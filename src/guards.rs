//! Refusals that protect the user's data, all checked before any output file is
//! created so a rejected run leaves nothing behind.

use std::io::IsTerminal;

use crate::config::Config;
use crate::io;

/// True iff writing `fmt`'s bytes to stdout would dump binary (BAM) or gzip
/// (FASTQ.gz) data into an interactive terminal: never useful output, and
/// almost always a forgotten `-o`/redirect. Plain FASTQ text is always fine.
/// Pure (no I/O) so it's trivial to unit-test without a real TTY.
pub(crate) fn binary_to_terminal(
    output_is_stdout: bool,
    fmt: io::Format,
    stdout_is_tty: bool,
) -> bool {
    output_is_stdout
        && stdout_is_tty
        && matches!(
            fmt,
            io::Format::Bam | io::Format::FastqGz | io::Format::FastqBgzf
        )
}

/// Reject binary output to an interactive terminal before creating a writer.
/// Report-only inference is exempt because it emits textual FASTA and exits
/// before workflow dispatch.
pub(crate) fn guard_stdout_binary(cfg: &Config, out_fmt: io::Format) -> anyhow::Result<()> {
    if cfg.adapter_infer.is_report() {
        return Ok(());
    }
    let stdout_is_tty = std::io::stdout().is_terminal();
    if binary_to_terminal(cfg.io.output.is_none(), out_fmt, stdout_is_tty) {
        let ext = match out_fmt {
            io::Format::Bam => "bam",
            io::Format::FastqGz => "fastq.gz",
            io::Format::FastqBgzf => "fastq.bgz",
            io::Format::Fastq => "fastq", // unreachable via binary_to_terminal, kept exhaustive
        };
        anyhow::bail!(
            "refusing to write {} to a terminal; redirect to a file/pipe (e.g. `> out.{ext}`) \
             or pass -o",
            out_fmt.label()
        );
    }
    Ok(())
}

/// Every file the run writes, checked against every file it reads.
///
/// `whittle` streams its input, so `File::create` truncating a file that is still
/// being read destroys it: a plain FASTQ run would emit an empty file and exit 0.
/// The same applies to `--summary-json`, which is written last and would replace
/// an input, the just-written output, or a folder member with JSON.
///
/// `extra_inputs` carries the folder-mode member files, which are not reachable
/// from `cfg` alone.
pub(crate) fn guard_output_collisions(
    cfg: &Config,
    extra_inputs: &[std::path::PathBuf],
) -> anyhow::Result<()> {
    let mut reads: Vec<(&str, &std::path::Path)> = Vec::new();
    if let Some(p) = cfg.io.input.as_deref() {
        reads.push(("the input", p));
    }
    if let Some(ac) = cfg.adapter_fasta.as_deref() {
        reads.push(("--adapter-fasta", ac));
    }
    for p in extra_inputs {
        reads.push(("an input file in the directory", p.as_path()));
    }

    for (what, dest) in [
        ("the output", cfg.io.output.as_deref()),
        ("--summary-json", cfg.summary_json.as_deref()),
    ] {
        let Some(dest) = dest else { continue };
        for (label, src) in &reads {
            if same_path(src, dest) {
                anyhow::bail!(
                    "{what} ({}) and {label} are the same file; whittle streams its input and \
                     would overwrite it, so write to a different path",
                    dest.display()
                );
            }
        }
        // With no `-i`, the input is stdin, which has no path to compare. Its
        // file descriptor still resolves to an inode, so a shell redirect from
        // the very file being written is caught here.
        if cfg.io.input.is_none() && stdin_is(dest) {
            anyhow::bail!(
                "{what} ({}) and the file being read on stdin are the same file; whittle \
                 streams its input and would truncate it mid-read, so write to a different path",
                dest.display()
            );
        }
    }

    // The summary is written after the output file, so it would clobber it.
    if let (Some(out), Some(sum)) = (cfg.io.output.as_deref(), cfg.summary_json.as_deref())
        && same_path(out, sum)
    {
        anyhow::bail!(
            "--summary-json ({}) and the output are the same file; the summary would replace \
             the trimmed reads with JSON",
            sum.display()
        );
    }
    Ok(())
}

/// Whether `path` names the same file the process has open on stdin.
///
/// A shell redirect (`whittle -o reads.fastq < reads.fastq`) leaves no path for
/// the same-file check to compare, but fd 0 still resolves to the inode.
#[cfg(unix)]
fn stdin_is(path: &std::path::Path) -> bool {
    use std::os::fd::AsFd;
    use std::os::unix::fs::MetadataExt;

    // `try_clone_to_owned` dups fd 0, so stdin itself is never closed here, and
    // it needs no `unsafe` (which this crate forbids).
    let Ok(dup) = std::io::stdin().as_fd().try_clone_to_owned() else {
        return false;
    };
    let Ok(m0) = std::fs::File::from(dup).metadata() else {
        return false;
    };
    let Ok(mp) = std::fs::metadata(path) else {
        return false;
    };
    // Only a regular file can collide; a pipe or tty shares no inode with a path.
    m0.is_file() && m0.dev() == mp.dev() && m0.ino() == mp.ino()
}

#[cfg(not(unix))]
fn stdin_is(_path: &std::path::Path) -> bool {
    false
}

/// Whether two paths resolve to the same file. Canonicalizes both so symlinks
/// and `./`-style aliasing are caught; the output usually does not exist yet, so
/// it falls back to canonicalizing the parent directory and re-joining the file
/// name. Conservative: any resolution failure yields `false` (don't block a run
/// on a path that cannot be resolved).
pub(crate) fn same_path(a: &std::path::Path, b: &std::path::Path) -> bool {
    // When both paths already exist, an inode+device match is definitive, and it
    // also catches hard links to one inode, which path canonicalization misses.
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(ma), Ok(mb)) = (std::fs::metadata(a), std::fs::metadata(b))
            && ma.dev() == mb.dev()
            && ma.ino() == mb.ino()
        {
            return true;
        }
    }
    fn resolve(p: &std::path::Path) -> Option<std::path::PathBuf> {
        if let Ok(c) = std::fs::canonicalize(p) {
            return Some(c);
        }
        let file = p.file_name()?;
        let parent = match p.parent() {
            Some(par) if !par.as_os_str().is_empty() => par,
            _ => std::path::Path::new("."),
        };
        std::fs::canonicalize(parent).ok().map(|c| c.join(file))
    }
    match (resolve(a), resolve(b)) {
        (Some(x), Some(y)) => x == y,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_to_terminal_flags_bam_on_a_tty_stdout() {
        assert!(binary_to_terminal(true, io::Format::Bam, true));
    }

    #[test]
    fn binary_to_terminal_flags_fastq_gz_on_a_tty_stdout() {
        assert!(binary_to_terminal(true, io::Format::FastqGz, true));
    }

    #[test]
    fn binary_to_terminal_allows_plain_fastq() {
        // Plain text FASTQ on a terminal is normal/expected output.
        assert!(!binary_to_terminal(true, io::Format::Fastq, true));
    }

    #[test]
    fn binary_to_terminal_allows_when_output_file_given() {
        // -o was given, so `output_is_stdout` is false regardless of format.
        assert!(!binary_to_terminal(false, io::Format::Bam, true));
    }

    #[test]
    fn binary_to_terminal_allows_when_not_a_tty() {
        // Redirected to a file/pipe: not a terminal, so it's fine.
        assert!(!binary_to_terminal(true, io::Format::Bam, false));
        assert!(!binary_to_terminal(true, io::Format::FastqGz, false));
    }
}
