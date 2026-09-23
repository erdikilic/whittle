//! Output collision guards and format- and tag-specific option validation.

use std::io::IsTerminal;
use std::path::Path;

use crate::config::Config;
use crate::io;

/// True iff writing `fmt`'s bytes to stdout would send binary (BAM) or gzip
/// (FASTQ.gz, FASTQ.bgz) data to an interactive terminal. Plain FASTQ text is
/// always allowed. Pure (no I/O), so it is unit-tested without a TTY.
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

/// Rejects binary output to an interactive terminal before creating a writer.
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
            io::Format::Fastq => "fastq", // unreachable: binary_to_terminal excludes plain FASTQ
        };
        anyhow::bail!(
            "refusing to write {} to a terminal; redirect to a file or pipe (e.g. `> out.{ext}`) \
             or pass -o",
            out_fmt.label()
        );
    }
    Ok(())
}

/// Rejects FASTQ input with BAM output: a FASTQ read carries no header from
/// which a BAM record could be built. `cli::parse` applies this to the formats
/// an explicit flag or a known extension names, and `run` applies it again
/// once detection has classified a stream or a directory, before any output is
/// announced or created.
pub(crate) fn guard_fastq_to_bam(in_fmt: io::Format, out_fmt: io::Format) -> anyhow::Result<()> {
    if in_fmt != io::Format::Bam && out_fmt == io::Format::Bam {
        anyhow::bail!(
            "FASTQ-to-BAM conversion is not supported (a FASTQ read carries no header for a BAM \
             record); write FASTQ output, or import with samtools import"
        );
    }
    Ok(())
}

/// Rejects `--update-moves` on any input format other than BAM: the ONT
/// signal tags are rewritten only on BAM-to-BAM output, so elsewhere the flag
/// would be accepted and silently do nothing. `cli::parse` applies this to the
/// format an explicit `--input-format` or a known extension names, and `run`
/// applies it again once detection has classified a stream.
pub(crate) fn guard_bam_only_flags(cfg: &Config, in_fmt: io::Format) -> anyhow::Result<()> {
    if in_fmt != io::Format::Bam && cfg.update_moves {
        anyhow::bail!(
            "--update-moves rewrites the ONT signal tags of BAM records and requires BAM input \
             (got {})",
            in_fmt.label()
        );
    }
    Ok(())
}

/// Rejects the flags that remove aux tags when the input carries none:
/// `--remove-tag` names tags. BAM records and tagged
/// FASTQ headers carry tags; a plain FASTQ does not, so there the flags would
/// be accepted and silently do nothing.
/// `run` applies this after the complete FASTQ stream has been inspected.
pub(crate) fn guard_tag_flags(
    cfg: &Config,
    in_fmt: io::Format,
    tagged: bool,
) -> anyhow::Result<()> {
    if in_fmt == io::Format::Bam || tagged {
        return Ok(());
    }
    let got = format!("{} without header tags", in_fmt.label());
    if !cfg.tag_filters.is_empty() {
        anyhow::bail!(
            "--tag-filter reads aux tags and requires BAM or tagged FASTQ input (got {got})"
        );
    }
    if !cfg.remove_tags.is_empty() {
        anyhow::bail!(
            "--remove-tag removes aux tags and requires BAM or tagged FASTQ input (got {got})"
        );
    }
    Ok(())
}

/// Checks every file the run writes against every file it reads and against
/// each other.
///
/// `whittle` streams its input, so `File::create` truncating a file that is still
/// being read destroys it: a plain FASTQ run would emit an empty file and exit 0.
/// The summary JSON is written after the reads, so it would replace an input,
/// the output written by the same run, or a folder member.
///
/// `extra_inputs` carries the folder-mode member files, which are not reachable
/// from `cfg` alone.
pub(crate) fn guard_output_collisions(
    cfg: &Config,
    extra_inputs: &[std::path::PathBuf],
) -> anyhow::Result<()> {
    let mut reads: Vec<(&str, &Path)> = Vec::new();
    if let Some(p) = cfg.io.input.as_deref() {
        reads.push(("the input", p));
    }
    if let Some(ac) = cfg.adapter_fasta.as_deref() {
        reads.push(("--adapter-fasta", ac));
    }
    for p in extra_inputs {
        reads.push(("an input file in the directory", p.as_path()));
    }

    let targets: Vec<(&str, &Path)> = cfg.write_targets().collect();
    for &(what, dest) in &targets {
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

    // Two artifacts on one path leave only the one written last. With no `-o`
    // the reads go to stdout, which a shell redirect can point at the very file
    // an artifact names, leaving no path for `same_path` to compare.
    for (i, &(later, dest)) in targets.iter().enumerate() {
        for &(earlier, other) in &targets[..i] {
            if same_path(other, dest) {
                anyhow::bail!(
                    "{later} ({}) and {earlier} are the same file; write each to its own path",
                    dest.display()
                );
            }
        }
        if cfg.io.output.is_none() && stdout_is(dest) {
            anyhow::bail!(
                "{later} ({}) and the output on stdout are the same file; the file would \
                 replace the trimmed reads",
                dest.display()
            );
        }
    }

    // A mistyped artifact path is caught during setup rather than after the
    // reads are written. Probing the parent directory instead of creating the
    // file leaves nothing behind on a run that writes no artifact
    // (`--adapter-report`). A permission failure still surfaces at write
    // time; a nonexistent directory does not.
    for &(flag, path) in &targets {
        if path.is_dir() {
            anyhow::bail!("{flag} {} is a directory", path.display());
        }
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        if !parent.is_dir() {
            anyhow::bail!(
                "{flag} {}: the directory {} does not exist",
                path.display(),
                parent.display()
            );
        }
    }
    Ok(())
}

/// Whether `path` names the same file the process has open on `fd`.
///
/// A shell redirect leaves no path for `same_path` to compare, but the
/// descriptor still resolves to an inode. Only a regular file can be clobbered,
/// so a pipe, tty, or `/dev/null` is never a match.
#[cfg(unix)]
fn redirects_to(fd: std::os::fd::BorrowedFd<'_>, path: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;

    // `try_clone_to_owned` dups the descriptor, so the original is never closed
    // here, and it needs no `unsafe` (which this crate forbids).
    let Ok(dup) = fd.try_clone_to_owned() else {
        return false;
    };
    let Ok(mine) = std::fs::File::from(dup).metadata() else {
        return false;
    };
    let Ok(theirs) = std::fs::metadata(path) else {
        return false;
    };
    mine.is_file() && mine.dev() == theirs.dev() && mine.ino() == theirs.ino()
}

/// Whether `path` is the file being read on stdin (`whittle -o x < x`).
#[cfg(unix)]
fn stdin_is(path: &Path) -> bool {
    use std::os::fd::AsFd;
    redirects_to(std::io::stdin().as_fd(), path)
}

/// Whether `path` is the file stdout is redirected to (`whittle --summary-json x > x`).
#[cfg(unix)]
fn stdout_is(path: &Path) -> bool {
    use std::os::fd::AsFd;
    redirects_to(std::io::stdout().as_fd(), path)
}

#[cfg(not(unix))]
fn stdin_is(_path: &Path) -> bool {
    false
}

#[cfg(not(unix))]
fn stdout_is(_path: &Path) -> bool {
    false
}

/// Whether two paths resolve to the same file. Canonicalizes both so symlinks
/// and `./`-style aliasing are caught; the output usually does not exist yet, so
/// it falls back to canonicalizing the parent directory and re-joining the file
/// name. Conservative: any resolution failure yields `false`, so an unresolvable
/// path never blocks a run.
pub(crate) fn same_path(a: &Path, b: &Path) -> bool {
    // Only a regular file has contents to destroy. `/dev/null` or a tty as both
    // endpoints is a discard, not a collision.
    for p in [a, b] {
        if let Ok(m) = std::fs::metadata(p)
            && !m.is_file()
        {
            return false;
        }
    }

    // When both paths exist, an inode and device match is definitive, and it
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
        // Plain-text FASTQ on a terminal is expected output.
        assert!(!binary_to_terminal(true, io::Format::Fastq, true));
    }

    #[test]
    fn binary_to_terminal_allows_when_output_file_given() {
        // -o was given, so `output_is_stdout` is false regardless of format.
        assert!(!binary_to_terminal(false, io::Format::Bam, true));
    }

    #[test]
    fn binary_to_terminal_allows_when_not_a_tty() {
        // Redirected to a file or pipe: not a terminal, so binary output is
        // allowed.
        assert!(!binary_to_terminal(true, io::Format::Bam, false));
        assert!(!binary_to_terminal(true, io::Format::FastqGz, false));
    }
}
