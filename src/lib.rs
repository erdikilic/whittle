pub mod cli;
pub mod config;
pub mod filter;
pub mod io;
pub mod mods;
pub mod pipeline;
pub mod qual;
pub mod record;
pub mod trim;

use std::io::{BufReader, BufWriter, Read, Write};

pub use config::Config;
use gzp::deflate::Gzip;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use gzp::{Compression, ZWriter};

/// Top-level entry point. Dispatches on the input: a directory triggers
/// folder-merge (all read files in it merged into one output); otherwise a
/// single file / stdin is trimmed. FASTQ and unaligned BAM are supported.
pub fn run(cfg: Config) -> anyhow::Result<()> {
    use config::EncodeKind;
    use io::Format;

    let mut cfg = cfg;

    // Scoped so the borrow of `cfg.io.input` ends before `run_folder` needs
    // `&mut cfg` — the directory path itself is cloned out first.
    if let Some(dir) = cfg
        .io
        .input
        .as_deref()
        .filter(|p| p.is_dir())
        .map(|p| p.to_path_buf())
    {
        return run_folder(&dir, &mut cfg);
    }

    // Refuse to read and write the same file: `whittle` streams the input, so
    // truncating it on `File::create` before it is fully read destroys the data
    // (a plain FASTQ run would silently emit an empty file with a success exit).
    if let (Some(inp), Some(outp)) = (cfg.io.input.as_deref(), cfg.io.output.as_deref())
        && same_path(inp, outp)
    {
        anyhow::bail!(
            "input and output are the same file ({}); whittle streams the input and \
             would truncate it before reading — write to a different path",
            outp.display()
        );
    }

    let in_path = cfg.io.input.as_deref();

    // Open the input (file or stdin) up front so format detection can sniff
    // its first bytes without losing them: any bytes consumed while sniffing
    // get prepended back via a Cursor+chain before the FASTQ reader is built.
    let raw: Box<dyn Read + Send> = match in_path {
        Some(p) => Box::new(std::fs::File::open(p)?),
        None => Box::new(std::io::stdin()),
    };
    let mut source: Box<dyn Read + Send> = Box::new(BufReader::new(raw));

    let in_fmt = match cfg.io.in_format {
        Some(f) => f,
        None => match in_path.and_then(io::from_extension) {
            Some(f) => f,
            None => {
                // Probe enough bytes to see a full BGZF block header (18 bytes),
                // so a BAM read from stdin/an unknown extension is told apart from
                // gzipped FASTQ. A single `read()` may return fewer; loop to fill.
                let mut probe = [0u8; 18];
                let mut n = 0;
                while n < probe.len() {
                    let r = source.read(&mut probe[n..])?;
                    if r == 0 {
                        break;
                    }
                    n += r;
                }
                let fmt = io::detect_input(in_path, &probe[..n])?;
                source = Box::new(std::io::Cursor::new(probe[..n].to_vec()).chain(source));
                fmt
            },
        },
    };
    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), in_fmt));

    // BAM dispatch happens before creating/truncating the output file, and so
    // do the FASTQ->BAM rejection and the BAM->FASTQ conversion, so a rejected
    // run never leaves a stray 0-byte file behind. Only the (Fastq*, Fastq*)
    // combinations fall through to the FASTQ path below.
    match (in_fmt, out_fmt) {
        (Format::Bam, Format::Bam) => {
            note_tags_ignored(&cfg, in_fmt, out_fmt);
            let b = config::thread_budget(cfg.threads, true, EncodeKind::Bgzf);
            // Read from `source` (not by re-opening `in_path`): for a stdin BAM the
            // sniff bytes were already consumed and chained back into `source`, so
            // re-opening stdin would drop the BGZF header. For a file, `source` is
            // the same handle positioned at the start.
            let (header, records) = io::bam::reader_from(source, b.decode)?;
            // Provenance: append our @PG line to a cloned header before writing.
            let out_header = provenance_header(header);
            let mut sink = io::bam::writer(
                cfg.io.output.as_deref(),
                &out_header,
                b.encode,
                cfg.compression_level,
            )?;
            cfg.render_workers = b.render;
            let stats = pipeline::run_bam(&out_header, records, &mut sink, &cfg)?;
            // Explicitly finish (final bgzf block + EOF marker) instead of relying
            // on `Drop`, whose `try_finish` error is silently discarded — an I/O
            // failure on final flush (e.g. ENOSPC) would otherwise yield a
            // truncated BAM with a success exit code.
            sink.finish()?;
            eprint_run_summary(&stats);
            return Ok(());
        },
        (Format::Bam, Format::Fastq | Format::FastqGz) => {
            let encode = if matches!(out_fmt, Format::FastqGz) {
                EncodeKind::Gzip
            } else {
                EncodeKind::None
            };
            let b = config::thread_budget(cfg.threads, true, encode);
            // See the note in the (Bam, Bam) arm: read from the chained `source`.
            let (_header, records) = io::bam::reader_from(source, b.decode)?;
            let mut writer = fastq_writer(&cfg, out_fmt, b.encode)?;
            cfg.render_workers = b.render;
            let stats = pipeline::run_bam_to_fastq(records, &mut writer, &cfg)?;
            writer.finish()?;
            eprint_run_summary(&stats);
            return Ok(());
        },
        (Format::Fastq | Format::FastqGz, Format::Bam) => {
            anyhow::bail!("cross-format FASTQ->BAM conversion is not supported")
        },
        _ => {},
    }

    note_tags_ignored(&cfg, in_fmt, out_fmt);
    let encode = if matches!(out_fmt, Format::FastqGz) {
        EncodeKind::Gzip
    } else {
        EncodeKind::None
    };
    let b = config::thread_budget(cfg.threads, false, encode);
    let mut writer = fastq_writer(&cfg, out_fmt, b.encode)?;
    cfg.render_workers = b.render;

    let gz_in = matches!(in_fmt, Format::FastqGz);
    let records = io::fastq::reader_from(source, gz_in);
    let stats = pipeline::run_fastq(records, &mut writer, &cfg)?;
    writer.finish()?;
    eprint_run_summary(&stats);
    Ok(())
}

/// FASTQ output writer: either a plain buffered writer, or a `gzp` parallel
/// gzip writer (used only when the output format is explicitly `FastqGz` —
/// see `io::resolve_output`, which no longer auto-compresses).
///
/// `gzp`'s `ParCompress` REQUIRES an explicit `finish()` call to flush its
/// final (possibly partial) compressed block plus the gzip footer/checksum:
/// its `Write` impl only ever hands off *full* buffered chunks to the
/// compressor threads, so the tail end only ever gets flushed by `finish()`,
/// never by `flush()`. `ParCompress`'s own `Drop` impl does call `finish()` as
/// a backstop if it's still live, but it `.unwrap()`s the result — any I/O
/// error at that point becomes an uncatchable panic instead of a propagated
/// `anyhow::Result::Err`. Calling `finish()` explicitly, as the single seam
/// both callers below go through instead of `flush()` + `drop()`, keeps that
/// failure mode as an ordinary `Err`.
enum FastqOut {
    Plain(BufWriter<Box<dyn Write + Send>>),
    Gz(ParCompress<'static, Gzip, Box<dyn Write + Send>>),
}

impl Write for FastqOut {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            FastqOut::Plain(w) => w.write(buf),
            FastqOut::Gz(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            FastqOut::Plain(w) => w.flush(),
            FastqOut::Gz(w) => w.flush(),
        }
    }
}

impl FastqOut {
    /// Finalize: for gz, flush the final block + gzip footer via gzp's
    /// `ZWriter::finish` (required — see the type's docs above); for plain,
    /// flush the `BufWriter`. Must be called before returning success.
    fn finish(self) -> anyhow::Result<()> {
        match self {
            FastqOut::Plain(mut w) => {
                w.flush()?;
                Ok(())
            },
            FastqOut::Gz(mut w) => {
                w.finish()?;
                Ok(())
            },
        }
    }
}

/// Build the FASTQ output writer: a file or stdout, wrapped in a parallel gzip
/// encoder (`gzp`, using `gz_workers` worker threads — the caller's
/// workload-aware ENCODE share of the `-t` budget) when the output format is
/// `FastqGz`, else a plain buffered writer. `gz_workers` is only read on the
/// `FastqGz` branch, so plain-output callers may pass any value.
fn fastq_writer(cfg: &Config, out_fmt: io::Format, gz_workers: usize) -> anyhow::Result<FastqOut> {
    let base: Box<dyn Write + Send> = match cfg.io.output.as_deref() {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    };
    if matches!(out_fmt, io::Format::FastqGz) {
        let w = ParCompressBuilder::<Gzip>::new()
            .num_threads(gz_workers)
            .unwrap()
            .compression_level(Compression::new(cfg.compression_level as u32))
            .from_writer(base);
        Ok(FastqOut::Gz(w))
    } else {
        Ok(FastqOut::Plain(BufWriter::new(base)))
    }
}

/// Folder-merge mode: `-i <dir>`. Classify the directory into one format family,
/// then merge all its read files into a single trimmed output using the same
/// pipelines as the single-file path.
fn run_folder(dir: &std::path::Path, cfg: &mut Config) -> anyhow::Result<()> {
    use config::EncodeKind;
    use io::Format;

    // Pass the output path so `classify` can hard-error if `-o` names a read file
    // inside `-i <dir>` — it could be a real input or a stale prior output, and
    // overwriting either while merging the rest is silent data loss. The merged
    // output must live outside the input directory.
    let (family, paths) = io::dir::classify(dir, cfg.io.output.as_deref())?;
    eprintln!(
        "Merging {} {:?} file(s) from {}",
        paths.len(),
        family,
        dir.display()
    );
    let family_fmt = match family {
        io::dir::Family::Fastq => Format::Fastq,
        io::dir::Family::Bam => Format::Bam,
    };
    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), family_fmt));

    match family {
        io::dir::Family::Fastq => {
            if matches!(out_fmt, Format::Bam) {
                anyhow::bail!(
                    "cross-format conversion (FASTQ folder to BAM) is not supported in v1"
                );
            }
            note_tags_ignored(cfg, family_fmt, out_fmt);
            let encode = if matches!(out_fmt, Format::FastqGz) {
                EncodeKind::Gzip
            } else {
                EncodeKind::None
            };
            let b = config::thread_budget(cfg.threads, false, encode);
            let mut writer = fastq_writer(cfg, out_fmt, b.encode)?;
            cfg.render_workers = b.render;
            let records = io::dir::fastq_records(&paths);
            let stats = pipeline::run_fastq(records, &mut writer, cfg)?;
            writer.finish()?;
            eprint_run_summary(&stats);
            Ok(())
        },
        io::dir::Family::Bam => match out_fmt {
            Format::Bam => {
                note_tags_ignored(cfg, family_fmt, out_fmt);
                // Only the first file's header is written; warn if the others
                // declare different read groups (relevant only for BAM output).
                io::dir::warn_on_bam_header_mismatch(&paths);
                let b = config::thread_budget(cfg.threads, true, EncodeKind::Bgzf);
                let (header, records) = io::dir::bam_reader(&paths, b.decode)?;
                let out_header = provenance_header(header);
                let mut sink = io::bam::writer(
                    cfg.io.output.as_deref(),
                    &out_header,
                    b.encode,
                    cfg.compression_level,
                )?;
                cfg.render_workers = b.render;
                let stats = pipeline::run_bam(&out_header, records, &mut sink, cfg)?;
                sink.finish()?;
                eprint_run_summary(&stats);
                Ok(())
            },
            Format::Fastq | Format::FastqGz => {
                let encode = if matches!(out_fmt, Format::FastqGz) {
                    EncodeKind::Gzip
                } else {
                    EncodeKind::None
                };
                let b = config::thread_budget(cfg.threads, true, encode);
                let (_header, records) = io::dir::bam_reader(&paths, b.decode)?;
                let mut writer = fastq_writer(cfg, out_fmt, b.encode)?;
                cfg.render_workers = b.render;
                let stats = pipeline::run_bam_to_fastq(records, &mut writer, cfg)?;
                writer.finish()?;
                eprint_run_summary(&stats);
                Ok(())
            },
        },
    }
}

/// Whether two paths resolve to the same file. Canonicalizes both so symlinks
/// and `./`-style aliasing are caught; the output usually does not exist yet, so
/// it falls back to canonicalizing the parent directory and re-joining the file
/// name. Conservative: any resolution failure yields `false` (don't block a run
/// on a path we can't resolve).
pub(crate) fn same_path(a: &std::path::Path, b: &std::path::Path) -> bool {
    // When both paths already exist, an inode+device match is definitive — and it
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

/// Print the end-of-run summary: the kept/total count, plus a one-line advisory
/// if any read carried a malformed per-base tag (see `has_malformed_perbase_tag`).
fn eprint_run_summary(stats: &pipeline::Stats) {
    eprintln!(
        "Kept {} reads out of {}",
        stats.output_reads, stats.input_reads
    );
    if stats.malformed_tag_reads > 0 {
        eprintln!(
            "note: {} read(s) carried a per-base kinetics tag (ip/pw/fi/fp/ri/rp) whose \
             length did not match the sequence; left unchanged",
            stats.malformed_tag_reads
        );
    }
}

/// `--fastq-tags` only affects BAM→FASTQ output. When the user set a non-default
/// value (`none`/an explicit list) on any other path, emit a one-line stderr note
/// rather than silently ignoring it. (An explicit `all` is the default and stays
/// silent.)
fn note_tags_ignored(cfg: &Config, in_fmt: io::Format, out_fmt: io::Format) {
    if !matches!(cfg.fastq_tags, config::FastqTags::All) {
        eprintln!(
            "note: --fastq-tags applies only to BAM->FASTQ output; ignored for {in_fmt:?}->{out_fmt:?}"
        );
    }
}

/// Append an `@PG` provenance record (`ID:whittle`, program name + version) to a
/// cloned header before writing. Best-effort: `Programs::add` can fail (e.g. on a
/// duplicate ID), in which case the header is written unchanged — the `@PG` line
/// is cosmetic and must never block record output.
fn provenance_header(mut header: noodles_sam::Header) -> noodles_sam::Header {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::Program;
    use noodles_sam::header::record::value::map::program::tag;

    // `Programs::add` walks the existing `@PG` chain via `Programs::leaves`,
    // which indexes the program map directly and panics if any program's `PP`
    // (previous-program) field names an ID that isn't itself a program in the
    // header. Real-world uBAMs can have exactly this: e.g. an ONT/dorado file
    // put through `samtools sort`/`view`/`reset` observed with
    // `@PG ID:samtools PP:basecaller` where no `ID:basecaller` record survived
    // into the header. Since the `@PG` line is cosmetic, skip adding it rather
    // than let a merely-untidy header crash the whole run.
    if has_dangling_program_chain(&header) {
        return header;
    }

    let program = Map::<Program>::builder()
        .insert(tag::NAME, "whittle")
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .build();

    if let Ok(program) = program {
        let _ = header.programs_mut().add("whittle", program);
    }

    header
}

/// True if any `@PG` record's `PP` field references an ID that is not itself a
/// program in the header. `Programs::leaves` (used internally by
/// `Programs::add`) panics on such a chain instead of returning an error, so
/// this must be checked before calling `add`.
fn has_dangling_program_chain(header: &noodles_sam::Header) -> bool {
    use noodles_sam::header::record::value::map::program::tag;

    let programs = header.programs().as_ref();
    programs.values().any(|program| {
        program
            .other_fields()
            .get(&tag::PREVIOUS_PROGRAM_ID)
            .is_some_and(|previous_id| !programs.contains_key(previous_id))
    })
}

#[cfg(test)]
mod tests {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::Program;
    use noodles_sam::header::record::value::map::program::tag;

    use super::*;

    /// Regression test for `d481c48`: a header with a dangling `@PG PP:` chain
    /// (a `PP` value that names a program ID not present in the header) used
    /// to panic inside `noodles_sam::header::Programs::add` — called via
    /// `provenance_header` — because `Programs::leaves` indexes the program
    /// map directly by the `PP` id without checking it exists first. Real
    /// ONT/samtools headers hit this in the wild (see `d481c48`'s commit
    /// message). `provenance_header` must detect the dangling reference via
    /// `has_dangling_program_chain` and return the header unchanged instead
    /// of calling `Programs::add`.
    #[test]
    fn provenance_header_does_not_panic_on_dangling_pp_chain() {
        // "pg1" claims a previous program "ghost", but "ghost" is never
        // added to the header — a genuinely dangling reference.
        let dangling_program = Map::<Program>::builder()
            .insert(tag::PREVIOUS_PROGRAM_ID, "ghost")
            .build()
            .expect("valid PP field");

        let header = noodles_sam::Header::builder()
            .add_program("pg1", dangling_program)
            .build();

        // Sanity-check that the header really is dangling (i.e. this test
        // isn't accidentally exercising the clean path).
        assert!(has_dangling_program_chain(&header));

        // Pre-fix, this call panicked inside `Programs::add` -> `leaves`
        // -> `has_cycle`, which indexes the program map with the `PP` id
        // and panics when that id isn't a key (`ghost` isn't present here).
        // Post-fix, `provenance_header` must return without panicking, and
        // since the chain is dangling it must skip adding the `whittle`
        // `@PG` line entirely.
        let out_header = provenance_header(header);

        assert!(
            !out_header.programs().as_ref().contains_key(&b"whittle"[..]),
            "expected no whittle @PG line to be added when the existing chain is dangling"
        );
    }

    /// Companion positive-path test: a plain header with no dangling `@PG`
    /// chain must still get the `whittle` provenance record added, so the
    /// dangling-chain guard doesn't accidentally suppress the common case.
    #[test]
    fn provenance_header_adds_whittle_program_on_clean_header() {
        let header = noodles_sam::Header::default();
        assert!(!has_dangling_program_chain(&header));

        let out_header = provenance_header(header);

        assert!(
            out_header
                .programs()
                .roots()
                .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"whittle"),
            "expected an @PG record with ID whittle in the output header, got {:?}",
            out_header.programs()
        );
    }
}
