pub mod cli;
pub mod config;
pub mod filter;
pub mod io;
pub mod mods;
pub mod pipeline;
pub mod qual;
pub mod record;
pub mod trim;

pub use config::Config;

use std::io::{BufReader, BufWriter, Read, Write};

use gzp::par::compress::{ParCompress, ParCompressBuilder};
use gzp::{Compression, ZWriter, deflate::Gzip};

/// Top-level entry point. Dispatches on the input: a directory triggers
/// folder-merge (all read files in it merged into one output); otherwise a
/// single file / stdin is trimmed. FASTQ and unaligned BAM are supported.
pub fn run(cfg: Config) -> anyhow::Result<()> {
    use io::Format;

    let in_path = cfg.io.input.as_deref();

    if let Some(p) = in_path
        && p.is_dir()
    {
        return run_folder(p, &cfg);
    }

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
                // A single `read()` may return fewer than 4 bytes; loop to fill.
                let mut probe = [0u8; 4];
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
            }
        },
    };
    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), in_fmt));

    // BAM dispatch, and the cross-format BAM<->FASTQ rejection, happen before
    // creating/truncating the output file so a rejected run never leaves a
    // stray 0-byte file behind. Only the (Fastq*, Fastq*) combinations fall
    // through to the FASTQ path below.
    match (in_fmt, out_fmt) {
        (Format::Bam, Format::Bam) => {
            let (header, records) = io::bam::reader(in_path)?;
            // Provenance: append our @PG line to a cloned header before writing.
            let out_header = provenance_header(header);
            let mut writer = io::bam::writer(cfg.io.output.as_deref(), &out_header)?;
            let stats = pipeline::run_bam(&out_header, records, &mut writer, &cfg)?;
            // Explicitly finish (final bgzf block + EOF marker) instead of relying
            // on `Drop`, whose `try_finish` error is silently discarded — an I/O
            // failure on final flush (e.g. ENOSPC) would otherwise yield a
            // truncated BAM with a success exit code.
            writer.try_finish()?;
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            return Ok(());
        }
        (Format::Bam, _) | (_, Format::Bam) => {
            anyhow::bail!("cross-format BAM<->FASTQ conversion is not supported in v1")
        }
        _ => {}
    }

    let mut writer = fastq_writer(&cfg, out_fmt)?;

    let gz_in = matches!(in_fmt, Format::FastqGz);
    let records = io::fastq::reader_from(source, gz_in);
    let stats = pipeline::run_fastq(records, &mut writer, &cfg)?;
    writer.finish()?;
    eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
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
            }
            FastqOut::Gz(mut w) => {
                w.finish()?;
                Ok(())
            }
        }
    }
}

/// Build the FASTQ output writer: a file or stdout, wrapped in a parallel gzip
/// encoder (`gzp`, using `cfg.threads` worker threads) when the output format
/// is `FastqGz`, else a plain buffered writer.
fn fastq_writer(cfg: &Config, out_fmt: io::Format) -> anyhow::Result<FastqOut> {
    let base: Box<dyn Write + Send> = match cfg.io.output.as_deref() {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    };
    if matches!(out_fmt, io::Format::FastqGz) {
        let w = ParCompressBuilder::<Gzip>::new()
            .num_threads(cfg.threads.max(1))
            .unwrap()
            .compression_level(Compression::new(6))
            .from_writer(base);
        Ok(FastqOut::Gz(w))
    } else {
        Ok(FastqOut::Plain(BufWriter::new(base)))
    }
}

/// Folder-merge mode: `-i <dir>`. Classify the directory into one format family,
/// then merge all its read files into a single trimmed output using the same
/// pipelines as the single-file path.
fn run_folder(dir: &std::path::Path, cfg: &Config) -> anyhow::Result<()> {
    use io::Format;

    let (family, paths) = io::dir::classify(dir)?;
    eprintln!("Merging {} {:?} file(s) from {}", paths.len(), family, dir.display());
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
            let mut writer = fastq_writer(cfg, out_fmt)?;
            let records = io::dir::fastq_records(&paths);
            let stats = pipeline::run_fastq(records, &mut writer, cfg)?;
            writer.finish()?;
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            Ok(())
        }
        io::dir::Family::Bam => {
            if !matches!(out_fmt, Format::Bam) {
                anyhow::bail!(
                    "cross-format conversion (BAM folder to FASTQ) is not supported in v1"
                );
            }
            let (header, records) = io::dir::bam_reader(&paths)?;
            let out_header = provenance_header(header);
            let mut writer = io::bam::writer(cfg.io.output.as_deref(), &out_header)?;
            let stats = pipeline::run_bam(&out_header, records, &mut writer, cfg)?;
            writer.try_finish()?;
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            Ok(())
        }
    }
}

/// Append an `@PG` provenance record (`ID:chopping`, program name + version) to a
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
        .insert(tag::NAME, "chopping")
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .build();

    if let Ok(program) = program {
        let _ = header.programs_mut().add("chopping", program);
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
    use super::*;
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::Program;
    use noodles_sam::header::record::value::map::program::tag;

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
        // since the chain is dangling it must skip adding the `chopping`
        // `@PG` line entirely.
        let out_header = provenance_header(header);

        assert!(
            !out_header.programs().as_ref().contains_key(&b"chopping"[..]),
            "expected no chopping @PG line to be added when the existing chain is dangling"
        );
    }

    /// Companion positive-path test: a plain header with no dangling `@PG`
    /// chain must still get the `chopping` provenance record added, so the
    /// dangling-chain guard doesn't accidentally suppress the common case.
    #[test]
    fn provenance_header_adds_chopping_program_on_clean_header() {
        let header = noodles_sam::Header::default();
        assert!(!has_dangling_program_chain(&header));

        let out_header = provenance_header(header);

        assert!(
            out_header
                .programs()
                .roots()
                .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"chopping"),
            "expected an @PG record with ID chopping in the output header, got {:?}",
            out_header.programs()
        );
    }
}
