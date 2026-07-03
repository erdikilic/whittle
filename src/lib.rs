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

/// Top-level entry point: dispatch on the resolved input/output formats and run
/// the matching pipeline. Plan 1 implements only the FASTQ path; Plan 2 adds BAM.
pub fn run(cfg: Config) -> anyhow::Result<()> {
    use io::Format;

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
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            return Ok(());
        }
        (Format::Bam, _) | (_, Format::Bam) => {
            anyhow::bail!("cross-format BAM<->FASTQ conversion is not supported in v1")
        }
        _ => {}
    }

    let base_writer: Box<dyn Write + Send> = match cfg.io.output.as_deref() {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    };
    let writer_inner: Box<dyn Write + Send> = if matches!(out_fmt, Format::FastqGz) {
        Box::new(flate2::write::GzEncoder::new(base_writer, flate2::Compression::default()))
    } else {
        base_writer
    };
    let mut writer = BufWriter::new(writer_inner);

    let gz_in = matches!(in_fmt, Format::FastqGz);
    let records = io::fastq::reader_from(source, gz_in);
    let stats = pipeline::run_fastq(records, &mut writer, &cfg)?;
    writer.flush()?;
    // Drop the writer (finishing the GzEncoder, if any, and flushing the
    // BufWriter) before returning so all bytes are on disk / stdout.
    drop(writer);
    eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
    Ok(())
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
