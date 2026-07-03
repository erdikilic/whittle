pub mod cli;
pub mod config;
pub mod filter;
pub mod io;
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

    // Reject BAM before creating/truncating the output file so a rejected run
    // never leaves a stray 0-byte file behind.
    if matches!(in_fmt, Format::Bam) || matches!(out_fmt, Format::Bam) {
        anyhow::bail!("BAM support arrives in Plan 2");
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
