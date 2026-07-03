pub mod cli;
pub mod config;
pub mod filter;
pub mod io;
pub mod pipeline;
pub mod qual;
pub mod record;
pub mod trim;

pub use config::Config;

use std::io::{BufWriter, Write};
use std::path::Path;

/// Top-level entry point: dispatch on the resolved input/output formats and run
/// the matching pipeline. Plan 1 implements only the FASTQ path; Plan 2 adds BAM.
pub fn run(cfg: Config) -> anyhow::Result<()> {
    use io::Format;

    let in_path = cfg.io.input.as_deref();
    // Sniff a few bytes for stdin/unknown-extension detection.
    let in_fmt = match cfg.io.in_format {
        Some(f) => f,
        None => {
            let sniff = peek_input(in_path)?;
            io::detect_input(in_path, &sniff)?
        }
    };
    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), in_fmt));

    let mut writer: BufWriter<Box<dyn Write>> = BufWriter::new(match cfg.io.output.as_deref() {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    });

    match (in_fmt, out_fmt) {
        (Format::Fastq, Format::Fastq)
        | (Format::FastqGz, Format::Fastq)
        | (Format::Fastq, Format::FastqGz)
        | (Format::FastqGz, Format::FastqGz) => {
            let gz_in = matches!(in_fmt, Format::FastqGz);
            let records = io::fastq::reader(in_path, gz_in)?;
            // gz output wrapping is added in a later task if out_fmt is FastqGz;
            // Plan 1 supports plain FASTQ output here, gz output in Task 9 note.
            let stats = pipeline::run_fastq(records, &mut writer, &cfg)?;
            writer.flush()?;
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            Ok(())
        }
        (Format::Bam, _) | (_, Format::Bam) => {
            anyhow::bail!("BAM support arrives in Plan 2")
        }
    }
}

fn peek_input(path: Option<&Path>) -> anyhow::Result<Vec<u8>> {
    use std::io::Read;
    let mut buf = vec![0u8; 4];
    if let Some(p) = path {
        let mut f = std::fs::File::open(p)?;
        let n = f.read(&mut buf)?;
        buf.truncate(n);
    }
    Ok(buf)
}
