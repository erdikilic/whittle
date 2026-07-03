pub mod cli;
pub mod config;
pub mod filter;
pub mod io;
pub mod pipeline;
pub mod qual;
pub mod record;
pub mod trim;

pub use config::Config;

/// Top-level entry point: dispatch on the resolved input/output formats and run
/// the matching pipeline. Plan 1 implements only the FASTQ path; Plan 2 adds BAM.
pub fn run(_cfg: Config) -> anyhow::Result<()> {
    anyhow::bail!("not yet implemented")
}
