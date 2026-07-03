use std::path::PathBuf;

use clap::Parser;

use crate::config::{Config, IoConfig};
use crate::filter::FilterConfig;
use crate::io::Format;
use crate::qual::QualMode;
use crate::trim::{QualityOp, TrimPlan};

#[derive(Parser, Debug)]
#[command(author, version, about = "Tag-aware long-read trimmer", long_about = None)]
struct Cli {
    #[arg(short = 'i', long, help_heading = "Setup")]
    input: Option<PathBuf>,
    #[arg(short = 'o', long, help_heading = "Setup")]
    output: Option<PathBuf>,
    #[arg(long, value_enum, help_heading = "Setup")]
    in_format: Option<FormatArg>,
    #[arg(long, value_enum, help_heading = "Setup")]
    out_format: Option<FormatArg>,
    #[arg(short = 't', long, default_value_t = 4, help_heading = "Setup")]
    threads: usize,

    #[arg(short = 'l', long, default_value_t = 1, help_heading = "Filtering")]
    min_length: usize,
    #[arg(short = 'L', long, help_heading = "Filtering")]
    max_length: Option<usize>,
    #[arg(short = 'q', long, default_value_t = 0.0, help_heading = "Filtering")]
    min_qual: f64,
    #[arg(short = 'Q', long, default_value_t = 1000.0, help_heading = "Filtering")]
    max_qual: f64,
    #[arg(short = 'g', long, help_heading = "Filtering")]
    min_gc: Option<f64>,
    #[arg(short = 'G', long, help_heading = "Filtering")]
    max_gc: Option<f64>,
    #[arg(short = 'm', long, value_enum, default_value_t = QualModeArg::Mean, help_heading = "Filtering")]
    qual_mode: QualModeArg,

    #[arg(short = 'H', long, default_value_t = 0, help_heading = "Trimming")]
    head_crop: usize,
    #[arg(short = 'T', long, default_value_t = 0, help_heading = "Trimming")]
    tail_crop: usize,
    #[arg(long, help_heading = "Trimming")]
    trim_qual: Option<u8>,
    #[arg(long, help_heading = "Trimming")]
    best_segment: Option<u8>,
    #[arg(long, help_heading = "Trimming")]
    split_qual: Option<u8>,
    #[arg(long, default_value_t = 1, help_heading = "Trimming")]
    split_window: usize,
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
enum FormatArg { Fastq, FastqGz, Bam }

impl From<FormatArg> for Format {
    fn from(f: FormatArg) -> Self {
        match f {
            FormatArg::Fastq => Format::Fastq,
            FormatArg::FastqGz => Format::FastqGz,
            FormatArg::Bam => Format::Bam,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
enum QualModeArg { Mean, Arithmetic, Median }

impl From<QualModeArg> for QualMode {
    fn from(m: QualModeArg) -> Self {
        match m {
            QualModeArg::Mean => QualMode::Mean,
            QualModeArg::Arithmetic => QualMode::Arithmetic,
            QualModeArg::Median => QualMode::Median,
        }
    }
}

pub fn parse() -> anyhow::Result<Config> {
    let c = Cli::parse();

    // Mutual exclusion of the three quality trim ops.
    let n_quality = [c.trim_qual.is_some(), c.best_segment.is_some(), c.split_qual.is_some()]
        .iter()
        .filter(|&&b| b)
        .count();
    if n_quality > 1 {
        anyhow::bail!("--trim-qual, --best-segment and --split-qual are mutually exclusive");
    }
    let quality = if let Some(q) = c.trim_qual {
        Some(QualityOp::TrimQual(q))
    } else if let Some(q) = c.best_segment {
        Some(QualityOp::BestSegment(q))
    } else if let Some(q) = c.split_qual {
        Some(QualityOp::Split { cutoff: q, window: c.split_window })
    } else {
        None
    };

    Ok(Config {
        io: IoConfig {
            input: c.input,
            output: c.output,
            in_format: c.in_format.map(Into::into),
            out_format: c.out_format.map(Into::into),
        },
        filter: FilterConfig {
            min_length: c.min_length,
            max_length: c.max_length.unwrap_or(usize::MAX),
            min_qual: c.min_qual,
            max_qual: c.max_qual,
            min_gc: c.min_gc,
            max_gc: c.max_gc,
            qual_mode: c.qual_mode.into(),
        },
        trim: TrimPlan { head: c.head_crop, tail: c.tail_crop, quality },
        threads: c.threads.max(1),
    })
}
