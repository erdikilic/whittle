use std::path::PathBuf;

use clap::Parser;

use crate::config::{Config, FastqTags, IoConfig};
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
    #[arg(long, default_value = "all", help_heading = "Setup")]
    fastq_tags: String,
    /// DEFLATE compression level for compressed output (bgzf for BAM, gzip for
    /// FASTQ.gz). Lower = faster/larger. Ignored for plain FASTQ.
    #[arg(short = 'c', long, default_value_t = 6, help_heading = "Setup")]
    compression_level: u8,

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
    /// Keep ONT signal tags consistent through trimming (slice `mv`, update
    /// `ts`/`ns`/`sp`/`pi`) for signal-aware tools (Remora, Clair3 v2), instead
    /// of dropping them. BAM→BAM only.
    #[arg(long, help_heading = "Trimming")]
    update_moves: bool,
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
    // bgzf (libdeflate) accepts up to 12 and gzip up to 9; cap at the common 0-9
    // so a single flag is valid for both compressed output formats.
    if c.compression_level > 9 {
        anyhow::bail!(
            "--compression-level must be between 0 and 9 (got {})",
            c.compression_level
        );
    }

    // Reject contradictory or out-of-domain filter bounds up front, rather than
    // silently keeping zero reads and exiting successfully.
    let max_length = c.max_length.unwrap_or(usize::MAX);
    if c.min_length > max_length {
        anyhow::bail!(
            "--min-length ({}) must not exceed --max-length ({max_length})",
            c.min_length
        );
    }
    if c.min_qual.is_nan() || c.max_qual.is_nan() {
        anyhow::bail!("--min-qual and --max-qual must be numbers (got NaN)");
    }
    if c.min_qual > c.max_qual {
        anyhow::bail!(
            "--min-qual ({}) must not exceed --max-qual ({})",
            c.min_qual,
            c.max_qual
        );
    }
    if let Some(g) = c.min_gc
        && !(0.0..=1.0).contains(&g)
    {
        anyhow::bail!("--min-gc ({g}) must be a fraction between 0 and 1");
    }
    if let Some(g) = c.max_gc
        && !(0.0..=1.0).contains(&g)
    {
        anyhow::bail!("--max-gc ({g}) must be a fraction between 0 and 1");
    }
    if let (Some(a), Some(b)) = (c.min_gc, c.max_gc)
        && a > b
    {
        anyhow::bail!("--min-gc ({a}) must not exceed --max-gc ({b})");
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
    let fastq_tags = FastqTags::parse(&c.fastq_tags)?;

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
        fastq_tags,
        render_workers: 0,
        compression_level: c.compression_level,
        update_moves: c.update_moves,
    })
}

/// Build a Config directly (used by integration tests). head/tail are fixed crops.
#[doc(hidden)]
pub fn config_for_test(
    input: &std::path::Path,
    output: &std::path::Path,
    head_crop: usize,
    tail_crop: usize,
) -> Config {
    config_for_test_threads(input, output, head_crop, tail_crop, 1)
}

/// Same as `config_for_test`, but with an explicit thread count (used by tests
/// that need to exercise the parallel BAM dispatch, e.g. a `t8` oracle run).
#[doc(hidden)]
pub fn config_for_test_threads(
    input: &std::path::Path,
    output: &std::path::Path,
    head_crop: usize,
    tail_crop: usize,
    threads: usize,
) -> Config {
    Config {
        io: IoConfig {
            input: Some(input.to_path_buf()),
            output: Some(output.to_path_buf()),
            in_format: Some(Format::Bam),
            out_format: Some(Format::Bam),
        },
        filter: FilterConfig {
            min_length: 1,
            max_length: usize::MAX,
            min_qual: 0.0,
            max_qual: 1000.0,
            min_gc: None,
            max_gc: None,
            qual_mode: QualMode::Mean,
        },
        trim: TrimPlan { head: head_crop, tail: tail_crop, quality: None },
        threads: threads.max(1),
        fastq_tags: FastqTags::All,
        render_workers: 0,
        compression_level: 6,
        update_moves: false,
    }
}
