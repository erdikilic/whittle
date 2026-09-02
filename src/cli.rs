//! Command-line parsing: the clap definition, cross-flag validation, and
//! construction of the resolved `Config`.

use std::path::PathBuf;

use clap::Parser;

use crate::config::{
    AdapterInfer, AdapterInferAction, AdapterInferPolicy, Advisory, Config, FastqTags, IoConfig,
    ProgressMode,
};
use crate::filter::FilterConfig;
use crate::io::Format;
use crate::qual::QualMode;
use crate::trim::{QualityOp, TrimPlan};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    disable_version_flag = true,
    about = "Tag-aware long-read trimmer",
    long_about = None
)]
struct Cli {
    /// Print version information and exit.
    #[arg(long, action = clap::ArgAction::Version, help_heading = "Setup")]
    version: Option<bool>,
    /// Input FASTQ-family file, unaligned BAM, or directory; - means stdin.
    /// Defaults to stdin.
    #[arg(short = 'i', long, help_heading = "Setup")]
    input: Option<PathBuf>,
    /// Output file, whose extension selects the format; - means stdout.
    /// Defaults to stdout.
    #[arg(short = 'o', long, help_heading = "Setup")]
    output: Option<PathBuf>,
    /// Force the input format instead of detecting it from the path or stream.
    #[arg(long, value_enum, help_heading = "Setup")]
    in_format: Option<FormatArg>,
    /// Force the output format instead of selecting it from the output path.
    #[arg(long, value_enum, help_heading = "Setup")]
    out_format: Option<FormatArg>,
    /// Worker threads, at least 1; values above the CPU count are clamped to
    /// it. Defaults to all detected CPUs.
    #[arg(
        short = 't',
        long,
        value_parser = clap::value_parser!(u64).range(1..),
        help_heading = "Setup"
    )]
    threads: Option<u64>,
    /// Write records in input order when running with more than one thread.
    /// Without it, records are written as they finish, which is faster and uses
    /// less memory but is not reproducible between runs.
    #[arg(long, help_heading = "Setup")]
    ordered: bool,
    /// BAM auxiliary tags copied into BAM-to-FASTQ headers: all, none, or a
    /// list such as MM,ML,RG. Defaults to all.
    #[arg(long, default_value = "all", help_heading = "Setup")]
    fastq_tags: String,
    /// DEFLATE compression level (0-9) for compressed output: bgzf for BAM and
    /// .bgz, gzip for FASTQ.gz. Lower levels are faster and produce larger
    /// files. Ignored for plain FASTQ. Defaults to 4 for gzip FASTQ and 6 for
    /// BGZF.
    // bgzf (libdeflate) accepts up to 12 and gzip up to 9; the cap is the
    // common 0-9 so a single flag is valid for both compressed output formats.
    #[arg(
        short = 'c',
        long,
        value_parser = clap::value_parser!(u8).range(0..=9),
        help_heading = "Setup"
    )]
    compression_level: Option<u8>,
    /// Write a machine-readable JSON run summary (counters plus the resolved
    /// settings) to this path. Written even under --quiet.
    #[arg(long, value_name = "PATH", help_heading = "Setup")]
    summary_json: Option<PathBuf>,

    /// Increase logging detail: -v is debug, -vv is trace (at most two).
    /// Overridden by WHITTLE_LOG.
    #[arg(short = 'v', long, action = clap::ArgAction::Count, help_heading = "Logging")]
    verbose: u8,
    /// Silence progress and info output; warnings and errors still print.
    /// Conflicts with -v and --progress.
    #[arg(long, conflicts_with = "verbose", help_heading = "Logging")]
    quiet: bool,
    /// Progress reporting, independent of the log level: auto selects a bar on
    /// a terminal and periodic lines otherwise; bar and plain force one form;
    /// none disables progress and keeps the banner and summary. A bar falls
    /// back to lines under -v or WHITTLE_LOG. Defaults to auto.
    #[arg(
        long,
        value_enum,
        default_value_t = ProgressArg::Auto,
        conflicts_with = "quiet",
        help_heading = "Logging"
    )]
    progress: ProgressArg,

    /// Minimum post-trim segment length. Defaults to 1.
    #[arg(short = 'l', long, default_value_t = 1, help_heading = "Filtering")]
    min_length: usize,
    /// Maximum post-trim segment length.
    #[arg(short = 'L', long, help_heading = "Filtering")]
    max_length: Option<usize>,
    /// Minimum post-trim read quality under --qual-mode. Defaults to 0.
    #[arg(short = 'q', long, default_value_t = 0.0, help_heading = "Filtering")]
    min_qual: f64,
    /// Maximum post-trim read quality under --qual-mode. Defaults to 1000.
    #[arg(
        short = 'Q',
        long,
        default_value_t = 1000.0,
        help_heading = "Filtering"
    )]
    max_qual: f64,
    /// Minimum post-trim GC fraction (0 to 1).
    #[arg(short = 'g', long, help_heading = "Filtering")]
    min_gc: Option<f64>,
    /// Maximum post-trim GC fraction (0 to 1).
    #[arg(short = 'G', long, help_heading = "Filtering")]
    max_gc: Option<f64>,
    /// Read-quality summary used by the quality filters. Defaults to mean.
    #[arg(short = 'm', long, value_enum, default_value_t = QualModeArg::Mean, help_heading = "Filtering")]
    qual_mode: QualModeArg,

    /// Remove this many bases from the 5' end before other trimming. Defaults
    /// to 0.
    #[arg(short = 'H', long, default_value_t = 0, help_heading = "Trimming")]
    head_crop: usize,
    /// Remove this many bases from the 3' end before other trimming. Defaults
    /// to 0.
    #[arg(short = 'T', long, default_value_t = 0, help_heading = "Trimming")]
    tail_crop: usize,
    /// Trim low-quality bases from both ends until each boundary reaches Q.
    #[arg(long, help_heading = "Trimming")]
    qual_trim: Option<u8>,
    /// Keep the longest contiguous segment whose bases are all at least Q.
    #[arg(long, help_heading = "Trimming")]
    qual_best_segment: Option<u8>,
    /// Split at low-quality runs below Q and keep the surviving segments.
    #[arg(long, help_heading = "Trimming")]
    qual_split: Option<u8>,
    /// Tolerate low-quality runs shorter than this many bases when splitting.
    /// Requires --qual-split. Defaults to 1.
    #[arg(
        long,
        default_value_t = 1,
        requires = "qual_split",
        help_heading = "Trimming"
    )]
    qual_split_window: usize,
    /// Keep ONT signal tags consistent through trimming (slice mv, update ts,
    /// ns, sp and pi) for signal-aware tools such as Remora and Clair3 v2,
    /// instead of dropping them. BAM-to-BAM only.
    #[arg(long, help_heading = "Trimming")]
    update_moves: bool,

    /// Adapter FASTA; each sequence must be at least 11 bp. Enables adapter
    /// trimming.
    #[arg(short = 'a', long, help_heading = "Adapter trimming")]
    adapter_fasta: Option<PathBuf>,
    /// Built-in ONT adapter catalog. Enables adapter trimming. Defaults to none.
    #[arg(long, value_enum, default_value_t = AdapterPresetArg::None, help_heading = "Adapter trimming")]
    adapter_preset: AdapterPresetArg,
    /// End-match tolerance as a fraction of adapter length; interior splits use
    /// half. Requires an adapter source. Defaults to 0.2.
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_error_rate: Option<f64>,
    /// Bases at each read end searched for a terminal adapter. Requires an
    /// adapter source. Defaults to 150.
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_end_size: Option<usize>,
    /// Trim adapters at read ends only; never split on interior adapters.
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_ends_only: bool,
    /// Reads sampled for preset detection or ab-initio inference; inference
    /// requires at least 100. Requires an adapter source. Defaults to 0 for
    /// preset detection and 40000 for inference.
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_sample: Option<usize>,
    /// Discover adapters de novo. Report prints the inferred FASTA and exits
    /// without writing read output. Defaults to trim when given no value.
    #[arg(
        long,
        value_enum,
        num_args = 0..=1,
        default_missing_value = "trim",
        help_heading = "Adapter trimming"
    )]
    adapter_infer: Option<AdapterInferActionArg>,
    /// Trust policy for inferred consensuses. Aggressive uses the complete
    /// consensus and permits splitting unless ends-only is set. Defaults to
    /// conservative.
    #[arg(
        long,
        value_enum,
        default_value_t = AdapterInferPolicyArg::Conservative,
        requires = "adapter_infer",
        help_heading = "Adapter trimming"
    )]
    adapter_infer_policy: AdapterInferPolicyArg,
}

/// The default `--adapter-error-rate`.
const DEFAULT_ADAPTER_ERROR_RATE: f64 = 0.2;
/// The default `--adapter-end-size`.
const DEFAULT_ADAPTER_END_SIZE: usize = 150;
/// The default `--adapter-sample` under `--adapter-infer`.
const DEFAULT_INFER_SAMPLE: usize = 40_000;

/// Returns the clap `Command` for the CLI. `examples/gen-man.rs` renders the
/// man page from it, so the page and the parser share one definition.
pub fn command() -> clap::Command {
    <Cli as clap::CommandFactory>::command()
}

#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
enum ProgressArg {
    /// A bar on a terminal, periodic lines otherwise.
    Auto,
    /// The animated bar, unless -v or WHITTLE_LOG asks for log lines.
    Bar,
    /// Always periodic lines, never a bar.
    Plain,
    /// No progress reporting; the banner and summary still print.
    None,
}

impl From<ProgressArg> for ProgressMode {
    fn from(value: ProgressArg) -> Self {
        match value {
            ProgressArg::Auto => ProgressMode::Auto,
            ProgressArg::Bar => ProgressMode::Bar,
            ProgressArg::Plain => ProgressMode::Plain,
            ProgressArg::None => ProgressMode::None,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
enum AdapterInferActionArg {
    /// Trim reads with the inferred sequences.
    Trim,
    /// Print inferred FASTA to stdout and do not write read output.
    Report,
}

impl From<AdapterInferActionArg> for AdapterInferAction {
    fn from(value: AdapterInferActionArg) -> Self {
        match value {
            AdapterInferActionArg::Trim => Self::Trim,
            AdapterInferActionArg::Report => Self::Report,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
enum AdapterInferPolicyArg {
    /// Use a short end-facing anchor and disable inferred interior splitting.
    Conservative,
    /// Use the complete recurrent consensus and allow interior splitting.
    Aggressive,
}

impl From<AdapterInferPolicyArg> for AdapterInferPolicy {
    fn from(value: AdapterInferPolicyArg) -> Self {
        match value {
            AdapterInferPolicyArg::Conservative => Self::Conservative,
            AdapterInferPolicyArg::Aggressive => Self::Aggressive,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
enum FormatArg {
    /// Plain FASTQ.
    Fastq,
    /// gzip-compressed FASTQ.
    FastqGz,
    /// BGZF-compressed FASTQ.
    #[value(name = "fastq-bgz", alias = "fastq-bgzf")]
    FastqBgzf,
    /// Unaligned BAM.
    Bam,
}

impl From<FormatArg> for Format {
    fn from(f: FormatArg) -> Self {
        match f {
            FormatArg::Fastq => Format::Fastq,
            FormatArg::FastqGz => Format::FastqGz,
            FormatArg::FastqBgzf => Format::FastqBgzf,
            FormatArg::Bam => Format::Bam,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
enum QualModeArg {
    /// Average error probabilities, then convert the result back to Phred Q.
    Mean,
    /// Take the arithmetic mean of the per-base Phred scores.
    Arithmetic,
    /// Take the median per-base Phred score.
    Median,
}

impl From<QualModeArg> for QualMode {
    fn from(m: QualModeArg) -> Self {
        match m {
            QualModeArg::Mean => QualMode::Mean,
            QualModeArg::Arithmetic => QualMode::Arithmetic,
            QualModeArg::Median => QualMode::Median,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
enum AdapterPresetArg {
    /// Do not load a built-in adapter catalog.
    None,
    /// Load the built-in Oxford Nanopore adapter and barcode catalog.
    Ont,
}

/// Parses the command line into a validated `Config`.
///
/// Diagnostics that are not errors are collected in `Config::advisories`
/// rather than printed: `parse` runs before the log subscriber exists, and
/// `run` emits them once it does. See `Advisory`.
pub fn parse() -> anyhow::Result<Config> {
    let mut c = Cli::parse();
    // `-` is the pipeline spelling of stdin and stdout, so it is never treated
    // as a file name.
    c.input = c.input.filter(|p| p.as_os_str() != "-");
    c.output = c.output.filter(|p| p.as_os_str() != "-");

    if c.verbose > 2 {
        anyhow::bail!("verbosity accepts at most -vv (debug with -v, trace with -vv)");
    }
    validate_filters(&c)?;
    let compression_level = compression_level_for(&c);
    let quality = quality_op_for(&c);
    let fastq_tags = FastqTags::parse(&c.fastq_tags)?;

    let mut advisories: Vec<Advisory> = Vec::new();
    let adapter_infer = resolve_infer(&c, &mut advisories)?;
    let adapters = resolve_adapters(&c, adapter_infer, &mut advisories)?;
    let adapter_sample = resolve_sample(&c, adapter_infer, &mut advisories)?;

    let ncpu = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    // `-t` is parsed as u64 so clap enforces the lower bound; a value beyond
    // usize saturates and is clamped to the CPU count like any other excess.
    let threads_requested = c.threads.map(|n| usize::try_from(n).unwrap_or(usize::MAX));
    let threads = crate::config::resolve_threads(threads_requested, ncpu);
    let threads_clamped = match threads_requested {
        Some(n) if n > ncpu => Some((n, ncpu)),
        _ => None,
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
        trim: TrimPlan {
            head: c.head_crop,
            tail: c.tail_crop,
            quality,
        },
        adapters,
        adapter_infer,
        threads,
        fastq_tags,
        render_workers: 0,
        adapter_sample,
        compression_level,
        update_moves: c.update_moves,
        ordered: c.ordered,
        verbosity: c.verbose,
        quiet: c.quiet,
        threads_clamped,
        summary_json: c.summary_json,
        advisories,
        progress: c.progress.into(),
        adapter_fasta: c.adapter_fasta,
        adapters_configured: None,
    })
}

/// Rejects contradictory or out-of-domain trim and filter settings before the
/// run, which would otherwise keep zero reads and exit successfully.
fn validate_filters(c: &Cli) -> anyhow::Result<()> {
    let n_quality = [
        c.qual_trim.is_some(),
        c.qual_best_segment.is_some(),
        c.qual_split.is_some(),
    ]
    .iter()
    .filter(|&&b| b)
    .count();
    if n_quality > 1 {
        anyhow::bail!("--qual-trim, --qual-best-segment and --qual-split are mutually exclusive");
    }
    let max_length = c.max_length.unwrap_or(usize::MAX);
    if c.min_length > max_length {
        anyhow::bail!(
            "--min-length ({}) must not exceed --max-length ({max_length})",
            c.min_length
        );
    }
    // NaN compares false against everything, so it slips past the ordering
    // check below and disables the filter; an infinite or negative bound is
    // outside the Phred domain.
    if c.min_qual.is_nan() || c.max_qual.is_nan() {
        anyhow::bail!("--min-qual and --max-qual must be numbers (got NaN)");
    }
    for (flag, value) in [("--min-qual", c.min_qual), ("--max-qual", c.max_qual)] {
        if !value.is_finite() || value < 0.0 {
            anyhow::bail!("{flag} ({value}) must be a finite quality of at least 0");
        }
    }
    if c.min_qual > c.max_qual {
        anyhow::bail!(
            "--min-qual ({}) must not exceed --max-qual ({})",
            c.min_qual,
            c.max_qual
        );
    }
    for (flag, value) in [("--min-gc", c.min_gc), ("--max-gc", c.max_gc)] {
        if let Some(g) = value
            && !(0.0..=1.0).contains(&g)
        {
            anyhow::bail!("{flag} ({g}) must be a fraction between 0 and 1");
        }
    }
    if let (Some(a), Some(b)) = (c.min_gc, c.max_gc)
        && a > b
    {
        anyhow::bail!("--min-gc ({a}) must not exceed --max-gc ({b})");
    }
    Ok(())
}

/// Resolves the compression level. An explicit `-c` wins; otherwise gzip FASTQ
/// output uses level 4, where libdeflate runs markedly faster than level 6 for
/// about 2% more output, and BGZF (BAM and `.bgz`) uses level 6, the BGZF
/// writer default.
fn compression_level_for(c: &Cli) -> u8 {
    let out_is_gz = match c.out_format {
        Some(FormatArg::FastqGz) => true,
        Some(_) => false,
        None => c.output.as_deref().and_then(crate::io::from_extension) == Some(Format::FastqGz),
    };
    c.compression_level.unwrap_or(if out_is_gz { 4 } else { 6 })
}

/// Returns the selected quality-trimming operation, if any; `validate_filters`
/// has already rejected a combination.
fn quality_op_for(c: &Cli) -> Option<QualityOp> {
    if let Some(q) = c.qual_trim {
        return Some(QualityOp::TrimQual(q));
    }
    if let Some(q) = c.qual_best_segment {
        return Some(QualityOp::BestSegment(q));
    }
    c.qual_split.map(|cutoff| QualityOp::Split {
        cutoff,
        window: c.qual_split_window,
    })
}

/// Resolves the ab-initio inference mode and checks it against the other
/// adapter sources.
fn resolve_infer(c: &Cli, advisories: &mut Vec<Advisory>) -> anyhow::Result<AdapterInfer> {
    let adapter_infer = c
        .adapter_infer
        .map_or(AdapterInfer::Off, |action| AdapterInfer::Enabled {
            action: action.into(),
            policy: c.adapter_infer_policy.into(),
        });

    // Trim mode excludes an explicit FASTA; report mode allows one so the
    // discoveries can be named against it.
    if matches!(
        adapter_infer,
        AdapterInfer::Enabled {
            action: AdapterInferAction::Trim,
            ..
        }
    ) && c.adapter_fasta.is_some()
    {
        anyhow::bail!(
            "--adapter-infer and --adapter-fasta are mutually exclusive (one discovers \
             the set, the other supplies it); --adapter-infer report --adapter-fasta <file> \
             names discovered adapters against a supplied FASTA"
        );
    }
    // Under inference the preset is not searched for trimming, since inference
    // builds its own set; it is retained only to name discovered adapters.
    if adapter_infer != AdapterInfer::Off && c.adapter_preset != AdapterPresetArg::None {
        advisories.push(Advisory::warn(
            "--adapter-preset is ignored for trimming under --adapter-infer \
             (used only for naming discovered adapters)",
        ));
    }
    // Report mode names discovered adapters against the union of the built-in
    // ONT catalog and the supplied FASTA (see `infer::discover`), so FASTA entry
    // names appear alongside catalog names.
    if adapter_infer.is_report() && c.adapter_fasta.is_some() {
        advisories.push(Advisory::info(
            "--adapter-infer report with --adapter-fasta: discovered adapters are named \
             against the built-in ONT catalog and the supplied FASTA",
        ));
    }
    Ok(adapter_infer)
}

/// Resolves the adapter set and its search settings, or `None` when no
/// adapter source is given.
///
/// The tuning flags are meaningful only with a source, so an explicit one
/// without a source is rejected rather than ignored.
fn resolve_adapters(
    c: &Cli,
    adapter_infer: AdapterInfer,
    advisories: &mut Vec<Advisory>,
) -> anyhow::Result<Option<crate::adapter::AdapterConfig>> {
    let mut adapter_seqs: Vec<crate::adapter::Adapter> = Vec::new();
    if c.adapter_preset == AdapterPresetArg::Ont {
        adapter_seqs.extend(crate::adapter::preset::preset_ont());
    }
    // Only the FASTA entries are carried onward as naming references under
    // inference; `infer::discover` looks up the built-in catalog itself.
    let mut fasta_adapters: Vec<crate::adapter::Adapter> = Vec::new();
    if let Some(path) = &c.adapter_fasta {
        let from_fasta = read_adapter_fasta(path, advisories)?;
        if from_fasta.is_empty() {
            anyhow::bail!(
                "--adapter-fasta {}: no usable adapters (all entries were empty, \
                 shorter than the {}-bp minimum, or non-nucleotide)",
                path.display(),
                crate::adapter::MIN_PATTERN_LEN
            );
        }
        fasta_adapters = from_fasta.clone();
        adapter_seqs.extend(from_fasta);
    }

    if adapter_seqs.is_empty() && adapter_infer == AdapterInfer::Off {
        require_adapter_source(c)?;
        if c.adapter_ends_only {
            advisories.push(Advisory::warn(
                "--adapter-ends-only has no effect without --adapter-fasta or --adapter-preset",
            ));
        }
        return Ok(None);
    }

    let error_rate = c.adapter_error_rate.unwrap_or(DEFAULT_ADAPTER_ERROR_RATE);
    if !(0.0..=1.0).contains(&error_rate) {
        anyhow::bail!("--adapter-error-rate ({error_rate}) must be between 0 and 1");
    }
    let end_size = c.adapter_end_size.unwrap_or(DEFAULT_ADAPTER_END_SIZE);
    if end_size == 0 {
        anyhow::bail!("--adapter-end-size must be >= 1");
    }
    // Under inference the trimming set is discovered later, so the preset
    // sequences are dropped here. A report-only FASTA is carried in this field
    // only as naming references for `infer::discover`: discovery replaces the
    // field before dispatch and report mode exits first, so the FASTA is never
    // trimmed against. Under `Trim` a FASTA is rejected by `resolve_infer`.
    let trim_adapters = if adapter_infer == AdapterInfer::Off {
        adapter_seqs
    } else {
        fasta_adapters
    };
    let infer_forces_ends_only =
        adapter_infer != AdapterInfer::Off && !adapter_infer.is_aggressive();
    if infer_forces_ends_only && !c.adapter_ends_only {
        advisories.push(Advisory::info(
            "Conservative adapter inference trims read ends only; use \
             --adapter-infer-policy aggressive to enable full-consensus interior splitting",
        ));
    }
    Ok(Some(crate::adapter::AdapterConfig {
        adapters: trim_adapters,
        error_rate,
        end_size,
        split: !c.adapter_ends_only && !infer_forces_ends_only,
        candidate_index: std::sync::OnceLock::new(),
    }))
}

/// Rejects an explicit adapter tuning flag given without an adapter source.
fn require_adapter_source(c: &Cli) -> anyhow::Result<()> {
    let explicit = [
        ("--adapter-error-rate", c.adapter_error_rate.is_some()),
        ("--adapter-end-size", c.adapter_end_size.is_some()),
        ("--adapter-sample", c.adapter_sample.is_some()),
    ];
    if let Some((flag, _)) = explicit.iter().find(|(_, given)| *given) {
        anyhow::bail!(
            "{flag} requires an adapter source (--adapter-fasta, --adapter-preset ont, or \
             --adapter-infer)"
        );
    }
    Ok(())
}

/// Resolves the sample size for presence detection or inference.
///
/// An omitted value means the mode default: 0 (detection off) with inference
/// off, 40000 with inference on. An explicit value must be 0 or at least
/// `MIN_SAMPLE_FOR_DETECTION`, and 0 is rejected under inference, which needs a
/// sample. Presence detection is preset-only: a user-supplied FASTA is a curated
/// set that is searched in full, since sampling could drop a rare custom
/// adapter, so detection is disabled whenever a FASTA is given and inference is
/// off.
fn resolve_sample(
    c: &Cli,
    adapter_infer: AdapterInfer,
    advisories: &mut Vec<Advisory>,
) -> anyhow::Result<usize> {
    let min = crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION;
    let requested = match c.adapter_sample {
        None if adapter_infer != AdapterInfer::Off => DEFAULT_INFER_SAMPLE,
        None => 0,
        Some(n) => {
            if n != 0 && n < min {
                anyhow::bail!(
                    "--adapter-sample ({n}) must be 0 (disable detection) or at least {min} \
                     (smaller samples are too few for reliable detection)"
                );
            }
            if n == 0 && adapter_infer != AdapterInfer::Off {
                anyhow::bail!(
                    "--adapter-sample 0 disables sampling, which --adapter-infer requires; \
                     omit it or pass >= {min}"
                );
            }
            n
        },
    };
    if adapter_infer != AdapterInfer::Off || c.adapter_fasta.is_none() {
        return Ok(requested);
    }
    if requested > 0 {
        advisories.push(Advisory::warn(
            "--adapter-sample is ignored with --adapter-fasta (presence detection is \
             preset-only)",
        ));
    }
    Ok(0)
}

/// Reads adapter sequences from a FASTA. Whitespace is removed, lowercase is
/// uppercased, and `U` is folded to `T`. IUPAC ambiguity codes are kept and
/// searched as the bases they stand for; an entry containing any other byte is
/// skipped with a warning advisory, as is an entry shorter than
/// `adapter::MIN_PATTERN_LEN`, the matcher's minimum pattern length. An entry
/// averaging two or more bases per position is kept with a warning advisory.
fn read_adapter_fasta(
    path: &std::path::Path,
    advisories: &mut Vec<crate::config::Advisory>,
) -> anyhow::Result<Vec<crate::adapter::Adapter>> {
    use seq_io::fasta::{Reader, Record};
    let mut reader = Reader::from_path(path)
        .map_err(|e| anyhow::anyhow!("--adapter-fasta {}: {e}", path.display()))?;
    let mut out = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| anyhow::anyhow!("--adapter-fasta {}: {e}", path.display()))?;
        let seq: Vec<u8> = rec
            .seq()
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(u8::to_ascii_uppercase)
            .collect();
        // `U` is folded to `T`: RNA primers are written with `U`, DNA reads
        // store `T`, and sassy treats `U` as a fifth base that matches nothing.
        let seq: Vec<u8> = seq
            .into_iter()
            .map(|b| if b == b'U' { b'T' } else { b })
            .collect();
        let name = String::from_utf8_lossy(rec.head()).into_owned();

        // IUPAC ambiguity codes are searched as the bases they stand for, as a
        // degenerate primer requires. A byte outside the nucleotide alphabet
        // marks a malformed record.
        let Some(degeneracy) = seq
            .iter()
            .map(|&b| crate::adapter::search::iupac_degeneracy(b).map(u32::from))
            .sum::<Option<u32>>()
        else {
            let bad: String = seq
                .iter()
                .filter(|&&b| crate::adapter::search::iupac_degeneracy(b).is_none())
                .map(|&b| b as char)
                .collect();
            advisories.push(crate::config::Advisory::warn(format!(
                "Adapter entry skipped, contains non-nucleotide characters: \
                 name={name:?}, chars={bad}"
            )));
            continue;
        };
        if seq.len() < crate::adapter::MIN_PATTERN_LEN {
            advisories.push(crate::config::Advisory::warn(format!(
                "Adapter entry skipped, shorter than the minimum match length: \
                 name={name:?}, len={}, min={}",
                seq.len(),
                crate::adapter::MIN_PATTERN_LEN
            )));
            continue;
        }
        // A fully degenerate stretch matches anywhere, so a pattern with more
        // ambiguity than specificity trims real insert. The pattern is still
        // searched; a warning advisory records the risk.
        if degeneracy as usize >= seq.len() * 2 {
            advisories.push(crate::config::Advisory::warn(format!(
                "Adapter entry is highly degenerate and matches almost anywhere, which may \
                 trim real sequence: name={name:?}, bases_per_position={:.1}",
                f64::from(degeneracy) / seq.len() as f64
            )));
        }

        out.push(crate::adapter::Adapter {
            name,
            seq,
            end: crate::adapter::End::Both,
        });
    }
    Ok(out)
}

/// Builds a `Config` directly for integration tests. `head_crop` and
/// `tail_crop` are fixed crops.
#[doc(hidden)]
pub fn config_for_test(
    input: &std::path::Path,
    output: &std::path::Path,
    head_crop: usize,
    tail_crop: usize,
) -> Config {
    config_for_test_threads(input, output, head_crop, tail_crop, 1)
}

/// Builds a `Config` as `config_for_test` does, with an explicit thread count
/// for tests that exercise the parallel BAM dispatch.
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
        trim: TrimPlan {
            head: head_crop,
            tail: tail_crop,
            quality: None,
        },
        adapters: None,
        adapter_infer: crate::config::AdapterInfer::Off,
        threads: threads.max(1),
        fastq_tags: FastqTags::All,
        render_workers: 0,
        adapter_sample: 0,
        compression_level: 6,
        update_moves: false,
        ordered: false,
        verbosity: 0,
        quiet: true,
        threads_clamped: None,
        summary_json: None,
        advisories: Vec::new(),
        adapter_fasta: None,
        progress: crate::config::ProgressMode::Auto,
        adapters_configured: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Help text quotes the matcher's minimum pattern length and the detection
    /// sample floor; the numbers stay in step with the constants.
    #[test]
    fn help_text_matches_the_adapter_constants() {
        let cmd = command();
        let help_for = |id: &str| {
            cmd.get_arguments()
                .find(|a| a.get_id().as_str() == id)
                .unwrap_or_else(|| panic!("Argument {id} exists"))
                .get_help()
                .expect("Argument has help text")
                .to_string()
        };
        assert!(
            help_for("adapter_fasta").contains(&format!("{} bp", crate::adapter::MIN_PATTERN_LEN)),
            "--adapter-fasta help quotes MIN_PATTERN_LEN"
        );
        assert!(
            help_for("adapter_sample").contains(&format!(
                "at least {}",
                crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION
            )),
            "--adapter-sample help quotes MIN_SAMPLE_FOR_DETECTION"
        );
    }

    /// FASTA loading and adapter search enforce the same minimum pattern length.
    #[test]
    fn read_adapter_fasta_skips_entries_below_min_pattern_len() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("adapters.fasta");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, ">kept_20bp").unwrap();
        writeln!(f, "ACGTACGTACGTACGTACGT").unwrap(); // 20 bp
        writeln!(f, ">skipped_8bp").unwrap();
        writeln!(f, "ACGTACGT").unwrap(); // 8 bp, below the 11-bp `MIN_PATTERN_LEN`
        drop(f);

        let mut advisories = Vec::new();
        let adapters = read_adapter_fasta(&path, &mut advisories).unwrap();

        assert_eq!(adapters.len(), 1);
        assert_eq!(adapters[0].seq, b"ACGTACGTACGTACGTACGT".to_vec());
    }

    /// IUPAC ambiguity codes are kept and searched as the bases they stand for,
    /// since degenerate primers are written with them. Only characters outside
    /// the nucleotide alphabet mark a malformed record.
    #[test]
    fn read_adapter_fasta_keeps_iupac_and_rejects_non_nucleotides() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("adapters.fasta");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, ">plain_20bp").unwrap();
        writeln!(f, "ACGTACGTACGTACGTACGT").unwrap();
        writeln!(f, ">degenerate_20bp").unwrap();
        writeln!(f, "ACGTACGTYCGTACGTACGN").unwrap();
        writeln!(f, ">rna_20bp").unwrap();
        writeln!(f, "ACGUACGUACGUACGUACGU").unwrap();
        writeln!(f, ">protein_20bp").unwrap();
        writeln!(f, "ACGTACGTZCGTACGTACGT").unwrap();
        drop(f);

        let mut advisories = Vec::new();
        let adapters = read_adapter_fasta(&path, &mut advisories).unwrap();

        let names: Vec<&str> = adapters.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, ["plain_20bp", "degenerate_20bp", "rna_20bp"]);
        assert_eq!(adapters[1].seq, b"ACGTACGTYCGTACGTACGN".to_vec());
        // `U` folds to `T`: a DNA read stores `T`, and sassy treats `U` as a
        // fifth base that matches nothing.
        assert_eq!(adapters[2].seq, b"ACGTACGTACGTACGTACGT".to_vec());
        assert!(
            advisories
                .iter()
                .any(|a| a.message.contains("non-nucleotide")),
            "The protein-alphabet entry is reported"
        );
    }

    /// A pattern with more ambiguity than specificity matches almost anywhere.
    /// It is still searched, and a warning advisory is recorded.
    #[test]
    fn read_adapter_fasta_warns_about_a_very_degenerate_pattern() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("adapters.fasta");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, ">mostly_n").unwrap();
        writeln!(f, "NNNNNNNNNNNNNNNACGTA").unwrap();
        drop(f);

        let mut advisories = Vec::new();
        let adapters = read_adapter_fasta(&path, &mut advisories).unwrap();
        assert_eq!(adapters.len(), 1, "Still searched");
        assert!(
            advisories.iter().any(|a| a.message.contains("degenerate")),
            "A near-wildcard pattern warns"
        );
    }
}
