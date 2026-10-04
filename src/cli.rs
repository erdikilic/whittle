//! Command-line parsing: the clap definition, cross-flag validation, and
//! construction of the resolved `Config`.

use std::path::PathBuf;

use clap::Parser;

use crate::config::{
    AdapterInfer, Advisory, Config, FastqTags, IoConfig, ProgressMode, TagRemoval,
};
use crate::filter::FilterConfig;
use crate::io::Format;
use crate::qual::QualMode;
use crate::split::{Action, KeyLevel, Require, Rules, Sheet, SplitOptions, Template};
use crate::trim::{QualityMethod, QualityOp, TrimPlan};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    disable_version_flag = true,
    about = "Coordinate-consistent trimming of long-read uBAM and tagged FASTQ",
    long_about = None,
    max_term_width = 100,
    after_help = EXAMPLES
)]
struct Cli {
    /// Print version information and exit.
    #[arg(long, action = clap::ArgAction::Version, help_heading = "Setup")]
    version: Option<bool>,
    /// Input FASTQ-family file, unaligned BAM, or directory; - means stdin.
    /// Defaults to stdin.
    #[arg(short = 'i', long, value_name = "PATH", help_heading = "Setup")]
    input: Option<PathBuf>,
    /// Output file, whose extension selects the format; - means stdout.
    /// Defaults to stdout. Under --split-by, a path holding {target},
    /// {group} or {barcode} is a template with one file per split key.
    #[arg(short = 'o', long, value_name = "PATH", help_heading = "Setup")]
    output: Option<PathBuf>,
    /// Force the input format instead of detecting it from the path or stream.
    #[arg(
        long = "input-format",
        value_enum,
        value_name = "FORMAT",
        help_heading = "Setup"
    )]
    in_format: Option<Format>,
    /// Force the output format instead of selecting it from the output path.
    #[arg(
        long = "output-format",
        value_enum,
        value_name = "FORMAT",
        help_heading = "Setup"
    )]
    out_format: Option<Format>,
    /// Worker threads, at least 1; values above the CPU count are clamped to
    /// it. Defaults to all detected CPUs.
    #[arg(
        short = 't',
        long,
        value_name = "N",
        value_parser = parse_threads,
        help_heading = "Setup"
    )]
    threads: Option<u64>,
    /// Write records in input order when running with more than one thread.
    /// Without it, records are written as they finish, which is faster and uses
    /// less memory but is not reproducible between runs.
    #[arg(long = "preserve-order", help_heading = "Setup")]
    ordered: bool,
    /// Aux tags written into FASTQ headers on BAM or tagged FASTQ input: all,
    /// none, or a list such as MM,ML,RG.
    #[arg(
        long,
        default_value = "all",
        value_name = "all|none|TAGS",
        help_heading = "Tags"
    )]
    fastq_tags: String,
    /// BGZF compression level (0-9) for BAM, .bgz and .gz output. Lower levels
    /// are faster and produce larger files. Ignored for plain FASTQ. Defaults
    /// to 4 for .gz and 6 for .bgz and BAM.
    // libdeflate accepts up to 12; the cap is the conventional gzip 0-9.
    #[arg(
        short = 'c',
        long,
        value_name = "0-9",
        value_parser = clap::value_parser!(u8).range(0..=9),
        help_heading = "Setup"
    )]
    compression_level: Option<u8>,
    /// Write a machine-readable JSON run summary (counters plus the resolved
    /// settings) to this path. Written even under --quiet.
    #[arg(long, value_name = "PATH", help_heading = "Setup")]
    summary_json: Option<PathBuf>,
    /// Write every input read or trimmed segment that does not reach the
    /// output to this path, in the output's format family, with a wr:Z tag
    /// naming the reason: tag_filter, trimmed_to_nothing, too_short,
    /// too_long, low_quality, high_quality, expected_errors or gc.
    #[arg(long = "rejected-output", value_name = "PATH", help_heading = "Setup")]
    rejected_output: Option<PathBuf>,

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
    /// back to lines under -v or WHITTLE_LOG.
    #[arg(
        long,
        value_enum,
        value_name = "MODE",
        default_value_t = ProgressMode::Auto,
        conflicts_with = "quiet",
        help_heading = "Logging"
    )]
    progress: ProgressMode,

    /// Minimum post-trim segment length.
    #[arg(
        short = 'l',
        long,
        value_name = "BASES",
        default_value_t = 1,
        help_heading = "Filtering"
    )]
    min_length: usize,
    /// Maximum post-trim segment length.
    #[arg(short = 'L', long, value_name = "BASES", help_heading = "Filtering")]
    max_length: Option<usize>,
    /// Minimum post-trim segment quality under --quality-mode.
    #[arg(
        short = 'q',
        long = "min-quality",
        value_name = "PHRED",
        default_value_t = 0.0,
        help_heading = "Filtering"
    )]
    min_qual: f64,
    /// Maximum post-trim segment quality under --quality-mode.
    #[arg(
        short = 'Q',
        long = "max-quality",
        value_name = "PHRED",
        default_value_t = 1000.0,
        help_heading = "Filtering"
    )]
    max_qual: f64,
    /// Maximum expected errors per output segment: the sum of its per-base
    /// error probabilities 10^(-Q/10) after trimming (Edgar and Flyvbjerg
    /// 2015). A total per segment rather than a rate, for amplicons of a fixed
    /// length; a non-negative number.
    #[arg(
        long = "max-expected-errors",
        value_name = "E",
        help_heading = "Filtering"
    )]
    max_expected_errors: Option<f64>,
    /// Minimum post-trim GC fraction (0 to 1; 0.4 means 40%).
    #[arg(short = 'g', long, value_name = "FRACTION", help_heading = "Filtering")]
    min_gc: Option<f64>,
    /// Keep only reads whose aux tags satisfy EXPR, a samtools-style
    /// expression over [tag] values, for example '[er]!="unblock_mux_change"'
    /// or '[qs]>=10 && [dx]==1'. Repeatable; every expression must hold.
    /// Applied before adapter discovery and trimming. BAM or tagged FASTQ
    /// input.
    #[arg(long = "tag-filter", value_name = "EXPR", help_heading = "Filtering")]
    tag_filter: Vec<String>,
    /// Maximum post-trim GC fraction (0 to 1; 0.4 means 40%).
    #[arg(short = 'G', long, value_name = "FRACTION", help_heading = "Filtering")]
    max_gc: Option<f64>,
    /// Quality calculation of the mean-quality filter, --min-quality and
    /// --max-quality, on each output segment. --quality-trim and
    /// --max-expected-errors use per-base qualities and are unaffected; a
    /// recomputed dorado qs tag is always the error-probability mean.
    #[arg(
        short = 'm',
        long = "quality-mode",
        value_enum,
        value_name = "MODE",
        default_value_t = QualMode::Mean,
        help_heading = "Filtering"
    )]
    qual_mode: QualMode,

    /// Remove this many bases from each adapter-derived segment's 5' end,
    /// after barcode restriction and before quality processing. Applied once;
    /// quality-split pieces are not cropped again.
    #[arg(
        short = 'H',
        long = "trim-front",
        visible_alias = "head-crop",
        value_name = "BASES",
        default_value_t = 0,
        help_heading = "Trimming"
    )]
    head_crop: usize,
    /// Remove this many bases from each adapter-derived segment's 3' end,
    /// after barcode restriction and before quality processing. Applied once;
    /// quality-split pieces are not cropped again.
    #[arg(
        short = 'T',
        long = "trim-tail",
        visible_alias = "tail-crop",
        value_name = "BASES",
        default_value_t = 0,
        help_heading = "Trimming"
    )]
    tail_crop: usize,
    /// Quality trimming of each adapter-derived segment after the fixed crop,
    /// scored against --quality-cutoff. The segments and runs methods can split
    /// a segment into pieces, numbered in original-read order.
    #[arg(
        long = "quality-trim",
        value_enum,
        value_name = "METHOD",
        help_heading = "Trimming"
    )]
    quality_trim: Option<QualityMethod>,
    /// Phred cutoff of --quality-trim. Required with it.
    #[arg(
        long = "quality-cutoff",
        value_name = "PHRED",
        help_heading = "Trimming"
    )]
    quality_cutoff: Option<u8>,
    /// Minimum run of consecutive bases below --quality-cutoff that splits a
    /// read under --quality-trim runs; shorter runs stay inside their piece.
    /// Defaults to 1.
    #[arg(
        long = "min-low-quality-run",
        value_name = "BASES",
        help_heading = "Trimming"
    )]
    min_low_quality_run: Option<usize>,
    /// Trim the ends of every piece of --quality-trim segments or runs, split
    /// or not, up to the first base at or above PHRED, so a read splits only at
    /// long weak interior regions while its ends are trimmed more strictly.
    /// Defaults to --quality-cutoff.
    #[arg(
        long = "quality-end-cutoff",
        value_name = "PHRED",
        help_heading = "Trimming"
    )]
    quality_end_cutoff: Option<u8>,
    /// Keep the ONT move table consistent through trimming (slice mv, update
    /// ts, ns, sp, pi and the poly-A tags) for signal-aware tools such as
    /// Remora and Clair3 v2, instead of dropping them. BAM-to-BAM only;
    /// requires a DNA or RNA basecall_model in the read-group description.
    #[arg(long = "update-moves", help_heading = "Tags")]
    update_moves: bool,
    /// Remove aux tags from every output record. Comma-separated and
    /// repeatable; an item is a two-character tag or a group: kinetics (ip pw
    /// fi fp ri rp sa sm sx), mods (MM ML MN), signal (mv ts ns sp pi). BAM or
    /// tagged FASTQ input.
    #[arg(long, value_name = "TAGS", help_heading = "Tags")]
    remove_tag: Vec<String>,

    /// Adapter FASTA; sequences may use IUPAC codes, and entries shorter than
    /// 11 bp are skipped. An entry whose header description contains the word
    /// primer or barcode is trimmed at read ends only; every other entry also
    /// splits reads at interior hits. Enables adapter trimming.
    #[arg(
        short = 'a',
        long,
        value_name = "FILE",
        help_heading = "Adapter trimming"
    )]
    adapter_fasta: Option<PathBuf>,
    /// Built-in kit presets, comma-separated: lsk114, rad114 (ulk114), rbk114,
    /// nbd114, pcb114 (pcs114), rpb114, mab114, rna004, pacbio, ont (every
    /// ONT kit), all. Enables adapter trimming. Defaults to none.
    #[arg(long, value_name = "KITS", help_heading = "Adapter trimming")]
    adapter_preset: Option<String>,
    /// End-match tolerance as a fraction of adapter length. Interior tolerance
    /// depends on the adapter sequence and read length and never exceeds this
    /// value. Requires an adapter source. Defaults to 0.2.
    #[arg(long, value_name = "FRACTION", help_heading = "Adapter trimming")]
    adapter_error_rate: Option<f64>,
    /// Bases at each read end searched for a terminal adapter. Requires an
    /// adapter source. Defaults to 150.
    #[arg(
        long = "adapter-end-search",
        value_name = "BASES",
        help_heading = "Adapter trimming"
    )]
    adapter_end_size: Option<usize>,
    /// Trim adapters at read ends only; never split on interior adapters.
    /// Does not disable quality splitting.
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_ends_only: bool,
    /// Reads inspected for preset presence or adapter discovery; at least 100.
    /// 0 disables presence detection and uses the full preset; discovery
    /// requires sampling. Ignored with --adapter-fasta unless discovering
    /// adapters. Defaults to 2000 with a preset and 40000 under
    /// --adapter-discover or --adapter-report.
    /// Sampling is also bounded by 256 MiB of payload and 64 Mi bases.
    #[arg(
        long = "adapter-sample-reads",
        value_name = "COUNT",
        help_heading = "Adapter trimming"
    )]
    adapter_sample: Option<usize>,
    /// Discover adapters, barcodes and primers de novo in layers from each
    /// read end and trim them. A preset or FASTA is trimmed first and
    /// discovery continues beyond it. Enables adapter trimming.
    #[arg(
        long = "adapter-discover",
        help_heading = "Adapter trimming",
        conflicts_with = "adapter_report"
    )]
    adapter_discover: bool,
    /// Run the same discovery as --adapter-discover, print the discovered
    /// FASTA to stdout and exit without read output or a JSON summary.
    #[arg(long = "adapter-report", help_heading = "Adapter trimming")]
    adapter_report: bool,

    /// Assign every output read or segment to a primer target from the
    /// primers at its ends and tag it wt:Z with the target, unassigned or
    /// ambiguous. Repeatable; each SPEC adds a source, merged into one
    /// sheet. SPEC is an existing sheet file (TSV or annotated FASTA), the
    /// preset mab114 (also enables --adapter-preset mab114 unless -a or
    /// --adapter-preset is given), or an inline target
    /// NAME:F:SEQ[,SEQ...]:R:SEQ[,SEQ...]. A target holds a primer mix:
    /// any of its forward primers pairs with any of its reverse primers.
    /// Further pairs after a comma (,F:...:R:...) form a pool of targets
    /// NAME.1, NAME.2, ... in group NAME. Split primers are
    /// trimmed with the adapters. An -o template such as
    /// out/{barcode}.{target}.fastq.gz writes each key to its own file:
    /// {target} or {group} names the key (unassigned and ambiguous
    /// included), {barcode} the BC:Z barcode call (unclassified without
    /// one); without a placeholder every read goes to the one output.
    #[arg(long = "split-by", value_name = "SPEC", help_heading = "Split")]
    split_by: Vec<String>,
    /// Which read ends need a located primer for an assignment.
    #[arg(
        long = "split-require",
        value_enum,
        value_name = "RULE",
        default_value_t = Require::Either,
        requires = "split_by",
        help_heading = "Split"
    )]
    split_require: Require,
    /// Whether located split primers are trimmed or kept in the output.
    #[arg(
        long = "split-action",
        value_enum,
        value_name = "ACTION",
        default_value_t = Action::Trim,
        requires = "split_by",
        help_heading = "Split"
    )]
    split_action: Action,
    /// Minimum edit-cost lead of the best target over the best different
    /// target; a smaller lead is ambiguous.
    #[arg(
        long = "split-lead",
        value_name = "N",
        default_value_t = 2,
        requires = "split_by",
        help_heading = "Split"
    )]
    split_lead: usize,
    /// Drop these bins instead of writing them, comma-separated: unassigned,
    /// ambiguous.
    #[arg(
        long = "split-discard",
        value_enum,
        value_name = "BINS",
        value_delimiter = ',',
        requires = "split_by",
        help_heading = "Split"
    )]
    split_discard: Vec<SplitBin>,
}

/// A `--split-discard` bin.
#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
enum SplitBin {
    /// Reads and segments no target is assigned to.
    Unassigned,
    /// Reads and segments with two targets within --split-lead.
    Ambiguous,
}

/// The examples block at the end of `--help`.
const EXAMPLES: &str = "\
Examples:
  whittle -i reads.fastq.gz -o trimmed.fastq.gz -H 20 --quality-trim ends --quality-cutoff 8 -l 500
  whittle -i reads.bam -o trimmed.bam --quality-trim segments --quality-cutoff 9 -l 1000
  whittle -i reads.bam -o trimmed.bam --adapter-preset lsk114 -l 500
  whittle -i reads.bam -o trimmed.bam --adapter-discover -l 500
  whittle -i 16s.fastq.gz -o trimmed.fastq.gz --adapter-preset mab114
  samtools fastq -T MM,ML,MN reads.bam | whittle -o trimmed.fastq.gz -H 10 -T 10
  whittle -i reads.bam -o reads.fastq.gz --quiet --summary-json qc.json";

/// The default `--adapter-error-rate`.
const DEFAULT_ADAPTER_ERROR_RATE: f64 = 0.2;
/// The default `--adapter-end-search`.
const DEFAULT_ADAPTER_END_SIZE: usize = 150;
/// The default `--adapter-sample-reads` under `--adapter-discover` and `--adapter-report`.
const DEFAULT_INFER_SAMPLE: usize = 40_000;
/// The default `--adapter-sample-reads` for preset presence detection.
const DEFAULT_PRESET_SAMPLE: usize = 2_000;

/// Returns the clap `Command` for the CLI. `examples/gen-man.rs` renders the
/// man page from it, so the page and the parser share one definition.
pub fn command() -> clap::Command {
    <Cli as clap::CommandFactory>::command()
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
    let fastq_tags = FastqTags::parse(&c.fastq_tags)?;
    let remove_tags = TagRemoval::parse(&c.remove_tag)?;
    let tag_filters = crate::tagfilter::TagFilters::parse(&c.tag_filter)?;
    let quality = resolve_quality(&c)?;

    let mut advisories: Vec<Advisory> = Vec::new();
    let split_opts = resolve_split(&mut c)?;
    let adapter_infer = resolve_infer(&c, &mut advisories)?;
    let adapters = resolve_adapters(&c, adapter_infer, &mut advisories)?;
    let adapter_sample = resolve_sample(&c, adapter_infer, &mut advisories)?;
    let split_search = match split_opts {
        Some(_) => Some(search_config(&c, Vec::new())?),
        None => None,
    };

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

    let cfg = Config {
        io: IoConfig {
            input: c.input,
            output: c.output,
            in_format: c.in_format,
            out_format: c.out_format,
        },
        filter: FilterConfig {
            min_length: c.min_length,
            max_length: c.max_length.unwrap_or(usize::MAX),
            min_qual: c.min_qual,
            max_qual: c.max_qual,
            max_expected_errors: c.max_expected_errors,
            min_gc: c.min_gc,
            max_gc: c.max_gc,
            qual_mode: c.qual_mode,
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
        rejected_output: c.rejected_output,
        advisories,
        progress: c.progress,
        adapter_fasta: c.adapter_fasta,
        adapters_configured: None,
        trim_classes: [false; 3],
        remove_tags,
        tag_filters,
        split_opts,
        split_search,
        splitter: None,
    };

    // Only an explicit `--input-format` or a known extension decides the input
    // format here; a stream or an extensionless path is classified by `run`,
    // which applies the same guard once detection has run.
    let in_fmt = cfg
        .io
        .in_format
        .or_else(|| cfg.io.input.as_deref().and_then(crate::io::from_extension));
    if let Some(fmt) = in_fmt {
        crate::guards::guard_bam_only_flags(&cfg, fmt)?;
    }
    let out_fmt = cfg
        .io
        .out_format
        .or_else(|| cfg.io.output.as_deref().and_then(crate::io::from_extension));
    if let (Some(in_fmt), Some(out_fmt)) = (in_fmt, out_fmt) {
        crate::guards::guard_fastq_to_bam(in_fmt, out_fmt)?;
    }
    Ok(cfg)
}

/// Parses `-t`: an integer of at least 1.
fn parse_threads(value: &str) -> Result<u64, String> {
    match value.parse::<u64>() {
        Ok(0) => Err("must be at least 1".to_string()),
        Ok(n) => Ok(n),
        Err(_) => Err("must be an integer of at least 1".to_string()),
    }
}

/// Rejects contradictory or out-of-domain trim and filter settings before the
/// run, which would otherwise keep zero reads and exit successfully.
fn validate_filters(c: &Cli) -> anyhow::Result<()> {
    let max_length = c.max_length.unwrap_or(usize::MAX);
    if c.min_length > max_length {
        anyhow::bail!(
            "--min-length ({}) must not exceed --max-length ({max_length})",
            c.min_length
        );
    }
    // NaN compares false against everything, so it would slip past the ordering
    // check below and disable the filter; `is_finite` rejects it together with
    // the infinite bounds, which are outside the Phred domain.
    for (flag, value) in [("--min-quality", c.min_qual), ("--max-quality", c.max_qual)] {
        if !value.is_finite() || value < 0.0 {
            anyhow::bail!("{flag} ({value}) must be a finite quality of at least 0");
        }
    }
    if let Some(e) = c.max_expected_errors
        && (!e.is_finite() || e < 0.0)
    {
        anyhow::bail!("--max-expected-errors ({e}) must be a finite number of at least 0");
    }
    if c.min_qual > c.max_qual {
        anyhow::bail!(
            "--min-quality ({}) must not exceed --max-quality ({})",
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
/// output uses level 4, which libdeflate compresses faster than level 6 for
/// about 2% more output, and BGZF (BAM and `.bgz`) uses level 6, the BGZF
/// writer default.
fn compression_level_for(c: &Cli) -> u8 {
    let out_is_gz = match c.out_format {
        Some(Format::FastqGz) => true,
        Some(_) => false,
        None => c.output.as_deref().and_then(crate::io::from_extension) == Some(Format::FastqGz),
    };
    c.compression_level.unwrap_or(if out_is_gz { 4 } else { 6 })
}

/// Resolves `--quality-trim`, `--quality-cutoff`, `--min-low-quality-run` and
/// `--quality-end-cutoff` into the quality operation. The cutoff is required
/// with a method, since no single value suits every platform and basecaller.
fn resolve_quality(c: &Cli) -> anyhow::Result<Option<QualityOp>> {
    if let Some(n) = c.min_low_quality_run {
        if c.quality_trim != Some(QualityMethod::Runs) {
            anyhow::bail!("--min-low-quality-run requires --quality-trim runs");
        }
        if n == 0 {
            anyhow::bail!("--min-low-quality-run must be at least 1");
        }
    }
    if c.quality_end_cutoff.is_some()
        && !matches!(
            c.quality_trim,
            Some(QualityMethod::Segments | QualityMethod::Runs)
        )
    {
        anyhow::bail!("--quality-end-cutoff requires --quality-trim segments or runs");
    }
    let Some(method) = c.quality_trim else {
        if c.quality_cutoff.is_some() {
            anyhow::bail!("--quality-cutoff requires --quality-trim");
        }
        return Ok(None);
    };
    let Some(cutoff) = c.quality_cutoff else {
        anyhow::bail!("--quality-trim requires --quality-cutoff");
    };
    Ok(Some(QualityOp {
        min_low_quality_run: c.min_low_quality_run.unwrap_or(1),
        end_cutoff: c.quality_end_cutoff.unwrap_or(cutoff),
        ..QualityOp::new(method, cutoff)
    }))
}

/// Resolves the de novo inference mode and checks it against the other
/// adapter sources.
fn resolve_infer(c: &Cli, advisories: &mut Vec<Advisory>) -> anyhow::Result<AdapterInfer> {
    let adapter_infer = if c.adapter_report {
        AdapterInfer::Report
    } else if c.adapter_discover {
        AdapterInfer::Discover
    } else {
        AdapterInfer::Off
    };

    // A FASTA or preset supplied with inference is trimmed first, and
    // inference continues from the boundary those sequences leave. Report
    // mode names discovered adapters against the built-in catalog and the
    // supplied FASTA (see `infer::discover`).
    if adapter_infer.is_report() && c.adapter_fasta.is_some() {
        advisories.push(Advisory::info(
            "--adapter-report with --adapter-fasta: discovered adapters are named \
             against the built-in adapter catalog and the supplied FASTA",
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
    // Preset and FASTA sequences form the known set. Under inference they are
    // trimmed first and also name the discoveries; `infer::discover` looks up
    // the built-in catalog itself. Report mode exits before any trimming.
    let mut adapter_seqs: Vec<crate::adapter::Adapter> = Vec::new();
    if let Some(kits) = preset_kits(c)? {
        adapter_seqs.extend(crate::adapter::preset::preset(&kits));
    }
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
        adapter_seqs.extend(from_fasta);
    }

    // `--split-by` searches its primers with the adapter settings, so the
    // tuning flags apply to it without another source.
    if adapter_seqs.is_empty() && adapter_infer == AdapterInfer::Off {
        if c.split_by.is_empty() {
            require_adapter_source(c)?;
            if c.adapter_ends_only {
                advisories.push(Advisory::warn(
                    "--adapter-ends-only has no effect without --adapter-fasta, --adapter-preset or \
                     --adapter-discover",
                ));
            }
        } else if c.adapter_sample.is_some() {
            advisories.push(Advisory::info(
                "--adapter-sample-reads has no effect with --split-by alone (presence detection \
                 applies to --adapter-preset, and sampling to --adapter-discover)",
            ));
        }
        return Ok(None);
    }
    search_config(c, adapter_seqs).map(Some)
}

/// Builds the adapter search settings for `adapters` from
/// `--adapter-error-rate`, `--adapter-end-search`, `--adapter-ends-only` and
/// `--min-length`, validating the first two.
fn search_config(
    c: &Cli,
    adapters: Vec<crate::adapter::Adapter>,
) -> anyhow::Result<crate::adapter::AdapterConfig> {
    let error_rate = c.adapter_error_rate.unwrap_or(DEFAULT_ADAPTER_ERROR_RATE);
    if !(0.0..=1.0).contains(&error_rate) {
        anyhow::bail!("--adapter-error-rate ({error_rate}) must be between 0 and 1");
    }
    let end_size = c.adapter_end_size.unwrap_or(DEFAULT_ADAPTER_END_SIZE);
    if end_size == 0 {
        anyhow::bail!("--adapter-end-search must be >= 1");
    }
    Ok(crate::adapter::AdapterConfig {
        adapters,
        error_rate,
        end_size,
        split: !c.adapter_ends_only,
        min_piece: c.min_length,
        candidate_index: std::sync::OnceLock::new(),
        amplicon: false,
        split_of: Vec::new(),
        split_opens: Vec::new(),
    })
}

/// Resolves `--split-by` and its options, or `None` without `--split-by`:
/// loads every given value, merges them into one sheet (`Sheet::load_all`),
/// and validates it against `--split-require`. A `mab114` source also
/// selects `--adapter-preset mab114` when neither `-a` nor
/// `--adapter-preset` is given; a sheet file or inline target enables no
/// adapters. A `-o` path holding placeholders is the output template, whose
/// key placeholder sets the key level (`{target}` over `{group}`) and whose
/// placeholders require every sheet name they take to be one path
/// component; without one the run is tag-only with one key per target and
/// names are unrestricted. Output placeholders without `--split-by` are an
/// error.
fn resolve_split(c: &mut Cli) -> anyhow::Result<Option<SplitOptions>> {
    if c.split_by.is_empty() {
        if let Some(out) = c.output.as_deref()
            && Template::has_placeholder(out)
        {
            anyhow::bail!(
                "-o {}: the placeholders {{target}}, {{group}} and {{barcode}} name split \
                 outputs and require --split-by",
                out.display()
            );
        }
        return Ok(None);
    }
    let specs = c.split_by.clone();
    let template = match c.output.as_deref() {
        Some(out) => Template::parse(out)?,
        None => None,
    };
    let sources = specs
        .iter()
        .map(std::path::PathBuf::from)
        .filter(|path| path.is_file())
        .collect();
    let sheet = Sheet::load_all(&specs)?;
    sheet.validate(c.split_require)?;
    if let Some(template) = &template {
        template.check_names(&sheet)?;
    }
    let preset = specs.iter().find(|spec| {
        !std::path::Path::new(spec.as_str()).is_file() && Sheet::preset(spec).is_some()
    });
    if let Some(preset) = preset
        && c.adapter_fasta.is_none()
        && c.adapter_preset.is_none()
    {
        c.adapter_preset = Some(preset.to_ascii_lowercase());
    }
    Ok(Some(SplitOptions {
        sheet,
        level: template.as_ref().map_or(KeyLevel::Target, Template::key),
        rules: Rules {
            require: c.split_require,
            lead: c.split_lead,
        },
        action: c.split_action,
        discard_unassigned: c.split_discard.contains(&SplitBin::Unassigned),
        discard_ambiguous: c.split_discard.contains(&SplitBin::Ambiguous),
        spec: specs,
        sources,
        template,
    }))
}

/// Parses `--adapter-preset` into the kits it names, or `None` when it is
/// absent or names nothing.
fn preset_kits(c: &Cli) -> anyhow::Result<Option<Vec<crate::adapter::preset::Kit>>> {
    let Some(spec) = &c.adapter_preset else {
        return Ok(None);
    };
    let kits = crate::adapter::preset::parse_presets(spec)
        .map_err(|e| anyhow::anyhow!("--adapter-preset: {e}"))?;
    Ok((!kits.is_empty()).then_some(kits))
}

/// Rejects an explicit adapter tuning flag given without an adapter source.
fn require_adapter_source(c: &Cli) -> anyhow::Result<()> {
    let explicit = [
        ("--adapter-error-rate", c.adapter_error_rate.is_some()),
        ("--adapter-end-search", c.adapter_end_size.is_some()),
        ("--adapter-sample-reads", c.adapter_sample.is_some()),
    ];
    if let Some((flag, _)) = explicit.iter().find(|(_, given)| *given) {
        anyhow::bail!(
            "{flag} requires an adapter source (--adapter-fasta, --adapter-preset, or \
             --adapter-discover)"
        );
    }
    Ok(())
}

/// Resolves the sample size for presence detection or inference.
///
/// An omitted value means the mode default: 2000 with inference off, 40000
/// with inference on. An explicit value must be 0 or at least
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
        None => DEFAULT_PRESET_SAMPLE,
        Some(n) => {
            if n != 0 && n < min {
                anyhow::bail!(
                    "--adapter-sample-reads ({n}) must be 0 (disable detection) or at least {min} \
                     (smaller samples are too few for reliable detection)"
                );
            }
            if n == 0 && adapter_infer != AdapterInfer::Off {
                anyhow::bail!(
                    "--adapter-sample-reads 0 disables sampling, which adapter discovery requires; \
                     omit it or pass >= {min}"
                );
            }
            n
        },
    };
    if adapter_infer != AdapterInfer::Off || c.adapter_fasta.is_none() {
        return Ok(requested);
    }
    if c.adapter_sample.is_some_and(|n| n > 0) {
        advisories.push(Advisory::warn(
            "--adapter-sample-reads is ignored with --adapter-fasta (presence detection is \
             preset-only)",
        ));
    }
    Ok(0)
}

/// Returns the role a FASTA header assigns: `primer` or `barcode` as a
/// whole word in the description after the name selects that role, and any
/// other header is an adapter.
fn fasta_role(head: &str) -> crate::adapter::Role {
    let description = head
        .split_once(char::is_whitespace)
        .map_or("", |(_, rest)| rest);
    let mut words = description
        .split(|c: char| c.is_whitespace() || matches!(c, '=' | ':' | ',' | ';'))
        .filter(|w| !w.is_empty());
    match words.find(|w| w.eq_ignore_ascii_case("primer") || w.eq_ignore_ascii_case("barcode")) {
        Some(w) if w.eq_ignore_ascii_case("primer") => crate::adapter::Role::Primer,
        Some(_) => crate::adapter::Role::Barcode,
        None => crate::adapter::Role::Adapter,
    }
}

/// Reads adapter sequences from a FASTA. Whitespace is removed, lowercase is
/// uppercased, and `U` is folded to `T`. IUPAC ambiguity codes are kept and
/// searched as the bases they stand for; an entry containing any other byte is
/// skipped with a warning advisory, as is an entry shorter than
/// `adapter::MIN_PATTERN_LEN`, the matcher's minimum pattern length. An entry
/// averaging two or more bases per position is kept with a warning advisory.
/// The header description selects the role; see `fasta_role`.
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
        let head = String::from_utf8_lossy(rec.head()).into_owned();
        let role = fasta_role(&head);
        let name = head
            .split_once(char::is_whitespace)
            .map_or(head.as_str(), |(name, _)| name)
            .to_string();

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

        out.push(crate::adapter::Adapter { name, seq, role });
    }
    Ok(out)
}

/// Builds a `Config` for integration tests: BAM in and out, the given crops,
/// one thread, quiet.
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
        trim: TrimPlan {
            head: head_crop,
            tail: tail_crop,
            quality: None,
        },
        threads: threads.max(1),
        quiet: true,
        ..Config::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn unsupported_parameter_names_are_rejected() {
        for flag in [
            "--in-format",
            "--out-format",
            "--ordered",
            "--min-qual",
            "--max-qual",
            "--qual-mode",
            "--qual-trim",
            "--qual-best-segment",
            "--qual-split",
            "--qual-split-window",
            "--trim-quality",
            "--best-quality-segment",
            "--split-quality",
            "--split-min-low-quality-bases",
            "--split-quality-segments",
            "--update-signal-tags",
            "--remove-kinetics",
            "--discover-adapters",
            "--strip-kinetics",
            "--adapter-end-size",
            "--adapter-sample",
            "--adapter-infer",
            "--adapter-infer-policy",
        ] {
            let error = Cli::try_parse_from(["whittle", flag]).unwrap_err();
            assert_eq!(
                error.kind(),
                clap::error::ErrorKind::UnknownArgument,
                "{flag}"
            );
        }
    }

    /// The two discovery modes are exclusive, so asking for both is a parse
    /// error rather than a silent choice.
    #[test]
    fn adapter_discover_and_report_conflict() {
        let error =
            Cli::try_parse_from(["whittle", "--adapter-discover", "--adapter-report"]).unwrap_err();
        assert_eq!(error.kind(), clap::error::ErrorKind::ArgumentConflict);
    }

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
            "--adapter-sample-reads help quotes MIN_SAMPLE_FOR_DETECTION"
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
