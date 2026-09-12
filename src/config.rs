//! Resolved run configuration: `Config`, tag carry-through policy,
//! adapter-inference settings, progress mode, and the thread budget.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use crate::filter::FilterConfig;
use crate::io::Format;
use crate::trim::TrimPlan;

/// Which aux tags to carry into FASTQ headers on BAM-to-FASTQ conversion.
/// `MM`/`ML`/`MN` are reconstructed (trim-aware); every other carried tag is
/// copied verbatim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FastqTags {
    /// Carry every aux tag from the source record.
    All,
    /// Carry no tags, emitting plain FASTQ.
    None,
    /// Carry only the listed 2-character SAM tags.
    Only(BTreeSet<[u8; 2]>),
}

impl FastqTags {
    /// Parses a `--fastq-tags` spec: `all`, `none`, or a comma list of
    /// two-character tags.
    pub fn parse(s: &str) -> anyhow::Result<Self> {
        match s {
            "all" => Ok(FastqTags::All),
            "none" => Ok(FastqTags::None),
            _ => {
                let mut set = BTreeSet::new();
                for tok in s.split(',') {
                    if tok.len() != 2 || !tok.bytes().all(|c| c.is_ascii_alphanumeric()) {
                        anyhow::bail!(
                            "--fastq-tags: invalid tag {tok:?} (SAM tags are exactly 2 \
                             characters); use `all`, `none`, or a comma list like `MM,ML,RG`"
                        );
                    }
                    let b = tok.as_bytes();
                    set.insert([b[0], b[1]]);
                }
                Ok(FastqTags::Only(set))
            },
        }
    }

    /// Whether a non-mod tag is carried.
    pub fn carries(&self, tag: &[u8; 2]) -> bool {
        match self {
            FastqTags::All => true,
            FastqTags::None => false,
            FastqTags::Only(s) => s.contains(tag),
        }
    }

    /// Whether the reconstructed `MM`/`ML`/`MN` block is carried. The block is a
    /// unit: on under `All`, or when an explicit list contains `MM` or `ML`.
    pub fn carries_mods(&self) -> bool {
        match self {
            FastqTags::All => true,
            FastqTags::None => false,
            FastqTags::Only(s) => s.contains(b"MM") || s.contains(b"ML"),
        }
    }
}

/// The per-base arrays `--strip-kinetics` removes: the PacBio kinetics
/// `ip`/`pw`/`fi`/`fp`, the reverse-strand `ri`/`rp`, the per-base aligned match
/// and mismatch counts `sm`/`sx`, and the run-length subread coverage `sa`.
/// Derived from the per-base tag constants, so the flag covers exactly the
/// arrays the BAM writer slices.
pub fn kinetics_tags() -> impl Iterator<Item = [u8; 2]> {
    crate::workflow::bam::KNOWN_PERBASE_TAGS
        .into_iter()
        .chain([crate::workflow::bam::RLE_COVERAGE_TAG])
}

/// Aux tags removed from every output record. `--remove-tag` names them one at a
/// time and `--strip-kinetics` folds in `kinetics_tags`, so both flags fill one
/// set and the writers have a single removal path. Removal runs after the
/// rewrite of a tag whittle maintains, so a removed `MM` or `mv` leaves the rest
/// of the record intact.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TagRemoval {
    /// The removed tags, sorted.
    tags: BTreeSet<[u8; 2]>,
    /// Whether `--remove-tag` named at least one tag.
    named: bool,
    /// Whether `--strip-kinetics` was given.
    kinetics: bool,
}

impl TagRemoval {
    /// Parses the `--remove-tag` values and folds in `--strip-kinetics`. Each
    /// value is exactly two ASCII alphanumeric characters, the shape of a SAM
    /// tag.
    pub fn parse(values: &[String], strip_kinetics: bool) -> anyhow::Result<Self> {
        let mut tags = BTreeSet::new();
        for value in values {
            if value.len() != 2 || !value.bytes().all(|c| c.is_ascii_alphanumeric()) {
                anyhow::bail!(
                    "--remove-tag: invalid tag {value:?} (SAM tags are exactly 2 alphanumeric \
                     characters, such as `ML` or `RG`)"
                );
            }
            let b = value.as_bytes();
            tags.insert([b[0], b[1]]);
        }
        let named = !tags.is_empty();
        if strip_kinetics {
            tags.extend(kinetics_tags());
        }
        Ok(TagRemoval {
            tags,
            named,
            kinetics: strip_kinetics,
        })
    }

    /// Whether nothing is removed, so every writer keeps its pass-through path.
    pub fn is_empty(&self) -> bool {
        self.tags.is_empty()
    }

    /// Whether `tag` is removed from output records.
    pub fn contains(&self, tag: &[u8; 2]) -> bool {
        !self.tags.is_empty() && self.tags.contains(tag)
    }

    /// The removed tags, sorted, for the run summary.
    pub fn tags(&self) -> impl Iterator<Item = &[u8; 2]> {
        self.tags.iter()
    }

    /// Whether `--strip-kinetics` was given, which the run summary records
    /// alongside the resolved set the flag expands to.
    pub fn strips_kinetics(&self) -> bool {
        self.kinetics
    }

    /// The flag or flags that configured the removal, for a diagnostic that has
    /// to name what the user wrote.
    pub fn flags(&self) -> &'static str {
        match (self.named, self.kinetics) {
            (true, true) => "--remove-tag and --strip-kinetics",
            (false, true) => "--strip-kinetics",
            _ => "--remove-tag",
        }
    }
}

/// What an enabled ab-initio adapter inference run does with its discoveries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum AdapterInferAction {
    /// Trim reads with the inferred sequences.
    Trim,
    /// Print inferred FASTA to stdout and do not write read output.
    Report,
}

impl AdapterInferAction {
    /// Lowercase label, as the banner and the summary print it.
    pub fn label(self) -> &'static str {
        match self {
            AdapterInferAction::Trim => "trim",
            AdapterInferAction::Report => "report",
        }
    }
}

/// How much of an inferred recurrent consensus is trusted as technical.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum AdapterInferPolicy {
    /// Use a short end-facing anchor.
    Conservative,
    /// Use the complete recurrent consensus.
    Aggressive,
}

impl AdapterInferPolicy {
    /// Lowercase label, as the banner and the summary print it.
    pub fn label(self) -> &'static str {
        match self {
            AdapterInferPolicy::Conservative => "conservative",
            AdapterInferPolicy::Aggressive => "aggressive",
        }
    }
}

/// Whether ab-initio adapter inference runs, and its independent action and
/// policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterInfer {
    /// Inference does not run.
    Off,
    /// Inference runs with the given action and policy.
    Enabled {
        /// What is done with the discoveries.
        action: AdapterInferAction,
        /// How much of each consensus is trusted.
        policy: AdapterInferPolicy,
    },
}

impl AdapterInfer {
    /// Whether inference is enabled in report mode.
    pub fn is_report(self) -> bool {
        matches!(
            self,
            Self::Enabled {
                action: AdapterInferAction::Report,
                ..
            }
        )
    }

    /// Whether inference is enabled with the aggressive policy.
    pub fn is_aggressive(self) -> bool {
        matches!(
            self,
            Self::Enabled {
                policy: AdapterInferPolicy::Aggressive,
                ..
            }
        )
    }
}

/// A parse-time diagnostic, held until the log subscriber exists.
///
/// `cli::parse` runs before `obs::init`, so a message printed there directly
/// would bypass the level filter (surviving `--quiet`), carry no
/// `[timestamp] [LEVEL]` prefix, and land ahead of the version and command lines
/// that open every run. The messages are collected instead, and `run` emits them
/// through tracing with the other deferred advisories.
#[derive(Debug, Clone)]
pub struct Advisory {
    /// True for a warning, false for informational.
    pub warn: bool,
    /// The message text, logged verbatim.
    pub message: String,
}

impl Advisory {
    /// Creates a warning-level advisory.
    pub fn warn(message: impl Into<String>) -> Self {
        Advisory {
            warn: true,
            message: message.into(),
        }
    }

    /// Creates an informational advisory.
    pub fn info(message: impl Into<String>) -> Self {
        Advisory {
            warn: false,
            message: message.into(),
        }
    }
}

/// How progress is reported, chosen independently of the log level.
///
/// `--quiet` conflicts with `--progress` at parse time and outranks this
/// setting. Progress and log level are separate so the summary can be kept
/// while in-flight progress lines or the bar are suppressed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum ProgressMode {
    /// A bar on a terminal, periodic lines otherwise.
    Auto,
    /// The animated bar, unless -v or WHITTLE_LOG asks for log lines.
    Bar,
    /// Always periodic lines, never a bar.
    Plain,
    /// No progress reporting; the banner and summary still print.
    None,
}

/// Input and output endpoints and any forced formats.
#[derive(Debug, Clone, Default)]
pub struct IoConfig {
    /// Input path; `None` reads stdin.
    pub input: Option<PathBuf>,
    /// Output path; `None` writes stdout.
    pub output: Option<PathBuf>,
    /// Input format forced by `--in-format`.
    pub in_format: Option<Format>,
    /// Output format forced by `--out-format`.
    pub out_format: Option<Format>,
}

/// The fully resolved settings for one run.
#[derive(Debug, Clone)]
pub struct Config {
    /// Input, output, and format selection.
    pub io: IoConfig,
    /// Post-trim segment filters.
    pub filter: FilterConfig,
    /// Fixed crops and the quality-trim operation.
    pub trim: TrimPlan,
    /// Adapter-trimming settings, or `None` when neither `--adapter-fasta` nor
    /// `--adapter-preset` was given (adapter trimming off, no per-read cost).
    pub adapters: Option<crate::adapter::AdapterConfig>,
    /// Whether ab-initio adapter inference runs and whether inferred adapters
    /// are also used for trimming.
    pub adapter_infer: AdapterInfer,
    /// Resolved worker-thread count.
    pub threads: usize,
    /// Aux tags carried into FASTQ headers on BAM-to-FASTQ output.
    pub fastq_tags: FastqTags,
    /// Resolved render-pool size for this dispatch; `0` means the workflow falls
    /// back to `threads`. Set by `settle` from `thread_budget(..).render`
    /// before the workflow runs.
    pub render_workers: usize,
    /// Reads to sample for adapter-presence detection before trimming the full
    /// dataset. `0` disables detection (trim against the full active set).
    /// Only meaningful when `adapters` is `Some`.
    pub adapter_sample: usize,
    /// DEFLATE compression level (0-9) for compressed output: bgzf for BAM and
    /// `.bgz`, gzip for FASTQ.gz. `cli::parse` defaults it to 4 for gzip FASTQ
    /// and 6 for BGZF and validates an explicit value to 0..=9. Plain FASTQ
    /// output ignores it.
    pub compression_level: u8,
    /// Whether ONT signal tags are kept consistent through trimming: the `mv`
    /// move table is sliced and `ts`/`ns`/`sp`/`pi` are updated (BAM-to-BAM
    /// only, see `workflow::bam`). When false, a trimmed read drops
    /// `mv`/`ts`/`ns`/`sp`/`pi`.
    pub update_moves: bool,
    /// Whether the barcode spans dorado recorded in the `bi` aux tag are
    /// removed, as the outermost trimming stage. BAM input only; `cli::parse`
    /// and `guards::guard_bam_only_flags` reject any other input format.
    pub trim_barcodes: bool,
    /// Aux tags removed from every output record (`--remove-tag`,
    /// `--strip-kinetics`). BAM input only; `cli::parse` and
    /// `guards::guard_bam_only_flags` reject any other input format.
    pub remove_tags: TagRemoval,
    /// Whether multithreaded runs write records in input order. When false,
    /// records are written in completion order.
    pub ordered: bool,
    /// Count of `-v` flags (0 to 2).
    pub verbosity: u8,
    /// Whether `--quiet` was given.
    pub quiet: bool,
    /// `Some((requested, ncpu))` when `-t` was clamped down; drives a warning in
    /// `run`.
    pub threads_clamped: Option<(usize, usize)>,
    /// Destination for the machine-readable run summary (`--summary-json`), or
    /// `None`. Written regardless of `--quiet` and the log level.
    pub summary_json: Option<PathBuf>,
    /// Diagnostics raised while parsing arguments, emitted by `run` once the log
    /// subscriber exists. See `Advisory`.
    pub advisories: Vec<Advisory>,
    /// How progress is reported. See `ProgressMode`.
    pub progress: ProgressMode,
    /// The `--adapter-fasta` path, kept so `run` can refuse to overwrite it.
    /// The sequences themselves are already resolved into `adapters`.
    pub adapter_fasta: Option<PathBuf>,
    /// How many adapters were configured before presence detection narrowed the
    /// set or inference replaced it, so the run summary can report both figures.
    /// Recorded by `settle`, the only thing that changes `adapters`. `None` when
    /// adapter trimming is off, and `0` under inference, where the set is
    /// discovered rather than configured.
    pub adapters_configured: Option<usize>,
}

impl Default for Config {
    /// The command line's defaults: stdin to stdout, no trimming or filtering,
    /// one thread, every tag, compression level 6, progress `auto`.
    fn default() -> Self {
        Config {
            io: IoConfig::default(),
            filter: crate::filter::FilterConfig::default(),
            trim: crate::trim::TrimPlan::default(),
            adapters: None,
            adapter_infer: AdapterInfer::Off,
            threads: 1,
            fastq_tags: FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            trim_barcodes: false,
            remove_tags: TagRemoval::default(),
            ordered: false,
            verbosity: 0,
            quiet: false,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            progress: ProgressMode::Auto,
            adapter_fasta: None,
            adapters_configured: None,
        }
    }
}

impl Config {
    /// Every file the run writes, each paired with the flag that named it. The
    /// overwrite guards and the report-only advisories both derive from this
    /// list, so an artifact flag added here is covered by both.
    pub fn write_targets(&self) -> impl Iterator<Item = (&'static str, &Path)> {
        [
            ("-o/--output", self.io.output.as_deref()),
            ("--summary-json", self.summary_json.as_deref()),
        ]
        .into_iter()
        .filter_map(|(flag, path)| path.map(|p| (flag, p)))
    }
}

/// How a `-t` total worker budget is spent. The render pool trims, rebuilds
/// tags and compresses the output blocks, so it holds the whole budget; BGZF
/// input adds decode workers that block whenever the pool is behind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThreadBudget {
    /// Workers for input decoding.
    pub decode: usize,
    /// Workers for the render pool.
    pub render: usize,
}

/// Returns the budget for `total` workers. `parallel_decode` names BGZF input,
/// whose blocks inflate in parallel: a quarter of the budget, at least one,
/// decodes ahead of the pool.
pub fn thread_budget(total: usize, parallel_decode: bool) -> ThreadBudget {
    let total = total.max(1);
    let decode = if parallel_decode && total > 1 {
        (total / 4).max(1)
    } else {
        1
    };
    ThreadBudget {
        decode,
        render: total,
    }
}

/// Resolves the worker-thread count. `None` (flag omitted) means all available
/// CPUs; `Some(n)` is clamped into `[1, ncpu]`. The caller warns when it clamped
/// down; `cli::parse` rejects 0 before this runs, so the floor covers library
/// callers only.
pub fn resolve_threads(requested: Option<usize>, ncpu: usize) -> usize {
    let ncpu = ncpu.max(1);
    match requested {
        None => ncpu,
        Some(n) => n.clamp(1, ncpu),
    }
}

#[cfg(test)]
mod resolve_threads_tests {
    use super::resolve_threads;

    #[test]
    fn auto_uses_all_cpus() {
        assert_eq!(resolve_threads(None, 8), 8);
    }
    #[test]
    fn in_range_is_unchanged() {
        assert_eq!(resolve_threads(Some(4), 8), 4);
    }
    #[test]
    fn over_spec_clamps_to_ncpu() {
        assert_eq!(resolve_threads(Some(32), 8), 8);
    }
    #[test]
    fn zero_floors_to_one() {
        assert_eq!(resolve_threads(Some(0), 8), 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The render pool holds the whole budget; BGZF input adds a quarter of it
    /// as decode workers, at least one, and a single worker stays sequential.
    #[test]
    fn thread_budget_gives_the_pool_every_worker() {
        assert_eq!(
            thread_budget(8, true),
            ThreadBudget {
                decode: 2,
                render: 8
            }
        );
        assert_eq!(
            thread_budget(8, false),
            ThreadBudget {
                decode: 1,
                render: 8
            }
        );
        assert_eq!(
            thread_budget(2, true),
            ThreadBudget {
                decode: 1,
                render: 2
            }
        );
        assert_eq!(
            thread_budget(1, true),
            ThreadBudget {
                decode: 1,
                render: 1
            }
        );
    }

    /// `--strip-kinetics` folds in exactly the nine per-base arrays the BAM
    /// writer slices, so the flag and the writer cannot drift apart.
    #[test]
    fn strip_kinetics_folds_in_the_nine_per_base_arrays() {
        let r = TagRemoval::parse(&[], true).unwrap();
        let names: Vec<String> = r
            .tags()
            .map(|t| String::from_utf8_lossy(t).into_owned())
            .collect();
        assert_eq!(
            names,
            ["fi", "fp", "ip", "pw", "ri", "rp", "sa", "sm", "sx"]
        );
        for tag in [
            b"ip", b"pw", b"fi", b"fp", b"ri", b"rp", b"sa", b"sm", b"sx",
        ] {
            assert!(r.contains(tag), "{}", String::from_utf8_lossy(tag));
        }
        assert!(!r.contains(b"MM"));
    }

    /// Both flags fill one set, so the writers have a single removal path.
    #[test]
    fn remove_tag_and_strip_kinetics_share_one_set() {
        let r = TagRemoval::parse(&["MM".to_string(), "RG".to_string()], true).unwrap();
        assert!(r.contains(b"MM") && r.contains(b"RG") && r.contains(b"ip"));
        assert_eq!(r.tags().count(), 11);
        assert_eq!(r.flags(), "--remove-tag and --strip-kinetics");
        assert_eq!(
            TagRemoval::parse(&[], true).unwrap().flags(),
            "--strip-kinetics"
        );
        assert_eq!(
            TagRemoval::parse(&["MM".to_string()], false)
                .unwrap()
                .flags(),
            "--remove-tag"
        );
    }

    #[test]
    fn no_flag_removes_nothing() {
        let r = TagRemoval::parse(&[], false).unwrap();
        assert!(r.is_empty());
        assert!(!r.contains(b"MM"));
        assert_eq!(r, TagRemoval::default());
    }

    /// A value that is not a two-character SAM tag is rejected, and the message
    /// names the flag.
    #[test]
    fn remove_tag_rejects_a_malformed_value() {
        for bad in ["M", "MMM", "M_", "", "\u{e9}"] {
            let err = TagRemoval::parse(&[bad.to_string()], false)
                .unwrap_err()
                .to_string();
            assert!(err.starts_with("--remove-tag:"), "{bad:?}: {err}");
        }
        assert!(TagRemoval::parse(&["M1".to_string()], false).is_ok());
    }

    #[test]
    fn parse_all_none() {
        assert_eq!(FastqTags::parse("all").unwrap(), FastqTags::All);
        assert_eq!(FastqTags::parse("none").unwrap(), FastqTags::None);
    }

    #[test]
    fn parse_list_collects_tags() {
        let t = FastqTags::parse("MM,ML,RG").unwrap();
        match t {
            FastqTags::Only(ref s) => {
                assert!(s.contains(b"MM") && s.contains(b"ML") && s.contains(b"RG"));
                assert_eq!(s.len(), 3);
            },
            other => panic!("Expected Only, got {other:?}"),
        }
    }

    #[test]
    fn parse_rejects_bad_token() {
        assert!(FastqTags::parse("MM,ABC").is_err()); // 3-char token
        assert!(FastqTags::parse("").is_err()); // empty -> one empty token
        assert!(FastqTags::parse("MM,").is_err()); // trailing empty token
    }

    #[test]
    fn parse_rejects_non_ascii_two_byte_token() {
        // A single non-ASCII code point encoded as two UTF-8 bytes would pass a
        // length-only check (`b.len() != 2`); it is rejected, while a two-byte
        // ASCII tag parses.
        assert!(FastqTags::parse("é").is_err());
        assert!(FastqTags::parse("RG").is_ok());
    }

    #[test]
    fn carries_rules() {
        assert!(FastqTags::All.carries(b"RG"));
        assert!(FastqTags::All.carries_mods());
        assert!(!FastqTags::None.carries(b"RG"));
        assert!(!FastqTags::None.carries_mods());

        let only = FastqTags::parse("ML,RG").unwrap();
        assert!(only.carries(b"RG"));
        assert!(!only.carries(b"XY"));
        // The mod block is carried when the list contains MM or ML.
        assert!(only.carries_mods());
        // MN alone does not enable the mod block.
        let mn_only = FastqTags::parse("MN").unwrap();
        assert!(!mn_only.carries_mods());
    }
}
