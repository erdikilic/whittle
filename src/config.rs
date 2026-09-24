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

/// The per-base arrays the `kinetics` removal group names: the PacBio kinetics
/// and alignment-count arrays a BAM writer slices with the sequence.
pub fn kinetics_tags() -> impl Iterator<Item = [u8; 2]> {
    [
        *b"ip", *b"pw", *b"fi", *b"fp", *b"ri", *b"rp", *b"sa", *b"sm", *b"sx",
    ]
    .into_iter()
}

/// The ONT signal tags the `signal` removal group names and the BAM writer
/// maintains under `--update-moves`.
pub const SIGNAL_TAGS: [[u8; 2]; 5] = [*b"mv", *b"ts", *b"ns", *b"sp", *b"pi"];

/// The modification block the `mods` removal group names.
pub const MOD_TAGS: [[u8; 2]; 3] = [*b"MM", *b"ML", *b"MN"];

/// Aux tags removed from every output record, filled by `--remove-tag` from
/// tags and group names. Removal runs after the rewrite of a tag whittle
/// maintains, so a removed `MM` or `mv` leaves the rest of the record intact.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TagRemoval {
    /// The removed tags, sorted.
    tags: BTreeSet<[u8; 2]>,
}

impl TagRemoval {
    /// Parses the `--remove-tag` values. Each value is a comma-separated list
    /// whose items are two ASCII alphanumeric characters, the shape of a SAM
    /// tag, or one of the group names `kinetics`, `mods` and `signal`.
    pub fn parse(values: &[String]) -> anyhow::Result<Self> {
        let mut tags = BTreeSet::new();
        for item in values.iter().flat_map(|v| v.split(',')) {
            let item = item.trim();
            match item {
                "kinetics" => tags.extend(kinetics_tags()),
                "mods" => tags.extend(MOD_TAGS),
                "signal" => tags.extend(SIGNAL_TAGS),
                _ if item.len() == 2 && item.bytes().all(|c| c.is_ascii_alphanumeric()) => {
                    let b = item.as_bytes();
                    tags.insert([b[0], b[1]]);
                },
                _ => anyhow::bail!(
                    "--remove-tag: invalid item {item:?} (a SAM tag is exactly 2 alphanumeric \
                     characters, such as `ML` or `RG`; the groups are kinetics, mods and signal)"
                ),
            }
        }
        Ok(TagRemoval { tags })
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
}

/// Whether de novo adapter discovery runs and what it does with its findings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterInfer {
    /// Discovery does not run.
    Off,
    /// Discovery runs and the run trims with the discovered sequences.
    Discover,
    /// Discovery runs, prints the discovered FASTA and the run ends.
    Report,
}

impl AdapterInfer {
    /// Whether discovery runs in report mode.
    pub fn is_report(self) -> bool {
        matches!(self, Self::Report)
    }

    /// The mode name used by the banner and the summary.
    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Discover => "discover",
            Self::Report => "report",
        }
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
    /// Input format forced by `--input-format`.
    pub in_format: Option<Format>,
    /// Output format forced by `--output-format`.
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
    /// Whether de novo adapter inference runs and whether inferred adapters
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
    /// Aux tags removed from every output record (`--remove-tag`,
    /// its groups). Requires BAM or tagged FASTQ;
    /// `guards::guard_tag_flags` rejects input without auxiliary tags.
    pub remove_tags: TagRemoval,
    /// The `--tag-filter` expressions; a read must satisfy every one to be
    /// processed. Requires BAM or tagged FASTQ input.
    pub tag_filters: crate::tagfilter::TagFilters,
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
    /// Destination for rejected reads and segments (`--rejected-output`), in the
    /// output's format family with a `wr:Z` reason tag, or `None`.
    pub rejected_output: Option<PathBuf>,
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
    /// Whether the resolved adapter set trims adapters, primers and barcodes,
    /// in dorado's `tm` token order. Recorded by `settle` and merged into the
    /// `@RG` and per-read `tm` fields.
    pub trim_classes: [bool; 3],
}

impl Default for Config {
    /// The library defaults: stdin to stdout, no trimming or filtering, one
    /// thread, every tag, compression level 6, progress `auto`. `cli::parse`
    /// overrides the thread count with the CPU count and the level with 4 for
    /// `.gz` output.
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
            remove_tags: TagRemoval::default(),
            tag_filters: crate::tagfilter::TagFilters::default(),
            ordered: false,
            verbosity: 0,
            quiet: false,
            threads_clamped: None,
            summary_json: None,
            rejected_output: None,
            advisories: Vec::new(),
            progress: ProgressMode::Auto,
            adapter_fasta: None,
            adapters_configured: None,
            trim_classes: [false; 3],
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
            ("--rejected-output", self.rejected_output.as_deref()),
            ("--summary-json", self.summary_json.as_deref()),
        ]
        .into_iter()
        .filter_map(|(flag, path)| path.map(|p| (flag, p)))
    }
}

/// How a `-t` total worker budget is spent. The render pool trims, rebuilds
/// tags and compresses the output blocks; BGZF input takes decode workers out
/// of the budget, which block whenever the pool is behind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThreadBudget {
    /// Workers for input decoding.
    pub decode: usize,
    /// Workers for the render pool.
    pub render: usize,
}

/// Returns the budget for `total` workers, which the two stages share.
/// `parallel_decode` names BGZF input, whose blocks inflate in parallel: a
/// quarter of the budget, at least one, decodes ahead of the pool and the rest
/// renders. Other input decodes on the pool, which then holds the whole budget.
pub fn thread_budget(total: usize, parallel_decode: bool) -> ThreadBudget {
    let total = total.max(1);
    if parallel_decode && total > 1 {
        let decode = (total / 4).max(1);
        ThreadBudget {
            decode,
            render: total - decode,
        }
    } else {
        ThreadBudget {
            decode: 1,
            render: total,
        }
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

    /// The stages share the budget: BGZF input takes a quarter of it, at least
    /// one, for decoding and the rest renders; other input renders on every
    /// worker, and a single worker stays sequential.
    #[test]
    fn thread_budget_shares_the_workers_between_stages() {
        assert_eq!(
            thread_budget(8, true),
            ThreadBudget {
                decode: 2,
                render: 6
            }
        );
        assert_eq!(
            thread_budget(32, true),
            ThreadBudget {
                decode: 8,
                render: 24
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
                render: 1
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

    /// The `kinetics` group folds in exactly the nine per-base arrays the BAM
    /// writer slices, so the group and the writer cannot drift apart.
    #[test]
    fn kinetics_group_folds_in_the_nine_per_base_arrays() {
        let r = TagRemoval::parse(&["kinetics".to_string()]).unwrap();
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

    /// Comma lists, repeated values and groups fill one set, so the writers
    /// have a single removal path.
    #[test]
    fn comma_lists_and_groups_share_one_set() {
        let r = TagRemoval::parse(&["MM,RG".to_string(), "kinetics".to_string()]).unwrap();
        assert!(r.contains(b"MM") && r.contains(b"RG") && r.contains(b"ip"));
        assert_eq!(r.tags().count(), 11);
        let r = TagRemoval::parse(&["mods,signal".to_string()]).unwrap();
        assert_eq!(r.tags().count(), 8);
        assert!(r.contains(b"ML") && r.contains(b"mv") && r.contains(b"pi"));
        assert!(TagRemoval::parse(&[]).unwrap().is_empty());
    }

    /// Group names are exact lowercase words; any other spelling is rejected
    /// rather than read as a tag or a group.
    #[test]
    fn group_names_are_lowercase_only() {
        for item in ["Kinetics", "KINETICS", "kinetic", "MODS", "xyz"] {
            let err = TagRemoval::parse(&[item.to_string()]).unwrap_err();
            assert!(err.to_string().contains("invalid item"), "{item}: {err}");
        }
    }

    #[test]
    fn no_flag_removes_nothing() {
        let r = TagRemoval::parse(&[]).unwrap();
        assert!(r.is_empty());
        assert!(!r.contains(b"MM"));
        assert_eq!(r, TagRemoval::default());
    }

    /// A value that is not a two-character SAM tag is rejected, and the message
    /// names the flag.
    #[test]
    fn remove_tag_rejects_a_malformed_value() {
        for bad in ["M", "MMM", "M_", "", "\u{e9}"] {
            let err = TagRemoval::parse(&[bad.to_string()])
                .unwrap_err()
                .to_string();
            assert!(err.starts_with("--remove-tag:"), "{bad:?}: {err}");
        }
        assert!(TagRemoval::parse(&["M1".to_string()]).is_ok());
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
