use std::collections::BTreeSet;
use std::path::PathBuf;

use crate::filter::FilterConfig;
use crate::io::Format;
use crate::trim::TrimPlan;

/// Which aux tags to carry into FASTQ headers on BAM→FASTQ conversion.
/// MM/ML/MN are reconstructed (trim-aware); every other carried tag is verbatim.
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
    /// Parse a `--fastq-tags` spec: `all`, `none`, or a comma list of 2-char tags.
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

    /// Whether the reconstructed MM/ML/MN block is carried. The block is a unit:
    /// on under `All`, or when an explicit list contains `MM` or `ML`.
    pub fn carries_mods(&self) -> bool {
        match self {
            FastqTags::All => true,
            FastqTags::None => false,
            FastqTags::Only(s) => s.contains(b"MM") || s.contains(b"ML"),
        }
    }
}

/// What an enabled ab-initio adapter inference run does with its discoveries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterInferAction {
    Trim,
    Report,
}

/// How much of an inferred recurrent consensus is trusted as technical.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterInferPolicy {
    Conservative,
    Aggressive,
}

/// Whether ab-initio adapter inference runs, and its independent action/policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AdapterInfer {
    #[default]
    Off,
    Enabled {
        action: AdapterInferAction,
        policy: AdapterInferPolicy,
    },
}

impl AdapterInfer {
    /// Whether ab-initio inference ran, so the adapter set was discovered
    /// rather than configured.
    pub fn is_enabled(self) -> bool {
        matches!(self, Self::Enabled { .. })
    }

    pub fn is_report(self) -> bool {
        matches!(
            self,
            Self::Enabled {
                action: AdapterInferAction::Report,
                ..
            }
        )
    }

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
/// `cli::parse` runs before `obs::init`, so anything it prints directly bypasses
/// the level filter (surviving `--quiet`), carries no `[timestamp] [LEVEL]`
/// prefix, and lands ahead of the version and command lines that are supposed to
/// open every run. Collecting the messages instead lets `run` emit them through
/// tracing alongside the other deferred advisories.
#[derive(Debug, Clone)]
pub struct Advisory {
    /// True for a warning, false for informational.
    pub warn: bool,
    pub message: String,
}

impl Advisory {
    pub fn warn(message: impl Into<String>) -> Self {
        Advisory {
            warn: true,
            message: message.into(),
        }
    }

    pub fn info(message: impl Into<String>) -> Self {
        Advisory {
            warn: false,
            message: message.into(),
        }
    }
}

/// How progress is reported, chosen independently of the log level.
///
/// `--quiet` silences everything and still wins over this. The point of the
/// separation is that a pipeline may want the run summary without a progress line
/// every thirty seconds in its log, and a terminal user may want the summary
/// without an animated bar.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ProgressMode {
    /// A bar on a terminal, periodic lines otherwise.
    #[default]
    Auto,
    /// Always the animated bar, even when output is redirected.
    Bar,
    /// Always periodic lines, never a bar.
    Plain,
    /// No progress reporting. The banner, warnings and summary still print.
    None,
}

#[derive(Debug, Clone)]
pub struct IoConfig {
    pub input: Option<PathBuf>,
    pub output: Option<PathBuf>,
    pub in_format: Option<Format>,
    pub out_format: Option<Format>,
}

#[derive(Debug, Clone)]
pub struct Config {
    pub io: IoConfig,
    pub filter: FilterConfig,
    pub trim: TrimPlan,
    /// Adapter-trimming settings, or `None` when neither `--adapter-fasta` nor
    /// `--adapter-preset ont` was given (adapter trimming off, no per-read cost).
    pub adapters: Option<crate::adapter::AdapterConfig>,
    /// Whether ab-initio adapter inference runs and whether inferred adapters
    /// are also used for trimming.
    pub adapter_infer: AdapterInfer,
    pub threads: usize,
    pub fastq_tags: FastqTags,
    /// Resolved render-pool size for this dispatch; `0` means "fall back to
    /// `threads`" (used by tests and any caller that hasn't computed a
    /// workload-aware budget). Set by `run`/`run_folder` from
    /// `thread_budget(..).render` before the workflow runs.
    pub render_workers: usize,
    /// Reads to sample for adapter-presence detection before trimming the full
    /// dataset. `0` disables detection (trim against the full active set).
    /// Only meaningful when `adapters` is `Some`.
    pub adapter_sample: usize,
    /// DEFLATE compression level (0-9) for compressed output: bgzf for BAM,
    /// gzip for FASTQ.gz. `6` is the bgzf/gzip default; lower it (e.g. 1-3) to
    /// trade ratio for speed on the compression-bound BAM path. Plain FASTQ
    /// output ignores it. Validated to 0..=9 by `cli::parse`.
    pub compression_level: u8,
    /// When true, keep ONT signal tags consistent through trimming instead of
    /// dropping them: slice the `mv` move table and update `ts`/`ns`/`sp`/`pi`
    /// (BAM→BAM only, see `workflow::bam`). Default false drops `mv`/`ts`/`ns`/
    /// `sp`/`pi` on any trimmed read.
    pub update_moves: bool,
    /// Whether multithreaded runs write records in input order. When false,
    /// records are written in completion order.
    pub ordered: bool,
    pub verbosity: u8,
    pub quiet: bool,
    /// `Some((requested, ncpu))` when `-t` was clamped down; drives a warning in `run`.
    pub threads_clamped: Option<(usize, usize)>,
    /// Where to write the machine-readable run summary (`--summary-json`), or
    /// `None` to write none. Written regardless of `--quiet` and the log level,
    /// since a caller that asked for it is a pipeline that needs the file.
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

/// How a `-t` total worker budget splits across the workflow stages. The split
/// is workload-aware (see `thread_budget`): decode never benefits from more
/// than 1 thread (serial inflate keeps up), while render (MM/ML
/// reconstruction, or the trim-only pass for FASTQ) and encode (bgzf/gzip
/// compression) are weighted against each other based on how heavy each stage
/// actually is for the dispatched (input, output) pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThreadBudget {
    pub decode: usize,
    pub render: usize,
    pub encode: usize,
}

impl ThreadBudget {
    /// Sum across all three stages: the resolved total worker count shown in
    /// the startup banner's `Threads: {total} total (...)` line. May exceed the
    /// requested `-t` value at very low counts, since `thread_budget` floors
    /// `render`/`encode` at >= 1 each even when the overall total is 1.
    pub fn total(&self) -> usize {
        self.decode + self.render + self.encode
    }
}

/// The output compression stage's weight, for thread budgeting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EncodeKind {
    /// No compression pool (plain FASTQ out).
    None,
    /// bgzf (BAM out): libdeflate, medium cost.
    Bgzf,
    /// gzip via gzp (FASTQ.gz out): heavier.
    Gzip,
}

/// Split a `-t` worker budget across decode, render, and encode. Parallel input
/// receives multiple decode workers; otherwise decoding stays serial and the
/// remaining workers are weighted toward the more expensive downstream stage.
pub fn thread_budget(
    total: usize,
    render_heavy: bool,
    parallel_decode: bool,
    encode: EncodeKind,
) -> ThreadBudget {
    let total = total.max(1);

    if parallel_decode && total == 2 {
        return ThreadBudget {
            decode: 2,
            render: 1,
            encode: 1,
        };
    }

    if parallel_decode && total >= 3 {
        let (decode, render, encode_n) = match encode {
            // The `min` matters at `total == 3`, where the `max(2)` floor would
            // otherwise take both remaining workers and leave render with none.
            // One worker is reserved for the writer thread and one for render, so
            // decode can claim at most `total - 2`. Parallel decode gives way to
            // render here rather than the other way round, because a render pool
            // of zero is not a slower configuration, it is a broken one.
            EncodeKind::None if render_heavy => {
                let decode = (total / 3).max(2).min(total - 2);
                (decode, total - decode - 1, 1)
            },
            EncodeKind::None => {
                let decode = (total * 2 / 3).max(2).min(total - 2);
                (decode, total - decode - 1, 1)
            },
            EncodeKind::Gzip | EncodeKind::Bgzf if render_heavy => {
                let render = (total * 2 / 5).max(1);
                let remaining = total - render;
                let decode = (remaining / 3).max(1);
                (decode, render, remaining - decode)
            },
            EncodeKind::Gzip | EncodeKind::Bgzf => {
                let render = 1;
                let remaining = total - render;
                let decode = (remaining / 3).max(1);
                (decode, render, remaining - decode)
            },
        };
        return ThreadBudget {
            decode: decode.max(1),
            render: render.max(1),
            encode: encode_n.max(1),
        };
    }

    let rest = total.saturating_sub(1).max(2); // >= 2 so both stages can get >= 1
    let (render, encode_n) = match (render_heavy, encode) {
        // No compression pool → render gets everything (encode field unused).
        (_, EncodeKind::None) => (rest, 1),
        // BAM in + bgzf out: split nearly evenly. Raw BAM field/tag conversion
        // now runs in the render pool alongside MM/ML reconstruction, while the
        // Noodles BGZF stage uses the other half for ordered block encoding.
        (true, EncodeKind::Bgzf) if rest <= 4 => (rest / 2, rest.div_ceil(2)),
        (true, EncodeKind::Bgzf) if rest <= 8 => (rest.div_ceil(2), rest / 2),
        (true, EncodeKind::Bgzf) => {
            let render = rest.div_ceil(2).max(1);
            (render, rest - render)
        },
        // BAM input with gzip output slightly favors encoding.
        (true, EncodeKind::Gzip) => (rest / 2, rest.div_ceil(2)),
        // FASTQ rendering is light, so compressed output favors encoding.
        (false, _) => {
            let r = (rest / 6).max(1);
            (r, rest - r)
        },
    };
    ThreadBudget {
        decode: 1,
        render: render.max(1),
        encode: encode_n.max(1),
    }
}

/// Resolve the worker-thread count. `None` (flag omitted) → all available CPUs;
/// `Some(n)` → clamp into `[1, ncpu]`. The caller warns when it clamped down.
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

/// The output compression stage's weight for a given output format: `Bgzf` for
/// BAM (always bgzf-compressed), `Gzip` for `FASTQ.gz`, `None` for plain FASTQ.
/// Paired with `render_heavy` (`in_fmt == Format::Bam`, or the folder-mode
/// equivalent), this is everything `config::thread_budget` needs; both call sites
/// (`run`, `run_folder`) resolve their budget from this exactly once, before the
/// startup banner, and reuse it for the actual workflow dispatch below.
pub(crate) fn encode_kind_for(out_fmt: crate::io::Format) -> EncodeKind {
    match out_fmt {
        crate::io::Format::Bam => EncodeKind::Bgzf,
        crate::io::Format::FastqGz => EncodeKind::Gzip,
        crate::io::Format::FastqBgzf => EncodeKind::Bgzf,
        crate::io::Format::Fastq => EncodeKind::None,
    }
}

/// Whether the render stage has substantial per-record work. BAM input remains
/// render-heavy even for a full-window output because the current parallel path
/// still clones owned `RecordBuf`s before handing them to the writer. FASTQ
/// input is normally trim-only (light), but adapter matching or ab-initio
/// inference runs an approximate search per read, which is heavy too, so it
/// gets a render-pool share rather than being starved as pure compression.
pub(crate) fn render_heavy_for(
    in_fmt: crate::io::Format,
    _out_fmt: crate::io::Format,
    cfg: &Config,
) -> bool {
    matches!(in_fmt, crate::io::Format::Bam)
        || cfg.adapters.is_some()
        || cfg.adapter_infer != AdapterInfer::Off
}

#[cfg(test)]
mod tests {

    /// Every stage must get at least one worker at every thread count, for every
    /// combination of inputs. A stage of zero is not a slower configuration, it is
    /// a broken one: the consumer falls back to its own default, so the banner
    /// reports a split the pipeline is not using.
    #[test]
    fn every_stage_gets_at_least_one_worker() {
        for total in 0..=64usize {
            for render_heavy in [false, true] {
                for parallel_decode in [false, true] {
                    for encode in [EncodeKind::None, EncodeKind::Gzip, EncodeKind::Bgzf] {
                        let b = thread_budget(total, render_heavy, parallel_decode, encode);
                        assert!(
                            b.decode >= 1 && b.render >= 1 && b.encode >= 1,
                            "total={total} render_heavy={render_heavy} \
                             parallel_decode={parallel_decode} encode={encode:?} gave {b:?}"
                        );
                    }
                }
            }
        }
    }

    /// With a real compression pool, the three stages must not claim more workers
    /// than were asked for.
    ///
    /// `EncodeKind::None` is excluded on purpose: plain FASTQ output has no encode
    /// pool, so its encode field counts the single writer thread rather than a
    /// share of the budget, and the sum can legitimately read one high. That is
    /// the case `ThreadBudget::total` documents.
    #[test]
    fn split_stays_within_the_requested_budget() {
        for total in 3..=64usize {
            for render_heavy in [false, true] {
                for parallel_decode in [false, true] {
                    for encode in [EncodeKind::Gzip, EncodeKind::Bgzf] {
                        let b = thread_budget(total, render_heavy, parallel_decode, encode);
                        assert!(
                            b.total() <= total,
                            "total={total} render_heavy={render_heavy} \
                             parallel_decode={parallel_decode} encode={encode:?} gave {b:?} \
                             summing to {}",
                            b.total()
                        );
                    }
                }
            }
        }
    }
    use super::*;

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
            other => panic!("expected Only, got {other:?}"),
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
        // "é" is a single codepoint encoded as 2 UTF-8 bytes, so a length-only
        // check (`b.len() != 2`) would wrongly accept it as a "tag". It must
        // be rejected, while a normal 2-byte ASCII tag still parses.
        assert!(FastqTags::parse("é").is_err());
        assert!(FastqTags::parse("RG").is_ok());
    }

    #[test]
    fn thread_budget_split() {
        use EncodeKind::*;
        assert_eq!(
            thread_budget(8, true, false, Bgzf),
            ThreadBudget {
                decode: 1,
                render: 4,
                encode: 3
            }
        );
        assert_eq!(
            thread_budget(16, true, false, Bgzf),
            ThreadBudget {
                decode: 1,
                render: 8,
                encode: 7
            }
        );
        assert_eq!(
            thread_budget(4, true, false, Bgzf),
            ThreadBudget {
                decode: 1,
                render: 1,
                encode: 2
            }
        );
        assert_eq!(
            thread_budget(8, true, false, Gzip),
            ThreadBudget {
                decode: 1,
                render: 3,
                encode: 4
            }
        );
        assert_eq!(
            thread_budget(8, false, false, Gzip),
            ThreadBudget {
                decode: 1,
                render: 1,
                encode: 6
            }
        );
        assert_eq!(
            thread_budget(8, true, false, None),
            ThreadBudget {
                decode: 1,
                render: 7,
                encode: 1
            }
        );
        for t in [1usize, 2, 3, 4, 16] {
            for rh in [true, false] {
                for e in [None, Bgzf, Gzip] {
                    let b = thread_budget(t, rh, false, e);
                    assert!(b.decode >= 1 && b.render >= 1 && b.encode >= 1);
                }
            }
        }
    }

    #[test]
    fn thread_budget_total_sums_all_three_stages() {
        use EncodeKind::*;
        for t in [1usize, 2, 8, 16] {
            for rh in [true, false] {
                for e in [None, Bgzf, Gzip] {
                    let b = thread_budget(t, rh, false, e);
                    assert_eq!(b.total(), b.decode + b.render + b.encode);
                }
            }
        }
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
        // mod block carried when the list has MM *or* ML:
        assert!(only.carries_mods());
        // MN alone does not turn on the mod block:
        let mn_only = FastqTags::parse("MN").unwrap();
        assert!(!mn_only.carries_mods());
    }

    #[test]
    fn bgzf_fastq_plain_output_favors_parallel_decode() {
        assert_eq!(
            thread_budget(16, false, true, EncodeKind::None),
            ThreadBudget {
                decode: 10,
                render: 5,
                encode: 1,
            }
        );
        assert_eq!(
            thread_budget(2, false, true, EncodeKind::None),
            ThreadBudget {
                decode: 2,
                render: 1,
                encode: 1,
            }
        );
    }
}
