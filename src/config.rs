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
    /// Carry no tags — emit plain FASTQ.
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
    /// `--adapter-preset ont` was given (adapter trimming off — no per-read cost).
    pub adapters: Option<crate::adapter::AdapterConfig>,
    pub threads: usize,
    pub fastq_tags: FastqTags,
    /// Resolved render-pool size for this dispatch; `0` means "fall back to
    /// `threads`" (used by tests and any caller that hasn't computed a
    /// workload-aware budget). Set by `run`/`run_folder` from
    /// `thread_budget(..).render` before the pipeline runs.
    pub render_workers: usize,
    /// DEFLATE compression level (0-9) for compressed output: bgzf for BAM,
    /// gzip for FASTQ.gz. `6` is the bgzf/gzip default; lower it (e.g. 1-3) to
    /// trade ratio for speed on the compression-bound BAM path. Plain FASTQ
    /// output ignores it. Validated to 0..=9 by `cli::parse`.
    pub compression_level: u8,
    /// When true, keep ONT signal tags consistent through trimming instead of
    /// dropping them: slice the `mv` move table and update `ts`/`ns`/`sp`/`pi`
    /// (BAM→BAM only — see `pipeline::bam`). Default false drops `mv`/`ts`/`ns`/
    /// `sp`/`pi` on any trimmed read.
    pub update_moves: bool,
    pub verbosity: u8,
    pub quiet: bool,
    /// `Some((requested, ncpu))` when `-t` was clamped down; drives a warning in `run`.
    pub threads_clamped: Option<(usize, usize)>,
}

/// How a `-t` total worker budget splits across the pipeline stages. The split
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
    /// Sum across all three stages — the resolved total worker count shown in
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
    /// bgzf (BAM out) — libdeflate, medium cost.
    Bgzf,
    /// gzip via gzp (FASTQ.gz out) — heavier.
    Gzip,
}

/// Split a `-t` total worker budget across decode/render/encode, given whether
/// RENDER is heavy (BAM input → MM/ML reconstruction) vs light (FASTQ input →
/// trim only), and the encode stage's kind. Empirically tuned (2026-07-03 sweep,
/// mid_eqbase): decode never benefits (serial inflate keeps up → 1); the rest is
/// split so the heavier stage gets more threads.
pub fn thread_budget(total: usize, render_heavy: bool, encode: EncodeKind) -> ThreadBudget {
    let total = total.max(1);
    let rest = total.saturating_sub(1).max(2); // >= 2 so both stages can get >= 1
    let (render, encode_n) = match (render_heavy, encode) {
        // No compression pool → render gets everything (encode field unused).
        (_, EncodeKind::None) => (rest, 1),
        // BAM in + bgzf out: render slightly favored (bgzf/libdeflate is fast). R4E3.
        (true, EncodeKind::Bgzf) => (rest.div_ceil(2), rest / 2),
        // BAM in + gzip out: both heavy, encode slightly favored. R3E4.
        (true, EncodeKind::Gzip) => (rest / 2, rest.div_ceil(2)),
        // FASTQ in (light render) + any compression: encode dominates. R1E6.
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

#[cfg(test)]
mod tests {
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
            thread_budget(8, true, Bgzf),
            ThreadBudget {
                decode: 1,
                render: 4,
                encode: 3
            }
        );
        assert_eq!(
            thread_budget(8, true, Gzip),
            ThreadBudget {
                decode: 1,
                render: 3,
                encode: 4
            }
        );
        assert_eq!(
            thread_budget(8, false, Gzip),
            ThreadBudget {
                decode: 1,
                render: 1,
                encode: 6
            }
        );
        assert_eq!(
            thread_budget(8, true, None),
            ThreadBudget {
                decode: 1,
                render: 7,
                encode: 1
            }
        );
        for t in [1usize, 2, 3, 4, 16] {
            for rh in [true, false] {
                for e in [None, Bgzf, Gzip] {
                    let b = thread_budget(t, rh, e);
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
                    let b = thread_budget(t, rh, e);
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
}
