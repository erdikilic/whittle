//! Live, thread-shared `--split-by` counters: one `KeyCounters` per sheet key,
//! the by-reason unassigned tally, and the ambiguous and discarded totals.
//! Built once the run's `Splitter` exists, sized by `Keys::names`, and
//! snapshotted into a plain `SplitStats` for `Stats` and the summary.

use std::sync::atomic::{AtomicU64, Ordering};

use super::classify::{Call, Ends, Strand, Unassigned};

/// Index into `SplitCounters::unassigned` and `SplitStats::unassigned`,
/// matching the summary's field order (`no_primer`, `require`,
/// `orientation`, `length`).
fn unassigned_index(reason: Unassigned) -> usize {
    match reason {
        Unassigned::NoPrimer => 0,
        Unassigned::Require => 1,
        Unassigned::Orientation => 2,
        Unassigned::Length => 3,
    }
}

/// Live per-key `--split-by` counters.
#[derive(Debug, Default)]
pub struct KeyCounters {
    pub reads: AtomicU64,
    pub bases: AtomicU64,
    pub both_ends: AtomicU64,
    pub five_only: AtomicU64,
    pub three_only: AtomicU64,
    pub plus: AtomicU64,
    pub minus: AtomicU64,
}

impl KeyCounters {
    fn snapshot(&self) -> KeyStats {
        KeyStats {
            reads: self.reads.load(Ordering::Relaxed),
            bases: self.bases.load(Ordering::Relaxed),
            both_ends: self.both_ends.load(Ordering::Relaxed),
            five_only: self.five_only.load(Ordering::Relaxed),
            three_only: self.three_only.load(Ordering::Relaxed),
            plus: self.plus.load(Ordering::Relaxed),
            minus: self.minus.load(Ordering::Relaxed),
        }
    }
}

/// Live, thread-shared `--split-by` counters for one run: one `KeyCounters` per
/// sheet key, in `Keys::names` order, the by-reason unassigned tally, and the
/// ambiguous and discarded totals.
#[derive(Debug)]
pub struct SplitCounters {
    keys: Vec<KeyCounters>,
    unassigned: [AtomicU64; 4],
    ambiguous: AtomicU64,
    discarded: AtomicU64,
}

impl SplitCounters {
    /// Builds zeroed counters for `n_keys` sheet keys.
    pub fn new(n_keys: usize) -> Self {
        SplitCounters {
            keys: (0..n_keys).map(|_| KeyCounters::default()).collect(),
            unassigned: [
                AtomicU64::new(0),
                AtomicU64::new(0),
                AtomicU64::new(0),
                AtomicU64::new(0),
            ],
            ambiguous: AtomicU64::new(0),
            discarded: AtomicU64::new(0),
        }
    }

    /// Records one classified segment that passed the post-trim filters:
    /// `call`'s key, unassigned reason, or ambiguous total is bumped, and an
    /// assigned call also adds `bases` and its end and strand breakdown.
    /// `discarded` additionally bumps the discarded total; the segment still
    /// counts toward its classification either way, so `assigned +
    /// unassigned + ambiguous` always equals `output segments + discarded`.
    pub fn record(&self, call: &Call, bases: u64, discarded: bool) {
        match call {
            Call::Assigned {
                key, strand, ends, ..
            } => {
                let k = &self.keys[*key];
                k.reads.fetch_add(1, Ordering::Relaxed);
                k.bases.fetch_add(bases, Ordering::Relaxed);
                match ends {
                    Ends::Both => &k.both_ends,
                    Ends::Five => &k.five_only,
                    Ends::Three => &k.three_only,
                }
                .fetch_add(1, Ordering::Relaxed);
                match strand {
                    Strand::Plus => &k.plus,
                    Strand::Minus => &k.minus,
                }
                .fetch_add(1, Ordering::Relaxed);
            },
            Call::Unassigned(reason) => {
                self.unassigned[unassigned_index(*reason)].fetch_add(1, Ordering::Relaxed);
            },
            Call::Ambiguous => {
                self.ambiguous.fetch_add(1, Ordering::Relaxed);
            },
        }
        if discarded {
            self.discarded.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Snapshots every counter into a plain `SplitStats`.
    pub fn snapshot(&self) -> SplitStats {
        SplitStats {
            keys: self.keys.iter().map(KeyCounters::snapshot).collect(),
            unassigned: [
                self.unassigned[0].load(Ordering::Relaxed),
                self.unassigned[1].load(Ordering::Relaxed),
                self.unassigned[2].load(Ordering::Relaxed),
                self.unassigned[3].load(Ordering::Relaxed),
            ],
            ambiguous: self.ambiguous.load(Ordering::Relaxed),
            discarded: self.discarded.load(Ordering::Relaxed),
        }
    }
}

/// A snapshot of one key's `KeyCounters`.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct KeyStats {
    pub reads: u64,
    pub bases: u64,
    pub both_ends: u64,
    pub five_only: u64,
    pub three_only: u64,
    pub plus: u64,
    pub minus: u64,
}

/// A snapshot of `SplitCounters`, keys in `Keys::names` order. `unassigned`
/// is indexed as `[no_primer, require, orientation, length]`.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SplitStats {
    pub keys: Vec<KeyStats>,
    pub unassigned: [u64; 4],
    pub ambiguous: u64,
    pub discarded: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assigned(key: usize, strand: Strand, ends: Ends) -> Call {
        Call::Assigned {
            key,
            target: 0,
            strand,
            ends,
        }
    }

    #[test]
    fn assigned_call_bumps_its_key_and_breakdown() {
        let counters = SplitCounters::new(2);
        counters.record(&assigned(0, Strand::Plus, Ends::Both), 100, false);
        counters.record(&assigned(0, Strand::Minus, Ends::Five), 40, false);
        counters.record(&assigned(1, Strand::Plus, Ends::Three), 30, false);

        let snap = counters.snapshot();
        assert_eq!(
            snap.keys[0],
            KeyStats {
                reads: 2,
                bases: 140,
                both_ends: 1,
                five_only: 1,
                three_only: 0,
                plus: 1,
                minus: 1,
            }
        );
        assert_eq!(
            snap.keys[1],
            KeyStats {
                reads: 1,
                bases: 30,
                both_ends: 0,
                five_only: 0,
                three_only: 1,
                plus: 1,
                minus: 0,
            }
        );
        assert_eq!(snap.unassigned, [0, 0, 0, 0]);
        assert_eq!(snap.ambiguous, 0);
        assert_eq!(snap.discarded, 0);
    }

    #[test]
    fn unassigned_call_bumps_its_reason_in_summary_field_order() {
        let counters = SplitCounters::new(1);
        counters.record(&Call::Unassigned(Unassigned::NoPrimer), 0, false);
        counters.record(&Call::Unassigned(Unassigned::Require), 0, false);
        counters.record(&Call::Unassigned(Unassigned::Orientation), 0, false);
        counters.record(&Call::Unassigned(Unassigned::Length), 0, false);
        counters.record(&Call::Unassigned(Unassigned::NoPrimer), 0, false);

        let snap = counters.snapshot();
        assert_eq!(snap.unassigned, [2, 1, 1, 1]);
        assert!(snap.keys.iter().all(|k| k.reads == 0));
    }

    #[test]
    fn ambiguous_call_bumps_the_ambiguous_total() {
        let counters = SplitCounters::new(1);
        counters.record(&Call::Ambiguous, 0, false);
        counters.record(&Call::Ambiguous, 0, false);
        assert_eq!(counters.snapshot().ambiguous, 2);
    }

    /// `discarded` is orthogonal to the classification: an assigned call is
    /// never discarded in practice (`Splitter::discards` never drops one),
    /// but the counter itself makes no such assumption.
    #[test]
    fn discarded_adds_to_the_discarded_total_without_changing_the_classification() {
        let counters = SplitCounters::new(1);
        counters.record(&Call::Unassigned(Unassigned::NoPrimer), 0, true);
        counters.record(&Call::Ambiguous, 0, true);
        counters.record(&Call::Ambiguous, 0, false);

        let snap = counters.snapshot();
        assert_eq!(snap.unassigned, [1, 0, 0, 0]);
        assert_eq!(snap.ambiguous, 2);
        assert_eq!(snap.discarded, 2);
    }

    /// A key with no reads still appears in the snapshot at its index, zeroed.
    #[test]
    fn zero_count_key_stays_in_the_snapshot() {
        let counters = SplitCounters::new(3);
        counters.record(&assigned(1, Strand::Plus, Ends::Both), 10, false);
        let snap = counters.snapshot();
        assert_eq!(snap.keys.len(), 3);
        assert_eq!(snap.keys[0], KeyStats::default());
        assert_eq!(snap.keys[2], KeyStats::default());
    }
}
