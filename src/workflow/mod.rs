mod bam;
mod fastq;

use std::sync::atomic::{AtomicU64, Ordering};

pub use bam::{reconstruct_mods, reconstruct_record, run_bam, run_bam_to_fastq, run_raw_bam};
pub use fastq::{run_fastq, run_fastq_seq};

use crate::filter::{DropReason, FilterConfig};

/// Keep parallel work units large enough to amortize scheduling/channel costs,
/// but small enough to balance unusually long reads across workers.
const BATCH_BASES: usize = 512 * 1024;
const BATCH_RECORDS: usize = 32;
const BAM_BATCH_BASES: usize = 256 * 1024;
const BAM_BATCH_RECORDS: usize = 4;

pub(crate) struct Batches<I, F> {
    records: I,
    weight: F,
    target_weight: usize,
    max_items: usize,
}

impl<I, F> Batches<I, F> {
    pub(crate) fn new(records: I, weight: F) -> Self {
        Self {
            records,
            weight,
            target_weight: BATCH_BASES,
            max_items: BATCH_RECORDS,
        }
    }

    pub(crate) fn bam(records: I, weight: F) -> Self {
        Self {
            records,
            weight,
            target_weight: BAM_BATCH_BASES,
            max_items: BAM_BATCH_RECORDS,
        }
    }
}

impl<I, F, T> Iterator for Batches<I, F>
where
    I: Iterator<Item = T>,
    F: Fn(&T) -> usize,
{
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut batch = Vec::with_capacity(self.max_items);
        let mut bases = 0usize;
        while batch.len() < self.max_items && bases < self.target_weight {
            let Some(record) = self.records.next() else {
                break;
            };
            bases = bases.saturating_add((self.weight)(&record));
            batch.push(record);
        }
        (!batch.is_empty()).then_some(batch)
    }
}

/// Live, thread-shared counters read by the progress ticker and finalized into `Stats`.
#[derive(Default)]
pub struct Counters {
    pub input_reads: AtomicU64,
    pub output_reads: AtomicU64,
    pub bytes_read: AtomicU64,
    /// Sum of SEQ lengths (bases) across every input read, regardless of
    /// whether it survives filtering/trimming.
    pub input_bases: AtomicU64,
    /// Sum of surviving segment lengths (bases) actually written to output.
    pub output_bases: AtomicU64,
    /// Input reads that produced at least one surviving output segment —
    /// bumped once per input read (not once per segment, unlike
    /// `output_reads`, which a `--qual-split` read can bump several times).
    /// Exists so `snapshot`'s `debug_assert_eq!` can check that every input
    /// read is accounted for by exactly one of the three read-level outcomes
    /// (the read-level third of the two-level counter model).
    pub reads_with_output: AtomicU64,
    /// Input reads that produced **zero** segments at all: `trim::apply`
    /// returned no intervals to even run the per-segment filter over (an
    /// empty read, a read fully consumed by adapter trimming, or an
    /// over-crop). Read-level, paired with `reads_with_output` and
    /// `reads_all_filtered` in the invariant below.
    pub reads_trimmed_to_nothing: AtomicU64,
    /// Input reads that produced **at least one** segment, but every one was
    /// rejected by post-trim `filter::check`. Read-level, paired with
    /// `reads_with_output` and `reads_trimmed_to_nothing` in the invariant
    /// below.
    pub reads_all_filtered: AtomicU64,
    /// Segment-level drop counters: one bump per **segment** (not read) that
    /// `filter::check` rejects, by reason, post-trim. A single input read can
    /// contribute to more than one of these (e.g. a `--qual-split` read whose
    /// several pieces are each judged independently) — these are NOT part of
    /// the read-level invariant.
    pub segments_dropped_short: AtomicU64,
    pub segments_dropped_long: AtomicU64,
    pub segments_dropped_low_qual: AtomicU64,
    pub segments_dropped_high_qual: AtomicU64,
    pub segments_dropped_gc: AtomicU64,
}

impl Counters {
    /// Bump the segment-level counter matching a `filter::check` failure
    /// reason. Called once per rejected **segment** (post-trim), not per read.
    pub fn record_segment_drop(&self, reason: DropReason) {
        let counter = match reason {
            DropReason::TooShort => &self.segments_dropped_short,
            DropReason::TooLong => &self.segments_dropped_long,
            DropReason::LowQuality => &self.segments_dropped_low_qual,
            DropReason::HighQuality => &self.segments_dropped_high_qual,
            DropReason::Gc => &self.segments_dropped_gc,
        };
        counter.fetch_add(1, Ordering::Relaxed);
    }

    /// Snapshot every counter into a `Stats` for end-of-run reporting.
    /// `malformed_tag_reads` is threaded through separately: only the BAM
    /// paths track it, and the parallel BAM path accumulates it in its own
    /// local atomic rather than in `Counters`.
    pub fn snapshot(&self, malformed_tag_reads: u64) -> Stats {
        let input_reads = self.input_reads.load(Ordering::Relaxed);
        let reads_with_output = self.reads_with_output.load(Ordering::Relaxed);
        let reads_trimmed_to_nothing = self.reads_trimmed_to_nothing.load(Ordering::Relaxed);
        let reads_all_filtered = self.reads_all_filtered.load(Ordering::Relaxed);
        let segments_dropped_short = self.segments_dropped_short.load(Ordering::Relaxed);
        let segments_dropped_long = self.segments_dropped_long.load(Ordering::Relaxed);
        let segments_dropped_low_qual = self.segments_dropped_low_qual.load(Ordering::Relaxed);
        let segments_dropped_high_qual = self.segments_dropped_high_qual.load(Ordering::Relaxed);
        let segments_dropped_gc = self.segments_dropped_gc.load(Ordering::Relaxed);

        // Every input read is exactly one of: produced at least one
        // surviving output segment (`reads_with_output`, regardless of how
        // many segments it split into), produced no segments at all
        // (`reads_trimmed_to_nothing` — an empty read, a fully-consumed
        // adapter read, or an over-crop), or produced at least one segment
        // but had every one of them rejected by `filter::check`
        // (`reads_all_filtered`) — never more than one of the three, never
        // none. Segment-level drops are intentionally excluded: a read can
        // shed several segments and still survive, so per-segment counts
        // don't belong in a read-level invariant. A future workflow path
        // that adds a `continue`/early return without bumping exactly one of
        // the three read-level counters would silently under/over-count the
        // summary; this catches that in debug builds instead of shipping a
        // quietly-wrong report.
        debug_assert_eq!(
            reads_with_output + reads_trimmed_to_nothing + reads_all_filtered,
            input_reads,
            "every input read must be exactly one of: produced output, trimmed to \
             nothing, or had every segment filtered"
        );

        Stats {
            input_reads,
            output_reads: self.output_reads.load(Ordering::Relaxed),
            input_bases: self.input_bases.load(Ordering::Relaxed),
            output_bases: self.output_bases.load(Ordering::Relaxed),
            malformed_tag_reads,
            reads_trimmed_to_nothing,
            reads_all_filtered,
            segments_dropped_short,
            segments_dropped_long,
            segments_dropped_low_qual,
            segments_dropped_high_qual,
            segments_dropped_gc,
        }
    }
}

/// Filter produced segments and update segment- and read-level counters for all
/// workflows. `seq` and `qual` contain the complete input read and `produced`
/// contains the ranges to evaluate. For each surviving segment, `render`
/// receives `(idx, total, start, end)`. A render error stops processing before
/// the read-level outcome counter is updated.
pub(crate) fn process_read_segments<Rn>(
    produced: &[(usize, usize)],
    seq: &[u8],
    qual: &[u8],
    filter_cfg: &FilterConfig,
    counters: &Counters,
    mut render: Rn,
) -> anyhow::Result<()>
where
    Rn: FnMut(usize, usize, usize, usize) -> anyhow::Result<()>,
{
    let total = produced.len();
    let mut survived = 0usize;
    for (idx, &(s, e)) in produced.iter().enumerate() {
        if let Some(reason) = crate::filter::check(&seq[s..e], &qual[s..e], filter_cfg) {
            counters.record_segment_drop(reason);
            continue;
        }
        render(idx, total, s, e)?;
        counters.output_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .output_bases
            .fetch_add((e - s) as u64, Ordering::Relaxed);
        survived += 1;
    }
    if produced.is_empty() {
        counters
            .reads_trimmed_to_nothing
            .fetch_add(1, Ordering::Relaxed);
    } else if survived == 0 {
        counters.reads_all_filtered.fetch_add(1, Ordering::Relaxed);
    } else {
        counters.reads_with_output.fetch_add(1, Ordering::Relaxed);
    }
    Ok(())
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Stats {
    pub input_reads: u64,
    pub output_reads: u64,
    /// Sum of SEQ lengths (bases) across every input read.
    pub input_bases: u64,
    /// Sum of surviving segment lengths (bases) actually written to output.
    pub output_bases: u64,
    /// Reads carrying a known per-base kinetics tag (ip/pw/…) whose array length
    /// did not match the sequence length — malformed and left untouched. Surfaced
    /// as a run-level advisory; not an error.
    pub malformed_tag_reads: u64,
    /// Read-level: input reads that produced zero segments at all (empty
    /// read, fully consumed by adapter trimming, or an over-crop) —
    /// `trim::apply` returned no intervals, so the per-segment filter loop
    /// never ran.
    pub reads_trimmed_to_nothing: u64,
    /// Read-level: input reads that produced at least one segment, but every
    /// one of them was rejected by post-trim `filter::check`.
    pub reads_all_filtered: u64,
    /// Segment-level: segments dropped by post-trim `filter::check` for being
    /// shorter than `min_length` (including empty segments).
    pub segments_dropped_short: u64,
    /// Segment-level: segments dropped by post-trim `filter::check` for exceeding `max_length`.
    pub segments_dropped_long: u64,
    /// Segment-level: segments dropped by post-trim `filter::check` for quality below `min_qual`.
    pub segments_dropped_low_qual: u64,
    /// Segment-level: segments dropped by post-trim `filter::check` for quality above `max_qual`.
    pub segments_dropped_high_qual: u64,
    /// Segment-level: segments dropped by post-trim `filter::check` for GC fraction
    /// outside `[min_gc, max_gc]`.
    pub segments_dropped_gc: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batches_stop_at_weight_or_record_limit() {
        let by_weight: Vec<Vec<usize>> =
            Batches::new(vec![200_000usize; 5].into_iter(), |n: &usize| *n).collect();
        assert_eq!(by_weight.iter().map(Vec::len).collect::<Vec<_>>(), [3, 2]);

        let bam: Vec<Vec<usize>> =
            Batches::bam(vec![1usize; 17].into_iter(), |n: &usize| *n).collect();
        assert_eq!(
            bam.iter().map(Vec::len).collect::<Vec<_>>(),
            [4, 4, 4, 4, 1]
        );
    }

    /// The three-way read-level invariant (`reads_with_output +
    /// reads_trimmed_to_nothing + reads_all_filtered == input_reads`) models
    /// three reads: (a) a read with 2 surviving segments (`reads_with_output`),
    /// (b) a read with 0 survivors whose 2 produced segments were both
    /// `TooShort` (`reads_all_filtered` — segments were produced, none
    /// survived), and (c) an empty read with no segments produced at all
    /// (`reads_trimmed_to_nothing` — distinct from (b): the per-segment
    /// filter loop never even ran). Segment level: `output_reads` bumped
    /// twice (the 2 survivors of read a); 2 short drops (both from read b).
    #[test]
    fn three_way_read_counters_hold_the_invariant() {
        let counters = Counters::default();

        // (a) a read that splits into 2 surviving segments.
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters.output_reads.fetch_add(2, Ordering::Relaxed);
        counters.reads_with_output.fetch_add(1, Ordering::Relaxed);

        // (b) a read whose 2 produced segments are both filtered TooShort ->
        // reads_all_filtered (segments were produced, but none survived).
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters.record_segment_drop(DropReason::TooShort);
        counters.record_segment_drop(DropReason::TooShort);
        counters.reads_all_filtered.fetch_add(1, Ordering::Relaxed);

        // (c) an empty input read: trim::apply produces no segments at all ->
        // reads_trimmed_to_nothing (no segment-level drop is recorded, since
        // the per-segment filter loop never runs).
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .reads_trimmed_to_nothing
            .fetch_add(1, Ordering::Relaxed);

        let stats = counters.snapshot(0);

        assert_eq!(stats.input_reads, 3);
        assert_eq!(stats.output_reads, 2);
        assert_eq!(
            stats.reads_all_filtered, 1,
            "read b produced segments, but every one was filtered"
        );
        assert_eq!(
            stats.reads_trimmed_to_nothing, 1,
            "read c produced no segments at all"
        );
        assert_eq!(counters.reads_with_output.load(Ordering::Relaxed), 1);
        // reads_with_output + reads_trimmed_to_nothing + reads_all_filtered ==
        // input_reads holds (also asserted internally by `snapshot`'s
        // debug_assert_eq!).
        assert_eq!(
            counters.reads_with_output.load(Ordering::Relaxed)
                + stats.reads_trimmed_to_nothing
                + stats.reads_all_filtered,
            stats.input_reads
        );
        assert_eq!(stats.segments_dropped_short, 2);
        assert_eq!(stats.segments_dropped_long, 0);
        assert_eq!(stats.segments_dropped_low_qual, 0);
        assert_eq!(stats.segments_dropped_high_qual, 0);
        assert_eq!(stats.segments_dropped_gc, 0);
    }

    /// Cover all read-level outcomes and the corresponding render arguments.
    #[test]
    fn process_read_segments_dispatches_and_counts_all_three_outcomes() {
        let filter_cfg = FilterConfig {
            min_length: 3,
            max_length: usize::MAX,
            min_qual: 0.0,
            max_qual: 1000.0,
            min_gc: None,
            max_gc: None,
            qual_mode: crate::qual::QualMode::Mean,
        };

        // Trimmed to nothing: no produced intervals at all -> render never
        // called, reads_trimmed_to_nothing bumped, no segment-level drop.
        {
            let counters = Counters::default();
            let mut calls: Vec<(usize, usize, usize, usize)> = Vec::new();
            process_read_segments(&[], b"", b"", &filter_cfg, &counters, |idx, total, s, e| {
                calls.push((idx, total, s, e));
                Ok(())
            })
            .unwrap();
            assert!(calls.is_empty());
            assert_eq!(counters.reads_trimmed_to_nothing.load(Ordering::Relaxed), 1);
            assert_eq!(counters.reads_all_filtered.load(Ordering::Relaxed), 0);
            assert_eq!(counters.reads_with_output.load(Ordering::Relaxed), 0);
            assert_eq!(counters.segments_dropped_short.load(Ordering::Relaxed), 0);
        }

        // All filtered: one produced segment, too short to pass -> render
        // never called for it, reads_all_filtered bumped, one segment drop.
        {
            let counters = Counters::default();
            let seq = b"AA";
            let qual = b"II";
            let mut calls: Vec<(usize, usize, usize, usize)> = Vec::new();
            process_read_segments(
                &[(0, 2)],
                seq,
                qual,
                &filter_cfg,
                &counters,
                |idx, total, s, e| {
                    calls.push((idx, total, s, e));
                    Ok(())
                },
            )
            .unwrap();
            assert!(calls.is_empty());
            assert_eq!(counters.reads_trimmed_to_nothing.load(Ordering::Relaxed), 0);
            assert_eq!(counters.reads_all_filtered.load(Ordering::Relaxed), 1);
            assert_eq!(counters.reads_with_output.load(Ordering::Relaxed), 0);
            assert_eq!(counters.segments_dropped_short.load(Ordering::Relaxed), 1);
        }

        // With output: two produced segments, both long enough -> render
        // called once per survivor with the correct (idx, total, s, e).
        {
            let counters = Counters::default();
            let seq = b"AAAAAA";
            let qual = b"IIIIII";
            let mut calls: Vec<(usize, usize, usize, usize)> = Vec::new();
            process_read_segments(
                &[(0, 3), (3, 6)],
                seq,
                qual,
                &filter_cfg,
                &counters,
                |idx, total, s, e| {
                    calls.push((idx, total, s, e));
                    Ok(())
                },
            )
            .unwrap();
            assert_eq!(calls, vec![(0, 2, 0, 3), (1, 2, 3, 6)]);
            assert_eq!(counters.reads_trimmed_to_nothing.load(Ordering::Relaxed), 0);
            assert_eq!(counters.reads_all_filtered.load(Ordering::Relaxed), 0);
            assert_eq!(counters.reads_with_output.load(Ordering::Relaxed), 1);
            assert_eq!(counters.output_reads.load(Ordering::Relaxed), 2);
            assert_eq!(counters.output_bases.load(Ordering::Relaxed), 6);
        }
    }
}
