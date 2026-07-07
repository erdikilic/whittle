mod bam;
mod fastq;

use std::sync::atomic::{AtomicU64, Ordering};

pub use bam::{reconstruct_mods, reconstruct_record, run_bam, run_bam_to_fastq};
pub use fastq::{run_fastq, run_fastq_seq};

use crate::filter::DropReason;

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
    /// read is accounted for by exactly one of "no output" or "produced
    /// output" (the read-level half of the two-level counter model).
    pub reads_with_output: AtomicU64,
    /// Input reads that produced **zero** surviving segments: an empty read,
    /// a read fully consumed by adapter trimming, or a read whose every
    /// produced segment was rejected by `filter::check`. Read-level, paired
    /// with `reads_with_output` in the invariant below.
    pub reads_no_output: AtomicU64,
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
    pub fn record_filter_drop(&self, reason: DropReason) {
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
        let reads_no_output = self.reads_no_output.load(Ordering::Relaxed);
        let segments_dropped_short = self.segments_dropped_short.load(Ordering::Relaxed);
        let segments_dropped_long = self.segments_dropped_long.load(Ordering::Relaxed);
        let segments_dropped_low_qual = self.segments_dropped_low_qual.load(Ordering::Relaxed);
        let segments_dropped_high_qual = self.segments_dropped_high_qual.load(Ordering::Relaxed);
        let segments_dropped_gc = self.segments_dropped_gc.load(Ordering::Relaxed);

        // Every input read either produced at least one surviving output
        // segment (counted once in `reads_with_output`, regardless of how
        // many segments it split into) or produced none at all (counted once
        // in `reads_no_output`) — never both, never neither. Segment-level
        // drops are intentionally excluded: a read can shed several segments
        // and still survive, so per-segment counts don't belong in a
        // read-level invariant. A future workflow path that adds a
        // `continue`/early return without bumping one of the two read-level
        // counters would silently under/over-count the summary; this catches
        // that in debug builds instead of shipping a quietly-wrong report.
        debug_assert_eq!(
            reads_with_output + reads_no_output,
            input_reads,
            "every input read must be either counted as no-output or have produced output"
        );

        Stats {
            input_reads,
            output_reads: self.output_reads.load(Ordering::Relaxed),
            input_bases: self.input_bases.load(Ordering::Relaxed),
            output_bases: self.output_bases.load(Ordering::Relaxed),
            malformed_tag_reads,
            reads_no_output,
            segments_dropped_short,
            segments_dropped_long,
            segments_dropped_low_qual,
            segments_dropped_high_qual,
            segments_dropped_gc,
        }
    }
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
    /// Read-level: input reads that produced zero surviving segments (empty
    /// read, fully consumed by adapter trimming, or every produced segment
    /// filtered post-trim).
    pub reads_no_output: u64,
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

    /// The two-level invariant (`reads_with_output + reads_no_output ==
    /// input_reads`) models three reads: (a) a read with 2 surviving segments,
    /// (b) a read with 0 survivors whose 2 produced segments were both
    /// `TooShort`, and (c) an empty read (no segments produced at all). Read
    /// level: 3 input reads, 1 with output, 2 with no output. Segment level:
    /// `output_reads` bumped twice (the 2 survivors of read a); 2 short drops.
    #[test]
    fn two_level_counters_hold_the_invariant() {
        let counters = Counters::default();

        // (a) a read that splits into 2 surviving segments.
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters.output_reads.fetch_add(2, Ordering::Relaxed);
        counters.reads_with_output.fetch_add(1, Ordering::Relaxed);

        // (b) a read whose 2 produced segments are both filtered TooShort.
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters.record_filter_drop(DropReason::TooShort);
        counters.record_filter_drop(DropReason::TooShort);
        counters.reads_no_output.fetch_add(1, Ordering::Relaxed);

        // (c) an empty input read: trim::apply produces no segments at all.
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters.reads_no_output.fetch_add(1, Ordering::Relaxed);

        let stats = counters.snapshot(0);

        assert_eq!(stats.input_reads, 3);
        assert_eq!(stats.output_reads, 2);
        assert_eq!(
            stats.reads_no_output, 2,
            "reads b and c both produced no surviving output"
        );
        assert_eq!(counters.reads_with_output.load(Ordering::Relaxed), 1);
        // reads_with_output + reads_no_output == input_reads holds (also
        // asserted internally by `snapshot`'s debug_assert_eq!).
        assert_eq!(
            counters.reads_with_output.load(Ordering::Relaxed) + stats.reads_no_output,
            stats.input_reads
        );
        assert_eq!(stats.segments_dropped_short, 2);
        assert_eq!(stats.segments_dropped_long, 0);
        assert_eq!(stats.segments_dropped_low_qual, 0);
        assert_eq!(stats.segments_dropped_high_qual, 0);
        assert_eq!(stats.segments_dropped_gc, 0);
    }
}
