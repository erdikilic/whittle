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
    /// Per-reason drop counters. Every input read that produces no output is
    /// counted under exactly one of these — the five `filter::check` reasons,
    /// plus `dropped_trimmed` for a read that passed the filter but whose
    /// trim intervals were empty (no window, or every segment below
    /// `min_length`).
    pub dropped_short: AtomicU64,
    pub dropped_long: AtomicU64,
    pub dropped_low_qual: AtomicU64,
    pub dropped_high_qual: AtomicU64,
    pub dropped_gc: AtomicU64,
    pub dropped_trimmed: AtomicU64,
}

impl Counters {
    /// Bump the counter matching a `filter::check` failure reason.
    pub fn record_filter_drop(&self, reason: DropReason) {
        let counter = match reason {
            DropReason::TooShort => &self.dropped_short,
            DropReason::TooLong => &self.dropped_long,
            DropReason::LowQuality => &self.dropped_low_qual,
            DropReason::HighQuality => &self.dropped_high_qual,
            DropReason::Gc => &self.dropped_gc,
        };
        counter.fetch_add(1, Ordering::Relaxed);
    }

    /// Snapshot every counter into a `Stats` for end-of-run reporting.
    /// `malformed_tag_reads` is threaded through separately: only the BAM
    /// paths track it, and the parallel BAM path accumulates it in its own
    /// local atomic rather than in `Counters`.
    pub fn snapshot(&self, malformed_tag_reads: u64) -> Stats {
        Stats {
            input_reads: self.input_reads.load(Ordering::Relaxed),
            output_reads: self.output_reads.load(Ordering::Relaxed),
            input_bases: self.input_bases.load(Ordering::Relaxed),
            output_bases: self.output_bases.load(Ordering::Relaxed),
            malformed_tag_reads,
            dropped_short: self.dropped_short.load(Ordering::Relaxed),
            dropped_long: self.dropped_long.load(Ordering::Relaxed),
            dropped_low_qual: self.dropped_low_qual.load(Ordering::Relaxed),
            dropped_high_qual: self.dropped_high_qual.load(Ordering::Relaxed),
            dropped_gc: self.dropped_gc.load(Ordering::Relaxed),
            dropped_trimmed: self.dropped_trimmed.load(Ordering::Relaxed),
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
    /// Input reads dropped by `filter::check` for being shorter than `min_length`
    /// (including empty reads).
    pub dropped_short: u64,
    /// Input reads dropped by `filter::check` for exceeding `max_length`.
    pub dropped_long: u64,
    /// Input reads dropped by `filter::check` for read quality below `min_qual`.
    pub dropped_low_qual: u64,
    /// Input reads dropped by `filter::check` for read quality above `max_qual`.
    pub dropped_high_qual: u64,
    /// Input reads dropped by `filter::check` for GC fraction outside `[min_gc, max_gc]`.
    pub dropped_gc: u64,
    /// Input reads that passed the filter but produced zero output segments
    /// after trimming (no interval, or every segment below `min_length`).
    pub dropped_trimmed: u64,
}
