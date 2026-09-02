pub(crate) mod bam;
mod fastq;

use std::collections::BTreeMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use rayon::prelude::*;

pub use bam::{reconstruct_mods, reconstruct_record, run_bam, run_bam_to_fastq, run_raw_bam};
pub use fastq::{run_fastq, run_fastq_seq};

use crate::config::Config;
use crate::filter::{DropReason, FilterConfig};

/// Batch sizing for one parallel workflow. `target_weight` bounds the summed
/// record weight (bases) per batch and `max_items` bounds the record count, so a
/// batch is large enough to amortize scheduling and channel costs and small
/// enough to balance unusually long reads across workers. `queue_per_worker`
/// sizes the bounded channel to the writer in batches per render worker.
#[derive(Clone, Copy)]
pub(crate) struct BatchPolicy {
    target_weight: usize,
    max_items: usize,
    queue_per_worker: usize,
}

/// FASTQ batches: owned records that render to a small buffer each.
pub(crate) const FASTQ_BATCH: BatchPolicy = BatchPolicy {
    target_weight: 512 * 1024,
    max_items: 32,
    queue_per_worker: 4,
};

/// BAM batches: records that decode to large owned buffers.
pub(crate) const BAM_BATCH: BatchPolicy = BatchPolicy {
    target_weight: 256 * 1024,
    max_items: 4,
    queue_per_worker: 1,
};

/// Groups an iterator's items into batches bounded by a weight sum and an item
/// count.
pub(crate) struct Batches<I, F> {
    records: I,
    weight: F,
    policy: BatchPolicy,
}

impl<I, F> Batches<I, F> {
    pub(crate) fn new(records: I, weight: F, policy: BatchPolicy) -> Self {
        Self {
            records,
            weight,
            policy,
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
        let mut batch = Vec::with_capacity(self.policy.max_items);
        let mut bases = 0usize;
        while batch.len() < self.policy.max_items && bases < self.policy.target_weight {
            let Some(record) = self.records.next() else {
                break;
            };
            bases = bases.saturating_add((self.weight)(&record));
            batch.push(record);
        }
        (!batch.is_empty()).then_some(batch)
    }
}

/// The output of rendering one input record.
pub(crate) struct Rendered<T> {
    /// Output items in the order they are written.
    pub items: Vec<T>,
    /// Whether the record carried a known per-base tag whose length disagrees
    /// with the sequence length.
    pub malformed_tags: bool,
}

/// The render-pool size for a run: the settled budget, or the thread count when
/// no budget was settled.
pub(crate) fn render_pool_size(cfg: &Config) -> usize {
    if cfg.render_workers >= 1 {
        cfg.render_workers
    } else {
        cfg.threads.max(1)
    }
}

/// Keeps the first error of a parallel run and raises the shared abort flag.
struct FirstError<E> {
    slot: Mutex<Option<E>>,
}

impl<E> FirstError<E> {
    fn new() -> Self {
        Self {
            slot: Mutex::new(None),
        }
    }

    fn record(&self, error: E, aborted: &AtomicBool) {
        aborted.store(true, Ordering::Relaxed);
        let mut slot = self.slot.lock().unwrap();
        if slot.is_none() {
            *slot = Some(error);
        }
    }

    fn take(&self) -> Option<E> {
        self.slot.lock().unwrap().take()
    }
}

/// A record stream that ends after its first `Err` (which is still delivered)
/// or once the run's abort flag is raised. The source is checked before every
/// poll, so it is never polled again after either event.
struct FuseOnError<'a, I> {
    inner: I,
    done: bool,
    aborted: &'a AtomicBool,
}

impl<I, R> Iterator for FuseOnError<'_, I>
where
    I: Iterator<Item = anyhow::Result<R>>,
{
    type Item = anyhow::Result<R>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done || self.aborted.load(Ordering::Relaxed) {
            return None;
        }
        let item = self.inner.next()?;
        self.done = item.is_err();
        Some(item)
    }
}

/// The parallel driver shared by every multithreaded workflow. Records are
/// batched under `policy`, rendered on a rayon pool of `render_pool_size(cfg)`
/// threads, and written by one thread. The record stream is fused on its first
/// `Err` and stops being read once any render or write error is recorded, so a
/// failing run neither re-polls a reader after an I/O error nor processes the
/// rest of the input. With `cfg.ordered` the writer emits batches in input
/// order; otherwise in completion order.
///
/// The read-level counters are updated inside `render` by
/// `process_read_segments`; this driver counts input reads and bases only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_parallel<R, T, S, Weight, Render, WriteOne>(
    records: impl Iterator<Item = anyhow::Result<R>> + Send,
    policy: BatchPolicy,
    weight: Weight,
    cfg: &Config,
    sink: &mut S,
    render: Render,
    write_one: WriteOne,
    counters: &Counters,
) -> anyhow::Result<Stats>
where
    R: Send,
    T: Send,
    S: Send,
    Weight: Fn(&R) -> usize + Sync,
    Render: Fn(R, &Config) -> anyhow::Result<Rendered<T>> + Sync,
    WriteOne: Fn(&mut S, &T) -> std::io::Result<()> + Send,
{
    let render_workers = render_pool_size(cfg);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(render_workers)
        .build()?;
    let queue = (render_workers * policy.queue_per_worker).max(2);
    let (tx, rx) = crossbeam_channel::bounded::<(usize, Vec<T>)>(queue);
    let ordered = cfg.ordered;
    let aborted = AtomicBool::new(false);
    let malformed = AtomicU64::new(0);
    let render_err: FirstError<anyhow::Error> = FirstError::new();
    let write_err: FirstError<std::io::Error> = FirstError::new();

    let aborted_ref = &aborted;
    let records = FuseOnError {
        inner: records,
        done: false,
        aborted: aborted_ref,
    };

    // The writer is a scoped OS thread; rendering runs on the local pool via
    // `pool.install`, so the nested `par_bridge` is bounded by `-t`. The writer
    // keeps draining after an error so a producer blocked on `tx.send` cannot
    // deadlock.
    std::thread::scope(|s| {
        let write_err = &write_err;
        s.spawn(move || {
            let mut next = 0usize;
            let mut pending: BTreeMap<usize, Vec<T>> = BTreeMap::new();
            let mut errored = false;
            let mut write_batch = |batch: &[T]| -> bool {
                for item in batch {
                    if let Err(e) = write_one(sink, item) {
                        write_err.record(e, aborted_ref);
                        return false;
                    }
                }
                true
            };
            for (idx, batch) in rx.iter() {
                if errored {
                    continue;
                }
                if ordered {
                    pending.insert(idx, batch);
                    while let Some(batch) = pending.remove(&next) {
                        if !write_batch(&batch) {
                            errored = true;
                            break;
                        }
                        next += 1;
                    }
                } else if !write_batch(&batch) {
                    errored = true;
                }
            }
        });

        pool.install(|| {
            let weight_of = |rec: &anyhow::Result<R>| rec.as_ref().map_or(0, &weight);
            Batches::new(records, weight_of, policy)
                .enumerate()
                .par_bridge()
                .for_each(|(idx, batch)| {
                    let mut out = Vec::with_capacity(batch.len());
                    let mut input_reads = 0u64;
                    let mut input_bases = 0u64;
                    let mut malformed_reads = 0u64;
                    for rec in batch {
                        if aborted.load(Ordering::Relaxed) {
                            break;
                        }
                        let rec = match rec {
                            Ok(r) => r,
                            Err(e) => {
                                render_err.record(e, &aborted);
                                break;
                            },
                        };
                        input_reads += 1;
                        input_bases += weight(&rec) as u64;
                        match render(rec, cfg) {
                            Ok(rendered) => {
                                malformed_reads += u64::from(rendered.malformed_tags);
                                out.extend(rendered.items);
                            },
                            Err(e) => {
                                render_err.record(e, &aborted);
                                break;
                            },
                        }
                    }
                    counters
                        .input_reads
                        .fetch_add(input_reads, Ordering::Relaxed);
                    counters
                        .input_bases
                        .fetch_add(input_bases, Ordering::Relaxed);
                    malformed.fetch_add(malformed_reads, Ordering::Relaxed);
                    // Every batch is sent, empty ones included, so the ordered
                    // writer can advance past it. A closed channel means the
                    // writer is gone; nothing more can be written.
                    if tx.send((idx, out)).is_err() {
                        aborted.store(true, Ordering::Relaxed);
                    }
                });
        });
        drop(tx);
    });

    if let Some(e) = render_err.take() {
        return Err(e);
    }
    if let Some(e) = write_err.take() {
        return Err(e.into());
    }
    Ok(counters.snapshot(malformed.load(Ordering::Relaxed)))
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
    /// Input reads whose `MM`/`ML`/`MN` block was malformed (an `MN` that
    /// disagrees with the sequence length, an `ML` whose length disagrees with
    /// `MM`, an `MM` that does not parse to its end, or a non-`B:C` `ML`) and was
    /// therefore removed from the output record.
    pub malformed_mod_reads: AtomicU64,
    /// Input reads that produced at least one surviving output segment,
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
    /// several pieces are each judged independently). These are NOT part of
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

        // Every input read lands in exactly one of the three read-level buckets: it
        // produced surviving segments, produced none at all, or produced some and
        // lost them all to `filter::check`. Segment-level drops are excluded, since a
        // read can shed segments and still survive. Catches a future early return
        // that forgets to bump one of the three.
        debug_assert_eq!(
            reads_with_output + reads_trimmed_to_nothing + reads_all_filtered,
            input_reads,
            "every input read must be exactly one of: produced output, trimmed to \
             nothing, or had every segment filtered"
        );

        Stats {
            input_reads,
            output_reads: self.output_reads.load(Ordering::Relaxed),
            malformed_mod_reads: self.malformed_mod_reads.load(Ordering::Relaxed),
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
/// A trace-level span naming the read, so the segment decisions logged while
/// processing it are attributable without repeating the name on every event.
///
/// Returns a disabled span when trace is off, which costs a level check rather
/// than any formatting, keeping the hot path unaffected at the default level.
pub(crate) fn read_span(name: &[u8]) -> tracing::Span {
    if tracing::enabled!(tracing::Level::TRACE) {
        tracing::trace_span!("read", name = %String::from_utf8_lossy(name))
    } else {
        tracing::Span::none()
    }
}

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
            // The per-segment verdict is the answer to "why is this read missing
            // from my output", which no run-level counter can give.
            tracing::trace!(
                segment = idx + 1,
                of = total,
                start = s,
                end = e,
                len = e - s,
                reason = reason.label(),
                "Segment dropped"
            );
            counters.record_segment_drop(reason);
            continue;
        }
        tracing::trace!(
            segment = idx + 1,
            of = total,
            start = s,
            end = e,
            len = e - s,
            "Segment kept"
        );
        render(idx, total, s, e)?;
        counters.output_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .output_bases
            .fetch_add((e - s) as u64, Ordering::Relaxed);
        survived += 1;
    }
    if produced.is_empty() {
        tracing::trace!("Read produced no segments");
        counters
            .reads_trimmed_to_nothing
            .fetch_add(1, Ordering::Relaxed);
    } else if survived == 0 {
        tracing::trace!(produced = total, "Every segment filtered");
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
    /// Reads carrying a known per-base kinetics tag (ip/pw/...) whose array length
    /// did not match the sequence length: malformed and left untouched. Surfaced
    /// as a run-level advisory; not an error.
    pub malformed_tag_reads: u64,
    /// Input reads whose modification block was malformed and removed; see
    /// `Counters::malformed_mod_reads`.
    pub malformed_mod_reads: u64,
    /// Read-level: input reads that produced zero segments at all (empty
    /// read, fully consumed by adapter trimming, or an over-crop).
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
            Batches::new(vec![200_000usize; 5].into_iter(), |n: &usize| *n, FASTQ_BATCH).collect();
        assert_eq!(by_weight.iter().map(Vec::len).collect::<Vec<_>>(), [3, 2]);

        let bam: Vec<Vec<usize>> =
            Batches::new(vec![1usize; 17].into_iter(), |n: &usize| *n, BAM_BATCH).collect();
        assert_eq!(
            bam.iter().map(Vec::len).collect::<Vec<_>>(),
            [4, 4, 4, 4, 1]
        );
    }

    /// The three-way read-level invariant (`reads_with_output +
    /// reads_trimmed_to_nothing + reads_all_filtered == input_reads`) over three
    /// reads: (a) 2 surviving segments, (b) 2 produced segments both `TooShort`,
    /// and (c) an empty read producing none. (c) differs from (b) in that its
    /// per-segment filter loop never runs.
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

    fn driver_cfg(threads: usize, ordered: bool) -> Config {
        let path = std::path::Path::new("/dev/null");
        let mut cfg = crate::cli::config_for_test_threads(path, path, 0, 0, threads);
        cfg.ordered = ordered;
        cfg
    }

    /// Rendering odd items slowly forces completion order to differ from input
    /// order, so the ordered writer is exercised rather than trivially satisfied.
    fn run_driver(ordered: bool) -> Vec<usize> {
        let cfg = driver_cfg(4, ordered);
        let mut sink: Vec<usize> = Vec::new();
        let counters = Counters::default();
        run_parallel(
            (0..200usize).map(anyhow::Ok),
            BatchPolicy {
                target_weight: 1,
                max_items: 1,
                queue_per_worker: 1,
            },
            |_: &usize| 1,
            &cfg,
            &mut sink,
            |n, _cfg| {
                if n % 2 == 1 {
                    std::thread::sleep(std::time::Duration::from_micros(200));
                }
                counters.reads_with_output.fetch_add(1, Ordering::Relaxed);
                Ok(Rendered {
                    items: vec![n],
                    malformed_tags: false,
                })
            },
            |sink, n: &usize| {
                sink.push(*n);
                Ok(())
            },
            &counters,
        )
        .unwrap();
        sink
    }

    #[test]
    fn ordered_driver_writes_in_input_order() {
        let out = run_driver(true);
        assert_eq!(out, (0..200).collect::<Vec<_>>());
    }

    #[test]
    fn unordered_driver_writes_every_item_once() {
        let mut out = run_driver(false);
        out.sort_unstable();
        assert_eq!(out, (0..200).collect::<Vec<_>>());
    }

    /// The record stream is fused on its first `Err`: the source is never
    /// polled again, so a reader left in an inconsistent state after an I/O
    /// error cannot panic.
    #[test]
    fn driver_stops_polling_the_source_after_its_first_error() {
        struct PoisonAfterError {
            n: usize,
            yielded_error: bool,
        }
        impl Iterator for PoisonAfterError {
            type Item = anyhow::Result<usize>;
            fn next(&mut self) -> Option<Self::Item> {
                assert!(!self.yielded_error, "source polled after it returned Err");
                if self.n == 50 {
                    self.yielded_error = true;
                    return Some(Err(anyhow::anyhow!("incomplete stream")));
                }
                self.n += 1;
                Some(Ok(self.n))
            }
        }
        let cfg = driver_cfg(4, false);
        let mut sink: Vec<usize> = Vec::new();
        let counters = Counters::default();
        let res = run_parallel(
            PoisonAfterError {
                n: 0,
                yielded_error: false,
            },
            FASTQ_BATCH,
            |_: &usize| 1,
            &cfg,
            &mut sink,
            |n, _cfg| {
                Ok(Rendered {
                    items: vec![n],
                    malformed_tags: false,
                })
            },
            |sink, n: &usize| {
                sink.push(*n);
                Ok(())
            },
            &counters,
        );
        assert_eq!(res.unwrap_err().to_string(), "incomplete stream");
    }

    /// A render error stops the run: records after the failing one are not
    /// rendered, so a failing run does not process the rest of its input.
    #[test]
    fn driver_stops_rendering_after_the_first_render_error() {
        use std::sync::atomic::AtomicUsize;

        let cfg = driver_cfg(2, false);
        let mut sink: Vec<usize> = Vec::new();
        let counters = Counters::default();
        let rendered = AtomicUsize::new(0);
        let res = run_parallel(
            (0..100_000usize).map(anyhow::Ok),
            BatchPolicy {
                target_weight: 1,
                max_items: 1,
                queue_per_worker: 1,
            },
            |_: &usize| 1,
            &cfg,
            &mut sink,
            |n, _cfg| {
                rendered.fetch_add(1, Ordering::Relaxed);
                if n == 10 {
                    anyhow::bail!("record 10 is malformed");
                }
                Ok(Rendered {
                    items: vec![n],
                    malformed_tags: false,
                })
            },
            |sink, n: &usize| {
                sink.push(*n);
                Ok(())
            },
            &counters,
        );
        assert!(res.is_err());
        assert!(
            rendered.load(Ordering::Relaxed) < 1_000,
            "rendering continued long after the first error: {} records",
            rendered.load(Ordering::Relaxed)
        );
    }
}
