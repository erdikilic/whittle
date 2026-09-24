//! Workflow drivers shared by the FASTQ and BAM paths: batching, the parallel driver, run counters and per-read segment filtering.

pub(crate) mod bam;
mod fastq;
#[cfg(feature = "paraseq")]
mod paraseq;
pub(crate) mod reject;

use std::collections::BTreeMap;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock};

use rayon::prelude::*;

pub(crate) use bam::run_bam_to_fastq;
pub use bam::run_raw_bam;

/// Wraps a record stream so that only reads `keep` accepts pass on. A rejected
/// read is counted here as input and as tag-filtered, since no workflow sees
/// it, and handed to `reject`, which writes it to the rejected output when one
/// is open; an evaluation error names the read.
pub(crate) fn filter_by_tags<R, I, K, W, N, J>(
    records: I,
    counters: Arc<Counters>,
    keep: K,
    weight: W,
    name: N,
    reject: J,
) -> Box<dyn Iterator<Item = anyhow::Result<R>> + Send>
where
    R: Send + 'static,
    I: Iterator<Item = anyhow::Result<R>> + Send + 'static,
    K: Fn(&R) -> anyhow::Result<bool> + Send + 'static,
    W: Fn(&R) -> usize + Send + 'static,
    N: Fn(&R) -> String + Send + 'static,
    J: Fn(&R) -> anyhow::Result<()> + Send + 'static,
{
    Box::new(records.filter_map(move |item| {
        let rec = match item {
            Ok(rec) => rec,
            Err(e) => return Some(Err(e)),
        };
        match keep(&rec) {
            Ok(true) => Some(Ok(rec)),
            Ok(false) => {
                counters.input_reads.fetch_add(1, Ordering::Relaxed);
                counters
                    .input_bases
                    .fetch_add(weight(&rec) as u64, Ordering::Relaxed);
                counters.tag_filtered_reads.fetch_add(1, Ordering::Relaxed);
                match reject(&rec) {
                    Ok(()) => None,
                    Err(e) => Some(Err(e.context(format!("read {}", name(&rec))))),
                }
            },
            Err(e) => Some(Err(e.context(format!("read {}", name(&rec))))),
        }
    }))
}
pub(crate) use fastq::{run_fastq, tag_filtered_fastq};
#[cfg(feature = "paraseq")]
pub(crate) use paraseq::{run_fastq_paraseq, selected as paraseq_selected};

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

/// FASTQ batches: owned records rendered into a shared batch buffer.
pub(crate) const FASTQ_BATCH: BatchPolicy = BatchPolicy {
    target_weight: 2 * 1024 * 1024,
    max_items: 128,
    queue_per_worker: 2,
};

/// BAM batches: records that decode to large owned buffers.
pub(crate) const BAM_BATCH: BatchPolicy = BatchPolicy {
    target_weight: 1024 * 1024,
    max_items: 16,
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
    /// Wraps `records`, weighing each item with `weight`, under `policy`.
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

/// A sink of rendered text batches: plain bytes, or BGZF blocks the render
/// workers compressed.
pub(crate) trait BatchSink: Write + Send {
    /// The BGZF level when the sink takes blocks compressed on the render
    /// pool; `None` for a sink that takes the rendered bytes as they are.
    fn block_level(&self) -> Option<u8> {
        None
    }
}

impl BatchSink for Vec<u8> {}

impl BatchSink for crate::io::fastq::FastqOut {
    fn block_level(&self) -> Option<u8> {
        crate::io::fastq::FastqOut::block_level(self)
    }
}

/// Renders directly into a batch buffer and compresses block output on the pool.
pub(crate) fn run_bytes_parallel<R, W, Weight, Render>(
    records: impl Iterator<Item = anyhow::Result<R>> + Send,
    policy: BatchPolicy,
    weight: Weight,
    cfg: &Config,
    writer: &mut W,
    render: Render,
    counters: &Counters,
) -> anyhow::Result<Stats>
where
    R: Send,
    W: BatchSink,
    Weight: Fn(&R) -> usize + Sync,
    Render: Fn(R, &Config, &mut Vec<u8>) -> anyhow::Result<()> + Sync,
{
    let level = writer.block_level();
    run_parallel(
        records,
        policy,
        weight,
        cfg,
        writer,
        render,
        |bytes: Vec<u8>| match level {
            Some(level) => crate::io::fastq::encode_blocks(level, &bytes),
            None => Ok(bytes),
        },
        |writer, bytes: &Vec<u8>| writer.write_all(bytes),
        counters,
    )
}

/// Returns the render-pool size for a run: the settled budget, or the thread
/// count when no budget was settled.
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

/// Bounds how far batch hand-out runs ahead of the ordered writer: batch `idx`
/// is handed out only once fewer than `size` batches before it are unwritten,
/// so a slow batch holds back at most `size` batches behind it.
struct ReorderWindow {
    /// Batches written so far.
    written: Mutex<usize>,
    /// Signaled when `written` advances or the run aborts.
    advanced: Condvar,
    /// Batches allowed in flight.
    size: usize,
}

impl ReorderWindow {
    fn new(size: usize) -> Self {
        Self {
            written: Mutex::new(0),
            advanced: Condvar::new(),
            size,
        }
    }

    /// Blocks until batch `idx` fits in the window or the run has aborted.
    fn wait_for(&self, idx: usize, aborted: &AtomicBool) {
        let mut written = self.written.lock().unwrap();
        while idx >= *written + self.size && !aborted.load(Ordering::Relaxed) {
            written = self.advanced.wait(written).unwrap();
        }
    }

    /// Records that `written` batches have been written.
    fn advance(&self, written: usize) {
        *self.written.lock().unwrap() = written;
        self.advanced.notify_all();
    }

    /// Wakes every waiter after the abort flag is raised. The lock orders the
    /// flag before a waiter's next check.
    fn release(&self) {
        let _written = self.written.lock().unwrap();
        self.advanced.notify_all();
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

/// Runs the parallel driver shared by every multithreaded workflow. Records are
/// batched under `policy`, rendered on a rayon pool of `render_pool_size(cfg)`
/// threads, and written by one thread. The record stream is fused on its first
/// `Err` and stops being read once any render or write error is recorded, so a
/// failing run neither re-polls a reader after an I/O error nor processes the
/// rest of the input. With `cfg.ordered` the writer emits batches in input
/// order, and a batch is handed out only while a bounded window of batches
/// separates it from the next one to write; otherwise in completion order.
///
/// `render` appends each record's output to the batch accumulator; `pack`
/// turns the accumulator into the unit the writer takes, on
/// the pool, so a compressing sink has its blocks compressed by the render
/// workers. The read-level counters are updated inside `render` by
/// `process_read_segments`; this driver counts input reads and bases only.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_parallel<R, T, P, S, Weight, Render, Pack, WriteOne>(
    records: impl Iterator<Item = anyhow::Result<R>> + Send,
    policy: BatchPolicy,
    weight: Weight,
    cfg: &Config,
    sink: &mut S,
    render: Render,
    pack: Pack,
    write_one: WriteOne,
    counters: &Counters,
) -> anyhow::Result<Stats>
where
    R: Send,
    T: Default + Send,
    P: Send,
    S: Send,
    Weight: Fn(&R) -> usize + Sync,
    Render: Fn(R, &Config, &mut T) -> anyhow::Result<()> + Sync,
    Pack: Fn(T) -> std::io::Result<P> + Sync,
    WriteOne: Fn(&mut S, &P) -> std::io::Result<()> + Send,
{
    let render_workers = render_pool_size(cfg);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(render_workers)
        .build()?;
    let queue = (render_workers * policy.queue_per_worker).max(2);
    let (tx, rx) = std::sync::mpsc::sync_channel::<(usize, P)>(queue);
    let ordered = cfg.ordered;
    let aborted = AtomicBool::new(false);
    let render_err: FirstError<anyhow::Error> = FirstError::new();
    let write_err: FirstError<std::io::Error> = FirstError::new();
    let window = ReorderWindow::new(queue + render_workers);

    let aborted_ref = &aborted;
    let window_ref = &window;
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
            let mut pending: BTreeMap<usize, P> = BTreeMap::new();
            let mut errored = false;
            let mut write_batch = |batch: &P| -> bool {
                if let Err(e) = write_one(sink, batch) {
                    write_err.record(e, aborted_ref);
                    window_ref.release();
                    return false;
                }
                true
            };
            for (idx, batch) in rx.iter() {
                if errored {
                    continue;
                }
                if ordered {
                    pending.insert(idx, batch);
                    let written = next;
                    while let Some(batch) = pending.remove(&next) {
                        if !write_batch(&batch) {
                            errored = true;
                            break;
                        }
                        next += 1;
                    }
                    if next > written {
                        window_ref.advance(next);
                    }
                } else if !write_batch(&batch) {
                    errored = true;
                }
            }
        });

        pool.install(|| {
            let weight_of = |rec: &anyhow::Result<R>| rec.as_ref().map_or(0, &weight);
            let mut batches = Batches::new(records, weight_of, policy);
            let mut handed_out = 0usize;
            // Under `ordered`, the next batch is read only once it fits in the
            // reorder window, which bounds the reordered output in memory.
            let batches = std::iter::from_fn(move || {
                if ordered {
                    window_ref.wait_for(handed_out, aborted_ref);
                }
                let batch = batches.next()?;
                handed_out += 1;
                Some((handed_out - 1, batch))
            });
            let render_batch = |(idx, batch): (usize, Vec<anyhow::Result<R>>)| {
                let mut out = T::default();
                let mut input_reads = 0u64;
                let mut input_bases = 0u64;
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
                    if let Err(e) = render(rec, cfg, &mut out) {
                        render_err.record(e, &aborted);
                        break;
                    }
                }
                if aborted.load(Ordering::Relaxed) {
                    window.release();
                }
                counters
                    .input_reads
                    .fetch_add(input_reads, Ordering::Relaxed);
                counters
                    .input_bases
                    .fetch_add(input_bases, Ordering::Relaxed);
                let packed = match pack(out) {
                    Ok(packed) => packed,
                    Err(e) => {
                        render_err.record(e.into(), &aborted);
                        window.release();
                        return;
                    },
                };
                // Every batch is sent, empty ones included, so the ordered
                // writer can advance past it. A closed channel means the
                // writer is gone; nothing more can be written.
                if tx.send((idx, packed)).is_err() {
                    aborted.store(true, Ordering::Relaxed);
                    window.release();
                }
            };
            batches.par_bridge().for_each(render_batch);
        });
        drop(tx);
    });

    if let Some(e) = render_err.take() {
        return Err(e);
    }
    if let Some(e) = write_err.take() {
        return Err(e.into());
    }
    Ok(counters.snapshot())
}

/// Live, thread-shared counters read by the progress ticker and finalized into `Stats`.
#[derive(Default)]
pub struct Counters {
    /// Whether any FASTQ header carries SAM auxiliary fields.
    pub(crate) tagged_fastq: AtomicBool,
    /// Input reads consumed.
    pub input_reads: AtomicU64,
    /// Output segments written (one per surviving segment).
    pub output_reads: AtomicU64,
    /// Compressed bytes pulled from the input source.
    pub bytes_read: AtomicU64,
    /// Sum of SEQ lengths (bases) across every input read, regardless of
    /// whether it survives filtering/trimming.
    pub input_bases: AtomicU64,
    /// Sum of surviving segment lengths (bases) written to output.
    pub output_bases: AtomicU64,
    /// Input reads carrying a known per-base tag (`ip`, `pw`, ...) whose array
    /// length disagrees with the sequence length, or an `sa` coverage array
    /// whose runs do not sum to it. The tag is left untouched and the count is
    /// surfaced as a run-level advisory. Tracked by the BAM paths only.
    pub malformed_tag_reads: AtomicU64,
    /// Input reads whose `MM`/`ML`/`MN` block was malformed (an `MN` that
    /// disagrees with the sequence length, an `ML` whose length disagrees with
    /// `MM`, an `MM` that does not parse to its end or exceeds the available
    /// counting-base occurrences, or a non-`B:C` `ML`) and was
    /// therefore removed from the output record.
    pub malformed_mod_reads: AtomicU64,
    /// Trimmed input reads whose PacBio undo blobs (`ds`, `ls`) were removed,
    /// since they describe the untrimmed read.
    pub undo_tags_dropped_reads: AtomicU64,
    /// Input reads whose `bi` barcode tag could not be read as a barcode window
    /// (not a seven-element float array, or an empty, inverted or out-of-range
    /// window) and were left untrimmed by that stage.
    pub barcode_tag_malformed_reads: AtomicU64,
    /// Input reads with a recorded barcode span at which no configured or
    /// named catalog barcode was found; that span was left untrimmed.
    pub barcode_tag_unverified_reads: AtomicU64,
    /// Input reads rejected by `--tag-filter` before any workflow saw them.
    /// Read-level, part of the invariant below.
    pub tag_filtered_reads: AtomicU64,
    /// The `--rejected-output` channel, set before the first record is read when
    /// the flag is given; unset otherwise.
    pub(crate) rejects: OnceLock<reject::Rejects>,
    /// Input reads that produced at least one surviving output segment,
    /// bumped once per input read (not once per segment, unlike
    /// `output_reads`, which a `--split-quality` read can bump several times).
    /// Exists so `snapshot` can check that every input read is accounted for
    /// by exactly one of the three read-level outcomes.
    pub reads_with_output: AtomicU64,
    /// Input reads that produced zero segments: `trim::apply` returned no
    /// intervals for the per-segment filter (an empty read, a read fully
    /// consumed by adapter trimming, or an over-crop). Read-level, paired with
    /// `reads_with_output` and `reads_all_filtered` in the invariant below.
    pub reads_trimmed_to_nothing: AtomicU64,
    /// Input reads that produced at least one segment, but every one was
    /// rejected by post-trim `filter::check`. Read-level, paired with
    /// `reads_with_output` and `reads_trimmed_to_nothing` in the invariant
    /// below.
    pub reads_all_filtered: AtomicU64,
    /// Segments dropped as `TooShort`. This and the four counters below are
    /// segment-level: one bump per segment (not read) that `filter::check`
    /// rejects, by reason, post-trim. A single input read can contribute to
    /// more than one of these (e.g. a `--split-quality` read whose several pieces
    /// are each judged independently). They are not part of the read-level
    /// invariant.
    pub segments_dropped_short: AtomicU64,
    /// Segments dropped as `TooLong`.
    pub segments_dropped_long: AtomicU64,
    /// Segments dropped as `LowQuality`.
    pub segments_dropped_low_qual: AtomicU64,
    /// Segments dropped as `HighQuality`.
    pub segments_dropped_high_qual: AtomicU64,
    /// Segments dropped as `Gc`.
    pub segments_dropped_gc: AtomicU64,
}

impl Counters {
    /// Whether rejected records are written, so producers can skip rendering them.
    pub(crate) fn wants_rejects(&self) -> bool {
        self.rejects.get().is_some()
    }

    /// Queues a rejected record when `--rejected-output` is set.
    pub(crate) fn reject(&self, item: reject::RejectItem) -> anyhow::Result<()> {
        match self.rejects.get() {
            Some(rejects) => rejects.send(item),
            None => Ok(()),
        }
    }

    /// Bumps the segment-level counter matching a `filter::check` failure
    /// reason. Called once per rejected segment (post-trim), not per read.
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

    /// Snapshots every counter into a `Stats` for end-of-run reporting.
    pub fn snapshot(&self) -> Stats {
        let input_reads = self.input_reads.load(Ordering::Relaxed);
        let reads_with_output = self.reads_with_output.load(Ordering::Relaxed);
        let reads_trimmed_to_nothing = self.reads_trimmed_to_nothing.load(Ordering::Relaxed);
        let reads_all_filtered = self.reads_all_filtered.load(Ordering::Relaxed);
        let reads_tag_filtered = self.tag_filtered_reads.load(Ordering::Relaxed);
        let segments_dropped_short = self.segments_dropped_short.load(Ordering::Relaxed);
        let segments_dropped_long = self.segments_dropped_long.load(Ordering::Relaxed);
        let segments_dropped_low_qual = self.segments_dropped_low_qual.load(Ordering::Relaxed);
        let segments_dropped_high_qual = self.segments_dropped_high_qual.load(Ordering::Relaxed);
        let segments_dropped_gc = self.segments_dropped_gc.load(Ordering::Relaxed);

        // Every input read lands in exactly one of the four read-level buckets: it
        // was rejected by the tag filter, produced surviving segments, produced
        // none at all, or produced some and lost them all to `filter::check`.
        // Segment-level drops are excluded, since a read can shed segments and
        // still survive. The assertion catches an early return that skips one.
        debug_assert_eq!(
            reads_with_output + reads_trimmed_to_nothing + reads_all_filtered + reads_tag_filtered,
            input_reads,
            "Every input read must be exactly one of: tag filtered, produced output, \
             trimmed to nothing, or had every segment filtered"
        );

        Stats {
            input_reads,
            output_reads: self.output_reads.load(Ordering::Relaxed),
            malformed_mod_reads: self.malformed_mod_reads.load(Ordering::Relaxed),
            undo_tags_dropped_reads: self.undo_tags_dropped_reads.load(Ordering::Relaxed),
            barcode_tag_malformed_reads: self.barcode_tag_malformed_reads.load(Ordering::Relaxed),
            barcode_tag_unverified_reads: self.barcode_tag_unverified_reads.load(Ordering::Relaxed),
            input_bases: self.input_bases.load(Ordering::Relaxed),
            output_bases: self.output_bases.load(Ordering::Relaxed),
            malformed_tag_reads: self.malformed_tag_reads.load(Ordering::Relaxed),
            reads_with_output,
            reads_trimmed_to_nothing,
            reads_all_filtered,
            reads_tag_filtered,
            segments_dropped_short,
            segments_dropped_long,
            segments_dropped_low_qual,
            segments_dropped_high_qual,
            segments_dropped_gc,
        }
    }
}

/// Returns a trace-level span naming the read, so the segment decisions logged
/// while processing it are attributable without repeating the name on every
/// event.
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

/// Filters produced segments and updates segment- and read-level counters for
/// all workflows. `seq` and `qual` contain the complete input read and
/// `produced` contains the final ranges from adapter and quality processing
/// in original-coordinate order. Segment numbers index this flattened list
/// before filtering; they do not restart at adapter boundaries. For each survivor,
/// `render` receives `(idx, total, start, end)`; for each rejected segment and
/// for a read that produced none, `reject` receives the `Rejection`. A render
/// error stops processing before the read-level outcome counter is updated.
pub(crate) fn process_read_segments<Rn, Rj>(
    produced: &[(usize, usize)],
    seq: &[u8],
    qual: &[u8],
    filter_cfg: &FilterConfig,
    counters: &Counters,
    mut render: Rn,
    mut reject: Rj,
) -> anyhow::Result<()>
where
    Rn: FnMut(usize, usize, usize, usize) -> anyhow::Result<()>,
    Rj: FnMut(Rejection) -> anyhow::Result<()>,
{
    let total = produced.len();
    let mut survived = 0usize;
    for (idx, &(s, e)) in produced.iter().enumerate() {
        if let Some(reason) = crate::filter::check(&seq[s..e], &qual[s..e], filter_cfg) {
            // The per-segment verdict attributes a missing read to a specific
            // segment and reason, which no run-level counter can do.
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
            reject(Rejection::Segment {
                idx,
                total,
                start: s,
                end: e,
                reason,
            })?;
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
        reject(Rejection::Whole)?;
    } else if survived == 0 {
        tracing::trace!(produced = total, "Every segment filtered");
        counters.reads_all_filtered.fetch_add(1, Ordering::Relaxed);
    } else {
        counters.reads_with_output.fetch_add(1, Ordering::Relaxed);
    }
    Ok(())
}

/// What `process_read_segments` rejected: one trimmed segment, or the whole
/// read because trimming produced no segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Rejection {
    /// A produced segment that failed a post-trim filter.
    Segment {
        /// 0-based segment index.
        idx: usize,
        /// Number of segments produced from the read.
        total: usize,
        /// First base of the segment, inclusive.
        start: usize,
        /// End of the segment, exclusive.
        end: usize,
        /// The failed filter.
        reason: DropReason,
    },
    /// The read produced no segment.
    Whole,
}

/// End-of-run counters, snapshotted from `Counters`.
#[derive(Debug, Default, Clone, Copy)]
pub struct Stats {
    /// Input reads consumed.
    pub input_reads: u64,
    /// Output segments written.
    pub output_reads: u64,
    /// Sum of SEQ lengths (bases) across every input read.
    pub input_bases: u64,
    /// Sum of surviving segment lengths (bases) written to output.
    pub output_bases: u64,
    /// Reads carrying a malformed per-base tag, left untouched; see
    /// `Counters::malformed_tag_reads`. Surfaced as a run-level advisory; not
    /// an error.
    pub malformed_tag_reads: u64,
    /// Input reads whose modification block was malformed and removed; see
    /// `Counters::malformed_mod_reads`.
    pub malformed_mod_reads: u64,
    /// Trimmed reads whose `ds`/`ls` undo blobs were removed.
    pub undo_tags_dropped_reads: u64,
    /// Reads whose `bi` barcode tag was unusable; see
    /// `Counters::barcode_tag_malformed_reads`.
    pub barcode_tag_malformed_reads: u64,
    /// Reads with an unverified barcode span; see
    /// `Counters::barcode_tag_unverified_reads`.
    pub barcode_tag_unverified_reads: u64,
    /// Read-level: input reads with at least one written segment.
    pub reads_with_output: u64,
    /// Read-level: input reads that produced zero segments at all (empty
    /// read, fully consumed by adapter trimming, or an over-crop).
    /// `trim::apply` returned no intervals, so the per-segment filter loop
    /// never ran.
    pub reads_trimmed_to_nothing: u64,
    /// Read-level: input reads that produced at least one segment, but every
    /// one of them was rejected by post-trim `filter::check`.
    pub reads_all_filtered: u64,
    /// Read-level: input reads rejected by `--tag-filter` before trimming.
    pub reads_tag_filtered: u64,
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
        let by_weight: Vec<Vec<usize>> = Batches::new(
            vec![800_000usize; 5].into_iter(),
            |n: &usize| *n,
            FASTQ_BATCH,
        )
        .collect();
        assert_eq!(by_weight.iter().map(Vec::len).collect::<Vec<_>>(), [3, 2]);

        let bam: Vec<Vec<usize>> =
            Batches::new(vec![1usize; 65].into_iter(), |n: &usize| *n, BAM_BATCH).collect();
        assert_eq!(
            bam.iter().map(Vec::len).collect::<Vec<_>>(),
            [16, 16, 16, 16, 1]
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

        // (b) a read whose 2 produced segments are both filtered `TooShort`:
        // `reads_all_filtered` (segments were produced, but none survived).
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters.record_segment_drop(DropReason::TooShort);
        counters.record_segment_drop(DropReason::TooShort);
        counters.reads_all_filtered.fetch_add(1, Ordering::Relaxed);

        // (c) an empty input read: `trim::apply` produces no segments, so
        // `reads_trimmed_to_nothing` is bumped and no segment-level drop is
        // recorded.
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .reads_trimmed_to_nothing
            .fetch_add(1, Ordering::Relaxed);

        let stats = counters.snapshot();

        assert_eq!(stats.input_reads, 3);
        assert_eq!(stats.output_reads, 2);
        assert_eq!(
            stats.reads_all_filtered, 1,
            "Read b produced segments, but every one was filtered"
        );
        assert_eq!(
            stats.reads_trimmed_to_nothing, 1,
            "Read c produced no segments"
        );
        assert_eq!(counters.reads_with_output.load(Ordering::Relaxed), 1);
        // `reads_with_output + reads_trimmed_to_nothing + reads_all_filtered ==
        // input_reads` holds (also asserted by `snapshot`).
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

    /// Covers all read-level outcomes and the corresponding render arguments.
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

        // Trimmed to nothing: no produced intervals, so `render` is never called,
        // `reads_trimmed_to_nothing` is bumped and no segment-level drop is
        // recorded.
        {
            let counters = Counters::default();
            let mut calls: Vec<(usize, usize, usize, usize)> = Vec::new();
            process_read_segments(
                &[],
                b"",
                b"",
                &filter_cfg,
                &counters,
                |idx, total, s, e| {
                    calls.push((idx, total, s, e));
                    Ok(())
                },
                |_| Ok(()),
            )
            .unwrap();
            assert!(calls.is_empty());
            assert_eq!(counters.reads_trimmed_to_nothing.load(Ordering::Relaxed), 1);
            assert_eq!(counters.reads_all_filtered.load(Ordering::Relaxed), 0);
            assert_eq!(counters.reads_with_output.load(Ordering::Relaxed), 0);
            assert_eq!(counters.segments_dropped_short.load(Ordering::Relaxed), 0);
        }

        // All filtered: one produced segment, too short to pass, so `render` is
        // not called, `reads_all_filtered` is bumped and one segment drop is
        // recorded.
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
                |_| Ok(()),
            )
            .unwrap();
            assert!(calls.is_empty());
            assert_eq!(counters.reads_trimmed_to_nothing.load(Ordering::Relaxed), 0);
            assert_eq!(counters.reads_all_filtered.load(Ordering::Relaxed), 1);
            assert_eq!(counters.reads_with_output.load(Ordering::Relaxed), 0);
            assert_eq!(counters.segments_dropped_short.load(Ordering::Relaxed), 1);
        }

        // With output: two produced segments, both long enough, so `render` is
        // called once per survivor with `(idx, total, s, e)`.
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
                |_| Ok(()),
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
            |n, _cfg, out: &mut Vec<usize>| {
                if n % 2 == 1 {
                    std::thread::sleep(std::time::Duration::from_micros(200));
                }
                counters.reads_with_output.fetch_add(1, Ordering::Relaxed);
                out.push(n);
                Ok(())
            },
            Ok,
            |sink, batch: &Vec<usize>| {
                sink.extend_from_slice(batch);
                Ok(())
            },
            &counters,
        )
        .unwrap();
        sink
    }

    #[test]
    fn ordered_driver_bounds_read_ahead_behind_a_slow_batch() {
        use std::sync::atomic::AtomicUsize;
        use std::time::Duration;
        let cfg = driver_cfg(4, true);
        let consumed = AtomicUsize::new(0);
        let (ready_tx, ready_rx) = std::sync::mpsc::channel();
        let (release_tx, release_rx) = std::sync::mpsc::channel();
        let release_rx = Mutex::new(release_rx);
        let mut output = Vec::new();
        let counters = Counters::default();
        std::thread::scope(|scope| {
            let consumed_ref = &consumed;
            let observed = scope.spawn(move || {
                for _ in 0..7 {
                    ready_rx.recv_timeout(Duration::from_secs(5)).unwrap();
                }
                std::thread::sleep(Duration::from_millis(20));
                let read_ahead = consumed_ref.load(Ordering::Relaxed);
                release_tx.send(()).unwrap();
                read_ahead
            });
            run_parallel(
                (0..1000usize).map(|n| {
                    consumed.fetch_add(1, Ordering::Relaxed);
                    Ok(n)
                }),
                BatchPolicy {
                    target_weight: 1,
                    max_items: 1,
                    queue_per_worker: 1,
                },
                |_| 1,
                &cfg,
                &mut output,
                |n, _, out: &mut Vec<usize>| {
                    if n == 0 {
                        release_rx
                            .lock()
                            .unwrap()
                            .recv_timeout(Duration::from_secs(5))
                            .unwrap();
                    } else if n < 8 {
                        ready_tx.send(()).unwrap();
                    }
                    counters.reads_with_output.fetch_add(1, Ordering::Relaxed);
                    out.push(n);
                    Ok(())
                },
                Ok,
                |out, batch: &Vec<usize>| {
                    out.extend_from_slice(batch);
                    Ok(())
                },
                &counters,
            )
            .unwrap();
            assert_eq!(observed.join().unwrap(), 8);
        });
        assert_eq!(output, (0..1000).collect::<Vec<_>>());
    }

    #[test]
    fn ordered_driver_propagates_write_errors() {
        let cfg = driver_cfg(4, true);
        let error = run_parallel(
            (0..1000usize).map(Ok),
            BatchPolicy {
                target_weight: 1,
                max_items: 1,
                queue_per_worker: 1,
            },
            |_| 1,
            &cfg,
            &mut (),
            |n, _, out: &mut Vec<usize>| {
                out.push(n);
                Ok(())
            },
            Ok,
            |_, _: &Vec<usize>| Err(std::io::Error::other("sink failed")),
            &Counters::default(),
        )
        .unwrap_err();
        assert_eq!(error.to_string(), "sink failed");
    }

    /// A batch that fails to pack is never written, so the ordered writer
    /// cannot advance past it; the run still ends with the error rather than
    /// waiting on the reorder window.
    #[test]
    fn ordered_driver_ends_on_a_pack_error() {
        let cfg = driver_cfg(4, true);
        let error = run_parallel(
            (0..1000usize).map(Ok),
            BatchPolicy {
                target_weight: 1,
                max_items: 1,
                queue_per_worker: 1,
            },
            |_| 1,
            &cfg,
            &mut Vec::new(),
            |n, _, out: &mut Vec<usize>| {
                out.push(n);
                Ok(())
            },
            |batch: Vec<usize>| {
                if batch == [3] {
                    return Err(std::io::Error::other("pack failed"));
                }
                Ok(batch)
            },
            |out: &mut Vec<usize>, batch: &Vec<usize>| {
                out.extend_from_slice(batch);
                Ok(())
            },
            &Counters::default(),
        )
        .unwrap_err();
        assert_eq!(error.to_string(), "pack failed");
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
                assert!(!self.yielded_error, "Source polled after it returned Err");
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
            |n, _cfg, out: &mut Vec<usize>| {
                out.push(n);
                Ok(())
            },
            Ok,
            |sink, batch: &Vec<usize>| {
                sink.extend_from_slice(batch);
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
            |n, _cfg, out: &mut Vec<usize>| {
                rendered.fetch_add(1, Ordering::Relaxed);
                if n == 10 {
                    anyhow::bail!("record 10 is malformed");
                }
                out.push(n);
                Ok(())
            },
            Ok,
            |sink, batch: &Vec<usize>| {
                sink.extend_from_slice(batch);
                Ok(())
            },
            &counters,
        );
        assert!(res.is_err());
        assert!(
            rendered.load(Ordering::Relaxed) < 1_000,
            "Rendering continued long after the first error: {} records",
            rendered.load(Ordering::Relaxed)
        );
    }
}
