//! FASTQ workflow over paraseq: each worker fills its own record set from the
//! shared byte stream, then renders, compresses and hands off the set, so
//! record copying leaves the serial reader.

use std::io::Read;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{SyncSender, sync_channel};
use std::sync::{Arc, Mutex};

use paraseq::prelude::*;

use super::fastq::render_record_ref;
use super::{BatchSink, Counters, FirstError, Parts, Stats};
use crate::config::Config;

/// Most records per set handed to one worker.
const RECORDS_PER_SET: usize = 128;

/// Bases per set the record count is sized for, from the first record's
/// length, so a set of long reads stays near the seq_io path's batch weight.
const BASES_PER_SET: usize = super::FASTQ_BATCH.target_weight;

/// Sets queued to the writer per worker.
const QUEUE_PER_WORKER: usize = 2;

/// Environment variable selecting the FASTQ parser: `paraseq` for this
/// path, anything else for the seq_io iterator.
pub(crate) const PARSER_VAR: &str = "WHITTLE_FASTQ_PARSER";

/// Whether the environment selects this path.
pub(crate) fn selected() -> bool {
    std::env::var_os(PARSER_VAR).is_some_and(|v| v == "paraseq")
}

/// Drained output buffers the writer returns for workers to render into.
type BufferPool = Arc<Mutex<Vec<Vec<u8>>>>;

/// One worker's state: its render parts (key 0 only, since a paraseq run has
/// no `--targets`), its raw Phred scratch buffer, its counts since the last
/// hand-off, the writer channel and the shared buffer pool.
#[derive(Clone)]
struct Processor {
    cfg: Arc<Config>,
    counters: Arc<Counters>,
    tx: SyncSender<Vec<u8>>,
    reuse: BufferPool,
    level: Option<u8>,
    ordered: bool,
    aborted: Arc<AtomicBool>,
    parts: Parts<u8>,
    qual: Vec<u8>,
    reads: u64,
    bases: u64,
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        if self.aborted.load(Ordering::Relaxed) {
            return Err(anyhow::anyhow!("run aborted").into());
        }
        let head = record.id();
        let seq = record.seq_raw();
        let id = crate::io::fastq::read_id(head);
        crate::io::fastq::phred_from_ascii_into(
            id,
            record.qual().unwrap_or_default(),
            &mut self.qual,
        )?;
        self.reads += 1;
        self.bases += seq.len() as u64;
        render_record_ref(
            head,
            seq,
            &self.qual,
            &self.cfg,
            &self.counters,
            &mut self.parts,
        )?;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.counters
            .input_reads
            .fetch_add(self.reads, Ordering::Relaxed);
        self.counters
            .input_bases
            .fetch_add(self.bases, Ordering::Relaxed);
        self.reads = 0;
        self.bases = 0;
        debug_assert!(
            self.parts.only_first(),
            "a paraseq run renders to key 0 only"
        );
        let packed = match self.level {
            Some(level) => {
                // Clearing in place keeps the buffer's capacity for the next set.
                let buf = self.parts.part(0);
                let out = crate::io::fastq::encode_blocks(level, buf)?;
                buf.clear();
                out
            },
            None => {
                let fresh = self
                    .reuse
                    .lock()
                    .ok()
                    .and_then(|mut pool| pool.pop())
                    .unwrap_or_default();
                std::mem::replace(self.parts.part(0), fresh)
            },
        };
        if self.tx.send(packed).is_err() {
            self.aborted.store(true, Ordering::Relaxed);
            return Err(anyhow::anyhow!("the output writer stopped").into());
        }
        Ok(())
    }

    fn requires_ordering(&self) -> bool {
        self.ordered
    }
}

/// Runs the FASTQ workflow over `input` on `render_pool_size(cfg)` paraseq
/// workers, writing rendered sets through one writer thread. Under
/// `cfg.ordered` the sets are handed off in input order. `level` is the
/// writer's BGZF level when it takes compressed blocks, `None` when it takes
/// the rendered bytes.
/// Parse errors, including a record cut short at the end of the input, name
/// `cfg.io.input` when it is a file.
pub(crate) fn run_fastq_paraseq<W: BatchSink>(
    input: Box<dyn Read + Send>,
    writer: &mut W,
    level: Option<u8>,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    let workers = super::render_pool_size(cfg);
    let queue = (workers * QUEUE_PER_WORKER).max(2);
    let (tx, rx) = sync_channel::<Vec<u8>>(queue);
    // Every buffer in flight is either queued or held by a worker.
    let pool_cap = queue + workers;
    let reuse: BufferPool = Arc::new(Mutex::new(Vec::with_capacity(pool_cap)));
    let aborted = Arc::new(AtomicBool::new(false));
    let write_err: FirstError<std::io::Error> = FirstError::new();
    let mut reader = paraseq::fastq::Reader::with_batch_size(input, RECORDS_PER_SET)?;
    if let Some(path) = cfg.io.input.as_deref() {
        reader.set_source_name(path.display().to_string());
    }
    reader.update_batch_size_in_bp(BASES_PER_SET)?;
    reader.set_batch_size(reader.batch_size().min(RECORDS_PER_SET))?;
    let mut processor = Processor {
        cfg: Arc::new(cfg.clone()),
        counters: Arc::clone(counters),
        tx,
        reuse: Arc::clone(&reuse),
        level,
        ordered: cfg.ordered,
        aborted: Arc::clone(&aborted),
        parts: Parts::default(),
        qual: Vec::new(),
        reads: 0,
        bases: 0,
    };

    let result = std::thread::scope(|s| {
        let write_err = &write_err;
        let aborted = &aborted;
        let reuse = &reuse;
        s.spawn(move || {
            let mut errored = false;
            for mut batch in &rx {
                if errored {
                    continue;
                }
                if let Err(e) = writer.write_all(&batch) {
                    write_err.record(e, aborted);
                    errored = true;
                }
                if level.is_none() {
                    batch.clear();
                    if let Ok(mut pool) = reuse.lock()
                        && pool.len() < pool_cap
                    {
                        pool.push(batch);
                    }
                }
            }
        });
        let result = reader.process_parallel(&mut processor, workers);
        drop(processor);
        result
    });

    if let Some(e) = write_err.take() {
        return Err(e.into());
    }
    result?;
    Ok(counters.snapshot())
}
