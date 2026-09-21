//! FASTQ workflow over paraseq: each worker fills its own record set from the
//! shared byte stream, then renders, compresses and hands off the set, so
//! record copying leaves the serial reader.

use std::io::Read;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{SyncSender, sync_channel};

use paraseq::prelude::*;

use super::fastq::render_record;
use super::{BatchSink, Counters, FirstError, Stats};
use crate::config::Config;
use crate::record::ReadRecord;

/// Records per set handed to one worker.
const RECORDS_PER_SET: usize = 128;

/// Sets queued to the writer per worker.
const QUEUE_PER_WORKER: usize = 2;

/// Environment variable selecting the FASTQ parser: `paraseq` for this
/// path, anything else for the seq_io iterator.
pub(crate) const PARSER_VAR: &str = "WHITTLE_FASTQ_PARSER";

/// Whether the environment selects this path.
pub(crate) fn selected() -> bool {
    std::env::var_os(PARSER_VAR).is_some_and(|v| v == "paraseq")
}

/// One worker's state: its render buffer, its counts since the last hand-off
/// and the writer channel.
#[derive(Clone)]
struct Processor {
    cfg: Arc<Config>,
    counters: Arc<Counters>,
    tx: SyncSender<Vec<u8>>,
    level: Option<u8>,
    ordered: bool,
    aborted: Arc<AtomicBool>,
    buf: Vec<u8>,
    reads: u64,
    bases: u64,
}

/// Highest raw Phred score a Phred+33 quality byte encodes.
const MAX_PHRED33: u8 = 126 - 33;

/// Copies a parsed record into an owned record with raw Phred qualities.
fn to_read_record<Rf: Record>(record: &Rf) -> anyhow::Result<ReadRecord> {
    let raw = record.qual().unwrap_or_default();
    let mut out_of_range = false;
    let qual: Vec<u8> = raw
        .iter()
        .map(|&b| {
            let q = b.wrapping_sub(33);
            out_of_range |= q > MAX_PHRED33;
            q
        })
        .collect();
    if out_of_range {
        return Err(crate::io::fastq::invalid_quality(record.id(), raw));
    }
    Ok(ReadRecord {
        name: record.id().to_vec(),
        seq: record.seq_raw().to_vec(),
        qual,
    })
}

impl<Rf: Record> ParallelProcessor<Rf> for Processor {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        if self.aborted.load(Ordering::Relaxed) {
            return Err(anyhow::anyhow!("run aborted").into());
        }
        let rec = to_read_record(&record)?;
        self.reads += 1;
        self.bases += rec.seq.len() as u64;
        render_record(rec, &self.cfg, &self.counters, &mut self.buf)?;
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
        let bytes = std::mem::take(&mut self.buf);
        let packed = match self.level {
            Some(level) => crate::io::fastq::encode_blocks(level, &bytes)?,
            None => bytes,
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
/// `cfg.ordered` the sets are handed off in input order.
pub(crate) fn run_fastq_paraseq<W: BatchSink>(
    input: Box<dyn Read + Send>,
    writer: &mut W,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    let workers = super::render_pool_size(cfg);
    let queue = (workers * QUEUE_PER_WORKER).max(2);
    let (tx, rx) = sync_channel::<Vec<u8>>(queue);
    let aborted = Arc::new(AtomicBool::new(false));
    let write_err: FirstError<std::io::Error> = FirstError::new();
    let reader = paraseq::fastq::Reader::with_batch_size(input, RECORDS_PER_SET)?;
    let mut processor = Processor {
        cfg: Arc::new(cfg.clone()),
        counters: Arc::clone(counters),
        tx,
        level: writer.block_level(),
        ordered: cfg.ordered,
        aborted: Arc::clone(&aborted),
        buf: Vec::new(),
        reads: 0,
        bases: 0,
    };

    let result = std::thread::scope(|s| {
        let write_err = &write_err;
        let aborted = &aborted;
        s.spawn(move || {
            let mut errored = false;
            for batch in rx.iter() {
                if errored {
                    continue;
                }
                if let Err(e) = writer.write_all(&batch) {
                    write_err.record(e, aborted);
                    errored = true;
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
