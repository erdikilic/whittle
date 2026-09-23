//! uBAM workflows: record reconstruction (sequence, quality, MM/ML/MN, per-base and signal tags) and the sequential, parallel and raw full-window drivers for BAM and FASTQ output.

use std::borrow::Cow;
use std::io;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};

use super::reject::{self, Reason, RejectItem};
use super::{
    BAM_BATCH, BatchSink, Counters, Rejection, Stats, process_read_segments, run_bytes_parallel,
    run_parallel,
};
use crate::config::{Config, FastqTags, TagRemoval};
use crate::io::fastq::{push_aux_field, push_mods_aux, push_record_body};
use crate::mods::reconstruct::IndexedMods;
use crate::{mods, trim};

mod barcode;
mod modblock;
mod naming;
mod raw;
mod render;
mod signal;
mod tags;
mod to_fastq;
pub(crate) use barcode::*;
pub(crate) use modblock::*;
pub(crate) use naming::*;
use raw::*;
pub(crate) use render::*;
pub(crate) use signal::*;
pub(crate) use tags::*;
pub(crate) use to_fastq::*;

/// ONT signal-mapping tags: the `mv` move table plus the `ts`/`ns` sample counts
/// and the `sp`/`pi` split linkage. On a trimmed read these are either rewritten
/// (`--update-moves`) or dropped (default), never left stale. Handled by
/// `signal_tag_updates`, not the per-base pass.
pub(crate) use crate::config::SIGNAL_TAGS;

/// Renders a raw record rejected by the tag filter for BAM output.
pub(crate) fn tag_filtered_bam(record: &bam::Record) -> anyhow::Result<RejectItem> {
    let mut rec = decode_raw_record(record)?;
    reject::tag_record(&mut rec, Reason::TagFilter);
    Ok(RejectItem::Bam(rec))
}

/// Renders a raw record rejected by the tag filter for FASTQ output: the whole
/// read with its selected tags and the reason tag.
pub(crate) fn tag_filtered_bam_fastq(
    record: &bam::Record,
    cfg: &Config,
) -> anyhow::Result<RejectItem> {
    let rec = decode_raw_record(record)?;
    let seq_len = rec.sequence().len();
    let mut out = Vec::new();
    render_fastq_window(
        &mut out,
        &rec,
        &[],
        Window {
            start: 0,
            end: seq_len,
            idx: 0,
            total: 1,
        },
        inspect_mod_block(&rec, seq_len),
        None,
        platform(&rec),
        &cfg.fastq_tags,
        &cfg.remove_tags,
        Some(Reason::TagFilter),
    );
    Ok(RejectItem::Fastq(out))
}

/// Runs the single-threaded uBAM workflow: refuses aligned reads, trims, filters
/// each produced segment and writes the reconstructed survivors.
fn run_bam_seq(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>>,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    for rec in records {
        let rec = decode_raw_record(&rec?)?;
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .input_bases
            .fetch_add(rec.sequence().as_ref().len() as u64, Ordering::Relaxed);
        render_bam_read(header, &rec, cfg, counters, |out| {
            sink.write_record(header, out.as_ref().unwrap_or(&rec))
        })?;
    }
    Ok(counters.snapshot())
}

/// Runs `workflow::run_parallel` for BAM input: decodes each raw record on the
/// pool and hands the decoded record to `render`, which appends output items
/// to the batch buffer. The per-segment filter and counters are updated inside
/// `render` by `process_read_segments`.
fn run_bam_parallel<T, P, S, Render, Pack, WriteOne>(
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    cfg: &Config,
    sink: &mut S,
    render: Render,
    pack: Pack,
    write_one: WriteOne,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats>
where
    T: Send,
    P: Send,
    S: Send,
    Render: Fn(&bam::Record, &RecordBuf, &Config, &mut Vec<T>) -> anyhow::Result<()> + Sync,
    Pack: Fn(Vec<T>) -> std::io::Result<P> + Sync,
    WriteOne: Fn(&mut S, &P) -> std::io::Result<()> + Send,
{
    run_parallel(
        records,
        BAM_BATCH,
        |record: &bam::Record| record.sequence().len(),
        cfg,
        sink,
        |rec, cfg, out| render(&rec, &decode_raw_record(&rec)?, cfg, out),
        pack,
        write_one,
        counters,
    )
}

/// Compresses a batch of output records into BGZF blocks at the sink's level.
fn pack_bam_blocks(
    header: &sam::Header,
    level: u8,
    records: Vec<BamOutputRecord>,
) -> std::io::Result<Vec<u8>> {
    use noodles_sam::alignment::io::Write as _;
    let mut w = bam::io::Writer::from(Vec::new());
    for rec in &records {
        match rec {
            BamOutputRecord::Raw(record) => w.write_record(header, record)?,
            BamOutputRecord::Decoded(record) => w.write_alignment_record(header, record)?,
        }
    }
    let mut blocks = Vec::new();
    crate::io::bgzf::encode(level, &w.into_inner(), &mut blocks)?;
    Ok(blocks)
}

/// Runs the uBAM workflow on raw records from a production reader. Full-window
/// runs filter and write unchanged records without building an owned
/// `RecordBuf`; any configuration that can alter sequence or tags is routed to
/// `run_bam`, tag removal included, since a record that would otherwise pass
/// through untouched still has to be rebuilt without the removed tags.
pub fn run_raw_bam(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    let full_window = cfg.trim.head == 0
        && cfg.trim.tail == 0
        && cfg.trim.quality.is_none()
        && cfg.adapters.is_none()
        && cfg.remove_tags.is_empty();
    if !full_window {
        return run_bam(header, records, sink, cfg, counters);
    }
    if cfg.threads <= 1 {
        run_raw_bam_full_window_seq(header, records, sink, cfg, counters)
    } else {
        run_raw_bam_full_window_parallel(header, records, sink, cfg, counters)
    }
}

/// Runs the uBAM workflow: decodes, refuses aligned reads, trims, filters and
/// reconstructs. Sequential for `cfg.threads <= 1`; otherwise renders on a
/// rayon pool and drains the `RecordBuf`s through `run_bam_parallel`'s bounded
/// channel to the writer, in input order under `cfg.ordered` and in completion
/// order otherwise.
pub(crate) fn run_bam(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    if cfg.threads <= 1 {
        return run_bam_seq(header, records, sink, cfg, counters);
    }
    let level = sink
        .block_level()
        .expect("A parallel run writes through a block sink");
    run_bam_parallel(
        records,
        cfg,
        sink,
        // Render: the survivors of one record. An untouched record is
        // written from its raw input without re-encoding.
        |raw, rec, cfg, items| {
            render_bam_read(header, rec, cfg, counters, |out| {
                items.push(match out {
                    Some(out) => BamOutputRecord::Decoded(out),
                    None => BamOutputRecord::Raw(raw.clone()),
                });
                Ok(())
            })
        },
        // Pack: encode and compress the batch on the pool.
        |records| pack_bam_blocks(header, level, records),
        // Write: the compressed blocks, on the writer thread.
        |sink, blocks: &Vec<u8>| sink.write_blocks(blocks),
        counters,
    )
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod barcode_tests;
