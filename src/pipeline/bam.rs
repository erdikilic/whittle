use std::io::Write;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};
use rayon::prelude::*;

use super::{Counters, Stats};
use crate::config::{Config, FastqTags};
use crate::io::fastq::{format_aux_field, format_mods_aux, write_segment, write_segment_tagged};
use crate::{filter, mods, trim};

/// PacBio per-base kinetics tags: one value per SEQ base (`B` arrays), so they
/// must be sliced in lockstep with the sequence when a read is trimmed. `ip`/`pw`
/// are single-strand IPD/pulse-width; `fi`/`fp`/`ri`/`rp` are the CCS forward/
/// reverse codec-V1 kinetics. Any *other* `B` array whose length equals the read
/// length is also treated as per-base (structural rule) so new/custom tags are
/// handled without a code change.
const KNOWN_PERBASE_TAGS: [[u8; 2]; 6] = [*b"ip", *b"pw", *b"fi", *b"fp", *b"ri", *b"rp"];

/// ONT signal-mapping tags: the `mv` move table plus the `ts`/`ns` sample counts
/// and the `sp`/`pi` split linkage. On a trimmed read these are either rewritten
/// (`--update-moves`) or dropped (default) — never left stale. Handled by
/// `signal_tag_updates`, not the per-base pass.
const SIGNAL_TAGS: [[u8; 2]; 5] = [*b"mv", *b"ts", *b"ns", *b"sp", *b"pi"];

/// Poly-A tail tags handled together with the move table: `pa` (signal boundaries,
/// stored in original-signal coordinates) and `pt` (tail length in bases). Under
/// `--update-moves` they are kept/shifted when the poly-A tail survives the trim
/// and dropped when it's cut; without it (or a malformed move table) they are
/// dropped, since we can't relate signal to sequence.
const POLYA_TAGS: [[u8; 2]; 2] = [*b"pa", *b"pt"];

/// `bi` (barcode info) embeds front/rear SEQUENCE positions that shift under a
/// crop and can't be reconstructed from the BAM, so it is dropped on any trimmed
/// read. The barcode call itself (`BC`/`bv`) is a per-read label and rides
/// through unchanged.
const DROP_ON_TRIM_TAGS: [[u8; 2]; 1] = [*b"bi"];

/// Tags dropped only when a read is SPLIT (not on a plain crop): `st` (read start
/// time) and `du` (duration) describe the whole parent read, but a split subread
/// starts later in the signal and spans less of it. Dorado recomputes both from
/// the sample rate, which isn't carried in the BAM, so we drop them rather than
/// ship a stale timestamp/duration. A head/tail crop keeps the same read identity,
/// so they stay valid there.
const DROP_ON_SPLIT_TAGS: [[u8; 2]; 2] = [*b"st", *b"du"];

fn array_len(a: &Array) -> usize {
    match a {
        Array::Int8(v) => v.len(),
        Array::UInt8(v) => v.len(),
        Array::Int16(v) => v.len(),
        Array::UInt16(v) => v.len(),
        Array::Int32(v) => v.len(),
        Array::UInt32(v) => v.len(),
        Array::Float(v) => v.len(),
    }
}

/// Slice every subtype of a `B` array to `[start, end)` (the element index is the
/// base index for a per-base tag). Subtype-agnostic, so `B:C` (codec-V1) and `B:S`
/// (raw frames) kinetics both work.
fn slice_array(a: &Array, start: usize, end: usize) -> Array {
    match a {
        Array::Int8(v) => Array::Int8(v[start..end].to_vec()),
        Array::UInt8(v) => Array::UInt8(v[start..end].to_vec()),
        Array::Int16(v) => Array::Int16(v[start..end].to_vec()),
        Array::UInt16(v) => Array::UInt16(v[start..end].to_vec()),
        Array::Int32(v) => Array::Int32(v[start..end].to_vec()),
        Array::UInt32(v) => Array::UInt32(v[start..end].to_vec()),
        Array::Float(v) => Array::Float(v[start..end].to_vec()),
    }
}

/// A per-base `B` array (any array whose length equals the read length — this
/// covers the known kinetics tags and any custom per-base tag) sliced to the
/// window `[start, end)`; `None` to leave the tag unchanged. Callers must already
/// have excluded MM/ML/MN and the signal tags. A known kinetics tag whose length
/// does NOT match is left unchanged and surfaced via `has_malformed_perbase_tag`.
fn perbase_slice(value: &Value, orig_len: usize, start: usize, end: usize) -> Option<Value> {
    match value {
        Value::Array(arr) if array_len(arr) == orig_len => {
            Some(Value::Array(slice_array(arr, start, end)))
        },
        _ => None,
    }
}

/// Parse an `mv` move table value into `(stride, moves)`. `None` unless it is a
/// `B:c` (Int8) array with a positive stride. `moves` excludes the stride prefix;
/// each entry corresponds to `stride` signal samples (1 = a base emitted here, so
/// the count of 1s equals the sequence length).
fn parse_move_table(value: &Value) -> Option<(i8, &[i8])> {
    match value {
        Value::Array(Array::Int8(a)) => {
            let (stride, moves) = a.split_first()?;
            if *stride > 0 {
                Some((*stride, moves))
            } else {
                None
            }
        },
        _ => None,
    }
}

/// Read an integer aux tag as `i64`, regardless of stored width.
fn signal_int(src: &RecordBuf, tag: &[u8; 2]) -> Option<i64> {
    match src.data().get(&Tag::new(tag[0], tag[1])) {
        Some(Value::Int8(n)) => Some(i64::from(*n)),
        Some(Value::UInt8(n)) => Some(i64::from(*n)),
        Some(Value::Int16(n)) => Some(i64::from(*n)),
        Some(Value::UInt16(n)) => Some(i64::from(*n)),
        Some(Value::Int32(n)) => Some(i64::from(*n)),
        Some(Value::UInt32(n)) => Some(i64::from(*n)),
        _ => None,
    }
}

/// The parent read id for a subread: the source's own `pi` if it already has one
/// (so `pi` always names the ultimate ancestor, matching dorado), else the source
/// read name.
fn parent_read_id(src: &RecordBuf) -> Vec<u8> {
    match src.data().get(&Tag::new(b'p', b'i')) {
        Some(Value::String(s)) => s.to_vec(),
        _ => src.name().map(|n| n.to_vec()).unwrap_or_default(),
    }
}

/// Trim-aware handling of the poly-A tags (`pa` signal boundaries, `pt` tail
/// length). `pa` holds absolute original-signal positions (`>= 0`; `-1`/`-2` are
/// dorado's not-found/not-enabled sentinels, left as-is). `[kept_start, kept_end)`
/// is the original-signal window the kept bases span. If every real `pa` position
/// falls inside that window the tail survived, so: on a split, shift `pa` into the
/// subread's own signal frame (its `ts` is 0) and keep `pt`; on a crop, keep both
/// unchanged (identity + POD5 signal are unchanged, so the absolute positions stay
/// valid). If any real position falls outside — the tail was (partly) trimmed — or
/// there's no poly-A array, drop `pa`/`pt`.
fn polya_updates(
    src: &RecordBuf,
    kept_start: i64,
    kept_end: i64,
    is_split: bool,
) -> Vec<(Tag, Option<Value>)> {
    let pa_tag = Tag::new(b'p', b'a');
    let pt_tag = Tag::new(b'p', b't');
    let drop_both = || vec![(pa_tag, None), (pt_tag, None)];

    let pa = match src.data().get(&pa_tag) {
        Some(Value::Array(Array::Int32(v))) => v,
        _ => return drop_both(),
    };
    // pa = [anchor, range0.start, range0.end, range1.start, range1.end]. dorado's
    // poly-A signal ranges are half-open `[start, end)`: the anchor and the range
    // starts are inclusive sample indices, so they must be `< kept_end`; the range
    // ENDS are exclusive and may equal `kept_end` (the window's own exclusive end).
    // Every real position must also be `>= kept_start`. Sentinels (`< 0`) skipped.
    let has_real = pa.iter().any(|&p| p >= 0);
    let survives = has_real
        && pa.iter().enumerate().all(|(i, &p)| {
            if p < 0 {
                return true; // sentinel (NOT_FOUND / NOT_ENABLED)
            }
            let p = i64::from(p);
            let within_upper = if i == 2 || i == 4 {
                p <= kept_end
            } else {
                p < kept_end
            };
            p >= kept_start && within_upper
        });
    if !survives {
        return drop_both();
    }
    if is_split {
        // Re-express into the subread's own frame (subread signal 0 == kept_start;
        // its ts is 0). Sentinels stay untouched. pt (base count) is unchanged.
        let shifted: Vec<i32> = pa
            .iter()
            .map(|&p| {
                if p >= 0 {
                    (i64::from(p) - kept_start) as i32
                } else {
                    p
                }
            })
            .collect();
        vec![(pa_tag, Some(Value::Array(Array::Int32(shifted))))]
    } else {
        // Crop: absolute original-signal positions are still valid; keep pa/pt.
        Vec::new()
    }
}

/// Trim-aware rewrite of the ONT signal tags for output window `[start, end)`.
/// Returns `(tag, Some(value))` to set or `(tag, None)` to remove; empty when the
/// read isn't trimmed. With `update_moves` off — or when the move table is
/// missing/malformed — all five signal tags are removed. With it on, the move
/// table is sliced by block range `moves[block_first .. block_second]`
/// (stride-aligned, following dorado `splitter::subread`) and:
///   - crop (`total == 1`, name kept): `mv` sliced, `ts += block_first*stride`,
///     `ns` unchanged (still matches the by-name POD5 signal length).
///   - split (`total > 1`, renamed): dorado subread encoding — `mv` sliced,
///     `ts = 0`, `ns = span*stride`, `sp = parent offset`, `pi = parent id`.
fn signal_tag_updates(
    src: &RecordBuf,
    seq_len: usize,
    start: usize,
    end: usize,
    total: usize,
    update_moves: bool,
) -> Vec<(Tag, Option<Value>)> {
    if start == 0 && end == seq_len {
        return Vec::new(); // untrimmed: leave everything
    }
    let drop_all = || -> Vec<(Tag, Option<Value>)> {
        SIGNAL_TAGS
            .iter()
            .chain(POLYA_TAGS.iter())
            .map(|t| (Tag::new(t[0], t[1]), None))
            .collect()
    };
    if !update_moves {
        return drop_all();
    }

    // A consistent move table (1-count == sequence length) is required to slice.
    let Some((stride, moves)) = src
        .data()
        .get(&Tag::new(b'm', b'v'))
        .and_then(parse_move_table)
    else {
        return drop_all();
    };
    let ones: Vec<usize> = moves
        .iter()
        .enumerate()
        .filter(|(_, m)| **m != 0)
        .map(|(i, _)| i)
        .collect();
    if ones.len() != seq_len {
        return drop_all(); // move table inconsistent with the sequence
    }

    let stride_n = stride as usize;
    let block_first = ones[start];
    let block_second = if end == seq_len {
        moves.len()
    } else {
        ones[end]
    };

    let mut new_mv = Vec::with_capacity(1 + block_second - block_first);
    new_mv.push(stride);
    new_mv.extend_from_slice(&moves[block_first..block_second]);
    let mut updates = vec![(
        Tag::new(b'm', b'v'),
        Some(Value::Array(Array::Int8(new_mv))),
    )];

    // Original-signal window the kept bases span: [ts0 + block_first*stride,
    // ts0 + block_second*stride). `ns = span + front trim` matches dorado's
    // `ns = raw_data_samples + num_trimmed_samples` (a tail crop shrinks ns, a
    // head-only crop leaves it unchanged, a split gets the subread span).
    let ts0 = signal_int(src, b"ts").unwrap_or(0);
    let kept_start = ts0 + (block_first * stride_n) as i64;
    let kept_end = ts0 + (block_second * stride_n) as i64;
    let span = ((block_second - block_first) * stride_n) as i64;

    if total > 1 {
        // Split -> dorado subread: renamed, front trim reset to 0, parent linkage.
        let sp = signal_int(src, b"sp").unwrap_or(0) + (block_first * stride_n) as i64;
        let pi = parent_read_id(src);
        updates.push((Tag::new(b't', b's'), Some(Value::Int32(0))));
        updates.push((Tag::new(b'n', b's'), Some(Value::Int32(span as i32))));
        updates.push((Tag::new(b's', b'p'), Some(Value::Int32(sp as i32))));
        updates.push((Tag::new(b'p', b'i'), Some(Value::String(pi.into()))));
        // Dorado marks split products with read_number -1.
        updates.push((Tag::new(b'r', b'n'), Some(Value::Int32(-1))));
    } else {
        // Head/tail crop in place: keep the read identity, advance the front trim.
        updates.push((Tag::new(b't', b's'), Some(Value::Int32(kept_start as i32))));
        updates.push((
            Tag::new(b'n', b's'),
            Some(Value::Int32((kept_start + span) as i32)),
        ));
    }
    updates.extend(polya_updates(src, kept_start, kept_end, total > 1));
    updates
}

/// True if the record carries a *known* per-base kinetics tag whose array length
/// disagrees with the sequence length — i.e. a malformed/unexpected per-base tag
/// that cannot be safely sliced. Used only to emit a run-level advisory.
pub fn has_malformed_perbase_tag(rec: &RecordBuf, seq_len: usize) -> bool {
    rec.data().iter().any(|(tag, value)| {
        let t = <[u8; 2]>::from(tag);
        KNOWN_PERBASE_TAGS.contains(&t)
            && matches!(value, Value::Array(a) if array_len(a) != seq_len)
    })
}

/// Build one output uBAM record for interval [start,end): slice SEQ/QUAL, rebuild
/// MM/ML/MN, slice per-base kinetics arrays, drop stale signal-space tags on trim,
/// suffix the name on splits. Remaining aux tags ride through unchanged.
pub fn reconstruct_record(
    src: &RecordBuf,
    start: usize,
    end: usize,
    total: usize,
    idx: usize,
    update_moves: bool,
) -> RecordBuf {
    let mut out = src.clone();

    // Slice sequence + quality.
    let seq = src.sequence().as_ref().to_vec();
    let qual = src.quality_scores().as_ref().to_vec();
    *out.sequence_mut() = seq[start..end].to_vec().into();
    *out.quality_scores_mut() = qual[start..end].to_vec().into();

    // Name suffix on splits.
    if total > 1 {
        let base = src.name().map(|n| n.to_vec()).unwrap_or_default();
        let mut name = base;
        name.extend_from_slice(format!("_segment_{}", idx + 1).as_bytes());
        *out.name_mut() = Some(name.into());
    }

    // Rebuild MM/ML/MN when the source carried modification tags. Only touch the
    // three tags when the source actually had `MM` (preserves prior behavior:
    // a source with ML/MN but no MM is left untouched).
    if matches!(
        src.data().get(&Tag::BASE_MODIFICATIONS),
        Some(Value::String(_))
    ) {
        let data = out.data_mut();
        match reconstruct_mods(src, &seq, start, end) {
            Some((mm_new, ml_new)) => {
                data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm_new.into()));
                match ml_new {
                    // Source had ML: write the sliced ML.
                    Some(ml) => {
                        data.insert(
                            Tag::BASE_MODIFICATION_PROBABILITIES,
                            Value::Array(Array::UInt8(ml)),
                        );
                    },
                    // Source was MM-only: drop the cloned ML so we stay MM-only
                    // rather than emit an empty (invalid) ML.
                    None => {
                        data.remove(&Tag::BASE_MODIFICATION_PROBABILITIES);
                    },
                }
                data.insert(
                    Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
                    Value::Int32((end - start) as i32),
                );
            },
            None => {
                data.remove(&Tag::BASE_MODIFICATIONS);
                data.remove(&Tag::BASE_MODIFICATION_PROBABILITIES);
                data.remove(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH);
            },
        }
    }

    // ONT signal tags (mv/ts/ns/sp/pi): rewrite for the window under
    // `--update-moves`, else drop on trim. Kept separate from the per-base pass
    // because it decodes the move table and updates several tags together.
    let orig_len = seq.len();
    for (tag, val) in signal_tag_updates(src, orig_len, start, end, total, update_moves) {
        match val {
            Some(v) => {
                out.data_mut().insert(tag, v);
            },
            None => {
                out.data_mut().remove(&tag);
            },
        }
    }

    if start != 0 || end != orig_len {
        // Drop position/signal tags we can't reconstruct (poly-A / barcode coords).
        for t in DROP_ON_TRIM_TAGS {
            out.data_mut().remove(&Tag::new(t[0], t[1]));
        }
        // Refresh qs (mean read qscore) from the trimmed quality — dorado
        // recomputes it per (sub)read — but only when the source carried one.
        if src.data().get(&Tag::new(b'q', b's')).is_some() {
            let qs = crate::qual::mean_prob_q(&qual[start..end]) as f32;
            out.data_mut()
                .insert(Tag::new(b'q', b's'), Value::Float(qs));
        }
    }

    // st/du describe the whole parent read; on a split they no longer fit the
    // subread (which starts later in the signal), so drop them.
    if total > 1 {
        for t in DROP_ON_SPLIT_TAGS {
            out.data_mut().remove(&Tag::new(t[0], t[1]));
        }
    }

    // Per-base arrays (PacBio ip/pw/fi/fp/ri/rp, or any read-length `B` array) are
    // sliced to the window so the trimmed record stays valid. MM/ML/MN, the signal
    // tags, and the dropped-on-trim tags are handled above and skipped here.
    if start != 0 || end != orig_len {
        let data = out.data_mut();
        let mut to_replace: Vec<(Tag, Value)> = Vec::new();
        for (tag, value) in data.iter() {
            let t = <[u8; 2]>::from(tag);
            // Skip every tag with dedicated handling above (mods, signal, poly-A,
            // and the drop-on-trim/split sets) so the structural per-base slicer
            // can't re-slice e.g. a kept `pa` array that happens to be read-length.
            if t == *b"MM"
                || t == *b"ML"
                || t == *b"MN"
                || SIGNAL_TAGS.contains(&t)
                || POLYA_TAGS.contains(&t)
                || DROP_ON_TRIM_TAGS.contains(&t)
                || DROP_ON_SPLIT_TAGS.contains(&t)
            {
                continue;
            }
            if let Some(v) = perbase_slice(value, orig_len, start, end) {
                to_replace.push((tag, v));
            }
        }
        for (tag, v) in to_replace {
            data.insert(tag, v);
        }
    }

    out
}

/// Slice MM/ML to the window `[start, end)` and re-serialize. Returns `None`
/// when the source has no `MM` tag, or when no modified position survives the
/// window (caller drops MM/ML/MN in that case). The inner `Option<Vec<u8>>` is
/// the sliced ML, or `None` when the source carried `MM` but no `ML` — ML is
/// optional per the SAM spec, so an MM-only record must stay MM-only rather than
/// gain a bogus empty ML. Shared by the BAM→BAM and BAM→FASTQ paths so they
/// cannot drift.
pub fn reconstruct_mods(
    src: &RecordBuf,
    seq: &[u8],
    start: usize,
    end: usize,
) -> Option<(Vec<u8>, Option<Vec<u8>>)> {
    let mm_raw = match src.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => s.to_vec(),
        _ => return None,
    };
    // Whether the source actually carried an ML array. If it didn't, we must not
    // emit one: declaring N modified positions with zero probabilities is an
    // invalid record that samtools/modkit reject.
    let ml_present = matches!(
        src.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES),
        Some(Value::Array(Array::UInt8(_)))
    );
    let ml_raw: Vec<u8> = match src.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        _ => Vec::new(),
    };
    let parsed = mods::parse(&mm_raw, &ml_raw);
    let sliced = mods::reconstruct(&parsed, seq, start, end);
    let (mm_new, ml_new) = mods::serialize(&sliced);
    if mm_new.is_empty() {
        None
    } else {
        Some((mm_new, ml_present.then_some(ml_new)))
    }
}

/// Single-threaded uBAM pipeline: refuse aligned reads, filter, trim, reconstruct.
fn run_bam_seq(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<RecordBuf>>,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    let mut malformed_tag_reads = 0u64;
    for rec in records {
        let rec = rec?;
        crate::io::bam::ensure_unaligned(&rec)?;
        counters.input_reads.fetch_add(1, Ordering::Relaxed);

        let seq = rec.sequence().as_ref().to_vec();
        counters
            .input_bases
            .fetch_add(seq.len() as u64, Ordering::Relaxed);
        let qual = rec.quality_scores().as_ref().to_vec();
        if qual.len() != seq.len() {
            let name = rec
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                .unwrap_or_else(|| "<unnamed>".to_string());
            anyhow::bail!(
                "read {name}: BAM record SEQ length {} != QUAL length {} \
                 (records without full per-base quality are not supported)",
                seq.len(),
                qual.len()
            );
        }
        if has_malformed_perbase_tag(&rec, seq.len()) {
            malformed_tag_reads += 1;
        }
        if !filter::passes(&seq, &qual, &cfg.filter) {
            continue;
        }
        let intervals = trim::apply(seq.len(), &qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        let mut out_bases = 0u64;
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            let out = reconstruct_record(&rec, s, e, total, idx, cfg.update_moves);
            sink.write_record(header, &out)?;
            counters.output_reads.fetch_add(1, Ordering::Relaxed);
            out_bases += (e - s) as u64;
        }
        counters
            .output_bases
            .fetch_add(out_bases, Ordering::Relaxed);
    }
    Ok(Stats {
        input_reads: counters.input_reads.load(Ordering::Relaxed),
        output_reads: counters.output_reads.load(Ordering::Relaxed),
        input_bases: counters.input_bases.load(Ordering::Relaxed),
        output_bases: counters.output_bases.load(Ordering::Relaxed),
        malformed_tag_reads,
    })
}

/// Shared parallel driver: reader iterator -> rayon pool (render) -> bounded
/// channel -> dedicated writer thread (write_one). Unordered. Mirrors
/// `run_fastq`'s error seam (see `pipeline/fastq.rs`): the first parse/render
/// error and the first write error are each captured in a `Mutex<Option<_>>`
/// slot; the writer task keeps draining `rx` after an error (never `break`) so
/// producer threads blocked on the bounded channel's `tx.send` can never
/// deadlock waiting for a writer that stopped consuming.
fn run_bam_parallel<T, S, Render, WriteOne>(
    records: impl Iterator<Item = anyhow::Result<RecordBuf>> + Send,
    cfg: &Config,
    sink: &mut S,
    render: Render,
    write_one: WriteOne,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats>
where
    T: Send,
    S: Send,
    // `Render` returns the surviving segments alongside their total base count
    // (segments can be `RecordBuf`s or rendered FASTQ byte buffers, so the base
    // count can't be recovered generically from `T` itself — the caller sums it
    // from the same intervals it renders from).
    Render: Fn(&RecordBuf, &Config) -> anyhow::Result<(Vec<T>, u64)> + Sync,
    WriteOne: Fn(&mut S, &T) -> std::io::Result<()> + Send,
{
    let render_workers = if cfg.render_workers >= 1 {
        cfg.render_workers
    } else {
        cfg.threads.max(1)
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(render_workers)
        .build()?;
    let (tx, rx) = crossbeam_channel::bounded::<T>(render_workers * 4);
    let malformed = AtomicU64::new(0);
    let proc_err: std::sync::Mutex<Option<anyhow::Error>> = std::sync::Mutex::new(None);
    let write_err: std::sync::Mutex<Option<std::io::Error>> = std::sync::Mutex::new(None);

    // Writer on a plain scoped OS thread; render on the budget-sized LOCAL rayon
    // pool via `pool.install` so the nested `par_bridge` uses THAT pool (not
    // rayon's global num_cpus pool) — this is what makes `-t` bound the render
    // threads. The writer keeps draining on a write error (never `break`) so
    // bounded-channel producers can't deadlock.
    std::thread::scope(|s| {
        let write_err_ref = &write_err;
        s.spawn(move || {
            let mut errored = false;
            for item in rx.iter() {
                if errored {
                    continue; // keep draining so bounded-channel producers never block
                }
                if let Err(e) = write_one(sink, &item) {
                    *write_err_ref.lock().unwrap() = Some(e);
                    errored = true;
                }
            }
        });

        pool.install(|| {
            records.par_bridge().for_each(|rec| {
                let rec = match rec {
                    Ok(r) => r,
                    Err(e) => {
                        let mut g = proc_err.lock().unwrap();
                        if g.is_none() {
                            *g = Some(e);
                        }
                        return;
                    },
                };
                counters.input_reads.fetch_add(1, Ordering::Relaxed);
                counters
                    .input_bases
                    .fetch_add(rec.sequence().as_ref().len() as u64, Ordering::Relaxed);
                if has_malformed_perbase_tag(&rec, rec.sequence().as_ref().len()) {
                    malformed.fetch_add(1, Ordering::Relaxed);
                }
                match render(&rec, cfg) {
                    Ok((items, out_bases)) => {
                        counters
                            .output_reads
                            .fetch_add(items.len() as u64, Ordering::Relaxed);
                        counters
                            .output_bases
                            .fetch_add(out_bases, Ordering::Relaxed);
                        for it in items {
                            let _ = tx.send(it);
                        }
                    },
                    Err(e) => {
                        let mut g = proc_err.lock().unwrap();
                        if g.is_none() {
                            *g = Some(e);
                        }
                    },
                }
            });
        });
        drop(tx);
    });

    if let Some(e) = proc_err.lock().unwrap().take() {
        return Err(e);
    }
    if let Some(e) = write_err.lock().unwrap().take() {
        return Err(e.into());
    }
    Ok(Stats {
        input_reads: counters.input_reads.load(Ordering::Relaxed),
        output_reads: counters.output_reads.load(Ordering::Relaxed),
        input_bases: counters.input_bases.load(Ordering::Relaxed),
        output_bases: counters.output_bases.load(Ordering::Relaxed),
        malformed_tag_reads: malformed.load(Ordering::Relaxed),
    })
}

/// Threads-aware uBAM pipeline entry point: refuse aligned reads, filter, trim,
/// reconstruct. Sequential for `cfg.threads <= 1`; otherwise renders each
/// record on a rayon work pool and drains the resulting `RecordBuf`s through
/// `run_bam_parallel`'s bounded channel onto a dedicated writer task. Output
/// order is unordered for `threads > 1` (records land in arrival order, not
/// input order) — the BAM format has no ordering requirement for uBAM.
pub fn run_bam(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<RecordBuf>> + Send,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    if cfg.threads <= 1 {
        return run_bam_seq(header, records, sink, cfg, counters);
    }
    run_bam_parallel(
        records,
        cfg,
        sink,
        // render: per-record guards + filter + trim + reconstruct -> (Vec<RecordBuf>, output bases)
        |rec, cfg| {
            crate::io::bam::ensure_unaligned(rec)?;
            let seq = rec.sequence().as_ref().to_vec();
            let qual = rec.quality_scores().as_ref().to_vec();
            if qual.len() != seq.len() {
                let name = rec
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_else(|| "<unnamed>".to_string());
                anyhow::bail!(
                    "read {name}: BAM record SEQ length {} != QUAL length {} \
                     (records without full per-base quality are not supported)",
                    seq.len(),
                    qual.len()
                );
            }
            if !filter::passes(&seq, &qual, &cfg.filter) {
                return Ok((Vec::new(), 0));
            }
            let intervals = trim::apply(seq.len(), &qual, &cfg.trim, cfg.filter.min_length);
            let total = intervals.len();
            let out_bases: u64 = intervals.iter().map(|&(s, e)| (e - s) as u64).sum();
            let items = intervals
                .into_iter()
                .enumerate()
                .map(|(idx, (s, e))| reconstruct_record(rec, s, e, total, idx, cfg.update_moves))
                .collect();
            Ok((items, out_bases))
        },
        // write_one: encode+write on the writer thread (bgzf compress is MT).
        |sink, rec| sink.write_record(header, rec),
        counters,
    )
}

/// Assemble the TAB-prefixed aux-tag block for one interval: carried non-mod
/// tags verbatim in source order, then the reconstructed MM/ML/MN block. Empty
/// when nothing is carried (caller writes a plain FASTQ record).
fn build_fastq_tags(
    src: &RecordBuf,
    seq: &[u8],
    start: usize,
    end: usize,
    total: usize,
    sel: &FastqTags,
) -> Vec<u8> {
    let mut tags = Vec::new();
    let orig_len = seq.len();
    let trimmed = start != 0 || end != orig_len;
    for (tag, value) in src.data().iter() {
        let t = <[u8; 2]>::from(tag);
        if t == *b"MM" || t == *b"ML" || t == *b"MN" {
            continue; // handled by the reconstructed block below
        }
        // On trim, drop the ONT signal tags (a sliced move table is impractical in
        // a FASTQ header, and signal-aware callers consume BAM — `--update-moves`
        // is BAM→BAM only) plus the poly-A and barcode-coordinate tags.
        if trimmed
            && (SIGNAL_TAGS.contains(&t)
                || POLYA_TAGS.contains(&t)
                || DROP_ON_TRIM_TAGS.contains(&t))
        {
            continue;
        }
        // On a split, st/du describe the parent read, not the subread.
        if total > 1 && DROP_ON_SPLIT_TAGS.contains(&t) {
            continue;
        }
        if !sel.carries(&t) {
            continue;
        }
        // Refresh qs from the trimmed quality (matches the BAM→BAM path).
        if t == *b"qs" && trimmed {
            let ql = src.quality_scores().as_ref();
            let qs = crate::qual::mean_prob_q(&ql[start..end]) as f32;
            tags.push(b'\t');
            tags.extend_from_slice(&format_aux_field(t, &Value::Float(qs)));
            continue;
        }
        // Per-base kinetics stay consistent with the trimmed sequence.
        tags.push(b'\t');
        match perbase_slice(value, orig_len, start, end) {
            Some(v) => tags.extend_from_slice(&format_aux_field(t, &v)),
            None => tags.extend_from_slice(&format_aux_field(t, value)),
        }
    }
    if sel.carries_mods()
        && let Some((mm, ml)) = reconstruct_mods(src, seq, start, end)
    {
        tags.push(b'\t');
        tags.extend_from_slice(&format_mods_aux(&mm, ml.as_deref(), end - start));
    }
    tags
}

/// Single-threaded uBAM→FASTQ pipeline: refuse aligned reads, filter, trim, then
/// write each surviving segment as FASTQ with the selected aux tags in the header
/// (MM/ML/MN reconstructed; others verbatim). gz compression, when requested, is
/// handled by the parallel `gzp` writer this drains into.
fn run_bam_to_fastq_seq<W>(
    records: impl Iterator<Item = anyhow::Result<RecordBuf>>,
    writer: &mut W,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats>
where
    W: Write,
{
    let mut malformed_tag_reads = 0u64;
    for rec in records {
        let rec = rec?;
        crate::io::bam::ensure_unaligned(&rec)?;
        counters.input_reads.fetch_add(1, Ordering::Relaxed);

        let seq = rec.sequence().as_ref().to_vec();
        counters
            .input_bases
            .fetch_add(seq.len() as u64, Ordering::Relaxed);
        let qual = rec.quality_scores().as_ref().to_vec();
        if qual.len() != seq.len() {
            let name = rec
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                .unwrap_or_else(|| "<unnamed>".to_string());
            anyhow::bail!(
                "read {name}: BAM record SEQ length {} != QUAL length {} \
                 (records without full per-base quality are not supported)",
                seq.len(),
                qual.len()
            );
        }
        if has_malformed_perbase_tag(&rec, seq.len()) {
            malformed_tag_reads += 1;
        }
        if !filter::passes(&seq, &qual, &cfg.filter) {
            continue;
        }
        let name = rec.name().map(|n| n.to_vec()).unwrap_or_default();
        let intervals = trim::apply(seq.len(), &qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        let mut out_bases = 0u64;
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            let tags = build_fastq_tags(&rec, &seq, s, e, total, &cfg.fastq_tags);
            if tags.is_empty() {
                write_segment(writer, &name, &seq[s..e], &qual[s..e], total, idx)?;
            } else {
                write_segment_tagged(writer, &name, &seq[s..e], &qual[s..e], total, idx, &tags)?;
            }
            counters.output_reads.fetch_add(1, Ordering::Relaxed);
            out_bases += (e - s) as u64;
        }
        counters
            .output_bases
            .fetch_add(out_bases, Ordering::Relaxed);
    }
    Ok(Stats {
        input_reads: counters.input_reads.load(Ordering::Relaxed),
        output_reads: counters.output_reads.load(Ordering::Relaxed),
        input_bases: counters.input_bases.load(Ordering::Relaxed),
        output_bases: counters.output_bases.load(Ordering::Relaxed),
        malformed_tag_reads,
    })
}

/// Threads-aware uBAM→FASTQ pipeline entry point: refuse aligned reads, filter,
/// trim, then write each surviving segment as FASTQ with the selected aux tags
/// in the header (MM/ML/MN reconstructed; others verbatim). Sequential for
/// `cfg.threads <= 1`; otherwise renders each record's FASTQ segments on a
/// rayon work pool and drains the resulting byte buffers through
/// `run_bam_parallel`'s bounded channel onto the writer, appended as-is.
/// Output order is unordered for `threads > 1` (records land in arrival order,
/// not input order).
pub fn run_bam_to_fastq<W: Write + Send>(
    records: impl Iterator<Item = anyhow::Result<RecordBuf>> + Send,
    writer: &mut W,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    if cfg.threads <= 1 {
        return run_bam_to_fastq_seq(records, writer, cfg, counters);
    }
    run_bam_parallel(
        records,
        cfg,
        writer,
        // render: guards + filter + trim -> (Vec<Vec<u8>>, output bases) (rendered FASTQ segments)
        |rec, cfg| {
            crate::io::bam::ensure_unaligned(rec)?;
            let seq = rec.sequence().as_ref().to_vec();
            let qual = rec.quality_scores().as_ref().to_vec();
            if qual.len() != seq.len() {
                let name = rec
                    .name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_else(|| "<unnamed>".to_string());
                anyhow::bail!(
                    "read {name}: BAM record SEQ length {} != QUAL length {} \
                     (records without full per-base quality are not supported)",
                    seq.len(),
                    qual.len()
                );
            }
            if !filter::passes(&seq, &qual, &cfg.filter) {
                return Ok((Vec::new(), 0));
            }
            let name = rec.name().map(|n| n.to_vec()).unwrap_or_default();
            let intervals = trim::apply(seq.len(), &qual, &cfg.trim, cfg.filter.min_length);
            let total = intervals.len();
            let mut out = Vec::with_capacity(total);
            let mut out_bases = 0u64;
            for (idx, (s, e)) in intervals.into_iter().enumerate() {
                let tags = build_fastq_tags(rec, &seq, s, e, total, &cfg.fastq_tags);
                let mut buf = Vec::new();
                if tags.is_empty() {
                    write_segment(&mut buf, &name, &seq[s..e], &qual[s..e], total, idx)?;
                } else {
                    write_segment_tagged(
                        &mut buf,
                        &name,
                        &seq[s..e],
                        &qual[s..e],
                        total,
                        idx,
                        &tags,
                    )?;
                }
                out.push(buf);
                out_bases += (e - s) as u64;
            }
            Ok((out, out_bases))
        },
        // write_one: append rendered bytes to the FastqOut writer.
        |w, buf| w.write_all(buf),
        counters,
    )
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    use super::*;

    fn ubam_with_mods(seq: &[u8], quals: Vec<u8>, mm: &[u8], ml: Vec<u8>) -> RecordBuf {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = seq.to_vec().into();
        *rec.quality_scores_mut() = quals.into();
        let data = rec.data_mut();
        data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.to_vec().into()));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(ml)),
        );
        rec
    }

    #[test]
    fn slices_seq_qual_and_rebuilds_tags() {
        // seq = C C A C ; C+m modified at C occ 0 and 2 -> pos 0 and 3; ML [10,20].
        let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0,1;", vec![10, 20]);
        // keep window [2,4): seq "AC", the modified C at pos3 survives (occ within window idx0).
        let out = reconstruct_record(&src, 2, 4, 1, 0, false);

        assert_eq!(out.sequence().as_ref(), b"AC");
        assert_eq!(out.quality_scores().as_ref(), &[32, 33]);

        let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => s.to_vec(),
            _ => panic!("no MM"),
        };
        assert_eq!(mm, b"C+m,0;");
        let ml = match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) => v.clone(),
            _ => panic!("no ML"),
        };
        assert_eq!(ml, vec![20]);
        // MN updated to the output length.
        let mn = match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(n)) => *n,
            _ => panic!("no MN"),
        };
        assert_eq!(mn, 2);
    }

    /// Regression test for a QUAL-absent uBAM record: `quality_scores()`
    /// decodes to empty while `sequence()` keeps its bases, so `seq.len() !=
    /// qual.len()`. Pre-fix, `run_bam` fed this straight into `trim::apply`,
    /// which slices `phred[start..end]` using `seq.len()`-derived bounds and
    /// panics (via `debug_assert_eq!` in debug builds, or an out-of-bounds
    /// slice panic in release) instead of returning an `Err`. Post-fix, the
    /// length-mismatch guard in `run_bam` must bail with an `Err` before any
    /// of that runs.
    #[test]
    fn qual_seq_length_mismatch_errors_without_panicking() {
        use crate::config::IoConfig;
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = b"ACGT".to_vec().into();
        // quality_scores left at its default (empty) -> SEQ/QUAL length mismatch.

        let header = sam::Header::default();
        let dir = tempfile::tempdir().unwrap();
        let mut sink =
            crate::io::bam::writer(Some(&dir.path().join("o.bam")), &header, 1, 6).unwrap();

        let cfg = Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
            render_workers: 0,
            compression_level: 6,
            update_moves: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        };

        let result = run_bam(
            &header,
            [Ok(rec)].into_iter(),
            &mut sink,
            &cfg,
            &Arc::new(Counters::default()),
        );
        assert!(
            result.is_err(),
            "SEQ/QUAL length mismatch must error, not panic"
        );
    }

    /// Regression test for the outer gate in `reconstruct_record`: a spec-invalid
    /// `MM` tag (typed as anything other than `Value::String`, e.g. `Int32`) must
    /// leave `MM`/`ML`/`MN` completely untouched. Pre-fix, the gate was a bare
    /// `.is_some()` check, which let this record enter the mod-rebuild block;
    /// `reconstruct_mods`'s `Some(Value::String(s)) => ... , _ => return None`
    /// match then falls through to `None` for a non-string MM, and the `None`
    /// branch in `reconstruct_record` REMOVES all three tags instead of leaving
    /// them alone. The fix narrows the gate to `matches!(.., Some(Value::String(_)))`
    /// so a spec-invalid MM never enters the rebuild block at all.
    #[test]
    fn reconstruct_record_leaves_non_string_mm_untouched() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        let data = src.data_mut();
        data.insert(Tag::BASE_MODIFICATIONS, Value::Int32(5)); // spec-invalid MM
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(vec![1, 2, 3])),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));

        let out = reconstruct_record(&src, 1, 4, 1, 0, false);

        match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::Int32(5)) => {},
            other => panic!("MM must be left untouched for a non-string value, got {other:?}"),
        }
        match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) if v == &[1u8, 2, 3] => {},
            other => panic!("ML must be left untouched, got {other:?}"),
        }
        match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(4)) => {},
            other => panic!("MN must be left untouched, got {other:?}"),
        }
    }

    #[test]
    fn split_suffixes_name_and_drops_empty_mods() {
        let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0;", vec![10]); // mod at pos0
        // segment [2,4) has no surviving C mod -> MM/ML removed entirely.
        let out = reconstruct_record(&src, 2, 4, 2, 1, false);
        // `.as_ref()` is ambiguous on `&BStr` (impls `AsRef<[u8]>` and `AsRef<BStr>`);
        // disambiguate with a turbofish (noodles-sam 0.85 / bstr 1.x adjustment).
        assert_eq!(AsRef::<[u8]>::as_ref(out.name().unwrap()), b"r1_segment_2");
        assert!(out.data().get(&Tag::BASE_MODIFICATIONS).is_none());
        assert!(
            out.data()
                .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
                .is_none()
        );
        assert!(
            out.data()
                .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
                .is_none()
        );
    }

    use crate::config::{FastqTags, IoConfig};
    use crate::filter::FilterConfig;
    use crate::qual::QualMode;
    use crate::trim::{QualityOp, TrimPlan};

    fn cfg_bam2fq(quality: Option<QualityOp>, head: usize, tags: FastqTags) -> Config {
        Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head,
                tail: 0,
                quality,
            },
            threads: 1,
            fastq_tags: tags,
            render_workers: 0,
            compression_level: 6,
            update_moves: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        }
    }

    // "CCACCCAC" C at seq idx 0,1,3,4,5,7; MM "C+m,0,1,0" -> occ 0,2,3 -> abs 0,3,4,
    // ML [10,20,30]. head-crop 2 -> window [2,8): keeps abs 3,4 renumbered -> "C+m,0,0;" ML [20,30] MN 6.
    fn read2_with_mods_and_rg() -> RecordBuf {
        let mut rec = ubam_with_mods(b"CCACCCAC", vec![35; 8], b"C+m,0,1,0;", vec![10, 20, 30]);
        rec.data_mut()
            .insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
        rec.data_mut()
            .insert(Tag::READ_GROUP, Value::String(b"grp1".as_slice().into()));
        rec
    }

    #[test]
    fn bam2fq_all_carries_rg_and_reconstructed_mods() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::All);
        let mut out = Vec::new();
        let stats = run_bam_to_fastq(
            [Ok(read2_with_mods_and_rg())].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
        let s = String::from_utf8(out).unwrap();
        // header carries RG verbatim + reconstructed mod block; seq head-cropped by 2.
        assert!(
            s.starts_with("@r1\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6\n"),
            "got: {s:?}"
        );
        assert!(s.contains("\nACCCAC\n+\n"), "cropped seq wrong: {s:?}");
    }

    #[test]
    fn bam2fq_only_mm_ml_drops_rg() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::parse("MM,ML").unwrap());
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(read2_with_mods_and_rg())].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(!s.contains("RG:Z"), "RG must be dropped: {s:?}");
        assert!(
            s.contains("MM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"),
            "mods missing: {s:?}"
        );
    }

    #[test]
    fn bam2fq_none_is_plain_fastq() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::None);
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(read2_with_mods_and_rg())].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(out, b"@r1\nACCCAC\n+\nDDDDDD\n"); // 35+33 = 'D'
    }

    #[test]
    fn bam2fq_split_suffixes_and_segments_mods() {
        // split at the low-qual base; each segment gets its own reconstructed mods.
        let cfg = cfg_bam2fq(
            Some(QualityOp::Split {
                cutoff: 20,
                window: 1,
            }),
            0,
            FastqTags::All,
        );
        // seq CCAC, C+m at occ 0 and 2 -> abs 0,3; qual: good good BAD good so split [0,2),[3,4)
        let mut rec = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0,1;", vec![100, 200]);
        rec.data_mut()
            .insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));
        let mut out = Vec::new();
        let stats = run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(stats.output_reads, 2);
        let s = String::from_utf8(out).unwrap();
        // segment 1 = [0,2) "CC" keeps abs-0 mod; segment 2 = [3,4) "C" keeps abs-3 mod.
        assert!(
            s.contains("@r1_segment_1\tMM:Z:C+m,0;\tML:B:C,100\tMN:i:2"),
            "seg1: {s:?}"
        );
        assert!(
            s.contains("@r1_segment_2\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:1"),
            "seg2: {s:?}"
        );
    }

    #[test]
    fn bam2fq_no_mods_read_is_plain() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"plain".into());
        *rec.sequence_mut() = b"ACGT".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 4].into();
        let cfg = cfg_bam2fq(None, 0, FastqTags::All);
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(out, b"@plain\nACGT\n+\nIIII\n");
    }

    /// Regression test: a source uBAM carrying `MM` but NO `ML` (valid — ML is
    /// optional per the SAM spec) must be rewritten as an MM-only record, NOT gain
    /// an empty `ML:B:C`. Pre-fix, `reconstruct_record` inserted `Array::UInt8([])`
    /// unconditionally, producing an MM that declares modified positions with zero
    /// probabilities — a record samtools/modkit reject as invalid.
    #[test]
    fn reconstruct_record_mm_without_ml_stays_mm_only() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"CCAC".to_vec().into(); // C at 0,1,3
        *src.quality_scores_mut() = vec![40; 4].into();
        src.data_mut().insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(b"C+m,0,1;".to_vec().into()),
        );
        // deliberately no ML, no MN.

        let out = reconstruct_record(&src, 0, 4, 1, 0, false);

        // MM retained (both modified Cs are in-window -> "C+m,0,1;").
        let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => s.to_vec(),
            other => panic!("expected MM retained, got {other:?}"),
        };
        assert_eq!(mm, b"C+m,0,1;");
        // ML must be ABSENT — never an empty array.
        assert!(
            out.data()
                .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
                .is_none(),
            "MM-only source must not gain an ML tag"
        );
        // MN set to the window length.
        match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(4)) => {},
            other => panic!("expected MN=4, got {other:?}"),
        }
    }

    /// Companion for the BAM→FASTQ path: an MM-only source must emit a FASTQ
    /// header with `MM` + `MN` but no `ML:B:C` field.
    #[test]
    fn bam2fq_mm_without_ml_omits_ml_field() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = b"CCAC".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 4].into();
        rec.data_mut().insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(b"C+m,0,1;".to_vec().into()),
        );

        let cfg = cfg_bam2fq(None, 0, FastqTags::All);
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(
            s.contains("MM:Z:C+m,0,1;\tMN:i:4"),
            "expected MM+MN, got: {s:?}"
        );
        assert!(
            !s.contains("ML:B"),
            "MM-only record must not emit an ML field: {s:?}"
        );
    }

    #[test]
    fn reconstruct_record_slices_kinetics_and_drops_mv() {
        // PacBio-style per-base kinetics (ip/pw, length == read length) must be
        // sliced with the sequence; ONT `mv` (signal-space) must be dropped on
        // trim; per-read RG must ride through unchanged.
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 6].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![10, 11, 12, 13, 14, 15])),
        );
        d.insert(
            Tag::new(b'p', b'w'),
            Value::Array(Array::UInt16(vec![20, 21, 22, 23, 24, 25])),
        );
        d.insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![5, 1, 0, 1, 0])),
        );
        d.insert(Tag::READ_GROUP, Value::String(b"grp".as_slice().into()));

        // window [2,5): "GTA" (head-crop 2, tail-crop 1).
        let out = reconstruct_record(&src, 2, 5, 1, 0, false);
        assert_eq!(out.sequence().as_ref(), b"GTA");
        match out.data().get(&Tag::new(b'i', b'p')) {
            Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[12, 13, 14]),
            other => panic!("ip should be sliced [2,5): {other:?}"),
        }
        match out.data().get(&Tag::new(b'p', b'w')) {
            Some(Value::Array(Array::UInt16(v))) => assert_eq!(v, &[22, 23, 24]),
            other => panic!("pw should be sliced [2,5): {other:?}"),
        }
        assert!(
            out.data().get(&Tag::new(b'm', b'v')).is_none(),
            "mv must be dropped on trim"
        );
        assert!(
            matches!(out.data().get(&Tag::READ_GROUP), Some(Value::String(_))),
            "RG kept"
        );
    }

    #[test]
    fn reconstruct_record_slices_unknown_read_length_array_but_not_others() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        // Unknown B-array whose length == read length -> sliced structurally.
        src.data_mut().insert(
            Tag::new(b'z', b'z'),
            Value::Array(Array::Int32(vec![1, 2, 3, 4])),
        );
        // B-array whose length != read length -> not per-base, left alone.
        src.data_mut()
            .insert(Tag::new(b'x', b'y'), Value::Array(Array::UInt8(vec![9, 9])));

        let out = reconstruct_record(&src, 1, 3, 1, 0, false); // window [1,3)
        match out.data().get(&Tag::new(b'z', b'z')) {
            Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 3]),
            other => panic!("zz sliced: {other:?}"),
        }
        match out.data().get(&Tag::new(b'x', b'y')) {
            Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[9, 9]),
            other => panic!("xy untouched: {other:?}"),
        }
    }

    #[test]
    fn reconstruct_record_untrimmed_keeps_kinetics_and_mv() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        src.data_mut().insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![1, 2, 3, 4])),
        );
        src.data_mut().insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![5, 1, 1])),
        );

        // Full window [0,4): nothing trimmed -> everything preserved verbatim.
        let out = reconstruct_record(&src, 0, 4, 1, 0, false);
        match out.data().get(&Tag::new(b'i', b'p')) {
            Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[1, 2, 3, 4]),
            other => panic!("ip unchanged: {other:?}"),
        }
        assert!(
            out.data().get(&Tag::new(b'm', b'v')).is_some(),
            "mv kept when untrimmed"
        );
    }

    #[test]
    fn bam2fq_slices_kinetics_in_header() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = b"ACGTAC".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 6].into();
        rec.data_mut().insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![10, 11, 12, 13, 14, 15])),
        );

        let cfg = cfg_bam2fq(None, 2, FastqTags::All); // head-crop 2 -> window [2,6)
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(
            s.contains("ip:B:C,12,13,14,15"),
            "kinetics not sliced in FASTQ header: {s:?}"
        );
    }

    #[test]
    fn malformed_perbase_tag_detected_and_left_untouched() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        // ip length 3 != read length 4 -> malformed known per-base tag.
        src.data_mut().insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![1, 2, 3])),
        );

        assert!(has_malformed_perbase_tag(&src, 4));
        // Can't safely slice it -> left exactly as-is.
        let out = reconstruct_record(&src, 1, 3, 1, 0, false);
        match out.data().get(&Tag::new(b'i', b'p')) {
            Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[1, 2, 3]),
            other => panic!("malformed ip left as-is: {other:?}"),
        }
    }

    // A synthetic move table: stride 2, 6 ones (one per base) at block indices
    // 0,1,3,4,6,7 -> 8 blocks total. Shared by the --update-moves tests below.
    fn ubam_with_moves() -> RecordBuf {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 6].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![2, 1, 1, 0, 1, 1, 0, 1, 1])),
        );
        d.insert(Tag::new(b't', b's'), Value::Int32(10));
        // Consistent: ts0 + n_blocks*stride = 10 + 8*2 = 26.
        d.insert(Tag::new(b'n', b's'), Value::Int32(26));
        d.insert(
            Tag::new(b's', b't'),
            Value::String(b"2024-06-21T10:00:00Z".as_slice().into()),
        );
        d.insert(Tag::new(b'd', b'u'), Value::Float(5.0));
        src
    }

    #[test]
    fn update_moves_head_crop_slices_mv_bumps_ts_keeps_ns() {
        // head-crop 2 -> window [2,6): block_first=ones[2]=3, block_second=8.
        let out = reconstruct_record(&ubam_with_moves(), 2, 6, 1, 0, true);
        assert_eq!(out.sequence().as_ref(), b"GTAC");
        assert_eq!(
            AsRef::<[u8]>::as_ref(out.name().unwrap()),
            b"r1",
            "crop keeps the read name"
        );
        // mv = [stride] + moves[3 .. 8] = [2] + [1,1,0,1,1].
        match out.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1, 1]),
            other => panic!("mv: {other:?}"),
        }
        // ts += block_first*stride = 10 + 3*2 = 16; ns = ts + span = 16 + (8-3)*2 = 26
        // (a head-only crop leaves ns unchanged).
        match out.data().get(&Tag::new(b't', b's')) {
            Some(Value::Int32(16)) => {},
            other => panic!("ts: {other:?}"),
        }
        match out.data().get(&Tag::new(b'n', b's')) {
            Some(Value::Int32(26)) => {},
            other => panic!("ns: {other:?}"),
        }
        assert!(out.data().get(&Tag::new(b's', b'p')).is_none());
        assert!(out.data().get(&Tag::new(b'p', b'i')).is_none());
        // A crop keeps the read identity, so st/du stay.
        assert!(
            out.data().get(&Tag::new(b's', b't')).is_some(),
            "st kept on crop"
        );
        assert!(
            out.data().get(&Tag::new(b'd', b'u')).is_some(),
            "du kept on crop"
        );
    }

    #[test]
    fn update_moves_tail_crop_shrinks_ns() {
        // tail-crop 2 -> window [0,4): block_first=ones[0]=0, block_second=ones[4]=6.
        let out = reconstruct_record(&ubam_with_moves(), 0, 4, 1, 0, true);
        // mv = [stride] + moves[0 .. 6] = [2] + [1,1,0,1,1,0].
        match out.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1, 1, 0]),
            other => panic!("mv: {other:?}"),
        }
        // ts unchanged (no head trim): 10; ns = ts + span = 10 + (6-0)*2 = 22 (< 26).
        match out.data().get(&Tag::new(b't', b's')) {
            Some(Value::Int32(10)) => {},
            other => panic!("ts: {other:?}"),
        }
        match out.data().get(&Tag::new(b'n', b's')) {
            Some(Value::Int32(22)) => {},
            other => panic!("ns must shrink on a tail crop (dorado ns = trim + span): {other:?}"),
        }
    }

    #[test]
    fn update_moves_split_emits_subread_tags() {
        // Split into [0,3) and [3,6): each is a dorado-style subread.
        let s1 = reconstruct_record(&ubam_with_moves(), 0, 3, 2, 0, true);
        assert_eq!(AsRef::<[u8]>::as_ref(s1.name().unwrap()), b"r1_segment_1");
        // mv = [2] + moves[ones[0]=0 .. ones[3]=4] = [2] + [1,1,0,1].
        match s1.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1]),
            other => panic!("s1 mv: {other:?}"),
        }
        match s1.data().get(&Tag::new(b't', b's')) {
            Some(Value::Int32(0)) => {},
            o => panic!("s1 ts should be 0: {o:?}"),
        }
        match s1.data().get(&Tag::new(b'n', b's')) {
            Some(Value::Int32(8)) => {}, // (block 4-0)*stride 2
            o => panic!("s1 ns: {o:?}"),
        }
        match s1.data().get(&Tag::new(b's', b'p')) {
            Some(Value::Int32(0)) => {}, // block_first 0 * stride
            o => panic!("s1 sp: {o:?}"),
        }
        match s1.data().get(&Tag::new(b'p', b'i')) {
            Some(Value::String(s)) => assert_eq!(s.to_vec(), b"r1"),
            o => panic!("s1 pi: {o:?}"),
        }
        // dorado marks split products with read_number -1.
        match s1.data().get(&Tag::new(b'r', b'n')) {
            Some(Value::Int32(-1)) => {},
            o => panic!("s1 rn should be -1: {o:?}"),
        }
        // st/du describe the parent read -> dropped on a split subread.
        assert!(
            s1.data().get(&Tag::new(b's', b't')).is_none(),
            "st dropped on split"
        );
        assert!(
            s1.data().get(&Tag::new(b'd', b'u')).is_none(),
            "du dropped on split"
        );

        let s2 = reconstruct_record(&ubam_with_moves(), 3, 6, 2, 1, true);
        assert_eq!(AsRef::<[u8]>::as_ref(s2.name().unwrap()), b"r1_segment_2");
        // mv = [2] + moves[ones[3]=4 .. 8] = [2] + [1,0,1,1].
        match s2.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 0, 1, 1]),
            other => panic!("s2 mv: {other:?}"),
        }
        match s2.data().get(&Tag::new(b'n', b's')) {
            Some(Value::Int32(8)) => {}, // (8-4)*2
            o => panic!("s2 ns: {o:?}"),
        }
        match s2.data().get(&Tag::new(b's', b'p')) {
            Some(Value::Int32(8)) => {}, // block_first 4 * stride 2
            o => panic!("s2 sp: {o:?}"),
        }
    }

    #[test]
    fn default_drops_all_signal_tags_on_trim() {
        let mut src = ubam_with_moves();
        src.data_mut().insert(Tag::new(b's', b'p'), Value::Int32(5));
        src.data_mut().insert(
            Tag::new(b'p', b'i'),
            Value::String(b"parent".as_slice().into()),
        );

        // update_moves = false + trimmed -> mv/ts/ns/sp/pi all removed.
        let out = reconstruct_record(&src, 2, 6, 1, 0, false);
        for t in [b"mv", b"ts", b"ns", b"sp", b"pi"] {
            assert!(
                out.data().get(&Tag::new(t[0], t[1])).is_none(),
                "{} must be dropped by default on trim",
                std::str::from_utf8(t).unwrap()
            );
        }
    }

    #[test]
    fn trim_drops_polya_barcode_tags_and_refreshes_qs() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        // First two bases low quality (phred 2), the rest Q40.
        *src.quality_scores_mut() = vec![2, 2, 40, 40, 40, 40].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![100, 200, 300, 400, 500])),
        );
        d.insert(Tag::new(b'p', b't'), Value::Int32(50));
        d.insert(
            Tag::new(b'b', b'i'),
            Value::Array(Array::Float(vec![0.9, 5.0, 20.0])),
        );
        d.insert(Tag::new(b'q', b's'), Value::Float(20.0)); // whole-read qs (stale after crop)
        d.insert(
            Tag::new(b'R', b'G'),
            Value::String(b"grp".as_slice().into()),
        );

        // head-crop 2 -> window [2,6): keeps only the Q40 bases.
        let out = reconstruct_record(&src, 2, 6, 1, 0, false);

        // Unreconstructable poly-A / barcode coordinate tags are dropped.
        for t in [b"pa", b"pt", b"bi"] {
            assert!(
                out.data().get(&Tag::new(t[0], t[1])).is_none(),
                "{} must be dropped on trim",
                std::str::from_utf8(t).unwrap()
            );
        }
        // qs is recomputed from the trimmed (all-Q40) quality, not left at 20.
        match out.data().get(&Tag::new(b'q', b's')) {
            Some(Value::Float(q)) => {
                let expected = crate::qual::mean_prob_q(&[40, 40, 40, 40]) as f32;
                assert!(
                    (q - expected).abs() < 1e-4,
                    "qs recomputed: got {q}, want {expected}"
                );
            },
            other => panic!("qs: {other:?}"),
        }
        // Per-read metadata (RG) is untouched.
        assert!(out.data().get(&Tag::new(b'R', b'G')).is_some());
    }

    // ubam_with_moves head-crop 2 spans original-signal window [ts0+3*2, ts0+8*2]
    // = [16, 26]; a split segment [3,6) spans [ts0+4*2, 26] = [18, 26].
    #[test]
    fn update_moves_crop_keeps_polya_when_tail_survives() {
        let mut src = ubam_with_moves();
        // anchor + boundary all inside [16,26]; the split range is a sentinel.
        src.data_mut().insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![20, 18, 24, -1, -1])),
        );
        src.data_mut()
            .insert(Tag::new(b'p', b't'), Value::Int32(30));

        let out = reconstruct_record(&src, 2, 6, 1, 0, true); // head-crop 2
        // A crop keeps the read identity + POD5 signal, so absolute pa stays valid.
        match out.data().get(&Tag::new(b'p', b'a')) {
            Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[20, 18, 24, -1, -1]),
            other => panic!("pa should be kept as-is on a crop: {other:?}"),
        }
        match out.data().get(&Tag::new(b'p', b't')) {
            Some(Value::Int32(30)) => {},
            other => panic!("pt: {other:?}"),
        }
    }

    #[test]
    fn update_moves_split_shifts_polya_into_subread_frame() {
        let mut src = ubam_with_moves();
        src.data_mut().insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![20, 18, 24, -1, -1])),
        );
        src.data_mut()
            .insert(Tag::new(b'p', b't'), Value::Int32(30));

        // split segment [3,6): kept signal window [18,26] -> shift real positions by -18.
        let out = reconstruct_record(&src, 3, 6, 2, 1, true);
        match out.data().get(&Tag::new(b'p', b'a')) {
            Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 0, 6, -1, -1]),
            other => panic!("pa should shift into the subread frame: {other:?}"),
        }
        match out.data().get(&Tag::new(b'p', b't')) {
            Some(Value::Int32(30)) => {}, // base count unchanged
            other => panic!("pt: {other:?}"),
        }
    }

    #[test]
    fn update_moves_drops_polya_when_tail_trimmed() {
        let mut src = ubam_with_moves();
        // anchor at 12 sits in the trimmed-off front signal (kept window is [16,26]).
        src.data_mut().insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![12, 10, 14, -1, -1])),
        );
        src.data_mut()
            .insert(Tag::new(b'p', b't'), Value::Int32(30));

        let out = reconstruct_record(&src, 2, 6, 1, 0, true); // head-crop 2
        assert!(
            out.data().get(&Tag::new(b'p', b'a')).is_none(),
            "pa dropped when tail trimmed"
        );
        assert!(
            out.data().get(&Tag::new(b'p', b't')).is_none(),
            "pt dropped when tail trimmed"
        );
    }

    #[test]
    fn update_moves_does_not_reslice_read_length_pa() {
        // Regression (review F1): a pa array whose length happens to equal the read
        // length must NOT be treated as a per-base array and sliced.
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTA".to_vec().into(); // 5 bases
        *src.quality_scores_mut() = vec![40; 5].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![2, 1, 1, 1, 1, 1])),
        ); // stride 2, 5 ones
        d.insert(Tag::new(b't', b's'), Value::Int32(0));
        d.insert(Tag::new(b'n', b's'), Value::Int32(10));
        // 5-element pa (== read length), all real positions inside the kept window.
        d.insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![4, 2, 6, -1, -1])),
        );

        // head-crop 1 -> window [1,5): kept signal window [2,10]; pa survives.
        let out = reconstruct_record(&src, 1, 5, 1, 0, true);
        match out.data().get(&Tag::new(b'p', b'a')) {
            Some(Value::Array(Array::Int32(v))) => {
                assert_eq!(v, &[4, 2, 6, -1, -1], "pa must not be re-sliced")
            },
            other => panic!("pa: {other:?}"),
        }
    }

    #[test]
    fn update_moves_without_move_table_drops_signal_and_polya() {
        // --update-moves but no move table -> can't relate signal to sequence, so
        // the signal + poly-A tags are dropped (parse_move_table -> None -> drop_all).
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 6].into();
        let d = src.data_mut();
        d.insert(Tag::new(b't', b's'), Value::Int32(10));
        d.insert(Tag::new(b'n', b's'), Value::Int32(100));
        d.insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![20, 18, 24, -1, -1])),
        );
        d.insert(Tag::new(b'p', b't'), Value::Int32(30));

        let out = reconstruct_record(&src, 2, 6, 1, 0, true);
        for t in [b"ts", b"ns", b"pa", b"pt"] {
            assert!(
                out.data().get(&Tag::new(t[0], t[1])).is_none(),
                "{} dropped when the move table is absent",
                std::str::from_utf8(t).unwrap()
            );
        }
    }

    #[test]
    fn update_moves_polya_boundary_end_inclusive_anchor_exclusive() {
        // split [3,6): kept window [18,26). A range END exactly at kept_end
        // (exclusive) survives; an anchor at kept_end is out of the window -> drop.
        let mk = |pa: Vec<i32>| {
            let mut src = ubam_with_moves();
            src.data_mut()
                .insert(Tag::new(b'p', b'a'), Value::Array(Array::Int32(pa)));
            src.data_mut()
                .insert(Tag::new(b'p', b't'), Value::Int32(30));
            src
        };
        // range end == kept_end (26) -> survives, shifted by -18.
        let kept = reconstruct_record(&mk(vec![20, 18, 26, -1, -1]), 3, 6, 2, 1, true);
        match kept.data().get(&Tag::new(b'p', b'a')) {
            Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 0, 8, -1, -1]),
            other => panic!("range-end at kept_end should survive: {other:?}"),
        }
        // anchor == kept_end (26) -> out of window -> dropped.
        let dropped = reconstruct_record(&mk(vec![26, 18, 24, -1, -1]), 3, 6, 2, 1, true);
        assert!(
            dropped.data().get(&Tag::new(b'p', b'a')).is_none(),
            "anchor at exclusive boundary drops"
        );
    }

    #[test]
    fn run_bam_parallel_matches_sequential_as_multiset() {
        use crate::config::{FastqTags, IoConfig};
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        let mk = |threads| Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 2,
                tail: 2,
                quality: None,
            },
            threads,
            fastq_tags: FastqTags::All,
            render_workers: 0,
            compression_level: 6,
            update_moves: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        };
        // 300 reads with mods so reconstruction runs on every one.
        let recs: Vec<RecordBuf> = (0..300)
            .map(|_| ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]))
            .collect();

        let header = sam::Header::default();
        let decode = |bytes: &[u8]| -> Vec<(Vec<u8>, Vec<u8>)> {
            // (seq, MM-bytes) pairs, sorted, as an order-independent fingerprint.
            let mut r = noodles_bam::io::Reader::new(bytes);
            let h = r.read_header().unwrap();
            let mut out = Vec::new();
            let mut buf = RecordBuf::default();
            while r.read_record_buf(&h, &mut buf).unwrap() != 0 {
                let seq = buf.sequence().as_ref().to_vec();
                let mm = match buf.data().get(&Tag::BASE_MODIFICATIONS) {
                    Some(Value::String(s)) => s.to_vec(),
                    _ => Vec::new(),
                };
                out.push((seq, mm));
            }
            out.sort();
            out
        };

        // t1 -> single-threaded bgzf sink, written to a tempfile.
        let dir = tempfile::tempdir().unwrap();
        let p1 = dir.path().join("t1.bam");
        let mut sink1 = crate::io::bam::writer(Some(&p1), &header, 1, 6).unwrap();
        run_bam(
            &header,
            recs.clone().into_iter().map(anyhow::Ok),
            &mut sink1,
            &mk(1),
            &Arc::new(Counters::default()),
        )
        .unwrap();
        sink1.finish().unwrap();
        let b1 = std::fs::read(&p1).unwrap();

        // t8 -> MT sink to a tempfile (MT writer needs an owned Write + Send).
        let p8 = dir.path().join("t8.bam");
        let mut sink8 = crate::io::bam::writer(Some(&p8), &header, 8, 6).unwrap();
        run_bam(
            &header,
            recs.into_iter().map(anyhow::Ok),
            &mut sink8,
            &mk(8),
            &Arc::new(Counters::default()),
        )
        .unwrap();
        sink8.finish().unwrap();
        let b8 = std::fs::read(&p8).unwrap();

        assert_eq!(
            decode(&b1),
            decode(&b8),
            "t1 and t8 must produce the same record set"
        );
    }

    /// Mirrors `pipeline::fastq`'s `parallel_surfaces_write_error_without_deadlock`,
    /// but drives `run_bam_parallel` directly with a stub sink whose `write_one`
    /// starts erroring after `limit` writes. Record count (3000) exceeds the
    /// bounded channel capacity (`threads * 4` = 16), so a pre-fix build that
    /// stops draining `rx` on the first write error would deadlock instead of
    /// returning.
    #[test]
    fn run_bam_parallel_surfaces_write_error_without_deadlock() {
        use std::io;

        use crate::config::IoConfig;
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        struct FailAfter {
            limit: usize,
            written: usize,
        }

        let cfg = Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            threads: 4,
            fastq_tags: crate::config::FastqTags::All,
            render_workers: 0,
            compression_level: 6,
            update_moves: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        };
        let recs: Vec<anyhow::Result<RecordBuf>> = (0..3000)
            .map(|_| anyhow::Ok(RecordBuf::default()))
            .collect();

        let mut sink = FailAfter {
            limit: 100,
            written: 0,
        };
        let res = run_bam_parallel(
            recs.into_iter(),
            &cfg,
            &mut sink,
            |_rec, _cfg| anyhow::Ok((vec![()], 0)),
            |sink, _item: &()| -> io::Result<()> {
                if sink.written >= sink.limit {
                    return Err(io::Error::new(io::ErrorKind::BrokenPipe, "boom"));
                }
                sink.written += 1;
                Ok(())
            },
            &Arc::new(Counters::default()),
        );
        assert!(
            res.is_err(),
            "write error must surface as Err, and must not hang"
        );
    }

    /// Mirrors `pipeline::fastq`'s `parallel_surfaces_parse_error_instead_of_dropping_it`,
    /// driving `run_bam_parallel` directly so a malformed upstream record (an
    /// `Err` item from the input iterator) is not silently swallowed.
    #[test]
    fn run_bam_parallel_surfaces_parse_error_instead_of_dropping_it() {
        use std::io;

        use crate::config::IoConfig;
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        struct NullSink;

        let cfg = Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            threads: 4,
            fastq_tags: crate::config::FastqTags::All,
            render_workers: 0,
            compression_level: 6,
            update_moves: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        };
        let good: Vec<anyhow::Result<RecordBuf>> =
            (0..5).map(|_| anyhow::Ok(RecordBuf::default())).collect();
        let recs = good
            .into_iter()
            .chain(std::iter::once(Err(anyhow::anyhow!("bad record"))));

        let mut sink = NullSink;
        let res = run_bam_parallel(
            recs,
            &cfg,
            &mut sink,
            |_rec, _cfg| anyhow::Ok((vec![()], 0)),
            |_sink: &mut NullSink, _item: &()| -> io::Result<()> { Ok(()) },
            &Arc::new(Counters::default()),
        );
        assert!(
            res.is_err(),
            "a malformed record must not be silently dropped on the parallel path"
        );
    }

    #[test]
    fn run_bam_to_fastq_parallel_matches_sequential_as_multiset() {
        use crate::config::{FastqTags, IoConfig};
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        let mk = |threads| Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 2,
                tail: 2,
                quality: None,
            },
            threads,
            fastq_tags: FastqTags::All,
            render_workers: 0,
            compression_level: 6,
            update_moves: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        };
        let recs: Vec<RecordBuf> = (0..300)
            .map(|_| ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]))
            .collect();

        let sorted_records = |bytes: &[u8]| {
            let s = String::from_utf8(bytes.to_vec()).unwrap();
            // Group every 4 consecutive lines into one record rather than
            // splitting on '@': a QUAL byte of Phred-31 (ASCII '@') would
            // corrupt an '@'-split. FASTQ records are exactly 4 lines each
            // here (no multiline), so this is a lossless re-chunking.
            let lines: Vec<&str> = s.lines().collect();
            assert_eq!(
                lines.len() % 4,
                0,
                "expected whole 4-line FASTQ records, got {} lines",
                lines.len()
            );
            let mut v: Vec<String> = lines.chunks(4).map(|c| c.join("\n")).collect();
            v.sort();
            v
        };

        let mut a = Vec::new();
        run_bam_to_fastq(
            recs.clone().into_iter().map(anyhow::Ok),
            &mut a,
            &mk(1),
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let mut b = Vec::new();
        run_bam_to_fastq(
            recs.into_iter().map(anyhow::Ok),
            &mut b,
            &mk(8),
            &Arc::new(Counters::default()),
        )
        .unwrap();

        assert_eq!(
            sorted_records(&a),
            sorted_records(&b),
            "t1 and t8 FASTQ must match as a multiset"
        );
    }
}
