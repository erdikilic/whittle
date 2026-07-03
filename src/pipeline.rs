use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

use noodles_sam::{self as sam, alignment::RecordBuf};
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

use crate::config::{Config, FastqTags};
use crate::io::fastq::{format_aux_field, format_mods_aux, write_segment, write_segment_tagged};
use crate::record::ReadRecord;
use crate::{filter, mods, trim};

#[derive(Debug, Default, Clone, Copy)]
pub struct Stats {
    pub input_reads: u64,
    pub output_reads: u64,
}

/// Single-threaded FASTQ pipeline: filter -> trim -> write each surviving segment.
pub fn run_fastq_seq<W: Write>(
    records: impl Iterator<Item = anyhow::Result<ReadRecord>>,
    writer: &mut W,
    cfg: &Config,
) -> anyhow::Result<Stats> {
    let mut stats = Stats::default();
    for rec in records {
        let rec = rec?;
        stats.input_reads += 1;
        if !filter::passes(&rec.seq, &rec.qual, &cfg.filter) {
            continue;
        }
        let intervals = trim::apply(rec.seq.len(), &rec.qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            write_segment(writer, &rec.name, &rec.seq[s..e], &rec.qual[s..e], total, idx)?;
            stats.output_reads += 1;
        }
    }
    Ok(stats)
}

/// Format the surviving segments of one record into an owned FASTQ byte buffer.
/// Returns the number of segments written alongside the buffer.
fn render_record(rec: &ReadRecord, cfg: &Config) -> (u64, Vec<u8>) {
    if !filter::passes(&rec.seq, &rec.qual, &cfg.filter) {
        return (0, Vec::new());
    }
    let intervals = trim::apply(rec.seq.len(), &rec.qual, &cfg.trim, cfg.filter.min_length);
    let total = intervals.len();
    let mut buf = Vec::new();
    let mut out = 0u64;
    for (idx, (s, e)) in intervals.into_iter().enumerate() {
        write_segment(&mut buf, &rec.name, &rec.seq[s..e], &rec.qual[s..e], total, idx).unwrap();
        out += 1;
    }
    (out, buf)
}

/// Threads-aware FASTQ pipeline entry point. Sequential (and output-order
/// deterministic) when `cfg.threads <= 1`; otherwise renders each record on a
/// rayon work pool and drains rendered buffers through a bounded channel on a
/// dedicated writer task. Output order is unordered for `threads > 1` since
/// buffers are written in arrival order, not input order.
///
/// Record-parse errors (`Err` items from `records`): both paths surface them
/// as an `Err` from this function. The sequential path aborts immediately via
/// `?`; the parallel path captures the first parse error in a shared slot and
/// keeps draining so producer threads never block on the bounded channel,
/// surfacing the error once the scope joins (matching the sequential path's
/// fail behavior, just at end-of-run instead of mid-stream).
pub fn run_fastq<W, I>(records: I, writer: &mut W, cfg: &Config) -> anyhow::Result<Stats>
where
    W: Write + Send,
    I: Iterator<Item = anyhow::Result<ReadRecord>> + Send,
{
    if cfg.threads <= 1 {
        return run_fastq_seq(records, writer, cfg);
    }

    let pool = rayon::ThreadPoolBuilder::new().num_threads(cfg.threads).build()?;
    let input_reads = AtomicU64::new(0);
    let output_reads = AtomicU64::new(0);
    let (tx, rx) = crossbeam_channel::bounded::<Vec<u8>>(cfg.threads * 4);

    let write_err: std::sync::Mutex<Option<std::io::Error>> = std::sync::Mutex::new(None);
    let parse_err: std::sync::Mutex<Option<anyhow::Error>> = std::sync::Mutex::new(None);

    pool.in_place_scope(|scope| {
        // Writer task drains rendered buffers in arrival order. On a write
        // error we must keep draining (not `break`): `rx` stays alive on this
        // stack frame, so if we stopped receiving, producers blocked on the
        // bounded `tx.send` in the `par_bridge().for_each` below would never
        // unblock, `in_place_scope` would never return, and the recorded
        // error would never surface — a deadlock instead of an `Err`.
        scope.spawn(|_| {
            let mut errored = false;
            for buf in rx.iter() {
                if errored {
                    continue; // keep draining so bounded-channel producers never block
                }
                if let Err(e) = writer.write_all(&buf) {
                    *write_err.lock().unwrap() = Some(e);
                    errored = true;
                }
            }
        });

        records.par_bridge().for_each(|rec| {
            let rec = match rec {
                Ok(r) => r,
                Err(e) => {
                    let mut g = parse_err.lock().unwrap();
                    if g.is_none() {
                        *g = Some(e);
                    }
                    return;
                }
            };
            input_reads.fetch_add(1, Ordering::Relaxed);
            let (out, buf) = render_record(&rec, cfg);
            if out > 0 {
                output_reads.fetch_add(out, Ordering::Relaxed);
                let _ = tx.send(buf);
            }
        });
        drop(tx);
    });

    if let Some(e) = parse_err.lock().unwrap().take() {
        return Err(e);
    }
    if let Some(e) = write_err.lock().unwrap().take() {
        return Err(e.into());
    }
    Ok(Stats {
        input_reads: input_reads.load(Ordering::Relaxed),
        output_reads: output_reads.load(Ordering::Relaxed),
    })
}

/// Build one output uBAM record for interval [start,end): slice SEQ/QUAL, rebuild
/// MM/ML/MN, suffix the name on splits. Non-mod aux tags are carried through
/// unchanged (they ride along in the cloned RecordBuf).
pub fn reconstruct_record(
    src: &RecordBuf,
    start: usize,
    end: usize,
    total: usize,
    idx: usize,
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
    if matches!(src.data().get(&Tag::BASE_MODIFICATIONS), Some(Value::String(_))) {
        let data = out.data_mut();
        match reconstruct_mods(src, &seq, start, end) {
            Some((mm_new, ml_new)) => {
                data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm_new.into()));
                data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(ml_new)));
                data.insert(
                    Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
                    Value::Int32((end - start) as i32),
                );
            }
            None => {
                data.remove(&Tag::BASE_MODIFICATIONS);
                data.remove(&Tag::BASE_MODIFICATION_PROBABILITIES);
                data.remove(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH);
            }
        }
    }

    out
}

/// Slice MM/ML to the window `[start, end)` and re-serialize. Returns `None`
/// when the source has no `MM` tag, or when no modified position survives the
/// window (caller drops MM/ML/MN in that case). Shared by the BAM→BAM and
/// BAM→FASTQ paths so they cannot drift.
pub fn reconstruct_mods(
    src: &RecordBuf,
    seq: &[u8],
    start: usize,
    end: usize,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let mm_raw = match src.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => s.to_vec(),
        _ => return None,
    };
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
        Some((mm_new, ml_new))
    }
}

/// Single-threaded uBAM pipeline: refuse aligned reads, filter, trim, reconstruct.
pub fn run_bam<W>(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<RecordBuf>>,
    writer: &mut noodles_bam::io::Writer<W>,
    cfg: &Config,
) -> anyhow::Result<Stats>
where
    W: std::io::Write,
{
    let mut stats = Stats::default();
    for rec in records {
        let rec = rec?;
        crate::io::bam::ensure_unaligned(&rec)?;
        stats.input_reads += 1;

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
            continue;
        }
        let intervals = trim::apply(seq.len(), &qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            let out = reconstruct_record(&rec, s, e, total, idx);
            writer.write_alignment_record(header, &out)?;
            stats.output_reads += 1;
        }
    }
    Ok(stats)
}

/// Assemble the TAB-prefixed aux-tag block for one interval: carried non-mod
/// tags verbatim in source order, then the reconstructed MM/ML/MN block. Empty
/// when nothing is carried (caller writes a plain FASTQ record).
fn build_fastq_tags(
    src: &RecordBuf,
    seq: &[u8],
    start: usize,
    end: usize,
    sel: &FastqTags,
) -> Vec<u8> {
    let mut tags = Vec::new();
    for (tag, value) in src.data().iter() {
        let t = <[u8; 2]>::from(tag);
        if t == *b"MM" || t == *b"ML" || t == *b"MN" {
            continue; // handled by the reconstructed block below
        }
        if sel.carries(&t) {
            tags.push(b'\t');
            tags.extend_from_slice(&format_aux_field(t, value));
        }
    }
    if sel.carries_mods()
        && let Some((mm, ml)) = reconstruct_mods(src, seq, start, end)
    {
        tags.push(b'\t');
        tags.extend_from_slice(&format_mods_aux(&mm, &ml, end - start));
    }
    tags
}

/// Single-threaded uBAM→FASTQ pipeline: refuse aligned reads, filter, trim, then
/// write each surviving segment as FASTQ with the selected aux tags in the header
/// (MM/ML/MN reconstructed; others verbatim). gz compression, when requested, is
/// handled by the parallel `gzp` writer this drains into.
pub fn run_bam_to_fastq<W>(
    records: impl Iterator<Item = anyhow::Result<RecordBuf>>,
    writer: &mut W,
    cfg: &Config,
) -> anyhow::Result<Stats>
where
    W: Write,
{
    let mut stats = Stats::default();
    for rec in records {
        let rec = rec?;
        crate::io::bam::ensure_unaligned(&rec)?;
        stats.input_reads += 1;

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
            continue;
        }
        let name = rec.name().map(|n| n.to_vec()).unwrap_or_default();
        let intervals = trim::apply(seq.len(), &qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            let tags = build_fastq_tags(&rec, &seq, s, e, &cfg.fastq_tags);
            if tags.is_empty() {
                write_segment(writer, &name, &seq[s..e], &qual[s..e], total, idx)?;
            } else {
                write_segment_tagged(writer, &name, &seq[s..e], &qual[s..e], total, idx, &tags)?;
            }
            stats.output_reads += 1;
        }
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::FilterConfig;
    use crate::qual::QualMode;
    use crate::record::ReadRecord;
    use crate::trim::{QualityOp, TrimPlan};

    fn rec(name: &str, seq: &[u8], phred: Vec<u8>) -> ReadRecord {
        ReadRecord { name: name.as_bytes().to_vec(), seq: seq.to_vec(), qual: phred }
    }

    fn base_filter() -> FilterConfig {
        FilterConfig {
            min_length: 1, max_length: usize::MAX, min_qual: 0.0, max_qual: 1000.0,
            min_gc: None, max_gc: None, qual_mode: QualMode::Mean,
        }
    }

    #[test]
    fn fixed_crop_writes_one_segment() {
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 1, tail: 1, quality: None },
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
        };
        let recs = vec![Ok(rec("r1", b"ACGT", vec![40, 40, 40, 40]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@r1\nCG\n+\nII\n");
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
    }

    #[test]
    fn split_writes_suffixed_segments() {
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 0, tail: 0, quality: Some(QualityOp::Split { cutoff: 10, window: 1 }) },
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
        };
        // good(3) bad(1) good(3): I I I # I I I  -> two segments (0,3),(4,7)
        let phred: Vec<u8> = b"III#III".iter().map(|&b| b - 33).collect();
        let recs = vec![Ok(rec("r1", b"AAATAAA", phred))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@r1_segment_1\nAAA\n+\nIII\n@r1_segment_2\nAAA\n+\nIII\n");
        assert_eq!((stats.input_reads, stats.output_reads), (1, 2));
    }

    #[test]
    fn filtered_read_produces_no_output() {
        let mut f = base_filter();
        f.min_length = 10;
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: f,
            trim: TrimPlan { head: 0, tail: 0, quality: None },
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
        };
        let recs = vec![Ok(rec("short", b"ACGT", vec![40; 4]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert!(out.is_empty());
        assert_eq!((stats.input_reads, stats.output_reads), (1, 0));
    }

    #[test]
    fn parallel_matches_sequential_as_multiset() {
        use crate::config::IoConfig;
        let mk = |threads| Config {
            io: IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 0, tail: 0, quality: Some(QualityOp::TrimQual(20)) },
            threads,
            fastq_tags: crate::config::FastqTags::All,
        };
        // Owned records (ReadRecord: Clone); wrap in Ok at iteration time so each run
        // gets a fresh Send iterator. anyhow::Error is not Clone, so we can't clone a
        // Vec<Result<..>> — clone the Vec<ReadRecord> and re-wrap instead.
        let recs: Vec<ReadRecord> = (0..500)
            .map(|i| rec(&format!("r{i}"), b"ACGTACGTAC", vec![40; 10]))
            .collect();

        let mut seq_out = Vec::new();
        run_fastq(recs.clone().into_iter().map(anyhow::Ok), &mut seq_out, &mk(1)).unwrap();

        let mut par_out = Vec::new();
        run_fastq(recs.into_iter().map(anyhow::Ok), &mut par_out, &mk(4)).unwrap();

        let sort_records = |bytes: &[u8]| {
            let mut v: Vec<Vec<u8>> = bytes
                .split(|&b| b == b'@')
                .filter(|s| !s.is_empty())
                .map(|s| s.to_vec())
                .collect();
            v.sort();
            v
        };
        assert_eq!(sort_records(&seq_out), sort_records(&par_out));
    }

    #[test]
    fn parallel_surfaces_write_error_without_deadlock() {
        use crate::config::IoConfig;
        use std::io::{self, Write};

        struct FailAfter { limit: usize, written: usize }
        impl Write for FailAfter {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                if self.written >= self.limit {
                    return Err(io::Error::new(io::ErrorKind::BrokenPipe, "boom"));
                }
                self.written += buf.len();
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> { Ok(()) }
        }

        let cfg = Config {
            io: IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 0, tail: 0, quality: None },
            threads: 4,
            fastq_tags: crate::config::FastqTags::All,
        };
        // Far more records than the bounded channel capacity (threads*4), so a
        // pre-fix build would deadlock instead of returning.
        let recs: Vec<ReadRecord> = (0..2000)
            .map(|i| rec(&format!("r{i}"), b"ACGTACGTAC", vec![40; 10]))
            .collect();
        let mut w = FailAfter { limit: 100, written: 0 };
        let res = run_fastq(recs.into_iter().map(anyhow::Ok), &mut w, &cfg);
        assert!(res.is_err(), "write error must surface as Err, and must not hang");
    }

    #[test]
    fn parallel_surfaces_parse_error_instead_of_dropping_it() {
        use crate::config::IoConfig;

        let cfg = Config {
            io: IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 0, tail: 0, quality: None },
            threads: 4,
            fastq_tags: crate::config::FastqTags::All,
        };
        let good: Vec<anyhow::Result<ReadRecord>> = (0..5)
            .map(|i| anyhow::Ok(rec(&format!("r{i}"), b"ACGTACGTAC", vec![40; 10])))
            .collect();
        let recs = good.into_iter().chain(std::iter::once(Err(anyhow::anyhow!("bad record"))));

        let mut out = Vec::new();
        let res = run_fastq(recs, &mut out, &cfg);
        assert!(res.is_err(), "a malformed record must not be silently dropped on the parallel path");
    }
}

#[cfg(test)]
mod bam_tests {
    use super::*;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    fn ubam_with_mods(seq: &[u8], quals: Vec<u8>, mm: &[u8], ml: Vec<u8>) -> RecordBuf {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = seq.to_vec().into();
        *rec.quality_scores_mut() = quals.into();
        let data = rec.data_mut();
        data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.to_vec().into()));
        data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(ml)));
        rec
    }

    #[test]
    fn slices_seq_qual_and_rebuilds_tags() {
        // seq = C C A C ; C+m modified at C occ 0 and 2 -> pos 0 and 3; ML [10,20].
        let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0,1;", vec![10, 20]);
        // keep window [2,4): seq "AC", the modified C at pos3 survives (occ within window idx0).
        let out = reconstruct_record(&src, 2, 4, 1, 0);

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
        let mut buf: Vec<u8> = Vec::new();
        let mut writer = noodles_bam::io::Writer::new(&mut buf);

        let cfg = Config {
            io: IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: FilterConfig {
                min_length: 1, max_length: usize::MAX, min_qual: 0.0, max_qual: 1000.0,
                min_gc: None, max_gc: None, qual_mode: QualMode::Mean,
            },
            trim: TrimPlan { head: 0, tail: 0, quality: None },
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
        };

        let result = run_bam(&header, [Ok(rec)].into_iter(), &mut writer, &cfg);
        assert!(result.is_err(), "SEQ/QUAL length mismatch must error, not panic");
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
        data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(vec![1, 2, 3])));
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));

        let out = reconstruct_record(&src, 1, 4, 1, 0);

        match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::Int32(5)) => {}
            other => panic!("MM must be left untouched for a non-string value, got {other:?}"),
        }
        match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) if v == &[1u8, 2, 3] => {}
            other => panic!("ML must be left untouched, got {other:?}"),
        }
        match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(4)) => {}
            other => panic!("MN must be left untouched, got {other:?}"),
        }
    }

    #[test]
    fn split_suffixes_name_and_drops_empty_mods() {
        let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0;", vec![10]); // mod at pos0
        // segment [2,4) has no surviving C mod -> MM/ML removed entirely.
        let out = reconstruct_record(&src, 2, 4, 2, 1);
        // `.as_ref()` is ambiguous on `&BStr` (impls `AsRef<[u8]>` and `AsRef<BStr>`);
        // disambiguate with a turbofish (noodles-sam 0.85 / bstr 1.x adjustment).
        assert_eq!(AsRef::<[u8]>::as_ref(out.name().unwrap()), b"r1_segment_2");
        assert!(out.data().get(&Tag::BASE_MODIFICATIONS).is_none());
        assert!(out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES).is_none());
        assert!(out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH).is_none());
    }

    use crate::config::{FastqTags, IoConfig};
    use crate::filter::FilterConfig;
    use crate::qual::QualMode;
    use crate::trim::{QualityOp, TrimPlan};

    fn cfg_bam2fq(quality: Option<QualityOp>, head: usize, tags: FastqTags) -> Config {
        Config {
            io: IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: FilterConfig {
                min_length: 1, max_length: usize::MAX, min_qual: 0.0, max_qual: 1000.0,
                min_gc: None, max_gc: None, qual_mode: QualMode::Mean,
            },
            trim: TrimPlan { head, tail: 0, quality },
            threads: 1,
            fastq_tags: tags,
        }
    }

    // "CCACCCAC" C at seq idx 0,1,3,4,5,7; MM "C+m,0,1,0" -> occ 0,2,3 -> abs 0,3,4,
    // ML [10,20,30]. head-crop 2 -> window [2,8): keeps abs 3,4 renumbered -> "C+m,0,0;" ML [20,30] MN 6.
    fn read2_with_mods_and_rg() -> RecordBuf {
        let mut rec = ubam_with_mods(b"CCACCCAC", vec![35; 8], b"C+m,0,1,0;", vec![10, 20, 30]);
        rec.data_mut().insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
        rec.data_mut().insert(
            Tag::READ_GROUP,
            Value::String(b"grp1".as_slice().into()),
        );
        rec
    }

    #[test]
    fn bam2fq_all_carries_rg_and_reconstructed_mods() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::All);
        let mut out = Vec::new();
        let stats = run_bam_to_fastq([Ok(read2_with_mods_and_rg())].into_iter(), &mut out, &cfg).unwrap();
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
        let s = String::from_utf8(out).unwrap();
        // header carries RG verbatim + reconstructed mod block; seq head-cropped by 2.
        assert!(s.starts_with("@r1\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6\n"), "got: {s:?}");
        assert!(s.contains("\nACCCAC\n+\n"), "cropped seq wrong: {s:?}");
    }

    #[test]
    fn bam2fq_only_mm_ml_drops_rg() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::parse("MM,ML").unwrap());
        let mut out = Vec::new();
        run_bam_to_fastq([Ok(read2_with_mods_and_rg())].into_iter(), &mut out, &cfg).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(!s.contains("RG:Z"), "RG must be dropped: {s:?}");
        assert!(s.contains("MM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"), "mods missing: {s:?}");
    }

    #[test]
    fn bam2fq_none_is_plain_fastq() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::None);
        let mut out = Vec::new();
        run_bam_to_fastq([Ok(read2_with_mods_and_rg())].into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@r1\nACCCAC\n+\nDDDDDD\n"); // 35+33 = 'D'
    }

    #[test]
    fn bam2fq_split_suffixes_and_segments_mods() {
        // split at the low-qual base; each segment gets its own reconstructed mods.
        let cfg = cfg_bam2fq(Some(QualityOp::Split { cutoff: 20, window: 1 }), 0, FastqTags::All);
        // seq CCAC, C+m at occ 0 and 2 -> abs 0,3; qual: good good BAD good so split [0,2),[3,4)
        let mut rec = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0,1;", vec![100, 200]);
        rec.data_mut().insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));
        let mut out = Vec::new();
        let stats = run_bam_to_fastq([Ok(rec)].into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(stats.output_reads, 2);
        let s = String::from_utf8(out).unwrap();
        // segment 1 = [0,2) "CC" keeps abs-0 mod; segment 2 = [3,4) "C" keeps abs-3 mod.
        assert!(s.contains("@r1_segment_1\tMM:Z:C+m,0;\tML:B:C,100\tMN:i:2"), "seg1: {s:?}");
        assert!(s.contains("@r1_segment_2\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:1"), "seg2: {s:?}");
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
        run_bam_to_fastq([Ok(rec)].into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@plain\nACGT\n+\nIIII\n");
    }
}
