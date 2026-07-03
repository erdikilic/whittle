use std::io::Write;

use noodles_sam::{self as sam, alignment::RecordBuf};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

use crate::config::{Config, FastqTags};
use crate::io::fastq::{format_aux_field, format_mods_aux, write_segment, write_segment_tagged};
use crate::{filter, mods, trim};

use super::Stats;

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
pub fn run_bam(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<RecordBuf>>,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
) -> anyhow::Result<Stats> {
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
            sink.write_record(header, &out)?;
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
        let dir = tempfile::tempdir().unwrap();
        let mut sink = crate::io::bam::writer(Some(&dir.path().join("o.bam")), &header, 1).unwrap();

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

        let result = run_bam(&header, [Ok(rec)].into_iter(), &mut sink, &cfg);
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
