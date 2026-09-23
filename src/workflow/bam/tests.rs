use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

use super::*;

#[test]
fn signal_direction_uses_the_record_group_or_a_shared_header_direction() {
    let header: sam::Header =
        "@RG\tID:dna\tDS:basecall_model=dna_r10.4.1\n@RG\tID:rna\tDS:basecall_model=rna004\n"
            .parse()
            .unwrap();
    let mut rec = RecordBuf::default();
    assert_eq!(signal_reversed(&header, &rec), None);
    for (id, expected) in [
        (b"dna".as_slice(), Some(false)),
        (b"rna".as_slice(), Some(true)),
        (b"absent".as_slice(), None),
    ] {
        rec.data_mut()
            .insert(Tag::READ_GROUP, Value::String(id.into()));
        assert_eq!(signal_reversed(&header, &rec), expected);
    }
    rec.data_mut().remove(&Tag::READ_GROUP);
    let header: sam::Header =
        "@RG\tID:a\tDS:basecall_model=rna004\n@RG\tID:b\tDS:runid=run;basecall_model=rna002\n"
            .parse()
            .unwrap();
    assert_eq!(signal_reversed(&header, &rec), Some(true));
    assert_eq!(signal_reversed(&sam::Header::default(), &rec), None);
}

/// Builds one output record for interval `[start, end)`, segment `idx` of
/// `total`, with the modification block classified from `src` and no tag
/// removal.
pub(super) fn reconstruct_record(
    src: &RecordBuf,
    start: usize,
    end: usize,
    total: usize,
    idx: usize,
    update_moves: bool,
) -> RecordBuf {
    let seq = src.sequence().as_ref();
    let mod_block = inspect_mod_block(src, seq.len());
    let window = Window {
        start,
        end,
        idx,
        total,
    };
    reconstruct_window_record(
        src,
        window,
        mod_block,
        None,
        update_moves
            .then(|| MoveIndex::new(src, false))
            .flatten()
            .as_ref(),
        &TagRemoval::default(),
    )
    .unwrap_or_else(|| src.clone())
}

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
    // Seq CCAC; `C+m` modified at C occurrences 0 and 2, positions 0 and 3;
    // ML [10, 20].
    let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0,1;", vec![10, 20]);
    // Window [2,4) keeps "AC"; the modified C at position 3 survives as
    // window occurrence 0.
    let out = reconstruct_record(&src, 2, 4, 1, 0, false);

    assert_eq!(out.sequence().as_ref(), b"AC");
    assert_eq!(out.quality_scores().as_ref(), &[32, 33]);

    let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => s.to_vec(),
        _ => panic!("No MM"),
    };
    assert_eq!(mm, b"C+m,0;");
    let ml = match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        _ => panic!("No ML"),
    };
    assert_eq!(ml, vec![20]);
    // MN updated to the output length.
    let mn = match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
        Some(Value::Int32(n)) => *n,
        _ => panic!("No MN"),
    };
    assert_eq!(mn, 2);
}

/// A BAM record with unequal SEQ and QUAL lengths returns an error.
#[test]
fn qual_seq_length_mismatch_errors_without_panicking() {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = b"ACGT".to_vec().into();
    // `quality_scores` is left at its default (empty), so SEQ and QUAL
    // lengths differ.

    let header = sam::Header::default();
    let dir = tempfile::tempdir().unwrap();
    let mut sink =
        crate::io::bam::writer(Some(&dir.path().join("o.bam")), &header, false, 6).unwrap();

    let cfg = Config {
        quiet: true,
        ..Config::default()
    };

    let result = run_bam(
        &header,
        [Ok(raw_record(&rec))].into_iter(),
        &mut sink,
        &cfg,
        &Arc::new(Counters::default()),
    );
    assert!(
        result.is_err(),
        "SEQ/QUAL length mismatch must error, not panic"
    );
}

/// A non-string MM value is outside the supported schema and remains untouched.
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

/// An empty assessment group remains present alongside populated groups.
#[test]
fn reconstruct_record_preserves_originally_empty_mm_group() {
    // Seq CACA: C at 0, 2; A at 1, 3. `A+a` modifies A occurrence 0
    // (position 1). `C+m` is empty.
    let src = ubam_with_mods(b"CACA", vec![30, 31, 32, 33], b"A+a,0;C+m;", vec![7]);
    // Missing MN requires reconstruction even for the complete sequence.
    let out = reconstruct_record(&src, 0, 4, 1, 0, false);
    let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => s.to_vec(),
        other => panic!("MM must survive, got {other:?}"),
    };
    assert_eq!(
        mm, b"A+a,0;C+m;",
        "The empty C+m assessment group must be preserved"
    );
}

/// A segment that loses every listed position keeps the group with no
/// positions: `C+m;` still declares the segment's C's canonical, which the
/// group's absence would turn into "no call".
#[test]
fn split_suffixes_name_and_keeps_empty_mod_group() {
    let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0;", vec![10]); // mod at position 0
    // Segment [2,4) has no surviving C modification.
    let out = reconstruct_record(&src, 2, 4, 2, 1, false);
    // `.as_ref()` is ambiguous on `&BStr` (it implements both `AsRef<[u8]>`
    // and `AsRef<BStr>`); the turbofish selects the byte view.
    assert_eq!(AsRef::<[u8]>::as_ref(out.name().unwrap()), b"r1_segment_2");
    match out.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => assert_eq!(s.to_vec(), b"C+m;"),
        other => panic!("Empty group must be kept, got {other:?}"),
    }
    match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
        Some(Value::Array(Array::UInt8(v))) => assert!(v.is_empty(), "ML must be empty"),
        other => panic!("ML must be an empty B:C array, got {other:?}"),
    }
    assert_eq!(
        out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH),
        Some(&Value::Int32(2))
    );
}

/// The FASTQ header spells an empty group as `MM:Z:C+m;` with a zero-length
/// `ML:B:C` array.
#[test]
fn bam2fq_keeps_empty_mod_group_after_a_crop() {
    let rec = ubam_with_mods(b"CCAC", vec![40; 4], b"C+m,0;", vec![10]);
    let cfg = cfg_bam2fq(None, 2, FastqTags::All);
    let mut out = Vec::new();
    run_bam_to_fastq(
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(
        s.starts_with("@r1\tMM:Z:C+m;\tML:B:C\tMN:i:2\n"),
        "Got: {s:?}"
    );
}

use crate::config::FastqTags;

use crate::trim::{QualityOp, TrimPlan};

pub(super) fn cfg_bam2fq(quality: Option<QualityOp>, head: usize, tags: FastqTags) -> Config {
    Config {
        trim: TrimPlan {
            head,
            tail: 0,
            quality,
        },
        fastq_tags: tags,
        quiet: true,
        ..Config::default()
    }
}

/// `CCACCCAC` has C at 0, 1, 3, 4, 5, 7; `C+m,0,1,0` marks occurrences 0, 2,
/// 3 (positions 0, 3, 4) with ML [10, 20, 30]. A head crop of 2 (window
/// [2,8)) keeps positions 3 and 4, renumbered to `C+m,0,0;` with ML [20, 30]
/// and MN 6.
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
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
    let s = String::from_utf8(out).unwrap();
    // The header carries RG unchanged and the reconstructed mod block; the
    // sequence is head-cropped by 2.
    assert!(
        s.starts_with("@r1\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6\n"),
        "Got: {s:?}"
    );
    assert!(s.contains("\nACCCAC\n+\n"), "Cropped sequence wrong: {s:?}");
}

#[test]
fn bam2fq_only_mm_ml_drops_rg() {
    let cfg = cfg_bam2fq(None, 2, FastqTags::parse("MM,ML").unwrap());
    let mut out = Vec::new();
    run_bam_to_fastq(
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(!s.contains("RG:Z"), "RG must be dropped: {s:?}");
    assert!(
        s.contains("MM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"),
        "Mods missing: {s:?}"
    );
}

#[test]
fn bam2fq_none_is_plain_fastq() {
    let cfg = cfg_bam2fq(None, 2, FastqTags::None);
    let mut out = Vec::new();
    run_bam_to_fastq(
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(out, b"@r1\nACCCAC\n+\nDDDDDD\n"); // 35+33 = 'D'
}

/// A split at the low-quality base gives each segment its own reconstructed
/// mods.
#[test]
fn bam2fq_split_suffixes_and_segments_mods() {
    let cfg = cfg_bam2fq(
        Some(QualityOp::Split {
            cutoff: 20,
            window: 1,
        }),
        0,
        FastqTags::All,
    );
    // Seq CCAC, `C+m` at occurrences 0 and 2 (positions 0 and 3); quality
    // good, good, bad, good, so the split is [0,2), [3,4).
    let mut rec = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0,1;", vec![100, 200]);
    rec.data_mut()
        .insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));
    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(stats.output_reads, 2);
    let s = String::from_utf8(out).unwrap();
    // Segment 1 = [0,2) "CC" keeps the position-0 mod; segment 2 = [3,4) "C"
    // keeps the position-3 mod.
    assert!(
        s.contains("@r1_segment_1\tMM:Z:C+m,0;\tML:B:C,100\tMN:i:2"),
        "Segment 1: {s:?}"
    );
    assert!(
        s.contains("@r1_segment_2\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:1"),
        "Segment 2: {s:?}"
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
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(out, b"@plain\nACGT\n+\nIIII\n");
}

/// MM without optional ML remains MM-only after reconstruction.
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
    // No ML and no MN.

    let out = reconstruct_record(&src, 0, 4, 1, 0, false);

    // MM is retained: both modified Cs are in the window, so `C+m,0,1;`.
    let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => s.to_vec(),
        other => panic!("Expected MM retained, got {other:?}"),
    };
    assert_eq!(mm, b"C+m,0,1;");
    // ML must be absent, never an empty array.
    assert!(
        out.data()
            .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
            .is_none(),
        "MM-only source must not gain an ML tag"
    );
    // MN set to the window length.
    match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
        Some(Value::Int32(4)) => {},
        other => panic!("Expected MN=4, got {other:?}"),
    }
}

/// Companion for the BAM-to-FASTQ path: an MM-only source must emit a FASTQ
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
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(
        s.contains("MM:Z:C+m,0,1;\tMN:i:4"),
        "Expected MM+MN, got: {s:?}"
    );
    assert!(
        !s.contains("ML:B"),
        "MM-only record must not emit an ML field: {s:?}"
    );
}

/// A record whose modification block cannot be placed on its sequence: the
/// variants `malformed_mod_blocks` feeds through both output paths.
fn malformed_mod_record(variant: &str) -> RecordBuf {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = b"CCCA".to_vec().into();
    *rec.quality_scores_mut() = vec![40; 4].into();
    let d = rec.data_mut();
    let (mm, ml, mn) = match variant {
        // MM declares 3 positions, ML has 1 byte.
        "ml_short" => (&b"C+m,0,0,0;"[..], Value::Array(Array::UInt8(vec![5])), 4),
        // MM stops parsing at the `x`.
        "mm_garbled" => (
            &b"C+m,0,0x,0;"[..],
            Value::Array(Array::UInt8(vec![5, 6, 7])),
            4,
        ),
        // ML at a subtype other than B:C.
        "ml_subtype" => (
            &b"C+m,0,0,0;"[..],
            Value::Array(Array::Int8(vec![5, 6, 7])),
            4,
        ),
        // MN disagrees with the sequence length.
        "mn_mismatch" => (
            &b"C+m,0,0,0;"[..],
            Value::Array(Array::UInt8(vec![5, 6, 7])),
            40,
        ),
        other => panic!("Unknown variant {other}"),
    };
    d.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.to_vec().into()));
    d.insert(Tag::BASE_MODIFICATION_PROBABILITIES, ml);
    d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(mn));
    d.insert(Tag::READ_GROUP, Value::String(b"grp".as_slice().into()));
    rec
}

const MALFORMED_MOD_VARIANTS: [&str; 4] = ["ml_short", "mm_garbled", "ml_subtype", "mn_mismatch"];

/// Each defect classifies as `Malformed`, on the full window and on a crop.
#[test]
fn malformed_mod_blocks_are_classified() {
    for variant in MALFORMED_MOD_VARIANTS {
        let rec = malformed_mod_record(variant);
        assert_eq!(inspect_mod_block(&rec, 4), ModBlock::Malformed, "{variant}");
    }
    let ok = ubam_with_mods(b"CCCA", vec![40; 4], b"C+m,0,0,0;", vec![5, 6, 7]);
    assert_eq!(inspect_mod_block(&ok, 4), ModBlock::MissingMn);
}

/// The BAM output drops the whole block, keeps every other tag, and the run
/// counts the read once, whether or not the read is trimmed or split.
#[test]
fn malformed_mod_block_is_removed_and_counted_on_bam_output() {
    for variant in MALFORMED_MOD_VARIANTS {
        for head in [0, 1] {
            let mut cfg = cfg_bam2fq(None, head, FastqTags::All);
            cfg.trim.quality = Some(QualityOp::Split {
                cutoff: 20,
                window: 1,
            });
            let mut rec = malformed_mod_record(variant);
            // A low-quality base in the middle splits the read in two.
            *rec.quality_scores_mut() = vec![40, 40, 1, 40].into();
            let header = sam::Header::default();
            let dir = tempfile::tempdir().unwrap();
            let path = dir.path().join("o.bam");
            let mut sink = crate::io::bam::writer(Some(&path), &header, false, 6).unwrap();
            let counters = Arc::new(Counters::default());
            let stats = run_bam(
                &header,
                [Ok(raw_record(&rec))].into_iter(),
                &mut sink,
                &cfg,
                &counters,
            )
            .unwrap();
            sink.finish().unwrap();
            assert_eq!(
                stats.malformed_mod_reads, 1,
                "{variant} head={head}: counted once per read"
            );

            let bytes = std::fs::read(&path).unwrap();
            let mut reader = noodles_bam::io::Reader::new(bytes.as_slice());
            let h = reader.read_header().unwrap();
            let mut buf = RecordBuf::default();
            let mut n = 0;
            while reader.read_record_buf(&h, &mut buf).unwrap() != 0 {
                n += 1;
                for t in MOD_TAGS {
                    assert!(
                        buf.data().get(&t).is_none(),
                        "{variant} head={head}: {t:?} must be removed"
                    );
                }
                assert!(
                    buf.data().get(&Tag::READ_GROUP).is_some(),
                    "{variant}: other tags are copied"
                );
            }
            assert_eq!(n, stats.output_reads as usize);
            assert!(n >= 1, "{variant} head={head}: segments written");
        }
    }
}

/// The FASTQ header carries no MM/ML/MN for a malformed block, and the run
/// counts the read.
#[test]
fn malformed_mod_block_is_omitted_and_counted_on_fastq_output() {
    for variant in MALFORMED_MOD_VARIANTS {
        for head in [0, 2] {
            let cfg = cfg_bam2fq(None, head, FastqTags::All);
            let mut out = Vec::new();
            let stats = run_bam_to_fastq(
                [Ok(raw_record(&malformed_mod_record(variant)))].into_iter(),
                &mut out,
                &cfg,
                &Arc::new(Counters::default()),
            )
            .unwrap();
            assert_eq!(stats.malformed_mod_reads, 1, "{variant} head={head}");
            let s = String::from_utf8(out).unwrap();
            assert!(
                s.starts_with("@r1\tRG:Z:grp\n"),
                "{variant} head={head}: mod block must be omitted, got {s:?}"
            );
        }
    }
}

/// A well-formed block is not counted, and an absent `MN` is added rather
/// than treated as a defect.
#[test]
fn well_formed_mod_block_is_not_counted_and_gains_mn() {
    let rec = ubam_with_mods(b"CCCA", vec![40; 4], b"C+m,0,0,0;", vec![5, 6, 7]);
    let cfg = cfg_bam2fq(None, 0, FastqTags::All);
    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(stats.malformed_mod_reads, 0);
    let s = String::from_utf8(out).unwrap();
    assert!(s.contains("MM:Z:C+m,0,0,0;\tML:B:C,5,6,7\tMN:i:4"), "{s:?}");
    let out = reconstruct_record(&rec, 0, 4, 1, 0, false);
    assert_eq!(
        out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH),
        Some(&Value::Int32(4)),
        "MN is added when absent"
    );
}

/// PacBio per-base kinetics (`ip`/`pw`, length equal to the read length) are
/// sliced with the sequence; the ONT `mv` (signal-space) tag is dropped on
/// trim; the per-read RG is copied unchanged.
#[test]
fn reconstruct_record_slices_kinetics_and_drops_mv() {
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

    // Window [2,5) is "GTA" (head crop 2, tail crop 1).
    let out = reconstruct_record(&src, 2, 5, 1, 0, false);
    assert_eq!(out.sequence().as_ref(), b"GTA");
    match out.data().get(&Tag::new(b'i', b'p')) {
        Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[12, 13, 14]),
        other => panic!("Expected ip sliced to [2,5): {other:?}"),
    }
    match out.data().get(&Tag::new(b'p', b'w')) {
        Some(Value::Array(Array::UInt16(v))) => assert_eq!(v, &[22, 23, 24]),
        other => panic!("Expected pw sliced to [2,5): {other:?}"),
    }
    assert!(
        out.data().get(&Tag::new(b'm', b'v')).is_none(),
        "The mv tag must be dropped on trim"
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
    // An unknown B array whose length equals the read length is sliced
    // structurally.
    src.data_mut().insert(
        Tag::new(b'z', b'z'),
        Value::Array(Array::Int32(vec![1, 2, 3, 4])),
    );
    // A B array whose length differs from the read length is not per-base
    // and is left alone.
    src.data_mut()
        .insert(Tag::new(b'x', b'y'), Value::Array(Array::UInt8(vec![9, 9])));

    let out = reconstruct_record(&src, 1, 3, 1, 0, false); // window [1,3)
    match out.data().get(&Tag::new(b'z', b'z')) {
        Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 3]),
        other => panic!("Expected zz sliced: {other:?}"),
    }
    match out.data().get(&Tag::new(b'x', b'y')) {
        Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[9, 9]),
        other => panic!("Expected xy untouched: {other:?}"),
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

    // Full window [0,4): nothing is trimmed, so everything is preserved.
    let out = reconstruct_record(&src, 0, 4, 1, 0, false);
    match out.data().get(&Tag::new(b'i', b'p')) {
        Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[1, 2, 3, 4]),
        other => panic!("Expected ip unchanged: {other:?}"),
    }
    assert!(
        out.data().get(&Tag::new(b'm', b'v')).is_some(),
        "The mv tag is kept when untrimmed"
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

    let cfg = cfg_bam2fq(None, 2, FastqTags::All); // head crop 2, window [2,6)
    let mut out = Vec::new();
    run_bam_to_fastq(
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(
        s.contains("ip:B:C,12,13,14,15"),
        "Kinetics not sliced in the FASTQ header: {s:?}"
    );
}

#[test]
fn malformed_perbase_tag_detected_and_left_untouched() {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(b"r1".into());
    *src.sequence_mut() = b"ACGT".to_vec().into();
    *src.quality_scores_mut() = vec![40; 4].into();
    // `ip` length 3 differs from read length 4: a malformed known per-base tag.
    src.data_mut().insert(
        Tag::new(b'i', b'p'),
        Value::Array(Array::UInt8(vec![1, 2, 3])),
    );

    assert!(has_malformed_perbase_tag(&src, 4));
    // It cannot be sliced, so it is left as is.
    let out = reconstruct_record(&src, 1, 3, 1, 0, false);
    match out.data().get(&Tag::new(b'i', b'p')) {
        Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[1, 2, 3]),
        other => panic!("Expected the malformed ip left as is: {other:?}"),
    }
}

/// A synthetic move table: stride 2, 6 ones (one per base) at block indexes
/// 0, 1, 3, 4, 6, 7, 8 blocks in total. Shared by the `--update-moves` tests.
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
    // Head crop 2 gives window [2,6): block_first = ones[2] = 3,
    // block_second = 8.
    let out = reconstruct_record(&ubam_with_moves(), 2, 6, 1, 0, true);
    assert_eq!(out.sequence().as_ref(), b"GTAC");
    assert_eq!(
        AsRef::<[u8]>::as_ref(out.name().unwrap()),
        b"r1",
        "A crop keeps the read name"
    );
    // `mv` = [stride] + moves[3..8] = [2] + [1,1,0,1,1].
    match out.data().get(&Tag::new(b'm', b'v')) {
        Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1, 1]),
        other => panic!("Unexpected mv: {other:?}"),
    }
    // `ts` becomes 10 + 3*2 = 16; `ns` = ts + span = 16 + (8-3)*2 = 26 (a
    // head-only crop leaves `ns` unchanged).
    match out.data().get(&Tag::new(b't', b's')) {
        Some(Value::Int32(16)) => {},
        other => panic!("Unexpected ts: {other:?}"),
    }
    match out.data().get(&Tag::new(b'n', b's')) {
        Some(Value::Int32(26)) => {},
        other => panic!("Unexpected ns: {other:?}"),
    }
    assert!(out.data().get(&Tag::new(b's', b'p')).is_none());
    assert!(out.data().get(&Tag::new(b'p', b'i')).is_none());
    // A crop keeps the read identity, so `st`/`du` stay.
    assert!(
        out.data().get(&Tag::new(b's', b't')).is_some(),
        "The st tag is kept on a crop"
    );
    assert!(
        out.data().get(&Tag::new(b'd', b'u')).is_some(),
        "The du tag is kept on a crop"
    );
}

#[test]
fn update_moves_large_signal_offsets_use_uint32_not_wrapped_i32() {
    let mut src = ubam_with_moves();
    src.data_mut()
        .insert(Tag::new(b't', b's'), Value::UInt32(2_147_483_645));
    src.data_mut()
        .insert(Tag::new(b'n', b's'), Value::UInt32(2_147_483_661));

    // Head crop 2 gives block_first = 3 with stride 2, so `ts` becomes
    // 2_147_483_651 (above `i32::MAX`) and `ns` becomes 2_147_483_661.
    let out = reconstruct_record(&src, 2, 6, 1, 0, true);

    match out.data().get(&Tag::new(b't', b's')) {
        Some(Value::UInt32(2_147_483_651)) => {},
        other => panic!("Large ts must stay positive as UInt32, got {other:?}"),
    }
    match out.data().get(&Tag::new(b'n', b's')) {
        Some(Value::UInt32(2_147_483_661)) => {},
        other => panic!("Large ns must stay positive as UInt32, got {other:?}"),
    }
}

#[test]
fn update_moves_tail_crop_shrinks_ns() {
    // Tail crop 2 gives window [0,4): block_first = ones[0] = 0,
    // block_second = ones[4] = 6.
    let out = reconstruct_record(&ubam_with_moves(), 0, 4, 1, 0, true);
    // `mv` = [stride] + moves[0..6] = [2] + [1,1,0,1,1,0].
    match out.data().get(&Tag::new(b'm', b'v')) {
        Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1, 1, 0]),
        other => panic!("Unexpected mv: {other:?}"),
    }
    // `ts` is unchanged (no head trim): 10; `ns` = 10 + (6-0)*2 = 22, below 26.
    match out.data().get(&Tag::new(b't', b's')) {
        Some(Value::Int32(10)) => {},
        other => panic!("Unexpected ts: {other:?}"),
    }
    match out.data().get(&Tag::new(b'n', b's')) {
        Some(Value::Int32(22)) => {},
        other => {
            panic!("The ns tag must shrink on a tail crop (dorado ns = trim + span): {other:?}")
        },
    }
}

#[test]
fn update_moves_split_emits_subread_tags() {
    // Split into [0,3) and [3,6): each is a dorado-style subread.
    let s1 = reconstruct_record(&ubam_with_moves(), 0, 3, 2, 0, true);
    assert_eq!(AsRef::<[u8]>::as_ref(s1.name().unwrap()), b"r1_segment_1");
    // `mv` = [2] + moves[ones[0]=0 .. ones[3]=4] = [2] + [1,1,0,1].
    match s1.data().get(&Tag::new(b'm', b'v')) {
        Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1]),
        other => panic!("Unexpected s1 mv: {other:?}"),
    }
    match s1.data().get(&Tag::new(b't', b's')) {
        Some(Value::Int32(0)) => {},
        o => panic!("Segment 1 ts should be 0: {o:?}"),
    }
    match s1.data().get(&Tag::new(b'n', b's')) {
        Some(Value::Int32(8)) => {}, // (block 4-0)*stride 2
        o => panic!("Unexpected s1 ns: {o:?}"),
    }
    // `sp` counts from the parent's POD5 signal start: ts0 10 + block 0 *
    // stride 2 (dorado adds `num_trimmed_samples` to the split point).
    match s1.data().get(&Tag::new(b's', b'p')) {
        Some(Value::Int32(10)) => {},
        o => panic!("Unexpected s1 sp: {o:?}"),
    }
    match s1.data().get(&Tag::new(b'p', b'i')) {
        Some(Value::String(s)) => assert_eq!(s.to_vec(), b"r1"),
        o => panic!("Unexpected s1 pi: {o:?}"),
    }
    // Dorado marks split products with read number -1.
    match s1.data().get(&Tag::new(b'r', b'n')) {
        Some(Value::Int32(-1)) => {},
        o => panic!("Segment 1 rn should be -1: {o:?}"),
    }
    // `st`/`du` describe the parent read; with the signal window known
    // they are recomputed for the subread (`update_moves_split_recomputes_st_and_du`).
    assert!(
        matches!(s1.data().get(&Tag::new(b's', b't')), Some(Value::String(s)) if s.to_vec() != b"2024-06-21T10:00:00Z"),
        "The st tag is recomputed on a split"
    );
    assert!(
        matches!(s1.data().get(&Tag::new(b'd', b'u')), Some(Value::Float(d)) if *d < 5.0),
        "The du tag is recomputed on a split"
    );

    let s2 = reconstruct_record(&ubam_with_moves(), 3, 6, 2, 1, true);
    assert_eq!(AsRef::<[u8]>::as_ref(s2.name().unwrap()), b"r1_segment_2");
    // `mv` = [2] + moves[ones[3]=4 .. 8] = [2] + [1,0,1,1].
    match s2.data().get(&Tag::new(b'm', b'v')) {
        Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 0, 1, 1]),
        other => panic!("Unexpected s2 mv: {other:?}"),
    }
    match s2.data().get(&Tag::new(b'n', b's')) {
        Some(Value::Int32(8)) => {}, // (8-4)*2
        o => panic!("Unexpected s2 ns: {o:?}"),
    }
    match s2.data().get(&Tag::new(b's', b'p')) {
        Some(Value::Int32(18)) => {}, // ts0 10 + block_first 4 * stride 2
        o => panic!("Unexpected s2 sp: {o:?}"),
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

    // `update_moves` off and trimmed: mv/ts/ns/sp/pi are all removed.
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

    // Head crop 2 gives window [2,6), which keeps only the Q40 bases.
    let out = reconstruct_record(&src, 2, 6, 1, 0, false);

    // The poly-A and barcode coordinate tags cannot be reconstructed and are
    // dropped.
    for t in [b"pa", b"pt", b"bi"] {
        assert!(
            out.data().get(&Tag::new(t[0], t[1])).is_none(),
            "{} must be dropped on trim",
            std::str::from_utf8(t).unwrap()
        );
    }
    // `qs` is recomputed from the trimmed (all-Q40) quality, not left at 20.
    match out.data().get(&Tag::new(b'q', b's')) {
        Some(Value::Float(q)) => {
            let expected = crate::qual::mean_prob_q(&[40, 40, 40, 40]) as f32;
            assert!(
                (q - expected).abs() < 1e-4,
                "Recomputed qs: got {q}, want {expected}"
            );
        },
        other => panic!("Unexpected qs: {other:?}"),
    }
    // Per-read metadata (RG) is untouched.
    assert!(out.data().get(&Tag::new(b'R', b'G')).is_some());
}

/// `ubam_with_moves` with head crop 2 spans the original-signal window
/// [ts0+3*2, ts0+8*2] = [16, 26]; a split segment [3,6) spans
/// [ts0+4*2, 26] = [18, 26]. The poly-A tags survive a crop that keeps the
/// tail.
#[test]
fn update_moves_crop_keeps_polya_when_tail_survives() {
    let mut src = ubam_with_moves();
    // Anchor and boundaries all inside [16,26]; the split range is a sentinel.
    src.data_mut().insert(
        Tag::new(b'p', b'a'),
        Value::Array(Array::Int32(vec![20, 18, 24, -1, -1])),
    );
    src.data_mut()
        .insert(Tag::new(b'p', b't'), Value::Int32(30));

    let out = reconstruct_record(&src, 2, 6, 1, 0, true); // head crop 2
    // A crop keeps the read identity and POD5 signal, so the absolute `pa`
    // stays valid.
    match out.data().get(&Tag::new(b'p', b'a')) {
        Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[20, 18, 24, -1, -1]),
        other => panic!("Expected pa kept as is on a crop: {other:?}"),
    }
    match out.data().get(&Tag::new(b'p', b't')) {
        Some(Value::Int32(30)) => {},
        other => panic!("Unexpected pt: {other:?}"),
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

    // Split segment [3,6): kept signal window [18,26], so real positions
    // shift by -18.
    let out = reconstruct_record(&src, 3, 6, 2, 1, true);
    match out.data().get(&Tag::new(b'p', b'a')) {
        Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 0, 6, -1, -1]),
        other => panic!("Expected pa shifted into the subread frame: {other:?}"),
    }
    match out.data().get(&Tag::new(b'p', b't')) {
        Some(Value::Int32(30)) => {}, // base count unchanged
        other => panic!("Unexpected pt: {other:?}"),
    }
}

#[test]
fn update_moves_drops_polya_when_tail_trimmed() {
    let mut src = ubam_with_moves();
    // Anchor at 12 sits in the trimmed-off front signal (kept window is [16,26]).
    src.data_mut().insert(
        Tag::new(b'p', b'a'),
        Value::Array(Array::Int32(vec![12, 10, 14, -1, -1])),
    );
    src.data_mut()
        .insert(Tag::new(b'p', b't'), Value::Int32(30));

    let out = reconstruct_record(&src, 2, 6, 1, 0, true); // head crop 2
    assert!(
        out.data().get(&Tag::new(b'p', b'a')).is_none(),
        "The pa tag is dropped when the tail is trimmed"
    );
    assert!(
        out.data().get(&Tag::new(b'p', b't')).is_none(),
        "The pt tag is dropped when the tail is trimmed"
    );
}

/// `pa` uses signal coordinates and is not a per-base array.
#[test]
fn update_moves_does_not_reslice_read_length_pa() {
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
    // A 5-element `pa` (equal to the read length) with all real positions
    // inside the kept window.
    d.insert(
        Tag::new(b'p', b'a'),
        Value::Array(Array::Int32(vec![4, 2, 6, -1, -1])),
    );

    // Head crop 1 gives window [1,5): kept signal window [2,10]; `pa` survives.
    let out = reconstruct_record(&src, 1, 5, 1, 0, true);
    match out.data().get(&Tag::new(b'p', b'a')) {
        Some(Value::Array(Array::Int32(v))) => {
            assert_eq!(v, &[4, 2, 6, -1, -1], "The pa tag must not be re-sliced")
        },
        other => panic!("Unexpected pa: {other:?}"),
    }
}

/// `--update-moves` without a move table cannot relate signal to sequence,
/// so the signal and poly-A tags are dropped (`parse_move_table` returns
/// `None`, which selects `drop_all`).
#[test]
fn update_moves_without_move_table_drops_signal_and_polya() {
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
    // Split [3,6): kept window [18,26). A range end exactly at `kept_end`
    // (exclusive) survives; an anchor at `kept_end` is outside the window and
    // drops the tags.
    let mk = |pa: Vec<i32>| {
        let mut src = ubam_with_moves();
        src.data_mut()
            .insert(Tag::new(b'p', b'a'), Value::Array(Array::Int32(pa)));
        src.data_mut()
            .insert(Tag::new(b'p', b't'), Value::Int32(30));
        src
    };
    // Range end equal to `kept_end` (26) survives, shifted by -18.
    let kept = reconstruct_record(&mk(vec![20, 18, 26, -1, -1]), 3, 6, 2, 1, true);
    match kept.data().get(&Tag::new(b'p', b'a')) {
        Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 0, 8, -1, -1]),
        other => panic!("Range end at kept_end should survive: {other:?}"),
    }
    // Anchor equal to `kept_end` (26) is outside the window, so the tags are
    // dropped.
    let dropped = reconstruct_record(&mk(vec![26, 18, 24, -1, -1]), 3, 6, 2, 1, true);
    assert!(
        dropped.data().get(&Tag::new(b'p', b'a')).is_none(),
        "Anchor at the exclusive boundary drops the tags"
    );
}

#[test]
fn run_bam_parallel_matches_sequential_as_multiset() {
    use crate::trim::TrimPlan;

    let mk = |threads| Config {
        trim: TrimPlan {
            head: 2,
            tail: 2,
            quality: None,
        },
        threads,
        quiet: true,
        ..Config::default()
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

    // t1: single-threaded BGZF sink, written to a temporary file.
    let dir = tempfile::tempdir().unwrap();
    let p1 = dir.path().join("t1.bam");
    let mut sink1 = crate::io::bam::writer(Some(&p1), &header, false, 6).unwrap();
    run_bam(
        &header,
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut sink1,
        &mk(1),
        &Arc::new(Counters::default()),
    )
    .unwrap();
    sink1.finish().unwrap();
    let b1 = std::fs::read(&p1).unwrap();

    // t8: multithreaded sink to a temporary file (the multithreaded writer
    // needs an owned `Write + Send`).
    let p8 = dir.path().join("t8.bam");
    let mut sink8 = crate::io::bam::writer(Some(&p8), &header, true, 6).unwrap();
    run_bam(
        &header,
        recs.iter().map(|r| Ok(raw_record(r))),
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
        "The t1 and t8 runs must produce the same record set"
    );
}

/// Writer errors remain observable after the bounded channel reaches capacity.
#[test]
fn run_bam_parallel_surfaces_write_error_without_deadlock() {
    use std::io;

    struct FailAfter {
        limit: usize,
        written: usize,
    }

    let cfg = Config {
        threads: 4,
        quiet: true,
        ..Config::default()
    };
    let recs: Vec<anyhow::Result<bam::Record>> = (0..3000)
        .map(|_| anyhow::Ok(bam::Record::default()))
        .collect();

    let mut sink = FailAfter {
        limit: 100,
        written: 0,
    };
    let res = run_bam_parallel(
        recs.into_iter(),
        &cfg,
        &mut sink,
        |_raw, _rec, _cfg, out: &mut Vec<()>| {
            out.push(());
            Ok(())
        },
        Ok,
        |sink, batch: &Vec<()>| -> io::Result<()> {
            if sink.written >= sink.limit {
                return Err(io::Error::new(io::ErrorKind::BrokenPipe, "boom"));
            }
            sink.written += batch.len();
            Ok(())
        },
        &Arc::new(Counters::default()),
    );
    assert!(
        res.is_err(),
        "Write error must surface as Err and must not hang"
    );
}

/// Mirrors `workflow::fastq`'s `parallel_surfaces_parse_error_instead_of_dropping_it`,
/// driving `run_bam_parallel` directly so a malformed upstream record (an
/// `Err` item from the input iterator) is not silently swallowed.
#[test]
fn run_bam_parallel_surfaces_parse_error_instead_of_dropping_it() {
    use std::io;

    struct NullSink;

    let cfg = Config {
        threads: 4,
        quiet: true,
        ..Config::default()
    };
    let good: Vec<anyhow::Result<bam::Record>> =
        (0..5).map(|_| anyhow::Ok(bam::Record::default())).collect();
    let recs = good
        .into_iter()
        .chain(std::iter::once(Err(anyhow::anyhow!("bad record"))));

    let mut sink = NullSink;
    let res = run_bam_parallel(
        recs,
        &cfg,
        &mut sink,
        |_raw, _rec, _cfg, out: &mut Vec<()>| {
            out.push(());
            Ok(())
        },
        Ok,
        |_sink: &mut NullSink, _batch: &Vec<()>| -> io::Result<()> { Ok(()) },
        &Arc::new(Counters::default()),
    );
    assert!(
        res.is_err(),
        "A malformed record must not be dropped on the parallel path"
    );
}

#[test]
fn run_bam_to_fastq_parallel_matches_sequential_as_multiset() {
    use crate::trim::TrimPlan;

    let mk = |threads| Config {
        trim: TrimPlan {
            head: 2,
            tail: 2,
            quality: None,
        },
        threads,
        quiet: true,
        ..Config::default()
    };
    let recs: Vec<RecordBuf> = (0..300)
        .map(|_| ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]))
        .collect();

    let sorted_records = |bytes: &[u8]| {
        let s = String::from_utf8(bytes.to_vec()).unwrap();
        // Records are grouped as 4 consecutive lines rather than split on
        // `@`: a QUAL byte of Phred 31 (ASCII `@`) would corrupt an `@`
        // split. FASTQ records are exactly 4 lines each here, so the
        // re-chunking is lossless.
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(
            lines.len() % 4,
            0,
            "Expected whole 4-line FASTQ records, got {} lines",
            lines.len()
        );
        let mut v: Vec<String> = lines.chunks(4).map(|c| c.join("\n")).collect();
        v.sort();
        v
    };

    let mut a = Vec::new();
    run_bam_to_fastq(
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut a,
        &mk(1),
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let mut b = Vec::new();
    run_bam_to_fastq(
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut b,
        &mk(8),
        &Arc::new(Counters::default()),
    )
    .unwrap();

    assert_eq!(
        sorted_records(&a),
        sorted_records(&b),
        "The t1 and t8 FASTQ outputs must match as a multiset"
    );
}

#[test]
fn interior_adapter_split_reconstructs_mods_per_segment() {
    use crate::adapter::{Adapter, AdapterConfig, Role};
    // Seq (64 bp): [flank1: C + 23 A][adapter GGGGTTTTGGGGTTTT (no C/A)][flank2: C + 23 A].
    // Only two Cs, at positions 0 and 40. `C+m,0,0;` marks both, with ML
    // [100, 200].
    let mut seq = b"CAAAAAAAAAAAAAAAAAAAAAAA".to_vec(); // C at 0
    seq.extend_from_slice(b"GGGGTTTTGGGGTTTT"); // interior adapter, 16 bp
    seq.extend_from_slice(b"CAAAAAAAAAAAAAAAAAAAAAAA"); // C at 40
    let quals = vec![40u8; seq.len()];
    let mut rec = ubam_with_mods(&seq, quals, b"C+m,0,0;", vec![100, 200]);
    rec.data_mut().insert(
        Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
        Value::Int32(seq.len() as i32),
    );

    // BAM-to-FASTQ path (renders MM/ML/MN as header text), adapters active,
    // split on.
    let mut cfg = cfg_bam2fq(None, 0, FastqTags::All);
    cfg.adapters = Some(AdapterConfig {
        adapters: vec![Adapter {
            name: "mid".into(),
            seq: b"GGGGTTTTGGGGTTTT".to_vec(),
            role: Role::Adapter,
        }],
        error_rate: 0.2,
        end_size: 8, // adapter at [24,40) is interior, more than 8 from both ends of 64 bp
        split: true,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    });

    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(
        stats.output_reads, 2,
        "Interior adapter splits into two subreads"
    );
    let s = String::from_utf8(out).unwrap();
    // Each segment keeps exactly its own C mod, renumbered to occurrence 0.
    assert!(
        s.contains("@r1_segment_1\tMM:Z:C+m,0;\tML:B:C,100\tMN:i:24"),
        "Segment 1 mods wrong: {s}"
    );
    assert!(
        s.contains("@r1_segment_2\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:24"),
        "Segment 2 mods wrong: {s}"
    );
}
/// `MN` is written at the smallest integer subtype that fits, so dorado emits
/// `MN:S` for an ordinary-length read. Accepting only `Int32` would make
/// every real mod-bearing record look inconsistent, forcing a full MM/ML
/// rebuild on an untrimmed record and rewriting `MN` as the wider `i`.
#[test]
fn mn_is_recognized_at_every_integer_subtype() {
    use noodles_sam::alignment::record_buf::data::field::Value;

    for v in [
        Value::UInt8(12),
        Value::Int8(12),
        Value::UInt16(12),
        Value::Int16(12),
        Value::UInt32(12),
        Value::Int32(12),
    ] {
        let mut rec = RecordBuf::default();
        *rec.sequence_mut() = b"ACGTACGTACGT".to_vec().into();
        let d = rec.data_mut();
        d.insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(b"C+m,0;".to_vec().into()),
        );
        d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, v.clone());
        assert_eq!(
            inspect_mod_block(&rec, 12),
            ModBlock::Consistent,
            "MN stored as {v:?} should be recognized"
        );
    }
}

/// An untrimmed record takes the pass-through path, so its tags come out
/// exactly as they went in, including `MN`'s storage width.
#[test]
fn untrimmed_mod_record_passes_through_byte_for_byte() {
    use noodles_sam::alignment::record_buf::data::field::Value;

    let mut rec = RecordBuf::default();
    *rec.flags_mut() = noodles_sam::alignment::record::Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = b"ACGTACGTACGT".to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; 12].into();
    let d = rec.data_mut();
    d.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1;".to_vec().into()),
    );
    d.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![200, 201])),
    );
    d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::UInt16(12));

    let out = reconstruct_record(&rec, 0, 12, 1, 0, false);
    assert_eq!(
        out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH),
        Some(&Value::UInt16(12)),
        "MN must keep the subtype it arrived with"
    );
    assert_eq!(out, rec, "An untrimmed record is passed through unchanged");
}

/// Encodes a decoded record as the raw BAM record a production reader yields.
fn raw_record(rec: &RecordBuf) -> bam::Record {
    use noodles_sam::alignment::io::Write as _;

    let header = sam::Header::default();
    let mut bytes = Vec::new();
    {
        let mut w = bam::io::Writer::new(&mut bytes);
        w.write_header(&header).unwrap();
        w.write_alignment_record(&header, rec).unwrap();
        w.try_finish().unwrap();
    }
    let mut r = bam::io::Reader::new(bytes.as_slice());
    r.read_header().unwrap();
    let mut raw = bam::Record::default();
    assert_ne!(r.read_record(&mut raw).unwrap(), 0);
    raw
}

/// A 10-base PacBio-style record: `fi[i]` belongs to base `i`, `ri[i]` to
/// base `9 - i` (reverse kinetics are stored last base first).
fn pacbio_kinetics_record() -> RecordBuf {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(b"r1".into());
    *src.sequence_mut() = b"ACGTACGTAC".to_vec().into();
    *src.quality_scores_mut() = vec![40; 10].into();
    let d = src.data_mut();
    d.insert(
        Tag::new(b'f', b'i'),
        Value::Array(Array::UInt8((0..10).collect())),
    );
    d.insert(
        Tag::new(b'r', b'i'),
        Value::Array(Array::UInt8((10..20).collect())),
    );
    src
}

fn u8_array(rec: &RecordBuf, tag: [u8; 2]) -> Vec<u8> {
    match rec.data().get(&Tag::new(tag[0], tag[1])) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        other => panic!("{}: {other:?}", std::str::from_utf8(&tag).unwrap()),
    }
}

/// Reverse-strand kinetics take array indexes `[len - end, len - start)`,
/// so a head crop removes entries from the end of `ri` and a tail crop from
/// its start, while `fi` is sliced in read order.
#[test]
fn reverse_strand_kinetics_are_sliced_from_the_other_end() {
    let src = pacbio_kinetics_record();
    // Head crop 2: bases 2..=9 keep fi[2..10] and ri[0..8] (bases 9..2).
    let head = reconstruct_record(&src, 2, 10, 1, 0, false);
    assert_eq!(u8_array(&head, *b"fi"), (2..10).collect::<Vec<u8>>());
    assert_eq!(u8_array(&head, *b"ri"), (10..18).collect::<Vec<u8>>());
    // Tail crop 3: bases 0..=6 keep fi[0..7] and ri[3..10] (bases 6..0).
    let tail = reconstruct_record(&src, 0, 7, 1, 0, false);
    assert_eq!(u8_array(&tail, *b"fi"), (0..7).collect::<Vec<u8>>());
    assert_eq!(u8_array(&tail, *b"ri"), (13..20).collect::<Vec<u8>>());
}

#[test]
fn bam2fq_slices_reverse_strand_kinetics_from_the_other_end() {
    let cfg = cfg_bam2fq(None, 2, FastqTags::All); // head crop 2
    let mut out = Vec::new();
    run_bam_to_fastq(
        [Ok(raw_record(&pacbio_kinetics_record()))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(s.contains("\tfi:B:C,2,3,4,5,6,7,8,9\t"), "{s:?}");
    assert!(s.contains("\tri:B:C,10,11,12,13,14,15,16,17\n"), "{s:?}");
}

/// PacBio's `qs:i`/`qe:i` are query coordinates, not a quality, and follow
/// the window in the original read's frame; only a float `qs` is a quality
/// and is recomputed.
#[test]
fn integer_qs_and_qe_follow_the_window() {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(b"r1".into());
    *src.sequence_mut() = b"ACGTAC".to_vec().into();
    *src.quality_scores_mut() = vec![2, 2, 40, 40, 40, 40].into();
    src.data_mut()
        .insert(Tag::new(b'q', b's'), Value::Int32(1200));
    src.data_mut()
        .insert(Tag::new(b'q', b'e'), Value::Int32(1206));

    let out = reconstruct_record(&src, 2, 6, 1, 0, false);
    assert_eq!(
        out.data().get(&Tag::new(b'q', b's')),
        Some(&Value::Int32(1202))
    );
    assert_eq!(
        out.data().get(&Tag::new(b'q', b'e')),
        Some(&Value::Int32(1206))
    );

    let cfg = cfg_bam2fq(None, 2, FastqTags::All);
    let mut fastq = Vec::new();
    run_bam_to_fastq(
        [Ok(raw_record(&src))].into_iter(),
        &mut fastq,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(fastq).unwrap();
    assert!(s.starts_with("@r1\tqs:i:1202\tqe:i:1206\n"), "{s:?}");

    // Dorado's float qs follows the trimmed quality on both outputs.
    src.data_mut()
        .insert(Tag::new(b'q', b's'), Value::Float(20.0));
    let expected = crate::qual::mean_prob_q(&[40, 40, 40, 40]) as f32;
    match reconstruct_record(&src, 2, 6, 1, 0, false)
        .data()
        .get(&Tag::new(b'q', b's'))
    {
        Some(Value::Float(q)) => assert!((q - expected).abs() < 1e-4),
        other => panic!("Unexpected qs: {other:?}"),
    }
    let mut fastq = Vec::new();
    run_bam_to_fastq(
        [Ok(raw_record(&src))].into_iter(),
        &mut fastq,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(fastq).unwrap();
    assert!(s.starts_with("@r1\tqs:f:"), "{s:?}");
    assert!(!s.contains("qs:f:20"), "{s:?}");
}

/// The untrimmed fast path removes a block whose `MN` disagrees with the
/// sequence instead of rewriting `MN`, adds a missing `MN`, and passes a
/// consistent record through raw.
#[test]
fn raw_full_window_removes_malformed_block_and_adds_missing_mn() {
    let cfg = cfg_bam2fq(None, 0, FastqTags::All);
    let counters = Arc::new(Counters::default());

    let out = process_raw_full_window(
        raw_record(&malformed_mod_record("mn_mismatch")),
        &cfg,
        &counters,
    )
    .unwrap();
    let Some(BamOutputRecord::Decoded(rec)) = out else {
        panic!("A malformed block forces a rebuild");
    };
    for t in MOD_TAGS {
        assert!(rec.data().get(&t).is_none(), "{t:?} must be removed");
    }
    assert!(rec.data().get(&Tag::READ_GROUP).is_some());
    assert_eq!(counters.malformed_mod_reads.load(Ordering::Relaxed), 1);

    let missing = ubam_with_mods(b"CCCA", vec![40; 4], b"C+m,0,0,0;", vec![5, 6, 7]);
    let out = process_raw_full_window(raw_record(&missing), &cfg, &counters).unwrap();
    let Some(BamOutputRecord::Decoded(rec)) = out else {
        panic!("A missing MN forces a rebuild");
    };
    assert_eq!(
        rec.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH),
        Some(&Value::Int32(4))
    );
    assert_eq!(
        counters.malformed_mod_reads.load(Ordering::Relaxed),
        1,
        "A missing MN is not a defect"
    );

    let out =
        process_raw_full_window(raw_record(&read2_with_mods_and_rg()), &cfg, &counters).unwrap();
    assert!(matches!(out, Some(BamOutputRecord::Raw(_))));
}

/// The move table resolves the signal end only to the stride, so a window
/// that runs to the last base keeps the source `ns`.
#[test]
fn update_moves_window_to_the_last_base_keeps_the_source_ns() {
    let mut src = ubam_with_moves();
    // One sample past the last full stride block (10 + 8*2 = 26).
    src.data_mut()
        .insert(Tag::new(b'n', b's'), Value::Int32(27));
    let ns = |rec: &RecordBuf| rec.data().get(&Tag::new(b'n', b's')).cloned();

    let head = reconstruct_record(&src, 2, 6, 1, 0, true);
    assert_eq!(ns(&head), Some(Value::Int32(27)), "Head crop keeps ns");
    let tail = reconstruct_record(&src, 0, 4, 1, 0, true);
    assert_eq!(
        ns(&tail),
        Some(Value::Int32(22)),
        "Tail crop ends at a block"
    );
    // The last split segment spans [18, 27).
    let last = reconstruct_record(&src, 3, 6, 2, 1, true);
    assert_eq!(ns(&last), Some(Value::Int32(9)));
}

/// An empty window at the sequence end has no start base; the signal tags
/// are dropped rather than the process aborted.
#[test]
fn update_moves_empty_window_at_the_end_drops_signal_tags() {
    let out = reconstruct_record(&ubam_with_moves(), 6, 6, 1, 0, true);
    assert!(out.sequence().as_ref().is_empty());
    for t in [b"mv", b"ts", b"ns"] {
        assert!(
            out.data().get(&Tag::new(t[0], t[1])).is_none(),
            "{} dropped",
            std::str::from_utf8(t).unwrap()
        );
    }
}

/// PacBio's fixed-size arrays are never per-base, even on a read whose
/// length equals their element count.
#[test]
fn fixed_size_pacbio_arrays_are_not_sliced() {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(b"r1".into());
    *src.sequence_mut() = b"ACGT".to_vec().into();
    *src.quality_scores_mut() = vec![40; 4].into();
    let d = src.data_mut();
    d.insert(
        Tag::new(b's', b'n'),
        Value::Array(Array::Float(vec![1.0, 2.0, 3.0, 4.0])),
    );
    d.insert(
        Tag::new(b'a', b'c'),
        Value::Array(Array::Int32(vec![0, 1, 1, 0])),
    );
    let out = reconstruct_record(&src, 1, 3, 1, 0, false);
    assert_eq!(
        out.data().get(&Tag::new(b's', b'n')),
        src.data().get(&Tag::new(b's', b'n'))
    );
    assert_eq!(
        out.data().get(&Tag::new(b'a', b'c')),
        src.data().get(&Tag::new(b'a', b'c'))
    );

    let mut two = RecordBuf::default();
    *two.flags_mut() = Flags::UNMAPPED;
    *two.name_mut() = Some(b"r2".into());
    *two.sequence_mut() = b"AC".to_vec().into();
    *two.quality_scores_mut() = vec![40; 2].into();
    two.data_mut().insert(
        Tag::new(b'b', b'c'),
        Value::Array(Array::UInt16(vec![3, 7])),
    );
    let cfg = cfg_bam2fq(None, 1, FastqTags::All);
    let mut out = Vec::new();
    run_bam_to_fastq(
        [Ok(raw_record(&two))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(s.contains("\tbc:B:S,3,7\n"), "{s:?}");
}

/// A split product is not the sequencer's read, so `rn` is -1 on both
/// outputs whether or not the move table is rewritten; a crop keeps it.
#[test]
fn split_sets_rn_to_minus_one_without_update_moves() {
    let mut src = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0;", vec![10]);
    src.data_mut().insert(Tag::new(b'r', b'n'), Value::Int32(7));
    let rn = |rec: &RecordBuf| rec.data().get(&Tag::new(b'r', b'n')).cloned();

    assert_eq!(
        rn(&reconstruct_record(&src, 0, 2, 2, 0, false)),
        Some(Value::Int32(-1))
    );
    assert_eq!(
        rn(&reconstruct_record(&src, 1, 4, 1, 0, false)),
        Some(Value::Int32(7))
    );

    let cfg = cfg_bam2fq(
        Some(QualityOp::Split {
            cutoff: 20,
            window: 1,
        }),
        0,
        FastqTags::All,
    );
    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        [Ok(raw_record(&src))].into_iter(),
        &mut out,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(stats.output_reads, 2);
    let s = String::from_utf8(out).unwrap();
    assert!(s.contains("@r1_segment_1\trn:i:-1\t"), "{s:?}");
    assert!(s.contains("@r1_segment_2\trn:i:-1\t"), "{s:?}");
}

/// The rebuilt record keeps its aux tags in source order: rewritten tags
/// stay in place, removed ones leave no hole, added ones are appended.
#[test]
fn rebuilt_record_keeps_aux_order() {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(b"r1".into());
    *src.sequence_mut() = b"CCAC".to_vec().into();
    *src.quality_scores_mut() = vec![40; 4].into();
    let d = src.data_mut();
    d.insert(Tag::READ_GROUP, Value::String(b"grp".as_slice().into()));
    d.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1;".to_vec().into()),
    );
    d.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![10, 20])),
    );
    d.insert(
        Tag::new(b'm', b'v'),
        Value::Array(Array::Int8(vec![5, 1, 1, 1, 1])),
    );
    d.insert(
        Tag::new(b'i', b'p'),
        Value::Array(Array::UInt8(vec![1, 2, 3, 4])),
    );
    d.insert(Tag::new(b'z', b'z'), Value::Int32(5));

    let out = reconstruct_record(&src, 2, 4, 1, 0, false);
    let order: Vec<[u8; 2]> = out.data().iter().map(|(t, _)| <[u8; 2]>::from(t)).collect();
    assert_eq!(
        order,
        [*b"RG", *b"MM", *b"ML", *b"ip", *b"zz", *b"MN"],
        "The mv tag is removed in place and MN appended"
    );
    assert_eq!(
        out.data().get(&Tag::new(b'z', b'z')),
        Some(&Value::Int32(5))
    );
}

const HIFI: &[u8] = b"m64011_190830_220126/123/ccs";

/// Catalog barcode BC01.
pub(super) const BC01: &[u8] = b"AAGAAAGTTGTCGGTGTCTTTGTG";

/// A 10-base PacBio HiFi record: `qs`/`qe` query coordinates 100..110,
/// `rn` passes, the `sa` coverage runs (5 over bases 0..4, 7 over 4..7,
/// 5 over 7..10), the `sm`/`sx` per-base counts and the `ds`/`ls` undo
/// blobs.
fn pacbio_hifi_record(name: &[u8]) -> RecordBuf {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(name.into());
    *src.sequence_mut() = b"ACGTACGTAC".to_vec().into();
    *src.quality_scores_mut() = vec![40; 10].into();
    let d = src.data_mut();
    d.insert(Tag::new(b'q', b's'), Value::Int32(100));
    d.insert(Tag::new(b'q', b'e'), Value::Int32(110));
    d.insert(Tag::new(b'r', b'n'), Value::Int32(3));
    d.insert(
        Tag::new(b's', b'a'),
        Value::Array(Array::UInt32(vec![4, 5, 3, 7, 3, 5])),
    );
    d.insert(
        Tag::new(b's', b'm'),
        Value::Array(Array::UInt8((0..10).collect())),
    );
    d.insert(
        Tag::new(b's', b'x'),
        Value::Array(Array::UInt8((10..20).collect())),
    );
    d.insert(
        Tag::new(b'd', b's'),
        Value::Array(Array::UInt8(vec![1, 2, 3])),
    );
    d.insert(
        Tag::new(b'l', b's'),
        Value::Array(Array::UInt8(vec![4, 5, 6])),
    );
    src
}

/// A dorado record: `ubam_with_moves` (`mv`, `ts` 10, `ns` 26, `st`, `du`
/// 5 s) plus `qs:f`, `me` and `er`, without `pi`.
fn ont_record() -> RecordBuf {
    let mut src = ubam_with_moves();
    let d = src.data_mut();
    d.insert(Tag::new(b'q', b's'), Value::Float(20.0));
    d.insert(Tag::new(b'm', b'e'), Value::Int32(12));
    d.insert(
        Tag::new(b'e', b'r'),
        Value::String(b"signal_positive".as_slice().into()),
    );
    src
}

fn tag(rec: &RecordBuf, t: [u8; 2]) -> Option<Value> {
    rec.data().get(&Tag::new(t[0], t[1])).cloned()
}

fn string_value(s: &[u8]) -> Option<Value> {
    Some(Value::String(s.into()))
}

fn name_of(rec: &RecordBuf) -> Vec<u8> {
    rec.name().map(|n| n.to_vec()).unwrap_or_default()
}

/// Splits at a low-quality base and carries every tag.
fn split_cfg() -> Config {
    cfg_bam2fq(
        Some(QualityOp::Split {
            cutoff: 20,
            window: 1,
        }),
        0,
        FastqTags::All,
    )
}

/// Runs the BAM-to-FASTQ workflow and returns the stats and the text.
fn bam2fq(recs: Vec<RecordBuf>, cfg: &Config) -> (Stats, String) {
    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut out,
        cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    (stats, String::from_utf8(out).unwrap())
}

/// Runs the BAM-to-BAM workflow and returns the stats and the decoded
/// output records.
fn bam2bam(recs: Vec<RecordBuf>, cfg: &Config) -> (Stats, Vec<RecordBuf>) {
    let header = sam::Header::default();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("o.bam");
    let mut sink = crate::io::bam::writer(Some(&path), &header, false, 6).unwrap();
    let stats = run_bam(
        &header,
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut sink,
        cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    sink.finish().unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let mut reader = noodles_bam::io::Reader::new(bytes.as_slice());
    let h = reader.read_header().unwrap();
    let mut out = Vec::new();
    let mut buf = RecordBuf::default();
    while reader.read_record_buf(&h, &mut buf).unwrap() != 0 {
        out.push(buf.clone());
    }
    (stats, out)
}

/// An integer `qs` or a PacBio read name classifies a record as PacBio;
/// anything else is ONT.
#[test]
fn platform_follows_the_qs_type_and_the_read_name() {
    assert_eq!(platform(&pacbio_hifi_record(HIFI)), Platform::PacBio);
    assert_eq!(
        platform(&pacbio_hifi_record(b"hifi_1")),
        Platform::PacBio,
        "An integer qs alone classifies"
    );
    assert_eq!(platform(&ont_record()), Platform::Ont);
    let named = |name: &[u8]| {
        let mut rec = RecordBuf::default();
        *rec.name_mut() = Some(name.into());
        platform(&rec)
    };
    let pacbio: [&[u8]; 6] = [
        b"m1/123/ccs",
        b"m1/123/ccs/fwd",
        b"m1/123/ccs/rev",
        b"m1/123/ccs/100_200",
        b"m1/123/ccs/fwd/100_200",
        b"m1/123/100_200",
    ];
    for name in pacbio {
        assert_eq!(
            named(name),
            Platform::PacBio,
            "{}",
            String::from_utf8_lossy(name)
        );
    }
    let ont: [&[u8]; 9] = [
        b"r1",
        b"m1/abc/ccs",
        b"m1/123",
        b"m1/123/ccs/x",
        b"/123/ccs",
        b"m1/123/100_",
        b"m1/123/ccs/fwd/rev",
        b"m1/123/ccs/100_200/fwd",
        b"0123e4-uuid_segment_1",
    ];
    for name in ont {
        assert_eq!(
            named(name),
            Platform::Ont,
            "{}",
            String::from_utf8_lossy(name)
        );
    }
}

/// One integer coordinate without the other leaves both unchanged.
#[test]
fn lone_integer_qs_is_left_unchanged() {
    let mut src = pacbio_hifi_record(HIFI);
    src.data_mut().remove(&Tag::new(b'q', b'e'));
    let out = reconstruct_record(&src, 2, 7, 1, 0, false);
    assert_eq!(tag(&out, *b"qs"), Some(Value::Int32(100)));
    assert!(tag(&out, *b"qe").is_none());
}

/// A crop keeps the name; a split names PacBio segments by their query
/// coordinates in the spec's forms and ONT segments `_segment_N`.
#[test]
fn split_names_follow_the_platform() {
    let w = |start, end, idx| Window {
        start,
        end,
        idx,
        total: 2,
    };
    let coords = Some((100, 104));
    let crop = Window {
        start: 0,
        end: 4,
        idx: 0,
        total: 1,
    };
    assert_eq!(segment_name(Platform::PacBio, HIFI, crop, coords), HIFI);
    let cases: [(&[u8], &[u8]); 5] = [
        (b"m1/123/ccs", b"m1/123/ccs/100_104"),
        (b"m1/123/ccs/100_110", b"m1/123/ccs/100_104"),
        (b"m1/123/ccs/fwd", b"m1/123/ccs/fwd/100_104"),
        (b"m1/123/ccs/rev/100_110", b"m1/123/ccs/rev/100_104"),
        (b"m1/123/100_110", b"m1/123/100_104"),
    ];
    for (name, want) in cases {
        assert_eq!(
            segment_name(Platform::PacBio, name, w(0, 4, 0), coords),
            want,
            "{}",
            String::from_utf8_lossy(name)
        );
    }
    // Without `qs`/`qe` the interval is offset from the name's own start.
    assert_eq!(
        segment_name(Platform::PacBio, b"m1/123/ccs/50_60", w(6, 10, 1), None),
        b"m1/123/ccs/56_60"
    );
    assert_eq!(
        segment_name(Platform::PacBio, b"m1/123/ccs", w(6, 10, 1), None),
        b"m1/123/ccs/6_10"
    );
    // A PacBio record without a PacBio name, and any ONT record, take the
    // suffix.
    assert_eq!(
        segment_name(Platform::PacBio, b"hifi_1", w(0, 4, 0), coords),
        b"hifi_1_segment_1"
    );
    assert_eq!(
        segment_name(Platform::Ont, b"m1/123/ccs", w(6, 10, 1), None),
        b"m1/123/ccs_segment_2"
    );
    let out = reconstruct_record(&pacbio_hifi_record(b"hifi_1"), 0, 4, 2, 0, false);
    assert_eq!(name_of(&out), b"hifi_1_segment_1");
}

/// Both output paths name a split PacBio read by its coordinates and
/// rewrite `qs`/`qe` per segment.
#[test]
fn pacbio_split_is_named_and_coordinated_on_both_outputs() {
    let mut src = pacbio_hifi_record(HIFI);
    // Base 2 is low quality, so the split is [0,2) and [3,10).
    let mut quals = vec![40; 10];
    quals[2] = 1;
    *src.quality_scores_mut() = quals.into();
    let cfg = split_cfg();
    let (stats, recs) = bam2bam(vec![src.clone()], &cfg);
    assert_eq!(stats.output_reads, 2);
    assert_eq!(name_of(&recs[0]), b"m64011_190830_220126/123/ccs/100_102");
    assert_eq!(name_of(&recs[1]), b"m64011_190830_220126/123/ccs/103_110");
    assert_eq!(tag(&recs[1], *b"qs"), Some(Value::Int32(103)));
    assert_eq!(tag(&recs[1], *b"qe"), Some(Value::Int32(110)));
    let (_, s) = bam2fq(vec![src], &cfg);
    assert!(
        s.contains("@m64011_190830_220126/123/ccs/100_102\tqs:i:100\tqe:i:102\t"),
        "{s:?}"
    );
    assert!(
        s.contains("@m64011_190830_220126/123/ccs/103_110\tqs:i:103\tqe:i:110\t"),
        "{s:?}"
    );
}

#[test]
fn adapter_and_quality_splits_retain_original_barcode_and_pacbio_coordinates() {
    use crate::adapter::{Adapter, AdapterConfig, Role};
    let adapter = b"GGGGTTTTGGGGTTTT";
    // Catalog BC01 at both ends, named by the `BC` call, so the recorded
    // `bi` spans verify: [0, 24) and [160, 184).
    let seq = [
        BC01.to_vec(),
        vec![b'C'; 60],
        adapter.to_vec(),
        vec![b'C'; 60],
        crate::adapter::reverse_complement(BC01),
    ]
    .concat();
    let mut qual = vec![40; seq.len()];
    qual[50..54].fill(2);
    qual[140..144].fill(2);
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(HIFI.into());
    *src.sequence_mut() = seq.into();
    *src.quality_scores_mut() = qual.into();
    src.data_mut()
        .insert(Tag::new(b'q', b's'), Value::Int32(100));
    src.data_mut()
        .insert(Tag::new(b'q', b'e'), Value::Int32(284));
    src.data_mut().insert(
        Tag::new(b'b', b'i'),
        Value::Array(Array::Float(vec![
            100.0, 0.0, 23.0, 100.0, 184.0, 24.0, 100.0,
        ])),
    );
    src.data_mut().insert(
        Tag::new(b'B', b'C'),
        Value::String(b"SQK-NBD114-24_barcode01".as_slice().into()),
    );
    let mut cfg = split_cfg();
    cfg.trim.head = 3;
    cfg.trim.tail = 5;
    cfg.trim.quality = Some(QualityOp::Split {
        cutoff: 9,
        window: 4,
    });
    cfg.filter.min_length = 20;
    cfg.adapters = Some(AdapterConfig {
        adapters: vec![Adapter {
            name: "junction".into(),
            seq: adapter.to_vec(),
            role: Role::Adapter,
        }],
        error_rate: 0.0,
        end_size: 8,
        split: true,
        min_piece: 20,
        candidate_index: std::sync::OnceLock::new(),
    });
    let (stats, recs) = bam2bam(vec![src.clone()], &cfg);
    assert_eq!(stats.output_reads, 3);
    assert_eq!(stats.segments_dropped_short, 1);
    let (_, fastq) = bam2fq(vec![src], &cfg);
    for (rec, (start, end)) in recs.iter().zip([(127, 150), (154, 179), (203, 240)]) {
        let name = format!("m64011_190830_220126/123/ccs/{start}_{end}");
        assert_eq!(name_of(rec), name.as_bytes());
        assert_eq!(tag(rec, *b"qs"), Some(Value::Int32(start)));
        assert_eq!(tag(rec, *b"qe"), Some(Value::Int32(end)));
        assert_eq!(rec.sequence().len(), (end - start) as usize);
        assert!(tag(rec, *b"bi").is_none());
        assert!(fastq.contains(&format!("@{name}\t")), "{fastq}");
    }
}

/// `rn` is a pass count on PacBio and passes through a split, as does
/// pbmarkdup's `du:Z`; an ONT split gets `rn` -1 and loses its float `du`.
#[test]
fn rn_and_du_follow_the_platform_on_a_split() {
    let mut src = pacbio_hifi_record(HIFI);
    src.data_mut().insert(
        Tag::new(b'd', b'u'),
        Value::String(b"dup-info".as_slice().into()),
    );
    let out = reconstruct_record(&src, 0, 4, 2, 0, false);
    assert_eq!(tag(&out, *b"rn"), Some(Value::Int32(3)));
    assert_eq!(tag(&out, *b"du"), string_value(b"dup-info"));
    assert!(tag(&out, *b"pi").is_none(), "pi is a dorado tag");
    let (_, s) = bam2fq(vec![src], &split_cfg());
    assert!(
        s.contains("\trn:i:3\t") && s.contains("\tdu:Z:dup-info"),
        "{s:?}"
    );

    let ont = reconstruct_record(&ont_record(), 0, 3, 2, 0, false);
    assert_eq!(tag(&ont, *b"rn"), Some(Value::Int32(-1)));
    assert!(tag(&ont, *b"du").is_none());
    assert!(tag(&ont, *b"st").is_none());
}

/// `sa` coverage runs are re-encoded for the window, `sm`/`sx` are sliced
/// per base, and dorado's scalar `sm:f` is untouched.
#[test]
fn sa_coverage_is_resliced_and_sm_sx_are_per_base() {
    let src = pacbio_hifi_record(HIFI);
    // Per-base coverage [5,5,5,5,7,7,7,5,5,5]; window [2,7) keeps
    // [5,5,7,7,7].
    let out = reconstruct_record(&src, 2, 7, 1, 0, false);
    assert_eq!(
        tag(&out, *b"sa"),
        Some(Value::Array(Array::UInt32(vec![2, 5, 3, 7])))
    );
    assert_eq!(u8_array(&out, *b"sm"), (2..7).collect::<Vec<u8>>());
    assert_eq!(u8_array(&out, *b"sx"), (12..17).collect::<Vec<u8>>());
    // Window [3,8) keeps [5,7,7,7,5]: equal coverage on either side of a
    // run stays two runs.
    let out = reconstruct_record(&src, 3, 8, 1, 0, false);
    assert_eq!(
        tag(&out, *b"sa"),
        Some(Value::Array(Array::UInt32(vec![1, 5, 3, 7, 1, 5])))
    );
    let (_, s) = bam2fq(vec![src], &cfg_bam2fq(None, 2, FastqTags::All));
    assert!(s.contains("\tsa:B:I,2,5,3,7,3,5\t"), "{s:?}");
    assert!(s.contains("\tsm:B:C,2,3,4,5,6,7,8,9\t"), "{s:?}");

    let mut ont = ont_record();
    ont.data_mut()
        .insert(Tag::new(b's', b'm'), Value::Float(1.5));
    assert!(!has_malformed_perbase_tag(&ont, 6));
    let out = reconstruct_record(&ont, 2, 6, 1, 0, false);
    assert_eq!(tag(&out, *b"sm"), Some(Value::Float(1.5)));
}

/// An `sa` whose runs do not sum to the read length cannot be sliced: it
/// is left unchanged and the read is counted as carrying a malformed
/// per-base tag, on the decoded and the raw path alike.
#[test]
fn malformed_sa_is_left_unchanged_and_counted() {
    let mut src = pacbio_hifi_record(HIFI);
    src.data_mut().insert(
        Tag::new(b's', b'a'),
        Value::Array(Array::UInt32(vec![4, 5, 3, 7])),
    );
    assert!(has_malformed_perbase_tag(&src, 10));
    let out = reconstruct_record(&src, 2, 7, 1, 0, false);
    assert_eq!(
        tag(&out, *b"sa"),
        Some(Value::Array(Array::UInt32(vec![4, 5, 3, 7])))
    );
    let (_, malformed) = raw_full_window_metadata(&raw_record(&src)).unwrap();
    assert!(malformed);
    let (stats, _) = bam2fq(vec![src], &cfg_bam2fq(None, 2, FastqTags::All));
    assert_eq!(stats.malformed_tag_reads, 1);

    let mut odd = pacbio_hifi_record(HIFI);
    odd.data_mut().insert(
        Tag::new(b's', b'a'),
        Value::Array(Array::UInt32(vec![10, 5, 3])),
    );
    assert!(has_malformed_perbase_tag(&odd, 10));

    let ok = pacbio_hifi_record(HIFI);
    assert!(!has_malformed_perbase_tag(&ok, 10));
    let (_, malformed) = raw_full_window_metadata(&raw_record(&ok)).unwrap();
    assert!(!malformed);
}

/// The undo blobs leave every output record of a trimmed read and the read
/// is counted once; an untrimmed read keeps them.
#[test]
fn undo_blobs_are_dropped_on_trim_and_counted_once_per_read() {
    let src = pacbio_hifi_record(HIFI);
    let out = reconstruct_record(&src, 2, 7, 1, 0, false);
    assert!(tag(&out, *b"ds").is_none() && tag(&out, *b"ls").is_none());
    let kept = reconstruct_record(&src, 0, 10, 1, 0, false);
    assert!(tag(&kept, *b"ds").is_some() && tag(&kept, *b"ls").is_some());

    // A cropped read with the blobs, one without, and a split read with
    // them: two reads counted, no output record carries them.
    let mut split = pacbio_hifi_record(b"m1/7/ccs");
    let mut quals = vec![40; 10];
    quals[5] = 1;
    *split.quality_scores_mut() = quals.into();
    let mut plain = pacbio_hifi_record(b"m1/8/ccs");
    plain.data_mut().remove(&Tag::new(b'd', b's'));
    plain.data_mut().remove(&Tag::new(b'l', b's'));
    let mut cfg = split_cfg();
    cfg.trim.head = 1;
    let reads = vec![src.clone(), plain, split];
    let (stats, recs) = bam2bam(reads.clone(), &cfg);
    assert_eq!(stats.output_reads, 4);
    assert_eq!(stats.undo_tags_dropped_reads, 2);
    assert!(
        recs.iter()
            .all(|r| tag(r, *b"ds").is_none() && tag(r, *b"ls").is_none())
    );
    let (stats, s) = bam2fq(reads, &cfg);
    assert_eq!(stats.undo_tags_dropped_reads, 2);
    assert!(!s.contains("ds:B") && !s.contains("ls:B"), "{s:?}");

    let (stats, recs) = bam2bam(vec![src], &cfg_bam2fq(None, 0, FastqTags::All));
    assert_eq!(stats.undo_tags_dropped_reads, 0);
    assert!(tag(&recs[0], *b"ds").is_some() && tag(&recs[0], *b"ls").is_some());
}

/// An ONT split carries the parent id, zero MinKNOW events and an unknown
/// end reason on every subread but the last, with or without
/// `--update-moves`; `sp` still needs the move table and `st`/`du` are
/// dropped without it.
#[test]
fn ont_split_marks_subreads_without_update_moves() {
    let src = ont_record();
    let first = reconstruct_record(&src, 0, 3, 2, 0, false);
    assert_eq!(tag(&first, *b"pi"), string_value(b"r1"));
    assert_eq!(tag(&first, *b"me"), Some(Value::Int32(0)));
    assert_eq!(tag(&first, *b"er"), string_value(b"unknown"));
    assert!(tag(&first, *b"sp").is_none());
    assert!(tag(&first, *b"st").is_none() && tag(&first, *b"du").is_none());
    let last = reconstruct_record(&src, 3, 6, 2, 1, false);
    assert_eq!(tag(&last, *b"pi"), string_value(b"r1"));
    assert_eq!(tag(&last, *b"me"), Some(Value::Int32(0)));
    assert_eq!(tag(&last, *b"er"), string_value(b"signal_positive"));

    // A crop is still the sequencer's read.
    let crop = reconstruct_record(&src, 2, 6, 1, 0, false);
    assert!(tag(&crop, *b"pi").is_none());
    assert_eq!(tag(&crop, *b"me"), Some(Value::Int32(12)));
    assert_eq!(tag(&crop, *b"er"), string_value(b"signal_positive"));

    // A record without `me`/`er` gains neither; `pi` is always set.
    let bare = reconstruct_record(&ubam_with_moves(), 0, 3, 2, 0, false);
    assert!(tag(&bare, *b"me").is_none() && tag(&bare, *b"er").is_none());
    assert_eq!(tag(&bare, *b"pi"), string_value(b"r1"));

    let mut fq = src;
    *fq.quality_scores_mut() = vec![40, 40, 1, 40, 40, 40].into();
    let (stats, s) = bam2fq(vec![fq], &split_cfg());
    assert_eq!(stats.output_reads, 2);
    assert!(
        s.contains("@r1_segment_1\tqs:f:") && s.contains("@r1_segment_2\tqs:f:"),
        "{s:?}"
    );
    assert!(
        s.contains("\tme:i:0\ter:Z:unknown\tpi:Z:r1\trn:i:-1\n"),
        "{s:?}"
    );
    assert!(
        s.contains("\tme:i:0\ter:Z:signal_positive\tpi:Z:r1\trn:i:-1\n"),
        "{s:?}"
    );
    assert!(!s.contains("\tst:Z:") && !s.contains("\tdu:f:"), "{s:?}");
}

/// Under `--update-moves` a split recomputes `du` from the subread's
/// samples and `st` from its start sample, at the rate `ns`/`du` gives
/// (26 samples over 5 s), in the source's offset form.
#[test]
fn update_moves_split_recomputes_st_and_du() {
    // Segment [0,3) spans samples [10,18): 8 samples at 5.2/s last 1.538 s
    // and start 1.923 s into the read.
    let first = reconstruct_record(&ont_record(), 0, 3, 2, 0, true);
    match tag(&first, *b"du") {
        Some(Value::Float(d)) => assert!((d - 1.538_461_5).abs() < 1e-5, "{d}"),
        other => panic!("Unexpected du: {other:?}"),
    }
    assert_eq!(
        tag(&first, *b"st"),
        string_value(b"2024-06-21T10:00:01.923Z")
    );
    // Segment [3,6) spans [18,26): 3.462 s in.
    let last = reconstruct_record(&ont_record(), 3, 6, 2, 1, true);
    assert_eq!(
        tag(&last, *b"st"),
        string_value(b"2024-06-21T10:00:03.462Z")
    );
    assert_eq!(tag(&last, *b"sp"), Some(Value::Int32(18)));
    // A crop keeps both.
    let crop = reconstruct_record(&ont_record(), 2, 6, 1, 0, true);
    assert_eq!(tag(&crop, *b"st"), string_value(b"2024-06-21T10:00:00Z"));
    assert_eq!(tag(&crop, *b"du"), Some(Value::Float(5.0)));

    let with_st = |st: &[u8]| {
        let mut r = ont_record();
        r.data_mut()
            .insert(Tag::new(b's', b't'), Value::String(st.into()));
        tag(&reconstruct_record(&r, 0, 3, 2, 0, true), *b"st")
    };
    assert_eq!(
        with_st(b"2024-06-21T10:00:00.000+00:00"),
        string_value(b"2024-06-21T10:00:01.923+00:00")
    );
    assert_eq!(
        with_st(b"2024-06-21T12:00:00.500+02:00"),
        string_value(b"2024-06-21T12:00:02.423+02:00")
    );
    assert_eq!(
        with_st(b"2024-06-21T10:00:00"),
        string_value(b"2024-06-21T10:00:01.923")
    );
    // An `st` that does not parse is left as is.
    assert_eq!(with_st(b"yesterday"), string_value(b"yesterday"));
}

/// The recompute needs a positive float `du` and an `ns`; otherwise both
/// tags are copied unchanged rather than dropped or guessed.
#[test]
fn update_moves_split_leaves_st_and_du_when_the_rate_is_unknown() {
    let variants: [(&[u8; 2], Option<Value>); 3] = [
        (b"du", Some(Value::Float(0.0))),
        (b"du", None),
        (b"ns", None),
    ];
    for (t, v) in variants {
        let mut r = ont_record();
        let t = Tag::new(t[0], t[1]);
        match v {
            Some(v) => {
                r.data_mut().insert(t, v);
            },
            None => {
                r.data_mut().remove(&t);
            },
        }
        let out = reconstruct_record(&r, 0, 3, 2, 0, true);
        assert_eq!(tag(&out, *b"st"), tag(&r, *b"st"), "{t:?}");
        assert_eq!(tag(&out, *b"du"), tag(&r, *b"du"), "{t:?}");
    }
}

/// Builds the removal set the flags produce.
fn removal(tags: &[&str]) -> crate::config::TagRemoval {
    let tags: Vec<String> = tags.iter().map(|t| (*t).to_string()).collect();
    crate::config::TagRemoval::parse(&tags).unwrap()
}

/// An 8-base record carrying a rewritten block (`MM`/`ML`/`MN`), a per-base
/// array (`ip`), a run-length coverage array (`sa`) and a copied scalar
/// (`RG`).
fn record_with_mixed_tags() -> RecordBuf {
    let mut rec = ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]);
    let data = rec.data_mut();
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
    data.insert(
        Tag::new(b'i', b'p'),
        Value::Array(Array::UInt8((0..8).collect())),
    );
    data.insert(
        Tag::new(b's', b'a'),
        Value::Array(Array::UInt32(vec![8, 3])),
    );
    data.insert(
        Tag::new(b'R', b'G'),
        Value::String(b"run1".as_slice().into()),
    );
    rec
}

/// A removed tag is left out of a trimmed record, whether whittle rewrites
/// it (`MM`) or copies it (`RG`). Removal runs after the rewrite, so the
/// rest of the modification block is still rebuilt against the window.
#[test]
fn removal_drops_a_rewritten_and_a_copied_tag_on_a_trimmed_record() {
    let rec = record_with_mixed_tags();
    let mut cfg = cfg_bam2fq(None, 2, FastqTags::All);
    cfg.remove_tags = removal(&["MM", "RG"]);
    let (_stats, out) = bam2bam(vec![rec], &cfg);

    assert_eq!(out.len(), 1);
    let out = &out[0];
    assert_eq!(out.sequence().as_ref(), b"ACCCAC");
    assert!(tag(out, *b"MM").is_none(), "MM was removed");
    assert!(tag(out, *b"RG").is_none(), "RG was removed");
    // The block is still rebuilt for the window: the call at base 0 falls
    // outside it, leaving the two at bases 3 and 4.
    assert_eq!(
        tag(out, *b"ML"),
        Some(Value::Array(Array::UInt8(vec![20, 30])))
    );
    assert_eq!(tag(out, *b"MN"), Some(Value::Int32(6)));
    // An untouched per-base array is still sliced to the window.
    assert_eq!(
        tag(out, *b"ip"),
        Some(Value::Array(Array::UInt8((2..8).collect())))
    );
}

/// Removal applies to a record no trimming would otherwise change, which is
/// what routes the raw full-window path through the decoded rebuild.
#[test]
fn removal_applies_to_an_untrimmed_record() {
    let rec = record_with_mixed_tags();
    let mut cfg = cfg_bam2fq(None, 0, FastqTags::All);
    cfg.remove_tags = removal(&["RG"]);
    let (_stats, out) = bam2bam(vec![rec.clone()], &cfg);

    assert_eq!(out.len(), 1);
    assert!(tag(&out[0], *b"RG").is_none(), "RG was removed");
    // Every other tag rides through unchanged, values included.
    for t in [*b"MM", *b"ML", *b"MN", *b"ip", *b"sa"] {
        assert_eq!(
            tag(&out[0], t),
            tag(&rec, t),
            "{}",
            String::from_utf8_lossy(&t)
        );
    }
    assert_eq!(out[0].sequence().as_ref(), rec.sequence().as_ref());
}

/// The `kinetics` group removes all nine per-base arrays and nothing else.
#[test]
fn kinetics_group_removes_every_per_base_array() {
    let mut rec = ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]);
    let names: [[u8; 2]; 9] = [
        *b"ip", *b"pw", *b"fi", *b"fp", *b"ri", *b"rp", *b"sm", *b"sx", *b"sa",
    ];
    {
        let data = rec.data_mut();
        for name in names {
            let value = if name == *b"sa" {
                Value::Array(Array::UInt32(vec![8, 3]))
            } else {
                Value::Array(Array::UInt8((0..8).collect()))
            };
            data.insert(Tag::new(name[0], name[1]), value);
        }
        data.insert(
            Tag::new(b'R', b'G'),
            Value::String(b"run1".as_slice().into()),
        );
    }

    let mut cfg = cfg_bam2fq(None, 0, FastqTags::All);
    cfg.remove_tags = removal(&["kinetics"]);
    let (_stats, out) = bam2bam(vec![rec], &cfg);

    assert_eq!(out.len(), 1);
    for name in names {
        assert!(
            tag(&out[0], name).is_none(),
            "{} was removed",
            String::from_utf8_lossy(&name)
        );
    }
    assert_eq!(
        tag(&out[0], *b"RG"),
        string_value(b"run1"),
        "An unnamed tag is untouched"
    );
    assert!(tag(&out[0], *b"MM").is_some());
}

/// On BAM-to-FASTQ the removal applies to the tags carried into the header,
/// per field of the modification block as it does on BAM output.
#[test]
fn removal_drops_the_tag_from_a_fastq_header() {
    let rec = record_with_mixed_tags();
    let mut cfg = cfg_bam2fq(None, 2, FastqTags::All);
    cfg.remove_tags = removal(&["MM", "RG"]);
    let (_stats, text) = bam2fq(vec![rec.clone()], &cfg);

    let header = text.lines().next().unwrap();
    assert!(!header.contains("MM:Z:"), "MM was removed: {header}");
    assert!(!header.contains("RG:Z:"), "RG was removed: {header}");
    assert!(header.contains("ML:B:C,20,30"), "{header}");
    assert!(header.contains("MN:i:6"), "{header}");
    assert!(header.contains("ip:B:C,2,3,4,5,6,7"), "{header}");

    // Without the flag the same run carries both.
    let mut keep = cfg.clone();
    keep.remove_tags = crate::config::TagRemoval::default();
    let (_stats, text) = bam2fq(vec![rec], &keep);
    let header = text.lines().next().unwrap();
    assert!(header.contains("MM:Z:"), "{header}");
    assert!(header.contains("RG:Z:run1"), "{header}");
}
