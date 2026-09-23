//! Sequence, quality and MM/ML/MN reconstruction of trimmed and split records.

use super::*;

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
