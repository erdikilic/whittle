//! Per-base arrays, run-length coverage, undo blobs, aux order and tag removal
//! on rebuilt records.

use super::*;

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
