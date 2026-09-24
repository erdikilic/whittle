//! ONT signal tags under trimming and splitting: move tables, sample counts,
//! poly(A) and subread tags.

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
fn default_drops_the_move_tags_and_keeps_the_parent_link_on_a_crop() {
    let mut src = ubam_with_moves();
    src.data_mut().insert(Tag::new(b's', b'p'), Value::Int32(5));
    src.data_mut().insert(
        Tag::new(b'p', b'i'),
        Value::String(b"parent".as_slice().into()),
    );

    // `update_moves` off and cropped: mv/ts/ns are removed, and sp/pi, which
    // place the unchanged raw signal in its parent, are kept.
    let out = reconstruct_record(&src, 2, 6, 1, 0, false);
    for t in [b"mv", b"ts", b"ns"] {
        assert!(
            out.data().get(&Tag::new(t[0], t[1])).is_none(),
            "{} must be dropped by default on trim",
            std::str::from_utf8(t).unwrap()
        );
    }
    assert_eq!(
        out.data().get(&Tag::new(b's', b'p')),
        Some(&Value::Int32(5))
    );
    assert_eq!(
        out.data().get(&Tag::new(b'p', b'i')),
        Some(&Value::String(b"parent".as_slice().into()))
    );
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

/// A tail crop under `update_moves` shortens `ns` and scales `du` with it, so
/// `ns` over `du` stays the sample rate; `st` is kept.
#[test]
fn update_moves_tail_crop_scales_du_with_ns() {
    let out = reconstruct_record(&ont_record(), 0, 4, 1, 0, true);
    assert_eq!(tag(&out, *b"ns"), Some(Value::Int32(22)));
    let Some(Value::Float(du)) = tag(&out, *b"du") else {
        panic!("du missing");
    };
    assert!((f64::from(du) - 5.0 * 22.0 / 26.0).abs() < 1e-6, "{du}");
    assert_eq!(tag(&out, *b"st"), string_value(b"2024-06-21T10:00:00Z"));
}
