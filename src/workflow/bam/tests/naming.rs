//! Platform detection and the names and coordinates of split segments.

use super::*;

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
