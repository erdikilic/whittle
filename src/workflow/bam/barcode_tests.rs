
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

use super::*;

/// A 20-base record carrying `bi` with the given seven floats.
fn record_with_bi(values: Vec<f32>) -> RecordBuf {
    let mut rec = plain_record();
    rec.data_mut()
        .insert(Tag::new(b'b', b'i'), Value::Array(Array::Float(values)));
    rec
}

/// A 20-base unmapped record with no barcode tag.
fn plain_record() -> RecordBuf {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = b"ACGTACGTACGTACGTACGT".to_vec().into();
    *rec.quality_scores_mut() = vec![40; 20].into();
    rec
}

/// Dorado's `bi` layout is `[barcode_score, front_start_index, front_len,
/// front_score, rear_end_index, rear_len, rear_score]`, and
/// `front_start_index + front_len` is the last base of the front barcode,
/// so the first kept base is one past it.
#[test]
fn window_follows_dorados_trim_interval() {
    let rec = record_with_bi(vec![90.0, 0.0, 4.0, 88.0, 18.0, 3.0, 87.0]);
    assert_eq!(
        barcode_window(&rec, 20),
        BarcodeSpan::Spans {
            front: Some((0, 5)),
            rear: Some((15, 18))
        }
    );
}

/// A front barcode that does not start at base 0 still ends at
/// `front_start_index + front_len`.
#[test]
fn window_honors_a_front_barcode_offset_from_the_read_start() {
    let rec = record_with_bi(vec![90.0, 2.0, 4.0, 88.0, 18.0, 3.0, 87.0]);
    assert_eq!(
        barcode_window(&rec, 20),
        BarcodeSpan::Spans {
            front: Some((2, 7)),
            rear: Some((15, 18))
        }
    );
}

/// Dorado writes `-1` positions for a barcode it did not find, which leaves
/// that end of the read where it was.
#[test]
fn a_missing_barcode_leaves_its_end_alone() {
    let front_only = record_with_bi(vec![90.0, 0.0, 4.0, 88.0, -1.0, 0.0, -1.0]);
    assert_eq!(
        barcode_window(&front_only, 20),
        BarcodeSpan::Spans {
            front: Some((0, 5)),
            rear: None
        }
    );

    let rear_only = record_with_bi(vec![90.0, -1.0, 0.0, -1.0, 18.0, 3.0, 87.0]);
    assert_eq!(
        barcode_window(&rear_only, 20),
        BarcodeSpan::Spans {
            front: None,
            rear: Some((15, 18))
        }
    );

    let neither = record_with_bi(vec![0.0, -1.0, 0.0, -1.0, -1.0, 0.0, -1.0]);
    assert_eq!(
        barcode_window(&neither, 20),
        BarcodeSpan::Spans {
            front: None,
            rear: None
        }
    );
}

#[test]
fn a_record_without_the_tag_is_absent() {
    assert_eq!(barcode_window(&plain_record(), 20), BarcodeSpan::Absent);
}

/// A tag that is not a seven-element float array, or that describes a
/// window outside the read, cannot place the barcodes on the sequence.
#[test]
fn unusable_tags_are_malformed() {
    let short = record_with_bi(vec![90.0, 0.0, 4.0, 88.0, 18.0, 3.0]);
    assert_eq!(barcode_window(&short, 20), BarcodeSpan::Malformed);

    let long = record_with_bi(vec![90.0, 0.0, 4.0, 88.0, 18.0, 3.0, 87.0, 1.0]);
    assert_eq!(barcode_window(&long, 20), BarcodeSpan::Malformed);

    let mut wrong_subtype = plain_record();
    wrong_subtype.data_mut().insert(
        Tag::new(b'b', b'i'),
        Value::Array(Array::Int32(vec![90, 0, 4, 88, 18, 3, 87])),
    );
    assert_eq!(barcode_window(&wrong_subtype, 20), BarcodeSpan::Malformed);

    let not_an_array = {
        let mut rec = plain_record();
        rec.data_mut().insert(Tag::new(b'b', b'i'), Value::Int32(7));
        rec
    };
    assert_eq!(barcode_window(&not_an_array, 20), BarcodeSpan::Malformed);

    // A rear barcode longer than its own end index puts the window end
    // before base 0.
    let negative_end = record_with_bi(vec![90.0, 0.0, 4.0, 88.0, 2.0, 5.0, 87.0]);
    assert_eq!(barcode_window(&negative_end, 20), BarcodeSpan::Malformed);

    // The front barcode ends past the rear one.
    let inverted = record_with_bi(vec![90.0, 0.0, 16.0, 88.0, 12.0, 3.0, 87.0]);
    assert_eq!(barcode_window(&inverted, 20), BarcodeSpan::Malformed);

    // The rear position names a base the read does not have.
    let past_end = record_with_bi(vec![90.0, 0.0, 4.0, 88.0, 40.0, 3.0, 87.0]);
    assert_eq!(barcode_window(&past_end, 20), BarcodeSpan::Malformed);

    // A window covering the whole read leaves nothing to keep.
    let empty = record_with_bi(vec![90.0, 0.0, 19.0, 88.0, 19.0, 0.0, 87.0]);
    assert_eq!(barcode_window(&empty, 20), BarcodeSpan::Malformed);
}

#[test]
fn a_non_finite_position_is_malformed() {
    let nan = record_with_bi(vec![90.0, 0.0, f32::NAN, 88.0, 18.0, 3.0, 87.0]);
    assert_eq!(barcode_window(&nan, 20), BarcodeSpan::Malformed);
}

/// Each read with an unusable `bi` is counted once and keeps every base for
/// the rest of the trim.
#[test]
fn a_malformed_tag_is_counted_once_and_leaves_the_read_untrimmed() {
    let counters = Counters::default();
    let cfg = cfg_barcodes();
    for values in [
        vec![90.0, 0.0, 4.0, 88.0, 18.0, 3.0],
        vec![90.0, 0.0, 4.0, 88.0, 2.0, 5.0, 87.0],
        vec![90.0, 0.0, 16.0, 88.0, 12.0, 3.0, 87.0],
        vec![90.0, 0.0, 4.0, 88.0, 40.0, 3.0, 87.0],
    ] {
        let rec = record_with_bi(values);
        let prepared = prepare_read(&rec, &cfg, &counters).unwrap();
        assert_eq!(prepared.barcode, None);
    }
    assert_eq!(
        counters.barcode_tag_malformed_reads.load(Ordering::Relaxed),
        4
    );

    // A record without the tag passes through without being counted.
    let rec = plain_record();
    assert_eq!(prepare_read(&rec, &cfg, &counters).unwrap().barcode, None);
    assert_eq!(
        counters.barcode_tag_malformed_reads.load(Ordering::Relaxed),
        4
    );
}

/// Barcode positions are read only with an adapter source, and a span is
/// trimmed only when a barcode sequence is found at it: the catalog entry
/// named by the `BC` call, or a barcode of the configured set.
#[test]
fn spans_are_trimmed_only_when_a_barcode_is_found_at_them() {
    let counters = Counters::default();
    let rec = barcoded_record(true);
    let mut off = cfg_barcodes();
    off.adapters = None;
    assert_eq!(prepare_read(&rec, &off, &counters).unwrap().barcode, None);

    let cfg = cfg_barcodes();
    assert_eq!(
        prepare_read(&rec, &cfg, &counters).unwrap().barcode,
        Some((24, 36))
    );
    assert_eq!(
        counters
            .barcode_tag_unverified_reads
            .load(Ordering::Relaxed),
        0
    );

    let mut unnamed = barcoded_record(true);
    unnamed.data_mut().remove(&Tag::new(b'B', b'C'));
    assert_eq!(
        prepare_read(&unnamed, &cfg, &counters).unwrap().barcode,
        None
    );
    assert_eq!(
        counters
            .barcode_tag_unverified_reads
            .load(Ordering::Relaxed),
        1
    );
    let mut configured = cfg_barcodes();
    configured
        .adapters
        .as_mut()
        .unwrap()
        .adapters
        .push(crate::adapter::Adapter {
            name: "BC01".into(),
            seq: BC01.to_vec(),
            role: crate::adapter::Role::Barcode,
        });
    assert_eq!(
        prepare_read(&unnamed, &configured, &counters)
            .unwrap()
            .barcode,
        Some((24, 36))
    );

    // Spans over sequence that holds no barcode, as on input dorado has
    // already trimmed, are left alone.
    let stale = barcoded_record(false);
    assert_eq!(prepare_read(&stale, &cfg, &counters).unwrap().barcode, None);
    assert_eq!(
        counters
            .barcode_tag_unverified_reads
            .load(Ordering::Relaxed),
        2
    );
}

use super::tests::BC01;

/// A 60-base record whose `bi` records 24-base barcodes at both ends and
/// whose `BC` names barcode 1. With `real`, the spans hold BC01 and its
/// reverse complement; otherwise they hold insert sequence.
fn barcoded_record(real: bool) -> RecordBuf {
    let ends: Vec<u8> = if real {
        BC01.to_vec()
    } else {
        b"ACGTTGCAACGTTGCAACGTTGCA".to_vec()
    };
    let seq = [
        ends.clone(),
        b"ACGTACGTACGT".to_vec(),
        crate::adapter::reverse_complement(&ends),
    ]
    .concat();
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.quality_scores_mut() = vec![40; seq.len()].into();
    *rec.sequence_mut() = seq.into();
    rec.data_mut().insert(
        Tag::new(b'b', b'i'),
        Value::Array(Array::Float(vec![90.0, 0.0, 23.0, 88.0, 60.0, 24.0, 87.0])),
    );
    rec.data_mut().insert(
        Tag::new(b'B', b'C'),
        Value::String(b"SQK-NBD114-24_barcode01".as_slice().into()),
    );
    rec
}

/// A base configuration for the decoded BAM workflows with an empty
/// adapter set, which enables the barcode stage.
fn cfg_barcodes() -> Config {
    let mut cfg = super::tests::cfg_bam2fq(None, 0, FastqTags::All);
    cfg.adapters = Some(crate::adapter::AdapterConfig {
        adapters: Vec::new(),
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 20,
        candidate_index: std::sync::OnceLock::new(),
    });
    cfg
}
