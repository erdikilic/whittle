//! Conversion of legacy tag forms.

use super::*;

/// `Mm`/`Ml` move to an absent `MM`/`ML` and are dropped beside a present
/// one, as htslib reads them.
#[test]
fn legacy_mod_tags_take_the_current_spelling() {
    let mut rec = RecordBuf::default();
    let data = rec.data_mut();
    data.insert(Tag::new(b'M', b'm'), Value::String(b"C+m,0;".into()));
    data.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![9])),
    );
    data.insert(Tag::new(b'M', b'l'), Value::Array(Array::UInt8(vec![200])));
    assert!(upgrade_legacy_tags(&mut rec));

    let tags: Vec<[u8; 2]> = rec.data().keys().map(<[u8; 2]>::from).collect();
    assert!(tags.contains(b"MM") && tags.contains(b"ML"), "{tags:?}");
    assert!(!tags.contains(b"Mm") && !tags.contains(b"Ml"), "{tags:?}");
    assert_eq!(tag(&rec, *b"MM"), string_value(b"C+m,0;"));
    assert_eq!(tag(&rec, *b"ML"), Some(Value::Array(Array::UInt8(vec![9]))));
    assert!(!upgrade_legacy_tags(&mut rec));
}

/// An integer `qs` without `qe` is a mean quality and becomes a float. A
/// `qs`/`qe` pair or a PacBio read name keeps an integer `qs` a query
/// coordinate.
#[test]
fn integer_ont_qs_becomes_a_float() {
    let mut ont = RecordBuf::default();
    *ont.name_mut() = Some(b"0a1b2c3d-ont-read".into());
    ont.data_mut()
        .insert(Tag::new(b'q', b's'), Value::UInt8(13));
    assert!(upgrade_legacy_tags(&mut ont));
    assert_eq!(tag(&ont, *b"qs"), Some(Value::Float(13.0)));
    assert!(!upgrade_legacy_tags(&mut ont));

    let mut pacbio = pacbio_hifi_record(b"hifi_1");
    assert!(!upgrade_legacy_tags(&mut pacbio));
    assert_eq!(tag(&pacbio, *b"qs"), Some(Value::Int32(100)));

    let mut named = pacbio_hifi_record(HIFI);
    named.data_mut().remove(&Tag::new(b'q', b'e'));
    assert!(!upgrade_legacy_tags(&mut named));
    assert_eq!(tag(&named, *b"qs"), Some(Value::Int32(100)));
}

/// HiFi kinetics that index the whole CCS read are cut to the read's bases
/// `[qs, qe)`: `fi`/`fp` from the first base of the CCS read, `ri`/`rp` from
/// its last. Arrays of the read length are left as they are.
#[test]
fn unclipped_pacbio_arrays_are_cut_to_the_query_coordinates() {
    // A 10-base read kept as [100, 110) of a 115-base read.
    let mut rec = pacbio_hifi_record(HIFI);
    let whole: Vec<u16> = (0..115).collect();
    for t in [*b"fi", *b"ri"] {
        rec.data_mut().insert(
            Tag::new(t[0], t[1]),
            Value::Array(Array::UInt16(whole.clone())),
        );
    }
    assert!(upgrade_legacy_tags(&mut rec));
    assert_eq!(
        tag(&rec, *b"fi"),
        Some(Value::Array(Array::UInt16((100..110).collect())))
    );
    assert_eq!(
        tag(&rec, *b"ri"),
        Some(Value::Array(Array::UInt16((5..15).collect())))
    );
    assert_eq!(
        tag(&rec, *b"sm"),
        Some(Value::Array(Array::UInt8((0..10).collect()))),
        "An array of the read length is left as it is"
    );
    assert!(!upgrade_legacy_tags(&mut rec));
}

/// Without coordinates that span the read, a longer array cannot be placed
/// and is left unchanged, to be counted as malformed.
#[test]
fn arrays_without_spanning_coordinates_are_left_unchanged() {
    let whole: Vec<u16> = (0..115).collect();
    let mut rec = pacbio_hifi_record(HIFI);
    rec.data_mut()
        .insert(Tag::new(b'q', b'e'), Value::Int32(112));
    rec.data_mut().insert(
        Tag::new(b'f', b'i'),
        Value::Array(Array::UInt16(whole.clone())),
    );
    assert!(!upgrade_legacy_tags(&mut rec));
    assert_eq!(tag(&rec, *b"fi"), Some(Value::Array(Array::UInt16(whole))));
}
