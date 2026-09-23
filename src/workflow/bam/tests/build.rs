//! The raw record builder: packed-sequence windows, aux copying, per-base
//! array slicing, tag replacement and the encoder-equivalent validation.

use super::*;

/// Encodes a decoded record with noodles, without the `block_size` prefix.
fn noodles_bytes(rec: &RecordBuf) -> Vec<u8> {
    use noodles_sam::alignment::io::Write as _;

    let mut w = bam::io::Writer::from(Vec::new());
    w.write_alignment_record(&sam::Header::default(), rec)
        .unwrap();
    w.into_inner()[4..].to_vec()
}

/// Reads raw record bytes, given without the `block_size` prefix, as the raw
/// record a production reader yields.
fn raw_from_bytes(bytes: &[u8]) -> bam::Record {
    let mut framed = u32::try_from(bytes.len()).unwrap().to_le_bytes().to_vec();
    framed.extend_from_slice(bytes);
    let mut reader = bam::io::Reader::from(framed.as_slice());
    let mut raw = bam::Record::default();
    assert_ne!(reader.read_record(&mut raw).unwrap(), 0);
    raw
}

/// A 7-base record with an aux field of every scalar type, text fields, a
/// per-base array of every subtype, reverse-strand kinetics and a fixed-size
/// array whose length equals the read length.
fn record_with_every_type() -> RecordBuf {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"read/1".into());
    *rec.sequence_mut() = b"ACGTNAC".to_vec().into();
    *rec.quality_scores_mut() = vec![3, 14, 15, 92, 65, 35, 89].into();
    let d = rec.data_mut();
    d.insert(Tag::READ_GROUP, Value::String(b"rg 1".as_slice().into()));
    d.insert(Tag::new(b'x', b'A'), Value::Character(b'Q'));
    d.insert(Tag::new(b'x', b'c'), Value::Int8(-3));
    d.insert(Tag::new(b'x', b'C'), Value::UInt8(250));
    d.insert(Tag::new(b'x', b's'), Value::Int16(-300));
    d.insert(Tag::new(b'x', b'S'), Value::UInt16(60000));
    d.insert(Tag::new(b'x', b'i'), Value::Int32(-70000));
    d.insert(Tag::new(b'x', b'I'), Value::UInt32(4_000_000_000));
    d.insert(Tag::new(b'x', b'f'), Value::Float(-1.25));
    d.insert(Tag::new(b'x', b'H'), Value::Hex(b"0AFF".as_slice().into()));
    d.insert(
        Tag::new(b'a', b'c'),
        Value::Array(Array::Int8(vec![-1, -2, -3, -4, -5, -6, -7])),
    );
    d.insert(
        Tag::new(b'a', b'C'),
        Value::Array(Array::UInt8(vec![1, 2, 3, 4, 5, 6, 7])),
    );
    d.insert(
        Tag::new(b'a', b's'),
        Value::Array(Array::Int16(vec![-10, -20, -30, -40, -50, -60, -70])),
    );
    d.insert(
        Tag::new(b'a', b'S'),
        Value::Array(Array::UInt16(vec![10, 20, 30, 40, 50, 60, 70])),
    );
    d.insert(
        Tag::new(b'a', b'i'),
        Value::Array(Array::Int32(vec![-1, 0, 1, 2, 3, 4, 5])),
    );
    d.insert(
        Tag::new(b'a', b'I'),
        Value::Array(Array::UInt32(vec![100, 200, 300, 400, 500, 600, 700])),
    );
    d.insert(
        Tag::new(b'a', b'f'),
        Value::Array(Array::Float(vec![0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])),
    );
    d.insert(
        Tag::new(b'r', b'p'),
        Value::Array(Array::UInt8(vec![70, 71, 72, 73, 74, 75, 76])),
    );
    d.insert(
        Tag::new(b's', b'n'),
        Value::Array(Array::Float(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])),
    );
    d.insert(
        Tag::new(b'a', b'x'),
        Value::Array(Array::UInt8(vec![9, 9, 9])),
    );
    rec
}

/// Every window of an odd-length sequence, at odd and even starts and ends,
/// packs to the noodles encoding of the sliced bases, the unused low nibble
/// of an odd-length result included.
#[test]
fn packed_bases_match_the_noodles_packing_at_every_window() {
    let bases = b"ACGTNACGTRYKMSWBDHV=";
    for len in [0, 1, 2, 7, 8, bases.len()] {
        let mut rec = RecordBuf::default();
        *rec.sequence_mut() = bases[..len].to_vec().into();
        *rec.quality_scores_mut() = vec![30; len].into();
        let raw = raw_record(&rec);
        let packed = raw.sequence().as_bytes();
        for start in 0..=len {
            for end in start..=len {
                let mut expected = RecordBuf::default();
                *expected.sequence_mut() = bases[start..end].to_vec().into();
                *expected.quality_scores_mut() = vec![30; end - start].into();
                let expected = raw_record(&expected);
                let mut out = Vec::new();
                push_packed_bases(&mut out, packed, start, end);
                assert_eq!(
                    out,
                    expected.sequence().as_bytes(),
                    "window [{start}, {end}) of {len} bases"
                );
            }
        }
    }
}

/// A packed sequence whose unused low nibble is set is written with it
/// cleared, at a full window and at a window ending on the last base.
#[test]
fn packed_bases_clear_the_unused_nibble() {
    let packed = [0x12, 0x48, 0x1f];
    let mut out = Vec::new();
    push_packed_bases(&mut out, &packed, 0, 5);
    assert_eq!(out, [0x12, 0x48, 0x10]);
    out.clear();
    push_packed_bases(&mut out, &packed, 1, 5);
    assert_eq!(out, [0x24, 0x81]);
    out.clear();
    push_packed_bases(&mut out, &packed, 2, 5);
    assert_eq!(out, [0x48, 0x10]);
}

/// Every window of a record carrying every aux type builds the bytes noodles
/// encodes for the same record decoded and cut by hand: scalars and text
/// copied, per-base arrays of every subtype sliced, `rp` sliced from the other
/// end, and the fixed-size `sn` and short `ax` arrays copied whole.
#[test]
fn built_windows_match_the_noodles_encoding_of_the_cut_record() {
    let src = record_with_every_type();
    let raw = raw_record(&src);
    let keep = TagRemoval::default();
    let len = src.sequence().len();
    for start in 0..len {
        for end in start + 1..=len {
            let edit = RecordEdit {
                start,
                end,
                ..RecordEdit::unchanged(len, &keep)
            };
            let built = build_record(&raw, edit).unwrap();

            let mut expected = src.clone();
            *expected.sequence_mut() = src.sequence().as_ref()[start..end].to_vec().into();
            *expected.quality_scores_mut() =
                src.quality_scores().as_ref()[start..end].to_vec().into();
            if start != 0 || end != len {
                let tags: Vec<Tag> = src.data().keys().collect();
                for tag in tags {
                    let value = expected.data_mut().get_mut(&tag).unwrap();
                    if let Some(sliced) = perbase_slice(tag.into(), value, len, start, end) {
                        *value = sliced;
                    }
                }
            }
            assert_eq!(built, noodles_bytes(&expected), "window [{start}, {end})");
        }
    }
    let rp = |rec: &RecordBuf| u8_array(rec, *b"rp");
    let edit = RecordEdit {
        start: 1,
        end: 3,
        ..RecordEdit::unchanged(len, &keep)
    };
    let out = decode_built(&build_record(&raw, edit).unwrap());
    assert_eq!(rp(&out), [74, 75], "rp is stored last base first");
    assert_eq!(u8_array(&out, *b"aC"), [2, 3]);
    assert_eq!(u8_array(&out, *b"ax"), [9, 9, 9]);
}

/// A record without aux fields builds its sliced core fields only, and a
/// rejection reason on it becomes the single aux field.
#[test]
fn empty_aux_builds_core_fields_and_takes_the_reason() {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(b"r1".into());
    *src.sequence_mut() = b"ACGTA".to_vec().into();
    *src.quality_scores_mut() = vec![10, 20, 30, 40, 50].into();
    let raw = raw_record(&src);
    let keep = TagRemoval::default();

    let edit = RecordEdit {
        start: 1,
        end: 4,
        name: Some(b"r1_segment_1".to_vec()),
        ..RecordEdit::unchanged(5, &keep)
    };
    let out = decode_built(&build_record(&raw, edit).unwrap());
    assert_eq!(out.sequence().as_ref(), b"CGT");
    assert_eq!(out.quality_scores().as_ref(), [20, 30, 40]);
    assert_eq!(name_of(&out), b"r1_segment_1");
    assert!(out.data().is_empty());

    let edit = RecordEdit {
        reason: Some(Reason::TagFilter),
        ..RecordEdit::unchanged(5, &keep)
    };
    let built = build_record(&raw, edit).unwrap();
    let mut expected = src.clone();
    expected.data_mut().insert(
        reject::REASON_TAG,
        Value::String(b"tag_filter".as_slice().into()),
    );
    assert_eq!(built, noodles_bytes(&expected));
}

/// An update replaces its tag in place, a removal update drops it, an update
/// for an absent tag is appended, `--remove-tag` wins over an update, and a
/// rejection reason replaces an existing `wr` in place.
#[test]
fn updates_replace_in_place_append_and_yield_to_removal() {
    let mut src = record_with_every_type();
    src.data_mut()
        .insert(reject::REASON_TAG, Value::String(b"old".as_slice().into()));
    src.data_mut().insert(Tag::new(b'z', b'z'), Value::UInt8(1));
    let raw = raw_record(&src);
    let remove = removal(&["xi"]);
    let edit = RecordEdit {
        updates: vec![
            (Tag::new(b'x', b'C'), Some(Value::Int32(-5))),
            (Tag::new(b'x', b'H'), None),
            (
                Tag::new(b'n', b'w'),
                Some(Value::String(b"new".as_slice().into())),
            ),
            (Tag::new(b'x', b'i'), Some(Value::Int32(1))),
            (
                Tag::new(b'a', b'C'),
                Some(Value::Array(Array::UInt8(vec![1, 2]))),
            ),
        ],
        reason: Some(Reason::TrimmedToNothing),
        ..RecordEdit::unchanged(7, &remove)
    };
    let out = decode_built(&build_record(&raw, edit).unwrap());
    let tags: Vec<[u8; 2]> = out.data().iter().map(|(t, _)| <[u8; 2]>::from(t)).collect();
    let mut expected: Vec<[u8; 2]> = src
        .data()
        .iter()
        .map(|(t, _)| <[u8; 2]>::from(t))
        .filter(|t| t != b"xH" && t != b"xi")
        .collect();
    expected.push(*b"nw");
    assert_eq!(tags, expected);
    assert_eq!(tag(&out, *b"xC"), Some(Value::Int32(-5)));
    assert_eq!(tag(&out, *b"nw"), string_value(b"new"));
    assert_eq!(
        tag(&out, *b"aC"),
        Some(Value::Array(Array::UInt8(vec![1, 2])))
    );
    assert_eq!(tag(&out, *b"wr"), string_value(b"trimmed_to_nothing"));
}

/// The CIGAR overflow tag is not written to a rebuilt record, as the noodles
/// encoder omits it from a decoded record's aux fields.
#[test]
fn cigar_overflow_tag_is_omitted() {
    let mut src = record_with_every_type();
    src.data_mut()
        .insert(Tag::CIGAR, Value::Array(Array::UInt32(vec![7 << 4 | 4])));
    let raw = raw_record(&src);
    let keep = TagRemoval::default();
    let built = build_record(&raw, RecordEdit::unchanged(7, &keep)).unwrap();
    assert_eq!(built, noodles_bytes(&src));
    assert!(decode_built(&built).data().get(&Tag::CIGAR).is_none());
}

/// Inputs the noodles decoder or encoder refuses are refused: an invalid
/// output name, a quality above 93, a copied string with a control character
/// and a duplicated tag.
#[test]
fn builder_refuses_what_noodles_refuses() {
    let keep = TagRemoval::default();
    let src = record_with_every_type();
    let raw = raw_record(&src);
    for name in [b"".as_slice(), b"*", b"a@b", b"a b"] {
        let edit = RecordEdit {
            start: 1,
            end: 3,
            name: Some(name.to_vec()),
            ..RecordEdit::unchanged(7, &keep)
        };
        assert!(build_record(&raw, edit).is_err(), "name {name:?}");
    }

    // QUAL follows the 32-byte fixed fields, the name `read/1` with its NUL
    // and the 4 bytes of 7 packed bases.
    let mut bytes = noodles_bytes(&src);
    bytes[32 + 7 + 4 + 4] = 94;
    let raw = raw_from_bytes(&bytes);
    assert!(build_record(&raw, RecordEdit::unchanged(7, &keep)).is_err());
    let edit = RecordEdit {
        start: 0,
        end: 4,
        ..RecordEdit::unchanged(7, &keep)
    };
    assert!(
        build_record(&raw, edit).is_ok(),
        "Only the written qualities are checked"
    );

    let mut bytes = noodles_bytes(&src);
    bytes.extend_from_slice(b"xZZa\tb\0");
    let raw = raw_from_bytes(&bytes);
    assert!(build_record(&raw, RecordEdit::unchanged(7, &keep)).is_err());
    assert!(
        build_record(&raw, RecordEdit::unchanged(7, &removal(&["xZ"]))).is_ok(),
        "A removed field is not written, so it is not checked"
    );

    let mut bytes = noodles_bytes(&src);
    bytes.extend_from_slice(b"xCC\x07");
    let raw = raw_from_bytes(&bytes);
    assert!(build_record(&raw, RecordEdit::unchanged(7, &keep)).is_err());
}

/// A record written by `write_record_bytes` carries its `block_size` and reads
/// back; a reference sequence ID absent from the header is refused.
#[test]
fn written_record_bytes_read_back_and_unknown_references_are_refused() {
    let src = record_with_every_type();
    let bytes = noodles_bytes(&src);
    let header = sam::Header::default();
    let mut out = Vec::new();
    crate::io::bam::write_record_bytes(&mut out, &header, &bytes).unwrap();
    assert_eq!(out[..4], u32::try_from(bytes.len()).unwrap().to_le_bytes());
    assert_eq!(out[4..], bytes);

    let mut placed = bytes.clone();
    placed[..4].copy_from_slice(&0i32.to_le_bytes());
    assert!(crate::io::bam::write_record_bytes(&mut Vec::new(), &header, &placed).is_err());
    let mut mate = bytes;
    mate[20..24].copy_from_slice(&3i32.to_le_bytes());
    assert!(crate::io::bam::write_record_bytes(&mut Vec::new(), &header, &mate).is_err());
}
