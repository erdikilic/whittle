// Decode-equivalence oracle: build a synthetic uBAM fixture with known C+m mods,
// head-crop it via `chopping::run`, then decode BOTH the original and our output
// with rust-htslib's `basemods_iter()` (an independent implementation of MM/ML
// decoding). Assert the output's per-position (canonical, modified, strand, qual)
// set equals the original's mods filtered to [start, len) and offset by `start`.
// This is independent of MM's multiple valid encodings — it only cares that the
// decoded modification calls match, not the byte-for-byte MM string.
use std::path::Path;

use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

use rust_htslib::bam::{self as hts, Read as _};

/// Build a one-read uBAM: seq with C's, C+m mods at chosen positions, ML bytes.
///
/// noodles-sam 0.85's `Header::builder()` has no bare `set_header` step that
/// takes `Default::default()` for the whole header; `Header::default()` (as
/// already used in `tests/bam_smoke.rs`) is sufficient since the htslib reader
/// only needs a valid, if minimal, header.
fn write_fixture(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"read1".into());
    // 12 bases, C at indices 1,4,7,10.
    *rec.sequence_mut() = b"ACGGCGGCGGCG".to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; 12].into();
    let data = rec.data_mut();
    // modify C occurrences 0,2,3 -> deltas 0,1,0 ; ML 3 bytes.
    data.insert(Tag::BASE_MODIFICATIONS, Value::String(b"C+m,0,1,0;".to_vec().into()));
    data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(vec![250, 5, 200])));
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(12));
    w.write_alignment_record(&header, &rec).unwrap();
}

/// Decode (0-based read pos, canonical, modified, strand, qual) with htslib.
fn hts_mods(path: &Path) -> Vec<(usize, char, char, i32, i32)> {
    let mut reader = hts::Reader::from_path(path).unwrap();
    let mut out = Vec::new();
    for rec in reader.records() {
        let rec = rec.unwrap();
        if let Ok(iter) = rec.basemods_iter() {
            for (pos, m) in iter.flatten() {
                out.push((
                    pos as usize,
                    m.canonical_base as u8 as char,
                    m.modified_base as u8 as char,
                    m.strand,
                    m.qual,
                ));
            }
        }
    }
    out
}

#[test]
fn trimmed_output_mods_match_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    write_fixture(&input);

    // Trim the first 3 bases (head-crop 3) via the library run().
    let cfg = chopping::cli::config_for_test(&input, &output, 3, 0);
    chopping::run(cfg).unwrap();

    let original = hts_mods(&input);
    let start = 3usize;
    let expected: Vec<_> = original
        .iter()
        .filter(|(pos, ..)| *pos >= start)
        .map(|&(pos, cb, mb, st, q)| (pos - start, cb, mb, st, q))
        .collect();

    let got = hts_mods(&output);
    let mut a = expected.clone();
    let mut b = got.clone();
    a.sort();
    b.sort();
    assert_eq!(a, b, "trimmed mod set must equal original filtered to [3, len) offset by 3");
}
