// Decode-equivalence oracle: build a synthetic uBAM fixture with known C+m mods,
// head-crop it via `chopping::run`, then decode BOTH the original and our output
// with rust-htslib's `basemods_iter()` (an independent implementation of MM/ML
// decoding). Assert the output's per-position (canonical, modified, strand, qual)
// set equals the original's mods filtered to [start, len) and offset by `start`.
// This is independent of MM's multiple valid encodings — it only cares that the
// decoded modification calls match, not the byte-for-byte MM string.
use std::collections::HashMap;
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

// --- Real-data sweep (Task 7) -----------------------------------------------
//
// Everything below is only exercised when `CHOPPING_UBAM` points at a real
// uBAM (e.g. one of the HG002 subsets under `data/`, or any unaligned BAM with
// MM/ML). It runs a fixed head=10/tail=10 crop — small enough that any real
// long read survives with exactly one output segment, so output read names
// are unchanged and can be matched 1:1 against the input by name — and checks
// every output read's htslib-decoded mods against the original read's mods
// filtered to the crop window and offset, exactly like the synthetic test
// above but per-read over a whole real file instead of one hand-built record.

/// One decoded base-modification call: (0-based read pos, canonical base,
/// modified base, strand, qual) — the same tuple shape `hts_mods` above
/// produces, just keyed per read instead of flattened across a whole file.
type ModCall = (usize, char, char, i32, i32);

/// Decode every read in a BAM/uBAM, keyed by read name, into its (SEQ length,
/// mod calls) — the same `basemods_iter()` decode as `hts_mods`, but split
/// per read so a multi-read real file can be compared read-by-read instead of
/// as one flattened bag.
fn hts_mods_by_read(path: &Path) -> HashMap<String, (usize, Vec<ModCall>)> {
    let mut reader = hts::Reader::from_path(path).unwrap();
    let mut out = HashMap::new();
    for rec in reader.records() {
        let rec = rec.unwrap();
        let name = String::from_utf8_lossy(rec.qname()).into_owned();
        let seq_len = rec.seq_len();
        let mut mods = Vec::new();
        if let Ok(iter) = rec.basemods_iter() {
            for (pos, m) in iter.flatten() {
                mods.push((
                    pos as usize,
                    m.canonical_base as u8 as char,
                    m.modified_base as u8 as char,
                    m.strand,
                    m.qual,
                ));
            }
        }
        out.insert(name, (seq_len, mods));
    }
    out
}

/// Sort a mod-call vector for order-independent comparison.
fn sorted(mods: &[ModCall]) -> Vec<ModCall> {
    let mut v = mods.to_vec();
    v.sort();
    v
}

/// Filter an original read's mod calls to the surviving window of a fixed
/// head/tail crop (`[head, orig_len - head)`, since head == tail == 10 here)
/// and offset positions back to the trimmed read's own coordinate frame.
fn filter_offset(mods: &[ModCall], head: usize, orig_len: usize) -> Vec<ModCall> {
    let tail_start = orig_len.saturating_sub(head);
    mods.iter()
        .filter(|(pos, ..)| *pos >= head && *pos < tail_start)
        .map(|&(pos, cb, mb, st, q)| (pos - head, cb, mb, st, q))
        .collect()
}

// Runs only when a real uBAM is provided, e.g.:
//   CHOPPING_UBAM=data/short_eqread/short_eqread.bam \
//     cargo test --test bam_mods_oracle -- --ignored
#[test]
#[ignore]
fn real_ubam_oracle_sweep() {
    let Some(path) = std::env::var_os("CHOPPING_UBAM") else { return };
    let input = std::path::PathBuf::from(path);
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.ubam");

    let cfg = chopping::cli::config_for_test(&input, &output, 10, 10);
    chopping::run(cfg).unwrap();

    let orig = hts_mods_by_read(&input);
    let out = hts_mods_by_read(&output);
    assert!(!out.is_empty(), "no output reads decoded from {}", output.display());

    for (name, (_, got_mods)) in &out {
        let (orig_len, orig_mods) = orig
            .get(name)
            .unwrap_or_else(|| panic!("output read {name} has no matching original read"));
        let expected = filter_offset(orig_mods, 10, *orig_len);
        assert_eq!(
            sorted(&expected),
            sorted(got_mods),
            "mod mismatch for read {name}"
        );
    }
}
