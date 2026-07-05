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
    data.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1,0;".to_vec().into()),
    );
    data.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![250, 5, 200])),
    );
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
    assert_eq!(
        a, b,
        "trimmed mod set must equal original filtered to [3, len) offset by 3"
    );
}

/// Multi-mod fixture: one read carrying THREE mod groups — `C+m`, `C+h`, `A+a`
/// (the real dorado shape) — with the C at abs 3 modified by BOTH `m` and `h`,
/// to exercise multiple groups, multiple fundamental bases, AND a same-position
/// double mod all reconstructed through a head+tail crop.
fn write_fixture_multimod(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"read1".into());
    // seq: C at 0,1,3,4,7,9 ; A at 2,6,8.
    *rec.sequence_mut() = b"CCACCGACAC".to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; 10].into();
    let data = rec.data_mut();
    // C+m at C-occ 0,2,5 -> abs 0,3,9 ; C+h at C-occ 2,4 -> abs 3,7 (abs3 shared) ;
    // A+a at A-occ 0,2 -> abs 2,8. ML concatenated in MM-group order.
    data.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1,2;C+h,2,1;A+a,0,1;".to_vec().into()),
    );
    data.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![200, 150, 100, 55, 66, 240, 10])),
    );
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(10));
    w.write_alignment_record(&header, &rec).unwrap();
}

#[test]
fn trimmed_output_multimod_mods_match_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    write_fixture_multimod(&input);

    // head-crop 2, tail-crop 2 -> surviving window [2, 8) on the length-10 read.
    let cfg = chopping::cli::config_for_test(&input, &output, 2, 2);
    chopping::run(cfg).unwrap();

    let (head, tail, len) = (2usize, 2usize, 10usize);
    let tail_start = len - tail;
    let original = hts_mods(&input);
    let expected: Vec<_> = original
        .iter()
        .filter(|(pos, ..)| *pos >= head && *pos < tail_start)
        .map(|&(pos, cb, mb, st, q)| (pos - head, cb, mb, st, q))
        .collect();

    let got = hts_mods(&output);
    let mut a = expected.clone();
    let mut b = got.clone();
    a.sort();
    b.sort();
    assert_eq!(
        a, b,
        "multi-mod trimmed set must equal original filtered to [2,8) offset by 2"
    );

    // Guard against a trivial pass: the surviving set must genuinely be multi-mod
    // (>1 distinct mod code) and non-empty, so this exercises real reconstruction.
    let codes: std::collections::HashSet<char> = b.iter().map(|t| t.2).collect();
    assert!(
        b.len() >= 3 && codes.len() >= 2,
        "expected a non-trivial multi-mod survivor set, got {b:?}"
    );
}

/// Threaded (t=8) variant of `trimmed_output_multimod_mods_match_oracle`: same
/// fixture, same head/tail crop, but driven through the parallel BAM dispatch
/// (`cfg.threads = 8`) via `config_for_test_threads`. `chopping::run` builds a
/// rayon pool regardless of record count, so this still exercises the
/// unordered render -> bounded-channel -> MT-bgzf-writer path end to end; the
/// assertion body is identical to the t1 oracle (order-independent by
/// construction — both sides are sorted before comparison), so this proves
/// MM/ML reconstruction is correct under parallelism, not just under t1.
#[test]
fn trimmed_output_multimod_mods_match_oracle_t8() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    write_fixture_multimod(&input);

    let cfg = chopping::cli::config_for_test_threads(&input, &output, 2, 2, 8);
    chopping::run(cfg).unwrap();

    let (head, tail, len) = (2usize, 2usize, 10usize);
    let tail_start = len - tail;
    let original = hts_mods(&input);
    let expected: Vec<_> = original
        .iter()
        .filter(|(pos, ..)| *pos >= head && *pos < tail_start)
        .map(|&(pos, cb, mb, st, q)| (pos - head, cb, mb, st, q))
        .collect();

    let got = hts_mods(&output);
    let mut a = expected.clone();
    let mut b = got.clone();
    a.sort();
    b.sort();
    assert_eq!(
        a, b,
        "multi-mod mods must survive parallel (t8) reconstruction"
    );
}

/// Multi-record cross-check at t=8: ~200 differently-named reads (each a copy
/// of the multi-mod fixture, but with a per-read-distinct ML payload — see
/// below), decoded and compared per-read via `hts_mods_by_read` rather than as
/// one flattened bag. `run_bam`'s parallel path is unordered (records land in
/// arrival order, not input order), so this specifically stresses that
/// per-read MM/ML reconstruction stays correct even when many records are in
/// flight across the rayon pool concurrently, not just that the aggregate
/// multiset matches.
///
/// Each read's 7 ML bytes are shifted by its own index `i` (`byte.wrapping_add(i
/// as u8)`, `i` in 0..200 so no wraparound), giving every read a distinct qual
/// fingerprint while keeping the same modified positions (MM unchanged) and
/// the same ML count (7, matching the 7 modified positions). Since the
/// per-read comparison below checks read `i`'s *own* decoded quals against
/// its *own* filtered original, a name<->payload mis-pairing bug (e.g. read
/// `i`'s output compared against read `j`'s original mods) would now produce
/// a qual mismatch and fail the test — with the old identical-payload fixture
/// such a mis-pairing was invisible because every read's mods were the same.
#[test]
fn trimmed_output_multimod_mods_match_oracle_t8_many_reads() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");

    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(&input).unwrap());
    w.write_header(&header).unwrap();
    for i in 0..200 {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(format!("read{i}").into_bytes().into());
        *rec.sequence_mut() = b"CCACCGACAC".to_vec().into();
        *rec.quality_scores_mut() = vec![40u8; 10].into();
        let data = rec.data_mut();
        data.insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(b"C+m,0,1,2;C+h,2,1;A+a,0,1;".to_vec().into()),
        );
        // Distinct ML payload per read: shift the base bytes by `i` so read
        // `i`'s decoded quals are unique (see doc comment above).
        let ml: Vec<u8> = [200u8, 150, 100, 55, 66, 240, 10]
            .into_iter()
            .map(|b| b.wrapping_add(i as u8))
            .collect();
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(ml)),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(10));
        w.write_alignment_record(&header, &rec).unwrap();
    }
    w.try_finish().unwrap();

    let cfg = chopping::cli::config_for_test_threads(&input, &output, 2, 2, 8);
    chopping::run(cfg).unwrap();

    let (head, orig_len) = (2usize, 10usize);
    let orig = hts_mods_by_read(&input);
    let out = hts_mods_by_read(&output);
    assert_eq!(
        out.len(),
        200,
        "all 200 reads must survive the t8 parallel run"
    );

    for (name, (_, got_mods)) in &out {
        let (_, orig_mods) = orig
            .get(name)
            .unwrap_or_else(|| panic!("output read {name} has no matching original read"));
        let expected = filter_offset(orig_mods, head, orig_len);
        assert_eq!(
            sorted(&expected),
            sorted(got_mods),
            "t8 mod mismatch for read {name}"
        );
    }
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
    let Some(path) = std::env::var_os("CHOPPING_UBAM") else {
        return;
    };
    let input = std::path::PathBuf::from(path);
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.ubam");

    let cfg = chopping::cli::config_for_test(&input, &output, 10, 10);
    chopping::run(cfg).unwrap();

    let orig = hts_mods_by_read(&input);
    let out = hts_mods_by_read(&output);
    assert!(
        !out.is_empty(),
        "no output reads decoded from {}",
        output.display()
    );

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

/// Threaded (t=8) companion to `real_ubam_oracle_sweep`: same real-data sweep,
/// same head=10/tail=10 crop, but driven through the parallel BAM dispatch.
/// Reads are matched by name (the parallel path is unordered), so this is a
/// real-world-scale spot-check that `-t 8` produces byte-valid, mod-correct
/// output on genuine ONT/dorado data, not just the small synthetic fixtures.
//   CHOPPING_UBAM=data/short_eqread/short_eqread.bam \
//     cargo test --test bam_mods_oracle -- --ignored
#[test]
#[ignore]
fn real_ubam_oracle_sweep_t8() {
    let Some(path) = std::env::var_os("CHOPPING_UBAM") else {
        return;
    };
    let input = std::path::PathBuf::from(path);
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.ubam");

    let cfg = chopping::cli::config_for_test_threads(&input, &output, 10, 10, 8);
    chopping::run(cfg).unwrap();

    let orig = hts_mods_by_read(&input);
    let out = hts_mods_by_read(&output);
    assert!(
        !out.is_empty(),
        "no output reads decoded from {}",
        output.display()
    );
    assert_eq!(
        out.len(),
        orig.len(),
        "t8 run must not drop or duplicate reads"
    );

    for (name, (_, got_mods)) in &out {
        let (orig_len, orig_mods) = orig
            .get(name)
            .unwrap_or_else(|| panic!("output read {name} has no matching original read"));
        let expected = filter_offset(orig_mods, 10, *orig_len);
        assert_eq!(
            sorted(&expected),
            sorted(got_mods),
            "t8 mod mismatch for read {name}"
        );
    }
}
