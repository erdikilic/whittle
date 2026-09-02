//! Decode-equivalence tests compare synthetic uBAM input and trimmed output with
//! rust-htslib's independent MM/ML decoder. Comparisons use decoded modification
//! calls because MM permits multiple equivalent encodings.

use std::collections::HashMap;
use std::path::Path;

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use rust_htslib::bam::{self as hts, Read as _};

/// Builds a one-read uBAM: a sequence with C bases, `C+m` calls at chosen
/// positions, and ML bytes. `Header::default()` is a valid minimal header for
/// the htslib reader.
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
    // Modified C occurrences 0,2,3 -> deltas 0,1,0; ML has 3 bytes.
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

/// The ASCII base an htslib base-modification code stands for.
fn base_char(code: i32) -> char {
    char::from(u8::try_from(code).expect("htslib base codes are ASCII"))
}

/// Decodes every call as (0-based read position, canonical, modified, strand,
/// qual) with htslib.
fn hts_mods(path: &Path) -> Vec<(usize, char, char, i32, i32)> {
    let mut reader = hts::Reader::from_path(path).unwrap();
    let mut out = Vec::new();
    for rec in reader.records() {
        let rec = rec.unwrap();
        if let Ok(iter) = rec.basemods_iter() {
            for (pos, m) in iter.flatten() {
                out.push((
                    pos as usize,
                    base_char(m.canonical_base),
                    base_char(m.modified_base),
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

    // Trim the first 3 bases (head-crop 3) through the library `run`.
    let cfg = whittle::cli::config_for_test(&input, &output, 3, 0);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).unwrap();

    let original = hts_mods(&input);
    let start = 3usize;
    let expected: Vec<_> = original
        .iter()
        .filter(|(pos, ..)| *pos >= start)
        .map(|&(pos, cb, mb, st, q)| (pos - start, cb, mb, st, q))
        .collect();

    let mut a = expected;
    let mut b = hts_mods(&output);
    a.sort();
    b.sort();
    assert_eq!(
        a, b,
        "Trimmed mod set equals the original filtered to [3, len) and offset by 3"
    );
}

/// Multi-mod fixture: one read carrying three mod groups (`C+m`, `C+h`, `A+a`,
/// the dorado shape), with the C at abs 3 modified by both `m` and `h`, to
/// exercise multiple groups, multiple fundamental bases, and a same-position
/// double mod, all reconstructed through a head and tail crop.
fn write_fixture_multimod(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"read1".into());
    // C at 0,1,3,4,7,9; A at 2,6,8.
    *rec.sequence_mut() = b"CCACCGACAC".to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; 10].into();
    let data = rec.data_mut();
    // C+m at C occurrences 0,2,5 -> abs 0,3,9; C+h at C occurrences 2,4 -> abs
    // 3,7 (abs 3 shared); A+a at A occurrences 0,2 -> abs 2,8. ML is
    // concatenated in MM group order.
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

    // Head-crop 2 and tail-crop 2 leave the window [2, 8) on the length-10 read.
    let cfg = whittle::cli::config_for_test(&input, &output, 2, 2);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).unwrap();

    let (head, tail, len) = (2usize, 2usize, 10usize);
    let tail_start = len - tail;
    let original = hts_mods(&input);
    let expected: Vec<_> = original
        .iter()
        .filter(|(pos, ..)| *pos >= head && *pos < tail_start)
        .map(|&(pos, cb, mb, st, q)| (pos - head, cb, mb, st, q))
        .collect();

    let mut a = expected;
    let mut b = hts_mods(&output);
    a.sort();
    b.sort();
    assert_eq!(
        a, b,
        "Multi-mod trimmed set equals the original filtered to [2,8) and offset by 2"
    );

    // Guard against a trivial pass: the surviving set contains more than one
    // distinct mod code and is non-empty.
    let codes: std::collections::HashSet<char> = b.iter().map(|t| t.2).collect();
    assert!(
        b.len() >= 3 && codes.len() >= 2,
        "Expected a non-trivial multi-mod survivor set, got {b:?}"
    );
}

/// Validates MM/ML reconstruction through the parallel render, channel, and
/// BGZF writer path. Sorting makes the comparison independent of output order.
#[test]
fn trimmed_output_multimod_mods_match_oracle_t8() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    write_fixture_multimod(&input);

    let cfg = whittle::cli::config_for_test_threads(&input, &output, 2, 2, 8);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).unwrap();

    let (head, tail, len) = (2usize, 2usize, 10usize);
    let tail_start = len - tail;
    let original = hts_mods(&input);
    let expected: Vec<_> = original
        .iter()
        .filter(|(pos, ..)| *pos >= head && *pos < tail_start)
        .map(|&(pos, cb, mb, st, q)| (pos - head, cb, mb, st, q))
        .collect();

    let mut a = expected;
    let mut b = hts_mods(&output);
    a.sort();
    b.sort();
    assert_eq!(a, b, "Multi-mod calls survive parallel (t8) reconstruction");
}

/// Compares many parallel records by name using distinct per-read ML payloads,
/// so a record and payload mismatch is observable despite unordered output.
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
        // Distinct ML payload per read: the base bytes are shifted by `i` so
        // read `i`'s decoded quals are unique.
        let shift = u8::try_from(i).expect("Fixture index fits in a byte");
        let ml: Vec<u8> = [200u8, 150, 100, 55, 66, 240, 10]
            .into_iter()
            .map(|b| b.wrapping_add(shift))
            .collect();
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(ml)),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(10));
        w.write_alignment_record(&header, &rec).unwrap();
    }
    w.try_finish().unwrap();

    let cfg = whittle::cli::config_for_test_threads(&input, &output, 2, 2, 8);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).unwrap();

    let (head, orig_len) = (2usize, 10usize);
    let orig = hts_mods_by_read(&input);
    let out = hts_mods_by_read(&output);
    assert_eq!(out.len(), 200, "All 200 reads survive the t8 parallel run");

    for (name, (_, got_mods)) in &out {
        let (_, orig_mods) = orig
            .get(name)
            .unwrap_or_else(|| panic!("Output read {name} has no matching original read"));
        let expected = filter_offset(orig_mods, head, orig_len);
        assert_eq!(
            sorted(&expected),
            sorted(got_mods),
            "Mod mismatch for read {name} at t8"
        );
    }
}

// Optional real-data sweep enabled by `WHITTLE_UBAM`.

/// The real uBAM the ignored sweeps run against. Running them without the
/// variable is a misconfigured invocation, so it fails rather than passing on
/// no data.
fn real_ubam_path() -> std::ffi::OsString {
    std::env::var_os("WHITTLE_UBAM")
        .expect("WHITTLE_UBAM must name a real uBAM to run the ignored oracle sweeps")
}

/// One decoded base-modification call: (0-based read position, canonical base,
/// modified base, strand, qual), the tuple shape `hts_mods` produces, keyed per
/// read rather than flattened across a whole file.
type ModCall = (usize, char, char, i32, i32);

/// Decodes every read in a BAM, keyed by read name, into its (SEQ length, mod
/// calls): the same `basemods_iter()` decode as `hts_mods`, split per read so a
/// multi-read file is compared read by read rather than as one flattened bag.
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
                    base_char(m.canonical_base),
                    base_char(m.modified_base),
                    m.strand,
                    m.qual,
                ));
            }
        }
        out.insert(name, (seq_len, mods));
    }
    out
}

/// Sorts a mod-call vector for order-independent comparison.
fn sorted(mods: &[ModCall]) -> Vec<ModCall> {
    let mut v = mods.to_vec();
    v.sort();
    v
}

/// Filters an original read's mod calls to the surviving window of a fixed
/// head and tail crop (`[head, orig_len - head)`, since `head == tail` here)
/// and offsets positions into the trimmed read's coordinate frame.
fn filter_offset(mods: &[ModCall], head: usize, orig_len: usize) -> Vec<ModCall> {
    let tail_start = orig_len.saturating_sub(head);
    mods.iter()
        .filter(|(pos, ..)| *pos >= head && *pos < tail_start)
        .map(|&(pos, cb, mb, st, q)| (pos - head, cb, mb, st, q))
        .collect()
}

/// Real-data sweep over the uBAM named by `WHITTLE_UBAM`:
///   WHITTLE_UBAM=/path/to/real.ubam cargo test --test bam_mods_oracle -- --ignored
#[test]
#[ignore = "needs WHITTLE_UBAM=<real uBAM>"]
fn real_ubam_oracle_sweep() {
    let path = real_ubam_path();
    let input = std::path::PathBuf::from(path);
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.ubam");

    let cfg = whittle::cli::config_for_test(&input, &output, 10, 10);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).unwrap();

    let orig = hts_mods_by_read(&input);
    let out = hts_mods_by_read(&output);
    assert!(
        !out.is_empty(),
        "No output reads decoded from {}",
        output.display()
    );

    for (name, (_, got_mods)) in &out {
        let (orig_len, orig_mods) = orig
            .get(name)
            .unwrap_or_else(|| panic!("Output read {name} has no matching original read"));
        let expected = filter_offset(orig_mods, 10, *orig_len);
        assert_eq!(
            sorted(&expected),
            sorted(got_mods),
            "Mod mismatch for read {name}"
        );
    }
}

/// Threaded (`-t 8`) companion to `real_ubam_oracle_sweep`: the same real-data
/// sweep and head 10, tail 10 crop through the parallel BAM dispatch. Reads are
/// matched by name, since that path is unordered.
#[test]
#[ignore = "needs WHITTLE_UBAM=<real uBAM>"]
fn real_ubam_oracle_sweep_t8() {
    let path = real_ubam_path();
    let input = std::path::PathBuf::from(path);
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.ubam");

    let cfg = whittle::cli::config_for_test_threads(&input, &output, 10, 10, 8);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).unwrap();

    let orig = hts_mods_by_read(&input);
    let out = hts_mods_by_read(&output);
    assert!(
        !out.is_empty(),
        "No output reads decoded from {}",
        output.display()
    );
    assert_eq!(
        out.len(),
        orig.len(),
        "A t8 run neither drops nor duplicates reads"
    );

    for (name, (_, got_mods)) in &out {
        let (orig_len, orig_mods) = orig
            .get(name)
            .unwrap_or_else(|| panic!("Output read {name} has no matching original read"));
        let expected = filter_offset(orig_mods, 10, *orig_len);
        assert_eq!(
            sorted(&expected),
            sorted(got_mods),
            "Mod mismatch for read {name} at t8"
        );
    }
}

// Adapter inference on uBAM preserves MM/ML coordinates after trimming.

/// The 28-bp adapter planted at the 5' end of every fixture read below: the same
/// SQK-NSK007/LSK109-neighborhood sequence `tests/adapter_cli.rs` plants as
/// `PLANTED_ADAPTER`. Duplicated because each integration-test file compiles as
/// its own binary.
const INFER_MM_ML_ADAPTER: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT";

/// Length of the per-read genomic tail appended after `INFER_MM_ML_ADAPTER`.
const INFER_MM_ML_TAIL_LEN: usize = 150;

/// Minimum modification offset, beyond the accepted inferred-adapter cut range.
const INFER_MM_ML_MOD_MIN_ABS: usize = 70;

/// The same splitmix64 bit mixer as `tests/adapter_cli.rs`'s `splitmix_tail`: a
/// deterministic, per-read, non-periodic ACGT background, so the tail region
/// carries no cross-read k-mer signal that discovery could report as a second
/// adapter.
fn splitmix_tail_infer(i: usize, len: usize) -> Vec<u8> {
    let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        out.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
    }
    out
}

/// 0-based occurrence indices (the Nth `canonical` base in read order) of
/// `canonical` bases at absolute position `>= min_abs`, up to `want` of them.
/// MM counts modified bases by occurrence of the canonical base across the whole
/// read (SAM spec: deltas count skipped occurrences from the read start), so an
/// MM string encodes occurrence indices, not absolute positions.
fn occurrence_indices_after(seq: &[u8], canonical: u8, min_abs: usize, want: usize) -> Vec<usize> {
    let mut occ = 0usize;
    let mut out = Vec::new();
    for (pos, &b) in seq.iter().enumerate() {
        if b == canonical {
            if pos >= min_abs && out.len() < want {
                out.push(occ);
            }
            occ += 1;
        }
    }
    out
}

/// Builds a delta-encoded MM segment (`"{canonical}+{code},{d0,d1,...};"`) and
/// matching ML bytes for the given 0-based occurrence indices.
fn mm_segment(canonical: u8, code: char, occ_indices: &[usize], ml_seed: u8) -> (String, Vec<u8>) {
    let mut deltas = Vec::with_capacity(occ_indices.len());
    let mut prev: Option<usize> = None;
    for &occ in occ_indices {
        deltas.push(match prev {
            None => occ,
            Some(p) => occ - p - 1,
        });
        prev = Some(occ);
    }
    let mm = format!(
        "{}+{},{};",
        canonical as char,
        code,
        deltas
            .iter()
            .map(|d| d.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );
    // Distinct ML values make record/payload mismatches observable.
    let ml: Vec<u8> = (0..occ_indices.len())
        .map(|k| {
            let k = u8::try_from(k).expect("Occurrence index fits in a byte");
            100u8.wrapping_add(ml_seed).wrapping_add(k.wrapping_mul(7))
        })
        .collect();
    (mm, ml)
}

/// Writes `n` uBAM records named `r0..r{n-1}`: each `INFER_MM_ML_ADAPTER`
/// followed by a per-read splitmix64 tail, carrying a `C+m` MM/ML tag anchored
/// at occurrence positions `>= INFER_MM_ML_MOD_MIN_ABS`.
fn write_infer_mm_ml_fixture(path: &Path, n: usize) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();
    for i in 0..n {
        let mut seq = INFER_MM_ML_ADAPTER.to_vec();
        seq.extend(splitmix_tail_infer(i, INFER_MM_ML_TAIL_LEN));

        let occ = occurrence_indices_after(&seq, b'C', INFER_MM_ML_MOD_MIN_ABS, 3);
        assert!(
            !occ.is_empty(),
            "Read {i}: the fixture has at least one C past position \
             {INFER_MM_ML_MOD_MIN_ABS} for a meaningful mod tag"
        );
        let ml_seed = u8::try_from(i).expect("Fixture index fits in a byte");
        let (mm, ml) = mm_segment(b'C', 'm', &occ, ml_seed);

        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(format!("r{i}").into_bytes().into());
        let seq_len = seq.len();
        *rec.sequence_mut() = seq.into();
        *rec.quality_scores_mut() = vec![40u8; seq_len].into();
        let data = rec.data_mut();
        data.insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(mm.into_bytes().into()),
        );
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(ml)),
        );
        data.insert(
            Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
            Value::Int32(seq_len as i32),
        );
        w.write_alignment_record(&header, &rec).unwrap();
    }
    w.try_finish().unwrap();
}

/// Decodes every record's (name, SEQ) from a BAM via noodles, independently of
/// the htslib oracle used for mods; only the raw bases are needed to check the
/// suffix and cut-length property.
fn read_bam_seqs(path: &Path) -> HashMap<String, Vec<u8>> {
    let mut r = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let hdr = r.read_header().unwrap();
    let mut out = HashMap::new();
    let mut buf = RecordBuf::default();
    while r.read_record_buf(&hdr, &mut buf).unwrap() != 0 {
        let name = String::from_utf8(buf.name().unwrap().to_vec()).unwrap();
        out.insert(name, buf.sequence().as_ref().to_vec());
    }
    out
}

#[test]
fn infer_on_ubam_preserves_mm_ml() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    let n = 150; // >= MIN_SAMPLE_FOR_DETECTION (100)
    write_infer_mm_ml_fixture(&input, n);

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--adapter-infer",
            "-t",
            "1",
        ])
        .assert()
        .success();

    // Output BAM parses and retains every long read.
    let in_seqs = read_bam_seqs(&input);
    let out_seqs = read_bam_seqs(&output);
    assert_eq!(
        out_seqs.len(),
        n,
        "All {n} reads survive (long inserts, default min-length 1)"
    );

    // Each output sequence is an adapter-trimmed suffix of its input.
    for (name, out_seq) in &out_seqs {
        let orig_seq = in_seqs.get(name).expect("Matching input read");
        assert!(
            orig_seq.ends_with(out_seq.as_slice()),
            "Output read {name:?} is an exact suffix of its original"
        );
        // `checked_sub` panics with a clear message if output exceeds input,
        // instead of an overflow panic.
        let cut = orig_seq
            .len()
            .checked_sub(out_seq.len())
            .expect("Output longer than input");
        assert!(
            (15..=60).contains(&cut),
            "Read {name:?}: cut length {cut} is not adapter-shaped (planted adapter is 28 bp)"
        );
    }

    // Decoded output modifications equal the input calls retained by each
    // read's cut interval.
    let orig_mods = hts_mods_by_read(&input);
    let out_mods = hts_mods_by_read(&output);
    assert_eq!(out_mods.len(), n, "All reads decode from the output BAM");

    let mut any_nonempty = false;
    for (name, (_out_len, got)) in &out_mods {
        let (orig_len, orig) = orig_mods
            .get(name)
            .unwrap_or_else(|| panic!("Output read {name} has no matching original"));
        let out_seq_len = out_seqs[name].len();
        // `checked_sub` panics with a clear message if output exceeds input,
        // instead of an overflow panic.
        let cut = orig_len
            .checked_sub(out_seq_len)
            .expect("Output longer than input");
        let expected: Vec<_> = orig
            .iter()
            .filter(|(pos, ..)| *pos >= cut)
            .map(|&(pos, cb, mb, st, q)| (pos - cut, cb, mb, st, q))
            .collect();
        assert_eq!(
            sorted(&expected),
            sorted(got),
            "Read {name}: MM/ML survives adapter-infer trimming in register"
        );
        if !got.is_empty() {
            any_nonempty = true;
        }
    }
    assert!(
        any_nonempty,
        "At least one read retains a non-empty mod call after trimming, so the \
         comparison above is not empty against empty"
    );
}

// Filtering one sibling segment does not change the retained segment's MM/ML.

/// Twelve-base retained flank with C at indices 1, 4, 7, and 10.
const KEPT_FLANK: &[u8] = b"ACGGCGGCGGCG";
/// 16-bp interior adapter, G/T only: it cannot match the all-A/C flanks, and it
/// is the adapter `tests/adapter_cli.rs`'s naming test uses.
const SPLIT_ADAPTER: &[u8] = b"GGGGTTTTGGGGTTTT";
/// 4-bp sibling segment, below `-l 5`, filtered after the split. It contains no
/// `C`, so it cannot perturb `KEPT_FLANK`'s C-occurrence indexing even though
/// MM counts occurrences over the whole original read.
const SHORT_FLANK: &[u8] = b"TTTT";

/// Writes a one-record uBAM: `seq`, quality 40 throughout, and `KEPT_FLANK`'s
/// mod tag (`C+m,0,1,0;`, ML `[250,5,200]`, MN `seq.len()`).
fn write_mods_fixture(path: &Path, seq: &[u8]) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".to_vec().into());
    let seq_len = seq.len();
    *rec.sequence_mut() = seq.to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; seq_len].into();
    let data = rec.data_mut();
    data.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1,0;".to_vec().into()),
    );
    data.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![250, 5, 200])),
    );
    data.insert(
        Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
        Value::Int32(seq_len as i32),
    );
    w.write_alignment_record(&header, &rec).unwrap();
    w.try_finish().unwrap();
}

/// A retained segment has the same decoded modifications with or without a
/// sibling that is removed by post-trim filtering.
#[test]
fn filtered_sibling_segment_does_not_corrupt_kept_segment_mods() {
    let dir = tempfile::tempdir().unwrap();

    // Split case: `KEPT_FLANK`, the interior adapter, and a 4-bp sibling that
    // `-l 5` filters after the split.
    let mut split_seq = KEPT_FLANK.to_vec();
    split_seq.extend_from_slice(SPLIT_ADAPTER);
    split_seq.extend_from_slice(SHORT_FLANK);
    let split_in = dir.path().join("split_in.ubam");
    let split_out = dir.path().join("split_out.ubam");
    write_mods_fixture(&split_in, &split_seq);

    // Reference: the same 12-bp kept flank with the same mod tag, run with no
    // adapter config, so the short sibling and the adapter are absent.
    let solo_in = dir.path().join("solo_in.ubam");
    let solo_out = dir.path().join("solo_out.ubam");
    write_mods_fixture(&solo_in, KEPT_FLANK);

    let mut fa = tempfile::NamedTempFile::new().unwrap();
    {
        use std::io::Write as _;
        writeln!(fa, ">mid\n{}", std::str::from_utf8(SPLIT_ADAPTER).unwrap()).unwrap();
    }

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "-i",
            split_in.to_str().unwrap(),
            "-o",
            split_out.to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "--adapter-error-rate",
            "0.1",
            "--adapter-end-size",
            "1",
            "-l",
            "5",
            "-t",
            "1",
        ])
        .assert()
        .success();

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "-i",
            solo_in.to_str().unwrap(),
            "-o",
            solo_out.to_str().unwrap(),
            "-l",
            "5",
            "-t",
            "1",
        ])
        .assert()
        .success();

    // The retained record keeps its produced index after its sibling is filtered.
    let split_names = read_bam_seqs(&split_out);
    assert_eq!(
        split_names.len(),
        1,
        "The 4-bp sibling is filtered, not written"
    );
    assert!(
        split_names.contains_key("r1_segment_1"),
        "Kept segment keeps its produced index: {:?}",
        split_names.keys().collect::<Vec<_>>()
    );

    let split_mods = hts_mods(&split_out);
    let solo_mods = hts_mods(&solo_out);
    let mut a = split_mods.clone();
    let mut b = solo_mods.clone();
    a.sort();
    b.sort();
    assert_eq!(
        a, b,
        "The kept segment's MM/ML equals the same 12-bp flank run without the \
         filtered sibling: split={split_mods:?} solo={solo_mods:?}"
    );
    assert!(
        a.len() >= 3,
        "All 3 modified positions decode, so the comparison is not empty against empty: {a:?}"
    );
}

/// Fixture covering the three MM shapes the `+`/ACGT fixtures above never reach:
/// a `-` strand group, an RNA `U+` group, and an `N+` group.
///
/// Each is a case where the counting base is not the literal MM byte: `G-m`
/// counts `G` rather than its complement, `U+m` counts `T` because BAM SEQ has
/// no `U`, and `N+n` counts every base rather than literal `N`.
fn write_fixture_strands(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"read1".into());
    // AGGTACGATCGT: G at 1,2,6,10; T at 3,8,11; length 12.
    *rec.sequence_mut() = b"AGGTACGATCGT".to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; 12].into();
    let data = rec.data_mut();
    // G-m at G-occurrences 0,2 -> abs 1,6.
    // U+m at T-occurrences 0,2 -> abs 3,11.
    // N+n at any-base occurrences 4,6 -> abs 4,6.
    data.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"G-m,0,1;U+m,0,1;N+n,4,1;".to_vec().into()),
    );
    data.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![200, 201, 150, 151, 90, 91])),
    );
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(12));
    w.write_alignment_record(&header, &rec).unwrap();
}

/// Untrimmed round trip: a no-op run carries all three group shapes through byte
/// for byte. This guards the pass-through path only; an untrimmed record never
/// reaches the rewrite, so `strand_and_ambiguity_groups_match_oracle_when_trimmed`
/// covers the counting-base rules.
#[test]
fn strand_and_ambiguity_groups_survive_untrimmed_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    write_fixture_strands(&input);

    let cfg = whittle::cli::config_for_test(&input, &output, 0, 0);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).unwrap();

    let mut before = hts_mods(&input);
    let mut after = hts_mods(&output);
    before.sort();
    after.sort();
    assert_eq!(
        before, after,
        "A no-op run neither moves, drops, nor relabels any modification call"
    );

    // Guard against a trivial pass: all three groups are present.
    assert_eq!(before.len(), 6, "Fixture decodes to six calls");
    let canonical: std::collections::HashSet<char> = before.iter().map(|t| t.1).collect();
    assert!(
        canonical.contains(&'G') && canonical.contains(&'T'),
        "Expected G and T canonical bases, got {canonical:?}"
    );
}

/// The same three shapes through a real head/tail crop, compared against htslib's
/// decode of the input restricted to the surviving window.
#[test]
fn strand_and_ambiguity_groups_match_oracle_when_trimmed() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    write_fixture_strands(&input);

    let (head, tail, len) = (2usize, 2usize, 12usize);
    let cfg = whittle::cli::config_for_test(&input, &output, head, tail);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).unwrap();

    let tail_start = len - tail;
    let mut expected: Vec<_> = hts_mods(&input)
        .iter()
        .filter(|(pos, ..)| *pos >= head && *pos < tail_start)
        .map(|&(pos, cb, mb, st, q)| (pos - head, cb, mb, st, q))
        .collect();
    let mut got = hts_mods(&output);
    expected.sort();
    got.sort();
    assert_eq!(
        expected, got,
        "Trimmed calls equal the original restricted to [{head},{tail_start})"
    );
    assert!(!got.is_empty(), "Window retains at least one call");
}

/// Builds a one-read uBAM with caller-chosen flags and tag spelling.
fn write_flagged(path: &Path, flags: Flags, mm_tag: &[u8; 2], ml_tag: &[u8; 2]) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut rec = RecordBuf::default();
    *rec.flags_mut() = flags;
    *rec.name_mut() = Some(b"read1".into());
    *rec.sequence_mut() = b"AGCTACGTNACGGTACCGTA".to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; 20].into();
    let data = rec.data_mut();
    data.insert(
        Tag::new(mm_tag[0], mm_tag[1]),
        Value::String(b"C+m,0,1;".to_vec().into()),
    );
    data.insert(
        Tag::new(ml_tag[0], ml_tag[1]),
        Value::Array(Array::UInt8(vec![200, 201])),
    );
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(20));
    w.write_alignment_record(&header, &rec).unwrap();
}

/// htslib reads the pre-spec lowercase `Mm`/`Ml` (emitted by early guppy and
/// megalodon). whittle rewrites only `MM`/`ML`, so a record carrying the
/// lowercase spelling is refused rather than copied onto a trimmed sequence with
/// every call pointing at the wrong base.
#[test]
fn legacy_lowercase_mod_tags_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("legacy.ubam");
    let output = dir.path().join("out.ubam");
    write_flagged(&input, Flags::UNMAPPED, b"Mm", b"Ml");

    let cfg = whittle::cli::config_for_test(&input, &output, 3, 4);
    let mut h = whittle::obs::ProgressHandle::disabled();
    let err = whittle::run(cfg, &mut h).expect_err("A legacy-tag record is refused");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("legacy") && msg.contains("Mm"),
        "Error names the legacy tag, got: {msg}"
    );
}

/// A record flagged reverse-complemented stores SEQ reverse-complemented and
/// htslib decodes its MM right to left with complemented bases. whittle trims in
/// read orientation, so it would crop the wrong ends and relocate every call.
/// Basecallers do not emit `0x4|0x10`, but `samtools view -f 4` of an aligned
/// file preserves it.
#[test]
fn unmapped_but_reverse_complemented_records_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rev.ubam");
    let output = dir.path().join("out.ubam");
    write_flagged(
        &input,
        Flags::UNMAPPED | Flags::REVERSE_COMPLEMENTED,
        b"MM",
        b"ML",
    );

    let cfg = whittle::cli::config_for_test(&input, &output, 2, 2);
    let mut h = whittle::obs::ProgressHandle::disabled();
    let err = whittle::run(cfg, &mut h).expect_err("A reverse-complemented record is refused");
    assert!(
        format!("{err:#}").contains("reverse-complemented"),
        "Error names the flag, got: {err:#}"
    );
}

/// The same record without either problem still round-trips, so the two
/// refusals above are specific.
#[test]
fn plain_unmapped_record_with_modern_tags_still_runs() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("ok.ubam");
    let output = dir.path().join("out.ubam");
    write_flagged(&input, Flags::UNMAPPED, b"MM", b"ML");

    let cfg = whittle::cli::config_for_test(&input, &output, 3, 4);
    let mut h = whittle::obs::ProgressHandle::disabled();
    whittle::run(cfg, &mut h).expect("A modern uBAM record trims normally");

    let (head, tail, len) = (3usize, 4usize, 20usize);
    let mut expected: Vec<_> = hts_mods(&input)
        .iter()
        .filter(|(pos, ..)| *pos >= head && *pos < len - tail)
        .map(|&(pos, cb, mb, st, q)| (pos - head, cb, mb, st, q))
        .collect();
    let mut got = hts_mods(&output);
    expected.sort();
    got.sort();
    assert_eq!(expected, got, "Surviving calls match htslib");
}
