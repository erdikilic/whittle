//! Conversion of legacy tag forms: an integer mean quality `qs` and HiFi
//! kinetics that index the whole CCS read. The `Mm`/`Ml` modification tags
//! are tested against htslib in `bam_mods_oracle.rs`.

use std::path::Path;

use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

/// Base qualities of the fixture read: a low-quality start and end.
const QUALS: [u8; 20] = [
    5, 8, 12, 30, 30, 30, 35, 35, 35, 40, 40, 40, 35, 35, 30, 30, 20, 12, 8, 5,
];

/// Builds a one-read uBAM whose mean quality is an integer `qs`, without
/// `qe`.
fn write_integer_qs(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"0a1b2c3d-ont-read".into());
    *rec.sequence_mut() = b"AGCTACGTAACGGTACCGTA".to_vec().into();
    *rec.quality_scores_mut() = QUALS.to_vec().into();
    rec.data_mut()
        .insert(Tag::new(b'q', b's'), Value::Int32(13));
    w.write_alignment_record(&header, &rec).unwrap();
}

/// The `qs` value of every record of the uBAM at `path`.
fn qs_values(path: &Path) -> Vec<Value> {
    let mut reader = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    reader
        .record_bufs(&header)
        .map(|rec| {
            rec.unwrap()
                .data()
                .get(&Tag::new(b'q', b's'))
                .cloned()
                .expect("qs is kept")
        })
        .collect()
}

/// An integer ONT `qs` is written as a float: unchanged in value when nothing
/// is trimmed, and as the mean quality of the kept bases after a trim, on the
/// single-threaded and the parallel path.
#[test]
fn integer_ont_qs_is_written_as_a_float() {
    for (head, tail) in [(0, 0), (3, 4)] {
        for threads in [1, 2] {
            let dir = tempfile::tempdir().unwrap();
            let input = dir.path().join("in.ubam");
            let output = dir.path().join("out.ubam");
            write_integer_qs(&input);
            let cfg = whittle::cli::config_for_test_threads(&input, &output, head, tail, threads);
            let mut h = whittle::obs::ProgressHandle::disabled();
            whittle::run(cfg, &mut h).unwrap();

            let expected = if head + tail == 0 {
                13.0
            } else {
                whittle::qual::mean_prob_q(&QUALS[head..QUALS.len() - tail]) as f32
            };
            match qs_values(&output).as_slice() {
                [Value::Float(q)] => assert!(
                    (q - expected).abs() < 1e-4,
                    "qs {q}, expected {expected} (head {head}, tail {tail}, {threads} threads)"
                ),
                other => panic!("Unexpected qs: {other:?}"),
            }
        }
    }
}

/// Length of the whole CCS read: 17 bases before the 20-base read and 8 after
/// it.
const WHOLE: usize = 45;

/// Builds a one-read uBAM of a HiFi read clipped to [17, 37) of a 45-base CCS
/// read (`qs`/`qe`), whose kinetic arrays index the whole CCS read.
fn write_unclipped_kinetics(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"m64043_220516_150157/1/ccs".into());
    *rec.sequence_mut() = b"AGCTACGTAACGGTACCGTA".to_vec().into();
    *rec.quality_scores_mut() = vec![40; 20].into();
    let data = rec.data_mut();
    data.insert(Tag::new(b'q', b's'), Value::Int32(17));
    data.insert(Tag::new(b'q', b'e'), Value::Int32(37));
    let whole: Vec<u8> = (0..WHOLE as u8).collect();
    for t in [b"fi", b"fp", b"ri", b"rp"] {
        data.insert(
            Tag::new(t[0], t[1]),
            Value::Array(Array::UInt8(whole.clone())),
        );
    }
    w.write_alignment_record(&header, &rec).unwrap();
}

/// The kinetic arrays of the single record of the uBAM at `path`.
fn kinetics(path: &Path) -> Vec<Vec<u8>> {
    let mut reader = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let rec = reader.record_bufs(&header).next().unwrap().unwrap();
    [b"fi", b"fp", b"ri", b"rp"]
        .iter()
        .map(|t| match rec.data().get(&Tag::new(t[0], t[1])) {
            Some(Value::Array(Array::UInt8(v))) => v.clone(),
            other => panic!("Unexpected {t:?}: {other:?}"),
        })
        .collect()
}

/// Kinetic arrays that index the whole CCS read are cut to the read by its
/// query coordinates and trimmed with it: forward arrays from the read's first
/// base, reverse arrays from its last, on both paths.
#[test]
fn unclipped_lima_kinetics_follow_the_read() {
    for (head, tail) in [(0, 0), (3, 4)] {
        for threads in [1, 2] {
            let dir = tempfile::tempdir().unwrap();
            let input = dir.path().join("in.ubam");
            let output = dir.path().join("out.ubam");
            write_unclipped_kinetics(&input);
            let cfg = whittle::cli::config_for_test_threads(&input, &output, head, tail, threads);
            let mut h = whittle::obs::ProgressHandle::disabled();
            whittle::run(cfg, &mut h).unwrap();

            let (start, end) = (17 + head, 37 - tail);
            let forward: Vec<u8> = (start as u8..end as u8).collect();
            let reverse: Vec<u8> = ((WHOLE - end) as u8..(WHOLE - start) as u8).collect();
            assert_eq!(
                kinetics(&output),
                vec![forward.clone(), forward, reverse.clone(), reverse],
                "head {head}, tail {tail}, {threads} threads"
            );
        }
    }
}

/// The run reports cut arrays as a legacy form, not as malformed tags.
#[test]
fn unclipped_lima_kinetics_are_counted_as_a_legacy_form() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    let json = dir.path().join("summary.json");
    write_unclipped_kinetics(&input);
    assert_cmd::Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .args(["--summary-json", json.to_str().unwrap(), "-H", "3"])
        .assert()
        .success();
    let summary: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(summary["warnings"]["legacy_tags_reads"], 1);
    assert_eq!(summary["warnings"]["malformed_tag_reads"], 0);
}
