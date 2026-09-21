//! `--tag-filter`: read selection by aux-tag expression on BAM, on tagged
//! FASTQ, and its refusal on plain FASTQ.

use std::path::Path;

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::{self as sam};
use predicates::prelude::*;

/// Three reads with ONT-style whole-read metadata: end reason, duplex state
/// and a dorado quality score stored as a float on two reads and as an
/// integer on the third.
fn write_fixture(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();
    let reads: [(&[u8], &[u8], i8, Value); 3] = [
        (b"kept", b"signal_positive", 1, Value::Float(12.5)),
        (
            b"unblocked",
            b"data_service_unblock_mux_change",
            0,
            Value::Float(9.0),
        ),
        (b"simplex", b"signal_positive", 0, Value::Int32(15)),
    ];
    for (name, er, dx, qs) in reads {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::UNMAPPED;
        *r.name_mut() = Some(name.into());
        *r.sequence_mut() = b"ACGTACGTACGTACGTACGT".to_vec().into();
        *r.quality_scores_mut() = vec![30; 20].into();
        let data = r.data_mut();
        data.insert(Tag::new(b'e', b'r'), Value::String(er.into()));
        data.insert(Tag::new(b'd', b'x'), Value::Int8(dx));
        data.insert(Tag::new(b'q', b's'), qs);
        w.write_alignment_record(&header, &r).unwrap();
    }
    w.try_finish().unwrap();
}

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

/// Reads every record name of a uBAM.
fn names(path: &Path) -> Vec<String> {
    let mut rdr = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let hdr = rdr.read_header().unwrap();
    let mut buf = RecordBuf::default();
    let mut out = Vec::new();
    while rdr.read_record_buf(&hdr, &mut buf).unwrap() > 0 {
        out.push(String::from_utf8_lossy(buf.name().unwrap().as_ref()).into_owned());
    }
    out
}

fn summary(path: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// The filter runs on the raw pass-through path (no trim flags) and on the
/// rebuilt path (`-H 1`), and the summary counts the rejected reads as
/// input that was tag-filtered.
#[test]
fn tag_filter_selects_reads_on_both_bam_paths() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.bam");
    write_fixture(&input);
    for (label, extra) in [("raw", vec![]), ("rebuilt", vec!["-H", "1"])] {
        let out = dir.path().join(format!("{label}.bam"));
        let json = dir.path().join(format!("{label}.json"));
        whittle()
            .args(["--tag-filter", "[er]!=\"data_service_unblock_mux_change\""])
            .args(["--tag-filter", "[qs]>=10"])
            .args(&extra)
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .args(["--summary-json", json.to_str().unwrap()])
            .assert()
            .success()
            .stderr(predicate::str::contains("Tag filter: [qs]>=10"))
            .stderr(predicate::str::contains(
                "Tag filtered: 1 input reads did not satisfy --tag-filter",
            ));
        assert_eq!(names(&out), ["kept", "simplex"], "{label}");
        let v = summary(&json);
        assert_eq!(v["reads"]["input"], 3, "{label}");
        assert_eq!(v["reads"]["tag_filtered"], 1, "{label}");
        assert_eq!(v["reads"]["with_output"], 2, "{label}");
        assert_eq!(v["params"]["tag_filter"][1], "[qs]>=10", "{label}");
    }
}

/// Integer and float storage of the same tag compare alike, a missing tag
/// fails every comparison, and presence tests work on their own.
#[test]
fn tag_filter_semantics_on_bam() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.bam");
    write_fixture(&input);
    for (expr, expected) in [
        ("[qs]>=12", vec!["kept", "simplex"]),
        ("[dx]==1", vec!["kept"]),
        ("[dx]==1 || ([dx]==0 && [qs]>=15)", vec!["kept", "simplex"]),
        ("[BC]!=\"barcode01\"", vec![]),
        ("![BC]", vec!["kept", "unblocked", "simplex"]),
        ("!([qs]>=12)", vec!["unblocked"]),
        ("exists([er])", vec!["kept", "unblocked", "simplex"]),
    ] {
        let out = dir.path().join("out.bam");
        whittle()
            .args(["--tag-filter", expr])
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .assert()
            .success();
        assert_eq!(names(&out), expected, "{expr}");
    }
}

/// A comparison between a number and a string is an error naming the read
/// and the tag, not a silent drop of every read.
#[test]
fn tag_filter_type_mismatch_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.bam");
    write_fixture(&input);
    whittle()
        .args(["--tag-filter", "[qs]==\"10\""])
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            dir.path().join("out.bam").to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("read kept"))
        .stderr(predicate::str::contains("[qs]"));
    whittle()
        .args(["--tag-filter", "[qs] >="])
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            dir.path().join("out.bam").to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--tag-filter \"[qs] >=\""));
}

/// Tagged FASTQ carries the same tags in its headers and is filtered the same
/// way; a plain FASTQ has no tags to read and the run is refused.
#[test]
fn tag_filter_on_tagged_and_plain_fastq() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.bam");
    write_fixture(&input);
    let tagged = dir.path().join("tagged.fastq");
    whittle()
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            tagged.to_str().unwrap(),
            "--quiet",
        ])
        .assert()
        .success();
    let out = dir.path().join("kept.fastq");
    whittle()
        .args(["--tag-filter", "[dx]==1 || [qs]>=15", "-H", "2"])
        .args(["-i", tagged.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .assert()
        .success()
        .stderr(predicate::str::contains("Tag filtered: 1 input reads"));
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(
        text.starts_with("@kept\t") && text.contains("\n@simplex\t"),
        "{text}"
    );
    assert!(!text.contains("unblocked"), "{text}");

    let plain = dir.path().join("plain.fastq");
    std::fs::write(&plain, "@r1 runid=x\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    whittle()
        .args(["--tag-filter", "[dx]==1"])
        .args([
            "-i",
            plain.to_str().unwrap(),
            "-o",
            dir.path().join("p.fastq").to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "--tag-filter reads aux tags and requires BAM or tagged FASTQ input",
        ));
}
