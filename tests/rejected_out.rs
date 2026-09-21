//! `--rejected-output`: every read or segment that does not reach the output is
//! written to a second file with a `wr:Z` reason tag.

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

/// Three reads: a 20-base duplex read, an adaptive-sampling reject and a
/// 12-base simplex read.
fn write_fixture(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();
    let reads: [(&[u8], &[u8], i8, usize); 3] = [
        (b"kept", b"signal_positive", 1, 20),
        (b"unblocked", b"data_service_unblock_mux_change", 0, 20),
        (b"simplex", b"signal_positive", 0, 12),
    ];
    for (name, er, dx, len) in reads {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::UNMAPPED;
        *r.name_mut() = Some(name.into());
        *r.sequence_mut() = b"ACGTACGTACGTACGTACGT"[..len].to_vec().into();
        *r.quality_scores_mut() = vec![30; len].into();
        let data = r.data_mut();
        data.insert(Tag::new(b'e', b'r'), Value::String(er.into()));
        data.insert(Tag::new(b'd', b'x'), Value::Int8(dx));
        w.write_alignment_record(&header, &r).unwrap();
    }
    w.try_finish().unwrap();
}

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

/// Every record of a uBAM as (name, sequence, wr tag).
fn records(path: &Path) -> Vec<(String, String, Option<String>)> {
    let mut rdr = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let hdr = rdr.read_header().unwrap();
    let mut buf = RecordBuf::default();
    let mut out = Vec::new();
    while rdr.read_record_buf(&hdr, &mut buf).unwrap() > 0 {
        let wr = match buf.data().get(&Tag::new(b'w', b'r')) {
            Some(Value::String(s)) => Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
            _ => None,
        };
        out.push((
            String::from_utf8_lossy(buf.name().unwrap().as_ref()).into_owned(),
            String::from_utf8_lossy(buf.sequence().as_ref()).into_owned(),
            wr,
        ));
    }
    out.sort();
    out
}

fn fastq_headers(text: &str) -> Vec<&str> {
    let mut v: Vec<&str> = text
        .lines()
        .step_by(4)
        .map(|h| h.trim_start_matches('@'))
        .collect();
    v.sort();
    v
}

/// On the rebuilt BAM path a tag-filtered read arrives untrimmed and a
/// too-short segment arrives trimmed; the main output holds the rest, and the
/// summary counts are those of the run without the flag.
#[test]
fn rejected_bam_holds_filtered_reads_and_dropped_segments() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.bam");
    write_fixture(&input);
    for threads in ["1", "4"] {
        let out = dir.path().join("out.bam");
        let rej = dir.path().join("rejected.bam");
        let json = dir.path().join("run.json");
        whittle()
            .args(["--tag-filter", "[er]!=\"data_service_unblock_mux_change\""])
            .args(["-H", "8", "-l", "5", "-t", threads])
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .args(["--rejected-output", rej.to_str().unwrap()])
            .args(["--summary-json", json.to_str().unwrap()])
            .assert()
            .success()
            .stderr(predicate::str::contains("Rejected: "));
        assert_eq!(
            records(&out),
            [("kept".to_string(), "ACGTACGTACGT".to_string(), None)],
            "-t {threads}"
        );
        assert_eq!(
            records(&rej),
            [
                (
                    "simplex".to_string(),
                    "ACGT".to_string(),
                    Some("too_short".to_string())
                ),
                (
                    "unblocked".to_string(),
                    "ACGTACGTACGTACGTACGT".to_string(),
                    Some("tag_filter".to_string())
                ),
            ],
            "-t {threads}"
        );
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
        assert_eq!(v["reads"]["input"], 3);
        assert_eq!(v["reads"]["tag_filtered"], 1);
        assert_eq!(v["reads"]["all_filtered"], 1);
        assert_eq!(v["reads"]["with_output"], 1);
        assert_eq!(v["params"]["rejected_output"], rej.to_str().unwrap());
    }
}

/// The raw pass-through path (no trim flags) rejects through the same file.
#[test]
fn rejected_bam_on_the_raw_path() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.bam");
    write_fixture(&input);
    let out = dir.path().join("out.bam");
    let rej = dir.path().join("rejected.bam");
    whittle()
        .args(["-l", "15", "-t", "2"])
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--rejected-output", rej.to_str().unwrap()])
        .assert()
        .success();
    assert_eq!(records(&out).len(), 2);
    assert_eq!(
        records(&rej),
        [(
            "simplex".to_string(),
            "ACGTACGTACGT".to_string(),
            Some("too_short".to_string())
        )]
    );
}

/// With FASTQ output the rejected file is FASTQ too, and the reason is a
/// header tag.
#[test]
fn rejected_fastq_from_bam_and_from_fastq_input() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.bam");
    write_fixture(&input);
    let out = dir.path().join("out.fastq");
    let rej = dir.path().join("rejected.fastq");
    whittle()
        .args(["--tag-filter", "[dx]==1 || [er]==\"signal_positive\""])
        .args(["-H", "8", "-l", "5"])
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--rejected-output", rej.to_str().unwrap()])
        .assert()
        .success();
    let text = std::fs::read_to_string(&rej).unwrap();
    let headers = fastq_headers(&text);
    assert_eq!(headers.len(), 2, "{text}");
    assert!(
        headers[0].starts_with("simplex") && headers[0].ends_with("\twr:Z:too_short"),
        "{text}"
    );
    assert!(
        headers[1].starts_with("unblocked") && headers[1].ends_with("\twr:Z:tag_filter"),
        "{text}"
    );
    assert!(
        text.contains("\nACGT\n+\n"),
        "the rejected segment is the trimmed one: {text}"
    );

    // Plain FASTQ input: a too-short read and a read cropped to nothing, into
    // a gzip rejected file.
    let plain = dir.path().join("plain.fastq");
    std::fs::write(
        &plain,
        "@a\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n@b\nACGTAC\n+\nIIIIII\n@c\nACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let out = dir.path().join("plain_out.fastq");
    let rej = dir.path().join("plain_rejected.fastq.gz");
    whittle()
        .args(["-H", "10", "-l", "12", "-t", "3"])
        .args(["-i", plain.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--rejected-output", rej.to_str().unwrap()])
        .assert()
        .success();
    let mut text = String::new();
    std::io::Read::read_to_string(
        &mut flate2::read::MultiGzDecoder::new(std::fs::File::open(&rej).unwrap()),
        &mut text,
    )
    .unwrap();
    assert_eq!(
        fastq_headers(&text),
        ["a\twr:Z:too_short", "b\twr:Z:trimmed_to_nothing"],
        "{text}"
    );
    assert_eq!(
        fastq_headers(&std::fs::read_to_string(&out).unwrap()),
        ["c"]
    );
}

/// The rejected file must share the output's format family.
#[test]
fn rejected_output_format_family_must_match_the_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.bam");
    write_fixture(&input);
    whittle()
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            dir.path().join("out.bam").to_str().unwrap(),
        ])
        .args([
            "--rejected-output",
            dir.path().join("rejected.fastq").to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--rejected-output"))
        .stderr(predicate::str::contains("format family"));
}
