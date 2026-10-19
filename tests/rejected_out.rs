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

/// Phred+33 text for `len` bases at quality `q`.
fn qual(q: u8, len: usize) -> String {
    char::from(q + 33).to_string().repeat(len)
}

/// `--max-expected-errors` judges each quality-split piece on its own and an
/// unsplit read as a whole; a rejected segment carries the `expected_errors`
/// reason and the summary counts it with the resolved bound.
#[test]
fn expected_errors_reject_each_segment_on_its_own() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fastq");
    // r1 splits into 100 Q30 bases (0.1 expected errors) and 100 Q15 bases
    // (3.2); r2 carries 0.5 and r3 2.0 expected errors.
    let r1 = [qual(30, 100), qual(3, 40), qual(15, 100)].concat();
    std::fs::write(
        &input,
        format!(
            "@r1\n{}\n+\n{r1}\n@r2\n{}\n+\n{}\n@r3\n{}\n+\n{}\n",
            "A".repeat(240),
            "C".repeat(50),
            qual(20, 50),
            "G".repeat(200),
            qual(20, 200)
        ),
    )
    .unwrap();
    let out = dir.path().join("out.fastq");
    let rej = dir.path().join("rejected.fastq");
    let json = dir.path().join("summary.json");
    whittle()
        .args(["--quality-trim", "runs", "--quality-cutoff", "10"])
        .args(["--min-low-quality-run", "10", "--max-expected-errors", "1"])
        .args(["-t", "1", "--quiet"])
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--rejected-output", rej.to_str().unwrap()])
        .args(["--summary-json", json.to_str().unwrap()])
        .assert()
        .success();
    assert_eq!(
        fastq_headers(&std::fs::read_to_string(&out).unwrap()),
        ["r1_segment_1", "r2"]
    );
    let rejected = std::fs::read_to_string(&rej).unwrap();
    assert_eq!(
        fastq_headers(&rejected),
        [
            "r1_segment_2\twr:Z:expected_errors",
            "r3\twr:Z:expected_errors"
        ],
        "{rejected}"
    );
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(v["params"]["max_expected_errors"], 1.0);
    assert_eq!(v["segments_dropped"]["expected_errors"], 2);
    assert_eq!(v["reads"]["with_output"], 2);
    assert_eq!(v["reads"]["all_filtered"], 1);
}

/// A failed write to the rejected output names the file and the system
/// error, at one thread and at several.
#[cfg(target_os = "linux")]
#[test]
fn rejected_output_write_failure_names_the_file_and_the_cause() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fastq");
    let read = format!("{}\n+\n{}\n", "ACGT".repeat(100), qual(30, 400));
    let text: String = (0..3000).map(|i| format!("@r{i}\n{read}")).collect();
    std::fs::write(&input, text).unwrap();
    for threads in ["1", "3"] {
        whittle()
            .args(["-l", "100000", "-t", threads])
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
            .args(["--rejected-output", "/dev/full"])
            .assert()
            .failure()
            .stderr(predicate::str::contains("rejected output /dev/full"))
            .stderr(predicate::str::contains("(os error 28)"));
    }
}

/// The `@HD` sort order of a BAM file.
fn sort_order(path: &Path) -> Option<String> {
    use noodles_sam::header::record::value::map::header::tag;
    let mut rdr = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let header = rdr.read_header().unwrap();
    header
        .header()
        .and_then(|hd| hd.other_fields().get(&tag::SORT_ORDER))
        .map(|v| String::from_utf8_lossy(v.as_ref()).into_owned())
}

/// Rejected records are written as they arrive, so under several threads
/// the rejected BAM header is unsorted even with `--preserve-order`, which
/// keeps the input order of the main output only. One thread keeps both.
#[test]
fn rejected_bam_header_is_unsorted_under_threads() {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::header::tag;

    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("sorted.bam");
    let hd = Map::<noodles_sam::header::record::value::map::Header>::builder()
        .insert(tag::SORT_ORDER, "queryname")
        .build()
        .unwrap();
    let header = sam::Header::builder().set_header(hd).build();
    let mut w = bam::io::Writer::new(std::fs::File::create(&input).unwrap());
    w.write_header(&header).unwrap();
    for i in 0..40 {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::UNMAPPED;
        *r.name_mut() = Some(format!("r{i:02}").into_bytes().into());
        let len = if i % 2 == 0 { 30 } else { 10 };
        *r.sequence_mut() = b"ACGT".repeat(10)[..len].to_vec().into();
        *r.quality_scores_mut() = vec![30; len].into();
        w.write_alignment_record(&header, &r).unwrap();
    }
    w.try_finish().unwrap();
    for (threads, rejected_order) in [("1", "queryname"), ("3", "unsorted")] {
        let out = dir.path().join(format!("out{threads}.bam"));
        let rej = dir.path().join(format!("rejected{threads}.bam"));
        whittle()
            .args(["-t", threads, "--preserve-order", "-l", "20"])
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .args(["--rejected-output", rej.to_str().unwrap()])
            .assert()
            .success();
        assert_eq!(sort_order(&out).as_deref(), Some("queryname"));
        assert_eq!(
            sort_order(&rej).as_deref(),
            Some(rejected_order),
            "-t {threads}"
        );
        assert_eq!(records(&rej).len(), 20);
    }
}

/// One record's modification tags: name, `MM`, `ML` length and `MN`.
type ModTags = (String, Option<String>, Option<usize>, Option<i64>);

/// The modification tags of every record of a uBAM, sorted by name.
fn mod_tags(path: &Path) -> Vec<ModTags> {
    use noodles_sam::alignment::record_buf::data::field::value::Array;
    let mut rdr = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let hdr = rdr.read_header().unwrap();
    let mut buf = RecordBuf::default();
    let mut out = Vec::new();
    while rdr.read_record_buf(&hdr, &mut buf).unwrap() > 0 {
        let data = buf.data();
        let mm = match data.get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
            _ => None,
        };
        let ml = match data.get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) => Some(v.len()),
            _ => None,
        };
        let mn = data
            .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
            .and_then(|v| v.as_int());
        out.push((
            String::from_utf8_lossy(buf.name().unwrap().as_ref()).into_owned(),
            mm,
            ml,
            mn,
        ));
    }
    out.sort();
    out
}

/// A rejected whole read gets the modification tags a trimmed record gets:
/// a malformed block is removed and a missing `MN` is added, on the raw
/// full-window path as on the decoded path that `--remove-tag` selects.
#[test]
fn rejected_reads_rewrite_their_modification_block_on_every_path() {
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("mods.bam");
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(&input).unwrap());
    w.write_header(&header).unwrap();
    for (name, mm) in [("malformed", b"C+m,5;".as_slice()), ("no_mn", b"C+m,0;")] {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::UNMAPPED;
        *r.name_mut() = Some(name.into());
        *r.sequence_mut() = b"ACGT".to_vec().into();
        *r.quality_scores_mut() = vec![30; 4].into();
        let data = r.data_mut();
        data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.into()));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(vec![200])),
        );
        w.write_alignment_record(&header, &r).unwrap();
    }
    w.try_finish().unwrap();
    let rejected = |extra: &[&str], name: &str| {
        let rej = dir.path().join(name);
        whittle()
            .args(["-l", "10", "-t", "1"])
            .args(extra)
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", dir.path().join("out.bam").to_str().unwrap()])
            .args(["--rejected-output", rej.to_str().unwrap()])
            .assert()
            .success();
        mod_tags(&rej)
    };
    let raw = rejected(&[], "raw.bam");
    assert_eq!(
        raw,
        [
            ("malformed".to_string(), None, None, None),
            (
                "no_mn".to_string(),
                Some("C+m,0;".to_string()),
                Some(1),
                Some(4)
            ),
        ]
    );
    assert_eq!(rejected(&["--remove-tag", "zz"], "decoded.bam"), raw);
}

/// A rejected FASTQ record carries one `wr` tag, the new reason, in place of
/// a `wr` its input header already holds, as a rejected BAM record does.
#[test]
fn rejected_fastq_replaces_an_input_reason_tag() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("rejected_before.fastq");
    std::fs::write(
        &input,
        "@short\twr:Z:too_short\nACGTACGT\n+\nIIIIIIII\n\
         @low\twr:Z:old\tqs:f:5\nACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let rej = dir.path().join("rejected_again.fastq");
    whittle()
        .args(["-l", "10", "-t", "1", "--tag-filter", "![qs] || [qs]>10"])
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args(["--rejected-output", rej.to_str().unwrap()])
        .assert()
        .success();
    let text = std::fs::read_to_string(&rej).unwrap();
    assert_eq!(
        fastq_headers(&text),
        ["low\tqs:f:5\twr:Z:tag_filter", "short\twr:Z:too_short"],
        "{text}"
    );
}

/// Without `-c`, the rejected file takes the default level of its own format,
/// 4 for `.fastq.gz`, whatever the main output's format; `-c` sets both.
#[test]
fn rejected_output_takes_the_default_level_of_its_own_format() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fastq");
    let mut state = 7u64;
    let text: String = (0..2000)
        .map(|i| {
            let seq: String = (0..200)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    char::from(b"ACGT"[(state >> 62) as usize])
                })
                .collect();
            format!("@r{i}\n{seq}\n+\n{}\n", qual(30, 200))
        })
        .collect();
    std::fs::write(&input, text).unwrap();
    let rejected = |name: &str, level: Option<&str>| {
        let rej = dir.path().join(format!("{name}.fastq.gz"));
        let mut cmd = whittle();
        cmd.args(["-l", "1000", "-t", "1"])
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
            .args(["--rejected-output", rej.to_str().unwrap()]);
        if let Some(level) = level {
            cmd.args(["-c", level]);
        }
        cmd.assert().success();
        std::fs::read(&rej).unwrap()
    };
    let default = rejected("default", None);
    assert_eq!(default, rejected("four", Some("4")));
    assert_ne!(default, rejected("six", Some("6")));
}
