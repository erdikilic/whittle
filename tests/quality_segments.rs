//! End-to-end `--quality-trim segments` over the compiled binary: segment
//! naming, rejected output and summary counters on FASTQ input, and the
//! modification calls of each kept segment on uBAM and tagged FASTQ input.

use std::path::Path;

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

/// Phred+33 text for `len` bases at quality `q`.
fn qual(q: u8, len: usize) -> String {
    char::from(q + 33).to_string().repeat(len)
}

/// A read of three Q30 regions of 100, 60 and 100 bases separated by 40
/// bases at Q3, and a read with no base above Q10.
fn fastq_reads() -> String {
    let seq = "ACGT".repeat(85);
    let q = [
        qual(30, 100),
        qual(3, 40),
        qual(30, 60),
        qual(3, 40),
        qual(30, 100),
    ]
    .concat();
    format!(
        "@r1\n{seq}\n+\n{q}\n@r2\n{}\n+\n{}\n",
        "ACGT".repeat(25),
        qual(8, 100)
    )
}

/// Each Q30 region becomes a segment named by its position among all three,
/// the one below `-l` is rejected as too short, and the read without a
/// segment is trimmed to nothing; the summary counts the same.
#[test]
fn fastq_segments_are_named_filtered_and_counted_like_split_pieces() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fastq");
    let out = dir.path().join("out.fastq");
    let rejected = dir.path().join("rejected.fastq");
    let json = dir.path().join("summary.json");
    std::fs::write(&input, fastq_reads()).unwrap();
    whittle()
        .args([
            "--quality-trim",
            "segments",
            "--quality-cutoff",
            "10",
            "-l",
            "80",
            "--quiet",
        ])
        .args(["-t", "1"])
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--rejected-output", rejected.to_str().unwrap()])
        .args(["--summary-json", json.to_str().unwrap()])
        .assert()
        .success();

    let text = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 8, "{text}");
    assert_eq!(lines[0], "@r1_segment_1");
    assert_eq!(lines[1], &"ACGT".repeat(85)[..100]);
    assert_eq!(lines[4], "@r1_segment_3");
    assert_eq!(lines[5], &"ACGT".repeat(85)[240..340]);
    assert_eq!(lines[7], qual(30, 100));

    let rejected = std::fs::read_to_string(&rejected).unwrap();
    let mut headers: Vec<&str> = rejected.lines().step_by(4).collect();
    headers.sort_unstable();
    assert_eq!(
        headers,
        [
            "@r1_segment_2\twr:Z:too_short",
            "@r2\twr:Z:trimmed_to_nothing"
        ],
        "{rejected}"
    );

    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    assert_eq!(v["params"]["quality_trim"]["method"], "segments");
    assert_eq!(v["params"]["quality_trim"]["cutoff"], 10);
    assert_eq!(v["reads"]["input"], 2);
    assert_eq!(v["reads"]["output"], 2);
    assert_eq!(v["reads"]["with_output"], 1);
    assert_eq!(v["reads"]["trimmed_to_nothing"], 1);
    assert_eq!(v["segments_dropped"]["too_short"], 1);
}

/// A 160-base uBAM record: two 60-base `CA` repeats at Q30 around 40 bases of
/// `T` at Q3. Its C bases lie at the even positions 0 to 58 and 100 to 158;
/// `C+m,1,30;` marks C occurrences 1 and 32, at positions 2 and 104.
fn write_mod_fixture(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();
    let seq = ["CA".repeat(30), "T".repeat(40), "CA".repeat(30)].concat();
    let quals = [vec![30u8; 60], vec![3u8; 40], vec![30u8; 60]].concat();
    let mut r = RecordBuf::default();
    *r.flags_mut() = Flags::UNMAPPED;
    *r.name_mut() = Some(b"r1".into());
    *r.sequence_mut() = seq.into_bytes().into();
    *r.quality_scores_mut() = quals.into();
    let d = r.data_mut();
    d.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,1,30;".to_vec().into()),
    );
    d.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![50, 250])),
    );
    d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(160));
    w.write_alignment_record(&header, &r).unwrap();
    w.try_finish().unwrap();
}

/// A record's name, `MM`, `ML` and `MN`.
type ModRecord = (String, Vec<u8>, Vec<u8>, Option<i64>);

/// Every record of a uBAM as its name and modification tags.
fn mod_records(path: &Path) -> Vec<ModRecord> {
    let mut rdr = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let hdr = rdr.read_header().unwrap();
    let mut buf = RecordBuf::default();
    let mut out = Vec::new();
    while rdr.read_record_buf(&hdr, &mut buf).unwrap() > 0 {
        let mm = match buf.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => s.to_vec(),
            other => panic!("MM {other:?}"),
        };
        let ml = match buf.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) => v.clone(),
            other => panic!("ML {other:?}"),
        };
        let mn = buf
            .data()
            .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
            .and_then(Value::as_int);
        out.push((
            String::from_utf8_lossy(buf.name().unwrap().as_ref()).into_owned(),
            mm,
            ml,
            mn,
        ));
    }
    out
}

/// Each flank keeps its own modification call, renumbered to the segment's
/// C occurrences, on BAM output and in tagged FASTQ output from uBAM and from
/// tagged FASTQ input.
#[test]
fn split_segments_shift_modification_calls_on_ubam_and_tagged_fastq() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.bam");
    write_mod_fixture(&input);

    let out = dir.path().join("out.bam");
    whittle()
        .args([
            "--quality-trim",
            "segments",
            "--quality-cutoff",
            "10",
            "--quiet",
        ])
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .assert()
        .success();
    assert_eq!(
        mod_records(&out),
        [
            (
                "r1_segment_1".into(),
                b"C+m,1;".to_vec(),
                vec![50],
                Some(60)
            ),
            (
                "r1_segment_2".into(),
                b"C+m,2;".to_vec(),
                vec![250],
                Some(60)
            ),
        ]
    );

    let tagged = dir.path().join("tagged.fastq");
    whittle()
        .args(["--quiet", "-i", input.to_str().unwrap()])
        .args(["-o", tagged.to_str().unwrap()])
        .assert()
        .success();
    let expected = [
        "@r1_segment_1\tMM:Z:C+m,1;\tML:B:C,50\tMN:i:60",
        "@r1_segment_2\tMM:Z:C+m,2;\tML:B:C,250\tMN:i:60",
    ];
    for source in [&input, &tagged] {
        let out = dir.path().join("out.fastq");
        whittle()
            .args([
                "--quality-trim",
                "segments",
                "--quality-cutoff",
                "10",
                "--fastq-tags",
                "MM,ML,MN",
                "--quiet",
            ])
            .args(["-i", source.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .assert()
            .success();
        let text = std::fs::read_to_string(&out).unwrap();
        let headers: Vec<&str> = text.lines().step_by(4).collect();
        assert_eq!(headers, expected, "{}", source.display());
        let seqs: Vec<&str> = text.lines().skip(1).step_by(4).collect();
        assert_eq!(seqs, ["CA".repeat(30), "CA".repeat(30)]);
    }
}
