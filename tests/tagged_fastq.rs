//! FASTQ input whose headers carry SAM aux tags (`samtools fastq -T`): the
//! tags are rewritten per output segment exactly as on the BAM-to-FASTQ path.

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

/// Writes a three-read uBAM: a plain read, a modified read with a per-base
/// array and a scalar, and a read with a `bi` barcode span.
fn write_fixture(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut r1 = RecordBuf::default();
    *r1.flags_mut() = Flags::UNMAPPED;
    *r1.name_mut() = Some(b"read1".into());
    *r1.sequence_mut() = b"ACGTACGTACGTACGTACGT".to_vec().into();
    *r1.quality_scores_mut() = vec![40; 20].into();
    w.write_alignment_record(&header, &r1).unwrap();

    // C at 0,1,3,4,5,7,8,9,11,12,13,15; MM occurrences 0,2,3,6 -> positions
    // 0,3,4,8.
    let mut r2 = RecordBuf::default();
    *r2.flags_mut() = Flags::UNMAPPED;
    *r2.name_mut() = Some(b"read2".into());
    *r2.sequence_mut() = b"CCACCCACCCACCCAC".to_vec().into();
    *r2.quality_scores_mut() = vec![35; 16].into();
    let d = r2.data_mut();
    d.insert(Tag::from(*b"RG"), Value::String(b"grp1".to_vec().into()));
    d.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1,0,2;".to_vec().into()),
    );
    d.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![10, 20, 30, 40])),
    );
    d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(16));
    d.insert(
        Tag::from(*b"pw"),
        Value::Array(Array::UInt8((1..=16).collect())),
    );
    d.insert(Tag::from(*b"qs"), Value::Float(35.0));
    w.write_alignment_record(&header, &r2).unwrap();

    // Front barcode over bases 0..4 and rear barcode over 16..20.
    let mut r3 = RecordBuf::default();
    *r3.flags_mut() = Flags::UNMAPPED;
    *r3.name_mut() = Some(b"read3".into());
    *r3.sequence_mut() = b"TTTTACGTACGTACGTAAAA".to_vec().into();
    *r3.quality_scores_mut() = vec![30; 20].into();
    r3.data_mut().insert(
        Tag::from(*b"bi"),
        Value::Array(Array::Float(vec![90.0, 0.0, 4.0, 90.0, 20.0, 4.0, 90.0])),
    );
    w.write_alignment_record(&header, &r3).unwrap();

    w.try_finish().unwrap();
}

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

/// Converts the fixture to tagged FASTQ through whittle's own BAM-to-FASTQ
/// path with no trimming, and returns both paths.
fn fixture(dir: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let bam = dir.join("reads.bam");
    let tagged = dir.join("tagged.fastq");
    write_fixture(&bam);
    whittle()
        .args([
            "-i",
            bam.to_str().unwrap(),
            "-o",
            tagged.to_str().unwrap(),
            "--quiet",
        ])
        .assert()
        .success();
    let text = std::fs::read_to_string(&tagged).unwrap();
    assert!(
        text.contains("@read2\tRG:Z:grp1\t") && text.contains("\tMM:Z:C+m,0,1,0,2;"),
        "{text}"
    );
    (bam, tagged)
}

/// Runs `args` over `input` and returns the output text.
fn run(args: &[&str], input: &Path, dir: &Path, name: &str) -> String {
    let out = dir.join(name);
    whittle()
        .args(args)
        .args([
            "-i",
            input.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--quiet",
        ])
        .assert()
        .success();
    std::fs::read_to_string(out).unwrap()
}

#[test]
fn tagged_fastq_matches_the_bam_path_on_every_trim() {
    let dir = tempfile::tempdir().unwrap();
    let (bam, tagged) = fixture(dir.path());
    for (i, args) in [
        vec!["-H", "3", "-T", "2"],
        vec![
            "--split-quality",
            "36",
            "--split-min-low-quality-bases",
            "1",
            "-l",
            "1",
        ],
        vec!["--best-quality-segment", "36"],
        vec!["--trim-barcodes", "-H", "1"],
        vec!["--remove-kinetics", "--remove-tag", "RG", "-T", "5"],
        vec!["--fastq-tags", "MM,ML,MN", "-H", "2"],
        vec!["--fastq-tags", "none", "-H", "2"],
    ]
    .iter()
    .enumerate()
    {
        let from_bam = run(args, &bam, dir.path(), &format!("bam{i}.fastq"));
        let from_tagged = run(args, &tagged, dir.path(), &format!("tagged{i}.fastq"));
        assert_eq!(from_tagged, from_bam, "args {args:?}");
    }
}

#[test]
fn tagged_fastq_rewrites_mods_and_slices_arrays() {
    let dir = tempfile::tempdir().unwrap();
    let (_bam, tagged) = fixture(dir.path());
    let out = run(&["-H", "3", "-T", "2"], &tagged, dir.path(), "out.fastq");
    let head = out.lines().find(|l| l.starts_with("@read2")).unwrap();
    // Window [3,14): C at 3,4,5,7,8,9,11,12,13 -> occurrences 0..9; the calls
    // at 3,4,8 survive as occurrences 0,1,4 with ML [20,30,40]; MN is 11.
    assert!(head.contains("\tMM:Z:C+m,0,0,2;"), "{head}");
    assert!(head.contains("\tML:B:C,20,30,40"), "{head}");
    assert!(head.contains("\tMN:i:11"), "{head}");
    assert!(
        head.contains("\tpw:B:C,4,5,6,7,8,9,10,11,12,13,14"),
        "{head}"
    );
    assert!(head.contains("\tRG:Z:grp1"), "{head}");
    assert!(head.contains("\tqs:f:35"), "{head}");
}

#[test]
fn tagged_fastq_over_stdin_and_gzip() {
    let dir = tempfile::tempdir().unwrap();
    let (_bam, tagged) = fixture(dir.path());
    let expected = run(&["-H", "3"], &tagged, dir.path(), "file.fastq");
    let text = std::fs::read(&tagged).unwrap();
    whittle()
        .args(["-H", "3", "--input-format", "fastq", "--quiet"])
        .write_stdin(text.clone())
        .assert()
        .success()
        .stdout(expected.clone());
    let gz = dir.path().join("tagged.fastq.gz");
    let mut enc = flate2::write::GzEncoder::new(
        std::fs::File::create(&gz).unwrap(),
        flate2::Compression::default(),
    );
    std::io::Write::write_all(&mut enc, &text).unwrap();
    enc.finish().unwrap();
    assert_eq!(run(&["-H", "3"], &gz, dir.path(), "gz.fastq"), expected);
}

#[test]
fn tag_flags_need_tags_in_the_input() {
    let dir = tempfile::tempdir().unwrap();
    let plain = dir.path().join("plain.fastq");
    std::fs::write(&plain, "@r1 runid=x\nACGT\n+\nIIII\n").unwrap();
    for (args, msg) in [
        (
            vec!["--trim-barcodes"],
            "--trim-barcodes reads barcode spans",
        ),
        (vec!["--remove-tag", "RG"], "--remove-tag removes aux tags"),
        (
            vec!["--remove-kinetics"],
            "--remove-kinetics removes aux tags",
        ),
    ] {
        whittle()
            .args(&args)
            .args(["-i", plain.to_str().unwrap(), "--quiet"])
            .assert()
            .failure()
            .stderr(predicates::str::contains(msg))
            .stderr(predicates::str::contains("FASTQ without header tags"));
    }
    // A tab that is not a SAM field leaves the header alone.
    let comment = dir.path().join("comment.fastq");
    std::fs::write(&comment, "@r1\tsome comment\nACGT\n+\nIIII\n").unwrap();
    whittle()
        .args(["-i", comment.to_str().unwrap(), "-H", "1", "--quiet"])
        .assert()
        .success()
        .stdout("@r1\tsome comment\nCGT\n+\nIII\n");
}

#[test]
fn malformed_tag_names_the_read() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("bad.fastq");
    std::fs::write(
        &input,
        "@r1\tMM:Z:C+m,0;\tML:B:C,5\nCCAC\n+\nIIII\n@r2\tML:B:C,x\nACGT\n+\nIIII\n",
    )
    .unwrap();
    whittle()
        .args(["-i", input.to_str().unwrap(), "-H", "1", "--quiet"])
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "read r2: header field \"ML:B:C,x\": array value \"x\" does not fit the subtype",
        ));
}

#[test]
fn late_tags_are_rewritten_in_files_and_directories() {
    let dir = tempfile::tempdir().unwrap();
    let folder = dir.path().join("inputs");
    std::fs::create_dir(&folder).unwrap();
    let plain: String = (0..160)
        .map(|i| format!("@plain{i}\nCCCC\n+\nIIII\n"))
        .collect();
    let tagged = "@tagged\tMM:Z:C+m,1;\tML:B:C,200\tMN:i:4\tRG:Z:group\nCCCC\n+\nIIII\n";
    std::fs::write(folder.join("1.fastq"), &plain).unwrap();
    std::fs::write(folder.join("2.fastq"), tagged).unwrap();
    let file = dir.path().join("all.fastq");
    std::fs::write(&file, plain + tagged).unwrap();
    for input in [&file, &folder] {
        for threads in ["1", "4"] {
            let output = run(
                &["-H", "1", "-t", threads, "--remove-tag", "RG"],
                input,
                dir.path(),
                "out.fastq",
            );
            let head = output
                .lines()
                .find(|line| line.starts_with("@tagged"))
                .unwrap();
            assert!(head.contains("\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:3"), "{head}");
            assert!(!head.contains("\tRG:"), "{head}");
            assert_eq!(output.lines().count(), 161 * 4);
        }
    }
}

#[test]
fn split_identifiers_precede_header_descriptions() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("descriptions.fastq");
    std::fs::write(&input, "@r1 description\tMM:Z:C+m,0;\tML:B:C,200\nCCCCC\n+\nII!II\n@plain another description\nCCCCC\n+\nII!II\n").unwrap();
    for threads in ["1", "4"] {
        let output = run(
            &["--split-quality", "10", "-t", threads, "--preserve-order"],
            &input,
            dir.path(),
            "out.fastq",
        );
        let heads: Vec<_> = output.lines().step_by(4).collect();
        assert!(heads[0].starts_with("@r1_segment_1 description\t"));
        assert!(heads[1].starts_with("@r1_segment_2 description\t"));
        assert!(heads[0].contains("\tpi:Z:r1\t"));
        assert_eq!(heads[2], "@plain_segment_1 another description");
        assert_eq!(heads[3], "@plain_segment_2 another description");
    }
}

#[test]
fn impossible_modification_positions_are_removed_and_counted() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("invalid.fastq");
    std::fs::write(&input, "@r1\tMM:Z:C+m,8;\tML:B:C,200\tMN:i:4\nCCCC\n+\nIIII\n@r2\tMM:Z:N+n,184467440737095516160;\tML:B:C,3\nACGT\n+\nIIII\n").unwrap();
    let summary = dir.path().join("summary.json");
    for threads in ["1", "4"] {
        for crop in ["0", "1"] {
            let output = run(
                &[
                    "-H",
                    crop,
                    "-t",
                    threads,
                    "--summary-json",
                    summary.to_str().unwrap(),
                ],
                &input,
                dir.path(),
                "out.fastq",
            );
            assert!(!output.contains("MM:"));
            assert!(!output.contains("ML:"));
            assert!(!output.contains("MN:"));
            let stats: serde_json::Value =
                serde_json::from_slice(&std::fs::read(&summary).unwrap()).unwrap();
            assert_eq!(stats["warnings"]["malformed_mod_reads"], 2);
        }
    }
}

#[test]
fn pacbio_interval_names_follow_repeated_crops() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("pacbio.fastq");
    std::fs::write(&input, "@movie/1/100_106\tqs:i:100\tqe:i:106\nCCCCCC\n+\nIIIIII\n@movie/2/ccs/100_106\tRG:Z:group\nCCCCCC\n+\nIIIIII\n").unwrap();
    let first = run(
        &["-H", "1", "--preserve-order"],
        &input,
        dir.path(),
        "first.fastq",
    );
    assert!(first.starts_with("@movie/1/101_106\tqs:i:101\tqe:i:106\n"));
    assert!(first.contains("@movie/2/ccs/101_106\t"));
    let second = run(
        &["-T", "1", "--preserve-order"],
        &dir.path().join("first.fastq"),
        dir.path(),
        "second.fastq",
    );
    assert!(second.starts_with("@movie/1/101_105\tqs:i:101\tqe:i:105\n"));
    assert!(second.contains("@movie/2/ccs/101_105\t"));
}
