use assert_cmd::Command;
use predicates::prelude::*;

fn chopping() -> Command {
    Command::cargo_bin("chopping").unwrap()
}

#[test]
fn head_tail_crop_over_stdin() {
    chopping()
        .args(["--head-crop", "1", "--tail-crop", "1", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .success()
        .stdout("@r1\nCG\n+\nII\n");
}

#[test]
fn mutually_exclusive_quality_ops_error() {
    chopping()
        .args(["--trim-qual", "10", "--best-segment", "10", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("mutually exclusive"));
}

#[test]
fn min_length_filters() {
    chopping()
        .args(["--min-length", "10", "--in-format", "fastq"])
        .write_stdin("@short\nACGT\n+\nIIII\n")
        .assert()
        .success()
        .stdout(""); // filtered out
}

#[test]
fn gz_output_roundtrips() {
    use std::io::Read;
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.fastq.gz");
    chopping()
        .args(["--in-format", "fastq", "-o"])
        .arg(&out)
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .success();
    let mut gz = flate2::read::MultiGzDecoder::new(std::fs::File::open(&out).unwrap());
    let mut s = String::new();
    gz.read_to_string(&mut s).unwrap();
    assert_eq!(s, "@r1\nACGT\n+\nIIII\n");
}
