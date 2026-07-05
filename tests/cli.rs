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
fn same_input_output_file_is_rejected_and_preserves_input() {
    // Streaming the input while truncating it on File::create destroys the data
    // (a plain FASTQ run silently emitted an empty file with a success exit).
    // The run must fail up front and leave the input untouched.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("reads.fastq");
    std::fs::write(&path, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n").unwrap();
    let before = std::fs::read(&path).unwrap();

    chopping()
        .arg("-i").arg(&path)
        .arg("-o").arg(&path)
        .assert()
        .failure()
        .stderr(predicate::str::contains("same file"));

    assert_eq!(std::fs::read(&path).unwrap(), before, "input must not be modified");
}

#[test]
fn contradictory_length_bounds_error() {
    chopping()
        .args(["-l", "10", "-L", "5", "--in-format", "fastq"])
        .write_stdin("@r1\nACGTACGTAC\n+\nIIIIIIIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("min-length"));
}

#[test]
fn contradictory_qual_bounds_error() {
    chopping()
        .args(["-q", "30", "-Q", "20", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("min-qual"));
}

#[test]
fn out_of_range_gc_bound_errors() {
    chopping()
        .args(["--min-gc", "2", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("min-gc"));
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
