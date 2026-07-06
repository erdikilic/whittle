use assert_cmd::Command;
use predicates::prelude::*;

fn whittle() -> Command {
    Command::cargo_bin("whittle").unwrap()
}

#[test]
fn head_tail_crop_over_stdin() {
    whittle()
        .args([
            "--head-crop",
            "1",
            "--tail-crop",
            "1",
            "--in-format",
            "fastq",
        ])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .success()
        .stdout("@r1\nCG\n+\nII\n");
}

#[test]
fn mutually_exclusive_quality_ops_error() {
    whittle()
        .args([
            "--trim-qual",
            "10",
            "--best-segment",
            "10",
            "--in-format",
            "fastq",
        ])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("mutually exclusive"));
}

#[test]
fn min_length_filters() {
    whittle()
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

    whittle()
        .arg("-i")
        .arg(&path)
        .arg("-o")
        .arg(&path)
        .assert()
        .failure()
        .stderr(predicate::str::contains("same file"));

    assert_eq!(
        std::fs::read(&path).unwrap(),
        before,
        "input must not be modified"
    );
}

#[test]
fn contradictory_length_bounds_error() {
    whittle()
        .args(["-l", "10", "-L", "5", "--in-format", "fastq"])
        .write_stdin("@r1\nACGTACGTAC\n+\nIIIIIIIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("min-length"));
}

#[test]
fn contradictory_qual_bounds_error() {
    whittle()
        .args(["-q", "30", "-Q", "20", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("min-qual"));
}

#[test]
fn out_of_range_gc_bound_errors() {
    whittle()
        .args(["--min-gc", "2", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("min-gc"));
}

#[test]
fn nan_quality_bound_errors() {
    // NaN slips past `min > max` (every NaN comparison is false) and would silently
    // disable quality filtering; it must be rejected explicitly.
    whittle()
        .args(["--min-qual", "nan", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("NaN"));
}

#[test]
#[cfg(unix)]
fn hard_linked_input_output_is_rejected_and_preserves_input() {
    // Two hard links to one inode canonicalize to distinct paths, so only the
    // inode+device check catches this — otherwise File::create truncates the input.
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fastq");
    let output = dir.path().join("out.fastq");
    std::fs::write(&input, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n").unwrap();
    std::fs::hard_link(&input, &output).unwrap();
    let before = std::fs::read(&input).unwrap();

    whittle()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .assert()
        .failure()
        .stderr(predicate::str::contains("same file"));

    assert_eq!(
        std::fs::read(&input).unwrap(),
        before,
        "hard-linked input must be preserved"
    );
}

#[test]
fn gz_output_roundtrips() {
    use std::io::Read;
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.fastq.gz");
    whittle()
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
