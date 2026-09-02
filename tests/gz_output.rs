//! Gzip FASTQ coverage. Unspecified output remains plain FASTQ, requested gzip
//! output is finalized with a complete footer, and damaged gzip input fails
//! with one clearly attributed error.

use std::io::{Read, Write};

use assert_cmd::Command;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use predicates::prelude::*;

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(bytes).unwrap();
    enc.finish().unwrap()
}

#[test]
fn truncated_gz_input_fails_once_with_record_context() {
    // Enough records that several parser buffers of records precede the
    // truncation point, so the error names the last record read.
    let mut fastq = String::new();
    for i in 0..20000 {
        fastq.push_str(&format!(
            "@r{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"
        ));
    }
    let gz = gzip(fastq.as_bytes());
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("trunc.fastq.gz");
    std::fs::write(&input, &gz[..gz.len() * 6 / 10]).unwrap();

    let assert = whittle()
        .arg("-i")
        .arg(&input)
        .args(["-t", "1"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("reading FASTQ record after r"));
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
    let failed = stderr
        .lines()
        .find(|l| l.contains("Failed after"))
        .unwrap_or_else(|| panic!("no failure line in: {stderr}"));
    let cause = failed.rsplit(": ").next().unwrap();
    assert_eq!(
        failed.matches(cause).count(),
        1,
        "the cause must be printed once: {failed}"
    );
}

#[test]
fn quality_byte_outside_phred33_range_is_a_hard_error() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("bad.fastq.gz");
    std::fs::write(&input, gzip(b"@r1\nACGTACGT\n+\nII I\x01III\n")).unwrap();

    whittle()
        .arg("-i")
        .arg(&input)
        .args(["-t", "1"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("record r1"))
        .stderr(predicate::str::contains("0x20"));
}

#[test]
fn plain_output_by_default_even_with_gz_input() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fastq.gz");

    // Build a small gzipped FASTQ input.
    let mut enc = GzEncoder::new(
        std::fs::File::create(&input).unwrap(),
        Compression::default(),
    );
    enc.write_all(b"@r1\nACGTACGTAC\n+\nIIIIIIIIII\n").unwrap();
    enc.finish().unwrap();

    // Input compression does not implicitly compress stdout.
    let assert = whittle()
        .arg("-i")
        .arg(&input)
        .args(["-H", "2", "-T", "2", "-t", "4"])
        .assert()
        .success();

    let stdout = assert.get_output().stdout.clone();
    assert_ne!(
        &stdout[..2.min(stdout.len())],
        &[0x1f, 0x8b][..],
        "stdout must be plain FASTQ, not gzip, when no output format is requested"
    );
    assert!(
        stdout.starts_with(b"@"),
        "expected plain FASTQ starting with '@', got {stdout:?}"
    );
    // ACGTACGTAC (10 bases), head-crop 2 + tail-crop 2 -> [2,8) = "GTACGT".
    assert_eq!(stdout, b"@r1\nGTACGT\n+\nIIIIII\n");
}

#[test]
fn explicit_gz_output_roundtrips_through_parallel_encoder() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fastq");
    std::fs::write(&input, "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n").unwrap();
    let out = dir.path().join("out.fastq.gz");

    // -t 4: exercise gzp's multi-threaded encoder, not just the trivial
    // single-thread case.
    whittle()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["-H", "2", "-T", "2", "-t", "4"])
        .assert()
        .success();

    // A missing `finish()` would leave this truncated/corrupt; decoding must
    // succeed and match the expected trimmed record exactly.
    let mut gz = MultiGzDecoder::new(std::fs::File::open(&out).unwrap());
    let mut s = String::new();
    gz.read_to_string(&mut s).unwrap();
    // ACGTACGTAC (10 bases), head-crop 2 + tail-crop 2 -> [2,8) = "GTACGT".
    assert_eq!(s, "@r1\nGTACGT\n+\nIIIIII\n");
}
