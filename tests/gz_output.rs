//! Gzip FASTQ coverage. Unspecified output remains plain FASTQ, requested gzip
//! output is finalized with a complete footer, and damaged gzip input fails
//! with one clearly attributed error.

use std::io::{Read, Write};

use assert_cmd::Command;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use predicates::prelude::*;

/// The binary with `WHITTLE_LOG` cleared.
fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

/// `bytes` as a single-member gzip stream.
fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(bytes).unwrap();
    enc.finish().unwrap()
}

/// A truncated gzip input fails with one error line that names the last record
/// read, whether it inflates inline (`-t 1`) or on its own thread.
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

    for threads in ["1", "4"] {
        let assert = whittle()
            .arg("-i")
            .arg(&input)
            .args(["-t", threads])
            .assert()
            .failure()
            .stderr(predicate::str::contains("reading FASTQ record after r"));
        let stderr = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
        let failed = stderr
            .lines()
            .find(|l| l.contains("Failed after"))
            .unwrap_or_else(|| panic!("No failure line in: {stderr}"));
        let cause = failed.rsplit(": ").next().unwrap();
        assert_eq!(
            failed.matches(cause).count(),
            1,
            "The cause is printed once: {failed}"
        );
    }
}

/// Multi-member gzip on stdin inflates on its own thread under `-t 4` and
/// yields every member's records in order.
#[test]
fn multi_member_gz_on_stdin_reads_every_member() {
    let mut plain = String::new();
    let mut gz = Vec::new();
    for member in 0..3 {
        let mut text = String::new();
        for i in 0..5000 {
            text.push_str(&format!("@m{member}r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n"));
        }
        gz.extend(gzip(text.as_bytes()));
        plain.push_str(&text);
    }
    let assert = whittle()
        .args(["-t", "4", "--preserve-order", "--quiet"])
        .write_stdin(gz)
        .assert()
        .success();
    assert_eq!(String::from_utf8_lossy(&assert.get_output().stdout), plain);
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
        "Stdout is plain FASTQ, not gzip, when no output format is requested"
    );
    assert!(
        stdout.starts_with(b"@"),
        "Expected plain FASTQ starting with '@', got {stdout:?}"
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

    // -t 4 takes the parallel path, where the render workers compress the
    // BGZF blocks; the output is multi-member gzip that a plain gzip decoder
    // reads back.
    whittle()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&out)
        .args(["-H", "2", "-T", "2", "-t", "4"])
        .assert()
        .success();

    // A missing `finish()` would leave the file truncated; decoding succeeds and
    // matches the trimmed record.
    let mut gz = MultiGzDecoder::new(std::fs::File::open(&out).unwrap());
    let mut s = String::new();
    gz.read_to_string(&mut s).unwrap();
    // ACGTACGTAC (10 bases), head-crop 2 + tail-crop 2 -> [2,8) = "GTACGT".
    assert_eq!(s, "@r1\nGTACGT\n+\nIIIIII\n");
}
