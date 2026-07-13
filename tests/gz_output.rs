//! Compressed FASTQ output coverage. Unspecified output remains plain FASTQ,
//! while requested gzip output is finalized with a complete footer.

use std::io::{Read, Write};

use assert_cmd::Command;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
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
