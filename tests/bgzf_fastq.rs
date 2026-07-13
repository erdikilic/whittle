use std::io::{Read, Write};

use assert_cmd::Command;

const INPUT: &[u8] = b"@r1\nACGTACGT\n+\nIIIIIIII\n@r2 note\nTTGGCCAA\n+\nHHHHHHHH\n";
const TRIMMED: &[u8] = b"@r1\nCGTACG\n+\nIIIIII\n@r2 note\nTGGCCA\n+\nHHHHHH\n";

fn bgzf_bytes() -> Vec<u8> {
    let mut writer = noodles_bgzf::io::Writer::new(Vec::new());
    writer.write_all(INPUT).unwrap();
    writer.finish().unwrap()
}

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

#[test]
fn bgzf_fastq_path_and_stdin_are_detected() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fastq.bgz");
    std::fs::write(&input, bgzf_bytes()).unwrap();

    whittle()
        .arg("-i")
        .arg(&input)
        .args(["-H", "1", "-T", "1", "-t", "4", "--quiet"])
        .assert()
        .success()
        .stdout(TRIMMED);

    whittle()
        .args(["-H", "1", "-T", "1", "-t", "4", "--quiet"])
        .write_stdin(bgzf_bytes())
        .assert()
        .success()
        .stdout(TRIMMED);
}

// The 28-byte BGZF EOF marker (an empty gzip block). samtools/bgzip require it
// to treat a .bgz stream as complete; a writer that omits it makes `bgzip -t`
// report truncation even when every data block is intact.
const BGZF_EOF: &[u8] = &[
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

#[test]
fn explicit_bgzf_output_roundtrips_and_finishes() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fastq");
    let output = dir.path().join("trimmed.fastq.bgz");
    std::fs::write(&input, INPUT).unwrap();

    whittle()
        .arg("-i")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["-H", "1", "-T", "1", "-t", "4", "--quiet"])
        .assert()
        .success();

    let raw = std::fs::read(&output).unwrap();
    assert!(
        raw.ends_with(BGZF_EOF),
        "bgzf output must end with the EOF marker, else samtools/bgzip sees a truncated file"
    );

    let mut reader = noodles_bgzf::io::Reader::new(raw.as_slice());
    let mut decoded = Vec::new();
    reader.read_to_end(&mut decoded).unwrap();
    assert_eq!(decoded, TRIMMED);
}
