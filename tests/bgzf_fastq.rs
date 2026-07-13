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

    let mut reader = noodles_bgzf::io::Reader::new(std::fs::File::open(output).unwrap());
    let mut decoded = Vec::new();
    reader.read_to_end(&mut decoded).unwrap();
    assert_eq!(decoded, TRIMMED);
}
