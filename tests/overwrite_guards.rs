//! Guards against whittle destroying a file it is also reading.
//!
//! Each case here was a real defect: the run exited 0 having replaced an input,
//! an adapter FASTA, or its own output with something else.

use std::process::{Command as StdCommand, Stdio};

use assert_cmd::Command;
use predicates::prelude::*;

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

const READS: &str = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";

#[test]
fn summary_json_may_not_overwrite_the_input() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fastq");
    std::fs::write(&input, READS).unwrap();

    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args(["--summary-json", input.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--summary-json"));

    assert_eq!(
        std::fs::read_to_string(&input).unwrap(),
        READS,
        "the input must be untouched"
    );
}

#[test]
fn summary_json_may_not_overwrite_the_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fastq");
    let out = dir.path().join("out.fastq");
    std::fs::write(&input, READS).unwrap();

    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", out.to_str().unwrap()])
        .args(["--summary-json", out.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--summary-json"));
}

#[test]
fn output_may_not_overwrite_the_adapter_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fastq");
    let fasta = dir.path().join("adapters.fasta");
    std::fs::write(&input, READS).unwrap();
    let fasta_text = ">ada1\nACGTACGTACGTACGT\n";
    std::fs::write(&fasta, fasta_text).unwrap();

    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-a", fasta.to_str().unwrap()])
        .args(["-o", fasta.to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--adapter-fasta"));

    assert_eq!(
        std::fs::read_to_string(&fasta).unwrap(),
        fasta_text,
        "the adapter FASTA must be untouched"
    );
}

#[test]
fn summary_json_may_not_overwrite_a_folder_member() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in");
    std::fs::create_dir(&input).unwrap();
    let member = input.join("a.fastq");
    std::fs::write(&member, READS).unwrap();
    std::fs::write(input.join("b.fastq"), READS).unwrap();

    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", dir.path().join("merged.fastq").to_str().unwrap()])
        .args(["--summary-json", member.to_str().unwrap()])
        .assert()
        .failure();

    assert_eq!(
        std::fs::read_to_string(&member).unwrap(),
        READS,
        "the folder member must be untouched"
    );
}

/// `whittle -o reads.fastq < reads.fastq` has no `-i` for the path comparison to
/// use, but fd 0 still resolves to the file's inode. Without this check the
/// output truncated the file mid-read and destroyed everything past the first
/// buffer.
#[test]
fn stdin_redirect_may_not_overwrite_its_own_source() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("reads.fastq");
    let body = READS.repeat(20_000);
    std::fs::write(&path, &body).unwrap();
    let before = std::fs::metadata(&path).unwrap().len();

    let out = StdCommand::new(assert_cmd::cargo::cargo_bin("whittle"))
        .args(["-o", path.to_str().unwrap()])
        .stdin(Stdio::from(std::fs::File::open(&path).unwrap()))
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .output()
        .unwrap();

    assert!(!out.status.success(), "the run must be refused");
    assert_eq!(
        std::fs::metadata(&path).unwrap().len(),
        before,
        "the source file must not be truncated"
    );
    assert_eq!(std::fs::read_to_string(&path).unwrap(), body);
}

/// A pipe on stdin shares no inode with any path, so the guard must not fire.
#[test]
fn piped_stdin_is_not_mistaken_for_the_output() {
    let dir = tempfile::tempdir().unwrap();
    whittle()
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args(["-l", "5"])
        .write_stdin(READS)
        .assert()
        .success();
}
