use std::io::Write;

use assert_cmd::Command;

#[test]
fn old_qual_flag_is_gone_new_one_listed() {
    Command::cargo_bin("whittle")
        .unwrap()
        .args(["--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("--qual-split"))
        .stdout(predicates::str::contains("--qual-trim"))
        .stdout(predicates::str::contains("--qual-best-segment"));
}

#[test]
fn adapter_help_lists_flags() {
    Command::cargo_bin("whittle")
        .unwrap()
        .args(["--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("--adapter-fasta"))
        .stdout(predicates::str::contains("--adapter-preset"))
        .stdout(predicates::str::contains("--adapter-error-rate"))
        .stdout(predicates::str::contains("--adapter-ends-only"));
}

#[test]
fn adapter_sample_flag_listed() {
    Command::cargo_bin("whittle")
        .unwrap()
        .args(["--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("--adapter-sample"));
}

#[test]
fn rejects_out_of_range_error_rate() {
    // error-rate is only validated when an adapter source is active.
    Command::cargo_bin("whittle")
        .unwrap()
        .args([
            "--adapter-preset",
            "ont",
            "--adapter-error-rate",
            "2.0",
            "-i",
            "/dev/null",
        ])
        .assert()
        .failure();
}

#[test]
fn fastq_end_adapter_is_trimmed() {
    let adapter = "ACGTACGTACGTACGTACGT"; // 20 bp
    let insert = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 40 bp
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">a\n{adapter}").unwrap();
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    writeln!(fq, "@r1\n{adapter}{insert}\n+\n{}", "I".repeat(60)).unwrap();

    let out = Command::cargo_bin("whittle")
        .unwrap()
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let out = String::from_utf8(out).unwrap();
    assert!(out.contains(insert), "insert kept");
    assert!(
        !out.contains(&format!("{adapter}{insert}")),
        "adapter trimmed off"
    );
}

#[test]
fn adapter_fasta_with_no_usable_entries_errors() {
    // only a too-short entry -> skipped -> zero usable adapters.
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">short\nACGT").unwrap();
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    writeln!(fq, "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII").unwrap();
    Command::cargo_bin("whittle")
        .unwrap()
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("no usable adapters"));
}

// Build a FASTQ where every read starts with adapter A (present) and none
// contains catalog barcode-ish B; with a custom 2-adapter FASTA, detection
// should keep only A. We assert the log line AND that trimming still works.
#[test]
fn detection_keeps_present_drops_absent_and_still_trims() {
    let present = "GGGGTTTTGGGGTTTTGGGG"; // 20bp present adapter
    let absent = "ACGACGACGACGACGACGAC"; // 20bp never in the reads
    let insert = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 40bp
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">present\n{present}\n>absent\n{absent}").unwrap();
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..200 {
        writeln!(fq, "@r{i}\n{present}{insert}\n+\n{}", "I".repeat(60)).unwrap();
    }
    let out = Command::cargo_bin("whittle")
        .unwrap()
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "-v",
        ]) // -v so the info detection line reaches stderr
        .assert()
        .success();
    let res = out.get_output();
    let stderr = String::from_utf8_lossy(&res.stderr);
    assert!(stderr.contains("kept 1 of 2 adapters"), "stderr: {stderr}");
    let stdout = String::from_utf8_lossy(&res.stdout);
    assert!(stdout.contains(insert), "insert kept");
    assert!(
        !stdout.contains(&format!("{present}{insert}")),
        "present adapter trimmed"
    );
}

#[test]
fn adapter_sample_zero_disables_detection() {
    let present = "GGGGTTTTGGGGTTTTGGGG";
    let insert = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">present\n{present}").unwrap();
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..200 {
        writeln!(fq, "@r{i}\n{present}{insert}\n+\n{}", "I".repeat(60)).unwrap();
    }
    let res = Command::cargo_bin("whittle")
        .unwrap()
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "--adapter-sample",
            "0",
            "-v",
        ])
        .assert()
        .success();
    let stderr = String::from_utf8_lossy(&res.get_output().stderr);
    assert!(
        !stderr.contains("Adapter presence"),
        "detection must be off: {stderr}"
    );
}

#[test]
fn tiny_input_skips_detection() {
    let present = "GGGGTTTTGGGGTTTTGGGG";
    let insert = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">present\n{present}").unwrap();
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..10 {
        writeln!(fq, "@r{i}\n{present}{insert}\n+\n{}", "I".repeat(60)).unwrap();
    }
    let res = Command::cargo_bin("whittle")
        .unwrap()
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "-v",
        ])
        .assert()
        .success();
    let stderr = String::from_utf8_lossy(&res.get_output().stderr);
    assert!(
        stderr.contains("using all"),
        "tiny input must skip detection: {stderr}"
    );
}

// The key correctness guarantee (Phase 1.5 spec): detection only ever *drops*
// adapters that don't act, so trimming a present adapter with detection ON
// must be byte-identical to trimming with detection OFF (`--adapter-sample
// 0`, i.e. the full set). Here one adapter is present and several absent
// adapters are C/G-rich 20mers that cannot match the G/T-only present
// adapter or the A-run insert -- neither forward nor reverse-complement --
// within the default edit-distance budget (k_end = floor(0.2 * 20) = 4):
// each absent sequence has >= 5 C's AND >= 5 G's, so both it and its
// revcomp need >= 5 edits to align anywhere in a read that contains no C at
// all, which exceeds k_end. That was confirmed empirically too: the
// detection run's stderr reports all 3 absent adapters dropped (see below).
#[test]
fn detection_output_equals_full_set_for_present_adapter() {
    let present = "GGGGTTTTGGGGTTTTGGGG"; // 20bp, G/T only
    let absent = [
        "CCCCGGGGCCCCGGGGCCCC", // 12 C, 8 G
        "ACGACGACGACGACGACGAC", // 7 C, 6 G
        "CCGGCCGGCCGGCCGGCCGG", // 10 C, 10 G
    ];
    let insert = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 40bp
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">present\n{present}").unwrap();
    for (i, seq) in absent.iter().enumerate() {
        writeln!(fa, ">absent{i}\n{seq}").unwrap();
    }
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..200 {
        writeln!(fq, "@r{i}\n{present}{insert}\n+\n{}", "I".repeat(60)).unwrap();
    }

    let run = |extra_args: &[&str]| {
        // `-t 1`: pin single-threaded so output order is deterministic (the
        // parallel writer path lands records in arrival order, not input
        // order, for `threads > 1` -- see `pipeline::fastq::run`). That
        // nondeterminism is orthogonal to what this test checks, so it must
        // be controlled for to get a meaningful byte comparison.
        let mut args = vec![
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "-t",
            "1",
        ];
        args.extend_from_slice(extra_args);
        Command::cargo_bin("whittle")
            .unwrap()
            .args(args)
            .assert()
            .success()
            .get_output()
            .stdout
            .clone()
    };

    // Detection ON (default sampling): drops the 3 absent adapters.
    let detect_on = run(&[]);
    // Detection OFF: trims against the full 4-adapter set.
    let detect_off = run(&["--adapter-sample", "0"]);

    assert!(
        !detect_on.is_empty(),
        "detection-on output must be non-empty"
    );
    let detect_on_str = String::from_utf8(detect_on.clone()).unwrap();
    assert!(
        detect_on_str.contains(insert),
        "detection-on output should keep the insert: {detect_on_str}"
    );
    assert_eq!(
        detect_on, detect_off,
        "trimming a present adapter must be byte-identical whether detection \
         is on (drops absent adapters) or off (uses the full set), since \
         detection only removes adapters that don't act"
    );
}

#[test]
fn no_adapter_flag_is_byte_identical() {
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    write!(fq, "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n").unwrap();
    let out = Command::cargo_bin("whittle")
        .unwrap()
        .args(["-i", fq.path().to_str().unwrap()])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(
        String::from_utf8(out).unwrap(),
        "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
    );
}
