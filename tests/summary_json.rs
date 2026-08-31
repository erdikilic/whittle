//! End-to-end `--summary-json` over the compiled binary. The unit tests in
//! `src/summary.rs` cover the schema shape; these cover the parts only a real
//! run can show: that the file is written on every dispatch path, that it
//! survives `--quiet`, and that its counters agree with the reads that actually
//! came out.

use assert_cmd::Command;
use predicates::prelude::*;

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

/// Four reads of 20 bp; two are all-Q40, two carry a low-quality tail.
fn reads() -> String {
    let mut s = String::new();
    for i in 1..=2 {
        s.push_str(&format!(
            "@hi{i}\n{}\n+\n{}\n",
            "ACGT".repeat(5),
            "I".repeat(20)
        ));
    }
    let low_qual = format!("{}{}", "I".repeat(10), "#".repeat(10));
    for i in 1..=2 {
        s.push_str(&format!("@lo{i}\n{}\n+\n{low_qual}\n", "ACGT".repeat(5)));
    }
    s
}

fn summary(dir: &std::path::Path) -> serde_json::Value {
    let text = std::fs::read_to_string(dir.join("summary.json")).unwrap();
    serde_json::from_str(&text).unwrap()
}

#[test]
fn summary_counts_match_the_run() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.fastq");
    let json = dir.path().join("summary.json");

    whittle()
        .args(["-o", out.to_str().unwrap()])
        .args(["--summary-json", json.to_str().unwrap()])
        .args(["-l", "15"])
        .args(["--qual-trim", "20"])
        .write_stdin(reads())
        .assert()
        .success();

    let v = summary(dir.path());
    assert_eq!(v["schema_version"], 1);
    assert_eq!(v["tool"], "whittle");
    assert_eq!(v["reads"]["input"], 4);

    // The two clean reads survive intact; the two low-quality ones trim to 10 bp
    // and are then dropped by -l 15.
    assert_eq!(v["reads"]["output"], 2);
    assert_eq!(v["reads"]["all_filtered"], 2);
    assert_eq!(v["reads"]["with_output"], 2);
    assert_eq!(v["segments_dropped"]["too_short"], 2);
    assert_eq!(v["bases"]["input"], 80);
    assert_eq!(v["bases"]["output"], 40);

    // The counters agree with what actually landed in the output file.
    let written = std::fs::read_to_string(&out).unwrap();
    assert_eq!(written.lines().count(), 8);
}

/// A pipeline that asks for the summary gets it even when it silences logging;
/// that is the whole point of the flag.
#[test]
fn summary_is_written_under_quiet() {
    let dir = tempfile::tempdir().unwrap();
    let json = dir.path().join("summary.json");

    whittle()
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args(["--summary-json", json.to_str().unwrap()])
        .args(["-l", "5"])
        .arg("--quiet")
        .write_stdin(reads())
        .assert()
        .success()
        // --quiet drops the banner, progress, and the human-readable summary;
        // only warnings and errors survive, and this run raises neither.
        .stderr(predicate::str::is_empty());

    assert_eq!(summary(dir.path())["reads"]["input"], 4);
}

#[test]
fn summary_records_the_resolved_parameters() {
    let dir = tempfile::tempdir().unwrap();
    let json = dir.path().join("summary.json");

    whittle()
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args(["--summary-json", json.to_str().unwrap()])
        .args(["-t", "2", "-H", "3", "-T", "4", "-l", "5", "-q", "7"])
        .args(["--qual-split", "9", "--qual-split-window", "6"])
        .write_stdin(reads())
        .assert()
        .success();

    let v = summary(dir.path());
    assert_eq!(v["params"]["threads"], 2);
    assert_eq!(v["params"]["head_crop"], 3);
    assert_eq!(v["params"]["tail_crop"], 4);
    assert_eq!(v["params"]["min_length"], 5);
    assert_eq!(v["params"]["min_qual"], 7.0);
    assert_eq!(v["params"]["quality_op"]["mode"], "split");
    assert_eq!(v["params"]["quality_op"]["threshold"], 9);
    assert_eq!(v["params"]["quality_op"]["window"], 6);
    assert_eq!(v["input"], "<stdin>");
    assert!(
        v["command"]
            .as_str()
            .unwrap()
            .contains("--qual-split-window")
    );
    assert!(v["elapsed_seconds"].as_f64().unwrap() >= 0.0);
}

/// Folder-merge mode goes through its own dispatch arm, so it needs its own
/// proof that the summary lands.
#[test]
fn summary_is_written_for_a_folder_run() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in");
    std::fs::create_dir(&input).unwrap();
    std::fs::write(input.join("a.fastq"), "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    std::fs::write(input.join("b.fastq"), "@r2\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    let json = dir.path().join("summary.json");

    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args(["--summary-json", json.to_str().unwrap()])
        .assert()
        .success();

    let v = summary(dir.path());
    assert_eq!(v["reads"]["input"], 2);
    assert_eq!(v["reads"]["output"], 2);
    assert_eq!(v["bases"]["output"], 16);
}

/// An unwritable summary path is a hard error: silently skipping it would let a
/// pipeline consume a stale file from a previous run.
#[test]
fn unwritable_summary_path_fails_the_run() {
    let dir = tempfile::tempdir().unwrap();

    whittle()
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args([
            "--summary-json",
            dir.path()
                .join("no-such-dir/summary.json")
                .to_str()
                .unwrap(),
        ])
        .write_stdin(reads())
        .assert()
        .failure()
        .stderr(predicate::str::contains("writing summary JSON"));
}

/// Without the flag, nothing is written and no stray file appears.
#[test]
fn no_summary_file_without_the_flag() {
    let dir = tempfile::tempdir().unwrap();

    whittle()
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .write_stdin(reads())
        .assert()
        .success();

    assert!(!dir.path().join("summary.json").exists());
}

/// `--quiet` silences the human-readable summary but must not blank out the
/// JSON's `elapsed_seconds`, which is the field a benchmarking pipeline reads.
#[test]
fn elapsed_seconds_survives_quiet() {
    let dir = tempfile::tempdir().unwrap();
    let json = dir.path().join("summary.json");

    whittle()
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args(["--summary-json", json.to_str().unwrap()])
        .args(["-l", "5"])
        .arg("--quiet")
        .write_stdin(reads())
        .assert()
        .success();

    let v = summary(dir.path());
    assert!(
        v["elapsed_seconds"].as_f64().is_some(),
        "elapsed_seconds must be a number under --quiet, got {}",
        v["elapsed_seconds"]
    );
}

/// Report mode writes no records, so it cannot write a summary either. Exiting 0
/// with neither file and no explanation would strand a pipeline waiting on one.
#[test]
fn report_mode_warns_that_summary_json_is_ignored() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fastq");
    let mut fq = String::new();
    for i in 0..200 {
        let body: String = (0..300)
            .map(|j| b"ACGT"[((i * 7 + j * 13 + j * j) % 4) as usize] as char)
            .collect();
        let s = format!("CCTGTACTTCGTTCAGTTACGTATTGC{body}");
        fq.push_str(&format!("@r{i}\n{s}\n+\n{}\n", "I".repeat(s.len())));
    }
    std::fs::write(&input, fq).unwrap();
    let json = dir.path().join("summary.json");

    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .args(["--summary-json", json.to_str().unwrap()])
        .args(["--adapter-infer", "report"])
        .assert()
        .success()
        .stderr(predicate::str::contains("--summary-json is ignored"))
        .stderr(predicate::str::contains("-o/--output is ignored"));

    assert!(!json.exists(), "report mode writes no summary file");
}

/// Argument-parsing diagnostics are collected and emitted through tracing, so the
/// level filter applies. Previously they were printed straight to stderr from
/// `cli::parse`, which bypassed the filter entirely: an informational advisory
/// showed up even under `--quiet`.
#[test]
fn parse_time_advisories_respect_the_level_filter() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fastq");
    std::fs::write(&input, "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n").unwrap();
    let out = dir.path().join("out.fastq");

    // Conservative inference raises an INFO advisory, which --quiet must drop.
    let args = [
        "-i",
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--adapter-preset",
        "ont",
        "--adapter-infer",
        "trim",
    ];

    whittle()
        .args(args)
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "conservative adapter inference trims read ends only",
        ));

    whittle()
        .args(args)
        .arg("--quiet")
        .assert()
        .success()
        // The INFO advisory is filtered out; the WARN one still belongs on stderr,
        // since --quiet keeps warnings.
        .stderr(
            predicate::str::contains("conservative adapter inference")
                .not()
                .and(predicate::str::contains("--adapter-preset is ignored")),
        );
}

/// Deferred advisories still carry the standard `[timestamp] [LEVEL]` prefix, and
/// still land after the version and command lines that open every run.
#[test]
fn parse_time_advisories_are_formatted_and_ordered() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fastq");
    std::fs::write(&input, "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n").unwrap();

    let out = whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .arg("--adapter-ends-only")
        .assert()
        .success()
        .get_output()
        .stderr
        .clone();
    let stderr = String::from_utf8_lossy(&out);

    let advisory = stderr
        .lines()
        .find(|l| l.contains("--adapter-ends-only has no effect"))
        .expect("advisory present");
    assert!(
        advisory.contains("[WARN]") && advisory.starts_with('['),
        "advisory must carry the standard prefix: {advisory}"
    );

    let version_at = stderr.find("whittle 0.").expect("version line present");
    let advisory_at = stderr
        .find("--adapter-ends-only has no effect")
        .expect("advisory present");
    assert!(
        version_at < advisory_at,
        "the version line must still open the run"
    );
}
