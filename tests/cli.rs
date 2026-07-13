use assert_cmd::Command;
use predicates::prelude::*;

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

#[test]
fn version_is_long_only() {
    whittle()
        .arg("--version")
        .assert()
        .success()
        .stdout(predicate::str::contains(env!("CARGO_PKG_VERSION")));

    whittle()
        .arg("-V")
        .assert()
        .failure()
        .stderr(predicate::str::contains("unexpected argument '-V'"));
}

#[test]
fn verbosity_above_trace_is_rejected() {
    whittle()
        .args(["-vvv", "-i", "/dev/null"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("accepts at most -vv"));
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
            "--qual-trim",
            "10",
            "--qual-best-segment",
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
fn quiet_suppresses_summary_but_keeps_stdout() {
    // A minimal FASTQ over stdin; default run prints the "Summary:" line to stderr.
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .arg("--quiet")
        .write_stdin(input)
        .assert()
        .success()
        .stdout(predicate::str::contains("@r1"))
        .stderr(predicate::str::contains("Summary:").not());
}

#[test]
fn default_run_prints_summary_to_stderr() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle().write_stdin(input).assert().success().stderr(
        predicate::str::contains("Summary:")
            .and(predicate::str::contains("input reads"))
            .and(predicate::str::contains("output reads")),
    );
}

#[test]
fn over_spec_threads_warns() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .args(["-t", "100000"])
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("exceeds").and(predicate::str::contains("using")));
}

#[test]
fn verbose_shows_phase_timing() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .arg("-v")
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("Processing")); // phase timing line appears at DEBUG
}

#[test]
fn default_hides_debug() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("Processing").not());
}

#[test]
fn quiet_beats_whittle_log() {
    // --quiet must win even when WHITTLE_LOG asks for verbose output.
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .env("WHITTLE_LOG", "debug")
        .arg("--quiet")
        .write_stdin(input)
        .assert()
        .success()
        .stdout(predicate::str::contains("@r1"))
        .stderr(predicate::str::contains("Summary:").not());
}

#[test]
fn whittle_log_overrides_verbosity_when_not_quiet() {
    // Without --quiet, WHITTLE_LOG still raises the level above the CLI's default.
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .env("WHITTLE_LOG", "debug")
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("Processing"));
}

#[test]
fn line_mode_banner_and_closer_appear_in_order() {
    // assert_cmd captures stderr to a pipe (non-tty), so this always runs in
    // line mode regardless of verbosity — the full startup banner plus the
    // Completed closer should appear, in order.
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle().write_stdin(input).assert().success().stderr(
        predicate::str::contains("whittle ")
            .and(predicate::str::contains("Command:"))
            .and(predicate::str::contains("Trimming"))
            .and(predicate::str::contains("Input: <stdin>"))
            .and(predicate::str::contains("Output: <stdout>"))
            .and(predicate::str::contains("Threads:"))
            .and(predicate::str::contains("Filters:"))
            .and(predicate::str::contains("Summary:"))
            .and(predicate::str::contains("Completed in")),
    );
}

#[test]
fn failure_path_prints_a_single_failed_after_line() {
    // The reworked `main.rs` failure path must print one clean "Failed after
    // ...: <message>" line via tracing (not a second, differently-formatted
    // anyhow dump from the default `fn main() -> anyhow::Result<()>` pattern).
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("reads.fastq");
    std::fs::write(&path, "@r1\nACGT\n+\nIIII\n").unwrap();

    whittle()
        .arg("-i")
        .arg(&path)
        .arg("-o")
        .arg(&path)
        .assert()
        .failure()
        .stderr(
            predicate::str::contains("Failed after").and(predicate::str::contains("same file")),
        );
}

#[test]
fn banner_version_and_command_come_first_in_line_mode() {
    // Version and command precede resolved configuration and diagnostics.
    let input = "@r1\nACGT\n+\nIIII\n";
    let assert = whittle().write_stdin(input).assert().success();
    let stderr = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    let version_pos = stderr.find("whittle ").expect("version line missing");
    let command_pos = stderr.find("Command:").expect("command line missing");
    let operation_pos = stderr.find("Trimming").expect("operation line missing");
    assert!(
        version_pos < command_pos && command_pos < operation_pos,
        "expected version, then Command:, then the operation line, in order: {stderr:?}"
    );
}

#[test]
fn non_tty_stderr_has_no_ansi_escapes() {
    // Captured stderr is non-interactive and must contain no ANSI escapes.
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("\u{1b}").not());
}

#[test]
fn all_dropped_run_warns() {
    // Every read fails an unreachable min-qual bound: nothing survives, but
    // the run itself still succeeds — the all-dropped guardrail WARN must
    // fire so this doesn't silently look like a clean empty-output run.
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .args(["-q", "50", "--in-format", "fastq"])
        .write_stdin(input)
        .assert()
        .success()
        .stderr(
            predicate::str::contains("No reads survived")
                .and(predicate::str::contains("input reads were dropped")),
        );
}

#[test]
fn empty_input_warns() {
    // Zero input reads (not an error) must still surface the empty-input
    // guardrail WARN rather than a silent, unremarkable "0 in, 0 out" summary.
    whittle()
        .args(["--in-format", "fastq"])
        .write_stdin("")
        .assert()
        .success()
        .stderr(predicate::str::contains("Input contained no reads"));
}

#[test]
fn sequential_threads_label_for_dash_t_1() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .args(["-t", "1"])
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("Threads: 1 (sequential)"));
}

#[test]
fn bam_to_fastq_conversion_phrasing() {
    use noodles_bam as bam;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::{self as sam};

    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.bam");
    let out = dir.path().join("o.fastq");

    // Build a minimal one-read uBAM in the tempdir so the test is hermetic
    // (it must not depend on a fixture outside the repository).
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(&inp).unwrap());
    w.write_header(&header).unwrap();
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = b"ACGTACGTAC".to_vec().into();
    *rec.quality_scores_mut() = vec![40; 10].into();
    w.write_alignment_record(&header, &rec).unwrap();
    w.try_finish().unwrap();

    whittle()
        .arg("-i")
        .arg(&inp)
        .arg("-o")
        .arg(&out)
        .assert()
        .success()
        .stderr(predicate::str::contains("Converting BAM to FASTQ"));
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
