//! CLI behavior over the compiled binary: argument validation, logging modes,
//! guardrail warnings, and the failure path.

use std::io::Read as _;
use std::path::Path;
use std::process::Stdio;

use assert_cmd::Command;
use predicates::prelude::*;

/// The binary with `WHITTLE_LOG` cleared, so the environment does not change
/// the log mode.
fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

/// One 20-bp Q40 read.
const READS: &str = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";

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
        .args(["-vvv", "--in-format", "fastq"])
        .write_stdin("")
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

/// Streaming the input while `File::create` truncates it destroys the data. The
/// run fails up front with one `Failed after` line and leaves the input
/// untouched.
#[test]
fn same_input_output_file_is_rejected_and_preserves_input() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("reads.fastq");
    std::fs::write(&path, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n").unwrap();
    let before = std::fs::read(&path).unwrap();

    let assert = whittle()
        .arg("-i")
        .arg(&path)
        .arg("-o")
        .arg(&path)
        .assert()
        .failure()
        .stderr(predicate::str::contains("same file"));
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr);
    assert_eq!(
        stderr.matches("Failed after").count(),
        1,
        "Exactly one failure line is printed: {stderr}"
    );

    assert_eq!(
        std::fs::read(&path).unwrap(),
        before,
        "Input must not be modified"
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

/// NaN slips past `min > max` (every NaN comparison is false) and would disable
/// quality filtering, so it is rejected explicitly.
#[test]
fn nan_quality_bound_errors() {
    whittle()
        .args(["--min-qual", "nan", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("NaN"));
}

/// How a run ends for one argument set.
enum Expect {
    /// Exit non-zero with the substring on stderr.
    Fails,
    /// Exit zero with the substring on stderr.
    Warns,
}

/// Every validation and advisory names the flag it concerns, so the offending
/// argument is identifiable from the message alone.
#[test]
fn every_validation_names_its_flag() {
    use std::io::Write as _;

    let dir = tempfile::tempdir().unwrap();
    let fastq = dir.path().join("reads.fastq");
    std::fs::write(&fastq, READS).unwrap();
    // gzip bytes behind a plain `.fastq` name, so `--in-format fastq-gz` is right
    // and the extension is the one that looks wrong.
    let misnamed_gz = dir.path().join("misnamed.fastq");
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    enc.write_all(READS.as_bytes()).unwrap();
    std::fs::write(&misnamed_gz, enc.finish().unwrap()).unwrap();
    let folder = dir.path().join("folder");
    std::fs::create_dir(&folder).unwrap();
    std::fs::write(folder.join("a.fastq"), READS).unwrap();
    let missing = dir.path().join("no-such-adapters.fasta");
    let p = |path: &Path| path.to_str().unwrap().to_string();

    let cases: Vec<(Vec<String>, Expect, &str)> = vec![
        (
            vec!["--compression-level".into(), "10".into()],
            Expect::Fails,
            "--compression-level must be between 0 and 9",
        ),
        (
            vec!["--max-gc".into(), "1.5".into()],
            Expect::Fails,
            "--max-gc (1.5) must be a fraction",
        ),
        (
            vec![
                "--min-gc".into(),
                "0.8".into(),
                "--max-gc".into(),
                "0.2".into(),
            ],
            Expect::Fails,
            "--min-gc (0.8) must not exceed --max-gc",
        ),
        (
            vec!["--max-qual".into(), "nan".into()],
            Expect::Fails,
            "NaN",
        ),
        (
            vec!["--min-qual=-5".into()],
            Expect::Fails,
            "--min-qual (-5) must be a finite quality",
        ),
        (
            vec!["--max-qual".into(), "inf".into()],
            Expect::Fails,
            "--max-qual (inf) must be a finite quality",
        ),
        (
            vec!["-t".into(), "0".into()],
            Expect::Fails,
            "-t/--threads must be at least 1",
        ),
        (
            vec!["--fastq-tags".into(), "MM,ABC".into()],
            Expect::Fails,
            "--fastq-tags: invalid tag",
        ),
        (
            vec!["--qual-split-window".into(), "5".into()],
            Expect::Fails,
            "--qual-split",
        ),
        (
            vec!["--progress".into(), "bar".into(), "--quiet".into()],
            Expect::Fails,
            "cannot be used with",
        ),
        (
            vec!["--adapter-fasta".into(), p(&missing)],
            Expect::Fails,
            "--adapter-fasta",
        ),
        (
            vec![
                "--adapter-preset".into(),
                "ont".into(),
                "--adapter-end-size".into(),
                "0".into(),
            ],
            Expect::Fails,
            "--adapter-end-size must be >= 1",
        ),
        (
            vec![
                "--adapter-infer".into(),
                "--adapter-sample".into(),
                "0".into(),
            ],
            Expect::Fails,
            "--adapter-sample 0 disables sampling",
        ),
        (
            vec!["--adapter-error-rate".into(), "0.5".into()],
            Expect::Fails,
            "--adapter-error-rate requires an adapter source",
        ),
        (
            vec!["--adapter-end-size".into(), "50".into()],
            Expect::Fails,
            "--adapter-end-size requires an adapter source",
        ),
        (
            vec!["--adapter-sample".into(), "50".into()],
            Expect::Fails,
            "--adapter-sample requires an adapter source",
        ),
        (
            vec![
                "-i".into(),
                p(&misnamed_gz),
                "--in-format".into(),
                "fastq-gz".into(),
                "-o".into(),
                p(&dir.path().join("out_in_mismatch.fastq")),
            ],
            Expect::Warns,
            "--in-format",
        ),
        (
            vec![
                "-i".into(),
                p(&fastq),
                "-o".into(),
                p(&dir.path().join("out_out_mismatch.fastq.gz")),
                "--out-format".into(),
                "fastq".into(),
            ],
            Expect::Warns,
            "--out-format",
        ),
        (
            vec![
                "-i".into(),
                p(&folder),
                "--in-format".into(),
                "fastq".into(),
                "-o".into(),
                p(&dir.path().join("out_folder.fastq")),
            ],
            Expect::Warns,
            "--in-format is ignored for a directory input",
        ),
        (
            vec![
                "-i".into(),
                p(&fastq),
                "-o".into(),
                p(&dir.path().join("out_noop.fastq")),
            ],
            Expect::Warns,
            "No trimming or filtering options set",
        ),
        (
            vec![
                "-i".into(),
                p(&fastq),
                "-o".into(),
                p(&dir.path().join("out.bam")),
            ],
            Expect::Fails,
            "FASTQ-to-BAM conversion is not supported",
        ),
        (
            vec![
                "-i".into(),
                p(&fastq),
                "-o".into(),
                p(&dir.path().join("out_sum.fastq")),
                "--summary-json".into(),
                p(dir.path()),
            ],
            Expect::Fails,
            "--summary-json",
        ),
    ];

    for (args, expect, needle) in cases {
        let mut cmd = whittle();
        cmd.args(&args);
        if !args.iter().any(|a| a == "-i") {
            cmd.args(["--in-format", "fastq"]).write_stdin(READS);
        }
        let output = cmd.output().unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        match expect {
            Expect::Fails => assert!(
                !output.status.success(),
                "Case {args:?} must fail: {stderr}"
            ),
            Expect::Warns => assert!(
                output.status.success(),
                "Case {args:?} must succeed: {stderr}"
            ),
        }
        assert!(
            stderr.contains(needle),
            "Case {args:?}: stderr lacks {needle:?}: {stderr}"
        );
    }
}

/// A missing input names its path, so a command line with several paths on it
/// says which one is wrong.
#[test]
fn missing_input_error_names_the_path() {
    let dir = tempfile::tempdir().unwrap();
    let missing = dir.path().join("nonexistent.fastq");
    whittle()
        .arg("-i")
        .arg(&missing)
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("opening input"))
        .stderr(predicate::str::contains("nonexistent.fastq"));
}

/// An unparseable `WHITTLE_LOG` does not disable logging: the level falls back
/// to the default, the failure line still prints, and the rejected directive is
/// named.
#[test]
fn invalid_whittle_log_falls_back_and_still_reports_the_failure() {
    let dir = tempfile::tempdir().unwrap();
    let missing = dir.path().join("nonexistent.fastq");
    whittle()
        .env("WHITTLE_LOG", "garbage=nope=1")
        .arg("-i")
        .arg(&missing)
        .args(["-o", dir.path().join("out.fastq").to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains("Failed after"))
        .stderr(predicate::str::contains("WHITTLE_LOG"))
        .stderr(predicate::str::contains("garbage=nope=1"));
}

/// `-` is the pipeline spelling of stdin and stdout; no file named `-` appears.
#[test]
fn dash_means_stdin_and_stdout() {
    let dir = tempfile::tempdir().unwrap();
    whittle()
        .current_dir(dir.path())
        .args(["-i", "-", "-o", "-", "--in-format", "fastq", "-l", "5"])
        .write_stdin(READS)
        .assert()
        .success()
        .stdout(READS)
        .stderr(predicate::str::contains("Input: <stdin>"))
        .stderr(predicate::str::contains("Output: <stdout>"));
    assert!(
        !dir.path().join("-").exists(),
        "No file named - may be created"
    );
}

/// A downstream reader that stops early (`whittle ... | head`) is not a failure
/// of this run: the process exits 0 without a `Failed after` line.
#[test]
fn downstream_closing_the_pipe_exits_quietly() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("reads.fastq");
    // Well past the pipe buffer, so writes are still in flight when the reader
    // goes away.
    std::fs::write(&input, READS.repeat(40_000)).unwrap();

    for threads in [&["-t", "1"][..], &[][..]] {
        let mut child = std::process::Command::new(assert_cmd::cargo::cargo_bin("whittle"))
            .env_remove("WHITTLE_LOG")
            .args(["-i", input.to_str().unwrap()])
            .args(threads)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let mut stdout = child.stdout.take().unwrap();
        let mut head = [0u8; 100];
        stdout.read_exact(&mut head).unwrap();
        drop(stdout);

        let output = child.wait_with_output().unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            output.status.success(),
            "A closed pipe exits 0 ({threads:?}): {stderr}"
        );
        assert!(
            !stderr.contains("Failed after"),
            "A closed pipe is not reported as a failure ({threads:?}): {stderr}"
        );
    }
}

/// Two hard links to one inode canonicalize to distinct paths, so only the inode
/// and device check catches this; otherwise `File::create` truncates the input.
#[test]
#[cfg(unix)]
fn hard_linked_input_output_is_rejected_and_preserves_input() {
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
        "Hard-linked input must be preserved"
    );
}

/// A default run prints the `Summary:` line to stderr; `--quiet` drops it and
/// keeps the reads on stdout.
#[test]
fn quiet_suppresses_summary_but_keeps_stdout() {
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
        .stderr(predicate::str::contains("Processing started")); // the phase timing line appears at DEBUG
}

#[test]
fn default_hides_debug() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("Processing started").not());
}

/// `--quiet` wins even when `WHITTLE_LOG` asks for verbose output.
#[test]
fn quiet_beats_whittle_log() {
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

/// Without `--quiet`, `WHITTLE_LOG` raises the level above the CLI default.
#[test]
fn whittle_log_overrides_verbosity_when_not_quiet() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .env("WHITTLE_LOG", "debug")
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("Processing started"));
}

/// assert_cmd captures stderr to a pipe (non-TTY), so the run is in line mode
/// regardless of verbosity. The full startup banner and the `Completed` closer
/// appear in order.
#[test]
fn line_mode_banner_and_closer_appear_in_order() {
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

/// The version and command lines precede the resolved configuration and the
/// diagnostics.
#[test]
fn banner_version_and_command_come_first_in_line_mode() {
    let input = "@r1\nACGT\n+\nIIII\n";
    let assert = whittle().write_stdin(input).assert().success();
    let stderr = String::from_utf8(assert.get_output().stderr.clone()).unwrap();
    let version_pos = stderr.find("whittle ").expect("Version line present");
    let command_pos = stderr.find("Command:").expect("Command line present");
    let operation_pos = stderr.find("Trimming").expect("Operation line present");
    assert!(
        version_pos < command_pos && command_pos < operation_pos,
        "Expected version, then Command:, then the operation line, in order: {stderr:?}"
    );
}

/// Captured stderr is non-interactive and contains no ANSI escapes.
#[test]
fn non_tty_stderr_has_no_ansi_escapes() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .write_stdin(input)
        .assert()
        .success()
        .stderr(predicate::str::contains("\u{1b}").not());
}

/// Every read fails an unreachable min-qual bound: nothing survives, and the run
/// succeeds. The all-dropped guardrail WARN fires so the run is distinguishable
/// from a clean empty-output run.
#[test]
fn all_dropped_run_warns() {
    let input = "@r1\nACGT\n+\nIIII\n";
    whittle()
        .args(["-q", "50", "--in-format", "fastq"])
        .write_stdin(input)
        .assert()
        .success()
        .stderr(
            predicate::str::contains("No reads survived")
                .and(predicate::str::contains("every input read was dropped")),
        );
}

/// Zero input reads is not an error, and the empty-input guardrail WARN still
/// fires.
#[test]
fn empty_input_warns() {
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

    // A minimal one-read uBAM built in the tempdir keeps the test hermetic.
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
