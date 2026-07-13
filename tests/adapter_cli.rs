use std::io::Write;

use assert_cmd::Command;

#[test]
fn old_qual_flag_is_gone_new_one_listed() {
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
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
        .env_remove("WHITTLE_LOG")
        .args(["--help"])
        .assert()
        .success()
        .stdout(predicates::str::contains("--adapter-fasta"))
        .stdout(predicates::str::contains("--adapter-preset"))
        .stdout(predicates::str::contains("--adapter-error-rate"))
        .stdout(predicates::str::contains("--adapter-ends-only"))
        .stdout(predicates::str::contains("--adapter-infer-policy"));
}

#[test]
fn adapter_sample_flag_listed() {
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
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
        .env_remove("WHITTLE_LOG")
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
fn adapter_sample_below_min_is_rejected() {
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "--adapter-preset",
            "ont",
            "--adapter-sample",
            "50",
            "-i",
            "/dev/null",
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("must be 0"));
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
        .env_remove("WHITTLE_LOG")
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
        .env_remove("WHITTLE_LOG")
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

// Custom FASTA entries remain active without presence-based reduction, including
// entries absent from the sampled reads.
#[test]
fn custom_fasta_never_reduces_even_with_an_absent_adapter() {
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
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "-v",
        ])
        .assert()
        .success();
    let res = out.get_output();
    let stderr = String::from_utf8_lossy(&res.stderr);
    assert!(
        !stderr.contains("Adapter presence"),
        "custom --adapter-fasta must disable detection outright, so no reduction is ever logged: {stderr}"
    );
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
        .env_remove("WHITTLE_LOG")
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

// Preset detection skips samples below the minimum discovery size.
#[test]
fn tiny_input_skips_detection() {
    let front = "CCTGTACTTCGTTCAGTTACGTATTGC"; // LSK114 front, real preset entry
    let insert = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..10 {
        writeln!(
            fq,
            "@r{i}\n{front}{insert}\n+\n{}",
            "I".repeat(front.len() + insert.len())
        )
        .unwrap();
    }
    let res = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-preset",
            "ont",
            "--adapter-sample",
            "10000",
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

// Presence detection is opt-in: --adapter-sample defaults to 0 (off), so a
// preset run with no --adapter-sample flag at all must never engage
// detection -- no "Adapter presence" line at any point -- while still
// trimming the preset adapter that's present in every read, since with
// detection off the full (un-reduced) catalog is what gets searched.
#[test]
fn default_does_not_run_detection() {
    let front = "CCTGTACTTCGTTCAGTTACGTATTGC"; // LSK114 front, real preset entry
    // Long insert -- see the comment on
    // `detection_output_equals_full_set_for_present_adapter` for why a short
    // insert can be consumed entirely by the catalog's paired front/rear
    // entries within the default 150bp end-zone (unrelated to detection).
    let insert = "A".repeat(300);
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..200 {
        writeln!(
            fq,
            "@r{i}\n{front}{insert}\n+\n{}",
            "I".repeat(front.len() + insert.len())
        )
        .unwrap();
    }
    let res = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-preset",
            "ont",
            "-v",
        ])
        .assert()
        .success();
    let out = res.get_output();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("Adapter presence"),
        "detection must be off by default (no --adapter-sample given): {stderr}"
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains(&"A".repeat(280)),
        "insert must still be kept: {stdout}"
    );
    assert!(
        !stdout.contains(&format!("{front}{insert}")),
        "preset adapter must still be trimmed off by default (full-set search): {stdout}"
    );
}

// Presence detection may remove inactive preset entries but must not change the
// output for a present adapter. The long insert keeps the paired rear entry
// outside the default end-search region.
#[test]
fn detection_output_equals_full_set_for_present_adapter() {
    let front = "CCTGTACTTCGTTCAGTTACGTATTGC"; // LSK114 front, 27bp
    let insert = "A".repeat(300);
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..200 {
        writeln!(
            fq,
            "@r{i}\n{front}{insert}\n+\n{}",
            "I".repeat(front.len() + insert.len())
        )
        .unwrap();
    }

    let run = |extra_args: &[&str]| {
        // `-t 1`: pin single-threaded so output order is deterministic (the
        // parallel writer path lands records in arrival order, not input
        // order, for `threads > 1` -- see `workflow::fastq::run`). That
        // nondeterminism is orthogonal to what this test checks, so it must
        // be controlled for to get a meaningful byte comparison.
        let mut args = vec![
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-preset",
            "ont",
            "-t",
            "1",
            "-v",
        ];
        args.extend_from_slice(extra_args);
        let out = Command::cargo_bin("whittle")
            .unwrap()
            .env_remove("WHITTLE_LOG")
            .args(args)
            .assert()
            .success();
        let res = out.get_output();
        (res.stdout.clone(), res.stderr.clone())
    };

    // Detection ON (explicit opt-in sampling, since detection defaults to
    // off): reduces the 124-entry catalog.
    let (detect_on, detect_on_err) = run(&["--adapter-sample", "10000"]);
    // Detection OFF (the default: no --adapter-sample flag needed, but pass
    // 0 explicitly for clarity): trims against the full 124-adapter set.
    let (detect_off, _) = run(&["--adapter-sample", "0"]);

    assert!(
        !detect_on.is_empty(),
        "detection-on output must be non-empty"
    );
    let stderr_on = String::from_utf8_lossy(&detect_on_err);
    assert!(
        stderr_on.contains("Adapter presence: sampled") && !stderr_on.contains("kept 124 of 124"),
        "detection must genuinely engage (kept < full 124): {stderr_on}"
    );
    let detect_on_str = String::from_utf8(detect_on.clone()).unwrap();
    // A long run of the insert's A's, not the full 300 -- fuzzy end-matching
    // can legitimately consume a base or two of the boundary between adapter
    // and insert (an equal-cost alignment choice), which isn't what this
    // test is checking; the real guarantee is the byte-identical comparison
    // below.
    assert!(
        detect_on_str.contains(&"A".repeat(280)),
        "detection-on output should keep (most of) the insert: {detect_on_str}"
    );
    assert_eq!(
        detect_on, detect_off,
        "trimming a present adapter must be byte-identical whether detection \
         is on (drops absent catalog adapters) or off (uses the full set), \
         since detection only removes adapters that don't act"
    );
}

// An adapter-free prefix must not disable a custom FASTA for later reads, even
// when an adapter sample size is explicitly supplied.
#[test]
fn custom_fasta_trims_adapters_after_a_clean_prefix() {
    let adapter = "GGGGTTTTGGGGTTTTGGGG"; // 20bp, G/T only
    let insert = "A".repeat(40);
    let clean = "C".repeat(60); // pure C: can't match a G/T-only adapter within budget
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">present\n{adapter}").unwrap();
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..10_000 {
        writeln!(fq, "@c{i}\n{clean}\n+\n{}", "I".repeat(clean.len())).unwrap();
    }
    for i in 0..10 {
        writeln!(
            fq,
            "@a{i}\n{adapter}{insert}\n+\n{}",
            "I".repeat(adapter.len() + insert.len())
        )
        .unwrap();
    }

    let out = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "--adapter-sample",
            "10000",
            "-t",
            "1",
        ])
        .assert()
        .success();
    let res = out.get_output();
    let stdout = String::from_utf8_lossy(&res.stdout);
    let stderr = String::from_utf8_lossy(&res.stderr);

    assert!(
        !stderr.contains("Adapter presence: sampled"),
        "custom --adapter-fasta must disable detection outright: {stderr}"
    );
    assert!(
        !stdout.contains(&format!("{adapter}{insert}")),
        "the 10 adapted reads must be trimmed, not left untouched: {stdout}"
    );
    let trimmed_lines = stdout.lines().filter(|l| *l == insert).count();
    assert_eq!(
        trimmed_lines, 10,
        "all 10 adapted reads' inserts must survive trimming: {stdout}"
    );
}

// When a sampled prefix contains no preset adapters, detection falls back to
// the full preset so adapters in later reads are still processed.
#[test]
fn preset_detection_falls_back_when_prefix_has_no_adapters() {
    let front = "CCTGTACTTCGTTCAGTTACGTATTGC"; // LSK114 front, real preset entry
    let insert = "A".repeat(40);
    let clean = "C".repeat(60);
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..100 {
        writeln!(fq, "@c{i}\n{clean}\n+\n{}", "I".repeat(clean.len())).unwrap();
    }
    for i in 0..10 {
        writeln!(
            fq,
            "@a{i}\n{front}{insert}\n+\n{}",
            "I".repeat(front.len() + insert.len())
        )
        .unwrap();
    }

    let res = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-preset",
            "ont",
            "--adapter-sample",
            "100",
            "-v",
            "-t",
            "1",
        ])
        .assert()
        .success();
    let out = res.get_output();
    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert!(
        stderr.contains("no adapters detected"),
        "must fall back to the full set with a warning: {stderr}"
    );
    assert!(
        !stdout.contains(&format!("{front}{insert}")),
        "the 10 adapted reads must be trimmed via the full-set fallback: {stdout}"
    );
}

#[test]
fn no_adapter_flag_is_byte_identical() {
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    write!(fq, "@r1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n").unwrap();
    let out = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
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

#[test]
fn infer_and_fasta_are_mutually_exclusive() {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        "x.fastq",
        "--adapter-infer",
        "--adapter-fasta",
        "a.fa",
    ]);
    cmd.assert()
        .failure()
        .stderr(predicates::str::contains("mutually exclusive"));
}

#[test]
fn adapter_sample_below_min_still_rejected_under_infer() {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args(["-i", "x.fastq", "--adapter-infer", "--adapter-sample", "50"]);
    cmd.assert()
        .failure()
        .stderr(predicates::str::contains("must be 0"));
}

#[test]
fn infer_policy_requires_an_inference_operation() {
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["-i", "x.fastq", "--adapter-infer-policy", "aggressive"])
        .assert()
        .failure()
        .stderr(predicates::str::contains(
            "--adapter-infer [<ADAPTER_INFER>]",
        ));
}

#[test]
fn infer_policy_rejects_unknown_value() {
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            "x.fastq",
            "--adapter-infer",
            "--adapter-infer-policy",
            "balanced",
        ])
        .assert()
        .failure()
        .stderr(predicates::str::contains("invalid value 'balanced'"));
}

// Report-only inference names discoveries against the built-in catalog and
// user-supplied FASTA entries.
#[test]
fn infer_report_with_fasta_notes_naming_includes_fasta() {
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">present\nACGTACGTACGTACGTACGT").unwrap();
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    for i in 0..10 {
        writeln!(fq, "@r{i}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII").unwrap();
    }
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-infer",
            "report",
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stderr(predicates::str::contains("plus your FASTA's adapters"));
}

// Inference fixtures plant an exact adapter before deterministic, non-periodic
// genomic backgrounds. Error-tolerant recovery is covered by unit tests.

/// The 28bp adapter planted at the 5' end of every synthetic read below (an
/// SQK-NSK007/LSK109-neighborhood front sequence -- same one used by
/// `discover_recovers_planted_adapter_under_error` in src/adapter/infer.rs).
const PLANTED_ADAPTER: &str = "AATGTACTTCGTTCAGTTACGTATTGCT";

/// Length of the per-read genomic tail appended after `PLANTED_ADAPTER`.
const TAIL_LEN: usize = 120;

/// Deterministic, non-periodic genomic background for read `i`: a
/// splitmix64-style bit-mix seeded from the read index, matching
/// `src/adapter/infer.rs`'s `discover_*` unit-test fixtures exactly. Distinct
/// per `i` (each read gets its own splitmix64 state), and not periodic (so it
/// carries no spurious cross-read k-mer signal for the discoverer to flag).
fn splitmix_tail(i: usize, len: usize) -> Vec<u8> {
    let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        out.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
    }
    out
}

/// The full (untrimmed) synthetic sequence for read `i`: the exact planted
/// adapter followed by its splitmix64 tail. Used both to write the input
/// fixture and, independently, to recompute what a genuinely-trimmed output
/// record must be a suffix of.
fn full_read_seq(i: usize) -> Vec<u8> {
    let mut seq = PLANTED_ADAPTER.as_bytes().to_vec();
    seq.extend(splitmix_tail(i, TAIL_LEN));
    seq
}

/// Write `n` synthetic reads (see fixture notes above) to `<dir>/adapted.fastq`
/// and return its path. Read `i`'s id is `@r{i}` (no description), so a test
/// can parse the trailing digits back into the same index `full_read_seq`
/// used to build it.
fn write_adapted_fastq(dir: &std::path::Path, n: usize) -> std::path::PathBuf {
    let path = dir.join("adapted.fastq");
    let mut f = std::fs::File::create(&path).unwrap();
    for i in 0..n {
        let seq = full_read_seq(i);
        let qual = "I".repeat(seq.len());
        writeln!(
            f,
            "@r{i}\n{}\n+\n{qual}",
            std::str::from_utf8(&seq).unwrap()
        )
        .unwrap();
    }
    path
}

#[test]
fn infer_action_and_policy_map_to_banner() {
    let dir = tempfile::tempdir().unwrap();
    let fq = write_adapted_fastq(dir.path(), 10);

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["-i", fq.to_str().unwrap(), "--adapter-infer", "-t", "1"])
        .assert()
        .success()
        .stderr(predicates::str::contains("infer trim · conservative"));

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.to_str().unwrap(),
            "--adapter-infer",
            "report",
            "--adapter-infer-policy",
            "aggressive",
            "-t",
            "1",
        ])
        .assert()
        .success()
        .stderr(predicates::str::contains("infer report · aggressive"));
}

#[test]
fn infer_report_prints_and_does_not_trim() {
    let dir = tempfile::tempdir().unwrap();
    let fq = write_adapted_fastq(dir.path(), 500);
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "--adapter-infer",
        "report",
        "-t",
        "1",
    ]);
    let assert = cmd.assert().success();
    let out = assert.get_output();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("inferred") || stderr.contains("support"),
        "report-only must log what it discovered: {stderr}"
    );
    // Report-only exits before dispatch: no FASTQ record header on stdout.
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        !stdout.contains('@'),
        "report-only must not write any trimmed FASTQ to stdout: {stdout}"
    );
}

// Report-only emits discovered adapters as FASTA records on stdout.
#[test]
fn infer_report_prints_sequence_to_stdout() {
    let dir = tempfile::tempdir().unwrap();
    let fq = write_adapted_fastq(dir.path(), 500);
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "--adapter-infer",
        "report",
        "-t",
        "1",
    ]);
    let assert = cmd.assert().success();
    let stdout = String::from_utf8_lossy(&assert.get_output().stdout).into_owned();
    assert!(
        stdout.contains(PLANTED_ADAPTER),
        "stdout must contain the discovered adapter's actual bases, not just \
         its name/support: {stdout}"
    );
    assert!(
        stdout.trim_start().starts_with('>'),
        "stdout must be FASTA (header line starting with '>'): {stdout}"
    );
    assert!(
        stdout.contains("boundary=")
            && stdout.contains("assembled_length=")
            && stdout.contains("uncertain_bases="),
        "FASTA header must expose boundary uncertainty: {stdout}"
    );
}

// Cross-naming considers both the ONT catalog and the user's FASTA. The custom
// name sorts first and deterministically wins an equal-identity tie.
#[test]
fn infer_report_cross_names_against_user_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let fq = write_adapted_fastq(dir.path(), 500);
    // Filename deliberately has no "MY_CUSTOM_ADAPTER" substring, so a stray
    // path/filename echo elsewhere in the log could never produce a false
    // pass -- the assertion below can only be satisfied by the discovered
    // adapter's own cross-name.
    let fa_path = dir.path().join("cross_name_refs.fa");
    std::fs::write(
        &fa_path,
        format!(">AAA_MY_CUSTOM_ADAPTER\n{PLANTED_ADAPTER}\n"),
    )
    .unwrap();

    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "--adapter-infer",
        "report",
        "--adapter-fasta",
        fa_path.to_str().unwrap(),
        "-t",
        "1",
    ]);
    let assert = cmd.assert().success();
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
    assert!(
        stderr.contains("MY_CUSTOM_ADAPTER"),
        "discovered adapter must be cross-named against the user's --adapter-fasta, \
         not just the built-in catalog: {stderr}"
    );
}

#[test]
fn infer_trims_planted_adapter() {
    let dir = tempfile::tempdir().unwrap();
    let n = 500;
    let fq = write_adapted_fastq(dir.path(), n);
    let out_path = dir.path().join("out.fastq");
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "--adapter-infer",
        "-t",
        "1",
    ]);
    let assert = cmd.assert().success();
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
    assert!(
        stderr.contains("inferred"),
        "stderr must show an inferred-adapter log line: {stderr}"
    );

    let trimmed = std::fs::read_to_string(&out_path).unwrap();
    assert!(
        !trimmed.contains(PLANTED_ADAPTER),
        "the planted adapter must not survive anywhere in the output: {trimmed}"
    );

    // Each surviving sequence must be a suffix of its independently rebuilt
    // input, with a cut length consistent with the planted adapter.
    let mut lines = trimmed.lines();
    let mut n_records = 0;
    while let Some(header) = lines.next() {
        assert!(header.starts_with("@r"), "unexpected header: {header}");
        let idx: usize = header[2..].parse().expect("header must be @r<index>");
        let seq_line = lines.next().expect("sequence line");
        let _plus = lines.next().expect("plus line");
        let _qual = lines.next().expect("quality line");

        let original = full_read_seq(idx);
        // `checked_sub`: a clear panic message if output ever somehow
        // exceeded input, instead of an underflow wraparound with a cryptic
        // "attempt to subtract with overflow" pointing at this line only.
        let cut = original
            .len()
            .checked_sub(seq_line.len())
            .expect("output longer than input");
        assert!(
            original.ends_with(seq_line.as_bytes()),
            "record {idx}'s output must be an exact suffix of its original read"
        );
        assert!(
            (20..=50).contains(&cut),
            "record {idx}: cut length {cut} is not adapter-shaped (planted adapter is 28bp)"
        );
        n_records += 1;
    }
    assert_eq!(n_records, n, "no reads were dropped by trimming");
}

#[test]
fn infer_on_tiny_input_warns_and_keeps_reads() {
    let dir = tempfile::tempdir().unwrap();
    let n = 10; // < MIN_SAMPLE_FOR_DETECTION (100)
    let fq = write_adapted_fastq(dir.path(), n);
    let out_path = dir.path().join("out.fastq");
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "--adapter-infer",
        "-t",
        "1",
    ]);
    let assert = cmd.assert().success();
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
    assert!(
        stderr.contains("too few") || stderr.contains("no adapters"),
        "must warn about the undersized sample: {stderr}"
    );

    // Untrimmed: every output record equals its full original (unstripped)
    // read exactly -- the planted adapter is still there, verbatim.
    let trimmed = std::fs::read_to_string(&out_path).unwrap();
    for i in 0..n {
        let expected = String::from_utf8(full_read_seq(i)).unwrap();
        assert!(
            trimmed.contains(&expected),
            "record {i} must be kept untrimmed: {trimmed}"
        );
    }
}

// `--adapter-infer report` must not write or modify record output.

/// Report-only writes no records when discovery is skipped for insufficient
/// input reads.
#[test]
fn infer_report_tiny_input_writes_no_output() {
    let dir = tempfile::tempdir().unwrap();
    let n = 10; // < MIN_SAMPLE_FOR_DETECTION (100)
    let fq = write_adapted_fastq(dir.path(), n);
    let out_path = dir.path().join("out.fastq");
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "--adapter-infer",
        "report",
        "-t",
        "1",
    ]);
    let assert = cmd.assert().success();
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
    assert!(
        stderr.contains("too few"),
        "must warn about the undersized sample: {stderr}"
    );

    // Report-only must write no records: either the `-o` file was never
    // created, or it exists but is empty / has no FASTQ record header.
    match std::fs::read(&out_path) {
        Ok(bytes) => assert!(
            bytes.is_empty() || !bytes.contains(&b'@'),
            "report-only must not write any output records to -o: {:?}",
            String::from_utf8_lossy(&bytes)
        ),
        Err(e) => assert_eq!(
            e.kind(),
            std::io::ErrorKind::NotFound,
            "unexpected error reading -o file: {e}"
        ),
    }
}

/// Report-only leaves an existing output path unchanged.
#[test]
fn infer_report_does_not_clobber_output_file() {
    let dir = tempfile::tempdir().unwrap();
    // Adequate (>= MIN_SAMPLE_FOR_DETECTION) planted-adapter input, so
    // discovery actually runs (not the too-few-reads path exercised above).
    let fq = write_adapted_fastq(dir.path(), 500);
    let out_path = dir.path().join("existing.txt");
    let sentinel = "SENTINEL: pre-existing file contents, must survive\n";
    std::fs::write(&out_path, sentinel).unwrap();

    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "--adapter-infer",
        "report",
        "-t",
        "1",
    ]);
    cmd.assert().success();

    let contents = std::fs::read_to_string(&out_path).unwrap();
    assert_eq!(
        contents, sentinel,
        "report-only must not touch a pre-existing -o file at all"
    );
}

// --- determinism ---------------------------------------------------------

/// Same input, run twice through `--adapter-infer` at `-t 1`, must produce
/// byte-identical output. Discovery itself (`infer::discover`) is pure over
/// its sampled slice with no RNG or hashmap-iteration-order dependence, and
/// `-t 1` pins the FASTQ dispatch to its sequential (order-preserving) path,
/// so this is a black-box seal on that guarantee rather than new logic.
#[test]
fn infer_is_deterministic() {
    let dir = tempfile::tempdir().unwrap();
    let fq = write_adapted_fastq(dir.path(), 500);
    let run = |name: &str| {
        let out = dir.path().join(name);
        let mut cmd = Command::cargo_bin("whittle").unwrap();
        cmd.env_remove("WHITTLE_LOG");
        cmd.args([
            "-i",
            fq.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--adapter-infer",
            "-t",
            "1",
        ]);
        cmd.assert().success();
        std::fs::read(&out).unwrap()
    };
    assert_eq!(
        run("a.fastq"),
        run("b.fastq"),
        "same input -> byte-identical output"
    );
}

// --- marginal-support warning ------------------------------------------

/// Fraction of reads that carry `PLANTED_ADAPTER` at all; the rest are pure
/// background (no adapter anywhere in the read), modelling a low-prevalence
/// / barcode-specific adapter rather than a per-read match-quality problem.
/// Support is now a whole-consensus PRESENCE fraction (see
/// `infer::assemble`'s doc comment), so a per-read *error rate* no longer
/// drags a genuine adapter's support down (a real, closely-matching
/// reconstruction now recovers at support ~1.0, see
/// `discover_recovers_planted_adapter_under_error`) -- what still lands an
/// adapter in the marginal band is being present in only a *minority* of
/// reads. 0.38 * 500 = 190 planted reads out of 500 puts support at ~0.38,
/// inside `[KEEP_SUPPORT, MARGINAL_SUPPORT)` = `[0.30, 0.45)` with headroom
/// on both sides.
const PLANTED_ADAPTER_PREVALENCE: f64 = 0.38;

/// Fixture for the marginal-support warning: `PLANTED_ADAPTER_PREVALENCE` of
/// `n` reads get an EXACT copy of `PLANTED_ADAPTER` (no injected
/// substitution error -- error-tolerant recovery is already covered by
/// `discover_recovers_planted_adapter_under_error`, this fixture targets
/// marginal *prevalence* instead) followed by a splitmix64 tail; the
/// remaining reads are pure splitmix64 background of the same total length,
/// carrying no adapter at all (same non-periodic bit-mix pattern used
/// throughout this file and in `src/adapter/infer.rs`'s own `discover_*`
/// unit tests, so it can't itself register as a spurious low-complexity
/// signal).
fn write_adapted_fastq_marginal(dir: &std::path::Path, n: usize) -> std::path::PathBuf {
    // NOTE: deliberately not named "*marginal*" -- the path is itself echoed
    // into the `[INFO] Input: ...` / `Command: ...` log lines, which would
    // make `stderr.contains("marginal")` a false positive unrelated to the
    // actual warning message under test.
    let path = dir.join("weak_adapter.fastq");
    let mut f = std::fs::File::create(&path).unwrap();
    let planted_n = (n as f64 * PLANTED_ADAPTER_PREVALENCE).round() as usize;
    for i in 0..n {
        let seq: Vec<u8> = if i < planted_n {
            let mut s = PLANTED_ADAPTER.as_bytes().to_vec();
            s.extend(splitmix_tail(i, TAIL_LEN));
            s
        } else {
            // pure background, no adapter -- same total read length as the
            // planted branch so both groups look alike apart from content.
            splitmix_tail(i, TAIL_LEN + PLANTED_ADAPTER.len())
        };
        let qual = "I".repeat(seq.len());
        writeln!(
            f,
            "@r{i}\n{}\n+\n{qual}",
            std::str::from_utf8(&seq).unwrap()
        )
        .unwrap();
    }
    path
}

/// A kept adapter whose support sits in `[KEEP_SUPPORT, MARGINAL_SUPPORT)`
/// (here ~0.38, see `write_adapted_fastq_marginal`) must get an explicit
/// `warn!` in addition to the plain per-adapter info line, so a marginal
/// discovery doesn't read the same as a confident one.
#[test]
fn infer_warns_on_marginal_support() {
    let dir = tempfile::tempdir().unwrap();
    let fq = write_adapted_fastq_marginal(dir.path(), 500);
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "--adapter-infer",
        "report",
        "-t",
        "1",
    ]);
    let assert = cmd.assert().success();
    let stderr = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
    assert!(
        stderr.contains("marginal"),
        "a support just above KEEP_SUPPORT must be flagged marginal: {stderr}"
    );
}

// BAM report-only mode emits FASTA text and completes without record output.
fn write_minimal_ubam(path: &std::path::Path, n: usize) {
    use noodles_bam as bam;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::{self as sam};

    let header = sam::Header::default();
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();
    for i in 0..n {
        let seq = full_read_seq(i);
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::UNMAPPED;
        *r.name_mut() = Some(format!("r{i}").into_bytes().into());
        let qual = vec![40u8; seq.len()];
        *r.sequence_mut() = seq.into();
        *r.quality_scores_mut() = qual.into();
        writer.write_alignment_record(&header, &r).unwrap();
    }
    writer.try_finish().unwrap();
}

#[test]
fn infer_report_on_bam_input_with_piped_stdout_succeeds() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("in.bam");
    write_minimal_ubam(&in_path, 500);

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            in_path.to_str().unwrap(),
            "--adapter-infer",
            "report",
            "-t",
            "1",
        ])
        .assert()
        .success();
}

// End-to-end trim-then-filter behavior and accounting.

/// Quality filtering evaluates the surviving segment. The complete read has
/// mean quality 24.8; cropping four Q2 bases leaves a Q40 segment.
#[test]
fn quality_filter_judges_the_trimmed_insert_not_the_raw_read() {
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    write!(fq, "@r1\nAAAAAAAAAA\n+\n####IIIIII\n").unwrap();

    let no_crop = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "-q",
            "30",
            "-m",
            "arithmetic",
            "-t",
            "1",
        ])
        .assert()
        .success();
    assert_eq!(
        no_crop.get_output().stdout,
        b"",
        "sanity: the raw whole-read mean (24.8) must fail -q 30 with no trim applied"
    );

    let cropped = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "-H",
            "4",
            "-q",
            "30",
            "-m",
            "arithmetic",
            "-t",
            "1",
        ])
        .assert()
        .success();
    assert_eq!(
        cropped.get_output().stdout,
        b"@r1\nAAAAAA\n+\nIIIIII\n",
        "trimmed insert (mean 40) must survive -q 30 once the bad flank is cropped away first"
    );

    let guard = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "-H",
            "4",
            "-q",
            "45",
            "-m",
            "arithmetic",
            "-t",
            "1",
        ])
        .assert()
        .success();
    let guard_out = guard.get_output();
    assert_eq!(
        guard_out.stdout, b"",
        "guard: the insert's own mean (40) must still fail -q 45"
    );
    let stderr = String::from_utf8_lossy(&guard_out.stderr);
    assert!(
        stderr.contains("No reads survived"),
        "guard run must report the read as fully dropped: {stderr}"
    );
}

/// Surviving segments retain their produced index rather than being
/// renumbered after filtering.
#[test]
fn produced_index_naming_end_to_end() {
    let adapter = "GGGGTTTTGGGGTTTT"; // 16bp, G/T only -> no accidental match in the A flanks
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">mid\n{adapter}").unwrap();

    let mut fq = tempfile::NamedTempFile::new().unwrap();
    writeln!(fq, "@one_seg\n{}\n+\n{}", "A".repeat(20), "I".repeat(20)).unwrap();
    writeln!(
        fq,
        "@two_seg\n{}{adapter}{}\n+\n{}",
        "A".repeat(10),
        "A".repeat(10),
        "I".repeat(36)
    )
    .unwrap();
    writeln!(
        fq,
        "@gap_seg\n{}{adapter}{}\n+\n{}",
        "A".repeat(3),
        "A".repeat(10),
        "I".repeat(29)
    )
    .unwrap();

    let out = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "--adapter-error-rate",
            "0.1",
            "--adapter-end-size",
            "1",
            "-l",
            "5",
            "-t",
            "1",
        ])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();

    let expected = format!(
        "@one_seg\n{}\n+\n{}\n@two_seg_segment_1\n{}\n+\n{}\n@two_seg_segment_2\n{}\n+\n{}\n@gap_seg_segment_2\n{}\n+\n{}\n",
        "A".repeat(20),
        "I".repeat(20),
        "A".repeat(10),
        "I".repeat(10),
        "A".repeat(10),
        "I".repeat(10),
        "A".repeat(10),
        "I".repeat(10),
    );
    assert_eq!(
        String::from_utf8(out).unwrap(),
        expected,
        "exact survivor names/order: unsplit unchanged, both-survive suffixed 1/2, \
         first-filtered leaves a lone _segment_2"
    );
}

/// Accounting/summary end-to-end: a mix of clean reads, a chimera with one
/// sub-`-l` half, an all-adapter read (fully consumed by terminal trimming, so
/// `trim::apply` produces zero segments), an empty read, and a read whose sole
/// produced segment is itself too short. Exercises the three-way read-level
/// counter split through the real binary's rendered summary
/// (`obs.rs`'s `Summary:`/`Trimmed to nothing:`/`All segments filtered:`/
/// `Segments dropped:` lines), not just the `Counters`/`Stats` unit tests in
/// `src/workflow/mod.rs`. In particular it distinguishes the all-adapter/empty
/// reads (`reads_trimmed_to_nothing` — no segments produced at all) from the
/// short-only read (`reads_all_filtered` — one segment produced, then
/// filtered).
#[test]
fn accounting_summary_end_to_end() {
    let adapter = "GGGGTTTTGGGGTTTTGGGG"; // 20bp, G/T only
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    writeln!(fa, ">mid\n{adapter}").unwrap();

    let mut fq = tempfile::NamedTempFile::new().unwrap();
    // Two clean reads: no adapter present, comfortably above -l 5.
    writeln!(fq, "@clean1\n{}\n+\n{}", "A".repeat(20), "I".repeat(20)).unwrap();
    writeln!(fq, "@clean2\n{}\n+\n{}", "A".repeat(20), "I".repeat(20)).unwrap();
    // Chimera: interior adapter splits into a 3bp flank (< -l 5, filtered
    // TooShort) and a 15bp flank (survives) -> 1 read with output, 1 segment
    // written, 1 segment dropped.
    writeln!(
        fq,
        "@chimera\n{}{adapter}{}\n+\n{}",
        "A".repeat(3),
        "A".repeat(15),
        "I".repeat(38)
    )
    .unwrap();
    // All-adapter: the read IS the adapter -- terminal trimming consumes the
    // whole thing, so `trim::apply` produces zero segments (no segment-level
    // drop; the per-segment filter loop never runs) -> reads_trimmed_to_nothing.
    writeln!(fq, "@alladapter\n{adapter}\n+\n{}", "I".repeat(20)).unwrap();
    // Empty read: 0-length SEQ/QUAL, not silently skipped from accounting ->
    // also reads_trimmed_to_nothing.
    write!(fq, "@empty\n\n+\n\n").unwrap();
    // Short-only: no adapter present, no split -- trim::apply returns the
    // whole (untrimmed) 3bp read as its single produced segment, which the
    // length filter then rejects -> reads_all_filtered (distinct from
    // reads_trimmed_to_nothing: a segment WAS produced, just not one that
    // survived).
    writeln!(fq, "@shortonly\n{}\n+\n{}", "A".repeat(3), "I".repeat(3)).unwrap();

    let res = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "--adapter-error-rate",
            "0.1",
            "--adapter-end-size",
            "1",
            "-l",
            "5",
            "-t",
            "1",
        ])
        .assert()
        .success();
    let stderr = String::from_utf8_lossy(&res.get_output().stderr);

    // Read level: 6 input reads; 3 produced output (2 clean + the chimera's
    // surviving half); 2 trimmed to nothing (the fully-consumed adapter read,
    // the empty read); 1 with every segment filtered (shortonly) -- the
    // three-way invariant (`reads_with_output + reads_trimmed_to_nothing +
    // reads_all_filtered == input_reads`) holds via `snapshot`'s debug assert.
    assert!(
        stderr.contains("Summary: 6 input reads, 3 output reads"),
        "input_reads=6, segments written=3: {stderr}"
    );
    assert!(
        stderr.contains("Trimmed to nothing: 2 input reads produced no segments at all"),
        "reads_trimmed_to_nothing=2 (alladapter + empty): {stderr}"
    );
    assert!(
        stderr.contains("All segments filtered: 1 input reads had every produced segment filtered"),
        "reads_all_filtered=1 (shortonly): {stderr}"
    );
    assert!(
        stderr.contains("Segments dropped: 2 (2 too short)"),
        "segments_dropped_short=2 (the chimera's short half + shortonly's sole segment): {stderr}"
    );
}

/// Quality splitting emits both high-quality pieces before length filtering.
/// The surviving second piece retains the `_segment_2` suffix.
#[test]
fn qual_split_emits_short_pieces_for_post_trim_filter_to_own() {
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    writeln!(fq, "@r1\nAAAATAAAAAA\n+\nIIII#IIIIII").unwrap();

    let res = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--qual-split",
            "10",
            "--qual-split-window",
            "1",
            "-l",
            "5",
        ])
        .assert()
        .success();
    let out = res.get_output();
    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert!(
        stdout.contains("r1_segment_2"),
        "the surviving AAAAAA piece [5,11) must be named _segment_2 (produced \
         as the 2nd of 2 pieces, even though the 1st was filtered): {stdout}"
    );
    assert!(
        !stdout.contains("@r1\n"),
        "the read must not appear unsuffixed -- it was split into 2 produced \
         pieces: {stdout}"
    );
    assert!(
        !stdout.contains("r1_segment_1"),
        "the short AAAA piece [0,4) is filtered TooShort by -l 5, not written: {stdout}"
    );
    assert!(
        stderr.contains("Segments dropped: 1 (1 too short)"),
        "the short high-quality piece must be a real produced-then-filtered \
         segment, bumping segments_dropped_short by exactly 1: {stderr}"
    );
}

/// Same read, `-l 4`: both produced pieces (len 4 and len 6) now clear the
/// length filter, so both survive and are numbered 1/2.
#[test]
fn qual_split_both_pieces_survive_at_lower_length_floor() {
    let mut fq = tempfile::NamedTempFile::new().unwrap();
    writeln!(fq, "@r1\nAAAATAAAAAA\n+\nIIII#IIIIII").unwrap();

    let res = Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            fq.path().to_str().unwrap(),
            "--qual-split",
            "10",
            "--qual-split-window",
            "1",
            "-l",
            "4",
        ])
        .assert()
        .success();
    let stdout = String::from_utf8_lossy(&res.get_output().stdout);

    assert!(
        stdout.contains("r1_segment_1"),
        "AAAA [0,4) now clears -l 4: {stdout}"
    );
    assert!(
        stdout.contains("r1_segment_2"),
        "AAAAAA [5,11) still survives: {stdout}"
    );
}
