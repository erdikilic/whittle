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
        .stdout(predicates::str::contains("--adapter-ends-only"));
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

// Build a FASTQ where every read starts with adapter A (present) and none
// contains adapter B (absent). With a *custom* FASTA, detection is disabled
// outright (Change 1: a curated FASTA should always be searched in full), so
// neither adapter is ever reduced -- both stay active, no "Adapter presence"
// log appears, and trimming on the present adapter still works exactly as
// without detection.
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

// Custom --adapter-fasta disables detection outright (no small-sample branch
// to hit at all), so this now exercises the small-sample path with
// --adapter-preset instead, which still runs detection for < 100 reads.
// Detection is opt-in now (default --adapter-sample is 0, i.e. off), so this
// must explicitly opt in to exercise the small-sample skip branch.
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

// The key correctness guarantee (Phase 1.5 spec): detection only ever *drops*
// adapters that don't act, so trimming a present adapter with detection ON
// must be byte-identical to trimming with detection OFF (`--adapter-sample
// 0`, i.e. the full set). Custom `--adapter-fasta` now disables detection
// outright (Change 1), so this equivalence can only be exercised against a
// preset. Uses `--adapter-preset ont` (124 real catalog entries) with the
// real LSK114 front ligation adapter present in every read; detection reduces
// 124 -> a handful (confirmed empirically: kept 4 -- LSK114_front/_rear and
// LSK109_front/_rear, whose sequences overlap the front adapter's tail).
//
// The insert is 300bp (not the usual 40bp): with a short ~67bp read, the
// catalog's paired front/rear entries (Y-adapter chemistry means a "rear"
// entry's sequence is a near-reverse-complement of the "front" one) can both
// match within the same read when `--adapter-end-size` (default 150) spans
// the whole thing, consuming it entirely (empirically: 0 output reads for
// *both* detection on and off -- technically "byte-identical" but violates
// the "non-empty" sanity check and isn't a meaningful comparison). A 300bp
// insert pushes the tail well outside the 150bp end-zone, so only the front
// adapter acts and the insert survives trimming, for both settings.
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
        // order, for `threads > 1` -- see `pipeline::fastq::run`). That
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

// Owner reproduction (a): 10000 clean reads followed by 10 adapted ones lost
// ALL adapter trimming on the tail, because presence detection sampled
// exactly the (adapter-free) first `adapter_sample` reads -- default 10000,
// matching this fixture's clean-read count precisely -- kept zero adapters,
// and reduced the active set to nothing, silently disabling trimming for
// every read that followed, including the 10 adapted ones. Fixed by Change 1:
// a custom --adapter-fasta is a curated set that should always be searched in
// full, so detection is now forced off unconditionally whenever a FASTA is
// given, regardless of --adapter-sample's value. Detection is opt-in now
// (default --adapter-sample is 0), so this passes --adapter-sample 10000
// explicitly -- a value that would otherwise enable detection -- to prove the
// fasta override still holds "regardless of --adapter-sample's value" rather
// than merely benefiting from the new off-by-default behavior.
//
// Confirmed RED under the pre-fix code (reverting Change 1 only): stderr
// logged "Adapter presence: sampled 10000 reads, kept 0 of 1 adapters" and
// the summary showed "10,010 input reads, 10,010 output reads ... 100.0%
// kept" -- i.e. the 10 adapted reads passed through completely untouched,
// their adapter+insert line appearing verbatim in stdout. GREEN after the
// fix: the banner reports "sample off", and the 10 adapted reads are
// trimmed down to their bare insert.
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

// Owner reproduction (b): a preset run where the sampled prefix happens to
// contain zero adapters (e.g. a run of clean reads ahead of the adapted
// ones) used to reduce the active set to *nothing*, silently disabling
// trimming for the rest of the file -- including the adapted reads outside
// the sample. Fixed by Change 2: when detection keeps zero adapters, fall
// back to the full configured set (with a WARN) instead of reducing to an
// empty one.
//
// Fixture: 100 clean reads (the entire sample, since --adapter-sample 100 ==
// MIN_SAMPLE_FOR_DETECTION, so detection runs rather than skipping), then 10
// reads carrying the real LSK114 front adapter. The clean reads are pure C,
// which cannot spuriously match any preset adapter (all are mixed-base)
// within the default edit budget -- confirmed empirically: detection keeps 0
// on the clean-only sample.
//
// Confirmed RED under the pre-fix code (reverting Change 2 only): stderr
// logged "Adapter presence: sampled 100 reads, kept 0 of 124 adapters" (no
// fallback), and the summary showed "110 input reads, 110 output reads ...
// 100.0% kept" -- the 10 adapted reads' untrimmed line appeared verbatim in
// stdout. GREEN after the fix: stderr carries the fallback WARN and the 10
// adapted reads are trimmed (dropped entirely here, since the short fixture
// read is fully consumed by the catalog's paired front/rear entries within
// the default 150bp end-zone -- see the comment on
// `detection_output_equals_full_set_for_present_adapter` for why).
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

// --adapter-infer-only + --adapter-fasta is allowed (unlike --adapter-infer,
// which rejects a FASTA outright), but v1 descopes cross-naming discovered
// sequences against the user's FASTA -- naming stays catalog-only. This just
// checks the one informational line fires so a user combining the two flags
// isn't left assuming the FASTA did something.
#[test]
fn infer_only_with_fasta_notes_naming_is_catalog_only() {
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
            "--adapter-infer-only",
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
        ])
        .assert()
        .success()
        .stderr(predicates::str::contains("catalog only"));
}

// --- ab-initio inference wiring (Task 11) -------------------------------
//
// Fixtures below plant an EXACT copy (no injected error -- error-tolerant
// recovery is already covered by `discover_recovers_planted_adapter_under_error`
// in src/adapter/infer.rs) of a real catalog-neighborhood adapter at the 5'
// end of every read, followed by a deterministic splitmix64-mixed genomic
// tail distinct per read index.
//
// IMPORTANT: a naive `(a*i + b*j) % 4` background generator is periodic
// (linear in `j` mod 4) and collapses into a phase-rotated ACGT tandem
// repeat -- a spurious, low-complexity-but-not-homopolymer signal that the
// k-mer discoverer picks up as a fake "adapter" of its own, breaking these
// tests. The splitmix64 bit-mix below is the same fixture pattern
// `src/adapter/infer.rs`'s own `discover_*` unit tests use, and does not
// have that defect.

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
fn infer_only_prints_and_does_not_trim() {
    let dir = tempfile::tempdir().unwrap();
    let fq = write_adapted_fastq(dir.path(), 500);
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd.args([
        "-i",
        fq.to_str().unwrap(),
        "--adapter-infer-only",
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

    // Genuine trimming check (not just "the adapter substring is gone"):
    // every surviving record's sequence must be an exact SUFFIX of the read
    // whittle actually read in (reconstructed independently via
    // `full_read_seq`, not re-derived from the output), and the amount cut
    // off the front must land in a sane window around the 28bp planted
    // adapter's length -- proving real per-read adapter-shaped trimming, not
    // a no-op, a fixed head-crop, or a whole-read wipe.
    let mut lines = trimmed.lines();
    let mut n_records = 0;
    while let Some(header) = lines.next() {
        assert!(header.starts_with("@r"), "unexpected header: {header}");
        let idx: usize = header[2..].parse().expect("header must be @r<index>");
        let seq_line = lines.next().expect("sequence line");
        let _plus = lines.next().expect("plus line");
        let _qual = lines.next().expect("quality line");

        let original = full_read_seq(idx);
        let cut = original.len() - seq_line.len();
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
