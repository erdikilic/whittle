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
// which rejects a FASTA outright): naming now covers the built-in catalog
// PLUS the user's FASTA (FU3 -- see `infer::discover`'s `name_refs`). This
// just checks the informational line reflects that (not the old "catalog
// only" wording), so a user combining the two flags isn't left assuming the
// FASTA did nothing. The actual cross-naming is proven end-to-end by
// `infer_only_cross_names_against_user_fasta` below.
#[test]
fn infer_only_with_fasta_notes_naming_includes_fasta() {
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
        .stderr(predicates::str::contains("plus your FASTA's adapters"));
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

// Bug 3 regression: CLI help promises `--adapter-infer-only` prints
// "sequences + support + catalog names", but pre-fix `log_discovered` only
// ever logged the sequence at `debug!` -- invisible at the default INFO
// level -- so plain `--adapter-infer-only` stdout was completely empty.
// Report-only must print each discovered adapter to stdout as a FASTA
// record (header + sequence line), so a user can redirect stdout straight
// into an adapter FASTA of their own.
#[test]
fn infer_only_prints_sequence_to_stdout() {
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
}

// FU3: report-only cross-names discovered adapters against the ONT catalog
// UNION the user's --adapter-fasta, not the catalog alone.
//
// `PLANTED_ADAPTER` is byte-identical to the catalog's own `LSK109_front`
// entry (see `src/adapter/ont_catalog.tsv`), so a discovered consensus that
// reconstructs it scores identically (same bytes compared, same edit-distance
// search) against BOTH the catalog entry and our own FASTA-supplied copy --
// an exact tie in `name_against`'s percent-identity, broken by its
// alphabetical (name asc) tie-break. The FASTA header is prefixed `AAA_` so
// it sorts before `LSK109_front` and therefore deterministically wins that
// tie, becoming `name_hits[0]` -- the only hit `log_discovered` prints -- no
// matter how `discover` actually reconstructs the consensus. That makes this
// a genuine proof that naming consulted the user's FASTA (not merely that it
// also happened to match the catalog): if `discover` still only checked the
// built-in catalog (pre-fix), the log would show `LSK109_front` instead and
// this assertion would fail.
#[test]
fn infer_only_cross_names_against_user_fasta() {
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
        "--adapter-infer-only",
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
        // `checked_sub` (M5): a clear panic message if output ever somehow
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

// --- Task 11 review-fix regressions: `--adapter-infer-only` must NEVER
// write or touch output -------------------------------------------------

/// HIGH bug: the too-few-reads branch of the infer path used to return
/// `Ok(Some(chain(sample, records)))` unconditionally, so `--adapter-infer-only`
/// on an undersized input warned, then still dispatched and wrote the full
/// (untrimmed) input back out through `-o`. `ReportOnly` must never write
/// output, no matter whether discovery itself ran or was skipped for too few
/// reads.
#[test]
fn infer_only_tiny_input_writes_no_output() {
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
        "--adapter-infer-only",
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

/// LOW bug: the FASTQ dispatch arm used to construct the output writer (a
/// truncating `File::create`) BEFORE the buffer-and-decide seam
/// (`maybe_reduce_adapters`), so `--adapter-infer-only -o existing.txt`
/// truncated `existing.txt` to zero bytes even though report-only writes no
/// records at all. The writer must only be created after the seam has had
/// its chance to return the "stop now, no dispatch" signal.
#[test]
fn infer_only_does_not_clobber_output_file() {
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
        "--adapter-infer-only",
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

// --- Task 12: determinism ------------------------------------------------

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
        "--adapter-infer-only",
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

// --- Bug 4: the binary-stdout TTY guard must not run before report-only's
// early exit -------------------------------------------------------------
//
// `guard_stdout_binary` used to run during output setup, BEFORE
// `maybe_reduce_adapters`'s `ReportOnly` early-exit, so `whittle -i
// reads.bam --adapter-infer-only` on a terminal was refused ("would write
// BAM to terminal") even though report-only writes no BAM at all -- it only
// ever prints a small FASTA text summary to stdout. The fix makes
// `guard_stdout_binary` itself exempt `AdapterInfer::ReportOnly`.
//
// This process's stdout is piped (not a terminal) under the test harness,
// so the guard's TTY check would never actually fire here regardless of the
// fix -- the real regression can only be exercised on a real terminal (see
// the fix's doc comment in `src/lib.rs` for the code-inspection argument).
// This is the "cheap assertion" fallback: a sanity check that report-only
// with a BAM input (which resolves `out_fmt` to `Bam`, the same format the
// guard used to hard-error on) and no `-o` still completes successfully end
// to end, rather than tripping over some other BAM-specific issue in the
// report-only path.
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
fn infer_only_on_bam_input_with_piped_stdout_succeeds() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("in.bam");
    write_minimal_ubam(&in_path, 500);

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "-i",
            in_path.to_str().unwrap(),
            "--adapter-infer-only",
            "-t",
            "1",
        ])
        .assert()
        .success();
}

// --- Task 2: trim-then-filter reorder end-to-end regressions --------------
//
// Task 1 moved `filter::check` from pre-trim (whole raw read) to post-trim
// (each produced segment), and split read-level accounting from segment-level
// drop counts. These pin that behavior change through the compiled binary
// end-to-end, on top of the internal `run_fastq_seq` unit tests already in
// `src/workflow/fastq.rs`.

/// The headline behavior change: quality is judged on the TRIMMED insert, not
/// the raw read. `r1` is 4 low-quality bases (phred 2, `'#'`) then 6
/// high-quality bases (phred 40, `'I'`); raw arithmetic mean = (2*4 + 40*6) /
/// 10 = 24.8.
///
/// Run 1 (sanity): no head-crop, `-q 30 -m arithmetic` -> the filter runs on
/// the whole raw read (its only produced segment, since there's no adapter/
/// crop to carve it up) and 24.8 < 30 fails -> no output. This reproduces
/// what the OLD pre-reorder (filter-before-trim) code measured, for real,
/// through the binary.
/// Run 2 (the fix): `-H 4` crops the bad flank away BEFORE filtering, SAME
/// `-q 30` -> the surviving 6-base insert's own mean (40) passes -> written.
/// Run 3 (guard): same crop, but `-q 45` (above the insert's own mean of 40)
/// -> dropped again, proving the insert's own quality is still enforced, not
/// that trimming silently disabled filtering.
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

/// Naming (produced-index, "option b"): a survivor keeps the index it was
/// PRODUCED at, not renumbered by how many segments survived. Three reads in
/// one run, exact output asserted: `one_seg` has no adapter match (produced
/// == 1) -> name unchanged; `two_seg` splits into two 10bp flanks that both
/// clear `-l 5` -> `_segment_1`/`_segment_2`; `gap_seg` splits into a 3bp
/// flank (< `-l 5`, filtered `TooShort`) and a 10bp flank (survives) -> the
/// survivor is a LONE `_segment_2`, signalling the split even though only the
/// second piece made it to output.
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
/// sub-`-l` half, an all-adapter read (fully consumed by terminal trimming,
/// so `trim::apply` produces zero segments), and an empty read. Exercises the
/// two-level counter model (Task 1) through the real binary's rendered
/// summary (`obs.rs`'s `Summary:`/`No output:`/`Segments dropped:` lines), not
/// just the `Counters`/`Stats` unit tests in `src/workflow/mod.rs`.
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
    // drop; the per-segment filter loop never runs).
    writeln!(fq, "@alladapter\n{adapter}\n+\n{}", "I".repeat(20)).unwrap();
    // Empty read: 0-length SEQ/QUAL, not silently skipped from accounting.
    write!(fq, "@empty\n\n+\n\n").unwrap();

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

    // Read level: 5 input reads; 3 produced output (2 clean + the chimera's
    // surviving half); 2 produced none (the fully-consumed adapter read, the
    // empty read) -- the two-level invariant (`reads_with_output +
    // reads_no_output == input_reads`) holds via `snapshot`'s debug assert.
    assert!(
        stderr.contains("Summary: 5 input reads, 3 output reads"),
        "input_reads=5, segments written=3: {stderr}"
    );
    assert!(
        stderr.contains("No output: 2 input reads produced no surviving segment"),
        "reads_no_output=2 (alladapter + empty): {stderr}"
    );
    assert!(
        stderr.contains("Segments dropped: 1 (1 too short)"),
        "segments_dropped_short=1 (the chimera's short half): {stderr}"
    );
}

/// Owner-reported regression: `--qual-split` must emit *every* high-quality
/// piece (even ones shorter than `-l`) and let the post-trim filter drop the
/// short ones -- it must not suppress them inside the trim stage itself.
/// `AAAATAAAAAA` / `IIII#IIIIII` with `--qual-split-window 1` splits on the
/// single low-quality base at index 4 into two produced pieces: `AAAA` at
/// [0,4) and `AAAAAA` at [5,11). With `-l 5` the first (len 4) is filtered
/// TooShort and the second (len 6) survives -- as the *second* produced
/// segment, so it must be named `_segment_2`, not `_segment_1` (which would
/// happen if the short piece were dropped pre-filter, shrinking the produced
/// count to 1).
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
