# chopping v1 — Plan 1: FASTQ trimmer core — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A complete, well-tested FASTQ(.gz) trimmer/splitter with chopper-parity filtering and trimming, and a shared record/pipeline/CLI foundation that Plan 2 (uBAM + MM/ML reconstruction) builds on.

**Architecture:** Approach A from the spec — filters and trimmers are pure functions over `(seq: &[u8], phred: &[u8])` slices, returning boolean / `Vec<(start,end)>` intervals; a format-specific reconstruction step slices the record per interval and writes it. Internal quality is stored as **raw Phred scores** (u8, 0-based); the FASTQ I/O layer converts to/from ASCII+33 at the boundary. Trimming never invents sequence — it only returns intervals — which is what makes Plan 2's tag reconstruction tractable.

**Tech Stack:** Rust 2024, `seq_io` 0.3.4 (FASTQ), `flate2` 1 (gz), `clap` 4 derive, `rayon` + `crossbeam-channel` (parallelism), `anyhow` (binary) + `thiserror` (library), `assert_cmd` + `predicates` (CLI tests).

Design source: the spec at `docs/superpowers/specs/2026-07-03-chopping-tag-aware-trimmer-design.md`. Ported algorithms come from `clone/chopper/src/{trimmers.rs,records.rs,main.rs}` (MIT).

## Global Constraints

- Rust **2024 edition**. `rust-version` not pinned.
- Internal quality representation is **raw Phred scores** (`u8`, e.g. Q10 == `10u8`), never ASCII. FASTQ reader subtracts 33; FASTQ writer adds 33.
- Dependencies (exact): `seq_io = "0.3.4"`, `flate2 = { version = "1", features = ["zlib-ng"] }`, `clap = { version = "4", features = ["derive"] }`, `rayon = "1"`, `crossbeam-channel = "0.5"`, `anyhow = "1"`, `thiserror = "2"`; dev: `assert_cmd = "2"`, `predicates = "3"`.
- Release profile: `lto = "fat"`, `codegen-units = 1` (matches chopper — CPU-bound workload).
- The three quality trim ops (`--trim-qual`, `--best-segment`, `--split-qual`) are **mutually exclusive**; supplying two is a hard error. Fixed crop (`--head-crop`/`--tail-crop`) composes and is applied first.
- Split segments and all output segments obey `--min-length`.
- `--qual-mode` default is `mean` (error-probability mean); it governs the `-q/-Q` read filter only.
- Every task ends green: `cargo test` passes and `cargo clippy -- -D warnings` is clean before commit.

## File Structure

```
Cargo.toml                     deps, release profile
src/
  main.rs                      clap entry -> chopping::run(cfg)
  lib.rs                       pub mods + run(cfg) orchestrator
  qual.rs                      Phred conversions + read-quality metrics + QualMode
  filter.rs                    FilterConfig + passes(seq, phred, cfg) -> bool
  trim/
    mod.rs                     TrimPlan, QualityOp, apply(...) -> Vec<(usize,usize)>
    strategies.rs              4 ported chopper algorithms over &[u8] phred
  record.rs                    ReadRecord { name, seq, qual } (qual = raw phred)
  io/
    mod.rs                     Format enum, detect_input, resolve_output
    fastq.rs                   seq_io reader (+ gz) -> ReadRecord; manual FASTQ writer
  pipeline.rs                  run_fastq(records, writer, cfg) single-thread + parallel
  cli.rs                       Cli (clap derive) -> Config, with validation
  config.rs                    Config { io, filter, trim, threads } shared struct
tests/
  cli.rs                       assert_cmd end-to-end + mutual-exclusion errors
test-data/
  basic.fastq                  small hand-built fixture
```

Plan 2 will add `src/mods/`, `src/io/bam.rs`, and extend `pipeline.rs` + `cli.rs`; nothing in Plan 1 should hard-code "FASTQ-only" in a way that blocks that (e.g. `Config` and `run()` dispatch on `Format`).

---

### Task 1: Project scaffold

**Files:**
- Create: `Cargo.toml`, `src/main.rs`, `src/lib.rs`

**Interfaces:**
- Produces: crate `chopping` (lib + bin); `pub fn chopping::run(cfg: Config) -> anyhow::Result<()>` (stub for now); binary calls it.

- [ ] **Step 1: Write `Cargo.toml`**

```toml
[package]
name = "chopping"
version = "0.1.0"
edition = "2024"
license = "MIT"
description = "Tag-aware long-read trimmer (ONT/PacBio): FASTQ/uBAM trimming + splitting with MM/ML reconstruction"

[[bin]]
name = "chopping"
path = "src/main.rs"

[lib]
name = "chopping"
path = "src/lib.rs"

[dependencies]
seq_io = "0.3.4"
flate2 = { version = "1", features = ["zlib-ng"], default-features = false }
clap = { version = "4", features = ["derive"] }
rayon = "1"
crossbeam-channel = "0.5"
anyhow = "1"
thiserror = "2"

[dev-dependencies]
assert_cmd = "2"
predicates = "3"

[profile.release]
lto = "fat"
codegen-units = 1
```

- [ ] **Step 2: Write `src/lib.rs`**

```rust
pub mod cli;
pub mod config;
pub mod filter;
pub mod io;
pub mod pipeline;
pub mod qual;
pub mod record;
pub mod trim;

pub use config::Config;

/// Top-level entry point: dispatch on the resolved input/output formats and run
/// the matching pipeline. Plan 1 implements only the FASTQ path; Plan 2 adds BAM.
pub fn run(_cfg: Config) -> anyhow::Result<()> {
    anyhow::bail!("not yet implemented")
}
```

- [ ] **Step 3: Write `src/main.rs`**

```rust
fn main() -> anyhow::Result<()> {
    let cfg = chopping::cli::parse()?;
    chopping::run(cfg)
}
```

- [ ] **Step 4: Create empty module files so it compiles**

Create each of `src/qual.rs`, `src/filter.rs`, `src/record.rs`, `src/pipeline.rs`, `src/config.rs`, `src/cli.rs`, `src/io/mod.rs`, `src/trim/mod.rs` with a single line `// implemented in a later task` — EXCEPT the ones a later step in *this* task fills. For Task 1, create `src/config.rs` and `src/cli.rs` with the minimal stubs below; the rest get real content in their tasks. To keep Task 1 compiling, temporarily stub the not-yet-written modules.

`src/config.rs` (minimal, expanded in Task 8):
```rust
#[derive(Debug, Clone)]
pub struct Config;
```

`src/cli.rs` (minimal, expanded in Task 8):
```rust
use crate::config::Config;

pub fn parse() -> anyhow::Result<Config> {
    Ok(Config)
}
```

Stub the remaining modules referenced by `lib.rs` with empty files containing `// filled in a later task`.

- [ ] **Step 5: Verify it builds and the empty test set passes**

Run: `cargo build && cargo test`
Expected: build succeeds; `test result: ok. 0 passed`.

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml Cargo.lock src/
git commit -m "feat: scaffold chopping crate (lib + bin, deps, run stub)"
```

---

### Task 2: Quality metrics (`src/qual.rs`)

**Files:**
- Create/replace: `src/qual.rs`

**Interfaces:**
- Produces:
  - `pub fn phred_to_prob(q: u8) -> f64`
  - `pub fn mean_prob_q(phred: &[u8]) -> f64`
  - `pub fn mean_arith_q(phred: &[u8]) -> f64`
  - `pub fn median_q(phred: &[u8]) -> f64`
  - `pub enum QualMode { Mean, Arithmetic, Median }` (derives `Debug, Clone, Copy, PartialEq, Eq`)
  - `pub fn read_quality(phred: &[u8], mode: QualMode) -> f64`

- [ ] **Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prob_matches_phred_definition() {
        assert!((phred_to_prob(20) - 0.01).abs() < 1e-12);
        assert!((phred_to_prob(30) - 0.001).abs() < 1e-12);
    }

    // Values ported from chopper's ave_qual test, but inputs are RAW phred (no +33).
    #[test]
    fn mean_prob_q_matches_chopper() {
        assert!((mean_prob_q(&[10]) - 10.0).abs() < 1e-9);
        assert!((mean_prob_q(&[10, 11, 12]) - 10.923583702678473).abs() < 1e-9);
        assert!((mean_prob_q(&[10, 11, 12, 20, 30, 40, 50]) - 14.408827647036087).abs() < 1e-9);
    }

    #[test]
    fn empty_is_zero() {
        assert_eq!(mean_prob_q(&[]), 0.0);
        assert_eq!(mean_arith_q(&[]), 0.0);
        assert_eq!(median_q(&[]), 0.0);
    }

    #[test]
    fn arithmetic_and_median() {
        assert!((mean_arith_q(&[10, 20, 30]) - 20.0).abs() < 1e-9);
        assert!((median_q(&[10, 20, 30]) - 20.0).abs() < 1e-9);
        // even count -> average of the two middle values
        assert!((median_q(&[10, 20, 30, 40]) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn read_quality_dispatches() {
        assert_eq!(read_quality(&[10, 20, 30], QualMode::Arithmetic), 20.0);
        assert_eq!(read_quality(&[10, 20, 30], QualMode::Median), 20.0);
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test qual::`
Expected: FAIL (functions/types undefined).

- [ ] **Step 3: Implement `src/qual.rs`**

```rust
use std::sync::LazyLock;

/// Precomputed 10^(-q/10) for every possible Phred byte. Sizing to the full u8
/// range means any quality byte indexes safely (ported from chopper's PHRED_LUT).
static PHRED_LUT: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut lut = [0.0f64; 256];
    for (i, v) in lut.iter_mut().enumerate() {
        *v = 10_f64.powf((i as f64) / -10.0);
    }
    lut
});

#[inline(always)]
pub fn phred_to_prob(q: u8) -> f64 {
    PHRED_LUT[q as usize]
}

/// Error-probability mean quality: the ONT-standard "read Q" (chopper's ave_qual).
pub fn mean_prob_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let sum: f64 = phred.iter().map(|&q| phred_to_prob(q)).sum();
    (sum / phred.len() as f64).log10() * -10.0
}

/// Plain arithmetic mean of the Phred integers.
pub fn mean_arith_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let sum: u64 = phred.iter().map(|&q| q as u64).sum();
    sum as f64 / phred.len() as f64
}

/// Median Phred via a 256-bucket histogram (O(n), no sort/alloc of the input).
pub fn median_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let mut hist = [0usize; 256];
    for &q in phred {
        hist[q as usize] += 1;
    }
    let n = phred.len();
    let mid = n / 2;
    // Walk buckets accumulating counts; find value(s) at the median rank(s).
    let value_at = |target: usize| -> usize {
        let mut cum = 0usize;
        for (v, &c) in hist.iter().enumerate() {
            cum += c;
            if cum > target {
                return v;
            }
        }
        255
    };
    if n % 2 == 1 {
        value_at(mid) as f64
    } else {
        (value_at(mid - 1) + value_at(mid)) as f64 / 2.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QualMode {
    Mean,
    Arithmetic,
    Median,
}

pub fn read_quality(phred: &[u8], mode: QualMode) -> f64 {
    match mode {
        QualMode::Mean => mean_prob_q(phred),
        QualMode::Arithmetic => mean_arith_q(phred),
        QualMode::Median => median_q(phred),
    }
}
```

Append the test module from Step 1 to this file.

- [ ] **Step 4: Run to verify pass**

Run: `cargo test qual::`
Expected: PASS (all 5 tests).

- [ ] **Step 5: Commit**

```bash
git add src/qual.rs
git commit -m "feat: quality metrics (prob/arith/median) with chopper-parity values"
```

---

### Task 3: Filters (`src/filter.rs`)

**Files:**
- Replace: `src/filter.rs`

**Interfaces:**
- Consumes: `crate::qual::{QualMode, read_quality}`
- Produces:
  - `pub struct FilterConfig { pub min_length: usize, pub max_length: usize, pub min_qual: f64, pub max_qual: f64, pub min_gc: Option<f64>, pub max_gc: Option<f64>, pub qual_mode: QualMode }`
  - `pub fn passes(seq: &[u8], phred: &[u8], cfg: &FilterConfig) -> bool`
  - `pub fn gc_fraction(seq: &[u8]) -> f64`

- [ ] **Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::qual::QualMode;

    fn base() -> FilterConfig {
        FilterConfig {
            min_length: 1,
            max_length: usize::MAX,
            min_qual: 0.0,
            max_qual: 1000.0,
            min_gc: None,
            max_gc: None,
            qual_mode: QualMode::Mean,
        }
    }

    #[test]
    fn length_bounds() {
        let mut c = base();
        c.min_length = 4;
        c.max_length = 8;
        assert!(!passes(b"ATG", &[30, 30, 30], &c)); // too short
        assert!(passes(b"ATGCG", &[30; 5], &c));
        assert!(!passes(b"ATGCGATGC", &[30; 9], &c)); // too long
    }

    #[test]
    fn quality_bound_uses_mode() {
        let mut c = base();
        c.min_qual = 15.0;
        // arithmetic mean of [10,20] = 15.0 -> passes at threshold
        c.qual_mode = QualMode::Arithmetic;
        assert!(passes(b"AT", &[10, 20], &c));
        // prob-mean of [10,20] < 15 -> fails
        c.qual_mode = QualMode::Mean;
        assert!(!passes(b"AT", &[10, 20], &c));
    }

    #[test]
    fn gc_fraction_and_filter() {
        assert!((gc_fraction(b"GGCC") - 1.0).abs() < 1e-12);
        assert!((gc_fraction(b"ATAT") - 0.0).abs() < 1e-12);
        let mut c = base();
        c.min_gc = Some(0.4);
        c.max_gc = Some(0.6);
        assert!(passes(b"ATGC", &[30; 4], &c)); // 0.5
        assert!(!passes(b"AAAT", &[30; 4], &c)); // 0.0
    }

    #[test]
    fn empty_seq_rejected() {
        assert!(!passes(b"", &[], &base()));
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test filter::`
Expected: FAIL.

- [ ] **Step 3: Implement `src/filter.rs`**

```rust
use crate::qual::{QualMode, read_quality};

#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub min_length: usize,
    pub max_length: usize,
    pub min_qual: f64,
    pub max_qual: f64,
    pub min_gc: Option<f64>,
    pub max_gc: Option<f64>,
    pub qual_mode: QualMode,
}

pub fn gc_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c'))
        .count();
    gc as f64 / seq.len() as f64
}

/// Cheapest-first, short-circuiting. Empty reads never pass.
pub fn passes(seq: &[u8], phred: &[u8], cfg: &FilterConfig) -> bool {
    let len = seq.len();
    if len == 0 || len < cfg.min_length || len > cfg.max_length {
        return false;
    }
    if cfg.min_qual > 0.0 || cfg.max_qual < 1000.0 {
        let q = read_quality(phred, cfg.qual_mode);
        if q < cfg.min_qual || q > cfg.max_qual {
            return false;
        }
    }
    if cfg.min_gc.is_some() || cfg.max_gc.is_some() {
        let gc = gc_fraction(seq);
        if gc < cfg.min_gc.unwrap_or(0.0) || gc > cfg.max_gc.unwrap_or(1.0) {
            return false;
        }
    }
    true
}
```

Append Step 1's tests.

- [ ] **Step 4: Run to verify pass**

Run: `cargo test filter::`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/filter.rs
git commit -m "feat: length/quality/GC filters with selectable quality mode"
```

---

### Task 4: Trimmers (`src/trim/strategies.rs` + `src/trim/mod.rs`)

**Files:**
- Create: `src/trim/strategies.rs`
- Replace: `src/trim/mod.rs`

**Interfaces:**
- Consumes: `crate::qual::phred_to_prob`
- Produces (in `strategies.rs`, all over raw phred slices, returning `Vec<(usize, usize)>` half-open intervals):
  - `pub fn trim_by_quality(phred: &[u8], cutoff: u8) -> Vec<(usize, usize)>`
  - `pub fn best_segment(phred: &[u8], cutoff_q: u8) -> Vec<(usize, usize)>`
  - `pub fn split_low_quality(phred: &[u8], cutoff: u8, min_length: usize, window: usize) -> Vec<(usize, usize)>`
- Produces (in `mod.rs`):
  - `pub enum QualityOp { TrimQual(u8), BestSegment(u8), Split { cutoff: u8, window: usize } }` (derive `Debug, Clone, PartialEq, Eq`)
  - `pub struct TrimPlan { pub head: usize, pub tail: usize, pub quality: Option<QualityOp> }`
  - `pub fn apply(seq_len: usize, phred: &[u8], plan: &TrimPlan, min_length: usize) -> Vec<(usize, usize)>`

The three ported algorithms mirror `clone/chopper/src/trimmers.rs` exactly, except inputs are **raw phred** (chopper did `qual[i] - 33`; we compare `phred[i]` directly), and `best_segment` takes a Q score and converts it to a probability cutoff internally (chopper converted at the call site).

- [ ] **Step 1: Write failing tests for `strategies.rs`**

These reuse chopper's `get_reads()` fixtures, converting each ASCII quality byte to raw phred (`b - 33`). Expected intervals are copied verbatim from chopper's passing tests.

```rust
#[cfg(test)]
mod tests {
    use super::*;

    // (seq, ascii_qual) — identical bytes to chopper's trimmers.rs get_reads().
    fn reads() -> [(Vec<u8>, Vec<u8>); 6] {
        let raw = [
            (b"AAAAAAAAAAAAAAATTTAA".to_vec(), b"&#3-G27C:(@G7B55+C4I".to_vec()),
            (b"TTTTTTTTTTTTTTTTTTTT".to_vec(), b"77%'24)FAF9@=94'%054".to_vec()),
            (b"AAAAAAAAAAAAAAATTTTA".to_vec(), b"'8$-BF2!C;+59->H@91#".to_vec()),
            (b"AAAAAAAAAAAAAAAAAAAA".to_vec(), b"%,42$CH*#0+0C6=0,*6/".to_vec()),
            (b"AAAAAAAAAAAAAAAAAAAT".to_vec(), b"-------------------J".to_vec()),
            (b"TAAAAAAAAAAAAAAAAAAA".to_vec(), b"I-------------------".to_vec()),
        ];
        raw.map(|(s, q)| (s, q.iter().map(|&b| b - 33).collect()))
    }

    #[test]
    fn trim_by_quality_matches_chopper() {
        let expected: [(u8, Vec<(usize, usize)>); 6] = [
            (20, vec![(4, 20)]),
            (7, vec![(0, 20)]),
            (15, vec![(1, 19)]),
            (40, vec![]),
            (40, vec![(19, 20)]),
            (40, vec![(0, 1)]),
        ];
        for ((cutoff, want), (_, phred)) in expected.iter().zip(reads()) {
            assert_eq!(trim_by_quality(&phred, *cutoff), *want);
        }
    }

    #[test]
    fn best_segment_matches_chopper() {
        // chopper cutoffs were probabilities; the equivalent Q scores are:
        // 0.01=Q20, 0.199..=Q7, 0.0316..=Q15, 0.0001=Q40.
        let expected: [(u8, Vec<(usize, usize)>); 6] = [
            (20, vec![(10, 16)]),
            (7, vec![(0, 20)]),
            (15, vec![(11, 19)]),
            (40, vec![]),
            (40, vec![(19, 20)]),
            (40, vec![(0, 1)]),
        ];
        for ((cutoff_q, want), (_, phred)) in expected.iter().zip(reads()) {
            assert_eq!(best_segment(&phred, *cutoff_q), *want);
        }
    }

    #[test]
    fn split_matches_chopper() {
        // (cutoff, min_length, expected) — from chopper split_by_low_quality_strategy_test, window=1.
        let cases: [(u8, usize, Vec<(usize, usize)>); 6] = [
            (20, 3, vec![(6, 9), (10, 16)]),
            (7, 3, vec![(4, 15), (17, 20)]),
            (15, 3, vec![(4, 7), (14, 19)]),
            (40, 3, vec![]),
            (40, 1, vec![(19, 20)]),
            (40, 1, vec![(0, 1)]),
        ];
        for ((cutoff, min_length, want), (_, phred)) in cases.iter().zip(reads()) {
            assert_eq!(split_low_quality(&phred, *cutoff, *min_length, 1), *want);
        }
    }

    #[test]
    fn split_window_tolerates_short_dips() {
        // chopper window_test: III#IIII###III with Q40=I, Q2=#
        let phred: Vec<u8> = b"III#IIII###III".iter().map(|&b| b - 33).collect();
        assert_eq!(split_low_quality(&phred, 10, 1, 1), vec![(0, 3), (4, 8), (11, 14)]);
        assert_eq!(split_low_quality(&phred, 10, 1, 3), vec![(0, 8), (11, 14)]);
        assert_eq!(split_low_quality(&phred, 10, 1, 4), vec![(0, 14)]);
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test trim::strategies`
Expected: FAIL (undefined).

- [ ] **Step 3: Implement `src/trim/strategies.rs`**

```rust
use crate::qual::phred_to_prob;

/// Trim low-quality bases from both ends until reaching a base with phred >= cutoff.
/// Ported from chopper TrimByQualityStrategy (inputs are raw phred here).
pub fn trim_by_quality(phred: &[u8], cutoff: u8) -> Vec<(usize, usize)> {
    let len = phred.len();
    let mut start = 0;
    while start < len && phred[start] < cutoff {
        start += 1;
    }
    let mut end = len;
    while end > start && phred[end - 1] < cutoff {
        end -= 1;
    }
    if end <= start { vec![] } else { vec![(start, end)] }
}

/// Modified Mott: the single segment with the lowest cumulative error probability.
/// Ported from chopper HighestQualityTrimStrategy; `cutoff_q` is converted to a
/// probability cutoff exactly as chopper did at its call site.
pub fn best_segment(phred: &[u8], cutoff_q: u8) -> Vec<(usize, usize)> {
    let cutoff = phred_to_prob(cutoff_q);
    let mut best_start = usize::MAX;
    let mut best_end = usize::MAX;
    let mut best_cumulative_error = 0.0;
    let mut best_length = 0usize;

    let mut current_start = 0usize;
    let mut current_cumulative_error = -1.0;
    for (i, &q) in phred.iter().enumerate() {
        let prob_error = cutoff - phred_to_prob(q);
        if current_cumulative_error < 0.0 {
            current_cumulative_error = 0.0;
            current_start = i;
        }
        current_cumulative_error += prob_error;
        if best_cumulative_error < current_cumulative_error
            || (best_cumulative_error == current_cumulative_error
                && best_length < i - current_start + 1)
        {
            best_start = current_start;
            best_end = i;
            best_cumulative_error = current_cumulative_error;
            best_length = i - current_start + 1;
        }
    }
    if best_start == usize::MAX {
        vec![]
    } else {
        vec![(best_start, best_end + 1)]
    }
}

/// Split into high-quality segments separated by runs of >= `window` low-quality
/// bases; drop segments shorter than `min_length`. Ported from chopper
/// SplitByLowQualityStrategy.
pub fn split_low_quality(phred: &[u8], cutoff: u8, min_length: usize, window: usize) -> Vec<(usize, usize)> {
    let window = window.max(1);
    let mut segments = Vec::new();
    let mut segment_start: Option<usize> = None;
    let mut last_good: Option<usize> = None;
    let mut bad_run = 0usize;

    let mut push = |start: usize, end: usize, out: &mut Vec<(usize, usize)>| {
        if end - start >= min_length {
            out.push((start, end));
        }
    };

    for (i, &q) in phred.iter().enumerate() {
        if q >= cutoff {
            if segment_start.is_none() {
                segment_start = Some(i);
            }
            last_good = Some(i);
            bad_run = 0;
        } else {
            bad_run += 1;
            if bad_run >= window {
                if let (Some(s), Some(lg)) = (segment_start, last_good) {
                    push(s, lg + 1, &mut segments);
                }
                segment_start = None;
                last_good = None;
            }
        }
    }
    if let (Some(s), Some(lg)) = (segment_start, last_good) {
        push(s, lg + 1, &mut segments);
    }
    segments
}
```

Append Step 1's tests.

- [ ] **Step 4: Run to verify pass**

Run: `cargo test trim::strategies`
Expected: PASS.

- [ ] **Step 5: Write failing test for `apply` (composition) in `src/trim/mod.rs`**

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_quality_op_is_fixed_crop() {
        let phred = vec![30u8; 20];
        let plan = TrimPlan { head: 5, tail: 3, quality: None };
        assert_eq!(apply(20, &phred, &plan, 1), vec![(5, 17)]);
    }

    #[test]
    fn crop_then_quality_offsets_back() {
        // 20 bases; crop head 2; then trim_qual on the remaining window.
        // phred: low low [good...] -> after crop the good region starts at 2.
        let mut phred = vec![40u8; 20];
        phred[0] = 2;
        phred[1] = 2;
        let plan = TrimPlan { head: 2, tail: 0, quality: Some(QualityOp::TrimQual(30)) };
        assert_eq!(apply(20, &phred, &plan, 1), vec![(2, 20)]);
    }

    #[test]
    fn min_length_drops_short_segments() {
        let phred = vec![40u8; 4];
        let plan = TrimPlan { head: 0, tail: 0, quality: None };
        assert_eq!(apply(4, &phred, &plan, 5), Vec::<(usize, usize)>::new());
    }

    #[test]
    fn empty_when_crop_exceeds_length() {
        let phred = vec![40u8; 4];
        let plan = TrimPlan { head: 3, tail: 3, quality: None };
        assert_eq!(apply(4, &phred, &plan, 1), Vec::<(usize, usize)>::new());
    }
}
```

- [ ] **Step 6: Implement `src/trim/mod.rs`**

```rust
pub mod strategies;

use strategies::{best_segment, split_low_quality, trim_by_quality};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QualityOp {
    TrimQual(u8),
    BestSegment(u8),
    Split { cutoff: u8, window: usize },
}

#[derive(Debug, Clone)]
pub struct TrimPlan {
    pub head: usize,
    pub tail: usize,
    pub quality: Option<QualityOp>,
}

/// Fixed crop first (positional), then the chosen quality op on the cropped
/// window, offsetting intervals back to original coordinates. Every returned
/// segment is >= `min_length`.
pub fn apply(seq_len: usize, phred: &[u8], plan: &TrimPlan, min_length: usize) -> Vec<(usize, usize)> {
    let start = plan.head.min(seq_len);
    let end = seq_len.saturating_sub(plan.tail).max(start);
    if start >= end {
        return vec![];
    }
    let window_phred = &phred[start..end];
    let inner = match &plan.quality {
        None => vec![(0, window_phred.len())],
        Some(QualityOp::TrimQual(q)) => trim_by_quality(window_phred, *q),
        Some(QualityOp::BestSegment(q)) => best_segment(window_phred, *q),
        Some(QualityOp::Split { cutoff, window }) => {
            split_low_quality(window_phred, *cutoff, min_length, *window)
        }
    };
    inner
        .into_iter()
        .map(|(s, e)| (s + start, e + start))
        .filter(|&(s, e)| e - s >= min_length)
        .collect()
}
```

Append Step 5's tests.

- [ ] **Step 7: Run and commit**

Run: `cargo test trim::`
Expected: PASS.

```bash
git add src/trim/
git commit -m "feat: port chopper trim algorithms + composable apply()"
```

---

### Task 5: Record model + FASTQ I/O (`src/record.rs`, `src/io/fastq.rs`)

**Files:**
- Replace: `src/record.rs`
- Create: `src/io/fastq.rs`
- Modify: `src/io/mod.rs` (add `pub mod fastq;`)

**Interfaces:**
- Produces:
  - `pub struct ReadRecord { pub name: Vec<u8>, pub seq: Vec<u8>, pub qual: Vec<u8> }` where `qual` is **raw phred** and `name` is the full FASTQ header line minus the leading `@`.
  - `pub fn crate::io::fastq::reader(input: Option<&Path>, gz: bool) -> anyhow::Result<Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send>>` (the `+ Send` matters: Task 9's parallel path needs a `Send` iterator)
  - `pub fn crate::io::fastq::write_segment<W: Write>(w: &mut W, name: &[u8], seq: &[u8], phred: &[u8], total_segments: usize, segment_idx: usize) -> io::Result<()>`

- [ ] **Step 1: Write `src/record.rs`**

```rust
/// Format-neutral read carrier. `qual` holds raw Phred scores (0-based), not ASCII.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadRecord {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}
```

- [ ] **Step 2: Write failing test for the FASTQ writer**

Add to `src/io/fastq.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_single_segment_verbatim_header() {
        let mut out = Vec::new();
        write_segment(&mut out, b"read1 desc", b"ACGT", &[40, 40, 40, 40], 1, 0).unwrap();
        assert_eq!(out, b"@read1 desc\nACGT\n+\nIIII\n");
    }

    #[test]
    fn split_segment_suffixes_id_before_desc() {
        let mut out = Vec::new();
        write_segment(&mut out, b"read1 desc", b"AC", &[40, 40], 2, 1, ).unwrap();
        assert_eq!(out, b"@read1_segment_2 desc\nAC\n+\nII\n");
    }

    #[test]
    fn roundtrip_reader_writer() {
        let fq = b"@r1\nACGT\n+\nIIII\n@r2 x\nTT\n+\n!!\n";
        let recs: Vec<ReadRecord> = reader_from_slice(fq).map(|r| r.unwrap()).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, b"r1");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[0].qual, vec![40, 40, 40, 40]); // 'I' = 73 - 33
        assert_eq!(recs[1].qual, vec![0, 0]); // '!' = 33 - 33
    }
}
```

- [ ] **Step 3: Implement `src/io/fastq.rs`**

```rust
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use seq_io::fastq::{Reader, Record};

use crate::record::ReadRecord;

/// Build a streaming FASTQ record iterator over a file (or stdin when `input`
/// is `None`), transparently decompressing gzip when `gz` is true.
pub fn reader(
    input: Option<&Path>,
    gz: bool,
) -> anyhow::Result<Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send>> {
    let raw: Box<dyn Read + Send> = match input {
        Some(p) => Box::new(File::open(p)?),
        None => Box::new(io::stdin()),
    };
    let buffered = BufReader::new(raw);
    let inner: Box<dyn Read + Send> = if gz {
        Box::new(MultiGzDecoder::new(buffered))
    } else {
        Box::new(buffered)
    };
    Ok(Box::new(RecordIter { reader: Reader::new(inner) }))
}

struct RecordIter<R: Read> {
    reader: Reader<R>,
}

impl<R: Read> Iterator for RecordIter<R> {
    type Item = anyhow::Result<ReadRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        let rec = self.reader.next()?;
        Some(rec.map_err(anyhow::Error::from).map(|r| ReadRecord {
            name: r.head().to_vec(),
            seq: r.seq().to_vec(),
            qual: r.qual().iter().map(|&b| b.saturating_sub(33)).collect(),
        }))
    }
}

/// Write one output segment as a FASTQ record. On splits (`total_segments > 1`)
/// the id gets a `_segment_N` suffix inserted before any description, matching
/// chopper's convention. `phred` is raw; ASCII is emitted by adding 33.
pub fn write_segment<W: Write>(
    w: &mut W,
    name: &[u8],
    seq: &[u8],
    phred: &[u8],
    total_segments: usize,
    segment_idx: usize,
) -> io::Result<()> {
    w.write_all(b"@")?;
    if total_segments > 1 {
        let (id, desc) = split_head(name);
        w.write_all(id)?;
        write!(w, "_segment_{}", segment_idx + 1)?;
        if let Some(d) = desc {
            w.write_all(b" ")?;
            w.write_all(d)?;
        }
    } else {
        w.write_all(name)?;
    }
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    let ascii: Vec<u8> = phred.iter().map(|&q| q + 33).collect();
    w.write_all(&ascii)?;
    w.write_all(b"\n")
}

fn split_head(name: &[u8]) -> (&[u8], Option<&[u8]>) {
    match name.iter().position(|&b| b == b' ') {
        Some(i) => (&name[..i], Some(&name[i + 1..])),
        None => (name, None),
    }
}

#[cfg(test)]
fn reader_from_slice(bytes: &'static [u8]) -> RecordIter<&'static [u8]> {
    RecordIter { reader: Reader::new(bytes) }
}
```

Add `pub mod fastq;` to `src/io/mod.rs` (create the file with just that line for now; Task 6 fills the rest).

- [ ] **Step 4: Run to verify pass**

Run: `cargo test io::fastq`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/record.rs src/io/
git commit -m "feat: ReadRecord + seq_io FASTQ reader (raw-phred) and segment writer"
```

---

### Task 6: Format detection (`src/io/mod.rs`)

**Files:**
- Replace: `src/io/mod.rs` (keep `pub mod fastq;`)

**Interfaces:**
- Produces:
  - `pub enum Format { Fastq, FastqGz, Bam }` (derive `Debug, Clone, Copy, PartialEq, Eq`)
  - `pub fn from_extension(path: &Path) -> Option<Format>`
  - `pub fn detect_input(path: Option<&Path>, sniff: &[u8]) -> anyhow::Result<Format>`
  - `pub fn resolve_output(path: Option<&Path>, input: Format) -> Format`

- [ ] **Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn extensions() {
        assert_eq!(from_extension(Path::new("x.fastq")), Some(Format::Fastq));
        assert_eq!(from_extension(Path::new("x.fq")), Some(Format::Fastq));
        assert_eq!(from_extension(Path::new("x.fastq.gz")), Some(Format::FastqGz));
        assert_eq!(from_extension(Path::new("x.fq.gz")), Some(Format::FastqGz));
        assert_eq!(from_extension(Path::new("x.bam")), Some(Format::Bam));
        assert_eq!(from_extension(Path::new("x.txt")), None);
    }

    #[test]
    fn stdin_sniff_falls_back_to_magic() {
        // no path -> sniff. gzip magic 1f 8b -> FastqGz; '@' -> Fastq; BAM magic -> Bam.
        assert_eq!(detect_input(None, &[0x1f, 0x8b, 0x08]).unwrap(), Format::FastqGz);
        assert_eq!(detect_input(None, b"@read").unwrap(), Format::Fastq);
        assert_eq!(detect_input(None, b"BAM\x01").unwrap(), Format::Bam);
    }

    #[test]
    fn output_mirrors_input_when_no_path() {
        assert_eq!(resolve_output(None, Format::Bam), Format::Bam);
        assert_eq!(resolve_output(None, Format::Fastq), Format::Fastq);
        assert_eq!(resolve_output(Some(Path::new("o.bam")), Format::Fastq), Format::Bam);
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test io::tests`
Expected: FAIL.

- [ ] **Step 3: Implement (prepend to `src/io/mod.rs`, keep `pub mod fastq;`)**

```rust
pub mod fastq;

use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fastq,
    FastqGz,
    Bam,
}

/// Extension-based detection. Recognises `.fastq`/`.fq`, the `.gz` variants, and `.bam`.
pub fn from_extension(path: &Path) -> Option<Format> {
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    if name.ends_with(".fastq.gz") || name.ends_with(".fq.gz") {
        Some(Format::FastqGz)
    } else if name.ends_with(".fastq") || name.ends_with(".fq") {
        Some(Format::Fastq)
    } else if name.ends_with(".bam") {
        Some(Format::Bam)
    } else {
        None
    }
}

/// Detect input format from the path extension; if unknown or reading stdin,
/// fall back to sniffing the first bytes.
pub fn detect_input(path: Option<&Path>, sniff: &[u8]) -> anyhow::Result<Format> {
    if let Some(p) = path {
        if let Some(f) = from_extension(p) {
            return Ok(f);
        }
    }
    if sniff.starts_with(&[0x1f, 0x8b]) {
        Ok(Format::FastqGz)
    } else if sniff.starts_with(b"BAM\x01") {
        Ok(Format::Bam)
    } else if sniff.first() == Some(&b'@') {
        Ok(Format::Fastq)
    } else {
        anyhow::bail!("cannot determine input format; pass --in-format")
    }
}

/// Output format from the path extension, else mirror the input format.
pub fn resolve_output(path: Option<&Path>, input: Format) -> Format {
    path.and_then(from_extension).unwrap_or(input)
}
```

- [ ] **Step 4: Run and commit**

Run: `cargo test io::tests`
Expected: PASS.

```bash
git add src/io/mod.rs
git commit -m "feat: format detection (extension + magic sniff + stdout mirror)"
```

---

### Task 7: Config + single-thread FASTQ pipeline (`src/config.rs`, `src/pipeline.rs`)

**Files:**
- Replace: `src/config.rs`, `src/pipeline.rs`
- Modify: `src/lib.rs` (`run` dispatches to the FASTQ pipeline)

**Interfaces:**
- Consumes: `filter::{FilterConfig, passes}`, `trim::{TrimPlan, apply}`, `io::fastq::{reader, write_segment}`, `io::{Format, resolve_output}`.
- Produces:
  - `pub struct IoConfig { pub input: Option<PathBuf>, pub output: Option<PathBuf>, pub in_format: Option<Format>, pub out_format: Option<Format> }`
  - `pub struct Config { pub io: IoConfig, pub filter: FilterConfig, pub trim: TrimPlan, pub threads: usize }`
  - `pub struct Stats { pub input_reads: u64, pub output_reads: u64 }`
  - `pub fn pipeline::run_fastq_seq(recs, writer, cfg) -> anyhow::Result<Stats>`

- [ ] **Step 1: Write `src/config.rs`**

```rust
use std::path::PathBuf;

use crate::filter::FilterConfig;
use crate::io::Format;
use crate::trim::TrimPlan;

#[derive(Debug, Clone)]
pub struct IoConfig {
    pub input: Option<PathBuf>,
    pub output: Option<PathBuf>,
    pub in_format: Option<Format>,
    pub out_format: Option<Format>,
}

#[derive(Debug, Clone)]
pub struct Config {
    pub io: IoConfig,
    pub filter: FilterConfig,
    pub trim: TrimPlan,
    pub threads: usize,
}
```

- [ ] **Step 2: Write failing test for `run_fastq_seq`**

Add to `src/pipeline.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::FilterConfig;
    use crate::qual::QualMode;
    use crate::record::ReadRecord;
    use crate::trim::{QualityOp, TrimPlan};

    fn rec(name: &str, seq: &[u8], phred: Vec<u8>) -> ReadRecord {
        ReadRecord { name: name.as_bytes().to_vec(), seq: seq.to_vec(), qual: phred }
    }

    fn base_filter() -> FilterConfig {
        FilterConfig {
            min_length: 1, max_length: usize::MAX, min_qual: 0.0, max_qual: 1000.0,
            min_gc: None, max_gc: None, qual_mode: QualMode::Mean,
        }
    }

    #[test]
    fn fixed_crop_writes_one_segment() {
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 1, tail: 1, quality: None },
            threads: 1,
        };
        let recs = vec![Ok(rec("r1", b"ACGT", vec![40, 40, 40, 40]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@r1\nCG\n+\nII\n");
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
    }

    #[test]
    fn split_writes_suffixed_segments() {
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 0, tail: 0, quality: Some(QualityOp::Split { cutoff: 10, window: 1 }) },
            threads: 1,
        };
        // good(3) bad(1) good(3): I I I # I I I  -> two segments (0,3),(4,7)
        let phred: Vec<u8> = b"III#III".iter().map(|&b| b - 33).collect();
        let recs = vec![Ok(rec("r1", b"AAATAAA", phred))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@r1_segment_1\nAAA\n+\nIII\n@r1_segment_2\nAAA\n+\nIII\n");
        assert_eq!((stats.input_reads, stats.output_reads), (1, 2));
    }

    #[test]
    fn filtered_read_produces_no_output() {
        let mut f = base_filter();
        f.min_length = 10;
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: f,
            trim: TrimPlan { head: 0, tail: 0, quality: None },
            threads: 1,
        };
        let recs = vec![Ok(rec("short", b"ACGT", vec![40; 4]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert!(out.is_empty());
        assert_eq!((stats.input_reads, stats.output_reads), (1, 0));
    }
}
```

- [ ] **Step 3: Implement `src/pipeline.rs`**

```rust
use std::io::Write;

use crate::config::Config;
use crate::io::fastq::write_segment;
use crate::record::ReadRecord;
use crate::{filter, trim};

#[derive(Debug, Default, Clone, Copy)]
pub struct Stats {
    pub input_reads: u64,
    pub output_reads: u64,
}

/// Single-threaded FASTQ pipeline: filter -> trim -> write each surviving segment.
pub fn run_fastq_seq<W: Write>(
    records: impl Iterator<Item = anyhow::Result<ReadRecord>>,
    writer: &mut W,
    cfg: &Config,
) -> anyhow::Result<Stats> {
    let mut stats = Stats::default();
    for rec in records {
        let rec = rec?;
        stats.input_reads += 1;
        if !filter::passes(&rec.seq, &rec.qual, &cfg.filter) {
            continue;
        }
        let intervals = trim::apply(rec.seq.len(), &rec.qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            write_segment(writer, &rec.name, &rec.seq[s..e], &rec.qual[s..e], total, idx)?;
            stats.output_reads += 1;
        }
    }
    Ok(stats)
}
```

- [ ] **Step 4: Wire `run()` in `src/lib.rs`**

Replace the `run` stub:

```rust
use std::io::{BufWriter, Write};
use std::path::Path;

/// Resolve formats, open reader/writer, dispatch to the matching pipeline.
pub fn run(cfg: Config) -> anyhow::Result<()> {
    use io::Format;

    let in_path = cfg.io.input.as_deref();
    // Sniff a few bytes for stdin/unknown-extension detection.
    let in_fmt = match cfg.io.in_format {
        Some(f) => f,
        None => {
            let sniff = peek_input(in_path)?;
            io::detect_input(in_path, &sniff)?
        }
    };
    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), in_fmt));

    let mut writer: BufWriter<Box<dyn Write>> = BufWriter::new(match cfg.io.output.as_deref() {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    });

    match (in_fmt, out_fmt) {
        (Format::Fastq, Format::Fastq)
        | (Format::FastqGz, Format::Fastq)
        | (Format::Fastq, Format::FastqGz)
        | (Format::FastqGz, Format::FastqGz) => {
            let gz_in = matches!(in_fmt, Format::FastqGz);
            let records = io::fastq::reader(in_path, gz_in)?;
            // gz output wrapping is added in a later task if out_fmt is FastqGz;
            // Plan 1 supports plain FASTQ output here, gz output in Task 9 note.
            let stats = pipeline::run_fastq(records, &mut writer, &cfg)?;
            writer.flush()?;
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            Ok(())
        }
        (Format::Bam, _) | (_, Format::Bam) => {
            anyhow::bail!("BAM support arrives in Plan 2")
        }
    }
}

fn peek_input(path: Option<&Path>) -> anyhow::Result<Vec<u8>> {
    use std::io::Read;
    let mut buf = vec![0u8; 4];
    if let Some(p) = path {
        let mut f = std::fs::File::open(p)?;
        let n = f.read(&mut buf)?;
        buf.truncate(n);
    }
    Ok(buf)
}
```

Note: `pipeline::run_fastq` is the parallel-capable entry added in Task 9; for now add a thin alias so `run` compiles:

```rust
// in src/pipeline.rs, temporary until Task 9:
pub fn run_fastq<W: Write>(
    records: Box<dyn Iterator<Item = anyhow::Result<ReadRecord>>>,
    writer: &mut W,
    cfg: &Config,
) -> anyhow::Result<Stats> {
    run_fastq_seq(records, writer, cfg)
}
```

- [ ] **Step 5: Run and commit**

Run: `cargo test pipeline::`
Expected: PASS (3 tests).

```bash
git add src/config.rs src/pipeline.rs src/lib.rs
git commit -m "feat: Config + single-thread FASTQ pipeline wired into run()"
```

---

### Task 8: CLI (`src/cli.rs`)

**Files:**
- Replace: `src/cli.rs`

**Interfaces:**
- Consumes: everything above.
- Produces: `pub fn parse() -> anyhow::Result<Config>` (parses argv, validates, builds `Config`).

- [ ] **Step 1: Write failing CLI integration tests in `tests/cli.rs`**

```rust
use assert_cmd::Command;
use predicates::prelude::*;

fn chopping() -> Command {
    Command::cargo_bin("chopping").unwrap()
}

#[test]
fn head_tail_crop_over_stdin() {
    chopping()
        .args(["--head-crop", "1", "--tail-crop", "1", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .success()
        .stdout("@r1\nCG\n+\nII\n");
}

#[test]
fn mutually_exclusive_quality_ops_error() {
    chopping()
        .args(["--trim-qual", "10", "--best-segment", "10", "--in-format", "fastq"])
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .failure()
        .stderr(predicate::str::contains("mutually exclusive"));
}

#[test]
fn min_length_filters() {
    chopping()
        .args(["--min-length", "10", "--in-format", "fastq"])
        .write_stdin("@short\nACGT\n+\nIIII\n")
        .assert()
        .success()
        .stdout(""); // filtered out
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test --test cli`
Expected: FAIL (`parse()` returns empty Config; flags unknown).

- [ ] **Step 3: Implement `src/cli.rs`**

```rust
use std::path::PathBuf;

use clap::Parser;

use crate::config::{Config, IoConfig};
use crate::filter::FilterConfig;
use crate::io::Format;
use crate::qual::QualMode;
use crate::trim::{QualityOp, TrimPlan};

#[derive(Parser, Debug)]
#[command(author, version, about = "Tag-aware long-read trimmer", long_about = None)]
struct Cli {
    #[arg(short = 'i', long, help_heading = "Setup")]
    input: Option<PathBuf>,
    #[arg(short = 'o', long, help_heading = "Setup")]
    output: Option<PathBuf>,
    #[arg(long, value_enum, help_heading = "Setup")]
    in_format: Option<FormatArg>,
    #[arg(long, value_enum, help_heading = "Setup")]
    out_format: Option<FormatArg>,
    #[arg(short = 't', long, default_value_t = 4, help_heading = "Setup")]
    threads: usize,

    #[arg(short = 'l', long, default_value_t = 1, help_heading = "Filtering")]
    min_length: usize,
    #[arg(short = 'L', long, help_heading = "Filtering")]
    max_length: Option<usize>,
    #[arg(short = 'q', long, default_value_t = 0.0, help_heading = "Filtering")]
    min_qual: f64,
    #[arg(short = 'Q', long, default_value_t = 1000.0, help_heading = "Filtering")]
    max_qual: f64,
    #[arg(short = 'g', long, help_heading = "Filtering")]
    min_gc: Option<f64>,
    #[arg(short = 'G', long, help_heading = "Filtering")]
    max_gc: Option<f64>,
    #[arg(short = 'm', long, value_enum, default_value_t = QualModeArg::Mean, help_heading = "Filtering")]
    qual_mode: QualModeArg,

    #[arg(short = 'H', long, default_value_t = 0, help_heading = "Trimming")]
    head_crop: usize,
    #[arg(short = 'T', long, default_value_t = 0, help_heading = "Trimming")]
    tail_crop: usize,
    #[arg(long, help_heading = "Trimming")]
    trim_qual: Option<u8>,
    #[arg(long, help_heading = "Trimming")]
    best_segment: Option<u8>,
    #[arg(long, help_heading = "Trimming")]
    split_qual: Option<u8>,
    #[arg(long, default_value_t = 1, help_heading = "Trimming")]
    split_window: usize,
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
enum FormatArg { Fastq, FastqGz, Bam }

impl From<FormatArg> for Format {
    fn from(f: FormatArg) -> Self {
        match f {
            FormatArg::Fastq => Format::Fastq,
            FormatArg::FastqGz => Format::FastqGz,
            FormatArg::Bam => Format::Bam,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
enum QualModeArg { Mean, Arithmetic, Median }

impl From<QualModeArg> for QualMode {
    fn from(m: QualModeArg) -> Self {
        match m {
            QualModeArg::Mean => QualMode::Mean,
            QualModeArg::Arithmetic => QualMode::Arithmetic,
            QualModeArg::Median => QualMode::Median,
        }
    }
}

pub fn parse() -> anyhow::Result<Config> {
    let c = Cli::parse();

    // Mutual exclusion of the three quality trim ops.
    let n_quality = [c.trim_qual.is_some(), c.best_segment.is_some(), c.split_qual.is_some()]
        .iter()
        .filter(|&&b| b)
        .count();
    if n_quality > 1 {
        anyhow::bail!("--trim-qual, --best-segment and --split-qual are mutually exclusive");
    }
    let quality = if let Some(q) = c.trim_qual {
        Some(QualityOp::TrimQual(q))
    } else if let Some(q) = c.best_segment {
        Some(QualityOp::BestSegment(q))
    } else if let Some(q) = c.split_qual {
        Some(QualityOp::Split { cutoff: q, window: c.split_window })
    } else {
        None
    };

    Ok(Config {
        io: IoConfig {
            input: c.input,
            output: c.output,
            in_format: c.in_format.map(Into::into),
            out_format: c.out_format.map(Into::into),
        },
        filter: FilterConfig {
            min_length: c.min_length,
            max_length: c.max_length.unwrap_or(usize::MAX),
            min_qual: c.min_qual,
            max_qual: c.max_qual,
            min_gc: c.min_gc,
            max_gc: c.max_gc,
            qual_mode: c.qual_mode.into(),
        },
        trim: TrimPlan { head: c.head_crop, tail: c.tail_crop, quality },
        threads: c.threads.max(1),
    })
}
```

- [ ] **Step 4: Run to verify pass**

Run: `cargo test --test cli`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/cli.rs tests/cli.rs
git commit -m "feat: clap CLI (kebab flags, min/max shorts, qual-mode, mutual-exclusion)"
```

---

### Task 9: Parallel pipeline + gz output

**Files:**
- Modify: `src/pipeline.rs` (replace the temporary `run_fastq` alias with a threads-aware version)
- Modify: `src/lib.rs` (wrap writer in gz when `out_fmt == FastqGz`)

**Interfaces:**
- Produces: `pub fn pipeline::run_fastq<W: Write + Send>(records, writer, cfg) -> anyhow::Result<Stats>` — sequential when `cfg.threads == 1`, otherwise a rayon work pool + dedicated writer thread (ported from chopper `main.rs`). Output order is unordered for `threads > 1`.

- [ ] **Step 1: Write failing test — parallel output equals sequential as a set**

Add to `src/pipeline.rs` tests:

```rust
#[test]
fn parallel_matches_sequential_as_multiset() {
    use crate::config::IoConfig;
    let mk = |threads| Config {
        io: IoConfig { input: None, output: None, in_format: None, out_format: None },
        filter: base_filter(),
        trim: TrimPlan { head: 0, tail: 0, quality: Some(QualityOp::TrimQual(20)) },
        threads,
    };
    // Owned records (ReadRecord: Clone); wrap in Ok at iteration time so each run
    // gets a fresh Send iterator. anyhow::Error is not Clone, so we can't clone a
    // Vec<Result<..>> — clone the Vec<ReadRecord> and re-wrap instead.
    let recs: Vec<ReadRecord> = (0..500)
        .map(|i| rec(&format!("r{i}"), b"ACGTACGTAC", vec![40; 10]))
        .collect();

    let mut seq_out = Vec::new();
    run_fastq(recs.clone().into_iter().map(anyhow::Ok), &mut seq_out, &mk(1)).unwrap();

    let mut par_out = Vec::new();
    run_fastq(recs.into_iter().map(anyhow::Ok), &mut par_out, &mk(4)).unwrap();

    let sort_records = |bytes: &[u8]| {
        let mut v: Vec<Vec<u8>> = bytes
            .split(|&b| b == b'@')
            .filter(|s| !s.is_empty())
            .map(|s| s.to_vec())
            .collect();
        v.sort();
        v
    };
    assert_eq!(sort_records(&seq_out), sort_records(&par_out));
}
```

`anyhow::Ok` fixes the error type to `anyhow::Error` so `run_fastq`'s `Item = anyhow::Result<ReadRecord>` bound is satisfied. `run_fastq` accepts `I: Iterator<Item = anyhow::Result<ReadRecord>> + Send`.

- [ ] **Step 2: Run to verify failure**

Run: `cargo test pipeline::parallel`
Expected: FAIL (current `run_fastq` is the single-thread alias; signature/behavior differs).

- [ ] **Step 3: Implement threads-aware `run_fastq`**

Replace the temporary alias in `src/pipeline.rs`:

```rust
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};

/// Format each surviving segment of one record into an owned FASTQ byte buffer.
fn render_record(rec: &ReadRecord, cfg: &Config) -> (u64, Vec<u8>) {
    if !filter::passes(&rec.seq, &rec.qual, &cfg.filter) {
        return (0, Vec::new());
    }
    let intervals = trim::apply(rec.seq.len(), &rec.qual, &cfg.trim, cfg.filter.min_length);
    let total = intervals.len();
    let mut buf = Vec::new();
    let mut out = 0u64;
    for (idx, (s, e)) in intervals.into_iter().enumerate() {
        write_segment(&mut buf, &rec.name, &rec.seq[s..e], &rec.qual[s..e], total, idx).unwrap();
        out += 1;
    }
    (out, buf)
}

pub fn run_fastq<W, I>(records: I, writer: &mut W, cfg: &Config) -> anyhow::Result<Stats>
where
    W: Write + Send,
    I: Iterator<Item = anyhow::Result<ReadRecord>> + Send,
{
    if cfg.threads <= 1 {
        return run_fastq_seq(records, writer, cfg);
    }

    let pool = rayon::ThreadPoolBuilder::new().num_threads(cfg.threads).build()?;
    let input_reads = AtomicU64::new(0);
    let output_reads = AtomicU64::new(0);
    let (tx, rx) = crossbeam_channel::bounded::<Vec<u8>>(cfg.threads * 4);

    let write_err: std::sync::Mutex<Option<std::io::Error>> = std::sync::Mutex::new(None);

    pool.in_place_scope(|scope| {
        // Writer task drains rendered buffers in arrival order.
        scope.spawn(|_| {
            for buf in rx.iter() {
                if let Err(e) = writer.write_all(&buf) {
                    *write_err.lock().unwrap() = Some(e);
                    break;
                }
            }
        });

        records.par_bridge().for_each(|rec| {
            let rec = match rec {
                Ok(r) => r,
                Err(_) => return,
            };
            input_reads.fetch_add(1, Ordering::Relaxed);
            let (out, buf) = render_record(&rec, cfg);
            if out > 0 {
                output_reads.fetch_add(out, Ordering::Relaxed);
                let _ = tx.send(buf);
            }
        });
        drop(tx);
    });

    if let Some(e) = write_err.lock().unwrap().take() {
        return Err(e.into());
    }
    Ok(Stats {
        input_reads: input_reads.load(Ordering::Relaxed),
        output_reads: output_reads.load(Ordering::Relaxed),
    })
}
```

Note on record-parse errors: the sequential path surfaces them via `?`; the parallel path currently drops them. Track a shared error slot the same way as `write_err` if strict parse-error propagation is required — for Plan 1, malformed FASTQ over multiple threads aborts the affected read only. Document this in the code comment.

- [ ] **Step 4: Add gz output wrapping in `src/lib.rs`**

In `run`, when `out_fmt == Format::FastqGz`, wrap the writer:

```rust
let base_writer: Box<dyn Write + Send> = match cfg.io.output.as_deref() {
    Some(p) => Box::new(std::fs::File::create(p)?),
    None => Box::new(std::io::stdout()),
};
let writer_inner: Box<dyn Write + Send> = if matches!(out_fmt, Format::FastqGz) {
    Box::new(flate2::write::GzEncoder::new(base_writer, flate2::Compression::default()))
} else {
    base_writer
};
let mut writer = BufWriter::new(writer_inner);
```

Replace the earlier writer construction accordingly, and ensure `writer.flush()?` plus (for gz) the encoder is finished by dropping `writer` before returning.

- [ ] **Step 5: Run the full suite and clippy**

Run: `cargo test && cargo clippy -- -D warnings`
Expected: all tests PASS; clippy clean.

- [ ] **Step 6: Add an end-to-end gz round-trip test in `tests/cli.rs`**

```rust
#[test]
fn gz_output_roundtrips() {
    use std::io::Read;
    let dir = tempfile::tempdir().unwrap(); // add tempfile to dev-deps
    let out = dir.path().join("out.fastq.gz");
    chopping()
        .args(["--in-format", "fastq", "-o"]).arg(&out)
        .write_stdin("@r1\nACGT\n+\nIIII\n")
        .assert()
        .success();
    let mut gz = flate2::read::MultiGzDecoder::new(std::fs::File::open(&out).unwrap());
    let mut s = String::new();
    gz.read_to_string(&mut s).unwrap();
    assert_eq!(s, "@r1\nACGT\n+\nIIII\n");
}
```

Add `tempfile = "3"` and `flate2 = { version = "1", features = ["zlib-ng"], default-features = false }` to `[dev-dependencies]` (flate2 is already a normal dep; the dev entry is only if the test needs the decoder directly — it can reuse the main dep, so this line is optional).

- [ ] **Step 7: Commit**

```bash
git add src/pipeline.rs src/lib.rs tests/cli.rs Cargo.toml Cargo.lock
git commit -m "feat: parallel FASTQ pipeline (rayon + writer thread) and gz output"
```

---

## Self-Review

**Spec coverage (FASTQ scope of the v1 spec):**
- Length/mean-quality/GC filtering → Task 3. ✅
- `--qual-mode` {mean|arithmetic|median} on `-q/-Q` only → Tasks 2, 3, 8. ✅
- Four trim ops as their own flags, no `--trim-approach`/`--cutoff` → Tasks 4, 8. ✅
- Fixed crop composes; three quality ops mutually exclusive → Tasks 4 (`apply`), 8 (validation). ✅
- `_segment_N` naming on splits → Task 5. ✅
- min/max short-flag scheme, `-H/-T`, removed `--contam`/`--inverse` → Task 8. ✅
- Format auto-detection (extension + magic + stdout mirror); `--in-format`/`--out-format` optional → Tasks 6, 7. ✅
- FASTQ(.gz) I/O via seq_io + flate2 → Tasks 5, 9. ✅
- Parallelism (parallelize per-read work, not parsing; `--threads 1` deterministic) → Task 9. ✅
- Goldens vs chopper for ported strategies → Task 4 uses chopper's exact expected intervals as unit goldens. ✅
- BAM path, MM/ML, aligned-read refusal, htslib oracle → **deferred to Plan 2** (out of scope here). Noted.

**Placeholder scan:** No "TBD"/"implement later" in deliverable code. Task 1 Step 4 intentionally creates temporary module stubs that later tasks replace — each such file is named and replaced by a specific task. The Task 7 `run_fastq` alias is explicitly temporary and replaced in Task 9. No vague "add error handling."

**Type consistency:** `ReadRecord{name,seq,qual}`, `FilterConfig`, `TrimPlan{head,tail,quality}`, `QualityOp`, `Config{io,filter,trim,threads}`, `Stats{input_reads,output_reads}`, `Format`, `QualMode` names are used identically across Tasks 2–9. `apply(seq_len, phred, plan, min_length)` signature matches every call site. `write_segment(w,name,seq,phred,total,idx)` matches pipeline + tests.

Two fixes applied during review: (1) `io::fastq::reader` returns `Box<dyn Iterator<..> + Send>` so Task 9's `par_bridge` accepts it; (2) the Task 9 multiset test wraps cloned `ReadRecord`s with `anyhow::Ok` (since `anyhow::Error` is not `Clone`, a `Vec<Result<..>>` can't be cloned). Both are reflected in the task code above.

---

## Plan 2 preview (separate document, written next)

**uBAM + MM/ML reconstruction** — builds directly on Plan 1's `Config`, `filter`, `trim`, `pipeline` seams:
1. `io/bam.rs`: noodles reader → `ReadRecord` + a `BamTags` sidecar (raw MM string, raw ML bytes, MN, passthrough aux, flags); **aligned-read refusal** (`flags().is_unmapped()` false → hard error). BAM writer via `RecordBuf`.
2. `mods/parse.rs`, `mods/reconstruct.rs`, `mods/serialize.rs`: the MM/ML/MN codec operating on intervals, bypassing noodles' typed parser (raw `data().get(b"MM")` / `Array::UInt8` per the bench crates).
3. Extend the record model so reconstruction rebuilds MM/ML/MN per interval; extend `pipeline` with a `run_bam` variant reusing `filter`/`trim`/`apply`.
4. `tests/bam_mods_oracle.rs`: `rust-htslib` dev-dependency; **decode-equivalence** oracle (decode our output with `basemods_iter()`, assert per-position set equals original filtered to the interval, offset by `start`) on a small HG002 uBAM fixture.
```
