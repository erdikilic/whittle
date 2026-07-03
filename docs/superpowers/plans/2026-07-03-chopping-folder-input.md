# chopping — folder-merge input — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `-i/--input` accept a directory; when given one, merge all read files in that flat folder (all FASTQ-family, or all BAM) into a single trimmed output.

**Architecture:** A new `src/io/dir.rs` classifies a directory into a format family + a sorted file list and builds a single chained record iterator over the files. `run()` gets one early branch: if the input path is a directory, route to a `run_folder` helper that reuses the existing `run_fastq`/`run_bam` pipelines and writer construction. Nothing on the single-file/stdin path changes behaviorally.

**Tech Stack:** Rust 2024; existing `io::fastq`/`io::bam` readers; `noodles-{bam,sam}`; `flate2`; dev: `assert_cmd`, `tempfile`.

Design source: `docs/superpowers/specs/2026-07-03-chopping-folder-input-design.md`. Current entry point: `src/lib.rs::run` (single-file/stdin today).

## Global Constraints

- Rust 2024; `cargo clippy --all-targets -- -D warnings` clean; pristine tests.
- Folder mode triggers only when `cfg.io.input` is `Some(path)` and `path.is_dir()`. Single file and stdin are unchanged.
- **Single flat folder, non-recursive** — immediate children only.
- Read files by extension: FASTQ-family (`.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`) or `.bam` (reuse `io::from_extension`). Non-read files/subdirs ignored.
- Files processed in **sorted filename order** (deterministic).
- Folder must be homogeneous: all FASTQ-family (plain+gz may coexist) → one FASTQ output; all `.bam` → one BAM output. Mixed fastq+bam → hard error; no read files → hard error.
- BAM merge uses the **first sorted file's header** (+ `@PG` via existing `provenance_header`); records from all files stream under it; `writer.try_finish()?` runs once.
- Internal quality stays raw Phred; the merge changes only where records come from, not how they're trimmed.

## File Structure

```
src/io/dir.rs     NEW: Family, classify(dir), fastq_records(paths), bam_reader(paths)
src/io/mod.rs     + `pub mod dir;`
src/lib.rs        run(): early is_dir() branch -> run_folder; extract fastq_writer() helper (DRY)
tests/folder.rs   NEW: end-to-end folder-merge tests (FASTQ + uBAM)
README.md         document `-i <dir>` folder-merge
```

---

### Task 1: `src/io/dir.rs` — classify + merged readers

**Files:**
- Create: `src/io/dir.rs`
- Modify: `src/io/mod.rs` (add `pub mod dir;`)

**Interfaces:**
- Consumes: `crate::io::{Format, from_extension}`; `crate::io::fastq::reader(Option<&Path>, bool)`; `crate::io::bam::reader(Option<&Path>)`; `crate::record::ReadRecord`; `noodles_sam::alignment::RecordBuf`; `noodles_sam as sam`.
- Produces:
  - `pub enum Family { Fastq, Bam }` (derive Debug, Clone, Copy, PartialEq, Eq)
  - `pub fn classify(dir: &Path) -> anyhow::Result<(Family, Vec<PathBuf>)>`
  - `pub fn fastq_records(paths: &[PathBuf]) -> Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send>`
  - `pub fn bam_reader(paths: &[PathBuf]) -> anyhow::Result<(sam::Header, Box<dyn Iterator<Item = anyhow::Result<RecordBuf>>>)>`

- [ ] **Step 1: Write failing tests for `classify`**

Add to `src/io/dir.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    fn touch(dir: &std::path::Path, name: &str) {
        std::fs::write(dir.join(name), b"x").unwrap();
    }

    #[test]
    fn fastq_only_folder_sorted() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "b_1.fastq.gz");
        touch(d.path(), "a_0.fastq");
        touch(d.path(), "sequencing_summary.txt"); // ignored
        std::fs::create_dir(d.path().join("subdir")).unwrap(); // ignored
        let (fam, paths) = classify(d.path()).unwrap();
        assert_eq!(fam, Family::Fastq);
        let names: Vec<_> = paths.iter().map(|p| p.file_name().unwrap().to_str().unwrap()).collect();
        assert_eq!(names, vec!["a_0.fastq", "b_1.fastq.gz"]); // sorted, .txt/subdir excluded
    }

    #[test]
    fn bam_only_folder() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "x.bam");
        touch(d.path(), "y.bam");
        let (fam, paths) = classify(d.path()).unwrap();
        assert_eq!(fam, Family::Bam);
        assert_eq!(paths.len(), 2);
    }

    #[test]
    fn mixed_formats_error() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "a.fastq");
        touch(d.path(), "b.bam");
        let err = classify(d.path()).unwrap_err().to_string();
        assert!(err.contains("mixes"));
    }

    #[test]
    fn empty_folder_error() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "notes.txt"); // no read files
        let err = classify(d.path()).unwrap_err().to_string();
        assert!(err.contains("no FASTQ or BAM"));
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test io::dir`
Expected: FAIL (module/functions undefined).

- [ ] **Step 3: Implement `src/io/dir.rs`**

```rust
use std::path::{Path, PathBuf};

use anyhow::Context;
use noodles_sam::{self as sam, alignment::RecordBuf};

use crate::io::{Format, from_extension};
use crate::record::ReadRecord;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    Fastq,
    Bam,
}

/// Classify a directory's immediate children into a single read-file family and
/// a sorted path list. Non-read files and subdirectories are ignored. Errors if
/// the folder mixes FASTQ and BAM, or contains no read files.
pub fn classify(dir: &Path) -> anyhow::Result<(Family, Vec<PathBuf>)> {
    let mut fastq = Vec::new();
    let mut bam = Vec::new();

    let entries = std::fs::read_dir(dir)
        .with_context(|| format!("reading directory {}", dir.display()))?;
    for entry in entries {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        match from_extension(&path) {
            Some(Format::Fastq | Format::FastqGz) => fastq.push(path),
            Some(Format::Bam) => bam.push(path),
            None => {} // ignore non-read files
        }
    }

    let (family, mut paths) = match (fastq.is_empty(), bam.is_empty()) {
        (false, true) => (Family::Fastq, fastq),
        (true, false) => (Family::Bam, bam),
        (true, true) => anyhow::bail!(
            "no FASTQ or BAM read files found in directory {}",
            dir.display()
        ),
        (false, false) => anyhow::bail!(
            "directory {} mixes FASTQ ({}) and BAM ({}) files; a folder must be \
             one format to merge into a single output",
            dir.display(),
            fastq.len(),
            bam.len()
        ),
    };
    paths.sort();
    Ok((family, paths))
}

/// One chained record stream over all FASTQ-family files, opened lazily as the
/// chain reaches each file (avoids exhausting file descriptors on large folders).
/// A file-open error surfaces as an `Err` item rather than aborting construction.
pub fn fastq_records(
    paths: &[PathBuf],
) -> Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send> {
    let paths = paths.to_vec();
    Box::new(paths.into_iter().flat_map(
        |p| -> Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send> {
            let gz = matches!(from_extension(&p), Some(Format::FastqGz));
            match crate::io::fastq::reader(Some(&p), gz) {
                Ok(reader) => reader,
                Err(e) => Box::new(std::iter::once(Err(e))),
            }
        },
    ))
}

/// The first file's header plus one chained record stream over all BAM files
/// (each file's own header is read and discarded; records stream under the
/// first header — `samtools cat` semantics for homogeneous uBAM).
pub fn bam_reader(
    paths: &[PathBuf],
) -> anyhow::Result<(sam::Header, Box<dyn Iterator<Item = anyhow::Result<RecordBuf>>>)> {
    let (header, _first_records) = crate::io::bam::reader(Some(&paths[0]))?;
    let paths = paths.to_vec();
    let records = paths.into_iter().flat_map(
        |p| -> Box<dyn Iterator<Item = anyhow::Result<RecordBuf>>> {
            match crate::io::bam::reader(Some(&p)) {
                Ok((_hdr, recs)) => recs,
                Err(e) => Box::new(std::iter::once(Err(e))),
            }
        },
    );
    Ok((header, Box::new(records)))
}
```

Add `pub mod dir;` to `src/io/mod.rs` (alongside the existing `pub mod fastq;` / `pub mod bam;`).

- [ ] **Step 4: Run to verify pass**

Run: `cargo test io::dir`
Expected: PASS (4 tests). Then `cargo clippy --all-targets -- -D warnings` clean.

- [ ] **Step 5: Commit**

```bash
git add src/io/dir.rs src/io/mod.rs
git commit -m "feat: io::dir — classify a folder + merged FASTQ/BAM record streams"
```

---

### Task 2: wire `run()` folder branch + integration tests + README

**Files:**
- Modify: `src/lib.rs` (extract `fastq_writer` helper; add `run_folder`; early `is_dir()` branch in `run`)
- Create: `tests/folder.rs`
- Modify: `README.md`

**Interfaces:**
- Consumes: `io::dir::{Family, classify, fastq_records, bam_reader}`; existing `pipeline::{run_fastq, run_bam}`, `io::bam::writer`, `provenance_header`, `io::resolve_output`.
- Produces: folder-merge behavior via `run(cfg)` when the input is a directory.

- [ ] **Step 1: Extract the FASTQ writer helper (behavior-preserving refactor)**

In `src/lib.rs`, add this helper and replace the inline writer construction in `run()` (currently the `base_writer`/`writer_inner`/`BufWriter::new` block for the FASTQ path) with a call to it:

```rust
/// Build the FASTQ output writer: a file or stdout, wrapped in a gzip encoder
/// when the output format is `FastqGz`, then buffered.
fn fastq_writer(
    cfg: &Config,
    out_fmt: io::Format,
) -> anyhow::Result<BufWriter<Box<dyn Write + Send>>> {
    let base_writer: Box<dyn Write + Send> = match cfg.io.output.as_deref() {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    };
    let writer_inner: Box<dyn Write + Send> = if matches!(out_fmt, io::Format::FastqGz) {
        Box::new(flate2::write::GzEncoder::new(base_writer, flate2::Compression::default()))
    } else {
        base_writer
    };
    Ok(BufWriter::new(writer_inner))
}
```

In the existing single-file FASTQ path, replace the inline block with:
```rust
let mut writer = fastq_writer(&cfg, out_fmt)?;
```
(Leave the rest of that path — `reader_from`, `run_fastq`, `flush`, `drop`, reporting — unchanged.)

- [ ] **Step 2: Add the `is_dir()` branch + `run_folder`**

At the very top of `run()`, after computing `in_path`, add:
```rust
if let Some(p) = in_path {
    if p.is_dir() {
        return run_folder(p, &cfg);
    }
}
```

Add `run_folder`:
```rust
/// Folder-merge mode: `-i <dir>`. Classify the directory into one format family,
/// then merge all its read files into a single trimmed output using the same
/// pipelines as the single-file path.
fn run_folder(dir: &std::path::Path, cfg: &Config) -> anyhow::Result<()> {
    use io::Format;

    let (family, paths) = io::dir::classify(dir)?;
    let family_fmt = match family {
        io::dir::Family::Fastq => Format::Fastq,
        io::dir::Family::Bam => Format::Bam,
    };
    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), family_fmt));

    match family {
        io::dir::Family::Fastq => {
            if matches!(out_fmt, Format::Bam) {
                anyhow::bail!(
                    "cross-format conversion (FASTQ folder to BAM) is not supported in v1"
                );
            }
            let mut writer = fastq_writer(cfg, out_fmt)?;
            let records = io::dir::fastq_records(&paths);
            let stats = pipeline::run_fastq(records, &mut writer, cfg)?;
            writer.flush()?;
            drop(writer);
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            Ok(())
        }
        io::dir::Family::Bam => {
            if !matches!(out_fmt, Format::Bam) {
                anyhow::bail!(
                    "cross-format conversion (BAM folder to FASTQ) is not supported in v1"
                );
            }
            let (header, records) = io::dir::bam_reader(&paths)?;
            let out_header = provenance_header(header);
            let mut writer = io::bam::writer(cfg.io.output.as_deref(), &out_header)?;
            let stats = pipeline::run_bam(&out_header, records, &mut writer, cfg)?;
            writer.try_finish()?;
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            Ok(())
        }
    }
}
```

- [ ] **Step 3: Write failing integration tests `tests/folder.rs`**

```rust
use std::fs::File;
use std::path::Path;

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;

fn chopping() -> Command {
    Command::cargo_bin("chopping").unwrap()
}

#[test]
fn folder_merge_fastq_sorted_and_ignores_non_read_files() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("a.fastq"), "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    std::fs::write(dir.path().join("b.fastq"), "@r2\nTTTTGGGG\n+\nIIIIIIII\n").unwrap();
    std::fs::write(dir.path().join("sequencing_summary.txt"), "junk\n").unwrap(); // ignored
    let out = dir.path().join("merged.fastq");

    chopping()
        .arg("-i").arg(dir.path())
        .arg("-o").arg(&out)
        .args(["-H", "2", "-T", "2", "-t", "1"]) // -t 1 => deterministic order
        .assert()
        .success();

    // Sorted: a.fastq then b.fastq. Head/tail-crop 2 => GTAC, TTGG.
    let got = std::fs::read_to_string(&out).unwrap();
    assert_eq!(got, "@r1\nGTAC\n+\nIIII\n@r2\nTTGG\n+\nIIII\n");
}

fn write_ubam(path: &Path, name: &[u8], seq: &[u8], quals: Vec<u8>) {
    let header = noodles_sam::Header::default();
    let mut w = bam::io::Writer::new(File::create(path).unwrap());
    w.write_header(&header).unwrap();
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(name.into());
    *rec.sequence_mut() = seq.to_vec().into();
    *rec.quality_scores_mut() = quals.into();
    w.write_alignment_record(&header, &rec).unwrap();
    w.try_finish().unwrap();
}

#[test]
fn folder_merge_bam_two_files() {
    let dir = tempfile::tempdir().unwrap();
    write_ubam(&dir.path().join("a.bam"), b"r1", b"ACGTACGT", vec![40; 8]);
    write_ubam(&dir.path().join("b.bam"), b"r2", b"TTTTGGGG", vec![40; 8]);
    let out = dir.path().join("merged.bam");

    chopping()
        .arg("-i").arg(dir.path())
        .arg("-o").arg(&out)
        .args(["-H", "2", "-T", "2", "-t", "1"])
        .assert()
        .success();

    // Read the merged BAM back: 2 records, @PG chopping present.
    let mut r = bam::io::Reader::new(File::open(&out).unwrap());
    let hdr = r.read_header().unwrap();
    assert!(
        hdr.programs().roots().any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"chopping"),
        "expected @PG chopping in merged header"
    );
    let mut count = 0usize;
    let mut buf = RecordBuf::default();
    while r.read_record_buf(&hdr, &mut buf).unwrap() != 0 {
        count += 1;
    }
    assert_eq!(count, 2);
}

#[test]
fn empty_folder_errors() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("notes.txt"), "x").unwrap();
    chopping()
        .arg("-i").arg(dir.path())
        .assert()
        .failure()
        .stderr(predicates::str::contains("no FASTQ or BAM"));
}
```

Add `predicates` if not already a dev-dependency (it is, from Plan 1).

- [ ] **Step 4: Run tests to verify (RED then GREEN)**

Run: `cargo test --test folder` — before Step 1/2 wiring it fails to compile / errors; after, it passes.
Then full `cargo test` and `cargo clippy --all-targets -- -D warnings` (clean).

- [ ] **Step 5: Update `README.md`**

Add a short "Folder input" subsection under usage: `-i` accepts a single file **or a directory**. When given a directory, chopping merges every read file directly in it (`.fastq`/`.fq`/`.fastq.gz`/`.fq.gz`, or `.bam`) — in sorted filename order — into one trimmed output (`-o` file or stdout). The folder must be one format (all FASTQ-family, or all BAM); non-recursive; a mixed or empty folder is an error. Example:
```
chopping -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --trim-qual 10
```

- [ ] **Step 6: Commit**

```bash
git add src/lib.rs tests/folder.rs README.md
git commit -m "feat: -i accepts a directory (folder-merge) + integration tests + README"
```

---

## Self-Review

**Spec coverage:**
- `-i <dir>` auto-detected via `is_dir()` → Task 2 Step 2. ✅
- Single flat folder, non-recursive, immediate children → `classify` (Task 1). ✅
- Read-file selection via `from_extension`, non-read ignored → `classify`. ✅
- Sorted order → `classify` `paths.sort()`; asserted in Task 2 FASTQ test. ✅
- Homogeneity + mixed/empty errors → `classify` match arms + tests (Tasks 1, 2). ✅
- FASTQ merge → chained `fastq_records` + `run_fastq` (Tasks 1, 2). ✅
- BAM merge, first-file header + `@PG`, `try_finish` → `bam_reader` + `run_folder` Bam arm (Tasks 1, 2). ✅
- Output format from `-o` else family → `run_folder` `resolve_output`. ✅
- DRY writer via `fastq_writer` shared with single-file path → Task 2 Step 1. ✅
- Tests: classify units + fastq/bam integration + empty error → Tasks 1, 2. ✅

**Placeholder scan:** none. All code blocks complete.

**Type consistency:** `Family{Fastq,Bam}`, `classify -> (Family, Vec<PathBuf>)`, `fastq_records(&[PathBuf]) -> Box<dyn Iterator<Item=anyhow::Result<ReadRecord>> + Send>`, `bam_reader(&[PathBuf]) -> anyhow::Result<(sam::Header, Box<dyn Iterator<Item=anyhow::Result<RecordBuf>>>)>` are used identically in Task 2. `fastq_writer(&Config, io::Format) -> anyhow::Result<BufWriter<Box<dyn Write + Send>>>` matches both call sites. `run_folder(&Path, &Config)` matches the branch call. The BAM `bam_reader` iterator is intentionally NOT `+ Send` (matches `run_bam`'s non-Send bound); the FASTQ one IS `+ Send` (matches `run_fastq`).
