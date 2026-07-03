# chopping — tag-aware long-read trimmer (v1 design)

- **Date:** 2026-07-03
- **Status:** Approved for v1
- **Repo:** https://github.com/erdikilic/chopping (private, empty — real initial commit + 0.1.0 release deferred)

## Goal

A well-tested CLI for trimming and splitting long reads (ONT/PacBio) that handles
both FASTQ and BAM. The differentiator: when a read is trimmed or split, its
base-modification tags (MM/ML — methylation and other epigenetic marks) are
**recomputed** so they stay correct on the shortened read, instead of being copied
through verbatim and silently corrupted.

## Background

The design harvests the algorithm core of `chopper` (MIT) but is built fresh.
Chopper's key idea worth carrying: **trimming returns coordinate intervals**
(`Vec<(start, end)>`) on the original read, never new sequences. That interval
model is what makes tag reconstruction tractable — a mod call's fate reduces to
"is its position inside a surviving interval, and what is its new offset."
Chopper itself is FASTQ/String-only and has no tag concept; it demonstrates the
bug we fix (it copies FASTQ-header tags through unchanged on trim).

## Scope

**In v1:**
- FASTQ(.gz) → FASTQ(.gz): length + mean-quality + GC filtering, trim/split.
- Unaligned BAM (uBAM) → uBAM: the same filtering + trim/split, **plus** MM/ML
  (and `MN`) recomputed on every output read.
- Same-format only (FASTQ→FASTQ, uBAM→uBAM).

**Deferred (explicitly out of v1):**
- Aligned BAM (CIGAR/POS recompute, hard/soft-clip, split-into-supplementary).
  v1 detects aligned records and refuses them with a clear error rather than
  silently mistreating them.
- minimap2 contamination filtering.
- Cross-format conversion (BAM→FASTQ, FASTQ→BAM).
- Reverse-complement bookkeeping for `-` strand mods beyond what uBAM (original
  basecall orientation) requires.

## Architecture — Approach A: format-neutral record + shared interval pipeline

One record model flows through pure, format-agnostic filter and trim stages that
produce intervals; a format-specific reconstruction step turns `(record, interval)`
into an output record. The MM/ML codec is an isolated module invoked only on the
BAM reconstruction path. Trim/filter logic is written once and shared by both
formats. Every unit is independently testable.

### Module layout

```
src/
  main.rs         clap entry, wires CLI -> pipeline
  cli.rs          arg struct + validation
  record.rs       ReadRecord { name, seq, qual, tags: Option<BamTags> }
  filter.rs       length / mean-qual / GC  -> bool  (format-agnostic)
  trim/
    mod.rs        TrimStrategy trait -> Vec<(start, end)>
    strategies.rs ported chopper algorithms
  mods/
    parse.rs      MM:Z / ML:B,C / MN  ->  typed model
    reconstruct.rs  (model, interval) -> sliced model
    serialize.rs  typed model -> MM:Z / ML:B,C / MN
  io/
    mod.rs        format detection, Reader/Writer enums, aligned-read refusal
    fastq.rs      seq_io reader/writer (+ gz)
    bam.rs        noodles reader/writer, raw-tag access
  pipeline.rs     filter -> trim -> reconstruct -> write; parallelism
```

### Record model & data flow

`ReadRecord { name, seq, qual, tags: Option<BamTags> }`. Per read:
1. Run filters cheapest-first, bail early.
2. If kept, run the selected `TrimStrategy` → `Vec<(start, end)>`.
3. For each interval, **reconstruct** an output record.

Reconstruction is the only format-aware step: FASTQ slices `seq`/`qual` and
suffixes the name on splits (`_segment_N`, chopper's convention); uBAM does the
same slice **and** rebuilds MM/ML/MN. No trimmer or filter knows the format.

### Trimmers & filters

Ported from chopper as pure functions over the quality/sequence slice, returning
intervals:
- `fixed-crop` (`--headcrop` / `--tailcrop`)
- `trim-by-quality` (trim ends below `--cutoff`)
- `best-read-segment` (Mott algorithm, `--cutoff`)
- `split-by-low-quality` (`--cutoff`, `--split-window`)

Filters: length (`--minlength`/`--maxlength`), mean quality
(`--minqual`/`--maxqual`), GC (`--mingc`/`--maxgc`). All single-pass, shared.

## MM/ML/MN reconstruction (the core)

We **do not** use noodles' typed base-mod parser — it drops whole-read mods on
certain "unspecified" mod codes (previously verified; htslib is the oracle).
Instead we read the **raw** `MM:Z` string and `ML:B,C` array directly off the
noodles record and run our own codec.

Tag model:
- `MM:Z` is a set of groups like `C+m?,5,12,0;` — a fundamental base (`C`), a
  strand (`+`/`-`), one or more mod codes (`m`), an optional `.`/`?` status flag,
  then skip-counts over occurrences of that fundamental base.
- `ML:B,C` holds one probability byte per listed modified position, in MM order
  (multiple mod codes on one position produce multiple ML bytes).
- `MN:i` records the SEQ length the tags were computed against.

Reconstruction for a surviving interval `[start, end)` on the read. uBAM SEQ is
in original basecall orientation, so **no reverse-complement bookkeeping in v1**:

1. For each MM group, walk occurrences of its fundamental base along the original
   SEQ, using skip-counts to mark which occurrences are *modified* (each consumes
   the next ML byte(s)).
2. Keep only modified occurrences whose absolute coordinate falls in
   `[start, end)`; drop the rest and their ML bytes.
3. Recompute skip-counts relative to the fundamental-base occurrences that lie
   **inside** the window (renumbered from the window start).
4. Emit the rebuilt `MM:Z`, the filtered `ML` array (same order), and
   `MN = end - start`.

Groups that end up empty inside the window are dropped. The `.`/`?` status flag
and mod codes are carried through per group untouched. Other aux tags are passed
through unchanged.

## I/O & aligned-read refusal

- FASTQ via `seq_io`; `.gz` via a gzip layer on reader/writer.
- BAM via `noodles-bam` + `noodles-bgzf` (libdeflate).
- On BAM read, any record whose unmapped flag is **clear** (aligned: has
  refID/CIGAR) → hard error naming the read and stating aligned BAM is not yet
  supported, rather than silently mistreating it.
- BAM output copies the input header and appends an `@PG` line for provenance.

## Pipeline & parallelism

Chopper's proven shape: one reader, a rayon work pool doing
filter/trim/reconstruct, a dedicated writer thread. Multi-threaded output is
**unordered** (reads are independent) for throughput; `--threads 1` gives
deterministic output, which golden tests use. Principle (from benchmarking):
don't parallelize parsing, parallelize the per-read work.

## CLI

`-i/--input` (or stdin), `-o/--output` (or stdout), `--in-format`/`--out-format`
overrides, `-t/--threads`. Filtering: `-l/--minlength`, `--maxlength`,
`-q/--minqual`, `--maxqual`, `--mingc`, `--maxgc`. Trimming:
`--trim-approach {fixed-crop|trim-by-quality|best-read-segment|split-by-low-quality}`
with `--cutoff`, `--headcrop`, `--tailcrop`, `--split-window`. Mod recompute is
automatic on BAM (no flag).

## Testing (first-class — this is a "well-tested" tool)

- **Trim/filter:** chopper's unit tests ported verbatim (proves algorithm
  parity) plus new edge cases.
- **MM/ML codec:** hand-built examples straight from the hts-specs; plus an
  **oracle test vs htslib** (`rust-htslib` as a dev-dependency, test-only) using
  **decode equivalence** rather than byte-equality — decode the mods from our
  trimmed output with htslib and assert the per-position set equals the original
  mods filtered to the interval and offset by `start`. Sidesteps MM's multiple
  valid encodings and directly catches a wrong codec.
- **Integration/golden:** real HG002 uBAM subset round-tripped at `--threads 1`;
  FASTQ goldens compared against chopper for the ported strategies.
- **CLI:** `assert_cmd`/`trycmd` for arg wiring, stdin/stdout, aligned-read
  refusal.

## Dependencies

`seq_io`, `noodles-{bam,sam,bgzf}`, `flate2` (gz), `clap` (derive),
`rayon` + `crossbeam`, `anyhow` (binary) + `thiserror` (library). `rust-htslib`
is a **dev-dependency** for the oracle only. Rust 2024 edition.

## Error handling

Fail fast with context. Malformed records and aligned BAM are hard errors, not
silent skips.
