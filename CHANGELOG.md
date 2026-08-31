# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Adapter and primer sequences may use the full IUPAC alphabet. A degenerate
  primer now matches every variant its wobble positions cover, instead of being
  skipped as "non-ACGT". `U` folds to `T`; non-nucleotide characters are still
  skipped with a warning, and a pattern averaging two or more bases per position
  is searched but flagged as near-wildcard. An ambiguity code in a *read* is
  treated as a mismatch instead, so it costs error budget rather than matching
  for free: a stray `N` in a real adapter still matches, a run of them is not
  excised as one.
- `--summary-json <PATH>`: writes a machine-readable JSON summary of the run,
  covering the resolved settings (`params`) and the read, base, and per-reason
  segment-drop counters. Written on every dispatch path, folder merges included,
  and regardless of `--quiet` or the log level, so a workflow manager always gets
  the file it asked for. `schema_version` is bumped only when an existing field
  changes meaning or disappears.
- A `whittle.1` man page, checked into `man/` and shipped in every release
  tarball. It is rendered from the live CLI definition by
  `cargo run --example gen-man`; `clap_mangen` is a dev-dependency, so the
  shipped binary is unchanged.

### Changed
- `-vv` now reports the per-read decisions it always claimed to: which adapter
  matched where and at what cost, what that made whittle do, and why each segment
  was kept or dropped, each line attributed to the read that produced it. There
  were previously no trace-level events at all, so `-vv` was indistinguishable
  from `-v`. `-v` gains the resolved thread budget and the run counters. Log
  events carry structured fields rather than preformatted prose. Measured at no
  cost to throughput at the default level.
- BAM to FASTQ no longer re-parses and re-serializes `MM`/`ML` for a record whose
  window is not being trimmed. Over the full window the reconstruction is the
  identity, so the source bytes are reused after an allocation-free `ML` length
  check. About 29% less CPU on an untrimmed conversion, with byte-identical
  output.

### Fixed
- Every base-modification record missed the untrimmed fast path, because the
  `MN` consistency check accepted only the `i` (Int32) subtype while SAM writes
  integers at the smallest width that fits, so dorado emits `MN:S`. A filter-only
  or pass-through BAM run therefore rebuilt `MM`/`ML` for every record and
  rewrote `MN` at the wider subtype. Recognizing every integer width cuts about
  27% of the CPU time from those runs and leaves `MN` at the width it arrived
  with. Decoded tag values are unchanged.
- Adapter trimming aborted the whole run with a panic, leaving a truncated
  output file, when any read contained an ambiguity code (`N`, `Y`, `R`, ...)
  near an end. Sassy's `Dna` profile panics during traceback on non-ACGT text,
  and the searcher fell back to exactly that profile whenever it detected one.
  Every searcher now uses the IUPAC profile, which is what tolerates ambiguity
  codes. Output on pure-ACGT input is byte-identical; adapter trimming is about
  19% slower.
- uBAM records flagged reverse-complemented (`0x4|0x10`) are refused rather than
  trimmed: htslib decodes their `MM` right to left with complemented bases and
  their stored `SEQ` is reverse-complemented, so trimming cropped the wrong ends
  and relocated every call.
- Records carrying the legacy lowercase `Mm`/`Ml` tags (old guppy, megalodon)
  are refused rather than copied through untouched onto a trimmed sequence.
  htslib still reads that spelling, so the calls decoded fine elsewhere while
  pointing at the wrong bases.
- `--summary-json` could destroy the run's own output when reads went to a
  redirected stdout (`whittle --summary-json out.fastq > out.fastq`): the
  collision guard only covered the `-o` form.
- Folder mode's progress bar was pinned at 0%, because nothing on that path
  counted input bytes. Folder readers now wrap their files the way the
  single-file path does.
- A mistyped `--summary-json` path now fails during setup instead of after the
  reads are written, and `-o /dev/null --summary-json /dev/null` is no longer
  refused as a self-overwrite.
- Folder mode (`-i <dir>`) was missing two advisories the single-file path
  emits: the `--out-format` extension mismatch and the no-trimming no-op
  warning. It also rendered a bare spinner instead of a progress bar with a
  percentage and ETA, despite having already measured the input.
- Reverse-strand `MM` groups were counted against the complement of the
  fundamental base, so every call on a `-` strand group was relocated onto a
  different base, some were dropped, and the `ML` bytes shifted relative to
  their positions. htslib counts the literal base; whittle now matches it.
  `U` is folded to `T` (BAM's SEQ encoding has no `U`) and `N` counts every
  base, both of which previously deleted the whole group.
- BGZF compression and decompression ran single-threaded regardless of `-t`.
  The `noodles-bgzf` 0.48 API ignored its worker count and used Rayon's global
  pool, which whittle configured; 0.51 takes an explicit count and builds a
  private pool, defaulting to one worker. Measured 3.5x faster BAM-to-BAM at
  `-t 16`.
- whittle no longer overwrites a file it is also reading: `--summary-json`
  pointing at the input, the output, or a folder member; `-o` pointing at
  `--adapter-fasta`; or `-o` naming the file redirected onto stdin, which the
  path-based check could not see.
- A cyclic `@PG` chain that does not return to its entry node (`pgA` to `pgB`
  to `pgC` to `pgB`) hung the run at 100% CPU. The chain walk now tracks
  visited IDs.
- `--summary-json` reported `elapsed_seconds: null` under `--quiet`, and was
  silently skipped under `--adapter-infer report`, which now says the flag is
  ignored.
- Argument-parsing diagnostics bypassed the log level filter, printing even
  under `--quiet` without the standard prefix and ahead of the version line.

### Changed
- `--summary-json` reports both adapter counts: `params.adapters.configured` is
  the set asked for, `params.adapters.count` the set actually trimmed against
  after presence detection or inference. The startup banner prints the former,
  so the two figures no longer look like a contradiction.
- Release tarballs now contain a versioned directory holding the binary, `man/`,
  and the README, CHANGELOG, and LICENSE, instead of a bare `whittle`
  executable.
- Documentation split: the README keeps the overview, install, and quick start;
  the full flag reference moved to `docs/cli.md`, the tag machinery to
  `docs/tags.md`, adapter trimming to `docs/adapters.md`, and the contributor
  policy to `CONTRIBUTING.md`.
- Comments, doc comments, and user-facing strings standardized to American
  English prose without em or en dashes.
- `noodles-bam`, `noodles-sam`, and `noodles-bgzf` updated to 0.95, 0.90, and
  0.51; `clap`, `anyhow`, `thiserror`, `crossbeam-channel`, `gzp`, `bstr`,
  `aho-corasick`, and `jiff` to their latest patch releases.

## [0.1.1] - 2026-07-14

### Added
- Startup AVX2 capability check: on x86-64 builds compiled with AVX2 (the
  default, via `target-cpu=x86-64-v3`), whittle now verifies at launch that the
  running CPU supports AVX2 and exits with a clear message instead of crashing
  with an illegal instruction on older CPUs.

## [0.1.0] - 2026-07-14

### Added
- Long-read (ONT / PacBio) trimming for FASTQ, gzipped FASTQ, and
  unaligned BAM: fixed head/tail crop, quality trimming (Mott, best-segment,
  quality-split), and adapter trimming with interior-chimera splitting.
- Length / quality / GC filtering, applied per surviving segment after trimming.
- Trim-aware rewriting of base-modification (`MM`/`ML`/`MN`) and per-base
  kinetics/signal tags, so every trim and split keeps its tags valid.
- BAM→FASTQ conversion with selectable aux-tag carry-through.
- Folder-merge mode, parallel processing with a workload-aware thread budget,
  and a progress/summary UI.

[Unreleased]: https://github.com/erdikilic/whittle/compare/0.1.1...HEAD
[0.1.1]: https://github.com/erdikilic/whittle/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/erdikilic/whittle/releases/tag/0.1.0
