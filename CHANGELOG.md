# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `--remove-tag <TAG>` and `--strip-kinetics`: remove auxiliary tags from every
  output record. `--remove-tag` names one two-character tag and is repeatable;
  `--strip-kinetics` removes the nine per-base kinetics and alignment-count
  arrays (`ip`, `pw`, `fi`, `fp`, `ri`, `rp`, `sa`, `sm`, `sx`) in one flag, and
  both fill the same removal set. Removal runs after the rewrite of the tags
  whittle keeps in register, so removing `MM` still leaves a rebuilt `ML`/`MN`
  and removing one per-base array still slices the others. It applies on every
  BAM output path, including the untrimmed fast path, whose records are rebuilt
  rather than passed through when tags are removed, and to the tags carried into
  a BAM-to-FASTQ header. The resolved set is reported under
  `params.remove_tags`. BAM input only.
- `--trim-barcodes`: removes the barcode spans dorado recorded in the `bi` aux
  tag. The positions are read rather than the sequences searched for, so the cut
  is the one `dorado demux` would have made, and it runs through the same
  trimming machinery as every other stage, keeping `MM`/`ML`/`MN`, per-base
  kinetics and the ONT move table in register. It is the outermost stage, so
  `--head-crop` counts from the first base after the front barcode. A barcode dorado
  did not find leaves that end of the read alone; a `bi` that does not describe a
  window inside the read leaves it untrimmed and is counted under
  `warnings.barcode_tag_malformed_reads`. BAM input only.
- `--progress {auto,bar,plain,none}` selects how progress is reported,
  independently of the log level. `--progress none` keeps the banner and the run
  summary while reporting nothing in flight, which suits a pipeline log;
  `--quiet` drops the summary as well and outranks it.
- Adapter and primer sequences may use the full IUPAC alphabet. A degenerate
  primer now matches every variant its wobble positions cover, instead of being
  skipped as "non-ACGT". `U` folds to `T`; non-nucleotide characters are still
  skipped with a warning, and a pattern averaging two or more bases per position
  is searched but flagged as near-wildcard. An ambiguity code in a read is
  treated as a mismatch instead, so it costs error budget rather than matching
  for free: a stray `N` in a real adapter still matches, a run of them is not
  excised as one.
- `--summary-json <PATH>`: writes a machine-readable JSON summary of the run,
  covering the resolved settings (`params`) and the read, base, and per-reason
  segment-drop counters. Written on every dispatch path, folder merges included,
  and regardless of `--quiet` or the log level. `schema_version` is bumped only
  when an existing field changes meaning or disappears.
- A `whittle.1` man page, checked into `man/` and shipped in every release
  tarball. It is rendered from the live CLI definition by
  `cargo run --example gen-man`; `clap_mangen` is a dev-dependency, so the
  shipped binary is unchanged.
- `--ordered` writes records in input order when running with more than one
  thread, so output is byte-identical between runs and to a single-threaded run.
  Without it records are written as they finish, which is faster and uses less
  memory; a BAM written that way carries `SO:unsorted` in its header.

### Changed
- The progress bar shows the output count beside the input count, so a filter
  discarding everything is visible while the run is going rather than only in
  the summary. The percentage is right-aligned so the bar no longer shifts
  sideways as it passes 9% and 99%, and the first frame carries the same fields
  as every later one. Still ASCII, so it renders the same over SSH, `screen` and
  a non-UTF-8 console.
- `-vv` reports the per-read decisions: which adapter matched where and at what
  cost, the resulting action, and why each segment was kept or dropped, each
  line attributed to its read. There were previously no trace-level events at
  all, so `-vv` was indistinguishable from `-v`. `-v` gains the resolved thread
  budget and the run counters. Log events carry structured fields rather than
  preformatted prose. Measured at no cost to throughput at the default level.
- BAM to FASTQ no longer re-parses and re-serializes `MM`/`ML` for a record whose
  window is not being trimmed. Over the full window the reconstruction is the
  identity, so the source bytes are reused after an allocation-free `ML` length
  check. About 29% less CPU on an untrimmed conversion, with byte-identical
  output.
- `--summary-json` reports both adapter counts: `params.adapters.configured` is
  the set asked for, `params.adapters.count` the set trimmed against after
  presence detection or inference. The startup banner prints the former, so the
  two figures no longer look like a contradiction.
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
- The gzip decoder is zlib-rs rather than zlib-ng, which reads compressed FASTQ
  about 7 percent faster and drops a C build dependency from the decode path.

### Fixed
- `--summary-json` wrote its `command` field with the banner's `Command: ` label
  glued to the front of the value, so a consumer re-running the recorded
  invocation had to strip it first. The label now belongs to the banner line and
  the JSON carries the bare, shell-quoted argv it always documented.
- `thread_budget` allocated zero render workers at exactly `-t 3` with parallel
  decode and uncompressed output, so the banner reported `trim 0` and the
  workflow silently fell back to its own default. Every other thread count is
  unchanged; a property test now checks that no stage is ever allocated zero.
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
  htslib still reads that spelling, so the calls decoded correctly elsewhere while
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
- PacBio reverse-strand kinetics (`ri`, `rp`) are sliced from the far end of the
  array, matching the last-base-first layout the PacBio BAM specification
  defines; trimmed HiFi reads previously carried shifted reverse IPD and pulse
  width values.
- The `qs` refresh applies only to a float mean-quality tag. PacBio `qs:i` and
  `qe:i` query coordinates pass through unchanged instead of being overwritten.
- Under `--adapter-preset ont`, reads shorter than about 175 bp that carry a
  normal front or rear adapter keep their insert. The rear catalog entries are
  reverse complements of the front entries, and a hit covered by both end zones
  was placed by its tag rather than its position, so the two trims met and the
  read was dropped.
- Adapters with IUPAC ambiguity codes split chimeras. The exact-seed filter
  expands ambiguity codes and the complement table covers the full IUPAC
  alphabet; previously such adapters trimmed at the ends but never matched an
  interior hit.
- A truncated or unreadable gzip input under `-t 2` or more exits with status 1
  and a message instead of panicking after writing partial output. A failing
  run also stops reading its input at the first error.
- An MM group whose positions all fall outside the kept window is emitted as an
  empty group rather than removed, so implicit-mode bases stay canonical. An
  `MN` that disagrees with the sequence length, an `ML` whose length disagrees
  with `MM`, an `MM` that does not parse to its end, or a non-`B:C` `ML` removes
  the modification block and is counted under `warnings.malformed_mod_reads` in
  the summary and in the end-of-run advisory, instead of being repaired or
  dropped silently.
- Quality bytes outside the Phred+33 range are an error naming the record and
  the byte, instead of being rewritten to Q0.
- A bgzip-compressed file named `.fastq.gz` is detected by its block header and
  read through the multithreaded BGZF decoder.
- Directory input skips hidden files and sorts members in natural order, so
  `run_2` precedes `run_10`.
- `-` names stdin for `-i` and stdout for `-o`. `-t 0` is rejected, quality
  bounds must be finite and non-negative, `--qual-split-window` and the adapter
  tuning flags require the flag they modify,
  and `--progress` conflicts with `--quiet`.
- A `WHITTLE_LOG` value that does not parse falls back to the verbosity level
  and is reported, instead of silencing every log line including errors.
- A missing input file error names the path; a closed downstream pipe exits
  quietly with status 0.
- The end-of-run `Completed` line prints after the summary JSON is written, so a
  failed write is not preceded by a success line.
- PacBio records are recognized by an integer `qs` or a PacBio read name, and
  follow the PacBio BAM specification: `qs`/`qe` are shifted to the kept
  window, split segments are named `{movie}/{zmw}/ccs/{qStart}_{qEnd}` (with a
  by-strand `fwd`/`rev` component kept before the interval) so pbbam parses
  them, `rn` (reverse passes) is left untouched, `du:Z` from pbmarkdup passes
  through, the run-length `sa` coverage array is re-sliced, `sm`/`sx` are
  sliced per base, and the `ds`/`ls` undo blobs are removed from trimmed reads
  and counted under `warnings.undo_tags_dropped_reads`, since `skera undo` and
  `lima-undo` would otherwise rebuild the wrong read.
- ONT split segments carry `pi` (the parent read id) on both outputs without
  `--update-moves`, `me:i:0` on every segment and `er:Z:unknown` on every
  segment but the last, matching dorado. Under `--update-moves` a split
  recomputes `st` and `du` from the parent's sample rate instead of dropping
  them, and `sp` includes the parent's trimmed-sample offset as dorado's
  splitter does.
- The `@PG` record carries the command line (`CL`); a repeat run gets a
  distinct ID and a `PP` link to the previous record.

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
- BAM-to-FASTQ conversion with selectable aux-tag carry-through.
- Folder-merge mode, parallel processing with a workload-aware thread budget,
  and a progress/summary UI.

[Unreleased]: https://github.com/erdikilic/whittle/compare/0.1.1...HEAD
[0.1.1]: https://github.com/erdikilic/whittle/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/erdikilic/whittle/releases/tag/0.1.0
