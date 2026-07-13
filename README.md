<div align="center">

# whittle

**A fast, tag-aware trimmer for long-read FASTQ and unaligned BAM.**

It rewrites `MM`/`ML` base-modification tags on every trim and split, so a trimmed read's methylation calls never drift out of register with its sequence.

[![CI](https://github.com/erdikilic/whittle/actions/workflows/ci.yml/badge.svg)](https://github.com/erdikilic/whittle/actions/workflows/ci.yml)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
![Rust](https://img.shields.io/badge/rust-2024%20edition-000000?logo=rust&logoColor=white)
![Long reads](https://img.shields.io/badge/long--reads-ONT%20%7C%20PacBio-1f6feb)

</div>

whittle filters and trims long reads (ONT, PacBio) in FASTQ, gzip/BGZF-compressed FASTQ, and unaligned BAM. It does the usual length/quality/GC filtering and head/tail/quality/adapter trimming — and, on uBAM, it reconstructs the base-modification (`MM`/`ML`/`MN`), per-base kinetics, and ONT signal tags so every output read stays valid instead of silently decoding to nonsense.

## Highlights

- **Correct modification tags.** `MM`/`ML`/`MN` are rebuilt for every trimmed or split uBAM read, verified against an independent `htslib` decoder.
- **Trim-aware tags.** Per-base kinetics (`ip`/`pw`/…) are sliced with the sequence; ONT signal tags (`mv`/`ts`/`ns`/…) are dropped, or rewritten dorado-style with `--update-moves`.
- **Adapter trimming.** Terminal trimming and interior chimera splitting, from a built-in ONT catalog, your own FASTA, or ab-initio discovery.
- **Formats.** FASTQ, gzip/BGZF-compressed FASTQ, and unaligned BAM, plus BAM→FASTQ conversion. Formats are auto-detected, including BGZF FASTQ or BAM piped over stdin.
- **Fast, self-contained.** Multithreaded throughout with a workload-aware thread budget, and no external `htslib` needed to build or run.

## Install

### Prebuilt binaries

Download a binary for your platform — Linux and macOS, x86-64 and arm64, glibc or static musl — from the [Releases](https://github.com/erdikilic/whittle/releases) page and put it on your `PATH`.

### From source

```bash
git clone https://github.com/erdikilic/whittle
cd whittle
cargo build --release   # -> target/release/whittle
```

### From crates.io

The adapter search ([`sassy`](https://crates.io/crates/sassy)) needs AVX2 on x86-64, and a `cargo install` doesn't inherit this repo's build config, so pass the flag yourself:

```bash
RUSTFLAGS="-C target-cpu=x86-64-v3" cargo install whittle
```

Building requires Rust 1.91 or newer. No external `htslib` is needed: BAM I/O uses the `libdeflate` backend of `noodles-bgzf` (an optional `rust-htslib` test dependency is dev-only).

## Quick start

Trim a FASTQ file — crop 20 bp off each end, quality-trim below Q8, drop reads under 500 bp or Q10, on 8 threads:

```bash
whittle -i reads.fastq.gz -o trimmed.fastq.gz -H 20 -T 20 --qual-trim 8 -l 500 -q 10 -t 8
```

Trim an unaligned BAM, splitting at low-quality runs, with modification tags rebuilt for every output read:

```bash
whittle -i reads.ubam.bam -o trimmed.ubam.bam -H 10 -T 10 -l 1000 --qual-split 9 --qual-split-window 50
```

## Usage

whittle reads from `-i`/`--input` (or stdin) and writes to `-o`/`--output` (or stdout). The format is taken from the file extension, sniffed from the first bytes of a stream, or forced with `--in-format`/`--out-format {fastq,fastq-gz,fastq-bgz,bam}`.

Output is plain FASTQ by default and is never auto-compressed — a `.gz` or `.bgz` input does not imply compressed output. Compressed output happens only when requested with a `.gz`/`.bgz` path or the corresponding format flag, and is written by a parallel encoder using `-t` threads. BGZF FASTQ input is also decompressed block-parallel; ordinary gzip remains a serial input format.

### How trimming works

Operations run in a fixed order, and the filters apply to whatever is left:

1. **Fixed crop** — `-H`/`--head-crop` and `-T`/`--tail-crop` remove a set number of bases from each end.
2. **Adapters** — terminal adapters are trimmed and interior adapters split the read (see [Adapter trimming](#adapter-trimming)).
3. **Quality** — one of `--qual-trim`, `--qual-best-segment`, or `--qual-split` (mutually exclusive).
4. **Filter** — each surviving segment must pass `-l`/`-L` (length), `-q`/`-Q` (quality), and `-g`/`-G` (GC).

When a read splits, each segment is filtered on its own and named `<read>_segment_N` (1-based), so `-l` is a post-trim, per-segment minimum.

### Formats

| input → | FASTQ | FASTQ.gz | FASTQ.bgz | BAM |
|---|:---:|:---:|:---:|:---:|
| FASTQ / FASTQ.gz / FASTQ.bgz | ✅ | ✅ | ✅ | ❌ |
| unaligned BAM | ✅ | ✅ | ✅ | ✅ |

With no output extension or `--out-format`, output mirrors the input, except compressed FASTQ input defaults to plain FASTQ. FASTQ→BAM is not supported — there's no header to build a BAM record from. BGZF streams are identified by their decompressed payload, so piped FASTQ.bgz and `samtools view -b … | whittle` need no hint.

On BAM→FASTQ, aux tags are written into the FASTQ header tab-delimited (the `samtools fastq -T` convention). `--fastq-tags` chooses which: `all` (default), `none`, or a list like `MM,ML,RG`. `MM`/`ML`/`MN` are reconstructed for the trimmed segment, per-base tags are sliced, and the rest are copied verbatim.

### Folder input

`-i` also accepts a directory: whittle merges every read file directly inside it, in sorted filename order, into one output. The folder must be a single format (all FASTQ-family or all BAM); subdirectories are ignored, and a mixed or empty folder is an error.

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --qual-trim 10
```

### Options

| Flag | Meaning |
|---|---|
| `-i, --input <PATH>` | Input file or directory (omit for stdin) |
| `-o, --output <PATH>` | Output file (omit for stdout) |
| `--in-format`, `--out-format {fastq,fastq-gz,fastq-bgz,bam}` | Force a format instead of detecting it |
| `--fastq-tags {all,none,LIST}` | Aux tags to carry into FASTQ headers on BAM→FASTQ (default `all`) |
| `-c, --compression-level <0-9>` | DEFLATE level for compressed output (default 6); ignored for plain FASTQ |
| `-t, --threads <N>` | Worker threads (default: all detected CPUs, clamped to that max) |
| `-l, --min-length <N>` | Minimum length to keep, per output segment (default 1) |
| `-L, --max-length <N>` | Maximum length to keep |
| `-q, --min-qual <F>` | Minimum read quality (default 0) |
| `-Q, --max-qual <F>` | Maximum read quality (default 1000) |
| `-g, --min-gc <F>`, `-G, --max-gc <F>` | GC-fraction bounds (0–1) |
| `-m, --qual-mode {mean,arithmetic,median}` | How read quality is summarized (default `mean`, the error-probability mean) |
| `-H, --head-crop <N>`, `-T, --tail-crop <N>` | Fixed crop from each end; always runs first |
| `--qual-trim <Q>` | Trim low-quality bases from both ends down to the first base ≥ Q |
| `--qual-best-segment <Q>` | Keep only the longest contiguous run of quality ≥ Q |
| `--qual-split <Q>` | Split at low-quality (< Q) runs, keeping each surviving segment |
| `--qual-split-window <N>` | Tolerate low-quality runs shorter than N without splitting (default 1) |
| `--update-moves` | Rewrite ONT signal tags through trimming instead of dropping them (BAM→BAM) |
| `-a, --adapter-fasta <FILE>` | Adapter/primer FASTA; enables adapter trimming |
| `--adapter-preset ont` | Use the built-in ONT catalog; enables adapter trimming |
| `--adapter-error-rate <F>` | End-match tolerance as a fraction of adapter length (default 0.2) |
| `--adapter-end-size <N>` | End-zone width searched for terminal adapters (default 150) |
| `--adapter-ends-only` | Trim ends only; never split on an interior adapter |
| `--adapter-sample <N>` | Sample the first N reads for presence detection (default 0 = off; N ≥ 100) |
| `--adapter-infer` | Discover adapters de novo, then trim with the discovered set |
| `--adapter-infer-only` | Discover and print adapters, then exit without trimming |
| `-v`, `-vv` | Increase log detail (debug, trace) |
| `--quiet` | Silence progress and the summary; warnings and errors still print |

`--qual-trim`, `--qual-best-segment`, and `--qual-split` are three strategies for the same step and are mutually exclusive — pass at most one. `-H`/`-T` are independent and compose with whichever you pick.

## Base-modification tags

Long-read modification calls (5mC, 6mA, …) live per-read in `MM` (which bases are modified, encoded as skip-counts over the sequence) and `ML` (their probabilities). Trim the `SEQ`/`QUAL` without touching these and the result decodes to nonsense — the skip-counts still index the original sequence, not the trimmed one.

whittle rebuilds them as part of trimming. For every output uBAM read — cropped, quality-trimmed, or split — `MM` and `ML` are reconstructed against the output window: skip-counts renumbered, probability bytes re-sliced, and `MN` updated to match, or all three dropped if trimming removed every modified base. Everything else in the record rides through unchanged.

This is checked by decode-equivalence tests that re-decode whittle's output with `rust-htslib`'s `basemods_iter()` — a different `MM`/`ML` implementation from the one whittle writes with — and compare against the original calls restricted to the surviving window. One always runs on a synthetic fixture; another sweeps a real uBAM when you point it at one:

```bash
WHITTLE_UBAM=/path/to/real.ubam cargo test --test bam_mods_oracle -- --ignored
```

### Other trim-aware tags

Every tag indexed by base position is kept consistent, on both BAM→BAM and BAM→FASTQ:

| Tag(s) | On a trimmed read |
|---|---|
| `MM` / `ML` / `MN` | Reconstructed for the output window |
| Per-base kinetics — `ip`/`pw`/`fi`/`fp`/`ri`/`rp`, and any read-length `B` array | Sliced in lockstep with the sequence |
| ONT signal — `mv`/`ts`/`ns`/`sp`/`pi` | Dropped, or rewritten with `--update-moves` |
| Poly-A — `pa`/`pt` | Kept/shifted with `--update-moves` if the tail survives, else dropped |
| `bi` (barcode positions) | Dropped — the positions shift under a crop |
| `qs` (mean qscore) | Recomputed from the trimmed quality |
| `st`/`du` (start time / duration) | Kept on a crop, dropped on a split |
| `RG`, `ch`, `mx`, `sm`/`sd`/`sv`, … | Copied verbatim |

Under `--update-moves`, a crop slices `mv` and advances `ts`/`ns`, while a split emits dorado-style subreads (`pi` parent id, `sp` parent-signal offset, `ns` subread span, `ts` 0, `rn` -1) so the renamed segment stays locatable in POD5 for tools like Remora. BAM→FASTQ always drops the signal tags on trim (a move table in a FASTQ header is impractical). A known per-base tag whose length doesn't match the sequence (malformed input) is left untouched, with a one-line advisory.

## Adapter trimming

Off by default. Enable it with `-a`/`--adapter-fasta <FILE>` (your own sequences, one per record, each ≥ 11 bp) and/or `--adapter-preset ont` (the built-in catalog); either alone is enough, and they can be combined. Every adapter is searched on both strands, so orientation doesn't matter, and each read gets two treatments:

- **Terminal trimming** — an adapter within `--adapter-end-size` bases of an end (default 150) is trimmed off.
- **Chimera splitting** — an adapter in the interior is treated as a junction: the read splits there, the adapter is excised, and both sides are kept. `--adapter-ends-only` disables this and searches only the two end-zones.

Interior hits use half the `--adapter-error-rate` budget (default 0.2) that terminal hits do, so a marginal end match still trims but only a tight interior match splits a read. Adapter trims flow through the same tag-rewrite machinery as every other trim, so `MM`/`ML`/`MN` and the per-base tags stay correct.

### Presence detection

A preset catalog (the ONT one has over a hundred sequences) usually holds far more adapters than any single run uses. `--adapter-sample <N>` (N ≥ 100) checks which adapters actually appear in the first N reads and trims the rest against only that set — faster, with no spurious trims from absent catalog entries.

It is off by default (`--adapter-sample 0`) and preset-only: a custom `--adapter-fasta` is always searched in full. If detection finds nothing — an ordered file with clean reads first can look adapter-free — whittle falls back to the full set with a warning rather than skipping trimming for the rest of the run.

### Ab-initio inference

`--adapter-infer` discovers adapters de novo from a sampled read prefix (Porechop_ABI-style k-mer assembly) and trims with the discovered set instead of a catalog. `--adapter-infer-only` discovers and prints them — sequence, support, and any catalog/FASTA cross-name — as FASTA, then exits without trimming or touching the output.

### Built-in ONT catalog

`--adapter-preset ont` loads a catalog assembled for whittle from ONT-published sources: dorado's `adapter_primer_kits.cpp`, Porechop's `adapters.py`, and qcat's kit definitions. It covers 124 sequences — ligation adapters (kit-14 and legacy), the rapid and direct-RNA adapters, PCR/cDNA and 10X primers, barcode flanks, and all 96 barcodes. Reverse-complement search handles both orientations, and fragments under 11 bp are never searched standalone (a pattern that short matches almost anywhere).

## Logging

The log level is `-v`/`-vv` (debug/trace) or `--quiet` (warnings and errors only), overridable with `WHITTLE_LOG` — a `RUST_LOG`-style filter such as `WHITTLE_LOG=whittle::workflow=trace` — where `--quiet` always wins. All logging goes to stderr, so stdout carries only read data. Progress renders as a live bar when stderr is a terminal, or as periodic lines (about every 30s) when it's redirected to a file or pipe.

## Limitations

- **Unaligned BAM only.** Aligned records are refused (whittle names the offending read); there's no CIGAR/POS adjustment for mapped reads.
- **No FASTQ→BAM.** There's no header to build a BAM record from a bare FASTQ read. The reverse, BAM→FASTQ, is supported.
- **No contamination screen.** There's no minimap2-based host filter; adapter discovery, however, is covered by `--adapter-infer`.
- **`--min-length` is post-trim**, applied per output segment rather than to the whole raw read.
- **One quality-trim strategy at a time.** `--qual-trim`, `--qual-best-segment`, and `--qual-split` are mutually exclusive; `-H`/`-T` compose with any of them.

## License

[Apache-2.0](LICENSE). Copyright 2026 Erdi Kılıç.
