<div align="center">

# whittle

**A fast, tag-aware trimmer for long-read FASTQ and unaligned BAM.**

It rewrites the position-indexed tags on every trim and split (`MM`/`ML` modification calls, per-base kinetics, and ONT signal), so a trimmed read stays in register with its sequence.

[![CI](https://github.com/erdikilic/whittle/actions/workflows/ci.yml/badge.svg)](https://github.com/erdikilic/whittle/actions/workflows/ci.yml)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
![Rust](https://img.shields.io/badge/rust-2024%20edition-000000?logo=rust&logoColor=white)
![Long reads](https://img.shields.io/badge/long--reads-ONT%20%7C%20PacBio-1f6feb)

</div>

whittle filters and trims long reads (ONT, PacBio) in FASTQ, gzip/BGZF-compressed FASTQ, and unaligned BAM. It handles the usual length/quality/GC filtering and head/tail/quality/adapter trimming. What most trimmers get wrong on uBAM is the tags: crop the sequence and the base-modification, kinetics, and signal tags now point at bases that are gone. whittle keeps them in sync, so every output read is still valid.

## Highlights

- **Correct modification tags.** `MM`/`ML`/`MN` are rebuilt for every trimmed or split uBAM read, and checked against an independent `htslib` decoder.
- **Trim-aware tags.** Per-base kinetics (`ip`/`pw`/…) are sliced along with the sequence. ONT signal tags (`mv`/`ts`/`ns`/…) are dropped, or rewritten dorado-style with `--update-moves`.
- **Adapter trimming.** Terminal trimming plus interior chimera splitting, driven by a built-in ONT catalog, your own FASTA, or ab-initio discovery.
- **Formats.** FASTQ, gzip/BGZF-compressed FASTQ, and unaligned BAM, plus BAM→FASTQ conversion. Formats are auto-detected, including BGZF FASTQ or BAM piped over stdin.
- **Fast and self-contained.** Multithreaded throughout, with a thread budget that adapts to the workload, and no external `htslib` to build or run.

## Install

### Prebuilt binaries

Download a binary for your platform from the [Releases](https://github.com/erdikilic/whittle/releases) page and put it on your `PATH`. Builds cover Linux and macOS, x86-64 and arm64, glibc and static musl.

### From source

```bash
git clone https://github.com/erdikilic/whittle
cd whittle
cargo build --release   # -> target/release/whittle
```

### From crates.io

The adapter search ([`sassy`](https://crates.io/crates/sassy)) needs AVX2 on x86-64. A `cargo install` won't inherit this repo's build config, so pass the flag yourself:

```bash
RUSTFLAGS="-C target-cpu=x86-64-v3" cargo install whittle
```

Building needs Rust 1.91 or newer. There's no external `htslib` dependency: BAM I/O goes through the `libdeflate` backend of `noodles-bgzf`, and the optional `rust-htslib` is a dev-only test dependency.

## Quick start

Trim a FASTQ file. Crop 20 bp off each end, quality-trim below Q8, drop reads under 500 bp or Q10, run on 8 threads:

```bash
whittle -i reads.fastq.gz -o trimmed.fastq.gz -H 20 -T 20 --qual-trim 8 -l 500 -q 10 -t 8
```

Trim an unaligned BAM, split at low-quality runs, and rebuild the modification tags on every output read:

```bash
whittle -i reads.bam -o trimmed.bam -H 10 -T 10 -l 1000 --qual-split 9 --qual-split-window 50
```

## Usage

whittle reads from `-i`/`--input` (or stdin) and writes to `-o`/`--output` (or stdout). It takes the format from the file extension, sniffs it from the first bytes of a stream, or accepts it from `--in-format`/`--out-format {fastq,fastq-gz,fastq-bgz,bam}`.

Output is plain FASTQ by default and is never compressed on its own. A `.gz` or `.bgz` input does not imply compressed output; you get that only by asking for it, with a `.gz`/`.bgz` output path or the matching format flag. Compressed output is written by a parallel encoder using `-t` threads. BGZF FASTQ input is decompressed block-parallel too; ordinary gzip stays a serial input format.

### How trimming works

The operations run in a fixed order, and the filters apply to whatever survives:

1. **Fixed crop.** `-H`/`--head-crop` and `-T`/`--tail-crop` remove a set number of bases from each end.
2. **Adapters.** Terminal adapters are trimmed, and interior adapters split the read (see [Adapter trimming](#adapter-trimming)).
3. **Quality.** One of `--qual-trim`, `--qual-best-segment`, or `--qual-split` (mutually exclusive).
4. **Filter.** Each surviving segment must pass `-l`/`-L` (length), `-q`/`-Q` (quality), and `-g`/`-G` (GC).

When a read splits, each segment is filtered on its own and named `<read>_segment_N` (1-based), so `-l` is a post-trim, per-segment minimum.

### Formats

| input → | FASTQ | FASTQ.gz | FASTQ.bgz | BAM |
|---|:---:|:---:|:---:|:---:|
| FASTQ / FASTQ.gz / FASTQ.bgz | ✅ | ✅ | ✅ | ❌ |
| unaligned BAM | ✅ | ✅ | ✅ | ✅ |

With no output extension or `--out-format`, the output format mirrors the input, except that compressed FASTQ input defaults back to plain FASTQ. FASTQ→BAM isn't supported; there's no header to build a BAM record from. BGZF streams are recognized by their decompressed payload, so piped FASTQ.bgz and `samtools view -b … | whittle` need no hint.

On BAM→FASTQ, aux tags go into the FASTQ header tab-delimited, following the `samtools fastq -T` convention. `--fastq-tags` picks which ones: `all` (default), `none`, or a list like `MM,ML,RG`. `MM`/`ML`/`MN` are reconstructed for the trimmed segment, per-base tags are sliced, and everything else is copied verbatim.

### Folder input

`-i` also takes a directory. whittle merges every read file directly inside it, in sorted filename order, into one output. The folder has to be a single format (all FASTQ-family or all BAM); subdirectories are ignored, and a mixed or empty folder is an error.

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --qual-trim 10
```

### Options

| Flag | Meaning |
|---|---|
| `--version` | Print the version and exit |
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
| `-g, --min-gc <F>`, `-G, --max-gc <F>` | GC-fraction bounds (0 to 1) |
| `-m, --qual-mode {mean,arithmetic,median}` | How read quality is summarized (default `mean`, the error-probability mean) |
| `-H, --head-crop <N>`, `-T, --tail-crop <N>` | Fixed crop from each end; always runs first |
| `--qual-trim <Q>` | Trim low-quality bases from both ends down to the first base ≥ Q |
| `--qual-best-segment <Q>` | Keep only the longest contiguous run of quality ≥ Q |
| `--qual-split <Q>` | Split at low-quality (< Q) runs, keeping each surviving segment |
| `--qual-split-window <N>` | Tolerate low-quality runs shorter than N without splitting (default 1) |
| `--update-moves` | Rewrite ONT signal tags through trimming instead of dropping them (BAM→BAM) |
| `-a, --adapter-fasta <FILE>` | Adapter/primer FASTA; enables adapter trimming |
| `--adapter-preset {none,ont}` | Built-in adapter catalog (default `none`; `ont` enables trimming) |
| `--adapter-error-rate <F>` | End-match tolerance as a fraction of adapter length (default 0.2) |
| `--adapter-end-size <N>` | End-zone width searched for terminal adapters (default 150) |
| `--adapter-ends-only` | Trim ends only; never split on an interior adapter |
| `--adapter-sample <N>` | Reads sampled for preset detection or inference (defaults `0` and `40000`, respectively) |
| `--adapter-infer [trim\|report]` | Discover adapters de novo; omitted value defaults to `trim` |
| `--adapter-infer-policy {conservative,aggressive}` | Trust policy for inferred adapters (default `conservative`) |
| `-v`, `-vv` | Increase log detail (debug, trace); higher counts are rejected |
| `--quiet` | Silence progress and the summary; warnings and errors still print |

`--qual-trim`, `--qual-best-segment`, and `--qual-split` are three strategies for the same step, so pass at most one. `-H`/`-T` are independent and compose with whichever you pick.

## Base-modification tags

Long-read modification calls (5mC, 6mA, and so on) live per-read in two tags: `MM`, which says which bases are modified as skip-counts over the sequence, and `ML`, which holds their probabilities. Trim `SEQ`/`QUAL` without touching these and the result decodes to nonsense, because the skip-counts still index the original sequence rather than the trimmed one.

whittle rebuilds them as part of trimming. For every output uBAM read, whether cropped, quality-trimmed, or split, `MM` and `ML` are reconstructed against the output window: skip-counts renumbered, probability bytes re-sliced, and `MN` updated to match. If trimming removed every modified base, all three are dropped. Everything else in the record rides through unchanged.

This is covered by decode-equivalence tests. They re-decode whittle's output with `rust-htslib`'s `basemods_iter()`, a different `MM`/`ML` implementation from the one whittle writes with, and compare against the original calls restricted to the surviving window. One test always runs on a synthetic fixture; another sweeps a real uBAM when you point it at one:

```bash
WHITTLE_UBAM=/path/to/real.ubam cargo test --test bam_mods_oracle -- --ignored
```

### Other trim-aware tags

Every tag indexed by base position is kept consistent, on both BAM→BAM and BAM→FASTQ:

| Tag(s) | On a trimmed read |
|---|---|
| `MM` / `ML` / `MN` | Reconstructed for the output window |
| Per-base kinetics (`ip`/`pw`/`fi`/`fp`/`ri`/`rp`), and any read-length `B` array | Sliced in lockstep with the sequence |
| ONT signal (`mv`/`ts`/`ns`/`sp`/`pi`) | Dropped, or rewritten with `--update-moves` |
| Poly-A (`pa`/`pt`) | Kept/shifted with `--update-moves` if the tail survives, else dropped |
| `bi` (barcode positions) | Dropped, since the positions shift under a crop |
| `qs` (mean qscore) | Recomputed from the trimmed quality |
| `st`/`du` (start time / duration) | Kept on a crop, dropped on a split |
| `RG`, `ch`, `mx`, `sm`/`sd`/`sv`, … | Copied verbatim |

With `--update-moves`, a crop slices `mv` and advances `ts`/`ns`, while a split emits dorado-style subreads (`pi` parent id, `sp` parent-signal offset, `ns` subread span, `ts` 0, `rn` -1) so the renamed segment stays locatable in POD5 for tools like Remora. BAM→FASTQ always drops the signal tags on a trim, since a move table in a FASTQ header is impractical. If a known per-base tag's length doesn't match the sequence (malformed input), whittle leaves it untouched and prints a one-line advisory.

## Adapter trimming

Off by default. Turn it on with `-a`/`--adapter-fasta <FILE>` (your own sequences, one per record, each at least 11 bp) and/or `--adapter-preset ont` (the built-in catalog). Either one alone is enough, and they combine. Every adapter is searched on both strands, so orientation doesn't matter, and each read gets two treatments:

- **Terminal trimming.** An adapter within `--adapter-end-size` bases of an end (default 150) is trimmed off.
- **Chimera splitting.** An adapter in the interior is treated as a junction. The read splits there, the adapter is excised, and both sides are kept. `--adapter-ends-only` turns this off and searches only the two end-zones.

Interior hits use half the `--adapter-error-rate` budget (default 0.2) that terminal hits do, so a marginal end match still trims but only a tight interior match splits a read. Adapter trims flow through the same tag-rewrite machinery as every other trim, so `MM`/`ML`/`MN` and the per-base tags stay correct.

### Presence detection

A preset catalog holds far more adapters than any single run uses (the ONT one has over a hundred). `--adapter-sample <N>` (N ≥ 100) checks which adapters actually turn up in the first N reads, then trims the rest against only that set. It's faster, and it avoids spurious trims from catalog entries that aren't present.

Detection is off by default (`--adapter-sample 0`) and preset-only; a custom `--adapter-fasta` is always searched in full. If detection finds nothing (an ordered file with clean reads first can look adapter-free), whittle warns and falls back to the full set rather than skipping trimming for the rest of the run.

### Ab-initio inference

`--adapter-infer` (the same as `--adapter-infer trim`) discovers recurrent read-end sequences de novo from a sampled read prefix, using Porechop_ABI-style k-mer assembly, then trims with what it found. By default it trims ends only, with a conservative anchor of at most 32 bp facing the physical end. Anything longer on the insert-facing side of the assembled consensus is reported as uncertain rather than assumed to be technical. This matters for amplicons: without a known primer or reference, a primer and a conserved marker-gene prefix can be statistically indistinguishable.

`--adapter-infer report` prints the recommended anchor, its support, the assembled length, the uncertain-base count, and any catalog/FASTA cross-name, all as FASTA, then exits without touching record output. Add `-v` to log the full review-only consensus. `--adapter-infer-policy aggressive` restores full-consensus trimming and allows interior splitting unless `--adapter-ends-only` is also set; reach for it only once you've ruled out overtrimming conserved biological sequence. The default policy is `conservative`.

### Built-in ONT catalog

`--adapter-preset ont` loads a catalog assembled for whittle from ONT-published sources: dorado's `adapter_primer_kits.cpp`, Porechop's `adapters.py`, and qcat's kit definitions. It has 124 sequences: ligation adapters (kit-14 and legacy), the rapid and direct-RNA adapters, PCR/cDNA and 10X primers, barcode flanks, and all 96 barcodes. Reverse-complement search covers both orientations, and fragments under 11 bp are never searched on their own, since a pattern that short matches almost anywhere.

## Logging

Set the log level with `-v`/`-vv` (debug/trace) or `--quiet` (warnings and errors only). `WHITTLE_LOG` overrides it with a `RUST_LOG`-style filter, for example `WHITTLE_LOG=whittle::workflow=trace`, and `--quiet` still wins over it. All logging goes to stderr, so stdout carries only read data. Progress shows as a live bar when stderr is a terminal, or as periodic lines (about every 30s) when it's redirected to a file or pipe.

## Limitations

- **Unaligned BAM only.** Aligned records are refused, with whittle naming the offending read; there's no CIGAR/POS adjustment for mapped reads.
- **No FASTQ→BAM.** There's no header to build a BAM record from a bare FASTQ read. The reverse, BAM→FASTQ, works.
- **`--min-length` is post-trim**, applied per output segment rather than to the whole raw read.
- **One quality-trim strategy at a time.** `--qual-trim`, `--qual-best-segment`, and `--qual-split` are mutually exclusive; `-H`/`-T` compose with any of them.

## Development

### Commit messages

Commit subjects follow [Conventional Commits](https://www.conventionalcommits.org/), for example `feat(adapter): add a trimming mode`, `fix: handle empty input`, or `perf: improve processing throughput`. Enable the repository's versioned validation hook once after cloning:

```bash
git config core.hooksPath .githooks
```

The hook rejects invalid subjects before a commit is created. Git won't enable repository-provided hooks on its own, for security reasons, so CI validates every pushed or pull-request commit as the repo-wide backstop. You can bypass the local hook with `--no-verify`; to make it non-optional, configure branch protection to require the `commit-message` CI job before merging.

Give manual merge commits a compliant subject instead of Git's default `Merge ...` message:

```bash
git merge --no-ff feature-branch -m "perf: integrate throughput improvements"
```

## License

[Apache-2.0](LICENSE). Copyright 2026 Erdi Kılıç.
