<div align="center">

# whittle

**A fast, tag-aware trimmer for long-read FASTQ and unaligned BAM.**

It rewrites the position-indexed tags on every trim and split (`MM`/`ML` modification calls, per-base kinetics, and ONT signal), so a trimmed read stays in register with its sequence.

[![CI](https://github.com/erdikilic/whittle/actions/workflows/ci.yml/badge.svg)](https://github.com/erdikilic/whittle/actions/workflows/ci.yml)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21355499.svg)](https://doi.org/10.5281/zenodo.21355499)
![Rust](https://img.shields.io/badge/rust-2024%20edition-000000?logo=rust&logoColor=white)
![Long reads](https://img.shields.io/badge/long--reads-ONT%20%7C%20PacBio-1f6feb)

</div>

whittle filters and trims long reads (ONT, PacBio) in FASTQ, gzip/BGZF-compressed FASTQ, and unaligned BAM. It handles the usual length/quality/GC filtering and head/tail/quality/adapter trimming. What most trimmers get wrong on uBAM is the tags: crop the sequence and the base-modification, kinetics, and signal tags now point at bases that are gone. whittle keeps them in sync, so every output read is still valid.

## Highlights

- **Correct modification tags.** `MM`/`ML`/`MN` are rebuilt for every trimmed or split uBAM read, and checked against an independent `htslib` decoder.
- **Trim-aware tags.** Per-base kinetics (`ip`/`pw`/…) are sliced along with the sequence. ONT signal tags (`mv`/`ts`/`ns`/…) are dropped, or rewritten dorado-style with `--update-moves`.
- **Adapter trimming.** Terminal trimming plus interior chimera splitting, driven by a built-in ONT catalog, your own FASTA (IUPAC codes included, so a degenerate primer matches every variant it covers), or ab-initio discovery.
- **Formats.** FASTQ, gzip/BGZF-compressed FASTQ, and unaligned BAM, plus BAM→FASTQ conversion. Formats are auto-detected, including BGZF FASTQ or BAM piped over stdin.
- **Pipeline-friendly.** `--summary-json` writes the run's counters and resolved settings as JSON, even under `--quiet`.
- **Fast and self-contained.** Multithreaded throughout, with a thread budget that adapts to the workload, and no external `htslib` to build or run.

## Install

**Prebuilt binaries.** Download one for your platform from the [Releases](https://github.com/erdikilic/whittle/releases) page and put it on your `PATH`. Builds cover Linux and macOS, x86-64 and arm64, glibc and static musl. Each tarball also carries the man page, which lives in [`man/`](man) here.

**From source.**

```bash
git clone https://github.com/erdikilic/whittle
cd whittle
cargo build --release   # -> target/release/whittle
```

**From crates.io.** The adapter search ([`sassy`](https://crates.io/crates/sassy)) needs AVX2 on x86-64. A `cargo install` won't inherit this repo's build config, so pass the flag yourself:

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

Merge and trim a whole basecalling folder:

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --qual-trim 10
```

## How trimming works

The operations run in a fixed order, and the filters apply to whatever survives:

1. **Fixed crop.** `-H`/`--head-crop` and `-T`/`--tail-crop` remove a set number of bases from each end.
2. **Adapters.** Terminal adapters are trimmed, and interior adapters split the read.
3. **Quality.** One of `--qual-trim`, `--qual-best-segment`, or `--qual-split` (mutually exclusive).
4. **Filter.** Each surviving segment must pass `-l`/`-L` (length), `-q`/`-Q` (quality), and `-g`/`-G` (GC).

When a read splits, each segment is filtered on its own and named `<read>_segment_N` (1-based), so `-l` is a post-trim, per-segment minimum.

## Formats

| input → | FASTQ | FASTQ.gz | FASTQ.bgz | BAM |
|---|:---:|:---:|:---:|:---:|
| FASTQ / FASTQ.gz / FASTQ.bgz | ✅ | ✅ | ✅ | ❌ |
| unaligned BAM | ✅ | ✅ | ✅ | ✅ |

Formats come from the path extension, a stream sniff, or `--in-format`/`--out-format`. Output is never compressed unless you ask for it, and FASTQ→BAM isn't supported because there's no header to build a BAM record from.

## Documentation

| Page | What's in it |
|---|---|
| [docs/cli.md](docs/cli.md) | Every flag, format selection, folder input, logging and progress |
| [docs/tags.md](docs/tags.md) | How `MM`/`ML`/`MN`, kinetics, and ONT signal tags are kept in register |
| [docs/adapters.md](docs/adapters.md) | Adapter trimming, presence detection, ab-initio inference, the ONT catalog |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Building, testing, commit conventions, style |

## Limitations

- **Unaligned BAM only.** Aligned records are refused, with whittle naming the offending read; there's no CIGAR/POS adjustment for mapped reads.
- **No FASTQ→BAM.** There's no header to build a BAM record from a bare FASTQ read. The reverse, BAM→FASTQ, works.
- **`--min-length` is post-trim**, applied per output segment rather than to the whole raw read.
- **One quality-trim strategy at a time.** `--qual-trim`, `--qual-best-segment`, and `--qual-split` are mutually exclusive; `-H`/`-T` compose with any of them.

## License

[Apache-2.0](LICENSE). Copyright 2026 Erdi Kılıç.
