<div align="center">

# whittle

**Tag-aware trimming of long-read FASTQ and unaligned BAM.**

[![CI](https://github.com/erdikilic/whittle/actions/workflows/ci.yml/badge.svg)](https://github.com/erdikilic/whittle/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/whittle.svg)](https://crates.io/crates/whittle)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/whittle.svg)](https://anaconda.org/bioconda/whittle)
[![License: Apache 2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21355499.svg)](https://doi.org/10.5281/zenodo.21355499)
![Rust](https://img.shields.io/badge/rust-2024%20edition-000000?logo=rust&logoColor=white)
![Long reads](https://img.shields.io/badge/long--reads-ONT%20%7C%20PacBio-1f6feb)

</div>

whittle filters, trims, and splits Oxford Nanopore and PacBio reads in FASTQ, compressed FASTQ, and unaligned BAM (uBAM). Unaligned BAM carries per-read tags that are indexed by base position: base-modification calls (`MM`, `ML`, `MN`), per-base kinetics, and the ONT move table (`mv`, `ns`, `ts`). whittle rewrites these tags on every trim and split, so each output read stays consistent with its sequence.

## Features

- **Base-modification tags.** `MM`, `ML`, and `MN` are reconstructed for every trimmed or split uBAM read. The test suite checks the result against an independent `htslib` decoder.
- **Kinetics and signal tags.** Per-base arrays (`ip`, `pw`, and related) are sliced with the sequence. ONT signal tags (`mv`, `ts`, `ns`, and related) are removed by default or rewritten with `--update-moves`.
- **Adapter and primer trimming.** Terminal adapters, adapters truncated by the read end, and interior adapters (chimera splitting with junction cleanup). Sequences come from built-in kit presets (ONT kit 14 ligation, rapid, barcoding, cDNA, and amplicon kits; RNA004; PacBio SMRTbell), a user FASTA with IUPAC codes, or ab-initio discovery. Presence detection restricts a preset to the sequences a library carries.
- **Quality trimming.** End trimming to a threshold, best-segment extraction, or splitting at low-quality runs. Every segment is filtered on its own.
- **Formats.** FASTQ, gzip and BGZF FASTQ, and unaligned BAM as input; the same, plus BAM-to-FASTQ, as output. Formats are detected from the path or the stream, including on stdin. A directory of files is merged in one run.
- **Pipeline integration.** `--summary-json` writes the resolved settings and all counters as JSON. `--ordered` keeps the input order under multithreading. Malformed tags are counted and reported rather than fatal.
- **Performance.** Multithreaded trimming, encoding, and decoding; `-t N` uses about N cores. Adapter search is SIMD bit-parallel ([sassy](https://github.com/RagnarGrootKoerkamp/sassy)). BAM and BGZF I/O use `noodles` with `libdeflate`. No `htslib` is required at build or run time.

## Installation

Bioconda:

```bash
conda install -c conda-forge -c bioconda whittle
```

Prebuilt binaries for Linux and macOS (x86-64 and arm64, glibc and static musl) are on the [Releases](https://github.com/erdikilic/whittle/releases) page. Each tarball includes the man page, which is also in [`man/`](man).

crates.io. The adapter search requires AVX2 on x86-64, and `cargo install` does not read the repository's `.cargo/config.toml`, so the target is passed explicitly:

```bash
RUSTFLAGS="-C target-cpu=x86-64-v3" cargo install whittle
```

From source (Rust 1.91 or newer):

```bash
git clone https://github.com/erdikilic/whittle
cd whittle
cargo build --release   # target/release/whittle
```

## Usage

Filter and trim FASTQ: crop 20 bp from each end, quality-trim below Q8, keep reads of at least 500 bp and Q10.

```bash
whittle -i reads.fastq.gz -o trimmed.fastq.gz -H 20 -T 20 --qual-trim 8 -l 500 -q 10 -t 8
```

Trim unaligned BAM, split at low-quality runs, and rewrite the modification tags of every output read.

```bash
whittle -i reads.bam -o trimmed.bam -H 10 -T 10 -l 1000 --qual-split 9 --qual-split-window 50
```

Trim adapters with a kit preset. Interior adapters split the read.

```bash
whittle -i reads.bam -o trimmed.bam --adapter-preset lsk114 -l 500
whittle -i reads.fastq.gz -o trimmed.fastq.gz --adapter-preset nbd114 -t 16
```

Trim amplicon primers with the `mab114` preset (degenerate 16S 27F/1492R and ITS1F/ITS4 primers) or a custom FASTA. IUPAC codes are accepted; `primer` in a header restricts the entry to the read ends.

```bash
whittle -i 16s.fastq.gz -o trimmed.fastq.gz --adapter-preset mab114
whittle -i 16s.fastq.gz -o trimmed.fastq.gz -a primers.fa
```

```
>27F primer
AGAGTTTGATYMTGGCTCAG
>1492R primer
TACGGYTACCTTGTTACGACTT
```

Trim PacBio SMRTbell adapters and split concatemers.

```bash
whittle -i hifi.bam -o trimmed.bam --adapter-preset pacbio
```

Discover adapters de novo: report the candidates, or trim with them directly.

```bash
whittle -i reads.fastq.gz --adapter-infer report
whittle -i reads.fastq.gz -o trimmed.fastq.gz --adapter-infer
```

Merge a directory, convert BAM to FASTQ, and write a machine-readable summary.

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.fastq.gz --qual-trim 10
whittle -i reads.bam -o reads.fastq.gz -l 500 --quiet --summary-json qc.json
```

`whittle --help` lists every option; [docs/cli.md](docs/cli.md) describes them.

## Trimming pipeline

The stages run in a fixed order, and each stage operates on what the previous one left:

1. **Barcodes.** `--trim-barcodes` removes the spans recorded in the `bi` tag.
2. **Fixed crop.** `-H`/`--head-crop` and `-T`/`--tail-crop`.
3. **Adapters.** Terminal adapters and primers are trimmed. An interior adapter splits the read, the junction is excised, and each side is re-trimmed at its new end.
4. **Quality.** One of `--qual-trim`, `--qual-best-segment`, or `--qual-split`.
5. **Filter.** Each surviving segment must pass `-l`/`-L` (length), `-q`/`-Q` (quality), and `-g`/`-G` (GC).

A split read yields segments named `<read>_segment_N`, each filtered independently, so `-l` is a per-segment minimum after trimming. Every stage is expressed as an interval on the original read, and the tags are rewritten once against the final interval ([docs/tags.md](docs/tags.md)).

## Formats

| input \ output | FASTQ | FASTQ.gz | FASTQ.bgz | BAM |
|---|:---:|:---:|:---:|:---:|
| FASTQ, FASTQ.gz, FASTQ.bgz | yes | yes | yes | no |
| unaligned BAM | yes | yes | yes | yes |

Formats are taken from the path extension, a stream sniff, or `--in-format`/`--out-format`. FASTQ-to-BAM is not supported, since a FASTQ read carries no header from which to build a BAM record.

## Documentation

| Page | Contents |
|---|---|
| [docs/cli.md](docs/cli.md) | Options, format selection, directory input, summary JSON, logging and progress |
| [docs/adapters.md](docs/adapters.md) | Adapter roles, partial adapters, chimera splitting, presence detection, ab-initio inference, kit catalog |
| [docs/tags.md](docs/tags.md) | Handling of `MM`/`ML`/`MN`, kinetics, and ONT signal tags |
| [CHANGELOG.md](CHANGELOG.md) | Release history |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Building, testing, commit conventions, style |

## Limitations

- **Unaligned BAM only.** Aligned records are refused; there is no CIGAR or POS adjustment for mapped reads.
- **No FASTQ-to-BAM.** BAM-to-FASTQ is supported.
- **`--min-length` applies after trimming**, per output segment, not to the raw read.
- **One quality-trim strategy per run.** `--qual-trim`, `--qual-best-segment`, and `--qual-split` are mutually exclusive; `-H`/`-T` combine with any of them.

## Citation

Kılıç E. whittle: tag-aware trimming of long-read FASTQ and unaligned BAM. Zenodo. https://doi.org/10.5281/zenodo.21355499

## License

[Apache-2.0](LICENSE). Copyright 2026 Erdi Kılıç.
