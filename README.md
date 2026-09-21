<div align="center">

# whittle

**Coordinate-consistent trimming of long-read uBAM and tagged FASTQ.**

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
- **Adapter and primer trimming.** Terminal adapters, adapters truncated by the read end, and interior adapters (chimera splitting with junction cleanup). Sequences come from built-in kit presets (ONT kit 14 ligation, rapid, barcoding, cDNA, and amplicon kits; RNA004; PacBio SMRTbell), a user FASTA with IUPAC codes, or de novo discovery. Presence detection restricts a preset to the sequences a library carries.
- **Quality trimming.** End trimming to a threshold, best-segment extraction, or splitting at consecutive low-quality bases. Every segment is filtered on its own.
- **Tagged FASTQ.** FASTQ whose headers carry SAM aux tags (`samtools fastq -T MM,ML,MN`) is trimmed with the same tag rewriting as uBAM: `MM`/`ML`/`MN` are rebuilt and per-base arrays sliced for every output segment.
- **Formats.** FASTQ, gzip and BGZF FASTQ, and unaligned BAM as input; the same, plus BAM-to-FASTQ, as output. Formats are detected from the path or the stream, including on stdin. A directory of files is merged in one run.
- **Pipeline integration.** `--summary-json` writes the resolved settings and all counters as JSON. `--preserve-order` keeps the input order under multithreading. Malformed tags are counted and reported rather than fatal.
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
whittle -i reads.fastq.gz -o trimmed.fastq.gz --trim-front 20 --trim-tail 20 \
  --trim-quality 8 --min-length 500 --min-quality 10 -t 8
```

Trim unaligned BAM, split at consecutive low-quality bases, and rewrite the modification tags of every output read.

```bash
whittle -i reads.bam -o trimmed.bam --trim-front 10 --trim-tail 10 \
  --split-quality 9 --split-min-low-quality-bases 50 --min-length 1000
```

Trim FASTQ exported with its tags; the tags are rewritten the same way.

```bash
samtools fastq -T MM,ML,MN reads.bam | whittle -o trimmed.fastq.gz -H 10 -T 10 --trim-quality 10
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
whittle -i reads.fastq.gz --discover-adapters report
whittle -i reads.fastq.gz -o trimmed.fastq.gz --discover-adapters
```

Discovery proceeds in layers from each read end: adapter, barcode flanks,
barcodes and primers, each accepted with its own read support and insert
boundary. Barcode panels are clustered from the reads between their flanks. A preset or FASTA given with `--discover-adapters` is trimmed first
and discovery continues beyond it, which recovers unknown primers behind a
known kit. It can recover multiple recurrent adapters and degenerate primers
from the same sample. Exact k-mer assembly and batched sassy alignments determine the sequences;
catalog matches provide names only. Minority families require independent read
support and enrichment near read ends; short repeats are excluded. Conserved
amplicon sequence is excluded when reads without a primer establish the insert
boundary. Candidates with unresolved boundaries are skipped. For known primers,
an explicit FASTA or kit preset remains the most direct choice; see
[adapter discovery](docs/adapters.md#adapter-discovery) for sampling
requirements and limitations.

Merge a directory, convert BAM to FASTQ, and write a machine-readable summary.

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.fastq.gz --trim-quality 10
whittle -i reads.bam -o reads.fastq.gz -l 500 --quiet --summary-json qc.json
```

`whittle --help` lists every option; [docs/cli.md](docs/cli.md) describes them.

## Quality filtering and trimming

`--min-quality` and `--max-quality` filter each output segment using the
calculation selected by `--quality-mode`: `mean` (average error probability
converted to Phred, the default), `arithmetic` (average Phred score), or
`median`. This setting does not affect trimming or splitting.

| Operation | Parameters | Behavior |
|---|---|---|
| Fixed crop | `--trim-front BASES`, `--trim-tail BASES` | Remove a fixed number of bases from each adapter-derived segment's 5' and 3' ends after barcode restriction |
| Quality end trimming | `--trim-quality PHRED` | Trim each end until reaching a base at or above the threshold |
| Best segment | `--best-quality-segment PHRED` | Select the highest-scoring contiguous segment using cumulative error probabilities; bases below the threshold can be retained |
| Quality splitting | `--split-quality PHRED`, `--split-min-low-quality-bases BASES` | Split at the specified number of consecutive bases below the threshold; retain shorter internal stretches |

The three quality operations are mutually exclusive and apply separately to
each segment produced by adapter processing. Filters apply after these
operations. `--min-length` sets the minimum retained segment length;
`--split-min-low-quality-bases` sets the number of low-quality bases required
to split.

`--head-crop` and `--tail-crop` are aliases for `--trim-front` and `--trim-tail`.
Both accept a base count and retain the short options `-H` and `-T`.

## Trimming pipeline

Adapter preparation loads a FASTA or preset, or discovers adapters from sampled
original reads. Sampled reads remain in the processing stream. Processing then
follows this order:

1. **Adapters.** Search the original read, trim terminal adapters and primers, and split at interior adapters. Clean each new end. Reads without a match continue as one segment.
2. **Barcodes.** Barcode spans recorded in the `bi` tag are removed where a barcode sequence is found at them.
3. **Fixed crop.** `--trim-front` and `--trim-tail` crop each retained adapter-derived segment once.
4. **Quality.** Apply `--trim-quality`, `--best-quality-segment`, or `--split-quality` to each cropped segment. Quality splitting can produce further segments; these are not cropped again.
5. **Filter.** Each final segment must pass the length, quality, and GC bounds.
6. **Output.** Rewrite tags against each surviving interval and write the records.

Final segments are numbered once, in their order along the original read.
Existing PacBio query intervals in read names are updated on crops as well as
splits. Tagged FASTQ headers are detected per record, including after plain
records in a merged directory.

Two adapter segments that each split into two quality segments produce
`<read>_segment_1` through `<read>_segment_4`. Filtering preserves the indices
of surviving segments. PacBio records use final query-coordinate names.
Every stage retains original-read coordinates, and tags are rewritten once
against each final interval ([docs/tags.md](docs/tags.md)).

## Formats

| input \ output | FASTQ | FASTQ.gz | FASTQ.bgz | BAM |
|---|:---:|:---:|:---:|:---:|
| FASTQ, FASTQ.gz, FASTQ.bgz | yes | yes | yes | no |
| unaligned BAM | yes | yes | yes | yes |

Formats are taken from the path extension, a stream sniff, or `--input-format`/`--output-format`. FASTQ-to-BAM is not supported, since a FASTQ read carries no header from which to build a BAM record.

## Documentation

| Page | Contents |
|---|---|
| [docs/cli.md](docs/cli.md) | Options, format selection, directory input, summary JSON, logging and progress |
| [docs/adapters.md](docs/adapters.md) | Adapter roles, partial adapters, chimera splitting, presence detection, de novo inference, kit catalog |
| [docs/tags.md](docs/tags.md) | Handling of `MM`/`ML`/`MN`, kinetics, and ONT signal tags |
| [CHANGELOG.md](CHANGELOG.md) | Release history |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Building, testing, commit conventions, style |

## Limitations

- **Unaligned BAM only.** Aligned records are refused; there is no CIGAR or POS adjustment for mapped reads.
- **No FASTQ-to-BAM.** BAM-to-FASTQ is supported.
- **`--min-length` applies after trimming**, per output segment, not to the raw read.
- **BAM folder output requires matching read groups.** Conflicting definitions are rejected.
- **Signal rewriting requires model direction.** `--update-moves` uses a DNA or RNA `basecall_model` in the BAM read-group description.
- **One quality-trim strategy per run.** `--trim-quality`, `--best-quality-segment`, and `--split-quality` are mutually exclusive; `-H`/`-T` combine with any of them.

## Citation

Kılıç E. whittle: coordinate-consistent trimming of long-read uBAM and tagged FASTQ. Zenodo. https://doi.org/10.5281/zenodo.21355499

## License

[Apache-2.0](LICENSE). Copyright 2026 Erdi Kılıç.
