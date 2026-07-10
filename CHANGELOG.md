# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/erdikilic/whittle/commits/main
