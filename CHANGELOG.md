# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
