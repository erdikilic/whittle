# whittle

[![License: Apache 2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
![Rust](https://img.shields.io/badge/rust-2024%20edition-000000?logo=rust&logoColor=white)
![Long reads](https://img.shields.io/badge/long--reads-ONT%20%7C%20PacBio-1f6feb)

A tag-aware long-read (ONT/PacBio) trimmer. `whittle` filters and trims
FASTQ and unaligned BAM (uBAM) reads the same way tools like `chopper` do —
but when the input is uBAM, it also **recomputes the `MM`/`ML`/`MN`
base-modification tags** so a trimmed or split read's methylation calls stay
correct instead of silently drifting out of register with the sequence.

## Install

```bash
cargo build --release
```

The binary is written to `target/release/whittle`. Building links against
`noodles-bgzf`'s `libdeflate` backend for BAM I/O, so no external `htslib` is
required to build or run (an optional `rust-htslib`-based test dependency is
dev-only, used by the oracle tests described below).

## Usage

`whittle` reads from `-i/--input` (or stdin) and writes to `-o/--output` (or
stdout). Format is detected from the file extension, or sniffed from the
first bytes when reading a stream; it can also be forced with
`--in-format`/`--out-format {fastq,fastq-gz,bam}`.

**Output is plain FASTQ by default — it is never auto-compressed.** A `.gz`
input does not imply gzipped output: with no `-o` extension and no
`--out-format`, `whittle` always writes plain FASTQ (this matters most on
stdout, where silently emitting gzip bytes would be surprising and is also
slower). Gzip output only happens when you ask for it explicitly, either via
an `-o` path ending in `.gz` (e.g. `out.fastq.gz`, `out.fq.gz`, or a bare
`out.gz`) or `--out-format fastq-gz`. When it
does, it's produced by a parallel gzip encoder (`gzp`) that uses `-t/--threads`
worker threads — which defaults to all detected CPUs, so gzip output is
parallel out of the box; `-t 8 -o out.fastq.gz` pins it to 8 threads instead
(values above the CPU count are clamped down to it).

### Format conversion

`whittle` auto-detects formats from file extensions (or `--in-format` /
`--out-format`). Supported conversions:

| input → output | FASTQ | FASTQ.gz | BAM |
|----------------|:-----:|:--------:|:---:|
| FASTQ / FASTQ.gz | ✅ | ✅ | ❌ |
| unaligned BAM   | ✅ | ✅ | ✅ |

When no `-o` extension / `--out-format` is given, output **mirrors** the input
(BAM stays BAM), except a `.gz` input defaults to plain FASTQ (never
auto-compressed). FASTQ→BAM is not supported (there is no header/tags to build a
BAM from).

#### BAM → FASTQ tags (`--fastq-tags`)

On BAM→FASTQ, aux tags are written into the FASTQ header, tab-delimited, in the
`samtools fastq -T` / `samtools import -T` convention
(`@read\tMM:Z:…\tML:B:C,…`). MM/ML/MN are **reconstructed** for the trimmed
segment; per-base tags are **sliced** (see below); remaining tags are copied
**verbatim**.

```text
--fastq-tags all     # default: carry every aux tag
--fastq-tags none    # plain FASTQ, no tags
--fastq-tags MM,ML   # only the (reconstructed) modification tags
--fastq-tags MM,ML,RG
```

### Trim-aware tag handling

Trimming/splitting rewrites tags that are indexed by base position so the output
stays internally consistent (applies to both BAM→BAM and BAM→FASTQ):

- **MM/ML/MN** (base modifications) — reconstructed for the window.
- **Per-base arrays** — sliced in lockstep with the sequence. This covers PacBio
  base kinetics (`ip`, `pw`, `fi`, `fp`, `ri`, `rp`) and, structurally, *any* `B`
  array whose length equals the read length. Without this a trimmed PacBio record
  would have kinetics arrays longer than its sequence — an invalid record.
- **ONT signal tags** (`mv` move table + `ts`/`ns`/`sp`/`pi`) map bases to the raw
  signal. By default they are **dropped** on a trimmed read (a stale move table
  misleads signal-aware tools). Pass **`--update-moves`** to instead keep them
  consistent (BAM→BAM), so a trimmed read stays usable by Remora / Clair3 v2:
  - a head/tail **crop** keeps the read name, slices `mv`, advances `ts` past the
    trimmed-off front signal, and sets `ns = ts + span` (dorado's `ns = trim +
    basecalled-span`, so a head-only crop leaves `ns` unchanged and a tail crop
    shrinks it);
  - a `--qual-split` **split** emits dorado-style subreads (`pi` = parent id,
    `sp` = offset into the parent signal, `ns` = subread span, `ts` = 0, `rn` = -1),
    so the renamed segment's signal is still locatable in POD5.

  The `mv`/`ts`/`ns`/`sp`/`pi`/`rn` values follow dorado's own `splitter::subread`
  and `generate_read_tags` (move blocks sliced by stride-aligned range; `#1`s stays
  equal to the sequence length). BAM→FASTQ always drops these tags on trim (a move
  table in a FASTQ header is impractical).
- **Poly-A tags** (`pa` signal boundaries, `pt` tail length in bases). `pa` holds
  absolute original-signal positions, so under `--update-moves` they're **kept**
  when the tail survives the trim: a crop leaves them valid as-is, a split shifts
  `pa` into the subread's own signal frame. If the trim cuts into the tail (any
  `pa` position falls outside the kept signal window), both are **dropped**.
  Without `--update-moves` (or a malformed move table) they're dropped.
- **`bi`** (barcode info) embeds front/rear sequence positions that shift under a
  crop, so it's **dropped** on trim. The barcode *call* (`BC`/`bv`) is a per-read
  label and is kept.
- **`qs`** (mean read qscore) is **recomputed** from the trimmed quality on a
  trimmed read (matching dorado's per-(sub)read `qs`).
- **`st` (start time) / `du` (duration)** are kept on a crop (same read identity)
  but **dropped on a split** — a subread starts later in the signal, and dorado
  recomputes these from the sample rate, which isn't carried in the BAM.
- **Signal-scaling scalars** (`sm`/`sd`/`sv`) and **per-read metadata** (`RG`,
  `ch`, `mx`, `dx`, `fn`, `BC`, …) are copied verbatim — base-trimming doesn't
  change them.

If a known per-base tag's length doesn't match the sequence (malformed input), it
is left untouched and the run prints a one-line advisory.

### FASTQ example

```bash
whittle -i reads.fastq.gz -o trimmed.fastq.gz \
  -l 500 -q 10 \
  -H 20 -T 20 \
  --qual-trim 8 \
  -t 8
```

The length/quality/GC filters (`-l`/`-q`) are applied to the whole read
first; then trimming runs — cropping 20 bases off each end (`-H`/`-T`) and
trimming any remaining low-quality edges below Q8 (`--qual-trim`); each
resulting output segment must still satisfy `-l`/`--min-length` — all using 8
worker threads (`-t`).

### Folder input

`-i` accepts a single file **or a directory**. When given a directory,
`whittle` merges every read file directly in it (`.fastq`/`.fq`/
`.fastq.gz`/`.fq.gz`, or `.bam`) — in sorted filename order — into one
trimmed output (`-o` file or stdout). The folder must be one format (all
FASTQ-family, or all BAM); non-recursive (subdirectories are ignored); a
mixed or empty folder is an error. In folder mode the format is detected
per file by extension, so `--in-format` has no effect on a directory input
(`--out-format`/`-o`'s extension still control the output).

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --qual-trim 10
```

### uBAM example

```bash
whittle -i reads.ubam.bam -o trimmed.ubam.bam \
  -l 1000 \
  -H 10 -T 10 \
  --qual-split 9 --qual-split-window 50
```

Same filtering/trimming model, but for unaligned BAM: `MM`, `ML`, and `MN`
are rebuilt for every output read so the modification calls still line up
with the (possibly shorter, possibly split-into-several) output sequence.
See **The MM/ML/ML guarantee** below.

### Flag reference

| Flag | Meaning |
|---|---|
| `-i`, `--input <PATH>` | Input file (omit for stdin) |
| `-o`, `--output <PATH>` | Output file (omit for stdout) |
| `--in-format`, `--out-format {fastq,fastq-gz,bam}` | Force format instead of detecting it |
| `--fastq-tags {all,none,LIST}` | Aux tags to carry into FASTQ headers on BAM→FASTQ (default `all`; MM/ML/MN reconstructed, per-base kinetics sliced, `mv` dropped on trim, rest verbatim) |
| `-c`, `--compression-level <0–9>` | DEFLATE level for compressed output — bgzf for BAM, gzip for FASTQ.gz (default 6). Lower is faster/larger; ignored for plain FASTQ |
| `-t`, `--threads <N>` | Total worker threads (default: all detected CPUs; values above the CPU count are clamped down to it), split workload-aware across the decode/render/encode stages; applies to both the FASTQ and uBAM pipelines |
| `-l`, `--min-length <N>` | Minimum read length to keep (default 1) — also the minimum length for a *split segment* to be kept, see below |
| `-L`, `--max-length <N>` | Maximum read length to keep |
| `-q`, `--min-qual <F>` | Minimum read quality to keep (default 0) |
| `-Q`, `--max-qual <F>` | Maximum read quality to keep (default 1000) |
| `-g`, `--min-gc <F>` | Minimum GC fraction (0–1) to keep |
| `-G`, `--max-gc <F>` | Maximum GC fraction (0–1) to keep |
| `-m`, `--qual-mode {mean,arithmetic,median}` | How a read's quality is summarized for `-q`/`-Q` (default `mean`, the ONT-standard error-probability mean; `arithmetic` averages Phred scores directly; `median` is the Phred median) |
| `-H`, `--head-crop <N>` | Trim N bases off the start of every read |
| `-T`, `--tail-crop <N>` | Trim N bases off the end of every read |
| `--qual-trim <Q>` | Trim low-quality bases off both ends down to the first/last base >= Q |
| `--qual-best-segment <Q>` | Keep only the single longest contiguous run of quality >= Q |
| `--qual-split <Q>` | Split the read at low-quality (< Q) runs, keeping each surviving segment as its own record |
| `--qual-split-window <N>` | Smoothing window for `--qual-split` (default 1): a low-quality run shorter than this is tolerated rather than causing a split |
| `--update-moves` | Keep ONT signal tags (`mv`/`ts`/`ns`/`sp`/`pi`) consistent through trimming for signal-aware tools (Remora, Clair3 v2) instead of dropping them. BAM→BAM only |
| `-a`, `--adapter-fasta <FILE>` | Custom adapter/primer FASTA (sequences shorter than the 11 bp minimum match length are skipped with a warning). Enables adapter trimming. See **Adapter trimming** below |
| `--adapter-preset {ont}` | Use the built-in ONT catalog instead of (or alongside) `--adapter-fasta`. Enables adapter trimming |
| `--adapter-error-rate <F>` | End-match error tolerance as a fraction of adapter length (default 0.2); interior/chimera-split hits use half this budget |
| `--adapter-end-size <N>` | Bases at each read end searched for a terminal adapter (default 150) |
| `--adapter-ends-only` | Trim adapters at read ends only; never split on an interior adapter |
| `--adapter-sample <N>` | Reads to sample (the first N, not random) for adapter presence detection with `--adapter-preset` (default 10000); trimming then runs only against adapters actually found in the sample (faster, fewer spurious trims), falling back to the full set if detection finds none. Always disabled (`0`) when `--adapter-fasta` is given — a custom FASTA is searched in full. `0` disables detection outright. Inputs under 100 reads always skip detection, regardless of this setting |
| `-v`, `-vv` | Increase logging detail: `-v` = debug, `-vv` = trace (default: info). See **Logging & progress** below |
| `--quiet` | Silence progress and the info-level summary; warnings and errors still print |

`-H`/`-T` are a positional fixed crop and always run first, before any
quality-based operation, on whatever remains of the read. `--qual-trim`,
`--qual-best-segment`, and `--qual-split` are three different quality-trimming
strategies and are **mutually exclusive** — pass at most one.

When a read is split into segments, each surviving segment's name gets a
`_segment_N` suffix (1-based); `-l/--min-length` filters out segments (not
just whole reads) that end up too short after trimming.

### Logging & progress

Logging level is `-v`/`-vv` (debug/trace) or `--quiet` (warnings/errors
only), with the `WHITTLE_LOG` environment variable available as a
`RUST_LOG`-style override (e.g. `WHITTLE_LOG=debug`, or a per-module
filter like `WHITTLE_LOG=whittle::pipeline=trace`). Precedence: `WHITTLE_LOG`
overrides `-v`/`-vv`, but `--quiet` always wins over `WHITTLE_LOG`. All of
this is on stderr — stdout carries only the read data. Progress itself
renders as a live bar/spinner when stderr is a terminal, or as periodic
log lines (every ~30s) when stderr is redirected to a file or pipe.

## The MM/ML/MN guarantee

Long-read base modification calls (5mC, 6mA, etc.) are stored per-read as
`MM` (which bases are modified, encoded as skip-counts) and `ML` (the
modification's probability). Naively trimming a read's `SEQ`/`QUAL` without
touching these tags produces a BAM that decodes to nonsense: `MM`'s
skip-counts still refer to base occurrences in the *original* sequence, not
the trimmed one.

`whittle` treats this as a first-class part of trimming: for **every**
output uBAM read — head/tail-cropped, quality-trimmed, or split into several
segments — the `MM` and `ML` tags are rebuilt from scratch against the
output window, renumbering skip-counts and re-slicing probability bytes so
they describe exactly the bases that survived. `MN` (the modification base
count) is updated to match, or dropped along with `MM`/`ML` if trimming
removed every modified base from a read. Everything else in the record
(other aux tags, flags, mapping-quality placeholder, etc.) rides through
unchanged.

This is validated by two decode-equivalence tests that independently
re-decode `whittle`'s output with `rust-htslib`'s `basemods_iter()` (a
separate MM/ML implementation from the one `whittle` uses to write) and
assert the decoded calls match the original read's calls filtered to the
surviving window:

- `tests/bam_mods_oracle.rs::trimmed_output_mods_match_oracle` — always runs,
  against a small synthetic fixture.
- `tests/bam_mods_oracle.rs::real_ubam_oracle_sweep` — `#[ignore]`d by
  default; opt in with a real uBAM:

  ```bash
  WHITTLE_UBAM=/path/to/real.ubam \
    cargo test --test bam_mods_oracle -- --ignored
  ```

  This runs a fixed head/tail crop over every read in the file and checks
  every output read's modification calls against the original, read by read.

## Adapter trimming

Adapter trimming is **opt-in and off by default** — every existing invocation
of `whittle` behaves exactly as before. Turn it on with `-a`/`--adapter-fasta
<FILE>` (your own adapter/primer sequences, one per FASTA record — each must
be >= the 11 bp minimum match length; shorter entries are skipped with a
warning) and/or `--adapter-preset ont` (the built-in ONT catalog, below); the
two can be combined, and either one alone is enough to enable the feature.

Once enabled, every adapter is searched for on **both strands** (each
sequence is also matched reverse-complemented, so it's found regardless of
read orientation), and `whittle` does two things per read:

- **Terminal trimming** — an adapter matching within `--adapter-end-size`
  bases of an end (default 150) is trimmed off, walking the keep-boundary
  inward past the match.
- **Chimera splitting** — an adapter matching in the read's interior (away
  from both ends) is treated as a chimera junction: the read is split there,
  the adapter excised, and each surviving side kept as its own segment. Pass
  `--adapter-ends-only` to disable this and only trim ends — this also skips
  searching the interior entirely, so each adapter lookup only scans the two
  end-zones (cheaper than the full-window scan chimera splitting requires).

Catalog entries tagged 5'/3'/both only gate which end is checked for a
*terminal* trim — any adapter can still trigger an interior split, since a
front/rear adapter sequence found mid-read is itself the chimera signal.
Custom `--adapter-fasta` sequences are always checked at both ends.

Interior/chimera hits use a **stricter, derived error budget** than terminal
hits: `--adapter-error-rate` (default 0.2, i.e. 20% of the adapter's length)
sets the terminal-match tolerance, while an interior match must fall within
*half* that fraction to trigger a split — so a marginal end-match still
trims, but only a tight match splits a read in two.

Trimming/splitting on adapters flows through the same tag-rewrite machinery
as quality- and `-H`/`-T`-trimming: on uBAM, `MM`/`ML`/`MN` are rebuilt for
every resulting window or segment (see **The MM/ML/MN guarantee** above),
and the rest of the trim-aware tag handling (per-base kinetics,
`--update-moves`, etc.) applies identically.

### Presence detection

Adapter/preset catalogs (especially `--adapter-preset ont`, over a hundred
sequences) usually contain far more entries than are actually present in a
given run. By default `whittle` samples the first `--adapter-sample <N>`
reads (default 10000) — the leading prefix, not a random draw — checks which
of the configured adapters actually match within that sample, and trims the
rest of the input against only that reduced set — same error-rate/end-size
settings, fewer adapters to search per read, and no spurious low-confidence
trims from catalog entries that were never in the data to begin with.

Detection is **preset-only**: it is automatically disabled (equivalent to
`--adapter-sample 0`) whenever a custom `--adapter-fasta` is given, since a
user-supplied FASTA is a curated set that should always be searched in full —
sampling could otherwise drop a rare custom adapter. Pass `--adapter-sample 0`
yourself to disable detection for a preset too and trim against the full
configured set, unconditionally.

Because the sample is a prefix rather than a random draw, ordered data (e.g. a
run of clean reads before any adapted ones) can look adapter-free during
detection even though adapters are present later in the file. To guard
against this, if detection keeps **zero** adapters, `whittle` falls back to
the full configured set (logged as a warning) instead of silently disabling
trimming for the rest of the run. Inputs under 100 reads always skip
detection (too few reads to sample reliably) and use the full set regardless
of `--adapter-sample`.

### The built-in ONT catalog

`--adapter-preset ont` loads `src/adapter/ont_catalog.tsv` — a catalog
assembled for `whittle` (not vendored from another project) from
ONT-published primary sources: dorado's `adapter_primer_kits.cpp` (kit-14
ligation adapters plus legacy chemistry), Porechop's `adapters.py` (legacy
adapters and the 96 barcode sequences), and qcat's kit definitions. It
covers ligation adapters (kit-14 and legacy), the rapid adapter, direct-RNA
adapters, PCR/cDNA and 10X primers, barcode flanking sequences, and all 96
native/PCR/rapid barcodes — 124 sequences total after de-duplication (one
pair of entries across sources shares an identical sequence and is folded
together). A barcode is one shared 24 bp oligo across kit families;
`whittle`'s reverse-complement search covers the native orientation without
needing a separate flank+revcomp catalog entry. Sequences shorter than 11 bp
(a handful of construction-only flank fragments) are never searched
standalone, since a pattern that short would match almost anywhere.

### Example

```bash
whittle -i reads.fastq.gz -o trimmed.fastq.gz --adapter-preset ont -t 8
```

Or supply your own sequences instead of (or alongside) the built-in catalog:

```bash
whittle -i reads.fastq.gz -o trimmed.fastq.gz \
  --adapter-fasta my_adapters.fasta --adapter-error-rate 0.15
```

### Build requirement

Adapter matching uses [`sassy`](https://crates.io/crates/sassy) for
approximate (edit-distance) search, which needs its AVX2 SIMD path on
x86-64: `.cargo/config.toml` sets `target-cpu=x86-64-v3` automatically for
`cfg(target_arch = "x86_64")` builds (aarch64 uses NEON by default, no flag
needed). Building `whittle` requires Rust >= 1.91.

## v1 limitations

- **uBAM only — aligned BAM is refused.** `whittle` checks the unmapped
  flag on every record and errors out (naming the offending read) if it
  finds an aligned one. There is no support for trimming alignments in
  place, adjusting CIGAR/POS, or otherwise handling mapped reads.
- **No contamination filtering or ab-initio adapter detection.** Adapter
  trimming (`-a`/`--adapter-preset`, see **Adapter trimming** above) matches
  a known catalog by approximate search; there is no minimap2-based
  host-contamination screen, and `whittle` does not infer unrecognized
  adapter sequences or auto-detect adapter presence (unlike e.g.
  Porechop_ABI).
- **No FASTQ→BAM conversion.** FASTQ in, BAM out is explicitly rejected —
  there's no header/tags to build a BAM record from a bare FASTQ read. The
  reverse, BAM→FASTQ, *is* supported; see [Format conversion](#format-conversion)
  above.
- **BAM over stdin auto-detects.** BAM files are BGZF-compressed, which shares
  gzip's magic bytes (`\x1f\x8b`) with gzipped FASTQ. `whittle` tells them apart
  by the BGZF block signature — the `FEXTRA` flag plus the mandatory `BC`
  subfield — which it checks before the plain-gzip fallback. So a raw BAM stream
  piped over stdin (e.g. `samtools view -b … | whittle`) is detected with no
  hint, just like a `.bam` path, plain FASTQ, or `.fastq.gz`. Pass
  `--in-format bam` only to force the interpretation of an unusual or headerless
  stream.
- **`--min-length` is dual-purpose.** It's both the whole-read minimum-length
  filter and the minimum length for a segment produced by `--qual-split` to
  be kept — there's no separate flag for the two.
- **The three quality-trim operations are mutually exclusive.** Pick one of
  `--qual-trim`, `--qual-best-segment`, `--qual-split` (or none, for filtering
  only). `-H`/`-T` fixed crop is independent of all three and always
  composes with whichever one you pick, running first.

## License

Licensed under the [Apache License, Version 2.0](LICENSE). Copyright 2026 Erdi Kılıç.
