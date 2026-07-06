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
worker threads, so `-t 8 -o out.fastq.gz` compresses using all 8 threads
instead of a single-threaded bottleneck.

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
  - a `--split-qual` **split** emits dorado-style subreads (`pi` = parent id,
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
  --trim-qual 8 \
  -t 8
```

The length/quality/GC filters (`-l`/`-q`) are applied to the whole read
first; then trimming runs — cropping 20 bases off each end (`-H`/`-T`) and
trimming any remaining low-quality edges below Q8 (`--trim-qual`); each
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
whittle -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --trim-qual 10
```

### uBAM example

```bash
whittle -i reads.ubam.bam -o trimmed.ubam.bam \
  -l 1000 \
  -H 10 -T 10 \
  --split-qual 9 --split-window 50
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
| `-t`, `--threads <N>` | Total worker threads (default 4), split workload-aware across the decode/render/encode stages; applies to both the FASTQ and uBAM pipelines |
| `-l`, `--min-length <N>` | Minimum read length to keep (default 1) — also the minimum length for a *split segment* to be kept, see below |
| `-L`, `--max-length <N>` | Maximum read length to keep |
| `-q`, `--min-qual <F>` | Minimum read quality to keep (default 0) |
| `-Q`, `--max-qual <F>` | Maximum read quality to keep (default 1000) |
| `-g`, `--min-gc <F>` | Minimum GC fraction (0–1) to keep |
| `-G`, `--max-gc <F>` | Maximum GC fraction (0–1) to keep |
| `-m`, `--qual-mode {mean,arithmetic,median}` | How a read's quality is summarized for `-q`/`-Q` (default `mean`, the ONT-standard error-probability mean; `arithmetic` averages Phred scores directly; `median` is the Phred median) |
| `-H`, `--head-crop <N>` | Trim N bases off the start of every read |
| `-T`, `--tail-crop <N>` | Trim N bases off the end of every read |
| `--trim-qual <Q>` | Trim low-quality bases off both ends down to the first/last base >= Q |
| `--best-segment <Q>` | Keep only the single longest contiguous run of quality >= Q |
| `--split-qual <Q>` | Split the read at low-quality (< Q) runs, keeping each surviving segment as its own record |
| `--split-window <N>` | Smoothing window for `--split-qual` (default 1): a low-quality run shorter than this is tolerated rather than causing a split |
| `--update-moves` | Keep ONT signal tags (`mv`/`ts`/`ns`/`sp`/`pi`) consistent through trimming for signal-aware tools (Remora, Clair3 v2) instead of dropping them. BAM→BAM only |

`-H`/`-T` are a positional fixed crop and always run first, before any
quality-based operation, on whatever remains of the read. `--trim-qual`,
`--best-segment`, and `--split-qual` are three different quality-trimming
strategies and are **mutually exclusive** — pass at most one.

When a read is split into segments, each surviving segment's name gets a
`_segment_N` suffix (1-based); `-l/--min-length` filters out segments (not
just whole reads) that end up too short after trimming.

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

## v1 limitations

- **uBAM only — aligned BAM is refused.** `whittle` checks the unmapped
  flag on every record and errors out (naming the offending read) if it
  finds an aligned one. There is no support for trimming alignments in
  place, adjusting CIGAR/POS, or otherwise handling mapped reads.
- **No contamination filtering.** There is no minimap2-based adapter/barcode/
  host-contamination screen (unlike e.g. Porechop_ABI). `whittle` only
  filters and trims by length, quality, and GC content.
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
  filter and the minimum length for a segment produced by `--split-qual` to
  be kept — there's no separate flag for the two.
- **The three quality-trim operations are mutually exclusive.** Pick one of
  `--trim-qual`, `--best-segment`, `--split-qual` (or none, for filtering
  only). `-H`/`-T` fixed crop is independent of all three and always
  composes with whichever one you pick, running first.

## License

Licensed under the [Apache License, Version 2.0](LICENSE). Copyright 2026 Erdi Kılıç.
