# chopping

A tag-aware long-read (ONT/PacBio) trimmer. `chopping` filters and trims
FASTQ and unaligned BAM (uBAM) reads the same way tools like `chopper` do —
but when the input is uBAM, it also **recomputes the `MM`/`ML`/`MN`
base-modification tags** so a trimmed or split read's methylation calls stay
correct instead of silently drifting out of register with the sequence.

## Install

```bash
cargo build --release
```

The binary is written to `target/release/chopping`. Building links against
`noodles-bgzf`'s `libdeflate` backend for BAM I/O, so no external `htslib` is
required to build or run (an optional `rust-htslib`-based test dependency is
dev-only, used by the oracle tests described below).

## Usage

`chopping` reads from `-i/--input` (or stdin) and writes to `-o/--output` (or
stdout). Format is detected from the file extension, or sniffed from the
first bytes when reading a stream; it can also be forced with
`--in-format`/`--out-format {fastq,fastq-gz,bam}`.

**Output is plain FASTQ by default — it is never auto-compressed.** A `.gz`
input does not imply gzipped output: with no `-o` extension and no
`--out-format`, `chopping` always writes plain FASTQ (this matters most on
stdout, where silently emitting gzip bytes would be surprising and is also
slower). Gzip output only happens when you ask for it explicitly, either via
an `-o` path ending in `.gz` (e.g. `out.fastq.gz`, `out.fq.gz`, or a bare
`out.gz`) or `--out-format fastq-gz`. When it
does, it's produced by a parallel gzip encoder (`gzp`) that uses `-t/--threads`
worker threads, so `-t 8 -o out.fastq.gz` compresses using all 8 threads
instead of a single-threaded bottleneck.

### FASTQ example

```bash
chopping -i reads.fastq.gz -o trimmed.fastq.gz \
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
`chopping` merges every read file directly in it (`.fastq`/`.fq`/
`.fastq.gz`/`.fq.gz`, or `.bam`) — in sorted filename order — into one
trimmed output (`-o` file or stdout). The folder must be one format (all
FASTQ-family, or all BAM); non-recursive (subdirectories are ignored); a
mixed or empty folder is an error. In folder mode the format is detected
per file by extension, so `--in-format` has no effect on a directory input
(`--out-format`/`-o`'s extension still control the output).

```bash
chopping -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --trim-qual 10
```

### uBAM example

```bash
chopping -i reads.ubam.bam -o trimmed.ubam.bam \
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
| `-t`, `--threads <N>` | Worker threads for the FASTQ pipeline (default 4; uBAM is single-threaded) |
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

`chopping` treats this as a first-class part of trimming: for **every**
output uBAM read — head/tail-cropped, quality-trimmed, or split into several
segments — the `MM` and `ML` tags are rebuilt from scratch against the
output window, renumbering skip-counts and re-slicing probability bytes so
they describe exactly the bases that survived. `MN` (the modification base
count) is updated to match, or dropped along with `MM`/`ML` if trimming
removed every modified base from a read. Everything else in the record
(other aux tags, flags, mapping-quality placeholder, etc.) rides through
unchanged.

This is validated by two decode-equivalence tests that independently
re-decode `chopping`'s output with `rust-htslib`'s `basemods_iter()` (a
separate MM/ML implementation from the one `chopping` uses to write) and
assert the decoded calls match the original read's calls filtered to the
surviving window:

- `tests/bam_mods_oracle.rs::trimmed_output_mods_match_oracle` — always runs,
  against a small synthetic fixture.
- `tests/bam_mods_oracle.rs::real_ubam_oracle_sweep` — `#[ignore]`d by
  default; opt in with a real uBAM:

  ```bash
  CHOPPING_UBAM=/path/to/real.ubam \
    cargo test --test bam_mods_oracle -- --ignored
  ```

  This runs a fixed head/tail crop over every read in the file and checks
  every output read's modification calls against the original, read by read.

## v1 limitations

- **uBAM only — aligned BAM is refused.** `chopping` checks the unmapped
  flag on every record and errors out (naming the offending read) if it
  finds an aligned one. There is no support for trimming alignments in
  place, adjusting CIGAR/POS, or otherwise handling mapped reads.
- **No contamination filtering.** There is no minimap2-based adapter/barcode/
  host-contamination screen (unlike e.g. Porechop_ABI). `chopping` only
  filters and trims by length, quality, and GC content.
- **No cross-format conversion.** BAM in, FASTQ out (or vice versa) is
  explicitly rejected — `chopping` trims BAM to BAM and FASTQ to FASTQ, not
  across formats.
- **BAM over stdin needs an explicit hint.** BAM files are BGZF-compressed,
  which shares gzip's magic bytes (`\x1f\x8b`) with the bytes `chopping` uses
  to auto-detect gzipped FASTQ. A `.bam` file *path* auto-detects fine (the
  extension is checked before any byte-sniffing), but **a raw BAM stream
  piped over stdin, with no path to check the extension of, will be
  misdetected as gzipped FASTQ** unless you pass `--in-format bam`
  explicitly. Plain FASTQ and `.fastq.gz` piped over stdin auto-detect
  correctly without any hint.
- **`--min-length` is dual-purpose.** It's both the whole-read minimum-length
  filter and the minimum length for a segment produced by `--split-qual` to
  be kept — there's no separate flag for the two.
- **The three quality-trim operations are mutually exclusive.** Pick one of
  `--trim-qual`, `--best-segment`, `--split-qual` (or none, for filtering
  only). `-H`/`-T` fixed crop is independent of all three and always
  composes with whichever one you pick, running first.
