# Command-line reference

Every flag whittle accepts, plus how it picks formats, merges folders, and reports
progress. `whittle --help` prints the same option list grouped by section.

## Input, output, and formats

whittle reads from `-i`/`--input` (or stdin) and writes to `-o`/`--output` (or
stdout). It takes the format from the file extension, sniffs it from the first
bytes of a stream, or accepts it from `--in-format`/`--out-format
{fastq,fastq-gz,fastq-bgz,bam}`.

Output is plain FASTQ by default and is never compressed on its own. A `.gz` or
`.bgz` input does not imply compressed output; you get that only by asking for
it, with a `.gz`/`.bgz` output path or the matching format flag. Compressed
output is written by a parallel encoder using `-t` threads. BGZF FASTQ input is
decompressed block-parallel too; ordinary gzip stays a serial input format.

With no output extension or `--out-format`, the output format mirrors the input,
except that compressed FASTQ input defaults back to plain FASTQ. FASTQ→BAM
isn't supported; there's no header to build a BAM record from. BGZF streams are
recognized by their decompressed payload, so piped FASTQ.bgz and `samtools view
-b … | whittle` need no hint.

On BAM→FASTQ, aux tags go into the FASTQ header tab-delimited, following the
`samtools fastq -T` convention. `--fastq-tags` picks which ones: `all` (default),
`none`, or a list like `MM,ML,RG`. `MM`/`ML`/`MN` are reconstructed for the
trimmed segment, per-base tags are sliced, and everything else is copied
verbatim.

## Folder input

`-i` also takes a directory. whittle merges every read file directly inside it,
in sorted filename order, into one output. The folder has to be a single format
(all FASTQ-family or all BAM); subdirectories are ignored, and a mixed or empty
folder is an error.

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --qual-trim 10
```

Naming the output inside the input directory is a hard error: whittle cannot
tell a real input from a stale prior output, and merging over either loses data.

## Options

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
| `--qual-trim <Q>` | Trim low-quality bases from both ends down to the first base >= Q |
| `--qual-best-segment <Q>` | Keep only the longest contiguous run of quality >= Q |
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

`--qual-trim`, `--qual-best-segment`, and `--qual-split` are three strategies for
the same step, so pass at most one. `-H`/`-T` are independent and compose with
whichever you pick.

## Logging and progress

Set the log level with `-v`/`-vv` (debug/trace) or `--quiet` (warnings and errors
only). `WHITTLE_LOG` overrides it with a `RUST_LOG`-style filter, for example
`WHITTLE_LOG=whittle::workflow=trace`, and `--quiet` still wins over it.

All logging goes to stderr, so stdout carries only read data. Progress shows as a
live bar when stderr is a terminal, or as periodic lines (about every 30s, or 10s
under `-v`) when it's redirected to a file or pipe.
