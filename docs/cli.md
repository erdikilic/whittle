# Command-line reference

Every flag whittle accepts, plus how it picks formats, merges folders, and reports
progress. `whittle --help` prints the same option list grouped by section.

## Input, output, and formats

whittle reads from `-i`/`--input` (or stdin) and writes to `-o`/`--output` (or
stdout). Either flag also takes `-` for the standard stream, the usual pipeline
spelling. It takes the format from the file extension, sniffs it from the first
bytes of a stream, or accepts it from `--in-format`/`--out-format
{fastq,fastq-gz,fastq-bgz,bam}`.

When a downstream reader closes the pipe early (`whittle ... | head`), whittle
stops writing and exits 0 without a message; the reader's decision to stop is not
a failure of the run.

Output is plain FASTQ by default and is never compressed on its own. A `.gz` or
`.bgz` input does not imply compressed output; you get that only by asking for
it, with a `.gz`/`.bgz` output path or the matching format flag. Compressed
output is written by a parallel encoder, which takes its share of the `-t`
budget alongside reading and trimming. BGZF FASTQ input is decompressed
block-parallel too, unless adapter trimming is on, in which case the budget goes
to trimming instead; ordinary gzip stays a serial input format.

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
| `-i, --input <PATH>` | Input file or directory (omit, or pass `-`, for stdin) |
| `-o, --output <PATH>` | Output file (omit, or pass `-`, for stdout) |
| `--in-format`, `--out-format {fastq,fastq-gz,fastq-bgz,bam}` | Force a format instead of detecting it |
| `--fastq-tags {all,none,LIST}` | Aux tags to carry into FASTQ headers on BAM→FASTQ (default `all`) |
| `-c, --compression-level <0-9>` | DEFLATE level for compressed output (default 4 for gzip FASTQ, 6 for BGZF and BAM); ignored for plain FASTQ |
| `--summary-json <PATH>` | Write a machine-readable run summary to PATH (see below) |
| `-t, --threads <N>` | Worker threads, at least 1 (default: all detected CPUs, clamped to that max) |
| `--ordered` | Write records in input order under `-t > 1`; without it output is in completion order (faster, not reproducible) |
| `-l, --min-length <N>` | Minimum length to keep, per output segment (default 1) |
| `-L, --max-length <N>` | Maximum length to keep |
| `-q, --min-qual <F>` | Minimum read quality, a finite value of at least 0 (default 0) |
| `-Q, --max-qual <F>` | Maximum read quality, a finite value of at least 0 (default 1000) |
| `-g, --min-gc <F>`, `-G, --max-gc <F>` | GC-fraction bounds (0 to 1) |
| `-m, --qual-mode {mean,arithmetic,median}` | How read quality is summarized (default `mean`, the error-probability mean) |
| `-H, --head-crop <N>`, `-T, --tail-crop <N>` | Fixed crop from each end; always runs first |
| `--qual-trim <Q>` | Trim low-quality bases from both ends down to the first base >= Q |
| `--qual-best-segment <Q>` | Keep only the longest contiguous run of quality >= Q |
| `--qual-split <Q>` | Split at low-quality (< Q) runs, keeping each surviving segment |
| `--qual-split-window <N>` | Tolerate low-quality runs shorter than N without splitting (default 1); requires `--qual-split` |
| `--update-moves` | Rewrite ONT signal tags through trimming instead of dropping them (BAM→BAM) |
| `-a, --adapter-fasta <FILE>` | Adapter/primer FASTA; enables adapter trimming |
| `--adapter-preset {none,ont}` | Built-in adapter catalog (default `none`; `ont` enables trimming) |
| `--adapter-error-rate <F>` | End-match tolerance as a fraction of adapter length (default 0.2); requires an adapter source |
| `--adapter-end-size <N>` | End-zone width searched for terminal adapters (default 150); requires an adapter source |
| `--adapter-ends-only` | Trim ends only; never split on an interior adapter |
| `--adapter-sample <N>` | Reads sampled for preset detection or inference (defaults `0` and `40000`, respectively); requires an adapter source |
| `--adapter-infer [trim\|report]` | Discover adapters de novo; omitted value defaults to `trim` |
| `--adapter-infer-policy {conservative,aggressive}` | Trust policy for inferred adapters (default `conservative`) |
| `-v`, `-vv` | Stage detail, then per-read decisions; higher counts are rejected |
| `--progress {auto,bar,plain,none}` | How to report progress, independently of the log level (default `auto`) |
| `--quiet` | Silence progress and the summary; warnings and errors still print. Conflicts with `-v` and `--progress` |

`--qual-trim`, `--qual-best-segment`, and `--qual-split` are three strategies for
the same step, so pass at most one. `-H`/`-T` are independent and compose with
whichever you pick.

An adapter source is `--adapter-fasta`, `--adapter-preset ont`, or
`--adapter-infer`. The tuning flags `--adapter-error-rate`, `--adapter-end-size`,
and `--adapter-sample` are rejected without one, since they would otherwise be
accepted and ignored.

## Machine-readable summary

`--summary-json <PATH>` writes one JSON object describing the run: the resolved
settings under `params`, and the counters under `reads`, `bases`, and
`segments_dropped`. It is written on every dispatch path, including folder merges,
and regardless of `--quiet` or the log level, so a workflow manager always gets
the file it asked for. A write failure fails the run rather than leaving a stale
file from a previous invocation in place.

```bash
whittle -i reads.bam -o trimmed.fastq.gz -l 500 --quiet --summary-json qc.json
```

```json
{
  "schema_version": 1,
  "tool": "whittle",
  "version": "0.1.1",
  "command": "whittle -i reads.bam -o trimmed.fastq.gz -l 500 --quiet --summary-json qc.json",
  "input": "reads.bam",
  "output": "trimmed.fastq.gz",
  "elapsed_seconds": 12.34,
  "params": { "threads": 8, "ordered": false, "min_length": 500, "qual_mode": "mean", "quality_op": null,
              "adapters": { "configured": 120, "count": 4, "sample": 500, "infer": "off" } },
  "reads": { "input": 1000, "output": 950, "with_output": 940, "trimmed_to_nothing": 30, "all_filtered": 30 },
  "bases": { "input": 10000000, "output": 9500000 },
  "segments_dropped": { "too_short": 12, "too_long": 0, "low_quality": 5, "high_quality": 0, "gc_out_of_range": 0 },
  "warnings": { "malformed_tag_reads": 0, "malformed_mod_reads": 0 }
}
```

`params` is abbreviated above; the real file carries every resolved setting,
including the ones left at their defaults. `params.ordered` records whether a
multithreaded run wrote records in input order.

Under `warnings`, `malformed_tag_reads` counts reads whose per-base kinetics tag
length disagreed with the sequence and was left untouched, and
`malformed_mod_reads` counts reads whose MM/ML/MN modification block could not be
parsed and was removed from the output. Both are also reported at the end of the
run on stderr.

`reads.output` counts output segments, not input reads, so under `--qual-split` it
can legitimately exceed `reads.input`. The three read-level buckets
(`with_output`, `trimmed_to_nothing`, `all_filtered`) do partition `reads.input`.

Under `params.adapters`, `configured` is the set asked for (the preset and/or
FASTA) and `count` is the set actually trimmed against, after presence detection
narrowed it or inference replaced it. They are equal when neither ran, and the
startup banner prints `configured`, since that is all that is known before reads
have been sampled. Under `--adapter-infer` nothing is configured up front, so
`configured` is `0`.

`schema_version` is bumped only when an existing field changes meaning or
disappears. New fields can appear without a bump, so parse leniently.

## Man page

Release tarballs ship `man/whittle.1` alongside the binary, and the same file is
checked into the repository, so you can install it without building:

```bash
install -Dm644 man/whittle.1 /usr/share/man/man1/whittle.1
```

## Logging and progress

Set the log level with `-v`/`-vv` or `--quiet` (warnings and errors only).

`-v` adds the resolved stage detail: the detected input format and how long
detection took, the thread budget actually handed to each stage, and the read and
base counts when processing finishes. Every field is emitted as structured data
rather than prose, so a filter or a log collector can read it.

Every line has the same shape: a timestamp, the level, the enclosing span with
the fields identifying it, a capitalized message, and then the structured fields.
The message is prose and the field values are data, so `Segment dropped
reason="too short"` is greppable on either half.

`-vv` adds the per-read decisions, each attributed to the read that produced it:

```text
[2026-09-01 12:13:45] [TRACE] [read{name=read_adapter}] Adapter hit adapter="LSK109_front" start=0 end=28 cost=0 action="trim 5'"
[2026-09-01 12:13:45] [TRACE] [read{name=read_adapter}] Segment kept segment=1 of=1 start=28 end=217 len=189
[2026-09-01 12:13:45] [TRACE] [read{name=read_short}] Segment dropped segment=1 of=1 start=0 end=30 len=30 reason="too short"
[2026-09-01 12:13:45] [TRACE] [read{name=read_short}] Every segment filtered produced=1
```

That is how to answer why a particular read was cut where it was, or why it is
missing from the output: which adapter matched, over what span, how far off an
exact match it was, and what that made whittle do. It is verbose by design, one
group of lines per read, so redirect it or filter it.

`WHITTLE_LOG` overrides the level with a `RUST_LOG`-style filter, for example
`WHITTLE_LOG=whittle::adapter=trace` to see adapter decisions without the
per-segment lines, and `--quiet` still wins over it. A value that does not parse
as a filter is reported as a warning and the level falls back to `-v`/default,
so a typo in the variable never silences the run.

All logging goes to stderr, so stdout carries only read data. By default progress
shows as a live bar when stderr is a terminal, and as periodic lines (about every
30s, or 10s under `-v`) when it is redirected to a file or a pipe.
`WHITTLE_PROGRESS_INTERVAL` overrides that cadence, in whole seconds. The bar is
never written to a non-terminal, so a redirected log holds no escape sequences or
carriage returns.

`--progress` chooses that independently of the log level:

| value | behavior |
|---|---|
| `auto` | A bar on a terminal, periodic lines otherwise. The default. |
| `bar` | The animated bar, even when redirected. Falls back to periodic lines under `-v`/`-vv` or `WHITTLE_LOG`, since the bar hides the multi-line banner a verbose run asks for and debug lines cannot share a terminal with it. |
| `plain` | Always periodic lines, never a bar. |
| `none` | No progress reporting. The banner, warnings and summary still print. |

`--progress none` is the one to reach for in a pipeline that wants the run
summary in its log without a progress line every interval. `--quiet` is stronger,
dropping the summary too, and conflicts with `--progress` rather than overriding
it.
