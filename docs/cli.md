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
`.bgz` input does not imply compressed output; compressed output is produced
only when requested, with a `.gz`/`.bgz` output path or the matching format
flag. Compressed output is written by a parallel encoder, which takes its share
of the `-t` budget alongside reading and trimming. BGZF FASTQ input is
decompressed block-parallel too, unless adapter trimming is on, in which case the
budget goes to trimming instead; ordinary gzip stays a serial input format.

With no output extension or `--out-format`, the output format mirrors the input,
except that compressed FASTQ input defaults back to plain FASTQ. FASTQ-to-BAM is
not supported: there is no header to build a BAM record from. BGZF streams are
recognized by their decompressed payload, so piped FASTQ.bgz and `samtools view
-b ... | whittle` need no hint.

On BAM-to-FASTQ, aux tags go into the FASTQ header tab-delimited, following the
`samtools fastq -T` convention. `--fastq-tags` picks which ones: `all` (default),
`none`, or a list like `MM,ML,RG`. `MM`/`ML`/`MN` are reconstructed for the
trimmed segment, per-base tags are sliced, and everything else is copied
verbatim.

## Folder input

`-i` also takes a directory. whittle merges every read file directly inside it,
in natural filename order (digit runs compared by value), into one output. The
folder has to be a single format (all FASTQ-family or all BAM); hidden files and
subdirectories are ignored, and a mixed or empty folder is an error.

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
| `--fastq-tags {all,none,LIST}` | Aux tags to carry into FASTQ headers on BAM-to-FASTQ (default `all`) |
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
| `-H, --head-crop <N>`, `-T, --tail-crop <N>` | Fixed crop from each end; runs before adapter and quality trimming |
| `--qual-trim <Q>` | Trim low-quality bases from both ends down to the first base >= Q |
| `--qual-best-segment <Q>` | Keep only the longest contiguous run of quality >= Q |
| `--qual-split <Q>` | Split at low-quality (< Q) runs, keeping each surviving segment |
| `--qual-split-window <N>` | Tolerate low-quality runs shorter than N without splitting (default 1); requires `--qual-split` |
| `--trim-barcodes` | Remove the barcode spans dorado recorded in the `bi` aux tag, before every other stage (BAM input) |
| `--update-moves` | Rewrite ONT signal tags through trimming instead of dropping them (BAM-to-BAM) |
| `--remove-tag <TAG>` | Remove this two-character aux tag from every output record; repeatable (BAM input) |
| `--strip-kinetics` | Remove the per-base kinetics and alignment-count arrays `ip pw fi fp ri rp sa sm sx` (BAM input) |
| `-a, --adapter-fasta <FILE>` | Adapter/primer FASTA; enables adapter trimming |
| `--adapter-preset {none,ont}` | Built-in adapter catalog (default `none`; `ont` enables trimming) |
| `--adapter-error-rate <F>` | End-match tolerance as a fraction of adapter length (default 0.2); requires an adapter source |
| `--adapter-end-size <N>` | End-zone width searched for terminal adapters (default 150); requires an adapter source |
| `--adapter-ends-only` | Trim ends only; never split on an interior adapter |
| `--adapter-sample <N>` | Reads sampled for preset detection or inference (defaults `0` and `40000`, respectively); requires an adapter source |
| `--adapter-infer [trim\|report]` | Discover adapters de novo; omitted value defaults to `trim` |
| `--adapter-infer-policy {conservative,aggressive}` | Trust policy for inferred adapters (default `conservative`); requires `--adapter-infer` |
| `-v`, `-vv` | Stage detail, then per-read decisions; higher counts are rejected |
| `--progress {auto,bar,plain,none}` | How to report progress, independently of the log level (default `auto`) |
| `--quiet` | Silence progress and the summary; warnings and errors print regardless. Conflicts with `-v` and `--progress` |

`--qual-trim`, `--qual-best-segment`, and `--qual-split` are three strategies for
the same step, so at most one is accepted. `-H`/`-T` are independent and compose
with any of them.

The trimming stages run in a fixed order: barcode removal, then the fixed crop,
then adapter trimming and chimera splitting, then the quality strategy. Each
stage sees what the previous one left, so `--head-crop` counts from the first
base after the front barcode, and the report separates the stages rather than
summing them.

## Barcode trimming

`--trim-barcodes` removes the barcode spans dorado recorded in the `bi` aux tag.
It reads the positions rather than searching for the sequences, so the cut is
exactly the one `dorado demux` would have made, and it runs through the same
trimming machinery as every other stage: `MM`/`ML`/`MN`, per-base kinetics and
the ONT move table stay in register with the trimmed sequence.

```bash
whittle -i barcoded.bam -o trimmed.bam --trim-barcodes --update-moves
```

The tag holds seven floats, of which four are positions: the front barcode's
start and length, and the rear barcode's end and length. A barcode dorado did
not find is written as a negative position and leaves that end of the read
alone, so a read barcoded at one end only is trimmed at that end only. A record
without `bi` passes through untouched.

`bi` itself is dropped from a trimmed read, since its positions index the
untrimmed sequence; the barcode call (`BC`, `bv`) is a per-read label and is
kept. A `bi` that is not a seven-element float array, or whose positions
describe an empty, inverted or out-of-range window, leaves the read untrimmed
and is counted under `warnings.barcode_tag_malformed_reads`.

The positions come from a BAM aux tag, so the flag requires BAM input and is
refused on FASTQ rather than accepted and ignored.

An adapter source is `--adapter-fasta`, `--adapter-preset ont`, or
`--adapter-infer`. The tuning flags `--adapter-error-rate`, `--adapter-end-size`,
and `--adapter-sample` are rejected without one, since they would otherwise be
accepted and ignored.

## Removing tags

`--remove-tag <TAG>` removes one two-character aux tag from every output record.
It is repeatable, and each value has to be exactly two alphanumeric characters,
checked before the run rather than silently matching nothing.

`--strip-kinetics` removes the per-base kinetics and alignment-count arrays in a
single flag: `ip`, `pw`, `fi`, `fp`, `ri`, `rp`, `sa`, `sm` and `sx`. It is
equivalent to naming each of them with `--remove-tag`, and both flags fill one
removal set, so they combine.

```bash
whittle -i reads.bam -o smaller.bam --strip-kinetics --remove-tag ML
```

Removal happens after whittle has rewritten the tags it keeps in register, so
removing one whittle maintains leaves the rest of the record intact: the record
is written without that tag rather than with a stale one. Removing `MM` keeps
the rebuilt `ML` and `MN`, and removing one per-base array still slices the
others to the trimmed window.

Removal applies on every BAM output path. A run that trims nothing normally
writes records back without decoding them; when tags are removed each record is
rebuilt instead, which costs the decode and changes nothing else about the
output. On BAM-to-FASTQ the removal applies to the tags carried into the header,
after `--fastq-tags` has chosen which tags are carried at all.

Both flags name BAM auxiliary tags, so they require BAM input and are refused on
FASTQ rather than accepted and ignored. Removing tags is a complete run on its
own: `whittle -i in.bam -o out.bam --strip-kinetics` with no trimming options is
valid and only strips tags. The resolved set is recorded under
`params.remove_tags` in the summary JSON, with `params.strip_kinetics` recording
which flag asked for it.

## Machine-readable summary

`--summary-json <PATH>` writes one JSON object describing the run: the resolved
settings under `params`, and the counters under `reads`, `bases`, and
`segments_dropped`. It is written on every dispatch path, including folder merges,
and regardless of `--quiet` or the log level. A write failure fails the run
rather than leaving a stale file from a previous invocation in place.

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
  "warnings": { "malformed_tag_reads": 0, "malformed_mod_reads": 0, "barcode_tag_malformed_reads": 0 }
}
```

`params` is abbreviated above; the real file carries every resolved setting,
including the ones left at their defaults. `params.ordered` records whether a
multithreaded run wrote records in input order.

Under `warnings`, `malformed_tag_reads` counts reads whose per-base kinetics tag
length disagreed with the sequence and was left untouched,
`malformed_mod_reads` counts reads whose MM/ML/MN modification block could not be
parsed and was removed from the output, and `barcode_tag_malformed_reads` counts
reads whose `bi` barcode positions did not describe a window inside the read
under `--trim-barcodes`. All are also reported at the end of the run on stderr.

`reads.output` counts output segments, not input reads, so under `--qual-split` it
can exceed `reads.input`. The three read-level buckets (`with_output`,
`trimmed_to_nothing`, `all_filtered`) partition `reads.input`.

Under `params.adapters`, `configured` is the set asked for (the preset and/or
FASTA) and `count` is the set trimmed against, after presence detection
narrowed it or inference replaced it. They are equal when neither ran, and the
startup banner prints `configured`, since that is all that is known before reads
have been sampled. Under `--adapter-infer` nothing is configured up front, so
`configured` is `0`.

`schema_version` is bumped only when an existing field changes meaning or
disappears. New fields can appear without a bump, so consumers tolerate unknown
fields.

## Man page

Release tarballs ship `man/whittle.1` alongside the binary, and the same file is
checked into the repository, so it can be installed without a build:

```bash
install -Dm644 man/whittle.1 /usr/share/man/man1/whittle.1
```

## Logging and progress

The log level is set with `-v`/`-vv` or `--quiet` (warnings and errors only).

`-v` adds the resolved stage detail: the detected input format and how long
detection took, the thread budget handed to each stage, and the read and
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

These lines answer why a read was cut where it was or why it is missing from the
output: which adapter matched, over what span, at what cost, and the resulting
action. The output is one group of lines per read, so it is normally redirected
or filtered.

`WHITTLE_LOG` overrides the level with a `RUST_LOG`-style filter, for example
`WHITTLE_LOG=whittle::adapter=trace` to see adapter decisions without the
per-segment lines; `--quiet` takes precedence over it. A value that does not parse
as a filter is reported as a warning and the level falls back to the `-v`/`-vv`
setting or the default, so a malformed variable does not silence the run.

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
| `bar` | The animated bar, even when redirected. Falls back to periodic lines under `-v`/`-vv` or `WHITTLE_LOG`, since the bar hides the multi-line banner of a verbose run and debug lines cannot share a terminal with it. |
| `plain` | Always periodic lines, never a bar. |
| `none` | No progress reporting. The banner, warnings and summary print regardless. |

`--progress none` suits a pipeline log that keeps the run summary without a
progress line every interval. `--quiet` is stronger, dropping the summary too,
and conflicts with `--progress` rather than overriding it.
