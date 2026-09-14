# Command-line reference

`whittle --help` prints the option list grouped by section. This page describes
every option, format selection, directory input, the summary file, and logging.

## Input and output

`-i`/`--input` names the input and `-o`/`--output` the output. Either flag may
be omitted or given as `-` for the standard stream. The format comes from the
file extension, from the first bytes of a stream, or from `--input-format` and
`--output-format` (`fastq`, `fastq-gz`, `fastq-bgz`, `bam`).

Without an output extension or `--output-format`, the output format is the input
format, except that compressed FASTQ input yields plain FASTQ. Compressed output
is produced only when requested by a `.gz`/`.bgz` extension or the format flag.
FASTQ-to-BAM is not supported: a FASTQ read carries no header from which to
build a BAM record. BGZF streams are recognized by their decompressed payload,
so BAM and FASTQ.bgz on stdin need no format flag.

When a downstream reader closes the pipe (`whittle ... | head`), whittle stops
writing and exits 0.

On BAM-to-FASTQ, auxiliary tags are appended to the FASTQ header, tab-delimited,
in the `samtools fastq -T` convention. `--fastq-tags` selects them: `all`
(default), `none`, or a list such as `MM,ML,RG`. `MM`/`ML`/`MN` are rewritten for
the trimmed segment, per-base arrays are sliced, and the remaining tags are
copied verbatim.

FASTQ input in the same convention (a tab after the read name, then
`TAG:TYPE:VALUE` fields) is tagged FASTQ. Any such field among the first 100
headers puts the whole input on the BAM-to-FASTQ path: every read's fields are
decoded as SAM aux tags and rewritten per output segment exactly as a uBAM
record's are ([tags.md](tags.md)), `--fastq-tags` selects the output tags, and
`--trim-barcodes`, `--remove-tag` and `--remove-kinetics` apply. A field that
does not parse as a SAM tag fails the run and names the read. A header whose
tab-delimited text is not in this form is copied verbatim.

### Threads

`-t`/`--threads` sets the worker count (default: every detected CPU). Trimming,
record serialization, and output compression run on a pool of that size. BGZF
input (BAM and FASTQ.bgz) is decompressed by an additional set of one quarter of
the workers, at least one; plain gzip input is decompressed serially. The
startup banner reports both figures.

Records are written in completion order under `-t > 1`. `--preserve-order` restores
the input order at a small throughput cost.

## Directory input

`-i` accepts a directory. Every read file directly inside it is merged into one
output in natural filename order (digit runs compared numerically). The files
must share one format family (all FASTQ or all BAM); hidden files and
subdirectories are ignored; a mixed or empty directory is an error. An output
path inside the input directory is rejected, since it cannot be distinguished
from an input on a later run.

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --trim-quality 10
```

## Options

| Flag | Meaning |
|---|---|
| `-h, --help` | Print the option list and examples |
| `--version` | Print the version and exit |
| `-i, --input <PATH>` | Input file or directory (omit, or pass `-`, for stdin) |
| `-o, --output <PATH>` | Output file (omit, or pass `-`, for stdout) |
| `--input-format`, `--output-format <FORMAT>` | Force a format instead of detecting it: `fastq`, `fastq-gz`, `fastq-bgz`, `bam` |
| `--fastq-tags <all\|none\|TAGS>` | Aux tags written into FASTQ headers on BAM or tagged FASTQ input (default `all`) |
| `-c, --compression-level <0-9>` | BGZF level for `.gz`, `.bgz` and BAM output (default 4 for `.gz`, 6 for `.bgz` and BAM); ignored for plain FASTQ |
| `--summary-json <PATH>` | Write a machine-readable run summary to PATH; ignored under `--discover-adapters report` |
| `-t, --threads <N>` | Worker threads, at least 1 (default: all detected CPUs, clamped to that maximum) |
| `--preserve-order` | Write records in input order under `-t > 1` |
| `-l, --min-length <BASES>` | Minimum length to keep, per output segment (default 1) |
| `-L, --max-length <BASES>` | Maximum length to keep |
| `-q, --min-quality <PHRED>` | Minimum post-trim segment quality, a finite value of at least 0 (default 0) |
| `-Q, --max-quality <PHRED>` | Maximum post-trim segment quality, a finite value of at least 0 (default 1000) |
| `-g, --min-gc <FRACTION>`, `-G, --max-gc <FRACTION>` | GC-fraction bounds (0 to 1; `0.4` means 40%) |
| `-m, --quality-mode <MODE>` | Quality calculation for `--min-quality`/`--max-quality` only: `mean` (mean error probability as a Phred score, the default), `arithmetic` (mean of the Phred scores), `median` |
| `-H, --trim-front <BASES>`, `-T, --trim-tail <BASES>` | Fixed crop from each adapter-derived segment after barcode restriction and before quality processing; applied once |
| `--trim-quality <PHRED>` | Trim both ends up to the first base of quality at least PHRED |
| `--best-quality-segment <PHRED>` | Keep the highest-scoring segment using cumulative base-error probabilities and the Phred cutoff (modified Mott); may retain bases below the cutoff |
| `--split-quality <PHRED>` | Split at consecutive bases below PHRED and keep each surviving segment |
| `--split-min-low-quality-bases <BASES>` | Minimum consecutive bases below the splitting threshold required to split; shorter internal stretches are retained and low-quality ends are trimmed (default 1); requires `--split-quality` |
| `--trim-barcodes` | Intersect adapter-derived segments with the retained interval from the original `bi` tag before cropping (BAM or tagged FASTQ input) |
| `--update-signal-tags` | Rewrite ONT signal tags through trimming instead of removing them (BAM-to-BAM) |
| `--remove-tag <TAG>` | Remove a two-character aux tag from every output record; repeatable (BAM or tagged FASTQ input) |
| `--remove-kinetics` | Remove the per-base kinetics and alignment-count arrays `ip pw fi fp ri rp sa sm sx` (BAM or tagged FASTQ input) |
| `-a, --adapter-fasta <FILE>` | Adapter and primer FASTA (IUPAC codes accepted; `primer` or `barcode` in a header description restricts the entry to the read ends); enables adapter trimming |
| `--adapter-preset <KITS>` | Built-in kit presets, comma-separated: `lsk114`, `rad114` (`ulk114`), `rbk114`, `nbd114`, `pcb114` (`pcs114`), `rpb114`, `mab114`, `rna004`, `pacbio`, `ont`, `all`; enables adapter trimming |
| `--adapter-error-rate <FRACTION>` | End-match tolerance as a fraction of adapter length (default 0.2); requires an adapter source |
| `--adapter-end-search <BASES>` | End-zone width searched for terminal adapters (default 150); requires an adapter source |
| `--adapter-ends-only` | Trim adapters at ends only; disable interior adapter splitting independently of quality splitting |
| `--adapter-sample-reads <COUNT>` | Reads inspected for preset presence or adapter discovery (defaults 2000 and 40000; at least 100); `0` disables preset detection and is rejected for discovery; ignored with `--adapter-fasta` unless discovering adapters |
| `--discover-adapters [<ACTION>]` | Discover adapters de novo: `trim` (the bare flag) or `report` (FASTA to stdout, then exit without read output or a JSON summary) |
| `--adapter-discovery-policy <POLICY>` | Trust policy for inferred adapters: `conservative` (default) or `aggressive`; requires `--discover-adapters` |
| `-v, --verbose` (repeatable) | Stage detail with `-v`, per-read decisions with `-vv` |
| `--progress <MODE>` | Progress reporting, independent of the log level: `auto` (default), `bar`, `plain`, `none` |
| `--quiet` | Silence progress and the summary; warnings and errors still print. Conflicts with `-v` and `--progress` |

`--trim-quality`, `--best-quality-segment`, and `--split-quality` are alternative
strategies for one stage, so at most one is accepted. `-H`/`-T` combine with any
of them.

An adapter source is `--adapter-fasta`, `--adapter-preset`, or
`--discover-adapters`. `--adapter-error-rate`, `--adapter-end-search`, and
`--adapter-sample-reads` are rejected without one. A forced `--input-format` that
disagrees with the stream, an output extension that names no format, and
FASTQ-to-BAM output are reported before any output is written. Adapter trimming is described in
[adapters.md](adapters.md).

## Quality filtering and trimming

`--quality-mode` determines the segment-level score used by `--min-quality`
and `--max-quality`. It does not affect trimming, best-segment selection, or
split locations.

| Mode | Calculation |
|---|---|
| `mean` (default) | Average per-base error probabilities, then convert the average to a Phred score |
| `arithmetic` | Average the numerical Phred scores directly |
| `median` | Take the median Phred score |

`--trim-quality PHRED` compares individual base scores with PHRED, removing
bases from each end until a base meets the threshold. `--split-quality PHRED`
removes stretches of consecutive bases below PHRED when they reach
`--split-min-low-quality-bases BASES`. Shorter internal stretches remain in
their segment; low-quality bases at segment ends are trimmed.

`--best-quality-segment PHRED` maximizes the cumulative score
`10^(-PHRED/10) - 10^(-base_quality/10)` over a contiguous segment. This
modified Mott calculation can retain bases below PHRED. Each segment from
adapter processing is evaluated separately, so an original read can still
produce multiple outputs.

```bash
whittle -i reads.fastq.gz -o trimmed.fastq.gz \
  --split-quality 9 --split-min-low-quality-bases 50 \
  --min-quality 12 --quality-mode median
```

The command splits at at least 50 consecutive bases below Q9, then keeps
segments with median quality of at least Q12. Length and GC filters also apply
to each output segment.

## Parameter aliases

`--head-crop` is an alias for `--trim-front`, and `--tail-crop` is an alias for
`--trim-tail`. Both appear in `--help`, accept a base count, and retain the
short options `-H` and `-T`. Other parameters use the names listed above.

## Stage order

Adapter preparation loads a FASTA or preset, or discovers adapters from a
sample of original reads. Sampled reads remain in the processing stream.
Each read then passes through these stages:

1. Adapter trimming and interior splitting on the original sequence, including
   terminal cleanup of the resulting segments.
2. Intersection of each segment with the retained barcode interval from the
   original `bi` tag, when `--trim-barcodes` is enabled.
3. Fixed cropping at both ends of each retained segment.
4. Quality end trimming, best-segment selection, or quality splitting.
5. Length, quality, and GC filtering of each final segment.
6. Tag reconstruction and output of surviving segments.

An unmatched read enters the barcode and crop stages as one full-length
segment. Cropping applies once per adapter-derived segment. Quality splitting
can divide that cropped segment again; its pieces are not cropped again.
Adapter matching uses the original sequence before barcode restriction and
fixed cropping.

Final intervals are collected in original-read order before filtering. FASTQ
and ONT names use one suffix per final interval, such as `read_segment_1`
through `read_segment_4` for two adapter segments each split into two quality
segments. Numbering does not restart at adapter boundaries or append nested
suffixes. Filtering preserves these indices: if only interval 3 passes, its
name remains `read_segment_3`. A read producing only one final interval keeps
its original name. PacBio names use final query coordinates according to the
[platform rules](tags.md#platform-rules).

## Barcode trimming

`--trim-barcodes` removes the barcode spans that dorado recorded in the `bi` aux
tag. The positions are read from the original tag and define the retained
interval. Each adapter-derived segment is intersected with that interval
before cropping. The trim uses the same tag-rewrite machinery as every other
stage: `MM`/`ML`/`MN`, per-base kinetics, and the ONT move table are rewritten
for the trimmed sequence.

```bash
whittle -i barcoded.bam -o trimmed.bam --trim-barcodes --update-signal-tags
```

The tag holds seven floats, four of which are positions: the front barcode's
start and length and the rear barcode's end and length. A barcode that dorado
did not find is stored as a negative position and leaves that end untouched, so
a read barcoded at one end is trimmed at that end only. A record without `bi`
passes through unchanged.

`bi` is dropped from a trimmed read, since its positions index the untrimmed
sequence; the barcode call (`BC`, `bv`) is a per-read label and is kept. A `bi`
that is not a seven-element float array, or whose positions describe an empty,
inverted, or out-of-range window, leaves the read untrimmed and is counted under
`warnings.barcode_tag_malformed_reads`.

The flag requires BAM or tagged FASTQ input carrying auxiliary tags.

## Tag removal

`--remove-tag <TAG>` removes one two-character aux tag from every output record.
It is repeatable, and each value must be exactly two alphanumeric characters,
validated before the run starts.

`--remove-kinetics` removes the per-base kinetics and alignment-count arrays in
one flag: `ip`, `pw`, `fi`, `fp`, `ri`, `rp`, `sa`, `sm`, and `sx`. It is
equivalent to naming each with `--remove-tag`, and the two flags fill one
removal set.

```bash
whittle -i reads.bam -o smaller.bam --remove-kinetics --remove-tag ML
```

Removal runs after the tags kept in register have been rewritten, so removing
one of them leaves the rest of the record intact: removing `MM` keeps the
rebuilt `ML` and `MN`, and removing one per-base array still slices the others
to the trimmed window.

Removal applies on every BAM output path. A run that trims nothing writes
records back without decoding them; with tag removal each record is rebuilt
instead, which costs the decode and changes nothing else. On BAM-to-FASTQ the
removal applies to the header tags selected by `--fastq-tags`.

Both flags require BAM or tagged FASTQ input. Tag removal is a complete run on its own:
`whittle -i in.bam -o out.bam --remove-kinetics` with no trimming options is
valid. The resolved set is recorded under `params.remove_tags` in the summary
JSON, and `params.strip_kinetics` records which flag requested it.

## Summary JSON

`--summary-json <PATH>` writes one JSON object describing the run: the resolved
settings under `params` and the counters under `reads`, `bases`, and
`segments_dropped`. It is written on every dispatch path, including directory
merges, regardless of `--quiet` or the log level. A write failure fails the run,
so a stale file from an earlier invocation is never left in place.

```bash
whittle -i reads.bam -o trimmed.fastq.gz -l 500 --quiet --summary-json qc.json
```

```json
{
  "schema_version": 1,
  "tool": "whittle",
  "version": "0.2.0",
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

`params` is abbreviated above; the file carries every resolved setting,
including defaults. `params.ordered` records whether a multithreaded run wrote
records in input order.

Under `warnings`, `malformed_tag_reads` counts reads whose per-base tag length
disagreed with the sequence and was left untouched, `malformed_mod_reads` counts
reads whose `MM`/`ML`/`MN` block could not be parsed and was removed from the
output, and `barcode_tag_malformed_reads` counts reads whose `bi` positions did
not describe a window inside the read under `--trim-barcodes`. All three are
also reported on stderr at the end of the run.

`reads.output` counts output segments, not input reads, so under `--split-quality`
or chimera splitting it can exceed `reads.input`. The read-level buckets
`with_output`, `trimmed_to_nothing`, and `all_filtered` partition `reads.input`.

Under `params.adapters`, `configured` is the set requested (preset and/or FASTA)
and `count` is the set searched, after presence detection narrowed it or
inference replaced it. The two are equal when neither ran. The startup banner
prints `configured`. Under `--discover-adapters` nothing is configured up front, so
`configured` is `0`.

`schema_version` is incremented only when an existing field changes meaning or
disappears. New fields may appear without an increment, so consumers should
tolerate unknown fields.

## Logging and progress

The log level is set with `-v`/`-vv` or `--quiet` (warnings and errors only).
All logging goes to stderr; stdout carries only read data.

`-v` adds stage detail: the detected input format and detection time, the
thread budget, and the read and base counts at the end. `-vv` adds per-read
decisions, each attributed to the read that produced it:

```text
[2026-09-01 12:13:45] [TRACE] [read{name=read_adapter}] Adapter hit adapter="LSK109_front" start=0 end=28 cost=0 action="trim 5'"
[2026-09-01 12:13:45] [TRACE] [read{name=read_adapter}] Segment kept segment=1 of=1 start=28 end=217 len=189
[2026-09-01 12:13:45] [TRACE] [read{name=read_short}] Segment dropped segment=1 of=1 start=0 end=30 len=30 reason="too short"
[2026-09-01 12:13:45] [TRACE] [read{name=read_short}] Every segment filtered produced=1
```

Every line has the same shape: timestamp, level, the enclosing span with its
identifying fields, a capitalized message, then structured fields. The message
is prose and the field values are data, so a line can be filtered on either.
The per-read lines record which adapter matched, over what span, at what cost,
and the resulting action, or why a segment was dropped.

`WHITTLE_LOG` overrides the level with a `RUST_LOG`-style filter, for example
`WHITTLE_LOG=whittle::adapter=trace` for adapter decisions without the
per-segment lines. `--quiet` takes precedence over it. A value that does not
parse is reported as a warning and the level falls back to the `-v`/`-vv`
setting.

Progress is shown as a live bar when stderr is a terminal and as periodic lines
(every 30 s, or 10 s under `-v`) when stderr is a file or a pipe. The bar is
never written to a non-terminal, so a redirected log contains no escape
sequences or carriage returns. `--progress` selects the mode independently of
the log level:

| Value | Behavior |
|---|---|
| `auto` | A bar on a terminal, periodic lines otherwise (default) |
| `bar` | The bar even when redirected; falls back to periodic lines under `-v`/`-vv` or `WHITTLE_LOG`, since debug lines and the bar cannot share a terminal |
| `plain` | Periodic lines, never a bar |
| `none` | No progress reporting; the banner, warnings, and summary still print |

`--quiet` also drops the summary and conflicts with `--progress`.

## Man page

Release tarballs ship `man/whittle.1` next to the binary, and the same file is
in the repository:

```bash
install -Dm644 man/whittle.1 /usr/share/man/man1/whittle.1
```
