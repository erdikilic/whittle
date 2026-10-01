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
`TAG:TYPE:VALUE` fields) is tagged FASTQ. Each header is inspected independently,
so tagged and plain reads can share a file or directory. Tagged records are
decoded as SAM aux tags and rewritten per output segment exactly as a uBAM
record's are ([tags.md](tags.md)), `--fastq-tags` selects the output tags, and
barcode positions (`bi`) and `--remove-tag` apply. A field that
does not parse as a SAM tag fails the run and names the read. A header whose
tab-delimited text is not in this form is copied verbatim.

### Threads

`-t`/`--threads` sets the worker count (default: every detected CPU) and bounds
the threads that do work. Trimming, record serialization, and output
compression run on a pool; BGZF input (BAM and FASTQ.bgz) is decompressed by one
quarter of the workers, at least one, and the pool takes the rest. Plain gzip
input is a single DEFLATE stream and is decompressed by one of the workers on
its own thread, ahead of the parser, so a gzip FASTQ run stops gaining from
more workers once that thread is saturated; BGZF
FASTQ input (`bgzip`, or whittle's own `.gz` output) has no such limit. The
startup banner reports the split.

Records are written in completion order under `-t > 1`. `--preserve-order`
restores the input order using a bounded window of batches. A slow batch limits
read-ahead until it is written.

## Directory input

`-i` accepts a directory. Every read file directly inside it is merged into one
output in natural filename order (digit runs compared numerically). The files
must share one format family (all FASTQ or all BAM); hidden files and
subdirectories are ignored; a mixed or empty directory is an error. An output
path inside the input directory is rejected, since it cannot be distinguished
from an input on a later run.

BAM output requires identical read-group definitions across input headers.
Conflicting definitions or different read-group sets fail the run. Only the
first header is written; FASTQ output does not require matching read groups.

```bash
whittle -i fastq_pass/barcode03/ -o barcode03.trimmed.fastq.gz --quality-trim ends --quality-cutoff 10
```

## Options

| Flag | Meaning |
|---|---|
| `-h, --help` | Print the option list and examples |
| `--version` | Print the version and exit |
| `-i, --input <PATH>` | Input file or directory (omit, or pass `-`, for stdin) |
| `-o, --output <PATH>` | Output file (omit, or pass `-`, for stdout); under `--split-by`, a path holding `{target}`, `{group}` or `{barcode}` is a template with one file per split key |
| `--input-format`, `--output-format <FORMAT>` | Force a format instead of detecting it: `fastq`, `fastq-gz`, `fastq-bgz`, `bam` |
| `--fastq-tags <all\|none\|TAGS>` | Aux tags written into FASTQ headers on BAM or tagged FASTQ input (default `all`) |
| `-c, --compression-level <0-9>` | BGZF level for `.gz`, `.bgz` and BAM output (default 4 for `.gz`, 6 for `.bgz` and BAM); ignored for plain FASTQ |
| `--summary-json <PATH>` | Write a machine-readable run summary to PATH; ignored under `--adapter-report` |
| `--rejected-output <PATH>` | Write every input read or trimmed segment that does not reach the output to PATH, in the output's format family, with a `wr:Z` tag naming the reason |
| `-t, --threads <N>` | Worker threads, at least 1 (default: all detected CPUs, clamped to that maximum) |
| `--preserve-order` | Write records in input order under `-t > 1` |
| `-l, --min-length <BASES>` | Minimum length to keep, per output segment (default 1) |
| `-L, --max-length <BASES>` | Maximum length to keep |
| `-q, --min-quality <PHRED>` | Minimum post-trim segment quality, a finite value of at least 0 (default 0) |
| `-Q, --max-quality <PHRED>` | Maximum post-trim segment quality, a finite value of at least 0 (default 1000) |
| `--max-expected-errors <E>` | Maximum expected errors per output segment, the sum of its per-base error probabilities; a finite value of at least 0 |
| `-g, --min-gc <FRACTION>`, `-G, --max-gc <FRACTION>` | GC-fraction bounds (0 to 1; `0.4` means 40%) |
| `-m, --quality-mode <MODE>` | Quality calculation of the mean-quality filter, `--min-quality`/`--max-quality`, only: `mean` (mean error probability as a Phred score, the default), `arithmetic` (mean of the Phred scores), `median` |
| `--tag-filter <EXPR>` | Keep only reads whose aux tags satisfy EXPR (samtools `-e` syntax over `[tag]` values); repeatable, every expression must hold; applied before adapter discovery and trimming (BAM or tagged FASTQ input) |
| `-H, --trim-front <BASES>`, `-T, --trim-tail <BASES>` | Fixed crop from each adapter-derived segment after barcode restriction and before quality processing; applied once |
| `--quality-trim <METHOD>` | Quality trimming of each adapter-derived segment: `ends`, `best`, `segments` or `runs` ([below](#quality-filtering-and-trimming)); requires `--quality-cutoff` |
| `--quality-cutoff <PHRED>` | Phred cutoff of `--quality-trim`; required with it and rejected without it |
| `--min-low-quality-run <BASES>` | Consecutive bases below the cutoff that split a read under `--quality-trim runs`; shorter runs stay inside their piece (default 1) |
| `--quality-end-cutoff <PHRED>` | Trim the ends of every piece of `--quality-trim segments` or `runs` up to the first base at or above PHRED (default `--quality-cutoff`) |
| `--update-moves` | Rewrite ONT signal tags through trimming instead of removing them (BAM-to-BAM; requires DNA or RNA model metadata in the read-group description) |
| `--remove-tag <TAGS>` | Remove aux tags from every output record; comma-separated and repeatable; an item is a two-character tag or a group: `kinetics` (`ip pw fi fp ri rp sa sm sx`), `mods` (`MM ML MN`), `signal` (`mv ts ns sp pi`) (BAM or tagged FASTQ input) |
| `-a, --adapter-fasta <FILE>` | Adapter and primer FASTA (IUPAC codes accepted; `primer` or `barcode` in a header description restricts the entry to the read ends); enables adapter trimming |
| `--adapter-preset <KITS>` | Built-in kit presets, comma-separated: `lsk114`, `rad114` (`ulk114`), `rbk114`, `nbd114`, `pcb114` (`pcs114`), `rpb114`, `mab114`, `rna004`, `pacbio`, `ont`, `all`; enables adapter trimming |
| `--adapter-error-rate <FRACTION>` | End-match tolerance as a fraction of adapter length (default 0.2); requires an adapter source |
| `--adapter-end-search <BASES>` | End-zone width searched for terminal adapters (default 150); requires an adapter source |
| `--adapter-ends-only` | Trim adapters at ends only; disable interior adapter splitting independently of quality splitting |
| `--adapter-sample-reads <COUNT>` | Reads inspected for preset presence or adapter discovery (defaults 2000 and 40000; at least 100); `0` disables preset detection and is rejected for discovery; ignored with `--adapter-fasta` unless discovering adapters |
| `--adapter-discover` | Discover adapters, barcodes and primers from the sampled reads and trim them; a preset or FASTA is trimmed first and discovery continues beyond it; enables adapter trimming |
| `--adapter-report` | Run the same discovery, print the discovered FASTA to stdout and exit without read output or a JSON summary; conflicts with `--adapter-discover` |
| `--split-by <SPEC>` (repeatable) | Assign every output read or segment to a primer target and tag it `wt:Z`; each `SPEC` adds a source, merged into one sheet: a sheet file (TSV or annotated FASTA), the preset `mab114`, or an inline target ([below](#primer-split)); a target holds one or more forward and reverse primers |
| `--split-require <RULE>` | Which end(s) need a located primer for an assignment: `either` (default), `both`, `fwd`, `rev` |
| `--split-action <trim\|retain>` | Whether a located split primer is trimmed (default) or kept by widening the piece back over it |
| `--split-lead <N>` | Minimum edit-cost lead of the best target over the best different target; a smaller lead is ambiguous (default 2) |
| `--split-discard <BINS>` | Drop these bins instead of writing them, comma-separated: `unassigned`, `ambiguous` |
| `-v, --verbose` (repeatable) | Stage detail with `-v`, per-read decisions with `-vv` |
| `--progress <MODE>` | Progress reporting, independent of the log level: `auto` (default), `bar`, `plain`, `none` |
| `--quiet` | Silence progress and the summary; warnings and errors still print. Conflicts with `-v` and `--progress` |

`--quality-trim` selects one method for the quality stage. `-H`/`-T` combine
with any of them.

An adapter source is `--adapter-fasta`, `--adapter-preset`, or
`--adapter-discover`. `--adapter-error-rate`, `--adapter-end-search`, and
`--adapter-sample-reads` are rejected without one. A forced `--input-format`
that disagrees with the stream, an output extension that names no format, and
FASTQ-to-BAM output are reported before any output is written. Adapter
trimming is described in [adapters.md](adapters.md).

Adapter sampling stops at the requested read count, 256 MiB of retained payload,
or 64 Mi bases, whichever is reached first. The last record is kept whole and
can exceed a limit. The log reports the sampled count and any payload limit;
sampled records remain in the processing stream.

## Quality filtering and trimming

### Mean quality

`--quality-mode` determines the segment-level score used by `--min-quality`
and `--max-quality`. It applies to this mean-quality filter only. The
`--quality-trim` methods and `--max-expected-errors` work on per-base
qualities or error probabilities and are not affected, and a recomputed
dorado `qs` tag is always the error-probability mean over the kept segment.

| Mode | Calculation |
|---|---|
| `mean` (default) | Average per-base error probabilities, then convert the average to a Phred score |
| `arithmetic` | Average the numerical Phred scores directly |
| `median` | Take the median Phred score |

The default `mean` averages the per-base error probabilities and converts the
average back to a Phred score, as dorado, chopper, NanoFilt and Filtlong do.
`--min-quality Q` in this mode keeps a segment whose expected error rate, its
expected errors divided by its length, is at most `10^(-Q/10)`: the test of
`vsearch --fastq_maxee_rate`. `arithmetic` averages the Phred values, as fastp
and fastplong do. It weighs each base by its score rather than by its error
probability, so on reads with many bases at the quality cap it gives much
higher values than `mean`, and the same `--min-quality` keeps more reads.

dorado computes its `qs` tag over the bases after the first 60, so
`--min-quality` and a read's `qs` can differ, most on reads whose first bases
are of low quality. `--tag-filter '[qs]>=10'` filters on dorado's own value,
before any trimming. A trimmed record's `qs` is recomputed over the whole kept
segment, the first 60 bases included.

On PacBio HiFi, base QVs are binned on current instruments and the read
accuracy is in the `rq` tag. Filter HiFi reads with
`--tag-filter '[rq]>=0.99'` rather than with a mean recomputed from binned
QVs. Quality trimming is rarely needed for HiFi reads.

### Expected errors

`--max-expected-errors E` rejects a segment whose expected errors, the sum of
its per-base error probabilities `10^(-Q/10)` (Edgar and Flyvbjerg 2015),
exceed E. It is computed on the final segment, after all trimming and
splitting, so each piece of a split read is judged on its own, and it is
reported as `expected_errors`. The mean-quality bound and this bound differ in
what they hold fixed: `--min-quality Q` in `mean` mode limits the expected
error rate, the expected errors divided by the length, to `10^(-Q/10)`
whatever the length, while `--max-expected-errors` limits the total per
segment, so a longer segment needs a higher per-base accuracy to pass. The
total suits reads of a fixed length, such as full-length 16S or ITS
amplicons, targeted panels and DADA2-style workflows. On genomic long reads,
whose lengths vary by orders of magnitude, it would reject long reads for
their length alone; `--min-quality` is the bound for them.

### Quality trimming

`--quality-trim METHOD` trims each segment from adapter processing, after the
fixed crop, against `--quality-cutoff PHRED`. The cutoff has no default: the
right value depends on the platform and basecaller, so it is always given.

| Method | Keeps |
|---|---|
| `ends` | The segment left after removing bases from each end up to the first base at or above PHRED |
| `best` | The contiguous segment with the highest modified Mott score; bases below PHRED can be kept |
| `segments` | Every maximal scoring segment under the modified Mott score |
| `runs` | The pieces between runs of at least `--min-low-quality-run` consecutive bases below PHRED |

`ends` compares individual base scores with PHRED, removing bases from each end
until a base meets the threshold. It stops at the first base at or above
PHRED, so a single good base inside a low-quality tail ends the trimming
there. `best` and `segments` score the whole tail, so an isolated good base
does not stop them. `runs` removes stretches of consecutive bases
below PHRED when they reach `--min-low-quality-run BASES`. Shorter internal
stretches remain in their piece; low-quality bases at piece ends are trimmed.

`best` maximizes the cumulative score `10^(-PHRED/10) - 10^(-base_quality/10)`
over a contiguous segment (modified Mott, as in phred). It can retain bases
below PHRED. Each segment from adapter processing is evaluated separately, so
an original read can still produce multiple outputs.

`segments` keeps every maximal scoring segment under the same score: the
maximal scoring subsequences of Ruzzo and Tompa (1999), computed in linear
time. The best segment is one of them; the others are the best segments of
the parts on either side, found recursively. Two high-quality regions stay
separate when the bases between them cost more than the smaller of their two
scores, so a read with a low-quality interior keeps both flanks. The split
follows the score rather than a run of consecutive low-quality bases, so it
also separates regions that mix bases below and above PHRED. Equal scores
favor the longer segment, as in `best`, and a read with one region above
PHRED yields the same segment under both methods when that segment reaches
the score floor.

A segment is kept when its score is at least that of 50 error-free bases,
`50 * 10^(-PHRED/10)`. An error-free base contributes the cutoff error
probability, the largest score one base can add, and a base at the cutoff
contributes nothing, so the floor is counted in error-free bases rather than in
bases at the cutoff. Counted this way, the floor needs the same bases at any
PHRED: 51 bases 20 Phred units above the cutoff, or 101 bases 3 units above
it, and never fewer than 51 bases. Short high-scoring stretches inside
low-quality regions fall below it and are removed with the low-quality bases.
`--min-length` then filters the kept segments like the pieces of `runs`; a
segment below it is reported as `too_short`.

`--quality-end-cutoff PHRED` separates the strictness of the split from that
of the ends under `segments` and `runs`: the method splits at
`--quality-cutoff`, then every piece, including a read that was not split,
has its ends trimmed as `ends` does at the end cutoff, up to the first base
at or above it. A piece with no such base is dropped, and `--min-length`
applies afterwards. A low `--quality-cutoff` with a higher end cutoff splits a
read only at long weak interior regions while trimming its ends more
strictly. The default equals `--quality-cutoff`, and a value at or below it
changes nothing, since every piece already starts and ends at a base at or
above `--quality-cutoff`. The end trim does not revisit the score floor of
`segments`.

```bash
whittle -i reads.fastq.gz -o trimmed.fastq.gz \
  --quality-trim runs --quality-cutoff 9 --min-low-quality-run 50 \
  --min-quality 12 --quality-mode median
```

### Choosing a quality operation

- **No quality trimming.** Most long-read uses need none: adapter trimming
  with a length and quality filter (`-l`, `-q`) keeps whole reads, which
  aligners and assemblers handle better than reads cut at local quality dips.
- **`best`** for reads with low-quality ends. It keeps one segment per
  adapter-derived segment and is not stopped by isolated good bases.
- **`segments`** for reads with low-quality interior regions. It keeps both
  flanks of a weak region, including regions that mix low and moderate bases,
  and `--quality-end-cutoff` trims the piece ends more strictly than the
  split.
- **`runs`** when a split must require a minimum run of consecutive bases
  below the cutoff, set with `--min-low-quality-run`.
- **`ends`** for plain end trimming to the first good base.
- **`--max-expected-errors`** for amplicons of a fixed length, next to or
  instead of `-q`.

The command splits at at least 50 consecutive bases below Q9, then keeps
segments with median quality of at least Q12. Length and GC filters also apply
to each output segment.

## Parameter aliases

`--head-crop` is an alias for `--trim-front`, and `--tail-crop` is an alias for
`--trim-tail`. Both appear in `--help`, accept a base count, and retain the
short options `-H` and `-T`. Other parameters use the names listed above.

## Filtering by aux tag

`--tag-filter <EXPR>` keeps only the reads whose aux tags satisfy an expression
and drops the rest before anything else happens: a rejected read is not
sampled for adapter discovery, not trimmed and not written. The flag is
repeatable and every expression must hold. It reads tags, so it requires BAM
or tagged FASTQ input.

The syntax is the aux-tag subset of `samtools view -e`: `[tag]` names a tag,
the comparisons are `==`, `!=`, `<`, `<=`, `>`, `>=`, the connectives are `!`,
`&&`, `||` and parentheses, string literals are double-quoted, and `[tag]`
alone (or `exists([tag])`) tests presence.

```bash
# Drop dorado's adaptive-sampling rejects.
whittle -i reads.bam -o kept.bam --tag-filter '[er]!="data_service_unblock_mux_change"'
# Duplex reads, or simplex reads with a dorado Q-score of at least 15.
whittle -i reads.bam -o kept.bam --tag-filter '[dx]==1 || ([dx]==0 && [qs]>=15)'
# HiFi reads by predicted accuracy and passes.
whittle -i hifi.bam -o kept.bam --tag-filter '[rq]>=0.99 && [np]>=3'
# One barcode, or unclassified reads.
whittle -i reads.bam -o kept.bam --tag-filter '[BC]=="barcode03" || ![BC]'
```

Numbers compare numerically whatever the tag's integer width or float type;
strings compare bytewise, so ISO 8601 timestamps such as `st` order correctly.
A comparison with a tag the read does not carry is false whatever the
operator, as in samtools, so `![tag]` or a negated comparison is the way to
accept absence. Comparing a number with a string, or comparing an array, is an
error that names the read and the tag, so a typo does not silently select or
drop every read. Rejected reads are counted as input and reported as `Tag
filtered` on stderr and as `reads.tag_filtered` in the summary JSON.

## Stage order

Adapter preparation loads a FASTA or preset, or discovers adapters from a
sample of original reads. Reads rejected by `--tag-filter` never reach the
sample. Sampled reads remain in the processing stream. Each read then passes
through these stages:

1. Adapter trimming and interior splitting on the original sequence, including
   terminal cleanup of the resulting segments.
2. Intersection of each segment with the verified barcode spans from the
   original `bi` tag, when an adapter source is given.
3. Fixed cropping at both ends of each retained segment.
4. Quality trimming by the `--quality-trim` method.
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

## Barcode positions

With an adapter source, whittle reads the barcode spans dorado recorded in
the `bi` aux tag and removes each span at which a barcode sequence is found:
a barcode of the configured set, or the catalog barcode named by the read's
`BC` call. A span holding no barcode sequence, such as stale positions on
input dorado has already trimmed, is left alone and counted under
`warnings.barcode_tag_unverified_reads`. There is no flag; the positions are
used whenever the input carries them. Each adapter-derived segment is
intersected with the retained interval before cropping, and the trim uses
the same tag-rewrite machinery as every other stage: `MM`/`ML`/`MN`, per-base
kinetics, and the ONT move table are rewritten for the trimmed sequence.

```bash
whittle -i barcoded.bam -o trimmed.bam --adapter-preset nbd114 --update-moves
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

Positions are read from BAM and tagged FASTQ input; plain FASTQ carries none.

## Tag removal

`--remove-tag <TAGS>` removes aux tags from every output record. The value is
a comma-separated list and the flag is repeatable. An item is either a
two-character alphanumeric tag or one of three group names, validated before
the run starts:

| group | tags |
|---|---|
| `kinetics` | `ip`, `pw`, `fi`, `fp`, `ri`, `rp`, `sa`, `sm`, `sx` (PacBio per-base kinetics and alignment counts) |
| `mods` | `MM`, `ML`, `MN` (base modifications) |
| `signal` | `mv`, `ts`, `ns`, `sp`, `pi` (ONT signal mapping) |

```bash
whittle -i reads.bam -o smaller.bam --remove-tag kinetics,ML
```

Removal runs after the tags kept in register have been rewritten, so removing
one of them leaves the rest of the record intact: removing `MM` keeps the
rebuilt `ML` and `MN`, and removing one per-base array still slices the others
to the trimmed window.

Removal applies on every BAM output path. A run that trims nothing writes
records back without decoding them; with tag removal each record is rebuilt
instead, which costs the decode and changes nothing else. On BAM-to-FASTQ the
removal applies to the header tags selected by `--fastq-tags`.

The flag requires BAM or tagged FASTQ input. Tag removal is a complete run on its own:
`whittle -i in.bam -o out.bam --remove-tag kinetics` with no trimming options is
valid. The resolved set, with groups expanded, is recorded under
`params.remove_tags` in the summary JSON.

## Rejected output

`--rejected-output <PATH>` writes everything that does not reach the main output
to a second file: reads rejected by `--tag-filter` (as read, untrimmed),
reads that trimming left without a segment (as read), and every trimmed
segment a filter dropped (as trimmed, with its `_segment_N` name and its
rewritten tags). Each record carries a `wr:Z` tag with one of `tag_filter`,
`trimmed_to_nothing`, `too_short`, `too_long`, `low_quality`, `high_quality`,
`expected_errors` or `gc`. In FASTQ output the tag is a header field, as for tagged FASTQ.

The file takes the output's format family, BAM for BAM output and FASTQ for
FASTQ output, with compression from its own extension (`.bam`, `.fastq`,
`.fastq.gz`, `.fastq.bgz`), and is written by one extra thread. Counts in the
summary are unchanged by the flag.

```bash
whittle -i reads.bam -o kept.bam --rejected-output rejected.bam -l 500 -q 10 --tag-filter '[dx]==1'
samtools view rejected.bam | cut -f 1,10 | head
```

## Primer split

`--split-by SPEC` assigns every output read, or every chimera-split
segment, to a target defined by a primer pair, a primer mix, or a pool of
either, and tags it `wt:Z` with the call. A target holds a list of forward
primers and a list of reverse primers; any forward primer pairs with any
reverse primer, as in a kit that ships several primer variants in one tube.
The flag is repeatable; each `SPEC` adds a source,
and every source merges into one sheet. `SPEC` resolves in this order: an
existing file path (a TSV or annotated FASTA sheet), the preset token
`mab114`, or, when neither matches, the inline form (below). Split primers
are located and trimmed alongside the other primers of the adapter search:
at each end, the outermost sheet primer in the orientation valid for that
end is located and trimmed, and a read is split at a sheet primer only
where two of them form a chimera junction (or one lies beside an interior
adapter or barcode), never at a lone primer inside it, so a nested or
overlapping panel keeps full-length amplicons whole. `--split-by` does not
force amplicon mode on.
[adapters.md](adapters.md#primer-split) describes how they join the engine.

### Sheet formats

A TSV sheet has a header row and requires the columns `target`, `fwd` and
`rev`; `group`, `min_len` and `max_len` are optional.

```
target	fwd	rev	group
16S_full	AGRGTTYGATYMTGGCTCAG	CGGTTACCTTGTTACGACTT	16S
16S_V34	CCTACGGGNGGCWGCAG	GACTACHVGGGTATCTAATCC	16S
ITS	CTTGGTCATTTAGAGGAAGTAA	TCCTCCGCTTATTGATATGC	ITS
```

`rev` is written 5' to 3', as ordered; whittle derives its reverse
complement. A `fwd` or `rev` cell holds one sequence or a comma-separated
list of them, the primer mix of that role:

```
target	fwd	rev
ITS	TCCGTAGGTGAACCTGCGG,TCCGTTGGTGAACCAGCGG,TCTGTAGGTGAACCTGCAG	TCCTCCGCTTATTGATATGC,TCCTCCGCTTATTAATATGC
```

Either primer cell may be empty, for a target with no primer of that role;
an empty sequence inside a list (a trailing or doubled comma) is an error.
A cell of one primer names it `<target>_fwd` or `<target>_rev`, a list
`<target>_fwd1`, `<target>_fwd2`, and so on.
`group` defaults to the target name and is the key an `-o` template names
under `{group}`. A primer sequence may appear in more than one row; whittle
stores it once and shares it among the targets that name it. `min_len` and
`max_len` bound the final segment length, inclusive; a segment outside the
window is unassigned (`length`).

An annotated FASTA extends the existing adapter-FASTA header parse with
`key=value` fields:

```
>27F primer target=16S end=fwd group=16S
AGRGTTYGATYMTGGCTCAG
>1492R primer target=16S end=rev group=16S
TACGGYTACCTTGTTACGACTT
```

`target=` marks an entry as a split primer and names its target; `end=fwd`
or `end=rev` places it. Several entries may share one `target=` and `end=`;
each adds a primer to that list. `group=` overrides the default group and
may appear on any entry of a target; two entries of one target giving
different `group=` values is an error. Entries with no `target=` field keep
their existing meaning as `-a` entries.

The preset `mab114` defines two targets, `16S` and `ITS`, holding the primer
mixes of the ONT Microbial Amplicon Barcoding kit (SQK-MAB114.24) under
their ONT names:

| Target | Role | Primer | Sequence |
|---|---|---|---|
| `16S` | forward | 16S_mix_F | `AGRGTTYGATYMTGGCTCAG` |
| `16S` | forward | 16S_Bor_F | `AGAGTTTGATCCTGGCTTAG` |
| `16S` | forward | 16S_Chl_F | `AGAATTTGATCTTRGTTCAG` |
| `16S` | forward | 16S_Ent_F | `AGAGTTTGATCATGGCTCAG` |
| `16S` | reverse | 16S_mix_R | `SGGYTACCTTGTTACGACTT` |
| `16S` | reverse | 16S_Bor_R | `CGGCTACCTTGTTACGACTT` |
| `16S` | reverse | 16S_Chl_R | `GGGCTACCTTGTTACGACTT` |
| `ITS` | forward | ITS1 | `TCCGTAGGTGAACCTGCGG` |
| `ITS` | forward | ITS1_Fus | `TCCGTTGGTGAACCAGCGG` |
| `ITS` | forward | ITS1_Mal | `TCTGTAGGTGAACCTGCAG` |
| `ITS` | reverse | ITS4 | `TCCTCCGCTTATTGATATGC` |
| `ITS` | reverse | ITS4_Pyt | `TCCTCCGCTTATTAATATGC` |

These are the primers `--adapter-preset mab114` searches; the split preset
takes their sequences from the same catalog entries. The token also enables
`--adapter-preset mab114` unless `-a` or `--adapter-preset` is given. A
sheet file enables no library adapters by itself.

Target and group names are checked at load time: `unassigned`, `ambiguous`
and `unclassified` are reserved, since the classifier writes them itself,
every target name must be unique, and every name must be printable ASCII
(space included), since it fills the `wt:Z` tag. Each primer must be at least 11 nt of
nucleotide or IUPAC alphabet. A primer sequence may appear in several
targets. A target with no primer of a role that
`--split-require` needs (no `rev` under `both`, say) is an error. whittle
also warns for any two primers of different targets whose IUPAC edit
distance is below `--split-lead`: those targets cannot be reliably told
apart at that lead. The primers of one target are never compared with each
other.

A target may use one sequence for both roles, which splits reads by a
barcode-like sequence flanking both ends: the 5' end of the read holds the
sequence as given and the 3' end its reverse complement, on either strand.
Such a target is assigned from one end or both under `--split-require
either`, and only from both under `both`. Its two strands look alike, so
its reads carry no orientation information: they are all counted under
`plus`, and the `plus`/`minus` counts of its key say nothing about strand.

### Inline targets

When a `--split-by` value is neither an existing file nor a preset token,
it is read as an inline target:

```
NAME:F:SEQ[,SEQ...]:R:SEQ[,SEQ...][,F:SEQ[,SEQ...]:R:SEQ[,SEQ...]...]
```

`NAME` is the target name. `F:` and `R:` each introduce the forward or
reverse primer list: one sequence, or several separated by commas. A pair
may give one tag or both, each once, but not neither. `SEQ` follows the
sheet's own normalization and checks: uppercase, `U` folded to `T`, IUPAC
alphabet, at least 11 nt. `R` sequences are written 5' to 3', as ordered,
like the sheet's `rev`.

After a comma, a field that is exactly the tag `F` or `R`, in either
case, followed by `:` begins a new primer pair; anything else is another
sequence of the current list. A sequence that starts with `F`-like bases or
holds the IUPAC code `R` is therefore never mistaken for a tag. An empty
sequence in a list (a trailing or doubled comma), a trailing `:`, and a tag
repeated within one pair are errors. So is a pool in which a pair gives
only one tag and the pair after it only the other, as `16S:F:A,R:B` does:
that is one pair written with `,R:` in place of `:R:`, and would otherwise
be read as two targets with one primer each. Such targets go in separate
`--split-by` values.

```bash
whittle -i pool.fastq.gz -o out/{target}.fastq.gz \
  --split-by '16S:F:AGRGTTYGATYMTGGCTCAG:R:CGGTTACCTTGTTACGACTT'
```

One pair gives target `NAME` in group `NAME`, with primers named `NAME_F`
and `NAME_R`; a list of several primers names them `NAME_F1`, `NAME_F2`,
... and `NAME_R1`, `NAME_R2`, ... by position. The primer mixes of the
MAB114 kit, written inline, are one target per mix:

```bash
--split-by '16S:F:AGRGTTYGATYMTGGCTCAG,AGAGTTTGATCCTGGCTTAG,AGAATTTGATCTTRGTTCAG,AGAGTTTGATCATGGCTCAG:R:SGGYTACCTTGTTACGACTT,CGGCTACCTTGTTACGACTT,GGGCTACCTTGTTACGACTT'
--split-by 'ITS:F:TCCGTAGGTGAACCTGCGG,TCCGTTGGTGAACCAGCGG,TCTGTAGGTGAACCTGCAG:R:TCCTCCGCTTATTGATATGC,TCCTCCGCTTATTAATATGC'
```

A read carrying any forward primer of `16S` and any of its reverse primers
is assigned to `16S`; the variants of one target never make a read
ambiguous.

Several pairs, each begun by `,F:` or `,R:`, form a pool of distinct
amplicons: targets `NAME.1`, `NAME.2`, ... in group `NAME`, with primers
`NAME.1_F`, `NAME.1_R`, `NAME.2_F`, and so on.

```bash
--split-by '16S:F:AGRGTTYGATYMTGGCTCAG:R:CGGTTACCTTGTTACGACTT,F:CCTACGGGNGGCWGCAG:R:GACTACHVGGGTATCTAATCC'
```

defines `16S.1` and `16S.2`, both in group `16S`, so an `out/{group}.fastq`
template collapses them into one file. Lists and pools combine:
`16S:F:a,b:R:x,F:c:R:y` gives `16S.1` with forward primers `a` and `b` and
reverse primer `x`, and `16S.2` with forward primer `c` and reverse primer
`y`.

### Merging sources

Every `--split-by` value loads independently and the results merge into one
sheet, in the order given. Target names must be unique across all sources;
a name repeated by two sources is an error. Primers are shared by sequence
across sources exactly as within one sheet, so a primer sequence common to
a sheet file and an inline target is stored once. Validation, the
`--split-lead` edit-distance warning, and the output-template name checks
all run on the merged sheet, not on each source alone.

### End rule, lead and action

`--split-require <either|both|fwd|rev>` (default `either`) sets which
end(s) of a segment need a located primer before it can be assigned; `fwd`
and `rev` mean any primer of the target's forward or reverse list. An end
that scores several primers of one list counts at the cheapest of them.
`--split-lead <N>` (default 2) is the minimum edit-cost lead the cheapest
consistent target needs over the next cheapest different target; a tie, or
a smaller lead, is `ambiguous` rather than a guess. `--split-action
<trim|retain>` (default `trim`) trims the located split primers like any
other primer; `retain` widens the kept span back over them, so the primer
sequence stays in the output.

### Classification

Each segment is classified independently, so a chimera split into several
pieces gives each piece its own call:

- **Assigned** to a target's key (the target name, or its group under a
  `{group}` template), when a strand-consistent target beats every other
  key by at least `--split-lead`.
- **Unassigned**, for one of four reasons: `no_primer` (no primer located at
  either end), `require` (a target matched, but not at the end(s)
  `--split-require` needs), `orientation` (a target's primer scored at both
  ends, so no strand is consistent), or `length` (the segment length falls
  outside the assigned target's `min_len`/`max_len` window).
- **Ambiguous**, when two keys tie for cheapest, or the cheapest does not
  beat the runner-up by `--split-lead`.

A target is matched at an end when a primer of the list its strand puts
there scored at that end. A primer of another key closer than
`--split-lead` to a scored primer counts as scored there at that cost plus
one, so near-identical primers of different keys are judged by the lead
rather than by which one cleared its budget; a target matched only through
such a stand-in is never assigned. An end that scores a primer of another
key and none of the target's own does not rule the target out when every
score there is within `--split-lead` of the penalty: the lowest budget
rescoring applied at that end plus one, the least any primer that did not
score there can cost. The missing primer then counts at the penalty, the
same for every target, and such a target competes only when no target is
matched at every scored end. A score that leads the penalty by
`--split-lead` or more rules out every target that does not list its
primer at that end, so a chimera with exact primers of two keys at its two
ends is assigned to neither. A penalised end does not count toward
`--split-require`.

### Output routing

Without a placeholder, `-o` writes every record to one file and tags it:
tag-only mode. A path holding `{target}`, `{group}` or `{barcode}` is a
template, and each key is written to the path it expands to. The key is the
target when `{target}` is present in the template, otherwise the group when
`{group}` is present; `{group}` may still appear alongside `{target}` to
nest the output by group, as in `out/{group}/{target}.fastq.gz`, without
changing which one names the key. Under `{group}`, ambiguity is judged
across groups only. `unassigned` and `ambiguous` calls substitute those
words for the key placeholder: `out/{barcode}.{target}.bam` gives
`out/BC03.unassigned.bam`. `{barcode}` takes the record's `BC:Z` barcode
call (BAM or tagged FASTQ input), or `unclassified` when the record carries
none; plain FASTQ carries no barcode call, so every record takes
`unclassified`.

A template needs a known output extension to fix the format. A target or
group name that is empty, `.` or `..`, or holds a path separator, is not one
path component and is refused at load time when the template uses it, as
are two target or group names the template uses that differ only in letter
case, since they would share one file on a case-insensitive file system.
Tag-only mode accepts any such name, since it never becomes a path. A
`BC:Z` barcode call under `{barcode}` must be one path component too; it is
checked per record as the record is routed, and a call that is not one path
component is an error. Placeholders may appear in directory and file components, and whittle
creates the directories. Each expanded path is
checked before its file is opened, against the input path, against
`--rejected-output`, and against every other key's path, so two keys
expanding to the same file are refused rather than one silently
overwriting the other. A run opening more than 512 files is warned once.
Files exist only for keys that received at least one record; the report
still lists every sheet key.

Filter drops still go to `--rejected-output` as usual. The `unassigned` and
`ambiguous` bins hold only segments that passed the length, quality and GC
filters. `--split-discard <unassigned,ambiguous>` drops those bins instead
of writing them, and counts the dropped segments as `discarded`.

### Tag and report

Every output record, in every mode, carries a `wt:Z` tag naming its call
(the key, `unassigned` or `ambiguous`); an existing `wt` tag in the input is
replaced. [tags.md](tags.md#split-target) has the tag's exact placement on
BAM and FASTQ output. A record written to `--rejected-output` carries no
`wt:Z` tag.

The end-of-run log prints a per-key table, and `--summary-json` carries the
same counts in an additive `split` object ([below](#summary-json)).

```bash
whittle -i pool.bam -o out/{barcode}.{target}.fastq.gz --split-by primers.tsv -t 8
whittle -i pool.fastq.gz -o trimmed.fastq.gz --split-by mab114
whittle -i pool.fastq.gz -o out/{group}.fastq.gz \
  --split-by mab114 --split-by 'V34:F:CCTACGGGNGGCWGCAG:R:GACTACHVGGGTATCTAATCC'
```

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
  "schema_version": 2,
  "tool": "whittle",
  "version": "0.2.0",
  "command": "whittle -i reads.bam -o trimmed.fastq.gz -l 500 --quiet --summary-json qc.json",
  "input": "reads.bam",
  "output": "trimmed.fastq.gz",
  "elapsed_seconds": 12.34,
  "params": { "threads": 8, "ordered": false, "min_length": 500, "qual_mode": "mean", "quality_trim": null,
              "adapters": { "configured": 120, "count": 4, "sample": 500, "infer": "off" } },
  "reads": { "input": 1000, "output": 950, "with_output": 940, "trimmed_to_nothing": 30, "all_filtered": 20, "tag_filtered": 10 },
  "bases": { "input": 10000000, "output": 9500000 },
  "segments_dropped": { "too_short": 12, "too_long": 0, "low_quality": 5, "high_quality": 0, "expected_errors": 0, "gc_out_of_range": 0 },
  "warnings": { "malformed_tag_reads": 0, "malformed_mod_reads": 0, "barcode_tag_malformed_reads": 0, "barcode_tag_unverified_reads": 0 }
}
```

`params` is abbreviated above; the file carries every resolved setting,
including defaults. `params.ordered` records whether a multithreaded run wrote
records in input order.

Under `warnings`, `malformed_tag_reads` counts reads whose per-base tag length
disagreed with the sequence and was left untouched, `malformed_mod_reads` counts
reads whose `MM`/`ML`/`MN` block could not be parsed and was removed from the
output, `barcode_tag_malformed_reads` counts reads whose `bi` positions did
not describe a window inside the read, and `barcode_tag_unverified_reads`
counts reads with a recorded barcode span at which no barcode sequence was
found. All four are also reported on stderr at the end of the run.

`reads.output` counts output segments, not input reads, so under
`--quality-trim segments` or `runs`, or chimera splitting, it can exceed
`reads.input`. The read-level buckets
`with_output`, `trimmed_to_nothing`, `all_filtered` and `tag_filtered` partition
`reads.input`; `params.tag_filter` lists the expressions as written.

Under `params.adapters`, `configured` is the set requested (preset and/or FASTA)
and `count` is the set searched, after presence detection narrowed it or
inference replaced it. The two are equal when neither ran. The startup banner
prints `configured`. Under `--adapter-discover` and `--adapter-report` nothing is configured up front, so
`configured` is `0`.

`schema_version` is incremented only when an existing field changes meaning or
disappears. New fields may appear without an increment, so consumers should
tolerate unknown fields.

Under `--split-by`, the summary carries an additive `split` object; it is
absent without `--split-by`.

```json
"split": {
  "spec": ["primers.tsv"], "require": "either", "action": "trim", "lead": 2,
  "keys": [
    { "key": "16S", "reads": 420, "bases": 168000, "both_ends": 400, "five_only": 15,
      "three_only": 5, "plus": 410, "minus": 10 }
  ],
  "unassigned": { "no_primer": 30, "require": 4, "orientation": 1, "length": 2 },
  "ambiguous": 6,
  "discarded": 0
}
```

`spec` is the `--split-by` values as given, in order; `require`, `action`
and `lead` are the resolved `--split-require`, `--split-action` and
`--split-lead` settings. `keys` lists every sheet key (one per target, or
one per group when the `-o` template holds `{group}`), in sheet order, zero
counts included; `both_ends`, `five_only` and `three_only` count how many
of a key's assigned segments had primer evidence at each end, and
`plus`/`minus` count its strand. A target that uses one sequence for both
roles has no strand to tell, so its segments all count under `plus`.
`unassigned` breaks unassigned segments down by reason. `discarded` counts
segments `--split-discard` dropped from the `unassigned` or `ambiguous` bins
rather than writing; a discarded segment is still counted under its key or
reason, so `sum(keys[].reads) + sum(unassigned.*) + ambiguous` always equals
`reads.output + discarded`. The same per-key counts are logged as a table at
the end of the run.

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
