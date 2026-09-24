# Position-indexed tags

Long-read records carry auxiliary tags that are indexed by base position.
Trimming the sequence without rewriting them leaves them referring to bases
that no longer exist. whittle rewrites every such tag on BAM-to-BAM and
BAM-to-FASTQ output, and on FASTQ input whose headers carry the tags in the
`samtools fastq -T` convention (tagged FASTQ, [cli.md](cli.md#input-and-output)),
which takes the BAM-to-FASTQ path.

## Base-modification tags

Modification calls are stored per read in `MM`, which lists modified bases as
skip counts over the sequence, and `ML`, which holds their probabilities.
Trimming `SEQ`/`QUAL` alone leaves a block whose skip counts index the original
sequence.

For every output uBAM read, whether cropped, quality-trimmed, or split, `MM` and
`ML` are reconstructed against the output window: skip counts renumbered,
probability bytes re-sliced, and `MN` updated (or added when the input had
none). A group whose listed positions all fall outside the window is kept with
no positions (`C+m;`), since a group with `.` or no status declares its unlisted
bases canonical and removing it would turn them into no-calls. Groups absent
from the input stay absent. The rest of the record is copied unchanged.

A block that cannot be placed on the sequence is removed from the output record
rather than repaired: an `MN` that disagrees with the sequence length, an `ML`
that is not a `B:C` array of the length `MM` declares, or an `MM` that does not
parse to its end, or modification coordinates beyond the available counting-base
occurrences. Such reads are counted and reported in the run summary.

### Decode-equivalence tests

The test suite re-decodes whittle's output with `rust-htslib`'s
`basemods_iter()`, an independent `MM`/`ML` implementation, and compares the
result with the original calls restricted to the surviving window. One test
runs on a synthetic fixture; another sweeps a uBAM named by `WHITTLE_UBAM`:

```bash
WHITTLE_UBAM=/path/to/reads.ubam cargo test --test bam_mods_oracle -- --ignored
```

## Tag handling

| Tag(s) | On a trimmed read |
|---|---|
| `MM` / `ML` / `MN` | Reconstructed for the output window; a group that loses every position is kept empty; a malformed block is removed and counted |
| Per-base kinetics (`ip`/`pw`/`fi`/`fp`), PacBio per-base match and mismatch counts (`sm`/`sx` as `B:C`), and any read-length `B` array | Sliced with the sequence; dorado's scalar `sm:f` is copied |
| Reverse-strand kinetics (`ri`/`rp`) | Sliced from the opposite end, since the PacBio BAM specification stores them last base first |
| `sa` (PacBio run-length subread coverage, `B:I`) | Decoded to per-base coverage, sliced, and re-encoded as `<length>,<coverage>` runs; runs that do not sum to the read length leave the tag unchanged and count as malformed |
| Fixed-size PacBio arrays (`sn`/`ac`/`bc`) | Copied verbatim |
| ONT signal (`mv`/`ts`/`ns`/`sp`) | Removed, or rewritten with `--update-moves`; `sp` is kept on a crop, since the raw signal it places is unchanged |
| `pi` (parent read id) | Set to the parent's name on every ONT split segment, with or without `--update-moves`; kept on a crop |
| Poly-A (`pa`/`pt`) | Kept and shifted with `--update-moves` when the tail survives, otherwise removed; `pa` positions are absolute POD5 sample indexes, the frame `ts` uses |
| `bi` (barcode positions) | Read with any adapter source to place the barcode trim where a barcode sequence is found, then removed, since the positions index the untrimmed read; a tag that is not a seven-element `B:f` array, or whose positions describe an empty, inverted, or out-of-range window, leaves the read untrimmed and is counted |
| `BC`/`bv` (barcode call and kit version) | Per-read labels, copied unchanged |
| `ds`/`ls` (PacBio undo blobs for `skera undo` and `lima-undo`) | Removed from every output record of a trimmed read, since they describe the untrimmed read; counted once per read and reported |
| `qs:f` (dorado mean qscore) | Recomputed from the trimmed quality |
| `qs:i`/`qe:i` (PacBio query coordinates) | Rewritten as `qs + start` and `qs + end` of the window, since the PacBio BAM specification keeps them relative to the original read; one without the other leaves both unchanged |
| Read name | A crop updates an existing PacBio query interval in the name; other names are kept. A split names ONT segments `{name}_segment_N`; PacBio segments (an integer `qs`, or a `{movie}/{zmw}/ccs[/fwd|/rev]` or `{movie}/{zmw}/{qStart}_{qEnd}` name) take the specification's `{stem}/{qStart}_{qEnd}` from the rewritten coordinates, replacing any existing interval |
| `rn` (read number) | Kept on a crop; `-1` on an ONT split (dorado's convention); PacBio's `rn` is a pass count and is copied |
| `st`/`du` (start time, duration) | Kept on a crop; with `--update-moves`, a crop that shortens `ns` scales `du` with it, so `ns` over `du` stays the sample rate. On a split, recomputed with `--update-moves`, otherwise removed; pbmarkdup's `du:Z` is not a duration and is copied |
| `me`/`er` (MinKNOW event count, end reason) | On an ONT split, `me` is 0 on every segment and `er` is `unknown` except on the segment retaining the parent signal end when moves are rewritten, or the last sequence segment otherwise; only when the source carries them |
| `RG`, `ch`, `mx`, `sd`/`sv`, and other scalar tags | Copied verbatim |
| `@RG` `tm` (dorado trim mode, BAM header) | The adapter, primer and barcode classes of the resolved adapter set are merged in, in dorado's `adapter,primer,barcode` order; a read group without it gains it; a value outside that grammar is left unchanged; quality trimming and crops alone leave it unchanged |

## Filtering by aux tag

`--tag-filter <EXPR>` selects reads by their aux tags (dorado's `er`, `dx`,
`qs`, `BC`, PacBio's `rq`, `np`, and any other scalar or string tag) before
discovery and trimming, on BAM and tagged FASTQ input. See
[cli.md](cli.md#filtering-by-aux-tag) for the syntax.

## Rejection reason

Records written to `--rejected-output` carry `wr:Z:<reason>`, in the aux data of
a BAM record or as a header field of a FASTQ record; see
[cli.md](cli.md#rejected-outputput) for the reasons.

## Tag removal

`--remove-tag <TAGS>` removes named aux tags from every output record; the
groups `kinetics` (`ip`, `pw`, `fi`, `fp`, `ri`, `rp`, `sa`, `sm`, `sx`),
`mods` (`MM`, `ML`, `MN`) and `signal` (`mv`, `ts`, `ns`, `sp`, `pi`) name
several at once. Removal runs after the rewrites in
the table above, so the remaining tags stay in register. It applies to BAM
output and to the tags carried into a BAM-to-FASTQ header, and requires BAM
input. See [cli.md](cli.md#tag-removal).

## Barcode positions

With an adapter source, the barcode spans recorded in `bi` are removed where a
barcode sequence is found at them ([cli.md](cli.md#barcode-positions)),
using an interval in the original read. Every position-indexed tag is
rewritten against each final interval. The tag is a `B:f` array of seven floats, `[barcode_score,
front_start_index, front_len, front_score, rear_end_index, rear_len,
rear_score]` (dorado `read_pipeline/base/messages.cpp`). `front_start_index +
front_len` is the last base of the front barcode and `rear_end_index -
rear_len` is the first base of the rear one, so the kept window is
`[front_start_index + front_len + 1, rear_end_index - rear_len)`, the interval
dorado's trimmer keeps (`demux/Trimmer.cpp`). A barcode that dorado did not
find is stored as a negative position, and each end is guarded on its own
value, so a read barcoded at one end is trimmed at that end only.

Each adapter-derived segment is intersected with the original barcode
interval before `--trim-front` and `--trim-tail` crop its ends. The barcode
coordinates are never interpreted relative to a split segment. BAM and tagged
FASTQ input support this operation.

Adapter and quality splitting share the original coordinate frame. Final
segments are numbered in original-read order before filtering, with a single
`_segment_N` suffix for ONT records. PacBio records use the final query
coordinates. Modification calls, kinetics, and signal tags are reconstructed
from the original record for each surviving interval.

## Platform rules

A record is treated as PacBio when it carries an integer `qs` (dorado's `qs` is
a float) or its name follows a PacBio convention: `{movie}/{zmw}/ccs`, with an
optional `/fwd` or `/rev` and an optional `/{qStart}_{qEnd}`, or the subread
form `{movie}/{zmw}/{qStart}_{qEnd}`. Every other record is treated as ONT. The
platform selects the split naming and the `rn`, `pi`, `me`, and `er` rules
above; the remaining rules are keyed on the tag type.

## ONT signal tags under `--update-moves`

Signal direction comes from `basecall_model=dna...` or `basecall_model=rna...`
in the record's `@RG` description. A record without `RG` can use a direction
shared by every read group in the header. Rewriting a trimmed read's move table
requires a resolvable direction; otherwise the run fails with the read name.
RNA sequence intervals are reversed into signal order before locating move
boundaries, since RNA basecalls and move tables have opposite directions.

With `--update-moves`, a crop slices `mv`, advances `ts` by the removed head
signal, and sets `ns` to the end of the kept signal. A DNA head-only crop and
an RNA tail-only crop preserve the original signal end. A split emits subreads
in dorado's convention: `pi` parent id, `sp` offset from the parent's POD5 signal start, `ns` subread span, `ts` 0, so each
renamed segment stays locatable in POD5. A split also recomputes `du` and `st`
as dorado does: the sample rate is the parent's `ns` over its `du`, the
subread's `du` is its sample count at that rate, and its `st` is the parent's
advanced by the subread's start sample, written at millisecond precision in the
parent's offset form (`Z`, a numeric offset, or none). A missing or unusable
`ns`, `du`, or `st` leaves the tag unchanged.

BAM-to-FASTQ always removes the signal tags on a trim. A known per-base tag
whose length does not match the sequence is left untouched and reported in a
one-line advisory; a malformed modification block is removed from the record
and reported the same way.
