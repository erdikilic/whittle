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
parse to its end. Such reads are counted and reported in the run summary.

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
| ONT signal (`mv`/`ts`/`ns`/`sp`) | Removed, or rewritten with `--update-signal-tags` |
| `pi` (parent read id) | Set to the parent's name on every ONT split segment, with or without `--update-signal-tags`; removed on a crop without it |
| Poly-A (`pa`/`pt`) | Kept and shifted with `--update-signal-tags` when the tail survives, otherwise removed; `pa` positions are absolute POD5 sample indexes, the frame `ts` uses |
| `bi` (barcode positions) | Read by `--trim-barcodes` to place the trim, then removed, since the positions index the untrimmed read; a tag that is not a seven-element `B:f` array, or whose positions describe an empty, inverted, or out-of-range window, leaves the read untrimmed and is counted |
| `BC`/`bv` (barcode call and kit version) | Per-read labels, copied unchanged |
| `ds`/`ls` (PacBio undo blobs for `skera undo` and `lima-undo`) | Removed from every output record of a trimmed read, since they describe the untrimmed read; counted once per read and reported |
| `qs:f` (dorado mean qscore) | Recomputed from the trimmed quality |
| `qs:i`/`qe:i` (PacBio query coordinates) | Rewritten as `qs + start` and `qs + end` of the window, since the PacBio BAM specification keeps them relative to the original read; one without the other leaves both unchanged |
| Read name | Kept on a crop. A split names ONT segments `{name}_segment_N`; PacBio segments (an integer `qs`, or a `{movie}/{zmw}/ccs[/fwd|/rev]` or `{movie}/{zmw}/{qStart}_{qEnd}` name) take the specification's `{stem}/{qStart}_{qEnd}` from the rewritten coordinates, replacing any existing interval |
| `rn` (read number) | Kept on a crop; `-1` on an ONT split (dorado's convention); PacBio's `rn` is a pass count and is copied |
| `st`/`du` (start time, duration) | Kept on a crop. On a split, recomputed with `--update-signal-tags`, otherwise removed; pbmarkdup's `du:Z` is not a duration and is copied |
| `me`/`er` (MinKNOW event count, end reason) | On an ONT split, `me` is 0 on every segment and `er` is `unknown` on every segment but the last, which ends where the read did; only when the source carries them |
| `RG`, `ch`, `mx`, `sd`/`sv`, and other scalar tags | Copied verbatim |

## Tag removal

`--remove-tag <TAG>` removes a named aux tag from every output record, and
`--remove-kinetics` removes the nine per-base arrays (`ip`, `pw`, `fi`, `fp`,
`ri`, `rp`, `sa`, `sm`, `sx`) in one flag. Removal runs after the rewrites in
the table above, so the remaining tags stay in register. It applies to BAM
output and to the tags carried into a BAM-to-FASTQ header, and requires BAM
input. See [cli.md](cli.md#tag-removal).

## Barcode positions under `--trim-barcodes`

`--trim-barcodes` removes the barcode spans recorded in `bi` rather than
searching for the sequences, so the cut equals the one `dorado demux` makes,
and every position-indexed tag is rewritten by the same machinery as any other
trim. The tag is a `B:f` array of seven floats, `[barcode_score,
front_start_index, front_len, front_score, rear_end_index, rear_len,
rear_score]` (dorado `read_pipeline/base/messages.cpp`). `front_start_index +
front_len` is the last base of the front barcode and `rear_end_index -
rear_len` is the first base of the rear one, so the kept window is
`[front_start_index + front_len + 1, rear_end_index - rear_len)`, the interval
dorado's trimmer keeps (`demux/Trimmer.cpp`). A barcode that dorado did not
find is stored as a negative position, and each end is guarded on its own
value, so a read barcoded at one end is trimmed at that end only.

Barcode removal is the outermost stage: it runs before `--trim-front` and
`--trim-tail`, which therefore count from the first base after the front
barcode. It requires BAM input, since no other format carries the tag.

## Platform rules

A record is treated as PacBio when it carries an integer `qs` (dorado's `qs` is
a float) or its name follows a PacBio convention: `{movie}/{zmw}/ccs`, with an
optional `/fwd` or `/rev` and an optional `/{qStart}_{qEnd}`, or the subread
form `{movie}/{zmw}/{qStart}_{qEnd}`. Every other record is treated as ONT. The
platform selects the split naming and the `rn`, `pi`, `me`, and `er` rules
above; the remaining rules are keyed on the tag type.

## ONT signal tags under `--update-signal-tags`

With `--update-signal-tags`, a crop slices `mv`, advances `ts` by the removed head
signal, and sets `ns` to the end of the kept signal (unchanged under a head-only
crop). A split emits subreads in dorado's convention: `pi` parent id, `sp`
offset from the parent's POD5 signal start, `ns` subread span, `ts` 0, so each
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
