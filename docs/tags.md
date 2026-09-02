# Trim-aware tags

Long-read records carry tags that are indexed by base position. Cropping the
sequence without touching them leaves them pointing at bases that are gone.
whittle keeps every one of them in register, on both BAM-to-BAM and BAM-to-FASTQ.

## Base-modification tags

Modification calls (5mC, 6mA, and so on) live per read in two tags: `MM`, which
says which bases are modified as skip-counts over the sequence, and `ML`, which
holds their probabilities. Trimming `SEQ`/`QUAL` without touching these leaves a
block that decodes to nonsense, because the skip-counts index the original
sequence rather than the trimmed one.

whittle rebuilds them as part of trimming. For every output uBAM read, whether
cropped, quality-trimmed, or split, `MM` and `ML` are reconstructed against the
output window: skip-counts renumbered, probability bytes re-sliced, and `MN`
updated to match (or added when the input had none). A group whose listed
positions all fall outside the window is kept with no positions (`C+m;`): a
group with `.` or no status declares its unlisted bases canonical, and removing
it would turn them into no-calls. Groups absent from the input stay absent.
Everything else in the record is copied unchanged.

A block that cannot be placed on the sequence is removed from the output record
rather than repaired: an `MN` that disagrees with the sequence length (the
sequence was altered after the calls were made), an `ML` that is not a `B:C`
array of the length `MM` declares, or an `MM` that does not parse to its end.
Such reads are counted and reported in the run summary.

### Decode-equivalence testing

This is covered by decode-equivalence tests. They re-decode whittle's output with
`rust-htslib`'s `basemods_iter()`, a different `MM`/`ML` implementation from the
one whittle writes with, and compare against the original calls restricted to the
surviving window. One test always runs on a synthetic fixture; another sweeps a
real uBAM when `WHITTLE_UBAM` names one:

```bash
WHITTLE_UBAM=/path/to/real.ubam cargo test --test bam_mods_oracle -- --ignored
```

## Every position-indexed tag

| Tag(s) | On a trimmed read |
|---|---|
| `MM` / `ML` / `MN` | Reconstructed for the output window; a group that loses every position is kept empty; a malformed block is removed and counted |
| Per-base kinetics (`ip`/`pw`/`fi`/`fp`), PacBio per-base match/mismatch counts (`sm`/`sx` as `B:C`), and any read-length `B` array | Sliced in lockstep with the sequence; dorado's scalar `sm:f` is not an array and is copied |
| Reverse-strand kinetics (`ri`/`rp`) | Sliced from the other end, since the PacBio BAM spec stores them last base first |
| `sa` (PacBio run-length subread coverage, `B:I`) | Decoded to per-base coverage, sliced, and re-encoded as `<length>,<coverage>` runs; runs that do not sum to the read length leave the tag unchanged and count as malformed |
| Fixed-size PacBio arrays (`sn`/`ac`/`bc`) | Copied verbatim, never treated as per-base |
| ONT signal (`mv`/`ts`/`ns`/`sp`) | Dropped, or rewritten with `--update-moves` |
| `pi` (parent read id) | Set to the parent's name on every ONT split segment, with or without `--update-moves`; dropped on a crop without it |
| Poly-A (`pa`/`pt`) | Kept/shifted with `--update-moves` if the tail survives, else dropped; `pa` positions are absolute POD5 sample indexes, the frame `ts` uses |
| `bi` (barcode positions) | Read by `--trim-barcodes` to place the trim, then dropped, since the positions index the untrimmed read; a tag that is not a seven-element `B:f` array, or whose positions describe an empty, inverted or out-of-range window, leaves the read untrimmed and is counted |
| `BC`/`bv` (barcode call and kit version) | Per-read labels, copied unchanged |
| `ds`/`ls` (PacBio undo blobs for `skera undo` and `lima-undo`) | Dropped from every output record of a trimmed read, since they describe the untrimmed read; counted once per read and reported |
| `qs:f` (dorado mean qscore) | Recomputed from the trimmed quality |
| `qs:i`/`qe:i` (PacBio query coordinates) | Rewritten as `qs + start` and `qs + end` of the window, since the PacBio BAM spec keeps them with respect to the original read; one without the other leaves both unchanged |
| Read name | Kept on a crop. A split names ONT segments `{name}_segment_N`; PacBio segments (an integer `qs`, or a `{movie}/{zmw}/ccs[/fwd|/rev]` or `{movie}/{zmw}/{qStart}_{qEnd}` name) take the spec's `{stem}/{qStart}_{qEnd}` from the rewritten coordinates, with any existing interval replaced |
| `rn` (read number) | Kept on a crop; `-1` on an ONT split (dorado's convention); PacBio's `rn` is a pass count and passes through |
| `st`/`du` (start time / duration) | Kept on a crop. On a split, recomputed with `--update-moves` (below), else dropped; pbmarkdup's `du:Z` is not a duration and passes through |
| `me`/`er` (MinKNOW event count / end reason) | On an ONT split, `me` is 0 on every segment and `er` is `unknown` on every segment but the last, which ends where the read did; both only when the source carries them |
| `RG`, `ch`, `mx`, `sd`/`sv`, and other scalar tags | Copied verbatim |

## Barcode positions under `--trim-barcodes`

`--trim-barcodes` removes the barcode spans dorado recorded in `bi` rather than
searching for the sequences, so the cut matches `dorado demux` and every
position-indexed tag is rewritten by the same machinery as any other trim. The
tag is a `B:f` array of seven floats, `[barcode_score, front_start_index,
front_len, front_score, rear_end_index, rear_len, rear_score]`
(`read_pipeline/base/messages.cpp`). `front_start_index + front_len` is the last
base of the front barcode and `rear_end_index - rear_len` is the first base of
the rear one, so the kept window is
`[front_start_index + front_len + 1, rear_end_index - rear_len)`, the interval
dorado's own trimmer keeps (`demux/Trimmer.cpp`). A barcode dorado did not find
is written as a negative position, and each end is guarded on its own value, so
a read barcoded at one end only is trimmed at that end only.

Barcode removal is the outermost stage: it runs before `--head-crop` and
`--tail-crop`, which therefore count from the first base after the front
barcode. It is BAM only, since no other input carries the tag.

## Platform rules

A record is treated as PacBio when it carries an integer `qs` (dorado's `qs` is
a float) or its name follows a PacBio convention: `{movie}/{zmw}/ccs`, with an
optional `/fwd` or `/rev` and an optional `/{qStart}_{qEnd}`, or the subread form
`{movie}/{zmw}/{qStart}_{qEnd}`. Every other record is treated as ONT. The
platform decides the split naming and the `rn`, `pi`, `me` and `er` rules above;
the remaining rules are keyed on the tag's type.

## ONT signal tags under `--update-moves`

With `--update-moves`, a crop slices `mv`, advances `ts` by the removed head
signal, and sets `ns` to the end of the kept signal (unchanged under a head-only
crop), while a split emits subreads in dorado's convention (`pi` parent id, `sp`
offset from the parent's POD5 signal start, `ns` subread span, `ts` 0) so the
renamed segment stays locatable in POD5 for tools like Remora. A split also
recomputes `du` and `st` the way dorado does: the sample rate is the parent's
`ns` over its `du`, the subread's `du` is its sample count at that rate, and its
`st` is the parent's advanced by the subread's start sample, written at
millisecond precision in the parent's offset form (`Z`, a numeric offset, or
none). A missing or unusable `ns`, `du` or `st` leaves the tag unchanged.

BAM-to-FASTQ always drops the signal tags on a trim, since a move table in a
FASTQ header is impractical. A known per-base tag whose length does not match the
sequence (malformed input) is left untouched and reported in a one-line advisory.
A malformed modification block is removed from the record and reported the same
way.
