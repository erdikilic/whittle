# Trim-aware tags

Long-read records carry tags that are indexed by base position. Crop the sequence
without touching them and they now point at bases that are gone. whittle keeps
every one of them in register, on both BAM→BAM and BAM→FASTQ.

## Base-modification tags

Modification calls (5mC, 6mA, and so on) live per read in two tags: `MM`, which
says which bases are modified as skip-counts over the sequence, and `ML`, which
holds their probabilities. Trim `SEQ`/`QUAL` without touching these and the result
decodes to nonsense, because the skip-counts still index the original sequence
rather than the trimmed one.

whittle rebuilds them as part of trimming. For every output uBAM read, whether
cropped, quality-trimmed, or split, `MM` and `ML` are reconstructed against the
output window: skip-counts renumbered, probability bytes re-sliced, and `MN`
updated to match (or added when the input had none). A group whose listed
positions all fall outside the window is kept with no positions (`C+m;`): a
group with `.` or no status declares its unlisted bases canonical, and removing
it would turn them into no-calls. Groups absent from the input stay absent.
Everything else in the record rides through unchanged.

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
real uBAM when you point it at one:

```bash
WHITTLE_UBAM=/path/to/real.ubam cargo test --test bam_mods_oracle -- --ignored
```

## Every position-indexed tag

| Tag(s) | On a trimmed read |
|---|---|
| `MM` / `ML` / `MN` | Reconstructed for the output window; a group that loses every position is kept empty; a malformed block is removed and counted |
| Per-base kinetics (`ip`/`pw`/`fi`/`fp`), and any read-length `B` array | Sliced in lockstep with the sequence |
| Reverse-strand kinetics (`ri`/`rp`) | Sliced from the other end, since the PacBio BAM spec stores them last base first |
| Fixed-size PacBio arrays (`sn`/`ac`/`bc`) | Copied verbatim, never treated as per-base |
| ONT signal (`mv`/`ts`/`ns`/`sp`/`pi`) | Dropped, or rewritten with `--update-moves` |
| Poly-A (`pa`/`pt`) | Kept/shifted with `--update-moves` if the tail survives, else dropped; `pa` positions are absolute POD5 sample indexes, the frame `ts` uses |
| `bi` (barcode positions) | Dropped, since the positions shift under a crop |
| `qs:f` (dorado mean qscore) | Recomputed from the trimmed quality; PacBio's `qs:i`/`qe:i` query coordinates are copied verbatim |
| `rn` (read number) | Kept on a crop; `-1` on a split, with or without `--update-moves` |
| `st`/`du` (start time / duration) | Kept on a crop, dropped on a split |
| `RG`, `ch`, `mx`, `sm`/`sd`/`sv`, … | Copied verbatim |

## ONT signal tags under `--update-moves`

With `--update-moves`, a crop slices `mv` and advances `ts` (a head-only crop
leaves `ns` unchanged), while a split emits dorado-style subreads (`pi` parent
id, `sp` parent-signal offset, `ns` subread span, `ts` 0) so the renamed segment
stays locatable in POD5 for tools like Remora.

BAM→FASTQ always drops the signal tags on a trim, since a move table in a
FASTQ header is impractical. If a known per-base tag's length doesn't match the
sequence (malformed input), whittle leaves it untouched and prints a one-line
advisory. A malformed modification block is removed from the record and reported
the same way.
