# Adapter trimming

Off by default. Turn it on with `-a`/`--adapter-fasta <FILE>` (your own sequences,
one per record, each at least 11 bp) and/or `--adapter-preset ont` (the built-in
catalog). Either one alone is enough, and they combine.

Adapter sequences may use the full IUPAC alphabet, so a degenerate primer is
written the way it is designed: `ACGTACGTYCGTACGRACGT` matches reads carrying
either base at each wobble position. `U` is folded to `T`, since a DNA read
stores `T`. A record holding anything outside the nucleotide alphabet is
malformed rather than degenerate, and is skipped with a warning; so is one
shorter than the 11-bp minimum. A pattern averaging two or more bases per
position is still searched, since you asked for it, but it warns that something
that degenerate matches almost anywhere.

An ambiguity code in a *read* is treated the other way round: as a mismatch, not
a free match. An uncalled base is evidence of nothing, so it costs error budget.
A stray `N` inside a real adapter still matches within `--adapter-error-rate`,
while a run of them never looks like an adapter and is not excised.

Every adapter is searched on both strands, so orientation doesn't matter, and
each read gets two treatments:

- **Terminal trimming.** An adapter within `--adapter-end-size` bases of an end
  (default 150) is trimmed off.
- **Chimera splitting.** An adapter in the interior is treated as a junction. The
  read splits there, the adapter is excised, and both sides are kept.
  `--adapter-ends-only` turns this off and searches only the two end-zones.

Interior hits use half the `--adapter-error-rate` budget (default 0.2) that
terminal hits do, so a marginal end match still trims but only a tight interior
match splits a read. Adapter trims flow through the same tag-rewrite machinery as
every other trim, so `MM`/`ML`/`MN` and the per-base tags stay correct (see
[tags.md](tags.md)).

## Presence detection

A preset catalog holds far more adapters than any single run uses (the ONT one has
over a hundred). `--adapter-sample <N>` (N >= 100) checks which adapters actually
turn up in the first N reads, then trims the rest against only that set. It's
faster, and it avoids spurious trims from catalog entries that aren't present.

Detection is off by default (`--adapter-sample 0`) and preset-only; a custom
`--adapter-fasta` is always searched in full. If detection finds nothing (an
ordered file with clean reads first can look adapter-free), whittle warns and
falls back to the full set rather than skipping trimming for the rest of the run.

## Ab-initio inference

`--adapter-infer` (the same as `--adapter-infer trim`) discovers recurrent
read-end sequences de novo from a sampled read prefix, using Porechop_ABI-style
k-mer assembly, then trims with what it found. By default it trims ends only, with
a conservative anchor of at most 32 bp facing the physical end. Anything longer on
the insert-facing side of the assembled consensus is reported as uncertain rather
than assumed to be technical.

This matters for amplicons: without a known primer or reference, a primer and a
conserved marker-gene prefix can be statistically indistinguishable.

`--adapter-infer report` prints the recommended anchor, its support, the assembled
length, the uncertain-base count, and any catalog/FASTA cross-name, all as FASTA,
then exits without touching record output. Add `-v` to log the full review-only
consensus.

`--adapter-infer-policy aggressive` restores full-consensus trimming and allows
interior splitting unless `--adapter-ends-only` is also set; reach for it only
once you've ruled out overtrimming conserved biological sequence. The default
policy is `conservative`.

## Built-in ONT catalog

`--adapter-preset ont` loads a catalog assembled for whittle from ONT-published
sources: dorado's `adapter_primer_kits.cpp`, Porechop's `adapters.py`, and qcat's
kit definitions. It has 120 sequences: ligation adapters (kit-14 and legacy), the
rapid and direct-RNA adapters, PCR/cDNA and 10X primers, barcode flanks, and all
96 barcodes.

Reverse-complement search covers both orientations. Flanks under 11 bp are left
out of the catalog, since a pattern that short matches almost anywhere and is
never searched on its own.
