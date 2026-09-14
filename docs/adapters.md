# Adapter trimming

Adapter trimming is off by default. It is enabled by an adapter source:
`-a`/`--adapter-fasta <FILE>` (user-supplied sequences), `--adapter-preset
<KITS>` (the built-in catalog, scoped to the named kits), or `--discover-adapters`
(ab-initio discovery). A FASTA and a preset combine into one search set.

Adapter search operates on the original read before barcode restriction and
fixed cropping. Each retained adapter segment then receives barcode
restriction, fixed cropping, and quality processing. Quality splitting can
produce several final segments from one adapter segment. Final segments are
numbered in original-read order, then filtered. Reads without an
adapter match continue through the same downstream stages.

## Sequences

A FASTA record holds one sequence of at least 11 bp. The full IUPAC alphabet is
accepted, so a degenerate primer is written as designed: `AGAGTTTGATYMTGGCTCAG`
matches either base at each wobble position. `U` is folded to `T`. A record with
a character outside the nucleotide alphabet, or shorter than 11 bp, is skipped
with a warning. A pattern averaging two or more bases per position is searched
with a warning, since it matches almost anywhere.

Ambiguity codes in a read are mismatches, not free matches: an uncalled base
costs error budget. A single `N` inside an adapter still matches within
`--adapter-error-rate`; a run of them does not.

Every sequence has a role. An **adapter** is trimmed at the read ends and
excised in the interior. A **primer** or **barcode** is trimmed at the ends
only, since a primer or barcode inside a read is as often part of the molecule
as a chimera signal. Catalog entries carry their role. A FASTA entry is an
adapter unless the word `primer` or `barcode` appears in its header description
(`>27F primer`). A preset made only of amplicon kits (`mab114`) promotes its
primers to adapters, since an amplicon library has no molecule with a primer
inside it.

## Search

Every sequence is searched on both strands, so orientation is irrelevant and a
rear adapter is the reverse complement of its front adapter. Each read receives
these treatments:

- **Terminal trimming.** A hit within `--adapter-end-search` bases of an end
  (default 150) trims that end together with everything outboard of it.
- **Partial adapters.** An adapter cut short by the read end (a truncated rear
  adapter, a front adapter missing its first bases) is trimmed when at least
  10 of its bases align flush with the end. Those 10 bases must match exactly;
  the error rate applies to the remainder.
- **Chimera splitting.** An interior adapter marks a junction. The read is split
  there, the adapter is excised, and both sides are kept. Each side is searched
  again at its new end, so a primer, barcode, or truncated adapter adjacent to
  the junction is trimmed as at a physical read end. Two excisions separated by
  fewer bases than `--min-length`, or by at most 11 bases, merge into one.
  `--adapter-ends-only` disables splitting and searches only the two end zones.

Interior hits are held to half the `--adapter-error-rate` budget (default 0.2)
of terminal hits, so a marginal end match trims while only a close interior
match splits a read. Adapter trims pass through the same tag-rewrite path as
every other trim, so `MM`/`ML`/`MN` and the per-base tags stay in register
([tags.md](tags.md)).

## Presence detection

A preset holds sequences a given library may not carry, such as the primers of
a kit's cDNA variant or most of the `ont` union. Presence detection runs the
trimming passes over the first `--adapter-sample-reads` reads (default 2000, minimum
100), keeps the entries that trimmed or split at least 0.2% of them (at least 3
reads), and searches the remaining input against that set only. This is faster
and avoids spurious trims from absent entries.

Detection applies to presets only. Supplying a FASTA disables presence
detection for the combined FASTA and preset set. `--adapter-sample-reads 0`
also disables detection and searches the full set. If detection finds nothing,
whittle warns and falls back to the full set.

## Adapter discovery

`--discover-adapters` (equivalent to `--discover-adapters trim`) discovers recurrent
read-end sequences de novo from the first `--adapter-sample-reads` reads (default
40000) by k-mer assembly in the manner of Porechop_ABI, then trims and splits
with the result. Sampling requires at least 100 reads; a sample count of 0 is
rejected. The default `conservative` policy anchors at most 32 bp facing
the physical read end; the insert-facing remainder of the assembled consensus is
reported as uncertain rather than assumed technical.

Without a known primer or reference, a primer and a conserved marker-gene
prefix can be statistically indistinguishable. A candidate without a sharp
boundary that lies inward of another candidate, in the reads carrying both, is
dropped as that candidate's insert-facing remainder, so a conserved gene start
is not trimmed from reads that carry no adapter.

`--discover-adapters report` prints the recommended anchor, its support, the
assembled length, the uncertain-base count, and any catalog or FASTA cross-name,
as FASTA to stdout, then exits without writing records or a JSON summary.
`-v` logs the full assembled consensus.

`--adapter-discovery-policy aggressive` trims the full consensus. It is appropriate
only when overtrimming of conserved biological sequence has been ruled out, and
unsuitable for amplicons, where the consensus extends into the conserved gene
start.

## Catalog and presets

`--adapter-preset <KITS>` loads the catalog entries of the named kits, given as
a comma-separated list of tokens:

| Token | Kit | Contents |
|---|---|---|
| `lsk114` | Ligation sequencing kit V14 | Ligation Y-adapter (kit 14 and legacy), cDNA and 10X primers |
| `rad114`, `ulk114` | Rapid and ultra-long kits V14 | Rapid adapter, ligation adapter |
| `rbk114` | Rapid barcoding kit V14 | Rapid adapter, rapid flank, 96 barcodes |
| `nbd114` | Native barcoding kit V14 | Ligation adapter, native flank, 96 barcodes |
| `pcb114`, `pcs114` | cDNA-PCR sequencing and barcoding kits V14 | Ligation adapter, PCS primers, PCR flanks, 24 barcodes |
| `rpb114` | Rapid PCR barcoding kit V14 | Rapid adapter, PCR primers, rapid ligation flank, 24 barcodes |
| `mab114` | Microbial amplicon barcoding kit V14 | Ligation and rapid adapters, degenerate 16S 27F/1492R and ITS1F/ITS4 primers, MAB flanks, 24 barcodes |
| `rna004` | Direct RNA sequencing kit | RNA adapter |
| `pacbio` | PacBio SMRTbell | SMRTbell adapter, C2 sequencing primer |
| `ont` | Every ONT kit above | |
| `all` | Every kit | |

The catalog is assembled from vendor-published sources: dorado's
`adapter_primer_kits.cpp` and `barcode_kits.cpp`, Porechop's `adapters.py`,
qcat's kit definitions, and the NCBI UniVec records of the PacBio sequences.
Every entry names its kits, so a kit preset searches only what its library can
contain; a union searches everything, at some cost in speed and precision.
Flanks shorter than 11 bp are excluded, since a pattern that short matches
almost anywhere.
