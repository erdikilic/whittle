# Adapter trimming

Adapter trimming is off by default. It is enabled by `-a`/`--adapter-fasta <FILE>`
(user-supplied sequences, one per record, each at least 11 bp) and/or
`--adapter-preset <KITS>` (the built-in catalog, scoped to the named kits).
Either one alone is enough, and they combine.

Adapter sequences may use the full IUPAC alphabet, so a degenerate primer is
written the way it is designed: `ACGTACGTYCGTACGRACGT` matches reads carrying
either base at each wobble position. `U` is folded to `T`, since a DNA read
stores `T`. A record holding anything outside the nucleotide alphabet is
malformed rather than degenerate, and is skipped with a warning; so is one
shorter than the 11-bp minimum. A pattern averaging two or more bases per
position is still searched, with a warning that a pattern that degenerate
matches almost anywhere.

An ambiguity code in a read is treated the other way round: as a mismatch, not
a free match. An uncalled base is evidence of nothing, so it costs error budget.
A stray `N` inside a real adapter still matches within `--adapter-error-rate`,
while a run of them never looks like an adapter and is not excised.

Every sequence has a role. An **adapter** is trimmed at the read ends and
excised in the interior; a **primer** or **barcode** is trimmed at the ends
only, since a primer or barcode inside a read is part of the molecule as often
as it is a chimera signal. Catalog entries carry their role; a FASTA entry is an
adapter unless the word `primer` or `barcode` appears in its header description
(`>27F primer`). A preset made only of amplicon kits (`mab114`) promotes its
primers to adapters, because an amplicon library has no molecule with a primer
inside it.

Every sequence is searched on both strands, so orientation does not matter and
a rear adapter is the reverse complement of its front adapter for every
purpose. Each read gets these treatments:

- **Terminal trimming.** A hit within `--adapter-end-size` bases of an end
  (default 150) trims that end, along with everything outboard of it.
- **Partial adapters.** An adapter cut short by the read end (a truncated rear
  adapter, a front adapter missing its first bases) is trimmed when at least
  10 of its bases align flush with the end. The first 10 bases must match
  exactly; the error rate applies to the rest.
- **Chimera splitting.** An adapter in the interior is treated as a junction.
  The read splits there, the adapter is excised, and both sides are kept.
  Each side is then searched again at its new end, so a primer, barcode or
  truncated adapter left next to the junction is trimmed as it would be at a
  physical read end. Two excisions separated by fewer bases than
  `--min-length` (or by at most 11 bases) merge into one. `--adapter-ends-only`
  turns splitting off and searches only the two end-zones.

Interior hits use half the `--adapter-error-rate` budget (default 0.2) that
terminal hits do, so a marginal end match still trims but only a tight interior
match splits a read. Adapter trims flow through the same tag-rewrite machinery as
every other trim, so `MM`/`ML`/`MN` and the per-base tags stay correct (see
[tags.md](tags.md)).

## Presence detection

A preset holds sequences a given library does not carry: the primers of a
kit's cDNA variant, or most of a union such as `ont`. Presence detection runs
the trimming passes over the first `--adapter-sample` reads (default 2000,
minimum 100) and keeps the entries that trimmed or split at least 0.2% of
them, then trims the rest of the input against that set. This is faster and
avoids spurious trims from absent entries.

Detection is preset-only; `--adapter-sample 0` turns it off, and a custom
`--adapter-fasta` is always searched in full. If detection finds nothing (an
ordered file with clean reads first can look adapter-free), whittle warns and
falls back to the full set rather than skipping trimming for the rest of the run.

## Ab-initio inference

`--adapter-infer` (the same as `--adapter-infer trim`) discovers recurrent
read-end sequences de novo from a sampled read prefix, using Porechop_ABI-style
k-mer assembly, then trims and splits with what it found. By default it uses a
conservative anchor of at most 32 bp facing the physical end. Anything longer on
the insert-facing side of the assembled consensus is reported as uncertain rather
than assumed to be technical.

This matters for amplicons: without a known primer or reference, a primer and a
conserved marker-gene prefix can be statistically indistinguishable. A
candidate without a sharp boundary that lies inward of another candidate in
the reads carrying both is dropped as that candidate's insert-facing
remainder, so a conserved gene start is not trimmed off reads that carry no
adapter.

`--adapter-infer report` prints the recommended anchor, its support, the assembled
length, the uncertain-base count, and any catalog/FASTA cross-name, all as FASTA,
then exits without touching record output. Add `-v` to log the full assembled
consensus.

`--adapter-infer-policy aggressive` trims the full consensus; it is appropriate
only after overtrimming of conserved biological sequence has been ruled out,
and it is unsuitable for amplicons, where the consensus grows into the
conserved gene start. The default policy is `conservative`.

## Built-in catalog and presets

`--adapter-preset <KITS>` loads the catalog entries of the named kits, as a
comma-separated list of tokens:

| Token | Kit | Contents |
|---|---|---|
| `lsk114` | Ligation sequencing kit V14 | ligation Y-adapter (kit-14 and legacy), cDNA and 10X primers |
| `rad114`, `ulk114` | Rapid and ultra-long kits V14 | rapid adapter, ligation adapter |
| `rbk114` | Rapid barcoding kit V14 | rapid adapter, rapid flank, 96 barcodes |
| `nbd114` | Native barcoding kit V14 | ligation adapter, native flank, 96 barcodes |
| `pcb114`, `pcs114` | cDNA-PCR sequencing and barcoding kits V14 | ligation adapter, PCS primers, PCR flanks, 24 barcodes |
| `rpb114` | Rapid PCR barcoding kit V14 | rapid adapter, PCR primers, rapid ligation flank, 24 barcodes |
| `mab114` | Microbial amplicon barcoding kit V14 | ligation and rapid adapters, degenerate 16S 27F/1492R and ITS1F/ITS4 primers, MAB flanks, 24 barcodes |
| `rna004` | Direct RNA sequencing kit | RNA adapter |
| `pacbio` | PacBio SMRTbell | SMRTbell adapter, C2 sequencing primer |
| `ont` | every ONT kit above | |
| `all` | every kit | |

The catalog is assembled for whittle from vendor-published sources: dorado's
`adapter_primer_kits.cpp` and `barcode_kits.cpp`, Porechop's `adapters.py`,
qcat's kit definitions, and the NCBI UniVec records of the PacBio sequences.
Every entry names its kits, so a kit preset searches only what its library can
contain; a union searches everything, at some cost in speed and precision.

Reverse-complement search covers both orientations. Flanks under 11 bp are left
out of the catalog, since a pattern that short matches almost anywhere and is
never searched on its own.
