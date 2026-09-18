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

Interior hits use a stricter budget than terminal hits. For each adapter and
read-length class (powers of two from 4 kb upward), the interior budget is the
largest edit count whose expected number of chance matches in a read of that
length, over both strands under independent uniform DNA, stays within `1e-4`;
it never exceeds the terminal budget, and exact matches are always accepted.
Long, informative patterns keep most of their tolerance in long reads; short or
IUPAC-rich patterns and very long reads receive fewer interior edits. This null
model is a specificity guard, not a calibrated false-split probability for
repetitive or composition-biased biological reads.
Adapter trims pass through the same tag-rewrite path as
every other trim, so `MM`/`ML`/`MN` and the per-base tags stay in register
([tags.md](tags.md)).

Unknown read bases consume edit budget and cannot form exact adapter matches.

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

Adapter samples are bounded by 256 MiB of retained payload and 64 Mi bases as
well as `--adapter-sample-reads`. The last sampled record is retained whole;
the actual sample size is reported, and all sampled reads are processed.

## Adapter discovery

`--discover-adapters` (equivalent to `--discover-adapters trim`) discovers recurrent
sequences from the first `--adapter-sample-reads` reads (default 40000), then trims
and splits with the supported sequences. Sampling requires at least 100 reads;
a count of 0 is rejected. Reporting and trimming use the same automatic boundary
rule. There is no conservative/aggressive policy switch.

Discovery counts exact 16-mers in the first and last 100 bases of sampled reads,
once per read-end window. Short tandem-repeat seeds are excluded.
A bounded graph assembly retains multiple paths, locates abrupt support changes
at insert boundaries, and corrects weak paths with aligned read evidence.
Sassy batches the approximate searches used to align supporting reads and
validate complete candidates. At most 4000 windows per end, distributed across
the sample and extending to 200 bases, participate in alignment validation.
Insert-facing termination is checked against both the original graph and
continuing bases in these longer windows, so the assembly window itself does
not supply evidence for an adapter boundary.

A candidate must occur in at least 1% of the usable validation windows at one
end and in at least 20 windows. At least 80% of its supporting alignments must
lie within 35 bases of the physical read end. Candidates dominated by short
approximate repeats are rejected; candidate prevalence must exceed that in
adjacent interior windows by more than fourfold. These checks allow minority
families without treating any recurrent sequence as an adapter. Rare families,
large barcode panels, distant adapters and insufficient samples can still be
missed. The bounded graph considers up to 12 paths per end; a weak fragment
does not terminate the search for other families.

Related candidates and reverse complements are merged after prioritizing
independently supported primer boundaries, with exact sequence support choosing
among their reconstructions. Other candidates are ranked by exact graph support
before complete-sequence support, so short error-derived fragments do not
displace a stronger assembly merely because they match more reads.

For amplicons, recurrent read starts can identify an insert boundary in reads
that lack a primer. Discovery aligns the upstream sequence in primer-bearing
reads, reconstructs supported IUPAC variants, and excludes the conserved insert.
Distinct primer families sharing an insert boundary are considered separately.
Candidates overlapping a validated insert start, or supported mainly on the
insert-facing side of a reconstructed primer, are excluded. Candidates whose
insert-facing boundaries remain unresolved are skipped.

Without reads exposing an insert boundary, an unknown primer and a conserved
gene prefix can be indistinguishable. Short primers, highly degenerate mixtures,
high sequencing error, adapter lengths approaching the sampling-window length,
or insufficient representation can prevent complete recovery. Use a known
primer FASTA or kit preset when available. Inferred sequences use the adapter
role: they trim ends and can excise interior junctions. Use
`--adapter-ends-only` to suppress adapter splitting, or an explicit FASTA with
`primer` in its headers to restrict primer matches to the ends.

`--discover-adapters report` prints the same sequences used for trimming,
together with their support, length, and catalog or supplied-FASTA annotations,
then exits without writing records or a JSON summary. Catalog entries provide
names only; their sequences do not choose or extend an inferred consensus.
`-v` logs inferred sequences and unresolved-boundary candidates.

The former `--adapter-discovery-policy` option and JSON `infer_policy` parameter
have been removed. Existing commands should omit that option.

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
