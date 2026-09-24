# Adapter trimming

Adapter trimming is off by default. It is enabled by an adapter source:
`-a`/`--adapter-fasta <FILE>` (user-supplied sequences), `--adapter-preset
<KITS>` (the built-in catalog, scoped to the named kits), or `--adapter-discover`
(de novo discovery). A FASTA and a preset combine into one search set.

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
- **Partial adapters.** An adapter, primer or barcode cut short by the read
  end (a truncated rear adapter, a front adapter or barcode missing its first
  bases) is trimmed when at least 10 of its bases align flush with the end. Those 10 bases must match exactly;
  the error rate applies to the remainder.
- **Chimera splitting.** An interior adapter marks a junction. The read is split
  there, the adapter is excised, and both sides are kept. Each side is searched
  again at its new end, so a primer, barcode, or truncated adapter adjacent to
  the junction is trimmed as at a physical read end. Two excisions separated by
  fewer bases than `--min-length`, or by at most 11 bases, merge into one.
  `--adapter-ends-only` disables splitting and searches only the two end zones.

A trim or excision boundary lies at the last well-aligned base of its hit.
Each end of the alignment is scored `+1` per match and `-2` per mismatch or
gap, and a run of end columns scoring below zero is left to the read. A catalog
entry that continues past the sequence in the read, its extra bases aligned as
edits against the insert, therefore trims where the shared sequence ends.

The native barcode kits carry a barcode construct, `NB_1st_FRONT`, the barcode
as 24 `N`, and the 8-base `NB_1st_REAR`, so a hit trims the barcode together
with the short flank on its insert side.

The edit budget of a hit is the error rate times the pattern length, rounded
down, bounded by a chance-match null model under independent uniform DNA. The
model sums alignment paths, so it overstates the chance rate and the bounds are
conservative.

A pattern's terminal budget is the largest edit count whose expected chance
matches over both end zones and both strands stay within `0.1` per read; at the
default error rate this lowers the budget of 11-, 12- and 15-base patterns by
one edit. The same bound also holds for the search set as a whole: a panel of
interchangeable sequences, such as 24 barcodes, multiplies the chance of a hit
at a read end by its size, so the distinct sequences of the set together stay
within `0.1` expected chance matches per read. The largest contributors lose an
edit first, and the members of a panel lose it together. The bound is applied
in two tiers. A hit anchored at the read end, or directly behind an accepted
hit, starting within 11 bases of it, keeps the budget bounded over that small
anchoring zone; a hit anywhere else in the end zone is held to the budget
bounded over the whole zone. A sequence and its reverse complement count once.
Every catalog entry keeps its per-pattern budget when anchored at the read end,
and the catalog adapters and 24-base ONT barcodes keep it anywhere in the end
zone.

Interior hits use a stricter budget. For each adapter and read-length class
(powers of two from 4 kb upward), the interior budget is the largest edit count
whose expected number of chance matches in a read of that length stays within
`1e-4`, for each splitting sequence alone and for the splitting sequences of
the set together. It never exceeds the terminal budget, and exact matches are
always accepted. Long, informative patterns keep most of their tolerance in
long reads; short or IUPAC-rich patterns, large panels and very long reads
receive fewer interior edits. Because exact matches are always accepted, a
splitting entry of about 12 bases or fewer can split long reads at chance
matches; give such entries the primer or barcode role in the FASTA header, or
use `--adapter-ends-only`. This null model is a specificity guard, not a
calibrated false-split probability for repetitive or composition-biased
biological reads.
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

`--adapter-discover` discovers recurrent sequences from the first
`--adapter-sample-reads` reads (default 40000), then trims and splits with the
supported sequences. Sampling requires at least 100 reads;
a count of 0 is rejected. Reporting and trimming use the same automatic boundary
rule. There is no conservative/aggressive policy switch.

Discovery proceeds in layers. Sequences from `--adapter-preset` or
`--adapter-fasta` are searched first and explain the outermost layers; a
preset given with `--adapter-discover` is therefore trimmed as usual and
discovery continues from the boundary it leaves, which recovers an unknown
primer behind a known kit. Each accepted layer moves the boundary inward and
the next layer is assembled from the unexplained sequence, so an adapter
remnant, a barcode flank and the barcodes of a rapid barcoding library are
found in turn. A layer ends where the support of adjacent k-mers changes
fourfold, such as where a barcode joins its shared flank, and the sequence
past that point belongs to the next layer. A variable layer is resolved from
the reads rather than the graph: when the best supported candidate lies at
least 12 bases past the boundary and every candidate at the boundary is
fourfold rarer, the stretch before it holds one member per read, and those
members are clustered by edit distance into families of at least 1% of the
windows and 20 reads. This recovers each barcode of a panel between its
flanks, from four to 96 barcodes, and reports them with the barcode role.
Layers stop when no candidate passes the checks below, after at most eight
per end. A discovered sequence flush with the physical read end takes the
adapter role and splits reads at interior hits; deeper sequences trim ends
only, like catalog barcodes and primers.

Discovery separates primers from the amplicon they bind by the far end of
the molecule. A conserved gene start that reads reach from the other side
appears there without the primer stack around it and is left in the read; a
primer whose reverse complement is absent from the far end, as in rapid
amplicon libraries, is trimmed. When reads run through the far primer, the
primer and the conserved start both recur at both ends and cannot be told
apart by structure: discovery then leaves the primers of one-sided libraries
and, in libraries barcoded at both ends, trims the conserved start with them.
Supply the primers with a kit preset or `--adapter-fasta` for such
libraries; discovery continues beyond them.

Discovery counts exact 16-mers in the first 100 unexplained bases at each end
of the sampled reads, once per read-end window. Short tandem-repeat seeds are
excluded.
A bounded graph assembly retains multiple paths, locates abrupt support changes
at insert boundaries, and corrects weak paths with aligned read evidence.
sassy batches the approximate searches used to align supporting reads and
validate complete candidates. At most 4000 windows per end, distributed across
the sample and extending to 200 bases, participate in alignment validation.
Insert-facing termination is checked against both the original graph and
continuing bases in these longer windows, so the assembly window itself does
not supply evidence for an adapter boundary.

The insert-facing end of each candidate is then placed by base conservation.
For each cut point near that end, the read base that follows an exact match of
the 11 candidate bases before it is tallied over the supporting windows, and
the candidate ends at the first cut point whose tally is not conserved. A
position is conserved when one base, or a pair of bases, holds at least the
midpoint between its share of the base composition of the windows and one: a
technical base is read at the platform's per-base accuracy and clears it, as
does a two-fold degenerate primer base on its pair, while an insert position
holds each base at its composition share. Judging against the composition keeps
the rule valid for AT- or GC-rich inserts up to the most extreme sequenced
genomes, around 80% AT. Assembly support alone cannot place this end for a
layer about one k-mer long: erosion at the physical read end removes the first
bases of the layer from part of the reads, which depresses the k-mer spanning
the whole layer toward the level of its insert continuations. The same tally
after the complete candidate counts as direct evidence of an insert boundary.

A candidate must occur in at least 1% of the usable validation windows at one
end and in at least 20 windows. At least 60% of its supporting alignments must
lie within 50 bases of the current boundary. A candidate whose inner segment
occurs as its reverse complement at the opposite physical read end at a
different depth is cut back to its outer part, or dropped when fewer than 11
bases remain: sequence that appears at the far end of the molecule without the
technical layers around it is insert, such as a conserved gene end that reads
reach from the other side. A mirror at the same depth, or none, is neutral, so
one-sided rapid libraries are unaffected. A candidate that occupies the same
reads as a stronger candidate of its layer is a sequencing variant and is
dropped. Candidates dominated by short
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
gene prefix can be indistinguishable, and a conserved gene end that mirrors at
the same depth in fully symmetric amplicon reads is kept as insert only when
unprimed reads establish the boundary. Short primers, highly degenerate mixtures,
high sequencing error, adapter lengths approaching the sampling-window length,
or insufficient representation can prevent complete recovery. Use a known
primer FASTA or kit preset when available. Inferred sequences use the adapter
role: they trim ends and can excise interior junctions. Use
`--adapter-ends-only` to suppress adapter splitting, or an explicit FASTA with
`primer` in its headers to restrict primer matches to the ends.

`--adapter-report` runs the same discovery and prints the sequences trimming would use,
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
