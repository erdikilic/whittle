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

Every sequence has a role: **adapter**, **primer** or **barcode**. Catalog
entries carry their role. A FASTA entry is an adapter unless the word `primer`
or `barcode` appears in its header description (`>27F primer`). Every role is
trimmed at the read ends and splits a read at an interior hit, with these
exceptions:

- The barcodes of a panel (equal-length barcode entries) do not split when the
  set carries barcode flanks or a barcode construct: the flanks at a barcode
  junction split it.
- The universal marker-gene primers (see
  [Marker-gene primers](#marker-gene-primers)) split only beside the rest of
  a junction stack in a selection that includes a genomic kit, since their
  sites lie inside every genomic read through an rRNA
  operon: another excision within 11 bases, or the other end's primer within
  the end zone, across the barcodes between them. A preset made only of
  amplicon kits (`mab114`) promotes them to adapters, which split alone.
  In an amplicon library, as the sampled reads show it (see
  [Amplicon libraries](#amplicon-libraries)), a marker primer in the primer
  role splits alone too, from a mixed preset, a FASTA or discovery, where
  the whole catalog primer aligns over its interior hit.
  Two of these primers also split a read where they lie side by side in the
  orientation of a chimera junction, whether they come from a preset, a
  FASTA or discovery: the first reverse complemented, as the end of one
  molecule holds it, and the second as synthesized, reading into the next
  molecule, the second starting within 11 bases of the end of the first or
  overlapping it by at most 11. Any two of them pair, including two copies
  of one primer. The two sites of a marker gene in a genome lie a gene
  apart and face each other, so they never form such a pair. Each primer
  of a pair may use up to its terminal budget, bounded jointly by the pair
  budget below; a primer's direction is that of the catalog primer it
  matches.
- An interior hit splits only when the whole sequence aligns: a hit whose
  alignment leaves 10 or more bases to clip at its ends is dropped, as at a
  genomic primer site matched by an adapter assembled with its primer.
- A sheet primer of `--split-by` splits a read only at a junction pair or
  beside an interior adapter or barcode, whatever its role and whatever the
  library; see [Primer split](#primer-split).
- The PCS114 UMI does not split: the SSP beside it splits a junction. It
  trims an end where it follows the SSP.
- A primer or barcode splits only reads short enough that its exact matches
  stay within the interior chance bound below; a short flank such as
  `NB_front` (14 bases) splits reads up to a few kilobases.

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
- **Chimera splitting.** An interior adapter, primer or barcode marks a
  junction. The read is split there, the sequence is excised, and both sides
  are kept. Each side is searched
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
bounded over the whole zone. A sequence and its reverse complement count once,
and so do the marker-gene primers of one binding site (see below), such as
the 16S forward mix and its variants: they count as the one of them with the
most chance matches, at one budget, the lowest of theirs, in the terminal,
interior and pair bounds alike, so the variants of a site leave the budgets
of every other entry where one primer would. Each primer of a site also
keeps no more than the budget it gets with every variant counted as a
sequence of its own.
Every catalog entry keeps its per-pattern budget when anchored at the read end,
a marker-gene primer the lowest of its site, and the catalog adapters and
24-base ONT barcodes keep it anywhere in the end zone.

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

The marker-gene primers of a pair share a pair budget. For each read-length
class, it is the largest edit count per primer, at most its terminal budget,
under which the expected chance pairs of the set in a read of that length
stay within `1e-4`: any primer of the set followed by any primer at one of
the 23 admitted offsets, so the chance rate of a pair at one position is the
square of the summed chance rates of the primers. The largest contributors
lose an edit first. Exact pairs are always accepted.
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

With `--adapter-discover`, discovery starts behind every preset entry, and
detection over the discovery sample then narrows the preset's barcode entries
in the set trimmed against. The barcode constructs and the flanks between a
barcode and the insert that the reads do not carry are dropped, since a
barcode panel leaves its junctions to them. When none of them remains, the
panel splits its own junctions, and its absent barcodes and the other absent
barcode flanks are dropped as well, since a splitting sequence takes part in
the interior chance bound of the set. The other adapters and primers of the
preset are kept. A FASTA disables this narrowing as well.

Adapter samples are bounded by 256 MiB of retained payload and 64 Mi bases as
well as `--adapter-sample-reads`. The last sampled record is retained whole;
the actual sample size is reported, and all sampled reads are processed.

### Amplicon libraries

When the set trimmed against holds a marker-gene primer in the primer role
(a preset that joins `mab114` to a genomic kit such as `mab114,lsk114`, `ont`
or `all`, a FASTA primer, or a discovered primer), the sample that presence
detection or discovery reads also decides whether the library is an amplicon
library. At most 4000 reads spread evenly over the sample are examined.

Every molecule of an amplicon library starts with one of its primers, so a
read keeps one at its 5' end unless it is a fragment, and the reads share the
lengths of a few targets. The primers of the library are the catalog marker
primers (see [Marker-gene primers](#marker-gene-primers)), the primer-role
entries of a preset
or FASTA that are marker primers or no catalog sequence, and every discovered
sequence except the members of a barcode layer and the sequences whose best
catalog match is another catalog entry (an adapter, barcode, flank or the PCR
handle of another kit). They may form any number of families, such as a 16S
pair and a gene-specific pair outside the catalog, and count together. A
read end is opened when a whole hit of a primer, on either strand within its
terminal budget, lies at the boundary that the terminal trims of the other
sequences of the set leave there:

- A marker primer opens the end when its hit starts within 11 bases of that
  boundary or before it.
- Any other primer opens the end when its hit starts within 11 bases of the
  boundary or before it and ends within 11 bases of it or after it. A primer
  hit that starts within 11 bases of the end of the previous one continues
  the stack, and the innermost primer opens the insert. Such a primer counts
  only when the inserts behind it start alike: one 8-base word, read at any
  of the first five offsets behind the primer and not dominated by one base,
  starts at least 3% of the ends it opens and at least 10 of them. The start
  of the targets of an amplicon recurs, while the random fragments behind a
  technical sequence of a genomic or cDNA library start anywhere.

The library is an amplicon library when primers open an end of at least 60%
of the examined reads and the median read shares its length, within 5%,
with at least 15% of them; both shares are logged. A genomic or cDNA library
fails both: its reads start at a primer site only where they happen to start
inside an rRNA operon, and its fragments spread over a wide range of
lengths, each sharing its length with a few percent of the others.

An amplicon holds no primer inside it, so in an amplicon library a marker
primer inside a read is a chimera junction, and an interior hit of a marker
primer in the primer role splits the read by itself, as under the `mab114`
preset alone. The hit must pass its interior budget, and the whole catalog
primer must align over it on the same strand within the primer's own
interior budget, leaving fewer than 10 bases to clip; the excision covers
both alignments. A primer form cut short on its outer side, as discovery
assembles it from eroded read ends, therefore splits only where the whole
primer is present. Primers outside the catalog split alone at a whole
interior hit in every library, within the interior chance bound, as the
other roles do. Outside an amplicon library the partner requirement above
stays, and `--adapter-report` is unaffected.

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
past that point belongs to the next layer. A layer also ends where its path
divides into two continuations that each keep at least a quarter of its
support for 24 k-mers, when each continuation recurs reverse complemented at
the other read end: a cDNA adapter core shared by both strands divides into
the strand-switching primer on one strand and the oligo(dT) primer on the
other, and each primer is found at the other end of the opposite strand. A
degenerate primer base opens a branch that rejoins the path within 16
k-mers, and variants of a conserved insert do not recur at the other end, so
neither divides a layer. A homopolymer run of at least 8 bases at the inner
end of a candidate, such as the poly(A) tail behind an oligo(dT) primer,
varies in length between reads and is left to the insert. A variable layer is resolved from
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

A random tag, such as a UMI, differs between reads and has no recurrent
k-mers, so it is read from base composition instead. Behind every known or
discovered layer other than a barcode, the bases that follow the layer in the
windows that hold it are tallied one position at a time, each round
realigning the windows to the layer and the tag so far. A position belongs to
the tag when one base holds the midpoint between its composition share and
one, or when a base falls fourfold below its composition share; it is written
with the ambiguity code of the bases present, such as `V` for a position that
excludes T. The tag ends with the run of fixed bases after its last
degenerate position; a further run, such as the G run a template switch adds,
varies in length between reads and is left to the insert. A tag needs 16
positions, at least 1% of the windows and 20 windows, and random degenerate
positions: no combination of their bases recurs in 20 windows, and at least
half of the windows hold distinct combinations. A barcode panel, the
conserved starts of an amplicon's species or two alternative primers
concentrate in a few combinations, and a genomic insert holds every base near
its composition share, whatever the genome's composition. A tag takes the
primer role: it trims read ends and, with its chance exact matches above the
interior bound, does not split reads. A tag that a known sequence already
describes is not reported.

Discovery separates primers from the amplicon they bind in two ways. The
far end of the molecule is one: a conserved gene start that reads reach from
the other side appears there without the primer stack around it and is left
in the read, and a primer whose reverse complement is absent from the far
end, as in rapid amplicon libraries, is trimmed. When reads run through the
far primer, the primer and the conserved start both recur at both ends, and
the variation of the templates tells them apart. Related templates share a
conserved start, but its variable bases follow the template: at such a
column, the reads of each of the two most frequent bases carry k-mers
further inside the insert, from 16 to 116 bases behind the candidate, that
at least 90% of the reads holding them share with that base, in at least a
quarter of that base's reads, and the reads of the two bases continue alike,
their most frequent bases agreeing in at least two thirds of the 16 columns
that follow. The degenerate base of a primer and a sequencing error are
carried by the reads of every template, and alternative technical sequences
at one position, such as the primers of the two strands behind a shared
adapter or the members of a barcode panel, continue differently. The second
base must hold at least 20 reads and one in 20 of the column.

A candidate is cut at the first of its columns that follows the template,
and dropped when that column lies within its first 11 bases or when one of
the 16 columns outboard of it follows the template: the candidate then lies
in the insert. A candidate whose first column behind it follows the
template, or that holds the inner end of a catalog marker primer (see
[Marker-gene primers](#marker-gene-primers)), ends at the insert boundary:
when accepted it
closes its read end, as a primer reconstructed from unprimed read starts
does, and its reverse complement marks the insert boundary at the other read
end, where the far primer of the molecule is read. A candidate that reads
past the inner end of such a sequence, in the reads holding both or with at
least 7 of its outer bases repeating that end, as a read eroded into the
primer at its physical end holds it, is cut there. A candidate that ends
inside such a primer in most of the reads holding both, as an outer layer
assembled together with the first bases of the primer does, is continued to
the primer's inner end. Where a candidate that follows the template from its
outer bases lies within 11 bases of the boundary, no other layer is taken at
that end. Behind the sequences of a preset or FASTA, the reads they trimmed
are tested column by column behind the boundary, as a whole and as two
families split by their most frequent far k-mer, as the two ends of a gene
read from the two strands are; when the template begins within 11 bases of
the boundary, discovery adds no layer at that end.

Conserved bases before the first variable column of the insert cannot be
told apart from the primer: a primer outside the catalog whose templates all
start with the same bases, or a library of one template, keeps them with the
primer. When the templates vary too rarely for the counts above, discovery
leaves the primers of one-sided libraries and, in libraries barcoded at both
ends, trims the conserved start with them. Supply the primers with a kit
preset or `--adapter-fasta` for such libraries.

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
one-sided rapid libraries are unaffected. A layer that divides into strand
primers, and the primer layers behind it, are exempt: their mirrors lie at a
depth set by how much outer adapter each end keeps, which differs between
the 5' and 3' ends of nanopore reads. A candidate that occupies the same
reads as a stronger candidate of its layer is a sequencing variant and is
dropped. The two ends of a library can read one adapter in forms that end at
different bases on the insert side, as the two strands of a Y adapter do.
When the other end's candidate of an accepted family, of at least 16 bases
and half the accepted sequence, aligns inside it short of its inner end, that
end's form is kept beside it: the candidate, extended outward with the
accepted sequence's bases for as long as most of the other end's windows, and
at least 20, carry them, since erosion at the physical end shortens it there.
An interior copy of either form then aligns without paying for bases it
lacks, and a junction holding the shorter form splits as one holding the
longer does. The other end can also read a family in a form of its own that
the accepted sequence does not trim, as a hairpin adapter read on the same
strand at both ends is: its first bases face the insert at the 3' end, where
a sequence assembled at the eroded 5' end lacks them. A candidate of the other
end of at least 16 bases in the family is a layer of that end when a trim with
the accepted sequence alone stops more than 3 bases short of the candidate's
inner edge in at least half of that end's windows holding it; discovery then
continues behind it at that end. Candidates dominated by short
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

Without reads exposing an insert boundary or variants that follow the
template, an unknown primer and a conserved gene prefix can be
indistinguishable, and a conserved gene end that mirrors at
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
| `pcb114`, `pcs114` | cDNA-PCR sequencing and barcoding kits V14 | Ligation adapter, PCS primers, PCS114 UMI, PCR flanks, 24 barcodes |
| `rpb114` | Rapid PCR barcoding kit V14 | Rapid adapter, PCR primers, rapid ligation flank, 24 barcodes |
| `mab114` | Microbial amplicon barcoding kit V14 | Ligation and rapid adapters, the kit's seven 16S and five ITS primers, MAB flanks, 24 barcodes |
| `rna004` | Direct RNA sequencing kit | RNA adapter |
| `pacbio` | PacBio SMRTbell | SMRTbell adapter, C2 sequencing primer |
| `ont` | Every ONT kit above | |
| `all` | Every kit | |

The catalog is assembled from vendor-published sources: dorado's
`adapter_primer_kits.cpp` and `barcode_kits.cpp`, Porechop's `adapters.py`,
qcat's kit definitions, ONT's primer sheet for SQK-MAB114.24, and the NCBI
UniVec records of the PacBio sequences.
Every entry names its kits, so a kit preset searches only what its library can
contain; a union searches everything, at some cost in speed and precision.
Flanks shorter than 11 bp are excluded, since a pattern that short matches
almost anywhere.

### Marker-gene primers

The `mab114` preset lists the primers of the kit as ONT ships them, seven
for 16S and five for ITS, under their ONT names:

| Primer | Sequence |
|---|---|
| 16S_mix_F | `AGRGTTYGATYMTGGCTCAG` |
| 16S_mix_R | `SGGYTACCTTGTTACGACTT` |
| 16S_Bor_F | `AGAGTTTGATCCTGGCTTAG` |
| 16S_Bor_R | `CGGCTACCTTGTTACGACTT` |
| 16S_Chl_F | `AGAATTTGATCTTRGTTCAG` |
| 16S_Chl_R | `GGGCTACCTTGTTACGACTT` |
| 16S_Ent_F | `AGAGTTTGATCATGGCTCAG` |
| ITS1 | `TCCGTAGGTGAACCTGCGG` |
| ITS1_Fus | `TCCGTTGGTGAACCAGCGG` |
| ITS1_Mal | `TCTGTAGGTGAACCTGCAG` |
| ITS4 | `TCCTCCGCTTATTGATATGC` |
| ITS4_Pyt | `TCCTCCGCTTATTAATATGC` |

The universal marker-gene primers are these twelve and two community
primers outside the kit, fungal ITS1F (`CTTGGTCATTTAGAGGAAGTAA`) and the
22-base 16S 1492R (`TACGGYTACCTTGTTACGACTT`). The classic 16S 27F
(`AGAGTTTGATYMTGGCTCAG`) is an instance of 16S_mix_F. A sequence from a
preset, a FASTA or discovery is a marker primer when it matches one of
them on either strand within its edit budget and differs from it in length
by fewer than 10 bases; a shorter sequence must reach an end of the primer.
Among primers that match at the same cost, the one nearest in length
decides, then the first in the order above. The binding sites group them:
16S forward (16S_mix_F, 16S_Bor_F, 16S_Chl_F, 16S_Ent_F), 16S reverse
(16S_mix_R, 16S_Bor_R, 16S_Chl_R, 1492R), ITS1 (ITS1, ITS1_Fus, ITS1_Mal),
ITS4 (ITS4, ITS4_Pyt), and ITS1F.

## Primer split

`--split-by` ([cli.md](cli.md#primer-split)) hands the adapter engine a
sheet of targets, each a list of forward and a list of reverse primers.
Each sheet primer joins the search set in the primer
role: it is matched by sequence, or by reverse complement, to an existing
preset or FASTA entry and reused under that entry's name, or appended as a
new entry when no match exists. An entry that matches a primer and the
reverse complement of another (a target whose reverse primer is the reverse
complement of its forward primer) reads into the insert in both
orientations. Sheet primers are attached to the search set
after presence detection or discovery has narrowed it, so they are always
searched and can never be narrowed away. `--split-by` leaves the
amplicon-library judgement as the sampled reads gave it: it still decides
whether a catalog marker primer outside the sheet splits alone, and it is
not forced on.

A sheet primer inside a read splits it only at a junction: a sheet primer in
its closing orientation (reverse complemented, as it ends an amplicon)
followed by a sheet primer in its opening orientation (as the next amplicon
starts), the second starting within 11 bases of the end of the first or
overlapping it by at most 11, or a sheet primer beside an interior adapter
or barcode. A lone sheet primer inside a read does not split it: a nested or
overlapping panel carries the primer sites of some targets inside the
amplicons of others, and a full-length amplicon stays one read. Overlapping
primers and two primers facing each other form no junction either.

At the read ends, a sheet primer acts only in the orientation a primer has
there: as given (reading into the insert) at the 5' end, reverse
complemented at the 3' end. A hit in the other orientation neither trims nor
locates a primer. Of the valid hits at one end, the outermost, nearest the
read end or the trimmed adapter, is the located primer, and trimming stops
at its inner edge: a sheet primer hit further inside the insert moves the
trim only when it overlaps the outermost one or starts within 11 bases of
it. Adapters, barcodes and other primers trim by their own rules. A valid
sheet primer hit that lies within or overlaps the trim of another entry at
that end (a catalog primer with other degenerate bases, or an entry that
fuses an adapter with the primer) is located there, and the trim stays
where that entry set it. When the terminal
search of an end leaves that end trimmed by another layer (an adapter,
barcode or other primer) and holds no sheet primer hit there, the bases just
inside the trim boundary it set are searched once more with the boundary as
a read end, so a primer that lost its first bases at the adapter junction
is located under the same partial-hit rules as a primer cut short by the
read end. When that end still holds no sheet primer hit, a whole sheet
primer hit in the orientation valid for that end whose outer edge lies at
most 3 bases outboard of the boundary or at most 11 bases inboard of it is
located there within the primer's anchored edit budget: the boundary fixes
its position, so the chance of a random match there is that of a few
positions, not of the whole end zone. The bases outboard of the boundary
are the flank the other layer matched, not random bases, so a primer
extends into them only by the few bases that layer's alignment may claim
from it. The anchored budget is the largest edit count, at most the
error-rate ceiling, whose expected chance matches over the 15 start
positions beside each of the two boundaries of a read, on the one strand
valid there, stay within the same per-read target as the terminal budgets,
for each sheet primer alone and for the sheet's primers together, with the
variants of one marker-primer site counted once as for the terminal
budgets. The search is measured from the boundary, not from the read end,
so a primer that the adapter, flank and barcode push past the end zone, or
that both end zones cover on a short read, is still the terminal primer of
its end.

Locating a split primer finds one hit per end, but the classifier needs to
know how every sheet primer, not only the one located, would score there.
Rescoring aligns each sheet primer against the located locus, widened into
a window sized to the longest sheet primer and clamped to the read, in the
orientation valid for that end: as given at the 5' end, reverse
complemented at the 3'. A primer scores only within its own terminal edit
budget, or, at a locus found whole at a trim boundary as above, a whole
primer hit scores within its anchored edit budget. A primer may hang off
the read end, and also off the outer edge of the located primer when that
primer was found hanging off a read end or a trim boundary: that edge then
counts as a read end for every sheet primer, so a primer the engine
located partially still scores there. A located primer that aligned whole
keeps a window widened on both sides. The two ends' scores go to the
classifier. For a target and strand, an end counts at the cheapest scored
primer of the list the strand puts there, so the primer variants of one
target never compete with each other. A primer of another key closer than
`--split-lead` edits to a scored primer stands in for it at one edit more.
An end that scored primers of other keys only, each within `--split-lead`
of the penalty, may hold the target's primer at that penalty: the lowest
budget rescoring applied at that end (terminal, or anchored at a locus
found whole at a trim boundary) plus one, the same for every target. A
primer that leads the penalty by the lead rules out every target that does
not list it there. The classifier assigns the cheapest strand-consistent
target only when it beats every other target by at least `--split-lead`
edits; a smaller margin, or a tie, is ambiguous rather than a guess. See
[cli.md](cli.md#classification) for the calls.
