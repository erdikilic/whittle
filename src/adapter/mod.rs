//! Adapter, primer and barcode trimming.
//!
//! Searches each read window for catalog sequences with sassy, classifies every
//! accepted hit as a terminal trim or an interior excision, re-trims the ends
//! that an excision creates, and returns the kept spans. Presence detection,
//! de novo inference and the built-in catalog live in the submodules.

pub mod catalog;
pub mod detect;
pub mod infer;
pub mod preset;
pub mod resolve;
pub mod search;

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::OnceLock;

use crate::split::Primer;

use search::{
    AmbiguousSearcher, EncodedAdapterBatch, Hit, MAX_TILED_PATTERN_LEN, PlainSearcher, Strands,
    encode_patterns, encoded_pattern_hits, for_each_hit, for_each_hit_in_texts,
    for_each_hit_on_strand, is_plain_acgt, iupac_bases, new_ambiguous_searcher,
    new_overhang_searcher, new_plain_searcher_fwd, new_searcher, new_searcher_fwd,
};

mod budget;
mod hits;
mod index;
mod passes;
mod refine;
pub(crate) use budget::*;
use hits::*;
pub(crate) use index::*;
pub(crate) use passes::*;

thread_local! {
    /// The searchers and per-read buffers of this thread, reused across reads
    /// so `adapter_segments` allocates neither a searcher nor its scratch on
    /// every call. Per-thread state keeps the parallel workflows free of
    /// sharing.
    static STATE: RefCell<ThreadState> = RefCell::new(ThreadState::new());
}

/// What a catalog sequence is. Every role is trimmed at the read ends and
/// splits a read at an interior hit, except a panel barcode of a set that
/// carries its flanks (see `CandidateIndex::splits`); only an adapter is
/// searched again with its end bases masked (`search_residue`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    /// A sequencing adapter.
    Adapter,
    /// A PCR or sequencing primer.
    Primer,
    /// A barcode, barcode flank or barcode construct.
    Barcode,
}

impl Role {
    /// Display label.
    pub fn label(self) -> &'static str {
        match self {
            Role::Adapter => "adapter",
            Role::Primer => "primer",
            Role::Barcode => "barcode",
        }
    }
}

/// One searchable adapter, primer, barcode or flank sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    /// Display name, used in logs and the report.
    pub name: String,
    /// Nucleotide sequence; may contain IUPAC ambiguity codes.
    pub seq: Vec<u8>,
    /// What the sequence is, which decides whether an interior hit splits.
    pub role: Role,
}

/// Resolved adapter-trimming settings for a run.
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    /// Sequences searched for, in configuration order.
    pub adapters: Vec<Adapter>,
    /// End-match tolerance as a fraction of adapter length (`k_end`), also the
    /// per-base cost of a pattern base hanging off a read end.
    pub error_rate: f64,
    /// Bases at each end within which a hit is terminal (trim) rather than
    /// interior (split).
    pub end_size: usize,
    /// Whether interior adapters split the read; `false` is ends-only
    /// (`--adapter-ends-only`).
    pub split: bool,
    /// Shortest segment worth keeping between two excisions. Excisions closer
    /// than this merge into one, since the bases between them would be
    /// discarded by the length filter anyway. Follows `--min-length`.
    pub min_piece: usize,
    /// Exact-seed index for lossless whole-read candidate filtering, built
    /// lazily once presence detection or inference has finalized `adapters`.
    pub(crate) candidate_index: OnceLock<CandidateIndex>,
    /// Whether the reads come from an amplicon library, as resolution judges
    /// from a sample (`detect::primer_share`). An amplicon holds no
    /// marker-gene primer inside it, so an interior hit of a marker primer in
    /// the primer role (`CandidateIndex::paired`) then splits the read by
    /// itself where the whole primer aligns, as a junction partner otherwise
    /// requires. Resolution sets it; `false` keeps the partner requirement.
    pub amplicon: bool,
    /// Per adapter: the index of the split sheet primer it is, or `None`.
    /// Empty when no split sheet is attached (`attach_split`).
    pub(crate) split_of: Vec<Option<usize>>,
    /// Per adapter with a `split_of` primer: the strands on which its hits
    /// read into the insert, as the sheet primer does. `Reversed` for an
    /// entry that is the reverse complement of its primer, `Both` for an
    /// entry that is a primer and the reverse complement of a primer. Empty
    /// when no split sheet is attached.
    pub(crate) split_opens: Vec<Opens>,
}

/// The strands on which a hit of a paired entry reads into the insert, as its
/// primer is synthesized, and the strands on which it reads out of it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Opens {
    /// A hit of the sequence as given reads into the insert; a hit of its
    /// reverse complement reads out of it.
    AsGiven,
    /// A hit of the reverse complement reads into the insert; a hit of the
    /// sequence as given reads out of it.
    Reversed,
    /// A hit on either strand reads into the insert and out of it: the entry
    /// is a primer and the reverse complement of a primer.
    Both,
}

impl Opens {
    /// Returns `AsGiven` when `as_given` and `Reversed` otherwise.
    pub(crate) fn from_given(as_given: bool) -> Self {
        if as_given {
            Opens::AsGiven
        } else {
            Opens::Reversed
        }
    }

    /// Whether a hit on the strand `rc` gives reads into the insert.
    pub(crate) fn reads_in(self, rc: bool) -> bool {
        match self {
            Opens::AsGiven => !rc,
            Opens::Reversed => rc,
            Opens::Both => true,
        }
    }

    /// Whether a hit on the strand `rc` gives reads out of the insert.
    pub(crate) fn reads_out(self, rc: bool) -> bool {
        match self {
            Opens::AsGiven => rc,
            Opens::Reversed => !rc,
            Opens::Both => true,
        }
    }

    /// Returns the sense of an entry that serves the primers of `self` and
    /// of `other`.
    fn merged(self, other: Opens) -> Opens {
        if self == other { self } else { Opens::Both }
    }
}

impl AdapterConfig {
    /// Replaces the adapter set, detaches the split sheet, and discards the
    /// candidate index built for the previous set.
    pub(crate) fn replace_adapters(&mut self, adapters: Vec<Adapter>) {
        self.adapters = adapters;
        self.split_of.clear();
        self.split_opens.clear();
        self.candidate_index = OnceLock::new();
    }

    /// Attaches the primers of a split sheet. Each primer maps to the first
    /// entry whose uppercase sequence equals the primer or its reverse
    /// complement, which keeps its name and role; a primer with no such entry
    /// is appended in the primer role. An entry maps to the first primer that
    /// matches it, and reads into the insert in the orientation of every
    /// primer that matches it (`split_opens`). Leaves `amplicon` as resolution
    /// judged it: a split entry splits a read only at a junction pair
    /// (`search_pairs`). Discards the candidate index built without the split
    /// sheet.
    pub fn attach_split(&mut self, primers: &[Primer]) {
        let mut split_of = vec![None; self.adapters.len()];
        let mut split_opens = vec![Opens::AsGiven; self.adapters.len()];
        for (primer_idx, primer) in primers.iter().enumerate() {
            let seq = primer.seq.to_ascii_uppercase();
            let rc = reverse_complement(&seq);
            let found = self.adapters.iter().enumerate().find_map(|(i, entry)| {
                let entry_seq = entry.seq.to_ascii_uppercase();
                match (entry_seq == seq, entry_seq == rc) {
                    (true, true) => Some((i, Opens::Both)),
                    (true, false) => Some((i, Opens::AsGiven)),
                    (false, true) => Some((i, Opens::Reversed)),
                    (false, false) => None,
                }
            });
            match found {
                Some((i, opens)) => {
                    if split_of[i].is_none() {
                        split_of[i] = Some(primer_idx);
                        split_opens[i] = opens;
                    } else {
                        split_opens[i] = split_opens[i].merged(opens);
                    }
                },
                None => {
                    let opens = if seq == rc {
                        Opens::Both
                    } else {
                        Opens::AsGiven
                    };
                    self.adapters.push(Adapter {
                        name: primer.name.clone(),
                        seq: primer.seq.clone(),
                        role: Role::Primer,
                    });
                    split_of.push(Some(primer_idx));
                    split_opens.push(opens);
                },
            }
        }
        self.split_of = split_of;
        self.split_opens = split_opens;
        self.candidate_index = OnceLock::new();
    }

    /// Returns whether the entry at `adapter_idx` is a split sheet primer.
    pub(crate) fn is_split(&self, adapter_idx: usize) -> bool {
        self.split_of.get(adapter_idx).copied().flatten().is_some()
    }

    /// Sets `amplicon` and discards the candidate index built without it.
    pub(crate) fn set_amplicon(&mut self, amplicon: bool) {
        self.amplicon = amplicon;
        self.candidate_index = OnceLock::new();
    }

    /// Returns whether a barcode sequence matches the read at `[start, end)`
    /// within the error budget: a barcode-role entry of the configured set,
    /// or a catalog barcode with the number of the barcode call `call`. The
    /// span is widened by `BARCODE_SPAN_SLACK` bases on each side.
    pub(crate) fn barcode_span_verified(
        &self,
        seq: &[u8],
        start: usize,
        end: usize,
        call: Option<&[u8]>,
    ) -> bool {
        let lo = start.saturating_sub(BARCODE_SPAN_SLACK);
        let hi = end.saturating_add(BARCODE_SPAN_SLACK).min(seq.len());
        if hi <= lo {
            return false;
        }
        let text = seq[lo..hi].to_ascii_uppercase();
        let number = call.and_then(barcode_number);
        let named = catalog_barcodes().iter().filter(|entry| {
            number.is_some_and(|n| barcode_number(entry.name.as_bytes()) == Some(n))
        });
        let mut searcher = search::new_ambiguous_searcher();
        self.adapters
            .iter()
            .filter(|entry| entry.role == Role::Barcode)
            .chain(named)
            .filter(|entry| entry.seq.len() >= MIN_PATTERN_LEN)
            .any(|entry| {
                let pattern = entry.seq.to_ascii_uppercase();
                let k = edit_budget(self.error_rate, pattern.len());
                !search::hits(&mut searcher, &pattern, &text, k).is_empty()
            })
    }
}

/// Bases by which a recorded barcode span is widened before verification.
const BARCODE_SPAN_SLACK: usize = 3;

/// Every barcode-role catalog entry, built once.
fn catalog_barcodes() -> &'static [Adapter] {
    static ENTRIES: OnceLock<Vec<Adapter>> = OnceLock::new();
    ENTRIES.get_or_init(|| {
        preset::preset(preset::Kit::ALL)
            .into_iter()
            .filter(|entry| entry.role == Role::Barcode)
            .collect()
    })
}

/// Returns the numeric value of the trailing digits of a barcode name, such
/// as `07` in `SQK-NBD114-24_barcode07` or `BC07`.
fn barcode_number(name: &[u8]) -> Option<u32> {
    let digits = name.iter().rev().take_while(|b| b.is_ascii_digit()).count();
    (digits > 0)
        .then(|| {
            std::str::from_utf8(&name[name.len() - digits..])
                .ok()?
                .parse()
                .ok()
        })
        .flatten()
}

/// Returns the IUPAC complement of one base or ambiguity code, case preserved.
/// `S`, `W` and `N` are their own complements; any other byte passes through.
fn complement(base: u8) -> u8 {
    let upper = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    };
    if base.is_ascii_lowercase() {
        upper.to_ascii_lowercase()
    } else {
        upper
    }
}

/// Returns the reverse complement of `seq`, code by code.
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Minimum searchable pattern length: a shorter pattern matches almost anywhere
/// under any error budget and is never searched standalone. Catalog flanks
/// below it are omitted from the catalog for the same reason.
pub const MIN_PATTERN_LEN: usize = 11;

/// Minimum aligned pattern length of a partial hit at a read end. A shorter
/// overlap matches a random read end too often to be evidence of an adapter.
pub const MIN_OVERLAP: usize = 10;

/// Outboard flank length at or below which a hit covered by both end zones is
/// a terminal trim rather than an excision, and at or below which the bases
/// between two excisions are merged into one. A flank this short is adapter
/// residue or junk and is not worth a read of its own.
const FLANK_SLACK: usize = MIN_PATTERN_LEN;

/// Computes the adapter keep segments for `window`: terminal hits within
/// `end_size` of an end trim that end inward, and interior hits (at the
/// stricter `k_mid`) excise and split. Every segment an excision creates is
/// then searched again at its new ends, so the residue next to a junction
/// (a truncated adapter, a primer, a barcode, junk) is trimmed as it would be
/// at a physical read end.
///
/// Under `--adapter-ends-only` (`cfg.split` false) only the two end zones are
/// searched, since no interior hit could be acted on.
///
/// Returns `[start, end)` spans in `window` coordinates.
pub fn adapter_segments(window: &[u8], cfg: &AdapterConfig) -> Vec<(usize, usize)> {
    spans(segments_tallied(window, cfg, None, cfg.min_piece))
}

/// `adapter_segments` with the split primer located at each end of every
/// segment (`Segment::five`, `Segment::three`). Without an attached split
/// sheet (`AdapterConfig::attach_split`) every locus is `None`.
pub fn adapter_segments_annotated(window: &[u8], cfg: &AdapterConfig) -> Vec<Segment> {
    segments_tallied(window, cfg, None, cfg.min_piece)
}

/// Locates segments without merging inserts by the output length threshold.
/// Retained primer bases contribute to the final length filter.
pub(crate) fn adapter_segments_retained(window: &[u8], cfg: &AdapterConfig) -> Vec<Segment> {
    segments_tallied(window, cfg, None, 0)
}

/// `adapter_segments` that also marks in `acted` every adapter whose hit
/// trimmed or excised part of `window`, for presence detection.
pub(crate) fn adapter_segments_tallied(
    window: &[u8],
    cfg: &AdapterConfig,
    acted: &mut [bool],
) -> Vec<(usize, usize)> {
    spans(segments_tallied(window, cfg, Some(acted), cfg.min_piece))
}

/// Returns the `[start, end)` span of each segment.
fn spans(segments: Vec<Segment>) -> Vec<(usize, usize)> {
    segments.into_iter().map(|s| (s.start, s.end)).collect()
}

/// The span of an accepted split primer hit, in window coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Locus {
    /// Hit start.
    pub start: usize,
    /// Hit end, exclusive.
    pub end: usize,
    /// Whether the outer edge of the locus (`start` at the 5' end, `end` at
    /// the 3' end) counts as a read end for rescoring: the hit that located
    /// it hung off the end of the searched span or off the trim boundary it
    /// was searched against, so that edge is a read end or that boundary.
    pub outer_open: bool,
    /// Whether the locus is a whole primer at an outer layer's trim boundary,
    /// located within its anchored budget (`anchored_budgets`), which
    /// rescoring then applies at this end.
    pub boundary: bool,
    /// The whole hit that located the locus and spans exactly `[start,
    /// end)`, or `None` for a locus whose hit overhangs or that an excision
    /// located.
    pub site: Option<PrimerSite>,
}

/// The whole entry hit behind a locus: the read holds the sequence of the
/// entry (its reverse complement when `rc`) over the locus span at `cost`
/// edits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PrimerSite {
    /// Index into `AdapterConfig::adapters`.
    pub entry: usize,
    /// Whether the entry matched as its reverse complement.
    pub rc: bool,
    /// Edit cost of the hit.
    pub cost: usize,
}

/// Whether `locus` keeps the convention of a site: without a site, or with
/// the sequence of its entry in `adapters` (its reverse complement when
/// `rc`) aligning end to end over `[start, end)` of `read` at no more than
/// the site's cost. Rescoring bounds a primer's cost at a locus from its
/// site on that convention (`Scorer::score_at_site`). A read byte matches a
/// code that stands for it after normalization (`normalize_base`), so a
/// byte outside ACGT matches nothing.
pub(crate) fn site_holds(read: &[u8], adapters: &[Adapter], locus: &Locus) -> bool {
    let Some(site) = locus.site else {
        return true;
    };
    let (Some(entry), Some(text)) = (adapters.get(site.entry), read.get(locus.start..locus.end))
    else {
        return false;
    };
    let pattern = if site.rc {
        reverse_complement(&entry.seq)
    } else {
        entry.seq.clone()
    };
    let matches = |code: u8, base: u8| {
        let base = normalize_base(base);
        iupac_bases(code).is_some_and(|bases| bases.contains(&base))
    };
    // `row[j]` is the least cost of aligning the pattern prefix so far with
    // the first `j` bases of the text.
    let mut row: Vec<usize> = (0..=text.len()).collect();
    for (i, &code) in pattern.iter().enumerate() {
        let mut diagonal = row[0];
        row[0] = i + 1;
        for (j, &base) in text.iter().enumerate() {
            let substitution = diagonal + usize::from(!matches(code, base));
            diagonal = row[j + 1];
            row[j + 1] = substitution.min(row[j + 1] + 1).min(row[j] + 1);
        }
    }
    row[text.len()] <= site.cost
}

/// One kept adapter segment with the split primer located at each of its
/// ends, in window coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Segment {
    /// Segment start.
    pub start: usize,
    /// Segment end, exclusive.
    pub end: usize,
    /// The split primer located at `start`: the outermost orientation-valid
    /// terminal hit of its end, or the excised junction hit beside it.
    pub five: Option<Locus>,
    /// The split primer located at `end`, as `five` is at `start`.
    pub three: Option<Locus>,
}

/// A located split primer and the span by which it is matched to a segment
/// end: the locus itself, or, where further hits abutting it moved the trim
/// boundary, the locus together with those hits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Located {
    /// The located primer.
    pub(crate) locus: Locus,
    /// Start of the matched span.
    pub(crate) start: usize,
    /// End of the matched span, exclusive.
    pub(crate) end: usize,
}

impl Located {
    /// Returns `self` moved `offset` bases to the right.
    fn shifted(self, offset: usize) -> Self {
        Located {
            locus: Locus {
                start: self.locus.start + offset,
                end: self.locus.end + offset,
                ..self.locus
            },
            start: self.start + offset,
            end: self.end + offset,
        }
    }
}

impl Segment {
    /// Returns the segment `[start, end)` with the loci of `hits` at its
    /// ends: at `start`, the locus whose matched span ends within
    /// `FLANK_SLACK` bases of it, the one ending last when several do; at
    /// `end`, the locus whose matched span starts within `FLANK_SLACK` bases
    /// of it, the one starting first when several do.
    fn located(start: usize, end: usize, hits: &[Located]) -> Self {
        let near = |pos: usize, at: usize| pos + FLANK_SLACK >= at && pos <= at + FLANK_SLACK;
        Segment {
            start,
            end,
            five: hits
                .iter()
                .filter(|h| near(h.end, start))
                .max_by_key(|h| h.end)
                .map(|h| h.locus),
            three: hits
                .iter()
                .filter(|h| near(h.start, end))
                .min_by_key(|h| h.start)
                .map(|h| h.locus),
        }
    }
}

#[cfg(test)]
mod segment_tests;
