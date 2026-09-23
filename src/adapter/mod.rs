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

use search::{
    AmbiguousSearcher, EncodedAdapterBatch, Hit, MAX_TILED_PATTERN_LEN, PlainSearcher, Strands,
    encode_patterns, encoded_pattern_hits, for_each_hit, for_each_hit_in_texts, is_plain_acgt,
    iupac_bases, new_ambiguous_searcher, new_overhang_searcher, new_searcher,
};

mod budget;
mod hits;
mod index;
mod passes;
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

/// What a catalog sequence is, which decides what a hit may do. Every role is
/// trimmed at the read ends; only an adapter splits a read at an interior hit,
/// since a primer or barcode inside a read is part of the molecule as often as
/// it is a chimera signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    /// A sequencing adapter: trimmed at the ends and excised in the interior.
    Adapter,
    /// A PCR or sequencing primer: trimmed at the ends only.
    Primer,
    /// A barcode or barcode flank: trimmed at the ends only.
    Barcode,
}

impl Role {
    /// Whether an interior hit of this role excises and splits the read.
    pub fn splits(self) -> bool {
        matches!(self, Role::Adapter)
    }

    /// Whether a partial hit hanging off a read end is accepted for this role.
    /// A barcode sits between its flanks, so a partial barcode at a read end is
    /// flank residue that the flank hit already trims.
    fn overhangs(self) -> bool {
        !matches!(self, Role::Barcode)
    }

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
}

impl AdapterConfig {
    /// Replaces the adapter set and discards the candidate index built for the
    /// previous set.
    pub(crate) fn replace_adapters(&mut self, adapters: Vec<Adapter>) {
        self.adapters = adapters;
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
    segments_tallied(window, cfg, None)
}

/// `adapter_segments` that also marks in `acted` every adapter whose hit
/// trimmed or excised part of `window`, for presence detection.
pub(crate) fn adapter_segments_tallied(
    window: &[u8],
    cfg: &AdapterConfig,
    acted: &mut [bool],
) -> Vec<(usize, usize)> {
    segments_tallied(window, cfg, Some(acted))
}

#[cfg(test)]
mod segment_tests;
