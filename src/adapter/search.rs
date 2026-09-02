//! Approximate sequence search over sassy's DNA and IUPAC profiles.
//!
//! The DNA profile is faster and panics on any byte outside A/C/G/T; the IUPAC
//! profile accepts ambiguity codes and is the only one with a batched pattern
//! search. Each entry point states which profile it uses and why.

use sassy::profiles::{Dna, Iupac, Profile};
use sassy::{CachedRev, EncodedPatterns, RcSearchAble, Searcher};

/// The fast searcher: sassy's DNA profile, for all-ACGT patterns against
/// all-ACGT text.
///
/// It panics during traceback on any other byte, in the pattern or the text, so
/// every use is gated on `is_plain_acgt` for both. The two profiles are kept
/// side by side instead of standardizing on IUPAC because, on a narrowed adapter
/// set where most patterns are searched one at a time rather than batched across
/// SIMD lanes, IUPAC costs several times more.
pub type PlainSearcher = Searcher<Dna>;

/// The general searcher: sassy's IUPAC profile, which handles ambiguity codes in
/// the pattern (a degenerate primer) and in the text (an `N` in a read).
pub type AmbiguousSearcher = Searcher<Iupac>;

/// The batched, pattern-parallel searcher. Sassy implements `search_patterns`
/// only for IUPAC, so the batched path is always ambiguity-safe.
pub type BatchedAdapterSearcher = Searcher<Iupac>;

/// Equal-length patterns encoded once for the tiled search, one encoding per
/// strand. Built by `encode_patterns`, searched by `encoded_pattern_hits`.
///
/// The reverse strand is searched as sassy's own two-strand search does it:
/// the complemented patterns over the reversed text, rather than the
/// reverse-complemented patterns over the forward text. The two are the same
/// alignments, but the local-minimum rule picks the rightmost end position of
/// a flat cost run and the traceback then fixes the start, so searching the
/// reversed text puts the tie on the same read position as `pattern_hits` and
/// `hits` and keeps every span identical.
#[derive(Debug, Clone)]
pub struct EncodedAdapterBatch {
    /// The patterns as given, for the forward text.
    forward: EncodedPatterns<Iupac>,
    /// The complemented patterns, for the reversed text.
    complement: EncodedPatterns<Iupac>,
}

/// Longest pattern the tiled search encodes: one pattern per 64-bit limb.
pub const MAX_TILED_PATTERN_LEN: usize = 64;

/// Returns whether every byte is an uppercase or lowercase A/C/G/T, so that
/// `PlainSearcher` can be used without its traceback panic.
pub fn is_plain_acgt(seq: &[u8]) -> bool {
    seq.iter()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

/// Returns the plain bases an IUPAC code stands for, or `None` if the byte is
/// not a nucleotide code.
///
/// `U` is absent: callers fold it to `T` before it reaches here, because sassy
/// treats `U` as a fifth base that matches nothing in a DNA read.
pub fn iupac_bases(code: u8) -> Option<&'static [u8]> {
    Some(match code.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => return None,
    })
}

/// Returns how many of the four bases an IUPAC code stands for, or `None` if the
/// byte is not a nucleotide code. See `iupac_bases` for the `U` rule.
pub fn iupac_degeneracy(code: u8) -> Option<u8> {
    iupac_bases(code).map(|bases| bases.len() as u8)
}

/// A text window with its reversal, both borrowed from per-read buffers, so a
/// two-strand search copies nothing. `reversed[i]` is
/// `forward[forward.len() - 1 - i]`.
#[derive(Debug, Clone, Copy)]
pub struct Strands<'a> {
    /// The window as read.
    pub forward: &'a [u8],
    /// The window reversed.
    pub reversed: &'a [u8],
}

impl RcSearchAble for Strands<'_> {
    fn text(&self) -> impl AsRef<[u8]> {
        self.forward
    }

    fn rev_text(&self) -> impl AsRef<[u8]> {
        self.reversed
    }
}

/// One approximate match of a pattern in the text. Strand is not exposed: a
/// reverse-complement hit occupies the same text span, which is all the trimmer
/// needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    /// Start of the span in the text, inclusive.
    pub start: usize,
    /// End of the span in the text, exclusive.
    pub end: usize,
    /// Edit distance of the match.
    pub cost: usize,
}

/// Returns a fresh DNA-profile searcher over both strands.
pub fn new_searcher() -> PlainSearcher {
    Searcher::<Dna>::new_rc()
}

/// Returns a fresh IUPAC-profile searcher over both strands.
pub fn new_ambiguous_searcher() -> AmbiguousSearcher {
    Searcher::<Iupac>::new_rc()
}

/// Returns a fresh IUPAC-profile searcher over the forward strand only. Used by
/// the inference k-mer recount, where each read-end window is already
/// strand-oriented and reverse-complement hits would inflate the per-window
/// presence count.
pub fn new_searcher_fwd() -> AmbiguousSearcher {
    Searcher::<Iupac>::new_fwd()
}

/// Returns a fresh pattern-batched searcher over both strands.
pub fn new_batched_searcher() -> BatchedAdapterSearcher {
    Searcher::<Iupac>::new_rc()
}

/// Encodes equal-length patterns for `encoded_pattern_hits`, or returns `None`
/// past `MAX_TILED_PATTERN_LEN`, where the caller keeps `pattern_hits`. The
/// encoding holds the bit profiles of every pattern on both strands, which
/// `pattern_hits` rebuilds on each call.
pub fn encode_patterns(patterns: &[Vec<u8>]) -> Option<EncodedAdapterBatch> {
    let len = patterns.first().map_or(0, Vec::len);
    if !(1..=MAX_TILED_PATTERN_LEN).contains(&len) {
        return None;
    }
    // A forward-only searcher encodes the given strand alone; the reverse
    // strand is its own encoding.
    let mut encoder = Searcher::<Iupac>::new_fwd();
    let complements: Vec<Vec<u8>> = patterns.iter().map(|p| Iupac::complement(p)).collect();
    Some(EncodedAdapterBatch {
        forward: encoder.encode_patterns(patterns),
        complement: encoder.encode_patterns(&complements),
    })
}

/// Searches a pre-encoded batch over both strands of `text`, one pattern per
/// SIMD lane, and calls `accept` with each hit's pattern index, text span and
/// cost. `reversed` is `text` reversed, which the caller keeps per read so the
/// reverse strand needs no copy. Hits are the rightmost local minima within
/// `k`, exactly as `pattern_hits` returns them.
pub fn encoded_pattern_hits(
    searcher: &mut BatchedAdapterSearcher,
    encoded: &EncodedAdapterBatch,
    text: &[u8],
    reversed: &[u8],
    k: usize,
    mut accept: impl FnMut(usize, usize, usize, usize),
) {
    debug_assert_eq!(text.len(), reversed.len());
    for m in searcher.search_encoded_patterns(&encoded.forward, text, k) {
        accept(m.pattern_idx, m.text_start, m.text_end, m.cost as usize);
    }
    let n = text.len();
    for m in searcher.search_encoded_patterns(&encoded.complement, reversed, k) {
        accept(
            m.pattern_idx,
            n - m.text_end,
            n - m.text_start,
            m.cost as usize,
        );
    }
}

/// Searches equal-length patterns together, packing them across SIMD lanes.
/// This retains `search`'s reverse-text semantics while avoiding the repeated
/// short-text setup of calling `search` once per pattern.
pub fn pattern_hits(
    searcher: &mut BatchedAdapterSearcher,
    patterns: &[Vec<u8>],
    text: &[u8],
    k: usize,
) -> Vec<sassy::Match> {
    // `search_patterns` processes SIMD-sized chunks internally; the cached
    // reversal keeps sassy from rebuilding the reversed text once per chunk.
    let cached_text = CachedRev::new(text, true);
    searcher.search_patterns(patterns, &cached_text, k)
}

/// Calls `accept` with every match of `pattern` in `text` within `k` edits,
/// as a text span. The strands searched depend on how `searcher` was built:
/// `new_searcher` matches both strands, `new_searcher_fwd` the forward strand
/// only. A plain slice reverses itself on each call; `Strands` borrows a
/// reversal the caller keeps. Reuses `searcher`'s internal buffers across
/// calls.
pub fn for_each_hit<P: Profile, T: RcSearchAble + ?Sized>(
    searcher: &mut Searcher<P>,
    pattern: &[u8],
    text: &T,
    k: usize,
    mut accept: impl FnMut(Hit),
) {
    for m in searcher.search(pattern, text, k) {
        accept(Hit {
            start: m.text_start,
            end: m.text_end,
            // Sassy's `Match::cost` is `pa_types::Cost`, an `i32` signed for other
            // algorithms in that crate; a returned match is within the
            // non-negative `k` budget, so the cast is lossless.
            cost: m.cost as usize,
        });
    }
}

/// Returns all matches of `pattern` in `text` within `k` edits, as text spans.
/// See `for_each_hit`.
pub fn hits<P: Profile>(
    searcher: &mut Searcher<P>,
    pattern: &[u8],
    text: &[u8],
    k: usize,
) -> Vec<Hit> {
    let mut out = Vec::new();
    for_each_hit(searcher, pattern, text, k, |hit| out.push(hit));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // `new_rc()` returns two same-span hits (forward and reverse complement) for
    // a reverse-complement-palindromic pattern. Count-based tests use a
    // non-palindromic pattern whose reverse complement is absent from the text,
    // so exactly one hit is returned. `adapter_segments` is unaffected: it
    // deduplicates terminal hits via max/min and merges interior ones.

    /// An exact forward occurrence is reported once with cost 0.
    #[test]
    fn exact_forward_match() {
        let mut s = new_searcher();
        // revcomp(AAAACCCCGGGG) = CCCCGGGGTTTT is absent from the text, so there
        // is one hit.
        let h = hits(&mut s, b"AAAACCCCGGGG", b"TTAAAACCCCGGGGTT", 0);
        assert_eq!(h.len(), 1);
        assert_eq!((h[0].start, h[0].end, h[0].cost), (2, 14, 0));
    }

    /// A reverse-complement occurrence is found by the both-strand searcher.
    #[test]
    fn finds_reverse_complement() {
        // Pattern AAAACCCC has reverse complement GGGGTTTT, which is embedded in
        // the text.
        let mut s = new_searcher();
        let h = hits(&mut s, b"AAAACCCC", b"TTGGGGTTTTAA", 0);
        assert_eq!(h.len(), 1);
        assert_eq!((h[0].start, h[0].end), (2, 10));
    }

    /// The forward-only searcher skips a reverse-complement occurrence and finds
    /// a forward one.
    #[test]
    fn forward_searcher_ignores_reverse_complement() {
        let mut s = new_searcher_fwd();
        // revcomp(AAAACCCC) = GGGGTTTT is in the text; forward-only skips it.
        assert_eq!(hits(&mut s, b"AAAACCCC", b"TTGGGGTTTTAA", 0).len(), 0);
        // The forward pattern is present and found.
        assert_eq!(hits(&mut s, b"AAAACCCC", b"TTAAAACCCCTT", 0).len(), 1);
    }

    /// Every IUPAC code maps to its bases in either case, and other bytes map
    /// to `None`.
    #[test]
    fn iupac_bases_cover_the_alphabet_and_reject_the_rest() {
        for (code, expected) in [
            (b'A', &b"A"[..]),
            (b'R', b"AG"),
            (b'y', b"CT"),
            (b'B', b"CGT"),
            (b'N', b"ACGT"),
        ] {
            assert_eq!(iupac_bases(code), Some(expected), "Code {}", code as char);
            assert_eq!(iupac_degeneracy(code), Some(expected.len() as u8));
        }
        for code in [b'U', b'X', b'.', b'-', b'0'] {
            assert_eq!(
                iupac_bases(code),
                None,
                "Code {} is not a nucleotide",
                code as char
            );
        }
    }

    /// One substitution is found at `k` 1 and not at `k` 0.
    #[test]
    fn tolerates_one_mismatch_within_budget() {
        let mut s = new_searcher();
        // One substitution (position 5, C to A) in AAAACCCCGGGG; the reverse
        // complement is absent from the text.
        assert_eq!(
            hits(&mut s, b"AAAACCCCGGGG", b"TTAAAACACCGGGGTT", 1).len(),
            1
        );
        assert_eq!(
            hits(&mut s, b"AAAACCCCGGGG", b"TTAAAACACCGGGGTT", 0).len(),
            0
        );
    }
}
