//! Approximate sequence search over sassy's DNA and IUPAC profiles.
//!
//! The DNA profile is faster and panics on any byte outside A/C/G/T; the IUPAC
//! profile accepts ambiguity codes and is the only one with a tiled pattern
//! search. Each entry point states which profile it uses and why.

use sassy::profiles::{Dna, Iupac, Profile};
use sassy::{EncodedPatterns, RcSearchAble, Searcher};

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
/// the pattern (a degenerate primer) and in the text (an `N` in a read). Sassy
/// implements the tiled pattern search for this profile only, so it also
/// serves `encoded_pattern_hits`.
pub type AmbiguousSearcher = Searcher<Iupac>;

/// Equal-length patterns encoded once for the tiled search, one encoding per
/// strand. Built by `encode_patterns`, searched by `encoded_pattern_hits`.
///
/// The reverse strand is searched as sassy's own two-strand search does it:
/// the complemented patterns over the reversed text, rather than the
/// reverse-complemented patterns over the forward text. The two are the same
/// alignments, but the local-minimum rule picks the rightmost end position of
/// a flat cost run and the traceback then fixes the start, so searching the
/// reversed text puts the tie on the same read position as `hits` and keeps
/// every span identical.
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
/// needs. The overhang fields are in text orientation and are zero unless the
/// searcher was built with an overhang cost.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    /// Start of the span in the text, inclusive.
    pub start: usize,
    /// End of the span in the text, exclusive.
    pub end: usize,
    /// Alignment cost: the edit distance of the aligned part plus the overhang
    /// cost, `floor(alpha * overhang)` per overhanging side.
    pub cost: usize,
    /// Pattern bases hanging off the text start.
    pub left_overhang: usize,
    /// Pattern bases hanging off the text end.
    pub right_overhang: usize,
    /// Whether the reverse complement of the pattern matched.
    pub reverse: bool,
}

impl Hit {
    /// Builds a hit from a sassy match of a `pattern_len`-base pattern. A
    /// reverse-complement match aligns the pattern end at the text start, so
    /// its pattern-side overhangs swap sides in text orientation.
    fn from_match(m: &sassy::Match, pattern_len: usize) -> Self {
        let (head, tail) = (m.pattern_start, pattern_len - m.pattern_end);
        let (left_overhang, right_overhang) = match m.strand {
            sassy::Strand::Fwd => (head, tail),
            sassy::Strand::Rc => (tail, head),
        };
        Hit {
            start: m.text_start,
            end: m.text_end,
            // Sassy's `Match::cost` is `pa_types::Cost`, an `i32` signed for other
            // algorithms in that crate; a returned match is within the
            // non-negative `k` budget, so the cast is lossless.
            cost: m.cost as usize,
            left_overhang,
            right_overhang,
            reverse: m.strand == sassy::Strand::Rc,
        }
    }
}

/// Returns a fresh DNA-profile searcher over both strands.
pub fn new_searcher() -> PlainSearcher {
    Searcher::<Dna>::new_rc()
}

/// Returns a fresh IUPAC-profile searcher over both strands.
pub fn new_ambiguous_searcher() -> AmbiguousSearcher {
    Searcher::<Iupac>::new_rc()
}

/// Returns a fresh IUPAC-profile searcher over both strands that also reports
/// partial matches hanging off either text end, each overhanging base costing
/// `alpha` (0 to 1). Sassy implements overhang alignment for the IUPAC profile
/// only.
pub fn new_overhang_searcher(alpha: f32) -> AmbiguousSearcher {
    Searcher::<Iupac>::new_rc_with_overhang(alpha.clamp(0.0, 1.0))
}

/// Returns a fresh IUPAC-profile searcher over the forward strand only. Used by
/// the inference k-mer recount, where each read-end window is already
/// strand-oriented and reverse-complement hits would inflate the per-window
/// presence count.
pub fn new_searcher_fwd() -> AmbiguousSearcher {
    Searcher::<Iupac>::new_fwd()
}

/// Encodes equal-length patterns of 1 to `MAX_TILED_PATTERN_LEN` bases for
/// `encoded_pattern_hits`. The encoding holds the bit profiles of every
/// pattern on both strands.
pub fn encode_patterns(patterns: &[Vec<u8>]) -> EncodedAdapterBatch {
    debug_assert!(
        patterns
            .first()
            .is_some_and(|p| (1..=MAX_TILED_PATTERN_LEN).contains(&p.len()))
    );
    // A forward-only searcher encodes the given strand alone; the reverse
    // strand is its own encoding.
    let mut encoder = Searcher::<Iupac>::new_fwd();
    let complements: Vec<Vec<u8>> = patterns.iter().map(|p| Iupac::complement(p)).collect();
    EncodedAdapterBatch {
        forward: encoder.encode_patterns(patterns),
        complement: encoder.encode_patterns(&complements),
    }
}

/// Searches a pre-encoded batch over both strands of `text`, one pattern per
/// SIMD lane, and calls `accept` with each hit's pattern index, text span,
/// cost, and whether the reverse complement matched. `reversed` is `text` reversed, which the caller keeps per read so the
/// reverse strand needs no copy. Hits are the rightmost local minima within
/// `k`, as `hits` returns them. The tiled search uses only the searcher's
/// pattern-tiling state, which its single-pattern searches never touch, so
/// one IUPAC searcher serves both.
pub fn encoded_pattern_hits(
    searcher: &mut AmbiguousSearcher,
    encoded: &EncodedAdapterBatch,
    text: &[u8],
    reversed: &[u8],
    k: usize,
    mut accept: impl FnMut(usize, usize, usize, usize, bool),
) {
    debug_assert_eq!(text.len(), reversed.len());
    for m in searcher.search_encoded_patterns(&encoded.forward, text, k) {
        accept(
            m.pattern_idx,
            m.text_start,
            m.text_end,
            m.cost as usize,
            false,
        );
    }
    let n = text.len();
    for m in searcher.search_encoded_patterns(&encoded.complement, reversed, k) {
        accept(
            m.pattern_idx,
            n - m.text_end,
            n - m.text_start,
            m.cost as usize,
            true,
        );
    }
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
        accept(Hit::from_match(&m, pattern.len()));
    }
}

/// Calls `accept` with the text index and every match of `pattern` in each of
/// `texts` within `k` edits. The texts share one pattern encoding and run in
/// parallel SIMD lanes, so two end windows cost about one search instead of
/// two; the spans are those `for_each_hit` reports on each text alone.
pub fn for_each_hit_in_texts<P: Profile, T: RcSearchAble>(
    searcher: &mut Searcher<P>,
    pattern: &[u8],
    texts: &[T],
    k: usize,
    mut accept: impl FnMut(usize, Hit),
) {
    for m in searcher.search_texts(pattern, texts, k) {
        accept(m.text_idx, Hit::from_match(&m, pattern.len()));
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
        for code in *b"UX.-0" {
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
