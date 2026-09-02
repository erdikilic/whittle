use sassy::profiles::{Dna, Iupac};
use sassy::{CachedRev, Searcher};

/// The fast searcher: sassy's DNA profile, for all-ACGT patterns against
/// all-ACGT text.
///
/// It PANICS during traceback on any other byte, in the pattern or the text, so
/// every use must be gated on `is_plain_acgt` for both. That is why the two
/// profiles are kept side by side instead of standardizing on IUPAC: on a
/// narrowed adapter set, where most patterns are searched one at a time rather
/// than batched across SIMD lanes, IUPAC costs several times more.
pub type PlainSearcher = Searcher<Dna>;

/// The general searcher: sassy's IUPAC profile, which handles ambiguity codes in
/// the pattern (a degenerate primer) and in the text (an `N` in a read).
pub type AmbiguousSearcher = Searcher<Iupac>;

/// The batched, pattern-parallel searcher. Sassy implements `search_patterns`
/// only for IUPAC, so the batched path is always ambiguity-safe.
pub type BatchedAdapterSearcher = Searcher<Iupac>;

/// Whether every byte is a plain uppercase-or-lowercase A/C/G/T, meaning
/// `PlainSearcher` can be used without risking its traceback panic.
pub fn is_plain_acgt(seq: &[u8]) -> bool {
    seq.iter()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

/// The plain bases an IUPAC code stands for, or `None` if the byte is not a
/// nucleotide code at all.
///
/// `U` is absent: callers fold it to `T` before it reaches here, because sassy
/// would otherwise treat `U` as a fifth base that matches nothing in a DNA read.
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

/// How many of the four bases an IUPAC code stands for, or `None` if the byte is
/// not a nucleotide code at all. See `iupac_bases` for the `U` rule.
pub fn iupac_degeneracy(code: u8) -> Option<u8> {
    iupac_bases(code).map(|bases| bases.len() as u8)
}

/// One approximate match of a pattern in the text: half-open `[start, end)` into
/// the text, with its edit `cost`. Strand is not exposed. A reverse-complement
/// hit still occupies the same text span, which is all the trimmer needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    pub start: usize,
    pub end: usize,
    pub cost: usize,
}

/// A fresh searcher configured to also match the reverse-complement strand.
pub fn new_searcher() -> PlainSearcher {
    Searcher::<Dna>::new_rc()
}

/// An ambiguity-tolerant searcher over both strands.
pub fn new_ambiguous_searcher() -> AmbiguousSearcher {
    Searcher::<Iupac>::new_rc()
}

/// A fresh searcher that matches the FORWARD strand only (no reverse-complement).
/// Used by inference's k-mer recount, where each read-end window is already
/// strand-oriented and RC hits would inflate the per-window presence count.
pub fn new_searcher_fwd() -> AmbiguousSearcher {
    Searcher::<Iupac>::new_fwd()
}

pub fn new_batched_searcher() -> BatchedAdapterSearcher {
    Searcher::<Iupac>::new_rc()
}

/// Search equal-length patterns together, packing them across SIMD lanes.
/// This retains `search`'s reverse-text semantics while avoiding the repeated
/// short-text setup of calling `search` once per pattern.
pub fn pattern_hits(
    searcher: &mut BatchedAdapterSearcher,
    patterns: &[Vec<u8>],
    text: &[u8],
    k: usize,
) -> Vec<sassy::Match> {
    // `search_patterns` processes SIMD-sized chunks internally. Without this
    // wrapper Sassy rebuilds the reversed text once per chunk.
    let cached_text = CachedRev::new(text, true);
    searcher.search_patterns(patterns, &cached_text, k)
}

/// All matches of `pattern` in `text` with edit distance <= `k`, as text
/// spans. Which strand(s) are searched depends on how `searcher` was built:
/// `new_searcher` (RC-enabled) matches both strands, `new_searcher_fwd`
/// matches the forward strand only. Reuses `searcher`'s internal buffers
/// across calls.
pub fn hits<P: sassy::profiles::Profile>(
    searcher: &mut Searcher<P>,
    pattern: &[u8],
    text: &[u8],
    k: usize,
) -> Vec<Hit> {
    searcher
        .search(pattern, text, k)
        .into_iter()
        .map(|m| Hit {
            start: m.text_start,
            end: m.text_end,
            // sassy's `Match::cost` is `pa_types::Cost` (`i32`, signed to support
            // other algorithms in that crate); an actual match is always within
            // the non-negative `k` budget, so the cast is lossless here.
            cost: m.cost as usize,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // NOTE: `new_rc()` returns TWO same-span hits for a reverse-complement-
    // palindromic pattern (Fwd + Rc). Count-based tests must use a NON-palindromic
    // pattern (revcomp absent from the text) to get exactly one hit. `adapter_segments`
    // is unaffected: it dedups terminal hits via max/min and merges interior ones.
    #[test]
    fn exact_forward_match() {
        let mut s = new_searcher();
        // revcomp(AAAACCCCGGGG) = CCCCGGGGTTTT, absent from the text -> one hit.
        let h = hits(&mut s, b"AAAACCCCGGGG", b"TTAAAACCCCGGGGTT", 0);
        assert_eq!(h.len(), 1);
        assert_eq!((h[0].start, h[0].end, h[0].cost), (2, 14, 0));
    }

    #[test]
    fn finds_reverse_complement() {
        // pattern AAAACCCC ; revcomp = GGGGTTTT ; embed GGGGTTTT in the text.
        let mut s = new_searcher();
        let h = hits(&mut s, b"AAAACCCC", b"TTGGGGTTTTAA", 0);
        assert_eq!(h.len(), 1);
        assert_eq!((h[0].start, h[0].end), (2, 10));
    }

    #[test]
    fn forward_searcher_ignores_reverse_complement() {
        let mut s = new_searcher_fwd();
        // revcomp(AAAACCCC) = GGGGTTTT is in the text, but forward-only must skip it.
        assert_eq!(hits(&mut s, b"AAAACCCC", b"TTGGGGTTTTAA", 0).len(), 0);
        // the forward pattern present -> found.
        assert_eq!(hits(&mut s, b"AAAACCCC", b"TTAAAACCCCTT", 0).len(), 1);
    }

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

    #[test]
    fn tolerates_one_mismatch_within_budget() {
        let mut s = new_searcher();
        // one substitution (pos 5, C->A) in AAAACCCCGGGG; revcomp absent from text.
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
