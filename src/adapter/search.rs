use sassy::profiles::Dna;
use sassy::Searcher;

/// Reusable DNA searcher (searches a pattern against both strands of the text).
pub type DnaSearcher = Searcher<Dna>;

/// One approximate match of a pattern in the text: half-open `[start, end)` into
/// the text, with its edit `cost`. Strand is not exposed — a reverse-complement
/// hit still occupies the same text span, which is all the trimmer needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit {
    pub start: usize,
    pub end: usize,
    pub cost: usize,
}

/// A fresh searcher configured to also match the reverse-complement strand.
pub fn new_searcher() -> DnaSearcher {
    Searcher::<Dna>::new_rc()
}

/// All matches of `pattern` in `text` with edit distance <= `k` (both strands),
/// as text spans. Reuses `searcher`'s internal buffers across calls.
pub fn hits(searcher: &mut DnaSearcher, pattern: &[u8], text: &[u8], k: usize) -> Vec<Hit> {
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
    // is unaffected — it dedups terminal hits via max/min and merges interior ones.
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
    fn tolerates_one_mismatch_within_budget() {
        let mut s = new_searcher();
        // one substitution (pos 5, C->A) in AAAACCCCGGGG; revcomp absent from text.
        assert_eq!(hits(&mut s, b"AAAACCCCGGGG", b"TTAAAACACCGGGGTT", 1).len(), 1);
        assert_eq!(hits(&mut s, b"AAAACCCCGGGG", b"TTAAAACACCGGGGTT", 0).len(), 0);
    }
}
