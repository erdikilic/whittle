use sassy::CachedRev;
use sassy::Searcher;
use sassy::profiles::Iupac;

/// Reusable adapter searcher (searches a pattern against both strands of the text).
///
/// The IUPAC profile, not `Dna`, and this is load-bearing: sassy's `Dna` profile
/// PANICS during traceback if the text holds any byte outside `ACGTacgt`, and
/// real ONT reads contain `N`. IUPAC treats an ambiguity code as matching the
/// bases it stands for, which is both the standard reading and the only one that
/// does not abort the run. On pure A/C/G/T text the two profiles are equivalent,
/// so this costs nothing on the common case. Sassy's pattern-batched API is
/// IUPAC-only anyway, so this makes one profile serve both paths.
pub type AdapterSearcher = Searcher<Iupac>;

/// The same profile, used through the pattern-batched API.
pub type BatchedAdapterSearcher = Searcher<Iupac>;

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
pub fn new_searcher() -> AdapterSearcher {
    Searcher::<Iupac>::new_rc()
}

/// A fresh searcher that matches the FORWARD strand only (no reverse-complement).
/// Used by inference's k-mer recount, where each read-end window is already
/// strand-oriented and RC hits would inflate the per-window presence count.
pub fn new_searcher_fwd() -> AdapterSearcher {
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
pub fn hits(searcher: &mut AdapterSearcher, pattern: &[u8], text: &[u8], k: usize) -> Vec<Hit> {
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
