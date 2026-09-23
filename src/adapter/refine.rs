//! Placement of a trim boundary inside an accepted hit.
//!
//! An edit-distance hit spans the whole pattern, so a pattern that runs past
//! the technical sequence into the insert (a catalog entry that shares its
//! core with a shorter one and continues beyond it) reports an end inside the
//! insert, the extra bases paid for as edits. The inner boundary of a hit is
//! therefore placed at the end of the best-scoring alignment of a pattern
//! prefix, scored `+1` per match and `-PENALTY` per mismatch or gap: the
//! aligned columns after that point score below zero and are left to the
//! insert.

use std::cell::RefCell;

/// Score of a mismatch or gap column; a match scores 1.
const PENALTY: i32 = 2;

thread_local! {
    /// Alignment rows of this thread, reused across hits.
    static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default());
}

/// Scratch rows for `aligned_prefix_len`, reused across calls.
#[derive(Debug, Default)]
struct Scratch {
    prev: Vec<i32>,
    cur: Vec<i32>,
    best: Vec<i32>,
    pattern: Vec<u8>,
    text: Vec<u8>,
}

/// The bases an IUPAC code stands for as a mask, A 1, C 2, G 4 and T 8; zero
/// for any other byte, including the `X` of an uncalled read base.
const BASE_MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    let codes: [(u8, u8); 15] = [
        (b'A', 1),
        (b'C', 2),
        (b'G', 4),
        (b'T', 8),
        (b'R', 5),
        (b'Y', 10),
        (b'S', 6),
        (b'W', 9),
        (b'K', 12),
        (b'M', 3),
        (b'B', 14),
        (b'D', 13),
        (b'H', 11),
        (b'V', 7),
        (b'N', 15),
    ];
    let mut i = 0;
    while i < codes.len() {
        t[codes[i].0 as usize] = codes[i].1;
        t[codes[i].0.to_ascii_lowercase() as usize] = codes[i].1;
        i += 1;
    }
    t
};

/// Returns whether read base `t` is one of the bases IUPAC code `p` stands
/// for. Read bases are uppercase A/C/G/T or `X`, which matches nothing.
fn matches(p: u8, t: u8) -> bool {
    BASE_MASK[usize::from(p)] & BASE_MASK[usize::from(t)] != 0
}

/// Returns the length of the text prefix that ends the best-scoring
/// alignment of `text` against a substring of `pattern` starting at the text
/// start, the later end on a tie. Pattern bases before the substring are
/// free, so a hit whose outer part hangs off the read is placed the same way.
fn aligned_prefix_len(pattern: &[u8], text: &[u8], s: &mut Scratch) -> usize {
    let n = text.len();
    let Scratch {
        prev, cur, best, ..
    } = s;
    prev.clear();
    prev.extend((0..=n).map(|j| -(j as i32) * PENALTY));
    best.clear();
    best.extend_from_slice(prev);
    for &p in pattern {
        cur.clear();
        cur.push(0);
        for j in 1..=n {
            let diagonal = prev[j - 1] + if matches(p, text[j - 1]) { 1 } else { -PENALTY };
            let score = diagonal.max(prev[j] - PENALTY).max(cur[j - 1] - PENALTY);
            cur.push(score);
            best[j] = best[j].max(score);
        }
        std::mem::swap(prev, cur);
    }
    let top = best.iter().copied().max().unwrap_or(0);
    best.iter().rposition(|&v| v == top).unwrap_or(0)
}

/// Returns the refined inner end of a hit at `text[start..end)` of `pattern`,
/// or of its reverse complement when `reverse`: the hit keeps its start and
/// ends where the best-scoring alignment ends. The result lies in
/// `start..=end`.
pub(super) fn refined_end(
    pattern: &[u8],
    reverse: bool,
    text: &[u8],
    start: usize,
    end: usize,
) -> usize {
    SCRATCH.with_borrow_mut(|s| {
        let mut oriented = std::mem::take(&mut s.pattern);
        oriented.clear();
        if reverse {
            oriented.extend(pattern.iter().rev().map(|&b| super::complement(b)));
        } else {
            oriented.extend_from_slice(pattern);
        }
        let len = aligned_prefix_len(&oriented, &text[start..end], s);
        s.pattern = oriented;
        start + len
    })
}

/// Returns the refined inner start of a hit at `text[start..end)`: the mirror
/// of `refined_end`, with the alignment anchored at the hit end.
pub(super) fn refined_start(
    pattern: &[u8],
    reverse: bool,
    text: &[u8],
    start: usize,
    end: usize,
) -> usize {
    SCRATCH.with_borrow_mut(|s| {
        let mut span = std::mem::take(&mut s.text);
        span.clear();
        span.extend(text[start..end].iter().rev());
        let mut oriented = std::mem::take(&mut s.pattern);
        oriented.clear();
        if reverse {
            oriented.extend(pattern.iter().map(|&b| super::complement(b)));
        } else {
            oriented.extend(pattern.iter().rev());
        }
        let len = aligned_prefix_len(&oriented, &span, s);
        s.pattern = oriented;
        s.text = span;
        end - len
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_hit_keeps_its_span() {
        let text = b"GGGGACGTACGTACGTTTTT";
        assert_eq!(refined_end(b"ACGTACGTACGT", false, text, 4, 16), 16);
        assert_eq!(refined_start(b"ACGTACGTACGT", false, text, 4, 16), 4);
    }

    #[test]
    fn mismatched_tail_is_left_to_the_insert() {
        // The pattern continues with GGGG where the read has TTTT.
        let text = b"ACGTTGCAACGTTTTTCCAG";
        assert_eq!(refined_end(b"ACGTTGCAACGTGGGG", false, text, 0, 16), 12);
    }

    #[test]
    fn mismatched_head_is_left_to_the_insert() {
        let text = b"CCAGTTTTACGTTGCAACGT";
        assert_eq!(refined_start(b"GGGGACGTTGCAACGT", false, text, 4, 20), 8);
    }

    #[test]
    fn reverse_complement_hit_is_refined_on_its_own_strand() {
        let pattern = b"ACGTTGCAACGTGGGG";
        let rc: Vec<u8> = pattern
            .iter()
            .rev()
            .map(|&b| crate::adapter::complement(b))
            .collect();
        let mut text = b"CATG".to_vec();
        text.extend_from_slice(b"AAAA");
        text.extend_from_slice(&rc[4..]);
        let n = text.len();
        assert_eq!(refined_start(pattern, true, &text, 4, n), 8);
    }

    #[test]
    fn degenerate_codes_match_their_bases() {
        assert_eq!(
            refined_end(b"ACGTNNACGTAC", false, b"ACGTGAACGTAC", 0, 12),
            12
        );
    }

    #[test]
    fn uncalled_bases_match_nothing() {
        assert_eq!(
            refined_end(b"ACGTACGTACGT", false, b"ACGTACGTAXXX", 0, 12),
            9
        );
    }
}
