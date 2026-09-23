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

/// Returns whether read base `t` is one of the bases IUPAC code `p` stands
/// for. Read bases are uppercase A/C/G/T or `X` for an uncalled base, which
/// matches nothing.
fn matches(p: u8, t: u8) -> bool {
    p == t || super::iupac_bases(p).is_some_and(|bases| bases.len() > 1 && bases.contains(&t))
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

/// Returns the refined inner end of a hit at `text[start..end)` whose pattern
/// is `pattern`: the hit keeps its start and ends where the best-scoring
/// alignment of either strand of the pattern ends. The result lies in
/// `start..=end`.
pub(super) fn refined_end(pattern: &[u8], text: &[u8], start: usize, end: usize) -> usize {
    SCRATCH.with_borrow_mut(|s| end_with(pattern, text, start, end, s))
}

/// Returns the refined inner start of a hit at `text[start..end)`: the mirror
/// of `refined_end`, with the alignment anchored at the hit end.
pub(super) fn refined_start(pattern: &[u8], text: &[u8], start: usize, end: usize) -> usize {
    SCRATCH.with_borrow_mut(|s| start_with(pattern, text, start, end, s))
}

fn end_with(pattern: &[u8], text: &[u8], start: usize, end: usize, s: &mut Scratch) -> usize {
    let span = &text[start..end];
    let forward = aligned_prefix_len(pattern, span, s);
    let mut rc = std::mem::take(&mut s.pattern);
    rc.clear();
    rc.extend(pattern.iter().rev().map(|&b| super::complement(b)));
    let reverse = aligned_prefix_len(&rc, span, s);
    s.pattern = rc;
    start + forward.max(reverse)
}

fn start_with(pattern: &[u8], text: &[u8], start: usize, end: usize, s: &mut Scratch) -> usize {
    let mut span = std::mem::take(&mut s.text);
    span.clear();
    span.extend(text[start..end].iter().rev());
    let mut reversed = std::mem::take(&mut s.pattern);
    reversed.clear();
    reversed.extend(pattern.iter().rev());
    let forward = aligned_prefix_len(&reversed, &span, s);
    reversed.clear();
    reversed.extend(pattern.iter().map(|&b| super::complement(b)));
    let reverse = aligned_prefix_len(&reversed, &span, s);
    s.pattern = reversed;
    s.text = span;
    end - forward.max(reverse)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_hit_keeps_its_span() {
        let text = b"GGGGACGTACGTACGTTTTT";
        assert_eq!(refined_end(b"ACGTACGTACGT", text, 4, 16), 16);
        assert_eq!(refined_start(b"ACGTACGTACGT", text, 4, 16), 4);
    }

    #[test]
    fn mismatched_tail_is_left_to_the_insert() {
        // The pattern continues with GGGG where the read has TTTT.
        let text = b"ACGTTGCAACGTTTTTCCAG";
        assert_eq!(refined_end(b"ACGTTGCAACGTGGGG", text, 0, 16), 12);
    }

    #[test]
    fn mismatched_head_is_left_to_the_insert() {
        let text = b"CCAGTTTTACGTTGCAACGT";
        assert_eq!(refined_start(b"GGGGACGTTGCAACGT", text, 4, 20), 8);
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
        assert_eq!(refined_start(pattern, &text, 4, n), 8);
    }

    #[test]
    fn degenerate_codes_match_their_bases() {
        assert_eq!(refined_end(b"ACGTNNACGTAC", b"ACGTGAACGTAC", 0, 12), 12);
    }

    #[test]
    fn uncalled_bases_match_nothing() {
        assert_eq!(refined_end(b"ACGTACGTACGT", b"ACGTACGTAXXX", 0, 12), 9);
    }
}
