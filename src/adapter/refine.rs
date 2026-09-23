//! Placement of a trim boundary inside an accepted hit.
//!
//! An edit-distance hit spans the whole pattern, so a pattern that runs past
//! the technical sequence into the insert (a catalog entry that shares its
//! core with a shorter one and continues beyond it) reports an end inside the
//! insert, the extra bases paid for as edits. Each end of a hit is therefore
//! moved inward past the stretch of its alignment that scores below zero at
//! `+1` per match and `-PENALTY` per mismatch or gap: those columns hold more
//! errors than matches and are left to the insert.

/// Score of a mismatch or gap column; a match scores 1.
const PENALTY: isize = 2;

/// Returns the text bases to drop from the end the columns of `ops` start
/// at: the text consumed by the lowest-scoring run of leading columns, or
/// zero when every run scores at least zero. Operations are sassy CIGAR
/// characters with their counts; the running score falls within a mismatch
/// or gap run and rises within a match run, so its minimum lies at a run
/// boundary.
fn clip_from(ops: impl Iterator<Item = (char, usize)>) -> usize {
    let (mut score, mut consumed) = (0isize, 0usize);
    let (mut lowest, mut clip) = (0isize, 0usize);
    for (op, count) in ops {
        let weight = if op == '=' { 1 } else { -PENALTY };
        score += weight * count as isize;
        if op != 'I' {
            consumed += count;
        }
        if score < lowest {
            lowest = score;
            clip = consumed;
        }
    }
    clip
}

/// Returns the text bases to drop at the text start and the text end of
/// match `m`. Sassy's CIGAR runs in pattern order, which is text order for a
/// forward match and the reverse of it for a reverse-complement match (`rc`).
pub(super) fn clips(m: &sassy::Match, rc: bool) -> (usize, usize) {
    let ops = || {
        m.cigar
            .ops
            .iter()
            .map(|e| (e.op.to_char(), usize::try_from(e.cnt).unwrap_or(0)))
    };
    let (head, tail) = (clip_from(ops()), clip_from(ops().rev()));
    if rc { (tail, head) } else { (head, tail) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::reverse_complement;
    use crate::adapter::search::{Hit, hits, new_searcher};

    fn only_hit(pattern: &[u8], text: &[u8], k: usize) -> Hit {
        let found = hits(&mut new_searcher(), pattern, text, k);
        assert_eq!(found.len(), 1, "{found:?}");
        found[0]
    }

    #[test]
    fn exact_hit_keeps_its_span() {
        let hit = only_hit(b"ACGTTGCAACGTAC", b"GGGGGACGTTGCAACGTACGGGGG", 0);
        assert_eq!((hit.clip_start, hit.clip_end), (0, 0));
    }

    #[test]
    fn mismatched_tail_is_left_to_the_insert() {
        // The pattern continues with GGGG where the read has TTTT.
        let hit = only_hit(b"CAGTTGCAACGTGGGG", b"AAAAACAGTTGCAACGTTTTTCCAAC", 4);
        assert_eq!(hit.end - hit.clip_end, 17);
    }

    #[test]
    fn mismatched_head_is_left_to_the_insert() {
        let hit = only_hit(b"GGGGCAGTTGCAACGT", b"AACCTTTTCAGTTGCAACGTAACC", 4);
        assert_eq!(hit.start + hit.clip_start, 8);
    }

    #[test]
    fn reverse_complement_hit_clips_its_own_ends() {
        let pattern = b"CAGTTGCAACGTGGGG";
        let text = reverse_complement(b"AAAAACAGTTGCAACGTTTTTCCAAC");
        let hit = only_hit(pattern, &text, 4);
        assert_eq!(hit.start + hit.clip_start, text.len() - 17);
        assert_eq!(hit.clip_end, 0);
    }

    #[test]
    fn clip_ends_at_the_lowest_scoring_run() {
        assert_eq!(clip_from([('X', 1), ('=', 5)].into_iter()), 1);
        assert_eq!(clip_from([('=', 1), ('X', 1), ('=', 5)].into_iter()), 2);
        assert_eq!(clip_from([('X', 1), ('=', 2), ('X', 1)].into_iter()), 1);
        assert_eq!(clip_from([('X', 1), ('=', 1), ('X', 1)].into_iter()), 3);
        assert_eq!(clip_from([('=', 2), ('X', 1), ('=', 5)].into_iter()), 0);
        assert_eq!(clip_from([('I', 2), ('=', 5)].into_iter()), 0);
        assert_eq!(clip_from([('D', 2), ('=', 5)].into_iter()), 2);
    }
}
