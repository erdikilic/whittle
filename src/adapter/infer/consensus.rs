//! Consensus refinement from aligned windows: polishing, upstream extension,
//! and the insert-facing end placed by base conservation against the read
//! composition.

use super::*;

/// Encodes a set of concrete DNA bases as an IUPAC symbol.
pub(super) fn ambiguity_code(mask: usize) -> u8 {
    b"-ACMGRSVTWYHKDBN"[mask]
}

/// Returns, for each base of `seq`, how often the best alignment of each
/// window within `edits` places A, C, G, T or a deletion against it, or
/// `None` when fewer than `MIN_SUPPORT_WINDOWS` windows align.
pub(super) fn column_counts(
    seq: &[u8],
    windows: &[&[u8]],
    edits: usize,
) -> Option<Vec<[usize; 5]>> {
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let mut best: Vec<Option<sassy::Match>> = vec![None; windows.len()];
    for hit in searcher.search_texts(seq, windows, edits) {
        let entry = &mut best[hit.text_idx];
        if entry
            .as_ref()
            .is_none_or(|old| (hit.cost, hit.text_start) < (old.cost, old.text_start))
        {
            *entry = Some(hit.clone());
        }
    }
    let mut counts = vec![[0usize; 5]; seq.len()];
    let mut aligned = 0;
    for (window, hit) in windows.iter().zip(best) {
        let Some(hit) = hit else { continue };
        aligned += 1;
        let path = hit.to_path();
        for (i, pos) in path.iter().enumerate() {
            let next = path
                .get(i + 1)
                .map(|p| (p.0, p.1))
                .unwrap_or((seq.len() as i32, hit.text_end as i32));
            if next.0 != pos.0 + 1 {
                continue;
            }
            let base = if next.1 == pos.1 + 1 {
                encode_kmer(&window[pos.1 as usize..pos.1 as usize + 1]).unwrap() as usize
            } else {
                4
            };
            counts[pos.0 as usize][base] += 1;
        }
    }
    (aligned >= MIN_SUPPORT_WINDOWS).then_some(counts)
}

/// Refines an assembled path with one best alignment per supporting window.
/// Majority substitutions and deletions correct graph branches caused by
/// sequencing errors without extending the assembly into unaligned sequence.
pub(super) fn polish_consensus(seq: &[u8], windows: &[&[u8]]) -> Vec<u8> {
    let Some(counts) = column_counts(seq, windows, edit_budget(0.25, seq.len())) else {
        return seq.to_vec();
    };
    seq.iter()
        .zip(counts)
        .filter_map(|(&base, counts)| {
            let total: usize = counts.iter().sum();
            if counts[4] * 2 > total {
                return None;
            }
            let best = (0..4).max_by_key(|&i| counts[i]).unwrap();
            Some(if counts[best] * 100 >= total * 70 {
                b"ACGT"[best]
            } else {
                base
            })
        })
        .collect()
}

/// Returns the fraction of A, C, G and T among the bases of `windows`, or a
/// uniform composition when they hold none.
pub(super) fn base_composition(windows: &[&[u8]]) -> [f64; 4] {
    let mut counts = [0usize; 4];
    for window in windows {
        for base in window.iter() {
            if let Some(code) = encode_kmer(std::slice::from_ref(base)) {
                counts[code as usize] += 1;
            }
        }
    }
    let total: usize = counts.iter().sum();
    if total == 0 {
        return [0.25; 4];
    }
    counts.map(|count| count as f64 / total as f64)
}

/// Returns whether a tally of read bases at one position is conserved: one
/// base, or a pair of bases, holds at least the midpoint between its share of
/// the `composition` and one. A technical base is read at the per-base
/// accuracy of the platform and clears the midpoint, and so does a two-fold
/// degenerate primer base on its pair; an insert position holds any base or
/// pair at its composition share, which the midpoint exceeds for any
/// composition.
pub(super) fn conserved(counts: [usize; 4], composition: [f64; 4]) -> bool {
    let total = counts.iter().sum::<usize>() as f64;
    if total == 0.0 {
        return false;
    }
    let clears = |count: usize, share: f64| 2.0 * count as f64 >= (1.0 + share) * total;
    (0..4).any(|i| {
        clears(counts[i], composition[i])
            || (i + 1..4).any(|j| clears(counts[i] + counts[j], composition[i] + composition[j]))
    })
}

/// The occurrences of every `MIN_PATTERN_LEN`-mer in the windows of one read
/// end, each with the base that follows it inward, for the conservation
/// tallies of `trim_unconserved_inner_end` and of the insert boundary. An
/// exact match does not collect windows that carry a similar sequence, such as
/// another member of a barcode panel, and the base after it lies outside the
/// match, so gap placement cannot bias it toward the consensus.
pub(super) struct FollowingBases {
    /// Every occurrence as the packed k-mer, the window, the distance of the
    /// k-mer's outer edge from the physical end, and the code of the
    /// following base, sorted.
    pub(super) occurrences: Vec<(u32, u32, u16, u8)>,
}

impl FollowingBases {
    /// Indexes the k-mers of `windows` whose outer edge lies within
    /// `ANCHOR_SLACK + LMAX` of the physical `end`.
    pub(super) fn new(windows: &[&[u8]], end: End) -> Self {
        let mut occurrences = Vec::new();
        for (idx, window) in windows.iter().enumerate() {
            let n = window.len();
            for distance in 0..n
                .saturating_sub(MIN_PATTERN_LEN)
                .min(ANCHOR_SLACK + LMAX + 1)
            {
                let (kmer, next) = match end {
                    End::Five => (
                        &window[distance..distance + MIN_PATTERN_LEN],
                        window[distance + MIN_PATTERN_LEN],
                    ),
                    End::Three => {
                        let stop = n - distance;
                        (
                            &window[stop - MIN_PATTERN_LEN..stop],
                            window[stop - MIN_PATTERN_LEN - 1],
                        )
                    },
                };
                if let (Some(code), Some(next)) =
                    (encode_kmer(kmer), encode_kmer(std::slice::from_ref(&next)))
                {
                    occurrences.push((code as u32, idx as u32, distance as u16, next as u8));
                }
            }
        }
        occurrences.sort_unstable();
        Self { occurrences }
    }

    /// Returns the counts of A, C, G and T in the base that follows `inner`,
    /// a `MIN_PATTERN_LEN`-base sequence with `outboard` further bases between
    /// it and the physical end, or `None` when fewer than
    /// `MIN_SUPPORT_WINDOWS` windows hold it. Each window contributes its
    /// occurrence nearest the end, when the sequence as a whole then starts
    /// within `ANCHOR_SLACK` of it.
    pub(super) fn counts(&self, inner: &[u8], outboard: usize) -> Option<[usize; 4]> {
        let code = encode_kmer(inner)? as u32;
        let first = self.occurrences.partition_point(|o| o.0 < code);
        let last = self.occurrences.partition_point(|o| o.0 <= code);
        let mut counts = [0usize; 4];
        let mut previous = None;
        // Sorted by window, then distance: the first occurrence of each
        // window is its nearest.
        for &(_, window, distance, next) in &self.occurrences[first..last] {
            if previous == Some(window) {
                continue;
            }
            previous = Some(window);
            if (distance as usize).saturating_sub(outboard) <= ANCHOR_SLACK {
                counts[next as usize] += 1;
            }
        }
        (counts.iter().sum::<usize>() >= MIN_SUPPORT_WINDOWS).then_some(counts)
    }

    /// Returns the counts of the base that follows the first `keep` bases of
    /// `seq` counted from the physical end; see `counts`.
    pub(super) fn after(&self, seq: &[u8], keep: usize, end: End) -> Option<[usize; 4]> {
        let outboard = keep.checked_sub(MIN_PATTERN_LEN)?;
        let inner = match end {
            End::Five => &seq[outboard..keep],
            End::Three => &seq[seq.len() - keep..seq.len() - outboard],
        };
        self.counts(inner, outboard)
    }
}

/// Removes the insert-facing bases of `seq` that the supporting windows do
/// not conserve. Cut points are tried from the outer edge of the last k-mer
/// inward; the first whose following base is not `conserved` against
/// `composition` ends the sequence. Cut points without a tally (an IUPAC code
/// or too few windows) are skipped. K-mer support cannot place this end for a
/// layer about one k-mer long, because erosion at the physical end weakens the
/// k-mer spanning the whole layer.
pub(super) fn trim_unconserved_inner_end(
    seq: &[u8],
    following: &FollowingBases,
    end: End,
    composition: [f64; 4],
) -> Vec<u8> {
    let len = seq.len();
    for keep in MIN_PATTERN_LEN.max(len.saturating_sub(KMER_K))..len {
        match following.after(seq, keep, end) {
            Some(counts) if !conserved(counts, composition) => {
                return match end {
                    End::Five => seq[..keep].to_vec(),
                    End::Three => seq[len - keep..].to_vec(),
                };
            },
            _ => {},
        }
    }
    seq.to_vec()
}

/// Extends a conserved insert anchor toward the physical read end. Each round
/// aligns the current consensus before voting on the preceding base, allowing
/// indels to shift the supporting reads without shifting the consensus.
pub(super) fn upstream_consensus(anchor: &[u8], windows: &[&[u8]]) -> Vec<u8> {
    let mut seq = anchor.to_vec();
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    for _ in 0..LMAX - anchor.len() {
        let mut best = vec![None; windows.len()];
        for hit in searcher.search_texts(&seq, windows, edit_budget(0.12, seq.len())) {
            let entry = &mut best[hit.text_idx];
            let key = (hit.cost, hit.text_start);
            if entry.is_none_or(|old| key < old) {
                *entry = Some(key);
            }
        }
        let mut counts = [0usize; 4];
        for (window, hit) in windows.iter().zip(best) {
            if let Some((_, start)) = hit.filter(|&(_, start)| start > 0)
                && let Some(base) = encode_kmer(&window[start - 1..start])
            {
                counts[base as usize] += 1;
            }
        }
        let total: usize = counts.iter().sum();
        if total < MIN_SUPPORT_WINDOWS.max((windows.len() as f64 * KEEP_SUPPORT).ceil() as usize) {
            break;
        }
        let mut ranked = [0, 1, 2, 3];
        ranked.sort_by_key(|&i| (std::cmp::Reverse(counts[i]), i));
        let first = ranked[0];
        let mask = if counts[first] * 100 >= total * 70 {
            1 << first
        } else if (counts[first] + counts[ranked[1]]) * 100 >= total * 80
            && counts[ranked[1]] * 100 >= total * 20
        {
            (1 << first) | (1 << ranked[1])
        } else {
            break;
        };
        seq.insert(0, ambiguity_code(mask));
        if seq.len() >= 2 * anchor.len() && is_repetitive(&seq) {
            return Vec::new();
        }
    }
    seq.truncate(seq.len() - anchor.len());
    seq
}
