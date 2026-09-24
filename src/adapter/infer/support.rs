//! Support of a candidate among read-end windows, sequence identity between
//! candidates, and catalog naming of discovered sequences.

use super::*;

/// Requires supporting alignments to concentrate near the physical read end.
pub(super) fn terminal_support(
    seq: &[u8],
    windows: &[&[u8]],
    end: End,
    edits: usize,
) -> (usize, usize) {
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let mut best = vec![None; windows.len()];
    for hit in searcher.search_texts(seq, windows, edits) {
        let distance = match end {
            End::Five => hit.text_start,
            End::Three => windows[hit.text_idx].len() - hit.text_end,
        };
        let key = (hit.cost, distance);
        let entry = &mut best[hit.text_idx];
        if entry.is_none_or(|old| key < old) {
            *entry = Some(key);
        }
    }
    let present = best.iter().filter(|hit| hit.is_some()).count();
    let anchored = best
        .iter()
        .flatten()
        .filter(|(_, distance)| *distance <= ANCHOR_SLACK)
        .count();
    (present, anchored)
}

/// Counts the distinct `windows` with at least one forward approximate
/// occurrence of `pattern` within `max_edits`. Each window counts at most once,
/// however often `pattern` occurs in it. `searcher` must be forward-only (see
/// `new_searcher_fwd`) so reverse-complement occurrences do not inflate the
/// count. Callers provide an already bounded window sample.
pub(super) fn windows_containing(
    searcher: &mut AmbiguousSearcher,
    pattern: &[u8],
    windows: &[&[u8]],
    max_edits: usize,
) -> u32 {
    let mut seen = vec![false; windows.len()];
    crate::adapter::search::for_each_hit_in_texts(
        searcher,
        pattern,
        windows,
        max_edits,
        |index, _| seen[index] = true,
    );
    seen.into_iter().filter(|&present| present).count() as u32
}

/// Returns, per window, whether `pattern` occurs within `max_edits`.
pub(super) fn windows_with(
    searcher: &mut AmbiguousSearcher,
    pattern: &[u8],
    windows: &[&[u8]],
    max_edits: usize,
) -> Vec<bool> {
    let mut seen = vec![false; windows.len()];
    crate::adapter::search::for_each_hit_in_texts(
        searcher,
        pattern,
        windows,
        max_edits,
        |index, _| seen[index] = true,
    );
    seen
}

/// Returns whether `a` and `b` are the same adapter within `error_rate`: an
/// approximate occurrence of the shorter in the longer on either strand (the
/// both-strand searcher covers the reverse-complement case).
pub(super) fn same_adapter(a: &[u8], b: &[u8], error_rate: f64) -> bool {
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if short.len() < MIN_PATTERN_LEN {
        return short == long;
    }
    let k = edit_budget(error_rate, short.len());
    let mut s = new_ambiguous_searcher();
    !hits(&mut s, short, long, k).is_empty()
}

/// Returns the best catalog matches for `seq` as `(name, percent_identity)`,
/// sorted by identity descending, at most three, and only at or above
/// `NAME_IDENTITY_MIN`. The result annotates an inferred adapter with the
/// catalog entry it corresponds to.
pub(super) fn name_against(seq: &[u8], refs: &[Adapter], error_rate: f64) -> Vec<(String, f32)> {
    let mut s = new_ambiguous_searcher();
    let mut named: Vec<(String, f32)> = Vec::new();
    // A barcode construct's `N` block matches any sequence and names nothing.
    for r in refs.iter().filter(|r| !crate::adapter::is_construct(r)) {
        let (short, long) = if seq.len() <= r.seq.len() {
            (seq, r.seq.as_slice())
        } else {
            (r.seq.as_slice(), seq)
        };
        if short.len() < MIN_PATTERN_LEN {
            continue;
        }
        let k = edit_budget(error_rate, short.len());
        if let Some(h) = hits(&mut s, short, long, k)
            .into_iter()
            .min_by_key(|h| h.cost)
        {
            let pct = 100.0 * (1.0 - h.cost as f32 / short.len() as f32);
            if pct >= NAME_IDENTITY_MIN {
                named.push((r.name.clone(), pct));
            }
        }
    }
    named.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap().then(a.0.cmp(&b.0)));
    named.truncate(3);
    named
}

/// Samples at most `cap` windows deterministically, spread across the whole
/// slice by stride rather than taken from its start, so the recount and
/// support frame in `assemble` is not order-biased. Returns every window in
/// order when `windows.len() <= cap`.
pub(super) fn stride_sample<'a>(windows: &[&'a [u8]], cap: usize) -> Vec<&'a [u8]> {
    let step = windows.len().div_ceil(cap.max(1)).max(1);
    windows.iter().step_by(step).copied().collect()
}

/// Returns whether two candidates reconstruct one family: the same adapter
/// within the error budget, or reconstructions whose longest common
/// substring covers `FAMILY_OVERLAP_PERCENT` of the shorter, as remnants of
/// one adapter that differ in length do.
pub(super) fn same_family(a: &[u8], b: &[u8], error_rate: f64) -> bool {
    same_adapter(a, b, error_rate)
        || shared_substring(a, b) * 100 >= a.len().min(b.len()) * FAMILY_OVERLAP_PERCENT
}

/// Returns whether the shorter of two candidates is a fragment of the
/// longer: one family by `same_family`, with less than `MIN_PATTERN_LEN` of
/// the shorter outside their longest common substring, as an end of the
/// longer extended by a few bases the assembly could not support is.
pub(super) fn fragment_of(a: &[u8], b: &[u8], error_rate: f64) -> bool {
    let short = a.len().min(b.len());
    same_adapter(a, b, error_rate) || {
        let shared = shared_substring(a, b);
        shared * 100 >= short * FAMILY_OVERLAP_PERCENT && short - shared < MIN_PATTERN_LEN
    }
}

/// Length of the longest substring `a` shares with `b` on either strand.
pub(super) fn shared_substring(a: &[u8], b: &[u8]) -> usize {
    longest_common_substring(a, b).max(longest_common_substring(
        &crate::adapter::reverse_complement(a),
        b,
    ))
}

/// Length of the longest common substring of `a` and `b`.
pub(super) fn longest_common_substring(a: &[u8], b: &[u8]) -> usize {
    let mut longest = 0;
    let mut previous = vec![0usize; b.len() + 1];
    for &x in a {
        let mut current = vec![0usize; b.len() + 1];
        for (j, &y) in b.iter().enumerate() {
            if x == y {
                current[j + 1] = previous[j] + 1;
                longest = longest.max(current[j + 1]);
            }
        }
        previous = current;
    }
    longest
}
