//! The layers of discovery: per-read boundaries, the windows of a layer,
//! boundary advance past accepted sequences, the symmetry cut, shared ends
//! between candidates, and variable layers such as barcode panels.

use super::*;

/// Per-read explained depth from each physical end, in bases.
pub(super) struct Boundaries {
    pub(super) five: Vec<usize>,
    pub(super) three: Vec<usize>,
}

impl Boundaries {
    pub(super) fn new(sample: &[&[u8]]) -> Self {
        Self {
            five: vec![0; sample.len()],
            three: sample.iter().map(|r| r.len()).collect(),
        }
    }
}

/// Windows of one read end, each paired with its sample index.
pub(super) type EndWindows<'a> = Vec<(usize, &'a [u8])>;

/// Returns the unexplained windows of at most `w` bases inward from each
/// boundary, paired with their sample indices. Reads without unexplained
/// bases contribute no window.
pub(super) fn layer_windows<'a>(
    sample: &[&'a [u8]],
    bounds: &Boundaries,
    w: usize,
) -> (EndWindows<'a>, EndWindows<'a>) {
    let mut five = Vec::new();
    let mut three = Vec::new();
    for (i, &read) in sample.iter().enumerate() {
        let (b5, b3) = (bounds.five[i], bounds.three[i]);
        if b3 <= b5 {
            continue;
        }
        let head = &read[b5..(b5 + w).min(b3)];
        let tail = &read[b3.saturating_sub(w).max(b5)..b3];
        if is_plain_acgt(head) {
            five.push((i, head));
        }
        if is_plain_acgt(tail) {
            three.push((i, tail));
        }
    }
    (five, three)
}

/// Returns the interior windows `w` bases past each boundary, for reads with
/// at least `4 * w` unexplained bases.
pub(super) fn background_windows<'a>(
    sample: &[&'a [u8]],
    bounds: &Boundaries,
    w: usize,
) -> Vec<&'a [u8]> {
    let mut out = Vec::new();
    for (i, &read) in sample.iter().enumerate() {
        let (b5, b3) = (bounds.five[i], bounds.three[i]);
        if b3 < b5 + 4 * w {
            continue;
        }
        for window in [&read[b5 + w..b5 + 2 * w], &read[b3 - 2 * w..b3 - w]] {
            if is_plain_acgt(window) {
                out.push(window);
            }
        }
    }
    out
}

/// Advances each boundary past the innermost hit of `patterns`, on either
/// strand, that starts within `ANCHOR_SLACK` bases of it, and with `repeat`
/// past every further hit the window holds. With `eroded`, a pattern may
/// hang off the boundary by up to its edit budget in overhang cost, as an
/// eroded adapter does at the physical read end, with at least
/// `MIN_OVERLAP` aligned bases. Only reads flagged in `active` are searched.
/// Returns the reads whose boundary moved, or an empty vector when none did.
pub(super) fn advance_boundaries(
    sample: &[&[u8]],
    bounds: &mut Boundaries,
    patterns: &[Vec<u8>],
    error_rate: f64,
    active: &[bool],
    repeat: bool,
    eroded: bool,
) -> Vec<bool> {
    let (five_w, three_w) = layer_windows(sample, bounds, 2 * WINDOW_LEN);
    let five_w: EndWindows = five_w.into_iter().filter(|(i, _)| active[*i]).collect();
    let three_w: EndWindows = three_w.into_iter().filter(|(i, _)| active[*i]).collect();
    let five_texts: Vec<&[u8]> = five_w.iter().map(|(_, w)| *w).collect();
    let three_texts: Vec<&[u8]> = three_w.iter().map(|(_, w)| *w).collect();
    let mut five_next = bounds.five.clone();
    let mut three_next = bounds.three.clone();
    // Hit spans per window, measured inward from the boundary. The inner edge
    // of each span drops the clipped bases of its alignment (`refine`), as a
    // trim does, so a layer the pattern does not truly cover is left for the
    // next layer's assembly.
    let mut five_hits: Vec<Vec<(usize, usize)>> = vec![Vec::new(); five_texts.len()];
    let mut three_hits: Vec<Vec<(usize, usize)>> = vec![Vec::new(); three_texts.len()];
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    // Patterns of one length, such as a barcode panel, share one tiled
    // search per window; other patterns are searched one strand at a time.
    // Windows without a hit are then searched with overhang when eroded.
    let mut singletons: Vec<Vec<u8>> = Vec::new();
    let mut by_length: std::collections::BTreeMap<usize, Vec<Vec<u8>>> =
        std::collections::BTreeMap::new();
    for pattern in patterns {
        by_length
            .entry(pattern.len())
            .or_default()
            .push(pattern.clone());
    }
    for (len, group) in by_length {
        let k = edit_budget(error_rate, len);
        if group.len() > 1 && len <= crate::adapter::search::MAX_TILED_PATTERN_LEN {
            let encoded = crate::adapter::search::encode_patterns(&group);
            for (text, hits) in five_texts.iter().zip(&mut five_hits) {
                let reversed: Vec<u8> = text.iter().rev().copied().collect();
                crate::adapter::search::encoded_pattern_hits(
                    &mut searcher,
                    &encoded,
                    text,
                    &reversed,
                    k,
                    |_, h| hits.push((h.start, h.end - h.clip_end)),
                );
            }
            for (text, hits) in three_texts.iter().zip(&mut three_hits) {
                let reversed: Vec<u8> = text.iter().rev().copied().collect();
                crate::adapter::search::encoded_pattern_hits(
                    &mut searcher,
                    &encoded,
                    text,
                    &reversed,
                    k,
                    |_, h| hits.push((text.len() - h.end, text.len() - h.start - h.clip_start)),
                );
            }
            continue;
        }
        for pattern in group {
            for strand in [
                pattern.clone(),
                crate::adapter::reverse_complement(&pattern),
            ] {
                for hit in searcher.search_texts(&strand, &five_texts, k) {
                    let (_, clip_end) = crate::adapter::refine::clips(&hit, false);
                    five_hits[hit.text_idx].push((hit.text_start, hit.text_end - clip_end));
                }
                for hit in searcher.search_texts(&strand, &three_texts, k) {
                    let len = three_texts[hit.text_idx].len();
                    let (clip_start, _) = crate::adapter::refine::clips(&hit, false);
                    three_hits[hit.text_idx]
                        .push((len - hit.text_end, len - hit.text_start - clip_start));
                }
            }
            singletons.push(pattern);
        }
    }
    if eroded && !singletons.is_empty() {
        let mut partial = crate::adapter::search::new_overhang_searcher(error_rate as f32);
        let missing = |hits: &[Vec<(usize, usize)>]| -> Vec<usize> {
            (0..hits.len()).filter(|&i| hits[i].is_empty()).collect()
        };
        let five_missing = missing(&five_hits);
        let three_missing = missing(&three_hits);
        let five_subset: Vec<&[u8]> = five_missing.iter().map(|&i| five_texts[i]).collect();
        let three_subset: Vec<&[u8]> = three_missing.iter().map(|&i| three_texts[i]).collect();
        for pattern in &singletons {
            let k = edit_budget(error_rate, pattern.len());
            for hit in partial.search_texts(pattern, &five_subset, k) {
                if hit.text_end - hit.text_start >= crate::adapter::MIN_OVERLAP {
                    let rc = hit.strand == sassy::Strand::Rc;
                    let (_, clip_end) = crate::adapter::refine::clips(&hit, rc);
                    five_hits[five_missing[hit.text_idx]]
                        .push((hit.text_start, hit.text_end - clip_end));
                }
            }
            for hit in partial.search_texts(pattern, &three_subset, k) {
                let len = three_subset[hit.text_idx].len();
                if hit.text_end - hit.text_start >= crate::adapter::MIN_OVERLAP {
                    let rc = hit.strand == sassy::Strand::Rc;
                    let (clip_start, _) = crate::adapter::refine::clips(&hit, rc);
                    three_hits[three_missing[hit.text_idx]]
                        .push((len - hit.text_end, len - hit.text_start - clip_start));
                }
            }
        }
    }
    let advance = |hits: &[(usize, usize)]| {
        let mut depth = 0;
        loop {
            let next = hits
                .iter()
                .filter(|&&(start, end)| start <= depth + ANCHOR_SLACK && end > depth)
                .map(|&(_, end)| end)
                .max();
            match next {
                Some(end) if repeat => depth = end,
                Some(end) => return end,
                None => return depth,
            }
        }
    };
    for ((i, _), hits) in five_w.iter().zip(&five_hits) {
        five_next[*i] = bounds.five[*i] + advance(hits);
    }
    for ((i, _), hits) in three_w.iter().zip(&three_hits) {
        three_next[*i] = bounds.three[*i] - advance(hits);
    }
    let moved: Vec<bool> = (0..sample.len())
        .map(|i| five_next[i] != bounds.five[i] || three_next[i] != bounds.three[i])
        .collect();
    bounds.five = five_next;
    bounds.three = three_next;
    if moved.contains(&true) {
        moved
    } else {
        Vec::new()
    }
}

/// Returns the median depth of the inner edge of `pattern` hits in physical
/// end windows of `end`, and the number of windows with a hit.
pub(super) fn inner_depth(
    searcher: &mut AmbiguousSearcher,
    pattern: &[u8],
    windows: &[&[u8]],
    end: End,
    error_rate: f64,
) -> (usize, usize) {
    let k = edit_budget(error_rate, pattern.len());
    let mut best: Vec<Option<(i32, usize)>> = vec![None; windows.len()];
    for hit in searcher.search_texts(pattern, windows, k) {
        let depth = match end {
            End::Five => hit.text_end,
            End::Three => windows[hit.text_idx].len() - hit.text_start,
        };
        let entry = &mut best[hit.text_idx];
        if entry.is_none_or(|(cost, _)| hit.cost < cost) {
            *entry = Some((hit.cost, depth));
        }
    }
    let mut depths: Vec<usize> = best.iter().flatten().map(|&(_, d)| d).collect();
    depths.sort_unstable();
    let present = depths.len();
    (depths.get(present / 2).copied().unwrap_or(0), present)
}

/// Returns the median depth of the outer edge of `pattern` hits in physical
/// end windows of `end`, or `usize::MAX` without hits.
pub(super) fn outer_depth(
    searcher: &mut AmbiguousSearcher,
    pattern: &[u8],
    windows: &[&[u8]],
    end: End,
    error_rate: f64,
) -> usize {
    let k = edit_budget(error_rate, pattern.len());
    let mut best: Vec<Option<(i32, usize)>> = vec![None; windows.len()];
    for hit in searcher.search_texts(pattern, windows, k) {
        let depth = match end {
            End::Five => hit.text_start,
            End::Three => windows[hit.text_idx].len() - hit.text_end,
        };
        let entry = &mut best[hit.text_idx];
        if entry.is_none_or(|(cost, _)| hit.cost < cost) {
            *entry = Some((hit.cost, depth));
        }
    }
    let mut depths: Vec<usize> = best.iter().flatten().map(|&(_, d)| d).collect();
    depths.sort_unstable();
    if depths.len() < MIN_SUPPORT_WINDOWS {
        return usize::MAX;
    }
    depths[depths.len() / 2]
}

/// Cuts the insert-facing part of `seq` whose reverse complement occurs at
/// the opposite physical end at a different depth. Returns the outer part,
/// or `None` when fewer than `MIN_PATTERN_LEN` outer bases remain, with the
/// removed insert stretch. A mirror at the same depth, or no mirror, leaves
/// `seq` unchanged.
pub(super) fn symmetry_cut(
    seq: &[u8],
    end: End,
    own: &[&[u8]],
    opposite: &[&[u8]],
    error_rate: f64,
) -> (Option<Vec<u8>>, Vec<u8>) {
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let (own_depth, present) = inner_depth(&mut searcher, seq, own, end, error_rate);
    if present < MIN_SUPPORT_WINDOWS {
        return (Some(seq.to_vec()), Vec::new());
    }
    let opposite_end = match end {
        End::Five => End::Three,
        End::Three => End::Five,
    };
    // Windows slide from the inner edge outward. A window whose mirror
    // lies at a different depth is insert, and the cut extends to its outer
    // edge; a window mirrored at the same depth is technical and ends the
    // scan. Windows without a recurrent mirror decide nothing, so a few
    // unmatched bases at the inner tip do not hide the insert behind them.
    let mut cut = 0;
    let mut offset = 0;
    while offset + MIRROR_SEGMENT <= seq.len() {
        let segment = match end {
            End::Five => &seq[seq.len() - offset - MIRROR_SEGMENT..seq.len() - offset],
            End::Three => &seq[offset..offset + MIRROR_SEGMENT],
        };
        let mirror = crate::adapter::reverse_complement(segment);
        let (depth, found) = inner_depth(
            &mut searcher,
            &mirror,
            opposite,
            opposite_end,
            MIRROR_ERROR_RATE * error_rate,
        );
        tracing::debug!(sequence = %String::from_utf8_lossy(seq), offset, own_depth, present, depth, found, "Mirror");
        if found >= MIN_SUPPORT_WINDOWS.max(present / MIRROR_SUPPORT_DIVISOR) {
            if depth.abs_diff(own_depth.saturating_sub(offset)) <= SYMMETRY_TOLERANCE {
                break;
            }
            cut = offset + MIRROR_SEGMENT;
        }
        offset += MIRROR_STEP;
    }
    if cut == 0 {
        return (Some(seq.to_vec()), Vec::new());
    }
    // The removed stretch is insert up to the last mirrored window, which
    // may straddle the junction with the technical sequence.
    let core = cut - MIRROR_SEGMENT;
    let (outer, insert) = match end {
        End::Five => (&seq[..seq.len() - cut], &seq[seq.len() - core..]),
        End::Three => (&seq[cut..], &seq[..core]),
    };
    let outer = (outer.len() >= MIN_PATTERN_LEN).then(|| outer.to_vec());
    (outer, insert.to_vec())
}

/// Removes from each candidate the prefix and suffix that occur in a
/// better supported candidate of another family in the same layer and end.
/// Sequence shared with a better supported neighbour, such as the flank
/// around a barcode, belongs to that neighbour's layer; a fragment of the
/// candidate itself is never stronger. Candidates left shorter than
/// `MIN_PATTERN_LEN` are dropped.
pub(super) fn strip_shared_ends(candidates: &mut Vec<Candidate>, error_rate: f64) {
    // Longest prefix of `a` that occurs anywhere in `b`.
    let shared_prefix = |a: &[u8], b: &[u8]| {
        (SHARED_END_MIN..=a.len())
            .rev()
            .find(|&n| b.windows(n).any(|w| w == &a[..n]))
            .unwrap_or(0)
    };
    let originals: Vec<(Vec<u8>, f64)> = candidates.iter().map(|c| (c.0.clone(), c.1)).collect();
    for (i, candidate) in candidates.iter_mut().enumerate() {
        let (seq, support) = &originals[i];
        let reversed: Vec<u8> = seq.iter().rev().copied().collect();
        let mut prefix = 0;
        let mut suffix = 0;
        for (j, (other, other_support)) in originals.iter().enumerate() {
            if j == i || other_support <= support || same_family(seq, other, error_rate) {
                continue;
            }
            prefix = prefix.max(shared_prefix(seq, other));
            let other_reversed: Vec<u8> = other.iter().rev().copied().collect();
            suffix = suffix.max(shared_prefix(&reversed, &other_reversed));
        }
        let start = if prefix >= SHARED_END_MIN { prefix } else { 0 };
        let end = if suffix >= SHARED_END_MIN {
            seq.len().saturating_sub(suffix).max(start)
        } else {
            seq.len()
        };
        candidate.0 = seq[start..end].to_vec();
    }
    candidates.retain(|c| c.0.len() >= MIN_PATTERN_LEN);
}

/// Resolves a variable layer in front of a constant one. When the best
/// supported candidate lies at least `MIN_PATTERN_LEN` bases past the
/// boundary and every candidate at the boundary is `RUN_JUMP` times rarer,
/// the stretch before it holds one member per read of a variable layer,
/// such as a barcode set between its flanks. The members are clustered
/// from the reads and replace the assembled fragments that started inside
/// the gap. Returns the candidates and the member sequences.
pub(super) fn with_variable_layer(
    searcher: &mut AmbiguousSearcher,
    candidates: Vec<Candidate>,
    windows: &[&[u8]],
    sample: &[&[u8]],
    end: End,
    error_rate: f64,
) -> (Vec<Candidate>, Vec<Vec<u8>>) {
    let depths: Vec<usize> = candidates
        .iter()
        .map(|c| outer_depth(searcher, &c.0, sample, end, error_rate))
        .collect();
    let shallow = candidates
        .iter()
        .zip(&depths)
        .filter(|(_, depth)| **depth < MIN_PATTERN_LEN)
        .map(|(c, _)| c.1)
        .fold(0.0, f64::max);
    let anchor = candidates
        .iter()
        .zip(&depths)
        .enumerate()
        .filter(|(_, (_, depth))| **depth >= MIN_PATTERN_LEN && **depth != usize::MAX)
        .max_by(|a, b| a.1.0.1.total_cmp(&b.1.0.1).then(b.0.cmp(&a.0)))
        .map(|(i, _)| i);
    let Some(anchor) = anchor.filter(|&i| candidates[i].1 >= shallow * f64::from(RUN_JUMP)) else {
        return (candidates, Vec::new());
    };
    let anchor_depth = depths[anchor];
    let members = variable_layer(searcher, &candidates[anchor].0, windows, end, error_rate);
    if members.is_empty() {
        return (candidates, Vec::new());
    }
    tracing::debug!(anchor = %String::from_utf8_lossy(&candidates[anchor].0), anchor_depth, members = members.len(), "Variable layer");
    let sequences: Vec<Vec<u8>> = members.iter().map(|c| c.0.clone()).collect();
    let kept: Vec<Candidate> = candidates
        .into_iter()
        .zip(depths)
        .filter(|(_, depth)| depth + MIN_PATTERN_LEN > anchor_depth)
        .map(|(c, _)| c)
        .chain(members)
        .collect();
    (kept, sequences)
}

/// Clusters the sequence between the boundary and the best `anchor` hit of
/// each window into supported families. Each family is a candidate with
/// its member count as support and weight. `windows` is every window of the
/// layer rather than the `RECOUNT_WINDOWS` validation sample: a member held
/// by `VARIABLE_MEMBER_SUPPORT` of the reads reaches the `MIN_SUPPORT_WINDOWS`
/// floor only in the complete sample.
pub(super) fn variable_layer(
    searcher: &mut AmbiguousSearcher,
    anchor: &[u8],
    windows: &[&[u8]],
    end: End,
    error_rate: f64,
) -> Vec<Candidate> {
    let k = edit_budget(error_rate, anchor.len());
    let mut best: Vec<Option<(i32, usize, usize, usize)>> = vec![None; windows.len()];
    for hit in searcher.search_texts(anchor, windows, k) {
        let depth = match end {
            End::Five => hit.text_start,
            End::Three => windows[hit.text_idx].len() - hit.text_end,
        };
        let entry = &mut best[hit.text_idx];
        if entry.is_none_or(|(cost, old, _, _)| (hit.cost, depth) < (cost, old)) {
            *entry = Some((hit.cost, depth, hit.text_start, hit.text_end));
        }
    }
    let gaps: Vec<&[u8]> = windows
        .iter()
        .zip(&best)
        .filter_map(|(window, hit)| {
            hit.map(|(_, _, start, stop)| match end {
                End::Five => &window[..start],
                End::Three => &window[stop..],
            })
        })
        .filter(|gap| gap.len() >= MIN_PATTERN_LEN)
        .collect();
    let floor =
        MIN_SUPPORT_WINDOWS.max((windows.len() as f64 * VARIABLE_MEMBER_SUPPORT).ceil() as usize);
    cluster_sequences(searcher, &gaps, error_rate, floor)
        .into_iter()
        .map(|(seq, members)| {
            let support = members as f64 / windows.len() as f64;
            let weight = members as u64 * (seq.len().saturating_sub(KMER_K) + 1) as u64;
            (seq, support, false, weight, false)
        })
        .collect()
}

/// Groups `sequences` into families by edit distance. Each family is seeded
/// by the most frequent exact sequence left, polished by its members, and
/// kept when it has at least `floor` members. Returns `(consensus, members)`.
pub(super) fn cluster_sequences(
    searcher: &mut AmbiguousSearcher,
    sequences: &[&[u8]],
    error_rate: f64,
    floor: usize,
) -> Vec<(Vec<u8>, usize)> {
    let mut unassigned: Vec<usize> = (0..sequences.len()).collect();
    let mut out = Vec::new();
    while unassigned.len() >= floor && out.len() < MAX_ADAPTERS_PER_END {
        let mut counts: std::collections::HashMap<&[u8], usize> = std::collections::HashMap::new();
        for &i in &unassigned {
            *counts.entry(sequences[i]).or_default() += 1;
        }
        let Some((&seed, _)) = counts.iter().max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0))) else {
            break;
        };
        let texts: Vec<&[u8]> = unassigned.iter().map(|&i| sequences[i]).collect();
        let hits = windows_with(searcher, seed, &texts, edit_budget(error_rate, seed.len()));
        let members: Vec<&[u8]> = texts
            .iter()
            .zip(&hits)
            .filter(|(_, hit)| **hit)
            .map(|(&text, _)| text)
            .collect();
        if members.len() >= floor {
            let consensus = polish_consensus(seed, &members);
            tracing::debug!(seed = %String::from_utf8_lossy(seed), consensus = %String::from_utf8_lossy(&consensus), members = members.len(), "Variable layer member");
            out.push((consensus, members.len()));
        }
        unassigned = unassigned
            .into_iter()
            .zip(hits)
            .filter(|(_, hit)| !hit)
            .map(|(i, _)| i)
            .collect();
    }
    out
}
