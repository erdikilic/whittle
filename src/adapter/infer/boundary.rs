//! Insert boundaries of assembled candidates: primer reconstruction at
//! conserved insert starts, insert-facing termination, and the support runs of
//! a path.

use super::*;

/// Primer families reconstructed at one conserved insert boundary.
pub(super) struct PrimerBoundary {
    /// Supported primer sequences, oriented as they occur in the read.
    pub(super) primers: Vec<Vec<u8>>,
    /// First `START_K` bases of the conserved insert, in read orientation.
    pub(super) insert_start: Vec<u8>,
}

/// Uses recurrent unprimed read starts to locate a conserved insert boundary.
/// A primer must be independently supported upstream of that boundary; the
/// conserved insert itself is excluded from the inferred trimming sequence.
/// `windows` read outward from the `end` boundary, reversed for the 3' end.
/// Read starts that open alternative layers behind the reconstructed sequence
/// (`opens_layers`) mark no boundary, and neither do read starts outward of
/// them, which lie in that sequence.
pub(super) fn contrast_boundary(
    consensus: &[u8],
    windows: &[&[u8]],
    end: End,
    error_rate: f64,
) -> Option<PrimerBoundary> {
    if consensus.len() < 2 * KMER_K {
        return None;
    }
    let reverse = end == End::Three;
    let mut seq: Vec<u8> = if reverse {
        consensus.iter().rev().copied().collect()
    } else {
        consensus.to_vec()
    };
    let mut starts = std::collections::HashMap::<u64, usize>::new();
    for window in windows {
        if window.len() >= START_K
            && let Some(code) = encode_kmer(&window[..START_K])
        {
            *starts.entry(code).or_default() += 1;
        }
    }
    let floor = 10.max(windows.len() / 32);
    if starts.values().all(|&count| count < floor) {
        return None;
    }
    let positions_in = |seq: &[u8]| -> Vec<(usize, usize)> {
        seq.windows(START_K)
            .enumerate()
            .take(seq.len().saturating_sub(2 * KMER_K) + usize::from(seq.len() >= 2 * KMER_K))
            .map(|(pos, word)| {
                (
                    pos,
                    encode_kmer(word)
                        .and_then(|code| starts.get(&code).copied())
                        .unwrap_or(0),
                )
            })
            .filter(|&(_, count)| count >= floor)
            .collect()
    };
    let mut positions = positions_in(&seq);
    if positions.is_empty() && seq.len() >= KMER_K {
        let mut prefix = upstream_consensus(&seq[..KMER_K], windows);
        if !prefix.is_empty() {
            prefix.extend_from_slice(&seq);
            seq = prefix;
            positions = positions_in(&seq);
        }
    }
    positions.sort_by_key(|&(pos, count)| (std::cmp::Reverse(count), pos));
    let mut tried = std::collections::HashSet::new();
    // Position of the innermost read start that opens layers. The read
    // starts outward of it lie in the sequence in front of those layers.
    let mut opened: Option<usize> = None;
    for (pos, _) in positions {
        if opened.is_some_and(|at| pos < at) {
            continue;
        }
        let anchor = &seq[pos..pos + KMER_K];
        if tried.len() >= MAX_BOUNDARY_ANCHORS || !tried.insert(anchor.to_vec()) {
            continue;
        }
        let primer = upstream_consensus(anchor, windows);
        let mut primers = if primer.len() >= MIN_PATTERN_LEN
            && primer
                .iter()
                .filter(|&&b| !matches!(b, b'A' | b'C' | b'G' | b'T'))
                .count()
                * 4
                <= primer.len()
        {
            vec![primer]
        } else {
            // Distinct primer families can share an insert boundary. Cluster
            // their adjacent sequence before extending each family separately.
            let mut searcher = crate::adapter::search::new_searcher_fwd();
            let mut best = vec![None; windows.len()];
            for hit in searcher.search_texts(anchor, windows, 2) {
                let key = (hit.cost, hit.text_start);
                let entry = &mut best[hit.text_idx];
                if entry.is_none_or(|old| key < old) {
                    *entry = Some(key);
                }
            }
            let mut groups = std::collections::BTreeMap::<u64, Vec<&[u8]>>::new();
            for (&window, hit) in windows.iter().zip(best) {
                if let Some((_, start)) = hit.filter(|&(_, start)| start >= 6) {
                    let key = encode_kmer(&window[start - 6..start]).unwrap();
                    groups.entry(key).or_default().push(window);
                }
            }
            let mut groups: Vec<_> = groups
                .into_iter()
                .filter(|(_, group)| {
                    group.len()
                        >= MIN_SUPPORT_WINDOWS
                            .max((windows.len() as f64 * KEEP_SUPPORT).ceil() as usize)
                })
                .collect();
            groups.sort_by_key(|(key, group)| (std::cmp::Reverse(group.len()), *key));
            groups
                .into_iter()
                .take(MAX_ADAPTERS_PER_END)
                .map(|(_, group)| upstream_consensus(anchor, &group))
                .filter(|p| {
                    p.len() >= MIN_PATTERN_LEN
                        && p.iter()
                            .filter(|&&b| !matches!(b, b'A' | b'C' | b'G' | b'T'))
                            .count()
                            * 4
                            <= p.len()
                })
                .collect()
        };
        // Reads that lack an outer adapter start at each of the layers behind
        // it, so a start that opens one of them is no insert boundary.
        if let Some(outer) = primers
            .iter()
            .find(|primer| opens_layers(primer, windows, &starts, floor, error_rate))
        {
            let outer: Vec<u8> = if reverse {
                outer.iter().rev().copied().collect()
            } else {
                outer.clone()
            };
            tracing::debug!(sequence = %String::from_utf8_lossy(&outer), "Read starts open layers behind a sequence");
            opened = Some(opened.map_or(pos, |at| at.max(pos)));
            continue;
        }
        if reverse {
            for primer in &mut primers {
                primer.reverse();
            }
        }
        if !primers.is_empty() {
            let mut insert_start = anchor[..START_K].to_vec();
            if reverse {
                insert_start.reverse();
            }
            return Some(PrimerBoundary {
                primers,
                insert_start,
            });
        }
    }
    None
}

/// Returns whether the read starts behind `primer`, a sequence reconstructed
/// outward of a recurrent read start, open alternative layers rather than one
/// insert: the windows that hold the inner `START_K` bases of `primer` within
/// its length of their outer end continue with two recurrent read starts
/// (`starts` holding them at least `floor` times) that are not one sequence
/// within `error_rate`, each in at least `MIN_SUPPORT_WINDOWS` windows. Reads
/// that lack an outer adapter start at the layer behind it, so an adapter
/// followed by either of two primers, as a cDNA adapter is by the primers of
/// the two strands, precedes two read starts; a primer precedes its insert.
fn opens_layers(
    primer: &[u8],
    windows: &[&[u8]],
    starts: &std::collections::HashMap<u64, usize>,
    floor: usize,
    error_rate: f64,
) -> bool {
    let Some(inner) = primer.len().checked_sub(START_K).map(|at| &primer[at..]) else {
        return false;
    };
    let mut following = std::collections::HashMap::<u64, usize>::new();
    for window in windows {
        let next = window
            .windows(2 * START_K)
            .take(primer.len() + 1)
            .find(|pair| &pair[..START_K] == inner)
            .and_then(|pair| encode_kmer(&pair[START_K..]))
            .filter(|code| starts.get(code).is_some_and(|&count| count >= floor));
        if let Some(code) = next {
            *following.entry(code).or_default() += 1;
        }
    }
    let mut layers: Vec<(usize, u64)> = following
        .into_iter()
        .filter(|&(_, count)| count >= MIN_SUPPORT_WINDOWS)
        .map(|(code, count)| (count, code))
        .collect();
    layers.sort_unstable_by(|a, b| b.cmp(a));
    let Some(&(_, first)) = layers.first() else {
        return false;
    };
    let first = decode_kmer(first, START_K);
    layers[1..]
        .iter()
        .any(|&(_, code)| !same_adapter(&first, &decode_kmer(code, START_K), error_rate))
}

/// Tests whether a candidate is predominantly on the insert-facing side of
/// independently reconstructed primers.
pub(super) fn predominantly_insert(
    seq: &[u8],
    windows: &[&[u8]],
    end: End,
    edits: usize,
    primer_edges: &[Option<(usize, usize)>],
    searcher: &mut AmbiguousSearcher,
) -> bool {
    let mut present = vec![false; windows.len()];
    let mut downstream = vec![false; windows.len()];
    for hit in searcher.search_texts(seq, windows, edits) {
        present[hit.text_idx] = true;
        if let Some((edge, overlap)) = primer_edges[hit.text_idx] {
            downstream[hit.text_idx] |= match end {
                End::Five => hit.text_start + overlap >= edge,
                End::Three => hit.text_end <= edge + overlap,
            };
        }
    }
    let n = present.iter().filter(|&&p| p).count();
    let paired = downstream.iter().filter(|&&p| p).count();
    paired >= MIN_SUPPORT_WINDOWS && paired * 2 >= n
}

/// Checks insert-facing continuation outside the bounded assembly windows.
/// Each window contributes at most once to each continuation count.
pub(super) fn supported_termination(word: &[u8], windows: &[&[u8]], end: End) -> bool {
    let mut present = 0;
    let mut extensions = [0; 4];
    for window in windows {
        let mut seen = false;
        let mut extended = [false; 4];
        for (pos, matched) in window.windows(word.len()).enumerate() {
            if matched != word {
                continue;
            }
            let next = match end {
                End::Five => window.get(pos + word.len()),
                End::Three => pos.checked_sub(1).and_then(|i| window.get(i)),
            };
            if let Some(code) = next.and_then(|b| encode_kmer(std::slice::from_ref(b))) {
                seen = true;
                extended[code as usize] = true;
            }
        }
        present += usize::from(seen);
        for (count, seen) in extensions.iter_mut().zip(extended) {
            *count += usize::from(seen);
        }
    }
    present > 0 && extensions.into_iter().max().unwrap_or(0) * 2 <= present
}

/// Ratio between the support plateaus on either side of a layer boundary,
/// such as a barcode joining its shared flank. A degenerate primer base
/// halves the support of the k-mers over it and does not reach it.
pub(super) const RUN_JUMP: u32 = 4;

/// K-mers of consistent support on each side of a boundary. A support ramp
/// without a plateau, such as an eroded adapter start, is one run.
pub(super) const PLATEAU: usize = 8;

/// K-mers spanning a boundary whose support ramps between the plateaus,
/// as k-mers overlapping both a barcode and its flank do.
pub(super) const RAMP: usize = 4;

/// K-mers over which the changed support level must persist. A
/// sequencing-error dip spans at most `KMER_K` k-mers; a barcode spans more.
pub(super) const RUN_PERSIST: usize = 24;

/// K-mers over which each continuation of a division keeps its support. A
/// branch opened by one base, such as a degenerate primer base, rejoins the
/// path after `KMER_K`.
pub(super) const DIVISION_PERSIST: usize = KMER_K + 1;

/// Returns the base span `[lo, hi)` of the outermost support run of a
/// consensus and whether the run ends at a rise in support. `weights` holds
/// the support of each k-mer by start base. Scanning inward from the read
/// end, the run ends at the first boundary between two flat plateaus whose
/// supports differ `RUN_JUMP`-fold, where the new level persists for
/// `RUN_PERSIST` k-mers or to the end. Within the run, k-mers below
/// `BOUNDARY_SUPPORT` of the run peak are trimmed from both sides.
pub(super) fn supported_span(weights: &[u32], end: End) -> (usize, usize, bool) {
    let n = weights.len();
    let order: Vec<usize> = match end {
        End::Five => (0..n).collect(),
        End::Three => (0..n).rev().collect(),
    };
    let first = order
        .iter()
        .position(|&j| weights[j] >= MIN_SUPPORT_WINDOWS as u32)
        .unwrap_or(0);
    let level = |idx: &[usize]| {
        let mut values: Vec<u32> = idx.iter().map(|&j| weights[j]).collect();
        values.sort_unstable();
        values[values.len() / 2]
    };
    let flat = |idx: &[usize]| {
        let (min, max) = idx
            .iter()
            .map(|&j| weights[j])
            .fold((u32::MAX, 0), |(lo, hi), w| (lo.min(w), hi.max(w)));
        max <= min.saturating_mul(2)
    };
    let mut last = order.len() - 1;
    let mut rises = false;
    for j in first + PLATEAU..order.len() {
        let left = &order[j - PLATEAU..j];
        let right_start = j + RAMP;
        let right_end = (right_start + PLATEAU).min(order.len());
        if right_end < right_start + RAMP {
            break;
        }
        let right = &order[right_start..right_end];
        if !flat(left) || !flat(right) {
            continue;
        }
        let (old, new) = (level(left), level(right));
        let drop = new.saturating_mul(RUN_JUMP) <= old;
        let rise = old.saturating_mul(RUN_JUMP) <= new;
        if !(drop || rise) {
            continue;
        }
        let persists = order[right_start..(right_start + RUN_PERSIST).min(order.len())]
            .iter()
            .all(|&k| {
                if drop {
                    weights[k].saturating_mul(RUN_JUMP) <= old
                } else {
                    weights[k] >= old.saturating_mul(RUN_JUMP)
                }
            });
        if !persists {
            continue;
        }
        // The old level ends at the last k-mer within twofold of it.
        let mut split = j;
        while split < right_start && {
            let w = weights[order[split]];
            if drop {
                w.saturating_mul(2) >= old
            } else {
                w <= old.saturating_mul(2)
            }
        } {
            split += 1;
        }
        last = split - 1;
        rises = rise;
        break;
    }
    // The run is never empty, and its peak clears the floor.
    let run = &order[first..=last];
    let peak = run.iter().map(|&j| weights[j]).max().unwrap_or(0);
    let floor = (peak as f64 * BOUNDARY_SUPPORT) as u32;
    let supported = run.iter().copied().filter(|&j| weights[j] >= floor);
    let (mut lo, mut hi) = (usize::MAX, 0);
    for j in supported {
        lo = lo.min(j);
        hi = hi.max(j + KMER_K);
    }
    (lo, hi, rises)
}

/// Returns the inner edges of the shared parts of a consensus whose path
/// divides into supported continuations, as a shared adapter divides into
/// the strand-specific primers behind it, with the first `KMER_K` bases of
/// two continuations, inward from the read end: the path's own and another
/// one, once for each other continuation. The edge is the end of the last
/// shared k-mer for the 5' end and its start for the 3' end.
///
/// Scanning inward over the k-mers of `[lo, hi)`, the innermost included,
/// whose continuations then lie past the span, a k-mer divides when it
/// lies within `RUN_JUMP` of the peak of the span, the path's own successor
/// and another successor each hold at least `1 / RUN_JUMP` of its support
/// and together at most twice it, and both open paths that keep `1 /
/// RUN_JUMP` of their opening support for `DIVISION_PERSIST` k-mers, the
/// other one off the consensus.
/// A degenerate primer base opens a branch that rejoins the consensus within
/// `KMER_K` k-mers and does not divide it. The shared part holds at least
/// the k-mer that divides, so it is longer than `MIN_PATTERN_LEN`. `weights`
/// holds the window counts of k-mers (`kmer_counts`) over windows that reach
/// past the assembly windows, so a division near the inner end of the
/// assembly windows still shows its continuations. What the path divides
/// into is read from the reads (`division_kind`).
pub(super) fn division_point(
    cons: &[u8],
    weights: &KmerMap<(u32, usize)>,
    lo: usize,
    hi: usize,
    end: End,
) -> Vec<(usize, [Vec<u8>; 2])> {
    let mut divisions = Vec::new();
    if hi < lo + KMER_K + 1 {
        return divisions;
    }
    let codes: Vec<u64> = cons
        .windows(KMER_K)
        .map(|w| encode_kmer(w).unwrap_or(u64::MAX))
        .collect();
    let on_path: std::collections::HashSet<u64> = codes.iter().copied().collect();
    let weight = |code: u64| weights.get(&code).map_or(0, |&(count, _)| count);
    let mask = (1u64 << (2 * KMER_K)) - 1;
    let step = |code: u64, base: u64| match end {
        End::Five => ((code << 2) | base) & mask,
        End::Three => (code >> 2) | (base << (2 * (KMER_K - 1))),
    };
    // Starts of the k-mers of the span, inward from the read end, that the
    // path continues past.
    let starts: Vec<usize> = match end {
        End::Five => (lo..=hi - KMER_K)
            .filter(|&j| j + 1 < codes.len())
            .collect(),
        End::Three => (lo.max(1)..=hi - KMER_K).rev().collect(),
    };
    let peak = starts.iter().map(|&j| weight(codes[j])).max().unwrap_or(0);
    for j in starts {
        let here = weight(codes[j]);
        let next = match end {
            End::Five => codes[j + 1],
            End::Three => codes[j - 1],
        };
        if here < MIN_SUPPORT_WINDOWS as u32
            || here.saturating_mul(RUN_JUMP) < peak
            || weight(next).saturating_mul(RUN_JUMP) < here
        {
            continue;
        }
        // A continuation persists when its heaviest path keeps at least
        // `1 / RUN_JUMP` of its opening support for `DIVISION_PERSIST`
        // k-mers, as each of two template starts that open with the same
        // base does once they part; a barcode panel behind a flank spreads
        // its support over more members instead.
        let persists = |start: u64, off_path: bool| {
            let floor = weight(start)
                .div_ceil(RUN_JUMP)
                .max(MIN_SUPPORT_WINDOWS as u32);
            let mut code = start;
            (0..DIVISION_PERSIST).all(|_| {
                let kept = weight(code) >= floor && !(off_path && on_path.contains(&code));
                code = (0..4)
                    .map(|base| step(code, base))
                    .max_by_key(|&c| weight(c))
                    .unwrap_or(code);
                kept
            })
        };
        if !persists(next, false) {
            continue;
        }
        // The first k-mer lying wholly in each continuation.
        let own = |start: u64| {
            let mut code = start;
            for _ in 1..KMER_K {
                code = (0..4)
                    .map(|base| step(code, base))
                    .max_by_key(|&c| weight(c))
                    .unwrap_or(code);
            }
            decode_kmer(code, KMER_K)
        };
        let edge = match end {
            End::Five => j + KMER_K,
            End::Three => j,
        };
        // Each other successor that opens a continuation is paired with the
        // path's own: behind a primer at both ends of a molecule, the starts
        // of several templates open more than two, and only a pair from two
        // templates tells them from the primers of two strands
        // (`division_kind`).
        for alt in (0..4).map(|base| step(codes[j], base)).filter(|&alt| {
            alt != next
                && weight(alt).saturating_mul(RUN_JUMP) >= here
                && weight(alt) + weight(next) <= here.saturating_mul(2)
                && persists(alt, true)
        }) {
            divisions.push((edge, [own(next), own(alt)]));
        }
    }
    divisions
}

/// What a path divides into at a division (`division_point`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Division {
    /// The primers of the two strands: the shared layer ends at the
    /// division, and each continuation opens a strand-specific layer.
    Strands,
    /// The starts of templates behind one primer: the insert begins at the
    /// division.
    Templates,
}

/// Reads what a path divides into from where each continuation recurs
/// reverse complemented at the other read end. The continuations open with
/// `words` behind the k-mer `shared`; `windows` are the windows of the end
/// assembled and `partners` the physical windows of the other end of the
/// same reads, empty where a read has none.
///
/// The primers of the two strands flank the same molecules, so the mirror
/// of each is read at the other end of the reads of the other
/// (`Division::Strands`). A template start behind one primer is read from
/// its other strand at the other end of reads that hold neither
/// continuation (`Division::Templates`): behind the other primer, or, with
/// one primer at both ends of a molecule, behind the same primer in reads
/// that open with a further continuation. A barcode read at both ends
/// recurs in its own reads and is neither. The reads holding each mirror,
/// at least `MIN_SUPPORT_WINDOWS` and `1 / MIRROR_SUPPORT_DIVISOR` of the
/// partners, must hold one kind `RUN_JUMP`-fold more often than the
/// others. The words must differ beyond their edit budget, as the variants
/// of one sequence do not. The primers of the two strands hold at least
/// `1 - 1 / RUN_JUMP` of the windows holding `shared`: the two orientations
/// of one of several templates behind a primer at both ends of a molecule
/// recur in each other's reads as well, but hold only that template's
/// share, and the starts of two of those templates divide the path into
/// templates. The primers of the two strands also open the layers of most
/// reads and continue in at least `1 / RUN_JUMP` of `windows`, which a path
/// specific to one barcode does not, while a shared layer that only some
/// reads keep whole, such as a rear adapter, does.
pub(super) fn division_kind(
    shared: &[u8],
    words: &[Vec<u8>; 2],
    windows: &[&[u8]],
    partners: &[&[u8]],
    error_rate: f64,
) -> Option<Division> {
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let edits = edit_budget(error_rate, KMER_K);
    if !hits(&mut searcher, &words[0], &words[1], edits).is_empty() {
        return None;
    }
    let jump = RUN_JUMP as usize;
    let at_shared = windows_with(&mut searcher, shared, windows, edits);
    let held = [
        windows_with(&mut searcher, &words[0], windows, edits),
        windows_with(&mut searcher, &words[1], windows, edits),
    ];
    let holding = at_shared.iter().filter(|&&at| at).count();
    let continued = (0..windows.len())
        .filter(|&j| at_shared[j] && (held[0][j] || held[1][j]))
        .count();
    // The primers of the two strands hold `1 - 1 / RUN_JUMP` of the windows
    // holding `shared` and `1 / RUN_JUMP` of all windows.
    let strands = continued * jump >= holding * (jump - 1) && continued * jump >= windows.len();
    let paired: Vec<usize> = (0..partners.len())
        .filter(|&j| !partners[j].is_empty())
        .collect();
    let texts: Vec<&[u8]> = paired.iter().map(|&j| partners[j]).collect();
    let floor = MIN_SUPPORT_WINDOWS.max(texts.len() / MIRROR_SUPPORT_DIVISOR);
    let mirror_edits = edit_budget(MIRROR_ERROR_RATE * error_rate, KMER_K);
    let mut kinds = [None; 2];
    for (index, word) in words.iter().enumerate() {
        let mirror = crate::adapter::reverse_complement(word);
        // Reads holding the mirror by what their own end holds: this
        // continuation, the other, neither, or both.
        let mirrored = windows_with(&mut searcher, &mirror, &texts, mirror_edits);
        let mut tally = [0usize; 4];
        for (&j, mirrored) in paired.iter().zip(mirrored) {
            if mirrored {
                tally[match (held[index][j], held[1 - index][j]) {
                    (true, false) => 0,
                    (false, true) => 1,
                    (false, false) => 2,
                    (true, true) => 3,
                }] += 1;
            }
        }
        let [own, other, neither, both] = tally;
        if own + both + other + neither < floor {
            return None;
        }
        kinds[index] = if other >= jump * (own + both + neither) && strands {
            Some(Division::Strands)
        } else if neither >= jump * (own + both + other) {
            Some(Division::Templates)
        } else {
            None
        };
    }
    match kinds {
        [Some(first), Some(second)] if first == second => Some(first),
        _ => None,
    }
}

/// Returns the length of the homopolymer run that ends `seq` on its
/// insert-facing side: the last bases for the 5' end, the first for the 3'
/// end.
pub(super) fn inner_run(seq: &[u8], end: End) -> usize {
    let mut bases: Box<dyn Iterator<Item = &u8>> = match end {
        End::Five => Box::new(seq.iter().rev()),
        End::Three => Box::new(seq.iter()),
    };
    let Some(&first) = bases.next() else {
        return 0;
    };
    1 + bases.take_while(|&&b| b == first).count()
}
