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
pub(super) fn contrast_boundary(
    consensus: &[u8],
    windows: &[&[u8]],
    end: End,
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
    for (pos, _) in positions {
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
    let run = &order[first..=last];
    let peak = run.iter().map(|&j| weights[j]).max().unwrap_or(0);
    let floor = (peak as f64 * BOUNDARY_SUPPORT) as u32;
    let supported = run.iter().copied().filter(|&j| weights[j] >= floor);
    let (mut lo, mut hi) = (usize::MAX, 0);
    for j in supported {
        lo = lo.min(j);
        hi = hi.max(j + KMER_K);
    }
    if lo == usize::MAX {
        (0, KMER_K.min(n + KMER_K - 1), rises)
    } else {
        (lo, hi, rises)
    }
}

/// Returns the inner edge of the shared part of a consensus whose path
/// divides into two supported continuations, as a shared adapter divides
/// into the strand-specific primers behind it, with the first `KMER_K`
/// bases of each continuation; `None` without a division. The edge is the
/// end of the last shared k-mer for the 5' end and its start for the 3' end.
///
/// Scanning inward over the k-mers of `[lo, hi)`, a k-mer divides when it
/// lies at the support level of the path (within `RUN_JUMP` of the peak of
/// the span, and in at least `1 / RUN_JUMP` of the `windows` counted), the
/// path's own successor and another successor each hold at least
/// `1 / RUN_JUMP` of its support and together at most twice it, and both
/// open paths that keep their support for `RUN_PERSIST` k-mers, the other
/// one off the consensus. A degenerate primer base opens a branch that
/// rejoins the consensus within `KMER_K` k-mers and does not divide it. The
/// shared part keeps at least `MIN_PATTERN_LEN` bases. `weights` counts
/// k-mers over windows that reach past the assembly windows (`kmer_counts`),
/// so a division near the inner end of the assembly windows still shows its
/// continuations.
pub(super) fn division_point(
    cons: &[u8],
    weights: &KmerMap<u32>,
    windows: usize,
    lo: usize,
    hi: usize,
    end: End,
) -> Option<(usize, [Vec<u8>; 2])> {
    if hi < lo + KMER_K + 1 {
        return None;
    }
    let codes: Vec<u64> = cons
        .windows(KMER_K)
        .map(|w| encode_kmer(w).unwrap_or(u64::MAX))
        .collect();
    let on_path: std::collections::HashSet<u64> = codes.iter().copied().collect();
    let weight = |code: u64| weights.get(&code).copied().unwrap_or(0);
    let mask = (1u64 << (2 * KMER_K)) - 1;
    let step = |code: u64, base: u64| match end {
        End::Five => ((code << 2) | base) & mask,
        End::Three => (code >> 2) | (base << (2 * (KMER_K - 1))),
    };
    // Starts of the k-mers of the span, inward from the read end.
    let starts: Vec<usize> = match end {
        End::Five => (lo..hi - KMER_K).collect(),
        End::Three => (lo + 1..=hi - KMER_K).rev().collect(),
    };
    let peak = starts.iter().map(|&j| weight(codes[j])).max().unwrap_or(0);
    for j in starts {
        let shared = match end {
            End::Five => j + KMER_K - lo,
            End::Three => hi - j,
        };
        if shared < MIN_PATTERN_LEN {
            continue;
        }
        let here = weight(codes[j]);
        let next = match end {
            End::Five => codes[j + 1],
            End::Three => codes[j - 1],
        };
        if here < MIN_SUPPORT_WINDOWS as u32
            || (here as usize).saturating_mul(RUN_JUMP as usize) < windows
            || here.saturating_mul(RUN_JUMP) < peak
            || weight(next).saturating_mul(RUN_JUMP) < here
        {
            continue;
        }
        // A continuation persists when its heaviest path keeps at least half
        // of its opening support for `RUN_PERSIST` k-mers; a barcode panel
        // behind a flank spreads its support over the members instead.
        let persists = |start: u64, off_path: bool| {
            let floor = weight(start).div_ceil(2).max(MIN_SUPPORT_WINDOWS as u32);
            let mut code = start;
            (0..RUN_PERSIST).all(|_| {
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
        let alt = (0..4).map(|base| step(codes[j], base)).find(|&alt| {
            alt != next
                && weight(alt).saturating_mul(RUN_JUMP) >= here
                && weight(alt) + weight(next) <= here.saturating_mul(2)
                && persists(alt, true)
        });
        if let Some(alt) = alt {
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
            return Some((edge, [own(next), own(alt)]));
        }
    }
    None
}

/// Returns the number of `windows` holding each k-mer of `KMER_K` plain
/// bases, as `top_kmers` counts them, without its ranking cut.
pub(super) fn kmer_counts(windows: &[&[u8]]) -> KmerMap<u32> {
    let mut counts: KmerMap<(u32, usize)> =
        KmerMap::with_capacity_and_hasher(windows.len() * 64, Default::default());
    let mask = u64::MAX >> (64 - 2 * KMER_K);
    for (window_index, window) in windows.iter().enumerate() {
        let (mut code, mut valid) = (0, 0);
        for &base in *window {
            if let Some(bits) = encode_kmer(&[base]) {
                code = ((code << 2) | bits) & mask;
                valid += 1;
                if valid >= KMER_K {
                    let entry = counts.entry(code).or_insert((0, usize::MAX));
                    if entry.1 != window_index {
                        entry.0 += 1;
                        entry.1 = window_index;
                    }
                }
            } else {
                (code, valid) = (0, 0);
            }
        }
    }
    counts
        .into_iter()
        .map(|(code, (count, _))| (code, count))
        .collect()
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
