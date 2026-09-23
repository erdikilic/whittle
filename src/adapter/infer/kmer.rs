//! Exact k-mers of read-end windows: packed encoding, the hasher and map keyed
//! by packed codes, repeat filters, the counts of the most frequent k-mers, and
//! the heaviest paths through their graph.

use super::*;

/// Encodes a k-mer at 2 bits per base (A=0, C=1, G=2, T=3). `None` when it
/// holds any byte other than uppercase ACGT (reads are uppercased upstream) or
/// is longer than 32 bases.
pub(super) fn encode_kmer(bytes: &[u8]) -> Option<u64> {
    if bytes.len() > 32 {
        return None;
    }
    let mut code = 0u64;
    for &b in bytes {
        let two = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        code = (code << 2) | two;
    }
    Some(code)
}

/// Decodes a 2-bit code back to its `k` bases; the inverse of `encode_kmer`.
pub(super) fn decode_kmer(mut code: u64, k: usize) -> Vec<u8> {
    let mut out = vec![0u8; k];
    for i in (0..k).rev() {
        out[i] = match code & 0b11 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
        code >>= 2;
    }
    out
}

/// Slices the first and last `w` bytes of each read into 5' and 3' window
/// lists. Empty reads are skipped, and so is any window holding a byte
/// outside ACGT: an uncalled base is evidence of nothing, and on the IUPAC
/// profile it would match every k-mer for free.
#[cfg(test)]
pub(super) fn end_windows<'a>(sample: &[&'a [u8]], w: usize) -> (Vec<&'a [u8]>, Vec<&'a [u8]>) {
    let (five, three) = layer_windows(sample, &Boundaries::new(sample), w);
    (
        five.into_iter().map(|(_, window)| window).collect(),
        three.into_iter().map(|(_, window)| window).collect(),
    )
}

/// Hashes packed k-mer codes by multiplicative mixing with a fold of the
/// high bits, which distributes the uniformly encoded keys across buckets.
#[derive(Default, Clone, Copy)]
pub(super) struct KmerHasher(u64);

impl std::hash::Hasher for KmerHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 = (self.0 ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
        }
    }

    fn write_u64(&mut self, code: u64) {
        let mixed = code.wrapping_mul(0x9e37_79b9_7f4a_7c15);
        self.0 = mixed ^ (mixed >> 32);
    }
}

/// Hash map keyed by packed k-mer codes.
pub(super) type KmerMap<V> =
    std::collections::HashMap<u64, V, std::hash::BuildHasherDefault<KmerHasher>>;

/// Returns whether a k-mer consists of a short tandem repeat.
pub(super) fn is_low_complexity(kmer: &[u8]) -> bool {
    (1..=8.min(kmer.len() / 2))
        .any(|period| kmer.iter().enumerate().all(|(i, &b)| b == kmer[i % period]))
}

/// Rejects consensuses dominated by a short approximate tandem repeat. Each
/// base is compared with the base one period earlier, so isolated insertions
/// and deletions shift the phase without hiding the repeat. Each period is
/// tested over at least ten comparisons so the threshold stays selective for
/// short candidates.
pub(super) fn is_repetitive(seq: &[u8]) -> bool {
    (1..=8.min(seq.len() / 2).min(seq.len().saturating_sub(10))).any(|period| {
        let matches = seq
            .iter()
            .zip(&seq[period..])
            .filter(|(a, b)| a == b)
            .count();
        matches * 100 >= (seq.len() - period) * 65
    })
}

/// Counts each exact k-mer once per window, excludes short tandem repeats,
/// sorts by count descending then code ascending, and retains at most `top`.
pub(super) fn top_kmers(windows: &[&[u8]], k: usize, top: usize) -> Vec<(u64, u32)> {
    let mut counts: KmerMap<(u32, usize)> =
        KmerMap::with_capacity_and_hasher(windows.len() * 8, Default::default());
    assert!((1..=32).contains(&k));
    let mask = u64::MAX >> (64 - 2 * k);
    for (window_index, &wnd) in windows.iter().enumerate() {
        let (mut code, mut valid) = (0, 0);
        for &base in wnd {
            if let Some(bits) = encode_kmer(&[base]) {
                code = ((code << 2) | bits) & mask;
                valid += 1;
                if valid >= k {
                    let entry = counts.entry(code).or_insert((0, usize::MAX));
                    if entry.1 != window_index {
                        entry.0 += 1;
                        entry.1 = window_index;
                    }
                }
            } else {
                code = 0;
                valid = 0;
            }
        }
    }
    let mut ranked: Vec<(u64, u32)> = counts
        .into_iter()
        .map(|(code, (count, _))| (code, count))
        .filter(|&(code, _)| {
            if k >= 4 {
                !(1..=8.min(k / 2))
                    .any(|period| code & (mask >> (2 * period)) == code >> (2 * period))
            } else {
                !is_low_complexity(&decode_kmer(code, k))
            }
        })
        .collect();
    // Count descending, then code ascending for a deterministic tie-break.
    let order = |a: &(u64, u32), b: &(u64, u32)| b.1.cmp(&a.1).then(a.0.cmp(&b.0));
    if ranked.len() > top {
        ranked.select_nth_unstable_by(top, order);
        ranked.truncate(top);
    }
    ranked.sort_unstable_by(order);
    ranked
}

/// Reconstructs a consensus adapter from weighted k-mer nodes by a cycle-safe
/// bidirectional greedy walk: seeds at the heaviest node, then extends both ways
/// through the heaviest unvisited neighbor. The visited set keeps this a simple
/// path; a length-bounded DP would re-traverse positive-weight cycles into a
/// long repetitive consensus. The walk is bidirectional because the heaviest
/// seed usually sits mid-adapter, so forward-only extension would recover only
/// the suffix.
///
/// Returns `(consensus, per-position weights, total weight)`, or `None` when
/// `nodes` is empty. `lmax` caps length, but at least one k-mer is always kept.
pub(super) fn bounded_heaviest_path(
    nodes: &[(u64, u32)],
    k: usize,
    lmax: usize,
) -> Option<(Vec<u8>, Vec<u32>, u64)> {
    use std::collections::HashMap;
    if nodes.is_empty() {
        return None;
    }
    let n = nodes.len();
    // Edge A to B exists when the last k-1 bases of A equal the first k-1 bases
    // of B; on 2-bit codes, `(A & suffix_mask) == (B >> 2)`.
    let suffix_mask: u64 = if k >= 1 {
        (1u64 << (2 * (k - 1))) - 1
    } else {
        0
    };
    // Successor index: (k-1)-prefix code to nodes whose prefix equals that code.
    // Predecessor index: (k-1)-suffix code to nodes whose suffix equals that code.
    let mut by_prefix: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut by_suffix: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, &(code, _)) in nodes.iter().enumerate() {
        by_prefix.entry(code >> 2).or_default().push(i);
        by_suffix.entry(code & suffix_mask).or_default().push(i);
    }

    // Greater support first, then the smaller code for deterministic ties.
    let weight_desc_code_asc = |&a: &usize, &b: &usize| {
        nodes[a]
            .1
            .cmp(&nodes[b].1)
            .then(nodes[b].0.cmp(&nodes[a].0))
    };

    // Deterministic pick: heaviest unvisited candidate, ties to the smaller code.
    let pick = |cands: Option<&Vec<usize>>, visited: &[bool]| -> Option<usize> {
        cands?
            .iter()
            .copied()
            .filter(|&i| !visited[i])
            .max_by(weight_desc_code_asc)
    };

    // Seed: the single heaviest node, ties to the smaller code.
    let seed = (0..n).max_by(weight_desc_code_asc).unwrap();
    let mut visited = vec![false; n];
    visited[seed] = true;

    // Forward extension: heaviest unvisited successor until none remains or
    // `lmax` is reached.
    let mut forward: Vec<usize> = Vec::new();
    let mut cur = seed;
    while k + forward.len() < lmax {
        match pick(by_prefix.get(&(nodes[cur].0 & suffix_mask)), &visited) {
            Some(v) => {
                visited[v] = true;
                forward.push(v);
                cur = v;
            },
            None => break,
        }
    }
    // Backward extension: heaviest unvisited predecessor.
    let mut backward: Vec<usize> = Vec::new();
    cur = seed;
    while k + forward.len() + backward.len() < lmax {
        match pick(by_suffix.get(&(nodes[cur].0 >> 2)), &visited) {
            Some(u) => {
                visited[u] = true;
                backward.push(u);
                cur = u;
            },
            None => break,
        }
    }

    // Full path: reverse(backward) ++ [seed] ++ forward.
    let mut chain: Vec<usize> = backward.iter().rev().copied().collect();
    chain.push(seed);
    chain.extend(forward.iter().copied());

    // The consensus: the first node emits k bases, each subsequent node its
    // last base.
    let mut cons = decode_kmer(nodes[chain[0]].0, k);
    let mut profile: Vec<u32> = vec![nodes[chain[0]].1; k];
    let mut weight: u64 = nodes[chain[0]].1 as u64;
    for &idx in &chain[1..] {
        cons.push(*decode_kmer(nodes[idx].0, k).last().unwrap());
        profile.push(nodes[idx].1);
        weight += nodes[idx].1 as u64;
    }
    Some((cons, profile, weight))
}

/// Peels up to `MAX_ADAPTERS_PER_END` distinct adapter consensuses out of one
/// end's weighted k-mer graph: each round runs `bounded_heaviest_path`, then
/// removes the k-mers of that path's outermost support run from `nodes` so
/// the next round is forced onto a different path. Sequence past a layer
/// boundary stays available to later paths. Every path is returned; support
/// floors apply to the validated candidates.
pub(super) fn peel_paths(
    mut nodes: Vec<(u64, u32)>,
    k: usize,
    end: End,
) -> Vec<(Vec<u8>, Vec<u32>)> {
    let mut out = Vec::new();
    while out.len() < MAX_ADAPTERS_PER_END {
        let Some((cons, _, _)) = bounded_heaviest_path(&nodes, k, LMAX) else {
            break;
        };
        let current: KmerMap<u32> = nodes.iter().copied().collect();
        let weights: Vec<u32> = cons
            .windows(k)
            .map(|w| {
                encode_kmer(w)
                    .and_then(|code| current.get(&code).copied())
                    .unwrap_or(0)
            })
            .collect();
        let (lo, hi, _) = supported_span(&weights, end);
        let mut used: std::collections::HashSet<u64> = cons[lo..hi.min(cons.len())]
            .windows(k)
            .filter_map(encode_kmer)
            .collect();
        if used.is_empty() {
            used = cons.windows(k).filter_map(encode_kmer).collect();
        }
        nodes.retain(|(code, _)| !used.contains(code));
        out.push((cons, weights));
        if nodes.is_empty() {
            break;
        }
    }
    out
}
