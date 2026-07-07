//! Ab-initio adapter inference (Phase 2). Discovers adapters de novo from a
//! read sample using Porechop_ABI's published method: read-end k-mer counting,
//! a weighted de Bruijn graph, length-bounded heaviest-path assembly, iterative
//! peeling, boundary drop-trim, and presence-fraction confidence. Implemented
//! from the paper (not translated from GPL source). Pure and format-neutral.

use crate::adapter::Adapter;
use crate::adapter::search::{DnaSearcher, hits};

/// Cap on the number of windows scanned per k-mer during the 2-error recount
/// (Task 9's confidence pass), bounding its cost on large samples.
const RECOUNT_WINDOWS: usize = 4000;

/// Max total emitted length of a single `bounded_heaviest_path` consensus,
/// used by `peel_paths` so no single peel can run away in length.
const LMAX: usize = 100;

/// Max number of adapters `peel_paths` will extract from one end's k-mer graph.
const MAX_ADAPTERS_PER_END: usize = 3;

/// A peeled path is kept only if its total weight is at least this fraction
/// of the first (heaviest) path's weight; below that it's noise, not a
/// distinct adapter.
const MIN_PATH_WEIGHT_FRAC: f64 = 0.25;

/// Neighbourhood size (in profile positions) `drop_trim` scans inward from
/// each end when looking for a sharp support drop.
const DROP_WINDOW: usize = 7;

/// Fraction of the profile's max weight added to the median-of-diffs baseline
/// to form `drop_trim`'s cut threshold.
const CUT_RATIO: f64 = 0.075;

/// One discovered adapter with inference metadata. Convert to a bare `Adapter`
/// (dropping `support`/`name_hits`) only when building the trim config.
#[derive(Debug, Clone)]
pub struct InferredAdapter {
    pub adapter: Adapter,
    pub support: f64,
    pub name_hits: Vec<(String, f32)>,
}

/// 2-bit-encode a k-mer (A=0,C=1,G=2,T=3). `None` if it contains any non-ACGT
/// base (e.g. `N`) or is longer than 32. Deterministic, case-sensitive to
/// uppercase (reads are uppercased upstream; lowercase/other -> None).
fn encode_kmer(bytes: &[u8]) -> Option<u64> {
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

/// Inverse of `encode_kmer` for a known length `k`.
fn decode_kmer(mut code: u64, k: usize) -> Vec<u8> {
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

/// Slices the first/last `w` bytes of each read into 5' and 3' window lists.
/// Returns `.0` = 5' windows (`&read[..min(w,len)]`), `.1` = 3' windows
/// (`&read[len-min(w,len)..]`). Empty reads are skipped.
// Not yet called from production code (only from the tests below); the
// Task 9 (discover) wires this in and this allow comes off then.
#[allow(dead_code)]
fn end_windows<'a>(sample: &[&'a [u8]], w: usize) -> (Vec<&'a [u8]>, Vec<&'a [u8]>) {
    let mut five = Vec::new();
    let mut three = Vec::new();
    for &read in sample {
        let n = read.len();
        if n == 0 {
            continue;
        }
        let take = w.min(n);
        five.push(&read[..take]);
        three.push(&read[n - take..]);
    }
    (five, three)
}

/// True if a k-mer is too low-complexity to be a useful adapter seed:
/// only 1 distinct base, or a period-1/period-2 repeat (homopolymer or
/// dinucleotide run). Deterministic, no allocation beyond the small set.
fn is_low_complexity(kmer: &[u8]) -> bool {
    if kmer.windows(2).all(|w| w[0] == w[1]) {
        return true; // homopolymer
    }
    // period-2 (e.g. ACACAC...)
    if kmer.len() >= 4 && kmer.iter().enumerate().all(|(i, &b)| b == kmer[i % 2]) {
        return true;
    }
    false
}

/// Exact k-mer counts across all windows, low-complexity k-mers dropped,
/// sorted by `(count desc, code asc)` and truncated to `top`.
// Not yet called from production code (only from the tests below); the
// Task 9 (discover) wires this in and this allow comes off then.
#[allow(dead_code)]
fn top_kmers(windows: &[&[u8]], k: usize, top: usize) -> Vec<(u64, u32)> {
    use std::collections::HashMap;
    let mut counts: HashMap<u64, u32> = HashMap::new();
    for &wnd in windows {
        if wnd.len() < k {
            continue;
        }
        for i in 0..=wnd.len() - k {
            let sub = &wnd[i..i + k];
            if is_low_complexity(sub) {
                continue;
            }
            if let Some(code) = encode_kmer(sub) {
                *counts.entry(code).or_insert(0) += 1;
            }
        }
    }
    let mut ranked: Vec<(u64, u32)> = counts.into_iter().collect();
    // count desc, then code asc for a deterministic tie-break.
    ranked.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    ranked.truncate(top);
    ranked
}

/// Number of distinct `windows` (capped at `RECOUNT_WINDOWS`) with >=1 forward
/// approximate occurrence of `kmer` (edit distance <= `max_edits`). Each window
/// counts at most once, even if `kmer` occurs in it multiple times. `searcher`
/// must be forward-only (see `new_searcher_fwd`) so a window's own reverse-
/// complement can't inflate the count.
// Not yet called from production code (only from the test below); the
// Task 9 (discover) confidence pass wires this in and this allow comes off then.
#[allow(dead_code)]
fn two_error_freq(
    searcher: &mut DnaSearcher,
    kmer: &[u8],
    windows: &[&[u8]],
    max_edits: usize,
) -> u32 {
    let mut present = 0u32;
    for &wnd in windows.iter().take(RECOUNT_WINDOWS) {
        if !hits(searcher, kmer, wnd, max_edits).is_empty() {
            present += 1; // per-window presence, counted once
        }
    }
    present
}

/// Reconstructs a consensus adapter sequence from weighted k-mer nodes via a
/// cycle-safe bidirectional greedy walk: seed at the single heaviest node,
/// then greedily extend forward (heaviest unvisited successor) and backward
/// (heaviest unvisited predecessor), marking each node visited so no node is
/// used twice. A length-bounded DP would traverse cycles repeatedly (weights
/// are positive, so revisiting a cycle only adds weight) and emit a long
/// repetitive consensus; the visited set makes this cycle-safe *and*
/// cycle-*correct* by construction (a simple path, not a loop). Bidirectional
/// extension is required because the heaviest seed k-mer usually sits in the
/// middle of the adapter -- forward-only would recover only the suffix.
/// Returns `(consensus bytes, per-position weight profile, total node-weight)`,
/// or `None` if `nodes` is empty. `lmax` caps the total emitted length; the
/// consensus is always at least one k-mer even if `lmax < k`.
fn bounded_heaviest_path(
    nodes: &[(u64, u32)],
    k: usize,
    lmax: usize,
) -> Option<(Vec<u8>, Vec<u32>, u64)> {
    use std::collections::HashMap;
    if nodes.is_empty() {
        return None;
    }
    let n = nodes.len();
    // (k-1)-overlap: node A -> node B iff last (k-1) bases of A == first (k-1)
    // bases of B. On 2-bit codes: (A & suffix_mask) == (B >> 2).
    let suffix_mask: u64 = if k >= 1 {
        (1u64 << (2 * (k - 1))) - 1
    } else {
        0
    };
    // successor index: (k-1)-prefix code -> nodes whose PREFIX == that code.
    // predecessor index: (k-1)-suffix code -> nodes whose SUFFIX == that code.
    let mut by_prefix: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut by_suffix: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, &(code, _)) in nodes.iter().enumerate() {
        by_prefix.entry(code >> 2).or_default().push(i);
        by_suffix.entry(code & suffix_mask).or_default().push(i);
    }

    // deterministic pick: heaviest unvisited candidate, tie -> smaller code.
    let pick = |cands: Option<&Vec<usize>>, visited: &[bool]| -> Option<usize> {
        cands?
            .iter()
            .copied()
            .filter(|&i| !visited[i])
            .max_by(|&a, &b| {
                nodes[a]
                    .1
                    .cmp(&nodes[b].1)
                    .then(nodes[b].0.cmp(&nodes[a].0))
            })
    };

    // seed = single heaviest node (tie -> smaller code).
    let seed = (0..n)
        .max_by(|&a, &b| {
            nodes[a]
                .1
                .cmp(&nodes[b].1)
                .then(nodes[b].0.cmp(&nodes[a].0))
        })
        .unwrap();
    let mut visited = vec![false; n];
    visited[seed] = true;

    // forward extension: heaviest unvisited successor, until none or lmax reached.
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
    // backward extension: heaviest unvisited predecessor.
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

    // full path: reverse(backward) ++ [seed] ++ forward
    let mut chain: Vec<usize> = backward.iter().rev().copied().collect();
    chain.push(seed);
    chain.extend(forward.iter().copied());

    // build consensus: first node emits k bases, each subsequent emits its last base.
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

fn median_u32(xs: &[u32]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let mut v = xs.to_vec();
    v.sort_unstable();
    let m = v.len() / 2;
    if v.len() % 2 == 1 {
        v[m] as f64
    } else {
        (v[m - 1] as f64 + v[m] as f64) / 2.0
    }
}

/// Presence-fraction confidence for one consensus: `median(profile) /
/// n_windows`, i.e. what share of the sampled end-windows actually
/// contained a k-mer from this path, clamped to `[0,1]`.
// Not yet called from production code (only from the tests below); Task 9
// (discover) wires this in and this allow comes off then.
#[allow(dead_code)]
fn support_from_profile(profile: &[u32], n_windows: usize) -> f64 {
    if n_windows == 0 {
        return 0.0;
    }
    (median_u32(profile) / n_windows as f64).clamp(0.0, 1.0)
}

fn median_f64(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.0;
    }
    let mut v = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = v.len() / 2;
    if v.len() % 2 == 1 {
        v[m]
    } else {
        (v[m - 1] + v[m]) / 2.0
    }
}

/// Trim low-support flanks. Walk from each end inward; cut at the first position
/// where the support jumps by more than the drop threshold relative to the
/// interior plateau. Threshold = median of |successive differences| + CUT_RATIO
/// * max(profile), evaluated over a DROP_WINDOW-sized neighbourhood.
// Not yet called from production code (only from the tests below); Task 9
// (discover) wires this in and this allow comes off then.
#[allow(dead_code)]
fn drop_trim(consensus: &[u8], profile: &[u32]) -> (Vec<u8>, Vec<u32>) {
    let n = profile.len();
    if n == 0 {
        return (consensus.to_vec(), profile.to_vec());
    }
    let maxp = *profile.iter().max().unwrap() as f64;
    let diffs: Vec<f64> = profile
        .windows(2)
        .map(|w| (w[0] as f64 - w[1] as f64).abs())
        .collect();
    let thresh = median_f64(&diffs) + CUT_RATIO * maxp;

    // left boundary: advance while a sharp drop-in from the edge is seen within
    // the first DROP_WINDOW positions.
    let mut lo = 0usize;
    while lo + 1 < n && lo < DROP_WINDOW {
        if (profile[lo] as f64) < maxp - thresh && (profile[lo + 1] as f64) >= maxp - thresh {
            lo += 1;
            break;
        }
        if (profile[lo] as f64) < maxp - thresh {
            lo += 1;
        } else {
            break;
        }
    }
    // right boundary: symmetric from the tail.
    let mut hi = n;
    while hi > lo + 1 && n - hi < DROP_WINDOW {
        if (profile[hi - 1] as f64) < maxp - thresh {
            hi -= 1;
        } else {
            break;
        }
    }
    if lo >= hi {
        return (consensus.to_vec(), profile.to_vec()); // never trim to nothing
    }
    (consensus[lo..hi].to_vec(), profile[lo..hi].to_vec())
}

/// Iteratively peels up to `MAX_ADAPTERS_PER_END` distinct adapter
/// consensuses out of a single end's weighted k-mer graph: each round runs
/// `bounded_heaviest_path`, then removes that path's k-mers from `nodes` so
/// the next round is forced onto a different (non-overlapping) path. Stops
/// early once a path's weight falls below `MIN_PATH_WEIGHT_FRAC` of the
/// first (heaviest) path's weight, or once no path / no nodes remain.
// Not yet called from production code (only from the tests below); Task 9
// (discover) wires this in and this allow comes off then.
#[allow(dead_code)]
fn peel_paths(mut nodes: Vec<(u64, u32)>, k: usize) -> Vec<(Vec<u8>, Vec<u32>)> {
    let mut out = Vec::new();
    let mut first_weight: Option<u64> = None;
    while out.len() < MAX_ADAPTERS_PER_END {
        let Some((cons, profile, weight)) = bounded_heaviest_path(&nodes, k, LMAX) else {
            break;
        };
        let fw = *first_weight.get_or_insert(weight);
        if (weight as f64) < MIN_PATH_WEIGHT_FRAC * fw as f64 {
            break;
        }
        // remove the nodes used by this path so the next peel finds a different one.
        let used: std::collections::HashSet<u64> =
            cons.windows(k).filter_map(encode_kmer).collect();
        nodes.retain(|(code, _)| !used.contains(code));
        out.push((cons, profile));
        if nodes.is_empty() {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_codec_roundtrips() {
        let k = b"ACGTACGTACGTACGT"; // 16bp
        let code = encode_kmer(k).unwrap();
        assert_eq!(decode_kmer(code, 16), k);
    }

    #[test]
    fn encode_rejects_non_acgt() {
        assert_eq!(encode_kmer(b"ACGTN"), None);
        assert_eq!(encode_kmer(b"acgt"), None); // lowercase not accepted
        assert_eq!(encode_kmer(&[b'A'; 33]), None); // > 32 bases rejected
    }

    #[test]
    fn end_windows_slices_both_ends() {
        let r1: &[u8] = b"AAAACCCCGGGGTTTTACGTACGT"; // 24bp
        let r2: &[u8] = b"TTTT"; // 4bp (< w) -> whole read both ends
        let sample: Vec<&[u8]> = vec![r1, r2, b""]; // empty skipped
        let (five, three) = end_windows(&sample, 8);
        assert_eq!(five, vec![&r1[..8], r2]); // first 8 / whole short read
        assert_eq!(three, vec![&r1[16..], r2]); // last 8 / whole short read
    }

    #[test]
    fn top_kmers_ranks_planted_over_background() {
        // A planted 16-mer appears in every window; each window also has unique
        // filler. The planted k-mer must rank first.
        let planted = b"ACGTACGTACGTACGT"; // 16bp, not low-complexity
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..50u8 {
            let mut wnd = planted.to_vec();
            // Varied filler; first byte cycles B..E (never 'A') so a window's
            // filler can never spell "ACGT" and accidentally reconstruct the
            // planted (period-4) k-mer at the trailing slide offset.
            wnd.extend_from_slice(&[b'B' + (i % 4), b'C', b'G', b'T']);
            owned.push(wnd);
        }
        let windows: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let ranked = top_kmers(&windows, 16, 500);
        assert_eq!(decode_kmer(ranked[0].0, 16), planted);
        assert_eq!(ranked[0].1, 50);
    }

    #[test]
    fn top_kmers_drops_homopolymer() {
        let windows: Vec<&[u8]> = vec![b"AAAAAAAAAAAAAAAA"]; // pure homopolymer, 16bp
        assert!(
            top_kmers(&windows, 16, 500).is_empty(),
            "low-complexity dropped"
        );
    }

    #[test]
    fn two_error_freq_counts_windows_once_and_ignores_rc() {
        use crate::adapter::search::new_searcher_fwd;
        // kmer chosen NOT self-reverse-complementary so the RC case is meaningful:
        // revcomp(AAAACCCCGGGGTATG) = CATACCCCGGGGTTTT (distinct from the kmer).
        let kmer = b"AAAACCCCGGGGTATG"; // 16bp
        let w0v = b"TTAAAACCCCGGGGTATGTT".to_vec(); // exact occurrence
        let w1v = b"TTAAAACACCGGGGTATGTT".to_vec(); // 1 substitution (C->A)
        let mut w2v = b"AAAACCCCGGGGTATG".to_vec(); // kmer twice -> counts ONCE
        w2v.extend_from_slice(b"GGGGAAAACCCCGGGGTATG");
        let w3v = b"TTCATACCCCGGGGTTTTTT".to_vec(); // reverse-complement only
        let windows: Vec<&[u8]> = vec![&w0v, &w1v, &w2v, &w3v];
        let mut s = new_searcher_fwd();
        // w0 (exact) + w1 (1 edit) + w2 (twice -> once) = 3; w3 (RC only) excluded.
        assert_eq!(two_error_freq(&mut s, kmer, &windows, 2), 3);
    }

    #[test]
    fn bounded_heaviest_path_reconstructs_known_consensus() {
        // Overlapping 4-mers that tile ACGTACGT with descending weights on the
        // intended path so the heaviest path is unambiguous.
        // ACGT(9) -> CGTA(8) -> GTAC(7) -> TACG(6) -> ACGT... use k=4.
        let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
        let nodes = vec![
            mk(b"ACGT", 9),
            mk(b"CGTA", 8),
            mk(b"GTAC", 7),
            mk(b"TACG", 6),
        ];
        let (cons, profile, weight) = bounded_heaviest_path(&nodes, 4, 100).unwrap();
        assert_eq!(cons, b"ACGTACG"); // ACGT + C + A + G  (4 nodes -> 4+3 = 7 nt)
        assert_eq!(profile.len(), cons.len());
        assert_eq!(weight, 9 + 8 + 7 + 6);
    }

    #[test]
    fn bounded_heaviest_path_terminates_on_cycle() {
        // A real 2-node cycle: ATAT -> TATA -> ATAT (k=4). The visited set must stop
        // the walk after each node is used once (no looping), so the consensus is a
        // short simple path, not a repeat filling lmax.
        let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
        let nodes = vec![mk(b"ATAT", 5), mk(b"TATA", 5)];
        let (cons, _profile, _w) = bounded_heaviest_path(&nodes, 4, 12).unwrap();
        assert!(cons.len() <= 12, "no loop: each node used at most once");
        assert!(cons.starts_with(b"ATAT") || cons.starts_with(b"TATA"));
    }

    #[test]
    fn peel_extracts_two_distinct_adapters() {
        // Two non-overlapping tilings, each internally chained, different bases.
        let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
        let nodes = vec![
            // adapter 1: ACGTACG... high weight
            mk(b"ACGT", 100),
            mk(b"CGTA", 99),
            mk(b"GTAC", 98),
            // adapter 2: TTGGTTG... lower but > 25% of 297
            mk(b"TTGG", 90),
            mk(b"TGGT", 89),
            mk(b"GGTT", 88),
        ];
        let paths = peel_paths(nodes, 4);
        assert_eq!(paths.len(), 2);
        assert!(paths[0].0.starts_with(b"ACGT"));
        assert!(paths[1].0.starts_with(b"TTGG"));
    }

    #[test]
    fn drop_trim_cuts_low_support_flank() {
        // high plateau then a sharp drop -> trailing low-support positions removed.
        let consensus = b"ACGTACGTACGTAAAA".to_vec(); // last 4 are the flank
        let profile = vec![
            100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 3, 3, 3, 3,
        ];
        let (trimmed, tprof) = drop_trim(&consensus, &profile);
        assert_eq!(trimmed, b"ACGTACGTACGT");
        assert_eq!(tprof.len(), trimmed.len());
    }

    #[test]
    fn support_is_median_over_windows() {
        let profile = vec![40u32, 50, 60];
        assert!((support_from_profile(&profile, 100) - 0.50).abs() < 1e-9);
    }
}
