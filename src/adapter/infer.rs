//! Ab-initio adapter inference from recurrent read-end sequences.
//!
//! Exact k-mer graphs establish consensus boundaries. Batched approximate
//! matching validates support and aligns primer extensions to conserved insert
//! starts. Catalog sequences annotate discoveries without choosing their bases.

use crate::adapter::search::{AmbiguousSearcher, hits, is_plain_acgt, new_ambiguous_searcher};
use crate::adapter::{Adapter, AdapterConfig, MIN_PATTERN_LEN, Role, edit_budget};

/// The physical read end from which a consensus was assembled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum End {
    /// Discovered in the 5' windows.
    Five,
    /// Discovered in the 3' windows.
    Three,
}

/// k-mer length used for end-window counting and assembly graph nodes.
const KMER_K: usize = 16;

/// Number of top exact k-mers retained per end for graph assembly.
const TOP_KMERS: usize = 500;

/// Length of the 5'/3' end window scanned per read for adapter discovery.
const WINDOW_LEN: usize = 100;

/// Minimum presence-fraction support required to keep a discovered adapter.
/// Support is the fraction of sampled end windows containing the consensus
/// within its length-scaled edit budget. The threshold retains common library
/// adapters while excluding sparse barcode-specific sequences and background.
const KEEP_SUPPORT: f64 = 0.15;

/// Maximum windows used for alignment and support validation.
const RECOUNT_WINDOWS: usize = 4000;

/// Minimum k-mer support relative to the path peak at an assembly boundary.
const BOUNDARY_SUPPORT: f64 = 0.45;

/// Exact prefix length used to locate recurrent unprimed read starts.
const START_K: usize = 11;

/// Max total emitted length of a single `bounded_heaviest_path` consensus,
/// used by `peel_paths` so no single peel can run away in length.
const LMAX: usize = 100;

/// Max number of adapters `peel_paths` will extract from one end's k-mer graph.
const MAX_ADAPTERS_PER_END: usize = 3;

/// Minimum fraction of the first (heaviest) path's weight a peeled path needs
/// to be kept; a lighter path is background rather than a distinct adapter.
const MIN_PATH_WEIGHT_FRAC: f64 = 0.25;

/// Minimum percent identity for a catalog entry to be reported as the match
/// of an inferred adapter. A 16 to 32 bp anchor searched against every catalog
/// entry on both strands names something spurious well above the 60 percent
/// that its trimming budget alone would allow.
const NAME_IDENTITY_MIN: f32 = 85.0;

/// One discovered adapter with inference metadata. The bare `Adapter` (without
/// `support` and `name_hits`) is extracted only when building the trim config.
#[derive(Debug, Clone)]
pub struct InferredAdapter {
    /// Sequence used for trimming (or printed as the recommendation), named
    /// `inferred_N` by presentation order.
    pub adapter: Adapter,
    /// Complete sequence retained after boundary validation.
    pub assembled_seq: Vec<u8>,
    /// Fraction of sampled end windows containing the consensus within its
    /// edit budget.
    pub support: f64,
    /// Catalog entries within `NAME_IDENTITY_MIN` of the consensus, best
    /// first, as `(name, percent identity)`. An annotation, not the name.
    pub name_hits: Vec<(String, f32)>,
}

impl InferredAdapter {
    /// Returns the number of assembled bases excluded from the trimming sequence.
    pub fn uncertain_bases(&self) -> usize {
        self.assembled_seq
            .len()
            .saturating_sub(self.adapter.seq.len())
    }
}

/// Encodes a k-mer at 2 bits per base (A=0, C=1, G=2, T=3). `None` when it
/// holds any byte other than uppercase ACGT (reads are uppercased upstream) or
/// is longer than 32 bases.
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

/// Decodes a 2-bit code back to its `k` bases; the inverse of `encode_kmer`.
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

/// Slices the first and last `w` bytes of each read into 5' and 3' window
/// lists. Returns the 5' windows (`&read[..min(w, len)]`) and the 3' windows
/// (`&read[len - min(w, len)..]`). Empty reads are skipped, and so is any window
/// holding a byte outside ACGT: an uncalled base is evidence of nothing, and on
/// the IUPAC profile it would match every k-mer for free.
fn end_windows<'a>(sample: &[&'a [u8]], w: usize) -> (Vec<&'a [u8]>, Vec<&'a [u8]>) {
    let mut five = Vec::new();
    let mut three = Vec::new();
    for &read in sample {
        let n = read.len();
        if n == 0 {
            continue;
        }
        let take = w.min(n);
        let head = &read[..take];
        let tail = &read[n - take..];
        if is_plain_acgt(head) {
            five.push(head);
        }
        if is_plain_acgt(tail) {
            three.push(tail);
        }
    }
    (five, three)
}

/// Returns whether a k-mer is too low-complexity to serve as an adapter seed: a
/// homopolymer or a dinucleotide repeat.
fn is_low_complexity(kmer: &[u8]) -> bool {
    if kmer.windows(2).all(|w| w[0] == w[1]) {
        return true; // homopolymer
    }
    // Period-2 repeat, such as ACACAC.
    if kmer.len() >= 4 && kmer.iter().enumerate().all(|(i, &b)| b == kmer[i % 2]) {
        return true;
    }
    false
}

/// Returns the exact k-mer counts across all windows, low-complexity k-mers
/// dropped, sorted by count descending then code ascending, and truncated to
/// `top`.
fn top_kmers(windows: &[&[u8]], k: usize, top: usize) -> Vec<(u64, u32)> {
    use std::collections::HashMap;
    let mut counts: HashMap<u64, u32> = HashMap::new();
    assert!((1..=32).contains(&k));
    let mask = u64::MAX >> (64 - 2 * k);
    for &wnd in windows {
        let (mut code, mut valid) = (0, 0);
        for &base in wnd {
            if let Some(bits) = encode_kmer(&[base]) {
                code = ((code << 2) | bits) & mask;
                valid += 1;
                if valid >= k {
                    *counts.entry(code).or_insert(0) += 1;
                }
            } else {
                code = 0;
                valid = 0;
            }
        }
    }
    let mut ranked: Vec<(u64, u32)> = counts
        .into_iter()
        .filter(|&(code, _)| {
            if k >= 4 {
                code & (mask >> 4) != code >> 4
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

/// Counts the distinct `windows` with at least one forward approximate
/// occurrence of `pattern` within `max_edits`. Each window counts at most once,
/// however often `pattern` occurs in it. `searcher` must be forward-only (see
/// `new_searcher_fwd`) so reverse-complement occurrences do not inflate the
/// count. Callers provide an already bounded window sample.
fn windows_containing(
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
/// removes that path's k-mers from `nodes` so the next round is forced onto a
/// different, non-overlapping path. Stops early once a path's weight falls
/// below `MIN_PATH_WEIGHT_FRAC` of the first (heaviest) path's weight, or once
/// no path or no nodes remain.
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
        // The nodes this path used are removed so the next peel finds a
        // different one.
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

/// Returns whether `a` and `b` are the same adapter within `error_rate`: an
/// approximate occurrence of the shorter in the longer on either strand (the
/// both-strand searcher covers the reverse-complement case).
fn same_adapter(a: &[u8], b: &[u8], error_rate: f64) -> bool {
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
fn name_against(seq: &[u8], refs: &[Adapter], error_rate: f64) -> Vec<(String, f32)> {
    let mut s = new_ambiguous_searcher();
    let mut named: Vec<(String, f32)> = Vec::new();
    for r in refs {
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
fn stride_sample<'a>(windows: &[&'a [u8]], cap: usize) -> Vec<&'a [u8]> {
    let step = windows.len().div_ceil(cap.max(1)).max(1);
    windows.iter().step_by(step).copied().collect()
}

/// Encodes a set of concrete DNA bases as an IUPAC symbol.
fn ambiguity_code(mask: usize) -> u8 {
    b"-ACMGRSVTWYHKDBN"[mask]
}

/// Refines an assembled path with one best alignment per supporting window.
/// Majority substitutions and deletions correct graph branches caused by
/// sequencing errors without extending the assembly into unaligned sequence.
fn polish_consensus(seq: &[u8], windows: &[&[u8]]) -> Vec<u8> {
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let mut best: Vec<Option<sassy::Match>> = vec![None; windows.len()];
    for hit in searcher.search_texts(seq, windows, edit_budget(0.25, seq.len())) {
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
    if aligned < 20 {
        return seq.to_vec();
    }
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

/// Extends a conserved insert anchor toward the physical read end. Each round
/// aligns the current consensus before voting on the preceding base, allowing
/// indels to shift the supporting reads without shifting the consensus.
fn upstream_consensus(anchor: &[u8], windows: &[&[u8]]) -> Vec<u8> {
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
        if total < 20.max(windows.len() / 5) {
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
    }
    seq.truncate(seq.len() - anchor.len());
    seq
}

/// Uses recurrent unprimed read starts to locate a conserved insert boundary.
/// A primer must be independently supported upstream of that boundary; the
/// conserved insert itself is excluded from the inferred trimming sequence.
fn contrast_boundary(consensus: &[u8], windows: &[&[u8]], end: End) -> Option<Vec<Vec<u8>>> {
    let reverse = end == End::Three;
    let oriented: Vec<Vec<u8>> = windows
        .iter()
        .map(|w| {
            if reverse {
                w.iter().rev().copied().collect()
            } else {
                w.to_vec()
            }
        })
        .collect();
    let windows: Vec<&[u8]> = oriented.iter().map(Vec::as_slice).collect();
    let mut seq: Vec<u8> = if reverse {
        consensus.iter().rev().copied().collect()
    } else {
        consensus.to_vec()
    };
    let mut starts = std::collections::HashMap::<u64, usize>::new();
    for window in &windows {
        if window.len() >= START_K
            && let Some(code) = encode_kmer(&window[..START_K])
        {
            *starts.entry(code).or_default() += 1;
        }
    }
    let floor = 10.max(windows.len() / 32);
    let positions_in = |seq: &[u8]| -> Vec<(usize, usize)> {
        seq.windows(START_K)
            .enumerate()
            .take(seq.len().saturating_sub(KMER_K) + 1)
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
        let mut prefix = upstream_consensus(&seq[..KMER_K], &windows);
        if !prefix.is_empty() {
            prefix.extend_from_slice(&seq);
            seq = prefix;
            positions = positions_in(&seq);
        }
    }
    positions.sort_by_key(|&(pos, count)| (std::cmp::Reverse(count), pos));
    for (pos, _) in positions {
        let anchor = &seq[pos..pos + KMER_K];
        let primer = upstream_consensus(anchor, &windows);
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
            for hit in searcher.search_texts(anchor, &windows, 2) {
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
                .filter(|(_, group)| group.len() >= 20.max(windows.len() / 10))
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
            return Some(primers);
        }
    }
    None
}

/// Assembles one end's candidates using exact k-mer support and validates each
/// complete sequence against a bounded, uniformly spaced window sample.
fn assemble(windows: &[&[u8]], base: &AdapterConfig, end: End) -> Vec<(Vec<u8>, f64)> {
    if windows.len() < 3 {
        return Vec::new();
    }
    // K-mer encoding and approximate matching operate on uppercase DNA.
    // Inference owns normalized copies and does not modify pipeline records.
    let upper: Vec<Vec<u8>> = windows.iter().map(|w| w.to_ascii_uppercase()).collect();
    let windows: Vec<&[u8]> = upper.iter().map(Vec::as_slice).collect();
    let windows = windows.as_slice();

    let exact = top_kmers(windows, KMER_K, TOP_KMERS);
    if exact.is_empty() {
        return Vec::new();
    }
    // Validation uses windows distributed across the complete sample.
    let recount = stride_sample(windows, RECOUNT_WINDOWS);
    let n_recount = recount.len();

    let mut fwd = crate::adapter::search::new_searcher_fwd();
    let weighted = exact;
    let mut out = Vec::new();
    for (cons, profile) in peel_paths(weighted, KMER_K) {
        let peak = *profile.iter().max().unwrap_or(&0);
        let weights = &profile[KMER_K - 1..];
        let floor = (peak as f64 * BOUNDARY_SUPPORT) as u32;
        let lo = weights.iter().position(|&w| w >= floor).unwrap_or(0);
        let hi = weights.iter().rposition(|&w| w >= floor).unwrap_or(0) + KMER_K;
        let trimmed = cons[lo..hi].to_vec();
        let trimmed = if peak as usize * 4 < windows.len() {
            polish_consensus(&trimmed, &recount)
        } else {
            trimmed
        };
        let sequences = if let Some(primers) = contrast_boundary(&cons, &recount, end) {
            primers
        } else {
            let bounded = match end {
                End::Five if hi < cons.len() => {
                    weights[hi - KMER_K + 1] * 2 <= weights[hi - KMER_K]
                },
                End::Three if lo > 0 => weights[lo - 1] * 2 <= weights[lo],
                _ => cons.len() < LMAX,
            };
            if !bounded {
                tracing::debug!(sequence = %String::from_utf8_lossy(&trimmed),
                    "Recurrent sequence has no supported insert boundary");
                continue;
            }
            vec![trimmed]
        };
        for trimmed in sequences {
            if trimmed.len() < MIN_PATTERN_LEN {
                continue;
            }
            // Presence counts each supporting window once, including reads
            // whose sequencing errors disrupted individual exact k-mers.
            let k_cons = edit_budget(base.error_rate, trimmed.len());
            let present = windows_containing(&mut fwd, &trimmed, &recount, k_cons);
            let support = present as f64 / n_recount as f64;
            out.push((trimmed, support));
        }
    }
    out
}

/// Discovers supported adapter sequences independently of the reference catalog.
/// Equivalent end assemblies share one trimming pattern. Catalog and supplied
/// FASTA entries provide names only after the inferred boundaries are fixed.
pub fn discover(sample: &[&[u8]], base: &AdapterConfig) -> Vec<InferredAdapter> {
    let (five_w, three_w) = end_windows(sample, WINDOW_LEN);
    let five = assemble(&five_w, base, End::Five);
    let three = assemble(&three_w, base, End::Three);

    let refs = crate::adapter::preset::preset(crate::adapter::preset::Kit::ALL);
    let name_refs: Vec<Adapter> = refs
        .into_iter()
        .chain(base.adapters.iter().cloned())
        .collect();

    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let mut candidates: Vec<(Vec<u8>, f64, u32)> = five
        .into_iter()
        .chain(three)
        .filter(|(seq, support)| seq.len() >= MIN_PATTERN_LEN && *support >= KEEP_SUPPORT)
        .map(|(seq, support)| {
            let exact = windows_containing(&mut searcher, &seq, &five_w, 0)
                + windows_containing(&mut searcher, &seq, &three_w, 0);
            (seq, support, exact)
        })
        .collect();
    // Exact support selects the reconstruction before fuzzy duplicates merge;
    // approximate support alone cannot distinguish a correct consensus from
    // several nearby error variants.
    candidates.sort_by(|a, b| {
        b.2.cmp(&a.2)
            .then(b.1.total_cmp(&a.1))
            .then(b.0.len().cmp(&a.0.len()))
            .then(a.0.cmp(&b.0))
    });
    let mut distinct: Vec<(Vec<u8>, f64)> = Vec::new();
    for (seq, support, _) in candidates {
        if let Some((_, previous)) = distinct
            .iter_mut()
            .find(|(other, _)| same_adapter(&seq, other, base.error_rate))
        {
            *previous = previous.max(support);
        } else {
            distinct.push((seq, support));
        }
    }
    distinct.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    distinct
        .into_iter()
        .enumerate()
        .map(|(i, (seq, support))| {
            let name_hits = name_against(&seq, &name_refs, base.error_rate);
            InferredAdapter {
                adapter: Adapter {
                    name: format!("inferred_{}", i + 1),
                    seq: seq.clone(),
                    role: Role::Adapter,
                },
                assembled_seq: seq,
                support,
                name_hits,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_bases(mut state: u64, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| {
                state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                z ^= z >> 31;
                b"ACGT"[(z >> 62) as usize]
            })
            .collect()
    }

    fn infer_owned(reads: &[Vec<u8>]) -> Vec<InferredAdapter> {
        let sample: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
        discover(
            &sample,
            &AdapterConfig {
                adapters: vec![],
                error_rate: 0.2,
                end_size: 150,
                split: true,
                min_piece: 20,
                candidate_index: std::sync::OnceLock::new(),
            },
        )
    }

    #[test]
    fn recovers_three_distinct_adapters_without_length_cap() {
        let adapters: Vec<Vec<u8>> = [27, 43, 71]
            .into_iter()
            .enumerate()
            .map(|(i, len)| random_bases(313 + i as u64, len))
            .collect();
        let reads: Vec<Vec<u8>> = (0..900)
            .map(|i| {
                let mut read = adapters[i % 3].clone();
                read.extend(random_bases(4321 + i as u64, 250));
                read
            })
            .collect();
        let found = infer_owned(&reads);
        for seq in &adapters {
            assert!(
                found.iter().any(|d| d.adapter.seq == *seq),
                "missing {}: {found:?}",
                String::from_utf8_lossy(seq)
            );
        }
        assert_eq!(found.len(), 3);
    }

    #[test]
    fn clean_conserved_inserts_are_not_adapters() {
        let prefix = random_bases(7877, 140);
        let suffix = random_bases(4512, 140);
        let reads: Vec<Vec<u8>> = (0..400)
            .map(|i| {
                let mut read = prefix.clone();
                read.extend(random_bases(915 + i, 200));
                read.extend_from_slice(&suffix);
                read
            })
            .collect();
        assert!(infer_owned(&reads).is_empty());
    }

    #[test]
    fn recovers_degenerate_primers_at_conserved_insert_boundaries() {
        let front = b"TCGATGARYCTACGTGACCT";
        let rear = b"GCTAGTACCGATGCTAGTCA";
        let prefix = random_bases(712, 140);
        let suffix = random_bases(815, 140);
        let reads: Vec<Vec<u8>> = (0..500usize)
            .map(|i| {
                let mut read = Vec::new();
                if i % 5 != 0 {
                    read.extend_from_slice(front);
                    read[7] = b"AG"[(i / 5) % 2];
                    read[8] = b"CT"[(i / 10) % 2];
                }
                read.extend_from_slice(&prefix);
                read.extend(random_bases(9712 + i as u64, 200));
                read.extend_from_slice(&suffix);
                if i % 5 != 0 {
                    read.extend_from_slice(rear);
                }
                read
            })
            .collect();
        let found = infer_owned(&reads);
        assert!(found.iter().any(|d| d.adapter.seq == front), "{found:?}");
        assert!(found.iter().any(|d| d.adapter.seq == rear), "{found:?}");
        assert_eq!(found.len(), 2);
    }

    #[test]
    fn ambiguity_codes_cover_exactly_the_voted_bases() {
        for mask in 1..16 {
            let bases = crate::adapter::search::iupac_bases(ambiguity_code(mask)).unwrap();
            let actual = bases
                .iter()
                .fold(0usize, |m, &b| m | (1 << encode_kmer(&[b]).unwrap()));
            assert_eq!(actual, mask);
        }
    }

    #[test]
    fn distinct_primers_at_one_insert_boundary_remain_distinct() {
        let primers = [random_bases(7651, 22), random_bases(2157, 24)];
        let insert = random_bases(2871, 140);
        let reads: Vec<Vec<u8>> = (0..600)
            .map(|i| {
                let mut read = Vec::new();
                if i % 5 != 0 {
                    read.extend_from_slice(&primers[i % 2]);
                }
                read.extend_from_slice(&insert);
                read.extend(random_bases(159 + i as u64, 200));
                read
            })
            .collect();
        let found = infer_owned(&reads);
        for primer in &primers {
            assert!(found.iter().any(|d| d.adapter.seq == *primer), "{found:?}");
        }
        assert_eq!(found.len(), 2);
    }

    #[test]
    fn truncated_insert_anchor_recovers_distinct_upstream_primers() {
        let primers = [random_bases(7651, 22), random_bases(2157, 24)];
        let insert = random_bases(2871, 140);
        let reads: Vec<Vec<u8>> = (0..600)
            .map(|i| {
                let mut read = Vec::new();
                if i % 5 != 0 {
                    read.extend_from_slice(&primers[i % 2]);
                }
                read.extend_from_slice(&insert);
                read.truncate(WINDOW_LEN);
                read
            })
            .collect();
        let windows: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
        let found = contrast_boundary(&insert[1..90], &windows, End::Five).unwrap();
        assert_eq!(found.len(), 2);
        for primer in &primers {
            assert!(found.contains(primer), "{found:?}");
        }
    }

    /// Encoding then decoding a k-mer is the identity.
    #[test]
    fn kmer_codec_roundtrips() {
        let k = b"ACGTACGTACGTACGT"; // 16bp
        let code = encode_kmer(k).unwrap();
        assert_eq!(decode_kmer(code, 16), k);
    }

    /// `encode_kmer` rejects ambiguity codes, lowercase bases and over-long
    /// k-mers.
    #[test]
    fn encode_rejects_non_acgt() {
        assert_eq!(encode_kmer(b"ACGTN"), None);
        assert_eq!(encode_kmer(b"acgt"), None); // lowercase not accepted
        assert_eq!(encode_kmer(&[b'A'; 33]), None); // > 32 bases rejected
    }

    /// A short read yields itself as both windows and an empty read yields
    /// nothing.
    #[test]
    fn end_windows_slices_both_ends() {
        let r1: &[u8] = b"AAAACCCCGGGGTTTTACGTACGT"; // 24bp
        let r2: &[u8] = b"TTTT"; // 4bp, below w: the whole read at both ends
        let sample: Vec<&[u8]> = vec![r1, r2, b""]; // empty skipped
        let (five, three) = end_windows(&sample, 8);
        assert_eq!(five, vec![&r1[..8], r2]); // first 8, then the whole short read
        assert_eq!(three, vec![&r1[16..], r2]); // last 8, then the whole short read
    }

    /// A window holding an ambiguity code is dropped from that end only.
    #[test]
    fn end_windows_drop_windows_holding_ambiguity_codes() {
        let r1: &[u8] = b"AAAANCCCGGGGTTTTACGTACGT"; // `N` in the 5' window only
        let r2: &[u8] = b"acgtacgtacgtacgtacgtacgn"; // `n` in the 3' window only
        let sample: Vec<&[u8]> = vec![r1, r2];
        let (five, three) = end_windows(&sample, 8);
        assert_eq!(five, vec![&r2[..8]]);
        assert_eq!(three, vec![&r1[16..]]);
    }

    /// A 16-mer planted in every window ranks first over unique filler.
    #[test]
    fn top_kmers_ranks_planted_over_background() {
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

    /// A homopolymer window contributes no k-mer.
    #[test]
    fn top_kmers_drops_homopolymer() {
        let windows: Vec<&[u8]> = vec![b"AAAAAAAAAAAAAAAA"]; // pure homopolymer, 16bp
        assert!(
            top_kmers(&windows, 16, 500).is_empty(),
            "Low-complexity k-mer dropped"
        );
    }

    /// Each window counts once however often the k-mer occurs, and a
    /// reverse-complement occurrence does not count.
    #[test]
    fn windows_containing_counts_windows_once_and_ignores_rc() {
        use crate::adapter::search::new_searcher_fwd;
        // The k-mer is not its own reverse complement, so the RC case is
        // meaningful: revcomp(AAAACCCCGGGGTATG) = CATACCCCGGGGTTTT.
        let kmer = b"AAAACCCCGGGGTATG"; // 16bp
        let w0v = b"TTAAAACCCCGGGGTATGTT".to_vec(); // exact occurrence
        let w1v = b"TTAAAACACCGGGGTATGTT".to_vec(); // 1 substitution (C to A)
        let mut w2v = b"AAAACCCCGGGGTATG".to_vec(); // k-mer twice; counted once
        w2v.extend_from_slice(b"GGGGAAAACCCCGGGGTATG");
        let w3v = b"TTCATACCCCGGGGTTTTTT".to_vec(); // reverse-complement only
        let windows: Vec<&[u8]> = vec![&w0v, &w1v, &w2v, &w3v];
        let mut s = new_searcher_fwd();
        // w0 (exact), w1 (1 edit) and w2 (twice, counted once) give 3; w3 (RC
        // only) is excluded.
        assert_eq!(windows_containing(&mut s, kmer, &windows, 2), 3);
    }

    /// Overlapping 4-mers that tile ACGTACG with descending weights along the
    /// intended path reconstruct it: ACGT(9), CGTA(8), GTAC(7), TACG(6).
    #[test]
    fn bounded_heaviest_path_reconstructs_known_consensus() {
        let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
        let nodes = vec![
            mk(b"ACGT", 9),
            mk(b"CGTA", 8),
            mk(b"GTAC", 7),
            mk(b"TACG", 6),
        ];
        let (cons, profile, weight) = bounded_heaviest_path(&nodes, 4, 100).unwrap();
        assert_eq!(cons, b"ACGTACG"); // ACGT + C + A + G: 4 nodes give 7 nt
        assert_eq!(profile.len(), cons.len());
        assert_eq!(weight, 9 + 8 + 7 + 6);
    }

    /// On the 2-node cycle ATAT, TATA, ATAT (k = 4) the visited set stops the
    /// walk after each node is used once, so the consensus is a short simple
    /// path rather than a repeat filling `lmax`.
    #[test]
    fn bounded_heaviest_path_terminates_on_cycle() {
        let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
        let nodes = vec![mk(b"ATAT", 5), mk(b"TATA", 5)];
        let (cons, _profile, _w) = bounded_heaviest_path(&nodes, 4, 12).unwrap();
        assert!(cons.len() <= 12, "No loop: each node used at most once");
        assert!(cons.starts_with(b"ATAT") || cons.starts_with(b"TATA"));
    }

    /// Two non-overlapping tilings with different bases peel as two adapters.
    #[test]
    fn peel_extracts_two_distinct_adapters() {
        let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
        let nodes = vec![
            // Adapter 1: ACGTACG..., high weight.
            mk(b"ACGT", 100),
            mk(b"CGTA", 99),
            mk(b"GTAC", 98),
            // Adapter 2: TTGGTTG..., lower weight but above 25% of 297.
            mk(b"TTGG", 90),
            mk(b"TGGT", 89),
            mk(b"GGTT", 88),
        ];
        let paths = peel_paths(nodes, 4);
        assert_eq!(paths.len(), 2);
        assert!(paths[0].0.starts_with(b"ACGT"));
        assert!(paths[1].0.starts_with(b"TTGG"));
    }

    /// An exact catalog sequence is named at 100 percent identity.
    #[test]
    fn name_against_matches_catalog_entry() {
        let refs = vec![Adapter {
            name: "SQK-TEST".into(),
            seq: b"ACGTACGTACGTACGT".to_vec(),
            role: Role::Adapter,
        }];
        let hits = name_against(b"ACGTACGTACGTACGT", &refs, 0.2);
        assert_eq!(hits[0].0, "SQK-TEST");
        assert!((hits[0].1 - 100.0).abs() < 1e-3);
    }

    /// On a 20 bp reference with a budget of floor(0.2 * 20) = 4 edits, two
    /// substitutions (90 percent) name it, three (85 percent) still do, and
    /// four (80 percent) do not.
    #[test]
    fn name_against_requires_high_identity() {
        let reference = b"GGGGTTTTGGGGTTTTGGGG";
        let refs = vec![Adapter {
            name: "REF".into(),
            seq: reference.to_vec(),
            role: Role::Adapter,
        }];
        let mutate = |count: usize| -> Vec<u8> {
            let mut seq = reference.to_vec();
            for i in 0..count {
                seq[3 + 5 * i] = b'C';
            }
            seq
        };
        assert_eq!(name_against(&mutate(2), &refs, 0.2).len(), 1);
        assert_eq!(name_against(&mutate(3), &refs, 0.2).len(), 1);
        assert!(
            name_against(&mutate(4), &refs, 0.2).is_empty(),
            "80 percent identity is below the naming floor"
        );
    }

    /// Reads starting with an exact catalog adapter (`LSK109_front`) yield an
    /// adapter named `inferred_1` that carries the catalog match separately.
    #[test]
    fn discovered_adapters_are_named_by_order_with_catalog_annotation() {
        let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT";
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..300usize {
            let mut read = adapter.to_vec();
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            for _ in 0..120usize {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            owned.push(read);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(!found.is_empty(), "The planted adapter is discovered");
        for (i, d) in found.iter().enumerate() {
            assert_eq!(d.adapter.name, format!("inferred_{}", i + 1));
        }
        assert_eq!(
            found[0].name_hits.first().map(|(name, _)| name.as_str()),
            Some("LSK109_front"),
            "The catalog match is an annotation: {:?}",
            found[0].name_hits
        );
    }

    /// Sixty `N`s then random bases: windows holding the run are dropped, so
    /// no poly-A-leading consensus is assembled from them.
    #[test]
    fn discover_finds_nothing_in_ambiguity_runs() {
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..200usize {
            let mut read = vec![b'N'; 60];
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            for _ in 0..100usize {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            owned.push(read);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(
            found.is_empty(),
            "An N run is not adapter evidence (got {found:?})"
        );
    }

    /// A catalog-like adapter planted at the 5' end of 500 synthetic reads with
    /// about 10 percent substitution error is recovered within a small edit
    /// distance. The noise is a fixed permutation of error positions per read
    /// index, with no RNG.
    #[test]
    fn discover_recovers_planted_adapter_under_error() {
        let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp, SQK-NSK007-like
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..500usize {
            let mut read = adapter.to_vec();
            // Deterministic genomic tail from a splitmix64-style mix. A formula
            // linear in the position modulo 4 collapses to a phase-rotated ACGT
            // tandem repeat, which is a spurious signal in 100% of reads that
            // crowds out the planted adapter.
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            for _ in 0..120usize {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            // Deterministic substitutions at roughly 10% of adapter positions.
            for p in (0..adapter.len()).step_by(10) {
                let q = (p + i) % adapter.len();
                read[q] = b"ACGT"[(read[q] as usize + 1) % 4];
            }
            owned.push(read);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(!found.is_empty(), "At least one adapter discovered");
        // The top candidate is a 5' or both-end adapter close to the planted
        // sequence.
        let top = &found[0];
        assert!(top.adapter.seq.len() >= MIN_PATTERN_LEN);
        // Near-match to the planted adapter; recovery is approximate.
        let mut s = new_ambiguous_searcher();
        let k = (0.25 * adapter.len() as f64).ceil() as usize;
        assert!(
            !hits(&mut s, &top.adapter.seq, adapter, k).is_empty()
                || !hits(&mut s, adapter, &top.adapter.seq, k).is_empty(),
            "Recovered adapter is within 25% edit distance of the planted one"
        );
    }

    /// `adapter` is planted at the 5' end with heavy substitutions (weak
    /// recovery) and its exact reverse complement at the 3' end (strong
    /// recovery) of every read, so `merge_both_ends` folds the two per-end
    /// discoveries into a single `End::Both` entry (per `same_adapter`). The
    /// noisy 5' and exact 3' assemblies are fuzzy-equivalent, so the merged
    /// adapter inherits the stronger 3' support.
    #[test]
    fn discover_dual_end_adapter_gets_max_support() {
        let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp
        let rc: Vec<u8> = adapter
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!("Adapter is pure ACGT"),
            })
            .collect();
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..200usize {
            // 5' copy: deterministic substitutions at every 6th (shifted)
            // position; weak but still independently recoverable.
            let mut read = adapter.to_vec();
            for p in (0..adapter.len()).step_by(6) {
                let q = (p + i) % adapter.len();
                read[q] = b"ACGT"[(read[q] as usize + 1) % 4];
            }
            // Deterministic non-periodic genomic middle, from the same
            // splitmix64 mix as the other `discover_*` fixtures.
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            for _ in 0..150usize {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            // 3' copy: exact reverse complement with no error, giving strong
            // recovery.
            read.extend_from_slice(&rc);
            owned.push(read);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);

        // The merged entry is the one near the planted adapter; recovery is
        // approximate, so the match is within 25% edit distance.
        let mut s = new_ambiguous_searcher();
        let k = (0.25 * adapter.len() as f64).ceil() as usize;
        let near: Vec<&InferredAdapter> = found
            .iter()
            .filter(|d| {
                !hits(&mut s, &d.adapter.seq, adapter, k).is_empty()
                    || !hits(&mut s, adapter, &d.adapter.seq, k).is_empty()
            })
            .collect();
        assert_eq!(
            near.len(),
            1,
            "The shared 5'/3' adapter is discovered as a single entry: {found:?}"
        );

        // The reported support reflects the stronger 3' end, not the weaker 5'
        // end alone (about 0.18). The unmerged `Five` entries this fixture also
        // produces carry that value and are dropped independently because
        // 0.18 < `KEEP_SUPPORT`.
        assert!(
            near[0].support > 0.7,
            "Merged adapter's support ({}) must reflect the max across ends \
             (the 3' end recovers at about 1.0), not the weaker 5' end alone (about 0.18)",
            near[0].support
        );
    }

    /// Deterministic non-periodic background (SplitMix64-derived bases) yields
    /// no adapter.
    #[test]
    fn discover_finds_nothing_in_clean_reads() {
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..300usize {
            let mut read = Vec::new();
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            for _ in 0..200usize {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            owned.push(read);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(
            found.is_empty(),
            "No spurious adapter in clean reads (got {found:?})"
        );
    }

    /// With `len <= cap` every window is returned in order (step 1).
    #[test]
    fn stride_sample_is_identity_when_within_cap() {
        let a: &[u8] = b"A";
        let b: &[u8] = b"C";
        let c: &[u8] = b"G";
        let windows: Vec<&[u8]> = vec![a, b, c];
        assert_eq!(stride_sample(&windows, 4), windows);
    }

    /// A four-element sample spans all 13 input positions rather than a prefix.
    #[test]
    fn stride_sample_spans_the_whole_range_not_just_a_prefix() {
        let bytes: Vec<u8> = (0..13u8).map(|i| b'A' + i).collect();
        let windows: Vec<&[u8]> = bytes.iter().map(std::slice::from_ref).collect();
        let sampled = stride_sample(&windows, 4);
        assert!(sampled.len() <= 4);
        // Expected indices: 0, 4, 8, 12.
        assert_eq!(
            sampled,
            vec![windows[0], windows[4], windows[8], windows[12]]
        );
        let last_idx = 12usize; // index of the last sampled window
        assert!(
            last_idx >= (13usize * 2).div_ceil(3),
            "Last sampled window must fall in the last third of the range, not a prefix"
        );
        assert_eq!(*sampled.last().unwrap(), windows[last_idx]);
    }

    /// The adapter occurs only in the latter half of an 8001-read sample, so
    /// the bounded recount has to cover the complete input range. Ignored by
    /// default; `cargo test --lib` runs it only when `--ignored` is passed
    /// through to the test binary.
    #[test]
    #[ignore]
    fn discover_is_not_order_biased_by_recount_window_cap() {
        let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp, same as the other discover_* fixtures
        let n_clean = RECOUNT_WINDOWS + 1;
        let n_planted = RECOUNT_WINDOWS; // 4000

        // Deterministic non-periodic background, from the same splitmix64 mix
        // as the other `discover_*` fixtures.
        let splitmix_tail = |i: usize, len: usize| -> Vec<u8> {
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            let mut out = Vec::with_capacity(len);
            for _ in 0..len {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                out.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            out
        };

        let mut owned: Vec<Vec<u8>> = Vec::with_capacity(n_clean + n_planted);
        for i in 0..n_clean {
            owned.push(splitmix_tail(i, 40)); // pure background, no adapter
        }
        for i in 0..n_planted {
            let mut read = adapter.to_vec();
            read.extend(splitmix_tail(n_clean + i, 12));
            owned.push(read);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(
            !found.is_empty(),
            "Adapter present in a clear majority of reads after the first \
             RECOUNT_WINDOWS must be discovered (got {found:?})"
        );
        let mut s = new_ambiguous_searcher();
        let k = (0.25 * adapter.len() as f64).ceil() as usize;
        assert!(
            found.iter().any(|d| {
                !hits(&mut s, &d.adapter.seq, adapter, k).is_empty()
                    || !hits(&mut s, adapter, &d.adapter.seq, k).is_empty()
            }),
            "Discovered adapters must include one within 25% edit distance \
             of the planted adapter: {found:?}"
        );
    }

    /// Lowercase reads produce the same inferred adapter as uppercase DNA, and
    /// the discovered sequence is uppercase.
    #[test]
    fn discover_recovers_planted_adapter_from_lowercase_reads() {
        let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..500usize {
            let mut read = adapter.to_vec();
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            for _ in 0..120usize {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            let lower: Vec<u8> = read.iter().map(u8::to_ascii_lowercase).collect();
            owned.push(lower);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(
            !found.is_empty(),
            "Lowercase reads must be inferable (got {found:?})"
        );
        let top = &found[0];
        let mut s = new_ambiguous_searcher();
        let k = (0.25 * adapter.len() as f64).ceil() as usize;
        assert!(
            !hits(&mut s, &top.adapter.seq, adapter, k).is_empty()
                || !hits(&mut s, adapter, &top.adapter.seq, k).is_empty(),
            "Discovered adapter (seq {:?}) must be within 25% edit distance \
             of the uppercase planted adapter",
            String::from_utf8_lossy(&top.adapter.seq)
        );
        // The discovered sequence is uppercase ACGT and carries no lowercase
        // byte through from the input.
        assert!(
            top.adapter.seq.iter().all(u8::is_ascii_uppercase),
            "Discovered sequence must be uppercase: {:?}",
            String::from_utf8_lossy(&top.adapter.seq)
        );
    }
}
