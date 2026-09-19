//! De novo adapter inference from recurrent read-end sequences.
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
/// within its length-scaled edit budget. Absolute support also limits sparse
/// discoveries in small samples.
const KEEP_SUPPORT: f64 = 0.01;

/// Minimum independently supporting read windows for a retained consensus.
const MIN_SUPPORT_WINDOWS: usize = 20;

/// Maximum windows used for alignment and support validation.
const RECOUNT_WINDOWS: usize = 4000;

/// Minimum k-mer support relative to the path peak at an assembly boundary.
const BOUNDARY_SUPPORT: f64 = 0.45;

/// Exact prefix length used to locate recurrent unprimed read starts.
const START_K: usize = 11;

/// Max total emitted length of a single `bounded_heaviest_path` consensus,
/// used by `peel_paths` so no single peel can run away in length.
const LMAX: usize = 100;

/// Maximum distinct insert-boundary anchors evaluated for one consensus.
/// Anchors are ranked by unprimed read-start support.
const MAX_BOUNDARY_ANCHORS: usize = 8;

/// Minimum fraction of a layer's windows that one member of a variable
/// layer must hold, below which a barcode panel would exceed 1000 members.
const VARIABLE_MEMBER_SUPPORT: f64 = 0.001;

/// Max number of adapters `peel_paths` will extract from one end's k-mer graph.
/// A barcode layer holds one family per barcode.
const MAX_ADAPTERS_PER_END: usize = 128;

/// Maximum discovery layers at one read end.
const MAX_LAYERS: usize = 8;

/// Bases from the physical read end searched for a candidate and its mirror.
const MIRROR_WINDOW: usize = 300;

/// Window tested for a mirror at the opposite end.
const MIRROR_SEGMENT: usize = KMER_K;

/// Depth difference between an inner segment and its mirror accepted as
/// symmetric, in bases.
const SYMMETRY_TOLERANCE: usize = 35;

/// Fraction of the terminal error rate allowed when locating a mirror, so
/// that a mirror does not extend past the sequence present at the read end.
const MIRROR_ERROR_RATE: f64 = 0.5;

/// Bases by which the tested window slides per mirror search.
const MIRROR_STEP: usize = 4;

/// Divisor of a candidate's support giving the fewest mirror occurrences
/// that count. A mirror occurs only in reads of the other orientation and is
/// searched at the stricter mirror error rate.
const MIRROR_SUPPORT_DIVISOR: usize = 8;

/// Percentage of a candidate's supporting windows in which it must start
/// within `ANCHOR_SLACK` of the boundary. A majority suffices once earlier
/// layers explain only part of the reads.
const ANCHORED_PERCENT: usize = 60;

/// Bases from the current boundary within which a hit counts as anchored at
/// it and advances it. The slack covers eroded or variable-length remnants
/// of the preceding layer.
const ANCHOR_SLACK: usize = 50;

/// Largest median distance from the physical read end at which an outermost
/// discovered sequence is a sequencing adapter, which splits reads at
/// interior hits; deeper sequences are end-only.
const ADAPTER_FLUSH: usize = 20;

/// Shortest end shared by candidates of one layer that is treated as a
/// neighbouring layer rather than part of the candidate.
const SHARED_END_MIN: usize = 12;

/// Percentage of the shorter of two candidates that their longest common
/// substring must cover for them to count as one family.
const FAMILY_OVERLAP_PERCENT: usize = 60;

/// Percentage of a candidate's supporting windows that a stronger candidate
/// of the same layer must share for the candidate to count as a sequencing
/// variant of it rather than a distinct family.
const VARIANT_OVERLAP_PERCENT: usize = 80;

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
    /// Discovery layer, counted from the read end after any known sequences.
    pub layer: usize,
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
/// lists. Empty reads are skipped, and so is any window holding a byte
/// outside ACGT: an uncalled base is evidence of nothing, and on the IUPAC
/// profile it would match every k-mer for free.
#[cfg(test)]
fn end_windows<'a>(sample: &[&'a [u8]], w: usize) -> (Vec<&'a [u8]>, Vec<&'a [u8]>) {
    let (five, three) = layer_windows(sample, &Boundaries::new(sample), w);
    (
        five.into_iter().map(|(_, window)| window).collect(),
        three.into_iter().map(|(_, window)| window).collect(),
    )
}

/// Hashes packed k-mer codes by multiplicative mixing with a fold of the
/// high bits, which distributes the uniformly encoded keys across buckets.
#[derive(Default, Clone, Copy)]
struct KmerHasher(u64);

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
type KmerMap<V> = std::collections::HashMap<u64, V, std::hash::BuildHasherDefault<KmerHasher>>;

/// Returns whether a k-mer consists of a short tandem repeat.
fn is_low_complexity(kmer: &[u8]) -> bool {
    (1..=8.min(kmer.len() / 2))
        .any(|period| kmer.iter().enumerate().all(|(i, &b)| b == kmer[i % period]))
}

/// Rejects consensuses dominated by a short approximate tandem repeat. Each
/// base is compared with the base one period earlier, so isolated insertions
/// and deletions shift the phase without hiding the repeat. Each period is
/// tested over at least ten comparisons so the threshold stays selective for
/// short candidates.
fn is_repetitive(seq: &[u8]) -> bool {
    (1..=8.min(seq.len() / 2).min(seq.len().saturating_sub(10))).any(|period| {
        let matches = seq
            .iter()
            .zip(&seq[period..])
            .filter(|(a, b)| a == b)
            .count();
        matches * 100 >= (seq.len() - period) * 65
    })
}

/// Requires supporting alignments to concentrate near the physical read end.
fn terminal_support(seq: &[u8], windows: &[&[u8]], end: End, edits: usize) -> (usize, usize) {
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

/// Counts each exact k-mer once per window, excludes short tandem repeats,
/// sorts by count descending then code ascending, and retains at most `top`.
fn top_kmers(windows: &[&[u8]], k: usize, top: usize) -> Vec<(u64, u32)> {
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

/// Returns, per window, whether `pattern` occurs within `max_edits`.
fn windows_with(
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
/// removes the k-mers of that path's outermost support run from `nodes` so
/// the next round is forced onto a different path. Sequence past a layer
/// boundary stays available to later paths. Every path is returned; support
/// floors apply to the validated candidates.
fn peel_paths(mut nodes: Vec<(u64, u32)>, k: usize, end: End) -> Vec<(Vec<u8>, Vec<u32>)> {
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

/// Primer families reconstructed at one conserved insert boundary.
struct PrimerBoundary {
    /// Supported primer sequences, oriented as they occur in the read.
    primers: Vec<Vec<u8>>,
    /// First `START_K` bases of the conserved insert, in read orientation.
    insert_start: Vec<u8>,
}

/// Uses recurrent unprimed read starts to locate a conserved insert boundary.
/// A primer must be independently supported upstream of that boundary; the
/// conserved insert itself is excluded from the inferred trimming sequence.
/// `windows` read outward from the `end` boundary, reversed for the 3' end.
fn contrast_boundary(consensus: &[u8], windows: &[&[u8]], end: End) -> Option<PrimerBoundary> {
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
fn predominantly_insert(
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
fn supported_termination(word: &[u8], windows: &[&[u8]], end: End) -> bool {
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
const RUN_JUMP: u32 = 4;

/// K-mers of consistent support on each side of a boundary. A support ramp
/// without a plateau, such as an eroded adapter start, is one run.
const PLATEAU: usize = 8;

/// K-mers spanning a boundary whose support ramps between the plateaus,
/// as k-mers overlapping both a barcode and its flank do.
const RAMP: usize = 4;

/// K-mers over which the changed support level must persist. A
/// sequencing-error dip spans at most `KMER_K` k-mers; a barcode spans more.
const RUN_PERSIST: usize = 24;

/// Returns the base span `[lo, hi)` of the outermost support run of a
/// consensus and whether the run ends at a rise in support. `weights` holds
/// the support of each k-mer by start base. Scanning inward from the read
/// end, the run ends at the first boundary between two flat plateaus whose
/// supports differ `RUN_JUMP`-fold, where the new level persists for
/// `RUN_PERSIST` k-mers or to the end. Within the run, k-mers below
/// `BOUNDARY_SUPPORT` of the run peak are trimmed from both sides.
fn supported_span(weights: &[u32], end: End) -> (usize, usize, bool) {
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

/// An assembled end candidate: sequence, support, whether an independent
/// insert boundary was found, the summed original k-mer support of the
/// retained span, and whether the assembly window could not bound the
/// insert-facing side.
type Candidate = (Vec<u8>, f64, bool, u64, bool);

/// Assembles one end's candidates and validates support and insert boundaries.
/// `contrast` enables primer reconstruction from recurrent unprimed window
/// starts, which describe an insert boundary only in the first discovery
/// layer.
fn assemble(windows: &[&[u8]], base: &AdapterConfig, end: End, contrast: bool) -> Vec<Candidate> {
    if windows.len() < 3 {
        return Vec::new();
    }
    // K-mer encoding and approximate matching operate on uppercase DNA.
    // Inference owns normalized copies and does not modify pipeline records.
    let upper: Vec<Vec<u8>> = windows.iter().map(|w| w.to_ascii_uppercase()).collect();
    let windows: Vec<&[u8]> = upper.iter().map(Vec::as_slice).collect();
    let windows = windows.as_slice();

    let assembly_windows: Vec<&[u8]> = windows
        .iter()
        .map(|w| match end {
            End::Five => &w[..WINDOW_LEN.min(w.len())],
            End::Three => &w[w.len().saturating_sub(WINDOW_LEN)..],
        })
        .collect();
    let exact = top_kmers(&assembly_windows, KMER_K, TOP_KMERS);
    if exact
        .first()
        .is_none_or(|&(_, count)| count < MIN_SUPPORT_WINDOWS as u32)
    {
        return Vec::new();
    }
    // Validation uses windows distributed across the complete sample.
    let recount = stride_sample(windows, RECOUNT_WINDOWS);
    let n_recount = recount.len();

    let original_weights: std::collections::HashMap<u64, u32> = exact.iter().copied().collect();
    let weighted = exact;
    let mut out = Vec::new();
    let mut insert_starts = Vec::new();
    let mut primer_edges = vec![None; recount.len()];
    let mut primer_searcher = crate::adapter::search::new_searcher_fwd();
    let mut known_primers = std::collections::HashSet::new();
    // Contrast reads the validation windows outward from the end boundary.
    let reversed: Vec<Vec<u8>> = if contrast && end == End::Three {
        recount
            .iter()
            .map(|w| w.iter().rev().copied().collect())
            .collect()
    } else {
        Vec::new()
    };
    let oriented: Vec<&[u8]> = if reversed.is_empty() {
        recount.clone()
    } else {
        reversed.iter().map(Vec::as_slice).collect()
    };
    for (cons, _) in peel_paths(weighted, KMER_K, end) {
        tracing::debug!(sequence = %String::from_utf8_lossy(&cons), "Assembled end candidate");
        // Original weights show the support of every k-mer of the path,
        // including k-mers earlier peels removed from the graph.
        let weights: Vec<u32> = cons
            .windows(KMER_K)
            .map(|w| {
                encode_kmer(w)
                    .and_then(|code| original_weights.get(&code).copied())
                    .unwrap_or(0)
            })
            .collect();
        let (lo, hi, rises) = supported_span(&weights, end);
        let span = &weights[lo..=hi - KMER_K];
        let peak = span.iter().copied().max().unwrap_or(0);
        // The path weight measures completeness: a fragment running into
        // the insert or a sequencing variant carries weak k-mers.
        let weight: u64 = span.iter().map(|&w| u64::from(w)).sum();
        if peak < MIN_SUPPORT_WINDOWS as u32 {
            continue;
        }
        let trimmed = cons[lo..hi].to_vec();
        if is_repetitive(&trimmed) {
            continue;
        }
        // A variant or fragment of a stronger validated candidate merges
        // into it and needs no validation of its own.
        if out.iter().any(|(seq, _, _, other, _): &Candidate| {
            *other >= weight && same_adapter(&trimmed, seq, base.error_rate)
        }) {
            continue;
        }
        if !known_primers.is_empty()
            && predominantly_insert(
                &trimmed,
                &recount,
                end,
                edit_budget(base.error_rate, trimmed.len()),
                &primer_edges,
                &mut primer_searcher,
            )
        {
            continue;
        }
        let trimmed = if peak as usize * 4 < windows.len() {
            polish_consensus(&trimmed, &recount)
        } else {
            trimmed
        };
        let contrast = contrast
            .then(|| contrast_boundary(&cons, &oriented, end))
            .flatten();
        let has_insert_boundary = contrast.is_some();
        let mut unbounded = false;
        let mut measured = None;
        let sequences = if let Some(boundary) = contrast {
            insert_starts.push(boundary.insert_start);
            boundary.primers
        } else {
            let word = match end {
                End::Five => &cons[hi - KMER_K..hi],
                End::Three => &cons[lo..lo + KMER_K],
            };
            let code = encode_kmer(word).unwrap();
            let boundary_weight = original_weights[&code];
            let mask = (1u64 << (2 * KMER_K)) - 1;
            let (continuation, next) = (0..4)
                .map(|base| {
                    let next = match end {
                        End::Five => ((code << 2) | base) & mask,
                        End::Three => (code >> 2) | (base << (2 * (KMER_K - 1))),
                    };
                    (original_weights.get(&next).copied().unwrap_or(0), next)
                })
                .max()
                .unwrap_or((0, 0));
            let observed_boundary = match end {
                End::Five => hi < cons.len(),
                End::Three => lo > 0,
            };
            // Support drops at an insert boundary. It rises where a layer
            // such as a barcode joins a shared downstream layer: between the
            // support plateaus of the path, or where the downstream k-mer
            // occurs in `RUN_JUMP` times more windows than contain the
            // candidate at all. Assembly windows shorter than the layer
            // stack depress the downstream k-mer count; a twofold excess in
            // those counts is confirmed on the validation windows, where a
            // variant path joining its own family never reaches twofold. A
            // boundary inside one technical sequence changes support little.
            let drop = continuation * 2 <= boundary_weight;
            let (present, anchored) = terminal_support(
                &trimmed,
                &recount,
                end,
                edit_budget(base.error_rate, trimmed.len()),
            );
            measured = Some((present, anchored));
            let rise = rises
                || present.saturating_mul(RUN_JUMP as usize) <= continuation as usize
                || (continuation >= 2 * boundary_weight && {
                    let downstream = windows_containing(
                        &mut primer_searcher,
                        &decode_kmer(next, KMER_K),
                        &recount,
                        edit_budget(base.error_rate, KMER_K),
                    );
                    downstream as usize >= 2 * present
                });
            let terminated = drop && supported_termination(word, &recount, end);
            let bounded = (observed_boundary || cons.len() < LMAX) && (rise || terminated);
            if !bounded {
                tracing::debug!(sequence = %String::from_utf8_lossy(&trimmed), observed_boundary,
                    boundary_weight, continuation, present, drop, terminated,
                    "Recurrent sequence has no supported insert boundary");
            }
            unbounded = !bounded;
            vec![trimmed]
        };
        for trimmed in sequences {
            if trimmed.len() < MIN_PATTERN_LEN || is_repetitive(&trimmed) {
                continue;
            }
            // Presence counts each supporting window once, including reads
            // whose sequencing errors disrupted individual exact k-mers.
            let k_cons = edit_budget(base.error_rate, trimmed.len());
            let (present, anchored) = measured
                .take()
                .unwrap_or_else(|| terminal_support(&trimmed, &recount, end, k_cons));
            let support = present as f64 / n_recount as f64;
            tracing::debug!(sequence = %String::from_utf8_lossy(&trimmed), present, anchored, peak, has_insert_boundary, "Validated end candidate");
            if present >= MIN_SUPPORT_WINDOWS
                && support >= KEEP_SUPPORT
                && anchored * 100 >= present * ANCHORED_PERCENT
            {
                if has_insert_boundary && known_primers.insert(trimmed.clone()) {
                    for hit in primer_searcher.search_texts(&trimmed, &recount, k_cons) {
                        let edge = match end {
                            End::Five => hit.text_end,
                            End::Three => hit.text_start,
                        };
                        primer_edges[hit.text_idx] = Some((edge, trimmed.len() / 2));
                    }
                }
                out.push((trimmed, support, has_insert_boundary, weight, unbounded));
            }
        }
    }
    out.retain(|(seq, _, boundary, _, _)| {
        *boundary
            || !insert_starts
                .iter()
                .any(|start| seq.windows(start.len()).any(|word| word == start))
    });
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    if primer_edges.iter().any(Option::is_some) {
        out.retain(|(seq, _, boundary, _, _)| {
            if *boundary {
                return true;
            }
            !predominantly_insert(
                seq,
                &recount,
                end,
                edit_budget(base.error_rate, seq.len()),
                &primer_edges,
                &mut searcher,
            )
        });
    }
    out
}

/// Per-read explained depth from each physical end, in bases.
struct Boundaries {
    five: Vec<usize>,
    three: Vec<usize>,
}

impl Boundaries {
    fn new(sample: &[&[u8]]) -> Self {
        Self {
            five: vec![0; sample.len()],
            three: sample.iter().map(|r| r.len()).collect(),
        }
    }
}

/// Windows of one read end, each paired with its sample index.
type EndWindows<'a> = Vec<(usize, &'a [u8])>;

/// Returns the unexplained windows of at most `w` bases inward from each
/// boundary, paired with their sample indices. Reads without unexplained
/// bases contribute no window.
fn layer_windows<'a>(
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
fn background_windows<'a>(sample: &[&'a [u8]], bounds: &Boundaries, w: usize) -> Vec<&'a [u8]> {
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
/// past every further hit the window holds. Only reads flagged in `active`
/// are searched. Returns the reads whose boundary moved, or an empty vector
/// when none did.
fn advance_boundaries(
    sample: &[&[u8]],
    bounds: &mut Boundaries,
    patterns: &[Vec<u8>],
    error_rate: f64,
    active: &[bool],
    repeat: bool,
) -> Vec<bool> {
    let (five_w, three_w) = layer_windows(sample, bounds, 2 * WINDOW_LEN);
    let five_w: EndWindows = five_w.into_iter().filter(|(i, _)| active[*i]).collect();
    let three_w: EndWindows = three_w.into_iter().filter(|(i, _)| active[*i]).collect();
    let five_texts: Vec<&[u8]> = five_w.iter().map(|(_, w)| *w).collect();
    let three_texts: Vec<&[u8]> = three_w.iter().map(|(_, w)| *w).collect();
    let mut five_next = bounds.five.clone();
    let mut three_next = bounds.three.clone();
    // Hit spans per window, measured inward from the boundary.
    let mut five_hits: Vec<Vec<(usize, usize)>> = vec![Vec::new(); five_texts.len()];
    let mut three_hits: Vec<Vec<(usize, usize)>> = vec![Vec::new(); three_texts.len()];
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    // Patterns of one length, such as a barcode panel, share one tiled
    // search per window; other patterns are searched one strand at a time.
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
                    |_, start, end, _| hits.push((start, end)),
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
                    |_, start, end, _| hits.push((text.len() - end, text.len() - start)),
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
                    five_hits[hit.text_idx].push((hit.text_start, hit.text_end));
                }
                for hit in searcher.search_texts(&strand, &three_texts, k) {
                    let len = three_texts[hit.text_idx].len();
                    three_hits[hit.text_idx].push((len - hit.text_end, len - hit.text_start));
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
fn inner_depth(
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
fn outer_depth(
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
fn symmetry_cut(
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

/// Returns whether two candidates reconstruct one family: the same adapter
/// within the error budget, or reconstructions whose longest common
/// substring covers `FAMILY_OVERLAP_PERCENT` of the shorter, as remnants of
/// one adapter that differ in length do.
fn same_family(a: &[u8], b: &[u8], error_rate: f64) -> bool {
    same_adapter(a, b, error_rate)
        || longest_common_substring(a, b) * 100 >= a.len().min(b.len()) * FAMILY_OVERLAP_PERCENT
}

/// Returns whether the shorter of two candidates is a fragment of the
/// longer: one family by `same_family`, with less than `MIN_PATTERN_LEN` of
/// the shorter outside their longest common substring, as an end of the
/// longer extended by a few bases the assembly could not support is.
fn fragment_of(a: &[u8], b: &[u8], error_rate: f64) -> bool {
    let short = a.len().min(b.len());
    same_adapter(a, b, error_rate) || {
        let shared = longest_common_substring(a, b);
        shared * 100 >= short * FAMILY_OVERLAP_PERCENT && short - shared < MIN_PATTERN_LEN
    }
}

/// Length of the longest common substring of `a` and `b`.
fn longest_common_substring(a: &[u8], b: &[u8]) -> usize {
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

/// Removes from each candidate the prefix and suffix that occur in a
/// better supported candidate of another family in the same layer and end.
/// Sequence shared with a better supported neighbour, such as the flank
/// around a barcode, belongs to that neighbour's layer; a fragment of the
/// candidate itself is never stronger. Candidates left shorter than
/// `MIN_PATTERN_LEN` are dropped.
fn strip_shared_ends(candidates: &mut Vec<Candidate>, error_rate: f64) {
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
fn with_variable_layer(
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
/// its member count as support and weight.
fn variable_layer(
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
fn cluster_sequences(
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

/// A merged layer family ready for variant suppression: sequence, support,
/// supporting windows, their count, whether it carries an insert boundary,
/// its end, and its path weight.
type Ranked = (Vec<u8>, f64, Vec<bool>, usize, bool, End, u64);

/// Discovers supported technical sequences in layers from each read end.
/// Known sequences in `base` explain the outermost layers first; each
/// accepted layer moves the boundary inward and the next layer is assembled
/// from the unexplained sequence. Equivalent assemblies share one trimming
/// pattern. Catalog and supplied FASTA entries provide names only after the
/// inferred boundaries are fixed.
pub fn discover(sample: &[&[u8]], base: &AdapterConfig) -> Vec<InferredAdapter> {
    let mut bounds = Boundaries::new(sample);
    let known: Vec<Vec<u8>> = base
        .adapters
        .iter()
        .map(|a| a.seq.to_ascii_uppercase())
        .collect();
    if !known.is_empty() {
        let mut active = vec![true; sample.len()];
        for _ in 0..MAX_LAYERS {
            active =
                advance_boundaries(sample, &mut bounds, &known, base.error_rate, &active, true);
            if active.is_empty() {
                break;
            }
        }
    }
    let physical = Boundaries::new(sample);
    let (five_p, three_p) = layer_windows(sample, &physical, MIRROR_WINDOW);
    let five_phys: Vec<&[u8]> = stride_sample(
        &five_p.iter().map(|(_, w)| *w).collect::<Vec<_>>(),
        RECOUNT_WINDOWS,
    );
    let three_phys: Vec<&[u8]> = stride_sample(
        &three_p.iter().map(|(_, w)| *w).collect::<Vec<_>>(),
        RECOUNT_WINDOWS,
    );

    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let mut distinct: Vec<(Vec<u8>, f64, usize, bool, bool)> = Vec::new();
    // An accepted primer marks the insert boundary at its end; no layer lies
    // beyond it.
    let mut open = [true, true];
    for layer in 0..MAX_LAYERS {
        let (mut five_w, mut three_w) = layer_windows(sample, &bounds, 2 * WINDOW_LEN);
        if !open[0] {
            five_w.clear();
        }
        if !open[1] {
            three_w.clear();
        }
        let five_texts: Vec<&[u8]> = five_w.iter().map(|(_, w)| *w).collect();
        let three_texts: Vec<&[u8]> = three_w.iter().map(|(_, w)| *w).collect();
        // Ranking statistics use windows distributed across the sample.
        let five_sample = stride_sample(&five_texts, RECOUNT_WINDOWS);
        let three_sample = stride_sample(&three_texts, RECOUNT_WINDOWS);
        let mut five = assemble(&five_texts, base, End::Five, layer == 0);
        let mut three = assemble(&three_texts, base, End::Three, layer == 0);
        strip_shared_ends(&mut five, base.error_rate);
        strip_shared_ends(&mut three, base.error_rate);
        let (five, five_variable) = with_variable_layer(
            &mut searcher,
            five,
            &five_texts,
            &five_sample,
            End::Five,
            base.error_rate,
        );
        let (three, three_variable) = with_variable_layer(
            &mut searcher,
            three,
            &three_texts,
            &three_sample,
            End::Three,
            base.error_rate,
        );
        let variable: Vec<Vec<u8>> = five_variable.into_iter().chain(three_variable).collect();
        let background = background_windows(sample, &bounds, WINDOW_LEN);
        let background = stride_sample(&background, RECOUNT_WINDOWS);

        // Insert stretches identified by their mirror also reject the graph
        // fragments that reconstruct part of them.
        let mut inserts: Vec<Vec<u8>> = Vec::new();
        let cut: Vec<(Vec<u8>, f64, bool, u64, End)> = five
            .into_iter()
            .map(|c| (c, End::Five))
            .chain(three.into_iter().map(|c| (c, End::Three)))
            .filter(|((seq, support, _, _, _), _)| {
                seq.len() >= MIN_PATTERN_LEN && (*support >= KEEP_SUPPORT || variable.contains(seq))
            })
            .filter_map(|((seq, support, boundary, weight, unbounded), end)| {
                let (own, opposite) = match end {
                    End::Five => (&five_phys, &three_phys),
                    End::Three => (&three_phys, &five_phys),
                };
                // A mirror at the opposite end is the only insert evidence
                // for a candidate the assembly window could not bound.
                let (cut, insert) = symmetry_cut(&seq, end, own, opposite, base.error_rate);
                if insert.len() >= MIN_PATTERN_LEN {
                    inserts.push(crate::adapter::reverse_complement(&insert));
                    inserts.push(insert);
                }
                let cut = cut?;
                if unbounded && cut.len() == seq.len() {
                    return None;
                }
                Some((cut, support, boundary, weight, end))
            })
            .collect();
        let mut candidates: Vec<(Vec<u8>, f64, u32, bool, u64, End)> = cut
            .into_iter()
            .filter(|(seq, _, _, _, _)| {
                !inserts
                    .iter()
                    .any(|insert| same_family(seq, insert, base.error_rate))
            })
            .filter_map(|(seq, support, boundary, weight, end)| {
                let count = windows_containing(
                    &mut searcher,
                    &seq,
                    &background,
                    edit_budget(base.error_rate, seq.len()),
                );
                if !background.is_empty() && count as f64 * 4.0 >= support * background.len() as f64
                {
                    return None;
                }
                let exact = windows_containing(&mut searcher, &seq, &five_sample, 0)
                    + windows_containing(&mut searcher, &seq, &three_sample, 0);
                tracing::debug!(sequence = %String::from_utf8_lossy(&seq), exact, "Layer candidate");
                Some((seq, support, exact, boundary, weight, end))
            })
            .collect();
        // Independently supported insert boundaries take precedence over graph
        // fragments that may include conserved insert sequence. Exact support
        // then distinguishes reconstructions from nearby sequencing-error variants.
        candidates.sort_by(|a, b| {
            b.3.cmp(&a.3)
                .then_with(|| {
                    if a.3 && b.3 {
                        b.2.cmp(&a.2)
                    } else {
                        std::cmp::Ordering::Equal
                    }
                })
                .then(b.4.cmp(&a.4))
                .then(b.2.cmp(&a.2))
                .then(b.1.total_cmp(&a.1))
                .then(b.0.len().cmp(&a.0.len()))
                .then(a.0.cmp(&b.0))
        });
        // Contained reconstructions and fragments merge into the heaviest
        // reconstruction of their family. A sequence supported `RUN_JUMP`
        // times better than the candidate containing it is the shared layer
        // of that candidate, not its fragment.
        let mut merged: Vec<(Vec<u8>, f64, u64, bool, End)> = Vec::new();
        for (seq, support, _, boundary, weight, end) in candidates {
            let matched = weight;
            // A family found again in a later layer, in reads its first
            // occurrence did not match, is not a new layer; a shorter
            // sequence contained in an earlier candidate is a layer that the
            // candidate fused with its neighbour. The heaviest
            // reconstruction represents a family with its own support.
            if distinct.iter().any(|(other, _, _, _, _)| {
                other.len() <= seq.len() + edit_budget(base.error_rate, seq.len())
                    && same_adapter(&seq, other, base.error_rate)
            }) {
                continue;
            } else if let Some((other, previous, best, _, _)) =
                merged.iter_mut().find(|(other, previous, _, _, _)| {
                    support < *previous * f64::from(RUN_JUMP)
                        && *previous < support * f64::from(RUN_JUMP)
                        && fragment_of(&seq, other, base.error_rate)
                })
            {
                if matched > *best {
                    *other = seq;
                    *best = matched;
                    *previous = support;
                }
            } else {
                merged.push((seq, support, matched, boundary, end));
            }
        }
        // Sequencing variants of a family occupy the same reads as the
        // family; distinct families occupy different reads.
        let layer_texts: Vec<&[u8]> = five_sample.iter().chain(&three_sample).copied().collect();
        let mut ranked: Vec<Ranked> = merged
            .into_iter()
            .map(|(seq, support, matched, boundary, end)| {
                let k = edit_budget(base.error_rate, seq.len());
                let mut windows = windows_with(&mut searcher, &seq, &layer_texts, k);
                let mirror = crate::adapter::reverse_complement(&seq);
                for (window, hit) in
                    windows
                        .iter_mut()
                        .zip(windows_with(&mut searcher, &mirror, &layer_texts, k))
                {
                    *window |= hit;
                }
                let own = windows.iter().filter(|&&w| w).count();
                (seq, support, windows, own, boundary, end, matched)
            })
            .collect();
        // Complete reconstructions rank above fragments that run into the
        // insert, which match more windows over fewer exact bases.
        ranked.sort_by(|a, b| {
            b.6.cmp(&a.6)
                .then(b.3.cmp(&a.3))
                .then(b.0.len().cmp(&a.0.len()))
                .then(a.0.cmp(&b.0))
        });
        let mut accepted: Vec<Vec<u8>> = Vec::new();
        let mut accepted_windows: Vec<Vec<bool>> = Vec::new();
        let mut accepted_ends: Vec<End> = Vec::new();
        for (seq, support, windows, own, boundary, end, matched) in ranked {
            // Members of a variable layer share their reads with the
            // constant layer behind them by construction.
            let member = variable
                .iter()
                .any(|v| same_adapter(&seq, v, base.error_rate));
            let shared = accepted_windows
                .iter()
                .map(|other| {
                    windows
                        .iter()
                        .zip(other)
                        .filter(|(a, b)| **a && **b)
                        .count()
                })
                .max()
                .unwrap_or(0);
            let variant = !member && shared * 100 >= own * VARIANT_OVERLAP_PERCENT;
            tracing::debug!(sequence = %String::from_utf8_lossy(&seq), matched, own, shared, variant, "Ranked candidate");
            if variant {
                continue;
            }
            let flush = layer == 0
                && known.is_empty()
                && [(End::Five, &five_phys), (End::Three, &three_phys)]
                    .iter()
                    .any(|(end, phys)| {
                        outer_depth(&mut searcher, &seq, phys, *end, base.error_rate)
                            <= ADAPTER_FLUSH
                    });
            if boundary {
                open[usize::from(end == End::Three)] = false;
            }
            accepted_ends.push(end);
            accepted.push(seq.clone());
            if !member {
                accepted_windows.push(windows);
            }
            distinct.push((seq, support, layer, flush, member));
        }
        tracing::debug!(
            layer = layer + 1,
            accepted = accepted.len(),
            "Discovery layer"
        );
        // An end without a layer here has none deeper either.
        for (index, end) in [End::Five, End::Three].into_iter().enumerate() {
            if !accepted_ends.contains(&end) {
                open[index] = false;
            }
        }
        if accepted.is_empty()
            || advance_boundaries(
                sample,
                &mut bounds,
                &accepted,
                base.error_rate,
                &vec![true; sample.len()],
                false,
            )
            .is_empty()
        {
            break;
        }
    }

    let refs = crate::adapter::preset::preset(crate::adapter::preset::Kit::ALL);
    let name_refs: Vec<Adapter> = refs
        .into_iter()
        .chain(base.adapters.iter().cloned())
        .collect();
    distinct.sort_by(|a, b| a.2.cmp(&b.2).then(b.1.total_cmp(&a.1)).then(a.0.cmp(&b.0)));
    distinct
        .into_iter()
        .enumerate()
        .map(|(i, (seq, support, layer, flush, member))| {
            let name_hits = name_against(&seq, &name_refs, base.error_rate);
            let role = if flush {
                Role::Adapter
            } else if member {
                Role::Barcode
            } else {
                Role::Primer
            };
            InferredAdapter {
                adapter: Adapter {
                    name: format!("inferred_{}", i + 1),
                    seq: seq.clone(),
                    role,
                },
                assembled_seq: seq,
                support,
                name_hits,
                layer,
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
        infer_with_known(reads, Vec::new())
    }

    fn infer_with_known(reads: &[Vec<u8>], known: Vec<Adapter>) -> Vec<InferredAdapter> {
        let sample: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
        discover(
            &sample,
            &AdapterConfig {
                adapters: known,
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
    fn recovers_minority_adapter_beside_a_dominant_family() {
        let common = random_bases(47319, 31);
        let rare = random_bases(18971, 43);
        let reads: Vec<Vec<u8>> = (0..2000)
            .map(|i| {
                let mut read = match i % 100 {
                    0..=69 => common.clone(),
                    70..=74 => rare.clone(),
                    _ => Vec::new(),
                };
                read.extend(random_bases(18231 + i, 420));
                read
            })
            .collect();
        let found = infer_owned(&reads);
        assert!(found.iter().any(|d| d.adapter.seq == common), "{found:?}");
        assert!(found.iter().any(|d| d.adapter.seq == rare), "{found:?}");
        assert_eq!(found.len(), 2, "{found:?}");
    }

    #[test]
    fn tandem_repeat_ends_do_not_produce_adapters() {
        let reads: Vec<Vec<u8>> = (0..600)
            .map(|i| {
                let mut read = b"TTAGGG".repeat(10);
                if i % 3 == 0 {
                    read[19] = b'C';
                }
                read.extend(random_bases(87123 + i, 400));
                read
            })
            .collect();
        let found = infer_owned(&reads);
        assert!(found.is_empty(), "{found:?}");
    }

    #[test]
    fn minority_primer_uses_unprimed_insert_boundaries() {
        let primer = random_bases(88971, 24);
        let anchor = random_bases(34723, 140);
        let reads: Vec<Vec<u8>> = (0..2000)
            .map(|i| {
                let mut read = if i % 20 == 0 {
                    primer.clone()
                } else {
                    Vec::new()
                };
                read.extend_from_slice(&anchor);
                read.extend(random_bases(8921 + i, 300));
                read
            })
            .collect();
        let found = infer_owned(&reads);
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].adapter.seq, primer);
    }

    #[test]
    fn kmer_support_counts_independent_windows() {
        let seed = b"ACGTCAGTGCATGACT";
        let repeated = seed.repeat(5);
        let windows = vec![repeated.as_slice(), seed.as_slice()];
        let counts = top_kmers(&windows, 16, 500);
        assert_eq!(
            counts
                .iter()
                .find(|(key, _)| *key == encode_kmer(seed).unwrap())
                .unwrap()
                .1,
            2
        );
        assert!(top_kmers(&[b"TTAGGGTTAGGGTTAGGG"], 16, 500).is_empty());
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
    fn variable_conserved_inserts_extend_beyond_assembly_windows() {
        let anchors: Vec<_> = (0..4)
            .map(|i| (random_bases(7877 + i, 120), random_bases(4512 + i, 120)))
            .collect();
        let reads: Vec<Vec<u8>> = (0..2000)
            .map(|i| {
                let (prefix, suffix) = &anchors[i % anchors.len()];
                let mut read = prefix.clone();
                read.extend(random_bases(915 + i as u64, 300));
                read.extend_from_slice(suffix);
                let mutations = random_bases(54371 + i as u64, read.len() * 5);
                let mut mutated = Vec::new();
                for (&base, event) in read.iter().zip(mutations.chunks_exact(5)) {
                    match encode_kmer(&event[..4]).unwrap() {
                        0..=4 => mutated.push(event[4]),
                        5..=6 => {},
                        7 => mutated.extend_from_slice(&[base, event[4]]),
                        _ => mutated.push(base),
                    }
                }
                mutated
            })
            .collect();
        let found = infer_owned(&reads);
        assert!(found.is_empty(), "{found:?}");
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
        assert_eq!(found.len(), 2, "{found:?}");
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
        let found = contrast_boundary(&insert[1..90], &windows, End::Five)
            .unwrap()
            .primers;
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
        let planted = b"ACGTCAGTGCATGACT";
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..50u8 {
            let mut wnd = planted.to_vec();
            // An uncalled separator prevents shifted copies of the seed.
            wnd.push(b'N');
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
        let paths = peel_paths(nodes, 4, End::Five);
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

    fn rc(seq: &[u8]) -> Vec<u8> {
        crate::adapter::reverse_complement(seq)
    }

    fn contains(hay: &[u8], needle: &[u8]) -> bool {
        needle.len() <= hay.len() && hay.windows(needle.len()).any(|w| w == needle)
    }

    #[test]
    fn known_adapter_precedes_a_discovered_primer_layer() {
        let adapter = random_bases(5011, 30);
        let primer = random_bases(5023, 22);
        let reads: Vec<Vec<u8>> = (0..600)
            .map(|i| {
                let mut read = adapter.clone();
                read.extend_from_slice(&primer);
                read.extend(random_bases(7000 + i, 400));
                read
            })
            .collect();
        let known = vec![Adapter {
            name: "known".into(),
            seq: adapter.clone(),
            role: Role::Adapter,
        }];
        let found = infer_with_known(&reads, known);
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].adapter.seq, primer);
        assert_eq!(found[0].adapter.role, Role::Primer);
        assert_eq!(found[0].layer, 0);
    }

    #[test]
    fn discovers_barcode_layer_between_flanks() {
        let adapter = random_bases(6011, 30);
        let flank1 = random_bases(6023, 16);
        let flank2 = random_bases(6031, 32);
        let barcodes: Vec<Vec<u8>> = (0..8).map(|i| random_bases(6100 + i, 24)).collect();
        let reads: Vec<Vec<u8>> = (0..1600)
            .map(|i| {
                let mut read = adapter.clone();
                read.extend_from_slice(&flank1);
                read.extend_from_slice(&barcodes[i % 8]);
                read.extend_from_slice(&flank2);
                read.extend(random_bases(9000 + i as u64, 400));
                read
            })
            .collect();
        let found = infer_owned(&reads);
        for d in &found {
            let technical = barcodes.iter().any(|b| {
                contains(
                    &[&adapter[..], &flank1, b, &flank2].concat(),
                    &d.adapter.seq,
                )
            });
            assert!(technical, "{:?}", String::from_utf8_lossy(&d.adapter.seq));
            let expected = if contains(&d.adapter.seq, &adapter[..16]) {
                Role::Adapter
            } else if barcodes.iter().any(|b| contains(&d.adapter.seq, &b[4..20])) {
                Role::Barcode
            } else {
                Role::Primer
            };
            assert_eq!(d.adapter.role, expected, "{d:?}");
        }
        assert!(
            found
                .iter()
                .any(|d| d.layer == 0 && contains(&d.adapter.seq, &adapter[..16]))
        );
        assert!(
            found
                .iter()
                .any(|d| d.layer > 0 && contains(&d.adapter.seq, &flank2)),
            "{found:?}"
        );
        let recovered = barcodes
            .iter()
            .filter(|b| {
                found
                    .iter()
                    .any(|d| d.adapter.role == Role::Barcode && contains(&d.adapter.seq, &b[4..20]))
            })
            .count();
        assert!(recovered >= 6, "{recovered} barcodes recovered: {found:?}");
        let cfg = AdapterConfig {
            adapters: found.into_iter().map(|d| d.adapter).collect(),
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 20,
            candidate_index: std::sync::OnceLock::new(),
        };
        for read in reads.iter().take(16) {
            let segments = crate::adapter::adapter_segments(read, &cfg);
            assert_eq!(segments, vec![(102, read.len())], "{segments:?}");
        }
    }

    #[test]
    fn conserved_insert_mirrored_at_truncated_read_ends_is_excluded() {
        let adapter = random_bases(7011, 30);
        let primer = random_bases(7023, 20);
        let conserved = random_bases(7031, 40);
        let reads: Vec<Vec<u8>> = (0..800)
            .map(|i| {
                let mut read = adapter.clone();
                read.extend_from_slice(&primer);
                read.extend_from_slice(&conserved);
                read.extend(random_bases(11000 + i as u64, 300));
                read.extend(rc(&conserved));
                read.extend_from_slice(&rc(&primer)[..7]);
                read
            })
            .collect();
        let found = infer_owned(&reads);
        assert!(!found.is_empty());
        let technical = [&adapter[..], &primer].concat();
        for d in &found {
            assert!(
                contains(&technical, &d.adapter.seq),
                "{:?}",
                String::from_utf8_lossy(&d.adapter.seq)
            );
        }
        assert!(
            found.iter().any(|d| contains(&d.adapter.seq, &adapter)),
            "{found:?}"
        );
    }
}
