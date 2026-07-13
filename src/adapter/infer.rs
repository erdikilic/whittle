//! Ab-initio adapter inference. Discovers adapters de novo from a
//! read sample using Porechop_ABI's published method: read-end k-mer counting,
//! a weighted de Bruijn graph, length-bounded heaviest-path assembly, iterative
//! peeling, boundary drop-trim, and presence-fraction confidence. Implemented
//! from the paper (not translated from GPL source). Pure and format-neutral.

use crate::adapter::search::{DnaSearcher, hits, new_searcher};
use crate::adapter::{Adapter, AdapterConfig, End, MIN_PATTERN_LEN};

/// k-mer length used for end-window counting and assembly graph nodes.
const KMER_K: usize = 16;

/// Number of top exact k-mers (by count) kept per end before reweighting.
const TOP_KMERS: usize = 500;

/// Length of the 5'/3' end window scanned per read for adapter discovery.
const WINDOW_LEN: usize = 100;

/// Edit-distance budget for the forward-only per-window presence recount.
const RECOUNT_EDITS: usize = 2;

/// Minimum presence-fraction support required to keep a discovered adapter.
/// Support is the fraction of sampled end windows containing the consensus
/// within its length-scaled edit budget. The threshold retains common library
/// adapters while excluding sparse barcode-specific sequences and background.
const KEEP_SUPPORT: f64 = 0.30;

/// Cap on the number of windows scanned per k-mer during the 2-error recount
/// (the confidence pass), bounding its cost on large samples.
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

/// Maximum length of the end-facing anchor used by conservative inference.
/// Two independent 16-mers are long enough to be specific in ordinary long
/// reads while avoiding the unidentifiable insert-facing tail of a recurrent
/// amplicon consensus. Terminal trimming still removes everything between the
/// physical read end and this anchor.
const CONSERVATIVE_ANCHOR_LEN: usize = 2 * KMER_K;

/// One discovered adapter with inference metadata. Convert to a bare `Adapter`
/// (dropping `support`/`name_hits`) only when building the trim config.
#[derive(Debug, Clone)]
pub struct InferredAdapter {
    /// Sequence actually used for trimming (or printed as the recommendation).
    pub adapter: Adapter,
    /// Complete recurrent consensus assembled before conservative anchoring.
    pub assembled_seq: Vec<u8>,
    pub support: f64,
    pub name_hits: Vec<(String, f32)>,
}

impl InferredAdapter {
    /// Number of insert-facing consensus bases deliberately excluded from the
    /// trimming anchor because their technical/biological status is unknown.
    pub fn uncertain_bases(&self) -> usize {
        self.assembled_seq
            .len()
            .saturating_sub(self.adapter.seq.len())
    }
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

/// Number of distinct `windows` with >=1 forward approximate occurrence of
/// `kmer` (edit distance <= `max_edits`). Each window counts at most once,
/// even if `kmer` occurs in it multiple times. `searcher` must be
/// forward-only (see `new_searcher_fwd`) so reverse-complement occurrences do
/// not inflate the count. Callers provide an already bounded window sample.
fn two_error_freq(
    searcher: &mut DnaSearcher,
    kmer: &[u8],
    windows: &[&[u8]],
    max_edits: usize,
) -> u32 {
    let mut present = 0u32;
    for &wnd in windows {
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

    // Prefer greater support, then the smaller code for deterministic ties.
    let weight_desc_code_asc = |&a: &usize, &b: &usize| {
        nodes[a]
            .1
            .cmp(&nodes[b].1)
            .then(nodes[b].0.cmp(&nodes[a].0))
    };

    // deterministic pick: heaviest unvisited candidate, tie -> smaller code.
    let pick = |cands: Option<&Vec<usize>>, visited: &[bool]| -> Option<usize> {
        cands?
            .iter()
            .copied()
            .filter(|&i| !visited[i])
            .max_by(weight_desc_code_asc)
    };

    // seed = single heaviest node (tie -> smaller code).
    let seed = (0..n).max_by(weight_desc_code_asc).unwrap();
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

/// True if `a` and `b` are the "same" adapter within `error_rate`: an
/// approximate occurrence of the shorter in the longer on either strand
/// (the both-strand searcher covers the RC case).
fn same_adapter(a: &[u8], b: &[u8], error_rate: f64) -> bool {
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if short.len() < MIN_PATTERN_LEN {
        return short == long;
    }
    let k = (error_rate * short.len() as f64).floor() as usize;
    let mut s = new_searcher();
    !hits(&mut s, short, long, k).is_empty()
}

/// Folds a sequence discovered at BOTH the 5' and 3' ends (per `same_adapter`)
/// into a single `End::Both` entry, so the matcher's nearest-end arbitration
/// (see `classify_terminal`) handles it rather than two independent
/// single-end entries. The rest keep their originating end tag.
fn merge_both_ends(
    five: Vec<Vec<u8>>,
    three: Vec<Vec<u8>>,
    error_rate: f64,
    aggressive: bool,
) -> Vec<(Vec<u8>, End)> {
    let mut out: Vec<(Vec<u8>, End)> = Vec::new();
    let mut three_used = vec![false; three.len()];
    for f in &five {
        if let Some(j) = three
            .iter()
            .enumerate()
            .position(|(j, t)| !three_used[j] && same_adapter(f, t, error_rate))
        {
            three_used[j] = true;
            // Conservative inference keeps the 5' representation so its
            // prefix remains the physical-end-facing side when we extract a
            // terminal anchor below. Aggressive inference preserves the old
            // full-consensus behaviour and keeps the longer reconstruction.
            let kept = if aggressive && three[j].len() > f.len() {
                three[j].clone()
            } else {
                f.clone()
            };
            out.push((kept, End::Both));
        } else {
            out.push((f.clone(), End::Five));
        }
    }
    for (j, t) in three.into_iter().enumerate() {
        if !three_used[j] {
            out.push((t, End::Three));
        }
    }
    out
}

/// Return only the physical-end-facing part of an assembled consensus. The
/// insert-facing extension is fundamentally ambiguous for reference-free
/// amplicon data: an unknown primer and a conserved marker-gene prefix can be
/// recurrent at exactly the same rate. For 5' (and merged candidates stored
/// in 5' orientation) the outer anchor is the prefix; for 3' it is the suffix.
fn conservative_terminal_anchor(seq: &[u8], end: End) -> Vec<u8> {
    if seq.len() <= CONSERVATIVE_ANCHOR_LEN {
        return seq.to_vec();
    }
    match end {
        End::Five | End::Both => seq[..CONSERVATIVE_ANCHOR_LEN].to_vec(),
        End::Three => seq[seq.len() - CONSERVATIVE_ANCHOR_LEN..].to_vec(),
    }
}

/// Best catalog matches for `seq` as `(name, percent_identity)`, sorted desc,
/// top 3, only >= 60%. Used to give an inferred adapter a human-readable name
/// when it corresponds to a known catalog entry.
fn name_against(seq: &[u8], refs: &[Adapter], error_rate: f64) -> Vec<(String, f32)> {
    let mut s = new_searcher();
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
        let k = (error_rate * short.len() as f64).ceil() as usize;
        if let Some(h) = hits(&mut s, short, long, k)
            .into_iter()
            .min_by_key(|h| h.cost)
        {
            let pct = 100.0 * (1.0 - h.cost as f32 / short.len() as f32);
            if pct >= 60.0 {
                named.push((r.name.clone(), pct));
            }
        }
    }
    named.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap().then(a.0.cmp(&b.0)));
    named.truncate(3);
    named
}

/// Deterministically sample <= `cap` windows spread across the WHOLE slice
/// (stride, not the first `cap`) so the recount/support frame isn't order-
/// biased (see `assemble`'s `recount` sample below). If `windows.len()
/// <= cap`, returns every window, in order (step == 1).
fn stride_sample<'a>(windows: &[&'a [u8]], cap: usize) -> Vec<&'a [u8]> {
    let step = windows.len().div_ceil(cap.max(1)).max(1);
    windows.iter().step_by(step).copied().collect()
}

/// One end's workflow: count -> reweight by 2-error freq -> peel -> drop-trim.
/// Returns (trimmed consensus, support) candidates for this end.
fn assemble(windows: &[&[u8]], base: &AdapterConfig) -> Vec<(Vec<u8>, f64)> {
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
    // Bound recount cost while sampling the complete input range uniformly.
    let recount = stride_sample(windows, RECOUNT_WINDOWS);
    let n_recount = recount.len();

    let mut fwd = crate::adapter::search::new_searcher_fwd();
    let weighted: Vec<(u64, u32)> = exact
        .iter()
        .map(|&(code, _)| {
            let kmer = decode_kmer(code, KMER_K);
            (
                code,
                two_error_freq(&mut fwd, &kmer, &recount, RECOUNT_EDITS),
            )
        })
        .filter(|&(_, w)| w > 0)
        .collect();
    let mut out = Vec::new();
    for (cons, profile) in peel_paths(weighted, KMER_K) {
        let (trimmed, _tprof) = drop_trim(&cons, &profile);
        if trimmed.len() < MIN_PATTERN_LEN {
            continue;
        }
        // Whole-consensus presence: what fraction of the same recount
        // sample actually contains this trimmed consensus (within an error
        // budget scaled to its own length), reusing the same forward
        // searcher and the same per-window presence counter
        // (`two_error_freq`) already used to reweight individual k-mers
        // above. Unlike a per-position profile statistic, this can't be
        // dragged down by an internal low-weight pocket inside an otherwise-
        // correct reconstruction.
        let k_cons = (base.error_rate * trimmed.len() as f64).floor() as usize;
        let present = two_error_freq(&mut fwd, &trimmed, &recount, k_cons);
        let support = present as f64 / n_recount as f64;
        out.push((trimmed, support));
    }
    out
}

/// Full ab-initio discovery workflow: per-end `assemble`, fold shared 5'/3'
/// discoveries into `End::Both` via `merge_both_ends`, drop anything too
/// short or too weakly supported, then name each survivor against the
/// built-in ONT catalog UNION `base.adapters` -- extra naming refs, e.g. the
/// user's `--adapter-fasta` entries under `AdapterInfer::ReportOnly` (see
/// `cli::parse`'s `trim_adapters`; empty under `Trim`, which rejects a FASTA
/// outright, and under a `ReportOnly` run with no FASTA). Deterministic
/// order: support desc, then sequence asc.
pub fn discover(sample: &[&[u8]], base: &AdapterConfig) -> Vec<InferredAdapter> {
    discover_with_policy(sample, base, false)
}

/// Discover adapters with an explicit boundary policy. Conservative mode
/// (the default exposed by [`discover`]) trims with a short physical-end-facing
/// anchor and never asserts that the complete recurrent consensus is
/// technical. Aggressive mode preserves the historical full-consensus output.
pub fn discover_with_policy(
    sample: &[&[u8]],
    base: &AdapterConfig,
    aggressive: bool,
) -> Vec<InferredAdapter> {
    let (five_w, three_w) = end_windows(sample, WINDOW_LEN);
    let five = assemble(&five_w, base);
    let three = assemble(&three_w, base);

    // A dual-end consensus inherits the strongest fuzzy-equivalent recovery
    // from either end, including reverse-complement representations.
    let support_of = |seq: &[u8]| -> f64 {
        five.iter()
            .chain(three.iter())
            .filter(|(s, _)| same_adapter(s, seq, base.error_rate))
            .map(|(_, sup)| *sup)
            .fold(0.0_f64, f64::max)
    };

    let merged = merge_both_ends(
        five.iter().map(|(s, _)| s.clone()).collect(),
        three.iter().map(|(s, _)| s.clone()).collect(),
        base.error_rate,
        aggressive,
    );

    // Use the ONT catalog and optional user entries only as naming references.
    let refs = crate::adapter::preset::preset_ont();
    let name_refs: Vec<Adapter> = if base.adapters.is_empty() {
        refs
    } else {
        refs.into_iter()
            .chain(base.adapters.iter().cloned())
            .collect()
    };

    // Candidate fields: sequence, end, support, and catalog matches.
    type Candidate = (Vec<u8>, End, f64, Vec<(String, f32)>);
    let mut candidates: Vec<Candidate> = Vec::new();
    for (assembled_seq, end) in merged.into_iter() {
        if assembled_seq.len() < MIN_PATTERN_LEN {
            continue;
        }
        let support = support_of(&assembled_seq);
        if support < KEEP_SUPPORT {
            continue;
        }
        let name_hits = name_against(&assembled_seq, &name_refs, base.error_rate);
        candidates.push((assembled_seq, end, support, name_hits));
    }
    // Stable presentation order: support descending, then sequence ascending.
    candidates.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap().then(a.0.cmp(&b.0)));
    // Assign fallback names after sorting so names match presentation order.
    candidates
        .into_iter()
        .enumerate()
        .map(|(i, (assembled_seq, end, support, name_hits))| {
            let name = name_hits
                .first()
                .map(|(n, _)| n.clone())
                .unwrap_or_else(|| format!("inferred_{}", i + 1));
            let seq = if aggressive {
                assembled_seq.clone()
            } else {
                conservative_terminal_anchor(&assembled_seq, end)
            };
            InferredAdapter {
                adapter: Adapter { name, seq, end },
                assembled_seq,
                support,
                name_hits,
            }
        })
        .collect()
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
    fn merge_folds_shared_sequence_to_both() {
        let a = b"ACGTACGTACGTACGT".to_vec();
        let five = vec![a.clone(), b"TTTTGGGGTTTTGGGG".to_vec()];
        let three = vec![a.clone()]; // same adapter seen at 3' too
        let merged = merge_both_ends(five, three, 0.2, false);
        // a -> Both; the 5'-only one stays Five; no 3'-only left.
        assert!(merged.iter().any(|(s, e)| s == &a && *e == End::Both));
        assert!(
            merged
                .iter()
                .any(|(s, e)| s == b"TTTTGGGGTTTTGGGG" && *e == End::Five)
        );
        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn aggressive_merge_keeps_longer_of_matched_both_end_pair() {
        // Aggressive fuzzy-equivalent end assemblies retain the longer consensus.
        let core = b"ACGTACGTACGTACGT".to_vec(); // 16bp truncated core (>= MIN_PATTERN_LEN)
        let mut longer = b"TT".to_vec();
        longer.extend_from_slice(&core);
        longer.extend_from_slice(b"TT"); // 20bp: `core` is an exact substring, so same_adapter-equal
        let five = vec![core.clone()];
        let three = vec![longer.clone()];
        let merged = merge_both_ends(five, three, 0.2, true);
        assert_eq!(
            merged,
            vec![(longer, End::Both)],
            "the longer (3') reconstruction must be kept, not the shorter 5' core"
        );
    }

    #[test]
    fn conservative_merge_keeps_five_prime_orientation() {
        let core = b"ACGTACGTACGTACGT".to_vec();
        let mut longer_three = b"TT".to_vec();
        longer_three.extend_from_slice(&core);
        longer_three.extend_from_slice(b"TT");
        let merged = merge_both_ends(vec![core.clone()], vec![longer_three], 0.2, false);
        assert_eq!(merged, vec![(core, End::Both)]);
    }

    #[test]
    fn conservative_anchor_uses_physical_end_facing_side() {
        let seq: Vec<u8> = (0..64).map(|i| b"ACGT"[i % 4]).collect();
        assert_eq!(
            conservative_terminal_anchor(&seq, End::Five),
            seq[..CONSERVATIVE_ANCHOR_LEN]
        );
        assert_eq!(
            conservative_terminal_anchor(&seq, End::Both),
            seq[..CONSERVATIVE_ANCHOR_LEN]
        );
        assert_eq!(
            conservative_terminal_anchor(&seq, End::Three),
            seq[seq.len() - CONSERVATIVE_ANCHOR_LEN..]
        );
    }

    #[test]
    fn conservative_anchor_does_not_pad_short_consensus() {
        let seq = b"AATGTACTTCGTTCAGTTACGTATTGCT";
        assert_eq!(conservative_terminal_anchor(seq, End::Five), seq);
        assert_eq!(conservative_terminal_anchor(seq, End::Three), seq);
    }

    #[test]
    fn name_against_matches_catalog_entry() {
        let refs = vec![Adapter {
            name: "SQK-TEST".into(),
            seq: b"ACGTACGTACGTACGT".to_vec(),
            end: End::Both,
        }];
        let hits = name_against(b"ACGTACGTACGTACGT", &refs, 0.2);
        assert_eq!(hits[0].0, "SQK-TEST");
        assert!((hits[0].1 - 100.0).abs() < 1e-3);
    }

    #[test]
    fn discover_recovers_planted_adapter_under_error() {
        // Plant a known catalog-like adapter at the 5' end of many synthetic reads,
        // inject ~10% substitution error, and require recovery within a small edit
        // distance + a catalog name hit. Deterministic pseudo-noise (no RNG): a
        // fixed permutation of error positions per read index.
        let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp (SQK-NSK007-like)
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..500usize {
            let mut read = adapter.to_vec();
            // deterministic genomic tail. A naive `(i*7 + j*13) % 4` formula
            // (as `i*31 + j*17` in the clean-reads test below) is linear in
            // `j` mod 4 and collapses to a phase-rotated ACGT tandem repeat:
            // a spurious signal present in 100% of reads that outweighs and
            // crowds out the real, noisy planted adapter. Use the same
            // splitmix64-style mix as the clean-reads test for a genuinely
            // non-periodic deterministic tail.
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            for _ in 0..120usize {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            // deterministic ~10% substitutions in the adapter region
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
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(!found.is_empty(), "at least one adapter discovered");
        // the top candidate should be a 5'/both adapter close to the planted seq.
        let top = &found[0];
        assert!(top.adapter.seq.len() >= MIN_PATTERN_LEN);
        // near-match to the planted adapter (fuzzy, since recovery is approximate)
        let mut s = new_searcher();
        let k = (0.25 * adapter.len() as f64).ceil() as usize;
        assert!(
            !hits(&mut s, &top.adapter.seq, adapter, k).is_empty()
                || !hits(&mut s, adapter, &top.adapter.seq, k).is_empty(),
            "recovered adapter is within ~25% edit distance of the planted one"
        );
    }

    #[test]
    fn discover_dual_end_adapter_gets_max_support() {
        // Plant `adapter` at the 5' end (with heavier substitutions -- weak
        // recovery) and its exact reverse complement at the 3' end (strong
        // recovery) of every read, so `merge_both_ends` folds the two
        // per-end discoveries into a single `End::Both` entry (per
        // The noisy 5' and exact 3' assemblies are fuzzy-equivalent. The merged
        // adapter must inherit the stronger 3' support.
        let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp
        let rc: Vec<u8> = adapter
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!("adapter is pure ACGT"),
            })
            .collect();
        let mut owned: Vec<Vec<u8>> = Vec::new();
        for i in 0..200usize {
            // 5' copy: deterministic substitutions at every 6th (shifted)
            // position -- weak but still independently recoverable.
            let mut read = adapter.to_vec();
            for p in (0..adapter.len()).step_by(6) {
                let q = (p + i) % adapter.len();
                read[q] = b"ACGT"[(read[q] as usize + 1) % 4];
            }
            // deterministic non-periodic genomic middle (same splitmix64
            // mix used by the other `discover_*` fixtures in this file).
            let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
            for _ in 0..150usize {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
            }
            // 3' copy: EXACT reverse complement, no error -- strong recovery.
            read.extend_from_slice(&rc);
            owned.push(read);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);

        let both = found
            .iter()
            .find(|d| d.adapter.end == End::Both)
            .expect("the shared 5'/3' adapter must be discovered as a single End::Both entry");

        // near-matches the planted adapter (fuzzy, since recovery is
        // approximate).
        let mut s = new_searcher();
        let k = (0.25 * adapter.len() as f64).ceil() as usize;
        assert!(
            !hits(&mut s, &both.adapter.seq, adapter, k).is_empty()
                || !hits(&mut s, adapter, &both.adapter.seq, k).is_empty(),
            "Both adapter (seq {:?}) must be within ~25% edit distance of the planted adapter",
            String::from_utf8_lossy(&both.adapter.seq)
        );

        // the reported support must reflect the stronger (3') end, not the
        // weaker 5' end alone (~0.18 -- see the sibling unmerged `Five`
        // entries this fixture also produces, at that same value, and
        // dropped independently since 0.18 < KEEP_SUPPORT).
        assert!(
            both.support > 0.7,
            "Both adapter's support ({}) must reflect the max across ends \
             (3' end recovers at ~1.0 here), not just the weaker 5' end alone (~0.18)",
            both.support
        );
    }

    #[test]
    fn discover_finds_nothing_in_clean_reads() {
        // SplitMix64-derived bases provide deterministic, nonperiodic background.
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
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(
            found.is_empty(),
            "no spurious adapter in clean reads (got {found:?})"
        );
    }

    #[test]
    fn stride_sample_is_identity_when_within_cap() {
        // len <= cap: every window is returned, in order (step == 1).
        let a: &[u8] = b"A";
        let b: &[u8] = b"C";
        let c: &[u8] = b"G";
        let windows: Vec<&[u8]> = vec![a, b, c];
        assert_eq!(stride_sample(&windows, 4), windows);
    }

    #[test]
    fn stride_sample_spans_the_whole_range_not_just_a_prefix() {
        // A four-element sample spans all 13 input positions.
        let bytes: Vec<u8> = (0..13u8).map(|i| b'A' + i).collect();
        let windows: Vec<&[u8]> = bytes.iter().map(std::slice::from_ref).collect();
        let sampled = stride_sample(&windows, 4);
        assert!(sampled.len() <= 4);
        // Expected indices: 0, 4, 8, 12.
        assert_eq!(
            sampled,
            vec![windows[0], windows[4], windows[8], windows[12]]
        );
        let last_idx = 12usize; // where the last sampled window actually lives
        assert!(
            last_idx >= (13usize * 2).div_ceil(3),
            "last sampled window must fall in the last third of the range, not a prefix"
        );
        assert_eq!(*sampled.last().unwrap(), windows[last_idx]);
    }

    // This 8001-read integration case is run on demand:
    //   cargo test --lib discover_is_not_order_biased_by_recount_window_cap -- --ignored
    #[test]
    #[ignore]
    fn discover_is_not_order_biased_by_recount_window_cap() {
        // The adapter occurs only in the latter half of the sample, requiring
        // the bounded recount to cover the complete input range.
        let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp, same as the other discover_* fixtures
        let n_clean = RECOUNT_WINDOWS + 1;
        let n_planted = RECOUNT_WINDOWS; // 4000

        // Deterministic non-periodic background, same splitmix64 mix used
        // throughout this file's other `discover_*` fixtures.
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
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(
            !found.is_empty(),
            "adapter present in a clear majority of reads after the first \
             RECOUNT_WINDOWS must still be discovered, not hidden by a \
             first-N window cap (got {found:?})"
        );
        let mut s = new_searcher();
        let k = (0.25 * adapter.len() as f64).ceil() as usize;
        assert!(
            found.iter().any(|d| {
                !hits(&mut s, &d.adapter.seq, adapter, k).is_empty()
                    || !hits(&mut s, adapter, &d.adapter.seq, k).is_empty()
            }),
            "discovered adapter(s) must include one within ~25% edit distance \
             of the planted adapter: {found:?}"
        );
    }

    #[test]
    fn discover_recovers_planted_adapter_from_lowercase_reads() {
        // Lowercase reads must produce the same inferred adapter as uppercase DNA.
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
            // Normalize the complete fixture to lowercase.
            let lower: Vec<u8> = read.iter().map(u8::to_ascii_lowercase).collect();
            owned.push(lower);
        }
        let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
        let base = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let found = discover(&sample, &base);
        assert!(
            !found.is_empty(),
            "lowercase reads must still be inferable (got {found:?})"
        );
        let top = &found[0];
        let mut s = new_searcher();
        let k = (0.25 * adapter.len() as f64).ceil() as usize;
        assert!(
            !hits(&mut s, &top.adapter.seq, adapter, k).is_empty()
                || !hits(&mut s, adapter, &top.adapter.seq, k).is_empty(),
            "discovered adapter (seq {:?}) must be within ~25% edit distance \
             of the (uppercase) planted adapter",
            String::from_utf8_lossy(&top.adapter.seq)
        );
        // The discovered sequence itself must be valid uppercase ACGT, not
        // carrying any lowercase byte through from the input.
        assert!(
            top.adapter.seq.iter().all(u8::is_ascii_uppercase),
            "discovered sequence must be uppercase: {:?}",
            String::from_utf8_lossy(&top.adapter.seq)
        );
    }
}
