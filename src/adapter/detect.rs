//! Adapter presence detection over a sampled read prefix.
//!
//! Narrows a configured adapter set to the entries that would act on at least
//! `presence_min` of the sampled reads, so a large catalog is not searched in
//! full on every read.

use super::search::{AmbiguousSearcher, hits, new_ambiguous_searcher};
use super::{Adapter, Budget, MIN_PATTERN_LEN, Terminal, classify_terminal, normalized_read};
use rayon::prelude::*;

/// Sample size below which presence detection is unreliable; callers skip it
/// and use the full adapter set.
pub const MIN_SAMPLE_FOR_DETECTION: usize = 100;

/// Returns the minimum sampled-read count for an adapter to be kept: 0.2% of
/// the sample, floored at 3 so that a single stray hit cannot promote an
/// adapter.
pub fn presence_min(sample_size: usize) -> usize {
    (sample_size / 500).max(3)
}

/// Returns whether `ad` would act on `window`: a terminal hit (trimmed) or,
/// when `split`, an interior hit (`cost <= k_mid`, split). Searches the same
/// normalized text with the same budgets as `adapter_segments`, so presence is
/// defined as having an effect on the read.
fn adapter_present_in(
    searcher: &mut AmbiguousSearcher,
    window: &[u8],
    ad: &Adapter,
    budget: Budget,
    end_size: usize,
    split: bool,
) -> bool {
    let n = window.len();
    if n == 0 || ad.seq.len() < MIN_PATTERN_LEN {
        return false;
    }
    let window = normalized_read(window);
    let end_size = end_size.min(n);
    for h in hits(searcher, &ad.seq, &window, budget.k_end) {
        match classify_terminal(h.start, h.end, n, end_size, ad.end) {
            // `Excise` acts either way: split when `split`, terminal-trim otherwise.
            Terminal::Five | Terminal::Three | Terminal::Excise => return true,
            Terminal::None => {
                if split && h.cost <= budget.k_mid {
                    return true;
                }
            },
        }
    }
    false
}

/// Retains the adapters that would act on at least `min_count` of the sampled
/// reads. Order is preserved.
pub fn present(
    sample: &[&[u8]],
    adapters: &[Adapter],
    error_rate: f64,
    end_size: usize,
    split: bool,
    min_count: usize,
    threads: usize,
) -> Vec<Adapter> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .expect("A positive Rayon worker count builds a pool");
    pool.install(|| {
        adapters
            .par_iter()
            .filter_map(|ad| {
                let mut searcher = new_ambiguous_searcher();
                let budget = Budget::new(ad.seq.len(), error_rate);
                let mut count = 0usize;
                for &seq in sample {
                    if adapter_present_in(&mut searcher, seq, ad, budget, end_size, split) {
                        count += 1;
                        if count >= min_count {
                            return Some(ad.clone());
                        }
                    }
                }
                None
            })
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::{AdapterConfig, End, adapter_segments};

    /// Builds an adapter from its parts.
    fn ad(name: &str, seq: &[u8], end: End) -> Adapter {
        Adapter {
            name: name.into(),
            seq: seq.to_vec(),
            end,
        }
    }

    /// Generates deterministic SplitMix64 bases with the same generator the
    /// inference fixtures use.
    fn splitmix_dna(seed: u64, len: usize) -> Vec<u8> {
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(seed);
        (0..len)
            .map(|_| {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                b"ACGT"[((z >> 62) & 0b11) as usize]
            })
            .collect()
    }

    /// The floor of 3 holds up to 1500 sampled reads; above that 0.2% applies.
    #[test]
    fn presence_min_boundaries() {
        assert_eq!(presence_min(0), 3);
        assert_eq!(presence_min(1000), 3); // 1000/500 = 2 -> max(3,2)=3
        assert_eq!(presence_min(10000), 20);
    }

    /// Over 200 reads that each start with adapter P and never contain adapter
    /// Q, P is kept and Q is dropped.
    #[test]
    fn keeps_present_drops_absent() {
        let p = b"GGGGTTTTGGGGTTTTGGGG"; // 20bp
        let q = b"ACGACGACGACGACGACGAC"; // 20bp, absent (and not P's revcomp)
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for _ in 0..200 {
            let mut r = p.to_vec();
            r.extend_from_slice(&[b'A'; 60]); // insert with no P/Q content
            reads.push(r);
        }
        let seqs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let adapters = vec![ad("P", p, End::Both), ad("Q", q, End::Both)];
        let kept = present(
            &seqs,
            &adapters,
            0.2,
            150,
            true,
            presence_min(seqs.len()),
            2,
        );
        let names: Vec<&str> = kept.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, vec!["P"], "Present kept, absent dropped");
    }

    /// A terminal hit counts, an absent adapter does not, and an interior hit
    /// counts only when splitting is on.
    #[test]
    fn adapter_present_in_terminal_and_interior() {
        let mut s = new_ambiguous_searcher();
        let a = ad("a", b"GGGGTTTTGGGGTTTTGGGG", End::Both);
        let budget = Budget::new(a.seq.len(), 0.2);
        // Terminal: adapter at the read start.
        let mut term = a.seq.clone();
        term.extend_from_slice(&[b'A'; 60]);
        assert!(adapter_present_in(&mut s, &term, &a, budget, 150, false));
        // Absent: pure-A read.
        assert!(!adapter_present_in(
            &mut s,
            &[b'A'; 80],
            &a,
            budget,
            150,
            true
        ));
        // Interior (deep, split on): adapter in the middle of a long read.
        let mut inter = vec![b'A'; 300];
        inter.splice(150..150, a.seq.iter().copied());
        assert!(
            adapter_present_in(&mut s, &inter, &a, budget, 20, true),
            "Interior found when split"
        );
        assert!(
            !adapter_present_in(&mut s, &inter, &a, budget, 20, false),
            "Interior ignored when ends-only"
        );
    }

    /// Sixty `N`s then random bases: no adapter is present, and detection agrees
    /// with the trimmer, which rewrites the run before searching.
    #[test]
    fn ambiguity_runs_in_reads_are_not_evidence() {
        let a = ad("a", b"GGGGTTTTGGGGTTTTGGGG", End::Both);
        let reads: Vec<Vec<u8>> = (0..200u64)
            .map(|i| {
                let mut r = vec![b'N'; 60];
                r.extend(splitmix_dna(i, 100));
                r
            })
            .collect();
        let seqs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let kept = present(
            &seqs,
            std::slice::from_ref(&a),
            0.2,
            150,
            true,
            presence_min(seqs.len()),
            2,
        );
        assert!(
            kept.is_empty(),
            "An N run is not adapter evidence: {kept:?}"
        );

        let cfg = AdapterConfig {
            adapters: vec![a.clone()],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let mut s = new_ambiguous_searcher();
        let budget = Budget::new(a.seq.len(), 0.2);
        for r in &reads {
            assert_eq!(adapter_segments(r, &cfg), vec![(0, r.len())]);
            assert!(!adapter_present_in(&mut s, r, &a, budget, 150, true));
        }

        // The adapter after the run: detection and the trimmer both act.
        let mut planted = vec![b'N'; 60];
        planted.extend_from_slice(&a.seq);
        planted.extend(splitmix_dna(9, 80));
        assert!(adapter_present_in(&mut s, &planted, &a, budget, 150, true));
        assert_ne!(adapter_segments(&planted, &cfg), vec![(0, planted.len())]);
    }
}
