//! Adapter presence detection over a sampled read prefix.
//!
//! Narrows a configured adapter set to the entries that would act on at least
//! `presence_min` of the sampled reads, so a large catalog is not searched in
//! full on every read.

use super::{Adapter, AdapterConfig, adapter_segments_tallied};
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

/// Retains the adapters of `cfg` that trim or excise at least `min_count` of
/// the sampled reads, running the same search passes the trimming pass runs.
/// Order is preserved.
pub fn present(
    sample: &[&[u8]],
    cfg: &AdapterConfig,
    min_count: usize,
    threads: usize,
) -> Vec<Adapter> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .expect("A positive Rayon worker count builds a pool");
    let n = cfg.adapters.len();
    let counts = pool.install(|| {
        sample
            .par_iter()
            .fold(
                || vec![0usize; n],
                |mut counts, &read| {
                    let mut acted = vec![false; n];
                    adapter_segments_tallied(read, cfg, &mut acted);
                    for (count, hit) in counts.iter_mut().zip(acted) {
                        *count += usize::from(hit);
                    }
                    counts
                },
            )
            .reduce(
                || vec![0usize; n],
                |mut a, b| {
                    for (x, y) in a.iter_mut().zip(b) {
                        *x += y;
                    }
                    a
                },
            )
    });
    cfg.adapters
        .iter()
        .zip(counts)
        .filter(|(_, count)| *count >= min_count)
        .map(|(adapter, _)| adapter.clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::{Role, adapter_segments};

    /// Builds an adapter from its parts.
    fn ad(name: &str, seq: &[u8], role: Role) -> Adapter {
        Adapter {
            name: name.into(),
            seq: seq.to_vec(),
            role,
        }
    }

    /// Builds a configuration at error rate 0.2 with the given end zone.
    fn cfg(adapters: Vec<Adapter>, end_size: usize, split: bool) -> AdapterConfig {
        AdapterConfig {
            adapters,
            error_rate: 0.2,
            end_size,
            split,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
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
        assert_eq!(presence_min(1000), 3);
        assert_eq!(presence_min(10000), 20);
    }

    /// Over 200 reads that each start with adapter P and never contain adapter
    /// Q, P is kept and Q is dropped.
    #[test]
    fn keeps_present_drops_absent() {
        let p = b"GGGGTTTTGGGGTTTTGGGG";
        let q = b"ACGACGACGACGACGACGAC";
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for _ in 0..200 {
            let mut r = p.to_vec();
            r.extend_from_slice(&[b'A'; 60]);
            reads.push(r);
        }
        let seqs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let c = cfg(
            vec![ad("P", p, Role::Adapter), ad("Q", q, Role::Adapter)],
            150,
            true,
        );
        let kept = present(&seqs, &c, presence_min(seqs.len()), 2);
        let names: Vec<&str> = kept.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, vec!["P"], "Present kept, absent dropped");
    }

    /// A terminal hit counts, an absent adapter does not, and an interior hit
    /// counts only when splitting is on.
    #[test]
    fn tally_marks_terminal_and_interior_hits() {
        let a = ad("a", b"GGGGTTTTGGGGTTTTGGGG", Role::Adapter);
        let mut term = a.seq.clone();
        term.extend_from_slice(&[b'A'; 60]);
        let mut acted = [false];
        adapter_segments_tallied(&term, &cfg(vec![a.clone()], 150, false), &mut acted);
        assert!(acted[0], "Terminal hit counts");

        let mut acted = [false];
        adapter_segments_tallied(&[b'A'; 80], &cfg(vec![a.clone()], 150, true), &mut acted);
        assert!(!acted[0], "An absent adapter does not count");

        let mut inter = vec![b'A'; 300];
        inter.splice(150..150, a.seq.iter().copied());
        let mut acted = [false];
        adapter_segments_tallied(&inter, &cfg(vec![a.clone()], 20, true), &mut acted);
        assert!(acted[0], "Interior found when split");
        let mut acted = [false];
        adapter_segments_tallied(&inter, &cfg(vec![a.clone()], 20, false), &mut acted);
        assert!(!acted[0], "Interior ignored when ends-only");
    }

    /// Sixty `N`s then random bases: no adapter is present, and detection agrees
    /// with the trimmer, which rewrites the run before searching.
    #[test]
    fn ambiguity_runs_in_reads_are_not_evidence() {
        let a = ad("a", b"GGGGTTTTGGGGTTTTGGGG", Role::Adapter);
        let reads: Vec<Vec<u8>> = (0..200u64)
            .map(|i| {
                let mut r = vec![b'N'; 60];
                r.extend(splitmix_dna(i, 100));
                r
            })
            .collect();
        let seqs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let c = cfg(vec![a.clone()], 150, true);
        let kept = present(&seqs, &c, presence_min(seqs.len()), 2);
        assert!(
            kept.is_empty(),
            "An N run is not adapter evidence: {kept:?}"
        );
        for r in &reads {
            assert_eq!(adapter_segments(r, &c), vec![(0, r.len())]);
        }

        // The adapter after the run: detection and the trimmer both act.
        let mut planted = vec![b'N'; 60];
        planted.extend_from_slice(&a.seq);
        planted.extend(splitmix_dna(9, 80));
        let mut acted = [false];
        adapter_segments_tallied(&planted, &c, &mut acted);
        assert!(acted[0]);
        assert_ne!(adapter_segments(&planted, &c), vec![(0, planted.len())]);
    }
}
