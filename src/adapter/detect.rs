use super::search::{DnaSearcher, hits, new_searcher};
use super::{Adapter, MIN_PATTERN_LEN, Terminal, classify_terminal};

/// Below this many sampled reads, presence detection is unreliable; callers
/// skip it and use the full adapter set.
pub const MIN_SAMPLE_FOR_DETECTION: usize = 100;

/// Minimum sampled-read count for an adapter to be kept: 0.2% of the sample,
/// floored at 3 (so a single stray hit can't promote an adapter).
pub fn presence_min(sample_size: usize) -> usize {
    (sample_size / 500).max(3)
}

/// Whether `ad` would actually act on `window` — a terminal hit (trimmed) or,
/// when `split`, an interior hit (`cost <= k_mid`, split). Mirrors the
/// per-adapter body of `adapter_segments` so "present" == "does something".
fn adapter_present_in(
    searcher: &mut DnaSearcher,
    window: &[u8],
    ad: &Adapter,
    error_rate: f64,
    end_size: usize,
    split: bool,
) -> bool {
    let n = window.len();
    let len = ad.seq.len();
    if n == 0 || len < MIN_PATTERN_LEN {
        return false;
    }
    let end_size = end_size.min(n);
    let k_end = (error_rate * len as f64).floor() as usize;
    let k_mid = (0.5 * error_rate * len as f64).floor() as usize;
    for h in hits(searcher, &ad.seq, window, k_end) {
        match classify_terminal(h.start, h.end, n, end_size, ad.end) {
            // `Excise` (adapter within end_size of both ends on a short read)
            // acts either way: split when `split`, terminal-trim otherwise.
            Terminal::Five | Terminal::Three | Terminal::Excise => return true,
            Terminal::None => {
                if split && h.cost <= k_mid {
                    return true;
                }
            },
        }
    }
    false
}

/// Retain the adapters found (would-act) in at least `min_count` of the sampled
/// reads. Order is preserved.
pub fn present(
    sample: &[&[u8]],
    adapters: &[Adapter],
    error_rate: f64,
    end_size: usize,
    split: bool,
    min_count: usize,
) -> Vec<Adapter> {
    let mut searcher = new_searcher();
    let mut counts = vec![0usize; adapters.len()];
    for &seq in sample {
        for (i, ad) in adapters.iter().enumerate() {
            if adapter_present_in(&mut searcher, seq, ad, error_rate, end_size, split) {
                counts[i] += 1;
            }
        }
    }
    adapters
        .iter()
        .zip(counts)
        .filter(|(_, c)| *c >= min_count)
        .map(|(a, _)| a.clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::End;

    fn ad(name: &str, seq: &[u8], end: End) -> Adapter {
        Adapter {
            name: name.into(),
            seq: seq.to_vec(),
            end,
        }
    }

    #[test]
    fn presence_min_boundaries() {
        assert_eq!(presence_min(0), 3);
        assert_eq!(presence_min(1000), 3); // 1000/500 = 2 -> max(3,2)=3
        assert_eq!(presence_min(10000), 20);
    }

    #[test]
    fn keeps_present_drops_absent() {
        // Build 200 reads each starting with adapter P (present) and never
        // containing adapter Q. Q must be dropped; P kept.
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
        let kept = present(&seqs, &adapters, 0.2, 150, true, presence_min(seqs.len()));
        let names: Vec<&str> = kept.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, vec!["P"], "present kept, absent dropped");
    }

    #[test]
    fn adapter_present_in_terminal_and_interior() {
        let mut s = new_searcher();
        let a = ad("a", b"GGGGTTTTGGGGTTTTGGGG", End::Both);
        // terminal: adapter at read start.
        let mut term = a.seq.clone();
        term.extend_from_slice(&[b'A'; 60]);
        assert!(adapter_present_in(&mut s, &term, &a, 0.2, 150, false));
        // absent: pure-A read.
        assert!(!adapter_present_in(&mut s, &[b'A'; 80], &a, 0.2, 150, true));
        // interior (deep, split on): adapter in the middle of a long read.
        let mut inter = vec![b'A'; 300];
        inter.splice(150..150, a.seq.iter().copied());
        assert!(
            adapter_present_in(&mut s, &inter, &a, 0.2, 20, true),
            "interior found when split"
        );
        assert!(
            !adapter_present_in(&mut s, &inter, &a, 0.2, 20, false),
            "interior ignored when ends-only"
        );
    }
}
