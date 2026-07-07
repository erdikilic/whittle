//! Ab-initio adapter inference (Phase 2). Discovers adapters de novo from a
//! read sample using Porechop_ABI's published method: read-end k-mer counting,
//! a weighted de Bruijn graph, length-bounded heaviest-path assembly, iterative
//! peeling, boundary drop-trim, and presence-fraction confidence. Implemented
//! from the paper (not translated from GPL source). Pure and format-neutral.

use crate::adapter::Adapter;

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
// Not yet called from production code (only from the tests below); the
// assembly task wires this in and this allow comes off then.
#[allow(dead_code)]
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
}
