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
// Not yet called from production code (only from the tests below); the
// k-mer counting task wires this in and this allow comes off then.
#[allow(dead_code)]
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
    }
}
