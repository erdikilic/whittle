use crate::qual::{QualMode, read_quality};

#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub min_length: usize,
    pub max_length: usize,
    pub min_qual: f64,
    pub max_qual: f64,
    pub min_gc: Option<f64>,
    pub max_gc: Option<f64>,
    pub qual_mode: QualMode,
}

pub fn gc_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c'))
        .count();
    gc as f64 / seq.len() as f64
}

/// Cheapest-first, short-circuiting. Empty reads never pass.
pub fn passes(seq: &[u8], phred: &[u8], cfg: &FilterConfig) -> bool {
    let len = seq.len();
    if len == 0 || len < cfg.min_length || len > cfg.max_length {
        return false;
    }
    if cfg.min_qual > 0.0 || cfg.max_qual < 1000.0 {
        let q = read_quality(phred, cfg.qual_mode);
        if q < cfg.min_qual || q > cfg.max_qual {
            return false;
        }
    }
    if cfg.min_gc.is_some() || cfg.max_gc.is_some() {
        let gc = gc_fraction(seq);
        if gc < cfg.min_gc.unwrap_or(0.0) || gc > cfg.max_gc.unwrap_or(1.0) {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qual::QualMode;

    fn base() -> FilterConfig {
        FilterConfig {
            min_length: 1,
            max_length: usize::MAX,
            min_qual: 0.0,
            max_qual: 1000.0,
            min_gc: None,
            max_gc: None,
            qual_mode: QualMode::Mean,
        }
    }

    #[test]
    fn length_bounds() {
        let mut c = base();
        c.min_length = 4;
        c.max_length = 8;
        assert!(!passes(b"ATG", &[30, 30, 30], &c)); // too short
        assert!(passes(b"ATGCG", &[30; 5], &c));
        assert!(!passes(b"ATGCGATGC", &[30; 9], &c)); // too long
    }

    #[test]
    fn quality_bound_uses_mode() {
        let mut c = base();
        c.min_qual = 15.0;
        // arithmetic mean of [10,20] = 15.0 -> passes at threshold
        c.qual_mode = QualMode::Arithmetic;
        assert!(passes(b"AT", &[10, 20], &c));
        // prob-mean of [10,20] < 15 -> fails
        c.qual_mode = QualMode::Mean;
        assert!(!passes(b"AT", &[10, 20], &c));
    }

    #[test]
    fn gc_fraction_and_filter() {
        assert!((gc_fraction(b"GGCC") - 1.0).abs() < 1e-12);
        assert!((gc_fraction(b"ATAT") - 0.0).abs() < 1e-12);
        let mut c = base();
        c.min_gc = Some(0.4);
        c.max_gc = Some(0.6);
        assert!(passes(b"ATGC", &[30; 4], &c)); // 0.5
        assert!(!passes(b"AAAT", &[30; 4], &c)); // 0.0
    }

    #[test]
    fn empty_seq_rejected() {
        assert!(!passes(b"", &[], &base()));
    }
}
