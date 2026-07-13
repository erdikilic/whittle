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

/// Why `check` dropped a segment. Both GC bounds (too low / too high) collapse
/// into `Gc` — the summary reports "GC out of range", not which side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    TooShort,
    TooLong,
    LowQuality,
    HighQuality,
    Gc,
}

/// Evaluate bounds cheapest-first and stop at the first rejection. `None`
/// indicates that the segment passes; empty segments are `TooShort` even when
/// `min_length` is zero.
///
/// Called once for each segment produced by trimming, so `seq` and `phred`
/// describe that segment rather than necessarily the complete input read.
pub fn check(seq: &[u8], phred: &[u8], cfg: &FilterConfig) -> Option<DropReason> {
    let gc = (cfg.min_gc.is_some() || cfg.max_gc.is_some()).then(|| gc_fraction(seq));
    check_metrics(seq.len(), phred, gc, cfg)
}

/// Filter from precomputed sequence metrics. The raw BAM fast path can obtain
/// length and GC directly from packed sequence views without materializing an
/// owned decoded sequence. `gc` is required only when a GC bound is active.
pub(crate) fn check_metrics(
    len: usize,
    phred: &[u8],
    gc: Option<f64>,
    cfg: &FilterConfig,
) -> Option<DropReason> {
    if len == 0 || len < cfg.min_length {
        return Some(DropReason::TooShort);
    }
    if len > cfg.max_length {
        return Some(DropReason::TooLong);
    }
    if cfg.min_qual > 0.0 || cfg.max_qual < 1000.0 {
        let q = read_quality(phred, cfg.qual_mode);
        if q < cfg.min_qual {
            return Some(DropReason::LowQuality);
        }
        if q > cfg.max_qual {
            return Some(DropReason::HighQuality);
        }
    }
    if cfg.min_gc.is_some() || cfg.max_gc.is_some() {
        let gc = gc.unwrap_or(0.0);
        if gc < cfg.min_gc.unwrap_or(0.0) || gc > cfg.max_gc.unwrap_or(1.0) {
            return Some(DropReason::Gc);
        }
    }
    None
}

/// Cheapest-first, short-circuiting. Empty reads never pass.
pub fn passes(seq: &[u8], phred: &[u8], cfg: &FilterConfig) -> bool {
    check(seq, phred, cfg).is_none()
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
        assert!(passes(b"ATGC", &[30; 4], &c)); // len == min_length, inclusive
        assert!(passes(b"ATGCGATG", &[30; 8], &c)); // len == max_length, inclusive
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
        assert!(passes(b"GCAAA", &[30; 5], &c)); // gc == min_gc (0.4), inclusive
    }

    #[test]
    fn empty_seq_rejected() {
        assert!(!passes(b"", &[], &base()));
    }

    #[test]
    fn check_reports_too_short() {
        let mut c = base();
        c.min_length = 4;
        assert_eq!(check(b"ATG", &[30, 30, 30], &c), Some(DropReason::TooShort));
        // Empty reads are TooShort regardless of min_length.
        let c0 = base();
        assert_eq!(check(b"", &[], &c0), Some(DropReason::TooShort));
    }

    #[test]
    fn check_reports_too_long() {
        let mut c = base();
        c.max_length = 4;
        assert_eq!(check(b"ATGCG", &[30; 5], &c), Some(DropReason::TooLong));
    }

    #[test]
    fn check_reports_low_and_high_quality() {
        let mut c = base();
        c.qual_mode = QualMode::Arithmetic;
        c.min_qual = 25.0;
        assert_eq!(check(b"AT", &[10, 20], &c), Some(DropReason::LowQuality));

        let mut c = base();
        c.qual_mode = QualMode::Arithmetic;
        c.max_qual = 12.0;
        assert_eq!(check(b"AT", &[10, 20], &c), Some(DropReason::HighQuality));
    }

    #[test]
    fn check_reports_gc_low_and_high() {
        let mut c = base();
        c.min_gc = Some(0.4);
        c.max_gc = Some(0.6);
        assert_eq!(check(b"AAAT", &[30; 4], &c), Some(DropReason::Gc)); // gc 0.0 < min
        assert_eq!(check(b"GGCC", &[30; 4], &c), Some(DropReason::Gc)); // gc 1.0 > max
    }

    #[test]
    fn check_none_when_passing() {
        assert_eq!(check(b"ACGT", &[30; 4], &base()), None);
    }

    #[test]
    fn passes_matches_check_is_none() {
        let mut c = base();
        c.min_length = 4;
        assert_eq!(
            passes(b"AT", &[30, 30], &c),
            check(b"AT", &[30, 30], &c).is_none()
        );
        assert_eq!(
            passes(b"ATGC", &[30; 4], &c),
            check(b"ATGC", &[30; 4], &c).is_none()
        );
    }
}
