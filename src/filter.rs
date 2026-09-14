//! Post-trim segment filtering by length, quality and GC fraction.

use crate::qual::{QualMode, read_quality};

/// Bounds applied to every produced segment.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// Minimum segment length, inclusive.
    pub min_length: usize,
    /// Maximum segment length, inclusive.
    pub max_length: usize,
    /// Minimum read quality, inclusive.
    pub min_qual: f64,
    /// Maximum read quality, inclusive.
    pub max_qual: f64,
    /// Minimum GC fraction, inclusive, when set.
    pub min_gc: Option<f64>,
    /// Maximum GC fraction, inclusive, when set.
    pub max_gc: Option<f64>,
    /// Quality summary used for the quality bounds.
    pub qual_mode: QualMode,
}

impl Default for FilterConfig {
    /// No bounds: every length, quality and GC content passes.
    fn default() -> Self {
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
}

/// Returns the fraction of `G`/`C` bases (either case) in `seq`; `0.0` for an
/// empty slice.
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

/// The reason `check` dropped a segment. Both GC bounds collapse into `Gc`:
/// the summary reports "GC out of range", not which side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    /// Shorter than `min_length`, or empty.
    TooShort,
    /// Longer than `max_length`.
    TooLong,
    /// Quality below `min_qual`.
    LowQuality,
    /// Quality above `max_qual`.
    HighQuality,
    /// GC fraction outside `[min_gc, max_gc]`.
    Gc,
}

impl DropReason {
    /// Returns the wording used for this reason wherever it is reported, so the
    /// end-of-run summary and a per-segment trace line name it identically.
    pub fn label(self) -> &'static str {
        match self {
            DropReason::TooShort => "too short",
            DropReason::TooLong => "too long",
            DropReason::LowQuality => "low quality",
            DropReason::HighQuality => "high quality",
            DropReason::Gc => "GC out of range",
        }
    }
}

/// Evaluates the bounds cheapest-first and stops at the first rejection. The
/// result names the bound the segment fails, or is `None` when the segment
/// passes; empty segments are `TooShort` even when `min_length` is zero.
///
/// Called once for each segment produced by trimming, so `seq` and `phred`
/// describe that segment rather than necessarily the complete input read.
pub fn check(seq: &[u8], phred: &[u8], cfg: &FilterConfig) -> Option<DropReason> {
    check_with_gc(seq.len(), phred, || gc_fraction(seq), cfg)
}

/// Evaluates sequence composition only after the length and quality bounds pass.
pub(crate) fn check_with_gc(
    len: usize,
    phred: &[u8],
    gc: impl FnOnce() -> f64,
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
        // Probability summation and logarithms introduce rounding at inclusive bounds.
        let tolerance = if cfg.qual_mode == QualMode::Mean {
            1e-10
        } else {
            0.0
        };
        if q < cfg.min_qual - tolerance {
            return Some(DropReason::LowQuality);
        }
        if q > cfg.max_qual + tolerance {
            return Some(DropReason::HighQuality);
        }
    }
    if cfg.min_gc.is_some() || cfg.max_gc.is_some() {
        let gc = gc();
        if gc < cfg.min_gc.unwrap_or(0.0) || gc > cfg.max_gc.unwrap_or(1.0) {
            return Some(DropReason::Gc);
        }
    }
    None
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

    fn passes(seq: &[u8], phred: &[u8], cfg: &FilterConfig) -> bool {
        check(seq, phred, cfg).is_none()
    }

    #[test]
    fn probability_quality_bounds_include_equal_scores() {
        for q in [0, 7, 10, 12, 20, 30, 40, 93] {
            for len in [1, 3, 30, 100, 511, 512, 513, 1000] {
                let seq = vec![b'C'; len];
                let phred = vec![q; len];
                let cfg = FilterConfig {
                    min_qual: f64::from(q),
                    max_qual: f64::from(q),
                    ..base()
                };
                assert_eq!(check(&seq, &phred, &cfg), None, "Q{q}, {len} bases");
                let above = FilterConfig {
                    min_qual: f64::from(q) + 1e-6,
                    max_qual: 1000.0,
                    ..base()
                };
                assert_eq!(check(&seq, &phred, &above), Some(DropReason::LowQuality));
            }
        }
    }

    #[test]
    fn gc_is_evaluated_only_after_earlier_filters_pass() {
        let cfg = FilterConfig {
            min_length: 10,
            min_qual: 20.0,
            min_gc: Some(0.5),
            ..base()
        };
        assert_eq!(
            check_with_gc(4, &[30; 4], || panic!("GC evaluated"), &cfg),
            Some(DropReason::TooShort)
        );
        assert_eq!(
            check_with_gc(12, &[10; 12], || panic!("GC evaluated"), &cfg),
            Some(DropReason::LowQuality)
        );
        assert_eq!(
            check_with_gc(12, &[30; 12], || 0.4, &cfg),
            Some(DropReason::Gc)
        );
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
        // The arithmetic mean of [10, 20] is 15.0, which passes at the threshold.
        c.qual_mode = QualMode::Arithmetic;
        assert!(passes(b"AT", &[10, 20], &c));
        // The probability mean of [10, 20] is below 15, which fails.
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
    fn gc_is_not_evaluated_without_a_bound() {
        assert_eq!(
            check_with_gc(4, &[30; 4], || panic!("GC evaluated"), &base()),
            None
        );
    }

    /// The bounds are evaluated cheapest-first, so a segment failing both
    /// length and quality reports the length verdict.
    #[test]
    fn length_is_evaluated_before_quality() {
        let mut c = base();
        c.qual_mode = QualMode::Arithmetic;
        c.min_qual = 30.0;
        c.min_length = 3;
        assert_eq!(check(b"AT", &[10, 20], &c), Some(DropReason::TooShort));
    }

    #[test]
    fn check_reports_too_short() {
        let mut c = base();
        c.min_length = 4;
        assert_eq!(check(b"ATG", &[30, 30, 30], &c), Some(DropReason::TooShort));
        // Empty reads are `TooShort` regardless of `min_length`.
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
}
