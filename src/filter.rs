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

/// Evaluates the bounds cheapest-first and stops at the first rejection. `None`
/// indicates that the segment passes; empty segments are `TooShort` even when
/// `min_length` is zero.
///
/// Called once for each segment produced by trimming, so `seq` and `phred`
/// describe that segment rather than necessarily the complete input read.
pub fn check(seq: &[u8], phred: &[u8], cfg: &FilterConfig) -> Option<DropReason> {
    let gc = (cfg.min_gc.is_some() || cfg.max_gc.is_some()).then(|| gc_fraction(seq));
    check_metrics(seq.len(), phred, gc, cfg)
}

/// Filters from precomputed sequence metrics. The raw BAM fast path obtains
/// length and GC from packed sequence views without materializing a decoded
/// sequence. `gc` is required whenever a GC bound is active; a missing value
/// panics rather than filtering on a fabricated fraction.
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
        let gc = gc.expect("The GC fraction is computed whenever a GC bound is active");
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
    fn check_metrics_uses_the_supplied_gc_when_a_bound_is_active() {
        let mut c = base();
        c.min_gc = Some(0.4);
        assert_eq!(check_metrics(4, &[30; 4], Some(0.5), &c), None);
        assert_eq!(
            check_metrics(4, &[30; 4], Some(0.1), &c),
            Some(DropReason::Gc)
        );
        // No bound active: a missing GC value is not consulted.
        assert_eq!(check_metrics(4, &[30; 4], None, &base()), None);
    }

    #[test]
    #[should_panic(expected = "GC")]
    fn check_metrics_panics_without_gc_when_a_bound_is_active() {
        let mut c = base();
        c.max_gc = Some(0.6);
        check_metrics(4, &[30; 4], None, &c);
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
