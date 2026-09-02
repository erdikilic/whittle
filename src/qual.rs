//! Phred-scale quality summaries for a read: error-probability mean, arithmetic mean and median.

use std::sync::LazyLock;

/// Precomputed `10^(-q/10)` for every Phred byte. The table covers the full
/// `u8` range, so any quality byte indexes it safely.
static PHRED_LUT: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut lut = [0.0f64; 256];
    for (i, v) in lut.iter_mut().enumerate() {
        *v = 10_f64.powf((i as f64) / -10.0);
    }
    lut
});

/// Returns the error probability `10^(-q/10)` for a raw Phred score, read from
/// `PHRED_LUT`.
#[inline(always)]
pub fn phred_to_prob(q: u8) -> f64 {
    PHRED_LUT[q as usize]
}

/// Returns the error-probability mean quality: the ONT read Q, the mean
/// per-base error probability converted back to a Phred score.
pub fn mean_prob_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let sum: f64 = phred.iter().map(|&q| phred_to_prob(q)).sum();
    (sum / phred.len() as f64).log10() * -10.0
}

/// Returns the arithmetic mean of the Phred integers.
pub fn mean_arith_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let sum: u64 = phred.iter().map(|&q| q as u64).sum();
    sum as f64 / phred.len() as f64
}

/// Returns the median Phred score via a 256-bucket histogram, in O(n) without
/// sorting or allocating.
pub fn median_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let mut hist = [0usize; 256];
    for &q in phred {
        hist[q as usize] += 1;
    }
    let n = phred.len();
    let mid = n / 2;
    // The value at rank `target`, found by accumulating bucket counts.
    let value_at = |target: usize| -> usize {
        let mut cum = 0usize;
        for (v, &c) in hist.iter().enumerate() {
            cum += c;
            if cum > target {
                return v;
            }
        }
        255
    };
    if n % 2 == 1 {
        value_at(mid) as f64
    } else {
        (value_at(mid - 1) + value_at(mid)) as f64 / 2.0
    }
}

/// The read-quality summary used by the quality filter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QualMode {
    /// Error-probability mean (the ONT read Q).
    Mean,
    /// Arithmetic mean of the Phred integers.
    Arithmetic,
    /// Median Phred score.
    Median,
}

/// Returns the read quality of `phred` under `mode`; `0.0` for an empty slice.
pub fn read_quality(phred: &[u8], mode: QualMode) -> f64 {
    match mode {
        QualMode::Mean => mean_prob_q(phred),
        QualMode::Arithmetic => mean_arith_q(phred),
        QualMode::Median => median_q(phred),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prob_matches_phred_definition() {
        assert!((phred_to_prob(20) - 0.01).abs() < 1e-12);
        assert!((phred_to_prob(30) - 0.001).abs() < 1e-12);
    }

    /// Known-good values for the error-probability mean; inputs are raw Phred
    /// scores without the +33 offset.
    #[test]
    fn mean_prob_q_known_values() {
        assert!((mean_prob_q(&[10]) - 10.0).abs() < 1e-9);
        assert!((mean_prob_q(&[10, 11, 12]) - 10.923583702678473).abs() < 1e-9);
        assert!((mean_prob_q(&[10, 11, 12, 20, 30, 40, 50]) - 14.408827647036087).abs() < 1e-9);
    }

    #[test]
    fn empty_is_zero() {
        assert_eq!(mean_prob_q(&[]), 0.0);
        assert_eq!(mean_arith_q(&[]), 0.0);
        assert_eq!(median_q(&[]), 0.0);
    }

    #[test]
    fn arithmetic_and_median() {
        assert!((mean_arith_q(&[10, 20, 30]) - 20.0).abs() < 1e-9);
        assert!((median_q(&[10, 20, 30]) - 20.0).abs() < 1e-9);
        // An even count averages the two middle values.
        assert!((median_q(&[10, 20, 30, 40]) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn read_quality_dispatches() {
        assert_eq!(read_quality(&[10, 20, 30], QualMode::Arithmetic), 20.0);
        assert_eq!(read_quality(&[10, 20, 30], QualMode::Median), 20.0);
        assert_eq!(read_quality(&[10], QualMode::Mean), mean_prob_q(&[10]));
    }
}
