use std::sync::LazyLock;

/// Precomputed 10^(-q/10) for every possible Phred byte. Sizing to the full u8
/// range means any quality byte indexes safely (ported from chopper's PHRED_LUT).
static PHRED_LUT: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut lut = [0.0f64; 256];
    for (i, v) in lut.iter_mut().enumerate() {
        *v = 10_f64.powf((i as f64) / -10.0);
    }
    lut
});

#[inline(always)]
pub fn phred_to_prob(q: u8) -> f64 {
    PHRED_LUT[q as usize]
}

/// Error-probability mean quality: the ONT-standard "read Q" (chopper's ave_qual).
pub fn mean_prob_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let sum: f64 = phred.iter().map(|&q| phred_to_prob(q)).sum();
    (sum / phred.len() as f64).log10() * -10.0
}

/// Plain arithmetic mean of the Phred integers.
pub fn mean_arith_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let sum: u64 = phred.iter().map(|&q| q as u64).sum();
    sum as f64 / phred.len() as f64
}

/// Median Phred via a 256-bucket histogram (O(n), no sort/alloc of the input).
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
    // Walk buckets accumulating counts; find value(s) at the median rank(s).
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QualMode {
    Mean,
    Arithmetic,
    Median,
}

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

    // Values ported from chopper's ave_qual test, but inputs are RAW phred (no +33).
    #[test]
    fn mean_prob_q_matches_chopper() {
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
        // even count -> average of the two middle values
        assert!((median_q(&[10, 20, 30, 40]) - 25.0).abs() < 1e-9);
    }

    #[test]
    fn read_quality_dispatches() {
        assert_eq!(read_quality(&[10, 20, 30], QualMode::Arithmetic), 20.0);
        assert_eq!(read_quality(&[10, 20, 30], QualMode::Median), 20.0);
        assert_eq!(read_quality(&[10], QualMode::Mean), mean_prob_q(&[10]));
    }
}
