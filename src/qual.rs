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

/// Length from which `mean_prob_q` sums through the quality histogram rather
/// than per base. The histogram's fixed cost (zeroing the bins and the
/// 256-term dot product) exceeds its per-base saving below this length.
const HISTOGRAM_MIN_LEN: usize = 512;

/// Returns the error-probability mean quality: the ONT read Q, the mean
/// per-base error probability converted back to a Phred score. Long reads sum
/// through `histogram_prob_sum`; shorter ones sum per base.
pub fn mean_prob_q(phred: &[u8]) -> f64 {
    if phred.is_empty() {
        return 0.0;
    }
    let sum: f64 = if phred.len() < HISTOGRAM_MIN_LEN {
        phred.iter().map(|&q| phred_to_prob(q)).sum()
    } else {
        histogram_prob_sum(phred)
    };
    (sum / phred.len() as f64).log10() * -10.0
}

/// Returns the sum of the per-base error probabilities of `phred` as the dot
/// product of a quality histogram with `PHRED_LUT`, so the per-base work is an
/// integer increment rather than a table lookup and a dependent float add.
/// The bytes are counted into four interleaved histograms, so a run of equal
/// qualities does not serialize on one counter, and the bins are merged once.
/// The histogram spans the `u8` range: FASTQ input is bounded to Phred 0-93,
/// BAM quality bytes are not. A `u32` bin holds 2^32 bases of one quality per
/// interleave, beyond the length of any record.
fn histogram_prob_sum(phred: &[u8]) -> f64 {
    let mut hist = [[0u32; 256]; 4];
    let mut chunks = phred.chunks_exact(4);
    for chunk in &mut chunks {
        hist[0][usize::from(chunk[0])] += 1;
        hist[1][usize::from(chunk[1])] += 1;
        hist[2][usize::from(chunk[2])] += 1;
        hist[3][usize::from(chunk[3])] += 1;
    }
    for &q in chunks.remainder() {
        hist[0][usize::from(q)] += 1;
    }
    let [h0, h1, h2, h3] = &hist;
    (0..256)
        .map(|q| {
            let count = u64::from(h0[q]) + u64::from(h1[q]) + u64::from(h2[q]) + u64::from(h3[q]);
            count as f64 * PHRED_LUT[q]
        })
        .sum()
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

    /// The histogram form agrees with the per-base sum to well within the
    /// rounding of either, across mixed and uniform quality strings and
    /// lengths that are not a multiple of the interleave; below the histogram
    /// length the two are the same computation.
    #[test]
    fn mean_prob_q_matches_the_per_base_sum() {
        let per_base = |phred: &[u8]| -> f64 {
            let sum: f64 = phred.iter().map(|&q| phred_to_prob(q)).sum();
            (sum / phred.len() as f64).log10() * -10.0
        };
        let mixed: Vec<u8> = (0..20_001u32).map(|i| ((i * 7919) % 94) as u8).collect();
        assert!((mean_prob_q(&mixed) - per_base(&mixed)).abs() < 1e-12);
        let uniform = vec![12u8; 5_003];
        assert!((mean_prob_q(&uniform) - per_base(&uniform)).abs() < 1e-12);
        assert!((mean_prob_q(&uniform) - 12.0).abs() < 1e-12);
        let short = &mixed[..HISTOGRAM_MIN_LEN - 1];
        assert_eq!(mean_prob_q(short), per_base(short));
        let at_threshold = &mixed[..HISTOGRAM_MIN_LEN];
        assert!((mean_prob_q(at_threshold) - per_base(at_threshold)).abs() < 1e-12);
        // A BAM quality byte above the FASTQ range indexes the table safely,
        // on both paths.
        assert!(mean_prob_q(&[200, 10]).is_finite());
        assert!(mean_prob_q(&[255u8; 1000]).is_finite());
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
