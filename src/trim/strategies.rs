use crate::qual::phred_to_prob;

/// A list of half-open `[start, end)` index ranges into a read.
type Segments = Vec<(usize, usize)>;

/// Trim low-quality bases from both ends until reaching a base with phred >= cutoff.
/// Ported from chopper TrimByQualityStrategy (inputs are raw phred here).
pub fn trim_by_quality(phred: &[u8], cutoff: u8) -> Vec<(usize, usize)> {
    let len = phred.len();
    let mut start = 0;
    while start < len && phred[start] < cutoff {
        start += 1;
    }
    let mut end = len;
    while end > start && phred[end - 1] < cutoff {
        end -= 1;
    }
    if end <= start {
        vec![]
    } else {
        vec![(start, end)]
    }
}

/// Modified Mott: the single segment with the lowest cumulative error probability.
/// Ported from chopper HighestQualityTrimStrategy; `cutoff_q` is converted to a
/// probability cutoff exactly as chopper did at its call site.
pub fn best_segment(phred: &[u8], cutoff_q: u8) -> Vec<(usize, usize)> {
    let cutoff = phred_to_prob(cutoff_q);
    let mut best_start = usize::MAX;
    let mut best_end = usize::MAX;
    let mut best_cumulative_error = 0.0;
    let mut best_length = 0usize;

    let mut current_start = 0usize;
    let mut current_cumulative_error = -1.0;
    for (i, &q) in phred.iter().enumerate() {
        let prob_error = cutoff - phred_to_prob(q);
        if current_cumulative_error < 0.0 {
            current_cumulative_error = 0.0;
            current_start = i;
        }
        current_cumulative_error += prob_error;
        if best_cumulative_error < current_cumulative_error
            || (best_cumulative_error == current_cumulative_error
                && best_length < i - current_start + 1)
        {
            best_start = current_start;
            best_end = i;
            best_cumulative_error = current_cumulative_error;
            best_length = i - current_start + 1;
        }
    }
    if best_start == usize::MAX {
        vec![]
    } else {
        vec![(best_start, best_end + 1)]
    }
}

/// Split into high-quality segments separated by runs of >= `window` low-quality
/// bases; drop segments shorter than `min_length`. Ported from chopper
/// SplitByLowQualityStrategy.
pub fn split_low_quality(
    phred: &[u8],
    cutoff: u8,
    min_length: usize,
    window: usize,
) -> Vec<(usize, usize)> {
    let window = window.max(1);
    let mut segments = Vec::new();
    let mut segment_start: Option<usize> = None;
    let mut last_good: Option<usize> = None;
    let mut bad_run = 0usize;

    let push = |start: usize, end: usize, out: &mut Segments| {
        if end - start >= min_length {
            out.push((start, end));
        }
    };

    for (i, &q) in phred.iter().enumerate() {
        if q >= cutoff {
            if segment_start.is_none() {
                segment_start = Some(i);
            }
            last_good = Some(i);
            bad_run = 0;
        } else {
            bad_run += 1;
            if bad_run >= window {
                if let (Some(s), Some(lg)) = (segment_start, last_good) {
                    push(s, lg + 1, &mut segments);
                }
                segment_start = None;
                last_good = None;
            }
        }
    }
    if let (Some(s), Some(lg)) = (segment_start, last_good) {
        push(s, lg + 1, &mut segments);
    }
    segments
}

#[cfg(test)]
mod tests {
    use super::*;

    // (seq, ascii_qual) — identical bytes to chopper's trimmers.rs get_reads().
    fn reads() -> [(Vec<u8>, Vec<u8>); 6] {
        let raw = [
            (
                b"AAAAAAAAAAAAAAATTTAA".to_vec(),
                b"&#3-G27C:(@G7B55+C4I".to_vec(),
            ),
            (
                b"TTTTTTTTTTTTTTTTTTTT".to_vec(),
                b"77%'24)FAF9@=94'%054".to_vec(),
            ),
            (
                b"AAAAAAAAAAAAAAATTTTA".to_vec(),
                b"'8$-BF2!C;+59->H@91#".to_vec(),
            ),
            (
                b"AAAAAAAAAAAAAAAAAAAA".to_vec(),
                b"%,42$CH*#0+0C6=0,*6/".to_vec(),
            ),
            (
                b"AAAAAAAAAAAAAAAAAAAT".to_vec(),
                b"-------------------J".to_vec(),
            ),
            (
                b"TAAAAAAAAAAAAAAAAAAA".to_vec(),
                b"I-------------------".to_vec(),
            ),
        ];
        raw.map(|(s, q)| (s, q.iter().map(|&b| b - 33).collect()))
    }

    #[test]
    fn trim_by_quality_matches_chopper() {
        let expected: [(u8, Vec<(usize, usize)>); 6] = [
            (20, vec![(4, 20)]),
            (7, vec![(0, 20)]),
            (15, vec![(1, 19)]),
            (40, vec![]),
            (40, vec![(19, 20)]),
            (40, vec![(0, 1)]),
        ];
        for ((cutoff, want), (_, phred)) in expected.iter().zip(reads()) {
            assert_eq!(trim_by_quality(&phred, *cutoff), *want);
        }
    }

    #[test]
    fn best_segment_matches_chopper() {
        // chopper cutoffs were probabilities; the equivalent Q scores are:
        // 0.01=Q20, 0.199..=Q7, 0.0316..=Q15, 0.0001=Q40.
        let expected: [(u8, Vec<(usize, usize)>); 6] = [
            (20, vec![(10, 16)]),
            (7, vec![(0, 20)]),
            (15, vec![(11, 19)]),
            (40, vec![]),
            (40, vec![(19, 20)]),
            (40, vec![(0, 1)]),
        ];
        for ((cutoff_q, want), (_, phred)) in expected.iter().zip(reads()) {
            assert_eq!(best_segment(&phred, *cutoff_q), *want);
        }
    }

    #[test]
    fn split_matches_chopper() {
        // (cutoff, min_length, expected) — from chopper split_by_low_quality_strategy_test, window=1.
        let cases: [(u8, usize, Segments); 6] = [
            (20, 3, vec![(6, 9), (10, 16)]),
            (7, 3, vec![(4, 15), (17, 20)]),
            (15, 3, vec![(4, 7), (14, 19)]),
            (40, 3, vec![]),
            (40, 1, vec![(19, 20)]),
            (40, 1, vec![(0, 1)]),
        ];
        for ((cutoff, min_length, want), (_, phred)) in cases.iter().zip(reads()) {
            assert_eq!(split_low_quality(&phred, *cutoff, *min_length, 1), *want);
        }
    }

    #[test]
    fn split_window_tolerates_short_dips() {
        // chopper window_test: III#IIII###III with Q40=I, Q2=#
        let phred: Vec<u8> = b"III#IIII###III".iter().map(|&b| b - 33).collect();
        assert_eq!(split_low_quality(&phred, 10, 1, 1), vec![
            (0, 3),
            (4, 8),
            (11, 14)
        ]);
        assert_eq!(split_low_quality(&phred, 10, 1, 3), vec![(0, 8), (11, 14)]);
        assert_eq!(split_low_quality(&phred, 10, 1, 4), vec![(0, 14)]);
    }
}
