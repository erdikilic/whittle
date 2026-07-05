pub mod strategies;

use strategies::{best_segment, split_low_quality, trim_by_quality};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QualityOp {
    TrimQual(u8),
    BestSegment(u8),
    Split { cutoff: u8, window: usize },
}

#[derive(Debug, Clone)]
pub struct TrimPlan {
    pub head: usize,
    pub tail: usize,
    pub quality: Option<QualityOp>,
}

/// Fixed crop first (positional), then the chosen quality op on the cropped
/// window, offsetting intervals back to original coordinates. Every returned
/// segment is >= `min_length`.
pub fn apply(
    seq_len: usize,
    phred: &[u8],
    plan: &TrimPlan,
    min_length: usize,
) -> Vec<(usize, usize)> {
    debug_assert_eq!(
        seq_len,
        phred.len(),
        "apply: seq_len must equal phred.len()"
    );
    let start = plan.head.min(seq_len);
    let end = seq_len.saturating_sub(plan.tail).max(start);
    if start >= end {
        return vec![];
    }
    let window_phred = &phred[start..end];
    let inner = match &plan.quality {
        None => vec![(0, window_phred.len())],
        Some(QualityOp::TrimQual(q)) => trim_by_quality(window_phred, *q),
        Some(QualityOp::BestSegment(q)) => best_segment(window_phred, *q),
        Some(QualityOp::Split { cutoff, window }) => {
            split_low_quality(window_phred, *cutoff, min_length, *window)
        }
    };
    inner
        .into_iter()
        .map(|(s, e)| (s + start, e + start))
        .filter(|&(s, e)| e - s >= min_length)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_quality_op_is_fixed_crop() {
        let phred = vec![30u8; 20];
        let plan = TrimPlan {
            head: 5,
            tail: 3,
            quality: None,
        };
        assert_eq!(apply(20, &phred, &plan, 1), vec![(5, 17)]);
    }

    #[test]
    fn crop_then_quality_offsets_back() {
        // 20 bases; crop head 2; then trim_qual on the remaining window.
        // phred: low low [good...] -> after crop the good region starts at 2.
        let mut phred = vec![40u8; 20];
        phred[0] = 2;
        phred[1] = 2;
        let plan = TrimPlan {
            head: 2,
            tail: 0,
            quality: Some(QualityOp::TrimQual(30)),
        };
        assert_eq!(apply(20, &phred, &plan, 1), vec![(2, 20)]);
    }

    #[test]
    fn min_length_drops_short_segments() {
        let phred = vec![40u8; 4];
        let plan = TrimPlan {
            head: 0,
            tail: 0,
            quality: None,
        };
        assert_eq!(apply(4, &phred, &plan, 5), Vec::<(usize, usize)>::new());
    }

    #[test]
    fn empty_when_crop_exceeds_length() {
        let phred = vec![40u8; 4];
        let plan = TrimPlan {
            head: 3,
            tail: 3,
            quality: None,
        };
        assert_eq!(apply(4, &phred, &plan, 1), Vec::<(usize, usize)>::new());
    }
}
