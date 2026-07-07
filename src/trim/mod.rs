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

/// Fixed crop first (positional), then the adapter stage on the cropped window
/// (when configured), then the chosen quality op within each adapter segment,
/// offsetting intervals back to original coordinates. Emits crop->adapter->quality
/// segments only; the caller filters (length/quality/GC) per segment. `min_length`
/// is not applied here as a filter — it is only forwarded to `split_low_quality`,
/// where it is the minimum quality-split piece size (a split-granularity control,
/// not a biological length filter).
pub fn apply(
    seq: &[u8],
    phred: &[u8],
    plan: &TrimPlan,
    adapters: Option<&crate::adapter::AdapterConfig>,
    min_length: usize,
) -> Vec<(usize, usize)> {
    debug_assert_eq!(
        seq.len(),
        phred.len(),
        "apply: seq.len() must equal phred.len()"
    );
    let seq_len = seq.len();
    let start = plan.head.min(seq_len);
    let end = seq_len.saturating_sub(plan.tail).max(start);
    if start >= end {
        return vec![];
    }

    // Adapter stage on the cropped window, mapped back to original coordinates.
    let adapter_segs: Vec<(usize, usize)> = match adapters {
        Some(cfg) => crate::adapter::adapter_segments(&seq[start..end], cfg)
            .into_iter()
            .map(|(s, e)| (s + start, e + start))
            .collect(),
        None => vec![(start, end)],
    };

    // Quality op within each adapter segment, offset back. No length filter here
    // — the caller filters each returned segment (length/quality/GC).
    let mut out = Vec::new();
    for (s, e) in adapter_segs {
        let window_phred = &phred[s..e];
        let inner = match &plan.quality {
            None => vec![(0, window_phred.len())],
            Some(QualityOp::TrimQual(q)) => trim_by_quality(window_phred, *q),
            Some(QualityOp::BestSegment(q)) => best_segment(window_phred, *q),
            Some(QualityOp::Split { cutoff, window }) => {
                split_low_quality(window_phred, *cutoff, min_length, *window)
            },
        };
        out.extend(inner.into_iter().map(|(is, ie)| (is + s, ie + s)));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_quality_op_is_fixed_crop() {
        let phred = vec![30u8; 20];
        let seq = vec![b'A'; 20];
        let plan = TrimPlan {
            head: 5,
            tail: 3,
            quality: None,
        };
        assert_eq!(apply(&seq, &phred, &plan, None, 1), vec![(5, 17)]);
    }

    #[test]
    fn crop_then_quality_offsets_back() {
        // 20 bases; crop head 2; then trim_qual on the remaining window.
        // phred: low low [good...] -> after crop the good region starts at 2.
        let mut phred = vec![40u8; 20];
        phred[0] = 2;
        phred[1] = 2;
        let seq = vec![b'A'; 20];
        let plan = TrimPlan {
            head: 2,
            tail: 0,
            quality: Some(QualityOp::TrimQual(30)),
        };
        assert_eq!(apply(&seq, &phred, &plan, None, 1), vec![(2, 20)]);
    }

    #[test]
    fn short_segments_are_emitted_not_filtered() {
        // `apply` no longer applies a `min_length` filter: filtering moved to the
        // caller (per-segment, post-trim). A segment shorter than `min_length` is
        // still RETURNED here.
        let phred = vec![40u8; 4];
        let seq = vec![b'A'; 4];
        let plan = TrimPlan {
            head: 0,
            tail: 0,
            quality: None,
        };
        assert_eq!(apply(&seq, &phred, &plan, None, 5), vec![(0, 4)]);
    }

    #[test]
    fn empty_when_crop_exceeds_length() {
        let phred = vec![40u8; 4];
        let seq = vec![b'A'; 4];
        let plan = TrimPlan {
            head: 3,
            tail: 3,
            quality: None,
        };
        assert_eq!(
            apply(&seq, &phred, &plan, None, 1),
            Vec::<(usize, usize)>::new()
        );
    }

    #[test]
    fn adapter_stage_runs_before_quality_op() {
        use crate::adapter::{Adapter, AdapterConfig, End};
        let adapter = b"ACGTACGTACGT";
        let mut seq = adapter.to_vec();
        seq.extend_from_slice(b"GGGGGGGGGGGG");
        let phred = vec![40u8; seq.len()];
        let plan = TrimPlan {
            head: 0,
            tail: 0,
            quality: None,
        };
        let ac = AdapterConfig {
            adapters: vec![Adapter {
                name: "a".into(),
                seq: adapter.to_vec(),
                end: End::Five,
            }],
            error_rate: 0.2,
            end_size: 20,
            split: false,
        };
        assert_eq!(apply(&seq, &phred, &plan, Some(&ac), 1), vec![(12, 24)]);
    }

    #[test]
    fn no_adapter_config_matches_old_behavior() {
        let phred = vec![30u8; 20];
        let seq = vec![b'A'; 20];
        let plan = TrimPlan {
            head: 5,
            tail: 3,
            quality: None,
        };
        assert_eq!(apply(&seq, &phred, &plan, None, 1), vec![(5, 17)]);
    }
}
