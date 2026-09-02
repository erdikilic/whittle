//! Per-read trimming: fixed crop, adapter stage and quality stage, producing kept intervals.

pub mod strategies;

use strategies::{best_segment, split_low_quality, trim_by_quality};

/// The quality-based operation applied within each adapter segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QualityOp {
    /// Trimming of both ends up to the first base at or above the cutoff.
    TrimQual(u8),
    /// The single highest-scoring segment (modified Mott).
    BestSegment(u8),
    /// A split at runs of at least `window` bases below `cutoff`.
    Split {
        /// Phred cutoff below which a base counts as low quality.
        cutoff: u8,
        /// Minimum run of low-quality bases that splits the read.
        window: usize,
    },
}

/// The per-read trim configuration.
#[derive(Debug, Clone)]
pub struct TrimPlan {
    /// Bases removed from the 5' end.
    pub head: usize,
    /// Bases removed from the 3' end.
    pub tail: usize,
    /// Quality operation applied after the crop and adapter stages, if any.
    pub quality: Option<QualityOp>,
}

/// Applies `plan` to one read and returns the kept intervals in read
/// coordinates: the fixed crop first, then the adapter stage on the cropped
/// window (when configured), then the quality operation within each adapter
/// segment. Every segment is returned, including short ones; the caller filters
/// each by length, quality and GC.
pub fn apply(
    seq: &[u8],
    phred: &[u8],
    plan: &TrimPlan,
    adapters: Option<&crate::adapter::AdapterConfig>,
) -> Vec<(usize, usize)> {
    debug_assert_eq!(
        seq.len(),
        phred.len(),
        "Sequence and quality lengths must be equal"
    );
    let seq_len = seq.len();
    let start = plan.head.min(seq_len);
    let end = seq_len.saturating_sub(plan.tail).max(start);
    if start >= end {
        return vec![];
    }

    // The quality op within one `[s, e)` segment, with results offset back to
    // read coordinates and appended to `out`. No length filter is applied; the
    // caller filters each returned segment.
    let quality_in = |s: usize, e: usize, out: &mut Vec<(usize, usize)>| {
        let wp = &phred[s..e];
        let offset = |v: Vec<(usize, usize)>, out: &mut Vec<(usize, usize)>| {
            out.extend(v.into_iter().map(|(is, ie)| (is + s, ie + s)));
        };
        match &plan.quality {
            None => out.push((s, e)),
            Some(QualityOp::TrimQual(q)) => offset(trim_by_quality(wp, *q), out),
            Some(QualityOp::BestSegment(q)) => offset(best_segment(wp, *q), out),
            Some(QualityOp::Split { cutoff, window }) => {
                offset(split_low_quality(wp, *cutoff, *window), out)
            },
        }
    };

    // The adapter stage on the cropped window, mapped back to read coordinates,
    // then the quality op within each segment. The no-adapter path goes directly
    // to the quality op with no intermediate segment vector.
    let mut out = Vec::new();
    match adapters {
        None => {
            quality_in(start, end, &mut out);
        },
        Some(cfg) => {
            for (s, e) in crate::adapter::adapter_segments(&seq[start..end], cfg) {
                quality_in(s + start, e + start, &mut out);
            }
        },
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
        assert_eq!(apply(&seq, &phred, &plan, None), vec![(5, 17)]);
    }

    #[test]
    fn crop_then_quality_offsets_back() {
        // 20 bases, head crop 2, then `TrimQual` on the remaining window. The
        // first two Phred values are low, so the good region starts at 2 after
        // the crop.
        let mut phred = vec![40u8; 20];
        phred[0] = 2;
        phred[1] = 2;
        let seq = vec![b'A'; 20];
        let plan = TrimPlan {
            head: 2,
            tail: 0,
            quality: Some(QualityOp::TrimQual(30)),
        };
        assert_eq!(apply(&seq, &phred, &plan, None), vec![(2, 20)]);
    }

    /// `apply` applies no length filter; the caller filters per segment after
    /// trimming, so a short segment is returned.
    #[test]
    fn short_segments_are_emitted_not_filtered() {
        let phred = vec![40u8; 4];
        let seq = vec![b'A'; 4];
        let plan = TrimPlan {
            head: 0,
            tail: 0,
            quality: None,
        };
        assert_eq!(apply(&seq, &phred, &plan, None), vec![(0, 4)]);
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
            apply(&seq, &phred, &plan, None),
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
            candidate_index: std::sync::OnceLock::new(),
        };
        assert_eq!(apply(&seq, &phred, &plan, Some(&ac)), vec![(12, 24)]);
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
        assert_eq!(apply(&seq, &phred, &plan, None), vec![(5, 17)]);
    }
}
