//! Adapter trimming and splitting, barcode restriction, per-segment cropping
//! and quality processing, expressed as intervals in the original read.

pub mod strategies;

use strategies::{best_segment, maximal_segments, split_low_quality, trim_by_quality};

/// The quality trimming method, selected by `--quality-trim`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum QualityMethod {
    /// Trim each end up to the first base at or above the cutoff
    Ends,
    /// Keep the single best segment under the modified Mott score
    Best,
    /// Keep every maximal segment under the modified Mott score
    Segments,
    /// Split at runs of --min-low-quality-run bases below the cutoff
    Runs,
}

impl QualityMethod {
    /// Lowercase label, as the command line, the banner and the summary spell
    /// it.
    pub fn label(self) -> &'static str {
        match self {
            QualityMethod::Ends => "ends",
            QualityMethod::Best => "best",
            QualityMethod::Segments => "segments",
            QualityMethod::Runs => "runs",
        }
    }
}

/// The quality trimming applied within each adapter segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QualityOp {
    /// The trimming method.
    pub method: QualityMethod,
    /// Phred cutoff the method scores or compares each base against.
    pub cutoff: u8,
    /// Minimum run of bases below `cutoff` that splits a read under
    /// `QualityMethod::Runs`; 1 for the other methods.
    pub min_low_quality_run: usize,
}

impl QualityOp {
    /// The operation `method` at `cutoff`, splitting at every low-quality base
    /// under `QualityMethod::Runs`.
    pub fn new(method: QualityMethod, cutoff: u8) -> Self {
        QualityOp {
            method,
            cutoff,
            min_low_quality_run: 1,
        }
    }

    /// The `QualityMethod::Runs` operation at `cutoff`, splitting at runs of at
    /// least `min_low_quality_run` bases below it.
    pub fn runs(cutoff: u8, min_low_quality_run: usize) -> Self {
        QualityOp {
            min_low_quality_run,
            ..QualityOp::new(QualityMethod::Runs, cutoff)
        }
    }

    /// Returns the intervals of `phred` the operation keeps, in its coordinates.
    fn apply(&self, phred: &[u8]) -> Vec<(usize, usize)> {
        match self.method {
            QualityMethod::Ends => trim_by_quality(phred, self.cutoff),
            QualityMethod::Best => best_segment(phred, self.cutoff),
            QualityMethod::Segments => maximal_segments(phred, self.cutoff),
            QualityMethod::Runs => split_low_quality(phred, self.cutoff, self.min_low_quality_run),
        }
    }
}

/// The per-read trim configuration.
#[derive(Debug, Clone, Default)]
pub struct TrimPlan {
    /// Bases removed from each adapter-derived segment's 5' end after barcode restriction.
    pub head: usize,
    /// Bases removed from each adapter-derived segment's 3' end after barcode restriction.
    pub tail: usize,
    /// Quality operation applied to each cropped segment.
    pub quality: Option<QualityOp>,
}

/// Searches the original read for adapters, intersects each retained segment
/// with `barcode`, crops each intersection once, and applies the quality
/// operation. Quality splitting can produce multiple intervals per adapter
/// segment. The result is flattened in original-coordinate order for final
/// numbering and filtering by length, quality and GC.
///
/// `barcode` is the retained interval resolved from the original record's
/// verified `bi` spans. `None` retains the whole read. An unmatched
/// read enters barcode restriction and cropping as one full-length segment.
pub fn apply(
    seq: &[u8],
    phred: &[u8],
    plan: &TrimPlan,
    adapters: Option<&crate::adapter::AdapterConfig>,
    barcode: Option<(usize, usize)>,
) -> Vec<(usize, usize)> {
    debug_assert_eq!(
        seq.len(),
        phred.len(),
        "Sequence and quality lengths must be equal"
    );
    let seq_len = seq.len();
    let (outer_start, outer_end) = match barcode {
        Some((s, e)) => (s.min(seq_len), e.clamp(s.min(seq_len), seq_len)),
        None => (0, seq_len),
    };
    if outer_start >= outer_end {
        return vec![];
    }

    let process_segment = |s: usize, e: usize, out: &mut Vec<(usize, usize)>| {
        let s = s.max(outer_start);
        let e = e.min(outer_end);
        if s >= e {
            return;
        }
        let s = s.saturating_add(plan.head).min(e);
        let e = e.saturating_sub(plan.tail).max(s);
        if s >= e {
            return;
        }
        let wp = &phred[s..e];
        let offset = |v: Vec<(usize, usize)>, out: &mut Vec<(usize, usize)>| {
            out.extend(v.into_iter().map(|(is, ie)| (is + s, ie + s)));
        };
        match &plan.quality {
            None => out.push((s, e)),
            Some(op) => offset(op.apply(wp), out),
        }
    };

    let mut out = Vec::new();
    match adapters {
        None => {
            process_segment(0, seq_len, &mut out);
        },
        Some(cfg) => {
            for (s, e) in crate::adapter::adapter_segments(seq, cfg) {
                process_segment(s, e, &mut out);
            }
        },
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn adapter_config() -> crate::adapter::AdapterConfig {
        crate::adapter::AdapterConfig {
            adapters: vec![crate::adapter::Adapter {
                name: "junction".into(),
                seq: b"GGGGTTTTGGGGTTTT".to_vec(),
                role: crate::adapter::Role::Adapter,
            }],
            error_rate: 0.0,
            end_size: 8,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
        }
    }

    #[test]
    fn adapter_segments_are_cropped_once_before_quality_processing() {
        let ac = adapter_config();
        let seq = [vec![b'C'; 64], ac.adapters[0].seq.clone(), vec![b'A'; 64]].concat();
        let mut phred = vec![40; seq.len()];
        phred[20..24].fill(2);
        phred[104..108].fill(2);
        let plan = TrimPlan {
            head: 3,
            tail: 5,
            quality: Some(QualityOp::runs(9, 4)),
        };
        assert_eq!(
            apply(&seq, &phred, &plan, Some(&ac), None),
            vec![(3, 20), (24, 59), (83, 104), (108, 139)]
        );
        assert_eq!(
            apply(&seq, &phred, &plan, Some(&ac), Some((10, 135))),
            vec![(13, 20), (24, 59), (83, 104), (108, 130)]
        );
        let end_trim = TrimPlan {
            quality: Some(QualityOp::new(QualityMethod::Ends, 20)),
            ..plan.clone()
        };
        assert_eq!(
            apply(&seq, &phred, &end_trim, Some(&ac), None),
            vec![(3, 59), (83, 139)]
        );
        let best = TrimPlan {
            quality: Some(QualityOp::new(QualityMethod::Best, 20)),
            ..plan
        };
        assert_eq!(
            apply(&seq, &phred, &best, Some(&ac), None),
            vec![(24, 59), (108, 139)]
        );
    }

    /// Each maximal segment of a cropped segment is offset back to read
    /// coordinates.
    #[test]
    fn split_segments_are_offset_to_read_coordinates() {
        let seq = vec![b'C'; 200];
        let mut phred = vec![40; seq.len()];
        phred[80..100].fill(2);
        let plan = TrimPlan {
            head: 3,
            tail: 5,
            quality: Some(QualityOp::new(QualityMethod::Segments, 20)),
        };
        assert_eq!(
            apply(&seq, &phred, &plan, None, None),
            vec![(3, 80), (100, 195)]
        );
    }

    #[test]
    fn terminal_adapter_is_recognized_before_barcode_restriction_and_crop() {
        let ac = adapter_config();
        let seq = [ac.adapters[0].seq.clone(), vec![b'C'; 64]].concat();
        let phred = vec![40; seq.len()];
        let plan = TrimPlan {
            head: 3,
            tail: 5,
            quality: None,
        };
        assert_eq!(apply(&seq, &phred, &plan, Some(&ac), None), vec![(19, 75)]);
        assert_eq!(
            apply(&seq, &phred, &plan, Some(&ac), Some((12, 70))),
            vec![(19, 65)]
        );
    }

    #[test]
    fn unmatched_read_is_cropped_and_quality_split() {
        let ac = adapter_config();
        let seq = vec![b'C'; 64];
        let mut phred = vec![40; seq.len()];
        phred[20..24].fill(2);
        let plan = TrimPlan {
            head: 3,
            tail: 5,
            quality: Some(QualityOp::runs(9, 4)),
        };
        assert_eq!(
            apply(&seq, &phred, &plan, Some(&ac), None),
            vec![(3, 20), (24, 59)]
        );
    }

    #[test]
    fn crop_can_consume_one_adapter_segment_without_consuming_its_sibling() {
        let ac = adapter_config();
        let seq = [vec![b'C'; 24], ac.adapters[0].seq.clone(), vec![b'A'; 64]].concat();
        let phred = vec![40; seq.len()];
        let plan = TrimPlan {
            head: 20,
            tail: 5,
            quality: None,
        };
        assert_eq!(apply(&seq, &phred, &plan, Some(&ac), None), vec![(60, 99)]);
        assert!(apply(&seq, &phred, &plan, Some(&ac), Some((0, 24))).is_empty());
    }

    #[test]
    fn no_quality_op_is_fixed_crop() {
        let phred = vec![30u8; 20];
        let seq = vec![b'A'; 20];
        let plan = TrimPlan {
            head: 5,
            tail: 3,
            quality: None,
        };
        assert_eq!(apply(&seq, &phred, &plan, None, None), vec![(5, 17)]);
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
            quality: Some(QualityOp::new(QualityMethod::Ends, 30)),
        };
        assert_eq!(apply(&seq, &phred, &plan, None, None), vec![(2, 20)]);
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
        assert_eq!(apply(&seq, &phred, &plan, None, None), vec![(0, 4)]);
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
            apply(&seq, &phred, &plan, None, None),
            Vec::<(usize, usize)>::new()
        );
    }

    #[test]
    fn adapter_stage_runs_before_quality_op() {
        use crate::adapter::{Adapter, AdapterConfig, Role};
        let adapter = b"ACGTACGTACGT";
        let mut seq = adapter.to_vec();
        seq.extend_from_slice(b"GGGGGGGGGGGG");
        let mut phred = vec![40u8; seq.len()];
        phred[12..15].fill(2);
        let plan = TrimPlan {
            head: 2,
            tail: 1,
            quality: Some(QualityOp::new(QualityMethod::Ends, 20)),
        };
        let ac = AdapterConfig {
            adapters: vec![Adapter {
                name: "a".into(),
                seq: adapter.to_vec(),
                role: Role::Adapter,
            }],
            error_rate: 0.2,
            end_size: 20,
            split: false,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
        };
        assert_eq!(apply(&seq, &phred, &plan, Some(&ac), None), vec![(15, 23)]);
    }

    #[test]
    fn crop_without_adapters_uses_the_full_read() {
        let phred = vec![30u8; 20];
        let seq = vec![b'A'; 20];
        let plan = TrimPlan {
            head: 5,
            tail: 3,
            quality: None,
        };
        assert_eq!(apply(&seq, &phred, &plan, None, None), vec![(5, 17)]);
    }
}

#[cfg(test)]
mod barcode_tests {
    use super::*;

    /// Fixed cropping operates within the original-coordinate barcode interval.
    #[test]
    fn barcode_window_precedes_the_crop() {
        let seq = vec![b'A'; 20];
        let phred = vec![40u8; 20];
        let plan = TrimPlan {
            head: 2,
            tail: 1,
            quality: None,
        };
        assert_eq!(
            apply(&seq, &phred, &plan, None, Some((5, 15))),
            vec![(7, 14)]
        );
    }

    /// A crop wider than the barcode window keeps nothing.
    #[test]
    fn crop_beyond_the_barcode_window_keeps_nothing() {
        let seq = vec![b'A'; 20];
        let phred = vec![40u8; 20];
        let plan = TrimPlan {
            head: 6,
            tail: 6,
            quality: None,
        };
        let kept = apply(&seq, &phred, &plan, None, Some((5, 15)));
        assert!(kept.is_empty());
    }
}
