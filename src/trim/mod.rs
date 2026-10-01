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
    /// Phred cutoff at which each piece of `QualityMethod::Segments` or
    /// `QualityMethod::Runs` has its ends trimmed after the split; equal to
    /// `cutoff` unless set.
    pub end_cutoff: u8,
}

impl QualityOp {
    /// The operation `method` at `cutoff`, splitting at every low-quality base
    /// under `QualityMethod::Runs`.
    pub fn new(method: QualityMethod, cutoff: u8) -> Self {
        QualityOp {
            method,
            cutoff,
            min_low_quality_run: 1,
            end_cutoff: cutoff,
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
    ///
    /// Every piece of `Segments` and `Runs` starts and ends at a base at or
    /// above `cutoff`, so the end trim changes a piece only when `end_cutoff`
    /// is higher; a piece with no base at or above it is dropped.
    fn apply(&self, phred: &[u8]) -> Vec<(usize, usize)> {
        let pieces = match self.method {
            QualityMethod::Ends => return trim_by_quality(phred, self.cutoff),
            QualityMethod::Best => return best_segment(phred, self.cutoff),
            QualityMethod::Segments => maximal_segments(phred, self.cutoff),
            QualityMethod::Runs => split_low_quality(phred, self.cutoff, self.min_low_quality_run),
        };
        if self.end_cutoff <= self.cutoff {
            return pieces;
        }
        pieces
            .into_iter()
            .flat_map(|(s, e)| {
                trim_by_quality(&phred[s..e], self.end_cutoff)
                    .into_iter()
                    .map(move |(ts, te)| (s + ts, s + te))
            })
            .collect()
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

/// One resulting piece of a read after adapter, barcode, crop and quality
/// processing: a kept `[start, end)` span in original read coordinates, with
/// the primer loci of the adapter segment it came from, if the split sheet
/// located any.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Piece {
    /// Piece start.
    pub start: usize,
    /// Piece end, exclusive.
    pub end: usize,
    /// The source adapter segment's 5' primer locus, if one was located.
    pub five: Option<crate::adapter::Locus>,
    /// The source adapter segment's 3' primer locus, if one was located.
    pub three: Option<crate::adapter::Locus>,
}

/// Shared core of `apply` and `apply_pieces`: searches the original read for
/// adapters, widens each segment under `retain`, intersects it with
/// `barcode`, crops the intersection once, applies the quality operation,
/// and emits every kept `[start, end)` interval to `emit` in
/// original-coordinate order, with the source segment's primer loci.
/// Quality splitting can emit multiple intervals per adapter segment, each
/// carrying that segment's loci unchanged.
///
/// `barcode` is the retained interval resolved from the original record's
/// verified `bi` spans. `None` retains the whole read. An unmatched read
/// enters barcode restriction and cropping as one full-length interval with
/// no loci.
///
/// With `retain`, each adapter segment is widened before barcode
/// restriction and cropping: its start moves back to its five locus's start
/// when one is located, and its end moves out to its three locus's end when
/// one is located. A read with no adapters is one segment with no loci, so
/// `retain` changes nothing for it.
fn apply_with(
    seq: &[u8],
    phred: &[u8],
    plan: &TrimPlan,
    adapters: Option<&crate::adapter::AdapterConfig>,
    barcode: Option<(usize, usize)>,
    retain: bool,
    mut emit: impl FnMut(usize, usize, Option<crate::adapter::Locus>, Option<crate::adapter::Locus>),
) {
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
        return;
    }

    let mut process_segment = |s: usize,
                               e: usize,
                               five: Option<crate::adapter::Locus>,
                               three: Option<crate::adapter::Locus>| {
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
        match &plan.quality {
            None => emit(s, e, five, three),
            Some(op) => {
                for (is, ie) in op.apply(wp) {
                    emit(is + s, ie + s, five, three);
                }
            },
        }
    };

    match adapters {
        None => process_segment(0, seq_len, None, None),
        Some(cfg) => {
            for seg in crate::adapter::adapter_segments_annotated(seq, cfg) {
                let (s, e) = if retain {
                    (
                        seg.five.map_or(seg.start, |l| l.start),
                        seg.three.map_or(seg.end, |l| l.end),
                    )
                } else {
                    (seg.start, seg.end)
                };
                process_segment(s, e, seg.five, seg.three);
            }
        },
    }
}

/// `apply_with`, collecting each emitted interval into a `Piece` that
/// carries its source segment's primer loci.
pub fn apply_pieces(
    seq: &[u8],
    phred: &[u8],
    plan: &TrimPlan,
    adapters: Option<&crate::adapter::AdapterConfig>,
    barcode: Option<(usize, usize)>,
    retain: bool,
) -> Vec<Piece> {
    let mut out = Vec::new();
    apply_with(
        seq,
        phred,
        plan,
        adapters,
        barcode,
        retain,
        |start, end, five, three| {
            out.push(Piece {
                start,
                end,
                five,
                three,
            });
        },
    );
    out
}

/// `apply_with` with `retain` false, collecting each emitted interval's
/// `[start, end)` span directly, without building a `Piece` for it.
pub fn apply(
    seq: &[u8],
    phred: &[u8],
    plan: &TrimPlan,
    adapters: Option<&crate::adapter::AdapterConfig>,
    barcode: Option<(usize, usize)>,
) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    apply_with(
        seq,
        phred,
        plan,
        adapters,
        barcode,
        false,
        |s, e, _five, _three| {
            out.push((s, e));
        },
    );
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
            split_of: Vec::new(),
            split_opens: Vec::new(),
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

    /// A read the method leaves whole has its ends trimmed at the stricter end
    /// cutoff under both splitting methods.
    #[test]
    fn end_cutoff_trims_the_ends_of_an_unsplit_read() {
        let phred = [vec![12u8; 5], vec![30; 100], vec![12; 5]].concat();
        let seq = vec![b'A'; phred.len()];
        for method in [QualityMethod::Segments, QualityMethod::Runs] {
            let op = QualityOp::new(method, 10);
            let plan = |quality| TrimPlan {
                quality: Some(quality),
                ..TrimPlan::default()
            };
            assert_eq!(
                apply(&seq, &phred, &plan(op.clone()), None, None),
                vec![(0, 110)],
                "{method:?}"
            );
            let strict = QualityOp {
                end_cutoff: 20,
                ..op
            };
            assert_eq!(
                apply(&seq, &phred, &plan(strict), None, None),
                vec![(5, 105)],
                "{method:?}"
            );
        }
    }

    /// Every piece of a split read has both of its ends trimmed at the end
    /// cutoff, including the ends facing the low-quality region, and a piece
    /// left with no base at or above the end cutoff is dropped.
    #[test]
    fn end_cutoff_trims_every_piece_of_a_split_read() {
        let phred = [
            vec![12u8; 5],
            vec![30; 80],
            vec![12; 3],
            vec![2; 30],
            vec![14; 4],
            vec![30; 80],
            vec![12; 5],
            vec![2; 30],
            vec![15; 60],
        ]
        .concat();
        let seq = vec![b'A'; phred.len()];
        let runs = TrimPlan {
            quality: Some(QualityOp {
                end_cutoff: 20,
                ..QualityOp::runs(10, 10)
            }),
            ..TrimPlan::default()
        };
        assert_eq!(
            apply(&seq, &phred, &runs, None, None),
            vec![(5, 85), (122, 202)]
        );
        let segments = TrimPlan {
            quality: Some(QualityOp {
                end_cutoff: 20,
                ..QualityOp::new(QualityMethod::Segments, 10)
            }),
            ..TrimPlan::default()
        };
        assert_eq!(
            apply(&seq, &phred, &segments, None, None),
            vec![(5, 85), (122, 202)]
        );
    }

    /// An end cutoff at or below the method cutoff changes nothing: every
    /// piece already starts and ends at a base at or above the method cutoff.
    #[test]
    fn end_cutoff_at_or_below_the_cutoff_keeps_the_method_output() {
        let mut state = 5u64;
        for _ in 0..500 {
            let len = 1 + (lcg(&mut state) % 400) as usize;
            let phred: Vec<u8> = (0..len).map(|_| (lcg(&mut state) % 40) as u8).collect();
            let seq = vec![b'A'; len];
            for op in [
                QualityOp::runs(12, 3),
                QualityOp::new(QualityMethod::Segments, 12),
            ] {
                let plan = |quality| TrimPlan {
                    quality: Some(quality),
                    ..TrimPlan::default()
                };
                let want = apply(&seq, &phred, &plan(op.clone()), None, None);
                for end_cutoff in [0, 7, 12] {
                    let lower = QualityOp {
                        end_cutoff,
                        ..op.clone()
                    };
                    assert_eq!(apply(&seq, &phred, &plan(lower), None, None), want);
                }
            }
        }
    }

    /// A deterministic linear congruential generator for the randomized tests.
    fn lcg(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        *state >> 33
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
            split_of: Vec::new(),
            split_opens: Vec::new(),
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

#[cfg(test)]
mod piece_tests {
    use super::*;
    use crate::adapter::{AdapterConfig, reverse_complement};
    use crate::split::Primer;

    /// Generates deterministic SplitMix64 bases.
    fn splitmix_dna(seed: u64, len: usize) -> Vec<u8> {
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(seed);
        (0..len)
            .map(|_| {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                b"ACGT"[((z >> 62) & 0b11) as usize]
            })
            .collect()
    }

    /// Builds `left`, then `mid`, then `right`.
    fn joined(left: &[u8], mid: &[u8], right: &[u8]) -> Vec<u8> {
        let mut out = left.to_vec();
        out.extend_from_slice(mid);
        out.extend_from_slice(right);
        out
    }

    /// A one-target split configuration naming `fwd` and `rev` as the
    /// target's forward and reverse primer, with no other adapters.
    fn split_config(fwd: &[u8], rev: &[u8]) -> AdapterConfig {
        let primers = vec![
            Primer {
                name: "fA".into(),
                seq: fwd.to_vec(),
            },
            Primer {
                name: "rA".into(),
                seq: rev.to_vec(),
            },
        ];
        let mut cfg = AdapterConfig {
            adapters: Vec::new(),
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
            split_of: Vec::new(),
            split_opens: Vec::new(),
        };
        cfg.attach_split(&primers);
        cfg
    }

    const FA: &[u8] = b"ACGGTTCAGCATTGACCGTA";
    const RA: &[u8] = b"TTGCACGGTAACCTGATCGA";

    /// `apply` and `apply_pieces(.., false)` agree on the spans they
    /// produce over planted reads whose quality splits a segment into
    /// several pieces, both with a split adapter configuration (the
    /// `Some(cfg)` branch of `apply_with`) and without one (the `None`
    /// branch).
    #[test]
    fn apply_matches_apply_pieces_spans() {
        let cfg = split_config(FA, RA);
        let mut insert = splitmix_dna(701, 300);
        insert[50..60].fill(b'A');
        insert[150..165].fill(b'A');
        let with_primers = joined(FA, &insert, &reverse_complement(RA));
        let mut phred_with_primers = vec![40u8; with_primers.len()];
        phred_with_primers[FA.len() + 50..FA.len() + 60].fill(2);
        phred_with_primers[FA.len() + 150..FA.len() + 165].fill(2);

        let without_primers = splitmix_dna(702, 300);
        let mut phred_without_primers = vec![40u8; without_primers.len()];
        phred_without_primers[50..60].fill(2);
        phred_without_primers[150..165].fill(2);

        let plan = TrimPlan {
            head: 2,
            tail: 3,
            quality: Some(QualityOp::runs(10, 5)),
        };

        let cases: [(&[u8], &[u8], Option<&crate::adapter::AdapterConfig>); 2] = [
            (&with_primers, &phred_with_primers, Some(&cfg)),
            (&without_primers, &phred_without_primers, None),
        ];
        for (seq, phred, adapters) in cases {
            let spans = apply(seq, phred, &plan, adapters, None);
            let piece_spans: Vec<(usize, usize)> =
                apply_pieces(seq, phred, &plan, adapters, None, false)
                    .into_iter()
                    .map(|p| (p.start, p.end))
                    .collect();
            assert!(
                spans.len() >= 3,
                "expected the two low-quality runs to split into at least 3 pieces: {spans:?}"
            );
            assert_eq!(spans, piece_spans);
        }
    }

    /// `fA` + insert + `rc(rA)`: trimming gives the span between the two
    /// primers; retaining widens it back to the whole read, over both
    /// located loci.
    #[test]
    fn retain_widens_over_primers() {
        let cfg = split_config(FA, RA);
        let insert = splitmix_dna(701, 400);
        let read = joined(FA, &insert, &reverse_complement(RA));
        let n = read.len();
        let phred = vec![40u8; n];
        let plan = TrimPlan::default();

        let trimmed = apply_pieces(&read, &phred, &plan, Some(&cfg), None, false);
        assert_eq!(trimmed.len(), 1, "{trimmed:?}");
        assert_eq!((trimmed[0].start, trimmed[0].end), (FA.len(), n - RA.len()));

        let retained = apply_pieces(&read, &phred, &plan, Some(&cfg), None, true);
        assert_eq!(retained.len(), 1, "{retained:?}");
        assert_eq!((retained[0].start, retained[0].end), (0, n));
        assert_eq!(retained[0].five, trimmed[0].five);
        assert_eq!(retained[0].three, trimmed[0].three);
    }

    /// Every piece a quality split produces from one adapter segment carries
    /// that segment's primer loci.
    #[test]
    fn quality_split_pieces_inherit_segment_loci() {
        let cfg = split_config(FA, RA);
        let mut insert = splitmix_dna(702, 200);
        insert[80..90].fill(b'A');
        let read = joined(FA, &insert, &reverse_complement(RA));
        let mut phred = vec![40u8; read.len()];
        phred[FA.len() + 80..FA.len() + 90].fill(2);
        let plan = TrimPlan {
            head: 0,
            tail: 0,
            quality: Some(QualityOp::runs(10, 5)),
        };

        let pieces = apply_pieces(&read, &phred, &plan, Some(&cfg), None, false);
        assert_eq!(pieces.len(), 2, "{pieces:?}");
        for piece in &pieces {
            assert_eq!(
                piece.five,
                Some(crate::adapter::Locus {
                    start: 0,
                    end: FA.len(),
                    outer_open: false,
                    boundary: false,
                    site: Some(crate::adapter::PrimerSite {
                        entry: 0,
                        rc: false,
                        cost: 0,
                    }),
                })
            );
            assert_eq!(
                piece.three,
                Some(crate::adapter::Locus {
                    start: read.len() - RA.len(),
                    end: read.len(),
                    outer_open: false,
                    boundary: false,
                    site: Some(crate::adapter::PrimerSite {
                        entry: 1,
                        rc: true,
                        cost: 0,
                    }),
                })
            );
        }
    }
}
