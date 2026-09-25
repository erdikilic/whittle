//! De novo adapter inference from recurrent read-end sequences.
//!
//! Exact k-mer graphs establish consensus boundaries. Batched approximate
//! matching validates support and aligns primer extensions to conserved insert
//! starts. Catalog sequences annotate discoveries without choosing their bases.

use crate::adapter::search::{AmbiguousSearcher, hits, is_plain_acgt, new_ambiguous_searcher};
use crate::adapter::{Adapter, AdapterConfig, MIN_PATTERN_LEN, Role, edit_budget};

mod assemble;
mod boundary;
mod consensus;
mod kmer;
mod layers;
mod support;
use assemble::*;
use boundary::*;
use consensus::*;
use kmer::*;
use layers::*;
use support::*;

/// The physical read end from which a consensus was assembled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum End {
    /// Discovered in the 5' windows.
    Five,
    /// Discovered in the 3' windows.
    Three,
}

/// k-mer length used for end-window counting and assembly graph nodes.
const KMER_K: usize = 16;

/// Number of top exact k-mers retained per end for graph assembly.
const TOP_KMERS: usize = 500;

/// Length of the 5'/3' end window scanned per read for adapter discovery.
const WINDOW_LEN: usize = 100;

/// Minimum presence-fraction support required to keep a discovered adapter.
/// Support is the fraction of sampled end windows containing the consensus
/// within its length-scaled edit budget. Absolute support also limits sparse
/// discoveries in small samples.
const KEEP_SUPPORT: f64 = 0.01;

/// Minimum independently supporting read windows for a retained consensus.
const MIN_SUPPORT_WINDOWS: usize = 20;

/// Maximum windows used for alignment and support validation.
const RECOUNT_WINDOWS: usize = 4000;

/// Minimum k-mer support relative to the path peak at an assembly boundary.
const BOUNDARY_SUPPORT: f64 = 0.45;

/// Exact prefix length used to locate recurrent unprimed read starts.
const START_K: usize = 11;

/// Max total emitted length of a single `bounded_heaviest_path` consensus,
/// used by `peel_paths` so no single peel can run away in length.
const LMAX: usize = 100;

/// Maximum distinct insert-boundary anchors evaluated for one consensus.
/// Anchors are ranked by unprimed read-start support.
const MAX_BOUNDARY_ANCHORS: usize = 8;

/// Minimum fraction of a layer's windows that one member of a variable
/// layer must hold, below which a barcode panel would exceed 1000 members.
const VARIABLE_MEMBER_SUPPORT: f64 = 0.001;

/// Max number of adapters `peel_paths` will extract from one end's k-mer graph.
/// A barcode layer holds one family per barcode.
const MAX_ADAPTERS_PER_END: usize = 128;

/// Maximum discovery layers at one read end.
const MAX_LAYERS: usize = 8;

/// Bases from the physical read end searched for a candidate and its mirror.
const MIRROR_WINDOW: usize = 300;

/// Window tested for a mirror at the opposite end.
const MIRROR_SEGMENT: usize = KMER_K;

/// Depth difference between an inner segment and its mirror accepted as
/// symmetric, in bases.
const SYMMETRY_TOLERANCE: usize = 35;

/// Fraction of the terminal error rate allowed when locating a mirror, so
/// that a mirror does not extend past the sequence present at the read end.
const MIRROR_ERROR_RATE: f64 = 0.5;

/// Bases by which the tested window slides per mirror search.
const MIRROR_STEP: usize = 4;

/// Divisor of a candidate's support giving the fewest mirror occurrences
/// that count. A mirror occurs only in reads of the other orientation and is
/// searched at the stricter mirror error rate.
const MIRROR_SUPPORT_DIVISOR: usize = 8;

/// Percentage of a candidate's supporting windows in which it must start
/// within `ANCHOR_SLACK` of the boundary. A majority suffices once earlier
/// layers explain only part of the reads.
const ANCHORED_PERCENT: usize = 60;

/// Bases from the current boundary within which a hit counts as anchored at
/// it and advances it. The slack covers eroded or variable-length remnants
/// of the preceding layer.
const ANCHOR_SLACK: usize = 50;

/// Largest median distance from the physical read end at which an outermost
/// discovered sequence is a sequencing adapter, which splits reads at
/// interior hits; deeper sequences are end-only.
const ADAPTER_FLUSH: usize = 20;

/// Shortest end shared by candidates of one layer that is treated as a
/// neighbouring layer rather than part of the candidate.
const SHARED_END_MIN: usize = 12;

/// Percentage of the shorter of two candidates that their longest common
/// substring must cover for them to count as one family.
const FAMILY_OVERLAP_PERCENT: usize = 60;

/// Percentage of a candidate's supporting windows that a stronger candidate
/// of the same layer must share for the candidate to count as a sequencing
/// variant of it rather than a distinct family.
const VARIANT_OVERLAP_PERCENT: usize = 80;

/// Minimum percent identity for a catalog entry to be reported as the match
/// of an inferred adapter. A 16 to 32 bp anchor searched against every catalog
/// entry on both strands names something spurious well above the 60 percent
/// that its trimming budget alone would allow.
const NAME_IDENTITY_MIN: f32 = 85.0;

/// One discovered adapter with inference metadata. The bare `Adapter` (without
/// `support` and `name_hits`) is extracted only when building the trim config.
#[derive(Debug, Clone)]
pub struct InferredAdapter {
    /// Sequence used for trimming (or printed as the recommendation), named
    /// `inferred_N` by presentation order.
    pub adapter: Adapter,
    /// Complete sequence retained after boundary validation.
    pub assembled_seq: Vec<u8>,
    /// Fraction of sampled end windows containing the consensus within its
    /// edit budget.
    pub support: f64,
    /// Catalog entries within `NAME_IDENTITY_MIN` of the consensus, best
    /// first, as `(name, percent identity)`. An annotation, not the name.
    pub name_hits: Vec<(String, f32)>,
    /// Discovery layer, counted from the read end after any known sequences.
    pub layer: usize,
}

impl InferredAdapter {
    /// Returns the number of assembled bases excluded from the trimming sequence.
    pub fn uncertain_bases(&self) -> usize {
        self.assembled_seq
            .len()
            .saturating_sub(self.adapter.seq.len())
    }
}

/// A merged layer family ready for variant suppression: sequence, support,
/// supporting windows, their count, whether it carries an insert boundary,
/// its end, and its path weight.
type Ranked = (Vec<u8>, f64, Vec<bool>, usize, bool, End, u64);

/// Discovers supported technical sequences in layers from each read end. Known
/// sequences in `base` explain the outermost layers first; each accepted layer
/// moves the boundary inward and the next layer is assembled from the
/// unexplained sequence. A barcode construct is not a known sequence: its `N`
/// block would explain any read end. Equivalent assemblies share one trimming
/// pattern. Catalog and supplied FASTA entries provide names only after the
/// inferred boundaries are fixed.
pub fn discover(sample: &[&[u8]], base: &AdapterConfig) -> Vec<InferredAdapter> {
    let mut bounds = Boundaries::new(sample);
    let known: Vec<Vec<u8>> = base
        .adapters
        .iter()
        .filter(|a| !crate::adapter::is_construct(a))
        .map(|a| a.seq.to_ascii_uppercase())
        .collect();
    if !known.is_empty() {
        let mut active = vec![true; sample.len()];
        for _ in 0..MAX_LAYERS {
            active = advance_boundaries(
                sample,
                &mut bounds,
                &known,
                base.error_rate,
                &active,
                true,
                true,
            );
            if active.is_empty() {
                break;
            }
        }
    }
    let physical = Boundaries::new(sample);
    let (five_p, three_p) = layer_windows(sample, &physical, MIRROR_WINDOW);
    let five_phys: Vec<&[u8]> = stride_sample(
        &five_p.iter().map(|(_, w)| *w).collect::<Vec<_>>(),
        RECOUNT_WINDOWS,
    );
    let three_phys: Vec<&[u8]> = stride_sample(
        &three_p.iter().map(|(_, w)| *w).collect::<Vec<_>>(),
        RECOUNT_WINDOWS,
    );

    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let mut strand_searcher = crate::adapter::search::new_searcher_fwd();
    let mut distinct: Vec<(Vec<u8>, f64, usize, bool, bool)> = Vec::new();
    // An accepted primer marks the insert boundary at its end; no layer lies
    // beyond it.
    let mut open = [true, true];
    // Opening k-mers of the strand-specific layers behind the divisions of
    // earlier layers.
    let mut strand_words: Vec<Vec<u8>> = Vec::new();
    for layer in 0..MAX_LAYERS {
        let (mut five_w, mut three_w) = layer_windows(sample, &bounds, 2 * WINDOW_LEN);
        if !open[0] {
            five_w.clear();
        }
        if !open[1] {
            three_w.clear();
        }
        let five_texts: Vec<&[u8]> = five_w.iter().map(|(_, w)| *w).collect();
        let three_texts: Vec<&[u8]> = three_w.iter().map(|(_, w)| *w).collect();
        // Ranking statistics use windows distributed across the sample.
        let five_sample = stride_sample(&five_texts, RECOUNT_WINDOWS);
        let three_sample = stride_sample(&three_texts, RECOUNT_WINDOWS);
        let mut five = assemble(&five_texts, &three_phys, base, End::Five, layer == 0);
        let mut three = assemble(&three_texts, &five_phys, base, End::Three, layer == 0);
        strip_shared_ends(&mut five, base.error_rate);
        strip_shared_ends(&mut three, base.error_rate);
        let (five, five_variable) = with_variable_layer(
            &mut searcher,
            five,
            &five_texts,
            &five_sample,
            End::Five,
            base.error_rate,
        );
        let (three, three_variable) = with_variable_layer(
            &mut searcher,
            three,
            &three_texts,
            &three_sample,
            End::Three,
            base.error_rate,
        );
        let variable: Vec<Vec<u8>> = five_variable.into_iter().chain(three_variable).collect();
        let background = background_windows(sample, &bounds, WINDOW_LEN);
        let background = stride_sample(&background, RECOUNT_WINDOWS);

        // Insert stretches identified by their mirror also reject the graph
        // fragments that reconstruct part of them.
        let mut inserts: Vec<Vec<u8>> = Vec::new();
        let mut next_words: Vec<Vec<u8>> = Vec::new();
        let cut: Vec<(Vec<u8>, f64, bool, u64, End)> = five
            .into_iter()
            .map(|c| (c, End::Five))
            .chain(three.into_iter().map(|c| (c, End::Three)))
            .filter(|((seq, support, _, _, _, _), _)| {
                seq.len() >= MIN_PATTERN_LEN && (*support >= KEEP_SUPPORT || variable.contains(seq))
            })
            .filter_map(
                |((seq, support, boundary, weight, unbounded, divided), end)| {
                    // A division into strand-specific layers marks the
                    // candidate and the layers it divides into technical;
                    // their mirrors lie at a depth that depends on how much
                    // outer adapter each end keeps.
                    let strand = strand_words.iter().any(|word| {
                        !hits(
                            &mut strand_searcher,
                            word,
                            &seq,
                            edit_budget(base.error_rate, word.len()),
                        )
                        .is_empty()
                    });
                    if !divided.is_empty() || strand {
                        next_words.extend(divided);
                        return Some((seq, support, boundary, weight, end));
                    }
                    let (own, opposite) = match end {
                        End::Five => (&five_phys, &three_phys),
                        End::Three => (&three_phys, &five_phys),
                    };
                    // A mirror at the opposite end is the only insert evidence
                    // for a candidate the assembly window could not bound.
                    let (cut, insert) = symmetry_cut(&seq, end, own, opposite, base.error_rate);
                    if insert.len() >= MIN_PATTERN_LEN {
                        inserts.push(crate::adapter::reverse_complement(&insert));
                        inserts.push(insert);
                    }
                    let cut = cut?;
                    if unbounded && cut.len() == seq.len() {
                        return None;
                    }
                    Some((cut, support, boundary, weight, end))
                },
            )
            .collect();
        strand_words.extend(next_words);
        let mut candidates: Vec<(Vec<u8>, f64, u32, bool, u64, End)> = cut
            .into_iter()
            .filter(|(seq, _, _, _, _)| {
                !inserts
                    .iter()
                    .any(|insert| same_family(seq, insert, base.error_rate))
            })
            .filter_map(|(seq, support, boundary, weight, end)| {
                let count = windows_containing(
                    &mut searcher,
                    &seq,
                    &background,
                    edit_budget(base.error_rate, seq.len()),
                );
                if !background.is_empty() && count as f64 * 4.0 >= support * background.len() as f64
                {
                    return None;
                }
                let exact = windows_containing(&mut searcher, &seq, &five_sample, 0)
                    + windows_containing(&mut searcher, &seq, &three_sample, 0);
                tracing::debug!(sequence = %String::from_utf8_lossy(&seq), exact, "Layer candidate");
                Some((seq, support, exact, boundary, weight, end))
            })
            .collect();
        // Independently supported insert boundaries take precedence over graph
        // fragments that may include conserved insert sequence. Exact support
        // then distinguishes reconstructions from nearby sequencing-error variants.
        candidates.sort_by(|a, b| {
            b.3.cmp(&a.3)
                .then_with(|| {
                    if a.3 && b.3 {
                        b.2.cmp(&a.2)
                    } else {
                        std::cmp::Ordering::Equal
                    }
                })
                .then(b.4.cmp(&a.4))
                .then(b.2.cmp(&a.2))
                .then(b.1.total_cmp(&a.1))
                .then(b.0.len().cmp(&a.0.len()))
                .then(a.0.cmp(&b.0))
        });
        // Contained reconstructions and fragments merge into the heaviest
        // reconstruction of their family. A sequence supported `RUN_JUMP`
        // times better than the candidate containing it is the shared layer
        // of that candidate, not its fragment.
        let mut merged: Vec<(Vec<u8>, f64, u64, bool, End)> = Vec::new();
        for (seq, support, _, boundary, weight, end) in candidates {
            let matched = weight;
            // A family found again in a later layer, in reads its first
            // occurrence did not match, is not a new layer; a shorter
            // sequence contained in an earlier candidate is a layer that the
            // candidate fused with its neighbour. The heaviest
            // reconstruction represents a family with its own support.
            if distinct.iter().any(|(other, _, _, _, _)| {
                other.len() <= seq.len() + edit_budget(base.error_rate, seq.len())
                    && same_adapter(&seq, other, base.error_rate)
            }) {
                continue;
            } else if let Some((other, previous, best, _, _)) =
                merged.iter_mut().find(|(other, previous, _, _, _)| {
                    support < *previous * f64::from(RUN_JUMP)
                        && *previous < support * f64::from(RUN_JUMP)
                        && fragment_of(&seq, other, base.error_rate)
                })
            {
                if matched > *best {
                    *other = seq;
                    *best = matched;
                    *previous = support;
                }
            } else {
                merged.push((seq, support, matched, boundary, end));
            }
        }
        // Sequencing variants of a family occupy the same reads as the
        // family; distinct families occupy different reads.
        let layer_texts: Vec<&[u8]> = five_sample.iter().chain(&three_sample).copied().collect();
        let mut ranked: Vec<Ranked> = merged
            .into_iter()
            .map(|(seq, support, matched, boundary, end)| {
                let k = edit_budget(base.error_rate, seq.len());
                let mut windows = windows_covered(&mut searcher, &seq, &layer_texts, k);
                let mirror = crate::adapter::reverse_complement(&seq);
                for (window, hit) in
                    windows
                        .iter_mut()
                        .zip(windows_covered(&mut searcher, &mirror, &layer_texts, k))
                {
                    *window |= hit;
                }
                let own = windows.iter().filter(|&&w| w).count();
                (seq, support, windows, own, boundary, end, matched)
            })
            .collect();
        // Complete reconstructions rank above fragments that run into the
        // insert, which match more windows over fewer exact bases.
        ranked.sort_by(|a, b| {
            b.6.cmp(&a.6)
                .then(b.3.cmp(&a.3))
                .then(b.0.len().cmp(&a.0.len()))
                .then(a.0.cmp(&b.0))
        });
        let mut accepted: Vec<Vec<u8>> = Vec::new();
        let mut accepted_windows: Vec<Vec<bool>> = Vec::new();
        let mut accepted_ends: Vec<End> = Vec::new();
        for (seq, support, windows, own, boundary, end, matched) in ranked {
            // Members of a variable layer share their reads with the
            // constant layer behind them by construction.
            let member = variable
                .iter()
                .any(|v| same_adapter(&seq, v, base.error_rate));
            let shared = accepted_windows
                .iter()
                .map(|other| {
                    windows
                        .iter()
                        .zip(other)
                        .filter(|(a, b)| **a && **b)
                        .count()
                })
                .max()
                .unwrap_or(0);
            let variant = !member && shared * 100 >= own * VARIANT_OVERLAP_PERCENT;
            tracing::debug!(sequence = %String::from_utf8_lossy(&seq), matched, own, shared, variant, "Ranked candidate");
            if variant {
                continue;
            }
            let flush = layer == 0
                && known.is_empty()
                && [(End::Five, &five_phys), (End::Three, &three_phys)]
                    .iter()
                    .any(|(end, phys)| {
                        outer_depth(&mut searcher, &seq, phys, *end, base.error_rate)
                            <= ADAPTER_FLUSH
                    });
            if boundary {
                open[usize::from(end == End::Three)] = false;
            }
            accepted_ends.push(end);
            accepted.push(seq.clone());
            if !member {
                accepted_windows.push(windows);
            }
            distinct.push((seq, support, layer, flush, member));
        }
        tracing::debug!(
            layer = layer + 1,
            accepted = accepted.len(),
            "Discovery layer"
        );
        // An end without a layer here has none deeper either.
        for (index, end) in [End::Five, End::Three].into_iter().enumerate() {
            if !accepted_ends.contains(&end) {
                open[index] = false;
            }
        }
        if accepted.is_empty()
            || advance_boundaries(
                sample,
                &mut bounds,
                &accepted,
                base.error_rate,
                &vec![true; sample.len()],
                false,
                layer == 0,
            )
            .is_empty()
        {
            break;
        }
    }

    let refs = crate::adapter::preset::preset(crate::adapter::preset::Kit::ALL);
    let name_refs: Vec<Adapter> = refs
        .into_iter()
        .chain(base.adapters.iter().cloned())
        .collect();
    distinct.sort_by(|a, b| a.2.cmp(&b.2).then(b.1.total_cmp(&a.1)).then(a.0.cmp(&b.0)));
    distinct
        .into_iter()
        .enumerate()
        .map(|(i, (seq, support, layer, flush, member))| {
            let name_hits = name_against(&seq, &name_refs, base.error_rate);
            // A flush layer is an adapter unless it is a marker-gene primer,
            // whose site also lies inside genomic reads.
            let role = if flush && !crate::adapter::matches_marker_primer(&seq, base.error_rate) {
                Role::Adapter
            } else if member {
                Role::Barcode
            } else {
                Role::Primer
            };
            InferredAdapter {
                adapter: Adapter {
                    name: format!("inferred_{}", i + 1),
                    seq: crate::adapter::with_marker_codes(&seq, base.error_rate),
                    role,
                },
                assembled_seq: seq,
                support,
                name_hits,
                layer,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests;
