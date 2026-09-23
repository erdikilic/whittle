//! Assembly of the candidates of one read end and one layer: graph paths,
//! boundary validation and support measurement.

use super::*;

/// An assembled end candidate: sequence, support, whether an independent
/// insert boundary was found, the summed original k-mer support of the
/// retained span, and whether the assembly window could not bound the
/// insert-facing side.
pub(super) type Candidate = (Vec<u8>, f64, bool, u64, bool);

/// Assembles one end's candidates and validates support and insert boundaries.
/// `contrast` enables primer reconstruction from recurrent unprimed window
/// starts, which describe an insert boundary only in the first discovery
/// layer.
pub(super) fn assemble(
    windows: &[&[u8]],
    base: &AdapterConfig,
    end: End,
    contrast: bool,
) -> Vec<Candidate> {
    if windows.len() < 3 {
        return Vec::new();
    }
    // K-mer encoding and approximate matching operate on uppercase DNA.
    // Inference owns normalized copies and does not modify pipeline records.
    let upper: Vec<Vec<u8>> = windows.iter().map(|w| w.to_ascii_uppercase()).collect();
    let windows: Vec<&[u8]> = upper.iter().map(Vec::as_slice).collect();
    let windows = windows.as_slice();

    let assembly_windows: Vec<&[u8]> = windows
        .iter()
        .map(|w| match end {
            End::Five => &w[..WINDOW_LEN.min(w.len())],
            End::Three => &w[w.len().saturating_sub(WINDOW_LEN)..],
        })
        .collect();
    let exact = top_kmers(&assembly_windows, KMER_K, TOP_KMERS);
    if exact
        .first()
        .is_none_or(|&(_, count)| count < MIN_SUPPORT_WINDOWS as u32)
    {
        return Vec::new();
    }
    // Validation uses windows distributed across the complete sample.
    let recount = stride_sample(windows, RECOUNT_WINDOWS);
    let n_recount = recount.len();
    let composition = base_composition(&recount);
    let following = FollowingBases::new(&recount, end);

    let original_weights: std::collections::HashMap<u64, u32> = exact.iter().copied().collect();
    let weighted = exact;
    let mut out = Vec::new();
    let mut insert_starts = Vec::new();
    let mut primer_edges = vec![None; recount.len()];
    let mut primer_searcher = crate::adapter::search::new_searcher_fwd();
    let mut known_primers = std::collections::HashSet::new();
    // Contrast reads the validation windows outward from the end boundary.
    let reversed: Vec<Vec<u8>> = if contrast && end == End::Three {
        recount
            .iter()
            .map(|w| w.iter().rev().copied().collect())
            .collect()
    } else {
        Vec::new()
    };
    let oriented: Vec<&[u8]> = if reversed.is_empty() {
        recount.clone()
    } else {
        reversed.iter().map(Vec::as_slice).collect()
    };
    for (cons, _) in peel_paths(weighted, KMER_K, end) {
        tracing::debug!(sequence = %String::from_utf8_lossy(&cons), "Assembled end candidate");
        // Original weights show the support of every k-mer of the path,
        // including k-mers earlier peels removed from the graph.
        let weights: Vec<u32> = cons
            .windows(KMER_K)
            .map(|w| {
                encode_kmer(w)
                    .and_then(|code| original_weights.get(&code).copied())
                    .unwrap_or(0)
            })
            .collect();
        let (lo, hi, rises) = supported_span(&weights, end);
        let span = &weights[lo..=hi - KMER_K];
        let peak = span.iter().copied().max().unwrap_or(0);
        // The path weight measures completeness: a fragment running into
        // the insert or a sequencing variant carries weak k-mers.
        let weight: u64 = span.iter().map(|&w| u64::from(w)).sum();
        if peak < MIN_SUPPORT_WINDOWS as u32 {
            continue;
        }
        let trimmed = cons[lo..hi].to_vec();
        if is_repetitive(&trimmed) {
            continue;
        }
        // A variant or fragment of a stronger validated candidate merges
        // into it and needs no validation of its own.
        if out.iter().any(|(seq, _, _, other, _): &Candidate| {
            *other >= weight && same_adapter(&trimmed, seq, base.error_rate)
        }) {
            continue;
        }
        if !known_primers.is_empty()
            && predominantly_insert(
                &trimmed,
                &recount,
                end,
                edit_budget(base.error_rate, trimmed.len()),
                &primer_edges,
                &mut primer_searcher,
            )
        {
            continue;
        }
        let trimmed = if peak as usize * 4 < windows.len() {
            polish_consensus(&trimmed, &recount)
        } else {
            trimmed
        };
        let trimmed = trim_unconserved_inner_end(&trimmed, &following, end, composition);
        let contrast = contrast
            .then(|| contrast_boundary(&cons, &oriented, end))
            .flatten();
        let has_insert_boundary = contrast.is_some();
        let mut unbounded = false;
        let mut measured = None;
        let sequences = if let Some(boundary) = contrast {
            insert_starts.push(boundary.insert_start);
            boundary.primers
        } else {
            let word = match end {
                End::Five => &cons[hi - KMER_K..hi],
                End::Three => &cons[lo..lo + KMER_K],
            };
            let code = encode_kmer(word).unwrap();
            let boundary_weight = original_weights[&code];
            let mask = (1u64 << (2 * KMER_K)) - 1;
            let (continuation, next) = (0..4)
                .map(|base| {
                    let next = match end {
                        End::Five => ((code << 2) | base) & mask,
                        End::Three => (code >> 2) | (base << (2 * (KMER_K - 1))),
                    };
                    (original_weights.get(&next).copied().unwrap_or(0), next)
                })
                .max()
                .unwrap_or((0, 0));
            let observed_boundary = match end {
                End::Five => hi < cons.len(),
                End::Three => lo > 0,
            };
            // Support drops at an insert boundary. It rises where a layer
            // such as a barcode joins a shared downstream layer: between the
            // support plateaus of the path, or where the downstream k-mer
            // occurs in `RUN_JUMP` times more windows than contain the
            // candidate at all. Assembly windows shorter than the layer
            // stack depress the downstream k-mer count; a twofold excess in
            // those counts is confirmed on the validation windows, where a
            // variant path joining its own family never reaches twofold. A
            // boundary inside one technical sequence changes support little.
            let drop = continuation * 2 <= boundary_weight;
            let (present, anchored) = terminal_support(
                &trimmed,
                &recount,
                end,
                edit_budget(base.error_rate, trimmed.len()),
            );
            measured = Some((present, anchored));
            let rise = rises
                || present.saturating_mul(RUN_JUMP as usize) <= continuation as usize
                || (continuation >= 2 * boundary_weight && {
                    let downstream = windows_containing(
                        &mut primer_searcher,
                        &decode_kmer(next, KMER_K),
                        &recount,
                        edit_budget(base.error_rate, KMER_K),
                    );
                    downstream as usize >= 2 * present
                });
            // The read base after the candidate marks an insert boundary
            // directly when it is not conserved (see `FollowingBases`);
            // erosion at the physical end depresses the boundary k-mer and can
            // hide the drop in k-mer support.
            let follows_insert = following
                .after(&trimmed, trimmed.len(), end)
                .is_some_and(|counts| !conserved(counts, composition));
            let terminated = (drop && supported_termination(word, &recount, end)) || follows_insert;
            let bounded = (observed_boundary || cons.len() < LMAX) && (rise || terminated);
            if !bounded {
                tracing::debug!(sequence = %String::from_utf8_lossy(&trimmed), observed_boundary,
                    boundary_weight, continuation, present, drop, terminated,
                    "Recurrent sequence has no supported insert boundary");
            }
            unbounded = !bounded;
            vec![trimmed]
        };
        for trimmed in sequences {
            if trimmed.len() < MIN_PATTERN_LEN || is_repetitive(&trimmed) {
                continue;
            }
            // Presence counts each supporting window once, including reads
            // whose sequencing errors disrupted individual exact k-mers.
            let k_cons = edit_budget(base.error_rate, trimmed.len());
            let (present, anchored) = measured
                .take()
                .unwrap_or_else(|| terminal_support(&trimmed, &recount, end, k_cons));
            let support = present as f64 / n_recount as f64;
            tracing::debug!(sequence = %String::from_utf8_lossy(&trimmed), present, anchored, peak, has_insert_boundary, "Validated end candidate");
            if present >= MIN_SUPPORT_WINDOWS
                && support >= KEEP_SUPPORT
                && anchored * 100 >= present * ANCHORED_PERCENT
            {
                if has_insert_boundary && known_primers.insert(trimmed.clone()) {
                    for hit in primer_searcher.search_texts(&trimmed, &recount, k_cons) {
                        let edge = match end {
                            End::Five => hit.text_end,
                            End::Three => hit.text_start,
                        };
                        primer_edges[hit.text_idx] = Some((edge, trimmed.len() / 2));
                    }
                }
                out.push((trimmed, support, has_insert_boundary, weight, unbounded));
            }
        }
    }
    out.retain(|(seq, _, boundary, _, _)| {
        *boundary
            || !insert_starts
                .iter()
                .any(|start| seq.windows(start.len()).any(|word| word == start))
    });
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    if primer_edges.iter().any(Option::is_some) {
        out.retain(|(seq, _, boundary, _, _)| {
            if *boundary {
                return true;
            }
            !predominantly_insert(
                seq,
                &recount,
                end,
                edit_budget(base.error_rate, seq.len()),
                &primer_edges,
                &mut searcher,
            )
        });
    }
    out
}
