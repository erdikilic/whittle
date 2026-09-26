//! Assembly of the candidates of one read end and one layer: graph paths,
//! boundary validation and support measurement.

use super::*;

/// An assembled end candidate: sequence, support, whether an independent
/// insert boundary was found, the summed original k-mer support of the
/// retained span, whether the assembly window could not bound the
/// insert-facing side, and the opening k-mers of the strand-specific layers
/// it divides into, empty unless it ends at a division (`division_point`).
pub(super) type Candidate = (Vec<u8>, f64, bool, u64, bool, Vec<Vec<u8>>);

/// Assembles one end's candidates and validates support and insert boundaries.
/// `contrast` enables primer reconstruction from recurrent unprimed window
/// starts, which describe an insert boundary only in the first discovery
/// layer. `opposite` holds physical windows of the other read end, where
/// the continuations of a division into strand-specific layers recur.
pub(super) fn assemble(
    windows: &[&[u8]],
    opposite: &[&[u8]],
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
    let mut deep = None;
    // Whether a candidate at the boundary follows the template, so that the
    // insert begins there for the reads holding it.
    let mut at_template = false;
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
        let (mut lo, mut hi, mut rises) = supported_span(&weights, end);
        // A shared layer ends where its path divides into the layers behind
        // it, a boundary like a rise in support. The division holds when both
        // continuations recur reverse complemented at the other read end, as
        // the primers of the two strands do; variants of a conserved insert
        // do not.
        let deep = deep.get_or_insert_with(|| kmer_counts(&recount));
        let mut divided = Vec::new();
        if let Some((edge, continuations)) = division_point(&cons, deep, recount.len(), lo, hi, end)
            && continuations.iter().all(|word| {
                let mirrored = windows_containing(
                    &mut primer_searcher,
                    &crate::adapter::reverse_complement(word),
                    opposite,
                    edit_budget(MIRROR_ERROR_RATE * base.error_rate, KMER_K),
                ) as usize;
                mirrored >= MIN_SUPPORT_WINDOWS
                    && mirrored * MIRROR_SUPPORT_DIVISOR >= opposite.len()
            })
        {
            tracing::debug!(sequence = %String::from_utf8_lossy(&cons), edge, "Path divides");
            match end {
                End::Five => hi = edge,
                End::Three => lo = edge,
            }
            rises = true;
            divided = continuations.to_vec();
        }
        // A homopolymer run at the inner end, such as the poly(A) tail of a
        // transcript behind a cDNA primer, varies in length between reads
        // and belongs to the insert; the layer ends where the run starts.
        let run = inner_run(&cons[lo..hi], end);
        if run >= PLATEAU && hi - lo >= run + KMER_K.max(MIN_PATTERN_LEN) {
            match end {
                End::Five => hi -= run,
                End::Three => lo += run,
            }
            rises = true;
        }
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
        if out.iter().any(|(seq, _, _, other, _, _): &Candidate| {
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
        let mut trimmed = trim_unconserved_inner_end(&trimmed, &following, end, composition);
        // The insert of an amplicon begins at the inner end of a marker-gene
        // primer, or at the first column whose bases follow the template. A
        // candidate that is insert from its outer bases is dropped; when too
        // few bases lie before it at the boundary to hold a technical layer,
        // the insert begins at the boundary itself. A division into strand
        // primers already accounts for the alternative bases behind a
        // candidate.
        let start = if !divided.is_empty() {
            (None, 0)
        } else if let Some(keep) = marker_primer_end(&trimmed, end, base.error_rate) {
            let start = if keep < trimmed.len() {
                TemplateStart::Inside(keep)
            } else {
                TemplateStart::Behind
            };
            (Some(start), usize::MAX)
        } else {
            template_start(&trimmed, &recount, end, base.error_rate)
        };
        let mut template = false;
        match start {
            (Some(TemplateStart::Before(distance)), depth) => {
                tracing::debug!(sequence = %String::from_utf8_lossy(&trimmed), distance, depth, "Candidate lies in the template");
                at_template |= depth < MIN_PATTERN_LEN + distance;
                continue;
            },
            (Some(TemplateStart::Inside(column)), depth) if column < MIN_PATTERN_LEN => {
                tracing::debug!(sequence = %String::from_utf8_lossy(&trimmed), column, depth, "Candidate follows the template");
                at_template |= depth + column < MIN_PATTERN_LEN;
                continue;
            },
            (Some(TemplateStart::Inside(column)), _) => {
                trimmed = outer_part(&trimmed, column, end);
                template = true;
            },
            (Some(TemplateStart::Behind), _) => template = true,
            (None, _) => {},
        }
        if template {
            tracing::debug!(sequence = %String::from_utf8_lossy(&trimmed), "Template begins behind candidate");
        }
        let contrast = (contrast && !template)
            .then(|| contrast_boundary(&cons, &oriented, end))
            .flatten();
        let has_insert_boundary = contrast.is_some() || template;
        let mut unbounded = false;
        let mut measured = None;
        let sequences = if let Some(boundary) = contrast {
            insert_starts.push(boundary.insert_start);
            boundary.primers
        } else if template {
            vec![trimmed]
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
                out.push((
                    trimmed,
                    support,
                    has_insert_boundary,
                    weight,
                    unbounded,
                    divided.clone(),
                ));
            }
        }
    }
    // Behind a boundary where the insert begins, only candidates that end at
    // an insert boundary themselves are technical.
    out.retain(|(seq, _, boundary, _, _, _)| {
        *boundary
            || (!at_template
                && !insert_starts
                    .iter()
                    .any(|start| seq.windows(start.len()).any(|word| word == start)))
    });
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    if primer_edges.iter().any(Option::is_some) {
        out.retain(|(seq, _, boundary, _, _, _)| {
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
    // A candidate that reads past the inner end of a primer in the reads
    // holding both continues into the insert that begins there; it ends
    // where the primer does. A candidate that ends inside such a primer, as
    // an outer layer assembled together with the primer's first bases does,
    // is continued to the primer's inner end.
    let spans: Vec<Vec<Option<(usize, usize)>>> = out
        .iter()
        .map(|(seq, _, _, _, _, _)| {
            hit_spans(seq, &recount, end, edit_budget(base.error_rate, seq.len()))
        })
        .collect();
    let flagged: Vec<bool> = out
        .iter()
        .map(|(_, _, boundary, _, _, _)| *boundary)
        .collect();
    let mut cut = vec![false; out.len()];
    for (i, (seq, _, boundary, _, unbounded, _)) in out.iter_mut().enumerate() {
        let excess = (0..spans.len())
            .filter(|&j| j != i && flagged[j])
            .map(|j| read_through(&spans[i], &spans[j]))
            .max()
            .unwrap_or(0);
        if excess > 0 {
            tracing::debug!(sequence = %String::from_utf8_lossy(seq), excess, "Candidate reads past a primer");
            *seq = outer_part(seq, seq.len().saturating_sub(excess), end);
            *boundary = true;
            *unbounded = false;
            cut[i] = true;
        }
    }
    let primers: Vec<(usize, Vec<u8>)> = out
        .iter()
        .enumerate()
        .filter(|&(j, _)| flagged[j] && !cut[j])
        .map(|(j, (seq, _, _, _, _, _))| (j, seq.clone()))
        .collect();
    for (i, (seq, _, boundary, _, unbounded, _)) in out.iter_mut().enumerate() {
        if flagged[i] || cut[i] {
            continue;
        }
        if let Some(joined) = primers
            .iter()
            .find_map(|(j, primer)| continued_into(seq, primer, &spans[i], &spans[*j], end))
        {
            tracing::debug!(sequence = %String::from_utf8_lossy(&joined), "Candidate continues into a primer");
            *seq = joined;
            *boundary = true;
            *unbounded = false;
        }
    }
    out.retain(|(seq, _, _, _, _, _)| seq.len() >= MIN_PATTERN_LEN);
    out
}
