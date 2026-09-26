//! Degenerate layers: random tags such as UMIs, whose bases differ between
//! reads while each position keeps a constrained composition. Exact k-mers
//! cannot assemble them, so they are read from the per-position base
//! composition of the reads behind an accepted layer.

use super::*;

/// Positions profiled behind an accepted layer for a degenerate layer.
const DEGENERATE_SPAN: usize = 3 * KMER_K;

/// Factor below its composition share at which a base is absent from a
/// position of a degenerate layer. A random tag excludes a base by design; a
/// genomic position holds every base near its share, whatever the genome's
/// composition.
const ABSENT_FACTOR: f64 = 4.0;

/// The code of one profiled position: a single conserved base, the set of
/// bases present at a position that excludes at least one base, or `None`
/// when every base is present.
pub(super) fn position_code(counts: [usize; 4], composition: [f64; 4]) -> Option<u8> {
    let total = counts.iter().sum::<usize>() as f64;
    if let Some(base) = (0..4).find(|&i| 2.0 * counts[i] as f64 >= (1.0 + composition[i]) * total) {
        return Some(b"ACGT"[base]);
    }
    let mask = (0..4)
        .filter(|&i| counts[i] as f64 * ABSENT_FACTOR > composition[i] * total)
        .fold(0, |mask, i| mask | (1 << i));
    (mask != 0b1111).then(|| ambiguity_code(mask))
}

/// Returns the degenerate layer behind `pattern` in `windows`, which read
/// inward from a physical read end, with the number of windows holding
/// `pattern`, or `None` without one.
///
/// The layer grows one position at a time. Each round aligns `pattern` and
/// the layer so far to the windows that hold `pattern`, so that indels shift
/// the supporting reads without shifting the layer, and tallies the base
/// that follows each alignment. A position is conserved when one base holds
/// the midpoint between its share of `composition` and one, as in
/// `conserved`, or when at least one base falls `ABSENT_FACTOR` below its
/// share, and it is written with the ambiguity code of the bases present;
/// the conserved positions end at the first position that is not. A random
/// tag closes with the homopolymer run of fixed bases after its last
/// degenerate position, and the layer ends there: a further fixed run, such
/// as the G run a template switch adds, varies in length between reads and
/// is left to the insert, as a poly(A) run is. The layer needs `KMER_K`
/// positions and support at the `KEEP_SUPPORT` and `MIN_SUPPORT_WINDOWS`
/// floors, and its degenerate positions must be random: no combination of
/// their bases recurs in `MIN_SUPPORT_WINDOWS` windows, and at least half of
/// the windows hold distinct combinations. A panel of barcodes, the conserved
/// starts of an amplicon's species or the primers of two strands concentrate
/// in a few combinations instead.
pub(super) fn degenerate_layer(
    searcher: &mut AmbiguousSearcher,
    pattern: &[u8],
    windows: &[&[u8]],
    composition: [f64; 4],
    error_rate: f64,
) -> Option<(Vec<u8>, usize)> {
    let floor = MIN_SUPPORT_WINDOWS.max((windows.len() as f64 * KEEP_SUPPORT).ceil() as usize);
    let mut best: Vec<Option<(i32, usize, usize)>> = vec![None; windows.len()];
    for hit in searcher.search_texts(pattern, windows, edit_budget(error_rate, pattern.len())) {
        let entry = &mut best[hit.text_idx];
        if entry.is_none_or(|(cost, start, _)| (hit.cost, hit.text_start) < (cost, start)) {
            *entry = Some((hit.cost, hit.text_start, hit.text_end));
        }
    }
    // Each text runs from the start of its window's hit to `DEGENERATE_SPAN`
    // bases past the hit, plus the edit budget of the longest query.
    let reach = DEGENERATE_SPAN + edit_budget(error_rate, pattern.len() + DEGENERATE_SPAN) + 1;
    let texts: Vec<Vec<u8>> = windows
        .iter()
        .zip(&best)
        .filter_map(|(window, hit)| {
            hit.map(|(_, start, end)| {
                window[start..window.len().min(end + reach)].to_ascii_uppercase()
            })
        })
        .collect();
    if texts.len() < floor {
        return None;
    }
    let texts: Vec<&[u8]> = texts.iter().map(Vec::as_slice).collect();
    let mut layer: Vec<u8> = Vec::new();
    // The base each text holds at each position of the layer.
    let mut observed: Vec<Vec<Option<u8>>> = vec![Vec::new(); texts.len()];
    let mut query = pattern.to_vec();
    for _ in 0..DEGENERATE_SPAN {
        let mut ends: Vec<Option<(i32, usize)>> = vec![None; texts.len()];
        for hit in searcher.search_texts(&query, &texts, edit_budget(error_rate, query.len())) {
            let entry = &mut ends[hit.text_idx];
            if entry.is_none_or(|old| (hit.cost, hit.text_end) < old) {
                *entry = Some((hit.cost, hit.text_end));
            }
        }
        let mut counts = [0usize; 4];
        for ((text, end), seen) in texts.iter().zip(&ends).zip(&mut observed) {
            let base = end.and_then(|(_, end)| text.get(end).copied());
            if let Some(code) = base.and_then(|b| encode_kmer(&[b])) {
                counts[code as usize] += 1;
            }
            seen.push(base);
        }
        if counts.iter().sum::<usize>() < floor {
            break;
        }
        match position_code(counts, composition) {
            Some(code) => {
                layer.push(code);
                query.push(code);
            },
            None => break,
        }
    }
    let degenerate: Vec<usize> = (0..layer.len())
        .filter(|&i| !matches!(layer[i], b'A' | b'C' | b'G' | b'T'))
        .collect();
    let &last = degenerate.last()?;
    if let Some(&spacer) = layer.get(last + 1) {
        let run = layer[last + 1..]
            .iter()
            .take_while(|&&b| b == spacer)
            .count();
        layer.truncate(last + 1 + run);
    }
    if layer.len() < KMER_K {
        return None;
    }
    let mut combinations = std::collections::HashMap::<Vec<u8>, usize>::new();
    for seen in &observed {
        if let Some(key) = degenerate
            .iter()
            .map(|&i| seen[i])
            .collect::<Option<Vec<u8>>>()
        {
            *combinations.entry(key).or_default() += 1;
        }
    }
    let held: usize = combinations.values().sum();
    let recurrent = combinations.values().copied().max().unwrap_or(0);
    tracing::debug!(pattern = %String::from_utf8_lossy(pattern), layer = %String::from_utf8_lossy(&layer),
        held, distinct = combinations.len(), recurrent, "Degenerate layer candidate");
    (held >= floor && recurrent < MIN_SUPPORT_WINDOWS && combinations.len() * 2 >= held)
        .then_some((layer, texts.len()))
}
