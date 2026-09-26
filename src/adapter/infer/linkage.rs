//! The insert boundary of an amplicon from variation that follows the
//! template.
//!
//! Every read of an amplicon library carries the same primer whatever
//! template it was copied from, while the sequence behind the primer is the
//! template's own. The conserved start that related templates share still
//! holds variable bases, and each variant is carried by the templates of one
//! lineage, whose sequence further inside the insert differs from that of the
//! other lineages. A variable base of a technical sequence, such as a
//! degenerate primer base or a recurrent sequencing error, is not tied to the
//! template. A column whose two most frequent bases are each identified by
//! k-mers far inside the insert, and whose reads of either base continue
//! alike, therefore lies in the insert, however conserved the column is.
//! Alternative technical sequences at one position, such as the primers of
//! the two strands, are identified by what follows them as well, but
//! continue differently.

use super::*;

/// Length of the k-mers that identify the template far inside the insert.
const FAR_K: usize = 12;

/// Bases of the insert, behind the tested columns, whose k-mers identify the
/// template.
const FAR_LEN: usize = 100;

/// Columns read on each side of a candidate: outboard of its outer edge and
/// past its inner end. The far k-mers begin behind the inner ones.
const FLANK: usize = 16;

/// Divisor of a column's aligned windows giving the fewest windows that its
/// second base must hold for the column to be tested.
const MINOR_BASE_DIVISOR: usize = 20;

/// Percentage of the windows holding a far k-mer, among the windows of the
/// two tested bases, that must hold one base for the k-mer to identify it.
const IDENTIFYING_PERCENT: usize = 90;

/// Fraction of the windows of each tested base that must hold a far k-mer
/// identifying that base.
const IDENTIFIED_FRACTION: f64 = 0.25;

/// Divisor of the windows behind a boundary giving the fewest windows of a
/// family of them that is tested on its own.
const FAMILY_DIVISOR: usize = 4;

/// Where the template of an amplicon begins relative to a candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum TemplateStart {
    /// A column this many bases outboard of the candidate's outer edge
    /// follows the template: the candidate lies in the insert.
    Before(usize),
    /// The column at this distance from the candidate's outer end, inside
    /// the candidate, follows the template; the outer part is technical.
    Inside(usize),
    /// The candidate is technical and the template begins behind it.
    Behind,
}

/// The bases that the windows holding a pattern read at the columns of the
/// pattern and at `FLANK` columns on either side, with the far k-mers of
/// each window, in inward orientation: from the physical read end toward
/// the insert. Column `FLANK` is the pattern's outer base.
struct Columns {
    /// Per window with a hit: the base code of each column, or 4 for a gap
    /// or a column outside the window.
    bases: Vec<Vec<u8>>,
    /// Per window: the distinct codes of the `FAR_K`-mers behind the inner
    /// flank, sorted.
    far: Vec<Vec<u64>>,
    /// Per window: the distance of the hit's outer edge from the window's
    /// outer end.
    outer: Vec<usize>,
}

/// Returns `seq` read inward from `end`: unchanged for the 5' end, reversed
/// for the 3' end, where the insert lies before the sequence.
fn inward(seq: &[u8], end: End) -> Vec<u8> {
    match end {
        End::Five => seq.to_vec(),
        End::Three => seq.iter().rev().copied().collect(),
    }
}

/// Returns the code of `base`, or 4 for an ambiguity code or a gap.
fn base_code(base: Option<&u8>) -> u8 {
    base.and_then(|b| encode_kmer(std::slice::from_ref(b)))
        .map_or(4, |code| code as u8)
}

/// Returns the sorted distinct codes of the `FAR_K`-mers that start in
/// `window[from..from + FAR_LEN]`.
fn far_kmers(window: &[u8], from: usize) -> Vec<u64> {
    let stop = (from + FAR_LEN).min(window.len().saturating_sub(FAR_K - 1));
    let mut codes: Vec<u64> = (from..stop)
        .filter_map(|start| encode_kmer(&window[start..start + FAR_K]))
        .collect();
    codes.sort_unstable();
    codes.dedup();
    codes
}

impl Columns {
    /// Aligns `seq` to each window of `end` within `edits` and reads the
    /// bases of its columns and of its flanks.
    fn of_candidate(seq: &[u8], windows: &[&[u8]], end: End, edits: usize) -> Self {
        let pattern = inward(seq, end);
        let texts: Vec<Vec<u8>> = windows.iter().map(|w| inward(w, end)).collect();
        let texts: Vec<&[u8]> = texts.iter().map(Vec::as_slice).collect();
        let mut searcher = crate::adapter::search::new_searcher_fwd();
        let mut best: Vec<Option<sassy::Match>> = vec![None; texts.len()];
        for hit in searcher.search_texts(&pattern, &texts, edits) {
            let entry = &mut best[hit.text_idx];
            if entry
                .as_ref()
                .is_none_or(|old| (hit.cost, hit.text_start) < (old.cost, old.text_start))
            {
                *entry = Some(hit.clone());
            }
        }
        let inner = FLANK + pattern.len();
        let mut columns = Self {
            bases: Vec::new(),
            far: Vec::new(),
            outer: Vec::new(),
        };
        for (text, hit) in texts.iter().zip(best) {
            let Some(hit) = hit else { continue };
            let mut bases = vec![4u8; inner + FLANK];
            for (offset, base) in bases[..FLANK].iter_mut().rev().enumerate() {
                *base = base_code(
                    hit.text_start
                        .checked_sub(offset + 1)
                        .and_then(|at| text.get(at)),
                );
            }
            let path = hit.to_path();
            for (i, pos) in path.iter().enumerate() {
                let next = path
                    .get(i + 1)
                    .map(|p| (p.0, p.1))
                    .unwrap_or((pattern.len() as i32, hit.text_end as i32));
                if next.0 == pos.0 + 1 && next.1 == pos.1 + 1 {
                    bases[FLANK + pos.0 as usize] = base_code(text.get(pos.1 as usize));
                }
            }
            for (offset, base) in bases[inner..].iter_mut().enumerate() {
                *base = base_code(text.get(hit.text_end + offset));
            }
            columns.bases.push(bases);
            columns.far.push(far_kmers(text, hit.text_end + FLANK));
            columns.outer.push(hit.text_start);
        }
        columns
    }

    /// Reads the first `FLANK` bases of each window of `end` as the columns
    /// behind a boundary, the far k-mers behind them.
    fn at_boundary(windows: &[&[u8]], end: End) -> Self {
        let mut columns = Self {
            bases: Vec::new(),
            far: Vec::new(),
            outer: Vec::new(),
        };
        for window in windows {
            let text = inward(window, end);
            columns
                .bases
                .push((0..FLANK).map(|at| base_code(text.get(at))).collect());
            columns.far.push(far_kmers(&text, FLANK));
            columns.outer.push(0);
        }
        columns
    }

    /// Returns the median distance of the hits' outer edges from the outer
    /// end of their windows, or `usize::MAX` without hits.
    fn median_outer(&self) -> usize {
        let mut outer = self.outer.clone();
        outer.sort_unstable();
        outer.get(outer.len() / 2).copied().unwrap_or(usize::MAX)
    }

    /// Returns the counts of A, C, G and T at `column`.
    fn counts(&self, column: usize) -> [usize; 4] {
        let mut counts = [0usize; 4];
        for bases in &self.bases {
            if let Some(&b) = bases.get(column).filter(|&&b| b < 4) {
                counts[b as usize] += 1;
            }
        }
        counts
    }

    /// Returns whether `column` follows the template: its two most frequent
    /// bases each hold far k-mers that identify them in
    /// `IDENTIFIED_FRACTION` of their windows, and the windows of the two
    /// continue alike (`continue_alike`). The second base must hold at least
    /// `MIN_SUPPORT_WINDOWS` windows and one in `MINOR_BASE_DIVISOR` of the
    /// column.
    fn follows_template(&self, column: usize) -> bool {
        let counts = self.counts(column);
        let total: usize = counts.iter().sum();
        let mut ranked = [0usize, 1, 2, 3];
        ranked.sort_by_key(|&b| (std::cmp::Reverse(counts[b]), b));
        let (first, second) = (ranked[0] as u8, ranked[1] as u8);
        if counts[second as usize] < MIN_SUPPORT_WINDOWS.max(total / MINOR_BASE_DIVISOR) {
            return false;
        }
        let allele = |bases: &[u8]| match bases[column] {
            b if b == first => Some(0),
            b if b == second => Some(1),
            _ => None,
        };
        let mut presence: KmerMap<[u32; 2]> = KmerMap::default();
        for (bases, far) in self.bases.iter().zip(&self.far) {
            if let Some(a) = allele(bases) {
                for &code in far {
                    presence.entry(code).or_default()[a] += 1;
                }
            }
        }
        let identifies = |code: &u64, a: usize| {
            presence.get(code).is_some_and(|held| {
                let n = (held[0] + held[1]) as usize;
                n >= MIN_SUPPORT_WINDOWS && held[a] as usize * 100 >= n * IDENTIFYING_PERCENT
            })
        };
        let mut identified = [0usize; 2];
        let mut windows = [0usize; 2];
        for (bases, far) in self.bases.iter().zip(&self.far) {
            if let Some(a) = allele(bases) {
                windows[a] += 1;
                identified[a] += usize::from(far.iter().any(|code| identifies(code, a)));
            }
        }
        (0..2).all(|a| identified[a] as f64 >= IDENTIFIED_FRACTION * windows[a] as f64)
            && self.continue_alike(column, first, second)
    }

    /// Returns whether the windows holding `first` and those holding
    /// `second` at `column` continue alike: the most frequent bases of the
    /// two groups agree in at least two thirds of the `FLANK` columns after it
    /// that both groups read in enough windows. Variants of a conserved
    /// insert are followed by the same conserved sequence; alternative
    /// technical sequences, such as the primers of the two strands or the
    /// members of a barcode panel, differ throughout.
    fn continue_alike(&self, column: usize, first: u8, second: u8) -> bool {
        let width = self.bases.first().map_or(0, Vec::len);
        let span = column + 1..(column + 1 + FLANK).min(width);
        let mut tallies = vec![[[0usize; 4]; 2]; span.len()];
        for bases in &self.bases {
            let group = match bases[column] {
                b if b == first => 0,
                b if b == second => 1,
                _ => continue,
            };
            for (tally, &b) in tallies.iter_mut().zip(&bases[span.clone()]) {
                if b < 4 {
                    tally[group][b as usize] += 1;
                }
            }
        }
        let top = |counts: &[usize; 4]| (0..4).max_by_key(|&b| (counts[b], std::cmp::Reverse(b)));
        let (mut read, mut agree) = (0, 0);
        for [a, b] in &tallies {
            let floor = MIN_SUPPORT_WINDOWS / 2;
            if a.iter().sum::<usize>() >= floor && b.iter().sum::<usize>() >= floor {
                read += 1;
                agree += usize::from(top(a) == top(b));
            }
        }
        read > 0 && 3 * agree >= 2 * read
    }
}

/// Places the start of the template of an amplicon relative to `seq`, a
/// candidate of `end` held by `windows`, or returns `None` when no column
/// follows the template, with the median depth of the candidate's hits from
/// the window ends. The `FLANK` columns outboard of the candidate are
/// tested, then every column of the candidate, then the first column past
/// its inner end, so that the candidate ends where the template begins.
pub(super) fn template_start(
    seq: &[u8],
    windows: &[&[u8]],
    end: End,
    error_rate: f64,
) -> (Option<TemplateStart>, usize) {
    let columns = Columns::of_candidate(seq, windows, end, edit_budget(error_rate, seq.len()));
    let depth = columns.median_outer();
    if columns.bases.len() < MIN_SUPPORT_WINDOWS {
        return (None, depth);
    }
    let inner = FLANK + seq.len();
    if let Some(column) = (0..FLANK).find(|&c| columns.follows_template(c)) {
        return (Some(TemplateStart::Before(FLANK - column)), depth);
    }
    if let Some(column) = (FLANK..inner).find(|&c| columns.follows_template(c)) {
        return (Some(TemplateStart::Inside(column - FLANK)), depth);
    }
    let behind = columns.follows_template(inner);
    (behind.then_some(TemplateStart::Behind), depth)
}

/// Returns, per window of `end`, the span of the best hit of `seq` within
/// `edits` as distances of its outer edge and its inner end from the
/// window's outer end, or `None` without a hit.
pub(super) fn hit_spans(
    seq: &[u8],
    windows: &[&[u8]],
    end: End,
    edits: usize,
) -> Vec<Option<(usize, usize)>> {
    let pattern = inward(seq, end);
    let texts: Vec<Vec<u8>> = windows.iter().map(|w| inward(w, end)).collect();
    let texts: Vec<&[u8]> = texts.iter().map(Vec::as_slice).collect();
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let mut best: Vec<Option<(i32, usize, usize)>> = vec![None; texts.len()];
    for hit in searcher.search_texts(&pattern, &texts, edits) {
        let entry = &mut best[hit.text_idx];
        if entry.is_none_or(|(cost, _, _)| hit.cost < cost) {
            *entry = Some((hit.cost, hit.text_start, hit.text_end));
        }
    }
    best.into_iter()
        .map(|hit| hit.map(|(_, outer, inner)| (outer, inner)))
        .collect()
}

/// Returns the bases by which a candidate with hits `spans` reads past
/// `primer`, the hit spans of a candidate that ends at the insert boundary,
/// in the windows that hold both: the median excess of the candidate's inner
/// end over the primer's, when it is positive in most of at least
/// `MIN_SUPPORT_WINDOWS` shared windows, or zero.
pub(super) fn read_through(
    spans: &[Option<(usize, usize)>],
    primer: &[Option<(usize, usize)>],
) -> usize {
    let mut excess: Vec<isize> = spans
        .iter()
        .zip(primer)
        .filter_map(|(own, primer)| {
            Some(own.as_ref()?.1.cast_signed() - primer.as_ref()?.1.cast_signed())
        })
        .collect();
    if excess.len() < MIN_SUPPORT_WINDOWS {
        return 0;
    }
    excess.sort_unstable();
    let median = excess[excess.len() / 2];
    if 2 * excess.iter().filter(|&&e| e > 0).count() > excess.len() {
        median.max(0).cast_unsigned()
    } else {
        0
    }
}

/// Returns `seq` continued to the inner end of `primer` when it ends inside
/// it, or `None`. `spans` and `primer_spans` are the hit spans of the two
/// candidates of `end`; `seq` ends inside the primer when, in most of at
/// least `MIN_SUPPORT_WINDOWS` windows holding both, its inner end lies past
/// the primer's outer edge and short of the primer's inner end. The
/// primer's bases past the median depth of that inner end within the
/// primer continue `seq`, which then ends at the insert boundary as the
/// primer does.
pub(super) fn continued_into(
    seq: &[u8],
    primer: &[u8],
    spans: &[Option<(usize, usize)>],
    primer_spans: &[Option<(usize, usize)>],
    end: End,
) -> Option<Vec<u8>> {
    let mut shared = 0;
    let mut offsets = Vec::new();
    for (own, other) in spans.iter().zip(primer_spans) {
        if let (Some((_, inner)), Some((outer, stop))) = (own, other) {
            shared += 1;
            if inner > outer && inner < stop {
                offsets.push(inner - outer);
            }
        }
    }
    if shared < MIN_SUPPORT_WINDOWS || 2 * offsets.len() <= shared {
        return None;
    }
    offsets.sort_unstable();
    let offset = offsets[offsets.len() / 2];
    let primer = inward(primer, end);
    let mut joined = inward(seq, end);
    joined.extend_from_slice(primer.get(offset..)?);
    Some(inward(&joined, end))
}

/// Returns whether the template of an amplicon begins at the boundary of
/// `windows`, which read inward from the boundary of `end`: one of the first
/// `MIN_PATTERN_LEN` columns follows the template, so that no technical
/// layer long enough to be discovered lies before it. Windows that read
/// different regions of the templates behind the boundary, such as the two
/// ends of a gene, share no conserved sequence to hold a column; each
/// family, split by the most frequent far k-mer and holding at least
/// `1 / FAMILY_DIVISOR` of the windows, is tested as well.
pub(super) fn template_at_boundary(windows: &[&[u8]], end: End) -> bool {
    if windows.len() < MIN_SUPPORT_WINDOWS {
        return false;
    }
    let columns = Columns::at_boundary(windows, end);
    if (0..MIN_PATTERN_LEN).any(|c| columns.follows_template(c)) {
        return true;
    }
    let mut counts: KmerMap<u32> = KmerMap::default();
    for far in &columns.far {
        for &code in far {
            *counts.entry(code).or_default() += 1;
        }
    }
    let Some(common) = counts
        .iter()
        .max_by_key(|&(code, count)| (*count, std::cmp::Reverse(*code)))
        .map(|(&code, _)| code)
    else {
        return false;
    };
    let mut families: [Vec<&[u8]>; 2] = [Vec::new(), Vec::new()];
    for (&window, far) in windows.iter().zip(&columns.far) {
        families[usize::from(far.binary_search(&common).is_err())].push(window);
    }
    families.iter().any(|family| {
        family.len() >= MIN_SUPPORT_WINDOWS.max(windows.len() / FAMILY_DIVISOR) && {
            let columns = Columns::at_boundary(family, end);
            (0..MIN_PATTERN_LEN).any(|c| columns.follows_template(c))
        }
    })
}

/// Fewest outer bases of a candidate that must repeat the inner end of a
/// sequence ending at the insert boundary for the candidate to read on past
/// that end, as a read eroded into the primer at its physical end does.
const OVERLAP_MIN: usize = 7;

/// Returns whether `base` is one of the bases of the IUPAC code `code`.
fn compatible(base: u8, code: u8) -> bool {
    let bases = |c| crate::adapter::search::iupac_bases(c).unwrap_or(&[]);
    bases(base).iter().all(|b| bases(code).contains(b))
}

/// Returns the number of outer bases of `seq`, a candidate of `end`, that
/// lie outboard of the inner end of `boundary`, a sequence of `end` whose
/// inner end is the insert boundary, when `seq` reads on past that end: it
/// holds the last `MIN_PATTERN_LEN` bases of `boundary` before its own inner
/// end, or its outer `OVERLAP_MIN` or more bases repeat the inner end of
/// `boundary`. `None` otherwise.
pub(super) fn past_boundary(seq: &[u8], boundary: &[u8], end: End) -> Option<usize> {
    let (seq, boundary) = (inward(seq, end), inward(boundary, end));
    let matches = |s: &[u8], b: &[u8]| s.iter().zip(b).all(|(&x, &y)| compatible(x, y));
    if boundary.len() >= MIN_PATTERN_LEN {
        let tail = &boundary[boundary.len() - MIN_PATTERN_LEN..];
        if let Some(at) = seq.windows(tail.len()).position(|w| matches(w, tail)) {
            return Some(at + tail.len()).filter(|&keep| keep < seq.len());
        }
    }
    (OVERLAP_MIN..seq.len().min(boundary.len() + 1))
        .rev()
        .find(|&n| matches(&seq[..n], &boundary[boundary.len() - n..]))
}

/// Ends each candidate of `end` that reads on past the inner end of one of
/// `boundaries` (`past_boundary`) at that end, marking it as ending at the
/// insert boundary, and drops the candidates left shorter than
/// `MIN_PATTERN_LEN`.
pub(super) fn end_at_insert_boundaries(
    candidates: &mut Vec<Candidate>,
    boundaries: &[Vec<u8>],
    end: End,
) {
    for (seq, _, boundary, _, unbounded, _) in candidates.iter_mut() {
        let keep = boundaries
            .iter()
            .filter(|b| *b != seq)
            .filter_map(|b| past_boundary(seq, b, end))
            .min();
        if let Some(keep) = keep {
            tracing::debug!(sequence = %String::from_utf8_lossy(seq), keep, "Candidate reads past an insert boundary");
            *seq = outer_part(seq, keep, end);
            *boundary = true;
            *unbounded = false;
        }
    }
    candidates.retain(|(seq, _, _, _, _, _)| seq.len() >= MIN_PATTERN_LEN);
}

/// Returns the `keep` bases of `seq` nearest the physical `end`.
pub(super) fn outer_part(seq: &[u8], keep: usize, end: End) -> Vec<u8> {
    match end {
        End::Five => seq[..keep].to_vec(),
        End::Three => seq[seq.len() - keep..].to_vec(),
    }
}

/// Returns the number of bases of `seq`, a candidate of `end`, that lie
/// outboard of the inner end of a marker-gene primer it holds reading into
/// the insert, or `None`. The primer's inner `KMER_K` bases must align
/// within their edit budget.
pub(super) fn marker_primer_end(seq: &[u8], end: End, error_rate: f64) -> Option<usize> {
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let seq = seq.to_ascii_uppercase();
    crate::adapter::catalog::MARKER_PRIMERS
        .iter()
        .filter_map(|&primer| {
            let pattern = match end {
                End::Five => primer[primer.len() - KMER_K..].to_vec(),
                End::Three => crate::adapter::reverse_complement(primer)[..KMER_K].to_vec(),
            };
            hits(
                &mut searcher,
                &pattern,
                &seq,
                edit_budget(error_rate, KMER_K),
            )
            .into_iter()
            .min_by_key(|hit| hit.cost)
            .map(|hit| match end {
                End::Five => hit.end,
                End::Three => seq.len() - hit.start,
            })
        })
        .filter(|&keep| keep >= MIN_PATTERN_LEN)
        .min()
}
