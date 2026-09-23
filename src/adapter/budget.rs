//! Edit budgets of adapter hits: the chance-match null model, the per-pattern
//! and set-wide terminal and interior bounds, and the budgets of partial hits.

use super::*;

/// Returns the edit budget for a `len`-base pattern at `rate`, rounded down.
/// The epsilon keeps an integral product whose double lands below its integer.
pub(crate) fn edit_budget(rate: f64, len: usize) -> usize {
    (rate * len as f64 + 1e-9).floor() as usize
}

/// Expected chance interior matches per read, over both strands, that an
/// interior edit budget may admit under an independent uniform base model.
pub(super) const INTERIOR_CHANCE_HITS_PER_READ: f64 = 1e-4;

/// Expected chance terminal matches per read, over both strands and both end
/// zones, that terminal edit budgets may admit under the same null model. The
/// model sums alignment paths and overstates the chance rate, so the bound is
/// conservative. It applies to each pattern alone, which at the default error
/// rate lowers the budget of 11-, 12- and 15-base patterns by one edit, and to
/// the distinct sequences of the set together; see `family_budgets`.
pub(super) const TERMINAL_CHANCE_HITS_PER_READ: f64 = 0.1;

/// Read-length class `c` holds reads shorter than `2^(INTERIOR_CLASS_BITS + c)`
/// bases. Class 0 also holds every shorter read.
pub(super) const INTERIOR_CLASS_BITS: u32 = 12;

/// Number of read-length classes. The last class holds every longer read.
pub(super) const INTERIOR_CLASSES: usize = 20;

/// Returns the read-length class of a read of `read_len` bases.
pub(super) fn interior_class(read_len: usize) -> usize {
    let bits = usize::BITS - read_len.leading_zeros();
    (bits.saturating_sub(INTERIOR_CLASS_BITS) as usize).min(INTERIOR_CLASSES - 1)
}

/// Interior alignment start positions, over both strands, in a read at the
/// ceiling of read-length class `class`.
pub(super) fn interior_positions(class: usize) -> f64 {
    2.0 * 2f64.powi(INTERIOR_CLASS_BITS as i32 + class as i32)
}

/// Returns, for each edit count up to `max_edits`, the probability under the
/// independent uniform DNA null model that `pattern` matches at one position
/// within that many edits. The recurrence sums alignment-path probabilities,
/// including substitutions, insertions and deletions, and therefore
/// overcounts sequences admitting multiple alignments. IUPAC ambiguity
/// increases the probability of a zero-cost match.
pub(super) fn chance_cumulative(pattern: &[u8], max_edits: usize) -> Vec<f64> {
    let mut previous = vec![1.0; max_edits + 1];
    for &base in pattern {
        let p = search::iupac_degeneracy(base).unwrap_or(4) as f64 / 4.0;
        let mut current = vec![0.0; max_edits + 1];
        current[0] = p * previous[0];
        for k in 1..=max_edits {
            current[k] = p * previous[k] + (2.0 - p) * previous[k - 1] + current[k - 1];
        }
        previous = current;
    }
    let mut cumulative = previous;
    for k in 1..=max_edits {
        cumulative[k] += cumulative[k - 1];
    }
    cumulative
}

/// Returns budgets, at most `caps`, under which the distinct sequences among
/// the entries flagged in `included` together admit at most `bound` expected
/// chance hits per read over `positions` alignment start positions. A panel
/// of interchangeable sequences multiplies the chance of a hit by its size,
/// which a bound per pattern does not see. Exact matches are always
/// admitted, so only the chance hits beyond them count. The largest
/// contributors lose one edit at a time, and sequences contributing equally
/// lose it together, so the members of a panel keep one budget. A sequence
/// and its reverse complement, both searched on both strands, are one
/// sequence. Entries not flagged keep their caps.
pub(super) fn family_budgets(
    adapters: &[Adapter],
    included: &[bool],
    caps: &[usize],
    positions: f64,
    bound: f64,
) -> Vec<usize> {
    let mut groups: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
    for (adapter_idx, adapter) in adapters.iter().enumerate() {
        if included[adapter_idx] {
            let forward = adapter.seq.to_ascii_uppercase();
            let reverse = reverse_complement(&forward);
            groups
                .entry(forward.min(reverse))
                .or_default()
                .push(adapter_idx);
        }
    }
    let members: Vec<Vec<usize>> = groups.into_values().collect();
    let mut edits: Vec<usize> = members
        .iter()
        .map(|group| group.iter().map(|&i| caps[i]).min().unwrap_or(0))
        .collect();
    let chance: Vec<Vec<f64>> = members
        .iter()
        .zip(&edits)
        .map(|(group, &k)| {
            chance_cumulative(&adapters[group[0]].seq, k)
                .into_iter()
                .map(|probability| probability * positions)
                .collect()
        })
        .collect();
    let excess = |g: usize, k: usize| chance[g][k] - chance[g][0];
    loop {
        let total: f64 = (0..members.len()).map(|g| excess(g, edits[g])).sum();
        if total <= bound {
            break;
        }
        let top = (0..members.len())
            .filter(|&g| edits[g] > 0)
            .map(|g| chance[g][edits[g]])
            .fold(0.0, f64::max);
        for g in 0..members.len() {
            if edits[g] > 0 && chance[g][edits[g]] >= top * (1.0 - 1e-9) {
                edits[g] -= 1;
            }
        }
    }
    let mut out = caps.to_vec();
    for (group, &k) in members.iter().zip(&edits) {
        for &adapter_idx in group {
            out[adapter_idx] = k;
        }
    }
    out
}

/// Edit budgets of one adapter: the configured terminal tolerance bounded by
/// the chance-match rate over the end zones of the pattern set, and an
/// interior tolerance bounded by the chance-match rate for each read length.
#[derive(Debug, Clone, Copy)]
pub(super) struct Budget {
    /// Pattern length in bases.
    pub(super) len: usize,
    /// Edit budget of the terminal search, which a hit anchored at the read
    /// end or at an accepted hit may use; see `Keep::settle`.
    pub(super) k_end: usize,
    /// Edit budget for a terminal hit anywhere in the end zone. At most
    /// `k_end`.
    pub(super) k_far: usize,
    /// Edit budget for interior hits, per read-length class. Budgets do not
    /// increase with the class.
    pub(super) k_mid: [usize; INTERIOR_CLASSES],
}

impl Budget {
    /// Computes the interior budgets from cumulative chance-match
    /// probabilities under an independent uniform DNA null model. The
    /// recurrence sums alignment-path probabilities, including substitutions,
    /// insertions and deletions, and therefore overcounts sequences admitting
    /// multiple alignments. IUPAC ambiguity increases the probability of a
    /// zero-cost match. Each class admits the largest edit count whose
    /// expected chance matches in a read at the class ceiling stay within
    /// `INTERIOR_CHANCE_HITS_PER_READ`; exact matches are always admitted.
    /// The terminal budget is the configured tolerance, lowered to the
    /// largest edit count whose expected chance hits over both end zones of
    /// `end_size` bases and both strands stay within
    /// `TERMINAL_CHANCE_HITS_PER_READ`.
    pub(super) fn new(pattern: &[u8], error_rate: f64, end_size: usize) -> Self {
        let len = pattern.len();
        let k_end = edit_budget(error_rate, len);
        let cumulative = chance_cumulative(pattern, k_end);
        let k_mid = std::array::from_fn(|class| {
            let positions = interior_positions(class);
            cumulative
                .iter()
                .take_while(|&&probability| {
                    probability * positions <= INTERIOR_CHANCE_HITS_PER_READ
                })
                .count()
                .saturating_sub(1)
        });
        let positions = 4.0 * (end_size + 1) as f64;
        let k_end = cumulative
            .iter()
            .take_while(|&&probability| probability * positions <= TERMINAL_CHANCE_HITS_PER_READ)
            .count()
            .saturating_sub(1)
            .min(k_end);
        Self {
            len,
            k_end,
            k_far: k_end,
            k_mid,
        }
    }

    /// Returns the interior edit budget for a read of `read_len` bases.
    pub(super) fn interior(&self, read_len: usize) -> usize {
        self.k_mid[interior_class(read_len)]
    }

    /// Returns the largest interior edit budget over all read lengths.
    pub(super) fn interior_max(&self) -> usize {
        self.k_mid[0]
    }
}

/// Returns the overhang cost the overhang searcher charges a hit at per-base
/// rate `alpha`: `floor(alpha * bases)` for each overhanging side, at the
/// searcher's `f32` precision.
pub(super) fn overhang_cost(alpha: f32, left: usize, right: usize) -> usize {
    let side = |bases: usize| (bases as f32 * alpha).floor() as usize;
    side(left) + side(right)
}

/// Returns the edit budget of a partial hit whose `overlap` bases aligned
/// inside the read. The first `MIN_OVERLAP` bases must match exactly and the
/// rate applies to the remainder, so the shortest accepted overlaps carry no
/// tolerance and a random read end is not mistaken for adapter residue.
pub(super) fn partial_budget(rate: f64, overlap: usize) -> usize {
    edit_budget(rate, overlap.saturating_sub(MIN_OVERLAP - 1))
}
