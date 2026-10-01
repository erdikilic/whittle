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
/// the chance families of the set together; see `family_budgets`.
pub(crate) const TERMINAL_CHANCE_HITS_PER_READ: f64 = 0.1;

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

/// Returns budgets, at most `caps`, under which the chance families flagged
/// in `included` together admit at most `bound` expected chance hits per
/// read over `positions` start positions; a panel multiplies the chance rate
/// by its size. Only chance hits beyond exact matches count. The largest
/// contributors lose one edit at a time, equal ones together. Unflagged
/// entries keep their caps. See `chance_families` for what one family is.
pub(super) fn family_budgets(
    adapters: &[Adapter],
    included: &[bool],
    sites: &[Option<usize>],
    caps: &[usize],
    positions: f64,
    bound: f64,
) -> Vec<usize> {
    let families = chance_families(adapters, included, sites);
    let mut edits = family_caps(&families, caps);
    let chance: Vec<Vec<Vec<f64>>> = families
        .iter()
        .zip(&edits)
        .map(|(family, &k)| {
            distinct_members(adapters, family)
                .into_iter()
                .map(|seq| {
                    chance_cumulative(seq, k)
                        .into_iter()
                        .map(|probability| probability * positions)
                        .collect()
                })
                .collect()
        })
        .collect();
    let excess = |f: usize, k: usize| {
        chance[f]
            .iter()
            .map(|member| member[k] - member[0])
            .fold(0.0, f64::max)
    };
    let rate = |f: usize, k: usize| chance[f].iter().map(|member| member[k]).fold(0.0, f64::max);
    loop {
        let total: f64 = (0..families.len()).map(|f| excess(f, edits[f])).sum();
        if total <= bound {
            break;
        }
        let top = (0..families.len())
            .filter(|&f| edits[f] > 0)
            .map(|f| rate(f, edits[f]))
            .fold(0.0, f64::max);
        for (f, k) in edits.iter_mut().enumerate() {
            if *k > 0 && rate(f, *k) >= top * (1.0 - 1e-9) {
                *k -= 1;
            }
        }
    }
    spread(&families, &edits, caps)
}

/// Returns the entries flagged in `included`, grouped into chance families:
/// the entries of one sequence, a sequence and its reverse complement
/// together, and the entries that bind one marker-primer site (`sites`),
/// whose variants match the same template bases. A family takes one budget,
/// the lowest cap of its members, and admits the chance hits of its member
/// that admits the most. The bounds of the other entries thus see a site as
/// one sequence however many variants a kit lists.
fn chance_families(
    adapters: &[Adapter],
    included: &[bool],
    sites: &[Option<usize>],
) -> Vec<Vec<usize>> {
    let mut families: BTreeMap<(Option<usize>, Vec<u8>), Vec<usize>> = BTreeMap::new();
    for (adapter_idx, adapter) in adapters.iter().enumerate() {
        if included[adapter_idx] {
            let key = match sites.get(adapter_idx).copied().flatten() {
                Some(site) => (Some(site), Vec::new()),
                None => {
                    let forward = adapter.seq.to_ascii_uppercase();
                    let reverse = reverse_complement(&forward);
                    (None, forward.min(reverse))
                },
            };
            families.entry(key).or_default().push(adapter_idx);
        }
    }
    families.into_values().collect()
}

/// Returns the lowest cap of each family's members.
fn family_caps(families: &[Vec<usize>], caps: &[usize]) -> Vec<usize> {
    families
        .iter()
        .map(|family| family.iter().map(|&i| caps[i]).min().unwrap_or(0))
        .collect()
}

/// Returns the distinct sequences of a family's members, a sequence and its
/// reverse complement once.
fn distinct_members<'a>(adapters: &'a [Adapter], family: &[usize]) -> Vec<&'a [u8]> {
    let mut seen: Vec<Vec<u8>> = Vec::new();
    let mut out = Vec::new();
    for &adapter_idx in family {
        let forward = adapters[adapter_idx].seq.to_ascii_uppercase();
        let canonical = reverse_complement(&forward).min(forward);
        if !seen.contains(&canonical) {
            seen.push(canonical);
            out.push(adapters[adapter_idx].seq.as_slice());
        }
    }
    out
}

/// Returns `caps` with the budget of each family in `edits` given to its
/// members.
fn spread(families: &[Vec<usize>], edits: &[usize], caps: &[usize]) -> Vec<usize> {
    let mut out = caps.to_vec();
    for (family, &k) in families.iter().zip(edits) {
        for &adapter_idx in family {
            out[adapter_idx] = k;
        }
    }
    out
}

/// Start positions of the second primer of a junction pair relative to the
/// end of the first: it may begin up to `FLANK_SLACK` bases after it or
/// overlap it by as many.
pub(super) const PAIR_OFFSETS: usize = 2 * FLANK_SLACK + 1;

/// Returns budgets, at most `caps`, under which chance junction pairs of the
/// entries flagged in `included` stay within `bound` expected per read over
/// `positions` start positions. A pair is any flagged sequence followed, at
/// one of `PAIR_OFFSETS` offsets, by any flagged sequence, so the chance rate
/// of a pair at one position is the square of the summed chance rates of the
/// chance families (`chance_families`), times the offsets. The largest
/// contributors lose one edit at a time, equal ones together, until the
/// bound holds or no edits remain. Unflagged entries keep their caps.
pub(super) fn pair_budgets(
    adapters: &[Adapter],
    included: &[bool],
    sites: &[Option<usize>],
    caps: &[usize],
    positions: f64,
    bound: f64,
) -> Vec<usize> {
    let families = chance_families(adapters, included, sites);
    let mut edits = family_caps(&families, caps);
    let chance: Vec<Vec<Vec<f64>>> = families
        .iter()
        .zip(&edits)
        .map(|(family, &k)| {
            distinct_members(adapters, family)
                .into_iter()
                .map(|seq| chance_cumulative(seq, k))
                .collect()
        })
        .collect();
    let rate = |f: usize, k: usize| chance[f].iter().map(|member| member[k]).fold(0.0, f64::max);
    let scale = positions * PAIR_OFFSETS as f64;
    loop {
        let total: f64 = (0..families.len()).map(|f| rate(f, edits[f])).sum();
        if total * total * scale <= bound {
            break;
        }
        let top = (0..families.len())
            .filter(|&f| edits[f] > 0)
            .map(|f| rate(f, edits[f]))
            .fold(0.0, f64::max);
        if top == 0.0 {
            break;
        }
        for (f, k) in edits.iter_mut().enumerate() {
            if *k > 0 && rate(f, *k) >= top * (1.0 - 1e-9) {
                *k -= 1;
            }
        }
    }
    spread(&families, &edits, caps)
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
    /// Edit budget of a whole hit of a split sheet primer anchored at an
    /// outer layer's trim boundary: the largest edit count, at most the
    /// error-rate ceiling, whose expected chance hits over
    /// `ANCHORED_POSITIONS` stay within `TERMINAL_CHANCE_HITS_PER_READ`,
    /// for the pattern alone and for the split entries of the set together
    /// (`set_budgets`).
    pub(super) k_anchor: usize,
    /// Edit budget for interior hits, per read-length class. Budgets do not
    /// increase with the class.
    pub(super) k_mid: [usize; INTERIOR_CLASSES],
    /// Edit budget for an interior hit of a `CandidateIndex::paired` entry
    /// beside its junction partner, per read-length class; see
    /// `pair_budgets`. Zero for other entries. Budgets do not increase with
    /// the class.
    pub(super) k_pair: [usize; INTERIOR_CLASSES],
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
    /// `TERMINAL_CHANCE_HITS_PER_READ`. The anchored budget bounds the same
    /// rate over the positions beside the trim boundaries of one read
    /// (`ANCHORED_POSITIONS`).
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
        let k_end = terminal_edits(&cumulative, 4.0 * (end_size + 1) as f64);
        Self {
            len,
            k_end,
            k_far: k_end,
            k_anchor: terminal_edits(&cumulative, ANCHORED_POSITIONS),
            k_mid,
            k_pair: [0; INTERIOR_CLASSES],
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

    /// Returns the pair budget for a read of `read_len` bases.
    pub(super) fn pair(&self, read_len: usize) -> usize {
        self.k_pair[interior_class(read_len)]
    }

    /// Returns the lower of each budget of `self` and `other`, two budgets
    /// of one pattern.
    pub(super) fn stricter(&self, other: &Budget) -> Budget {
        Budget {
            len: self.len,
            k_end: self.k_end.min(other.k_end),
            k_far: self.k_far.min(other.k_far),
            k_anchor: self.k_anchor.min(other.k_anchor),
            k_mid: std::array::from_fn(|c| self.k_mid[c].min(other.k_mid[c])),
            k_pair: std::array::from_fn(|c| self.k_pair[c].min(other.k_pair[c])),
        }
    }
}

/// Returns the terminal edit budget of `pattern` at `error_rate`, bounded by
/// the chance-match rate over both end zones of `end_size` bases; see
/// `Budget::new`.
pub(crate) fn terminal_budget(pattern: &[u8], error_rate: f64, end_size: usize) -> usize {
    Budget::new(pattern, error_rate, end_size).k_end
}

/// Bases by which the outer edge of a whole hit anchored at an outer layer's
/// trim boundary may lie outboard of that boundary. The bases outboard of the
/// boundary are the flank the layer matched, not random bases, so a primer
/// extends into them only by the few bases the layer's alignment may claim
/// from the primer at their junction.
pub(super) const BOUNDARY_OUTER_SLACK: usize = 3;

/// Start positions, per read, of a whole hit anchored at an outer layer's
/// trim boundary: its outer edge lies at most `BOUNDARY_OUTER_SLACK` bases
/// outboard of the boundary or at most `FLANK_SLACK` bases inboard of it, at
/// each of the two ends, on the one strand that faces the insert there.
pub(super) const ANCHORED_POSITIONS: f64 = 2.0 * (FLANK_SLACK + BOUNDARY_OUTER_SLACK + 1) as f64;

/// Returns the largest edit count, at most the last of `cumulative`, whose
/// expected chance hits over `positions` start positions stay within
/// `TERMINAL_CHANCE_HITS_PER_READ`; exact matches are always admitted.
fn terminal_edits(cumulative: &[f64], positions: f64) -> usize {
    cumulative
        .iter()
        .take_while(|&&probability| probability * positions <= TERMINAL_CHANCE_HITS_PER_READ)
        .count()
        .saturating_sub(1)
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

/// Returns whether a hit whose `overhang` sides cost `overhang_cost(rate as
/// f32, left_overhang, right_overhang)` stays, once that charge is
/// discounted from `cost`, within the partial budget of its `overlap`
/// aligned bases. An overhanging hit must clear this in addition to its
/// whole-pattern terminal budget, or a hit could pass on overhang discount
/// alone with an unsupported aligned portion.
pub(crate) fn residue_within_budget(
    rate: f64,
    cost: usize,
    left_overhang: usize,
    right_overhang: usize,
    overlap: usize,
) -> bool {
    let charged = overhang_cost(rate as f32, left_overhang, right_overhang);
    cost.saturating_sub(charged) <= partial_budget(rate, overlap)
}
