//! Primer split classifier: assigns a read end pair to a sheet target from
//! the primers rescored there. Pure logic, no I/O and no engine types; the
//! engine hands it the rescoring output (`Score`) for each end.

use std::collections::HashMap;

use super::sheet::{Require, Sheet};

/// The granularity a key names: one key per target, or one key per group
/// (several targets sharing a group collapse into one key).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum KeyLevel {
    /// Each target is its own key.
    Target,
    /// Targets sharing a `group` collapse into one key.
    Group,
}

/// The keys a sheet classifies into at a given `KeyLevel`, and the key each
/// target belongs to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Keys {
    /// Key names, in first-seen sheet order.
    pub names: Vec<String>,
    /// `of_target[i]` is the index into `names` of `sheet.targets[i]`'s key.
    pub of_target: Vec<usize>,
}

impl Keys {
    /// Builds the keys of `sheet` at `level`. At `KeyLevel::Target`, `names`
    /// is the target names in sheet order. At `KeyLevel::Group`, `names` is
    /// the distinct group names in first-seen order, and several targets
    /// sharing a group point at the same key.
    pub fn new(sheet: &Sheet, level: KeyLevel) -> Keys {
        match level {
            KeyLevel::Target => Keys {
                names: sheet.targets.iter().map(|t| t.name.clone()).collect(),
                of_target: (0..sheet.targets.len()).collect(),
            },
            KeyLevel::Group => {
                let mut names: Vec<String> = Vec::new();
                let mut index: HashMap<&str, usize> = HashMap::new();
                let of_target = sheet
                    .targets
                    .iter()
                    .map(|t| {
                        *index.entry(t.group.as_str()).or_insert_with(|| {
                            names.push(t.group.clone());
                            names.len() - 1
                        })
                    })
                    .collect();
                Keys { names, of_target }
            },
        }
    }
}

/// One rescored primer at a read end: `primer` indexes `Sheet::primers`, and
/// `cost` is its alignment cost there. When the same primer appears twice at
/// one end, the lowest cost counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Score {
    /// Index into `Sheet::primers`.
    pub primer: usize,
    /// Alignment cost of this primer at this end.
    pub cost: usize,
}

/// The read strand a target was found on: which physical end (5' or 3') the
/// target's forward primer sits at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// The forward primer sits at the 5' end, the reverse primer at 3'.
    Plus,
    /// The reverse primer sits at the 5' end, the forward primer at 3'.
    Minus,
}

/// Which read end(s) carried the primer evidence behind a call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ends {
    /// Both the 5' and 3' ends had a score.
    Both,
    /// Only the 5' end had a score.
    Five,
    /// Only the 3' end had a score.
    Three,
}

/// Why a segment was not assigned to a key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unassigned {
    /// Neither end had a score.
    NoPrimer,
    /// A target matched, but not at the end(s) `Rules::require` needs.
    Require,
    /// A target's primer scored at both ends, so no strand is consistent.
    Orientation,
    /// The assigned target's `len` window excludes the segment length.
    Length,
}

impl Unassigned {
    /// The wording used for this reason in the report and the `wt:Z` tag.
    pub fn label(self) -> &'static str {
        match self {
            Unassigned::NoPrimer => "no_primer",
            Unassigned::Require => "require",
            Unassigned::Orientation => "orientation",
            Unassigned::Length => "length",
        }
    }
}

/// The classifier's outcome for one segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Call {
    /// Assigned to `key`, whose cheapest consistent target is `target`.
    Assigned {
        /// Index into `Keys::names`.
        key: usize,
        /// Index into `Sheet::targets`: the cheapest of `key`'s consistent
        /// targets. Equal ranks prefer a compatible length window when a
        /// length is supplied, then the lower target index.
        target: usize,
        /// The strand the target was found on.
        strand: Strand,
        /// The end(s) that carried the primer evidence.
        ends: Ends,
    },
    /// Not assigned, for the given reason.
    Unassigned(Unassigned),
    /// Two or more keys tied for cheapest, or none beat the runner-up by
    /// `Rules::lead`.
    Ambiguous,
}

/// The end rule and lead margin `classify` applies.
#[derive(Clone, Copy, Debug)]
pub struct Rules {
    /// Which end(s) a target needs a located primer at.
    pub require: Require,
    /// Minimum cost lead the best key needs over the best different key to
    /// be assigned rather than ambiguous.
    pub lead: usize,
}

/// Facts `classify` weighs a used end with besides its scores: the keys
/// that list each primer, its close primers of other keys, and the penalty
/// of a role primer that did not score.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bounds {
    /// `keys[p]` lists the keys whose targets name primer `p`.
    keys: Vec<Vec<usize>>,
    /// `close[p]` lists each primer of a different key whose IUPAC edit
    /// distance to `p` is below the lead (`Sheet::close_pairs_by`).
    close: Vec<Vec<usize>>,
    /// The cost a role primer that did not score counts at an end rescored
    /// at the terminal budgets: the lowest terminal budget of the sheet's
    /// primers plus one, the least any primer that did not score there can
    /// cost.
    penalty: usize,
    /// The same cost at a boundary locus (`Locus::boundary`), from the
    /// anchored budgets rescoring applies there.
    penalty_anchored: usize,
}

impl Bounds {
    /// Builds the bounds of `sheet` classified into `keys`, with `budget`
    /// and `anchored` the terminal and anchored rescoring budgets of each
    /// sheet primer (`Scorer::budgets`, `Scorer::anchored_budgets`) and
    /// `lead` the lead the classifier applies.
    pub fn new(
        sheet: &Sheet,
        keys: &Keys,
        lead: usize,
        budget: &[usize],
        anchored: &[usize],
    ) -> Bounds {
        let n = sheet.primers.len();
        let mut of_primer: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (ti, target) in sheet.targets.iter().enumerate() {
            let key = keys.of_target[ti];
            for &p in target.fwd.iter().chain(&target.rev) {
                if !of_primer[p].contains(&key) {
                    of_primer[p].push(key);
                }
            }
        }
        let mut close: Vec<Vec<usize>> = vec![Vec::new(); n];
        let pairs = sheet.close_pairs_by(lead, |i, j| keys.of_target[i] == keys.of_target[j]);
        for (a, b, _) in pairs {
            close[a].push(b);
            close[b].push(a);
        }
        let penalty = |budgets: &[usize]| budgets.iter().min().map_or(0, |&k| k + 1);
        Bounds {
            keys: of_primer,
            close,
            penalty: penalty(budget),
            penalty_anchored: penalty(anchored),
        }
    }
}

/// The primer evidence at one read end: the scores the scorer returned
/// there, the close primers of the scored ones, and the penalty of a role
/// primer that did not score there.
struct Evidence<'a> {
    /// The scores returned at this end, each within its primer's budget.
    scored: &'a [Score],
    /// For each scored primer, every close primer of a different key at the
    /// scored cost plus one: the two cannot be told apart at the lead, so the
    /// hit of either stands for both, and the primer that scored keeps a
    /// one-edit preference.
    close: Vec<Score>,
    /// `Bounds::penalty`, or `Bounds::penalty_anchored` at a boundary locus.
    penalty: usize,
}

impl<'a> Evidence<'a> {
    fn new(scored: &'a [Score], bounds: &Bounds, boundary: bool) -> Evidence<'a> {
        let close = scored
            .iter()
            .flat_map(|s| {
                bounds.close[s.primer].iter().map(move |&primer| Score {
                    primer,
                    cost: s.cost + 1,
                })
            })
            .collect();
        let penalty = if boundary {
            bounds.penalty_anchored
        } else {
            bounds.penalty
        };
        Evidence {
            scored,
            close,
            penalty,
        }
    }

    fn used(&self) -> bool {
        !self.scored.is_empty()
    }
}

/// How one end fits a target's role there.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Fit {
    /// The end holds no scores and makes no demand.
    Unused,
    /// A primer of the role scored at `cost`, or a close primer of another
    /// key stands for one; `scored` tells whether the role primer itself
    /// scored.
    Matched { cost: usize, scored: bool },
    /// The end scored a primer of another key, none of the target's own,
    /// and nothing that leads a primer that did not score by the lead, so
    /// the role primer may sit there beyond its budget: it counts at the
    /// end's penalty (`Evidence::penalty`), which every target shares.
    Penalised(usize),
    /// The end rules the target out on this strand.
    Vetoed,
}

impl Fit {
    fn cost(self) -> usize {
        match self {
            Fit::Unused | Fit::Vetoed => 0,
            Fit::Matched { cost, .. } | Fit::Penalised(cost) => cost,
        }
    }

    fn matched(self) -> bool {
        matches!(self, Fit::Matched { .. })
    }
}

/// How `end` fits `role`, one primer list of a target of key `key`, whose
/// other list is `other`. The end matches at the lowest cost among the
/// primers of `role` that scored there or that a close primer stands for,
/// so the variants of one list never compete with each other. Otherwise the
/// end is penalised when `role` is not empty, the end scored a primer of
/// another key and no primer of `other`, and every score there is within
/// `lead` of the end's penalty (`cost + lead > penalty`): no scored primer
/// leads a role primer that did not score by the lead, so the end may hold
/// the role primer beside a near-identical primer of another key. Every
/// other case vetoes the target on this strand: a scored primer of `other`
/// places the target's own primer at the wrong end, a primer that leads by
/// the lead names another target there, and scores of the target's key
/// alone leave the end to the targets that list them.
fn fit(
    end: &Evidence,
    role: &[usize],
    other: &[usize],
    key: usize,
    bounds: &Bounds,
    lead: usize,
) -> Fit {
    if !end.used() {
        return Fit::Unused;
    }
    let lowest = |scores: &[Score]| {
        scores
            .iter()
            .filter(|s| role.contains(&s.primer))
            .map(|s| s.cost)
            .min()
    };
    let scored = lowest(end.scored);
    if let Some(cost) = scored.into_iter().chain(lowest(&end.close)).min() {
        return Fit::Matched {
            cost,
            scored: scored.is_some(),
        };
    }
    let own = end.scored.iter().any(|s| other.contains(&s.primer));
    let foreign = end
        .scored
        .iter()
        .any(|s| bounds.keys[s.primer].iter().any(|&k| k != key));
    let within_lead = end.scored.iter().all(|s| s.cost + lead > end.penalty);
    if !role.is_empty() && !own && foreign && within_lead {
        Fit::Penalised(end.penalty)
    } else {
        Fit::Vetoed
    }
}

/// Whether any primer of `role` was scored at `scores`, for the orientation
/// check, which asks about raw primer placement rather than the
/// strand-consistent fit `fit` performs.
fn scored(scores: &[Score], role: &[usize]) -> bool {
    scores.iter().any(|s| role.contains(&s.primer))
}

/// One target/strand that `classify` keeps for its key.
#[derive(Clone, Copy, Debug)]
struct Candidate {
    /// Whether an end is penalised rather than matched.
    penalised: bool,
    /// The cost summed over the used ends.
    cost: usize,
    /// Whether a role primer scored at a matched end, rather than only a
    /// close primer standing for one.
    scored: bool,
    /// Index into `Sheet::targets`.
    target: usize,
    /// The strand the target was found on.
    strand: Strand,
}

impl Candidate {
    /// The order of candidates within one key: matched at every used end
    /// first, then cheaper, then with a scored role primer. Equal ranks keep
    /// the first candidate tried unless a length window breaks the tie.
    fn rank(&self) -> (bool, usize, bool) {
        (self.penalised, self.cost, !self.scored)
    }
}

/// Whether a `(strand, require)` combination is satisfied by which ends
/// matched a primer of the target. `require` names a primer role
/// (`Fwd`/`Rev`) rather than a physical end, so which end it demands depends
/// on the strand.
fn passes_require(require: Require, strand: Strand, five: bool, three: bool) -> bool {
    match require {
        Require::Both => five && three,
        Require::Either => five || three,
        Require::Fwd => match strand {
            Strand::Plus => five,
            Strand::Minus => three,
        },
        Require::Rev => match strand {
            Strand::Plus => three,
            Strand::Minus => five,
        },
    }
}

/// Classifies one segment from the primers rescored at its ends, both
/// rescored at the terminal budgets; see `classify_at`.
pub fn classify(
    sheet: &Sheet,
    keys: &Keys,
    bounds: &Bounds,
    rules: Rules,
    five: &[Score],
    three: &[Score],
) -> Call {
    classify_at(sheet, keys, bounds, rules, five, three, [false, false])
}

/// Classifies one segment from the primers rescored at its ends. An empty
/// `five` or `three` slice means that end has no primer evidence.
/// `boundary` tells, for the 5' and the 3' end, whether the end was
/// rescored at a boundary locus (`Locus::boundary`), which sets the penalty
/// of a role primer that did not score there.
///
/// For each target and strand, each end is fitted to the list that strand's
/// role puts there (`Plus`: forward at 5', reverse at 3'; `Minus`: reverse at
/// 5', forward at 3'), as `fit` describes: an end without scores makes no
/// demand, a matched end counts at its cost, a penalised end at the end's
/// penalty, the same for every target, and a vetoed end rules the target out
/// on that strand. A target/strand with no vetoed end and at least one
/// matched end is a candidate at the cost summed over its used ends, kept
/// only when its matched ends satisfy `rules.require`. Each key keeps its
/// best candidate by `Candidate::rank`. Keys with a penalised candidate
/// compete only when no key has a candidate matched at every used end. The
/// cheapest competing key is assigned when its candidate has a role primer
/// that scored, not only one a close primer stands for, and it beats every
/// other competing key by at least `rules.lead`; two keys tying for cheapest
/// are always ambiguous, whatever the lead. With no kept candidate at all: a
/// candidate whose used ends all matched, one of them by a scored role
/// primer, that failed `require` gives `Unassigned::Require`; otherwise a
/// target with primers of one role list scored at both ends gives
/// `Unassigned::Orientation`; otherwise the call is `Ambiguous`.
///
/// A target that lists one primer in both role lists is consistent on both
/// strands at the same cost whichever ends score that primer, and takes
/// `Strand::Plus`, the first strand tried.
pub fn classify_at(
    sheet: &Sheet,
    keys: &Keys,
    bounds: &Bounds,
    rules: Rules,
    five: &[Score],
    three: &[Score],
    boundary: [bool; 2],
) -> Call {
    classify_at_length(sheet, keys, bounds, rules, five, three, boundary, None)
}

/// Classifies primer evidence with a final segment length. Equal-rank
/// targets within one key prefer a compatible length window. The selected
/// target's length window is checked after competition between keys.
#[allow(clippy::too_many_arguments)]
pub(crate) fn classify_at_length(
    sheet: &Sheet,
    keys: &Keys,
    bounds: &Bounds,
    rules: Rules,
    five: &[Score],
    three: &[Score],
    boundary: [bool; 2],
    len: Option<usize>,
) -> Call {
    let ends = match (!five.is_empty(), !three.is_empty()) {
        (true, true) => Ends::Both,
        (true, false) => Ends::Five,
        (false, true) => Ends::Three,
        (false, false) => return Call::Unassigned(Unassigned::NoPrimer),
    };
    let five = Evidence::new(five, bounds, boundary[0]);
    let three = Evidence::new(three, bounds, boundary[1]);

    let length_fits = |target: usize| {
        len.is_none_or(|len| {
            sheet.targets[target]
                .len
                .is_none_or(|(min, max)| (min..=max).contains(&len))
        })
    };
    let mut structurally_consistent = false;
    let mut best: Vec<Option<Candidate>> = vec![None; keys.names.len()];

    for (ti, target) in sheet.targets.iter().enumerate() {
        let key = keys.of_target[ti];
        for strand in [Strand::Plus, Strand::Minus] {
            let (five_role, three_role) = match strand {
                Strand::Plus => (&target.fwd, &target.rev),
                Strand::Minus => (&target.rev, &target.fwd),
            };
            let five_fit = fit(&five, five_role, three_role, key, bounds, rules.lead);
            let three_fit = fit(&three, three_role, five_role, key, bounds, rules.lead);
            let fits = [five_fit, three_fit];
            if fits.contains(&Fit::Vetoed) || !fits.iter().any(|f| f.matched()) {
                continue;
            }
            let scored = fits
                .iter()
                .any(|f| matches!(f, Fit::Matched { scored: true, .. }));
            let penalised = fits.iter().any(|f| matches!(f, Fit::Penalised(_)));
            structurally_consistent |= scored && !penalised;
            if !passes_require(
                rules.require,
                strand,
                five_fit.matched(),
                three_fit.matched(),
            ) {
                continue;
            }
            let candidate = Candidate {
                penalised,
                cost: five_fit.cost() + three_fit.cost(),
                scored,
                target: ti,
                strand,
            };
            let slot = &mut best[key];
            if slot.is_none_or(|current| {
                (candidate.rank(), !length_fits(candidate.target))
                    < (current.rank(), !length_fits(current.target))
            }) {
                *slot = Some(candidate);
            }
        }
    }

    let matched_everywhere = best.iter().flatten().any(|c| !c.penalised);
    let mut entries: Vec<(usize, Candidate)> = best
        .into_iter()
        .enumerate()
        .filter_map(|(key, slot)| slot.map(|c| (key, c)))
        .filter(|(_, c)| !(matched_everywhere && c.penalised))
        .collect();

    if entries.is_empty() {
        if structurally_consistent {
            return Call::Unassigned(Unassigned::Require);
        }
        for target in &sheet.targets {
            let fwd_both = scored(five.scored, &target.fwd) && scored(three.scored, &target.fwd);
            let rev_both = scored(five.scored, &target.rev) && scored(three.scored, &target.rev);
            if fwd_both || rev_both {
                return Call::Unassigned(Unassigned::Orientation);
            }
        }
        return Call::Ambiguous;
    }

    entries.sort_by_key(|(_, c)| c.cost);
    let (key, winner) = entries[0];
    if entries
        .iter()
        .filter(|(_, c)| c.cost == winner.cost)
        .count()
        > 1
    {
        return Call::Ambiguous;
    }
    let leads = match entries.get(1) {
        Some((_, second)) => winner.cost + rules.lead <= second.cost,
        None => true,
    };
    if !winner.scored || !leads {
        return Call::Ambiguous;
    }
    let call = Call::Assigned {
        key,
        target: winner.target,
        strand: winner.strand,
        ends,
    };
    len.map_or(call, |len| check_length(sheet, call, len))
}

/// Narrows an `Assigned` call to `Unassigned(Length)` when `len` falls
/// outside its target's `len` window. Passes every other call through
/// unchanged, including an `Assigned` call whose target has no window.
pub fn check_length(sheet: &Sheet, call: Call, len: usize) -> Call {
    if let Call::Assigned { target, .. } = call
        && let Some((min, max)) = sheet.targets[target].len
        && !(min..=max).contains(&len)
    {
        return Call::Unassigned(Unassigned::Length);
    }
    call
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Targets `A` (fwd `fA`, rev `rA`), `B` (fwd `fB`, rev `rB`) and `C`
    /// (fwd `fA`, rev `rC`): `C` shares its forward primer with `A`. Sheet
    /// order gives primer indices `fA=0, rA=1, fB=2, rB=3, rC=4` and target
    /// indices `A=0, B=1, C=2`. `A` and `C` share group `G`; `B` is its own
    /// group.
    fn fixture() -> Sheet {
        Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\n\
             A\tAAAAAAAAAAA\tCCCCCCCCCCC\tG\n\
             B\tGGGGGGGGGGG\tTTTTTTTTTTT\tB\n\
             C\tAAAAAAAAAAA\tACACACACACA\tG\n",
        )
        .unwrap()
    }

    fn rules(require: Require, lead: usize) -> Rules {
        Rules { require, lead }
    }

    fn score(primer: usize, cost: usize) -> Score {
        Score { primer, cost }
    }

    /// The bounds of `sheet` at lead 2, with every primer at budget 3.
    fn bounds(sheet: &Sheet, keys: &Keys) -> Bounds {
        let budgets = vec![3; sheet.primers.len()];
        Bounds::new(sheet, keys, 2, &budgets, &budgets)
    }

    #[test]
    fn plus_both_ends() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(0, 0)],
            &[score(1, 1)],
        );
        assert_eq!(
            call,
            Call::Assigned {
                key: 0,
                target: 0,
                strand: Strand::Plus,
                ends: Ends::Both,
            }
        );
    }

    #[test]
    fn minus_both_ends() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(1, 0)],
            &[score(0, 0)],
        );
        assert_eq!(
            call,
            Call::Assigned {
                key: 0,
                target: 0,
                strand: Strand::Minus,
                ends: Ends::Both,
            }
        );
    }

    #[test]
    fn either_five_only() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(2, 1)],
            &[],
        );
        assert_eq!(
            call,
            Call::Assigned {
                key: 1,
                target: 1,
                strand: Strand::Plus,
                ends: Ends::Five,
            }
        );
    }

    #[test]
    fn both_rule_rejects_one_end() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Both, 2),
            &[score(2, 0)],
            &[],
        );
        assert_eq!(call, Call::Unassigned(Unassigned::Require));
    }

    #[test]
    fn fwd_rule_needs_fwd_primer() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        // five=[rB:0] alone is only consistent as B on Minus (reverse primer
        // at 5'), whose forward primer would sit at the unused 3' end.
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Fwd, 2),
            &[score(3, 0)],
            &[],
        );
        assert_eq!(call, Call::Unassigned(Unassigned::Require));
    }

    #[test]
    fn f_at_both_ends_is_orientation() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(2, 0)],
            &[score(2, 0)],
        );
        assert_eq!(call, Call::Unassigned(Unassigned::Orientation));
    }

    #[test]
    fn conflicting_ends_are_ambiguous() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        // five scores A's forward primer, three scores B's reverse primer,
        // both exactly. Each exact primer leads a primer that did not score
        // by the lead, so it rules out A and C at 3' and B at 5'.
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(0, 0)],
            &[score(3, 0)],
        );
        assert_eq!(call, Call::Ambiguous);
    }

    #[test]
    fn shared_primer_one_end_is_ambiguous() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        // A and C both have fA as their forward primer; with only the 5' end
        // scored, both tie at cost 0 and neither target's identity is
        // resolved.
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(0, 0)],
            &[],
        );
        assert_eq!(call, Call::Ambiguous);
    }

    #[test]
    fn shared_primer_resolved_by_other_end() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        // The 3' end scores C's reverse primer (rC) exactly, which A does
        // not have, so A is ruled out there and C, matched at both ends, is
        // assigned.
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(0, 0)],
            &[score(4, 0)],
        );
        assert_eq!(
            call,
            Call::Assigned {
                key: 2,
                target: 2,
                strand: Strand::Plus,
                ends: Ends::Both,
            }
        );
    }

    #[test]
    fn lead_boundary() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        // rA (Minus, target A) and fB (Plus, target B) do not share a
        // primer with any other target, so only A and B compete.
        let assigned = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(1, 1), score(2, 3)],
            &[],
        );
        assert_eq!(
            assigned,
            Call::Assigned {
                key: 0,
                target: 0,
                strand: Strand::Minus,
                ends: Ends::Five,
            },
            "cost 1 leads cost 3 by 2: assigned"
        );

        let ambiguous = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(1, 2), score(2, 3)],
            &[],
        );
        assert_eq!(
            ambiguous,
            Call::Ambiguous,
            "cost 2 leads cost 3 by only 1, short of lead 2"
        );
    }

    #[test]
    fn group_level_collapses_same_group() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Group);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(0, 0)],
            &[],
        );
        // A and C are both consistent at cost 0 in the same group key; they
        // do not compete, and A wins the within-key tie by lower target
        // index.
        assert_eq!(
            call,
            Call::Assigned {
                key: keys.of_target[0],
                target: 0,
                strand: Strand::Plus,
                ends: Ends::Five,
            }
        );
    }

    #[test]
    fn no_scores_is_no_primer() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[],
            &[],
        );
        assert_eq!(call, Call::Unassigned(Unassigned::NoPrimer));
    }

    /// Targets `M` (forward `f1`, `f2`; reverse `r1`, `r2`), `N` (forward
    /// `f3`; reverse `r3`) and `O` (forward `f4`, no reverse). Sheet order
    /// gives primer indices `f1=0, f2=1, r1=2, r2=3, f3=4, r3=5, f4=6` and
    /// target indices `M=0, N=1, O=2`.
    fn mix_fixture() -> Sheet {
        Sheet::parse_tsv(
            "target\tfwd\trev\n\
             M\tAAAAAAAAAAA,AAAAAAAAAAC\tCCCCCCCCCCC,CCCCCCCCCCA\n\
             N\tGGGGGGGGGGG\tTTTTTTTTTTT\n\
             O\tACACACACACA\t\n",
        )
        .unwrap()
    }

    fn assigned(target: usize, strand: Strand, ends: Ends) -> Call {
        Call::Assigned {
            key: target,
            target,
            strand,
            ends,
        }
    }

    #[test]
    fn second_forward_variant_alone_is_assigned() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(1, 1)],
                &[]
            ),
            assigned(0, Strand::Plus, Ends::Five)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[],
                &[score(1, 0)]
            ),
            assigned(0, Strand::Minus, Ends::Three)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(3, 0)],
                &[score(1, 1)]
            ),
            assigned(0, Strand::Minus, Ends::Both),
            "any reverse variant pairs with any forward variant"
        );
    }

    #[test]
    fn equal_cost_variants_of_one_target_are_assigned() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(0, 1), score(1, 1)],
            &[score(2, 0), score(3, 0)],
        );
        assert_eq!(call, assigned(0, Strand::Plus, Ends::Both));
    }

    /// An end counts at the lowest cost among the scored variants of the
    /// role, so the lead is measured from the best variant.
    #[test]
    fn lowest_cost_variant_sets_the_end_cost() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(0, 3), score(1, 0), score(4, 2)],
                &[]
            ),
            assigned(0, Strand::Plus, Ends::Five),
            "variant f2 at cost 0 leads N at cost 2"
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(0, 3), score(1, 1), score(4, 2)],
                &[]
            ),
            Call::Ambiguous,
            "variant f2 at cost 1 leads N at cost 2 by 1, short of lead 2"
        );
    }

    #[test]
    fn variants_of_different_keys_within_the_lead_are_ambiguous() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(1, 0), score(4, 1)],
                &[]
            ),
            Call::Ambiguous
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(1, 0), score(4, 0)],
                &[]
            ),
            Call::Ambiguous
        );
    }

    #[test]
    fn both_rule_takes_any_variant_at_each_end() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let both = rules(Require::Both, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                both,
                &[score(1, 0)],
                &[score(3, 1)]
            ),
            assigned(0, Strand::Plus, Ends::Both)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                both,
                &[score(1, 0)],
                &[]
            ),
            Call::Unassigned(Unassigned::Require)
        );
    }

    #[test]
    fn role_rules_take_any_variant_of_the_role_list() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                rules(Require::Fwd, 2),
                &[],
                &[score(1, 0)]
            ),
            assigned(0, Strand::Minus, Ends::Three)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                rules(Require::Rev, 2),
                &[],
                &[score(1, 0)]
            ),
            Call::Unassigned(Unassigned::Require)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                rules(Require::Rev, 2),
                &[score(3, 0)],
                &[]
            ),
            assigned(0, Strand::Minus, Ends::Five)
        );
    }

    #[test]
    fn empty_reverse_list_under_either() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(6, 0)],
                &[]
            ),
            assigned(2, Strand::Plus, Ends::Five)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[],
                &[score(6, 1)]
            ),
            assigned(2, Strand::Minus, Ends::Three)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(6, 0)],
                &[score(5, 0)]
            ),
            Call::Ambiguous,
            "O has no reverse list to place at 3', and O's exact forward primer leads N's \
             forward primer that did not score by the lead, so it rules N out at 5'"
        );
    }

    /// Two variants of one role list at the two ends place the same role at
    /// both ends, so no strand is consistent.
    #[test]
    fn forward_variants_at_both_ends_are_orientation() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Either, 2),
            &[score(0, 0)],
            &[score(1, 0)],
        );
        assert_eq!(call, Call::Unassigned(Unassigned::Orientation));
    }

    /// Targets `BC01` and `BC02`, each with one sequence in both role lists.
    /// Sheet order gives primer indices `BC01=0, BC02=1`.
    fn symmetric_fixture() -> Sheet {
        Sheet::parse_tsv(
            "target\tfwd\trev\n\
             BC01\tAAAAAAAAAAA\tAAAAAAAAAAA\n\
             BC02\tGGGGGGGGGGG\tGGGGGGGGGGG\n",
        )
        .unwrap()
    }

    #[test]
    fn one_primer_at_both_ends_of_a_symmetric_target_is_assigned() {
        let sheet = symmetric_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        for require in [Require::Either, Require::Both, Require::Fwd, Require::Rev] {
            assert_eq!(
                classify(
                    &sheet,
                    &keys,
                    &bounds(&sheet, &keys),
                    rules(require, 2),
                    &[score(0, 0)],
                    &[score(0, 1)]
                ),
                assigned(0, Strand::Plus, Ends::Both),
                "{require:?}"
            );
        }
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                rules(Require::Either, 2),
                &[score(0, 0)],
                &[score(1, 0)]
            ),
            Call::Ambiguous,
            "the two ends score different targets"
        );
    }

    #[test]
    fn one_end_of_a_symmetric_target_follows_the_end_rule() {
        let sheet = symmetric_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[score(1, 0)],
                &[]
            ),
            assigned(1, Strand::Plus, Ends::Five)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                either,
                &[],
                &[score(1, 1)]
            ),
            assigned(1, Strand::Plus, Ends::Three)
        );
        let both = rules(Require::Both, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                both,
                &[score(1, 0)],
                &[]
            ),
            Call::Unassigned(Unassigned::Require)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                both,
                &[],
                &[score(1, 0)]
            ),
            Call::Unassigned(Unassigned::Require)
        );
    }

    /// Targets `V34` (forward 341F, reverse 785R) and `V4` (forward 515F,
    /// reverse 806R), whose reverse primers are one edit apart. Sheet order
    /// gives primer indices `341F=0, 785R=1, 515F=2, 806R=3` and target
    /// indices `V34=0, V4=1`. Budgets are 3 for 341F and 515F and 4 for the
    /// reverse primers.
    fn amplicon_fixture() -> (Sheet, Keys, Bounds) {
        let sheet = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             V34\tCCTACGGGNGGCWGCAG\tGACTACHVGGGTATCTAATCC\n\
             V4\tGTGYCAGCMGCCGCGGTAA\tGGACTACNVGGGTWTCTAAT\n",
        )
        .unwrap();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let bounds = Bounds::new(&sheet, &keys, 2, &[3, 4, 3, 4], &[3, 4, 3, 4]);
        (sheet, keys, bounds)
    }

    /// A V4 read whose clipped 806R did not score while 785R did at the
    /// same end: 806R counts at 785R's cost plus one, so V4 is matched at
    /// both ends and V34, penalised at the 5' end, does not compete.
    #[test]
    fn close_primer_of_another_key_does_not_veto_the_target() {
        let (sheet, keys, bounds) = amplicon_fixture();
        let either = rules(Require::Either, 2);
        for cost_515f in 0..=3 {
            assert_eq!(
                classify(
                    &sheet,
                    &keys,
                    &bounds,
                    either,
                    &[score(2, cost_515f)],
                    &[score(1, 2)]
                ),
                assigned(1, Strand::Plus, Ends::Both),
                "515F at {cost_515f}"
            );
        }
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                either,
                &[score(1, 2)],
                &[score(2, 1)]
            ),
            assigned(1, Strand::Minus, Ends::Both),
            "the same read on the minus strand"
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                either,
                &[score(2, 1)],
                &[score(1, 0), score(3, 4)]
            ),
            assigned(1, Strand::Plus, Ends::Both),
            "806R counts at the lower of its own cost and 785R's cost plus one"
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                rules(Require::Both, 2),
                &[score(2, 1)],
                &[score(1, 2)]
            ),
            assigned(1, Strand::Plus, Ends::Both),
            "a close primer satisfies the end rule for its key"
        );
    }

    /// A V4 read whose 806R is clipped while 785R scores within its budget
    /// is never assigned V34, at any cost of 785R or 515F: at lead 2 the
    /// close 785R stands in for 806R, and at lead 1, where the primers are
    /// not close, 515F at the 5' end rules V34 out unless it is within the
    /// lead of the penalty.
    #[test]
    fn clipped_806r_is_never_assigned_v34() {
        let (sheet, keys, _) = amplicon_fixture();
        for lead in [1, 2] {
            let budgets = [3, 4, 3, 4];
            let bounds = Bounds::new(&sheet, &keys, lead, &budgets, &budgets);
            for cost_515f in 0..=3 {
                for cost_785r in 0..=4 {
                    for (five, three, strand) in [
                        (score(2, cost_515f), score(1, cost_785r), Strand::Plus),
                        (score(1, cost_785r), score(2, cost_515f), Strand::Minus),
                    ] {
                        let call = classify(
                            &sheet,
                            &keys,
                            &bounds,
                            rules(Require::Either, lead),
                            &[five],
                            &[three],
                        );
                        let case = format!("lead {lead}, 515F {cost_515f}, 785R {cost_785r}");
                        assert!(
                            !matches!(call, Call::Assigned { key: 0, .. }),
                            "{case}: {call:?}"
                        );
                        if lead == 2 {
                            assert_eq!(call, assigned(1, strand, Ends::Both), "{case}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn clear_v34_read_stays_v34() {
        let (sheet, keys, bounds) = amplicon_fixture();
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                either,
                &[score(0, 0)],
                &[score(1, 1)]
            ),
            assigned(0, Strand::Plus, Ends::Both)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                either,
                &[score(0, 1)],
                &[score(1, 0), score(3, 2)]
            ),
            assigned(0, Strand::Plus, Ends::Both)
        );
    }

    /// With only the shared end scored, near-identical primers of two keys
    /// are judged by the lead, whichever of them cleared its budget.
    #[test]
    fn close_primers_at_one_end_are_judged_by_the_lead() {
        let (sheet, keys, bounds) = amplicon_fixture();
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(&sheet, &keys, &bounds, either, &[], &[score(1, 0)]),
            Call::Ambiguous
        );
        assert_eq!(
            classify(&sheet, &keys, &bounds, either, &[], &[score(3, 1)]),
            Call::Ambiguous
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                rules(Require::Either, 1),
                &[],
                &[score(1, 0)]
            ),
            assigned(0, Strand::Plus, Ends::Three),
            "at lead 1 the primers are not close, and 785R alone names V34"
        );
    }

    /// A penalised end needs a scored primer of another key: an end whose
    /// only score is the target's own primer of the other role still vetoes
    /// the target, and a target no role primer matches is no candidate.
    #[test]
    fn penalty_needs_another_keys_primer() {
        let (sheet, keys, bounds) = amplicon_fixture();
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                either,
                &[score(2, 0)],
                &[score(2, 0)]
            ),
            Call::Unassigned(Unassigned::Orientation)
        );
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                rules(Require::Both, 2),
                &[score(0, 0)],
                &[]
            ),
            Call::Unassigned(Unassigned::Require)
        );
    }

    /// Bounds of `fixture` at lead 2 with `rA` at budget 1 and every other
    /// primer at budget 3, so the penalty is 2 at every end.
    fn unequal_bounds(sheet: &Sheet, keys: &Keys) -> Bounds {
        let budgets = [3, 1, 3, 3, 3];
        Bounds::new(sheet, keys, 2, &budgets, &budgets)
    }

    /// A target matched at every used end outranks a cheaper target with a
    /// penalised end. Sheet order gives primer indices `fA=0, rA=1, fB=2,
    /// rB=3, rC=4`; `A`, penalised at 3', costs the penalty 2.
    #[test]
    fn a_target_matched_at_every_end_outranks_a_penalised_one() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let bounds = unequal_bounds(&sheet, &keys);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                rules(Require::Either, 2),
                &[score(0, 0)],
                &[score(4, 3)]
            ),
            Call::Assigned {
                key: 2,
                target: 2,
                strand: Strand::Plus,
                ends: Ends::Both,
            }
        );
    }

    /// Exact primers of two keys at the two ends, a chimera, are not
    /// assigned to either key, whatever the budgets of the primers that did
    /// not score: each exact primer leads the penalty by the lead and rules
    /// out the targets that do not list it. With the penalty taken from the
    /// budget of the missing primer, `A`, missing `rA` at budget 1, led `B`
    /// and `C`, missing primers at budget 3, and was assigned.
    #[test]
    fn exact_primers_of_two_keys_at_the_two_ends_are_not_assigned() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let bounds = unequal_bounds(&sheet, &keys);
        for require in [Require::Either, Require::Both] {
            assert_eq!(
                classify(
                    &sheet,
                    &keys,
                    &bounds,
                    rules(require, 2),
                    &[score(0, 0)],
                    &[score(3, 0)]
                ),
                Call::Ambiguous,
                "{require:?}"
            );
        }
    }

    /// Penalised targets compete at one penalty, so no target wins by its
    /// missing primer having the lower budget: `A` and `C` both match `fA`
    /// at 5' and are penalised at 3', where `rB` scored at its budget, and
    /// tie, although `rA` has the lower budget. At the group level they
    /// share a key, which is assigned.
    #[test]
    fn penalised_targets_compete_at_one_penalty() {
        let sheet = fixture();
        let bounds = |keys: &Keys| unequal_bounds(&sheet, keys);
        let targets = Keys::new(&sheet, KeyLevel::Target);
        let either = rules(Require::Either, 2);
        assert_eq!(
            classify(
                &sheet,
                &targets,
                &bounds(&targets),
                either,
                &[score(0, 0)],
                &[score(3, 3)]
            ),
            Call::Ambiguous
        );
        let groups = Keys::new(&sheet, KeyLevel::Group);
        assert_eq!(
            classify(
                &sheet,
                &groups,
                &bounds(&groups),
                either,
                &[score(0, 0)],
                &[score(3, 3)]
            ),
            Call::Assigned {
                key: groups.of_target[0],
                target: 0,
                strand: Strand::Plus,
                ends: Ends::Both,
            }
        );
        assert_eq!(
            classify(
                &sheet,
                &targets,
                &bounds(&targets),
                rules(Require::Both, 2),
                &[score(0, 0)],
                &[score(3, 3)]
            ),
            Call::Ambiguous,
            "a penalised end does not satisfy the both rule"
        );
    }

    /// One exact primer at both ends, beside a weak primer of another key at
    /// the 3' end, is an orientation conflict: the exact `fB` rules out every
    /// target that does not list it at either end, and `B` cannot hold its
    /// forward primer at both ends.
    #[test]
    fn exact_primer_at_both_ends_beside_a_weak_one_is_orientation() {
        let sheet = fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds(&sheet, &keys),
                rules(Require::Either, 2),
                &[score(2, 0)],
                &[score(2, 0), score(1, 3)]
            ),
            Call::Unassigned(Unassigned::Orientation)
        );
    }

    /// At a boundary locus the penalty follows the anchored budgets that
    /// rescoring applied there: with terminal budgets 3 and anchored budgets
    /// 4, `N`'s reverse primer at cost 3 is within the lead of the terminal
    /// penalty 4, so `M`, matched at 5', is penalised at 3' and assigned, and
    /// not within the lead of the anchored penalty 5, so it rules `M` out.
    #[test]
    fn boundary_locus_takes_the_anchored_penalty() {
        let sheet = mix_fixture();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let n = sheet.primers.len();
        let bounds = Bounds::new(&sheet, &keys, 2, &vec![3; n], &vec![4; n]);
        let either = rules(Require::Either, 2);
        let call = |boundary| {
            classify_at(
                &sheet,
                &keys,
                &bounds,
                either,
                &[score(0, 0)],
                &[score(5, 3)],
                boundary,
            )
        };
        assert_eq!(call([false, false]), assigned(0, Strand::Plus, Ends::Both));
        assert_eq!(call([true, false]), assigned(0, Strand::Plus, Ends::Both));
        assert_eq!(call([false, true]), Call::Ambiguous);
    }

    /// A key whose every matched end is only stood for by close primers of
    /// other keys is not assigned, however cheap.
    #[test]
    fn a_key_matched_only_by_close_primers_is_not_assigned() {
        let sheet = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             T1\tAGRGTTYGATYMTGGCTCAG\tTCCTCCGCTTATTGATATGC\n\
             T2\tGTACACACCGCCCGTC\tCGGTTACCTTGTTACGACTT\n\
             T3\tAGRGTTYGATYMTGGCTCAC\tCGGTTACCTTGTTACGACTA\n",
        )
        .unwrap();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let budgets = vec![3; sheet.primers.len()];
        let bounds = Bounds::new(&sheet, &keys, 2, &budgets, &budgets);
        let (f1, r2) = (sheet.targets[0].fwd[0], sheet.targets[1].rev[0]);
        // The exact primers rule out T1 at 3' and T2 at 5'; T3 costs 1 + 1
        // from close primers only.
        assert_eq!(
            classify(
                &sheet,
                &keys,
                &bounds,
                rules(Require::Either, 2),
                &[score(f1, 0)],
                &[score(r2, 0)]
            ),
            Call::Ambiguous
        );
    }

    #[test]
    fn length_window() {
        let sheet = Sheet::parse_tsv(
            "target\tfwd\trev\tmin_len\tmax_len\n\
             A\tAAAAAAAAAAA\tCCCCCCCCCCC\t300\t1000\n",
        )
        .unwrap();
        let assigned = Call::Assigned {
            key: 0,
            target: 0,
            strand: Strand::Plus,
            ends: Ends::Both,
        };
        assert_eq!(
            check_length(&sheet, assigned, 200),
            Call::Unassigned(Unassigned::Length)
        );
        assert_eq!(check_length(&sheet, assigned, 500), assigned);
    }

    #[test]
    fn length_windows_do_not_override_stronger_primer_evidence() {
        let mut sheet = fixture();
        sheet.targets[0].len = Some((100, 200));
        sheet.targets[2].len = Some((300, 500));
        let keys = Keys::new(&sheet, KeyLevel::Group);
        let call = classify_at_length(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Both, 2),
            &[score(0, 0)],
            &[score(1, 0), score(4, 1)],
            [false; 2],
            Some(400),
        );
        assert_eq!(call, Call::Unassigned(Unassigned::Length));
    }

    #[test]
    fn length_windows_do_not_resolve_competing_keys() {
        let sheet = Sheet::parse_tsv(
            "target\tfwd\trev\tmin_len\tmax_len\n\
             A\tAAAAAAAAAAA\tCCCCCCCCCCC\t100\t200\n\
             B\tAAAAAAAAAAA\tCCCCCCCCCCC\t300\t500\n",
        )
        .unwrap();
        let keys = Keys::new(&sheet, KeyLevel::Target);
        let call = classify_at_length(
            &sheet,
            &keys,
            &bounds(&sheet, &keys),
            rules(Require::Both, 2),
            &[score(0, 0)],
            &[score(1, 0)],
            [false; 2],
            Some(400),
        );
        assert_eq!(call, Call::Ambiguous);
    }
}
