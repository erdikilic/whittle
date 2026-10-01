//! Rescoring of sheet primers at a located primer locus.
//!
//! `Scorer` holds every sheet primer's forward and reverse-complement
//! patterns and terminal edit budget. `Scorer::score` widens a located locus
//! into a search window, aligns every primer there in the orientation valid
//! for the end, and returns the minimum-cost hit per primer within its
//! budget. Primers that share a length and budgets are aligned together,
//! one per SIMD lane. `Scorer::score_at_site` first drops the primers that
//! the whole hit behind a locus rules out (`site_bound`).

use std::borrow::Cow;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::OnceLock;

use super::classify::Score;
use super::sheet::Sheet;
use crate::adapter::search::{self, AmbiguousSearcher};
use crate::adapter::{
    Adapter, MIN_OVERLAP, PrimerSite, anchored_budgets, normalize_into, residue_within_budget,
    reverse_complement, terminal_budget,
};

/// Which end of a read a locus was located at, and so the orientation every
/// sheet primer is aligned in there: as given at the 5' end, reverse
/// complemented at the 3' end.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum End {
    /// The 5' end: primers are aligned as stored in the sheet.
    Five,
    /// The 3' end: primers are aligned reverse complemented.
    Three,
}

/// One sheet primer's search patterns and edit budgets.
#[derive(Debug, Clone)]
struct Entry {
    /// Pattern searched at `End::Five`, as stored in the sheet.
    forward: Vec<u8>,
    /// Pattern searched at `End::Three`: the reverse complement of `forward`.
    reverse: Vec<u8>,
    /// Minimum cost at or below which a hit of this primer scores.
    budget: usize,
    /// Minimum cost at or below which a whole hit of this primer scores at
    /// a boundary locus (`Locus::boundary`): its anchored budget over the
    /// sheet's primers together (`anchored_budgets`), the budget that located
    /// the locus.
    anchored: usize,
}

/// Sheet primers that share a length and both budgets, so one search covers
/// them: the window and the budget are those of each member's own search.
#[derive(Debug, Clone)]
struct Group {
    /// Indices into `Scorer::entries`, ascending.
    members: Vec<usize>,
    /// Length of every member, in bases.
    len: usize,
    /// Budget of every member.
    budget: usize,
    /// Anchored budget of every member.
    anchored: usize,
}

/// Rescores sheet primers at a located primer locus.
///
/// Holds every sheet primer's search patterns and budget and holds no
/// searcher, so it is `Send + Sync` and shared by every render worker
/// through an `Arc`. Each thread keeps its own searcher in thread-local
/// state, built lazily by `score`.
#[derive(Debug, Clone)]
pub struct Scorer {
    entries: Vec<Entry>,
    /// The entries grouped by length and budgets; every entry is in one
    /// group.
    groups: Vec<Group>,
    /// Longest sheet primer, in bases; sets how far a locus widens into a
    /// search window.
    max_primer_len: usize,
    /// The run's adapter error rate, at `f64` precision: the partial budget
    /// of an overhanging hit's aligned portion is computed at this
    /// precision, as at the engine's own terminal search.
    error_rate: f64,
    /// Overhang cost per hanging base: `error_rate` at the searcher's `f32`
    /// precision.
    alpha: f32,
    /// The uppercase sequence of every configured adapter entry, indexed
    /// like `AdapterConfig::adapters`, which a locus site names; empty
    /// unless set by `with_sites`.
    sites: Vec<Vec<u8>>,
    /// Every sheet primer's `site_bound` per tabled site shape
    /// (`SiteShape::slot`), each computed on first use and then read
    /// without a lock by every thread; reset by `with_sites`.
    bounds: Vec<OnceLock<Box<[usize]>>>,
}

/// The shape of a rescoring window around a locus site (`site_bound`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SiteShape {
    /// Index into `Scorer::sites`.
    entry: usize,
    /// Whether the read holds the reverse complement of the entry.
    rc: bool,
    /// The end whose orientation the primers are aligned in.
    end: End,
    /// Wildcards left of the site, at most `Scorer::max_primer_len`.
    left: usize,
    /// Wildcards right of the site, at most `Scorer::max_primer_len`.
    right: usize,
}

impl SiteShape {
    /// Shapes tabled per entry, orientation and end, for flanks of at most
    /// `max` wildcards: both flanks at `max`, one flank at `max` and the
    /// other below it, or both flanks equal and below `max`.
    fn slots_per_key(max: usize) -> usize {
        3 * max + 1
    }

    /// Returns the index of this shape in a table of `slots_per_key(max)`
    /// slots per entry, orientation and end, or `None` for a shape outside
    /// it. A window side that reaches a read end takes `max` wildcards and
    /// every other side the widening of the locus, the same on both sides,
    /// so every shape of a rescoring window is tabled.
    fn slot(&self, max: usize) -> Option<usize> {
        let within = match (self.left, self.right) {
            (l, r) if l == max && r == max => 0,
            (l, r) if l == max => 1 + r,
            (l, r) if r == max => 1 + max + l,
            (l, r) if l == r => 1 + 2 * max + l,
            _ => return None,
        };
        let key = (self.entry * 2 + usize::from(self.rc)) * 2 + usize::from(self.end == End::Five);
        Some(key * Self::slots_per_key(max) + within)
    }
}

/// Returns a 4-bit mask of the bases an IUPAC code stands for, every base
/// for any other byte, so an unknown code never raises `site_bound`.
fn base_mask(code: u8) -> u8 {
    search::iupac_bases(code).map_or(0b1111, |bases| {
        bases.iter().fold(0, |mask, &base| {
            mask | match base {
                b'A' => 1,
                b'C' => 2,
                b'G' => 4,
                _ => 8,
            }
        })
    })
}

/// Returns the semi-global edit distance of `pattern` against the text of
/// `left` wildcards, `site` and `right` wildcards: the least cost of
/// aligning the whole pattern with any substring of that text, where two
/// IUPAC codes match when they share a base and a wildcard matches every
/// code.
///
/// The distance bounds from below the cost of a scoring primer hit in a
/// rescoring window around a whole hit of `site`. Read bases are normalized
/// to A, C, G, T or a symbol that matches nothing (`normalize_into`), so a
/// primer code and a site code that both match one read base share that
/// base. Composing the alignment of a primer hit in the window with the
/// alignment of the site over the locus, the window bases outside the locus
/// taking wildcards, thus aligns the primer with this text at no more than
/// the two costs together: at most `e + c` for a hit whose aligned part
/// costs `e` and a site of cost `c`. A scoring hit overhangs only a side of
/// the window that counts as a read end, which then takes further wildcards
/// for the hanging bases. An optimal alignment uses no more wildcards on a
/// side than the pattern has bases, as each further one is a deletion, so a
/// flank of `pattern.len()` wildcards gives the bound of any longer flank.
/// A hit costs at least its aligned edits, and a scoring hit at most the
/// budget `k` of its search, so a primer whose bound exceeds `c + k` has no
/// scoring hit in the window.
fn site_bound(pattern: &[u8], site: &[u8], left: usize, right: usize) -> usize {
    let text: Vec<u8> = std::iter::repeat_n(0b1111, left)
        .chain(site.iter().map(|&code| base_mask(code)))
        .chain(std::iter::repeat_n(0b1111, right))
        .collect();
    // `row[j]` is the least cost of aligning the pattern prefix so far with
    // a text substring ending at `j`; a substring may start anywhere.
    let mut row = vec![0usize; text.len() + 1];
    let mut next = vec![0usize; text.len() + 1];
    for (i, &code) in pattern.iter().enumerate() {
        let mask = base_mask(code);
        next[0] = i + 1;
        for j in 1..=text.len() {
            let substitution = usize::from(mask & text[j - 1] == 0);
            next[j] = (row[j - 1] + substitution)
                .min(row[j] + 1)
                .min(next[j - 1] + 1);
        }
        std::mem::swap(&mut row, &mut next);
    }
    row.into_iter().min().unwrap_or(0)
}

/// This thread's rescoring state: the overhang-aware, forward-only searcher,
/// rebuilt when its alpha differs from the scorer's, and the normalized read
/// buffer, the per-primer minimum costs and the group members left to
/// search, reused across calls.
struct ThreadState {
    searcher: Option<(f32, AmbiguousSearcher)>,
    normalized: Vec<u8>,
    best: Vec<Option<usize>>,
    kept: Vec<usize>,
}

thread_local! {
    static STATE: RefCell<ThreadState> = const {
        RefCell::new(ThreadState {
            searcher: None,
            normalized: Vec::new(),
            best: Vec::new(),
            kept: Vec::new(),
        })
    };
}

impl Scorer {
    /// Builds a scorer for every primer in `sheet`. `error_rate` sets each
    /// primer's terminal edit budget (`terminal_budget`) and the per-base
    /// overhang cost; `end_size` shapes the chance-match bound behind the
    /// budget, as at the engine's own terminal search. `error_rate` also
    /// sets each primer's anchored budget (`anchored_budgets`).
    pub fn new(sheet: &Sheet, error_rate: f64, end_size: usize) -> Scorer {
        let seqs: Vec<&[u8]> = sheet.primers.iter().map(|p| p.seq.as_slice()).collect();
        let entries: Vec<Entry> = sheet
            .primers
            .iter()
            .zip(anchored_budgets(&seqs, error_rate))
            .map(|(primer, anchored)| Entry {
                forward: primer.seq.clone(),
                reverse: reverse_complement(&primer.seq),
                budget: terminal_budget(&primer.seq, error_rate, end_size),
                anchored,
            })
            .collect();
        let max_primer_len = entries.iter().map(|e| e.forward.len()).max().unwrap_or(0);
        let mut by_window: BTreeMap<(usize, usize, usize), Vec<usize>> = BTreeMap::new();
        for (primer, entry) in entries.iter().enumerate() {
            by_window
                .entry((entry.forward.len(), entry.budget, entry.anchored))
                .or_default()
                .push(primer);
        }
        let groups = by_window
            .into_iter()
            .map(|((len, budget, anchored), members)| Group {
                members,
                len,
                budget,
                anchored,
            })
            .collect();
        Scorer {
            entries,
            groups,
            max_primer_len,
            error_rate,
            alpha: error_rate as f32,
            sites: Vec::new(),
            bounds: Vec::new(),
        }
    }

    /// Records the sequences of `adapters`, the configured entries that a
    /// locus site (`PrimerSite::entry`) indexes, for `score_at_site`, and
    /// empties the table of their bounds.
    pub fn with_sites(mut self, adapters: &[Adapter]) -> Scorer {
        self.sites = adapters
            .iter()
            .map(|a| a.seq.to_ascii_uppercase())
            .collect();
        let slots = self.sites.len() * 4 * SiteShape::slots_per_key(self.max_primer_len);
        self.bounds = std::iter::repeat_with(OnceLock::new).take(slots).collect();
        self
    }

    /// Returns every sheet primer's `site_bound` in the orientation `end`
    /// gives, around the whole hit of entry `entry` (its reverse complement
    /// when `rc`) with `left` and `right` wildcards beside it, each capped
    /// at the longest primer, or `None` for an entry without a recorded
    /// sequence. A tabled shape is computed once and borrowed from the
    /// table after; any other shape is computed on each call.
    fn site_bounds(
        &self,
        (entry, rc): (usize, bool),
        end: End,
        left: usize,
        right: usize,
    ) -> Option<Cow<'_, [usize]>> {
        let seq = self.sites.get(entry)?;
        let shape = SiteShape {
            entry,
            rc,
            end,
            left: left.min(self.max_primer_len),
            right: right.min(self.max_primer_len),
        };
        let compute = || -> Box<[usize]> {
            let site = if rc {
                reverse_complement(seq)
            } else {
                seq.clone()
            };
            self.entries
                .iter()
                .map(|e| {
                    let pattern = match end {
                        End::Five => &e.forward,
                        End::Three => &e.reverse,
                    };
                    site_bound(pattern, &site, shape.left, shape.right)
                })
                .collect()
        };
        let cell = shape
            .slot(self.max_primer_len)
            .and_then(|slot| self.bounds.get(slot));
        Some(match cell {
            Some(cell) => Cow::Borrowed(cell.get_or_init(compute)),
            None => Cow::Owned(compute().into_vec()),
        })
    }

    /// Returns every sheet primer's terminal edit budget, indexed like
    /// `Sheet::primers`.
    pub fn budgets(&self) -> Vec<usize> {
        self.entries.iter().map(|e| e.budget).collect()
    }

    /// Returns every sheet primer's anchored edit budget, the budget of a
    /// whole hit at a boundary locus, indexed like `Sheet::primers`.
    pub fn anchored_budgets(&self) -> Vec<usize> {
        self.entries.iter().map(|e| e.anchored).collect()
    }

    /// Returns the minimum-cost hit of every sheet primer within its budget,
    /// searched in the orientation `end` gives: at most one `Score` per
    /// primer, sorted ascending by cost, ties by primer index.
    ///
    /// The search window is `locus` widened by `max_primer_len - (locus.1 -
    /// locus.0) + 3` bases on each side, saturating, clamped to the read. A
    /// hit's overhang counts only on the side where the window reaches the
    /// matching read end; a hit overhanging the other side is discarded,
    /// since that side of the window is bounded by widening rather than by
    /// the read, and an overhang there would not reflect a genuine clip at
    /// the read end. An overhanging hit also needs at least `MIN_OVERLAP`
    /// aligned bases in the read, and the edits among those aligned
    /// bases (the hit's cost with the overhang charge discounted) must stay
    /// within the partial budget of that overlap, so a primer cannot score
    /// from an overhang alone with an unsupported aligned portion; see
    /// `residue_within_budget`.
    ///
    /// With `outer_open`, the outer edge of the locus (its start at
    /// `End::Five`, its end at `End::Three`) counts as a read end: the
    /// window stops there on that side, and a hit may overhang it under the
    /// same rules, as a primer that lost bases at a trimmed adapter junction
    /// does.
    ///
    /// With `boundary`, the locus is a whole primer at an outer layer's trim
    /// boundary (`Locus::boundary`): a hit without overhang scores within
    /// the primer's anchored budget instead of its terminal budget. A hit
    /// with overhang keeps the terminal budget, and the search runs at the
    /// larger of the two.
    ///
    /// A group of two or more primers is aligned in one search where the
    /// window and the longest overhang fit one scan block of the searcher
    /// (`search::fits_one_block`), and one primer at a time otherwise; both
    /// give every primer the hits of its own search.
    pub fn score(
        &self,
        read: &[u8],
        locus: (usize, usize),
        end: End,
        outer_open: bool,
        boundary: bool,
    ) -> Vec<Score> {
        self.score_at_site(read, locus, end, outer_open, boundary, None)
    }

    /// `score` with `site`, the whole hit spanning `locus` exactly, which
    /// skips every primer whose `site_bound` over the window exceeds the
    /// site's cost plus the primer's search budget: such a primer has no
    /// scoring hit there. The other primers are searched as `score` searches
    /// them, so the result is that of `score`. A site whose entry has no
    /// recorded sequence (`with_sites`) skips none.
    pub fn score_at_site(
        &self,
        read: &[u8],
        locus: (usize, usize),
        end: End,
        outer_open: bool,
        boundary: bool,
        site: Option<PrimerSite>,
    ) -> Vec<Score> {
        let (lo, hi) = locus;
        let widen = self.max_primer_len.saturating_sub(hi.saturating_sub(lo)) + 3;
        let start = match end {
            End::Five if outer_open => lo.min(read.len()),
            _ => lo.saturating_sub(widen),
        };
        let stop = match end {
            End::Three if outer_open => hi.min(read.len()),
            _ => hi.saturating_add(widen).min(read.len()),
        };
        let touches_start = start == 0 || (outer_open && end == End::Five);
        let touches_end = stop == read.len() || (outer_open && end == End::Three);
        if start >= stop {
            return Vec::new();
        }
        let window = &read[start..stop];
        let flank = |touches: bool, bases: usize| {
            if touches {
                bases + self.max_primer_len
            } else {
                bases
            }
        };
        let pruning = site.and_then(|site| {
            let bounds = self.site_bounds(
                (site.entry, site.rc),
                end,
                flank(touches_start, lo.saturating_sub(start)),
                flank(touches_end, stop.saturating_sub(hi)),
            )?;
            Some((bounds, site.cost))
        });

        STATE.with_borrow_mut(|state| {
            let ThreadState {
                searcher,
                normalized,
                best,
                kept,
            } = state;
            let (window, _) = normalize_into(window, normalized);
            let searcher = match searcher {
                Some((a, s)) if *a == self.alpha => s,
                slot => {
                    &mut slot
                        .insert((self.alpha, search::new_overhang_searcher_fwd(self.alpha)))
                        .1
                },
            };
            let pattern_of = |primer: usize| -> &[u8] {
                match end {
                    End::Five => &self.entries[primer].forward,
                    End::Three => &self.entries[primer].reverse,
                }
            };
            let scores = |primer: usize, len: usize, hit: &search::Hit| {
                let entry = &self.entries[primer];
                let overhang = hit.left_overhang + hit.right_overhang;
                let whole = if boundary {
                    entry.anchored
                } else {
                    entry.budget
                };
                if overhang == 0 {
                    return hit.cost <= whole;
                }
                let overlap = len - overhang;
                hit.cost <= entry.budget
                    && (hit.left_overhang == 0 || touches_start)
                    && (hit.right_overhang == 0 || touches_end)
                    && overlap >= MIN_OVERLAP
                    && residue_within_budget(
                        self.error_rate,
                        hit.cost,
                        hit.left_overhang,
                        hit.right_overhang,
                        overlap,
                    )
            };
            best.clear();
            best.resize(self.entries.len(), None);
            let mut keep = |primer: usize, len: usize, hit: search::Hit| {
                if scores(primer, len, &hit) && best[primer].is_none_or(|cost| hit.cost < cost) {
                    best[primer] = Some(hit.cost);
                }
            };
            for group in &self.groups {
                let Group {
                    members,
                    len,
                    budget,
                    anchored,
                } = group;
                let budget = if boundary {
                    *anchored.max(budget)
                } else {
                    *budget
                };
                let members: &[usize] = match &pruning {
                    Some((bounds, cost)) => {
                        kept.clear();
                        kept.extend(
                            members
                                .iter()
                                .copied()
                                .filter(|&primer| bounds[primer] <= cost + budget),
                        );
                        kept
                    },
                    None => members,
                };
                if members.len() >= 2 && search::fits_one_block(*len, window.len()) {
                    for lanes in members.chunks(PATTERN_LANES) {
                        let mut patterns: [&[u8]; PATTERN_LANES] = [&[]; PATTERN_LANES];
                        for (slot, &primer) in patterns.iter_mut().zip(lanes) {
                            *slot = pattern_of(primer);
                        }
                        search::for_each_hit_in_patterns(
                            searcher,
                            &patterns[..lanes.len()],
                            window,
                            budget,
                            |lane, hit| keep(lanes[lane], *len, hit),
                        );
                    }
                } else {
                    for &primer in members {
                        search::for_each_hit(searcher, pattern_of(primer), window, budget, |hit| {
                            keep(primer, *len, hit)
                        });
                    }
                }
            }
            let mut out: Vec<Score> = best
                .iter()
                .enumerate()
                .filter_map(|(primer, cost)| cost.map(|cost| Score { primer, cost }))
                .collect();
            out.sort_by_key(|s| (s.cost, s.primer));
            out
        })
    }
}

/// Primers aligned per search of a group: a search runs one pattern per SIMD
/// lane, and four lanes are available on every build.
const PATTERN_LANES: usize = 4;

#[cfg(test)]
mod tests {
    use super::*;

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

    /// Substitutes the base at each of `positions` for a different base.
    fn substituted(seq: &[u8], positions: &[usize]) -> Vec<u8> {
        let mut out = seq.to_vec();
        for &i in positions {
            out[i] = if out[i] == b'A' { b'C' } else { b'A' };
        }
        out
    }

    /// A one-target sheet naming `fwd` and `rev` as the target's forward and
    /// reverse primers.
    fn sheet(fwd: &[u8], rev: &[u8]) -> Sheet {
        let text = format!(
            "target\tfwd\trev\nA\t{}\t{}\n",
            String::from_utf8_lossy(fwd),
            String::from_utf8_lossy(rev),
        );
        Sheet::parse_tsv(&text).unwrap()
    }

    // Plain ACGT fixtures: a read is literal DNA, never an ambiguity code, so
    // these stand in for `fA`/`rA` everywhere except the degenerate-pattern
    // test, which needs a primer with wobble positions.
    const FA: &[u8] = b"ACGGTTCAGCATTGACCGTA";
    const RA: &[u8] = b"TTGCACGGTAACCTGATCGA";

    #[test]
    fn exact_primer_at_five_prime_scores_zero() {
        let s = sheet(FA, RA);
        let scorer = Scorer::new(&s, 0.2, 150);
        let fa_idx = s.targets[0].fwd[0];

        let mut read = FA.to_vec();
        read.extend(splitmix_dna(1, 400));

        let scores = scorer.score(&read, (0, FA.len()), End::Five, false, false);
        assert!(
            scores.contains(&Score {
                primer: fa_idx,
                cost: 0
            }),
            "{scores:?}"
        );
    }

    #[test]
    fn three_prime_uses_reverse_complement() {
        let s = sheet(FA, RA);
        let scorer = Scorer::new(&s, 0.2, 150);
        let fa_idx = s.targets[0].fwd[0];
        let ra_idx = s.targets[0].rev[0];

        let mut read = splitmix_dna(2, 400);
        let start = read.len();
        read.extend(reverse_complement(RA));

        let scores = scorer.score(&read, (start, read.len()), End::Three, false, false);
        assert!(
            scores.contains(&Score {
                primer: ra_idx,
                cost: 0
            }),
            "{scores:?}"
        );
        assert!(
            !scores.iter().any(|s| s.primer == fa_idx),
            "fA's reverse complement is not in this window: {scores:?}"
        );
    }

    #[test]
    fn degenerate_positions_cost_nothing() {
        // 27F, with wobble at R/Y/Y/M.
        let primer_27f: &[u8] = b"AGRGTTYGATYMTGGCTCAG";
        let s = sheet(primer_27f, RA);
        let scorer = Scorer::new(&s, 0.2, 150);
        let fa_idx = s.targets[0].fwd[0];

        // Every degenerate position (R, Y, Y, M) resolved to one of its
        // represented bases; an exact IUPAC-consistent read costs nothing.
        let mut read = b"AGAGTTTGATCATGGCTCAG".to_vec();
        read.extend(splitmix_dna(3, 400));

        let scores = scorer.score(&read, (0, 20), End::Five, false, false);
        assert!(
            scores.contains(&Score {
                primer: fa_idx,
                cost: 0
            }),
            "{scores:?}"
        );
    }

    /// A whole hit above the terminal budget and within the anchored budget
    /// scores at a boundary locus only.
    #[test]
    fn boundary_locus_scores_at_the_anchored_budget() {
        let primer_27f: &[u8] = b"AGRGTTYGATYMTGGCTCAG";
        let s = sheet(primer_27f, RA);
        let scorer = Scorer::new(&s, 0.2, 150);
        let fa_idx = s.targets[0].fwd[0];
        assert_eq!(terminal_budget(primer_27f, 0.2, 150), 3);
        assert_eq!(anchored_budgets(&[primer_27f, RA], 0.2)[0], 4);

        let site = substituted(b"AGAGTTTGATCATGGCTCAG", &[4, 8, 13, 17]);
        let read = [splitmix_dna(14, 60), site, splitmix_dna(15, 400)].concat();
        let locus = (60, 80);
        let scores = scorer.score(&read, locus, End::Five, false, false);
        assert!(!scores.iter().any(|s| s.primer == fa_idx), "{scores:?}");
        let scores = scorer.score(&read, locus, End::Five, false, true);
        assert!(
            scores.contains(&Score {
                primer: fa_idx,
                cost: 4
            }),
            "{scores:?}"
        );
    }

    #[test]
    fn edits_within_budget_score_beyond_are_dropped() {
        let s = sheet(FA, RA);
        let scorer = Scorer::new(&s, 0.2, 150);
        let fa_idx = s.targets[0].fwd[0];
        let budget = terminal_budget(FA, 0.2, 150);
        assert!(budget >= 2, "fixture needs a budget of at least 2 edits");

        let within = substituted(FA, &(0..budget).map(|i| 2 + i * 3).collect::<Vec<_>>());
        let mut read = within.clone();
        read.extend(splitmix_dna(4, 400));
        let scores = scorer.score(&read, (0, FA.len()), End::Five, false, false);
        assert!(
            scores.contains(&Score {
                primer: fa_idx,
                cost: budget
            }),
            "{scores:?}"
        );

        let beyond = substituted(FA, &(0..budget + 1).map(|i| 2 + i * 3).collect::<Vec<_>>());
        let mut read = beyond;
        read.extend(splitmix_dna(5, 400));
        let scores = scorer.score(&read, (0, FA.len()), End::Five, false, false);
        assert!(
            !scores.iter().any(|s| s.primer == fa_idx),
            "beyond budget {budget}: {scores:?}"
        );
    }

    #[test]
    fn clipped_primer_at_read_start_scores_via_overhang() {
        let s = sheet(FA, RA);
        let scorer = Scorer::new(&s, 0.2, 150);
        let fa_idx = s.targets[0].fwd[0];

        let clipped = &FA[5..];
        let mut read = clipped.to_vec();
        read.extend(splitmix_dna(6, 400));

        let scores = scorer.score(&read, (0, clipped.len()), End::Five, false, false);
        assert!(
            scores.iter().any(|s| s.primer == fa_idx),
            "a clip flush with the read start scores via overhang: {scores:?}"
        );
    }

    /// The partial budget of an overlap of exactly `MIN_OVERLAP` bases is 0
    /// (`partial_budget`), so an overhanging hit whose aligned portion is
    /// exactly that long tolerates no edits there, even though its total
    /// cost (aligned edits plus overhang charge) is still within the
    /// primer's whole-pattern terminal budget.
    #[test]
    fn overhang_hit_needs_the_aligned_portion_within_budget() {
        let s = sheet(FA, RA);
        let scorer = Scorer::new(&s, 0.2, 150);
        let fa_idx = s.targets[0].fwd[0];

        // Exactly MIN_OVERLAP bases of FA at the read start; the rest hangs
        // off the window as overhang.
        let overlap = &FA[FA.len() - MIN_OVERLAP..];

        let mut exact = overlap.to_vec();
        exact.extend(splitmix_dna(10, 400));
        let scores = scorer.score(&exact, (0, overlap.len()), End::Five, false, false);
        assert!(
            scores.iter().any(|s| s.primer == fa_idx),
            "an exact MIN_OVERLAP-base overlap scores via overhang: {scores:?}"
        );

        let mut two_edits = substituted(overlap, &[2, 7]);
        two_edits.extend(splitmix_dna(11, 400));
        let scores = scorer.score(&two_edits, (0, overlap.len()), End::Five, false, false);
        assert!(
            !scores.iter().any(|s| s.primer == fa_idx),
            "two edits inside a MIN_OVERLAP overlap exceed its partial budget: {scores:?}"
        );
    }

    #[test]
    fn clipped_primer_away_from_read_end_is_not_overhang() {
        // A lower error rate than the other tests keeps the budget tight
        // enough that no coincidental, overhang-free alignment of `FA`
        // scores anywhere else in this wider window either.
        let s = sheet(FA, RA);
        let scorer = Scorer::new(&s, 0.1, 150);
        let fa_idx = s.targets[0].fwd[0];

        let clipped = &FA[5..];
        let mut read = splitmix_dna(7, 50);
        let start = read.len();
        read.extend(clipped);
        read.extend(splitmix_dna(8, 400));

        let scores = scorer.score(
            &read,
            (start, start + clipped.len()),
            End::Five,
            false,
            false,
        );
        assert!(
            !scores.iter().any(|s| s.primer == fa_idx),
            "the same clip away from the read start is not overhang: {scores:?}"
        );
    }

    /// An open outer edge counts as a read end: the primer clipped there
    /// scores via overhang at either end, under the same rules as at the
    /// physical read end.
    #[test]
    fn clipped_primer_at_an_open_outer_edge_scores_via_overhang() {
        let s = sheet(FA, RA);
        let scorer = Scorer::new(&s, 0.1, 150);
        let fa_idx = s.targets[0].fwd[0];
        let ra_idx = s.targets[0].rev[0];

        let clipped = &FA[5..];
        let mut read = splitmix_dna(7, 50);
        let start = read.len();
        read.extend(clipped);
        read.extend(splitmix_dna(8, 400));
        let scores = scorer.score(
            &read,
            (start, start + clipped.len()),
            End::Five,
            true,
            false,
        );
        assert!(scores.iter().any(|s| s.primer == fa_idx), "{scores:?}");

        let ra_rc = reverse_complement(RA);
        let clipped = &ra_rc[..ra_rc.len() - 5];
        let mut read = splitmix_dna(12, 400);
        let start = read.len();
        read.extend(clipped);
        let end = read.len();
        read.extend(splitmix_dna(13, 50));
        let scores = scorer.score(&read, (start, end), End::Three, false, false);
        assert!(!scores.iter().any(|s| s.primer == ra_idx), "{scores:?}");
        let scores = scorer.score(&read, (start, end), End::Three, true, false);
        assert!(scores.iter().any(|s| s.primer == ra_idx), "{scores:?}");
    }

    /// The scores of one search per primer, the reference for the grouped
    /// search of `Scorer::score`; `boundary` is as there.
    fn scored_singly(
        scorer: &Scorer,
        read: &[u8],
        (start, stop): (usize, usize),
        end: End,
        (touches_start, touches_end): (bool, bool),
        boundary: bool,
    ) -> Vec<Score> {
        let mut searcher = search::new_overhang_searcher_fwd(scorer.alpha);
        let window = &read[start..stop];
        let mut out: Vec<Score> = scorer
            .entries
            .iter()
            .enumerate()
            .filter_map(|(primer, entry)| {
                let pattern = match end {
                    End::Five => &entry.forward,
                    End::Three => &entry.reverse,
                };
                let whole = if boundary {
                    entry.anchored
                } else {
                    entry.budget
                };
                search::hits(&mut searcher, pattern, window, whole.max(entry.budget))
                    .into_iter()
                    .filter(|hit| {
                        let overhang = hit.left_overhang + hit.right_overhang;
                        let overlap = pattern.len() - overhang;
                        (overhang == 0 && hit.cost <= whole)
                            || (overhang > 0
                                && hit.cost <= entry.budget
                                && (hit.left_overhang == 0 || touches_start)
                                && (hit.right_overhang == 0 || touches_end)
                                && overlap >= MIN_OVERLAP
                                && residue_within_budget(
                                    scorer.error_rate,
                                    hit.cost,
                                    hit.left_overhang,
                                    hit.right_overhang,
                                    overlap,
                                ))
                    })
                    .map(|hit| hit.cost)
                    .min()
                    .map(|cost| Score { primer, cost })
            })
            .collect();
        out.sort_by_key(|s| (s.cost, s.primer));
        out
    }

    #[test]
    fn grouped_search_scores_as_one_search_per_primer() {
        // Five forward primers of one length, one of them degenerate and
        // three a few substitutions apart, so several score in one window,
        // and a shorter sixth primer in a group of its own.
        let forwards: Vec<Vec<u8>> = vec![
            FA.to_vec(),
            substituted(FA, &[3]),
            substituted(FA, &[3, 15]),
            b"ACGGTTYAGCATTGRCCGTA".to_vec(),
            splitmix_dna(40, FA.len()),
            splitmix_dna(41, FA.len() - 1),
        ];
        let mut text = String::from("target\tfwd\trev\n");
        for (i, fwd) in forwards.iter().enumerate() {
            let rev = splitmix_dna(50 + i as u64, RA.len());
            text.push_str(&format!(
                "T{i}\t{}\t{}\n",
                String::from_utf8_lossy(fwd),
                String::from_utf8_lossy(&rev),
            ));
        }
        let s = Sheet::parse_tsv(&text).unwrap();
        let scorer = Scorer::new(&s, 0.2, 150);
        assert!(scorer.groups.iter().any(|g| g.members.len() >= 4));

        let mut compared = 0;
        for seed in 0..200u64 {
            let edits = (seed % 5) as usize;
            let positions: Vec<usize> = (0..edits)
                .map(|e| (3 + 4 * e + seed as usize) % 20)
                .collect();
            let primer = substituted(&forwards[(seed % 3) as usize], &positions);
            let cut = if seed % 3 == 0 {
                (seed % 9) as usize
            } else {
                0
            };
            for end in [End::Five, End::Three] {
                let insert = splitmix_dna(100 + seed, 300);
                let (read, locus) = match end {
                    End::Five => {
                        let mut read = primer[cut..].to_vec();
                        read.extend(&insert);
                        (read, (0, primer.len() - cut))
                    },
                    End::Three => {
                        let site = reverse_complement(&primer[cut..]);
                        let mut read = insert.clone();
                        read.extend(&site);
                        let n = read.len();
                        (read, (n - site.len(), n))
                    },
                };
                for (outer_open, boundary) in [(false, false), (true, false), (false, true)] {
                    let widen =
                        scorer.max_primer_len - (locus.1 - locus.0).min(scorer.max_primer_len) + 3;
                    let start = match end {
                        End::Five if outer_open => locus.0,
                        _ => locus.0.saturating_sub(widen),
                    };
                    let stop = match end {
                        End::Three if outer_open => locus.1,
                        _ => (locus.1 + widen).min(read.len()),
                    };
                    let touches = (
                        start == 0 || (outer_open && end == End::Five),
                        stop == read.len() || (outer_open && end == End::Three),
                    );
                    let expected =
                        scored_singly(&scorer, &read, (start, stop), end, touches, boundary);
                    let scores = scorer.score(&read, locus, end, outer_open, boundary);
                    assert_eq!(
                        scores, expected,
                        "seed {seed} {end:?} open {outer_open} boundary {boundary}"
                    );
                    compared += usize::from(!expected.is_empty());
                }
            }
        }
        assert!(compared > 300, "{compared}");
    }

    /// Linear congruential generator for the pruning fixtures.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, n: usize) -> usize {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 33) as usize) % n
        }

        fn dna(&mut self, n: usize) -> Vec<u8> {
            (0..n).map(|_| b"ACGT"[self.below(4)]).collect()
        }
    }

    /// One planted copy of `primer`: each ambiguity code resolved to one of
    /// its bases, then up to three substitutions, deletions or insertions,
    /// some beside a homopolymer run.
    fn planted(primer: &[u8], rng: &mut Lcg) -> Vec<u8> {
        let mut seq: Vec<u8> = primer
            .iter()
            .map(|&code| {
                let bases = search::iupac_bases(code).unwrap();
                bases[rng.below(bases.len())]
            })
            .collect();
        for _ in 0..rng.below(4) {
            let at = rng.below(seq.len());
            match rng.below(4) {
                0 => seq[at] = b"ACGT"[rng.below(4)],
                1 => {
                    seq.remove(at);
                },
                2 => seq.insert(at, b"ACGT"[rng.below(4)]),
                _ => seq.insert(at, seq[at]),
            }
        }
        seq
    }

    /// A sheet of IUPAC and plain primers of several lengths: five random
    /// targets with up to four ambiguity codes per primer, and two targets
    /// whose primers are a few substitutions apart.
    fn iupac_sheet() -> Sheet {
        let mut text = String::from("target\tfwd\trev\n");
        let mut rng = Lcg(0x7072_756e);
        for t in 0..5u64 {
            let len = 17 + rng.below(8);
            let mut fwd = splitmix_dna(300 + t, len);
            let mut rev = splitmix_dna(400 + t, len + 1);
            for _ in 0..t {
                fwd[rng.below(len)] = b"RYSWKMBDHVN"[rng.below(11)];
                rev[rng.below(len)] = b"RYKM"[rng.below(4)];
            }
            text.push_str(&format!(
                "T{t}\t{}\t{}\n",
                String::from_utf8_lossy(&fwd),
                String::from_utf8_lossy(&rev),
            ));
        }
        let variants = substituted(FA, &[3, 15]);
        text.push_str(&format!(
            "V\t{}\t{}\nW\t{}\t{}\n",
            String::from_utf8_lossy(FA),
            String::from_utf8_lossy(RA),
            String::from_utf8_lossy(&variants),
            String::from_utf8_lossy(&substituted(RA, &[7])),
        ));
        Sheet::parse_tsv(&text).unwrap()
    }

    /// `score_at_site` gives the scores of `score` at every whole hit of an
    /// entry, at both ends, plain, open and at a boundary: under the MAB114
    /// sheet, whose primers are degenerate and a few edits apart, and under
    /// a sheet of IUPAC and plain primers of several lengths, with entries
    /// that are primers and reverse complements of primers. The reads hold
    /// planted primers with edits, flush with a read end, behind a few bases
    /// or behind another primer cut short at the end, which scores there by
    /// overhang, and some hold `N` or lowercase bases. Every fourth read is
    /// one planted primer between a few bases, shorter than two primers, so
    /// its window is clamped at both read ends. More than one and a half
    /// primers per site and end are skipped on average.
    #[test]
    fn site_pruning_keeps_the_scores_of_the_full_search() {
        let mut rng = Lcg(0x7072_756e);
        let sheets = [Sheet::preset("mab114").unwrap(), iupac_sheet()];
        let (mut compared, mut scored, mut pruned) = (0usize, 0usize, 0usize);
        let mut clamped = 0usize;
        for sheet in &sheets {
            let mut adapters: Vec<Adapter> = sheet
                .primers
                .iter()
                .map(|p| Adapter {
                    name: p.name.clone(),
                    seq: p.seq.clone(),
                    role: crate::adapter::Role::Primer,
                })
                .collect();
            adapters.extend(sheet.primers.iter().step_by(3).map(|p| Adapter {
                name: format!("{}_rc", p.name),
                seq: reverse_complement(&p.seq),
                role: crate::adapter::Role::Primer,
            }));
            let scorer = Scorer::new(sheet, 0.2, 150).with_sites(&adapters);
            let primers: Vec<&[u8]> = sheet.primers.iter().map(|p| p.seq.as_slice()).collect();
            let mut searcher = search::new_ambiguous_searcher();
            for case in 0..400 {
                let mut pick = |rng: &mut Lcg| planted(primers[rng.below(primers.len())], rng);
                let outboard =
                    |rng: &mut Lcg, pick: &mut dyn FnMut(&mut Lcg) -> Vec<u8>| match rng.below(3) {
                        0 => Vec::new(),
                        1 => {
                            let len = 1 + rng.below(6);
                            rng.dna(len)
                        },
                        _ => {
                            let other = pick(rng);
                            other[other.len() - (5 + rng.below(10)).min(other.len())..].to_vec()
                        },
                    };
                let head = [outboard(&mut rng, &mut pick), pick(&mut rng)].concat();
                let tail = [outboard(&mut rng, &mut pick), pick(&mut rng)].concat();
                let insert = splitmix_dna(5_000 + case, 40 + rng.below(200));
                let mut read = if case % 4 == 3 {
                    let (left, right) = (rng.below(6), rng.below(6));
                    [rng.dna(left), pick(&mut rng), rng.dna(right)].concat()
                } else {
                    [head, insert, reverse_complement(&tail)].concat()
                };
                if rng.below(4) == 0 {
                    let at = rng.below(read.len());
                    read[at] = b'N';
                }
                if rng.below(5) == 0 {
                    let at = rng.below(read.len());
                    let to = (at + 30).min(read.len());
                    read[at..to].make_ascii_lowercase();
                }
                let mut buf = Vec::new();
                let (normalized, _) = normalize_into(&read, &mut buf);
                for (entry, adapter) in adapters.iter().enumerate() {
                    for hit in search::hits(&mut searcher, &adapter.seq, normalized, 5) {
                        let site = PrimerSite {
                            entry,
                            rc: hit.rc,
                            cost: hit.cost,
                        };
                        for end in [End::Five, End::Three] {
                            for (outer_open, boundary) in
                                [(false, false), (false, true), (true, false)]
                            {
                                let locus = (hit.start, hit.end);
                                let expected =
                                    scorer.score(&read, locus, end, outer_open, boundary);
                                let at_site = scorer.score_at_site(
                                    &read,
                                    locus,
                                    end,
                                    outer_open,
                                    boundary,
                                    Some(site),
                                );
                                assert_eq!(
                                    at_site, expected,
                                    "case {case} {site:?} {locus:?} {end:?} \
                                     open {outer_open} boundary {boundary}"
                                );
                                compared += 1;
                                scored += usize::from(!expected.is_empty());
                            }
                            let widen =
                                scorer.max_primer_len.saturating_sub(hit.end - hit.start) + 3;
                            let start = hit.start.saturating_sub(widen);
                            let stop = (hit.end + widen).min(read.len());
                            clamped += usize::from(start == 0 && stop == read.len());
                            let flank = |touches: bool, bases: usize| {
                                bases + if touches { scorer.max_primer_len } else { 0 }
                            };
                            let bounds = scorer
                                .site_bounds(
                                    (entry, hit.rc),
                                    end,
                                    flank(start == 0, hit.start - start),
                                    flank(stop == read.len(), stop - hit.end),
                                )
                                .unwrap();
                            pruned += scorer
                                .entries
                                .iter()
                                .zip(bounds.iter())
                                .filter(|&(e, &b)| b > hit.cost + e.budget.max(e.anchored))
                                .count();
                        }
                    }
                }
            }
        }
        assert!(
            compared > 20_000 && scored > 10_000 && 2 * pruned > compared && clamped > 500,
            "{compared} loci, {scored} scored, {pruned} primers skipped, \
             {clamped} sites clamped at both ends"
        );
    }

    /// Every site `adapter_segments_annotated` records aligns its entry over
    /// its locus within its cost (`site_holds`), and `score_at_site` gives
    /// the scores of `score` there: under the MAB114 kit with its sheet, and
    /// under the ligation adapters with `iupac_sheet`. The reads are
    /// amplicons of planted primers behind whole, cut or absent outer
    /// layers, so some primers sit at a trim boundary, chimeras of two
    /// amplicons in either orientation, and short reads of one primer
    /// between a few bases; some hold `N` bases.
    #[test]
    fn recorded_sites_align_over_their_loci() {
        use crate::adapter::preset::{Kit, preset};
        use crate::adapter::{AdapterConfig, Locus, adapter_segments_annotated, site_holds};

        let exact = Locus {
            start: 2,
            end: 2 + FA.len(),
            outer_open: false,
            boundary: false,
            site: Some(PrimerSite {
                entry: 0,
                rc: false,
                cost: 0,
            }),
        };
        let entry = [Adapter {
            name: "fA".into(),
            seq: FA.to_vec(),
            role: crate::adapter::Role::Primer,
        }];
        assert!(site_holds(&[b"GG", FA].concat(), &entry, &exact));
        assert!(!site_holds(
            &[b"GG".as_slice(), &substituted(FA, &[5])].concat(),
            &entry,
            &exact
        ));

        let named = |adapters: &[Adapter], name: &str| {
            adapters
                .iter()
                .find(|a| a.name == name)
                .unwrap()
                .seq
                .clone()
        };
        let mab = preset(&[Kit::Mab114]);
        let mab_layers: Vec<Vec<u8>> = ["TP03", "TP17"]
            .iter()
            .map(|bc| {
                [
                    named(&mab, "RBK4_front"),
                    named(&mab, bc),
                    named(&mab, "MAB_rear"),
                ]
                .concat()
            })
            .collect();
        let mab_tails: Vec<Vec<u8>> = mab_layers.iter().map(|l| reverse_complement(l)).collect();
        let lsk = preset(&[Kit::Lsk114]);
        let lsk_layers = vec![named(&lsk, "LSK114_front"), named(&lsk, "LSK109_front")];
        let lsk_tails = vec![named(&lsk, "LSK114_rear"), named(&lsk, "LSK109_rear")];
        let sets = [
            (
                "mab114",
                mab,
                Sheet::preset("mab114").unwrap(),
                mab_layers,
                mab_tails,
                true,
            ),
            ("iupac", lsk, iupac_sheet(), lsk_layers, lsk_tails, false),
        ];

        let mut rng = Lcg(0x7369_7465);
        for (label, adapters, sheet, heads, tails, amplicon) in sets {
            let mut cfg = AdapterConfig {
                adapters,
                error_rate: 0.2,
                end_size: 150,
                split: true,
                min_piece: 1,
                candidate_index: OnceLock::new(),
                amplicon,
                split_of: Vec::new(),
                split_opens: Vec::new(),
            };
            cfg.attach_split(&sheet.primers);
            let scorer = Scorer::new(&sheet, 0.2, 150).with_sites(&cfg.adapters);
            // A third of the primers carry as many substitutions as their
            // anchored budget, more than their terminal budget, so behind an
            // outer layer they locate only at its trim boundary.
            let anchored = scorer.anchored_budgets();
            let primer = |rng: &mut Lcg, list: &[usize]| {
                let index = list[rng.below(list.len())];
                let codes = &sheet.primers[index].seq;
                if rng.below(3) != 0 {
                    return planted(codes, rng);
                }
                let mut seq: Vec<u8> = codes
                    .iter()
                    .map(|&code| {
                        let bases = search::iupac_bases(code).unwrap();
                        bases[rng.below(bases.len())]
                    })
                    .collect();
                let step = seq.len() / (anchored[index] + 1);
                for at in (1..=anchored[index]).map(|i| i * step) {
                    seq[at] = if seq[at] == b'A' { b'C' } else { b'A' };
                }
                seq
            };
            let layer = |rng: &mut Lcg, layers: &[Vec<u8>], five: bool| {
                let layer = &layers[rng.below(layers.len())];
                let keep = 5 + rng.below(layer.len() - 5);
                match rng.below(3) {
                    0 => Vec::new(),
                    1 => layer.clone(),
                    _ if five => layer[layer.len() - keep..].to_vec(),
                    _ => layer[..keep].to_vec(),
                }
            };
            let amplicon = |rng: &mut Lcg, seed: u64| {
                let target = &sheet.targets[rng.below(sheet.targets.len())];
                let read = [
                    layer(rng, &heads, true),
                    primer(rng, &target.fwd),
                    splitmix_dna(seed, 150 + rng.below(300)),
                    reverse_complement(&primer(rng, &target.rev)),
                    layer(rng, &tails, false),
                ]
                .concat();
                if rng.below(2) == 0 {
                    reverse_complement(&read)
                } else {
                    read
                }
            };
            let (mut sites, mut boundary, mut chimeric, mut short, mut with_n) =
                (0usize, 0usize, 0usize, 0usize, 0usize);
            for case in 0..500u64 {
                let kind = case % 4;
                let mut read = match kind {
                    0 | 1 => amplicon(&mut rng, 9_000 + case),
                    2 => [
                        amplicon(&mut rng, 9_000 + case),
                        amplicon(&mut rng, 9_500 + case),
                    ]
                    .concat(),
                    _ => {
                        let target = &sheet.targets[rng.below(sheet.targets.len())];
                        let (left, right) = (rng.below(6), rng.below(6));
                        [rng.dna(left), primer(&mut rng, &target.fwd), rng.dna(right)].concat()
                    },
                };
                let has_n = rng.below(3) == 0;
                if has_n {
                    for _ in 0..1 + rng.below(3) {
                        let at = if rng.below(2) == 0 {
                            rng.below(read.len().min(80))
                        } else {
                            read.len() - 1 - rng.below(read.len().min(80))
                        };
                        read[at] = b'N';
                    }
                }
                let segments = adapter_segments_annotated(&read, &cfg);
                for seg in &segments {
                    for (locus, end) in [(seg.five, End::Five), (seg.three, End::Three)] {
                        let Some(locus) = locus else {
                            continue;
                        };
                        assert!(
                            site_holds(&read, &cfg.adapters, &locus),
                            "{label} case {case} {locus:?}"
                        );
                        let Some(site) = locus.site else {
                            continue;
                        };
                        let span = (locus.start, locus.end);
                        let (open, at) = (locus.outer_open, locus.boundary);
                        assert_eq!(
                            scorer.score_at_site(&read, span, end, open, at, Some(site)),
                            scorer.score(&read, span, end, open, at),
                            "{label} case {case} {locus:?} {end:?}"
                        );
                        sites += 1;
                        boundary += usize::from(locus.boundary);
                        chimeric += usize::from(kind == 2 && segments.len() > 1);
                        short += usize::from(kind == 3);
                        with_n += usize::from(has_n);
                    }
                }
            }
            assert!(
                sites > 500 && boundary > 5 && chimeric > 100 && short > 25 && with_n > 100,
                "{label}: {sites} sites, {boundary} at a boundary, {chimeric} in split \
                 chimeras, {short} in short reads, {with_n} in reads with N"
            );
        }
    }

    #[test]
    fn scores_sorted_by_cost() {
        // A second target whose forward primer is `FA` with two
        // substitutions, so both primers are found in the same window at
        // different costs.
        let fb = substituted(FA, &[3, 15]);
        let text = format!(
            "target\tfwd\trev\n\
             A\t{}\t{}\n\
             B\t{}\t{}\n",
            String::from_utf8_lossy(FA),
            String::from_utf8_lossy(RA),
            String::from_utf8_lossy(&fb),
            String::from_utf8_lossy(RA),
        );
        let s = Sheet::parse_tsv(&text).unwrap();
        let scorer = Scorer::new(&s, 0.2, 150);
        let fa_idx = s.targets[0].fwd[0];
        let fb_idx = s.targets[1].fwd[0];

        let mut read = FA.to_vec();
        read.extend(splitmix_dna(9, 400));

        let scores = scorer.score(&read, (0, FA.len()), End::Five, false, false);
        let costs: Vec<usize> = scores.iter().map(|s| s.cost).collect();
        let mut sorted = costs.clone();
        sorted.sort_unstable();
        assert_eq!(costs, sorted, "{scores:?}");
        assert_eq!(scores.iter().find(|s| s.primer == fa_idx).unwrap().cost, 0);
        assert_eq!(scores.iter().find(|s| s.primer == fb_idx).unwrap().cost, 2);
    }
}
