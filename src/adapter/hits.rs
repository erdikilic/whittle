//! Classification and accumulation of adapter hits: terminal trims, interior
//! excisions, the anchoring of deferred terminal hits, and the cuts of one
//! window.

use super::*;

/// Terminal classification of a hit: which end, if any, it trims.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Terminal {
    /// Trims the 5' end.
    Five,
    /// Trims the 3' end.
    Three,
    /// Covered by both end zones with a real flank on each side: excise the
    /// adapter span and keep both flanks. See `classify_terminal`.
    Excise,
    /// Trims neither end.
    None,
}

/// Emits a trace event for one adapter hit: the sequence, its span, its edit
/// cost, and the action taken (a terminal trim, an excision, or none).
pub(super) fn trace_hit(
    name: &str,
    start: usize,
    end: usize,
    cost: usize,
    action: Option<HitAction>,
) {
    tracing::trace!(
        adapter = name,
        start,
        end,
        cost,
        action = match action {
            Some(HitAction::TrimFivePrime) => "trim 5'",
            Some(HitAction::TrimThreePrime) => "trim 3'",
            Some(HitAction::Excise) => "excise and split",
            None => "no action",
        },
        "Adapter hit"
    );
}

/// What the trimmer did with an accepted adapter hit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum HitAction {
    /// Terminal hit at the 5' end: the keep-boundary moved inward past it.
    TrimFivePrime,
    /// Terminal hit at the 3' end.
    TrimThreePrime,
    /// Interior hit: the span is cut out and both flanks are kept.
    Excise,
}

/// Returns the end nearer to a hit at `[start, end)` in a length-`n` window.
pub(super) fn nearer_end(start: usize, end: usize, n: usize) -> Terminal {
    if start <= n - end {
        Terminal::Five
    } else {
        Terminal::Three
    }
}

/// Classifies a hit at window coordinates `[start, end)` in a length-`n`
/// window by geometry alone.
///
/// A hit inside one end zone trims that end. Every search is
/// reverse-complement aware and a rear adapter is the reverse complement of
/// its front adapter, so which catalog entry matched says nothing about which
/// end the hit is at. A hit covered by both end zones (`n <= end_size + hit
/// length`) trims toward an end whose outboard flank is at most `FLANK_SLACK`
/// (the nearer one when both are), and otherwise becomes `Excise`: cut out the
/// adapter, keep both flanks, as a central chimera junction needs.
pub(super) fn classify_terminal(start: usize, end: usize, n: usize, end_size: usize) -> Terminal {
    let in_head = start <= end_size;
    let in_tail = end >= n.saturating_sub(end_size);
    match (in_head, in_tail) {
        (true, true) => match (start <= FLANK_SLACK, n - end <= FLANK_SLACK) {
            (false, false) => Terminal::Excise,
            (true, false) => Terminal::Five,
            (false, true) => Terminal::Three,
            (true, true) => nearer_end(start, end, n),
        },
        (true, false) => Terminal::Five,
        (false, true) => Terminal::Three,
        (false, false) => Terminal::None,
    }
}

/// Classifies a hit for ends-only mode: splitting is disabled, so an `Excise`
/// outcome of `classify_terminal` resolves to a terminal trim toward the nearer
/// end, and every other outcome is unchanged.
pub(super) fn ends_only_terminal(start: usize, end: usize, n: usize, end_size: usize) -> Terminal {
    match classify_terminal(start, end, n, end_size) {
        Terminal::Excise => nearer_end(start, end, n),
        other => other,
    }
}

/// Where a hit was found. Each site acts on the outcomes it owns, so a hit
/// that the head and tail windows both contain (they overlap on a read shorter
/// than twice their reach) is applied and reported once.
///
/// Every 5' trim and every excision lies inside the head window and every 3'
/// trim inside the tail window: a hit reaching an end zone starts or ends within
/// `end_size` of that end and spans at most `len + k_end` bases.
#[derive(Debug, Clone, Copy)]
pub(super) enum Site {
    /// `[0, end_size + len + k_end)`: owns 5' trims and excisions.
    Head,
    /// `[n - (end_size + len + k_end), n)`: owns 3' trims. `head_end` is where
    /// the head window stopped, so hits before it are traced there only.
    Tail { head_end: usize },
    /// A candidate window searched at `k_mid`: owns interior excisions.
    Interior,
}

/// Returns the head end and tail start for a `len`-base adapter at budget
/// `k_end`. Each window covers `end_size` bases plus the longest alignment,
/// `len + k_end`, so it holds every hit that can reach its end zone.
pub(super) fn terminal_windows(
    n: usize,
    end_size: usize,
    len: usize,
    k_end: usize,
) -> (usize, usize) {
    let reach = end_size + len + k_end;
    (reach.min(n), n.saturating_sub(reach))
}

/// Accumulator for the accepted hits of one window: the keep boundaries and
/// interior cuts. Coordinates are those of the window, `[0, n)`.
pub(super) struct Keep<'a> {
    /// The configured adapters, for roles and names.
    pub(super) adapters: &'a [Adapter],
    /// Per-adapter edit budgets, for the anchoring of terminal hits.
    pub(super) budgets: &'a [Budget],
    /// The run's error rate, which scales the budget of a partial hit.
    pub(super) error_rate: f64,
    /// Window length.
    pub(super) n: usize,
    /// Terminal zone depth, capped at `n`.
    pub(super) end_size: usize,
    /// Whether interior hits split the read.
    pub(super) split: bool,
    /// 5' keep boundary; advances inward on 5' trims.
    pub(super) lo: usize,
    /// 3' keep boundary; retreats inward on 3' trims.
    pub(super) hi: usize,
    /// Accepted excisions, merged by `into_cuts`.
    pub(super) interior: Vec<(usize, usize)>,
    /// The hit behind each excision in `interior`, for `refine`.
    excised: Vec<Applied>,
    /// The hits behind the applied 5' trims, for `refine`.
    five: Vec<Applied>,
    /// The hits behind the applied 3' trims, for `refine`.
    three: Vec<Applied>,
    /// 5' boundary set by trims that `refine` leaves as found.
    fixed_lo: usize,
    /// 3' boundary set by trims that `refine` leaves as found.
    fixed_hi: usize,
    /// The adapters whose hits trimmed or excised, for presence detection.
    pub(super) acted: Vec<usize>,
    /// Terminal trims above the `k_far` budget of their adapter, applied by
    /// `settle` once anchored.
    pub(super) deferred: Vec<Deferred>,
}

/// An applied hit, kept for `refine`: its span and the clips of its
/// alignment ends.
#[derive(Debug, Clone, Copy)]
pub(super) struct Applied {
    /// Hit start in window coordinates.
    start: usize,
    /// Hit end in window coordinates.
    end: usize,
    /// See `Hit::clip_start`.
    clip_start: usize,
    /// See `Hit::clip_end`.
    clip_end: usize,
}

/// A terminal trim held until it is anchored at the read end or at an
/// accepted hit.
#[derive(Debug, Clone, Copy)]
pub(super) struct Deferred {
    /// Index into the configured adapters.
    pub(super) adapter_idx: usize,
    /// Hit start in window coordinates.
    pub(super) start: usize,
    /// Hit end in window coordinates.
    pub(super) end: usize,
    /// Edit cost of the hit.
    pub(super) cost: usize,
    /// See `Hit::clip_start`.
    pub(super) clip_start: usize,
    /// See `Hit::clip_end`.
    pub(super) clip_end: usize,
    /// `TrimFivePrime` or `TrimThreePrime`.
    pub(super) action: HitAction,
}

impl<'a> Keep<'a> {
    /// Creates an accumulator that keeps the whole `[0, n)` window. `split`
    /// selects the classification: with it, a hit covered by both end zones
    /// may excise; without it, every hit trims an end.
    pub(super) fn new(
        cfg: &'a AdapterConfig,
        index: &'a CandidateIndex,
        n: usize,
        split: bool,
    ) -> Self {
        Self {
            adapters: &cfg.adapters,
            budgets: &index.budgets,
            error_rate: cfg.error_rate,
            n,
            end_size: cfg.end_size.min(n),
            split,
            lo: 0,
            hi: n,
            interior: Vec::new(),
            excised: Vec::new(),
            five: Vec::new(),
            three: Vec::new(),
            fixed_lo: 0,
            fixed_hi: n,
            acted: Vec::new(),
            deferred: Vec::new(),
        }
    }

    /// Returns whether a partial hit is acceptable: enough of the pattern
    /// aligned, flush with the window end it hangs off, and within the edit
    /// budget of the aligned part alone. A whole-pattern hit always passes.
    pub(super) fn partial_hit_is_valid(&self, adapter_idx: usize, hit: Hit) -> bool {
        let overhang = hit.left_overhang + hit.right_overhang;
        if overhang == 0 {
            return true;
        }
        let overlap = self.adapters[adapter_idx].seq.len() - overhang;
        overlap >= MIN_OVERLAP
            && (hit.left_overhang == 0 || hit.start == 0)
            && (hit.right_overhang == 0 || hit.end == self.n)
            && self.residue_within_budget(hit, overlap)
    }

    /// Returns whether a hit whose `overlap` bases lie in the text is within
    /// the partial budget of that overlap, once the cost of the pattern bases
    /// beyond the text end is discounted.
    pub(super) fn residue_within_budget(&self, hit: Hit, overlap: usize) -> bool {
        let charged = overhang_cost(
            self.error_rate as f32,
            hit.left_overhang,
            hit.right_overhang,
        );
        hit.cost.saturating_sub(charged) <= partial_budget(self.error_rate, overlap)
    }

    /// Classifies one hit and applies it when `site` owns the outcome.
    pub(super) fn accept(&mut self, site: Site, adapter_idx: usize, hit: Hit) {
        let adapter = &self.adapters[adapter_idx];
        let Hit {
            start, end, cost, ..
        } = hit;
        if !self.partial_hit_is_valid(adapter_idx, hit) {
            return;
        }
        let terminal = if self.split {
            classify_terminal(start, end, self.n, self.end_size)
        } else {
            ends_only_terminal(start, end, self.n, self.end_size)
        };
        let action = match (site, terminal) {
            (Site::Head, Terminal::Five) => HitAction::TrimFivePrime,
            (Site::Head, Terminal::Excise) if adapter.role.splits() => HitAction::Excise,
            (Site::Head, Terminal::Excise) => match nearer_end(start, end, self.n) {
                Terminal::Five => HitAction::TrimFivePrime,
                _ => HitAction::TrimThreePrime,
            },
            (Site::Tail { .. }, Terminal::Three) => HitAction::TrimThreePrime,
            // An interior hit within the flank slack of an end is residue of
            // that end rather than a junction.
            (Site::Interior, _) if start <= FLANK_SLACK => HitAction::TrimFivePrime,
            (Site::Interior, _) if self.n - end <= FLANK_SLACK => HitAction::TrimThreePrime,
            (Site::Interior, Terminal::None) if adapter.role.splits() => HitAction::Excise,
            (Site::Interior, Terminal::None) => {
                trace_hit(&adapter.name, start, end, cost, None);
                return;
            },
            (Site::Head, Terminal::None) => {
                trace_hit(&adapter.name, start, end, cost, None);
                return;
            },
            (Site::Tail { head_end }, Terminal::None) => {
                if start >= head_end {
                    trace_hit(&adapter.name, start, end, cost, None);
                }
                return;
            },
            _ => return,
        };
        let whole = hit.left_overhang + hit.right_overhang == 0;
        let terminal = matches!(site, Site::Head | Site::Tail { .. })
            && matches!(action, HitAction::TrimFivePrime | HitAction::TrimThreePrime);
        if whole && terminal && cost > self.budgets[adapter_idx].k_far {
            self.deferred.push(Deferred {
                adapter_idx,
                start,
                end,
                cost,
                clip_start: hit.clip_start,
                clip_end: hit.clip_end,
                action,
            });
            return;
        }
        self.apply(
            adapter_idx,
            start,
            end,
            cost,
            (hit.clip_start, hit.clip_end),
            action,
        );
    }

    /// Applies an accepted hit to the keep boundaries or the excisions.
    pub(super) fn apply(
        &mut self,
        adapter_idx: usize,
        start: usize,
        end: usize,
        cost: usize,
        (clip_start, clip_end): (usize, usize),
        action: HitAction,
    ) {
        trace_hit(
            &self.adapters[adapter_idx].name,
            start,
            end,
            cost,
            Some(action),
        );
        self.acted.push(adapter_idx);
        let applied = Applied {
            start,
            end,
            clip_start,
            clip_end,
        };
        match action {
            HitAction::TrimFivePrime => {
                self.lo = self.lo.max(end);
                self.five.push(applied);
            },
            HitAction::TrimThreePrime => {
                self.hi = self.hi.min(start);
                self.three.push(applied);
            },
            HitAction::Excise => {
                self.interior.push((start, end));
                self.excised.push(applied);
            },
        }
    }

    /// Trims the 5' end to `end` at a boundary `refine` keeps as found, for
    /// hits aligned against a masked copy of the window.
    pub(super) fn trim_five_fixed(&mut self, end: usize) {
        self.lo = self.lo.max(end);
        self.fixed_lo = self.fixed_lo.max(end);
    }

    /// Trims the 3' end to `start` at a boundary `refine` keeps as found.
    pub(super) fn trim_three_fixed(&mut self, start: usize) {
        self.hi = self.hi.min(start);
        self.fixed_hi = self.fixed_hi.min(start);
    }

    /// Moves every boundary to the inner end of the alignment that set it,
    /// dropping the clipped bases of the hit behind it; see `refine`. A
    /// boundary never moves outward. The 5' trims are visited in decreasing end
    /// order until a raw end falls at or below the best clipped one, and the 3'
    /// trims likewise.
    pub(super) fn refine(&mut self) {
        if !self.five.is_empty() {
            self.five.sort_unstable_by_key(|a| std::cmp::Reverse(a.end));
            let mut lo = self.fixed_lo;
            for a in &self.five {
                if a.end <= lo {
                    break;
                }
                lo = lo.max(a.end - a.clip_end);
            }
            self.lo = lo;
        }
        if !self.three.is_empty() {
            self.three.sort_unstable_by_key(|a| a.start);
            let mut hi = self.fixed_hi;
            for a in &self.three {
                if a.start >= hi {
                    break;
                }
                hi = hi.min(a.start + a.clip_start);
            }
            self.hi = hi;
        }
        for (cut, a) in self.interior.iter_mut().zip(&self.excised) {
            *cut = (a.start + a.clip_start, a.end - a.clip_end);
        }
    }

    /// Applies each deferred trim that is anchored: a 5' hit starting within
    /// `FLANK_SLACK` of the 5' keep boundary, or a 3' hit ending within
    /// `FLANK_SLACK` of the 3' boundary. An applied trim moves the boundary,
    /// which may anchor further deferred trims behind it.
    pub(super) fn settle(&mut self) {
        loop {
            let (lo, hi) = (self.lo, self.hi);
            let anchored = self.deferred.iter().position(|d| match d.action {
                HitAction::TrimFivePrime => d.start <= lo + FLANK_SLACK && d.end > lo,
                HitAction::TrimThreePrime => d.end + FLANK_SLACK >= hi && d.start < hi,
                HitAction::Excise => false,
            });
            let Some(i) = anchored else {
                return;
            };
            let d = self.deferred.swap_remove(i);
            self.apply(
                d.adapter_idx,
                d.start,
                d.end,
                d.cost,
                (d.clip_start, d.clip_end),
                d.action,
            );
        }
    }

    /// Returns the keep boundaries and the excisions, refined and clipped to
    /// the boundaries, merged: two excisions overlapping,
    /// touching, or separated by at most `FLANK_SLACK` bases or fewer than
    /// `min_piece` bases become one, since the bases between them are junction
    /// residue or a piece the length filter would discard.
    pub(super) fn into_cuts(mut self, min_piece: usize) -> (usize, usize, Vec<(usize, usize)>) {
        self.settle();
        self.refine();
        for d in &self.deferred {
            trace_hit(
                &self.adapters[d.adapter_idx].name,
                d.start,
                d.end,
                d.cost,
                None,
            );
        }
        let Keep {
            lo, hi, interior, ..
        } = self;
        if lo >= hi {
            return (lo, hi, Vec::new());
        }
        let mut cuts: Vec<(usize, usize)> = interior
            .into_iter()
            .filter_map(|(s, e)| {
                let s = s.max(lo);
                let e = e.min(hi);
                (s < e).then_some((s, e))
            })
            .collect();
        cuts.sort_unstable();
        let mut merged: Vec<(usize, usize)> = Vec::new();
        for (s, e) in cuts {
            if let Some(last) = merged.last_mut() {
                let gap = s.saturating_sub(last.1);
                if s <= last.1 || gap <= FLANK_SLACK || gap < min_piece {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }
        (lo, hi, merged)
    }
}
