//! The candidate index of an adapter set: per-adapter budgets, exact partition
//! seeds that bound the interior search, end seeds that gate the partial
//! search, and the equal-length batches of the terminal search.

use super::*;

/// Upper bound on the plain strings one seed piece may expand to. A piece past
/// it marks its adapter `unfiltered`, and the interior search covers the whole
/// read for that adapter instead of candidate windows.
pub(super) const MAX_SEED_EXPANSIONS: usize = 256;

/// Exact-seed index over the adapter set. Partition seeds, looked up through a
/// prefix table, bound the interior search to candidate windows, and
/// equal-length barcode entries are grouped into SIMD batches for the terminal
/// search.
#[derive(Debug, Clone)]
pub(crate) struct CandidateIndex {
    /// Table over the interior seeds of every splitting adapter; `None` when
    /// no adapter has seeds.
    pub(super) seeds: Option<SeedTable>,
    /// Per-adapter edit budgets, computed once per adapter set.
    pub(super) budgets: Vec<Budget>,
    /// Per-adapter `is_plain_acgt`, which selects the search profile.
    pub(super) plain: Vec<bool>,
    /// Adapters with no usable seeds (see `MAX_SEED_EXPANSIONS`).
    pub(super) unfiltered: Vec<bool>,
    /// Equal-length adapter groups searched together over the end windows.
    pub(super) terminal_batches: Vec<TerminalBatch>,
    /// Adapters searched one pattern at a time over the end windows: those of
    /// searchable length that no batch covers.
    pub(super) singletons: Vec<bool>,
    /// Table over every `END_SEED_LEN`-mer of the entries eligible for
    /// partial matching, both strands; `None` when there are none. Gates the
    /// overhang search of an end to the entries with an exact seed in it.
    pub(super) end_seeds: Option<SeedTable>,
    /// The longest end-window reach over the partial-matching entries.
    pub(super) end_reach: usize,
}

/// Equal-length adapters searched together through sassy's pattern-parallel
/// API, which matches the forward and reverse-complement strands for the batch.
#[derive(Debug, Clone)]
pub(super) struct TerminalBatch {
    /// Indices into `AdapterConfig::adapters`, in pattern order.
    pub(super) adapter_indices: Vec<usize>,
    /// The adapter sequences encoded once for the tiled search.
    pub(super) encoded: EncodedAdapterBatch,
    /// The shared pattern length.
    pub(super) len: usize,
    /// The largest terminal edit budget of the batch; a hit above the budget
    /// of its own entry is discarded.
    pub(super) k_end: usize,
}

impl CandidateIndex {
    /// Builds the index for `adapters` searched over end zones of `end_size`
    /// bases; interior seeds are built only when `include_interior`, and
    /// only for roles that split.
    pub(super) fn new(
        adapters: &[Adapter],
        error_rate: f64,
        end_size: usize,
        include_interior: bool,
    ) -> Self {
        // A pattern below `MIN_PATTERN_LEN` takes part in no search: it gets
        // no seeds, no batch and no singleton search.
        let searchable: Vec<bool> = adapters
            .iter()
            .map(|adapter| adapter.seq.len() >= MIN_PATTERN_LEN)
            .collect();
        let mut budgets: Vec<Budget> = adapters
            .iter()
            .map(|adapter| Budget::new(&adapter.seq, error_rate, end_size))
            .collect();
        // A hit anchored at the read end or at an accepted hit starts within
        // `FLANK_SLACK` of it; a hit anywhere in the zone may start at any of
        // its `end_size + 1` positions. Each count covers both strands and
        // both ends.
        let caps: Vec<usize> = budgets.iter().map(|b| b.k_end).collect();
        let near = family_budgets(
            adapters,
            &searchable,
            &caps,
            4.0 * (FLANK_SLACK + 1) as f64,
            TERMINAL_CHANCE_HITS_PER_READ,
        );
        let far = family_budgets(
            adapters,
            &searchable,
            &near,
            4.0 * (end_size + 1) as f64,
            TERMINAL_CHANCE_HITS_PER_READ,
        );
        for ((budget, k_end), k_far) in budgets.iter_mut().zip(near).zip(far) {
            budget.k_end = k_end;
            budget.k_far = k_far;
        }
        // Interior hits split only for splitting roles, which are the set the
        // interior chance bound covers, per read-length class. Budgets do not
        // increase with the class.
        let splitting: Vec<bool> = adapters
            .iter()
            .zip(&searchable)
            .map(|(adapter, &searchable)| searchable && adapter.role.splits())
            .collect();
        for class in 0..INTERIOR_CLASSES {
            let caps: Vec<usize> = budgets.iter().map(|b| b.k_mid[class]).collect();
            let edits = family_budgets(
                adapters,
                &splitting,
                &caps,
                interior_positions(class),
                INTERIOR_CHANCE_HITS_PER_READ,
            );
            for (budget, k) in budgets.iter_mut().zip(edits) {
                let previous = class.checked_sub(1).map_or(usize::MAX, |c| budget.k_mid[c]);
                budget.k_mid[class] = k.min(previous);
            }
        }
        let plain: Vec<bool> = adapters
            .iter()
            .map(|adapter| is_plain_acgt(&adapter.seq))
            .collect();
        let mut unfiltered = vec![false; adapters.len()];
        let seeds = if include_interior {
            let mut seeds: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
            for (adapter_idx, adapter) in adapters.iter().enumerate() {
                let k_mid = budgets[adapter_idx].interior_max();
                if !searchable[adapter_idx] || !adapter.role.splits() {
                    continue;
                }
                let pattern = adapter.seq.to_ascii_uppercase();
                let forward = partition_seeds(&pattern, k_mid);
                let reverse = partition_seeds(&reverse_complement(&pattern), k_mid);
                match (forward, reverse) {
                    (Some(forward), Some(reverse)) => {
                        for seed in forward.into_iter().chain(reverse) {
                            let owners = seeds.entry(seed).or_default();
                            if owners.last() != Some(&adapter_idx) {
                                owners.push(adapter_idx);
                            }
                        }
                    },
                    _ => unfiltered[adapter_idx] = true,
                }
            }

            SeedTable::new(seeds)
        } else {
            None
        };

        // Sassy packs equal-length patterns across SIMD lanes, one pattern per
        // 64-bit limb. Only the roles without overhang alignment are batched:
        // the tiled search is a whole-pattern search, and the barcode sets are
        // where equal lengths occur in numbers. Singletons stay on the
        // ordinary search path: a batch of one has no pattern-level
        // parallelism and is slower over these terminal windows.
        let mut by_len: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (adapter_idx, adapter) in adapters.iter().enumerate() {
            let len = adapter.seq.len();
            if searchable[adapter_idx] && len <= MAX_TILED_PATTERN_LEN && !adapter.role.overhangs()
            {
                by_len.entry(len).or_default().push(adapter_idx);
            }
        }
        let mut terminal_batches = Vec::new();
        let mut singletons = searchable.clone();
        for (len, adapter_indices) in by_len {
            if adapter_indices.len() >= 2 {
                let patterns: Vec<Vec<u8>> = adapter_indices
                    .iter()
                    .map(|&idx| adapters[idx].seq.clone())
                    .collect();
                for &adapter_idx in &adapter_indices {
                    singletons[adapter_idx] = false;
                }
                terminal_batches.push(TerminalBatch {
                    k_end: adapter_indices
                        .iter()
                        .map(|&idx| budgets[idx].k_end)
                        .max()
                        .unwrap_or(0),
                    adapter_indices,
                    encoded: encode_patterns(&patterns),
                    len,
                });
            }
        }

        // End seeds: every k-mer of every partial-matching entry on both
        // strands. A partial hit within its budget keeps at least one intact
        // k-mer (see `END_SEED_LEN`), so an end window without a seed of an
        // entry cannot hold a partial hit of it and skips the overhang search.
        let mut end_seeds: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
        let mut end_reach = 0;
        for (adapter_idx, adapter) in adapters.iter().enumerate() {
            let Budget { len, k_end, .. } = budgets[adapter_idx];
            if !searchable[adapter_idx] || !adapter.role.overhangs() {
                continue;
            }
            end_reach = end_reach.max(len + k_end);
            let pattern = adapter.seq.to_ascii_uppercase();
            for strand in [pattern.clone(), reverse_complement(&pattern)] {
                for kmer in strand.windows(END_SEED_LEN) {
                    let Some(expansions) = expand_iupac(kmer) else {
                        continue;
                    };
                    for seed in expansions {
                        let owners = end_seeds.entry(seed).or_default();
                        if owners.last() != Some(&adapter_idx) {
                            owners.push(adapter_idx);
                        }
                    }
                }
            }
        }
        let end_seeds = SeedTable::new(end_seeds);

        Self {
            seeds,
            budgets,
            plain,
            unfiltered,
            terminal_batches,
            singletons,
            end_seeds,
            end_reach,
        }
    }

    /// Marks in `head` and `tail` the entries with an exact end seed inside
    /// the head or tail window of `window[ws..we]`, each window reaching
    /// `end_size + end_reach` bases in from its end.
    pub(super) fn end_candidates(
        &self,
        window: &[u8],
        ws: usize,
        we: usize,
        end_size: usize,
        head: &mut Vec<bool>,
        tail: &mut Vec<bool>,
    ) {
        let adapters = self.budgets.len();
        head.clear();
        head.resize(adapters, false);
        tail.clear();
        tail.resize(adapters, false);
        let Some(table) = &self.end_seeds else {
            return;
        };
        let n = we - ws;
        let reach = (end_size + self.end_reach).min(n);
        table.scan(&window[ws..ws + reach], |adapter_idx, _| {
            head[adapter_idx] = true
        });
        table.scan(&window[we - reach..we], |adapter_idx, _| {
            tail[adapter_idx] = true
        });
    }

    /// Fills `windows` with the text spans that can hold an interior hit, as
    /// `(adapter, start, end)` sorted by adapter then start: a radius around
    /// every exact seed occurrence, merged per adapter, or the whole text for
    /// an `unfiltered` adapter. `text` is normalized (see `normalize_into`);
    /// the seed table encodes uppercase bases only.
    pub(super) fn candidate_windows(&self, text: &[u8], windows: &mut Vec<(usize, usize, usize)>) {
        windows.clear();
        for (adapter_idx, &whole) in self.unfiltered.iter().enumerate() {
            if whole {
                windows.push((adapter_idx, 0, text.len()));
            }
        }
        if let Some(table) = &self.seeds {
            table.scan(text, |adapter_idx, (start, end)| {
                let Budget { len, k_end, .. } = self.budgets[adapter_idx];
                // The exact seed lies inside the `<= k_mid` alignment. A
                // radius of pattern length + `k_end` on each side contains
                // that entire alignment and enough context for the
                // full-window `k_end` search to reproduce its span and tie
                // behavior.
                let radius = len + k_end;
                windows.push((
                    adapter_idx,
                    start.saturating_sub(radius),
                    end.saturating_add(radius).min(text.len()),
                ));
            });
        }

        // Merges overlapping or touching windows of one adapter in place; the
        // merged prefix never reaches the window being read.
        windows.sort_unstable();
        let mut merged = 0;
        for i in 0..windows.len() {
            let (adapter_idx, start, end) = windows[i];
            if merged > 0 && windows[merged - 1].0 == adapter_idx && start <= windows[merged - 1].2
            {
                windows[merged - 1].2 = windows[merged - 1].2.max(end);
            } else {
                windows[merged] = (adapter_idx, start, end);
                merged += 1;
            }
        }
        windows.truncate(merged);
    }
}

/// Longest prefix the seed table indexes. Its table holds `4^len` entries of
/// two bytes; the cap keeps it inside the second-level cache.
pub(super) const MAX_SEED_PREFIX_LEN: usize = 10;

/// Exact-seed index over a table of every prefix code. The scan encodes the
/// text two bits per base in one rolling word and looks each prefix up, one
/// table read per base; a hit is confirmed against the whole seed before it
/// is reported, so a longer seed matches no more often than an exact search
/// for the whole seed would.
#[derive(Debug, Clone)]
pub(super) struct SeedTable {
    /// Prefix length: the shortest seed, capped at `MAX_SEED_PREFIX_LEN`,
    /// and shortened until the distinct prefixes fit the `u16` slots.
    pub(super) prefix: usize,
    /// Per prefix code, one plus the index into `lists`; zero for no seed.
    pub(super) slots: Vec<u16>,
    /// Per slot, the `(adapter, seed)` pairs whose seed begins with the prefix.
    pub(super) lists: Vec<Vec<(usize, Vec<u8>)>>,
}

impl SeedTable {
    /// Builds the table over `seeds`, each mapped to the adapters owning it.
    /// `None` when there are no seeds.
    pub(super) fn new(seeds: BTreeMap<Vec<u8>, Vec<usize>>) -> Option<Self> {
        let code_of = |seed: &[u8], prefix: usize| {
            seed[..prefix].iter().fold(0usize, |code, &b| {
                (code << 2) | usize::from(BASE_CODE[usize::from(b)])
            })
        };
        // A prefix of `SLOT_PREFIX_FLOOR` bases has at most 4^7 codes, which
        // always fit.
        const SLOT_PREFIX_FLOOR: usize = 7;
        let longest = seeds.keys().map(Vec::len).min()?.min(MAX_SEED_PREFIX_LEN);
        let prefix = (SLOT_PREFIX_FLOOR.min(longest)..=longest)
            .rev()
            .find(|&prefix| {
                let mut seen = vec![false; 1 << (2 * prefix)];
                let distinct = seeds
                    .keys()
                    .filter(|seed| !std::mem::replace(&mut seen[code_of(seed, prefix)], true))
                    .count();
                distinct <= usize::from(u16::MAX)
            })
            .unwrap_or(longest.min(SLOT_PREFIX_FLOOR));
        let mut slots = vec![0u16; 1 << (2 * prefix)];
        let mut lists: Vec<Vec<(usize, Vec<u8>)>> = Vec::new();
        for (seed, owners) in seeds {
            let code = code_of(&seed, prefix);
            if slots[code] == 0 {
                lists.push(Vec::new());
                slots[code] = u16::try_from(lists.len()).expect("The prefix bounds the slot count");
            }
            let list = &mut lists[usize::from(slots[code]) - 1];
            list.extend(owners.into_iter().map(|adapter| (adapter, seed.clone())));
        }
        Some(SeedTable {
            prefix,
            slots,
            lists,
        })
    }

    /// Calls `hit` with the owning adapter and the `[start, end)` span of every
    /// seed occurrence in `text`, which is normalized ACGT (see
    /// `normalize_into`). The loop body is a code lookup, a shift and a slot
    /// load per base; the slot slice is cut to exactly the code range so the
    /// index needs no bounds check.
    pub(super) fn scan(&self, text: &[u8], mut hit: impl FnMut(usize, (usize, usize))) {
        let prefix = self.prefix;
        if text.len() < prefix {
            return;
        }
        let mask = (1usize << (2 * prefix)) - 1;
        let slots = &self.slots[..=mask];
        let mut code = text[..prefix - 1].iter().fold(0usize, |code, &b| {
            (code << 2) | usize::from(BASE_CODE[usize::from(b)])
        });
        for (i, &b) in text[prefix - 1..].iter().enumerate() {
            code = ((code << 2) | usize::from(BASE_CODE[usize::from(b)])) & mask;
            let slot = slots[code];
            if slot != 0 {
                for (adapter, seed) in &self.lists[usize::from(slot) - 1] {
                    if text[i..].starts_with(seed) {
                        hit(*adapter, (i, i + seed.len()));
                    }
                }
            }
        }
    }
}

/// The two-bit code of a normalized base: A 0, C 1, G 2, T 3.
pub(super) const BASE_CODE: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
};

/// Returns the exact seeds of one strand of `pattern`: `max_edits + 1` pieces,
/// each expanded over its ambiguity codes. At most `max_edits` edits leave one
/// piece untouched, so a read carrying the pattern holds one expansion of one
/// piece verbatim. `None` when a piece expands past `MAX_SEED_EXPANSIONS` or
/// holds a byte outside the nucleotide alphabet.
pub(super) fn partition_seeds(pattern: &[u8], max_edits: usize) -> Option<Vec<Vec<u8>>> {
    let parts = (max_edits + 1).min(pattern.len());
    let mut seeds = Vec::new();
    for i in 0..parts {
        let start = i * pattern.len() / parts;
        let end = (i + 1) * pattern.len() / parts;
        seeds.extend(expand_iupac(&pattern[start..end])?);
    }
    Some(seeds)
}

/// Returns every plain ACGT string `piece` stands for, or `None` past
/// `MAX_SEED_EXPANSIONS` or on a byte outside the nucleotide alphabet.
pub(super) fn expand_iupac(piece: &[u8]) -> Option<Vec<Vec<u8>>> {
    let mut expansions: Vec<Vec<u8>> = vec![Vec::with_capacity(piece.len())];
    for &code in piece {
        let bases = iupac_bases(code)?;
        if expansions.len() * bases.len() > MAX_SEED_EXPANSIONS {
            return None;
        }
        let mut next = Vec::with_capacity(expansions.len() * bases.len());
        for prefix in &expansions {
            for &base in bases {
                let mut expansion = prefix.clone();
                expansion.push(base);
                next.push(expansion);
            }
        }
        expansions = next;
    }
    Some(expansions)
}

/// Length of the exact seeds that gate the overhang search. An exact partial
/// hit keeps at least `MIN_OVERLAP` intact bases and always passes the gate. A
/// partial hit with edits keeps an intact stretch of at least
/// `(o - e) / (e + 1)` bases for `o` aligned bases and `e` edits, which under
/// `partial_budget` falls below this length only for the shortest overlaps
/// with an edit near their middle; those are given up for the gate, which
/// spares the overhang search on most read ends without an adapter.
pub(super) const END_SEED_LEN: usize = 9;
