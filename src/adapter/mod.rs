//! Adapter, primer and barcode trimming.
//!
//! Searches each read window for catalog sequences with sassy, classifies every
//! accepted hit as a terminal trim or an interior excision, re-trims the ends
//! that an excision creates, and returns the kept spans. Presence detection,
//! de novo inference and the built-in catalog live in the submodules.

pub mod catalog;
pub mod detect;
pub mod infer;
pub mod preset;
pub mod resolve;
pub mod search;

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::OnceLock;

use search::{
    AmbiguousSearcher, EncodedAdapterBatch, Hit, MAX_TILED_PATTERN_LEN, PlainSearcher, Strands,
    encode_patterns, encoded_pattern_hits, for_each_hit, for_each_hit_in_texts, is_plain_acgt,
    iupac_bases, new_ambiguous_searcher, new_overhang_searcher, new_searcher,
};

thread_local! {
    /// The searchers and per-read buffers of this thread, reused across reads
    /// so `adapter_segments` allocates neither a searcher nor its scratch on
    /// every call. Per-thread state keeps the parallel workflows free of
    /// sharing.
    static STATE: RefCell<ThreadState> = RefCell::new(ThreadState::new());
}

/// One thread's searchers and per-read buffers. Each buffer keeps its
/// capacity across reads, so the adapter stage itself allocates only the
/// segments it returns; the remaining per-read allocations are sassy's own,
/// inside each search call.
struct ThreadState {
    /// The fast all-ACGT searcher. Used only when both the pattern and the
    /// searched text are plain ACGT; see `is_plain_acgt`.
    plain: PlainSearcher,
    /// The ambiguity-tolerant searcher, for a degenerate primer and for the
    /// tiled terminal batches, which sassy implements for the IUPAC profile
    /// only.
    ambiguous: AmbiguousSearcher,
    /// The overhang-aware searcher for the read ends, keyed by the overhang
    /// cost it was built with so a run at another error rate rebuilds it.
    overhang: Option<(f32, AmbiguousSearcher)>,
    /// The normalized read, when the input is not its own normalization.
    normalized: Vec<u8>,
    /// The normalized read reversed, for the reverse strand of every search.
    reversed: Vec<u8>,
    /// Candidate windows as `(adapter, start, end)`; see `candidate_windows`.
    windows: Vec<(usize, usize, usize)>,
    /// Masked end windows; see `search_residue`.
    mask: MaskScratch,
    /// Per-adapter end-seed flags for the head window; see `end_candidates`.
    head_flags: Vec<bool>,
    /// Per-adapter end-seed flags for the tail window.
    tail_flags: Vec<bool>,
}

impl ThreadState {
    /// Creates the searchers with empty buffers.
    fn new() -> Self {
        Self {
            plain: new_searcher(),
            ambiguous: new_ambiguous_searcher(),
            overhang: None,
            normalized: Vec::new(),
            reversed: Vec::new(),
            windows: Vec::new(),
            mask: MaskScratch::default(),
            head_flags: Vec::new(),
            tail_flags: Vec::new(),
        }
    }
}

/// What a catalog sequence is, which decides what a hit may do. Every role is
/// trimmed at the read ends; only an adapter splits a read at an interior hit,
/// since a primer or barcode inside a read is part of the molecule as often as
/// it is a chimera signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    /// A sequencing adapter: trimmed at the ends and excised in the interior.
    Adapter,
    /// A PCR or sequencing primer: trimmed at the ends only.
    Primer,
    /// A barcode or barcode flank: trimmed at the ends only.
    Barcode,
}

impl Role {
    /// Whether an interior hit of this role excises and splits the read.
    pub fn splits(self) -> bool {
        matches!(self, Role::Adapter)
    }

    /// Whether a partial hit hanging off a read end is accepted for this role.
    /// A barcode sits between its flanks, so a partial barcode at a read end is
    /// flank residue that the flank hit already trims.
    fn overhangs(self) -> bool {
        !matches!(self, Role::Barcode)
    }

    /// Display label.
    pub fn label(self) -> &'static str {
        match self {
            Role::Adapter => "adapter",
            Role::Primer => "primer",
            Role::Barcode => "barcode",
        }
    }
}

/// One searchable adapter, primer, barcode or flank sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    /// Display name, used in logs and the report.
    pub name: String,
    /// Nucleotide sequence; may contain IUPAC ambiguity codes.
    pub seq: Vec<u8>,
    /// What the sequence is, which decides whether an interior hit splits.
    pub role: Role,
}

/// Resolved adapter-trimming settings for a run.
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    /// Sequences searched for, in configuration order.
    pub adapters: Vec<Adapter>,
    /// End-match tolerance as a fraction of adapter length (`k_end`), also the
    /// per-base cost of a pattern base hanging off a read end.
    pub error_rate: f64,
    /// Bases at each end within which a hit is terminal (trim) rather than
    /// interior (split).
    pub end_size: usize,
    /// Whether interior adapters split the read; `false` is ends-only
    /// (`--adapter-ends-only`).
    pub split: bool,
    /// Shortest segment worth keeping between two excisions. Excisions closer
    /// than this merge into one, since the bases between them would be
    /// discarded by the length filter anyway. Follows `--min-length`.
    pub min_piece: usize,
    /// Exact-seed index for lossless whole-read candidate filtering, built
    /// lazily once presence detection or inference has finalized `adapters`.
    pub(crate) candidate_index: OnceLock<CandidateIndex>,
}

impl AdapterConfig {
    /// Replaces the adapter set and discards the candidate index built for the
    /// previous set.
    pub(crate) fn replace_adapters(&mut self, adapters: Vec<Adapter>) {
        self.adapters = adapters;
        self.candidate_index = OnceLock::new();
    }

    /// Returns whether a barcode sequence matches the read at `[start, end)`
    /// within the error budget: a barcode-role entry of the configured set,
    /// or a catalog barcode with the number of the barcode call `call`. The
    /// span is widened by `BARCODE_SPAN_SLACK` bases on each side.
    pub(crate) fn barcode_span_verified(
        &self,
        seq: &[u8],
        start: usize,
        end: usize,
        call: Option<&[u8]>,
    ) -> bool {
        let lo = start.saturating_sub(BARCODE_SPAN_SLACK);
        let hi = end.saturating_add(BARCODE_SPAN_SLACK).min(seq.len());
        if hi <= lo {
            return false;
        }
        let text = seq[lo..hi].to_ascii_uppercase();
        let number = call.and_then(barcode_number);
        let named = catalog_barcodes().iter().filter(|entry| {
            number.is_some_and(|n| barcode_number(entry.name.as_bytes()) == Some(n))
        });
        let mut searcher = search::new_ambiguous_searcher();
        self.adapters
            .iter()
            .filter(|entry| entry.role == Role::Barcode)
            .chain(named)
            .filter(|entry| entry.seq.len() >= MIN_PATTERN_LEN)
            .any(|entry| {
                let pattern = entry.seq.to_ascii_uppercase();
                let k = edit_budget(self.error_rate, pattern.len());
                !search::hits(&mut searcher, &pattern, &text, k).is_empty()
            })
    }
}

/// Bases by which a recorded barcode span is widened before verification.
const BARCODE_SPAN_SLACK: usize = 3;

/// Every barcode-role catalog entry, built once.
fn catalog_barcodes() -> &'static [Adapter] {
    static ENTRIES: OnceLock<Vec<Adapter>> = OnceLock::new();
    ENTRIES.get_or_init(|| {
        preset::preset(preset::Kit::ALL)
            .into_iter()
            .filter(|entry| entry.role == Role::Barcode)
            .collect()
    })
}

/// Returns the numeric value of the trailing digits of a barcode name, such
/// as `07` in `SQK-NBD114-24_barcode07` or `BC07`.
fn barcode_number(name: &[u8]) -> Option<u32> {
    let digits = name.iter().rev().take_while(|b| b.is_ascii_digit()).count();
    (digits > 0)
        .then(|| {
            std::str::from_utf8(&name[name.len() - digits..])
                .ok()?
                .parse()
                .ok()
        })
        .flatten()
}

/// Returns the edit budget for a `len`-base pattern at `rate`, rounded down.
/// The epsilon keeps an integral product whose double lands below its integer.
pub(crate) fn edit_budget(rate: f64, len: usize) -> usize {
    (rate * len as f64 + 1e-9).floor() as usize
}

/// Expected chance interior matches per read, over both strands, that an
/// interior edit budget may admit under an independent uniform base model.
const INTERIOR_CHANCE_HITS_PER_READ: f64 = 1e-4;

/// Expected chance terminal matches per read, over both strands and both end
/// zones, that terminal edit budgets may admit under the same null model. The
/// model sums alignment paths and overstates the chance rate, so the bound is
/// conservative. It applies to each pattern alone, which at the default error
/// rate lowers the budget of 11-, 12- and 15-base patterns by one edit, and to
/// the distinct sequences of the set together; see `family_budgets`.
const TERMINAL_CHANCE_HITS_PER_READ: f64 = 0.1;

/// Read-length class `c` holds reads shorter than `2^(INTERIOR_CLASS_BITS + c)`
/// bases. Class 0 also holds every shorter read.
const INTERIOR_CLASS_BITS: u32 = 12;

/// Number of read-length classes. The last class holds every longer read.
const INTERIOR_CLASSES: usize = 20;

/// Returns the read-length class of a read of `read_len` bases.
fn interior_class(read_len: usize) -> usize {
    let bits = usize::BITS - read_len.leading_zeros();
    (bits.saturating_sub(INTERIOR_CLASS_BITS) as usize).min(INTERIOR_CLASSES - 1)
}

/// Interior alignment start positions, over both strands, in a read at the
/// ceiling of read-length class `class`.
fn interior_positions(class: usize) -> f64 {
    2.0 * 2f64.powi(INTERIOR_CLASS_BITS as i32 + class as i32)
}

/// Returns, for each edit count up to `max_edits`, the probability under the
/// independent uniform DNA null model that `pattern` matches at one position
/// within that many edits. The recurrence sums alignment-path probabilities,
/// including substitutions, insertions and deletions, and therefore
/// overcounts sequences admitting multiple alignments. IUPAC ambiguity
/// increases the probability of a zero-cost match.
fn chance_cumulative(pattern: &[u8], max_edits: usize) -> Vec<f64> {
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
fn family_budgets(
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
struct Budget {
    /// Pattern length in bases.
    len: usize,
    /// Edit budget of the terminal search, which a hit anchored at the read
    /// end or at an accepted hit may use; see `Keep::settle`.
    k_end: usize,
    /// Edit budget for a terminal hit anywhere in the end zone. At most
    /// `k_end`.
    k_far: usize,
    /// Edit budget for interior hits, per read-length class. Budgets do not
    /// increase with the class.
    k_mid: [usize; INTERIOR_CLASSES],
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
    fn new(pattern: &[u8], error_rate: f64, end_size: usize) -> Self {
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
    fn interior(&self, read_len: usize) -> usize {
        self.k_mid[interior_class(read_len)]
    }

    /// Returns the largest interior edit budget over all read lengths.
    fn interior_max(&self) -> usize {
        self.k_mid[0]
    }
}

/// Upper bound on the plain strings one seed piece may expand to. A piece past
/// it marks its adapter `unfiltered`, and the interior search covers the whole
/// read for that adapter instead of candidate windows.
const MAX_SEED_EXPANSIONS: usize = 256;

/// Exact-seed index over the adapter set. Partition seeds, looked up through a
/// prefix table, bound the interior search to candidate windows, and
/// equal-length barcode entries are grouped into SIMD batches for the terminal
/// search.
#[derive(Debug, Clone)]
pub(crate) struct CandidateIndex {
    /// Table over the interior seeds of every splitting adapter; `None` when
    /// no adapter has seeds.
    seeds: Option<SeedTable>,
    /// Per-adapter edit budgets, computed once per adapter set.
    budgets: Vec<Budget>,
    /// Per-adapter `is_plain_acgt`, which selects the search profile.
    plain: Vec<bool>,
    /// Adapters with no usable seeds (see `MAX_SEED_EXPANSIONS`).
    unfiltered: Vec<bool>,
    /// Equal-length adapter groups searched together over the end windows.
    terminal_batches: Vec<TerminalBatch>,
    /// Adapters searched one pattern at a time over the end windows: those of
    /// searchable length that no batch covers.
    singletons: Vec<bool>,
    /// Table over every `END_SEED_LEN`-mer of the entries eligible for
    /// partial matching, both strands; `None` when there are none. Gates the
    /// overhang search of an end to the entries with an exact seed in it.
    end_seeds: Option<SeedTable>,
    /// The longest end-window reach over the partial-matching entries.
    end_reach: usize,
}

/// Equal-length adapters searched together through sassy's pattern-parallel
/// API, which matches the forward and reverse-complement strands for the batch.
#[derive(Debug, Clone)]
struct TerminalBatch {
    /// Indices into `AdapterConfig::adapters`, in pattern order.
    adapter_indices: Vec<usize>,
    /// The adapter sequences encoded once for the tiled search.
    encoded: EncodedAdapterBatch,
    /// The shared pattern length.
    len: usize,
    /// The largest terminal edit budget of the batch; a hit above the budget
    /// of its own entry is discarded.
    k_end: usize,
}

impl CandidateIndex {
    /// Builds the index for `adapters` searched over end zones of `end_size`
    /// bases; interior seeds are built only when `include_interior`, and
    /// only for roles that split.
    fn new(adapters: &[Adapter], error_rate: f64, end_size: usize, include_interior: bool) -> Self {
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
    fn end_candidates(
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
    fn candidate_windows(&self, text: &[u8], windows: &mut Vec<(usize, usize, usize)>) {
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
const MAX_SEED_PREFIX_LEN: usize = 10;

/// Exact-seed index over a table of every prefix code. The scan encodes the
/// text two bits per base in one rolling word and looks each prefix up, one
/// table read per base; a hit is confirmed against the whole seed before it
/// is reported, so a longer seed matches no more often than an exact search
/// for the whole seed would.
#[derive(Debug, Clone)]
struct SeedTable {
    /// Prefix length: the shortest seed, capped at `MAX_SEED_PREFIX_LEN`,
    /// and shortened until the distinct prefixes fit the `u16` slots.
    prefix: usize,
    /// Per prefix code, one plus the index into `lists`; zero for no seed.
    slots: Vec<u16>,
    /// Per slot, the `(adapter, seed)` pairs whose seed begins with the prefix.
    lists: Vec<Vec<(usize, Vec<u8>)>>,
}

impl SeedTable {
    /// Builds the table over `seeds`, each mapped to the adapters owning it.
    /// `None` when there are no seeds.
    fn new(seeds: BTreeMap<Vec<u8>, Vec<usize>>) -> Option<Self> {
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
    fn scan(&self, text: &[u8], mut hit: impl FnMut(usize, (usize, usize))) {
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
const BASE_CODE: [u8; 256] = {
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
fn partition_seeds(pattern: &[u8], max_edits: usize) -> Option<Vec<Vec<u8>>> {
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
fn expand_iupac(piece: &[u8]) -> Option<Vec<Vec<u8>>> {
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

/// Returns the IUPAC complement of one base or ambiguity code, case preserved.
/// `S`, `W` and `N` are their own complements; any other byte passes through.
fn complement(base: u8) -> u8 {
    let upper = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    };
    if base.is_ascii_lowercase() {
        upper.to_ascii_lowercase()
    } else {
        upper
    }
}

/// Returns the reverse complement of `seq`, code by code.
pub(crate) fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Minimum searchable pattern length: a shorter pattern matches almost anywhere
/// under any error budget and is never searched standalone. Catalog flanks
/// below it are omitted from the catalog for the same reason.
pub const MIN_PATTERN_LEN: usize = 11;

/// Minimum aligned pattern length of a partial hit at a read end. A shorter
/// overlap matches a random read end too often to be evidence of an adapter.
pub const MIN_OVERLAP: usize = 10;

/// Length of the exact seeds that gate the overhang search. An exact partial
/// hit keeps at least `MIN_OVERLAP` intact bases and always passes the gate. A
/// partial hit with edits keeps an intact stretch of at least
/// `(o - e) / (e + 1)` bases for `o` aligned bases and `e` edits, which under
/// `partial_budget` falls below this length only for the shortest overlaps
/// with an edit near their middle; those are given up for the gate, which
/// spares the overhang search on most read ends without an adapter.
const END_SEED_LEN: usize = 9;

/// Returns the overhang cost the overhang searcher charges a hit at per-base
/// rate `alpha`: `floor(alpha * bases)` for each overhanging side, at the
/// searcher's `f32` precision.
fn overhang_cost(alpha: f32, left: usize, right: usize) -> usize {
    let side = |bases: usize| (bases as f32 * alpha).floor() as usize;
    side(left) + side(right)
}

/// Returns the edit budget of a partial hit whose `overlap` bases aligned
/// inside the read. The first `MIN_OVERLAP` bases must match exactly and the
/// rate applies to the remainder, so the shortest accepted overlaps carry no
/// tolerance and a random read end is not mistaken for adapter residue.
fn partial_budget(rate: f64, overlap: usize) -> usize {
    edit_budget(rate, overlap.saturating_sub(MIN_OVERLAP - 1))
}

/// Outboard flank length at or below which a hit covered by both end zones is
/// a terminal trim rather than an excision, and at or below which the bases
/// between two excisions are merged into one. A flank this short is adapter
/// residue or junk and is not worth a read of its own.
const FLANK_SLACK: usize = MIN_PATTERN_LEN;

/// Terminal classification of a hit: which end, if any, it trims.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Terminal {
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

/// Searches the adapter at `adapter_idx` in `text` and passes each hit to
/// `accept`, choosing the profile by the pattern's alphabet (`plain`,
/// precomputed per adapter).
///
/// A degenerate primer needs the IUPAC profile so its wobble positions match the
/// bases they stand for. A plain ACGT pattern takes the faster DNA profile,
/// which matters because a narrowed adapter set is searched one pattern at a
/// time rather than batched across SIMD lanes.
///
/// Sassy's per-pattern search rebuilds the pattern profile and the
/// complemented pattern on every call, a few allocations each. The tiled
/// searcher would avoid them, but its column-major scan has no early exit and
/// fills one lane of eight with a single pattern, which costs more on these
/// short windows than the allocations do.
///
/// The DNA profile requires a plain pattern and a plain read.
fn search(
    engine: &mut Engine<'_>,
    index: &CandidateIndex,
    adapter_idx: usize,
    pattern: &[u8],
    text: Strands<'_>,
    k: usize,
    accept: impl FnMut(Hit),
) {
    if engine.plain_read && index.plain[adapter_idx] {
        for_each_hit(engine.plain, pattern, &text, k, accept);
    } else {
        for_each_hit(engine.ambiguous, pattern, &text, k, accept);
    }
}

/// The IUPAC profile's nonmatching symbol for uncalled read bases.
const AMBIGUOUS_READ_BASE: u8 = b'X';

/// Returns the read as every searcher sees it: uppercase, with each byte
/// outside ACGT rewritten to `AMBIGUOUS_READ_BASE`. An uppercase plain read,
/// the common case, is returned as is; any other read is rewritten into
/// `buf`, which keeps its capacity across calls. Sassy's profiles fold case
/// themselves; the seed table does not, so the text is folded once here
/// rather than on every lookup. Also returns whether the normalized read
/// holds no `AMBIGUOUS_READ_BASE`.
pub(crate) fn normalize_into<'a>(window: &'a [u8], buf: &'a mut Vec<u8>) -> (&'a [u8], bool) {
    if is_upper_acgt(window) {
        return (window, true);
    }
    buf.clear();
    let mut plain = true;
    buf.extend(window.iter().map(|&b| {
        let normalized = normalize_base(b);
        plain &= normalized != AMBIGUOUS_READ_BASE;
        normalized
    }));
    (buf, plain)
}

/// Returns the normalized form of one read byte: its uppercase base, or
/// `AMBIGUOUS_READ_BASE` for a byte outside ACGT.
#[inline]
fn normalize_base(b: u8) -> u8 {
    match b {
        b'A' | b'C' | b'G' | b'T' => b,
        b'a' | b'c' | b'g' | b't' => b.to_ascii_uppercase(),
        _ => AMBIGUOUS_READ_BASE,
    }
}

/// Returns whether every byte is an uppercase A/C/G/T, so that the window is
/// its own normalization. Folded 32 bytes at a time without an early exit
/// and with the four comparisons or-ed rather than matched, which lets the
/// scan vectorize; it runs over every base of every read.
fn is_upper_acgt(seq: &[u8]) -> bool {
    let upper_acgt = |b: u8| (b == b'A') | (b == b'C') | (b == b'G') | (b == b'T');
    let mut chunks = seq.chunks_exact(32);
    let body = chunks.all(|chunk| chunk.iter().fold(true, |ok, &b| ok & upper_acgt(b)));
    body && chunks.remainder().iter().all(|&b| upper_acgt(b))
}

/// Emits a trace event for one adapter hit: the sequence, its span, its edit
/// cost, and the action taken (a terminal trim, an excision, or none).
fn trace_hit(name: &str, start: usize, end: usize, cost: usize, action: Option<HitAction>) {
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
enum HitAction {
    /// Terminal hit at the 5' end: the keep-boundary moved inward past it.
    TrimFivePrime,
    /// Terminal hit at the 3' end.
    TrimThreePrime,
    /// Interior hit: the span is cut out and both flanks are kept.
    Excise,
}

/// Returns the end nearer to a hit at `[start, end)` in a length-`n` window.
fn nearer_end(start: usize, end: usize, n: usize) -> Terminal {
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
fn classify_terminal(start: usize, end: usize, n: usize, end_size: usize) -> Terminal {
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
fn ends_only_terminal(start: usize, end: usize, n: usize, end_size: usize) -> Terminal {
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
enum Site {
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
fn terminal_windows(n: usize, end_size: usize, len: usize, k_end: usize) -> (usize, usize) {
    let reach = end_size + len + k_end;
    (reach.min(n), n.saturating_sub(reach))
}

/// Accumulator for the accepted hits of one window: the keep boundaries and
/// interior cuts. Coordinates are those of the window, `[0, n)`.
struct Keep<'a> {
    /// The configured adapters, for roles and names.
    adapters: &'a [Adapter],
    /// Per-adapter edit budgets, for the anchoring of terminal hits.
    budgets: &'a [Budget],
    /// The run's error rate, which scales the budget of a partial hit.
    error_rate: f64,
    /// Window length.
    n: usize,
    /// Terminal zone depth, capped at `n`.
    end_size: usize,
    /// Whether interior hits split the read.
    split: bool,
    /// 5' keep boundary; advances inward on 5' trims.
    lo: usize,
    /// 3' keep boundary; retreats inward on 3' trims.
    hi: usize,
    /// Accepted excisions, merged by `into_cuts`.
    interior: Vec<(usize, usize)>,
    /// The adapters whose hits trimmed or excised, for presence detection.
    acted: Vec<usize>,
    /// Terminal trims above the `k_far` budget of their adapter, applied by
    /// `settle` once anchored.
    deferred: Vec<Deferred>,
}

/// A terminal trim held until it is anchored at the read end or at an
/// accepted hit.
#[derive(Debug, Clone, Copy)]
struct Deferred {
    /// Index into the configured adapters.
    adapter_idx: usize,
    /// Hit start in window coordinates.
    start: usize,
    /// Hit end in window coordinates.
    end: usize,
    /// Edit cost of the hit.
    cost: usize,
    /// `TrimFivePrime` or `TrimThreePrime`.
    action: HitAction,
}

impl<'a> Keep<'a> {
    /// Creates an accumulator that keeps the whole `[0, n)` window. `split`
    /// selects the classification: with it, a hit covered by both end zones
    /// may excise; without it, every hit trims an end.
    fn new(cfg: &'a AdapterConfig, index: &'a CandidateIndex, n: usize, split: bool) -> Self {
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
            acted: Vec::new(),
            deferred: Vec::new(),
        }
    }

    /// Returns whether a partial hit is acceptable: enough of the pattern
    /// aligned, flush with the window end it hangs off, and within the edit
    /// budget of the aligned part alone. A whole-pattern hit always passes.
    fn partial_hit_is_valid(&self, adapter_idx: usize, hit: Hit) -> bool {
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
    fn residue_within_budget(&self, hit: Hit, overlap: usize) -> bool {
        let charged = overhang_cost(
            self.error_rate as f32,
            hit.left_overhang,
            hit.right_overhang,
        );
        hit.cost.saturating_sub(charged) <= partial_budget(self.error_rate, overlap)
    }

    /// Classifies one hit and applies it when `site` owns the outcome.
    fn accept(&mut self, site: Site, adapter_idx: usize, hit: Hit) {
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
                action,
            });
            return;
        }
        self.apply(adapter_idx, start, end, cost, action);
    }

    /// Applies an accepted hit to the keep boundaries or the excisions.
    fn apply(
        &mut self,
        adapter_idx: usize,
        start: usize,
        end: usize,
        cost: usize,
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
        match action {
            HitAction::TrimFivePrime => self.lo = self.lo.max(end),
            HitAction::TrimThreePrime => self.hi = self.hi.min(start),
            HitAction::Excise => self.interior.push((start, end)),
        }
    }

    /// Applies each deferred trim that is anchored: a 5' hit starting within
    /// `FLANK_SLACK` of the 5' keep boundary, or a 3' hit ending within
    /// `FLANK_SLACK` of the 3' boundary. An applied trim moves the boundary,
    /// which may anchor further deferred trims behind it.
    fn settle(&mut self) {
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
            self.apply(d.adapter_idx, d.start, d.end, d.cost, d.action);
        }
    }

    /// Returns the keep boundaries and the excisions clipped to them, merged:
    /// two excisions overlapping, touching, or separated by at most
    /// `FLANK_SLACK` bases or fewer than `min_piece` bases become one, since
    /// the bases between them are junction residue or a piece the length
    /// filter would discard.
    fn into_cuts(mut self, min_piece: usize) -> (usize, usize, Vec<(usize, usize)>) {
        self.settle();
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

/// A normalized read and its reversal, borrowed from per-thread buffers so
/// every two-strand search copies nothing.
#[derive(Debug, Clone, Copy)]
struct Read<'a> {
    /// The read as every searcher sees it; see `normalize_into`.
    window: &'a [u8],
    /// `window` reversed.
    reversed: &'a [u8],
}

impl<'a> Read<'a> {
    /// Returns the strands of `window[start..end]`.
    fn strands(&self, start: usize, end: usize) -> Strands<'a> {
        let n = self.window.len();
        Strands {
            forward: &self.window[start..end],
            reversed: &self.reversed[n - end..n - start],
        }
    }
}

/// One read's search context: the configuration, its index, and the read.
#[derive(Clone, Copy)]
struct Context<'a> {
    /// The run's adapter settings.
    cfg: &'a AdapterConfig,
    /// The index built for `cfg.adapters`.
    index: &'a CandidateIndex,
    /// The read under search.
    read: Read<'a>,
}

/// The per-thread searchers and buffers one read is processed with, borrowed
/// from the thread's `ThreadState` for the duration of `adapter_segments`.
struct Engine<'a> {
    /// Whether the complete read contains only ACGT bases.
    plain_read: bool,
    /// The DNA-profile searcher, for a plain pattern.
    plain: &'a mut PlainSearcher,
    /// The IUPAC-profile searcher, for a degenerate pattern and for the
    /// terminal batches.
    ambiguous: &'a mut AmbiguousSearcher,
    /// The IUPAC-profile searcher with overhang alignment, for the read ends.
    overhang: &'a mut AmbiguousSearcher,
    /// Candidate windows of the interior search; see `candidate_windows`.
    windows: &'a mut Vec<(usize, usize, usize)>,
    /// Masked end windows of the residue search.
    mask: &'a mut MaskScratch,
    /// Per-adapter end-seed flags for the head window; see `end_candidates`.
    head_flags: &'a mut Vec<bool>,
    /// Per-adapter end-seed flags for the tail window.
    tail_flags: &'a mut Vec<bool>,
}

/// A span `[start, end)` of the read that is searched as a read of its own.
type Span = (usize, usize);

/// Searches every equal-length batch over the two end windows of the span.
/// All adapters in a batch share a length and budget, so the windows are
/// shared too; this collapses a kit's equal-length barcode searches into one
/// SIMD pattern search per end. Hits are passed to `keep` in span
/// coordinates.
fn search_batched(
    ctx: Context<'_>,
    span: Span,
    searcher: &mut AmbiguousSearcher,
    keep: &mut Keep<'_>,
) {
    let (ws, we) = span;
    let n = we - ws;
    for batch in &ctx.index.terminal_batches {
        let (head_end, tail_start) = terminal_windows(n, keep.end_size, batch.len, batch.k_end);
        let head = ctx.read.strands(ws, ws + head_end);
        let tail = ctx.read.strands(ws + tail_start, we);
        accept_batch_hits(batch, searcher, head, 0, Site::Head, keep);
        accept_batch_hits(
            batch,
            searcher,
            tail,
            tail_start,
            Site::Tail { head_end },
            keep,
        );
    }
}

/// Searches one batch over `text`, a window starting at `offset` in the span,
/// and passes every hit to `keep` at `site` in span coordinates.
fn accept_batch_hits(
    batch: &TerminalBatch,
    searcher: &mut AmbiguousSearcher,
    text: Strands<'_>,
    offset: usize,
    site: Site,
    keep: &mut Keep<'_>,
) {
    let accept = |pattern_idx: usize, start: usize, end: usize, cost: usize| {
        let adapter_idx = batch.adapter_indices[pattern_idx];
        if cost > keep.budgets[adapter_idx].k_end {
            return;
        }
        keep.accept(
            site,
            adapter_idx,
            Hit {
                start: offset + start,
                end: offset + end,
                cost,
                left_overhang: 0,
                right_overhang: 0,
            },
        );
    };
    encoded_pattern_hits(
        searcher,
        &batch.encoded,
        text.forward,
        text.reversed,
        batch.k_end,
        accept,
    );
}

/// Searches every adapter without an equal-length partner over the two end
/// windows of the span, one pattern at a time, for whole-pattern hits.
fn search_singletons(ctx: Context<'_>, span: Span, engine: &mut Engine<'_>, keep: &mut Keep<'_>) {
    let (ws, we) = span;
    let n = we - ws;
    for (adapter_idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
        if !ctx.index.singletons[adapter_idx] {
            continue;
        }
        let Budget { len, k_end, .. } = ctx.index.budgets[adapter_idx];
        let (head_end, tail_start) = terminal_windows(n, keep.end_size, len, k_end);
        let windows = [
            ctx.read.strands(ws, ws + head_end),
            ctx.read.strands(ws + tail_start, we),
        ];
        let accept = |text_idx: usize, h: Hit| {
            if text_idx == 0 {
                keep.accept(Site::Head, adapter_idx, h);
            } else {
                keep.accept(Site::Tail { head_end }, adapter_idx, shifted(h, tail_start));
            }
        };
        if engine.plain_read && ctx.index.plain[adapter_idx] {
            for_each_hit_in_texts(engine.plain, &adapter.seq, &windows, k_end, accept);
        } else {
            for_each_hit_in_texts(engine.ambiguous, &adapter.seq, &windows, k_end, accept);
        }
    }
}

/// Searches the partial-matching entries with overhang alignment over each
/// end window of the span that the whole-pattern pass left untrimmed, for
/// the entries whose end seeds occur in that window. The end-seed flags are
/// set only for the partial-matching entries (see `CandidateIndex::new`), so
/// they gate the role as well as the seed. A partial hit flush with the read
/// end trims it; see `Keep::partial_hit_is_valid`.
fn search_partial(ctx: Context<'_>, span: Span, engine: &mut Engine<'_>, keep: &mut Keep<'_>) {
    let (ws, we) = span;
    let n = we - ws;
    let head_open = keep.lo == 0;
    let tail_open = keep.hi == n;
    for (adapter_idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
        let Budget { len, k_end, .. } = ctx.index.budgets[adapter_idx];
        let (head_end, tail_start) = terminal_windows(n, keep.end_size, len, k_end);
        if head_open && engine.head_flags[adapter_idx] {
            let head = ctx.read.strands(ws, ws + head_end);
            for_each_hit(engine.overhang, &adapter.seq, &head, k_end, |h| {
                if h.left_overhang > 0 {
                    keep.accept(Site::Head, adapter_idx, h);
                }
            });
        }
        if tail_open && engine.tail_flags[adapter_idx] {
            let tail = ctx.read.strands(ws + tail_start, we);
            for_each_hit(engine.overhang, &adapter.seq, &tail, k_end, |h| {
                if h.right_overhang > 0 {
                    keep.accept(Site::Tail { head_end }, adapter_idx, shifted(h, tail_start));
                }
            });
        }
    }
}

/// Mask lengths of the residue search: the outboard bases of an end window
/// rewritten to `N` so a partial adapter followed by that many unalignable
/// bases aligns as if flush with the read end. A remnant of `r` bases followed
/// by `j` bases is found by a mask `m` when `m - r + MIN_OVERLAP <= j <= m`;
/// the steps keep every junk length up to the last mask covered for remnants
/// of about half the adapter length or more.
const RESIDUE_MASKS: &[usize] = &[12, 24, 36];

/// Masked copies of one end window for the residue search, reused across reads.
#[derive(Debug, Default)]
struct MaskScratch {
    /// The window with its outboard bases rewritten to `N`.
    text: Vec<u8>,
    /// `text` reversed.
    reversed: Vec<u8>,
}

impl MaskScratch {
    /// Fills the buffers with `window`, masking its first `mask` bases when
    /// `head` and its last `mask` bases otherwise, and returns the strands.
    fn fill(&mut self, window: &[u8], mask: usize, head: bool) -> Strands<'_> {
        self.text.clear();
        self.text.extend_from_slice(window);
        let n = self.text.len();
        let masked = if head { 0..mask } else { n - mask..n };
        self.text[masked].fill(b'N');
        self.reversed.clear();
        self.reversed.extend(self.text.iter().rev());
        Strands {
            forward: &self.text,
            reversed: &self.reversed,
        }
    }
}

/// Searches the adapter entries once more at each end of the span that the
/// plain pass left untrimmed, over the end window with its outboard bases
/// masked to `N` at each of `RESIDUE_MASKS`. The IUPAC profile matches `N`
/// in the text at no cost, so a partial adapter followed by unalignable bases
/// aligns through the mask as if it were flush with the read end. A hit must
/// align at least `MIN_OVERLAP` unmasked bases within the partial budget of
/// that overlap and lie inside the current keep boundaries, so a hit the
/// plain pass already trimmed at the opposite end is not read as residue of
/// this one; the masked stretch itself is trimmed with the hit. The end-seed
/// flags gate the entries as in `search_partial`.
fn search_residue(ctx: Context<'_>, span: Span, engine: &mut Engine<'_>, keep: &mut Keep<'_>) {
    let (ws, we) = span;
    let n = we - ws;
    let window = ctx.read.window;
    let retry_head = keep.lo == 0;
    let retry_tail = keep.hi == n;
    for (adapter_idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
        let retry_head = retry_head && engine.head_flags[adapter_idx];
        let retry_tail = retry_tail && engine.tail_flags[adapter_idx];
        if !adapter.role.splits() || (!retry_head && !retry_tail) {
            continue;
        }
        let Budget { len, k_end, .. } = ctx.index.budgets[adapter_idx];
        let reach = (keep.end_size + len + k_end).min(n);
        for &masked in RESIDUE_MASKS {
            if masked + MIN_OVERLAP > reach {
                break;
            }
            if retry_head {
                let text = engine.mask.fill(&window[ws..ws + reach], masked, true);
                for_each_hit(engine.overhang, &adapter.seq, &text, k_end, |h| {
                    let overlap = h.end.saturating_sub(masked.max(h.start));
                    if h.start <= masked
                        && h.end <= keep.hi
                        && h.right_overhang == 0
                        && overlap >= MIN_OVERLAP
                        && keep.residue_within_budget(h, overlap)
                    {
                        trace_hit(
                            &adapter.name,
                            h.start,
                            h.end,
                            h.cost,
                            Some(HitAction::TrimFivePrime),
                        );
                        keep.lo = keep.lo.max(h.end);
                    }
                });
            }
            if retry_tail {
                let start = n - reach;
                let text = engine.mask.fill(&window[ws + start..we], masked, false);
                let unmasked = reach - masked;
                for_each_hit(engine.overhang, &adapter.seq, &text, k_end, |h| {
                    let overlap = unmasked.min(h.end).saturating_sub(h.start);
                    if h.end >= unmasked
                        && start + h.start >= keep.lo
                        && h.left_overhang == 0
                        && overlap >= MIN_OVERLAP
                        && keep.residue_within_budget(h, overlap)
                    {
                        let (s, e) = (start + h.start, start + h.end.min(reach));
                        trace_hit(&adapter.name, s, e, h.cost, Some(HitAction::TrimThreePrime));
                        keep.hi = keep.hi.min(s);
                    }
                });
            }
        }
    }
}

/// Returns `hit` moved `offset` bases to the right.
fn shifted(hit: Hit, offset: usize) -> Hit {
    Hit {
        start: hit.start + offset,
        end: hit.end + offset,
        ..hit
    }
}

/// Searches every adapter's candidate windows at the interior budget of the
/// read's length class. Exact partition seeds cover the largest interior
/// budget, so they identify every possible interior match, and the search
/// runs at the class budget rather than the looser end budget. An adapter
/// below `MIN_PATTERN_LEN` or without a splitting role has no seeds and no
/// windows.
fn search_interior(ctx: Context<'_>, engine: &mut Engine<'_>, keep: &mut Keep<'_>) {
    ctx.index.candidate_windows(ctx.read.window, engine.windows);
    for i in 0..engine.windows.len() {
        let (adapter_idx, start, end) = engine.windows[i];
        let k_mid = ctx.index.budgets[adapter_idx].interior(ctx.read.window.len());
        search(
            engine,
            ctx.index,
            adapter_idx,
            &ctx.cfg.adapters[adapter_idx].seq,
            ctx.read.strands(start, end),
            k_mid,
            |h| keep.accept(Site::Interior, adapter_idx, shifted(h, start)),
        );
    }
}

/// Runs the terminal passes over the span: the batched and singleton
/// whole-pattern searches, then the partial and residue searches over the
/// ends they left untrimmed.
fn search_terminal(ctx: Context<'_>, span: Span, engine: &mut Engine<'_>, keep: &mut Keep<'_>) {
    if !ctx.index.terminal_batches.is_empty() {
        search_batched(ctx, span, engine.ambiguous, keep);
    }
    search_singletons(ctx, span, engine, keep);
    keep.settle();
    let (ws, we) = span;
    if keep.lo != 0 && keep.hi != we - ws {
        return;
    }
    ctx.index.end_candidates(
        ctx.read.window,
        ws,
        we,
        keep.end_size,
        engine.head_flags,
        engine.tail_flags,
    );
    search_partial(ctx, span, engine, keep);
    search_residue(ctx, span, engine, keep);
    keep.settle();
}

/// Computes the adapter keep segments for `window`: terminal hits within
/// `end_size` of an end trim that end inward, and interior hits (at the
/// stricter `k_mid`) excise and split. Every segment an excision creates is
/// then searched again at its new ends, so the residue next to a junction
/// (a truncated adapter, a primer, a barcode, junk) is trimmed as it would be
/// at a physical read end.
///
/// Under `--adapter-ends-only` (`cfg.split` false) only the two end zones are
/// searched, since no interior hit could be acted on.
///
/// Returns `[start, end)` spans in `window` coordinates.
pub fn adapter_segments(window: &[u8], cfg: &AdapterConfig) -> Vec<(usize, usize)> {
    segments_tallied(window, cfg, None)
}

/// `adapter_segments` that also marks in `acted` every adapter whose hit
/// trimmed or excised part of `window`, for presence detection.
pub(crate) fn adapter_segments_tallied(
    window: &[u8],
    cfg: &AdapterConfig,
    acted: &mut [bool],
) -> Vec<(usize, usize)> {
    segments_tallied(window, cfg, Some(acted))
}

/// Runs the search passes over per-thread state; see `adapter_segments`.
fn segments_tallied(
    window: &[u8],
    cfg: &AdapterConfig,
    acted: Option<&mut [bool]>,
) -> Vec<(usize, usize)> {
    let n = window.len();
    if n == 0 {
        return vec![];
    }
    if cfg.adapters.is_empty() {
        return vec![(0, n)];
    }
    let index = cfg.candidate_index.get_or_init(|| {
        CandidateIndex::new(&cfg.adapters, cfg.error_rate, cfg.end_size, cfg.split)
    });
    // The overhang cost per base of the terminal search is the error rate, so
    // a partial adapter costs what its missing part would have been allowed
    // in edits.
    let alpha = cfg.error_rate as f32;
    STATE.with_borrow_mut(|state| {
        let ThreadState {
            plain,
            ambiguous,
            overhang,
            normalized,
            reversed,
            windows,
            mask,
            head_flags,
            tail_flags,
        } = state;
        let (window, plain_read) = normalize_into(window, normalized);
        reversed.clear();
        reversed.extend_from_slice(window);
        reversed.reverse();
        let overhang = match overhang {
            Some((a, s)) if *a == alpha => s,
            slot => &mut slot.insert((alpha, new_overhang_searcher(alpha))).1,
        };
        let ctx = Context {
            cfg,
            index,
            read: Read { window, reversed },
        };
        let mut engine = Engine {
            plain_read,
            plain,
            ambiguous,
            overhang,
            windows,
            mask,
            head_flags,
            tail_flags,
        };
        segments_with(ctx, &mut engine, acted)
    })
}

/// The search passes behind `adapter_segments`, over per-thread searchers.
/// `acted` receives the adapters that trimmed or excised, when given.
fn segments_with(
    ctx: Context<'_>,
    engine: &mut Engine<'_>,
    mut acted: Option<&mut [bool]>,
) -> Vec<(usize, usize)> {
    let cfg = ctx.cfg;
    let n = ctx.read.window.len();
    let mut tally = |keep: &Keep<'_>| {
        if let Some(acted) = acted.as_deref_mut() {
            for &adapter_idx in &keep.acted {
                acted[adapter_idx] = true;
            }
        }
    };
    let mut keep = Keep::new(cfg, ctx.index, n, cfg.split);
    search_terminal(ctx, (0, n), engine, &mut keep);
    if cfg.split {
        search_interior(ctx, engine, &mut keep);
    }
    // A trim near an end found by the interior search can anchor a deferred
    // terminal hit, so deferred hits are settled before the tally counts the
    // adapters that acted.
    keep.settle();
    tally(&keep);
    let (lo, hi, cuts) = keep.into_cuts(cfg.min_piece);
    if lo >= hi {
        return vec![];
    }
    if cuts.is_empty() {
        return vec![(lo, hi)];
    }

    // Each piece between cuts gets the terminal search again over its own
    // span, in ends-only mode. The end a cut created is a read end for every
    // purpose, and the opposite end has already been trimmed, so the second
    // pass changes it only when a further hit lies within `end_size` of the
    // new boundary.
    let mut segs = Vec::with_capacity(cuts.len() + 1);
    let mut cursor = lo;
    let mut push_piece = |s: usize, e: usize, segs: &mut Vec<(usize, usize)>| {
        if s >= e {
            return;
        }
        let mut keep = Keep::new(cfg, ctx.index, e - s, false);
        search_terminal(ctx, (s, e), engine, &mut keep);
        tally(&keep);
        if keep.lo < keep.hi {
            segs.push((s + keep.lo, s + keep.hi));
        }
    };
    for (s, e) in cuts {
        push_piece(cursor, s, &mut segs);
        cursor = cursor.max(e);
    }
    push_piece(cursor, hi, &mut segs);
    segs
}

#[cfg(test)]
mod segment_tests;
