//! Adapter, primer and barcode trimming.
//!
//! Searches each read window for catalog sequences with sassy, classifies every
//! accepted hit as a terminal trim or an interior excision, and returns the kept
//! spans. Presence detection, ab-initio inference and the ONT catalog live in
//! the submodules.

pub mod detect;
pub mod infer;
mod ont_catalog;
pub mod preset;
pub mod resolve;
pub mod search;

use std::borrow::Cow;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::OnceLock;

use aho_corasick::{AhoCorasick, AhoCorasickKind};
use search::{
    AmbiguousSearcher, BatchedAdapterSearcher, EncodedAdapterBatch, Hit, PlainSearcher, Strands,
    encode_patterns, encoded_pattern_hits, for_each_hit, is_plain_acgt, iupac_bases,
    new_ambiguous_searcher, new_batched_searcher, new_searcher, pattern_hits,
};

// One searcher of each kind per thread, reused across reads so
// `adapter_segments` does not allocate a searcher and its scratch buffers on
// every call. Per-thread state keeps the parallel workflows free of sharing.
thread_local! {
    /// The fast all-ACGT searcher. Used only when both the pattern and the
    /// searched text are plain ACGT; see `is_plain_acgt`.
    static RC_SEARCHER: RefCell<PlainSearcher> = RefCell::new(new_searcher());

    /// The ambiguity-tolerant searcher, used for a degenerate primer.
    static RC_AMBIGUOUS: RefCell<AmbiguousSearcher> = RefCell::new(new_ambiguous_searcher());

    /// The tiled searcher, for the pre-encoded terminal batches. Sassy
    /// implements pattern tiling only for the IUPAC profile, which on A/C/G/T
    /// input is equivalent to the DNA profile.
    static BATCH_SEARCHER: RefCell<BatchedAdapterSearcher> = RefCell::new(new_batched_searcher());

    /// Per-read buffers, reused across reads.
    static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default());
}

/// Per-read buffers of one thread. Each keeps its capacity across reads, so
/// the adapter stage itself allocates only the segments it returns; the
/// remaining per-read allocations are sassy's own, inside each search call.
#[derive(Debug, Default)]
struct Scratch {
    /// The normalized read, when the input is not its own normalization.
    normalized: Vec<u8>,
    /// The normalized read reversed, for the reverse strand of every search.
    reversed: Vec<u8>,
    /// Candidate windows; see `candidate_windows`.
    windows: WindowScratch,
}

/// Candidate windows of one read as `(adapter, start, end)` triples.
#[derive(Debug, Default)]
struct WindowScratch {
    /// Windows in seed-emission order.
    emitted: Vec<(usize, usize, usize)>,
    /// Windows grouped by adapter and merged: the interior search input.
    grouped: Vec<(usize, usize, usize)>,
    /// Per-adapter offsets into `grouped` while grouping.
    offsets: Vec<usize>,
}

/// The read end a catalog sequence is expected at. The tag gates terminal
/// trimming only: interior chimera splitting, when enabled, accepts any adapter
/// that matches in the read interior regardless of the tag, since a front or
/// rear adapter appearing mid-read is itself the chimera signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum End {
    /// Expected at the 5' end.
    Five,
    /// Expected at the 3' end.
    Three,
    /// Expected at either end.
    Both,
}

/// One searchable adapter, primer, barcode or flank sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    /// Display name, used in logs and the report.
    pub name: String,
    /// Nucleotide sequence; may contain IUPAC ambiguity codes.
    pub seq: Vec<u8>,
    /// Read end the sequence is expected at.
    pub end: End,
}

/// Resolved adapter-trimming settings for a run.
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    /// Sequences searched for, in configuration order.
    pub adapters: Vec<Adapter>,
    /// End-match tolerance as a fraction of adapter length (`k_end`).
    pub error_rate: f64,
    /// Bases at each end within which a hit is terminal (trim) rather than
    /// interior (split).
    pub end_size: usize,
    /// Whether interior adapters split the read; `false` is ends-only
    /// (`--adapter-ends-only`).
    pub split: bool,
    /// Exact-seed automaton for lossless whole-read candidate filtering, built
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
}

/// Returns the edit budget for a `len`-base pattern at `rate`, rounded down.
/// The epsilon keeps an integral product whose double lands below its integer.
pub(crate) fn edit_budget(rate: f64, len: usize) -> usize {
    (rate * len as f64 + 1e-9).floor() as usize
}

/// Edit budgets of one adapter: `k_end` for terminal hits and `k_mid`, half
/// of it, for interior hits.
#[derive(Debug, Clone, Copy)]
struct Budget {
    /// Pattern length in bases.
    len: usize,
    /// Edit budget for terminal hits.
    k_end: usize,
    /// Edit budget for interior hits.
    k_mid: usize,
}

impl Budget {
    /// Computes the budgets for a `len`-base pattern at `error_rate`.
    fn new(len: usize, error_rate: f64) -> Self {
        Self {
            len,
            k_end: edit_budget(error_rate, len),
            k_mid: edit_budget(0.5 * error_rate, len),
        }
    }
}

/// Upper bound on the plain strings one seed piece may expand to. A piece past
/// it marks its adapter `unfiltered`, and the interior search covers the whole
/// read for that adapter instead of candidate windows.
const MAX_SEED_EXPANSIONS: usize = 256;

/// Upper bound on the seed automaton's DFA size. A DFA past it is rebuilt as
/// a contiguous NFA, which is a few times slower per byte but grows with the
/// seed count rather than the state count times the alphabet stride.
const MAX_SEED_DFA_BYTES: usize = 64 << 20;

/// Exact-seed index over the adapter set. Aho-Corasick partition seeds bound
/// the interior search to candidate windows, and equal-length adapters are
/// grouped into SIMD batches for the terminal search.
#[derive(Debug, Clone)]
pub(crate) struct CandidateIndex {
    /// Automaton over every seed expansion (see `seed_automaton`); `None`
    /// when no adapter has seeds.
    matcher: Option<AhoCorasick>,
    /// Adapters owning each seed, indexed by automaton pattern id.
    seed_adapters: Vec<Vec<usize>>,
    /// Per-adapter edit budgets, computed once per adapter set.
    budgets: Vec<Budget>,
    /// Per-adapter `is_plain_acgt`, which selects the search profile.
    plain: Vec<bool>,
    /// Adapters with no usable seeds (see `MAX_SEED_EXPANSIONS`).
    unfiltered: Vec<bool>,
    /// Equal-length adapter groups searched together over the end windows.
    terminal_batches: Vec<TerminalBatch>,
    /// Adapters covered by a batch and skipped by the singleton search.
    batched_adapters: Vec<bool>,
}

/// The per-thread single-pattern searchers `search` chooses between.
struct Searchers<'a> {
    /// The DNA-profile searcher, for a plain pattern.
    plain: &'a mut PlainSearcher,
    /// The IUPAC-profile searcher, for a degenerate pattern.
    ambiguous: &'a mut AmbiguousSearcher,
}

/// Equal-length adapters searched together through sassy's pattern-parallel
/// API, which matches the forward and reverse-complement strands for the batch.
#[derive(Debug, Clone)]
struct TerminalBatch {
    /// Indices into `AdapterConfig::adapters`, in pattern order.
    adapter_indices: Vec<usize>,
    /// The adapter sequences, in `adapter_indices` order.
    patterns: Vec<Vec<u8>>,
    /// `patterns` encoded once for the tiled search; `None` past
    /// `MAX_TILED_PATTERN_LEN`, where the batch is searched from `patterns`.
    encoded: Option<EncodedAdapterBatch>,
    /// The shared pattern length.
    len: usize,
    /// The shared terminal edit budget.
    k_end: usize,
}

impl CandidateIndex {
    /// Builds the index for `adapters`; interior seeds are built only when
    /// `include_interior`.
    fn new(adapters: &[Adapter], error_rate: f64, include_interior: bool) -> Self {
        let budgets: Vec<Budget> = adapters
            .iter()
            .map(|adapter| Budget::new(adapter.seq.len(), error_rate))
            .collect();
        let plain: Vec<bool> = adapters
            .iter()
            .map(|adapter| is_plain_acgt(&adapter.seq))
            .collect();
        let mut unfiltered = vec![false; adapters.len()];
        let (matcher, seed_adapters) = if include_interior {
            let mut seeds: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
            for (adapter_idx, adapter) in adapters.iter().enumerate() {
                let Budget { len, k_mid, .. } = budgets[adapter_idx];
                if len < MIN_PATTERN_LEN {
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

            let patterns: Vec<Vec<u8>> = seeds.keys().cloned().collect();
            let seed_adapters: Vec<Vec<usize>> = seeds.into_values().collect();
            let matcher = (!patterns.is_empty()).then(|| seed_automaton(&patterns));
            (matcher, seed_adapters)
        } else {
            (None, Vec::new())
        };

        // Sassy packs equal-length patterns across SIMD lanes. Singletons stay
        // on the ordinary search path: a batch of one has no pattern-level
        // parallelism and is slower over these terminal windows.
        let mut by_len: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (adapter_idx, adapter) in adapters.iter().enumerate() {
            if adapter.seq.len() >= MIN_PATTERN_LEN {
                by_len
                    .entry(adapter.seq.len())
                    .or_default()
                    .push(adapter_idx);
            }
        }
        let mut terminal_batches = Vec::new();
        let mut batched_adapters = vec![false; adapters.len()];
        for (len, adapter_indices) in by_len {
            if adapter_indices.len() >= 2 {
                let batch_patterns: Vec<Vec<u8>> = adapter_indices
                    .iter()
                    .map(|&idx| adapters[idx].seq.clone())
                    .collect();
                for &adapter_idx in &adapter_indices {
                    batched_adapters[adapter_idx] = true;
                }
                let k_end = budgets[adapter_indices[0]].k_end;
                let encoded = encode_patterns(&batch_patterns);
                terminal_batches.push(TerminalBatch {
                    adapter_indices,
                    patterns: batch_patterns,
                    encoded,
                    len,
                    k_end,
                });
            }
        }
        Self {
            matcher,
            seed_adapters,
            budgets,
            plain,
            unfiltered,
            terminal_batches,
            batched_adapters,
        }
    }

    /// Fills `scratch.grouped` with the text spans that can hold an interior
    /// hit, as `(adapter, start, end)` sorted by adapter then start: a radius
    /// around every exact seed occurrence, merged per adapter, or the whole
    /// text for an `unfiltered` adapter. `text` is a `normalized_read`; the
    /// seed automaton matches uppercase bases only.
    fn candidate_windows(&self, text: &[u8], scratch: &mut WindowScratch) {
        let WindowScratch {
            emitted,
            grouped,
            offsets,
        } = scratch;
        emitted.clear();
        for (adapter_idx, &whole) in self.unfiltered.iter().enumerate() {
            if whole {
                emitted.push((adapter_idx, 0, text.len()));
            }
        }
        if let Some(matcher) = &self.matcher {
            for m in matcher.find_overlapping_iter(text) {
                for &adapter_idx in &self.seed_adapters[m.pattern().as_usize()] {
                    let Budget { len, k_end, .. } = self.budgets[adapter_idx];
                    // The exact seed lies inside the `<= k_mid` alignment. A
                    // radius of pattern length + `k_end` on each side contains
                    // that entire alignment and enough context for the
                    // full-window `k_end` search to reproduce its span and tie
                    // behavior.
                    let radius = len + k_end;
                    emitted.push((
                        adapter_idx,
                        m.start().saturating_sub(radius),
                        m.end().saturating_add(radius).min(text.len()),
                    ));
                }
            }
        }

        // Groups the windows per adapter with a counting scatter and sorts
        // each run on its own. The automaton emits in text order, so a run
        // is nearly sorted already and costs a short insertion pass, where
        // one sort of every window costs a comparison sort of hundreds.
        let adapters = self.budgets.len();
        offsets.clear();
        offsets.resize(adapters + 1, 0);
        for &(adapter_idx, _, _) in emitted.iter() {
            offsets[adapter_idx + 1] += 1;
        }
        for adapter_idx in 0..adapters {
            offsets[adapter_idx + 1] += offsets[adapter_idx];
        }
        grouped.clear();
        grouped.resize(emitted.len(), (0, 0, 0));
        for &window in emitted.iter() {
            let slot = &mut offsets[window.0];
            grouped[*slot] = window;
            *slot += 1;
        }

        // Merges overlapping or touching windows of one adapter in place.
        // After the scatter `offsets[a]` is the end of run `a`, and the
        // merged prefix never reaches the run being read.
        let mut merged = 0;
        let mut run_start = 0;
        for &run_end in offsets[..adapters].iter() {
            grouped[run_start..run_end].sort_unstable();
            let first = merged;
            for i in run_start..run_end {
                let (_, start, end) = grouped[i];
                if merged > first && start <= grouped[merged - 1].2 {
                    grouped[merged - 1].2 = grouped[merged - 1].2.max(end);
                } else {
                    grouped[merged] = grouped[i];
                    merged += 1;
                }
            }
            run_start = run_end;
        }
        grouped.truncate(merged);
    }
}

/// Builds the overlapping-match automaton over `seeds`. The seed scan runs
/// over every base of every read, so the automaton is a DFA: one table lookup
/// per byte, against the contiguous NFA's per-state transition scan. Byte
/// classes fold the bytes outside ACGT into one column and keep the table
/// small. The automaton is case-sensitive over uppercase seeds: the scanned
/// text is always a `normalized_read`, and a case-insensitive alphabet would
/// double the table stride for nothing. A DFA past `MAX_SEED_DFA_BYTES`, or
/// one the builder rejects, gives way to the contiguous NFA.
fn seed_automaton(seeds: &[Vec<u8>]) -> AhoCorasick {
    let build = |kind: AhoCorasickKind| {
        AhoCorasick::builder()
            .kind(Some(kind))
            .byte_classes(true)
            .prefilter(true)
            .build(seeds)
    };
    match build(AhoCorasickKind::DFA) {
        Ok(dfa) if dfa.memory_usage() <= MAX_SEED_DFA_BYTES => dfa,
        _ => build(AhoCorasickKind::ContiguousNFA)
            .expect("Adapter seeds are nonempty ASCII DNA patterns"),
    }
}

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
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Minimum searchable pattern length: a shorter pattern matches almost anywhere
/// under any error budget and is never searched standalone. Catalog flanks
/// below it are omitted from the catalog for the same reason.
pub const MIN_PATTERN_LEN: usize = 11;

/// Outboard flank length at or below which a hit covered by both end zones is
/// a terminal trim rather than an excision. A flank this short is adapter
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
/// fills one lane of eight with a single pattern, and it measured slower on
/// these short windows than the allocations cost.
///
/// `text` is plain ACGT on every call; see `normalized_read`.
fn search(
    searchers: &mut Searchers<'_>,
    index: &CandidateIndex,
    adapter_idx: usize,
    pattern: &[u8],
    text: Strands<'_>,
    k: usize,
    accept: impl FnMut(Hit),
) {
    if index.plain[adapter_idx] {
        for_each_hit(searchers.plain, pattern, &text, k, accept);
    } else {
        for_each_hit(searchers.ambiguous, pattern, &text, k, accept);
    }
}

/// Returns the strands of `read[start..end]`, given `reversed` as `read`
/// reversed.
fn strands<'a>(read: &'a [u8], reversed: &'a [u8], start: usize, end: usize) -> Strands<'a> {
    let n = read.len();
    Strands {
        forward: &read[start..end],
        reversed: &reversed[n - end..n - start],
    }
}

/// The base an ambiguity code in a read is rewritten to before searching.
///
/// An uncalled base is evidence of nothing, so it consumes error budget rather
/// than matching for free: a short `N` inside a real adapter still matches
/// within `--adapter-error-rate`, while a run of them never looks like an
/// adapter. Leaving the codes in place and searching on the IUPAC profile would
/// do the opposite and excise the whole run as adapter.
///
/// Any ACGT byte serves; the requirement is one fixed base, so that no real
/// adapter matches a homopolymer run of it within its budget. The rewrite also
/// keeps the fast profile usable, which panics during traceback on any other
/// byte, and keeps the batched path consistent, since sassy implements batching
/// only for IUPAC.
const AMBIGUOUS_READ_BASE: u8 = b'A';

/// Returns the read as every searcher sees it: uppercase, with each byte
/// outside ACGT rewritten to `AMBIGUOUS_READ_BASE`. An uppercase plain read,
/// the common case, borrows unchanged and allocates nothing. Sassy's profiles
/// fold case themselves; the seed automaton does not, so the text is folded
/// once here rather than on every seed transition.
pub(crate) fn normalized_read(window: &[u8]) -> Cow<'_, [u8]> {
    if is_upper_acgt(window) {
        return Cow::Borrowed(window);
    }
    Cow::Owned(window.iter().map(|&b| normalize_base(b)).collect())
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
/// its own `normalized_read`. Folded 32 bytes at a time without an early exit
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
pub enum HitAction {
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

/// Classifies a hit at window coordinates `[start, end)` in a length-`n` window.
///
/// A hit inside one end zone only trims that end when the adapter's tag allows
/// it. A hit covered by both end zones (`n <= end_size + hit length`) is placed
/// by geometry alone: every search is reverse-complement aware and catalog rear
/// entries are reverse complements of front entries, so the same span carries
/// hits of both tags there and the tag says nothing about which end it is at.
/// Such a hit trims toward an end whose outboard flank is at most `FLANK_SLACK`
/// (the nearer one when both are), and otherwise becomes `Excise`: cut out the
/// adapter, keep both flanks, as a central chimera junction needs.
fn classify_terminal(start: usize, end: usize, n: usize, end_size: usize, tag: End) -> Terminal {
    let in_head = start <= end_size;
    let in_tail = end >= n.saturating_sub(end_size);
    match (in_head, in_tail) {
        (true, true) => match (start <= FLANK_SLACK, n - end <= FLANK_SLACK) {
            (false, false) => Terminal::Excise,
            (true, false) => Terminal::Five,
            (false, true) => Terminal::Three,
            (true, true) => nearer_end(start, end, n),
        },
        (true, false) if matches!(tag, End::Five | End::Both) => Terminal::Five,
        (false, true) if matches!(tag, End::Three | End::Both) => Terminal::Three,
        _ => Terminal::None,
    }
}

/// Classifies a hit for ends-only mode: splitting is disabled, so an `Excise`
/// outcome of `classify_terminal` resolves to a terminal trim toward the nearer
/// end, and every other outcome is unchanged.
fn ends_only_terminal(start: usize, end: usize, n: usize, end_size: usize, tag: End) -> Terminal {
    match classify_terminal(start, end, n, end_size, tag) {
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
/// interior cuts.
struct Keep<'a> {
    /// The configured adapters, for tags and names.
    adapters: &'a [Adapter],
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
    /// Accepted excisions, merged by `into_segments`.
    interior: Vec<(usize, usize)>,
}

impl<'a> Keep<'a> {
    /// Creates an accumulator that keeps the whole `[0, n)` window.
    fn new(adapters: &'a [Adapter], n: usize, end_size: usize, split: bool) -> Self {
        Self {
            adapters,
            n,
            end_size,
            split,
            lo: 0,
            hi: n,
            interior: Vec::new(),
        }
    }

    /// Classifies one hit and applies it when `site` owns the outcome.
    fn accept(&mut self, site: Site, adapter_idx: usize, start: usize, end: usize, cost: usize) {
        let adapter = &self.adapters[adapter_idx];
        let terminal = if self.split {
            classify_terminal(start, end, self.n, self.end_size, adapter.end)
        } else {
            ends_only_terminal(start, end, self.n, self.end_size, adapter.end)
        };
        let action = match (site, terminal) {
            (Site::Head, Terminal::Five) => HitAction::TrimFivePrime,
            (Site::Head, Terminal::Excise) => HitAction::Excise,
            (Site::Tail { .. }, Terminal::Three) => HitAction::TrimThreePrime,
            (Site::Interior, Terminal::None) => HitAction::Excise,
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
        trace_hit(&adapter.name, start, end, cost, Some(action));
        match action {
            HitAction::TrimFivePrime => self.lo = self.lo.max(end),
            HitAction::TrimThreePrime => self.hi = self.hi.min(start),
            HitAction::Excise => self.interior.push((start, end)),
        }
    }

    /// Returns the keep segments: `[lo, hi)` with the merged interior cuts
    /// carved out.
    fn into_segments(self) -> Vec<(usize, usize)> {
        let Keep {
            lo, hi, interior, ..
        } = self;
        if lo >= hi {
            return vec![];
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
            if let Some(last) = merged.last_mut()
                && s <= last.1
            {
                last.1 = last.1.max(e);
                continue;
            }
            merged.push((s, e));
        }
        let mut segs = Vec::new();
        let mut cursor = lo;
        for (s, e) in merged {
            if s > cursor {
                segs.push((cursor, s));
            }
            cursor = cursor.max(e);
        }
        if cursor < hi {
            segs.push((cursor, hi));
        }
        segs
    }
}

/// Searches every equal-length batch over the two end windows. All adapters in
/// a batch share a length and budget, so the windows are shared too; this
/// collapses the ONT catalog's 96 equal-length barcode searches into one SIMD
/// pattern search per end. `reversed` is `window` reversed.
fn search_batched(
    index: &CandidateIndex,
    window: &[u8],
    reversed: &[u8],
    searcher: &mut BatchedAdapterSearcher,
    keep: &mut Keep<'_>,
) {
    let n = window.len();
    for batch in &index.terminal_batches {
        let (head_end, tail_start) = terminal_windows(n, keep.end_size, batch.len, batch.k_end);
        accept_batch_hits(
            batch,
            searcher,
            &window[..head_end],
            &reversed[n - head_end..],
            0,
            Site::Head,
            keep,
        );
        accept_batch_hits(
            batch,
            searcher,
            &window[tail_start..],
            &reversed[..n - tail_start],
            tail_start,
            Site::Tail { head_end },
            keep,
        );
    }
}

/// Searches one batch over `text`, a window starting at `offset` in the read
/// whose reversal is `reversed`, and passes every hit to `keep` at `site` in
/// read coordinates.
fn accept_batch_hits(
    batch: &TerminalBatch,
    searcher: &mut BatchedAdapterSearcher,
    text: &[u8],
    reversed: &[u8],
    offset: usize,
    site: Site,
    keep: &mut Keep<'_>,
) {
    let mut accept = |pattern_idx: usize, start: usize, end: usize, cost: usize| {
        keep.accept(
            site,
            batch.adapter_indices[pattern_idx],
            offset + start,
            offset + end,
            cost,
        );
    };
    match &batch.encoded {
        Some(encoded) => {
            encoded_pattern_hits(searcher, encoded, text, reversed, batch.k_end, accept);
        },
        None => {
            for h in pattern_hits(searcher, &batch.patterns, text, batch.k_end) {
                accept(h.pattern_idx, h.text_start, h.text_end, h.cost as usize);
            }
        },
    }
}

/// Searches every adapter without an equal-length partner over the two end
/// windows, one pattern at a time. `reversed` is `window` reversed.
fn search_singletons(
    cfg: &AdapterConfig,
    index: &CandidateIndex,
    window: &[u8],
    reversed: &[u8],
    searchers: &mut Searchers<'_>,
    keep: &mut Keep<'_>,
) {
    let n = window.len();
    for (adapter_idx, adapter) in cfg.adapters.iter().enumerate() {
        let Budget { len, k_end, .. } = index.budgets[adapter_idx];
        if len < MIN_PATTERN_LEN || index.batched_adapters[adapter_idx] {
            continue;
        }
        let (head_end, tail_start) = terminal_windows(n, keep.end_size, len, k_end);
        search(
            searchers,
            index,
            adapter_idx,
            &adapter.seq,
            strands(window, reversed, 0, head_end),
            k_end,
            |h| keep.accept(Site::Head, adapter_idx, h.start, h.end, h.cost),
        );
        search(
            searchers,
            index,
            adapter_idx,
            &adapter.seq,
            strands(window, reversed, tail_start, n),
            k_end,
            |h| {
                keep.accept(
                    Site::Tail { head_end },
                    adapter_idx,
                    tail_start + h.start,
                    tail_start + h.end,
                    h.cost,
                );
            },
        );
    }
}

/// Searches every adapter's candidate windows at `k_mid`. Exact partition
/// seeds identify every possible interior match, and interior hits are
/// accepted only up to `k_mid`, so the search runs at that limit rather than
/// the looser end budget. `reversed` is `window` reversed. An adapter below
/// `MIN_PATTERN_LEN` has no seeds and no windows.
fn search_interior(
    cfg: &AdapterConfig,
    index: &CandidateIndex,
    window: &[u8],
    reversed: &[u8],
    windows: &mut WindowScratch,
    searchers: &mut Searchers<'_>,
    keep: &mut Keep<'_>,
) {
    index.candidate_windows(window, windows);
    for &(adapter_idx, start, end) in windows.grouped.iter() {
        let Budget { k_mid, .. } = index.budgets[adapter_idx];
        search(
            searchers,
            index,
            adapter_idx,
            &cfg.adapters[adapter_idx].seq,
            strands(window, reversed, start, end),
            k_mid,
            |h| {
                keep.accept(
                    Site::Interior,
                    adapter_idx,
                    start + h.start,
                    start + h.end,
                    h.cost,
                );
            },
        );
    }
}

/// Computes the adapter keep segments for `window`: terminal hits within
/// `end_size` of an end trim that end inward, and interior hits (at the
/// stricter `k_mid`) excise and split.
///
/// Under `--adapter-ends-only` (`cfg.split` false) only the two end zones are
/// searched, since no interior hit could be acted on.
///
/// Returns `[start, end)` spans in `window` coordinates.
pub fn adapter_segments(window: &[u8], cfg: &AdapterConfig) -> Vec<(usize, usize)> {
    let n = window.len();
    if n == 0 {
        return vec![];
    }
    if cfg.adapters.is_empty() {
        return vec![(0, n)];
    }
    let index = cfg
        .candidate_index
        .get_or_init(|| CandidateIndex::new(&cfg.adapters, cfg.error_rate, cfg.split));
    let mut keep = Keep::new(&cfg.adapters, n, cfg.end_size.min(n), cfg.split);
    SCRATCH.with_borrow_mut(|scratch| {
        let Scratch {
            normalized,
            reversed,
            windows,
        } = scratch;
        let window: &[u8] = if is_upper_acgt(window) {
            window
        } else {
            normalized.clear();
            normalized.extend(window.iter().map(|&b| normalize_base(b)));
            normalized
        };
        reversed.clear();
        reversed.extend(window.iter().rev());
        if !index.terminal_batches.is_empty() {
            BATCH_SEARCHER.with_borrow_mut(|tiled| {
                search_batched(index, window, reversed, tiled, &mut keep);
            });
        }
        RC_SEARCHER.with_borrow_mut(|plain| {
            RC_AMBIGUOUS.with_borrow_mut(|ambiguous| {
                let mut searchers = Searchers { plain, ambiguous };
                search_singletons(cfg, index, window, reversed, &mut searchers, &mut keep);
                if cfg.split {
                    search_interior(
                        cfg,
                        index,
                        window,
                        reversed,
                        windows,
                        &mut searchers,
                        &mut keep,
                    );
                }
            });
        });
    });

    keep.into_segments()
}

#[cfg(test)]
mod segment_tests {
    use super::preset::preset_ont;
    use super::*;

    /// Builds a configuration at error rate 0.2 and `end_size` 20.
    fn cfg(adapters: Vec<Adapter>, split: bool) -> AdapterConfig {
        AdapterConfig {
            adapters,
            error_rate: 0.2,
            end_size: 20,
            split,
            candidate_index: std::sync::OnceLock::new(),
        }
    }

    /// Builds an adapter from its parts.
    fn ad(name: &str, seq: &[u8], end: End) -> Adapter {
        Adapter {
            name: name.into(),
            seq: seq.to_vec(),
            end,
        }
    }

    /// Generates deterministic SplitMix64 bases with the same generator the
    /// inference fixtures use.
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

    /// Searches `pattern` over a plain slice through the profile its alphabet
    /// selects, as the reference for `search`.
    fn reference_search(
        plain: &mut PlainSearcher,
        ambiguous: &mut AmbiguousSearcher,
        pattern: &[u8],
        text: &[u8],
        k: usize,
    ) -> Vec<Hit> {
        if is_plain_acgt(pattern) {
            search::hits(plain, pattern, text, k)
        } else {
            search::hits(ambiguous, pattern, text, k)
        }
    }

    /// Returns `candidate_windows` grouped per adapter.
    fn windows_by_adapter(index: &CandidateIndex, text: &[u8]) -> Vec<Vec<(usize, usize)>> {
        let mut scratch = WindowScratch::default();
        index.candidate_windows(text, &mut scratch);
        let mut by_adapter = vec![Vec::new(); index.budgets.len()];
        for (adapter_idx, start, end) in scratch.grouped {
            by_adapter[adapter_idx].push((start, end));
        }
        by_adapter
    }

    /// Computes the segments by exhaustive full-window search, as the reference
    /// for the candidate filter.
    fn reference_segments(window: &[u8], cfg: &AdapterConfig) -> Vec<(usize, usize)> {
        let n = window.len();
        if n == 0 {
            return vec![];
        }
        if cfg.adapters.is_empty() {
            return vec![(0, n)];
        }
        let window = normalized_read(window);
        let end_size = cfg.end_size.min(n);
        let mut lo = 0usize;
        let mut hi = n;
        let mut interior = Vec::new();
        let mut plain = new_searcher();
        let mut ambiguous = new_ambiguous_searcher();

        for adapter in &cfg.adapters {
            let len = adapter.seq.len();
            if len < MIN_PATTERN_LEN {
                continue;
            }
            let Budget { k_end, k_mid, .. } = Budget::new(len, cfg.error_rate);
            if cfg.split {
                for hit in
                    reference_search(&mut plain, &mut ambiguous, &adapter.seq, &window, k_end)
                {
                    match classify_terminal(hit.start, hit.end, n, end_size, adapter.end) {
                        Terminal::Five => lo = lo.max(hit.end),
                        Terminal::Three => hi = hi.min(hit.start),
                        Terminal::Excise => interior.push((hit.start, hit.end)),
                        Terminal::None if hit.cost <= k_mid => {
                            interior.push((hit.start, hit.end));
                        },
                        Terminal::None => {},
                    }
                }
            } else {
                let head_end = (end_size + len + k_end).min(n);
                for hit in reference_search(
                    &mut plain,
                    &mut ambiguous,
                    &adapter.seq,
                    &window[..head_end],
                    k_end,
                ) {
                    if ends_only_terminal(hit.start, hit.end, n, end_size, adapter.end)
                        == Terminal::Five
                    {
                        lo = lo.max(hit.end);
                    }
                }
                let tail_start = n.saturating_sub(end_size + len + k_end);
                for hit in reference_search(
                    &mut plain,
                    &mut ambiguous,
                    &adapter.seq,
                    &window[tail_start..],
                    k_end,
                ) {
                    let (start, end) = (tail_start + hit.start, tail_start + hit.end);
                    if ends_only_terminal(start, end, n, end_size, adapter.end) == Terminal::Three {
                        hi = hi.min(start);
                    }
                }
            }
        }

        if lo >= hi {
            return vec![];
        }
        let mut cuts: Vec<_> = interior
            .into_iter()
            .filter_map(|(start, end)| {
                let start = start.max(lo);
                let end = end.min(hi);
                (start < end).then_some((start, end))
            })
            .collect();
        cuts.sort_unstable();
        let mut merged: Vec<(usize, usize)> = Vec::new();
        for (start, end) in cuts {
            if let Some(last) = merged.last_mut()
                && start <= last.1
            {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        }
        let mut segments = Vec::new();
        let mut cursor = lo;
        for (start, end) in merged {
            if start > cursor {
                segments.push((cursor, start));
            }
            cursor = cursor.max(end);
        }
        if cursor < hi {
            segments.push((cursor, hi));
        }
        segments
    }

    /// Linear congruential generator for deterministic fixtures.
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> usize {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 32) as usize
        }
        fn below(&mut self, n: usize) -> usize {
            self.next() % n
        }
        fn dna(&mut self, n: usize) -> Vec<u8> {
            (0..n).map(|_| b"ACGT"[self.below(4)]).collect()
        }
    }

    /// Plants one adapter with up to `k_mid` edits at an end or in the
    /// interior of a random window and checks the candidate search against the
    /// full-window reference. `degenerate` rewrites a share of adapter
    /// positions to ambiguity codes; the planted copy then carries one base
    /// each code stands for.
    fn check_candidate_search_randomized(seed: u64, degenerate: bool) {
        const CODES: &[u8] = b"RYSWKMBDHVN";
        let mut rng = Lcg(seed);
        for case in 0..400 {
            let adapters: Vec<Adapter> = (0..(1 + rng.below(10)))
                .map(|i| {
                    let len = 11 + rng.below(40);
                    let mut seq = rng.dna(len);
                    if degenerate {
                        for _ in 0..(1 + len / 8) {
                            let p = rng.below(len);
                            seq[p] = CODES[rng.below(CODES.len())];
                        }
                    }
                    Adapter {
                        name: format!("a{i}"),
                        seq,
                        end: match rng.below(3) {
                            0 => End::Five,
                            1 => End::Three,
                            _ => End::Both,
                        },
                    }
                })
                .collect();
            let window_len = 80 + rng.below(660);
            let mut window = rng.dna(window_len);

            let planted = rng.below(adapters.len());
            let pattern = &adapters[planted].seq;
            if pattern.len() <= window.len() {
                let max_edits = (0.1 * pattern.len() as f64).floor() as usize;
                let mut copy: Vec<u8> = pattern
                    .iter()
                    .map(|&code| {
                        let bases = iupac_bases(code).expect("adapter bytes are nucleotide codes");
                        bases[rng.below(bases.len())]
                    })
                    .collect();
                for _ in 0..rng.below(max_edits + 1) {
                    match rng.below(3) {
                        0 => {
                            let p = rng.below(copy.len());
                            let old = copy[p];
                            copy[p] =
                                b"ACGT"[(b"ACGT".iter().position(|&b| b == old).unwrap() + 1) % 4];
                        },
                        1 => {
                            let p = rng.below(copy.len() + 1);
                            copy.insert(p, b"ACGT"[rng.below(4)]);
                        },
                        _ => {
                            let p = rng.below(copy.len());
                            copy.remove(p);
                        },
                    }
                }
                let planted_len = copy.len();
                let pos = match rng.below(3) {
                    0 => rng.below(8.min(window.len() - planted_len + 1)),
                    1 => {
                        window.len()
                            - planted_len
                            - rng.below(8.min(window.len() - planted_len + 1))
                    },
                    _ => rng.below(window.len() - planted_len + 1),
                };
                window[pos..pos + planted_len].copy_from_slice(&copy);
                if case % 7 == 0 {
                    window.make_ascii_lowercase();
                }
            }

            let cfg = AdapterConfig {
                adapters,
                error_rate: 0.2,
                end_size: 1 + rng.below(180),
                split: true,
                candidate_index: std::sync::OnceLock::new(),
            };
            assert_eq!(
                adapter_segments(&window, &cfg),
                reference_segments(&window, &cfg),
                "Candidate/reference mismatch in randomized case {case} (degenerate: {degenerate})"
            );
        }
    }

    /// The candidate search matches the full-window reference on random plain
    /// adapters.
    #[test]
    fn candidate_search_matches_full_search_randomized() {
        check_candidate_search_randomized(0x4e4f_4f44_4c45_5301, false);
    }

    /// The candidate search matches the full-window reference on random
    /// degenerate adapters.
    #[test]
    fn candidate_search_matches_full_search_randomized_degenerate() {
        check_candidate_search_randomized(0x4445_4745_4e45_5241, true);
    }

    /// A window holding an `N` run matches the reference after normalization.
    #[test]
    fn non_acgt_window_falls_back_to_scalar_search() {
        let cfg = cfg(
            vec![
                ad("a", b"ACGTACGTACGT", End::Five),
                ad("b", b"TTTTGGGGCCCC", End::Three),
            ],
            true,
        );
        let window = b"ACGTACGTACGTNNNNNNNNNNNNNNNNNNNNTTTTGGGGCCCC";
        assert_eq!(
            adapter_segments(window, &cfg),
            reference_segments(window, &cfg)
        );
    }

    /// A copy within the edit budget always retains one exact seed piece.
    #[test]
    fn partition_seeds_survive_random_indels_and_substitutions() {
        struct Lcg(u64);
        impl Lcg {
            fn next(&mut self) -> usize {
                self.0 = self
                    .0
                    .wrapping_mul(2862933555777941757)
                    .wrapping_add(3037000493);
                (self.0 >> 32) as usize
            }
            fn below(&mut self, n: usize) -> usize {
                self.next() % n
            }
            fn base(&mut self) -> u8 {
                b"ACGT"[self.below(4)]
            }
        }

        let mut rng = Lcg(0x5049_4745_4f4e_484f);
        for case in 0..1000 {
            let pattern: Vec<u8> = (0..(11 + rng.below(50))).map(|_| rng.base()).collect();
            let k = (0.1 * pattern.len() as f64).floor() as usize;
            let mut mutated = pattern.clone();
            for _ in 0..rng.below(k + 1) {
                match rng.below(3) {
                    0 => {
                        let p = rng.below(mutated.len());
                        mutated[p] = rng.base();
                    },
                    1 => {
                        let p = rng.below(mutated.len() + 1);
                        mutated.insert(p, rng.base());
                    },
                    _ if mutated.len() > 1 => {
                        let p = rng.below(mutated.len());
                        mutated.remove(p);
                    },
                    _ => {},
                }
            }
            let adapter = Adapter {
                name: "a".into(),
                seq: pattern,
                end: End::Both,
            };
            let index = CandidateIndex::new(&[adapter], 0.2, true);
            let mut text: Vec<u8> = (0..17).map(|_| rng.base()).collect();
            text.extend_from_slice(&mutated);
            text.extend((0..19).map(|_| rng.base()));
            if case % 2 == 0 {
                text.make_ascii_lowercase();
            }
            assert!(
                !windows_by_adapter(&index, &normalized_read(&text))[0].is_empty(),
                "Lossless seed filter rejected <=k edit case {case}"
            );
        }
    }

    /// `expand_iupac` enumerates every concrete string and refuses to expand
    /// past the cap or over a non-nucleotide byte.
    #[test]
    fn expand_iupac_enumerates_every_base_and_caps() {
        assert_eq!(expand_iupac(b"AC"), Some(vec![b"AC".to_vec()]));
        let mut r = expand_iupac(b"RY").unwrap();
        r.sort();
        assert_eq!(
            r,
            vec![
                b"AC".to_vec(),
                b"AT".to_vec(),
                b"GC".to_vec(),
                b"GT".to_vec()
            ]
        );
        assert_eq!(expand_iupac(b"NNNN").map(|v| v.len()), Some(256));
        assert_eq!(expand_iupac(b"NNNNN"), None, "Five N's expand past the cap");
        assert_eq!(
            expand_iupac(b"ACXT"),
            None,
            "A non-nucleotide byte has no expansion"
        );
    }

    /// Reverse complementing twice is the identity for every IUPAC code, with
    /// case preserved.
    #[test]
    fn reverse_complement_uses_the_full_iupac_table() {
        for &code in b"ACGTRYSWKMBDHVNacgtryswkmbdhvn" {
            let once = reverse_complement(&[code]);
            assert_eq!(
                reverse_complement(&once),
                vec![code],
                "Code {}",
                code as char
            );
            assert_eq!(
                once[0].is_ascii_lowercase(),
                code.is_ascii_lowercase(),
                "Case is preserved for {}",
                code as char
            );
        }
        assert_eq!(reverse_complement(b"RYKMBDHV"), b"BDHVKMRY");
        assert_eq!(reverse_complement(b"SWN"), b"NWS");
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
    }

    /// `edit_budget` does not round an integral product down through
    /// floating-point error.
    #[test]
    fn edit_budget_keeps_an_integral_product() {
        assert_eq!(edit_budget(0.29, 100), 29);
        assert_eq!(edit_budget(0.57, 100), 57);
        assert_eq!(edit_budget(0.2, 22), 4);
        assert_eq!(edit_budget(0.1, 12), 1);
        assert_eq!(edit_budget(0.0, 50), 0);
    }

    /// Both seed pieces of the 12-mer hold an `N`; the planted copy is one
    /// concrete instance and still splits the read.
    #[test]
    fn degenerate_adapter_splits_interior_chimera() {
        let adapter = b"GTNGTTGGNTGT";
        let mut w = vec![b'A'; 40];
        w.extend_from_slice(b"GTGGTTGGGTGT");
        w.extend_from_slice(&[b'C'; 40]);
        let c = AdapterConfig {
            adapters: vec![ad("deg", adapter, End::Both)],
            error_rate: 0.2,
            end_size: 10,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        assert_eq!(adapter_segments(&w, &c), vec![(0, 40), (52, 92)]);
    }

    /// The `N` sits in the first piece and the substitution (C to A) in the
    /// second, so only an expanded first-piece seed finds the copy.
    #[test]
    fn degenerate_adapter_splits_with_one_substitution_in_the_plain_piece() {
        let adapter = b"GTNGTTGGCTGT";
        let mut w = vec![b'A'; 40];
        w.extend_from_slice(b"GTGGTTGGATGT");
        w.extend_from_slice(&[b'C'; 40]);
        let c = AdapterConfig {
            adapters: vec![ad("deg", adapter, End::Both)],
            error_rate: 0.2,
            end_size: 10,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        assert_eq!(adapter_segments(&w, &c), vec![(0, 40), (52, 92)]);
    }

    /// Each 8-base piece holds five `N`s (1024 expansions), past the cap, so
    /// the adapter is searched over the whole window and still splits it.
    #[test]
    fn overly_degenerate_adapter_is_searched_unfiltered() {
        let adapter = b"NNNNNGGTTGGNNNNN";
        let mut w = vec![b'A'; 35];
        w.extend_from_slice(b"CACGTGGTTGGACGTC");
        w.extend_from_slice(&[b'C'; 41]);
        let c = AdapterConfig {
            adapters: vec![ad("deg", adapter, End::Both)],
            error_rate: 0.2,
            end_size: 10,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let index = CandidateIndex::new(&c.adapters, c.error_rate, true);
        assert_eq!(index.unfiltered, vec![true]);
        assert!(
            index.matcher.is_none(),
            "No seed is built for an unfiltered adapter"
        );
        assert_eq!(windows_by_adapter(&index, &w), vec![vec![(0, w.len())]]);
        let segs = adapter_segments(&w, &c);
        assert_eq!(segs, reference_segments(&w, &c));
        assert_eq!(
            segs.len(),
            2,
            "The unfiltered adapter still splits the read"
        );
    }

    /// An empty adapter set keeps the whole window.
    #[test]
    fn no_adapters_is_identity() {
        let w = b"ACGTACGTACGTACGT";
        assert_eq!(adapter_segments(w, &cfg(vec![], true)), vec![(0, w.len())]);
    }

    /// A 5' adapter at the read start is trimmed.
    #[test]
    fn trims_5prime_adapter_and_outboard() {
        let adapter = b"ACGTACGTACGT"; // 12 bp
        let mut w = adapter.to_vec();
        w.extend_from_slice(b"AAAAAAAAAAAA");
        let c = cfg(vec![ad("a", adapter, End::Five)], false);
        assert_eq!(adapter_segments(&w, &c), vec![(12, 24)]);
    }

    /// With the default `end_size` of 150 the two end zones overlap for any
    /// read of at most 300 bp. A chimera-junction adapter within `end_size` of
    /// both ends splits the read and keeps both inserts; treating it as a
    /// terminal adapter would discard the entire outboard arm, up to `end_size`
    /// bases of real insert.
    #[test]
    fn central_chimera_on_short_read_splits_both_arms() {
        let adapter = b"GGGGTTTTGGGGTTTTGGGG"; // 20bp, G/T only (no A/C to collide)
        let mut w = vec![b'A'; 115]; // insert1
        let cut = w.len();
        w.extend_from_slice(adapter); // junction adapter at [115,135)
        w.extend_from_slice(&[b'C'; 115]); // insert2 -> n=250
        let c = AdapterConfig {
            adapters: vec![ad("mid", adapter, End::Both)],
            error_rate: 0.2,
            end_size: 150, // default: end-zones overlap on this 250bp read
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let segs = adapter_segments(&w, &c);
        assert_eq!(
            segs,
            vec![(0, cut), (cut + adapter.len(), w.len())],
            "Central chimera must split into both arms, not lose insert1"
        );
    }

    /// Three junk bases, the adapter, then the insert: the 3 bp flank is
    /// trimmed with the adapter instead of surviving as its own segment.
    #[test]
    fn near_terminal_excision_folds_into_a_trim() {
        let adapter = b"GGGGTTTTGGGGTTTTGGGG";
        let mut w = b"AAA".to_vec();
        w.extend_from_slice(adapter);
        w.extend_from_slice(&[b'C'; 37]);
        let c = AdapterConfig {
            adapters: vec![ad("a", adapter, End::Both)],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        assert_eq!(adapter_segments(&w, &c), vec![(23, 60)]);

        // Mirror: insert, adapter, three junk bases.
        let mut w = vec![b'C'; 37];
        w.extend_from_slice(adapter);
        w.extend_from_slice(b"AAA");
        assert_eq!(adapter_segments(&w, &c), vec![(0, 37)]);
    }

    /// Asserts one kept segment that reaches the far end, with the trim boundary
    /// at the adapter or at most `slack` bases into the insert.
    fn assert_insert_kept(
        segs: &[(usize, usize)],
        n: usize,
        adapter: usize,
        slack: usize,
        front: bool,
    ) {
        assert_eq!(segs.len(), 1, "Read length {n}: {segs:?}");
        let (start, end) = segs[0];
        if front {
            assert_eq!(end, n, "Read length {n}: {segs:?}");
            assert!(
                (adapter..=adapter + slack).contains(&start),
                "Read length {n}: {segs:?}"
            );
        } else {
            assert_eq!(start, 0, "Read length {n}: {segs:?}");
            let insert = n - adapter;
            assert!(
                (insert - slack..=insert).contains(&end),
                "Read length {n}: {segs:?}"
            );
        }
    }

    /// The catalog holds longer entries that begin with `PCR1_front`
    /// (`cDNA_rear` adds `TTT`), so a fuzzy hit can carry the trim up to three
    /// bases into the insert. The same insert padded to 200 bp, where the end
    /// zones do not overlap, keeps the insert as well.
    #[test]
    fn short_read_with_front_adapter_keeps_insert_under_ont_preset() {
        // PCR1_front (5') plus a 140 bp insert: PCR1_rear (3') is its reverse
        // complement, so the same span also carries a 3'-tagged hit, and both
        // end zones cover it on a 162 bp read.
        let mut short = b"ACTTGCCTGTCGCTCTATCTTC".to_vec();
        short.extend_from_slice(b"GGGG");
        short.extend(splitmix_dna(0, 136));
        let mut long = short.clone();
        long.extend(splitmix_dna(5, 38));
        for split in [true, false] {
            let c = AdapterConfig {
                adapters: preset_ont(),
                error_rate: 0.2,
                end_size: 150,
                split,
                candidate_index: std::sync::OnceLock::new(),
            };
            assert_insert_kept(&adapter_segments(&short, &c), 162, 22, 3, true);
            assert_insert_kept(&adapter_segments(&long, &c), 200, 22, 3, true);
        }
    }

    /// Mirror: `cDNA_front` and `PCS110_front` end in `PCR2_front`, so their
    /// reverse complements are `PCR2_rear` with three or four leading bases.
    #[test]
    fn short_read_with_rear_adapter_keeps_insert_under_ont_preset() {
        let mut short = splitmix_dna(0, 136);
        short.extend_from_slice(b"GGGG");
        short.extend_from_slice(b"GCAATATCAGCACCAACAGAAA");
        let mut long = splitmix_dna(5, 38);
        long.extend_from_slice(&short);
        for split in [true, false] {
            let c = AdapterConfig {
                adapters: preset_ont(),
                error_rate: 0.2,
                end_size: 150,
                split,
                candidate_index: std::sync::OnceLock::new(),
            };
            assert_insert_kept(&adapter_segments(&short, &c), 162, 22, 4, false);
            assert_insert_kept(&adapter_segments(&long, &c), 200, 22, 4, false);
        }
    }

    /// A minimal user FASTA: `f` tagged 5' and its reverse complement tagged
    /// 3', with the adapter and its insert on a 162 bp read. Equal lengths take
    /// the batched path; the shortened rear entry takes the singleton path.
    #[test]
    fn front_and_reverse_complement_rear_pair_keep_insert_on_short_read() {
        let f = b"ACTTGCCTGTCGCTCTATCTTC";
        let r = reverse_complement(f);
        let mut w = f.to_vec();
        w.extend(splitmix_dna(3, 140));
        for rear in [r.as_slice(), &r[1..]] {
            for split in [true, false] {
                let c = AdapterConfig {
                    adapters: vec![ad("f", f, End::Five), ad("r", rear, End::Three)],
                    error_rate: 0.2,
                    end_size: 150,
                    split,
                    candidate_index: std::sync::OnceLock::new(),
                };
                assert_eq!(
                    adapter_segments(&w, &c),
                    vec![(22, 162)],
                    "Split mode {split}, rear length {}",
                    rear.len()
                );
            }
        }
    }

    /// Inside the overlap of both end zones the outcome is geometric and the
    /// tag is ignored; with one zone only, the tag decides.
    #[test]
    fn classify_terminal_is_geometric_in_the_overlap() {
        // Both zones cover every hit on a 60 bp window with `end_size` 60.
        for tag in [End::Five, End::Three, End::Both] {
            assert_eq!(classify_terminal(0, 20, 60, 60, tag), Terminal::Five);
            assert_eq!(classify_terminal(40, 60, 60, 60, tag), Terminal::Three);
            assert_eq!(classify_terminal(20, 40, 60, 60, tag), Terminal::Excise);
            assert_eq!(classify_terminal(11, 31, 60, 60, tag), Terminal::Five);
            assert_eq!(classify_terminal(12, 32, 60, 60, tag), Terminal::Excise);
            assert_eq!(classify_terminal(29, 49, 60, 60, tag), Terminal::Three);
            // Both flanks within slack: the nearer end.
            assert_eq!(classify_terminal(8, 30, 35, 35, tag), Terminal::Three);
            assert_eq!(classify_terminal(4, 30, 35, 35, tag), Terminal::Five);
        }
        // One zone only: the tag decides.
        assert_eq!(
            classify_terminal(0, 20, 400, 150, End::Five),
            Terminal::Five
        );
        assert_eq!(
            classify_terminal(0, 20, 400, 150, End::Three),
            Terminal::None
        );
        assert_eq!(
            classify_terminal(380, 400, 400, 150, End::Three),
            Terminal::Three
        );
        assert_eq!(
            classify_terminal(380, 400, 400, 150, End::Five),
            Terminal::None
        );
        assert_eq!(
            classify_terminal(200, 220, 400, 150, End::Both),
            Terminal::None
        );
    }

    /// An interior adapter splits the read into its two flanks.
    #[test]
    fn splits_on_interior_adapter() {
        let adapter = b"GGGGTTTTGGGGTTTT"; // 16 bp, no C/A so it cannot match the flanks
        let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec(); // 24 bp lead (> end_size 20)
        let cut_start = w.len();
        w.extend_from_slice(adapter);
        w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC"); // 24 bp tail
        let c = cfg(vec![ad("mid", adapter, End::Both)], true);
        let segs = adapter_segments(&w, &c);
        assert_eq!(segs.len(), 2, "Interior adapter splits the read");
        assert_eq!(segs[0], (0, cut_start));
        assert_eq!(segs[1], (cut_start + adapter.len(), w.len()));
    }

    /// Ends-only mode leaves an interior adapter in place.
    #[test]
    fn ends_only_suppresses_interior_split() {
        let adapter = b"GGGGTTTTGGGGTTTT";
        let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
        w.extend_from_slice(adapter);
        w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC");
        let c = cfg(vec![ad("mid", adapter, End::Both)], false); // ends-only
        assert_eq!(adapter_segments(&w, &c), vec![(0, w.len())]);
    }

    /// A 5' adapter, an insert and a 3' adapter in ends-only mode: both ends
    /// trim to the insert although only the two end zones are searched.
    #[test]
    fn ends_only_trims_both_terminal_adapters() {
        let adapter5 = b"ACGTACGTACGT"; // 12 bp
        let adapter3 = b"TTTTGGGGCCCC"; // 12 bp, distinct from adapter5
        let insert = b"AAAAAAAAAAAA"; // 12 bp
        let mut w = adapter5.to_vec();
        w.extend_from_slice(insert);
        w.extend_from_slice(adapter3);
        let c = cfg(
            vec![
                ad("five", adapter5, End::Five),
                ad("three", adapter3, End::Three),
            ],
            false, // ends-only
        );
        assert_eq!(adapter_segments(&w, &c), vec![(12, 24)]);
    }

    /// A terminal 5' adapter that starts inside `end_size` but ends beyond it:
    /// with `end_size` 4, a 12 bp adapter at position 2 spans [2, 14). A head
    /// zone of `window[..end_size]` (4 bytes) cannot contain a 12-byte match;
    /// the `end_size + len` sizing gives `window[..16]`, which does.
    #[test]
    fn ends_only_trims_adapter_straddling_end_size() {
        let adapter = b"ACGTACGTACGT"; // 12 bp
        let mut w = b"AA".to_vec(); // 2 bp prefix -> adapter starts at position 2
        w.extend_from_slice(adapter); // adapter occupies [2, 14)
        w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCC"); // 20 bp tail
        let c = AdapterConfig {
            adapters: vec![ad("five", adapter, End::Five)],
            error_rate: 0.2,
            end_size: 4,
            split: false, // ends-only
            candidate_index: std::sync::OnceLock::new(),
        };
        let segs = adapter_segments(&w, &c);
        assert_eq!(segs, vec![(14, w.len())]);
    }

    /// A pattern below `MIN_PATTERN_LEN` is never searched.
    #[test]
    fn short_pattern_is_skipped() {
        let short = b"GGTGCTG"; // 7 bp < MIN_PATTERN_LEN
        let w = b"GGTGCTGAAAAAAAAAAAAAAAA";
        let c = cfg(vec![ad("flank", short, End::Five)], true);
        assert_eq!(adapter_segments(w, &c), vec![(0, w.len())]);
    }

    /// An empty window yields no segments.
    #[test]
    fn empty_window_returns_empty() {
        let c = cfg(vec![ad("a", b"ACGTACGTACGT", End::Both)], true);
        assert_eq!(adapter_segments(b"", &c), vec![]);
    }

    /// The window is the adapter, matched `End::Both` at the start: the 5'
    /// trim advances `lo` to `n`, so `lo >= hi` and the whole window is
    /// consumed.
    #[test]
    fn whole_window_consumed_returns_empty() {
        let adapter = b"ACGTACGTACGT"; // 12 bp
        let c = cfg(vec![ad("a", adapter, End::Both)], true);
        assert_eq!(adapter_segments(adapter, &c), vec![]);
    }

    /// Mirror of `trims_5prime_adapter_and_outboard` with the adapter at the
    /// 3' end: insert first, adapter last.
    #[test]
    fn trims_3prime_adapter() {
        let adapter = b"ACGTACGTACGT"; // 12 bp
        let mut w = b"AAAAAAAAAAAA".to_vec();
        w.extend_from_slice(adapter);
        let c = cfg(vec![ad("a", adapter, End::Three)], false);
        assert_eq!(adapter_segments(&w, &c), vec![(0, 12)]);
    }

    /// Two distinct interior adapters whose hits overlap by 6 bp: `a` matches
    /// [24, 40) and `b` matches [34, 50), constructed so their shared 6 bp
    /// region ("TGTGTG", the tail of `a` and head of `b`) is the same window
    /// bytes, giving both an exact (cost 0) hit. The overlap merges into one
    /// excision, leaving exactly 2 segments.
    #[test]
    fn overlapping_interior_cuts_merge() {
        let a = b"GGGGTTTTTGTGTGTG"; // 16 bp
        let b = b"TGTGTGTGTTTTGGGG"; // 16 bp, shares a's last 6 bp as its first 6 bp
        let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec(); // 24 bp lead
        w.extend_from_slice(a); // a occupies [24, 40)
        w.extend_from_slice(&b[6..]); // appends b's non-overlapping tail; b occupies [34, 50)
        w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC"); // 24 bp tail
        let c = cfg(vec![ad("a", a, End::Both), ad("b", b, End::Both)], true);
        let segs = adapter_segments(&w, &c);
        assert_eq!(
            segs.len(),
            2,
            "Overlapping interior cuts merge into one excision"
        );
        assert_eq!(segs[0], (0, 24));
        assert_eq!(segs[1], (50, w.len()));
    }

    /// The terminal hit [0, 16) overlaps the interior hit [10, 26); clipping the
    /// interior interval to the keep window still excises [16, 26).
    #[test]
    fn straddling_cut_is_clipped_not_leaked() {
        let t_prefix = b"GGTGTGGTTT"; // 10 bp
        let overlap = b"GTTGGT"; // 6 bp, shared
        let s_suffix = b"TGGTGTTGGG"; // 10 bp
        let mut t = t_prefix.to_vec();
        t.extend_from_slice(overlap); // t = 16 bp, occupies [0, 16)
        let mut s = overlap.to_vec();
        s.extend_from_slice(s_suffix); // s = 16 bp, occupies [10, 26)

        let mut w = t.clone();
        w.extend_from_slice(s_suffix); // s's non-overlapping tail; s occupies [10, 26)
        w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"); // 34 bp tail -> n=60

        let c = AdapterConfig {
            adapters: vec![ad("t", &t, End::Five), ad("s", &s, End::Both)],
            error_rate: 0.2,
            end_size: 9,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let segs = adapter_segments(&w, &c);
        assert_eq!(
            segs,
            vec![(26, 60)],
            "Adapter s is fully excised, no leaked bases before 26"
        );
        // No retained segment may contain the interior adapter.
        for &(seg_start, seg_end) in &segs {
            assert!(
                !w[seg_start..seg_end]
                    .windows(s.len())
                    .any(|win| win == s.as_slice())
            );
        }
    }

    /// A cost-4 hit within `k_end` of 6 but above the interior `k_mid` of 3
    /// does not split the read.
    #[test]
    fn interior_above_k_mid_does_not_split() {
        let adapter = b"GGTTGGTTGGTT"; // 12 bp
        let mut mutated = adapter.to_vec();
        for &i in &[1usize, 4, 7, 10] {
            mutated[i] = match mutated[i] {
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            };
        }
        let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec(); // 24 bp lead
        w.extend_from_slice(&mutated); // interior copy at [24, 36), cost 4 vs `adapter`
        w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC"); // 24 bp tail
        // `end_size` of 10 keeps the intended hit and sassy's other fuzzy hits
        // (found under the wide `k_end` of 6) away from the near-end terminal
        // checks.
        let c = AdapterConfig {
            adapters: vec![ad("mid", adapter, End::Both)],
            error_rate: 0.5,
            end_size: 10,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        assert_eq!(
            adapter_segments(&w, &c),
            vec![(0, w.len())],
            "Cost 4 hit is above k_mid=3 and must not split the read"
        );
    }

    /// A six-base insertion expands the terminal alignment from 20 to 26 bases.
    /// The terminal search window includes `k_end` additional bases, so
    /// ends-only and split modes select the same [2, 28) alignment.
    #[test]
    fn ends_only_equals_split_on_indel_terminal_adapter() {
        let adapter = b"AAAACCCCGGGGTTTTACGT"; // 20 bp
        let extra = b"CTGACT"; // 6 bp splice, foreign bases -> forces insertion
        let mut copy = adapter[..10].to_vec();
        copy.extend_from_slice(extra);
        copy.extend_from_slice(&adapter[10..]); // copy = 26 bp

        let mut w = b"AA".to_vec(); // 2 bp prefix -> copy occupies [2, 28)
        w.extend_from_slice(&copy);
        w.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"); // 30 bp clean insert tail

        let c_split = AdapterConfig {
            adapters: vec![ad("five", adapter, End::Five)],
            error_rate: 0.3,
            end_size: 4,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let c_ends_only = AdapterConfig {
            split: false,
            ..c_split.clone()
        };

        let split_segs = adapter_segments(&w, &c_split);
        let ends_only_segs = adapter_segments(&w, &c_ends_only);

        assert_eq!(
            split_segs,
            vec![(28, w.len())],
            "Split mode finds the full 26bp indel-bearing hit and trims to 28"
        );
        assert_eq!(
            ends_only_segs, split_segs,
            "Ends-only must match split mode exactly: the end zone must be wide \
             enough (end_size + len + k_end) to contain the full indel-lengthened hit"
        );
        // Adapter bases are removed, rather than both results being equal and empty.
        assert_eq!(ends_only_segs[0].0, 28);
    }

    /// A 40 bp insert plus a 20 bp adapter at the 3' end, tagged `End::Both`,
    /// with `end_size >= n` so both zones overlap. The insert [0, 40) is kept
    /// and the read is not dropped.
    #[test]
    fn three_prime_both_adapter_on_short_read_trims_tail_not_whole_read() {
        let adapter = b"GGGGTTTTGGGGTTTTGGGG"; // 20bp, G/T only (no A/C to collide with insert)
        let mut w = vec![b'A'; 40];
        w.extend_from_slice(adapter);
        let split = AdapterConfig {
            adapters: vec![ad("a", adapter, End::Both)],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let ends = AdapterConfig {
            split: false,
            ..split.clone()
        };
        assert_eq!(adapter_segments(&w, &split), vec![(0, 40)], "Split mode");
        assert_eq!(adapter_segments(&w, &ends), vec![(0, 40)], "Ends-only mode");
    }

    /// A 5' `End::Both` adapter on a short read trims the head and keeps the
    /// insert.
    #[test]
    fn five_prime_both_adapter_on_short_read_trims_head() {
        let adapter = b"GGGGTTTTGGGGTTTTGGGG";
        let mut w = adapter.to_vec();
        w.extend_from_slice(&[b'A'; 40]); // adapter [0,20) + 40bp insert
        let split = AdapterConfig {
            adapters: vec![ad("a", adapter, End::Both)],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        assert_eq!(adapter_segments(&w, &split), vec![(20, 60)]);
    }

    /// The two terminal patterns and their reverse complements are distinct,
    /// leaving the 40-base insert as the only retained segment.
    #[test]
    fn both_adapters_at_both_ends_keep_middle() {
        let a5 = b"GGGGTTTTGGGGTTTTGGGG";
        let a3 = b"AAAAGGGGAAAAGGGGAAAA"; // A/G only (purine): not self-complementary, not revcomp(a5)
        let mut w = a5.to_vec();
        w.extend_from_slice(&[b'T'; 40]); // insert bytes do not match either adapter's revcomp
        w.extend_from_slice(a3);
        let cfg = AdapterConfig {
            adapters: vec![ad("a5", a5, End::Both), ad("a3", a3, End::Both)],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let segs = adapter_segments(&w, &cfg);
        assert_eq!(segs, vec![(20, 60)]);
    }

    /// Distinct end-specific adapters preserve the insert when the end search
    /// regions overlap on a short read.
    #[test]
    fn inferred_single_end_adapters_on_short_read_keep_insert() {
        let a5 = b"GGGGTTTTGGGGTTTTGGGG"; // 20bp, G/T only
        let a3 = b"AAAAGGGGAAAAGGGGAAAA"; // 20bp, A/G only: not a5, not revcomp(a5)
        let mut w = a5.to_vec();
        w.extend_from_slice(&[b'C'; 40]); // 40bp insert, no match to either adapter/revcomp
        w.extend_from_slice(a3);
        let c = AdapterConfig {
            adapters: vec![ad("five", a5, End::Five), ad("three", a3, End::Three)],
            error_rate: 0.2,
            end_size: 150, // >= n, zones overlap
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        assert_eq!(adapter_segments(&w, &c), vec![(20, 60)], "Insert survives");
    }
}
