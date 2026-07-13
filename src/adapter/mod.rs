pub mod detect;
pub mod infer;
pub mod preset;
pub mod search;

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::OnceLock;

use aho_corasick::AhoCorasick;
use search::{
    BatchedDnaSearcher, DnaSearcher, hits, new_batched_searcher, new_searcher, pattern_hits,
};

thread_local! {
    /// One reverse-complement searcher per thread, reused across reads so
    /// `adapter_segments` doesn't allocate a fresh searcher (and its scratch
    /// buffers) on every call. Per-thread, so the parallel workflows stay
    /// data-race-free without sharing.
    static RC_SEARCHER: RefCell<DnaSearcher> = RefCell::new(new_searcher());

    /// Pattern-batched searcher. Sassy's `Dna` profile does not implement
    /// `search_patterns`; on A/C/G/T input its IUPAC profile is equivalent.
    static BATCH_SEARCHER: RefCell<BatchedDnaSearcher> = RefCell::new(new_batched_searcher());
}

/// Which read end a catalog sequence is expected at — this gates TERMINAL
/// trimming only: `Five` trims at the 5' end, `Three` at the 3' end, `Both` at
/// either. Interior chimera-splitting (when enabled) considers any adapter that
/// matches in the read interior regardless of this tag, since a front/rear
/// adapter appearing mid-read is itself the chimera signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum End {
    Five,
    Three,
    Both,
}

/// One searchable adapter/primer/barcode/flank sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Adapter {
    pub name: String,
    pub seq: Vec<u8>,
    pub end: End,
}

/// Resolved adapter-trimming settings for a run.
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    pub adapters: Vec<Adapter>,
    /// End-match tolerance as a fraction of adapter length (`k_end`).
    pub error_rate: f64,
    /// Bases at each end classified as "terminal" (trim) vs interior (split).
    pub end_size: usize,
    /// Split on interior adapters. False = ends-only (`--adapter-ends-only`).
    pub split: bool,
    /// Exact-seed automaton for lossless whole-read candidate filtering. Built
    /// lazily after presence detection/inference has finalized `adapters`.
    pub(crate) candidate_index: OnceLock<CandidateIndex>,
}

impl AdapterConfig {
    pub(crate) fn replace_adapters(&mut self, adapters: Vec<Adapter>) {
        self.adapters = adapters;
        self.candidate_index = OnceLock::new();
    }
}

#[derive(Debug, Clone)]
pub(crate) struct CandidateIndex {
    matcher: Option<AhoCorasick>,
    seed_adapters: Vec<Vec<usize>>,
    adapter_lens: Vec<usize>,
    error_rate: f64,
    terminal_batches: Vec<TerminalBatch>,
    batched_adapters: Vec<bool>,
}

/// Equal-length adapters searched together through Sassy's pattern-parallel
/// API. Sassy handles forward and reverse-complement matching for the batch.
#[derive(Debug, Clone)]
struct TerminalBatch {
    adapter_indices: Vec<usize>,
    patterns: Vec<Vec<u8>>,
    len: usize,
    k_end: usize,
}

impl CandidateIndex {
    fn new(adapters: &[Adapter], error_rate: f64, include_interior: bool) -> Self {
        let (matcher, seed_adapters) = if include_interior {
            let mut seeds: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
            for (adapter_idx, adapter) in adapters.iter().enumerate() {
                let len = adapter.seq.len();
                if len < MIN_PATTERN_LEN {
                    continue;
                }
                let k_mid = (0.5 * error_rate * len as f64).floor() as usize;
                add_partition_seeds(&mut seeds, adapter_idx, &adapter.seq, k_mid);
                let rc = reverse_complement(&adapter.seq);
                add_partition_seeds(&mut seeds, adapter_idx, &rc, k_mid);
            }

            let patterns: Vec<Vec<u8>> = seeds.keys().cloned().collect();
            let seed_adapters: Vec<Vec<usize>> = seeds.into_values().collect();
            let matcher = (!patterns.is_empty()).then(|| {
                AhoCorasick::builder()
                    .ascii_case_insensitive(true)
                    .build(&patterns)
                    .expect("adapter seeds are nonempty ASCII DNA patterns")
            });
            (matcher, seed_adapters)
        } else {
            (None, Vec::new())
        };

        // Sassy can pack equal-length patterns across SIMD lanes. Keep
        // singletons on the ordinary search path: batching them cannot use
        // pattern-level parallelism and is slower for these terminal windows.
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
                terminal_batches.push(TerminalBatch {
                    adapter_indices,
                    patterns: batch_patterns,
                    len,
                    k_end: (error_rate * len as f64).floor() as usize,
                });
            }
        }
        Self {
            matcher,
            seed_adapters,
            adapter_lens: adapters.iter().map(|adapter| adapter.seq.len()).collect(),
            error_rate,
            terminal_batches,
            batched_adapters,
        }
    }

    fn candidate_windows(&self, text: &[u8], adapter_count: usize) -> Vec<Vec<(usize, usize)>> {
        let mut windows = vec![Vec::new(); adapter_count];
        let Some(matcher) = &self.matcher else {
            return windows;
        };
        for m in matcher.find_overlapping_iter(text) {
            for &adapter_idx in &self.seed_adapters[m.pattern().as_usize()] {
                let len = self.adapter_lens[adapter_idx];
                let k_end = (self.error_rate * len as f64).floor() as usize;
                // The exact seed lies inside the <=k_mid alignment. A radius of
                // pattern length + k_end on each side necessarily contains that
                // entire alignment and enough context for the original k_end
                // Sassy search to reproduce its span/tie behavior.
                let radius = len + k_end;
                windows[adapter_idx].push((
                    m.start().saturating_sub(radius),
                    m.end().saturating_add(radius).min(text.len()),
                ));
            }
        }
        for adapter_windows in &mut windows {
            adapter_windows.sort_unstable();
            let mut merged: Vec<(usize, usize)> = Vec::with_capacity(adapter_windows.len());
            for &(start, end) in adapter_windows.iter() {
                if let Some(last) = merged.last_mut()
                    && start <= last.1
                {
                    last.1 = last.1.max(end);
                } else {
                    merged.push((start, end));
                }
            }
            *adapter_windows = merged;
        }
        windows
    }
}

fn add_partition_seeds(
    seeds: &mut BTreeMap<Vec<u8>, Vec<usize>>,
    adapter_idx: usize,
    pattern: &[u8],
    max_edits: usize,
) {
    let parts = (max_edits + 1).min(pattern.len());
    for i in 0..parts {
        let start = i * pattern.len() / parts;
        let end = (i + 1) * pattern.len() / parts;
        let owners = seeds.entry(pattern[start..end].to_vec()).or_default();
        if owners.last() != Some(&adapter_idx) {
            owners.push(adapter_idx);
        }
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'a' => b't',
            b'C' => b'G',
            b'c' => b'g',
            b'G' => b'C',
            b'g' => b'c',
            b'T' => b'A',
            b't' => b'a',
            _ => b,
        })
        .collect()
}

/// Sequences shorter than this are never searched standalone — a <11 bp pattern
/// matches almost anywhere under any error budget. The 7 bp catalog flanks are
/// construction anchors, not standalone patterns.
pub const MIN_PATTERN_LEN: usize = 11;

/// Terminal classification of a hit: which end (if any) it trims.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Terminal {
    Five,
    Three,
    /// Eligible for BOTH ends (short read, overlapping end-zones): excise the
    /// adapter span and keep BOTH flanks. See `classify_terminal`.
    Excise,
    None,
}

/// Classify a hit at window coords `[start, end)` in a length-`n` window.
///
/// A hit eligible for BOTH ends means the read is short enough that the two
/// end-zones overlap (`n <= 2*end_size`), so the adapter sits within `end_size`
/// of each end. Trimming toward the nearer end here would delete the whole
/// outboard arm — which, for a central chimera-junction adapter, is real
/// insert. Instead classify it as `Excise`: cut out just the adapter and keep
/// both flanks. For a genuinely terminal adapter (abutting an end) the outboard
/// flank is empty, so excising is identical to trimming that end; only a central
/// adapter — which cannot be a terminal leader — is actually split.
fn classify_terminal(start: usize, end: usize, n: usize, end_size: usize, tag: End) -> Terminal {
    let near5 = start <= end_size && matches!(tag, End::Five | End::Both);
    let near3 = end >= n.saturating_sub(end_size) && matches!(tag, End::Three | End::Both);
    match (near5, near3) {
        (true, true) => Terminal::Excise,
        (true, false) => Terminal::Five,
        (false, true) => Terminal::Three,
        (false, false) => Terminal::None,
    }
}

/// Ends-only variant of `classify_terminal`: splitting is disabled, so an
/// `Excise` hit can't keep both flanks. Resolve it back to a terminal trim
/// toward the nearer end (`start` vs `n - end`), leaving every other outcome
/// untouched.
fn ends_only_terminal(start: usize, end: usize, n: usize, end_size: usize, tag: End) -> Terminal {
    match classify_terminal(start, end, n, end_size, tag) {
        Terminal::Excise => {
            if start <= n - end {
                Terminal::Five
            } else {
                Terminal::Three
            }
        },
        other => other,
    }
}

/// Compute adapter keep-segments for `window`:
///   - terminal hits within `end_size` of an end trim that end inward;
///   - interior hits (stricter `k_mid`) excise and split.
///
/// When `cfg.split` is false (`--adapter-ends-only`), only the two end-zones
/// are searched at all — the interior is never scanned, since no hit found
/// there could ever be actioned.
///
/// Returns `[start,end)` spans in `window` coordinates.
pub fn adapter_segments(window: &[u8], cfg: &AdapterConfig) -> Vec<(usize, usize)> {
    let n = window.len();
    if n == 0 {
        return vec![];
    }
    if cfg.adapters.is_empty() {
        return vec![(0, n)];
    }
    let end_size = cfg.end_size.min(n);

    let mut lo = 0usize; // 5' keep-boundary (advances inward on terminal hits)
    let mut hi = n; // 3' keep-boundary (retreats inward on terminal hits)
    let mut interior: Vec<(usize, usize)> = Vec::new();
    let search_index = cfg
        .candidate_index
        .get_or_init(|| CandidateIndex::new(&cfg.adapters, cfg.error_rate, cfg.split));
    let interior_windows = cfg
        .split
        .then(|| search_index.candidate_windows(window, cfg.adapters.len()));
    // `search_patterns` is only available through Sassy's IUPAC profile. It is
    // exactly equivalent to the existing DNA path on A/C/G/T text; preserve
    // prior behavior when either terminal search region has another symbol.
    // Only inspect the bounded terminal regions, not the potentially huge
    // interior of a long read that no terminal search will touch.
    let batch_span = search_index
        .terminal_batches
        .iter()
        .map(|batch| end_size + batch.len + batch.k_end)
        .max()
        .unwrap_or(0)
        .min(n);
    let use_batches = window[..batch_span]
        .iter()
        .chain(&window[n - batch_span..])
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'));

    RC_SEARCHER.with_borrow_mut(|searcher| {
        // All adapters in a batch have the same length and edit budget, so the
        // head/tail windows are shared. This collapses the ONT catalog's 96
        // equal-length barcode searches into one SIMD pattern search per end.
        if use_batches {
            BATCH_SEARCHER.with_borrow_mut(|batch_searcher| {
                for batch in &search_index.terminal_batches {
                    let head_end = (end_size + batch.len + batch.k_end).min(n);
                    for h in pattern_hits(
                        batch_searcher,
                        &batch.patterns,
                        &window[..head_end],
                        batch.k_end,
                    ) {
                        let ad = &cfg.adapters[batch.adapter_indices[h.pattern_idx]];
                        let terminal = if cfg.split {
                            classify_terminal(h.text_start, h.text_end, n, end_size, ad.end)
                        } else {
                            ends_only_terminal(h.text_start, h.text_end, n, end_size, ad.end)
                        };
                        match terminal {
                            Terminal::Five => lo = lo.max(h.text_end),
                            Terminal::Excise => interior.push((h.text_start, h.text_end)),
                            Terminal::Three | Terminal::None => {},
                        }
                    }

                    let tail_start = n.saturating_sub(end_size + batch.len + batch.k_end);
                    for h in pattern_hits(
                        batch_searcher,
                        &batch.patterns,
                        &window[tail_start..],
                        batch.k_end,
                    ) {
                        let ad = &cfg.adapters[batch.adapter_indices[h.pattern_idx]];
                        let (start, end) = (tail_start + h.text_start, tail_start + h.text_end);
                        let terminal = if cfg.split {
                            classify_terminal(start, end, n, end_size, ad.end)
                        } else {
                            ends_only_terminal(start, end, n, end_size, ad.end)
                        };
                        match terminal {
                            Terminal::Three => hi = hi.min(start),
                            Terminal::Excise => interior.push((start, end)),
                            Terminal::Five | Terminal::None => {},
                        }
                    }
                }
            });
        }

        // Singleton length groups retain the original one-pattern search. On a
        // non-ACGT window every adapter takes this path for exact compatibility.
        for (adapter_idx, ad) in cfg.adapters.iter().enumerate() {
            let len = ad.seq.len();
            if len < MIN_PATTERN_LEN || (use_batches && search_index.batched_adapters[adapter_idx])
            {
                continue;
            }
            let k_end = (cfg.error_rate * len as f64).floor() as usize;
            let head_end = (end_size + len + k_end).min(n);
            for h in hits(searcher, &ad.seq, &window[..head_end], k_end) {
                if cfg.split {
                    match classify_terminal(h.start, h.end, n, end_size, ad.end) {
                        Terminal::Five => lo = lo.max(h.end),
                        Terminal::Excise => interior.push((h.start, h.end)),
                        Terminal::Three | Terminal::None => {},
                    }
                } else if ends_only_terminal(h.start, h.end, n, end_size, ad.end) == Terminal::Five
                {
                    lo = lo.max(h.end);
                }
            }
            let tail_start = n.saturating_sub(end_size + len + k_end);
            for h in hits(searcher, &ad.seq, &window[tail_start..], k_end) {
                let (s, e) = (tail_start + h.start, tail_start + h.end);
                if cfg.split {
                    match classify_terminal(s, e, n, end_size, ad.end) {
                        Terminal::Three => hi = hi.min(s),
                        Terminal::Excise => interior.push((s, e)),
                        Terminal::Five | Terminal::None => {},
                    }
                } else if ends_only_terminal(s, e, n, end_size, ad.end) == Terminal::Three {
                    hi = hi.min(s);
                }
            }
        }

        // Exact partition seeds identify every possible interior match.
        // Interior hits are accepted only up to `k_mid`, so Sassy searches the
        // candidate window at that limit rather than the looser end budget.
        if let Some(candidate_windows) = &interior_windows {
            for (adapter_idx, ad) in cfg.adapters.iter().enumerate() {
                if ad.seq.len() < MIN_PATTERN_LEN {
                    continue;
                }
                let k_mid = (0.5 * cfg.error_rate * ad.seq.len() as f64).floor() as usize;
                for &(candidate_start, candidate_end) in &candidate_windows[adapter_idx] {
                    for h in hits(
                        searcher,
                        &ad.seq,
                        &window[candidate_start..candidate_end],
                        k_mid,
                    ) {
                        let start = candidate_start + h.start;
                        let end = candidate_start + h.end;
                        if classify_terminal(start, end, n, end_size, ad.end) == Terminal::None {
                            interior.push((start, end));
                        }
                    }
                }
            }
        }
    });

    if lo >= hi {
        return vec![]; // whole window consumed by terminal adapters
    }

    // Merge interior cuts strictly inside (lo, hi), then carve gaps.
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

#[cfg(test)]
mod segment_tests {
    use super::*;

    fn cfg(adapters: Vec<Adapter>, split: bool) -> AdapterConfig {
        AdapterConfig {
            adapters,
            error_rate: 0.2,
            end_size: 20,
            split,
            candidate_index: std::sync::OnceLock::new(),
        }
    }
    fn ad(name: &str, seq: &[u8], end: End) -> Adapter {
        Adapter {
            name: name.into(),
            seq: seq.to_vec(),
            end,
        }
    }

    /// Exhaustive reference implementation for candidate-filter comparisons.
    fn reference_segments(window: &[u8], cfg: &AdapterConfig) -> Vec<(usize, usize)> {
        let n = window.len();
        if n == 0 {
            return vec![];
        }
        if cfg.adapters.is_empty() {
            return vec![(0, n)];
        }
        let end_size = cfg.end_size.min(n);
        let mut lo = 0usize;
        let mut hi = n;
        let mut interior = Vec::new();
        let mut searcher = new_searcher();

        for adapter in &cfg.adapters {
            let len = adapter.seq.len();
            if len < MIN_PATTERN_LEN {
                continue;
            }
            let k_end = (cfg.error_rate * len as f64).floor() as usize;
            let k_mid = (0.5 * cfg.error_rate * len as f64).floor() as usize;
            if cfg.split {
                for hit in hits(&mut searcher, &adapter.seq, window, k_end) {
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
                for hit in hits(&mut searcher, &adapter.seq, &window[..head_end], k_end) {
                    if ends_only_terminal(hit.start, hit.end, n, end_size, adapter.end)
                        == Terminal::Five
                    {
                        lo = lo.max(hit.end);
                    }
                }
                let tail_start = n.saturating_sub(end_size + len + k_end);
                for hit in hits(&mut searcher, &adapter.seq, &window[tail_start..], k_end) {
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

    #[test]
    fn candidate_search_matches_full_search_randomized() {
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

        let mut rng = Lcg(0x4e4f_4f44_4c45_5301);
        for case in 0..400 {
            let adapters: Vec<Adapter> = (0..(1 + rng.below(10)))
                .map(|i| {
                    let len = 11 + rng.below(40);
                    Adapter {
                        name: format!("a{i}"),
                        seq: rng.dna(len),
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

            // Plant one adapter with up to k_mid substitutions at an end or in
            // the interior. Random background also exercises false candidates.
            let planted = rng.below(adapters.len());
            let pattern = &adapters[planted].seq;
            if pattern.len() <= window.len() {
                let max_edits = (0.1 * pattern.len() as f64).floor() as usize;
                let mut copy = pattern.clone();
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
                "candidate/reference mismatch in randomized case {case}"
            );
        }
    }

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
                !index.candidate_windows(&text, 1)[0].is_empty(),
                "lossless seed filter rejected <=k edit case {case}"
            );
        }
    }

    #[test]
    fn no_adapters_is_identity() {
        let w = b"ACGTACGTACGTACGT";
        assert_eq!(adapter_segments(w, &cfg(vec![], true)), vec![(0, w.len())]);
    }

    #[test]
    fn trims_5prime_adapter_and_outboard() {
        let adapter = b"ACGTACGTACGT"; // 12 bp
        let mut w = adapter.to_vec();
        w.extend_from_slice(b"AAAAAAAAAAAA");
        let c = cfg(vec![ad("a", adapter, End::Five)], false);
        assert_eq!(adapter_segments(&w, &c), vec![(12, 24)]);
    }

    #[test]
    fn central_chimera_on_short_read_splits_both_arms() {
        // With the default end_size=150, both end-zones overlap for any read
        // <= 2*end_size (300bp). A chimera-junction adapter sitting within
        // end_size of BOTH ends must SPLIT the read (keep both inserts), not be
        // treated as a terminal adapter — which discarded the entire outboard
        // arm (up to ~end_size bases of real insert).
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
            "central chimera must split into both arms, not lose insert1"
        );
    }

    #[test]
    fn splits_on_interior_adapter() {
        let adapter = b"GGGGTTTTGGGGTTTT"; // 16 bp, no C/A so it can't match the flanks
        let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec(); // 24 bp lead (> end_size 20)
        let cut_start = w.len();
        w.extend_from_slice(adapter);
        w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC"); // 24 bp tail
        let c = cfg(vec![ad("mid", adapter, End::Both)], true);
        let segs = adapter_segments(&w, &c);
        assert_eq!(segs.len(), 2, "interior adapter splits the read");
        assert_eq!(segs[0], (0, cut_start));
        assert_eq!(segs[1], (cut_start + adapter.len(), w.len()));
    }

    #[test]
    fn ends_only_suppresses_interior_split() {
        let adapter = b"GGGGTTTTGGGGTTTT";
        let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
        w.extend_from_slice(adapter);
        w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC");
        let c = cfg(vec![ad("mid", adapter, End::Both)], false); // ends-only
        assert_eq!(adapter_segments(&w, &c), vec![(0, w.len())]);
    }

    #[test]
    fn ends_only_trims_both_terminal_adapters() {
        // [5' adapter][insert][3' adapter], ends-only mode: both ends must
        // still trim to the insert even though only the two end-zones (not
        // the whole window) are searched.
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

    #[test]
    fn ends_only_trims_adapter_straddling_end_size() {
        // A terminal 5' adapter whose match STARTS within end_size but ENDS
        // beyond it: end_size=4, a 12bp adapter starting at position 2, so
        // it spans [2,14) crossing the end_size=4 boundary. If the ends-only
        // head zone were naively sized as `window[..end_size]` (4 bytes),
        // this adapter (needing ~12 bytes of text) could never be found, and
        // the read would come back untrimmed. With the correct
        // `end_size + len` zone sizing, the head zone is `window[..16]`,
        // which fully contains the match.
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

    #[test]
    fn short_pattern_is_skipped() {
        let short = b"GGTGCTG"; // 7 bp < MIN_PATTERN_LEN
        let w = b"GGTGCTGAAAAAAAAAAAAAAAA";
        let c = cfg(vec![ad("flank", short, End::Five)], true);
        assert_eq!(adapter_segments(w, &c), vec![(0, w.len())]);
    }

    #[test]
    fn empty_window_returns_empty() {
        let c = cfg(vec![ad("a", b"ACGTACGTACGT", End::Both)], true);
        assert_eq!(adapter_segments(b"", &c), vec![]);
    }

    #[test]
    fn whole_window_consumed_returns_empty() {
        // Window IS the adapter, matched End::Both at the very start: the terminal-5'
        // branch advances `lo` to `n`, so `lo >= hi` and the whole window is consumed.
        let adapter = b"ACGTACGTACGT"; // 12 bp
        let c = cfg(vec![ad("a", adapter, End::Both)], true);
        assert_eq!(adapter_segments(adapter, &c), vec![]);
    }

    #[test]
    fn trims_3prime_adapter() {
        // Mirror of trims_5prime_adapter_and_outboard, but the adapter sits at the
        // 3' end: insert first, adapter last.
        let adapter = b"ACGTACGTACGT"; // 12 bp
        let mut w = b"AAAAAAAAAAAA".to_vec();
        w.extend_from_slice(adapter);
        let c = cfg(vec![ad("a", adapter, End::Three)], false);
        assert_eq!(adapter_segments(&w, &c), vec![(0, 12)]);
    }

    #[test]
    fn overlapping_interior_cuts_merge() {
        // Two DIFFERENT interior adapters whose hits overlap by 6 bp: `a` matches
        // [24,40), `b` matches [34,50) — constructed so their shared 6 bp region
        // ("TGTGTG", the tail of `a` / head of `b`) is literally the same window
        // bytes, giving both an exact (cost 0) hit. The overlap must merge into one
        // excision, leaving exactly 2 segments (not 3).
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
            "overlapping interior cuts merge into one excision"
        );
        assert_eq!(segs[0], (0, 24));
        assert_eq!(segs[1], (50, w.len()));
    }

    #[test]
    fn straddling_cut_is_clipped_not_leaked() {
        // The terminal hit [0,16) overlaps the interior hit [10,26).
        // Clipping the interior interval to the keep window must still excise [16,26).
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
            "s is fully excised, no leaked bases before 26"
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

    #[test]
    fn interior_above_k_mid_does_not_split() {
        // The cost-4 hit is within k_end=6 but exceeds the interior k_mid=3 limit.
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
        // end_size=10 keeps the intended hit AND sassy's other fuzzy hits (found
        // under the wide k_end=6 budget) away from the near-end terminal checks.
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
            "cost 4 hit is above k_mid=3 and must not split the read"
        );
    }

    #[test]
    fn ends_only_equals_split_on_indel_terminal_adapter() {
        // A six-base insertion expands the terminal alignment from 20 to 26 bases.
        // The terminal search window includes `k_end` additional bases so ends-only
        // and split modes select the same [2,28) alignment.
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
            "split-mode finds the full 26bp indel-bearing hit and trims to 28"
        );
        assert_eq!(
            ends_only_segs, split_segs,
            "ends-only must match split-mode exactly: the end-zone must be wide \
             enough (end_size + len + k_end) to contain the full indel-lengthened hit"
        );
        // Adapter bases actually removed, not just nominally "equal but empty".
        assert_eq!(ends_only_segs[0].0, 28);
    }

    #[test]
    fn three_prime_both_adapter_on_short_read_trims_tail_not_whole_read() {
        // 40bp insert + 20bp adapter at the 3' end; End::Both; end_size >= n so both
        // zones overlap. Must keep the insert [0,40), NOT drop the read.
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
        assert_eq!(adapter_segments(&w, &split), vec![(0, 40)], "split mode");
        assert_eq!(adapter_segments(&w, &ends), vec![(0, 40)], "ends-only mode");
    }

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

    #[test]
    fn both_adapters_at_both_ends_keep_middle() {
        // The two terminal patterns and their reverse complements are distinct,
        // leaving the 40-base insert as the only retained segment.
        let a5 = b"GGGGTTTTGGGGTTTTGGGG";
        let a3 = b"AAAAGGGGAAAAGGGGAAAA"; // A/G only (purine): NOT self-complementary, NOT revcomp(a5)
        let mut w = a5.to_vec();
        w.extend_from_slice(&[b'T'; 40]); // insert bytes don't match either adapter's revcomp
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

    #[test]
    fn inferred_single_end_adapters_on_short_read_keep_insert() {
        // Distinct end-specific adapters must preserve the insert when the end
        // search regions overlap on a short read.
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
        assert_eq!(adapter_segments(&w, &c), vec![(20, 60)], "insert survives");
    }
}
