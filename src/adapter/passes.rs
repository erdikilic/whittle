//! The search passes over one read: per-thread searchers, read normalization,
//! the batched, singleton, partial, residue and interior searches, and the
//! re-search of the pieces an excision creates.

use super::*;

/// One thread's searchers and per-read buffers. Each buffer keeps its
/// capacity across reads, so the adapter stage itself allocates only the
/// segments it returns; the remaining per-read allocations are sassy's own,
/// inside each search call.
pub(super) struct ThreadState {
    /// The fast all-ACGT searcher. Used only when both the pattern and the
    /// searched text are plain ACGT; see `is_plain_acgt`.
    pub(super) plain: PlainSearcher,
    /// The ambiguity-tolerant searcher, for a degenerate primer and for the
    /// tiled terminal batches, which sassy implements for the IUPAC profile
    /// only.
    pub(super) ambiguous: AmbiguousSearcher,
    /// The overhang-aware searcher for the read ends, keyed by the overhang
    /// cost it was built with so a run at another error rate rebuilds it.
    pub(super) overhang: Option<(f32, AmbiguousSearcher)>,
    /// The normalized read, when the input is not its own normalization.
    pub(super) normalized: Vec<u8>,
    /// The normalized read reversed, for the reverse strand of every search.
    pub(super) reversed: Vec<u8>,
    /// Candidate windows as `(adapter, start, end)`; see `candidate_windows`.
    pub(super) windows: Vec<(usize, usize, usize)>,
    /// Masked end windows; see `search_residue`.
    pub(super) mask: MaskScratch,
    /// Per-adapter end-seed flags for the head window; see `end_candidates`.
    pub(super) head_flags: Vec<bool>,
    /// Per-adapter end-seed flags for the tail window.
    pub(super) tail_flags: Vec<bool>,
}

impl ThreadState {
    /// Creates the searchers with empty buffers.
    pub(super) fn new() -> Self {
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
pub(super) fn search(
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
pub(super) const AMBIGUOUS_READ_BASE: u8 = b'X';

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
pub(super) fn normalize_base(b: u8) -> u8 {
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
pub(super) fn is_upper_acgt(seq: &[u8]) -> bool {
    let upper_acgt = |b: u8| (b == b'A') | (b == b'C') | (b == b'G') | (b == b'T');
    let (chunks, remainder) = seq.as_chunks::<32>();
    let body = chunks
        .iter()
        .all(|chunk| chunk.iter().fold(true, |ok, &b| ok & upper_acgt(b)));
    body && remainder.iter().all(|&b| upper_acgt(b))
}

/// A normalized read and its reversal, borrowed from per-thread buffers so
/// every two-strand search copies nothing.
#[derive(Debug, Clone, Copy)]
pub(super) struct Read<'a> {
    /// The read as every searcher sees it; see `normalize_into`.
    pub(super) window: &'a [u8],
    /// `window` reversed.
    pub(super) reversed: &'a [u8],
}

impl<'a> Read<'a> {
    /// Returns the strands of `window[start..end]`.
    pub(super) fn strands(&self, start: usize, end: usize) -> Strands<'a> {
        let n = self.window.len();
        Strands {
            forward: &self.window[start..end],
            reversed: &self.reversed[n - end..n - start],
        }
    }
}

/// One read's search context: the configuration, its index, and the read.
#[derive(Clone, Copy)]
pub(super) struct Context<'a> {
    /// The run's adapter settings.
    pub(super) cfg: &'a AdapterConfig,
    /// The index built for `cfg.adapters`.
    pub(super) index: &'a CandidateIndex,
    /// The read under search.
    pub(super) read: Read<'a>,
}

/// The per-thread searchers and buffers one read is processed with, borrowed
/// from the thread's `ThreadState` for the duration of `adapter_segments`.
pub(super) struct Engine<'a> {
    /// Whether the complete read contains only ACGT bases.
    pub(super) plain_read: bool,
    /// The DNA-profile searcher, for a plain pattern.
    pub(super) plain: &'a mut PlainSearcher,
    /// The IUPAC-profile searcher, for a degenerate pattern and for the
    /// terminal batches.
    pub(super) ambiguous: &'a mut AmbiguousSearcher,
    /// The IUPAC-profile searcher with overhang alignment, for the read ends.
    pub(super) overhang: &'a mut AmbiguousSearcher,
    /// Candidate windows of the interior search; see `candidate_windows`.
    pub(super) windows: &'a mut Vec<(usize, usize, usize)>,
    /// Masked end windows of the residue search.
    pub(super) mask: &'a mut MaskScratch,
    /// Per-adapter end-seed flags for the head window; see `end_candidates`.
    pub(super) head_flags: &'a mut Vec<bool>,
    /// Per-adapter end-seed flags for the tail window.
    pub(super) tail_flags: &'a mut Vec<bool>,
}

/// A span `[start, end)` of the read that is searched as a read of its own.
pub(super) type Span = (usize, usize);

/// Searches every equal-length batch over the two end windows of the span.
/// All adapters in a batch share a length and budget, so the windows are
/// shared too; this collapses a kit's equal-length barcode searches into one
/// SIMD pattern search per end. An end trimmed through an entry that bounds
/// its barcode (see `bounds_barcode`) is not searched. Hits are passed to
/// `keep` in span coordinates.
pub(super) fn search_batched(
    ctx: Context<'_>,
    span: Span,
    searcher: &mut AmbiguousSearcher,
    keep: &mut Keep<'_>,
) {
    let (ws, we) = span;
    let n = we - ws;
    for batch in &ctx.index.terminal_batches {
        let (head_end, tail_start) = terminal_windows(n, keep.end_size, batch.len, batch.k_end);
        if !(keep.gate_panels && keep.five_bounded) {
            let head = ctx.read.strands(ws, ws + head_end);
            accept_batch_hits(batch, searcher, head, 0, Site::Head, keep);
        }
        if !(keep.gate_panels && keep.three_bounded) {
            let tail = ctx.read.strands(ws + tail_start, we);
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
}

/// Searches one batch over `text`, a window starting at `offset` in the span,
/// and passes every hit to `keep` at `site` in span coordinates.
pub(super) fn accept_batch_hits(
    batch: &TerminalBatch,
    searcher: &mut AmbiguousSearcher,
    text: Strands<'_>,
    offset: usize,
    site: Site,
    keep: &mut Keep<'_>,
) {
    let accept = |pattern_idx: usize, hit: Hit| {
        let adapter_idx = batch.adapter_indices[pattern_idx];
        if hit.cost > keep.budgets[adapter_idx].k_end {
            return;
        }
        keep.accept(site, adapter_idx, shifted(hit, offset));
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
pub(super) fn search_singletons(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
) {
    let (ws, we) = span;
    let n = we - ws;
    for (adapter_idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
        if !ctx.index.singletons[adapter_idx] {
            continue;
        }
        let Budget { len, k_end, .. } = ctx.index.budgets[adapter_idx];
        let (head_end, tail_start) = terminal_windows(n, keep.end_size, len, k_end);
        // An end window of at least four alignment lengths is cut into two
        // texts at a split point: the first owns the hits ending at or before
        // it, the second those ending after it. An alignment spans at most
        // `reach` bases, so each text extends `reach` bases past its owned
        // ends, and every owned end sees the costs the whole window gives it.
        // The four texts fill sassy's four lanes.
        let reach = len + k_end;
        let empty = ctx.read.strands(ws, ws);
        let mut windows = [empty; 4];
        let mut owned = [(0, Site::Head, 0, 0); 4];
        let mut count = 0;
        for (start, end, site) in [
            (0, head_end, Site::Head),
            (tail_start, n, Site::Tail { head_end }),
        ] {
            let parts = if end - start >= 4 * reach {
                let split = start + (end - start) / 2;
                [
                    (start, split + reach + 1, start, split + 1),
                    (split - reach, end, split + 1, end + 1),
                ]
            } else {
                [(start, end, start, end + 1), (0, 0, 0, 0)]
            };
            for (text_start, text_end, first_end, last_end) in parts {
                if text_end > text_start {
                    windows[count] = ctx.read.strands(ws + text_start, ws + text_end);
                    owned[count] = (text_start, site, first_end, last_end);
                    count += 1;
                }
            }
        }
        let accept = |text_idx: usize, h: Hit| {
            let (offset, site, first_end, last_end) = owned[text_idx];
            let h = shifted(h, offset);
            if (first_end..last_end).contains(&h.end) {
                keep.accept(site, adapter_idx, h);
            }
        };
        let windows = &windows[..count];
        if engine.plain_read && ctx.index.plain[adapter_idx] {
            for_each_hit_in_texts(engine.plain, &adapter.seq, windows, k_end, accept);
        } else {
            for_each_hit_in_texts(engine.ambiguous, &adapter.seq, windows, k_end, accept);
        }
    }
}

/// Searches the partial-matching entries with overhang alignment over each
/// end window of the span that the whole-pattern pass left untrimmed, for
/// the entries whose end seeds occur in that window. The end-seed flags are
/// set only for the partial-matching entries (see `CandidateIndex::new`), so
/// they gate the role as well as the seed. A partial hit flush with the read
/// end trims it; see `Keep::partial_hit_is_valid`.
pub(super) fn search_partial(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
) {
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
pub(super) const RESIDUE_MASKS: &[usize] = &[12, 24, 36];

/// Masked copies of one end window for the residue search, reused across reads.
#[derive(Debug, Default)]
pub(super) struct MaskScratch {
    /// The window with its outboard bases rewritten to `N`.
    pub(super) text: Vec<u8>,
    /// `text` reversed.
    pub(super) reversed: Vec<u8>,
}

impl MaskScratch {
    /// Fills the buffers with `window`, masking its first `mask` bases when
    /// `head` and its last `mask` bases otherwise, and returns the strands.
    pub(super) fn fill(&mut self, window: &[u8], mask: usize, head: bool) -> Strands<'_> {
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
pub(super) fn search_residue(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
) {
    let (ws, we) = span;
    let n = we - ws;
    let window = ctx.read.window;
    let retry_head = keep.lo == 0;
    let retry_tail = keep.hi == n;
    for (adapter_idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
        let retry_head = retry_head && engine.head_flags[adapter_idx];
        let retry_tail = retry_tail && engine.tail_flags[adapter_idx];
        if adapter.role != Role::Adapter || (!retry_head && !retry_tail) {
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
                        keep.trim_five_fixed(h.end);
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
                        keep.trim_three_fixed(s);
                    }
                });
            }
        }
    }
}

/// Returns `hit` moved `offset` bases to the right.
pub(super) fn shifted(hit: Hit, offset: usize) -> Hit {
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
pub(super) fn search_interior(ctx: Context<'_>, engine: &mut Engine<'_>, keep: &mut Keep<'_>) {
    let n = ctx.read.window.len();
    ctx.index.candidate_windows(ctx.read.window, engine.windows);
    for i in 0..engine.windows.len() {
        let (adapter_idx, start, end) = engine.windows[i];
        if !keep.splits(adapter_idx) {
            continue;
        }
        let Budget { len, k_end, .. } = ctx.index.budgets[adapter_idx];
        // An interior hit acts only outside both end zones, where the terminal
        // search does not reach: it starts after `end_size` and ends before
        // `n - end_size`. The window keeps `len + k_end` bases of context past
        // that region, as `candidate_windows` does around a seed.
        let reach = len + k_end;
        let start = start.max((keep.end_size + 1).saturating_sub(reach));
        let end = end.min(n.saturating_sub(keep.end_size + 1) + reach);
        if end <= start || end - start < len.saturating_sub(k_end) {
            continue;
        }
        let k_mid = ctx.index.budgets[adapter_idx].interior(n);
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
pub(super) fn search_terminal(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
) {
    search_singletons(ctx, span, engine, keep);
    keep.settle();
    if !ctx.index.terminal_batches.is_empty() {
        search_batched(ctx, span, engine.ambiguous, keep);
        keep.settle();
    }
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

/// Runs the search passes over per-thread state; see `adapter_segments`.
pub(super) fn segments_tallied(
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
pub(super) fn segments_with(
    ctx: Context<'_>,
    engine: &mut Engine<'_>,
    mut acted: Option<&mut [bool]>,
) -> Vec<(usize, usize)> {
    let cfg = ctx.cfg;
    let n = ctx.read.window.len();
    let gate_panels = acted.is_none();
    let mut tally = |keep: &Keep<'_>| {
        if let Some(acted) = acted.as_deref_mut() {
            for &adapter_idx in &keep.acted {
                acted[adapter_idx] = true;
            }
        }
    };
    let mut keep = Keep::new(cfg, ctx.index, n, cfg.split);
    keep.gate_panels = gate_panels;
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
        keep.gate_panels = gate_panels;
        search_terminal(ctx, (s, e), engine, &mut keep);
        keep.refine();
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
