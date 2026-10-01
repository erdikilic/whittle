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
    /// The all-ACGT searcher over the forward strand only, for the search
    /// of one strand at a time in `search_pairs`.
    pub(super) plain_fwd: PlainSearcher,
    /// The ambiguity-tolerant searcher over the forward strand only.
    pub(super) ambiguous_fwd: AmbiguousSearcher,
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
    /// Hits of the tiled singleton batches; see `search_singleton_batches`.
    pub(super) pooled: Vec<Pooled>,
}

impl ThreadState {
    /// Creates the searchers with empty buffers.
    pub(super) fn new() -> Self {
        Self {
            plain: new_searcher(),
            ambiguous: new_ambiguous_searcher(),
            overhang: None,
            plain_fwd: new_plain_searcher_fwd(),
            ambiguous_fwd: new_searcher_fwd(),
            normalized: Vec::new(),
            reversed: Vec::new(),
            windows: Vec::new(),
            mask: MaskScratch::default(),
            head_flags: Vec::new(),
            tail_flags: Vec::new(),
            pooled: Vec::new(),
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
    /// The DNA-profile searcher over the forward strand only.
    pub(super) plain_fwd: &'a mut PlainSearcher,
    /// The IUPAC-profile searcher over the forward strand only.
    pub(super) ambiguous_fwd: &'a mut AmbiguousSearcher,
    /// Candidate windows of the interior search; see `candidate_windows`.
    pub(super) windows: &'a mut Vec<(usize, usize, usize)>,
    /// Masked end windows of the residue search.
    pub(super) mask: &'a mut MaskScratch,
    /// Per-adapter end-seed flags for the head window; see `end_candidates`.
    pub(super) head_flags: &'a mut Vec<bool>,
    /// Per-adapter end-seed flags for the tail window.
    pub(super) tail_flags: &'a mut Vec<bool>,
    /// Hits of the tiled singleton batches; see `search_singleton_batches`.
    pub(super) pooled: &'a mut Vec<Pooled>,
}

/// A span `[start, end)` of the read that is searched as a read of its own.
pub(super) type Span = (usize, usize);

/// Searches every equal-length barcode batch over the two end windows of the
/// span. All adapters in a batch share a length and budget, so the windows are
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

/// The texts one whole-pattern search covers at the two ends of a span, and
/// the hit ends each text owns.
pub(super) struct EndTexts<'a> {
    /// The texts; the first `count` are set.
    windows: [Strands<'a>; 4],
    /// Per text: its offset in the span, the site of its end window, and the
    /// range of hit ends, in span coordinates, that it owns.
    owned: [(usize, Site, usize, usize); 4],
    /// The number of texts.
    count: usize,
}

impl<'a> EndTexts<'a> {
    /// Returns the texts for a pattern of `len` bases at budget `k_end` over
    /// the end windows of the span, with zones of `end_size` bases.
    ///
    /// An end window of at least four alignment lengths is cut into two
    /// texts at a split point: the first owns the hits ending at or before
    /// it, the second those ending after it. An alignment spans at most
    /// `reach` bases, so each text extends `reach` bases past its owned
    /// ends, and every owned end sees the costs the whole window gives it.
    /// The four texts fill sassy's four lanes.
    fn new(ctx: Context<'a>, span: Span, end_size: usize, len: usize, k_end: usize) -> Self {
        let (ws, we) = span;
        let n = we - ws;
        let (head_end, tail_start) = terminal_windows(n, end_size, len, k_end);
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
        Self {
            windows,
            owned,
            count,
        }
    }

    /// Returns the site of text `text_idx` and `hit`, found in that text, in
    /// span coordinates, or `None` when the text does not own the hit.
    fn owned_hit(&self, text_idx: usize, hit: Hit) -> Option<(Site, Hit)> {
        let (offset, site, first_end, last_end) = self.owned[text_idx];
        let hit = shifted(hit, offset);
        (first_end..last_end)
            .contains(&hit.end)
            .then_some((site, hit))
    }
}

/// A hit of the tiled search of a singleton batch, held until
/// `search_singletons` reaches its adapter.
#[derive(Debug, Clone, Copy)]
pub(super) struct Pooled {
    /// Index into the configured adapters.
    adapter_idx: usize,
    /// The index of the text the hit was found in; see `EndTexts`.
    text_idx: usize,
    /// The site of that text.
    site: Site,
    /// The hit in span coordinates.
    hit: Hit,
}

impl Pooled {
    /// The position of the hit among the hits of its adapter as the
    /// single-pattern search reports them: the forward strand before the
    /// reverse, each strand text by text, and each text in scan order, by
    /// end on the forward strand and by descending start on the reverse
    /// strand, which is scanned over the reversed text.
    fn order(&self) -> (usize, bool, usize, usize) {
        let scan = if self.hit.rc {
            usize::MAX - self.hit.start
        } else {
            self.hit.end
        };
        (self.adapter_idx, self.hit.rc, self.text_idx, scan)
    }
}

/// Searches the tiled singleton batches (`TerminalBatch::tiled`) over the
/// texts of their members and fills `engine.pooled` with the owned hits,
/// sorted by adapter in the order of `Pooled::order`. The members of a batch
/// share their length and budget, so the texts and the budget are those of
/// each member's own search.
fn search_singleton_batches(
    ctx: Context<'_>,
    span: Span,
    end_size: usize,
    engine: &mut Engine<'_>,
) {
    engine.pooled.clear();
    for batch in &ctx.index.singleton_batches {
        if !batch.tiled(engine.plain_read) {
            continue;
        }
        let texts = EndTexts::new(ctx, span, end_size, batch.len, batch.k_end);
        for (text_idx, text) in texts.windows[..texts.count].iter().enumerate() {
            let pooled = &mut *engine.pooled;
            encoded_pattern_hits(
                engine.ambiguous,
                &batch.encoded,
                text.forward,
                text.reversed,
                batch.k_end,
                |pattern_idx, hit| {
                    if let Some((site, hit)) = texts.owned_hit(text_idx, hit) {
                        pooled.push(Pooled {
                            adapter_idx: batch.adapter_indices[pattern_idx],
                            text_idx,
                            site,
                            hit,
                        });
                    }
                },
            );
        }
    }
    engine.pooled.sort_unstable_by_key(Pooled::order);
}

/// Searches every adapter outside the barcode batches over the two end
/// windows of the span for whole-pattern hits, in adapter order, and passes
/// `keep` the hits of each (`singleton_hits`) that no cheaper overlapping
/// hit of the same adapter dominates.
pub(super) fn search_singletons(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
) {
    singleton_hits(ctx, span, keep.end_size, engine, |adapter_idx, found| {
        // Overlapping hits of one pattern are placements of one occurrence,
        // as a pattern with a repeating unit aligns one unit apart; the
        // cheapest applies.
        for &(site, h) in found {
            let dominated = found
                .iter()
                .any(|(_, g)| g.cost < h.cost && g.start < h.end && h.start < g.end);
            if !dominated {
                keep.accept(site, adapter_idx, h);
            }
        }
    });
}

/// Calls `visit` with each adapter outside the barcode batches, in adapter
/// order, and its whole-pattern hits over the two end windows of the span
/// with zones of `end_size` bases. An adapter of a tiled singleton batch
/// takes its hits from the tiled search of the batch
/// (`search_singleton_batches`); every other adapter is searched one pattern
/// at a time.
///
/// The two searches are separate sassy code paths, each with its own
/// traceback, so an equal-cost choice of hit start, as an indel beside a
/// homopolymer run allows, is not the same by construction. Both report the
/// rightmost local minima within the budget over the same texts, and
/// `segment_tests::tiled_singleton_batches_report_the_hits_of_the_one_by_one_search`
/// checks that they give every adapter the same sites and hits in the same
/// order, on planted entries with edits and homopolymer indels, reads with
/// `N`, spans shorter than a pattern, chimeras and random reads, under
/// sassy 0.2.6, the exact version `Cargo.toml` pins.
pub(super) fn singleton_hits(
    ctx: Context<'_>,
    span: Span,
    end_size: usize,
    engine: &mut Engine<'_>,
    mut visit: impl FnMut(usize, &[(Site, Hit)]),
) {
    search_singleton_batches(ctx, span, end_size, engine);
    let mut pooled = 0;
    let mut found: Vec<(Site, Hit)> = Vec::new();
    for (adapter_idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
        if !ctx.index.singletons[adapter_idx] {
            continue;
        }
        found.clear();
        let tiled = ctx.index.singleton_batch_of[adapter_idx]
            .is_some_and(|batch| ctx.index.singleton_batches[batch].tiled(engine.plain_read));
        if tiled {
            while let Some(p) = engine.pooled.get(pooled)
                && p.adapter_idx == adapter_idx
            {
                found.push((p.site, p.hit));
                pooled += 1;
            }
        } else {
            let Budget { len, k_end, .. } = ctx.index.budgets[adapter_idx];
            let texts = EndTexts::new(ctx, span, end_size, len, k_end);
            let accept = |text_idx: usize, h: Hit| found.extend(texts.owned_hit(text_idx, h));
            let windows = &texts.windows[..texts.count];
            if engine.plain_read && ctx.index.plain[adapter_idx] {
                for_each_hit_in_texts(engine.plain, &adapter.seq, windows, k_end, accept);
            } else {
                for_each_hit_in_texts(engine.ambiguous, &adapter.seq, windows, k_end, accept);
            }
        }
        visit(adapter_idx, &found);
    }
}

/// Searches the partial-matching entries with overhang alignment over each
/// end window of the span that the whole-pattern pass left untrimmed, for
/// the entries whose end seeds occur in that window. The end-seed flags are
/// set only for the partial-matching entries (see `CandidateIndex::new`), so
/// they gate the role as well as the seed. A partial hit flush with the read
/// end trims it; see `Keep::partial_hit_is_valid`. An end is untrimmed as
/// `Keep::five_open` and `Keep::three_open` judge it.
pub(super) fn search_partial(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
) {
    let (ws, we) = span;
    let n = we - ws;
    let head_open = keep.five_open();
    let tail_open = keep.three_open();
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
/// flags gate the entries as in `search_partial`. A hit of a split sheet
/// primer is held for `Keep::place_split` (`Keep::accept_residue`).
pub(super) fn search_residue(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
) {
    let (ws, we) = span;
    let n = we - ws;
    let window = ctx.read.window;
    let retry_head = keep.five_open();
    let retry_tail = keep.three_open();
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
                        let unmasked = (masked.max(h.start), h.end);
                        keep.accept_residue(adapter_idx, h, unmasked, HitAction::TrimFivePrime);
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
                        let hit = Hit {
                            start: s,
                            end: e,
                            ..h
                        };
                        let aligned = (s, start + unmasked.min(h.end));
                        keep.accept_residue(adapter_idx, hit, aligned, HitAction::TrimThreePrime);
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
        let pattern = &ctx.cfg.adapters[adapter_idx].seq;
        let text = ctx.read.strands(start, end);
        let Some(whole) = &ctx.index.whole_primers[adapter_idx] else {
            search(engine, ctx.index, adapter_idx, pattern, text, k_mid, |h| {
                keep.accept(Site::Interior, adapter_idx, shifted(h, start))
            });
            continue;
        };
        let mut found: Vec<Hit> = Vec::new();
        search(engine, ctx.index, adapter_idx, pattern, text, k_mid, |h| {
            found.push(shifted(h, start))
        });
        for hit in found {
            match whole_primer_over(ctx, engine, whole, hit) {
                Some(spanned) => keep.accept_standalone(adapter_idx, spanned),
                None => keep.accept(Site::Interior, adapter_idx, hit),
            }
        }
    }
}

/// Returns `hit` widened to the whole marker primer `whole` where that primer
/// aligns over it on the same strand within its interior budget for the
/// read, leaving fewer than `MIN_OVERLAP` bases to clip, or `None` where it
/// does not. The cheapest such alignment applies. A primer cut short by the
/// read end or by the search window does not align whole.
fn whole_primer_over(
    ctx: Context<'_>,
    engine: &mut Engine<'_>,
    whole: &WholePrimer,
    hit: Hit,
) -> Option<Hit> {
    let n = ctx.read.window.len();
    let k = whole.k_mid[interior_class(n)];
    let reach = whole.seq.len() + k;
    let start = hit.start.saturating_sub(reach);
    let end = (hit.end + reach).min(n);
    let mut best: Option<Hit> = None;
    for_each_hit(
        engine.ambiguous,
        &whole.seq,
        &ctx.read.strands(start, end),
        k,
        |g| {
            let g = shifted(g, start);
            if g.rc == hit.rc
                && g.start < hit.end
                && hit.start < g.end
                && g.clip_start + g.clip_end < MIN_OVERLAP
                && best.is_none_or(|b| g.cost < b.cost)
            {
                best = Some(g);
            }
        },
    );
    let g = best?;
    let (lo, hi) = (hit.start.min(g.start), hit.end.max(g.end));
    let inner_start = (hit.start + hit.clip_start).min(g.start + g.clip_start);
    let inner_end = (hit.end - hit.clip_end).max(g.end - g.clip_end);
    Some(Hit {
        start: lo,
        end: hi,
        clip_start: inner_start - lo,
        clip_end: hi - inner_end,
        ..hit
    })
}

/// Searches the `CandidateIndex::paired` entries over the interior of the
/// read at their pair budget and excises every junction pair: a hit that
/// reads out of the insert before it, followed by a hit that reads into the
/// insert after it, starting within `FLANK_SLACK` bases of its end or
/// overlapping it by at most as many. A chimera junction joins the end of one
/// molecule to the start of the next, so it holds the primers of both beside
/// each other in that orientation, whichever primers they are; the sites of a
/// marker gene in a genome lie a gene apart and face each other. Of the
/// overlapping hits of one entry, the cheapest stands for the occurrence.
/// The hits are those of a two-strand search of each entry, in its order.
pub(super) fn search_pairs(ctx: Context<'_>, engine: &mut Engine<'_>, keep: &mut Keep<'_>) {
    let n = ctx.read.window.len();
    // A pair needs a hit that reads out of an insert, so the strand on which
    // each entry reads out is searched first, and the other strand only in
    // a read that holds such a hit.
    let mut found: Vec<(usize, Hit)> = Vec::new();
    for closing in [true, false] {
        for (adapter_idx, &paired) in ctx.index.paired.iter().enumerate() {
            if !paired {
                continue;
            }
            let Budget { len, k_end, .. } = ctx.index.budgets[adapter_idx];
            let reach = len + k_end;
            let start = (keep.end_size + 1).saturating_sub(reach);
            let end = (n.saturating_sub(keep.end_size + 1) + reach).min(n);
            if end <= start || end - start < len {
                continue;
            }
            let opens = ctx.index.opens[adapter_idx];
            let pattern = &ctx.cfg.adapters[adapter_idx].seq;
            let text = ctx.read.strands(start, end);
            let k = ctx.index.budgets[adapter_idx].pair(n);
            let accept = |h: Hit| found.push((adapter_idx, shifted(h, start)));
            if opens == Opens::Both {
                if closing {
                    search(engine, ctx.index, adapter_idx, pattern, text, k, accept);
                }
                continue;
            }
            let rc = opens.reads_out(true) == closing;
            if engine.plain_read && ctx.index.plain[adapter_idx] {
                for_each_hit_on_strand(engine.plain_fwd, pattern, &text, rc, k, accept);
            } else {
                for_each_hit_on_strand(engine.ambiguous_fwd, pattern, &text, rc, k, accept);
            }
        }
        if found.is_empty() {
            return;
        }
    }
    // The hits of each entry in the order of a two-strand search, the
    // forward strand first; the sort is stable, so each strand keeps its
    // scan order.
    found.sort_by_key(|&(adapter_idx, hit)| (adapter_idx, hit.rc));
    let found: Vec<(usize, Hit)> = found
        .chunk_by(|a, b| a.0 == b.0)
        .flat_map(|hits| {
            hits.iter().copied().filter(|&(_, h)| {
                !hits.iter().any(|&(_, g)| {
                    g.rc == h.rc && g.cost < h.cost && g.start < h.end && h.start < g.end
                })
            })
        })
        .collect();
    if found.len() < 2 {
        return;
    }
    let opens =
        |&&(adapter_idx, hit): &&(usize, Hit)| ctx.index.opens[adapter_idx].reads_in(hit.rc);
    let closes =
        |&&(adapter_idx, hit): &&(usize, Hit)| ctx.index.opens[adapter_idx].reads_out(hit.rc);
    for closing in found.iter().filter(closes) {
        for opening in found.iter().filter(opens) {
            let (c, o) = (closing.1, opening.1);
            if o.start > c.start
                && o.end > c.end
                && o.start <= c.end + FLANK_SLACK
                && o.start + FLANK_SLACK >= c.end
            {
                keep.accept_pair(*closing, *opening);
            }
        }
    }
}

/// Runs the terminal passes over the span: the singleton and the barcode
/// batch whole-pattern searches, then the partial and residue searches over the
/// ends they left untrimmed, and, with a split sheet attached, the search of
/// the ends that another layer trimmed (`search_trimmed_ends`).
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
    if keep.five_open() || keep.three_open() {
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
    if !ctx.cfg.split_of.is_empty() {
        search_trimmed_ends(ctx, span, engine, keep);
    }
}

/// Searches the split sheet primers at each end of the span whose keep
/// boundary another layer moved and where no split primer hit is held
/// (`Keep::holds`), with that boundary standing in for the read end. A
/// primer that lost bases at its junction with the layer hangs off the
/// boundary, as one cut short by the read end hangs off the end: each end is
/// searched with overhang alignment over the `end_reach` bases inward of its
/// boundary, for the entries with an end seed there, and a hit is held under
/// the rules of a partial hit at a read end (`Keep::hold_at_boundary`). An
/// end that still holds no split primer hit is then searched for a whole
/// primer at the boundary (`search_boundary_primers`).
pub(super) fn search_trimmed_ends(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
) {
    let Some(table) = &ctx.index.end_seeds else {
        return;
    };
    let (ws, we) = span;
    let n = we - ws;
    let reach = ctx.index.end_reach;
    let ends = [
        (HitAction::TrimFivePrime, keep.lo > 0 && keep.lo < keep.hi),
        (HitAction::TrimThreePrime, keep.hi < n && keep.lo < keep.hi),
    ];
    for (action, trimmed) in ends {
        if !trimmed || keep.holds(action) {
            continue;
        }
        let (start, end) = match action {
            HitAction::TrimFivePrime => (keep.lo, (keep.lo + reach).min(keep.hi)),
            _ => (keep.hi.saturating_sub(reach).max(keep.lo), keep.hi),
        };
        if end - start < MIN_OVERLAP {
            continue;
        }
        let flags = &mut *engine.head_flags;
        flags.clear();
        flags.resize(ctx.cfg.adapters.len(), false);
        table.scan(&ctx.read.window[ws + start..ws + end], |adapter_idx, _| {
            flags[adapter_idx] = true
        });
        for (adapter_idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
            if !engine.head_flags[adapter_idx] || !keep.is_split(adapter_idx) {
                continue;
            }
            let k_end = ctx.index.budgets[adapter_idx].k_end;
            let text = ctx.read.strands(ws + start, ws + end);
            for_each_hit(engine.overhang, &adapter.seq, &text, k_end, |h| {
                keep.hold_at_boundary(adapter_idx, shifted(h, start), action);
            });
        }
        if !keep.holds(action) {
            search_boundary_primers(ctx, span, engine, keep, action);
        }
    }
}

/// Searches every split sheet primer for whole hits at the keep boundary of
/// the end `action` trims (`Keep::fits_boundary`): the outer edge, the start
/// at the 5' end or the end at the 3' end, lies at most
/// `BOUNDARY_OUTER_SLACK` bases outboard of the boundary or at most
/// `FLANK_SLACK` bases inboard of it, within the primer's anchored budget
/// (`Budget::k_anchor`). Of the hits that fit, each that no cheaper
/// overlapping fitting hit of the same primer dominates is held
/// (`Keep::hold_whole_at_boundary`), so a hit that does not fit, such as
/// one on the strand that faces away from the insert, suppresses none. The
/// boundary fixes where such a hit lies, so its chance of a random match is
/// that of the few positions beside the boundaries on one strand, which the
/// anchored budget bounds, not that of a whole end zone.
fn search_boundary_primers(
    ctx: Context<'_>,
    span: Span,
    engine: &mut Engine<'_>,
    keep: &mut Keep<'_>,
    action: HitAction,
) {
    let (ws, we) = span;
    let n = we - ws;
    let mut found: Vec<Hit> = Vec::new();
    for (adapter_idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
        if !keep.is_split(adapter_idx) {
            continue;
        }
        let Budget { len, k_anchor, .. } = ctx.index.budgets[adapter_idx];
        let reach = FLANK_SLACK + len + k_anchor;
        let (start, end) = match action {
            HitAction::TrimFivePrime => (
                keep.lo.saturating_sub(BOUNDARY_OUTER_SLACK),
                (keep.lo + reach).min(keep.hi),
            ),
            _ => (
                keep.hi.saturating_sub(reach).max(keep.lo),
                (keep.hi + BOUNDARY_OUTER_SLACK).min(n),
            ),
        };
        if end < start + len {
            continue;
        }
        found.clear();
        let text = ctx.read.strands(ws + start, ws + end);
        search(
            engine,
            ctx.index,
            adapter_idx,
            &adapter.seq,
            text,
            k_anchor,
            |h| found.push(shifted(h, start)),
        );
        found.retain(|h| keep.fits_boundary(adapter_idx, h, action));
        for &h in &found {
            let dominated = found
                .iter()
                .any(|g| g.cost < h.cost && g.start < h.end && h.start < g.end);
            if !dominated {
                keep.hold_whole_at_boundary(adapter_idx, h, action);
            }
        }
    }
}

/// Runs the search passes over per-thread state; see
/// `adapter_segments_annotated`.
pub(super) fn segments_tallied(
    window: &[u8],
    cfg: &AdapterConfig,
    acted: Option<&mut [bool]>,
) -> Vec<Segment> {
    let n = window.len();
    if n == 0 {
        return vec![];
    }
    if cfg.adapters.is_empty() {
        return vec![Segment::located(0, n, &[])];
    }
    with_engine(window, cfg, |ctx, engine| segments_with(ctx, engine, acted))
}

/// Runs `f` over the context of `window` under `cfg` and this thread's
/// searchers and buffers.
pub(super) fn with_engine<T>(
    window: &[u8],
    cfg: &AdapterConfig,
    f: impl FnOnce(Context<'_>, &mut Engine<'_>) -> T,
) -> T {
    let index = cfg
        .candidate_index
        .get_or_init(|| CandidateIndex::for_config(cfg));
    // The overhang cost per base of the terminal search is the error rate, so
    // a partial adapter costs what its missing part would have been allowed
    // in edits.
    let alpha = cfg.error_rate as f32;
    STATE.with_borrow_mut(|state| {
        let ThreadState {
            plain,
            ambiguous,
            overhang,
            plain_fwd,
            ambiguous_fwd,
            normalized,
            reversed,
            windows,
            mask,
            head_flags,
            tail_flags,
            pooled,
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
            plain_fwd,
            ambiguous_fwd,
            windows,
            mask,
            head_flags,
            tail_flags,
            pooled,
        };
        f(ctx, &mut engine)
    })
}

/// Asserts in debug builds that every primer `place_split` located over
/// `window` keeps the convention of its site (`site_holds`).
fn debug_assert_sites(window: &[u8], cfg: &AdapterConfig, located: &[Located]) {
    debug_assert!(
        located
            .iter()
            .all(|h| site_holds(window, &cfg.adapters, &h.locus)),
        "a located primer's site does not align over its locus: {located:?}"
    );
}

/// The search passes behind `adapter_segments_annotated`, over per-thread
/// searchers. `acted` receives the adapters that trimmed or excised, when
/// given. The loci of each segment are taken from the split primer hits of
/// the whole-window pass and of every piece pass.
pub(super) fn segments_with(
    ctx: Context<'_>,
    engine: &mut Engine<'_>,
    mut acted: Option<&mut [bool]>,
) -> Vec<Segment> {
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
        search_pairs(ctx, engine, &mut keep);
    }
    // A trim near an end found by the interior search, or a placed split
    // primer, can anchor a deferred terminal hit, so deferred hits are
    // settled before the tally counts the adapters that acted.
    keep.place_split();
    debug_assert_sites(ctx.read.window, cfg, &keep.primer_hits);
    keep.settle();
    tally(&keep);
    let (lo, hi, cuts, mut hits) = keep.into_cuts(cfg.min_piece);
    if lo >= hi {
        return vec![];
    }
    if cuts.is_empty() {
        return vec![Segment::located(lo, hi, &hits)];
    }

    // Each piece between cuts gets the terminal search again over its own
    // span, in ends-only mode. The end a cut created is a read end for every
    // purpose, and the opposite end has already been trimmed, so the second
    // pass changes it only when a further hit lies within `end_size` of the
    // new boundary.
    let mut segs = Vec::with_capacity(cuts.len() + 1);
    let mut cursor = lo;
    let mut push_piece = |s: usize, e: usize, segs: &mut Vec<Segment>| {
        if s >= e {
            return;
        }
        let mut keep = Keep::new(cfg, ctx.index, e - s, false);
        keep.gate_panels = gate_panels;
        search_terminal(ctx, (s, e), engine, &mut keep);
        keep.place_split();
        debug_assert_sites(&ctx.read.window[s..e], cfg, &keep.primer_hits);
        keep.settle();
        keep.refine();
        tally(&keep);
        hits.extend(keep.primer_hits.iter().map(|h| h.shifted(s)));
        if keep.lo < keep.hi {
            segs.push(Segment::located(s + keep.lo, s + keep.hi, &[]));
        }
    };
    for (s, e) in cuts {
        push_piece(cursor, s, &mut segs);
        cursor = cursor.max(e);
    }
    push_piece(cursor, hi, &mut segs);
    if !hits.is_empty() {
        for seg in &mut segs {
            *seg = Segment::located(seg.start, seg.end, &hits);
        }
    }
    segs
}
