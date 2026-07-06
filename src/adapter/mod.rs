pub mod preset;
pub mod search;

use search::{hits, new_searcher};

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
}

/// Sequences shorter than this are never searched standalone — a <11 bp pattern
/// matches almost anywhere under any error budget. The 7 bp catalog flanks are
/// construction anchors, not standalone patterns.
pub const MIN_PATTERN_LEN: usize = 11;

/// Compute adapter keep-segments for `window`:
///   - terminal hits within `end_size` of an end trim that end inward;
///   - interior hits (stricter `k_mid`) excise and split.
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
    let mut searcher = new_searcher();

    let mut lo = 0usize; // 5' keep-boundary (advances inward on terminal hits)
    let mut hi = n; // 3' keep-boundary (retreats inward on terminal hits)
    let mut interior: Vec<(usize, usize)> = Vec::new();

    for ad in &cfg.adapters {
        let len = ad.seq.len();
        if len < MIN_PATTERN_LEN {
            continue;
        }
        let k_end = (cfg.error_rate * len as f64).floor() as usize;
        let k_mid = (0.5 * cfg.error_rate * len as f64).floor() as usize;
        for h in hits(&mut searcher, &ad.seq, window, k_end) {
            let near5 = h.start <= end_size && matches!(ad.end, End::Five | End::Both);
            let near3 =
                h.end >= n.saturating_sub(end_size) && matches!(ad.end, End::Three | End::Both);
            if near5 {
                lo = lo.max(h.end);
            } else if near3 {
                hi = hi.min(h.start);
            } else if cfg.split && h.cost <= k_mid {
                interior.push((h.start, h.end));
            }
        }
    }

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
        if let Some(last) = merged.last_mut() {
            if s <= last.1 {
                last.1 = last.1.max(e);
                continue;
            }
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
        }
    }
    fn ad(name: &str, seq: &[u8], end: End) -> Adapter {
        Adapter {
            name: name.into(),
            seq: seq.to_vec(),
            end,
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
        // `t` (End::Five, matches at [0,16)) is a terminal hit that advances `lo` to
        // 16. `s` (End::Both) is a DIFFERENT adapter whose only hit is [10,26) —
        // straddling `lo`: it starts inside t's terminal span (< lo) and ends well
        // past it. `t` and `s` share a 6 bp overlap ("GTTGGT", t's tail / s's head)
        // so both get an exact (cost 0) hit from the same physical bytes.
        // end_size=9 keeps s's hit (start=10) out of the near-5' terminal check
        // (10 > end_size, so `h.start <= end_size` is false) so it is classified
        // interior, not terminal.
        //
        // Pre-fix, the interior filter required the WHOLE cut inside [lo, hi) and
        // dropped [10,26) entirely (10 < lo=16), leaking bytes [16,26) — which
        // belong to `s` — into the kept segment as (16, 60). Post-fix, the cut is
        // clipped to (16, 26) and excised, so the kept segment starts at 26.
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
        };
        let segs = adapter_segments(&w, &c);
        assert_eq!(
            segs,
            vec![(26, 60)],
            "s is fully excised, no leaked bases before 26"
        );
        // Explicit invariant: no kept segment contains s's sequence.
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
        // len=12, error_rate=0.5 -> k_end = floor(0.5*12) = 6, k_mid = floor(0.25*12) = 3.
        // A copy with exactly 4 substitutions (positions 1,4,7,10) sits in the
        // interior. Verified empirically (see probe run in review) that sassy finds
        // it at cost 4 — strictly between k_mid(3) and k_end(6) — so the hit is
        // found but must NOT be actioned as an interior cut: the read stays whole.
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
        };
        assert_eq!(
            adapter_segments(&w, &c),
            vec![(0, w.len())],
            "cost 4 hit is above k_mid=3 and must not split the read"
        );
    }
}
