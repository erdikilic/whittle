pub mod preset;
pub mod search;

use search::{hits, new_searcher};

/// Which read end a catalog sequence is expected at. `Both` is searched at both
/// ends and (when splitting) in the interior.
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
    if n == 0 || cfg.adapters.is_empty() {
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
        .filter(|&(s, e)| s >= lo && e <= hi && s < e)
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
}
