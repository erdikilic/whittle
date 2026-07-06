pub mod preset;
pub mod search;

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
