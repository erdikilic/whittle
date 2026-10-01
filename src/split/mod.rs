//! Primer split: read demultiplexing by primer or primer pool.
//!
//! `sheet` holds the target model and its TSV, FASTA and preset parsers.
//! `score` rescores sheet primers at a located primer locus. `classify` is
//! the pure logic that turns rescored primer evidence at a segment's ends
//! into a `Call`. `splitter` builds a sheet's keys and rescoring state and
//! calls each trimmed piece against them. `route` expands the `-o` template
//! into the output path of each key. `counters` tallies calls per key for the
//! run summary and the end-of-run report.

pub mod classify;
pub mod counters;
pub mod route;
pub mod score;
pub mod sheet;
pub mod splitter;

pub use classify::{
    Bounds, Call, Ends, KeyLevel, Keys, Rules, Score, Strand, Unassigned, check_length, classify,
    classify_at,
};
pub use counters::{KeyCounters, KeyStats, SplitCounters, SplitStats};
pub use route::{KeyTable, Owner, Template};
pub use score::{End, Scorer};
pub use sheet::{Primer, RESERVED, Require, Sheet, Target};
pub use splitter::{Action, SplitOptions, Splitter};
