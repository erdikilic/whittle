//! Observability: leveled logging (tracing) and progress reporting (indicatif).

use tracing::level_filters::LevelFilter;

/// Map the CLI verbosity/quiet flags to a tracing level. `WHITTLE_LOG`, when set,
/// is applied separately (in `init`) and takes precedence over this.
pub fn level_from(verbosity: u8, quiet: bool) -> LevelFilter {
    if quiet {
        LevelFilter::WARN
    } else {
        match verbosity {
            0 => LevelFilter::INFO,
            1 => LevelFilter::DEBUG,
            _ => LevelFilter::TRACE,
        }
    }
}

use std::io::{self, IsTerminal, Write};

use indicatif::MultiProgress;
use tracing_subscriber::fmt::MakeWriter;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, fmt};

use crate::pipeline::Stats;

/// A `MakeWriter` that routes each fmt write through `MultiProgress::suspend`, so log
/// lines are printed cleanly above the live progress bar (and plainly when no bar exists).
#[derive(Clone)]
struct MpWriter {
    multi: MultiProgress,
}

impl Write for MpWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.multi.suspend(|| {
            let mut err = io::stderr().lock();
            err.write_all(buf)?;
            Ok(buf.len())
        })
    }
    fn flush(&mut self) -> io::Result<()> {
        io::stderr().flush()
    }
}

impl<'a> MakeWriter<'a> for MpWriter {
    type Writer = MpWriter;
    fn make_writer(&'a self) -> Self::Writer {
        self.clone()
    }
}

/// Owns the progress `MultiProgress` and (later) the ticker. Created in the binary.
///
/// `multi`/`enabled`/`tty` aren't read anywhere yet — they're the seam a later task
/// wires up (per-file progress bars gated on `enabled`/`tty`); `#[allow(dead_code)]`
/// keeps clippy clean in the meantime without dropping the fields early.
#[allow(dead_code)]
pub struct ProgressHandle {
    pub(crate) multi: MultiProgress,
    pub(crate) enabled: bool,
    pub(crate) tty: bool,
}

impl ProgressHandle {
    /// A no-op handle for tests / library callers that install nothing.
    pub fn disabled() -> Self {
        ProgressHandle {
            multi: MultiProgress::new(),
            enabled: false,
            tty: false,
        }
    }

    /// Print the end-of-run summary through tracing (subject to the level filter).
    pub fn finish(&self, stats: &Stats) {
        tracing::info!(
            "Kept {} reads out of {}",
            stats.output_reads,
            stats.input_reads
        );
        if stats.malformed_tag_reads > 0 {
            tracing::warn!(
                "note: {} read(s) carried a per-base kinetics tag (ip/pw/fi/fp/ri/rp) whose \
                 length did not match the sequence; left unchanged",
                stats.malformed_tag_reads
            );
        }
    }
}

/// Install the global subscriber and return the progress handle. Call once, in the binary.
pub fn init(verbosity: u8, quiet: bool) -> ProgressHandle {
    let filter = match std::env::var("WHITTLE_LOG") {
        Ok(s) if !s.is_empty() => EnvFilter::new(s),
        _ => EnvFilter::new(level_from(verbosity, quiet).to_string()),
    };
    let multi = MultiProgress::new();
    tracing_subscriber::registry()
        .with(filter)
        .with(
            fmt::layer()
                .without_time()
                .with_target(false)
                .with_level(false)
                .with_writer(MpWriter {
                    multi: multi.clone(),
                }),
        )
        .init();
    ProgressHandle {
        multi,
        enabled: !quiet,
        tty: io::stderr().is_terminal(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn level_mapping() {
        assert_eq!(level_from(0, false), LevelFilter::INFO);
        assert_eq!(level_from(1, false), LevelFilter::DEBUG);
        assert_eq!(level_from(2, false), LevelFilter::TRACE);
        assert_eq!(level_from(5, false), LevelFilter::TRACE);
        // quiet wins over verbosity
        assert_eq!(level_from(0, true), LevelFilter::WARN);
        assert_eq!(level_from(3, true), LevelFilter::WARN);
    }
}
