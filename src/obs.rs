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
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use tracing_subscriber::fmt::MakeWriter;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, fmt};

use crate::pipeline::{Counters, Stats};

/// How often the ticker thread polls the shared counters and refreshes the
/// bar's message/position (or, off-TTY, checks whether it's time to log).
const TICK_INTERVAL: Duration = Duration::from_millis(250);
/// Steady-tick interval for the indicatif spinner shown when `total` is
/// unknown (no byte count to drive a determinate bar).
const SPINNER_TICK: Duration = Duration::from_millis(120);
/// Minimum gap between off-TTY throttled progress log lines.
const LOG_THROTTLE: Duration = Duration::from_secs(30);

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

/// Owns the progress `MultiProgress`, the live ticker thread, and (on a TTY) the
/// bar/spinner it drives. Created in the binary.
pub struct ProgressHandle {
    pub(crate) multi: MultiProgress,
    pub(crate) enabled: bool,
    pub(crate) tty: bool,
    ticker: Option<(Arc<AtomicBool>, JoinHandle<()>)>,
    bar: Option<ProgressBar>,
}

impl ProgressHandle {
    /// A no-op handle for tests / library callers that install nothing.
    pub fn disabled() -> Self {
        ProgressHandle {
            multi: MultiProgress::new(),
            enabled: false,
            tty: false,
            ticker: None,
            bar: None,
        }
    }

    /// Begin live progress once the input is open. TTY → animated bar/spinner driven by
    /// a ticker thread that polls the shared counters; non-TTY → the ticker emits a
    /// throttled INFO line every ~30s instead. `total` (input byte count) is `None`
    /// until byte counting lands; a `None` total renders a spinner rather than a bar.
    /// No-op when `enabled` is false (quiet mode).
    pub fn start(&mut self, total: Option<u64>, counters: Arc<Counters>) {
        if !self.enabled {
            return;
        }
        debug_assert!(
            self.ticker.is_none(),
            "start() called twice without finish()"
        );
        let bar = if self.tty {
            let pb = match total {
                Some(t) => {
                    let pb = self.multi.add(ProgressBar::new(t));
                    pb.set_style(
                        ProgressStyle::with_template("{bar:30} {msg}")
                            .unwrap()
                            .progress_chars("=>-"),
                    );
                    pb
                },
                None => {
                    let pb = self.multi.add(ProgressBar::new_spinner());
                    pb.enable_steady_tick(SPINNER_TICK);
                    pb
                },
            };
            Some(pb)
        } else {
            None
        };

        let stop = Arc::new(AtomicBool::new(false));
        let stop_t = stop.clone();
        let bar_t = bar.clone();
        let tty = self.tty;
        let start = Instant::now();
        let handle = std::thread::spawn(move || {
            let mut last_log = Instant::now();
            while !stop_t.load(Ordering::Relaxed) {
                std::thread::sleep(TICK_INTERVAL);
                let ir = counters.input_reads.load(Ordering::Relaxed);
                let or = counters.output_reads.load(Ordering::Relaxed);
                let by = counters.bytes_read.load(Ordering::Relaxed);
                let msg = format_progress(ir, or, by, total, start.elapsed());
                if let Some(pb) = &bar_t {
                    if total.is_some() {
                        pb.set_position(by);
                    }
                    pb.set_message(msg);
                } else if !tty && last_log.elapsed() >= LOG_THROTTLE {
                    tracing::info!("{msg}");
                    last_log = Instant::now();
                }
            }
        });
        self.ticker = Some((stop, handle));
        self.bar = bar;
    }

    /// Stop the ticker (signal + join) and clear the bar, if either is live. Idempotent:
    /// both fields are `.take()`n, so a second call — including the one implicit in
    /// `Drop` after an explicit `finish()` — is a no-op. Shared by `finish()` (which
    /// follows it with the end-of-run summary) and `Drop` (which must clean up silently
    /// on an early `?`/`bail!` return, before the summary would otherwise ever print).
    fn stop_ticker(&mut self) {
        if let Some((stop, handle)) = self.ticker.take() {
            stop.store(true, Ordering::Relaxed);
            let _ = handle.join();
        }
        if let Some(pb) = self.bar.take() {
            pb.finish_and_clear();
        }
    }

    /// Stop the ticker and clear the bar, then print the end-of-run summary through
    /// tracing (subject to the level filter). The ticker is joined and the bar cleared
    /// *before* logging so no stale bar/spinner frame is left behind the summary line.
    pub fn finish(&mut self, stats: &Stats) {
        self.stop_ticker();
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

/// RAII backstop for early `?`/`bail!` returns from `run`/`run_folder` after
/// `start()` but before `finish()`: without this, the ticker thread and (on a
/// TTY) the steady-tick spinner keep running after an error propagates, and
/// the spinner can overwrite the fatal error message `main` prints. Stops and
/// joins the ticker and clears the bar — no summary logging, since an error
/// path has no `Stats` to summarize. A no-op after an explicit `finish()`
/// (both fields are already `None`).
impl Drop for ProgressHandle {
    fn drop(&mut self) {
        self.stop_ticker();
    }
}

/// Install the global subscriber and return the progress handle. Call once, in the binary.
pub fn init(verbosity: u8, quiet: bool) -> ProgressHandle {
    let filter = match std::env::var("WHITTLE_LOG") {
        Ok(s) if !s.is_empty() => EnvFilter::new(s),
        _ => EnvFilter::new(level_from(verbosity, quiet).to_string()),
    };
    let multi = MultiProgress::new();
    let tty = io::stderr().is_terminal();
    tracing_subscriber::registry()
        .with(filter)
        .with(
            fmt::layer()
                .without_time()
                .with_target(false)
                .with_level(false)
                .with_ansi(tty)
                .with_writer(MpWriter {
                    multi: multi.clone(),
                }),
        )
        .init();
    ProgressHandle {
        multi,
        enabled: !quiet,
        tty,
        ticker: None,
        bar: None,
    }
}

fn human_count(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.0}k", n as f64 / 1_000.0)
    } else {
        n.to_string()
    }
}

fn fmt_hms(d: Duration) -> String {
    let s = d.as_secs();
    format!("{}:{:02}:{:02}", s / 3600, (s % 3600) / 60, s % 60)
}

/// One-line progress summary shared by the bar message and the non-TTY log line.
/// `total` is the input byte count when known (adds a %-complete suffix).
pub fn format_progress(
    input_reads: u64,
    output_reads: u64,
    bytes: u64,
    total: Option<u64>,
    elapsed: Duration,
) -> String {
    let secs = elapsed.as_secs_f64().max(1e-3);
    let rps = input_reads as f64 / secs;
    let mbps = (bytes as f64 / 1_000_000.0) / secs;
    let kept_pct = if input_reads > 0 {
        100.0 * output_reads as f64 / input_reads as f64
    } else {
        0.0
    };
    let mut s = format!(
        "[{}] {} in → {} kept ({:.0}%) · {:.0}k reads/s",
        fmt_hms(elapsed),
        human_count(input_reads),
        human_count(output_reads),
        kept_pct,
        rps / 1000.0,
    );
    if bytes > 0 {
        s.push_str(&format!(" · {mbps:.0} MB/s"));
    }
    if let Some(t) = total.filter(|&t| t > 0) {
        let pct = (100.0 * bytes as f64 / t as f64).min(100.0);
        s.push_str(&format!(" · {pct:.0}%"));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression test for the RAII `Drop` cleanup: dropping a handle whose ticker
    /// thread is still running (e.g. via an early `?`/`bail!` return in `run` that
    /// never calls `finish()`) must stop and join that thread rather than leaking
    /// it or hanging the process. Built non-TTY so no real spinner/terminal is
    /// needed; `enabled: true` so `start()` actually spawns the ticker.
    #[test]
    fn dropping_started_handle_stops_ticker_without_hanging() {
        let mut h = ProgressHandle {
            multi: MultiProgress::new(),
            enabled: true,
            tty: false,
            ticker: None,
            bar: None,
        };
        h.start(None, Arc::new(Counters::default()));
        assert!(h.ticker.is_some(), "start() should have spawned a ticker");
        drop(h); // must join the ticker thread and return, not hang
    }

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

    #[test]
    fn progress_line_without_total_is_a_spinner_string() {
        let s = format_progress(
            1_200_000,
            1_050_000,
            0,
            None,
            std::time::Duration::from_secs(10),
        );
        assert!(s.contains("1.2M in"));
        assert!(s.contains("1.0M kept") || s.contains("1.1M kept"));
        assert!(s.contains("reads/s"));
        assert!(!s.contains('%') || s.contains("88%") || s.contains("87%")); // kept% present, no ETA%
        assert!(
            !s.contains("MB/s"),
            "bytes=0 (untracked, e.g. folder-merge mode) must not render a misleading MB/s field: {s}"
        );
    }

    #[test]
    fn progress_line_with_total_adds_percent() {
        let s = format_progress(
            500,
            400,
            42_000_000,
            Some(100_000_000),
            std::time::Duration::from_secs(2),
        );
        assert!(s.contains("42%")); // 42MB / 100MB
        assert!(s.contains("MB/s"));
    }
}
