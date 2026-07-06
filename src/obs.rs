//! Observability: leveled logging (tracing) and progress reporting (indicatif).

use tracing::level_filters::LevelFilter;

/// Map the CLI verbosity/quiet flags to a tracing level. `WHITTLE_LOG`, when set, is
/// applied separately (in `init`) and overrides this — unless `quiet` is set, in which
/// case `quiet` always wins (WARN) regardless of `WHITTLE_LOG`.
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
use tracing_subscriber::fmt::format::Writer;
use tracing_subscriber::fmt::time::FormatTime;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, fmt};

use crate::pipeline::{Counters, Stats};

/// How often the ticker thread refreshes the bar's position/message in `Mode::Bar`.
const TICK_INTERVAL: Duration = Duration::from_millis(250);
/// Steady-tick interval for the indicatif spinner shown when `total` is
/// unknown (no byte count to drive a determinate bar).
const SPINNER_TICK: Duration = Duration::from_millis(120);
/// Cadence of the periodic progress log line in `Mode::Line`. The ticker thread
/// sleeps in much shorter `TICK_INTERVAL` steps so `stop_ticker`'s join stays
/// prompt; it only actually logs once this much time has elapsed.
const LOG_INTERVAL: Duration = Duration::from_secs(10);

/// Local wall-clock timestamp `[YYYY-MM-DD HH:MM:SS]`, via `jiff`.
struct LocalStamp;

impl FormatTime for LocalStamp {
    fn format_time(&self, w: &mut Writer<'_>) -> std::fmt::Result {
        write!(w, "[{}]", jiff::Zoned::now().strftime("%Y-%m-%d %H:%M:%S"))
    }
}

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

/// The output mode for a run, computed once in `init` from `quiet`/`tty`/`verbosity`.
/// Exactly one applies — bar and line-log output never coexist.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Mode {
    /// `--quiet`: warnings/errors only. No bar, no progress line, no summary.
    Off,
    /// Default level on a real terminal: an animated bar/spinner, warnings/errors
    /// (suspended above it), and the final summary. No periodic log lines, no debug.
    Bar,
    /// `-v`/`-vv` on a TTY, or any non-TTY run: the start line, a periodic progress
    /// line every `LOG_INTERVAL`, debug/trace output (per level), and the summary.
    /// No bar.
    Line,
}

/// Owns the progress `MultiProgress`, the live ticker thread, and (in `Mode::Bar`)
/// the bar/spinner it drives. Created in the binary.
pub struct ProgressHandle {
    pub(crate) multi: MultiProgress,
    pub(crate) mode: Mode,
    ticker: Option<(Arc<AtomicBool>, JoinHandle<()>)>,
    bar: Option<ProgressBar>,
    /// Wall-clock start, set by `start()`; consumed by `finish()` to compute the
    /// summary's "in <dur> (<rate>)" tail. `None` if `start()` was never called
    /// (or after `finish()` has already consumed it).
    start: Option<Instant>,
}

impl ProgressHandle {
    /// A no-op handle for tests / library callers that install nothing.
    pub fn disabled() -> Self {
        ProgressHandle {
            multi: MultiProgress::new(),
            mode: Mode::Off,
            ticker: None,
            bar: None,
            start: None,
        }
    }

    /// True iff this run is in line mode (the periodic-log, no-bar mode) — either
    /// `-v`/`-vv` on a TTY, or any non-TTY run. Used to gate output that must not
    /// appear over a live bar (e.g. the one-time start line in `lib.rs`): bar mode
    /// stays clean of everything but the bar itself, warnings/errors, and the final
    /// summary. False in both `Mode::Bar` and `Mode::Off`.
    pub fn shows_lines(&self) -> bool {
        matches!(self.mode, Mode::Line)
    }

    /// Begin live progress once the input is open. `Mode::Bar` → animated bar/spinner
    /// driven by a ticker thread that polls the shared counters every `TICK_INTERVAL`;
    /// `Mode::Line` → no bar, the ticker instead emits a periodic INFO line every
    /// `LOG_INTERVAL`. `total` (input byte count) is `None` until byte counting lands;
    /// a `None` total renders a spinner rather than a bar. No-op in `Mode::Off`.
    pub fn start(&mut self, total: Option<u64>, counters: Arc<Counters>) {
        if matches!(self.mode, Mode::Off) {
            return;
        }
        debug_assert!(
            self.ticker.is_none(),
            "start() called twice without finish()"
        );
        let bar = if matches!(self.mode, Mode::Bar) {
            let pb = match total {
                Some(t) => {
                    let pb = self.multi.add(ProgressBar::new(t));
                    pb.set_style(
                        ProgressStyle::with_template(
                            "{elapsed_precise} [{bar:20}] {percent}% {msg} ETA {eta_precise}",
                        )
                        .unwrap()
                        .progress_chars("=>-"),
                    );
                    pb
                },
                None => {
                    let pb = self.multi.add(ProgressBar::new_spinner());
                    pb.set_style(
                        ProgressStyle::with_template("{elapsed_precise} {spinner} {msg}").unwrap(),
                    );
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
        let mode = self.mode;
        let start = Instant::now();
        self.start = Some(start);
        let handle = std::thread::spawn(move || {
            let mut last_log = start;
            while !stop_t.load(Ordering::Relaxed) {
                std::thread::sleep(TICK_INTERVAL);
                let ir = counters.input_reads.load(Ordering::Relaxed);
                let by = counters.bytes_read.load(Ordering::Relaxed);
                match mode {
                    Mode::Bar => {
                        if let Some(pb) = &bar_t {
                            if total.is_some() {
                                pb.set_position(by);
                            }
                            pb.set_message(bar_message(ir, by, total, start.elapsed()));
                        }
                    },
                    Mode::Line => {
                        if last_log.elapsed() >= LOG_INTERVAL {
                            tracing::info!("{}", periodic_line(ir, by, total, start.elapsed()));
                            last_log = Instant::now();
                        }
                    },
                    Mode::Off => break,
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

        let pct = if stats.input_reads > 0 {
            100.0 * stats.output_reads as f64 / stats.input_reads as f64
        } else {
            0.0
        };
        let mut msg = format!(
            "Kept {} of {} reads ({pct:.1}%)",
            commas(stats.output_reads),
            commas(stats.input_reads),
        );
        if let Some(start) = self.start.take() {
            let elapsed = start.elapsed();
            let rate = stats.input_reads as f64 / elapsed.as_secs_f64().max(1e-9);
            msg.push_str(&format!(
                " in {} ({} reads/s)",
                human_dur(elapsed),
                human_count(rate.round() as u64)
            ));
        }
        tracing::info!("{msg}");

        if stats.malformed_tag_reads > 0 {
            tracing::warn!(
                "Note: {} read(s) carried a per-base kinetics tag (ip/pw/fi/fp/ri/rp) whose \
                 length did not match the sequence; left unchanged",
                stats.malformed_tag_reads
            );
        }
    }
}

/// RAII backstop for early `?`/`bail!` returns from `run`/`run_folder` after
/// `start()` but before `finish()`: without this, the ticker thread and (in
/// `Mode::Bar`) the steady-tick spinner keep running after an error propagates, and
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
///
/// Precedence: `--quiet` always wins (WARN, regardless of `WHITTLE_LOG`); otherwise a
/// non-empty `WHITTLE_LOG` overrides `-v`/`-vv`; otherwise the level is derived from
/// `verbosity`.
///
/// Mode selection (never both bar and line log in the same run):
/// `quiet -> Off`; `!quiet && tty && verbosity==0 -> Bar`; otherwise (`-v`/`-vv` on a
/// TTY, or any non-TTY run) `-> Line`.
pub fn init(verbosity: u8, quiet: bool) -> ProgressHandle {
    let filter = if quiet {
        EnvFilter::new(level_from(verbosity, true).to_string())
    } else {
        match std::env::var("WHITTLE_LOG") {
            Ok(s) if !s.is_empty() => EnvFilter::new(s),
            _ => EnvFilter::new(level_from(verbosity, false).to_string()),
        }
    };
    let multi = MultiProgress::new();
    let tty = io::stderr().is_terminal();
    let mode = if quiet {
        Mode::Off
    } else if tty && verbosity == 0 {
        Mode::Bar
    } else {
        Mode::Line
    };
    tracing_subscriber::registry()
        .with(filter)
        .with(
            fmt::layer()
                .with_timer(LocalStamp)
                .with_level(true)
                .with_target(false)
                .with_ansi(tty)
                .with_writer(MpWriter {
                    multi: multi.clone(),
                }),
        )
        .init();
    ProgressHandle {
        multi,
        mode,
        ticker: None,
        bar: None,
        start: None,
    }
}

/// Compact magnitude for live progress fields: `750`, `145k`, `1.2M`.
fn human_count(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.0}k", n as f64 / 1_000.0)
    } else {
        n.to_string()
    }
}

/// Full thousands-separated integer, for the summary: `3,050,000`.
fn commas(n: u64) -> String {
    let digits = n.to_string();
    let bytes = digits.as_bytes();
    let mut out = String::with_capacity(bytes.len() + bytes.len() / 3);
    for (i, &b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(b as char);
    }
    out
}

/// Human-readable duration for the summary/debug lines: `420ms`, `1.42s`, `1m08s`, `1h02m`.
pub(crate) fn human_dur(d: Duration) -> String {
    let secs = d.as_secs_f64();
    if secs < 1.0 {
        format!("{}ms", d.as_millis())
    } else if secs < 60.0 {
        format!("{secs:.2}s")
    } else if secs < 3600.0 {
        let total = d.as_secs();
        format!("{}m{:02}s", total / 60, total % 60)
    } else {
        let total = d.as_secs();
        format!("{}h{:02}m", total / 3600, (total % 3600) / 60)
    }
}

/// `HH:MM:SS`-style duration, for the periodic line's ETA field (indicatif draws its
/// own ETA off the bar template; this covers the line-mode text equivalent).
fn fmt_hms(d: Duration) -> String {
    let s = d.as_secs();
    format!("{}:{:02}:{:02}", s / 3600, (s % 3600) / 60, s % 60)
}

/// Bar-mode message (the bar/spinner itself already draws elapsed, `%`, and ETA off
/// its own template — see `ProgressStyle` in `start()` — so this covers only the data
/// fields): `145k of 154k reads, 53 MB/s`. The "of <estimated total>" clause is only
/// shown when `total` (input bytes) is known, by extrapolating from the read/byte
/// ratio seen so far; a `None`/zero total (e.g. stdin, or folder-merge mode where byte
/// counting isn't wired up) drops it. `bytes == 0` likewise drops the MB/s field
/// rather than render a misleading rate.
fn bar_message(input_reads: u64, bytes: u64, total: Option<u64>, elapsed: Duration) -> String {
    let mut s = match total.filter(|&t| t > 0 && bytes > 0) {
        Some(t) => {
            let est_total = input_reads as f64 * (t as f64 / bytes as f64);
            format!(
                "{} of {} reads",
                human_count(input_reads),
                human_count(est_total.round() as u64)
            )
        },
        None => format!("{} reads", human_count(input_reads)),
    };
    if bytes > 0 {
        let secs = elapsed.as_secs_f64().max(1e-3);
        let mbps = (bytes as f64 / 1_000_000.0) / secs;
        s.push_str(&format!(", {mbps:.0} MB/s"));
    }
    s
}

/// Line-mode periodic progress log, emitted at INFO every `LOG_INTERVAL`:
/// `1,200,000 reads, 42%, 45k reads/s, 380 MB/s, ETA 00:00:40`. Fields, in order:
/// full-precision read count, `%` complete (if `total` bytes known), reads/s, MB/s
/// (if any bytes have been read), ETA (if `total` known).
fn periodic_line(input_reads: u64, bytes: u64, total: Option<u64>, elapsed: Duration) -> String {
    let secs = elapsed.as_secs_f64().max(1e-3);
    let rps = input_reads as f64 / secs;

    let mut s = format!("{} reads", commas(input_reads));
    if let Some(t) = total.filter(|&t| t > 0) {
        let pct = (100.0 * bytes as f64 / t as f64).min(100.0);
        s.push_str(&format!(", {pct:.0}%"));
    }
    s.push_str(&format!(", {} reads/s", human_count(rps.round() as u64)));
    if bytes > 0 {
        let mbps = (bytes as f64 / 1_000_000.0) / secs;
        s.push_str(&format!(", {mbps:.0} MB/s"));
    }
    if let Some(t) = total.filter(|&t| t > 0 && bytes > 0) {
        let bps = bytes as f64 / secs;
        let eta = Duration::from_secs_f64(((t.saturating_sub(bytes)) as f64 / bps).max(0.0));
        s.push_str(&format!(", ETA {}", fmt_hms(eta)));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression test for the RAII `Drop` cleanup: dropping a handle whose ticker
    /// thread is still running (e.g. via an early `?`/`bail!` return in `run` that
    /// never calls `finish()`) must stop and join that thread rather than leaking
    /// it or hanging the process. Built non-TTY / `Mode::Line` so no real
    /// spinner/terminal is needed; the ticker still spawns and must join cleanly.
    #[test]
    fn dropping_started_handle_stops_ticker_without_hanging() {
        let mut h = ProgressHandle {
            multi: MultiProgress::new(),
            mode: Mode::Line,
            ticker: None,
            bar: None,
            start: None,
        };
        h.start(None, Arc::new(Counters::default()));
        assert!(h.ticker.is_some(), "start() should have spawned a ticker");
        drop(h); // must join the ticker thread and return, not hang
    }

    /// Same regression as above, exercised through the `Mode::Bar` ticker branch
    /// (which additionally owns a live indicatif bar) rather than `Mode::Line`'s —
    /// both branches must stop/join/clear without hanging on drop.
    #[test]
    fn dropping_started_bar_handle_stops_ticker_without_hanging() {
        let mut h = ProgressHandle {
            multi: MultiProgress::new(),
            mode: Mode::Bar,
            ticker: None,
            bar: None,
            start: None,
        };
        h.start(Some(1_000), Arc::new(Counters::default()));
        assert!(h.ticker.is_some(), "start() should have spawned a ticker");
        assert!(h.bar.is_some(), "Mode::Bar should create a live bar");
        drop(h); // must join the ticker thread and clear the bar, not hang
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
    fn human_count_formats_magnitudes() {
        assert_eq!(human_count(999), "999");
        assert_eq!(human_count(1_000), "1k");
        assert_eq!(human_count(145_000), "145k");
        assert_eq!(human_count(1_200_000), "1.2M");
    }

    #[test]
    fn commas_inserts_thousands_separators() {
        assert_eq!(commas(0), "0");
        assert_eq!(commas(999), "999");
        assert_eq!(commas(1_000), "1,000");
        assert_eq!(commas(3_050_000), "3,050,000");
    }

    #[test]
    fn human_dur_formats_tiers() {
        assert_eq!(human_dur(Duration::from_millis(420)), "420ms");
        assert_eq!(human_dur(Duration::from_millis(1_420)), "1.42s");
        assert_eq!(human_dur(Duration::from_secs(68)), "1m08s");
        assert_eq!(human_dur(Duration::from_secs(3_720)), "1h02m");
    }

    #[test]
    fn periodic_line_without_total_has_no_percent_or_eta() {
        let s = periodic_line(1_200_000, 0, None, Duration::from_secs(10));
        assert!(s.contains("1,200,000 reads"));
        assert!(s.contains("reads/s"));
        assert!(!s.contains('%'));
        assert!(!s.contains("ETA"));
        assert!(
            !s.contains("MB/s"),
            "bytes=0 (untracked, e.g. folder-merge mode) must not render a misleading MB/s field: {s}"
        );
        assert!(
            !s.contains("->") && !s.contains('\u{b7}'),
            "plain ASCII only: {s}"
        );
    }

    #[test]
    fn periodic_line_with_total_adds_percent_and_eta() {
        let s = periodic_line(500, 42_000_000, Some(100_000_000), Duration::from_secs(2));
        assert!(s.contains("42%")); // 42MB / 100MB
        assert!(s.contains("MB/s"));
        assert!(s.contains("ETA"));
    }

    #[test]
    fn bar_message_without_total_has_no_of_clause() {
        let s = bar_message(145_000, 0, None, Duration::from_secs(60));
        assert_eq!(s, "145k reads");
    }

    #[test]
    fn bar_message_with_total_shows_estimated_total_and_rate() {
        let s = bar_message(
            145_000,
            50_000_000,
            Some(100_000_000),
            Duration::from_secs(60),
        );
        assert!(s.contains("145k of"));
        assert!(s.contains("reads"));
        assert!(s.contains("MB/s"));
        assert!(!s.contains('%'), "bar draws % itself via the template: {s}");
    }
}
