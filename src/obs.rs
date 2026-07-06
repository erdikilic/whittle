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
use tracing::{Event, Level, Subscriber};
use tracing_subscriber::fmt::format::Writer;
use tracing_subscriber::fmt::{FmtContext, FormatEvent, FormatFields, MakeWriter};
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::registry::LookupSpan;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, fmt};

use crate::pipeline::{Counters, Stats};

/// How often the ticker thread refreshes the bar's position/message in `Mode::Bar`.
const TICK_INTERVAL: Duration = Duration::from_millis(250);
/// Steady-tick interval for the indicatif spinner shown when `total` is
/// unknown (no byte count to drive a determinate bar).
const SPINNER_TICK: Duration = Duration::from_millis(120);

/// Resolve the periodic-log cadence for `Mode::Line`: 30s by default, 10s at
/// `-v`/`-vv` (a verbose run wants more frequent feedback), overridable either
/// way via `WHITTLE_PROGRESS_INTERVAL` (integer seconds). Pure — takes the env
/// var's value as a parameter — so it's unit-testable without mutating real
/// process env (which would race across parallel test threads); the real
/// entry point, `resolve_log_interval`, is a thin wrapper reading the actual var.
fn log_interval_from(verbosity: u8, env_override: Option<&str>) -> Duration {
    if let Some(secs) = env_override.and_then(|s| s.parse::<u64>().ok()) {
        return Duration::from_secs(secs);
    }
    if verbosity >= 1 {
        Duration::from_secs(10)
    } else {
        Duration::from_secs(30)
    }
}

/// The ticker thread sleeps in much shorter `TICK_INTERVAL` steps so
/// `stop_ticker`'s join stays prompt; it only actually logs once
/// `log_interval_from`'s duration has elapsed. Malformed
/// `WHITTLE_PROGRESS_INTERVAL` values (unset, empty, non-numeric) are ignored.
fn resolve_log_interval(verbosity: u8) -> Duration {
    log_interval_from(
        verbosity,
        std::env::var("WHITTLE_PROGRESS_INTERVAL").ok().as_deref(),
    )
}

/// Custom event formatter: `[YYYY-MM-DD HH:MM:SS] [LEVEL] Message`, with a
/// bracketed local wall-clock timestamp (via `jiff`) AND a bracketed level.
/// The stock formatter's ` INFO`-padded, unbracketed level is what this
/// replaces. `Level`'s `Display` yields `INFO`/`WARN`/`DEBUG`/`TRACE`/`ERROR`,
/// so the level renders as `[INFO]` etc. `color` (set once in `init` from
/// whether stderr is a TTY) gates ANSI coloring of the `[LEVEL]` token only —
/// the timestamp and message always stay plain; when `false`, output carries
/// zero escape bytes, which matters for redirected/non-TTY runs.
struct WhittleFormat {
    color: bool,
}

/// ANSI color codes for the bracketed `[LEVEL]` token, raw (no new dependency):
/// ERROR bold red, WARN yellow, INFO green, DEBUG/TRACE dim. `color == false`
/// (non-TTY) yields the plain `[LEVEL]` token with no escape bytes at all.
fn format_level(level: &Level, color: bool) -> String {
    if !color {
        return format!("[{level}]");
    }
    let code = match *level {
        Level::ERROR => "\x1b[1;31m",
        Level::WARN => "\x1b[33m",
        Level::INFO => "\x1b[32m",
        Level::DEBUG | Level::TRACE => "\x1b[2m",
    };
    format!("{code}[{level}]\x1b[0m")
}

impl<S, N> FormatEvent<S, N> for WhittleFormat
where
    S: Subscriber + for<'a> LookupSpan<'a>,
    N: for<'a> FormatFields<'a> + 'static,
{
    fn format_event(
        &self,
        ctx: &FmtContext<'_, S, N>,
        mut w: Writer<'_>,
        event: &Event<'_>,
    ) -> std::fmt::Result {
        write!(
            w,
            "[{}] {} ",
            jiff::Zoned::now().strftime("%Y-%m-%d %H:%M:%S"),
            format_level(event.metadata().level(), self.color)
        )?;
        ctx.field_format().format_fields(w.by_ref(), event)?;
        writeln!(w)
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
    /// Default level on a real terminal: a one-line start banner, an animated
    /// bar/spinner, warnings/errors (suspended above it), and the final summary.
    /// No periodic log lines, no debug.
    Bar,
    /// `-v`/`-vv` on a TTY, or any non-TTY run: the full multi-line start banner, a
    /// periodic progress line every `log_interval` (see `resolve_log_interval`),
    /// debug/trace output (per level), and the summary. No bar.
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
    /// `Mode::Line` periodic-log cadence, resolved once in `init()` from
    /// verbosity/`WHITTLE_PROGRESS_INTERVAL` (see `resolve_log_interval`). Unused
    /// outside `Mode::Line`.
    log_interval: Duration,
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
            log_interval: Duration::from_secs(30),
        }
    }

    /// True iff this run is in line mode (the periodic-log, no-bar mode) — either
    /// `-v`/`-vv` on a TTY, or any non-TTY run. Used to gate output that must not
    /// appear over a live bar (e.g. the startup banner in `lib.rs`): bar mode
    /// stays clean of everything but its own one-line start, the bar itself,
    /// warnings/errors, and the final summary. False in both `Mode::Bar` and
    /// `Mode::Off`.
    pub fn shows_lines(&self) -> bool {
        matches!(self.mode, Mode::Line)
    }

    /// True iff this run is in bar mode (the animated bar/spinner, default level on
    /// a real terminal). Used to gate the one-line bar-mode start banner in
    /// `lib.rs` — unlike line mode's full multi-line banner, bar mode gets exactly
    /// one line so the bar stays clean. False in both `Mode::Line` and `Mode::Off`.
    pub fn is_bar(&self) -> bool {
        matches!(self.mode, Mode::Bar)
    }

    /// Begin live progress once the input is open. `Mode::Bar` → animated bar/spinner
    /// driven by a ticker thread that polls the shared counters every `TICK_INTERVAL`;
    /// `Mode::Line` → no bar, the ticker instead emits a periodic INFO line every
    /// `log_interval`. `total` (input byte count) is `None` until byte counting lands;
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
        let log_interval = self.log_interval;
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
                            pb.set_message(bar_message(ir, by, start.elapsed()));
                        }
                    },
                    Mode::Line => {
                        if last_log.elapsed() >= log_interval {
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
    /// `output` is the output path (or `<stdout>`) shown in the closing `Completed`
    /// line — the end-of-run counterpart to the startup banner's `Output:` line.
    /// `Completed` is always the last thing logged (after the malformed-tag note, if
    /// any) so it closes out the run symmetrically with the banner that opened it.
    /// Omitted when `elapsed` is unknown (a library caller using
    /// `ProgressHandle::disabled()`, which never calls `start()`).
    pub fn finish(&mut self, stats: &Stats, output: &str) {
        // Snapshot elapsed *before* `stop_ticker()`, not after: `stop_ticker` blocks
        // on `handle.join()`, and the ticker thread only wakes from its
        // `TICK_INTERVAL` (250ms) sleep to notice the stop flag — so measuring
        // afterward could add up to a full tick's worth of pure join-wait onto a
        // genuinely fast run, reporting e.g. "250ms" for a run that actually took
        // single-digit milliseconds. Capturing here first makes the reported
        // duration true wall-clock processing time, not processing time plus
        // ticker-shutdown latency.
        let elapsed = self.start.take().map(|start| start.elapsed());
        self.stop_ticker();

        tracing::info!("{}", summary_line(stats, elapsed));

        if let Some(line) = bases_line(stats) {
            tracing::info!("{}", line);
        }

        if let Some(line) = dropped_line(stats) {
            tracing::info!("{}", line);
        }

        if stats.malformed_tag_reads > 0 {
            tracing::warn!(
                "Note: {} read(s) carried a per-base kinetics tag (ip/pw/fi/fp/ri/rp) whose \
                 length did not match the sequence; left unchanged",
                stats.malformed_tag_reads
            );
        }

        if let Some(d) = elapsed {
            tracing::info!("{}", completed_line(d, output));
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
                .event_format(WhittleFormat { color: tty })
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
        log_interval: resolve_log_interval(verbosity),
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

/// Human-readable byte count for the startup banner's `Input:`/`Output:` fields:
/// `5.4 GB`, `183 MB`, `512 B`. Decimal (SI, 1000-based) units — consistent with
/// the MB/s figures already computed elsewhere in this module off `1_000_000.0`.
/// Bytes render as a bare integer (no fractional byte makes sense); above that,
/// values under 10 in their unit get one decimal place (`5.4 GB`), 10 and over
/// round to a whole number (`183 MB`).
pub(crate) fn human_bytes(n: u64) -> String {
    const UNITS: [&str; 6] = ["B", "KB", "MB", "GB", "TB", "PB"];
    let mut val = n as f64;
    let mut unit = 0usize;
    while val >= 1000.0 && unit + 1 < UNITS.len() {
        val /= 1000.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{n} {}", UNITS[0])
    } else if val < 10.0 {
        format!("{val:.1} {}", UNITS[unit])
    } else {
        format!("{val:.0} {}", UNITS[unit])
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

/// The end-of-run summary line: `Summary: 1 input reads, 3 output reads in 2.00s`.
/// Deliberately split-safe — no "kept X%" figure — because `--split-qual` can
/// turn one input read into several output segments, so `output_reads` can
/// legitimately exceed `input_reads` (a naive percentage would then read
/// "300%"). `elapsed` is `None` when the caller never started a timer (e.g. a
/// library caller using `ProgressHandle::disabled()`), in which case the
/// trailing "in <dur>" clause is omitted.
fn summary_line(stats: &Stats, elapsed: Option<Duration>) -> String {
    let mut msg = format!(
        "Summary: {} input reads, {} output reads",
        commas(stats.input_reads),
        commas(stats.output_reads),
    );
    if let Some(d) = elapsed {
        msg.push_str(&format!(" in {}", human_dur(d)));
    }
    msg
}

/// Human-readable base count for the yield summary's `Bases:` line: `12.4 Gbp`,
/// `460.0 Mbp`, `8.2 kbp`, `500 bp`. Decimal (1000-based) tiers, always one
/// decimal place above the `bp` tier — unlike `human_bytes`, there's no
/// "round to a whole number above 10" step, since a fixed one-decimal figure
/// reads more consistently across the Gbp-scale totals this line is built for.
fn human_bases(n: u64) -> String {
    if n >= 1_000_000_000 {
        format!("{:.1} Gbp", n as f64 / 1_000_000_000.0)
    } else if n >= 1_000_000 {
        format!("{:.1} Mbp", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1} kbp", n as f64 / 1_000.0)
    } else {
        format!("{n} bp")
    }
}

/// The end-of-run yield line: `Bases: 12.4 Gbp in, 11.9 Gbp out (95.8% kept)`.
/// Sits between `summary_line` and the malformed-tag/`Completed` lines (see
/// `finish`). `None` when `input_bases` is 0 (no bases were ever counted, e.g.
/// a library caller that never wired up the byte-level counters) — there is no
/// meaningful kept-percentage to report in that case.
fn bases_line(stats: &Stats) -> Option<String> {
    if stats.input_bases == 0 {
        return None;
    }
    let pct = 100.0 * stats.output_bases as f64 / stats.input_bases as f64;
    Some(format!(
        "Bases: {} in, {} out ({pct:.1}% kept)",
        human_bases(stats.input_bases),
        human_bases(stats.output_bases),
    ))
}

/// The end-of-run "why reads were dropped" line, shown right after `Bases:`:
/// `Dropped: 3,200 input reads (2,100 too short, 1,100 low quality)`. Only the
/// non-zero reasons appear, in this fixed order: too short, too long, low
/// quality, high quality, GC out of range, trimmed away (passed the filter but
/// every trimmed segment was empty/sub-`min_length`). `None` when nothing was
/// dropped, so a clean run gets no extra line.
fn dropped_line(stats: &Stats) -> Option<String> {
    let total = stats.dropped_short
        + stats.dropped_long
        + stats.dropped_low_qual
        + stats.dropped_high_qual
        + stats.dropped_gc
        + stats.dropped_trimmed;
    if total == 0 {
        return None;
    }
    let mut parts = Vec::new();
    if stats.dropped_short > 0 {
        parts.push(format!("{} too short", commas(stats.dropped_short)));
    }
    if stats.dropped_long > 0 {
        parts.push(format!("{} too long", commas(stats.dropped_long)));
    }
    if stats.dropped_low_qual > 0 {
        parts.push(format!("{} low quality", commas(stats.dropped_low_qual)));
    }
    if stats.dropped_high_qual > 0 {
        parts.push(format!("{} high quality", commas(stats.dropped_high_qual)));
    }
    if stats.dropped_gc > 0 {
        parts.push(format!("{} GC out of range", commas(stats.dropped_gc)));
    }
    if stats.dropped_trimmed > 0 {
        parts.push(format!("{} trimmed away", commas(stats.dropped_trimmed)));
    }
    Some(format!(
        "Dropped: {} input reads ({})",
        commas(total),
        parts.join(", ")
    ))
}

/// Human-readable duration for the summary/debug/closer lines: `420ms`, `1.42s`,
/// `1m08s`, `1h02m`. `pub`: `main.rs`'s failure path calls this directly to render
/// the "Failed after ..." elapsed time before any run-scoped state exists.
pub fn human_dur(d: Duration) -> String {
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

/// The end-of-run closer, emitted after the summary (and after the malformed-tag
/// note, if any) so it's always the true last line of a run — the end-of-run
/// counterpart to the startup banner's `Output:` line: `Completed in 2.00s,
/// output /path/to/out.fastq.gz`.
fn completed_line(elapsed: Duration, output: &str) -> String {
    format!("Completed in {}, output {output}", human_dur(elapsed))
}

/// `HH:MM:SS`-style duration, for the periodic line's ETA field (indicatif draws its
/// own ETA off the bar template; this covers the line-mode text equivalent).
fn fmt_hms(d: Duration) -> String {
    let s = d.as_secs();
    format!("{}:{:02}:{:02}", s / 3600, (s % 3600) / 60, s % 60)
}

/// Bar-mode message (the bar/spinner itself already draws elapsed, `%`, and ETA off
/// its own template — see `ProgressStyle` in `start()` — so this covers only the data
/// fields): `145k reads, 53 MB/s`. Shows just the processed read count — no invented
/// "of <total>" clause, since only total *bytes* are known up front, never total
/// *reads*; the byte-based `%`/ETA (drawn by the bar template) convey real progress.
/// `bytes == 0` (e.g. folder-merge mode where byte counting isn't wired up) drops the
/// MB/s field rather than render a misleading rate.
fn bar_message(input_reads: u64, bytes: u64, elapsed: Duration) -> String {
    let mut s = format!("{} reads", human_count(input_reads));
    if bytes > 0 {
        let secs = elapsed.as_secs_f64().max(1e-3);
        let mbps = (bytes as f64 / 1_000_000.0) / secs;
        s.push_str(&format!(", {mbps:.0} MB/s"));
    }
    s
}

/// Line-mode periodic progress log, emitted at INFO every `log_interval` (see
/// `resolve_log_interval`): `Processed 1,200,000 input reads, 42%, 45k reads/s,
/// 380 MB/s, ETA 00:00:40`. Fields, in order: full-precision *input* read count
/// (explicit — this is reads consumed, not reads emitted, which can legitimately
/// differ under `--split-qual`), `%` complete (if `total` bytes known), reads/s,
/// MB/s (if any bytes have been read), ETA (if `total` known).
fn periodic_line(input_reads: u64, bytes: u64, total: Option<u64>, elapsed: Duration) -> String {
    let secs = elapsed.as_secs_f64().max(1e-3);
    let rps = input_reads as f64 / secs;

    let mut s = format!("Processed {} input reads", commas(input_reads));
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
            log_interval: Duration::from_secs(30),
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
            log_interval: Duration::from_secs(30),
        };
        h.start(Some(1_000), Arc::new(Counters::default()));
        assert!(h.ticker.is_some(), "start() should have spawned a ticker");
        assert!(h.bar.is_some(), "Mode::Bar should create a live bar");
        drop(h); // must join the ticker thread and clear the bar, not hang
    }

    #[test]
    fn format_level_plain_has_no_escape_bytes() {
        for level in [
            Level::ERROR,
            Level::WARN,
            Level::INFO,
            Level::DEBUG,
            Level::TRACE,
        ] {
            let s = format_level(&level, false);
            assert!(
                !s.contains('\x1b'),
                "non-color output must carry zero escape bytes: {s:?}"
            );
        }
        assert_eq!(format_level(&Level::INFO, false), "[INFO]");
        assert_eq!(format_level(&Level::ERROR, false), "[ERROR]");
    }

    #[test]
    fn format_level_color_wraps_each_level_with_its_own_code_and_a_reset() {
        assert_eq!(
            format_level(&Level::ERROR, true),
            "\x1b[1;31m[ERROR]\x1b[0m"
        );
        assert_eq!(format_level(&Level::WARN, true), "\x1b[33m[WARN]\x1b[0m");
        assert_eq!(format_level(&Level::INFO, true), "\x1b[32m[INFO]\x1b[0m");
        assert_eq!(format_level(&Level::DEBUG, true), "\x1b[2m[DEBUG]\x1b[0m");
        assert_eq!(format_level(&Level::TRACE, true), "\x1b[2m[TRACE]\x1b[0m");
        for level in [
            Level::ERROR,
            Level::WARN,
            Level::INFO,
            Level::DEBUG,
            Level::TRACE,
        ] {
            let s = format_level(&level, true);
            assert!(s.contains('\x1b'), "color output must escape: {s:?}");
            assert!(s.ends_with("\x1b[0m"), "must always reset: {s:?}");
        }
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
    fn human_bytes_formats_magnitudes() {
        assert_eq!(human_bytes(0), "0 B");
        assert_eq!(human_bytes(512), "512 B");
        assert_eq!(human_bytes(183_000_000), "183 MB");
        assert_eq!(human_bytes(5_400_000_000), "5.4 GB");
        // 1000 bytes rolls over to the next unit rather than staying "1000 B".
        assert_eq!(human_bytes(1_000), "1.0 KB");
    }

    #[test]
    fn log_interval_defaults_to_30s_and_10s_when_verbose() {
        assert_eq!(log_interval_from(0, None), Duration::from_secs(30));
        assert_eq!(log_interval_from(1, None), Duration::from_secs(10));
        assert_eq!(log_interval_from(2, None), Duration::from_secs(10));
    }

    #[test]
    fn log_interval_env_override_wins_either_way() {
        assert_eq!(log_interval_from(0, Some("5")), Duration::from_secs(5));
        assert_eq!(log_interval_from(1, Some("60")), Duration::from_secs(60));
    }

    #[test]
    fn log_interval_ignores_unparseable_env_override() {
        assert_eq!(
            log_interval_from(0, Some("not-a-number")),
            Duration::from_secs(30)
        );
        assert_eq!(log_interval_from(0, Some("")), Duration::from_secs(30));
        assert_eq!(log_interval_from(1, Some("")), Duration::from_secs(10));
    }

    #[test]
    fn completed_line_formats_elapsed_and_output() {
        assert_eq!(
            completed_line(Duration::from_secs(2), "/tmp/out.fastq.gz"),
            "Completed in 2.00s, output /tmp/out.fastq.gz"
        );
        assert_eq!(
            completed_line(Duration::from_millis(420), "<stdout>"),
            "Completed in 420ms, output <stdout>"
        );
    }

    #[test]
    fn summary_line_is_split_safe_with_no_percentage() {
        // Regression: --split-qual can turn one input read into several output
        // segments, so output_reads > input_reads is legitimate. The summary
        // must not compute a "kept X%" figure off these counts.
        let stats = Stats {
            input_reads: 1,
            output_reads: 3,
            ..Default::default()
        };
        let s = summary_line(&stats, Some(Duration::from_secs(2)));
        assert_eq!(s, "Summary: 1 input reads, 3 output reads in 2.00s");
        assert!(!s.contains('%'));
        assert!(!s.contains("Kept"));
    }

    #[test]
    fn summary_line_omits_duration_when_elapsed_unknown() {
        let stats = Stats {
            input_reads: 5,
            output_reads: 5,
            ..Default::default()
        };
        assert_eq!(
            summary_line(&stats, None),
            "Summary: 5 input reads, 5 output reads"
        );
    }

    #[test]
    fn human_bases_formats_magnitudes() {
        assert_eq!(human_bases(12_400_000_000), "12.4 Gbp");
        assert_eq!(human_bases(460_000_000), "460.0 Mbp");
        assert_eq!(human_bases(8_240), "8.2 kbp");
        assert_eq!(human_bases(500), "500 bp");
    }

    #[test]
    fn bases_line_reports_kept_percentage() {
        let stats = Stats {
            input_reads: 1,
            output_reads: 1,
            input_bases: 12_400_000_000,
            output_bases: 11_900_000_000,
            ..Default::default()
        };
        assert_eq!(
            bases_line(&stats).unwrap(),
            "Bases: 12.4 Gbp in, 11.9 Gbp out (96.0% kept)"
        );
    }

    #[test]
    fn bases_line_omitted_when_input_bases_zero() {
        let stats = Stats {
            input_reads: 0,
            output_reads: 0,
            ..Default::default()
        };
        assert_eq!(bases_line(&stats), None);
    }

    #[test]
    fn dropped_line_lists_only_nonzero_reasons_in_fixed_order() {
        let stats = Stats {
            dropped_short: 2_100,
            dropped_low_qual: 1_100,
            ..Default::default()
        };
        assert_eq!(
            dropped_line(&stats).unwrap(),
            "Dropped: 3,200 input reads (2,100 too short, 1,100 low quality)"
        );
    }

    #[test]
    fn dropped_line_covers_every_reason_in_order() {
        let stats = Stats {
            dropped_short: 1,
            dropped_long: 2,
            dropped_low_qual: 3,
            dropped_high_qual: 4,
            dropped_gc: 5,
            dropped_trimmed: 6,
            ..Default::default()
        };
        assert_eq!(
            dropped_line(&stats).unwrap(),
            "Dropped: 21 input reads (1 too short, 2 too long, 3 low quality, \
             4 high quality, 5 GC out of range, 6 trimmed away)"
        );
    }

    #[test]
    fn dropped_line_omitted_when_total_zero() {
        assert_eq!(dropped_line(&Stats::default()), None);
    }

    #[test]
    fn periodic_line_without_total_has_no_percent_or_eta() {
        let s = periodic_line(1_200_000, 0, None, Duration::from_secs(10));
        assert!(s.contains("Processed 1,200,000 input reads"));
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
    fn bar_message_without_bytes_is_just_the_read_count() {
        let s = bar_message(145_000, 0, Duration::from_secs(60));
        assert_eq!(s, "145k reads");
    }

    #[test]
    fn bar_message_with_bytes_adds_rate_but_never_a_total() {
        let s = bar_message(145_000, 50_000_000, Duration::from_secs(60));
        assert!(s.starts_with("145k reads"));
        assert!(s.contains("MB/s"));
        assert!(
            !s.contains(" of "),
            "must not invent a total-reads figure: {s}"
        );
        assert!(!s.contains('%'), "bar draws % itself via the template: {s}");
    }
}
