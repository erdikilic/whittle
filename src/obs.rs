//! Observability: leveled logging (tracing) and progress reporting (indicatif).

use tracing::level_filters::LevelFilter;

use crate::config::{Advisory, Config, ProgressMode};
use tracing_subscriber::fmt::FormattedFields;

/// Maps the CLI verbosity and quiet flags to a tracing level. `WHITTLE_LOG`,
/// when set, is applied separately in `init` and overrides this unless `quiet`
/// is set; `quiet` yields WARN regardless of `WHITTLE_LOG`.
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

use crate::workflow::{Counters, Stats};

/// How often the ticker thread refreshes the bar's position and message in
/// `Mode::Bar`.
const TICK_INTERVAL: Duration = Duration::from_millis(250);
/// Steady-tick interval for the indicatif spinner shown when `total` is
/// unknown (no byte count to drive a determinate bar).
const SPINNER_TICK: Duration = Duration::from_millis(120);

/// Resolves the periodic-log cadence for `Mode::Line`: 30 s by default, 10 s at
/// `-v`/`-vv`; `WHITTLE_PROGRESS_INTERVAL` (integer seconds) overrides either.
/// Pure: the environment value is a parameter, so tests run without mutating
/// process environment, which races across parallel test threads.
/// `resolve_log_interval` reads the variable and delegates here.
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

/// Reads `WHITTLE_PROGRESS_INTERVAL` and resolves the periodic-log cadence
/// through `log_interval_from`. An unset, empty, or non-numeric value is
/// ignored. The ticker sleeps in `TICK_INTERVAL` steps so `stop_ticker` joins
/// promptly, and logs only once this cadence has elapsed.
fn resolve_log_interval(verbosity: u8) -> Duration {
    log_interval_from(
        verbosity,
        std::env::var("WHITTLE_PROGRESS_INTERVAL").ok().as_deref(),
    )
}

/// Custom event formatter: `[YYYY-MM-DD HH:MM:SS] [LEVEL] Message`, replacing
/// the stock formatter's ` INFO`-padded, unbracketed level. The timestamp is the
/// local wall clock via `jiff`. `color` (set once in `init` from whether stderr
/// is a TTY) gates ANSI coloring of the `[LEVEL]` token only, so a non-TTY run
/// carries no escape bytes.
struct WhittleFormat {
    /// Whether the `[LEVEL]` token is ANSI-colored.
    color: bool,
}

/// The bracketed `[LEVEL]` token, with ANSI color codes written inline: ERROR
/// bold red, WARN yellow, INFO green, DEBUG and TRACE dim. `color == false`
/// (non-TTY) yields the plain token with no escape bytes.
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

        // Enclosing span names, outermost first, so a line from deep in the
        // pipeline says which stage produced it. Nothing is printed for an event
        // outside every span, which is where the banner and summary lines sit.
        if let Some(scope) = ctx.event_scope() {
            let mut spans = scope.from_root().peekable();
            if spans.peek().is_some() {
                w.write_char('[')?;
                let mut first = true;
                for span in spans {
                    if !first {
                        w.write_char(':')?;
                    }
                    write!(w, "{}", span.name())?;
                    // A span's own fields identify the instance that produced
                    // the line: `read{name=...}` rather than a bare `read`.
                    let ext = span.extensions();
                    if let Some(fields) = ext.get::<FormattedFields<N>>()
                        && !fields.is_empty()
                    {
                        write!(w, "{{{fields}}}")?;
                    }
                    first = false;
                }
                w.write_str("] ")?;
            }
        }

        ctx.field_format().format_fields(w.by_ref(), event)?;
        writeln!(w)
    }
}

/// A `MakeWriter` that routes each fmt write through `MultiProgress::suspend`,
/// so log lines are printed cleanly above the live progress bar (and plainly
/// when no bar exists).
#[derive(Clone)]
struct MpWriter {
    /// The `MultiProgress` that writes are suspended around.
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

/// The output mode for a run, computed once in `init` from `quiet`, the TTY
/// state, and `verbosity`. Exactly one applies: bar and line-log output never
/// coexist.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Mode {
    /// `--quiet`: warnings and errors only. No bar, no progress line, no
    /// summary.
    Off,
    /// `--progress none`: the full banner, warnings and summary, but nothing
    /// reporting progress while the run is in flight.
    Silent,
    /// Default level on a terminal: a one-line start banner, an animated bar or
    /// spinner, warnings and errors (suspended above it), and the final summary.
    /// No periodic log lines, no debug.
    Bar,
    /// `-v`/`-vv` on a TTY, or any non-TTY run: the full multi-line start
    /// banner, a periodic progress line every `log_interval` (see
    /// `resolve_log_interval`), debug and trace output (per level), and the
    /// summary. No bar.
    Line,
}

/// Owns the progress `MultiProgress`, the live ticker thread, and (in
/// `Mode::Bar`) the bar or spinner it drives. Created in the binary.
pub struct ProgressHandle {
    /// The indicatif `MultiProgress` that log writes are suspended around.
    pub(crate) multi: MultiProgress,
    /// The output mode selected by `init`.
    pub(crate) mode: Mode,
    /// The ticker thread and its stop flag, while live.
    ticker: Option<(Arc<AtomicBool>, JoinHandle<()>)>,
    /// The live bar or spinner in `Mode::Bar`.
    bar: Option<ProgressBar>,
    /// Wall-clock start, set by `start`; consumed by `finish` to compute the
    /// summary's `in <dur>` tail. `None` if `start` was never called, or after
    /// `finish` has consumed it.
    start: Option<Instant>,
    /// `Mode::Line` periodic-log cadence, resolved once in `init` from the
    /// verbosity and `WHITTLE_PROGRESS_INTERVAL` (see `resolve_log_interval`).
    /// Unused outside `Mode::Line`.
    log_interval: Duration,
}

impl ProgressHandle {
    /// Returns a no-op handle for tests and library callers that install no
    /// subscriber.
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

    /// True in `Mode::Line` and `Mode::Silent`: the modes with a multi-line
    /// banner and no live bar. Gates output that must not appear over a bar,
    /// such as the startup banner in `lib.rs`. False in `Mode::Bar` and
    /// `Mode::Off`.
    pub fn shows_lines(&self) -> bool {
        matches!(self.mode, Mode::Line | Mode::Silent)
    }

    /// True in `Mode::Bar` only. Gates the one-line bar-mode start banner in
    /// `lib.rs`: bar mode prints one line so the bar stays clean, where line
    /// mode prints the full banner.
    pub fn is_bar(&self) -> bool {
        matches!(self.mode, Mode::Bar)
    }

    /// Begins live progress once the input is open. In `Mode::Bar` a ticker
    /// thread polls the shared counters every `TICK_INTERVAL` and drives a bar,
    /// or a spinner when `total` is `None` (stdin has no byte count). In
    /// `Mode::Line` the ticker emits a periodic INFO line every `log_interval`
    /// instead. In `Mode::Off` and `Mode::Silent` only the run clock starts.
    pub fn start(&mut self, total: Option<u64>, counters: Arc<Counters>) {
        // The timer runs in every mode, including `Off`: `--quiet` silences the
        // human-readable summary, but `elapsed_seconds` in `--summary-json` is
        // still reported.
        let start = Instant::now();
        self.start = Some(start);
        if matches!(self.mode, Mode::Off | Mode::Silent) {
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
                            "{elapsed_precise} [{bar:20}] {percent:>3}% {msg} ETA {eta_precise}",
                        )
                        .unwrap()
                        .progress_chars("=>-"),
                    );
                    // Seeded so the first frame carries the same fields as every
                    // later one; the ticker replaces it on its first pass.
                    pb.set_message(bar_message(0, 0, 0, Duration::ZERO));
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
        let handle = std::thread::spawn(move || {
            let mut last_log = start;
            while !stop_t.load(Ordering::Relaxed) {
                std::thread::sleep(TICK_INTERVAL);
                let ir = counters.input_reads.load(Ordering::Relaxed);
                let or = counters.output_reads.load(Ordering::Relaxed);
                let by = counters.bytes_read.load(Ordering::Relaxed);
                match mode {
                    Mode::Bar => {
                        if let Some(pb) = &bar_t {
                            if total.is_some() {
                                pb.set_position(by);
                            }
                            pb.set_message(bar_message(ir, or, by, start.elapsed()));
                        }
                    },
                    Mode::Line => {
                        if last_log.elapsed() >= log_interval {
                            tracing::info!("{}", periodic_line(ir, by, total, start.elapsed()));
                            last_log = Instant::now();
                        }
                    },
                    // Neither mode spawns this thread; stopping is the safe
                    // response to a mode changed underneath it.
                    Mode::Off | Mode::Silent => break,
                }
            }
        });
        self.ticker = Some((stop, handle));
        self.bar = bar;
    }

    /// Stops the ticker (signal and join) and clears the bar, if either is live.
    /// Idempotent: both fields are taken, so a second call (including the one in
    /// `Drop` after an explicit `finish`) is a no-op. Shared by `finish`, which
    /// follows it with the end-of-run summary, and `Drop`, which cleans up
    /// silently on an early error return.
    fn stop_ticker(&mut self) {
        if let Some((stop, handle)) = self.ticker.take() {
            stop.store(true, Ordering::Relaxed);
            // The default panic hook has already printed a ticker panic; the
            // assertion documents that a clean join is the only expected outcome.
            let joined = handle.join();
            debug_assert!(joined.is_ok(), "Progress ticker thread panicked");
        }
        if let Some(pb) = self.bar.take() {
            pb.finish_and_clear();
        }
    }

    /// Stops the ticker, clears the bar, then logs the end-of-run summary.
    /// Clearing happens first so no stale bar frame is left behind the summary.
    /// The closing `Completed` line is separate (`complete`), so the caller can
    /// write its artifacts between the two and a failed write is never reported
    /// after a success line.
    ///
    /// Returns the elapsed duration it reported, so `summary::Summary` quotes
    /// the same number.
    pub fn finish(&mut self, stats: &Stats) -> Option<Duration> {
        // Elapsed is taken before `stop_ticker`, which joins the ticker thread;
        // the thread notices the stop flag only when it wakes from its
        // `TICK_INTERVAL` sleep, so measuring afterward would charge up to a full
        // tick to a fast run.
        let elapsed = self.start.take().map(|start| start.elapsed());
        self.stop_ticker();

        tracing::info!("{}", summary_line(stats, elapsed));

        if let Some(line) = bases_line(stats) {
            tracing::info!("{}", line);
        }

        if let Some(line) = trimmed_to_nothing_line(stats) {
            tracing::info!("{}", line);
        }

        if let Some(line) = all_filtered_line(stats) {
            tracing::info!("{}", line);
        }

        if let Some(line) = segments_dropped_line(stats) {
            tracing::info!("{}", line);
        }

        if stats.malformed_tag_reads > 0 {
            tracing::warn!(
                reads = stats.malformed_tag_reads,
                "Per-base tags (ip/pw/fi/fp/ri/rp/sm/sx or a malformed sa) whose length did not match \
                 the sequence were left unchanged"
            );
        }
        if stats.malformed_mod_reads > 0 {
            tracing::warn!(
                reads = stats.malformed_mod_reads,
                "Malformed MM/ML/MN modification blocks were removed from the output"
            );
        }
        if stats.undo_tags_dropped_reads > 0 {
            tracing::warn!(
                reads = stats.undo_tags_dropped_reads,
                "PacBio undo blobs (ds/ls) were removed from trimmed reads"
            );
        }

        // Guardrail warnings: an empty input and an all-dropped run both exit
        // successfully, and both are reported at WARN so they are
        // distinguishable from a normal run in a log.
        if stats.input_reads == 0 {
            tracing::warn!("Input contained no reads");
        } else if stats.output_reads == 0 {
            tracing::warn!(
                input_reads = stats.input_reads,
                "No reads survived; every input read was dropped"
            );
        }

        elapsed
    }

    /// Logs the closing `Completed` line, the last line of a run. `output` is
    /// the path (or `<stdout>`) the reads went to; `elapsed` is the value
    /// `finish` returned, and the line is omitted when it is `None`.
    pub fn complete(&self, elapsed: Option<Duration>, output: &str) {
        if let Some(d) = elapsed {
            tracing::info!("{}", completed_line(d, output));
        }
    }
}

/// Backstop for an early error return between `start` and `finish`: without it
/// the ticker thread and the bar-mode spinner keep running past the error and
/// can overwrite the fatal message `main` prints. Stops the ticker and clears
/// the bar without a summary, since an error path has no `Stats`. A no-op after
/// an explicit `finish`.
impl Drop for ProgressHandle {
    fn drop(&mut self) {
        self.stop_ticker();
    }
}

/// Pure mode selection (see `init`): the TTY and `WHITTLE_LOG` state are
/// parameters, so tests run without mutating process environment, which races
/// across parallel test threads. `init` reads the real state and delegates here.
fn select_mode(
    quiet: bool,
    tty: bool,
    verbosity: u8,
    whittle_log_set: bool,
    progress: ProgressMode,
) -> Mode {
    // `--quiet` silences everything, including the summary, so it outranks any
    // progress preference.
    if quiet {
        return Mode::Off;
    }
    // Debug output and a live bar cannot share a terminal, so verbosity and a
    // `WHITTLE_LOG` filter both fall back to periodic lines, under an explicit
    // `--progress bar` as well as under `auto`. `auto` additionally needs a
    // terminal to redraw on.
    let bar_fits = verbosity == 0 && !whittle_log_set;
    match progress {
        ProgressMode::None => Mode::Silent,
        ProgressMode::Bar if bar_fits => Mode::Bar,
        ProgressMode::Auto if bar_fits && tty => Mode::Bar,
        ProgressMode::Bar | ProgressMode::Auto | ProgressMode::Plain => Mode::Line,
    }
}

/// Resolves the level filter from `--quiet`, `WHITTLE_LOG`, and the verbosity.
///
/// `--quiet` always wins (WARN); otherwise a non-empty `WHITTLE_LOG` overrides
/// `-v`/`-vv`; otherwise the level follows `verbosity`. A `WHITTLE_LOG` that does
/// not parse falls back to the verbosity level and returns an advisory naming it,
/// since a lossy parse would enable nothing and hide even the ERROR line of a
/// failing run. The returned flag says whether `WHITTLE_LOG` took effect, which
/// `select_mode` uses to keep its lines off a bar.
fn log_filter(
    whittle_log: Option<&str>,
    verbosity: u8,
    quiet: bool,
) -> (EnvFilter, bool, Option<Advisory>) {
    let fallback = || EnvFilter::new(level_from(verbosity, quiet).to_string());
    if quiet {
        return (fallback(), false, None);
    }
    let Some(spec) = whittle_log else {
        return (fallback(), false, None);
    };
    match EnvFilter::builder().parse(spec) {
        Ok(filter) => (filter, true, None),
        Err(e) => {
            let level = level_from(verbosity, false);
            let advisory = Advisory::warn(format!(
                "WHITTLE_LOG is not a valid log filter and is ignored: \
                 value={spec:?}, error={e}, level={level}"
            ));
            (fallback(), false, Some(advisory))
        },
    }
}

/// Installs the global subscriber and returns the progress handle. Called once,
/// in the binary.
///
/// The level follows `log_filter`; a rejected `WHITTLE_LOG` lands in
/// `cfg.advisories`, which `run` prints once the banner is up.
///
/// Mode is never both bar and line log: `quiet` gives `Off`, a default-verbosity
/// TTY with no `WHITTLE_LOG` gives `Bar`, everything else gives `Line`.
/// `WHITTLE_LOG` forces `Line` so its debug lines cannot interleave with a bar.
pub fn init(cfg: &mut Config) -> ProgressHandle {
    let (verbosity, quiet, progress) = (cfg.verbosity, cfg.quiet, cfg.progress);
    let whittle_log = std::env::var("WHITTLE_LOG").ok().filter(|s| !s.is_empty());
    let (filter, whittle_log_set, advisory) = log_filter(whittle_log.as_deref(), verbosity, quiet);
    cfg.advisories.extend(advisory);
    let multi = MultiProgress::new();
    let tty = io::stderr().is_terminal();
    let mode = select_mode(quiet, tty, verbosity, whittle_log_set, progress);
    tracing_subscriber::registry()
        .with(filter)
        .with(
            fmt::layer()
                .event_format(WhittleFormat { color: tty })
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
        log_interval: resolve_log_interval(verbosity),
    }
}

/// Compact magnitude for live progress fields: `750`, `145k`, `1.2M`.
fn human_count(n: u64) -> String {
    // Promote values that would round to 1000k so the result stays normalized.
    if n >= 999_500 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.0}k", n as f64 / 1_000.0)
    } else {
        n.to_string()
    }
}

/// Human-readable byte count for the startup banner's `Input:`/`Output:`
/// fields: `5.4 GB`, `183 MB`, `512 B`. Decimal (SI, 1000-based) units,
/// matching the MB/s figures in this module. Bytes render as a bare integer;
/// above that, values under 10 in their unit get one decimal place (`5.4 GB`)
/// and 10 and over round to a whole number (`183 MB`).
pub(crate) fn human_bytes(n: u64) -> String {
    const UNITS: [&str; 6] = ["B", "KB", "MB", "GB", "TB", "PB"];
    let mut val = n as f64;
    let mut unit = 0usize;
    // The threshold is `999.5`, not `1000.0`: a value that `{:.0}` would round
    // up to `1000` in its unit is promoted, so it reads `1.0 MB` rather than
    // `1000 KB`.
    while val >= 999.5 && unit + 1 < UNITS.len() {
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

/// The end-of-run summary line: `Summary: 1 input reads, 3 output reads in
/// 2.00s`. It carries no kept percentage, since `--qual-split` can turn one
/// input read into several segments and a read-count percentage would exceed
/// 100%. The trailing `in <dur>` clause is omitted when `elapsed` is `None`.
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
/// `460.0 Mbp`, `8.2 kbp`, `500 bp`. Uses decimal units and one decimal place
/// for kbp and larger values.
fn human_bases(n: u64) -> String {
    // Promote values that would round to 1000.0 in the current unit.
    if n >= 999_950_000 {
        format!("{:.1} Gbp", n as f64 / 1_000_000_000.0)
    } else if n >= 999_950 {
        format!("{:.1} Mbp", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1} kbp", n as f64 / 1_000.0)
    } else {
        format!("{n} bp")
    }
}

/// The end-of-run yield line: `Bases: 12.4 Gbp in, 11.9 Gbp out (95.8% kept)`.
/// Sits between `summary_line` and the malformed-tag and `Completed` lines (see
/// `finish`). `None` when `input_bases` is 0, where no kept percentage exists.
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

/// The end-of-run read-level line: `Trimmed to nothing: 1,234 input reads produced
/// no segments at all`. Covers reads for which `trim::apply` returned no
/// intervals: empty, fully consumed by adapter trimming, or over-cropped. Distinct
/// from `all_filtered_line`, which covers reads that did produce segments. `None`
/// when no read was trimmed to nothing.
fn trimmed_to_nothing_line(stats: &Stats) -> Option<String> {
    if stats.reads_trimmed_to_nothing == 0 {
        return None;
    }
    Some(format!(
        "Trimmed to nothing: {} input reads produced no segments at all",
        commas(stats.reads_trimmed_to_nothing)
    ))
}

/// The end-of-run read-level "every segment filtered" line, shown after
/// `trimmed_to_nothing_line`: `All segments filtered: 567 input reads had every
/// produced segment filtered`. Covers reads that produced at least one segment
/// (unlike `trimmed_to_nothing_line`) but had every one rejected by post-trim
/// `filter::check`. `None` when no input read had all of its segments filtered.
fn all_filtered_line(stats: &Stats) -> Option<String> {
    if stats.reads_all_filtered == 0 {
        return None;
    }
    Some(format!(
        "All segments filtered: {} input reads had every produced segment filtered",
        commas(stats.reads_all_filtered)
    ))
}

/// The end-of-run drop-reason line: `Segments dropped: 3,200 (2,100 too short,
/// 1,100 low quality)`. Only non-zero reasons appear, in a fixed order: too short,
/// too long, low quality, high quality, GC out of range. Counts segments, not
/// reads, since a split read can contribute several and still survive. `None` when
/// nothing was dropped.
fn segments_dropped_line(stats: &Stats) -> Option<String> {
    use crate::filter::DropReason;

    // Fixed order, and the wording comes from `DropReason::label` so this line
    // and the per-segment trace event cannot describe the same rejection
    // differently.
    let by_reason = [
        (DropReason::TooShort, stats.segments_dropped_short),
        (DropReason::TooLong, stats.segments_dropped_long),
        (DropReason::LowQuality, stats.segments_dropped_low_qual),
        (DropReason::HighQuality, stats.segments_dropped_high_qual),
        (DropReason::Gc, stats.segments_dropped_gc),
    ];
    let total: u64 = by_reason.iter().map(|(_, n)| n).sum();
    if total == 0 {
        return None;
    }
    let parts: Vec<String> = by_reason
        .iter()
        .filter(|(_, n)| *n > 0)
        .map(|(reason, n)| format!("{} {}", commas(*n), reason.label()))
        .collect();
    Some(format!(
        "Segments dropped: {} ({})",
        commas(total),
        parts.join(", ")
    ))
}

/// Human-readable duration for the summary, debug, and closer lines: `420ms`,
/// `1.42s`, `1m08s`, `1h02m`. `pub` because `main.rs`'s failure path renders
/// the `Failed after ...` elapsed time before any run-scoped state exists.
pub fn human_dur(d: Duration) -> String {
    let secs = d.as_secs_f64();
    if secs < 1.0 {
        format!("{}ms", d.as_millis())
    } else if secs < 59.995 {
        // Below the value `{:.2}s` would round up to "60.00s".
        format!("{secs:.2}s")
    } else if secs < 3599.5 {
        // Round to the nearest second so 59.996s reads "1m00s".
        let total = secs.round() as u64;
        format!("{}m{:02}s", total / 60, total % 60)
    } else {
        let total = secs.round() as u64;
        format!("{}h{:02}m", total / 3600, (total % 3600) / 60)
    }
}

/// The end-of-run closer, emitted after the summary (and after the
/// malformed-tag note, if any) so it is the last line of a run; the counterpart
/// to the startup banner's `Output:` line: `Completed in 2.00s, output
/// /path/to/out.fastq.gz`.
fn completed_line(elapsed: Duration, output: &str) -> String {
    format!("Completed in {}, output {output}", human_dur(elapsed))
}

/// `HH:MM:SS`-style duration for the periodic line's ETA field (indicatif draws
/// its own ETA from the bar template; this covers the line-mode equivalent).
fn fmt_hms(d: Duration) -> String {
    let s = d.as_secs();
    format!("{}:{:02}:{:02}", s / 3600, (s % 3600) / 60, s % 60)
}

/// Bar-mode message: `145k reads, 53 MB/s`. Covers only the data fields, since
/// the bar template draws elapsed, percent, and ETA. The processed read count
/// carries no `of <total>` because only total bytes are known up front.
/// `bytes == 0` omits the MB/s field.
fn bar_message(input_reads: u64, output_reads: u64, bytes: u64, elapsed: Duration) -> String {
    // Reads consumed and segments emitted, so a filter discarding everything
    // shows while the run is going rather than only in the summary. Labeled
    // `out` rather than `kept` because a split read emits several segments, so
    // the second figure can exceed the first. The two counters are also sampled
    // a moment apart, so the pair drifts until the run settles.
    let mut s = format!(
        "{} reads, {} out",
        human_count(input_reads),
        human_count(output_reads)
    );
    if bytes > 0 {
        let secs = elapsed.as_secs_f64().max(1e-3);
        let mbps = (bytes as f64 / 1_000_000.0) / secs;
        s.push_str(&format!(", {mbps:.0} MB/s"));
    }
    s
}

/// Line-mode periodic progress log, emitted at INFO every `log_interval` (see
/// `resolve_log_interval`): `Processed 1,200,000 input reads, 42%, 45k reads/s,
/// 380 MB/s, ETA 00:00:40`. Fields, in order: full-precision input read count
/// (reads consumed, not reads emitted, which differ under `--qual-split`),
/// percent complete (if `total` bytes are known), reads/s, MB/s (if any bytes
/// have been read), ETA (if `total` is known).
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

    /// Dropping an active line-mode handle joins its ticker thread.
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
        assert!(h.ticker.is_some(), "start() spawns a ticker");
        drop(h); // joins the ticker thread and returns
    }

    /// Dropping an active bar-mode handle joins its ticker and clears the bar.
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
        assert!(h.ticker.is_some(), "start() spawns a ticker");
        assert!(h.bar.is_some(), "Mode::Bar creates a live bar");
        drop(h); // joins the ticker thread and clears the bar
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
                "Non-color output carries no escape bytes: {s:?}"
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
            assert!(
                s.contains('\x1b'),
                "Color output carries escape bytes: {s:?}"
            );
            assert!(
                s.ends_with("\x1b[0m"),
                "Color output ends with a reset: {s:?}"
            );
        }
    }

    #[test]
    fn level_mapping() {
        assert_eq!(level_from(0, false), LevelFilter::INFO);
        assert_eq!(level_from(1, false), LevelFilter::DEBUG);
        assert_eq!(level_from(2, false), LevelFilter::TRACE);
        assert_eq!(level_from(5, false), LevelFilter::TRACE);
        // Quiet outranks verbosity.
        assert_eq!(level_from(0, true), LevelFilter::WARN);
        assert_eq!(level_from(3, true), LevelFilter::WARN);
    }

    /// `--progress` is independent of the log level: the summary is kept while
    /// the periodic line or the bar is suppressed.
    #[test]
    fn progress_mode_overrides_the_terminal_default() {
        // A terminal would otherwise get a bar.
        assert_eq!(
            select_mode(false, true, 0, false, ProgressMode::None),
            Mode::Silent
        );
        assert_eq!(
            select_mode(false, true, 0, false, ProgressMode::Plain),
            Mode::Line
        );
        // Redirected output would otherwise get periodic lines.
        assert_eq!(
            select_mode(false, false, 0, false, ProgressMode::Bar),
            Mode::Bar
        );
        assert_eq!(
            select_mode(false, false, 0, false, ProgressMode::None),
            Mode::Silent
        );
    }

    /// `--quiet` drops the summary as well, so it outranks any progress choice.
    #[test]
    fn quiet_outranks_every_progress_mode() {
        for p in [
            ProgressMode::Auto,
            ProgressMode::Bar,
            ProgressMode::Plain,
            ProgressMode::None,
        ] {
            assert_eq!(select_mode(true, true, 0, false, p), Mode::Off);
        }
    }

    /// Silent keeps the multi-line banner and the summary; only the in-flight
    /// progress reporting is gone.
    #[test]
    fn silent_still_shows_the_banner() {
        let h = ProgressHandle {
            mode: Mode::Silent,
            multi: MultiProgress::new(),
            bar: None,
            ticker: None,
            start: None,
            log_interval: Duration::from_secs(30),
        };
        assert!(h.shows_lines());
        assert!(!h.is_bar());
    }

    #[test]
    fn select_mode_quiet_always_off() {
        // Quiet outranks TTY state, verbosity, and WHITTLE_LOG.
        assert_eq!(
            select_mode(true, true, 0, false, ProgressMode::Auto),
            Mode::Off
        );
        assert_eq!(
            select_mode(true, false, 2, true, ProgressMode::Auto),
            Mode::Off
        );
    }

    #[test]
    fn select_mode_default_tty_is_bar() {
        assert_eq!(
            select_mode(false, true, 0, false, ProgressMode::Auto),
            Mode::Bar
        );
    }

    #[test]
    fn select_mode_non_tty_is_always_line() {
        assert_eq!(
            select_mode(false, false, 0, false, ProgressMode::Auto),
            Mode::Line
        );
    }

    #[test]
    fn select_mode_verbose_tty_is_line() {
        assert_eq!(
            select_mode(false, true, 1, false, ProgressMode::Auto),
            Mode::Line
        );
    }

    #[test]
    fn select_mode_whittle_log_forces_line_even_on_a_bar_eligible_tty() {
        // A non-empty `WHITTLE_LOG` forces line mode even at the default
        // verbosity on a TTY; its debug and trace lines would otherwise
        // interleave with a live bar.
        assert_eq!(
            select_mode(false, true, 0, true, ProgressMode::Auto),
            Mode::Line
        );
    }

    /// An explicit bar gives way to lines under `-v` or `WHITTLE_LOG`: bar mode
    /// hides the multi-line banner that a verbose run is asking for, and debug
    /// lines would interleave with the bar.
    #[test]
    fn select_mode_explicit_bar_downgrades_to_line_when_verbose() {
        assert_eq!(
            select_mode(false, true, 1, false, ProgressMode::Bar),
            Mode::Line
        );
        assert_eq!(
            select_mode(false, false, 2, false, ProgressMode::Bar),
            Mode::Line
        );
        assert_eq!(
            select_mode(false, true, 0, true, ProgressMode::Bar),
            Mode::Line
        );
        assert_eq!(
            select_mode(false, false, 0, false, ProgressMode::Bar),
            Mode::Bar
        );
    }

    /// An unparseable `WHITTLE_LOG` does not disable logging: the filter falls
    /// back to the verbosity level and the rejected directive is reported.
    #[test]
    fn log_filter_falls_back_and_warns_on_an_invalid_whittle_log() {
        let (filter, applied, advisory) = log_filter(Some("garbage=nope=1"), 0, false);
        assert!(
            !applied,
            "An invalid filter does not count as WHITTLE_LOG set"
        );
        assert_eq!(filter.to_string(), "info");
        let advisory = advisory.expect("An invalid WHITTLE_LOG raises an advisory");
        assert!(advisory.warn);
        assert!(
            advisory.message.contains("garbage=nope=1"),
            "The advisory names the directive: {}",
            advisory.message
        );
    }

    #[test]
    fn log_filter_applies_a_valid_whittle_log_unless_quiet() {
        let (filter, applied, advisory) = log_filter(Some("debug"), 0, false);
        assert!(applied);
        assert!(advisory.is_none());
        assert_eq!(filter.to_string(), "debug");

        let (filter, applied, advisory) = log_filter(Some("debug"), 0, true);
        assert!(!applied, "Quiet outranks WHITTLE_LOG");
        assert!(advisory.is_none());
        assert_eq!(filter.to_string(), "warn");
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
        // 1000 bytes rolls over to the next unit rather than staying `1000 B`.
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
        // Split reads can produce more output records than input records.
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

    /// Boundary values remain normalized after rounding.
    #[test]
    fn human_count_rolls_k_to_m_at_boundary() {
        assert_eq!(human_count(999_500), "1.0M");
        assert_eq!(human_count(999_499), "999k");
    }

    #[test]
    fn human_bytes_rolls_over_at_unit_boundary() {
        assert_eq!(human_bytes(999_999), "1.0 MB");
        assert_eq!(human_bytes(999_999_999), "1.0 GB");
        assert_eq!(human_bytes(999_499), "999 KB"); // immediately below the boundary; stays KB
    }

    #[test]
    fn human_bases_rolls_over_at_unit_boundary() {
        assert_eq!(human_bases(999_999_999), "1.0 Gbp");
        assert_eq!(human_bases(999_950), "1.0 Mbp");
    }

    #[test]
    fn human_dur_rolls_seconds_to_minutes_at_boundary() {
        assert_eq!(human_dur(Duration::from_millis(59_996)), "1m00s");
        assert_eq!(human_dur(Duration::from_millis(59_990)), "59.99s");
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
    fn segments_dropped_line_lists_only_nonzero_reasons_in_fixed_order() {
        let stats = Stats {
            segments_dropped_short: 2_100,
            segments_dropped_low_qual: 1_100,
            ..Default::default()
        };
        assert_eq!(
            segments_dropped_line(&stats).unwrap(),
            "Segments dropped: 3,200 (2,100 too short, 1,100 low quality)"
        );
    }

    #[test]
    fn segments_dropped_line_covers_every_reason_in_order() {
        let stats = Stats {
            segments_dropped_short: 1,
            segments_dropped_long: 2,
            segments_dropped_low_qual: 3,
            segments_dropped_high_qual: 4,
            segments_dropped_gc: 5,
            ..Default::default()
        };
        assert_eq!(
            segments_dropped_line(&stats).unwrap(),
            "Segments dropped: 15 (1 too short, 2 too long, 3 low quality, \
             4 high quality, 5 GC out of range)"
        );
    }

    #[test]
    fn segments_dropped_line_omitted_when_total_zero() {
        assert_eq!(segments_dropped_line(&Stats::default()), None);
    }

    #[test]
    fn trimmed_to_nothing_line_reports_read_level_count() {
        let stats = Stats {
            reads_trimmed_to_nothing: 42,
            ..Default::default()
        };
        assert_eq!(
            trimmed_to_nothing_line(&stats).unwrap(),
            "Trimmed to nothing: 42 input reads produced no segments at all"
        );
    }

    #[test]
    fn trimmed_to_nothing_line_omitted_when_zero() {
        assert_eq!(trimmed_to_nothing_line(&Stats::default()), None);
    }

    #[test]
    fn all_filtered_line_reports_read_level_count() {
        let stats = Stats {
            reads_all_filtered: 7,
            ..Default::default()
        };
        assert_eq!(
            all_filtered_line(&stats).unwrap(),
            "All segments filtered: 7 input reads had every produced segment filtered"
        );
    }

    #[test]
    fn all_filtered_line_omitted_when_zero() {
        assert_eq!(all_filtered_line(&Stats::default()), None);
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
            "Zero bytes must not render a misleading MB/s field: {s}"
        );
        assert!(
            !s.contains("->") && !s.contains('\u{b7}'),
            "Plain ASCII only: {s}"
        );
    }

    #[test]
    fn periodic_line_with_total_adds_percent_and_eta() {
        let s = periodic_line(500, 42_000_000, Some(100_000_000), Duration::from_secs(2));
        assert!(s.contains("42%")); // 42 MB of 100 MB
        assert!(s.contains("MB/s"));
        assert!(s.contains("ETA"));
    }

    #[test]
    fn bar_message_without_bytes_omits_the_rate() {
        let s = bar_message(145_000, 140_000, 0, Duration::from_secs(60));
        assert_eq!(s, "145k reads, 140k out");
    }

    #[test]
    fn bar_message_with_bytes_adds_rate_but_never_a_total() {
        let s = bar_message(145_000, 145_000, 50_000_000, Duration::from_secs(60));
        assert!(s.starts_with("145k reads, 145k out"));
        assert!(s.contains("MB/s"));
        assert!(
            !s.contains(" of "),
            "No total-reads figure is invented: {s}"
        );
        assert!(
            !s.contains('%'),
            "The bar template draws the percentage: {s}"
        );
    }
}
