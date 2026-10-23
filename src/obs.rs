//! Observability: leveled logging (tracing) and periodic progress lines.

use tracing::level_filters::LevelFilter;

use crate::config::{Advisory, Config};
use tracing_subscriber::fmt::FormattedFields;

/// Maps the CLI verbosity and quiet flags to a tracing level. `WHITTLE_LOG`,
/// when set, is applied separately in `init` and overrides this unless `quiet`
/// is set; `quiet` yields WARN regardless of `WHITTLE_LOG`.
fn level_from(verbosity: u8, quiet: bool) -> LevelFilter {
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

use std::io::{self, IsTerminal};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use tracing::{Event, Level, Subscriber};
use tracing_subscriber::filter::Targets;
use tracing_subscriber::fmt;
use tracing_subscriber::fmt::format::Writer;
use tracing_subscriber::fmt::{FmtContext, FormatEvent, FormatFields};
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::registry::LookupSpan;
use tracing_subscriber::util::SubscriberInitExt;

use crate::workflow::{Counters, Stats};

/// How often the ticker thread wakes to check whether the periodic line is
/// due. The thread waits this long via `park_timeout`, but `stop_ticker` cuts
/// the wait short with `unpark`, so a finished run does not idle for it.
const TICK_INTERVAL: Duration = Duration::from_millis(250);

/// Returns the periodic-log cadence of `Mode::Line` for the effective log
/// level: 30 s at info, 10 s at debug and trace (`-v`, `-vv`, or a
/// `WHITTLE_LOG` filter that enables them). The ticker wakes at most every
/// `TICK_INTERVAL`, or sooner if `stop_ticker` unparks it, and logs only once
/// this cadence has elapsed.
fn log_interval(level: LevelFilter) -> Duration {
    if level >= LevelFilter::DEBUG {
        Duration::from_secs(10)
    } else {
        Duration::from_secs(30)
    }
}

/// Custom event formatter: `[YYYY-MM-DD HH:MM:SS] [LEVEL] Message`, replacing
/// the stock formatter's ` INFO`-padded, unbracketed level. The timestamp is the
/// local wall clock via `jiff`. `color` (set once in `init`: stderr is a
/// terminal and `NO_COLOR` is unset) gates ANSI coloring of the `[LEVEL]` token
/// only, so a redirected run carries no escape bytes.
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

/// The output mode for a run, computed once in `init` from `quiet`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    /// `--quiet`: warnings and errors only. No banner, no progress line, no
    /// summary.
    Off,
    /// The default: the multi-line start banner, a periodic progress line
    /// every `log_interval`, debug and trace output (per level), and the
    /// summary.
    Line,
}

/// Owns the live ticker thread that logs the periodic progress line, and the
/// run clock. Created in the binary.
pub struct ProgressHandle {
    /// The output mode selected by `init`.
    mode: Mode,
    /// The ticker thread, its stop flag, and its join handle, while live. The
    /// join handle also serves `stop_ticker`'s `unpark` call that wakes the
    /// thread out of its `park_timeout` wait.
    ticker: Option<(Arc<AtomicBool>, JoinHandle<()>)>,
    /// Wall-clock start: construction time until `start` resets it once the
    /// input is open; `finish` reports the elapsed time from it.
    start: Instant,
    /// `Mode::Line` periodic-log cadence, resolved once in `init` from the
    /// verbosity. Unused outside `Mode::Line`.
    log_interval: Duration,
}

impl ProgressHandle {
    /// Returns a no-op handle for tests and library callers that install no
    /// subscriber.
    pub fn disabled() -> Self {
        ProgressHandle {
            mode: Mode::Off,
            ticker: None,
            start: Instant::now(),
            log_interval: Duration::from_secs(30),
        }
    }

    /// True in `Mode::Line`, which prints the multi-line banner. False in
    /// `Mode::Off`, where INFO output is filtered.
    pub fn shows_lines(&self) -> bool {
        self.mode == Mode::Line
    }

    /// Begins live progress once the input is open. In `Mode::Line` a ticker
    /// thread polls the shared counters every `TICK_INTERVAL` and logs a
    /// periodic INFO line every `log_interval`; `total` adds the percentage
    /// and ETA when the input size is known. In `Mode::Off` only the run clock
    /// starts. A ticker left running by an earlier `start` is stopped first.
    pub fn start(&mut self, total: Option<u64>, counters: Arc<Counters>) {
        // The timer runs in every mode, including `Off`: `--quiet` silences the
        // human-readable summary, but `elapsed_seconds` in `--summary-json` is
        // still reported.
        let start = Instant::now();
        self.start = start;
        if self.mode != Mode::Line {
            return;
        }
        self.stop_ticker();

        let stop = Arc::new(AtomicBool::new(false));
        let stop_t = stop.clone();
        let log_interval = self.log_interval;
        let handle = std::thread::spawn(move || {
            let mut last_log = start;
            while !stop_t.load(Ordering::Relaxed) {
                std::thread::park_timeout(TICK_INTERVAL);
                if last_log.elapsed() >= log_interval {
                    let ir = counters.input_reads.load(Ordering::Relaxed);
                    let or = counters.output_reads.load(Ordering::Relaxed);
                    let by = counters.bytes_read.load(Ordering::Relaxed);
                    tracing::info!("{}", periodic_line(ir, or, by, total, start.elapsed()));
                    last_log = Instant::now();
                }
            }
        });
        self.ticker = Some((stop, handle));
    }

    /// Stops the ticker (signal, unpark, and join), if it is live. The unpark
    /// cuts short whatever `park_timeout` wait the thread is in, so the join
    /// returns immediately instead of waiting for the next tick. Idempotent:
    /// the field is taken, so a second call (including the one in `Drop` after
    /// an explicit `finish`) is a no-op. Shared by `finish`, which follows it
    /// with the end-of-run summary, and `Drop`, which cleans up silently on an
    /// early error return.
    fn stop_ticker(&mut self) {
        if let Some((stop, handle)) = self.ticker.take() {
            stop.store(true, Ordering::Relaxed);
            handle.thread().unpark();
            // The default panic hook has already printed a ticker panic; the
            // assertion documents that a clean join is the only expected outcome.
            let joined = handle.join();
            debug_assert!(joined.is_ok(), "Progress ticker thread panicked");
        }
    }

    /// Stops the ticker, then logs the end-of-run summary. Stopping happens
    /// first so no periodic line interleaves with the summary. The closing
    /// `Completed` line is separate (`complete`), so the caller can write its
    /// artifacts between the two and a failed write is never reported after a
    /// success line.
    ///
    /// Returns the elapsed duration it reported, so `summary::Summary` quotes
    /// the same number.
    pub fn finish(&mut self, stats: &Stats, cfg: &Config) -> Duration {
        // Elapsed is taken before `stop_ticker`, so the join of the ticker
        // thread is not charged to the run.
        let elapsed = self.start.elapsed();
        self.stop_ticker();

        tracing::info!("{}", summary_line(stats, elapsed));

        if let Some(line) = bases_line(stats) {
            tracing::info!("{}", line);
        }

        if let Some(line) = tag_filtered_line(stats) {
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

        let key_names = cfg.splitter.as_deref().map(|s| s.keys.names.as_slice());
        for line in target_lines(stats, key_names) {
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
        if stats.footprint_tags_dropped_reads > 0 {
            tracing::warn!(
                reads = stats.footprint_tags_dropped_reads,
                "fibertools footprints (ns/nl/as/al) were removed from trimmed or split reads"
            );
        }
        if stats.barcode_tag_malformed_reads > 0 {
            tracing::warn!(
                reads = stats.barcode_tag_malformed_reads,
                "Barcode positions (bi) that did not describe a window inside the read were \
                 left untrimmed"
            );
        }
        if stats.barcode_tag_unverified_reads > 0 {
            tracing::warn!(
                reads = stats.barcode_tag_unverified_reads,
                "Barcode positions (bi) at which no barcode sequence was found were left \
                 untrimmed"
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
    /// `finish` returned.
    pub fn complete(&self, elapsed: Duration, output: &str) {
        tracing::info!("{}", completed_line(elapsed, output));
    }
}

/// Backstop for an early error return between `start` and `finish`: without it
/// the ticker thread keeps running past the error and can log a progress line
/// after the fatal message `main` prints. Stops the ticker without a summary,
/// since an error path has no `Stats`. A no-op after an explicit `finish`.
impl Drop for ProgressHandle {
    fn drop(&mut self) {
        self.stop_ticker();
    }
}

/// Resolves the level filter from `--quiet`, `WHITTLE_LOG`, and the verbosity.
///
/// `--quiet` always wins (WARN); otherwise a non-empty `WHITTLE_LOG` overrides
/// `-v`/`-vv`; otherwise the level follows `verbosity`. A `WHITTLE_LOG` that does
/// not parse falls back to the verbosity level and returns an advisory naming it,
/// since a lossy parse would enable nothing and hide even the ERROR line of a
/// failing run. A `WHITTLE_LOG` that names only targets takes the verbosity
/// level for every other target, which it would otherwise disable.
fn log_filter(
    whittle_log: Option<&str>,
    verbosity: u8,
    quiet: bool,
) -> (Targets, Option<Advisory>) {
    let fallback = || Targets::new().with_default(level_from(verbosity, quiet));
    if quiet {
        return (fallback(), None);
    }
    let Some(spec) = whittle_log else {
        return (fallback(), None);
    };
    match spec.parse::<Targets>() {
        Ok(filter) if filter.default_level().is_none() => {
            (filter.with_default(level_from(verbosity, false)), None)
        },
        Ok(filter) => (filter, None),
        Err(e) => {
            let level = level_from(verbosity, false);
            let advisory = Advisory::warn(format!(
                "WHITTLE_LOG is not a valid log filter and is ignored: \
                 value={spec:?}, error={e}, level={level}"
            ));
            (fallback(), Some(advisory))
        },
    }
}

/// Installs the global subscriber and returns the progress handle. Called once,
/// in the binary.
///
/// The level follows `log_filter`; a rejected `WHITTLE_LOG` lands in
/// `cfg.advisories`, which `run` prints once the banner is up. Every line goes
/// to stderr. The `[LEVEL]` token is colored only when stderr is a terminal
/// and `NO_COLOR` is unset. `--quiet` selects `Mode::Off`; every other run
/// logs the periodic progress line at the cadence of the effective level.
pub fn init(cfg: &mut Config) -> ProgressHandle {
    let (verbosity, quiet) = (cfg.verbosity, cfg.quiet);
    let whittle_log = std::env::var("WHITTLE_LOG").ok().filter(|s| !s.is_empty());
    let (filter, advisory) = log_filter(whittle_log.as_deref(), verbosity, quiet);
    cfg.advisories.extend(advisory);
    let color = io::stderr().is_terminal() && std::env::var_os("NO_COLOR").is_none();
    tracing_subscriber::registry()
        .with(filter)
        .with(
            fmt::layer()
                .event_format(WhittleFormat { color })
                .with_ansi(false)
                .with_writer(io::stderr),
        )
        .init();
    ProgressHandle {
        mode: if quiet { Mode::Off } else { Mode::Line },
        ticker: None,
        start: Instant::now(),
        log_interval: log_interval(LevelFilter::current()),
    }
}

/// Compact magnitude for live progress fields: `750`, `145k`, `1.2M`.
fn human_count(n: u64) -> String {
    // A value that would round to 1000k is promoted so the result stays
    // normalized.
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
/// 2.00s`. It carries no kept percentage, since quality splitting can turn one
/// input read into several segments and a read-count percentage would exceed
/// 100%.
fn summary_line(stats: &Stats, elapsed: Duration) -> String {
    format!(
        "Summary: {} input reads, {} output reads in {}",
        commas(stats.input_reads),
        commas(stats.output_reads),
        human_dur(elapsed)
    )
}

/// Human-readable base count for the yield summary's `Bases:` line: `12.4 Gbp`,
/// `460.0 Mbp`, `8.2 kbp`, `500 bp`. Uses decimal units and one decimal place
/// for kbp and larger values.
fn human_bases(n: u64) -> String {
    // A value that would round to 1000.0 in its unit is promoted to the next.
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

/// The end-of-run tag-filter line: `Tag filtered: 1,234 input reads did not
/// satisfy --tag-filter`. `None` when no read was rejected by the filter.
fn tag_filtered_line(stats: &Stats) -> Option<String> {
    if stats.reads_tag_filtered == 0 {
        return None;
    }
    Some(format!(
        "Tag filtered: {} input reads did not satisfy --tag-filter",
        commas(stats.reads_tag_filtered)
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
/// too long, low quality, high quality, too many expected errors, GC out of
/// range. Counts segments, not
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
        (
            DropReason::ExpectedErrors,
            stats.segments_dropped_expected_errors,
        ),
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

/// The end-of-run `--targets` table: one line per sheet key (`Target <key>: 120
/// segments, 45.2 kbp (60.0% of assigned)`), then the unassigned reasons on one
/// line, then ambiguous and discarded totals. Percentages are of the total
/// assigned reads across every key, `0.0%` when that total is zero. Empty
/// without `--targets` (`key_names` is `None`) or before `stats` carries a
/// split snapshot. `key_names` is `Keys::names`, so it and `split.keys` are
/// the same length and order.
fn target_lines(stats: &Stats, key_names: Option<&[String]>) -> Vec<String> {
    use crate::split::Unassigned;

    let (Some(key_names), Some(split)) = (key_names, stats.split.as_ref()) else {
        return Vec::new();
    };
    let assigned_total: u64 = split.keys.iter().map(|k| k.reads).sum();
    let mut lines: Vec<String> = key_names
        .iter()
        .zip(&split.keys)
        .map(|(name, k)| {
            let pct = if assigned_total == 0 {
                0.0
            } else {
                100.0 * k.reads as f64 / assigned_total as f64
            };
            format!(
                "Target {name}: {} segments, {} ({pct:.1}% of assigned)",
                commas(k.reads),
                human_bases(k.bases),
            )
        })
        .collect();
    let unassigned_total: u64 = split.unassigned.iter().sum();
    // `SplitStats::unassigned` counts the reasons in this order.
    let reasons = [
        Unassigned::NoMatch,
        Unassigned::MissingEnd,
        Unassigned::Orientation,
        Unassigned::Length,
    ];
    let by_reason: Vec<String> = reasons
        .iter()
        .zip(&split.unassigned)
        .map(|(reason, n)| format!("{} {}", reason.label(), commas(*n)))
        .collect();
    lines.push(format!(
        "Target unassigned: {} segments ({})",
        commas(unassigned_total),
        by_reason.join(", ")
    ));
    lines.push(format!(
        "Target ambiguous: {} segments",
        commas(split.ambiguous)
    ));
    lines.push(format!(
        "Target discarded: {} segments",
        commas(split.discarded)
    ));
    lines
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
        // Rounded to the nearest second, so 59.996s reads "1m00s".
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

/// `HH:MM:SS` duration, rounded to the second, for the periodic line's ETA
/// field.
fn fmt_hms(d: Duration) -> String {
    let s = d.as_secs_f64().round() as u64;
    format!("{:02}:{:02}:{:02}", s / 3600, (s % 3600) / 60, s % 60)
}

/// The periodic progress line, emitted at INFO every `log_interval`:
/// `Processed 1,200,000 input reads, 1,150,000 output reads, 42%, 45k
/// reads/s, 380 MB/s, ETA 00:00:40`. Fields, in order: input reads consumed,
/// output reads written (a split read writes several, and a filter that
/// discards every read shows as a zero here while the run is in flight),
/// percent complete (if `total` bytes are known), reads/s, MB/s (if any bytes
/// have been read), ETA (if `total` is known). The two counts are sampled a
/// moment apart, so their ratio is approximate until the run ends.
fn periodic_line(
    input_reads: u64,
    output_reads: u64,
    bytes: u64,
    total: Option<u64>,
    elapsed: Duration,
) -> String {
    let secs = elapsed.as_secs_f64().max(1e-3);
    let rps = input_reads as f64 / secs;

    let mut s = format!(
        "Processed {} input reads, {} output reads",
        commas(input_reads),
        commas(output_reads)
    );
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

    /// A line-mode handle with a 30 s cadence, as `init` builds at info.
    fn line_handle() -> ProgressHandle {
        ProgressHandle {
            mode: Mode::Line,
            ticker: None,
            start: Instant::now(),
            log_interval: Duration::from_secs(30),
        }
    }

    /// Dropping an active line-mode handle joins its ticker thread.
    #[test]
    fn dropping_started_handle_stops_ticker_without_hanging() {
        let mut h = line_handle();
        h.start(None, Arc::new(Counters::default()));
        assert!(h.ticker.is_some(), "start() spawns a ticker");
        drop(h); // joins the ticker thread and returns
    }

    /// A second `start` without `finish` stops the first ticker rather than
    /// leaving its thread running.
    #[test]
    fn restarting_a_handle_stops_the_previous_ticker() {
        let mut h = line_handle();
        h.start(None, Arc::new(Counters::default()));
        let first = Arc::clone(&h.ticker.as_ref().expect("a ticker runs").0);
        h.start(None, Arc::new(Counters::default()));
        assert!(first.load(std::sync::atomic::Ordering::Relaxed));
        assert!(h.ticker.is_some());
    }

    /// Stopping a live ticker returns without waiting for its next
    /// `TICK_INTERVAL` tick, since `stop_ticker` unparks it. The bound
    /// (200 ms) stays below the tick interval (250 ms), so a stop that waits
    /// for the tick fails, and leaves room for scheduling delay.
    #[test]
    fn stop_ticker_wakes_immediately_instead_of_waiting_for_the_next_tick() {
        let mut h = line_handle();
        h.start(None, Arc::new(Counters::default()));
        let began = Instant::now();
        h.stop_ticker();
        let took = began.elapsed();
        assert!(
            took < Duration::from_millis(200),
            "stop_ticker took {took:?}, expected well under the {TICK_INTERVAL:?} tick interval"
        );
    }

    /// `--quiet` starts only the run clock: no ticker and no banner lines.
    #[test]
    fn off_mode_starts_no_ticker_and_shows_no_lines() {
        let mut h = ProgressHandle::disabled();
        h.start(Some(1_000), Arc::new(Counters::default()));
        assert!(h.ticker.is_none());
        assert!(!h.shows_lines());
        assert!(line_handle().shows_lines());
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

    /// An unparseable `WHITTLE_LOG` does not disable logging: the filter falls
    /// back to the verbosity level and the rejected directive is reported.
    #[test]
    fn log_filter_falls_back_and_warns_on_an_invalid_whittle_log() {
        let (filter, advisory) = log_filter(Some("garbage=nope=1"), 0, false);
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
        let (filter, advisory) = log_filter(Some("debug"), 0, false);
        assert!(advisory.is_none());
        assert_eq!(filter.to_string(), "debug");

        // Quiet outranks WHITTLE_LOG.
        let (filter, advisory) = log_filter(Some("debug"), 0, true);
        assert!(advisory.is_none());
        assert_eq!(filter.to_string(), "warn");
    }

    /// A `WHITTLE_LOG` that names only targets keeps the verbosity level for
    /// every other target, so errors, warnings and the banner still print; a
    /// bare level in the value stays the default.
    #[test]
    fn log_filter_gives_a_target_only_whittle_log_the_verbosity_default() {
        let (filter, advisory) = log_filter(Some("whittle::adapter=trace"), 0, false);
        assert!(advisory.is_none());
        assert_eq!(filter.default_level(), Some(LevelFilter::INFO));
        assert!(
            filter
                .iter()
                .any(|(target, level)| target == "whittle::adapter" && level == LevelFilter::TRACE)
        );

        let (filter, _) = log_filter(Some("whittle::adapter=trace"), 1, false);
        assert_eq!(filter.default_level(), Some(LevelFilter::DEBUG));

        let (filter, _) = log_filter(Some("warn,whittle::adapter=trace"), 2, false);
        assert_eq!(filter.default_level(), Some(LevelFilter::WARN));
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

    /// The periodic-log cadence is 30 s at info and 10 s at debug and trace.
    #[test]
    fn log_interval_is_30s_and_10s_when_verbose() {
        assert_eq!(log_interval(LevelFilter::WARN), Duration::from_secs(30));
        assert_eq!(log_interval(LevelFilter::INFO), Duration::from_secs(30));
        assert_eq!(log_interval(LevelFilter::DEBUG), Duration::from_secs(10));
        assert_eq!(log_interval(LevelFilter::TRACE), Duration::from_secs(10));
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
        let s = summary_line(&stats, Duration::from_secs(2));
        assert_eq!(s, "Summary: 1 input reads, 3 output reads in 2.00s");
        assert!(!s.contains('%'));
        assert!(!s.contains("Kept"));
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
            segments_dropped_expected_errors: 6,
            segments_dropped_gc: 5,
            ..Default::default()
        };
        assert_eq!(
            segments_dropped_line(&stats).unwrap(),
            "Segments dropped: 21 (1 too short, 2 too long, 3 low quality, \
             4 high quality, 6 too many expected errors, 5 GC out of range)"
        );
    }

    #[test]
    fn segments_dropped_line_omitted_when_total_zero() {
        assert_eq!(segments_dropped_line(&Stats::default()), None);
    }

    fn key_stats(reads: u64, bases: u64) -> crate::split::KeyStats {
        crate::split::KeyStats {
            reads,
            bases,
            both_ends: reads,
            five_only: 0,
            three_only: 0,
            plus: reads,
            minus: 0,
        }
    }

    #[test]
    fn target_lines_reports_one_line_per_key_then_unassigned_ambiguous_discarded() {
        let names = ["16S".to_string(), "ITS".to_string()];
        let stats = Stats {
            split: Some(crate::split::SplitStats {
                keys: vec![key_stats(30, 12_000), key_stats(10, 4_000)],
                unassigned: [5, 2, 1, 0],
                ambiguous: 3,
                discarded: 4,
            }),
            ..Default::default()
        };
        let lines = target_lines(&stats, Some(&names));
        assert_eq!(
            lines,
            vec![
                "Target 16S: 30 segments, 12.0 kbp (75.0% of assigned)",
                "Target ITS: 10 segments, 4.0 kbp (25.0% of assigned)",
                "Target unassigned: 8 segments (no_match 5, missing_end 2, orientation 1, length 0)",
                "Target ambiguous: 3 segments",
                "Target discarded: 4 segments",
            ]
        );
    }

    /// A key with no reads reports `0.0%`, not a NaN from dividing by an
    /// assigned total of zero.
    #[test]
    fn target_lines_reports_zero_percent_when_nothing_was_assigned() {
        let names = ["16S".to_string()];
        let stats = Stats {
            split: Some(crate::split::SplitStats {
                keys: vec![key_stats(0, 0)],
                unassigned: [7, 0, 0, 0],
                ambiguous: 0,
                discarded: 0,
            }),
            ..Default::default()
        };
        let lines = target_lines(&stats, Some(&names));
        assert_eq!(lines[0], "Target 16S: 0 segments, 0 bp (0.0% of assigned)");
    }

    #[test]
    fn target_lines_empty_without_targets() {
        assert_eq!(target_lines(&Stats::default(), None), Vec::<String>::new());
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
        let s = periodic_line(1_200_000, 1_150_000, 0, None, Duration::from_secs(10));
        assert!(s.contains("Processed 1,200,000 input reads, 1,150,000 output reads"));
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
        let s = periodic_line(
            500,
            450,
            42_000_000,
            Some(100_000_000),
            Duration::from_secs(2),
        );
        assert!(s.contains("42%")); // 42 MB of 100 MB
        assert!(s.contains("MB/s"));
        assert!(s.contains("ETA 00:00:03"), "{s}");
    }
}
