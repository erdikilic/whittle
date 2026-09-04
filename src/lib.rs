//! whittle: a tag-aware trimmer for long-read FASTQ and unaligned BAM.
//!
//! `run` is the library entry point; `cli::parse` builds its `Config`.

pub mod adapter;
pub(crate) mod banner;
pub mod cli;
pub mod config;
pub mod filter;
pub(crate) mod guards;
pub mod io;
pub mod mods;
pub mod obs;
pub mod qual;
pub mod record;
pub mod summary;
pub mod trim;
pub mod workflow;

pub use banner::command_line;

use std::borrow::Cow;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::Context as _;

use config::AdapterInfer;
pub use config::Config;
use io::Format;

/// Runs one trimming job. A directory input selects folder-merge mode (every
/// read file in it is merged into one output); any other input is a single file
/// or stdin. FASTQ and unaligned BAM are supported.
///
/// `obs` drives progress and end-of-run output; library callers pass
/// `obs::ProgressHandle::disabled()`.
pub fn run(cfg: Config, obs: &mut obs::ProgressHandle) -> anyhow::Result<()> {
    let mut cfg = cfg;
    let dir = cfg
        .io
        .input
        .as_deref()
        .filter(|p| p.is_dir())
        .map(Path::to_path_buf);
    let result = match dir {
        Some(dir) => run_folder(&dir, &mut cfg, obs),
        None => run_single(&mut cfg, obs),
    };
    // `announce` drains the parse-time advisories as it prints them. A setup
    // failure ahead of the banner leaves them pending, and they print here so an
    // early error never hides a skipped FASTA entry or a rejected `WHITTLE_LOG`.
    emit_advisories(&mut cfg.advisories);
    result
}

/// Runs single-file (or stdin) mode: detects the format, announces, and
/// dispatches.
fn run_single(cfg: &mut Config, obs: &mut obs::ProgressHandle) -> anyhow::Result<()> {
    let setup_start = Instant::now();
    guards::guard_output_collisions(cfg, &[])?;

    let in_path = cfg.io.input.clone();
    let in_path = in_path.as_deref();

    // The total input byte count, known for a file, drives a determinate
    // progress bar with percent and ETA; stdin has no size, so `None` selects a
    // spinner (see `obs::ProgressHandle::start`).
    let total: Option<u64> = in_path
        .and_then(|p| std::fs::metadata(p).ok())
        .map(|m| m.len());

    // Created before the reader so one `Arc` is shared by `CountingReader`, the
    // workflow call, and `obs.start`.
    let counters = Arc::new(workflow::Counters::default());

    // The input (file or stdin) is opened before format detection so the sniffed
    // bytes are kept: whatever detection consumes is chained back in front of
    // the stream through a `Cursor` before the reader is built. `CountingReader`
    // sits innermost, so sniffed bytes are counted once.
    let raw: Box<dyn Read + Send> = match in_path {
        Some(p) => Box::new(
            std::fs::File::open(p).with_context(|| format!("opening input {}", p.display()))?,
        ),
        None => Box::new(std::io::stdin()),
    };
    let raw: Box<dyn Read + Send> =
        Box::new(io::counting::CountingReader::new(raw, counters.clone()));
    let mut source: Box<dyn Read + Send> = Box::new(BufReader::new(raw));

    let in_fmt = match cfg.io.in_format {
        Some(f) => f,
        None => match in_path.and_then(io::from_extension) {
            // A `.gz` extension covers both plain gzip and BGZF; the first block
            // header decides, and BGZF gets a block-parallel decode share.
            Some(Format::FastqGz) => io::probe_gz(&mut source)?,
            Some(f) => f,
            None => {
                let (fmt, replayed) = detect_format(in_path, source)?;
                source = replayed;
                fmt
            },
        },
    };

    // Advisory only: an explicit `--in-format` or `--out-format` decides the
    // format, and a disagreement with the path's extension (`--out-format fastq`
    // on an `out.fastq.gz` path) is reported as a warning. Skipped for stdin,
    // stdout, and paths without an extension. Only the check runs here; both
    // warnings are logged after the banner with the other advisories.
    let mismatch_warn = io::format_mismatch_warning("--in-format", cfg.io.in_format, in_path);
    let out_mismatch_warn =
        io::format_mismatch_warning("--out-format", cfg.io.out_format, cfg.io.output.as_deref());

    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), in_fmt));

    // Hard error before any writer or output file is created: BAM or gzipped
    // bytes are refused on an interactive terminal.
    guards::guard_stdout_binary(cfg, out_fmt)?;
    // Detection has classified the stream, so a piped or extensionless FASTQ
    // reaches the same refusal `cli::parse` applies to a named one.
    guards::guard_barcode_input(cfg, in_fmt)?;
    guards::guard_remove_tag_input(cfg, in_fmt)?;

    tracing::debug!(
        stage = "setup",
        format = in_fmt.label(),
        elapsed_ms = setup_start.elapsed().as_millis() as u64,
        "Input format detected"
    );

    let budget = plan_budget(
        cfg,
        in_fmt,
        out_fmt,
        matches!(in_fmt, Format::FastqBgzf | Format::Bam),
    );

    let mut warnings: Vec<String> = Vec::new();
    warnings.extend(mismatch_warn);
    warnings.extend(out_mismatch_warn);
    if is_no_op(cfg, in_fmt == out_fmt) {
        warnings.push(NO_OP_WARNING.to_string());
    }
    announce(
        cfg,
        obs,
        counters.clone(),
        &Startup {
            input_line: match (in_path, total) {
                (Some(p), Some(size)) => {
                    format!("Input: {} ({})", p.display(), obs::human_bytes(size))
                },
                (Some(p), None) => format!("Input: {}", p.display()),
                (None, _) => "Input: <stdin>".to_string(),
            },
            total,
            in_fmt,
            out_fmt,
            budget,
            warnings,
        },
    );

    let session = Session::begin(cfg, budget, counters, in_fmt, out_fmt);
    session.dispatch(cfg, obs, in_fmt, out_fmt, Source::Stream(source))
}

/// Runs folder-merge mode (`-i <dir>`): classifies the directory into one
/// format family, then merges its read files into a single trimmed output
/// through the same dispatch as the single-file path.
fn run_folder(dir: &Path, cfg: &mut Config, obs: &mut obs::ProgressHandle) -> anyhow::Result<()> {
    // `--in-format` is inert for a directory, whose family is decided per file
    // by extension, so a warning is queued below.
    let folder_in_format_ignored = cfg.io.in_format.is_some();

    // `classify` hard-errors when `-o` names a read file inside `-i <dir>`:
    // overwriting real input or stale prior output while merging the rest is
    // silent data loss. The guard then checks every other written file against
    // every member.
    let (family, paths) = io::dir::classify(dir, cfg.io.output.as_deref())?;
    guards::guard_output_collisions(cfg, &paths)?;
    let family_fmt = match family {
        io::dir::Family::Fastq => Format::Fastq,
        io::dir::Family::Bam => Format::Bam,
    };
    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), family_fmt));
    guards::guard_stdout_binary(cfg, out_fmt)?;
    guards::guard_barcode_input(cfg, family_fmt)?;
    guards::guard_remove_tag_input(cfg, family_fmt)?;

    // A `.gz` member is BGZF when its first block header says so.
    let bgzf_input = family_fmt == Format::Bam
        || paths.iter().any(|p| match io::from_extension(p) {
            Some(Format::FastqBgzf) => true,
            Some(Format::FastqGz) => io::is_bgzf_file(p),
            _ => false,
        });
    let budget = plan_budget(cfg, family_fmt, out_fmt, bgzf_input);
    let counters = Arc::new(workflow::Counters::default());

    // Summed unconditionally, not only when the banner prints: it also drives
    // the determinate progress bar.
    let total_bytes: u64 = paths
        .iter()
        .filter_map(|p| std::fs::metadata(p).ok())
        .map(|m| m.len())
        .sum();
    tracing::debug!(
        stage = "setup",
        format = family_fmt.label(),
        files = paths.len(),
        "Folder classified"
    );

    let mut warnings: Vec<String> = Vec::new();
    if folder_in_format_ignored {
        warnings.push(
            "--in-format is ignored for a directory input; folder files are classified by \
             extension per file"
                .to_string(),
        );
    }
    warnings.extend(io::format_mismatch_warning(
        "--out-format",
        cfg.io.out_format,
        cfg.io.output.as_deref(),
    ));
    // Every member file's precise format, not the collapsed family: a folder of
    // `.fastq.gz` merged into plain `.fastq` is a conversion.
    let same_format = paths.iter().all(|p| io::from_extension(p) == Some(out_fmt));
    if is_no_op(cfg, same_format) {
        warnings.push(NO_OP_WARNING.to_string());
    }
    announce(
        cfg,
        obs,
        counters.clone(),
        &Startup {
            input_line: format!(
                "Input: {} {} files, {}",
                paths.len(),
                family_fmt.label(),
                obs::human_bytes(total_bytes)
            ),
            total: Some(total_bytes),
            in_fmt: family_fmt,
            out_fmt,
            budget,
            warnings,
        },
    );

    let session = Session::begin(cfg, budget, counters, family_fmt, out_fmt);
    session.dispatch(cfg, obs, family_fmt, out_fmt, Source::Folder(paths))
}

/// Sniffs the input format from the stream's first bytes when neither
/// `--in-format` nor the path extension decides it. Returns the format and the
/// stream with the consumed bytes chained back in front, so the reader built
/// next sees the input from its start.
fn detect_format(
    in_path: Option<&Path>,
    mut source: Box<dyn Read + Send>,
) -> anyhow::Result<(Format, Box<dyn Read + Send>)> {
    // The probe covers a full BGZF block header (18 bytes), which tells a BAM on
    // stdin or under an unknown extension apart from gzipped FASTQ. A single
    // `read()` may return fewer bytes, so the loop fills the buffer.
    let mut probe = [0u8; 18];
    let mut n = 0;
    while n < probe.len() {
        let r = source.read(&mut probe[n..])?;
        if r == 0 {
            break;
        }
        n += r;
    }
    let mut replay = probe[..n].to_vec();
    let fmt = if io::is_bgzf(&replay) {
        let block_size = usize::from(u16::from_le_bytes([replay[16], replay[17]])) + 1;
        if block_size < replay.len() {
            anyhow::bail!("invalid BGZF block size {block_size}");
        }
        replay.resize(block_size, 0);
        source.read_exact(&mut replay[n..]).with_context(|| {
            format!(
                "truncated BGZF block in {}",
                in_path.map_or_else(|| "<stdin>".to_string(), |p| p.display().to_string())
            )
        })?;
        io::detect_bgzf_block(&replay)?
    } else {
        io::detect_input(in_path, &replay)?
    };
    Ok((fmt, Box::new(std::io::Cursor::new(replay).chain(source))))
}

/// Splits the `-t` worker budget for a dispatch once, so the banner's
/// `Threads:` line and the workflow use the same numbers.
///
/// BGZF containers (BAM and bgzf FASTQ) decode block-parallel, but only when
/// render is light: adapter search (preset, FASTA, or inference) makes render
/// the bottleneck instead, and the decode share goes there.
fn plan_budget(
    cfg: &Config,
    in_fmt: Format,
    out_fmt: Format,
    bgzf_input: bool,
) -> config::ThreadBudget {
    let parallel_decode =
        bgzf_input && cfg.adapters.is_none() && cfg.adapter_infer == AdapterInfer::Off;
    config::thread_budget(
        cfg.threads,
        config::render_heavy_for(in_fmt, cfg),
        parallel_decode,
        config::encode_kind_for(out_fmt),
    )
}

/// Where the records come from: one stream (a file or stdin, with any sniffed
/// bytes chained back in), or the member files of a folder.
enum Source {
    /// A single byte stream.
    Stream(Box<dyn Read + Send>),
    /// The read files of a directory, in merge order.
    Folder(Vec<PathBuf>),
}

/// The per-run state every dispatch arm shares: the thread split, the live
/// counters, the output label, and the processing-phase clock. One `dispatch`
/// serves both entry points, so every arm settles the configuration and writes
/// the summary file.
struct Session {
    /// The per-stage worker split.
    budget: config::ThreadBudget,
    /// Live counters shared with the reader and the progress ticker.
    counters: Arc<workflow::Counters>,
    /// The output path or `<stdout>`, as shown in the banner and the closer.
    out_desc: String,
    /// Start of the processing phase. Stages run concurrently (read, trim, and
    /// write overlap across threads), so this marks a phase boundary, not a
    /// CPU-time split.
    t0: Instant,
}

impl Session {
    /// Opens the processing phase, after the banner has printed.
    fn begin(
        cfg: &Config,
        budget: config::ThreadBudget,
        counters: Arc<workflow::Counters>,
        in_fmt: Format,
        out_fmt: Format,
    ) -> Self {
        tracing::debug!(
            stage = "dispatch",
            input = in_fmt.label(),
            output = out_fmt.label(),
            threads = cfg.threads,
            decode = budget.decode,
            render = budget.render,
            encode = budget.encode,
            "Processing started"
        );
        Session {
            budget,
            counters,
            out_desc: banner::output_desc(cfg.io.output.as_deref()),
            t0: Instant::now(),
        }
    }

    /// Runs the workflow for one (input, output) format pair.
    ///
    /// Every arm settles the adapter set before creating the output file: an
    /// inference-report exit (`Ok(None)` from `settle`) returns before any output
    /// is touched, since building the writer first would truncate a pre-existing
    /// `-o` file even though report mode writes no records. The FASTQ-to-BAM
    /// rejection sits before its writer for the same reason, so a rejected run
    /// leaves no 0-byte file behind.
    fn dispatch(
        self,
        cfg: &mut Config,
        obs: &mut obs::ProgressHandle,
        in_fmt: Format,
        out_fmt: Format,
        source: Source,
    ) -> anyhow::Result<()> {
        match (in_fmt, out_fmt) {
            (Format::Bam, Format::Bam) => {
                note_tags_ignored(cfg, in_fmt, out_fmt);
                // Only the first file's header is written, so differing read
                // groups in the other files are reported for BAM output.
                if let Source::Folder(paths) = &source {
                    io::dir::warn_on_bam_header_mismatch(paths);
                }
                let (header, records) = self.bam_reader(source)?;
                let Some(records) = settle(records, cfg, self.budget, adapter::resolve::bam_seq)?
                else {
                    return Ok(());
                };
                let out_header = io::bam::provenance_header(
                    header,
                    cfg.threads <= 1 || cfg.ordered,
                    &command_line(std::env::args_os()),
                );
                let mut sink = io::bam::writer(
                    cfg.io.output.as_deref(),
                    &out_header,
                    self.budget.encode,
                    cfg.compression_level,
                )?;
                let stats =
                    workflow::run_raw_bam(&out_header, records, &mut sink, cfg, &self.counters)?;
                // Explicit finish (final bgzf block and EOF marker) rather than
                // `Drop`, which discards a `try_finish` error: an I/O failure on
                // the final flush (ENOSPC) would otherwise yield a truncated BAM
                // and a success exit code.
                sink.finish()?;
                self.finish(obs, &stats, cfg)
            },
            (Format::Bam, Format::Fastq | Format::FastqGz | Format::FastqBgzf) => {
                let (_header, records) = self.bam_reader(source)?;
                let Some(records) = settle(records, cfg, self.budget, adapter::resolve::bam_seq)?
                else {
                    return Ok(());
                };
                let mut writer = io::fastq::writer(cfg, out_fmt, self.budget.encode)?;
                let stats = workflow::run_bam_to_fastq(records, &mut writer, cfg, &self.counters)?;
                writer.finish()?;
                self.finish(obs, &stats, cfg)
            },
            (Format::Fastq | Format::FastqGz | Format::FastqBgzf, Format::Bam) => {
                anyhow::bail!("FASTQ-to-BAM conversion is not supported")
            },
            (Format::Fastq | Format::FastqGz | Format::FastqBgzf, _) => {
                note_tags_ignored(cfg, in_fmt, out_fmt);
                let records = self.fastq_reader(source, in_fmt)?;
                let Some(records) = settle(records, cfg, self.budget, |r| {
                    Cow::Borrowed(r.seq.as_slice())
                })?
                else {
                    return Ok(());
                };
                let mut writer = io::fastq::writer(cfg, out_fmt, self.budget.encode)?;
                let stats = workflow::run_fastq(records, &mut writer, cfg, &self.counters)?;
                writer.finish()?;
                self.finish(obs, &stats, cfg)
            },
        }
    }

    /// Opens the BAM record stream. A stdin BAM's sniffed bytes are chained back
    /// into the stream, so it is read as is rather than reopened.
    fn bam_reader(
        &self,
        source: Source,
    ) -> anyhow::Result<(noodles_sam::Header, io::bam::RawRecordIter)> {
        match source {
            Source::Stream(src) => io::bam::reader_from(src, self.budget.decode),
            Source::Folder(paths) => {
                io::dir::bam_reader(&paths, self.budget.decode, self.counters.clone())
            },
        }
    }

    /// Opens the FASTQ-family record stream for `in_fmt`.
    fn fastq_reader(
        &self,
        source: Source,
        in_fmt: Format,
    ) -> anyhow::Result<Box<dyn Iterator<Item = anyhow::Result<record::ReadRecord>> + Send>> {
        match source {
            Source::Stream(src) => match in_fmt {
                Format::Fastq => Ok(io::fastq::reader_from(src, false)),
                Format::FastqGz => Ok(io::fastq::reader_from(src, true)),
                Format::FastqBgzf => io::fastq::reader_from_bgzf(src, self.budget.decode),
                Format::Bam => unreachable!("BAM input is dispatched to bam_reader"),
            },
            Source::Folder(paths) => Ok(io::dir::fastq_records(
                &paths,
                self.budget.decode,
                self.counters.clone(),
            )),
        }
    }

    /// Closes out a finished dispatch: logs the phase duration and the
    /// end-of-run summary, writes `--summary-json`, then logs the closer. The
    /// artifact precedes the `Completed` line so a failed write is never
    /// reported after a success line.
    fn finish(
        &self,
        obs: &mut obs::ProgressHandle,
        stats: &workflow::Stats,
        cfg: &Config,
    ) -> anyhow::Result<()> {
        tracing::debug!(
            stage = "dispatch",
            elapsed_ms = self.t0.elapsed().as_millis() as u64,
            reads_in = stats.input_reads,
            reads_out = stats.output_reads,
            bases_in = stats.input_bases,
            bases_out = stats.output_bases,
            "Processing finished"
        );
        let elapsed = obs.finish(stats);
        if let Some(path) = cfg.summary_json.as_deref() {
            summary::Summary::new(
                cfg,
                stats,
                command_line(std::env::args_os()),
                self.out_desc.clone(),
                elapsed,
            )
            .write(path)?;
        }
        obs.complete(elapsed, &self.out_desc);
        Ok(())
    }
}

/// Settles the configuration for dispatch and returns the record stream.
///
/// Two fields are only knowable here. The adapter set depends on the reads:
/// presence detection samples a prefix and inference discovers from one, so the
/// set is final only after the banner has printed the configured set. The
/// render-pool size comes from the thread budget. Both are assigned in this one
/// place so every dispatch arm sees the same narrowed set and pool size.
///
/// `Ok(None)` means the run is over without writing records: that is
/// `--adapter-infer report`, which prints the inferred FASTA and stops.
fn settle<R, I, F>(
    records: I,
    cfg: &mut Config,
    budget: config::ThreadBudget,
    seq_of: F,
) -> anyhow::Result<Option<Box<dyn Iterator<Item = anyhow::Result<R>> + Send>>>
where
    I: Iterator<Item = anyhow::Result<R>> + Send + 'static,
    R: Send + 'static,
    F: for<'a> Fn(&'a R) -> std::borrow::Cow<'a, [u8]>,
{
    let Some(resolved) = adapter::resolve::resolve(records, cfg, seq_of)? else {
        note_report_only_ignores(cfg);
        return Ok(None);
    };
    // Captured before the overwrite: the banner already printed this figure, and
    // the summary reports it alongside what resolution settled on.
    cfg.adapters_configured = cfg.adapters.as_ref().map(|a| a.adapters.len());
    cfg.adapters = resolved.adapters;
    cfg.render_workers = budget.render;
    Ok(Some(resolved.records))
}

/// Warns for every artifact flag that report-only inference leaves unwritten;
/// the run exits 0 without creating any of them.
fn note_report_only_ignores(cfg: &Config) {
    for (flag, _) in cfg.write_targets() {
        tracing::warn!("{flag} is ignored under --adapter-infer report, which writes no records");
    }
}

/// The resolved facts both entry points announce before processing starts.
struct Startup {
    /// The `Input:` line body: a path and size, `<stdin>`, or a folder's file
    /// count and total size.
    input_line: String,
    /// Input bytes, when known. Drives a determinate bar; `None` gives a spinner.
    total: Option<u64>,
    /// The input format, or a folder's format family.
    in_fmt: Format,
    /// The resolved output format.
    out_fmt: Format,
    /// The per-stage worker split shown in the `Threads:` line.
    budget: config::ThreadBudget,
    /// Advisories specific to this entry point, logged after the shared ones.
    warnings: Vec<String>,
}

/// Prints the startup banner and the deferred advisories, then starts progress.
///
/// The order is the output contract: `whittle {version}` and `Command:` (from
/// `main`), then the resolved config, then every advisory, then live progress,
/// so the run's provenance is at the top and nothing interleaves with the bar.
/// Both entry points call this function, so an advisory is emitted on both
/// paths. Starting progress here rather than at the call sites makes the
/// ordering a property of this function.
fn announce(
    cfg: &mut Config,
    obs: &mut obs::ProgressHandle,
    counters: Arc<workflow::Counters>,
    s: &Startup,
) {
    if obs.shows_lines() {
        tracing::info!("{}", banner::operation_line(s.in_fmt, s.out_fmt));
        tracing::info!("{}", s.input_line);
        tracing::info!(
            "{}",
            banner::output_banner_line(
                cfg.io.output.as_deref(),
                s.out_fmt,
                cfg.compression_level,
                s.budget.encode
            )
        );
        tracing::info!("{}", banner::threads_banner_line(cfg.threads, s.budget));
        tracing::info!("{}", banner::filters_and_trim_line(&cfg.filter, &cfg.trim));
        if let Some(line) = banner::adapter_banner_line(
            cfg.adapters.as_ref(),
            cfg.adapter_sample,
            cfg.adapter_infer,
        ) {
            tracing::info!("{line}");
        }
    } else if obs.is_bar() {
        // Bar mode prints one line so the live bar stays clean.
        tracing::info!(
            "{} ({} threads)",
            banner::operation_line(s.in_fmt, s.out_fmt),
            cfg.threads
        );
    }

    emit_advisories(&mut cfg.advisories);
    if let Some((requested, ncpu)) = cfg.threads_clamped {
        tracing::warn!(
            requested,
            ncpu,
            "Requested -t exceeds the CPU count; using the CPU count"
        );
    }
    for w in &s.warnings {
        tracing::warn!("{w}");
    }

    obs.start(s.total, counters);
}

/// Logs and drains the parse-time diagnostics held by `cli::parse` and
/// `obs::init` until the subscriber exists, so `--quiet` silences them and they
/// carry the standard prefix. Draining makes a second call a no-op, so `run`
/// prints only what the banner did not reach.
fn emit_advisories(advisories: &mut Vec<config::Advisory>) {
    for a in advisories.drain(..) {
        if a.warn {
            tracing::warn!("{}", a.message);
        } else {
            tracing::info!("{}", a.message);
        }
    }
}

/// True when the run neither trims nor filters, so the output mirrors the input
/// apart from format.
///
/// `same_format` is supplied by the caller because the two entry points answer
/// it differently: folder mode's family collapses `.fastq`, `.fastq.gz`, and
/// `.fastq.bgz` into one value, so comparing families there would classify a
/// decompression run as a no-op.
fn is_no_op(cfg: &Config, same_format: bool) -> bool {
    let no_trim = cfg.trim.head == 0 && cfg.trim.tail == 0 && cfg.trim.quality.is_none();
    let pass_through_filter = cfg.filter.min_length <= 1
        && cfg.filter.max_length == usize::MAX
        && cfg.filter.min_qual <= 0.0
        && cfg.filter.max_qual >= 1000.0
        && cfg.filter.min_gc.is_none()
        && cfg.filter.max_gc.is_none();
    no_trim
        && pass_through_filter
        && cfg.adapters.is_none()
        && !cfg.trim_barcodes
        && cfg.remove_tags.is_empty()
        && same_format
}

/// Warning for a run that neither trims nor filters.
const NO_OP_WARNING: &str =
    "No trimming or filtering options set; output will mostly mirror the input";

/// Warns when `--fastq-tags` was set to a non-default value (`none` or an
/// explicit list) on a path other than BAM-to-FASTQ, where it has no effect. An
/// explicit `all` equals the default and is silent.
fn note_tags_ignored(cfg: &Config, in_fmt: Format, out_fmt: Format) {
    if !matches!(cfg.fastq_tags, config::FastqTags::All) {
        tracing::warn!(
            input = in_fmt.label(),
            output = out_fmt.label(),
            "--fastq-tags applies only to BAM-to-FASTQ output and is ignored"
        );
    }
}

#[cfg(test)]
mod tests {
    fn base_filter() -> filter::FilterConfig {
        filter::FilterConfig {
            min_length: 1,
            max_length: usize::MAX,
            min_qual: 0.0,
            max_qual: 1000.0,
            min_gc: None,
            max_gc: None,
            qual_mode: qual::QualMode::Mean,
        }
    }

    fn base_trim() -> trim::TrimPlan {
        trim::TrimPlan {
            head: 0,
            tail: 0,
            quality: None,
        }
    }

    use super::*;

    #[test]
    fn encode_kind_for_maps_output_format() {
        assert_eq!(
            config::encode_kind_for(io::Format::Bam),
            config::EncodeKind::Bgzf
        );
        assert_eq!(
            config::encode_kind_for(io::Format::FastqGz),
            config::EncodeKind::Gzip
        );
        assert_eq!(
            config::encode_kind_for(io::Format::Fastq),
            config::EncodeKind::None
        );
    }

    fn a_read() -> crate::record::ReadRecord {
        crate::record::ReadRecord {
            name: b"r1".to_vec(),
            seq: b"ACGTACGTACGT".to_vec(),
            qual: vec![40; 12],
        }
    }

    /// `settle` is the one place both dispatch-time fields are written. A
    /// missing `render_workers` assignment would oversubscribe the render pool
    /// or run it single-threaded on the BAM full-window path.
    #[test]
    fn settle_sets_both_dispatch_fields() {
        let mut cfg = base_config();
        cfg.render_workers = 0;
        cfg.adapters = None;
        let budget = config::ThreadBudget {
            decode: 1,
            render: 5,
            encode: 2,
        };

        let records = vec![Ok(a_read())].into_iter();
        let got = settle(records, &mut cfg, budget, |r| {
            Cow::Borrowed(r.seq.as_slice())
        })
        .expect("Settle succeeds")
        .expect("Records are returned when not in report mode");

        assert_eq!(
            cfg.render_workers, 5,
            "Render pool size comes from the budget"
        );
        assert!(cfg.adapters.is_none(), "No adapters configured stays none");
        assert_eq!(got.count(), 1, "The record stream is handed back intact");
    }

    /// With no adapter work to do, the stream passes through untouched rather
    /// than being buffered.
    #[test]
    fn settle_passes_records_through_when_adapters_are_off() {
        let mut cfg = base_config();
        let budget = config::ThreadBudget {
            decode: 1,
            render: 1,
            encode: 1,
        };
        let records = (0..7).map(|_| Ok(a_read()));
        let got = settle(records, &mut cfg, budget, |r| {
            Cow::Borrowed(r.seq.as_slice())
        })
        .unwrap()
        .unwrap();
        assert_eq!(got.count(), 7);
    }

    fn base_config() -> Config {
        Config {
            io: config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: base_trim(),
            adapters: None,
            adapter_infer: config::AdapterInfer::Off,
            threads: 8,
            fastq_tags: config::FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            ordered: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
            progress: crate::config::ProgressMode::Auto,
            adapters_configured: None,
            trim_barcodes: false,
            remove_tags: crate::config::TagRemoval::default(),
        }
    }

    #[test]
    fn render_heavy_for_treats_bam_as_heavy() {
        let cfg = base_config();
        assert!(!config::render_heavy_for(io::Format::Fastq, &cfg));
        assert!(config::render_heavy_for(io::Format::Bam, &cfg));
    }
}
