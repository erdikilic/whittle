pub mod adapter;
pub mod banner;
pub mod cli;
pub mod config;
pub mod filter;
pub mod guards;
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

use config::AdapterInfer;
pub use config::Config;

/// Top-level entry point. Dispatches on the input: a directory triggers
/// folder-merge (all read files in it merged into one output); otherwise a
/// single file / stdin is trimmed. FASTQ and unaligned BAM are supported.
///
/// `obs` drives progress + end-of-run output; library callers pass `ProgressHandle::disabled()`.
pub fn run(cfg: Config, obs: &mut obs::ProgressHandle) -> anyhow::Result<()> {
    use io::Format;

    let mut cfg = cfg;
    let setup_start = std::time::Instant::now();

    // Scoped so the borrow of `cfg.io.input` ends before `run_folder` needs
    // `&mut cfg`. The directory path itself is cloned out first.
    if let Some(dir) = cfg
        .io
        .input
        .as_deref()
        .filter(|p| p.is_dir())
        .map(|p| p.to_path_buf())
    {
        return run_folder(&dir, &mut cfg, obs);
    }

    guards::guard_output_collisions(&cfg, &[])?;

    let in_path = cfg.io.input.as_deref();

    // Total input bytes, when known (a real file), drives a determinate
    // progress bar with %/ETA; stdin has no metadata, so it stays `None` and
    // renders a spinner instead (see `obs::ProgressHandle::start`).
    let total: Option<u64> = in_path
        .and_then(|p| std::fs::metadata(p).ok())
        .map(|m| m.len());

    // Created here (before the reader) so the same `Arc` can be shared into
    // `CountingReader` below, then cloned again for the workflow call and
    // `obs.start`.
    let counters = std::sync::Arc::new(workflow::Counters::default());

    // Open the input (file or stdin) up front so format detection can sniff its
    // first bytes without losing them: anything consumed gets prepended back via a
    // Cursor+chain before the FASTQ reader is built. `CountingReader` sits
    // innermost, so sniff bytes are counted once and re-serving them is free.
    let raw: Box<dyn Read + Send> = match in_path {
        Some(p) => Box::new(std::fs::File::open(p)?),
        None => Box::new(std::io::stdin()),
    };
    let raw: Box<dyn Read + Send> =
        Box::new(io::counting::CountingReader::new(raw, counters.clone()));
    let mut source: Box<dyn Read + Send> = Box::new(BufReader::new(raw));

    let in_fmt = match cfg.io.in_format {
        Some(f) => f,
        None => match in_path.and_then(io::from_extension) {
            Some(f) => f,
            None => {
                // Probe enough bytes to see a full BGZF block header (18 bytes),
                // so a BAM read from stdin/an unknown extension is told apart from
                // gzipped FASTQ. A single `read()` may return fewer; loop to fill.
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
                    source.read_exact(&mut replay[n..])?;
                    io::detect_bgzf_block(&replay)?
                } else {
                    io::detect_input(in_path, &replay)?
                };
                source = Box::new(std::io::Cursor::new(replay).chain(source));
                fmt
            },
        },
    };

    // Advisory only: an explicit --in-format/--out-format always wins for
    // detection, but usually signals a mistake when it disagrees with the path's
    // extension, e.g. `--out-format fastq` on an `out.fastq.gz` path. Skipped for
    // stdin/stdout and paths without an extension. Only the detection runs here;
    // both warnings fire after the banner, with the other advisories.
    let mismatch_warn = io::format_mismatch_warning("--in-format", cfg.io.in_format, in_path);
    let out_mismatch_warn =
        io::format_mismatch_warning("--out-format", cfg.io.out_format, cfg.io.output.as_deref());

    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), in_fmt));

    // Hard-error before any writer/output file is created: dumping BAM or
    // gzipped bytes into an interactive terminal is never useful and almost
    // always means the user forgot `-o`/a redirect.
    guards::guard_stdout_binary(&cfg, out_fmt)?;

    tracing::debug!(
        "Detected {} input in {}",
        in_fmt.label(),
        obs::human_dur(setup_start.elapsed())
    );

    // Resolved once here so the banner's Threads line and the dispatch arm below
    // agree; recomputing per arm risked showing one number and running another.
    // BAM and bgzf-FASTQ inputs are BGZF containers whose decode can go
    // multithreaded, but grant that budget only when render is light: adapter
    // search (preset, FASTA, or inference) makes render the bottleneck instead.
    let parallel_decode = matches!(in_fmt, Format::FastqBgzf | Format::Bam)
        && cfg.adapters.is_none()
        && cfg.adapter_infer == AdapterInfer::Off;
    let budget = config::thread_budget(
        cfg.threads,
        config::render_heavy_for(in_fmt, out_fmt, &cfg),
        parallel_decode,
        config::encode_kind_for(out_fmt),
    );
    let out_desc = banner::output_desc(cfg.io.output.as_deref());

    let mut warnings: Vec<String> = Vec::new();
    warnings.extend(mismatch_warn);
    warnings.extend(out_mismatch_warn);
    if is_no_op(&cfg, in_fmt, out_fmt) {
        warnings.push(NO_OP_WARNING.to_string());
    }
    announce(
        &cfg,
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

    // Coarse wall-clock timer for the processing phase (dispatch below); each
    // arm logs elapsed time from this point just before its own `obs.finish`.
    // Stages run concurrently internally (read/trim/write overlap across
    // threads), so this is a phase boundary, not a CPU-time split.
    let t0 = std::time::Instant::now();
    tracing::debug!("Processing {}, {} threads", in_fmt.label(), cfg.threads);

    // BAM dispatch happens before creating/truncating the output file, and so
    // do the FASTQ->BAM rejection and the BAM->FASTQ conversion, so a rejected
    // run never leaves a stray 0-byte file behind. Only the (Fastq*, Fastq*)
    // combinations fall through to the FASTQ path below.
    match (in_fmt, out_fmt) {
        (Format::Bam, Format::Bam) => {
            note_tags_ignored(&cfg, in_fmt, out_fmt);
            // Read from `source` (not by re-opening `in_path`): for a stdin BAM the
            // sniff bytes were already consumed and chained back into `source`, so
            // re-opening stdin would drop the BGZF header. For a file, `source` is
            // the same handle positioned at the start.
            let (header, records) = io::bam::reader_from(source, budget.decode)?;
            let Some(records) = adapter::resolve::maybe_reduce_adapters(
                records,
                &mut cfg,
                adapter::resolve::bam_seq,
            )?
            else {
                return Ok(());
            };
            // Append the invocation's @PG provenance line before writing.
            let out_header = io::bam::provenance_header(header);
            let mut sink = io::bam::writer(
                cfg.io.output.as_deref(),
                &out_header,
                budget.encode,
                cfg.compression_level,
            )?;
            cfg.render_workers = budget.render;
            let stats = workflow::run_raw_bam(&out_header, records, &mut sink, &cfg, &counters)?;
            // Explicitly finish (final bgzf block + EOF marker) instead of relying
            // on `Drop`, whose `try_finish` error is silently discarded. An I/O
            // failure on final flush (e.g. ENOSPC) would otherwise yield a
            // truncated BAM with a success exit code.
            sink.finish()?;
            finish_run(obs, &stats, &out_desc, &cfg, t0)?;
            return Ok(());
        },
        (Format::Bam, Format::Fastq | Format::FastqGz | Format::FastqBgzf) => {
            // See the note in the (Bam, Bam) arm: read from the chained `source`.
            let (_header, records) = io::bam::reader_from(source, budget.decode)?;
            let Some(records) = adapter::resolve::maybe_reduce_adapters(
                records,
                &mut cfg,
                adapter::resolve::bam_seq,
            )?
            else {
                return Ok(());
            };
            let mut writer = io::fastq::writer(&cfg, out_fmt, budget.encode)?;
            cfg.render_workers = budget.render;
            let stats = workflow::run_bam_to_fastq(records, &mut writer, &cfg, &counters)?;
            writer.finish()?;
            finish_run(obs, &stats, &out_desc, &cfg, t0)?;
            return Ok(());
        },
        (Format::Fastq | Format::FastqGz | Format::FastqBgzf, Format::Bam) => {
            anyhow::bail!("cross-format FASTQ->BAM conversion is not supported")
        },
        _ => {},
    }

    note_tags_ignored(&cfg, in_fmt, out_fmt);

    // Writer construction (a `File::create`, which eagerly truncates any
    // existing `-o` target) happens AFTER `adapter::resolve::maybe_reduce_adapters`, not before,
    // matching the BAM arms above. An inference-report early exit (`Ok(None)`)
    // must return before any output file is touched; building the writer
    // first would truncate a pre-existing `-o` file even though report-only
    // writes no records at all.
    let records = match in_fmt {
        Format::Fastq => io::fastq::reader_from(source, false),
        Format::FastqGz => io::fastq::reader_from(source, true),
        Format::FastqBgzf => io::fastq::reader_from_bgzf(source, budget.decode)?,
        Format::Bam => unreachable!("BAM dispatch returned above"),
    };
    let Some(records) = adapter::resolve::maybe_reduce_adapters(records, &mut cfg, |r| {
        Cow::Borrowed(r.seq.as_slice())
    })?
    else {
        return Ok(());
    };
    let mut writer = io::fastq::writer(&cfg, out_fmt, budget.encode)?;
    cfg.render_workers = budget.render;
    let stats = workflow::run_fastq(records, &mut writer, &cfg, &counters)?;
    writer.finish()?;
    finish_run(obs, &stats, &out_desc, &cfg, t0)?;
    Ok(())
}

/// The resolved facts both entry points announce before processing starts.
struct Startup {
    /// The `Input:` line body: a path and size, `<stdin>`, or a folder's file
    /// count and total size.
    input_line: String,
    /// Input bytes, when known. Drives a determinate bar; `None` gives a spinner.
    total: Option<u64>,
    in_fmt: io::Format,
    out_fmt: io::Format,
    budget: config::ThreadBudget,
    /// Advisories specific to this entry point, logged after the shared ones.
    warnings: Vec<String>,
}

/// Print the startup banner and the deferred advisories.
///
/// The order is the output contract: `whittle {version}` and `Command:` (from
/// `main`), then the resolved config, then every advisory, then live progress. A
/// reader can always find what ran at the top, and nothing interleaves with the
/// bar. Both entry points go through here, so an advisory added to one cannot
/// silently skip the other, which is how the folder path came to be missing two.
///
/// Starting progress is the last step, and belongs here rather than at the call
/// sites: it is what makes "nothing interleaves with the bar" a property of this
/// function instead of a convention two callers have to remember.
fn announce(
    cfg: &Config,
    obs: &mut obs::ProgressHandle,
    counters: std::sync::Arc<workflow::Counters>,
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
        // Bar mode gets exactly one line so the live bar stays clean.
        tracing::info!(
            "{} ({} threads)",
            banner::operation_line(s.in_fmt, s.out_fmt),
            cfg.threads
        );
    }

    // Parse-time diagnostics, held by `cli::parse` until the subscriber exists so
    // `--quiet` can silence them and they carry the standard prefix.
    for a in &cfg.advisories {
        if a.warn {
            tracing::warn!("{}", a.message);
        } else {
            tracing::info!("{}", a.message);
        }
    }
    if let Some((requested, ncpu)) = cfg.threads_clamped {
        tracing::warn!("Requested -t {requested} exceeds {ncpu} CPUs; using {ncpu}");
    }
    for w in &s.warnings {
        tracing::warn!("{w}");
    }

    obs.start(s.total, counters);
}

/// True when the run neither trims, filters, nor converts, so it just re-emits
/// (almost) what it read. Legitimate for a conversion-only run, which is why the
/// format comparison is part of the test.
fn is_no_op(cfg: &Config, in_fmt: io::Format, out_fmt: io::Format) -> bool {
    let no_trim = cfg.trim.head == 0 && cfg.trim.tail == 0 && cfg.trim.quality.is_none();
    let pass_through_filter = cfg.filter.min_length <= 1
        && cfg.filter.max_length == usize::MAX
        && cfg.filter.min_qual <= 0.0
        && cfg.filter.max_qual >= 1000.0
        && cfg.filter.min_gc.is_none()
        && cfg.filter.max_gc.is_none();
    no_trim && pass_through_filter && cfg.adapters.is_none() && in_fmt == out_fmt
}

const NO_OP_WARNING: &str =
    "No trimming or filtering options set; output will mostly mirror the input";

/// Close out a finished dispatch: log the phase duration, print the end-of-run
/// summary, and write `--summary-json` when one was requested. The one seam every
/// dispatch arm goes through, so a new arm cannot forget the summary file.
fn finish_run(
    obs: &mut obs::ProgressHandle,
    stats: &workflow::Stats,
    out_desc: &str,
    cfg: &Config,
    t0: std::time::Instant,
) -> anyhow::Result<()> {
    tracing::debug!("Processing finished in {}", obs::human_dur(t0.elapsed()));
    let elapsed = obs.finish(stats, out_desc);
    if let Some(path) = cfg.summary_json.as_deref() {
        summary::Summary::new(
            cfg,
            stats,
            command_line(std::env::args_os()),
            out_desc.to_string(),
            elapsed,
        )
        .write(path)?;
    }
    Ok(())
}

/// Folder-merge mode: `-i <dir>`. Classify the directory into one format family,
/// then merge its read files into a single trimmed output through the same
/// workflows as the single-file path.
fn run_folder(
    dir: &std::path::Path,
    cfg: &mut Config,
    obs: &mut obs::ProgressHandle,
) -> anyhow::Result<()> {
    use io::Format;

    // Pass the output path so `classify` can hard-error when `-o` names a read
    // file inside `-i <dir>`: real input or stale prior output, overwriting either
    // while merging the rest is silent data loss. `--in-format` is inert here,
    // since a directory's family is decided per file by extension, so it gets a
    // warning below rather than being silently ignored.
    let folder_in_format_ignored = cfg.io.in_format.is_some();

    let (family, paths) = io::dir::classify(dir, cfg.io.output.as_deref())?;
    // `classify` already refuses an output inside the directory; this also covers
    // --summary-json and --adapter-fasta against every member file.
    guards::guard_output_collisions(cfg, &paths)?;
    let family_fmt = match family {
        io::dir::Family::Fastq => Format::Fastq,
        io::dir::Family::Bam => Format::Bam,
    };
    let out_fmt = cfg
        .io
        .out_format
        .unwrap_or_else(|| io::resolve_output(cfg.io.output.as_deref(), family_fmt));

    // Hard-error before any writer/output file is created (see `run`'s
    // matching guard for the single-file path).
    guards::guard_stdout_binary(cfg, out_fmt)?;

    // Resolved once, here, so the banner's Threads line and the actual dispatch
    // arm below agree on the same split (see the matching comment in `run`).
    let bgzf_input = family_fmt == Format::Bam
        || paths
            .iter()
            .any(|p| io::from_extension(p) == Some(Format::FastqBgzf));
    // See `run`: parallel decode only when render is light (no adapter search).
    let parallel_decode =
        bgzf_input && cfg.adapters.is_none() && cfg.adapter_infer == AdapterInfer::Off;
    let budget = config::thread_budget(
        cfg.threads,
        config::render_heavy_for(family_fmt, out_fmt, cfg),
        parallel_decode,
        config::encode_kind_for(out_fmt),
    );
    let out_desc = banner::output_desc(cfg.io.output.as_deref());
    let counters = std::sync::Arc::new(workflow::Counters::default());

    // Summed unconditionally, not just when the banner prints: it also drives a
    // determinate progress bar, which folder mode previously never got.
    let total_bytes: u64 = paths
        .iter()
        .filter_map(|p| std::fs::metadata(p).ok())
        .map(|m| m.len())
        .sum();

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
    if is_no_op(cfg, family_fmt, out_fmt) {
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

    let t0 = std::time::Instant::now();
    tracing::debug!(
        "Processing folder ({}), {} threads",
        family_fmt.label(),
        cfg.threads
    );

    match family {
        io::dir::Family::Fastq => {
            if matches!(out_fmt, Format::Bam) {
                anyhow::bail!(
                    "cross-format conversion (FASTQ folder to BAM) is not supported in v1"
                );
            }
            note_tags_ignored(cfg, family_fmt, out_fmt);
            // Resolve report-only mode before creating the output file.
            let records = io::dir::fastq_records(&paths, budget.decode);
            let Some(records) = adapter::resolve::maybe_reduce_adapters(records, cfg, |r| {
                Cow::Borrowed(r.seq.as_slice())
            })?
            else {
                return Ok(());
            };
            let mut writer = io::fastq::writer(cfg, out_fmt, budget.encode)?;
            cfg.render_workers = budget.render;
            let stats = workflow::run_fastq(records, &mut writer, cfg, &counters)?;
            writer.finish()?;
            finish_run(obs, &stats, &out_desc, cfg, t0)?;
            Ok(())
        },
        io::dir::Family::Bam => match out_fmt {
            Format::Bam => {
                note_tags_ignored(cfg, family_fmt, out_fmt);
                // Only the first file's header is written; warn if the others
                // declare different read groups (relevant only for BAM output).
                io::dir::warn_on_bam_header_mismatch(&paths);
                let (header, records) = io::dir::bam_reader(&paths, budget.decode)?;
                let Some(records) = adapter::resolve::maybe_reduce_adapters(
                    records,
                    cfg,
                    adapter::resolve::bam_seq,
                )?
                else {
                    return Ok(());
                };
                let out_header = io::bam::provenance_header(header);
                let mut sink = io::bam::writer(
                    cfg.io.output.as_deref(),
                    &out_header,
                    budget.encode,
                    cfg.compression_level,
                )?;
                cfg.render_workers = budget.render;
                let stats = workflow::run_raw_bam(&out_header, records, &mut sink, cfg, &counters)?;
                sink.finish()?;
                finish_run(obs, &stats, &out_desc, cfg, t0)?;
                Ok(())
            },
            Format::Fastq | Format::FastqGz | Format::FastqBgzf => {
                let (_header, records) = io::dir::bam_reader(&paths, budget.decode)?;
                let Some(records) = adapter::resolve::maybe_reduce_adapters(
                    records,
                    cfg,
                    adapter::resolve::bam_seq,
                )?
                else {
                    return Ok(());
                };
                let mut writer = io::fastq::writer(cfg, out_fmt, budget.encode)?;
                cfg.render_workers = budget.render;
                let stats = workflow::run_bam_to_fastq(records, &mut writer, cfg, &counters)?;
                writer.finish()?;
                finish_run(obs, &stats, &out_desc, cfg, t0)?;
                Ok(())
            },
        },
    }
}

/// `--fastq-tags` only affects BAM→FASTQ output. When the user set a non-default
/// value (`none`/an explicit list) on any other path, emit a one-line stderr note
/// rather than silently ignoring it. (An explicit `all` is the default and stays
/// silent.)
fn note_tags_ignored(cfg: &Config, in_fmt: io::Format, out_fmt: io::Format) {
    if !matches!(cfg.fastq_tags, config::FastqTags::All) {
        tracing::warn!(
            "--fastq-tags applies only to BAM-to-FASTQ output; ignored for {} to {}",
            in_fmt.label(),
            out_fmt.label()
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
        }
    }

    #[test]
    fn render_heavy_for_treats_bam_as_heavy() {
        let cfg = base_config();
        assert!(!config::render_heavy_for(
            io::Format::Fastq,
            io::Format::FastqGz,
            &cfg
        ));
        assert!(config::render_heavy_for(
            io::Format::Bam,
            io::Format::Bam,
            &cfg
        ));
        assert!(config::render_heavy_for(
            io::Format::Bam,
            io::Format::Fastq,
            &cfg
        ));
    }
}
