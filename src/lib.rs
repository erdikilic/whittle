pub mod adapter;
pub mod banner;
pub mod cli;
pub mod config;
pub mod filter;
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
use std::io::{BufReader, IsTerminal, Read};

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

    guard_output_collisions(&cfg, &[])?;

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
    guard_stdout_binary(&cfg, out_fmt)?;

    // Advisory only: no trimming, a pass-through filter, and no format
    // conversion means the run just re-emits (almost) the same reads it read,
    // usually not what was intended. Skipped for a conversion-only run
    // (in_fmt != out_fmt), which is legitimate on its own. Warning deferred to
    // the consolidated block below, same as `mismatch_warn` above.
    let no_trim = cfg.trim.head == 0 && cfg.trim.tail == 0 && cfg.trim.quality.is_none();
    let pass_through_filter = cfg.filter.min_length <= 1
        && cfg.filter.max_length == usize::MAX
        && cfg.filter.min_qual <= 0.0
        && cfg.filter.max_qual >= 1000.0
        && cfg.filter.min_gc.is_none()
        && cfg.filter.max_gc.is_none();
    let no_op_warn = no_trim && pass_through_filter && cfg.adapters.is_none() && in_fmt == out_fmt;

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
        render_heavy_for(in_fmt, out_fmt, &cfg),
        parallel_decode,
        encode_kind_for(out_fmt),
    );
    let out_desc = banner::output_desc(cfg.io.output.as_deref());

    if obs.shows_lines() {
        tracing::info!("{}", banner::operation_line(in_fmt, out_fmt));
        match (in_path, total) {
            (Some(p), Some(size)) => {
                tracing::info!("Input: {} ({})", p.display(), obs::human_bytes(size));
            },
            (Some(p), None) => tracing::info!("Input: {}", p.display()),
            (None, _) => tracing::info!("Input: <stdin>"),
        }
        tracing::info!(
            "{}",
            banner::output_banner_line(
                cfg.io.output.as_deref(),
                out_fmt,
                cfg.compression_level,
                budget.encode
            )
        );
        tracing::info!("{}", banner::threads_banner_line(cfg.threads, budget));
        tracing::info!("{}", banner::filters_and_trim_line(&cfg.filter, &cfg.trim));
        if let Some(line) = banner::adapter_banner_line(
            cfg.adapters.as_ref(),
            cfg.adapter_sample,
            cfg.adapter_infer,
        ) {
            tracing::info!("{line}");
        }
    } else if obs.is_bar() {
        tracing::info!(
            "{} ({} threads)",
            banner::operation_line(in_fmt, out_fmt),
            cfg.threads
        );
    }

    // Warnings fire after the resolved-config banner (not before it, and not
    // interleaved with it): `whittle {version}`/`Command: ...` (printed by
    // `main` before `run` is even called) and the banner above are meant to be
    // the first things a reader sees; only then do clamp/mismatch/no-op
    // advisories follow, ahead of the live progress/summary.
    // Parse-time diagnostics, held by `cli::parse` until this subscriber exists
    // so `--quiet` can silence them and they carry the standard prefix.
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
    if let Some(msg) = mismatch_warn {
        tracing::warn!("{msg}");
    }
    if let Some(msg) = out_mismatch_warn {
        tracing::warn!("{msg}");
    }
    if no_op_warn {
        tracing::warn!("No trimming or filtering options set; output will mostly mirror the input");
    }

    obs.start(total, counters.clone());

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
            let out_header = provenance_header(header);
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

/// True iff writing `fmt`'s bytes to stdout would dump binary (BAM) or gzip
/// (FASTQ.gz) data into an interactive terminal: never useful output, and
/// almost always a forgotten `-o`/redirect. Plain FASTQ text is always fine.
/// Pure (no I/O) so it's trivial to unit-test without a real TTY.
fn binary_to_terminal(output_is_stdout: bool, fmt: io::Format, stdout_is_tty: bool) -> bool {
    output_is_stdout
        && stdout_is_tty
        && matches!(
            fmt,
            io::Format::Bam | io::Format::FastqGz | io::Format::FastqBgzf
        )
}

/// Reject binary output to an interactive terminal before creating a writer.
/// Report-only inference is exempt because it emits textual FASTA and exits
/// before workflow dispatch.
fn guard_stdout_binary(cfg: &Config, out_fmt: io::Format) -> anyhow::Result<()> {
    if cfg.adapter_infer.is_report() {
        return Ok(());
    }
    let stdout_is_tty = std::io::stdout().is_terminal();
    if binary_to_terminal(cfg.io.output.is_none(), out_fmt, stdout_is_tty) {
        let ext = match out_fmt {
            io::Format::Bam => "bam",
            io::Format::FastqGz => "fastq.gz",
            io::Format::FastqBgzf => "fastq.bgz",
            io::Format::Fastq => "fastq", // unreachable via binary_to_terminal, kept exhaustive
        };
        anyhow::bail!(
            "refusing to write {} to a terminal; redirect to a file/pipe (e.g. `> out.{ext}`) \
             or pass -o",
            out_fmt.label()
        );
    }
    Ok(())
}

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

/// Every file the run writes, checked against every file it reads.
///
/// `whittle` streams its input, so `File::create` truncating a file that is still
/// being read destroys it: a plain FASTQ run would emit an empty file and exit 0.
/// The same applies to `--summary-json`, which is written last and would replace
/// an input, the just-written output, or a folder member with JSON.
///
/// `extra_inputs` carries the folder-mode member files, which are not reachable
/// from `cfg` alone.
fn guard_output_collisions(
    cfg: &Config,
    extra_inputs: &[std::path::PathBuf],
) -> anyhow::Result<()> {
    let mut reads: Vec<(&str, &std::path::Path)> = Vec::new();
    if let Some(p) = cfg.io.input.as_deref() {
        reads.push(("the input", p));
    }
    if let Some(ac) = cfg.adapter_fasta.as_deref() {
        reads.push(("--adapter-fasta", ac));
    }
    for p in extra_inputs {
        reads.push(("an input file in the directory", p.as_path()));
    }

    for (what, dest) in [
        ("the output", cfg.io.output.as_deref()),
        ("--summary-json", cfg.summary_json.as_deref()),
    ] {
        let Some(dest) = dest else { continue };
        for (label, src) in &reads {
            if same_path(src, dest) {
                anyhow::bail!(
                    "{what} ({}) and {label} are the same file; whittle streams its input and \
                     would overwrite it, so write to a different path",
                    dest.display()
                );
            }
        }
        // With no `-i`, the input is stdin, which has no path to compare. Its
        // file descriptor still resolves to an inode, so a shell redirect from
        // the very file being written is caught here.
        if cfg.io.input.is_none() && stdin_is(dest) {
            anyhow::bail!(
                "{what} ({}) and the file being read on stdin are the same file; whittle \
                 streams its input and would truncate it mid-read, so write to a different path",
                dest.display()
            );
        }
    }

    // The summary is written after the output file, so it would clobber it.
    if let (Some(out), Some(sum)) = (cfg.io.output.as_deref(), cfg.summary_json.as_deref())
        && same_path(out, sum)
    {
        anyhow::bail!(
            "--summary-json ({}) and the output are the same file; the summary would replace \
             the trimmed reads with JSON",
            sum.display()
        );
    }
    Ok(())
}

/// Whether `path` names the same file the process has open on stdin.
///
/// A shell redirect (`whittle -o reads.fastq < reads.fastq`) leaves no path for
/// the same-file check to compare, but fd 0 still resolves to the inode.
#[cfg(unix)]
fn stdin_is(path: &std::path::Path) -> bool {
    use std::os::fd::AsFd;
    use std::os::unix::fs::MetadataExt;

    // `try_clone_to_owned` dups fd 0, so stdin itself is never closed here, and
    // it needs no `unsafe` (which this crate forbids).
    let Ok(dup) = std::io::stdin().as_fd().try_clone_to_owned() else {
        return false;
    };
    let Ok(m0) = std::fs::File::from(dup).metadata() else {
        return false;
    };
    let Ok(mp) = std::fs::metadata(path) else {
        return false;
    };
    // Only a regular file can collide; a pipe or tty shares no inode with a path.
    m0.is_file() && m0.dev() == mp.dev() && m0.ino() == mp.ino()
}

#[cfg(not(unix))]
fn stdin_is(_path: &std::path::Path) -> bool {
    false
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
    guard_output_collisions(cfg, &paths)?;
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
    guard_stdout_binary(cfg, out_fmt)?;

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
        render_heavy_for(family_fmt, out_fmt, cfg),
        parallel_decode,
        encode_kind_for(out_fmt),
    );
    let out_desc = banner::output_desc(cfg.io.output.as_deref());

    if obs.shows_lines() {
        tracing::info!("{}", banner::operation_line(family_fmt, out_fmt));
        let total_bytes: u64 = paths
            .iter()
            .filter_map(|p| std::fs::metadata(p).ok())
            .map(|m| m.len())
            .sum();
        tracing::info!(
            "Input: {} {} files, {}",
            paths.len(),
            family_fmt.label(),
            obs::human_bytes(total_bytes)
        );
        tracing::info!(
            "{}",
            banner::output_banner_line(
                cfg.io.output.as_deref(),
                out_fmt,
                cfg.compression_level,
                budget.encode
            )
        );
        tracing::info!("{}", banner::threads_banner_line(cfg.threads, budget));
        tracing::info!("{}", banner::filters_and_trim_line(&cfg.filter, &cfg.trim));
        if let Some(line) = banner::adapter_banner_line(
            cfg.adapters.as_ref(),
            cfg.adapter_sample,
            cfg.adapter_infer,
        ) {
            tracing::info!("{line}");
        }
    } else if obs.is_bar() {
        tracing::info!(
            "{} ({} threads)",
            banner::operation_line(family_fmt, out_fmt),
            cfg.threads
        );
    }

    // Parse-time diagnostics, held by `cli::parse` until this subscriber exists
    // so `--quiet` can silence them and they carry the standard prefix.
    for a in &cfg.advisories {
        if a.warn {
            tracing::warn!("{}", a.message);
        } else {
            tracing::info!("{}", a.message);
        }
    }

    // See the matching comment in `run`: the clamp warning fires after the
    // banner, not before it.
    if let Some((requested, ncpu)) = cfg.threads_clamped {
        tracing::warn!("Requested -t {requested} exceeds {ncpu} CPUs; using {ncpu}");
    }
    if folder_in_format_ignored {
        tracing::warn!(
            "--in-format is ignored for a directory input; folder files are classified \
             by extension per file"
        );
    }

    let counters = std::sync::Arc::new(workflow::Counters::default());
    obs.start(None, counters.clone());

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
                let out_header = provenance_header(header);
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

/// The output compression stage's weight for a given output format: `Bgzf` for
/// BAM (always bgzf-compressed), `Gzip` for `FASTQ.gz`, `None` for plain FASTQ.
/// Paired with `render_heavy` (`in_fmt == Format::Bam`, or the folder-mode
/// equivalent), this is everything `config::thread_budget` needs; both call sites
/// (`run`, `run_folder`) resolve their budget from this exactly once, before the
/// startup banner, and reuse it for the actual workflow dispatch below.
fn encode_kind_for(out_fmt: io::Format) -> config::EncodeKind {
    match out_fmt {
        io::Format::Bam => config::EncodeKind::Bgzf,
        io::Format::FastqGz => config::EncodeKind::Gzip,
        io::Format::FastqBgzf => config::EncodeKind::Bgzf,
        io::Format::Fastq => config::EncodeKind::None,
    }
}

/// Whether the render stage has substantial per-record work. BAM input remains
/// render-heavy even for a full-window output because the current parallel path
/// still clones owned `RecordBuf`s before handing them to the writer. FASTQ
/// input is normally trim-only (light), but adapter matching or ab-initio
/// inference runs an approximate search per read, which is heavy too, so it
/// gets a render-pool share rather than being starved as pure compression.
fn render_heavy_for(in_fmt: io::Format, _out_fmt: io::Format, cfg: &Config) -> bool {
    matches!(in_fmt, io::Format::Bam)
        || cfg.adapters.is_some()
        || cfg.adapter_infer != AdapterInfer::Off
}

/// Whether two paths resolve to the same file. Canonicalizes both so symlinks
/// and `./`-style aliasing are caught; the output usually does not exist yet, so
/// it falls back to canonicalizing the parent directory and re-joining the file
/// name. Conservative: any resolution failure yields `false` (don't block a run
/// on a path that cannot be resolved).
pub(crate) fn same_path(a: &std::path::Path, b: &std::path::Path) -> bool {
    // When both paths already exist, an inode+device match is definitive, and it
    // also catches hard links to one inode, which path canonicalization misses.
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(ma), Ok(mb)) = (std::fs::metadata(a), std::fs::metadata(b))
            && ma.dev() == mb.dev()
            && ma.ino() == mb.ino()
        {
            return true;
        }
    }
    fn resolve(p: &std::path::Path) -> Option<std::path::PathBuf> {
        if let Ok(c) = std::fs::canonicalize(p) {
            return Some(c);
        }
        let file = p.file_name()?;
        let parent = match p.parent() {
            Some(par) if !par.as_os_str().is_empty() => par,
            _ => std::path::Path::new("."),
        };
        std::fs::canonicalize(parent).ok().map(|c| c.join(file))
    }
    match (resolve(a), resolve(b)) {
        (Some(x), Some(y)) => x == y,
        _ => false,
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

/// Append an `@PG` provenance record (`ID:whittle`, program name + version) to a
/// cloned header before writing. Best-effort: `Programs::add` can fail (e.g. on a
/// duplicate ID), in which case the header is written unchanged. The `@PG` line
/// is cosmetic and must never block record output.
fn provenance_header(mut header: noodles_sam::Header) -> noodles_sam::Header {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::Program;
    use noodles_sam::header::record::value::map::program::tag;

    // `Programs::add` walks the `@PG` chain via `Programs::leaves`, which indexes
    // the program map directly and panics when a `PP` names an ID no longer in the
    // header. Real uBAMs hit this: a dorado file through `samtools reset` can keep
    // `@PG ID:samtools PP:basecaller` with no `ID:basecaller` record. The `@PG`
    // line is cosmetic, so skip it rather than crash on an untidy header.
    if has_dangling_program_chain(&header) {
        return header;
    }

    let program = Map::<Program>::builder()
        .insert(tag::NAME, "whittle")
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .build();

    if let Ok(program) = program {
        let _ = header.programs_mut().add("whittle", program);
    }

    header
}

/// True if the header's `@PG` chain is one `Programs::add` cannot walk safely.
///
/// `Programs::add` calls `Programs::leaves`, which indexes the program map
/// directly and panics when a `PP` names an absent ID, and which only terminates
/// a cycle that returns to the node it started from. A rho-shaped chain
/// (`pgA -> pgB -> pgC -> pgB`) has every ID present and never revisits `pgA`, so
/// it loops forever. Both shapes are rejected here by walking each chain with a
/// visited set.
fn has_dangling_program_chain(header: &noodles_sam::Header) -> bool {
    use std::collections::HashSet;

    use noodles_sam::header::record::value::map::program::tag;

    let programs = header.programs().as_ref();
    programs.keys().any(|start| {
        let mut seen: HashSet<&[u8]> = HashSet::new();
        let mut id: &[u8] = start.as_ref();
        loop {
            if !seen.insert(id) {
                return true; // revisited a node: cyclic
            }
            let Some(program) = programs.get(id) else {
                return true; // PP names an ID that is not a program: dangling
            };
            match program.other_fields().get(&tag::PREVIOUS_PROGRAM_ID) {
                Some(previous) => id = previous.as_ref(),
                None => return false, // reached the root of this chain
            }
        }
    })
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
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::Program;
    use noodles_sam::header::record::value::map::program::tag;

    use super::*;

    #[test]
    fn binary_to_terminal_flags_bam_on_a_tty_stdout() {
        assert!(binary_to_terminal(true, io::Format::Bam, true));
    }

    #[test]
    fn binary_to_terminal_flags_fastq_gz_on_a_tty_stdout() {
        assert!(binary_to_terminal(true, io::Format::FastqGz, true));
    }

    #[test]
    fn binary_to_terminal_allows_plain_fastq() {
        // Plain text FASTQ on a terminal is normal/expected output.
        assert!(!binary_to_terminal(true, io::Format::Fastq, true));
    }

    #[test]
    fn binary_to_terminal_allows_when_output_file_given() {
        // -o was given, so `output_is_stdout` is false regardless of format.
        assert!(!binary_to_terminal(false, io::Format::Bam, true));
    }

    #[test]
    fn binary_to_terminal_allows_when_not_a_tty() {
        // Redirected to a file/pipe: not a terminal, so it's fine.
        assert!(!binary_to_terminal(true, io::Format::Bam, false));
        assert!(!binary_to_terminal(true, io::Format::FastqGz, false));
    }

    /// A dangling `@PG PP:` reference must leave the header unchanged because
    /// Noodles requires every parent program ID to exist.
    #[test]
    fn provenance_header_does_not_panic_on_dangling_pp_chain() {
        // `pg1` references a parent that is absent from the header.
        let dangling_program = Map::<Program>::builder()
            .insert(tag::PREVIOUS_PROGRAM_ID, "ghost")
            .build()
            .expect("valid PP field");

        let header = noodles_sam::Header::builder()
            .add_program("pg1", dangling_program)
            .build();

        assert!(has_dangling_program_chain(&header));

        let out_header = provenance_header(header);

        assert!(
            !out_header.programs().as_ref().contains_key(&b"whittle"[..]),
            "expected no whittle @PG line to be added when the existing chain is dangling"
        );
    }

    /// A rho-shaped chain (`pgA -> pgB -> pgC -> pgB`) has no absent ID, so the
    /// old dangling-only check passed it through to `Programs::add`, whose
    /// `leaves()` walk only terminates on a cycle that returns to its start node.
    /// Walking from `pgA` never revisits `pgA`, so it looped forever at 100% CPU.
    #[test]
    fn provenance_header_rejects_a_cycle_that_excludes_the_entry_node() {
        fn with_pp(previous: &str) -> Map<Program> {
            Map::<Program>::builder()
                .insert(tag::PREVIOUS_PROGRAM_ID, previous)
                .build()
                .expect("valid PP field")
        }

        let header = noodles_sam::Header::builder()
            .add_program("pgA", with_pp("pgB"))
            .add_program("pgB", with_pp("pgC"))
            .add_program("pgC", with_pp("pgB"))
            .build();

        assert!(
            has_dangling_program_chain(&header),
            "a rho-shaped chain must be rejected before `Programs::add` sees it"
        );

        // Reaching this line at all is the assertion: the old code hung here.
        let out_header = provenance_header(header);
        assert!(
            !out_header.programs().as_ref().contains_key(&b"whittle"[..]),
            "no @PG line should be added when the existing chain cannot be walked"
        );
    }

    /// A self-referential record (`pgA -> pgA`) is the degenerate cycle.
    #[test]
    fn provenance_header_rejects_a_self_referential_program() {
        let header = noodles_sam::Header::builder()
            .add_program(
                "pgA",
                Map::<Program>::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pgA")
                    .build()
                    .expect("valid PP field"),
            )
            .build();
        assert!(has_dangling_program_chain(&header));
    }

    /// A valid program chain receives the `whittle` provenance record.
    #[test]
    fn provenance_header_adds_whittle_program_on_clean_header() {
        let header = noodles_sam::Header::default();
        assert!(!has_dangling_program_chain(&header));

        let out_header = provenance_header(header);

        assert!(
            out_header
                .programs()
                .roots()
                .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"whittle"),
            "expected an @PG record with ID whittle in the output header, got {:?}",
            out_header.programs()
        );
    }

    #[test]
    fn encode_kind_for_maps_output_format() {
        assert_eq!(encode_kind_for(io::Format::Bam), config::EncodeKind::Bgzf);
        assert_eq!(
            encode_kind_for(io::Format::FastqGz),
            config::EncodeKind::Gzip
        );
        assert_eq!(encode_kind_for(io::Format::Fastq), config::EncodeKind::None);
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
        assert!(!render_heavy_for(
            io::Format::Fastq,
            io::Format::FastqGz,
            &cfg
        ));
        assert!(render_heavy_for(io::Format::Bam, io::Format::Bam, &cfg));
        assert!(render_heavy_for(io::Format::Bam, io::Format::Fastq, &cfg));
    }
}
