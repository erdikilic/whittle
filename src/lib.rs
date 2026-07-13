pub mod adapter;
pub mod cli;
pub mod config;
pub mod filter;
pub mod io;
pub mod mods;
pub mod obs;
pub mod qual;
pub mod record;
pub mod trim;
pub mod workflow;

use std::borrow::Cow;
use std::io::{BufReader, BufWriter, IsTerminal, Read, Write};

use config::AdapterInfer;
pub use config::Config;
use gzp::deflate::Gzip;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use gzp::{Compression, ZWriter};

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
    // `&mut cfg` — the directory path itself is cloned out first.
    if let Some(dir) = cfg
        .io
        .input
        .as_deref()
        .filter(|p| p.is_dir())
        .map(|p| p.to_path_buf())
    {
        return run_folder(&dir, &mut cfg, obs);
    }

    // Refuse to read and write the same file: `whittle` streams the input, so
    // truncating it on `File::create` before it is fully read destroys the data
    // (a plain FASTQ run would silently emit an empty file with a success exit).
    if let (Some(inp), Some(outp)) = (cfg.io.input.as_deref(), cfg.io.output.as_deref())
        && same_path(inp, outp)
    {
        anyhow::bail!(
            "input and output are the same file ({}); whittle streams the input and \
             would truncate it before reading — write to a different path",
            outp.display()
        );
    }

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

    // Open the input (file or stdin) up front so format detection can sniff
    // its first bytes without losing them: any bytes consumed while sniffing
    // get prepended back via a Cursor+chain before the FASTQ reader is built.
    // Wrapped in `CountingReader` here, innermost, so it counts actual bytes
    // pulled from the file/stdin — the sniff bytes are counted once when
    // first read; re-serving them from the in-memory `Cursor` below does not
    // double-count them.
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
    // actual detection (this never changes behavior), but it usually signals a
    // mistake when it disagrees with the path's own extension — e.g.
    // `--in-format bam` on a `.fastq` file, or `--out-format fastq` on an
    // `out.fastq.gz` path (which would write a plain FASTQ into a .gz name).
    // Extension-only check: skipped for stdin/stdout / no-extension. Both
    // warnings fire later, after the banner (see the consolidated warnings
    // block below) — only the detection runs here.
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
    // conversion means the run just re-emits (almost) the same reads it read —
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

    // Resolved once, here, so the banner's Threads line and the actual dispatch
    // arm below agree on the same split — recomputing per arm risked the banner
    // showing one number and the workflow running another.
    let budget = config::thread_budget(
        cfg.threads,
        render_heavy_for(in_fmt, out_fmt, &cfg),
        matches!(in_fmt, Format::FastqBgzf),
        encode_kind_for(out_fmt),
    );
    configure_shared_bgzf_pool(
        matches!(in_fmt, Format::FastqBgzf | Format::Bam),
        matches!(out_fmt, Format::FastqBgzf | Format::Bam),
        budget,
    )?;
    let out_desc = output_desc(cfg.io.output.as_deref());

    if obs.shows_lines() {
        tracing::info!("{}", operation_line(in_fmt, out_fmt));
        match (in_path, total) {
            (Some(p), Some(size)) => {
                tracing::info!("Input: {} ({})", p.display(), obs::human_bytes(size));
            },
            (Some(p), None) => tracing::info!("Input: {}", p.display()),
            (None, _) => tracing::info!("Input: <stdin>"),
        }
        tracing::info!(
            "{}",
            output_banner_line(
                cfg.io.output.as_deref(),
                out_fmt,
                cfg.compression_level,
                budget.encode
            )
        );
        tracing::info!("{}", threads_banner_line(cfg.threads, budget));
        tracing::info!("{}", filters_and_trim_line(&cfg.filter, &cfg.trim));
        if let Some(line) =
            adapter_banner_line(cfg.adapters.as_ref(), cfg.adapter_sample, cfg.adapter_infer)
        {
            tracing::info!("{line}");
        }
    } else if obs.is_bar() {
        tracing::info!(
            "{} ({} threads)",
            operation_line(in_fmt, out_fmt),
            cfg.threads
        );
    }

    // Warnings fire after the resolved-config banner (not before it, and not
    // interleaved with it) — `whittle {version}`/`Command: ...` (printed by
    // `main` before `run` is even called) and the banner above are meant to be
    // the first things a reader sees; only then do clamp/mismatch/no-op
    // advisories follow, ahead of the live progress/summary.
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
            let Some(records) = maybe_reduce_adapters(records, &mut cfg, bam_seq)? else {
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
            // on `Drop`, whose `try_finish` error is silently discarded — an I/O
            // failure on final flush (e.g. ENOSPC) would otherwise yield a
            // truncated BAM with a success exit code.
            sink.finish()?;
            tracing::debug!("Processing finished in {}", obs::human_dur(t0.elapsed()));
            obs.finish(&stats, &out_desc);
            return Ok(());
        },
        (Format::Bam, Format::Fastq | Format::FastqGz | Format::FastqBgzf) => {
            // See the note in the (Bam, Bam) arm: read from the chained `source`.
            let (_header, records) = io::bam::reader_from(source, budget.decode)?;
            let Some(records) = maybe_reduce_adapters(records, &mut cfg, bam_seq)? else {
                return Ok(());
            };
            let mut writer = fastq_writer(&cfg, out_fmt, budget.encode)?;
            cfg.render_workers = budget.render;
            let stats = workflow::run_bam_to_fastq(records, &mut writer, &cfg, &counters)?;
            writer.finish()?;
            tracing::debug!("Processing finished in {}", obs::human_dur(t0.elapsed()));
            obs.finish(&stats, &out_desc);
            return Ok(());
        },
        (Format::Fastq | Format::FastqGz | Format::FastqBgzf, Format::Bam) => {
            anyhow::bail!("cross-format FASTQ->BAM conversion is not supported")
        },
        _ => {},
    }

    note_tags_ignored(&cfg, in_fmt, out_fmt);

    // Writer construction (a `File::create`, which eagerly truncates any
    // existing `-o` target) happens AFTER `maybe_reduce_adapters`, not before
    // — matching the BAM arms above. A `ReportOnly` early-exit (`Ok(None)`)
    // must return before any output file is touched; building the writer
    // first would truncate a pre-existing `-o` file even though report-only
    // writes no records at all.
    let records = match in_fmt {
        Format::Fastq => io::fastq::reader_from(source, false),
        Format::FastqGz => io::fastq::reader_from(source, true),
        Format::FastqBgzf => io::fastq::reader_from_bgzf(source, budget.decode)?,
        Format::Bam => unreachable!("BAM dispatch returned above"),
    };
    let Some(records) =
        maybe_reduce_adapters(records, &mut cfg, |r| Cow::Borrowed(r.seq.as_slice()))?
    else {
        return Ok(());
    };
    let mut writer = fastq_writer(&cfg, out_fmt, budget.encode)?;
    cfg.render_workers = budget.render;
    let stats = workflow::run_fastq(records, &mut writer, &cfg, &counters)?;
    writer.finish()?;
    tracing::debug!("Processing finished in {}", obs::human_dur(t0.elapsed()));
    obs.finish(&stats, &out_desc);
    Ok(())
}

/// Decode the packed SEQ of a lazy raw BAM record only when adapter sampling
/// needs it. Normal workflow records stay packed until a render worker converts
/// them to `RecordBuf`.
fn bam_seq(rec: &noodles_bam::Record) -> Cow<'_, [u8]> {
    Cow::Owned(rec.sequence().iter().collect())
}

/// A kept adapter's support below this is close enough to `infer::KEEP_SUPPORT`
/// (0.30, a whole-consensus presence fraction — see its doc comment) that
/// it's worth an explicit warning rather than trusting the plain info line:
/// a genuinely marginal discovery (e.g. an adapter only present in a fraction
/// of reads, such as a barcode-specific sequence) can still clear the keep
/// floor while remaining far from a confident, near-all-reads presence, so a
/// kept adapter this weak deserves scrutiny before trusting it, not silent
/// trust. ~1.5x `KEEP_SUPPORT` gives headroom above the floor while staying
/// well below a genuine, high-prevalence adapter's typical near-1.0 presence.
const MARGINAL_SUPPORT: f64 = 0.45;

/// Log each ab-initio discovery (`infer::discover` output) at
/// `info!`: `inferred_N ≈ NAME (pct%) · support X.XX` when the discovered
/// sequence cross-names against the built-in ONT catalog, else `inferred_N
/// (no catalog match) · support X.XX`. `N` is the 1-based position in
/// `discovered`'s own order (support desc, then sequence asc — see
/// `infer::discover`), which now agrees with any `inferred_{N}` fallback in
/// `InferredAdapter::adapter.name` (both derive from the same post-sort order).
/// The raw sequence is logged separately at `debug!` — too noisy for the
/// default INFO level, but useful with `-v` when checking a discovery by eye.
/// Support below `MARGINAL_SUPPORT` additionally gets a `warn!`, since it's
/// close enough to the `KEEP_SUPPORT` floor to warrant double-checking.
fn log_discovered(discovered: &[crate::adapter::infer::InferredAdapter], n_sampled: usize) {
    tracing::info!(
        "Adapter inference: sampled {n_sampled} reads, discovered {} adapter{}",
        discovered.len(),
        if discovered.len() == 1 { "" } else { "s" }
    );
    for (i, d) in discovered.iter().enumerate() {
        let n = i + 1;
        match d.name_hits.first() {
            Some((name, pct)) => {
                tracing::info!(
                    "inferred_{n} \u{2248} {name} ({pct:.0}%) \u{b7} support {:.2}",
                    d.support
                );
            },
            None => {
                tracing::info!(
                    "inferred_{n} (no catalog match) \u{b7} support {:.2}",
                    d.support
                );
            },
        }
        if d.support < MARGINAL_SUPPORT {
            tracing::warn!(
                "adapter '{}' support {:.2} is marginal (near the KEEP_SUPPORT floor); \
                 verify with --adapter-infer-only",
                d.adapter.name,
                d.support
            );
        }
        tracing::debug!(
            "inferred_{n} sequence: {}",
            String::from_utf8_lossy(&d.adapter.seq)
        );
    }
}

/// Print inferred adapters as FASTA with support and the best catalog match.
/// Numbering follows the final discovery order used by the status log.
fn print_discovered_fasta(discovered: &[crate::adapter::infer::InferredAdapter]) {
    for (i, d) in discovered.iter().enumerate() {
        let n = i + 1;
        let name_suffix = match d.name_hits.first() {
            Some((name, pct)) => format!(" [\u{2248} {name} ({pct:.0}%)]"),
            None => String::new(),
        };
        println!(">inferred_{n} support={:.2}{name_suffix}", d.support);
        println!("{}", String::from_utf8_lossy(&d.adapter.seq));
    }
}

/// Buffer at most `n` records, stopping when the input is exhausted.
fn buffer_prefix<R>(
    records: &mut impl Iterator<Item = anyhow::Result<R>>,
    n: usize,
) -> anyhow::Result<Vec<R>> {
    let mut sample = Vec::new();
    for _ in 0..n {
        match records.next() {
            Some(Ok(r)) => sample.push(r),
            Some(Err(e)) => return Err(e),
            None => break,
        }
    }
    Ok(sample)
}

/// Resolve adapter inference or presence sampling before workflow dispatch.
/// Report-only mode returns `None` so callers do not create an output writer;
/// trimming modes return the buffered prefix chained with the remaining input.
/// `seq_of` exposes a record's sequence without constraining its storage type.
fn maybe_reduce_adapters<R, I, F>(
    mut records: I,
    cfg: &mut Config,
    // The returned sequence view borrows the record passed to `seq_of`.
    seq_of: F,
) -> anyhow::Result<Option<Box<dyn Iterator<Item = anyhow::Result<R>> + Send>>>
where
    // Workflow iterators are boxed and may cross worker-thread boundaries.
    I: Iterator<Item = anyhow::Result<R>> + Send + 'static,
    R: Send + 'static,
    F: for<'a> Fn(&'a R) -> Cow<'a, [u8]>,
{
    if cfg.adapter_infer != AdapterInfer::Off {
        // Inference mode stores an empty configuration until discovery completes.
        let base = cfg
            .adapters
            .clone()
            .expect("adapter_infer != Off implies cfg.adapters is Some (see cli::parse)");

        let sample: Vec<R> = buffer_prefix(&mut records, cfg.adapter_sample)?;
        let s = sample.len();
        let chain =
            |sample: Vec<R>, records: I| -> Box<dyn Iterator<Item = anyhow::Result<R>> + Send> {
                Box::new(sample.into_iter().map(anyhow::Ok).chain(records))
            };
        if s < crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION {
            // Report-only mode must not create output when the sample is too small.
            tracing::warn!(
                "adapter inference: too few reads ({s}, need >= {}) to infer reliably; \
                 keeping reads untrimmed",
                crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION
            );
            if cfg.adapter_infer == AdapterInfer::ReportOnly {
                return Ok(None);
            }
            let mut reduced = base;
            reduced.replace_adapters(Vec::new());
            cfg.adapters = Some(reduced);
            return Ok(Some(chain(sample, records)));
        }

        let seq_storage: Vec<Cow<'_, [u8]>> = sample.iter().map(&seq_of).collect();
        let seqs: Vec<&[u8]> = seq_storage.iter().map(|s| s.as_ref()).collect();
        let discovered = crate::adapter::infer::discover(&seqs, &base);
        log_discovered(&discovered, s);

        if cfg.adapter_infer == AdapterInfer::ReportOnly {
            print_discovered_fasta(&discovered);
            return Ok(None);
        }

        if discovered.is_empty() {
            tracing::warn!(
                "adapter inference: no adapters inferred from the first {s} reads; keeping \
                 reads untrimmed"
            );
        }
        let mut reduced = base;
        reduced.replace_adapters(discovered.into_iter().map(|d| d.adapter).collect());
        cfg.adapters = Some(reduced);
        return Ok(Some(chain(sample, records)));
    }

    // Avoid buffering when neither inference nor presence sampling is active.
    if cfg.adapters.is_none() || cfg.adapter_sample == 0 {
        return Ok(Some(Box::new(records)));
    }

    // Reduce configured adapters using the sampled prefix.
    let mut sample: Vec<R> = Vec::new();
    if let Some(ac) = cfg.adapters.clone()
        && cfg.adapter_sample > 0
    {
        sample = buffer_prefix(&mut records, cfg.adapter_sample)?;
        let s = sample.len();
        let full = ac.adapters.len();
        let kept = if s < crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION {
            tracing::info!(
                "Adapter presence: only {s} reads (< {}); using all {full} adapters",
                crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION
            );
            ac.adapters.clone()
        } else {
            let seq_storage: Vec<Cow<'_, [u8]>> = sample.iter().map(&seq_of).collect();
            let seqs: Vec<&[u8]> = seq_storage.iter().map(|s| s.as_ref()).collect();
            let detected = crate::adapter::detect::present(
                &seqs,
                &ac.adapters,
                ac.error_rate,
                ac.end_size,
                ac.split,
                crate::adapter::detect::presence_min(s),
                cfg.threads,
            );
            if detected.is_empty() {
                tracing::warn!(
                    "Adapter presence: no adapters detected in the first {s} sampled reads; using all {full} \
                     (the sampled prefix may be unrepresentative — pass --adapter-sample 0 to always use the full set)"
                );
                ac.adapters.clone()
            } else {
                let names: Vec<&str> = detected.iter().take(12).map(|a| a.name.as_str()).collect();
                let more = detected.len().saturating_sub(names.len());
                tracing::info!(
                    "Adapter presence: sampled {s} reads, kept {} of {full} adapters{}{}",
                    detected.len(),
                    if names.is_empty() {
                        String::new()
                    } else {
                        format!(" ({})", names.join(", "))
                    },
                    if more > 0 {
                        format!(" +{more} more")
                    } else {
                        String::new()
                    },
                );
                detected
            }
        };
        let mut reduced = ac;
        reduced.replace_adapters(kept);
        cfg.adapters = Some(reduced);
    }
    Ok(Some(Box::new(
        sample.into_iter().map(anyhow::Ok).chain(records),
    )))
}

/// FASTQ output writer: either a plain buffered writer, or a `gzp` parallel
/// gzip writer (used only when the output format is explicitly `FastqGz` —
/// see `io::resolve_output`, which no longer auto-compresses).
///
/// `gzp`'s `ParCompress` REQUIRES an explicit `finish()` call to flush its
/// final (possibly partial) compressed block plus the gzip footer/checksum:
/// its `Write` impl only ever hands off *full* buffered chunks to the
/// compressor threads, so the tail end only ever gets flushed by `finish()`,
/// never by `flush()`. `ParCompress`'s own `Drop` impl does call `finish()` as
/// a backstop if it's still live, but it `.unwrap()`s the result — any I/O
/// error at that point becomes an uncatchable panic instead of a propagated
/// `anyhow::Result::Err`. Calling `finish()` explicitly, as the single seam
/// both callers below go through instead of `flush()` + `drop()`, keeps that
/// failure mode as an ordinary `Err`.
enum FastqOut {
    Plain(BufWriter<Box<dyn Write + Send>>),
    Gz(ParCompress<'static, Gzip, Box<dyn Write + Send>>),
    Bgzf(noodles_bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>),
}

impl Write for FastqOut {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            FastqOut::Plain(w) => w.write(buf),
            FastqOut::Gz(w) => w.write(buf),
            FastqOut::Bgzf(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            FastqOut::Plain(w) => w.flush(),
            FastqOut::Gz(w) => w.flush(),
            FastqOut::Bgzf(w) => w.flush(),
        }
    }
}

impl FastqOut {
    /// Finalize: for gz, flush the final block + gzip footer via gzp's
    /// `ZWriter::finish` (required — see the type's docs above); for plain,
    /// flush the `BufWriter`. Must be called before returning success.
    fn finish(self) -> anyhow::Result<()> {
        match self {
            FastqOut::Plain(mut w) => {
                w.flush()?;
                Ok(())
            },
            FastqOut::Gz(mut w) => {
                w.finish()?;
                Ok(())
            },
            FastqOut::Bgzf(mut w) => {
                w.finish()?;
                Ok(())
            },
        }
    }
}

/// Build the FASTQ output writer: a file or stdout, wrapped in a parallel gzip
/// encoder (`gzp`, using `gz_workers` worker threads — the caller's
/// workload-aware ENCODE share of the `-t` budget) when the output format is
/// `FastqGz`, else a plain buffered writer. `gz_workers` is only read on the
/// `FastqGz` branch, so plain-output callers may pass any value.
fn fastq_writer(cfg: &Config, out_fmt: io::Format, gz_workers: usize) -> anyhow::Result<FastqOut> {
    let base: Box<dyn Write + Send> = match cfg.io.output.as_deref() {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    };
    match out_fmt {
        io::Format::FastqGz => {
            let w = ParCompressBuilder::<Gzip>::new()
                .num_threads(gz_workers)
                .unwrap()
                .compression_level(Compression::new(cfg.compression_level as u32))
                .from_writer(base);
            Ok(FastqOut::Gz(w))
        },
        io::Format::FastqBgzf => {
            io::bam::configure_bgzf_pool(gz_workers)?;
            let level = noodles_bgzf::io::writer::CompressionLevel::new(cfg.compression_level)
                .ok_or_else(|| anyhow::anyhow!("invalid BGZF compression level"))?;
            let w = noodles_bgzf::io::multithreaded_writer::Builder::default()
                .set_compression_level(level)
                .build_from_writer(base);
            Ok(FastqOut::Bgzf(w))
        },
        io::Format::Fastq | io::Format::Bam => Ok(FastqOut::Plain(BufWriter::new(base))),
    }
}

/// True iff writing `fmt`'s bytes to stdout would dump binary (BAM) or gzip
/// (FASTQ.gz) data into an interactive terminal — never useful output, and
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
    if cfg.adapter_infer == AdapterInfer::ReportOnly {
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
            "refusing to write {} to a terminal — redirect to a file/pipe (e.g. `> out.{ext}`) \
             or pass -o",
            out_fmt.label()
        );
    }
    Ok(())
}

/// Folder-merge mode: `-i <dir>`. Classify the directory into one format family,
/// then merge all its read files into a single trimmed output using the same
/// workflows as the single-file path.
fn run_folder(
    dir: &std::path::Path,
    cfg: &mut Config,
    obs: &mut obs::ProgressHandle,
) -> anyhow::Result<()> {
    use io::Format;

    // Pass the output path so `classify` can hard-error if `-o` names a read file
    // inside `-i <dir>` — it could be a real input or a stale prior output, and
    // overwriting either while merging the rest is silent data loss. The merged
    // output must live outside the input directory.
    // --in-format has no effect here: a directory's family is decided per file
    // by extension (see io::dir::classify), so a forced input format is inert.
    // Warn rather than silently ignore it (the warning proper fires below, with
    // the other advisories, after the banner).
    let folder_in_format_ignored = cfg.io.in_format.is_some();

    let (family, paths) = io::dir::classify(dir, cfg.io.output.as_deref())?;
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
    let budget = config::thread_budget(
        cfg.threads,
        render_heavy_for(family_fmt, out_fmt, cfg),
        paths
            .iter()
            .any(|p| io::from_extension(p) == Some(Format::FastqBgzf)),
        encode_kind_for(out_fmt),
    );
    configure_shared_bgzf_pool(
        matches!(family_fmt, Format::Bam)
            || paths
                .iter()
                .any(|p| io::from_extension(p) == Some(Format::FastqBgzf)),
        matches!(out_fmt, Format::FastqBgzf | Format::Bam),
        budget,
    )?;
    let out_desc = output_desc(cfg.io.output.as_deref());

    if obs.shows_lines() {
        tracing::info!("{}", operation_line(family_fmt, out_fmt));
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
            output_banner_line(
                cfg.io.output.as_deref(),
                out_fmt,
                cfg.compression_level,
                budget.encode
            )
        );
        tracing::info!("{}", threads_banner_line(cfg.threads, budget));
        tracing::info!("{}", filters_and_trim_line(&cfg.filter, &cfg.trim));
        if let Some(line) =
            adapter_banner_line(cfg.adapters.as_ref(), cfg.adapter_sample, cfg.adapter_infer)
        {
            tracing::info!("{line}");
        }
    } else if obs.is_bar() {
        tracing::info!(
            "{} ({} threads)",
            operation_line(family_fmt, out_fmt),
            cfg.threads
        );
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
            let Some(records) =
                maybe_reduce_adapters(records, cfg, |r| Cow::Borrowed(r.seq.as_slice()))?
            else {
                return Ok(());
            };
            let mut writer = fastq_writer(cfg, out_fmt, budget.encode)?;
            cfg.render_workers = budget.render;
            let stats = workflow::run_fastq(records, &mut writer, cfg, &counters)?;
            writer.finish()?;
            tracing::debug!("Processing finished in {}", obs::human_dur(t0.elapsed()));
            obs.finish(&stats, &out_desc);
            Ok(())
        },
        io::dir::Family::Bam => match out_fmt {
            Format::Bam => {
                note_tags_ignored(cfg, family_fmt, out_fmt);
                // Only the first file's header is written; warn if the others
                // declare different read groups (relevant only for BAM output).
                io::dir::warn_on_bam_header_mismatch(&paths);
                let (header, records) = io::dir::bam_reader(&paths, budget.decode)?;
                let Some(records) = maybe_reduce_adapters(records, cfg, bam_seq)? else {
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
                tracing::debug!("Processing finished in {}", obs::human_dur(t0.elapsed()));
                obs.finish(&stats, &out_desc);
                Ok(())
            },
            Format::Fastq | Format::FastqGz | Format::FastqBgzf => {
                let (_header, records) = io::dir::bam_reader(&paths, budget.decode)?;
                let Some(records) = maybe_reduce_adapters(records, cfg, bam_seq)? else {
                    return Ok(());
                };
                let mut writer = fastq_writer(cfg, out_fmt, budget.encode)?;
                cfg.render_workers = budget.render;
                let stats = workflow::run_bam_to_fastq(records, &mut writer, cfg, &counters)?;
                writer.finish()?;
                tracing::debug!("Processing finished in {}", obs::human_dur(t0.elapsed()));
                obs.finish(&stats, &out_desc);
                Ok(())
            },
        },
    }
}

/// The output compression stage's weight for a given output format — `Bgzf` for
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

/// Noodles decode and encode jobs share Rayon's one global registry. When both
/// sides are BGZF, configure it before either reader or writer can lock in only
/// its own smaller stage allocation.
fn configure_shared_bgzf_pool(
    bgzf_input: bool,
    bgzf_output: bool,
    budget: config::ThreadBudget,
) -> anyhow::Result<()> {
    let workers = usize::from(bgzf_input && budget.decode > 1) * budget.decode
        + usize::from(bgzf_output && budget.encode > 1) * budget.encode;
    if workers > 1 {
        io::bam::configure_bgzf_pool(workers)?;
    }
    Ok(())
}

/// Whether the render stage has substantial per-record work. BAM input remains
/// render-heavy even for a full-window output because the current parallel path
/// still clones owned `RecordBuf`s before handing them to the writer. FASTQ
/// input is normally trim-only (light), but adapter matching or ab-initio
/// inference runs an approximate search per read, which is heavy too — so it
/// gets a render-pool share rather than being starved as pure compression.
fn render_heavy_for(in_fmt: io::Format, _out_fmt: io::Format, cfg: &Config) -> bool {
    matches!(in_fmt, io::Format::Bam)
        || cfg.adapters.is_some()
        || cfg.adapter_infer != AdapterInfer::Off
}

/// The startup banner's operation line (LINE mode's item 3 / BAR mode's own
/// one-liner build on the same wording): `Trimming FASTQ` when input and output
/// share a `Format::family` — including a `FASTQ` -> `FASTQ.gz` run, which is a
/// compression change, not a format conversion — else `Converting {in_label} to
/// {out_label}` (e.g. `Converting BAM to FASTQ`) for a genuine cross-family
/// conversion.
fn operation_line(in_fmt: io::Format, out_fmt: io::Format) -> String {
    if in_fmt.family() == out_fmt.family() {
        format!("Trimming {}", in_fmt.family())
    } else {
        format!("Converting {} to {}", in_fmt.label(), out_fmt.label())
    }
}

/// The startup banner's `Output: ...` line: `Output: <stdout>` when writing to
/// stdout (no compression detail), else `Output: {path}`, with
/// `(gzip|bgzf level {level}, {encode_workers} workers)` appended for compressed
/// output formats (gzip for `FASTQ.gz`, bgzf for BAM; plain FASTQ gets no suffix).
fn output_banner_line(
    output: Option<&std::path::Path>,
    out_fmt: io::Format,
    level: u8,
    encode_workers: usize,
) -> String {
    let Some(path) = output else {
        return "Output: <stdout>".to_string();
    };
    let mut line = format!("Output: {}", path.display());
    match out_fmt {
        io::Format::Bam => {
            line.push_str(&format!(" (bgzf level {level}, {encode_workers} workers)"));
        },
        io::Format::FastqGz => {
            line.push_str(&format!(" (gzip level {level}, {encode_workers} workers)"));
        },
        io::Format::FastqBgzf => {
            line.push_str(&format!(" (bgzf level {level}, {encode_workers} workers)"));
        },
        io::Format::Fastq => {},
    }
    line
}

/// The startup banner's `Threads: ...` line: the resolved `-t`/auto worker
/// count (`threads`) as the header, with the per-stage split (mapping the
/// `ThreadBudget`'s internal stage names — decode/render/encode — onto the
/// workflow-stage vocabulary shown to the user: read/trim/write) in
/// parentheses: `Threads: 8 (read 1, trim 4, write 3)`. Deliberately *not*
/// `b.total()`: that per-stage sum can exceed `threads` (each stage is floored
/// at >= 1 even when the overall total is 1 — see `ThreadBudget::total`'s
/// doc), which read as a confusing second, larger "total" next to the `-t`
/// value the user actually asked for.
///
/// `threads <= 1` instead prints `Threads: 1 (sequential)`: `thread_budget`
/// still floors `render`/`encode` at >= 1 each even at a total of 1, so the
/// read/trim/write split would show e.g. `(read 1, trim 1, write 1)` — three
/// threads' worth of detail for a run that is, in fact, single-threaded.
fn threads_banner_line(threads: usize, b: config::ThreadBudget) -> String {
    if threads <= 1 {
        return "Threads: 1 (sequential)".to_string();
    }
    format!(
        "Threads: {threads} (read {}, trim {}, write {})",
        b.decode, b.render, b.encode
    )
}

/// Lowercase label for a `QualMode`, used only in the startup banner's Filters
/// line (`{qual_mode} quality >=...`).
fn qual_mode_label(mode: qual::QualMode) -> &'static str {
    match mode {
        qual::QualMode::Mean => "mean",
        qual::QualMode::Arithmetic => "arithmetic",
        qual::QualMode::Median => "median",
    }
}

/// The startup banner's `Filters: ...; trim: ...` line, built from the resolved
/// `FilterConfig` + `TrimPlan`. Pure (no I/O), so it's unit-testable directly.
/// Shows only *active* (non-default) clauses/ops — a fresh-defaults run (no
/// filters, no trim) reads as `Filters: none; trim: none` rather than
/// spelling out every no-op threshold (e.g. `mean quality >=0`).
///
/// Filters clause: `length >={min}` only if `min_length > 1`, plus ` <={max}`
/// only if `max_length != usize::MAX`; `{qual_mode} quality >={min}` only if
/// `min_qual > 0.0`, plus ` <={max}` only if `max_qual < 1000.0`; `GC
/// {min}-{max}` only if either GC bound was set. `none` if nothing above fired.
///
/// Trim clause: `head {N}, tail {N}` only if either crop is non-zero, plus the
/// configured quality op's own wording, joined with a comma; `none` if neither
/// a crop nor a quality op is set.
fn filters_and_trim_line(filter: &filter::FilterConfig, trim: &trim::TrimPlan) -> String {
    let mut filters = Vec::new();

    let length_active = filter.min_length > 1 || filter.max_length != usize::MAX;
    if length_active {
        let mut length = String::new();
        if filter.min_length > 1 {
            length.push_str(&format!(">={}", filter.min_length));
        }
        if filter.max_length != usize::MAX {
            if !length.is_empty() {
                length.push(' ');
            }
            length.push_str(&format!("<={}", filter.max_length));
        }
        filters.push(format!("length {length}"));
    }

    let qual_active = filter.min_qual > 0.0 || filter.max_qual < 1000.0;
    if qual_active {
        let mut quality = format!("{} quality", qual_mode_label(filter.qual_mode));
        if filter.min_qual > 0.0 {
            quality.push_str(&format!(" >={}", filter.min_qual));
        }
        if filter.max_qual < 1000.0 {
            quality.push_str(&format!(" <={}", filter.max_qual));
        }
        filters.push(quality);
    }

    if filter.min_gc.is_some() || filter.max_gc.is_some() {
        filters.push(format!(
            "GC {}-{}",
            filter.min_gc.unwrap_or(0.0),
            filter.max_gc.unwrap_or(1.0)
        ));
    }

    let filters_str = if filters.is_empty() {
        "none".to_string()
    } else {
        filters.join("; ")
    };

    let mut trim_parts = Vec::new();
    if trim.head > 0 || trim.tail > 0 {
        trim_parts.push(format!("head {}, tail {}", trim.head, trim.tail));
    }
    if let Some(op) = &trim.quality {
        trim_parts.push(match op {
            trim::QualityOp::TrimQual(q) => format!("trim quality <{q}"),
            trim::QualityOp::BestSegment(q) => format!("best segment >={q}"),
            trim::QualityOp::Split { cutoff, .. } => format!("split quality <{cutoff}"),
        });
    }
    let trim_str = if trim_parts.is_empty() {
        "none".to_string()
    } else {
        trim_parts.join(", ")
    };

    format!("Filters: {filters_str}; trim: {trim_str}")
}

/// The startup banner's `Adapters: ...` line, shown only when adapter trimming
/// is active — `None` when `cfg.adapters` is unset, so callers can skip the
/// line entirely for an off run (same convention as the other banner-line
/// helpers being pure/unit-testable). Reports the adapter count, `trim +
/// split` vs `ends-only` mode (`AdapterConfig::split`), the end-match error
/// rate, the end-zone size in bp, and (via `adapter_sample`, i.e.
/// `cfg.adapter_sample`) whether presence detection will sample the input —
/// `sample {N}` when active, `sample off` when `N == 0` disables detection.
///
/// Under `AdapterInfer::Trim`/`ReportOnly`, the count printed here is always
/// `0` (discovery hasn't run yet — it replaces `cfg.adapters` only once the
/// buffer-then-decide seam runs, after this banner prints), so a `· infer` /
/// `· infer-only` suffix is appended to make clear the set is about to be
/// discovered, not that trimming is configured with zero adapters. This is
/// forced to `0` explicitly (rather than read off `a.adapters.len()`) because
/// under `ReportOnly` with a `--adapter-fasta`, `a.adapters` may itself hold
/// the user's FASTA entries (see `cli::parse`'s `trim_adapters`) -- carried
/// through purely as extra naming refs for `infer::discover`, never as a
/// trimming set, so they must not be counted here as if they were one.
fn adapter_banner_line(
    adapters: Option<&crate::adapter::AdapterConfig>,
    adapter_sample: usize,
    adapter_infer: AdapterInfer,
) -> Option<String> {
    let a = adapters?;
    let mode = if a.split { "trim + split" } else { "ends-only" };
    let sample = if adapter_sample > 0 {
        format!("sample {adapter_sample}")
    } else {
        "sample off".to_string()
    };
    let infer_suffix = match adapter_infer {
        AdapterInfer::Off => "",
        AdapterInfer::Trim => " \u{b7} infer",
        AdapterInfer::ReportOnly => " \u{b7} infer-only",
    };
    let n_adapters = if adapter_infer == AdapterInfer::Off {
        a.adapters.len()
    } else {
        0
    };
    Some(format!(
        "Adapters: {} sequences · {mode} · error {:.2} · end-zone {} bp · {sample}{infer_suffix}",
        n_adapters, a.error_rate, a.end_size
    ))
}

/// Shell-quote a single argument the way Python's `shlex.quote` does: bare when
/// non-empty and every character is in the POSIX-shell-safe set
/// (`[A-Za-z0-9_@%+=:,./-]`); otherwise wrapped in single quotes, with any
/// embedded single quote escaped as `'\''` (close the quote, an escaped literal
/// quote, reopen the quote). An empty argument is never safe bare (it would
/// vanish when re-run), so it renders as `''`.
pub(crate) fn shell_quote(arg: &str) -> String {
    let is_safe = |c: char| c.is_ascii_alphanumeric() || "_@%+=:,./-".contains(c);
    if !arg.is_empty() && arg.chars().all(is_safe) {
        arg.to_string()
    } else {
        format!("'{}'", arg.replace('\'', r"'\''"))
    }
}

/// The startup banner's `Command: ...` line: the real process argv, space-joined
/// and each argument shell-quoted via `shell_quote` so the line can be copied
/// back out and re-run verbatim. Takes `OsStr`-like items (the caller passes
/// `std::env::args_os()`, NOT `args()` — the latter panics on non-UTF-8 argv) and
/// lossily converts each to `str` here, at the one seam that must never panic on
/// a malformed argv. Generic over the argument iterator so it's unit-testable
/// without touching the real process argv.
pub fn command_line<I, S>(args: I) -> String
where
    I: IntoIterator<Item = S>,
    S: AsRef<std::ffi::OsStr>,
{
    let joined = args
        .into_iter()
        .map(|a| shell_quote(&a.as_ref().to_string_lossy()))
        .collect::<Vec<_>>()
        .join(" ");
    format!("Command: {joined}")
}

/// The output path (or `<stdout>`) shown in both the startup banner's `Output:`
/// line and the end-of-run `Completed`/closer line — the two bookend on the same
/// text so a reader can match them up at a glance.
fn output_desc(output: Option<&std::path::Path>) -> String {
    match output {
        Some(p) => p.display().to_string(),
        None => "<stdout>".to_string(),
    }
}

/// Whether two paths resolve to the same file. Canonicalizes both so symlinks
/// and `./`-style aliasing are caught; the output usually does not exist yet, so
/// it falls back to canonicalizing the parent directory and re-joining the file
/// name. Conservative: any resolution failure yields `false` (don't block a run
/// on a path we can't resolve).
pub(crate) fn same_path(a: &std::path::Path, b: &std::path::Path) -> bool {
    // When both paths already exist, an inode+device match is definitive — and it
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
/// duplicate ID), in which case the header is written unchanged — the `@PG` line
/// is cosmetic and must never block record output.
fn provenance_header(mut header: noodles_sam::Header) -> noodles_sam::Header {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::Program;
    use noodles_sam::header::record::value::map::program::tag;

    // `Programs::add` walks the existing `@PG` chain via `Programs::leaves`,
    // which indexes the program map directly and panics if any program's `PP`
    // (previous-program) field names an ID that isn't itself a program in the
    // header. Real-world uBAMs can have exactly this: e.g. an ONT/dorado file
    // put through `samtools sort`/`view`/`reset` observed with
    // `@PG ID:samtools PP:basecaller` where no `ID:basecaller` record survived
    // into the header. Since the `@PG` line is cosmetic, skip adding it rather
    // than let a merely-untidy header crash the whole run.
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

/// True if any `@PG` record's `PP` field references an ID that is not itself a
/// program in the header. `Programs::leaves` (used internally by
/// `Programs::add`) panics on such a chain instead of returning an error, so
/// this must be checked before calling `add`.
fn has_dangling_program_chain(header: &noodles_sam::Header) -> bool {
    use noodles_sam::header::record::value::map::program::tag;

    let programs = header.programs().as_ref();
    programs.values().any(|program| {
        program
            .other_fields()
            .get(&tag::PREVIOUS_PROGRAM_ID)
            .is_some_and(|previous_id| !programs.contains_key(previous_id))
    })
}

#[cfg(test)]
mod tests {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::Program;
    use noodles_sam::header::record::value::map::program::tag;

    use super::*;
    use crate::adapter::{Adapter, AdapterConfig, End};

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

    #[test]
    fn operation_line_collapses_matching_families() {
        assert_eq!(
            operation_line(io::Format::Fastq, io::Format::Fastq),
            "Trimming FASTQ"
        );
        // FASTQ -> FASTQ.gz shares the FASTQ family (a compression change, not
        // a format conversion), so it collapses too rather than reading as an
        // "X to X" conversion.
        assert_eq!(
            operation_line(io::Format::Fastq, io::Format::FastqGz),
            "Trimming FASTQ"
        );
    }

    #[test]
    fn operation_line_converting_wording_for_cross_family() {
        assert_eq!(
            operation_line(io::Format::Bam, io::Format::Fastq),
            "Converting BAM to FASTQ"
        );
    }

    #[test]
    fn output_banner_line_stdout_has_no_compression_detail() {
        // Even for a format that would otherwise show a compression suffix.
        assert_eq!(
            output_banner_line(None, io::Format::Bam, 6, 3),
            "Output: <stdout>"
        );
    }

    #[test]
    fn output_banner_line_plain_fastq_has_no_suffix() {
        let p = std::path::Path::new("/tmp/out.fastq");
        assert_eq!(
            output_banner_line(Some(p), io::Format::Fastq, 6, 3),
            "Output: /tmp/out.fastq"
        );
    }

    #[test]
    fn output_banner_line_appends_compression_detail() {
        let p = std::path::Path::new("/tmp/out.fastq.gz");
        assert_eq!(
            output_banner_line(Some(p), io::Format::FastqGz, 6, 4),
            "Output: /tmp/out.fastq.gz (gzip level 6, 4 workers)"
        );
        let p = std::path::Path::new("/tmp/out.bam");
        assert_eq!(
            output_banner_line(Some(p), io::Format::Bam, 3, 5),
            "Output: /tmp/out.bam (bgzf level 3, 5 workers)"
        );
    }

    #[test]
    fn output_desc_stdout_vs_path() {
        assert_eq!(output_desc(None), "<stdout>");
        assert_eq!(
            output_desc(Some(std::path::Path::new("/tmp/out.fastq"))),
            "/tmp/out.fastq"
        );
    }

    #[test]
    fn threads_banner_line_shows_requested_threads_not_the_stage_sum() {
        let b = config::thread_budget(8, true, false, config::EncodeKind::Bgzf);
        assert_eq!(
            threads_banner_line(8, b),
            format!(
                "Threads: 8 (read {}, trim {}, write {})",
                b.decode, b.render, b.encode
            )
        );
        // Concrete figure too, so a change in `thread_budget`'s split is noticed here.
        assert_eq!(
            threads_banner_line(8, b),
            "Threads: 8 (read 1, trim 4, write 3)"
        );
    }

    #[test]
    fn threads_banner_line_header_is_requested_even_when_stage_sum_differs() {
        // The banner reports the requested limit, not the sum of stage fields.
        let b = config::thread_budget(8, true, false, config::EncodeKind::None);
        assert_eq!(b.total(), 9);
        assert_eq!(
            threads_banner_line(8, b),
            "Threads: 8 (read 1, trim 7, write 1)"
        );
    }

    #[test]
    fn threads_banner_line_sequential_for_one_or_fewer() {
        // `-t 1` (or `-t 0`, which `resolve_threads` floors to 1): the
        // read/trim/write split would otherwise show e.g. "(read 1, trim 1,
        // write 1)" for what is actually a single-threaded run — collapse it
        // to a plain "sequential" label instead.
        let b = config::thread_budget(1, true, false, config::EncodeKind::Bgzf);
        assert_eq!(threads_banner_line(1, b), "Threads: 1 (sequential)");
    }

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

    #[test]
    fn filters_and_trim_line_defaults() {
        // All-default filter/trim: no active clause, so it reads "none" rather
        // than spelling out no-op thresholds like "mean quality >=0".
        assert_eq!(
            filters_and_trim_line(&base_filter(), &base_trim()),
            "Filters: none; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_min_length_active() {
        let mut f = base_filter();
        f.min_length = 500;
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: length >=500; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_max_length_active() {
        let mut f = base_filter();
        f.max_length = 10_000;
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: length <=10000; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_min_qual_active() {
        let mut f = base_filter();
        f.min_qual = 10.0;
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: mean quality >=10; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_max_qual_active() {
        let mut f = base_filter();
        f.max_qual = 40.0;
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: mean quality <=40; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_gc_active() {
        let mut f = base_filter();
        f.min_gc = Some(0.3);
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: GC 0.3-1; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_trim_active() {
        let mut t = base_trim();
        t.head = 5;
        assert_eq!(
            filters_and_trim_line(&base_filter(), &t),
            "Filters: none; trim: head 5, tail 0"
        );
    }

    #[test]
    fn filters_and_trim_line_all_bounds_set() {
        let mut f = base_filter();
        f.min_length = 200;
        f.max_length = 10_000;
        f.min_qual = 8.0;
        f.max_qual = 30.0;
        f.min_gc = Some(0.4);
        f.max_gc = Some(0.6);
        f.qual_mode = qual::QualMode::Median;

        let mut t = base_trim();
        t.head = 10;
        t.tail = 5;
        t.quality = Some(trim::QualityOp::TrimQual(12));

        assert_eq!(
            filters_and_trim_line(&f, &t),
            "Filters: length >=200 <=10000; median quality >=8 <=30; GC 0.4-0.6; \
             trim: head 10, tail 5, trim quality <12"
        );
    }

    #[test]
    fn filters_and_trim_line_quality_ops() {
        let f = base_filter();
        let mut t = base_trim();

        t.quality = Some(trim::QualityOp::BestSegment(20));
        assert!(filters_and_trim_line(&f, &t).ends_with("trim: best segment >=20"));

        t.quality = Some(trim::QualityOp::Split {
            cutoff: 15,
            window: 50,
        });
        assert!(filters_and_trim_line(&f, &t).ends_with("trim: split quality <15"));

        // head/tail-only (no quality op): no trailing quality-op clause.
        t.quality = None;
        t.head = 3;
        t.tail = 0;
        assert!(filters_and_trim_line(&f, &t).ends_with("trim: head 3, tail 0"));
    }

    #[test]
    fn adapter_banner_line_none_when_off_and_describes_when_on() {
        assert!(adapter_banner_line(None, 10000, AdapterInfer::Off).is_none());
        let cfg = AdapterConfig {
            adapters: vec![Adapter {
                name: "a".into(),
                seq: b"ACGTACGTACGT".to_vec(),
                end: End::Both,
            }],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let line = adapter_banner_line(Some(&cfg), 10000, AdapterInfer::Off).unwrap();
        assert!(line.contains("1 sequences"));
        assert!(line.contains("trim + split"));
        assert!(line.contains("error 0.20"));
        assert!(line.contains("end-zone 150 bp"));
        assert!(line.contains("sample 10000"));
        assert!(!line.contains("infer"), "no infer suffix when off: {line}");

        let off_line = adapter_banner_line(Some(&cfg), 0, AdapterInfer::Off).unwrap();
        assert!(off_line.contains("sample off"));
    }

    #[test]
    fn adapter_banner_line_ends_only_when_split_disabled() {
        let cfg = AdapterConfig {
            adapters: vec![Adapter {
                name: "a".into(),
                seq: b"ACGTACGTACGT".to_vec(),
                end: End::Both,
            }],
            error_rate: 0.2,
            end_size: 150,
            split: false,
            candidate_index: std::sync::OnceLock::new(),
        };
        assert!(
            adapter_banner_line(Some(&cfg), 10000, AdapterInfer::Off)
                .unwrap()
                .contains("ends-only")
        );
    }

    #[test]
    fn adapter_banner_line_notes_infer_mode() {
        let cfg = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let trim_line = adapter_banner_line(Some(&cfg), 40000, AdapterInfer::Trim).unwrap();
        assert!(trim_line.ends_with("infer"), "{trim_line}");
        assert!(!trim_line.ends_with("infer-only"), "{trim_line}");

        let report_line = adapter_banner_line(Some(&cfg), 40000, AdapterInfer::ReportOnly).unwrap();
        assert!(report_line.ends_with("infer-only"), "{report_line}");
    }

    #[test]
    fn command_line_quotes_only_unsafe_args() {
        assert_eq!(
            command_line(["whittle", "-i", "in.fastq", "-o", "out.fastq"]),
            "Command: whittle -i in.fastq -o out.fastq"
        );
        assert_eq!(
            command_line(["whittle", "-i", "my reads.fastq"]),
            "Command: whittle -i 'my reads.fastq'"
        );
    }

    #[test]
    fn shell_quote_leaves_plain_args_bare() {
        assert_eq!(shell_quote("whittle"), "whittle");
        assert_eq!(shell_quote("-i"), "-i");
        assert_eq!(shell_quote("in.fastq"), "in.fastq");
        assert_eq!(
            shell_quote("path/to/file_1.0.fq.gz"),
            "path/to/file_1.0.fq.gz"
        );
    }

    #[test]
    fn shell_quote_wraps_args_with_spaces() {
        assert_eq!(shell_quote("my reads.fastq"), "'my reads.fastq'");
    }

    #[test]
    fn shell_quote_escapes_embedded_single_quotes() {
        assert_eq!(shell_quote("it's here.fastq"), r"'it'\''s here.fastq'");
    }

    #[test]
    fn shell_quote_wraps_shell_metacharacters() {
        assert_eq!(shell_quote("$HOME"), "'$HOME'");
        assert_eq!(shell_quote("a;b"), "'a;b'");
    }

    #[test]
    fn shell_quote_wraps_empty_string() {
        // Bare would vanish entirely when the line is re-run.
        assert_eq!(shell_quote(""), "''");
    }
}
