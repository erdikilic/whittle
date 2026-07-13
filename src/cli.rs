use std::path::PathBuf;

use clap::Parser;

use crate::config::{AdapterInfer, Config, FastqTags, IoConfig};
use crate::filter::FilterConfig;
use crate::io::Format;
use crate::qual::QualMode;
use crate::trim::{QualityOp, TrimPlan};

#[derive(Parser, Debug)]
#[command(author, version, about = "Tag-aware long-read trimmer", long_about = None)]
struct Cli {
    #[arg(short = 'i', long, help_heading = "Setup")]
    input: Option<PathBuf>,
    #[arg(short = 'o', long, help_heading = "Setup")]
    output: Option<PathBuf>,
    #[arg(long, value_enum, help_heading = "Setup")]
    in_format: Option<FormatArg>,
    #[arg(long, value_enum, help_heading = "Setup")]
    out_format: Option<FormatArg>,
    /// Worker threads (default: all detected CPUs; values above the CPU count are clamped).
    #[arg(short = 't', long, help_heading = "Setup")]
    threads: Option<usize>,
    #[arg(long, default_value = "all", help_heading = "Setup")]
    fastq_tags: String,
    /// DEFLATE compression level for compressed output (bgzf for BAM, gzip for
    /// FASTQ.gz). Lower = faster/larger. Ignored for plain FASTQ.
    #[arg(short = 'c', long, default_value_t = 6, help_heading = "Setup")]
    compression_level: u8,

    /// Increase logging detail: -v = debug, -vv = trace. Overridden by WHITTLE_LOG.
    #[arg(short = 'v', long, action = clap::ArgAction::Count, help_heading = "Logging")]
    verbose: u8,
    /// Silence progress and info output; warnings and errors still print.
    #[arg(long, conflicts_with = "verbose", help_heading = "Logging")]
    quiet: bool,

    #[arg(short = 'l', long, default_value_t = 1, help_heading = "Filtering")]
    min_length: usize,
    #[arg(short = 'L', long, help_heading = "Filtering")]
    max_length: Option<usize>,
    #[arg(short = 'q', long, default_value_t = 0.0, help_heading = "Filtering")]
    min_qual: f64,
    #[arg(
        short = 'Q',
        long,
        default_value_t = 1000.0,
        help_heading = "Filtering"
    )]
    max_qual: f64,
    #[arg(short = 'g', long, help_heading = "Filtering")]
    min_gc: Option<f64>,
    #[arg(short = 'G', long, help_heading = "Filtering")]
    max_gc: Option<f64>,
    #[arg(short = 'm', long, value_enum, default_value_t = QualModeArg::Mean, help_heading = "Filtering")]
    qual_mode: QualModeArg,

    #[arg(short = 'H', long, default_value_t = 0, help_heading = "Trimming")]
    head_crop: usize,
    #[arg(short = 'T', long, default_value_t = 0, help_heading = "Trimming")]
    tail_crop: usize,
    #[arg(long, help_heading = "Trimming")]
    qual_trim: Option<u8>,
    #[arg(long, help_heading = "Trimming")]
    qual_best_segment: Option<u8>,
    #[arg(long, help_heading = "Trimming")]
    qual_split: Option<u8>,
    #[arg(long, default_value_t = 1, help_heading = "Trimming")]
    qual_split_window: usize,
    /// Keep ONT signal tags consistent through trimming (slice `mv`, update
    /// `ts`/`ns`/`sp`/`pi`) for signal-aware tools (Remora, Clair3 v2), instead
    /// of dropping them. BAM→BAM only.
    #[arg(long, help_heading = "Trimming")]
    update_moves: bool,

    /// Adapter FASTA (each sequence >= the 11-bp minimum match length).
    /// Enables adapter trimming.
    #[arg(short = 'a', long, help_heading = "Adapter trimming")]
    adapter_fasta: Option<PathBuf>,
    /// Built-in ONT adapter catalog. Enables adapter trimming.
    #[arg(long, value_enum, default_value_t = AdapterPresetArg::None, help_heading = "Adapter trimming")]
    adapter_preset: AdapterPresetArg,
    /// End-match tolerance (fraction of adapter length). Interior splits use half.
    #[arg(long, default_value_t = 0.2, help_heading = "Adapter trimming")]
    adapter_error_rate: f64,
    /// Bases at each read end searched for a terminal adapter.
    #[arg(long, default_value_t = 150, help_heading = "Adapter trimming")]
    adapter_end_size: usize,
    /// Trim adapters at read ends only; never split on interior adapters.
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_ends_only: bool,
    /// Reads to sample. Detection (>=100) or, under --adapter-infer, the inference
    /// buffer (default 40000). Omitted = mode default; explicit 0 = off (no infer).
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_sample: Option<usize>,
    /// Discover adapters de novo from a read sample, then trim with only the
    /// discovered set (ignores --adapter-preset for trimming; conflicts with
    /// --adapter-fasta). Off by default.
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_infer: bool,
    /// Discover adapters and print them (sequences + support + catalog names),
    /// then exit without trimming. Implies --adapter-infer. May be combined with
    /// --adapter-fasta (naming covers the built-in catalog plus your FASTA's
    /// adapters -- see the printed note).
    #[arg(long, help_heading = "Adapter trimming")]
    adapter_infer_only: bool,
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
enum FormatArg {
    Fastq,
    FastqGz,
    #[value(name = "fastq-bgz", alias = "fastq-bgzf")]
    FastqBgzf,
    Bam,
}

impl From<FormatArg> for Format {
    fn from(f: FormatArg) -> Self {
        match f {
            FormatArg::Fastq => Format::Fastq,
            FormatArg::FastqGz => Format::FastqGz,
            FormatArg::FastqBgzf => Format::FastqBgzf,
            FormatArg::Bam => Format::Bam,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
enum QualModeArg {
    Mean,
    Arithmetic,
    Median,
}

impl From<QualModeArg> for QualMode {
    fn from(m: QualModeArg) -> Self {
        match m {
            QualModeArg::Mean => QualMode::Mean,
            QualModeArg::Arithmetic => QualMode::Arithmetic,
            QualModeArg::Median => QualMode::Median,
        }
    }
}

#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
enum AdapterPresetArg {
    None,
    Ont,
}

pub fn parse() -> anyhow::Result<Config> {
    let c = Cli::parse();

    // Mutual exclusion of the three quality trim ops.
    let n_quality = [
        c.qual_trim.is_some(),
        c.qual_best_segment.is_some(),
        c.qual_split.is_some(),
    ]
    .iter()
    .filter(|&&b| b)
    .count();
    if n_quality > 1 {
        anyhow::bail!("--qual-trim, --qual-best-segment and --qual-split are mutually exclusive");
    }
    // bgzf (libdeflate) accepts up to 12 and gzip up to 9; cap at the common 0-9
    // so a single flag is valid for both compressed output formats.
    if c.compression_level > 9 {
        anyhow::bail!(
            "--compression-level must be between 0 and 9 (got {})",
            c.compression_level
        );
    }

    // Reject contradictory or out-of-domain filter bounds up front, rather than
    // silently keeping zero reads and exiting successfully.
    let max_length = c.max_length.unwrap_or(usize::MAX);
    if c.min_length > max_length {
        anyhow::bail!(
            "--min-length ({}) must not exceed --max-length ({max_length})",
            c.min_length
        );
    }
    if c.min_qual.is_nan() || c.max_qual.is_nan() {
        anyhow::bail!("--min-qual and --max-qual must be numbers (got NaN)");
    }
    if c.min_qual > c.max_qual {
        anyhow::bail!(
            "--min-qual ({}) must not exceed --max-qual ({})",
            c.min_qual,
            c.max_qual
        );
    }
    if let Some(g) = c.min_gc
        && !(0.0..=1.0).contains(&g)
    {
        anyhow::bail!("--min-gc ({g}) must be a fraction between 0 and 1");
    }
    if let Some(g) = c.max_gc
        && !(0.0..=1.0).contains(&g)
    {
        anyhow::bail!("--max-gc ({g}) must be a fraction between 0 and 1");
    }
    if let (Some(a), Some(b)) = (c.min_gc, c.max_gc)
        && a > b
    {
        anyhow::bail!("--min-gc ({a}) must not exceed --max-gc ({b})");
    }

    let quality = if let Some(q) = c.qual_trim {
        Some(QualityOp::TrimQual(q))
    } else if let Some(q) = c.qual_best_segment {
        Some(QualityOp::BestSegment(q))
    } else if let Some(q) = c.qual_split {
        Some(QualityOp::Split {
            cutoff: q,
            window: c.qual_split_window,
        })
    } else {
        None
    };
    let fastq_tags = FastqTags::parse(&c.fastq_tags)?;

    let adapter_infer = if c.adapter_infer_only {
        AdapterInfer::ReportOnly
    } else if c.adapter_infer {
        AdapterInfer::Trim
    } else {
        AdapterInfer::Off
    };

    // Mutual exclusion with an explicit FASTA (trim mode only; report-only
    // allows a FASTA so it can be cross-named against what inference finds).
    if adapter_infer == AdapterInfer::Trim && c.adapter_fasta.is_some() {
        anyhow::bail!(
            "--adapter-infer and --adapter-fasta are mutually exclusive (one discovers \
             the set, the other supplies it). To trim with your FASTA and also see what \
             inference finds, run --adapter-infer-only --adapter-fasta <file> first."
        );
    }
    // The preset is redundant for trimming under infer (inference builds its
    // own set) -- kept only so its names can be cross-referenced.
    if adapter_infer != AdapterInfer::Off && c.adapter_preset != AdapterPresetArg::None {
        eprintln!(
            "[WARN] --adapter-preset is ignored for trimming under --adapter-infer \
             (used only for naming discovered adapters)"
        );
    }
    // --adapter-infer-only + --adapter-fasta is allowed (unlike --adapter-infer,
    // which rejects a FASTA outright above): report-only names discovered
    // adapters against the built-in ONT catalog UNION the user's FASTA (see
    // `infer::discover`), so a user combining the two flags gets their own
    // adapter names surfaced too, not just catalog matches.
    if adapter_infer == AdapterInfer::ReportOnly && c.adapter_fasta.is_some() {
        eprintln!(
            "[INFO] --adapter-infer-only with --adapter-fasta: discovered adapters are named \
             against the built-in ONT catalog plus your FASTA's adapters"
        );
    }

    let mut adapter_seqs: Vec<crate::adapter::Adapter> = Vec::new();
    if c.adapter_preset == AdapterPresetArg::Ont {
        adapter_seqs.extend(crate::adapter::preset::preset_ont());
    }
    // Tracked separately from `adapter_seqs` (which also mixes in the preset
    // above): only the user's own FASTA entries are carried onward as extra
    // naming refs (see the `trim_adapters` comment below), not the preset
    // (already covered by `infer::discover`'s own built-in catalog lookup).
    let mut fasta_adapters: Vec<crate::adapter::Adapter> = Vec::new();
    if let Some(path) = &c.adapter_fasta {
        let from_fasta = read_adapter_fasta(path)?;
        if from_fasta.is_empty() {
            anyhow::bail!(
                "--adapter-fasta {}: no usable adapters (all entries were empty, \
                 shorter than the {}-bp minimum, or non-ACGT)",
                path.display(),
                crate::adapter::MIN_PATTERN_LEN
            );
        }
        fasta_adapters = from_fasta.clone();
        adapter_seqs.extend(from_fasta);
    }
    let adapters = if adapter_seqs.is_empty() && adapter_infer == AdapterInfer::Off {
        if c.adapter_ends_only {
            eprintln!(
                "[WARN] --adapter-ends-only has no effect without --adapter-fasta or --adapter-preset"
            );
        }
        None
    } else {
        if !(0.0..=1.0).contains(&c.adapter_error_rate) {
            anyhow::bail!(
                "--adapter-error-rate ({}) must be between 0 and 1",
                c.adapter_error_rate
            );
        }
        if c.adapter_end_size == 0 {
            anyhow::bail!("--adapter-end-size must be >= 1");
        }
        // Under infer, the trimming set is discovered later (discovery fills
        // it in); any preset sequences gathered above are ignored here (the
        // preset WARN above already told the user). A report-only FASTA is
        // carried through as `fasta_adapters` instead of being dropped: it's
        // never trimmed against under infer (discovery always replaces this
        // field before any dispatch -- see `maybe_reduce_adapters`, and
        // report-only exits before dispatch entirely), so reusing this field
        // to ferry the FASTA refs to `infer::discover` for cross-naming is
        // safe. Under `Trim`, a FASTA is rejected above, so this is always
        // empty there.
        let trim_adapters = if adapter_infer == AdapterInfer::Off {
            adapter_seqs
        } else {
            fasta_adapters
        };
        Some(crate::adapter::AdapterConfig {
            adapters: trim_adapters,
            error_rate: c.adapter_error_rate,
            end_size: c.adapter_end_size,
            split: !c.adapter_ends_only,
            candidate_index: std::sync::OnceLock::new(),
        })
    };

    // Resolve the sample size (distinguish unset from explicit 0): omitted
    // means "the mode default" (0 = off normally, 40000 under infer); an
    // explicit value is validated against the existing detection-floor rule,
    // plus a 0-under-infer rejection (0 would starve inference entirely).
    let adapter_sample = match c.adapter_sample {
        None => {
            if adapter_infer != AdapterInfer::Off {
                40_000
            } else {
                0
            }
        },
        Some(n) => {
            if n != 0 && n < crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION {
                anyhow::bail!(
                    "--adapter-sample ({}) must be 0 (disable detection) or at least {} \
                     (smaller samples are too few for reliable detection)",
                    n,
                    crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION
                );
            }
            if n == 0 && adapter_infer != AdapterInfer::Off {
                anyhow::bail!(
                    "--adapter-sample 0 disables sampling, which --adapter-infer requires; \
                     omit it or pass >= {}",
                    crate::adapter::detect::MIN_SAMPLE_FOR_DETECTION
                );
            }
            n
        },
    };

    let ncpu = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let threads = crate::config::resolve_threads(c.threads, ncpu);
    let threads_clamped = match c.threads {
        Some(n) if n > ncpu => Some((n, ncpu)),
        _ => None,
    };

    // Presence detection is preset-only: a user-supplied --adapter-fasta is a
    // curated set that should all be searched, and sampling could wrongly drop a
    // rare custom adapter. So detection is disabled whenever a FASTA is provided
    // and we're not inferring (under infer, --adapter-sample means the
    // inference buffer, not presence detection, so it's left alone).
    let adapter_sample = if adapter_infer == AdapterInfer::Off && c.adapter_fasta.is_some() {
        if adapter_sample > 0 {
            eprintln!(
                "[WARN] --adapter-sample is ignored with --adapter-fasta (presence detection is preset-only)"
            );
        }
        0
    } else {
        adapter_sample
    };

    Ok(Config {
        io: IoConfig {
            input: c.input,
            output: c.output,
            in_format: c.in_format.map(Into::into),
            out_format: c.out_format.map(Into::into),
        },
        filter: FilterConfig {
            min_length: c.min_length,
            max_length: c.max_length.unwrap_or(usize::MAX),
            min_qual: c.min_qual,
            max_qual: c.max_qual,
            min_gc: c.min_gc,
            max_gc: c.max_gc,
            qual_mode: c.qual_mode.into(),
        },
        trim: TrimPlan {
            head: c.head_crop,
            tail: c.tail_crop,
            quality,
        },
        adapters,
        adapter_infer,
        threads,
        fastq_tags,
        render_workers: 0,
        adapter_sample,
        compression_level: c.compression_level,
        update_moves: c.update_moves,
        verbosity: c.verbose,
        quiet: c.quiet,
        threads_clamped,
    })
}

/// Read adapter sequences from a FASTA. Lowercase acgt is uppercased and
/// accepted; entries with any other non-ACGT byte (e.g. IUPAC ambiguity
/// codes) are skipped with a warning, matching `preset::parse_catalog`'s
/// ACGT-only rule. Entries shorter than `adapter::MIN_PATTERN_LEN` (the
/// matcher's minimum pattern length) are also skipped with a warning, since a
/// shorter pattern would never be matched anyway.
fn read_adapter_fasta(path: &std::path::Path) -> anyhow::Result<Vec<crate::adapter::Adapter>> {
    use seq_io::fasta::{Reader, Record};
    let mut reader = Reader::from_path(path)
        .map_err(|e| anyhow::anyhow!("--adapter-fasta {}: {e}", path.display()))?;
    let mut out = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| anyhow::anyhow!("--adapter-fasta {}: {e}", path.display()))?;
        let seq: Vec<u8> = rec
            .seq()
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(u8::to_ascii_uppercase)
            .collect();
        let name = String::from_utf8_lossy(rec.head()).into_owned();
        if !seq.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            eprintln!("[WARN] adapter {name:?} has non-ACGT bases; skipped");
            continue;
        }
        if seq.len() < crate::adapter::MIN_PATTERN_LEN {
            eprintln!(
                "[WARN] adapter {name:?} is {} bp; shorter than the {}-bp minimum match length, skipped",
                seq.len(),
                crate::adapter::MIN_PATTERN_LEN
            );
            continue;
        }
        out.push(crate::adapter::Adapter {
            name,
            seq,
            end: crate::adapter::End::Both,
        });
    }
    Ok(out)
}

/// Build a Config directly (used by integration tests). head/tail are fixed crops.
#[doc(hidden)]
pub fn config_for_test(
    input: &std::path::Path,
    output: &std::path::Path,
    head_crop: usize,
    tail_crop: usize,
) -> Config {
    config_for_test_threads(input, output, head_crop, tail_crop, 1)
}

/// Same as `config_for_test`, but with an explicit thread count (used by tests
/// that need to exercise the parallel BAM dispatch, e.g. a `t8` oracle run).
#[doc(hidden)]
pub fn config_for_test_threads(
    input: &std::path::Path,
    output: &std::path::Path,
    head_crop: usize,
    tail_crop: usize,
    threads: usize,
) -> Config {
    Config {
        io: IoConfig {
            input: Some(input.to_path_buf()),
            output: Some(output.to_path_buf()),
            in_format: Some(Format::Bam),
            out_format: Some(Format::Bam),
        },
        filter: FilterConfig {
            min_length: 1,
            max_length: usize::MAX,
            min_qual: 0.0,
            max_qual: 1000.0,
            min_gc: None,
            max_gc: None,
            qual_mode: QualMode::Mean,
        },
        trim: TrimPlan {
            head: head_crop,
            tail: tail_crop,
            quality: None,
        },
        adapters: None,
        adapter_infer: crate::config::AdapterInfer::Off,
        threads: threads.max(1),
        fastq_tags: FastqTags::All,
        render_workers: 0,
        adapter_sample: 0,
        compression_level: 6,
        update_moves: false,
        verbosity: 0,
        quiet: true,
        threads_clamped: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    // FASTA loading and adapter search must enforce the same minimum length.
    #[test]
    fn read_adapter_fasta_skips_entries_below_min_pattern_len() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("adapters.fasta");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, ">kept_20bp").unwrap();
        writeln!(f, "ACGTACGTACGTACGTACGT").unwrap(); // 20 bp
        writeln!(f, ">skipped_8bp").unwrap();
        writeln!(f, "ACGTACGT").unwrap(); // 8 bp, below the 11bp MIN_PATTERN_LEN
        drop(f);

        let adapters = read_adapter_fasta(&path).unwrap();

        assert_eq!(adapters.len(), 1);
        assert_eq!(adapters[0].seq, b"ACGTACGTACGTACGTACGT".to_vec());
    }

    // Parity with `preset::parse_catalog`: entries with non-ACGT bases (IUPAC
    // ambiguity codes like N) are rejected, not silently searched as-is.
    #[test]
    fn read_adapter_fasta_skips_non_acgt_entries() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("adapters.fasta");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, ">kept_valid_20bp").unwrap();
        writeln!(f, "ACGTACGTACGTACGTACGT").unwrap(); // 20 bp, valid ACGT
        writeln!(f, ">skipped_n_20bp").unwrap();
        writeln!(f, "ACGTACGTACGTACGTACGN").unwrap(); // 20 bp, but has an N
        drop(f);

        let adapters = read_adapter_fasta(&path).unwrap();

        assert_eq!(adapters.len(), 1);
        assert_eq!(adapters[0].name, "kept_valid_20bp");
        assert_eq!(adapters[0].seq, b"ACGTACGTACGTACGTACGT".to_vec());
    }
}
