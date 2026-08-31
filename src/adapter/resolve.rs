//! Resolving the adapter set a run actually trims against.
//!
//! Presence detection and ab-initio inference both need to see reads before the
//! set is final, so this owns the buffer-then-decide seam: sample a prefix,
//! narrow or discover the set, then hand the record stream back to the caller
//! with the sampled prefix chained in front of it.

use std::borrow::Cow;

use super::{detect, infer};
use crate::config::{AdapterInfer, Config};

/// Decode the packed SEQ of a lazy raw BAM record only when adapter sampling
/// needs it. Normal workflow records stay packed until a render worker converts
/// them to `RecordBuf`.
pub(crate) fn bam_seq(rec: &noodles_bam::Record) -> Cow<'_, [u8]> {
    Cow::Owned(rec.sequence().iter().collect())
}

/// A kept adapter's support below this is close enough to `infer::KEEP_SUPPORT`
/// (0.30) to warrant a warning rather than a plain info line: a barcode-specific
/// sequence present in only a fraction of reads can clear the keep floor while
/// staying far from a confident near-1.0 presence. ~1.5x the floor gives headroom
/// without reaching a genuine high-prevalence adapter's typical support.
pub(crate) const MARGINAL_SUPPORT: f64 = 0.45;

/// Log each ab-initio discovery at `info!`: `inferred_N ≈ NAME (pct%) · support
/// X.XX`, or `(no catalog match)` when the sequence cross-names against nothing
/// in the ONT catalog. `N` is the 1-based position in `discovered`'s own order,
/// which agrees with the `inferred_{N}` name fallback. The raw sequence goes to
/// `debug!` instead, too noisy for INFO. Support below `MARGINAL_SUPPORT` also
/// gets a `warn!`, being close enough to the `KEEP_SUPPORT` floor to re-check.
pub(crate) fn log_discovered(discovered: &[infer::InferredAdapter], n_sampled: usize) {
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
                "Adapter '{}' support {:.2} is marginal (near the KEEP_SUPPORT floor); \
                 verify with --adapter-infer report",
                d.adapter.name,
                d.support
            );
        }
        if d.uncertain_bases() > 0 {
            tracing::warn!(
                "Adapter '{}' uses a conservative {} bp terminal anchor; {} bp of the \
                 {} bp recurrent consensus remain uncertain and will not be trimmed \
                 (--adapter-infer-policy aggressive opts into the full consensus)",
                d.adapter.name,
                d.adapter.seq.len(),
                d.uncertain_bases(),
                d.assembled_seq.len(),
            );
        }
        tracing::debug!(
            "inferred_{n} trimming sequence: {}",
            String::from_utf8_lossy(&d.adapter.seq)
        );
        if d.uncertain_bases() > 0 {
            tracing::debug!(
                "inferred_{n} full recurrent consensus (review only): {}",
                String::from_utf8_lossy(&d.assembled_seq)
            );
        }
    }
}

/// Print inferred adapters as FASTA with support and the best catalog match.
/// Numbering follows the final discovery order used by the status log.
pub(crate) fn print_discovered_fasta(discovered: &[infer::InferredAdapter]) {
    for (i, d) in discovered.iter().enumerate() {
        let n = i + 1;
        let name_suffix = match d.name_hits.first() {
            Some((name, pct)) => format!(" [\u{2248} {name} ({pct:.0}%)]"),
            None => String::new(),
        };
        println!(
            ">inferred_{n} support={:.2} boundary={} assembled_length={} uncertain_bases={}{name_suffix}",
            d.support,
            if d.uncertain_bases() == 0 {
                "full"
            } else {
                "conservative"
            },
            d.assembled_seq.len(),
            d.uncertain_bases(),
        );
        println!("{}", String::from_utf8_lossy(&d.adapter.seq));
    }
}

/// Buffer at most `n` records, stopping when the input is exhausted.
pub(crate) fn buffer_prefix<R>(
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

/// What resolution decided: the stream to process, with any sampled prefix
/// chained back in front of it, and the adapter set to trim it against.
pub(crate) struct Resolved<R> {
    pub records: Box<dyn Iterator<Item = anyhow::Result<R>> + Send>,
    /// The set actually trimmed against, after presence detection narrowed the
    /// configured set or inference replaced it. `None` when trimming is off.
    pub adapters: Option<super::AdapterConfig>,
}

/// Decide the adapter set for a run, reading a prefix of the stream when it has
/// to, and return it with the stream intact.
///
/// `Ok(None)` means the run is over without writing records: that is
/// `--adapter-infer report`, which prints the inferred FASTA and stops.
///
/// Takes `&Config` and returns the outcome rather than writing back into the
/// config: the set is only final after reads have been seen, which is well after
/// the banner has printed the configured one, and an in-place overwrite made
/// those two silently different views of the same field.
pub(crate) fn resolve<R, I, F>(
    mut records: I,
    cfg: &Config,
    // The returned sequence view borrows the record passed to `seq_of`.
    seq_of: F,
) -> anyhow::Result<Option<Resolved<R>>>
where
    // Workflow iterators are boxed and may cross worker-thread boundaries.
    I: Iterator<Item = anyhow::Result<R>> + Send + 'static,
    R: Send + 'static,
    F: for<'a> Fn(&'a R) -> Cow<'a, [u8]>,
{
    if cfg.adapter_infer != AdapterInfer::Off {
        // Inference mode stores an empty configuration until discovery completes.
        // `cli::parse` always pairs inference with an (initially empty) adapter
        // config, but `run` is public, so a library caller can hand over a
        // mismatched pair. An error beats a panic at a public boundary.
        let Some(base) = cfg.adapters.clone() else {
            anyhow::bail!(
                "adapter inference was requested without an adapter configuration; \
                 this is a caller error, not a bad input file"
            );
        };

        let sample: Vec<R> = buffer_prefix(&mut records, cfg.adapter_sample)?;
        let s = sample.len();
        let chain =
            |sample: Vec<R>, records: I| -> Box<dyn Iterator<Item = anyhow::Result<R>> + Send> {
                Box::new(sample.into_iter().map(anyhow::Ok).chain(records))
            };
        if s < detect::MIN_SAMPLE_FOR_DETECTION {
            // Report-only mode must not create output when the sample is too small.
            tracing::warn!(
                "Adapter inference: too few reads ({s}, need >= {}) to infer reliably; \
                 keeping reads untrimmed",
                detect::MIN_SAMPLE_FOR_DETECTION
            );
            if cfg.adapter_infer.is_report() {
                return Ok(None);
            }
            let mut reduced = base;
            reduced.replace_adapters(Vec::new());
            return Ok(Some(Resolved {
                records: chain(sample, records),
                adapters: Some(reduced),
            }));
        }

        let seq_storage: Vec<Cow<'_, [u8]>> = sample.iter().map(&seq_of).collect();
        let seqs: Vec<&[u8]> = seq_storage.iter().map(|s| s.as_ref()).collect();
        let discovered =
            infer::discover_with_policy(&seqs, &base, cfg.adapter_infer.is_aggressive());
        log_discovered(&discovered, s);

        if cfg.adapter_infer.is_report() {
            // Report mode prints the inferred FASTA and writes no records, so
            // there is nothing for `-o` to hold and no counters worth
            // summarizing. Both flags are named explicitly: exiting 0 having
            // silently created neither file strands a pipeline that expected one.
            for (flag, given) in [
                ("-o/--output", cfg.io.output.is_some()),
                ("--summary-json", cfg.summary_json.is_some()),
            ] {
                if given {
                    tracing::warn!(
                        "{flag} is ignored under --adapter-infer report, which writes no records"
                    );
                }
            }
            print_discovered_fasta(&discovered);
            return Ok(None);
        }

        if discovered.is_empty() {
            tracing::warn!(
                "Adapter inference: no adapters inferred from the first {s} reads; keeping \
                 reads untrimmed"
            );
        }
        let mut reduced = base;
        reduced.replace_adapters(discovered.into_iter().map(|d| d.adapter).collect());
        return Ok(Some(Resolved {
            records: chain(sample, records),
            adapters: Some(reduced),
        }));
    }

    // Avoid buffering when neither inference nor presence sampling is active.
    if cfg.adapters.is_none() || cfg.adapter_sample == 0 {
        return Ok(Some(Resolved {
            records: Box::new(records),
            adapters: cfg.adapters.clone(),
        }));
    }

    // Narrow the configured set to what the sampled prefix actually contains.
    let mut sample: Vec<R> = Vec::new();
    let mut adapters = cfg.adapters.clone();
    if let Some(ac) = cfg.adapters.clone() {
        sample = buffer_prefix(&mut records, cfg.adapter_sample)?;
        let s = sample.len();
        let full = ac.adapters.len();
        let kept = if s < detect::MIN_SAMPLE_FOR_DETECTION {
            tracing::info!(
                "Adapter presence: only {s} reads (< {}); using all {full} adapters",
                detect::MIN_SAMPLE_FOR_DETECTION
            );
            ac.adapters.clone()
        } else {
            let seq_storage: Vec<Cow<'_, [u8]>> = sample.iter().map(&seq_of).collect();
            let seqs: Vec<&[u8]> = seq_storage.iter().map(|s| s.as_ref()).collect();
            let detected = detect::present(
                &seqs,
                &ac.adapters,
                ac.error_rate,
                ac.end_size,
                ac.split,
                detect::presence_min(s),
                cfg.threads,
            );
            if detected.is_empty() {
                tracing::warn!(
                    "Adapter presence: no adapters detected in the first {s} sampled reads; using all {full} \
                     (the sampled prefix may be unrepresentative; pass --adapter-sample 0 to always use the full set)"
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
        adapters = Some(reduced);
    }
    Ok(Some(Resolved {
        records: Box::new(sample.into_iter().map(anyhow::Ok).chain(records)),
        adapters,
    }))
}
