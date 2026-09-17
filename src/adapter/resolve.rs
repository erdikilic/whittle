//! Resolution of the adapter set a run trims against.
//!
//! Presence detection and ab-initio inference both read a sample before the set
//! is final. This module buffers a prefix, narrows or discovers the set, and
//! returns the record stream with the sampled prefix chained in front of it.

use std::borrow::Cow;

use super::{detect, infer};
use crate::config::{AdapterInfer, Config};

/// Decodes the packed `SEQ` of a lazy raw BAM record for adapter sampling.
/// Workflow records otherwise stay packed until a render worker converts them
/// to `RecordBuf`.
pub(crate) fn bam_seq(rec: &noodles_bam::Record) -> Cow<'_, [u8]> {
    Cow::Owned(rec.sequence().iter().collect())
}

/// Support below which a kept adapter is logged with a warning rather than a
/// plain info line. It is three times `infer::KEEP_SUPPORT` (0.15): a
/// barcode-specific sequence present in a fraction of reads can clear the keep
/// floor while staying far from the near-1.0 support of a library adapter.
pub(crate) const MARGINAL_SUPPORT: f64 = 0.45;

/// Logs each ab-initio discovery: one `info!` line per adapter with its support
/// and best catalog match (an annotation; `inferred_N` is the name), a `warn!`
/// when the support is below `MARGINAL_SUPPORT`,
/// and the sequences at `debug!`.
pub(crate) fn log_discovered(discovered: &[infer::InferredAdapter], n_sampled: usize) {
    tracing::info!(
        reads = n_sampled,
        discovered = discovered.len(),
        "Adapter inference: sampled prefix scanned"
    );
    for d in discovered {
        let support = format!("{:.2}", d.support);
        match d.name_hits.first() {
            Some((name, pct)) => {
                let identity_pct = format!("{pct:.0}");
                tracing::info!(
                    adapter = %d.adapter.name,
                    catalog_match = %name,
                    identity_pct = %identity_pct,
                    support = %support,
                    "Inferred adapter"
                );
            },
            None => {
                tracing::info!(
                    adapter = %d.adapter.name,
                    support = %support,
                    "Inferred adapter with no catalog match"
                );
            },
        }
        if d.support < MARGINAL_SUPPORT {
            tracing::warn!(
                adapter = %d.adapter.name,
                support = %support,
                floor = MARGINAL_SUPPORT,
                "Inferred adapter support is marginal; verify with --discover-adapters report"
            );
        }
        if d.uncertain_bases() > 0 {
            tracing::warn!(
                adapter = %d.adapter.name,
                anchor_bp = d.adapter.seq.len(),
                uncertain_bp = d.uncertain_bases(),
                consensus_bp = d.assembled_seq.len(),
                "Inferred adapter excludes unsupported consensus bases"
            );
        }
        let sequence = String::from_utf8_lossy(&d.adapter.seq);
        tracing::debug!(
            adapter = %d.adapter.name,
            sequence = %sequence,
            "Inferred adapter trimming sequence"
        );
        if d.uncertain_bases() > 0 {
            let consensus = String::from_utf8_lossy(&d.assembled_seq);
            tracing::debug!(
                adapter = %d.adapter.name,
                consensus = %consensus,
                "Inferred adapter full recurrent consensus"
            );
        }
    }
}

/// Prints inferred adapters as FASTA with support and the best catalog match.
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
                "bounded"
            },
            d.assembled_seq.len(),
            d.uncertain_bases(),
        );
        println!("{}", String::from_utf8_lossy(&d.adapter.seq));
    }
}

/// Runs `f` over the sampled sequences as plain slices. The decoded views are
/// materialized once here, so detection and inference share one borrow shape.
fn with_sequences<R, F, T>(sample: &[R], seq_of: &F, f: impl FnOnce(&[&[u8]]) -> T) -> T
where
    F: for<'a> Fn(&'a R) -> Cow<'a, [u8]>,
{
    let storage: Vec<Cow<'_, [u8]>> = sample.iter().map(seq_of).collect();
    let seqs: Vec<&[u8]> = storage.iter().map(|s| s.as_ref()).collect();
    f(&seqs)
}

/// Maximum retained sample payload, including decoded BAM sampling sequences.
const SAMPLE_BYTES: usize = 256 * 1024 * 1024;
/// Maximum bases retained for adapter sampling.
const SAMPLE_BASES: usize = 64 * 1024 * 1024;

/// Estimates a raw BAM record's retained payload and decoded sampling sequence.
pub(crate) fn bam_sample_weight(rec: &noodles_bam::Record) -> (usize, usize) {
    let bases = rec.sequence().len();
    let bytes = std::mem::size_of::<noodles_bam::Record>()
        + rec.name().map_or(0, |n| n.len() + 1)
        + rec.cigar().len() * 4
        + bases.div_ceil(2)
        + rec.quality_scores().len()
        + rec.data().as_bytes().len()
        + bases;
    (bytes, bases)
}

/// Buffers a prefix bounded by record count, payload bytes and bases. The last
/// record is retained whole, so one record can exceed a payload limit.
fn buffer_prefix<R>(
    records: &mut impl Iterator<Item = anyhow::Result<R>>,
    n: usize,
    weight: &impl Fn(&R) -> (usize, usize),
    max_bytes: usize,
    max_bases: usize,
) -> anyhow::Result<Vec<R>> {
    let mut sample = Vec::new();
    let (mut bytes, mut bases) = (0usize, 0usize);
    for _ in 0..n {
        match records.next() {
            Some(Ok(r)) => {
                let (record_bytes, record_bases) = weight(&r);
                bytes = bytes.saturating_add(record_bytes);
                bases = bases.saturating_add(record_bases);
                sample.push(r);
                if bytes >= max_bytes || bases >= max_bases {
                    tracing::info!(
                        reads = sample.len(),
                        requested = n,
                        bytes,
                        bases,
                        "Adapter sample reached its payload limit"
                    );
                    break;
                }
            },
            Some(Err(e)) => return Err(e),
            None => break,
        }
    }
    Ok(sample)
}

/// The outcome of resolution: the stream to process, with any sampled prefix
/// chained back in front of it, and the adapter set to trim it against.
pub(crate) struct Resolved<R> {
    /// The record stream, with any sampled prefix chained back in front of it.
    pub records: Box<dyn Iterator<Item = anyhow::Result<R>> + Send>,
    /// The set trimmed against, after presence detection narrowed the
    /// configured set or inference replaced it. `None` when trimming is off.
    pub adapters: Option<super::AdapterConfig>,
}

/// Decides the adapter set for a run, reading a prefix of the stream when
/// needed, and returns it with the stream intact.
///
/// `Ok(None)` means the run is over without writing records: that is
/// `--discover-adapters report`, which prints the inferred FASTA and stops.
///
/// Takes `&Config` and returns the outcome rather than writing back into the
/// config: the set is final only after reads have been seen, which is after the
/// banner has printed the configured set, and an in-place overwrite would leave
/// the banner and the run describing different sets.
pub(crate) fn resolve<R, I, F>(
    mut records: I,
    cfg: &Config,
    // The returned sequence view borrows the record passed to `seq_of`.
    seq_of: F,
    weight: impl Fn(&R) -> (usize, usize),
) -> anyhow::Result<Option<Resolved<R>>>
where
    // Workflow iterators are boxed and may cross worker-thread boundaries.
    I: Iterator<Item = anyhow::Result<R>> + Send + 'static,
    R: Send + 'static,
    F: for<'a> Fn(&'a R) -> Cow<'a, [u8]>,
{
    if cfg.adapter_infer != AdapterInfer::Off {
        // `cli::parse` pairs inference with an initially empty adapter config,
        // but `run` is public, so a library caller can omit it. The mismatch is
        // reported as an error rather than a panic.
        let Some(base) = cfg.adapters.clone() else {
            anyhow::bail!("adapter inference requires an adapter configuration");
        };

        let sample: Vec<R> = buffer_prefix(
            &mut records,
            cfg.adapter_sample,
            &weight,
            SAMPLE_BYTES,
            SAMPLE_BASES,
        )?;
        let s = sample.len();
        let chain =
            |sample: Vec<R>, records: I| -> Box<dyn Iterator<Item = anyhow::Result<R>> + Send> {
                Box::new(sample.into_iter().map(anyhow::Ok).chain(records))
            };
        if s < detect::MIN_SAMPLE_FOR_DETECTION {
            // Report-only mode writes no output when the sample is too small.
            tracing::warn!(
                reads = s,
                minimum = detect::MIN_SAMPLE_FOR_DETECTION,
                "Adapter inference: too few reads to infer reliably; keeping reads untrimmed"
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

        let discovered = with_sequences(&sample, &seq_of, |seqs| infer::discover(seqs, &base));
        log_discovered(&discovered, s);

        if cfg.adapter_infer.is_report() {
            // Report mode prints the inferred FASTA and writes no records.
            // `lib::settle` warns about every write target this leaves unused.
            print_discovered_fasta(&discovered);
            return Ok(None);
        }

        if discovered.is_empty() {
            tracing::warn!(
                reads = s,
                "Adapter inference: no adapters inferred from the sampled prefix; keeping \
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

    // No buffering when neither inference nor presence sampling is active.
    let Some(ac) = cfg.adapters.clone().filter(|_| cfg.adapter_sample > 0) else {
        return Ok(Some(Resolved {
            records: Box::new(records),
            adapters: cfg.adapters.clone(),
        }));
    };

    // Presence detection narrows the configured set to what the sampled prefix
    // contains.
    let sample: Vec<R> = buffer_prefix(
        &mut records,
        cfg.adapter_sample,
        &weight,
        SAMPLE_BYTES,
        SAMPLE_BASES,
    )?;
    let s = sample.len();
    let full = ac.adapters.len();
    let kept = if s < detect::MIN_SAMPLE_FOR_DETECTION {
        tracing::info!(
            reads = s,
            minimum = detect::MIN_SAMPLE_FOR_DETECTION,
            configured = full,
            "Adapter presence: sample too small; using all configured adapters"
        );
        ac.adapters.clone()
    } else {
        let detected = with_sequences(&sample, &seq_of, |seqs| {
            detect::present(seqs, &ac, detect::presence_min(s), cfg.threads)
        });
        if detected.is_empty() {
            tracing::warn!(
                reads = s,
                configured = full,
                "Adapter presence: no adapters detected in the sampled prefix; using all \
                 configured adapters (the prefix may be unrepresentative; --adapter-sample-reads 0 \
                 skips sampling)"
            );
            ac.adapters.clone()
        } else {
            let names: Vec<&str> = detected.iter().take(12).map(|a| a.name.as_str()).collect();
            let listed = names.join(", ");
            let unlisted = detected.len().saturating_sub(names.len());
            tracing::info!(
                reads = s,
                kept = detected.len(),
                configured = full,
                adapters = %listed,
                unlisted,
                "Adapter presence: sampled prefix narrowed the adapter set"
            );
            detected
        }
    };
    let mut reduced = ac;
    reduced.replace_adapters(kept);
    Ok(Some(Resolved {
        records: Box::new(sample.into_iter().map(anyhow::Ok).chain(records)),
        adapters: Some(reduced),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_limits_preserve_the_complete_record_stream() {
        for (max_bytes, max_bases, expected) in [(60, 100, 3), (100, 12, 3), (1, 1, 1)] {
            let mut records = (0..10usize).map(Ok);
            let sample =
                buffer_prefix(&mut records, 10, &|_| (20, 4), max_bytes, max_bases).unwrap();
            assert_eq!(sample.len(), expected);
            let replayed: Vec<_> = sample
                .into_iter()
                .map(Ok)
                .chain(records)
                .collect::<anyhow::Result<_>>()
                .unwrap();
            assert_eq!(replayed, (0..10).collect::<Vec<_>>());
        }
    }

    #[test]
    fn sample_limit_accounts_for_large_auxiliary_payloads() {
        let rec = crate::record::ReadRecord {
            name: vec![b'A'; 1000],
            seq: vec![b'C'; 4],
            qual: vec![40; 4],
        };
        let mut records = std::iter::repeat_with(|| Ok(rec.clone())).take(20);
        let sample = buffer_prefix(
            &mut records,
            20,
            &|r| (r.name.len() + r.seq.len() + r.qual.len(), r.seq.len()),
            2000,
            1000,
        )
        .unwrap();
        assert_eq!(sample.len(), 2);
        assert_eq!(records.count(), 18);
    }
}
