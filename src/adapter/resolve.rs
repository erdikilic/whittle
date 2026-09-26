//! Resolution of the adapter set a run trims against.
//!
//! Presence detection and de novo inference both read a sample before the set
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

/// Support below which a discovered sequencing adapter is logged with a
/// warning rather than a plain info line. Sparse families warrant review even
/// when their consensus clears the minimum discovery support. End-only
/// layers such as barcodes are sparse by design and are not warned about.
pub(crate) const MARGINAL_SUPPORT: f64 = 0.45;

/// Logs each de novo discovery: one `info!` line per adapter with its support
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
                    layer = d.layer + 1,
                    role = d.adapter.role.label(),
                    catalog_match = %name,
                    identity_pct = %identity_pct,
                    support = %support,
                    "Inferred adapter"
                );
            },
            None => {
                tracing::info!(
                    adapter = %d.adapter.name,
                    layer = d.layer + 1,
                    role = d.adapter.role.label(),
                    support = %support,
                    "Inferred adapter with no catalog match"
                );
            },
        }
        if d.adapter.role == crate::adapter::Role::Adapter && d.support < MARGINAL_SUPPORT {
            tracing::warn!(
                adapter = %d.adapter.name,
                support = %support,
                floor = MARGINAL_SUPPORT,
                "Inferred adapter support is marginal; verify with --adapter-report"
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
            ">inferred_{n} layer={} role={} support={:.2} boundary={} assembled_length={} uncertain_bases={}{name_suffix}",
            d.layer + 1,
            d.adapter.role.label(),
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
/// `--adapter-report`, which prints the inferred FASTA and stops.
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
            if base.adapters.is_empty() {
                tracing::warn!(
                    reads = s,
                    "Adapter inference: no adapters inferred from the sampled prefix; keeping \
                     reads untrimmed"
                );
            } else {
                tracing::info!(
                    reads = s,
                    "Adapter inference: no sequences beyond the configured set"
                );
            }
        }
        // Known sequences stay in the set; discoveries extend it. Discovery
        // starts behind every preset entry, while the preset's barcode
        // entries that the sampled reads do not carry leave the set trimmed
        // against. A supplied FASTA is searched in full.
        let mut reduced = if cfg.adapter_fasta.is_none() {
            without_absent_barcodes(&sample, &seq_of, base, cfg.threads)
        } else {
            base
        };
        let primers = library_primers(&reduced.adapters, &discovered, reduced.error_rate);
        let mut adapters = reduced.adapters.clone();
        adapters.extend(discovered.into_iter().map(|d| d.adapter));
        reduced.replace_adapters(adapters);
        judge_amplicon(&sample, &seq_of, &mut reduced, &primers, cfg.threads);
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
    let primers = library_primers(&kept, &[], ac.error_rate);
    let mut reduced = ac;
    reduced.replace_adapters(kept);
    judge_amplicon(&sample, &seq_of, &mut reduced, &primers, cfg.threads);
    Ok(Some(Resolved {
        records: Box::new(sample.into_iter().map(anyhow::Ok).chain(records)),
        adapters: Some(reduced),
    }))
}

/// Returns the sequences that may be primers of the library in
/// `detect::primer_share`: the catalog marker-gene primers, the primer-role
/// entries of `known` that are marker primers or no catalog sequence, and
/// the sequences of `discovered`, other than the members of a barcode layer,
/// whose best catalog match, if any, is a marker primer. The catalog's other
/// sequences (adapters, barcodes, flanks and the PCR handles of other kits)
/// are technical layers that a library of any kind carries at its read ends,
/// and a barcode that discovery finds as a layer of its own is named after
/// its catalog entry. A discovered sequence of the adapter role is included,
/// as a primer that opens the reads flush with their ends is one.
fn library_primers(
    known: &[super::Adapter],
    discovered: &[infer::InferredAdapter],
    error_rate: f64,
) -> Vec<Vec<u8>> {
    let catalog = super::preset::preset(super::preset::Kit::ALL);
    let marker = |seq: &[u8]| super::matches_marker_primer(seq, error_rate);
    let catalog_seq = |seq: &[u8]| {
        catalog
            .iter()
            .any(|entry| entry.seq == seq || entry.seq == super::reverse_complement(seq))
    };
    let mut primers: Vec<Vec<u8>> = super::catalog::MARKER_PRIMERS
        .iter()
        .map(|primer| primer.to_vec())
        .collect();
    for entry in known.iter().filter(|a| a.role == super::Role::Primer) {
        let seq = entry.seq.to_ascii_uppercase();
        if marker(&seq) || !catalog_seq(&seq) {
            primers.push(seq);
        }
    }
    for found in discovered
        .iter()
        .filter(|d| d.adapter.role != super::Role::Barcode)
    {
        let named = found
            .name_hits
            .first()
            .and_then(|(name, _)| catalog.iter().find(|entry| entry.name == *name));
        if named.is_none_or(|entry| marker(&entry.seq)) {
            primers.push(found.adapter.seq.to_ascii_uppercase());
        }
    }
    primers.sort();
    primers.dedup();
    primers
}

/// Marks `ac` as an amplicon library (`AdapterConfig::amplicon`) when
/// `primers` (`library_primers`) open at least `detect::AMPLICON_PERCENT` of
/// the sampled reads (`detect::primer_share`) and the median read shares its
/// length with at least `detect::AMPLICON_LENGTH_PERCENT` of them
/// (`detect::length_share`), and logs both shares whenever they are
/// measured. Only the reads the shares examine, spread evenly over the
/// sample, are decoded.
fn judge_amplicon<R, F>(
    sample: &[R],
    seq_of: &F,
    ac: &mut super::AdapterConfig,
    primers: &[Vec<u8>],
    threads: usize,
) where
    F: for<'a> Fn(&'a R) -> Cow<'a, [u8]>,
{
    let step = sample.len().div_ceil(detect::AMPLICON_SAMPLE_READS).max(1);
    let storage: Vec<Cow<'_, [u8]>> = sample.iter().step_by(step).map(seq_of).collect();
    let seqs: Vec<&[u8]> = storage.iter().map(|s| s.as_ref()).collect();
    let Some(share) = detect::primer_share(&seqs, ac, primers, threads) else {
        return;
    };
    let lengths = detect::length_share(&seqs);
    let amplicon = share * 100.0 >= detect::AMPLICON_PERCENT as f64
        && lengths * 100.0 >= detect::AMPLICON_LENGTH_PERCENT as f64;
    let percent = format!("{:.1}", share * 100.0);
    let length_percent = format!("{:.1}", lengths * 100.0);
    if amplicon {
        tracing::info!(
            reads = sample.len(),
            percent = %percent,
            length_percent = %length_percent,
            "Adapter presence: primers open most sampled reads, which share a few lengths; \
             interior marker-gene primers split reads as amplicon junctions"
        );
    } else {
        tracing::info!(
            reads = sample.len(),
            percent = %percent,
            length_percent = %length_percent,
            "Adapter presence: no amplicon library; interior marker-gene primers split reads \
             only beside a junction partner"
        );
    }
    ac.set_amplicon(amplicon);
}

/// Removes from `ac` the barcode entries that presence detection
/// (`detect::present`) does not find in the sampled reads, in two steps.
/// First the absent entries that bound a barcode (`bounds_barcode`) go: the
/// barcode constructs and the flanks between a barcode and the insert. A
/// barcode panel leaves its junctions to them, so one that the reads do not
/// carry would keep the panel from splitting. When none of them remains, the
/// panel splits its own junctions, and its absent members and the other
/// absent barcode-role entries go too, since a splitting entry takes part in
/// the interior chance bound of the set. The other adapters and primers are
/// kept. A sample below `detect::MIN_SAMPLE_FOR_DETECTION`, or a set without
/// barcode entries, is returned as given.
fn without_absent_barcodes<R, F>(
    sample: &[R],
    seq_of: &F,
    ac: super::AdapterConfig,
    threads: usize,
) -> super::AdapterConfig
where
    F: for<'a> Fn(&'a R) -> Cow<'a, [u8]>,
{
    let s = sample.len();
    let barcoded = ac
        .adapters
        .iter()
        .any(|a| a.role == super::Role::Barcode || super::bounds_barcode(a));
    if s < detect::MIN_SAMPLE_FOR_DETECTION || !barcoded {
        return ac;
    }
    let detected = with_sequences(sample, seq_of, |seqs| {
        detect::present(seqs, &ac, detect::presence_min(s), threads)
    });
    let mut kept: Vec<super::Adapter> = ac
        .adapters
        .iter()
        .filter(|a| !super::bounds_barcode(a) || detected.contains(a))
        .cloned()
        .collect();
    if !kept.iter().any(super::bounds_barcode) {
        kept.retain(|a| a.role != super::Role::Barcode || detected.contains(a));
    }
    if kept.len() < ac.adapters.len() {
        tracing::info!(
            reads = s,
            dropped = ac.adapters.len() - kept.len(),
            configured = ac.adapters.len(),
            "Adapter presence: barcode entries absent from the sampled prefix are not searched"
        );
    }
    let mut reduced = ac;
    reduced.replace_adapters(kept);
    reduced
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

    /// Deterministic bases from a SplitMix64 stream.
    fn bases(mut state: u64, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| {
                state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                b"ACGT"[((z ^ (z >> 31)) >> 62) as usize]
            })
            .collect()
    }

    /// Builds a set at error rate 0.2 from `(name, sequence, role)` entries.
    fn set(entries: Vec<(&str, Vec<u8>, crate::adapter::Role)>) -> crate::adapter::AdapterConfig {
        crate::adapter::AdapterConfig {
            adapters: entries
                .into_iter()
                .map(|(name, seq, role)| crate::adapter::Adapter {
                    name: name.into(),
                    seq,
                    role,
                })
                .collect(),
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 20,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
        }
    }

    /// Returns the names of the entries of `ac` kept for 200 reads that each
    /// start with `lead`.
    fn kept_names(ac: crate::adapter::AdapterConfig, lead: &[u8]) -> Vec<String> {
        let reads: Vec<Vec<u8>> = (0..200)
            .map(|i| [lead, &bases(100 + i, 600)].concat())
            .collect();
        without_absent_barcodes(&reads, &|r: &Vec<u8>| Cow::Borrowed(&r[..]), ac, 1)
            .adapters
            .into_iter()
            .map(|a| a.name)
            .collect()
    }

    /// A barcode construct the reads do not carry is dropped, and with no
    /// flank or construct left the absent barcodes go too; every other
    /// entry, present or not, is kept.
    #[test]
    fn absent_barcode_bounds_and_then_absent_barcodes_are_dropped() {
        use crate::adapter::Role;
        let adapter = bases(1, 30);
        let present = bases(2, 24);
        let mut construct = bases(3, 20);
        construct.extend(std::iter::repeat_n(b'N', 24));
        construct.extend(bases(4, 8));
        let ac = set(vec![
            ("adapter", adapter.clone(), Role::Adapter),
            ("absent_adapter", bases(5, 30), Role::Adapter),
            ("construct", construct, Role::Barcode),
            ("present", present.clone(), Role::Barcode),
            ("absent", bases(6, 24), Role::Barcode),
        ]);
        assert_eq!(
            kept_names(ac, &[adapter, present].concat()),
            ["adapter", "absent_adapter", "present"]
        );
    }

    /// While a barcode flank the reads carry remains, the barcodes leave
    /// their junctions to it and the absent ones are kept; an absent
    /// construct is still dropped.
    #[test]
    fn absent_barcodes_stay_beside_a_present_flank() {
        use crate::adapter::Role;
        let flank = b"CCATATCCGTGTCGCCCTT".to_vec();
        let present = bases(2, 24);
        let mut construct = bases(3, 20);
        construct.extend(std::iter::repeat_n(b'N', 24));
        construct.extend(bases(4, 8));
        let ac = set(vec![
            ("construct", construct, Role::Barcode),
            ("present", present.clone(), Role::Barcode),
            ("absent", bases(6, 24), Role::Barcode),
            ("flank", flank.clone(), Role::Barcode),
        ]);
        assert_eq!(
            kept_names(ac, &[present, flank].concat()),
            ["present", "absent", "flank"]
        );
    }

    /// The 16S 27F primer, resolved at its degenerate positions.
    const FORWARD: &[u8] = b"AGAGTTTGATCCTGGCTCAG";
    /// The 16S 1492R primer, resolved at its degenerate position.
    const REVERSE: &[u8] = b"TACGGTTACCTTGTTACGACTT";
    /// A ligation adapter at the 5' end of every genomic read.
    const LIGATION: &[u8] = b"CCTGTACTTCGTTCAGTTACGTATTGC";

    /// Returns `seq` with about one base in sixty-four substituted.
    fn with_errors(seq: &[u8], seed: u64) -> Vec<u8> {
        let draws = bases(seed, seq.len() * 3);
        seq.iter()
            .zip(draws.chunks(3))
            .map(|(&base, draw)| match (draw, base) {
                (b"AAA", b'A') => b'C',
                (b"AAA", _) => b'A',
                _ => base,
            })
            .collect()
    }

    /// Reads of both strands of an amplicon library: 27F, a 1400-base
    /// insert and the reverse complement of 1492R, eroded at the 5' end by
    /// up to two bases. The first `amplicons` reads of 1500 are amplicons;
    /// the rest are genomic reads of 500 to 6000 bases that start with
    /// `LIGATION`.
    fn library(amplicons: u64) -> Vec<Vec<u8>> {
        (0..1500u64)
            .map(|i| {
                let read = if i < amplicons {
                    [
                        &FORWARD[(i % 3) as usize..],
                        &bases(7_000 + i, 1400)[..],
                        &crate::adapter::reverse_complement(REVERSE)[..],
                    ]
                    .concat()
                } else {
                    [LIGATION, &bases(9_000 + i, genomic_length(i))[..]].concat()
                };
                let read = with_errors(&read, 11_000 + i);
                if i % 2 == 1 && i < amplicons {
                    crate::adapter::reverse_complement(&read)
                } else {
                    read
                }
            })
            .collect()
    }

    /// Returns a genomic fragment length between 500 and 6000 bases for read
    /// `i`, spread as the fragments of a sheared genome are.
    fn genomic_length(i: u64) -> usize {
        500 + (i * 7919 % 5501) as usize
    }

    /// Resolves the set for `reads` under `infer` from the known entries of
    /// `base`, sampling every read.
    fn resolved(
        reads: Vec<Vec<u8>>,
        base: Vec<crate::adapter::Adapter>,
        infer: AdapterInfer,
    ) -> crate::adapter::AdapterConfig {
        let mut adapters = set(Vec::new());
        adapters.adapters = base;
        let cfg = Config {
            adapters: Some(adapters),
            adapter_infer: infer,
            adapter_sample: reads.len(),
            ..Config::default()
        };
        let records = reads.into_iter().map(anyhow::Ok);
        resolve(
            records,
            &cfg,
            |r: &Vec<u8>| Cow::Borrowed(&r[..]),
            |r| (r.len(), r.len()),
        )
        .unwrap()
        .unwrap()
        .adapters
        .unwrap()
    }

    /// Returns two amplicons joined by the 27F primer of the second, the
    /// first ending without its 1492R site, with the position of the
    /// junction primer, and a genomic read holding the 27F site alone.
    fn single_primer_reads() -> (Vec<u8>, usize, Vec<u8>) {
        let first = [FORWARD, &bases(21_000, 1400)[..]].concat();
        let second = [
            FORWARD,
            &bases(22_000, 1400)[..],
            &crate::adapter::reverse_complement(REVERSE)[..],
        ]
        .concat();
        let site = [
            LIGATION,
            &bases(23_000, 900)[..],
            FORWARD,
            &bases(24_000, 900)[..],
        ]
        .concat();
        let junction = first.len();
        let chimera = [first, second].concat();
        (chimera, junction, site)
    }

    /// Returns whether `read` splits into two segments under `ac` at the
    /// 27F primer that starts at `at`, give or take 30 bases.
    fn split_near(read: &[u8], ac: &crate::adapter::AdapterConfig, at: usize) -> bool {
        let segments = crate::adapter::adapter_segments(read, ac);
        segments.len() == 2
            && segments[0].1.abs_diff(at) <= 30
            && segments[1].0.abs_diff(at + FORWARD.len()) <= 30
    }

    /// Discovery on an amplicon library finds the marker primers at most
    /// read ends and splits a chimera at the single primer that joins two
    /// amplicons.
    #[test]
    fn discovery_on_an_amplicon_library_splits_at_a_single_primer() {
        let ac = resolved(library(1500), Vec::new(), AdapterInfer::Discover);
        assert!(ac.amplicon, "{:?}", ac.adapters);
        let (chimera, junction, _) = single_primer_reads();
        assert!(
            split_near(&chimera, &ac, junction),
            "{:?} {:?}",
            crate::adapter::adapter_segments(&chimera, &ac),
            ac.adapters
        );
    }

    /// A library in which marker primers open only a minority of the reads
    /// is no amplicon library: the discovered primers trim read ends, and a
    /// primer site inside a genomic read stays whole.
    #[test]
    fn discovery_on_a_genomic_library_keeps_a_single_primer_site() {
        let ac = resolved(library(300), Vec::new(), AdapterInfer::Discover);
        assert!(
            ac.adapters.iter().any(|a| {
                a.role == crate::adapter::Role::Primer
                    && crate::adapter::matches_marker_primer(&a.seq, 0.2)
            }),
            "{:?}",
            ac.adapters
        );
        assert!(!ac.amplicon);
        let (_, _, site) = single_primer_reads();
        let segments = crate::adapter::adapter_segments(&site, &ac);
        assert_eq!(segments.len(), 1, "{segments:?}");
    }

    /// A kit set that joins the amplicon kit to a genomic kit keeps the
    /// primer site of a genomic read whole, with and without discovery, and
    /// splits an amplicon chimera once the reads show an amplicon library.
    #[test]
    fn mixed_kit_set_splits_single_primers_only_in_an_amplicon_library() {
        use super::super::preset::{Kit, preset};
        let kits = preset(&[Kit::Mab114, Kit::Lsk114]);
        let (chimera, junction, site) = single_primer_reads();
        for infer in [AdapterInfer::Off, AdapterInfer::Discover] {
            let mut genomic = library(0);
            // Presence detection keeps an entry the reads carry.
            genomic.extend(std::iter::repeat_n(site.clone(), 20));
            let ac = resolved(genomic, kits.clone(), infer);
            assert!(!ac.amplicon);
            let segments = crate::adapter::adapter_segments(&site, &ac);
            assert_eq!(segments.len(), 1, "{infer:?} {segments:?}");

            let ac = resolved(library(1500), kits.clone(), infer);
            assert!(ac.amplicon, "{infer:?}");
            assert!(split_near(&chimera, &ac, junction), "{infer:?}");
        }
    }

    /// A primer pair outside the catalog, as a gene-specific amplicon uses.
    const CUSTOM_FORWARD: &[u8] = b"GGTCAACAAATCATAAAGATATTGG";
    /// The reverse primer of the custom pair.
    const CUSTOM_REVERSE: &[u8] = b"TAAACTTCAGGGTGACCAAAAAATCA";

    /// Returns variant `variant` of the 850-base target of the custom primer
    /// pair, as the few species of a sample carry the gene.
    fn custom_target(variant: u64) -> Vec<u8> {
        bases(32_000 + variant, 850)
    }

    /// Reads of both strands of a library amplified with the custom primer
    /// pair: the forward primer, one of four target variants and the reverse
    /// complement of the reverse primer, eroded at the 5' end by up to two
    /// bases and read with errors. `lead` precedes the forward primer of
    /// every third read.
    fn custom_reads(count: u64, seed: u64, lead: &[u8]) -> Vec<Vec<u8>> {
        (0..count)
            .map(|i| {
                let lead = if i % 3 == 0 { lead } else { &[][..] };
                let read = [
                    lead,
                    &CUSTOM_FORWARD[(i % 3) as usize..],
                    &custom_target(i % 12)[..],
                    &crate::adapter::reverse_complement(CUSTOM_REVERSE)[..],
                ]
                .concat();
                let read = with_errors(&read, seed + i);
                if i % 2 == 1 {
                    crate::adapter::reverse_complement(&read)
                } else {
                    read
                }
            })
            .collect()
    }

    /// Returns two amplicons of the custom pair joined by the forward primer
    /// of the second, the first ending without its reverse primer site, with
    /// the position of the junction primer.
    fn custom_chimera() -> (Vec<u8>, usize) {
        let first = [CUSTOM_FORWARD, &with_errors(&custom_target(0), 41_000)[..]].concat();
        let second = [
            CUSTOM_FORWARD,
            &with_errors(&custom_target(1), 41_001)[..],
            &crate::adapter::reverse_complement(CUSTOM_REVERSE)[..],
        ]
        .concat();
        let junction = first.len();
        ([first, second].concat(), junction)
    }

    /// Returns whether `read` splits into two segments under `ac` at the
    /// custom forward primer that starts at `at`, give or take 30 bases.
    fn custom_split_near(read: &[u8], ac: &crate::adapter::AdapterConfig, at: usize) -> bool {
        let segments = crate::adapter::adapter_segments(read, ac);
        segments.len() == 2
            && segments[0].1.abs_diff(at) <= 30
            && segments[1].0.abs_diff(at + CUSTOM_FORWARD.len()) <= 30
    }

    /// A library of two primer pairs, the 16S pair and a pair outside the
    /// catalog, is one amplicon library: the primers open most reads
    /// together, though neither pair does alone, and a chimera of either kind
    /// splits at the single primer that joins its amplicons.
    #[test]
    fn two_primer_pairs_make_one_amplicon_library() {
        let mut reads = library(750);
        reads.truncate(750);
        reads.extend(custom_reads(750, 51_000, b""));
        let ac = resolved(reads.clone(), Vec::new(), AdapterInfer::Discover);
        assert!(ac.amplicon, "{:?}", ac.adapters);
        let (chimera, junction, _) = single_primer_reads();
        assert!(split_near(&chimera, &ac, junction), "{:?}", ac.adapters);
        let (custom, at) = custom_chimera();
        assert!(custom_split_near(&custom, &ac, at), "{:?}", ac.adapters);
        use super::super::preset::{Kit, preset};
        let ac = resolved(
            reads,
            preset(&[Kit::Mab114, Kit::Lsk114]),
            AdapterInfer::Discover,
        );
        assert!(ac.amplicon);
    }

    /// A library amplified with a primer pair outside the catalog splits a
    /// chimera at a single interior copy of its discovered primer.
    #[test]
    fn a_custom_primer_splits_its_chimeras() {
        let ac = resolved(
            custom_reads(1500, 61_000, b""),
            Vec::new(),
            AdapterInfer::Discover,
        );
        let (custom, at) = custom_chimera();
        assert!(custom_split_near(&custom, &ac, at), "{:?}", ac.adapters);
    }

    /// Adapter remnants that trimming left before the primer at some read
    /// ends are trimmed with it, and the library still counts as an amplicon
    /// library.
    #[test]
    fn adapter_remnants_leave_the_amplicon_library_intact() {
        let remnant = &LIGATION[LIGATION.len() - 14..];
        let mut reads = library(900);
        reads.truncate(900);
        reads.extend(custom_reads(600, 71_000, remnant));
        let reads: Vec<Vec<u8>> = reads
            .into_iter()
            .enumerate()
            .map(|(i, read)| {
                if i % 4 == 0 {
                    [remnant, &read[..]].concat()
                } else {
                    read
                }
            })
            .collect();
        let ac = resolved(reads, Vec::new(), AdapterInfer::Discover);
        assert!(ac.amplicon, "{:?}", ac.adapters);
        let (chimera, junction, _) = single_primer_reads();
        let read = [remnant, &chimera[..]].concat();
        let segments = crate::adapter::adapter_segments(&read, &ac);
        assert_eq!(segments.len(), 2, "{segments:?} {:?}", ac.adapters);
        assert!(
            segments[0].0 >= remnant.len() + FORWARD.len() - 3,
            "{segments:?}"
        );
        assert!(
            segments[1]
                .0
                .abs_diff(remnant.len() + junction + FORWARD.len())
                <= 30
        );
    }

    /// Genomic reads that hold the sites of the 16S and the custom primers
    /// inside them keep them whole, with a kit set that joins the amplicon
    /// kit to a genomic kit and with discovery: the primers open no reads.
    #[test]
    fn genomic_reads_keep_every_primer_site_whole() {
        use super::super::preset::{Kit, preset};
        let site = [
            LIGATION,
            &bases(81_000, 1200)[..],
            FORWARD,
            &bases(81_001, 1300)[..],
            CUSTOM_FORWARD,
            &bases(81_002, 1100)[..],
        ]
        .concat();
        let mut genomic = library(0);
        genomic.extend(std::iter::repeat_n(site.clone(), 20));
        let known = preset(&[Kit::Mab114, Kit::Lsk114]);
        for infer in [AdapterInfer::Off, AdapterInfer::Discover] {
            let ac = resolved(genomic.clone(), known.clone(), infer);
            assert!(!ac.amplicon);
            let segments = crate::adapter::adapter_segments(&site, &ac);
            assert_eq!(segments, vec![(LIGATION.len(), site.len())], "{infer:?}");
        }
    }
}
