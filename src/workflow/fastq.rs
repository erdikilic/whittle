use std::io::Write;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use super::{Counters, FASTQ_BATCH, Rendered, Stats, process_read_segments, run_parallel};
use crate::config::Config;
use crate::io::fastq::write_segment;
use crate::record::ReadRecord;
use crate::trim;

/// Single-threaded FASTQ workflow: trim -> filter each produced segment -> write
/// survivors.
pub fn run_fastq_seq<W: Write>(
    records: impl Iterator<Item = anyhow::Result<ReadRecord>>,
    writer: &mut W,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    for rec in records {
        let rec = rec?;
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .input_bases
            .fetch_add(rec.seq.len() as u64, Ordering::Relaxed);
        let _read = super::read_span(&rec.name);
        let _read = _read.enter();
        let produced = trim::apply(&rec.seq, &rec.qual, &cfg.trim, cfg.adapters.as_ref());
        process_read_segments(
            &produced,
            &rec.seq,
            &rec.qual,
            &cfg.filter,
            counters,
            |idx, total, s, e| {
                write_segment(
                    writer,
                    &rec.name,
                    &rec.seq[s..e],
                    &rec.qual[s..e],
                    total,
                    idx,
                )?;
                Ok(())
            },
        )?;
    }
    Ok(counters.snapshot(0))
}

/// Trim one record, filter each produced segment through `process_read_segments`,
/// and render the survivors into an owned FASTQ buffer. Writing into an in-memory
/// `Vec<u8>` cannot fail, so the `.expect` below is an assertion, not error
/// handling: the parallel caller runs inside a plain `for_each` with no `Result`
/// seam (see `run_fastq`).
fn render_record(rec: &ReadRecord, cfg: &Config, counters: &Counters, buf: &mut Vec<u8>) {
    let _read = super::read_span(&rec.name);
    let _read = _read.enter();
    let produced = trim::apply(&rec.seq, &rec.qual, &cfg.trim, cfg.adapters.as_ref());
    process_read_segments(
        &produced,
        &rec.seq,
        &rec.qual,
        &cfg.filter,
        counters,
        |idx, total, s, e| {
            write_segment(
                &mut *buf,
                &rec.name,
                &rec.seq[s..e],
                &rec.qual[s..e],
                total,
                idx,
            )?;
            Ok(())
        },
    )
    .expect("writing FASTQ segments into an in-memory Vec<u8> cannot fail");
}

/// Threads-aware FASTQ workflow entry point. Sequential when `cfg.threads <= 1`;
/// otherwise records render on a rayon pool and drain through `run_parallel`,
/// in input order under `cfg.ordered` and in completion order otherwise.
pub fn run_fastq<W, I>(
    records: I,
    writer: &mut W,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats>
where
    W: Write + Send,
    I: Iterator<Item = anyhow::Result<ReadRecord>> + Send,
{
    if cfg.threads <= 1 {
        return run_fastq_seq(records, writer, cfg, counters);
    }
    run_parallel(
        records,
        FASTQ_BATCH,
        |rec: &ReadRecord| rec.seq.len(),
        cfg,
        writer,
        |rec, cfg| {
            let mut buf = Vec::with_capacity(rec.seq.len().saturating_mul(2) + rec.name.len() + 6);
            render_record(&rec, cfg, counters, &mut buf);
            Ok(Rendered {
                items: if buf.is_empty() { Vec::new() } else { vec![buf] },
                malformed_tags: false,
            })
        },
        |writer, buf: &Vec<u8>| writer.write_all(buf),
        counters,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::{self, FilterConfig};
    use crate::qual::QualMode;
    use crate::record::ReadRecord;
    use crate::trim::{QualityOp, TrimPlan};

    fn rec(name: &str, seq: &[u8], phred: Vec<u8>) -> ReadRecord {
        ReadRecord {
            name: name.as_bytes().to_vec(),
            seq: seq.to_vec(),
            qual: phred,
        }
    }

    fn base_filter() -> FilterConfig {
        FilterConfig {
            min_length: 1,
            max_length: usize::MAX,
            min_qual: 0.0,
            max_qual: 1000.0,
            min_gc: None,
            max_gc: None,
            qual_mode: QualMode::Mean,
        }
    }

    #[test]
    fn shared_counters_reflect_totals() {
        use std::sync::Arc;

        use crate::workflow::Counters;
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: TrimPlan {
                head: 1,
                tail: 1,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        let recs = vec![Ok(rec("r1", b"ACGT", vec![40, 40, 40, 40]))];
        let mut out = Vec::new();
        let counters = Arc::new(Counters::default());
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg, &counters).unwrap();
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
        assert_eq!(
            counters
                .input_reads
                .load(std::sync::atomic::Ordering::Relaxed),
            1
        );
        assert_eq!(
            counters
                .output_reads
                .load(std::sync::atomic::Ordering::Relaxed),
            1
        );
    }

    #[test]
    fn fixed_crop_writes_one_segment() {
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: TrimPlan {
                head: 1,
                tail: 1,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        let recs = vec![Ok(rec("r1", b"ACGT", vec![40, 40, 40, 40]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(
            recs.into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(out, b"@r1\nCG\n+\nII\n");
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
    }

    #[test]
    fn split_writes_suffixed_segments() {
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: Some(QualityOp::Split {
                    cutoff: 10,
                    window: 1,
                }),
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        // good(3) bad(1) good(3): I I I # I I I  -> two segments (0,3),(4,7)
        let phred: Vec<u8> = b"III#III".iter().map(|&b| b - 33).collect();
        let recs = vec![Ok(rec("r1", b"AAATAAA", phred))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(
            recs.into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(
            out,
            b"@r1_segment_1\nAAA\n+\nIII\n@r1_segment_2\nAAA\n+\nIII\n"
        );
        assert_eq!((stats.input_reads, stats.output_reads), (1, 2));
    }

    #[test]
    fn filtered_read_produces_no_output() {
        let mut f = base_filter();
        f.min_length = 10;
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: f,
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        let recs = vec![Ok(rec("short", b"ACGT", vec![40; 4]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(
            recs.into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert!(out.is_empty());
        assert_eq!((stats.input_reads, stats.output_reads), (1, 0));
    }

    #[test]
    fn too_short_segment_bumps_segments_dropped_short_counter() {
        // One produced segment rejected by length counts as all-filtered, not
        // trimmed-to-nothing.
        let mut f = base_filter();
        f.min_length = 10;
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: f,
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        let recs = vec![Ok(rec("short", b"ACGT", vec![40; 4]))];
        let mut out = Vec::new();
        let counters = Arc::new(Counters::default());
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg, &counters).unwrap();
        assert_eq!(stats.segments_dropped_short, 1);
        assert_eq!(stats.reads_all_filtered, 1);
        assert_eq!(stats.reads_trimmed_to_nothing, 0);
        assert_eq!(
            counters
                .segments_dropped_short
                .load(std::sync::atomic::Ordering::Relaxed),
            1
        );
    }

    #[test]
    fn trimmed_to_nothing_bumps_reads_trimmed_to_nothing_counter() {
        // A crop that removes the complete read produces no segments.
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: TrimPlan {
                head: 10,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        let recs = vec![Ok(rec("r1", b"ACGT", vec![40; 4]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(
            recs.into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert!(out.is_empty());
        assert_eq!(stats.reads_trimmed_to_nothing, 1);
        assert_eq!(stats.reads_all_filtered, 0);
        assert_eq!(stats.segments_dropped_short, 0);
    }

    /// Filtering uses the cropped sequence rather than the original read.
    #[test]
    fn quality_below_raw_mean_but_above_trimmed_insert_survives() {
        let mut f = base_filter();
        f.qual_mode = QualMode::Arithmetic;
        f.min_qual = 30.0;
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: f,
            trim: TrimPlan {
                head: 4,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        // Original mean: 24.8. Cropping four Q2 bases leaves six Q40 bases.
        let mut phred = vec![2u8; 4];
        phred.extend(std::iter::repeat_n(40u8, 6));
        assert!(
            filter::check(b"AAAAAAAAAA", &phred, &cfg.filter).is_some(),
            "the complete input read must fail the quality filter"
        );
        let recs = vec![Ok(rec("r1", b"AAAAAAAAAA", phred))];
        let mut out = Vec::new();
        let counters = Arc::new(Counters::default());
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg, &counters).unwrap();
        assert_eq!(out, b"@r1\nAAAAAA\n+\nIIIIII\n");
        assert_eq!(stats.output_reads, 1);
        assert_eq!(
            counters
                .reads_with_output
                .load(std::sync::atomic::Ordering::Relaxed),
            1
        );
        assert_eq!(stats.reads_trimmed_to_nothing, 0);
        assert_eq!(stats.reads_all_filtered, 0);
        assert_eq!(stats.segments_dropped_low_qual, 0);
    }

    /// An interior adapter splits a read
    /// into a long insert (survives length filtering) and a short insert
    /// (rejected `TooShort`). The survivor keeps its PRODUCED index: it
    /// is named `_segment_1` (not renamed to look unsplit), even
    /// though its sibling `_segment_2` never made it to output. A lone
    /// suffix with a gap correctly signals "this read was split".
    #[test]
    fn split_produces_long_survivor_and_short_segment_drop() {
        use crate::adapter::{Adapter, AdapterConfig, End};

        let adapter = b"GGGGTTTTGGGGTTTT"; // 16 bp, no A/C so it can't match the flanks
        let mut seq = vec![b'A'; 24]; // long flank -> survives length filter
        seq.extend_from_slice(adapter);
        seq.extend_from_slice(&[b'C'; 4]); // short flank -> TooShort
        let phred = vec![40u8; seq.len()];

        let mut f = base_filter();
        f.min_length = 5;
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: f,
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: Some(AdapterConfig {
                adapters: vec![Adapter {
                    name: "mid".into(),
                    seq: adapter.to_vec(),
                    end: End::Both,
                }],
                error_rate: 0.1,
                // end_size=1: both flanks (distance 24 and 4 from the match)
                // sit outside end_size, so the adapter classifies as interior
                // and the read splits rather than being terminal-trimmed.
                end_size: 1,
                split: true,
                candidate_index: std::sync::OnceLock::new(),
            }),
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        let recs = vec![Ok(rec("r1", &seq, phred))];
        let mut out = Vec::new();
        let counters = Arc::new(Counters::default());
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg, &counters).unwrap();

        assert_eq!(stats.output_reads, 1, "only the long flank survives");
        assert_eq!(
            stats.segments_dropped_short, 1,
            "the short flank is dropped"
        );
        assert_eq!(
            counters
                .reads_with_output
                .load(std::sync::atomic::Ordering::Relaxed),
            1
        );
        assert_eq!(stats.reads_trimmed_to_nothing, 0);
        assert_eq!(stats.reads_all_filtered, 0);
        let s = String::from_utf8(out).unwrap();
        assert!(
            s.starts_with("@r1_segment_1\n"),
            "survivor keeps its PRODUCED index (1 of 2), not renamed unsplit: {s:?}"
        );
        assert!(
            !s.contains("_segment_2"),
            "the dropped short segment must not appear in output: {s:?}"
        );
    }

    /// An empty input read produces no trim
    /// intervals at all, so it bumps the read-level `reads_trimmed_to_nothing`
    /// counter with NO segment-level drop (the per-segment filter loop never
    /// runs, since there is nothing to iterate).
    #[test]
    fn empty_read_bumps_reads_trimmed_to_nothing_with_no_segment_drop() {
        let cfg = Config {
            io: crate::config::IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        let recs = vec![Ok(rec("empty", b"", vec![]))];
        let mut out = Vec::new();
        let counters = Arc::new(Counters::default());
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg, &counters).unwrap();
        assert!(out.is_empty());
        assert_eq!(stats.input_reads, 1);
        assert_eq!(stats.output_reads, 0);
        assert_eq!(stats.reads_trimmed_to_nothing, 1);
        assert_eq!(stats.reads_all_filtered, 0);
        assert_eq!(stats.segments_dropped_short, 0);
        assert_eq!(stats.segments_dropped_long, 0);
        assert_eq!(stats.segments_dropped_low_qual, 0);
        assert_eq!(stats.segments_dropped_high_qual, 0);
        assert_eq!(stats.segments_dropped_gc, 0);
    }

    #[test]
    fn parallel_matches_sequential_as_multiset() {
        use crate::config::IoConfig;
        let mk = |threads| Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: Some(QualityOp::TrimQual(20)),
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        // Owned records (ReadRecord: Clone); wrap in Ok at iteration time so each run
        // gets a fresh Send iterator. anyhow::Error is not Clone, so a
        // Vec<Result<..>> cannot be cloned; clone the Vec<ReadRecord> and re-wrap instead.
        let recs: Vec<ReadRecord> = (0..500)
            .map(|i| rec(&format!("r{i}"), b"ACGTACGTAC", vec![40; 10]))
            .collect();

        let mut seq_out = Vec::new();
        run_fastq(
            recs.clone().into_iter().map(anyhow::Ok),
            &mut seq_out,
            &mk(1),
            &Arc::new(Counters::default()),
        )
        .unwrap();

        let mut par_out = Vec::new();
        run_fastq(
            recs.into_iter().map(anyhow::Ok),
            &mut par_out,
            &mk(4),
            &Arc::new(Counters::default()),
        )
        .unwrap();

        let sort_records = |bytes: &[u8]| {
            let mut v: Vec<Vec<u8>> = bytes
                .split(|&b| b == b'@')
                .filter(|s| !s.is_empty())
                .map(|s| s.to_vec())
                .collect();
            v.sort();
            v
        };
        assert_eq!(sort_records(&seq_out), sort_records(&par_out));
    }

    #[test]
    fn parallel_surfaces_write_error_without_deadlock() {
        use std::io::{self, Write};

        use crate::config::IoConfig;

        struct FailAfter {
            limit: usize,
            written: usize,
        }
        impl Write for FailAfter {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                if self.written >= self.limit {
                    return Err(io::Error::new(io::ErrorKind::BrokenPipe, "boom"));
                }
                self.written += buf.len();
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let cfg = Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 4,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        // Exceed the bounded channel capacity before the writer fails.
        let recs: Vec<ReadRecord> = (0..2000)
            .map(|i| rec(&format!("r{i}"), b"ACGTACGTAC", vec![40; 10]))
            .collect();
        let mut w = FailAfter {
            limit: 100,
            written: 0,
        };
        let res = run_fastq(
            recs.into_iter().map(anyhow::Ok),
            &mut w,
            &cfg,
            &Arc::new(Counters::default()),
        );
        assert!(
            res.is_err(),
            "write error must surface as Err, and must not hang"
        );
    }

    #[test]
    fn parallel_surfaces_parse_error_instead_of_dropping_it() {
        use crate::config::IoConfig;

        let cfg = Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: base_filter(),
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 4,
            fastq_tags: crate::config::FastqTags::All,
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
        };
        let good: Vec<anyhow::Result<ReadRecord>> = (0..5)
            .map(|i| anyhow::Ok(rec(&format!("r{i}"), b"ACGTACGTAC", vec![40; 10])))
            .collect();
        let recs = good
            .into_iter()
            .chain(std::iter::once(Err(anyhow::anyhow!("bad record"))));

        let mut out = Vec::new();
        let res = run_fastq(recs, &mut out, &cfg, &Arc::new(Counters::default()));
        assert!(
            res.is_err(),
            "a malformed record must not be silently dropped on the parallel path"
        );
    }
}
