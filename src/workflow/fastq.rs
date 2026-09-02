//! FASTQ workflows: sequential and parallel trim, filter and write drivers.

use std::io::Write;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use super::{Counters, FASTQ_BATCH, Rendered, Stats, process_read_segments, run_parallel};
use crate::config::Config;
use crate::io::fastq::write_segment;
use crate::record::ReadRecord;
use crate::trim;

/// Runs the single-threaded FASTQ workflow: trims, filters each produced segment
/// and writes the survivors.
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

/// Trims one record, filters each produced segment through
/// `process_read_segments`, and renders the survivors into `buf`. Writing into
/// an in-memory `Vec<u8>` cannot fail, so the `expect` is an assertion rather
/// than error handling.
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
    .expect("Writing FASTQ segments into an in-memory Vec<u8> cannot fail");
}

/// Runs the FASTQ workflow: sequential when `cfg.threads <= 1`; otherwise
/// records render on a rayon pool and drain through `run_parallel`, in input
/// order under `cfg.ordered` and in completion order otherwise.
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
                items: if buf.is_empty() {
                    Vec::new()
                } else {
                    vec![buf]
                },
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
    use crate::filter;
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

    /// A pass-through config (no trimming, no filtering, no adapters) with the
    /// given thread count; each test overrides the fields it exercises.
    fn test_cfg(threads: usize) -> Config {
        let null = std::path::Path::new("/dev/null");
        crate::cli::config_for_test_threads(null, null, 0, 0, threads)
    }

    #[test]
    fn shared_counters_reflect_totals() {
        use std::sync::Arc;

        use crate::workflow::Counters;
        let mut cfg = test_cfg(1);
        cfg.trim = TrimPlan {
            head: 1,
            tail: 1,
            quality: None,
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
        let mut cfg = test_cfg(1);
        cfg.trim = TrimPlan {
            head: 1,
            tail: 1,
            quality: None,
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
        let mut cfg = test_cfg(1);
        cfg.trim.quality = Some(QualityOp::Split {
            cutoff: 10,
            window: 1,
        });
        // Three good, one bad, three good (`III#III`) gives two segments, (0,3)
        // and (4,7).
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
        let mut cfg = test_cfg(1);
        cfg.filter.min_length = 10;
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

    /// One produced segment rejected by length counts as all-filtered, not
    /// trimmed-to-nothing.
    #[test]
    fn too_short_segment_bumps_segments_dropped_short_counter() {
        let mut cfg = test_cfg(1);
        cfg.filter.min_length = 10;
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

    /// A crop that removes the complete read produces no segments.
    #[test]
    fn trimmed_to_nothing_bumps_reads_trimmed_to_nothing_counter() {
        let mut cfg = test_cfg(1);
        cfg.trim.head = 10;
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
        let mut cfg = test_cfg(1);
        cfg.filter.qual_mode = QualMode::Arithmetic;
        cfg.filter.min_qual = 30.0;
        cfg.trim.head = 4;
        // Original mean: 24.8. Cropping four Q2 bases leaves six Q40 bases.
        let mut phred = vec![2u8; 4];
        phred.extend(std::iter::repeat_n(40u8, 6));
        assert!(
            filter::check(b"AAAAAAAAAA", &phred, &cfg.filter).is_some(),
            "The complete input read must fail the quality filter"
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

    /// An interior adapter splits a read into a long insert (survives length
    /// filtering) and a short insert (rejected `TooShort`). The survivor keeps
    /// its produced index: it is named `_segment_1`, not renamed to look
    /// unsplit, although its sibling `_segment_2` is not written. A lone suffix
    /// with a gap signals that the read was split.
    #[test]
    fn split_produces_long_survivor_and_short_segment_drop() {
        use crate::adapter::{Adapter, AdapterConfig, End};

        let adapter = b"GGGGTTTTGGGGTTTT"; // 16 bp, no A/C, so it cannot match the flanks
        let mut seq = vec![b'A'; 24]; // long flank, survives the length filter
        seq.extend_from_slice(adapter);
        seq.extend_from_slice(&[b'C'; 4]); // short flank, TooShort
        let phred = vec![40u8; seq.len()];

        let mut cfg = test_cfg(1);
        cfg.filter.min_length = 5;
        cfg.adapters = Some(AdapterConfig {
            adapters: vec![Adapter {
                name: "mid".into(),
                seq: adapter.to_vec(),
                end: End::Both,
            }],
            error_rate: 0.1,
            // With `end_size` 1, both flanks (24 and 4 bases from the match)
            // sit outside the end zone, so the adapter is interior and the
            // read splits rather than being terminal-trimmed.
            end_size: 1,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        });
        let recs = vec![Ok(rec("r1", &seq, phred))];
        let mut out = Vec::new();
        let counters = Arc::new(Counters::default());
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg, &counters).unwrap();

        assert_eq!(stats.output_reads, 1, "Only the long flank survives");
        assert_eq!(
            stats.segments_dropped_short, 1,
            "The short flank is dropped"
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
            "Survivor keeps its produced index (1 of 2): {s:?}"
        );
        assert!(
            !s.contains("_segment_2"),
            "The dropped short segment must not appear in output: {s:?}"
        );
    }

    /// An empty input read produces no trim intervals, so it bumps the
    /// read-level `reads_trimmed_to_nothing` counter with no segment-level
    /// drop: the per-segment filter loop never runs.
    #[test]
    fn empty_read_bumps_reads_trimmed_to_nothing_with_no_segment_drop() {
        let cfg = test_cfg(1);
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
        let mk = |threads| {
            let mut cfg = test_cfg(threads);
            cfg.trim.quality = Some(QualityOp::TrimQual(20));
            cfg
        };
        // Owned records (`ReadRecord: Clone`), wrapped in `Ok` at iteration time
        // so each run gets a fresh `Send` iterator. `anyhow::Error` is not
        // `Clone`, so a `Vec<Result<..>>` cannot be cloned; the `Vec<ReadRecord>`
        // is cloned and re-wrapped instead.
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

        let cfg = test_cfg(4);
        // Enough records to exceed the bounded channel capacity before the
        // writer fails.
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
            "Write error must surface as Err and must not hang"
        );
    }

    #[test]
    fn parallel_surfaces_parse_error_instead_of_dropping_it() {
        let cfg = test_cfg(4);
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
            "A malformed record must not be dropped on the parallel path"
        );
    }
}
