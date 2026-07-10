use std::io::Write;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use rayon::prelude::*;

use super::{Counters, Stats, process_read_segments};
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

/// Trim one record, then filter each produced segment via the shared
/// `process_read_segments` helper, rendering the survivors into an owned
/// FASTQ byte buffer. Writing into an in-memory `Vec<u8>` cannot fail, so the
/// `.expect` below is an assertion, not real error handling — the parallel
/// caller runs inside a plain `for_each` with no `Result` propagation seam
/// (see `run_fastq`), matching the pre-refactor `.unwrap()` on the same
/// write.
fn render_record(rec: &ReadRecord, cfg: &Config, counters: &Counters) -> Vec<u8> {
    let produced = trim::apply(&rec.seq, &rec.qual, &cfg.trim, cfg.adapters.as_ref());
    let mut buf = Vec::new();
    process_read_segments(
        &produced,
        &rec.seq,
        &rec.qual,
        &cfg.filter,
        counters,
        |idx, total, s, e| {
            write_segment(
                &mut buf,
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
    buf
}

/// Threads-aware FASTQ workflow entry point. Sequential (and output-order
/// deterministic) when `cfg.threads <= 1`; otherwise renders each record on a
/// rayon work pool and drains rendered buffers through a bounded channel on a
/// dedicated writer task. Output order is unordered for `threads > 1` since
/// buffers are written in arrival order, not input order.
///
/// Record-parse errors (`Err` items from `records`): both paths surface them
/// as an `Err` from this function. The sequential path aborts immediately via
/// `?`; the parallel path captures the first parse error in a shared slot and
/// keeps draining so producer threads never block on the bounded channel,
/// surfacing the error once the scope joins (matching the sequential path's
/// fail behavior, just at end-of-run instead of mid-stream).
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

    let render_workers = if cfg.render_workers >= 1 {
        cfg.render_workers
    } else {
        cfg.threads.max(1)
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(render_workers)
        .build()?;
    let (tx, rx) = crossbeam_channel::bounded::<Vec<u8>>(render_workers * 4);

    let write_err: std::sync::Mutex<Option<std::io::Error>> = std::sync::Mutex::new(None);
    let parse_err: std::sync::Mutex<Option<anyhow::Error>> = std::sync::Mutex::new(None);

    // Writer on a plain scoped OS thread; render on the budget-sized LOCAL rayon
    // pool via `pool.install` so the nested `par_bridge` uses THAT pool (not
    // rayon's global num_cpus pool) — this is what makes `-t` bound the render
    // threads. The writer keeps draining on a write error (never `break`) so
    // bounded-channel producers can't deadlock.
    std::thread::scope(|s| {
        let write_err_ref = &write_err;
        s.spawn(move || {
            let mut errored = false;
            for buf in rx.iter() {
                if errored {
                    continue; // keep draining so bounded-channel producers never block
                }
                if let Err(e) = writer.write_all(&buf) {
                    *write_err_ref.lock().unwrap() = Some(e);
                    errored = true;
                }
            }
        });

        pool.install(|| {
            records.par_bridge().for_each(|rec| {
                let rec = match rec {
                    Ok(r) => r,
                    Err(e) => {
                        let mut g = parse_err.lock().unwrap();
                        if g.is_none() {
                            *g = Some(e);
                        }
                        return;
                    },
                };
                counters.input_reads.fetch_add(1, Ordering::Relaxed);
                counters
                    .input_bases
                    .fetch_add(rec.seq.len() as u64, Ordering::Relaxed);
                let buf = render_record(&rec, cfg, counters);
                if !buf.is_empty() {
                    let _ = tx.send(buf);
                }
            });
        });
        drop(tx);
    });

    if let Some(e) = parse_err.lock().unwrap().take() {
        return Err(e);
    }
    if let Some(e) = write_err.lock().unwrap().take() {
        return Err(e.into());
    }
    // FASTQ input carries no BAM per-base tags.
    Ok(counters.snapshot(0))
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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
        // Reorder regression: `filter::check` now runs POST-trim, per produced
        // segment (previously it ran pre-trim on the whole raw read). With no
        // adapters/quality-op configured, `trim::apply` still returns the whole
        // (untrimmed) read as its one produced segment, so the length filter
        // rejects that segment rather than the raw read up front — same
        // observable drop, but now counted at the segment level. Since exactly
        // one segment was produced (not zero) and it was filtered, this is
        // `reads_all_filtered`, not `reads_trimmed_to_nothing`.
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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
        // Reorder regression: `trim::apply` producing zero segments (a head-crop
        // of 10 exceeds the 4-base read length, so no window survives) never
        // enters the per-segment loop at all, bumping the read-level
        // `reads_trimmed_to_nothing` counter (distinct from `reads_all_filtered`,
        // which requires at least one produced segment that was then filtered).
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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

    /// A read's RAW mean quality is
    /// dragged below `-q` only by a low-quality head flank; once that flank is
    /// cropped away (crop runs before the filter in the new order), the
    /// trimmed insert's own mean passes. Pre-fix (filter-before-trim), this
    /// read would have been rejected on the raw mean; post-fix it SURVIVES.
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        };
        // 4 low-quality bases (phred 2) then 6 high-quality bases (phred 40).
        // Raw arithmetic mean = (2*4 + 40*6) / 10 = 24.8 < 30 (would fail the
        // OLD pre-trim whole-read filter). After a head-crop of 4, the
        // trimmed insert's mean is 40 >= 30 -> passes.
        let mut phred = vec![2u8; 4];
        phred.extend(std::iter::repeat_n(40u8, 6));
        assert!(
            filter::check(b"AAAAAAAAAA", &phred, &cfg.filter).is_some(),
            "sanity: the RAW whole read must fail the filter on its own"
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
    /// though its sibling `_segment_2` never made it to output — a lone
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
            }),
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        };
        // Owned records (ReadRecord: Clone); wrap in Ok at iteration time so each run
        // gets a fresh Send iterator. anyhow::Error is not Clone, so we can't clone a
        // Vec<Result<..>> — clone the Vec<ReadRecord> and re-wrap instead.
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
        };
        // Far more records than the bounded channel capacity (threads*4), so a
        // pre-fix build would deadlock instead of returning.
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
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
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
