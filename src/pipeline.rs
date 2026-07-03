use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

use crate::config::Config;
use crate::io::fastq::write_segment;
use crate::record::ReadRecord;
use crate::{filter, trim};

#[derive(Debug, Default, Clone, Copy)]
pub struct Stats {
    pub input_reads: u64,
    pub output_reads: u64,
}

/// Single-threaded FASTQ pipeline: filter -> trim -> write each surviving segment.
pub fn run_fastq_seq<W: Write>(
    records: impl Iterator<Item = anyhow::Result<ReadRecord>>,
    writer: &mut W,
    cfg: &Config,
) -> anyhow::Result<Stats> {
    let mut stats = Stats::default();
    for rec in records {
        let rec = rec?;
        stats.input_reads += 1;
        if !filter::passes(&rec.seq, &rec.qual, &cfg.filter) {
            continue;
        }
        let intervals = trim::apply(rec.seq.len(), &rec.qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            write_segment(writer, &rec.name, &rec.seq[s..e], &rec.qual[s..e], total, idx)?;
            stats.output_reads += 1;
        }
    }
    Ok(stats)
}

/// Format the surviving segments of one record into an owned FASTQ byte buffer.
/// Returns the number of segments written alongside the buffer.
fn render_record(rec: &ReadRecord, cfg: &Config) -> (u64, Vec<u8>) {
    if !filter::passes(&rec.seq, &rec.qual, &cfg.filter) {
        return (0, Vec::new());
    }
    let intervals = trim::apply(rec.seq.len(), &rec.qual, &cfg.trim, cfg.filter.min_length);
    let total = intervals.len();
    let mut buf = Vec::new();
    let mut out = 0u64;
    for (idx, (s, e)) in intervals.into_iter().enumerate() {
        write_segment(&mut buf, &rec.name, &rec.seq[s..e], &rec.qual[s..e], total, idx).unwrap();
        out += 1;
    }
    (out, buf)
}

/// Threads-aware FASTQ pipeline entry point. Sequential (and output-order
/// deterministic) when `cfg.threads <= 1`; otherwise renders each record on a
/// rayon work pool and drains rendered buffers through a bounded channel on a
/// dedicated writer task. Output order is unordered for `threads > 1` since
/// buffers are written in arrival order, not input order.
///
/// Record-parse errors (`Err` items from `records`): the sequential path
/// surfaces them via `?` and aborts the whole run. The parallel path currently
/// drops the offending item and continues — for Plan 1, malformed FASTQ over
/// multiple threads only aborts the affected read, not the whole run. Strict
/// propagation would need a shared error slot mirroring `write_err` below.
pub fn run_fastq<W, I>(records: I, writer: &mut W, cfg: &Config) -> anyhow::Result<Stats>
where
    W: Write + Send,
    I: Iterator<Item = anyhow::Result<ReadRecord>> + Send,
{
    if cfg.threads <= 1 {
        return run_fastq_seq(records, writer, cfg);
    }

    let pool = rayon::ThreadPoolBuilder::new().num_threads(cfg.threads).build()?;
    let input_reads = AtomicU64::new(0);
    let output_reads = AtomicU64::new(0);
    let (tx, rx) = crossbeam_channel::bounded::<Vec<u8>>(cfg.threads * 4);

    let write_err: std::sync::Mutex<Option<std::io::Error>> = std::sync::Mutex::new(None);

    pool.in_place_scope(|scope| {
        // Writer task drains rendered buffers in arrival order.
        scope.spawn(|_| {
            for buf in rx.iter() {
                if let Err(e) = writer.write_all(&buf) {
                    *write_err.lock().unwrap() = Some(e);
                    break;
                }
            }
        });

        records.par_bridge().for_each(|rec| {
            let rec = match rec {
                Ok(r) => r,
                Err(_) => return,
            };
            input_reads.fetch_add(1, Ordering::Relaxed);
            let (out, buf) = render_record(&rec, cfg);
            if out > 0 {
                output_reads.fetch_add(out, Ordering::Relaxed);
                let _ = tx.send(buf);
            }
        });
        drop(tx);
    });

    if let Some(e) = write_err.lock().unwrap().take() {
        return Err(e.into());
    }
    Ok(Stats {
        input_reads: input_reads.load(Ordering::Relaxed),
        output_reads: output_reads.load(Ordering::Relaxed),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::filter::FilterConfig;
    use crate::qual::QualMode;
    use crate::record::ReadRecord;
    use crate::trim::{QualityOp, TrimPlan};

    fn rec(name: &str, seq: &[u8], phred: Vec<u8>) -> ReadRecord {
        ReadRecord { name: name.as_bytes().to_vec(), seq: seq.to_vec(), qual: phred }
    }

    fn base_filter() -> FilterConfig {
        FilterConfig {
            min_length: 1, max_length: usize::MAX, min_qual: 0.0, max_qual: 1000.0,
            min_gc: None, max_gc: None, qual_mode: QualMode::Mean,
        }
    }

    #[test]
    fn fixed_crop_writes_one_segment() {
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 1, tail: 1, quality: None },
            threads: 1,
        };
        let recs = vec![Ok(rec("r1", b"ACGT", vec![40, 40, 40, 40]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@r1\nCG\n+\nII\n");
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
    }

    #[test]
    fn split_writes_suffixed_segments() {
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 0, tail: 0, quality: Some(QualityOp::Split { cutoff: 10, window: 1 }) },
            threads: 1,
        };
        // good(3) bad(1) good(3): I I I # I I I  -> two segments (0,3),(4,7)
        let phred: Vec<u8> = b"III#III".iter().map(|&b| b - 33).collect();
        let recs = vec![Ok(rec("r1", b"AAATAAA", phred))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@r1_segment_1\nAAA\n+\nIII\n@r1_segment_2\nAAA\n+\nIII\n");
        assert_eq!((stats.input_reads, stats.output_reads), (1, 2));
    }

    #[test]
    fn filtered_read_produces_no_output() {
        let mut f = base_filter();
        f.min_length = 10;
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: f,
            trim: TrimPlan { head: 0, tail: 0, quality: None },
            threads: 1,
        };
        let recs = vec![Ok(rec("short", b"ACGT", vec![40; 4]))];
        let mut out = Vec::new();
        let stats = run_fastq_seq(recs.into_iter(), &mut out, &cfg).unwrap();
        assert!(out.is_empty());
        assert_eq!((stats.input_reads, stats.output_reads), (1, 0));
    }

    #[test]
    fn parallel_matches_sequential_as_multiset() {
        use crate::config::IoConfig;
        let mk = |threads| Config {
            io: IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 0, tail: 0, quality: Some(QualityOp::TrimQual(20)) },
            threads,
        };
        // Owned records (ReadRecord: Clone); wrap in Ok at iteration time so each run
        // gets a fresh Send iterator. anyhow::Error is not Clone, so we can't clone a
        // Vec<Result<..>> — clone the Vec<ReadRecord> and re-wrap instead.
        let recs: Vec<ReadRecord> = (0..500)
            .map(|i| rec(&format!("r{i}"), b"ACGTACGTAC", vec![40; 10]))
            .collect();

        let mut seq_out = Vec::new();
        run_fastq(recs.clone().into_iter().map(anyhow::Ok), &mut seq_out, &mk(1)).unwrap();

        let mut par_out = Vec::new();
        run_fastq(recs.into_iter().map(anyhow::Ok), &mut par_out, &mk(4)).unwrap();

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
}
