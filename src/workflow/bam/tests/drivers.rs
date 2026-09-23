//! Parallel BAM drivers against their sequential counterparts, and error
//! propagation.

use super::*;

#[test]
fn run_bam_parallel_matches_sequential_as_multiset() {
    use crate::trim::TrimPlan;

    let mk = |threads| Config {
        trim: TrimPlan {
            head: 2,
            tail: 2,
            quality: None,
        },
        threads,
        quiet: true,
        ..Config::default()
    };
    // 300 reads with mods so reconstruction runs on every one.
    let recs: Vec<RecordBuf> = (0..300)
        .map(|_| ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]))
        .collect();

    let header = sam::Header::default();
    let decode = |bytes: &[u8]| -> Vec<(Vec<u8>, Vec<u8>)> {
        // (seq, MM-bytes) pairs, sorted, as an order-independent fingerprint.
        let mut r = noodles_bam::io::Reader::new(bytes);
        let h = r.read_header().unwrap();
        let mut out = Vec::new();
        let mut buf = RecordBuf::default();
        while r.read_record_buf(&h, &mut buf).unwrap() != 0 {
            let seq = buf.sequence().as_ref().to_vec();
            let mm = match buf.data().get(&Tag::BASE_MODIFICATIONS) {
                Some(Value::String(s)) => s.to_vec(),
                _ => Vec::new(),
            };
            out.push((seq, mm));
        }
        out.sort();
        out
    };

    // t1: single-threaded BGZF sink, written to a temporary file.
    let dir = tempfile::tempdir().unwrap();
    let p1 = dir.path().join("t1.bam");
    let mut sink1 = crate::io::bam::writer(Some(&p1), &header, false, 6).unwrap();
    run_bam(
        &header,
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut sink1,
        &mk(1),
        &Arc::new(Counters::default()),
    )
    .unwrap();
    sink1.finish().unwrap();
    let b1 = std::fs::read(&p1).unwrap();

    // t8: multithreaded sink to a temporary file (the multithreaded writer
    // needs an owned `Write + Send`).
    let p8 = dir.path().join("t8.bam");
    let mut sink8 = crate::io::bam::writer(Some(&p8), &header, true, 6).unwrap();
    run_bam(
        &header,
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut sink8,
        &mk(8),
        &Arc::new(Counters::default()),
    )
    .unwrap();
    sink8.finish().unwrap();
    let b8 = std::fs::read(&p8).unwrap();

    assert_eq!(
        decode(&b1),
        decode(&b8),
        "The t1 and t8 runs must produce the same record set"
    );
}

/// Writer errors remain observable after the bounded channel reaches capacity.
#[test]
fn run_bam_parallel_surfaces_write_error_without_deadlock() {
    use std::io;

    struct FailAfter {
        limit: usize,
        written: usize,
    }

    let cfg = Config {
        threads: 4,
        quiet: true,
        ..Config::default()
    };
    let recs: Vec<anyhow::Result<bam::Record>> = (0..3000)
        .map(|_| anyhow::Ok(bam::Record::default()))
        .collect();

    let mut sink = FailAfter {
        limit: 100,
        written: 0,
    };
    let res = run_bam_parallel(
        recs.into_iter(),
        &cfg,
        &mut sink,
        |_raw, _rec, _cfg, out: &mut Vec<()>| {
            out.push(());
            Ok(())
        },
        Ok,
        |sink, batch: &Vec<()>| -> io::Result<()> {
            if sink.written >= sink.limit {
                return Err(io::Error::new(io::ErrorKind::BrokenPipe, "boom"));
            }
            sink.written += batch.len();
            Ok(())
        },
        &Arc::new(Counters::default()),
    );
    assert!(
        res.is_err(),
        "Write error must surface as Err and must not hang"
    );
}

/// Mirrors `workflow::fastq`'s `parallel_surfaces_parse_error_instead_of_dropping_it`,
/// driving `run_bam_parallel` directly so a malformed upstream record (an
/// `Err` item from the input iterator) is not silently swallowed.
#[test]
fn run_bam_parallel_surfaces_parse_error_instead_of_dropping_it() {
    use std::io;

    struct NullSink;

    let cfg = Config {
        threads: 4,
        quiet: true,
        ..Config::default()
    };
    let good: Vec<anyhow::Result<bam::Record>> =
        (0..5).map(|_| anyhow::Ok(bam::Record::default())).collect();
    let recs = good
        .into_iter()
        .chain(std::iter::once(Err(anyhow::anyhow!("bad record"))));

    let mut sink = NullSink;
    let res = run_bam_parallel(
        recs,
        &cfg,
        &mut sink,
        |_raw, _rec, _cfg, out: &mut Vec<()>| {
            out.push(());
            Ok(())
        },
        Ok,
        |_sink: &mut NullSink, _batch: &Vec<()>| -> io::Result<()> { Ok(()) },
        &Arc::new(Counters::default()),
    );
    assert!(
        res.is_err(),
        "A malformed record must not be dropped on the parallel path"
    );
}

#[test]
fn run_bam_to_fastq_parallel_matches_sequential_as_multiset() {
    use crate::trim::TrimPlan;

    let mk = |threads| Config {
        trim: TrimPlan {
            head: 2,
            tail: 2,
            quality: None,
        },
        threads,
        quiet: true,
        ..Config::default()
    };
    let recs: Vec<RecordBuf> = (0..300)
        .map(|_| ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]))
        .collect();

    let sorted_records = |bytes: &[u8]| {
        let s = String::from_utf8(bytes.to_vec()).unwrap();
        // Records are grouped as 4 consecutive lines rather than split on
        // `@`: a QUAL byte of Phred 31 (ASCII `@`) would corrupt an `@`
        // split. FASTQ records are exactly 4 lines each here, so the
        // re-chunking is lossless.
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(
            lines.len() % 4,
            0,
            "Expected whole 4-line FASTQ records, got {} lines",
            lines.len()
        );
        let mut v: Vec<String> = lines.chunks(4).map(|c| c.join("\n")).collect();
        v.sort();
        v
    };

    let mut a = Vec::new();
    run_bam_to_fastq(
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut a,
        &mk(1),
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let mut b = Vec::new();
    run_bam_to_fastq(
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut b,
        &mk(8),
        &Arc::new(Counters::default()),
    )
    .unwrap();

    assert_eq!(
        sorted_records(&a),
        sorted_records(&b),
        "The t1 and t8 FASTQ outputs must match as a multiset"
    );
}
