//! Writes the synthetic corpus the profile-guided release build trains on, into
//! the directory named by the first argument (`target/training-data` by
//! default).
//!
//! The corpus covers the paths `scripts/pgo-train.sh` exercises: plain and gzip
//! FASTQ carrying ONT header fields and preset adapters, and an unaligned BAM
//! carrying modification calls, per-base kinetics, a move table, and the ONT run
//! fields. A fixed seed makes the corpus byte-identical on every run, so a
//! rebuilt profile differs only when the code does.
//!
//! An example, not a binary, so it stays out of the shipped binary.

use std::io::Write as _;
use std::path::{Path, PathBuf};

use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};
use whittle::adapter::{Adapter, End, preset};

/// Seed of the corpus generator. Any fixed value serves; changing it reshuffles
/// the whole corpus.
const SEED: u64 = 0x1f2e_3d4c_5b6a_7988;

/// FASTQ reads written to both the plain and the gzip file.
const FASTQ_READS: usize = 3000;

/// BAM records written.
const BAM_RECORDS: usize = 2500;

/// Signal samples per move-table block.
const MOVE_STRIDE: i8 = 5;

/// Signal samples the basecaller trimmed from the front of every read, the `ts`
/// tag.
const SIGNAL_TRIM: i32 = 10;

/// Signal sampling rate, used to derive `du` from the sample count.
const SAMPLE_RATE: f32 = 5000.0;

/// Sequencing start as seconds since the Unix epoch, 2024-06-22T10:00:00Z. Read
/// timestamps are spread forward from it so the report's timeline has more than
/// one bucket.
const RUN_START: i64 = 1_719_050_400;

fn main() -> std::io::Result<()> {
    let dir = std::env::args_os()
        .nth(1)
        .map_or_else(|| PathBuf::from("target/training-data"), PathBuf::from);
    std::fs::create_dir_all(&dir)?;

    // Both FASTQ files hold the same reads, so the gzip path trims the same work
    // as the plain one.
    let fastq = build_fastq();
    let plain = dir.join("reads.fastq");
    std::fs::write(&plain, &fastq)?;

    let gz = dir.join("reads.fastq.gz");
    let mut encoder =
        flate2::write::GzEncoder::new(std::fs::File::create(&gz)?, flate2::Compression::new(4));
    encoder.write_all(&fastq)?;
    encoder.finish()?;

    let bam = dir.join("reads.bam");
    write_bam(&bam)?;

    for path in [&plain, &gz, &bam] {
        println!("{}", path.display());
    }
    Ok(())
}

/// SplitMix64, so the corpus is reproducible without a dependency on a random
/// number generator.
struct Rng(u64);

impl Rng {
    /// Creates a generator from a seed.
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// Returns the next value in the sequence.
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Returns a value in `[0, n)`.
    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

/// Returns the first preset sequence expected at `end`.
fn adapter_seq(adapters: &[Adapter], end: End) -> Vec<u8> {
    adapters
        .iter()
        .find(|a| a.end == end)
        .expect("the ONT preset carries a sequence for each end")
        .seq
        .clone()
}

/// Draws a FASTQ read length. Most reads are one to three kilobases, a tenth are
/// short enough for `--min-length` to drop, and a fifth run long, so the length
/// panels and the filters both see a spread.
fn fastq_length(rng: &mut Rng) -> usize {
    match rng.below(10) {
        0 => 120 + rng.below(280) as usize,
        1..=7 => 600 + rng.below(2400) as usize,
        _ => 3000 + rng.below(5000) as usize,
    }
}

/// Draws a BAM read length. Shorter than the FASTQ distribution, because every
/// record also carries per-base arrays.
fn bam_length(rng: &mut Rng) -> usize {
    match rng.below(10) {
        0 => 150 + rng.below(250) as usize,
        _ => 500 + rng.below(2500) as usize,
    }
}

/// Builds a sequence with a per-read GC content between 0.30 and 0.60, so the GC
/// panel and the `--min-gc`/`--max-gc` filters see a spread.
fn sequence(rng: &mut Rng, len: usize) -> Vec<u8> {
    let gc = 30 + rng.below(31);
    let mut seq = Vec::with_capacity(len);
    for _ in 0..len {
        let strong = rng.below(100) < gc;
        let second = rng.below(2) == 0;
        seq.push(match (strong, second) {
            (true, true) => b'G',
            (true, false) => b'C',
            (false, true) => b'A',
            (false, false) => b'T',
        });
    }
    seq
}

/// Builds raw Phred scores: a per-read mean with per-base jitter, degraded over
/// the first and last 60 bases so the end-trimming stages have something to cut.
fn quality(rng: &mut Rng, len: usize) -> Vec<u8> {
    let mean = 8 + rng.below(22) as i32;
    let mut qual = Vec::with_capacity(len);
    for i in 0..len {
        let edge = i.min(len - 1 - i);
        let penalty = if edge < 60 { (60 - edge) as i32 / 5 } else { 0 };
        let jitter = rng.below(7) as i32 - 3;
        qual.push((mean - penalty + jitter).clamp(2, 40) as u8);
    }
    qual
}

/// Returns the mean of a read's Phred scores, written into the `qs` field.
fn mean_quality(qual: &[u8]) -> f64 {
    let sum: u64 = qual.iter().map(|&q| u64::from(q)).sum();
    sum as f64 / qual.len() as f64
}

/// Inserts an adapter copy at `at`, substituting one base in eight, so the
/// approximate matcher resolves it rather than an exact one.
fn splice(seq: &mut Vec<u8>, at: usize, adapter: &[u8], rng: &mut Rng) {
    let mut copy = adapter.to_vec();
    for _ in 0..adapter.len() / 8 {
        let i = rng.below(copy.len() as u64) as usize;
        copy[i] = b"ACGT"[rng.below(4) as usize];
    }
    let tail = seq.split_off(at);
    seq.extend_from_slice(&copy);
    seq.extend_from_slice(&tail);
}

/// Formats seconds since the Unix epoch as the RFC 3339 timestamp ONT writes.
fn iso_time(epoch_seconds: i64) -> String {
    jiff::Timestamp::from_second(epoch_seconds)
        .expect("training timestamps are in range")
        .to_string()
}

/// Builds the FASTQ text both read files hold.
fn build_fastq() -> Vec<u8> {
    let mut rng = Rng::new(SEED);
    let adapters = preset::preset_ont();
    let five = adapter_seq(&adapters, End::Five);
    let three = adapter_seq(&adapters, End::Three);

    let mut out = Vec::with_capacity(16 << 20);
    for i in 0..FASTQ_READS {
        let len = fastq_length(&mut rng);
        let mut seq = sequence(&mut rng, len);
        // Half the reads carry an adapter, at the 5' end, at the 3' end, or in
        // the interior, so terminal trimming and chimera splitting both run.
        match i % 6 {
            0 => splice(&mut seq, 0, &five, &mut rng),
            1 => {
                let at = seq.len();
                splice(&mut seq, at, &three, &mut rng);
            },
            2 => {
                let at = seq.len() / 2;
                splice(&mut seq, at, &five, &mut rng);
            },
            _ => {},
        }
        let qual = quality(&mut rng, seq.len());

        // Half the headers carry the MinKNOW fields the run panel reads, so both
        // the present and the absent branch run.
        out.push(b'@');
        if i % 2 == 0 {
            let header = format!(
                "read{i} runid=8f0c2ac4 ch={} start_time={} qs={:.1}",
                1 + i % 512,
                iso_time(RUN_START + i as i64 * 7),
                mean_quality(&qual)
            );
            out.extend_from_slice(header.as_bytes());
        } else {
            out.extend_from_slice(format!("read{i}").as_bytes());
        }
        out.push(b'\n');
        out.extend_from_slice(&seq);
        out.extend_from_slice(b"\n+\n");
        out.extend(qual.iter().map(|q| q + b'!'));
        out.push(b'\n');
    }
    out
}

/// Builds a two-group `MM`/`ML` pair: 5mC on every third C and 6mA on every
/// fifth A. Likelihoods are bimodal, so the report's distribution has mass on
/// both sides of its threshold.
fn modifications(seq: &[u8], rng: &mut Rng) -> (Vec<u8>, Vec<u8>) {
    let mut mm = Vec::new();
    let mut ml = Vec::new();
    for (base, code, stride) in [(b'C', &b"C+m"[..], 3usize), (b'A', &b"A+a"[..], 5usize)] {
        mm.extend_from_slice(code);
        mm.push(b'?');
        let mut skipped = 0usize;
        let mut seen = 0usize;
        for &b in seq {
            if b != base {
                continue;
            }
            if seen.is_multiple_of(stride) {
                mm.extend_from_slice(format!(",{skipped}").as_bytes());
                ml.push(if rng.below(2) == 0 {
                    20 + rng.below(60) as u8
                } else {
                    180 + rng.below(70) as u8
                });
                skipped = 0;
            } else {
                skipped += 1;
            }
            seen += 1;
        }
        mm.push(b';');
    }
    (mm, ml)
}

/// Builds a per-base kinetics array, one entry per base.
fn kinetics(rng: &mut Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| 10 + rng.below(200) as u8).collect()
}

/// Builds a move table: the stride, then one block per base plus an idle block
/// every third base. The number of 1s equals the sequence length, which is the
/// invariant the move-table rewrite checks.
fn move_table(len: usize) -> Vec<i8> {
    let mut mv = Vec::with_capacity(len + len / 3 + 1);
    mv.push(MOVE_STRIDE);
    for i in 0..len {
        mv.push(1);
        if i % 3 == 0 {
            mv.push(0);
        }
    }
    mv
}

/// Writes the unaligned BAM.
fn write_bam(path: &Path) -> std::io::Result<()> {
    let mut rng = Rng::new(SEED ^ 0x9e37_79b9_7f4a_7c15);
    let header = sam::Header::default();
    let mut writer = bam::io::Writer::new(std::fs::File::create(path)?);
    writer.write_header(&header)?;

    for i in 0..BAM_RECORDS {
        let len = bam_length(&mut rng);
        let seq = sequence(&mut rng, len);
        let qual = quality(&mut rng, len);

        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(format!("read{i}").into_bytes().into());
        *rec.sequence_mut() = seq.clone().into();
        *rec.quality_scores_mut() = qual.clone().into();

        let (mm, ml) = modifications(&seq, &mut rng);
        let data = rec.data_mut();
        data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.into()));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(ml)),
        );
        data.insert(
            Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
            Value::Int32(len as i32),
        );

        // Per-base kinetics, sliced with the sequence on every trim.
        data.insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(kinetics(&mut rng, len))),
        );
        data.insert(
            Tag::new(b'p', b'w'),
            Value::Array(Array::UInt8(kinetics(&mut rng, len))),
        );

        // `ts` and `ns` agree with the block count, so `--update-moves`
        // recomputes the signal window instead of dropping the tags.
        let moves = move_table(len);
        let blocks = (moves.len() - 1) as i32;
        let samples = SIGNAL_TRIM + blocks * i32::from(MOVE_STRIDE);
        data.insert(Tag::new(b'm', b'v'), Value::Array(Array::Int8(moves)));
        data.insert(Tag::new(b't', b's'), Value::Int32(SIGNAL_TRIM));
        data.insert(Tag::new(b'n', b's'), Value::Int32(samples));

        // A quarter of the records carry the ONT run fields, so the run panel
        // sees reads with and without them. `du` accompanies `st`, which the
        // split path needs to place a segment in time.
        if i % 4 == 0 {
            data.insert(Tag::new(b'c', b'h'), Value::Int32(1 + (i % 512) as i32));
            data.insert(
                Tag::new(b'q', b's'),
                Value::Float(mean_quality(&qual) as f32),
            );
            data.insert(
                Tag::new(b's', b't'),
                Value::String(iso_time(RUN_START + i as i64 * 5).into_bytes().into()),
            );
            data.insert(
                Tag::new(b'd', b'u'),
                Value::Float(samples as f32 / SAMPLE_RATE),
            );
        }
        // An ordinary tag, copied through rather than rewritten.
        data.insert(
            Tag::new(b'R', b'G'),
            Value::String(b"train1".to_vec().into()),
        );

        writer.write_alignment_record(&header, &rec)?;
    }
    writer.try_finish()
}
