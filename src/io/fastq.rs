//! FASTQ reading and writing: streaming record iterators over plain, gzip and BGZF input, and segment writers with optional SAM-style header tags.

use std::io::{self, BufWriter, Read, Write};

use flate2::read::MultiGzDecoder;
use gzp::deflate::Mgzip;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use gzp::{Compression, ZWriter};
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use seq_io::fastq::{Reader, Record};

use crate::config::Config;
use crate::record::ReadRecord;

/// Builds a streaming FASTQ record iterator over an already-open source (e.g. a
/// peeked-and-chained stdin stream), decompressing gzip when `gz` is true.
pub fn reader_from(
    inner: Box<dyn Read + Send>,
    gz: bool,
) -> Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send> {
    let inner: Box<dyn Read + Send> = if gz {
        Box::new(MultiGzDecoder::new(inner))
    } else {
        inner
    };
    Box::new(RecordIter::new(inner))
}

/// Builds a FASTQ iterator over BGZF-compressed input. Unlike ordinary gzip,
/// BGZF is independently framed, so blocks inflate on a private pool of
/// `workers` threads while the parser consumes completed blocks in order.
pub fn reader_from_bgzf(
    inner: Box<dyn Read + Send>,
    workers: usize,
) -> anyhow::Result<Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send>> {
    let inner: Box<dyn Read + Send> = if workers > 1 {
        Box::new(noodles_bgzf::io::MultithreadedReader::with_worker_count(
            std::num::NonZero::new(workers).unwrap_or(std::num::NonZero::<usize>::MIN),
            inner,
        ))
    } else {
        Box::new(noodles_bgzf::io::Reader::new(inner))
    };
    Ok(Box::new(RecordIter::new(inner)))
}

/// Highest raw Phred score a Phred+33 quality byte encodes (`~`, ASCII 126).
const MAX_PHRED33: u8 = 126 - 33;

struct RecordIter<R: Read> {
    reader: Reader<R>,
    /// Records yielded so far.
    count: u64,
    /// Id of the last record yielded, kept for error context.
    last_id: Vec<u8>,
}

impl<R: Read> RecordIter<R> {
    fn new(inner: R) -> Self {
        RecordIter {
            reader: Reader::new(inner),
            count: 0,
            last_id: Vec::new(),
        }
    }

    /// Converts a reader error. An `Io` error is unwrapped to the inner
    /// `io::Error` and given stream-position context: seq_io's wrapper both
    /// displays and sources the same error, which prints its message twice
    /// under anyhow's alternate formatting. Parse errors already carry their
    /// record and line.
    fn describe(&self, e: seq_io::fastq::Error) -> anyhow::Error {
        match e {
            seq_io::fastq::Error::Io(e) => {
                let context = if self.count == 0 {
                    "reading the first FASTQ record".to_string()
                } else {
                    format!(
                        "reading FASTQ record after {}",
                        String::from_utf8_lossy(&self.last_id)
                    )
                };
                anyhow::Error::new(e).context(context)
            },
            other => anyhow::Error::new(other),
        }
    }
}

/// Names the first quality byte outside the Phred+33 range. `qual` holds at
/// least one such byte.
fn invalid_quality(id: &[u8], qual: &[u8]) -> anyhow::Error {
    let (pos, &byte) = qual
        .iter()
        .enumerate()
        .find(|&(_, &b)| b.wrapping_sub(33) > MAX_PHRED33)
        .expect("The caller found a quality byte outside the Phred+33 range");
    anyhow::anyhow!(
        "record {}: quality byte 0x{byte:02x} at position {} is outside the Phred+33 range \
         (ASCII 33..=126)",
        String::from_utf8_lossy(id),
        pos + 1
    )
}

impl<R: Read> Iterator for RecordIter<R> {
    type Item = anyhow::Result<ReadRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        let rec = match self.reader.next()? {
            Ok(rec) => rec,
            Err(e) => return Some(Err(self.describe(e))),
        };
        // One pass converts and validates: a byte below 33 wraps above the
        // maximum, so a single compare per byte flags both ends of the range.
        let raw = rec.qual();
        let mut out_of_range = false;
        let qual: Vec<u8> = raw
            .iter()
            .map(|&b| {
                let q = b.wrapping_sub(33);
                out_of_range |= q > MAX_PHRED33;
                q
            })
            .collect();
        if out_of_range {
            return Some(Err(invalid_quality(rec.id_bytes(), raw)));
        }
        self.count += 1;
        self.last_id.clear();
        self.last_id.extend_from_slice(rec.id_bytes());
        Some(Ok(ReadRecord {
            name: rec.head().to_vec(),
            seq: rec.seq().to_vec(),
            qual,
        }))
    }
}

/// Writes the `@`-prefixed header id for a segment (no trailing newline, no tags).
/// On splits (`total_segments > 1`) the id gets a `_segment_N` suffix inserted
/// before any space-separated description.
fn write_head<W: Write>(
    w: &mut W,
    name: &[u8],
    total_segments: usize,
    segment_idx: usize,
) -> io::Result<()> {
    w.write_all(b"@")?;
    if total_segments > 1 {
        // The suffix follows the read ID, preserving the original delimiter,
        // description, and tab-delimited tags.
        match name.iter().position(|&b| b == b' ' || b == b'\t') {
            Some(i) => {
                w.write_all(&name[..i])?;
                write!(w, "_segment_{}", segment_idx + 1)?;
                w.write_all(&name[i..])?;
            },
            None => {
                w.write_all(name)?;
                write!(w, "_segment_{}", segment_idx + 1)?;
            },
        }
    } else {
        w.write_all(name)?;
    }
    Ok(())
}

/// Writes one output segment as a plain FASTQ record. `phred` is raw; ASCII is
/// emitted by adding 33. Thin wrapper over `write_segment_tagged` with no tags,
/// so the record layout lives in one place.
pub fn write_segment<W: Write>(
    w: &mut W,
    name: &[u8],
    seq: &[u8],
    phred: &[u8],
    total_segments: usize,
    segment_idx: usize,
) -> io::Result<()> {
    write_segment_tagged(w, name, seq, phred, total_segments, segment_idx, b"")
}

/// Writes one output segment like `write_segment`, inserting `tags` (already
/// TAB-prefixed per field, or empty) between the header id and the newline:
/// `@<id>[_segment_N]<tags>`.
pub fn write_segment_tagged<W: Write>(
    w: &mut W,
    name: &[u8],
    seq: &[u8],
    phred: &[u8],
    total_segments: usize,
    segment_idx: usize,
    tags: &[u8],
) -> io::Result<()> {
    write_head(w, name, total_segments, segment_idx)?;
    w.write_all(tags)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    // Phred to ASCII in fixed stack chunks, avoiding a per-segment heap
    // allocation.
    let mut ascii = [0u8; 1024];
    for chunk in phred.chunks(ascii.len()) {
        for (dst, &q) in ascii.iter_mut().zip(chunk) {
            *dst = q.saturating_add(33);
        }
        w.write_all(&ascii[..chunk.len()])?;
    }
    w.write_all(b"\n")
}

/// Appends `n` as ASCII decimal without invoking `core::fmt`. The aux `B` arrays
/// (notably `ML` and the per-base kinetics tags) can hold tens of thousands of
/// integers per record; a stack buffer and a digit loop avoid the per-call
/// overhead of `write!` on a `Vec` at that volume.
#[inline]
pub(crate) fn push_u64(out: &mut Vec<u8>, mut n: u64) {
    let mut buf = [0u8; 20];
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    out.extend_from_slice(&buf[i..]);
}

/// Appends `n` as signed ASCII decimal, the counterpart of [`push_u64`].
/// `unsigned_abs` yields the correct magnitude even for `i64::MIN`, where `-n`
/// would overflow.
#[inline]
pub(crate) fn push_i64(out: &mut Vec<u8>, n: i64) {
    if n < 0 {
        out.push(b'-');
    }
    push_u64(out, n.unsigned_abs());
}

/// Formats one SAM aux field as text `XX:T:VALUE` (no leading TAB). Integers of
/// any source width serialize with SAM type code `i`; `B` arrays keep their
/// subtype.
pub fn format_aux_field(tag: [u8; 2], value: &Value) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&tag);
    out.push(b':');
    match value {
        Value::Character(c) => {
            out.extend_from_slice(b"A:");
            out.push(*c);
        },
        Value::Int8(n) => {
            out.extend_from_slice(b"i:");
            push_i64(&mut out, i64::from(*n));
        },
        Value::UInt8(n) => {
            out.extend_from_slice(b"i:");
            push_u64(&mut out, u64::from(*n));
        },
        Value::Int16(n) => {
            out.extend_from_slice(b"i:");
            push_i64(&mut out, i64::from(*n));
        },
        Value::UInt16(n) => {
            out.extend_from_slice(b"i:");
            push_u64(&mut out, u64::from(*n));
        },
        Value::Int32(n) => {
            out.extend_from_slice(b"i:");
            push_i64(&mut out, i64::from(*n));
        },
        Value::UInt32(n) => {
            out.extend_from_slice(b"i:");
            push_u64(&mut out, u64::from(*n));
        },
        Value::Float(x) => write!(out, "f:{x}").unwrap(),
        Value::String(s) => {
            out.extend_from_slice(b"Z:");
            out.extend_from_slice(AsRef::<[u8]>::as_ref(s));
        },
        Value::Hex(s) => {
            out.extend_from_slice(b"H:");
            out.extend_from_slice(AsRef::<[u8]>::as_ref(s));
        },
        Value::Array(a) => {
            out.extend_from_slice(b"B:");
            write_array(&mut out, a);
        },
    }
    out
}

fn write_array(out: &mut Vec<u8>, a: &Array) {
    match a {
        Array::Int8(v) => {
            out.push(b'c');
            for &x in v {
                out.push(b',');
                push_i64(out, i64::from(x));
            }
        },
        Array::UInt8(v) => {
            out.push(b'C');
            for &x in v {
                out.push(b',');
                push_u64(out, u64::from(x));
            }
        },
        Array::Int16(v) => {
            out.push(b's');
            for &x in v {
                out.push(b',');
                push_i64(out, i64::from(x));
            }
        },
        Array::UInt16(v) => {
            out.push(b'S');
            for &x in v {
                out.push(b',');
                push_u64(out, u64::from(x));
            }
        },
        Array::Int32(v) => {
            out.push(b'i');
            for &x in v {
                out.push(b',');
                push_i64(out, i64::from(x));
            }
        },
        Array::UInt32(v) => {
            out.push(b'I');
            for &x in v {
                out.push(b',');
                push_u64(out, u64::from(x));
            }
        },
        Array::Float(v) => {
            out.push(b'f');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        },
    }
}

/// Formats the reconstructed MM/ML/MN block as SAM aux text (no leading TAB):
/// `MM:Z:<mm>\tML:B:C,<ml...>\tMN:i:<mn>`. `ml` is `None` for an MM-only source
/// record (ML is optional per the SAM spec), in which case the `ML:B:C` field is
/// omitted rather than emitted empty.
pub fn format_mods_aux(mm: &[u8], ml: Option<&[u8]>, mn: usize) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"MM:Z:");
    out.extend_from_slice(mm);
    if let Some(ml) = ml {
        out.extend_from_slice(b"\tML:B:C");
        for &b in ml {
            out.push(b',');
            push_u64(&mut out, u64::from(b));
        }
    }
    out.extend_from_slice(b"\tMN:i:");
    push_u64(&mut out, mn as u64);
    out
}

/// FASTQ output writer: a plain buffered writer, a `gzp` parallel gzip writer
/// for `FastqGz`, or a multithreaded BGZF writer for `FastqBgzf`.
///
/// `gzp`'s `ParCompress` requires an explicit `finish()`: its `Write` impl hands
/// only full chunks to the compressor threads, so the tail block and gzip footer
/// are never flushed by `flush()`. Its `Drop` calls `finish()` as a backstop but
/// `.unwrap()`s the result, turning an I/O error into a panic; calling it
/// explicitly keeps that failure an ordinary `Err`.
pub(crate) enum FastqOut {
    /// Plain buffered output.
    Plain(BufWriter<Box<dyn Write + Send>>),
    /// Parallel gzip (`gzp` `Mgzip`) output.
    Gz(ParCompress<'static, Mgzip, Box<dyn Write + Send>>),
    /// Multithreaded BGZF output.
    Bgzf(noodles_bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>),
}

impl Write for FastqOut {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            FastqOut::Plain(w) => w.write(buf),
            FastqOut::Gz(w) => w.write(buf),
            FastqOut::Bgzf(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            FastqOut::Plain(w) => w.flush(),
            FastqOut::Gz(w) => w.flush(),
            FastqOut::Bgzf(w) => w.flush(),
        }
    }
}

impl FastqOut {
    /// Finalizes the writer: `Gz` flushes the final block and gzip footer
    /// through `ZWriter::finish`, `Bgzf` writes the BGZF EOF block through
    /// `finish`, and `Plain` flushes the `BufWriter`. Must be called before
    /// returning success.
    pub(crate) fn finish(self) -> anyhow::Result<()> {
        match self {
            FastqOut::Plain(mut w) => {
                w.flush()?;
                Ok(())
            },
            FastqOut::Gz(mut w) => {
                w.finish()?;
                Ok(())
            },
            FastqOut::Bgzf(mut w) => {
                w.finish()?;
                Ok(())
            },
        }
    }
}

/// Builds the FASTQ output writer over a file or stdout: a parallel gzip
/// encoder (`gzp`) for `FastqGz`, a multithreaded BGZF writer for `FastqBgzf`,
/// and a plain buffered writer otherwise. `gz_workers` is the caller's encode
/// share of the `-t` budget; both compressed formats clamp it to at least one
/// thread and plain output ignores it.
pub(crate) fn writer(
    cfg: &Config,
    out_fmt: crate::io::Format,
    gz_workers: usize,
) -> anyhow::Result<FastqOut> {
    let base: Box<dyn Write + Send> = match cfg.io.output.as_deref() {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    };
    let workers = std::num::NonZero::new(gz_workers).unwrap_or(std::num::NonZero::<usize>::MIN);
    match out_fmt {
        crate::io::Format::FastqGz => {
            // gzp's `Mgzip` (libdeflate-backed blocked gzip) rather than `Gzip`
            // (flate2/zlib-ng). libdeflater is already linked, and the output is
            // a valid multi-member gzip stream that `MultiGzDecoder` and standard
            // gzip tools decode.
            let w = ParCompressBuilder::<Mgzip>::new()
                .num_threads(workers.get())?
                .compression_level(Compression::new(cfg.compression_level as u32))
                .from_writer(base);
            Ok(FastqOut::Gz(w))
        },
        crate::io::Format::FastqBgzf => {
            let level = noodles_bgzf::io::writer::CompressionLevel::new(cfg.compression_level)
                .ok_or_else(|| anyhow::anyhow!("invalid BGZF compression level"))?;
            let w = noodles_bgzf::io::multithreaded_writer::Builder::default()
                .set_compression_level(level)
                .set_worker_count(workers)
                .build_from_writer(base);
            Ok(FastqOut::Bgzf(w))
        },
        crate::io::Format::Fastq | crate::io::Format::Bam => {
            Ok(FastqOut::Plain(BufWriter::new(base)))
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_single_segment_verbatim_header() {
        let mut out = Vec::new();
        write_segment(&mut out, b"read1 desc", b"ACGT", &[40, 40, 40, 40], 1, 0).unwrap();
        assert_eq!(out, b"@read1 desc\nACGT\n+\nIIII\n");
    }

    #[test]
    fn split_segment_suffixes_id_before_desc() {
        let mut out = Vec::new();
        write_segment(&mut out, b"read1 desc", b"AC", &[40, 40], 2, 1).unwrap();
        assert_eq!(out, b"@read1_segment_2 desc\nAC\n+\nII\n");
    }

    /// A `samtools fastq -T`-style header carries SAM tags after a TAB. The
    /// `_segment_N` suffix lands directly after the read id, keeping the TAB
    /// and the tag value intact.
    #[test]
    fn split_segment_preserves_tab_delimited_tags() {
        let mut out = Vec::new();
        write_segment(&mut out, b"r1\tRG:Z:grp", b"AC", &[40, 40], 2, 0).unwrap();
        assert_eq!(out, b"@r1_segment_1\tRG:Z:grp\nAC\n+\nII\n");
    }

    #[test]
    fn roundtrip_reader_writer() {
        let fq = b"@r1\nACGT\n+\nIIII\n@r2 x\nTT\n+\n!!\n";
        let recs: Vec<ReadRecord> = RecordIter::new(&fq[..]).map(|r| r.unwrap()).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, b"r1");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[0].qual, vec![40, 40, 40, 40]); // 'I' = 73 - 33
        assert_eq!(recs[1].qual, vec![0, 0]); // '!' = 33 - 33
    }

    #[test]
    fn quality_range_bounds_are_accepted() {
        let fq = b"@r1\nAC\n+\n!~\n";
        let recs: Vec<ReadRecord> = RecordIter::new(&fq[..]).map(|r| r.unwrap()).collect();
        assert_eq!(recs[0].qual, vec![0, 93]);
    }

    #[test]
    fn quality_byte_outside_phred33_range_is_an_error() {
        let fq = b"@r1 desc\nACGTACGT\n+\nII I\x01III\n";
        let err = RecordIter::new(&fq[..])
            .next()
            .unwrap()
            .unwrap_err()
            .to_string();
        assert!(err.contains("record r1"), "Names the record id: {err}");
        assert!(
            err.contains("0x20"),
            "Names the first offending byte: {err}"
        );
        assert!(err.contains("Phred+33"), "Names the expected range: {err}");

        let fq = b"@r2\nAC\n+\nI\x7f\n";
        let err = RecordIter::new(&fq[..])
            .next()
            .unwrap()
            .unwrap_err()
            .to_string();
        assert!(err.contains("record r2") && err.contains("0x7f"), "{err}");
    }

    /// Serves `data` and fails every read past its end with `msg`.
    struct FailAfter {
        data: Vec<u8>,
        pos: usize,
        msg: &'static str,
    }

    impl Read for FailAfter {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.pos == self.data.len() {
                return Err(io::Error::other(self.msg));
            }
            let n = (self.data.len() - self.pos).min(buf.len());
            buf[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
            self.pos += n;
            Ok(n)
        }
    }

    #[test]
    fn io_error_carries_record_context_and_prints_the_cause_once() {
        // More than one parser buffer of records, so the failing read follows
        // records that were yielded.
        let mut data = Vec::new();
        for i in 0..8000 {
            data.extend_from_slice(format!("@r{i} desc\nACGT\n+\nIIII\n").as_bytes());
        }
        let mut it = RecordIter::new(FailAfter {
            data,
            pos: 0,
            msg: "incomplete deflate stream",
        });
        let mut last_id = None;
        let err = loop {
            match it.next().unwrap() {
                Ok(rec) => last_id = Some(rec.name[..rec.name.len() - 5].to_vec()),
                Err(e) => break e,
            }
        };
        let last_id = String::from_utf8(last_id.expect("Records precede the error")).unwrap();
        assert_eq!(
            format!("{err:#}"),
            format!("reading FASTQ record after {last_id}: incomplete deflate stream")
        );

        let mut it = RecordIter::new(FailAfter {
            data: Vec::new(),
            pos: 0,
            msg: "boom",
        });
        let err = it.next().unwrap().unwrap_err();
        assert_eq!(format!("{err:#}"), "reading the first FASTQ record: boom");
    }

    #[test]
    fn parse_error_message_is_not_duplicated() {
        let fq = b"@r1\nACGT\n+\nII\n";
        let err = RecordIter::new(&fq[..]).next().unwrap().unwrap_err();
        let msg = format!("{err:#}");
        assert_eq!(msg.matches("FASTQ parse error").count(), 1, "{msg}");
    }

    #[test]
    fn gz_writer_clamps_zero_workers() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("o.fastq.gz");
        let mut cfg = crate::cli::config_for_test(&path, &path, 0, 0);
        cfg.io.output = Some(path);
        let w = writer(&cfg, crate::io::Format::FastqGz, 0).unwrap();
        w.finish().unwrap();
    }

    #[test]
    fn bgzf_reader_roundtrips_fastq() {
        let fq = b"@r1\nACGT\n+\nIIII\n@r2 x\nTT\n+\n!!\n";
        let mut writer = noodles_bgzf::io::Writer::new(Vec::new());
        writer.write_all(fq).unwrap();
        let compressed = writer.finish().unwrap();

        let records: Vec<_> = reader_from_bgzf(Box::new(std::io::Cursor::new(compressed)), 1)
            .unwrap()
            .map(Result::unwrap)
            .collect();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, b"r1");
        assert_eq!(records[0].seq, b"ACGT");
        assert_eq!(records[1].name, b"r2 x");
    }

    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    #[test]
    fn aux_scalar_types() {
        assert_eq!(
            format_aux_field(*b"RG", &Value::String(b"grp1".as_slice().into())),
            b"RG:Z:grp1"
        );
        assert_eq!(format_aux_field(*b"NM", &Value::Int32(-3)), b"NM:i:-3");
        assert_eq!(format_aux_field(*b"Uq", &Value::UInt8(200)), b"Uq:i:200");
        assert_eq!(format_aux_field(*b"pa", &Value::Float(0.5)), b"pa:f:0.5");
        assert_eq!(format_aux_field(*b"bc", &Value::Character(b'K')), b"bc:A:K");
        assert_eq!(
            format_aux_field(*b"H2", &Value::Hex(b"1AE3".as_slice().into())),
            b"H2:H:1AE3"
        );
        // Every integer width serializes with SAM type code `i`, regardless of
        // signedness or size.
        assert_eq!(format_aux_field(*b"i1", &Value::Int8(-5)), b"i1:i:-5");
        assert_eq!(format_aux_field(*b"i2", &Value::Int16(-300)), b"i2:i:-300");
        assert_eq!(format_aux_field(*b"i3", &Value::UInt16(400)), b"i3:i:400");
        assert_eq!(
            format_aux_field(*b"i4", &Value::UInt32(70000)),
            b"i4:i:70000"
        );
    }

    #[test]
    fn aux_array_subtypes() {
        assert_eq!(
            format_aux_field(*b"a1", &Value::Array(Array::UInt8(vec![1, 2, 3]))),
            b"a1:B:C,1,2,3"
        );
        assert_eq!(
            format_aux_field(*b"a2", &Value::Array(Array::Int8(vec![-1, 2]))),
            b"a2:B:c,-1,2"
        );
        assert_eq!(
            format_aux_field(*b"a3", &Value::Array(Array::Int16(vec![-5]))),
            b"a3:B:s,-5"
        );
        assert_eq!(
            format_aux_field(*b"a4", &Value::Array(Array::UInt16(vec![5]))),
            b"a4:B:S,5"
        );
        assert_eq!(
            format_aux_field(*b"a5", &Value::Array(Array::Int32(vec![7]))),
            b"a5:B:i,7"
        );
        assert_eq!(
            format_aux_field(*b"a6", &Value::Array(Array::UInt32(vec![8]))),
            b"a6:B:I,8"
        );
        assert_eq!(
            format_aux_field(*b"a7", &Value::Array(Array::Float(vec![1.5]))),
            b"a7:B:f,1.5"
        );
    }

    #[test]
    fn mods_aux_layout() {
        assert_eq!(
            format_mods_aux(b"C+m,0;", Some(&[10, 20]), 6),
            b"MM:Z:C+m,0;\tML:B:C,10,20\tMN:i:6"
        );
        // ML present but empty (e.g. all mods sliced away yet MM retained) yields
        // a zero-length `B:C` array.
        assert_eq!(
            format_mods_aux(b"C+m;", Some(&[]), 4),
            b"MM:Z:C+m;\tML:B:C\tMN:i:4"
        );
        // ML absent (MM-only source record): the ML field is omitted rather than
        // emitted empty, so the record stays valid.
        assert_eq!(format_mods_aux(b"C+m,0;", None, 4), b"MM:Z:C+m,0;\tMN:i:4");
    }

    #[test]
    fn tagged_writer_appends_tags_after_id() {
        let mut out = Vec::new();
        write_segment_tagged(
            &mut out,
            b"read2",
            b"AC",
            &[40, 40],
            1,
            0,
            b"\tRG:Z:grp1\tMM:Z:C+m,0;\tML:B:C,20\tMN:i:2",
        )
        .unwrap();
        assert_eq!(
            out,
            b"@read2\tRG:Z:grp1\tMM:Z:C+m,0;\tML:B:C,20\tMN:i:2\nAC\n+\nII\n"
        );
    }

    #[test]
    fn tagged_writer_empty_tags_is_plain_record() {
        let mut a = Vec::new();
        write_segment_tagged(&mut a, b"read1", b"ACGT", &[40, 40, 40, 40], 1, 0, b"").unwrap();
        let mut b = Vec::new();
        write_segment(&mut b, b"read1", b"ACGT", &[40, 40, 40, 40], 1, 0).unwrap();
        assert_eq!(a, b);
        assert_eq!(a, b"@read1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn tagged_writer_split_suffix_then_tags() {
        let mut out = Vec::new();
        write_segment_tagged(&mut out, b"read2", b"AC", &[40, 40], 2, 1, b"\tMN:i:2").unwrap();
        assert_eq!(out, b"@read2_segment_2\tMN:i:2\nAC\n+\nII\n");
    }
}
