//! FASTQ reading and writing: streaming record iterators over plain, gzip and BGZF input, and segment writers with optional SAM-style header tags.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

use anyhow::Context;

use flate2::bufread::MultiGzDecoder;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use seq_io::fastq::{Reader, Record};

use crate::config::Config;
use crate::record::ReadRecord;

/// Capacity of the buffer between a compressed source and its decoder, which
/// pulls the compressed stream through it in one `read` syscall per megabyte.
///
/// The parser's own buffer keeps the seq_io default: the parser fills it to
/// capacity before yielding, so an I/O error during a fill discards every
/// record in it, and a larger buffer would widen the span of records an error
/// on damaged input cannot be attributed to.
const INPUT_BUFFER_CAPACITY: usize = 1 << 20;

/// Builds a streaming FASTQ record iterator over an already-open source (e.g. a
/// peeked-and-chained stdin stream), decompressing gzip when `gz` is true.
pub fn reader_from(
    inner: Box<dyn Read + Send>,
    gz: bool,
) -> Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send> {
    let inner: Box<dyn Read + Send> = if gz {
        // `bufread::MultiGzDecoder` over an explicit buffer; the `read` variant
        // wraps its source in an 8 KiB `BufReader`.
        Box::new(MultiGzDecoder::new(BufReader::with_capacity(
            INPUT_BUFFER_CAPACITY,
            inner,
        )))
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
    // The block reader issues a `read` per frame part (header, payload,
    // trailer); the buffer in front of the source coalesces them.
    let inner = BufReader::with_capacity(INPUT_BUFFER_CAPACITY, inner);
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
/// emitted by adding 33.
pub fn write_segment<W: Write>(
    w: &mut W,
    name: &[u8],
    seq: &[u8],
    phred: &[u8],
    total_segments: usize,
    segment_idx: usize,
) -> io::Result<()> {
    write_head(w, name, total_segments, segment_idx)?;
    write_body(w, seq, phred)
}

/// Appends the rest of a record after its header id and tags to `out`: the
/// newline, the sequence, the `+` line and the Phred+33 qualities. The `Vec`
/// counterpart of `write_body` for callers that assemble the header in place.
pub(crate) fn push_record_body(out: &mut Vec<u8>, seq: &[u8], phred: &[u8]) {
    write_body(out, seq, phred).expect("Writing to a Vec cannot fail");
}

/// Writes the rest of a record after its header line: the newline, the
/// sequence, the `+` line and the Phred+33 qualities.
fn write_body<W: Write>(w: &mut W, seq: &[u8], phred: &[u8]) -> io::Result<()> {
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

/// The decimal text of `0..=99` as digit pairs, for the two-digit fast path of
/// [`push_u64`]: `MM` skip counts are mostly one or two digits.
const DIGIT_PAIRS: [[u8; 2]; 100] = {
    let mut t = [[0u8; 2]; 100];
    let mut i = 0;
    while i < 100 {
        t[i] = [b'0' + (i / 10) as u8, b'0' + (i % 10) as u8];
        i += 1;
    }
    t
};

/// Appends `n` as ASCII decimal without invoking `core::fmt`. The aux `B` arrays
/// (notably `ML` and the per-base kinetics tags) can hold tens of thousands of
/// integers per record; a stack buffer and a digit loop avoid the per-call
/// overhead of `write!` on a `Vec` at that volume.
#[inline]
pub(crate) fn push_u64(out: &mut Vec<u8>, mut n: u64) {
    if n < 10 {
        out.push(b'0' + n as u8);
        return;
    }
    if n < 100 {
        out.extend_from_slice(&DIGIT_PAIRS[n as usize]);
        return;
    }
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

/// Appends one SAM aux field as text `XX:T:VALUE` (no leading TAB) to `out`.
/// Integers of any source width serialize with SAM type code `i`; `B` arrays
/// keep their subtype.
pub fn push_aux_field(out: &mut Vec<u8>, tag: [u8; 2], value: &Value) {
    out.extend_from_slice(&tag);
    out.push(b':');
    match value {
        Value::Character(c) => {
            out.extend_from_slice(b"A:");
            out.push(*c);
        },
        Value::Int8(n) => {
            out.extend_from_slice(b"i:");
            push_i64(out, i64::from(*n));
        },
        Value::UInt8(n) => {
            out.extend_from_slice(b"i:");
            push_u64(out, u64::from(*n));
        },
        Value::Int16(n) => {
            out.extend_from_slice(b"i:");
            push_i64(out, i64::from(*n));
        },
        Value::UInt16(n) => {
            out.extend_from_slice(b"i:");
            push_u64(out, u64::from(*n));
        },
        Value::Int32(n) => {
            out.extend_from_slice(b"i:");
            push_i64(out, i64::from(*n));
        },
        Value::UInt32(n) => {
            out.extend_from_slice(b"i:");
            push_u64(out, u64::from(*n));
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
            write_array(out, a);
        },
    }
}

/// One `B:C` array element as text: `,` and the decimal digits of the value,
/// with the length of the used prefix. Indexed by value. `ML` and the
/// per-base kinetics arrays hold tens of thousands of elements per record, so
/// each element is a fixed-width table copy rather than a digit loop.
const COMMA_U8: [([u8; 4], usize); 256] = {
    let mut t = [([0u8; 4], 0usize); 256];
    let mut i = 0;
    while i < 256 {
        let n = i as u8;
        let mut b = [b',', 0, 0, 0];
        let len = if n >= 100 {
            b[1] = b'0' + n / 100;
            b[2] = b'0' + (n / 10) % 10;
            b[3] = b'0' + n % 10;
            4
        } else if n >= 10 {
            b[1] = b'0' + n / 10;
            b[2] = b'0' + n % 10;
            3
        } else {
            b[1] = b'0' + n;
            2
        };
        t[i] = (b, len);
        i += 1;
    }
    t
};

/// One `B:c` array element as text, the signed counterpart of [`COMMA_U8`],
/// indexed by the value's bit pattern (`x as u8`). The move table `mv` is a
/// `B:c` array with one element per signal stride.
const COMMA_I8: [([u8; 5], usize); 256] = {
    let mut t = [([0u8; 5], 0usize); 256];
    let mut i = 0;
    while i < 256 {
        let x = i as u8 as i8;
        let m = x.unsigned_abs();
        let mut b = [b',', 0, 0, 0, 0];
        let mut len = 1;
        if x < 0 {
            b[len] = b'-';
            len += 1;
        }
        if m >= 100 {
            b[len] = b'0' + m / 100;
            b[len + 1] = b'0' + (m / 10) % 10;
            b[len + 2] = b'0' + m % 10;
            len += 3;
        } else if m >= 10 {
            b[len] = b'0' + m / 10;
            b[len + 1] = b'0' + m % 10;
            len += 2;
        } else {
            b[len] = b'0' + m;
            len += 1;
        }
        t[i] = (b, len);
        i += 1;
    }
    t
};

/// Appends `,` and the decimal text of every element of a `B:C` array. The
/// capacity is reserved once for the widest element, so each element is a
/// fixed four-byte copy followed by a truncate rather than a variable-length
/// copy.
fn push_comma_u8s(out: &mut Vec<u8>, v: &[u8]) {
    out.reserve(v.len() * 4);
    for &x in v {
        let (bytes, len) = COMMA_U8[usize::from(x)];
        out.extend_from_slice(&bytes);
        out.truncate(out.len() - (4 - len));
    }
}

/// Appends `,` and the decimal text of every element of a `B:c` array; see
/// [`push_comma_u8s`].
fn push_comma_i8s(out: &mut Vec<u8>, v: &[i8]) {
    out.reserve(v.len() * 5);
    for &x in v {
        let (bytes, len) = COMMA_I8[usize::from(x as u8)];
        out.extend_from_slice(&bytes);
        out.truncate(out.len() - (5 - len));
    }
}

/// Appends `,` and the decimal text of every element of a wider unsigned `B`
/// array. The capacity is reserved at two bytes per element, the shortest text
/// an element can take.
fn push_comma_u64s(out: &mut Vec<u8>, v: impl ExactSizeIterator<Item = u64>) {
    out.reserve(v.len() * 2);
    for x in v {
        out.push(b',');
        push_u64(out, x);
    }
}

/// Appends `,` and the decimal text of every element of a wider signed `B`
/// array; see [`push_comma_u64s`].
fn push_comma_i64s(out: &mut Vec<u8>, v: impl ExactSizeIterator<Item = i64>) {
    out.reserve(v.len() * 2);
    for x in v {
        out.push(b',');
        push_i64(out, x);
    }
}

/// Appends a `B` array's subtype code and comma-prefixed elements.
fn write_array(out: &mut Vec<u8>, a: &Array) {
    match a {
        Array::Int8(v) => {
            out.push(b'c');
            push_comma_i8s(out, v);
        },
        Array::UInt8(v) => {
            out.push(b'C');
            push_comma_u8s(out, v);
        },
        Array::Int16(v) => {
            out.push(b's');
            push_comma_i64s(out, v.iter().map(|&x| i64::from(x)));
        },
        Array::UInt16(v) => {
            out.push(b'S');
            push_comma_u64s(out, v.iter().map(|&x| u64::from(x)));
        },
        Array::Int32(v) => {
            out.push(b'i');
            push_comma_i64s(out, v.iter().map(|&x| i64::from(x)));
        },
        Array::UInt32(v) => {
            out.push(b'I');
            push_comma_u64s(out, v.iter().map(|&x| u64::from(x)));
        },
        Array::Float(v) => {
            out.push(b'f');
            out.reserve(v.len() * 2);
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        },
    }
}

/// Appends the reconstructed MM/ML/MN block as SAM aux text to `out`, each
/// field TAB-prefixed: `\tMM:Z:<mm>\tML:B:C,<ml...>\tMN:i:<mn>`. `ml` is
/// `None` for an MM-only source record (ML is optional per the SAM spec), in
/// which case the `ML:B:C` field is omitted rather than emitted empty.
///
/// A field named by `remove` is left out and the others still emitted, which is
/// what BAM output does with the same setting. Every field TAB-prefixed rather
/// than TAB-separated is what lets any one of them be omitted.
pub fn push_mods_aux(
    out: &mut Vec<u8>,
    mm: &[u8],
    ml: Option<&[u8]>,
    mn: usize,
    remove: &crate::config::TagRemoval,
) {
    out.reserve(mm.len() + 40);
    if !remove.contains(b"MM") {
        out.extend_from_slice(b"\tMM:Z:");
        out.extend_from_slice(mm);
    }
    if let Some(ml) = ml
        && !remove.contains(b"ML")
    {
        out.extend_from_slice(b"\tML:B:C");
        push_comma_u8s(out, ml);
    }
    if !remove.contains(b"MN") {
        out.extend_from_slice(b"\tMN:i:");
        push_u64(out, mn as u64);
    }
}

/// Capacity of the buffer in front of the output file or stdout. Plain output
/// is written record by record and the compressed encoders emit one write per
/// block; the buffer coalesces either into one `write` syscall per megabyte.
const OUTPUT_BUFFER_CAPACITY: usize = 1 << 20;

/// The buffered destination of a FASTQ writer.
type BufferedOutput = BufWriter<Box<dyn Write + Send>>;

/// FASTQ output writer: a plain buffered writer, or BGZF for the compressed
/// formats. `FastqGz` and `FastqBgzf` share the BGZF framing: BGZF is a valid
/// multi-member gzip stream, so `.fastq.gz` output decodes with every gzip
/// tool and `MultiGzDecoder`, and also inflates block-parallel when read
/// back. `Bgzf` compresses in the writing thread for a sequential run;
/// `Blocks` takes blocks the render workers compressed (see
/// `workflow::run_parallel`).
pub(crate) enum FastqOut {
    /// Plain buffered output.
    Plain(BufferedOutput),
    /// BGZF output compressed in the writing thread.
    Bgzf(noodles_bgzf::io::Writer<BufferedOutput>),
    /// Pre-compressed BGZF blocks.
    Blocks {
        /// The output.
        inner: BufferedOutput,
        /// The BGZF DEFLATE level the render workers compress at.
        level: u8,
    },
}

impl Write for FastqOut {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            FastqOut::Plain(w) | FastqOut::Blocks { inner: w, .. } => w.write(buf),
            FastqOut::Bgzf(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            FastqOut::Plain(w) | FastqOut::Blocks { inner: w, .. } => w.flush(),
            FastqOut::Bgzf(w) => w.flush(),
        }
    }
}

impl FastqOut {
    /// The BGZF level the render workers compress at, for a `Blocks` sink.
    pub(crate) fn block_level(&self) -> Option<u8> {
        match self {
            FastqOut::Blocks { level, .. } => Some(*level),
            FastqOut::Plain(_) | FastqOut::Bgzf(_) => None,
        }
    }

    /// Finalizes the writer: the BGZF variants write their last block and the
    /// BGZF EOF block, and every variant then flushes the output buffer, whose
    /// write error surfaces here. Must be called before returning success;
    /// the encoder's `Drop` swallows errors.
    pub(crate) fn finish(self) -> anyhow::Result<()> {
        let mut inner = match self {
            FastqOut::Plain(w) => w,
            FastqOut::Bgzf(w) => w.finish()?,
            FastqOut::Blocks { mut inner, .. } => {
                inner.write_all(&crate::io::bam::eof_block())?;
                inner
            },
        };
        inner.flush()?;
        Ok(())
    }
}

/// Compresses `data` into BGZF blocks at `level`: the bytes of a BGZF stream
/// fragment with complete blocks and no EOF block, which a `FastqOut::Blocks`
/// writes through.
pub(crate) fn encode_blocks(level: u8, data: &[u8]) -> std::io::Result<Vec<u8>> {
    let clevel = crate::io::bam::compression_level(level).map_err(std::io::Error::other)?;
    let mut w = noodles_bgzf::io::writer::Builder::default()
        .set_compression_level(clevel)
        .build_from_writer(Vec::new());
    w.write_all(data)?;
    w.flush()?;
    Ok(w.into_inner())
}

/// Builds the FASTQ output writer over a file or stdout: for `FastqGz` and
/// `FastqBgzf`, a BGZF writer at `cfg.compression_level` compressing in the
/// writing thread, or a `Blocks` sink when `parallel`; a plain buffered
/// writer for `Fastq`.
pub(crate) fn writer(
    cfg: &Config,
    out_fmt: crate::io::Format,
    parallel: bool,
) -> anyhow::Result<FastqOut> {
    let inner: Box<dyn Write + Send> = match cfg.io.output.as_deref() {
        Some(p) => {
            Box::new(File::create(p).with_context(|| format!("creating output {}", p.display()))?)
        },
        None => Box::new(io::stdout()),
    };
    let base = BufWriter::with_capacity(OUTPUT_BUFFER_CAPACITY, inner);
    match out_fmt {
        crate::io::Format::FastqGz | crate::io::Format::FastqBgzf => {
            let level = cfg.compression_level;
            if parallel {
                return Ok(FastqOut::Blocks { inner: base, level });
            }
            let w = noodles_bgzf::io::writer::Builder::default()
                .set_compression_level(crate::io::bam::compression_level(level)?)
                .build_from_writer(base);
            Ok(FastqOut::Bgzf(w))
        },
        crate::io::Format::Fastq => Ok(FastqOut::Plain(base)),
        crate::io::Format::Bam => unreachable!("BAM output is written by `io::bam::writer`"),
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

    /// The sequential `.fastq.gz` writer emits BGZF-framed gzip that a plain
    /// gzip decoder reads back.
    #[test]
    fn gz_writer_clamps_zero_workers_and_writes_bgzf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("o.fastq.gz");
        let mut cfg = crate::cli::config_for_test(&path, &path, 0, 0);
        cfg.io.output = Some(path.clone());
        let mut w = writer(&cfg, crate::io::Format::FastqGz, false).unwrap();
        w.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        w.finish().unwrap();

        assert!(crate::io::is_bgzf_file(&path));
        let mut text = String::new();
        MultiGzDecoder::new(BufReader::new(std::fs::File::open(&path).unwrap()))
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text, "@r1\nACGT\n+\nIIII\n");
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

    /// One aux field in a buffer of its own, so an expectation reads as the
    /// field text rather than as the state of a shared buffer.
    fn aux(tag: [u8; 2], value: &Value) -> Vec<u8> {
        let mut out = Vec::new();
        push_aux_field(&mut out, tag, value);
        out
    }

    /// The MM/ML/MN block in a buffer of its own, for the same reason.
    fn mods_aux(
        mm: &[u8],
        ml: Option<&[u8]>,
        mn: usize,
        remove: &crate::config::TagRemoval,
    ) -> Vec<u8> {
        let mut out = Vec::new();
        push_mods_aux(&mut out, mm, ml, mn, remove);
        out
    }

    #[test]
    fn aux_scalar_types() {
        assert_eq!(
            aux(*b"RG", &Value::String(b"grp1".as_slice().into())),
            b"RG:Z:grp1"
        );
        assert_eq!(aux(*b"NM", &Value::Int32(-3)), b"NM:i:-3");
        assert_eq!(aux(*b"Uq", &Value::UInt8(200)), b"Uq:i:200");
        assert_eq!(aux(*b"pa", &Value::Float(0.5)), b"pa:f:0.5");
        assert_eq!(aux(*b"bc", &Value::Character(b'K')), b"bc:A:K");
        assert_eq!(
            aux(*b"H2", &Value::Hex(b"1AE3".as_slice().into())),
            b"H2:H:1AE3"
        );
        // Every integer width serializes with SAM type code `i`, regardless of
        // signedness or size.
        assert_eq!(aux(*b"i1", &Value::Int8(-5)), b"i1:i:-5");
        assert_eq!(aux(*b"i2", &Value::Int16(-300)), b"i2:i:-300");
        assert_eq!(aux(*b"i3", &Value::UInt16(400)), b"i3:i:400");
        assert_eq!(aux(*b"i4", &Value::UInt32(70000)), b"i4:i:70000");
    }

    #[test]
    fn aux_array_subtypes() {
        assert_eq!(
            aux(*b"a1", &Value::Array(Array::UInt8(vec![1, 2, 3]))),
            b"a1:B:C,1,2,3"
        );
        assert_eq!(
            aux(*b"a2", &Value::Array(Array::Int8(vec![-1, 2]))),
            b"a2:B:c,-1,2"
        );
        assert_eq!(
            aux(*b"a3", &Value::Array(Array::Int16(vec![-5]))),
            b"a3:B:s,-5"
        );
        assert_eq!(
            aux(*b"a4", &Value::Array(Array::UInt16(vec![5]))),
            b"a4:B:S,5"
        );
        assert_eq!(
            aux(*b"a5", &Value::Array(Array::Int32(vec![7]))),
            b"a5:B:i,7"
        );
        assert_eq!(
            aux(*b"a6", &Value::Array(Array::UInt32(vec![8]))),
            b"a6:B:I,8"
        );
        assert_eq!(
            aux(*b"a7", &Value::Array(Array::Float(vec![1.5]))),
            b"a7:B:f,1.5"
        );
    }

    /// `push_u64` agrees with `core::fmt` across the fast paths and the loop.
    /// The byte-element tables print every value as the digit loop does.
    #[test]
    fn comma_element_tables_match_the_digit_loop() {
        let all: Vec<u8> = (0..=255).collect();
        let mut fast = Vec::new();
        push_comma_u8s(&mut fast, &all);
        let mut slow = Vec::new();
        push_comma_u64s(&mut slow, all.iter().map(|&x| u64::from(x)));
        assert_eq!(fast, slow);

        let all: Vec<i8> = (-128..=127).collect();
        let mut fast = Vec::new();
        push_comma_i8s(&mut fast, &all);
        let mut slow = Vec::new();
        push_comma_i64s(&mut slow, all.iter().map(|&x| i64::from(x)));
        assert_eq!(fast, slow);
    }

    #[test]
    fn push_u64_matches_fmt() {
        for n in (0..1200).chain([9_999, 10_000, u64::from(u32::MAX), u64::MAX]) {
            let mut out = Vec::new();
            push_u64(&mut out, n);
            assert_eq!(out, n.to_string().as_bytes(), "{n}");
        }
        let mut out = Vec::new();
        push_i64(&mut out, i64::MIN);
        assert_eq!(out, i64::MIN.to_string().as_bytes());
    }

    #[test]
    fn mods_aux_layout() {
        let keep = crate::config::TagRemoval::default();
        assert_eq!(
            mods_aux(b"C+m,0;", Some(&[10, 20]), 6, &keep),
            b"\tMM:Z:C+m,0;\tML:B:C,10,20\tMN:i:6"
        );
        // ML present but empty (e.g. all mods sliced away yet MM retained) yields
        // a zero-length `B:C` array.
        assert_eq!(
            mods_aux(b"C+m;", Some(&[]), 4, &keep),
            b"\tMM:Z:C+m;\tML:B:C\tMN:i:4"
        );
        // ML absent (MM-only source record): the ML field is omitted rather than
        // emitted empty, so the record stays valid.
        assert_eq!(
            mods_aux(b"C+m,0;", None, 4, &keep),
            b"\tMM:Z:C+m,0;\tMN:i:4"
        );
        // Each removed field goes on its own, leaving the rest of the block.
        let removal =
            |tag: &str| crate::config::TagRemoval::parse(&[tag.to_string()], false).unwrap();
        assert_eq!(
            mods_aux(b"C+m,0;", Some(&[10]), 6, &removal("ML")),
            b"\tMM:Z:C+m,0;\tMN:i:6"
        );
        assert_eq!(
            mods_aux(b"C+m,0;", Some(&[10]), 6, &removal("MN")),
            b"\tMM:Z:C+m,0;\tML:B:C,10"
        );
        assert_eq!(
            mods_aux(b"C+m,0;", Some(&[10]), 6, &removal("MM")),
            b"\tML:B:C,10\tMN:i:6"
        );
    }
}
