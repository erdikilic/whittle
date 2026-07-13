use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use seq_io::fastq::{Reader, Record};

use crate::record::ReadRecord;

/// Build a streaming FASTQ record iterator over a file (or stdin when `input`
/// is `None`), transparently decompressing gzip when `gz` is true.
pub fn reader(
    input: Option<&Path>,
    gz: bool,
) -> anyhow::Result<Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send>> {
    let raw: Box<dyn Read + Send> = match input {
        Some(p) => Box::new(File::open(p)?),
        None => Box::new(io::stdin()),
    };
    let buffered: Box<dyn Read + Send> = Box::new(BufReader::new(raw));
    Ok(reader_from(buffered, gz))
}

/// Build a streaming FASTQ record iterator over an already-open source (e.g. a
/// peeked-and-chained stdin stream), transparently decompressing gzip when
/// `gz` is true.
pub fn reader_from(
    inner: Box<dyn Read + Send>,
    gz: bool,
) -> Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send> {
    let inner: Box<dyn Read + Send> = if gz {
        Box::new(MultiGzDecoder::new(inner))
    } else {
        inner
    };
    Box::new(RecordIter {
        reader: Reader::new(inner),
    })
}

/// Build a FASTQ iterator over BGZF-compressed input. Unlike ordinary gzip,
/// BGZF is independently framed and can inflate blocks on the shared Rayon
/// codec pool while the parser consumes completed blocks in order.
pub fn reader_from_bgzf(
    inner: Box<dyn Read + Send>,
    workers: usize,
) -> anyhow::Result<Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send>> {
    let inner: Box<dyn Read + Send> = if workers > 1 {
        crate::io::bam::configure_bgzf_pool(workers)?;
        Box::new(noodles_bgzf::io::MultithreadedReader::new(inner))
    } else {
        Box::new(noodles_bgzf::io::Reader::new(inner))
    };
    Ok(Box::new(RecordIter {
        reader: Reader::new(inner),
    }))
}

struct RecordIter<R: Read> {
    reader: Reader<R>,
}

impl<R: Read> Iterator for RecordIter<R> {
    type Item = anyhow::Result<ReadRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        let rec = self.reader.next()?;
        Some(rec.map_err(anyhow::Error::from).map(|r| ReadRecord {
            name: r.head().to_vec(),
            seq: r.seq().to_vec(),
            qual: r.qual().iter().map(|&b| b.saturating_sub(33)).collect(),
        }))
    }
}

/// Write the `@`-prefixed header id for a segment (no trailing newline, no tags).
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
        // Insert the suffix after the read ID while preserving the original
        // delimiter, description, and tab-delimited tags.
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

/// Write one output segment as a plain FASTQ record. `phred` is raw; ASCII is
/// emitted by adding 33. Thin wrapper over `write_segment_tagged` with no tags,
/// so the record layout lives in exactly one place.
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

/// Like `write_segment`, but inserts `tags` (already TAB-prefixed per field, or
/// empty) between the header id and the newline: `@<id>[_segment_N]<tags>`.
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
    // Encode phred -> ASCII in fixed stack chunks, avoiding a per-segment heap
    // allocation (this runs once per output segment).
    let mut ascii = [0u8; 1024];
    for chunk in phred.chunks(ascii.len()) {
        for (dst, &q) in ascii.iter_mut().zip(chunk) {
            *dst = q.saturating_add(33);
        }
        w.write_all(&ascii[..chunk.len()])?;
    }
    w.write_all(b"\n")
}

/// One SAM aux field as text `XX:T:VALUE` (no leading TAB). Integers of any
/// source width serialize with SAM type code `i`; `B` arrays keep their subtype.
pub fn format_aux_field(tag: [u8; 2], value: &Value) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&tag);
    out.push(b':');
    match value {
        Value::Character(c) => {
            out.extend_from_slice(b"A:");
            out.push(*c);
        },
        Value::Int8(n) => write!(out, "i:{n}").unwrap(),
        Value::UInt8(n) => write!(out, "i:{n}").unwrap(),
        Value::Int16(n) => write!(out, "i:{n}").unwrap(),
        Value::UInt16(n) => write!(out, "i:{n}").unwrap(),
        Value::Int32(n) => write!(out, "i:{n}").unwrap(),
        Value::UInt32(n) => write!(out, "i:{n}").unwrap(),
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
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        },
        Array::UInt8(v) => {
            out.push(b'C');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        },
        Array::Int16(v) => {
            out.push(b's');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        },
        Array::UInt16(v) => {
            out.push(b'S');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        },
        Array::Int32(v) => {
            out.push(b'i');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        },
        Array::UInt32(v) => {
            out.push(b'I');
            for x in v {
                write!(out, ",{x}").unwrap();
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

/// The reconstructed MM/ML/MN block as SAM aux text (no leading TAB):
/// `MM:Z:<mm>\tML:B:C,<ml…>\tMN:i:<mn>`. `ml` is `None` for an MM-only source
/// record (ML is optional per the SAM spec), in which case the `ML:B:C` field is
/// omitted entirely rather than emitted empty.
pub fn format_mods_aux(mm: &[u8], ml: Option<&[u8]>, mn: usize) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"MM:Z:");
    out.extend_from_slice(mm);
    if let Some(ml) = ml {
        out.extend_from_slice(b"\tML:B:C");
        for b in ml {
            write!(out, ",{b}").unwrap();
        }
    }
    write!(out, "\tMN:i:{mn}").unwrap();
    out
}

#[cfg(test)]
fn reader_from_slice(bytes: &'static [u8]) -> RecordIter<&'static [u8]> {
    RecordIter {
        reader: Reader::new(bytes),
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

    #[test]
    fn split_segment_preserves_tab_delimited_tags() {
        // A `samtools fastq -T`-style header carries SAM tags after a TAB. The
        // `_segment_N` suffix must land right after the read id, keeping the TAB
        // and the tag value intact — not append past the tag (which would mutate
        // both the id and the RG value).
        let mut out = Vec::new();
        write_segment(&mut out, b"r1\tRG:Z:grp", b"AC", &[40, 40], 2, 0).unwrap();
        assert_eq!(out, b"@r1_segment_1\tRG:Z:grp\nAC\n+\nII\n");
    }

    #[test]
    fn roundtrip_reader_writer() {
        let fq = b"@r1\nACGT\n+\nIIII\n@r2 x\nTT\n+\n!!\n";
        let recs: Vec<ReadRecord> = reader_from_slice(fq).map(|r| r.unwrap()).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, b"r1");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[0].qual, vec![40, 40, 40, 40]); // 'I' = 73 - 33
        assert_eq!(recs[1].qual, vec![0, 0]); // '!' = 33 - 33
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
        // ML present but empty (e.g. all mods sliced away yet MM retained) -> zero-length B:C array
        assert_eq!(
            format_mods_aux(b"C+m;", Some(&[]), 4),
            b"MM:Z:C+m;\tML:B:C\tMN:i:4"
        );
        // ML absent (MM-only source record) -> the ML field is omitted entirely,
        // never emitted empty, so the record stays valid.
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
