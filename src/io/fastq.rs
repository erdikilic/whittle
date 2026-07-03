use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
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
    let buffered = BufReader::new(raw);
    let inner: Box<dyn Read + Send> = if gz {
        Box::new(MultiGzDecoder::new(buffered))
    } else {
        Box::new(buffered)
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

/// Write one output segment as a FASTQ record. On splits (`total_segments > 1`)
/// the id gets a `_segment_N` suffix inserted before any description, matching
/// chopper's convention. `phred` is raw; ASCII is emitted by adding 33.
pub fn write_segment<W: Write>(
    w: &mut W,
    name: &[u8],
    seq: &[u8],
    phred: &[u8],
    total_segments: usize,
    segment_idx: usize,
) -> io::Result<()> {
    w.write_all(b"@")?;
    if total_segments > 1 {
        let (id, desc) = split_head(name);
        w.write_all(id)?;
        write!(w, "_segment_{}", segment_idx + 1)?;
        if let Some(d) = desc {
            w.write_all(b" ")?;
            w.write_all(d)?;
        }
    } else {
        w.write_all(name)?;
    }
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    let ascii: Vec<u8> = phred.iter().map(|&q| q + 33).collect();
    w.write_all(&ascii)?;
    w.write_all(b"\n")
}

fn split_head(name: &[u8]) -> (&[u8], Option<&[u8]>) {
    match name.iter().position(|&b| b == b' ') {
        Some(i) => (&name[..i], Some(&name[i + 1..])),
        None => (name, None),
    }
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
    fn roundtrip_reader_writer() {
        let fq = b"@r1\nACGT\n+\nIIII\n@r2 x\nTT\n+\n!!\n";
        let recs: Vec<ReadRecord> = reader_from_slice(fq).map(|r| r.unwrap()).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, b"r1");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[0].qual, vec![40, 40, 40, 40]); // 'I' = 73 - 33
        assert_eq!(recs[1].qual, vec![0, 0]); // '!' = 33 - 33
    }
}
