use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::RecordBuf};

/// A boxed, owning iterator over decoded `RecordBuf`s (or per-record errors).
type RecordBufIterBox = Box<dyn Iterator<Item = anyhow::Result<RecordBuf>>>;

/// Error (naming the read) if the record is aligned. uBAM only in v1.
pub fn ensure_unaligned(rec: &RecordBuf) -> anyhow::Result<()> {
    if rec.flags().is_unmapped() {
        return Ok(());
    }
    let name = rec
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .unwrap_or_else(|| "<unnamed>".to_string());
    anyhow::bail!(
        "read {name} is aligned (mapped); chopping v1 supports unaligned BAM (uBAM) only"
    )
}

/// Open a BAM reader; return the header and an owning `RecordBuf` iterator.
pub fn reader(input: Option<&Path>) -> anyhow::Result<(sam::Header, RecordBufIterBox)> {
    let inner: Box<dyn io::Read> = match input {
        Some(p) => Box::new(File::open(p)?),
        None => Box::new(io::stdin()),
    };
    let mut r = bam::io::Reader::new(inner);
    let header = r.read_header()?;
    let header_for_iter = header.clone();
    let iter = RecordBufIter { reader: r, header: header_for_iter };
    Ok((header, Box::new(iter)))
}

struct RecordBufIter<R: io::Read> {
    reader: bam::io::Reader<R>,
    header: sam::Header,
}

impl<R: io::Read> Iterator for RecordBufIter<R> {
    type Item = anyhow::Result<RecordBuf>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = RecordBuf::default();
        match self.reader.read_record_buf(&self.header, &mut buf) {
            Ok(0) => None,
            Ok(_) => Some(Ok(buf)),
            Err(e) => Some(Err(e.into())),
        }
    }
}

/// A BAM writer with the (provenance-annotated) header already written.
pub fn writer(
    output: Option<&Path>,
    header: &sam::Header,
) -> anyhow::Result<bam::io::Writer<noodles_bgzf::io::Writer<Box<dyn Write>>>> {
    let inner: Box<dyn Write> = match output {
        Some(p) => Box::new(File::create(p)?),
        None => Box::new(io::stdout()),
    };
    let mut w = bam::io::Writer::new(inner);
    w.write_header(header)?;
    Ok(w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;

    #[test]
    fn unmapped_ok_mapped_rejected() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        assert!(ensure_unaligned(&rec).is_ok());

        *rec.flags_mut() = Flags::empty(); // mapped
        let err = ensure_unaligned(&rec).unwrap_err().to_string();
        assert!(err.contains("r1"));
        assert!(err.contains("aligned"));
    }
}
