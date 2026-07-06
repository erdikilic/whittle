use std::fs::File;
use std::io::{self, Write};
use std::num::NonZero;
use std::path::Path;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _; // write_header / write_alignment_record
use noodles_sam::{self as sam};

/// A boxed, owning iterator over decoded `RecordBuf`s (or per-record errors).
type RecordBufIterBox = Box<dyn Iterator<Item = anyhow::Result<RecordBuf>> + Send>;

/// Error (naming the read) if the record is aligned. uBAM only in v1.
pub fn ensure_unaligned(rec: &RecordBuf) -> anyhow::Result<()> {
    if rec.flags().is_unmapped() {
        return Ok(());
    }
    let name = rec
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .unwrap_or_else(|| "<unnamed>".to_string());
    anyhow::bail!("read {name} is aligned (mapped); whittle v1 supports unaligned BAM (uBAM) only")
}

/// Open a BAM reader; MT-bgzf when `workers > 1`. Returns the header and a Send
/// owning `RecordBuf` iterator.
pub fn reader(
    input: Option<&Path>,
    workers: usize,
) -> anyhow::Result<(sam::Header, RecordBufIterBox)> {
    let inner: Box<dyn io::Read + Send> = match input {
        Some(p) => Box::new(File::open(p)?),
        None => Box::new(io::stdin()),
    };
    reader_from(inner, workers)
}

/// Like `reader`, but over an already-open stream rather than a path/stdin. Used
/// by the single-file dispatch so a stdin BAM whose first bytes were consumed for
/// format sniffing (and chained back into `inner`) is read from the true start —
/// re-opening `io::stdin()` would drop those already-consumed bytes. MT-bgzf when
/// `workers > 1`.
pub fn reader_from(
    inner: Box<dyn io::Read + Send>,
    workers: usize,
) -> anyhow::Result<(sam::Header, RecordBufIterBox)> {
    if workers > 1 {
        let mt =
            bgzf::io::MultithreadedReader::with_worker_count(NonZero::new(workers).unwrap(), inner);
        let mut r = bam::io::Reader::from(mt);
        let header = r.read_header()?;
        let hc = header.clone();
        Ok((
            header,
            Box::new(RecordBufIter {
                reader: r,
                header: hc,
            }),
        ))
    } else {
        let mut r = bam::io::Reader::new(inner);
        let header = r.read_header()?;
        let hc = header.clone();
        Ok((
            header,
            Box::new(RecordBufIter {
                reader: r,
                header: hc,
            }),
        ))
    }
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

/// A BAM output sink: single-threaded bgzf (t1) or multithreaded bgzf (t>1).
pub enum BamSink {
    Single(bam::io::Writer<bgzf::io::Writer<Box<dyn Write + Send>>>),
    Multi(bam::io::Writer<bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>>),
}

/// Build the sink (header written), MT-bgzf when `workers > 1`. `level` is the
/// bgzf DEFLATE compression level (0-9 per the CLI, though libdeflate accepts up
/// to 12); it is applied to both the single- and multi-threaded encoders.
pub fn writer(
    output: Option<&Path>,
    header: &sam::Header,
    workers: usize,
    level: u8,
) -> anyhow::Result<BamSink> {
    let clevel = bgzf::io::writer::CompressionLevel::new(level)
        .ok_or_else(|| anyhow::anyhow!("invalid bgzf compression level {level} (expected 0-12)"))?;
    let inner: Box<dyn Write + Send> = match output {
        Some(p) => Box::new(File::create(p)?),
        None => Box::new(io::stdout()),
    };
    if workers > 1 {
        let mt = bgzf::io::multithreaded_writer::Builder::default()
            .set_compression_level(clevel)
            .set_worker_count(NonZero::new(workers).unwrap())
            .build_from_writer(inner);
        let mut w = bam::io::Writer::from(mt);
        w.write_header(header)?;
        Ok(BamSink::Multi(w))
    } else {
        // Build the single-threaded bgzf writer explicitly (rather than
        // `bam::io::Writer::new`, which would force the default level) so `level`
        // takes effect.
        let bgzf_w = bgzf::io::writer::Builder::default()
            .set_compression_level(clevel)
            .build_from_writer(inner);
        let mut w = bam::io::Writer::from(bgzf_w);
        w.write_header(header)?;
        Ok(BamSink::Single(w))
    }
}

impl BamSink {
    pub fn write_record(&mut self, header: &sam::Header, rec: &RecordBuf) -> io::Result<()> {
        match self {
            BamSink::Single(w) => w.write_alignment_record(header, rec),
            BamSink::Multi(w) => w.write_alignment_record(header, rec),
        }
    }

    /// Flush + finalize (bgzf EOF block). Single: `try_finish`; Multi:
    /// `into_inner().finish()` (its `Drop` swallows errors — must be explicit).
    pub fn finish(self) -> anyhow::Result<()> {
        match self {
            BamSink::Single(mut w) => {
                w.try_finish()?;
                Ok(())
            },
            BamSink::Multi(w) => {
                w.into_inner().finish()?;
                Ok(())
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;

    use super::*;

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

    #[test]
    fn mt_writer_roundtrips_through_mt_reader() {
        use noodles_sam::alignment::record::Flags;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mt.bam");

        // Write two unmapped records through a 4-worker MT BamSink.
        let header = sam::Header::default();
        let mut sink = writer(Some(&path), &header, 4, 6).unwrap();
        for name in [b"r1".as_slice(), b"r2".as_slice()] {
            let mut rec = RecordBuf::default();
            *rec.flags_mut() = Flags::UNMAPPED;
            *rec.name_mut() = Some(name.into());
            *rec.sequence_mut() = b"ACGT".to_vec().into();
            *rec.quality_scores_mut() = vec![40u8; 4].into();
            sink.write_record(&header, &rec).unwrap();
        }
        sink.finish().unwrap();

        // Read back through a 4-worker MT reader.
        let (_h, records) = reader(Some(&path), 4).unwrap();
        let names: Vec<Vec<u8>> = records
            .map(|r| r.unwrap().name().map(|n| n.to_vec()).unwrap_or_default())
            .collect();
        assert_eq!(names.len(), 2);
        assert!(names.contains(&b"r1".to_vec()) && names.contains(&b"r2".to_vec()));
    }
}
