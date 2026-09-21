//! BGZF block encoding over libdeflate: one compressor per thread reused across
//! blocks, the block framing, and a streaming writer for sequential output.

use std::cell::RefCell;
use std::io::{self, Write};

use libdeflater::{CompressionLvl, Compressor, Crc};

/// Largest uncompressed payload of one block, the htslib block size. A block
/// of this size fits the 16-bit `BSIZE` field even when stored uncompressed.
pub(crate) const BLOCK_PAYLOAD: usize = 0xff00;

/// The gzip header with the `BC` extra subfield.
const HEADER_SIZE: usize = 18;

/// CRC-32 and uncompressed size.
const TRAILER_SIZE: usize = 8;

/// The empty block every BGZF stream ends with.
pub(crate) const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Returns the DEFLATE level for `level`, or an error naming the accepted range.
pub(crate) fn compression_level(level: u8) -> anyhow::Result<CompressionLvl> {
    if level > 9 {
        anyhow::bail!("invalid BGZF compression level {level} (expected 0-9)");
    }
    CompressionLvl::new(i32::from(level))
        .map_err(|e| anyhow::anyhow!("invalid BGZF compression level {level}: {e:?}"))
}

thread_local! {
    /// The calling thread's compressor and the level it was built for.
    static COMPRESSOR: RefCell<Option<(CompressionLvl, Compressor)>> = const { RefCell::new(None) };
}

/// Runs `f` with the calling thread's compressor at `level`, building one on
/// first use or when the level changes.
fn with_compressor<T>(level: CompressionLvl, f: impl FnOnce(&mut Compressor) -> T) -> T {
    COMPRESSOR.with(|slot| {
        let mut slot = slot.borrow_mut();
        if !matches!(slot.as_ref(), Some((l, _)) if *l == level) {
            *slot = Some((level, Compressor::new(level)));
        }
        let (_, compressor) = slot.as_mut().expect("The compressor was set above");
        f(compressor)
    })
}

/// Appends `data` to `out` as complete blocks at `level`. Empty `data`
/// appends nothing.
pub(crate) fn encode(level: u8, data: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    let level = compression_level(level).map_err(io::Error::other)?;
    with_compressor(level, |compressor| {
        for payload in data.chunks(BLOCK_PAYLOAD) {
            encode_block(compressor, payload, out)?;
        }
        Ok(())
    })
}

/// Appends one block holding `payload` to `out`.
fn encode_block(compressor: &mut Compressor, payload: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    let start = out.len();
    let data_start = start + HEADER_SIZE;
    let bound = compressor.deflate_compress_bound(payload.len());
    out.resize(data_start + bound, 0);
    let compressed = compressor
        .deflate_compress(payload, &mut out[data_start..])
        .map_err(|e| io::Error::other(format!("BGZF block compression failed: {e:?}")))?;
    out.truncate(data_start + compressed);
    let block_size = HEADER_SIZE + compressed + TRAILER_SIZE;
    let bsize = u16::try_from(block_size - 1)
        .map_err(|_| io::Error::other("BGZF block exceeds 64 KiB"))?
        .to_le_bytes();
    out[start..data_start].copy_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02,
        0x00, bsize[0], bsize[1],
    ]);
    let mut crc = Crc::new();
    crc.update(payload);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    Ok(())
}

/// A streaming BGZF writer for a sequential run: bytes stage up to one block
/// payload and are compressed on the writing thread.
pub struct Writer<W: Write> {
    inner: W,
    level: u8,
    staged: Vec<u8>,
    block: Vec<u8>,
}

impl<W: Write> Writer<W> {
    /// Wraps `inner`, compressing at `level`.
    pub(crate) fn new(inner: W, level: u8) -> Self {
        Self {
            inner,
            level,
            staged: Vec::with_capacity(BLOCK_PAYLOAD),
            block: Vec::new(),
        }
    }

    /// Compresses the staged bytes, if any, into one block and writes it.
    fn write_block(&mut self) -> io::Result<()> {
        if self.staged.is_empty() {
            return Ok(());
        }
        self.block.clear();
        encode(self.level, &self.staged, &mut self.block)?;
        self.inner.write_all(&self.block)?;
        self.staged.clear();
        Ok(())
    }

    /// Writes the staged block and returns the destination without an EOF
    /// block, for output that continues with blocks compressed elsewhere.
    pub(crate) fn into_inner(mut self) -> io::Result<W> {
        self.write_block()?;
        Ok(self.inner)
    }

    /// Writes the staged block and the EOF block, flushes the destination and
    /// returns it.
    pub(crate) fn finish(mut self) -> io::Result<W> {
        self.write_block()?;
        self.inner.write_all(&EOF_BLOCK)?;
        self.inner.flush()?;
        Ok(self.inner)
    }
}

impl<W: Write> Write for Writer<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let room = BLOCK_PAYLOAD - self.staged.len();
        let n = buf.len().min(room);
        self.staged.extend_from_slice(&buf[..n]);
        if self.staged.len() == BLOCK_PAYLOAD {
            self.write_block()?;
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.write_block()?;
        self.inner.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    fn inflate(bgzf: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        noodles_bgzf::io::Reader::new(bgzf)
            .read_to_end(&mut out)
            .unwrap();
        out
    }

    fn sample(len: usize) -> Vec<u8> {
        (0..len).map(|i| b"ACGT"[(i * 7 + i / 13) % 4]).collect()
    }

    #[test]
    fn eof_block_matches_the_reference_stream_end() {
        let reference = noodles_bgzf::io::Writer::new(Vec::new()).finish().unwrap();
        assert_eq!(reference, EOF_BLOCK);
    }

    #[test]
    fn encode_splits_at_the_block_payload_and_roundtrips() {
        let data = sample(3 * BLOCK_PAYLOAD + 1);
        let mut out = Vec::new();
        encode(4, &data, &mut out).unwrap();
        out.extend_from_slice(&EOF_BLOCK);
        assert_eq!(inflate(&out), data);
        let mut blocks = 0;
        let mut pos = 0;
        while pos < out.len() {
            assert_eq!(&out[pos..pos + 4], &[0x1f, 0x8b, 0x08, 0x04]);
            let bsize = u16::from_le_bytes([out[pos + 16], out[pos + 17]]) as usize + 1;
            pos += bsize;
            blocks += 1;
        }
        assert_eq!(blocks, 5, "four data blocks and the EOF block");
    }

    #[test]
    fn encode_appends_nothing_for_empty_data_and_stores_at_level_zero() {
        let mut out = Vec::new();
        encode(4, &[], &mut out).unwrap();
        assert!(out.is_empty());
        let data = sample(BLOCK_PAYLOAD);
        encode(0, &data, &mut out).unwrap();
        out.extend_from_slice(&EOF_BLOCK);
        assert_eq!(inflate(&out), data);
    }

    #[test]
    fn compression_level_rejects_values_above_nine() {
        assert!(compression_level(9).is_ok());
        assert!(compression_level(10).is_err());
    }

    #[test]
    fn writer_stages_to_full_blocks_and_finishes_with_eof() {
        let data = sample(2 * BLOCK_PAYLOAD + 100);
        let mut w = Writer::new(Vec::new(), 1);
        for chunk in data.chunks(1000) {
            w.write_all(chunk).unwrap();
        }
        let out = w.finish().unwrap();
        assert!(out.ends_with(&EOF_BLOCK));
        assert_eq!(inflate(&out), data);
    }

    #[test]
    fn writer_into_inner_omits_the_eof_block() {
        let mut w = Writer::new(Vec::new(), 1);
        w.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        let out = w.into_inner().unwrap();
        assert!(!out.ends_with(&EOF_BLOCK));
        let mut with_eof = out.clone();
        with_eof.extend_from_slice(&EOF_BLOCK);
        assert_eq!(inflate(&with_eof), b"@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn output_is_detected_as_bgzf() {
        let mut out = Vec::new();
        encode(6, b"@r1\nACGT\n+\nIIII\n", &mut out).unwrap();
        assert!(crate::io::is_bgzf(&out));
    }
}
