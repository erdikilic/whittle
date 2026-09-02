//! A `Read` adapter that counts the bytes pulled from its inner reader into `Counters::bytes_read`.

use std::io::Read;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use crate::workflow::Counters;

/// A `Read` wrapper that adds the byte count of each `read` call to
/// `counters.bytes_read`, which drives the progress bar.
pub struct CountingReader<R> {
    inner: R,
    counters: Arc<Counters>,
}

impl<R: Read> CountingReader<R> {
    /// Wraps `inner` so that every `read` adds its byte count to `counters.bytes_read`.
    pub fn new(inner: R, counters: Arc<Counters>) -> Self {
        Self { inner, counters }
    }
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.counters
            .bytes_read
            .fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::sync::atomic::Ordering;

    use super::*;

    #[test]
    fn counts_all_bytes_read() {
        let data = b"hello world, this is 30 bytes!";
        let counters = Arc::new(Counters::default());
        let mut r = CountingReader::new(Cursor::new(&data[..]), counters.clone());
        let mut out = Vec::new();
        r.read_to_end(&mut out).unwrap();
        assert_eq!(out, data);
        assert_eq!(
            counters.bytes_read.load(Ordering::Relaxed),
            data.len() as u64
        );
    }
}
