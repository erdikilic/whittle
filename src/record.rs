/// Format-neutral read carrier. `qual` holds raw Phred scores (0-based), not ASCII.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadRecord {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    /// Raw Phred scores. Values must stay <= 222 for lossless FASTQ round-trips
    /// (ASCII emission adds 33; printable ASCII tops out at 126).
    pub qual: Vec<u8>,
}
