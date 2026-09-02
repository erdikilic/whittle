/// Format-neutral read carrier. `qual` holds raw Phred scores (0-based), not ASCII.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadRecord {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    /// Raw Phred scores. The FASTQ reader rejects quality bytes outside ASCII
    /// 33..=126, so FASTQ-sourced values lie in 0..=93; ASCII emission adds 33
    /// with saturation.
    pub qual: Vec<u8>,
}
