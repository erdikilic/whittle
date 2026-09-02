//! The format-neutral read record shared by the FASTQ and BAM workflows.

/// A format-neutral read. `qual` holds raw Phred scores (0-based), not ASCII.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadRecord {
    /// Read identifier, including any description or tab-delimited tags that
    /// followed it in the source header.
    pub name: Vec<u8>,
    /// Nucleotide sequence as ASCII bytes.
    pub seq: Vec<u8>,
    /// Raw Phred scores. The FASTQ reader rejects quality bytes outside ASCII
    /// 33..=126, so FASTQ-sourced values lie in 0..=93; ASCII emission adds 33
    /// with saturation.
    pub qual: Vec<u8>,
}
