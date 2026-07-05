mod bam;
mod fastq;

pub use bam::{reconstruct_mods, reconstruct_record, run_bam, run_bam_to_fastq};
pub use fastq::{run_fastq, run_fastq_seq};

#[derive(Debug, Default, Clone, Copy)]
pub struct Stats {
    pub input_reads: u64,
    pub output_reads: u64,
    /// Reads carrying a known per-base kinetics tag (ip/pw/…) whose array length
    /// did not match the sequence length — malformed and left untouched. Surfaced
    /// as a run-level advisory; not an error.
    pub malformed_tag_reads: u64,
}
