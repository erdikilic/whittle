pub mod parse;
pub mod reconstruct;
pub mod serialize;

pub use parse::parse;
pub use reconstruct::reconstruct;
pub use serialize::serialize;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModCode {
    Char(u8),
    Chebi(u32),
}

/// One MM group, e.g. `C+m?,5,12` with its slice of ML bytes.
/// `ml.len() == deltas.len() * codes.len()` (position-major, see plan header).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmGroup {
    pub base: u8,
    pub strand: u8,
    pub codes: Vec<ModCode>,
    pub status: Option<u8>,
    pub deltas: Vec<usize>,
    pub ml: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Mods {
    pub groups: Vec<MmGroup>,
}

pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        b'U' | b'u' => b'A',
        _ => b'N',
    }
}

/// The SEQ base whose occurrences the MM skip-counts index: the fundamental base
/// for `+`, its complement for `-` (the mods sit on the opposite strand). Slicing
/// only needs to count the SAME base the encoder counted — the htslib oracle
/// confirms this matches real data.
pub fn counting_base(base: u8, strand: u8) -> u8 {
    if strand == b'-' {
        complement(base)
    } else {
        base.to_ascii_uppercase()
    }
}
