pub mod parse;
pub mod reconstruct;
pub mod serialize;

pub use parse::{MalformedMm, expected_ml_len, parse, parse_checked};
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

/// The SEQ base whose occurrences an MM group's skip-counts index: the
/// fundamental base as written in MM, with `U` folded to `T`.
///
/// The strand is deliberately NOT applied. htslib counts the literal base for
/// both `+` and `-` groups: in `sam_mods.c` the parsed base goes into
/// `canonical[]` while the strand is stored separately, and the only
/// complementing (`seqi_rc`) is gated on the record's `BAM_FREVERSE` flag, which
/// an unaligned BAM never carries. `U` is folded because BAM's 4-bit SEQ
/// encoding has no `U`, so an RNA group would otherwise match nothing.
pub fn counting_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'U' => b'T',
        b => b,
    }
}

/// Whether `seq_base` is one of the bases a group with counting base `cbase`
/// indexes. `N` in MM means "any base", matching every position in SEQ.
pub fn counts(seq_base: u8, cbase: u8) -> bool {
    cbase == b'N' || seq_base.to_ascii_uppercase() == cbase
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pinned against htslib: `MM:Z:G-m,0,2;` decodes to G positions, not C.
    /// Complementing here silently relocated every reverse-strand call.
    #[test]
    fn minus_strand_counts_the_literal_base() {
        assert_eq!(counting_base(b'G'), b'G');
        assert_eq!(counting_base(b'C'), b'C');
    }

    /// BAM SEQ has no `U`, so an RNA group must count `T` or it matches nothing
    /// and the whole group is dropped.
    #[test]
    fn uracil_folds_to_thymine() {
        assert_eq!(counting_base(b'U'), b'T');
        assert_eq!(counting_base(b'u'), b'T');
    }

    #[test]
    fn n_counts_every_base() {
        for b in [b'A', b'C', b'G', b'T', b'N'] {
            assert!(counts(b, b'N'), "N must count {}", b as char);
        }
        assert!(counts(b'C', b'C'));
        assert!(!counts(b'G', b'C'));
        // Lowercase SEQ bytes still count.
        assert!(counts(b'c', b'C'));
    }
}
