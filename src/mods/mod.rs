//! Base-modification (`MM`/`ML`) parsing, windowing and serialization.

pub mod parse;
pub mod reconstruct;
pub mod serialize;

pub use parse::{MalformedMm, expected_ml_len, parse, parse_checked};
pub use reconstruct::reconstruct;
pub use serialize::serialize;

/// One modification code in an `MM` group: a single-letter code or a numeric
/// ChEBI id.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModCode {
    /// Single-letter code, e.g. `m`.
    Char(u8),
    /// Numeric ChEBI identifier.
    Chebi(u32),
}

/// One `MM` group, e.g. `C+m?,5,12`, with its slice of `ML` bytes.
/// `ml.len() == deltas.len() * codes.len()`: ML bytes are position-major, one
/// byte per code for each modified position.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmGroup {
    /// Fundamental base as written in `MM` (`C`, `A`, `N`, ...).
    pub base: u8,
    /// Strand character, `+` or `-`.
    pub strand: u8,
    /// Modification codes for this group.
    pub codes: Vec<ModCode>,
    /// Optional status flag, `?` or `.`.
    pub status: Option<u8>,
    /// Skip-counts between successive modified occurrences of the counting base.
    pub deltas: Vec<usize>,
    /// `ML` probability bytes for this group, position-major.
    pub ml: Vec<u8>,
}

/// The parsed `MM`/`ML` content of one record.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Mods {
    /// Groups in source order.
    pub groups: Vec<MmGroup>,
}

/// Returns the SEQ base whose occurrences an `MM` group's skip-counts index:
/// the fundamental base as written in `MM`, with `U` folded to `T`.
///
/// The strand is not applied. htslib counts the literal base for
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

/// Reports whether `seq_base` is a base that a group with counting base `cbase`
/// indexes. `N` in `MM` means any base and matches every position in SEQ.
pub fn counts(seq_base: u8, cbase: u8) -> bool {
    cbase == b'N' || seq_base.to_ascii_uppercase() == cbase
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pinned against htslib: `MM:Z:G-m,0,2;` decodes to G positions, not C.
    /// Complementing the base would relocate every reverse-strand call.
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
