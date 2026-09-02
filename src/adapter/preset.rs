//! The built-in adapter presets, built from the ONT catalog.

use super::ont_catalog::{CATALOG, Entry};
use super::{Adapter, End};

/// Builds the searchable adapter set from catalog entries.
///
/// Identical sequences collapse to one entry, keeping the first name. When a
/// duplicate carries a different end tag the survivor is upgraded to
/// `End::Both`, so a sequence that appears at either end (for example a primer
/// that is also a barcode flank) stays searchable at both.
fn build(entries: &[Entry]) -> Vec<Adapter> {
    let mut out: Vec<Adapter> = Vec::with_capacity(entries.len());
    let mut idx: std::collections::HashMap<&[u8], usize> = std::collections::HashMap::new();
    for &(name, end, seq) in entries {
        match idx.get(seq) {
            Some(&i) => {
                if out[i].end != end {
                    out[i].end = End::Both;
                }
            },
            None => {
                idx.insert(seq, out.len());
                out.push(Adapter {
                    name: name.to_string(),
                    seq: seq.to_vec(),
                    end,
                });
            },
        }
    }
    out
}

/// Returns the built-in ONT catalog, deduplicated and ready to search.
pub fn preset_ont() -> Vec<Adapter> {
    build(CATALOG)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Identical sequences collapse to the first name, and differing end tags
    /// merge to `End::Both`.
    #[test]
    fn duplicate_sequences_collapse_and_merge_ends() {
        let entries: &[Entry] = &[
            ("A", End::Five, b"ACGTACGTACGT"),
            ("B", End::Three, b"TTTTGGGGCCCC"),
            ("Dup", End::Three, b"ACGTACGTACGT"),
        ];
        let v = build(entries);
        assert_eq!(v.len(), 2, "Duplicate sequence collapsed");
        assert_eq!(v[0].name, "A", "First name kept");
        assert_eq!(
            v[0].end,
            End::Both,
            "5' + 3' of the same sequence merges to Both"
        );
        assert_eq!(v[1].end, End::Three, "Unique entry keeps its own end");
    }

    /// A zero-length pattern matches everywhere, and a byte outside the
    /// nucleotide alphabet would panic the searcher. The catalog is a
    /// compile-time literal with no parse-time validation, so this test
    /// enforces both invariants.
    ///
    /// Ambiguity codes are permitted (the searcher handles them), though every
    /// catalog entry is plain ACGT.
    #[test]
    fn entries_are_valid_nucleotide_sequences() {
        for &(name, _, seq) in CATALOG {
            assert!(
                !seq.is_empty(),
                "Catalog entry {name} has an empty sequence"
            );
            for &b in seq {
                assert!(
                    crate::adapter::search::iupac_degeneracy(b).is_some(),
                    "Catalog entry {name} has a non-nucleotide byte {:?}: {}",
                    b as char,
                    String::from_utf8_lossy(seq)
                );
                assert_eq!(
                    b,
                    b.to_ascii_uppercase(),
                    "Catalog entry {name} must be uppercase"
                );
            }
        }
    }

    /// A pattern below `MIN_PATTERN_LEN` is skipped by every search loop, so a
    /// catalog entry that short would be counted as configured and never act.
    #[test]
    fn entries_meet_the_minimum_pattern_length() {
        for &(name, _, seq) in CATALOG {
            assert!(
                seq.len() >= crate::adapter::MIN_PATTERN_LEN,
                "{name} is {} bp, below the {} bp searchable minimum",
                seq.len(),
                crate::adapter::MIN_PATTERN_LEN
            );
        }
    }

    /// No two catalog entries share a display name.
    #[test]
    fn entry_names_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for &(name, _, _) in CATALOG {
            assert!(seen.insert(name), "Duplicate catalog name {name}");
        }
    }

    /// `PCR1_front` (5') and `LWB_rear1` (3') share one sequence, which is
    /// searchable at both ends.
    #[test]
    fn preset_merges_pcr1_lwb_shared_sequence_to_both() {
        let v = preset_ont();
        let e = v
            .iter()
            .find(|a| a.seq == b"ACTTGCCTGTCGCTCTATCTTC")
            .expect("Shared sequence is present");
        assert_eq!(
            e.end,
            End::Both,
            "Shared 5'/3' sequence must be searchable at both ends"
        );
    }

    /// The catalog holds 121 entries and 120 after the one exact duplicate
    /// (`PCR1_front` and `LWB_rear1`) collapses, including all 96 barcodes.
    #[test]
    fn preset_has_the_expected_shape() {
        assert_eq!(CATALOG.len(), 121, "Catalog entries before dedup");
        let v = preset_ont();
        assert_eq!(v.len(), 120, "Catalog entry count after dedup");
        assert!(v.iter().any(|a| a.name == "LSK114_front"));
        assert_eq!(v.iter().filter(|a| a.name.starts_with("BC")).count(), 96);
    }
}
