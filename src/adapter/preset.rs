use super::ont_catalog::{CATALOG, Entry};
use super::{Adapter, End};

/// Build the searchable adapter set from catalog entries.
///
/// Identical sequences collapse to one entry, keeping the first name. When a
/// duplicate carries a different end tag the survivor is upgraded to
/// `End::Both`, so a sequence that legitimately appears at either end (a primer
/// that is also a barcode flank, say) stays searchable at both.
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

/// The built-in ONT catalog, deduplicated and ready to search.
pub fn preset_ont() -> Vec<Adapter> {
    build(CATALOG)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn duplicate_sequences_collapse_and_merge_ends() {
        let entries: &[Entry] = &[
            ("A", End::Five, b"ACGTACGTACGT"),
            ("B", End::Three, b"TTTTGGGGCCCC"),
            ("Dup", End::Three, b"ACGTACGTACGT"),
        ];
        let v = build(entries);
        assert_eq!(v.len(), 2, "duplicate sequence collapsed");
        assert_eq!(v[0].name, "A", "first name kept");
        assert_eq!(
            v[0].end,
            End::Both,
            "5' + 3' of the same sequence merges to Both"
        );
        assert_eq!(v[1].end, End::Three, "unique entry keeps its own end");
    }

    /// The search engine only handles uppercase ACGT, and a zero-length pattern
    /// would match everywhere. The old TSV filtered both out at parse time; as
    /// literals they cannot be filtered, so this guards future catalog edits.
    #[test]
    fn entries_are_uppercase_acgt() {
        for &(name, _, seq) in CATALOG {
            assert!(!seq.is_empty(), "{name} has an empty sequence");
            assert!(
                seq.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')),
                "{name} has a non-ACGT base: {}",
                String::from_utf8_lossy(seq)
            );
        }
    }

    #[test]
    fn entry_names_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for &(name, _, _) in CATALOG {
            assert!(seen.insert(name), "duplicate catalog name {name}");
        }
    }

    #[test]
    fn preset_merges_pcr1_lwb_shared_sequence_to_both() {
        let v = preset_ont();
        // PCR1_front (5') and LWB_rear1 (3') share ACTTGCCTGTCGCTCTATCTTC.
        let e = v
            .iter()
            .find(|a| a.seq == b"ACTTGCCTGTCGCTCTATCTTC")
            .expect("shared seq present");
        assert_eq!(
            e.end,
            End::Both,
            "shared 5'/3' sequence must be searchable at both ends"
        );
    }

    #[test]
    fn preset_has_the_expected_shape() {
        assert_eq!(CATALOG.len(), 125, "catalog entries before dedup");
        let v = preset_ont();
        // 96 barcodes plus adapters/primers/flanks, minus the one exact-duplicate
        // sequence (PCR1_front == LWB_rear1). Expect 124 after dedup.
        assert_eq!(v.len(), 124, "catalog entry count after dedup");
        assert!(v.iter().any(|a| a.name == "LSK114_front"));
        assert_eq!(v.iter().filter(|a| a.name.starts_with("BC")).count(), 96);
    }
}
