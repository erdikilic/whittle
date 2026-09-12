//! The built-in adapter presets: one per sequencing kit family, built from the
//! catalog.

use std::fmt;

use super::catalog::{CATALOG, Entry};
use super::{Adapter, Role};

/// A sequencing kit family the catalog knows. Each variant selects the catalog
/// entries a library made with that kit can contain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Kit {
    /// Ligation sequencing kit V14 (SQK-LSK114), including its cDNA and 10X
    /// primer variants.
    Lsk114,
    /// Rapid sequencing kit V14 (SQK-RAD114) and the ultra-long kit (SQK-ULK114).
    Rad114,
    /// Rapid barcoding kit V14 (SQK-RBK114, 24 or 96 barcodes).
    Rbk114,
    /// Native barcoding kit V14 (SQK-NBD114, 24 or 96 barcodes).
    Nbd114,
    /// cDNA-PCR sequencing and barcoding kits V14 (SQK-PCS114, SQK-PCB114).
    Pcb114,
    /// Rapid PCR barcoding kit V14 (SQK-RPB114).
    Rpb114,
    /// Microbial amplicon barcoding kit V14 (SQK-MAB114): 16S and ITS primers
    /// and its own 24 barcodes.
    Mab114,
    /// Direct RNA sequencing kit (SQK-RNA004).
    Rna004,
    /// PacBio SMRTbell libraries.
    PacBio,
}

impl Kit {
    /// Every kit, in display order.
    pub const ALL: &[Kit] = &[
        Kit::Lsk114,
        Kit::Rad114,
        Kit::Rbk114,
        Kit::Nbd114,
        Kit::Pcb114,
        Kit::Rpb114,
        Kit::Mab114,
        Kit::Rna004,
        Kit::PacBio,
    ];

    /// The Oxford Nanopore kits: everything but PacBio.
    pub const ONT: &[Kit] = &[
        Kit::Lsk114,
        Kit::Rad114,
        Kit::Rbk114,
        Kit::Nbd114,
        Kit::Pcb114,
        Kit::Rpb114,
        Kit::Mab114,
        Kit::Rna004,
    ];

    /// Whether every library of this kit is a PCR amplicon, so that a primer
    /// inside a read marks a chimeric junction rather than part of the molecule.
    pub fn amplicon(self) -> bool {
        matches!(self, Kit::Mab114)
    }
}

/// The preset tokens `--adapter-preset` accepts, with their expansions.
pub const PRESET_TOKENS: &[(&str, &[Kit])] = &[
    ("lsk114", &[Kit::Lsk114]),
    ("rad114", &[Kit::Rad114]),
    ("ulk114", &[Kit::Rad114]),
    ("rbk114", &[Kit::Rbk114]),
    ("nbd114", &[Kit::Nbd114]),
    ("pcb114", &[Kit::Pcb114]),
    ("pcs114", &[Kit::Pcb114]),
    ("rpb114", &[Kit::Rpb114]),
    ("mab114", &[Kit::Mab114]),
    ("rna004", &[Kit::Rna004]),
    ("pacbio", &[Kit::PacBio]),
    ("ont", Kit::ONT),
    ("all", Kit::ALL),
];

/// A preset token that names no kit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnknownPreset(pub String);

impl fmt::Display for UnknownPreset {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let known: Vec<&str> = PRESET_TOKENS.iter().map(|(token, _)| *token).collect();
        write!(
            f,
            "unknown adapter preset {:?}; expected one of {}",
            self.0,
            known.join(", ")
        )
    }
}

/// Parses a comma-separated preset list into the kits it names, deduplicated
/// and in catalog order. Tokens are case-insensitive; `none` and empty tokens
/// name nothing.
pub fn parse_presets(spec: &str) -> Result<Vec<Kit>, UnknownPreset> {
    let mut kits: Vec<Kit> = Vec::new();
    for token in spec.split(',') {
        let token = token.trim();
        if token.is_empty() || token.eq_ignore_ascii_case("none") {
            continue;
        }
        let expansion = PRESET_TOKENS
            .iter()
            .find(|(name, _)| name.eq_ignore_ascii_case(token))
            .map(|(_, kits)| *kits)
            .ok_or_else(|| UnknownPreset(token.to_string()))?;
        kits.extend_from_slice(expansion);
    }
    kits.sort();
    kits.dedup();
    Ok(kits)
}

/// Builds the searchable adapter set from catalog entries.
///
/// Identical sequences collapse to one entry, keeping the first name. When a
/// duplicate carries a different role the survivor keeps the more permissive
/// one, so a sequence that is an adapter in one kit and a flank in another
/// still splits.
fn build(entries: &[Entry]) -> Vec<Adapter> {
    let mut out: Vec<Adapter> = Vec::with_capacity(entries.len());
    let mut idx: std::collections::HashMap<&[u8], usize> = std::collections::HashMap::new();
    for &(name, role, _, seq) in entries {
        match idx.get(seq) {
            Some(&i) => {
                if role.splits() {
                    out[i].role = role;
                }
            },
            None => {
                idx.insert(seq, out.len());
                out.push(Adapter {
                    name: name.to_string(),
                    seq: seq.to_vec(),
                    role,
                });
            },
        }
    }
    out
}

/// Returns the catalog entries of `kits`, deduplicated and ready to search.
///
/// When every selected kit is an amplicon kit, its primers take the adapter
/// role and split reads at interior hits: an amplicon library has no molecule
/// with a primer inside it, so an interior primer is a chimeric junction. In a
/// selection that includes a genomic kit the primers keep their terminal-only
/// role, since the primer site of a marker gene is part of a genomic read.
pub fn preset(kits: &[Kit]) -> Vec<Adapter> {
    let amplicon_only = !kits.is_empty() && kits.iter().all(|kit| kit.amplicon());
    let entries: Vec<Entry> = CATALOG
        .iter()
        .filter(|(_, _, members, _)| members.iter().any(|kit| kits.contains(kit)))
        .map(|&(name, role, members, seq)| {
            let role = if amplicon_only && role == Role::Primer {
                Role::Adapter
            } else {
                role
            };
            (name, role, members, seq)
        })
        .collect();
    build(&entries)
}

/// Returns the Oxford Nanopore catalog, deduplicated.
pub fn preset_ont() -> Vec<Adapter> {
    preset(Kit::ONT)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Identical sequences collapse to the first name, and a splitting role wins
    /// over a terminal-only one.
    #[test]
    fn duplicate_sequences_collapse_and_keep_the_splitting_role() {
        let entries: &[Entry] = &[
            ("A", Role::Primer, &[Kit::Lsk114], b"ACGTACGTACGT"),
            ("B", Role::Barcode, &[Kit::Lsk114], b"TTTTGGGGCCCC"),
            ("Dup", Role::Adapter, &[Kit::Lsk114], b"ACGTACGTACGT"),
        ];
        let v = build(entries);
        assert_eq!(v.len(), 2, "Duplicate sequence collapsed");
        assert_eq!(v[0].name, "A", "First name kept");
        assert_eq!(v[0].role, Role::Adapter, "Splitting role wins");
        assert_eq!(v[1].role, Role::Barcode, "Unique entry keeps its role");
    }

    /// A zero-length pattern matches everywhere, and a byte outside the
    /// nucleotide alphabet would panic the searcher. The catalog is a
    /// compile-time literal with no parse-time validation, so this test
    /// enforces both invariants. Ambiguity codes are permitted.
    #[test]
    fn entries_are_valid_nucleotide_sequences() {
        for &(name, _, _, seq) in CATALOG {
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
        for &(name, _, _, seq) in CATALOG {
            assert!(
                seq.len() >= crate::adapter::MIN_PATTERN_LEN,
                "{name} is {} bp, below the {} bp searchable minimum",
                seq.len(),
                crate::adapter::MIN_PATTERN_LEN
            );
        }
    }

    /// No two catalog entries share a display name, and every entry belongs to
    /// at least one kit.
    #[test]
    fn entry_names_are_unique_and_kits_are_nonempty() {
        let mut seen = std::collections::HashSet::new();
        for &(name, _, kits, _) in CATALOG {
            assert!(seen.insert(name), "Duplicate catalog name {name}");
            assert!(!kits.is_empty(), "{name} belongs to no kit");
        }
    }

    /// Every kit selects its ligation or rapid adapter, and the amplicon kit
    /// carries its degenerate primers and barcodes.
    #[test]
    fn presets_select_their_kit_sequences() {
        let lsk = preset(&[Kit::Lsk114]);
        assert!(lsk.iter().any(|a| a.name == "LSK114_front"));
        assert!(lsk.iter().all(|a| !a.name.starts_with("BC")));
        assert!(lsk.iter().all(|a| a.name != "RAD"));

        let mab = preset(&[Kit::Mab114]);
        assert!(
            mab.iter()
                .any(|a| a.name == "16S_27F" && a.role == Role::Adapter)
        );
        assert_eq!(mab.iter().filter(|a| a.name.starts_with("TP")).count(), 24);
        assert!(
            mab.iter()
                .any(|a| a.name == "RAD" && a.role == Role::Adapter)
        );
        let mixed = preset(&[Kit::Mab114, Kit::Lsk114]);
        assert!(
            mixed
                .iter()
                .any(|a| a.name == "16S_27F" && a.role == Role::Primer)
        );

        let pb = preset(&[Kit::PacBio]);
        assert_eq!(pb.len(), 2);
        assert!(pb.iter().any(|a| a.name == "SMRTbell"));
    }

    /// The whole catalog holds every barcode, and the ONT union excludes PacBio.
    #[test]
    fn unions_have_the_expected_shape() {
        let all = preset(Kit::ALL);
        assert_eq!(all.iter().filter(|a| a.name.starts_with("BC")).count(), 96);
        assert_eq!(all.iter().filter(|a| a.name.starts_with("TP")).count(), 24);
        assert_eq!(CATALOG.len(), all.len(), "Catalog sequences are unique");
        let ont = preset_ont();
        assert!(ont.iter().all(|a| a.name != "SMRTbell"));
        assert_eq!(ont.len(), all.len() - 2);
    }

    /// Preset lists parse case-insensitively, expand aliases and unions, and
    /// reject unknown tokens.
    #[test]
    fn preset_lists_parse_and_reject_unknown_tokens() {
        assert_eq!(
            parse_presets("LSK114,mab114").unwrap(),
            vec![Kit::Lsk114, Kit::Mab114]
        );
        assert_eq!(parse_presets("ulk114").unwrap(), vec![Kit::Rad114]);
        assert_eq!(parse_presets("ont").unwrap(), Kit::ONT.to_vec());
        assert_eq!(parse_presets("none").unwrap(), Vec::<Kit>::new());
        assert_eq!(
            parse_presets("lsk114,16s").unwrap_err(),
            UnknownPreset("16s".into())
        );
    }
}
