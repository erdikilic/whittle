use super::{Adapter, End};

/// The built-in ONT catalog, embedded at compile time (lives next to this file).
const CATALOG_TSV: &str = include_str!("ont_catalog.tsv");

/// Parse the catalog TSV: skip blank and `#` lines, take columns
/// `id, category, end, sequence, ...`. Column 3 (`end`) is `5`/`3`/`both`.
/// Identical sequences are deduplicated (first name kept). Non-ACGT rows are
/// skipped defensively.
pub fn parse_catalog(text: &str) -> Vec<Adapter> {
    let mut out: Vec<Adapter> = Vec::new();
    let mut seen: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            continue;
        }
        let name = cols[0].to_string();
        let end = match cols[2] {
            "5" => End::Five,
            "3" => End::Three,
            _ => End::Both,
        };
        let seq = cols[3].as_bytes().to_vec();
        if seq.is_empty() || !seq.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        if seen.insert(seq.clone()) {
            out.push(Adapter { name, seq, end });
        }
    }
    out
}

/// The parsed built-in ONT catalog.
pub fn preset_ont() -> Vec<Adapter> {
    parse_catalog(CATALOG_TSV)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_rows_skips_comments_and_dedups() {
        let tsv = "# header comment\n\
                   id\tcategory\tend\tsequence\tkits\tsource\n\
                   A\tligation-adapter\t5\tACGTACGTACGT\tk\ts\n\
                   B\tflank\t3\tTTTTGGGGCCCC\tk\ts\n\
                   Dup\tbarcode\tboth\tACGTACGTACGT\tk\ts\n";
        // The literal "id...sequence..." header row has a non-ACGT seq column
        // ("sequence") so it is dropped by the ACGT filter; the comment is skipped;
        // "Dup" duplicates A's sequence.
        let v = parse_catalog(tsv);
        assert_eq!(v.len(), 2, "header + duplicate must be dropped");
        assert_eq!(v[0].name, "A");
        assert_eq!(v[0].end, End::Five);
        assert_eq!(v[1].end, End::Three);
    }

    #[test]
    fn preset_has_the_expected_shape() {
        let v = preset_ont();
        // 96 barcodes + adapters/primers/flanks, minus the one exact-duplicate
        // sequence (PCR1_front == LWB flank). Expect 124 after dedup.
        assert_eq!(v.len(), 124, "catalog entry count after dedup");
        assert!(v.iter().any(|a| a.name == "LSK114_front"));
        assert_eq!(v.iter().filter(|a| a.name.starts_with("BC")).count(), 96);
    }
}
