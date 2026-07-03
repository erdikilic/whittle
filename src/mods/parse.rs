use super::{MmGroup, ModCode, Mods};

/// Parse a raw MM:Z string plus its ML:B,C array into groups. Malformed tails are
/// tolerated (best-effort): parsing a group stops at the first unexpected byte.
pub fn parse(mm: &[u8], ml: &[u8]) -> Mods {
    let mut groups = Vec::new();
    let mut ml_pos = 0usize;

    for token in mm.split(|&b| b == b';') {
        if token.len() < 2 {
            continue; // empty (trailing ';') or malformed
        }
        let base = token[0];
        let strand = token[1];
        let mut i = 2;

        // Codes: either a run of letters (each one code) or a numeric ChEBI id.
        let mut codes = Vec::new();
        if i < token.len() && token[i].is_ascii_digit() {
            let mut id = 0u32;
            while i < token.len() && token[i].is_ascii_digit() {
                id = id * 10 + (token[i] - b'0') as u32;
                i += 1;
            }
            codes.push(ModCode::Chebi(id));
        } else {
            while i < token.len() && token[i].is_ascii_alphabetic() {
                codes.push(ModCode::Char(token[i]));
                i += 1;
            }
        }

        // Optional status flag.
        let mut status = None;
        if i < token.len() && (token[i] == b'.' || token[i] == b'?') {
            status = Some(token[i]);
            i += 1;
        }

        // Skip-count deltas: (',' number)*
        let mut deltas = Vec::new();
        while i < token.len() {
            if token[i] != b',' {
                break;
            }
            i += 1;
            let mut n = 0usize;
            let mut saw = false;
            while i < token.len() && token[i].is_ascii_digit() {
                n = n * 10 + (token[i] - b'0') as usize;
                i += 1;
                saw = true;
            }
            if saw {
                deltas.push(n);
            }
        }

        // Claim this group's ML bytes: positions * codes, position-major.
        let want = deltas.len() * codes.len().max(1);
        let end = (ml_pos + want).min(ml.len());
        let group_ml = ml[ml_pos..end].to_vec();
        ml_pos = end;

        groups.push(MmGroup { base, strand, codes, status, deltas, ml: group_ml });
    }

    Mods { groups }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::ModCode;

    #[test]
    fn single_group_single_code() {
        let m = parse(b"C+m?,5,12,0;", &[200, 10, 128]);
        assert_eq!(m.groups.len(), 1);
        let g = &m.groups[0];
        assert_eq!((g.base, g.strand), (b'C', b'+'));
        assert_eq!(g.codes, vec![ModCode::Char(b'm')]);
        assert_eq!(g.status, Some(b'?'));
        assert_eq!(g.deltas, vec![5, 12, 0]);
        assert_eq!(g.ml, vec![200, 10, 128]);
    }

    #[test]
    fn multi_code_group_takes_two_ml_per_position() {
        // C+mh with 2 positions -> 4 ML bytes, position-major.
        let m = parse(b"C+mh,1,3;", &[10, 20, 30, 40]);
        let g = &m.groups[0];
        assert_eq!(g.codes, vec![ModCode::Char(b'm'), ModCode::Char(b'h')]);
        assert_eq!(g.deltas, vec![1, 3]);
        assert_eq!(g.ml, vec![10, 20, 30, 40]);
    }

    #[test]
    fn chebi_numeric_code() {
        let m = parse(b"C+16061,2;", &[99]);
        assert_eq!(m.groups[0].codes, vec![ModCode::Chebi(16061)]);
        assert_eq!(m.groups[0].deltas, vec![2]);
    }

    #[test]
    fn two_groups_split_ml() {
        let m = parse(b"C+m,0;A+a,1,4;", &[1, 2, 3]);
        assert_eq!(m.groups.len(), 2);
        assert_eq!(m.groups[0].ml, vec![1]); // 1 position
        assert_eq!(m.groups[1].ml, vec![2, 3]); // 2 positions
        assert_eq!(m.groups[1].base, b'A');
    }

    #[test]
    fn no_status_and_empty_positions() {
        let m = parse(b"C+m;", &[]);
        let g = &m.groups[0];
        assert_eq!(g.status, None);
        assert!(g.deltas.is_empty());
        assert!(g.ml.is_empty());
    }
}
