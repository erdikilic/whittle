use super::{MmGroup, ModCode, Mods};

/// The number of `ML` bytes a well-formed record carrying this `MM` string must
/// have: one per modified position per mod code, summed over groups.
///
/// Counts without allocating, so a caller can check `ML` consistency on the hot
/// path without building the parsed representation. Tokenization mirrors `parse`,
/// including its tolerance of a malformed tail.
pub fn expected_ml_len(mm: &[u8]) -> usize {
    let mut total = 0usize;
    for token in mm.split(|&b| b == b';') {
        if token.len() < 2 {
            continue;
        }
        let mut i = 2;
        let mut codes = 0usize;
        if i < token.len() && token[i].is_ascii_digit() {
            while i < token.len() && token[i].is_ascii_digit() {
                i += 1;
            }
            codes = 1; // a ChEBI id is one code
        } else {
            while i < token.len() && token[i].is_ascii_alphabetic() {
                codes += 1;
                i += 1;
            }
        }
        if i < token.len() && (token[i] == b'.' || token[i] == b'?') {
            i += 1;
        }
        let mut deltas = 0usize;
        while i < token.len() && token[i] == b',' {
            i += 1;
            let mut saw = false;
            while i < token.len() && token[i].is_ascii_digit() {
                i += 1;
                saw = true;
            }
            // An empty field (`,,`) contributes no delta but does not end the
            // list, matching `parse`, which skips it and keeps reading.
            if saw {
                deltas += 1;
            }
        }
        total += deltas * codes.max(1);
    }
    total
}

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
                // Saturating: a corrupt over-long id must clamp, never overflow
                // (which would panic in debug / silently wrap in release).
                id = id
                    .saturating_mul(10)
                    .saturating_add((token[i] - b'0') as u32);
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
                // Saturating for the same reason as the ChEBI id above; a delta
                // this large is unreachable and gets dropped in reconstruct.
                n = n
                    .saturating_mul(10)
                    .saturating_add((token[i] - b'0') as usize);
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

        groups.push(MmGroup {
            base,
            strand,
            codes,
            status,
            deltas,
            ml: group_ml,
        });
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

    #[test]
    fn over_long_numeric_fields_saturate_without_panicking() {
        // A corrupt 20-digit ChEBI id and delta overflow u32/usize with the naive
        // `n*10 + d`; saturating arithmetic must clamp instead of panicking (this
        // test is a debug build, where overflow panics).
        let m = parse(b"C+99999999999999999999,88888888888888888888;", &[1]);
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].codes, vec![ModCode::Chebi(u32::MAX)]);
        assert_eq!(m.groups[0].deltas, vec![usize::MAX]);
    }
    /// The counting scan must agree with the parsed representation, since it is
    /// what lets the full-window shortcut check ML consistency without parsing.
    #[test]
    fn expected_ml_len_matches_the_parsed_length() {
        for mm in [
            &b"C+m,0,1,2;"[..],
            b"C+m,0,1;C+h,2;",
            b"C+mh,0,1;",
            b"A+a?,3;",
            b"C+16061,0,1;",
            b"C+m;",
            b"C+m,0,1,2;A+a,0;N+n,4,1;",
            b"",
            b";",
            b"C+m,0,,1;",
        ] {
            let parsed = super::parse(mm, &[]);
            let want: usize = parsed
                .groups
                .iter()
                .map(|g| g.deltas.len() * g.codes.len().max(1))
                .sum();
            assert_eq!(
                super::expected_ml_len(mm),
                want,
                "counting scan disagreed for {}",
                String::from_utf8_lossy(mm)
            );
        }
    }
}
