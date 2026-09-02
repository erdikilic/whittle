use super::{MmGroup, ModCode, Mods};

/// An `MM` string holding a byte the group grammar does not accept.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MalformedMm {
    /// Byte offset of the first unexpected byte.
    pub offset: usize,
}

/// One group token read in place: its header fields, its code and delta
/// counts, and the offset of the first byte the grammar rejects, which is
/// `token.len()` for a well-formed group.
struct GroupScan {
    base: u8,
    strand: u8,
    status: Option<u8>,
    codes: usize,
    deltas: usize,
    stop: usize,
}

/// Reads one group token without allocating. Codes and deltas are handed to
/// `on_code` and `on_delta` as they are read, so `parse` collects them while
/// the counting scan discards them. The grammar is the SAM tags one:
/// `[A-Za-z][+-]([a-z]+|[0-9]+)[.?]?(,[0-9]+)*`.
fn scan_group(
    token: &[u8],
    mut on_code: impl FnMut(ModCode),
    mut on_delta: impl FnMut(usize),
) -> GroupScan {
    let mut scan = GroupScan {
        base: 0,
        strand: 0,
        status: None,
        codes: 0,
        deltas: 0,
        stop: 0,
    };
    let Some(&base) = token.first().filter(|b| b.is_ascii_alphabetic()) else {
        return scan;
    };
    scan.base = base;
    scan.stop = 1;
    let Some(&strand) = token.get(1).filter(|s| matches!(s, b'+' | b'-')) else {
        return scan;
    };
    scan.strand = strand;

    let mut i = 2;
    if i < token.len() && token[i].is_ascii_digit() {
        // Saturating: a corrupt over-long id clamps instead of overflowing.
        let mut id = 0u32;
        while i < token.len() && token[i].is_ascii_digit() {
            id = id
                .saturating_mul(10)
                .saturating_add(u32::from(token[i] - b'0'));
            i += 1;
        }
        on_code(ModCode::Chebi(id));
        scan.codes = 1;
    } else {
        while i < token.len() && token[i].is_ascii_alphabetic() {
            on_code(ModCode::Char(token[i]));
            scan.codes += 1;
            i += 1;
        }
    }
    scan.stop = i;
    if scan.codes == 0 {
        return scan;
    }

    if i < token.len() && matches!(token[i], b'.' | b'?') {
        scan.status = Some(token[i]);
        i += 1;
    }

    while i < token.len() && token[i] == b',' {
        i += 1;
        let digits = i;
        let mut n = 0usize;
        while i < token.len() && token[i].is_ascii_digit() {
            // Saturating for the same reason as the ChEBI id; a delta this
            // large lies outside any window and is dropped by `reconstruct`.
            n = n
                .saturating_mul(10)
                .saturating_add(usize::from(token[i] - b'0'));
            i += 1;
        }
        if i == digits {
            break;
        }
        on_delta(n);
        scan.deltas += 1;
    }
    scan.stop = i;
    scan
}

/// The group tokens of `mm` with their byte offsets. The empty remainder after
/// a final `;` is not a token; an empty token anywhere else is one, and fails
/// the grammar like any group without a code.
fn group_tokens(mm: &[u8]) -> impl Iterator<Item = (usize, &[u8])> {
    let len = mm.len();
    let mut next = 0usize;
    mm.split(|&b| b == b';').filter_map(move |token| {
        let start = next;
        next += token.len() + 1;
        (start < len).then_some((start, token))
    })
}

/// The number of `ML` bytes a well-formed record carrying this `MM` string must
/// have: one per listed position per mod code, summed over groups. Counts
/// without allocating, so a caller can check `ML` on the hot path without
/// building the groups. `Err` when `mm` does not conform to the grammar.
pub fn expected_ml_len(mm: &[u8]) -> Result<usize, MalformedMm> {
    let mut total = 0usize;
    for (start, token) in group_tokens(mm) {
        let scan = scan_group(token, |_| {}, |_| {});
        if scan.codes == 0 || scan.stop != token.len() {
            return Err(MalformedMm {
                offset: start + scan.stop,
            });
        }
        total += scan.deltas * scan.codes;
    }
    Ok(total)
}

/// Parses a raw `MM:Z` string plus its `ML:B,C` array into groups. A group is
/// read up to its first unexpected byte and the remainder of that group is
/// skipped; the groups after it are still read. `parse_checked` refuses such
/// strings instead.
pub fn parse(mm: &[u8], ml: &[u8]) -> Mods {
    parse_inner(mm, ml).0
}

/// Parses like `parse` but refuses an `MM` string that does not conform to the
/// group grammar to its end.
pub fn parse_checked(mm: &[u8], ml: &[u8]) -> Result<Mods, MalformedMm> {
    match parse_inner(mm, ml) {
        (mods, None) => Ok(mods),
        (_, Some(offset)) => Err(MalformedMm { offset }),
    }
}

/// The parsed groups and the offset of the first unexpected byte, if any.
fn parse_inner(mm: &[u8], ml: &[u8]) -> (Mods, Option<usize>) {
    let mut groups = Vec::new();
    let mut ml_pos = 0usize;
    let mut malformed = None;

    for (start, token) in group_tokens(mm) {
        let mut codes = Vec::new();
        let mut deltas = Vec::new();
        let scan = scan_group(token, |c| codes.push(c), |d| deltas.push(d));
        if (scan.codes == 0 || scan.stop != token.len()) && malformed.is_none() {
            malformed = Some(start + scan.stop);
        }
        if scan.codes == 0 {
            continue;
        }

        // Claim this group's ML bytes: positions * codes, position-major.
        let want = deltas.len() * scan.codes;
        let end = (ml_pos + want).min(ml.len());
        let group_ml = ml[ml_pos..end].to_vec();
        ml_pos = end;

        groups.push(MmGroup {
            base: scan.base,
            strand: scan.strand,
            codes,
            status: scan.status,
            deltas,
            ml: group_ml,
        });
    }

    (Mods { groups }, malformed)
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
            b"C+m,0",
        ] {
            let parsed = parse_checked(mm, &[]).unwrap();
            let want: usize = parsed
                .groups
                .iter()
                .map(|g| g.deltas.len() * g.codes.len())
                .sum();
            assert_eq!(
                expected_ml_len(mm),
                Ok(want),
                "counting scan disagreed for {}",
                String::from_utf8_lossy(mm)
            );
        }
    }

    /// Every departure from the grammar is reported at the offending byte by
    /// both the checked parse and the counting scan.
    #[test]
    fn malformed_strings_report_the_first_unexpected_byte() {
        for (mm, offset) in [
            (&b"C+m,5,1x,7;"[..], 7),
            (b"C+m,0;;A+a,1;", 6),
            (b"C+m,,1;", 4),
            (b"C,5;", 1),
            (b"C+;", 2),
            (b"+m,0;", 0),
            (b";", 0),
            (b"C+m,0,1;C+h,2 ;", 13),
        ] {
            let want = Some(MalformedMm { offset });
            assert_eq!(
                parse_checked(mm, &[]).err(),
                want,
                "{}",
                String::from_utf8_lossy(mm)
            );
            assert_eq!(
                expected_ml_len(mm).err(),
                want,
                "{}",
                String::from_utf8_lossy(mm)
            );
        }
    }

    /// The lenient parse keeps what precedes the unexpected byte and still
    /// reads the groups after it.
    #[test]
    fn lenient_parse_keeps_the_readable_prefix_of_a_malformed_group() {
        let m = parse(b"C+m,5,1x,7;A+a,2;", &[1, 2, 3]);
        assert_eq!(m.groups.len(), 2);
        assert_eq!(m.groups[0].deltas, vec![5, 1]);
        assert_eq!(m.groups[0].ml, vec![1, 2]);
        assert_eq!(m.groups[1].deltas, vec![2]);
        assert_eq!(m.groups[1].ml, vec![3]);
        // A group without a usable header contributes nothing.
        assert_eq!(parse(b"C+;A+a,2;", &[3]).groups.len(), 1);
    }
}
