//! Parsing of `MM` text and `ML` bytes into `Mods`.

use super::{MmGroup, ModCode, Mods};

/// The header of one group token, read in place: `[A-Za-z][+-]` then either
/// a run of single-letter codes or one numeric ChEBI id, then an optional
/// `.`/`?` status. `tail` is the rest of the token, the `(,[0-9]+)*` delta
/// list.
struct GroupHead<'a> {
    base: u8,
    strand: u8,
    /// The code text: letters, or the digits of a ChEBI id when `chebi`.
    codes: &'a [u8],
    chebi: bool,
    status: Option<u8>,
    tail: &'a [u8],
}

impl GroupHead<'_> {
    /// The number of modification codes the group lists. A ChEBI id is one.
    fn code_count(&self) -> usize {
        if self.chebi { 1 } else { self.codes.len() }
    }

    /// The group's codes as `ModCode`s.
    fn mod_codes(&self) -> Vec<ModCode> {
        if self.chebi {
            // Saturating: a corrupt over-long id clamps instead of overflowing.
            let id = self.codes.iter().fold(0u32, |id, &d| {
                id.saturating_mul(10).saturating_add(u32::from(d - b'0'))
            });
            vec![ModCode::Chebi(id)]
        } else {
            self.codes.iter().map(|&c| ModCode::Char(c)).collect()
        }
    }
}

/// Reads a group token's header. `None` when the token does not start with
/// `[A-Za-z][+-]([a-z]+|[0-9]+)`.
fn group_head(token: &[u8]) -> Option<GroupHead<'_>> {
    let (&base, rest) = token.split_first()?;
    if !base.is_ascii_alphabetic() {
        return None;
    }
    let (&strand, rest) = rest.split_first()?;
    if !matches!(strand, b'+' | b'-') {
        return None;
    }
    let chebi = rest.first().is_some_and(u8::is_ascii_digit);
    let code_len = rest
        .iter()
        .take_while(|b| {
            if chebi {
                b.is_ascii_digit()
            } else {
                b.is_ascii_alphabetic()
            }
        })
        .count();
    if code_len == 0 {
        return None;
    }
    let (codes, mut tail) = rest.split_at(code_len);
    let mut status = None;
    if let Some((&s, rest)) = tail.split_first()
        && matches!(s, b'.' | b'?')
    {
        status = Some(s);
        tail = rest;
    }
    Some(GroupHead {
        base,
        strand,
        codes,
        chebi,
        status,
        tail,
    })
}

/// Reads the deltas of a group's tail, `(,[0-9]+)*`, in place: `Some(delta)`
/// per listed position, then `None` once when the tail holds a byte the
/// grammar rejects. A single trailing comma lists no position and is accepted.
fn deltas(tail: &[u8]) -> impl Iterator<Item = Option<usize>> + '_ {
    let mut rest = tail;
    let mut done = false;
    std::iter::from_fn(move || {
        if done {
            return None;
        }
        let Some((&b',', digits)) = rest.split_first() else {
            done = true;
            return (!rest.is_empty()).then_some(None);
        };
        let n = digits.iter().take_while(|b| b.is_ascii_digit()).count();
        if n == 0 {
            done = true;
            return (!digits.is_empty()).then_some(None);
        }
        // Saturating: a delta this large lies outside any window and is dropped
        // by `reconstruct`.
        let delta = digits[..n].iter().fold(0usize, |d, &b| {
            d.saturating_mul(10).saturating_add(usize::from(b - b'0'))
        });
        rest = &digits[n..];
        Some(Some(delta))
    })
}

/// Returns the group tokens of `mm`. The empty remainder after a final `;` is
/// not a token; an empty token anywhere else is one, and fails the grammar
/// like any group without a code.
fn group_tokens(mm: &[u8]) -> impl Iterator<Item = &[u8]> {
    let len = mm.len();
    let mut next = 0usize;
    mm.split(|&b| b == b';').filter(move |token| {
        let start = next;
        next += token.len() + 1;
        start < len
    })
}

/// Returns the number of `ML` bytes a well-formed record carrying this `MM`
/// string must have: one per listed position per mod code, summed over groups.
/// Counts without allocating, so a caller can check `ML` on the hot path
/// without building the groups. `None` when `mm` does not conform to the
/// grammar `[A-Za-z][+-]([a-z]+|[0-9]+)[.?]?(,[0-9]+)*` to its end.
pub fn expected_ml_len(mm: &[u8]) -> Option<usize> {
    let mut total = 0usize;
    for token in group_tokens(mm) {
        let head = group_head(token)?;
        let positions = deltas(head.tail).try_fold(0usize, |n, delta| delta.map(|_| n + 1))?;
        total += positions * head.code_count();
    }
    Some(total)
}

/// Validates every cumulative modification occurrence against the sequence.
pub(crate) fn positions_valid(mm: &[u8], seq: impl Iterator<Item = u8>) -> bool {
    let mut totals = [0usize; 256];
    let mut len = 0;
    for base in seq {
        totals[usize::from(base.to_ascii_uppercase())] += 1;
        len += 1;
    }
    totals[usize::from(b'N')] = len;
    group_tokens(mm).all(|token| {
        let Some(head) = group_head(token) else {
            return false;
        };
        let total = totals[usize::from(super::counting_base(head.base))];
        let mut next = 0usize;
        deltas(head.tail).all(|delta| {
            let Some(position) = delta.and_then(|d| next.checked_add(d)) else {
                return false;
            };
            if position >= total {
                return false;
            }
            next = position + 1;
            true
        })
    })
}

/// Parses a raw `MM:Z` string plus its `ML:B,C` array into groups. A group
/// without a usable header contributes nothing, and a group is read up to its
/// first unexpected byte; the groups after it are still read. The workflows
/// validate `MM` with `expected_ml_len` first, so a production parse sees only
/// well-formed strings.
pub fn parse(mm: &[u8], ml: &[u8]) -> Mods {
    let mut groups = Vec::new();
    let mut ml_pos = 0usize;

    for token in group_tokens(mm) {
        let Some(head) = group_head(token) else {
            continue;
        };
        let codes = head.mod_codes();
        let deltas: Vec<usize> = deltas(head.tail).map_while(|delta| delta).collect();

        // This group's ML bytes: positions * codes, position-major.
        let want = deltas.len() * codes.len();
        let end = (ml_pos + want).min(ml.len());
        let group_ml = ml[ml_pos..end].to_vec();
        ml_pos = end;

        groups.push(MmGroup {
            base: head.base,
            strand: head.strand,
            codes,
            status: head.status,
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
    fn cumulative_positions_must_fit_the_counting_base() {
        let seq = b"ACCT";
        for mm in [b"C+m,0,0;".as_slice(), b"N+n,3;", b"U+m,0;", b"G-m;"] {
            assert!(positions_valid(mm, seq.iter().copied()));
        }
        for mm in [
            b"C+m,2;".as_slice(),
            b"C+m,1,0;",
            b"N+n,4;",
            b"G-m,0;",
            b"C+m,99999999999999999999999999;",
        ] {
            assert!(!positions_valid(mm, seq.iter().copied()));
        }
    }

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
        // `C+mh` with 2 positions takes 4 ML bytes, position-major.
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

    /// A corrupt 20-digit ChEBI id or delta overflows `u32`/`usize` under
    /// `n * 10 + d`; saturating arithmetic clamps instead of panicking (debug
    /// builds panic on overflow).
    #[test]
    fn over_long_numeric_fields_saturate_without_panicking() {
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
            let parsed = parse(mm, &[]);
            let want: usize = parsed
                .groups
                .iter()
                .map(|g| g.deltas.len() * g.codes.len())
                .sum();
            assert_eq!(
                expected_ml_len(mm),
                Some(want),
                "Counting scan disagreed for {}",
                String::from_utf8_lossy(mm)
            );
        }
    }

    /// The validator accepts well-formed strings and the grammar's accepted
    /// oddities (a trailing comma, uppercase codes, an empty string) with the
    /// declared `ML` length, and rejects every malformed shape.
    #[test]
    fn expected_ml_len_rejects_every_malformed_shape() {
        let cases: &[(&[u8], Option<usize>)] = &[
            (b"", Some(0)),
            (b"C+m,0,5;", Some(2)),
            (b"C+m,0;A+a,1,2;", Some(3)),
            (b"C+mh,1,2;", Some(4)),
            (b"C+MH,1,2;", Some(4)),
            (b"C+12345,3;", Some(1)),
            (b"C+12345;", Some(0)),
            (b"C+m?,1;", Some(1)),
            (b"C+m.,1;", Some(1)),
            (b"C+m;", Some(0)),
            (b"C+m", Some(0)),
            (b"C+m,", Some(0)),
            (b"C+m,1,", Some(1)),
            (b"C+m,1,;A+a,2", Some(2)),
            (b"N-x,0", Some(1)),
            (b";", None),
            (b"C+m?.,1;", None),
            (b"C+m,,1;", None),
            (b"C+;", None),
            (b"C+", None),
            (b"+m,1;", None),
            (b"Cm,1;", None),
            (b"C*m,1;", None),
            (b"C+m,a;", None),
            (b"C+m,1 ;", None),
            (b"C+m,1;;C+h;", None),
            (b"C+m1;", None),
            (b"C+m,1;junk", None),
            (b"C+m,1;C+h,2,x", None),
            (b"C+m,-1;", None),
            (b"C+m,1.5;", None),
            (b"9+m,1;", None),
            (b"C+m,5,1x,7;", None),
            (b"C,5;", None),
        ];
        for &(mm, want) in cases {
            assert_eq!(expected_ml_len(mm), want, "{}", String::from_utf8_lossy(mm));
        }
    }
}
