use super::{ModCode, Mods};

pub fn serialize(mods: &Mods) -> (Vec<u8>, Vec<u8>) {
    let mut mm = Vec::new();
    let mut ml = Vec::new();

    for g in &mods.groups {
        if g.deltas.is_empty() {
            continue;
        }
        mm.push(g.base);
        mm.push(g.strand);
        for code in &g.codes {
            match code {
                ModCode::Char(c) => mm.push(*c),
                ModCode::Chebi(id) => mm.extend_from_slice(id.to_string().as_bytes()),
            }
        }
        if let Some(s) = g.status {
            mm.push(s);
        }
        for d in &g.deltas {
            mm.push(b',');
            mm.extend_from_slice(d.to_string().as_bytes());
        }
        mm.push(b';');
        ml.extend_from_slice(&g.ml);
    }

    (mm, ml)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::parse;

    #[test]
    fn roundtrip_single_group() {
        let (mm, ml) = serialize(&parse(b"C+m?,5,12,0;", &[200, 10, 128]));
        assert_eq!(mm, b"C+m?,5,12,0;");
        assert_eq!(ml, vec![200, 10, 128]);
    }

    #[test]
    fn roundtrip_multi_group_and_chebi() {
        let input = b"C+mh,1,3;A+16061,2;".as_slice();
        let (mm, ml) = serialize(&parse(input, &[1, 2, 3, 4, 5]));
        assert_eq!(mm, input);
        assert_eq!(ml, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn skips_empty_groups() {
        // A group with no deltas should not be emitted.
        let mut mods = parse(b"C+m,0;", &[7]);
        mods.groups.push(crate::mods::MmGroup {
            base: b'A', strand: b'+', codes: vec![crate::mods::ModCode::Char(b'a')],
            status: None, deltas: vec![], ml: vec![],
        });
        let (mm, _) = serialize(&mods);
        assert_eq!(mm, b"C+m,0;");
    }
}
