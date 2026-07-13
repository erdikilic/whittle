use super::{ModCode, Mods};
use crate::io::fastq::push_u64;

pub fn serialize(mods: &Mods) -> (Vec<u8>, Vec<u8>) {
    let mut mm = Vec::new();
    let mut ml = Vec::new();

    for g in &mods.groups {
        // Empty-delta groups ARE emitted: `reconstruct` only ever keeps an
        // empty group that was empty in the source (a valid "assessed, none
        // found" record), and dropping it here would silently erase that
        // information. A group that lost all its positions to the window is
        // already excluded upstream, so it never reaches here.
        mm.push(g.base);
        mm.push(g.strand);
        for code in &g.codes {
            match code {
                ModCode::Char(c) => mm.push(*c),
                ModCode::Chebi(id) => push_u64(&mut mm, u64::from(*id)),
            }
        }
        if let Some(s) = g.status {
            mm.push(s);
        }
        for d in &g.deltas {
            mm.push(b',');
            push_u64(&mut mm, *d as u64);
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
    fn emits_empty_groups() {
        // A zero-position group is valid SAM ("assessed, none found") and must
        // be re-emitted, not dropped. Contributes no ML bytes.
        let mut mods = parse(b"C+m,0;", &[7]);
        mods.groups.push(crate::mods::MmGroup {
            base: b'A',
            strand: b'+',
            codes: vec![crate::mods::ModCode::Char(b'a')],
            status: None,
            deltas: vec![],
            ml: vec![],
        });
        let (mm, ml) = serialize(&mods);
        assert_eq!(mm, b"C+m,0;A+a;");
        assert_eq!(ml, vec![7]);
    }
}
