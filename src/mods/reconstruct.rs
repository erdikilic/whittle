use super::{counting_base, MmGroup, Mods};

pub fn reconstruct(mods: &Mods, seq: &[u8], start: usize, end: usize) -> Mods {
    let mut out = Vec::new();

    for g in &mods.groups {
        let ncodes = g.codes.len().max(1);
        let cbase = counting_base(g.base, g.strand);

        // All positions of the counting base along the whole SEQ (ascending).
        let positions: Vec<usize> = seq
            .iter()
            .enumerate()
            .filter(|&(_, &b)| b.to_ascii_uppercase() == cbase)
            .map(|(i, _)| i)
            .collect();

        // Positions inside the output window (ascending), for renumbering.
        let window: Vec<usize> = positions
            .iter()
            .copied()
            .filter(|&p| p >= start && p < end)
            .collect();

        // Walk the group's deltas to recover each modified absolute position and
        // its ML byte run, then keep the ones inside the window.
        let mut new_deltas = Vec::new();
        let mut new_ml = Vec::new();
        let mut prev_widx: isize = -1;
        let mut cursor = 0usize; // index into `positions`

        for (k, &d) in g.deltas.iter().enumerate() {
            cursor += d;
            if cursor >= positions.len() {
                break; // malformed / past end
            }
            let abs = positions[cursor];
            cursor += 1;

            if abs < start || abs >= end {
                continue;
            }
            // Index of `abs` within the window (present by construction).
            let widx = window.partition_point(|&p| p < abs) as isize;
            new_deltas.push((widx - prev_widx - 1) as usize);
            prev_widx = widx;

            let ml_start = (k * ncodes).min(g.ml.len());
            let ml_end = (ml_start + ncodes).min(g.ml.len());
            new_ml.extend_from_slice(&g.ml[ml_start..ml_end]);
        }

        if !new_deltas.is_empty() {
            out.push(MmGroup {
                base: g.base,
                strand: g.strand,
                codes: g.codes.clone(),
                status: g.status,
                deltas: new_deltas,
                ml: new_ml,
            });
        }
    }

    Mods { groups: out }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::parse;

    // Helper: build from raw MM/ML then reconstruct.
    fn recon(mm: &[u8], ml: &[u8], seq: &[u8], start: usize, end: usize) -> Mods {
        reconstruct(&parse(mm, ml), seq, start, end)
    }

    #[test]
    fn keeps_only_in_window_and_renumbers() {
        // seq C at indices 0,2,5,8. MM "C+m,0,1,0" -> modified Cs at occ 0,2,3 => pos 0,5,8.
        // ML one byte per position: [11,22,33].
        let seq = b"CACATCTTC"; // 'C' at 0,2,5,8
        let m = recon(b"C+m,0,1,0;", &[11, 22, 33], seq, 3, 9);
        // window [3,9): C occurrences at 5,8. Surviving modified: pos5 (occ idx1), pos8 (occ idx2).
        // window C positions = [5,8] -> renumber. pos5 -> delta 0; pos8 -> delta 0.
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].deltas, vec![0, 0]);
        assert_eq!(m.groups[0].ml, vec![22, 33]);
    }

    #[test]
    fn drops_group_with_no_survivors() {
        let seq = b"CCCC";
        let m = recon(b"C+m,0;", &[50], seq, 2, 4); // modified C at pos0, outside [2,4)
        assert!(m.groups.is_empty());
    }

    #[test]
    fn multi_code_keeps_both_ml_bytes_per_position() {
        let seq = b"CC";
        // C+mh with 2 positions (0 and next) -> ML [a0,b0,a1,b1]. Keep both in [0,2).
        let m = recon(b"C+mh,0,0;", &[1, 2, 3, 4], seq, 0, 2);
        assert_eq!(m.groups[0].ml, vec![1, 2, 3, 4]);
        assert_eq!(m.groups[0].deltas, vec![0, 0]);
    }

    #[test]
    fn minus_strand_counts_complement() {
        // G-m: counting base = complement(G) = C. seq C at 1,3. "G-m,0,0" -> modified at C occ 0,1 => pos1,3.
        let seq = b"ACAC";
        let m = recon(b"G-m,0,0;", &[7, 8], seq, 2, 4); // window keeps pos3 only
        assert_eq!(m.groups[0].deltas, vec![0]);
        assert_eq!(m.groups[0].ml, vec![8]);
        assert_eq!((m.groups[0].base, m.groups[0].strand), (b'G', b'-'));
    }

    #[test]
    fn truncated_ml_does_not_panic() {
        // MM lists 3 modified C's but ML has only 1 byte (malformed/truncated).
        let m = reconstruct(&crate::mods::parse(b"C+m,0,0,0;", &[5]), b"CCC", 0, 3);
        // All three C's are in-window; deltas renumber to [0,0,0]; only the ML bytes
        // that actually exist are kept (no panic).
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].deltas, vec![0, 0, 0]);
        assert_eq!(m.groups[0].ml, vec![5]);
    }

    #[test]
    fn renumber_with_nonzero_gap_and_offset_first() {
        // seq CCCCC (C at 0,1,2,3,4). MM C+m,0,1,0 -> modified at occ 0,2,3 => abs 0,2,3, ML [a,b,c].
        // Window [1,5): drop abs 0 (outside); keep abs 2 and 3.
        // window C positions = [1,2,3,4]; widx(2)=1 -> delta 1; widx(3)=2 -> delta 0.
        let m = reconstruct(&crate::mods::parse(b"C+m,0,1,0;", &[10, 20, 30]), b"CCCCC", 1, 5);
        assert_eq!(m.groups[0].deltas, vec![1, 0]);
        assert_eq!(m.groups[0].ml, vec![20, 30]);
    }
}
