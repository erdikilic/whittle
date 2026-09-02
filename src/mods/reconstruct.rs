//! Windowing of parsed modifications to a trimmed sequence interval.

use super::{MmGroup, Mods, counting_base, counts};

/// Restricts every group to the window `[start, end)`: skip-counts are
/// renumbered against the window's own counting-base occurrences and each
/// surviving position keeps its `ML` bytes. Every source group is emitted,
/// with no positions when none survive, because a group declares how the
/// bases it does not list are read (canonical under `.` or no status, no call
/// under `?`) and removing it would change that reading for the whole window.
pub fn reconstruct(mods: &Mods, seq: &[u8], start: usize, end: usize) -> Mods {
    let mut out = Vec::new();
    // MM deltas index occurrences of the canonical base, so reconstructing a
    // window only needs its occurrence-index bounds. The bounds are cached for
    // groups sharing a counting base (e.g. `C+h` and `C+m`) without
    // materializing every matching SEQ position.
    let mut bounds_cache: Vec<(u8, usize, usize)> = Vec::new();

    for g in &mods.groups {
        let ncodes = g.codes.len().max(1);
        let cbase = counting_base(g.base);

        // Occurrence indexes in `[before, window_end)` are inside `[start, end)`.
        // An index lookup avoids holding a borrow across the conditional insert.
        let cache_idx = match bounds_cache.iter().position(|(b, _, _)| *b == cbase) {
            Some(i) => i,
            None => {
                let before = seq[..start].iter().filter(|&&b| counts(b, cbase)).count();
                let in_window = seq[start..end]
                    .iter()
                    .filter(|&&b| counts(b, cbase))
                    .count();
                bounds_cache.push((cbase, before, before + in_window));
                bounds_cache.len() - 1
            },
        };
        let (_, before, window_end) = bounds_cache[cache_idx];

        // Each delta yields a modified absolute position and its ML byte run;
        // those inside the window are kept.
        let mut new_deltas = Vec::new();
        let mut new_ml = Vec::new();
        let mut prev_widx: isize = -1;
        let mut cursor = 0usize; // occurrence index of the counting base

        for (k, &d) in g.deltas.iter().enumerate() {
            // Saturating so a corrupt (clamped) delta cannot overflow the running
            // cursor; it remains outside the finite window bounds.
            cursor = cursor.saturating_add(d);
            let occurrence = cursor;
            cursor = cursor.saturating_add(1);

            if occurrence < before || occurrence >= window_end {
                continue;
            }
            let widx = (occurrence - before) as isize;
            new_deltas.push((widx - prev_widx - 1) as usize);
            prev_widx = widx;

            let ml_start = (k * ncodes).min(g.ml.len());
            let ml_end = (ml_start + ncodes).min(g.ml.len());
            new_ml.extend_from_slice(&g.ml[ml_start..ml_end]);
        }

        out.push(MmGroup {
            base: g.base,
            strand: g.strand,
            codes: g.codes.clone(),
            status: g.status,
            deltas: new_deltas,
            ml: new_ml,
        });
    }

    Mods { groups: out }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::parse;

    /// Parses raw MM/ML and reconstructs the window.
    fn recon(mm: &[u8], ml: &[u8], seq: &[u8], start: usize, end: usize) -> Mods {
        reconstruct(&parse(mm, ml), seq, start, end)
    }

    #[test]
    fn keeps_only_in_window_and_renumbers() {
        // C at indices 0, 2, 5, 8. `C+m,0,1,0` marks occurrences 0, 2, 3, which
        // are positions 0, 5, 8. ML has one byte per position: [11, 22, 33].
        let seq = b"CACATCTTC"; // 'C' at 0,2,5,8
        let m = recon(b"C+m,0,1,0;", &[11, 22, 33], seq, 3, 9);
        // Window [3,9) holds C occurrences 5 and 8; the survivors are position 5
        // (occurrence 1) and position 8 (occurrence 2). Renumbered against the
        // window's C positions [5, 8], both deltas are 0.
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].deltas, vec![0, 0]);
        assert_eq!(m.groups[0].ml, vec![22, 33]);
    }

    /// An MM group with zero positions (e.g. `C+m?;`, assessed, none found) is
    /// valid SAM and survives windowing.
    #[test]
    fn preserves_originally_empty_group() {
        let seq = b"ACAC"; // A at 0,2 ; C at 1,3
        // `A+a` has a modified A at occurrence 0 (position 0, in window); `C+m?`
        // is empty.
        let m = recon(b"A+a,0;C+m?;", &[9], seq, 0, 4);
        assert_eq!(m.groups.len(), 2, "Empty C+m? group is kept alongside A+a");
        let cm = m
            .groups
            .iter()
            .find(|g| g.base == b'C')
            .expect("Empty C+m? group is preserved");
        assert!(cm.deltas.is_empty());
        assert!(cm.ml.is_empty());
        assert_eq!(cm.status, Some(b'?'));
    }

    /// A group whose listed positions all fall outside the window is emitted
    /// with no positions: under `.` (or no status) it still declares the
    /// window's unlisted bases canonical, which its absence would not.
    #[test]
    fn keeps_group_with_no_survivors_as_empty() {
        let seq = b"CCCC";
        let m = recon(b"C+m.,0;", &[50], seq, 2, 4); // modified C at pos0, outside [2,4)
        assert_eq!(m.groups.len(), 1);
        let g = &m.groups[0];
        assert!(g.deltas.is_empty());
        assert!(g.ml.is_empty());
        assert_eq!(g.status, Some(b'.'));
        assert_eq!((g.base, g.strand), (b'C', b'+'));
    }

    #[test]
    fn multi_code_keeps_both_ml_bytes_per_position() {
        let seq = b"CC";
        // `C+mh` with 2 positions has ML [a0, b0, a1, b1]; both are kept in [0,2).
        let m = recon(b"C+mh,0,0;", &[1, 2, 3, 4], seq, 0, 2);
        assert_eq!(m.groups[0].ml, vec![1, 2, 3, 4]);
        assert_eq!(m.groups[0].deltas, vec![0, 0]);
    }

    /// Pinned against htslib, which counts the literal fundamental base on both
    /// strands. Counting the complement would move every reverse-strand call
    /// onto a different base.
    #[test]
    fn minus_strand_counts_the_literal_base() {
        // G at 1 and 3. `G-m,0,0` marks G occurrences 0 and 1, positions 1 and 3.
        let seq = b"AGAG";
        let m = recon(b"G-m,0,0;", &[7, 8], seq, 2, 4); // window [2,4) keeps abs 3 only
        assert_eq!(m.groups[0].deltas, vec![0]);
        assert_eq!(m.groups[0].ml, vec![8]);
        assert_eq!((m.groups[0].base, m.groups[0].strand), (b'G', b'-'));
    }

    /// BAM SEQ stores `T`, never `U`, so an RNA group must count `T`.
    #[test]
    fn uracil_group_counts_thymine_positions() {
        let m = recon(b"U+m,0,0;", &[7, 8], b"ATAT", 2, 4);
        assert_eq!(m.groups[0].deltas, vec![0]);
        assert_eq!(m.groups[0].ml, vec![8]);
    }

    /// `N` in MM means "any base", so its skip-counts index every position.
    #[test]
    fn n_group_counts_every_base() {
        let m = recon(b"N+n,2,1;", &[7, 8], b"ACGTACGT", 0, 8);
        assert_eq!(m.groups[0].deltas, vec![2, 1]);
        assert_eq!(m.groups[0].ml, vec![7, 8]);

        // Window [3,8) drops position 2 and keeps position 4, which becomes
        // window index 1.
        let m = recon(b"N+n,2,1;", &[7, 8], b"ACGTACGT", 3, 8);
        assert_eq!(m.groups[0].deltas, vec![1]);
        assert_eq!(m.groups[0].ml, vec![8]);
    }

    #[test]
    fn truncated_ml_does_not_panic() {
        // MM lists 3 modified Cs but ML has only 1 byte (truncated).
        let m = reconstruct(&crate::mods::parse(b"C+m,0,0,0;", &[5]), b"CCC", 0, 3);
        // All three Cs are in the window; deltas renumber to [0, 0, 0]; only the
        // ML bytes present are kept.
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].deltas, vec![0, 0, 0]);
        assert_eq!(m.groups[0].ml, vec![5]);
    }

    #[test]
    fn renumber_with_nonzero_gap_and_offset_first() {
        // CCCCC has C at 0..=4. `C+m,0,1,0` marks occurrences 0, 2, 3, which are
        // positions 0, 2, 3, with ML [a, b, c]. Window [1,5) drops position 0 and
        // keeps positions 2 and 3. Against window C positions [1, 2, 3, 4],
        // position 2 has window index 1 (delta 1) and position 3 has window
        // index 2 (delta 0).
        let m = reconstruct(
            &crate::mods::parse(b"C+m,0,1,0;", &[10, 20, 30]),
            b"CCCCC",
            1,
            5,
        );
        assert_eq!(m.groups[0].deltas, vec![1, 0]);
        assert_eq!(m.groups[0].ml, vec![20, 30]);
    }

    #[test]
    fn mod_exactly_at_window_start_is_kept() {
        // C at 0..=4, all modified. In window [2,5), position 2 (the start) and
        // position 4 (end - 1) survive; the half-open end excludes position 5.
        let m = recon(b"C+m,0,0,0,0,0;", &[1, 2, 3, 4, 5], b"CCCCC", 2, 5);
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].deltas, vec![0, 0, 0]); // abs 2,3,4 renumbered
        assert_eq!(m.groups[0].ml, vec![3, 4, 5]);
    }

    #[test]
    fn mod_exactly_at_window_end_is_excluded() {
        // Window [0,3) keeps position 0 (the start) and drops position 3 (the
        // exclusive end).
        let m = recon(b"C+m,0,0,0,0,0;", &[1, 2, 3, 4, 5], b"CCCCC", 0, 3);
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].deltas, vec![0, 0, 0]); // abs 0,1,2
        assert_eq!(m.groups[0].ml, vec![1, 2, 3]);
    }

    /// Random windows preserve the ML bytes associated with surviving positions.
    #[test]
    fn ml_stays_byte_aligned_over_random_windows() {
        // A deterministic LCG keeps the test reproducible without an external RNG.
        struct Lcg(u64);
        impl Lcg {
            fn next_u64(&mut self) -> u64 {
                self.0 = self
                    .0
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                self.0
            }
            fn below(&mut self, n: usize) -> usize {
                ((self.next_u64() >> 33) as usize) % n
            }
        }
        let mut rng = Lcg(0x0123_4567_89ab_cdef);

        for _ in 0..3000 {
            let n = 5 + rng.below(40);
            let seq: Vec<u8> = (0..n).map(|_| b"ACGT"[rng.below(4)]).collect();
            let c_pos: Vec<usize> = seq
                .iter()
                .enumerate()
                .filter(|&(_, &b)| b == b'C')
                .map(|(i, _)| i)
                .collect();
            if c_pos.is_empty() {
                continue;
            }
            // Random subset of C occurrences to modify (ascending), one ML byte each.
            let mut occ = Vec::new();
            for o in 0..c_pos.len() {
                if rng.below(2) == 0 {
                    occ.push(o);
                }
            }
            if occ.is_empty() {
                continue;
            }
            let mut deltas = Vec::new();
            let mut prev: i64 = -1;
            for &o in &occ {
                deltas.push((o as i64 - prev - 1) as usize);
                prev = o as i64;
            }
            let ml: Vec<u8> = (0..occ.len()).map(|_| rng.below(256) as u8).collect();
            let mut mm = b"C+m".to_vec();
            for d in &deltas {
                mm.extend_from_slice(format!(",{d}").as_bytes());
            }
            mm.push(b';');
            let a = rng.below(n + 1);
            let b = rng.below(n + 1);
            let (start, end) = if a <= b { (a, b) } else { (b, a) };

            let out = recon(&mm, &ml, &seq, start, end);

            // Independent expected survivors (position order preserved).
            let mut exp_ml = Vec::new();
            for (k, &o) in occ.iter().enumerate() {
                let abs = c_pos[o];
                if abs >= start && abs < end {
                    exp_ml.push(ml[k]);
                }
            }

            assert_eq!(
                out.groups.len(),
                1,
                "One group expected: seq={seq:?} window=[{start},{end})"
            );
            if exp_ml.is_empty() {
                assert!(
                    out.groups[0].deltas.is_empty() && out.groups[0].ml.is_empty(),
                    "No survivors in [{start},{end}) but positions were emitted: seq={seq:?}"
                );
            } else {
                let g = &out.groups[0];
                // With a single code there is exactly one ML byte per surviving delta.
                assert_eq!(
                    g.ml.len(),
                    g.deltas.len(),
                    "ML/deltas length mismatch: seq={seq:?} window=[{start},{end})"
                );
                assert_eq!(
                    g.ml, exp_ml,
                    "Kept ML must equal the in-window positions' bytes: seq={seq:?} window=[{start},{end})"
                );
            }
        }
    }

    /// A corrupt second delta clamps to `usize::MAX` in `parse`; the running
    /// cursor saturates rather than overflowing, and the unreachable position
    /// is dropped.
    #[test]
    fn huge_delta_does_not_overflow_cursor() {
        let m = recon(b"C+m,1,99999999999999999999;", &[5, 6], b"CCCC", 0, 4);
        // Only the first modified C (occurrence 1, position 1) survives.
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].deltas, vec![1]);
        assert_eq!(m.groups[0].ml, vec![5]);
    }

    /// `ml_stays_byte_aligned_over_random_windows` extended to multi-code groups
    /// (`ncodes` of 1 or 2, e.g. `C+m` and `C+mh`). ML is position-major, so
    /// after slicing its length equals `surviving_positions * ncodes` and its
    /// bytes are the in-window positions' `ncodes`-byte runs in order. The
    /// single-code test cannot reach this misalignment class.
    #[test]
    fn ml_stays_byte_aligned_multicode_over_random_windows() {
        struct Lcg(u64);
        impl Lcg {
            fn next_u64(&mut self) -> u64 {
                self.0 = self
                    .0
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                self.0
            }
            fn below(&mut self, n: usize) -> usize {
                ((self.next_u64() >> 33) as usize) % n
            }
        }
        let mut rng = Lcg(0xdead_beef_cafe_babe);

        for _ in 0..3000 {
            let ncodes = 1 + rng.below(2); // 1 or 2
            let codes: &[u8] = if ncodes == 1 { b"m" } else { b"mh" };

            let n = 5 + rng.below(40);
            let seq: Vec<u8> = (0..n).map(|_| b"ACGT"[rng.below(4)]).collect();
            let c_pos: Vec<usize> = seq
                .iter()
                .enumerate()
                .filter(|&(_, &b)| b == b'C')
                .map(|(i, _)| i)
                .collect();
            if c_pos.is_empty() {
                continue;
            }

            // Random subset of C occurrences to modify (ascending).
            let mut occ = Vec::new();
            for o in 0..c_pos.len() {
                if rng.below(2) == 0 {
                    occ.push(o);
                }
            }
            if occ.is_empty() {
                continue;
            }
            let mut deltas = Vec::new();
            let mut prev: i64 = -1;
            for &o in &occ {
                deltas.push((o as i64 - prev - 1) as usize);
                prev = o as i64;
            }
            // ML holds `ncodes` bytes per modified position, position-major.
            let ml: Vec<u8> = (0..occ.len() * ncodes)
                .map(|_| rng.below(256) as u8)
                .collect();

            let mut mm = b"C+".to_vec();
            mm.extend_from_slice(codes);
            for d in &deltas {
                mm.extend_from_slice(format!(",{d}").as_bytes());
            }
            mm.push(b';');

            let a = rng.below(n + 1);
            let b = rng.below(n + 1);
            let (start, end) = if a <= b { (a, b) } else { (b, a) };

            let out = recon(&mm, &ml, &seq, start, end);

            // Independent expected survivors: each in-window position's
            // `ncodes`-byte ML run, concatenated in order.
            let mut exp_ml = Vec::new();
            for (k, &o) in occ.iter().enumerate() {
                let abs = c_pos[o];
                if abs >= start && abs < end {
                    exp_ml.extend_from_slice(&ml[k * ncodes..k * ncodes + ncodes]);
                }
            }

            assert_eq!(
                out.groups.len(),
                1,
                "One group expected: seq={seq:?} window=[{start},{end}) ncodes={ncodes}"
            );
            if exp_ml.is_empty() {
                assert!(
                    out.groups[0].deltas.is_empty() && out.groups[0].ml.is_empty(),
                    "No survivors but positions were emitted: seq={seq:?} ncodes={ncodes}"
                );
            } else {
                let g = &out.groups[0];
                assert_eq!(g.codes.len(), ncodes, "Code count must be preserved");
                assert_eq!(
                    g.ml.len(),
                    g.deltas.len() * ncodes,
                    "ML must stay position-major: seq={seq:?} window=[{start},{end}) ncodes={ncodes}"
                );
                assert_eq!(
                    g.ml, exp_ml,
                    "Kept ML must equal the in-window positions' runs: seq={seq:?} window=[{start},{end}) ncodes={ncodes}"
                );
            }
        }
    }
}
