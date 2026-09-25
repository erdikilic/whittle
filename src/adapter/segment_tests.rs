use super::preset::preset_ont;
use super::*;

#[test]
fn interior_tolerance_depends_on_pattern_and_read_length() {
    let lsk = Budget::new(b"AATGTACTTCGTTCAGTTACGTATTGCT", 0.2, 150);
    assert_eq!(lsk.k_end, 5);
    assert_eq!(lsk.interior(2_000), 4);
    assert_eq!(lsk.interior(30_000), 4);
    assert_eq!(lsk.interior(200_000), 3);
    let short = Budget::new(b"TGGTTAGACTACGTATTGCTG", 0.2, 150);
    assert_eq!(short.interior(2_000), 2);
    assert_eq!(short.interior(30_000), 1);
    let ambiguous = Budget::new(&[b'N'; 40], 0.2, 150);
    assert_eq!(ambiguous.interior_max(), 0);
    for length in 11..=100 {
        let budget = Budget::new(&vec![b'A'; length], 0.2, 150);
        assert!(budget.interior_max() <= budget.k_end);
        assert!(budget.k_mid.windows(2).all(|pair| pair[0] >= pair[1]));
        assert_eq!(budget.interior(0), budget.interior_max());
        assert_eq!(
            budget.interior(usize::MAX),
            budget.k_mid[INTERIOR_CLASSES - 1]
        );
    }
}

#[test]
fn terminal_tolerance_shrinks_with_short_patterns() {
    assert_eq!(Budget::new(b"AAGAAAGTTGTCGGTGTCTTTGTG", 0.2, 150).k_end, 4);
    assert_eq!(Budget::new(b"AGAGTTTGATYMTGGCTCAG", 0.2, 150).k_end, 4);
    assert_eq!(Budget::new(b"TGGTTAGACTACGTAT", 0.2, 150).k_end, 3);
    assert_eq!(Budget::new(b"TGGTTAGACTACGTA", 0.2, 150).k_end, 2);
    assert_eq!(Budget::new(b"GCTTGGGTGTT", 0.2, 150).k_end, 1);
    assert_eq!(
        Budget::new(b"AATGTACTTCGTTCAGTTACGTATTGCT", 0.2, 150).k_end,
        5
    );
}

/// Builds a configuration at error rate 0.2 and `end_size` 20.
fn cfg(adapters: Vec<Adapter>, split: bool) -> AdapterConfig {
    AdapterConfig {
        adapters,
        error_rate: 0.2,
        end_size: 20,
        split,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    }
}

/// Builds a configuration from its parts with `min_piece` 1.
fn cfg_with(
    adapters: Vec<Adapter>,
    error_rate: f64,
    end_size: usize,
    split: bool,
) -> AdapterConfig {
    AdapterConfig {
        adapters,
        error_rate,
        end_size,
        split,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    }
}

/// Builds an adapter-role entry.
fn ad(name: &str, seq: &[u8]) -> Adapter {
    Adapter {
        name: name.into(),
        seq: seq.to_vec(),
        role: Role::Adapter,
    }
}

/// Builds an entry with an explicit role.
fn entry(name: &str, seq: &[u8], role: Role) -> Adapter {
    Adapter {
        name: name.into(),
        seq: seq.to_vec(),
        role,
    }
}

/// Generates deterministic SplitMix64 bases with the same generator the
/// inference fixtures use.
fn splitmix_dna(seed: u64, len: usize) -> Vec<u8> {
    let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(seed);
    (0..len)
        .map(|_| {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            b"ACGT"[((z >> 62) & 0b11) as usize]
        })
        .collect()
}

/// Returns `candidate_windows` grouped per adapter.
fn windows_by_adapter(index: &CandidateIndex, text: &[u8]) -> Vec<Vec<(usize, usize)>> {
    let mut windows = Vec::new();
    index.candidate_windows(text, &mut windows);
    let mut by_adapter = vec![Vec::new(); index.budgets.len()];
    for (adapter_idx, start, end) in windows {
        by_adapter[adapter_idx].push((start, end));
    }
    by_adapter
}

/// Computes the segments with every splitting adapter searched over the
/// whole window instead of its candidate windows, as the reference for the
/// seed filter. Every other pass is shared with `adapter_segments`.
fn reference_segments(window: &[u8], cfg: &AdapterConfig) -> Vec<(usize, usize)> {
    let mut index = CandidateIndex::new(&cfg.adapters, cfg.error_rate, cfg.end_size, cfg.split);
    index.seeds = None;
    for adapter_idx in 0..cfg.adapters.len() {
        index.unfiltered[adapter_idx] = cfg.split && index.split_classes[adapter_idx] > 0;
    }
    let exhaustive = AdapterConfig {
        candidate_index: std::sync::OnceLock::from(index),
        ..cfg.clone()
    };
    adapter_segments(window, &exhaustive)
}

/// Linear congruential generator for deterministic fixtures.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 32) as usize
    }
    fn below(&mut self, n: usize) -> usize {
        self.next() % n
    }
    fn dna(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[self.below(4)]).collect()
    }
}

/// Plants one adapter with up to the interior budget of edits at an end or in the
/// interior of a random window and checks the candidate search against the
/// full-window reference. `degenerate` rewrites a share of adapter
/// positions to ambiguity codes; the planted copy then carries one base
/// each code stands for.
fn check_candidate_search_randomized(seed: u64, degenerate: bool) {
    const CODES: &[u8] = b"RYSWKMBDHVN";
    let mut rng = Lcg(seed);
    for case in 0..400 {
        let adapters: Vec<Adapter> = (0..(1 + rng.below(10)))
            .map(|i| {
                let len = 11 + rng.below(40);
                let mut seq = rng.dna(len);
                if degenerate {
                    for _ in 0..(1 + len / 8) {
                        let p = rng.below(len);
                        seq[p] = CODES[rng.below(CODES.len())];
                    }
                }
                let role = match rng.below(4) {
                    0 => Role::Primer,
                    1 => Role::Barcode,
                    _ => Role::Adapter,
                };
                entry(&format!("a{i}"), &seq, role)
            })
            .collect();
        let window_len = 80 + rng.below(660);
        let mut window = rng.dna(window_len);

        let planted = rng.below(adapters.len());
        let pattern = &adapters[planted].seq;
        if pattern.len() <= window.len() {
            let max_edits = (0.1 * pattern.len() as f64).floor() as usize;
            let mut copy: Vec<u8> = pattern
                .iter()
                .map(|&code| {
                    let bases = iupac_bases(code).expect("adapter bytes are nucleotide codes");
                    bases[rng.below(bases.len())]
                })
                .collect();
            for _ in 0..rng.below(max_edits + 1) {
                match rng.below(3) {
                    0 => {
                        let p = rng.below(copy.len());
                        let old = copy[p];
                        copy[p] =
                            b"ACGT"[(b"ACGT".iter().position(|&b| b == old).unwrap() + 1) % 4];
                    },
                    1 => {
                        let p = rng.below(copy.len() + 1);
                        copy.insert(p, b"ACGT"[rng.below(4)]);
                    },
                    _ => {
                        let p = rng.below(copy.len());
                        copy.remove(p);
                    },
                }
            }
            let planted_len = copy.len();
            let pos = match rng.below(3) {
                0 => rng.below(8.min(window.len() - planted_len + 1)),
                1 => window.len() - planted_len - rng.below(8.min(window.len() - planted_len + 1)),
                _ => rng.below(window.len() - planted_len + 1),
            };
            window[pos..pos + planted_len].copy_from_slice(&copy);
            if case % 7 == 0 {
                window.make_ascii_lowercase();
            }
        }

        let cfg = AdapterConfig {
            adapters,
            error_rate: 0.2,
            end_size: 1 + rng.below(180),
            split: true,
            min_piece: 1 + rng.below(60),
            candidate_index: std::sync::OnceLock::new(),
        };
        assert_eq!(
            adapter_segments(&window, &cfg),
            reference_segments(&window, &cfg),
            "Candidate/reference mismatch in randomized case {case} (degenerate: {degenerate})"
        );
    }
}

#[test]
fn unknown_read_bases_consume_adapter_error_budget() {
    let adapter = b"ACGTAAGTCAGTACGATCAG";
    let mut text = adapter.to_vec();
    text[5] = b'N';
    text.extend_from_slice(&[b'C'; 80]);
    let exact = cfg_with(vec![ad("a", adapter)], 0.0, 8, true);
    assert_eq!(adapter_segments(&text, &exact), vec![(0, text.len())]);
    let tolerant = cfg_with(vec![ad("a", adapter)], 0.1, 8, true);
    assert_ne!(adapter_segments(&text, &tolerant), vec![(0, text.len())]);
    let mut unknown = vec![b'N'; 40];
    unknown.extend_from_slice(&[b'C'; 80]);
    let homopolymer = cfg_with(vec![ad("poly_a", &[b'A'; 20])], 0.0, 8, true);
    assert_eq!(
        adapter_segments(&unknown, &homopolymer),
        vec![(0, unknown.len())]
    );
}

/// The candidate search matches the full-window reference on random plain
/// adapters.
#[test]
fn candidate_search_matches_full_search_randomized() {
    check_candidate_search_randomized(0x4e4f_4f44_4c45_5301, false);
}

/// The candidate search matches the full-window reference on random
/// degenerate adapters.
#[test]
fn candidate_search_matches_full_search_randomized_degenerate() {
    check_candidate_search_randomized(0x4445_4745_4e45_5241, true);
}

/// A window holding an `N` run matches the reference after normalization.
#[test]
fn non_acgt_window_falls_back_to_scalar_search() {
    let cfg = cfg(
        vec![ad("a", b"ACGTACGTACGT"), ad("b", b"TTTTGGGGCCCC")],
        true,
    );
    let window = b"ACGTACGTACGTNNNNNNNNNNNNNNNNNNNNTTTTGGGGCCCC";
    assert_eq!(
        adapter_segments(window, &cfg),
        reference_segments(window, &cfg)
    );
}

/// A copy within the edit budget always retains one exact seed piece.
#[test]
fn partition_seeds_survive_random_indels_and_substitutions() {
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> usize {
            self.0 = self
                .0
                .wrapping_mul(2862933555777941757)
                .wrapping_add(3037000493);
            (self.0 >> 32) as usize
        }
        fn below(&mut self, n: usize) -> usize {
            self.next() % n
        }
        fn base(&mut self) -> u8 {
            b"ACGT"[self.below(4)]
        }
    }

    let mut rng = Lcg(0x5049_4745_4f4e_484f);
    for case in 0..1000 {
        let pattern: Vec<u8> = (0..(11 + rng.below(50))).map(|_| rng.base()).collect();
        let k = Budget::new(&pattern, 0.2, 150).interior_max();
        let mut mutated = pattern.clone();
        for _ in 0..rng.below(k + 1) {
            match rng.below(3) {
                0 => {
                    let p = rng.below(mutated.len());
                    mutated[p] = rng.base();
                },
                1 => {
                    let p = rng.below(mutated.len() + 1);
                    mutated.insert(p, rng.base());
                },
                _ if mutated.len() > 1 => {
                    let p = rng.below(mutated.len());
                    mutated.remove(p);
                },
                _ => {},
            }
        }
        let index = CandidateIndex::new(&[ad("a", &pattern)], 0.2, 150, true);
        let mut text: Vec<u8> = (0..17).map(|_| rng.base()).collect();
        text.extend_from_slice(&mutated);
        text.extend((0..19).map(|_| rng.base()));
        if case % 2 == 0 {
            text.make_ascii_lowercase();
        }
        let mut buf = Vec::new();
        assert!(
            !windows_by_adapter(&index, normalize_into(&text, &mut buf).0)[0].is_empty(),
            "Lossless seed filter rejected <=k edit case {case}"
        );
    }
}

/// An exact partial hit always keeps an intact `END_SEED_LEN`-mer, so the
/// end-seed gate never drops one; a hit with edits keeps one whenever its
/// longest intact stretch reaches the seed length.
#[test]
fn end_seeds_survive_partial_hits_within_budget() {
    let mut rng = Lcg(0x454e_4453_4545_4453);
    for case in 0..2000 {
        let len = 11 + rng.below(45);
        let pattern = rng.dna(len);
        let overlap = MIN_OVERLAP + rng.below(len - MIN_OVERLAP + 1);
        let budget = partial_budget(0.2, overlap);
        let mut copy = pattern[..overlap].to_vec();
        let mut edited = std::collections::BTreeSet::new();
        for _ in 0..budget {
            let p = rng.below(copy.len());
            copy[p] = b"ACGT"[(b"ACGT".iter().position(|&b| b == copy[p]).unwrap() + 1) % 4];
            edited.insert(p);
        }
        let mut longest = 0;
        let mut run = 0;
        for p in 0..overlap {
            run = if edited.contains(&p) { 0 } else { run + 1 };
            longest = longest.max(run);
        }
        if longest < END_SEED_LEN {
            continue;
        }
        let index = CandidateIndex::new(&[ad("a", &pattern)], 0.2, 150, true);
        let mut text = rng.dna(60);
        text.extend_from_slice(&copy);
        let (mut head, mut tail) = (Vec::new(), Vec::new());
        index.end_candidates(&text, 0, text.len(), 150, &mut head, &mut tail);
        assert!(tail[0], "End seed missing in case {case}");
    }
}

/// `expand_iupac` enumerates every concrete string and refuses to expand
/// past the cap or over a non-nucleotide byte.
#[test]
fn expand_iupac_enumerates_every_base_and_caps() {
    assert_eq!(expand_iupac(b"AC"), Some(vec![b"AC".to_vec()]));
    let mut r = expand_iupac(b"RY").unwrap();
    r.sort();
    assert_eq!(
        r,
        vec![
            b"AC".to_vec(),
            b"AT".to_vec(),
            b"GC".to_vec(),
            b"GT".to_vec()
        ]
    );
    assert_eq!(expand_iupac(b"NNNN").map(|v| v.len()), Some(256));
    assert_eq!(expand_iupac(b"NNNNN"), None, "Five N's expand past the cap");
    assert_eq!(
        expand_iupac(b"ACXT"),
        None,
        "A non-nucleotide byte has no expansion"
    );
}

/// Reverse complementing twice is the identity for every IUPAC code, with
/// case preserved.
#[test]
fn reverse_complement_uses_the_full_iupac_table() {
    for &code in b"ACGTRYSWKMBDHVNacgtryswkmbdhvn" {
        let once = reverse_complement(&[code]);
        assert_eq!(
            reverse_complement(&once),
            vec![code],
            "Code {}",
            code as char
        );
        assert_eq!(
            once[0].is_ascii_lowercase(),
            code.is_ascii_lowercase(),
            "Case is preserved for {}",
            code as char
        );
    }
    assert_eq!(reverse_complement(b"RYKMBDHV"), b"BDHVKMRY");
    assert_eq!(reverse_complement(b"SWN"), b"NWS");
    assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
    assert_eq!(reverse_complement(b"AACG"), b"CGTT");
}

/// `edit_budget` does not round an integral product down through
/// floating-point error, and `partial_budget` allows no edit at the
/// shortest overlaps.
/// More distinct seed prefixes than a `u16` slot numbers shorten the
/// prefix instead of failing, and every seed is still found.
#[test]
fn seed_table_shortens_the_prefix_past_u16_slots() {
    let mut seeds: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
    for i in 0..70_000u64 {
        let seed: Vec<u8> = (0..12)
            .map(|j| b"ACGT"[((i >> (2 * j)) & 3) as usize])
            .collect();
        seeds.insert(seed, vec![0]);
    }
    let probe = seeds.keys().nth(69_999).unwrap().clone();
    let table = SeedTable::new(seeds).expect("seeds are present");
    assert!(table.prefix < MAX_SEED_PREFIX_LEN);
    let mut text = splitmix_dna(4242, 40);
    text.extend_from_slice(&probe);
    let mut found = false;
    table.scan(&text, |_, (start, _)| found |= start == 40);
    assert!(found);
}

/// The overhang discount is the cost the overhang searcher charged: each
/// overhanging side is floored on its own, at the searcher's `f32` rate.
#[test]
fn overhang_discount_matches_the_charged_cost() {
    assert_eq!(overhang_cost(0.2, 3, 3), 0);
    assert_eq!(overhang_cost(0.2, 5, 0), 1);
    assert_eq!(overhang_cost(0.2, 5, 5), 2);
    assert_eq!(overhang_cost(0.2, 4, 6), 1);
}

#[test]
fn budgets_are_integral_and_partial_budget_is_strict() {
    assert_eq!(edit_budget(0.29, 100), 29);
    assert_eq!(edit_budget(0.57, 100), 57);
    assert_eq!(edit_budget(0.2, 22), 4);
    assert_eq!(edit_budget(0.1, 12), 1);
    assert_eq!(edit_budget(0.0, 50), 0);
    assert_eq!(partial_budget(0.2, MIN_OVERLAP), 0);
    assert_eq!(partial_budget(0.2, 13), 0);
    assert_eq!(partial_budget(0.2, 14), 1);
    assert_eq!(partial_budget(0.2, 19), 2);
    assert_eq!(partial_budget(0.2, 28), 3);
}

/// Both seed pieces of the 12-mer hold an `N`; the planted copy is one
/// concrete instance and still splits the read.
#[test]
fn degenerate_adapter_splits_interior_chimera() {
    let adapter = b"GTNGTTGGNTGT";
    let mut w = vec![b'A'; 40];
    w.extend_from_slice(b"GTGGTTGGGTGT");
    w.extend_from_slice(&[b'C'; 40]);
    let c = cfg_with(vec![ad("deg", adapter)], 0.2, 10, true);
    assert_eq!(adapter_segments(&w, &c), vec![(0, 40), (52, 92)]);
}

/// The `N` sits in the first piece and the substitution (C to A) in the
/// second, so only an expanded first-piece seed finds the copy.
#[test]
fn degenerate_adapter_splits_with_one_substitution_in_the_plain_piece() {
    let adapter = b"GTNGTTGGCTGTACCGATCA";
    let mut w = vec![b'A'; 40];
    w.extend_from_slice(b"GTGGTTGGATGTACCGATCA");
    w.extend_from_slice(&[b'C'; 40]);
    let c = cfg_with(vec![ad("deg", adapter)], 0.2, 10, true);
    assert_eq!(adapter_segments(&w, &c), vec![(0, 40), (60, 100)]);
}

/// Each 8-base piece holds five `N`s (1024 expansions), past the cap, so
/// the adapter is searched over the whole window and still splits it.
#[test]
fn overly_degenerate_adapter_is_searched_unfiltered() {
    let adapter = b"NNNNNGGTTGGNNNNN";
    let mut w = vec![b'A'; 35];
    w.extend_from_slice(b"CACGTGGTTGGACGTC");
    w.extend_from_slice(&[b'C'; 41]);
    let c = cfg_with(vec![ad("deg", adapter)], 0.2, 10, true);
    let index = CandidateIndex::new(&c.adapters, c.error_rate, c.end_size, true);
    assert_eq!(index.unfiltered, vec![true]);
    assert!(
        index.seeds.is_none(),
        "No seed is built for an unfiltered adapter"
    );
    assert_eq!(windows_by_adapter(&index, &w), vec![vec![(0, w.len())]]);
    let segs = adapter_segments(&w, &c);
    assert_eq!(segs, reference_segments(&w, &c));
    assert_eq!(
        segs.len(),
        2,
        "The unfiltered adapter still splits the read"
    );
}

/// An empty adapter set keeps the whole window.
#[test]
fn no_adapters_is_identity() {
    let w = b"ACGTACGTACGTACGT";
    assert_eq!(adapter_segments(w, &cfg(vec![], true)), vec![(0, w.len())]);
}

/// A 5' adapter at the read start is trimmed.
#[test]
fn trims_5prime_adapter_and_outboard() {
    let adapter = b"ACGTACGTACGT";
    let mut w = adapter.to_vec();
    w.extend_from_slice(b"AAAAAAAAAAAA");
    let c = cfg(vec![ad("a", adapter)], false);
    assert_eq!(adapter_segments(&w, &c), vec![(12, 24)]);
}

/// Bases before a 5' adapter are trimmed with it, whichever catalog entry
/// matched: the reverse complement of a rear entry at the read start is
/// a front adapter.
#[test]
fn leading_junk_is_trimmed_through_a_reverse_complement_hit() {
    let front = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let rear = reverse_complement(front);
    let mut w = b"TGACCGTAG".to_vec();
    w.extend_from_slice(front);
    w.extend(splitmix_dna(1, 300));
    let c = cfg_with(vec![ad("rear", &rear)], 0.2, 150, true);
    assert_eq!(adapter_segments(&w, &c), vec![(9 + front.len(), w.len())]);
}

/// A longer entry sharing its core with the primer in the read, its extra
/// bases matched as edits against the insert, trims at the end of the core.
#[test]
fn extended_entry_trims_at_the_end_of_its_shared_core() {
    let primer = b"TTTCTGTTGGTGCTGATATTGC";
    let extended = b"TTTCTGTTGGTGCTGATATTGCTTT";
    let mut w = primer.to_vec();
    w.extend(splitmix_dna(7, 300));
    w.extend(reverse_complement(primer));
    let c = cfg_with(
        vec![
            entry("primer", primer, Role::Primer),
            entry("extended", extended, Role::Primer),
        ],
        0.2,
        150,
        true,
    );
    assert_eq!(
        adapter_segments(&w, &c),
        vec![(primer.len(), w.len() - primer.len())]
    );
}

/// A rear adapter cut short by the read end is trimmed from its first
/// aligned base when at least `MIN_OVERLAP` bases align, and left in place
/// below that.
#[test]
fn truncated_rear_adapter_at_the_read_end_is_trimmed() {
    let front = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let rear = reverse_complement(front);
    let insert = splitmix_dna(2, 400);
    let c = cfg_with(vec![ad("front", front)], 0.2, 150, true);
    for remnant in [MIN_OVERLAP, 14, 20, 27] {
        let mut w = insert.clone();
        w.extend_from_slice(&rear[..remnant]);
        assert_eq!(
            adapter_segments(&w, &c),
            vec![(0, insert.len())],
            "Remnant of {remnant} bases"
        );
    }
    let mut w = insert.clone();
    w.extend_from_slice(&rear[..MIN_OVERLAP - 1]);
    assert_eq!(adapter_segments(&w, &c), vec![(0, w.len())]);
}

/// A front adapter missing its leading bases is trimmed from the read
/// start, with the tolerance `partial_budget` allows for the aligned part.
#[test]
fn truncated_front_adapter_at_the_read_start_is_trimmed() {
    let front = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let insert = splitmix_dna(3, 400);
    let c = cfg_with(vec![ad("front", front)], 0.2, 150, true);
    for missing in [1, 6, 12, front.len() - MIN_OVERLAP] {
        let mut w = front[missing..].to_vec();
        w.extend_from_slice(&insert);
        assert_eq!(
            adapter_segments(&w, &c),
            vec![(front.len() - missing, w.len())],
            "Missing {missing} leading bases"
        );
    }
    // One substitution inside a 20-base remnant is within budget; inside
    // a 12-base remnant it is not.
    let mut long = front[8..].to_vec();
    long[5] = if long[5] == b'A' { b'C' } else { b'A' };
    long.extend_from_slice(&insert);
    assert_eq!(adapter_segments(&long, &c), vec![(20, long.len())]);
    let mut short = front[16..].to_vec();
    short[5] = if short[5] == b'A' { b'C' } else { b'A' };
    short.extend_from_slice(&insert);
    assert_eq!(adapter_segments(&short, &c), vec![(0, short.len())]);
}

/// An adapter prefix inside the read, away from either end, is not a
/// partial hit: the overhang search accepts a partial alignment only flush
/// with the read end it hangs off.
#[test]
fn interior_adapter_prefix_is_not_a_partial_hit() {
    let front = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let mut w = splitmix_dna(4, 60);
    w.extend_from_slice(&front[..14]);
    w.extend(splitmix_dna(5, 300));
    let c = cfg_with(vec![ad("front", front)], 0.2, 150, true);
    assert_eq!(adapter_segments(&w, &c), vec![(0, w.len())]);
}

/// Every role trims the read ends and splits at an interior hit.
#[test]
fn every_role_trims_ends_and_splits() {
    let seq = b"GGGGTTTTGGGGTTTTGGGG";
    let mut w = seq.to_vec();
    w.extend_from_slice(&[b'A'; 60]);
    w.extend_from_slice(seq);
    w.extend_from_slice(&[b'C'; 60]);
    w.extend_from_slice(seq);
    for role in [Role::Adapter, Role::Primer, Role::Barcode] {
        let c = cfg_with(vec![entry("p", seq, role)], 0.2, 30, true);
        assert_eq!(
            adapter_segments(&w, &c),
            vec![(20, 80), (100, 160)],
            "{role:?}"
        );
    }
}

/// A panel barcode splits a read when the set has no barcode flanks, and
/// leaves the split to the flanks when the set carries them.
#[test]
fn panel_barcodes_split_only_without_flanks() {
    let panel = [
        entry("bc1", b"AAGAAAGTTGTCGGTGTCTTTGTG", Role::Barcode),
        entry("bc2", b"TCGATTCCGTTTGTAGTCGTCTGT", Role::Barcode),
    ];
    let flank = entry("rear", b"TTAACCTTTCTGTTGGTGCTGATATTGC", Role::Barcode);
    let bare = CandidateIndex::new(&panel, 0.2, 150, true);
    assert!(bare.split_classes.iter().all(|&c| c > 0));
    let mut flanked = panel.to_vec();
    flanked.push(flank);
    let index = CandidateIndex::new(&flanked, 0.2, 150, true);
    assert_eq!(index.split_classes[..2], [0, 0]);
    assert!(index.split_classes[2] > 0);
}

/// The piece on each side of an excision is trimmed at its new end: a
/// primer left before the junction adapter and one after it are removed
/// with the split.
#[test]
fn split_pieces_are_trimmed_at_the_junction() {
    let adapter = b"GGGGTTTTGGGGTTTTGGGG";
    let primer = b"CACACAGAGAGACACACAGAGA";
    let primer_rc = reverse_complement(primer);
    let mut w = vec![b'A'; 80];
    w.extend_from_slice(&primer_rc);
    w.extend_from_slice(adapter);
    w.extend_from_slice(primer);
    w.extend_from_slice(&[b'C'; 80]);
    let c = cfg_with(
        vec![ad("a", adapter), entry("p", primer, Role::Primer)],
        0.2,
        30,
        true,
    );
    assert_eq!(adapter_segments(&w, &c), vec![(0, 80), (144, 224)]);
}

/// A truncated rear adapter followed by unalignable bases before the
/// junction adapter is trimmed by the residue search of the left piece.
#[test]
fn junction_residue_behind_junk_is_trimmed() {
    let front = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let rear = reverse_complement(front);
    let insert1 = splitmix_dna(6, 300);
    let insert2 = splitmix_dna(7, 300);
    let junk = splitmix_dna(8, 20);
    let mut w = insert1.clone();
    w.extend_from_slice(&rear[..16]);
    w.extend_from_slice(&junk);
    w.extend_from_slice(front);
    w.extend_from_slice(&insert2);
    let c = cfg_with(vec![ad("front", front)], 0.2, 150, true);
    let cut = insert1.len() + 16 + junk.len() + front.len();
    assert_eq!(
        adapter_segments(&w, &c),
        vec![(0, insert1.len()), (cut, w.len())]
    );
}

/// With the default `end_size` of 150 the two end zones overlap for any
/// read of at most 300 bp. A chimera-junction adapter within `end_size` of
/// both ends splits the read and keeps both inserts; treating it as a
/// terminal adapter would discard the entire outboard arm, up to `end_size`
/// bases of real insert.
#[test]
fn central_chimera_on_short_read_splits_both_arms() {
    let adapter = b"GGGGTTTTGGGGTTTTGGGG";
    let mut w = vec![b'A'; 115];
    let cut = w.len();
    w.extend_from_slice(adapter);
    w.extend_from_slice(&[b'C'; 115]);
    let c = cfg_with(vec![ad("mid", adapter)], 0.2, 150, true);
    let segs = adapter_segments(&w, &c);
    assert_eq!(
        segs,
        vec![(0, cut), (cut + adapter.len(), w.len())],
        "Central chimera must split into both arms, not lose insert1"
    );
}

/// Three junk bases, the adapter, then the insert: the 3 bp flank is
/// trimmed with the adapter instead of surviving as its own segment.
#[test]
fn near_terminal_excision_folds_into_a_trim() {
    let adapter = b"GGGGTTTTGGGGTTTTGGGG";
    let mut w = b"AAA".to_vec();
    w.extend_from_slice(adapter);
    w.extend_from_slice(&[b'C'; 37]);
    let c = cfg_with(vec![ad("a", adapter)], 0.2, 150, true);
    assert_eq!(adapter_segments(&w, &c), vec![(23, 60)]);

    // Mirror: insert, adapter, three junk bases.
    let mut w = vec![b'C'; 37];
    w.extend_from_slice(adapter);
    w.extend_from_slice(b"AAA");
    assert_eq!(adapter_segments(&w, &c), vec![(0, 37)]);
}

/// Asserts one kept segment that reaches the far end, with the trim boundary
/// at the adapter or at most `slack` bases into the insert.
fn assert_insert_kept(
    segs: &[(usize, usize)],
    n: usize,
    adapter: usize,
    slack: usize,
    front: bool,
) {
    assert_eq!(segs.len(), 1, "Read length {n}: {segs:?}");
    let (start, end) = segs[0];
    if front {
        assert_eq!(end, n, "Read length {n}: {segs:?}");
        assert!(
            (adapter..=adapter + slack).contains(&start),
            "Read length {n}: {segs:?}"
        );
    } else {
        assert_eq!(start, 0, "Read length {n}: {segs:?}");
        let insert = n - adapter;
        assert!(
            (insert - slack..=insert).contains(&end),
            "Read length {n}: {segs:?}"
        );
    }
}

/// The catalog holds longer entries that begin with `PCR1_front`
/// (`cDNA_rear` adds `TTT`); their extra bases do not carry the trim into the
/// insert. The same insert padded to 200 bp, where the end zones do not
/// overlap, keeps the insert as well.
#[test]
fn short_read_with_front_adapter_keeps_insert_under_ont_preset() {
    let mut short = b"ACTTGCCTGTCGCTCTATCTTC".to_vec();
    short.extend_from_slice(b"GGGG");
    short.extend(splitmix_dna(0, 136));
    let mut long = short.clone();
    long.extend(splitmix_dna(5, 38));
    for split in [true, false] {
        let c = cfg_with(preset_ont(), 0.2, 150, split);
        assert_insert_kept(&adapter_segments(&short, &c), 162, 22, 0, true);
        assert_insert_kept(&adapter_segments(&long, &c), 200, 22, 0, true);
    }
}

/// Mirror: `cDNA_front` and `PCS110_front` end in `PCR2_front`, so their
/// reverse complements are `PCR2_rear` with three or four leading bases.
#[test]
fn short_read_with_rear_adapter_keeps_insert_under_ont_preset() {
    let mut short = splitmix_dna(0, 136);
    short.extend_from_slice(b"GGGG");
    short.extend_from_slice(b"GCAATATCAGCACCAACAGAAA");
    let mut long = splitmix_dna(5, 38);
    long.extend_from_slice(&short);
    for split in [true, false] {
        let c = cfg_with(preset_ont(), 0.2, 150, split);
        assert_insert_kept(&adapter_segments(&short, &c), 162, 22, 0, false);
        assert_insert_kept(&adapter_segments(&long, &c), 200, 22, 0, false);
    }
}

/// A native-barcoded read end is trimmed through the 8 bp inner flank that
/// follows the barcode, at either end.
#[test]
fn native_barcode_construct_is_trimmed_through_its_inner_flank() {
    let mut construct = b"ATTGCTAAGGTTAA".to_vec();
    construct.extend(reverse_complement(b"AAGAAAGTTGTCGGTGTCTTTGTG"));
    construct.extend_from_slice(b"CAGCACCT");
    let insert = splitmix_dna(3, 400);
    let mut w = construct.clone();
    w.extend_from_slice(&insert);
    w.extend(reverse_complement(&construct));
    let c = cfg_with(
        super::preset::preset(&[super::preset::Kit::Nbd114]),
        0.2,
        150,
        true,
    );
    assert_eq!(
        adapter_segments(&w, &c),
        vec![(construct.len(), construct.len() + insert.len())]
    );
}

/// A minimal user FASTA: `f` and its reverse complement, with the adapter
/// and its insert on a 162 bp read. The shortened rear entry exercises
/// two entries of different lengths.
#[test]
fn front_and_reverse_complement_rear_pair_keep_insert_on_short_read() {
    let f = b"ACTTGCCTGTCGCTCTATCTTC";
    let r = reverse_complement(f);
    let mut w = f.to_vec();
    w.extend(splitmix_dna(3, 140));
    for rear in [r.as_slice(), &r[1..]] {
        for split in [true, false] {
            let c = cfg_with(vec![ad("f", f), ad("r", rear)], 0.2, 150, split);
            assert_eq!(
                adapter_segments(&w, &c),
                vec![(22, 162)],
                "Split mode {split}, rear length {}",
                rear.len()
            );
        }
    }
}

/// Classification is geometric: inside the overlap of both end zones the
/// flank slack decides between a trim and an excision, and with one zone
/// only the zone decides.
#[test]
fn classify_terminal_is_geometric() {
    assert_eq!(classify_terminal(0, 20, 60, 60), Terminal::Five);
    assert_eq!(classify_terminal(40, 60, 60, 60), Terminal::Three);
    assert_eq!(classify_terminal(20, 40, 60, 60), Terminal::Excise);
    assert_eq!(classify_terminal(11, 31, 60, 60), Terminal::Five);
    assert_eq!(classify_terminal(12, 32, 60, 60), Terminal::Excise);
    assert_eq!(classify_terminal(29, 49, 60, 60), Terminal::Three);
    // Both flanks within slack: the nearer end.
    assert_eq!(classify_terminal(8, 30, 35, 35), Terminal::Three);
    assert_eq!(classify_terminal(4, 30, 35, 35), Terminal::Five);
    // One zone only.
    assert_eq!(classify_terminal(0, 20, 400, 150), Terminal::Five);
    assert_eq!(classify_terminal(380, 400, 400, 150), Terminal::Three);
    assert_eq!(classify_terminal(200, 220, 400, 150), Terminal::None);
}

/// An interior adapter splits the read into its two flanks.
#[test]
fn splits_on_interior_adapter() {
    let adapter = b"GGGGTTTTGGGGTTTT";
    let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
    let cut_start = w.len();
    w.extend_from_slice(adapter);
    w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC");
    let c = cfg(vec![ad("mid", adapter)], true);
    let segs = adapter_segments(&w, &c);
    assert_eq!(segs.len(), 2, "Interior adapter splits the read");
    assert_eq!(segs[0], (0, cut_start));
    assert_eq!(segs[1], (cut_start + adapter.len(), w.len()));
}

/// Ends-only mode leaves an interior adapter in place.
#[test]
fn ends_only_suppresses_interior_split() {
    let adapter = b"GGGGTTTTGGGGTTTT";
    let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
    w.extend_from_slice(adapter);
    w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC");
    let c = cfg(vec![ad("mid", adapter)], false);
    assert_eq!(adapter_segments(&w, &c), vec![(0, w.len())]);
}

/// A 5' adapter, an insert and a 3' adapter in ends-only mode: both ends
/// trim to the insert although only the two end zones are searched.
#[test]
fn ends_only_trims_both_terminal_adapters() {
    let adapter5 = b"ACGTACGTACGT";
    let adapter3 = b"TTTTGGGGCCCC";
    let insert = b"AAAAAAAAAAAA";
    let mut w = adapter5.to_vec();
    w.extend_from_slice(insert);
    w.extend_from_slice(adapter3);
    let c = cfg(vec![ad("five", adapter5), ad("three", adapter3)], false);
    assert_eq!(adapter_segments(&w, &c), vec![(12, 24)]);
}

/// A terminal 5' adapter that starts inside `end_size` but ends beyond it:
/// with `end_size` 4, a 12 bp adapter at position 2 spans [2, 14). A head
/// zone of `window[..end_size]` (4 bytes) cannot contain a 12-byte match;
/// the `end_size + len` sizing gives `window[..16]`, which does.
#[test]
fn ends_only_trims_adapter_straddling_end_size() {
    let adapter = b"ACGTACGTACGT";
    let mut w = b"AA".to_vec();
    w.extend_from_slice(adapter);
    w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCC");
    let c = cfg_with(vec![ad("five", adapter)], 0.2, 4, false);
    assert_eq!(adapter_segments(&w, &c), vec![(14, w.len())]);
}

/// A pattern below `MIN_PATTERN_LEN` is never searched.
#[test]
fn short_pattern_is_skipped() {
    let short = b"GGTGCTG";
    let w = b"GGTGCTGAAAAAAAAAAAAAAAA";
    let c = cfg(vec![ad("flank", short)], true);
    assert_eq!(adapter_segments(w, &c), vec![(0, w.len())]);
}

/// An empty window yields no segments.
#[test]
fn empty_window_returns_empty() {
    let c = cfg(vec![ad("a", b"ACGTACGTACGT")], true);
    assert_eq!(adapter_segments(b"", &c), vec![]);
}

/// The window is the adapter: the 5' trim advances `lo` to `n`, so
/// `lo >= hi` and the whole window is consumed.
#[test]
fn whole_window_consumed_returns_empty() {
    let adapter = b"ACGTACGTACGT";
    let c = cfg(vec![ad("a", adapter)], true);
    assert_eq!(adapter_segments(adapter, &c), vec![]);
}

/// Mirror of `trims_5prime_adapter_and_outboard` with the adapter at the
/// 3' end: insert first, adapter last.
#[test]
fn trims_3prime_adapter() {
    let adapter = b"ACGTACGTACGT";
    let mut w = b"AAAAAAAAAAAA".to_vec();
    w.extend_from_slice(adapter);
    let c = cfg(vec![ad("a", adapter)], false);
    assert_eq!(adapter_segments(&w, &c), vec![(0, 12)]);
}

/// Two distinct interior adapters whose hits overlap by 6 bp: `a` matches
/// [24, 40) and `b` matches [34, 50), constructed so their shared 6 bp
/// region is the same window bytes, giving both an exact hit. The overlap
/// merges into one excision, leaving exactly 2 segments.
#[test]
fn overlapping_interior_cuts_merge() {
    let a = b"GGGGTTTTTGTGTGTG";
    let b = b"TGTGTGTGTTTTGGGG";
    let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
    w.extend_from_slice(a);
    w.extend_from_slice(&b[6..]);
    w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC");
    let c = cfg(vec![ad("a", a), ad("b", b)], true);
    let segs = adapter_segments(&w, &c);
    assert_eq!(
        segs.len(),
        2,
        "Overlapping interior cuts merge into one excision"
    );
    assert_eq!(segs[0], (0, 24));
    assert_eq!(segs[1], (50, w.len()));
}

/// Two excisions separated by a gap within `FLANK_SLACK`, or shorter than
/// `min_piece`, merge into one; a longer gap stays a segment of its own.
#[test]
fn nearby_interior_cuts_merge_by_slack_and_min_piece() {
    let a = b"GGGGTTTTGGGGTTTTGGGG";
    let build = |gap: usize| {
        let mut w = vec![b'A'; 60];
        w.extend_from_slice(a);
        w.extend(splitmix_dna(9, gap));
        w.extend_from_slice(a);
        w.extend_from_slice(&[b'C'; 60]);
        w
    };
    let slack = build(FLANK_SLACK);
    let c = cfg_with(vec![ad("a", a)], 0.2, 30, true);
    assert_eq!(
        adapter_segments(&slack, &c),
        vec![(0, 60), (slack.len() - 60, slack.len())]
    );

    let wide = build(40);
    assert_eq!(
        adapter_segments(&wide, &c),
        vec![(0, 60), (80, 120), (wide.len() - 60, wide.len())]
    );
    let c = AdapterConfig {
        min_piece: 41,
        ..cfg_with(vec![ad("a", a)], 0.2, 30, true)
    };
    assert_eq!(
        adapter_segments(&wide, &c),
        vec![(0, 60), (wide.len() - 60, wide.len())]
    );
}

/// The terminal hit [0, 16) overlaps the interior hit [10, 26); clipping the
/// interior interval to the keep window still excises [16, 26).
#[test]
fn straddling_cut_is_clipped_not_leaked() {
    let t_prefix = b"GGTGTGGTTT";
    let overlap = b"GTTGGT";
    let s_suffix = b"TGGTGTTGGG";
    let mut t = t_prefix.to_vec();
    t.extend_from_slice(overlap);
    let mut s = overlap.to_vec();
    s.extend_from_slice(s_suffix);

    let mut w = t.clone();
    w.extend_from_slice(s_suffix);
    w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");

    let c = cfg_with(vec![ad("t", &t), ad("s", &s)], 0.2, 9, true);
    let segs = adapter_segments(&w, &c);
    assert_eq!(
        segs,
        vec![(26, 60)],
        "Adapter s is fully excised, no leaked bases before 26"
    );
    for &(seg_start, seg_end) in &segs {
        assert!(
            !w[seg_start..seg_end]
                .windows(s.len())
                .any(|win| win == s.as_slice())
        );
    }
}

/// A cost-4 hit within `k_end` of 6 but above the interior budget does
/// not split the read.
#[test]
fn interior_above_k_mid_does_not_split() {
    let adapter = b"GGTTGGTTGGTT";
    let mut mutated = adapter.to_vec();
    for &i in &[1usize, 4, 7, 10] {
        mutated[i] = match mutated[i] {
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        };
    }
    let mut w = b"AAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
    w.extend_from_slice(&mutated);
    w.extend_from_slice(b"CCCCCCCCCCCCCCCCCCCCCCCC");
    let c = cfg_with(vec![ad("mid", adapter)], 0.5, 10, true);
    assert_eq!(
        adapter_segments(&w, &c),
        vec![(0, w.len())],
        "Cost 4 hit is above the interior budget and must not split the read"
    );
}

/// A six-base insertion expands the terminal alignment from 40 to 46 bases.
/// The terminal search window includes `k_end` additional bases, so
/// ends-only and split modes select the same [2, 48) alignment.
#[test]
fn ends_only_equals_split_on_indel_terminal_adapter() {
    let adapter = b"AAAACCCCGGGGTTTTACGTTGCATCAGTCCAGTGACTGA";
    let extra = b"CTGACT";
    let mut copy = adapter[..20].to_vec();
    copy.extend_from_slice(extra);
    copy.extend_from_slice(&adapter[20..]);

    let mut w = b"AA".to_vec();
    w.extend_from_slice(&copy);
    w.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

    let c_split = cfg_with(vec![ad("five", adapter)], 0.15, 4, true);
    let c_ends_only = AdapterConfig {
        split: false,
        ..c_split.clone()
    };

    let split_segs = adapter_segments(&w, &c_split);
    let ends_only_segs = adapter_segments(&w, &c_ends_only);

    assert_eq!(
        split_segs,
        vec![(48, w.len())],
        "Split mode finds the full 46bp indel-bearing hit and trims to 48"
    );
    assert_eq!(
        ends_only_segs, split_segs,
        "Ends-only must match split mode exactly: the end zone must be wide \
             enough (end_size + len + k_end) to contain the full indel-lengthened hit"
    );
    assert_eq!(ends_only_segs[0].0, 48);
}

/// A 40 bp insert plus a 20 bp adapter at the 3' end with `end_size >= n`,
/// so both zones overlap. The insert [0, 40) is kept and the read is not
/// dropped.
#[test]
fn three_prime_adapter_on_short_read_trims_tail_not_whole_read() {
    let adapter = b"GGGGTTTTGGGGTTTTGGGG";
    let mut w = vec![b'A'; 40];
    w.extend_from_slice(adapter);
    let split = cfg_with(vec![ad("a", adapter)], 0.2, 150, true);
    let ends = AdapterConfig {
        split: false,
        ..split.clone()
    };
    assert_eq!(adapter_segments(&w, &split), vec![(0, 40)], "Split mode");
    assert_eq!(adapter_segments(&w, &ends), vec![(0, 40)], "Ends-only mode");
}

/// A 5' adapter on a short read trims the head and keeps the insert.
#[test]
fn five_prime_adapter_on_short_read_trims_head() {
    let adapter = b"GGGGTTTTGGGGTTTTGGGG";
    let mut w = adapter.to_vec();
    w.extend_from_slice(&[b'A'; 40]);
    let split = cfg_with(vec![ad("a", adapter)], 0.2, 150, true);
    assert_eq!(adapter_segments(&w, &split), vec![(20, 60)]);
}

/// The two terminal patterns and their reverse complements are distinct,
/// leaving the 40-base insert as the only retained segment.
#[test]
fn both_adapters_at_both_ends_keep_middle() {
    let a5 = b"GGGGTTTTGGGGTTTTGGGG";
    let a3 = b"AAAAGGGGAAAAGGGGAAAA";
    let mut w = a5.to_vec();
    w.extend_from_slice(&[b'T'; 40]);
    w.extend_from_slice(a3);
    let c = cfg_with(vec![ad("a5", a5), ad("a3", a3)], 0.2, 150, true);
    assert_eq!(adapter_segments(&w, &c), vec![(20, 60)]);
}

/// Returns a panel of `size` random barcodes of `len` bases.
fn barcode_panel(size: u64, len: usize) -> Vec<Adapter> {
    (1..=size)
        .map(|i| entry(&format!("bc{i}"), &splitmix_dna(i, len), Role::Barcode))
        .collect()
}

/// Returns `seq` with a different base at each of `positions`.
fn substituted(seq: &[u8], positions: &[usize]) -> Vec<u8> {
    let mut out = seq.to_vec();
    for &i in positions {
        out[i] = if out[i] == b'A' { b'C' } else { b'A' };
    }
    out
}

/// A 24-member panel of 16-base barcodes shares one chance bound: a hit
/// deep in the end zone is held to two edits, while a hit anchored at
/// the read end keeps the three edits a lone barcode of that length has.
#[test]
fn panel_members_share_one_terminal_chance_bound() {
    let panel = barcode_panel(24, 16);
    let index = CandidateIndex::new(&panel, 0.2, 150, true);
    for budget in &index.budgets {
        assert_eq!(budget.k_end, 3);
        assert_eq!(budget.k_far, 2);
    }
    let lone = CandidateIndex::new(&panel[..1], 0.2, 150, true);
    assert_eq!(lone.budgets[0].k_end, 3);
    assert_eq!(lone.budgets[0].k_far, 3);
}

/// Identical entries, and entries that are reverse complements of each
/// other, are one sequence to the chance bound.
#[test]
fn duplicate_entries_count_once_toward_the_chance_bound() {
    let panel = barcode_panel(24, 16);
    let mut doubled = panel.clone();
    doubled.extend(panel.iter().map(|a| {
        entry(
            &format!("{}_rc", a.name),
            &reverse_complement(&a.seq),
            a.role,
        )
    }));
    let single = CandidateIndex::new(&panel, 0.2, 150, true);
    let double = CandidateIndex::new(&doubled, 0.2, 150, true);
    for (a, b) in single.budgets.iter().zip(&double.budgets) {
        assert_eq!((a.k_end, a.k_far), (b.k_end, b.k_far));
    }
}

/// The whole catalog keeps the per-pattern terminal budget of every
/// adapter and primer at the read end, and deep in the end zone for every
/// entry of 24 bases or more.
#[test]
fn catalog_budgets_keep_their_per_pattern_tolerance() {
    let catalog = super::preset::preset(super::preset::Kit::ALL);
    let index = CandidateIndex::new(&catalog, 0.2, 150, true);
    for (adapter, budget) in catalog.iter().zip(&index.budgets) {
        let own = Budget::new(&adapter.seq, 0.2, 150);
        assert_eq!(budget.k_end, own.k_end, "{}", adapter.name);
        if adapter.seq.len() >= 24 {
            assert_eq!(budget.k_far, own.k_end, "{}", adapter.name);
        }
    }
}

/// A panel barcode with three substitutions trims at the read end, and
/// directly behind an adapter hit, but not 80 bases into the read, where
/// the panel admits two edits.
#[test]
fn marginal_panel_hit_trims_only_when_anchored() {
    let panel = barcode_panel(24, 16);
    let adapter = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let mut set = panel.clone();
    set.push(ad("lsk", adapter));
    let c = cfg_with(set, 0.2, 150, true);
    let insert = splitmix_dna(9001, 2000);
    let marginal = substituted(&panel[4].seq, &[3, 8, 13]);
    let n = insert.len() + marginal.len();

    assert_eq!(adapter_segments(&insert, &c), vec![(0, insert.len())]);

    let mut at_end = marginal.clone();
    at_end.extend_from_slice(&insert);
    assert_eq!(adapter_segments(&at_end, &c), vec![(16, n)]);

    let mut behind_adapter = adapter.to_vec();
    behind_adapter.extend_from_slice(&marginal);
    behind_adapter.extend_from_slice(&insert);
    assert_eq!(
        adapter_segments(&behind_adapter, &c),
        vec![(44, n + adapter.len())]
    );

    let mut deep = insert[..80].to_vec();
    deep.extend_from_slice(&marginal);
    deep.extend_from_slice(&insert[80..]);
    assert_eq!(adapter_segments(&deep, &c), vec![(0, n)]);

    let mut tail = insert.clone();
    tail.extend_from_slice(&reverse_complement(&marginal));
    assert_eq!(adapter_segments(&tail, &c), vec![(0, insert.len())]);

    let mut tail_deep = insert[..insert.len() - 80].to_vec();
    tail_deep.extend_from_slice(&reverse_complement(&marginal));
    tail_deep.extend_from_slice(&insert[insert.len() - 80..]);
    assert_eq!(adapter_segments(&tail_deep, &c), vec![(0, n)]);
}

/// A panel barcode within two edits trims wherever it lies in the end zone.
/// The alignment may absorb one insert base at equal cost.
#[test]
fn confident_panel_hit_trims_anywhere_in_the_end_zone() {
    let panel = barcode_panel(24, 16);
    let c = cfg_with(panel.clone(), 0.2, 150, true);
    let insert = splitmix_dna(9001, 2000);
    let confident = substituted(&panel[4].seq, &[3, 13]);
    let mut deep = insert[..80].to_vec();
    deep.extend_from_slice(&confident);
    deep.extend_from_slice(&insert[80..]);
    let segments = adapter_segments(&deep, &c);
    assert_eq!(segments.len(), 1);
    assert!((96..=97).contains(&segments[0].0), "{segments:?}");
    assert_eq!(segments[0].1, deep.len());
}

/// A 48-member panel of 24-base adapters shares one interior chance
/// bound: its members search the interior of a long read with fewer edits
/// than one of them alone, and a long adapter in the same set keeps its
/// budget.
#[test]
fn panel_members_share_one_interior_chance_bound() {
    let panel: Vec<Adapter> = (1..=48)
        .map(|i| ad(&format!("p{i}"), &splitmix_dna(900 + i, 24)))
        .collect();
    let alone = CandidateIndex::new(&panel[..1], 0.2, 150, true);
    let mut set = panel.clone();
    set.push(ad(
        "long",
        b"AATGTACTTCGTTCAGTTACGTATTGCTGGTTTTCGCATTTATCGTGAAACGCTTTC",
    ));
    let index = CandidateIndex::new(&set, 0.2, 150, true);
    for budget in &index.budgets[..48] {
        assert!(budget.interior(30_000) < alone.budgets[0].interior(30_000));
    }
    let own = Budget::new(&set[48].seq, 0.2, 150);
    assert_eq!(index.budgets[48].interior(30_000), own.interior(30_000));
}

/// An adapter whose seeds would open windows over most of a read is searched
/// over the whole read; an adapter with rare seeds keeps its windows.
#[test]
fn dense_seeds_select_the_whole_read_search() {
    let adapters = preset_ont();
    let index = CandidateIndex::new(&adapters, 0.2, 150, true);
    let position = |name: &str| adapters.iter().position(|a| a.name == name).unwrap();
    assert!(index.unfiltered[position("RAD")]);
    assert!(!index.unfiltered[position("LSK114_rear")]);
}

/// Constructs and insert-side flanks bound their barcode on either strand;
/// outer flanks and ligation adapters do not.
#[test]
fn barcode_bounding_entries_are_constructs_and_insert_side_flanks() {
    let adapters = preset_ont();
    let bounded = |name: &str| bounds_barcode(adapters.iter().find(|a| a.name == name).unwrap());
    for name in [
        "NB_construct",
        "PBC_rear",
        "RAD",
        "RLB_rear",
        "MAB_rear",
        "PCR1_front",
    ] {
        assert!(bounded(name), "{name}");
    }
    for name in [
        "NB_front",
        "RBK4_front",
        "LSK114_front",
        "PCR2_front",
        "BC01",
    ] {
        assert!(!bounded(name), "{name}");
    }
    let rear = entry(
        "rc",
        &reverse_complement(b"CCATATCCGTGTCGCCCTT"),
        Role::Barcode,
    );
    assert!(bounds_barcode(&rear));
}

/// Presence detection tallies a panel barcode at an end that a construct
/// trims, although the trimming pass skips the panel search there.
#[test]
fn presence_detection_tallies_the_barcode_behind_a_construct() {
    let adapters = super::preset::preset(&[super::preset::Kit::Nbd114]);
    let mut construct = b"ATTGCTAAGGTTAA".to_vec();
    construct.extend(reverse_complement(b"AAGAAAGTTGTCGGTGTCTTTGTG"));
    construct.extend_from_slice(b"CAGCACCT");
    let mut w = construct.clone();
    w.extend(splitmix_dna(9, 400));
    let c = cfg_with(adapters.clone(), 0.2, 150, true);
    let mut acted = vec![false; adapters.len()];
    let tallied = adapter_segments_tallied(&w, &c, &mut acted);
    assert_eq!(tallied, adapter_segments(&w, &c));
    let bc01 = adapters.iter().position(|a| a.name == "BC01").unwrap();
    assert!(acted[bc01]);
}

/// A terminal hit is found wherever it lies in the end zone, including across
/// the point where the singleton search cuts an end window in two.
#[test]
fn terminal_hit_is_found_across_the_window_split() {
    let adapter = b"GATCGGAAGAGCACACGTCTGAACTCCAGTC";
    let c = cfg_with(vec![ad("a", adapter)], 0.1, 150, false);
    for pos in (0..140).step_by(3) {
        let mut w = splitmix_dna(11, pos);
        w.extend_from_slice(adapter);
        w.extend(splitmix_dna(12, 600));
        let segs = adapter_segments(&w, &c);
        assert_eq!(
            segs,
            vec![(pos + adapter.len(), w.len())],
            "adapter at {pos}"
        );
    }
}

/// A barcode cut short by the read end is trimmed like a truncated adapter
/// when enough of it aligns flush with the end; a barcode construct, whose `N`
/// block aligns at no cost, is not matched partially.
#[test]
fn truncated_barcode_at_the_read_start_is_trimmed() {
    let barcode = b"TCCGCCATACCTCCATAGGT";
    let mut w = barcode[6..].to_vec();
    w.extend(splitmix_dna(13, 400));
    let c = cfg_with(vec![entry("bc", barcode, Role::Barcode)], 0.1, 150, true);
    assert_eq!(adapter_segments(&w, &c), vec![(barcode.len() - 6, w.len())]);

    let construct = entry(
        "construct",
        b"ATTGCTAAGGTTAANNNNNNNNNNNNNNNNNNNNNNNNCAGCACCT",
        Role::Barcode,
    );
    let index = CandidateIndex::new(&[construct], 0.1, 150, true);
    assert!(index.end_seeds.is_none());
}

/// A marker-gene primer splits freely where an amplicon-only preset makes it
/// an adapter; in a selection with a genomic kit it splits only beside the
/// rest of a junction stack.
#[test]
fn marker_primers_pair_outside_amplicon_selections() {
    use super::preset::{Kit, preset};
    let paired_of = |kits: &[Kit]| {
        let adapters = preset(kits);
        let index = CandidateIndex::new(&adapters, 0.2, 150, true);
        let i = adapters.iter().position(|a| a.name == "16S_27F").unwrap();
        assert!(index.split_classes[i] > 0);
        index.paired[i]
    };
    assert!(!paired_of(&[Kit::Mab114]));
    assert!(paired_of(&[Kit::Mab114, Kit::Lsk114]));
    let pcr = preset(&[Kit::Pcb114]);
    let index = CandidateIndex::new(&pcr, 0.2, 150, true);
    let i = pcr.iter().position(|a| a.name == "PCR2_front").unwrap();
    assert!(index.split_classes[i] > 0 && !index.paired[i]);
}

/// A 16S primer resolved from degenerate reads, as discovery assembles it,
/// counts as a marker primer and splits only in pairs in the primer role; an
/// unrelated primer splits alone.
#[test]
fn resolved_marker_primer_variants_split_in_pairs() {
    let resolved = entry("variant", b"AGAGTTTGATCCTGGCTCAG", Role::Primer);
    let pcr = entry("pcr", b"TTTCTGTTGGTGCTGATATTGC", Role::Primer);
    let index = CandidateIndex::new(&[resolved, pcr], 0.15, 150, true);
    assert_eq!(index.paired, [true, false]);
}

/// A marker primer inside a read is a genomic site and stays; beside the
/// other end's primer it is a junction and splits the read.
#[test]
fn marker_primer_splits_only_beside_a_junction_partner() {
    let forward = b"AGAGTTTGATCCTGGCTCAG";
    let reverse = b"TACGGTTACCTTGTTACGACTT";
    let c = cfg_with(
        vec![
            entry("27F", forward, Role::Primer),
            entry("1492R", reverse, Role::Primer),
        ],
        0.1,
        150,
        true,
    );
    let (left, right) = (splitmix_dna(21, 700), splitmix_dna(22, 700));
    let mut site = left.clone();
    site.extend_from_slice(forward);
    site.extend_from_slice(&right);
    assert_eq!(adapter_segments(&site, &c), vec![(0, site.len())]);

    let mut junction = left.clone();
    junction.extend(reverse_complement(reverse));
    junction.extend_from_slice(forward);
    junction.extend_from_slice(&right);
    let cut = left.len() + reverse.len() + forward.len();
    assert_eq!(
        adapter_segments(&junction, &c),
        vec![(0, left.len()), (cut, junction.len())]
    );
}

/// A resolved 16S primer, alone or with a few bases of its neighbour, is a
/// marker primer; an adapter assembled together with the primer behind it,
/// or an unrelated primer, is not.
#[test]
fn marker_primer_match_requires_the_sequence_to_be_the_primer() {
    assert!(matches_marker_primer(b"AGAGTTTGATCCTGGCTCAG", 0.15));
    assert!(matches_marker_primer(b"AGATAGAGTTTGATTCTGGCTCAG", 0.15));
    let mut with_adapter = b"ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT".to_vec();
    with_adapter.extend_from_slice(b"AGAGTTTGATCCTGGCTCAG");
    assert!(!matches_marker_primer(&with_adapter, 0.15));
    assert!(!matches_marker_primer(b"TTTCTGTTGGTGCTGATATTGC", 0.15));
    assert!(!matches_marker_primer(b"TGGTCCAGGATCAACA", 0.2));
}

/// A primer cut short at either end, as erosion or a layer boundary leaves
/// it, is a marker primer; a stretch from inside a primer is not.
#[test]
fn eroded_marker_primers_match() {
    assert!(matches_marker_primer(b"GAGTTTGATCATGGCTCAG", 0.2));
    assert!(matches_marker_primer(b"AAGTCGTAACAAGGTAAC", 0.2));
    assert!(matches_marker_primer(b"AGTCGTAACAAGGTAACCGTA", 0.2));
    assert!(!matches_marker_primer(b"GTTTGATCATGGCTC", 0.0));
}

/// A resolved marker primer takes the primer's ambiguity codes at the bases
/// they align to, on either strand; other sequences are unchanged.
#[test]
fn marker_codes_restore_degenerate_positions() {
    assert_eq!(
        with_marker_codes(b"AGAGTTTGATCCTGGCTCAG", 0.2),
        b"AGAGTTTGATYMTGGCTCAG"
    );
    assert_eq!(
        with_marker_codes(b"AAGTCGTAACAAGGTAAC", 0.2),
        b"AAGTCGTAACAAGGTARC"
    );
    assert_eq!(
        with_marker_codes(b"TTTCTGTTGGTGCTGATATTGC", 0.2),
        b"TTTCTGTTGGTGCTGATATTGC"
    );
}

/// An interior hit splits only when the whole adapter aligns: a site that
/// holds the inner part of an adapter and pays for the rest in edits, as a
/// genomic primer site does for an adapter assembled with its primer, stays.
#[test]
fn interior_split_needs_the_whole_adapter() {
    let primer = b"AAGTCGTAACAAGGTAGCCGTA";
    let flank = b"CGGTCCGAACGTG";
    let mut adapter = primer.to_vec();
    adapter.extend_from_slice(flank);
    let c = cfg_with(vec![entry("ad", &adapter, Role::Adapter)], 0.2, 150, true);
    let (left, right) = (splitmix_dna(31, 700), splitmix_dna(32, 700));
    let mut junction = left.clone();
    junction.extend_from_slice(&adapter);
    junction.extend_from_slice(&right);
    assert_eq!(
        adapter_segments(&junction, &c),
        vec![
            (0, left.len()),
            (left.len() + adapter.len(), junction.len())
        ]
    );
    let mut site = left.clone();
    site.extend_from_slice(primer);
    site.extend_from_slice(b"AGCTCAGTAGCTG");
    site.extend_from_slice(&right);
    assert_eq!(adapter_segments(&site, &c), vec![(0, site.len())]);
}

/// A PCS114 read end is trimmed through the UMI that follows the SSP, at
/// either end, and the `GGG` after it stays, also where the UMI pattern
/// aligns one repeat unit further at a higher cost; an end without a UMI is
/// trimmed at the end of the SSP.
#[test]
fn pcs114_ssp_is_trimmed_through_its_umi() {
    let ssp = b"TTTCTGTTGGTGCTGATATTGCTTT";
    let umi = b"ACGGTTCAGCTTGGAATTAGCCTTT";
    let insert = splitmix_dna(7, 400);
    let c = cfg_with(
        super::preset::preset(&[super::preset::Kit::Pcb114]),
        0.2,
        150,
        true,
    );
    let mut tagged = ssp.to_vec();
    tagged.extend_from_slice(umi);
    let mut forward = tagged.clone();
    forward.extend_from_slice(b"GGG");
    forward.extend_from_slice(&insert);
    assert_eq!(
        adapter_segments(&forward, &c),
        vec![(tagged.len(), forward.len())]
    );

    let mut reverse = insert.clone();
    reverse.extend_from_slice(b"CCC");
    reverse.extend(reverse_complement(&tagged));
    assert_eq!(adapter_segments(&reverse, &c), vec![(0, insert.len() + 3)]);

    let mut shifted = tagged.clone();
    shifted.extend_from_slice(b"GGGATTT");
    shifted.extend_from_slice(&insert);
    assert_eq!(
        adapter_segments(&shifted, &c),
        vec![(tagged.len(), shifted.len())]
    );

    let mut bare = ssp.to_vec();
    bare.extend_from_slice(b"GGG");
    bare.extend_from_slice(&insert);
    assert_eq!(adapter_segments(&bare, &c), vec![(ssp.len(), bare.len())]);
}
