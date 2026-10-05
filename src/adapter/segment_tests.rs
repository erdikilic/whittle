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
        amplicon: false,
        split_of: Vec::new(),
        split_opens: Vec::new(),
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
        amplicon: false,
        split_of: Vec::new(),
        split_opens: Vec::new(),
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
    let mut index = CandidateIndex::for_config(cfg);
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
            amplicon: false,
            split_of: Vec::new(),
            split_opens: Vec::new(),
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
/// adapter and primer at the read end, the lowest of its site for a marker
/// primer, and deep in the end zone for every entry of 24 bases or more.
#[test]
fn catalog_budgets_keep_their_per_pattern_tolerance() {
    let catalog = super::preset::preset(super::preset::Kit::ALL);
    let index = CandidateIndex::new(&catalog, 0.2, 150, true);
    let site = |adapter: &Adapter| {
        marker_primer_of(&adapter.seq, 0.2)
            .filter(|_| !is_umi(adapter))
            .map(|(primer, _)| super::catalog::marker_site(primer))
    };
    for (adapter, budget) in catalog.iter().zip(&index.budgets) {
        let own = Budget::new(&adapter.seq, 0.2, 150);
        let expected = match site(adapter) {
            Some(s) => catalog
                .iter()
                .filter(|other| site(other) == Some(s))
                .map(|other| Budget::new(&other.seq, 0.2, 150).k_end)
                .min()
                .unwrap(),
            None => own.k_end,
        };
        assert_eq!(budget.k_end, expected, "{}", adapter.name);
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
        let i = adapters.iter().position(|a| a.name == "16S_mix_F").unwrap();
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

/// Builds `left`, then each of `parts`, then `right`.
fn joined(left: &[u8], parts: &[&[u8]], right: &[u8]) -> Vec<u8> {
    let mut out = left.to_vec();
    for part in parts {
        out.extend_from_slice(part);
    }
    out.extend_from_slice(right);
    out
}

/// Two marker primers adjacent in the orientation of a junction, the first
/// reverse complemented as it closes one molecule and the second reading
/// into the next, split the read even where each carries more edits than the
/// interior budget allows a primer alone. Either primer alone stays.
#[test]
fn adjacent_marker_primer_pair_splits_beyond_the_interior_budget() {
    let forward = b"AGAGTTTGATCCTGGCTCAG";
    let reverse = b"TACGGTTACCTTGTTACGACTT";
    let c = cfg_with(
        vec![
            entry("27F", forward, Role::Primer),
            entry("1492R", reverse, Role::Primer),
        ],
        0.2,
        150,
        true,
    );
    let (left, right) = (splitmix_dna(41, 900), splitmix_dna(42, 900));
    let closing = substituted(&reverse_complement(reverse), &[5, 11, 16]);
    let opening = substituted(forward, &[4, 9, 14]);
    let junction = joined(&left, &[&closing, &opening], &right);
    let n = junction.len();
    let index = CandidateIndex::new(&c.adapters, c.error_rate, c.end_size, true);
    assert!(index.budgets.iter().all(|b| b.interior(n) < 3));
    let cut = left.len() + closing.len() + opening.len();
    assert_eq!(
        adapter_segments(&junction, &c),
        vec![(0, left.len()), (cut, n)]
    );
    for lone in [
        joined(&left, &[&closing], &right),
        joined(&left, &[&opening], &right),
    ] {
        assert_eq!(adapter_segments(&lone, &c), vec![(0, lone.len())]);
    }
}

/// Two marker primers facing each other, as the two sites of an amplicon lie
/// in a genome, do not form a junction, however close they are.
#[test]
fn marker_primers_facing_each_other_do_not_split() {
    let forward = b"AGAGTTTGATCCTGGCTCAG";
    let reverse = b"TACGGTTACCTTGTTACGACTT";
    let c = cfg_with(
        vec![
            entry("27F", forward, Role::Primer),
            entry("1492R", reverse, Role::Primer),
        ],
        0.2,
        150,
        true,
    );
    let (left, right) = (splitmix_dna(43, 900), splitmix_dna(44, 900));
    let opening = substituted(forward, &[4, 9, 14]);
    let closing = substituted(&reverse_complement(reverse), &[5, 11, 16]);
    let inward = joined(&left, &[&opening, b"ACGTA", &closing], &right);
    assert_eq!(adapter_segments(&inward, &c), vec![(0, inward.len())]);
}

/// One marker primer closing a molecule and again opening the next, as where
/// a read continues into the reverse strand of another copy, is a junction.
#[test]
fn inverted_copies_of_one_marker_primer_split() {
    let forward = b"AGAGTTTGATCCTGGCTCAG";
    let c = cfg_with(vec![entry("27F", forward, Role::Primer)], 0.2, 150, true);
    let (left, right) = (splitmix_dna(45, 900), splitmix_dna(46, 900));
    let closing = substituted(&reverse_complement(forward), &[4, 10, 15]);
    let opening = substituted(forward, &[3, 8, 13]);
    let junction = joined(&left, &[&closing, &opening], &right);
    let cut = left.len() + closing.len() + opening.len();
    assert_eq!(
        adapter_segments(&junction, &c),
        vec![(0, left.len()), (cut, junction.len())]
    );
}

/// Marker primers cut short on their outer side, as discovery assembles them
/// from eroded read ends, pair across the bases they lack at a junction that
/// holds both primers whole.
#[test]
fn eroded_marker_primer_forms_pair_across_the_bases_they_lack() {
    let forward = b"AGAGTTTGATCCTGGCTCAG";
    let closing_whole = reverse_complement(b"TACGGTTACCTTGTTACGACTT");
    let c = cfg_with(
        vec![
            entry("five", &forward[1..], Role::Primer),
            entry("three", &closing_whole[..17], Role::Primer),
        ],
        0.2,
        150,
        true,
    );
    let (left, right) = (splitmix_dna(47, 900), splitmix_dna(48, 900));
    let closing = substituted(&closing_whole, &[5, 11]);
    let opening = substituted(forward, &[5, 10, 15]);
    let junction = joined(&left, &[&closing, &opening], &right);
    let n = junction.len();
    let index = CandidateIndex::new(&c.adapters, c.error_rate, c.end_size, true);
    assert!(index.budgets[0].interior(n) < 3 && index.budgets[1].interior(n) < 2);
    let cut = left.len() + closing.len() + opening.len();
    assert_eq!(
        adapter_segments(&junction, &c),
        vec![(0, left.len()), (cut, n)]
    );
}

/// Builds the marker primers of a 16S library in the primer role, as a set
/// of `entries`, at error rate 0.2 and end zone 150, for an amplicon library
/// when `amplicon`.
fn marker_set(entries: Vec<Adapter>, amplicon: bool) -> AdapterConfig {
    let mut c = cfg_with(entries, 0.2, 150, true);
    c.amplicon = amplicon;
    c
}

/// In an amplicon library one marker primer inside a read is a chimera
/// junction and splits the read, on either strand; outside one it is a
/// genomic site and the read stays whole.
#[test]
fn amplicon_library_splits_at_a_single_marker_primer() {
    let forward = b"AGAGTTTGATCCTGGCTCAG";
    let reverse = b"TACGGTTACCTTGTTACGACTT";
    let entries = || {
        vec![
            entry("27F", forward, Role::Primer),
            entry("1492R", reverse, Role::Primer),
        ]
    };
    let (left, right) = (splitmix_dna(51, 900), splitmix_dna(52, 900));
    for primer in [
        substituted(forward, &[9]),
        reverse_complement(reverse),
        reverse.to_vec(),
    ] {
        let read = joined(&left, &[&primer], &right);
        let cut = left.len() + primer.len();
        assert_eq!(
            adapter_segments(&read, &marker_set(entries(), true)),
            vec![(0, left.len()), (cut, read.len())]
        );
        assert_eq!(
            adapter_segments(&read, &marker_set(entries(), false)),
            vec![(0, read.len())]
        );
    }
}

/// A marker primer form cut short on its outer side, as discovery assembles
/// it, splits an amplicon library where the whole primer lies inside the
/// read, and the excision covers the whole primer; where only the bases of
/// the form are present the primer is not whole and the read stays.
#[test]
fn amplicon_splits_need_the_whole_marker_primer() {
    let forward = b"AGAGTTTGATCCTGGCTCAG";
    let c = marker_set(vec![entry("five", &forward[4..], Role::Primer)], true);
    let (left, right) = (splitmix_dna(53, 900), splitmix_dna(54, 900));
    let whole = joined(&left, &[forward], &right);
    let cut = left.len() + forward.len();
    assert_eq!(
        adapter_segments(&whole, &c),
        vec![(0, left.len()), (cut, whole.len())]
    );
    let mut lead = left.clone();
    let tail = lead.len() - 4;
    for (base, primer_base) in lead[tail..].iter_mut().zip(&forward[..4]) {
        *base = complement(*primer_base);
    }
    let partial = joined(&lead, &[&forward[4..]], &right);
    assert_eq!(adapter_segments(&partial, &c), vec![(0, partial.len())]);
}

/// The pair budget of a marker primer never exceeds its terminal budget,
/// does not grow with the read-length class, and keeps chance pairs of the
/// set within the interior bound wherever an edit is admitted.
#[test]
fn pair_budgets_respect_the_interior_chance_bound() {
    use super::preset::{Kit, preset};
    let adapters = preset(&[Kit::Mab114, Kit::Lsk114]);
    let index = CandidateIndex::new(&adapters, 0.2, 150, true);
    let paired: Vec<usize> = (0..adapters.len()).filter(|&i| index.paired[i]).collect();
    assert_eq!(paired.len(), 12);
    for &i in &paired {
        let budget = &index.budgets[i];
        assert!(budget.k_pair[0] >= budget.interior_max());
        assert!(budget.k_pair[0] <= budget.k_end);
        assert!(budget.k_pair.windows(2).all(|pair| pair[0] >= pair[1]));
    }
    // The primers of one marker-primer site count as their member of
    // highest chance.
    let site = |i: usize| {
        let (primer, _) = marker_primer_of(&adapters[i].seq, 0.2).unwrap();
        super::catalog::marker_site(primer)
    };
    for class in 0..INTERIOR_CLASSES {
        let mut by_site: BTreeMap<usize, f64> = BTreeMap::new();
        for &i in &paired {
            let k = index.budgets[i].k_pair[class];
            let rate = by_site.entry(site(i)).or_default();
            *rate = rate.max(chance_cumulative(&adapters[i].seq, k)[k]);
        }
        let rate: f64 = by_site.values().sum();
        let edited = paired.iter().any(|&i| index.budgets[i].k_pair[class] > 0);
        let expected = rate * rate * PAIR_OFFSETS as f64 * interior_positions(class);
        assert!(!edited || expected <= INTERIOR_CHANCE_HITS_PER_READ);
    }
    let unpaired = (0..adapters.len()).find(|&i| !index.paired[i]).unwrap();
    assert_eq!(index.budgets[unpaired].k_pair, [0; INTERIOR_CLASSES]);
}

/// Every marker primer is listed under exactly one site, and every site
/// lists only marker primers.
#[test]
fn marker_sites_partition_the_marker_primers() {
    use super::catalog::{MARKER_PRIMERS, MARKER_SITES};
    for primer in MARKER_PRIMERS {
        let listed = MARKER_SITES
            .iter()
            .filter(|site| site.contains(primer))
            .count();
        assert_eq!(listed, 1, "{}", String::from_utf8_lossy(primer));
    }
    let sited: usize = MARKER_SITES.iter().map(|site| site.len()).sum();
    assert_eq!(sited, MARKER_PRIMERS.len());
}

/// Under the `mab114` preset alone, its twelve primers take the adapter
/// role (`preset::preset`'s amplicon-only conversion), so none of them
/// joins another entry's chance family or pairs; only the set-wide
/// terminal and interior chance bounds can still move another entry's
/// budget when the twelve kit primers stand in for the four legacy
/// primers 16S 27F and 1492R, ITS1F and ITS4. `LSK109_front`, `RAD` and
/// `MAB_rear` move by one edit in one or two read-length classes; every
/// other entry of the preset, including `TP01` to `TP24` and the barcode
/// flank `RBK4_front`, keeps its budgets exactly.
#[test]
fn kit_primer_variants_move_named_budgets_under_mab114_alone() {
    use super::catalog::CATALOG;
    use super::preset::{Kit, preset};
    let one_per_site: [(&str, &[u8]); 4] = [
        ("16S_27F", b"AGAGTTTGATYMTGGCTCAG"),
        ("16S_1492R", b"TACGGYTACCTTGTTACGACTT"),
        ("ITS1F", b"CTTGGTCATTTAGAGGAAGTAA"),
        ("ITS4", b"TCCTCCGCTTATTGATATGC"),
    ];
    let is_kit_primer = |name: &str| {
        CATALOG.iter().any(|&(entry, role, kits, _)| {
            entry == name && kits == [Kit::Mab114] && role == Role::Primer
        })
    };
    let kit = preset(&[Kit::Mab114]);
    let mut single: Vec<Adapter> = kit
        .iter()
        .filter(|a| !is_kit_primer(&a.name))
        .cloned()
        .collect();
    single.extend(one_per_site.iter().map(|&(name, seq)| Adapter {
        name: name.to_string(),
        seq: seq.to_vec(),
        // preset() gives the mab114 primers the adapter role since the kit
        // is amplicon-only; the stand-ins take the same role.
        role: Role::Adapter,
    }));
    let variants = CandidateIndex::new(&kit, 0.2, 150, true);
    let baseline = CandidateIndex::new(&single, 0.2, 150, true);

    #[rustfmt::skip]
    let expected_move: &[(&str, [usize; INTERIOR_CLASSES], [usize; INTERIOR_CLASSES])] = &[
        (
            "LSK109_front",
            [4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1],
            [4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1],
        ),
        (
            "RAD",
            [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8],
            [10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8],
        ),
        (
            "MAB_rear",
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ),
    ];

    let mut compared = 0;
    let mut moved = 0;
    for (i, adapter) in kit.iter().enumerate() {
        if is_kit_primer(&adapter.name) {
            continue;
        }
        let j = single.iter().position(|a| a.name == adapter.name).unwrap();
        let (a, b) = (&variants.budgets[i], &baseline.budgets[j]);
        assert_eq!(a.k_end, b.k_end, "{}", adapter.name);
        assert_eq!(a.k_far, b.k_far, "{}", adapter.name);
        assert_eq!(a.k_pair, b.k_pair, "{}", adapter.name);
        match expected_move
            .iter()
            .find(|(name, _, _)| *name == adapter.name)
        {
            Some((_, want_variants, want_baseline)) => {
                assert_eq!(&a.k_mid, want_variants, "{}", adapter.name);
                assert_eq!(&b.k_mid, want_baseline, "{}", adapter.name);
                moved += 1;
            },
            None => assert_eq!(a.k_mid, b.k_mid, "{}", adapter.name),
        }
        compared += 1;
    }
    assert_eq!(moved, expected_move.len());
    assert_eq!(compared, 31, "{compared}");
}

/// A barcode entry whose sequence matches a marker primer within its edit
/// budget keeps its own chance family: only a primer-role entry or a split
/// primer joins a marker-primer site.
#[test]
fn a_barcode_within_a_marker_primers_budget_keeps_its_own_family() {
    let primer = Adapter {
        name: "16S_27F".to_string(),
        role: Role::Primer,
        seq: b"AGAGTTTGATCCTGGCTCAG".to_vec(),
    };
    // An eroded variant of the same primer, within its edit budget, tagged
    // as a barcode instead of a primer.
    let barcode = Adapter {
        name: "fake_barcode".to_string(),
        role: Role::Barcode,
        seq: b"GAGTTTGATCATGGCTCAG".to_vec(),
    };
    assert!(matches_marker_primer(&barcode.seq, 0.2));

    // The same variant, but tagged as another primer, does join the site
    // and is bounded together with it.
    let mut as_primer = barcode.clone();
    as_primer.role = Role::Primer;
    let grouped = CandidateIndex::new(&[primer.clone(), as_primer], 0.2, 150, true);
    let with_barcode = CandidateIndex::new(&[primer.clone(), barcode.clone()], 0.2, 150, true);
    assert_ne!(
        (grouped.budgets[0].k_end, grouped.budgets[0].k_mid),
        (with_barcode.budgets[0].k_end, with_barcode.budgets[0].k_mid),
        "a same-role variant and a barcode should not bound the primer the same way"
    );

    // The barcode keeps the budget a plain entry of its own sequence gets,
    // unaffected by sitting beside a marker primer it happens to match.
    let barcode_alone = CandidateIndex::new(std::slice::from_ref(&barcode), 0.2, 150, true);
    let (a, b) = (&with_barcode.budgets[1], &barcode_alone.budgets[0]);
    assert_eq!(
        (a.k_end, a.k_far, a.k_mid, a.k_pair),
        (b.k_end, b.k_far, b.k_mid, b.k_pair)
    );
}

/// The kit primers of one marker-primer site form one chance family, so the
/// variants of a site leave the budgets of every other entry where the
/// catalog of one primer per site put them: under the `ont` and `all`
/// presets, each entry outside the amplicon kit keeps its terminal and
/// interior budgets when the twelve kit primers stand in for the four
/// primers 16S 27F and 1492R, ITS1F and ITS4. The pair budget of the UMI,
/// the one such entry that pairs with the primers, follows theirs.
#[test]
fn kit_primer_variants_leave_the_budgets_of_other_entries() {
    use super::catalog::CATALOG;
    use super::preset::{Kit, preset};
    let one_per_site: [(&str, &[u8]); 4] = [
        ("16S_27F", b"AGAGTTTGATYMTGGCTCAG"),
        ("16S_1492R", b"TACGGYTACCTTGTTACGACTT"),
        ("ITS1F", b"CTTGGTCATTTAGAGGAAGTAA"),
        ("ITS4", b"TCCTCCGCTTATTGATATGC"),
    ];
    let amplicon_only = |name: &str| {
        CATALOG
            .iter()
            .any(|&(entry, _, kits, _)| entry == name && kits == [Kit::Mab114])
    };
    for kits in [Kit::ONT, Kit::ALL] {
        let kit = preset(kits);
        let mut single: Vec<Adapter> = kit
            .iter()
            .filter(|a| !amplicon_only(&a.name) || a.role != Role::Primer)
            .cloned()
            .collect();
        single.extend(one_per_site.iter().map(|&(name, seq)| Adapter {
            name: name.to_string(),
            seq: seq.to_vec(),
            role: Role::Primer,
        }));
        let variants = CandidateIndex::new(&kit, 0.2, 150, true);
        let sites = CandidateIndex::new(&single, 0.2, 150, true);
        let mut compared = 0;
        for (i, adapter) in kit.iter().enumerate() {
            if amplicon_only(&adapter.name) {
                continue;
            }
            let j = single.iter().position(|a| a.name == adapter.name).unwrap();
            let (a, b) = (&variants.budgets[i], &sites.budgets[j]);
            assert_eq!(
                (a.k_end, a.k_far, a.k_mid),
                (b.k_end, b.k_far, b.k_mid),
                "{}",
                adapter.name
            );
            if !is_umi(adapter) {
                assert_eq!(a.k_pair, b.k_pair, "{}", adapter.name);
            }
            compared += 1;
        }
        assert!(compared > 100, "{compared}");
    }
}

/// A marker primer's direction is read from the catalog primer it matches:
/// as synthesized it reads into the insert, and its reverse complement, as a
/// read's 3' end holds it, reads out of it.
#[test]
fn marker_primer_direction_follows_the_catalog_primer() {
    assert_eq!(
        marker_primer_opens(b"AGAGTTTGATCCTGGCTCAG", 0.2),
        Some(true)
    );
    assert_eq!(marker_primer_opens(b"GAGTTTGATCATGGCTCAG", 0.2), Some(true));
    assert_eq!(marker_primer_opens(b"AAGTCGTAACAAGGTAAC", 0.2), Some(false));
    assert_eq!(
        marker_primer_opens(b"TACGGTTACCTTGTTACGACTT", 0.2),
        Some(true)
    );
    assert_eq!(marker_primer_opens(b"TTTCTGTTGGTGCTGATATTGC", 0.2), None);
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
        b"AGRGTTYGATYMTGGCTCAG"
    );
    assert_eq!(
        with_marker_codes(b"AGAGTTTGATCCTGGCTTAG", 0.2),
        b"AGAGTTTGATCCTGGCTTAG"
    );
    assert_eq!(
        with_marker_codes(b"TACGGTTACCTTGTTACGACTT", 0.2),
        b"TACGGYTACCTTGTTACGACTT"
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

/// Every primer of the amplicon kit is a marker primer that reads into the
/// insert as listed, and pairs in the primer role.
#[test]
fn kit_primers_are_marker_primers_as_listed() {
    use super::preset::{Kit, preset};
    let adapters = preset(&[Kit::Mab114, Kit::Lsk114]);
    let index = CandidateIndex::new(&adapters, 0.2, 150, true);
    let mut primers = 0;
    for (i, adapter) in adapters.iter().enumerate() {
        if !(adapter.name.starts_with("16S_") || adapter.name.starts_with("ITS")) {
            continue;
        }
        primers += 1;
        assert_eq!(
            marker_primer_of(&adapter.seq, 0.2).map(|(_, opens)| opens),
            Some(true),
            "{}",
            adapter.name
        );
        assert!(
            index.paired[i] && index.opens[i] == Opens::AsGiven,
            "{}",
            adapter.name
        );
    }
    assert_eq!(primers, 12);
}

/// The community primers outside the amplicon kit are marker primers as a
/// FASTA gives them. The classic 16S 27F is an instance of the kit's forward
/// mix; the 22-base 1492R, which the kit's reverse mix also aligns within,
/// and fungal ITS1F match their own entries.
#[test]
fn community_marker_primers_are_recognised() {
    let classic: [(&str, &[u8], &[u8]); 3] = [
        ("27F", b"AGAGTTTGATYMTGGCTCAG", b"AGRGTTYGATYMTGGCTCAG"),
        (
            "1492R",
            b"TACGGYTACCTTGTTACGACTT",
            b"TACGGYTACCTTGTTACGACTT",
        ),
        (
            "ITS1F",
            b"CTTGGTCATTTAGAGGAAGTAA",
            b"CTTGGTCATTTAGAGGAAGTAA",
        ),
    ];
    for (name, seq, primer) in classic {
        assert_eq!(marker_primer_of(seq, 0.2), Some((primer, true)), "{name}");
        assert_eq!(
            marker_primer_of(&reverse_complement(seq), 0.2),
            Some((primer, false)),
            "{name}"
        );
    }
    let entries: Vec<Adapter> = classic
        .iter()
        .map(|(name, seq, _)| entry(name, seq, Role::Primer))
        .collect();
    let index = CandidateIndex::new(&entries, 0.2, 150, true);
    for (i, (name, _, _)) in classic.iter().enumerate() {
        assert!(
            index.paired[i] && index.opens[i] == Opens::AsGiven,
            "{name}"
        );
    }
    let mut c = cfg_with(entries, 0.2, 150, true);
    c.set_amplicon(true);
    let index = CandidateIndex::for_config(&c);
    for (i, (name, _, primer)) in classic.iter().enumerate() {
        assert_eq!(
            index.whole_primers[i].as_ref().map(|w| w.seq.as_slice()),
            Some(*primer),
            "{name}"
        );
    }
}

/// An MAB114 amplicon that carries variant primers of the kit's mixes is
/// trimmed to its insert under the amplicon preset: 16S_Bor_F with
/// 16S_Chl_R, and ITS1_Mal with ITS4_Pyt, behind the barcode and its flanks.
#[test]
fn mab114_variant_primer_pairs_trim_to_the_insert() {
    use super::preset::{Kit, preset};
    let c = cfg_with(preset(&[Kit::Mab114]), 0.2, 150, true);
    let head = |primer: &[u8]| {
        [
            b"GCTTGGGTGTTTAACC".as_slice(),
            b"GCACCTGGAACTTGTGCCTTCCAC",
            b"CCATATCCGTGTCGCCCTT",
            primer,
        ]
        .concat()
    };
    let pairs: [(&[u8], &[u8], usize); 2] = [
        (b"AGAGTTTGATCCTGGCTTAG", b"GGGCTACCTTGTTACGACTT", 1500),
        (b"TCTGTAGGTGAACCTGCAG", b"TCCTCCGCTTATTAATATGC", 600),
    ];
    for (i, (forward, reverse, len)) in pairs.into_iter().enumerate() {
        let (left, right) = (head(forward), reverse_complement(&head(reverse)));
        let insert = splitmix_dna(700 + i as u64, len);
        let read = [left.as_slice(), &insert, &right].concat();
        let span = vec![(left.len(), left.len() + len)];
        assert_eq!(adapter_segments(&read, &c), span);
        let flipped = reverse_complement(&read);
        let span = vec![(right.len(), right.len() + len)];
        assert_eq!(adapter_segments(&flipped, &c), span);
    }
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

/// The sheet primers of two targets, `A` and `B`, as SplitMix sequences:
/// forward and reverse primer of each, as ordered.
fn split_primers() -> Vec<crate::split::Primer> {
    [
        ("fA", 111, 20),
        ("rA", 112, 22),
        ("fB", 113, 21),
        ("rB", 114, 20),
    ]
    .into_iter()
    .map(|(name, seed, len)| crate::split::Primer {
        name: name.into(),
        seq: splitmix_dna(seed, len),
    })
    .collect()
}

/// A split configuration: one sequencing adapter at error rate 0.2 and end
/// zone 150, with the two-target sheet of `split_primers` attached.
fn split_cfg() -> AdapterConfig {
    let mut c = cfg_with(
        vec![ad("lsk", b"AATGTACTTCGTTCAGTTACGTATTGCT")],
        0.2,
        150,
        true,
    );
    c.attach_split(&split_primers());
    c
}

/// Searches both strands of every paired entry over the full read and
/// applies every undominated pair in adapter and strand order.
fn full_read_pairs(ctx: Context<'_>, engine: &mut Engine<'_>, keep: &mut Keep<'_>) {
    let n = ctx.read.window.len();
    let mut hits = Vec::new();
    for (idx, adapter) in ctx.cfg.adapters.iter().enumerate() {
        if ctx.index.paired[idx] && n >= adapter.seq.len() {
            search(
                engine,
                ctx.index,
                idx,
                &adapter.seq,
                ctx.read.strands(0, n),
                ctx.index.budgets[idx].pair(n),
                |hit| hits.push((idx, hit)),
            );
        }
    }
    let best: Vec<_> = hits
        .iter()
        .copied()
        .filter(|&(idx, h)| {
            !hits.iter().any(|&(other, g)| {
                idx == other
                    && h.rc == g.rc
                    && g.cost < h.cost
                    && g.start < h.end
                    && h.start < g.end
            })
        })
        .collect();
    for &(ci, c) in &best {
        for &(oi, o) in &best {
            if ctx.index.opens[ci].reads_out(c.rc)
                && ctx.index.opens[oi].reads_in(o.rc)
                && o.start > c.start
                && o.end > c.end
                && o.start.abs_diff(c.end) <= FLANK_SLACK
            {
                keep.accept_pair((ci, c), (oi, o));
            }
        }
    }
}

#[test]
fn junction_search_matches_full_read_search() {
    use super::preset::{Kit, preset};
    let mut kit = cfg_with(preset(&[Kit::Mab114]), 0.2, 150, true);
    kit.attach_split(&crate::split::Sheet::preset("mab114").unwrap().primers);
    let mut primers = split_primers();
    primers.push(Primer {
        name: "rB_rc".into(),
        seq: reverse_complement(&primers[3].seq),
    });
    let mut reversed = cfg_with(
        vec![ad("fA_rc", &reverse_complement(&primers[0].seq))],
        0.2,
        150,
        true,
    );
    reversed.attach_split(&primers);
    let index = CandidateIndex::for_config(&reversed);
    assert!(index.opens.contains(&Opens::Reversed));
    assert!(index.opens.contains(&Opens::Both));
    let mut rng = Lcg(0x6a75_6e63_7469_6f6e);
    for c in [split_cfg(), kit, reversed] {
        let index = CandidateIndex::for_config(&c);
        let paired: Vec<_> = (0..c.adapters.len()).filter(|&i| index.paired[i]).collect();
        let mut reads = planted_primer_reads(100);
        for case in 0..400 {
            let head = match case % 4 {
                0 => FLANK_SLACK,
                1 => FLANK_SLACK + 1,
                _ => rng.below(2000),
            };
            let mut read = rng.dna(head);
            for closing in [true, false] {
                let idx = paired[rng.below(paired.len())];
                let mut seq: Vec<_> = c.adapters[idx]
                    .seq
                    .iter()
                    .map(|&b| {
                        let bases = iupac_bases(b).unwrap();
                        bases[rng.below(bases.len())]
                    })
                    .collect();
                for _ in 0..rng.below(4) {
                    let at = rng.below(seq.len());
                    match rng.below(3) {
                        0 => seq[at] = b"ACGT"[rng.below(4)],
                        1 => {
                            seq.remove(at);
                        },
                        _ => seq.insert(at, b"ACGT"[rng.below(4)]),
                    }
                }
                if index.opens[idx].reads_out(false) != closing {
                    seq = reverse_complement(&seq);
                }
                read.extend(seq);
                if closing {
                    let gap = rng.below(2 * FLANK_SLACK + 3);
                    if gap < FLANK_SLACK {
                        read.truncate(read.len() - gap);
                    } else {
                        read.extend(rng.dna(gap - FLANK_SLACK));
                    }
                }
            }
            let tail = match case % 5 {
                0 => 0,
                1 => FLANK_SLACK,
                2 => FLANK_SLACK + 1,
                _ => rng.below(2000),
            };
            read.extend(rng.dna(tail));
            if case % 3 == 0 {
                let at = rng.below(read.len());
                read[at] = b'N';
            }
            reads.push(read.clone());
            reads.push(reverse_complement(&read));
        }
        let mut excised = 0;
        for (i, read) in reads.iter().enumerate() {
            let run = |reference| {
                with_engine(read, &c, |ctx, engine| {
                    let mut keep = Keep::new(&c, ctx.index, read.len(), true);
                    search_terminal(ctx, (0, read.len()), engine, &mut keep);
                    search_interior(ctx, engine, &mut keep);
                    if reference {
                        full_read_pairs(ctx, engine, &mut keep);
                    } else {
                        search_pairs(ctx, engine, &mut keep);
                    }
                    keep.place_split();
                    keep.settle();
                    (keep.acted.clone(), keep.into_cuts(1))
                })
            };
            let expected = run(true);
            excised += usize::from(!expected.1.2.is_empty());
            assert_eq!(run(false), expected, "read {i}");
        }
        assert!(excised > 50, "{excised} reads excised");
    }
}

/// The anchored budgets rescoring takes for a sheet (`anchored_budgets`)
/// are those the candidate index gives the sheet's split entries, whatever
/// else the adapter set holds: the MAB114 preset with the `mab114` sheet,
/// and an adapter with a panel of twenty random pairs of 15-base primers,
/// whose set-wide bound lowers every anchored budget below the error-rate
/// ceiling.
#[test]
fn sheet_anchored_budgets_match_the_index() {
    let mab114 = crate::split::Sheet::preset("mab114").unwrap().primers;
    let panel: Vec<crate::split::Primer> = (0..40)
        .map(|i| crate::split::Primer {
            name: format!("p{i}"),
            seq: splitmix_dna(7_000 + i, 15),
        })
        .collect();
    let sets = [
        (
            super::preset::preset(&[super::preset::Kit::Mab114]),
            mab114,
            false,
        ),
        (
            vec![ad("lsk", b"AATGTACTTCGTTCAGTTACGTATTGCT")],
            panel,
            true,
        ),
    ];
    for (adapters, primers, below_ceiling) in sets {
        let mut c = cfg_with(adapters, 0.2, 150, true);
        c.attach_split(&primers);
        let index = CandidateIndex::for_config(&c);
        let seqs: Vec<&[u8]> = primers.iter().map(|p| p.seq.as_slice()).collect();
        let anchored = anchored_budgets(&seqs, 0.2);
        for (adapter_idx, primer) in c.split_of.iter().enumerate() {
            if let Some(primer) = *primer {
                assert_eq!(
                    index.budgets[adapter_idx].k_anchor, anchored[primer],
                    "{}",
                    c.adapters[adapter_idx].name
                );
            }
        }
        if below_ceiling {
            assert!(anchored.iter().all(|&k| k < edit_budget(0.2, 15)));
        }
    }
}

/// Returns a segment at `[start, end)` with its primer loci.
fn seg(start: usize, end: usize, five: Option<Locus>, three: Option<Locus>) -> Segment {
    Segment {
        start,
        end,
        five,
        three,
    }
}

/// Returns the locus `[start, end)`, its outer edge no read end, located by
/// an exact whole hit of the entry at `entry`, as its reverse complement
/// when `rc`.
fn locus(start: usize, end: usize, entry: usize, rc: bool) -> Option<Locus> {
    sited(start, end, PrimerSite { entry, rc, cost: 0 })
}

/// Returns the locus `[start, end)`, its outer edge no read end, located by
/// the whole hit `site`.
fn sited(start: usize, end: usize, site: PrimerSite) -> Option<Locus> {
    Some(Locus {
        start,
        end,
        outer_open: false,
        boundary: false,
        site: Some(site),
    })
}

/// Returns the locus `[start, end)` located by an excision, which carries
/// no site.
fn excised(start: usize, end: usize) -> Option<Locus> {
    Some(Locus {
        start,
        end,
        outer_open: false,
        boundary: false,
        site: None,
    })
}

/// Without split primers the segments carry no loci, and a configuration
/// whose `split_of` marks no entry gives the spans of one without it.
#[test]
fn split_off_output_unchanged() {
    let adapter = b"GGGGTTTTGGGGTTTTGGGG";
    let primer = b"CACACAGAGAGACACACAGAGA";
    let mut w = vec![b'A'; 80];
    w.extend_from_slice(&reverse_complement(primer));
    w.extend_from_slice(adapter);
    w.extend_from_slice(primer);
    w.extend_from_slice(&[b'C'; 80]);
    let c = cfg_with(
        vec![ad("a", adapter), entry("p", primer, Role::Primer)],
        0.2,
        30,
        true,
    );
    assert!(c.split_of.is_empty());
    let expected = vec![(0, 80), (144, 224)];
    assert_eq!(adapter_segments(&w, &c), expected);
    assert_eq!(
        adapter_segments_annotated(&w, &c),
        vec![seg(0, 80, None, None), seg(144, 224, None, None)]
    );
    let mut unmarked = c.clone();
    unmarked.split_of = vec![None; c.adapters.len()];
    unmarked.split_opens = vec![Opens::AsGiven; c.adapters.len()];
    assert_eq!(adapter_segments(&w, &unmarked), expected);
}

/// A plus-strand amplicon reports the forward primer at its 5' end and the
/// reverse primer, reverse complemented, at its 3' end.
#[test]
fn loci_at_both_ends() {
    let p = split_primers();
    let read = joined(
        &p[0].seq,
        &[&splitmix_dna(121, 800)],
        &reverse_complement(&p[1].seq),
    );
    let n = read.len();
    let (f, r) = (p[0].seq.len(), p[1].seq.len());
    assert_eq!(
        adapter_segments_annotated(&read, &split_cfg()),
        vec![seg(
            f,
            n - r,
            locus(0, f, 1, false),
            locus(n - r, n, 2, true)
        )]
    );
}

/// A minus-strand amplicon reports the reverse primer at its 5' end and the
/// forward primer, reverse complemented, at its 3' end.
#[test]
fn loci_on_minus_strand_read() {
    let p = split_primers();
    let read = joined(
        &p[1].seq,
        &[&splitmix_dna(122, 800)],
        &reverse_complement(&p[0].seq),
    );
    let n = read.len();
    let (f, r) = (p[0].seq.len(), p[1].seq.len());
    assert_eq!(
        adapter_segments_annotated(&read, &split_cfg()),
        vec![seg(
            r,
            n - f,
            locus(0, r, 2, false),
            locus(n - f, n, 1, true)
        )]
    );
}

/// A chimera of two amplicons splits at the junction, and each piece
/// reports the primer on its side of the junction.
#[test]
fn junction_pieces_get_loci() {
    let p = split_primers();
    let (ins1, ins2) = (splitmix_dna(123, 800), splitmix_dna(124, 800));
    let (fa, ra_rc) = (p[0].seq.clone(), reverse_complement(&p[1].seq));
    let (fb, rb_rc) = (p[2].seq.clone(), reverse_complement(&p[3].seq));
    let read = joined(&fa, &[&ins1, &ra_rc, &fb, &ins2], &rb_rc);
    let n = read.len();
    let junction = fa.len() + ins1.len();
    let second = junction + ra_rc.len() + fb.len();
    assert_eq!(
        adapter_segments_annotated(&read, &split_cfg()),
        vec![
            seg(
                fa.len(),
                junction,
                locus(0, fa.len(), 1, false),
                excised(junction, junction + ra_rc.len())
            ),
            seg(
                second,
                n - rb_rc.len(),
                excised(second - fb.len(), second),
                locus(n - rb_rc.len(), n, 4, true)
            ),
        ]
    );
}

/// Two amplicons joined through an adapter split at the adapter, and the
/// sheet primers beside it, backed by its excision, are the loci of the
/// pieces on either side.
#[test]
fn adapter_junction_splits_beside_split_primers() {
    let p = split_primers();
    let lsk = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let (ins1, ins2) = (splitmix_dna(127, 800), splitmix_dna(128, 800));
    let (fa, ra_rc) = (p[0].seq.clone(), reverse_complement(&p[1].seq));
    let (fb, rb_rc) = (p[2].seq.clone(), reverse_complement(&p[3].seq));
    let read = joined(&fa, &[&ins1, &ra_rc, lsk, &fb, &ins2], &rb_rc);
    let n = read.len();
    let junction = fa.len() + ins1.len();
    let second = junction + ra_rc.len() + lsk.len() + fb.len();
    assert_eq!(
        adapter_segments_annotated(&read, &split_cfg()),
        vec![
            seg(
                fa.len(),
                junction,
                locus(0, fa.len(), 1, false),
                excised(junction, junction + ra_rc.len())
            ),
            seg(
                second,
                n - rb_rc.len(),
                excised(second - fb.len(), second),
                locus(n - rb_rc.len(), n, 4, true)
            ),
        ]
    );
}

/// A lone sheet primer inside a read and a catalog marker primer beside it,
/// within the end zone of each other, back neither excision: the read stays
/// whole under `--split-by`.
#[test]
fn nested_sheet_primer_does_not_back_a_marker_primer() {
    let p = split_primers();
    let marker = b"AGAGTTTGATCATGGCTCAG";
    let entries = vec![
        ad("lsk", b"AATGTACTTCGTTCAGTTACGTATTGCT"),
        entry("27F", b"AGAGTTTGATYMTGGCTCAG", Role::Primer),
    ];
    let mut c = cfg_with(entries, 0.2, 150, true);
    c.attach_split(&p);
    assert!(CandidateIndex::for_config(&c).paired[1]);
    let (fa, ra_rc, fb) = (
        p[0].seq.clone(),
        reverse_complement(&p[1].seq),
        p[2].seq.clone(),
    );
    let read = joined(
        &fa,
        &[
            &splitmix_dna(129, 500),
            &fb,
            &splitmix_dna(130, 50),
            marker,
            &splitmix_dna(131, 500),
        ],
        &ra_rc,
    );
    let n = read.len();
    assert_eq!(
        adapter_segments_annotated(&read, &c),
        vec![seg(
            fa.len(),
            n - ra_rc.len(),
            locus(0, fa.len(), 2, false),
            locus(n - ra_rc.len(), n, 3, true)
        )]
    );
}

/// A sheet primer that maps to the reverse complement of an entry
/// (`split_opens` `Reversed`) locates a primer only in the orientation valid for
/// its end: as ordered at the 5' end it trims, and the entry's own sequence
/// there, which closes an amplicon, neither trims nor becomes a locus.
#[test]
fn reverse_complement_entry_is_gated_by_orientation() {
    let p = split_primers();
    let ra_rc = reverse_complement(&p[1].seq);
    let mut c = cfg_with(vec![entry("rA_rc", &ra_rc, Role::Barcode)], 0.2, 150, true);
    c.attach_split(&p);
    assert_eq!(c.split_opens[0], Opens::Reversed);
    let fa_rc = reverse_complement(&p[0].seq);
    let insert = splitmix_dna(132, 800);
    let (r, f) = (ra_rc.len(), fa_rc.len());

    let valid = joined(&p[1].seq, &[&insert], &fa_rc);
    let n = valid.len();
    assert_eq!(
        adapter_segments_annotated(&valid, &c),
        vec![seg(
            r,
            n - f,
            locus(0, r, 0, true),
            locus(n - f, n, 1, true)
        )]
    );

    let wrong = joined(&ra_rc, &[&insert], &fa_rc);
    let n = wrong.len();
    assert_eq!(
        adapter_segments_annotated(&wrong, &c),
        vec![seg(0, n - f, None, locus(n - f, n, 1, true))]
    );
}

/// An entry that fuses an adapter with a sheet primer trims the 5' end as
/// given, and the sheet primer within its trim is the locus there.
#[test]
fn fused_adapter_primer_entry_yields_the_sheet_primer_locus() {
    let p = split_primers();
    let lsk = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let fused = [lsk.as_slice(), &p[0].seq].concat();
    let mut c = cfg_with(vec![ad("lsk_fA", &fused)], 0.2, 150, true);
    c.attach_split(&p);
    assert!(!c.is_split(0));
    let ra_rc = reverse_complement(&p[1].seq);
    let read = joined(&fused, &[&splitmix_dna(135, 800)], &ra_rc);
    let n = read.len();
    let (a, f, r) = (lsk.len(), p[0].seq.len(), ra_rc.len());
    assert_eq!(
        adapter_segments_annotated(&read, &c),
        vec![seg(
            a + f,
            n - r,
            locus(a, a + f, 1, false),
            locus(n - r, n, 2, true)
        )]
    );
}

/// A sheet primer that overlaps a longer entry of the set at the 3' end is
/// the locus there, and the trim stays where the longer entry set it.
#[test]
fn overlapping_catalog_entry_keeps_its_trim_and_yields_the_locus() {
    let p = split_primers();
    let ra_rc = reverse_complement(&p[1].seq);
    let longer = [ra_rc.as_slice(), &splitmix_dna(136, 6)].concat();
    let mut c = cfg_with(
        vec![entry("rA_long", &longer, Role::Primer)],
        0.2,
        150,
        true,
    );
    c.attach_split(&p);
    let read = joined(&p[0].seq, &[&splitmix_dna(137, 800)], &longer);
    let n = read.len();
    let (f, l) = (p[0].seq.len(), longer.len());
    assert_eq!(
        adapter_segments_annotated(&read, &c),
        vec![seg(
            f,
            n - l,
            locus(0, f, 1, false),
            locus(n - l, n - l + ra_rc.len(), 2, true)
        )]
    );
}

/// A remnant of a sheet primer in the adapter role behind unalignable bases,
/// which only the residue search finds, trims the end and is the locus
/// there, its outer edge a read end.
#[test]
fn residue_trimmed_split_adapter_yields_an_open_locus() {
    let p = split_primers();
    let mut c = cfg_with(vec![ad("fA", &p[0].seq)], 0.2, 150, true);
    c.attach_split(&p);
    assert!(c.is_split(0));
    let junk = splitmix_dna(139, 12);
    let remnant = &p[0].seq[8..];
    let ra_rc = reverse_complement(&p[1].seq);
    let read = joined(&junk, &[remnant, &splitmix_dna(142, 800)], &ra_rc);
    let n = read.len();
    let (j, r) = (junk.len(), ra_rc.len());
    let open = Some(Locus {
        start: j,
        end: j + remnant.len(),
        outer_open: true,
        boundary: false,
        site: None,
    });
    assert_eq!(
        adapter_segments_annotated(&read, &c),
        vec![seg(
            j + remnant.len(),
            n - r,
            open,
            locus(n - r, n, 1, true)
        )]
    );
}

/// An entry that serves a primer and its reverse complement, as a sheet
/// whose reverse primer is the reverse complement of its forward primer
/// gives, reads into the insert on both strands: its hits are located at
/// both ends of a read on either strand.
#[test]
fn entry_of_a_primer_and_its_reverse_complement_opens_both_ways() {
    let px = splitmix_dna(140, 22);
    let primers = vec![
        crate::split::Primer {
            name: "fX".into(),
            seq: px.clone(),
        },
        crate::split::Primer {
            name: "rX".into(),
            seq: reverse_complement(&px),
        },
    ];
    let mut c = cfg_with(
        vec![ad("lsk", b"AATGTACTTCGTTCAGTTACGTATTGCT")],
        0.2,
        150,
        true,
    );
    c.attach_split(&primers);
    assert_eq!(c.adapters.len(), 2);
    assert_eq!(c.split_opens[1], Opens::Both);
    assert_eq!(CandidateIndex::for_config(&c).opens[1], Opens::Both);
    let insert = splitmix_dna(141, 800);
    let x = px.len();
    for (ends, rc) in [(px.clone(), false), (reverse_complement(&px), true)] {
        let read = joined(&ends, &[&insert], &ends);
        let n = read.len();
        assert_eq!(
            adapter_segments_annotated(&read, &c),
            vec![seg(x, n - x, locus(0, x, 1, rc), locus(n - x, n, 1, rc))]
        );
    }
}

/// A valid sheet primer hit that overlaps the outermost one, or starts
/// within `FLANK_SLACK` bases of its end, moves the trim to its own end;
/// the outermost hit stays the locus. One starting further away stays in
/// the read.
#[test]
fn abutting_inner_primer_moves_the_trim_but_not_the_locus() {
    let p = split_primers();
    let (fa, fb) = (p[0].seq.clone(), p[2].seq.clone());
    let ra_rc = reverse_complement(&p[1].seq);
    let insert = splitmix_dna(133, 800);
    for gap in [0, 5, FLANK_SLACK, FLANK_SLACK + 1] {
        let read = joined(&fa, &[&splitmix_dna(134, gap), &fb, &insert], &ra_rc);
        let n = read.len();
        let start = if gap <= FLANK_SLACK {
            fa.len() + gap + fb.len()
        } else {
            fa.len()
        };
        assert_eq!(
            adapter_segments_annotated(&read, &split_cfg()),
            vec![seg(
                start,
                n - ra_rc.len(),
                locus(0, fa.len(), 1, false),
                locus(n - ra_rc.len(), n, 2, true)
            )],
            "gap {gap}"
        );
    }
}

/// A whole sheet primer hit above the `k_far` budget of its entry is placed
/// only when anchored: at the read start, or directly behind an adapter, it
/// trims and is the locus; 80 bases into the read it neither trims nor
/// locates a primer. A sheet of 24 primers of 16 bases shares one terminal
/// chance bound, which admits three edits anchored and two anywhere in the
/// end zone.
#[test]
fn far_primer_hit_is_placed_only_when_anchored() {
    let primers: Vec<crate::split::Primer> = (1..=24)
        .map(|i| crate::split::Primer {
            name: format!("p{i}"),
            seq: splitmix_dna(i, 16),
        })
        .collect();
    let adapter = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let mut c = cfg_with(vec![ad("lsk", adapter)], 0.2, 150, true);
    c.attach_split(&primers);
    let Budget { k_end, k_far, .. } = CandidateIndex::for_config(&c).budgets[5];
    assert_eq!((k_end, k_far), (3, 2));
    let marginal = substituted(&primers[4].seq, &[3, 8, 13]);
    let marginal_site = PrimerSite {
        entry: 5,
        rc: false,
        cost: 3,
    };
    let insert = splitmix_dna(9001, 2000);
    let n = insert.len() + marginal.len();

    let at_end = [marginal.clone(), insert.clone()].concat();
    assert_eq!(
        adapter_segments_annotated(&at_end, &c),
        vec![seg(16, n, sited(0, 16, marginal_site), None)]
    );

    let behind_adapter = [adapter.to_vec(), marginal.clone(), insert.clone()].concat();
    assert_eq!(
        adapter_segments_annotated(&behind_adapter, &c),
        vec![seg(
            44,
            n + adapter.len(),
            sited(28, 44, marginal_site),
            None
        )]
    );

    let deep = [insert[..80].to_vec(), marginal, insert[80..].to_vec()].concat();
    assert_eq!(
        adapter_segments_annotated(&deep, &c),
        vec![seg(0, n, None, None)]
    );
}

/// A held sheet primer hit with fewer than `MIN_OVERLAP` bases outboard of
/// it holds the read end: with `MIN_OVERLAP - 1` bases of another primer
/// before `fA`, `fA` is the locus. With `MIN_OVERLAP` bases, the read end is
/// left to the partial search, which locates the other primer cut short
/// there as an open locus.
#[test]
fn held_primer_holds_the_read_end_below_the_minimum_overlap() {
    let p = split_primers();
    let (fa, fb) = (p[0].seq.clone(), p[2].seq.clone());
    let ra_rc = reverse_complement(&p[1].seq);
    for outboard in [MIN_OVERLAP - 1, MIN_OVERLAP] {
        let cut = &fb[fb.len() - outboard..];
        let read = joined(cut, &[&fa, &splitmix_dna(155, 800)], &ra_rc);
        let segments = adapter_segments_annotated(&read, &split_cfg());
        assert_eq!(segments.len(), 1, "{outboard}: {segments:?}");
        let five = segments[0].five.expect("a primer is located");
        if outboard < MIN_OVERLAP {
            assert_eq!(
                (five.start, five.end, five.outer_open),
                (outboard, outboard + fa.len(), false),
                "{outboard}: {segments:?}"
            );
        } else {
            assert_eq!(
                (five.start, five.end, five.outer_open),
                (0, outboard, true),
                "{outboard}: {segments:?}"
            );
        }
    }
}

/// A sheet primer site deep in an end zone leaves the read end to the
/// partial search: the primer cut short at each read end is the locus there
/// and trims, and the whole site of the other target's primer further inward
/// stays in the read.
#[test]
fn deep_sheet_primer_site_leaves_the_read_end_to_the_partial_search() {
    let p = split_primers();
    let (fa, fb) = (p[0].seq.clone(), p[2].seq.clone());
    let (ra_rc, rb_rc) = (reverse_complement(&p[1].seq), reverse_complement(&p[3].seq));
    let cut = 6;
    let (fa_cut, ra_cut) = (&fa[cut..], &ra_rc[..ra_rc.len() - cut]);
    let read = joined(
        fa_cut,
        &[
            &splitmix_dna(150, 100),
            &fb,
            &splitmix_dna(151, 800),
            &rb_rc,
            &splitmix_dna(152, 100),
        ],
        ra_cut,
    );
    let n = read.len();
    let (f, r) = (fa_cut.len(), ra_cut.len());
    let open = |start, end| {
        Some(Locus {
            start,
            end,
            outer_open: true,
            boundary: false,
            site: None,
        })
    };
    assert_eq!(
        adapter_segments_annotated(&read, &split_cfg()),
        vec![seg(f, n - r, open(0, f), open(n - r, n))]
    );
}

/// An adapter-role entry of a primer variant, as an amplicon kit lists it,
/// cut short at the 3' read end trims there even when a whole sheet primer
/// site lies deeper in the end zone, and the sheet primer it covers is the
/// locus of that end.
#[test]
fn kit_primer_cut_short_at_the_read_end_covers_the_sheet_primer() {
    let p = split_primers();
    let ra_rc = reverse_complement(&p[1].seq);
    let rb_rc = reverse_complement(&p[3].seq);
    let kit = substituted(&p[1].seq, &[0]);
    let mut c = cfg_with(
        vec![
            ad("lsk", b"AATGTACTTCGTTCAGTTACGTATTGCT"),
            ad("kit_rA", &kit),
        ],
        0.2,
        150,
        true,
    );
    c.attach_split(&p);
    assert!(!c.is_split(1));
    let cut = 6;
    let ra_cut = &ra_rc[..ra_rc.len() - cut];
    let fa = &p[0].seq;
    let read = joined(
        fa,
        &[&splitmix_dna(153, 800), &rb_rc, &splitmix_dna(154, 100)],
        ra_cut,
    );
    let n = read.len();
    let (f, r) = (fa.len(), ra_cut.len());
    assert_eq!(
        adapter_segments_annotated(&read, &c),
        vec![seg(
            f,
            n - r,
            locus(0, f, 2, false),
            Some(Locus {
                start: n - r,
                end: n,
                outer_open: true,
                boundary: false,
                site: None,
            })
        )]
    );
}

/// A read with no primer keeps its whole span and reports no locus.
#[test]
fn no_locus_without_primer() {
    let read = splitmix_dna(125, 800);
    assert_eq!(
        adapter_segments_annotated(&read, &split_cfg()),
        vec![seg(0, read.len(), None, None)]
    );
}

/// Sheet primers already in the adapter set, as the MAB114 preset holds
/// them, map to their entries: the set keeps its size, names and roles.
/// Attaching a sheet leaves the amplicon judgement as it was, and a split
/// entry gets no whole primer even in an amplicon library, since it splits a
/// read only at a junction pair.
#[test]
fn attach_split_maps_preset_duplicate() {
    let mut c = cfg_with(
        super::preset::preset(&[super::preset::Kit::Mab114]),
        0.2,
        150,
        true,
    );
    let before = c.adapters.clone();
    let sheet = crate::split::Sheet::parse_fasta(
        b">16S_mix_F target=16S end=fwd\nAGRGTTYGATYMTGGCTCAG\n\
          >16S_mix_R target=16S end=rev\nSGGYTACCTTGTTACGACTT\n\
          >ITS1 target=ITS end=fwd\nTCCGTAGGTGAACCTGCGG\n\
          >ITS4 target=ITS end=rev\nTCCTCCGCTTATTGATATGC\n",
    )
    .unwrap();
    c.attach_split(&sheet.primers);
    assert_eq!(c.adapters, before);
    assert!(!c.amplicon);
    for (primer_idx, primer) in sheet.primers.iter().enumerate() {
        let i = c
            .adapters
            .iter()
            .position(|a| a.name == primer.name)
            .unwrap();
        assert_eq!(c.split_of[i], Some(primer_idx), "{}", primer.name);
    }
    assert_eq!(c.split_of.iter().flatten().count(), sheet.primers.len());
    let index = CandidateIndex::for_config(&c);
    let i = c
        .adapters
        .iter()
        .position(|a| a.name == "16S_mix_F")
        .unwrap();
    assert!(index.paired[i] && index.opens[i] == Opens::AsGiven);
    assert!(index.whole_primers[i].is_none());
    c.set_amplicon(true);
    let index = CandidateIndex::for_config(&c);
    assert!(index.whole_primers.iter().all(Option::is_none));
}

/// A sheet primer absent from the adapter set is appended in the primer
/// role under its sheet name.
#[test]
fn attach_split_appends_new_primer() {
    let c = split_cfg();
    let p = split_primers();
    assert_eq!(c.adapters.len(), 1 + p.len());
    assert_eq!(c.split_of[0], None);
    for (primer_idx, primer) in p.iter().enumerate() {
        let entry = &c.adapters[1 + primer_idx];
        assert_eq!(
            (entry.name.as_str(), entry.seq.as_slice(), entry.role),
            (primer.name.as_str(), primer.seq.as_slice(), Role::Primer)
        );
        assert_eq!(c.split_of[1 + primer_idx], Some(primer_idx));
    }
    let index = CandidateIndex::for_config(&c);
    assert!((1..c.adapters.len()).all(|i| index.paired[i] && index.opens[i] == Opens::AsGiven));
    assert!(!index.paired[0]);
}

/// A sheet primer whose reverse complement is already in the set maps to
/// that entry, which reads out of the insert, and still locates the primer.
#[test]
fn attach_split_maps_reverse_complement_entry() {
    let p = split_primers();
    let ra_rc = reverse_complement(&p[1].seq);
    let mut c = cfg_with(vec![entry("rA_rc", &ra_rc, Role::Barcode)], 0.2, 150, true);
    c.attach_split(&p);
    assert_eq!(c.adapters.len(), p.len());
    assert_eq!(
        (c.adapters[0].name.as_str(), c.adapters[0].role),
        ("rA_rc", Role::Barcode)
    );
    assert_eq!(c.split_of[0], Some(1));
    let index = CandidateIndex::for_config(&c);
    assert!(index.paired[0] && index.opens[0] == Opens::Reversed);
    let read = joined(&p[0].seq, &[&splitmix_dna(126, 800)], &ra_rc);
    let n = read.len();
    let (f, r) = (p[0].seq.len(), ra_rc.len());
    assert_eq!(
        adapter_segments_annotated(&read, &c),
        vec![seg(
            f,
            n - r,
            locus(0, f, 1, false),
            locus(n - r, n, 0, false)
        )]
    );
}

/// Replacing the adapter set detaches the split primers.
#[test]
fn replace_adapters_clears_split() {
    let mut c = split_cfg();
    c.replace_adapters(vec![ad("lsk", b"AATGTACTTCGTTCAGTTACGTATTGCT")]);
    assert!(c.split_of.is_empty() && c.split_opens.is_empty());
}

/// Returns `cfg` with an index whose singleton batches are all searched
/// tiled (`tiled`) or all one pattern at a time, whatever the read.
fn with_singleton_batches(cfg: &AdapterConfig, tiled: bool) -> AdapterConfig {
    let mut index = CandidateIndex::for_config(cfg);
    for batch in &mut index.singleton_batches {
        batch.tiled_on_plain_read = tiled;
        batch.tiled_on_ambiguous_read = tiled;
    }
    AdapterConfig {
        candidate_index: std::sync::OnceLock::from(index),
        ..cfg.clone()
    }
}

/// Amplicon reads of the MAB114 kit for the batch comparison: SplitMix
/// inserts between every forward and reverse primer variant, with the
/// ambiguity codes resolved, behind a barcode construct or bare, with edits
/// in the primers, in both orientations, some cut short at an end, some
/// joined into chimeras and some holding an `N`.
fn planted_primer_reads(count: usize) -> Vec<Vec<u8>> {
    use super::catalog::sequence_of;
    let forward = [
        "16S_mix_F",
        "16S_Bor_F",
        "16S_Chl_F",
        "16S_Ent_F",
        "ITS1",
        "ITS1_Fus",
        "ITS1_Mal",
    ];
    let reverse = ["16S_mix_R", "16S_Bor_R", "16S_Chl_R", "ITS4", "ITS4_Pyt"];
    let mut rng = Lcg(114);
    let planted = |name: &str, rng: &mut Lcg| -> Vec<u8> {
        let mut seq: Vec<u8> = sequence_of(name)
            .unwrap()
            .iter()
            .map(|&code| {
                let bases = iupac_bases(code).unwrap();
                bases[rng.below(bases.len())]
            })
            .collect();
        for _ in 0..rng.below(4) {
            let at = rng.below(seq.len());
            match rng.below(3) {
                0 => seq[at] = b"ACGT"[rng.below(4)],
                1 => {
                    seq.remove(at);
                },
                _ => seq.insert(at, b"ACGT"[rng.below(4)]),
            }
        }
        seq
    };
    let amplicon = |seed: u64, rng: &mut Lcg| -> Vec<u8> {
        let construct: &[u8] = if rng.below(3) == 0 {
            b""
        } else {
            b"GCTTGGGTGTTTAACCGCACCTGGAACTTGTGCCTTCCACCCATATCCGTGTCGCCCTT"
        };
        let head = [construct, &planted(forward[rng.below(forward.len())], rng)].concat();
        let tail = [construct, &planted(reverse[rng.below(reverse.len())], rng)].concat();
        let insert = splitmix_dna(40_000 + seed, 200 + rng.below(400));
        let read = [head, insert, reverse_complement(&tail)].concat();
        if rng.below(2) == 0 {
            read
        } else {
            reverse_complement(&read)
        }
    };
    (0..count)
        .map(|i| {
            let mut read = amplicon(2 * i as u64, &mut rng);
            if rng.below(8) == 0 {
                read.extend(amplicon(2 * i as u64 + 1, &mut rng));
            }
            if rng.below(6) == 0 {
                let cut = rng.below(25);
                read.drain(..cut);
                read.truncate(read.len() - rng.below(25));
            }
            if rng.below(7) == 0 {
                let at = rng.below(read.len());
                read[at] = b'N';
            }
            read
        })
        .collect()
}

/// Some batch takes the one-at-a-time path (`TerminalBatch::tiled`) on a
/// plain read, under the kit, amplicon, split and ends-only configs.
#[test]
fn tiled_singleton_batches_match_the_one_by_one_search() {
    use super::preset::{Kit, preset};
    let kit = cfg_with(preset(&[Kit::Mab114]), 0.2, 150, true);
    let mut amplicon = kit.clone();
    amplicon.set_amplicon(true);
    let mut split = kit.clone();
    split.attach_split(&crate::split::Sheet::preset("mab114").unwrap().primers);
    let mut ends_only = kit.clone();
    ends_only.split = false;
    ends_only.candidate_index = std::sync::OnceLock::new();
    let reads = planted_primer_reads(300);
    for (label, c) in [
        ("kit", kit),
        ("amplicon", amplicon),
        ("split", split),
        ("ends only", ends_only),
    ] {
        let index = CandidateIndex::for_config(&c);
        assert!(
            index.singleton_batches.iter().any(|b| b.tiled(false)),
            "{label}: no batch is tiled on a read with an N"
        );
        assert!(
            index.singleton_batches.iter().any(|b| !b.tiled(true)),
            "{label}: no batch keeps the DNA profile on a plain read"
        );
        let one_by_one = with_singleton_batches(&c, false);
        let tiled = with_singleton_batches(&c, true);
        let run = |cfg: &AdapterConfig, read: &[u8]| {
            let mut acted = vec![false; cfg.adapters.len()];
            adapter_segments_tallied(read, cfg, &mut acted);
            (adapter_segments_annotated(read, cfg), acted)
        };
        let mut trimmed = 0;
        for (i, read) in reads.iter().enumerate() {
            let expected = run(&one_by_one, read);
            assert_eq!(run(&tiled, read), expected, "{label}: read {i}, tiled");
            assert_eq!(run(&c, read), expected, "{label}: read {i}, default");
            trimmed += usize::from(expected.0 != [Segment::located(0, read.len(), &[])]);
        }
        assert!(
            trimmed > reads.len() / 2,
            "{label}: {trimmed} reads trimmed"
        );
    }
}

/// Returns every singleton's whole-pattern hits over the end windows of
/// `span` of `read` under `cfg`, as `singleton_hits` reports them.
fn singleton_hit_lists(
    read: &[u8],
    cfg: &AdapterConfig,
    span: (usize, usize),
) -> Vec<(usize, Vec<(Site, search::Hit)>)> {
    let mut lists = Vec::new();
    super::passes::with_engine(read, cfg, |ctx, engine| {
        super::passes::singleton_hits(ctx, span, cfg.end_size, engine, |adapter_idx, found| {
            lists.push((adapter_idx, found.to_vec()))
        });
    });
    lists
}

/// Reads for the hit-level comparison of the tiled and one-by-one terminal
/// searches under `cfg`: the entries of the set that a singleton batch holds,
/// planted with edits at the ends and inside, with indels beside
/// homopolymer runs, with `N` bases, cut to spans shorter than a pattern,
/// joined into chimeras, and SplitMix reads without a planted entry.
fn batch_comparison_reads(cfg: &AdapterConfig, random: usize) -> Vec<Vec<u8>> {
    let index = CandidateIndex::for_config(cfg);
    let batched: Vec<&[u8]> = index
        .singleton_batches
        .iter()
        .flat_map(|batch| batch.adapter_indices.iter())
        .map(|&idx| cfg.adapters[idx].seq.as_slice())
        .collect();
    assert!(!batched.is_empty());
    let mut rng = Lcg(0x7469_6c65_6421);
    let resolved = |seq: &[u8], rng: &mut Lcg| -> Vec<u8> {
        seq.iter()
            .map(|&code| {
                let bases = iupac_bases(code).unwrap();
                bases[rng.below(bases.len())]
            })
            .collect()
    };
    let edited = |mut seq: Vec<u8>, rng: &mut Lcg| -> Vec<u8> {
        for _ in 0..rng.below(4) {
            let at = rng.below(seq.len());
            match rng.below(3) {
                0 => seq[at] = b"ACGT"[rng.below(4)],
                1 => {
                    seq.remove(at);
                },
                _ => seq.insert(at, b"ACGT"[rng.below(4)]),
            }
        }
        seq
    };
    // An insertion or deletion of one base of a homopolymer run, where the
    // alignment may place the edit anywhere along the run.
    let homopolymer_indel = |mut seq: Vec<u8>, rng: &mut Lcg| -> Vec<u8> {
        let runs: Vec<usize> = (1..seq.len()).filter(|&i| seq[i] == seq[i - 1]).collect();
        let at = if runs.is_empty() {
            let at = rng.below(seq.len());
            seq.insert(at, seq[at]);
            at
        } else {
            runs[rng.below(runs.len())]
        };
        if rng.below(2) == 0 {
            seq.insert(at, seq[at]);
        } else {
            seq.remove(at);
        }
        seq
    };
    let planted = |rng: &mut Lcg| -> Vec<u8> {
        let entry = resolved(batched[rng.below(batched.len())], rng);
        let entry = match rng.below(3) {
            0 => entry,
            1 => edited(entry, rng),
            _ => homopolymer_indel(entry, rng),
        };
        if rng.below(2) == 0 {
            entry
        } else {
            reverse_complement(&entry)
        }
    };
    let flank = |rng: &mut Lcg| -> Vec<u8> {
        let len = rng.below(40);
        if rng.below(3) == 0 {
            vec![b"ACGT"[rng.below(4)]; len]
        } else {
            rng.dna(len)
        }
    };
    let mut reads = Vec::new();
    for i in 0..600u64 {
        let mut read = [
            flank(&mut rng),
            planted(&mut rng),
            splitmix_dna(70_000 + i, 20 + rng.below(300)),
            planted(&mut rng),
            flank(&mut rng),
        ]
        .concat();
        if rng.below(4) == 0 {
            let at = 150 + rng.below(read.len().saturating_sub(150).max(1));
            let inner = planted(&mut rng);
            let at = at.min(read.len());
            read.splice(at..at, inner);
        }
        if rng.below(5) == 0 {
            let other = [
                planted(&mut rng),
                splitmix_dna(80_000 + i, 200),
                planted(&mut rng),
            ]
            .concat();
            read.extend(other);
        }
        if rng.below(4) == 0 {
            for _ in 0..1 + rng.below(3) {
                let at = rng.below(read.len());
                read[at] = b'N';
            }
        }
        reads.push(read);
    }
    for _ in 0..200 {
        let entry = planted(&mut rng);
        let len = 1 + rng.below(entry.len() + 8);
        let mut read = [entry, rng.dna(8)].concat();
        read.truncate(len);
        reads.push(read);
    }
    for i in 0..random as u64 {
        reads.push(splitmix_dna(90_000 + i, 1 + rng.below(500)));
    }
    reads
}

/// The tiled search of a singleton batch reports, for every adapter it
/// holds, the hits and sites of the one-by-one search in the same order:
/// under the `all` and `ont` presets, in an amplicon library, over whole
/// reads and over spans of every length, on the reads of
/// `batch_comparison_reads`.
#[test]
fn tiled_singleton_batches_report_the_hits_of_the_one_by_one_search() {
    use super::preset::{Kit, preset};
    let mut amplicon = cfg_with(preset(&[Kit::Mab114]), 0.2, 150, true);
    amplicon.set_amplicon(true);
    for (label, c) in [
        ("all", cfg_with(preset(Kit::ALL), 0.2, 150, true)),
        ("ont", cfg_with(preset(Kit::ONT), 0.2, 150, true)),
        ("mab114", amplicon),
    ] {
        let one_by_one = with_singleton_batches(&c, false);
        let tiled = with_singleton_batches(&c, true);
        let reads = batch_comparison_reads(&c, 3000);
        let batched = CandidateIndex::for_config(&c).singleton_batch_of;
        let mut rng = Lcg(0x7370_616e);
        let (mut spans, mut hits) = (0, 0);
        for (i, read) in reads.iter().enumerate() {
            let n = read.len();
            let mut read_spans = vec![(0, n)];
            for _ in 0..2 {
                let start = rng.below(n + 1);
                let end = start + rng.below(n - start + 1);
                read_spans.push((start, end));
            }
            for span in read_spans {
                if span.0 == span.1 {
                    continue;
                }
                let expected = singleton_hit_lists(read, &one_by_one, span);
                assert_eq!(
                    singleton_hit_lists(read, &tiled, span),
                    expected,
                    "{label}: read {i}, span {span:?}"
                );
                spans += 1;
                hits += expected
                    .iter()
                    .filter(|(adapter_idx, _)| batched[*adapter_idx].is_some())
                    .map(|(_, found)| found.len())
                    .sum::<usize>();
            }
        }
        assert!(
            spans > reads.len() && hits > reads.len(),
            "{label}: {hits} hits of batched entries"
        );
    }
}
