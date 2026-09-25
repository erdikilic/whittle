use super::*;

fn random_bases(mut state: u64, len: usize) -> Vec<u8> {
    (0..len)
        .map(|_| {
            state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^= z >> 31;
            b"ACGT"[(z >> 62) as usize]
        })
        .collect()
}

fn infer_owned(reads: &[Vec<u8>]) -> Vec<InferredAdapter> {
    infer_with_known(reads, Vec::new())
}

fn infer_with_known(reads: &[Vec<u8>], known: Vec<Adapter>) -> Vec<InferredAdapter> {
    let sample: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    discover(
        &sample,
        &AdapterConfig {
            adapters: known,
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 20,
            candidate_index: std::sync::OnceLock::new(),
        },
    )
}

#[test]
fn recovers_three_distinct_adapters_without_length_cap() {
    let adapters: Vec<Vec<u8>> = [27, 43, 71]
        .into_iter()
        .enumerate()
        .map(|(i, len)| random_bases(313 + i as u64, len))
        .collect();
    let reads: Vec<Vec<u8>> = (0..900)
        .map(|i| {
            let mut read = adapters[i % 3].clone();
            read.extend(random_bases(4321 + i as u64, 250));
            read
        })
        .collect();
    let found = infer_owned(&reads);
    for seq in &adapters {
        assert!(
            found.iter().any(|d| d.adapter.seq == *seq),
            "missing {}: {found:?}",
            String::from_utf8_lossy(seq)
        );
    }
    assert_eq!(found.len(), 3);
}

#[test]
fn recovers_minority_adapter_beside_a_dominant_family() {
    let common = random_bases(47319, 31);
    let rare = random_bases(18971, 43);
    let reads: Vec<Vec<u8>> = (0..2000)
        .map(|i| {
            let mut read = match i % 100 {
                0..=69 => common.clone(),
                70..=74 => rare.clone(),
                _ => Vec::new(),
            };
            read.extend(random_bases(18231 + i, 420));
            read
        })
        .collect();
    let found = infer_owned(&reads);
    assert!(found.iter().any(|d| d.adapter.seq == common), "{found:?}");
    assert!(found.iter().any(|d| d.adapter.seq == rare), "{found:?}");
    assert_eq!(found.len(), 2, "{found:?}");
}

#[test]
fn tandem_repeat_ends_do_not_produce_adapters() {
    let reads: Vec<Vec<u8>> = (0..600)
        .map(|i| {
            let mut read = b"TTAGGG".repeat(10);
            if i % 3 == 0 {
                read[19] = b'C';
            }
            read.extend(random_bases(87123 + i, 400));
            read
        })
        .collect();
    let found = infer_owned(&reads);
    assert!(found.is_empty(), "{found:?}");
}

#[test]
fn minority_primer_uses_unprimed_insert_boundaries() {
    let primer = random_bases(88971, 24);
    let anchor = random_bases(34723, 140);
    let reads: Vec<Vec<u8>> = (0..2000)
        .map(|i| {
            let mut read = if i % 20 == 0 {
                primer.clone()
            } else {
                Vec::new()
            };
            read.extend_from_slice(&anchor);
            read.extend(random_bases(8921 + i, 300));
            read
        })
        .collect();
    let found = infer_owned(&reads);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0].adapter.seq, primer);
}

#[test]
fn kmer_support_counts_independent_windows() {
    let seed = b"ACGTCAGTGCATGACT";
    let repeated = seed.repeat(5);
    let windows = vec![repeated.as_slice(), seed.as_slice()];
    let counts = top_kmers(&windows, 16, 500);
    assert_eq!(
        counts
            .iter()
            .find(|(key, _)| *key == encode_kmer(seed).unwrap())
            .unwrap()
            .1,
        2
    );
    assert!(top_kmers(&[b"TTAGGGTTAGGGTTAGGG"], 16, 500).is_empty());
}

#[test]
fn clean_conserved_inserts_are_not_adapters() {
    let prefix = random_bases(7877, 140);
    let suffix = random_bases(4512, 140);
    let reads: Vec<Vec<u8>> = (0..400)
        .map(|i| {
            let mut read = prefix.clone();
            read.extend(random_bases(915 + i, 200));
            read.extend_from_slice(&suffix);
            read
        })
        .collect();
    assert!(infer_owned(&reads).is_empty());
}

#[test]
fn variable_conserved_inserts_extend_beyond_assembly_windows() {
    let anchors: Vec<_> = (0..4)
        .map(|i| (random_bases(7877 + i, 120), random_bases(4512 + i, 120)))
        .collect();
    let reads: Vec<Vec<u8>> = (0..2000)
        .map(|i| {
            let (prefix, suffix) = &anchors[i % anchors.len()];
            let mut read = prefix.clone();
            read.extend(random_bases(915 + i as u64, 300));
            read.extend_from_slice(suffix);
            let mutations = random_bases(54371 + i as u64, read.len() * 5);
            let mut mutated = Vec::new();
            for (&base, event) in read.iter().zip(mutations.as_chunks::<5>().0) {
                match encode_kmer(&event[..4]).unwrap() {
                    0..=4 => mutated.push(event[4]),
                    5..=6 => {},
                    7 => mutated.extend_from_slice(&[base, event[4]]),
                    _ => mutated.push(base),
                }
            }
            mutated
        })
        .collect();
    let found = infer_owned(&reads);
    assert!(found.is_empty(), "{found:?}");
}

#[test]
fn recovers_degenerate_primers_at_conserved_insert_boundaries() {
    let front = b"TCGATGARYCTACGTGACCT";
    let rear = b"GCTAGTACCGATGCTAGTCA";
    let prefix = random_bases(712, 140);
    let suffix = random_bases(815, 140);
    let reads: Vec<Vec<u8>> = (0..500usize)
        .map(|i| {
            let mut read = Vec::new();
            if i % 5 != 0 {
                read.extend_from_slice(front);
                read[7] = b"AG"[(i / 5) % 2];
                read[8] = b"CT"[(i / 10) % 2];
            }
            read.extend_from_slice(&prefix);
            read.extend(random_bases(9712 + i as u64, 200));
            read.extend_from_slice(&suffix);
            if i % 5 != 0 {
                read.extend_from_slice(rear);
            }
            read
        })
        .collect();
    let found = infer_owned(&reads);
    assert!(found.iter().any(|d| d.adapter.seq == front), "{found:?}");
    assert!(found.iter().any(|d| d.adapter.seq == rear), "{found:?}");
    assert_eq!(found.len(), 2, "{found:?}");
}

#[test]
fn ambiguity_codes_cover_exactly_the_voted_bases() {
    for mask in 1..16 {
        let bases = crate::adapter::search::iupac_bases(ambiguity_code(mask)).unwrap();
        let actual = bases
            .iter()
            .fold(0usize, |m, &b| m | (1 << encode_kmer(&[b]).unwrap()));
        assert_eq!(actual, mask);
    }
}

#[test]
fn distinct_primers_at_one_insert_boundary_remain_distinct() {
    let primers = [random_bases(7651, 22), random_bases(2157, 24)];
    let insert = random_bases(2871, 140);
    let reads: Vec<Vec<u8>> = (0..600)
        .map(|i| {
            let mut read = Vec::new();
            if i % 5 != 0 {
                read.extend_from_slice(&primers[i % 2]);
            }
            read.extend_from_slice(&insert);
            read.extend(random_bases(159 + i as u64, 200));
            read
        })
        .collect();
    let found = infer_owned(&reads);
    for primer in &primers {
        assert!(found.iter().any(|d| d.adapter.seq == *primer), "{found:?}");
    }
    assert_eq!(found.len(), 2);
}

#[test]
fn truncated_insert_anchor_recovers_distinct_upstream_primers() {
    let primers = [random_bases(7651, 22), random_bases(2157, 24)];
    let insert = random_bases(2871, 140);
    let reads: Vec<Vec<u8>> = (0..600)
        .map(|i| {
            let mut read = Vec::new();
            if i % 5 != 0 {
                read.extend_from_slice(&primers[i % 2]);
            }
            read.extend_from_slice(&insert);
            read.truncate(WINDOW_LEN);
            read
        })
        .collect();
    let windows: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let found = contrast_boundary(&insert[1..90], &windows, End::Five)
        .unwrap()
        .primers;
    assert_eq!(found.len(), 2);
    for primer in &primers {
        assert!(found.contains(primer), "{found:?}");
    }
}

/// Encoding then decoding a k-mer is the identity.
#[test]
fn kmer_codec_roundtrips() {
    let k = b"ACGTACGTACGTACGT"; // 16bp
    let code = encode_kmer(k).unwrap();
    assert_eq!(decode_kmer(code, 16), k);
}

/// `encode_kmer` rejects ambiguity codes, lowercase bases and over-long
/// k-mers.
#[test]
fn encode_rejects_non_acgt() {
    assert_eq!(encode_kmer(b"ACGTN"), None);
    assert_eq!(encode_kmer(b"acgt"), None); // lowercase not accepted
    assert_eq!(encode_kmer(&[b'A'; 33]), None); // > 32 bases rejected
}

/// A short read yields itself as both windows and an empty read yields
/// nothing.
#[test]
fn end_windows_slices_both_ends() {
    let r1: &[u8] = b"AAAACCCCGGGGTTTTACGTACGT"; // 24bp
    let r2: &[u8] = b"TTTT"; // 4bp, below w: the whole read at both ends
    let sample: Vec<&[u8]> = vec![r1, r2, b""]; // empty skipped
    let (five, three) = end_windows(&sample, 8);
    assert_eq!(five, vec![&r1[..8], r2]); // first 8, then the whole short read
    assert_eq!(three, vec![&r1[16..], r2]); // last 8, then the whole short read
}

/// A window holding an ambiguity code is dropped from that end only.
#[test]
fn end_windows_drop_windows_holding_ambiguity_codes() {
    let r1: &[u8] = b"AAAANCCCGGGGTTTTACGTACGT"; // `N` in the 5' window only
    let r2: &[u8] = b"acgtacgtacgtacgtacgtacgn"; // `n` in the 3' window only
    let sample: Vec<&[u8]> = vec![r1, r2];
    let (five, three) = end_windows(&sample, 8);
    assert_eq!(five, vec![&r2[..8]]);
    assert_eq!(three, vec![&r1[16..]]);
}

/// A 16-mer planted in every window ranks first over unique filler.
#[test]
fn top_kmers_ranks_planted_over_background() {
    let planted = b"ACGTCAGTGCATGACT";
    let mut owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..50u8 {
        let mut wnd = planted.to_vec();
        // An uncalled separator prevents shifted copies of the seed.
        wnd.push(b'N');
        wnd.extend_from_slice(&[b'B' + (i % 4), b'C', b'G', b'T']);
        owned.push(wnd);
    }
    let windows: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let ranked = top_kmers(&windows, 16, 500);
    assert_eq!(decode_kmer(ranked[0].0, 16), planted);
    assert_eq!(ranked[0].1, 50);
}

/// A homopolymer window contributes no k-mer.
#[test]
fn top_kmers_drops_homopolymer() {
    let windows: Vec<&[u8]> = vec![b"AAAAAAAAAAAAAAAA"]; // pure homopolymer, 16bp
    assert!(
        top_kmers(&windows, 16, 500).is_empty(),
        "Low-complexity k-mer dropped"
    );
}

/// Each window counts once however often the k-mer occurs, and a
/// reverse-complement occurrence does not count.
#[test]
fn windows_containing_counts_windows_once_and_ignores_rc() {
    use crate::adapter::search::new_searcher_fwd;
    // The k-mer is not its own reverse complement, so the RC case is
    // meaningful: revcomp(AAAACCCCGGGGTATG) = CATACCCCGGGGTTTT.
    let kmer = b"AAAACCCCGGGGTATG"; // 16bp
    let w0v = b"TTAAAACCCCGGGGTATGTT".to_vec(); // exact occurrence
    let w1v = b"TTAAAACACCGGGGTATGTT".to_vec(); // 1 substitution (C to A)
    let mut w2v = b"AAAACCCCGGGGTATG".to_vec(); // k-mer twice; counted once
    w2v.extend_from_slice(b"GGGGAAAACCCCGGGGTATG");
    let w3v = b"TTCATACCCCGGGGTTTTTT".to_vec(); // reverse-complement only
    let windows: Vec<&[u8]> = vec![&w0v, &w1v, &w2v, &w3v];
    let mut s = new_searcher_fwd();
    // w0 (exact), w1 (1 edit) and w2 (twice, counted once) give 3; w3 (RC
    // only) is excluded.
    assert_eq!(windows_containing(&mut s, kmer, &windows, 2), 3);
}

/// Overlapping 4-mers that tile ACGTACG with descending weights along the
/// intended path reconstruct it: ACGT(9), CGTA(8), GTAC(7), TACG(6).
#[test]
fn bounded_heaviest_path_reconstructs_known_consensus() {
    let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
    let nodes = vec![
        mk(b"ACGT", 9),
        mk(b"CGTA", 8),
        mk(b"GTAC", 7),
        mk(b"TACG", 6),
    ];
    let (cons, profile, weight) = bounded_heaviest_path(&nodes, 4, 100).unwrap();
    assert_eq!(cons, b"ACGTACG"); // ACGT + C + A + G: 4 nodes give 7 nt
    assert_eq!(profile.len(), cons.len());
    assert_eq!(weight, 9 + 8 + 7 + 6);
}

/// On the 2-node cycle ATAT, TATA, ATAT (k = 4) the visited set stops the
/// walk after each node is used once, so the consensus is a short simple
/// path rather than a repeat filling `lmax`.
#[test]
fn bounded_heaviest_path_terminates_on_cycle() {
    let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
    let nodes = vec![mk(b"ATAT", 5), mk(b"TATA", 5)];
    let (cons, _profile, _w) = bounded_heaviest_path(&nodes, 4, 12).unwrap();
    assert!(cons.len() <= 12, "No loop: each node used at most once");
    assert!(cons.starts_with(b"ATAT") || cons.starts_with(b"TATA"));
}

/// Two non-overlapping tilings with different bases peel as two adapters.
#[test]
fn peel_extracts_two_distinct_adapters() {
    let mk = |s: &[u8], w: u32| (encode_kmer(s).unwrap(), w);
    let nodes = vec![
        // Adapter 1: ACGTACG..., high weight.
        mk(b"ACGT", 100),
        mk(b"CGTA", 99),
        mk(b"GTAC", 98),
        // Adapter 2: TTGGTTG..., lower weight but above 25% of 297.
        mk(b"TTGG", 90),
        mk(b"TGGT", 89),
        mk(b"GGTT", 88),
    ];
    let paths = peel_paths(nodes, 4, End::Five);
    assert_eq!(paths.len(), 2);
    assert!(paths[0].0.starts_with(b"ACGT"));
    assert!(paths[1].0.starts_with(b"TTGG"));
}

/// An exact catalog sequence is named at 100 percent identity.
#[test]
fn name_against_matches_catalog_entry() {
    let refs = vec![Adapter {
        name: "SQK-TEST".into(),
        seq: b"ACGTACGTACGTACGT".to_vec(),
        role: Role::Adapter,
    }];
    let hits = name_against(b"ACGTACGTACGTACGT", &refs, 0.2);
    assert_eq!(hits[0].0, "SQK-TEST");
    assert!((hits[0].1 - 100.0).abs() < 1e-3);
}

/// On a 20 bp reference with a budget of floor(0.2 * 20) = 4 edits, two
/// substitutions (90 percent) name it, three (85 percent) still do, and
/// four (80 percent) do not.
#[test]
fn name_against_requires_high_identity() {
    let reference = b"GGGGTTTTGGGGTTTTGGGG";
    let refs = vec![Adapter {
        name: "REF".into(),
        seq: reference.to_vec(),
        role: Role::Adapter,
    }];
    let mutate = |count: usize| -> Vec<u8> {
        let mut seq = reference.to_vec();
        for i in 0..count {
            seq[3 + 5 * i] = b'C';
        }
        seq
    };
    assert_eq!(name_against(&mutate(2), &refs, 0.2).len(), 1);
    assert_eq!(name_against(&mutate(3), &refs, 0.2).len(), 1);
    assert!(
        name_against(&mutate(4), &refs, 0.2).is_empty(),
        "80 percent identity is below the naming floor"
    );
}

/// Reads starting with an exact catalog adapter (`LSK109_front`) yield an
/// adapter named `inferred_1` that carries the catalog match separately.
#[test]
fn discovered_adapters_are_named_by_order_with_catalog_annotation() {
    let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT";
    let mut owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..300usize {
        let mut read = adapter.to_vec();
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
        for _ in 0..120usize {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
        }
        owned.push(read);
    }
    let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let base = AdapterConfig {
        adapters: vec![],
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    };
    let found = discover(&sample, &base);
    assert!(!found.is_empty(), "The planted adapter is discovered");
    for (i, d) in found.iter().enumerate() {
        assert_eq!(d.adapter.name, format!("inferred_{}", i + 1));
    }
    assert_eq!(
        found[0].name_hits.first().map(|(name, _)| name.as_str()),
        Some("LSK109_front"),
        "The catalog match is an annotation: {:?}",
        found[0].name_hits
    );
}

/// Sixty `N`s then random bases: windows holding the run are dropped, so
/// no poly-A-leading consensus is assembled from them.
#[test]
fn discover_finds_nothing_in_ambiguity_runs() {
    let mut owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..200usize {
        let mut read = vec![b'N'; 60];
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
        for _ in 0..100usize {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
        }
        owned.push(read);
    }
    let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let base = AdapterConfig {
        adapters: vec![],
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    };
    let found = discover(&sample, &base);
    assert!(
        found.is_empty(),
        "An N run is not adapter evidence (got {found:?})"
    );
}

/// A catalog-like adapter planted at the 5' end of 500 synthetic reads with
/// about 10 percent substitution error is recovered within a small edit
/// distance. The noise is a fixed permutation of error positions per read
/// index, with no RNG.
#[test]
fn discover_recovers_planted_adapter_under_error() {
    let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp, SQK-NSK007-like
    let mut owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..500usize {
        let mut read = adapter.to_vec();
        // Deterministic genomic tail from a splitmix64-style mix. A formula
        // linear in the position modulo 4 collapses to a phase-rotated ACGT
        // tandem repeat, which is a spurious signal in 100% of reads that
        // crowds out the planted adapter.
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
        for _ in 0..120usize {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
        }
        // Deterministic substitutions at roughly 10% of adapter positions.
        for p in (0..adapter.len()).step_by(10) {
            let q = (p + i) % adapter.len();
            read[q] = b"ACGT"[(read[q] as usize + 1) % 4];
        }
        owned.push(read);
    }
    let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let base = AdapterConfig {
        adapters: vec![],
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    };
    let found = discover(&sample, &base);
    assert!(!found.is_empty(), "At least one adapter discovered");
    // The top candidate is a 5' or both-end adapter close to the planted
    // sequence.
    let top = &found[0];
    assert!(top.adapter.seq.len() >= MIN_PATTERN_LEN);
    // Near-match to the planted adapter; recovery is approximate.
    let mut s = new_ambiguous_searcher();
    let k = (0.25 * adapter.len() as f64).ceil() as usize;
    assert!(
        !hits(&mut s, &top.adapter.seq, adapter, k).is_empty()
            || !hits(&mut s, adapter, &top.adapter.seq, k).is_empty(),
        "Recovered adapter is within 25% edit distance of the planted one"
    );
}

/// `adapter` is planted at the 5' end with heavy substitutions (weak
/// recovery) and its exact reverse complement at the 3' end (strong
/// recovery) of every read, so `merge_both_ends` folds the two per-end
/// discoveries into a single `End::Both` entry (per `same_adapter`). The
/// noisy 5' and exact 3' assemblies are fuzzy-equivalent, so the merged
/// adapter inherits the stronger 3' support.
#[test]
fn discover_dual_end_adapter_gets_max_support() {
    let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp
    let rc: Vec<u8> = adapter
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!("Adapter is pure ACGT"),
        })
        .collect();
    let mut owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..200usize {
        // 5' copy: deterministic substitutions at every 6th (shifted)
        // position; weak but still independently recoverable.
        let mut read = adapter.to_vec();
        for p in (0..adapter.len()).step_by(6) {
            let q = (p + i) % adapter.len();
            read[q] = b"ACGT"[(read[q] as usize + 1) % 4];
        }
        // Deterministic non-periodic genomic middle, from the same
        // splitmix64 mix as the other `discover_*` fixtures.
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
        for _ in 0..150usize {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
        }
        // 3' copy: exact reverse complement with no error, giving strong
        // recovery.
        read.extend_from_slice(&rc);
        owned.push(read);
    }
    let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let base = AdapterConfig {
        adapters: vec![],
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    };
    let found = discover(&sample, &base);

    // The merged entry is the one near the planted adapter; recovery is
    // approximate, so the match is within 25% edit distance.
    let mut s = new_ambiguous_searcher();
    let k = (0.25 * adapter.len() as f64).ceil() as usize;
    let near: Vec<&InferredAdapter> = found
        .iter()
        .filter(|d| {
            !hits(&mut s, &d.adapter.seq, adapter, k).is_empty()
                || !hits(&mut s, adapter, &d.adapter.seq, k).is_empty()
        })
        .collect();
    assert_eq!(
        near.len(),
        1,
        "The shared 5'/3' adapter is discovered as a single entry: {found:?}"
    );

    // The reported support reflects the stronger 3' end, not the weaker 5'
    // end alone (about 0.18). The unmerged `Five` entries this fixture also
    // produces carry that value and are dropped independently because
    // 0.18 < `KEEP_SUPPORT`.
    assert!(
        near[0].support > 0.7,
        "Merged adapter's support ({}) must reflect the max across ends \
             (the 3' end recovers at about 1.0), not the weaker 5' end alone (about 0.18)",
        near[0].support
    );
}

/// Deterministic non-periodic background (SplitMix64-derived bases) yields
/// no adapter.
#[test]
fn discover_finds_nothing_in_clean_reads() {
    let mut owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..300usize {
        let mut read = Vec::new();
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
        for _ in 0..200usize {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
        }
        owned.push(read);
    }
    let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let base = AdapterConfig {
        adapters: vec![],
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    };
    let found = discover(&sample, &base);
    assert!(
        found.is_empty(),
        "No spurious adapter in clean reads (got {found:?})"
    );
}

/// With `len <= cap` every window is returned in order (step 1).
#[test]
fn stride_sample_is_identity_when_within_cap() {
    let a: &[u8] = b"A";
    let b: &[u8] = b"C";
    let c: &[u8] = b"G";
    let windows: Vec<&[u8]> = vec![a, b, c];
    assert_eq!(stride_sample(&windows, 4), windows);
}

/// A four-element sample spans all 13 input positions rather than a prefix.
#[test]
fn stride_sample_spans_the_whole_range_not_just_a_prefix() {
    let bytes: Vec<u8> = (0..13u8).map(|i| b'A' + i).collect();
    let windows: Vec<&[u8]> = bytes.iter().map(std::slice::from_ref).collect();
    let sampled = stride_sample(&windows, 4);
    assert!(sampled.len() <= 4);
    // Expected indices: 0, 4, 8, 12.
    assert_eq!(
        sampled,
        vec![windows[0], windows[4], windows[8], windows[12]]
    );
    let last_idx = 12usize; // index of the last sampled window
    assert!(
        last_idx >= (13usize * 2).div_ceil(3),
        "Last sampled window must fall in the last third of the range, not a prefix"
    );
    assert_eq!(*sampled.last().unwrap(), windows[last_idx]);
}

/// The adapter occurs only in the latter half of an 8001-read sample, so
/// the bounded recount has to cover the complete input range. Ignored by
/// default; `cargo test --lib` runs it only when `--ignored` is passed
/// through to the test binary.
#[test]
#[ignore]
fn discover_is_not_order_biased_by_recount_window_cap() {
    let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp, same as the other discover_* fixtures
    let n_clean = RECOUNT_WINDOWS + 1;
    let n_planted = RECOUNT_WINDOWS; // 4000

    // Deterministic non-periodic background, from the same splitmix64 mix
    // as the other `discover_*` fixtures.
    let splitmix_tail = |i: usize, len: usize| -> Vec<u8> {
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            out.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
        }
        out
    };

    let mut owned: Vec<Vec<u8>> = Vec::with_capacity(n_clean + n_planted);
    for i in 0..n_clean {
        owned.push(splitmix_tail(i, 40)); // pure background, no adapter
    }
    for i in 0..n_planted {
        let mut read = adapter.to_vec();
        read.extend(splitmix_tail(n_clean + i, 12));
        owned.push(read);
    }
    let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let base = AdapterConfig {
        adapters: vec![],
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    };
    let found = discover(&sample, &base);
    assert!(
        !found.is_empty(),
        "Adapter present in a clear majority of reads after the first \
             RECOUNT_WINDOWS must be discovered (got {found:?})"
    );
    let mut s = new_ambiguous_searcher();
    let k = (0.25 * adapter.len() as f64).ceil() as usize;
    assert!(
        found.iter().any(|d| {
            !hits(&mut s, &d.adapter.seq, adapter, k).is_empty()
                || !hits(&mut s, adapter, &d.adapter.seq, k).is_empty()
        }),
        "Discovered adapters must include one within 25% edit distance \
             of the planted adapter: {found:?}"
    );
}

/// Lowercase reads produce the same inferred adapter as uppercase DNA, and
/// the discovered sequence is uppercase.
#[test]
fn discover_recovers_planted_adapter_from_lowercase_reads() {
    let adapter: &[u8] = b"AATGTACTTCGTTCAGTTACGTATTGCT"; // 28bp
    let mut owned: Vec<Vec<u8>> = Vec::new();
    for i in 0..500usize {
        let mut read = adapter.to_vec();
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(i as u64);
        for _ in 0..120usize {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            read.push(b"ACGT"[((z >> 62) & 0b11) as usize]);
        }
        let lower: Vec<u8> = read.iter().map(u8::to_ascii_lowercase).collect();
        owned.push(lower);
    }
    let sample: Vec<&[u8]> = owned.iter().map(|v| v.as_slice()).collect();
    let base = AdapterConfig {
        adapters: vec![],
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 1,
        candidate_index: std::sync::OnceLock::new(),
    };
    let found = discover(&sample, &base);
    assert!(
        !found.is_empty(),
        "Lowercase reads must be inferable (got {found:?})"
    );
    let top = &found[0];
    let mut s = new_ambiguous_searcher();
    let k = (0.25 * adapter.len() as f64).ceil() as usize;
    assert!(
        !hits(&mut s, &top.adapter.seq, adapter, k).is_empty()
            || !hits(&mut s, adapter, &top.adapter.seq, k).is_empty(),
        "Discovered adapter (seq {:?}) must be within 25% edit distance \
             of the uppercase planted adapter",
        String::from_utf8_lossy(&top.adapter.seq)
    );
    // The discovered sequence is uppercase ACGT and carries no lowercase
    // byte through from the input.
    assert!(
        top.adapter.seq.iter().all(u8::is_ascii_uppercase),
        "Discovered sequence must be uppercase: {:?}",
        String::from_utf8_lossy(&top.adapter.seq)
    );
}

fn rc(seq: &[u8]) -> Vec<u8> {
    crate::adapter::reverse_complement(seq)
}

fn contains(hay: &[u8], needle: &[u8]) -> bool {
    needle.len() <= hay.len() && hay.windows(needle.len()).any(|w| w == needle)
}

#[test]
fn known_adapter_precedes_a_discovered_primer_layer() {
    let adapter = random_bases(5011, 30);
    let primer = random_bases(5023, 22);
    let reads: Vec<Vec<u8>> = (0..600)
        .map(|i| {
            let mut read = adapter.clone();
            read.extend_from_slice(&primer);
            read.extend(random_bases(7000 + i, 400));
            read
        })
        .collect();
    let known = vec![Adapter {
        name: "known".into(),
        seq: adapter.clone(),
        role: Role::Adapter,
    }];
    let found = infer_with_known(&reads, known);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0].adapter.seq, primer);
    assert_eq!(found[0].adapter.role, Role::Primer);
    assert_eq!(found[0].layer, 0);
}

#[test]
fn discovers_barcode_layer_between_flanks() {
    let adapter = random_bases(6011, 30);
    let flank1 = random_bases(6023, 16);
    let flank2 = random_bases(6031, 32);
    let barcodes: Vec<Vec<u8>> = (0..8).map(|i| random_bases(6100 + i, 24)).collect();
    let reads: Vec<Vec<u8>> = (0..1600)
        .map(|i| {
            let mut read = adapter.clone();
            read.extend_from_slice(&flank1);
            read.extend_from_slice(&barcodes[i % 8]);
            read.extend_from_slice(&flank2);
            read.extend(random_bases(9000 + i as u64, 400));
            read
        })
        .collect();
    let found = infer_owned(&reads);
    for d in &found {
        let technical = barcodes.iter().any(|b| {
            contains(
                &[&adapter[..], &flank1, b, &flank2].concat(),
                &d.adapter.seq,
            )
        });
        assert!(technical, "{:?}", String::from_utf8_lossy(&d.adapter.seq));
        let expected = if contains(&d.adapter.seq, &adapter[..16]) {
            Role::Adapter
        } else if barcodes.iter().any(|b| contains(&d.adapter.seq, &b[4..20])) {
            Role::Barcode
        } else {
            Role::Primer
        };
        assert_eq!(d.adapter.role, expected, "{d:?}");
    }
    assert!(
        found
            .iter()
            .any(|d| d.layer == 0 && contains(&d.adapter.seq, &adapter[..16]))
    );
    assert!(
        found
            .iter()
            .any(|d| d.layer > 0 && contains(&d.adapter.seq, &flank2)),
        "{found:?}"
    );
    let recovered = barcodes
        .iter()
        .filter(|b| {
            found
                .iter()
                .any(|d| d.adapter.role == Role::Barcode && contains(&d.adapter.seq, &b[4..20]))
        })
        .count();
    assert!(recovered >= 6, "{recovered} barcodes recovered: {found:?}");
    let cfg = AdapterConfig {
        adapters: found.into_iter().map(|d| d.adapter).collect(),
        error_rate: 0.2,
        end_size: 150,
        split: true,
        min_piece: 20,
        candidate_index: std::sync::OnceLock::new(),
    };
    for read in reads.iter().take(16) {
        let segments = crate::adapter::adapter_segments(read, &cfg);
        assert_eq!(segments, vec![(102, read.len())], "{segments:?}");
    }
}

#[test]
fn eroded_short_barcodes_end_at_the_insert_boundary() {
    let barcodes: Vec<Vec<u8>> = (0..12).map(|i| random_bases(8100 + i, 16)).collect();
    let reads: Vec<Vec<u8>> = (0..2400)
        .map(|i| {
            let barcode = &barcodes[i % 12];
            let mut read = barcode[(i / 12) % 2..].to_vec();
            read.extend(random_bases(12000 + i as u64, 400));
            read
        })
        .collect();
    let found = infer_owned(&reads);
    for d in &found {
        assert!(
            barcodes.iter().any(|b| contains(b, &d.adapter.seq)),
            "{:?}",
            String::from_utf8_lossy(&d.adapter.seq)
        );
    }
    let recovered = barcodes
        .iter()
        .filter(|b| found.iter().any(|d| contains(&d.adapter.seq, &b[1..])))
        .count();
    assert!(recovered >= 10, "{recovered} barcodes recovered: {found:?}");
}

#[test]
fn conservation_is_judged_against_the_background_composition() {
    let uniform = [0.25; 4];
    let at_rich = [0.4, 0.1, 0.1, 0.4];
    assert!(conserved([97, 1, 1, 1], uniform));
    assert!(conserved([3, 90, 3, 4], at_rich));
    assert!(
        conserved([1, 48, 1, 50], uniform),
        "two-fold degenerate base"
    );
    assert!(
        conserved([1, 48, 1, 50], at_rich),
        "two-fold degenerate base"
    );
    assert!(!conserved([25, 25, 25, 25], uniform));
    assert!(!conserved([30, 20, 20, 30], uniform));
    assert!(!conserved([40, 10, 10, 40], at_rich));
    assert!(!conserved([47, 8, 5, 40], at_rich));
}

/// Returns `len` bases of which about 78% are A or T, the composition of
/// the most AT-rich sequenced genomes.
fn at_rich_bases(seed: u64, len: usize) -> Vec<u8> {
    let first = random_bases(seed ^ 0x5555, len);
    let second = random_bases(seed ^ 0xAAAA, len);
    random_bases(seed, len)
        .into_iter()
        .zip(first.into_iter().zip(second))
        .map(|(base, draws)| match (base, draws) {
            (b'C' | b'G', (b'A' | b'G' | b'T', b'A' | b'G')) => b'A',
            (b'C' | b'G', (b'A' | b'G' | b'T', b'T')) => b'T',
            _ => base,
        })
        .collect()
}

#[test]
fn eroded_short_barcodes_end_at_an_at_rich_insert_boundary() {
    let barcodes: Vec<Vec<u8>> = (0..12).map(|i| random_bases(8300 + i, 16)).collect();
    let reads: Vec<Vec<u8>> = (0..2400)
        .map(|i| {
            let barcode = &barcodes[i % 12];
            let mut read = barcode[(i / 12) % 2..].to_vec();
            read.extend(at_rich_bases(14000 + i as u64, 400));
            read
        })
        .collect();
    let found = infer_owned(&reads);
    for d in &found {
        assert!(
            barcodes.iter().any(|b| contains(b, &d.adapter.seq)),
            "{:?}",
            String::from_utf8_lossy(&d.adapter.seq)
        );
    }
    let recovered = barcodes
        .iter()
        .filter(|b| found.iter().any(|d| contains(&d.adapter.seq, &b[1..])))
        .count();
    assert!(recovered >= 10, "{recovered} barcodes recovered: {found:?}");
}

#[test]
fn conserved_insert_mirrored_at_truncated_read_ends_is_excluded() {
    let adapter = random_bases(7011, 30);
    let primer = random_bases(7023, 20);
    let conserved = random_bases(7031, 40);
    let reads: Vec<Vec<u8>> = (0..800)
        .map(|i| {
            let mut read = adapter.clone();
            read.extend_from_slice(&primer);
            read.extend_from_slice(&conserved);
            read.extend(random_bases(11000 + i as u64, 300));
            read.extend(rc(&conserved));
            read.extend_from_slice(&rc(&primer)[..7]);
            read
        })
        .collect();
    let found = infer_owned(&reads);
    assert!(!found.is_empty());
    let technical = [&adapter[..], &primer].concat();
    for d in &found {
        assert!(
            contains(&technical, &d.adapter.seq),
            "{:?}",
            String::from_utf8_lossy(&d.adapter.seq)
        );
    }
    assert!(
        found.iter().any(|d| contains(&d.adapter.seq, &adapter)),
        "{found:?}"
    );
}

/// A barcode construct names no inferred sequence: its `N` block would match
/// any sequence at full identity.
#[test]
fn barcode_construct_names_no_inferred_sequence() {
    let construct = Adapter {
        name: "NB_construct".into(),
        seq: b"ATTGCTAAGGTTAANNNNNNNNNNNNNNNNNNNNNNNNCAGCACCT".to_vec(),
        role: Role::Barcode,
    };
    assert!(name_against(b"TGACTCCTCGCTTTCGA", &[construct], 0.15).is_empty());
}

/// A pattern that shares its outer part with a window and pays for a
/// different inner part in edits does not cover the window; the same pattern
/// covers a window holding all of it.
#[test]
fn windows_covered_requires_the_whole_pattern() {
    let shared = b"ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT";
    let mut pattern = shared.to_vec();
    pattern.extend_from_slice(b"TACGGCTACCTTGTTACGACTT");
    let mut other = shared.to_vec();
    other.extend_from_slice(b"AGAGTTTGATCCTGGCTCAGTTACCGATGG");
    let mut own = pattern.clone();
    own.extend_from_slice(b"GGATCCATTAC");
    let mut searcher = crate::adapter::search::new_searcher_fwd();
    let k = crate::adapter::edit_budget(0.2, pattern.len());
    let covered = windows_covered(&mut searcher, &pattern, &[&other, &own], k);
    assert_eq!(covered, [false, true]);
}

/// Returns a splitmix64 draw from `state`, advanced in place.
fn draw(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Copies `seq` with `per_mille` errors per thousand bases, a third each
/// substitutions, deletions and insertions.
fn with_errors(seq: &[u8], seed: u64, per_mille: u64) -> Vec<u8> {
    let mut state = seed;
    let mut out = Vec::with_capacity(seq.len() + 8);
    for &base in seq {
        let roll = draw(&mut state) % 3000;
        let pick = b"ACGT"[(draw(&mut state) % 4) as usize];
        if roll < per_mille {
            out.push(if pick == base {
                b"ACGT"[(pick as usize + 1) % 4]
            } else {
                pick
            });
        } else if roll < 2 * per_mille {
        } else if roll < 3 * per_mille {
            out.extend_from_slice(&[base, pick]);
        } else {
            out.push(base);
        }
    }
    out
}

/// Reads of both strands of a cDNA library with 3% errors: the outer
/// adapter, noisier and eroded by a geometric number of bases, a shared core,
/// then one strand primer at the 5' end and the other, reverse complemented,
/// at the 3' end, where less outer adapter remains. The reverse strand primer
/// is followed by a poly(T) run of varying length.
fn strand_primer_reads() -> (Vec<Vec<u8>>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let outer = random_bases(8101, 45);
    let core = random_bases(8111, 40);
    let forward = random_bases(8121, 25);
    let reverse = random_bases(8131, 30);
    let reads = (0..1500u64)
        .map(|i| {
            let mut state = 30_000 + i;
            let mut erosion = 0;
            while erosion < 30 && !draw(&mut state).is_multiple_of(5) {
                erosion += 1;
            }
            let mut read = with_errors(&outer[erosion..], 40_000 + i, 120);
            let mut rest = core.clone();
            rest.extend_from_slice(&forward);
            rest.extend(random_bases(20_000 + i, 400));
            rest.extend(std::iter::repeat_n(b'A', 12 + (i % 9) as usize));
            rest.extend(rc(&reverse));
            rest.extend(rc(&core));
            read.extend(with_errors(&rest, 50_000 + i, 30));
            read.extend(with_errors(
                &rc(&outer)[..5 + (i % 7) as usize],
                60_000 + i,
                120,
            ));
            if i % 2 == 1 { rc(&read) } else { read }
        })
        .collect();
    (reads, core, forward, reverse)
}

/// A shared adapter core that divides into the primers of the two strands
/// is one layer, and each primer is the next; the poly(T) behind the reverse
/// primer is left to the insert. The core lies deeper at the 5' end than its
/// mirror at the 3' end, which keeps less outer adapter.
#[test]
fn shared_core_divides_into_strand_primers() {
    let (reads, core, forward, reverse) = strand_primer_reads();
    let found = infer_owned(&reads);
    let seqs: Vec<Vec<u8>> = found.iter().map(|d| d.adapter.seq.clone()).collect();
    let overlap = |part: &[u8], seq: &[u8]| {
        longest_common_substring(part, seq).max(longest_common_substring(part, &rc(seq)))
    };
    for (part, least) in [(&core, 30), (&forward, 22), (&reverse, 25)] {
        assert!(
            seqs.iter().any(|s| overlap(part, s) >= least),
            "{:?}",
            String::from_utf8_lossy(part)
        );
    }
    for s in &seqs {
        assert!(
            !contains(s, b"TTTTTTTT") && !contains(s, b"AAAAAAAA"),
            "{:?}",
            String::from_utf8_lossy(s)
        );
    }
}

/// A path divides where two continuations keep their support; a degenerate
/// base opens a bubble that rejoins the path and does not divide it.
#[test]
fn division_needs_continuations_that_do_not_rejoin() {
    let shared = random_bases(8201, 30);
    let left = random_bases(8211, 40);
    let right = random_bases(8221, 40);
    let tail = random_bases(8231, 40);
    let mut windows: Vec<Vec<u8>> = Vec::new();
    for i in 0..200 {
        let mut divided = shared.clone();
        divided.extend_from_slice(if i % 2 == 0 { &left } else { &right });
        windows.push(divided);
        let mut bubble = shared.clone();
        bubble.push(if i % 2 == 0 { b'A' } else { b'C' });
        bubble.extend_from_slice(&tail);
        windows.push(bubble);
    }
    let (divided, bubbled): (Vec<&[u8]>, Vec<&[u8]>) = (
        windows.iter().step_by(2).map(Vec::as_slice).collect(),
        windows
            .iter()
            .skip(1)
            .step_by(2)
            .map(Vec::as_slice)
            .collect(),
    );
    let mut cons = shared.clone();
    cons.extend_from_slice(&left);
    let counts = kmer_counts(&divided);
    let edge = division_point(&cons, &counts, divided.len(), 0, cons.len(), End::Five);
    assert_eq!(edge.map(|(edge, _)| edge), Some(shared.len()));

    let mut cons = shared.clone();
    cons.push(b'A');
    cons.extend_from_slice(&tail);
    let counts = kmer_counts(&bubbled);
    assert!(division_point(&cons, &counts, bubbled.len(), 0, cons.len(), End::Five).is_none());
}

/// The homopolymer run on the insert-facing side of a sequence is measured
/// from the inner end.
#[test]
fn inner_run_measures_the_insert_facing_homopolymer() {
    assert_eq!(inner_run(b"ACGTTTTTT", End::Five), 6);
    assert_eq!(inner_run(b"AAAAACGT", End::Three), 5);
    assert_eq!(inner_run(b"ACGT", End::Five), 1);
}
