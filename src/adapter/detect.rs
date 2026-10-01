//! Adapter presence detection over a sampled read prefix.
//!
//! Narrows a configured adapter set to the entries that would act on at least
//! `presence_min` of the sampled reads, so a large catalog is not searched in
//! full on every read.

use super::{Adapter, AdapterConfig, adapter_segments_tallied};
use rayon::prelude::*;

/// Sample size below which presence detection is unreliable; callers skip it
/// and use the full adapter set.
pub const MIN_SAMPLE_FOR_DETECTION: usize = 100;

/// Returns the minimum sampled-read count for an adapter to be kept: 0.2% of
/// the sample, floored at 3 so that a single stray hit cannot promote an
/// adapter.
pub fn presence_min(sample_size: usize) -> usize {
    (sample_size / 500).max(3)
}

/// Retains the adapters of `cfg` that trim or excise at least `min_count` of
/// the sampled reads, running the same search passes the trimming pass runs.
/// Order is preserved.
pub fn present(
    sample: &[&[u8]],
    cfg: &AdapterConfig,
    min_count: usize,
    threads: usize,
) -> Vec<Adapter> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .expect("A positive Rayon worker count builds a pool");
    let n = cfg.adapters.len();
    let counts = pool.install(|| {
        sample
            .par_iter()
            .fold(
                || vec![0usize; n],
                |mut counts, &read| {
                    let mut acted = vec![false; n];
                    adapter_segments_tallied(read, cfg, &mut acted);
                    for (count, hit) in counts.iter_mut().zip(acted) {
                        *count += usize::from(hit);
                    }
                    counts
                },
            )
            .reduce(
                || vec![0usize; n],
                |mut a, b| {
                    for (x, y) in a.iter_mut().zip(b) {
                        *x += y;
                    }
                    a
                },
            )
    });
    cfg.adapters
        .iter()
        .zip(counts)
        .filter(|(_, count)| *count >= min_count)
        .map(|(adapter, _)| adapter.clone())
        .collect()
}

/// Percentage of the sampled reads that the primers of the set must open for
/// the library to count as an amplicon library. Every molecule of an amplicon
/// library starts with one of its primers, which a read keeps at its 5' end
/// unless the read is a fragment; a genomic library holds a primer site at a
/// read end only where the read happens to start at one.
pub const AMPLICON_PERCENT: usize = 60;

/// Bases beyond the boundary left by the outer layers within which a primer
/// still opens the read end: the anchoring slack of a terminal trim.
const AMPLICON_ANCHOR_SLACK: usize = super::FLANK_SLACK;

/// Most sampled reads, spread evenly over the sample, that `primer_share`
/// examines. The share is a proportion, which this many reads estimate to
/// within about a percentage point.
pub(crate) const AMPLICON_SAMPLE_READS: usize = 4000;

/// Returns whether `seq` can stand for a primer in `primer_share`: it is of
/// searchable length, no UMI pattern (`catalog::UMIS`), holds ambiguity
/// codes at no more than a quarter of its positions, as a degenerate primer
/// does and a random tag does not, and is at least as specific as a plain
/// sequence of `MIN_PATTERN_LEN` bases, so that it does not open read ends
/// by chance.
fn usable_primer(seq: &[u8]) -> bool {
    let umi = super::Adapter {
        name: String::new(),
        seq: seq.to_vec(),
        role: super::Role::Primer,
    };
    let ambiguous = seq
        .iter()
        .filter(|&&b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
        .count();
    seq.len() >= super::MIN_PATTERN_LEN
        && !super::is_umi(&umi)
        && 4 * ambiguous <= seq.len()
        && super::chance_cumulative(seq, 0)[0] <= 0.25f64.powi(super::MIN_PATTERN_LEN as i32)
}

/// Bases by which a read length may differ from another, as a percentage of
/// it, for `length_share` to count them alike.
const LENGTH_WINDOW_PERCENT: usize = 5;

/// Percentage of the sampled reads that the median read must share its
/// length with (`length_share`) for the library to count as an amplicon
/// library. The reads of an amplicon library share the lengths of a few
/// targets; the fragments of a genome or transcriptome spread over a wide
/// range, and each shares its length with a few percent of the others.
pub const AMPLICON_LENGTH_PERCENT: usize = 15;

/// Returns, over `sample`, the median share of the reads whose length lies
/// within `LENGTH_WINDOW_PERCENT` of a read's own length.
pub fn length_share(sample: &[&[u8]]) -> f64 {
    let mut lengths: Vec<usize> = sample.iter().map(|read| read.len()).collect();
    lengths.sort_unstable();
    let n = lengths.len();
    if n == 0 {
        return 0.0;
    }
    let mut shares: Vec<f64> = lengths
        .iter()
        .map(|&length| {
            let slack = length * LENGTH_WINDOW_PERCENT / 100;
            let lo = lengths.partition_point(|&l| l < length.saturating_sub(slack));
            let hi = lengths.partition_point(|&l| l <= length + slack);
            (hi - lo) as f64 / n as f64
        })
        .collect();
    shares.sort_by(f64::total_cmp);
    shares[n / 2]
}

/// Length of the insert words that `primer_share` compares across the read
/// ends a primer opens.
const INSERT_WORD: usize = 8;

/// Offsets from the insert boundary at which insert words are read, so that
/// a boundary placed a few bases off by an edit still yields the word.
const INSERT_WORD_OFFSETS: usize = 5;

/// Percentage of the read ends a primer opens that one insert word must
/// start, and the fewest such ends, for the primer to open amplicons: the
/// start of a limited set of target sequences recurs, while random fragments
/// of a genome or transcriptome start anywhere.
const CONSERVED_INSERT_PERCENT: usize = 3;

/// Fewest read ends that must share the insert word of a primer; see
/// `CONSERVED_INSERT_PERCENT`.
const CONSERVED_INSERT_ENDS: usize = 10;

/// Returns the share of `sample` in which a primer of the library opens a
/// read end, over at most `AMPLICON_SAMPLE_READS` reads spread evenly over
/// it, or `None` when the question does not arise: `cfg` splits no read,
/// holds no marker-gene primer in the primer role (the only entries whose
/// interior split depends on the answer; an amplicon-only preset gives its
/// primers the adapter role), or the sample is smaller than
/// `MIN_SAMPLE_FOR_DETECTION`.
///
/// A read end is opened by one of `primers` (`opening_primers`) when a whole
/// hit of it lies at the boundary that the outer layers leave there. The
/// outer layers are the entries of `cfg` other than `primers` and the marker
/// primers. A marker primer counts by itself. Any other primer counts only
/// when the inserts behind it start alike (`conserved_insert`): in an
/// amplicon library a primer is followed by the start of its targets, while
/// a technical sequence that a genomic or cDNA library carries at its read
/// ends is followed by random fragments. The primers that count are taken
/// together, whichever families they form, so a library of several primer
/// pairs is judged as a whole.
pub fn primer_share(
    sample: &[&[u8]],
    cfg: &AdapterConfig,
    primers: &[Vec<u8>],
    threads: usize,
) -> Option<f64> {
    let guarded = cfg.adapters.iter().any(|a| {
        a.role == super::Role::Primer && super::matches_marker_primer(&a.seq, cfg.error_rate)
    });
    if !cfg.split || !guarded || sample.len() < MIN_SAMPLE_FOR_DETECTION {
        return None;
    }
    let mut outer = cfg.clone();
    outer.split = false;
    outer.set_amplicon(false);
    outer.replace_adapters(
        cfg.adapters
            .iter()
            .filter(|a| {
                let seq = a.seq.to_ascii_uppercase();
                !primers.contains(&seq) && !super::matches_marker_primer(&seq, cfg.error_rate)
            })
            .cloned()
            .collect(),
    );
    let primers: Vec<LibraryPrimer<'_>> = primers
        .iter()
        .filter(|seq| usable_primer(seq))
        .map(|seq| LibraryPrimer {
            seq,
            k: super::Budget::new(seq, cfg.error_rate, cfg.end_size).k_end,
            marker: super::matches_marker_primer(seq, cfg.error_rate),
        })
        .collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .expect("A positive Rayon worker count builds a pool");
    let step = sample.len().div_ceil(AMPLICON_SAMPLE_READS).max(1);
    let examined: Vec<&[u8]> = sample.iter().step_by(step).copied().collect();
    let opened: Vec<[Option<(usize, InsertWords)>; 2]> = pool.install(|| {
        examined
            .par_iter()
            .map_init(
                || (super::search::new_ambiguous_searcher(), Vec::new()),
                |(searcher, buf), &read| opening_primers(read, &outer, &primers, searcher, buf),
            )
            .collect()
    });
    let conserved: Vec<bool> = (0..primers.len())
        .map(|primer| {
            primers[primer].marker
                || conserved_insert(
                    opened
                        .iter()
                        .flatten()
                        .flatten()
                        .filter(|(p, _)| *p == primer)
                        .map(|(_, words)| words),
                )
        })
        .collect();
    let reads = opened
        .iter()
        .filter(|ends| ends.iter().flatten().any(|(primer, _)| conserved[*primer]))
        .count();
    Some(reads as f64 / examined.len() as f64)
}

/// A primer of the library as `primer_share` searches it.
struct LibraryPrimer<'a> {
    /// The primer sequence.
    seq: &'a [u8],
    /// Its terminal edit budget.
    k: usize,
    /// Whether it is a marker-gene primer (`matches_marker_primer`).
    marker: bool,
}

/// The insert words at `INSERT_WORD_OFFSETS` successive offsets from an
/// insert boundary, read into the insert, two bits per base; `None` for a
/// word that holds a base outside ACGT, runs past the read, or is of low
/// complexity (`low_complexity`).
type InsertWords = [Option<u16>; INSERT_WORD_OFFSETS];

/// Returns whether `word` is of low complexity: one base fills three
/// quarters of it or more, as in the poly(A) or poly(T) tail of a cDNA read.
fn low_complexity(word: &[u8]) -> bool {
    b"ACGT"
        .iter()
        .any(|&base| 4 * word.iter().filter(|&&b| b == base).count() >= 3 * word.len())
}

/// Returns the insert words of `insert`, the bases from an insert boundary
/// into the insert.
fn insert_words(insert: &[u8]) -> InsertWords {
    std::array::from_fn(|offset| {
        let word = insert.get(offset..offset + INSERT_WORD)?;
        if low_complexity(word) {
            return None;
        }
        word.iter().try_fold(0u16, |code, &base| {
            let bits = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            Some((code << 2) | bits)
        })
    })
}

/// Returns whether the inserts behind one primer start alike: one insert
/// word occurs among the words of at least `CONSERVED_INSERT_PERCENT` of the
/// read ends in `ends`, and of at least `CONSERVED_INSERT_ENDS` of them.
fn conserved_insert<'a>(ends: impl Iterator<Item = &'a InsertWords>) -> bool {
    let mut counts: std::collections::HashMap<u16, usize> = std::collections::HashMap::new();
    let mut total = 0usize;
    for words in ends {
        total += 1;
        let mut seen: Vec<u16> = words.iter().flatten().copied().collect();
        seen.sort_unstable();
        seen.dedup();
        for word in seen {
            *counts.entry(word).or_default() += 1;
        }
    }
    let top = counts.values().copied().max().unwrap_or(0);
    top >= CONSERVED_INSERT_ENDS && top * 100 >= total * CONSERVED_INSERT_PERCENT
}

/// Returns, for the 5' and the 3' end of `read`, the primer that opens the
/// insert there and the insert words behind it (`opening_at`). The 3' end
/// is read on the reverse complement, inward from its end, like the 5' end.
/// `buf` holds the normalized read (`normalize_into`).
fn opening_primers(
    read: &[u8],
    outer: &AdapterConfig,
    primers: &[LibraryPrimer<'_>],
    searcher: &mut super::search::AmbiguousSearcher,
    buf: &mut Vec<u8>,
) -> [Option<(usize, InsertWords)>; 2] {
    let segments = super::adapter_segments(read, outer);
    let (Some(&(lo, _)), Some(&(_, hi))) = (segments.first(), segments.last()) else {
        return [None, None];
    };
    let (text, _) = super::normalize_into(read, buf);
    let reversed = super::reverse_complement(text);
    let reach = outer.end_size;
    [
        opening_at(text, lo, reach, primers, searcher),
        opening_at(&reversed, read.len() - hi, reach, primers, searcher),
    ]
}

/// Returns the primer that opens the insert of `text` read from its start,
/// behind the boundary `lo` that the outer layers leave, with the insert
/// words behind it. Hits are whole, on either strand, as discovery keeps the
/// form each read end shows, and within the terminal edit budget paired with
/// the primer; they are searched within `reach` bases beyond `lo`. A marker
/// primer opens the read end with a hit that starts within
/// `AMPLICON_ANCHOR_SLACK` of `lo` or before it, as inside an outer layer
/// that discovery assembled together with it. Any other primer opens it with
/// a hit that starts within as many bases of `lo` or before it and ends
/// within as many bases of it or after it, so that no outer layer lies
/// behind it; a primer hit that starts within as many bases of the end of
/// the last one continues the stack, as a primer behind a barcode that
/// discovery also reads as a primer does, and the innermost primer of the
/// stack opens the insert.
fn opening_at(
    text: &[u8],
    lo: usize,
    reach: usize,
    primers: &[LibraryPrimer<'_>],
    searcher: &mut super::search::AmbiguousSearcher,
) -> Option<(usize, InsertWords)> {
    let window = &text[..(lo + reach).min(text.len())];
    let hits: Vec<(usize, usize, usize)> = primers
        .iter()
        .enumerate()
        .flat_map(|(primer, entry)| {
            super::search::hits(searcher, entry.seq, window, entry.k)
                .into_iter()
                .map(move |h| (primer, h.start, h.end))
        })
        .collect();
    let slack = AMPLICON_ANCHOR_SLACK;
    if let Some(&(primer, _, end)) = hits
        .iter()
        .filter(|&&(primer, start, _)| primers[primer].marker && start <= lo + slack)
        .min_by_key(|&&(primer, start, _)| (start, primer))
    {
        return Some((primer, insert_words(&text[end..])));
    }
    let mut current = hits
        .iter()
        .filter(|&&(primer, start, end)| {
            !primers[primer].marker && start <= lo + slack && end + slack >= lo
        })
        .min_by_key(|&&(primer, start, _)| (start, primer))
        .copied()?;
    while let Some(&next) = hits
        .iter()
        .filter(|&&(primer, start, end)| {
            !primers[primer].marker
                && start + slack >= current.2
                && start <= current.2 + slack
                && end > current.2
        })
        .max_by_key(|&&(primer, _, end)| (end, std::cmp::Reverse(primer)))
    {
        current = next;
    }
    Some((current.0, insert_words(&text[current.2..])))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::{Role, adapter_segments};

    /// Builds an adapter from its parts.
    fn ad(name: &str, seq: &[u8], role: Role) -> Adapter {
        Adapter {
            name: name.into(),
            seq: seq.to_vec(),
            role,
        }
    }

    /// Builds a configuration at error rate 0.2 with the given end zone.
    fn cfg(adapters: Vec<Adapter>, end_size: usize, split: bool) -> AdapterConfig {
        AdapterConfig {
            adapters,
            error_rate: 0.2,
            end_size,
            split,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
            split_of: Vec::new(),
            split_opens: Vec::new(),
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

    /// The floor of 3 holds up to 1500 sampled reads; above that 0.2% applies.
    #[test]
    fn presence_min_boundaries() {
        assert_eq!(presence_min(0), 3);
        assert_eq!(presence_min(1000), 3);
        assert_eq!(presence_min(10000), 20);
    }

    /// Over 200 reads that each start with adapter P and never contain adapter
    /// Q, P is kept and Q is dropped.
    #[test]
    fn keeps_present_drops_absent() {
        let p = b"GGGGTTTTGGGGTTTTGGGG";
        let q = b"ACGACGACGACGACGACGAC";
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for _ in 0..200 {
            let mut r = p.to_vec();
            r.extend_from_slice(&[b'A'; 60]);
            reads.push(r);
        }
        let seqs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let c = cfg(
            vec![ad("P", p, Role::Adapter), ad("Q", q, Role::Adapter)],
            150,
            true,
        );
        let kept = present(&seqs, &c, presence_min(seqs.len()), 2);
        let names: Vec<&str> = kept.iter().map(|a| a.name.as_str()).collect();
        assert_eq!(names, vec!["P"], "Present kept, absent dropped");
    }

    /// A terminal hit counts, an absent adapter does not, and an interior hit
    /// counts only when splitting is on.
    #[test]
    fn tally_marks_terminal_and_interior_hits() {
        let a = ad("a", b"GGGGTTTTGGGGTTTTGGGG", Role::Adapter);
        let mut term = a.seq.clone();
        term.extend_from_slice(&[b'A'; 60]);
        let mut acted = [false];
        adapter_segments_tallied(&term, &cfg(vec![a.clone()], 150, false), &mut acted);
        assert!(acted[0], "Terminal hit counts");

        let mut acted = [false];
        adapter_segments_tallied(&[b'A'; 80], &cfg(vec![a.clone()], 150, true), &mut acted);
        assert!(!acted[0], "An absent adapter does not count");

        let mut inter = vec![b'A'; 300];
        inter.splice(150..150, a.seq.iter().copied());
        let mut acted = [false];
        adapter_segments_tallied(&inter, &cfg(vec![a.clone()], 20, true), &mut acted);
        assert!(acted[0], "Interior found when split");
        let mut acted = [false];
        adapter_segments_tallied(&inter, &cfg(vec![a.clone()], 20, false), &mut acted);
        assert!(!acted[0], "Interior ignored when ends-only");
    }

    /// Sixty `N`s then random bases: no adapter is present, and detection agrees
    /// with the trimmer, which rewrites the run before searching.
    #[test]
    fn ambiguity_runs_in_reads_are_not_evidence() {
        let a = ad("a", b"GGGGTTTTGGGGTTTTGGGG", Role::Adapter);
        let reads: Vec<Vec<u8>> = (0..200u64)
            .map(|i| {
                let mut r = vec![b'N'; 60];
                r.extend(splitmix_dna(i, 100));
                r
            })
            .collect();
        let seqs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let c = cfg(vec![a.clone()], 150, true);
        let kept = present(&seqs, &c, presence_min(seqs.len()), 2);
        assert!(
            kept.is_empty(),
            "An N run is not adapter evidence: {kept:?}"
        );
        for r in &reads {
            assert_eq!(adapter_segments(r, &c), vec![(0, r.len())]);
        }

        // The adapter after the run: detection and the trimmer both act.
        let mut planted = vec![b'N'; 60];
        planted.extend_from_slice(&a.seq);
        planted.extend(splitmix_dna(9, 80));
        let mut acted = [false];
        adapter_segments_tallied(&planted, &c, &mut acted);
        assert!(acted[0]);
        assert_ne!(adapter_segments(&planted, &c), vec![(0, planted.len())]);
    }

    /// The share counts the reads that a whole primer opens behind the outer
    /// layers of the set, at either end and on either strand: a marker
    /// primer by itself, and another primer where the inserts behind it start
    /// alike. A read whose primer lies deeper is not opened, and a primer
    /// followed by random inserts, as a technical sequence of a genomic
    /// library is, opens none. The share is not measured where no marker
    /// primer holds the primer role or reads are not split.
    #[test]
    fn primer_share_counts_reads_that_primers_open() {
        use crate::adapter::preset::{Kit, preset};
        use crate::adapter::reverse_complement;
        let forward = b"AGAGTTTGATCCTGGCTCAG".as_slice();
        let adapter = b"CCTGTACTTCGTTCAGTTACGTATTGC".as_slice();
        let custom = b"GGTCAACAAATCATAAAGATATTGG".as_slice();
        let target = splitmix_dna(290, 700);
        let reads: Vec<Vec<u8>> = (0..200u64)
            .map(|i| {
                let insert = splitmix_dna(300 + i, 800);
                match i % 5 {
                    0 => [adapter, forward, &insert].concat(),
                    1 => [&insert[..], &reverse_complement(forward)].concat(),
                    2 => [custom, &target].concat(),
                    3 => [adapter, &insert[..40], forward, &insert[40..]].concat(),
                    _ => [&reverse_complement(forward)[..], &insert].concat(),
                }
            })
            .collect();
        let seqs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let mut primers: Vec<Vec<u8>> = crate::adapter::catalog::MARKER_PRIMERS
            .iter()
            .map(|p| p.to_vec())
            .collect();
        primers.push(custom.to_vec());
        let mixed = cfg(preset(&[Kit::Mab114, Kit::Lsk114]), 150, true);
        let share = primer_share(&seqs, &mixed, &primers, 2).unwrap();
        assert!((share - 0.8).abs() < 1e-9, "{share}");

        let technical: Vec<Vec<u8>> = (0..200u64)
            .map(|i| [custom, &splitmix_dna(600 + i, 800)].concat())
            .collect();
        let seqs: Vec<&[u8]> = technical.iter().map(|r| r.as_slice()).collect();
        assert_eq!(primer_share(&seqs, &mixed, &primers, 2), Some(0.0));

        let amplicon_only = cfg(preset(&[Kit::Mab114]), 150, true);
        assert_eq!(primer_share(&seqs, &amplicon_only, &primers, 2), None);
        let ends_only = cfg(preset(&[Kit::Mab114, Kit::Lsk114]), 150, false);
        assert_eq!(primer_share(&seqs, &ends_only, &primers, 2), None);
    }

    /// Reads of a few shared lengths, as amplicons are, share their length
    /// with most others; reads of spread lengths, as genomic fragments are,
    /// with few.
    #[test]
    fn length_share_separates_amplicons_from_fragments() {
        let amplicons: Vec<Vec<u8>> = (0..300usize)
            .map(|i| {
                vec![
                    b'A';
                    if i % 2 == 0 {
                        1500 + i % 20
                    } else {
                        650 + i % 15
                    }
                ]
            })
            .collect();
        let seqs: Vec<&[u8]> = amplicons.iter().map(|r| r.as_slice()).collect();
        assert!(length_share(&seqs) >= 0.45);
        let fragments: Vec<Vec<u8>> = (0..300usize)
            .map(|i| vec![b'A'; 500 + i * 7919 % 9500])
            .collect();
        let seqs: Vec<&[u8]> = fragments.iter().map(|r| r.as_slice()).collect();
        assert!(length_share(&seqs) * 100.0 < AMPLICON_LENGTH_PERCENT as f64);
        assert_eq!(length_share(&[]), 0.0);
    }
}
