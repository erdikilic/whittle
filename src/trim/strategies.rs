//! Quality-based trimming strategies over raw Phred scores.

use crate::qual::phred_to_prob;

/// A list of half-open `[start, end)` index ranges into a read.
type Segments = Vec<(usize, usize)>;

/// Trims low-quality bases from both ends up to the first base with a Phred
/// score at or above `cutoff`. Inputs are raw Phred scores.
pub fn trim_by_quality(phred: &[u8], cutoff: u8) -> Segments {
    let len = phred.len();
    let mut start = 0;
    while start < len && phred[start] < cutoff {
        start += 1;
    }
    let mut end = len;
    while end > start && phred[end - 1] < cutoff {
        end -= 1;
    }
    if end <= start {
        vec![]
    } else {
        vec![(start, end)]
    }
}

/// Returns the segment maximizing the cumulative difference between the
/// cutoff error probability and each base's error probability (modified Mott).
/// Equal scores favor the longer segment. Bases below `cutoff_q` can be retained.
pub fn best_segment(phred: &[u8], cutoff_q: u8) -> Segments {
    let cutoff = phred_to_prob(cutoff_q);
    let mut best_start = usize::MAX;
    let mut best_end = usize::MAX;
    let mut best_cumulative_error = 0.0;
    let mut best_length = 0usize;

    let mut current_start = 0usize;
    let mut current_cumulative_error = -1.0;
    for (i, &q) in phred.iter().enumerate() {
        let prob_error = cutoff - phred_to_prob(q);
        if current_cumulative_error < 0.0 {
            current_cumulative_error = 0.0;
            current_start = i;
        }
        current_cumulative_error += prob_error;
        if best_cumulative_error < current_cumulative_error
            || (best_cumulative_error == current_cumulative_error
                && best_length < i - current_start + 1)
        {
            best_start = current_start;
            best_end = i;
            best_cumulative_error = current_cumulative_error;
            best_length = i - current_start + 1;
        }
    }
    if best_start == usize::MAX {
        vec![]
    } else {
        vec![(best_start, best_end + 1)]
    }
}

/// Score floor of a maximal segment, in error-free bases: a segment is kept
/// when its modified Mott score is at least this many times the cutoff error
/// probability, the score an error-free base contributes. The floor depends on
/// the qualities relative to the cutoff and not on the cutoff itself, so it
/// holds the same number of bases at any threshold.
pub const MAXIMAL_SEGMENT_MIN_BASES: u32 = 50;

/// Returns every maximal scoring segment under the modified Mott score of
/// `best_segment`, in read order, keeping those whose score reaches
/// `MAXIMAL_SEGMENT_MIN_BASES` error-free bases. The segments are the maximal
/// scoring subsequences of Ruzzo and Tompa (1999): the best segment of the
/// read, then recursively the best segments of the parts on either side of it.
/// Two high-quality regions stay separate when the bases between them cost
/// more than the smaller of their two scores.
pub fn maximal_segments(phred: &[u8], cutoff_q: u8) -> Segments {
    let cutoff = phred_to_prob(cutoff_q);
    let min_score = f64::from(MAXIMAL_SEGMENT_MIN_BASES) * cutoff;
    maximal_scoring_segments(phred, |&q| cutoff - phred_to_prob(q), min_score)
}

/// A candidate segment of the Ruzzo-Tompa list: the half-open base range, the
/// cumulative score before its first base and after its last, and the index of
/// the nearest earlier candidate whose starting cumulative score is at most
/// this one's.
struct Candidate {
    start: usize,
    end: usize,
    left: f64,
    right: f64,
    prev: Option<usize>,
}

/// Returns the maximal scoring subsequences of the per-element `score` of
/// `values` whose score is at least `min_score`, in linear time (Ruzzo and
/// Tompa 1999).
///
/// Each run of consecutive non-negative scores enters the candidate list as
/// one candidate: the candidates of its single bases would merge into one
/// another, and a longer run only reaches a higher cumulative score, so the
/// run makes every merge its bases would make. Merging on equal cumulative
/// scores resolves ties toward the longer segment, as `best_segment` does, and
/// includes bases at the cutoff, whose score is zero, at segment edges.
fn maximal_scoring_segments<T>(
    values: &[T],
    score: impl Fn(&T) -> f64,
    min_score: f64,
) -> Segments {
    let mut list: Vec<Candidate> = Vec::new();
    let mut cumulative = 0.0;
    let mut i = 0;
    while i < values.len() {
        let s = score(&values[i]);
        if s < 0.0 {
            cumulative += s;
            i += 1;
            continue;
        }
        let (start, left) = (i, cumulative);
        cumulative += s;
        i += 1;
        while let Some(v) = values.get(i) {
            let s = score(v);
            if s < 0.0 {
                break;
            }
            cumulative += s;
            i += 1;
        }
        insert_candidate(&mut list, start, i, left, cumulative);
    }
    list.into_iter()
        .filter(|c| c.right - c.left >= min_score)
        .map(|c| (c.start, c.end))
        .collect()
}

/// Appends the candidate `[start, end)` spanning cumulative scores `left` to
/// `right`, after it absorbs the earlier candidates the Ruzzo-Tompa rule
/// merges into it: while the rightmost earlier candidate that starts at a
/// cumulative score no higher than the new one ends at one no higher, the new
/// candidate extends back over it and every candidate after it. The `prev`
/// links skip candidates that start higher, so the list is walked in amortized
/// constant time per candidate. Kept out of line, so the per-base loop of
/// `maximal_scoring_segments` holds its running sum in registers.
#[inline(never)]
fn insert_candidate(list: &mut Vec<Candidate>, start: usize, end: usize, left: f64, right: f64) {
    let mut new = Candidate {
        start,
        end,
        left,
        right,
        prev: None,
    };
    let mut j = list.len().checked_sub(1);
    loop {
        while let Some(idx) = j
            && list[idx].left > new.left
        {
            j = list[idx].prev;
        }
        match j {
            Some(idx) if list[idx].right <= new.right => {
                new.start = list[idx].start;
                new.left = list[idx].left;
                j = list[idx].prev;
                list.truncate(idx);
            },
            _ => {
                new.prev = j;
                list.push(new);
                return;
            },
        }
    }
}

/// Splits the read into high-quality segments separated by runs of at least
/// `window` bases below `cutoff`. Every high-quality run is emitted regardless
/// of length; the caller's length filter (`-l`) drops short pieces after
/// trimming.
pub fn split_low_quality(phred: &[u8], cutoff: u8, window: usize) -> Segments {
    let window = window.max(1);
    let mut segments = Vec::new();
    let mut segment_start: Option<usize> = None;
    let mut last_good: Option<usize> = None;
    let mut bad_run = 0usize;

    let push = |start: usize, end: usize, out: &mut Segments| {
        if end > start {
            out.push((start, end));
        }
    };

    for (i, &q) in phred.iter().enumerate() {
        if q >= cutoff {
            if segment_start.is_none() {
                segment_start = Some(i);
            }
            last_good = Some(i);
            bad_run = 0;
        } else {
            bad_run += 1;
            if bad_run >= window {
                if let (Some(s), Some(lg)) = (segment_start, last_good) {
                    push(s, lg + 1, &mut segments);
                }
                segment_start = None;
                last_good = None;
            }
        }
    }
    if let (Some(s), Some(lg)) = (segment_start, last_good) {
        push(s, lg + 1, &mut segments);
    }
    segments
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Shared `(seq, phred)` test vectors for the strategy tests.
    fn reads() -> [(Vec<u8>, Vec<u8>); 6] {
        let raw = [
            (
                b"AAAAAAAAAAAAAAATTTAA".to_vec(),
                b"&#3-G27C:(@G7B55+C4I".to_vec(),
            ),
            (
                b"TTTTTTTTTTTTTTTTTTTT".to_vec(),
                b"77%'24)FAF9@=94'%054".to_vec(),
            ),
            (
                b"AAAAAAAAAAAAAAATTTTA".to_vec(),
                b"'8$-BF2!C;+59->H@91#".to_vec(),
            ),
            (
                b"AAAAAAAAAAAAAAAAAAAA".to_vec(),
                b"%,42$CH*#0+0C6=0,*6/".to_vec(),
            ),
            (
                b"AAAAAAAAAAAAAAAAAAAT".to_vec(),
                b"-------------------J".to_vec(),
            ),
            (
                b"TAAAAAAAAAAAAAAAAAAA".to_vec(),
                b"I-------------------".to_vec(),
            ),
        ];
        raw.map(|(s, q)| (s, q.iter().map(|&b| b - 33).collect()))
    }

    #[test]
    fn trim_by_quality_expected_ranges() {
        let expected: [(u8, Segments); 6] = [
            (20, vec![(4, 20)]),
            (7, vec![(0, 20)]),
            (15, vec![(1, 19)]),
            (40, vec![]),
            (40, vec![(19, 20)]),
            (40, vec![(0, 1)]),
        ];
        for ((cutoff, want), (_, phred)) in expected.iter().zip(reads()) {
            assert_eq!(trim_by_quality(&phred, *cutoff), *want);
        }
    }

    #[test]
    fn best_segment_expected_ranges() {
        // The Q-score cutoffs correspond to probability thresholds: 0.01 = Q20,
        // 0.199 = Q7, 0.0316 = Q15, 0.0001 = Q40.
        let expected: [(u8, Segments); 6] = [
            (20, vec![(10, 16)]),
            (7, vec![(0, 20)]),
            (15, vec![(11, 19)]),
            (40, vec![]),
            (40, vec![(19, 20)]),
            (40, vec![(0, 1)]),
        ];
        for ((cutoff_q, want), (_, phred)) in expected.iter().zip(reads()) {
            assert_eq!(best_segment(&phred, *cutoff_q), *want);
        }
    }

    #[test]
    fn split_expected_segments() {
        // (cutoff, expected) with window 1. Every high-quality run is emitted,
        // including short ones; length filtering belongs to the caller.
        let cases: [(u8, Segments); 6] = [
            (20, vec![(4, 5), (6, 9), (10, 16), (17, 18), (19, 20)]),
            (7, vec![(0, 2), (4, 15), (17, 20)]),
            (15, vec![(1, 2), (4, 7), (8, 10), (11, 13), (14, 19)]),
            (40, vec![]),
            (40, vec![(19, 20)]),
            (40, vec![(0, 1)]),
        ];
        for ((cutoff, want), (_, phred)) in cases.iter().zip(reads()) {
            assert_eq!(split_low_quality(&phred, *cutoff, 1), *want);
        }
    }

    /// Raw Phred scores from a Phred+33 string.
    fn phred(q: &[u8]) -> Vec<u8> {
        q.iter().map(|&b| b - 33).collect()
    }

    /// The recursive definition of the maximal scoring subsequences: the
    /// highest-scoring non-negative segment of the range, the longest and then
    /// the leftmost among equal scores, followed by those of the ranges on
    /// either side. Cubic time; a reference for small inputs.
    fn maximal_reference(scores: &[f64], lo: usize, hi: usize, out: &mut Segments) {
        let mut best: Option<(f64, usize, usize)> = None;
        for s in lo..hi {
            let mut sum = 0.0;
            for e in s + 1..=hi {
                sum += scores[e - 1];
                let better = match best {
                    None => sum >= 0.0,
                    Some((b, bs, be)) => sum > b || (sum == b && e - s > be - bs),
                };
                if better {
                    best = Some((sum, s, e));
                }
            }
        }
        if let Some((_, s, e)) = best {
            maximal_reference(scores, lo, s, out);
            out.push((s, e));
            maximal_reference(scores, e, hi, out);
        }
    }

    /// The maximal scoring subsequences of `scores` with no score floor.
    fn scored_segments(scores: &[f64]) -> Segments {
        maximal_scoring_segments(scores, |&s| s, 0.0)
    }

    /// A deterministic linear congruential generator for the randomized tests.
    fn lcg(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        *state >> 33
    }

    /// Integer scores make every tie exact, so the linear algorithm and the
    /// recursive definition agree on each tie as well as on the segments.
    #[test]
    fn maximal_scoring_segments_match_the_recursive_definition() {
        let mut state = 7;
        for _ in 0..20_000 {
            let len = (lcg(&mut state) % 14) as usize;
            let scores: Vec<f64> = (0..len)
                .map(|_| (lcg(&mut state) % 7) as f64 - 3.0)
                .collect();
            let mut want = Vec::new();
            maximal_reference(&scores, 0, len, &mut want);
            let got = scored_segments(&scores);
            assert_eq!(got, want, "{scores:?}");
        }
    }

    /// A read with one region above the cutoff, scoring above the floor,
    /// yields the segment `best_segment` selects, including the interior dips
    /// it tolerates.
    #[test]
    fn maximal_segments_match_best_segment_on_one_positive_region() {
        let read = [
            phred(b"%%%%&&"),
            vec![40; 60],
            phred(b"+#5"),
            vec![30; 60],
            phred(b"%%%%"),
        ]
        .concat();
        for cutoff in [7, 10, 12, 15] {
            assert_eq!(
                maximal_segments(&read, cutoff),
                best_segment(&read, cutoff),
                "Q{cutoff}"
            );
        }
    }

    /// The best segment is one of the maximal segments of any read, and no
    /// maximal segment scores higher.
    #[test]
    fn best_segment_is_one_of_the_maximal_segments() {
        let mut state = 11;
        for _ in 0..5_000 {
            let len = 1 + (lcg(&mut state) % 200) as usize;
            let q: Vec<u8> = (0..len).map(|_| (lcg(&mut state) % 42) as u8).collect();
            let cutoff = 5 + (lcg(&mut state) % 20) as u8;
            let p = phred_to_prob(cutoff);
            let all = maximal_scoring_segments(&q, |&b| p - phred_to_prob(b), 0.0);
            let best = best_segment(&q, cutoff);
            assert!(best.iter().all(|seg| all.contains(seg)), "{q:?} Q{cutoff}");
        }
    }

    /// Two high-quality flanks around a low-quality interior are both kept,
    /// where `best_segment` keeps only the better flank.
    #[test]
    fn maximal_segments_keep_both_flanks_of_a_low_quality_interior() {
        let q = [vec![30u8; 60], vec![3u8; 40], vec![25u8; 80]].concat();
        assert_eq!(maximal_segments(&q, 10), vec![(0, 60), (100, 180)]);
        assert_eq!(best_segment(&q, 10), vec![(100, 180)]);
    }

    /// An interior region of alternating low and moderate qualities, with no
    /// two consecutive bases below the cutoff, scores negative overall and
    /// separates the flanks; a split at runs of low-quality bases keeps it.
    #[test]
    fn maximal_segments_split_at_a_mixed_low_quality_region() {
        let mixed: Vec<u8> = (0..81).map(|i| if i % 2 == 0 { 4 } else { 11 }).collect();
        let q = [vec![30u8; 60], mixed, vec![30u8; 60]].concat();
        assert_eq!(maximal_segments(&q, 10), vec![(0, 60), (141, 201)]);
        assert_eq!(split_low_quality(&q, 10, 2), vec![(0, 201)]);
    }

    #[test]
    fn maximal_segments_are_empty_for_an_all_low_quality_read() {
        assert!(maximal_segments(&[5u8; 100], 10).is_empty());
        assert!(maximal_segments(&[9u8; 100], 10).is_empty());
    }

    /// A low-quality gap that costs exactly what the flank before it gains is
    /// absorbed, as `best_segment` resolves equal scores toward the longer
    /// segment; bases at the cutoff score zero and extend the segment edges.
    #[test]
    fn maximal_segments_resolve_ties_toward_the_longer_segment() {
        assert_eq!(scored_segments(&[1.0, -1.0, 1.0]), vec![(0, 3)]);
        assert_eq!(scored_segments(&[1.0, -2.0, 1.0]), vec![(0, 1), (2, 3)]);
        let q = [
            vec![5u8; 10],
            vec![10u8; 3],
            vec![40u8; 60],
            vec![10u8; 2],
            vec![5u8; 10],
        ]
        .concat();
        assert_eq!(maximal_segments(&q, 10), vec![(10, 75)]);
        assert_eq!(best_segment(&q, 10), vec![(10, 75)]);
    }

    /// A segment scoring below the floor is dropped and its neighbors are
    /// kept unchanged.
    #[test]
    fn maximal_segments_below_the_score_floor_are_dropped() {
        let floor = f64::from(MAXIMAL_SEGMENT_MIN_BASES) * phred_to_prob(10);
        let bases = (floor / (phred_to_prob(10) - phred_to_prob(20))).ceil() as usize;
        let (short, long) = (bases - 1, bases + 1);
        let gap = vec![2u8; 30];
        let q = [
            vec![20u8; 100],
            gap.clone(),
            vec![20u8; short],
            gap.clone(),
            vec![20u8; long],
            gap,
        ]
        .concat();
        let long_start = 100 + 30 + short + 30;
        assert_eq!(
            maximal_segments(&q, 10),
            vec![(0, 100), (long_start, long_start + long)]
        );
    }

    /// Qualities at the FASTQ cap and at the BAM byte limit score as
    /// near-error-free bases; at a cutoff equal to the cap, bases at the cap
    /// score zero and form no segment above the floor.
    #[test]
    fn maximal_segments_at_the_quality_cap() {
        let q = [vec![93u8; 60], vec![0u8; 60], vec![255u8; 60]].concat();
        assert_eq!(maximal_segments(&q, 20), vec![(0, 60), (120, 180)]);
        assert_eq!(maximal_segments(&[93u8; 500], 93), Vec::new());
    }

    #[test]
    fn maximal_segments_of_empty_and_one_base_reads() {
        assert!(maximal_segments(&[], 10).is_empty());
        assert!(maximal_segments(&[40], 10).is_empty());
        assert!(maximal_segments(&[2], 10).is_empty());
        assert_eq!(scored_segments(&[0.5]), vec![(0, 1)]);
        assert!(scored_segments(&[-0.5]).is_empty());
        assert!(scored_segments(&[]).is_empty());
    }

    #[test]
    fn split_window_tolerates_short_dips() {
        // `III#IIII###III` with `I` = Q40 and `#` = Q2.
        let phred: Vec<u8> = b"III#IIII###III".iter().map(|&b| b - 33).collect();
        assert_eq!(
            split_low_quality(&phred, 10, 1),
            vec![(0, 3), (4, 8), (11, 14)]
        );
        assert_eq!(split_low_quality(&phred, 10, 3), vec![(0, 8), (11, 14)]);
        assert_eq!(split_low_quality(&phred, 10, 4), vec![(0, 14)]);
    }
}
