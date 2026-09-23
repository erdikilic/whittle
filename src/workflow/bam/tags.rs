//! Aux tag rules for trimmed records: the tag classes, per-base array and
//! run-length coverage slicing, and the tag updates of one output window.

use super::*;

/// Per-base `B` arrays, one value per SEQ base, sliced in lockstep with the
/// sequence when a read is trimmed: the PacBio kinetics `ip`/`pw` (single
/// strand) and `fi`/`fp`/`ri`/`rp` (CCS forward/reverse codec-V1), and the
/// PacBio `sm`/`sx` per-base aligned match and mismatch counts. Dorado's scalar
/// `sm:f` (scaling midpoint) shares the `sm` name; the array match leaves it
/// untouched. Any other `B` array whose length equals the read length is also
/// treated as per-base (structural rule), so custom tags need no dedicated
/// handling.
pub(crate) const KNOWN_PERBASE_TAGS: [[u8; 2]; 8] = [
    *b"ip", *b"pw", *b"fi", *b"fp", *b"ri", *b"rp", *b"sm", *b"sx",
];

/// PacBio `sa:B:I`: run-length encoded per-base coverage by subread alignments,
/// stored as `<length>,<coverage>` pairs. Decoded to per-base coverage, sliced
/// to the window and re-encoded (`slice_rle_coverage`). A tag whose run lengths
/// do not sum to the read length is left unchanged and surfaced via
/// `has_malformed_perbase_tag`.
pub(crate) const RLE_COVERAGE_TAG: [u8; 2] = *b"sa";

/// PacBio undo blobs: `ds` (segmented-read reconstitution data for `skera
/// undo`) and `ls` (clipped data for `lima-undo`). Both describe the untrimmed
/// read and are removed from every output record of a trimmed read, counted
/// once per read in `Counters::undo_tags_dropped_reads`.
pub(crate) const UNDO_TAGS: [[u8; 2]; 2] = [*b"ds", *b"ls"];

/// PacBio reverse-strand kinetics: the PacBio BAM spec stores them from the
/// last base to the first, so the window `[start, end)` maps to array indexes
/// `[len - end, len - start)`.
pub(crate) const REVERSED_PERBASE_TAGS: [[u8; 2]; 2] = [*b"ri", *b"rp"];

/// PacBio `B` arrays with a fixed element count unrelated to the read length:
/// `sn` (SNR per channel, 4), `ac` (adapter counts, 4), `bc` (barcode
/// indexes, 2). Excluded from the structural per-base rule, which would
/// otherwise slice them on a read whose length equals their size.
pub(crate) const FIXED_ARRAY_TAGS: [[u8; 2]; 3] = [*b"sn", *b"ac", *b"bc"];

/// Poly-A tail tags handled together with the move table: `pa` (signal
/// boundaries) and `pt` (tail length in bases). `pa` positions are absolute
/// POD5 sample indexes, the frame `ts` uses: dorado adds `num_trimmed_samples`
/// to the anchor and to both boundary ranges before writing the tag
/// (`PolyACalculatorNode.cpp`, `poly_tail_calculator.cpp`). Under
/// `--update-moves` they are kept or shifted when the poly-A tail survives the
/// trim and dropped when it is cut; without it (or with a malformed move table)
/// they are dropped, since signal cannot be related to sequence.
pub(crate) const POLYA_TAGS: [[u8; 2]; 2] = [*b"pa", *b"pt"];

/// Tags dropped on any trimmed read: `bi` (barcode info) embeds front and rear
/// sequence positions that index the untrimmed read and cannot be reconstructed
/// from the BAM, and the PacBio undo blobs `ds`/`ls` describe the untrimmed
/// read. The barcode stage reads `bi` to place the trim (`barcode_window`)
/// before it is dropped. The barcode call itself (`BC`/`bv`) is a per-read label
/// and is copied unchanged.
pub(crate) const DROP_ON_TRIM_TAGS: [[u8; 2]; 3] = [*b"bi", *b"ds", *b"ls"];

/// Tags that describe the whole parent read, not a split subread: `st` (read
/// start time) and `du` (duration). A head/tail crop keeps the same read
/// identity, so they stay valid there. On a split they are recomputed from the
/// sample rate when `--update-moves` resolves the subread's signal window
/// (`split_time_updates`) and dropped otherwise. A non-float `du` (pbmarkdup's
/// `du:Z`) is not a duration and is copied.
pub(crate) const DROP_ON_SPLIT_TAGS: [[u8; 2]; 2] = [*b"st", *b"du"];

/// Returns the element count of a `B` array of any subtype.
pub(crate) fn array_len(a: &Array) -> usize {
    match a {
        Array::Int8(v) => v.len(),
        Array::UInt8(v) => v.len(),
        Array::Int16(v) => v.len(),
        Array::UInt16(v) => v.len(),
        Array::Int32(v) => v.len(),
        Array::UInt32(v) => v.len(),
        Array::Float(v) => v.len(),
    }
}

/// Slices a `B` array of any subtype to `[start, end)` (the element index is the
/// base index for a per-base tag). Subtype-agnostic, so `B:C` (codec-V1) and
/// `B:S` (raw frames) kinetics both work.
pub(super) fn slice_array(a: &Array, start: usize, end: usize) -> Array {
    match a {
        Array::Int8(v) => Array::Int8(v[start..end].to_vec()),
        Array::UInt8(v) => Array::UInt8(v[start..end].to_vec()),
        Array::Int16(v) => Array::Int16(v[start..end].to_vec()),
        Array::UInt16(v) => Array::UInt16(v[start..end].to_vec()),
        Array::Int32(v) => Array::Int32(v[start..end].to_vec()),
        Array::UInt32(v) => Array::UInt32(v[start..end].to_vec()),
        Array::Float(v) => Array::Float(v[start..end].to_vec()),
    }
}

/// Returns a per-base `B` array (any array whose length equals the read length,
/// which covers the known kinetics tags and any custom per-base tag) sliced to
/// the window `[start, end)`, or `None` to leave the tag unchanged.
/// Reverse-strand kinetics are stored last base first and are sliced from the
/// other end. Callers must already have excluded MM/ML/MN and the signal tags.
/// A known kinetics tag whose length does not match is left unchanged and
/// surfaced via `has_malformed_perbase_tag`.
pub(super) fn perbase_slice(
    tag: [u8; 2],
    value: &Value,
    orig_len: usize,
    start: usize,
    end: usize,
) -> Option<Value> {
    if FIXED_ARRAY_TAGS.contains(&tag) {
        return None;
    }
    match value {
        Value::Array(arr) if array_len(arr) == orig_len => {
            let (s, e) = if REVERSED_PERBASE_TAGS.contains(&tag) {
                (orig_len - end, orig_len - start)
            } else {
                (start, end)
            };
            Some(Value::Array(slice_array(arr, s, e)))
        },
        _ => None,
    }
}

/// Returns the integers a `B` array holds, `None` for a float array.
pub(super) fn array_integers(a: &Array) -> Option<Vec<i64>> {
    Some(match a {
        Array::Int8(v) => v.iter().map(|&n| i64::from(n)).collect(),
        Array::UInt8(v) => v.iter().map(|&n| i64::from(n)).collect(),
        Array::Int16(v) => v.iter().map(|&n| i64::from(n)).collect(),
        Array::UInt16(v) => v.iter().map(|&n| i64::from(n)).collect(),
        Array::Int32(v) => v.iter().map(|&n| i64::from(n)).collect(),
        Array::UInt32(v) => v.iter().map(|&n| i64::from(n)).collect(),
        Array::Float(_) => return None,
    })
}

/// Returns `values` stored at the subtype of `template`, `None` when a value
/// does not fit that subtype or the template is a float array.
pub(super) fn array_at_subtype(template: &Array, values: &[i64]) -> Option<Array> {
    fn narrow<T: TryFrom<i64>>(values: &[i64]) -> Option<Vec<T>> {
        values.iter().map(|&n| T::try_from(n).ok()).collect()
    }
    Some(match template {
        Array::Int8(_) => Array::Int8(narrow(values)?),
        Array::UInt8(_) => Array::UInt8(narrow(values)?),
        Array::Int16(_) => Array::Int16(narrow(values)?),
        Array::UInt16(_) => Array::UInt16(narrow(values)?),
        Array::Int32(_) => Array::Int32(narrow(values)?),
        Array::UInt32(_) => Array::UInt32(narrow(values)?),
        Array::Float(_) => return None,
    })
}

/// Returns the number of bases the `<length>,<coverage>` runs of an `sa` array
/// cover. `None` for an odd element count or a negative or overflowing length.
pub(super) fn rle_runs_len(runs: &[i64]) -> Option<usize> {
    if !runs.len().is_multiple_of(2) {
        return None;
    }
    runs.iter().step_by(2).try_fold(0usize, |sum, &len| {
        sum.checked_add(usize::try_from(len).ok()?)
    })
}

/// Returns the number of bases an `sa` run-length coverage array covers,
/// `None` when it is not a well-formed run list (`rle_runs_len`).
pub(super) fn rle_coverage_len(a: &Array) -> Option<usize> {
    rle_runs_len(&array_integers(a)?)
}

/// Slices an `sa` run-length coverage array to the window `[start, end)`: each
/// run is clipped to the window and adjacent runs of equal coverage are merged.
/// `None` when the runs do not cover exactly `orig_len` bases, which leaves the
/// tag unchanged.
pub(super) fn slice_rle_coverage(
    a: &Array,
    orig_len: usize,
    start: usize,
    end: usize,
) -> Option<Array> {
    let runs = array_integers(a)?;
    if rle_runs_len(&runs)? != orig_len {
        return None;
    }
    let mut out: Vec<i64> = Vec::new();
    let mut pos = 0usize;
    for &[len, coverage] in runs.as_chunks::<2>().0 {
        // `rle_runs_len` has checked that every length is a non-negative
        // `usize` and that the sum fits.
        let len = len as usize;
        let run_start = pos;
        pos += len;
        let kept = run_start.max(start)..pos.min(end);
        if kept.is_empty() {
            continue;
        }
        match out.as_mut_slice() {
            [.., last_len, last_coverage] if *last_coverage == coverage => {
                *last_len += kept.len() as i64;
            },
            _ => out.extend([kept.len() as i64, coverage]),
        }
    }
    array_at_subtype(a, &out)
}

/// Returns a tag's value for the window `[start, end)`: the `sa` coverage runs
/// re-encoded, or a per-base array sliced (`perbase_slice`). `None` leaves the
/// tag unchanged.
pub(super) fn windowed_value(
    tag: [u8; 2],
    value: &Value,
    orig_len: usize,
    start: usize,
    end: usize,
) -> Option<Value> {
    if tag == RLE_COVERAGE_TAG {
        return match value {
            Value::Array(a) => slice_rle_coverage(a, orig_len, start, end).map(Value::Array),
            _ => None,
        };
    }
    perbase_slice(tag, value, orig_len, start, end)
}

/// Returns the integer an aux value holds, whatever width it was stored at.
///
/// SAM integer tags are written at the smallest subtype that fits, so a tag is
/// `C` below 256, `S` below 65536 and `I` above that. Matching on one subtype
/// therefore fails on most real records.
pub(super) fn aux_integer(value: &Value) -> Option<i64> {
    Some(match value {
        Value::UInt8(n) => i64::from(*n),
        Value::Int8(n) => i64::from(*n),
        Value::UInt16(n) => i64::from(*n),
        Value::Int16(n) => i64::from(*n),
        Value::UInt32(n) => i64::from(*n),
        Value::Int32(n) => i64::from(*n),
        _ => return None,
    })
}

/// Aux tag updates for one output window: `(tag, Some(value))` replaces the
/// tag in place or appends it when the source lacks it, `(tag, None)` removes
/// it.
pub(super) type TagUpdates = Vec<(Tag, Option<Value>)>;

/// Returns true if the record carries a known per-base tag whose array length
/// disagrees with the sequence length, or an `sa` coverage array whose runs do
/// not sum to it: a malformed per-base tag that cannot be sliced. Used only to
/// emit a run-level advisory.
pub(crate) fn has_malformed_perbase_tag(rec: &RecordBuf, seq_len: usize) -> bool {
    rec.data().iter().any(|(tag, value)| {
        let t = <[u8; 2]>::from(tag);
        match value {
            Value::Array(a) if KNOWN_PERBASE_TAGS.contains(&t) => array_len(a) != seq_len,
            Value::Array(a) if t == RLE_COVERAGE_TAG => rle_coverage_len(a) != Some(seq_len),
            _ => false,
        }
    })
}

/// Returns whether `t` is handled by a dedicated rule rather than the structural
/// per-base slice: the modification block, the signal and poly-A tags, and the
/// tags dropped on a trim or a split.
pub(super) fn has_dedicated_rule(t: [u8; 2]) -> bool {
    matches!(&t, b"MM" | b"ML" | b"MN")
        || SIGNAL_TAGS.contains(&t)
        || POLYA_TAGS.contains(&t)
        || DROP_ON_TRIM_TAGS.contains(&t)
        || DROP_ON_SPLIT_TAGS.contains(&t)
}

/// Replaces the update for `tag` or appends one.
pub(super) fn set_update(updates: &mut TagUpdates, tag: Tag, value: Option<Value>) {
    match updates.iter_mut().find(|(t, _)| *t == tag) {
        Some(slot) => slot.1 = value,
        None => updates.push((tag, value)),
    }
}

/// Computes the aux tag updates for output window `window` of `src`, shared by
/// the BAM and FASTQ output paths: `(tag, Some(value))` replaces the tag or
/// adds it, `(tag, None)` removes it. Empty for a full, unsplit window. Covers
/// the signal and poly-A tags, the tags dropped on a trim, the `qs:f` refresh,
/// the PacBio `qs`/`qe` coordinates, `st`/`du` on a split, and the dorado
/// subread tags `rn`, `pi`, `me` and `er` on an ONT split. The modification
/// block and the per-base arrays are the callers' own.
pub(super) fn window_tag_updates(
    src: &RecordBuf,
    qual: &[u8],
    window: Window,
    platform: Platform,
    moves: Option<&MoveIndex<'_>>,
) -> TagUpdates {
    let Window {
        start,
        end,
        idx,
        total,
    } = window;
    let orig_len = qual.len();
    let trimmed = start != 0 || end != orig_len;
    let split = total > 1;
    if !trimmed && !split {
        return Vec::new();
    }

    let (mut updates, signal) = signal_tag_updates(src, orig_len, start, end, total, moves);
    if trimmed {
        updates.extend(DROP_ON_TRIM_TAGS.map(|t| (Tag::new(t[0], t[1]), None)));
        let qs = Tag::new(b'q', b's');
        // Dorado's `qs:f` is the mean read qscore and follows the trimmed
        // quality.
        if matches!(src.data().get(&qs), Some(Value::Float(_))) {
            let mean = crate::qual::mean_prob_q(&qual[start..end]) as f32;
            updates.push((qs, Some(Value::Float(mean))));
        }
        // PacBio's `qs:i`/`qe:i` are with respect to the original read and
        // follow the window; one without the other is left as is.
        if let Some((qs0, _)) = query_span(src) {
            let (qs_new, qe_new) = window_coords(qs0, start, end);
            if let (Some(qs_value), Some(qe_value)) =
                (signal_int_value(qs_new), signal_int_value(qe_new))
            {
                updates.push((qs, Some(qs_value)));
                updates.push((Tag::new(b'q', b'e'), Some(qe_value)));
            }
        }
    }
    if split {
        let du = Tag::new(b'd', b'u');
        match signal {
            Some(window) => updates.extend(split_time_updates(src, window)),
            None => {
                updates.push((Tag::new(b's', b't'), None));
                if matches!(src.data().get(&du), Some(Value::Float(_))) {
                    updates.push((du, None));
                }
            },
        }
    }
    if split && platform == Platform::Ont {
        // Dorado's subread convention (`splitter_utils.cpp`): read number -1,
        // the parent read id, zero MinKNOW events, and an unknown end reason on
        // subreads that do not retain the end of the parent signal. Without
        // a move index, the last sequence segment retains the source end reason.
        let me = Tag::new(b'm', b'e');
        let er = Tag::new(b'e', b'r');
        updates.push((Tag::new(b'r', b'n'), Some(Value::Int32(-1))));
        set_update(
            &mut updates,
            Tag::new(b'p', b'i'),
            Some(Value::String(parent_read_id(src).into())),
        );
        if src.data().get(&me).is_some() {
            updates.push((me, Some(Value::Int32(0))));
        }
        let retains_end = moves.map_or(idx + 1 == total, |index| {
            if index.reversed {
                start == 0
            } else {
                end == orig_len
            }
        });
        if !retains_end && src.data().get(&er).is_some() {
            updates.push((er, Some(Value::String(b"unknown".as_slice().into()))));
        }
    }
    updates
}

/// Bumps `Counters::undo_tags_dropped_reads` when a read carrying a PacBio
/// undo blob (`ds`/`ls`) loses it in the output: at least one written window
/// does not span the whole read.
pub(super) fn count_undo_tags_dropped(
    counters: &Counters,
    src: &RecordBuf,
    orig_len: usize,
    survivors: &[(usize, usize)],
) {
    let trimmed = survivors.iter().any(|&(s, e)| s != 0 || e != orig_len);
    let carries_undo = UNDO_TAGS
        .iter()
        .any(|t| src.data().get(&Tag::new(t[0], t[1])).is_some());
    if trimmed && carries_undo {
        counters
            .undo_tags_dropped_reads
            .fetch_add(1, Ordering::Relaxed);
    }
}
