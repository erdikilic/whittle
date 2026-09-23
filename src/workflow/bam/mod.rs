//! uBAM workflows: record reconstruction (sequence, quality, MM/ML/MN, per-base and signal tags) and the sequential, parallel and raw full-window drivers for BAM and FASTQ output.

use std::borrow::Cow;
use std::io;
use std::sync::Arc;
use std::sync::atomic::Ordering;

use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};

use super::reject::{self, Reason, RejectItem};
use super::{
    BAM_BATCH, BatchSink, Counters, Rejection, Stats, process_read_segments, run_bytes_parallel,
    run_parallel,
};
use crate::config::{Config, FastqTags, TagRemoval};
use crate::io::fastq::{push_aux_field, push_mods_aux, push_record_body};
use crate::mods::reconstruct::IndexedMods;
use crate::{mods, trim};

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

/// ONT signal-mapping tags: the `mv` move table plus the `ts`/`ns` sample counts
/// and the `sp`/`pi` split linkage. On a trimmed read these are either rewritten
/// (`--update-moves`) or dropped (default), never left stale. Handled by
/// `signal_tag_updates`, not the per-base pass.
pub(crate) use crate::config::SIGNAL_TAGS;

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

/// The base-modification block, in the order it is emitted.
const MOD_TAGS: [Tag; 3] = [
    Tag::BASE_MODIFICATIONS,
    Tag::BASE_MODIFICATION_PROBABILITIES,
    Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
];

/// Converts a raw record to a `RecordBuf` on the render worker without routing
/// sequence, quality and every aux value through the generic SAM trait
/// iterators. The concrete noodles views have bulk conversions for these large
/// fields and reduce conversion overhead on long reads.
///
/// BAM's `CG:B:I` overflow representation of a CIGAR longer than 65535
/// operations is not expanded: the workflows accept unaligned records only,
/// whose CIGAR is empty.
pub(crate) fn decode_raw_record(src: &bam::Record) -> std::io::Result<RecordBuf> {
    let mut dst = RecordBuf::default();
    *dst.name_mut() = src.name().map(Into::into);
    *dst.flags_mut() = src.flags();
    *dst.reference_sequence_id_mut() = src.reference_sequence_id().transpose()?;
    *dst.alignment_start_mut() = src.alignment_start().transpose()?;
    *dst.mapping_quality_mut() = src.mapping_quality();

    let cigar = dst.cigar_mut().as_mut();
    cigar.clear();
    for result in src.cigar().iter() {
        cigar.push(result?);
    }

    *dst.mate_reference_sequence_id_mut() = src.mate_reference_sequence_id().transpose()?;
    *dst.mate_alignment_start_mut() = src.mate_alignment_start().transpose()?;
    *dst.template_length_mut() = src.template_length();
    *dst.sequence_mut() = src.sequence().into();
    *dst.quality_scores_mut() = src.quality_scores().into();
    *dst.data_mut() = src.data().try_into()?;
    Ok(dst)
}

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
fn slice_array(a: &Array, start: usize, end: usize) -> Array {
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
fn perbase_slice(
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
fn array_integers(a: &Array) -> Option<Vec<i64>> {
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
fn array_at_subtype(template: &Array, values: &[i64]) -> Option<Array> {
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
fn rle_runs_len(runs: &[i64]) -> Option<usize> {
    if !runs.len().is_multiple_of(2) {
        return None;
    }
    runs.iter().step_by(2).try_fold(0usize, |sum, &len| {
        sum.checked_add(usize::try_from(len).ok()?)
    })
}

/// Returns the number of bases an `sa` run-length coverage array covers,
/// `None` when it is not a well-formed run list (`rle_runs_len`).
fn rle_coverage_len(a: &Array) -> Option<usize> {
    rle_runs_len(&array_integers(a)?)
}

/// Slices an `sa` run-length coverage array to the window `[start, end)`: each
/// run is clipped to the window and adjacent runs of equal coverage are merged.
/// `None` when the runs do not cover exactly `orig_len` bases, which leaves the
/// tag unchanged.
fn slice_rle_coverage(a: &Array, orig_len: usize, start: usize, end: usize) -> Option<Array> {
    let runs = array_integers(a)?;
    if rle_runs_len(&runs)? != orig_len {
        return None;
    }
    let mut out: Vec<i64> = Vec::new();
    let mut pos = 0usize;
    for run in runs.chunks_exact(2) {
        // `rle_runs_len` has checked that every length is a non-negative
        // `usize` and that the sum fits.
        let (len, coverage) = (run[0] as usize, run[1]);
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
fn windowed_value(
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
fn aux_integer(value: &Value) -> Option<i64> {
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

/// Dorado's `bi` barcode-info tag: a `B:f` array of exactly seven floats,
/// `[barcode_score, front_start_index, front_len, front_score, rear_end_index,
/// rear_len, rear_score]` (`read_pipeline/base/messages.cpp`).
const BARCODE_TAG: [u8; 2] = *b"bi";

/// The barcode spans a record's `bi` positions describe.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BarcodeSpan {
    /// The record carries no `bi` tag, so the barcode stage keeps every base.
    Absent,
    /// The recorded front and rear barcode spans, each `[start, end)` in read
    /// coordinates; a barcode dorado did not find is `None`.
    Spans {
        front: Option<(usize, usize)>,
        rear: Option<(usize, usize)>,
    },
    /// The tag is not a seven-element float array, or its positions describe an
    /// empty, inverted or out-of-range window. The read is left untrimmed by
    /// this stage and counted in `Counters::barcode_tag_malformed_reads`.
    Malformed,
}

/// Returns the integer a `bi` position holds. `f32` represents every base index
/// a read can hold exactly and the cast saturates, so only a non-finite value
/// is rejected.
fn barcode_position(value: f32) -> Option<i64> {
    value.is_finite().then_some(value as i64)
}

/// Returns the barcode spans a `bi` position pair describes over a
/// `seq_len`-base read.
///
/// `front_start + front_len` is the last base of the front barcode and
/// `rear_end - rear_len` is the first base of the rear one, so dorado's own
/// trimmer keeps `[front_start + front_len + 1, rear_end - rear_len)`
/// (`demux/Trimmer.cpp::determine_trim_interval`). Each end is guarded on its
/// own raw position: a barcode dorado did not find is written as `-1`, which
/// leaves that end at the read's boundary.
fn barcode_interval(
    front_start: f32,
    front_len: f32,
    rear_end: f32,
    rear_len: f32,
    seq_len: usize,
) -> BarcodeSpan {
    let (Some(front_start), Some(front_len), Some(rear_end), Some(rear_len)) = (
        barcode_position(front_start),
        barcode_position(front_len),
        barcode_position(rear_end),
        barcode_position(rear_len),
    ) else {
        return BarcodeSpan::Malformed;
    };
    let len = i64::try_from(seq_len).unwrap_or(i64::MAX);
    let front = (front_start >= 0).then(|| {
        (
            front_start,
            front_start.saturating_add(front_len).saturating_add(1),
        )
    });
    let rear = (rear_end >= 0).then(|| (rear_end.saturating_sub(rear_len), rear_end));
    let start = front.map_or(0, |(_, end)| end);
    let end = rear.map_or(len, |(start, _)| start);
    let valid = |(s, e): (i64, i64)| s >= 0 && s < e && e <= len;
    if start >= end || !front.is_none_or(valid) || !rear.is_none_or(valid) {
        return BarcodeSpan::Malformed;
    }
    let cast = |(s, e): (i64, i64)| (s as usize, e as usize);
    BarcodeSpan::Spans {
        front: front.map(cast),
        rear: rear.map(cast),
    }
}

/// Returns the record's barcode call (`BC`), the value dorado writes as
/// `<kit>_barcodeNN`.
fn barcode_call(rec: &RecordBuf) -> Option<&[u8]> {
    match rec.data().get(&Tag::new(b'B', b'C'))? {
        Value::String(value) => Some(value.as_ref()),
        _ => None,
    }
}

/// Resolves the retained window from a record's verified barcode spans: a
/// span is trimmed only when a barcode sequence is found at it. Returns the
/// window and whether any recorded span failed verification.
fn verified_barcode_window(
    rec: &RecordBuf,
    seq: &[u8],
    adapters: &crate::adapter::AdapterConfig,
    front: Option<(usize, usize)>,
    rear: Option<(usize, usize)>,
) -> (Option<(usize, usize)>, bool) {
    let call = barcode_call(rec);
    let verify = |span: Option<(usize, usize)>| {
        span.map(|(s, e)| (s, e, adapters.barcode_span_verified(seq, s, e, call)))
    };
    let (front, rear) = (verify(front), verify(rear));
    let unverified = front.is_some_and(|(_, _, ok)| !ok) || rear.is_some_and(|(_, _, ok)| !ok);
    let start = front.filter(|f| f.2).map_or(0, |(_, e, _)| e);
    let end = rear.filter(|r| r.2).map_or(seq.len(), |(s, _, _)| s);
    (
        (start < end && (start > 0 || end < seq.len())).then_some((start, end)),
        unverified,
    )
}

/// Resolves a record's `bi` barcode positions into barcode spans over a
/// `seq_len`-base sequence.
pub(crate) fn barcode_window(rec: &RecordBuf, seq_len: usize) -> BarcodeSpan {
    let Some(value) = rec.data().get(&Tag::new(BARCODE_TAG[0], BARCODE_TAG[1])) else {
        return BarcodeSpan::Absent;
    };
    let Value::Array(Array::Float(values)) = value else {
        return BarcodeSpan::Malformed;
    };
    let &[
        _score,
        front_start,
        front_len,
        _front_score,
        rear_end,
        rear_len,
        _rear_score,
    ] = values.as_slice()
    else {
        return BarcodeSpan::Malformed;
    };
    barcode_interval(front_start, front_len, rear_end, rear_len, seq_len)
}

/// The state of a record's `MM`/`ML`/`MN` block relative to its sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ModBlock {
    /// No `MM:Z` tag. An `ML` or `MN` present on its own is copied verbatim.
    Absent,
    /// `MM` parses to its end, `ML` (when present) is a `B:C` array of the
    /// length `MM` declares, all calls fit the sequence, and `MN` equals its length.
    Consistent,
    /// As `Consistent`, with `MN` absent; the output gains one.
    MissingMn,
    /// `MM` does not parse to its end, `ML` is not a `B:C` array or has the
    /// wrong length, `MN` disagrees with the sequence length, or a call exceeds
    /// the available counting-base occurrences. The block is removed from the
    /// output and the read is counted in `Counters::malformed_mod_reads`.
    Malformed,
}

/// Checks modification syntax and tag lengths; sequence positions are checked
/// separately. `ml` is `None` when the tag
/// is absent and `Some(None)` when it is present with a subtype other than
/// `B:C`; `mn` is `None` when absent and `Some(None)` when not an integer.
fn classify_mod_block(
    mm: &[u8],
    ml: Option<Option<usize>>,
    mn: Option<Option<i64>>,
    seq_len: usize,
) -> ModBlock {
    let Some(expected) = mods::expected_ml_len(mm) else {
        return ModBlock::Malformed;
    };
    match ml {
        Some(None) => return ModBlock::Malformed,
        Some(Some(len)) if len != expected => return ModBlock::Malformed,
        _ => {},
    }
    match mn {
        None => ModBlock::MissingMn,
        Some(Some(n)) if i64::try_from(seq_len).ok() == Some(n) => ModBlock::Consistent,
        Some(_) => ModBlock::Malformed,
    }
}

/// Returns the `MM` bytes of a record and, when it carries a `B:C` array, its
/// `ML` bytes. `None` when `MM` is absent or not a string.
fn mod_tags(src: &RecordBuf) -> Option<(&[u8], Option<&[u8]>)> {
    let mm: &[u8] = match src.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => AsRef::<[u8]>::as_ref(s),
        _ => return None,
    };
    let ml = match src.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
        Some(Value::Array(Array::UInt8(v))) => Some(v.as_slice()),
        _ => None,
    };
    Some((mm, ml))
}

/// Classifies the modification block of a decoded record whose sequence has
/// `seq_len` bases.
pub(crate) fn inspect_mod_block(src: &RecordBuf, seq_len: usize) -> ModBlock {
    let Some((mm, _)) = mod_tags(src) else {
        return ModBlock::Absent;
    };
    let ml = src
        .data()
        .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
        .map(|v| match v {
            Value::Array(Array::UInt8(a)) => Some(a.len()),
            _ => None,
        });
    let mn = src
        .data()
        .get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH)
        .map(aux_integer);
    let block = classify_mod_block(mm, ml, mn, seq_len);
    if block != ModBlock::Malformed
        && !mods::parse::positions_valid(mm, src.sequence().as_ref().iter().copied())
    {
        ModBlock::Malformed
    } else {
        block
    }
}

/// Parses an `mv` move table value into `(stride, moves)`. `None` unless it is a
/// `B:c` (Int8) array with a positive stride. `moves` excludes the stride prefix;
/// each entry corresponds to `stride` signal samples (1 = a base emitted here, so
/// the count of 1s equals the sequence length).
pub(crate) fn parse_move_table(value: &Value) -> Option<(i8, &[i8])> {
    match value {
        Value::Array(Array::Int8(a)) => {
            let (stride, moves) = a.split_first()?;
            if *stride > 0 {
                Some((*stride, moves))
            } else {
                None
            }
        },
        _ => None,
    }
}

/// Signal direction resolved from the record's basecalling model.
fn signal_reversed(header: &sam::Header, rec: &RecordBuf) -> Option<bool> {
    use sam::header::record::value::map::read_group::tag::DESCRIPTION;
    let direction =
        |group: &sam::header::record::value::Map<sam::header::record::value::map::ReadGroup>| {
            let description = group.other_fields().get(&DESCRIPTION)?;
            description
                .split(|b| b.is_ascii_whitespace() || *b == b';')
                .find_map(|field| {
                    let model = field.strip_prefix(b"basecall_model=")?;
                    if model.starts_with(b"rna") {
                        Some(true)
                    } else if model.starts_with(b"dna") {
                        Some(false)
                    } else {
                        None
                    }
                })
        };
    match rec.data().get(&Tag::READ_GROUP) {
        Some(Value::String(id)) => direction(header.read_groups().get(AsRef::<[u8]>::as_ref(id))?),
        None => {
            let mut groups = header.read_groups().values();
            let first = direction(groups.next()?)?;
            groups
                .all(|group| direction(group) == Some(first))
                .then_some(first)
        },
        _ => None,
    }
}

/// Emitted-base boundaries in signal block order, indexed once per read.
struct MoveIndex<'a> {
    stride: i8,
    moves: &'a [i8],
    boundaries: Vec<usize>,
    reversed: bool,
}

impl<'a> MoveIndex<'a> {
    /// Accepts binary move tables with one emission per sequence base.
    fn new(src: &'a RecordBuf, reversed: bool) -> Option<Self> {
        let (stride, moves) = src
            .data()
            .get(&Tag::new(b'm', b'v'))
            .and_then(parse_move_table)?;
        let mut boundaries = Vec::new();
        for (i, &m) in moves.iter().enumerate() {
            match m {
                0 => {},
                1 => boundaries.push(i),
                _ => return None,
            }
        }
        if boundaries.len() != src.sequence().len() {
            return None;
        }
        boundaries.push(moves.len());
        Some(Self {
            stride,
            moves,
            boundaries,
            reversed,
        })
    }
}

/// Reads an integer aux tag as `i64`, regardless of stored width.
fn signal_int(src: &RecordBuf, tag: &[u8; 2]) -> Option<i64> {
    src.data()
        .get(&Tag::new(tag[0], tag[1]))
        .and_then(aux_integer)
}

fn signal_offset(blocks: usize, stride: usize) -> Option<i64> {
    blocks
        .checked_mul(stride)
        .and_then(|n| i64::try_from(n).ok())
}

fn signal_int_value(n: i64) -> Option<Value> {
    if let Ok(n) = i32::try_from(n) {
        Some(Value::Int32(n))
    } else if let Ok(n) = u32::try_from(n) {
        Some(Value::UInt32(n))
    } else {
        None
    }
}

/// Returns the parent read id for a subread: the source's own `pi` if it has
/// one (so `pi` always names the ultimate ancestor, matching dorado), else the
/// source read name.
fn parent_read_id(src: &RecordBuf) -> Vec<u8> {
    match src.data().get(&Tag::new(b'p', b'i')) {
        Some(Value::String(s)) => s.to_vec(),
        _ => src.name().map(|n| n.to_vec()).unwrap_or_default(),
    }
}

/// The platform whose tag conventions a record follows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Platform {
    /// PacBio: integer `qs`/`qe` query coordinates, `{movie}/{zmw}/...` read
    /// names, `rn` as a pass count, `du:Z` from pbmarkdup.
    PacBio,
    /// ONT (dorado): float `qs`, the `mv`/`ts`/`ns` signal tags, `st`/`du`
    /// read timing, `_segment_N` split names.
    Ont,
}

/// Classifies a record: `PacBio` when it carries an integer `qs` (dorado's `qs`
/// is a float) or its name follows a PacBio convention (`parse_pacbio_name`),
/// `Ont` otherwise.
pub(crate) fn platform(rec: &RecordBuf) -> Platform {
    let integer_qs = rec
        .data()
        .get(&Tag::new(b'q', b's'))
        .and_then(aux_integer)
        .is_some();
    let pacbio_name = rec
        .name()
        .is_some_and(|n| parse_pacbio_name(AsRef::<[u8]>::as_ref(n)).is_some());
    if integer_qs || pacbio_name {
        Platform::PacBio
    } else {
        Platform::Ont
    }
}

/// The parts of a PacBio read name a split segment name is built from.
struct PacBioName<'a> {
    /// The name without any `/{qStart}_{qEnd}` interval: `{movie}/{zmw}` for
    /// a subread, `{movie}/{zmw}/ccs` for a CCS read, with `/fwd` or `/rev`
    /// for a by-strand read.
    stem: &'a [u8],
    /// The `qStart` of a `{qStart}_{qEnd}` interval in the name.
    query_start: Option<i64>,
}

/// Parses a read name against the PacBio BAM conventions: `{movie}/{zmw}/ccs`
/// with an optional `/fwd` or `/rev` and an optional `/{qStart}_{qEnd}`
/// (segmented reads), or the subread form `{movie}/{zmw}/{qStart}_{qEnd}`.
/// `None` for any other name.
fn parse_pacbio_name(name: &[u8]) -> Option<PacBioName<'_>> {
    fn is_digits(s: &[u8]) -> bool {
        !s.is_empty() && s.iter().all(u8::is_ascii_digit)
    }
    fn interval_start(s: &[u8]) -> Option<i64> {
        let (start, end) = s.split_at(s.iter().position(|&c| c == b'_')?);
        if !is_digits(start) || !is_digits(&end[1..]) {
            return None;
        }
        std::str::from_utf8(start).ok()?.parse().ok()
    }

    let parts: Vec<&[u8]> = name.split(|&c| c == b'/').collect();
    let [movie, zmw, rest @ ..] = parts.as_slice() else {
        return None;
    };
    if movie.is_empty() || !is_digits(zmw) {
        return None;
    }
    match rest {
        [b"ccs"] | [b"ccs", b"fwd" | b"rev"] => Some(PacBioName {
            stem: name,
            query_start: None,
        }),
        [interval] | [b"ccs", interval] | [b"ccs", b"fwd" | b"rev", interval] => {
            let query_start = interval_start(interval)?;
            Some(PacBioName {
                stem: &name[..name.len() - interval.len() - 1],
                query_start: Some(query_start),
            })
        },
        _ => None,
    }
}

/// Returns a record's integer `qs`/`qe` query coordinates, `None` unless both
/// are present as integers.
fn query_span(src: &RecordBuf) -> Option<(i64, i64)> {
    Some((signal_int(src, b"qs")?, signal_int(src, b"qe")?))
}

/// Returns the query coordinates of window `[start, end)` in the frame of the
/// original PacBio read whose query starts at `qs0`: the PacBio BAM spec keeps
/// `qs`/`qe` with respect to the original read through clipping.
fn window_coords(qs0: i64, start: usize, end: usize) -> (i64, i64) {
    (qs0 + start as i64, qs0 + end as i64)
}

/// Returns the output name of `window`, updating existing PacBio query intervals.
/// A split names an ONT segment `{name}_segment_{n}` (1-based) and a
/// PacBio segment `{stem}/{qStart}_{qEnd}` from `coords`, the segment's
/// rewritten `qs`/`qe`; without them the interval is offset from the name's
/// own `qStart` (0 when the name has none). A PacBio record whose name follows
/// no PacBio convention takes the ONT suffix.
pub(crate) fn segment_name(
    platform: Platform,
    name: &[u8],
    window: Window,
    coords: Option<(i64, i64)>,
) -> Vec<u8> {
    if platform == Platform::PacBio
        && let Some(parts) = parse_pacbio_name(name)
        && (window.total > 1 || parts.query_start.is_some())
    {
        let (qs, qe) = coords.unwrap_or_else(|| {
            window_coords(parts.query_start.unwrap_or(0), window.start, window.end)
        });
        let mut out = parts.stem.to_vec();
        out.extend_from_slice(format!("/{qs}_{qe}").as_bytes());
        return out;
    }
    if window.total <= 1 {
        return name.to_vec();
    }
    let mut out = name.to_vec();
    out.extend_from_slice(format!("_segment_{}", window.idx + 1).as_bytes());
    out
}

/// Computes the poly-A tag updates (`pa` signal boundaries, `pt` tail length)
/// for a trimmed read. `pa` holds absolute original-signal positions, the frame
/// `ts` and `ns` use; `-1`/`-2` are dorado's not-found/not-enabled sentinels and
/// are left as is. When every real position falls inside
/// `[kept_start, kept_end)` the tail survived: a split shifts `pa` into the
/// subread's own signal frame, a crop keeps both unchanged. Otherwise, or with
/// no poly-A array, `pa`/`pt` are dropped.
fn polya_updates(
    src: &RecordBuf,
    kept_start: i64,
    kept_end: i64,
    is_split: bool,
) -> Vec<(Tag, Option<Value>)> {
    let pa_tag = Tag::new(b'p', b'a');
    let pt_tag = Tag::new(b'p', b't');
    let drop_both = || vec![(pa_tag, None), (pt_tag, None)];

    let pa = match src.data().get(&pa_tag) {
        Some(Value::Array(Array::Int32(v))) => v,
        _ => return drop_both(),
    };
    // `pa` = [anchor, range0.start, range0.end, range1.start, range1.end].
    // Dorado's poly-A signal ranges are half-open `[start, end)`: the anchor and
    // the range starts are inclusive sample indexes and must be `< kept_end`; the
    // range ends are exclusive and may equal `kept_end`. Every real position must
    // also be `>= kept_start`. Sentinels (`< 0`) are skipped.
    let has_real = pa.iter().any(|&p| p >= 0);
    let survives = has_real
        && pa.iter().enumerate().all(|(i, &p)| {
            if p < 0 {
                return true; // sentinel (NOT_FOUND / NOT_ENABLED)
            }
            let p = i64::from(p);
            let within_upper = if i == 2 || i == 4 {
                p <= kept_end
            } else {
                p < kept_end
            };
            p >= kept_start && within_upper
        });
    if !survives {
        return drop_both();
    }
    if is_split {
        // Shifted into the subread's own frame (subread signal 0 is `kept_start`;
        // its `ts` is 0). Sentinels stay unchanged, as does `pt` (a base count).
        let mut shifted = Vec::with_capacity(pa.len());
        for &p in pa {
            if p >= 0 {
                let Some(q) = i64::from(p)
                    .checked_sub(kept_start)
                    .and_then(|n| i32::try_from(n).ok())
                else {
                    return drop_both();
                };
                shifted.push(q);
            } else {
                shifted.push(p);
            }
        }
        vec![(pa_tag, Some(Value::Array(Array::Int32(shifted))))]
    } else {
        // A crop keeps `pa`/`pt`: absolute original-signal positions remain valid.
        Vec::new()
    }
}

/// The original-signal window a trimmed read's kept bases span, in the frame
/// `ts`/`ns` use: samples `[kept_start, kept_end)`.
#[derive(Debug, Clone, Copy)]
struct SignalWindow {
    /// First kept sample, inclusive.
    kept_start: i64,
    /// End of the kept signal, exclusive.
    kept_end: i64,
}

/// Aux tag updates for one output window: `(tag, Some(value))` replaces the
/// tag in place or appends it when the source lacks it, `(tag, None)` removes
/// it.
type TagUpdates = Vec<(Tag, Option<Value>)>;

/// Computes the ONT signal tag updates for output window `[start, end)`.
/// Returns `(tag, Some(value))` to set or `(tag, None)` to remove, with the
/// kept signal window when the tags were rewritten; empty when the read is not
/// trimmed. With `update_moves` off, or a missing or malformed move table, the
/// five signal tags and both poly-A tags are removed. With it on, `mv` is
/// sliced by block range (stride-aligned, following dorado
/// `splitter::subread`) and:
///   - crop (`total == 1`, name kept): `ts += block_first*stride`; `ns` is the
///     kept signal's end, which is the source `ns` when the window runs to the
///     last base.
///   - split (`total > 1`, renamed): `ts = 0`, `ns = span`,
///     `sp = parent offset`, `pi = parent id`.
fn signal_tag_updates(
    src: &RecordBuf,
    seq_len: usize,
    start: usize,
    end: usize,
    total: usize,
    moves: Option<&MoveIndex<'_>>,
) -> (TagUpdates, Option<SignalWindow>) {
    if start == 0 && end == seq_len {
        return (Vec::new(), None);
    }
    if let Some(moves) = moves
        && let Some((updates, window)) = signal_rewrite(src, seq_len, start, end, total, moves)
    {
        return (updates, Some(window));
    }
    let dropped = SIGNAL_TAGS
        .iter()
        .chain(POLYA_TAGS.iter())
        .map(|t| (Tag::new(t[0], t[1]), None))
        .collect();
    (dropped, None)
}

/// Rewrites the signal and poly-A tags of a trimmed window from the move
/// table. `None` when the table is missing or malformed, its base count
/// disagrees with the sequence, the window has no start base, or a signal
/// offset does not fit its tag; the caller then removes the tags.
fn signal_rewrite(
    src: &RecordBuf,
    seq_len: usize,
    start: usize,
    end: usize,
    total: usize,
    index: &MoveIndex<'_>,
) -> Option<(TagUpdates, SignalWindow)> {
    if start >= end || end > seq_len {
        return None;
    }
    let (start, end) = if index.reversed {
        (seq_len - end, seq_len - start)
    } else {
        (start, end)
    };
    let stride = index.stride;
    let moves = index.moves;
    let stride_n = stride as usize;
    let block_first = *index.boundaries.get(start)?;
    let block_second = *index.boundaries.get(end)?;

    let mut new_mv = Vec::with_capacity(1 + block_second - block_first);
    new_mv.push(stride);
    new_mv.extend_from_slice(&moves[block_first..block_second]);
    let mut updates = vec![(
        Tag::new(b'm', b'v'),
        Some(Value::Array(Array::Int8(new_mv))),
    )];

    // Original-signal window the kept bases span: [ts0 + block_first*stride,
    // ts0 + block_second*stride). `ns = span + front trim` matches dorado's
    // `ns = raw_data_samples + num_trimmed_samples` (a tail crop shrinks ns, a
    // head-only crop leaves it unchanged, a split gets the subread span).
    let ts0 = signal_int(src, b"ts").unwrap_or(0);
    let kept_start = ts0.checked_add(signal_offset(block_first, stride_n)?)?;
    let block_end = ts0.checked_add(signal_offset(block_second, stride_n)?)?;
    // The move table resolves the signal end only to the stride; the source
    // `ns` names it exactly when the window runs to the last base.
    let kept_end = match signal_int(src, b"ns") {
        Some(ns0) if end == seq_len && ns0 > block_end => ns0,
        _ => block_end,
    };

    if total > 1 {
        // A split yields a dorado subread: renamed, front trim reset to 0, parent
        // linkage set. Dorado's `split_point` is the parent's own plus the
        // parent's trimmed samples plus the signal offset, so `sp` counts from
        // the parent's POD5 signal start, not from its first basecalled sample.
        let sp = signal_int(src, b"sp")
            .unwrap_or(0)
            .checked_add(kept_start)?;
        let ns_value = signal_int_value(kept_end - kept_start)?;
        let sp_value = signal_int_value(sp)?;
        let pi = parent_read_id(src);
        updates.push((Tag::new(b't', b's'), Some(Value::Int32(0))));
        updates.push((Tag::new(b'n', b's'), Some(ns_value)));
        updates.push((Tag::new(b's', b'p'), Some(sp_value)));
        updates.push((Tag::new(b'p', b'i'), Some(Value::String(pi.into()))));
    } else {
        // A head or tail crop keeps the read identity and advances the front trim.
        let ts_value = signal_int_value(kept_start)?;
        let ns_value = signal_int_value(kept_end)?;
        updates.push((Tag::new(b't', b's'), Some(ts_value)));
        updates.push((Tag::new(b'n', b's'), Some(ns_value)));
    }
    updates.extend(polya_updates(src, kept_start, kept_end, total > 1));
    Some((
        updates,
        SignalWindow {
            kept_start,
            kept_end,
        },
    ))
}

/// Recomputes `st`/`du` for a split subread whose signal window is known,
/// following dorado `splitter_utils.cpp`: the sample rate is the source `ns`
/// over its `du`, the subread duration is its sample count over that rate, and
/// its start time is the source `st` advanced by `kept_start` samples. A `du`
/// that is not a positive float, a missing `ns`, or an `st` that does not
/// parse leaves the tag unchanged.
fn split_time_updates(src: &RecordBuf, window: SignalWindow) -> TagUpdates {
    let du_tag = Tag::new(b'd', b'u');
    let st_tag = Tag::new(b's', b't');
    let mut updates = Vec::new();
    let Some(Value::Float(du0)) = src.data().get(&du_tag) else {
        return updates;
    };
    let Some(ns0) = signal_int(src, b"ns") else {
        return updates;
    };
    if ns0 <= 0 || !du0.is_finite() || *du0 <= 0.0 {
        return updates;
    }
    let rate = ns0 as f64 / f64::from(*du0);
    let samples = (window.kept_end - window.kept_start) as f64;
    updates.push((du_tag, Some(Value::Float((samples / rate) as f32))));
    if let Some(Value::String(st0)) = src.data().get(&st_tag)
        && let Some(shifted) = shift_timestamp(st0.as_ref(), window.kept_start as f64 / rate)
    {
        updates.push((st_tag, Some(Value::String(shifted.into()))));
    }
    updates
}

/// Returns an RFC 3339 `st` value advanced by `seconds`, at millisecond
/// precision (dorado's own) and in the offset form the source uses: `Z`, a
/// numeric offset, or none (a civil time read as UTC). `None` when the source
/// does not parse.
fn shift_timestamp(value: &[u8], seconds: f64) -> Option<Vec<u8>> {
    use jiff::fmt::temporal::{Pieces, PiecesOffset};
    use jiff::tz::Offset;

    let text = std::str::from_utf8(value).ok()?;
    let pieces = Pieces::parse(text).ok()?;
    let datetime = pieces.date().to_datetime(pieces.time().unwrap_or_default());
    let offset = pieces.offset();
    let instant = offset
        .as_ref()
        .map_or(Offset::UTC, PiecesOffset::to_numeric_offset)
        .to_timestamp(datetime)
        .ok()?;
    let millis = (seconds * 1000.0).round();
    if !millis.is_finite() {
        return None;
    }
    let shifted = instant
        .checked_add(jiff::SignedDuration::from_millis(millis as i64))
        .ok()?;
    let out = match offset {
        Some(PiecesOffset::Zulu) => format!("{shifted:.3}"),
        Some(o) => format!("{:.3}", shifted.display_with_offset(o.to_numeric_offset())),
        None => format!("{:.3}", Offset::UTC.to_datetime(shifted)),
    };
    Some(out.into_bytes())
}

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

/// One output window of a read: bases `[start, end)`, segment `idx` (0-based)
/// of `total`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Window {
    /// First base of the window, inclusive.
    pub start: usize,
    /// End of the window, exclusive.
    pub end: usize,
    /// 0-based segment index.
    pub idx: usize,
    /// Number of segments produced from the read.
    pub total: usize,
}

/// Returns whether `t` is handled by a dedicated rule rather than the structural
/// per-base slice: the modification block, the signal and poly-A tags, and the
/// tags dropped on a trim or a split.
fn has_dedicated_rule(t: [u8; 2]) -> bool {
    matches!(&t, b"MM" | b"ML" | b"MN")
        || SIGNAL_TAGS.contains(&t)
        || POLYA_TAGS.contains(&t)
        || DROP_ON_TRIM_TAGS.contains(&t)
        || DROP_ON_SPLIT_TAGS.contains(&t)
}

/// Replaces the update for `tag` or appends one.
fn set_update(updates: &mut TagUpdates, tag: Tag, value: Option<Value>) {
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
fn window_tag_updates(
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
fn count_undo_tags_dropped(
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

/// Builds one output uBAM record for `window`: SEQ/QUAL sliced, `MM`/`ML`/`MN`
/// rebuilt, per-base kinetics sliced, stale signal-space tags rewritten or
/// dropped, and the name updated for splits and PacBio interval crops. The record
/// is assembled field by field: aux tags are copied in source order with the
/// rewritten ones replaced in place, removed ones skipped and added ones
/// appended. A `Malformed` block is removed. An untrimmed, unsplit record with
/// an `Absent` or `Consistent` block and nothing to remove is returned as
/// `None`: its output is its input.
///
/// `remove` names the tags `--remove-tag` drops. Removal
/// is applied last, to the rewritten tag set, so a removed tag whittle
/// maintains (`MM`, the move table, a per-base array) is absent from the
/// output rather than left stale.
fn reconstruct_window_record(
    src: &RecordBuf,
    window: Window,
    mod_block: ModBlock,
    indexed: Option<&IndexedMods>,
    moves: Option<&MoveIndex<'_>>,
    remove: &TagRemoval,
) -> Option<RecordBuf> {
    let Window {
        start, end, total, ..
    } = window;
    let seq = src.sequence().as_ref();
    let qual = src.quality_scores().as_ref();
    let orig_len = seq.len();
    let trimmed = start != 0 || end != orig_len;
    let split = total > 1;
    if !trimmed
        && !split
        && remove.is_empty()
        && matches!(mod_block, ModBlock::Absent | ModBlock::Consistent)
    {
        return None;
    }

    let platform = platform(src);
    let mut out = RecordBuf::default();
    *out.name_mut() = if trimmed || split {
        let name: &[u8] = src.name().map(AsRef::as_ref).unwrap_or_default();
        let coords = query_span(src).map(|(qs0, _)| window_coords(qs0, start, end));
        Some(segment_name(platform, name, window, coords).into())
    } else {
        src.name().map(Into::into)
    };
    *out.flags_mut() = src.flags();
    *out.reference_sequence_id_mut() = src.reference_sequence_id();
    *out.alignment_start_mut() = src.alignment_start();
    *out.mapping_quality_mut() = src.mapping_quality();
    *out.cigar_mut() = src.cigar().clone();
    *out.mate_reference_sequence_id_mut() = src.mate_reference_sequence_id();
    *out.mate_alignment_start_mut() = src.mate_alignment_start();
    *out.template_length_mut() = src.template_length();
    *out.sequence_mut() = seq[start..end].to_vec().into();
    *out.quality_scores_mut() = qual[start..end].to_vec().into();

    // Tags with dedicated handling: `Some` replaces the source value in place,
    // or is appended when the source lacks the tag; `None` removes it.
    let mut updates: Vec<(Tag, Option<Value>)> = Vec::new();
    match mod_block {
        ModBlock::Absent => {},
        ModBlock::Malformed => updates.extend(MOD_TAGS.map(|t| (t, None))),
        ModBlock::Consistent | ModBlock::MissingMn => {
            if let Some((mm, ml)) = mod_tags(src) {
                let (mm, ml) = rebuild_mods(mm, ml, seq, start, end, indexed);
                updates.push((Tag::BASE_MODIFICATIONS, Some(Value::String(mm.into()))));
                updates.push((
                    Tag::BASE_MODIFICATION_PROBABILITIES,
                    ml.map(|ml| Value::Array(Array::UInt8(ml))),
                ));
                updates.push((
                    Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
                    Some(Value::Int32((end - start) as i32)),
                ));
            }
        },
    }
    updates.extend(window_tag_updates(src, qual, window, platform, moves));
    // A removed tag is dropped from the rewrite list as well, so neither the
    // copy loop below nor the append loop after it can put it back.
    if !remove.is_empty() {
        updates.retain(|(t, _)| !remove.contains(&<[u8; 2]>::from(*t)));
    }

    let data = out.data_mut();
    for (tag, value) in src.data().iter() {
        if remove.contains(&<[u8; 2]>::from(tag)) {
            continue;
        }
        if let Some(i) = updates.iter().position(|(t, _)| *t == tag) {
            if let (_, Some(v)) = updates.remove(i) {
                data.insert(tag, v);
            }
            continue;
        }
        let t = <[u8; 2]>::from(tag);
        let sliced = if trimmed && !has_dedicated_rule(t) {
            windowed_value(t, value, orig_len, start, end)
        } else {
            None
        };
        data.insert(tag, sliced.unwrap_or_else(|| value.clone()));
    }
    for (tag, value) in updates {
        if let Some(v) = value {
            data.insert(tag, v);
        }
    }

    Some(out)
}

/// Rebuilds the `MM`/`ML` block of a `Consistent` or `MissingMn` record for the
/// window `[start, end)`: skip-counts renumbered, `ML` re-sliced. `ml` is
/// `None` when the source carries no `ML`, and so is the result's.
fn rebuild_mods(
    mm: &[u8],
    ml: Option<&[u8]>,
    seq: &[u8],
    start: usize,
    end: usize,
    indexed: Option<&IndexedMods>,
) -> (Vec<u8>, Option<Vec<u8>>) {
    // Over the full window the rebuild is the identity, so the source bytes are
    // returned and the parse, slice and re-serialize are skipped; that work is
    // the dominant cost of an untrimmed BAM-to-FASTQ run.
    if start == 0 && end == seq.len() {
        return (mm.to_vec(), ml.map(<[u8]>::to_vec));
    }
    let sliced = match indexed {
        Some(indexed) => indexed.window(start, end),
        None => mods::reconstruct(&mods::parse(mm, ml.unwrap_or(&[])), seq, start, end),
    };
    let (mm, probabilities) = mods::serialize(&sliced);
    (mm, ml.map(|_| probabilities))
}

/// The per-read state the decoded BAM workflows share.
struct PreparedRead<'a> {
    /// The record's bases.
    seq: &'a [u8],
    /// The record's per-base qualities.
    qual: &'a [u8],
    /// The state of the record's modification block.
    mod_block: ModBlock,
    /// The original-coordinate interval retained by barcode restriction, `None`
    /// without an adapter source or when the record carries no verified `bi`.
    barcode: Option<(usize, usize)>,
}

/// Runs the per-read guards and bookkeeping shared by the decoded workflows:
/// refuses aligned reads and legacy mod tags, requires full per-base quality,
/// classifies the modification block, counting a malformed one, counts a
/// malformed per-base tag, and resolves the barcode window from the verified
/// `bi` spans, counting an unusable or unverified `bi`.
fn prepare_read<'a>(
    rec: &'a RecordBuf,
    cfg: &Config,
    counters: &Counters,
) -> anyhow::Result<PreparedRead<'a>> {
    crate::io::bam::ensure_trimmable(rec)?;
    let seq = rec.sequence().as_ref();
    let qual = rec.quality_scores().as_ref();
    if qual.len() != seq.len() || crate::io::bam::quality_absent(qual) {
        anyhow::bail!(
            "read {}: BAM record SEQ length {} != QUAL length {} \
             (records without per-base quality are not supported)",
            crate::io::bam::display_name(rec.name().map(AsRef::as_ref)),
            seq.len(),
            qual.len()
        );
    }
    let mod_block = inspect_mod_block(rec, seq.len());
    if mod_block == ModBlock::Malformed {
        counters.malformed_mod_reads.fetch_add(1, Ordering::Relaxed);
    }
    if has_malformed_perbase_tag(rec, seq.len()) {
        counters.malformed_tag_reads.fetch_add(1, Ordering::Relaxed);
    }
    let barcode = match cfg.adapters.as_ref() {
        Some(adapters) => match barcode_window(rec, seq.len()) {
            BarcodeSpan::Spans { front, rear } => {
                let (window, unverified) = verified_barcode_window(rec, seq, adapters, front, rear);
                if unverified {
                    counters
                        .barcode_tag_unverified_reads
                        .fetch_add(1, Ordering::Relaxed);
                }
                window
            },
            BarcodeSpan::Absent => None,
            BarcodeSpan::Malformed => {
                counters
                    .barcode_tag_malformed_reads
                    .fetch_add(1, Ordering::Relaxed);
                None
            },
        },
        None => None,
    };
    Ok(PreparedRead {
        seq,
        qual,
        mod_block,
        barcode,
    })
}

/// Runs the per-read guards and the trim on a decoded record, filters each
/// produced segment and calls `render` with every window and the record's
/// modification block: `None` for a survivor, or the reason for a rejected
/// segment (a read that produced no segment is one rejected full window).
/// Rejected windows are rendered only while a rejected output is open.
/// Counts a dropped undo blob once the survivors are known. Shared by the BAM
/// and FASTQ output paths.
fn render_windows(
    rec: &RecordBuf,
    cfg: &Config,
    counters: &Counters,
    render: impl FnMut(Window, ModBlock, Option<&IndexedMods>, Option<Reason>) -> anyhow::Result<()>,
) -> anyhow::Result<()> {
    let PreparedRead {
        seq,
        qual,
        mod_block,
        barcode,
    } = prepare_read(rec, cfg, counters)?;
    let _read = crate::workflow::read_span(rec.name().map(|n| n.as_ref()).unwrap_or(b"<unnamed>"));
    let _read = _read.enter();
    let produced = trim::apply(seq, qual, &cfg.trim, cfg.adapters.as_ref(), barcode);
    let indexed =
        if matches!(mod_block, ModBlock::Consistent | ModBlock::MissingMn) && produced.len() > 1 {
            mod_tags(rec).map(|(mm, ml)| IndexedMods::new(mods::parse(mm, ml.unwrap_or(&[])), seq))
        } else {
            None
        };
    let mut survivors: Vec<(usize, usize)> = Vec::new();
    // Both callbacks render, so the closure is shared through a cell.
    let render = std::cell::RefCell::new(render);
    process_read_segments(
        &produced,
        seq,
        qual,
        &cfg.filter,
        counters,
        |idx, total, start, end| {
            survivors.push((start, end));
            render.borrow_mut()(
                Window {
                    start,
                    end,
                    idx,
                    total,
                },
                mod_block,
                indexed.as_ref(),
                None,
            )
        },
        |rejection| {
            if !counters.wants_rejects() {
                return Ok(());
            }
            let (window, reason) = match rejection {
                Rejection::Segment {
                    idx,
                    total,
                    start,
                    end,
                    reason,
                } => (
                    Window {
                        start,
                        end,
                        idx,
                        total,
                    },
                    Reason::Dropped(reason),
                ),
                Rejection::Whole => (
                    Window {
                        start: 0,
                        end: seq.len(),
                        idx: 0,
                        total: 1,
                    },
                    Reason::TrimmedToNothing,
                ),
            };
            render.borrow_mut()(window, mod_block, indexed.as_ref(), Some(reason))
        },
    )?;
    count_undo_tags_dropped(counters, rec, seq.len(), &survivors);
    Ok(())
}

/// Renders one decoded record for BAM output: every surviving window is
/// rebuilt into an output record and handed to `emit`. Shared by the
/// sequential and parallel drivers.
/// `emit` receives `None` for a window whose output record is the input
/// record.
fn render_bam_read(
    header: &sam::Header,
    rec: &RecordBuf,
    cfg: &Config,
    counters: &Counters,
    mut emit: impl FnMut(Option<RecordBuf>) -> io::Result<()>,
) -> anyhow::Result<()> {
    let direction = if cfg.update_moves {
        signal_reversed(header, rec)
    } else {
        None
    };
    // The move table is indexed on the first window that needs it.
    let moves: std::cell::OnceCell<Option<MoveIndex<'_>>> = std::cell::OnceCell::new();
    let seq = rec.sequence().as_ref();
    render_windows(rec, cfg, counters, |window, mod_block, indexed, reason| {
        let partial = window.start != 0 || window.end != seq.len();
        if cfg.update_moves
            && direction.is_none()
            && partial
            && rec.data().get(&Tag::new(b'm', b'v')).is_some()
        {
            anyhow::bail!(
                "read {}: --update-moves requires a DNA or RNA basecall_model in the @RG description",
                crate::io::bam::display_name(rec.name().map(AsRef::as_ref))
            );
        }
        let moves = if partial {
            moves
                .get_or_init(|| direction.and_then(|reverse| MoveIndex::new(rec, reverse)))
                .as_ref()
        } else {
            None
        };
        let out =
            reconstruct_window_record(rec, window, mod_block, indexed, moves, &cfg.remove_tags);
        match reason {
            None => Ok(emit(out)?),
            Some(reason) => {
                let mut rejected = out.unwrap_or_else(|| rec.clone());
                reject::tag_record(&mut rejected, reason);
                counters.reject(RejectItem::Bam(rejected))
            },
        }
    })
}

/// Renders a raw record rejected by the tag filter for BAM output.
pub(crate) fn tag_filtered_bam(record: &bam::Record) -> anyhow::Result<RejectItem> {
    let mut rec = decode_raw_record(record)?;
    reject::tag_record(&mut rec, Reason::TagFilter);
    Ok(RejectItem::Bam(rec))
}

/// Renders a raw record rejected by the tag filter for FASTQ output: the whole
/// read with its selected tags and the reason tag.
pub(crate) fn tag_filtered_bam_fastq(
    record: &bam::Record,
    cfg: &Config,
) -> anyhow::Result<RejectItem> {
    let rec = decode_raw_record(record)?;
    let seq_len = rec.sequence().len();
    let mut out = Vec::new();
    render_fastq_window(
        &mut out,
        &rec,
        &[],
        Window {
            start: 0,
            end: seq_len,
            idx: 0,
            total: 1,
        },
        inspect_mod_block(&rec, seq_len),
        None,
        platform(&rec),
        &cfg.fastq_tags,
        &cfg.remove_tags,
        Some(Reason::TagFilter),
    );
    Ok(RejectItem::Fastq(out))
}

/// Runs the single-threaded uBAM workflow: refuses aligned reads, trims, filters
/// each produced segment and writes the reconstructed survivors.
fn run_bam_seq(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>>,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    for rec in records {
        let rec = decode_raw_record(&rec?)?;
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .input_bases
            .fetch_add(rec.sequence().as_ref().len() as u64, Ordering::Relaxed);
        render_bam_read(header, &rec, cfg, counters, |out| {
            sink.write_record(header, out.as_ref().unwrap_or(&rec))
        })?;
    }
    Ok(counters.snapshot())
}

/// Runs `workflow::run_parallel` for BAM input: decodes each raw record on the
/// pool and hands the decoded record to `render`, which appends output items
/// to the batch buffer. The per-segment filter and counters are updated inside
/// `render` by `process_read_segments`.
fn run_bam_parallel<T, P, S, Render, Pack, WriteOne>(
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    cfg: &Config,
    sink: &mut S,
    render: Render,
    pack: Pack,
    write_one: WriteOne,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats>
where
    T: Send,
    P: Send,
    S: Send,
    Render: Fn(&bam::Record, &RecordBuf, &Config, &mut Vec<T>) -> anyhow::Result<()> + Sync,
    Pack: Fn(Vec<T>) -> std::io::Result<P> + Sync,
    WriteOne: Fn(&mut S, &P) -> std::io::Result<()> + Send,
{
    run_parallel(
        records,
        BAM_BATCH,
        |record: &bam::Record| record.sequence().len(),
        cfg,
        sink,
        |rec, cfg, out| render(&rec, &decode_raw_record(&rec)?, cfg, out),
        pack,
        write_one,
        counters,
    )
}

/// Compresses a batch of output records into BGZF blocks at the sink's level.
fn pack_bam_blocks(
    header: &sam::Header,
    level: u8,
    records: Vec<BamOutputRecord>,
) -> std::io::Result<Vec<u8>> {
    use noodles_sam::alignment::io::Write as _;
    let mut w = bam::io::Writer::from(Vec::new());
    for rec in &records {
        match rec {
            BamOutputRecord::Raw(record) => w.write_record(header, record)?,
            BamOutputRecord::Decoded(record) => w.write_alignment_record(header, record)?,
        }
    }
    let mut blocks = Vec::new();
    crate::io::bgzf::encode(level, &w.into_inner(), &mut blocks)?;
    Ok(blocks)
}

/// A record ready to write: the untouched raw input or a rebuilt decoded record.
enum BamOutputRecord {
    /// The raw input record, written without decoding.
    Raw(bam::Record),
    /// A rebuilt record.
    Decoded(RecordBuf),
}

fn raw_array_len(value: &noodles_sam::alignment::record::data::field::Value<'_>) -> Option<usize> {
    use noodles_sam::alignment::record::data::field::Value as RawValue;
    use noodles_sam::alignment::record::data::field::value::Array as RawArray;

    match value {
        RawValue::Array(RawArray::Int8(v)) => Some(v.len()),
        RawValue::Array(RawArray::UInt8(v)) => Some(v.len()),
        RawValue::Array(RawArray::Int16(v)) => Some(v.len()),
        RawValue::Array(RawArray::UInt16(v)) => Some(v.len()),
        RawValue::Array(RawArray::Int32(v)) => Some(v.len()),
        RawValue::Array(RawArray::UInt32(v)) => Some(v.len()),
        RawValue::Array(RawArray::Float(v)) => Some(v.len()),
        _ => None,
    }
}

/// The number of bases a raw `sa` coverage array covers; the borrowed
/// counterpart of `rle_coverage_len`.
fn raw_rle_coverage_len(
    value: &noodles_sam::alignment::record::data::field::Value<'_>,
) -> Option<usize> {
    use noodles_sam::alignment::record::data::field::Value as RawValue;
    use noodles_sam::alignment::record::data::field::value::Array as RawArray;

    fn collect<N: Into<i64>>(values: impl Iterator<Item = io::Result<N>>) -> Option<Vec<i64>> {
        values.map(|n| n.ok().map(Into::into)).collect()
    }
    let runs = match value {
        RawValue::Array(RawArray::Int8(v)) => collect(v.iter())?,
        RawValue::Array(RawArray::UInt8(v)) => collect(v.iter())?,
        RawValue::Array(RawArray::Int16(v)) => collect(v.iter())?,
        RawValue::Array(RawArray::UInt16(v)) => collect(v.iter())?,
        RawValue::Array(RawArray::Int32(v)) => collect(v.iter())?,
        RawValue::Array(RawArray::UInt32(v)) => collect(v.iter())?,
        _ => return None,
    };
    rle_runs_len(&runs)
}

/// The integer a raw aux value holds, whatever width it was stored at; the
/// borrowed counterpart of `aux_integer`.
fn raw_integer(value: &noodles_sam::alignment::record::data::field::Value<'_>) -> Option<i64> {
    use noodles_sam::alignment::record::data::field::Value as RawValue;

    Some(match value {
        RawValue::UInt8(n) => i64::from(*n),
        RawValue::Int8(n) => i64::from(*n),
        RawValue::UInt16(n) => i64::from(*n),
        RawValue::Int16(n) => i64::from(*n),
        RawValue::UInt32(n) => i64::from(*n),
        RawValue::Int32(n) => i64::from(*n),
        _ => return None,
    })
}

/// Inspects only the aux metadata that can change or affect advisories on an
/// otherwise full-window record. Returns the modification block's state and
/// whether a known per-base tag is malformed, without allocating owned tag
/// values.
fn raw_full_window_metadata(record: &bam::Record) -> std::io::Result<(ModBlock, bool)> {
    use noodles_sam::alignment::record::data::field::Value as RawValue;
    use noodles_sam::alignment::record::data::field::value::Array as RawArray;

    let seq_len = record.sequence().len();
    let data = record.data();
    let mut mm: Option<&[u8]> = None;
    let mut ml: Option<Option<usize>> = None;
    let mut mn: Option<Option<i64>> = None;
    let mut malformed_perbase = false;

    for result in data.iter() {
        let (tag, value) = result?;
        if tag == Tag::BASE_MODIFICATIONS {
            if let RawValue::String(s) = &value {
                mm = Some(AsRef::<[u8]>::as_ref(*s));
            }
        } else if tag == Tag::BASE_MODIFICATION_PROBABILITIES {
            ml = Some(match &value {
                RawValue::Array(RawArray::UInt8(v)) => Some(v.len()),
                _ => None,
            });
        } else if tag == Tag::BASE_MODIFICATION_SEQUENCE_LENGTH {
            mn = Some(raw_integer(&value));
        }

        let tag_bytes = <[u8; 2]>::from(tag);
        if KNOWN_PERBASE_TAGS.contains(&tag_bytes)
            && raw_array_len(&value).is_some_and(|len| len != seq_len)
        {
            malformed_perbase = true;
        }
        if tag_bytes == RLE_COVERAGE_TAG
            && raw_array_len(&value).is_some()
            && raw_rle_coverage_len(&value) != Some(seq_len)
        {
            malformed_perbase = true;
        }
    }

    let block = match mm {
        None => ModBlock::Absent,
        Some(mm) => {
            let block = classify_mod_block(mm, ml, mn, seq_len);
            if block != ModBlock::Malformed
                && !mods::parse::positions_valid(mm, record.sequence().iter())
            {
                ModBlock::Malformed
            } else {
                block
            }
        },
    };
    Ok((block, malformed_perbase))
}

/// Applies the aligned, reverse-complement and legacy-tag guards to a raw
/// record; the counterpart of `io::bam::ensure_trimmable`.
fn ensure_raw_trimmable(record: &bam::Record) -> anyhow::Result<()> {
    let legacy_tag = crate::io::bam::LEGACY_MOD_TAGS
        .into_iter()
        .find(|t| record.data().get(&Tag::new(t[0], t[1])).is_some());
    crate::io::bam::refuse_untrimmable(record.flags(), legacy_tag, || {
        crate::io::bam::display_name(record.name().map(AsRef::as_ref))
    })
}

/// Returns the GC fraction of a raw record's sequence, counted over its packed
/// bases without decoding them into a buffer, by the rule of
/// `filter::gc_fraction`.
fn raw_gc_fraction(record: &bam::Record) -> f64 {
    let sequence = record.sequence();
    if sequence.is_empty() {
        return 0.0;
    }
    let gc = sequence.iter().filter(|&b| crate::filter::is_gc(b)).count();
    gc as f64 / sequence.len() as f64
}

/// Filters one raw record over its full window and decides its output: the raw
/// record itself when nothing changes, a decoded rebuild when `MN` is missing
/// or the modification block is malformed, nothing when the filter drops it.
/// A malformed modification block or per-base tag is counted.
///
/// Tag removal never reaches here: `run_raw_bam` excludes it from the
/// full-window shortcut, so a run that removes tags rebuilds every record
/// through `run_bam`.
fn process_raw_full_window(
    record: bam::Record,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Option<BamOutputRecord>> {
    let seq_len = record.sequence().len();
    let qualities = record.quality_scores();
    let qual = qualities.as_ref();
    if qual.len() != seq_len || crate::io::bam::quality_absent(qual) {
        anyhow::bail!(
            "read {}: BAM record SEQ length {} != QUAL length {} \
             (records without per-base quality are not supported)",
            crate::io::bam::display_name(record.name().map(AsRef::as_ref)),
            seq_len,
            qual.len()
        );
    }

    let (mod_block, malformed_perbase) = raw_full_window_metadata(&record)?;
    if mod_block == ModBlock::Malformed {
        counters.malformed_mod_reads.fetch_add(1, Ordering::Relaxed);
    }
    if malformed_perbase {
        counters.malformed_tag_reads.fetch_add(1, Ordering::Relaxed);
    }

    let rejected = if seq_len == 0 {
        counters
            .reads_trimmed_to_nothing
            .fetch_add(1, Ordering::Relaxed);
        Some(Reason::TrimmedToNothing)
    } else {
        match crate::filter::check_with_gc(seq_len, qual, || raw_gc_fraction(&record), &cfg.filter)
        {
            Some(reason) => {
                counters.record_segment_drop(reason);
                counters.reads_all_filtered.fetch_add(1, Ordering::Relaxed);
                Some(Reason::Dropped(reason))
            },
            None => None,
        }
    };
    if let Some(reason) = rejected {
        if counters.wants_rejects() {
            let mut rec = decode_raw_record(&record)?;
            reject::tag_record(&mut rec, reason);
            counters.reject(RejectItem::Bam(rec))?;
        }
        return Ok(None);
    }

    // The window spans the whole record, so a survivor's output is its input.
    counters.output_reads.fetch_add(1, Ordering::Relaxed);
    counters
        .output_bases
        .fetch_add(seq_len as u64, Ordering::Relaxed);
    counters.reads_with_output.fetch_add(1, Ordering::Relaxed);

    let output = match mod_block {
        ModBlock::Absent | ModBlock::Consistent => BamOutputRecord::Raw(record),
        ModBlock::MissingMn | ModBlock::Malformed => {
            let decoded = decode_raw_record(&record)?;
            let window = Window {
                start: 0,
                end: seq_len,
                idx: 0,
                total: 1,
            };
            match reconstruct_window_record(
                &decoded,
                window,
                mod_block,
                None,
                None,
                &cfg.remove_tags,
            ) {
                Some(rebuilt) => BamOutputRecord::Decoded(rebuilt),
                None => BamOutputRecord::Raw(record),
            }
        },
    };
    Ok(Some(output))
}

fn run_raw_bam_full_window_seq(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>>,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    for record in records {
        let record = record?;
        ensure_raw_trimmable(&record)?;
        let seq_len = record.sequence().len();
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .input_bases
            .fetch_add(seq_len as u64, Ordering::Relaxed);
        match process_raw_full_window(record, cfg, counters)? {
            Some(BamOutputRecord::Raw(record)) => sink.write_raw_record(header, &record)?,
            Some(BamOutputRecord::Decoded(record)) => sink.write_record(header, &record)?,
            None => {},
        }
    }
    Ok(counters.snapshot())
}

fn run_raw_bam_full_window_parallel(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    let level = sink
        .block_level()
        .expect("A parallel run writes through a block sink");
    run_parallel(
        records,
        BAM_BATCH,
        |record: &bam::Record| record.sequence().len(),
        cfg,
        sink,
        |record, cfg, out: &mut Vec<BamOutputRecord>| {
            ensure_raw_trimmable(&record)?;
            out.extend(process_raw_full_window(record, cfg, counters)?);
            Ok(())
        },
        |records| pack_bam_blocks(header, level, records),
        |sink, blocks: &Vec<u8>| sink.write_blocks(blocks),
        counters,
    )
}

/// Runs the uBAM workflow on raw records from a production reader. Full-window
/// runs filter and write unchanged records without building an owned
/// `RecordBuf`; any configuration that can alter sequence or tags is routed to
/// `run_bam`, tag removal included, since a record that would otherwise pass
/// through untouched still has to be rebuilt without the removed tags.
pub fn run_raw_bam(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    let full_window = cfg.trim.head == 0
        && cfg.trim.tail == 0
        && cfg.trim.quality.is_none()
        && cfg.adapters.is_none()
        && cfg.remove_tags.is_empty();
    if !full_window {
        return run_bam(header, records, sink, cfg, counters);
    }
    if cfg.threads <= 1 {
        run_raw_bam_full_window_seq(header, records, sink, cfg, counters)
    } else {
        run_raw_bam_full_window_parallel(header, records, sink, cfg, counters)
    }
}

/// Runs the uBAM workflow: decodes, refuses aligned reads, trims, filters and
/// reconstructs. Sequential for `cfg.threads <= 1`; otherwise renders on a
/// rayon pool and drains the `RecordBuf`s through `run_bam_parallel`'s bounded
/// channel to the writer, in input order under `cfg.ordered` and in completion
/// order otherwise.
pub(crate) fn run_bam(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    if cfg.threads <= 1 {
        return run_bam_seq(header, records, sink, cfg, counters);
    }
    let level = sink
        .block_level()
        .expect("A parallel run writes through a block sink");
    run_bam_parallel(
        records,
        cfg,
        sink,
        // Render: the survivors of one record. An untouched record is
        // written from its raw input without re-encoding.
        |raw, rec, cfg, items| {
            render_bam_read(header, rec, cfg, counters, |out| {
                items.push(match out {
                    Some(out) => BamOutputRecord::Decoded(out),
                    None => BamOutputRecord::Raw(raw.clone()),
                });
                Ok(())
            })
        },
        // Pack: encode and compress the batch on the pool.
        |records| pack_bam_blocks(header, level, records),
        // Write: the compressed blocks, on the writer thread.
        |sink, blocks: &Vec<u8>| sink.write_blocks(blocks),
        counters,
    )
}

/// Appends the TAB-prefixed aux-tag block for one window to `tags`: carried
/// non-mod tags in source order with the `window_tag_updates` rewrites applied
/// in place and the removed ones skipped, per-base arrays sliced, then the
/// rebuilt MM/ML/MN block, then the added tags. Nothing is appended when
/// nothing is carried (the record then has a plain header). A `Malformed`
/// block is omitted. A tag named by `remove` is left out of the header, after
/// the rewrite, exactly as on BAM output.
#[allow(clippy::too_many_arguments)]
fn push_fastq_tags(
    tags: &mut Vec<u8>,
    src: &RecordBuf,
    seq: &[u8],
    window: Window,
    mod_block: ModBlock,
    indexed: Option<&IndexedMods>,
    sel: &FastqTags,
    platform: Platform,
    remove: &TagRemoval,
) {
    let Window { start, end, .. } = window;
    let orig_len = seq.len();
    let trimmed = start != 0 || end != orig_len;
    // BAM-to-FASTQ never rewrites the move table (a sliced one is impractical
    // in a FASTQ header, and signal-aware consumers read BAM), so a trim drops
    // the signal and poly-A tags.
    let mut updates =
        window_tag_updates(src, src.quality_scores().as_ref(), window, platform, None);
    if !remove.is_empty() {
        updates.retain(|(t, _)| !remove.contains(&<[u8; 2]>::from(*t)));
    }
    for (tag, value) in src.data().iter() {
        let t = <[u8; 2]>::from(tag);
        if matches!(&t, b"MM" | b"ML" | b"MN") {
            continue; // handled by the rebuilt block below
        }
        let rewritten = updates
            .iter()
            .position(|(u, _)| *u == tag)
            .map(|i| updates.remove(i).1);
        if !sel.carries(&t) || remove.contains(&t) {
            continue;
        }
        let value: Cow<Value> = match rewritten {
            Some(None) => continue,
            Some(Some(v)) => Cow::Owned(v),
            None => match trimmed
                .then(|| windowed_value(t, value, orig_len, start, end))
                .flatten()
            {
                Some(v) => Cow::Owned(v),
                None => Cow::Borrowed(value),
            },
        };
        tags.push(b'\t');
        push_aux_field(tags, t, &value);
    }
    if sel.carries_mods()
        && matches!(mod_block, ModBlock::Consistent | ModBlock::MissingMn)
        && let Some((mm, ml)) = mod_tags(src)
    {
        let (mm, ml) = rebuild_mods(mm, ml, seq, start, end, indexed);
        push_mods_aux(tags, &mm, ml.as_deref(), end - start, remove);
    }
    for (tag, value) in updates {
        let t = <[u8; 2]>::from(tag);
        if let Some(v) = value
            && sel.carries(&t)
        {
            tags.push(b'\t');
            push_aux_field(tags, t, &v);
        }
    }
}

/// Appends one surviving window of a decoded record to `out` as a FASTQ
/// record: the platform's segment name (`segment_name`), the selected aux
/// tags, then the sliced bases and qualities. The header is assembled in
/// place, so the tag text is formatted once, into the output buffer.
#[allow(clippy::too_many_arguments)]
fn render_fastq_window(
    out: &mut Vec<u8>,
    rec: &RecordBuf,
    description: &[u8],
    window: Window,
    mod_block: ModBlock,
    indexed: Option<&IndexedMods>,
    platform: Platform,
    sel: &FastqTags,
    remove: &TagRemoval,
    reason: Option<Reason>,
) {
    let Window { start, end, .. } = window;
    let seq = rec.sequence().as_ref();
    let qual = rec.quality_scores().as_ref();
    let coords = query_span(rec).map(|(qs0, _)| window_coords(qs0, start, end));
    let name = rec.name().map(|n| n.as_ref()).unwrap_or_default();
    out.push(b'@');
    if window.total > 1 || start != 0 || end != seq.len() {
        out.extend_from_slice(&segment_name(platform, name, window, coords));
    } else {
        out.extend_from_slice(name);
    }
    out.extend_from_slice(description);
    push_fastq_tags(
        out, rec, seq, window, mod_block, indexed, sel, platform, remove,
    );
    if let Some(reason) = reason {
        reject::push_fastq_tag(out, reason);
    }
    push_record_body(out, &seq[start..end], &qual[start..end]);
}

/// Renders one decoded record for FASTQ output: every surviving window is
/// appended to `buf`. The buffer is shared across records within a batch.
fn render_bam_fastq_read(
    rec: &RecordBuf,
    cfg: &Config,
    counters: &Counters,
    buf: &mut Vec<u8>,
    description: &[u8],
) -> anyhow::Result<()> {
    let platform = platform(rec);
    render_windows(rec, cfg, counters, |window, mod_block, indexed, reason| {
        let mut rejected = Vec::new();
        let out = if reason.is_some() {
            &mut rejected
        } else {
            &mut *buf
        };
        render_fastq_window(
            out,
            rec,
            description,
            window,
            mod_block,
            indexed,
            platform,
            &cfg.fastq_tags,
            &cfg.remove_tags,
            reason,
        );
        if reason.is_some() {
            counters.reject(RejectItem::Fastq(rejected))?;
        }
        Ok(())
    })
}

/// Decodes and renders BAM records as FASTQ, reusing a buffer per batch.
pub(crate) fn run_bam_to_fastq<W: BatchSink>(
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    writer: &mut W,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    if cfg.threads <= 1 {
        let mut buf = Vec::new();
        for rec in records {
            let rec = decode_raw_record(&rec?)?;
            counters.input_reads.fetch_add(1, Ordering::Relaxed);
            counters
                .input_bases
                .fetch_add(rec.sequence().len() as u64, Ordering::Relaxed);
            buf.clear();
            render_bam_fastq_read(&rec, cfg, counters, &mut buf, &[])?;
            writer.write_all(&buf)?;
        }
        return Ok(counters.snapshot());
    }
    run_bytes_parallel(
        records,
        BAM_BATCH,
        |rec| rec.sequence().len(),
        cfg,
        writer,
        |rec, cfg, buf| render_bam_fastq_read(&decode_raw_record(&rec)?, cfg, counters, buf, &[]),
        counters,
    )
}

/// Appends a tagged FASTQ read's output while preserving its header description.
pub(crate) fn render_tagged_fastq_read(
    rec: crate::record::ReadRecord,
    cfg: &Config,
    counters: &Counters,
    buf: &mut Vec<u8>,
) -> anyhow::Result<()> {
    let head_end = rec
        .name
        .iter()
        .position(|&b| b == b'\t')
        .unwrap_or(rec.name.len());
    let id_end = rec.name[..head_end]
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(head_end);
    let description = rec.name[id_end..head_end].to_vec();
    let rec = crate::io::tagged::record_from_tagged(rec)?;
    render_bam_fastq_read(&rec, cfg, counters, buf, &description)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod barcode_tests;
