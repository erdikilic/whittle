//! The barcode interval of a record, from its `bi` tag and barcode call, which
//! bounds adapter trimming.

use super::*;

/// Dorado's `bi` barcode-info tag: a `B:f` array of exactly seven floats,
/// `[barcode_score, front_start_index, front_len, front_score, rear_end_index,
/// rear_len, rear_score]` (`read_pipeline/base/messages.cpp`).
pub(super) const BARCODE_TAG: [u8; 2] = *b"bi";

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
pub(super) fn barcode_position(value: f32) -> Option<i64> {
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
pub(super) fn barcode_interval(
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
pub(super) fn barcode_call(rec: &RecordBuf) -> Option<&[u8]> {
    match rec.data().get(&Tag::new(b'B', b'C'))? {
        Value::String(value) => Some(value.as_ref()),
        _ => None,
    }
}

/// Resolves the retained window from a record's verified barcode spans: a
/// span is trimmed only when a barcode sequence is found at it. Returns the
/// window and whether any recorded span failed verification.
pub(super) fn verified_barcode_window(
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
