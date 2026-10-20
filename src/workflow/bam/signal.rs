//! ONT signal tags of a trimmed record: the move table, the sample counts and
//! split linkage, poly(A) coordinates and read timestamps.

use super::*;

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

/// Signal direction named by a dorado read-group id,
/// `<run_id>_<basecall_model>[_<barcode>]`: true for an RNA model (`rna0...`),
/// false for a DNA model (`dna_r...`), found at the start of the id or after an
/// underscore.
fn group_id_reversed(id: &[u8]) -> Option<bool> {
    (0..id.len())
        .filter(|&i| i == 0 || id[i - 1] == b'_')
        .find_map(|i| {
            let model = &id[i..];
            if model.starts_with(b"rna0") {
                Some(true)
            } else if model.starts_with(b"dna_r") {
                Some(false)
            } else {
                None
            }
        })
}

/// Signal direction named by `basecall_model=` in a read group's description:
/// true for an RNA model, false for a DNA model.
fn described_reversed(
    group: &sam::header::record::value::Map<sam::header::record::value::map::ReadGroup>,
) -> Option<bool> {
    use sam::header::record::value::map::read_group::tag::DESCRIPTION;
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
}

/// Signal direction resolved from the record's basecalling model: the
/// `basecall_model` in the description of its `@RG` line, or else the model a
/// dorado read-group id names, which also covers an id the header lacks.
/// Without `header` (tagged FASTQ input), only the id in the record's `RG` tag
/// is read. A BAM record without `RG` takes the direction every read group of
/// the header shares.
pub(super) fn signal_reversed(header: Option<&sam::Header>, rec: &RecordBuf) -> Option<bool> {
    let direction = |id: &[u8], group| described_reversed(group).or_else(|| group_id_reversed(id));
    match (rec.data().get(&Tag::READ_GROUP), header) {
        (Some(Value::String(id)), Some(header)) => {
            let id: &[u8] = id.as_ref();
            match header.read_groups().get(id) {
                Some(group) => direction(id, group),
                None => group_id_reversed(id),
            }
        },
        (Some(Value::String(id)), None) => group_id_reversed(id.as_ref()),
        (None, Some(header)) => {
            let mut groups = header.read_groups().iter();
            let (id, group) = groups.next()?;
            let first = direction(id.as_ref(), group)?;
            groups
                .all(|(id, group)| direction(id.as_ref(), group) == Some(first))
                .then_some(first)
        },
        _ => None,
    }
}

/// Emitted-base boundaries in signal block order, indexed once per read.
pub(super) struct MoveIndex<'a> {
    pub(super) stride: i8,
    pub(super) moves: &'a [i8],
    pub(super) boundaries: Vec<usize>,
    pub(super) reversed: bool,
}

impl<'a> MoveIndex<'a> {
    /// Accepts binary move tables with one emission per sequence base.
    pub(super) fn new(src: &'a RecordBuf, reversed: bool) -> Option<Self> {
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

/// The move table of one record under `--update-moves`, indexed on the first
/// trimmed window that needs it. Shared by the BAM and FASTQ output paths, so
/// both rewrite the signal tags of a window alike.
pub(super) struct ReadMoves<'a> {
    /// The record whose move table is indexed.
    rec: &'a RecordBuf,
    /// Whether the signal tags are rewritten (`--update-moves`).
    update: bool,
    /// The signal direction (`signal_reversed`), `None` when the model is
    /// unknown or `update` is off.
    direction: Option<bool>,
    /// The move index, built on first use; `None` inside for a missing or
    /// malformed table.
    index: std::cell::OnceCell<Option<MoveIndex<'a>>>,
}

impl<'a> ReadMoves<'a> {
    /// Resolves the record's signal direction when `update` is set; `header`
    /// is `None` for tagged FASTQ input.
    pub(super) fn new(rec: &'a RecordBuf, header: Option<&sam::Header>, update: bool) -> Self {
        Self {
            rec,
            update,
            direction: update.then(|| signal_reversed(header, rec)).flatten(),
            index: std::cell::OnceCell::new(),
        }
    }

    /// Returns the move index for `window`: `None` for a full window, without
    /// `--update-moves`, or for a missing or malformed table, which leaves the
    /// signal tags to be removed. A trimmed window of a read with a move table
    /// and no resolvable model is an error naming the read.
    pub(super) fn for_window(&self, window: Window) -> anyhow::Result<Option<&MoveIndex<'a>>> {
        let partial = window.start != 0 || window.end != self.rec.sequence().len();
        if !self.update || !partial {
            return Ok(None);
        }
        let Some(reverse) = self.direction else {
            if self.rec.data().get(&Tag::new(b'm', b'v')).is_some() {
                anyhow::bail!(
                    "read {}: --update-moves requires a DNA or RNA basecall model, named by \
                     basecall_model in the @RG description or by a dorado read-group id in RG",
                    crate::io::bam::display_name(self.rec.name().map(AsRef::as_ref))
                );
            }
            return Ok(None);
        };
        Ok(self
            .index
            .get_or_init(|| MoveIndex::new(self.rec, reverse))
            .as_ref())
    }
}

/// Reads an integer aux tag as `i64`, regardless of stored width.
pub(super) fn signal_int(src: &RecordBuf, tag: &[u8; 2]) -> Option<i64> {
    src.data()
        .get(&Tag::new(tag[0], tag[1]))
        .and_then(aux_integer)
}

pub(super) fn signal_offset(blocks: usize, stride: usize) -> Option<i64> {
    blocks
        .checked_mul(stride)
        .and_then(|n| i64::try_from(n).ok())
}

pub(super) fn signal_int_value(n: i64) -> Option<Value> {
    if let Ok(n) = i32::try_from(n) {
        Some(Value::Int32(n))
    } else if let Ok(n) = u32::try_from(n) {
        Some(Value::UInt32(n))
    } else {
        None
    }
}

/// Computes the poly-A tag updates (`pa` signal boundaries, `pt` tail length)
/// for a trimmed read. `pa` holds absolute original-signal positions, the frame
/// `ts` and `ns` use; `-1`/`-2` are dorado's not-found/not-enabled sentinels and
/// are left as is. When every real position falls inside
/// `[kept_start, kept_end)` the tail survived: a split shifts `pa` into the
/// subread's own signal frame, a crop keeps both unchanged. Otherwise, or with
/// no poly-A array, `pa`/`pt` are dropped.
pub(super) fn polya_updates(
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
pub(super) struct SignalWindow {
    /// First kept sample, inclusive.
    pub(super) kept_start: i64,
    /// End of the kept signal, exclusive.
    pub(super) kept_end: i64,
}

/// Computes the ONT signal tag updates for output window `[start, end)`.
/// Returns `(tag, Some(value))` to set or `(tag, None)` to remove, with the
/// kept signal window when the tags were rewritten; empty when the read is not
/// trimmed. With `update_moves` off, or a missing or malformed move table, the
/// signal tags and both poly-A tags are removed; a crop keeps `sp` and `pi`,
/// which place the read's raw signal in its parent and stay valid when only
/// bases are removed. With it on, `mv` is
/// sliced by block range (stride-aligned, following dorado
/// `splitter::subread`) and:
///   - crop (`total == 1`, name kept): `ts += block_first*stride`; `ns` is the
///     kept signal's end, which is the source `ns` when the window runs to the
///     last base.
///   - split (`total > 1`, renamed): `ts = 0`, `ns = span`,
///     `sp = parent offset`, `pi = parent id`.
pub(super) fn signal_tag_updates(
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
        .filter(|t| total > 1 || !matches!(*t, b"sp" | b"pi"))
        .chain(POLYA_TAGS.iter())
        .map(|t| (Tag::new(t[0], t[1]), None))
        .collect();
    (dropped, None)
}

/// Rewrites the signal and poly-A tags of a trimmed window from the move
/// table. `None` when the table is missing or malformed, its base count
/// disagrees with the sequence, the window has no start base, or a signal
/// offset does not fit its tag; the caller then removes the tags.
pub(super) fn signal_rewrite(
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
pub(super) fn split_time_updates(src: &RecordBuf, window: SignalWindow) -> TagUpdates {
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

/// Rescales `du` for a crop whose signal window is known, so that `ns` over
/// `du` stays the sample rate: a tail crop shortens `ns`, and dorado's own
/// trimming writes `du` as `ns` over the sample rate. Empty when the crop keeps
/// the source `ns`, or `du` or `ns` is missing or not positive.
pub(super) fn crop_time_updates(src: &RecordBuf, window: SignalWindow) -> TagUpdates {
    let du_tag = Tag::new(b'd', b'u');
    let (Some(Value::Float(du0)), Some(ns0)) = (src.data().get(&du_tag), signal_int(src, b"ns"))
    else {
        return Vec::new();
    };
    if ns0 <= 0 || !du0.is_finite() || *du0 <= 0.0 || window.kept_end == ns0 {
        return Vec::new();
    }
    let du = f64::from(*du0) * window.kept_end as f64 / ns0 as f64;
    vec![(du_tag, Some(Value::Float(du as f32)))]
}

/// Returns an RFC 3339 `st` value advanced by `seconds`, at millisecond
/// precision (dorado's own) and in the offset form the source uses: `Z`, a
/// numeric offset, or none (a civil time read as UTC). `None` when the source
/// does not parse.
pub(super) fn shift_timestamp(value: &[u8], seconds: f64) -> Option<Vec<u8>> {
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
