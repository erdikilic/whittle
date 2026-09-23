//! Platform detection and the names of output segments: parent read
//! identifiers, PacBio movie coordinates and split suffixes.

use super::*;

/// Returns the parent read id for a subread: the source's own `pi` if it has
/// one (so `pi` always names the ultimate ancestor, matching dorado), else the
/// source read name.
pub(super) fn parent_read_id(src: &RecordBuf) -> Vec<u8> {
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
pub(super) struct PacBioName<'a> {
    /// The name without any `/{qStart}_{qEnd}` interval: `{movie}/{zmw}` for
    /// a subread, `{movie}/{zmw}/ccs` for a CCS read, with `/fwd` or `/rev`
    /// for a by-strand read.
    pub(super) stem: &'a [u8],
    /// The `qStart` of a `{qStart}_{qEnd}` interval in the name.
    pub(super) query_start: Option<i64>,
}

/// Parses a read name against the PacBio BAM conventions: `{movie}/{zmw}/ccs`
/// with an optional `/fwd` or `/rev` and an optional `/{qStart}_{qEnd}`
/// (segmented reads), or the subread form `{movie}/{zmw}/{qStart}_{qEnd}`.
/// `None` for any other name.
pub(super) fn parse_pacbio_name(name: &[u8]) -> Option<PacBioName<'_>> {
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
pub(super) fn query_span(src: &RecordBuf) -> Option<(i64, i64)> {
    Some((signal_int(src, b"qs")?, signal_int(src, b"qe")?))
}

/// Returns the query coordinates of window `[start, end)` in the frame of the
/// original PacBio read whose query starts at `qs0`: the PacBio BAM spec keeps
/// `qs`/`qe` with respect to the original read through clipping.
pub(super) fn window_coords(qs0: i64, start: usize, end: usize) -> (i64, i64) {
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
