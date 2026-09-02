//! uBAM workflows: record reconstruction (sequence, quality, MM/ML/MN, per-base and signal tags) and the sequential, parallel and raw full-window drivers for BAM and FASTQ output.

use std::borrow::Cow;
use std::io::{self, Write};
use std::sync::Arc;
use std::sync::atomic::Ordering;

use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::{Op, op::Kind};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};

use super::{BAM_BATCH, Counters, Rendered, Stats, process_read_segments, run_parallel};
use crate::config::{Config, FastqTags};
use crate::io::fastq::{format_aux_field, format_mods_aux, write_named_record};
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
pub(crate) const SIGNAL_TAGS: [[u8; 2]; 5] = [*b"mv", *b"ts", *b"ns", *b"sp", *b"pi"];

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
/// sequence positions that shift under a crop and cannot be reconstructed from
/// the BAM, and the PacBio undo blobs `ds`/`ls` describe the untrimmed read.
/// The barcode call itself (`BC`/`bv`) is a per-read label and is copied
/// unchanged.
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

/// Input accepted by the BAM workflows. Production readers yield lazy raw
/// `bam::Record`s, so structured field/tag decoding happens on a render worker;
/// tests and library callers can still provide an already-decoded `RecordBuf`.
pub trait InputRecord: Send {
    /// Length of the record's sequence, available without decoding.
    fn sequence_len(&self) -> usize;
    /// Decodes into an owned `RecordBuf`.
    fn decode(self) -> std::io::Result<RecordBuf>;
}

impl InputRecord for bam::Record {
    fn sequence_len(&self) -> usize {
        self.sequence().len()
    }

    fn decode(self) -> std::io::Result<RecordBuf> {
        decode_raw_record(&self)
    }
}

/// Converts a raw record to a `RecordBuf` on the render worker without routing
/// sequence, quality and every aux value through the generic SAM trait
/// iterators. The concrete noodles views have bulk conversions for these large
/// fields and reduce conversion overhead on long reads.
fn decode_raw_record(src: &bam::Record) -> std::io::Result<RecordBuf> {
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
    resolve_long_cigar(&mut dst)?;
    Ok(dst)
}

/// Expands BAM's `CG:B:I` overflow representation for records whose CIGAR does
/// not fit in the 16-bit operation count. This mirrors the resolution performed
/// by `noodles_bam::io::Reader::read_record_buf` after decoding a buffered
/// record, which the lazy raw-record API leaves to its caller.
fn resolve_long_cigar(record: &mut RecordBuf) -> io::Result<()> {
    let is_overflow_placeholder = match record.cigar().as_ref() {
        [op_0, op_1] => {
            *op_0 == Op::new(Kind::SoftClip, record.sequence().len()) && op_1.kind() == Kind::Skip
        },
        _ => false,
    };

    if !is_overflow_placeholder {
        return Ok(());
    }

    let Some((_, value)) = record.data_mut().remove(&Tag::CIGAR) else {
        return Ok(());
    };
    let Value::Array(Array::UInt32(values)) = value else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CG data field type",
        ));
    };

    let cigar = record.cigar_mut().as_mut();
    cigar.clear();
    for n in values {
        let kind = match n & 0x0f {
            0 => Kind::Match,
            1 => Kind::Insertion,
            2 => Kind::Deletion,
            3 => Kind::Skip,
            4 => Kind::SoftClip,
            5 => Kind::HardClip,
            6 => Kind::Pad,
            7 => Kind::SequenceMatch,
            8 => Kind::SequenceMismatch,
            actual => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid CG CIGAR operation kind: {actual}"),
                ));
            },
        };
        cigar.push(Op::new(kind, (n >> 4) as usize));
    }

    Ok(())
}

impl InputRecord for RecordBuf {
    fn sequence_len(&self) -> usize {
        self.sequence().as_ref().len()
    }

    fn decode(self) -> std::io::Result<RecordBuf> {
        Ok(self)
    }
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

/// The state of a record's `MM`/`ML`/`MN` block relative to its sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModBlock {
    /// No `MM:Z` tag. An `ML` or `MN` present on its own is copied verbatim.
    Absent,
    /// `MM` parses to its end, `ML` (when present) is a `B:C` array of the
    /// length `MM` declares, and `MN` equals the sequence length.
    Consistent,
    /// As `Consistent`, with `MN` absent; the output gains one.
    MissingMn,
    /// `MM` does not parse to its end, `ML` is not a `B:C` array or has the
    /// wrong length, or `MN` disagrees with the sequence length. The calls
    /// cannot be placed on the sequence, so the block is removed from the
    /// output and the read is counted in `Counters::malformed_mod_reads`.
    Malformed,
}

/// Classifies a modification block from its parts. `ml` is `None` when the tag
/// is absent and `Some(None)` when it is present with a subtype other than
/// `B:C`; `mn` is `None` when absent and `Some(None)` when not an integer.
fn classify_mod_block(
    mm: &[u8],
    ml: Option<Option<usize>>,
    mn: Option<Option<i64>>,
    seq_len: usize,
) -> ModBlock {
    let Ok(expected) = mods::expected_ml_len(mm) else {
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
pub fn inspect_mod_block(src: &RecordBuf, seq_len: usize) -> ModBlock {
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
    classify_mod_block(mm, ml, mn, seq_len)
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
pub enum Platform {
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
pub fn platform(rec: &RecordBuf) -> Platform {
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

/// Returns the output name of `window` for a read named `name`. A crop keeps
/// the name. A split names an ONT segment `{name}_segment_{n}` (1-based) and a
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
    if window.total <= 1 {
        return name.to_vec();
    }
    if platform == Platform::PacBio
        && let Some(parts) = parse_pacbio_name(name)
    {
        let (qs, qe) = coords.unwrap_or_else(|| {
            window_coords(parts.query_start.unwrap_or(0), window.start, window.end)
        });
        let mut out = parts.stem.to_vec();
        out.extend_from_slice(format!("/{qs}_{qe}").as_bytes());
        return out;
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
    update_moves: bool,
) -> (TagUpdates, Option<SignalWindow>) {
    if start == 0 && end == seq_len {
        return (Vec::new(), None);
    }
    if update_moves && let Some((updates, window)) = signal_rewrite(src, seq_len, start, end, total)
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
) -> Option<(TagUpdates, SignalWindow)> {
    let (stride, moves) = src
        .data()
        .get(&Tag::new(b'm', b'v'))
        .and_then(parse_move_table)?;
    // The move index of the `start`-th and `end`-th base (each `1` in `moves`
    // is one emitted base), and the total base count, found in a single pass,
    // without materializing the whole positions list.
    let mut ones_seen = 0usize;
    let mut block_first = None;
    let mut block_second = None;
    for (i, &m) in moves.iter().enumerate() {
        if m != 0 {
            if ones_seen == start {
                block_first = Some(i);
            }
            if ones_seen == end {
                block_second = Some(i);
            }
            ones_seen += 1;
        }
    }
    if ones_seen != seq_len {
        return None;
    }

    let stride_n = stride as usize;
    // The end-th base exists only when `end < seq_len`; otherwise the window
    // runs to the table end. An empty or out-of-range window has no start base.
    let block_first = block_first?;
    let block_second = if end == seq_len {
        moves.len()
    } else {
        block_second?
    };

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
pub fn has_malformed_perbase_tag(rec: &RecordBuf, seq_len: usize) -> bool {
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
    update_moves: bool,
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

    let (mut updates, signal) = signal_tag_updates(src, orig_len, start, end, total, update_moves);
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
        // every subread but the last, which ends where the read did.
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
        if idx + 1 < total && src.data().get(&er).is_some() {
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

/// Builds one output uBAM record for interval `[start, end)`, segment `idx` of
/// `total`: SEQ/QUAL sliced, `MM`/`ML`/`MN` rebuilt, per-base kinetics sliced,
/// stale signal-space tags rewritten or dropped, the name set per platform on a
/// split (`segment_name`).
/// Remaining aux tags are copied unchanged.
pub fn reconstruct_record(
    src: &RecordBuf,
    start: usize,
    end: usize,
    total: usize,
    idx: usize,
    update_moves: bool,
) -> RecordBuf {
    let seq = src.sequence().as_ref();
    let qual = src.quality_scores().as_ref();
    let mod_block = inspect_mod_block(src, seq.len());
    let window = Window {
        start,
        end,
        idx,
        total,
    };
    reconstruct_record_with_bases(src, seq, qual, window, mod_block, update_moves)
}

/// Builds one output record for `window`. The record is assembled field by
/// field: SEQ/QUAL are sliced, aux tags are copied in source order with the
/// rewritten ones replaced in place, removed ones skipped and added ones
/// appended. A `Malformed` block is removed. An untrimmed, unsplit record with
/// an `Absent` or `Consistent` block is cloned as is.
fn reconstruct_record_with_bases(
    src: &RecordBuf,
    seq: &[u8],
    qual: &[u8],
    window: Window,
    mod_block: ModBlock,
    update_moves: bool,
) -> RecordBuf {
    let Window {
        start, end, total, ..
    } = window;
    let orig_len = seq.len();
    debug_assert_eq!(src.sequence().as_ref(), seq);
    debug_assert_eq!(src.quality_scores().as_ref(), qual);
    let trimmed = start != 0 || end != orig_len;
    let split = total > 1;
    if !trimmed && !split && matches!(mod_block, ModBlock::Absent | ModBlock::Consistent) {
        return src.clone();
    }

    let platform = platform(src);
    let mut out = RecordBuf::default();
    *out.name_mut() = if split {
        let name = src.name().map(|n| n.to_vec()).unwrap_or_default();
        let coords = query_span(src).map(|(qs0, _)| window_coords(qs0, start, end));
        Some(segment_name(platform, &name, window, coords).into())
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
                let (mm, ml) = rebuild_mods(mm, ml, seq, start, end);
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
    updates.extend(window_tag_updates(
        src,
        qual,
        window,
        platform,
        update_moves,
    ));

    let data = out.data_mut();
    for (tag, value) in src.data().iter() {
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

    out
}

/// Slices a record's `MM`/`ML` to the window `[start, end)` and re-serializes
/// them. `None` when the record carries no `MM:Z` or its block is `Malformed`
/// (see `ModBlock`). The inner `Option` is the sliced `ML`, `None` when the
/// source has none: `ML` is optional per the SAM spec, so such a record must not
/// gain one. Shared by the BAM-to-BAM and BAM-to-FASTQ paths.
pub fn reconstruct_mods(
    src: &RecordBuf,
    seq: &[u8],
    start: usize,
    end: usize,
) -> Option<(Vec<u8>, Option<Vec<u8>>)> {
    match inspect_mod_block(src, seq.len()) {
        ModBlock::Consistent | ModBlock::MissingMn => {
            let (mm, ml) = mod_tags(src)?;
            Some(rebuild_mods(mm, ml, seq, start, end))
        },
        ModBlock::Absent | ModBlock::Malformed => None,
    }
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
) -> (Vec<u8>, Option<Vec<u8>>) {
    // Over the full window the rebuild is the identity, so the source bytes are
    // returned and the parse, slice and re-serialize are skipped; that work is
    // the dominant cost of an untrimmed BAM-to-FASTQ run.
    if start == 0 && end == seq.len() {
        let fast = (mm.to_vec(), ml.map(<[u8]>::to_vec));
        debug_assert_eq!(
            fast,
            rebuild_mods_windowed(mm, ml, seq, start, end),
            "The full-window shortcut must agree with the general path"
        );
        return fast;
    }
    rebuild_mods_windowed(mm, ml, seq, start, end)
}

/// Rebuilds the MM/ML block for any window the full-window shortcut does not cover.
fn rebuild_mods_windowed(
    mm: &[u8],
    ml: Option<&[u8]>,
    seq: &[u8],
    start: usize,
    end: usize,
) -> (Vec<u8>, Option<Vec<u8>>) {
    let parsed = mods::parse(mm, ml.unwrap_or(&[]));
    let sliced = mods::reconstruct(&parsed, seq, start, end);
    let (mm_new, ml_new) = mods::serialize(&sliced);
    (mm_new, ml.map(|_| ml_new))
}

/// Runs the per-read guards and bookkeeping shared by the decoded workflows:
/// refuses aligned reads and legacy mod tags, requires full per-base quality, and
/// classifies the modification block, counting a malformed one. Returns the
/// record's bases, qualities and block.
fn prepare_read<'a>(
    rec: &'a RecordBuf,
    counters: &Counters,
) -> anyhow::Result<(&'a [u8], &'a [u8], ModBlock)> {
    crate::io::bam::ensure_unaligned(rec)?;
    crate::io::bam::ensure_modern_mod_tags(rec)?;
    let seq = rec.sequence().as_ref();
    let qual = rec.quality_scores().as_ref();
    if qual.len() != seq.len() {
        let name = rec
            .name()
            .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
            .unwrap_or_else(|| "<unnamed>".to_string());
        anyhow::bail!(
            "read {name}: BAM record SEQ length {} != QUAL length {} \
             (records without full per-base quality are not supported)",
            seq.len(),
            qual.len()
        );
    }
    let mod_block = inspect_mod_block(rec, seq.len());
    if mod_block == ModBlock::Malformed {
        counters.malformed_mod_reads.fetch_add(1, Ordering::Relaxed);
    }
    Ok((seq, qual, mod_block))
}

/// Runs the single-threaded uBAM workflow: refuses aligned reads, trims, filters
/// each produced segment and reconstructs the survivors.
fn run_bam_seq<R: InputRecord>(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<R>>,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    let mut malformed_tag_reads = 0u64;
    for rec in records {
        let rec = rec?.decode()?;
        let (seq, qual, mod_block) = prepare_read(&rec, counters)?;
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .input_bases
            .fetch_add(seq.len() as u64, Ordering::Relaxed);
        if has_malformed_perbase_tag(&rec, seq.len()) {
            malformed_tag_reads += 1;
        }
        let _read =
            crate::workflow::read_span(rec.name().map(|n| n.as_ref()).unwrap_or(b"<unnamed>"));
        let _read = _read.enter();
        let produced = trim::apply(seq, qual, &cfg.trim, cfg.adapters.as_ref());
        let mut survivors: Vec<(usize, usize)> = Vec::new();
        process_read_segments(
            &produced,
            seq,
            qual,
            &cfg.filter,
            counters,
            |idx, total, s, e| {
                survivors.push((s, e));
                let window = Window {
                    start: s,
                    end: e,
                    idx,
                    total,
                };
                let out = reconstruct_record_with_bases(
                    &rec,
                    seq,
                    qual,
                    window,
                    mod_block,
                    cfg.update_moves,
                );
                sink.write_record(header, &out)?;
                Ok(())
            },
        )?;
        count_undo_tags_dropped(counters, &rec, seq.len(), &survivors);
    }
    Ok(counters.snapshot(malformed_tag_reads))
}

/// Runs `workflow::run_parallel` for BAM input: decodes each raw record on the
/// pool, notes a malformed per-base tag, and hands the decoded record to
/// `render`. `render` returns the surviving segments only; the per-segment
/// filter and counters are updated inside it by `process_read_segments`.
fn run_bam_parallel<R, T, S, Render, WriteOne>(
    records: impl Iterator<Item = anyhow::Result<R>> + Send,
    cfg: &Config,
    sink: &mut S,
    render: Render,
    write_one: WriteOne,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats>
where
    R: InputRecord + Send,
    T: Send,
    S: Send,
    Render: Fn(&RecordBuf, &Config) -> anyhow::Result<Vec<T>> + Sync,
    WriteOne: Fn(&mut S, &T) -> std::io::Result<()> + Send,
{
    run_parallel(
        records,
        BAM_BATCH,
        InputRecord::sequence_len,
        cfg,
        sink,
        |rec, cfg| {
            let rec = rec.decode()?;
            let seq_len = rec.sequence().as_ref().len();
            let malformed_tags = has_malformed_perbase_tag(&rec, seq_len);
            let items = render(&rec, cfg)?;
            Ok(Rendered {
                items,
                malformed_tags,
            })
        },
        write_one,
        counters,
    )
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
        Some(mm) => classify_mod_block(mm, ml, mn, seq_len),
    };
    Ok((block, malformed_perbase))
}

/// Applies the aligned, legacy-tag and reverse-complement guards to a raw
/// record; the counterpart of `io::bam::ensure_unaligned` and
/// `ensure_modern_mod_tags`.
fn ensure_raw_unaligned(record: &bam::Record) -> anyhow::Result<()> {
    let flags = record.flags();
    let name = || {
        record
            .name()
            .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
            .unwrap_or_else(|| "<unnamed>".to_string())
    };
    if !flags.is_unmapped() {
        anyhow::bail!(
            "read {} is aligned (mapped); only unaligned BAM (uBAM) input is supported",
            name()
        );
    }
    for lt in crate::io::bam::LEGACY_MOD_TAGS {
        if record.data().get(&Tag::new(lt[0], lt[1])).is_some() {
            anyhow::bail!(
                "read {} carries the legacy `{}` base-modification tag; whittle rewrites only \
                 the current `MM`/`ML` spelling, so trimming this record would leave its \
                 modification calls pointing at the wrong bases",
                name(),
                String::from_utf8_lossy(&lt)
            );
        }
    }
    // See `io::bam::ensure_unaligned` for why a reverse-complemented record is
    // refused rather than trimmed.
    if flags.is_reverse_complemented() {
        anyhow::bail!(
            "read {} is flagged reverse-complemented; whittle trims in read orientation and \
             cannot keep position-indexed tags correct for it",
            name()
        );
    }
    Ok(())
}

fn raw_gc_fraction(record: &bam::Record) -> f64 {
    let sequence = record.sequence();
    if sequence.is_empty() {
        return 0.0;
    }
    let gc = sequence
        .iter()
        .filter(|&b| matches!(b, b'G' | b'g' | b'C' | b'c'))
        .count();
    gc as f64 / sequence.len() as f64
}

/// Filters one raw record over its full window and decides its output: the raw
/// record itself when nothing changes, a decoded rebuild when `MN` is missing
/// or the modification block is malformed, nothing when the filter drops it.
/// The QC observations follow the trimmed path: the input read and its tags
/// are observed for every record, the tags against the surviving window or
/// against none when the read is dropped. Returns the output and whether a
/// known per-base tag is malformed.
fn process_raw_full_window(
    record: bam::Record,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<(Option<BamOutputRecord>, bool)> {
    let seq_len = record.sequence().len();
    let qualities = record.quality_scores();
    let qual = qualities.as_ref();
    if qual.len() != seq_len {
        let name = record
            .name()
            .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
            .unwrap_or_else(|| "<unnamed>".to_string());
        anyhow::bail!(
            "read {name}: BAM record SEQ length {} != QUAL length {} \
             (records without full per-base quality are not supported)",
            seq_len,
            qual.len()
        );
    }

    let (mod_block, malformed_perbase) = raw_full_window_metadata(&record)?;
    if mod_block == ModBlock::Malformed {
        counters.malformed_mod_reads.fetch_add(1, Ordering::Relaxed);
    }

    let mut tag_source: Option<RecordBuf> = None;

    let dropped = if seq_len == 0 {
        counters
            .reads_trimmed_to_nothing
            .fetch_add(1, Ordering::Relaxed);
        true
    } else {
        let gc = (cfg.filter.min_gc.is_some() || cfg.filter.max_gc.is_some())
            .then(|| raw_gc_fraction(&record));
        match crate::filter::check_metrics(seq_len, qual, gc, &cfg.filter) {
            Some(reason) => {
                counters.record_segment_drop(reason);
                counters.reads_all_filtered.fetch_add(1, Ordering::Relaxed);
                true
            },
            None => false,
        }
    };
    if dropped {
        return Ok((None, malformed_perbase));
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
            let decoded = tag_source
                .take()
                .map_or_else(|| decode_raw_record(&record), Ok)?;
            let seq = decoded.sequence().as_ref();
            let window = Window {
                start: 0,
                end: seq_len,
                idx: 0,
                total: 1,
            };
            BamOutputRecord::Decoded(reconstruct_record_with_bases(
                &decoded,
                seq,
                qual,
                window,
                mod_block,
                cfg.update_moves,
            ))
        },
    };
    Ok((Some(output), malformed_perbase))
}

fn run_raw_bam_full_window_seq(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>>,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    let mut malformed = 0;
    for record in records {
        let record = record?;
        ensure_raw_unaligned(&record)?;
        let seq_len = record.sequence().len();
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .input_bases
            .fetch_add(seq_len as u64, Ordering::Relaxed);
        let (output, is_malformed) = process_raw_full_window(record, cfg, counters)?;
        malformed += u64::from(is_malformed);
        match output {
            Some(BamOutputRecord::Raw(record)) => sink.write_raw_record(header, &record)?,
            Some(BamOutputRecord::Decoded(record)) => sink.write_record(header, &record)?,
            None => {},
        }
    }
    Ok(counters.snapshot(malformed))
}

fn run_raw_bam_full_window_parallel(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    run_parallel(
        records,
        BAM_BATCH,
        |record: &bam::Record| record.sequence().len(),
        cfg,
        sink,
        |record, cfg| {
            ensure_raw_unaligned(&record)?;
            let (output, malformed_tags) = process_raw_full_window(record, cfg, counters)?;
            Ok(Rendered {
                items: output.into_iter().collect(),
                malformed_tags,
            })
        },
        |sink, output| match output {
            BamOutputRecord::Raw(record) => sink.write_raw_record(header, record),
            BamOutputRecord::Decoded(record) => sink.write_record(header, record),
        },
        counters,
    )
}

/// Runs the uBAM workflow on raw records from a production reader. Full-window
/// runs filter and write unchanged records without building an owned
/// `RecordBuf`; any configuration that can alter sequence or tags is routed to
/// `run_bam`.
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
        && cfg.adapters.is_none();
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
pub fn run_bam<R: InputRecord>(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<R>> + Send,
    sink: &mut crate::io::bam::BamSink,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    if cfg.threads <= 1 {
        return run_bam_seq(header, records, sink, cfg, counters);
    }
    run_bam_parallel(
        records,
        cfg,
        sink,
        // Render: per-record guards and trim, then `process_read_segments` filters
        // per segment and reconstructs survivors into `Vec<RecordBuf>`.
        |rec, cfg| {
            let (seq, qual, mod_block) = prepare_read(rec, counters)?;
            let _read =
                crate::workflow::read_span(rec.name().map(|n| n.as_ref()).unwrap_or(b"<unnamed>"));
            let _read = _read.enter();
            let produced = trim::apply(seq, qual, &cfg.trim, cfg.adapters.as_ref());
            let mut items = Vec::with_capacity(produced.len());
            let mut survivors: Vec<(usize, usize)> = Vec::new();
            process_read_segments(
                &produced,
                seq,
                qual,
                &cfg.filter,
                counters,
                |idx, total, s, e| {
                    survivors.push((s, e));
                    let window = Window {
                        start: s,
                        end: e,
                        idx,
                        total,
                    };
                    items.push(reconstruct_record_with_bases(
                        rec,
                        seq,
                        qual,
                        window,
                        mod_block,
                        cfg.update_moves,
                    ));
                    Ok(())
                },
            )?;
            count_undo_tags_dropped(counters, rec, seq.len(), &survivors);
            Ok(items)
        },
        // Write: encode and write on the writer thread (BGZF compression is
        // multithreaded).
        |sink, rec| sink.write_record(header, rec),
        counters,
    )
}

/// Assembles the TAB-prefixed aux-tag block for one window: carried non-mod
/// tags in source order with the `window_tag_updates` rewrites applied in
/// place and the removed ones skipped, per-base arrays sliced, then the rebuilt
/// MM/ML/MN block, then the added tags. Empty when nothing is carried (the
/// caller writes a plain FASTQ record). A `Malformed` block is omitted.
fn build_fastq_tags(
    src: &RecordBuf,
    seq: &[u8],
    window: Window,
    mod_block: ModBlock,
    sel: &FastqTags,
    platform: Platform,
) -> Vec<u8> {
    let Window { start, end, .. } = window;
    let mut tags = Vec::new();
    let orig_len = seq.len();
    let trimmed = start != 0 || end != orig_len;
    // BAM-to-FASTQ never rewrites the move table (a sliced one is impractical
    // in a FASTQ header, and signal-aware consumers read BAM), so a trim drops
    // the signal and poly-A tags.
    let mut updates =
        window_tag_updates(src, src.quality_scores().as_ref(), window, platform, false);
    for (tag, value) in src.data().iter() {
        let t = <[u8; 2]>::from(tag);
        if matches!(&t, b"MM" | b"ML" | b"MN") {
            continue; // handled by the rebuilt block below
        }
        let rewritten = updates
            .iter()
            .position(|(u, _)| *u == tag)
            .map(|i| updates.remove(i).1);
        if !sel.carries(&t) {
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
        tags.extend_from_slice(&format_aux_field(t, &value));
    }
    if sel.carries_mods()
        && matches!(mod_block, ModBlock::Consistent | ModBlock::MissingMn)
        && let Some((mm, ml)) = mod_tags(src)
    {
        let (mm, ml) = rebuild_mods(mm, ml, seq, start, end);
        tags.push(b'\t');
        tags.extend_from_slice(&format_mods_aux(&mm, ml.as_deref(), end - start));
    }
    for (tag, value) in updates {
        let t = <[u8; 2]>::from(tag);
        if let Some(v) = value
            && sel.carries(&t)
        {
            tags.push(b'\t');
            tags.extend_from_slice(&format_aux_field(t, &v));
        }
    }
    tags
}

/// Writes one surviving window of a decoded record as a FASTQ record: the
/// platform's segment name (`segment_name`), the selected aux tags, then the
/// sliced bases and qualities.
fn write_fastq_window<W: Write>(
    w: &mut W,
    rec: &RecordBuf,
    window: Window,
    mod_block: ModBlock,
    platform: Platform,
    sel: &FastqTags,
) -> io::Result<()> {
    let Window { start, end, .. } = window;
    let seq = rec.sequence().as_ref();
    let qual = rec.quality_scores().as_ref();
    let coords = query_span(rec).map(|(qs0, _)| window_coords(qs0, start, end));
    let name = rec.name().map(|n| n.to_vec()).unwrap_or_default();
    let name = segment_name(platform, &name, window, coords);
    let tags = build_fastq_tags(rec, seq, window, mod_block, sel, platform);
    write_named_record(w, &name, &seq[start..end], &qual[start..end], &tags)
}

/// Runs the single-threaded uBAM-to-FASTQ workflow: refuses aligned reads,
/// trims, filters each produced segment, then writes each surviving segment as
/// FASTQ with the selected aux tags in the header (MM/ML/MN reconstructed,
/// per-base arrays sliced, other tags copied). gzip compression, when
/// requested, is handled by the parallel `gzp` writer this drains into.
fn run_bam_to_fastq_seq<R, W>(
    records: impl Iterator<Item = anyhow::Result<R>>,
    writer: &mut W,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats>
where
    R: InputRecord,
    W: Write,
{
    let mut malformed_tag_reads = 0u64;
    for rec in records {
        let rec = rec?.decode()?;
        let (seq, qual, mod_block) = prepare_read(&rec, counters)?;
        counters.input_reads.fetch_add(1, Ordering::Relaxed);
        counters
            .input_bases
            .fetch_add(seq.len() as u64, Ordering::Relaxed);
        if has_malformed_perbase_tag(&rec, seq.len()) {
            malformed_tag_reads += 1;
        }
        let platform = platform(&rec);
        let _read =
            crate::workflow::read_span(rec.name().map(|n| n.as_ref()).unwrap_or(b"<unnamed>"));
        let _read = _read.enter();
        let produced = trim::apply(seq, qual, &cfg.trim, cfg.adapters.as_ref());
        let mut survivors: Vec<(usize, usize)> = Vec::new();
        process_read_segments(
            &produced,
            seq,
            qual,
            &cfg.filter,
            counters,
            |idx, total, s, e| {
                survivors.push((s, e));
                let window = Window {
                    start: s,
                    end: e,
                    idx,
                    total,
                };
                write_fastq_window(writer, &rec, window, mod_block, platform, &cfg.fastq_tags)?;
                Ok(())
            },
        )?;
        count_undo_tags_dropped(counters, &rec, seq.len(), &survivors);
    }
    Ok(counters.snapshot(malformed_tag_reads))
}

/// Runs the uBAM-to-FASTQ workflow: decodes, refuses aligned reads, trims,
/// filters, then writes each surviving segment as FASTQ with the selected aux
/// tags in the header (MM/ML/MN reconstructed, per-base arrays sliced, other
/// tags copied). Sequential for `cfg.threads <= 1`; otherwise renders on a
/// rayon pool and drains through `run_bam_parallel`'s bounded channel, in
/// input order under `cfg.ordered` and in completion order otherwise.
pub fn run_bam_to_fastq<R: InputRecord, W: Write + Send>(
    records: impl Iterator<Item = anyhow::Result<R>> + Send,
    writer: &mut W,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    if cfg.threads <= 1 {
        return run_bam_to_fastq_seq(records, writer, cfg, counters);
    }
    run_bam_parallel(
        records,
        cfg,
        writer,
        // Render: guards and trim, then `process_read_segments` filters per
        // segment into `Vec<Vec<u8>>` (rendered FASTQ segments, survivors only).
        |rec, cfg| {
            let (seq, qual, mod_block) = prepare_read(rec, counters)?;
            let platform = platform(rec);
            let _read =
                crate::workflow::read_span(rec.name().map(|n| n.as_ref()).unwrap_or(b"<unnamed>"));
            let _read = _read.enter();
            let produced = trim::apply(seq, qual, &cfg.trim, cfg.adapters.as_ref());
            let mut out = Vec::with_capacity(produced.len());
            let mut survivors: Vec<(usize, usize)> = Vec::new();
            process_read_segments(
                &produced,
                seq,
                qual,
                &cfg.filter,
                counters,
                |idx, total, s, e| {
                    survivors.push((s, e));
                    let window = Window {
                        start: s,
                        end: e,
                        idx,
                        total,
                    };
                    let mut buf = Vec::new();
                    write_fastq_window(
                        &mut buf,
                        rec,
                        window,
                        mod_block,
                        platform,
                        &cfg.fastq_tags,
                    )?;
                    out.push(buf);
                    Ok(())
                },
            )?;
            count_undo_tags_dropped(counters, rec, seq.len(), &survivors);
            Ok(out)
        },
        // Write: append the rendered bytes to the `FastqOut` writer.
        |w, buf| w.write_all(buf),
        counters,
    )
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::cigar::{Op, op::Kind};
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    use super::*;

    #[test]
    fn resolves_long_cigar_overflow_from_cg_tag() {
        let mut record = RecordBuf::default();
        *record.sequence_mut() = b"ACGT".to_vec().into();
        record
            .cigar_mut()
            .as_mut()
            .extend([Op::new(Kind::SoftClip, 4), Op::new(Kind::Skip, 4)]);
        record.data_mut().insert(
            Tag::CIGAR,
            Value::Array(Array::UInt32(vec![2 << 4, (1 << 4) | 1, (1 << 4) | 7])),
        );

        resolve_long_cigar(&mut record).unwrap();

        assert_eq!(
            record.cigar().as_ref(),
            [
                Op::new(Kind::Match, 2),
                Op::new(Kind::Insertion, 1),
                Op::new(Kind::SequenceMatch, 1),
            ]
        );
        assert!(record.data().get(&Tag::CIGAR).is_none());
    }

    fn ubam_with_mods(seq: &[u8], quals: Vec<u8>, mm: &[u8], ml: Vec<u8>) -> RecordBuf {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = seq.to_vec().into();
        *rec.quality_scores_mut() = quals.into();
        let data = rec.data_mut();
        data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.to_vec().into()));
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(ml)),
        );
        rec
    }

    #[test]
    fn slices_seq_qual_and_rebuilds_tags() {
        // Seq CCAC; `C+m` modified at C occurrences 0 and 2, positions 0 and 3;
        // ML [10, 20].
        let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0,1;", vec![10, 20]);
        // Window [2,4) keeps "AC"; the modified C at position 3 survives as
        // window occurrence 0.
        let out = reconstruct_record(&src, 2, 4, 1, 0, false);

        assert_eq!(out.sequence().as_ref(), b"AC");
        assert_eq!(out.quality_scores().as_ref(), &[32, 33]);

        let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => s.to_vec(),
            _ => panic!("No MM"),
        };
        assert_eq!(mm, b"C+m,0;");
        let ml = match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) => v.clone(),
            _ => panic!("No ML"),
        };
        assert_eq!(ml, vec![20]);
        // MN updated to the output length.
        let mn = match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(n)) => *n,
            _ => panic!("No MN"),
        };
        assert_eq!(mn, 2);
    }

    /// A BAM record with unequal SEQ and QUAL lengths returns an error.
    #[test]
    fn qual_seq_length_mismatch_errors_without_panicking() {
        use crate::config::IoConfig;
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = b"ACGT".to_vec().into();
        // `quality_scores` is left at its default (empty), so SEQ and QUAL
        // lengths differ.

        let header = sam::Header::default();
        let dir = tempfile::tempdir().unwrap();
        let mut sink =
            crate::io::bam::writer(Some(&dir.path().join("o.bam")), &header, 1, 6).unwrap();

        let cfg = Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            ordered: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
            progress: crate::config::ProgressMode::Auto,
            adapters_configured: None,
        };

        let result = run_bam(
            &header,
            [Ok(rec)].into_iter(),
            &mut sink,
            &cfg,
            &Arc::new(Counters::default()),
        );
        assert!(
            result.is_err(),
            "SEQ/QUAL length mismatch must error, not panic"
        );
    }

    /// A non-string MM value is outside the supported schema and remains untouched.
    #[test]
    fn reconstruct_record_leaves_non_string_mm_untouched() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        let data = src.data_mut();
        data.insert(Tag::BASE_MODIFICATIONS, Value::Int32(5)); // spec-invalid MM
        data.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(vec![1, 2, 3])),
        );
        data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));

        let out = reconstruct_record(&src, 1, 4, 1, 0, false);

        match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::Int32(5)) => {},
            other => panic!("MM must be left untouched for a non-string value, got {other:?}"),
        }
        match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) if v == &[1u8, 2, 3] => {},
            other => panic!("ML must be left untouched, got {other:?}"),
        }
        match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(4)) => {},
            other => panic!("MN must be left untouched, got {other:?}"),
        }
    }

    /// An empty assessment group remains present alongside populated groups.
    #[test]
    fn reconstruct_record_preserves_originally_empty_mm_group() {
        // Seq CACA: C at 0, 2; A at 1, 3. `A+a` modifies A occurrence 0
        // (position 1). `C+m` is empty.
        let src = ubam_with_mods(b"CACA", vec![30, 31, 32, 33], b"A+a,0;C+m;", vec![7]);
        // Missing MN requires reconstruction even for the complete sequence.
        let out = reconstruct_record(&src, 0, 4, 1, 0, false);
        let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => s.to_vec(),
            other => panic!("MM must survive, got {other:?}"),
        };
        assert_eq!(
            mm, b"A+a,0;C+m;",
            "The empty C+m assessment group must be preserved"
        );
    }

    /// A segment that loses every listed position keeps the group with no
    /// positions: `C+m;` still declares the segment's C's canonical, which the
    /// group's absence would turn into "no call".
    #[test]
    fn split_suffixes_name_and_keeps_empty_mod_group() {
        let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0;", vec![10]); // mod at position 0
        // Segment [2,4) has no surviving C modification.
        let out = reconstruct_record(&src, 2, 4, 2, 1, false);
        // `.as_ref()` is ambiguous on `&BStr` (it implements both `AsRef<[u8]>`
        // and `AsRef<BStr>`); the turbofish selects the byte view.
        assert_eq!(AsRef::<[u8]>::as_ref(out.name().unwrap()), b"r1_segment_2");
        match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => assert_eq!(s.to_vec(), b"C+m;"),
            other => panic!("Empty group must be kept, got {other:?}"),
        }
        match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) => assert!(v.is_empty(), "ML must be empty"),
            other => panic!("ML must be an empty B:C array, got {other:?}"),
        }
        assert_eq!(
            out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH),
            Some(&Value::Int32(2))
        );
    }

    /// The FASTQ header spells an empty group as `MM:Z:C+m;` with a zero-length
    /// `ML:B:C` array.
    #[test]
    fn bam2fq_keeps_empty_mod_group_after_a_crop() {
        let rec = ubam_with_mods(b"CCAC", vec![40; 4], b"C+m,0;", vec![10]);
        let cfg = cfg_bam2fq(None, 2, FastqTags::All);
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(
            s.starts_with("@r1\tMM:Z:C+m;\tML:B:C\tMN:i:2\n"),
            "Got: {s:?}"
        );
    }

    use crate::config::{FastqTags, IoConfig};
    use crate::filter::FilterConfig;
    use crate::qual::QualMode;
    use crate::trim::{QualityOp, TrimPlan};

    fn cfg_bam2fq(quality: Option<QualityOp>, head: usize, tags: FastqTags) -> Config {
        Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head,
                tail: 0,
                quality,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 1,
            fastq_tags: tags,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            ordered: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
            progress: crate::config::ProgressMode::Auto,
            adapters_configured: None,
        }
    }

    /// `CCACCCAC` has C at 0, 1, 3, 4, 5, 7; `C+m,0,1,0` marks occurrences 0, 2,
    /// 3 (positions 0, 3, 4) with ML [10, 20, 30]. A head crop of 2 (window
    /// [2,8)) keeps positions 3 and 4, renumbered to `C+m,0,0;` with ML [20, 30]
    /// and MN 6.
    fn read2_with_mods_and_rg() -> RecordBuf {
        let mut rec = ubam_with_mods(b"CCACCCAC", vec![35; 8], b"C+m,0,1,0;", vec![10, 20, 30]);
        rec.data_mut()
            .insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
        rec.data_mut()
            .insert(Tag::READ_GROUP, Value::String(b"grp1".as_slice().into()));
        rec
    }

    #[test]
    fn bam2fq_all_carries_rg_and_reconstructed_mods() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::All);
        let mut out = Vec::new();
        let stats = run_bam_to_fastq(
            [Ok(read2_with_mods_and_rg())].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
        let s = String::from_utf8(out).unwrap();
        // The header carries RG unchanged and the reconstructed mod block; the
        // sequence is head-cropped by 2.
        assert!(
            s.starts_with("@r1\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6\n"),
            "Got: {s:?}"
        );
        assert!(s.contains("\nACCCAC\n+\n"), "Cropped sequence wrong: {s:?}");
    }

    #[test]
    fn bam2fq_only_mm_ml_drops_rg() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::parse("MM,ML").unwrap());
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(read2_with_mods_and_rg())].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(!s.contains("RG:Z"), "RG must be dropped: {s:?}");
        assert!(
            s.contains("MM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"),
            "Mods missing: {s:?}"
        );
    }

    #[test]
    fn bam2fq_none_is_plain_fastq() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::None);
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(read2_with_mods_and_rg())].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(out, b"@r1\nACCCAC\n+\nDDDDDD\n"); // 35+33 = 'D'
    }

    /// A split at the low-quality base gives each segment its own reconstructed
    /// mods.
    #[test]
    fn bam2fq_split_suffixes_and_segments_mods() {
        let cfg = cfg_bam2fq(
            Some(QualityOp::Split {
                cutoff: 20,
                window: 1,
            }),
            0,
            FastqTags::All,
        );
        // Seq CCAC, `C+m` at occurrences 0 and 2 (positions 0 and 3); quality
        // good, good, bad, good, so the split is [0,2), [3,4).
        let mut rec = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0,1;", vec![100, 200]);
        rec.data_mut()
            .insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));
        let mut out = Vec::new();
        let stats = run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(stats.output_reads, 2);
        let s = String::from_utf8(out).unwrap();
        // Segment 1 = [0,2) "CC" keeps the position-0 mod; segment 2 = [3,4) "C"
        // keeps the position-3 mod.
        assert!(
            s.contains("@r1_segment_1\tMM:Z:C+m,0;\tML:B:C,100\tMN:i:2"),
            "Segment 1: {s:?}"
        );
        assert!(
            s.contains("@r1_segment_2\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:1"),
            "Segment 2: {s:?}"
        );
    }

    #[test]
    fn bam2fq_no_mods_read_is_plain() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"plain".into());
        *rec.sequence_mut() = b"ACGT".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 4].into();
        let cfg = cfg_bam2fq(None, 0, FastqTags::All);
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(out, b"@plain\nACGT\n+\nIIII\n");
    }

    /// MM without optional ML remains MM-only after reconstruction.
    #[test]
    fn reconstruct_record_mm_without_ml_stays_mm_only() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"CCAC".to_vec().into(); // C at 0,1,3
        *src.quality_scores_mut() = vec![40; 4].into();
        src.data_mut().insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(b"C+m,0,1;".to_vec().into()),
        );
        // No ML and no MN.

        let out = reconstruct_record(&src, 0, 4, 1, 0, false);

        // MM is retained: both modified Cs are in the window, so `C+m,0,1;`.
        let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => s.to_vec(),
            other => panic!("Expected MM retained, got {other:?}"),
        };
        assert_eq!(mm, b"C+m,0,1;");
        // ML must be absent, never an empty array.
        assert!(
            out.data()
                .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
                .is_none(),
            "MM-only source must not gain an ML tag"
        );
        // MN set to the window length.
        match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(4)) => {},
            other => panic!("Expected MN=4, got {other:?}"),
        }
    }

    /// Companion for the BAM-to-FASTQ path: an MM-only source must emit a FASTQ
    /// header with `MM` + `MN` but no `ML:B:C` field.
    #[test]
    fn bam2fq_mm_without_ml_omits_ml_field() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = b"CCAC".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 4].into();
        rec.data_mut().insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(b"C+m,0,1;".to_vec().into()),
        );

        let cfg = cfg_bam2fq(None, 0, FastqTags::All);
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(
            s.contains("MM:Z:C+m,0,1;\tMN:i:4"),
            "Expected MM+MN, got: {s:?}"
        );
        assert!(
            !s.contains("ML:B"),
            "MM-only record must not emit an ML field: {s:?}"
        );
    }

    /// A record whose modification block cannot be placed on its sequence: the
    /// variants `malformed_mod_blocks` feeds through both output paths.
    fn malformed_mod_record(variant: &str) -> RecordBuf {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = b"CCCA".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 4].into();
        let d = rec.data_mut();
        let (mm, ml, mn) = match variant {
            // MM declares 3 positions, ML has 1 byte.
            "ml_short" => (&b"C+m,0,0,0;"[..], Value::Array(Array::UInt8(vec![5])), 4),
            // MM stops parsing at the `x`.
            "mm_garbled" => (
                &b"C+m,0,0x,0;"[..],
                Value::Array(Array::UInt8(vec![5, 6, 7])),
                4,
            ),
            // ML at a subtype other than B:C.
            "ml_subtype" => (
                &b"C+m,0,0,0;"[..],
                Value::Array(Array::Int8(vec![5, 6, 7])),
                4,
            ),
            // MN disagrees with the sequence length.
            "mn_mismatch" => (
                &b"C+m,0,0,0;"[..],
                Value::Array(Array::UInt8(vec![5, 6, 7])),
                40,
            ),
            other => panic!("Unknown variant {other}"),
        };
        d.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.to_vec().into()));
        d.insert(Tag::BASE_MODIFICATION_PROBABILITIES, ml);
        d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(mn));
        d.insert(Tag::READ_GROUP, Value::String(b"grp".as_slice().into()));
        rec
    }

    const MALFORMED_MOD_VARIANTS: [&str; 4] =
        ["ml_short", "mm_garbled", "ml_subtype", "mn_mismatch"];

    /// Each defect classifies as `Malformed`, on the full window and on a crop.
    #[test]
    fn malformed_mod_blocks_are_classified() {
        for variant in MALFORMED_MOD_VARIANTS {
            let rec = malformed_mod_record(variant);
            assert_eq!(inspect_mod_block(&rec, 4), ModBlock::Malformed, "{variant}");
        }
        let ok = ubam_with_mods(b"CCCA", vec![40; 4], b"C+m,0,0,0;", vec![5, 6, 7]);
        assert_eq!(inspect_mod_block(&ok, 4), ModBlock::MissingMn);
    }

    /// The BAM output drops the whole block, keeps every other tag, and the run
    /// counts the read once, whether or not the read is trimmed or split.
    #[test]
    fn malformed_mod_block_is_removed_and_counted_on_bam_output() {
        for variant in MALFORMED_MOD_VARIANTS {
            for head in [0, 1] {
                let mut cfg = cfg_bam2fq(None, head, FastqTags::All);
                cfg.trim.quality = Some(QualityOp::Split {
                    cutoff: 20,
                    window: 1,
                });
                let mut rec = malformed_mod_record(variant);
                // A low-quality base in the middle splits the read in two.
                *rec.quality_scores_mut() = vec![40, 40, 1, 40].into();
                let header = sam::Header::default();
                let dir = tempfile::tempdir().unwrap();
                let path = dir.path().join("o.bam");
                let mut sink = crate::io::bam::writer(Some(&path), &header, 1, 6).unwrap();
                let counters = Arc::new(Counters::default());
                let stats =
                    run_bam(&header, [Ok(rec)].into_iter(), &mut sink, &cfg, &counters).unwrap();
                sink.finish().unwrap();
                assert_eq!(
                    stats.malformed_mod_reads, 1,
                    "{variant} head={head}: counted once per read"
                );

                let bytes = std::fs::read(&path).unwrap();
                let mut reader = noodles_bam::io::Reader::new(bytes.as_slice());
                let h = reader.read_header().unwrap();
                let mut buf = RecordBuf::default();
                let mut n = 0;
                while reader.read_record_buf(&h, &mut buf).unwrap() != 0 {
                    n += 1;
                    for t in MOD_TAGS {
                        assert!(
                            buf.data().get(&t).is_none(),
                            "{variant} head={head}: {t:?} must be removed"
                        );
                    }
                    assert!(
                        buf.data().get(&Tag::READ_GROUP).is_some(),
                        "{variant}: other tags are copied"
                    );
                }
                assert_eq!(n, stats.output_reads as usize);
                assert!(n >= 1, "{variant} head={head}: segments written");
            }
        }
    }

    /// The FASTQ header carries no MM/ML/MN for a malformed block, and the run
    /// counts the read.
    #[test]
    fn malformed_mod_block_is_omitted_and_counted_on_fastq_output() {
        for variant in MALFORMED_MOD_VARIANTS {
            for head in [0, 2] {
                let cfg = cfg_bam2fq(None, head, FastqTags::All);
                let mut out = Vec::new();
                let stats = run_bam_to_fastq(
                    [Ok(malformed_mod_record(variant))].into_iter(),
                    &mut out,
                    &cfg,
                    &Arc::new(Counters::default()),
                )
                .unwrap();
                assert_eq!(stats.malformed_mod_reads, 1, "{variant} head={head}");
                let s = String::from_utf8(out).unwrap();
                assert!(
                    s.starts_with("@r1\tRG:Z:grp\n"),
                    "{variant} head={head}: mod block must be omitted, got {s:?}"
                );
            }
        }
    }

    /// A well-formed block is not counted, and an absent `MN` is added rather
    /// than treated as a defect.
    #[test]
    fn well_formed_mod_block_is_not_counted_and_gains_mn() {
        let rec = ubam_with_mods(b"CCCA", vec![40; 4], b"C+m,0,0,0;", vec![5, 6, 7]);
        let cfg = cfg_bam2fq(None, 0, FastqTags::All);
        let mut out = Vec::new();
        let stats = run_bam_to_fastq(
            [Ok(rec.clone())].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(stats.malformed_mod_reads, 0);
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("MM:Z:C+m,0,0,0;\tML:B:C,5,6,7\tMN:i:4"), "{s:?}");
        let out = reconstruct_record(&rec, 0, 4, 1, 0, false);
        assert_eq!(
            out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH),
            Some(&Value::Int32(4)),
            "MN is added when absent"
        );
    }

    /// PacBio per-base kinetics (`ip`/`pw`, length equal to the read length) are
    /// sliced with the sequence; the ONT `mv` (signal-space) tag is dropped on
    /// trim; the per-read RG is copied unchanged.
    #[test]
    fn reconstruct_record_slices_kinetics_and_drops_mv() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 6].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![10, 11, 12, 13, 14, 15])),
        );
        d.insert(
            Tag::new(b'p', b'w'),
            Value::Array(Array::UInt16(vec![20, 21, 22, 23, 24, 25])),
        );
        d.insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![5, 1, 0, 1, 0])),
        );
        d.insert(Tag::READ_GROUP, Value::String(b"grp".as_slice().into()));

        // Window [2,5) is "GTA" (head crop 2, tail crop 1).
        let out = reconstruct_record(&src, 2, 5, 1, 0, false);
        assert_eq!(out.sequence().as_ref(), b"GTA");
        match out.data().get(&Tag::new(b'i', b'p')) {
            Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[12, 13, 14]),
            other => panic!("Expected ip sliced to [2,5): {other:?}"),
        }
        match out.data().get(&Tag::new(b'p', b'w')) {
            Some(Value::Array(Array::UInt16(v))) => assert_eq!(v, &[22, 23, 24]),
            other => panic!("Expected pw sliced to [2,5): {other:?}"),
        }
        assert!(
            out.data().get(&Tag::new(b'm', b'v')).is_none(),
            "The mv tag must be dropped on trim"
        );
        assert!(
            matches!(out.data().get(&Tag::READ_GROUP), Some(Value::String(_))),
            "RG kept"
        );
    }

    #[test]
    fn reconstruct_record_slices_unknown_read_length_array_but_not_others() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        // An unknown B array whose length equals the read length is sliced
        // structurally.
        src.data_mut().insert(
            Tag::new(b'z', b'z'),
            Value::Array(Array::Int32(vec![1, 2, 3, 4])),
        );
        // A B array whose length differs from the read length is not per-base
        // and is left alone.
        src.data_mut()
            .insert(Tag::new(b'x', b'y'), Value::Array(Array::UInt8(vec![9, 9])));

        let out = reconstruct_record(&src, 1, 3, 1, 0, false); // window [1,3)
        match out.data().get(&Tag::new(b'z', b'z')) {
            Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 3]),
            other => panic!("Expected zz sliced: {other:?}"),
        }
        match out.data().get(&Tag::new(b'x', b'y')) {
            Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[9, 9]),
            other => panic!("Expected xy untouched: {other:?}"),
        }
    }

    #[test]
    fn reconstruct_record_untrimmed_keeps_kinetics_and_mv() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        src.data_mut().insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![1, 2, 3, 4])),
        );
        src.data_mut().insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![5, 1, 1])),
        );

        // Full window [0,4): nothing is trimmed, so everything is preserved.
        let out = reconstruct_record(&src, 0, 4, 1, 0, false);
        match out.data().get(&Tag::new(b'i', b'p')) {
            Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[1, 2, 3, 4]),
            other => panic!("Expected ip unchanged: {other:?}"),
        }
        assert!(
            out.data().get(&Tag::new(b'm', b'v')).is_some(),
            "The mv tag is kept when untrimmed"
        );
    }

    #[test]
    fn bam2fq_slices_kinetics_in_header() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = b"ACGTAC".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 6].into();
        rec.data_mut().insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![10, 11, 12, 13, 14, 15])),
        );

        let cfg = cfg_bam2fq(None, 2, FastqTags::All); // head crop 2, window [2,6)
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(
            s.contains("ip:B:C,12,13,14,15"),
            "Kinetics not sliced in the FASTQ header: {s:?}"
        );
    }

    #[test]
    fn malformed_perbase_tag_detected_and_left_untouched() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        // `ip` length 3 differs from read length 4: a malformed known per-base tag.
        src.data_mut().insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![1, 2, 3])),
        );

        assert!(has_malformed_perbase_tag(&src, 4));
        // It cannot be sliced, so it is left as is.
        let out = reconstruct_record(&src, 1, 3, 1, 0, false);
        match out.data().get(&Tag::new(b'i', b'p')) {
            Some(Value::Array(Array::UInt8(v))) => assert_eq!(v, &[1, 2, 3]),
            other => panic!("Expected the malformed ip left as is: {other:?}"),
        }
    }

    /// A synthetic move table: stride 2, 6 ones (one per base) at block indexes
    /// 0, 1, 3, 4, 6, 7, 8 blocks in total. Shared by the `--update-moves` tests.
    fn ubam_with_moves() -> RecordBuf {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 6].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![2, 1, 1, 0, 1, 1, 0, 1, 1])),
        );
        d.insert(Tag::new(b't', b's'), Value::Int32(10));
        // Consistent: ts0 + n_blocks*stride = 10 + 8*2 = 26.
        d.insert(Tag::new(b'n', b's'), Value::Int32(26));
        d.insert(
            Tag::new(b's', b't'),
            Value::String(b"2024-06-21T10:00:00Z".as_slice().into()),
        );
        d.insert(Tag::new(b'd', b'u'), Value::Float(5.0));
        src
    }

    #[test]
    fn update_moves_head_crop_slices_mv_bumps_ts_keeps_ns() {
        // Head crop 2 gives window [2,6): block_first = ones[2] = 3,
        // block_second = 8.
        let out = reconstruct_record(&ubam_with_moves(), 2, 6, 1, 0, true);
        assert_eq!(out.sequence().as_ref(), b"GTAC");
        assert_eq!(
            AsRef::<[u8]>::as_ref(out.name().unwrap()),
            b"r1",
            "A crop keeps the read name"
        );
        // `mv` = [stride] + moves[3..8] = [2] + [1,1,0,1,1].
        match out.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1, 1]),
            other => panic!("Unexpected mv: {other:?}"),
        }
        // `ts` becomes 10 + 3*2 = 16; `ns` = ts + span = 16 + (8-3)*2 = 26 (a
        // head-only crop leaves `ns` unchanged).
        match out.data().get(&Tag::new(b't', b's')) {
            Some(Value::Int32(16)) => {},
            other => panic!("Unexpected ts: {other:?}"),
        }
        match out.data().get(&Tag::new(b'n', b's')) {
            Some(Value::Int32(26)) => {},
            other => panic!("Unexpected ns: {other:?}"),
        }
        assert!(out.data().get(&Tag::new(b's', b'p')).is_none());
        assert!(out.data().get(&Tag::new(b'p', b'i')).is_none());
        // A crop keeps the read identity, so `st`/`du` stay.
        assert!(
            out.data().get(&Tag::new(b's', b't')).is_some(),
            "The st tag is kept on a crop"
        );
        assert!(
            out.data().get(&Tag::new(b'd', b'u')).is_some(),
            "The du tag is kept on a crop"
        );
    }

    #[test]
    fn update_moves_large_signal_offsets_use_uint32_not_wrapped_i32() {
        let mut src = ubam_with_moves();
        src.data_mut()
            .insert(Tag::new(b't', b's'), Value::UInt32(2_147_483_645));
        src.data_mut()
            .insert(Tag::new(b'n', b's'), Value::UInt32(2_147_483_661));

        // Head crop 2 gives block_first = 3 with stride 2, so `ts` becomes
        // 2_147_483_651 (above `i32::MAX`) and `ns` becomes 2_147_483_661.
        let out = reconstruct_record(&src, 2, 6, 1, 0, true);

        match out.data().get(&Tag::new(b't', b's')) {
            Some(Value::UInt32(2_147_483_651)) => {},
            other => panic!("Large ts must stay positive as UInt32, got {other:?}"),
        }
        match out.data().get(&Tag::new(b'n', b's')) {
            Some(Value::UInt32(2_147_483_661)) => {},
            other => panic!("Large ns must stay positive as UInt32, got {other:?}"),
        }
    }

    #[test]
    fn update_moves_tail_crop_shrinks_ns() {
        // Tail crop 2 gives window [0,4): block_first = ones[0] = 0,
        // block_second = ones[4] = 6.
        let out = reconstruct_record(&ubam_with_moves(), 0, 4, 1, 0, true);
        // `mv` = [stride] + moves[0..6] = [2] + [1,1,0,1,1,0].
        match out.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1, 1, 0]),
            other => panic!("Unexpected mv: {other:?}"),
        }
        // `ts` is unchanged (no head trim): 10; `ns` = 10 + (6-0)*2 = 22, below 26.
        match out.data().get(&Tag::new(b't', b's')) {
            Some(Value::Int32(10)) => {},
            other => panic!("Unexpected ts: {other:?}"),
        }
        match out.data().get(&Tag::new(b'n', b's')) {
            Some(Value::Int32(22)) => {},
            other => {
                panic!("The ns tag must shrink on a tail crop (dorado ns = trim + span): {other:?}")
            },
        }
    }

    #[test]
    fn update_moves_split_emits_subread_tags() {
        // Split into [0,3) and [3,6): each is a dorado-style subread.
        let s1 = reconstruct_record(&ubam_with_moves(), 0, 3, 2, 0, true);
        assert_eq!(AsRef::<[u8]>::as_ref(s1.name().unwrap()), b"r1_segment_1");
        // `mv` = [2] + moves[ones[0]=0 .. ones[3]=4] = [2] + [1,1,0,1].
        match s1.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1]),
            other => panic!("Unexpected s1 mv: {other:?}"),
        }
        match s1.data().get(&Tag::new(b't', b's')) {
            Some(Value::Int32(0)) => {},
            o => panic!("Segment 1 ts should be 0: {o:?}"),
        }
        match s1.data().get(&Tag::new(b'n', b's')) {
            Some(Value::Int32(8)) => {}, // (block 4-0)*stride 2
            o => panic!("Unexpected s1 ns: {o:?}"),
        }
        // `sp` counts from the parent's POD5 signal start: ts0 10 + block 0 *
        // stride 2 (dorado adds `num_trimmed_samples` to the split point).
        match s1.data().get(&Tag::new(b's', b'p')) {
            Some(Value::Int32(10)) => {},
            o => panic!("Unexpected s1 sp: {o:?}"),
        }
        match s1.data().get(&Tag::new(b'p', b'i')) {
            Some(Value::String(s)) => assert_eq!(s.to_vec(), b"r1"),
            o => panic!("Unexpected s1 pi: {o:?}"),
        }
        // Dorado marks split products with read number -1.
        match s1.data().get(&Tag::new(b'r', b'n')) {
            Some(Value::Int32(-1)) => {},
            o => panic!("Segment 1 rn should be -1: {o:?}"),
        }
        // `st`/`du` describe the parent read; with the signal window known
        // they are recomputed for the subread (`update_moves_split_recomputes_st_and_du`).
        assert!(
            matches!(s1.data().get(&Tag::new(b's', b't')), Some(Value::String(s)) if s.to_vec() != b"2024-06-21T10:00:00Z"),
            "The st tag is recomputed on a split"
        );
        assert!(
            matches!(s1.data().get(&Tag::new(b'd', b'u')), Some(Value::Float(d)) if *d < 5.0),
            "The du tag is recomputed on a split"
        );

        let s2 = reconstruct_record(&ubam_with_moves(), 3, 6, 2, 1, true);
        assert_eq!(AsRef::<[u8]>::as_ref(s2.name().unwrap()), b"r1_segment_2");
        // `mv` = [2] + moves[ones[3]=4 .. 8] = [2] + [1,0,1,1].
        match s2.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 0, 1, 1]),
            other => panic!("Unexpected s2 mv: {other:?}"),
        }
        match s2.data().get(&Tag::new(b'n', b's')) {
            Some(Value::Int32(8)) => {}, // (8-4)*2
            o => panic!("Unexpected s2 ns: {o:?}"),
        }
        match s2.data().get(&Tag::new(b's', b'p')) {
            Some(Value::Int32(18)) => {}, // ts0 10 + block_first 4 * stride 2
            o => panic!("Unexpected s2 sp: {o:?}"),
        }
    }

    #[test]
    fn default_drops_all_signal_tags_on_trim() {
        let mut src = ubam_with_moves();
        src.data_mut().insert(Tag::new(b's', b'p'), Value::Int32(5));
        src.data_mut().insert(
            Tag::new(b'p', b'i'),
            Value::String(b"parent".as_slice().into()),
        );

        // `update_moves` off and trimmed: mv/ts/ns/sp/pi are all removed.
        let out = reconstruct_record(&src, 2, 6, 1, 0, false);
        for t in [b"mv", b"ts", b"ns", b"sp", b"pi"] {
            assert!(
                out.data().get(&Tag::new(t[0], t[1])).is_none(),
                "{} must be dropped by default on trim",
                std::str::from_utf8(t).unwrap()
            );
        }
    }

    #[test]
    fn trim_drops_polya_barcode_tags_and_refreshes_qs() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        // First two bases low quality (phred 2), the rest Q40.
        *src.quality_scores_mut() = vec![2, 2, 40, 40, 40, 40].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![100, 200, 300, 400, 500])),
        );
        d.insert(Tag::new(b'p', b't'), Value::Int32(50));
        d.insert(
            Tag::new(b'b', b'i'),
            Value::Array(Array::Float(vec![0.9, 5.0, 20.0])),
        );
        d.insert(Tag::new(b'q', b's'), Value::Float(20.0)); // whole-read qs (stale after crop)
        d.insert(
            Tag::new(b'R', b'G'),
            Value::String(b"grp".as_slice().into()),
        );

        // Head crop 2 gives window [2,6), which keeps only the Q40 bases.
        let out = reconstruct_record(&src, 2, 6, 1, 0, false);

        // The poly-A and barcode coordinate tags cannot be reconstructed and are
        // dropped.
        for t in [b"pa", b"pt", b"bi"] {
            assert!(
                out.data().get(&Tag::new(t[0], t[1])).is_none(),
                "{} must be dropped on trim",
                std::str::from_utf8(t).unwrap()
            );
        }
        // `qs` is recomputed from the trimmed (all-Q40) quality, not left at 20.
        match out.data().get(&Tag::new(b'q', b's')) {
            Some(Value::Float(q)) => {
                let expected = crate::qual::mean_prob_q(&[40, 40, 40, 40]) as f32;
                assert!(
                    (q - expected).abs() < 1e-4,
                    "Recomputed qs: got {q}, want {expected}"
                );
            },
            other => panic!("Unexpected qs: {other:?}"),
        }
        // Per-read metadata (RG) is untouched.
        assert!(out.data().get(&Tag::new(b'R', b'G')).is_some());
    }

    /// `ubam_with_moves` with head crop 2 spans the original-signal window
    /// [ts0+3*2, ts0+8*2] = [16, 26]; a split segment [3,6) spans
    /// [ts0+4*2, 26] = [18, 26]. The poly-A tags survive a crop that keeps the
    /// tail.
    #[test]
    fn update_moves_crop_keeps_polya_when_tail_survives() {
        let mut src = ubam_with_moves();
        // Anchor and boundaries all inside [16,26]; the split range is a sentinel.
        src.data_mut().insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![20, 18, 24, -1, -1])),
        );
        src.data_mut()
            .insert(Tag::new(b'p', b't'), Value::Int32(30));

        let out = reconstruct_record(&src, 2, 6, 1, 0, true); // head crop 2
        // A crop keeps the read identity and POD5 signal, so the absolute `pa`
        // stays valid.
        match out.data().get(&Tag::new(b'p', b'a')) {
            Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[20, 18, 24, -1, -1]),
            other => panic!("Expected pa kept as is on a crop: {other:?}"),
        }
        match out.data().get(&Tag::new(b'p', b't')) {
            Some(Value::Int32(30)) => {},
            other => panic!("Unexpected pt: {other:?}"),
        }
    }

    #[test]
    fn update_moves_split_shifts_polya_into_subread_frame() {
        let mut src = ubam_with_moves();
        src.data_mut().insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![20, 18, 24, -1, -1])),
        );
        src.data_mut()
            .insert(Tag::new(b'p', b't'), Value::Int32(30));

        // Split segment [3,6): kept signal window [18,26], so real positions
        // shift by -18.
        let out = reconstruct_record(&src, 3, 6, 2, 1, true);
        match out.data().get(&Tag::new(b'p', b'a')) {
            Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 0, 6, -1, -1]),
            other => panic!("Expected pa shifted into the subread frame: {other:?}"),
        }
        match out.data().get(&Tag::new(b'p', b't')) {
            Some(Value::Int32(30)) => {}, // base count unchanged
            other => panic!("Unexpected pt: {other:?}"),
        }
    }

    #[test]
    fn update_moves_drops_polya_when_tail_trimmed() {
        let mut src = ubam_with_moves();
        // Anchor at 12 sits in the trimmed-off front signal (kept window is [16,26]).
        src.data_mut().insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![12, 10, 14, -1, -1])),
        );
        src.data_mut()
            .insert(Tag::new(b'p', b't'), Value::Int32(30));

        let out = reconstruct_record(&src, 2, 6, 1, 0, true); // head crop 2
        assert!(
            out.data().get(&Tag::new(b'p', b'a')).is_none(),
            "The pa tag is dropped when the tail is trimmed"
        );
        assert!(
            out.data().get(&Tag::new(b'p', b't')).is_none(),
            "The pt tag is dropped when the tail is trimmed"
        );
    }

    /// `pa` uses signal coordinates and is not a per-base array.
    #[test]
    fn update_moves_does_not_reslice_read_length_pa() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTA".to_vec().into(); // 5 bases
        *src.quality_scores_mut() = vec![40; 5].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![2, 1, 1, 1, 1, 1])),
        ); // stride 2, 5 ones
        d.insert(Tag::new(b't', b's'), Value::Int32(0));
        d.insert(Tag::new(b'n', b's'), Value::Int32(10));
        // A 5-element `pa` (equal to the read length) with all real positions
        // inside the kept window.
        d.insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![4, 2, 6, -1, -1])),
        );

        // Head crop 1 gives window [1,5): kept signal window [2,10]; `pa` survives.
        let out = reconstruct_record(&src, 1, 5, 1, 0, true);
        match out.data().get(&Tag::new(b'p', b'a')) {
            Some(Value::Array(Array::Int32(v))) => {
                assert_eq!(v, &[4, 2, 6, -1, -1], "The pa tag must not be re-sliced")
            },
            other => panic!("Unexpected pa: {other:?}"),
        }
    }

    /// `--update-moves` without a move table cannot relate signal to sequence,
    /// so the signal and poly-A tags are dropped (`parse_move_table` returns
    /// `None`, which selects `drop_all`).
    #[test]
    fn update_moves_without_move_table_drops_signal_and_polya() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 6].into();
        let d = src.data_mut();
        d.insert(Tag::new(b't', b's'), Value::Int32(10));
        d.insert(Tag::new(b'n', b's'), Value::Int32(100));
        d.insert(
            Tag::new(b'p', b'a'),
            Value::Array(Array::Int32(vec![20, 18, 24, -1, -1])),
        );
        d.insert(Tag::new(b'p', b't'), Value::Int32(30));

        let out = reconstruct_record(&src, 2, 6, 1, 0, true);
        for t in [b"ts", b"ns", b"pa", b"pt"] {
            assert!(
                out.data().get(&Tag::new(t[0], t[1])).is_none(),
                "{} dropped when the move table is absent",
                std::str::from_utf8(t).unwrap()
            );
        }
    }

    #[test]
    fn update_moves_polya_boundary_end_inclusive_anchor_exclusive() {
        // Split [3,6): kept window [18,26). A range end exactly at `kept_end`
        // (exclusive) survives; an anchor at `kept_end` is outside the window and
        // drops the tags.
        let mk = |pa: Vec<i32>| {
            let mut src = ubam_with_moves();
            src.data_mut()
                .insert(Tag::new(b'p', b'a'), Value::Array(Array::Int32(pa)));
            src.data_mut()
                .insert(Tag::new(b'p', b't'), Value::Int32(30));
            src
        };
        // Range end equal to `kept_end` (26) survives, shifted by -18.
        let kept = reconstruct_record(&mk(vec![20, 18, 26, -1, -1]), 3, 6, 2, 1, true);
        match kept.data().get(&Tag::new(b'p', b'a')) {
            Some(Value::Array(Array::Int32(v))) => assert_eq!(v, &[2, 0, 8, -1, -1]),
            other => panic!("Range end at kept_end should survive: {other:?}"),
        }
        // Anchor equal to `kept_end` (26) is outside the window, so the tags are
        // dropped.
        let dropped = reconstruct_record(&mk(vec![26, 18, 24, -1, -1]), 3, 6, 2, 1, true);
        assert!(
            dropped.data().get(&Tag::new(b'p', b'a')).is_none(),
            "Anchor at the exclusive boundary drops the tags"
        );
    }

    #[test]
    fn run_bam_parallel_matches_sequential_as_multiset() {
        use crate::config::{FastqTags, IoConfig};
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        let mk = |threads| Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 2,
                tail: 2,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads,
            fastq_tags: FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            ordered: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
            progress: crate::config::ProgressMode::Auto,
            adapters_configured: None,
        };
        // 300 reads with mods so reconstruction runs on every one.
        let recs: Vec<RecordBuf> = (0..300)
            .map(|_| ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]))
            .collect();

        let header = sam::Header::default();
        let decode = |bytes: &[u8]| -> Vec<(Vec<u8>, Vec<u8>)> {
            // (seq, MM-bytes) pairs, sorted, as an order-independent fingerprint.
            let mut r = noodles_bam::io::Reader::new(bytes);
            let h = r.read_header().unwrap();
            let mut out = Vec::new();
            let mut buf = RecordBuf::default();
            while r.read_record_buf(&h, &mut buf).unwrap() != 0 {
                let seq = buf.sequence().as_ref().to_vec();
                let mm = match buf.data().get(&Tag::BASE_MODIFICATIONS) {
                    Some(Value::String(s)) => s.to_vec(),
                    _ => Vec::new(),
                };
                out.push((seq, mm));
            }
            out.sort();
            out
        };

        // t1: single-threaded BGZF sink, written to a temporary file.
        let dir = tempfile::tempdir().unwrap();
        let p1 = dir.path().join("t1.bam");
        let mut sink1 = crate::io::bam::writer(Some(&p1), &header, 1, 6).unwrap();
        run_bam(
            &header,
            recs.clone().into_iter().map(anyhow::Ok),
            &mut sink1,
            &mk(1),
            &Arc::new(Counters::default()),
        )
        .unwrap();
        sink1.finish().unwrap();
        let b1 = std::fs::read(&p1).unwrap();

        // t8: multithreaded sink to a temporary file (the multithreaded writer
        // needs an owned `Write + Send`).
        let p8 = dir.path().join("t8.bam");
        let mut sink8 = crate::io::bam::writer(Some(&p8), &header, 8, 6).unwrap();
        run_bam(
            &header,
            recs.into_iter().map(anyhow::Ok),
            &mut sink8,
            &mk(8),
            &Arc::new(Counters::default()),
        )
        .unwrap();
        sink8.finish().unwrap();
        let b8 = std::fs::read(&p8).unwrap();

        assert_eq!(
            decode(&b1),
            decode(&b8),
            "The t1 and t8 runs must produce the same record set"
        );
    }

    /// Writer errors remain observable after the bounded channel reaches capacity.
    #[test]
    fn run_bam_parallel_surfaces_write_error_without_deadlock() {
        use std::io;

        use crate::config::IoConfig;
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        struct FailAfter {
            limit: usize,
            written: usize,
        }

        let cfg = Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 4,
            fastq_tags: crate::config::FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            ordered: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
            progress: crate::config::ProgressMode::Auto,
            adapters_configured: None,
        };
        let recs: Vec<anyhow::Result<RecordBuf>> = (0..3000)
            .map(|_| anyhow::Ok(RecordBuf::default()))
            .collect();

        let mut sink = FailAfter {
            limit: 100,
            written: 0,
        };
        let res = run_bam_parallel(
            recs.into_iter(),
            &cfg,
            &mut sink,
            |_rec, _cfg| anyhow::Ok(vec![()]),
            |sink, _item: &()| -> io::Result<()> {
                if sink.written >= sink.limit {
                    return Err(io::Error::new(io::ErrorKind::BrokenPipe, "boom"));
                }
                sink.written += 1;
                Ok(())
            },
            &Arc::new(Counters::default()),
        );
        assert!(
            res.is_err(),
            "Write error must surface as Err and must not hang"
        );
    }

    /// Mirrors `workflow::fastq`'s `parallel_surfaces_parse_error_instead_of_dropping_it`,
    /// driving `run_bam_parallel` directly so a malformed upstream record (an
    /// `Err` item from the input iterator) is not silently swallowed.
    #[test]
    fn run_bam_parallel_surfaces_parse_error_instead_of_dropping_it() {
        use std::io;

        use crate::config::IoConfig;
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        struct NullSink;

        let cfg = Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 0,
                tail: 0,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads: 4,
            fastq_tags: crate::config::FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            ordered: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
            progress: crate::config::ProgressMode::Auto,
            adapters_configured: None,
        };
        let good: Vec<anyhow::Result<RecordBuf>> =
            (0..5).map(|_| anyhow::Ok(RecordBuf::default())).collect();
        let recs = good
            .into_iter()
            .chain(std::iter::once(Err(anyhow::anyhow!("bad record"))));

        let mut sink = NullSink;
        let res = run_bam_parallel(
            recs,
            &cfg,
            &mut sink,
            |_rec, _cfg| anyhow::Ok(vec![()]),
            |_sink: &mut NullSink, _item: &()| -> io::Result<()> { Ok(()) },
            &Arc::new(Counters::default()),
        );
        assert!(
            res.is_err(),
            "A malformed record must not be dropped on the parallel path"
        );
    }

    #[test]
    fn run_bam_to_fastq_parallel_matches_sequential_as_multiset() {
        use crate::config::{FastqTags, IoConfig};
        use crate::filter::FilterConfig;
        use crate::qual::QualMode;
        use crate::trim::TrimPlan;

        let mk = |threads| Config {
            io: IoConfig {
                input: None,
                output: None,
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 1,
                max_length: usize::MAX,
                min_qual: 0.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 2,
                tail: 2,
                quality: None,
            },
            adapters: None,
            adapter_infer: crate::config::AdapterInfer::Off,
            threads,
            fastq_tags: FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            ordered: false,
            verbosity: 0,
            quiet: true,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
            progress: crate::config::ProgressMode::Auto,
            adapters_configured: None,
        };
        let recs: Vec<RecordBuf> = (0..300)
            .map(|_| ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]))
            .collect();

        let sorted_records = |bytes: &[u8]| {
            let s = String::from_utf8(bytes.to_vec()).unwrap();
            // Records are grouped as 4 consecutive lines rather than split on
            // `@`: a QUAL byte of Phred 31 (ASCII `@`) would corrupt an `@`
            // split. FASTQ records are exactly 4 lines each here, so the
            // re-chunking is lossless.
            let lines: Vec<&str> = s.lines().collect();
            assert_eq!(
                lines.len() % 4,
                0,
                "Expected whole 4-line FASTQ records, got {} lines",
                lines.len()
            );
            let mut v: Vec<String> = lines.chunks(4).map(|c| c.join("\n")).collect();
            v.sort();
            v
        };

        let mut a = Vec::new();
        run_bam_to_fastq(
            recs.clone().into_iter().map(anyhow::Ok),
            &mut a,
            &mk(1),
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let mut b = Vec::new();
        run_bam_to_fastq(
            recs.into_iter().map(anyhow::Ok),
            &mut b,
            &mk(8),
            &Arc::new(Counters::default()),
        )
        .unwrap();

        assert_eq!(
            sorted_records(&a),
            sorted_records(&b),
            "The t1 and t8 FASTQ outputs must match as a multiset"
        );
    }

    #[test]
    fn interior_adapter_split_reconstructs_mods_per_segment() {
        use crate::adapter::{Adapter, AdapterConfig, End};
        // Seq (64 bp): [flank1: C + 23 A][adapter GGGGTTTTGGGGTTTT (no C/A)][flank2: C + 23 A].
        // Only two Cs, at positions 0 and 40. `C+m,0,0;` marks both, with ML
        // [100, 200].
        let mut seq = b"CAAAAAAAAAAAAAAAAAAAAAAA".to_vec(); // C at 0
        seq.extend_from_slice(b"GGGGTTTTGGGGTTTT"); // interior adapter, 16 bp
        seq.extend_from_slice(b"CAAAAAAAAAAAAAAAAAAAAAAA"); // C at 40
        let quals = vec![40u8; seq.len()];
        let mut rec = ubam_with_mods(&seq, quals, b"C+m,0,0;", vec![100, 200]);
        rec.data_mut().insert(
            Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
            Value::Int32(seq.len() as i32),
        );

        // BAM-to-FASTQ path (renders MM/ML/MN as header text), adapters active,
        // split on.
        let mut cfg = cfg_bam2fq(None, 0, FastqTags::All);
        cfg.adapters = Some(AdapterConfig {
            adapters: vec![Adapter {
                name: "mid".into(),
                seq: b"GGGGTTTTGGGGTTTT".to_vec(),
                end: End::Both,
            }],
            error_rate: 0.2,
            end_size: 8, // adapter at [24,40) is interior, more than 8 from both ends of 64 bp
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        });

        let mut out = Vec::new();
        let stats = run_bam_to_fastq(
            [Ok(rec)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(
            stats.output_reads, 2,
            "Interior adapter splits into two subreads"
        );
        let s = String::from_utf8(out).unwrap();
        // Each segment keeps exactly its own C mod, renumbered to occurrence 0.
        assert!(
            s.contains("@r1_segment_1\tMM:Z:C+m,0;\tML:B:C,100\tMN:i:24"),
            "Segment 1 mods wrong: {s}"
        );
        assert!(
            s.contains("@r1_segment_2\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:24"),
            "Segment 2 mods wrong: {s}"
        );
    }
    /// `MN` is written at the smallest integer subtype that fits, so dorado emits
    /// `MN:S` for an ordinary-length read. Accepting only `Int32` would make
    /// every real mod-bearing record look inconsistent, forcing a full MM/ML
    /// rebuild on an untrimmed record and rewriting `MN` as the wider `i`.
    #[test]
    fn mn_is_recognized_at_every_integer_subtype() {
        use noodles_sam::alignment::record_buf::data::field::Value;

        for v in [
            Value::UInt8(12),
            Value::Int8(12),
            Value::UInt16(12),
            Value::Int16(12),
            Value::UInt32(12),
            Value::Int32(12),
        ] {
            let mut rec = RecordBuf::default();
            *rec.sequence_mut() = b"ACGTACGTACGT".to_vec().into();
            let d = rec.data_mut();
            d.insert(
                Tag::BASE_MODIFICATIONS,
                Value::String(b"C+m,0;".to_vec().into()),
            );
            d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, v.clone());
            assert_eq!(
                inspect_mod_block(&rec, 12),
                ModBlock::Consistent,
                "MN stored as {v:?} should be recognized"
            );
        }
    }

    /// An untrimmed record takes the pass-through path, so its tags come out
    /// exactly as they went in, including `MN`'s storage width.
    #[test]
    fn untrimmed_mod_record_passes_through_byte_for_byte() {
        use noodles_sam::alignment::record_buf::data::field::Value;

        let mut rec = RecordBuf::default();
        *rec.flags_mut() = noodles_sam::alignment::record::Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = b"ACGTACGTACGT".to_vec().into();
        *rec.quality_scores_mut() = vec![40u8; 12].into();
        let d = rec.data_mut();
        d.insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(b"C+m,0,1;".to_vec().into()),
        );
        d.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(vec![200, 201])),
        );
        d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::UInt16(12));

        let out = reconstruct_record(&rec, 0, 12, 1, 0, false);
        assert_eq!(
            out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH),
            Some(&Value::UInt16(12)),
            "MN must keep the subtype it arrived with"
        );
        assert_eq!(out, rec, "An untrimmed record is passed through unchanged");
    }

    /// Encodes a decoded record as the raw BAM record a production reader yields.
    fn raw_record(rec: &RecordBuf) -> bam::Record {
        use noodles_sam::alignment::io::Write as _;

        let header = sam::Header::default();
        let mut bytes = Vec::new();
        {
            let mut w = bam::io::Writer::new(&mut bytes);
            w.write_header(&header).unwrap();
            w.write_alignment_record(&header, rec).unwrap();
            w.try_finish().unwrap();
        }
        let mut r = bam::io::Reader::new(bytes.as_slice());
        r.read_header().unwrap();
        let mut raw = bam::Record::default();
        assert_ne!(r.read_record(&mut raw).unwrap(), 0);
        raw
    }

    /// A 10-base PacBio-style record: `fi[i]` belongs to base `i`, `ri[i]` to
    /// base `9 - i` (reverse kinetics are stored last base first).
    fn pacbio_kinetics_record() -> RecordBuf {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 10].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b'f', b'i'),
            Value::Array(Array::UInt8((0..10).collect())),
        );
        d.insert(
            Tag::new(b'r', b'i'),
            Value::Array(Array::UInt8((10..20).collect())),
        );
        src
    }

    fn u8_array(rec: &RecordBuf, tag: [u8; 2]) -> Vec<u8> {
        match rec.data().get(&Tag::new(tag[0], tag[1])) {
            Some(Value::Array(Array::UInt8(v))) => v.clone(),
            other => panic!("{}: {other:?}", std::str::from_utf8(&tag).unwrap()),
        }
    }

    /// Reverse-strand kinetics take array indexes `[len - end, len - start)`,
    /// so a head crop removes entries from the end of `ri` and a tail crop from
    /// its start, while `fi` is sliced in read order.
    #[test]
    fn reverse_strand_kinetics_are_sliced_from_the_other_end() {
        let src = pacbio_kinetics_record();
        // Head crop 2: bases 2..=9 keep fi[2..10] and ri[0..8] (bases 9..2).
        let head = reconstruct_record(&src, 2, 10, 1, 0, false);
        assert_eq!(u8_array(&head, *b"fi"), (2..10).collect::<Vec<u8>>());
        assert_eq!(u8_array(&head, *b"ri"), (10..18).collect::<Vec<u8>>());
        // Tail crop 3: bases 0..=6 keep fi[0..7] and ri[3..10] (bases 6..0).
        let tail = reconstruct_record(&src, 0, 7, 1, 0, false);
        assert_eq!(u8_array(&tail, *b"fi"), (0..7).collect::<Vec<u8>>());
        assert_eq!(u8_array(&tail, *b"ri"), (13..20).collect::<Vec<u8>>());
    }

    #[test]
    fn bam2fq_slices_reverse_strand_kinetics_from_the_other_end() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::All); // head crop 2
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(pacbio_kinetics_record())].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("\tfi:B:C,2,3,4,5,6,7,8,9\t"), "{s:?}");
        assert!(s.contains("\tri:B:C,10,11,12,13,14,15,16,17\n"), "{s:?}");
    }

    /// PacBio's `qs:i`/`qe:i` are query coordinates, not a quality, and follow
    /// the window in the original read's frame; only a float `qs` is a quality
    /// and is recomputed.
    #[test]
    fn integer_qs_and_qe_follow_the_window() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![2, 2, 40, 40, 40, 40].into();
        src.data_mut()
            .insert(Tag::new(b'q', b's'), Value::Int32(1200));
        src.data_mut()
            .insert(Tag::new(b'q', b'e'), Value::Int32(1206));

        let out = reconstruct_record(&src, 2, 6, 1, 0, false);
        assert_eq!(
            out.data().get(&Tag::new(b'q', b's')),
            Some(&Value::Int32(1202))
        );
        assert_eq!(
            out.data().get(&Tag::new(b'q', b'e')),
            Some(&Value::Int32(1206))
        );

        let cfg = cfg_bam2fq(None, 2, FastqTags::All);
        let mut fastq = Vec::new();
        run_bam_to_fastq(
            [Ok(src.clone())].into_iter(),
            &mut fastq,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(fastq).unwrap();
        assert!(s.starts_with("@r1\tqs:i:1202\tqe:i:1206\n"), "{s:?}");

        // Dorado's float qs follows the trimmed quality on both outputs.
        src.data_mut()
            .insert(Tag::new(b'q', b's'), Value::Float(20.0));
        let expected = crate::qual::mean_prob_q(&[40, 40, 40, 40]) as f32;
        match reconstruct_record(&src, 2, 6, 1, 0, false)
            .data()
            .get(&Tag::new(b'q', b's'))
        {
            Some(Value::Float(q)) => assert!((q - expected).abs() < 1e-4),
            other => panic!("Unexpected qs: {other:?}"),
        }
        let mut fastq = Vec::new();
        run_bam_to_fastq(
            [Ok(src)].into_iter(),
            &mut fastq,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(fastq).unwrap();
        assert!(s.starts_with("@r1\tqs:f:"), "{s:?}");
        assert!(!s.contains("qs:f:20"), "{s:?}");
    }

    /// The untrimmed fast path removes a block whose `MN` disagrees with the
    /// sequence instead of rewriting `MN`, adds a missing `MN`, and passes a
    /// consistent record through raw.
    #[test]
    fn raw_full_window_removes_malformed_block_and_adds_missing_mn() {
        let cfg = cfg_bam2fq(None, 0, FastqTags::All);
        let counters = Arc::new(Counters::default());

        let (out, _) = process_raw_full_window(
            raw_record(&malformed_mod_record("mn_mismatch")),
            &cfg,
            &counters,
        )
        .unwrap();
        let Some(BamOutputRecord::Decoded(rec)) = out else {
            panic!("A malformed block forces a rebuild");
        };
        for t in MOD_TAGS {
            assert!(rec.data().get(&t).is_none(), "{t:?} must be removed");
        }
        assert!(rec.data().get(&Tag::READ_GROUP).is_some());
        assert_eq!(counters.malformed_mod_reads.load(Ordering::Relaxed), 1);

        let missing = ubam_with_mods(b"CCCA", vec![40; 4], b"C+m,0,0,0;", vec![5, 6, 7]);
        let (out, _) = process_raw_full_window(raw_record(&missing), &cfg, &counters).unwrap();
        let Some(BamOutputRecord::Decoded(rec)) = out else {
            panic!("A missing MN forces a rebuild");
        };
        assert_eq!(
            rec.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH),
            Some(&Value::Int32(4))
        );
        assert_eq!(
            counters.malformed_mod_reads.load(Ordering::Relaxed),
            1,
            "A missing MN is not a defect"
        );

        let (out, _) =
            process_raw_full_window(raw_record(&read2_with_mods_and_rg()), &cfg, &counters)
                .unwrap();
        assert!(matches!(out, Some(BamOutputRecord::Raw(_))));
    }

    /// The move table resolves the signal end only to the stride, so a window
    /// that runs to the last base keeps the source `ns`.
    #[test]
    fn update_moves_window_to_the_last_base_keeps_the_source_ns() {
        let mut src = ubam_with_moves();
        // One sample past the last full stride block (10 + 8*2 = 26).
        src.data_mut()
            .insert(Tag::new(b'n', b's'), Value::Int32(27));
        let ns = |rec: &RecordBuf| rec.data().get(&Tag::new(b'n', b's')).cloned();

        let head = reconstruct_record(&src, 2, 6, 1, 0, true);
        assert_eq!(ns(&head), Some(Value::Int32(27)), "Head crop keeps ns");
        let tail = reconstruct_record(&src, 0, 4, 1, 0, true);
        assert_eq!(
            ns(&tail),
            Some(Value::Int32(22)),
            "Tail crop ends at a block"
        );
        // The last split segment spans [18, 27).
        let last = reconstruct_record(&src, 3, 6, 2, 1, true);
        assert_eq!(ns(&last), Some(Value::Int32(9)));
    }

    /// An empty window at the sequence end has no start base; the signal tags
    /// are dropped rather than the process aborted.
    #[test]
    fn update_moves_empty_window_at_the_end_drops_signal_tags() {
        let out = reconstruct_record(&ubam_with_moves(), 6, 6, 1, 0, true);
        assert!(out.sequence().as_ref().is_empty());
        for t in [b"mv", b"ts", b"ns"] {
            assert!(
                out.data().get(&Tag::new(t[0], t[1])).is_none(),
                "{} dropped",
                std::str::from_utf8(t).unwrap()
            );
        }
    }

    /// PacBio's fixed-size arrays are never per-base, even on a read whose
    /// length equals their element count.
    #[test]
    fn fixed_size_pacbio_arrays_are_not_sliced() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"ACGT".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        let d = src.data_mut();
        d.insert(
            Tag::new(b's', b'n'),
            Value::Array(Array::Float(vec![1.0, 2.0, 3.0, 4.0])),
        );
        d.insert(
            Tag::new(b'a', b'c'),
            Value::Array(Array::Int32(vec![0, 1, 1, 0])),
        );
        let out = reconstruct_record(&src, 1, 3, 1, 0, false);
        assert_eq!(
            out.data().get(&Tag::new(b's', b'n')),
            src.data().get(&Tag::new(b's', b'n'))
        );
        assert_eq!(
            out.data().get(&Tag::new(b'a', b'c')),
            src.data().get(&Tag::new(b'a', b'c'))
        );

        let mut two = RecordBuf::default();
        *two.flags_mut() = Flags::UNMAPPED;
        *two.name_mut() = Some(b"r2".into());
        *two.sequence_mut() = b"AC".to_vec().into();
        *two.quality_scores_mut() = vec![40; 2].into();
        two.data_mut().insert(
            Tag::new(b'b', b'c'),
            Value::Array(Array::UInt16(vec![3, 7])),
        );
        let cfg = cfg_bam2fq(None, 1, FastqTags::All);
        let mut out = Vec::new();
        run_bam_to_fastq(
            [Ok(two)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("\tbc:B:S,3,7\n"), "{s:?}");
    }

    /// A split product is not the sequencer's read, so `rn` is -1 on both
    /// outputs whether or not the move table is rewritten; a crop keeps it.
    #[test]
    fn split_sets_rn_to_minus_one_without_update_moves() {
        let mut src = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0;", vec![10]);
        src.data_mut().insert(Tag::new(b'r', b'n'), Value::Int32(7));
        let rn = |rec: &RecordBuf| rec.data().get(&Tag::new(b'r', b'n')).cloned();

        assert_eq!(
            rn(&reconstruct_record(&src, 0, 2, 2, 0, false)),
            Some(Value::Int32(-1))
        );
        assert_eq!(
            rn(&reconstruct_record(&src, 1, 4, 1, 0, false)),
            Some(Value::Int32(7))
        );

        let cfg = cfg_bam2fq(
            Some(QualityOp::Split {
                cutoff: 20,
                window: 1,
            }),
            0,
            FastqTags::All,
        );
        let mut out = Vec::new();
        let stats = run_bam_to_fastq(
            [Ok(src)].into_iter(),
            &mut out,
            &cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        assert_eq!(stats.output_reads, 2);
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("@r1_segment_1\trn:i:-1\t"), "{s:?}");
        assert!(s.contains("@r1_segment_2\trn:i:-1\t"), "{s:?}");
    }

    /// The rebuilt record keeps its aux tags in source order: rewritten tags
    /// stay in place, removed ones leave no hole, added ones are appended.
    #[test]
    fn rebuilt_record_keeps_aux_order() {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(b"r1".into());
        *src.sequence_mut() = b"CCAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 4].into();
        let d = src.data_mut();
        d.insert(Tag::READ_GROUP, Value::String(b"grp".as_slice().into()));
        d.insert(
            Tag::BASE_MODIFICATIONS,
            Value::String(b"C+m,0,1;".to_vec().into()),
        );
        d.insert(
            Tag::BASE_MODIFICATION_PROBABILITIES,
            Value::Array(Array::UInt8(vec![10, 20])),
        );
        d.insert(
            Tag::new(b'm', b'v'),
            Value::Array(Array::Int8(vec![5, 1, 1, 1, 1])),
        );
        d.insert(
            Tag::new(b'i', b'p'),
            Value::Array(Array::UInt8(vec![1, 2, 3, 4])),
        );
        d.insert(Tag::new(b'z', b'z'), Value::Int32(5));

        let out = reconstruct_record(&src, 2, 4, 1, 0, false);
        let order: Vec<[u8; 2]> = out.data().iter().map(|(t, _)| <[u8; 2]>::from(t)).collect();
        assert_eq!(
            order,
            [*b"RG", *b"MM", *b"ML", *b"ip", *b"zz", *b"MN"],
            "The mv tag is removed in place and MN appended"
        );
        assert_eq!(
            out.data().get(&Tag::new(b'z', b'z')),
            Some(&Value::Int32(5))
        );
    }

    const HIFI: &[u8] = b"m64011_190830_220126/123/ccs";

    /// A 10-base PacBio HiFi record: `qs`/`qe` query coordinates 100..110,
    /// `rn` passes, the `sa` coverage runs (5 over bases 0..4, 7 over 4..7,
    /// 5 over 7..10), the `sm`/`sx` per-base counts and the `ds`/`ls` undo
    /// blobs.
    fn pacbio_hifi_record(name: &[u8]) -> RecordBuf {
        let mut src = RecordBuf::default();
        *src.flags_mut() = Flags::UNMAPPED;
        *src.name_mut() = Some(name.into());
        *src.sequence_mut() = b"ACGTACGTAC".to_vec().into();
        *src.quality_scores_mut() = vec![40; 10].into();
        let d = src.data_mut();
        d.insert(Tag::new(b'q', b's'), Value::Int32(100));
        d.insert(Tag::new(b'q', b'e'), Value::Int32(110));
        d.insert(Tag::new(b'r', b'n'), Value::Int32(3));
        d.insert(
            Tag::new(b's', b'a'),
            Value::Array(Array::UInt32(vec![4, 5, 3, 7, 3, 5])),
        );
        d.insert(
            Tag::new(b's', b'm'),
            Value::Array(Array::UInt8((0..10).collect())),
        );
        d.insert(
            Tag::new(b's', b'x'),
            Value::Array(Array::UInt8((10..20).collect())),
        );
        d.insert(
            Tag::new(b'd', b's'),
            Value::Array(Array::UInt8(vec![1, 2, 3])),
        );
        d.insert(
            Tag::new(b'l', b's'),
            Value::Array(Array::UInt8(vec![4, 5, 6])),
        );
        src
    }

    /// A dorado record: `ubam_with_moves` (`mv`, `ts` 10, `ns` 26, `st`, `du`
    /// 5 s) plus `qs:f`, `me` and `er`, without `pi`.
    fn ont_record() -> RecordBuf {
        let mut src = ubam_with_moves();
        let d = src.data_mut();
        d.insert(Tag::new(b'q', b's'), Value::Float(20.0));
        d.insert(Tag::new(b'm', b'e'), Value::Int32(12));
        d.insert(
            Tag::new(b'e', b'r'),
            Value::String(b"signal_positive".as_slice().into()),
        );
        src
    }

    fn tag(rec: &RecordBuf, t: [u8; 2]) -> Option<Value> {
        rec.data().get(&Tag::new(t[0], t[1])).cloned()
    }

    fn string_value(s: &[u8]) -> Option<Value> {
        Some(Value::String(s.into()))
    }

    fn name_of(rec: &RecordBuf) -> Vec<u8> {
        rec.name().map(|n| n.to_vec()).unwrap_or_default()
    }

    /// Splits at a low-quality base and carries every tag.
    fn split_cfg() -> Config {
        cfg_bam2fq(
            Some(QualityOp::Split {
                cutoff: 20,
                window: 1,
            }),
            0,
            FastqTags::All,
        )
    }

    /// Runs the BAM-to-FASTQ workflow and returns the stats and the text.
    fn bam2fq(recs: Vec<RecordBuf>, cfg: &Config) -> (Stats, String) {
        let mut out = Vec::new();
        let stats = run_bam_to_fastq(
            recs.into_iter().map(anyhow::Ok),
            &mut out,
            cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        (stats, String::from_utf8(out).unwrap())
    }

    /// Runs the BAM-to-BAM workflow and returns the stats and the decoded
    /// output records.
    fn bam2bam(recs: Vec<RecordBuf>, cfg: &Config) -> (Stats, Vec<RecordBuf>) {
        let header = sam::Header::default();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("o.bam");
        let mut sink = crate::io::bam::writer(Some(&path), &header, 1, 6).unwrap();
        let stats = run_bam(
            &header,
            recs.into_iter().map(anyhow::Ok),
            &mut sink,
            cfg,
            &Arc::new(Counters::default()),
        )
        .unwrap();
        sink.finish().unwrap();
        let bytes = std::fs::read(&path).unwrap();
        let mut reader = noodles_bam::io::Reader::new(bytes.as_slice());
        let h = reader.read_header().unwrap();
        let mut out = Vec::new();
        let mut buf = RecordBuf::default();
        while reader.read_record_buf(&h, &mut buf).unwrap() != 0 {
            out.push(buf.clone());
        }
        (stats, out)
    }

    /// An integer `qs` or a PacBio read name classifies a record as PacBio;
    /// anything else is ONT.
    #[test]
    fn platform_follows_the_qs_type_and_the_read_name() {
        assert_eq!(platform(&pacbio_hifi_record(HIFI)), Platform::PacBio);
        assert_eq!(
            platform(&pacbio_hifi_record(b"hifi_1")),
            Platform::PacBio,
            "An integer qs alone classifies"
        );
        assert_eq!(platform(&ont_record()), Platform::Ont);
        let named = |name: &[u8]| {
            let mut rec = RecordBuf::default();
            *rec.name_mut() = Some(name.into());
            platform(&rec)
        };
        let pacbio: [&[u8]; 6] = [
            b"m1/123/ccs",
            b"m1/123/ccs/fwd",
            b"m1/123/ccs/rev",
            b"m1/123/ccs/100_200",
            b"m1/123/ccs/fwd/100_200",
            b"m1/123/100_200",
        ];
        for name in pacbio {
            assert_eq!(
                named(name),
                Platform::PacBio,
                "{}",
                String::from_utf8_lossy(name)
            );
        }
        let ont: [&[u8]; 9] = [
            b"r1",
            b"m1/abc/ccs",
            b"m1/123",
            b"m1/123/ccs/x",
            b"/123/ccs",
            b"m1/123/100_",
            b"m1/123/ccs/fwd/rev",
            b"m1/123/ccs/100_200/fwd",
            b"0123e4-uuid_segment_1",
        ];
        for name in ont {
            assert_eq!(
                named(name),
                Platform::Ont,
                "{}",
                String::from_utf8_lossy(name)
            );
        }
    }

    /// One integer coordinate without the other leaves both unchanged.
    #[test]
    fn lone_integer_qs_is_left_unchanged() {
        let mut src = pacbio_hifi_record(HIFI);
        src.data_mut().remove(&Tag::new(b'q', b'e'));
        let out = reconstruct_record(&src, 2, 7, 1, 0, false);
        assert_eq!(tag(&out, *b"qs"), Some(Value::Int32(100)));
        assert!(tag(&out, *b"qe").is_none());
    }

    /// A crop keeps the name; a split names PacBio segments by their query
    /// coordinates in the spec's forms and ONT segments `_segment_N`.
    #[test]
    fn split_names_follow_the_platform() {
        let w = |start, end, idx| Window {
            start,
            end,
            idx,
            total: 2,
        };
        let coords = Some((100, 104));
        let crop = Window {
            start: 0,
            end: 4,
            idx: 0,
            total: 1,
        };
        assert_eq!(segment_name(Platform::PacBio, HIFI, crop, coords), HIFI);
        let cases: [(&[u8], &[u8]); 5] = [
            (b"m1/123/ccs", b"m1/123/ccs/100_104"),
            (b"m1/123/ccs/100_110", b"m1/123/ccs/100_104"),
            (b"m1/123/ccs/fwd", b"m1/123/ccs/fwd/100_104"),
            (b"m1/123/ccs/rev/100_110", b"m1/123/ccs/rev/100_104"),
            (b"m1/123/100_110", b"m1/123/100_104"),
        ];
        for (name, want) in cases {
            assert_eq!(
                segment_name(Platform::PacBio, name, w(0, 4, 0), coords),
                want,
                "{}",
                String::from_utf8_lossy(name)
            );
        }
        // Without `qs`/`qe` the interval is offset from the name's own start.
        assert_eq!(
            segment_name(Platform::PacBio, b"m1/123/ccs/50_60", w(6, 10, 1), None),
            b"m1/123/ccs/56_60"
        );
        assert_eq!(
            segment_name(Platform::PacBio, b"m1/123/ccs", w(6, 10, 1), None),
            b"m1/123/ccs/6_10"
        );
        // A PacBio record without a PacBio name, and any ONT record, take the
        // suffix.
        assert_eq!(
            segment_name(Platform::PacBio, b"hifi_1", w(0, 4, 0), coords),
            b"hifi_1_segment_1"
        );
        assert_eq!(
            segment_name(Platform::Ont, b"m1/123/ccs", w(6, 10, 1), None),
            b"m1/123/ccs_segment_2"
        );
        let out = reconstruct_record(&pacbio_hifi_record(b"hifi_1"), 0, 4, 2, 0, false);
        assert_eq!(name_of(&out), b"hifi_1_segment_1");
    }

    /// Both output paths name a split PacBio read by its coordinates and
    /// rewrite `qs`/`qe` per segment.
    #[test]
    fn pacbio_split_is_named_and_coordinated_on_both_outputs() {
        let mut src = pacbio_hifi_record(HIFI);
        // Base 2 is low quality, so the split is [0,2) and [3,10).
        let mut quals = vec![40; 10];
        quals[2] = 1;
        *src.quality_scores_mut() = quals.into();
        let cfg = split_cfg();
        let (stats, recs) = bam2bam(vec![src.clone()], &cfg);
        assert_eq!(stats.output_reads, 2);
        assert_eq!(name_of(&recs[0]), b"m64011_190830_220126/123/ccs/100_102");
        assert_eq!(name_of(&recs[1]), b"m64011_190830_220126/123/ccs/103_110");
        assert_eq!(tag(&recs[1], *b"qs"), Some(Value::Int32(103)));
        assert_eq!(tag(&recs[1], *b"qe"), Some(Value::Int32(110)));
        let (_, s) = bam2fq(vec![src], &cfg);
        assert!(
            s.contains("@m64011_190830_220126/123/ccs/100_102\tqs:i:100\tqe:i:102\t"),
            "{s:?}"
        );
        assert!(
            s.contains("@m64011_190830_220126/123/ccs/103_110\tqs:i:103\tqe:i:110\t"),
            "{s:?}"
        );
    }

    /// `rn` is a pass count on PacBio and passes through a split, as does
    /// pbmarkdup's `du:Z`; an ONT split gets `rn` -1 and loses its float `du`.
    #[test]
    fn rn_and_du_follow_the_platform_on_a_split() {
        let mut src = pacbio_hifi_record(HIFI);
        src.data_mut().insert(
            Tag::new(b'd', b'u'),
            Value::String(b"dup-info".as_slice().into()),
        );
        let out = reconstruct_record(&src, 0, 4, 2, 0, false);
        assert_eq!(tag(&out, *b"rn"), Some(Value::Int32(3)));
        assert_eq!(tag(&out, *b"du"), string_value(b"dup-info"));
        assert!(tag(&out, *b"pi").is_none(), "pi is a dorado tag");
        let (_, s) = bam2fq(vec![src], &split_cfg());
        assert!(
            s.contains("\trn:i:3\t") && s.contains("\tdu:Z:dup-info"),
            "{s:?}"
        );

        let ont = reconstruct_record(&ont_record(), 0, 3, 2, 0, false);
        assert_eq!(tag(&ont, *b"rn"), Some(Value::Int32(-1)));
        assert!(tag(&ont, *b"du").is_none());
        assert!(tag(&ont, *b"st").is_none());
    }

    /// `sa` coverage runs are re-encoded for the window, `sm`/`sx` are sliced
    /// per base, and dorado's scalar `sm:f` is untouched.
    #[test]
    fn sa_coverage_is_resliced_and_sm_sx_are_per_base() {
        let src = pacbio_hifi_record(HIFI);
        // Per-base coverage [5,5,5,5,7,7,7,5,5,5]; window [2,7) keeps
        // [5,5,7,7,7].
        let out = reconstruct_record(&src, 2, 7, 1, 0, false);
        assert_eq!(
            tag(&out, *b"sa"),
            Some(Value::Array(Array::UInt32(vec![2, 5, 3, 7])))
        );
        assert_eq!(u8_array(&out, *b"sm"), (2..7).collect::<Vec<u8>>());
        assert_eq!(u8_array(&out, *b"sx"), (12..17).collect::<Vec<u8>>());
        // Window [3,8) keeps [5,7,7,7,5]: equal coverage on either side of a
        // run stays two runs.
        let out = reconstruct_record(&src, 3, 8, 1, 0, false);
        assert_eq!(
            tag(&out, *b"sa"),
            Some(Value::Array(Array::UInt32(vec![1, 5, 3, 7, 1, 5])))
        );
        let (_, s) = bam2fq(vec![src], &cfg_bam2fq(None, 2, FastqTags::All));
        assert!(s.contains("\tsa:B:I,2,5,3,7,3,5\t"), "{s:?}");
        assert!(s.contains("\tsm:B:C,2,3,4,5,6,7,8,9\t"), "{s:?}");

        let mut ont = ont_record();
        ont.data_mut()
            .insert(Tag::new(b's', b'm'), Value::Float(1.5));
        assert!(!has_malformed_perbase_tag(&ont, 6));
        let out = reconstruct_record(&ont, 2, 6, 1, 0, false);
        assert_eq!(tag(&out, *b"sm"), Some(Value::Float(1.5)));
    }

    /// An `sa` whose runs do not sum to the read length cannot be sliced: it
    /// is left unchanged and the read is counted as carrying a malformed
    /// per-base tag, on the decoded and the raw path alike.
    #[test]
    fn malformed_sa_is_left_unchanged_and_counted() {
        let mut src = pacbio_hifi_record(HIFI);
        src.data_mut().insert(
            Tag::new(b's', b'a'),
            Value::Array(Array::UInt32(vec![4, 5, 3, 7])),
        );
        assert!(has_malformed_perbase_tag(&src, 10));
        let out = reconstruct_record(&src, 2, 7, 1, 0, false);
        assert_eq!(
            tag(&out, *b"sa"),
            Some(Value::Array(Array::UInt32(vec![4, 5, 3, 7])))
        );
        let (_, malformed) = raw_full_window_metadata(&raw_record(&src)).unwrap();
        assert!(malformed);
        let (stats, _) = bam2fq(vec![src], &cfg_bam2fq(None, 2, FastqTags::All));
        assert_eq!(stats.malformed_tag_reads, 1);

        let mut odd = pacbio_hifi_record(HIFI);
        odd.data_mut().insert(
            Tag::new(b's', b'a'),
            Value::Array(Array::UInt32(vec![10, 5, 3])),
        );
        assert!(has_malformed_perbase_tag(&odd, 10));

        let ok = pacbio_hifi_record(HIFI);
        assert!(!has_malformed_perbase_tag(&ok, 10));
        let (_, malformed) = raw_full_window_metadata(&raw_record(&ok)).unwrap();
        assert!(!malformed);
    }

    /// The undo blobs leave every output record of a trimmed read and the read
    /// is counted once; an untrimmed read keeps them.
    #[test]
    fn undo_blobs_are_dropped_on_trim_and_counted_once_per_read() {
        let src = pacbio_hifi_record(HIFI);
        let out = reconstruct_record(&src, 2, 7, 1, 0, false);
        assert!(tag(&out, *b"ds").is_none() && tag(&out, *b"ls").is_none());
        let kept = reconstruct_record(&src, 0, 10, 1, 0, false);
        assert!(tag(&kept, *b"ds").is_some() && tag(&kept, *b"ls").is_some());

        // A cropped read with the blobs, one without, and a split read with
        // them: two reads counted, no output record carries them.
        let mut split = pacbio_hifi_record(b"m1/7/ccs");
        let mut quals = vec![40; 10];
        quals[5] = 1;
        *split.quality_scores_mut() = quals.into();
        let mut plain = pacbio_hifi_record(b"m1/8/ccs");
        plain.data_mut().remove(&Tag::new(b'd', b's'));
        plain.data_mut().remove(&Tag::new(b'l', b's'));
        let mut cfg = split_cfg();
        cfg.trim.head = 1;
        let reads = vec![src.clone(), plain, split];
        let (stats, recs) = bam2bam(reads.clone(), &cfg);
        assert_eq!(stats.output_reads, 4);
        assert_eq!(stats.undo_tags_dropped_reads, 2);
        assert!(
            recs.iter()
                .all(|r| tag(r, *b"ds").is_none() && tag(r, *b"ls").is_none())
        );
        let (stats, s) = bam2fq(reads, &cfg);
        assert_eq!(stats.undo_tags_dropped_reads, 2);
        assert!(!s.contains("ds:B") && !s.contains("ls:B"), "{s:?}");

        let (stats, recs) = bam2bam(vec![src], &cfg_bam2fq(None, 0, FastqTags::All));
        assert_eq!(stats.undo_tags_dropped_reads, 0);
        assert!(tag(&recs[0], *b"ds").is_some() && tag(&recs[0], *b"ls").is_some());
    }

    /// An ONT split carries the parent id, zero MinKNOW events and an unknown
    /// end reason on every subread but the last, with or without
    /// `--update-moves`; `sp` still needs the move table and `st`/`du` are
    /// dropped without it.
    #[test]
    fn ont_split_marks_subreads_without_update_moves() {
        let src = ont_record();
        let first = reconstruct_record(&src, 0, 3, 2, 0, false);
        assert_eq!(tag(&first, *b"pi"), string_value(b"r1"));
        assert_eq!(tag(&first, *b"me"), Some(Value::Int32(0)));
        assert_eq!(tag(&first, *b"er"), string_value(b"unknown"));
        assert!(tag(&first, *b"sp").is_none());
        assert!(tag(&first, *b"st").is_none() && tag(&first, *b"du").is_none());
        let last = reconstruct_record(&src, 3, 6, 2, 1, false);
        assert_eq!(tag(&last, *b"pi"), string_value(b"r1"));
        assert_eq!(tag(&last, *b"me"), Some(Value::Int32(0)));
        assert_eq!(tag(&last, *b"er"), string_value(b"signal_positive"));

        // A crop is still the sequencer's read.
        let crop = reconstruct_record(&src, 2, 6, 1, 0, false);
        assert!(tag(&crop, *b"pi").is_none());
        assert_eq!(tag(&crop, *b"me"), Some(Value::Int32(12)));
        assert_eq!(tag(&crop, *b"er"), string_value(b"signal_positive"));

        // A record without `me`/`er` gains neither; `pi` is always set.
        let bare = reconstruct_record(&ubam_with_moves(), 0, 3, 2, 0, false);
        assert!(tag(&bare, *b"me").is_none() && tag(&bare, *b"er").is_none());
        assert_eq!(tag(&bare, *b"pi"), string_value(b"r1"));

        let mut fq = src;
        *fq.quality_scores_mut() = vec![40, 40, 1, 40, 40, 40].into();
        let (stats, s) = bam2fq(vec![fq], &split_cfg());
        assert_eq!(stats.output_reads, 2);
        assert!(
            s.contains("@r1_segment_1\tqs:f:") && s.contains("@r1_segment_2\tqs:f:"),
            "{s:?}"
        );
        assert!(
            s.contains("\tme:i:0\ter:Z:unknown\tpi:Z:r1\trn:i:-1\n"),
            "{s:?}"
        );
        assert!(
            s.contains("\tme:i:0\ter:Z:signal_positive\tpi:Z:r1\trn:i:-1\n"),
            "{s:?}"
        );
        assert!(!s.contains("\tst:Z:") && !s.contains("\tdu:f:"), "{s:?}");
    }

    /// Under `--update-moves` a split recomputes `du` from the subread's
    /// samples and `st` from its start sample, at the rate `ns`/`du` gives
    /// (26 samples over 5 s), in the source's offset form.
    #[test]
    fn update_moves_split_recomputes_st_and_du() {
        // Segment [0,3) spans samples [10,18): 8 samples at 5.2/s last 1.538 s
        // and start 1.923 s into the read.
        let first = reconstruct_record(&ont_record(), 0, 3, 2, 0, true);
        match tag(&first, *b"du") {
            Some(Value::Float(d)) => assert!((d - 1.538_461_5).abs() < 1e-5, "{d}"),
            other => panic!("Unexpected du: {other:?}"),
        }
        assert_eq!(
            tag(&first, *b"st"),
            string_value(b"2024-06-21T10:00:01.923Z")
        );
        // Segment [3,6) spans [18,26): 3.462 s in.
        let last = reconstruct_record(&ont_record(), 3, 6, 2, 1, true);
        assert_eq!(
            tag(&last, *b"st"),
            string_value(b"2024-06-21T10:00:03.462Z")
        );
        assert_eq!(tag(&last, *b"sp"), Some(Value::Int32(18)));
        // A crop keeps both.
        let crop = reconstruct_record(&ont_record(), 2, 6, 1, 0, true);
        assert_eq!(tag(&crop, *b"st"), string_value(b"2024-06-21T10:00:00Z"));
        assert_eq!(tag(&crop, *b"du"), Some(Value::Float(5.0)));

        let with_st = |st: &[u8]| {
            let mut r = ont_record();
            r.data_mut()
                .insert(Tag::new(b's', b't'), Value::String(st.into()));
            tag(&reconstruct_record(&r, 0, 3, 2, 0, true), *b"st")
        };
        assert_eq!(
            with_st(b"2024-06-21T10:00:00.000+00:00"),
            string_value(b"2024-06-21T10:00:01.923+00:00")
        );
        assert_eq!(
            with_st(b"2024-06-21T12:00:00.500+02:00"),
            string_value(b"2024-06-21T12:00:02.423+02:00")
        );
        assert_eq!(
            with_st(b"2024-06-21T10:00:00"),
            string_value(b"2024-06-21T10:00:01.923")
        );
        // An `st` that does not parse is left as is.
        assert_eq!(with_st(b"yesterday"), string_value(b"yesterday"));
    }

    /// The recompute needs a positive float `du` and an `ns`; otherwise both
    /// tags are copied unchanged rather than dropped or guessed.
    #[test]
    fn update_moves_split_leaves_st_and_du_when_the_rate_is_unknown() {
        let variants: [(&[u8; 2], Option<Value>); 3] = [
            (b"du", Some(Value::Float(0.0))),
            (b"du", None),
            (b"ns", None),
        ];
        for (t, v) in variants {
            let mut r = ont_record();
            let t = Tag::new(t[0], t[1]);
            match v {
                Some(v) => {
                    r.data_mut().insert(t, v);
                },
                None => {
                    r.data_mut().remove(&t);
                },
            }
            let out = reconstruct_record(&r, 0, 3, 2, 0, true);
            assert_eq!(tag(&out, *b"st"), tag(&r, *b"st"), "{t:?}");
            assert_eq!(tag(&out, *b"du"), tag(&r, *b"du"), "{t:?}");
        }
    }
}
