//! Output uBAM records assembled from the raw bytes of their input record: the
//! fixed fields and read name, the window's slice of the packed sequence and
//! qualities, and the aux fields copied verbatim, cut to the window or
//! re-encoded.
//!
//! The output is byte-identical to the noodles encoding of the decoded record
//! with the same edits applied, and the builder refuses the inputs that
//! encoder or the record decoder refuses.

use std::ops::Range;

use noodles_bam::record::codec::encoder::data::field::{ty, write_value};
use noodles_sam::alignment::record::data::field::Value as FieldValue;

use super::*;

/// The bin of a record without a position, `reg2bin(-1, 0)` (SAM spec
/// section 4.2.1).
const UNMAPPED_BIN: u16 = 4680;

/// The highest base quality a BAM record stores (SAM spec section 4.2.3).
const MAX_QUALITY: u8 = 93;

/// The longest read name a BAM record stores, NUL excluded.
const MAX_NAME_LENGTH: usize = 254;

/// The overflow store of a CIGAR longer than 65535 operations. The noodles
/// encoder does not write it from the aux fields of a decoded record, so a
/// rebuilt record omits it.
const CIGAR_TAG: [u8; 2] = *b"CG";

/// The length of a `B` array field's header: tag, type, subtype and count.
const ARRAY_HEADER_LEN: usize = 8;

/// The changes one output record makes to its input record.
pub(super) struct RecordEdit<'a> {
    /// First base of the window, inclusive.
    pub(super) start: usize,
    /// End of the window, exclusive.
    pub(super) end: usize,
    /// The output read name; `None` keeps the input name.
    pub(super) name: Option<Vec<u8>>,
    /// Tag rewrites: `Some` replaces the input value in place or is appended
    /// when the input lacks the tag, `None` removes the tag.
    pub(super) updates: TagUpdates,
    /// The tags `--remove-tag` drops, applied after the rewrites.
    pub(super) remove: &'a TagRemoval,
    /// The rejection reason, written as the `wr:Z` tag.
    pub(super) reason: Option<Reason>,
}

impl<'a> RecordEdit<'a> {
    /// The edit that keeps every base and name of a `len`-base record and
    /// drops only the tags `remove` names.
    pub(super) fn unchanged(len: usize, remove: &'a TagRemoval) -> Self {
        Self {
            start: 0,
            end: len,
            name: None,
            updates: Vec::new(),
            remove,
            reason: None,
        }
    }
}

/// The source of one output aux field's bytes.
enum FieldBytes {
    /// A field of the input aux block, copied verbatim.
    Source(Range<usize>),
    /// A `B` array of the input aux block cut to the window: the rewritten
    /// field header and the byte range of the kept elements.
    Sliced([u8; ARRAY_HEADER_LEN], Range<usize>),
    /// A field encoded into the scratch buffer.
    Encoded(Range<usize>),
}

/// One aux field of the input aux block.
struct RawField {
    /// The field's tag.
    tag: [u8; 2],
    /// The field's type code.
    ty: u8,
    /// The whole field, tag included.
    bytes: Range<usize>,
    /// The element size and count of a `B` array.
    array: Option<(usize, usize)>,
}

fn invalid_input(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidInput, message.to_string())
}

fn invalid_data(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_string())
}

/// Locates the aux field that starts at byte `pos` of `aux`. A truncated
/// field, an unknown type or array subtype, or a string without its NUL is
/// refused, as the noodles record decoder refuses it.
fn next_field(aux: &[u8], pos: usize) -> io::Result<RawField> {
    let &[t0, t1, ty, ..] = &aux[pos..] else {
        return Err(invalid_data("truncated aux field"));
    };
    let value = pos + 3;
    let (end, array) = match ty {
        b'A' | b'c' | b'C' => (value + 1, None),
        b's' | b'S' => (value + 2, None),
        b'i' | b'I' | b'f' => (value + 4, None),
        b'Z' | b'H' => {
            let nul = aux[value..]
                .iter()
                .position(|&b| b == 0)
                .ok_or_else(|| invalid_data("aux string is not NUL terminated"))?;
            (value + nul + 1, None)
        },
        b'B' => {
            let &[subtype, c0, c1, c2, c3, ..] = &aux[value.min(aux.len())..] else {
                return Err(invalid_data("truncated aux array"));
            };
            let size = match subtype {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return Err(invalid_data("invalid aux array subtype")),
            };
            let count = usize::try_from(u32::from_le_bytes([c0, c1, c2, c3]))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let end = count
                .checked_mul(size)
                .and_then(|n| n.checked_add(value + 5))
                .ok_or_else(|| invalid_data("truncated aux array"))?;
            (end, Some((size, count)))
        },
        _ => return Err(invalid_data("invalid aux type")),
    };
    if end > aux.len() {
        return Err(invalid_data("truncated aux field"));
    }
    Ok(RawField {
        tag: [t0, t1],
        ty,
        bytes: pos..end,
        array,
    })
}

/// Checks a copied `Z` or `H` value, NUL excluded, as the noodles encoder
/// does: printable ASCII for a string, an even count of uppercase hex digits
/// for a hex value.
fn check_text(ty: u8, value: &[u8]) -> io::Result<()> {
    match ty {
        b'Z' if !value.iter().all(|b| matches!(b, b' '..=b'~')) => {
            Err(invalid_input("invalid aux string"))
        },
        b'H' if !value.len().is_multiple_of(2)
            || !value.iter().all(|b| matches!(b, b'0'..=b'9' | b'A'..=b'F')) =>
        {
            Err(invalid_input("invalid aux hex value"))
        },
        _ => Ok(()),
    }
}

/// Returns whether a per-base array under `tag` is cut to the window by the
/// structural rule: a tag without a dedicated rule, other than the `sa`
/// coverage runs and the fixed-size PacBio arrays.
fn slices_per_base(tag: [u8; 2]) -> bool {
    !has_dedicated_rule(tag) && tag != RLE_COVERAGE_TAG && !FIXED_ARRAY_TAGS.contains(&tag)
}

/// Appends a `B` array value, type code included.
fn push_array(dst: &mut Vec<u8>, array: &Array) -> io::Result<()> {
    fn header(dst: &mut Vec<u8>, subtype: u8, len: usize) -> io::Result<()> {
        let count =
            u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        dst.extend_from_slice(&[b'B', subtype]);
        dst.extend_from_slice(&count.to_le_bytes());
        Ok(())
    }
    match array {
        Array::Int8(v) => {
            header(dst, b'c', v.len())?;
            dst.extend(v.iter().map(|&n| n as u8));
        },
        Array::UInt8(v) => {
            header(dst, b'C', v.len())?;
            dst.extend_from_slice(v);
        },
        Array::Int16(v) => {
            header(dst, b's', v.len())?;
            dst.extend(v.iter().flat_map(|n| n.to_le_bytes()));
        },
        Array::UInt16(v) => {
            header(dst, b'S', v.len())?;
            dst.extend(v.iter().flat_map(|n| n.to_le_bytes()));
        },
        Array::Int32(v) => {
            header(dst, b'i', v.len())?;
            dst.extend(v.iter().flat_map(|n| n.to_le_bytes()));
        },
        Array::UInt32(v) => {
            header(dst, b'I', v.len())?;
            dst.extend(v.iter().flat_map(|n| n.to_le_bytes()));
        },
        Array::Float(v) => {
            header(dst, b'f', v.len())?;
            dst.extend(v.iter().flat_map(|n| n.to_le_bytes()));
        },
    }
    Ok(())
}

/// Encodes one aux field into `scratch` and returns its byte range. Arrays are
/// written in bulk; every other value goes through the noodles value encoder,
/// which also validates strings and hex values.
fn encode_field(scratch: &mut Vec<u8>, tag: [u8; 2], value: &Value) -> io::Result<Range<usize>> {
    let begin = scratch.len();
    scratch.extend_from_slice(&tag);
    match value {
        Value::Array(array) => push_array(scratch, array)?,
        _ => {
            let value = FieldValue::from(value);
            scratch.push(ty::encode(value.ty()));
            write_value(scratch, &value)?;
        },
    }
    Ok(begin..scratch.len())
}

/// Places `bytes` under `tag`: in place of a field already written under
/// `tag`, otherwise after the last field. This is the insertion rule of a
/// decoded record's aux data.
fn insert(fields: &mut Vec<([u8; 2], FieldBytes)>, tag: [u8; 2], bytes: FieldBytes) {
    match fields.iter_mut().find(|(t, _)| *t == tag) {
        Some(slot) => slot.1 = bytes,
        None => fields.push((tag, bytes)),
    }
}

/// Collects the output aux fields of `edit` over the input aux block `aux`:
/// input fields in order, each removed, replaced by its update, cut to the
/// window (a per-base array of a trimmed read) or copied verbatim, then the
/// updates the input lacks and the rejection reason.
fn edit_fields(
    aux: &[u8],
    orig_len: usize,
    edit: &mut RecordEdit<'_>,
    scratch: &mut Vec<u8>,
) -> io::Result<Vec<([u8; 2], FieldBytes)>> {
    let (start, end) = (edit.start, edit.end);
    let trimmed = start != 0 || end != orig_len;
    let remove = edit.remove;
    let updates = &mut edit.updates;
    if !remove.is_empty() {
        updates.retain(|(t, _)| !remove.contains(&<[u8; 2]>::from(*t)));
    }

    let mut fields = Vec::new();
    let mut seen: Vec<[u8; 2]> = Vec::new();
    let mut pos = 0;
    while pos < aux.len() {
        let field = next_field(aux, pos)?;
        pos = field.bytes.end;
        let tag = field.tag;
        if seen.contains(&tag) {
            return Err(invalid_data("duplicate aux tag"));
        }
        seen.push(tag);
        if tag == CIGAR_TAG || remove.contains(&tag) {
            continue;
        }
        if let Some(i) = updates.iter().position(|(t, _)| <[u8; 2]>::from(*t) == tag) {
            if let (_, Some(value)) = updates.remove(i) {
                let bytes = encode_field(scratch, tag, &value)?;
                insert(&mut fields, tag, FieldBytes::Encoded(bytes));
            }
            continue;
        }
        let bytes = match field.array {
            Some((size, count)) if trimmed && count == orig_len && slices_per_base(tag) => {
                // Reverse-strand kinetics are stored last base first.
                let (s, e) = if REVERSED_PERBASE_TAGS.contains(&tag) {
                    (orig_len - end, orig_len - start)
                } else {
                    (start, end)
                };
                let first = field.bytes.start;
                let elements = first + ARRAY_HEADER_LEN;
                let mut head = [0; ARRAY_HEADER_LEN];
                head[..4].copy_from_slice(&aux[first..first + 4]);
                // `e - s` is at most the source count, which fits a `u32`.
                head[4..].copy_from_slice(&((e - s) as u32).to_le_bytes());
                FieldBytes::Sliced(head, elements + s * size..elements + e * size)
            },
            _ => {
                check_text(
                    field.ty,
                    &aux[field.bytes.start + 3..field.bytes.end.saturating_sub(1)],
                )?;
                FieldBytes::Source(field.bytes)
            },
        };
        insert(&mut fields, tag, bytes);
    }
    for (tag, value) in updates.drain(..) {
        if let Some(value) = value {
            let tag = <[u8; 2]>::from(tag);
            let bytes = encode_field(scratch, tag, &value)?;
            insert(&mut fields, tag, FieldBytes::Encoded(bytes));
        }
    }
    if let Some(reason) = edit.reason {
        let tag = <[u8; 2]>::from(reject::REASON_TAG);
        let value = Value::String(reason.label().into());
        let bytes = encode_field(scratch, tag, &value)?;
        insert(&mut fields, tag, FieldBytes::Encoded(bytes));
    }
    Ok(fields)
}

/// Checks a read name as the noodles encoder does: 1 to 254 printable ASCII
/// characters other than `@`, and not the missing-name marker `*`.
fn check_name(name: &[u8]) -> io::Result<()> {
    let valid = (1..=MAX_NAME_LENGTH).contains(&name.len())
        && name != b"*"
        && name.iter().all(|&b| b.is_ascii_graphic() && b != b'@');
    if valid {
        Ok(())
    } else {
        Err(invalid_input("invalid read name"))
    }
}

/// Returns the BAM bin of a record placed at 1-based `start` whose alignment
/// covers `span` reference bases (SAM spec section 5.3, `reg2bin`): the
/// smallest bin level whose 16 kbp to 64 Mbp window holds both ends. An empty
/// span ends the record at its start, and a bin past `u16` is truncated, as
/// the noodles encoder computes it.
fn region_bin(start: usize, span: usize) -> u16 {
    /// The shift of each bin level's window size and the level's first bin.
    const LEVELS: [(u32, usize); 5] = [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)];
    let beg = start - 1;
    let end = if span == 0 { beg } else { beg + span - 1 };
    let bin = LEVELS
        .iter()
        .find(|&&(shift, _)| beg >> shift == end >> shift)
        .map_or(0, |&(shift, first)| first + (beg >> shift));
    bin as u16
}

/// Returns the reference span and the read length of a packed CIGAR, refusing
/// an unknown operation.
fn cigar_lengths(cigar: &[u8]) -> io::Result<(usize, usize)> {
    let mut span = 0;
    let mut read_len = 0;
    for op in cigar.as_chunks::<4>().0 {
        let n = u32::from_le_bytes(*op);
        let len = (n >> 4) as usize;
        match n & 0x0f {
            // M, =, X
            0 | 7 | 8 => {
                span += len;
                read_len += len;
            },
            // D, N
            2 | 3 => span += len,
            // I, S
            1 | 4 => read_len += len,
            // H, P
            5 | 6 => {},
            _ => return Err(invalid_data("invalid CIGAR operation")),
        }
    }
    Ok((span, read_len))
}

/// Appends bases `[start, end)` of a 4-bit packed sequence, which stores the
/// first base of each byte in its high nibble. An odd `start` shifts every
/// base by one nibble; an odd base count leaves the low nibble of the last
/// byte zero, as the SAM spec recommends.
pub(super) fn push_packed_bases(dst: &mut Vec<u8>, packed: &[u8], start: usize, end: usize) {
    let n = end - start;
    if n == 0 {
        return;
    }
    let src = &packed[start / 2..end.div_ceil(2)];
    if start.is_multiple_of(2) {
        dst.extend_from_slice(src);
        if !n.is_multiple_of(2)
            && let Some(last) = dst.last_mut()
        {
            *last &= 0xf0;
        }
    } else {
        dst.extend(src.windows(2).map(|w| (w[0] << 4) | (w[1] >> 4)));
        if !n.is_multiple_of(2) {
            dst.push(src[src.len() - 1] << 4);
        }
    }
}

/// Returns the BAM encoding of an optional 0-based field: -1 when absent.
fn encoded_index(index: Option<usize>) -> io::Result<i32> {
    index.map_or(Ok(-1), |i| {
        i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })
}

/// Builds the output record `edit` describes from the raw input record `src`,
/// without the `block_size` prefix: the fixed fields with the window's lengths
/// and a recomputed bin, the output name, the input CIGAR, the window's bases
/// and qualities, and the edited aux fields (`edit_fields`).
///
/// The reference sequence IDs are checked against the output header when the
/// record is written (`io::bam::write_record_bytes`).
pub(super) fn build_record(src: &bam::Record, mut edit: RecordEdit<'_>) -> io::Result<Vec<u8>> {
    let sequence = src.sequence();
    let orig_len = sequence.len();
    let (start, end) = (edit.start, edit.end);
    if start > end || end > orig_len {
        return Err(invalid_input("output window exceeds the read"));
    }
    let base_count = end - start;

    let data = src.data();
    let aux = data.as_bytes();
    let mut scratch = Vec::new();
    let fields = edit_fields(aux, orig_len, &mut edit, &mut scratch)?;

    let name: &[u8] = match (&edit.name, src.name()) {
        (Some(name), _) => name,
        (None, Some(name)) => name,
        (None, None) => b"*",
    };
    if edit.name.is_some() || src.name().is_some() {
        check_name(name)?;
    }

    let cigar = src.cigar();
    let cigar = cigar.as_bytes();
    let (span, read_len) = cigar_lengths(cigar)?;
    if base_count > 0 && read_len > 0 && read_len != base_count {
        return Err(invalid_input("read length-sequence length mismatch"));
    }

    let qual = &src.quality_scores().as_bytes()[start..end];
    if qual.iter().any(|&q| q > MAX_QUALITY) {
        return Err(invalid_input("invalid base quality"));
    }

    let alignment_start = src.alignment_start().transpose()?.map(usize::from);
    let bin = alignment_start.map_or(UNMAPPED_BIN, |start| region_bin(start, span));
    let aux_len: usize = fields
        .iter()
        .map(|(_, bytes)| match bytes {
            FieldBytes::Source(r) | FieldBytes::Encoded(r) => r.len(),
            FieldBytes::Sliced(head, r) => head.len() + r.len(),
        })
        .sum();

    let mut dst = Vec::with_capacity(
        32 + name.len() + 1 + cigar.len() + base_count.div_ceil(2) + base_count + aux_len,
    );
    let reference_id = encoded_index(src.reference_sequence_id().transpose()?)?;
    dst.extend_from_slice(&reference_id.to_le_bytes());
    let position = encoded_index(alignment_start.map(|p| p - 1))?;
    dst.extend_from_slice(&position.to_le_bytes());
    // `check_name` bounds the name to 254 bytes, so the length with its NUL
    // fits a byte.
    dst.push((name.len() + 1) as u8);
    dst.push(src.mapping_quality().map_or(0xff, u8::from));
    dst.extend_from_slice(&bin.to_le_bytes());
    // A packed CIGAR holds at most `u16::MAX` operations.
    dst.extend_from_slice(&((cigar.len() / 4) as u16).to_le_bytes());
    dst.extend_from_slice(&u16::from(src.flags()).to_le_bytes());
    // The window is no longer than the input sequence, whose length is a `u32`.
    dst.extend_from_slice(&(base_count as u32).to_le_bytes());
    let mate_reference_id = encoded_index(src.mate_reference_sequence_id().transpose()?)?;
    dst.extend_from_slice(&mate_reference_id.to_le_bytes());
    let mate_position = encoded_index(
        src.mate_alignment_start()
            .transpose()?
            .map(|p| usize::from(p) - 1),
    )?;
    dst.extend_from_slice(&mate_position.to_le_bytes());
    dst.extend_from_slice(&src.template_length().to_le_bytes());
    dst.extend_from_slice(name);
    dst.push(0);
    dst.extend_from_slice(cigar);
    push_packed_bases(&mut dst, sequence.as_bytes(), start, end);
    dst.extend_from_slice(qual);
    for (_, bytes) in &fields {
        match bytes {
            FieldBytes::Source(r) => dst.extend_from_slice(&aux[r.clone()]),
            FieldBytes::Sliced(head, r) => {
                dst.extend_from_slice(head);
                dst.extend_from_slice(&aux[r.clone()]);
            },
            FieldBytes::Encoded(r) => dst.extend_from_slice(&scratch[r.clone()]),
        }
    }
    Ok(dst)
}
