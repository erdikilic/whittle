//! Base-modification tags of a record: classification of its MM/ML block and
//! the rebuild for an output window.

use super::*;

/// The base-modification block, in the order it is emitted.
pub(super) const MOD_TAGS: [Tag; 3] = [
    Tag::BASE_MODIFICATIONS,
    Tag::BASE_MODIFICATION_PROBABILITIES,
    Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
];

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
pub(super) fn classify_mod_block(
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
pub(super) fn mod_tags(src: &RecordBuf) -> Option<(&[u8], Option<&[u8]>)> {
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

/// Rebuilds the `MM`/`ML` block of a `Consistent` or `MissingMn` record for the
/// window `[start, end)`: skip-counts renumbered, `ML` re-sliced. `ml` is
/// `None` when the source carries no `ML`, and so is the result's.
pub(super) fn rebuild_mods(
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
