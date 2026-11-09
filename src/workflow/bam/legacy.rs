//! Conversion of legacy tag forms on input: the lowercase `Mm`/`Ml`
//! modification tags, an integer mean quality `qs` and HiFi kinetics that
//! index the whole CCS read.

use super::*;

/// The lowercase spellings of the base-modification tags, each with the tag
/// it stands for. htslib reads `Mm` in place of an absent `MM` and `Ml` in
/// place of an absent `ML` (`sam_mods.c`).
const LEGACY_MOD_TAGS: [([u8; 2], [u8; 2]); 2] = [(*b"Mm", *b"MM"), (*b"Ml", *b"ML")];

/// The HiFi kinetics arrays: `fi`/`fp` in base order, `ri`/`rp` last base
/// first (`REVERSED_PERBASE_TAGS`).
const HIFI_KINETICS_TAGS: [[u8; 2]; 4] = [*b"fi", *b"fp", *b"ri", *b"rp"];

/// The tags `upgrade_legacy_tags` can change.
const UPGRADED_TAGS: [[u8; 2]; 9] = [
    *b"Mm", *b"Ml", *b"MM", *b"ML", *b"qs", *b"fi", *b"fp", *b"ri", *b"rp",
];

/// Returns whether an integer `qs` is a mean read quality rather than a
/// query coordinate: the record carries no `qe` and its name follows no
/// PacBio convention (`parse_pacbio_name`). PacBio query coordinates are a
/// `qs`/`qe` pair.
fn ont_integer_qs(name: Option<&[u8]>, has_qe: bool) -> bool {
    !has_qe && !name.is_some_and(|n| parse_pacbio_name(n).is_some())
}

/// Returns the query interval `[qs, qe)` of a read of `len` bases, when it
/// spans exactly `len` bases.
fn kept_span(qs: i64, qe: i64, len: usize) -> Option<(usize, usize)> {
    let (qs, qe) = (usize::try_from(qs).ok()?, usize::try_from(qe).ok()?);
    (qe.checked_sub(qs)? == len).then_some((qs, qe))
}

/// Returns the entries of a HiFi kinetics array of `whole` values that hold
/// the bases `[qs, qe)` of a read of `len` bases, for an array that indexes
/// the whole CCS read: `[qs, qe)` in base order, `[whole - qe, whole - qs)`
/// last base first. `None` for an array of the read length or one too short
/// to hold the interval.
fn unclipped_window(
    tag: [u8; 2],
    whole: usize,
    (qs, qe): (usize, usize),
    len: usize,
) -> Option<(usize, usize)> {
    if whole == len || whole < qe {
        return None;
    }
    Some(if REVERSED_PERBASE_TAGS.contains(&tag) {
        (whole - qe, whole - qs)
    } else {
        (qs, qe)
    })
}

/// Converts the legacy tag forms of `rec`. `Mm`/`Ml` values move to an absent
/// `MM`/`ML` and are dropped beside a present one, as htslib reads them. An
/// integer mean quality `qs` (`ont_integer_qs`) becomes a float. HiFi
/// kinetics that index the whole CCS read are cut to the bases of the read
/// (`unclipped_window`). Returns whether `rec` held a legacy form.
pub(super) fn upgrade_legacy_tags(rec: &mut RecordBuf) -> bool {
    let qs = Tag::new(b'q', b's');
    let ont_qs = ont_integer_qs(
        rec.name().map(AsRef::as_ref),
        rec.data().get(&Tag::new(b'q', b'e')).is_some(),
    );
    let len = rec.sequence().len();
    let kept = query_span(rec).and_then(|(qs, qe)| kept_span(qs, qe, len));
    let data = rec.data_mut();
    let mut found = false;
    for (legacy, current) in LEGACY_MOD_TAGS {
        if let Some((_, value)) = data.remove(&Tag::new(legacy[0], legacy[1])) {
            found = true;
            let current = Tag::new(current[0], current[1]);
            if data.get(&current).is_none() {
                data.insert(current, value);
            }
        }
    }
    if ont_qs && let Some(mean) = data.get(&qs).and_then(|value| value.as_int()) {
        data.insert(qs, Value::Float(mean as f32));
        found = true;
    }
    if let Some(span) = kept {
        for t in HIFI_KINETICS_TAGS {
            let tag = Tag::new(t[0], t[1]);
            let cut = match data.get(&tag) {
                Some(Value::Array(array)) => unclipped_window(t, array_len(array), span, len)
                    .map(|(start, end)| slice_array(array, start, end)),
                _ => None,
            };
            if let Some(cut) = cut {
                data.insert(tag, Value::Array(cut));
                found = true;
            }
        }
    }
    found
}

/// Returns `record` with its legacy tag forms converted (`upgraded_tags`); a
/// record without one is returned as is.
pub(super) fn current_tags(
    record: bam::Record,
    counters: &Counters,
) -> anyhow::Result<bam::Record> {
    Ok(upgraded_tags(&record, counters)?.unwrap_or(record))
}

/// Rebuilds a raw record that carries a legacy tag form with the tags of
/// `upgrade_legacy_tags` and counts it in `Counters::legacy_tags_reads`. The
/// workflows process the rebuilt record as any other. `None` for a record
/// without a legacy form.
pub(super) fn upgraded_tags(
    record: &bam::Record,
    counters: &Counters,
) -> anyhow::Result<Option<bam::Record>> {
    let data = record.data();
    let has = |t: [u8; 2]| data.get(&Tag::new(t[0], t[1])).is_some();
    let integer_qs = data
        .get(&Tag::new(b'q', b's'))
        .and_then(Result::ok)
        .is_some_and(|value| raw_integer(&value).is_some())
        && ont_integer_qs(record.name().map(AsRef::as_ref), has(*b"qe"));
    let integer = |t: &[u8; 2]| {
        data.get(&Tag::new(t[0], t[1]))
            .and_then(Result::ok)
            .and_then(|value| raw_integer(&value))
    };
    let len = record.sequence().len();
    let unclipped = integer(b"qs")
        .zip(integer(b"qe"))
        .and_then(|(qs, qe)| kept_span(qs, qe, len))
        .is_some_and(|span| {
            HIFI_KINETICS_TAGS.into_iter().any(|t| {
                data.get(&Tag::new(t[0], t[1]))
                    .and_then(Result::ok)
                    .and_then(|value| raw_array_len(&value))
                    .is_some_and(|whole| unclipped_window(t, whole, span, len).is_some())
            })
        });
    if !integer_qs && !unclipped && !LEGACY_MOD_TAGS.iter().any(|&(t, _)| has(t)) {
        return Ok(None);
    }
    let mut decoded = decode_raw_record(record).with_context(read_context(record))?;
    upgrade_legacy_tags(&mut decoded);
    let updates = UPGRADED_TAGS
        .into_iter()
        .map(|t| {
            let tag = Tag::new(t[0], t[1]);
            (tag, decoded.data().get(&tag).cloned())
        })
        .collect();
    let keep = TagRemoval::default();
    let edit = RecordEdit {
        updates,
        ..RecordEdit::unchanged(record.sequence().len(), &keep)
    };
    let bytes = build_record(record, edit).with_context(read_context(record))?;
    counters.legacy_tags_reads.fetch_add(1, Ordering::Relaxed);
    let upgraded = raw_record_from_bytes(&bytes).with_context(read_context(record))?;
    Ok(Some(upgraded))
}
