//! The full-window path for BAM input that no step trims: records pass through
//! as raw bytes after metadata checks.

use super::*;

/// A record ready to write: the untouched raw input or a rebuilt record.
pub(super) enum BamOutputRecord {
    /// The raw input record, written without decoding.
    Raw(bam::Record),
    /// A rebuilt record's bytes, without the `block_size` prefix.
    Built(Vec<u8>),
}

pub(super) fn raw_array_len(
    value: &noodles_sam::alignment::record::data::field::Value<'_>,
) -> Option<usize> {
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
pub(super) fn raw_rle_coverage_len(
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
pub(super) fn raw_integer(
    value: &noodles_sam::alignment::record::data::field::Value<'_>,
) -> Option<i64> {
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
pub(super) fn raw_full_window_metadata(record: &bam::Record) -> std::io::Result<(ModBlock, bool)> {
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
pub(super) fn ensure_raw_trimmable(record: &bam::Record) -> anyhow::Result<()> {
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
pub(super) fn raw_gc_fraction(record: &bam::Record) -> f64 {
    let sequence = record.sequence();
    if sequence.is_empty() {
        return 0.0;
    }
    let gc = sequence.iter().filter(|&b| crate::filter::is_gc(b)).count();
    gc as f64 / sequence.len() as f64
}

/// Returns the edit of a raw record's full window that the decoded workflow
/// applies to its modification block (`window_edit`): `MN` added to a block
/// without it, a malformed block removed. `None` for any other block, whose
/// record is its own output.
fn full_window_mod_edit<'a>(
    record: &bam::Record,
    mod_block: ModBlock,
    remove: &'a TagRemoval,
) -> anyhow::Result<Option<RecordEdit<'a>>> {
    if !matches!(mod_block, ModBlock::MissingMn | ModBlock::Malformed) {
        return Ok(None);
    }
    let decoded = decode_raw_record(record).with_context(read_context(record))?;
    let window = Window {
        start: 0,
        end: record.sequence().len(),
        idx: 0,
        total: 1,
    };
    Ok(window_edit(&decoded, window, mod_block, None, None, remove))
}

/// Filters one raw record over its full window and decides its output: the raw
/// record itself when nothing changes, a rebuild when `MN` is missing or the
/// modification block is malformed, nothing when the filter drops it. A
/// rejected record gets the same rebuild, with its reason.
/// A malformed modification block or per-base tag is counted. The decision
/// is traced as `process_read_segments` traces it on the decoded path.
///
/// Tag removal never reaches here: `run_raw_bam` excludes it from the
/// full-window shortcut, so a run that removes tags rebuilds every record
/// through `run_bam`.
pub(super) fn process_raw_full_window(
    record: bam::Record,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Option<BamOutputRecord>> {
    let seq_len = record.sequence().len();
    let qualities = record.quality_scores();
    let qual = qualities.as_ref();
    require_qualities(seq_len, qual).map_err(|e| {
        e.context(format!(
            "read {}",
            crate::io::bam::display_name(record.name().map(AsRef::as_ref))
        ))
    })?;

    let (mod_block, malformed_perbase) =
        raw_full_window_metadata(&record).with_context(read_context(&record))?;
    if mod_block == ModBlock::Malformed {
        counters.malformed_mod_reads.fetch_add(1, Ordering::Relaxed);
    }
    if malformed_perbase {
        counters.malformed_tag_reads.fetch_add(1, Ordering::Relaxed);
    }

    let read = crate::workflow::read_span(record.name().map_or(b"<unnamed>", |n| n.as_ref()));
    let _read = read.enter();
    let rejected = if seq_len == 0 {
        tracing::trace!("Read produced no segments");
        counters
            .reads_trimmed_to_nothing
            .fetch_add(1, Ordering::Relaxed);
        Some(Reason::TrimmedToNothing)
    } else {
        match crate::filter::check_with_gc(seq_len, qual, || raw_gc_fraction(&record), &cfg.filter)
        {
            Some(reason) => {
                tracing::trace!(
                    segment = 1,
                    of = 1,
                    start = 0,
                    end = seq_len,
                    len = seq_len,
                    reason = reason.label(),
                    "Segment dropped"
                );
                tracing::trace!(produced = 1, "Every segment filtered");
                counters.record_segment_drop(reason);
                counters.reads_all_filtered.fetch_add(1, Ordering::Relaxed);
                Some(Reason::Dropped(reason))
            },
            None => {
                tracing::trace!(
                    segment = 1,
                    of = 1,
                    start = 0,
                    end = seq_len,
                    len = seq_len,
                    "Segment kept"
                );
                None
            },
        }
    };
    if let Some(reason) = rejected {
        if counters.wants_rejects() {
            let mut edit = full_window_mod_edit(&record, mod_block, &cfg.remove_tags)?
                .unwrap_or_else(|| RecordEdit::unchanged(seq_len, &cfg.remove_tags));
            edit.reason = Some(reason);
            let rejected = build_record(&record, edit).with_context(read_context(&record))?;
            counters.reject(RejectItem::Bam(rejected))?;
        }
        return Ok(None);
    }

    // The window spans the whole record, so a survivor's output is its input.
    counters.output_reads.fetch_add(1, Ordering::Relaxed);
    counters
        .output_bases
        .fetch_add(seq_len as u64, Ordering::Relaxed);
    counters.reads_with_output.fetch_add(1, Ordering::Relaxed);

    let output = match full_window_mod_edit(&record, mod_block, &cfg.remove_tags)? {
        Some(edit) => {
            BamOutputRecord::Built(build_record(&record, edit).with_context(read_context(&record))?)
        },
        None => BamOutputRecord::Raw(record),
    };
    Ok(Some(output))
}

pub(super) fn run_raw_bam_full_window_seq(
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
            Some(BamOutputRecord::Built(bytes)) => sink.write_record_bytes(header, &bytes)?,
            None => {},
        }
    }
    Ok(counters.snapshot())
}

pub(super) fn run_raw_bam_full_window_parallel(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<bam::Record>> + Send,
    sink: &mut crate::io::bam::BamSink,
    level: u8,
    cfg: &Config,
    counters: &Arc<Counters>,
) -> anyhow::Result<Stats> {
    run_parallel(
        records,
        BAM_BATCH,
        |record: &bam::Record| record.sequence().len(),
        cfg,
        &mut KeyedSinks::single(sink),
        |record, cfg, out: &mut Parts<BamOutputRecord>| {
            ensure_raw_trimmable(&record)?;
            out.part(0)
                .extend(process_raw_full_window(record, cfg, counters)?);
            Ok(())
        },
        |records| pack_bam_blocks(header, level, records),
        |sink, _key, blocks: &Vec<u8>| sink.write_blocks(blocks),
        counters,
    )
}
