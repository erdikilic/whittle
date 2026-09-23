//! Reconstruction of BAM output records: decoding, the kept windows of a read,
//! and the rebuilt record of each window.

use super::*;

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

/// Computes the edit that builds output record `window` of `src`: SEQ/QUAL
/// cut to the window, `MM`/`ML`/`MN` rebuilt, `sa` coverage runs re-encoded,
/// stale signal-space tags rewritten or dropped, and the name updated for
/// splits and PacBio interval crops. `build_record` applies it to the raw
/// record: aux tags are copied in source order with the rewritten ones
/// replaced in place, per-base arrays cut to the window, removed ones skipped
/// and added ones appended. A `Malformed` block is removed. An untrimmed,
/// unsplit record with an `Absent` or `Consistent` block and nothing to remove
/// yields `None`: its output is its input.
///
/// `remove` names the tags `--remove-tag` drops. Removal
/// is applied last, to the rewritten tag set, so a removed tag whittle
/// maintains (`MM`, the move table, a per-base array) is absent from the
/// output rather than left stale.
pub(super) fn window_edit<'a>(
    src: &RecordBuf,
    window: Window,
    mod_block: ModBlock,
    indexed: Option<&IndexedMods>,
    moves: Option<&MoveIndex<'_>>,
    remove: &'a TagRemoval,
) -> Option<RecordEdit<'a>> {
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
    let name = (trimmed || split).then(|| {
        let name: &[u8] = src.name().map(AsRef::as_ref).unwrap_or_default();
        let coords = query_span(src).map(|(qs0, _)| window_coords(qs0, start, end));
        segment_name(platform, name, window, coords)
    });

    // Tags with dedicated handling: `Some` replaces the source value in place,
    // or is appended when the source lacks the tag; `None` removes it.
    let mut updates: TagUpdates = Vec::new();
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
    // The `sa` coverage runs are re-encoded for the window; runs that do not
    // cover the read leave the tag unchanged.
    let sa = Tag::new(RLE_COVERAGE_TAG[0], RLE_COVERAGE_TAG[1]);
    if trimmed
        && let Some(value) = src.data().get(&sa)
        && let Some(sliced) = windowed_value(RLE_COVERAGE_TAG, value, orig_len, start, end)
    {
        updates.push((sa, Some(sliced)));
    }

    Some(RecordEdit {
        start,
        end,
        name,
        updates,
        remove,
        reason: None,
    })
}

/// The per-read state the decoded BAM workflows share.
pub(super) struct PreparedRead<'a> {
    /// The record's bases.
    pub(super) seq: &'a [u8],
    /// The record's per-base qualities.
    pub(super) qual: &'a [u8],
    /// The state of the record's modification block.
    pub(super) mod_block: ModBlock,
    /// The original-coordinate interval retained by barcode restriction, `None`
    /// without an adapter source or when the record carries no verified `bi`.
    pub(super) barcode: Option<(usize, usize)>,
}

/// Runs the per-read guards and bookkeeping shared by the decoded workflows:
/// refuses aligned reads and legacy mod tags, requires full per-base quality,
/// classifies the modification block, counting a malformed one, counts a
/// malformed per-base tag, and resolves the barcode window from the verified
/// `bi` spans, counting an unusable or unverified `bi`.
pub(super) fn prepare_read<'a>(
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
pub(super) fn render_windows(
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

/// Renders one record for BAM output: every surviving window is built from
/// the raw record `raw` with the edit its decoded form `rec` determines and
/// handed to `emit`, and every rejected window is sent to the rejected output.
/// Shared by the sequential and parallel drivers. `emit` receives `None` for
/// a window whose output record is the input record.
pub(super) fn render_bam_read(
    header: &sam::Header,
    raw: &bam::Record,
    rec: &RecordBuf,
    cfg: &Config,
    counters: &Counters,
    mut emit: impl FnMut(Option<Vec<u8>>) -> io::Result<()>,
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
        let edit = window_edit(rec, window, mod_block, indexed, moves, &cfg.remove_tags);
        match (edit, reason) {
            (None, None) => Ok(emit(None)?),
            (Some(edit), None) => Ok(emit(Some(build_record(raw, edit)?))?),
            (edit, Some(reason)) => {
                let mut edit =
                    edit.unwrap_or_else(|| RecordEdit::unchanged(seq.len(), &cfg.remove_tags));
                edit.reason = Some(reason);
                counters.reject(RejectItem::Bam(build_record(raw, edit)?))
            },
        }
    })
}
