//! FASTQ output from BAM input: selected aux tags in the header and the record
//! body of each output window.

use super::*;

/// Appends the TAB-prefixed aux-tag block for one window to `tags`: carried
/// non-mod tags in source order with the `window_tag_updates` rewrites applied
/// in place and the removed ones skipped, per-base arrays sliced, then the
/// rebuilt MM/ML/MN block, then the added tags. Nothing is appended when
/// nothing is carried (the record then has a plain header). A `Malformed`
/// block is omitted. A tag named by `remove` is left out of the header, after
/// the rewrite, exactly as on BAM output.
#[allow(clippy::too_many_arguments)]
pub(super) fn push_fastq_tags(
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
pub(super) fn render_fastq_window(
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
pub(super) fn render_bam_fastq_read(
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
