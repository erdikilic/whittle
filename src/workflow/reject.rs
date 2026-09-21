//! Rejected records: every input read or trimmed segment that does not reach
//! the main output, written to `--rejected-output` in the output's format family
//! with a `wr:Z` tag naming the reason, by one writer thread fed from the
//! reader and the render workers.

use std::path::Path;
use std::sync::mpsc::{SyncSender, sync_channel};
use std::thread::JoinHandle;

use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::{self as sam};

use crate::filter::DropReason;
use crate::io::Format;

/// The aux tag carrying the rejection reason.
pub(crate) const REASON_TAG: Tag = Tag::new(b'w', b'r');

/// Items queued to the writer per worker.
const QUEUE_PER_WORKER: usize = 4;

/// Why a read or segment was rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Reason {
    /// The read did not satisfy `--tag-filter`.
    TagFilter,
    /// Trimming left no segment.
    TrimmedToNothing,
    /// A trimmed segment failed a post-trim filter.
    Dropped(DropReason),
}

impl Reason {
    /// The tag value naming the reason.
    pub(crate) fn label(self) -> &'static str {
        match self {
            Reason::TagFilter => "tag_filter",
            Reason::TrimmedToNothing => "trimmed_to_nothing",
            Reason::Dropped(DropReason::TooShort) => "too_short",
            Reason::Dropped(DropReason::TooLong) => "too_long",
            Reason::Dropped(DropReason::LowQuality) => "low_quality",
            Reason::Dropped(DropReason::HighQuality) => "high_quality",
            Reason::Dropped(DropReason::Gc) => "gc",
        }
    }
}

/// Sets the reason tag on a rejected BAM record.
pub(crate) fn tag_record(rec: &mut RecordBuf, reason: Reason) {
    rec.data_mut()
        .insert(REASON_TAG, Value::String(reason.label().into()));
}

/// Appends the reason tag as a tagged FASTQ header field.
pub(crate) fn push_fastq_tag(out: &mut Vec<u8>, reason: Reason) {
    out.extend_from_slice(b"\twr:Z:");
    out.extend_from_slice(reason.label().as_bytes());
}

/// One rejected record, rendered for the output format family.
pub(crate) enum RejectItem {
    /// A BAM record carrying the reason tag.
    Bam(RecordBuf),
    /// A complete FASTQ record, reason tag included.
    Fastq(Vec<u8>),
    /// The last item; the writer finishes its file.
    End,
}

/// The sending side, cloned into every producer.
#[derive(Clone)]
pub(crate) struct Rejects {
    tx: SyncSender<RejectItem>,
}

impl Rejects {
    /// Queues one rejected record; an error means the writer stopped.
    pub(crate) fn send(&self, item: RejectItem) -> anyhow::Result<()> {
        self.tx
            .send(item)
            .map_err(|_| anyhow::anyhow!("the rejected-record writer stopped"))
    }
}

/// The writer thread and the channel end that closes it.
pub(crate) struct RejectWriter {
    tx: SyncSender<RejectItem>,
    handle: JoinHandle<anyhow::Result<()>>,
}

impl RejectWriter {
    /// Opens `path` in `format` and starts the writer. `header` is the BAM
    /// header for BAM output and `None` for the FASTQ family.
    pub(crate) fn start(
        path: &Path,
        format: Format,
        level: u8,
        header: Option<sam::Header>,
        workers: usize,
    ) -> anyhow::Result<(Rejects, RejectWriter)> {
        let (tx, rx) = sync_channel::<RejectItem>((workers * QUEUE_PER_WORKER).max(4));
        let handle = match format {
            Format::Bam => {
                let header = header.expect("BAM rejected output is given its header");
                let mut sink = crate::io::bam::writer(Some(path), &header, false, level)?;
                std::thread::Builder::new()
                    .name("rejected-writer".into())
                    .spawn(move || {
                        for item in rx.iter() {
                            match item {
                                RejectItem::Bam(rec) => sink.write_record(&header, &rec)?,
                                RejectItem::Fastq(_) => {
                                    anyhow::bail!("a FASTQ record reached the BAM rejected output")
                                },
                                RejectItem::End => break,
                            }
                        }
                        sink.finish()
                    })?
            },
            Format::Fastq | Format::FastqGz | Format::FastqBgzf => {
                let mut sink = crate::io::fastq::writer_to(Some(path), format, level, false)?;
                std::thread::Builder::new()
                    .name("rejected-writer".into())
                    .spawn(move || {
                        use std::io::Write;
                        for item in rx.iter() {
                            match item {
                                RejectItem::Fastq(bytes) => sink.write_all(&bytes)?,
                                RejectItem::Bam(_) => {
                                    anyhow::bail!("a BAM record reached the FASTQ rejected output")
                                },
                                RejectItem::End => break,
                            }
                        }
                        sink.finish()
                    })?
            },
        };
        Ok((Rejects { tx: tx.clone() }, RejectWriter { tx, handle }))
    }

    /// Sends the end marker, waits for the writer and returns its result.
    pub(crate) fn finish(self) -> anyhow::Result<()> {
        // A send failure means the writer already stopped; its error follows.
        let _ = self.tx.send(RejectItem::End);
        self.handle
            .join()
            .map_err(|_| anyhow::anyhow!("the rejected-record writer panicked"))?
    }
}

/// Resolves the rejected output format from `path` and checks that it shares
/// the main output's format family.
pub(crate) fn resolve_format(path: &Path, out_fmt: Format) -> anyhow::Result<Format> {
    let format = crate::io::resolve_output(Some(path), out_fmt);
    if format.family() != out_fmt.family() {
        anyhow::bail!(
            "--rejected-output {} is {} but the output is {}; rejected records are written in the \
             output's format family",
            path.display(),
            format.family(),
            out_fmt.family()
        );
    }
    Ok(format)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reasons_have_stable_labels() {
        assert_eq!(Reason::TagFilter.label(), "tag_filter");
        assert_eq!(Reason::TrimmedToNothing.label(), "trimmed_to_nothing");
        assert_eq!(Reason::Dropped(DropReason::TooShort).label(), "too_short");
        assert_eq!(Reason::Dropped(DropReason::Gc).label(), "gc");
    }

    #[test]
    fn rejected_format_follows_the_output_family() {
        assert_eq!(
            resolve_format(Path::new("r.fastq.gz"), Format::Fastq).unwrap(),
            Format::FastqGz
        );
        assert_eq!(
            resolve_format(Path::new("r.bam"), Format::Bam).unwrap(),
            Format::Bam
        );
        assert!(resolve_format(Path::new("r.bam"), Format::Fastq).is_err());
        assert!(resolve_format(Path::new("r.fastq"), Format::Bam).is_err());
    }

    #[test]
    fn tags_name_the_reason_in_both_families() {
        let mut rec = RecordBuf::default();
        tag_record(&mut rec, Reason::Dropped(DropReason::LowQuality));
        assert_eq!(
            rec.data().get(&REASON_TAG),
            Some(&Value::String(b"low_quality".as_slice().into()))
        );
        let mut out = b"@r1".to_vec();
        push_fastq_tag(&mut out, Reason::TagFilter);
        assert_eq!(out, b"@r1\twr:Z:tag_filter");
    }
}
