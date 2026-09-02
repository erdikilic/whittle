use std::fs::File;
use std::io::{self, Write};
use std::num::NonZero;
use std::path::Path;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _; // write_header / write_alignment_record
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::{self as sam};

/// A boxed, owning iterator over raw BAM records (or per-record errors). A
/// `bam::Record` owns the validated BAM record bytes but decodes fields lazily;
/// the workflow converts it to `RecordBuf` on a render worker instead of doing
/// all structured decoding on the serial reader thread.
pub type RawRecordIter = Box<dyn Iterator<Item = anyhow::Result<bam::Record>> + Send>;

/// Worker count for a noodles BGZF reader or writer.
///
/// noodles 0.51 builds a private Rayon pool per reader/writer from this count.
/// Up to 0.48 the count was ignored and jobs went to Rayon's global registry,
/// so a caller configured that registry instead; passing the count explicitly is
/// what replaced it. Getting this wrong is silent: the type checks either way,
/// the codec just runs single-threaded.
fn workers_nonzero(workers: usize) -> NonZero<usize> {
    NonZero::new(workers.max(1)).unwrap_or(NonZero::<usize>::MIN)
}

/// Error (naming the read) if the record is aligned. uBAM only in v1.
pub fn ensure_unaligned(rec: &RecordBuf) -> anyhow::Result<()> {
    let flags = rec.flags();
    let name = || {
        rec.name()
            .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
            .unwrap_or_else(|| "<unnamed>".to_string())
    };
    if !flags.is_unmapped() {
        anyhow::bail!(
            "read {} is aligned (mapped); whittle v1 supports unaligned BAM (uBAM) only",
            name()
        );
    }
    // A record flagged reverse stores SEQ as the reverse complement of the read,
    // and htslib decodes its MM right to left with complemented bases
    // (`sam_mods.c`, the `BAM_FREVERSE` branches). whittle trims and renumbers
    // left to right, so it would crop the wrong ends and relocate every call.
    // Basecallers do not emit `0x4|0x10`, but the SAM spec does not forbid it and
    // `samtools view -f 4` of an aligned file preserves it, so refuse rather than
    // silently disagree with every other tool.
    if flags.is_reverse_complemented() {
        anyhow::bail!(
            "read {} is flagged reverse-complemented; whittle trims in read \
             orientation and cannot keep position-indexed tags correct for it",
            name()
        );
    }
    Ok(())
}

/// The pre-spec lowercase spellings of the base-modification tags.
///
/// htslib still reads them (`sam_mods.c` falls back to `Mm` when `MM` is absent,
/// and to `Ml` for `Ml`), so old guppy and megalodon output decodes fine
/// elsewhere while whittle, which looks only for the uppercase tags, would copy
/// them through untouched onto a trimmed sequence and silently relocate every
/// call.
pub(crate) const LEGACY_MOD_TAGS: [[u8; 2]; 2] = [*b"Mm", *b"Ml"];

/// Error if the record carries legacy `Mm`/`Ml` rather than `MM`/`ML`.
///
/// Refused rather than rewritten: supporting both spellings means choosing which
/// to emit, and quietly modernizing a tag is its own surprise. Loudly declining
/// beats the silent corruption that copying them through produces.
pub fn ensure_modern_mod_tags(rec: &RecordBuf) -> anyhow::Result<()> {
    for t in LEGACY_MOD_TAGS {
        if rec.data().get(&Tag::new(t[0], t[1])).is_some() {
            anyhow::bail!(
                "read {} carries the legacy `{}` base-modification tag; whittle rewrites only \
                 the current `MM`/`ML` spelling, so trimming this record would leave its \
                 modification calls pointing at the wrong bases",
                rec.name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_else(|| "<unnamed>".to_string()),
                String::from_utf8_lossy(&t)
            );
        }
    }
    Ok(())
}

/// Open a BAM reader; MT-bgzf when `workers > 1`. Returns the header and a Send
/// owning raw-record iterator.
pub fn reader(
    input: Option<&Path>,
    workers: usize,
) -> anyhow::Result<(sam::Header, RawRecordIter)> {
    let inner: Box<dyn io::Read + Send> = match input {
        Some(p) => Box::new(File::open(p)?),
        None => Box::new(io::stdin()),
    };
    reader_from(inner, workers)
}

/// Like `reader`, but over an already-open stream rather than a path/stdin. Used
/// by the single-file dispatch so a stdin BAM whose first bytes were consumed for
/// format sniffing (and chained back into `inner`) is read from the true start.
/// Re-opening `io::stdin()` would drop those already-consumed bytes. MT-bgzf when
/// `workers > 1`.
pub fn reader_from(
    inner: Box<dyn io::Read + Send>,
    workers: usize,
) -> anyhow::Result<(sam::Header, RawRecordIter)> {
    if workers > 1 {
        let mt = bgzf::io::MultithreadedReader::with_worker_count(workers_nonzero(workers), inner);
        let mut r = bam::io::Reader::from(mt);
        let header = r.read_header()?;
        Ok((header, Box::new(RawRecordIterImpl { reader: r })))
    } else {
        let mut r = bam::io::Reader::new(inner);
        let header = r.read_header()?;
        Ok((header, Box::new(RawRecordIterImpl { reader: r })))
    }
}

struct RawRecordIterImpl<R: io::Read> {
    reader: bam::io::Reader<R>,
}

impl<R: io::Read> Iterator for RawRecordIterImpl<R> {
    type Item = anyhow::Result<bam::Record>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut record = bam::Record::default();
        match self.reader.read_record(&mut record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(record)),
            Err(e) => Some(Err(e.into())),
        }
    }
}

/// A BAM output sink: single-threaded bgzf (t1) or multithreaded bgzf (t>1).
pub enum BamSink {
    Single(bam::io::Writer<bgzf::io::Writer<Box<dyn Write + Send>>>),
    Multi(bam::io::Writer<bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>>),
}

/// Build the sink (header written), MT-bgzf when `workers > 1`. `level` is the
/// bgzf DEFLATE compression level (0-9 per the CLI, though libdeflate accepts up
/// to 12); it is applied to both the single- and multi-threaded encoders.
pub fn writer(
    output: Option<&Path>,
    header: &sam::Header,
    workers: usize,
    level: u8,
) -> anyhow::Result<BamSink> {
    let clevel = bgzf::io::writer::CompressionLevel::new(level)
        .ok_or_else(|| anyhow::anyhow!("invalid bgzf compression level {level} (expected 0-12)"))?;
    let inner: Box<dyn Write + Send> = match output {
        Some(p) => Box::new(File::create(p)?),
        None => Box::new(io::stdout()),
    };
    if workers > 1 {
        let mt = bgzf::io::multithreaded_writer::Builder::default()
            .set_compression_level(clevel)
            .set_worker_count(workers_nonzero(workers))
            .build_from_writer(inner);
        let mut w = bam::io::Writer::from(mt);
        w.write_header(header)?;
        Ok(BamSink::Multi(w))
    } else {
        // Build the single-threaded bgzf writer explicitly (rather than
        // `bam::io::Writer::new`, which would force the default level) so `level`
        // takes effect.
        let bgzf_w = bgzf::io::writer::Builder::default()
            .set_compression_level(clevel)
            .build_from_writer(inner);
        let mut w = bam::io::Writer::from(bgzf_w);
        w.write_header(header)?;
        Ok(BamSink::Single(w))
    }
}

impl BamSink {
    pub fn write_record(&mut self, header: &sam::Header, rec: &RecordBuf) -> io::Result<()> {
        match self {
            BamSink::Single(w) => w.write_alignment_record(header, rec),
            BamSink::Multi(w) => w.write_alignment_record(header, rec),
        }
    }

    pub fn write_raw_record(&mut self, header: &sam::Header, rec: &bam::Record) -> io::Result<()> {
        match self {
            BamSink::Single(w) => w.write_record(header, rec),
            BamSink::Multi(w) => w.write_record(header, rec),
        }
    }

    /// Flush + finalize (bgzf EOF block). Single: `try_finish`; Multi:
    /// `into_inner().finish()` (its `Drop` swallows errors, so this must be explicit).
    pub fn finish(self) -> anyhow::Result<()> {
        match self {
            BamSink::Single(mut w) => {
                w.try_finish()?;
                Ok(())
            },
            BamSink::Multi(w) => {
                w.into_inner().finish()?;
                Ok(())
            },
        }
    }
}

/// The output header: the input header with an `@PG` provenance record
/// (`ID:whittle`, program name and version) appended, and with `@HD SO:` set to
/// `unsorted` (and `GO`/`SS` removed) when `order_kept` is false, since a
/// multithreaded run without `--ordered` writes records in completion order.
///
/// The `@PG` record is best-effort: `Programs::add` fails on a duplicate ID and
/// cannot walk a dangling `PP` chain (`samtools reset` leaves
/// `@PG ID:samtools PP:basecaller` without an `ID:basecaller` record), in which
/// case the programs are left unchanged. The `@PG` line never blocks record
/// output.
pub(crate) fn provenance_header(mut header: sam::Header, order_kept: bool) -> sam::Header {
    use sam::header::record::value::Map;
    use sam::header::record::value::map::Program;
    use sam::header::record::value::map::header::tag as header_tag;
    use sam::header::record::value::map::program::tag;

    if let (false, Some(hd)) = (order_kept, header.header_mut()) {
        let fields = hd.other_fields_mut();
        fields.insert(header_tag::SORT_ORDER, "unsorted".into());
        fields.shift_remove(&header_tag::GROUP_ORDER);
        fields.shift_remove(&header_tag::SUBSORT_ORDER);
    }

    if has_dangling_program_chain(&header) {
        return header;
    }

    let program = Map::<Program>::builder()
        .insert(tag::NAME, "whittle")
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .build();

    if let Ok(program) = program {
        let _ = header.programs_mut().add("whittle", program);
    }

    header
}

/// True if the header's `@PG` chain is one `Programs::add` cannot walk safely.
///
/// `Programs::add` calls `Programs::leaves`, which indexes the program map
/// directly and panics when a `PP` names an absent ID, and which only terminates
/// a cycle that returns to the node it started from. A rho-shaped chain
/// (`pgA -> pgB -> pgC -> pgB`) has every ID present and never revisits `pgA`, so
/// it loops forever. Both shapes are rejected here by walking each chain with a
/// visited set.
fn has_dangling_program_chain(header: &sam::Header) -> bool {
    use std::collections::HashSet;

    use sam::header::record::value::map::program::tag;

    let programs = header.programs().as_ref();
    programs.keys().any(|start| {
        let mut seen: HashSet<&[u8]> = HashSet::new();
        let mut id: &[u8] = start.as_ref();
        loop {
            if !seen.insert(id) {
                return true; // revisited a node: cyclic
            }
            let Some(program) = programs.get(id) else {
                return true; // PP names an ID that is not a program: dangling
            };
            match program.other_fields().get(&tag::PREVIOUS_PROGRAM_ID) {
                Some(previous) => id = previous.as_ref(),
                None => return false, // reached the root of this chain
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::Program;
    use noodles_sam::header::record::value::map::program::tag;

    /// A dangling `@PG PP:` reference must leave the header unchanged because
    /// Noodles requires every parent program ID to exist.
    #[test]
    fn provenance_header_does_not_panic_on_dangling_pp_chain() {
        // `pg1` references a parent that is absent from the header.
        let dangling_program = Map::<Program>::builder()
            .insert(tag::PREVIOUS_PROGRAM_ID, "ghost")
            .build()
            .expect("valid PP field");

        let header = sam::Header::builder()
            .add_program("pg1", dangling_program)
            .build();

        assert!(has_dangling_program_chain(&header));

        let out_header = provenance_header(header, true);

        assert!(
            !out_header.programs().as_ref().contains_key(&b"whittle"[..]),
            "expected no whittle @PG line to be added when the existing chain is dangling"
        );
    }

    /// A rho-shaped chain (`pgA -> pgB -> pgC -> pgB`) has no absent ID, so the
    /// old dangling-only check passed it through to `Programs::add`, whose
    /// `leaves()` walk only terminates on a cycle that returns to its start node.
    /// Walking from `pgA` never revisits `pgA`, so it looped forever at 100% CPU.
    #[test]
    fn provenance_header_rejects_a_cycle_that_excludes_the_entry_node() {
        fn with_pp(previous: &str) -> Map<Program> {
            Map::<Program>::builder()
                .insert(tag::PREVIOUS_PROGRAM_ID, previous)
                .build()
                .expect("valid PP field")
        }

        let header = sam::Header::builder()
            .add_program("pgA", with_pp("pgB"))
            .add_program("pgB", with_pp("pgC"))
            .add_program("pgC", with_pp("pgB"))
            .build();

        assert!(
            has_dangling_program_chain(&header),
            "a rho-shaped chain must be rejected before `Programs::add` sees it"
        );

        // Reaching this line at all is the assertion: the old code hung here.
        let out_header = provenance_header(header, true);
        assert!(
            !out_header.programs().as_ref().contains_key(&b"whittle"[..]),
            "no @PG line should be added when the existing chain cannot be walked"
        );
    }

    /// A self-referential record (`pgA -> pgA`) is the degenerate cycle.
    #[test]
    fn provenance_header_rejects_a_self_referential_program() {
        let header = sam::Header::builder()
            .add_program(
                "pgA",
                Map::<Program>::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pgA")
                    .build()
                    .expect("valid PP field"),
            )
            .build();
        assert!(has_dangling_program_chain(&header));
    }

    /// A valid program chain receives the `whittle` provenance record.
    #[test]
    fn provenance_header_adds_whittle_program_on_clean_header() {
        let header = sam::Header::default();
        assert!(!has_dangling_program_chain(&header));

        let out_header = provenance_header(header, true);

        assert!(
            out_header
                .programs()
                .roots()
                .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"whittle"),
            "expected an @PG record with ID whittle in the output header, got {:?}",
            out_header.programs()
        );
    }
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;

    use super::*;

    #[test]
    fn unmapped_ok_mapped_rejected() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        assert!(ensure_unaligned(&rec).is_ok());

        *rec.flags_mut() = Flags::empty(); // mapped
        let err = ensure_unaligned(&rec).unwrap_err().to_string();
        assert!(err.contains("r1"));
        assert!(err.contains("aligned"));
    }

    #[test]
    fn mt_writer_roundtrips_through_mt_reader() {
        use noodles_sam::alignment::record::Flags;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mt.bam");

        // Write two unmapped records through a 4-worker MT BamSink.
        let header = sam::Header::default();
        let mut sink = writer(Some(&path), &header, 4, 6).unwrap();
        for name in [b"r1".as_slice(), b"r2".as_slice()] {
            let mut rec = RecordBuf::default();
            *rec.flags_mut() = Flags::UNMAPPED;
            *rec.name_mut() = Some(name.into());
            *rec.sequence_mut() = b"ACGT".to_vec().into();
            *rec.quality_scores_mut() = vec![40u8; 4].into();
            sink.write_record(&header, &rec).unwrap();
        }
        sink.finish().unwrap();

        // Read back through a 4-worker MT reader.
        let (_h, records) = reader(Some(&path), 4).unwrap();
        let names: Vec<Vec<u8>> = records
            .map(|r| r.unwrap().name().map(|n| n.to_vec()).unwrap_or_default())
            .collect();
        assert_eq!(names.len(), 2);
        assert!(names.contains(&b"r1".to_vec()) && names.contains(&b"r2".to_vec()));
    }

    #[test]
    fn provenance_header_marks_unordered_output_unsorted() {
        use sam::header::record::value::Map;
        use sam::header::record::value::map::header::tag;

        let hd = Map::<sam::header::record::value::map::Header>::builder()
            .insert(tag::SORT_ORDER, "queryname")
            .insert(tag::GROUP_ORDER, "query")
            .build()
            .unwrap();
        let header = sam::Header::builder().set_header(hd).build();

        let kept = provenance_header(header.clone(), true);
        let fields = kept.header().unwrap().other_fields();
        assert_eq!(fields.get(&tag::SORT_ORDER).map(|v| v.as_slice()), Some(&b"queryname"[..]));
        assert!(fields.contains_key(&tag::GROUP_ORDER));

        let unordered = provenance_header(header, false);
        let fields = unordered.header().unwrap().other_fields();
        assert_eq!(fields.get(&tag::SORT_ORDER).map(|v| v.as_slice()), Some(&b"unsorted"[..]));
        assert!(!fields.contains_key(&tag::GROUP_ORDER));
    }
}
