//! BAM reading and writing over noodles: raw-record readers, single- and multithreaded BGZF sinks, and record-level input guards.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
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

/// Returns the worker count for a noodles BGZF reader or writer.
///
/// noodles builds a private Rayon pool per reader or writer from this count;
/// the global Rayon registry is not consulted. A count of zero is clamped to
/// one. An incorrect count is not a type error: the codec runs single-threaded.
fn workers_nonzero(workers: usize) -> NonZero<usize> {
    NonZero::new(workers.max(1)).unwrap_or(NonZero::<usize>::MIN)
}

/// Errors (naming the read) if the record is aligned or flagged
/// reverse-complemented; only unaligned BAM (uBAM) input is supported.
pub fn ensure_unaligned(rec: &RecordBuf) -> anyhow::Result<()> {
    let flags = rec.flags();
    let name = || {
        rec.name()
            .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
            .unwrap_or_else(|| "<unnamed>".to_string())
    };
    if !flags.is_unmapped() {
        anyhow::bail!(
            "read {} is aligned (mapped); only unaligned BAM (uBAM) input is supported",
            name()
        );
    }
    // A record flagged reverse stores SEQ as the reverse complement of the read,
    // and htslib decodes its MM right to left with complemented bases
    // (`sam_mods.c`, the `BAM_FREVERSE` branches). whittle trims and renumbers
    // left to right, so it would crop the wrong ends and relocate every call.
    // Basecallers do not emit `0x4|0x10`, but the SAM spec does not forbid it and
    // `samtools view -f 4` of an aligned file preserves it, so such a record is
    // refused rather than trimmed in the opposite orientation from every
    // htslib-based consumer.
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
/// and to `Ml` when `ML` is absent), so guppy and megalodon output decodes
/// correctly in htslib-based tools, while whittle, which reads only the
/// uppercase tags, would copy them through unchanged onto a trimmed sequence
/// and relocate every call.
pub(crate) const LEGACY_MOD_TAGS: [[u8; 2]; 2] = [*b"Mm", *b"Ml"];

/// Errors if the record carries legacy `Mm`/`Ml` rather than `MM`/`ML`.
///
/// Refused rather than rewritten: supporting both spellings would require
/// choosing which to emit, and rewriting the tag changes the record's schema.
/// Refusing avoids the corruption that copying the tags through would produce.
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

/// Opens a BAM reader (multithreaded BGZF when `workers > 1`) and returns the
/// header and a `Send` owning raw-record iterator.
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

/// Opens a BAM reader like `reader`, but over an already-open stream rather
/// than a path or stdin. Used by the single-file dispatch so a stdin BAM whose
/// first bytes were consumed for format sniffing (and chained back into
/// `inner`) is read from the true start; reopening `io::stdin()` would drop
/// those bytes. Multithreaded BGZF when `workers > 1`.
pub fn reader_from(
    inner: Box<dyn io::Read + Send>,
    workers: usize,
) -> anyhow::Result<(sam::Header, RawRecordIter)> {
    // The block reader issues a `read` per frame part (header, payload,
    // trailer); the buffer in front of the source coalesces them.
    let inner = BufReader::with_capacity(INPUT_BUFFER_CAPACITY, inner);
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

/// Capacity of the buffer between a BGZF block reader and its file or stdin.
const INPUT_BUFFER_CAPACITY: usize = 1 << 20;

/// Capacity of the buffer between a BGZF encoder and its file or stdout.
/// The encoder emits each block frame as several small writes (header, payload,
/// trailer); the buffer coalesces them into one `write` syscall per megabyte.
const OUTPUT_BUFFER_CAPACITY: usize = 1 << 20;

/// The buffered destination of a BAM sink.
type BufferedOutput = BufWriter<Box<dyn Write + Send>>;

/// A BAM output sink: single-threaded BGZF (`-t 1`) or multithreaded BGZF.
pub enum BamSink {
    /// Single-threaded BGZF writer.
    Single(bam::io::Writer<bgzf::io::Writer<BufferedOutput>>),
    /// Multithreaded BGZF writer.
    Multi(bam::io::Writer<bgzf::io::MultithreadedWriter<BufferedOutput>>),
}

/// Builds the sink with the header written; multithreaded BGZF when
/// `workers > 1`. `level` is the BGZF DEFLATE compression level (0-9 per the
/// CLI, though libdeflate accepts up to 12) and is applied to both encoders.
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
    let inner = BufWriter::with_capacity(OUTPUT_BUFFER_CAPACITY, inner);
    if workers > 1 {
        let mt = bgzf::io::multithreaded_writer::Builder::default()
            .set_compression_level(clevel)
            .set_worker_count(workers_nonzero(workers))
            .build_from_writer(inner);
        let mut w = bam::io::Writer::from(mt);
        w.write_header(header)?;
        Ok(BamSink::Multi(w))
    } else {
        // The single-threaded BGZF writer is built explicitly rather than through
        // `bam::io::Writer::new`, which would force the default level.
        let bgzf_w = bgzf::io::writer::Builder::default()
            .set_compression_level(clevel)
            .build_from_writer(inner);
        let mut w = bam::io::Writer::from(bgzf_w);
        w.write_header(header)?;
        Ok(BamSink::Single(w))
    }
}

impl BamSink {
    /// Writes one decoded record under `header`.
    pub fn write_record(&mut self, header: &sam::Header, rec: &RecordBuf) -> io::Result<()> {
        match self {
            BamSink::Single(w) => w.write_alignment_record(header, rec),
            BamSink::Multi(w) => w.write_alignment_record(header, rec),
        }
    }

    /// Writes one raw record under `header` without decoding it.
    pub fn write_raw_record(&mut self, header: &sam::Header, rec: &bam::Record) -> io::Result<()> {
        match self {
            BamSink::Single(w) => w.write_record(header, rec),
            BamSink::Multi(w) => w.write_record(header, rec),
        }
    }

    /// Finalizes the stream: the BGZF encoder flushes its last block and
    /// writes the EOF block, then the output buffer is flushed to the file.
    /// Both encoders hand back the buffer from `finish`, so the flush error
    /// surfaces here; the encoders' `Drop` impls swallow errors, so the call
    /// must be explicit.
    pub fn finish(self) -> anyhow::Result<()> {
        let mut inner = match self {
            BamSink::Single(w) => w.into_inner().finish()?,
            BamSink::Multi(w) => w.into_inner().finish()?,
        };
        inner.flush()?;
        Ok(())
    }
}

/// Returns the output header: the input header with an `@PG` provenance record
/// (`ID:whittle`, program name and version) appended, and with `@HD SO:` set to
/// `unsorted` (and `GO`/`SS` removed) when `order_kept` is false, since a
/// multithreaded run without `--ordered` writes records in completion order.
///
/// The `@PG` record is best-effort: `Programs::add` fails on a duplicate ID and
/// cannot walk a dangling `PP` chain (`samtools reset` leaves
/// `@PG ID:samtools PP:basecaller` without an `ID:basecaller` record), in which
/// case the programs are left unchanged. The `@PG` line never blocks record
/// output.
pub(crate) fn provenance_header(
    mut header: sam::Header,
    order_kept: bool,
    command_line: &str,
) -> sam::Header {
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
        .insert(tag::COMMAND_LINE, command_line)
        .build();

    // `Programs::add` links `PP` to each chain leaf and suffixes the ID when
    // `whittle` is already present.
    if let Ok(program) = program
        && let Err(e) = header.programs_mut().add("whittle", program)
    {
        tracing::warn!(error = %e, "The @PG provenance record was not added");
    }

    header
}

/// Returns true if the header's `@PG` chain is one `Programs::add` cannot walk
/// safely.
///
/// `Programs::add` calls `Programs::leaves`, which indexes the program map
/// directly and panics when a `PP` names an absent ID, and which only terminates
/// a cycle that returns to the node it started from. A rho-shaped chain
/// (`pgA -> pgB -> pgC -> pgB`) has every ID present and never revisits `pgA`, so
/// the walk does not terminate. Both shapes are rejected here by walking each
/// chain with a visited set.
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

    /// A dangling `@PG PP:` reference leaves the header unchanged because
    /// noodles requires every parent program ID to exist.
    #[test]
    fn provenance_header_does_not_panic_on_dangling_pp_chain() {
        // `pg1` references a parent that is absent from the header.
        let dangling_program = Map::<Program>::builder()
            .insert(tag::PREVIOUS_PROGRAM_ID, "ghost")
            .build()
            .expect("Valid PP field");

        let header = sam::Header::builder()
            .add_program("pg1", dangling_program)
            .build();

        assert!(has_dangling_program_chain(&header));

        let out_header = provenance_header(header, true, "whittle");

        assert!(
            !out_header.programs().as_ref().contains_key(&b"whittle"[..]),
            "Expected no whittle @PG line when the existing chain is dangling"
        );
    }

    /// A rho-shaped chain (`pgA -> pgB -> pgC -> pgB`) has no absent ID, so a
    /// dangling-only check would pass it to `Programs::add`, whose `leaves()`
    /// walk terminates only on a cycle that returns to its start node. A walk
    /// from `pgA` never revisits `pgA` and does not terminate.
    #[test]
    fn provenance_header_rejects_a_cycle_that_excludes_the_entry_node() {
        fn with_pp(previous: &str) -> Map<Program> {
            Map::<Program>::builder()
                .insert(tag::PREVIOUS_PROGRAM_ID, previous)
                .build()
                .expect("Valid PP field")
        }

        let header = sam::Header::builder()
            .add_program("pgA", with_pp("pgB"))
            .add_program("pgB", with_pp("pgC"))
            .add_program("pgC", with_pp("pgB"))
            .build();

        assert!(
            has_dangling_program_chain(&header),
            "A rho-shaped chain must be rejected before `Programs::add` sees it"
        );

        // Returning from `provenance_header` is the assertion: an unwalkable
        // chain must not loop.
        let out_header = provenance_header(header, true, "whittle");
        assert!(
            !out_header.programs().as_ref().contains_key(&b"whittle"[..]),
            "No @PG line should be added when the existing chain cannot be walked"
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
                    .expect("Valid PP field"),
            )
            .build();
        assert!(has_dangling_program_chain(&header));
    }

    /// A valid program chain receives the `whittle` provenance record.
    #[test]
    fn provenance_header_adds_whittle_program_on_clean_header() {
        let header = sam::Header::default();
        assert!(!has_dangling_program_chain(&header));

        let out_header = provenance_header(header, true, "whittle");

        assert!(
            out_header
                .programs()
                .roots()
                .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"whittle"),
            "Expected an @PG record with ID whittle in the output header, got {:?}",
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

        // Two unmapped records go through a 4-worker multithreaded `BamSink`.
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

        // The records are read back through a 4-worker multithreaded reader.
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

        let kept = provenance_header(header.clone(), true, "whittle");
        let fields = kept.header().unwrap().other_fields();
        assert_eq!(
            fields.get(&tag::SORT_ORDER).map(|v| v.as_slice()),
            Some(&b"queryname"[..])
        );
        assert!(fields.contains_key(&tag::GROUP_ORDER));

        let unordered = provenance_header(header, false, "whittle");
        let fields = unordered.header().unwrap().other_fields();
        assert_eq!(
            fields.get(&tag::SORT_ORDER).map(|v| v.as_slice()),
            Some(&b"unsorted"[..])
        );
        assert!(!fields.contains_key(&tag::GROUP_ORDER));
    }

    /// The provenance record carries the command line and links to the chain
    /// leaf; a second run gets a distinct ID.
    #[test]
    fn provenance_header_records_the_command_line_and_links_the_chain() {
        use sam::header::record::value::Map;
        use sam::header::record::value::map::Program;
        use sam::header::record::value::map::program::tag;

        let header = sam::Header::builder()
            .add_program("dorado", Map::<Program>::default())
            .build();
        let once = provenance_header(header, true, "whittle -i a.bam -o b.bam");
        let pg = once.programs().as_ref().get(&b"whittle"[..]).unwrap();
        assert_eq!(
            pg.other_fields()
                .get(&tag::COMMAND_LINE)
                .map(|v| v.as_slice()),
            Some(&b"whittle -i a.bam -o b.bam"[..])
        );
        assert_eq!(
            pg.other_fields()
                .get(&tag::PREVIOUS_PROGRAM_ID)
                .map(|v| v.as_slice()),
            Some(&b"dorado"[..])
        );

        let twice = provenance_header(once, true, "whittle -i b.bam -o c.bam");
        let ids: Vec<&[u8]> = twice
            .programs()
            .as_ref()
            .keys()
            .map(|k| k.as_ref())
            .collect();
        assert_eq!(
            ids.len(),
            3,
            "Two whittle records coexist with the dorado one: {ids:?}"
        );
        let second = twice
            .programs()
            .as_ref()
            .get(&b"whittle-whittle"[..])
            .unwrap();
        assert_eq!(
            second
                .other_fields()
                .get(&tag::PREVIOUS_PROGRAM_ID)
                .map(|v| v.as_slice()),
            Some(&b"whittle"[..])
        );
    }
}
