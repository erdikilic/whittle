//! BAM reading and writing over noodles: raw-record readers, single- and multithreaded BGZF sinks, and record-level input guards.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _; // write_header / write_alignment_record
use noodles_sam::alignment::record::Flags;
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

/// The read name for a message, `<unnamed>` for a record without one.
pub(crate) fn display_name(name: Option<&[u8]>) -> String {
    name.map_or_else(
        || "<unnamed>".to_string(),
        |n| String::from_utf8_lossy(n).into_owned(),
    )
}

/// The pre-spec lowercase spellings of the base-modification tags.
///
/// htslib reads them (`sam_mods.c` falls back to `Mm` when `MM` is absent,
/// and to `Ml` when `ML` is absent), so guppy and megalodon output decodes
/// correctly in htslib-based tools, while whittle, which reads only the
/// uppercase tags, would copy them through unchanged onto a trimmed sequence
/// and relocate every call.
pub(crate) const LEGACY_MOD_TAGS: [[u8; 2]; 2] = [*b"Mm", *b"Ml"];

/// Refuses a record the workflows cannot trim, naming the read: one that is
/// aligned, one flagged reverse-complemented, or one carrying a legacy
/// `Mm`/`Ml` tag (`legacy_tag`, the first such tag it carries). Shared by the
/// decoded guard (`ensure_trimmable`) and the raw-record guard in
/// `workflow::bam`, so both refuse with the same message. `name` renders the
/// read name for the message.
pub(crate) fn refuse_untrimmable(
    flags: Flags,
    legacy_tag: Option<[u8; 2]>,
    name: impl Fn() -> String,
) -> anyhow::Result<()> {
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
    // A legacy tag is refused rather than rewritten: supporting both spellings
    // would require choosing which to emit, and rewriting the tag changes the
    // record's schema. Refusing avoids the corruption that copying the tag
    // through would produce.
    if let Some(t) = legacy_tag {
        anyhow::bail!(
            "read {} carries the legacy `{}` base-modification tag; whittle rewrites only \
             the current `MM`/`ML` spelling, so trimming this record would leave its \
             modification calls pointing at the wrong bases",
            name(),
            String::from_utf8_lossy(&t)
        );
    }
    Ok(())
}

/// Errors (naming the read) if a decoded record is aligned, flagged
/// reverse-complemented or carries a legacy `Mm`/`Ml` tag; see
/// `refuse_untrimmable`. Only unaligned BAM (uBAM) input is supported.
pub fn ensure_trimmable(rec: &RecordBuf) -> anyhow::Result<()> {
    let legacy_tag = LEGACY_MOD_TAGS
        .into_iter()
        .find(|t| rec.data().get(&Tag::new(t[0], t[1])).is_some());
    refuse_untrimmable(rec.flags(), legacy_tag, || {
        display_name(rec.name().map(AsRef::as_ref))
    })
}

/// Opens a BAM reader over an already-open stream and returns the header and a
/// `Send` owning raw-record iterator. The single-file dispatch hands over a
/// stream whose sniffed first bytes are chained back in front, so a stdin BAM
/// is read from its true start. Multithreaded BGZF when `workers > 1`.
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

/// A BAM output sink. `Single` compresses in the writing thread for a
/// sequential run; `Blocks` takes BGZF blocks the render workers compressed
/// (see `workflow::run_parallel`) and writes the bytes through.
pub enum BamSink {
    /// Single-threaded BGZF writer.
    Single(bam::io::Writer<bgzf::io::Writer<BufferedOutput>>),
    /// Pre-compressed BGZF blocks, header already written.
    Blocks {
        /// The output, positioned after the header blocks.
        inner: BufferedOutput,
        /// The BGZF DEFLATE level the render workers compress at.
        level: u8,
    },
}

/// Builds the sink with the header written: `Blocks` when `parallel`,
/// `Single` otherwise. `level` is the BGZF DEFLATE compression level (0-9 per
/// the CLI, though libdeflate accepts up to 12).
pub fn writer(
    output: Option<&Path>,
    header: &sam::Header,
    parallel: bool,
    level: u8,
) -> anyhow::Result<BamSink> {
    let clevel = compression_level(level)?;
    let inner: Box<dyn Write + Send> = match output {
        Some(p) => Box::new(File::create(p)?),
        None => Box::new(io::stdout()),
    };
    let inner = BufWriter::with_capacity(OUTPUT_BUFFER_CAPACITY, inner);
    // The single-threaded BGZF writer is built explicitly rather than through
    // `bam::io::Writer::new`, which would force the default level.
    let bgzf_w = bgzf::io::writer::Builder::default()
        .set_compression_level(clevel)
        .build_from_writer(inner);
    let mut w = bam::io::Writer::from(bgzf_w);
    w.write_header(header)?;
    if !parallel {
        return Ok(BamSink::Single(w));
    }
    let mut bgzf_w = w.into_inner();
    bgzf_w.flush()?;
    Ok(BamSink::Blocks {
        inner: bgzf_w.into_inner(),
        level,
    })
}

/// Returns the BGZF compression level for `level`, or an error naming the
/// accepted range.
pub(crate) fn compression_level(level: u8) -> anyhow::Result<bgzf::io::writer::CompressionLevel> {
    bgzf::io::writer::CompressionLevel::new(level)
        .ok_or_else(|| anyhow::anyhow!("invalid bgzf compression level {level} (expected 0-12)"))
}

/// Encodes `records` under `header` into BGZF blocks at `level`: the bytes of
/// a BGZF stream fragment with complete blocks and no EOF block, which a
/// `BamSink::Blocks` writes through.
pub fn encode_blocks<'a>(
    header: &sam::Header,
    level: u8,
    records: impl IntoIterator<Item = &'a RecordBuf>,
) -> io::Result<Vec<u8>> {
    let clevel = compression_level(level).map_err(io::Error::other)?;
    let bgzf_w = bgzf::io::writer::Builder::default()
        .set_compression_level(clevel)
        .build_from_writer(Vec::new());
    let mut w = bam::io::Writer::from(bgzf_w);
    for rec in records {
        w.write_alignment_record(header, rec)?;
    }
    let mut bgzf_w = w.into_inner();
    bgzf_w.flush()?;
    Ok(bgzf_w.into_inner())
}

/// The BGZF EOF block: the empty block every BGZF stream ends with.
pub(crate) fn eof_block() -> Vec<u8> {
    bgzf::io::Writer::new(Vec::new())
        .finish()
        .expect("An in-memory BGZF writer finishes without I/O errors")
}

impl BamSink {
    /// Writes one decoded record under `header`; `Single` only.
    pub fn write_record(&mut self, header: &sam::Header, rec: &RecordBuf) -> io::Result<()> {
        match self {
            BamSink::Single(w) => w.write_alignment_record(header, rec),
            BamSink::Blocks { .. } => Err(io::Error::other(
                "a block sink takes compressed blocks, not records",
            )),
        }
    }

    /// Writes one raw record under `header` without decoding it; `Single` only.
    pub fn write_raw_record(&mut self, header: &sam::Header, rec: &bam::Record) -> io::Result<()> {
        match self {
            BamSink::Single(w) => w.write_record(header, rec),
            BamSink::Blocks { .. } => Err(io::Error::other(
                "a block sink takes compressed blocks, not records",
            )),
        }
    }

    /// Writes compressed BGZF blocks through; `Blocks` only.
    pub fn write_blocks(&mut self, blocks: &[u8]) -> io::Result<()> {
        match self {
            BamSink::Blocks { inner, .. } => inner.write_all(blocks),
            BamSink::Single(_) => Err(io::Error::other(
                "a record sink takes records, not compressed blocks",
            )),
        }
    }

    /// The BGZF level the render workers compress at, for a `Blocks` sink.
    pub fn block_level(&self) -> Option<u8> {
        match self {
            BamSink::Blocks { level, .. } => Some(*level),
            BamSink::Single(_) => None,
        }
    }

    /// Finalizes the stream: the last block and the EOF block are written,
    /// then the output buffer is flushed to the file. The encoder's `Drop`
    /// swallows errors, so the call must be explicit.
    pub fn finish(self) -> anyhow::Result<()> {
        let mut inner = match self {
            BamSink::Single(w) => w.into_inner().finish()?,
            BamSink::Blocks { mut inner, .. } => {
                inner.write_all(&eof_block())?;
                inner
            },
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
        assert!(ensure_trimmable(&rec).is_ok());

        *rec.flags_mut() = Flags::empty(); // mapped
        let err = ensure_trimmable(&rec).unwrap_err().to_string();
        assert!(err.contains("r1"));
        assert!(err.contains("aligned"));
    }

    /// The reverse-complement and legacy-tag refusals name the read and the
    /// offending fact.
    #[test]
    fn reverse_and_legacy_tag_records_are_refused() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED | Flags::REVERSE_COMPLEMENTED;
        *rec.name_mut() = Some(b"r1".into());
        let err = ensure_trimmable(&rec).unwrap_err().to_string();
        assert!(
            err.contains("r1") && err.contains("reverse-complemented"),
            "{err}"
        );

        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        rec.data_mut().insert(
            Tag::new(b'M', b'l'),
            noodles_sam::alignment::record_buf::data::field::Value::Array(
                noodles_sam::alignment::record_buf::data::field::value::Array::UInt8(vec![1]),
            ),
        );
        let err = ensure_trimmable(&rec).unwrap_err().to_string();
        assert!(
            err.contains("<unnamed>") && err.contains("legacy `Ml`"),
            "{err}"
        );
    }

    #[test]
    fn mt_writer_roundtrips_through_mt_reader() {
        use noodles_sam::alignment::record::Flags;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mt.bam");

        // Two unmapped records go through a block `BamSink`, compressed as one
        // batch the way a render worker does it.
        let header = sam::Header::default();
        let mut sink = writer(Some(&path), &header, true, 6).unwrap();
        let mut records = Vec::new();
        for name in [b"r1".as_slice(), b"r2".as_slice()] {
            let mut rec = RecordBuf::default();
            *rec.flags_mut() = Flags::UNMAPPED;
            *rec.name_mut() = Some(name.into());
            *rec.sequence_mut() = b"ACGT".to_vec().into();
            *rec.quality_scores_mut() = vec![40u8; 4].into();
            records.push(rec);
        }
        let blocks = encode_blocks(&header, 6, &records).unwrap();
        sink.write_blocks(&blocks).unwrap();
        sink.finish().unwrap();

        // The records are read back through a 4-worker multithreaded reader.
        let (_h, records) = reader_from(Box::new(File::open(&path).unwrap()), 4).unwrap();
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
