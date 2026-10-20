//! FASTQ output from BAM input: header tags, modification tags and per-base
//! arrays.

use super::*;

/// The FASTQ header spells an empty group as `MM:Z:C+m;` with a zero-length
/// `ML:B:C` array.
#[test]
fn bam2fq_keeps_empty_mod_group_after_a_crop() {
    let rec = ubam_with_mods(b"CCAC", vec![40; 4], b"C+m,0;", vec![10]);
    let cfg = cfg_bam2fq(None, 2, FastqTags::All);
    let mut out = Vec::new();
    run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&rec))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(
        s.starts_with("@r1\tMM:Z:C+m;\tML:B:C\tMN:i:2\n"),
        "Got: {s:?}"
    );
}

#[test]
fn bam2fq_all_carries_rg_and_reconstructed_mods() {
    let cfg = cfg_bam2fq(None, 2, FastqTags::All);
    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
    let s = String::from_utf8(out).unwrap();
    // The header carries RG unchanged and the reconstructed mod block; the
    // sequence is head-cropped by 2.
    assert!(
        s.starts_with("@r1\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6\n"),
        "Got: {s:?}"
    );
    assert!(s.contains("\nACCCAC\n+\n"), "Cropped sequence wrong: {s:?}");
}

#[test]
fn bam2fq_only_mm_ml_drops_rg() {
    let cfg = cfg_bam2fq(None, 2, FastqTags::parse("MM,ML").unwrap());
    let mut out = Vec::new();
    run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(!s.contains("RG:Z"), "RG must be dropped: {s:?}");
    assert!(
        s.contains("MM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"),
        "Mods missing: {s:?}"
    );
}

#[test]
fn bam2fq_none_is_plain_fastq() {
    let cfg = cfg_bam2fq(None, 2, FastqTags::None);
    let mut out = Vec::new();
    run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(out, b"@r1\nACCCAC\n+\nDDDDDD\n"); // 35+33 = 'D'
}

/// A split at the low-quality base gives each segment its own reconstructed
/// mods.
#[test]
fn bam2fq_split_suffixes_and_segments_mods() {
    let cfg = cfg_bam2fq(Some(QualityOp::runs(20, 1)), 0, FastqTags::All);
    // Seq CCAC, `C+m` at occurrences 0 and 2 (positions 0 and 3); quality
    // good, good, bad, good, so the split is [0,2), [3,4).
    let mut rec = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0,1;", vec![100, 200]);
    rec.data_mut()
        .insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));
    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&rec))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(stats.output_reads, 2);
    let s = String::from_utf8(out).unwrap();
    // Segment 1 = [0,2) "CC" keeps the position-0 mod; segment 2 = [3,4) "C"
    // keeps the position-3 mod.
    assert!(
        s.contains("@r1_segment_1\tMM:Z:C+m,0;\tML:B:C,100\tMN:i:2"),
        "Segment 1: {s:?}"
    );
    assert!(
        s.contains("@r1_segment_2\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:1"),
        "Segment 2: {s:?}"
    );
}

#[test]
fn bam2fq_no_mods_read_is_plain() {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"plain".into());
    *rec.sequence_mut() = b"ACGT".to_vec().into();
    *rec.quality_scores_mut() = vec![40; 4].into();
    let cfg = cfg_bam2fq(None, 0, FastqTags::All);
    let mut out = Vec::new();
    run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&rec))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    assert_eq!(out, b"@plain\nACGT\n+\nIIII\n");
}

/// Companion for the BAM-to-FASTQ path: an MM-only source must emit a FASTQ
/// header with `MM` + `MN` but no `ML:B:C` field.
#[test]
fn bam2fq_mm_without_ml_omits_ml_field() {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = b"CCAC".to_vec().into();
    *rec.quality_scores_mut() = vec![40; 4].into();
    rec.data_mut().insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1;".to_vec().into()),
    );

    let cfg = cfg_bam2fq(None, 0, FastqTags::All);
    let mut out = Vec::new();
    run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&rec))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(
        s.contains("MM:Z:C+m,0,1;\tMN:i:4"),
        "Expected MM+MN, got: {s:?}"
    );
    assert!(
        !s.contains("ML:B"),
        "MM-only record must not emit an ML field: {s:?}"
    );
}

/// The FASTQ header carries no MM/ML/MN for a malformed block, and the run
/// counts the read.
#[test]
fn malformed_mod_block_is_omitted_and_counted_on_fastq_output() {
    for variant in MALFORMED_MOD_VARIANTS {
        for head in [0, 2] {
            let cfg = cfg_bam2fq(None, head, FastqTags::All);
            let mut out = Vec::new();
            let stats = run_bam_to_fastq(
                &sam::Header::default(),
                [Ok(raw_record(&malformed_mod_record(variant)))].into_iter(),
                &mut KeyedSinks::single(&mut out),
                None,
                &cfg,
                &Arc::new(Counters::default()),
            )
            .unwrap();
            assert_eq!(stats.malformed_mod_reads, 1, "{variant} head={head}");
            let s = String::from_utf8(out).unwrap();
            assert!(
                s.starts_with("@r1\tRG:Z:grp\n"),
                "{variant} head={head}: mod block must be omitted, got {s:?}"
            );
        }
    }
}

#[test]
fn bam2fq_slices_kinetics_in_header() {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = b"ACGTAC".to_vec().into();
    *rec.quality_scores_mut() = vec![40; 6].into();
    rec.data_mut().insert(
        Tag::new(b'i', b'p'),
        Value::Array(Array::UInt8(vec![10, 11, 12, 13, 14, 15])),
    );

    let cfg = cfg_bam2fq(None, 2, FastqTags::All); // head crop 2, window [2,6)
    let mut out = Vec::new();
    run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&rec))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(
        s.contains("ip:B:C,12,13,14,15"),
        "Kinetics not sliced in the FASTQ header: {s:?}"
    );
}

#[test]
fn bam2fq_slices_reverse_strand_kinetics_from_the_other_end() {
    let cfg = cfg_bam2fq(None, 2, FastqTags::All); // head crop 2
    let mut out = Vec::new();
    run_bam_to_fastq(
        &sam::Header::default(),
        [Ok(raw_record(&pacbio_kinetics_record()))].into_iter(),
        &mut KeyedSinks::single(&mut out),
        None,
        &cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(s.contains("\tfi:B:C,2,3,4,5,6,7,8,9\t"), "{s:?}");
    assert!(s.contains("\tri:B:C,10,11,12,13,14,15,16,17\n"), "{s:?}");
}

/// On BAM-to-FASTQ the removal applies to the tags carried into the header,
/// per field of the modification block as it does on BAM output.
#[test]
fn removal_drops_the_tag_from_a_fastq_header() {
    let rec = record_with_mixed_tags();
    let mut cfg = cfg_bam2fq(None, 2, FastqTags::All);
    cfg.remove_tags = removal(&["MM", "RG"]);
    let (_stats, text) = bam2fq(vec![rec.clone()], &cfg);

    let header = text.lines().next().unwrap();
    assert!(!header.contains("MM:Z:"), "MM was removed: {header}");
    assert!(!header.contains("RG:Z:"), "RG was removed: {header}");
    assert!(header.contains("ML:B:C,20,30"), "{header}");
    assert!(header.contains("MN:i:6"), "{header}");
    assert!(header.contains("ip:B:C,2,3,4,5,6,7"), "{header}");

    // Without the flag the same run carries both.
    let mut keep = cfg.clone();
    keep.remove_tags = crate::config::TagRemoval::default();
    let (_stats, text) = bam2fq(vec![rec], &keep);
    let header = text.lines().next().unwrap();
    assert!(header.contains("MM:Z:"), "{header}");
    assert!(header.contains("RG:Z:run1"), "{header}");
}

/// Without an `MM:Z` block, a lone `MN` and the tags of a non-string `MM`
/// are copied into the FASTQ header as BAM output copies them, and a
/// selection that leaves them out drops them.
#[test]
fn bam2fq_copies_mod_tags_without_a_string_mm() {
    let record = |name: &[u8], tags: Vec<(Tag, Value)>| {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(name.into());
        *rec.sequence_mut() = b"ACGT".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 4].into();
        for (tag, value) in tags {
            rec.data_mut().insert(tag, value);
        }
        rec
    };
    let lone = record(
        b"r1",
        vec![(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4))],
    );
    let odd = record(
        b"r2",
        vec![
            (Tag::BASE_MODIFICATIONS, Value::Int32(5)),
            (
                Tag::BASE_MODIFICATION_PROBABILITIES,
                Value::Array(Array::UInt8(vec![1, 2, 3])),
            ),
            (Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4)),
        ],
    );
    let cfg = cfg_bam2fq(None, 1, FastqTags::All);
    let (_stats, text) = bam2fq(vec![lone.clone(), odd.clone()], &cfg);
    let headers: Vec<&str> = text.lines().step_by(4).collect();
    assert_eq!(
        headers,
        ["@r1\tMN:i:4", "@r2\tMM:i:5\tML:B:C,1,2,3\tMN:i:4"],
        "{text}"
    );

    let cfg = cfg_bam2fq(None, 1, FastqTags::parse("ML").unwrap());
    let (_stats, text) = bam2fq(vec![lone, odd], &cfg);
    let headers: Vec<&str> = text.lines().step_by(4).collect();
    assert_eq!(headers, ["@r1", "@r2\tML:B:C,1,2,3"], "{text}");
}

/// A dorado read of 30 bases with every ONT tag the rules touch: a move table
/// of stride 5 (a stay after every third base), the mods block, the sample
/// counts, poly-A, timing, MinKNOW and barcode tags, and a read-length custom
/// array. A low-quality base at 14 lets a quality split cut it in two.
fn dorado_read(name: &[u8], group: &[u8]) -> RecordBuf {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(name.into());
    *rec.sequence_mut() = b"ACGTCCGTACCATGCAGTCCATGACGTCAC".to_vec().into();
    let mut qual = vec![30; 30];
    qual[14] = 5;
    *rec.quality_scores_mut() = qual.into();
    let mut mv = vec![5];
    for i in 0..30 {
        mv.push(1);
        if i % 3 == 0 {
            mv.push(0);
        }
    }
    let d = rec.data_mut();
    d.insert(Tag::READ_GROUP, Value::String(group.into()));
    d.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m?,0,2,1;".as_slice().into()),
    );
    d.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![200, 100, 50])),
    );
    d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(30));
    d.insert(Tag::new(b'm', b'v'), Value::Array(Array::Int8(mv)));
    // 40 blocks of 5 samples after a front trim of 100.
    d.insert(Tag::new(b't', b's'), Value::Int32(100));
    d.insert(Tag::new(b'n', b's'), Value::Int32(300));
    d.insert(
        Tag::new(b'p', b'a'),
        Value::Array(Array::Int32(vec![280, 260, 290, -1, -1])),
    );
    d.insert(Tag::new(b'p', b't'), Value::Int32(5));
    d.insert(
        Tag::new(b's', b't'),
        Value::String(b"2024-06-21T10:00:00.000+00:00".as_slice().into()),
    );
    d.insert(Tag::new(b'd', b'u'), Value::Float(0.06));
    d.insert(Tag::new(b'q', b's'), Value::Float(25.0));
    d.insert(Tag::new(b'r', b'n'), Value::Int32(42));
    d.insert(Tag::new(b'c', b'h'), Value::Int32(7));
    d.insert(Tag::new(b'm', b'e'), Value::Int32(12));
    d.insert(
        Tag::new(b'e', b'r'),
        Value::String(b"signal_positive".as_slice().into()),
    );
    d.insert(
        Tag::new(b'B', b'C'),
        Value::String(b"barcode01".as_slice().into()),
    );
    d.insert(
        Tag::new(b'x', b'c'),
        Value::Array(Array::UInt8((0..30).collect())),
    );
    rec
}

/// The reads of the BAM and FASTQ tag parity check: DNA and RNA dorado reads
/// whose read-group ids name their models, a subread of an earlier dorado
/// split, a PacBio HiFi read with kinetics and fixed-size arrays, and a mods
/// block without `MN`.
fn parity_reads() -> Vec<RecordBuf> {
    let dna = dorado_read(b"dna1", b"4f2a_dna_r10.4.1_e8.2_400bps_sup@v5.0.0");
    let rna = dorado_read(b"rna1", b"4f2a_rna004_130bps_sup@v5.1.0");
    let mut subread = dorado_read(b"dna2", b"4f2a_dna_r10.4.1_e8.2_400bps_sup@v5.0.0");
    let d = subread.data_mut();
    d.insert(
        Tag::new(b'p', b'i'),
        Value::String(b"parent-read".as_slice().into()),
    );
    d.insert(Tag::new(b's', b'p'), Value::Int32(1000));
    d.remove(&Tag::new(b'p', b'a'));
    let mut hifi = pacbio_hifi_record(HIFI);
    let mut qual = vec![40; 10];
    qual[5] = 5;
    *hifi.quality_scores_mut() = qual.into();
    let d = hifi.data_mut();
    for (t, offset) in [(*b"fi", 0u8), (*b"ri", 10), (*b"ip", 20), (*b"pw", 30)] {
        d.insert(
            Tag::new(t[0], t[1]),
            Value::Array(Array::UInt8((offset..offset + 10).collect())),
        );
    }
    d.insert(
        Tag::new(b's', b'n'),
        Value::Array(Array::Float(vec![5.5, 9.0, 4.25, 7.75])),
    );
    d.insert(Tag::new(b'r', b'q'), Value::Float(0.999));
    let mut no_mn = record_with_mixed_tags();
    *no_mn.name_mut() = Some(b"mixed".into());
    no_mn
        .data_mut()
        .remove(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH);
    vec![dna, rna, subread, hifi, no_mn]
}

/// Splits FASTQ text into reads whose names hold the full header lines.
fn fastq_reads(text: &str) -> Vec<crate::record::ReadRecord> {
    text.lines()
        .collect::<Vec<_>>()
        .chunks(4)
        .map(|r| crate::record::ReadRecord {
            name: r[0].as_bytes()[1..].to_vec(),
            seq: r[1].as_bytes().to_vec(),
            qual: r[3].bytes().map(|q| q - 33).collect(),
        })
        .collect()
}

/// Decodes FASTQ text into records whose aux data are the header tags.
fn fastq_records(text: &str) -> Vec<RecordBuf> {
    fastq_reads(text)
        .into_iter()
        .map(|read| crate::io::tagged::record_from_tagged(read).unwrap())
        .collect()
}

/// A record's name, bases, qualities and aux tags, with integers compared by
/// value, since BAM and SAM text store them at different widths.
#[derive(Debug, PartialEq)]
struct Comparable {
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    tags: Vec<(String, String)>,
}

impl Comparable {
    fn of(rec: &RecordBuf) -> Self {
        let mut tags: Vec<(String, String)> = rec
            .data()
            .iter()
            .map(|(t, v)| {
                let value = aux_integer(v).map_or_else(|| format!("{v:?}"), |n| format!("int {n}"));
                (String::from_utf8_lossy(t.as_ref()).into_owned(), value)
            })
            .collect();
        tags.sort();
        Self {
            name: name_of(rec),
            seq: rec.sequence().as_ref().to_vec(),
            qual: rec.quality_scores().as_ref().to_vec(),
            tags,
        }
    }
}

/// Runs tagged FASTQ text through the FASTQ workflow.
fn tagged2fq(text: &str, cfg: &Config) -> String {
    let mut out = Vec::new();
    crate::workflow::run_fastq(
        fastq_reads(text).into_iter().map(Ok),
        &mut KeyedSinks::single(&mut out),
        None,
        cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    String::from_utf8(out).unwrap()
}

/// Every tag rule applies alike to BAM output, to BAM-to-FASTQ output and to
/// tagged FASTQ input: for crops, a quality split and a removal, with and
/// without `--update-moves`, each output read carries the same name, bases,
/// qualities and aux tags on all three paths.
#[test]
fn fastq_output_tags_match_bam_output_tags() {
    let reads = parity_reads();
    let (_, tagged) = bam2fq(reads.clone(), &cfg_bam2fq(None, 0, FastqTags::All));
    let head = cfg_bam2fq(None, 3, FastqTags::All);
    let mut tail = cfg_bam2fq(None, 0, FastqTags::All);
    tail.trim.tail = 4;
    let split = cfg_bam2fq(Some(QualityOp::runs(20, 1)), 0, FastqTags::All);
    let mut removed = cfg_bam2fq(Some(QualityOp::runs(20, 1)), 2, FastqTags::All);
    removed.remove_tags = removal(&["kinetics", "ML"]);
    for (label, cfg) in [
        ("head crop", head),
        ("tail crop", tail),
        ("quality split", split),
        ("crop, split and removal", removed),
    ] {
        for update_moves in [false, true] {
            let cfg = Config {
                update_moves,
                ..cfg.clone()
            };
            let (_, bam) = bam2bam(reads.clone(), &cfg);
            let (_, fastq) = bam2fq(reads.clone(), &cfg);
            let fastq = fastq_records(&fastq);
            let from_tagged = fastq_records(&tagged2fq(&tagged, &cfg));
            let bam: Vec<_> = bam.iter().map(Comparable::of).collect();
            let fastq: Vec<_> = fastq.iter().map(Comparable::of).collect();
            let from_tagged: Vec<_> = from_tagged.iter().map(Comparable::of).collect();
            assert_eq!(fastq, bam, "{label}, update_moves {update_moves}");
            assert_eq!(from_tagged, bam, "{label}, update_moves {update_moves}");
            // The rewritten move table reaches the FASTQ header.
            let carries_mv = fastq
                .iter()
                .filter(|read| read.name.starts_with(b"dna1"))
                .all(|read| read.tags.iter().any(|(t, _)| t == "mv"));
            assert_eq!(carries_mv, update_moves, "{label}");
        }
    }
}
