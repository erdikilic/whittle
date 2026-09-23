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
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
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
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut out,
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
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut out,
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
        [Ok(raw_record(&read2_with_mods_and_rg()))].into_iter(),
        &mut out,
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
    let cfg = cfg_bam2fq(
        Some(QualityOp::Split {
            cutoff: 20,
            window: 1,
        }),
        0,
        FastqTags::All,
    );
    // Seq CCAC, `C+m` at occurrences 0 and 2 (positions 0 and 3); quality
    // good, good, bad, good, so the split is [0,2), [3,4).
    let mut rec = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0,1;", vec![100, 200]);
    rec.data_mut()
        .insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));
    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
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
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
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
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
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
                [Ok(raw_record(&malformed_mod_record(variant)))].into_iter(),
                &mut out,
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
        [Ok(raw_record(&rec))].into_iter(),
        &mut out,
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
        [Ok(raw_record(&pacbio_kinetics_record()))].into_iter(),
        &mut out,
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
