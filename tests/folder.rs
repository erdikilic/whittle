use std::fs::File;
use std::path::Path;

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use predicates::prelude::*;

fn chopping() -> Command {
    Command::cargo_bin("chopping").unwrap()
}

#[test]
fn folder_merge_fastq_sorted_and_ignores_non_read_files() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("a.fastq"), "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    std::fs::write(dir.path().join("b.fastq"), "@r2\nTTTTGGGG\n+\nIIIIIIII\n").unwrap();
    std::fs::write(dir.path().join("sequencing_summary.txt"), "junk\n").unwrap(); // ignored
    let out = dir.path().join("merged.fastq");

    chopping()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-H", "2", "-T", "2", "-t", "1"]) // -t 1 => deterministic order
        .assert()
        .success();

    // Sorted: a.fastq then b.fastq. Head/tail-crop 2 => GTAC, TTGG.
    let got = std::fs::read_to_string(&out).unwrap();
    assert_eq!(got, "@r1\nGTAC\n+\nIIII\n@r2\nTTGG\n+\nIIII\n");
}

fn write_ubam(path: &Path, name: &[u8], seq: &[u8], quals: Vec<u8>) {
    let header = noodles_sam::Header::default();
    let mut w = bam::io::Writer::new(File::create(path).unwrap());
    w.write_header(&header).unwrap();
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(name.into());
    *rec.sequence_mut() = seq.to_vec().into();
    *rec.quality_scores_mut() = quals.into();
    w.write_alignment_record(&header, &rec).unwrap();
    w.try_finish().unwrap();
}

#[test]
fn folder_merge_bam_two_files() {
    let dir = tempfile::tempdir().unwrap();
    write_ubam(&dir.path().join("a.bam"), b"r1", b"ACGTACGT", vec![40; 8]);
    write_ubam(&dir.path().join("b.bam"), b"r2", b"TTTTGGGG", vec![40; 8]);
    let out = dir.path().join("merged.bam");

    chopping()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-H", "2", "-T", "2", "-t", "1"])
        .assert()
        .success();

    // Read the merged BAM back: 2 records, @PG chopping present.
    let mut r = bam::io::Reader::new(File::open(&out).unwrap());
    let hdr = r.read_header().unwrap();
    assert!(
        hdr.programs()
            .roots()
            .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"chopping"),
        "expected @PG chopping in merged header"
    );
    let mut count = 0usize;
    let mut buf = RecordBuf::default();
    while r.read_record_buf(&hdr, &mut buf).unwrap() != 0 {
        count += 1;
    }
    assert_eq!(count, 2);
}

#[test]
fn empty_folder_errors() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("notes.txt"), "x").unwrap();
    chopping()
        .arg("-i")
        .arg(dir.path())
        .assert()
        .failure()
        .stderr(predicates::str::contains("no FASTQ or BAM"));
}

#[test]
fn folder_rerun_does_not_reingest_its_own_output() {
    // When `-o` lands inside the `-i <dir>`, a rerun must exclude the previously
    // written output from the input set rather than merge (and truncate) it into
    // itself.
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("a.fastq"), "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    let out = dir.path().join("merged.fastq");

    chopping()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-t", "1"])
        .assert()
        .success();
    let first = std::fs::read_to_string(&out).unwrap();

    chopping()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-t", "1"])
        .assert()
        .success()
        .stderr(predicate::str::contains("excluding the output file"));
    let second = std::fs::read_to_string(&out).unwrap();

    assert_eq!(
        first, second,
        "rerun must produce the same output, not re-ingest it"
    );
}

#[test]
fn folder_bam_to_fastq_rerun_excludes_cross_format_output() {
    // A BAM folder producing a FASTQ output *inside itself*: on rerun the stale
    // merged.fastq must be excluded before the mixed-family check — otherwise the
    // folder (a.bam + merged.fastq) would wrongly look "mixed" and error.
    let dir = tempfile::tempdir().unwrap();
    write_ubam(&dir.path().join("a.bam"), b"r1", b"ACGTACGT", vec![40; 8]);
    let out = dir.path().join("merged.fastq");

    // First run creates merged.fastq inside the BAM folder.
    chopping()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["--out-format", "fastq", "-t", "1"])
        .assert()
        .success();
    let first = std::fs::read_to_string(&out).unwrap();

    // Rerun must succeed (not error as "mixed"), exclude the output, same result.
    chopping()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["--out-format", "fastq", "-t", "1"])
        .assert()
        .success()
        .stderr(predicate::str::contains("excluding the output file"));
    assert_eq!(std::fs::read_to_string(&out).unwrap(), first);
}

fn write_ubam_with_rg(path: &Path, name: &[u8], rg: &str) {
    use noodles_sam::header::record::value::Map;
    use noodles_sam::header::record::value::map::ReadGroup;
    let header = noodles_sam::Header::builder()
        .add_read_group(rg, Map::<ReadGroup>::default())
        .build();
    let mut w = bam::io::Writer::new(File::create(path).unwrap());
    w.write_header(&header).unwrap();
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(name.into());
    *rec.sequence_mut() = b"ACGTACGT".to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; 8].into();
    w.write_alignment_record(&header, &rec).unwrap();
    w.try_finish().unwrap();
}

#[test]
fn folder_merge_bam_warns_on_differing_read_groups() {
    // Folder merge keeps only the first header, so records from a file declaring a
    // different @RG would reference a read group missing from the merged output.
    let dir = tempfile::tempdir().unwrap();
    write_ubam_with_rg(&dir.path().join("a.bam"), b"r1", "rg_a");
    write_ubam_with_rg(&dir.path().join("b.bam"), b"r2", "rg_b");
    let out = dir.path().join("merged.bam");

    chopping()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-t", "1"])
        .assert()
        .success()
        .stderr(predicate::str::contains("different @RG"));
}
