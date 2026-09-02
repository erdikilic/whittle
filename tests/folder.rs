//! Folder-merge mode (`-i <dir>`) over the compiled binary.

use std::fs::File;
use std::path::Path;

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use predicates::prelude::*;

/// The binary with `WHITTLE_LOG` cleared.
fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

#[test]
fn folder_merge_fastq_sorted_and_ignores_non_read_files() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("a.fastq"), "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    std::fs::write(dir.path().join("b.fastq"), "@r2\nTTTTGGGG\n+\nIIIIIIII\n").unwrap();
    std::fs::write(dir.path().join("sequencing_summary.txt"), "junk\n").unwrap(); // ignored
    let out = dir.path().join("merged.fastq");

    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-H", "2", "-T", "2", "-t", "1"]) // -t 1 gives deterministic order
        .assert()
        .success();

    // Sorted: a.fastq then b.fastq; head and tail crop 2 give GTAC, TTGG.
    let got = std::fs::read_to_string(&out).unwrap();
    assert_eq!(got, "@r1\nGTAC\n+\nIIII\n@r2\nTTGG\n+\nIIII\n");
}

#[test]
fn folder_merge_skips_hidden_files_and_orders_members_naturally() {
    let dir = tempfile::tempdir().unwrap();
    for (name, id) in [
        ("run_10.fastq", "r10"),
        ("run_2.fastq", "r2"),
        ("run_1.fastq", "r1"),
    ] {
        std::fs::write(dir.path().join(name), format!("@{id}\nACGT\n+\nIIII\n")).unwrap();
    }
    // An AppleDouble sidecar is not FASTQ and fails to parse if ingested.
    std::fs::write(dir.path().join("._run_1.fastq"), b"\x00\x05\x16\x07junk").unwrap();
    std::fs::write(dir.path().join(".hidden.fastq"), "@hidden\nACGT\n+\nIIII\n").unwrap();
    let out = dir.path().join("merged.fastq");

    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-t", "1"])
        .assert()
        .success();

    let got = std::fs::read_to_string(&out).unwrap();
    assert_eq!(
        got,
        "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n@r10\nACGT\n+\nIIII\n"
    );
}

/// Writes a one-read uBAM with the given name, sequence, and qualities.
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

    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-H", "2", "-T", "2", "-t", "1"])
        .assert()
        .success();

    // Read the merged BAM back: 2 records, @PG whittle present.
    let mut r = bam::io::Reader::new(File::open(&out).unwrap());
    let hdr = r.read_header().unwrap();
    assert!(
        hdr.programs()
            .roots()
            .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"whittle"),
        "Expected @PG whittle in the merged header"
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
    whittle()
        .arg("-i")
        .arg(dir.path())
        .assert()
        .failure()
        .stderr(predicates::str::contains("no FASTQ or BAM"));
}

/// A folder merge cannot overwrite one of its input files.
#[test]
fn folder_output_matching_a_real_input_is_rejected_and_preserves_it() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("a.fastq"), "@a\nACGT\n+\nIIII\n").unwrap();
    std::fs::write(dir.path().join("b.fastq"), "@b\nTTTT\n+\nIIII\n").unwrap();
    let a = dir.path().join("a.fastq");
    let before = std::fs::read(&a).unwrap();

    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&a)
        .args(["-t", "1"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("refusing to overwrite"));
    assert_eq!(
        std::fs::read(&a).unwrap(),
        before,
        "The real input file is untouched"
    );
}

/// A rerun whose `-o` is inside `-i <dir>` (the previous output is a read file,
/// indistinguishable from real input) hard-errors rather than overwriting.
#[test]
fn folder_rerun_with_output_inside_dir_hard_errors() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("a.fastq"), "@r1\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    let out = dir.path().join("merged.fastq");

    // First run: merged.fastq does not exist; the run succeeds and creates it.
    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-t", "1"])
        .assert()
        .success();
    let first = std::fs::read_to_string(&out).unwrap();

    // Rerun: merged.fastq exists in the directory; hard error, prior output kept.
    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-t", "1"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("refusing to overwrite"));
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        first,
        "Prior output preserved"
    );
}

/// A BAM folder producing a FASTQ output inside itself: the first run succeeds,
/// and the rerun (merged.fastq is then a read file in the folder) hard-errors
/// rather than overwriting.
#[test]
fn folder_bam_to_fastq_rerun_with_output_inside_dir_hard_errors() {
    let dir = tempfile::tempdir().unwrap();
    write_ubam(&dir.path().join("a.bam"), b"r1", b"ACGTACGT", vec![40; 8]);
    let out = dir.path().join("merged.fastq");

    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["--out-format", "fastq", "-t", "1"])
        .assert()
        .success();
    let first = std::fs::read_to_string(&out).unwrap();

    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["--out-format", "fastq", "-t", "1"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("refusing to overwrite"));
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        first,
        "Prior output preserved"
    );
}

/// Writes a one-read uBAM whose header declares the read group `rg`.
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

/// Folder merge preserves custom-FASTA semantics across the combined stream.
#[test]
fn folder_merge_custom_fasta_never_reduces_still_trims() {
    let present = "GGGGTTTTGGGGTTTTGGGG"; // 20 bp, present adapter
    let absent = "ACGACGACGACGACGACGAC"; // 20 bp, never in the reads
    let insert = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 40 bp
    let dir = tempfile::tempdir().unwrap();
    let mut fa = tempfile::NamedTempFile::new().unwrap();
    std::io::Write::write_fmt(
        &mut fa,
        format_args!(">present\n{present}\n>absent\n{absent}\n"),
    )
    .unwrap();

    let mut a = String::new();
    for i in 0..100 {
        a.push_str(&format!(
            "@r{i}\n{present}{insert}\n+\n{}\n",
            "I".repeat(60)
        ));
    }
    std::fs::write(dir.path().join("a.fastq"), a).unwrap();
    let mut b = String::new();
    for i in 100..200 {
        b.push_str(&format!(
            "@r{i}\n{present}{insert}\n+\n{}\n",
            "I".repeat(60)
        ));
    }
    std::fs::write(dir.path().join("b.fastq"), b).unwrap();

    let out = dir.path().join("merged.fastq");
    let res = whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args([
            "--adapter-fasta",
            fa.path().to_str().unwrap(),
            "-v",
            "-t",
            "1",
        ])
        .assert()
        .success();

    let stderr = String::from_utf8_lossy(&res.get_output().stderr).into_owned();
    assert!(
        !stderr.contains("Adapter presence"),
        "Custom --adapter-fasta disables detection outright in folder-merge mode as well: {stderr}"
    );

    let got = std::fs::read_to_string(&out).unwrap();
    assert!(got.contains(insert), "Insert kept");
    assert!(
        !got.contains(&format!("{present}{insert}")),
        "Present adapter trimmed off in the merged output"
    );
    assert_eq!(
        got.matches("@r").count(),
        200,
        "All 200 merged reads present"
    );
}

/// Folder merge keeps only the first header, so records from a file declaring a
/// different `@RG` would reference a read group missing from the merged output.
#[test]
fn folder_merge_bam_warns_on_differing_read_groups() {
    let dir = tempfile::tempdir().unwrap();
    write_ubam_with_rg(&dir.path().join("a.bam"), b"r1", "rg_a");
    write_ubam_with_rg(&dir.path().join("b.bam"), b"r2", "rg_b");
    let out = dir.path().join("merged.bam");

    whittle()
        .arg("-i")
        .arg(dir.path())
        .arg("-o")
        .arg(&out)
        .args(["-t", "1"])
        .assert()
        .success()
        .stderr(predicate::str::contains("different @RG"));
}
