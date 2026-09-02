//! End-to-end BAM-to-FASTQ and BAM-to-FASTQ.gz conversion over the compiled
//! binary. Builds a small uBAM fixture (a plain read, and a read with `RG` and
//! `MM`/`ML`/`MN` mods), converts it, and checks header tags. The load-bearing
//! check is `cross_check_fastq_header_mods_equal_bam_path`: the FASTQ-header
//! `MM`/`ML`/`MN` equals what the BAM-to-BAM path writes, which
//! `tests/bam_mods_oracle.rs` verifies against htslib.

use std::io::Read;
use std::path::Path;

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};

/// Writes the two-read uBAM fixture: `read1` plain, `read2` with `RG` and mods.
fn write_fixture(path: &Path) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    // read1: plain, no tags.
    let mut r1 = RecordBuf::default();
    *r1.flags_mut() = Flags::UNMAPPED;
    *r1.name_mut() = Some(b"read1".into());
    *r1.sequence_mut() = b"ACGTACGTAC".to_vec().into();
    *r1.quality_scores_mut() = vec![40; 10].into();
    w.write_alignment_record(&header, &r1).unwrap();

    // read2: RG and mods. C at seq idx 0,1,3,4,5,7; MM occ 0,2,3 -> abs 0,3,4; ML [10,20,30].
    let mut r2 = RecordBuf::default();
    *r2.flags_mut() = Flags::UNMAPPED;
    *r2.name_mut() = Some(b"read2".into());
    *r2.sequence_mut() = b"CCACCCAC".to_vec().into();
    *r2.quality_scores_mut() = vec![35; 8].into();
    let d = r2.data_mut();
    d.insert(Tag::from(*b"RG"), Value::String(b"grp1".as_slice().into()));
    d.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1,0;".to_vec().into()),
    );
    d.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![10, 20, 30])),
    );
    d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
    w.write_alignment_record(&header, &r2).unwrap();

    w.try_finish().unwrap();
}

/// Runs the binary with `args` plus `-i input -o output` and asserts success.
fn run(args: &[&str], input: &Path, output: &Path) {
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(args)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output)
        .assert()
        .success();
}

/// The `@read2` header line of a FASTQ text.
fn read2_header_line(fastq: &str) -> &str {
    fastq
        .lines()
        .find(|l| l.starts_with("@read2"))
        .expect("No read2 header in output")
}

#[test]
fn bam_to_fastq_all_carries_rg_and_mods() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.bam");
    let out = dir.path().join("out.fastq");
    write_fixture(&inp);

    run(&["--out-format", "fastq", "--head-crop", "2"], &inp, &out);

    let s = std::fs::read_to_string(&out).unwrap();
    // read1: plain, no tags.
    assert!(s.contains("@read1\nGTACGTAC\n+\n"), "read1 wrong: {s:?}");
    // read2: RG verbatim and the reconstructed mod block; window [2,8) -> "C+m,0,0;" ML 20,30 MN 6.
    assert_eq!(
        read2_header_line(&s),
        "@read2\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"
    );
}

#[test]
fn bam_to_fastq_none_is_plain() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.bam");
    let out = dir.path().join("out.fastq");
    write_fixture(&inp);

    run(
        &[
            "--out-format",
            "fastq",
            "--head-crop",
            "2",
            "--fastq-tags",
            "none",
        ],
        &inp,
        &out,
    );

    let s = std::fs::read_to_string(&out).unwrap();
    assert_eq!(read2_header_line(&s), "@read2"); // no tags
    assert!(!s.contains("MM:Z"), "Mods are dropped under none: {s:?}");
}

#[test]
fn bam_to_fastq_only_mm_ml_drops_rg() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.bam");
    let out = dir.path().join("out.fastq");
    write_fixture(&inp);

    run(
        &[
            "--out-format",
            "fastq",
            "--head-crop",
            "2",
            "--fastq-tags",
            "MM,ML",
        ],
        &inp,
        &out,
    );

    let s = std::fs::read_to_string(&out).unwrap();
    assert_eq!(
        read2_header_line(&s),
        "@read2\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"
    );
}

#[test]
fn bam_to_fastq_gz_roundtrips() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.bam");
    let out = dir.path().join("out.fastq.gz");
    write_fixture(&inp);

    run(
        &["--out-format", "fastq-gz", "--head-crop", "2", "-t", "4"],
        &inp,
        &out,
    );

    // Decode the gz output and compare it to the plain conversion.
    let mut gz = flate2::read::MultiGzDecoder::new(std::fs::File::open(&out).unwrap());
    let mut decoded = String::new();
    gz.read_to_string(&mut decoded).unwrap();
    assert_eq!(
        read2_header_line(&decoded),
        "@read2\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"
    );
}

/// The FASTQ-header `MM`/`ML`/`MN` is byte-identical to the BAM-to-BAM
/// output's, transitively inheriting the htslib oracle guarantee from
/// `tests/bam_mods_oracle.rs`.
#[test]
fn cross_check_fastq_header_mods_equal_bam_path() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.bam");
    let fq = dir.path().join("out.fastq");
    let ba = dir.path().join("out.bam");
    write_fixture(&inp);

    run(&["--out-format", "fastq", "--head-crop", "2"], &inp, &fq);
    run(&["--out-format", "bam", "--head-crop", "2"], &inp, &ba);

    // Extract MM/ML/MN from the BAM read2.
    let mut reader = bam::io::Reader::new(std::fs::File::open(&ba).unwrap());
    let header = reader.read_header().unwrap();
    let mut buf = RecordBuf::default();
    let mut mm_bam = None;
    while reader.read_record_buf(&header, &mut buf).unwrap() != 0 {
        if AsRef::<[u8]>::as_ref(buf.name().unwrap()) == b"read2" {
            let mm = match buf.data().get(&Tag::BASE_MODIFICATIONS) {
                Some(Value::String(s)) => s.to_vec(),
                other => panic!("No MM in the BAM output: {other:?}"),
            };
            let ml = match buf.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
                Some(Value::Array(Array::UInt8(v))) => v.clone(),
                other => panic!("No ML in the BAM output: {other:?}"),
            };
            let mn = match buf.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
                Some(Value::Int32(n)) => *n,
                other => panic!("No MN in the BAM output: {other:?}"),
            };
            // Render the same SAM-text block the FASTQ path would.
            let mut expect = format!("MM:Z:{}", String::from_utf8(mm).unwrap());
            expect.push_str("\tML:B:C");
            for b in &ml {
                expect.push_str(&format!(",{b}"));
            }
            expect.push_str(&format!("\tMN:i:{mn}"));
            mm_bam = Some(expect);
        }
    }
    let mm_bam = mm_bam.expect("Read2 missing from the BAM output");

    let s = std::fs::read_to_string(&fq).unwrap();
    let header_line = read2_header_line(&s);
    assert!(
        header_line.ends_with(&mm_bam),
        "FASTQ header mods {header_line:?} end with the BAM-path mods {mm_bam:?}"
    );
}

/// The single-file BAM-to-FASTQ path is covered above; this test exercises
/// `run_folder`'s BAM-family-to-FASTQ arm by pointing `-i` at a directory.
#[test]
fn folder_dispatch_bam_to_fastq() {
    let dir = tempfile::tempdir().unwrap();
    let sub = dir.path().join("barcode01");
    std::fs::create_dir_all(&sub).unwrap();
    let inp = sub.join("in.bam");
    let out = dir.path().join("out.fastq");
    write_fixture(&inp);

    run(&["--out-format", "fastq", "--head-crop", "2"], &sub, &out);

    let s = std::fs::read_to_string(&out).unwrap();
    assert_eq!(
        read2_header_line(&s),
        "@read2\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"
    );
}

#[test]
fn fastq_tags_on_fastq_input_prints_ignored_note() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.fastq");
    let out = dir.path().join("out.fastq");
    std::fs::write(&inp, b"@r\nACGT\n+\nIIII\n").unwrap();

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["--fastq-tags", "none", "-i"])
        .arg(&inp)
        .arg("-o")
        .arg(&out)
        .assert()
        .success()
        .stderr(predicates::str::contains(
            "--fastq-tags applies only to BAM-to-FASTQ",
        ));
}
