// End-to-end BAM→FASTQ/.gz conversion over the compiled binary. Builds a small
// uBAM fixture (a plain read, and a read with RG + MM/ML/MN mods), converts it,
// and checks header tags. The load-bearing correctness check is `cross_check`:
// the FASTQ-header MM/ML/MN must equal what the BAM→BAM path writes (which is
// itself htslib-oracle-verified by tests/bam_mods_oracle.rs).
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

    // read2: RG + mods. C at seq idx 0,1,3,4,5,7; MM occ 0,2,3 -> abs 0,3,4; ML [10,20,30].
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

fn run(args: &[&str], input: &Path, output: &Path) {
    Command::cargo_bin("whittle")
        .unwrap()
        .args(args)
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output)
        .assert()
        .success();
}

fn read2_header_line(fastq: &str) -> &str {
    fastq
        .lines()
        .find(|l| l.starts_with("@read2"))
        .expect("no read2 header in output")
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
    // read2: RG verbatim + reconstructed mod block; window [2,8) -> "C+m,0,0;" ML 20,30 MN 6.
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
    assert!(
        !s.contains("MM:Z"),
        "mods must be dropped under none: {s:?}"
    );
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

    // decode the gz and compare to the plain conversion.
    let mut gz = flate2::read::MultiGzDecoder::new(std::fs::File::open(&out).unwrap());
    let mut decoded = String::new();
    gz.read_to_string(&mut decoded).unwrap();
    assert_eq!(
        read2_header_line(&decoded),
        "@read2\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"
    );
}

#[test]
fn cross_check_fastq_header_mods_equal_bam_path() {
    // The FASTQ-header MM/ML/MN must be byte-identical to the BAM→BAM output's,
    // transitively inheriting the htslib oracle guarantee from bam_mods_oracle.rs.
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
                other => panic!("no MM in bam: {other:?}"),
            };
            let ml = match buf.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
                Some(Value::Array(Array::UInt8(v))) => v.clone(),
                other => panic!("no ML in bam: {other:?}"),
            };
            let mn = match buf.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
                Some(Value::Int32(n)) => *n,
                other => panic!("no MN in bam: {other:?}"),
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
    let mm_bam = mm_bam.expect("read2 missing from bam output");

    let s = std::fs::read_to_string(&fq).unwrap();
    let header_line = read2_header_line(&s);
    assert!(
        header_line.ends_with(&mm_bam),
        "fastq header mods {header_line:?} must end with bam-path mods {mm_bam:?}"
    );
}

#[test]
fn folder_dispatch_bam_to_fastq() {
    // The single-file BAM->FASTQ path is covered above; this exercises the
    // FOLDER dispatch arm (`run_folder`'s Bam-family -> FASTQ branch in
    // src/lib.rs) by pointing `-i` at a directory instead of a file.
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
        .args(["--fastq-tags", "none", "-i"])
        .arg(&inp)
        .arg("-o")
        .arg(&out)
        .assert()
        .success()
        .stderr(predicates::str::contains(
            "--fastq-tags applies only to BAM->FASTQ",
        ));
}
