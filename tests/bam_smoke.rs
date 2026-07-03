// End-to-end regression test for the uBAM pipeline: drives the compiled `chopping`
// binary over a real on-disk BAM file (reader -> provenance_header -> run_bam ->
// writer), rather than exercising `reconstruct_record` directly against synthetic
// `RecordBuf`s the way `pipeline::bam_tests` does. Catches wiring bugs that unit
// tests can't (header handling, writer generics, CLI flag plumbing).
use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::RecordBuf};
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

fn write_fixture(path: &std::path::Path) {
    let header = sam::Header::default();
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).unwrap();

    // Record 1: plain unmapped read, no mods.
    let mut r1 = RecordBuf::default();
    *r1.flags_mut() = Flags::UNMAPPED;
    *r1.name_mut() = Some(b"read1".into());
    *r1.sequence_mut() = b"ACGTACGTAC".to_vec().into();
    *r1.quality_scores_mut() = vec![40; 10].into();
    writer.write_alignment_record(&header, &r1).unwrap();

    // Record 2: unmapped read with MM/ML/MN mods on Cs.
    let mut r2 = RecordBuf::default();
    *r2.flags_mut() = Flags::UNMAPPED;
    *r2.name_mut() = Some(b"read2".into());
    *r2.sequence_mut() = b"CCACCCAC".to_vec().into();
    *r2.quality_scores_mut() = vec![35; 8].into();
    let data = r2.data_mut();
    data.insert(Tag::BASE_MODIFICATIONS, Value::String(b"C+m,0,1,0;".to_vec().into()));
    data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(vec![10, 20, 30])));
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
    writer.write_alignment_record(&header, &r2).unwrap();

    writer.try_finish().unwrap();
}

#[test]
fn bam_to_bam_end_to_end() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("in.bam");
    let out_path = dir.path().join("out.bam");
    write_fixture(&in_path);

    Command::cargo_bin("chopping")
        .unwrap()
        .args(["--in-format", "bam", "--out-format", "bam", "--head-crop", "2", "-i"])
        .arg(&in_path)
        .arg("-o")
        .arg(&out_path)
        .assert()
        .success();

    // Read back the output BAM: check @PG provenance + reconstructed records.
    let mut reader = bam::io::Reader::new(std::fs::File::open(&out_path).unwrap());
    let header = reader.read_header().unwrap();

    assert!(
        header.programs().roots().any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"chopping"),
        "expected an @PG record with ID chopping in the output header, got {:?}",
        header.programs()
    );

    let mut out_records = Vec::new();
    let mut buf = RecordBuf::default();
    loop {
        match reader.read_record_buf(&header, &mut buf).unwrap() {
            0 => break,
            _ => out_records.push(buf.clone()),
        }
    }
    assert_eq!(out_records.len(), 2, "both reads should survive --head-crop 2 (len 10, 8 > min-length)");

    // read1: "ACGTACGTAC" head-cropped by 2 -> "GTACGTAC" (8 bases), no mods.
    let o1 = &out_records[0];
    assert_eq!(o1.name().unwrap().to_vec(), b"read1");
    assert_eq!(AsRef::<[u8]>::as_ref(o1.sequence()), b"GTACGTAC");
    assert!(o1.data().get(&Tag::BASE_MODIFICATIONS).is_none());

    // read2: "CCACCCAC" (C at seq idx 0,1,3,4,5,7) head-cropped by 2 -> "ACCCAC" (6 bases).
    // MM "C+m,0,1,0" walks the C occurrences with skip-counts 0,1,0 -> occurrences
    // 0,2,3 -> abs positions 0,3,4 (ML [10,20,30] one-per-position). Window [2,8)
    // keeps abs 3,4 (drops abs 0); renumbered against the window's C positions
    // (3,4,5,7) that's occurrences 0,1 -> deltas [0,0], ML [20,30].
    let o2 = &out_records[1];
    assert_eq!(o2.name().unwrap().to_vec(), b"read2");
    assert_eq!(AsRef::<[u8]>::as_ref(o2.sequence()), b"ACCCAC");
    let mm = match o2.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => s.to_vec(),
        other => panic!("expected MM tag, got {other:?}"),
    };
    assert_eq!(mm, b"C+m,0,0;");
    let ml = match o2.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        other => panic!("expected ML tag, got {other:?}"),
    };
    assert_eq!(ml, vec![20, 30]);
    let mn = match o2.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
        Some(Value::Int32(n)) => *n,
        other => panic!("expected MN tag, got {other:?}"),
    };
    assert_eq!(mn, 6, "MN should equal the output segment length");
}
