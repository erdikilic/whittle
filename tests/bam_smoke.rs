// End-to-end uBAM coverage through the compiled binary, including header,
// reader, workflow, writer, and CLI integration.
use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};

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
    data.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,0,1,0;".to_vec().into()),
    );
    data.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![10, 20, 30])),
    );
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

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "--head-crop",
            "2",
            "-i",
        ])
        .arg(&in_path)
        .arg("-o")
        .arg(&out_path)
        .assert()
        .success();

    // Read back the output BAM: check @PG provenance + reconstructed records.
    let mut reader = bam::io::Reader::new(std::fs::File::open(&out_path).unwrap());
    let header = reader.read_header().unwrap();

    assert!(
        header
            .programs()
            .roots()
            .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"whittle"),
        "expected an @PG record with ID whittle in the output header, got {:?}",
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
    assert_eq!(
        out_records.len(),
        2,
        "both reads should survive --head-crop 2 (len 10, 8 > min-length)"
    );

    // Default threads (4) run the BAM workflow unordered, so look records up by
    // name instead of assuming input order is preserved on output.
    let by_name: std::collections::HashMap<Vec<u8>, RecordBuf> = out_records
        .iter()
        .map(|r| (r.name().unwrap().to_vec(), r.clone()))
        .collect();

    // read1: "ACGTACGTAC" head-cropped by 2 -> "GTACGTAC" (8 bases), no mods.
    let o1 = &by_name[&b"read1".to_vec()];
    assert_eq!(o1.name().unwrap().to_vec(), b"read1");
    assert_eq!(AsRef::<[u8]>::as_ref(o1.sequence()), b"GTACGTAC");
    assert!(o1.data().get(&Tag::BASE_MODIFICATIONS).is_none());

    // read2: "CCACCCAC" (C at seq idx 0,1,3,4,5,7) head-cropped by 2 -> "ACCCAC" (6 bases).
    // MM "C+m,0,1,0" walks the C occurrences with skip-counts 0,1,0 -> occurrences
    // 0,2,3 -> abs positions 0,3,4 (ML [10,20,30] one-per-position). Window [2,8)
    // keeps abs 3,4 (drops abs 0); renumbered against the window's C positions
    // (3,4,5,7) that's occurrences 0,1 -> deltas [0,0], ML [20,30].
    let o2 = &by_name[&b"read2".to_vec()];
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

#[test]
fn bam_raw_full_window_path_filters_without_rebuilding_records() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("raw_in.bam");
    let out_path = dir.path().join("raw_out.bam");
    write_fixture(&in_path);

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["--min-length", "9", "-t", "4", "-i"])
        .arg(&in_path)
        .arg("-o")
        .arg(&out_path)
        .assert()
        .success();

    let mut reader = bam::io::Reader::new(std::fs::File::open(out_path).unwrap());
    let header = reader.read_header().unwrap();
    let records: Vec<_> = reader
        .record_bufs(&header)
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].name().unwrap().to_vec(), b"read1");
    assert_eq!(records[0].sequence().as_ref(), b"ACGTACGTAC");
    assert_eq!(records[0].quality_scores().as_ref(), &[40; 10]);
}

/// End-to-end: a PacBio-style uBAM with per-base kinetics (`ip`/`pw`, one value
/// per base) must have those arrays sliced in lockstep with the sequence when the
/// read is trimmed — otherwise the output record is invalid (array length != SEQ
/// length) and breaks kinetics/methylation callers.
#[test]
fn bam_to_bam_slices_pacbio_kinetics() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("pb.bam");
    let out_path = dir.path().join("pb_out.bam");

    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(&in_path).unwrap());
    w.write_header(&header).unwrap();
    let mut r = RecordBuf::default();
    *r.flags_mut() = Flags::UNMAPPED;
    *r.name_mut() = Some(b"ccs1".into());
    *r.sequence_mut() = b"ACGTACGTAC".to_vec().into(); // 10 bases
    *r.quality_scores_mut() = vec![40; 10].into();
    r.data_mut().insert(
        Tag::new(b'i', b'p'),
        Value::Array(Array::UInt8((0..10).collect())),
    );
    r.data_mut().insert(
        Tag::new(b'p', b'w'),
        Value::Array(Array::UInt8((100..110).collect())),
    );
    w.write_alignment_record(&header, &r).unwrap();
    w.try_finish().unwrap();

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "--head-crop",
            "3",
            "--tail-crop",
            "2",
            "-i",
        ])
        .arg(&in_path)
        .arg("-o")
        .arg(&out_path)
        .assert()
        .success();

    // window [3,8): seq "TACGT"; ip -> [3,4,5,6,7]; pw -> [103,104,105,106,107].
    let mut rdr = bam::io::Reader::new(std::fs::File::open(&out_path).unwrap());
    let hdr = rdr.read_header().unwrap();
    let mut buf = RecordBuf::default();
    rdr.read_record_buf(&hdr, &mut buf).unwrap();

    assert_eq!(buf.sequence().as_ref(), b"TACGT");
    let ip = match buf.data().get(&Tag::new(b'i', b'p')) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        other => panic!("ip: {other:?}"),
    };
    assert_eq!(
        ip,
        vec![3, 4, 5, 6, 7],
        "ip must be sliced to the trimmed window"
    );
    let pw = match buf.data().get(&Tag::new(b'p', b'w')) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        other => panic!("pw: {other:?}"),
    };
    assert_eq!(pw, vec![103, 104, 105, 106, 107], "pw must be sliced too");
}

/// End-to-end: `--update-moves` slices the ONT `mv` move table and bumps `ts`
/// through the actual binary, so a trimmed read stays signal-mappable for
/// Remora/Clair3 v2 instead of dropping the move table.
#[test]
fn bam_update_moves_slices_move_table() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("mv.bam");
    let out_path = dir.path().join("mv_out.bam");

    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(&in_path).unwrap());
    w.write_header(&header).unwrap();
    let mut r = RecordBuf::default();
    *r.flags_mut() = Flags::UNMAPPED;
    *r.name_mut() = Some(b"read1".into());
    *r.sequence_mut() = b"ACGTAC".to_vec().into(); // 6 bases
    *r.quality_scores_mut() = vec![40; 6].into();
    // stride 2; 6 ones (one per base) at block indices 0,1,3,4,6,7.
    r.data_mut().insert(
        Tag::new(b'm', b'v'),
        Value::Array(Array::Int8(vec![2, 1, 1, 0, 1, 1, 0, 1, 1])),
    );
    r.data_mut().insert(Tag::new(b't', b's'), Value::Int32(10));
    r.data_mut().insert(Tag::new(b'n', b's'), Value::Int32(100));
    w.write_alignment_record(&header, &r).unwrap();
    w.try_finish().unwrap();

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "--update-moves",
            "--head-crop",
            "2",
            "-t",
            "1",
            "-i",
        ])
        .arg(&in_path)
        .arg("-o")
        .arg(&out_path)
        .assert()
        .success();

    let mut rdr = bam::io::Reader::new(std::fs::File::open(&out_path).unwrap());
    let hdr = rdr.read_header().unwrap();
    let mut buf = RecordBuf::default();
    rdr.read_record_buf(&hdr, &mut buf).unwrap();

    assert_eq!(buf.sequence().as_ref(), b"GTAC");
    // mv sliced to blocks [3,8): [stride] + [1,1,0,1,1].
    match buf.data().get(&Tag::new(b'm', b'v')) {
        Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1, 1]),
        other => panic!("mv not sliced: {other:?}"),
    }
    // ts bumped by block_first*stride = 3*2 = 6 -> 16 (any integer width).
    let ts = match buf.data().get(&Tag::new(b't', b's')) {
        Some(Value::Int8(n)) => i64::from(*n),
        Some(Value::Int16(n)) => i64::from(*n),
        Some(Value::Int32(n)) => i64::from(*n),
        other => panic!("ts: {other:?}"),
    };
    assert_eq!(ts, 16, "ts must advance past the trimmed head signal");
}

/// A BGZF-framed BAM piped through stdin is detected as BAM rather than generic
/// gzip and is converted through the chained probe stream.
#[test]
fn bam_on_stdin_without_in_format_is_detected() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("in.bam");
    write_fixture(&in_path);
    let bam_bytes = std::fs::read(&in_path).unwrap();

    // BGZF uses the gzip magic prefix.
    assert_eq!(
        &bam_bytes[..2],
        &[0x1f, 0x8b],
        "BAM fixture must be BGZF (gzip magic)"
    );

    let out_path = dir.path().join("out.fastq");
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        // No --in-format: detection must come from the piped bytes alone.
        .args(["--out-format", "fastq", "-o"])
        .arg(&out_path)
        .write_stdin(bam_bytes)
        .assert()
        .success();

    // Both fixture reads should convert to FASTQ (4 lines each).
    let out = std::fs::read_to_string(&out_path).unwrap();
    let lines = out.lines().count();
    assert_eq!(
        lines, 8,
        "expected 2 FASTQ records (8 lines) from the piped BAM, got:\n{out}"
    );
    assert!(
        out.contains("@read1") && out.contains("@read2"),
        "missing reads: {out}"
    );
}
