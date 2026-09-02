//! End-to-end uBAM coverage through the compiled binary, including header,
//! reader, workflow, writer, and CLI integration.

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};

/// Writes a two-read uBAM: a plain read and a read with `MM`/`ML`/`MN` mods.
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

    // Read back the output BAM: check @PG provenance and reconstructed records.
    let mut reader = bam::io::Reader::new(std::fs::File::open(&out_path).unwrap());
    let header = reader.read_header().unwrap();

    assert!(
        header
            .programs()
            .roots()
            .any(|(id, _)| AsRef::<[u8]>::as_ref(id) == b"whittle"),
        "Expected an @PG record with ID whittle in the output header, got {:?}",
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
        "Both reads survive --head-crop 2 (lengths 10 and 8, above min-length)"
    );

    // The BAM workflow emits records unordered at `threads > 1`, so records are
    // looked up by name rather than by input order.
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
    // MM "C+m,0,1,0" walks the C occurrences with skip counts 0,1,0 -> occurrences
    // 0,2,3 -> abs positions 0,3,4 (ML [10,20,30], one per position). Window [2,8)
    // keeps abs 3,4 (drops abs 0); renumbered against the window's C positions
    // (3,4,5,7) that is occurrences 0,1 -> deltas [0,0], ML [20,30].
    let o2 = &by_name[&b"read2".to_vec()];
    assert_eq!(o2.name().unwrap().to_vec(), b"read2");
    assert_eq!(AsRef::<[u8]>::as_ref(o2.sequence()), b"ACCCAC");
    let mm = match o2.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => s.to_vec(),
        other => panic!("Expected MM tag, got {other:?}"),
    };
    assert_eq!(mm, b"C+m,0,0;");
    let ml = match o2.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        other => panic!("Expected ML tag, got {other:?}"),
    };
    assert_eq!(ml, vec![20, 30]);
    let mn = match o2.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
        Some(Value::Int32(n)) => *n,
        other => panic!("Expected MN tag, got {other:?}"),
    };
    assert_eq!(mn, 6, "MN equals the output segment length");
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

/// A PacBio-style uBAM with per-base kinetics (`ip`/`pw`, one value per base)
/// has those arrays sliced in lockstep with the sequence when the read is
/// trimmed. Otherwise the output record is invalid (array length differs from
/// SEQ length) and breaks kinetics and methylation callers.
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

    // Window [3,8): seq "TACGT"; ip -> [3,4,5,6,7]; pw -> [103,104,105,106,107].
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

/// `--update-moves` slices the ONT `mv` move table and advances `ts` through
/// the compiled binary, so a trimmed read stays signal-mappable for Remora and
/// Clair3 v2 instead of dropping the move table.
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
    // Stride 2; 6 ones (one per base) at block indices 0,1,3,4,6,7.
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
    // `mv` sliced to blocks [3,8): the stride, then [1,1,0,1,1].
    match buf.data().get(&Tag::new(b'm', b'v')) {
        Some(Value::Array(Array::Int8(v))) => assert_eq!(v, &[2, 1, 1, 0, 1, 1]),
        other => panic!("mv not sliced: {other:?}"),
    }
    // `ts` advances by block_first * stride = 3 * 2 = 6, to 16 (any integer
    // width).
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

    // Both fixture reads convert to FASTQ (4 lines each).
    let out = std::fs::read_to_string(&out_path).unwrap();
    let lines = out.lines().count();
    assert_eq!(
        lines, 8,
        "Expected 2 FASTQ records (8 lines) from the piped BAM, got:\n{out}"
    );
    assert!(
        out.contains("@read1") && out.contains("@read2"),
        "Missing reads: {out}"
    );
}

/// Writes a one-read uBAM with `MM`/`ML`/`MN`, per-base kinetics and the aux
/// tags a barcoded dorado record carries. `bi` is written only when `barcoded`
/// is set, so the same read can be trimmed by `--trim-barcodes` and by the
/// equivalent fixed crop.
///
/// Sequence: 20 bases, `CCACCCACGTCCACCCACGT`, C at 0, 1, 3, 4, 5, 7, 10, 11,
/// 13, 14, 15, 17. `C+m,4,1,2;` marks C occurrences 4, 6 and 9, positions 5, 10
/// and 14, all inside the barcode window.
fn write_barcode_fixture(path: &std::path::Path, barcoded: bool) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut r = RecordBuf::default();
    *r.flags_mut() = Flags::UNMAPPED;
    *r.name_mut() = Some(b"read1".into());
    *r.sequence_mut() = b"CCACCCACGTCCACCCACGT".to_vec().into();
    *r.quality_scores_mut() = vec![40; 20].into();
    let data = r.data_mut();
    data.insert(
        Tag::BASE_MODIFICATIONS,
        Value::String(b"C+m,4,1,2;".to_vec().into()),
    );
    data.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(vec![10, 20, 30])),
    );
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(20));
    data.insert(
        Tag::new(b'i', b'p'),
        Value::Array(Array::UInt8((0..20).collect())),
    );
    data.insert(
        Tag::new(b'p', b'w'),
        Value::Array(Array::UInt8((100..120).collect())),
    );
    data.insert(
        Tag::new(b'B', b'C'),
        Value::String(b"barcode07".as_slice().into()),
    );
    data.insert(Tag::new(b'b', b'v'), Value::String(b"v5".as_slice().into()));
    if barcoded {
        // Front barcode over bases [0, 4], rear barcode over [15, 18): dorado
        // writes [score, front_start, front_len, front_score, rear_end,
        // rear_len, rear_score], so the kept window is [0 + 4 + 1, 18 - 3).
        data.insert(
            Tag::new(b'b', b'i'),
            Value::Array(Array::Float(vec![90.0, 0.0, 4.0, 88.0, 18.0, 3.0, 87.0])),
        );
    }
    w.write_alignment_record(&header, &r).unwrap();
    w.try_finish().unwrap();
}

/// Reads the single record of a uBAM.
fn read_one(path: &std::path::Path) -> RecordBuf {
    let mut rdr = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let hdr = rdr.read_header().unwrap();
    let mut buf = RecordBuf::default();
    assert!(rdr.read_record_buf(&hdr, &mut buf).unwrap() > 0);
    buf
}

/// `--trim-barcodes` removes the spans `bi` records through the same machinery
/// as a fixed crop, so the modification calls and the per-base kinetics land on
/// exactly the bases an equivalent `--head-crop`/`--tail-crop` run leaves. `bi`
/// itself goes, since its positions index the untrimmed read; the barcode call
/// (`BC`, `bv`) is a per-read label and stays.
#[test]
fn trim_barcodes_matches_the_equivalent_crop_and_keeps_tags_in_register() {
    let dir = tempfile::tempdir().unwrap();
    let barcoded = dir.path().join("barcoded.bam");
    let plain = dir.path().join("plain.bam");
    write_barcode_fixture(&barcoded, true);
    write_barcode_fixture(&plain, false);

    let by_tag = dir.path().join("by_tag.bam");
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["--trim-barcodes", "-t", "1", "-i"])
        .arg(&barcoded)
        .arg("-o")
        .arg(&by_tag)
        .assert()
        .success();

    // The window `bi` describes is [5, 15), the same as cropping 5 from the
    // head and 5 from the tail of a 20-base read.
    let by_crop = dir.path().join("by_crop.bam");
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["--head-crop", "5", "--tail-crop", "5", "-t", "1", "-i"])
        .arg(&plain)
        .arg("-o")
        .arg(&by_crop)
        .assert()
        .success();

    let tagged = read_one(&by_tag);
    let cropped = read_one(&by_crop);

    assert_eq!(tagged.sequence().as_ref(), b"CACGTCCACC");
    assert_eq!(
        tagged.sequence().as_ref(),
        cropped.sequence().as_ref(),
        "The barcode window must cut where the equivalent crop cuts"
    );
    assert_eq!(
        tagged.quality_scores().as_ref(),
        cropped.quality_scores().as_ref()
    );
    for t in [
        Tag::BASE_MODIFICATIONS,
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
        Tag::new(b'i', b'p'),
        Tag::new(b'p', b'w'),
    ] {
        assert_eq!(
            tagged.data().get(&t),
            cropped.data().get(&t),
            "Tag {t:?} must land on the same bases as under the equivalent crop"
        );
    }
    // The window's C positions are 0, 2, 5, 6, 8 and 9; the three modified
    // bases land on occurrences 0, 2 and 5, renumbered to deltas 0, 1 and 2.
    assert_eq!(
        tagged.data().get(&Tag::BASE_MODIFICATIONS),
        Some(&Value::String(b"C+m,0,1,2;".to_vec().into()))
    );
    assert_eq!(
        tagged.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES),
        Some(&Value::Array(Array::UInt8(vec![10, 20, 30])))
    );
    assert_eq!(
        tagged.data().get(&Tag::new(b'i', b'p')),
        Some(&Value::Array(Array::UInt8((5..15).collect())))
    );

    assert!(
        tagged.data().get(&Tag::new(b'b', b'i')).is_none(),
        "bi indexes the untrimmed read and must be dropped"
    );
    assert_eq!(
        tagged.data().get(&Tag::new(b'B', b'C')),
        Some(&Value::String(b"barcode07".as_slice().into()))
    );
    assert_eq!(
        tagged.data().get(&Tag::new(b'b', b'v')),
        Some(&Value::String(b"v5".as_slice().into()))
    );
}

/// `--head-crop` counts from the first base after the front barcode, since the
/// barcode stage is the outermost one, and the JSON summary records the
/// resulting base counts.
#[test]
fn trim_barcodes_runs_before_the_crop_and_updates_json_counts() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("barcoded.bam");
    let out_path = dir.path().join("out.bam");
    let summary = dir.path().join("qc.json");
    write_barcode_fixture(&in_path, true);

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["--trim-barcodes", "--head-crop", "2", "-t", "1", "-i"])
        .arg(&in_path)
        .arg("-o")
        .arg(&out_path)
        .arg("--summary-json")
        .arg(&summary)
        .assert()
        .success();

    // Window [5, 15) less a head crop of 2 leaves [7, 15).
    assert_eq!(read_one(&out_path).sequence().as_ref(), b"CGTCCACC");

    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&summary).unwrap()).unwrap();
    assert_eq!(v["params"]["trim_barcodes"], true);
    assert_eq!(v["warnings"]["barcode_tag_malformed_reads"], 0);
    assert_eq!(v["bases"]["input"].as_u64().unwrap(), 20);
    assert_eq!(v["bases"]["output"].as_u64().unwrap(), 8);
}

/// A `bi` that does not describe a window inside the read leaves the read
/// untrimmed and is counted once, rather than failing the run.
#[test]
fn a_malformed_barcode_tag_is_counted_and_the_read_survives() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("bad.bam");
    let out_path = dir.path().join("out.bam");
    let summary = dir.path().join("qc.json");

    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(&in_path).unwrap());
    w.write_header(&header).unwrap();
    let mut r = RecordBuf::default();
    *r.flags_mut() = Flags::UNMAPPED;
    *r.name_mut() = Some(b"read1".into());
    *r.sequence_mut() = b"ACGTACGTAC".to_vec().into();
    *r.quality_scores_mut() = vec![40; 10].into();
    // A rear end past the sequence.
    r.data_mut().insert(
        Tag::new(b'b', b'i'),
        Value::Array(Array::Float(vec![90.0, 0.0, 2.0, 88.0, 40.0, 3.0, 87.0])),
    );
    w.write_alignment_record(&header, &r).unwrap();
    w.try_finish().unwrap();

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["--trim-barcodes", "-t", "1", "-i"])
        .arg(&in_path)
        .arg("-o")
        .arg(&out_path)
        .arg("--summary-json")
        .arg(&summary)
        .assert()
        .success();

    let out = read_one(&out_path);
    assert_eq!(out.sequence().as_ref(), b"ACGTACGTAC");
    assert!(
        out.data().get(&Tag::new(b'b', b'i')).is_some(),
        "An untrimmed read keeps its bi tag"
    );
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&summary).unwrap()).unwrap();
    assert_eq!(v["warnings"]["barcode_tag_malformed_reads"], 1);
}

/// The barcode spans come from a BAM aux tag, so the flag is refused on FASTQ
/// input rather than accepted and ignored.
#[test]
fn trim_barcodes_is_rejected_on_fastq_input() {
    let dir = tempfile::tempdir().unwrap();
    let in_path = dir.path().join("reads.fastq");
    std::fs::write(&in_path, "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n").unwrap();

    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["--trim-barcodes", "-i"])
        .arg(&in_path)
        .arg("-o")
        .arg(dir.path().join("out.fastq"))
        .assert()
        .failure()
        .stderr(predicates::str::contains("--trim-barcodes"));

    // A stream carries no extension, so the refusal comes from the detected
    // format instead.
    Command::cargo_bin("whittle")
        .unwrap()
        .env_remove("WHITTLE_LOG")
        .args(["--trim-barcodes", "--out-format", "fastq", "-o"])
        .arg(dir.path().join("out2.fastq"))
        .write_stdin("@r1\nACGTACGTAC\n+\nIIIIIIIIII\n")
        .assert()
        .failure()
        .stderr(predicates::str::contains("--trim-barcodes"));
}
