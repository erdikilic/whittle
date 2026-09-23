use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

use super::*;
use crate::config::FastqTags;
use crate::trim::{QualityOp, TrimPlan};

mod build;
mod drivers;
mod mods;
mod naming;
mod signal;
mod tags;
mod to_fastq;

/// Builds one output record for interval `[start, end)`, segment `idx` of
/// `total`, with the modification block classified from `src` and no tag
/// removal: the record is encoded, rebuilt from its raw bytes and decoded.
pub(super) fn reconstruct_record(
    src: &RecordBuf,
    start: usize,
    end: usize,
    total: usize,
    idx: usize,
    update_moves: bool,
) -> RecordBuf {
    let seq = src.sequence().as_ref();
    let mod_block = inspect_mod_block(src, seq.len());
    let window = Window {
        start,
        end,
        idx,
        total,
    };
    let keep = TagRemoval::default();
    let moves = update_moves.then(|| MoveIndex::new(src, false)).flatten();
    match window_edit(src, window, mod_block, None, moves.as_ref(), &keep) {
        Some(edit) => decode_built(&build_record(&raw_record(src), edit).unwrap()),
        None => src.clone(),
    }
}

/// Decodes a record `build_record` produced.
pub(super) fn decode_built(bytes: &[u8]) -> RecordBuf {
    let mut framed = u32::try_from(bytes.len()).unwrap().to_le_bytes().to_vec();
    framed.extend_from_slice(bytes);
    let mut reader = bam::io::Reader::from(framed.as_slice());
    let mut raw = bam::Record::default();
    assert_ne!(reader.read_record(&mut raw).unwrap(), 0);
    decode_raw_record(&raw).unwrap()
}

fn ubam_with_mods(seq: &[u8], quals: Vec<u8>, mm: &[u8], ml: Vec<u8>) -> RecordBuf {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = seq.to_vec().into();
    *rec.quality_scores_mut() = quals.into();
    let data = rec.data_mut();
    data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.to_vec().into()));
    data.insert(
        Tag::BASE_MODIFICATION_PROBABILITIES,
        Value::Array(Array::UInt8(ml)),
    );
    rec
}

pub(super) fn cfg_bam2fq(quality: Option<QualityOp>, head: usize, tags: FastqTags) -> Config {
    Config {
        trim: TrimPlan {
            head,
            tail: 0,
            quality,
        },
        fastq_tags: tags,
        quiet: true,
        ..Config::default()
    }
}

/// `CCACCCAC` has C at 0, 1, 3, 4, 5, 7; `C+m,0,1,0` marks occurrences 0, 2,
/// 3 (positions 0, 3, 4) with ML [10, 20, 30]. A head crop of 2 (window
/// [2,8)) keeps positions 3 and 4, renumbered to `C+m,0,0;` with ML [20, 30]
/// and MN 6.
fn read2_with_mods_and_rg() -> RecordBuf {
    let mut rec = ubam_with_mods(b"CCACCCAC", vec![35; 8], b"C+m,0,1,0;", vec![10, 20, 30]);
    rec.data_mut()
        .insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
    rec.data_mut()
        .insert(Tag::READ_GROUP, Value::String(b"grp1".as_slice().into()));
    rec
}

/// A record whose modification block cannot be placed on its sequence: the
/// variants `malformed_mod_blocks` feeds through both output paths.
fn malformed_mod_record(variant: &str) -> RecordBuf {
    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"r1".into());
    *rec.sequence_mut() = b"CCCA".to_vec().into();
    *rec.quality_scores_mut() = vec![40; 4].into();
    let d = rec.data_mut();
    let (mm, ml, mn) = match variant {
        // MM declares 3 positions, ML has 1 byte.
        "ml_short" => (&b"C+m,0,0,0;"[..], Value::Array(Array::UInt8(vec![5])), 4),
        // MM stops parsing at the `x`.
        "mm_garbled" => (
            &b"C+m,0,0x,0;"[..],
            Value::Array(Array::UInt8(vec![5, 6, 7])),
            4,
        ),
        // ML at a subtype other than B:C.
        "ml_subtype" => (
            &b"C+m,0,0,0;"[..],
            Value::Array(Array::Int8(vec![5, 6, 7])),
            4,
        ),
        // MN disagrees with the sequence length.
        "mn_mismatch" => (
            &b"C+m,0,0,0;"[..],
            Value::Array(Array::UInt8(vec![5, 6, 7])),
            40,
        ),
        other => panic!("Unknown variant {other}"),
    };
    d.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.to_vec().into()));
    d.insert(Tag::BASE_MODIFICATION_PROBABILITIES, ml);
    d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(mn));
    d.insert(Tag::READ_GROUP, Value::String(b"grp".as_slice().into()));
    rec
}

const MALFORMED_MOD_VARIANTS: [&str; 4] = ["ml_short", "mm_garbled", "ml_subtype", "mn_mismatch"];

/// A synthetic move table: stride 2, 6 ones (one per base) at block indexes
/// 0, 1, 3, 4, 6, 7, 8 blocks in total. Shared by the `--update-moves` tests.
fn ubam_with_moves() -> RecordBuf {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(b"r1".into());
    *src.sequence_mut() = b"ACGTAC".to_vec().into();
    *src.quality_scores_mut() = vec![40; 6].into();
    let d = src.data_mut();
    d.insert(
        Tag::new(b'm', b'v'),
        Value::Array(Array::Int8(vec![2, 1, 1, 0, 1, 1, 0, 1, 1])),
    );
    d.insert(Tag::new(b't', b's'), Value::Int32(10));
    // Consistent: ts0 + n_blocks*stride = 10 + 8*2 = 26.
    d.insert(Tag::new(b'n', b's'), Value::Int32(26));
    d.insert(
        Tag::new(b's', b't'),
        Value::String(b"2024-06-21T10:00:00Z".as_slice().into()),
    );
    d.insert(Tag::new(b'd', b'u'), Value::Float(5.0));
    src
}

/// Encodes a decoded record as the raw BAM record a production reader yields.
fn raw_record(rec: &RecordBuf) -> bam::Record {
    use noodles_sam::alignment::io::Write as _;

    let header = sam::Header::default();
    let mut bytes = Vec::new();
    {
        let mut w = bam::io::Writer::new(&mut bytes);
        w.write_header(&header).unwrap();
        w.write_alignment_record(&header, rec).unwrap();
        w.try_finish().unwrap();
    }
    let mut r = bam::io::Reader::new(bytes.as_slice());
    r.read_header().unwrap();
    let mut raw = bam::Record::default();
    assert_ne!(r.read_record(&mut raw).unwrap(), 0);
    raw
}

/// A 10-base PacBio-style record: `fi[i]` belongs to base `i`, `ri[i]` to
/// base `9 - i` (reverse kinetics are stored last base first).
fn pacbio_kinetics_record() -> RecordBuf {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(b"r1".into());
    *src.sequence_mut() = b"ACGTACGTAC".to_vec().into();
    *src.quality_scores_mut() = vec![40; 10].into();
    let d = src.data_mut();
    d.insert(
        Tag::new(b'f', b'i'),
        Value::Array(Array::UInt8((0..10).collect())),
    );
    d.insert(
        Tag::new(b'r', b'i'),
        Value::Array(Array::UInt8((10..20).collect())),
    );
    src
}

fn u8_array(rec: &RecordBuf, tag: [u8; 2]) -> Vec<u8> {
    match rec.data().get(&Tag::new(tag[0], tag[1])) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        other => panic!("{}: {other:?}", std::str::from_utf8(&tag).unwrap()),
    }
}

const HIFI: &[u8] = b"m64011_190830_220126/123/ccs";

/// Catalog barcode BC01.
pub(super) const BC01: &[u8] = b"AAGAAAGTTGTCGGTGTCTTTGTG";

/// A 10-base PacBio HiFi record: `qs`/`qe` query coordinates 100..110,
/// `rn` passes, the `sa` coverage runs (5 over bases 0..4, 7 over 4..7,
/// 5 over 7..10), the `sm`/`sx` per-base counts and the `ds`/`ls` undo
/// blobs.
fn pacbio_hifi_record(name: &[u8]) -> RecordBuf {
    let mut src = RecordBuf::default();
    *src.flags_mut() = Flags::UNMAPPED;
    *src.name_mut() = Some(name.into());
    *src.sequence_mut() = b"ACGTACGTAC".to_vec().into();
    *src.quality_scores_mut() = vec![40; 10].into();
    let d = src.data_mut();
    d.insert(Tag::new(b'q', b's'), Value::Int32(100));
    d.insert(Tag::new(b'q', b'e'), Value::Int32(110));
    d.insert(Tag::new(b'r', b'n'), Value::Int32(3));
    d.insert(
        Tag::new(b's', b'a'),
        Value::Array(Array::UInt32(vec![4, 5, 3, 7, 3, 5])),
    );
    d.insert(
        Tag::new(b's', b'm'),
        Value::Array(Array::UInt8((0..10).collect())),
    );
    d.insert(
        Tag::new(b's', b'x'),
        Value::Array(Array::UInt8((10..20).collect())),
    );
    d.insert(
        Tag::new(b'd', b's'),
        Value::Array(Array::UInt8(vec![1, 2, 3])),
    );
    d.insert(
        Tag::new(b'l', b's'),
        Value::Array(Array::UInt8(vec![4, 5, 6])),
    );
    src
}

/// A dorado record: `ubam_with_moves` (`mv`, `ts` 10, `ns` 26, `st`, `du`
/// 5 s) plus `qs:f`, `me` and `er`, without `pi`.
fn ont_record() -> RecordBuf {
    let mut src = ubam_with_moves();
    let d = src.data_mut();
    d.insert(Tag::new(b'q', b's'), Value::Float(20.0));
    d.insert(Tag::new(b'm', b'e'), Value::Int32(12));
    d.insert(
        Tag::new(b'e', b'r'),
        Value::String(b"signal_positive".as_slice().into()),
    );
    src
}

fn tag(rec: &RecordBuf, t: [u8; 2]) -> Option<Value> {
    rec.data().get(&Tag::new(t[0], t[1])).cloned()
}

fn string_value(s: &[u8]) -> Option<Value> {
    Some(Value::String(s.into()))
}

fn name_of(rec: &RecordBuf) -> Vec<u8> {
    rec.name().map(|n| n.to_vec()).unwrap_or_default()
}

/// Splits at a low-quality base and carries every tag.
fn split_cfg() -> Config {
    cfg_bam2fq(
        Some(QualityOp::Split {
            cutoff: 20,
            window: 1,
        }),
        0,
        FastqTags::All,
    )
}

/// Runs the BAM-to-FASTQ workflow and returns the stats and the text.
fn bam2fq(recs: Vec<RecordBuf>, cfg: &Config) -> (Stats, String) {
    let mut out = Vec::new();
    let stats = run_bam_to_fastq(
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut out,
        cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    (stats, String::from_utf8(out).unwrap())
}

/// Runs the BAM-to-BAM workflow and returns the stats and the decoded
/// output records.
fn bam2bam(recs: Vec<RecordBuf>, cfg: &Config) -> (Stats, Vec<RecordBuf>) {
    let header = sam::Header::default();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("o.bam");
    let mut sink = crate::io::bam::writer(Some(&path), &header, false, 6).unwrap();
    let stats = run_bam(
        &header,
        recs.iter().map(|r| Ok(raw_record(r))),
        &mut sink,
        cfg,
        &Arc::new(Counters::default()),
    )
    .unwrap();
    sink.finish().unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let mut reader = noodles_bam::io::Reader::new(bytes.as_slice());
    let h = reader.read_header().unwrap();
    let mut out = Vec::new();
    let mut buf = RecordBuf::default();
    while reader.read_record_buf(&h, &mut buf).unwrap() != 0 {
        out.push(buf.clone());
    }
    (stats, out)
}

/// Builds the removal set the flags produce.
fn removal(tags: &[&str]) -> crate::config::TagRemoval {
    let tags: Vec<String> = tags.iter().map(|t| (*t).to_string()).collect();
    crate::config::TagRemoval::parse(&tags).unwrap()
}

/// An 8-base record carrying a rewritten block (`MM`/`ML`/`MN`), a per-base
/// array (`ip`), a run-length coverage array (`sa`) and a copied scalar
/// (`RG`).
fn record_with_mixed_tags() -> RecordBuf {
    let mut rec = ubam_with_mods(b"CCACCCAC", vec![40; 8], b"C+m,0,1,0;", vec![10, 20, 30]);
    let data = rec.data_mut();
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
    data.insert(
        Tag::new(b'i', b'p'),
        Value::Array(Array::UInt8((0..8).collect())),
    );
    data.insert(
        Tag::new(b's', b'a'),
        Value::Array(Array::UInt32(vec![8, 3])),
    );
    data.insert(
        Tag::new(b'R', b'G'),
        Value::String(b"run1".as_slice().into()),
    );
    rec
}
