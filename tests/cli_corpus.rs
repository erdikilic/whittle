use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

use assert_cmd::Command;
use flate2::Compression;
use flate2::write::GzEncoder;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::{self as sam};

const READS: usize = 1024;
const ADAPTER: &[u8] = b"GGGGTTTTGGGGTTTT";

#[derive(Debug, Clone)]
struct SourceRead {
    id: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    rg: Option<String>,
    mods: Option<ModFixture>,
    ip: Option<Vec<u8>>,
    moves: bool,
    ts: i32,
    adapter: Option<AdapterCase>,
}

#[derive(Debug, Clone)]
struct ModFixture {
    abs: Vec<usize>,
    probs: Vec<u8>,
}

#[derive(Debug, Clone, Copy)]
enum AdapterCase {
    FivePrime,
    ThreePrime,
    Interior { start: usize },
}

#[derive(Debug, Clone, Copy)]
enum QualityOp {
    None,
    Trim(u8),
    Split { cutoff: u8, window: usize },
}

#[derive(Debug, Clone, Copy)]
struct ExpectCfg {
    head: usize,
    tail: usize,
    min_len: usize,
    max_len: usize,
    min_qual: f64,
    min_gc: Option<f64>,
    max_gc: Option<f64>,
    quality: QualityOp,
    adapters: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct FastqRecord {
    head: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct BamRecord {
    name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    rg: Option<String>,
    mm: Option<String>,
    ml: Option<Vec<u8>>,
    mn: Option<i32>,
    ip: Option<Vec<u8>>,
    mv: Option<Vec<i8>>,
    ts: Option<i64>,
    ns: Option<i64>,
}

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

fn splitmix_base(read_idx: usize, pos: usize) -> u8 {
    let mut z = 0x9E37_79B9_7F4A_7C15u64
        .wrapping_add(read_idx as u64)
        .wrapping_add((pos as u64) << 32);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    b"ACGT"[(z & 0b11) as usize]
}

fn build_seq(i: usize, len: usize) -> Vec<u8> {
    match i % 16 {
        0 => (0..len)
            .map(|j| if j % 5 == 0 { b'T' } else { b'A' })
            .collect(),
        1 => (0..len)
            .map(|j| if j % 2 == 0 { b'C' } else { b'G' })
            .collect(),
        _ => (0..len).map(|j| splitmix_base(i, j)).collect(),
    }
}

fn plant_adapter(i: usize, mut seq: Vec<u8>) -> (Vec<u8>, Option<AdapterCase>) {
    match i % 32 {
        8 => {
            let mut out = ADAPTER.to_vec();
            out.extend_from_slice(&seq);
            (out, Some(AdapterCase::FivePrime))
        },
        9 => {
            seq.extend_from_slice(ADAPTER);
            (seq, Some(AdapterCase::ThreePrime))
        },
        10 => {
            let start = 32 + (i % 5);
            seq.splice(start..start, ADAPTER.iter().copied());
            (seq, Some(AdapterCase::Interior { start }))
        },
        _ => (seq, None),
    }
}

fn build_qual(i: usize, len: usize) -> Vec<u8> {
    let mut qual = vec![35u8; len];
    match i % 16 {
        2 => qual.fill(12),
        3 => {
            for q in qual.iter_mut().take(8) {
                *q = 5;
            }
            for q in qual.iter_mut().skip(len.saturating_sub(8)) {
                *q = 5;
            }
        },
        4 => {
            for q in qual.iter_mut().take(35).skip(30) {
                *q = 5;
            }
            for q in qual.iter_mut().take(72).skip(65) {
                *q = 5;
            }
        },
        _ => {},
    }
    qual
}

fn c_positions(seq: &[u8]) -> Vec<usize> {
    seq.iter()
        .enumerate()
        .filter(|(_, b)| matches!(b, b'C' | b'c'))
        .map(|(i, _)| i)
        .collect()
}

fn build_mod_fixture(i: usize, seq: &[u8]) -> Option<ModFixture> {
    if !i.is_multiple_of(3) {
        return None;
    }
    let cs = c_positions(seq);
    if cs.len() < 3 {
        return None;
    }
    let occ = [0usize, cs.len() / 2, cs.len() - 1];
    Some(ModFixture {
        abs: occ.into_iter().map(|idx| cs[idx]).collect(),
        probs: vec![
            20 + (i % 40) as u8,
            90 + (i % 50) as u8,
            170 + (i % 60) as u8,
        ],
    })
}

fn corpus() -> Vec<SourceRead> {
    (0..READS)
        .map(|i| {
            let len = 90 + ((i * 37) % 41);
            let (seq, adapter) = plant_adapter(i, build_seq(i, len));
            let qual = build_qual(i, seq.len());
            let mods = build_mod_fixture(i, &seq);
            let ip = (i % 5 == 0).then(|| (0..seq.len()).map(|j| ((i + j) & 0xff) as u8).collect());
            SourceRead {
                id: format!("r{i:04}"),
                seq,
                qual,
                rg: (i % 4 == 0).then(|| format!("rg{}", i % 3)),
                mods,
                ip,
                moves: i % 7 == 0,
                ts: 1000 + i as i32,
                adapter,
            }
        })
        .collect()
}

fn mod_mm_for_source(read: &SourceRead) -> Option<String> {
    let mods = read.mods.as_ref()?;
    let cs = c_positions(&read.seq);
    let mut prev = -1isize;
    let mut mm = String::from("C+m");
    for &abs in &mods.abs {
        let occ = cs
            .iter()
            .position(|&p| p == abs)
            .expect("fixture mod must land on a C") as isize;
        let delta = occ - prev - 1;
        mm.push_str(&format!(",{delta}"));
        prev = occ;
    }
    mm.push(';');
    Some(mm)
}

fn write_fastq_gz(path: &Path, reads: &[SourceRead]) {
    let mut enc = GzEncoder::new(std::fs::File::create(path).unwrap(), Compression::default());
    for read in reads {
        enc.write_all(b"@").unwrap();
        enc.write_all(read.id.as_bytes()).unwrap();
        enc.write_all(b"\n").unwrap();
        enc.write_all(&read.seq).unwrap();
        enc.write_all(b"\n+\n").unwrap();
        let ascii: Vec<u8> = read.qual.iter().map(|q| q + 33).collect();
        enc.write_all(&ascii).unwrap();
        enc.write_all(b"\n").unwrap();
    }
    enc.finish().unwrap();
}

fn write_bam(path: &Path, reads: &[SourceRead]) {
    let header = sam::Header::default();
    let mut writer = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for read in reads {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(read.id.as_bytes().into());
        *rec.sequence_mut() = read.seq.clone().into();
        *rec.quality_scores_mut() = read.qual.clone().into();

        if let Some(rg) = &read.rg {
            rec.data_mut()
                .insert(Tag::new(b'R', b'G'), Value::String(rg.as_bytes().into()));
        }
        if let Some(mm) = mod_mm_for_source(read) {
            let mods = read.mods.as_ref().unwrap();
            rec.data_mut().insert(
                Tag::BASE_MODIFICATIONS,
                Value::String(mm.into_bytes().into()),
            );
            rec.data_mut().insert(
                Tag::BASE_MODIFICATION_PROBABILITIES,
                Value::Array(Array::UInt8(mods.probs.clone())),
            );
            rec.data_mut().insert(
                Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
                Value::Int32(read.seq.len() as i32),
            );
        }
        if let Some(ip) = &read.ip {
            rec.data_mut()
                .insert(Tag::new(b'i', b'p'), Value::Array(Array::UInt8(ip.clone())));
        }
        if read.moves {
            let mut mv = Vec::with_capacity(read.seq.len() + 1);
            mv.push(2);
            mv.extend(std::iter::repeat_n(1, read.seq.len()));
            rec.data_mut()
                .insert(Tag::new(b'm', b'v'), Value::Array(Array::Int8(mv)));
            rec.data_mut()
                .insert(Tag::new(b't', b's'), Value::Int32(read.ts));
            rec.data_mut().insert(
                Tag::new(b'n', b's'),
                Value::Int32(read.ts + (read.seq.len() as i32 * 2)),
            );
        }

        writer.write_alignment_record(&header, &rec).unwrap();
    }
    writer.try_finish().unwrap();
}

fn trim_edge(qual: &[u8], cutoff: u8) -> Vec<(usize, usize)> {
    let mut start = 0usize;
    while start < qual.len() && qual[start] < cutoff {
        start += 1;
    }
    let mut end = qual.len();
    while end > start && qual[end - 1] < cutoff {
        end -= 1;
    }
    if start < end {
        vec![(start, end)]
    } else {
        Vec::new()
    }
}

fn split_low_quality(qual: &[u8], cutoff: u8, window: usize) -> Vec<(usize, usize)> {
    let window = window.max(1);
    let mut out = Vec::new();
    let mut segment_start = None;
    let mut last_good = None;
    let mut bad_run = 0usize;

    for (i, &q) in qual.iter().enumerate() {
        if q >= cutoff {
            segment_start.get_or_insert(i);
            last_good = Some(i);
            bad_run = 0;
        } else {
            bad_run += 1;
            if bad_run >= window {
                if let (Some(s), Some(e)) = (segment_start, last_good)
                    && e + 1 > s
                {
                    out.push((s, e + 1));
                }
                segment_start = None;
                last_good = None;
            }
        }
    }
    if let (Some(s), Some(e)) = (segment_start, last_good)
        && e + 1 > s
    {
        out.push((s, e + 1));
    }
    out
}

fn adapter_windows(read: &SourceRead, start: usize, end: usize) -> Vec<(usize, usize)> {
    let Some(adapter) = read.adapter else {
        return vec![(start, end)];
    };

    let ad_len = ADAPTER.len();
    match adapter {
        AdapterCase::FivePrime if start == 0 && end > ad_len => vec![(ad_len, end)],
        AdapterCase::ThreePrime if end == read.seq.len() && end.saturating_sub(start) > ad_len => {
            vec![(start, end - ad_len)]
        },
        AdapterCase::Interior { start: cut_start }
            if cut_start >= start && cut_start + ad_len <= end =>
        {
            let cut_end = cut_start + ad_len;
            [(start, cut_start), (cut_end, end)]
                .into_iter()
                .filter(|(s, e)| e > s)
                .collect()
        },
        _ => vec![(start, end)],
    }
}

fn produced_windows(read: &SourceRead, cfg: ExpectCfg) -> Vec<(usize, usize)> {
    let start = cfg.head.min(read.seq.len());
    let end = read.seq.len().saturating_sub(cfg.tail).max(start);
    if start >= end {
        return Vec::new();
    }

    let adapter_segs = if cfg.adapters {
        adapter_windows(read, start, end)
    } else {
        vec![(start, end)]
    };

    let mut out = Vec::new();
    for (seg_start, seg_end) in adapter_segs {
        let inner = match cfg.quality {
            QualityOp::None => vec![(0, seg_end - seg_start)],
            QualityOp::Trim(cutoff) => trim_edge(&read.qual[seg_start..seg_end], cutoff),
            QualityOp::Split { cutoff, window } => {
                split_low_quality(&read.qual[seg_start..seg_end], cutoff, window)
            },
        };
        out.extend(
            inner
                .into_iter()
                .map(|(s, e)| (s + seg_start, e + seg_start)),
        );
    }
    out
}

fn mean_prob_q(qual: &[u8]) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum: f64 = qual
        .iter()
        .map(|q| 10_f64.powf(f64::from(*q) / -10.0))
        .sum();
    (sum / qual.len() as f64).log10() * -10.0
}

fn gc_fraction(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|b| matches!(b, b'G' | b'g' | b'C' | b'c'))
        .count();
    gc as f64 / seq.len() as f64
}

fn passes(seq: &[u8], qual: &[u8], cfg: ExpectCfg) -> bool {
    if seq.is_empty() || seq.len() < cfg.min_len || seq.len() > cfg.max_len {
        return false;
    }
    if cfg.min_qual > 0.0 && mean_prob_q(qual) < cfg.min_qual {
        return false;
    }
    if let Some(min_gc) = cfg.min_gc
        && gc_fraction(seq) < min_gc
    {
        return false;
    }
    if let Some(max_gc) = cfg.max_gc
        && gc_fraction(seq) > max_gc
    {
        return false;
    }
    true
}

fn output_name(id: &str, total: usize, idx: usize) -> String {
    if total > 1 {
        format!("{id}_segment_{}", idx + 1)
    } else {
        id.to_string()
    }
}

fn expected_mods(read: &SourceRead, start: usize, end: usize) -> Option<(String, Vec<u8>, i32)> {
    let mods = read.mods.as_ref()?;
    let window_cs: Vec<_> = c_positions(&read.seq)
        .into_iter()
        .filter(|p| *p >= start && *p < end)
        .collect();
    let mut prev = -1isize;
    let mut deltas = Vec::new();
    let mut probs = Vec::new();
    for (&abs, &prob) in mods.abs.iter().zip(&mods.probs) {
        if abs < start || abs >= end {
            continue;
        }
        let widx = window_cs
            .iter()
            .position(|&p| p == abs)
            .expect("surviving mod must land on a C in the output window")
            as isize;
        deltas.push((widx - prev - 1) as usize);
        probs.push(prob);
        prev = widx;
    }
    if deltas.is_empty() {
        return None;
    }
    let mut mm = String::from("C+m");
    for d in deltas {
        mm.push_str(&format!(",{d}"));
    }
    mm.push(';');
    Some((mm, probs, (end - start) as i32))
}

fn expected_fastq(reads: &[SourceRead], cfg: ExpectCfg) -> Vec<FastqRecord> {
    let mut out = Vec::new();
    for read in reads {
        let produced = produced_windows(read, cfg);
        let total = produced.len();
        for (idx, (s, e)) in produced.into_iter().enumerate() {
            if passes(&read.seq[s..e], &read.qual[s..e], cfg) {
                out.push(FastqRecord {
                    head: output_name(&read.id, total, idx),
                    seq: read.seq[s..e].to_vec(),
                    qual: read.qual[s..e].to_vec(),
                });
            }
        }
    }
    out.sort();
    out
}

fn expected_fastq_from_bam_with_tags(reads: &[SourceRead], cfg: ExpectCfg) -> Vec<FastqRecord> {
    let mut out = Vec::new();
    for read in reads {
        let produced = produced_windows(read, cfg);
        let total = produced.len();
        for (idx, (s, e)) in produced.into_iter().enumerate() {
            if !passes(&read.seq[s..e], &read.qual[s..e], cfg) {
                continue;
            }
            let mut head = output_name(&read.id, total, idx);
            if let Some(rg) = &read.rg {
                head.push_str(&format!("\tRG:Z:{rg}"));
            }
            if let Some((mm, ml, mn)) = expected_mods(read, s, e) {
                head.push_str(&format!("\tMM:Z:{mm}\tML:B:C"));
                for p in ml {
                    head.push_str(&format!(",{p}"));
                }
                head.push_str(&format!("\tMN:i:{mn}"));
            }
            out.push(FastqRecord {
                head,
                seq: read.seq[s..e].to_vec(),
                qual: read.qual[s..e].to_vec(),
            });
        }
    }
    out.sort();
    out
}

fn expected_bam(reads: &[SourceRead], cfg: ExpectCfg, update_moves: bool) -> Vec<BamRecord> {
    let mut out = Vec::new();
    for read in reads {
        let produced = produced_windows(read, cfg);
        let total = produced.len();
        for (idx, (s, e)) in produced.into_iter().enumerate() {
            if !passes(&read.seq[s..e], &read.qual[s..e], cfg) {
                continue;
            }
            let mods = expected_mods(read, s, e);
            let (mm, ml, mn) = match mods {
                Some((mm, ml, mn)) => (Some(mm), Some(ml), Some(mn)),
                None => (None, None, None),
            };
            let trimmed = s != 0 || e != read.seq.len();
            let (mv, ts, ns) = if read.moves && !trimmed {
                (
                    Some({
                        let mut mv = Vec::with_capacity(read.seq.len() + 1);
                        mv.push(2);
                        mv.extend(std::iter::repeat_n(1, read.seq.len()));
                        mv
                    }),
                    Some(i64::from(read.ts)),
                    Some(i64::from(read.ts) + (read.seq.len() as i64 * 2)),
                )
            } else if update_moves && read.moves {
                let kept_len = e - s;
                (
                    Some({
                        let mut mv = Vec::with_capacity(kept_len + 1);
                        mv.push(2);
                        mv.extend(std::iter::repeat_n(1, kept_len));
                        mv
                    }),
                    Some(i64::from(read.ts) + (s as i64 * 2)),
                    Some(i64::from(read.ts) + (e as i64 * 2)),
                )
            } else {
                (None, None, None)
            };
            out.push(BamRecord {
                name: output_name(&read.id, total, idx),
                seq: read.seq[s..e].to_vec(),
                qual: read.qual[s..e].to_vec(),
                rg: read.rg.clone(),
                mm,
                ml,
                mn,
                ip: read.ip.as_ref().map(|ip| ip[s..e].to_vec()),
                mv,
                ts,
                ns,
            });
        }
    }
    out.sort();
    out
}

fn read_fastq(path: &Path) -> Vec<FastqRecord> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let mut out = Vec::new();
    while let Some(head) = lines.next() {
        let seq = lines
            .next()
            .expect("FASTQ sequence line")
            .as_bytes()
            .to_vec();
        let plus = lines.next().expect("FASTQ plus line");
        let qual = lines.next().expect("FASTQ quality line");
        assert_eq!(plus, "+");
        out.push(FastqRecord {
            head: head.strip_prefix('@').unwrap().to_string(),
            seq,
            qual: qual.as_bytes().iter().map(|q| q - 33).collect(),
        });
    }
    out.sort();
    out
}

fn aux_string(rec: &RecordBuf, tag: Tag) -> Option<String> {
    match rec.data().get(&tag) {
        Some(Value::String(s)) => Some(String::from_utf8(s.to_vec()).unwrap()),
        _ => None,
    }
}

fn aux_i64(rec: &RecordBuf, tag: Tag) -> Option<i64> {
    match rec.data().get(&tag) {
        Some(Value::Int8(n)) => Some(i64::from(*n)),
        Some(Value::UInt8(n)) => Some(i64::from(*n)),
        Some(Value::Int16(n)) => Some(i64::from(*n)),
        Some(Value::UInt16(n)) => Some(i64::from(*n)),
        Some(Value::Int32(n)) => Some(i64::from(*n)),
        Some(Value::UInt32(n)) => Some(i64::from(*n)),
        _ => None,
    }
}

fn read_bam(path: &Path) -> Vec<BamRecord> {
    let mut reader = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let mut buf = RecordBuf::default();
    let mut out = Vec::new();
    while reader.read_record_buf(&header, &mut buf).unwrap() != 0 {
        let mm = aux_string(&buf, Tag::BASE_MODIFICATIONS);
        let ml = match buf.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) => Some(v.clone()),
            _ => None,
        };
        let mn = match buf.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(n)) => Some(*n),
            _ => None,
        };
        let ip = match buf.data().get(&Tag::new(b'i', b'p')) {
            Some(Value::Array(Array::UInt8(v))) => Some(v.clone()),
            _ => None,
        };
        let mv = match buf.data().get(&Tag::new(b'm', b'v')) {
            Some(Value::Array(Array::Int8(v))) => Some(v.clone()),
            _ => None,
        };
        out.push(BamRecord {
            name: String::from_utf8(buf.name().unwrap().to_vec()).unwrap(),
            seq: buf.sequence().as_ref().to_vec(),
            qual: buf.quality_scores().as_ref().to_vec(),
            rg: aux_string(&buf, Tag::new(b'R', b'G')),
            mm,
            ml,
            mn,
            ip,
            mv,
            ts: aux_i64(&buf, Tag::new(b't', b's')),
            ns: aux_i64(&buf, Tag::new(b'n', b's')),
        });
    }
    out.sort();
    out
}

fn assert_same<T>(label: &str, actual: &[T], expected: &[T])
where
    T: std::fmt::Debug + PartialEq,
{
    if actual == expected {
        return;
    }
    let first = actual
        .iter()
        .zip(expected)
        .position(|(a, e)| a != e)
        .unwrap_or_else(|| actual.len().min(expected.len()));
    panic!(
        "{label} mismatch: actual {} records, expected {} records, first mismatch at {first}\nactual: {:?}\nexpected: {:?}",
        actual.len(),
        expected.len(),
        actual.get(first),
        expected.get(first)
    );
}

fn fixed_crop_filter_cfg() -> ExpectCfg {
    ExpectCfg {
        head: 5,
        tail: 7,
        min_len: 64,
        max_len: 110,
        min_qual: 20.0,
        min_gc: Some(0.25),
        max_gc: Some(0.75),
        quality: QualityOp::None,
        adapters: false,
    }
}

fn qual_trim_cfg() -> ExpectCfg {
    ExpectCfg {
        head: 0,
        tail: 0,
        min_len: 50,
        max_len: usize::MAX,
        min_qual: 0.0,
        min_gc: None,
        max_gc: None,
        quality: QualityOp::Trim(20),
        adapters: false,
    }
}

fn qual_split_cfg() -> ExpectCfg {
    ExpectCfg {
        head: 0,
        tail: 0,
        min_len: 20,
        max_len: usize::MAX,
        min_qual: 0.0,
        min_gc: None,
        max_gc: None,
        quality: QualityOp::Split {
            cutoff: 20,
            window: 3,
        },
        adapters: false,
    }
}

fn update_moves_cfg() -> ExpectCfg {
    ExpectCfg {
        head: 5,
        tail: 7,
        min_len: 1,
        max_len: usize::MAX,
        min_qual: 0.0,
        min_gc: None,
        max_gc: None,
        quality: QualityOp::None,
        adapters: false,
    }
}

fn adapter_cfg() -> ExpectCfg {
    ExpectCfg {
        head: 0,
        tail: 0,
        min_len: 20,
        max_len: usize::MAX,
        min_qual: 0.0,
        min_gc: None,
        max_gc: None,
        quality: QualityOp::None,
        adapters: true,
    }
}

fn write_inputs(dir: &Path) -> (Vec<SourceRead>, std::path::PathBuf, std::path::PathBuf) {
    let reads = corpus();
    let fastq_gz = dir.join("reads.fastq.gz");
    let bam = dir.join("reads.bam");
    write_fastq_gz(&fastq_gz, &reads);
    write_bam(&bam, &reads);
    (reads, fastq_gz, bam)
}

fn write_adapter_fasta(dir: &Path) -> std::path::PathBuf {
    let path = dir.join("adapter.fa");
    std::fs::write(
        &path,
        format!(">planted\n{}\n", std::str::from_utf8(ADAPTER).unwrap()),
    )
    .unwrap();
    path
}

#[test]
fn fastq_gz_corpus_fixed_crop_filter_matches_expected() {
    let dir = tempfile::tempdir().unwrap();
    let (reads, fastq_gz, _) = write_inputs(dir.path());
    let out = dir.path().join("out.fastq");

    whittle()
        .arg("-i")
        .arg(&fastq_gz)
        .arg("-o")
        .arg(&out)
        .args([
            "-H", "5", "-T", "7", "-l", "64", "-L", "110", "-g", "0.25", "-G", "0.75", "-q", "20",
            "-t", "4", "--quiet",
        ])
        .assert()
        .success();

    assert_same(
        "FASTQ.gz fixed crop/filter",
        &read_fastq(&out),
        &expected_fastq(&reads, fixed_crop_filter_cfg()),
    );
}

#[test]
fn fastq_gz_corpus_quality_trim_matches_expected() {
    let dir = tempfile::tempdir().unwrap();
    let (reads, fastq_gz, _) = write_inputs(dir.path());
    let out = dir.path().join("out.fastq");

    whittle()
        .arg("-i")
        .arg(&fastq_gz)
        .arg("-o")
        .arg(&out)
        .args(["--qual-trim", "20", "-l", "50", "-t", "4", "--quiet"])
        .assert()
        .success();

    assert_same(
        "FASTQ.gz quality trim",
        &read_fastq(&out),
        &expected_fastq(&reads, qual_trim_cfg()),
    );
}

#[test]
fn fastq_gz_corpus_quality_split_matches_expected() {
    let dir = tempfile::tempdir().unwrap();
    let (reads, fastq_gz, _) = write_inputs(dir.path());
    let out = dir.path().join("out.fastq");

    whittle()
        .arg("-i")
        .arg(&fastq_gz)
        .arg("-o")
        .arg(&out)
        .args([
            "--qual-split",
            "20",
            "--qual-split-window",
            "3",
            "-l",
            "20",
            "-t",
            "4",
            "--quiet",
        ])
        .assert()
        .success();

    assert_same(
        "FASTQ.gz quality split",
        &read_fastq(&out),
        &expected_fastq(&reads, qual_split_cfg()),
    );
}

#[test]
fn fastq_gz_corpus_adapter_trim_and_split_matches_expected() {
    let dir = tempfile::tempdir().unwrap();
    let (reads, fastq_gz, _) = write_inputs(dir.path());
    let adapters = write_adapter_fasta(dir.path());
    let out = dir.path().join("out.fastq");

    whittle()
        .arg("-i")
        .arg(&fastq_gz)
        .arg("-o")
        .arg(&out)
        .arg("--adapter-fasta")
        .arg(&adapters)
        .args([
            "--adapter-error-rate",
            "0",
            "--adapter-end-size",
            "20",
            "-l",
            "20",
            "-t",
            "4",
            "--quiet",
        ])
        .assert()
        .success();

    let actual = read_fastq(&out);
    let expected = expected_fastq(&reads, adapter_cfg());
    assert!(
        actual.iter().any(|r| r.head.ends_with("_segment_2")),
        "adapter fixture must exercise interior split suffixes"
    );
    assert_same("FASTQ.gz adapter trim/split", &actual, &expected);
}

#[test]
fn bam_corpus_fixed_crop_filter_matches_expected() {
    let dir = tempfile::tempdir().unwrap();
    let (reads, _, bam) = write_inputs(dir.path());
    let out = dir.path().join("out.bam");

    whittle()
        .arg("-i")
        .arg(&bam)
        .arg("-o")
        .arg(&out)
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "-H",
            "5",
            "-T",
            "7",
            "-l",
            "64",
            "-L",
            "110",
            "-g",
            "0.25",
            "-G",
            "0.75",
            "-q",
            "20",
            "-t",
            "4",
            "--quiet",
        ])
        .assert()
        .success();

    assert_same(
        "BAM fixed crop/filter",
        &read_bam(&out),
        &expected_bam(&reads, fixed_crop_filter_cfg(), false),
    );
}

#[test]
fn bam_corpus_adapter_trim_and_split_matches_expected() {
    let dir = tempfile::tempdir().unwrap();
    let (reads, _, bam) = write_inputs(dir.path());
    let adapters = write_adapter_fasta(dir.path());
    let out = dir.path().join("out.bam");

    whittle()
        .arg("-i")
        .arg(&bam)
        .arg("-o")
        .arg(&out)
        .arg("--adapter-fasta")
        .arg(&adapters)
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "--adapter-error-rate",
            "0",
            "--adapter-end-size",
            "20",
            "-l",
            "20",
            "-t",
            "4",
            "--quiet",
        ])
        .assert()
        .success();

    let actual = read_bam(&out);
    let expected = expected_bam(&reads, adapter_cfg(), false);
    assert!(
        actual.iter().any(|r| r.name.ends_with("_segment_2")),
        "adapter fixture must exercise interior split suffixes"
    );
    assert_same("BAM adapter trim/split", &actual, &expected);
}

#[test]
fn bam_corpus_to_fastq_tags_match_expected() {
    let dir = tempfile::tempdir().unwrap();
    let (reads, _, bam) = write_inputs(dir.path());
    let out = dir.path().join("out.fastq");

    whittle()
        .arg("-i")
        .arg(&bam)
        .arg("-o")
        .arg(&out)
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "fastq",
            "-H",
            "5",
            "-T",
            "7",
            "--fastq-tags",
            "MM,ML,RG",
            "-t",
            "4",
            "--quiet",
        ])
        .assert()
        .success();

    assert_same(
        "BAM to FASTQ selected tags",
        &read_fastq(&out),
        &expected_fastq_from_bam_with_tags(&reads, update_moves_cfg()),
    );
}

#[test]
fn bam_corpus_update_moves_matches_expected() {
    let dir = tempfile::tempdir().unwrap();
    let (reads, _, bam) = write_inputs(dir.path());
    let out = dir.path().join("out.bam");

    whittle()
        .arg("-i")
        .arg(&bam)
        .arg("-o")
        .arg(&out)
        .args([
            "--in-format",
            "bam",
            "--out-format",
            "bam",
            "--update-moves",
            "-H",
            "5",
            "-T",
            "7",
            "-t",
            "4",
            "--quiet",
        ])
        .assert()
        .success();

    let actual_by_name: BTreeMap<_, _> = read_bam(&out)
        .into_iter()
        .map(|r| (r.name.clone(), r))
        .collect();
    let expected_by_name: BTreeMap<_, _> = expected_bam(&reads, update_moves_cfg(), true)
        .into_iter()
        .map(|r| (r.name.clone(), r))
        .collect();
    assert_eq!(actual_by_name.len(), expected_by_name.len());

    for (name, expected) in expected_by_name {
        let actual = actual_by_name
            .get(&name)
            .unwrap_or_else(|| panic!("missing update-moves output record {name}"));
        assert_eq!(
            (&actual.seq, &actual.qual, &actual.mv, actual.ts, actual.ns),
            (
                &expected.seq,
                &expected.qual,
                &expected.mv,
                expected.ts,
                expected.ns
            ),
            "signal tags mismatch for {name}"
        );
    }
}
