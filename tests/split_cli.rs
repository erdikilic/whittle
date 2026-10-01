//! `--split-by`: primer split calls in tag-only mode, over the compiled binary.
//! Every output record carries a `wt:Z` tag naming its target, `unassigned`
//! or `ambiguous`.

use std::io::Write;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::{self as sam};
use predicates::prelude::*;

/// Generates deterministic SplitMix64 bases.
fn splitmix_dna(seed: u64, len: usize) -> Vec<u8> {
    let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(seed);
    (0..len)
        .map(|_| {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            b"ACGT"[((z >> 62) & 0b11) as usize]
        })
        .collect()
}

/// Reverse complement of an ACGT sequence.
fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// Concatenates sequence parts.
fn cat(parts: &[&[u8]]) -> Vec<u8> {
    parts.concat()
}

/// The primers of the two-target sheet: 16S and ITS forward and reverse,
/// each written 5' to 3'.
struct Primers {
    f16: Vec<u8>,
    r16: Vec<u8>,
    fits: Vec<u8>,
    rits: Vec<u8>,
}

fn primers() -> Primers {
    Primers {
        f16: splitmix_dna(901, 20),
        r16: splitmix_dna(902, 22),
        fits: splitmix_dna(903, 21),
        rits: splitmix_dna(904, 20),
    }
}

/// Writes the two-target TSV sheet to `dir/sheet.tsv`.
fn write_sheet(dir: &Path, p: &Primers) -> PathBuf {
    let path = dir.join("sheet.tsv");
    let s = |v: &[u8]| String::from_utf8(v.to_vec()).unwrap();
    std::fs::write(
        &path,
        format!(
            "target\tfwd\trev\n16S\t{}\t{}\nITS\t{}\t{}\n",
            s(&p.f16),
            s(&p.r16),
            s(&p.fits),
            s(&p.rits)
        ),
    )
    .unwrap();
    path
}

/// The insert of read `i`.
fn insert(i: u64) -> Vec<u8> {
    splitmix_dna(1000 + i, 400)
}

/// Three reads: `plus16` (16S on the plus strand), `minusits` (ITS on the
/// minus strand) and `none` (no primer).
fn three_reads(p: &Primers) -> Vec<(String, Vec<u8>)> {
    vec![
        ("plus16".into(), cat(&[&p.f16, &insert(1), &rc(&p.r16)])),
        ("minusits".into(), cat(&[&p.rits, &insert(2), &rc(&p.fits)])),
        ("none".into(), insert(3)),
    ]
}

fn write_fastq(path: &Path, reads: &[(String, Vec<u8>)]) {
    let mut f = std::fs::File::create(path).unwrap();
    for (name, seq) in reads {
        writeln!(
            f,
            "@{name}\n{}\n+\n{}",
            String::from_utf8_lossy(seq),
            "I".repeat(seq.len())
        )
        .unwrap();
    }
}

/// Writes `reads` as unmapped records, each with an `xx:Z:keep` tag; the
/// read named `none` also carries a stale `wt:Z:old` tag.
fn write_ubam(path: &Path, reads: &[(String, Vec<u8>)]) {
    write_ubam_with_header(path, reads, sam::Header::default());
}

/// Writes `reads` as `write_ubam` does, under `header`.
fn write_ubam_with_header(path: &Path, reads: &[(String, Vec<u8>)], header: sam::Header) {
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();
    for (name, seq) in reads {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::UNMAPPED;
        *r.name_mut() = Some(name.as_bytes().into());
        *r.sequence_mut() = seq.clone().into();
        *r.quality_scores_mut() = vec![40; seq.len()].into();
        let data = r.data_mut();
        data.insert(Tag::new(b'x', b'x'), Value::String("keep".into()));
        if name == "none" {
            data.insert(Tag::new(b'w', b't'), Value::String("old".into()));
        }
        w.write_alignment_record(&header, &r).unwrap();
    }
    w.try_finish().unwrap();
}

fn whittle() -> Command {
    let mut cmd = Command::cargo_bin("whittle").unwrap();
    cmd.env_remove("WHITTLE_LOG");
    cmd
}

/// Every FASTQ record of `text` as (header, sequence), in file order.
fn fastq_records(text: &str) -> Vec<(String, String)> {
    let lines: Vec<&str> = text.lines().collect();
    lines
        .chunks(4)
        .map(|r| (r[0].trim_start_matches('@').to_string(), r[1].to_string()))
        .collect()
}

/// The record of `records` whose read id is `name`.
fn by_name<'a>(records: &'a [(String, String)], name: &str) -> &'a (String, String) {
    records
        .iter()
        .find(|(h, _)| h.split(['\t', ' ']).next() == Some(name))
        .unwrap_or_else(|| panic!("record {name} missing: {records:?}"))
}

/// One uBAM record: name, sequence, and its string aux tags as `(tag, value)`
/// in record order.
type BamRecord = (String, String, Vec<(String, String)>);

/// Every record of a uBAM, in file order.
fn bam_records(path: &Path) -> Vec<BamRecord> {
    let mut rdr = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let hdr = rdr.read_header().unwrap();
    let mut buf = RecordBuf::default();
    let mut out = Vec::new();
    while rdr.read_record_buf(&hdr, &mut buf).unwrap() > 0 {
        let tags = buf
            .data()
            .iter()
            .filter_map(|(tag, value)| match value {
                Value::String(s) => Some((
                    String::from_utf8_lossy(tag.as_ref()).into_owned(),
                    String::from_utf8_lossy(s.as_ref()).into_owned(),
                )),
                _ => None,
            })
            .collect();
        out.push((
            String::from_utf8_lossy(buf.name().unwrap().as_ref()).into_owned(),
            String::from_utf8_lossy(buf.sequence().as_ref()).into_owned(),
            tags,
        ));
    }
    out
}

#[test]
fn tag_only_fastq_labels_every_read() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    for threads in ["1", "2"] {
        let out = dir.path().join(format!("out{threads}.fastq"));
        whittle()
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .assert()
            .success();
        let text = std::fs::read_to_string(&out).unwrap();
        let recs = fastq_records(&text);
        assert_eq!(recs.len(), 3, "{text}");
        let (head, seq) = by_name(&recs, "plus16");
        assert_eq!(head, "plus16\twt:Z:16S");
        assert_eq!(seq.as_bytes(), insert(1), "Both 16S primers are trimmed");
        let (head, seq) = by_name(&recs, "minusits");
        assert_eq!(head, "minusits\twt:Z:ITS");
        assert_eq!(seq.as_bytes(), insert(2), "Both ITS primers are trimmed");
        let (head, seq) = by_name(&recs, "none");
        assert_eq!(head, "none\twt:Z:unassigned");
        assert_eq!(seq.as_bytes(), insert(3));
    }
}

/// A tagged FASTQ header keeps its tags; the target tag follows them.
#[test]
fn tagged_fastq_appends_the_target_after_the_input_tags() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.fastq");
    let reads: Vec<(String, Vec<u8>)> = three_reads(&p)
        .into_iter()
        .map(|(name, seq)| (format!("{name}\txx:Z:keep"), seq))
        .collect();
    write_fastq(&input, &reads);
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    assert_eq!(by_name(&recs, "plus16").0, "plus16\txx:Z:keep\twt:Z:16S");
    assert_eq!(by_name(&recs, "none").0, "none\txx:Z:keep\twt:Z:unassigned");
}

#[test]
fn tag_only_bam_writes_wt_tag() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.bam");
    write_ubam(&input, &three_reads(&p));
    for threads in ["1", "2"] {
        let out = dir.path().join(format!("out{threads}.bam"));
        whittle()
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .args(["--preserve-order"])
            .assert()
            .success();
        let recs = bam_records(&out);
        let tags = |name: &str| {
            recs.iter()
                .find(|(n, _, _)| n == name)
                .unwrap_or_else(|| panic!("{name} missing: {recs:?}"))
                .2
                .clone()
        };
        let pair = |t: &str, v: &str| (t.to_string(), v.to_string());
        assert_eq!(tags("plus16"), [pair("xx", "keep"), pair("wt", "16S")]);
        assert_eq!(tags("minusits"), [pair("xx", "keep"), pair("wt", "ITS")]);
        // The input `wt` is replaced in place, not duplicated.
        assert_eq!(tags("none"), [pair("xx", "keep"), pair("wt", "unassigned")]);
        let seq = |name: &str| recs.iter().find(|(n, _, _)| n == name).unwrap().1.clone();
        assert_eq!(seq("plus16").as_bytes(), insert(1));
        assert_eq!(seq("none").as_bytes(), insert(3));
    }
}

/// BAM to FASTQ writes the target tag after the carried tags.
#[test]
fn bam_to_fastq_appends_the_target_tag() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.bam");
    write_ubam(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    assert!(
        by_name(&recs, "plus16")
            .0
            .ends_with("\txx:Z:keep\twt:Z:16S"),
        "{recs:?}"
    );
    assert_eq!(by_name(&recs, "none").0, "none\txx:Z:keep\twt:Z:unassigned");
}

#[test]
fn retain_keeps_primer_bases() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .args(["--split-action", "retain"])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    let (head, seq) = by_name(&recs, "plus16");
    assert_eq!(head, "plus16\twt:Z:16S");
    assert!(
        seq.as_bytes().starts_with(&p.f16),
        "The forward primer is retained: {seq}"
    );
    assert!(seq.as_bytes().ends_with(&rc(&p.r16)));
}

/// Primer trimming and the target call work from the full read after a
/// fixed crop moves the piece start.
#[test]
fn head_crop_keeps_the_call() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1", "-H", "5"])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    let (head, seq) = by_name(&recs, "plus16");
    assert_eq!(head, "plus16\twt:Z:16S");
    assert_eq!(seq.as_bytes(), &insert(1)[5..]);
}

#[test]
fn split_primers_survive_presence_sampling() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let lsk_front = b"CCTGTACTTCGTTCAGTTACGTATTGC";
    let reads: Vec<(String, Vec<u8>)> = (0..150)
        .map(|i| {
            let seq = cat(&[
                lsk_front,
                &splitmix_dna(2000 + i, 5),
                &p.f16,
                &insert(10 + i),
                &rc(&p.r16),
            ]);
            (format!("r{i}"), seq)
        })
        .collect();
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "2"])
        .args(["--adapter-preset", "lsk114"])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    assert_eq!(recs.len(), 150);
    for (head, seq) in &recs {
        assert!(head.ends_with("\twt:Z:16S"), "{head}");
        assert!(
            !seq.contains(std::str::from_utf8(lsk_front).unwrap()),
            "The adapter is trimmed: {head}"
        );
    }
}

#[test]
fn split_flags_require_split() {
    for args in [
        &["--split-lead", "3"][..],
        &["--split-require", "both"],
        &["--split-action", "retain"],
        &["--split-discard", "unassigned"],
    ] {
        whittle()
            .args(args)
            .args(["--input-format", "fastq"])
            .write_stdin("")
            .assert()
            .failure()
            .stderr(predicates::str::contains("--split-by <SPEC>"));
    }
}

#[test]
fn split_mab114_enables_preset() {
    let dir = tempfile::tempdir().unwrap();
    // The MAB114 construct: the rapid barcode flank and a barcode ahead of
    // the kit primers, the degenerate ones resolved to concrete bases here.
    let rbk_front = b"GCTTGGGTGTTTAACC";
    let tp01 = b"GCACCTGGAACTTGTGCCTTCCAC";
    let s16_f = b"AGAGTTTGATCCTGGCTCAG";
    let s16_r = b"CGGTTACCTTGTTACGACTT";
    let its1 = b"TCCGTAGGTGAACCTGCGG";
    let its4 = b"TCCTCCGCTTATTGATATGC";
    let reads = vec![
        (
            "s16".to_string(),
            cat(&[rbk_front, tp01, s16_f, &insert(1), &rc(s16_r)]),
        ),
        (
            "sits".to_string(),
            cat(&[rbk_front, tp01, its1, &insert(2), &rc(its4)]),
        ),
    ];
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", "mab114", "-t", "1"])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    let (head, seq) = by_name(&recs, "s16");
    assert_eq!(head, "s16\twt:Z:16S");
    assert_eq!(
        seq.as_bytes(),
        insert(1),
        "Barcode flank and primers trimmed"
    );
    let (head, seq) = by_name(&recs, "sits");
    assert_eq!(head, "sits\twt:Z:ITS");
    assert_eq!(seq.as_bytes(), insert(2));
}

/// Reads carrying different primer variants of the MAB114 16S and ITS
/// mixes, on both strands, plus a primer-free read: `(name, sequence,
/// target)`.
fn mab114_mix_pool() -> Vec<(String, Vec<u8>, &'static str)> {
    let bor_f = b"AGAGTTTGATCCTGGCTTAG";
    let ent_f = b"AGAGTTTGATCATGGCTCAG";
    let chl_f = b"AGAATTTGATCTTAGTTCAG";
    let mix_r = b"CGGTTACCTTGTTACGACTT";
    let bor_r = b"CGGCTACCTTGTTACGACTT";
    let chl_r = b"GGGCTACCTTGTTACGACTT";
    let its1 = b"TCCGTAGGTGAACCTGCGG";
    let its1_fus = b"TCCGTTGGTGAACCAGCGG";
    let its1_mal = b"TCTGTAGGTGAACCTGCAG";
    let its4 = b"TCCTCCGCTTATTGATATGC";
    let its4_pyt = b"TCCTCCGCTTATTAATATGC";
    let plus = |f: &[u8], i: u64, r: &[u8]| cat(&[f, &insert(i), &rc(r)]);
    let minus = |f: &[u8], i: u64, r: &[u8]| cat(&[r, &insert(i), &rc(f)]);
    vec![
        ("bor_chl".to_string(), plus(bor_f, 1, chl_r), "16S"),
        ("ent_mix".to_string(), minus(ent_f, 2, mix_r), "16S"),
        ("chl_bor".to_string(), plus(chl_f, 3, bor_r), "16S"),
        ("mal_pyt".to_string(), plus(its1_mal, 4, its4_pyt), "ITS"),
        ("fus_its4".to_string(), minus(its1_fus, 5, its4), "ITS"),
        ("its1_pyt".to_string(), plus(its1, 6, its4_pyt), "ITS"),
        ("none".to_string(), insert(7), "unassigned"),
    ]
}

/// Runs `specs` over `mab114_mix_pool` under `--split-require both` with a
/// `{target}` template, and checks that each target's file holds exactly
/// its planted reads, whichever variants they carry.
fn assert_mab114_mix_split(specs: &[&str]) {
    let dir = tempfile::tempdir().unwrap();
    let pool = mab114_mix_pool();
    let reads: Vec<(String, Vec<u8>)> = pool
        .iter()
        .map(|(n, s, _)| (n.clone(), s.clone()))
        .collect();
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    let out = dir.path().join("out");
    let mut cmd = whittle();
    cmd.args(["-i", input.to_str().unwrap()])
        .args(["-o", out.join("{target}.fastq").to_str().unwrap()])
        .args(["--split-require", "both", "--preserve-order", "-t", "1"]);
    for spec in specs {
        cmd.args(["--split-by", spec]);
    }
    cmd.assert().success();
    assert_eq!(
        files_under(&out),
        ["16S.fastq", "ITS.fastq", "unassigned.fastq"]
    );
    for bin in ["16S", "ITS", "unassigned"] {
        let text = std::fs::read_to_string(out.join(format!("{bin}.fastq"))).unwrap();
        let expected: Vec<&str> = pool
            .iter()
            .filter(|(_, _, target)| *target == bin)
            .map(|(name, _, _)| name.as_str())
            .collect();
        assert_eq!(ids_with_label(&text, bin), expected, "{bin}");
    }
}

/// Inline targets holding the MAB114 primer mixes as lists assign reads
/// that carry different variants to one file per target.
#[test]
fn inline_primer_mixes_assign_every_variant_to_its_target() {
    assert_mab114_mix_split(&[
        "16S:F:AGRGTTYGATYMTGGCTCAG,AGAGTTTGATCCTGGCTTAG,AGAATTTGATCTTRGTTCAG,\
         AGAGTTTGATCATGGCTCAG:R:SGGYTACCTTGTTACGACTT,CGGCTACCTTGTTACGACTT,GGGCTACCTTGTTACGACTT",
        "ITS:F:TCCGTAGGTGAACCTGCGG,TCCGTTGGTGAACCAGCGG,TCTGTAGGTGAACCTGCAG\
         :R:TCCTCCGCTTATTGATATGC,TCCTCCGCTTATTAATATGC",
    ]);
}

/// The `mab114` preset holds the same mixes as the inline lists.
#[test]
fn mab114_preset_assigns_every_variant_to_its_target() {
    assert_mab114_mix_split(&["mab114"]);
}

/// A TSV sheet with list cells and a FASTA sheet with repeated ends give
/// the same targets as the inline lists.
#[test]
fn sheet_files_with_primer_lists_assign_every_variant() {
    let dir = tempfile::tempdir().unwrap();
    let tsv = dir.path().join("mix.tsv");
    std::fs::write(
        &tsv,
        "target\tfwd\trev\n\
         16S\tAGRGTTYGATYMTGGCTCAG,AGAGTTTGATCCTGGCTTAG,AGAATTTGATCTTRGTTCAG,\
         AGAGTTTGATCATGGCTCAG\tSGGYTACCTTGTTACGACTT,CGGCTACCTTGTTACGACTT,GGGCTACCTTGTTACGACTT\n",
    )
    .unwrap();
    let fasta = dir.path().join("mix.fa");
    std::fs::write(
        &fasta,
        ">ITS1 target=ITS end=fwd\nTCCGTAGGTGAACCTGCGG\n\
         >ITS1_Fus target=ITS end=fwd\nTCCGTTGGTGAACCAGCGG\n\
         >ITS1_Mal target=ITS end=fwd\nTCTGTAGGTGAACCTGCAG\n\
         >ITS4 target=ITS end=rev\nTCCTCCGCTTATTGATATGC\n\
         >ITS4_Pyt target=ITS end=rev\nTCCTCCGCTTATTAATATGC\n",
    )
    .unwrap();
    assert_mab114_mix_split(&[tsv.to_str().unwrap(), fasta.to_str().unwrap()]);
}

/// An inline-only run, with no sheet file or preset, assigns the planted 16S
/// and ITS reads from two `--split-by` targets.
#[test]
fn inline_split_by_assigns_planted_reads() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let s = |v: &[u8]| String::from_utf8(v.to_vec()).unwrap();
    let spec16 = format!("16S:F:{}:R:{}", s(&p.f16), s(&p.r16));
    let spec_its = format!("ITS:F:{}:R:{}", s(&p.fits), s(&p.rits));
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", &spec16])
        .args(["--split-by", &spec_its])
        .args(["-t", "1"])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    assert_eq!(recs.len(), 3, "{text}");
    let (head, seq) = by_name(&recs, "plus16");
    assert_eq!(head, "plus16\twt:Z:16S");
    assert_eq!(seq.as_bytes(), insert(1));
    let (head, seq) = by_name(&recs, "minusits");
    assert_eq!(head, "minusits\twt:Z:ITS");
    assert_eq!(seq.as_bytes(), insert(2));
    let (head, _) = by_name(&recs, "none");
    assert_eq!(head, "none\twt:Z:unassigned");
}

/// `--split-by` given twice, both times an inline target, merges into one
/// sheet: a plain target and a pool, routed through a `{group}` template.
/// The pool's two targets share a group and collapse into one file.
#[test]
fn split_by_given_twice_as_inline_targets_pool_collapses_to_one_file() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let s = |v: &[u8]| String::from_utf8(v.to_vec()).unwrap();
    let spec16 = format!("16S:F:{}:R:{}", s(&p.f16), s(&p.r16));
    let (f34a, r34a) = (splitmix_dna(910, 20), splitmix_dna(911, 20));
    let (f34b, r34b) = (splitmix_dna(912, 20), splitmix_dna(913, 20));
    let spec_v34 = format!(
        "V34:F:{}:R:{},F:{}:R:{}",
        s(&f34a),
        s(&r34a),
        s(&f34b),
        s(&r34b)
    );
    let reads = vec![
        ("s16".to_string(), cat(&[&p.f16, &insert(1), &rc(&p.r16)])),
        ("v34a".to_string(), cat(&[&f34a, &insert(2), &rc(&r34a)])),
        ("v34b".to_string(), cat(&[&f34b, &insert(3), &rc(&r34b)])),
    ];
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    let out = dir.path().join("out");
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", out.join("{group}.fastq").to_str().unwrap()])
        .args(["--split-by", &spec16])
        .args(["--split-by", &spec_v34])
        .args(["-t", "1"])
        .assert()
        .success();
    assert_eq!(files_under(&out), ["16S.fastq", "V34.fastq"]);
    let text = std::fs::read_to_string(out.join("16S.fastq")).unwrap();
    assert_eq!(ids_with_label(&text, "16S"), ["s16"]);
    let text = std::fs::read_to_string(out.join("V34.fastq")).unwrap();
    assert_eq!(
        ids_with_label(&text, "V34"),
        ["v34a", "v34b"],
        "the pool's two targets collapse into the one V34 group file"
    );
}

/// Three barcode-style targets, each with one sequence in both role lists,
/// split reads that carry the barcode at both ends in either orientation or
/// at one end only. Both located ends are trimmed. `--split-require both`
/// assigns only the reads with the barcode at both ends.
#[test]
fn symmetric_targets_split_by_a_barcode_at_both_ends() {
    let dir = tempfile::tempdir().unwrap();
    let s = |v: &[u8]| String::from_utf8(v.to_vec()).unwrap();
    let names = ["BC01", "BC02", "BC03"];
    let barcodes: Vec<Vec<u8>> = (0..3).map(|i| splitmix_dna(920 + i, 24)).collect();
    let mut reads: Vec<(String, Vec<u8>)> = Vec::new();
    for (i, (name, bc)) in names.iter().zip(&barcodes).enumerate() {
        let ins = |j: u64| insert(10 * i as u64 + j);
        let whole = cat(&[bc, &ins(0), &rc(bc)]);
        reads.push((format!("{name}_plus"), whole));
        reads.push((format!("{name}_minus"), rc(&cat(&[bc, &ins(1), &rc(bc)]))));
        reads.push((format!("{name}_five"), cat(&[bc, &ins(2)])));
        reads.push((format!("{name}_three"), cat(&[&ins(3), &rc(bc)])));
    }
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    let run = |out: &Path, require: &str| {
        let mut cmd = whittle();
        cmd.args(["-i", input.to_str().unwrap()])
            .args(["-o", out.join("{target}.fastq").to_str().unwrap()])
            .args(["--split-require", require, "-t", "1", "--preserve-order"]);
        for (name, bc) in names.iter().zip(&barcodes) {
            cmd.args(["--split-by", &format!("{name}:F:{0}:R:{0}", s(bc))]);
        }
        cmd.assert().success();
    };

    let either = dir.path().join("either");
    run(&either, "either");
    assert_eq!(
        files_under(&either),
        ["BC01.fastq", "BC02.fastq", "BC03.fastq"]
    );
    for (i, name) in names.iter().enumerate() {
        let text = std::fs::read_to_string(either.join(format!("{name}.fastq"))).unwrap();
        let ids = ["plus", "minus", "five", "three"].map(|kind| format!("{name}_{kind}"));
        assert_eq!(ids_with_label(&text, name), ids);
        let recs = fastq_records(&text);
        let ins = |j: u64| insert(10 * i as u64 + j);
        let seq = |kind: &str| {
            by_name(&recs, &format!("{name}_{kind}"))
                .1
                .as_bytes()
                .to_vec()
        };
        assert_eq!(seq("plus"), ins(0), "both ends trimmed");
        assert_eq!(seq("minus"), rc(&ins(1)), "both ends trimmed");
        assert_eq!(seq("five"), ins(2));
        assert_eq!(seq("three"), ins(3));
    }

    let both = dir.path().join("both");
    run(&both, "both");
    assert_eq!(
        files_under(&both),
        ["BC01.fastq", "BC02.fastq", "BC03.fastq", "unassigned.fastq"]
    );
    for name in names {
        let text = std::fs::read_to_string(both.join(format!("{name}.fastq"))).unwrap();
        let ids = ["plus", "minus"].map(|kind| format!("{name}_{kind}"));
        assert_eq!(ids_with_label(&text, name), ids);
    }
    let text = std::fs::read_to_string(both.join("unassigned.fastq")).unwrap();
    assert_eq!(ids_with_label(&text, "unassigned").len(), 6);
}

#[test]
fn discard_unassigned_in_tag_only_mode() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .args(["--split-discard", "unassigned"])
        .assert()
        .success()
        .stderr(predicates::str::contains("No trimming or filtering").not());
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    let names: Vec<&str> = recs.iter().map(|(h, _)| h.as_str()).collect();
    assert_eq!(names, ["plus16\twt:Z:16S", "minusits\twt:Z:ITS"]);
}

/// A run without `--split-by` writes no target tag.
#[test]
fn no_split_writes_no_target_tag() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["-t", "1"])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(!text.contains("wt:Z"));
}

/// The end-of-run log carries one line per sheet key, then unassigned by
/// reason, ambiguous and discarded, at the same level as the rest of the
/// summary.
#[test]
fn end_of_run_log_reports_the_split_table() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    let output = whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .success()
        .get_output()
        .stderr
        .clone();
    let stderr = String::from_utf8_lossy(&output);
    assert!(
        stderr.contains("Split 16S: 1 segments, 400 bp (50.0% of assigned)"),
        "{stderr}"
    );
    assert!(
        stderr.contains("Split ITS: 1 segments, 400 bp (50.0% of assigned)"),
        "{stderr}"
    );
    assert!(
        stderr.contains(
            "Split unassigned: 1 segments (no_primer 1, require 0, orientation 0, length 0)"
        ),
        "{stderr}"
    );
    assert!(stderr.contains("Split ambiguous: 0 segments"), "{stderr}");
    assert!(stderr.contains("Split discarded: 0 segments"), "{stderr}");
}

/// Without `--split-by`, the end-of-run log has no split table at all.
#[test]
fn end_of_run_log_has_no_split_table_without_split() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    let output = whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["-t", "1"])
        .assert()
        .success()
        .get_output()
        .stderr
        .clone();
    let stderr = String::from_utf8_lossy(&output);
    assert!(!stderr.contains("Split "), "{stderr}");
}

/// The `tm` field of the `@RG` record of a uBAM header, if any.
fn read_group_trim_mode(path: &Path) -> Option<String> {
    let mut rdr = bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let header = rdr.read_header().unwrap();
    let tm = sam::header::record::value::map::tag::Other::try_from(*b"tm").unwrap();
    header
        .read_groups()
        .values()
        .find_map(|group| group.other_fields().get(&tm))
        .map(|v| String::from_utf8_lossy(v.as_ref()).into_owned())
}

/// Split primers count as trimmed primers in `@RG tm` under `trim`, and not
/// under `retain`, which keeps them.
#[test]
fn retain_leaves_primer_out_of_the_trim_mode() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.bam");
    let header: sam::Header = "@HD\tVN:1.6\n@RG\tID:rg1\n".parse().unwrap();
    write_ubam_with_header(&input, &three_reads(&p), header);
    let mode = |action: &str| {
        let out = dir.path().join(format!("{action}.bam"));
        whittle()
            .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
            .args(["--split-action", action])
            .assert()
            .success();
        read_group_trim_mode(&out)
    };
    let trimmed = mode("trim");
    assert!(
        trimmed.as_deref().is_some_and(|tm| tm.contains("primer")),
        "{trimmed:?}"
    );
    let retained = mode("retain");
    assert!(
        !retained.as_deref().is_some_and(|tm| tm.contains("primer")),
        "{retained:?}"
    );
}

/// A split-only run has no presence detection or sampling, so
/// `--adapter-sample-reads` is reported as having no effect.
#[test]
fn split_only_notes_that_sampling_has_no_effect() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .args(["--adapter-sample-reads", "500"])
        .assert()
        .success()
        .stderr(predicates::str::contains(
            "--adapter-sample-reads has no effect with --split-by alone",
        ));
}

/// Without `--split-by`, the ends-only advisory names the three adapter sources.
#[test]
fn ends_only_advisory_without_split_names_the_adapter_sources() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &three_reads(&p));
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--adapter-ends-only", "-t", "1"])
        .assert()
        .success()
        .stderr(predicates::str::contains(
            "--adapter-ends-only has no effect without --adapter-fasta, --adapter-preset or \
             --adapter-discover",
        ));
}

/// Writes a TSV sheet of `rows` (target, fwd, rev, group) to `dir/name`.
fn write_rows(dir: &Path, name: &str, rows: &[(&str, &[u8], &[u8], &str)]) -> PathBuf {
    let path = dir.join(name);
    let mut text = String::from("target\tfwd\trev\tgroup\n");
    for (target, fwd, rev, group) in rows {
        text.push_str(&format!(
            "{target}\t{}\t{}\t{group}\n",
            String::from_utf8_lossy(fwd),
            String::from_utf8_lossy(rev)
        ));
    }
    std::fs::write(&path, text).unwrap();
    path
}

/// A read with its planted bin.
type Planted = (String, Vec<u8>, &'static str);

/// A pool with planted truth: `n16` 16S reads on alternating strands,
/// `nits` ITS reads and `nnone` primer-free reads, interleaved.
fn pool(p: &Primers, n16: u64, nits: u64, nnone: u64) -> Vec<Planted> {
    let mut left = [n16, nits, nnone];
    let mut out = Vec::new();
    let mut i = 0u64;
    while left.iter().any(|&n| n > 0) {
        let bin = (0..3)
            .map(|k| (i as usize + k) % 3)
            .find(|&k| left[k] > 0)
            .unwrap();
        left[bin] -= 1;
        let body = insert(100 + i);
        let planted = match bin {
            0 if i.is_multiple_of(2) => (cat(&[&p.f16, &body, &rc(&p.r16)]), "16S"),
            0 => (cat(&[&p.r16, &body, &rc(&p.f16)]), "16S"),
            1 => (cat(&[&p.fits, &body, &rc(&p.rits)]), "ITS"),
            _ => (body, "unassigned"),
        };
        out.push((format!("r{i}"), planted.0, planted.1));
        i += 1;
    }
    out
}

/// The (name, sequence) pairs of a planted pool, for the FASTQ and BAM
/// writers.
fn pool_reads(pool: &[Planted]) -> Vec<(String, Vec<u8>)> {
    pool.iter()
        .map(|(n, s, _)| (n.clone(), s.clone()))
        .collect()
}

/// The names of a planted pool's reads in `bin`, in input order.
fn truth(pool: &[Planted], bin: &str) -> Vec<String> {
    pool.iter()
        .filter(|(_, _, b)| *b == bin)
        .map(|(n, _, _)| n.clone())
        .collect()
}

/// Every file under `dir`, as paths relative to it, sorted.
fn files_under(dir: &Path) -> Vec<String> {
    fn walk(root: &Path, dir: &Path, out: &mut Vec<String>) {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                walk(root, &path, out);
            } else {
                out.push(path.strip_prefix(root).unwrap().display().to_string());
            }
        }
    }
    let mut out = Vec::new();
    walk(dir, dir, &mut out);
    out.sort();
    out
}

/// Reads a gzip or BGZF FASTQ file to text.
fn read_gz(path: &Path) -> String {
    use std::io::Read;
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// The read ids of FASTQ `text`, in file order, each checked to carry the
/// target tag `label`.
fn ids_with_label(text: &str, label: &str) -> Vec<String> {
    fastq_records(text)
        .into_iter()
        .map(|(head, _)| {
            assert!(head.ends_with(&format!("\twt:Z:{label}")), "{head}");
            head.split(['\t', ' ']).next().unwrap().to_string()
        })
        .collect()
}

#[test]
fn fastq_gz_split_writes_one_file_per_key() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let reads = pool(&p, 9, 6, 4);
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &pool_reads(&reads));
    for threads in ["1", "3"] {
        let out = dir.path().join(format!("t{threads}"));
        let template = out.join("{target}.fastq.gz");
        whittle()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", template.to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .args(["--preserve-order"])
            .assert()
            .success();
        assert_eq!(
            files_under(&out),
            ["16S.fastq.gz", "ITS.fastq.gz", "unassigned.fastq.gz"]
        );
        for bin in ["16S", "ITS", "unassigned"] {
            let text = read_gz(&out.join(format!("{bin}.fastq.gz")));
            assert_eq!(ids_with_label(&text, bin), truth(&reads, bin), "{bin}");
        }
    }
}

#[test]
fn bam_split_files_have_headers() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let reads = pool(&p, 5, 4, 3);
    let input = dir.path().join("in.bam");
    let header: sam::Header = "@HD\tVN:1.6\n@RG\tID:rg1\n".parse().unwrap();
    write_ubam_with_header(&input, &pool_reads(&reads), header);
    for threads in ["1", "2"] {
        let out = dir.path().join(format!("t{threads}"));
        whittle()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", out.join("{target}.bam").to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .args(["--preserve-order"])
            .assert()
            .success();
        assert_eq!(files_under(&out), ["16S.bam", "ITS.bam", "unassigned.bam"]);
        for bin in ["16S", "ITS", "unassigned"] {
            let path = out.join(format!("{bin}.bam"));
            let mut rdr = bam::io::Reader::new(std::fs::File::open(&path).unwrap());
            let hdr = rdr.read_header().unwrap();
            assert!(
                hdr.programs().as_ref().contains_key(&b"whittle"[..]),
                "{bin}: @PG whittle missing"
            );
            let tm = read_group_trim_mode(&path);
            assert!(
                tm.as_deref().is_some_and(|tm| tm.contains("primer")),
                "{tm:?}"
            );
            let recs = bam_records(&path);
            let names: Vec<String> = recs.iter().map(|r| r.0.clone()).collect();
            assert_eq!(names, truth(&reads, bin), "{bin}");
            for (_, _, tags) in &recs {
                assert!(
                    tags.contains(&("wt".to_string(), bin.to_string())),
                    "{tags:?}"
                );
            }
        }
    }
}

#[test]
fn group_template_merges_targets() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let (f34, r34) = (splitmix_dna(905, 18), splitmix_dna(906, 21));
    let sheet = write_rows(
        dir.path(),
        "groups.tsv",
        &[
            ("16S_full", &p.f16, &p.r16, "16S"),
            ("16S_V34", &f34, &r34, "16S"),
            ("ITS", &p.fits, &p.rits, "ITS"),
        ],
    );
    let reads = vec![
        ("full".to_string(), cat(&[&p.f16, &insert(1), &rc(&p.r16)])),
        ("v34".to_string(), cat(&[&f34, &insert(2), &rc(&r34)])),
        ("its".to_string(), cat(&[&p.fits, &insert(3), &rc(&p.rits)])),
    ];
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    let out = dir.path().join("out");
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", out.join("{group}.fq").to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .success();
    assert_eq!(files_under(&out), ["16S.fq", "ITS.fq"]);
    let text = std::fs::read_to_string(out.join("16S.fq")).unwrap();
    assert_eq!(ids_with_label(&text, "16S"), ["full", "v34"]);
    let text = std::fs::read_to_string(out.join("ITS.fq")).unwrap();
    assert_eq!(ids_with_label(&text, "ITS"), ["its"]);

    // Beside `{target}`, `{group}` names the assigned target's group.
    let out = dir.path().join("nested");
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", out.join("{group}/{target}.fq").to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .success();
    assert_eq!(
        files_under(&out),
        ["16S/16S_V34.fq", "16S/16S_full.fq", "ITS/ITS.fq"]
    );
}

/// A barcoded uBAM read: name, sequence and `BC:Z` barcode call.
type Barcoded = (String, Vec<u8>, Option<&'static str>);

/// Writes `reads` as unmapped records, each with a `BC:Z` barcode call when
/// one is given.
fn write_ubam_barcoded(path: &Path, reads: &[Barcoded]) {
    let header = sam::Header::default();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();
    for (name, seq, bc) in reads {
        let mut r = RecordBuf::default();
        *r.flags_mut() = Flags::UNMAPPED;
        *r.name_mut() = Some(name.as_bytes().into());
        *r.sequence_mut() = seq.clone().into();
        *r.quality_scores_mut() = vec![40; seq.len()].into();
        if let Some(bc) = bc {
            r.data_mut()
                .insert(Tag::new(b'B', b'C'), Value::String((*bc).into()));
        }
        w.write_alignment_record(&header, &r).unwrap();
    }
    w.try_finish().unwrap();
}

#[test]
fn barcode_placeholder_from_bc_tag() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let reads: Vec<Barcoded> = vec![
        (
            "a".into(),
            cat(&[&p.f16, &insert(1), &rc(&p.r16)]),
            Some("barcode01"),
        ),
        (
            "b".into(),
            cat(&[&p.fits, &insert(2), &rc(&p.rits)]),
            Some("barcode02"),
        ),
        ("c".into(), cat(&[&p.f16, &insert(3), &rc(&p.r16)]), None),
        ("d".into(), insert(4), Some("barcode02")),
    ];
    let input = dir.path().join("in.bam");
    write_ubam_barcoded(&input, &reads);
    let expected = [
        "barcode01.16S.bam",
        "barcode02.ITS.bam",
        "barcode02.unassigned.bam",
        "unclassified.16S.bam",
    ];
    for threads in ["1", "2"] {
        let out = dir.path().join(format!("bam{threads}"));
        whittle()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", out.join("{barcode}.{target}.bam").to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .assert()
            .success();
        assert_eq!(files_under(&out), expected);
        let names = |file: &str| -> Vec<String> {
            bam_records(&out.join(file))
                .into_iter()
                .map(|r| r.0)
                .collect()
        };
        assert_eq!(names("barcode01.16S.bam"), ["a"]);
        assert_eq!(names("barcode02.ITS.bam"), ["b"]);
        assert_eq!(names("unclassified.16S.bam"), ["c"]);
        assert_eq!(names("barcode02.unassigned.bam"), ["d"]);
    }

    // BAM to FASTQ and tagged FASTQ read the same tag; plain FASTQ has none.
    let out = dir.path().join("fq");
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", out.join("{barcode}/{target}.fq").to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .success();
    assert_eq!(
        files_under(&out),
        [
            "barcode01/16S.fq",
            "barcode02/ITS.fq",
            "barcode02/unassigned.fq",
            "unclassified/16S.fq"
        ]
    );
    let tagged = dir.path().join("tagged.fastq");
    let tagged_reads: Vec<(String, Vec<u8>)> = reads
        .iter()
        .map(|(n, s, bc)| match bc {
            Some(bc) => (format!("{n}\tBC:Z:{bc}"), s.clone()),
            None => (n.clone(), s.clone()),
        })
        .collect();
    write_fastq(&tagged, &tagged_reads);
    let out = dir.path().join("tagged");
    whittle()
        .args(["-i", tagged.to_str().unwrap()])
        .args(["-o", out.join("{barcode}.{target}.fq").to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "2"])
        .assert()
        .success();
    assert_eq!(
        files_under(&out),
        [
            "barcode01.16S.fq",
            "barcode02.ITS.fq",
            "barcode02.unassigned.fq",
            "unclassified.16S.fq"
        ]
    );
    let plain = dir.path().join("plain.fastq");
    write_fastq(&plain, &three_reads(&p));
    let out = dir.path().join("plain");
    whittle()
        .args(["-i", plain.to_str().unwrap()])
        .args(["-o", out.join("{barcode}.{target}.fq").to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .success();
    assert_eq!(
        files_under(&out),
        [
            "unclassified.16S.fq",
            "unclassified.ITS.fq",
            "unclassified.unassigned.fq"
        ]
    );
}

/// A sheet whose `alt` target shares the 16S forward primer, and an input
/// whose read `amb` carries only that primer, which is ambiguous between
/// 16S and `alt`.
fn ambiguous_setup(dir: &Path, p: &Primers) -> (PathBuf, PathBuf) {
    let ralt = splitmix_dna(907, 22);
    let sheet = write_rows(
        dir,
        "amb.tsv",
        &[
            ("16S", &p.f16, &p.r16, "16S"),
            ("ITS", &p.fits, &p.rits, "ITS"),
            ("alt", &p.f16, &ralt, "alt"),
        ],
    );
    let mut reads = three_reads(p);
    reads.push(("amb".to_string(), cat(&[&p.f16, &insert(9)])));
    let input = dir.join("amb.fastq");
    write_fastq(&input, &reads);
    (sheet, input)
}

#[test]
fn discard_ambiguous_writes_no_file() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let (sheet, input) = ambiguous_setup(dir.path(), &p);
    let kept = dir.path().join("kept");
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", kept.join("{target}.fq").to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .success();
    assert_eq!(
        files_under(&kept),
        ["16S.fq", "ITS.fq", "ambiguous.fq", "unassigned.fq"]
    );
    let text = std::fs::read_to_string(kept.join("ambiguous.fq")).unwrap();
    assert_eq!(ids_with_label(&text, "ambiguous"), ["amb"]);
    for threads in ["1", "2"] {
        let out = dir.path().join(format!("discard{threads}"));
        whittle()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", out.join("{target}.fq").to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .args(["--split-discard", "ambiguous"])
            .assert()
            .success();
        assert_eq!(files_under(&out), ["16S.fq", "ITS.fq", "unassigned.fq"]);
    }
}

/// A key opens its file only once a record routes to it: a sheet target
/// without reads, and the `unassigned` bin of a run where every read is
/// assigned, leave no file.
#[test]
fn unused_keys_create_no_file() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let (fv, rv) = (splitmix_dna(908, 20), splitmix_dna(909, 20));
    let sheet = write_rows(
        dir.path(),
        "three.tsv",
        &[
            ("16S", &p.f16, &p.r16, "16S"),
            ("ITS", &p.fits, &p.rits, "ITS"),
            ("V9", &fv, &rv, "V9"),
        ],
    );
    let reads: Vec<(String, Vec<u8>)> = three_reads(&p).into_iter().take(2).collect();
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    for threads in ["1", "2"] {
        let out = dir.path().join(format!("t{threads}"));
        whittle()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", out.join("{target}.fastq.gz").to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .assert()
            .success();
        assert_eq!(files_under(&out), ["16S.fastq.gz", "ITS.fastq.gz"]);
    }
}

#[test]
fn template_hitting_input_is_error() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_rows(
        dir.path(),
        "in.tsv",
        &[
            ("in", &p.f16, &p.r16, "in"),
            ("ITS", &p.fits, &p.rits, "ITS"),
        ],
    );
    let input = dir.path().join("in.fq");
    write_fastq(&input, &three_reads(&p));
    let before = std::fs::read(&input).unwrap();
    for threads in ["1", "2"] {
        whittle()
            .current_dir(dir.path())
            .args(["-i", "in.fq", "-o", "{target}.fq"])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .assert()
            .failure()
            .stderr(predicates::str::contains("in.fq").and(predicates::str::contains("input")));
        assert_eq!(
            std::fs::read(&input).unwrap(),
            before,
            "The input is intact"
        );
        // The collision is found before any input is read, so no other
        // key's file is written either.
        assert_eq!(files_under(dir.path()), ["in.fq", "in.tsv"]);
    }
}

/// With a template, a target name that is not one path component is
/// refused at load time; tag-only mode takes any name.
#[test]
fn template_refuses_names_that_are_not_path_components() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let input = dir.path().join("in.fq");
    write_fastq(&input, &three_reads(&p));
    for (i, bad) in ["../x", "a/b"].into_iter().enumerate() {
        let sheet = write_rows(
            dir.path(),
            &format!("bad{i}.tsv"),
            &[(bad, &p.f16, &p.r16, "g"), ("ITS", &p.fits, &p.rits, "ITS")],
        );
        let out = dir.path().join(format!("out{i}"));
        whittle()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", out.join("{target}.fq").to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
            .assert()
            .failure()
            .stderr(predicates::str::contains(format!("target {bad:?}")));
        assert!(!out.exists(), "Nothing is written");
        let tagged = dir.path().join(format!("tagged{i}.fq"));
        whittle()
            .args([
                "-i",
                input.to_str().unwrap(),
                "-o",
                tagged.to_str().unwrap(),
            ])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
            .assert()
            .success();
        let text = std::fs::read_to_string(&tagged).unwrap();
        assert_eq!(
            by_name(&fastq_records(&text), "plus16").0,
            format!("plus16\twt:Z:{bad}")
        );
    }
}

/// With a template, two target names that differ only in letter case are
/// refused at load time, before anything is written.
#[test]
fn template_refuses_names_that_differ_only_in_case() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let input = dir.path().join("in.fq");
    write_fastq(&input, &three_reads(&p));
    let sheet = write_rows(
        dir.path(),
        "case.tsv",
        &[("its", &p.f16, &p.r16, "a"), ("ITS", &p.fits, &p.rits, "b")],
    );
    let out = dir.path().join("out");
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", out.join("{target}.fq").to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .failure()
        .stderr(
            predicates::str::contains("\"its\" and \"ITS\"")
                .and(predicates::str::contains("letter case")),
        );
    assert!(!out.exists(), "Nothing is written");
}

/// Two different bins and barcodes whose paths are the same text are
/// refused, naming the path and both owners.
#[test]
fn two_owners_on_one_path_is_error() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_rows(
        dir.path(),
        "owners.tsv",
        &[("A", &p.f16, &p.r16, "A"), ("Ab", &p.fits, &p.rits, "Ab")],
    );
    let reads: Vec<Barcoded> = vec![
        (
            "a".into(),
            cat(&[&p.f16, &insert(1), &rc(&p.r16)]),
            Some("b1"),
        ),
        (
            "b".into(),
            cat(&[&p.fits, &insert(2), &rc(&p.rits)]),
            Some("1"),
        ),
    ];
    let input = dir.path().join("in.bam");
    write_ubam_barcoded(&input, &reads);
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args([
            "-o",
            dir.path().join("{target}{barcode}.bam").to_str().unwrap(),
        ])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .failure()
        .stderr(
            predicates::str::contains("Ab1.bam")
                .and(predicates::str::contains("A with barcode b1"))
                .and(predicates::str::contains("Ab with barcode 1")),
        );
}

#[test]
fn template_hitting_rejected_output_is_error() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_rows(
        dir.path(),
        "rej.tsv",
        &[
            ("rej", &p.f16, &p.r16, "rej"),
            ("ITS", &p.fits, &p.rits, "ITS"),
        ],
    );
    let input = dir.path().join("in.fq");
    write_fastq(&input, &three_reads(&p));
    let rejected = dir.path().join("rej.fq");
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", dir.path().join("{target}.fq").to_str().unwrap()])
        .args(["--rejected-output", rejected.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .failure()
        .stderr(
            predicates::str::contains("rej.fq").and(predicates::str::contains("--rejected-output")),
        );
}

/// Two keys whose paths differ as text but name one file are refused.
#[test]
fn two_keys_on_one_file_is_error() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let reads: Vec<Barcoded> = vec![
        (
            "a".into(),
            cat(&[&p.f16, &insert(1), &rc(&p.r16)]),
            Some("b1"),
        ),
        (
            "b".into(),
            cat(&[&p.f16, &insert(2), &rc(&p.r16)]),
            Some("b2"),
        ),
    ];
    let input = dir.path().join("in.bam");
    write_ubam_barcoded(&input, &reads);
    let template = dir.path().join("out/{barcode}/../{target}.bam");
    whittle()
        .args(["-i", input.to_str().unwrap()])
        .args(["-o", template.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap(), "-t", "1"])
        .assert()
        .failure()
        .stderr(predicates::str::contains("16S.bam").and(predicates::str::contains("same file")));
}

#[test]
fn template_errors_name_the_path() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let input = dir.path().join("in.fq");
    write_fastq(&input, &three_reads(&p));
    for (template, needle) in [
        ("out/{barcode}.fq", "{target}"),
        ("out/{sample}.fq", "{sample}"),
        ("out/{target}", "extension"),
    ] {
        whittle()
            .args(["-i", input.to_str().unwrap(), "-o", template])
            .args(["--split-by", sheet.to_str().unwrap()])
            .assert()
            .failure()
            .stderr(predicates::str::contains(template).and(predicates::str::contains(needle)));
    }
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", "out/{target}.fq"])
        .assert()
        .failure()
        .stderr(predicates::str::contains("--split-by"));
}

#[test]
fn preserve_order_same_across_threads() {
    let dir = tempfile::tempdir().unwrap();
    let p = primers();
    let sheet = write_sheet(dir.path(), &p);
    let reads = pool(&p, 400, 300, 200);
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &pool_reads(&reads));
    let run = |threads: &str| {
        let out = dir.path().join(format!("t{threads}"));
        whittle()
            .args(["-i", input.to_str().unwrap()])
            .args(["-o", out.join("{target}.fastq").to_str().unwrap()])
            .args(["--split-by", sheet.to_str().unwrap(), "-t", threads])
            .args(["--preserve-order"])
            .assert()
            .success();
        out
    };
    let one = run("1");
    let four = run("4");
    assert_eq!(files_under(&one), files_under(&four));
    for file in files_under(&one) {
        let a = std::fs::read(one.join(&file)).unwrap();
        let b = std::fs::read(four.join(&file)).unwrap();
        assert!(a == b, "{file} differs between -t 1 and -t 4");
    }
    // Every read is written once, and each file holds its reads in input
    // order.
    let order: Vec<String> = reads.iter().map(|r| r.0.clone()).collect();
    let mut total = 0;
    for file in files_under(&one) {
        let bin = file.trim_end_matches(".fastq");
        let text = std::fs::read_to_string(one.join(&file)).unwrap();
        let ids = ids_with_label(&text, bin);
        let positions: Vec<usize> = ids
            .iter()
            .map(|id| order.iter().position(|n| n == id).unwrap())
            .collect();
        assert!(positions.is_sorted(), "{file} is out of input order");
        total += ids.len();
    }
    assert_eq!(total, reads.len());
}

/// The TSV sheet of the degenerate 16S and ITS primer pairs, whose 16S
/// primers overlap the catalog 16S primers without equaling them.
fn write_degenerate_sheet(dir: &Path) -> PathBuf {
    let path = dir.join("degenerate.tsv");
    std::fs::write(
        &path,
        "target\tfwd\trev\tgroup\n\
         16S_full\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\t16S\n\
         ITS\tCTTGGTCATTTAGAGGAAGTAA\tTCCTCCGCTTATTGATATGC\tITS\n",
    )
    .unwrap();
    path
}

/// A pooled 16S and ITS amplicon run behind the LSK114 adapter: `n` reads
/// of each target, alternating, with 8 random bases at each end. The 16S
/// primer sites match both the catalog entries and the degenerate sheet
/// primers. Returns each read with its insert.
fn degenerate_pool(n: u64) -> Vec<(String, Vec<u8>, Vec<u8>)> {
    let lsk = b"CCTGTACTTCGTTCAGTTACGTATTGC";
    let f16 = b"AGAGTTTGATCATGGCTCAG";
    let r16 = b"TACGGTTACCTTGTTACGACTT";
    let fits = b"CTTGGTCATTTAGAGGAAGTAA";
    let rits = b"TCCTCCGCTTATTGATATGC";
    let mut reads = Vec::new();
    for i in 0..2 * n {
        let (name, body, ins) = if i % 2 == 0 {
            let ins = splitmix_dna(5000 + i, 1400);
            (format!("r{i}_16S"), cat(&[f16, &ins, &rc(r16)]), ins)
        } else {
            let ins = splitmix_dna(5000 + i, 600);
            (format!("r{i}_ITS"), cat(&[fits, &ins, &rc(rits)]), ins)
        };
        let head = splitmix_dna(7000 + i, 8);
        let tail = splitmix_dna(8000 + i, 8);
        reads.push((name, cat(&[&head, lsk, &body, &tail]), ins));
    }
    reads
}

/// Runs the degenerate sheet over `degenerate_pool(n)` with the MAB114
/// preset and `extra` flags, and checks that every read is assigned its
/// target, and with `exact`, that both primers are trimmed and the insert
/// kept whole.
fn assert_degenerate_pool_assigned(n: u64, extra: &[&str], exact: bool) {
    let dir = tempfile::tempdir().unwrap();
    let sheet = write_degenerate_sheet(dir.path());
    let pool = degenerate_pool(n);
    let reads: Vec<(String, Vec<u8>)> = pool
        .iter()
        .map(|(n, s, _)| (n.clone(), s.clone()))
        .collect();
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    let out = dir.path().join("out.fastq");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args(["--split-by", sheet.to_str().unwrap()])
        .args(["--adapter-preset", "mab114", "-t", "1"])
        .args(extra)
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    assert_eq!(recs.len(), reads.len(), "{text}");
    for (name, _, ins) in &pool {
        let key = if name.ends_with("_16S") {
            "16S_full"
        } else {
            "ITS"
        };
        let (head, seq) = by_name(&recs, name);
        assert_eq!(head, &format!("{name}\twt:Z:{key}"));
        if exact {
            assert_eq!(seq.as_bytes(), ins.as_slice(), "{name}");
        }
    }
}

/// Sheet primers that overlap non-identical catalog entries of the preset
/// are located within the catalog trims and assign every read.
#[test]
fn degenerate_sheet_under_the_mab114_preset_assigns_both_targets() {
    assert_degenerate_pool_assigned(20, &[], true);
}

/// Adapter discovery on top of the preset leaves the sheet primers located.
#[test]
fn degenerate_sheet_under_the_mab114_preset_with_discovery_assigns_both_targets() {
    assert_degenerate_pool_assigned(100, &["--adapter-discover"], false);
}

/// A target whose reverse primer is the reverse complement of its forward
/// primer is located at both ends on both strands and assigned under
/// `--split-require both`, beside an ordinary target.
#[test]
fn reverse_primer_equal_to_reverse_complemented_forward_opens_both_ends() {
    let dir = tempfile::tempdir().unwrap();
    let px = splitmix_dna(951, 22);
    let (fy, ry) = (splitmix_dna(952, 21), splitmix_dna(953, 20));
    let sheet = write_rows(
        dir.path(),
        "rc.tsv",
        &[("X", &px, &rc(&px), "X"), ("Y", &fy, &ry, "Y")],
    );
    let mut reads = Vec::new();
    for i in 0..6u64 {
        let body = if i % 2 == 0 {
            cat(&[&px, &insert(200 + i), &px])
        } else {
            cat(&[&rc(&px), &insert(200 + i), &rc(&px)])
        };
        reads.push((format!("x{i}"), body));
        reads.push((format!("y{i}"), cat(&[&fy, &insert(300 + i), &rc(&ry)])));
    }
    let input = dir.path().join("in.fastq");
    write_fastq(&input, &reads);
    let out = dir.path().join("out.fastq");
    let json = dir.path().join("summary.json");
    whittle()
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .args([
            "--split-by",
            sheet.to_str().unwrap(),
            "--split-require",
            "both",
        ])
        .args(["-t", "1", "--summary-json", json.to_str().unwrap()])
        .assert()
        .success();
    let text = std::fs::read_to_string(&out).unwrap();
    let recs = fastq_records(&text);
    for i in 0..6u64 {
        let (head, seq) = by_name(&recs, &format!("x{i}"));
        assert_eq!(head, &format!("x{i}\twt:Z:X"));
        assert_eq!(seq.as_bytes(), insert(200 + i));
        let (head, _) = by_name(&recs, &format!("y{i}"));
        assert_eq!(head, &format!("y{i}\twt:Z:Y"));
    }
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
    let x = v["split"]["keys"]
        .as_array()
        .unwrap()
        .iter()
        .find(|k| k["key"] == "X")
        .unwrap()
        .clone();
    assert_eq!(x["reads"], 6);
    assert_eq!(x["both_ends"], 6);
    assert_eq!(x["plus"], 3);
    assert_eq!(x["minus"], 3);
}
