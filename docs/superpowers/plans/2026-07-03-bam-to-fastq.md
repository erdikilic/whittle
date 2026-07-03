# BAM → FASTQ/.gz Conversion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add unaligned-BAM → FASTQ/.gz conversion (single file, stdin/stdout, and BAM folder), reconstructing MM/ML/MN into the FASTQ header alongside other aux tags, controlled by a new `--fastq-tags all|none|LIST` flag.

**Architecture:** Reuse the existing format-neutral pipeline. A new single-threaded `run_bam_to_fastq` mirrors `run_bam`'s per-record guards but writes tagged FASTQ; the MM/ML slice/serialize block is extracted from `reconstruct_record` into a shared `reconstruct_mods` so the BAM and FASTQ paths cannot drift; `io/fastq.rs` gains pure SAM-text aux serializers and a tagged writer; `lib.rs` gains new dispatch arms. gz output still parallelizes through the existing `gzp` writer.

**Tech Stack:** Rust 2024, `noodles-bam`/`noodles-sam` (BAM read + `Value` types), `gzp` (parallel gz output), `flate2` (gz decode, and gz read in tests), `clap` derive, existing `mods` codec, `rust-htslib` (dev-dep, oracle only), `assert_cmd` + `tempfile` (integration tests).

## Global Constraints

- Rust 2024 edition; every task must leave `cargo clippy --all-targets -- -D warnings` clean.
- No new runtime or dev dependencies.
- Internal quality is raw Phred (`u8`); the FASTQ boundary adds/subtracts 33.
- MM/ML/MN are the only trim-aware tags; all other aux tags are copied verbatim.
- `--fastq-tags` default is `all`. Values: `all`, `none`, or a comma list of exactly-2-char SAM tags.
- FASTQ→BAM stays a hard error. Aligned BAM stays refused via the existing `ensure_unaligned`.
- SAM aux textual form is `XX:T:VALUE`; integers of any width serialize with type code `i`; `B` arrays keep subtype `c/C/s/S/i/I/f`.
- Header tags follow the read id, TAB-separated: `@<qname>[_segment_N]\t<tag>\t<tag>…`. A read with no carried tags is written as an ordinary plain FASTQ record (no trailing TAB).
- Spec: `docs/superpowers/specs/2026-07-03-chopping-bam-to-fastq-design.md`.

---

### Task 1: `FastqTags` config type + `--fastq-tags` CLI parsing

**Files:**
- Modify: `src/config.rs` (add `FastqTags` enum, its methods, and a `fastq_tags` field on `Config`)
- Modify: `src/cli.rs` (add the `--fastq-tags` arg; wire it through `parse()` and `config_for_test`)
- Modify: `src/pipeline.rs` (add `fastq_tags: FastqTags::All` to the 7 test `Config { … }` literals so the crate compiles)

**Interfaces:**
- Produces:
  - `pub enum FastqTags { All, None, Only(std::collections::BTreeSet<[u8; 2]>) }` (in `crate::config`)
  - `impl FastqTags { pub fn parse(s: &str) -> anyhow::Result<Self>; pub fn carries(&self, tag: &[u8; 2]) -> bool; pub fn carries_mods(&self) -> bool }`
  - `Config` gains `pub fastq_tags: FastqTags`
- Consumes: nothing from earlier tasks.

- [ ] **Step 1: Write the failing tests for `FastqTags::parse` + carry rules**

Add to the bottom of `src/config.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_all_none() {
        assert_eq!(FastqTags::parse("all").unwrap(), FastqTags::All);
        assert_eq!(FastqTags::parse("none").unwrap(), FastqTags::None);
    }

    #[test]
    fn parse_list_collects_tags() {
        let t = FastqTags::parse("MM,ML,RG").unwrap();
        match t {
            FastqTags::Only(ref s) => {
                assert!(s.contains(b"MM") && s.contains(b"ML") && s.contains(b"RG"));
                assert_eq!(s.len(), 3);
            }
            other => panic!("expected Only, got {other:?}"),
        }
    }

    #[test]
    fn parse_rejects_bad_token() {
        assert!(FastqTags::parse("MM,ABC").is_err()); // 3-char token
        assert!(FastqTags::parse("").is_err()); // empty -> one empty token
        assert!(FastqTags::parse("MM,").is_err()); // trailing empty token
    }

    #[test]
    fn carries_rules() {
        assert!(FastqTags::All.carries(b"RG"));
        assert!(FastqTags::All.carries_mods());
        assert!(!FastqTags::None.carries(b"RG"));
        assert!(!FastqTags::None.carries_mods());

        let only = FastqTags::parse("ML,RG").unwrap();
        assert!(only.carries(b"RG"));
        assert!(!only.carries(b"XY"));
        // mod block carried when the list has MM *or* ML:
        assert!(only.carries_mods());
        // MN alone does not turn on the mod block:
        let mn_only = FastqTags::parse("MN").unwrap();
        assert!(!mn_only.carries_mods());
    }
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --lib config::tests 2>&1 | tail -20`
Expected: FAIL — `FastqTags` does not exist yet (compile error).

- [ ] **Step 3: Add `FastqTags` and the `Config` field**

In `src/config.rs`, add the import and type (top of file, after existing `use`s), and the field on `Config`:

```rust
use std::collections::BTreeSet;
use std::path::PathBuf;

use crate::filter::FilterConfig;
use crate::io::Format;
use crate::trim::TrimPlan;

/// Which aux tags to carry into FASTQ headers on BAM→FASTQ conversion.
/// MM/ML/MN are reconstructed (trim-aware); every other carried tag is verbatim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FastqTags {
    /// Carry every aux tag from the source record.
    All,
    /// Carry no tags — emit plain FASTQ.
    None,
    /// Carry only the listed 2-character SAM tags.
    Only(BTreeSet<[u8; 2]>),
}

impl FastqTags {
    /// Parse a `--fastq-tags` spec: `all`, `none`, or a comma list of 2-char tags.
    pub fn parse(s: &str) -> anyhow::Result<Self> {
        match s {
            "all" => Ok(FastqTags::All),
            "none" => Ok(FastqTags::None),
            _ => {
                let mut set = BTreeSet::new();
                for tok in s.split(',') {
                    let b = tok.as_bytes();
                    if b.len() != 2 {
                        anyhow::bail!(
                            "--fastq-tags: invalid tag {tok:?} (SAM tags are exactly 2 \
                             characters); use `all`, `none`, or a comma list like `MM,ML,RG`"
                        );
                    }
                    set.insert([b[0], b[1]]);
                }
                Ok(FastqTags::Only(set))
            }
        }
    }

    /// Whether a non-mod tag is carried.
    pub fn carries(&self, tag: &[u8; 2]) -> bool {
        match self {
            FastqTags::All => true,
            FastqTags::None => false,
            FastqTags::Only(s) => s.contains(tag),
        }
    }

    /// Whether the reconstructed MM/ML/MN block is carried. The block is a unit:
    /// on under `All`, or when an explicit list contains `MM` or `ML`.
    pub fn carries_mods(&self) -> bool {
        match self {
            FastqTags::All => true,
            FastqTags::None => false,
            FastqTags::Only(s) => s.contains(b"MM") || s.contains(b"ML"),
        }
    }
}
```

Then add the field to `Config`:

```rust
#[derive(Debug, Clone)]
pub struct Config {
    pub io: IoConfig,
    pub filter: FilterConfig,
    pub trim: TrimPlan,
    pub threads: usize,
    pub fastq_tags: FastqTags,
}
```

- [ ] **Step 4: Run the config tests to verify they pass**

Run: `cargo test --lib config::tests 2>&1 | tail -20`
Expected: the 4 `config::tests` pass. (The crate as a whole will NOT build yet — every `Config { … }` literal is now missing a field; that's fixed in the next steps.)

- [ ] **Step 5: Wire the `--fastq-tags` flag into `cli.rs`**

In `src/cli.rs`, add the import and the arg. Add to the `use` block:

```rust
use crate::config::{Config, FastqTags, IoConfig};
```

Add this arg to the `Cli` struct (in the `Setup` help heading group, after `threads`):

```rust
    #[arg(long, default_value = "all", help_heading = "Setup")]
    fastq_tags: String,
```

In `parse()`, before the `Ok(Config { … })`, add:

```rust
    let fastq_tags = FastqTags::parse(&c.fastq_tags)?;
```

Add `fastq_tags,` to the `Config { … }` returned by `parse()`:

```rust
        trim: TrimPlan { head: c.head_crop, tail: c.tail_crop, quality },
        threads: c.threads.max(1),
        fastq_tags,
    })
```

In `config_for_test`, add the field to its `Config { … }` literal:

```rust
        trim: TrimPlan { head: head_crop, tail: tail_crop, quality: None },
        threads: 1,
        fastq_tags: FastqTags::All,
    }
```

- [ ] **Step 6: Fix the 7 `Config { … }` literals in `src/pipeline.rs`**

Each of these test literals (at lines ~269, ~284, ~303, the `mk` closure ~319, ~367, ~387, ~478) ends with `threads: <n>,` inside `Config { … }`. Add `fastq_tags: crate::config::FastqTags::All,` immediately after the `threads:` line in each. Example for the first one:

```rust
        let cfg = Config {
            io: crate::config::IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: base_filter(),
            trim: TrimPlan { head: 1, tail: 1, quality: None },
            threads: 1,
            fastq_tags: crate::config::FastqTags::All,
        };
```

Apply the same one-line addition to all 7 literals.

- [ ] **Step 7: Verify the whole crate builds, the flag exists, and a bad value errors**

Run: `cargo build 2>&1 | tail -5`
Expected: builds clean.

Run: `cargo run -- --help 2>&1 | grep -A1 fastq-tags`
Expected: shows `--fastq-tags <FASTQ_TAGS>` with default `all`.

Run: `printf '@r\nACGT\n+\nIIII\n' | cargo run -- --fastq-tags MMM 2>&1 | tail -3`
Expected: exits non-zero with the "invalid tag" message.

- [ ] **Step 8: Run the full test suite + clippy**

Run: `cargo test 2>&1 | grep -E "test result:"`
Expected: all existing tests + the 4 new `config::tests` pass.

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -3`
Expected: no warnings.

- [ ] **Step 9: Commit**

```bash
git add src/config.rs src/cli.rs src/pipeline.rs
git commit -m "feat: add --fastq-tags selector and FastqTags config type"
```

---

### Task 2: SAM-text aux serializers + tagged FASTQ writer

**Files:**
- Modify: `src/io/fastq.rs` (add `format_aux_field`, `format_mods_aux`, `write_head`, `write_segment_tagged`; refactor `write_segment` to share `write_head`)

**Interfaces:**
- Consumes: nothing (pure formatting; uses `noodles_sam` `Value`/`Array`).
- Produces:
  - `pub fn format_aux_field(tag: [u8; 2], value: &Value) -> Vec<u8>` — one SAM aux field `XX:T:val`, no leading TAB.
  - `pub fn format_mods_aux(mm: &[u8], ml: &[u8], mn: usize) -> Vec<u8>` — `MM:Z:<mm>\tML:B:C,<ml…>\tMN:i:<mn>`, no leading TAB.
  - `pub fn write_segment_tagged<W: Write>(w, name, seq, phred, total_segments, segment_idx, tags: &[u8]) -> io::Result<()>` — like `write_segment` but inserts `tags` (empty, or a block already prefixed with TABs) between the header id and the newline.

- [ ] **Step 1: Write the failing tests**

Add these tests to the `tests` module at the bottom of `src/io/fastq.rs`:

```rust
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    #[test]
    fn aux_scalar_types() {
        assert_eq!(format_aux_field(*b"RG", &Value::String(b"grp1".as_slice().into())), b"RG:Z:grp1");
        assert_eq!(format_aux_field(*b"NM", &Value::Int32(-3)), b"NM:i:-3");
        assert_eq!(format_aux_field(*b"Uq", &Value::UInt8(200)), b"Uq:i:200");
        assert_eq!(format_aux_field(*b"pa", &Value::Float(0.5)), b"pa:f:0.5");
        assert_eq!(format_aux_field(*b"bc", &Value::Character(b'K')), b"bc:A:K");
        assert_eq!(format_aux_field(*b"H2", &Value::Hex(b"1AE3".as_slice().into())), b"H2:H:1AE3");
    }

    #[test]
    fn aux_array_subtypes() {
        assert_eq!(format_aux_field(*b"a1", &Value::Array(Array::UInt8(vec![1, 2, 3]))), b"a1:B:C,1,2,3");
        assert_eq!(format_aux_field(*b"a2", &Value::Array(Array::Int8(vec![-1, 2]))), b"a2:B:c,-1,2");
        assert_eq!(format_aux_field(*b"a3", &Value::Array(Array::Int16(vec![-5]))), b"a3:B:s,-5");
        assert_eq!(format_aux_field(*b"a4", &Value::Array(Array::UInt16(vec![5]))), b"a4:B:S,5");
        assert_eq!(format_aux_field(*b"a5", &Value::Array(Array::Int32(vec![7]))), b"a5:B:i,7");
        assert_eq!(format_aux_field(*b"a6", &Value::Array(Array::UInt32(vec![8]))), b"a6:B:I,8");
        assert_eq!(format_aux_field(*b"a7", &Value::Array(Array::Float(vec![1.5]))), b"a7:B:f,1.5");
    }

    #[test]
    fn mods_aux_layout() {
        assert_eq!(format_mods_aux(b"C+m,0;", &[10, 20], 6), b"MM:Z:C+m,0;\tML:B:C,10,20\tMN:i:6");
        // empty ML -> zero-length B:C array
        assert_eq!(format_mods_aux(b"C+m;", &[], 4), b"MM:Z:C+m;\tML:B:C\tMN:i:4");
    }

    #[test]
    fn tagged_writer_appends_tags_after_id() {
        let mut out = Vec::new();
        write_segment_tagged(&mut out, b"read2", b"AC", &[40, 40], 1, 0, b"\tRG:Z:grp1\tMM:Z:C+m,0;\tML:B:C,20\tMN:i:2").unwrap();
        assert_eq!(out, b"@read2\tRG:Z:grp1\tMM:Z:C+m,0;\tML:B:C,20\tMN:i:2\nAC\n+\nII\n");
    }

    #[test]
    fn tagged_writer_empty_tags_is_plain_record() {
        let mut a = Vec::new();
        write_segment_tagged(&mut a, b"read1", b"ACGT", &[40, 40, 40, 40], 1, 0, b"").unwrap();
        let mut b = Vec::new();
        write_segment(&mut b, b"read1", b"ACGT", &[40, 40, 40, 40], 1, 0).unwrap();
        assert_eq!(a, b);
        assert_eq!(a, b"@read1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn tagged_writer_split_suffix_then_tags() {
        let mut out = Vec::new();
        write_segment_tagged(&mut out, b"read2", b"AC", &[40, 40], 2, 1, b"\tMN:i:2").unwrap();
        assert_eq!(out, b"@read2_segment_2\tMN:i:2\nAC\n+\nII\n");
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --lib io::fastq 2>&1 | tail -20`
Expected: FAIL — `format_aux_field`, `format_mods_aux`, `write_segment_tagged` are undefined.

- [ ] **Step 3: Implement the serializers and refactor the writer**

In `src/io/fastq.rs`, add these imports near the top (after the existing `use` lines):

```rust
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
```

Replace the existing `write_segment` function with a shared `write_head` plus both writers:

```rust
/// Write the `@`-prefixed header id for a segment (no trailing newline, no tags).
/// On splits (`total_segments > 1`) the id gets a `_segment_N` suffix inserted
/// before any space-separated description, matching chopper's convention.
fn write_head<W: Write>(
    w: &mut W,
    name: &[u8],
    total_segments: usize,
    segment_idx: usize,
) -> io::Result<()> {
    w.write_all(b"@")?;
    if total_segments > 1 {
        let (id, desc) = split_head(name);
        w.write_all(id)?;
        write!(w, "_segment_{}", segment_idx + 1)?;
        if let Some(d) = desc {
            w.write_all(b" ")?;
            w.write_all(d)?;
        }
    } else {
        w.write_all(name)?;
    }
    Ok(())
}

/// Write one output segment as a plain FASTQ record. `phred` is raw; ASCII is
/// emitted by adding 33. Thin wrapper over `write_segment_tagged` with no tags,
/// so the record layout lives in exactly one place.
pub fn write_segment<W: Write>(
    w: &mut W,
    name: &[u8],
    seq: &[u8],
    phred: &[u8],
    total_segments: usize,
    segment_idx: usize,
) -> io::Result<()> {
    write_segment_tagged(w, name, seq, phred, total_segments, segment_idx, b"")
}

/// Like `write_segment`, but inserts `tags` (already TAB-prefixed per field, or
/// empty) between the header id and the newline: `@<id>[_segment_N]<tags>`.
pub fn write_segment_tagged<W: Write>(
    w: &mut W,
    name: &[u8],
    seq: &[u8],
    phred: &[u8],
    total_segments: usize,
    segment_idx: usize,
    tags: &[u8],
) -> io::Result<()> {
    write_head(w, name, total_segments, segment_idx)?;
    w.write_all(tags)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    let ascii: Vec<u8> = phred.iter().map(|&q| q.saturating_add(33)).collect();
    w.write_all(&ascii)?;
    w.write_all(b"\n")
}

/// One SAM aux field as text `XX:T:VALUE` (no leading TAB). Integers of any
/// source width serialize with SAM type code `i`; `B` arrays keep their subtype.
pub fn format_aux_field(tag: [u8; 2], value: &Value) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&tag);
    out.push(b':');
    match value {
        Value::Character(c) => {
            out.extend_from_slice(b"A:");
            out.push(*c);
        }
        Value::Int8(n) => write!(out, "i:{n}").unwrap(),
        Value::UInt8(n) => write!(out, "i:{n}").unwrap(),
        Value::Int16(n) => write!(out, "i:{n}").unwrap(),
        Value::UInt16(n) => write!(out, "i:{n}").unwrap(),
        Value::Int32(n) => write!(out, "i:{n}").unwrap(),
        Value::UInt32(n) => write!(out, "i:{n}").unwrap(),
        Value::Float(x) => write!(out, "f:{x}").unwrap(),
        Value::String(s) => {
            out.extend_from_slice(b"Z:");
            out.extend_from_slice(AsRef::<[u8]>::as_ref(s));
        }
        Value::Hex(s) => {
            out.extend_from_slice(b"H:");
            out.extend_from_slice(AsRef::<[u8]>::as_ref(s));
        }
        Value::Array(a) => {
            out.extend_from_slice(b"B:");
            write_array(&mut out, a);
        }
    }
    out
}

fn write_array(out: &mut Vec<u8>, a: &Array) {
    match a {
        Array::Int8(v) => {
            out.push(b'c');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        }
        Array::UInt8(v) => {
            out.push(b'C');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        }
        Array::Int16(v) => {
            out.push(b's');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        }
        Array::UInt16(v) => {
            out.push(b'S');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        }
        Array::Int32(v) => {
            out.push(b'i');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        }
        Array::UInt32(v) => {
            out.push(b'I');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        }
        Array::Float(v) => {
            out.push(b'f');
            for x in v {
                write!(out, ",{x}").unwrap();
            }
        }
    }
}

/// The reconstructed MM/ML/MN block as SAM aux text (no leading TAB):
/// `MM:Z:<mm>\tML:B:C,<ml…>\tMN:i:<mn>`.
pub fn format_mods_aux(mm: &[u8], ml: &[u8], mn: usize) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"MM:Z:");
    out.extend_from_slice(mm);
    out.extend_from_slice(b"\tML:B:C");
    for b in ml {
        write!(out, ",{b}").unwrap();
    }
    write!(out, "\tMN:i:{mn}").unwrap();
    out
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test --lib io::fastq 2>&1 | tail -20`
Expected: all `io::fastq` tests pass (the new ones plus the pre-existing `writes_single_segment_verbatim_header`, `split_segment_suffixes_id_before_desc`, `roundtrip_reader_writer`).

- [ ] **Step 5: Clippy**

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -3`
Expected: no warnings.

- [ ] **Step 6: Commit**

```bash
git add src/io/fastq.rs
git commit -m "feat: SAM-text aux serializers and tagged FASTQ writer"
```

---

### Task 3: `reconstruct_mods` extraction + `run_bam_to_fastq`

**Files:**
- Modify: `src/pipeline.rs` (extract `reconstruct_mods` from `reconstruct_record`; add private `build_fastq_tags`; add `run_bam_to_fastq`)

**Interfaces:**
- Consumes: `crate::config::FastqTags` + `carries`/`carries_mods` (Task 1); `crate::io::fastq::{format_aux_field, format_mods_aux, write_segment_tagged}` (Task 2).
- Produces:
  - `pub fn reconstruct_mods(src: &RecordBuf, seq: &[u8], start: usize, end: usize) -> Option<(Vec<u8>, Vec<u8>)>` — sliced `(mm, ml)`, or `None` when the source has no `MM` or nothing survives.
  - `pub fn run_bam_to_fastq<W: Write>(records: impl Iterator<Item = anyhow::Result<RecordBuf>>, writer: &mut W, cfg: &Config) -> anyhow::Result<Stats>`

- [ ] **Step 1: Write the failing tests**

Add to the `bam_tests` module at the bottom of `src/pipeline.rs`:

```rust
    use crate::config::{FastqTags, IoConfig};
    use crate::filter::FilterConfig;
    use crate::qual::QualMode;
    use crate::trim::{QualityOp, TrimPlan};

    fn cfg_bam2fq(quality: Option<QualityOp>, head: usize, tags: FastqTags) -> Config {
        Config {
            io: IoConfig { input: None, output: None, in_format: None, out_format: None },
            filter: FilterConfig {
                min_length: 1, max_length: usize::MAX, min_qual: 0.0, max_qual: 1000.0,
                min_gc: None, max_gc: None, qual_mode: QualMode::Mean,
            },
            trim: TrimPlan { head, tail: 0, quality },
            threads: 1,
            fastq_tags: tags,
        }
    }

    // "CCACCCAC" C at seq idx 0,1,3,4,5,7; MM "C+m,0,1,0" -> occ 0,2,3 -> abs 0,3,4,
    // ML [10,20,30]. head-crop 2 -> window [2,8): keeps abs 3,4 renumbered -> "C+m,0,0;" ML [20,30] MN 6.
    fn read2_with_mods_and_rg() -> RecordBuf {
        let mut rec = ubam_with_mods(b"CCACCCAC", vec![35; 8], b"C+m,0,1,0;", vec![10, 20, 30]);
        rec.data_mut().insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
        rec.data_mut().insert(
            Tag::READ_GROUP,
            Value::String(b"grp1".as_slice().into()),
        );
        rec
    }

    #[test]
    fn bam2fq_all_carries_rg_and_reconstructed_mods() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::All);
        let mut out = Vec::new();
        let stats = run_bam_to_fastq([Ok(read2_with_mods_and_rg())].into_iter(), &mut out, &cfg).unwrap();
        assert_eq!((stats.input_reads, stats.output_reads), (1, 1));
        let s = String::from_utf8(out).unwrap();
        // header carries RG verbatim + reconstructed mod block; seq head-cropped by 2.
        assert!(s.starts_with("@read1\tRG:Z:grp1\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6\n"), "got: {s:?}");
        assert!(s.contains("\nACCCAC\n+\n"), "cropped seq wrong: {s:?}");
    }

    #[test]
    fn bam2fq_only_mm_ml_drops_rg() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::parse("MM,ML").unwrap());
        let mut out = Vec::new();
        run_bam_to_fastq([Ok(read2_with_mods_and_rg())].into_iter(), &mut out, &cfg).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(!s.contains("RG:Z"), "RG must be dropped: {s:?}");
        assert!(s.contains("MM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6"), "mods missing: {s:?}");
    }

    #[test]
    fn bam2fq_none_is_plain_fastq() {
        let cfg = cfg_bam2fq(None, 2, FastqTags::None);
        let mut out = Vec::new();
        run_bam_to_fastq([Ok(read2_with_mods_and_rg())].into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@read1\nACCCAC\n+\nDDDDDD\n"); // 35+33 = 'D'
    }

    #[test]
    fn bam2fq_split_suffixes_and_segments_mods() {
        // split at the low-qual base; each segment gets its own reconstructed mods.
        let cfg = cfg_bam2fq(Some(QualityOp::Split { cutoff: 20, window: 1 }), 0, FastqTags::All);
        // seq CCAC, C+m at occ 0 and 2 -> abs 0,3; qual: good good BAD good so split [0,2),[3,4)
        let mut rec = ubam_with_mods(b"CCAC", vec![40, 40, 1, 40], b"C+m,0,1;", vec![100, 200]);
        rec.data_mut().insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(4));
        let mut out = Vec::new();
        let stats = run_bam_to_fastq([Ok(rec)].into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(stats.output_reads, 2);
        let s = String::from_utf8(out).unwrap();
        // segment 1 = [0,2) "CC" keeps abs-0 mod; segment 2 = [3,4) "C" keeps abs-3 mod.
        assert!(s.contains("@read1_segment_1\tMM:Z:C+m,0;\tML:B:C,100\tMN:i:2"), "seg1: {s:?}");
        assert!(s.contains("@read1_segment_2\tMM:Z:C+m,0;\tML:B:C,200\tMN:i:1"), "seg2: {s:?}");
    }

    #[test]
    fn bam2fq_no_mods_read_is_plain() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"plain".into());
        *rec.sequence_mut() = b"ACGT".to_vec().into();
        *rec.quality_scores_mut() = vec![40; 4].into();
        let cfg = cfg_bam2fq(None, 0, FastqTags::All);
        let mut out = Vec::new();
        run_bam_to_fastq([Ok(rec)].into_iter(), &mut out, &cfg).unwrap();
        assert_eq!(out, b"@plain\nACGT\n+\nIIII\n");
    }
```

Note: `Tag::READ_GROUP` is the noodles constant for `RG`. If it is not exported under that name, replace `Tag::READ_GROUP` with `Tag::from(*b"RG")` (both compile against noodles-sam 0.85; verify the exact constant during Step 3).

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test --lib pipeline::bam_tests 2>&1 | tail -20`
Expected: FAIL — `run_bam_to_fastq` is undefined.

- [ ] **Step 3: Extract `reconstruct_mods` and refactor `reconstruct_record`**

In `src/pipeline.rs`, add `reconstruct_mods` and change `reconstruct_record` to call it. Replace the current MM/ML rebuild block inside `reconstruct_record` (the part starting at `// Rebuild MM/ML/MN when the source carried modification tags.` through the end of that `if let Some(mm_raw) = mm_raw { … }`) with a call:

```rust
    // Rebuild MM/ML/MN when the source carried modification tags. Only touch the
    // three tags when the source actually had `MM` (preserves prior behavior:
    // a source with ML/MN but no MM is left untouched).
    if src.data().get(&Tag::BASE_MODIFICATIONS).is_some() {
        let data = out.data_mut();
        match reconstruct_mods(src, &seq, start, end) {
            Some((mm_new, ml_new)) => {
                data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm_new.into()));
                data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(ml_new)));
                data.insert(
                    Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
                    Value::Int32((end - start) as i32),
                );
            }
            None => {
                data.remove(&Tag::BASE_MODIFICATIONS);
                data.remove(&Tag::BASE_MODIFICATION_PROBABILITIES);
                data.remove(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH);
            }
        }
    }

    out
}

/// Slice MM/ML to the window `[start, end)` and re-serialize. Returns `None`
/// when the source has no `MM` tag, or when no modified position survives the
/// window (caller drops MM/ML/MN in that case). Shared by the BAM→BAM and
/// BAM→FASTQ paths so they cannot drift.
pub fn reconstruct_mods(
    src: &RecordBuf,
    seq: &[u8],
    start: usize,
    end: usize,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let mm_raw = match src.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => s.to_vec(),
        _ => return None,
    };
    let ml_raw: Vec<u8> = match src.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
        Some(Value::Array(Array::UInt8(v))) => v.clone(),
        _ => Vec::new(),
    };
    let parsed = mods::parse(&mm_raw, &ml_raw);
    let sliced = mods::reconstruct(&parsed, seq, start, end);
    let (mm_new, ml_new) = mods::serialize(&sliced);
    if mm_new.is_empty() {
        None
    } else {
        Some((mm_new, ml_new))
    }
}
```

(Keep the earlier part of `reconstruct_record` — the `src.clone()`, seq/qual slicing, and name-suffix logic — exactly as-is. The block above replaces only the MM/ML rebuild tail and closes the function, then defines `reconstruct_mods`.)

- [ ] **Step 4: Verify the existing BAM tests still pass (behavior preserved)**

Run: `cargo test --lib pipeline::bam_tests::slices_seq_qual_and_rebuilds_tags pipeline::bam_tests::split_suffixes_name_and_drops_empty_mods 2>&1 | tail -10`
Expected: both pass — the extraction is behavior-preserving.

- [ ] **Step 5: Add `build_fastq_tags` and `run_bam_to_fastq`**

Add to `src/pipeline.rs` (after `run_bam`). First extend the imports at the top of the file:

```rust
use crate::config::{Config, FastqTags};
use crate::io::fastq::{format_aux_field, format_mods_aux, write_segment, write_segment_tagged};
```

(Replace the existing `use crate::config::Config;` and `use crate::io::fastq::write_segment;` lines with the two above.)

Then add:

```rust
/// Assemble the TAB-prefixed aux-tag block for one interval: carried non-mod
/// tags verbatim in source order, then the reconstructed MM/ML/MN block. Empty
/// when nothing is carried (caller writes a plain FASTQ record).
fn build_fastq_tags(
    src: &RecordBuf,
    seq: &[u8],
    start: usize,
    end: usize,
    sel: &FastqTags,
) -> Vec<u8> {
    let mut tags = Vec::new();
    for (tag, value) in src.data().iter() {
        let t = <[u8; 2]>::from(tag);
        if t == *b"MM" || t == *b"ML" || t == *b"MN" {
            continue; // handled by the reconstructed block below
        }
        if sel.carries(&t) {
            tags.push(b'\t');
            tags.extend_from_slice(&format_aux_field(t, value));
        }
    }
    if sel.carries_mods()
        && let Some((mm, ml)) = reconstruct_mods(src, seq, start, end)
    {
        tags.push(b'\t');
        tags.extend_from_slice(&format_mods_aux(&mm, &ml, end - start));
    }
    tags
}

/// Single-threaded uBAM→FASTQ pipeline: refuse aligned reads, filter, trim, then
/// write each surviving segment as FASTQ with the selected aux tags in the header
/// (MM/ML/MN reconstructed; others verbatim). gz compression, when requested, is
/// handled by the parallel `gzp` writer this drains into.
pub fn run_bam_to_fastq<W>(
    records: impl Iterator<Item = anyhow::Result<RecordBuf>>,
    writer: &mut W,
    cfg: &Config,
) -> anyhow::Result<Stats>
where
    W: Write,
{
    let mut stats = Stats::default();
    for rec in records {
        let rec = rec?;
        crate::io::bam::ensure_unaligned(&rec)?;
        stats.input_reads += 1;

        let seq = rec.sequence().as_ref().to_vec();
        let qual = rec.quality_scores().as_ref().to_vec();
        if qual.len() != seq.len() {
            let name = rec
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                .unwrap_or_else(|| "<unnamed>".to_string());
            anyhow::bail!(
                "read {name}: BAM record SEQ length {} != QUAL length {} \
                 (records without full per-base quality are not supported)",
                seq.len(),
                qual.len()
            );
        }
        if !filter::passes(&seq, &qual, &cfg.filter) {
            continue;
        }
        let name = rec.name().map(|n| n.to_vec()).unwrap_or_default();
        let intervals = trim::apply(seq.len(), &qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            let tags = build_fastq_tags(&rec, &seq, s, e, &cfg.fastq_tags);
            if tags.is_empty() {
                write_segment(writer, &name, &seq[s..e], &qual[s..e], total, idx)?;
            } else {
                write_segment_tagged(writer, &name, &seq[s..e], &qual[s..e], total, idx, &tags)?;
            }
            stats.output_reads += 1;
        }
    }
    Ok(stats)
}
```

Note: `build_fastq_tags` uses the `let-chain` (`if a && let Some(x) = b`) which is stable on Rust 2024 (already used elsewhere in this crate, e.g. `lib.rs` `run()`).

- [ ] **Step 6: Run the new BAM→FASTQ tests**

Run: `cargo test --lib pipeline::bam_tests 2>&1 | tail -20`
Expected: the 5 new `bam2fq_*` tests pass alongside the existing ones.

- [ ] **Step 7: Full suite + clippy**

Run: `cargo test 2>&1 | grep -E "test result:"`
Expected: all pass.

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -3`
Expected: no warnings.

- [ ] **Step 8: Commit**

```bash
git add src/pipeline.rs
git commit -m "feat: run_bam_to_fastq + shared reconstruct_mods extraction"
```

---

### Task 4: Dispatch wiring in `lib.rs` (`run` + `run_folder`) + ignored-note

**Files:**
- Modify: `src/lib.rs` (`run()` match arms; `run_folder()` Bam arm; add `note_tags_ignored` helper)

**Interfaces:**
- Consumes: `pipeline::run_bam_to_fastq` (Task 3); `config::FastqTags` (Task 1); existing `io::bam::reader`, `io::dir::bam_reader`, `fastq_writer`, `FastqOut::finish`.
- Produces: end-to-end BAM→FASTQ/.gz for single file and folder.

- [ ] **Step 1: Replace the cross-format bail in `run()`**

In `src/lib.rs`, replace the `match (in_fmt, out_fmt) { … }` block (currently lines ~71–90) with:

```rust
    match (in_fmt, out_fmt) {
        (Format::Bam, Format::Bam) => {
            note_tags_ignored(&cfg, in_fmt, out_fmt);
            let (header, records) = io::bam::reader(in_path)?;
            let out_header = provenance_header(header);
            let mut writer = io::bam::writer(cfg.io.output.as_deref(), &out_header)?;
            let stats = pipeline::run_bam(&out_header, records, &mut writer, &cfg)?;
            writer.try_finish()?;
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            return Ok(());
        }
        (Format::Bam, Format::Fastq | Format::FastqGz) => {
            let (_header, records) = io::bam::reader(in_path)?;
            let mut writer = fastq_writer(&cfg, out_fmt)?;
            let stats = pipeline::run_bam_to_fastq(records, &mut writer, &cfg)?;
            writer.finish()?;
            eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
            return Ok(());
        }
        (Format::Fastq | Format::FastqGz, Format::Bam) => {
            anyhow::bail!("cross-format FASTQ->BAM conversion is not supported")
        }
        _ => {}
    }

    note_tags_ignored(&cfg, in_fmt, out_fmt);
    let mut writer = fastq_writer(&cfg, out_fmt)?;
```

(The final `let mut writer = fastq_writer(&cfg, out_fmt)?;` already exists just below the match — do not duplicate it; the snippet shows it only to anchor where `note_tags_ignored` goes. Place the `note_tags_ignored(&cfg, in_fmt, out_fmt);` call immediately before the existing `let mut writer = fastq_writer(&cfg, out_fmt)?;` line.)

- [ ] **Step 2: Add the `note_tags_ignored` helper**

Add near `provenance_header` in `src/lib.rs`:

```rust
/// `--fastq-tags` only affects BAM→FASTQ output. When the user set a non-default
/// value (`none`/an explicit list) on any other path, emit a one-line stderr note
/// rather than silently ignoring it. (An explicit `all` is the default and stays
/// silent.)
fn note_tags_ignored(cfg: &Config, in_fmt: io::Format, out_fmt: io::Format) {
    if !matches!(cfg.fastq_tags, config::FastqTags::All) {
        eprintln!(
            "note: --fastq-tags applies only to BAM->FASTQ output; ignored for {in_fmt:?}->{out_fmt:?}"
        );
    }
}
```

- [ ] **Step 3: Wire the BAM folder → FASTQ path in `run_folder()`**

In `src/lib.rs`, replace the `io::dir::Family::Bam => { … }` arm of `run_folder` with:

```rust
        io::dir::Family::Bam => match out_fmt {
            Format::Bam => {
                note_tags_ignored(cfg, family_fmt, out_fmt);
                let (header, records) = io::dir::bam_reader(&paths)?;
                let out_header = provenance_header(header);
                let mut writer = io::bam::writer(cfg.io.output.as_deref(), &out_header)?;
                let stats = pipeline::run_bam(&out_header, records, &mut writer, cfg)?;
                writer.try_finish()?;
                eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
                Ok(())
            }
            Format::Fastq | Format::FastqGz => {
                let (_header, records) = io::dir::bam_reader(&paths)?;
                let mut writer = fastq_writer(cfg, out_fmt)?;
                let stats = pipeline::run_bam_to_fastq(records, &mut writer, cfg)?;
                writer.finish()?;
                eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
                Ok(())
            }
        },
```

Also add `note_tags_ignored(cfg, family_fmt, out_fmt);` at the top of the `io::dir::Family::Fastq =>` arm's success path (just before `let mut writer = fastq_writer(cfg, out_fmt)?;`), so a FASTQ folder with a stray `--fastq-tags` also notes it.

- [ ] **Step 4: Manual end-to-end smoke test**

Build, then convert a hand-built uBAM to FASTQ and gz:

Run:
```bash
cargo build 2>&1 | tail -1
cat > /tmp/mk_ubam.rs <<'EOF'
// (informal) — instead use the existing test fixture path below.
EOF
# Use the integration fixture instead (added in Task 5); for now just check dispatch compiles + errors correctly:
printf '@r\nACGT\n+\nIIII\n' | cargo run -- --out-format bam 2>&1 | tail -2
```
Expected: the FASTQ→BAM run prints `cross-format FASTQ->BAM conversion is not supported` and exits non-zero. (Full BAM→FASTQ end-to-end is covered by Task 5's integration tests.)

- [ ] **Step 5: Full suite + clippy**

Run: `cargo test 2>&1 | grep -E "test result:"`
Expected: all pass.

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -3`
Expected: no warnings.

- [ ] **Step 6: Commit**

```bash
git add src/lib.rs
git commit -m "feat: dispatch BAM->FASTQ/.gz in run() and run_folder()"
```

---

### Task 5: Integration tests (assert_cmd, hermetic)

**Files:**
- Create: `tests/bam_to_fastq.rs`

**Interfaces:**
- Consumes: the compiled `chopping` binary (via `assert_cmd`), the same in-tmp BAM-fixture pattern used by `tests/bam_smoke.rs`, and `flate2` for reading `.gz` output.

- [ ] **Step 1: Write the integration test file**

Create `tests/bam_to_fastq.rs`:

```rust
// End-to-end BAM→FASTQ/.gz conversion over the compiled binary. Builds a small
// uBAM fixture (a plain read, and a read with RG + MM/ML/MN mods), converts it,
// and checks header tags. The load-bearing correctness check is `cross_check`:
// the FASTQ-header MM/ML/MN must equal what the BAM→BAM path writes (which is
// itself htslib-oracle-verified by tests/bam_mods_oracle.rs).
use std::io::Read;
use std::path::Path;

use assert_cmd::Command;
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::RecordBuf};
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

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
    d.insert(Tag::BASE_MODIFICATIONS, Value::String(b"C+m,0,1,0;".to_vec().into()));
    d.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(vec![10, 20, 30])));
    d.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(8));
    w.write_alignment_record(&header, &r2).unwrap();

    w.try_finish().unwrap();
}

fn run(args: &[&str], input: &Path, output: &Path) {
    Command::cargo_bin("chopping")
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

    run(&["--out-format", "fastq", "--head-crop", "2", "--fastq-tags", "none"], &inp, &out);

    let s = std::fs::read_to_string(&out).unwrap();
    assert_eq!(read2_header_line(&s), "@read2"); // no tags
    assert!(!s.contains("MM:Z"), "mods must be dropped under none: {s:?}");
}

#[test]
fn bam_to_fastq_only_mm_ml_drops_rg() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.bam");
    let out = dir.path().join("out.fastq");
    write_fixture(&inp);

    run(&["--out-format", "fastq", "--head-crop", "2", "--fastq-tags", "MM,ML"], &inp, &out);

    let s = std::fs::read_to_string(&out).unwrap();
    assert_eq!(read2_header_line(&s), "@read2\tMM:Z:C+m,0,0;\tML:B:C,20,30\tMN:i:6");
}

#[test]
fn bam_to_fastq_gz_roundtrips() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.bam");
    let out = dir.path().join("out.fastq.gz");
    write_fixture(&inp);

    run(&["--out-format", "fastq-gz", "--head-crop", "2", "-t", "4"], &inp, &out);

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
fn fastq_tags_on_fastq_input_prints_ignored_note() {
    let dir = tempfile::tempdir().unwrap();
    let inp = dir.path().join("in.fastq");
    let out = dir.path().join("out.fastq");
    std::fs::write(&inp, b"@r\nACGT\n+\nIIII\n").unwrap();

    Command::cargo_bin("chopping")
        .unwrap()
        .args(["--fastq-tags", "none", "-i"])
        .arg(&inp)
        .arg("-o")
        .arg(&out)
        .assert()
        .success()
        .stderr(predicates::str::contains("--fastq-tags applies only to BAM->FASTQ"));
}
```

Note on `predicates`: `assert_cmd` re-exports it, but the crate may not be a direct dev-dependency. If `predicates::str` does not resolve, replace the last assertion with a manual check: capture output via `.assert().success()` then read `String::from_utf8(cmd.output().unwrap().stderr)`. Verify `predicates` availability during Step 2 and adjust.

- [ ] **Step 2: Run the integration tests**

Run: `cargo test --test bam_to_fastq 2>&1 | tail -30`
Expected: all 6 tests pass. If `predicates` is unavailable, apply the fallback noted above and re-run.

- [ ] **Step 3: Full suite + clippy**

Run: `cargo test 2>&1 | grep -E "test result:"`
Expected: all pass.

Run: `cargo clippy --all-targets -- -D warnings 2>&1 | tail -3`
Expected: no warnings.

- [ ] **Step 4: Commit**

```bash
git add tests/bam_to_fastq.rs
git commit -m "test: end-to-end BAM->FASTQ/.gz conversion + cross-check vs BAM path"
```

---

### Task 6: Documentation (README)

**Files:**
- Modify: `README.md` (support matrix / conversion section, `--fastq-tags` docs, stale-tag caveat)

**Interfaces:** none (docs only).

- [ ] **Step 1: Read the current README I/O section**

Run: `grep -n -iE 'bam|fastq|convert|out-format|support' README.md | head -40`
Identify the section describing input/output formats and where a support matrix belongs.

- [ ] **Step 2: Add a conversion / support-matrix subsection**

Add a subsection (place it near the existing format/output documentation) with this content:

````markdown
### Format conversion

`chopping` auto-detects formats from file extensions (or `--in-format` /
`--out-format`). Supported conversions:

| input → output | FASTQ | FASTQ.gz | BAM |
|----------------|:-----:|:--------:|:---:|
| FASTQ / FASTQ.gz | ✅ | ✅ | ❌ |
| unaligned BAM   | ✅ | ✅ | ✅ |

When no `-o` extension / `--out-format` is given, output **mirrors** the input
(BAM stays BAM), except a `.gz` input defaults to plain FASTQ (never
auto-compressed). FASTQ→BAM is not supported (there is no header/tags to build a
BAM from).

#### BAM → FASTQ tags (`--fastq-tags`)

On BAM→FASTQ, aux tags are written into the FASTQ header, tab-delimited, in the
`samtools fastq -T` / `samtools import -T` convention
(`@read\tMM:Z:…\tML:B:C,…`). MM/ML/MN are **reconstructed** for the trimmed
segment; every other tag is copied **verbatim**.

```
--fastq-tags all     # default: carry every aux tag
--fastq-tags none    # plain FASTQ, no tags
--fastq-tags MM,ML   # only the (reconstructed) modification tags
--fastq-tags MM,ML,RG
```

> **Caveat:** only MM/ML/MN are trim-aware. Some other aux tags are themselves
> position-dependent (e.g. dorado's `mv` move table, `ts`/`ns` signal-trim
> counts); under `--fastq-tags all` they are copied verbatim and will be **stale**
> after trimming — the same behaviour as `samtools fastq -T`. Exclude them with an
> explicit `--fastq-tags` list if that matters for your downstream.
````

- [ ] **Step 3: Verify no stale claims remain**

Run: `grep -n -iE 'cross-format|not supported|BAM.*FASTQ|FASTQ.*BAM' README.md`
Expected: any statement that BAM↔FASTQ conversion is unsupported is either removed or narrowed to FASTQ→BAM only. Fix any that still claim BAM→FASTQ is unsupported.

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "docs: document BAM->FASTQ conversion and --fastq-tags"
```

---

## Final verification (after all tasks)

- [ ] `cargo test 2>&1 | grep -E "test result:"` — all green (≈65 prior + ~20 new).
- [ ] `cargo clippy --all-targets -- -D warnings` — clean.
- [ ] `cargo fmt --check` — clean (run `cargo fmt` if not).
- [ ] Manual: build a uBAM with mods, `chopping -i x.bam -o y.fastq.gz -H 5 -t 8`, `zcat y.fastq.gz | head` shows reconstructed `MM:Z:`/`ML:B:C,`/`MN:i:` in headers.
- [ ] Then use superpowers:finishing-a-development-branch to merge `feat/bam-to-fastq` → main and push.
