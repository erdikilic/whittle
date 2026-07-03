# chopping v1 — Plan 2: uBAM + MM/ML reconstruction — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. **Prerequisite: Plan 1 (FASTQ core) is implemented and green** — this plan reuses `filter`, `trim`, `Config`, and the pipeline seams unchanged.

**Goal:** Add unaligned-BAM (uBAM) input/output that runs the same filtering + trim/split as the FASTQ path, and **recomputes the MM/ML/MN base-modification tags** for every trimmed or split output read, verified against the htslib oracle.

**Architecture:** A standalone MM/ML/MN codec (`src/mods/`) — parse → reconstruct-over-interval → serialize — that operates on the **raw** `MM:Z` string and `ML:B,C` array (bypassing noodles' typed base-mod parser, which drops mods on some "unspecified" codes). The BAM I/O layer reads owned `RecordBuf`s, refuses aligned records, and the pipeline reuses Plan 1's `filter::passes` + `trim::apply` to get intervals, then reconstructs each output `RecordBuf` (sliced SEQ/QUAL + rebuilt tags). Correctness is pinned by a decode-equivalence test against `rust-htslib::basemods_iter`.

**Tech Stack (added to Plan 1):** `noodles-bam` 0.90, `noodles-sam` 0.85, `noodles-bgzf` 0.47 (`libdeflate`), `bstr` 1; dev: `rust-htslib` 1.0 (`default-features = false`, `libdeflate`).

Design source: spec `docs/superpowers/specs/2026-07-03-chopping-tag-aware-trimmer-design.md`, §"MM/ML/MN reconstruction". API usage verified against the bench crates (`bench/crates/bam-noodles`, `bench/crates/bam-mods-cmp`) and the noodles source in the cargo cache.

## Global Constraints

- Everything from Plan 1's Global Constraints still holds (Rust 2024; internal quality = raw Phred; three quality ops mutually exclusive; etc.).
- **uBAM only.** Any record with the unmapped flag **clear** (aligned) is a hard error naming the read: aligned BAM is deferred.
- Internal quality for BAM is already raw Phred: `RecordBuf::quality_scores().as_ref()` yields raw bytes — **no ±33 conversion** on the BAM path.
- BAM SEQ from `RecordBuf::sequence().as_ref()` is decoded uppercase `ACGTN` bytes.
- MM/ML are read/written as **raw tags** via `Data::get`/`insert` with `Value::String` / `Value::Array(Array::UInt8(..))`; the typed `BaseModifications::parse` path is never used.
- `MN` (`Tag::BASE_MODIFICATION_SEQUENCE_LENGTH`) is set to the **output** segment length on every emitted read.
- `bam::io::Writer::new(w)` wraps `w` in a BGZF encoder (libdeflate via the `noodles-bgzf` feature). Write records with `writer.write_alignment_record(&header, &record_buf)` (bring `use noodles_sam::alignment::io::Write as _;` into scope).
- Every task ends green: `cargo test` + `cargo clippy -- -D warnings`.

## ML ordering contract (load-bearing — the whole codec depends on it)

For one MM group with `C` codes and `P` listed (modified) positions, the group owns `P × C` `ML` bytes, **position-major**: `[pos0·code0, pos0·code1, …, pos1·code0, …]`. Keeping/dropping a position keeps/drops its whole `C`-byte run together. The htslib oracle (Task 6) is what proves this ordering is right on real data.

## File Structure (additions to Plan 1)

```
src/
  mods/
    mod.rs           types (ModCode, MmGroup, Mods) + counting_base + re-exports
    parse.rs         parse(mm: &[u8], ml: &[u8]) -> Mods
    reconstruct.rs   reconstruct(&Mods, seq, start, end) -> Mods
    serialize.rs     serialize(&Mods) -> (Vec<u8> /*MM*/, Vec<u8> /*ML*/)
  io/
    bam.rs           reader (Header + RecordBuf iter, aligned refusal) + writer
pipeline.rs          + run_bam(header, records, writer, cfg)
lib.rs               run() dispatches BAM->BAM to run_bam
tests/
  bam_mods_oracle.rs htslib decode-equivalence on a uBAM fixture
test-data/
  mods_small.ubam    tiny synthetic uBAM with known C+m mods (built in Task 6)
```

---

### Task 1: Mods model + parser (`src/mods/mod.rs`, `src/mods/parse.rs`)

**Files:**
- Create: `src/mods/mod.rs`, `src/mods/parse.rs`
- Modify: `src/lib.rs` (add `pub mod mods;`)
- Modify: `Cargo.toml` (add `bstr = "1"`, noodles deps — see below)

**Interfaces:**
- Produces:
  - `pub enum ModCode { Char(u8), Chebi(u32) }` (derive `Debug, Clone, PartialEq, Eq`)
  - `pub struct MmGroup { pub base: u8, pub strand: u8, pub codes: Vec<ModCode>, pub status: Option<u8>, pub deltas: Vec<usize>, pub ml: Vec<u8> }` (derive `Debug, Clone, PartialEq, Eq`)
  - `pub struct Mods { pub groups: Vec<MmGroup> }` (derive `Debug, Clone, PartialEq, Eq`)
  - `pub fn complement(base: u8) -> u8`
  - `pub fn counting_base(base: u8, strand: u8) -> u8`
  - `pub fn parse(mm: &[u8], ml: &[u8]) -> Mods`

- [ ] **Step 1: Add deps to `Cargo.toml`**

```toml
# [dependencies] additions
noodles-bam = "0.90.0"
noodles-sam = "0.85.0"
noodles-bgzf = { version = "0.47.0", features = ["libdeflate"] }
bstr = "1"

# [dev-dependencies] additions
rust-htslib = { version = "1.0", default-features = false, features = ["libdeflate"] }
```

Add `pub mod mods;` to `src/lib.rs`.

- [ ] **Step 2: Write `src/mods/mod.rs` (types + helpers)**

```rust
pub mod parse;
pub mod reconstruct;
pub mod serialize;

pub use parse::parse;
pub use reconstruct::reconstruct;
pub use serialize::serialize;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModCode {
    Char(u8),
    Chebi(u32),
}

/// One MM group, e.g. `C+m?,5,12` with its slice of ML bytes.
/// `ml.len() == deltas.len() * codes.len()` (position-major, see plan header).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmGroup {
    pub base: u8,
    pub strand: u8,
    pub codes: Vec<ModCode>,
    pub status: Option<u8>,
    pub deltas: Vec<usize>,
    pub ml: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Mods {
    pub groups: Vec<MmGroup>,
}

pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        b'U' | b'u' => b'A',
        _ => b'N',
    }
}

/// The SEQ base whose occurrences the MM skip-counts index: the fundamental base
/// for `+`, its complement for `-` (the mods sit on the opposite strand). Slicing
/// only needs to count the SAME base the encoder counted — the htslib oracle
/// confirms this matches real data.
pub fn counting_base(base: u8, strand: u8) -> u8 {
    if strand == b'-' { complement(base) } else { base.to_ascii_uppercase() }
}
```

- [ ] **Step 3: Write failing tests in `src/mods/parse.rs`**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::ModCode;

    #[test]
    fn single_group_single_code() {
        let m = parse(b"C+m?,5,12,0;", &[200, 10, 128]);
        assert_eq!(m.groups.len(), 1);
        let g = &m.groups[0];
        assert_eq!((g.base, g.strand), (b'C', b'+'));
        assert_eq!(g.codes, vec![ModCode::Char(b'm')]);
        assert_eq!(g.status, Some(b'?'));
        assert_eq!(g.deltas, vec![5, 12, 0]);
        assert_eq!(g.ml, vec![200, 10, 128]);
    }

    #[test]
    fn multi_code_group_takes_two_ml_per_position() {
        // C+mh with 2 positions -> 4 ML bytes, position-major.
        let m = parse(b"C+mh,1,3;", &[10, 20, 30, 40]);
        let g = &m.groups[0];
        assert_eq!(g.codes, vec![ModCode::Char(b'm'), ModCode::Char(b'h')]);
        assert_eq!(g.deltas, vec![1, 3]);
        assert_eq!(g.ml, vec![10, 20, 30, 40]);
    }

    #[test]
    fn chebi_numeric_code() {
        let m = parse(b"C+16061,2;", &[99]);
        assert_eq!(m.groups[0].codes, vec![ModCode::Chebi(16061)]);
        assert_eq!(m.groups[0].deltas, vec![2]);
    }

    #[test]
    fn two_groups_split_ml() {
        let m = parse(b"C+m,0;A+a,1,4;", &[1, 2, 3]);
        assert_eq!(m.groups.len(), 2);
        assert_eq!(m.groups[0].ml, vec![1]); // 1 position
        assert_eq!(m.groups[1].ml, vec![2, 3]); // 2 positions
        assert_eq!(m.groups[1].base, b'A');
    }

    #[test]
    fn no_status_and_empty_positions() {
        let m = parse(b"C+m;", &[]);
        let g = &m.groups[0];
        assert_eq!(g.status, None);
        assert!(g.deltas.is_empty());
        assert!(g.ml.is_empty());
    }
}
```

- [ ] **Step 4: Run to verify failure**

Run: `cargo test mods::parse`
Expected: FAIL.

- [ ] **Step 5: Implement `src/mods/parse.rs`**

```rust
use super::{MmGroup, ModCode, Mods};

/// Parse a raw MM:Z string plus its ML:B,C array into groups. Malformed tails are
/// tolerated (best-effort): parsing a group stops at the first unexpected byte.
pub fn parse(mm: &[u8], ml: &[u8]) -> Mods {
    let mut groups = Vec::new();
    let mut ml_pos = 0usize;

    for token in mm.split(|&b| b == b';') {
        if token.len() < 2 {
            continue; // empty (trailing ';') or malformed
        }
        let base = token[0];
        let strand = token[1];
        let mut i = 2;

        // Codes: either a run of letters (each one code) or a numeric ChEBI id.
        let mut codes = Vec::new();
        if i < token.len() && token[i].is_ascii_digit() {
            let mut id = 0u32;
            while i < token.len() && token[i].is_ascii_digit() {
                id = id * 10 + (token[i] - b'0') as u32;
                i += 1;
            }
            codes.push(ModCode::Chebi(id));
        } else {
            while i < token.len() && token[i].is_ascii_alphabetic() {
                codes.push(ModCode::Char(token[i]));
                i += 1;
            }
        }

        // Optional status flag.
        let mut status = None;
        if i < token.len() && (token[i] == b'.' || token[i] == b'?') {
            status = Some(token[i]);
            i += 1;
        }

        // Skip-count deltas: (',' number)*
        let mut deltas = Vec::new();
        while i < token.len() {
            if token[i] != b',' {
                break;
            }
            i += 1;
            let mut n = 0usize;
            let mut saw = false;
            while i < token.len() && token[i].is_ascii_digit() {
                n = n * 10 + (token[i] - b'0') as usize;
                i += 1;
                saw = true;
            }
            if saw {
                deltas.push(n);
            }
        }

        // Claim this group's ML bytes: positions * codes, position-major.
        let want = deltas.len() * codes.len().max(1);
        let end = (ml_pos + want).min(ml.len());
        let group_ml = ml[ml_pos..end].to_vec();
        ml_pos = end;

        groups.push(MmGroup { base, strand, codes, status, deltas, ml: group_ml });
    }

    Mods { groups }
}
```

Append Step 3's tests.

- [ ] **Step 6: Run and commit**

Run: `cargo test mods::parse`
Expected: PASS (5 tests).

```bash
git add Cargo.toml Cargo.lock src/lib.rs src/mods/mod.rs src/mods/parse.rs
git commit -m "feat: MM/ML tag model + raw parser (bypasses noodles typed parser)"
```

---

### Task 2: Reconstruct over interval (`src/mods/reconstruct.rs`)

**Files:**
- Create: `src/mods/reconstruct.rs`

**Interfaces:**
- Consumes: `super::{MmGroup, Mods, counting_base}`
- Produces: `pub fn reconstruct(mods: &Mods, seq: &[u8], start: usize, end: usize) -> Mods` — the mods for the sub-read `seq[start..end)`. Groups with no surviving positions are dropped.

- [ ] **Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::{ModCode, MmGroup, Mods, parse};

    // Helper: build from raw MM/ML then reconstruct.
    fn recon(mm: &[u8], ml: &[u8], seq: &[u8], start: usize, end: usize) -> Mods {
        reconstruct(&parse(mm, ml), seq, start, end)
    }

    #[test]
    fn keeps_only_in_window_and_renumbers() {
        // seq C at indices 0,2,5,8. MM "C+m,0,1,0" -> modified Cs at occ 0,2,3 => pos 0,5,8.
        // ML one byte per position: [11,22,33].
        let seq = b"CACatCttC"; // upper only counts: positions of 'C' = 0,2,5,8 (index2 'C', 't' lower ignored? use uppercase seq)
        let seq = b"CACATCTTC"; // 'C' at 0,2,5,8
        let m = recon(b"C+m,0,1,0;", &[11, 22, 33], seq, 3, 9);
        // window [3,9): C occurrences at 5,8. Surviving modified: pos5 (occ idx1), pos8 (occ idx2).
        // window C positions = [5,8] -> renumber. pos5 -> delta 0; pos8 -> delta 0.
        assert_eq!(m.groups.len(), 1);
        assert_eq!(m.groups[0].deltas, vec![0, 0]);
        assert_eq!(m.groups[0].ml, vec![22, 33]);
    }

    #[test]
    fn drops_group_with_no_survivors() {
        let seq = b"CCCC";
        let m = recon(b"C+m,0;", &[50], seq, 2, 4); // modified C at pos0, outside [2,4)
        assert!(m.groups.is_empty());
    }

    #[test]
    fn multi_code_keeps_both_ml_bytes_per_position() {
        let seq = b"CC";
        // C+mh with 2 positions (0 and next) -> ML [a0,b0,a1,b1]. Keep both in [0,2).
        let m = recon(b"C+mh,0,0;", &[1, 2, 3, 4], seq, 0, 2);
        assert_eq!(m.groups[0].ml, vec![1, 2, 3, 4]);
        assert_eq!(m.groups[0].deltas, vec![0, 0]);
    }

    #[test]
    fn minus_strand_counts_complement() {
        // G-m: counting base = complement(G) = C. seq C at 1,3. "G-m,0,0" -> modified at C occ 0,1 => pos1,3.
        let seq = b"ACAC";
        let m = recon(b"G-m,0,0;", &[7, 8], seq, 2, 4); // window keeps pos3 only
        assert_eq!(m.groups[0].deltas, vec![0]);
        assert_eq!(m.groups[0].ml, vec![8]);
        assert_eq!((m.groups[0].base, m.groups[0].strand), (b'G', b'-'));
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test mods::reconstruct`
Expected: FAIL.

- [ ] **Step 3: Implement `src/mods/reconstruct.rs`**

```rust
use super::{MmGroup, Mods, counting_base};

pub fn reconstruct(mods: &Mods, seq: &[u8], start: usize, end: usize) -> Mods {
    let mut out = Vec::new();

    for g in &mods.groups {
        let ncodes = g.codes.len().max(1);
        let cbase = counting_base(g.base, g.strand);

        // All positions of the counting base along the whole SEQ (ascending).
        let positions: Vec<usize> = seq
            .iter()
            .enumerate()
            .filter(|(_, &b)| b.to_ascii_uppercase() == cbase)
            .map(|(i, _)| i)
            .collect();

        // Positions inside the output window (ascending), for renumbering.
        let window: Vec<usize> = positions
            .iter()
            .copied()
            .filter(|&p| p >= start && p < end)
            .collect();

        // Walk the group's deltas to recover each modified absolute position and
        // its ML byte run, then keep the ones inside the window.
        let mut new_deltas = Vec::new();
        let mut new_ml = Vec::new();
        let mut prev_widx: isize = -1;
        let mut cursor = 0usize; // index into `positions`

        for (k, &d) in g.deltas.iter().enumerate() {
            cursor += d;
            if cursor >= positions.len() {
                break; // malformed / past end
            }
            let abs = positions[cursor];
            cursor += 1;

            if abs < start || abs >= end {
                continue;
            }
            // Index of `abs` within the window (present by construction).
            let widx = window.partition_point(|&p| p < abs) as isize;
            new_deltas.push((widx - prev_widx - 1) as usize);
            prev_widx = widx;

            // Clamp BOTH ends: a truncated/malformed ML (shorter than
            // deltas.len()*ncodes) must degrade to an empty run, never an
            // inverted slice range (which would panic).
            let ml_start = (k * ncodes).min(g.ml.len());
            let ml_end = (ml_start + ncodes).min(g.ml.len());
            new_ml.extend_from_slice(&g.ml[ml_start..ml_end]);
        }

        if !new_deltas.is_empty() {
            out.push(MmGroup {
                base: g.base,
                strand: g.strand,
                codes: g.codes.clone(),
                status: g.status,
                deltas: new_deltas,
                ml: new_ml,
            });
        }
    }

    Mods { groups: out }
}
```

Append Step 1's tests.

- [ ] **Step 4: Run and commit**

Run: `cargo test mods::reconstruct`
Expected: PASS (4 tests).

```bash
git add src/mods/reconstruct.rs
git commit -m "feat: MM/ML reconstruction over a trim interval (with minus-strand counting)"
```

---

### Task 3: Serialize (`src/mods/serialize.rs`)

**Files:**
- Create: `src/mods/serialize.rs`

**Interfaces:**
- Consumes: `super::{ModCode, Mods}`
- Produces: `pub fn serialize(mods: &Mods) -> (Vec<u8>, Vec<u8>)` — `(MM bytes, ML bytes)`. Empty-delta groups are skipped.

- [ ] **Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::parse;

    #[test]
    fn roundtrip_single_group() {
        let (mm, ml) = serialize(&parse(b"C+m?,5,12,0;", &[200, 10, 128]));
        assert_eq!(mm, b"C+m?,5,12,0;");
        assert_eq!(ml, vec![200, 10, 128]);
    }

    #[test]
    fn roundtrip_multi_group_and_chebi() {
        let input = b"C+mh,1,3;A+16061,2;".as_slice();
        let (mm, ml) = serialize(&parse(input, &[1, 2, 3, 4, 5]));
        assert_eq!(mm, input);
        assert_eq!(ml, vec![1, 2, 3, 4, 5]);
    }

    #[test]
    fn skips_empty_groups() {
        // A group with no deltas should not be emitted.
        let mut mods = parse(b"C+m,0;", &[7]);
        mods.groups.push(crate::mods::MmGroup {
            base: b'A', strand: b'+', codes: vec![crate::mods::ModCode::Char(b'a')],
            status: None, deltas: vec![], ml: vec![],
        });
        let (mm, _) = serialize(&mods);
        assert_eq!(mm, b"C+m,0;");
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test mods::serialize`
Expected: FAIL.

- [ ] **Step 3: Implement `src/mods/serialize.rs`**

```rust
use super::{ModCode, Mods};

pub fn serialize(mods: &Mods) -> (Vec<u8>, Vec<u8>) {
    let mut mm = Vec::new();
    let mut ml = Vec::new();

    for g in &mods.groups {
        if g.deltas.is_empty() {
            continue;
        }
        mm.push(g.base);
        mm.push(g.strand);
        for code in &g.codes {
            match code {
                ModCode::Char(c) => mm.push(*c),
                ModCode::Chebi(id) => mm.extend_from_slice(id.to_string().as_bytes()),
            }
        }
        if let Some(s) = g.status {
            mm.push(s);
        }
        for d in &g.deltas {
            mm.push(b',');
            mm.extend_from_slice(d.to_string().as_bytes());
        }
        mm.push(b';');
        ml.extend_from_slice(&g.ml);
    }

    (mm, ml)
}
```

Append Step 1's tests.

- [ ] **Step 4: Run and commit**

Run: `cargo test mods::`
Expected: PASS (whole codec).

```bash
git add src/mods/serialize.rs
git commit -m "feat: MM/ML serializer (round-trips parse; skips empty groups)"
```

---

### Task 4: BAM I/O + aligned-read refusal (`src/io/bam.rs`)

**Files:**
- Create: `src/io/bam.rs`
- Modify: `src/io/mod.rs` (add `pub mod bam;`)

**Interfaces:**
- Produces:
  - `pub fn reader(input: Option<&Path>) -> anyhow::Result<(sam::Header, Box<dyn Iterator<Item = anyhow::Result<RecordBuf>>>)>`
  - `pub fn writer(output: Option<&Path>, header: &sam::Header) -> anyhow::Result<bam::io::Writer<...>>` — returns a writer with the header already written and a `@PG` line appended.
  - `pub fn ensure_unaligned(rec: &RecordBuf) -> anyhow::Result<()>` — errors (naming the read) if the record is aligned.
- Where `RecordBuf = noodles_sam::alignment::RecordBuf`, `sam = noodles_sam`, `bam = noodles_bam`.

- [ ] **Step 1: Write failing test for `ensure_unaligned`**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;

    #[test]
    fn unmapped_ok_mapped_rejected() {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        assert!(ensure_unaligned(&rec).is_ok());

        *rec.flags_mut() = Flags::empty(); // mapped
        let err = ensure_unaligned(&rec).unwrap_err().to_string();
        assert!(err.contains("r1"));
        assert!(err.contains("aligned"));
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test io::bam`
Expected: FAIL.

- [ ] **Step 3: Implement `src/io/bam.rs`**

```rust
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::RecordBuf};

/// Error (naming the read) if the record is aligned. uBAM only in v1.
pub fn ensure_unaligned(rec: &RecordBuf) -> anyhow::Result<()> {
    if rec.flags().is_unmapped() {
        return Ok(());
    }
    let name = rec
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .unwrap_or_else(|| "<unnamed>".to_string());
    anyhow::bail!(
        "read {name} is aligned (mapped); chopping v1 supports unaligned BAM (uBAM) only"
    )
}

/// Open a BAM reader; return the header and an owning `RecordBuf` iterator.
pub fn reader(
    input: Option<&Path>,
) -> anyhow::Result<(sam::Header, Box<dyn Iterator<Item = anyhow::Result<RecordBuf>>>)> {
    let inner: Box<dyn io::Read> = match input {
        Some(p) => Box::new(File::open(p)?),
        None => Box::new(io::stdin()),
    };
    let mut r = bam::io::Reader::new(inner);
    let header = r.read_header()?;
    let header_for_iter = header.clone();
    let iter = RecordBufIter { reader: r, header: header_for_iter };
    Ok((header, Box::new(iter)))
}

struct RecordBufIter<R: io::Read> {
    reader: bam::io::Reader<R>,
    header: sam::Header,
}

impl<R: io::Read> Iterator for RecordBufIter<R> {
    type Item = anyhow::Result<RecordBuf>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = RecordBuf::default();
        match self.reader.read_record_buf(&self.header, &mut buf) {
            Ok(0) => None,
            Ok(_) => Some(Ok(buf)),
            Err(e) => Some(Err(e.into())),
        }
    }
}

/// A BAM writer with the (provenance-annotated) header already written.
pub fn writer(
    output: Option<&Path>,
    header: &sam::Header,
) -> anyhow::Result<bam::io::Writer<noodles_bgzf::io::Writer<Box<dyn Write>>>> {
    let inner: Box<dyn Write> = match output {
        Some(p) => Box::new(File::create(p)?),
        None => Box::new(io::stdout()),
    };
    let mut w = bam::io::Writer::new(inner);
    w.write_header(header)?;
    Ok(w)
}
```

Add `pub mod bam;` to `src/io/mod.rs`.

Note: `Flags::UNMAPPED` and `Flags::is_unmapped()` come from `noodles_sam::alignment::record::Flags`. `RecordBuf::name()` returns `Option<&BStr>`.

- [ ] **Step 4: Run and commit**

Run: `cargo test io::bam`
Expected: PASS.

```bash
git add src/io/bam.rs src/io/mod.rs
git commit -m "feat: uBAM reader/writer (RecordBuf) + aligned-read refusal"
```

---

### Task 5: BAM reconstruction + pipeline (`src/pipeline.rs`, `src/lib.rs`)

**Files:**
- Modify: `src/pipeline.rs` (add `run_bam` + `reconstruct_record`)
- Modify: `src/lib.rs` (`run` dispatches BAM→BAM)

**Interfaces:**
- Consumes: `filter::passes`, `trim::apply`, `mods::{parse, reconstruct, serialize}`, `io::bam`.
- Produces:
  - `pub fn pipeline::reconstruct_record(src: &RecordBuf, start: usize, end: usize, total: usize, idx: usize) -> RecordBuf`
  - `pub fn pipeline::run_bam(header: &sam::Header, records, writer, cfg) -> anyhow::Result<Stats>`

- [ ] **Step 1: Write failing test for `reconstruct_record`**

```rust
#[cfg(test)]
mod bam_tests {
    use super::*;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    fn ubam_with_mods(seq: &[u8], quals: Vec<u8>, mm: &[u8], ml: Vec<u8>) -> RecordBuf {
        let mut rec = RecordBuf::default();
        *rec.flags_mut() = Flags::UNMAPPED;
        *rec.name_mut() = Some(b"r1".into());
        *rec.sequence_mut() = seq.to_vec().into();
        *rec.quality_scores_mut() = quals.into();
        let data = rec.data_mut();
        data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm.to_vec().into()));
        data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(ml)));
        rec
    }

    #[test]
    fn slices_seq_qual_and_rebuilds_tags() {
        // seq = C C A C ; C+m modified at C occ 0 and 2 -> pos 0 and 3; ML [10,20].
        let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0,1;", vec![10, 20]);
        // keep window [2,4): seq "AC", the modified C at pos3 survives (occ within window idx0).
        let out = reconstruct_record(&src, 2, 4, 1, 0);

        assert_eq!(out.sequence().as_ref(), b"AC");
        assert_eq!(out.quality_scores().as_ref(), &[32, 33]);

        let mm = match out.data().get(&Tag::BASE_MODIFICATIONS) {
            Some(Value::String(s)) => s.to_vec(),
            _ => panic!("no MM"),
        };
        assert_eq!(mm, b"C+m,0;");
        let ml = match out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) => v.clone(),
            _ => panic!("no ML"),
        };
        assert_eq!(ml, vec![20]);
        // MN updated to the output length.
        let mn = match out.data().get(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH) {
            Some(Value::Int32(n)) => *n,
            _ => panic!("no MN"),
        };
        assert_eq!(mn, 2);
    }

    #[test]
    fn split_suffixes_name_and_drops_empty_mods() {
        let src = ubam_with_mods(b"CCAC", vec![30, 31, 32, 33], b"C+m,0;", vec![10]); // mod at pos0
        // segment [2,4) has no surviving C mod -> MM/ML removed entirely.
        let out = reconstruct_record(&src, 2, 4, 2, 1);
        assert_eq!(out.name().unwrap().as_ref(), b"r1_segment_2");
        assert!(out.data().get(&Tag::BASE_MODIFICATIONS).is_none());
        assert!(out.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES).is_none());
    }
}
```

- [ ] **Step 2: Run to verify failure**

Run: `cargo test pipeline::bam_tests`
Expected: FAIL.

- [ ] **Step 3: Implement `reconstruct_record` + `run_bam` in `src/pipeline.rs`**

```rust
use noodles_sam::{self as sam, alignment::RecordBuf};
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;

use crate::mods;

/// Build one output uBAM record for interval [start,end): slice SEQ/QUAL, rebuild
/// MM/ML/MN, suffix the name on splits. Non-mod aux tags are carried through
/// unchanged (they ride along in the cloned RecordBuf).
pub fn reconstruct_record(
    src: &RecordBuf,
    start: usize,
    end: usize,
    total: usize,
    idx: usize,
) -> RecordBuf {
    let mut out = src.clone();

    // Slice sequence + quality.
    let seq = src.sequence().as_ref().to_vec();
    let qual = src.quality_scores().as_ref().to_vec();
    *out.sequence_mut() = seq[start..end].to_vec().into();
    *out.quality_scores_mut() = qual[start..end].to_vec().into();

    // Name suffix on splits.
    if total > 1 {
        let base = src.name().map(|n| n.to_vec()).unwrap_or_default();
        let mut name = base;
        name.extend_from_slice(format!("_segment_{}", idx + 1).as_bytes());
        *out.name_mut() = Some(name.into());
    }

    // Rebuild MM/ML/MN when the source carried modification tags.
    let mm_raw = match src.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => Some(s.to_vec()),
        _ => None,
    };
    if let Some(mm_raw) = mm_raw {
        let ml_raw: Vec<u8> = match src.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES) {
            Some(Value::Array(Array::UInt8(v))) => v.clone(),
            _ => Vec::new(),
        };
        let parsed = mods::parse(&mm_raw, &ml_raw);
        let sliced = mods::reconstruct(&parsed, &seq, start, end);
        let (mm_new, ml_new) = mods::serialize(&sliced);

        let data = out.data_mut();
        if mm_new.is_empty() {
            // No mods survive in this segment: drop MM/ML AND MN (an MN with no
            // MM/ML is an orphaned, meaningless tag that confuses downstream tools).
            data.remove(&Tag::BASE_MODIFICATIONS);
            data.remove(&Tag::BASE_MODIFICATION_PROBABILITIES);
            data.remove(&Tag::BASE_MODIFICATION_SEQUENCE_LENGTH);
        } else {
            data.insert(Tag::BASE_MODIFICATIONS, Value::String(mm_new.into()));
            data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(ml_new)));
            // MN reflects the output segment length.
            data.insert(
                Tag::BASE_MODIFICATION_SEQUENCE_LENGTH,
                Value::Int32((end - start) as i32),
            );
        }
    }

    out
}

/// Single-threaded uBAM pipeline: refuse aligned reads, filter, trim, reconstruct.
pub fn run_bam<W>(
    header: &sam::Header,
    records: impl Iterator<Item = anyhow::Result<RecordBuf>>,
    writer: &mut noodles_bam::io::Writer<W>,
    cfg: &Config,
) -> anyhow::Result<Stats>
where
    W: std::io::Write,
{
    let mut stats = Stats::default();
    for rec in records {
        let rec = rec?;
        crate::io::bam::ensure_unaligned(&rec)?;
        stats.input_reads += 1;

        let seq = rec.sequence().as_ref().to_vec();
        let qual = rec.quality_scores().as_ref().to_vec();
        if !filter::passes(&seq, &qual, &cfg.filter) {
            continue;
        }
        let intervals = trim::apply(seq.len(), &qual, &cfg.trim, cfg.filter.min_length);
        let total = intervals.len();
        for (idx, (s, e)) in intervals.into_iter().enumerate() {
            let out = reconstruct_record(&rec, s, e, total, idx);
            writer.write_alignment_record(header, &out)?;
            stats.output_reads += 1;
        }
    }
    Ok(stats)
}
```

- [ ] **Step 4: Wire BAM dispatch in `src/lib.rs` `run()`**

Replace the `(Format::Bam, _) | (_, Format::Bam)` arm:

```rust
(Format::Bam, Format::Bam) => {
    let (header, records) = io::bam::reader(in_path)?;
    // Provenance: append our @PG line to a cloned header before writing.
    let out_header = provenance_header(header);
    let mut writer = io::bam::writer(cfg.io.output.as_deref(), &out_header)?;
    let stats = pipeline::run_bam(&out_header, records, &mut writer, &cfg)?;
    eprintln!("Kept {} reads out of {}", stats.output_reads, stats.input_reads);
    Ok(())
}
(Format::Bam, _) | (_, Format::Bam) => {
    anyhow::bail!("cross-format BAM<->FASTQ conversion is not supported in v1")
}
```

Add a small helper in `lib.rs`:

```rust
fn provenance_header(mut header: noodles_sam::Header) -> noodles_sam::Header {
    use noodles_sam::header::record::value::{map::Program, Map};
    let program = Map::<Program>::default();
    // Best-effort: insert a program record keyed "chopping"; ignore duplicates.
    let _ = header
        .programs_mut()
        .insert(b"chopping".into(), program);
    header
}
```

Note: if the exact `programs_mut().insert` signature differs in noodles 0.85, fall back to writing the header unchanged — the `@PG` line is cosmetic and must not block record output. Confirm the `Programs::insert` signature against the cargo source when implementing (`header/programs.rs`).

- [ ] **Step 5: Run and commit**

Run: `cargo test pipeline::bam_tests && cargo clippy -- -D warnings`
Expected: PASS; clean.

```bash
git add src/pipeline.rs src/lib.rs
git commit -m "feat: uBAM pipeline — slice + MM/ML/MN reconstruction, run() dispatch"
```

---

### Task 6: htslib oracle integration test (`tests/bam_mods_oracle.rs`)

**Files:**
- Create: `tests/bam_mods_oracle.rs`
- Create: a tiny builder that writes `test-data/mods_small.ubam` (a Rust test helper, so no external tooling needed).

**Interfaces:**
- Consumes: `chopping` public API (`run` or `pipeline::run_bam`) + `rust-htslib::bam::Read` / `basemods_iter`.

The oracle asserts **decode equivalence**: decode the mods of our trimmed output with htslib, and assert the per-position `(canonical, code, strand, ml)` set equals the original mods decoded by htslib, filtered to the trim interval and offset by `start`. This is independent of MM's multiple valid encodings.

- [ ] **Step 1: Write the fixture builder + oracle test**

```rust
use std::path::Path;

use noodles_bam as bam;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam as sam;

use rust_htslib::bam::{self as hts, Read as _};

/// Build a one-read uBAM: seq with C's, C+m mods at chosen positions, ML bytes.
fn write_fixture(path: &Path) {
    let header = sam::Header::builder()
        .set_header(Default::default())
        .build();
    let mut w = bam::io::Writer::new(std::fs::File::create(path).unwrap());
    w.write_header(&header).unwrap();

    let mut rec = RecordBuf::default();
    *rec.flags_mut() = Flags::UNMAPPED;
    *rec.name_mut() = Some(b"read1".into());
    // 12 bases, C at indices 1,4,7,10.
    *rec.sequence_mut() = b"ACGGCGGCGGCG".to_vec().into();
    *rec.quality_scores_mut() = vec![40u8; 12].into();
    let data = rec.data_mut();
    // modify C occurrences 0,2,3 -> deltas 0,1,0 ; ML 3 bytes.
    data.insert(Tag::BASE_MODIFICATIONS, Value::String(b"C+m,0,1,0;".to_vec().into()));
    data.insert(Tag::BASE_MODIFICATION_PROBABILITIES, Value::Array(Array::UInt8(vec![250, 5, 200])));
    data.insert(Tag::BASE_MODIFICATION_SEQUENCE_LENGTH, Value::Int32(12));
    w.write_alignment_record(&header, &rec).unwrap();
}

/// Decode (0-based read pos, canonical, modified, strand, qual) with htslib.
fn hts_mods(path: &Path) -> Vec<(usize, char, char, i32, i32)> {
    let mut reader = hts::Reader::from_path(path).unwrap();
    let mut out = Vec::new();
    for rec in reader.records() {
        let rec = rec.unwrap();
        if let Ok(iter) = rec.basemods_iter() {
            for (pos, m) in iter.flatten() {
                out.push((
                    pos as usize,
                    m.canonical_base as u8 as char,
                    m.modified_base as u8 as char,
                    m.strand,
                    m.qual,
                ));
            }
        }
    }
    out
}

#[test]
fn trimmed_output_mods_match_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.ubam");
    let output = dir.path().join("out.ubam");
    write_fixture(&input);

    // Trim the first 3 bases (head-crop 3) via the library run().
    let cfg = chopping::cli::config_for_test(&input, &output, 3, 0);
    chopping::run(cfg).unwrap();

    let original = hts_mods(&input);
    let start = 3usize;
    let expected: Vec<_> = original
        .iter()
        .filter(|(pos, ..)| *pos >= start)
        .map(|&(pos, cb, mb, st, q)| (pos - start, cb, mb, st, q))
        .collect();

    let got = hts_mods(&output);
    let mut a = expected.clone();
    let mut b = got.clone();
    a.sort();
    b.sort();
    assert_eq!(a, b, "trimmed mod set must equal original filtered to [3, len) offset by 3");
}
```

- [ ] **Step 2: Add a test-only config constructor**

Because `cli::parse()` reads real argv, add a small helper in `src/cli.rs` guarded for tests and integration use:

```rust
/// Build a Config directly (used by integration tests). head/tail are fixed crops.
pub fn config_for_test(
    input: &std::path::Path,
    output: &std::path::Path,
    head_crop: usize,
    tail_crop: usize,
) -> Config {
    Config {
        io: IoConfig {
            input: Some(input.to_path_buf()),
            output: Some(output.to_path_buf()),
            in_format: Some(Format::Bam),
            out_format: Some(Format::Bam),
        },
        filter: FilterConfig {
            min_length: 1, max_length: usize::MAX, min_qual: 0.0, max_qual: 1000.0,
            min_gc: None, max_gc: None, qual_mode: QualMode::Mean,
        },
        trim: TrimPlan { head: head_crop, tail: tail_crop, quality: None },
        threads: 1,
    }
}
```

Add `tempfile = "3"` to `[dev-dependencies]` if not already present from Plan 1.

- [ ] **Step 3: Run the oracle test**

Run: `cargo test --test bam_mods_oracle`
Expected: PASS. (First build compiles bundled htslib — needs a C toolchain; slow once.)

- [ ] **Step 4: Commit**

```bash
git add tests/bam_mods_oracle.rs src/cli.rs Cargo.toml Cargo.lock
git commit -m "test: htslib decode-equivalence oracle for trimmed MM/ML"
```

---

### Task 7: Real-data oracle sweep + docs

**Files:**
- Modify: `tests/bam_mods_oracle.rs` (add an ignored, data-gated sweep over a real HG002 uBAM subset)
- Create: `README.md` (usage, formats, uBAM mod-recompute note, aligned-BAM limitation)

**Interfaces:** none new.

- [ ] **Step 1: Add a `#[ignore]` real-data test**

```rust
// Runs only when a real uBAM is provided, e.g.:
//   CHOPPING_UBAM=/path/hg002.subset.ubam cargo test --test bam_mods_oracle -- --ignored
#[test]
#[ignore]
fn real_ubam_oracle_sweep() {
    let Some(path) = std::env::var_os("CHOPPING_UBAM") else { return };
    let input = std::path::PathBuf::from(path);
    let dir = tempfile::tempdir().unwrap();
    let output = dir.path().join("out.ubam");

    let cfg = chopping::cli::config_for_test(&input, &output, 10, 10);
    chopping::run(cfg).unwrap();

    // For every read, our output mods (htslib-decoded) must equal the original
    // decoded mods filtered to [10, len-10) offset by 10. Compare per read name.
    let orig = hts_mods_by_read(&input);
    let out = hts_mods_by_read(&output);
    for (name, got) in &out {
        let base = name.strip_suffix("").unwrap_or(name); // no split with fixed crop
        let expected = orig.get(base).cloned().unwrap_or_default();
        // (expected filtering/offset applied inside hts_mods_by_read's caller — see helper)
        assert_eq!(sorted(got), sorted(&filter_offset(&expected, 10, len_of(&out, name))));
    }
}
```

Implement the small `hts_mods_by_read`, `sorted`, `filter_offset`, `len_of` helpers alongside (decode per read name into a map). Keep the assertion identical in spirit to Task 6.

- [ ] **Step 2: Write `README.md`**

Cover: what chopping is; install (`cargo build --release`); FASTQ and uBAM usage examples with the real flags (`-q/-Q`, `-l/-L`, `-H/-T`, `--trim-qual`, `--best-segment`, `--split-qual`/`--split-window`, `-m/--qual-mode`); the guarantee that **uBAM MM/ML/MN are recomputed on trim/split**; and the explicit v1 limitations (aligned BAM refused; no contamination filter; no cross-format).

- [ ] **Step 3: Run full suite + clippy, then commit**

Run: `cargo test && cargo clippy -- -D warnings`
Expected: PASS (the real-data sweep stays ignored without `CHOPPING_UBAM`).

```bash
git add tests/bam_mods_oracle.rs README.md
git commit -m "test: real-data oracle sweep (gated) + README"
```

---

## Self-Review

**Spec coverage (uBAM scope):**
- uBAM → uBAM, same filter/trim as FASTQ → Task 5 reuses `filter`/`trim`. ✅
- MM/ML recomputed + MN updated on every output read → Tasks 1–3 (codec), 5 (`reconstruct_record`). ✅
- Bypass noodles typed parser; use raw MM/ML → Task 1 parses raw; Task 5 reads `Value::String`/`Array::UInt8`. ✅
- Aligned-read refusal (unmapped flag clear → hard error naming read) → Task 4 `ensure_unaligned`, enforced in Task 5. ✅
- `@PG` provenance on output header → Task 5 `provenance_header` (best-effort, non-blocking). ✅
- Splits produce `_segment_N` names, empty-in-window mods dropped → Task 5 tests. ✅
- htslib **decode-equivalence** oracle (dev-dependency, test-only) on synthetic + real HG002 subset → Tasks 6, 7. ✅
- Cross-format BAM↔FASTQ explicitly rejected in v1 → Task 5 dispatch arm. ✅

**Placeholder scan:** Codec Tasks 1–3 are fully concrete. Two implementation-time confirmations are flagged, not hidden: (a) `Programs::insert` exact signature in noodles 0.85 (Task 5, with a non-blocking fallback), and (b) `sam::Header::builder().set_header(..)` for the fixture (Task 6 — if the exact builder differs, `sam::Header::default()` suffices, since the reader only needs a valid empty header). The Task 7 real-data helpers (`hts_mods_by_read`, etc.) are described with their exact contract; implement them as thin decode-into-map functions.

**Type consistency:** `Mods{groups}`, `MmGroup{base,strand,codes,status,deltas,ml}`, `ModCode`, `parse(mm,ml)`, `reconstruct(&Mods,seq,start,end)`, `serialize(&Mods)->(mm,ml)` are used identically in Tasks 1–3, 5, 6. `reconstruct_record(src,start,end,total,idx)` and `run_bam(header,records,writer,cfg)` match their call sites in `lib.rs`. `Tag::BASE_MODIFICATIONS` / `BASE_MODIFICATION_PROBABILITIES` / `BASE_MODIFICATION_SEQUENCE_LENGTH`, `Value::String` / `Value::Array(Array::UInt8)` / `Value::Int32` are the verified noodles 0.85 names.

**One risk noted for execution:** minus-strand mod handling (Task 2 `counting_base` complement rule) is the least-certain semantic. It is covered by a unit test and, decisively, by the htslib oracle sweep (Task 7) on real HG002 reads — if a real `-` group diverges, the oracle fails loudly and the rule gets fixed there rather than shipping silently wrong.
```
