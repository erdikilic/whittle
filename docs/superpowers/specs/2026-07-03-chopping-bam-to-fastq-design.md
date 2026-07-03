# chopping — BAM → FASTQ/.gz conversion (design)

- **Date:** 2026-07-03
- **Status:** Approved
- **Repo:** https://github.com/erdikilic/chopping (private dev repo)
- **Builds on:** [2026-07-03-chopping-tag-aware-trimmer-design.md](2026-07-03-chopping-tag-aware-trimmer-design.md)
  (v1: FASTQ↔FASTQ and uBAM→uBAM), which explicitly deferred cross-format conversion.

## Goal

Add the one missing conversion direction: **unaligned BAM → FASTQ (plain or gzip)**,
for both a single file and a folder of BAMs. When a read is trimmed or split, its
MM/ML/MN modification tags are reconstructed exactly as they already are on the
BAM→BAM path, then written into the FASTQ header alongside the read's other aux
tags — the `samtools fastq -T` / `samtools import -T` convention — so the mods
survive the format change and can be round-tripped back into BAM.

## Scope

**In:**
- `BAM → FASTQ` and `BAM → FASTQ.gz`, single file and stdin/stdout.
- `BAM folder → FASTQ/.gz` (merge), mirroring the existing folder-merge path.
- A `--fastq-tags` selector controlling which aux tags land in the FASTQ header.

**Out (unchanged from v1):**
- `FASTQ → BAM` (any direction): we cannot invent a header or aux tags from FASTQ,
  so this stays a hard error.
- Aligned BAM input: still refused with the existing aligned-read error.
- BAM-path parallelism: the BAM read/trim loop stays single-threaded (as `run_bam`
  today). Gzip *output* is still parallel via `gzp`.

## Support matrix (after this change)

| input ↓ / output → | FASTQ | FASTQ.gz | BAM |
|--------------------|:-----:|:--------:|:---:|
| FASTQ / FASTQ.gz   |  ✅   |    ✅    | ❌ (bail) |
| uBAM               | ✅ new | ✅ new  | ✅ |

Format is auto-detected from extensions / `--out-format`, exactly as today. The
default when `-o` is omitted still **mirrors** the input (BAM→BAM), so existing
invocations are unaffected; BAM→FASTQ happens only when the output is explicitly
FASTQ (via `-o *.fastq`/`*.fastq.gz`/`*.gz` or `--out-format fastq|fastq-gz`).

## CLI — one new flag

```
--fastq-tags <all|none|LIST>   aux tags to carry into FASTQ headers on BAM→FASTQ
                               [default: all]
```

- `all` — carry every aux tag from the source record.
- `none` — carry no tags; emit plain FASTQ.
- `LIST` — a comma-separated list of 2-character SAM tags, e.g. `MM,ML,RG`; carry
  only those.

**Validation:** each listed token must be exactly two characters (SAM tag shape);
anything else is a hard error naming the offending token. `all`/`none` are
reserved keywords and are matched case-sensitively (lowercase).

**Applicability:** the flag is consulted **only** on the BAM→FASTQ path. If it is
set on a path where it cannot apply (FASTQ input, or BAM→BAM), chopping prints a
one-line note to stderr (`--fastq-tags is ignored for <in>→<out>`) and proceeds —
not a hard error. (Precedent: `--in-format` is a no-op in folder mode.)

## Tag handling

MM/ML/MN are **position-dependent** and are reconstructed (sliced to the surviving
interval) by the existing `mods` codec. Every other aux tag is **position-agnostic
from this tool's perspective** and is copied **verbatim**.

> **Caveat (documented, by design):** a handful of non-MM/ML aux tags are in fact
> position-dependent (e.g. dorado's `mv` move table, `ts`/`ns` signal-trim counts).
> Under `--fastq-tags all` these are copied verbatim and will be **stale** after
> trimming — identical to `samtools fastq -T`, which likewise copies without
> adjusting. Only MM/ML/MN are trim-aware. Users who care can exclude such tags via
> an explicit `--fastq-tags` list.

### Per-interval algorithm

For each surviving interval `[start, end)` of a kept read:

1. Compute reconstructed mods once: `reconstruct_mods(src, seq, start, end) ->
   Option<(mm, ml)>` (returns `None` when the source has no `MM`, or when every
   modified position was trimmed away).
2. Walk the source record's aux fields **in source order**. For each field whose
   tag is carried (per the rules below) **and** is not one of `MM`/`ML`/`MN`,
   append its SAM-text form `XX:T:val`.
3. If mods are carried **and** step 1 returned `Some((mm, ml))`, append the
   trim-aware block `MM:Z:<mm>`, `ML:B:C,<b0>,<b1>,…`, `MN:i:<end-start>`.

Header order is irrelevant to correctness — `samtools import` parses aux by tag
name — so appending the reconstructed mod block last is fine and keeps the walk
deterministic.

**"Carried" rules** for a given `--fastq-tags` value:
- `all`: every tag is carried.
- `none`: no tag is carried (header is just the read id; result is plain FASTQ).
- `LIST`: a non-mod tag is carried iff it appears in the list. The **mod block**
  (MM/ML/MN as a unit) is carried iff the list contains `MM` or `ML`. `MN` is
  bound to the mod block and is never independently selectable — listing `MN`
  alone carries nothing.

### Header format

```
@<qname>[_segment_N]\t<tag>\t<tag>…
<sliced seq>
+
<sliced qual>
```

- `_segment_N` suffix on splits, matching the existing FASTQ writer (`total > 1`).
- Tags are TAB-separated and follow the (possibly suffixed) read id.
- A read that ends up with no carried tags is written as an ordinary plain FASTQ
  record (no trailing TAB).

### SAM-text aux serialization (`format_aux_field`)

A small self-contained serializer maps each noodles `RecordBuf` data `Value` to
its SAM textual form `XX:T:VALUE`. We write our own rather than depend on noodles'
internal SAM record writer — consistent with the project already owning its MM/ML
codec.

| noodles `Value`            | SAM type | example        |
|----------------------------|:--------:|----------------|
| `Character(c)`             | `A`      | `XX:A:c`       |
| `Int8/Int16/Int32/UInt8/UInt16/UInt32` | `i` | `XX:i:-3` |
| `Float(f)`                 | `f`      | `XX:f:0.5`     |
| `String(s)`                | `Z`      | `XX:Z:abc`     |
| `Hex(h)`                   | `H`      | `XX:H:1AE3`    |
| `Array(Int8)`              | `B:c`    | `XX:B:c,-1,2`  |
| `Array(UInt8)`             | `B:C`    | `XX:B:C,1,2`   |
| `Array(Int16)`             | `B:s`    | `XX:B:s,…`     |
| `Array(UInt16)`            | `B:S`    | `XX:B:S,…`     |
| `Array(Int32)`             | `B:i`    | `XX:B:i,…`     |
| `Array(UInt32)`            | `B:I`    | `XX:B:I,…`     |
| `Array(Float)`             | `B:f`    | `XX:B:f,…`     |

Integers all serialize with SAM type code `i` on output (SAM's textual integer
type), regardless of the source width. The reconstructed `MM`/`ML`/`MN` block is
produced directly by `format_mods_aux(mm, ml, mn)` (not routed through
`format_aux_field`), since MM is a `Z` string, ML a `B:C` array, and MN an `i`.

## Architecture

Reuses the existing format-neutral pipeline. New/changed units:

```
src/
  cli.rs        + --fastq-tags arg + parse/validate -> FastqTags
  config.rs     + FastqTags { All, None, Only(BTreeSet<[u8;2]>) }; field on Config
  pipeline.rs   + reconstruct_mods() extracted from reconstruct_record()
                + run_bam_to_fastq()
  io/fastq.rs   + write_segment_tagged(), format_aux_field(), format_mods_aux()
  lib.rs        run()/run_folder(): new (Bam, Fastq|FastqGz) dispatch arm
```

### `pipeline.rs`

- **`reconstruct_mods(src: &RecordBuf, seq: &[u8], start: usize, end: usize) ->
  Option<(Vec<u8>, Vec<u8>)>`** — the MM/ML parse→reconstruct→serialize block
  currently inlined in `reconstruct_record` is lifted into this function. Returns
  `None` when there is no `MM` tag or the sliced result is empty. `reconstruct_record`
  is refactored to call it (behavior identical; its existing tests must still pass).

- **`run_bam_to_fastq<W: Write>(records, writer, cfg) -> anyhow::Result<Stats>`** —
  mirrors `run_bam`'s per-record guards (`ensure_unaligned`, SEQ/QUAL-length check),
  runs `filter::passes` + `trim::apply`, and for each interval calls
  `write_segment_tagged`. Single-threaded. The header is not needed for FASTQ
  output. Parse/IO errors propagate via `?`, as in `run_bam`.

### `io/fastq.rs`

- **`write_segment_tagged(w, name, seq, phred, total, idx, tags: &[u8])`** — like
  `write_segment` but inserts the pre-built `tags` byte block (each element already
  prefixed with a TAB, or empty) between the header id and the newline. The shared
  id-building logic (`@`, `_segment_N` suffix, optional description) is factored so
  `write_segment` and `write_segment_tagged` cannot drift.
- **`format_aux_field(tag: [u8;2], value: &Value) -> Vec<u8>`** — the table above.
- **`format_mods_aux(mm: &[u8], ml: &[u8], mn: usize) -> Vec<u8>`** — emits
  `\tMM:Z:<mm>\tML:B:C,<ml…>\tMN:i:<mn>`.

The assembly of the full `tags` block for one interval (walking source fields,
applying the carry rules, appending the reconstructed mod block) lives in
`run_bam_to_fastq` (or a private helper next to it in `pipeline.rs`), so
`io/fastq.rs` stays a pure formatter with no knowledge of `--fastq-tags`.

### `lib.rs` dispatch

In both `run()` and `run_folder()`:
- `(Format::Bam, Format::Fastq)` and `(Format::Bam, Format::FastqGz)` → open the
  BAM reader (`io::bam::reader` / `io::dir::bam_reader`), build the FASTQ writer
  (`fastq_writer(cfg, out_fmt)`), call `run_bam_to_fastq`, then `writer.finish()?`
  (the same gz-finalize seam the FASTQ paths use).
- `(Format::Fastq | Format::FastqGz, Format::Bam)` → keep the existing hard bail
  (`FASTQ→BAM not supported`).

## Config plumbing

`Config` gains `fastq_tags: FastqTags`. `cli.rs` parses the `--fastq-tags` string:
`all` → `FastqTags::All`, `none` → `FastqTags::None`, otherwise split on `,`,
validate each token is 2 chars, collect into `FastqTags::Only`. All existing
`Config { … }` literals in tests must be updated to set the new field (default
`FastqTags::All`).

## Error handling

- Malformed `--fastq-tags` token → hard error at CLI parse time.
- Aligned BAM record → existing aligned-read error (unchanged).
- SEQ/QUAL length mismatch → existing guard error (unchanged), reused verbatim in
  `run_bam_to_fastq`.
- FASTQ→BAM → existing hard bail (unchanged).

## Testing

**Unit:**
- `format_aux_field`: one case per `Value` variant and per array subtype (`c/C/s/S/i/I/f`).
- `format_mods_aux`: `MM:Z`, `ML:B:C`, `MN:i` assembly, including empty `ml`.
- `--fastq-tags` parsing: `all`, `none`, `MM,ML,RG`, and rejection of a bad token
  (`ABC`, empty).
- `write_segment_tagged`: header id + `_segment_N` + TAB-joined tags; empty-tags
  falls back to a plain record byte-identical to `write_segment`.
- `run_bam_to_fastq` over synthetic uBAM: (a) read with mods → header carries
  reconstructed `MM/ML/MN` sliced to the window; (b) split → `_segment_N` + per-
  segment mods; (c) `--fastq-tags none` → plain FASTQ; (d) segment with all mods
  trimmed → plain FASTQ; (e) a non-mod tag (`RG`) carried verbatim under `all` and
  dropped under `MM,ML`.

**Cross-check:** for the same record+interval, the `MM/ML/MN` rendered into the
FASTQ header must equal the values `reconstruct_record` writes into the BAM record
(both go through `reconstruct_mods`) — guards against the two paths drifting.

**Correctness (hermetic, no external tools):** the cross-check above is the
load-bearing guarantee — the BAM→FASTQ `MM/ML/MN` bytes must equal the BAM→BAM
output's, which is already htslib-oracle-verified (`tests/bam_mods_oracle.rs`,
including a 57k-read real-data sweep). Because both paths call the same
`reconstruct_mods`, equality transitively proves the FASTQ-header mods are
oracle-correct. No `samtools` dependency in the suite. (An optional, gated
`samtools import -T '*'` interop check may be added later purely to confirm the
real tool accepts our header format — never load-bearing.)

**CLI (`assert_cmd`):** `-i tests/…/mods.bam -o out.fastq` and `-o out.fastq.gz`
produce the expected header tags; `--out-format fastq` on a `.bam` input works;
`--fastq-tags` on a FASTQ input prints the ignored-note and still succeeds.

## Dependencies

No new dependencies. Reuses `noodles-bam`/`noodles-sam` (read + `Value` types),
`gzp` (parallel gz output), the existing `mods` codec, and `rust-htslib` (dev-dep,
oracle only).
