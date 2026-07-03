# chopping — folder-merge input (design)

- **Date:** 2026-07-03
- **Status:** Approved
- **Builds on:** v1 (Plan 1 FASTQ + Plan 2 uBAM), merged at `main` @ cec929f.

## Goal

Let `-i/--input` accept a **directory** as well as a single file. When given a
directory, chopping merges all read files in that directory into one trimmed
output — the ONT barcode-folder workflow (a `barcode03/` folder of `.fastq.gz`
chunks, or a folder of dorado `.bam` files, trimmed into one clean output).

## Scope

**In:**
- `-i <path>`: if `<path>` is a directory, enter folder-merge mode (auto-detected
  via `is_dir()` — no new flag). Single file and stdin behave exactly as today.
- **Single flat folder only** (non-recursive): the directory's immediate children.
- Merge → **one** output (`-o <file>` or stdout).

**Deferred (not this feature):** recursion; batch-per-subfolder
(`fastq_pass/` → one output per `barcode*/`); per-file mirroring; cross-`@RG`
header merging.

## Behavior

### File selection
Immediate children of the directory that are read files by extension:
FASTQ-family (`.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`) or `.bam`. Everything else
(`sequencing_summary.txt`, indexes, subdirectories, …) is ignored. Selected files
are processed in **sorted filename order** for deterministic output.

### Homogeneity
The folder must resolve to one format family:
- **All FASTQ-family** (plain and `.gz` may coexist; each file is gz-decoded by its
  own extension) → one FASTQ output.
- **All `.bam`** → one BAM output.

A folder mixing `.bam` with `.fastq*` is a hard error (they cannot merge into a
single output). A folder with **no** read files is a hard error. Both name the
directory and explain.

### FASTQ merge
Chain the per-file record iterators into one stream (each file opened and
gz-decoded per its own extension), fed to the existing `pipeline::run_fastq`.
Trimming/filtering/splitting and `--threads` behave exactly as for a single file.
Output: `-o` file/stdout, gz when the output extension is `.gz`.

### BAM merge
Use the **first** (sorted) file's header — `samtools cat` semantics; a single
barcode's uBAMs share their `@HD`/`@RG`/`@PG`. Append our `@PG` provenance
(`provenance_header`, already dangling-`PP`-safe). Stream every file's records
under that header into `pipeline::run_bam` → one merged BAM. Aligned-read refusal
(`ensure_unaligned`) still fires per record across all files. The explicit
`writer.try_finish()?` still runs once at the end.

### Output
Unchanged: `-o <file>` or stdout; format from the `-o` extension, else the
folder's format family (FASTQ-family folder → FASTQ; BAM folder → BAM).

## Code shape

New module `src/io/dir.rs`:
- `classify(dir: &Path) -> anyhow::Result<(Family, Vec<PathBuf>)>` — enumerate
  immediate children, keep read files, sort, determine `Family` (`Fastq` | `Bam`),
  error on mixed/empty. `enum Family`.
- `fastq_records(paths: &[PathBuf]) -> anyhow::Result<Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send>>`
  — chain each file's `io::fastq::reader(path, is_gz(path))`.
- `bam_reader(paths: &[PathBuf]) -> anyhow::Result<(sam::Header, Box<dyn Iterator<Item = anyhow::Result<RecordBuf>>>)>`
  — header from `paths[0]`; chain each file's `io::bam::reader(path).1`.

`run()` in `src/lib.rs`: at the top, if `in_path` is a directory, `classify` it and
route to a folder-merge branch that reuses the existing writer construction +
`run_fastq`/`run_bam`. The single-file/stdin path below is left untouched. (Factor
the writer-open + `run_bam` finish + "Kept …" reporting so both the single-file and
folder branches share them rather than duplicating.)

## Errors
Fail fast with context: mixed-format folder, empty folder, unreadable directory.
A malformed/aligned record mid-merge surfaces as today (nonzero exit).

## Testing
- `dir::classify` unit tests: FASTQ-only folder, BAM-only folder, gz+plain FASTQ
  mix (→ Fastq), mixed fastq+bam (→ error), empty/no-read-files (→ error),
  non-read files and subdirs ignored, sorted order.
- Integration (temp dir): two small FASTQ files → one merged trimmed output with
  the expected record set.
- Integration (temp dir): two small uBAM files → one merged BAM that is readable,
  has the expected record count, and carries the `@PG chopping` line.
