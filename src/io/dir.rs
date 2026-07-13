use std::path::{Path, PathBuf};

use anyhow::Context;
use noodles_sam::{self as sam};

use crate::io::{Format, from_extension};
use crate::record::ReadRecord;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    Fastq,
    Bam,
}

/// Classify a directory's immediate children into a single read-file family and
/// a sorted path list. Non-read files and subdirectories are ignored.
/// Errors if the folder mixes FASTQ and BAM, or contains no read files.
///
/// If `output` names a read file *inside* the directory, this is a hard error: we
/// can't tell a real input from a stale prior output, and either way merging the
/// rest and overwriting it silently loses data (this covers both `-o` pointing at
/// a genuine input file and a rerun re-ingesting its own output). The merged
/// output must be written to a path outside the input directory.
pub fn classify(dir: &Path, output: Option<&Path>) -> anyhow::Result<(Family, Vec<PathBuf>)> {
    let mut fastq = Vec::new();
    let mut bam = Vec::new();

    let entries =
        std::fs::read_dir(dir).with_context(|| format!("reading directory {}", dir.display()))?;
    for entry in entries {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        let format = from_extension(&path);
        if matches!(format, Some(Format::Fastq | Format::FastqGz | Format::Bam))
            && output.is_some_and(|o| crate::same_path(&path, o))
        {
            anyhow::bail!(
                "output {} is a read file inside the input directory {}; refusing to \
                 overwrite it (it may be one of the inputs) — write the merged output \
                 to a path outside {}",
                path.display(),
                dir.display(),
                dir.display()
            );
        }
        match format {
            Some(Format::Fastq | Format::FastqGz) => fastq.push(path),
            Some(Format::Bam) => bam.push(path),
            None => {}, // ignore non-read files
        }
    }

    let (family, mut paths) = match (fastq.is_empty(), bam.is_empty()) {
        (false, true) => (Family::Fastq, fastq),
        (true, false) => (Family::Bam, bam),
        (true, true) => anyhow::bail!(
            "no FASTQ or BAM read files found in directory {}",
            dir.display()
        ),
        (false, false) => anyhow::bail!(
            "directory {} mixes FASTQ ({}) and BAM ({}) files; a folder must be \
             one format to merge into a single output",
            dir.display(),
            fastq.len(),
            bam.len()
        ),
    };
    paths.sort();
    Ok((family, paths))
}

/// One chained record stream over all FASTQ-family files, opened lazily as the
/// chain reaches each file (avoids exhausting file descriptors on large folders).
/// A file-open error surfaces as an `Err` item rather than aborting construction.
pub fn fastq_records(
    paths: &[PathBuf],
) -> Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send> {
    let paths = paths.to_vec();
    Box::new(paths.into_iter().flat_map(
        |p| -> Box<dyn Iterator<Item = anyhow::Result<ReadRecord>> + Send> {
            let gz = matches!(from_extension(&p), Some(Format::FastqGz));
            match crate::io::fastq::reader(Some(&p), gz) {
                Ok(reader) => reader,
                Err(e) => Box::new(std::iter::once(Err(e))),
            }
        },
    ))
}

/// A boxed, owning iterator over lazy raw BAM records (or per-record errors).
/// Named to satisfy `clippy::type_complexity` on the `bam_reader` signature below.
/// `Send` so it can feed `workflow::run_bam`'s parallel path.
type BamRecordIter = crate::io::bam::RawRecordIter;

/// The first file's header plus one chained record stream over all BAM files
/// (each file's own header is read and discarded past the first; records
/// stream as lazy raw records under the first header — `samtools cat` semantics
/// for homogeneous uBAM).
///
/// `workers` is the MT-bgzf decode worker count (same knob as the single-file
/// `io::bam::reader`), passed through to every per-file reader. This is safe
/// from oversubscription despite chaining N files: the first file's reader is
/// opened eagerly right here, but each `rest` file's reader is opened lazily,
/// one at a time, only once `flat_map` actually reaches it — and `Iterator`'s
/// `Chain`/`FlatMap` drop the exhausted inner iterator before advancing to the
/// next, closing its MT reader (and joining its worker threads) first. So at
/// most one file's `workers` bgzf threads are ever live at once, never N×.
///
/// Returns an `Err` if `paths` is empty rather than panicking. Each file is
/// opened exactly once: the first file's record iterator (obtained alongside
/// its header) is reused via `chain` rather than reopening that file.
pub fn bam_reader(
    paths: &[PathBuf],
    workers: usize,
) -> anyhow::Result<(sam::Header, BamRecordIter)> {
    let (first, rest) = paths
        .split_first()
        .ok_or_else(|| anyhow::anyhow!("bam_reader called with no BAM files"))?;
    let (header, first_records) = crate::io::bam::reader(Some(first), workers)?;
    let rest = rest.to_vec();
    let rest_records = rest.into_iter().flat_map(move |p| -> BamRecordIter {
        match crate::io::bam::reader(Some(&p), workers) {
            Ok((_hdr, recs)) => recs,
            Err(e) => Box::new(std::iter::once(Err(e))),
        }
    });
    let records: BamRecordIter = Box::new(first_records.chain(rest_records));
    Ok((header, records))
}

/// The set of `@RG` ids declared in a BAM file's header (empty if none / unreadable).
fn bam_read_group_ids(path: &Path) -> Option<std::collections::BTreeSet<Vec<u8>>> {
    let (header, _records) = crate::io::bam::reader(Some(path), 1).ok()?;
    Some(header.read_groups().keys().map(|k| k.to_vec()).collect())
}

/// Name the first BAM whose `@RG` id set differs from the first file's, if any.
/// Best-effort: files that fail to open are skipped. Used to warn that folder
/// merge keeps only the first header.
fn first_rg_mismatch(paths: &[PathBuf]) -> Option<(PathBuf, PathBuf)> {
    let mut iter = paths.iter();
    let first = iter.next()?;
    let first_rgs = bam_read_group_ids(first)?;
    for p in iter {
        if let Some(rgs) = bam_read_group_ids(p)
            && rgs != first_rgs
        {
            return Some((first.clone(), p.clone()));
        }
    }
    None
}

/// Warn (once) if the folder's BAM files don't all share the same `@RG` set.
/// Folder merge writes only the first file's header (samtools-cat semantics), so
/// records from later files that carry `RG` tags for read groups absent from that
/// header would reference `@RG` lines missing from the output. Homogeneous
/// single-run uBAM — the intended input — shares one header and stays silent.
/// This is a header-only pre-pass (cheap for uBAM); a future "better method"
/// would union the `@RG`/`@PG` records into the output header instead.
pub fn warn_on_bam_header_mismatch(paths: &[PathBuf]) {
    if let Some((first, offender)) = first_rg_mismatch(paths) {
        tracing::warn!(
            "warning: the folder's BAM files have different @RG sets ({} vs {}); \
             only the first file's header is written, so records from other files \
             may reference read groups missing from the merged output header",
            first.display(),
            offender.display()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn touch(dir: &std::path::Path, name: &str) {
        std::fs::write(dir.join(name), b"x").unwrap();
    }

    #[test]
    fn fastq_only_folder_sorted() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "b_1.fastq.gz");
        touch(d.path(), "a_0.fastq");
        touch(d.path(), "sequencing_summary.txt"); // ignored
        std::fs::create_dir(d.path().join("subdir")).unwrap(); // ignored
        let (fam, paths) = classify(d.path(), None).unwrap();
        assert_eq!(fam, Family::Fastq);
        let names: Vec<_> = paths
            .iter()
            .map(|p| p.file_name().unwrap().to_str().unwrap())
            .collect();
        assert_eq!(names, vec!["a_0.fastq", "b_1.fastq.gz"]); // sorted, .txt/subdir excluded
    }

    #[test]
    fn bam_only_folder() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "x.bam");
        touch(d.path(), "y.bam");
        let (fam, paths) = classify(d.path(), None).unwrap();
        assert_eq!(fam, Family::Bam);
        assert_eq!(paths.len(), 2);
    }

    #[test]
    fn mixed_formats_error() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "a.fastq");
        touch(d.path(), "b.bam");
        let err = classify(d.path(), None).unwrap_err().to_string();
        assert!(err.contains("mixes"));
    }

    #[test]
    fn errors_when_output_is_a_read_file_inside_the_dir() {
        // `-o` pointing at a read file inside the input dir (a real input, or a
        // stale prior output — indistinguishable) must hard-error, not silently
        // exclude+overwrite it.
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "a.fastq");
        touch(d.path(), "b.fastq");
        let a = d.path().join("a.fastq");
        let err = classify(d.path(), Some(a.as_path()))
            .unwrap_err()
            .to_string();
        assert!(err.contains("refusing to overwrite"), "got: {err}");

        // A non-read output (or one outside the dir) does not trip it.
        let (fam, paths) = classify(d.path(), Some(d.path().join("notes.txt").as_path())).unwrap();
        assert_eq!(fam, Family::Fastq);
        assert_eq!(paths.len(), 2);
    }

    #[test]
    fn output_collision_errors_before_the_mixed_check() {
        // A BAM folder with a stale `.fastq` output would otherwise look "mixed";
        // the output-collision error must fire first (and name the overwrite).
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "a.bam");
        touch(d.path(), "merged.fastq");
        let out = d.path().join("merged.fastq");
        let err = classify(d.path(), Some(out.as_path()))
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("refusing to overwrite"),
            "should not be the 'mixes' error: {err}"
        );
    }

    #[test]
    fn empty_folder_error() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "notes.txt"); // no read files
        let err = classify(d.path(), None).unwrap_err().to_string();
        assert!(err.contains("no FASTQ or BAM"));
    }
}
