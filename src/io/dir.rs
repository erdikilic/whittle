use std::path::{Path, PathBuf};

use anyhow::Context;
use noodles_sam::{self as sam, alignment::RecordBuf};

use crate::io::{Format, from_extension};
use crate::record::ReadRecord;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    Fastq,
    Bam,
}

/// Classify a directory's immediate children into a single read-file family and
/// a sorted path list. Non-read files and subdirectories are ignored. Errors if
/// the folder mixes FASTQ and BAM, or contains no read files.
pub fn classify(dir: &Path) -> anyhow::Result<(Family, Vec<PathBuf>)> {
    let mut fastq = Vec::new();
    let mut bam = Vec::new();

    let entries = std::fs::read_dir(dir)
        .with_context(|| format!("reading directory {}", dir.display()))?;
    for entry in entries {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        match from_extension(&path) {
            Some(Format::Fastq | Format::FastqGz) => fastq.push(path),
            Some(Format::Bam) => bam.push(path),
            None => {} // ignore non-read files
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

/// A boxed, owning iterator over decoded `RecordBuf`s (or per-record errors).
/// Named to satisfy `clippy::type_complexity` on the `bam_reader` signature below.
type BamRecordIter = Box<dyn Iterator<Item = anyhow::Result<RecordBuf>>>;

/// The first file's header plus one chained record stream over all BAM files
/// (each file's own header is read and discarded past the first; records
/// stream under the first header — `samtools cat` semantics for homogeneous
/// uBAM).
///
/// Returns an `Err` if `paths` is empty rather than panicking. Each file is
/// opened exactly once: the first file's record iterator (obtained alongside
/// its header) is reused via `chain` rather than reopening that file.
pub fn bam_reader(paths: &[PathBuf]) -> anyhow::Result<(sam::Header, BamRecordIter)> {
    let (first, rest) = paths
        .split_first()
        .ok_or_else(|| anyhow::anyhow!("bam_reader called with no BAM files"))?;
    let (header, first_records) = crate::io::bam::reader(Some(first))?;
    let rest = rest.to_vec();
    let rest_records = rest.into_iter().flat_map(|p| -> BamRecordIter {
        match crate::io::bam::reader(Some(&p)) {
            Ok((_hdr, recs)) => recs,
            Err(e) => Box::new(std::iter::once(Err(e))),
        }
    });
    let records: BamRecordIter = Box::new(first_records.chain(rest_records));
    Ok((header, records))
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
        let (fam, paths) = classify(d.path()).unwrap();
        assert_eq!(fam, Family::Fastq);
        let names: Vec<_> = paths.iter().map(|p| p.file_name().unwrap().to_str().unwrap()).collect();
        assert_eq!(names, vec!["a_0.fastq", "b_1.fastq.gz"]); // sorted, .txt/subdir excluded
    }

    #[test]
    fn bam_only_folder() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "x.bam");
        touch(d.path(), "y.bam");
        let (fam, paths) = classify(d.path()).unwrap();
        assert_eq!(fam, Family::Bam);
        assert_eq!(paths.len(), 2);
    }

    #[test]
    fn mixed_formats_error() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "a.fastq");
        touch(d.path(), "b.bam");
        let err = classify(d.path()).unwrap_err().to_string();
        assert!(err.contains("mixes"));
    }

    #[test]
    fn empty_folder_error() {
        let d = tempfile::tempdir().unwrap();
        touch(d.path(), "notes.txt"); // no read files
        let err = classify(d.path()).unwrap_err().to_string();
        assert!(err.contains("no FASTQ or BAM"));
    }
}
