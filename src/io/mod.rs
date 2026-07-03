pub mod bam;
pub mod dir;
pub mod fastq;

use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fastq,
    FastqGz,
    Bam,
}

/// Extension-based detection. Recognises `.fastq`/`.fq`, the `.gz` variants, and `.bam`.
pub fn from_extension(path: &Path) -> Option<Format> {
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    if name.ends_with(".fastq.gz") || name.ends_with(".fq.gz") {
        Some(Format::FastqGz)
    } else if name.ends_with(".fastq") || name.ends_with(".fq") {
        Some(Format::Fastq)
    } else if name.ends_with(".bam") {
        Some(Format::Bam)
    } else {
        None
    }
}

/// Detect input format from the path extension; if unknown or reading stdin,
/// fall back to sniffing the first bytes.
pub fn detect_input(path: Option<&Path>, sniff: &[u8]) -> anyhow::Result<Format> {
    if let Some(f) = path.and_then(from_extension) {
        return Ok(f);
    }
    if sniff.starts_with(&[0x1f, 0x8b]) {
        Ok(Format::FastqGz)
    } else if sniff.starts_with(b"BAM\x01") {
        Ok(Format::Bam)
    } else if sniff.first() == Some(&b'@') {
        Ok(Format::Fastq)
    } else {
        anyhow::bail!("cannot determine input format; pass --in-format")
    }
}

/// Output format from the path extension, else mirror the input format — with
/// one exception: never auto-compress. A `.gz` (`FastqGz`) input with no output
/// extension defaults to plain `Fastq`, so gzip output only ever happens when the
/// caller explicitly asks (`-o *.gz` / `--out-format fastq-gz`).
pub fn resolve_output(path: Option<&Path>, input: Format) -> Format {
    // An explicit output extension always wins.
    if let Some(f) = path.and_then(from_extension) {
        return f;
    }
    // Otherwise mirror the input EXCEPT never auto-compress: a .gz input
    // defaults to PLAIN fastq output (gz only when explicitly requested via
    // `-o *.gz` above or `--out-format fastq-gz`, which is handled upstream).
    match input {
        Format::FastqGz => Format::Fastq,
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn extensions() {
        assert_eq!(from_extension(Path::new("x.fastq")), Some(Format::Fastq));
        assert_eq!(from_extension(Path::new("x.fq")), Some(Format::Fastq));
        assert_eq!(from_extension(Path::new("x.fastq.gz")), Some(Format::FastqGz));
        assert_eq!(from_extension(Path::new("x.fq.gz")), Some(Format::FastqGz));
        assert_eq!(from_extension(Path::new("x.bam")), Some(Format::Bam));
        assert_eq!(from_extension(Path::new("x.txt")), None);
    }

    #[test]
    fn stdin_sniff_falls_back_to_magic() {
        // no path -> sniff. gzip magic 1f 8b -> FastqGz; '@' -> Fastq; BAM magic -> Bam.
        assert_eq!(detect_input(None, &[0x1f, 0x8b, 0x08]).unwrap(), Format::FastqGz);
        assert_eq!(detect_input(None, b"@read").unwrap(), Format::Fastq);
        assert_eq!(detect_input(None, b"BAM\x01").unwrap(), Format::Bam);
    }

    #[test]
    fn output_mirrors_input_when_no_path() {
        assert_eq!(resolve_output(None, Format::Bam), Format::Bam);
        assert_eq!(resolve_output(None, Format::Fastq), Format::Fastq);
        assert_eq!(resolve_output(Some(Path::new("o.bam")), Format::Fastq), Format::Bam);
    }

    #[test]
    fn output_never_auto_compresses_gz_input() {
        // A .gz input with no output path/format defaults to PLAIN fastq —
        // auto-compressing on stdout would be silent and surprising.
        assert_eq!(resolve_output(None, Format::FastqGz), Format::Fastq);
        // gz output is still available when explicitly requested via a `.gz`
        // output path extension.
        assert_eq!(
            resolve_output(Some(Path::new("o.fastq.gz")), Format::Fastq),
            Format::FastqGz
        );
    }
}
