pub mod bam;
pub mod counting;
pub mod dir;
pub mod fastq;

use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fastq,
    FastqGz,
    Bam,
}

impl Format {
    /// Human-facing label used in log/summary output: `FASTQ`, `FASTQ.gz`, `BAM`.
    pub fn label(&self) -> &'static str {
        match self {
            Format::Fastq => "FASTQ",
            Format::FastqGz => "FASTQ.gz",
            Format::Bam => "BAM",
        }
    }
}

/// Extension-based detection. Recognises `.fastq`/`.fq`, the `.gz` variants, and `.bam`.
pub fn from_extension(path: &Path) -> Option<Format> {
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    // Any trailing `.gz` is gzipped FASTQ — the only gz format this tool handles —
    // so a bare `out.gz` counts as fastq-gz, not just `.fastq.gz`/`.fq.gz`.
    if name.ends_with(".gz") {
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
    // A BAM file is BGZF-compressed, which shares gzip's `1f 8b` magic, so BGZF
    // must be recognised BEFORE the plain-gzip branch — otherwise every real BAM
    // read from stdin or an unrecognised extension is misdetected as gzipped
    // FASTQ (the `BAM\x01` magic only appears *after* bgzf decompression).
    if is_bgzf(sniff) {
        Ok(Format::Bam)
    } else if sniff.starts_with(&[0x1f, 0x8b]) {
        Ok(Format::FastqGz)
    } else if sniff.starts_with(b"BAM\x01") {
        // A naked (non-BGZF) BAM stream — unusual, but accept it rather than fail.
        Ok(Format::Bam)
    } else if sniff.first() == Some(&b'@') {
        Ok(Format::Fastq)
    } else {
        anyhow::bail!("cannot determine input format; pass --in-format")
    }
}

/// True if `sniff` begins with a BGZF block header: gzip magic + deflate method +
/// the `FEXTRA` flag carrying the mandatory `BC` subfield. This distinguishes a
/// (BGZF) BAM from a plain-gzip FASTQ, which shares the leading `1f 8b` but sets
/// neither `FEXTRA` nor `BC`. Requires the full 18-byte block header to be present.
fn is_bgzf(sniff: &[u8]) -> bool {
    sniff.len() >= 18
        && sniff[0] == 0x1f
        && sniff[1] == 0x8b
        && sniff[2] == 0x08          // CM = DEFLATE
        && (sniff[3] & 0x04) != 0    // FLG.FEXTRA set
        && sniff[12] == b'B'         // first extra subfield SI1
        && sniff[13] == b'C' //                      SI2 -> "BC" = BGZF
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
    use std::path::Path;

    use super::*;

    #[test]
    fn format_labels() {
        assert_eq!(Format::Fastq.label(), "FASTQ");
        assert_eq!(Format::FastqGz.label(), "FASTQ.gz");
        assert_eq!(Format::Bam.label(), "BAM");
    }

    #[test]
    fn extensions() {
        assert_eq!(from_extension(Path::new("x.fastq")), Some(Format::Fastq));
        assert_eq!(from_extension(Path::new("x.fq")), Some(Format::Fastq));
        assert_eq!(
            from_extension(Path::new("x.fastq.gz")),
            Some(Format::FastqGz)
        );
        assert_eq!(from_extension(Path::new("x.fq.gz")), Some(Format::FastqGz));
        assert_eq!(from_extension(Path::new("x.gz")), Some(Format::FastqGz)); // bare .gz
        assert_eq!(from_extension(Path::new("x.bam")), Some(Format::Bam));
        assert_eq!(from_extension(Path::new("x.txt")), None);
    }

    #[test]
    fn stdin_sniff_falls_back_to_magic() {
        // no path -> sniff. gzip magic 1f 8b -> FastqGz; '@' -> Fastq; BAM magic -> Bam.
        assert_eq!(
            detect_input(None, &[0x1f, 0x8b, 0x08]).unwrap(),
            Format::FastqGz
        );
        assert_eq!(detect_input(None, b"@read").unwrap(), Format::Fastq);
        assert_eq!(detect_input(None, b"BAM\x01").unwrap(), Format::Bam);
    }

    #[test]
    fn bgzf_header_sniffs_as_bam_not_gz() {
        // A real BAM is BGZF, which starts with the gzip magic but sets FLG.FEXTRA
        // and carries a "BC" subfield. It must be detected as Bam, not FastqGz.
        let mut bgzf = vec![
            0x1f, 0x8b, 0x08, 0x04, // magic, CM=deflate, FLG=FEXTRA
            0x00, 0x00, 0x00, 0x00, // MTIME
            0x00, 0xff, // XFL, OS
            0x06, 0x00, // XLEN = 6
            b'B', b'C', 0x02, 0x00, // "BC" subfield, SLEN=2
            0x1b, 0x00, // BSIZE
        ];
        assert_eq!(detect_input(None, &bgzf).unwrap(), Format::Bam);

        // A plain-gzip FASTQ stream (FLG=0, no BC) must still be FastqGz even with
        // a full-length header present.
        bgzf[3] = 0x00; // clear FEXTRA
        assert_eq!(detect_input(None, &bgzf).unwrap(), Format::FastqGz);

        // Too-short gzip-magic buffer can't be BGZF -> defaults to FastqGz.
        assert_eq!(
            detect_input(None, &[0x1f, 0x8b, 0x08, 0x04]).unwrap(),
            Format::FastqGz
        );
    }

    #[test]
    fn output_mirrors_input_when_no_path() {
        assert_eq!(resolve_output(None, Format::Bam), Format::Bam);
        assert_eq!(resolve_output(None, Format::Fastq), Format::Fastq);
        assert_eq!(
            resolve_output(Some(Path::new("o.bam")), Format::Fastq),
            Format::Bam
        );
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
