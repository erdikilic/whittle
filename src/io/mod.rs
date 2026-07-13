pub mod bam;
pub mod counting;
pub mod dir;
pub mod fastq;

use std::io::Read;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fastq,
    FastqGz,
    FastqBgzf,
    Bam,
}

impl Format {
    /// Human-facing label used in log/summary output: `FASTQ`, `FASTQ.gz`, `BAM`.
    pub fn label(&self) -> &'static str {
        match self {
            Format::Fastq => "FASTQ",
            Format::FastqGz => "FASTQ.gz",
            Format::FastqBgzf => "FASTQ.bgz",
            Format::Bam => "BAM",
        }
    }

    /// Coarse format family, used to decide whether an (in, out) pair reads as
    /// a "conversion" in the startup banner's operation line: the two FASTQ
    /// variants collapse to the same `"FASTQ"` family (a `Fastq` ->
    /// `FastqGz` run is a compression change, not a format conversion), while
    /// `Bam` is its own family.
    pub fn family(&self) -> &'static str {
        match self {
            Format::Fastq | Format::FastqGz | Format::FastqBgzf => "FASTQ",
            Format::Bam => "BAM",
        }
    }
}

/// Extension-based detection. Recognises `.fastq`/`.fq`, the `.gz` variants, and `.bam`.
pub fn from_extension(path: &Path) -> Option<Format> {
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    // Any trailing `.gz` is gzipped FASTQ — the only gz format this tool handles —
    // so a bare `out.gz` counts as fastq-gz, not just `.fastq.gz`/`.fq.gz`.
    if name.ends_with(".bgz") || name.ends_with(".bgzf") {
        Some(Format::FastqBgzf)
    } else if name.ends_with(".gz") {
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
    // BAM is BGZF-compressed and shares gzip's magic, so test BGZF first. The
    // `BAM\x01` magic becomes visible only after decompression.
    if is_bgzf(sniff) {
        Ok(Format::Bam)
    } else if sniff.starts_with(&[0x1f, 0x8b]) {
        Ok(Format::FastqGz)
    } else if sniff.starts_with(b"BAM\x01") {
        // A naked (non-BGZF) BAM stream: the reader always wraps input in a
        // bgzf decoder, so this could never actually be read — fail with a
        // precise message instead of surfacing an opaque bgzf framing error.
        anyhow::bail!(
            "input looks like an uncompressed (non-BGZF) BAM stream; a BGZF-compressed BAM \
             is required (re-compress with `samtools view -b`)"
        )
    } else if sniff.first() == Some(&b'@') {
        Ok(Format::Fastq)
    } else {
        anyhow::bail!("cannot determine input format; pass --in-format")
    }
}

/// Advisory text when an explicit `--in-format`/`--out-format` (`forced`)
/// disagrees with what `path`'s extension suggests — e.g. `--out-format fastq`
/// on an `out.fastq.gz` path. `None` when there's no forced format, the path
/// has no recognized extension, it's stdin/stdout (`path` is `None`), or the
/// two agree. `flag` names the CLI flag for the message.
pub fn format_mismatch_warning(
    flag: &str,
    forced: Option<Format>,
    path: Option<&Path>,
) -> Option<String> {
    let forced = forced?;
    let detected = from_extension(path?)?;
    (detected != forced).then(|| {
        format!(
            "{flag} {} but the file extension looks like {}",
            forced.label(),
            detected.label()
        )
    })
}

/// True if `sniff` begins with a BGZF block header: gzip magic + deflate method +
/// the `FEXTRA` flag carrying the mandatory `BC` subfield. This distinguishes a
/// (BGZF) BAM from a plain-gzip FASTQ, which shares the leading `1f 8b` but sets
/// neither `FEXTRA` nor `BC`. Requires the full 18-byte block header to be present.
pub(crate) fn is_bgzf(sniff: &[u8]) -> bool {
    sniff.len() >= 18
        && sniff[0] == 0x1f
        && sniff[1] == 0x8b
        && sniff[2] == 0x08          // CM = DEFLATE
        && (sniff[3] & 0x04) != 0    // FLG.FEXTRA set
        && sniff[12] == b'B'         // first extra subfield SI1
        && sniff[13] == b'C' //                      SI2 -> "BC" = BGZF
}

/// Identify the payload carried by one complete BGZF block. BAM begins with
/// `BAM\x01`; FASTQ begins with `@`. The caller replays the original compressed
/// block into the selected reader after this probe.
pub(crate) fn detect_bgzf_block(block: &[u8]) -> anyhow::Result<Format> {
    let mut reader = noodles_bgzf::io::Reader::new(std::io::Cursor::new(block));
    let mut probe = [0u8; 4];
    reader.read_exact(&mut probe)?;
    if probe.starts_with(b"BAM\x01") {
        Ok(Format::Bam)
    } else if probe.first() == Some(&b'@') {
        Ok(Format::FastqBgzf)
    } else {
        anyhow::bail!("BGZF input is neither BAM nor FASTQ")
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
        Format::FastqGz | Format::FastqBgzf => Format::Fastq,
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
        assert_eq!(Format::FastqBgzf.label(), "FASTQ.bgz");
        assert_eq!(Format::Bam.label(), "BAM");
    }

    #[test]
    fn format_families() {
        assert_eq!(Format::Fastq.family(), "FASTQ");
        assert_eq!(Format::FastqGz.family(), "FASTQ");
        assert_eq!(Format::FastqBgzf.family(), "FASTQ");
        assert_eq!(Format::Bam.family(), "BAM");
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
        assert_eq!(
            from_extension(Path::new("x.fastq.bgz")),
            Some(Format::FastqBgzf)
        );
        assert_eq!(
            from_extension(Path::new("x.fq.bgzf")),
            Some(Format::FastqBgzf)
        );
        assert_eq!(from_extension(Path::new("x.gz")), Some(Format::FastqGz)); // bare .gz
        assert_eq!(from_extension(Path::new("x.bam")), Some(Format::Bam));
        assert_eq!(from_extension(Path::new("x.txt")), None);
    }

    #[test]
    fn stdin_sniff_falls_back_to_magic() {
        // no path -> sniff. gzip magic 1f 8b -> FastqGz; '@' -> Fastq.
        assert_eq!(
            detect_input(None, &[0x1f, 0x8b, 0x08]).unwrap(),
            Format::FastqGz
        );
        assert_eq!(detect_input(None, b"@read").unwrap(), Format::Fastq);
    }

    #[test]
    fn naked_non_bgzf_bam_is_rejected() {
        // A bare `BAM\x01` stream (no BGZF framing) cannot be read by the
        // bgzf-wrapping reader, so detection must fail with a clear message
        // rather than claim Format::Bam and then surface an opaque bgzf error.
        let err = detect_input(None, b"BAM\x01rest").unwrap_err().to_string();
        assert!(
            err.to_ascii_lowercase().contains("bgzf"),
            "message should name BGZF, got: {err}"
        );
    }

    #[test]
    fn format_mismatch_warning_fires_on_disagreement() {
        let w = format_mismatch_warning(
            "--out-format",
            Some(Format::Fastq),
            Some(Path::new("out.fastq.gz")),
        );
        assert_eq!(
            w.as_deref(),
            Some("--out-format FASTQ but the file extension looks like FASTQ.gz")
        );
    }

    #[test]
    fn format_mismatch_warning_silent_when_absent_or_agreeing() {
        // no forced format, matching extension, unknown extension, and stdin all stay silent.
        assert_eq!(
            format_mismatch_warning("--in-format", None, Some(Path::new("x.bam"))),
            None
        );
        assert_eq!(
            format_mismatch_warning("--in-format", Some(Format::Bam), Some(Path::new("x.bam"))),
            None
        );
        assert_eq!(
            format_mismatch_warning("--in-format", Some(Format::Bam), Some(Path::new("x.txt"))),
            None
        );
        assert_eq!(
            format_mismatch_warning("--in-format", Some(Format::Bam), None),
            None
        );
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
    fn complete_bgzf_fastq_block_is_identified() {
        use std::io::Write;

        let mut writer = noodles_bgzf::io::Writer::new(Vec::new());
        writer.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
        let compressed = writer.finish().unwrap();
        let first_block_size =
            usize::from(u16::from_le_bytes([compressed[16], compressed[17]])) + 1;
        assert_eq!(
            detect_bgzf_block(&compressed[..first_block_size]).unwrap(),
            Format::FastqBgzf
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
