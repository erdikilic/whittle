//! The startup banner and the `Command:` line.
//!
//! Every line whittle prints before processing starts is built here, as a pure
//! `String` from resolved config, so each one is unit-testable without a run.
//! `obs` owns the end-of-run counterparts; the two bookend on the same
//! `output_desc` text.

use crate::config::{AdapterInfer, AdapterInferAction, AdapterInferPolicy};
use crate::{config, filter, io, qual, trim};

/// The startup banner's operation line (LINE mode's item 3 / BAR mode's own
/// one-liner build on the same wording): `Trimming FASTQ` when input and output
/// share a `Format::family` (including a `FASTQ` -> `FASTQ.gz` run, which is a
/// compression change, not a format conversion), else `Converting {in_label} to
/// {out_label}` (e.g. `Converting BAM to FASTQ`) for a genuine cross-family
/// conversion.
pub(crate) fn operation_line(in_fmt: io::Format, out_fmt: io::Format) -> String {
    if in_fmt.family() == out_fmt.family() {
        format!("Trimming {}", in_fmt.family())
    } else {
        format!("Converting {} to {}", in_fmt.label(), out_fmt.label())
    }
}

/// The startup banner's `Output: ...` line: `Output: <stdout>` when writing to
/// stdout (no compression detail), else `Output: {path}`, with
/// `(gzip|bgzf level {level}, {encode_workers} workers)` appended for compressed
/// output formats (gzip for `FASTQ.gz`, bgzf for BAM; plain FASTQ gets no suffix).
pub(crate) fn output_banner_line(
    output: Option<&std::path::Path>,
    out_fmt: io::Format,
    level: u8,
    encode_workers: usize,
) -> String {
    let Some(path) = output else {
        return "Output: <stdout>".to_string();
    };
    let mut line = format!("Output: {}", path.display());
    match out_fmt {
        io::Format::Bam => {
            line.push_str(&format!(" (bgzf level {level}, {encode_workers} workers)"));
        },
        io::Format::FastqGz => {
            line.push_str(&format!(" (gzip level {level}, {encode_workers} workers)"));
        },
        io::Format::FastqBgzf => {
            line.push_str(&format!(" (bgzf level {level}, {encode_workers} workers)"));
        },
        io::Format::Fastq => {},
    }
    line
}

/// The startup banner's `Threads: ...` line: the resolved worker count, then the
/// per-stage split with `ThreadBudget`'s decode/render/encode renamed to the
/// user-facing read/trim/write, as `Threads: 8 (read 1, trim 4, write 3)`.
///
/// Deliberately not `b.total()`, which floors each stage at >= 1 and so can
/// exceed `-t`, reading as a confusing second, larger total. For the same reason
/// `threads <= 1` prints `Threads: 1 (sequential)` instead of a three-thread
/// split for a run that is single-threaded.
pub(crate) fn threads_banner_line(threads: usize, b: config::ThreadBudget) -> String {
    if threads <= 1 {
        return "Threads: 1 (sequential)".to_string();
    }
    format!(
        "Threads: {threads} (read {}, trim {}, write {})",
        b.decode, b.render, b.encode
    )
}

/// Lowercase label for a `QualMode`, used only in the startup banner's Filters
/// line (`{qual_mode} quality >=...`).
pub(crate) fn qual_mode_label(mode: qual::QualMode) -> &'static str {
    match mode {
        qual::QualMode::Mean => "mean",
        qual::QualMode::Arithmetic => "arithmetic",
        qual::QualMode::Median => "median",
    }
}

/// The startup banner's `Filters: ...; trim: ...` line. Pure, so it is
/// unit-testable directly. Shows only active clauses, so a fresh-defaults run
/// reads `Filters: none; trim: none` rather than spelling out every no-op
/// threshold. A bound appears only when it differs from its default: length from
/// `min_length > 1` / `max_length != usize::MAX`, quality from `min_qual > 0.0` /
/// `max_qual < 1000.0`, GC from either bound being set.
pub(crate) fn filters_and_trim_line(
    filter: &filter::FilterConfig,
    trim: &trim::TrimPlan,
) -> String {
    let mut filters = Vec::new();

    let length_active = filter.min_length > 1 || filter.max_length != usize::MAX;
    if length_active {
        let mut length = String::new();
        if filter.min_length > 1 {
            length.push_str(&format!(">={}", filter.min_length));
        }
        if filter.max_length != usize::MAX {
            if !length.is_empty() {
                length.push(' ');
            }
            length.push_str(&format!("<={}", filter.max_length));
        }
        filters.push(format!("length {length}"));
    }

    let qual_active = filter.min_qual > 0.0 || filter.max_qual < 1000.0;
    if qual_active {
        let mut quality = format!("{} quality", qual_mode_label(filter.qual_mode));
        if filter.min_qual > 0.0 {
            quality.push_str(&format!(" >={}", filter.min_qual));
        }
        if filter.max_qual < 1000.0 {
            quality.push_str(&format!(" <={}", filter.max_qual));
        }
        filters.push(quality);
    }

    if filter.min_gc.is_some() || filter.max_gc.is_some() {
        filters.push(format!(
            "GC {}-{}",
            filter.min_gc.unwrap_or(0.0),
            filter.max_gc.unwrap_or(1.0)
        ));
    }

    let filters_str = if filters.is_empty() {
        "none".to_string()
    } else {
        filters.join("; ")
    };

    let mut trim_parts = Vec::new();
    if trim.head > 0 || trim.tail > 0 {
        trim_parts.push(format!("head {}, tail {}", trim.head, trim.tail));
    }
    if let Some(op) = &trim.quality {
        trim_parts.push(match op {
            trim::QualityOp::TrimQual(q) => format!("trim quality <{q}"),
            trim::QualityOp::BestSegment(q) => format!("best segment >={q}"),
            trim::QualityOp::Split { cutoff, .. } => format!("split quality <{cutoff}"),
        });
    }
    let trim_str = if trim_parts.is_empty() {
        "none".to_string()
    } else {
        trim_parts.join(", ")
    };

    format!("Filters: {filters_str}; trim: {trim_str}")
}

/// The startup banner's `Adapters: ...` line: adapter count, `trim + split` vs
/// `ends-only`, error rate, end-zone size, and whether presence detection will
/// sample. `None` when adapter trimming is off, so the caller skips the line.
///
/// Under inference the count is forced to `0` rather than read off `a.adapters`:
/// in report mode that field may hold the user's FASTA purely as naming refs for
/// `infer::discover`, never as a trimming set.
pub(crate) fn adapter_banner_line(
    adapters: Option<&crate::adapter::AdapterConfig>,
    adapter_sample: usize,
    adapter_infer: AdapterInfer,
) -> Option<String> {
    let a = adapters?;
    let mode = if a.split { "trim + split" } else { "ends-only" };
    let sample = if adapter_sample > 0 {
        format!("sample {adapter_sample}")
    } else {
        "sample off".to_string()
    };
    let infer_suffix = match adapter_infer {
        AdapterInfer::Off => "",
        AdapterInfer::Enabled {
            action: AdapterInferAction::Trim,
            policy: AdapterInferPolicy::Conservative,
        } => " \u{b7} infer trim \u{b7} conservative",
        AdapterInfer::Enabled {
            action: AdapterInferAction::Trim,
            policy: AdapterInferPolicy::Aggressive,
        } => " \u{b7} infer trim \u{b7} aggressive",
        AdapterInfer::Enabled {
            action: AdapterInferAction::Report,
            policy: AdapterInferPolicy::Conservative,
        } => " \u{b7} infer report \u{b7} conservative",
        AdapterInfer::Enabled {
            action: AdapterInferAction::Report,
            policy: AdapterInferPolicy::Aggressive,
        } => " \u{b7} infer report \u{b7} aggressive",
    };
    let n_adapters = if adapter_infer == AdapterInfer::Off {
        a.adapters.len()
    } else {
        0
    };
    Some(format!(
        "Adapters: {} sequences · {mode} · error {:.2} · end-zone {} bp · {sample}{infer_suffix}",
        n_adapters, a.error_rate, a.end_size
    ))
}

/// Shell-quote a single argument the way Python's `shlex.quote` does: bare when
/// non-empty and every character is in the POSIX-shell-safe set
/// (`[A-Za-z0-9_@%+=:,./-]`); otherwise wrapped in single quotes, with any
/// embedded single quote escaped as `'\''` (close the quote, an escaped literal
/// quote, reopen the quote). An empty argument is never safe bare (it would
/// vanish when re-run), so it renders as `''`.
pub(crate) fn shell_quote(arg: &str) -> String {
    let is_safe = |c: char| c.is_ascii_alphanumeric() || "_@%+=:,./-".contains(c);
    if !arg.is_empty() && arg.chars().all(is_safe) {
        arg.to_string()
    } else {
        format!("'{}'", arg.replace('\'', r"'\''"))
    }
}

/// The real process argv, space-joined and shell-quoted so it can be copied back
/// out and re-run. Takes `OsStr`-like items (`args_os()`, not `args()`, which
/// panics on non-UTF-8 argv) and converts lossily here, at the one seam that must
/// never panic on a malformed argv. Generic over the iterator so it is
/// unit-testable without the real argv.
///
/// The value carries no label: the banner prefixes its own `Command: `, while
/// `--summary-json` stores the bare line under its `command` key.
pub fn command_line<I, S>(args: I) -> String
where
    I: IntoIterator<Item = S>,
    S: AsRef<std::ffi::OsStr>,
{
    args.into_iter()
        .map(|a| shell_quote(&a.as_ref().to_string_lossy()))
        .collect::<Vec<_>>()
        .join(" ")
}

/// The output path (or `<stdout>`) shown in both the startup banner's `Output:`
/// line and the end-of-run `Completed`/closer line; the two bookend on the same
/// text so a reader can match them up at a glance.
pub(crate) fn output_desc(output: Option<&std::path::Path>) -> String {
    match output {
        Some(p) => p.display().to_string(),
        None => "<stdout>".to_string(),
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    fn base_filter() -> filter::FilterConfig {
        filter::FilterConfig {
            min_length: 1,
            max_length: usize::MAX,
            min_qual: 0.0,
            max_qual: 1000.0,
            min_gc: None,
            max_gc: None,
            qual_mode: qual::QualMode::Mean,
        }
    }

    fn base_trim() -> trim::TrimPlan {
        trim::TrimPlan {
            head: 0,
            tail: 0,
            quality: None,
        }
    }
    use crate::adapter::{Adapter, AdapterConfig, End};

    #[test]
    fn operation_line_collapses_matching_families() {
        assert_eq!(
            operation_line(io::Format::Fastq, io::Format::Fastq),
            "Trimming FASTQ"
        );
        // FASTQ -> FASTQ.gz shares the FASTQ family (a compression change, not
        // a format conversion), so it collapses too rather than reading as an
        // "X to X" conversion.
        assert_eq!(
            operation_line(io::Format::Fastq, io::Format::FastqGz),
            "Trimming FASTQ"
        );
    }

    #[test]
    fn operation_line_converting_wording_for_cross_family() {
        assert_eq!(
            operation_line(io::Format::Bam, io::Format::Fastq),
            "Converting BAM to FASTQ"
        );
    }

    #[test]
    fn output_banner_line_plain_fastq_has_no_suffix() {
        let p = std::path::Path::new("/tmp/out.fastq");
        assert_eq!(
            output_banner_line(Some(p), io::Format::Fastq, 6, 3),
            "Output: /tmp/out.fastq"
        );
    }

    #[test]
    fn output_banner_line_appends_compression_detail() {
        let p = std::path::Path::new("/tmp/out.fastq.gz");
        assert_eq!(
            output_banner_line(Some(p), io::Format::FastqGz, 6, 4),
            "Output: /tmp/out.fastq.gz (gzip level 6, 4 workers)"
        );
        let p = std::path::Path::new("/tmp/out.bam");
        assert_eq!(
            output_banner_line(Some(p), io::Format::Bam, 3, 5),
            "Output: /tmp/out.bam (bgzf level 3, 5 workers)"
        );
    }

    #[test]
    fn output_banner_line_stdout_has_no_compression_detail() {
        // Even for a format that would otherwise show a compression suffix.
        assert_eq!(
            output_banner_line(None, io::Format::Bam, 6, 3),
            "Output: <stdout>"
        );
    }

    #[test]
    fn threads_banner_line_shows_requested_threads_not_the_stage_sum() {
        let b = config::thread_budget(8, true, false, config::EncodeKind::Bgzf);
        assert_eq!(
            threads_banner_line(8, b),
            format!(
                "Threads: 8 (read {}, trim {}, write {})",
                b.decode, b.render, b.encode
            )
        );
        // Concrete figure too, so a change in `thread_budget`'s split is noticed here.
        assert_eq!(
            threads_banner_line(8, b),
            "Threads: 8 (read 1, trim 4, write 3)"
        );
    }

    #[test]
    fn threads_banner_line_header_is_requested_even_when_stage_sum_differs() {
        // The banner reports the requested limit, not the sum of stage fields.
        let b = config::thread_budget(8, true, false, config::EncodeKind::None);
        assert_eq!(b.total(), 9);
        assert_eq!(
            threads_banner_line(8, b),
            "Threads: 8 (read 1, trim 7, write 1)"
        );
    }

    #[test]
    fn threads_banner_line_sequential_for_one_or_fewer() {
        // `-t 1` (or `-t 0`, which `resolve_threads` floors to 1): the
        // read/trim/write split would otherwise show e.g. "(read 1, trim 1,
        // write 1)" for what is actually a single-threaded run; collapse it
        // to a plain "sequential" label instead.
        let b = config::thread_budget(1, true, false, config::EncodeKind::Bgzf);
        assert_eq!(threads_banner_line(1, b), "Threads: 1 (sequential)");
    }

    #[test]
    fn filters_and_trim_line_defaults() {
        // All-default filter/trim: no active clause, so it reads "none" rather
        // than spelling out no-op thresholds like "mean quality >=0".
        assert_eq!(
            filters_and_trim_line(&base_filter(), &base_trim()),
            "Filters: none; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_min_length_active() {
        let mut f = base_filter();
        f.min_length = 500;
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: length >=500; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_max_length_active() {
        let mut f = base_filter();
        f.max_length = 10_000;
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: length <=10000; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_min_qual_active() {
        let mut f = base_filter();
        f.min_qual = 10.0;
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: mean quality >=10; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_max_qual_active() {
        let mut f = base_filter();
        f.max_qual = 40.0;
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: mean quality <=40; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_gc_active() {
        let mut f = base_filter();
        f.min_gc = Some(0.3);
        assert_eq!(
            filters_and_trim_line(&f, &base_trim()),
            "Filters: GC 0.3-1; trim: none"
        );
    }

    #[test]
    fn filters_and_trim_line_only_trim_active() {
        let mut t = base_trim();
        t.head = 5;
        assert_eq!(
            filters_and_trim_line(&base_filter(), &t),
            "Filters: none; trim: head 5, tail 0"
        );
    }

    #[test]
    fn filters_and_trim_line_quality_ops() {
        let f = base_filter();
        let mut t = base_trim();

        t.quality = Some(trim::QualityOp::BestSegment(20));
        assert!(filters_and_trim_line(&f, &t).ends_with("trim: best segment >=20"));

        t.quality = Some(trim::QualityOp::Split {
            cutoff: 15,
            window: 50,
        });
        assert!(filters_and_trim_line(&f, &t).ends_with("trim: split quality <15"));

        // head/tail-only (no quality op): no trailing quality-op clause.
        t.quality = None;
        t.head = 3;
        t.tail = 0;
        assert!(filters_and_trim_line(&f, &t).ends_with("trim: head 3, tail 0"));
    }

    #[test]
    fn filters_and_trim_line_all_bounds_set() {
        let mut f = base_filter();
        f.min_length = 200;
        f.max_length = 10_000;
        f.min_qual = 8.0;
        f.max_qual = 30.0;
        f.min_gc = Some(0.4);
        f.max_gc = Some(0.6);
        f.qual_mode = qual::QualMode::Median;

        let mut t = base_trim();
        t.head = 10;
        t.tail = 5;
        t.quality = Some(trim::QualityOp::TrimQual(12));

        assert_eq!(
            filters_and_trim_line(&f, &t),
            "Filters: length >=200 <=10000; median quality >=8 <=30; GC 0.4-0.6; \
             trim: head 10, tail 5, trim quality <12"
        );
    }

    #[test]
    fn adapter_banner_line_none_when_off_and_describes_when_on() {
        assert!(adapter_banner_line(None, 10000, AdapterInfer::Off).is_none());
        let cfg = AdapterConfig {
            adapters: vec![Adapter {
                name: "a".into(),
                seq: b"ACGTACGTACGT".to_vec(),
                end: End::Both,
            }],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let line = adapter_banner_line(Some(&cfg), 10000, AdapterInfer::Off).unwrap();
        assert!(line.contains("1 sequences"));
        assert!(line.contains("trim + split"));
        assert!(line.contains("error 0.20"));
        assert!(line.contains("end-zone 150 bp"));
        assert!(line.contains("sample 10000"));
        assert!(!line.contains("infer"), "no infer suffix when off: {line}");

        let off_line = adapter_banner_line(Some(&cfg), 0, AdapterInfer::Off).unwrap();
        assert!(off_line.contains("sample off"));
    }

    #[test]
    fn adapter_banner_line_ends_only_when_split_disabled() {
        let cfg = AdapterConfig {
            adapters: vec![Adapter {
                name: "a".into(),
                seq: b"ACGTACGTACGT".to_vec(),
                end: End::Both,
            }],
            error_rate: 0.2,
            end_size: 150,
            split: false,
            candidate_index: std::sync::OnceLock::new(),
        };
        assert!(
            adapter_banner_line(Some(&cfg), 10000, AdapterInfer::Off)
                .unwrap()
                .contains("ends-only")
        );
    }

    #[test]
    fn adapter_banner_line_notes_infer_mode() {
        let cfg = AdapterConfig {
            adapters: vec![],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        };
        let trim_line = adapter_banner_line(
            Some(&cfg),
            40000,
            AdapterInfer::Enabled {
                action: AdapterInferAction::Trim,
                policy: AdapterInferPolicy::Conservative,
            },
        )
        .unwrap();
        assert!(
            trim_line.ends_with("infer trim · conservative"),
            "{trim_line}"
        );

        let report_line = adapter_banner_line(
            Some(&cfg),
            40000,
            AdapterInfer::Enabled {
                action: AdapterInferAction::Report,
                policy: AdapterInferPolicy::Conservative,
            },
        )
        .unwrap();
        assert!(
            report_line.ends_with("infer report · conservative"),
            "{report_line}"
        );
    }

    #[test]
    fn command_line_quotes_only_unsafe_args() {
        assert_eq!(
            command_line(["whittle", "-i", "in.fastq", "-o", "out.fastq"]),
            "whittle -i in.fastq -o out.fastq"
        );
        assert_eq!(
            command_line(["whittle", "-i", "my reads.fastq"]),
            "whittle -i 'my reads.fastq'"
        );
    }

    #[test]
    fn shell_quote_leaves_plain_args_bare() {
        assert_eq!(shell_quote("whittle"), "whittle");
        assert_eq!(shell_quote("-i"), "-i");
        assert_eq!(shell_quote("in.fastq"), "in.fastq");
        assert_eq!(
            shell_quote("path/to/file_1.0.fq.gz"),
            "path/to/file_1.0.fq.gz"
        );
    }

    #[test]
    fn shell_quote_wraps_args_with_spaces() {
        assert_eq!(shell_quote("my reads.fastq"), "'my reads.fastq'");
    }

    #[test]
    fn shell_quote_escapes_embedded_single_quotes() {
        assert_eq!(shell_quote("it's here.fastq"), r"'it'\''s here.fastq'");
    }

    #[test]
    fn shell_quote_wraps_shell_metacharacters() {
        assert_eq!(shell_quote("$HOME"), "'$HOME'");
        assert_eq!(shell_quote("a;b"), "'a;b'");
    }

    #[test]
    fn shell_quote_wraps_empty_string() {
        // Bare would vanish entirely when the line is re-run.
        assert_eq!(shell_quote(""), "''");
    }

    #[test]
    fn output_desc_stdout_vs_path() {
        assert_eq!(output_desc(None), "<stdout>");
        assert_eq!(
            output_desc(Some(std::path::Path::new("/tmp/out.fastq"))),
            "/tmp/out.fastq"
        );
    }
}
