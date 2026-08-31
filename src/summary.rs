//! Machine-readable run summary (`--summary-json`).
//!
//! One JSON object per run: the resolved settings (`params`), what was consumed
//! and produced (`reads`, `bases`), and why anything was dropped
//! (`segments_dropped`). Written regardless of `--quiet` or the log level. The
//! field names are interface: `schema_version` is bumped when an existing field
//! changes meaning or disappears, but not when one is added.

use std::path::Path;

use serde::Serialize;

use crate::config::{AdapterInfer, AdapterInferAction, AdapterInferPolicy, Config, FastqTags};
use crate::qual::QualMode;
use crate::trim::QualityOp;
use crate::workflow::Stats;

/// Incremented only when an existing field changes meaning or is removed.
const SCHEMA_VERSION: u32 = 1;

#[derive(Debug, Serialize)]
pub struct Summary {
    schema_version: u32,
    tool: &'static str,
    version: &'static str,
    /// The invocation, shell-quoted so it can be copied back out and re-run.
    command: String,
    /// Input path, or `<stdin>`.
    input: String,
    /// Output path, or `<stdout>`.
    output: String,
    /// Wall-clock processing time. `None` when the caller never started a timer.
    elapsed_seconds: Option<f64>,
    params: Params,
    reads: Reads,
    bases: Bases,
    segments_dropped: SegmentsDropped,
    warnings: Warnings,
}

/// The run's resolved settings, after defaults and clamping. Reading these back
/// rather than re-deriving them from `command` is what makes the file a
/// self-contained provenance record.
#[derive(Debug, Serialize)]
struct Params {
    threads: usize,
    compression_level: u8,
    min_length: usize,
    /// `None` when `--max-length` was not given.
    max_length: Option<usize>,
    min_qual: f64,
    max_qual: f64,
    min_gc: Option<f64>,
    max_gc: Option<f64>,
    qual_mode: &'static str,
    head_crop: usize,
    tail_crop: usize,
    /// `None` when no quality-trimming strategy was selected.
    quality_op: Option<QualityOpParams>,
    update_moves: bool,
    /// `all`, `none`, or the comma-joined tag list.
    fastq_tags: String,
    /// `None` when adapter trimming is off.
    adapters: Option<AdapterParams>,
}

#[derive(Debug, Serialize)]
struct QualityOpParams {
    /// `trim`, `best_segment`, or `split`.
    mode: &'static str,
    threshold: u8,
    /// Only meaningful for `split`; `None` otherwise.
    window: Option<usize>,
}

#[derive(Debug, Serialize)]
struct AdapterParams {
    /// Adapters configured before resolution: the preset and/or FASTA asked for.
    /// This is the figure the startup banner prints, since it is all that is
    /// known before reads have been sampled. `0` under inference, where the set
    /// is discovered rather than configured.
    configured: usize,
    /// Adapters actually trimmed against, after presence detection narrowed the
    /// configured set or inference replaced it. Equal to `configured` when
    /// neither ran.
    count: usize,
    error_rate: f64,
    end_size: usize,
    /// False under `--adapter-ends-only`, or under conservative inference.
    split: bool,
    /// Reads sampled for presence detection or inference; `0` disables detection.
    sample: usize,
    /// `off`, `trim`, or `report`.
    infer: &'static str,
    /// `None` when inference is off.
    infer_policy: Option<&'static str>,
}

#[derive(Debug, Serialize)]
struct Reads {
    input: u64,
    /// Output segments written. A `--qual-split` read can contribute several,
    /// so this can legitimately exceed `input`.
    output: u64,
    /// Input reads that produced at least one surviving segment.
    with_output: u64,
    /// Input reads that produced no segments at all (empty, fully consumed by
    /// adapter trimming, or over-cropped).
    trimmed_to_nothing: u64,
    /// Input reads that produced segments, all of which the filters rejected.
    all_filtered: u64,
}

#[derive(Debug, Serialize)]
struct Bases {
    input: u64,
    output: u64,
}

/// Segment-level filter rejections, by reason. One input read can contribute to
/// more than one of these when it splits, so these do not sum to a read count.
#[derive(Debug, Serialize)]
struct SegmentsDropped {
    too_short: u64,
    too_long: u64,
    low_quality: u64,
    high_quality: u64,
    gc_out_of_range: u64,
}

#[derive(Debug, Serialize)]
struct Warnings {
    /// Reads whose per-base kinetics tag length disagreed with the sequence and
    /// was left untouched.
    malformed_tag_reads: u64,
}

impl Summary {
    /// Build the summary from a finished run. `command` is the shell-quoted
    /// invocation (see `crate::command_line`), `output` the output path or
    /// `<stdout>` label already computed for the banner.
    pub fn new(
        cfg: &Config,
        stats: &Stats,
        command: String,
        output: String,
        elapsed: Option<std::time::Duration>,
    ) -> Self {
        Summary {
            schema_version: SCHEMA_VERSION,
            tool: "whittle",
            version: env!("CARGO_PKG_VERSION"),
            command,
            input: cfg
                .io
                .input
                .as_deref()
                .map_or_else(|| "<stdin>".to_string(), |p| p.display().to_string()),
            output,
            elapsed_seconds: elapsed.map(|d| d.as_secs_f64()),
            params: Params::from_config(cfg),
            reads: Reads {
                input: stats.input_reads,
                output: stats.output_reads,
                with_output: stats
                    .input_reads
                    .saturating_sub(stats.reads_trimmed_to_nothing)
                    .saturating_sub(stats.reads_all_filtered),
                trimmed_to_nothing: stats.reads_trimmed_to_nothing,
                all_filtered: stats.reads_all_filtered,
            },
            bases: Bases {
                input: stats.input_bases,
                output: stats.output_bases,
            },
            segments_dropped: SegmentsDropped {
                too_short: stats.segments_dropped_short,
                too_long: stats.segments_dropped_long,
                low_quality: stats.segments_dropped_low_qual,
                high_quality: stats.segments_dropped_high_qual,
                gc_out_of_range: stats.segments_dropped_gc,
            },
            warnings: Warnings {
                malformed_tag_reads: stats.malformed_tag_reads,
            },
        }
    }

    /// Serialize to pretty-printed JSON with a trailing newline, so the file is
    /// both diff-friendly and safe to `cat` in a shell.
    pub fn to_json(&self) -> anyhow::Result<String> {
        let mut s = serde_json::to_string_pretty(self)?;
        s.push('\n');
        Ok(s)
    }

    /// Write the summary to `path`. A failure here is a hard error: a caller
    /// that passed `--summary-json` is a pipeline that needs the file.
    pub fn write(&self, path: &Path) -> anyhow::Result<()> {
        let json = self.to_json()?;
        std::fs::write(path, json)
            .map_err(|e| anyhow::anyhow!("writing summary JSON to {}: {e}", path.display()))
    }
}

impl Params {
    fn from_config(cfg: &Config) -> Self {
        Params {
            threads: cfg.threads,
            compression_level: cfg.compression_level,
            min_length: cfg.filter.min_length,
            max_length: (cfg.filter.max_length != usize::MAX).then_some(cfg.filter.max_length),
            min_qual: cfg.filter.min_qual,
            max_qual: cfg.filter.max_qual,
            min_gc: cfg.filter.min_gc,
            max_gc: cfg.filter.max_gc,
            qual_mode: match cfg.filter.qual_mode {
                QualMode::Mean => "mean",
                QualMode::Arithmetic => "arithmetic",
                QualMode::Median => "median",
            },
            head_crop: cfg.trim.head,
            tail_crop: cfg.trim.tail,
            quality_op: cfg.trim.quality.as_ref().map(|op| match op {
                QualityOp::TrimQual(q) => QualityOpParams {
                    mode: "trim",
                    threshold: *q,
                    window: None,
                },
                QualityOp::BestSegment(q) => QualityOpParams {
                    mode: "best_segment",
                    threshold: *q,
                    window: None,
                },
                QualityOp::Split { cutoff, window } => QualityOpParams {
                    mode: "split",
                    threshold: *cutoff,
                    window: Some(*window),
                },
            }),
            update_moves: cfg.update_moves,
            fastq_tags: match &cfg.fastq_tags {
                FastqTags::All => "all".to_string(),
                FastqTags::None => "none".to_string(),
                FastqTags::Only(tags) => tags
                    .iter()
                    .map(|t| String::from_utf8_lossy(t).into_owned())
                    .collect::<Vec<_>>()
                    .join(","),
            },
            adapters: cfg.adapters.as_ref().map(|ac| AdapterParams {
                // Falls back to the resolved count for a caller that drove the
                // workflow directly and never went through `settle`.
                configured: cfg.adapters_configured.unwrap_or(ac.adapters.len()),
                count: ac.adapters.len(),
                error_rate: ac.error_rate,
                end_size: ac.end_size,
                split: ac.split,
                sample: cfg.adapter_sample,
                infer: match cfg.adapter_infer {
                    AdapterInfer::Off => "off",
                    AdapterInfer::Enabled {
                        action: AdapterInferAction::Trim,
                        ..
                    } => "trim",
                    AdapterInfer::Enabled {
                        action: AdapterInferAction::Report,
                        ..
                    } => "report",
                },
                infer_policy: match cfg.adapter_infer {
                    AdapterInfer::Off => None,
                    AdapterInfer::Enabled {
                        policy: AdapterInferPolicy::Conservative,
                        ..
                    } => Some("conservative"),
                    AdapterInfer::Enabled {
                        policy: AdapterInferPolicy::Aggressive,
                        ..
                    } => Some("aggressive"),
                },
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;
    use crate::config::IoConfig;
    use crate::filter::FilterConfig;
    use crate::trim::TrimPlan;

    fn cfg() -> Config {
        Config {
            io: IoConfig {
                input: Some("reads.bam".into()),
                output: Some("out.fastq".into()),
                in_format: None,
                out_format: None,
            },
            filter: FilterConfig {
                min_length: 500,
                max_length: usize::MAX,
                min_qual: 10.0,
                max_qual: 1000.0,
                min_gc: None,
                max_gc: None,
                qual_mode: QualMode::Mean,
            },
            trim: TrimPlan {
                head: 20,
                tail: 20,
                quality: Some(QualityOp::Split {
                    cutoff: 9,
                    window: 50,
                }),
            },
            adapters: None,
            adapter_infer: AdapterInfer::Off,
            threads: 8,
            fastq_tags: FastqTags::All,
            render_workers: 0,
            adapter_sample: 0,
            compression_level: 6,
            update_moves: false,
            verbosity: 0,
            quiet: false,
            threads_clamped: None,
            summary_json: None,
            advisories: Vec::new(),
            adapter_fasta: None,
            adapters_configured: None,
        }
    }

    fn stats() -> Stats {
        Stats {
            input_reads: 100,
            output_reads: 110,
            input_bases: 10_000,
            output_bases: 9_000,
            malformed_tag_reads: 2,
            reads_trimmed_to_nothing: 5,
            reads_all_filtered: 3,
            segments_dropped_short: 7,
            segments_dropped_long: 0,
            segments_dropped_low_qual: 1,
            segments_dropped_high_qual: 0,
            segments_dropped_gc: 0,
        }
    }

    fn value(summary: &Summary) -> serde_json::Value {
        serde_json::from_str(&summary.to_json().unwrap()).unwrap()
    }

    /// The three read-level buckets partition `input_reads`, so `with_output` is
    /// derived rather than counted, and must not go negative on a partial run.
    #[test]
    fn read_buckets_partition_the_input() {
        let s = Summary::new(
            &cfg(),
            &stats(),
            "whittle -i reads.bam".into(),
            "out.fastq".into(),
            None,
        );
        let v = value(&s);
        assert_eq!(v["reads"]["input"], 100);
        assert_eq!(v["reads"]["with_output"], 92);
        assert_eq!(v["reads"]["trimmed_to_nothing"], 5);
        assert_eq!(v["reads"]["all_filtered"], 3);
        // Output segments can exceed input reads under --qual-split.
        assert_eq!(v["reads"]["output"], 110);
    }

    #[test]
    fn derived_read_bucket_saturates_instead_of_underflowing() {
        let mut st = stats();
        st.input_reads = 1;
        st.reads_trimmed_to_nothing = 5;
        st.reads_all_filtered = 3;
        let s = Summary::new(&cfg(), &st, String::new(), String::new(), None);
        assert_eq!(value(&s)["reads"]["with_output"], 0);
    }

    #[test]
    fn params_record_the_resolved_run() {
        let s = Summary::new(&cfg(), &stats(), String::new(), String::new(), None);
        let v = value(&s);
        assert_eq!(v["schema_version"], 1);
        assert_eq!(v["tool"], "whittle");
        assert_eq!(v["params"]["threads"], 8);
        assert_eq!(v["params"]["min_length"], 500);
        assert_eq!(v["params"]["qual_mode"], "mean");
        assert_eq!(v["params"]["head_crop"], 20);
        assert_eq!(v["params"]["quality_op"]["mode"], "split");
        assert_eq!(v["params"]["quality_op"]["threshold"], 9);
        assert_eq!(v["params"]["quality_op"]["window"], 50);
        assert_eq!(v["params"]["fastq_tags"], "all");
        // An unset --max-length is null, not usize::MAX leaking into the file.
        assert!(v["params"]["max_length"].is_null());
        // Adapter trimming off is null, not a zeroed block.
        assert!(v["params"]["adapters"].is_null());
    }

    #[test]
    fn adapter_block_present_when_trimming_is_on() {
        let mut c = cfg();
        c.adapters = Some(crate::adapter::AdapterConfig {
            adapters: Vec::new(),
            error_rate: 0.2,
            end_size: 150,
            split: true,
            candidate_index: std::sync::OnceLock::new(),
        });
        c.adapter_sample = 5_000;
        c.adapters_configured = Some(124);
        let s = Summary::new(&c, &stats(), String::new(), String::new(), None);
        let v = value(&s);
        // Both figures are reported: detection narrowed 124 down to 0 here, and a
        // reader can see that rather than guessing which number they are looking at.
        assert_eq!(v["params"]["adapters"]["configured"], 124);
        assert_eq!(v["params"]["adapters"]["count"], 0);
        assert_eq!(v["params"]["adapters"]["error_rate"], 0.2);
        assert_eq!(v["params"]["adapters"]["end_size"], 150);
        assert_eq!(v["params"]["adapters"]["split"], true);
        assert_eq!(v["params"]["adapters"]["sample"], 5_000);
        assert_eq!(v["params"]["adapters"]["infer"], "off");
        assert!(v["params"]["adapters"]["infer_policy"].is_null());
    }

    #[test]
    fn fastq_tag_list_renders_as_a_comma_list() {
        let mut c = cfg();
        c.fastq_tags = FastqTags::Only(BTreeSet::from([*b"MM", *b"ML", *b"RG"]));
        let s = Summary::new(&c, &stats(), String::new(), String::new(), None);
        assert_eq!(value(&s)["params"]["fastq_tags"], "ML,MM,RG");
    }

    #[test]
    fn elapsed_is_seconds_and_omitted_when_unknown() {
        let s = Summary::new(
            &cfg(),
            &stats(),
            String::new(),
            String::new(),
            Some(std::time::Duration::from_millis(2500)),
        );
        assert_eq!(value(&s)["elapsed_seconds"], 2.5);

        let s = Summary::new(&cfg(), &stats(), String::new(), String::new(), None);
        assert!(value(&s)["elapsed_seconds"].is_null());
    }

    /// stdin/stdout runs still name their endpoints, so a summary file is never
    /// ambiguous about where the reads came from.
    #[test]
    fn stdin_and_stdout_get_explicit_labels() {
        let mut c = cfg();
        c.io.input = None;
        let s = Summary::new(&c, &stats(), String::new(), "<stdout>".into(), None);
        let v = value(&s);
        assert_eq!(v["input"], "<stdin>");
        assert_eq!(v["output"], "<stdout>");
    }

    #[test]
    fn write_creates_a_parseable_file_ending_in_a_newline() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("summary.json");
        Summary::new(&cfg(), &stats(), String::new(), String::new(), None)
            .write(&path)
            .unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.ends_with("}\n"));
        let v: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(v["bases"]["input"], 10_000);
        assert_eq!(v["segments_dropped"]["too_short"], 7);
        assert_eq!(v["warnings"]["malformed_tag_reads"], 2);
    }
}
