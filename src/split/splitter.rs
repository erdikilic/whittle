//! Splitter: builds the keys and rescoring state of a sheet, then calls each
//! trimmed piece by rescoring the sheet's primers at its located loci and
//! names the output key it is written to.

use crate::adapter::{AdapterConfig, Locus};
use crate::trim::Piece;

use super::classify::{Bounds, Call, KeyLevel, Keys, Rules, check_length, classify_at};
use super::route::{KeyTable, Owner, Template, UNCLASSIFIED, barcode_name};
use super::score::{End, Scorer};
use super::sheet::Sheet;
use crate::workflow::KeyId;

/// What `--split-by` does with a classified piece.
#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum Action {
    /// Keep the piece as trimmed.
    Trim,
    /// Keep the piece widened over its located primer loci.
    Retain,
}

impl Action {
    /// The CLI spelling of this action, as printed in the summary.
    pub fn label(self) -> &'static str {
        match self {
            Action::Trim => "trim",
            Action::Retain => "retain",
        }
    }
}

/// Resolved `--split-by` configuration for a run.
#[derive(Debug, Clone)]
pub struct SplitOptions {
    /// The parsed primer split sheet.
    pub sheet: Sheet,
    /// Key granularity: one key per target or one key per group.
    pub level: KeyLevel,
    /// The end rule and lead margin the classifier applies.
    pub rules: Rules,
    /// Whether an assigned piece is trimmed or retained.
    pub action: Action,
    /// Drops reads classified `Unassigned` instead of writing them to the
    /// `unassigned` bin.
    pub discard_unassigned: bool,
    /// Drops reads classified `Ambiguous` instead of writing them to the
    /// `ambiguous` bin.
    pub discard_ambiguous: bool,
    /// The `--split-by` values as given, in order, carried through to the
    /// report.
    pub spec: Vec<String>,
    /// The `-o` output template, or `None` in tag-only mode, where every
    /// record goes to the one output.
    pub template: Option<Template>,
}

/// Builds a sheet's keys and rescoring state once and classifies every piece
/// against them.
#[derive(Debug)]
pub struct Splitter {
    /// The resolved `--split-by` configuration.
    pub opts: SplitOptions,
    /// The sheet's keys at `opts.level`.
    pub keys: Keys,
    /// Rescores the sheet's primers at a piece's located loci.
    pub scorer: Scorer,
    /// The per-primer budgets and close primers the classifier weighs.
    pub bounds: Bounds,
    /// The output paths `opts.template` has expanded to, one per output key.
    pub table: KeyTable,
    /// How a call's output key is found.
    routes: Routes,
}

/// How `Splitter::output_key` finds a call's output key.
#[derive(Debug)]
enum Routes {
    /// Tag-only mode: every call goes to key 0.
    Single,
    /// A template without `{barcode}`: the key of each bin (each key, then
    /// `unassigned` and `ambiguous`), interned once at construction.
    Fixed(Vec<KeyId>),
    /// A template with `{barcode}`: each call's path is expanded and interned.
    Barcoded,
}

impl Splitter {
    /// Builds the keys and scorer of `opts.sheet`, scoring with `adapters`'
    /// error rate and end size, and records the sequences of `adapters`,
    /// which locus sites index (`Scorer::with_sites`). For a template
    /// without `{barcode}`, the path of every bin is interned here; two bins
    /// expanding to one path are an error.
    pub fn new(opts: SplitOptions, adapters: &AdapterConfig) -> anyhow::Result<Splitter> {
        let keys = Keys::new(&opts.sheet, opts.level);
        let scorer = Scorer::new(&opts.sheet, adapters.error_rate, adapters.end_size)
            .with_sites(&adapters.adapters);
        let bounds = Bounds::new(
            &opts.sheet,
            &keys,
            opts.rules.lead,
            &scorer.budgets(),
            &scorer.anchored_budgets(),
        );
        let table = KeyTable::default();
        let routes = match &opts.template {
            None => Routes::Single,
            Some(template) if template.barcode() => Routes::Barcoded,
            Some(template) => Routes::Fixed(template.intern_bins(&opts.sheet, &keys, &table)?),
        };
        Ok(Splitter {
            opts,
            keys,
            scorer,
            bounds,
            table,
            routes,
        })
    }

    /// Whether output keys depend on the record's `BC:Z` barcode call: the
    /// template holds `{barcode}`.
    pub fn reads_barcode(&self) -> bool {
        matches!(self.routes, Routes::Barcoded)
    }

    /// The index of `call`'s bin, in the order `Template::intern_bins` uses.
    fn bin(&self, call: &Call) -> usize {
        let keys = self.keys.names.len();
        match call {
            Call::Assigned { key, .. } => *key,
            Call::Unassigned(_) => keys,
            Call::Ambiguous => keys + 1,
        }
    }

    /// The output key of a written `call`: 0 in tag-only mode; the bin's
    /// precomputed key for a template without `{barcode}`; otherwise the
    /// interned path the template expands to for its bin, its group and the
    /// record's `BC:Z` barcode call `barcode` (`unclassified` when `None`).
    /// A barcode call that is not one path component, and a path already
    /// interned for a different bin or barcode, are errors.
    pub fn output_key(&self, call: &Call, barcode: Option<&[u8]>) -> anyhow::Result<KeyId> {
        let template = match &self.routes {
            Routes::Single => return Ok(0),
            Routes::Fixed(ids) => return Ok(ids[self.bin(call)]),
            Routes::Barcoded => self
                .opts
                .template
                .as_ref()
                .expect("A barcoded route has a template"),
        };
        let barcode = barcode
            .map(barcode_name)
            .transpose()?
            .unwrap_or(UNCLASSIFIED);
        let label = self.label(call);
        let group = match call {
            Call::Assigned { target, .. } => self.opts.sheet.targets[*target].group.as_str(),
            Call::Unassigned(_) | Call::Ambiguous => label,
        };
        let owner = Owner {
            key: label,
            group,
            barcode: Some(barcode),
        };
        self.table.intern(
            &template.expand_in_group(label, group, Some(barcode)),
            owner,
        )
    }

    /// Classifies `piece`: rescores the sheet's primers at each located
    /// locus and its site (`Scorer::score_at_site`), with its outer edge as
    /// a read end where the locus is `outer_open` (an end without a locus
    /// scores an empty slice) and at the anchored budget where it is
    /// `boundary`, classifies the result with the penalty of each end taken
    /// from the budgets rescoring used there, and applies the length check
    /// to `piece.end - piece.start`.
    pub fn call(&self, seq: &[u8], piece: &Piece) -> Call {
        let score = |locus: Option<Locus>, end: End| match locus {
            Some(l) => self.scorer.score_at_site(
                seq,
                (l.start, l.end),
                end,
                l.outer_open,
                l.boundary,
                l.site,
            ),
            None => Vec::new(),
        };
        let five = score(piece.five, End::Five);
        let three = score(piece.three, End::Three);
        let boundary = |locus: Option<Locus>| locus.is_some_and(|l| l.boundary);
        let call = classify_at(
            &self.opts.sheet,
            &self.keys,
            &self.bounds,
            self.opts.rules,
            &five,
            &three,
            [boundary(piece.five), boundary(piece.three)],
        );
        check_length(&self.opts.sheet, call, piece.end - piece.start)
    }

    /// The bin name of `call`: the assigned key's name, `"unassigned"`, or
    /// `"ambiguous"`.
    pub fn label(&self, call: &Call) -> &str {
        match call {
            Call::Assigned { key, .. } => self.keys.names[*key].as_str(),
            Call::Unassigned(_) => "unassigned",
            Call::Ambiguous => "ambiguous",
        }
    }

    /// Whether `call` is dropped rather than written, under
    /// `opts.discard_unassigned` and `opts.discard_ambiguous`.
    pub fn discards(&self, call: &Call) -> bool {
        match call {
            Call::Assigned { .. } => false,
            Call::Unassigned(_) => self.opts.discard_unassigned,
            Call::Ambiguous => self.opts.discard_ambiguous,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adapter::{Segment, adapter_segments_annotated, reverse_complement};
    use crate::split::classify::{Ends, Strand, Unassigned};
    use crate::split::sheet::Require;

    /// Generates deterministic SplitMix64 bases.
    fn splitmix_dna(seed: u64, len: usize) -> Vec<u8> {
        let mut state = 0x9E37_79B9_7F4A_7C15u64.wrapping_add(seed);
        (0..len)
            .map(|_| {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                b"ACGT"[((z >> 62) & 0b11) as usize]
            })
            .collect()
    }

    /// Builds `left`, then `mid`, then `right`.
    fn joined(left: &[u8], mid: &[u8], right: &[u8]) -> Vec<u8> {
        let mut out = left.to_vec();
        out.extend_from_slice(mid);
        out.extend_from_slice(right);
        out
    }

    /// A two-target sheet: `A` (fwd `fA`, rev `rA`) and `B` (fwd `fB`, rev
    /// `rB`), naming no group and no length window.
    fn two_target_sheet() -> Sheet {
        let (fa, ra) = (splitmix_dna(801, 20), splitmix_dna(802, 22));
        let (fb, rb) = (splitmix_dna(803, 21), splitmix_dna(804, 20));
        let text = format!(
            "target\tfwd\trev\n\
             A\t{}\t{}\n\
             B\t{}\t{}\n",
            String::from_utf8_lossy(&fa),
            String::from_utf8_lossy(&ra),
            String::from_utf8_lossy(&fb),
            String::from_utf8_lossy(&rb),
        );
        Sheet::parse_tsv(&text).unwrap()
    }

    /// An adapter configuration with `sheet`'s primers attached and no other
    /// adapters, at error rate 0.2 and end zone 150.
    fn adapter_config(sheet: &Sheet) -> AdapterConfig {
        let mut cfg = AdapterConfig {
            adapters: Vec::new(),
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
            split_of: Vec::new(),
            split_opens: Vec::new(),
        };
        cfg.attach_split(&sheet.primers);
        cfg
    }

    /// Resolved `--split-by` options with default action and no discards.
    fn options(sheet: Sheet, require: Require, lead: usize) -> SplitOptions {
        SplitOptions {
            sheet,
            level: KeyLevel::Target,
            rules: Rules { require, lead },
            action: Action::Trim,
            discard_unassigned: false,
            discard_ambiguous: false,
            spec: vec!["sheet.tsv".into()],
            template: None,
        }
    }

    /// Returns `read`'s single adapter segment as a `Piece`, panicking if
    /// `adapter_segments_annotated` finds a different count.
    fn single_piece(read: &[u8], cfg: &AdapterConfig) -> Piece {
        let segments = adapter_segments_annotated(read, cfg);
        let [
            Segment {
                start,
                end,
                five,
                three,
            },
        ] = segments[..]
        else {
            panic!("expected one segment, got {segments:?}");
        };
        Piece {
            start,
            end,
            five,
            three,
        }
    }

    #[test]
    fn call_assigns_plus_and_minus_reads() {
        let sheet = two_target_sheet();
        let fa = sheet.primers[sheet.targets[0].fwd[0]].seq.clone();
        let ra = sheet.primers[sheet.targets[0].rev[0]].seq.clone();
        let cfg = adapter_config(&sheet);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();

        let plus = joined(&fa, &splitmix_dna(811, 400), &reverse_complement(&ra));
        let piece = single_piece(&plus, &cfg);
        assert_eq!(
            splitter.call(&plus, &piece),
            Call::Assigned {
                key: 0,
                target: 0,
                strand: Strand::Plus,
                ends: Ends::Both,
            }
        );

        let minus = joined(&ra, &splitmix_dna(812, 400), &reverse_complement(&fa));
        let piece = single_piece(&minus, &cfg);
        assert_eq!(
            splitter.call(&minus, &piece),
            Call::Assigned {
                key: 0,
                target: 0,
                strand: Strand::Minus,
                ends: Ends::Both,
            }
        );
    }

    /// The pieces `trim::apply_pieces` returns carry loci in read
    /// coordinates, so each is called on the full read even after a fixed
    /// crop and a quality split move its bounds.
    #[test]
    fn apply_pieces_output_calls_on_the_full_read() {
        use crate::trim::{QualityOp, TrimPlan, apply_pieces};

        let sheet = two_target_sheet();
        let fb = sheet.primers[sheet.targets[1].fwd[0]].seq.clone();
        let rb = sheet.primers[sheet.targets[1].rev[0]].seq.clone();
        let cfg = adapter_config(&sheet);
        let splitter = Splitter::new(options(sheet, Require::Both, 2), &cfg).unwrap();

        let insert = splitmix_dna(814, 600);
        let read = joined(&fb, &insert, &reverse_complement(&rb));
        let mut phred = vec![40u8; read.len()];
        let low = fb.len() + 300;
        phred[low..low + 5].fill(2);
        let plan = TrimPlan {
            head: 7,
            tail: 0,
            quality: Some(QualityOp::runs(10, 5)),
        };
        let pieces = apply_pieces(&read, &phred, &plan, Some(&cfg), None, false);
        assert_eq!(pieces.len(), 2, "{pieces:?}");
        assert_eq!(pieces[0].start, fb.len() + 7);
        let expected = Call::Assigned {
            key: 1,
            target: 1,
            strand: Strand::Plus,
            ends: Ends::Both,
        };
        for piece in &pieces {
            assert_eq!(splitter.call(&read, piece), expected, "{piece:?}");
        }
    }

    #[test]
    fn call_unassigned_without_primer() {
        let sheet = two_target_sheet();
        let cfg = adapter_config(&sheet);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();

        let read = splitmix_dna(813, 400);
        let piece = single_piece(&read, &cfg);
        assert_eq!(piece.five, None);
        assert_eq!(piece.three, None);
        assert_eq!(
            splitter.call(&read, &piece),
            Call::Unassigned(Unassigned::NoPrimer)
        );
    }

    #[test]
    fn label_names() {
        let sheet = Sheet::parse_tsv(&format!(
            "target\tfwd\trev\n16S\t{}\t{}\n",
            String::from_utf8_lossy(&splitmix_dna(821, 20)),
            String::from_utf8_lossy(&splitmix_dna(822, 20)),
        ))
        .unwrap();
        let cfg = adapter_config(&sheet);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();

        let assigned = Call::Assigned {
            key: 0,
            target: 0,
            strand: Strand::Plus,
            ends: Ends::Both,
        };
        assert_eq!(splitter.label(&assigned), "16S");
        assert_eq!(
            splitter.label(&Call::Unassigned(Unassigned::NoPrimer)),
            "unassigned"
        );
        assert_eq!(splitter.label(&Call::Ambiguous), "ambiguous");
    }

    #[test]
    fn output_key_interns_the_expanded_path() {
        let sheet = two_target_sheet();
        let cfg = adapter_config(&sheet);
        let tag_only = Splitter::new(options(sheet.clone(), Require::Either, 2), &cfg).unwrap();
        let a = Call::Assigned {
            key: 0,
            target: 0,
            strand: Strand::Plus,
            ends: Ends::Both,
        };
        assert_eq!(tag_only.output_key(&a, Some(b"bc1")).unwrap(), 0);
        assert!(!tag_only.reads_barcode());

        // Without `{barcode}`, every bin's key is interned up front.
        let mut opts = options(sheet.clone(), Require::Either, 2);
        opts.template = Template::parse(std::path::Path::new("out/{target}.fq")).unwrap();
        let fixed = Splitter::new(opts, &cfg).unwrap();
        assert!(!fixed.reads_barcode());
        let paths: Vec<_> = (0..4).map(|id| fixed.table.path(id)).collect();
        assert_eq!(
            paths,
            [
                "out/A.fq",
                "out/B.fq",
                "out/unassigned.fq",
                "out/ambiguous.fq"
            ]
            .map(std::path::PathBuf::from)
        );
        assert_eq!(fixed.output_key(&a, Some(b"../x")).unwrap(), 0);
        assert_eq!(fixed.output_key(&Call::Ambiguous, None).unwrap(), 3);

        let mut opts = options(sheet, Require::Either, 2);
        opts.template = Template::parse(std::path::Path::new("{barcode}/{target}.fq")).unwrap();
        let splitter = Splitter::new(opts, &cfg).unwrap();
        assert!(splitter.reads_barcode());
        let none = Call::Unassigned(Unassigned::NoPrimer);
        let first = splitter.output_key(&none, Some(b"bc1")).unwrap();
        let second = splitter.output_key(&a, None).unwrap();
        assert_ne!(first, second);
        assert_eq!(splitter.output_key(&a, None).unwrap(), second);
        assert_eq!(
            splitter.table.path(first),
            std::path::PathBuf::from("bc1/unassigned.fq")
        );
        assert_eq!(
            splitter.table.path(second),
            std::path::PathBuf::from("unclassified/A.fq")
        );
        assert!(splitter.output_key(&a, Some(b"../x")).is_err());
    }

    /// Returns `read`'s adapter segments as pieces.
    fn pieces(read: &[u8], cfg: &AdapterConfig) -> Vec<Piece> {
        adapter_segments_annotated(read, cfg)
            .into_iter()
            .map(|s| Piece {
                start: s.start,
                end: s.end,
                five: s.five,
                three: s.three,
            })
            .collect()
    }

    /// The sequences of the two-target sheet: `fA`, `rA`, `fB`, `rB`.
    fn primer_seqs(sheet: &Sheet) -> [Vec<u8>; 4] {
        let seq = |list: &[usize]| sheet.primers[list[0]].seq.clone();
        [
            seq(&sheet.targets[0].fwd),
            seq(&sheet.targets[0].rev),
            seq(&sheet.targets[1].fwd),
            seq(&sheet.targets[1].rev),
        ]
    }

    /// The call of a plus-strand read of target `target` with primers at
    /// both ends.
    fn plus_both(target: usize) -> Call {
        Call::Assigned {
            key: target,
            target,
            strand: Strand::Plus,
            ends: Ends::Both,
        }
    }

    /// An amplicon of target `A` whose insert holds the whole amplicon of
    /// target `B`, both primers intact, is one amplicon: a sheet primer inside
    /// a read splits it only at a junction pair.
    #[test]
    fn nested_site_does_not_split() {
        let sheet = two_target_sheet();
        let [fa, ra, fb, rb] = primer_seqs(&sheet);
        let cfg = adapter_config(&sheet);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();

        for inner in [100, 400] {
            let read = joined(
                &fa,
                &[
                    splitmix_dna(831, 400),
                    fb.clone(),
                    splitmix_dna(832, inner),
                    reverse_complement(&rb),
                    splitmix_dna(833, 400),
                ]
                .concat(),
                &reverse_complement(&ra),
            );
            let pieces = pieces(&read, &cfg);
            assert_eq!(pieces.len(), 1, "inner {inner}: {pieces:?}");
            assert_eq!(
                (pieces[0].start, pieces[0].end),
                (fa.len(), read.len() - ra.len())
            );
            assert_eq!(splitter.call(&read, &pieces[0]), plus_both(0));
        }
    }

    /// Two amplicons joined directly hold a junction pair, the closing primer
    /// of one beside the opening primer of the next, and split there into
    /// one piece per target.
    #[test]
    fn junction_pair_still_splits() {
        let sheet = two_target_sheet();
        let [fa, ra, fb, rb] = primer_seqs(&sheet);
        let cfg = adapter_config(&sheet);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();

        let read = joined(
            &fa,
            &[
                splitmix_dna(834, 500),
                reverse_complement(&ra),
                fb.clone(),
                splitmix_dna(835, 500),
            ]
            .concat(),
            &reverse_complement(&rb),
        );
        let pieces = pieces(&read, &cfg);
        assert_eq!(pieces.len(), 2, "{pieces:?}");
        assert_eq!(splitter.call(&read, &pieces[0]), plus_both(0));
        assert_eq!(splitter.call(&read, &pieces[1]), plus_both(1));
    }

    /// A sheet primer in the closing orientation inside the 5' end zone, just
    /// behind the opening primer, neither trims nor becomes the locus.
    #[test]
    fn wrong_orientation_inner_hit_ignored() {
        let sheet = two_target_sheet();
        let [fa, ra, _, rb] = primer_seqs(&sheet);
        let cfg = adapter_config(&sheet);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();

        for gap in 0..=5 {
            let read = joined(
                &fa,
                &[
                    splitmix_dna(836, gap),
                    reverse_complement(&rb),
                    splitmix_dna(837, 400),
                ]
                .concat(),
                &reverse_complement(&ra),
            );
            let pieces = pieces(&read, &cfg);
            assert_eq!(pieces.len(), 1, "gap {gap}: {pieces:?}");
            assert_eq!(pieces[0].start, fa.len(), "gap {gap}: {pieces:?}");
            let five = pieces[0].five.unwrap();
            assert_eq!((five.start, five.end), (0, fa.len()), "gap {gap}");
            assert_eq!(splitter.call(&read, &pieces[0]), plus_both(0), "gap {gap}");
        }
    }

    /// Of two valid 3' primers, the outermost is the locus and sets the trim:
    /// the primer 20 bases further inside and the bases between stay.
    #[test]
    fn outermost_primer_is_the_locus() {
        let sheet = two_target_sheet();
        let [fa, ra, _, rb] = primer_seqs(&sheet);
        let cfg = adapter_config(&sheet);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();

        let read = joined(
            &fa,
            &[
                splitmix_dna(838, 400),
                reverse_complement(&rb),
                splitmix_dna(839, 20),
            ]
            .concat(),
            &reverse_complement(&ra),
        );
        let n = read.len();
        let pieces = pieces(&read, &cfg);
        assert_eq!(pieces.len(), 1, "{pieces:?}");
        assert_eq!(pieces[0].end, n - ra.len(), "{pieces:?}");
        let three = pieces[0].three.unwrap();
        assert_eq!((three.start, three.end), (n - ra.len(), n));
        assert_eq!(splitter.call(&read, &pieces[0]), plus_both(0));
    }

    /// A primer that lost its first bases at the adapter junction is located
    /// at the adapter's trim boundary and scored there, as a primer cut short
    /// by the read start is. The adapter ends in bases that the lost primer
    /// bases do not match, so the whole primer does not align within its
    /// budget there.
    #[test]
    fn clipped_primer_behind_adapter_scores() {
        let sheet = two_target_sheet();
        let [fa, ra, _, _] = primer_seqs(&sheet);
        let adapter = [splitmix_dna(842, 23), b"AAAAA".to_vec()].concat();
        let mut cfg = AdapterConfig {
            adapters: vec![crate::adapter::Adapter {
                name: "adapter".into(),
                seq: adapter.clone(),
                role: crate::adapter::Role::Adapter,
            }],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
            split_of: Vec::new(),
            split_opens: Vec::new(),
        };
        cfg.attach_split(&sheet.primers);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();

        let head = [splitmix_dna(840, 30), adapter].concat();
        let read = joined(
            &head,
            &[fa[5..].to_vec(), splitmix_dna(841, 400)].concat(),
            &reverse_complement(&ra),
        );
        let budget = crate::adapter::terminal_budget(&fa, 0.2, 150);
        let mut searcher = crate::adapter::search::new_ambiguous_searcher();
        let whole = crate::adapter::search::hits(
            &mut searcher,
            &fa,
            &read[..head.len() + fa.len()],
            budget,
        );
        assert!(whole.is_empty(), "{whole:?}");

        let pieces = pieces(&read, &cfg);
        assert_eq!(pieces.len(), 1, "{pieces:?}");
        let five = pieces[0].five.expect("the clipped primer is located");
        assert_eq!(
            (five.start, five.end),
            (head.len(), head.len() + fa.len() - 5)
        );
        assert!(five.outer_open, "{five:?}");
        let fa_idx = splitter.opts.sheet.targets[0].fwd[0];
        let scores = splitter
            .scorer
            .score(&read, (five.start, five.end), End::Five, true, false);
        assert!(scores.iter().any(|s| s.primer == fa_idx), "{scores:?}");
        assert_eq!(splitter.call(&read, &pieces[0]), plus_both(0));
    }

    /// Two forward primers on one site, the neighbour one base longer on its
    /// outer side, one edit apart: the true primer, whole behind an adapter
    /// with four edits in its first bases, is called for its own target. It
    /// aligns whole, so its outer edge stays closed: neither primer may hang
    /// its first bases off the adapter boundary for free, and the neighbour
    /// pays for the adapter base it lacks, beyond its budget. The two primers
    /// are closer than the default lead, so the one-edit lead decides at
    /// lead 1.
    #[test]
    fn whole_primer_behind_adapter_keeps_its_lead_over_a_shifted_neighbour() {
        let site = splitmix_dna(850, 21);
        let (neighbour, fwd) = (site.clone(), site[1..].to_vec());
        let text = format!(
            "target\tfwd\trev\n\
             A\t{}\t{}\n\
             B\t{}\t{}\n",
            String::from_utf8_lossy(&fwd),
            String::from_utf8_lossy(&splitmix_dna(851, 20)),
            String::from_utf8_lossy(&neighbour),
            String::from_utf8_lossy(&splitmix_dna(852, 20)),
        );
        let sheet = Sheet::parse_tsv(&text).unwrap();
        let mut adapter = splitmix_dna(853, 28);
        adapter[27] = if site[0] == b'A' { b'C' } else { b'A' };
        let mut cfg = AdapterConfig {
            adapters: vec![crate::adapter::Adapter {
                name: "adapter".into(),
                seq: adapter.clone(),
                role: crate::adapter::Role::Adapter,
            }],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
            split_of: Vec::new(),
            split_opens: Vec::new(),
        };
        cfg.attach_split(&sheet.primers);
        let splitter = Splitter::new(options(sheet, Require::Either, 1), &cfg).unwrap();

        let mut edited = fwd.clone();
        for i in [0, 2, 4, 6] {
            edited[i] = if edited[i] == b'A' { b'C' } else { b'A' };
        }
        let head = [splitmix_dna(854, 30), adapter].concat();
        let read = joined(&head, &edited, &splitmix_dna(855, 400));
        let pieces = pieces(&read, &cfg);
        assert_eq!(pieces.len(), 1, "{pieces:?}");
        let five = pieces[0].five.expect("the primer is located");
        assert_eq!(
            splitter.call(&read, &pieces[0]),
            Call::Assigned {
                key: 0,
                target: 0,
                strand: Strand::Plus,
                ends: Ends::Five,
            },
            "{five:?}"
        );
        assert!(!five.outer_open, "{five:?}");
    }

    /// The MAB114 adapter preset with the `mab114` sheet attached, at error
    /// rate 0.2 and end zone 150, and a splitter over that sheet.
    fn mab114_config() -> (AdapterConfig, Splitter) {
        let sheet = Sheet::preset("mab114").unwrap();
        let mut cfg = AdapterConfig {
            adapters: crate::adapter::preset::preset(&[crate::adapter::preset::Kit::Mab114]),
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
            split_of: Vec::new(),
            split_opens: Vec::new(),
        };
        cfg.attach_split(&sheet.primers);
        let splitter = Splitter::new(options(sheet, Require::Either, 2), &cfg).unwrap();
        (cfg, splitter)
    }

    /// Returns the catalog sequence of `name`.
    fn kit(name: &str) -> Vec<u8> {
        crate::adapter::catalog::sequence_of(name).unwrap().to_vec()
    }

    /// Returns the anchored budget of the `mab114` sheet primer `name`, over
    /// the sheet's primers together (`anchored_budgets`).
    fn mab114_anchored(name: &str) -> usize {
        let sheet = Sheet::preset("mab114").unwrap();
        let seqs: Vec<&[u8]> = sheet.primers.iter().map(|p| p.seq.as_slice()).collect();
        let at = sheet.primers.iter().position(|p| p.name == name).unwrap();
        crate::adapter::anchored_budgets(&seqs, 0.2)[at]
    }

    /// Returns `seq` with a different base at each of `positions`.
    fn substituted(seq: &[u8], positions: &[usize]) -> Vec<u8> {
        let mut out = seq.to_vec();
        for &i in positions {
            out[i] = if out[i] == b'A' { b'C' } else { b'A' };
        }
        out
    }

    /// The MAB114 barcode layers that read into the insert: the rapid
    /// adapter behind `lead` bases, the front flank, barcode 5 and the rear
    /// flank.
    fn mab114_head(lead: usize) -> Vec<u8> {
        [
            splitmix_dna(870, lead),
            kit("RAD"),
            kit("RBK4_front"),
            kit("TP05"),
            kit("MAB_rear"),
        ]
        .concat()
    }

    /// The call of a read of the 16S target located at `ends` on `strand`.
    fn call_16s(splitter: &Splitter, strand: Strand, ends: Ends) -> Call {
        let target = splitter
            .opts
            .sheet
            .targets
            .iter()
            .position(|t| t.name == "16S")
            .unwrap();
        Call::Assigned {
            key: target,
            target,
            strand,
            ends,
        }
    }

    /// A kit primer with four edits, within its own terminal budget but
    /// above the set-wide budget of its site, flush behind the rear flank of
    /// the barcode layers is located at that trim boundary and assigned, at
    /// the 5' end and, reverse complemented, at the 3' end.
    #[test]
    fn kit_primer_at_its_own_budget_behind_the_flank_is_assigned() {
        let (cfg, splitter) = mab114_config();
        let fwd = substituted(&kit("16S_Bor_F"), &[3, 8, 12, 17]);
        let rev = substituted(&kit("16S_Bor_R"), &[2, 7, 11, 16]);
        assert_eq!(crate::adapter::terminal_budget(&fwd, 0.2, 150), 4);

        let head = mab114_head(20);
        let read = [head.clone(), fwd.clone(), splitmix_dna(871, 1400)].concat();
        let found = pieces(&read, &cfg);
        assert_eq!(found.len(), 1, "{found:?}");
        let five = found[0].five.expect("the primer is located");
        assert_eq!((five.start, five.end), (head.len(), head.len() + fwd.len()));
        assert_eq!(found[0].start, head.len() + fwd.len());
        assert_eq!(
            splitter.call(&read, &found[0]),
            call_16s(&splitter, Strand::Plus, Ends::Five)
        );

        let tail = reverse_complement(&mab114_head(20));
        let read = [splitmix_dna(872, 1400), reverse_complement(&rev), tail].concat();
        let found = pieces(&read, &cfg);
        assert_eq!(found.len(), 1, "{found:?}");
        let three = found[0].three.expect("the primer is located");
        assert_eq!((three.start, three.end), (1400, 1400 + rev.len()));
        assert_eq!(
            splitter.call(&read, &found[0]),
            call_16s(&splitter, Strand::Plus, Ends::Three)
        );
    }

    /// The same four-edit primer inside the 5' end zone is not located where
    /// no trim boundary sits beside it: at a read without adapters, and
    /// thirty bases behind the rear flank.
    #[test]
    fn kit_primer_at_its_own_budget_away_from_a_boundary_is_not_located() {
        let (cfg, splitter) = mab114_config();
        let fwd = substituted(&kit("16S_Bor_F"), &[3, 8, 12, 17]);
        let heads = [
            splitmix_dna(873, 60),
            [mab114_head(20), splitmix_dna(874, 30)].concat(),
        ];
        for head in heads {
            let read = [head, fwd.clone(), splitmix_dna(875, 1400)].concat();
            let pieces = pieces(&read, &cfg);
            assert_eq!(pieces.len(), 1, "{pieces:?}");
            assert_eq!(pieces[0].five, None, "{pieces:?}");
            assert_eq!(
                splitter.call(&read, &pieces[0]),
                Call::Unassigned(Unassigned::NoPrimer)
            );
        }
    }

    /// A kit primer that the barcode layers push past the 5' end zone, and
    /// one on a read short enough that both end zones cover it, is the
    /// terminal primer of the trim boundary it sits on and is assigned.
    #[test]
    fn kit_primer_pushed_past_the_end_zone_is_assigned() {
        let (cfg, splitter) = mab114_config();
        let fwd = kit("16S_Ent_F");
        let deep = mab114_head(70);
        assert!(deep.len() > 150);
        let short = [
            splitmix_dna(876, 32),
            kit("RBK4_front"),
            kit("TP05"),
            kit("MAB_rear"),
        ]
        .concat();
        for (head, insert) in [(deep, 1400), (short, 143)] {
            let read = [head.clone(), fwd.clone(), splitmix_dna(877, insert)].concat();
            let pieces = pieces(&read, &cfg);
            assert_eq!(pieces.len(), 1, "{pieces:?}");
            let five = pieces[0].five.expect("the primer is located");
            assert_eq!((five.start, five.end), (head.len(), head.len() + fwd.len()));
            assert_eq!(
                splitter.call(&read, &pieces[0]),
                call_16s(&splitter, Strand::Plus, Ends::Five)
            );
        }
    }

    /// A lone kit primer inside the insert of a MAB114 read, away from every
    /// trim boundary, neither splits the read nor becomes a locus.
    #[test]
    fn lone_interior_kit_primer_does_not_split() {
        let (cfg, splitter) = mab114_config();
        let head = mab114_head(20);
        let fwd = kit("16S_Ent_F");
        for offset in [150, 600] {
            let read = [
                head.clone(),
                fwd.clone(),
                splitmix_dna(878, offset),
                kit("16S_Bor_F"),
                splitmix_dna(879, 1200 - offset),
            ]
            .concat();
            let pieces = pieces(&read, &cfg);
            assert_eq!(pieces.len(), 1, "offset {offset}: {pieces:?}");
            assert_eq!(pieces[0].three, None, "offset {offset}: {pieces:?}");
            assert_eq!(
                splitter.call(&read, &pieces[0]),
                call_16s(&splitter, Strand::Plus, Ends::Five)
            );
        }
    }

    /// The degenerate 16S_mix_F primer, resolved to plain bases, with four
    /// substitutions: above its own terminal budget, within its anchored
    /// budget.
    fn mix_primer_four_edits() -> Vec<u8> {
        let primer = kit("16S_mix_F");
        assert_eq!(crate::adapter::terminal_budget(&primer, 0.2, 150), 3);
        assert_eq!(mab114_anchored("16S_mix_F"), 4);
        substituted(b"AGAGTTTGATCATGGCTCAG", &[4, 8, 13, 17])
    }

    /// A primer above its own terminal budget and within its anchored budget,
    /// flush behind the rear flank of the barcode layers, is located at that
    /// trim boundary as a boundary locus and assigned, at the 5' end and,
    /// reverse complemented, at the 3' end.
    #[test]
    fn primer_at_its_anchored_budget_behind_the_flank_is_assigned() {
        let (cfg, splitter) = mab114_config();
        let fwd = mix_primer_four_edits();

        let head = mab114_head(20);
        let read = [head.clone(), fwd.clone(), splitmix_dna(880, 1400)].concat();
        let found = pieces(&read, &cfg);
        assert_eq!(found.len(), 1, "{found:?}");
        let five = found[0].five.expect("the primer is located");
        assert_eq!((five.start, five.end), (head.len(), head.len() + fwd.len()));
        assert!(five.boundary, "{five:?}");
        assert_eq!(
            splitter.call(&read, &found[0]),
            call_16s(&splitter, Strand::Plus, Ends::Five)
        );

        let tail = reverse_complement(&mab114_head(20));
        let read = [splitmix_dna(881, 1400), reverse_complement(&fwd), tail].concat();
        let found = pieces(&read, &cfg);
        assert_eq!(found.len(), 1, "{found:?}");
        let three = found[0].three.expect("the primer is located");
        assert_eq!((three.start, three.end), (1400, 1400 + fwd.len()));
        assert!(three.boundary, "{three:?}");
        assert_eq!(
            splitter.call(&read, &found[0]),
            call_16s(&splitter, Strand::Minus, Ends::Three)
        );
    }

    /// The same primer forty bases behind the rear flank, away from every
    /// trim boundary, is not located, and the read is unassigned.
    #[test]
    fn primer_at_its_anchored_budget_away_from_a_boundary_is_not_located() {
        let (cfg, splitter) = mab114_config();
        let head = [mab114_head(20), splitmix_dna(882, 40)].concat();
        let read = [head, mix_primer_four_edits(), splitmix_dna(883, 1400)].concat();
        let pieces = pieces(&read, &cfg);
        assert_eq!(pieces.len(), 1, "{pieces:?}");
        assert_eq!(pieces[0].five, None, "{pieces:?}");
        assert_eq!(
            splitter.call(&read, &pieces[0]),
            Call::Unassigned(Unassigned::NoPrimer)
        );
    }

    /// Reads of random bases between `head` and `tail`, drawn per sheet.
    const CHANCE_TRIALS: u64 = 1_000;

    /// Asserts that the anchored budgets of `sheet` admit at most the
    /// per-read terminal chance target of whole chance hits at the trim
    /// boundaries (`anchored_chance_per_read`), and that the boundary loci
    /// (`Locus::boundary`) `cfg` locates at both ends of `CHANCE_TRIALS`
    /// reads of random bases between `head` and `tail` stay within three
    /// times that expectation.
    fn assert_boundary_loci_within_chance(
        sheet: &Sheet,
        cfg: &AdapterConfig,
        head: &[u8],
        tail: &[u8],
    ) {
        let seqs: Vec<&[u8]> = sheet.primers.iter().map(|p| p.seq.as_slice()).collect();
        let expected = crate::adapter::anchored_chance_per_read(&seqs, cfg.error_rate);
        assert!(
            expected <= crate::adapter::TERMINAL_CHANCE_HITS_PER_READ,
            "{expected}"
        );
        let found: usize = (0..CHANCE_TRIALS)
            .map(|seed| {
                let read = [head, &splitmix_dna(10_000 + seed, 220), tail].concat();
                pieces(&read, cfg)
                    .iter()
                    .flat_map(|p| [p.five, p.three])
                    .flatten()
                    .filter(|l| l.boundary)
                    .count()
            })
            .sum();
        let bound = 3.0 * expected * CHANCE_TRIALS as f64;
        assert!(
            found as f64 <= bound,
            "{found} boundary loci in {CHANCE_TRIALS} reads, bound {bound:.2}"
        );
    }

    /// Random bases between the MAB114 barcode layers at both ends hold
    /// boundary loci within the chance the anchored budgets admit.
    #[test]
    fn random_bases_between_mab114_flanks_stay_within_the_chance_target() {
        let (cfg, splitter) = mab114_config();
        let head = mab114_head(20);
        let tail = reverse_complement(&head);
        assert_boundary_loci_within_chance(&splitter.opts.sheet, &cfg, &head, &tail);
    }

    /// Random bases between adapters at both ends hold boundary loci within
    /// the chance the anchored budgets of a panel of twenty random primer
    /// pairs admit. Each primer alone keeps its error-rate ceiling within the
    /// chance target; the set-wide bound lowers the panel's budgets.
    #[test]
    fn random_bases_between_adapters_stay_within_the_chance_target_of_a_panel() {
        const PAIRS: u64 = 20;
        const LEN: usize = 15;
        let mut text = String::from("target\tfwd\trev\n");
        for pair in 0..PAIRS {
            let (fwd, rev) = (
                splitmix_dna(5_000 + 2 * pair, LEN),
                splitmix_dna(5_001 + 2 * pair, LEN),
            );
            text.push_str(&format!(
                "T{pair}\t{}\t{}\n",
                String::from_utf8_lossy(&fwd),
                String::from_utf8_lossy(&rev)
            ));
        }
        let sheet = Sheet::parse_tsv(&text).unwrap();
        let seqs: Vec<&[u8]> = sheet.primers.iter().map(|p| p.seq.as_slice()).collect();
        let edit_ceiling = crate::adapter::edit_budget(0.2, LEN);
        assert_eq!(
            crate::adapter::anchored_budgets(&seqs[..1], 0.2),
            [edit_ceiling]
        );
        assert!(
            crate::adapter::anchored_budgets(&seqs, 0.2)
                .iter()
                .all(|&k| k < edit_ceiling)
        );

        let adapter = b"AATGTACTTCGTTCAGTTACGTATTGCT";
        let mut cfg = AdapterConfig {
            adapters: vec![crate::adapter::Adapter {
                name: "adapter".into(),
                seq: adapter.to_vec(),
                role: crate::adapter::Role::Adapter,
            }],
            error_rate: 0.2,
            end_size: 150,
            split: true,
            min_piece: 1,
            candidate_index: std::sync::OnceLock::new(),
            amplicon: false,
            split_of: Vec::new(),
            split_opens: Vec::new(),
        };
        cfg.attach_split(&sheet.primers);
        let head = [splitmix_dna(905, 20), adapter.to_vec()].concat();
        let tail = reverse_complement(&head);

        let at_ceiling = substituted(
            seqs[0],
            &(0..edit_ceiling).map(|i| 2 + 5 * i).collect::<Vec<_>>(),
        );
        let read = [&head[..], &at_ceiling, &splitmix_dna(906, 220), &tail].concat();
        let found = pieces(&read, &cfg);
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(found[0].five, None, "{found:?}");

        assert_boundary_loci_within_chance(&sheet, &cfg, &head, &tail);
    }

    #[test]
    fn discards_follow_options() {
        let sheet = two_target_sheet();
        let cfg = adapter_config(&sheet);
        let mut opts = options(sheet, Require::Either, 2);
        opts.discard_unassigned = true;
        opts.discard_ambiguous = false;
        let splitter = Splitter::new(opts, &cfg).unwrap();

        let assigned = Call::Assigned {
            key: 0,
            target: 0,
            strand: Strand::Plus,
            ends: Ends::Both,
        };
        assert!(!splitter.discards(&assigned));
        assert!(splitter.discards(&Call::Unassigned(Unassigned::NoPrimer)));
        assert!(!splitter.discards(&Call::Ambiguous));

        let sheet = two_target_sheet();
        let cfg = adapter_config(&sheet);
        let mut opts = options(sheet, Require::Either, 2);
        opts.discard_unassigned = false;
        opts.discard_ambiguous = true;
        let splitter = Splitter::new(opts, &cfg).unwrap();
        assert!(!splitter.discards(&assigned));
        assert!(!splitter.discards(&Call::Unassigned(Unassigned::NoPrimer)));
        assert!(splitter.discards(&Call::Ambiguous));
    }
}
