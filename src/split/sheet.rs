//! Primer split sheet: the target model plus its TSV, FASTA and preset
//! parsers.
//!
//! A sheet names the primer targets `--split-by` demultiplexes reads into. Each
//! target is a primer mix (a list of forward primers and a list of reverse
//! primers, either of which may be empty) plus a group and an optional
//! expected length range. Any forward primer of a target pairs with any of
//! its reverse primers. Primers are stored once and shared by sequence, so
//! two targets naming the same primer point at one entry.

use std::collections::{HashMap, HashSet};

use anyhow::{Context, bail};
use clap::ValueEnum;

use crate::adapter::MIN_PATTERN_LEN;
use crate::adapter::search;

/// Reserved target and group key names. A sheet cannot define a target with
/// one of these names, since they are the classifier's own output bins.
pub const RESERVED: [&str; 3] = ["unassigned", "ambiguous", "unclassified"];

/// One primer sequence, stored once and shared by every target that names it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Primer {
    /// Display name.
    pub name: String,
    /// Nucleotide sequence, uppercase, `U` folded to `T`. May contain IUPAC
    /// ambiguity codes.
    pub seq: Vec<u8>,
}

/// One primer target: a name, a reporting group, and its forward and reverse
/// primer lists.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Target {
    /// Target name, unique within the sheet.
    pub name: String,
    /// Reporting group; defaults to the target name.
    pub group: String,
    /// Indices into `Sheet::primers` of the forward primers, each listed
    /// once; empty when the target has none.
    pub fwd: Vec<usize>,
    /// Indices into `Sheet::primers` of the reverse primers, each listed
    /// once; empty when the target has none. Each is stored as ordered, 5'
    /// to 3', as the primer itself is written.
    pub rev: Vec<usize>,
    /// Expected amplicon length range (min, max), inclusive, when the sheet
    /// gives a `min_len` or `max_len` column.
    pub len: Option<(usize, usize)>,
}

/// A parsed primer split sheet: every distinct primer and the targets that
/// reference them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sheet {
    /// Every distinct primer sequence in the sheet.
    pub primers: Vec<Primer>,
    /// Every target, in sheet order.
    pub targets: Vec<Target>,
}

/// Which end(s) of a read a target needs a located primer at to be assigned.
#[derive(Clone, Copy, PartialEq, Eq, Debug, clap::ValueEnum)]
pub enum Require {
    /// Either end alone is enough.
    Either,
    /// Both ends must carry a located primer.
    Both,
    /// The forward end must carry a located primer.
    Fwd,
    /// The reverse end must carry a located primer.
    Rev,
}

impl Require {
    /// The CLI spelling of this rule, as printed in the summary.
    pub fn label(self) -> &'static str {
        match self {
            Require::Either => "either",
            Require::Both => "both",
            Require::Fwd => "fwd",
            Require::Rev => "rev",
        }
    }
}

/// Uppercases a raw sequence and strips whitespace, folding `U` to `T`: RNA
/// primers are written with `U`, DNA reads store `T`.
fn normalize_seq(raw: &[u8]) -> Vec<u8> {
    raw.iter()
        .filter(|b| !b.is_ascii_whitespace())
        .map(|&b| {
            let upper = b.to_ascii_uppercase();
            if upper == b'U' { b'T' } else { upper }
        })
        .collect()
}

/// Accumulates a sheet's deduplicated primers while parsing. A sequence
/// already present returns its existing index instead of a new entry.
#[derive(Default)]
struct PrimerTable {
    primers: Vec<Primer>,
    by_seq: HashMap<Vec<u8>, usize>,
}

impl PrimerTable {
    /// Adds `raw` under `name`, validated against the minimum pattern length
    /// and the nucleotide alphabet, and returns its index.
    fn add(&mut self, name: &str, raw: &[u8]) -> anyhow::Result<usize> {
        let seq = normalize_seq(raw);
        if seq.len() < MIN_PATTERN_LEN {
            bail!(
                "primer {:?} is {} nt, shorter than the minimum of {MIN_PATTERN_LEN} nt",
                name,
                seq.len()
            );
        }
        if let Some(&bad) = seq.iter().find(|&&b| search::iupac_degeneracy(b).is_none()) {
            bail!(
                "primer {:?} contains a non-nucleotide character {:?}",
                name,
                bad as char
            );
        }
        if let Some(&idx) = self.by_seq.get(&seq) {
            return Ok(idx);
        }
        let idx = self.primers.len();
        self.by_seq.insert(seq.clone(), idx);
        self.primers.push(Primer {
            name: name.to_string(),
            seq,
        });
        Ok(idx)
    }
}

/// Adds every sequence of `seqs` to `table` and returns their indices in
/// order, each listed once. A single sequence is named `stem`; several are
/// named `stem1`, `stem2`, ... by position.
fn add_list(table: &mut PrimerTable, stem: &str, seqs: &[&str]) -> anyhow::Result<Vec<usize>> {
    let mut list = Vec::with_capacity(seqs.len());
    for (n, seq) in seqs.iter().enumerate() {
        let name = if seqs.len() == 1 {
            stem.to_string()
        } else {
            format!("{stem}{}", n + 1)
        };
        push_unique(&mut list, table.add(&name, seq.as_bytes())?);
    }
    Ok(list)
}

/// Appends `idx` to `list` unless the list already holds it.
fn push_unique(list: &mut Vec<usize>, idx: usize) {
    if !list.contains(&idx) {
        list.push(idx);
    }
}

/// Splits a TSV primer cell into its comma-separated sequences. An empty
/// cell gives no sequence; an empty element of a list is an error.
fn cell_list<'a>(cell: &'a str, target: &str, column: &str) -> anyhow::Result<Vec<&'a str>> {
    if cell.is_empty() {
        return Ok(Vec::new());
    }
    let seqs: Vec<&str> = cell.split(',').map(str::trim).collect();
    if seqs.iter().any(|seq| seq.is_empty()) {
        bail!("target {target:?}: empty sequence in the {column} list");
    }
    Ok(seqs)
}

/// Whether every byte of `name` is printable ASCII (`[ -~]`, space
/// included), as the `Z` type of a SAM tag value requires.
fn is_printable_ascii(name: &str) -> bool {
    name.bytes().all(|b| (b' '..=b'~').contains(&b))
}

/// Rejects a target or group name that is not printable ASCII
/// (`is_printable_ascii`), a reserved target or group name, and a duplicate
/// target name.
fn check_target_names(targets: &[Target]) -> anyhow::Result<()> {
    let mut seen: HashSet<&str> = HashSet::new();
    for t in targets {
        if !is_printable_ascii(&t.name) {
            bail!(
                "target {:?} holds a character outside printable ASCII; target names fill the \
                 wt:Z tag, which takes printable ASCII only",
                t.name
            );
        }
        if !is_printable_ascii(&t.group) {
            bail!(
                "group {:?} of target {:?} holds a character outside printable ASCII; group \
                 names fill the wt:Z tag, which takes printable ASCII only",
                t.group,
                t.name
            );
        }
        if RESERVED.iter().any(|r| r.eq_ignore_ascii_case(&t.name)) {
            bail!(
                "target {:?} is a reserved name (one of {:?})",
                t.name,
                RESERVED
            );
        }
        if RESERVED.iter().any(|r| r.eq_ignore_ascii_case(&t.group)) {
            bail!(
                "group {:?} of target {:?} is a reserved name (one of {:?})",
                t.group,
                t.name,
                RESERVED
            );
        }
        if !seen.insert(t.name.as_str()) {
            bail!("target {:?} is defined more than once", t.name);
        }
    }
    Ok(())
}

/// A preset target: name, forward primer mix, reverse primer mix. Each
/// primer is named by its catalog entry (`adapter::catalog::CATALOG`), which
/// holds its sequence.
type PresetTarget = (
    &'static str,
    &'static [&'static str],
    &'static [&'static str],
);

/// The primer mixes of the ONT Microbial Amplicon Barcoding kit
/// SQK-MAB114.24, under their ONT names: one 16S mix and one ITS mix, in
/// which any forward primer pairs with any reverse primer.
const MAB114: [PresetTarget; 2] = [
    (
        "16S",
        &["16S_mix_F", "16S_Bor_F", "16S_Chl_F", "16S_Ent_F"],
        &["16S_mix_R", "16S_Bor_R", "16S_Chl_R"],
    ),
    (
        "ITS",
        &["ITS1", "ITS1_Fus", "ITS1_Mal"],
        &["ITS4", "ITS4_Pyt"],
    ),
];

/// Minimum IUPAC edit distance of the shorter sequence aligned within the
/// longer, searched up to `k` edits, or `None` when no alignment scores
/// within that budget.
fn edit_distance(
    searcher: &mut search::AmbiguousSearcher,
    a: &[u8],
    b: &[u8],
    k: usize,
) -> Option<usize> {
    let (pattern, text) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    search::hits(searcher, pattern, text, k)
        .into_iter()
        .map(|hit| hit.cost)
        .min()
}

/// The three `--split-by` forms, with an inline example, for resolution
/// errors.
const SPEC_FORMS: &str = "SPEC is a sheet file, a preset (mab114), or an inline target \
                          NAME:F:SEQ[,SEQ...]:R:SEQ[,SEQ...][,F:...] such as \
                          16S:F:AGAGTTTGATYMTGGCTCAG:R:TACGGYTACCTTGTTACGACTT";

impl Sheet {
    /// Builds a sheet from parsed primers and targets, checking the target
    /// names (`check_target_names`). A sequence may be in both role lists of
    /// one target: the 5' end of such a read holds the sequence as given and
    /// the 3' end its reverse complement, as a barcode flanking both ends
    /// does.
    fn checked(primers: Vec<Primer>, targets: Vec<Target>) -> anyhow::Result<Sheet> {
        check_target_names(&targets)?;
        Ok(Sheet { primers, targets })
    }
}

impl Sheet {
    /// Loads a sheet from one `--split-by` value: an existing file path first
    /// (FASTA when its first non-blank byte is `>`, TSV otherwise), then a
    /// preset token, then the inline form (`parse_inline`).
    pub fn load(spec: &str) -> anyhow::Result<Sheet> {
        let path = std::path::Path::new(spec);
        if path.is_file() {
            let bytes = std::fs::read(path)
                .with_context(|| format!("reading --split-by sheet {}", path.display()))?;
            return if bytes.iter().find(|b| !b.is_ascii_whitespace()) == Some(&b'>') {
                Sheet::parse_fasta(&bytes)
            } else {
                let text = String::from_utf8(bytes).with_context(|| {
                    format!("--split-by sheet {} is not valid UTF-8", path.display())
                })?;
                Sheet::parse_tsv(&text)
            };
        }
        if let Some(sheet) = Sheet::preset(spec) {
            return Ok(sheet);
        }
        if !spec.contains(':') {
            bail!(
                "--split-by {spec:?} is not an existing file, a known preset or an inline \
                 target; {SPEC_FORMS}"
            );
        }
        Sheet::parse_inline(spec).with_context(|| {
            format!(
                "--split-by {spec:?} is not an existing file or a known preset, and not a valid \
                 inline target; {SPEC_FORMS}"
            )
        })
    }

    /// Loads every `--split-by` value (`load`) and merges them into one
    /// sheet (`merge`).
    pub fn load_all(specs: &[String]) -> anyhow::Result<Sheet> {
        let sheets = specs
            .iter()
            .map(|spec| Sheet::load(spec))
            .collect::<anyhow::Result<Vec<_>>>()?;
        Sheet::merge(sheets)
    }

    /// Merges `sheets` into one, in order. Primers are shared by sequence
    /// across sheets as within one, keeping the first sheet's name for a
    /// shared primer; target names must be unique across all sheets. Two
    /// primers of one list that the merged table stores as one entry are
    /// listed once.
    pub fn merge(sheets: Vec<Sheet>) -> anyhow::Result<Sheet> {
        let mut table = PrimerTable::default();
        let mut targets = Vec::new();
        for sheet in sheets {
            let index = sheet
                .primers
                .iter()
                .map(|p| table.add(&p.name, &p.seq))
                .collect::<anyhow::Result<Vec<usize>>>()?;
            let remap = |list: &[usize]| {
                let mut out = Vec::with_capacity(list.len());
                for &i in list {
                    push_unique(&mut out, index[i]);
                }
                out
            };
            targets.extend(sheet.targets.into_iter().map(|t| Target {
                fwd: remap(&t.fwd),
                rev: remap(&t.rev),
                ..t
            }));
        }
        Sheet::checked(table.primers, targets)
    }

    /// Parses an inline target,
    /// `NAME:F:SEQ[,SEQ...]:R:SEQ[,SEQ...][,F:...]`. `NAME` runs to the first
    /// `:`; the rest is one or more primer pairs. A pair is `:`-separated tag
    /// and list fields, `F`, `R`, or both, once each, and a list is one or
    /// more comma-separated sequences; `R` sequences are written 5' to 3' as
    /// ordered. After a comma, a field that is exactly a tag followed by `:`
    /// begins a new pair, and anything else is another sequence of the
    /// current list, so a sequence is never taken for a tag. A pair that
    /// gives one tag followed by a pair that gives only the other is an
    /// error: it is one pair with `,` written for `:`. One pair is target
    /// `NAME`; several are targets `NAME.1`, `NAME.2`, ... in group `NAME`.
    /// A list of one primer is named `<target>_F` or `<target>_R`; a longer
    /// list `<target>_F1`, `<target>_F2`, ... by position.
    pub fn parse_inline(spec: &str) -> anyhow::Result<Sheet> {
        /// Tag spellings by role index: forward, reverse.
        const TAGS: [&str; 2] = ["F", "R"];

        /// The role index a tag field names, or `None` for any other field.
        fn role_of(field: &str) -> Option<usize> {
            TAGS.iter().position(|tag| tag.eq_ignore_ascii_case(field))
        }

        /// The sequences one primer pair gives per role, and whether each
        /// tag has appeared in it.
        #[derive(Default)]
        struct Pair<'a> {
            seqs: [Vec<&'a str>; 2],
            tagged: [bool; 2],
        }

        let (name, rest) = spec.split_once(':').unwrap_or((spec, ""));
        let name = name.trim();
        if name.is_empty() {
            bail!("inline target has an empty target name");
        }
        if rest.trim().is_empty() {
            bail!("inline target {name:?} has no primer; give F:SEQ, R:SEQ or both");
        }

        let mut pairs: Vec<Pair> = Vec::new();
        // Role of the list the latest sequence joined.
        let mut current = 0;
        for chunk in rest.split(',') {
            let fields: Vec<&str> = chunk.split(':').map(str::trim).collect();
            let opens_pair = pairs.is_empty() || (fields.len() > 1 && role_of(fields[0]).is_some());
            if opens_pair {
                pairs.push(Pair::default());
            }
            let n = pairs.len();
            let pair = &mut pairs[n - 1];
            let mut fields = fields.into_iter();
            if !opens_pair {
                let seq = fields.next().unwrap_or("");
                if seq.is_empty() {
                    bail!(
                        "inline target {name:?}: empty sequence in the {} list of primer pair \
                         {n}; a list takes no trailing or doubled comma",
                        TAGS[current]
                    );
                }
                pair.seqs[current].push(seq);
            }
            while let Some(tag) = fields.next() {
                if tag.is_empty() {
                    let what = if fields.len() == 0 {
                        "trailing"
                    } else {
                        "doubled"
                    };
                    bail!("inline target {name:?}: {what} ':' in primer pair {n}");
                }
                let Some(role) = role_of(tag) else {
                    bail!(
                        "inline target {name:?}: unknown tag {tag:?} in primer pair {n}; tags \
                         are F and R"
                    );
                };
                let tag = TAGS[role];
                if pair.tagged[role] {
                    bail!("inline target {name:?}: primer pair {n} gives {tag} more than once");
                }
                let seq = fields.next().unwrap_or("");
                if seq.is_empty() {
                    bail!("inline target {name:?}: tag {tag} has no sequence in primer pair {n}");
                }
                pair.tagged[role] = true;
                pair.seqs[role].push(seq);
                current = role;
            }
        }

        // A pair of one tag followed by a pair of only the other is one pair
        // whose second tag was written after a comma in place of a colon.
        for (n, twin) in pairs.windows(2).enumerate() {
            let only = |pair: &Pair| match pair.tagged {
                [true, false] => Some(0),
                [false, true] => Some(1),
                _ => None,
            };
            if let (Some(first), Some(second)) = (only(&twin[0]), only(&twin[1]))
                && first != second
            {
                let tag = TAGS[second];
                bail!(
                    "inline target {name:?}: primer pair {} gives only {tag} after a comma, and \
                     pair {} only {}; write \":{tag}:\" to add {} primers to pair {}, or give \
                     each target its own --split-by",
                    n + 2,
                    n + 1,
                    TAGS[first],
                    ["forward", "reverse"][second],
                    n + 1,
                );
            }
        }

        let mut table = PrimerTable::default();
        let mut targets = Vec::with_capacity(pairs.len());
        for (n, pair) in pairs.iter().enumerate() {
            let target = if pairs.len() == 1 {
                name.to_string()
            } else {
                format!("{name}.{}", n + 1)
            };
            let fwd = add_list(&mut table, &format!("{target}_F"), &pair.seqs[0])?;
            let rev = add_list(&mut table, &format!("{target}_R"), &pair.seqs[1])?;
            targets.push(Target {
                name: target,
                group: name.to_string(),
                fwd,
                rev,
                len: None,
            });
        }
        Sheet::checked(table.primers, targets)
    }

    /// Parses a tab-separated sheet. The header is matched case-insensitively
    /// and requires `target`, `fwd` and `rev`; `group`, `min_len` and
    /// `max_len` are optional. Lines starting with `#`, and blank lines, are
    /// skipped. A primer cell holds one sequence or a comma-separated list;
    /// either cell may be empty (no primer of that role). A cell of one
    /// primer is named `<target>_fwd` or `<target>_rev`; a list
    /// `<target>_fwd1`, `<target>_fwd2`, ... by position.
    pub fn parse_tsv(text: &str) -> anyhow::Result<Sheet> {
        let mut lines = text
            .lines()
            .filter(|l| !l.trim().is_empty() && !l.trim_start().starts_with('#'));
        let header = lines
            .next()
            .ok_or_else(|| anyhow::anyhow!("primer sheet is empty"))?;
        let columns: Vec<String> = header
            .split('\t')
            .map(|c| c.trim().to_ascii_lowercase())
            .collect();
        let col = |name: &str| -> anyhow::Result<usize> {
            columns.iter().position(|c| c == name).ok_or_else(|| {
                anyhow::anyhow!("primer sheet header is missing the required column {name:?}")
            })
        };
        let target_col = col("target")?;
        let fwd_col = col("fwd")?;
        let rev_col = col("rev")?;
        let group_col = columns.iter().position(|c| c == "group");
        let min_len_col = columns.iter().position(|c| c == "min_len");
        let max_len_col = columns.iter().position(|c| c == "max_len");

        let mut table = PrimerTable::default();
        let mut targets = Vec::new();
        for line in lines {
            let fields: Vec<&str> = line.split('\t').collect();
            let field = |idx: usize| fields.get(idx).map_or("", |s| s.trim());

            let name = field(target_col);
            if name.is_empty() {
                bail!("a primer sheet row has an empty target name");
            }
            let fwd = add_list(
                &mut table,
                &format!("{name}_fwd"),
                &cell_list(field(fwd_col), name, "fwd")?,
            )?;
            let rev = add_list(
                &mut table,
                &format!("{name}_rev"),
                &cell_list(field(rev_col), name, "rev")?,
            )?;
            let group = group_col
                .map(field)
                .filter(|g| !g.is_empty())
                .unwrap_or(name)
                .to_string();
            let len = if min_len_col.is_some() || max_len_col.is_some() {
                let min = min_len_col
                    .map(field)
                    .filter(|s| !s.is_empty())
                    .map(|s| {
                        s.parse::<usize>()
                            .with_context(|| format!("target {name:?}: min_len {s:?}"))
                    })
                    .transpose()?
                    .unwrap_or(0);
                let max = max_len_col
                    .map(field)
                    .filter(|s| !s.is_empty())
                    .map(|s| {
                        s.parse::<usize>()
                            .with_context(|| format!("target {name:?}: max_len {s:?}"))
                    })
                    .transpose()?
                    .unwrap_or(usize::MAX);
                Some((min, max))
            } else {
                None
            };
            targets.push(Target {
                name: name.to_string(),
                group,
                fwd,
                rev,
                len,
            });
        }
        Sheet::checked(table.primers, targets)
    }

    /// Parses an annotated FASTA: entries carrying a `target=` header field
    /// are split primers, named as the record head's first word; `end=fwd` or
    /// `end=rev` places the primer. Several entries may share one `target=`
    /// and `end=`; each adds a primer to that list. `group=` overrides the
    /// default group (the
    /// target name); it may appear on any entry of the target, and two
    /// entries of one target giving different `group=` values is an error.
    /// Entries without a `target=` field are ignored, and belong to the
    /// caller's other adapter-FASTA handling instead.
    pub fn parse_fasta(bytes: &[u8]) -> anyhow::Result<Sheet> {
        use seq_io::fasta::{Reader, Record};

        /// A target under construction: `group` is `None` until some entry
        /// gives an explicit `group=`, and defaults to the target name only
        /// once parsing is done.
        struct Building {
            name: String,
            group: Option<String>,
            fwd: Vec<usize>,
            rev: Vec<usize>,
        }

        let mut reader = Reader::new(bytes);
        let mut table = PrimerTable::default();
        let mut targets: Vec<Building> = Vec::new();
        let mut index: HashMap<String, usize> = HashMap::new();

        while let Some(rec) = reader.next() {
            let rec = rec.map_err(|e| anyhow::anyhow!("primer FASTA: {e}"))?;
            let head = String::from_utf8_lossy(rec.head()).into_owned();
            let (name, description) = head
                .split_once(char::is_whitespace)
                .unwrap_or((head.as_str(), ""));
            let name = name.to_string();

            let mut target_field: Option<&str> = None;
            let mut end_field: Option<&str> = None;
            let mut group_field: Option<&str> = None;
            for tok in description.split_whitespace() {
                if let Some((k, v)) = tok.split_once('=') {
                    match k.to_ascii_lowercase().as_str() {
                        "target" => target_field = Some(v),
                        "end" => end_field = Some(v),
                        "group" => group_field = Some(v),
                        _ => {},
                    }
                }
            }
            let Some(target_name) = target_field else {
                continue;
            };
            let end = end_field.ok_or_else(|| {
                anyhow::anyhow!("primer {name:?}: target={target_name:?} needs end=fwd or end=rev")
            })?;

            let idx = table.add(&name, rec.seq())?;
            let target_idx = *index.entry(target_name.to_string()).or_insert_with(|| {
                targets.push(Building {
                    name: target_name.to_string(),
                    group: None,
                    fwd: Vec::new(),
                    rev: Vec::new(),
                });
                targets.len() - 1
            });
            if let Some(g) = group_field {
                match &targets[target_idx].group {
                    Some(existing) if existing != g => bail!(
                        "primer {name:?}: target={target_name:?} group={g:?} conflicts with the \
                         group {existing:?} already given for this target",
                    ),
                    _ => targets[target_idx].group = Some(g.to_string()),
                }
            }
            match end.to_ascii_lowercase().as_str() {
                "fwd" => push_unique(&mut targets[target_idx].fwd, idx),
                "rev" => push_unique(&mut targets[target_idx].rev, idx),
                other => bail!("primer {name:?}: end={other:?} must be fwd or rev"),
            }
        }

        if targets.is_empty() {
            bail!("primer FASTA has no entry with a target= field");
        }
        let targets: Vec<Target> = targets
            .into_iter()
            .map(|t| Target {
                group: t.group.unwrap_or_else(|| t.name.clone()),
                name: t.name,
                fwd: t.fwd,
                rev: t.rev,
                len: None,
            })
            .collect();
        Sheet::checked(table.primers, targets)
    }

    /// Returns the built-in sheet for `token`, or `None` when it names no
    /// known preset. `mab114` defines targets `16S` and `ITS` from the kit's
    /// primer mixes (`MAB114`), with the sequences of the catalog entries of
    /// those names.
    pub fn preset(token: &str) -> Option<Sheet> {
        if !token.eq_ignore_ascii_case("mab114") {
            return None;
        }
        let mut table = PrimerTable::default();
        let mut mix = |primers: &[&str]| -> Vec<usize> {
            let mut list = Vec::with_capacity(primers.len());
            for name in primers {
                let seq = crate::adapter::catalog::sequence_of(name)
                    .expect("preset primers name catalog entries");
                let idx = table
                    .add(name, seq)
                    .expect("preset primers meet the sheet's own length and alphabet rules");
                push_unique(&mut list, idx);
            }
            list
        };
        let targets = MAB114
            .iter()
            .map(|(name, fwd, rev)| Target {
                name: name.to_string(),
                group: name.to_string(),
                fwd: mix(fwd),
                rev: mix(rev),
                len: None,
            })
            .collect();
        Some(Sheet {
            primers: table.primers,
            targets,
        })
    }

    /// Checks that every target has a primer in each role list `require`
    /// needs. Load-time parsing allows an empty list; this is the separate
    /// check for a specific `--split-require` value.
    pub fn validate(&self, require: Require) -> anyhow::Result<()> {
        for t in &self.targets {
            let ok = match require {
                Require::Both => !t.fwd.is_empty() && !t.rev.is_empty(),
                Require::Either => !t.fwd.is_empty() || !t.rev.is_empty(),
                Require::Fwd => !t.fwd.is_empty(),
                Require::Rev => !t.rev.is_empty(),
            };
            if !ok {
                let spelling = require
                    .to_possible_value()
                    .expect("Require has a possible value for every variant")
                    .get_name()
                    .to_string();
                bail!(
                    "target {:?} has no primer for --split-require {}",
                    t.name,
                    spelling
                );
            }
        }
        Ok(())
    }

    /// Returns every pair of primers of two targets of different keys (target
    /// names by default, or groups when `by_group`) whose IUPAC edit distance
    /// is below `lead`, regardless of which role list each primer is in:
    /// primer index, primer index (lower first), distance. Rescoring aligns
    /// every sheet primer in the orientation valid for each end, so a forward
    /// primer of one key and a reverse primer of another compete at the same
    /// end too; those keys cannot be reliably told apart at `--split-lead
    /// lead`. A pair is reported once even when several role combinations or
    /// target pairs reach the same two primers. Primers of one key, the
    /// variants of one target among them, are never compared.
    pub fn close_pairs(&self, lead: usize, by_group: bool) -> Vec<(usize, usize, usize)> {
        self.close_pairs_by(lead, |i, j| {
            let (a, b) = (&self.targets[i], &self.targets[j]);
            if by_group {
                a.group == b.group
            } else {
                a.name == b.name
            }
        })
    }

    /// `close_pairs` with the key relation given by `same_key`, which tells
    /// whether two target indices share a key.
    pub fn close_pairs_by(
        &self,
        lead: usize,
        same_key: impl Fn(usize, usize) -> bool,
    ) -> Vec<(usize, usize, usize)> {
        fn primers(t: &Target) -> impl Iterator<Item = usize> + '_ {
            t.fwd.iter().chain(&t.rev).copied()
        }

        let mut searcher = search::new_searcher_fwd();
        let mut seen: HashSet<(usize, usize)> = HashSet::new();
        let mut out = Vec::new();
        for i in 0..self.targets.len() {
            for j in (i + 1)..self.targets.len() {
                if same_key(i, j) {
                    continue;
                }
                for a in primers(&self.targets[i]) {
                    for b in primers(&self.targets[j]) {
                        if a == b {
                            continue;
                        }
                        let pair = (a.min(b), a.max(b));
                        if !seen.insert(pair) {
                            continue;
                        }
                        let dist = edit_distance(
                            &mut searcher,
                            &self.primers[pair.0].seq,
                            &self.primers[pair.1].seq,
                            lead,
                        );
                        if let Some(dist) = dist
                            && dist < lead
                        {
                            out.push((pair.0, pair.1, dist));
                        }
                    }
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tsv_three_targets_two_groups() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\n\
             16S_full\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\t16S\n\
             16S_V34\tCCTACGGGNGGCWGCAG\tGACTACHVGGGTATCTAATCC\t16S\n\
             ITS\tCTTGGTCATTTAGAGGAAGTAA\tTCCTCCGCTTATTGATATGC\t\n",
        )
        .unwrap();
        assert_eq!(s.targets.len(), 3);
        assert_eq!(s.targets[2].group, "ITS"); // empty group defaults to target
        assert_eq!(s.primers.len(), 6);
        assert_eq!(s.primers[s.targets[0].rev[0]].seq, b"CGGTTACCTTGTTACGACTT");
    }

    #[test]
    fn tsv_shared_primer_is_one_entry() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\n\
             A\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\t\n\
             B\tAGRGTTYGATYMTGGCTCAG\tTCCTCCGCTTATTGATATGC\t\n",
        )
        .unwrap();
        assert_eq!(s.primers.len(), 3, "the shared fwd primer is one entry");
        assert_eq!(
            s.targets[0].fwd, s.targets[1].fwd,
            "both targets point at the same shared primer"
        );
    }

    #[test]
    fn tsv_optional_length_columns() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\tmin_len\tmax_len\n\
             16S\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\t16S\t300\t1000\n",
        )
        .unwrap();
        assert_eq!(s.targets[0].len, Some((300, 1000)));
    }

    #[test]
    fn tsv_rejects_missing_header_column() {
        let err = Sheet::parse_tsv(
            "target\tfwd\tgroup\n\
             16S\tAGRGTTYGATYMTGGCTCAG\t16S\n",
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("rev"),
            "error names the missing column: {err}"
        );
    }

    #[test]
    fn tsv_rejects_reserved_and_duplicate_names() {
        let err = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             unassigned\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("unassigned"),
            "error names the reserved target: {err}"
        );

        let err = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             16S\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\n\
             16S\tCCTACGGGNGGCWGCAG\tGACTACHVGGGTATCTAATCC\n",
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("16S"),
            "error names the duplicate target: {err}"
        );

        let err = Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\n\
             16S\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\tAmbiguous\n",
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("group \"Ambiguous\""),
            "error names the reserved group: {err}"
        );
    }

    #[test]
    fn names_outside_printable_ascii_are_refused() {
        let err = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             16S_\u{e9}\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("target \"16S_\u{e9}\""), "{err}");
        assert!(err.contains("printable ASCII"), "{err}");

        let err = Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\n\
             16S\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\tg\u{1}\n",
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("group \"g\\u{1}\""), "{err}");

        let err = Sheet::parse_fasta(
            ">27F primer target=16S\u{b5} end=fwd\nAGRGTTYGATYMTGGCTCAG\n".as_bytes(),
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("printable ASCII"), "{err}");

        let s = Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\n\
             16S full ~1\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\tmarker genes\n",
        )
        .unwrap();
        assert_eq!(
            (s.targets[0].name.as_str(), s.targets[0].group.as_str()),
            ("16S full ~1", "marker genes")
        );
    }

    #[test]
    fn tsv_rejects_short_and_non_nucleotide_primer() {
        let err = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             A\tACGTACGTAC\tCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("11"),
            "error names the minimum length: {err}"
        );

        let err = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             A\tACGTXACGTACGT\tCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("non-nucleotide"),
            "error names the alphabet problem: {err}"
        );
        assert!(
            err.to_string().contains('X'),
            "error names the offending character: {err}"
        );
    }

    /// Uppercasing and the `U` to `T` fold apply to every primer the sheet
    /// stores, not only ones already written in FASTA/DNA convention.
    #[test]
    fn lowercase_primer_with_u_is_stored_uppercase_with_t() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             A\tacguacguacgu\tCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap();
        assert_eq!(s.primers[s.targets[0].fwd[0]].seq, b"ACGTACGTACGT");
    }

    #[test]
    fn fasta_headers_define_targets() {
        let s = Sheet::parse_fasta(
            b">27F primer target=16S end=fwd\nAGRGTTYGATYMTGGCTCAG\n\
              >1492R primer target=16S end=rev\nCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap();
        assert_eq!(s.targets.len(), 1);
        assert_eq!((s.targets[0].fwd.len(), s.targets[0].rev.len()), (1, 1));
        assert_eq!(s.primers[0].name, "27F");
    }

    #[test]
    fn fasta_without_target_key_is_error() {
        let err = Sheet::parse_fasta(
            b">27F primer\nAGRGTTYGATYMTGGCTCAG\n\
              >1492R primer\nCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("target"),
            "error names the missing key: {err}"
        );
    }

    /// `group=` may appear on any entry of a target, not only the first one
    /// parsed.
    #[test]
    fn fasta_group_field_on_any_entry_sets_the_target_group() {
        let s = Sheet::parse_fasta(
            b">27F primer target=16S end=fwd\nAGRGTTYGATYMTGGCTCAG\n\
              >1492R primer target=16S end=rev group=bacteria\nCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap();
        assert_eq!(s.targets.len(), 1);
        assert_eq!(s.targets[0].group, "bacteria");
    }

    #[test]
    fn fasta_conflicting_group_fields_is_error() {
        let err = Sheet::parse_fasta(
            b">27F primer target=16S end=fwd group=A\nAGRGTTYGATYMTGGCTCAG\n\
              >1492R primer target=16S end=rev group=B\nCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap_err();
        assert!(
            err.to_string().contains("group"),
            "error names the conflicting field: {err}"
        );
    }

    #[test]
    fn preset_mab114_holds_the_kit_primer_mixes() {
        let s = Sheet::preset("mab114").unwrap();
        let names: Vec<&str> = s.targets.iter().map(|t| t.name.as_str()).collect();
        assert_eq!(names, ["16S", "ITS"]);
        assert_eq!(s.primers.len(), 12);
        assert_eq!(
            primers_of(&s, &s.targets[0].fwd),
            [
                ("16S_mix_F", "AGRGTTYGATYMTGGCTCAG"),
                ("16S_Bor_F", "AGAGTTTGATCCTGGCTTAG"),
                ("16S_Chl_F", "AGAATTTGATCTTRGTTCAG"),
                ("16S_Ent_F", "AGAGTTTGATCATGGCTCAG"),
            ]
        );
        assert_eq!(
            primers_of(&s, &s.targets[0].rev),
            [
                ("16S_mix_R", "SGGYTACCTTGTTACGACTT"),
                ("16S_Bor_R", "CGGCTACCTTGTTACGACTT"),
                ("16S_Chl_R", "GGGCTACCTTGTTACGACTT"),
            ]
        );
        assert_eq!(
            primers_of(&s, &s.targets[1].fwd),
            [
                ("ITS1", "TCCGTAGGTGAACCTGCGG"),
                ("ITS1_Fus", "TCCGTTGGTGAACCAGCGG"),
                ("ITS1_Mal", "TCTGTAGGTGAACCTGCAG"),
            ]
        );
        assert_eq!(
            primers_of(&s, &s.targets[1].rev),
            [
                ("ITS4", "TCCTCCGCTTATTGATATGC"),
                ("ITS4_Pyt", "TCCTCCGCTTATTAATATGC"),
            ]
        );
        assert!(s.validate(Require::Both).is_ok());
    }

    #[test]
    fn tsv_cells_hold_comma_separated_lists() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             16S\tAGRGTTYGATYMTGGCTCAG, AGAGTTTGATCCTGGCTTAG\tSGGYTACCTTGTTACGACTT\n\
             ITS\t\tTCCTCCGCTTATTGATATGC,TCCTCCGCTTATTAATATGC\n",
        )
        .unwrap();
        assert_eq!(
            primers_of(&s, &s.targets[0].fwd),
            [
                ("16S_fwd1", "AGRGTTYGATYMTGGCTCAG"),
                ("16S_fwd2", "AGAGTTTGATCCTGGCTTAG"),
            ]
        );
        assert_eq!(
            primers_of(&s, &s.targets[0].rev),
            [("16S_rev", "SGGYTACCTTGTTACGACTT")]
        );
        assert!(s.targets[1].fwd.is_empty());
        assert_eq!(
            primers_of(&s, &s.targets[1].rev),
            [
                ("ITS_rev1", "TCCTCCGCTTATTGATATGC"),
                ("ITS_rev2", "TCCTCCGCTTATTAATATGC"),
            ]
        );
    }

    #[test]
    fn tsv_list_with_an_empty_sequence_is_error() {
        for cell in [
            "AGRGTTYGATYMTGGCTCAG,",
            ",AGRGTTYGATYMTGGCTCAG",
            "AGRGTTYGATYMTGGCTCAG,,AGAGTTTGATCCTGGCTTAG",
        ] {
            let err = Sheet::parse_tsv(&format!("target\tfwd\trev\nA\t{cell}\t\n"))
                .unwrap_err()
                .to_string();
            assert!(err.contains("empty sequence"), "{cell}: {err}");
            assert!(err.contains("fwd"), "{cell}: {err}");
        }
    }

    #[test]
    fn fasta_repeated_ends_add_to_the_list() {
        let s = Sheet::parse_fasta(
            b">16S_mix_F primer target=16S end=fwd\nAGRGTTYGATYMTGGCTCAG\n\
              >16S_mix_R primer target=16S end=rev\nSGGYTACCTTGTTACGACTT\n\
              >16S_Bor_F primer target=16S end=fwd\nAGAGTTTGATCCTGGCTTAG\n\
              >16S_Bor_R primer target=16S end=rev\nCGGCTACCTTGTTACGACTT\n\
              >16S_Chl_F primer target=16S end=fwd\nAGAATTTGATCTTRGTTCAG\n",
        )
        .unwrap();
        assert_eq!(s.targets.len(), 1);
        let names = |list: &[usize]| -> Vec<&str> {
            primers_of(&s, list).into_iter().map(|p| p.0).collect()
        };
        assert_eq!(
            names(&s.targets[0].fwd),
            ["16S_mix_F", "16S_Bor_F", "16S_Chl_F"]
        );
        assert_eq!(names(&s.targets[0].rev), ["16S_mix_R", "16S_Bor_R"]);
    }

    /// A sequence in both role lists of one target loads in every sheet
    /// form and is stored once, shared by the two lists. A reverse primer
    /// that is the reverse complement of the forward primer is a different
    /// sequence.
    #[test]
    fn one_sequence_in_both_role_lists_of_a_target_loads() {
        let seq = "AGAGTTTGATCCTGGCTTAG";
        let inline = Sheet::parse_inline(&format!("A:F:{seq}:R:{seq}")).unwrap();
        let tsv = Sheet::parse_tsv(&format!("target\tfwd\trev\nA\t{seq}\t{seq}\n")).unwrap();
        let fasta = Sheet::parse_fasta(
            format!(">a_f target=A end=fwd\n{seq}\n>a_r target=A end=rev\n{seq}\n").as_bytes(),
        )
        .unwrap();
        for (form, s) in [("inline", inline), ("tsv", tsv), ("fasta", fasta)] {
            assert_eq!(s.primers.len(), 1, "{form}");
            assert_eq!(s.primers[0].seq, seq.as_bytes(), "{form}");
            assert_eq!(s.targets[0].fwd, [0], "{form}");
            assert_eq!(s.targets[0].rev, [0], "{form}");
            for require in [Require::Either, Require::Both, Require::Fwd, Require::Rev] {
                assert!(s.validate(require).is_ok(), "{form} {require:?}");
            }
        }

        let s = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             A\tAGRGTTYGATYMTGGCTCAG,AGAGTTTGATCCTGGCTTAG\tagagtttgatcctggcttag\n",
        )
        .unwrap();
        assert_eq!(s.primers.len(), 2);
        assert_eq!(s.targets[0].fwd, [0, 1]);
        assert_eq!(s.targets[0].rev, [1]);

        let s = Sheet::parse_inline("A:F:AGAGTTTGATCCTGGCTTAG:R:CTAAGCCAGGATCAAACTCT").unwrap();
        assert_eq!(s.primers.len(), 2);
        let s = Sheet::parse_inline("A:F:ACGTACGTACGT:R:ACGTACGTACGT").unwrap();
        assert_eq!(s.targets[0].fwd, s.targets[0].rev);

        let s = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             A\tAGAGTTTGATCCTGGCTTAG\tCGGTTACCTTGTTACGACTT\n\
             B\tTCCTCCGCTTATTGATATGC\tAGAGTTTGATCCTGGCTTAG\n",
        )
        .unwrap();
        assert_eq!(
            s.targets[0].fwd, s.targets[1].rev,
            "a sequence may take different roles in different targets"
        );
    }

    #[test]
    fn validate_checks_the_role_lists_the_rule_needs() {
        let fwd_only =
            Sheet::parse_inline("A:F:AGRGTTYGATYMTGGCTCAG,AGAGTTTGATCCTGGCTTAG").unwrap();
        assert!(fwd_only.validate(Require::Either).is_ok());
        assert!(fwd_only.validate(Require::Fwd).is_ok());
        assert!(fwd_only.validate(Require::Rev).is_err());
        assert!(fwd_only.validate(Require::Both).is_err());
        let mix = Sheet::parse_inline(
            "A:F:AGRGTTYGATYMTGGCTCAG,AGAGTTTGATCCTGGCTTAG:R:SGGYTACCTTGTTACGACTT",
        )
        .unwrap();
        for require in [Require::Either, Require::Both, Require::Fwd, Require::Rev] {
            assert!(mix.validate(require).is_ok(), "{require:?}");
        }
    }

    /// Variants of one target are never compared with each other; every
    /// variant of one key is compared with every variant of another.
    #[test]
    fn close_pairs_compares_lists_of_different_keys_only() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             A\tAGAGTTTGATCCTGGCTTAG,AGAGTTTGATCCTGGCTTAC\tCGGTTACCTTGTTACGACTT\n\
             B\tTCCGTAGGTGAACCTGCGG,AGAGTTTGATCCTGGCTTAT\tTCCTCCGCTTATTGATATGC\n",
        )
        .unwrap();
        let pairs = s.close_pairs(2, false);
        let flagged: HashSet<(usize, usize)> = pairs.iter().map(|p| (p.0, p.1)).collect();
        let (a, b) = (&s.targets[0].fwd, &s.targets[1].fwd);
        let expected: HashSet<(usize, usize)> = [
            (a[0].min(b[1]), a[0].max(b[1])),
            (a[1].min(b[1]), a[1].max(b[1])),
        ]
        .into();
        assert_eq!(flagged, expected);
    }

    #[test]
    fn validate_both_requires_rev() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             A\tAGRGTTYGATYMTGGCTCAG\t\n",
        )
        .unwrap();
        assert!(s.validate(Require::Both).is_err());
        assert!(s.validate(Require::Either).is_ok());
    }

    /// The error names the rule in its CLI spelling (`both`), not the
    /// `Debug` form (`Both`), so it reads the way the user typed the flag.
    #[test]
    fn validate_error_uses_the_cli_spelling() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\n\
             A\tAGRGTTYGATYMTGGCTCAG\t\n",
        )
        .unwrap();
        let err = s.validate(Require::Both).unwrap_err();
        assert!(err.to_string().contains("both"), "{err}");
        assert!(!err.to_string().contains("Both"), "{err}");
    }

    #[test]
    fn close_pairs_flags_one_edit_primers() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\n\
             A\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\tG\n\
             B\tAGRGTTYGATYMTGGCTCAC\tTCCTCCGCTTATTGATATGC\tG\n",
        )
        .unwrap();

        let pairs = s.close_pairs(2, false);
        assert_eq!(pairs.len(), 1, "different targets: one close pair");
        assert_eq!(pairs[0].2, 1, "one substitution apart");
        let flagged: std::collections::HashSet<usize> = [pairs[0].0, pairs[0].1].into();
        let expected: std::collections::HashSet<usize> =
            [s.targets[0].fwd[0], s.targets[1].fwd[0]].into();
        assert_eq!(flagged, expected);

        let grouped = s.close_pairs(2, true);
        assert!(
            grouped.is_empty(),
            "same group under by_group: no pair reported"
        );
    }

    /// A forward primer of one key and a reverse primer of another key
    /// compete at the same end too (rescoring aligns every sheet primer at
    /// every end), so `close_pairs` must flag a close pair across roles, not
    /// only fwd-vs-fwd or rev-vs-rev.
    #[test]
    fn close_pairs_flags_cross_role_pairs() {
        let s = Sheet::parse_tsv(
            "target\tfwd\trev\tgroup\n\
             A\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\tGA\n\
             B\t\tAGRGTTYGATYMTGGCTCAC\tGB\n",
        )
        .unwrap();

        let pairs = s.close_pairs(2, false);
        assert_eq!(pairs.len(), 1, "one cross-role pair: {pairs:?}");
        assert_eq!(pairs[0].2, 1, "one substitution apart");
        let flagged: std::collections::HashSet<usize> = [pairs[0].0, pairs[0].1].into();
        let expected: std::collections::HashSet<usize> =
            [s.targets[0].fwd[0], s.targets[1].rev[0]].into();
        assert_eq!(flagged, expected);
    }

    #[test]
    fn load_prefers_existing_file_then_preset() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mab114");
        std::fs::write(
            &path,
            "target\tfwd\trev\n16S\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\n",
        )
        .unwrap();

        let s = Sheet::load(path.to_str().unwrap()).unwrap();
        assert_eq!(s.targets.len(), 1);
        assert_eq!(s.targets[0].name, "16S");

        let s = Sheet::load("mab114").unwrap();
        assert_eq!(s.targets.len(), 2, "falls back to the mab114 preset");
    }

    /// Returns the (name, sequence) of the only primer of `list`.
    fn primer<'a>(s: &'a Sheet, list: &[usize]) -> (&'a str, &'a [u8]) {
        assert_eq!(list.len(), 1, "the list holds one primer");
        let p = &s.primers[list[0]];
        (p.name.as_str(), p.seq.as_slice())
    }

    /// Returns the (name, sequence) of every primer of `list`, in order.
    fn primers_of<'a>(s: &'a Sheet, list: &[usize]) -> Vec<(&'a str, &'a str)> {
        list.iter()
            .map(|&i| {
                let p = &s.primers[i];
                (p.name.as_str(), std::str::from_utf8(&p.seq).unwrap())
            })
            .collect()
    }

    #[test]
    fn inline_single_pair() {
        let s = Sheet::parse_inline("16S:F:AGRGTTYGATYMTGGCTCAG:R:CGGTTACCTTGTTACGACTT").unwrap();
        assert_eq!(s.targets.len(), 1);
        let t = &s.targets[0];
        assert_eq!((t.name.as_str(), t.group.as_str()), ("16S", "16S"));
        assert_eq!(t.len, None);
        assert_eq!(primer(&s, &t.fwd), ("16S_F", &b"AGRGTTYGATYMTGGCTCAG"[..]));
        assert_eq!(primer(&s, &t.rev), ("16S_R", &b"CGGTTACCTTGTTACGACTT"[..]));
    }

    #[test]
    fn inline_pool_numbers_targets_in_one_group() {
        let s = Sheet::parse_inline(
            "16S:F:AGRGTTYGATYMTGGCTCAG:R:CGGTTACCTTGTTACGACTT,\
             F:CCTACGGGNGGCWGCAG:R:GACTACHVGGGTATCTAATCC",
        )
        .unwrap();
        let names: Vec<(&str, &str)> = s
            .targets
            .iter()
            .map(|t| (t.name.as_str(), t.group.as_str()))
            .collect();
        assert_eq!(names, [("16S.1", "16S"), ("16S.2", "16S")]);
        assert_eq!(primer(&s, &s.targets[0].fwd).0, "16S.1_F");
        assert_eq!(primer(&s, &s.targets[1].rev).0, "16S.2_R");
        assert_eq!(primer(&s, &s.targets[1].fwd).1, b"CCTACGGGNGGCWGCAG");
    }

    #[test]
    fn inline_forward_only_and_reverse_only_targets() {
        let s = Sheet::parse_inline("ITS:F:CTTGGTCATTTAGAGGAAGTAA").unwrap();
        assert_eq!(s.targets[0].name, "ITS");
        assert_eq!(primer(&s, &s.targets[0].fwd).0, "ITS_F");
        assert!(s.targets[0].rev.is_empty());

        let s = Sheet::parse_inline("ITS:R:TCCTCCGCTTATTGATATGC").unwrap();
        assert!(s.targets[0].fwd.is_empty());
        assert_eq!(primer(&s, &s.targets[0].rev).0, "ITS_R");

        let s = Sheet::parse_inline("ITS:R:TCCTCCGCTTATTGATATGC:F:CTTGGTCATTTAGAGGAAGTAA").unwrap();
        assert_eq!(
            primer(&s, &s.targets[0].fwd).1,
            b"CTTGGTCATTTAGAGGAAGTAA",
            "tags place primers in either order"
        );
    }

    #[test]
    fn inline_normalizes_case_and_u() {
        let s = Sheet::parse_inline("A:f:acguacguacgu:r:cggttaccuugttacgacuu").unwrap();
        assert_eq!(primer(&s, &s.targets[0].fwd).1, b"ACGTACGTACGT");
        assert_eq!(primer(&s, &s.targets[0].rev).1, b"CGGTTACCTTGTTACGACTT");
    }

    /// The tag is the first field of each tag and sequence pair, so a
    /// sequence of IUPAC `R` codes, or one starting with `R` or `F`-like
    /// bases, is never read as a tag.
    #[test]
    fn inline_reverse_sequence_holding_iupac_r_parses_by_position() {
        let s = Sheet::parse_inline("X:R:RRGTACGTACGR:F:RACGTACGTACG").unwrap();
        assert_eq!(primer(&s, &s.targets[0].rev).1, b"RRGTACGTACGR");
        assert_eq!(primer(&s, &s.targets[0].fwd).1, b"RACGTACGTACG");

        let s = Sheet::parse_inline("X:R:RRRRRRRRRRRR").unwrap();
        assert_eq!(primer(&s, &s.targets[0].rev).1, b"RRRRRRRRRRRR");
        assert!(s.targets[0].fwd.is_empty());
    }

    #[test]
    fn inline_lists_give_one_target_several_primers() {
        let s = Sheet::parse_inline(
            "16S:F:AGRGTTYGATYMTGGCTCAG,AGAGTTTGATCCTGGCTTAG,AGAATTTGATCTTRGTTCAG\
             :R:SGGYTACCTTGTTACGACTT,CGGCTACCTTGTTACGACTT",
        )
        .unwrap();
        assert_eq!(s.targets.len(), 1);
        let t = &s.targets[0];
        assert_eq!((t.name.as_str(), t.group.as_str()), ("16S", "16S"));
        assert_eq!(
            primers_of(&s, &t.fwd),
            [
                ("16S_F1", "AGRGTTYGATYMTGGCTCAG"),
                ("16S_F2", "AGAGTTTGATCCTGGCTTAG"),
                ("16S_F3", "AGAATTTGATCTTRGTTCAG"),
            ]
        );
        assert_eq!(
            primers_of(&s, &t.rev),
            [
                ("16S_R1", "SGGYTACCTTGTTACGACTT"),
                ("16S_R2", "CGGCTACCTTGTTACGACTT"),
            ]
        );
    }

    #[test]
    fn inline_pool_pairs_hold_lists() {
        let s = Sheet::parse_inline(
            "16S:F:AGRGTTYGATYMTGGCTCAG,AGAGTTTGATCCTGGCTTAG:R:SGGYTACCTTGTTACGACTT,\
             F:CCTACGGGNGGCWGCAG:R:GACTACHVGGGTATCTAATCC",
        )
        .unwrap();
        let names: Vec<(&str, &str)> = s
            .targets
            .iter()
            .map(|t| (t.name.as_str(), t.group.as_str()))
            .collect();
        assert_eq!(names, [("16S.1", "16S"), ("16S.2", "16S")]);
        assert_eq!(
            primers_of(&s, &s.targets[0].fwd),
            [
                ("16S.1_F1", "AGRGTTYGATYMTGGCTCAG"),
                ("16S.1_F2", "AGAGTTTGATCCTGGCTTAG"),
            ]
        );
        assert_eq!(
            primers_of(&s, &s.targets[0].rev),
            [("16S.1_R", "SGGYTACCTTGTTACGACTT")]
        );
        assert_eq!(
            primers_of(&s, &s.targets[1].fwd),
            [("16S.2_F", "CCTACGGGNGGCWGCAG")]
        );
        assert_eq!(
            primers_of(&s, &s.targets[1].rev),
            [("16S.2_R", "GACTACHVGGGTATCTAATCC")]
        );
    }

    /// After a comma, only a field that is exactly a tag followed by `:`
    /// begins a new pair; a sequence starting with `R` or `F`-like bases, or
    /// holding the IUPAC code `R`, stays in the current list.
    #[test]
    fn inline_list_sequence_starting_with_a_tag_letter_is_a_sequence() {
        let s = Sheet::parse_inline("16S:F:AGAGTTTGATCATGGCTCAG,RGGTTACCTTGTTACGACTT").unwrap();
        assert_eq!(s.targets.len(), 1);
        assert_eq!(
            primers_of(&s, &s.targets[0].fwd),
            [
                ("16S_F1", "AGAGTTTGATCATGGCTCAG"),
                ("16S_F2", "RGGTTACCTTGTTACGACTT"),
            ]
        );
        assert!(s.targets[0].rev.is_empty());

        let s = Sheet::parse_inline("X:R:RRGTACGTACGR,RRRRRRRRRRRR,r:RACGTACGTACG").unwrap();
        assert_eq!(s.targets.len(), 2);
        assert_eq!(
            primers_of(&s, &s.targets[0].rev),
            [("X.1_R1", "RRGTACGTACGR"), ("X.1_R2", "RRRRRRRRRRRR")]
        );
        assert_eq!(
            primers_of(&s, &s.targets[1].rev),
            [("X.2_R", "RACGTACGTACG")]
        );
    }

    /// A sequence repeated within one list is stored and listed once.
    #[test]
    fn inline_list_repeating_a_sequence_lists_it_once() {
        let s = Sheet::parse_inline("A:F:AGRGTTYGATYMTGGCTCAG,agrgttygatymtggctcag").unwrap();
        assert_eq!(s.targets[0].fwd.len(), 1);
    }

    #[test]
    fn inline_malformed_fields_are_errors() {
        const F: &str = "AGRGTTYGATYMTGGCTCAG";
        const R: &str = "CGGTTACCTTGTTACGACTT";
        for (spec, needle) in [
            (format!(":F:{F}"), "empty target name"),
            ("16S".to_string(), "no primer"),
            ("16S:".to_string(), "no primer"),
            (format!("16S:X:{F}"), "\"X\""),
            (format!("16S:F:{F}:Q:{R}"), "\"Q\""),
            (format!("16S:F:{F}:R"), "R has no sequence"),
            (format!("16S:F::R:{R}"), "F has no sequence"),
            (format!("16S:F:{F}:F:{R}"), "F more than once"),
            (format!("16S:R:{R}:R:{F}"), "R more than once"),
            (
                format!("16S:F:{F},"),
                "empty sequence in the F list of primer pair 1",
            ),
            (format!("16S:F:{F},,R:{R}"), "empty sequence in the F list"),
            (format!("16S:F:{F},{R},"), "empty sequence in the F list"),
            (
                format!("16S:F:{F}:R:{R},"),
                "empty sequence in the R list of primer pair 1",
            ),
            (
                format!("16S:F:{F}:R:{R},F:{R},,{F}"),
                "F list of primer pair 2",
            ),
            (format!("16S:F:{F},:R:{R}"), "empty sequence in the F list"),
            (format!("16S:F:{F},{R}:F:{R}"), "F more than once"),
            (
                format!("16S:F:{F},F:"),
                "F has no sequence in primer pair 2",
            ),
            (format!("16S:F:{F}:"), "trailing ':' in primer pair 1"),
            (format!("16S:F:{F}:R:{R}:"), "trailing ':' in primer pair 1"),
            (format!("16S:F:{F}::R:{R}"), "doubled ':' in primer pair 1"),
            (format!("16S:F:{F},ACGT"), "11"),
            (format!("16S:F:{F}:R:{F}GG{R}X"), "non-nucleotide"),
            ("16S:F:ACGTACGTAC".to_string(), "11"),
            (format!("unassigned:F:{F}"), "reserved"),
        ] {
            let err = format!("{:#}", Sheet::parse_inline(&spec).unwrap_err());
            assert!(err.contains(needle), "{spec}: {err}");
        }
    }

    /// A pool in which a pair gives only one tag and the next pair only the
    /// other is one pair with a comma in place of a colon, and is rejected
    /// with the fix.
    #[test]
    fn inline_pool_of_complementary_one_tag_pairs_is_an_error() {
        const F: &str = "AGRGTTYGATYMTGGCTCAG";
        const R: &str = "CGGTTACCTTGTTACGACTT";
        const P: &str = "CCTACGGGNGGCWGCAG";
        for (spec, needle) in [
            (
                format!("16S:F:{F},R:{R}"),
                "primer pair 2 gives only R after a comma, and pair 1 only F; write \":R:\" \
                 to add reverse primers to pair 1",
            ),
            (
                format!("16S:R:{R},F:{F}"),
                "write \":F:\" to add forward primers to pair 1",
            ),
            (format!("16S:F:{F},{P},R:{R}"), "primer pair 2 gives only R"),
            (
                format!("16S:F:{F}:R:{R},F:{P},R:{R}"),
                "primer pair 3 gives only R after a comma, and pair 2 only F",
            ),
        ] {
            let err = format!("{:#}", Sheet::parse_inline(&spec).unwrap_err());
            assert!(err.contains(needle), "{spec}: {err}");
        }
    }

    /// Pools whose pairs all give both tags, pools of one-tag pairs of one
    /// role, and single-pair targets of one tag stay valid.
    #[test]
    fn inline_pools_of_full_or_like_pairs_stay_valid() {
        const F: &str = "AGRGTTYGATYMTGGCTCAG";
        const R: &str = "CGGTTACCTTGTTACGACTT";
        const P: &str = "CCTACGGGNGGCWGCAG";
        for (spec, targets) in [
            (format!("16S:F:{F}:R:{R},F:{P}:R:{R}"), 2),
            (format!("16S:R:{R}:F:{F},R:{R}:F:{P}"), 2),
            (format!("16S:F:{F},F:{P}"), 2),
            (format!("16S:R:{R},R:{P}"), 2),
            (format!("16S:F:{F}:R:{R},R:{P}"), 2),
            (format!("16S:F:{F},{P}"), 1),
            (format!("16S:F:{F}"), 1),
            (format!("16S:R:{R}"), 1),
        ] {
            let s = Sheet::parse_inline(&spec).unwrap_or_else(|e| panic!("{spec}: {e:#}"));
            assert_eq!(s.targets.len(), targets, "{spec}");
        }
    }

    #[test]
    fn load_resolves_inline_last_and_names_every_form_on_failure() {
        let s = Sheet::load("16S:F:AGRGTTYGATYMTGGCTCAG").unwrap();
        assert_eq!(s.targets[0].name, "16S");

        let err = format!("{:#}", Sheet::load("missing_sheet.tsv").unwrap_err());
        for needle in [
            "missing_sheet.tsv",
            "file",
            "mab114",
            "NAME:F:SEQ[,SEQ...]:R:SEQ[,SEQ...]",
            "16S:F:",
        ] {
            assert!(err.contains(needle), "{needle}: {err}");
        }
        let err = format!("{:#}", Sheet::load("16S:X:ACGT").unwrap_err());
        for needle in ["\"X\"", "file", "mab114", "NAME:F:SEQ[,SEQ...]"] {
            assert!(err.contains(needle), "{needle}: {err}");
        }
    }

    #[test]
    fn load_all_merges_a_sheet_file_and_an_inline_target() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sheet.tsv");
        std::fs::write(
            &path,
            "target\tfwd\trev\tmin_len\n16S\tAGRGTTYGATYMTGGCTCAG\tCGGTTACCTTGTTACGACTT\t300\n",
        )
        .unwrap();
        let s = Sheet::load_all(&[
            path.to_str().unwrap().to_string(),
            "ITS:F:CTTGGTCATTTAGAGGAAGTAA:R:TCCTCCGCTTATTGATATGC".to_string(),
        ])
        .unwrap();
        let names: Vec<&str> = s.targets.iter().map(|t| t.name.as_str()).collect();
        assert_eq!(names, ["16S", "ITS"]);
        assert_eq!(s.primers.len(), 4);
        assert_eq!(s.targets[0].len, Some((300, usize::MAX)));
        assert_eq!(primer(&s, &s.targets[0].rev).1, b"CGGTTACCTTGTTACGACTT");
        assert_eq!(
            primer(&s, &s.targets[1].fwd),
            ("ITS_F", &b"CTTGGTCATTTAGAGGAAGTAA"[..])
        );
        assert_eq!(primer(&s, &s.targets[1].rev).1, b"TCCTCCGCTTATTGATATGC");
    }

    #[test]
    fn load_all_merges_a_preset_and_an_inline_pool() {
        let s = Sheet::load_all(&[
            "mab114".to_string(),
            "V34:F:CCTACGGGNGGCWGCAG:R:GACTACHVGGGTATCTAATCC,F:GTGCCAGCMGCCGCGGTAA".to_string(),
        ])
        .unwrap();
        let names: Vec<(&str, &str)> = s
            .targets
            .iter()
            .map(|t| (t.name.as_str(), t.group.as_str()))
            .collect();
        assert_eq!(
            names,
            [
                ("16S", "16S"),
                ("ITS", "ITS"),
                ("V34.1", "V34"),
                ("V34.2", "V34")
            ]
        );
        assert_eq!(s.primers.len(), 15);
        assert_eq!(primer(&s, &s.targets[3].fwd).1, b"GTGCCAGCMGCCGCGGTAA");
        assert!(s.targets[3].rev.is_empty());
    }

    #[test]
    fn load_all_refuses_a_target_name_given_by_two_sources() {
        let err = Sheet::load_all(&[
            "mab114".to_string(),
            "ITS:F:CTTGGTCATTTAGAGGAAGTAC".to_string(),
        ])
        .unwrap_err()
        .to_string();
        assert!(err.contains("\"ITS\""), "{err}");
        assert!(err.contains("more than once"), "{err}");
    }

    #[test]
    fn load_all_shares_a_primer_named_by_two_sources() {
        let s = Sheet::load_all(&[
            "A:F:AGRGTTYGATYMTGGCTCAG:R:CGGTTACCTTGTTACGACTT".to_string(),
            "B:F:agrgttygatymtggctcag:R:TCCTCCGCTTATTGATATGC".to_string(),
        ])
        .unwrap();
        assert_eq!(s.primers.len(), 3, "the shared forward primer is one entry");
        assert_eq!(s.targets[0].fwd, s.targets[1].fwd);
        assert_eq!(primer(&s, &s.targets[1].fwd).0, "A_F");
        assert_eq!(primer(&s, &s.targets[1].rev).0, "B_R");
    }
}
