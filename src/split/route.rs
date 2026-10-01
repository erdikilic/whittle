//! Output routing for `--split-by`: the `-o` path template, the table of
//! expanded output paths the render workers route segments to, and the
//! checks each path passes before its file is created.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::RwLock;

use anyhow::Context as _;

use super::classify::{KeyLevel, Keys};
use super::sheet::Sheet;
use crate::workflow::KeyId;

/// The `{barcode}` value of a record without a `BC:Z` barcode call.
pub const UNCLASSIFIED: &str = "unclassified";

/// Opened output files above which a run is warned, once, about the
/// open-file limit.
pub(crate) const OPEN_FILES_WARNING: usize = 512;

/// A `-o` path holding at least one of the placeholders `{target}`,
/// `{group}` and `{barcode}`; each split key is written to the path the
/// template expands to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Template {
    /// The path as given.
    raw: String,
    /// The key the path names: `Target` when `{target}` is present,
    /// otherwise `Group`.
    key: KeyLevel,
    /// Whether the path holds `{barcode}`.
    barcode: bool,
}

/// One piece of a template: literal text, or the name inside a `{...}`
/// placeholder.
#[derive(Clone, Copy)]
enum Token<'a> {
    Text(&'a str),
    Field(&'a str),
}

/// Splits `raw` into tokens, in order. A `{` opens a placeholder when a `}`
/// follows it before any other `{` or `/`; otherwise it is text.
fn scan<'a>(raw: &'a str, mut on: impl FnMut(Token<'a>)) {
    let mut rest = raw;
    while let Some(open) = rest.find('{') {
        let after = &rest[open + 1..];
        match after.find(['}', '{', '/']) {
            Some(close) if after.as_bytes()[close] == b'}' => {
                on(Token::Text(&rest[..open]));
                on(Token::Field(&after[..close]));
                rest = &after[close + 1..];
            },
            _ => {
                on(Token::Text(&rest[..=open]));
                rest = after;
            },
        }
    }
    on(Token::Text(rest));
}

impl Template {
    /// Parses `path` as an output template: `None` when it holds no
    /// placeholder. A `{barcode}` without `{target}` or `{group}`, an unknown
    /// placeholder, a path whose extension names no output format, and a
    /// template that is not UTF-8 are errors.
    pub fn parse(path: &Path) -> anyhow::Result<Option<Template>> {
        let Some(raw) = path.to_str() else {
            if Template::has_placeholder(path) {
                anyhow::bail!(
                    "-o {}: an output template must be valid UTF-8",
                    path.display()
                );
            }
            return Ok(None);
        };
        let (mut target, mut group, mut barcode) = (false, false, false);
        let mut unknown = None;
        scan(raw, |token| {
            if let Token::Field(name) = token {
                match name {
                    "target" => target = true,
                    "group" => group = true,
                    "barcode" => barcode = true,
                    other => {
                        unknown.get_or_insert(other);
                    },
                }
            }
        });
        if let Some(name) = unknown {
            anyhow::bail!(
                "-o {raw}: unknown placeholder {{{name}}}; the output placeholders are \
                 {{target}}, {{group}} and {{barcode}}"
            );
        }
        if !(target || group || barcode) {
            return Ok(None);
        }
        if !(target || group) {
            anyhow::bail!(
                "-o {raw}: {{barcode}} needs {{target}} or {{group}} in the path to name the \
                 split key"
            );
        }
        if crate::io::from_extension(path).is_none() {
            anyhow::bail!(
                "-o {raw}: an output template needs a known extension (.fastq, .fq, .fastq.gz, \
                 .fq.gz, .fastq.bgz or .bam) to fix the output format"
            );
        }
        let key = if target {
            KeyLevel::Target
        } else {
            KeyLevel::Group
        };
        Ok(Some(Template {
            raw: raw.to_string(),
            key,
            barcode,
        }))
    }

    /// Whether `path` holds one of the output placeholders.
    pub fn has_placeholder(path: &Path) -> bool {
        let raw = path.to_string_lossy();
        let mut found = false;
        scan(&raw, |token| {
            found |= matches!(token, Token::Field("target" | "group" | "barcode"));
        });
        found
    }

    /// The template as given.
    pub fn raw(&self) -> &str {
        &self.raw
    }

    /// The key granularity the template names.
    pub fn key(&self) -> KeyLevel {
        self.key
    }

    /// Whether the template holds `{barcode}`.
    pub fn barcode(&self) -> bool {
        self.barcode
    }

    /// Whether the template holds the placeholder `{name}`.
    fn holds(&self, name: &str) -> bool {
        let mut found = false;
        scan(&self.raw, |token| {
            found |= matches!(token, Token::Field(field) if field == name);
        });
        found
    }

    /// Rejects every sheet name the template would put into a path that is
    /// not one path component (`is_path_component`), and two distinct names
    /// that differ only in letter case, which name one file on a
    /// case-insensitive file system: target names when it holds `{target}`,
    /// group names when it holds `{group}`.
    pub fn check_names(&self, sheet: &Sheet) -> anyhow::Result<()> {
        for (field, on) in [
            ("target", self.holds("target")),
            ("group", self.holds("group")),
        ] {
            if !on {
                continue;
            }
            let mut seen: Vec<&str> = Vec::with_capacity(sheet.targets.len());
            for target in &sheet.targets {
                let name = if field == "target" {
                    &target.name
                } else {
                    &target.group
                };
                if !is_path_component(name) {
                    anyhow::bail!(
                        "-o {}: {field} {name:?} cannot fill {{{field}}} in an output path; \
                         {COMPONENT_RULE}",
                        self.raw
                    );
                }
                if let Some(other) = seen
                    .iter()
                    .find(|other| **other != name && other.eq_ignore_ascii_case(name))
                {
                    anyhow::bail!(
                        "-o {}: {field}s {other:?} and {name:?} differ only in letter case and \
                         would share one output file on a case-insensitive file system; rename \
                         one of them",
                        self.raw
                    );
                }
                seen.push(name);
            }
        }
        Ok(())
    }

    /// Interns the path of every bin of `keys` (each key, then `unassigned`
    /// and `ambiguous`) into `table`, and returns their ids in that order.
    /// For a template without `{barcode}`, whose paths depend on the bin
    /// alone. Two bins expanding to one path are an error.
    pub fn intern_bins(
        &self,
        sheet: &Sheet,
        keys: &Keys,
        table: &KeyTable,
    ) -> anyhow::Result<Vec<KeyId>> {
        bins(sheet, keys)
            .into_iter()
            .map(|(key, group)| {
                let owner = Owner {
                    key,
                    group,
                    barcode: None,
                };
                table.intern(&self.expand_in_group(key, group, None), owner)
            })
            .collect()
    }

    /// The path of every bin of `sheet` at the template's key level, as
    /// `intern_bins` interns them. For a template without `{barcode}`.
    pub fn bin_paths(&self, sheet: &Sheet) -> anyhow::Result<Vec<PathBuf>> {
        let table = KeyTable::default();
        let ids = self.intern_bins(sheet, &Keys::new(sheet, self.key), &table)?;
        Ok(ids.into_iter().map(|id| table.path(id)).collect())
    }

    /// Expands the template for `key` and the barcode call `barcode`
    /// (`unclassified` when `None`). Both `{target}` and `{group}` take
    /// `key`.
    pub fn expand(&self, key: &str, barcode: Option<&str>) -> PathBuf {
        self.expand_in_group(key, key, barcode)
    }

    /// Expands the template as `expand` does, except that a `{group}` beside
    /// a `{target}` key takes `group`, the key's group.
    pub fn expand_in_group(&self, key: &str, group: &str, barcode: Option<&str>) -> PathBuf {
        let group = match self.key {
            KeyLevel::Target => group,
            KeyLevel::Group => key,
        };
        let barcode = barcode.unwrap_or(UNCLASSIFIED);
        let mut out = String::with_capacity(self.raw.len() + 2 * key.len() + barcode.len());
        scan(&self.raw, |token| {
            out.push_str(match token {
                Token::Text(text) => text,
                Token::Field("target") => key,
                Token::Field("group") => group,
                Token::Field(_) => barcode,
            })
        });
        PathBuf::from(out)
    }

    /// The template with every placeholder as `*`: the path the up-front
    /// output guards check.
    pub fn glob(&self) -> PathBuf {
        let mut out = String::with_capacity(self.raw.len());
        scan(&self.raw, |token| match token {
            Token::Text(text) => out.push_str(text),
            Token::Field(_) => out.push('*'),
        });
        PathBuf::from(out)
    }
}

/// The rule `is_path_component` applies, as error messages state it.
const COMPONENT_RULE: &str = "a name in an output path must be one path component (not empty, . or .., and without /, \\ \
     or NUL)";

/// Whether `name` can stand in a path as exactly one component: not empty,
/// `.` or `..`, and without a path separator or NUL. Applied to every sheet
/// name and barcode call an output template substitutes.
pub(crate) fn is_path_component(name: &str) -> bool {
    !matches!(name, "" | "." | "..") && !name.contains(['/', '\\', '\0'])
}

/// Reads a `BC:Z` barcode call as the `{barcode}` value: UTF-8 text that
/// is one path component (`is_path_component`).
pub(crate) fn barcode_name(value: &[u8]) -> anyhow::Result<&str> {
    std::str::from_utf8(value)
        .ok()
        .filter(|v| is_path_component(v))
        .ok_or_else(|| {
            anyhow::anyhow!(
                "the BC:Z barcode {:?} cannot fill {{barcode}} in an output path; \
                 {COMPONENT_RULE}",
                String::from_utf8_lossy(value)
            )
        })
}

/// The bins of `keys`, in bin order, each as its name and group: every key
/// (its group is its first target's), then `unassigned` and `ambiguous`,
/// which are their own groups.
pub(crate) fn bins<'a>(sheet: &'a Sheet, keys: &'a Keys) -> Vec<(&'a str, &'a str)> {
    let mut out: Vec<(&str, &str)> = keys
        .names
        .iter()
        .enumerate()
        .map(|(key, name)| {
            let target = keys
                .of_target
                .iter()
                .position(|&k| k == key)
                .expect("Every key has a target");
            (name.as_str(), sheet.targets[target].group.as_str())
        })
        .collect();
    out.push(("unassigned", "unassigned"));
    out.push(("ambiguous", "ambiguous"));
    out
}

/// What an output path was expanded for: the bin name, its group, and the
/// `{barcode}` value (`None` for a template without `{barcode}`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Owner<'a> {
    /// The bin name: a key, `unassigned` or `ambiguous`.
    pub key: &'a str,
    /// The bin's group.
    pub group: &'a str,
    /// The `{barcode}` value.
    pub barcode: Option<&'a str>,
}

impl std::fmt::Display for Owner<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.key)?;
        if self.group != self.key {
            write!(f, " (group {})", self.group)?;
        }
        if let Some(barcode) = self.barcode {
            write!(f, " with barcode {barcode}")?;
        }
        Ok(())
    }
}

/// An `Owner` held by the table.
#[derive(Debug)]
struct OwnedOwner {
    key: String,
    group: String,
    barcode: Option<String>,
}

impl OwnedOwner {
    fn new(owner: Owner<'_>) -> OwnedOwner {
        OwnedOwner {
            key: owner.key.to_string(),
            group: owner.group.to_string(),
            barcode: owner.barcode.map(str::to_string),
        }
    }

    fn get(&self) -> Owner<'_> {
        Owner {
            key: &self.key,
            group: &self.group,
            barcode: self.barcode.as_deref(),
        }
    }
}

/// The expanded output paths of a run, each interned to a `KeyId` in
/// first-seen order with the owner it was first expanded for. Shared by the
/// render workers: a path already interned is found under the read lock.
#[derive(Debug, Default)]
pub struct KeyTable {
    inner: RwLock<Interned>,
}

/// The interned paths and their ids.
#[derive(Debug, Default)]
struct Interned {
    /// The id of each interned path.
    ids: HashMap<PathBuf, KeyId>,
    /// The path and owner of each id, at index `id`.
    paths: Vec<(PathBuf, OwnedOwner)>,
}

impl Interned {
    /// The id of `path` when interned for `owner`; an error when it was
    /// interned for a different owner.
    fn find(&self, path: &Path, owner: Owner<'_>) -> anyhow::Result<Option<KeyId>> {
        let Some(&id) = self.ids.get(path) else {
            return Ok(None);
        };
        let first = self.paths[id as usize].1.get();
        if first != owner {
            anyhow::bail!(
                "the -o template expands both {first} and {owner} to {}; add placeholders or \
                 separators that tell them apart",
                path.display()
            );
        }
        Ok(Some(id))
    }
}

impl KeyTable {
    /// Returns the id of `path`, interning it for `owner` on first sight. A
    /// path already interned for a different owner is an error naming the
    /// path and both owners.
    pub fn intern(&self, path: &Path, owner: Owner<'_>) -> anyhow::Result<KeyId> {
        if let Some(id) = self.inner.read().unwrap().find(path, owner)? {
            return Ok(id);
        }
        let mut inner = self.inner.write().unwrap();
        if let Some(id) = inner.find(path, owner)? {
            return Ok(id);
        }
        let id = KeyId::try_from(inner.paths.len()).expect("A run has fewer than 2^32 outputs");
        inner
            .paths
            .push((path.to_path_buf(), OwnedOwner::new(owner)));
        inner.ids.insert(path.to_path_buf(), id);
        Ok(id)
    }

    /// Returns the path interned as `id`.
    pub fn path(&self, id: KeyId) -> PathBuf {
        self.inner.read().unwrap().paths[id as usize].0.clone()
    }
}

/// Checks each expanded output path before its file is created: its
/// directories are created, and it must not name a file the run reads or
/// writes otherwise, nor a file an earlier key already opened. Warns once
/// when the opened count passes `OPEN_FILES_WARNING`.
pub(crate) struct OpenGuard {
    /// Files the run reads or writes besides the split outputs, each with
    /// the description an error names it by.
    protected: Vec<(String, PathBuf)>,
    /// Whether the input is stdin, whose file is compared by descriptor.
    stdin: bool,
    /// Every admitted path, resolved (canonical parent and file name), with
    /// the path as expanded.
    opened: Vec<(PathBuf, PathBuf)>,
}

impl OpenGuard {
    /// A guard refusing every path in `protected`, and the file read on
    /// stdin when `stdin`.
    pub(crate) fn new(protected: Vec<(String, PathBuf)>, stdin: bool) -> OpenGuard {
        OpenGuard {
            protected,
            stdin,
            opened: Vec::new(),
        }
    }

    /// Admits `path` as the next output file, creating its directories
    /// before `check`ing it.
    pub(crate) fn admit(&mut self, path: &Path) -> anyhow::Result<()> {
        if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating the output directory {}", parent.display()))?;
        }
        self.check(path)?;
        if self.opened.len() == OPEN_FILES_WARNING + 1 {
            tracing::warn!(
                "The -o template has opened more than {OPEN_FILES_WARNING} output files; raise \
                 the open-file limit (ulimit -n) if opening more fails"
            );
        }
        Ok(())
    }

    /// Checks `path` against the protected files and every path checked
    /// before, and records it. Paths are compared as the file system
    /// resolves them (`guards::same_path`, `guards::resolve`), which reads
    /// file metadata and canonicalizes the path or its parent directory; a
    /// path that does not resolve is compared as written. A directory at
    /// `path` is an error.
    pub(crate) fn check(&mut self, path: &Path) -> anyhow::Result<()> {
        if path.is_dir() {
            anyhow::bail!(
                "the -o template expands to {}, which is a directory",
                path.display()
            );
        }
        for (what, other) in &self.protected {
            if crate::guards::same_path(other, path) {
                anyhow::bail!(
                    "the -o template expands to {}, which is {what}; whittle would overwrite it, \
                     so change the template or the target names",
                    path.display()
                );
            }
        }
        if self.stdin && crate::guards::stdin_is(path) {
            anyhow::bail!(
                "the -o template expands to {}, which is the file being read on stdin; whittle \
                 would overwrite it, so change the template or the target names",
                path.display()
            );
        }
        let resolved = crate::guards::resolve(path).unwrap_or_else(|| path.to_path_buf());
        if let Some((_, earlier)) = self.opened.iter().find(|(r, _)| *r == resolved) {
            anyhow::bail!(
                "the -o template expands two keys to the same file ({} and {}); make each \
                 key's path distinct",
                earlier.display(),
                path.display()
            );
        }
        self.opened.push((resolved, path.to_path_buf()));
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(raw: &str) -> anyhow::Result<Option<Template>> {
        Template::parse(Path::new(raw))
    }

    #[test]
    fn parse_detects_placeholders() {
        let template = parse("out/{target}.bam").unwrap().unwrap();
        assert_eq!(template.key(), KeyLevel::Target);
        assert!(!template.barcode);
        let template = parse("out/{barcode}/{group}.fq").unwrap().unwrap();
        assert_eq!(template.key(), KeyLevel::Group);
        assert!(template.barcode);
        assert_eq!(
            parse("{group}.{target}.fq").unwrap().unwrap().key(),
            KeyLevel::Target
        );
        assert_eq!(parse("out.bam").unwrap(), None);
    }

    #[test]
    fn barcode_alone_is_error() {
        let err = parse("out/{barcode}.bam").unwrap_err().to_string();
        assert!(err.contains("{target}"), "{err}");
        assert!(err.contains("out/{barcode}.bam"), "{err}");
    }

    #[test]
    fn unknown_placeholder_is_error() {
        let err = parse("{sample}.fq").unwrap_err().to_string();
        assert!(err.contains("{sample}"), "{err}");
    }

    #[test]
    fn unknown_extension_is_error() {
        let err = parse("out/{target}.txt").unwrap_err().to_string();
        assert!(err.contains("extension"), "{err}");
    }

    #[test]
    fn expand_substitutes_and_defaults_barcode() {
        let template = parse("{barcode}.{target}.fq.gz").unwrap().unwrap();
        assert_eq!(
            template.expand("16S", None),
            PathBuf::from("unclassified.16S.fq.gz")
        );
        assert_eq!(
            template.expand("unassigned", Some("barcode03")),
            PathBuf::from("barcode03.unassigned.fq.gz")
        );
        let template = parse("out/{group}/{group}.bam").unwrap().unwrap();
        assert_eq!(
            template.expand("ITS", None),
            PathBuf::from("out/ITS/ITS.bam")
        );
    }

    #[test]
    fn expand_fills_group_beside_target() {
        let template = parse("{group}/{target}.fq").unwrap().unwrap();
        assert_eq!(
            template.expand_in_group("V34", "16S", None),
            PathBuf::from("16S/V34.fq")
        );
    }

    #[test]
    fn unmatched_braces_are_text() {
        let template = parse("out{/{target}.fq").unwrap().unwrap();
        assert_eq!(template.expand("A", None), PathBuf::from("out{/A.fq"));
    }

    #[test]
    fn glob_marks_placeholders() {
        let template = parse("out/{barcode}.{target}.bam").unwrap().unwrap();
        assert_eq!(template.glob(), PathBuf::from("out/*.*.bam"));
    }

    /// An owner whose group is its key, without a barcode.
    fn owner(key: &str) -> Owner<'_> {
        Owner {
            key,
            group: key,
            barcode: None,
        }
    }

    #[test]
    fn intern_is_stable() {
        let table = KeyTable::default();
        let a = table.intern(Path::new("out/16S.fq"), owner("16S")).unwrap();
        let b = table.intern(Path::new("out/ITS.fq"), owner("ITS")).unwrap();
        assert_ne!(a, b);
        assert_eq!(
            table.intern(Path::new("out/16S.fq"), owner("16S")).unwrap(),
            a
        );
        assert_eq!(table.path(a), PathBuf::from("out/16S.fq"));
        assert_eq!(table.path(b), PathBuf::from("out/ITS.fq"));
    }

    #[test]
    fn intern_refuses_a_second_owner() {
        let table = KeyTable::default();
        let first = Owner {
            key: "A",
            group: "A",
            barcode: Some("b1"),
        };
        let second = Owner {
            key: "Ab",
            group: "Ab",
            barcode: Some("1"),
        };
        table.intern(Path::new("Ab1.bam"), first).unwrap();
        let err = table
            .intern(Path::new("Ab1.bam"), second)
            .unwrap_err()
            .to_string();
        assert!(err.contains("Ab1.bam"), "{err}");
        assert!(err.contains("A with barcode b1"), "{err}");
        assert!(err.contains("Ab with barcode 1"), "{err}");
    }

    /// A two-target sheet with the given names and groups.
    fn sheet(rows: [(&str, &str); 2]) -> Sheet {
        let [(t1, g1), (t2, g2)] = rows;
        Sheet::parse_tsv(&format!(
            "target\tfwd\trev\tgroup\n\
             {t1}\tACGTACGTACGTACGTAC\tTTGCATTGCATTGCATTG\t{g1}\n\
             {t2}\tGGCCAAGGCCAAGGCCAA\tCATGCATGCATGCATGCA\t{g2}\n"
        ))
        .unwrap()
    }

    #[test]
    fn check_names_refuses_names_that_are_not_one_component() {
        let target = parse("out/{target}.fq").unwrap().unwrap();
        let group = parse("out/{group}.fq").unwrap().unwrap();
        for bad in ["../x", "a/b"] {
            let err = target
                .check_names(&sheet([(bad, "g"), ("ok", "g")]))
                .unwrap_err()
                .to_string();
            assert!(err.contains(&format!("target {bad:?}")), "{err}");
            // A group template puts group names in paths, not target names.
            group
                .check_names(&sheet([(bad, "g"), ("ok", "g")]))
                .unwrap();
            let err = group
                .check_names(&sheet([("a", bad), ("ok", "g")]))
                .unwrap_err()
                .to_string();
            assert!(err.contains(&format!("group {bad:?}")), "{err}");
        }
        target
            .check_names(&sheet([("a", "a/b"), ("ok", "g")]))
            .unwrap();
    }

    #[test]
    fn check_names_refuses_names_that_differ_only_in_case() {
        let target = parse("out/{target}.fq").unwrap().unwrap();
        let group = parse("out/{group}.fq").unwrap().unwrap();
        let err = target
            .check_names(&sheet([("its", "g"), ("ITS", "g")]))
            .unwrap_err()
            .to_string();
        assert!(err.contains("\"its\" and \"ITS\""), "{err}");
        assert!(err.contains("letter case"), "{err}");
        // A group template puts group names in paths, not target names.
        group
            .check_names(&sheet([("its", "g"), ("ITS", "g")]))
            .unwrap();
        let err = group
            .check_names(&sheet([("a", "g16"), ("b", "G16")]))
            .unwrap_err()
            .to_string();
        assert!(err.contains("groups \"g16\" and \"G16\""), "{err}");
        target
            .check_names(&sheet([("a", "g16"), ("b", "G16")]))
            .unwrap();
    }

    #[test]
    fn bin_paths_cover_every_bin() {
        let template = parse("out/{group}/{target}.fq").unwrap().unwrap();
        let paths = template
            .bin_paths(&sheet([("A", "g1"), ("B", "g1")]))
            .unwrap();
        assert_eq!(
            paths,
            [
                "out/g1/A.fq",
                "out/g1/B.fq",
                "out/unassigned/unassigned.fq",
                "out/ambiguous/ambiguous.fq"
            ]
            .map(PathBuf::from)
        );
        let merged = parse("{group}.fq").unwrap().unwrap();
        assert_eq!(
            merged.bin_paths(&sheet([("A", "g"), ("B", "g")])).unwrap(),
            ["g.fq", "unassigned.fq", "ambiguous.fq"].map(PathBuf::from)
        );
    }

    #[test]
    fn bin_paths_refuse_two_bins_on_one_path() {
        let template = parse("{target}{group}.fq").unwrap().unwrap();
        let err = template
            .bin_paths(&sheet([("a", "bc"), ("ab", "c")]))
            .unwrap_err()
            .to_string();
        assert!(err.contains("abc.fq"), "{err}");
        assert!(err.contains("a (group bc)"), "{err}");
        assert!(err.contains("ab (group c)"), "{err}");
    }

    #[test]
    fn intern_is_shared_across_threads() {
        let table = KeyTable::default();
        let ids: Vec<Vec<KeyId>> = std::thread::scope(|s| {
            let handles: Vec<_> = (0..4)
                .map(|_| {
                    s.spawn(|| {
                        (0..50)
                            .map(|i| {
                                let key = format!("k{}", i % 7);
                                table.intern(Path::new(&key), owner(&key)).unwrap()
                            })
                            .collect()
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        for run in &ids {
            assert_eq!(run, &ids[0]);
        }
        assert_eq!(table.path(ids[0][3]), PathBuf::from("k3"));
    }

    #[test]
    fn open_guard_refuses_protected_and_aliased_paths() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.fq");
        std::fs::write(&input, "").unwrap();
        let mut guard = OpenGuard::new(vec![("the input".into(), input.clone())], false);
        let err = guard.admit(&input).unwrap_err().to_string();
        assert!(err.contains("the input"), "{err}");
        assert!(err.contains("in.fq"), "{err}");

        let first = dir.path().join("sub/a/../16S.fq");
        guard.admit(&first).unwrap();
        assert!(
            dir.path().join("sub/a").is_dir(),
            "Parent directories are created"
        );
        let alias = dir.path().join("sub/16S.fq");
        let err = guard.admit(&alias).unwrap_err().to_string();
        assert!(err.contains("16S.fq"), "{err}");
        guard.admit(&dir.path().join("sub/ITS.fq")).unwrap();

        // `check` compares without creating directories.
        let mut guard = OpenGuard::new(vec![("the input".into(), input.clone())], false);
        guard.check(&dir.path().join("new/16S.fq")).unwrap();
        assert!(!dir.path().join("new").exists());
        assert!(guard.check(&input).is_err());
    }

    #[test]
    fn barcode_values_that_are_not_file_names_are_refused() {
        for bad in ["", ".", "..", "a/b", "a\\b", "a\0b"] {
            assert!(barcode_name(bad.as_bytes()).is_err(), "{bad:?}");
            assert!(!is_path_component(bad), "{bad:?}");
        }
        assert!(barcode_name(&[0xff]).is_err());
        assert_eq!(barcode_name(b"SQK_barcode01").unwrap(), "SQK_barcode01");
    }
}
