//! The pedigree rules of step 38 and the Nextflow FAMILY workflow.
//!
//! A families file is a pedigree (sample, sex, father, mother, partner, tab
//! separated) with an optional leading family column, so one file can
//! declare several families. This crate reads it as `bin/family_report.py`
//! did, with the same checks and the same messages, word for word:
//!
//! - [`read_pedigree`], [`split_families`], [`family_rows`]: the rows of the
//!   file and of one family;
//! - [`build_members`]: every person the rows name, and what is wrong;
//! - [`read_members`]: the `members.tsv` that `og family prepare` writes;
//! - [`Kinship`], [`declared_relationship`], [`degree`]: the relatedness a
//!   pedigree implies, binned as KING bins it;
//! - [`Par`]: the pseudoautosomal regions of a BED file.
//!
//! A problem is a [`Die`]: `ERROR: ` lines on stderr and an exit code.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufReader, Write};

use indexmap::IndexMap;
use og_io::{Lines, Utf8};
use pycompat::text::{lower, strip};

/// A families file's columns; a pedigree is the same without `family`.
pub const COLUMNS: [&str; 6] = ["family", "sample", "sex", "father", "mother", "partner"];
/// The pedigree's own columns.
pub const PED_COLUMNS: [&str; 5] = ["sample", "sex", "father", "mother", "partner"];
/// The family of a file without the family column.
pub const DEFAULT_FAMILY: &str = "family";
/// Relatedness (twice the kinship coefficient) bins, as KING draws them:
/// the geometric middle between two degrees.
pub const BINS: [(f64, &str); 4] = [
    (0.707, "same person"),
    (0.354, "first degree"),
    (0.177, "second degree"),
    (0.0884, "third degree"),
];

/// A stop: each line is printed as `ERROR: <line>` on stderr, then the
/// tool exits with `code`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Die {
    pub lines: Vec<String>,
    pub code: u8,
}

impl Die {
    /// A stop with exit 2, the code for a problem in the input.
    pub fn new(lines: Vec<String>) -> Die {
        Die { lines, code: 2 }
    }

    /// One line, exit 2.
    pub fn one(line: impl Into<String>) -> Die {
        Die::new(vec![line.into()])
    }

    /// A file that cannot be read or written.
    pub fn io(path: &str, e: &io::Error) -> Die {
        Die::one(format!("cannot read {path}: {e}"))
    }

    /// Prints the lines on stderr and returns the exit code.
    pub fn report(&self) -> u8 {
        let mut err = io::stderr().lock();
        for line in &self.lines {
            let _ = writeln!(err, "ERROR: {line}");
        }
        self.code
    }
}

/// Values that mean "not given".
fn is_missing_word(v: &str) -> bool {
    matches!(v, "" | "." | "0" | "-9")
}

/// `missing(v)`: `None`, or a value that strips to a missing word.
pub fn missing(v: Option<&str>) -> bool {
    v.is_none_or(|v| is_missing_word(strip(v, None)))
}

/// The sex words a pedigree may use: `Some(Some(sex))`, `Some(None)` for
/// unknown, `None` when the word is not one.
pub fn sex_word(v: &str) -> Option<Option<&'static str>> {
    match v {
        "male" | "m" | "1" => Some(Some("male")),
        "female" | "f" | "2" => Some(Some("female")),
        "unknown" | "" | "." | "0" | "-9" => Some(None),
        _ => None,
    }
}

/// A sample or family name: letters, digits, `.`, `_` and `-`, not
/// starting with `-`.
pub fn bad_name(v: &str) -> bool {
    v.is_empty()
        || v.starts_with('-')
        || !v
            .bytes()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'.' | b'_' | b'-'))
}

/// Text lines of a file as Python's `open(path)` reads them: universal
/// newlines, strict UTF-8, each line without its ending.
pub fn text_lines(path: &str) -> Result<Vec<String>, Die> {
    let f = File::open(path).map_err(|e| Die::io(path, &e))?;
    Lines::new(BufReader::new(f), Utf8::Strict)
        .collect::<io::Result<Vec<String>>>()
        .map_err(|e| Die::io(path, &e))
}

/// One row of a pedigree: the header's columns and the line number.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub line: usize,
    /// Column name to stripped value, in header order (a column named
    /// twice keeps its last value).
    pub vals: IndexMap<String, String>,
}

impl Row {
    /// `r.get(k)`: the value, `None` when the header has no such column.
    pub fn get(&self, k: &str) -> Option<&str> {
        self.vals.get(k).map(String::as_str)
    }

    /// `r.get(k, "").strip()`.
    pub fn field(&self, k: &str) -> &str {
        strip(self.get(k).unwrap_or(""), None)
    }
}

/// The rows of a pedigree file, and the problems found. The first
/// non-blank line is the header (leading `#` dropped); later lines
/// starting with `#` are comments.
pub fn read_pedigree(path: &str) -> Result<(Vec<Row>, Vec<String>), Die> {
    Ok(parse_pedigree(path, &text_lines(path)?))
}

/// [`read_pedigree`] on lines already read; `path` names the file in messages.
pub fn parse_pedigree(path: &str, lines: &[String]) -> (Vec<Row>, Vec<String>) {
    let (mut rows, mut errors) = (Vec::new(), Vec::new());
    let mut header: Option<Vec<String>> = None;
    for (i, line) in lines.iter().enumerate() {
        let n = i + 1;
        let line = line.trim_end_matches(['\r', '\n']);
        if strip(line, None).is_empty() {
            continue;
        }
        let Some(head) = &header else {
            let head: Vec<String> = line
                .trim_start_matches('#')
                .split('\t')
                .map(|c| lower(strip(c, None)))
                .collect();
            if !head.iter().any(|c| c == "sample") {
                return (
                    Vec::new(),
                    vec![format!(
                        "{path}: the first line must be a header with a 'sample' column \
                         (columns: {}, tab separated); it is '{line}'",
                        COLUMNS.join(", ")
                    )],
                );
            }
            let bad: Vec<&str> = head
                .iter()
                .map(String::as_str)
                .filter(|c| !COLUMNS.contains(c))
                .collect();
            if !bad.is_empty() {
                errors.push(format!(
                    "{path}: unknown column{} {}; the columns are {}",
                    if bad.len() > 1 { "s" } else { "" },
                    bad.join(", "),
                    COLUMNS.join(", ")
                ));
            }
            header = Some(head);
            continue;
        };
        if line.starts_with('#') {
            continue;
        }
        let vals: Vec<&str> = line.split('\t').map(|v| strip(v, None)).collect();
        if vals.len() > head.len() {
            errors.push(format!(
                "{path} line {n}: {} fields, the header has {} (fields are tab separated)",
                vals.len(),
                head.len()
            ));
            continue;
        }
        let mut rec = IndexMap::new();
        for (k, col) in head.iter().enumerate() {
            rec.insert(col.clone(), vals.get(k).copied().unwrap_or("").to_owned());
        }
        rows.push(Row { line: n, vals: rec });
    }
    if header.is_none() {
        errors.push(format!("{path} is empty"));
    }
    (rows, errors)
}

/// One person of a family.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Member {
    pub sex: Option<String>,
    pub father: Option<String>,
    pub mother: Option<String>,
    pub partner: Option<String>,
    /// Has a row of their own in the pedigree.
    pub row: bool,
    /// Has a gVCF in the joint call.
    pub sequenced: bool,
}

impl Member {
    /// The father, then the mother, those that are named.
    pub fn parents(&self) -> impl Iterator<Item = &str> {
        [&self.father, &self.mother]
            .into_iter()
            .filter_map(|p| p.as_deref())
    }

    fn rel(&self, k: &str) -> Option<&str> {
        match k {
            "father" => self.father.as_deref(),
            "mother" => self.mother.as_deref(),
            _ => self.partner.as_deref(),
        }
    }
}

/// Every member by name, in the order the pedigree names them.
pub type Members = IndexMap<String, Member>;

/// Every person the rows name: the rows, and the parents and partners they
/// name without a row of their own (not sequenced, no parents known), with
/// the problems found. `declared` is the sex given elsewhere (`source`
/// names where, in messages).
pub fn build_members(
    rows: &[Row],
    declared: &HashMap<String, String>,
    path: &str,
    source: &str,
) -> (Members, Vec<String>) {
    let mut errors = Vec::new();
    let mut members = Members::new();
    for r in rows {
        let where_ = format!("{path} line {}", r.line);
        let s = r.field("sample");
        if missing(Some(s)) {
            errors.push(format!("{where_}: no sample name"));
            continue;
        }
        if bad_name(s) {
            errors.push(format!(
                "{where_}: sample name '{s}' may hold only letters, digits, '.', '_' and '-' \
                 and may not start with '-'"
            ));
            continue;
        }
        if members.contains_key(s) {
            errors.push(format!(
                "{where_}: sample '{s}' has a second row; each person has one"
            ));
            continue;
        }
        let sex_raw = lower(r.field("sex"));
        let sex = match sex_word(&sex_raw) {
            Some(x) => x,
            None => {
                errors.push(format!(
                    "{where_}: sex '{}' of '{s}' is not male, female or unknown",
                    r.get("sex").unwrap_or("None")
                ));
                None
            }
        };
        let mut m = Member {
            sex: sex.map(str::to_owned),
            row: true,
            ..Member::default()
        };
        for k in ["father", "mother", "partner"] {
            let v = r.field(k);
            let val = if missing(Some(v)) {
                None
            } else if bad_name(v) {
                errors.push(format!(
                    "{where_}: {k} '{v}' of '{s}' is not a valid sample name"
                ));
                None
            } else if v == s {
                errors.push(format!("{where_}: '{s}' is named as its own {k}"));
                None
            } else {
                Some(v.to_owned())
            };
            match k {
                "father" => m.father = val,
                "mother" => m.mother = val,
                _ => m.partner = val,
            }
        }
        if m.father.is_some() && m.father == m.mother {
            errors.push(format!(
                "{where_}: '{s}' has the same person, '{}', as father and mother",
                m.father.as_deref().unwrap_or_default()
            ));
        }
        if let Some(d) = declared.get(s).filter(|d| !d.is_empty()) {
            if let Some(own) = &m.sex
                && own != d
            {
                errors.push(format!(
                    "{where_}: the pedigree says '{s}' is {own}, {source} says {d}"
                ));
            }
            if m.sex.is_none() {
                m.sex = Some(d.clone());
            }
        }
        members.insert(s.to_owned(), m);
    }

    // People named but without a row: parents and partners not in the file.
    let named: Vec<Member> = members.values().cloned().collect();
    for m in &named {
        for (k, sex) in [
            ("father", Some("male")),
            ("mother", Some("female")),
            ("partner", None),
        ] {
            if let Some(p) = m.rel(k)
                && !members.contains_key(p)
            {
                members.insert(
                    p.to_owned(),
                    Member {
                        sex: sex.map(str::to_owned),
                        ..Member::default()
                    },
                );
            }
        }
    }

    let sex_of = |p: &str| members.get(p).and_then(|m| m.sex.as_deref());
    for (s, m) in &members {
        if let Some(f) = &m.father
            && sex_of(f) == Some("female")
        {
            errors.push(format!(
                "'{f}' is the father of '{s}' but its sex is female"
            ));
        }
        if let Some(mo) = &m.mother
            && sex_of(mo) == Some("male")
        {
            errors.push(format!("'{mo}' is the mother of '{s}' but its sex is male"));
        }
    }
    let fathers: HashSet<&str> = members
        .values()
        .filter_map(|m| m.father.as_deref())
        .collect();
    let mothers: HashSet<&str> = members
        .values()
        .filter_map(|m| m.mother.as_deref())
        .collect();
    let mut both: Vec<&str> = fathers.intersection(&mothers).copied().collect();
    both.sort_unstable();
    for p in both {
        errors.push(format!("'{p}' is named as a father and as a mother"));
    }

    for (s, m) in &members {
        if let Some(p) = &m.partner
            && let Some(pm) = members.get(p)
            && pm.row
            && pm.partner.as_deref().is_some_and(|q| q != s)
        {
            errors.push(format!(
                "'{s}' names '{p}' as partner, but '{p}' names '{}'",
                pm.partner.as_deref().unwrap_or_default()
            ));
        }
    }

    // No one is their own ancestor.
    for s in members.keys() {
        let mut seen: HashSet<&str> = HashSet::new();
        let mut todo: Vec<&str> = vec![s];
        'walk: while let Some(x) = todo.pop() {
            let Some(mx) = members.get(x) else { continue };
            for p in mx.parents() {
                if p == s {
                    errors.push(format!(
                        "'{s}' is among its own ancestors: the parent columns form a loop"
                    ));
                    break 'walk;
                }
                if seen.insert(p) {
                    todo.push(p);
                }
            }
        }
    }
    (members, errors)
}

/// The families of a file: one family keyed `None` for a file without
/// the family column, else each family's rows in file order.
pub type Families = IndexMap<Option<String>, Vec<Row>>;

/// The families of a file's rows, and the problems found. The family column
/// is filled on every row or on none; a person is in one family; the
/// parents and partner a row names are in its own family; a family name
/// follows the sample-name rule and is no person's name.
pub fn split_families(rows: &[Row], path: &str) -> (Families, Vec<String>) {
    let single = |rows: &[Row]| -> Families {
        let mut f = Families::new();
        if !rows.is_empty() {
            f.insert(None, rows.to_vec());
        }
        f
    };
    if !rows.iter().any(|r| r.vals.contains_key("family")) {
        return (single(rows), Vec::new());
    }
    let mut errors = Vec::new();
    let mut fams: IndexMap<String, Vec<Row>> = IndexMap::new();
    let mut where_of: HashMap<String, (String, usize)> = HashMap::new();
    let filled: Vec<&Row> = rows.iter().filter(|r| !missing(r.get("family"))).collect();
    for r in rows {
        let fam = r.field("family");
        let where_ = format!("{path} line {}", r.line);
        if missing(Some(fam)) {
            if let Some(first) = filled.first() {
                errors.push(format!(
                    "{where_}: no family name, but other rows have one (line {}: '{}'); fill the family \
                     column on every row, or leave it out of the file",
                    first.line,
                    first.field("family")
                ));
            }
            continue;
        }
        if bad_name(fam) {
            errors.push(format!(
                "{where_}: family name '{fam}' may hold only letters, digits, '.', '_' and '-' and may not start \
                 with '-' (it names the family's folder)"
            ));
            continue;
        }
        let s = r.field("sample");
        if !s.is_empty() && !missing(Some(s)) {
            if let Some((f0, l0)) = where_of.get(s)
                && f0 != fam
            {
                errors.push(format!(
                    "{where_}: sample '{s}' is in family '{fam}' here and in family '{f0}' on line {l0}; a person \
                     belongs to one family"
                ));
                continue;
            }
            where_of
                .entry(s.to_owned())
                .or_insert_with(|| (fam.to_owned(), r.line));
        }
        fams.entry(fam.to_owned()).or_default().push(r.clone());
    }
    if filled.is_empty() {
        return (single(rows), errors);
    }
    let mut people: HashMap<&str, usize> = HashMap::new();
    for r in rows {
        for k in ["sample", "father", "mother", "partner"] {
            let v = r.field(k);
            if !v.is_empty() && !missing(Some(v)) {
                people.entry(v).or_insert(r.line);
            }
        }
    }
    let mut names: Vec<&String> = fams.keys().collect();
    names.sort();
    for fam in names {
        if let Some(l) = people.get(fam.as_str()) {
            errors.push(format!(
                "{path}: the family name '{fam}' is also a person's name (line {l}); its folder would be the \
                 sample's. Give the family another name"
            ));
        }
    }
    for (fam, rs) in &fams {
        for r in rs {
            let s = r.field("sample");
            for k in ["father", "mother", "partner"] {
                let v = r.field(k);
                if !v.is_empty()
                    && !missing(Some(v))
                    && let Some((f0, l0)) = where_of.get(v)
                    && f0 != fam
                {
                    errors.push(format!(
                        "{path} line {}: the {k} of '{s}', '{v}', is in family '{f0}' (line {l0}), '{s}' in family \
                         '{fam}'; a family lists its own parents and partners",
                        r.line
                    ));
                }
            }
        }
    }
    let out = fams.into_iter().map(|(k, v)| (Some(k), v)).collect();
    (out, errors)
}

/// One family of a families file (or a pedigree): the family named, else
/// the one `sample` is in, else the file's only family. The name is
/// [`DEFAULT_FAMILY`] (or the `family` given) for a file without the
/// family column; `None` with no rows when `sample` is in no family.
pub fn family_rows(
    path: &str,
    family: Option<&str>,
    sample: Option<&str>,
) -> Result<(Option<String>, Vec<Row>), Die> {
    let (rows, errors) = read_pedigree(path)?;
    if !errors.is_empty() {
        return Err(Die::new(errors));
    }
    let (mut fams, errors) = split_families(&rows, path);
    if !errors.is_empty() {
        return Err(Die::new(errors));
    }
    let has_sample = |rs: &[Row], s: &str| rs.iter().any(|r| r.field("sample") == s);
    if let Some(rs) = fams.shift_remove(&None) {
        let name = family.unwrap_or(DEFAULT_FAMILY);
        if bad_name(name) {
            return Err(Die::one(format!(
                "family name '{name}' may hold only letters, digits, '.', '_' and '-' and may not start with '-' \
                 (it names the family's folder)"
            )));
        }
        for r in &rs {
            if ["sample", "father", "mother", "partner"]
                .iter()
                .any(|k| r.field(k) == name)
            {
                return Err(Die::one(format!(
                    "{path} line {}: the family name '{name}' is also a person's name; its folder would be the \
                     sample's. Give the family another name",
                    r.line
                )));
            }
        }
        if let Some(s) = sample
            && family.is_none()
            && !has_sample(&rs, s)
        {
            return Ok((Some(name.to_owned()), Vec::new()));
        }
        return Ok((Some(name.to_owned()), rs));
    }
    let mut names: Vec<String> = fams.keys().flatten().cloned().collect();
    names.sort();
    if let Some(f) = family {
        let key = Some(f.to_owned());
        return match fams.shift_remove(&key) {
            Some(rs) => Ok((key, rs)),
            None => Err(Die::one(format!(
                "{path} has no family '{f}'; its families are: {}",
                if names.is_empty() {
                    "none".to_owned()
                } else {
                    names.join(", ")
                }
            ))),
        };
    }
    if let Some(s) = sample {
        for name in &names {
            let key = Some(name.clone());
            if has_sample(&fams[&key], s) {
                let rs = fams.shift_remove(&key).unwrap_or_default();
                return Ok((key, rs));
            }
        }
        return Ok((None, Vec::new()));
    }
    if names.len() > 1 {
        return Err(Die::one(format!(
            "{path} holds {} families ({}); name the one to run (--family, or --family_name in Nextflow)",
            names.len(),
            names.join(", ")
        )));
    }
    match names.into_iter().next() {
        Some(name) => {
            let key = Some(name);
            let rs = fams.shift_remove(&key).unwrap_or_default();
            Ok((key, rs))
        }
        None => Ok((Some(DEFAULT_FAMILY.to_owned()), Vec::new())),
    }
}

/// A family's rows as a pedigree of its own: the pedigree columns, no
/// family column, `.` for an empty value.
pub fn pedigree_text(rows: &[Row]) -> String {
    let mut out = PED_COLUMNS.join("\t");
    out.push('\n');
    for r in rows {
        let cells: Vec<&str> = PED_COLUMNS
            .iter()
            .map(|k| match r.field(k) {
                "" => ".",
                v => v,
            })
            .collect();
        out.push_str(&cells.join("\t"));
        out.push('\n');
    }
    out
}

/// A `\t`-separated table read as Python's `dict(zip(head, line.split("\t")))`
/// per line: a short line lacks the last columns.
pub fn read_table(path: &str, strip_hash: bool) -> Result<Vec<IndexMap<String, String>>, Die> {
    let lines = text_lines(path)?;
    let mut it = lines.into_iter();
    let head_line = it.next().unwrap_or_default();
    let head_line = if strip_hash {
        head_line.trim_start_matches('#').to_owned()
    } else {
        head_line
    };
    let head: Vec<String> = head_line.split('\t').map(str::to_owned).collect();
    Ok(it
        .map(|l| {
            let mut r = IndexMap::new();
            for (k, v) in head.iter().zip(l.split('\t')) {
                r.insert(k.clone(), v.to_owned());
            }
            r
        })
        .collect())
}

/// The value of `k` in a table row, or a stop naming the file.
pub fn cell<'a>(r: &'a IndexMap<String, String>, k: &str, path: &str) -> Result<&'a str, Die> {
    r.get(k)
        .map(String::as_str)
        .ok_or_else(|| Die::one(format!("{path}: a row has no '{k}' column")))
}

/// The `members.tsv` that `og family prepare` writes.
pub fn read_members(path: &str) -> Result<Members, Die> {
    let mut members = Members::new();
    for r in read_table(path, false)? {
        let opt = |k: &str| -> Result<Option<String>, Die> {
            let v = cell(&r, k, path)?;
            Ok((!v.is_empty()).then(|| v.to_owned()))
        };
        let m = Member {
            sex: opt("sex")?,
            father: opt("father")?,
            mother: opt("mother")?,
            partner: opt("partner")?,
            row: cell(&r, "in_pedigree", path)? == "yes",
            sequenced: cell(&r, "sequenced", path)? == "yes",
        };
        members.insert(cell(&r, "sample", path)?.to_owned(), m);
    }
    Ok(members)
}

/// phi(a, b), the kinship coefficient of two members, from the pedigree.
/// Founders (no parent known) are unrelated; an unknown parent is a
/// founder. The same floating-point operations as the Python, memoised.
pub struct Kinship<'a> {
    members: &'a Members,
    depth: HashMap<String, u32>,
    phi: HashMap<(String, String), f64>,
}

impl<'a> Kinship<'a> {
    pub fn new(members: &'a Members) -> Self {
        Kinship {
            members,
            depth: HashMap::new(),
            phi: HashMap::new(),
        }
    }

    fn parents(&self, x: &str) -> (Option<&'a str>, Option<&'a str>) {
        self.members
            .get(x)
            .map_or((None, None), |m| (m.father.as_deref(), m.mother.as_deref()))
    }

    fn depth(&mut self, x: &str) -> u32 {
        if let Some(&d) = self.depth.get(x) {
            return d;
        }
        let (f, m) = self.parents(x);
        let d = [f, m]
            .into_iter()
            .flatten()
            .map(|p| self.depth(p))
            .max()
            .map_or(0, |d| d + 1);
        self.depth.insert(x.to_owned(), d);
        d
    }

    /// The kinship coefficient of `a` and `b`.
    pub fn phi(&mut self, a: &str, b: &str) -> f64 {
        if a == b {
            return match self.parents(a) {
                (Some(f), Some(m)) => 0.5 * (1.0 + self.phi(f, m)),
                _ => 0.5,
            };
        }
        let key = if a < b {
            (a.to_owned(), b.to_owned())
        } else {
            (b.to_owned(), a.to_owned())
        };
        if let Some(&v) = self.phi.get(&key) {
            return v;
        }
        let (a, b) = if self.depth(a) < self.depth(b) {
            (b, a)
        } else {
            (a, b)
        };
        let v = match self.parents(a) {
            (None, None) => 0.0,
            (f, m) => {
                let terms: Vec<f64> = [f, m]
                    .into_iter()
                    .flatten()
                    .map(|p| self.phi(p, b))
                    .collect();
                0.5 * pycompat::num::sum(terms)
            }
        };
        self.phi.insert(key, v);
        v
    }
}

/// The KING degree of a relatedness.
pub fn degree(r: f64) -> &'static str {
    BINS.iter()
        .find(|(lo, _)| r >= *lo)
        .map_or("unrelated", |(_, name)| name)
}

/// The relationship the pedigree declares between `a` and `b`, with
/// `expected` (twice their kinship) for the pairs it does not name.
pub fn declared_relationship(a: &str, b: &str, members: &Members, expected: f64) -> String {
    let none = Member::default();
    let ma = members.get(a).unwrap_or(&none);
    let mb = members.get(b).unwrap_or(&none);
    if mb.parents().any(|p| p == a) || ma.parents().any(|p| p == b) {
        return "parent-child".to_owned();
    }
    let pa: HashSet<&str> = ma.parents().collect();
    let shared = mb
        .parents()
        .collect::<HashSet<&str>>()
        .intersection(&pa)
        .count();
    match shared {
        2 => return "full siblings".to_owned(),
        1 => return "half siblings".to_owned(),
        _ => {}
    }
    if ma.partner.as_deref() == Some(b) || mb.partner.as_deref() == Some(a) {
        return if expected == 0.0 {
            "partners".to_owned()
        } else {
            format!("partners, {}", degree(expected))
        };
    }
    if expected > 0.0 {
        degree(expected).to_owned()
    } else {
        "unrelated".to_owned()
    }
}

/// The pseudoautosomal regions of a BED file, per contig.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Par(HashMap<String, Vec<(i64, i64)>>);

impl Par {
    /// Reads `path` (no path: no regions). A line with fewer than three
    /// tab-separated fields is skipped.
    pub fn read(path: Option<&str>) -> Result<Par, Die> {
        let mut par: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        let Some(path) = path else {
            return Ok(Par(par));
        };
        for (i, line) in text_lines(path)?.iter().enumerate() {
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 3 {
                continue;
            }
            let num = |s: &str| {
                pycompat::num::parse_int(s)
                    .ok()
                    .and_then(|v| i64::try_from(v).ok())
                    .ok_or_else(|| {
                        Die::one(format!(
                            "{path} line {}: '{s}' is not a whole number",
                            i + 1
                        ))
                    })
            };
            let (s, e) = (num(c[1])?, num(c[2])?);
            par.entry(c[0].to_owned()).or_default().push((s, e));
        }
        Ok(Par(par))
    }

    /// No region at all (no file, or an empty one).
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Whether the 1-based `pos` on `chrom` is in a region (BED: 0-based
    /// start, end exclusive).
    pub fn contains(&self, chrom: &str, pos: i64) -> bool {
        self.0
            .get(chrom)
            .is_some_and(|v| v.iter().any(|&(s, e)| s < pos && pos <= e))
    }
}

#[cfg(test)]
mod tests;
