//! `og live family`: the live site's family page, from the files the
//! family steps wrote.
//!
//! A family folder under GENOME_DIR (one holding `members.tsv` or
//! `pedigree.tsv`, as steps 38 and 39 and the Nextflow family workflow
//! leave it) is a run of kind `family`. Its page shows the pedigree (drawn
//! from that file at request time, [`og_pedigree::layout`]), each child's
//! chromosome map (`<child>_segments.svg`), the relationship check, the
//! shared segments and the couple carrier screen, as the steps wrote them.
//!
//! ```text
//! og live family state DIR [--name NAME]
//!     state.json's plugins.family of the family in DIR: its members, each
//!     child with its parents and map, the relationship and shared-segment
//!     rows, each couple's counts, the plan and the files the page may open.
//! og live family page DIR [--name NAME] [--root ROOT] [--input FILE|-]
//!     family.json: the state, the pedigree's layout (so each person is a
//!     link), the relationship and shared-segment rows with their measured
//!     numbers, each couple's gene rows, each child's inherited stretches
//!     (<child>_segments.tsv), the de novo and compound-het candidates
//!     (de_novo.tsv, compound_het.tsv) with whether each child was assessed
//!     and why not, and GRCh38's chromosome lengths. INPUT (JSON) holds what
//!     only the site knows: "runs" (each member run's entry, copied into
//!     family.json as it is), "candidates" (each child's run folder: the
//!     compound-het line then says whether step 23 wrote the candidates step
//!     38 pairs) and "summaries" (each member's run folder: the members side
//!     by side from each summary.json). A run folder is read only inside
//!     ROOT. Exit 1 when a file of the family cannot be read.
//! og live family families ROOT
//!     every family folder under ROOT as /runs.json lists it: name, dir,
//!     kind, members ([] while the pedigree cannot be read) and unreadable
//!     (while family.json cannot be built from the folder).
//! og live family files DIR
//!     the files of the family folder the page may open, as a JSON list.
//!     Exit 2 when the folder holds no pedigree that names anyone.
//! og live family svg DIR [--title TEXT]
//!     the pedigree as an SVG (default title: the folder's name). Exit 2 when
//!     the folder holds no pedigree that names anyone.
//! ```
//!
//! It computes nothing about genetics: it reads the files the steps wrote,
//! counts their rows and draws the pedigree the file declares. A file of
//! the family folder that links out of it is never read, and the page never
//! shows where it points.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{self, BufReader, Read, Write};
use std::path::{Component, Path, PathBuf};
use std::time::SystemTime;

use indexmap::IndexMap;
use og_io::{Lines, Utf8};
use og_pedigree::layout::{self, People};
use pycompat::json::{Dump, Value, dumps, loads};
use pycompat::value::{get, py_str, truthy};

const USAGE: &str = "usage: og live family state DIR [--name NAME]\n       \
     og live family page DIR [--name NAME] [--root ROOT] [--input FILE|-]\n       \
     og live family families ROOT\n       \
     og live family files DIR\n       \
     og live family svg DIR [--title TEXT]\n";

/// The family's result files the page may open, besides the drawn pedigree
/// and the per-child and per-couple files (named from the members).
const FIXED: [&str; 7] = [
    "family_report.md",
    "ibd_report.md",
    "relationships.tsv",
    "ibd_pairs.tsv",
    "family_plan.tsv",
    "members.tsv",
    "pedigree.tsv",
];
pub const TRIO_FILES: [&str; 2] = ["de_novo.tsv", "compound_het.tsv"];
const SEGMENT_KEYS: [&str; 10] = [
    "side",
    "parent",
    "chrom",
    "start",
    "end",
    "length_bp",
    "sites",
    "discordant",
    "ancestor",
    "basis",
];
const DE_NOVO_KEYS: [&str; 12] = [
    "child",
    "chrom",
    "pos",
    "ref",
    "alt",
    "child_gt",
    "father_gt",
    "mother_gt",
    "child_dp",
    "child_gq",
    "child_ab",
    "quality",
];
const COMP_KEYS: [&str; 7] = [
    "child",
    "gene",
    "variant_a",
    "variant_b",
    "a_from",
    "b_from",
    "phase",
];
const REL_KEYS: [&str; 10] = [
    "sample_a",
    "sample_b",
    "declared",
    "expected_relatedness",
    "relatedness",
    "ibs0",
    "ibs2",
    "sites",
    "verdict",
    "reason",
];
const IBD_KEYS: [&str; 9] = [
    "sample_a",
    "sample_b",
    "declared",
    "segments",
    "shared_cm",
    "ibd2_cm",
    "longest_cm",
    "relatedness",
    "degree",
];
const COUPLE_KEYS: [&str; 10] = [
    "gene",
    "inheritance",
    "condition",
    "result",
    "offspring_risk",
    "partner_a",
    "findings_a",
    "partner_b",
    "findings_b",
    "explanation",
];
/// GRCh38's chr1 to chr22 and chrX, which the painting of each child's
/// stretches scales by (UCSC hg38 chrom.sizes).
pub const GRCH38: [(&str, i64); 23] = [
    ("chr1", 248956422),
    ("chr2", 242193529),
    ("chr3", 198295559),
    ("chr4", 190214555),
    ("chr5", 181538259),
    ("chr6", 170805979),
    ("chr7", 159345973),
    ("chr8", 145138636),
    ("chr9", 138394717),
    ("chr10", 133797422),
    ("chr11", 135086622),
    ("chr12", 133275309),
    ("chr13", 114364328),
    ("chr14", 107043718),
    ("chr15", 101991189),
    ("chr16", 90338345),
    ("chr17", 83257441),
    ("chr18", 80373285),
    ("chr19", 58617616),
    ("chr20", 64444167),
    ("chr21", 46709983),
    ("chr22", 50818468),
    ("chrX", 156040895),
];
/// The groups of the side-by-side table, in its order.
const GROUPS: [&str; 4] = ["Medicines", "HLA", "Haplogroups", "Blood groups"];

/// A file of the family that exists but cannot be read: family.json is
/// then not built. The text names the kind of error, never a path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Unreadable(pub String);

impl Unreadable {
    fn of(e: &io::Error) -> Unreadable {
        Unreadable(
            match e.kind() {
                io::ErrorKind::PermissionDenied => "PermissionError",
                io::ErrorKind::NotFound => "FileNotFoundError",
                io::ErrorKind::IsADirectory => "IsADirectoryError",
                _ => "OSError",
            }
            .into(),
        )
    }
}

// --- paths ------------------------------------------------------------------------------

/// `os.path.realpath`: links followed as far as they go, `..` taken after
/// them; a part that does not exist is kept as written.
pub fn realpath(p: &Path) -> PathBuf {
    fn resolve(p: &Path, depth: usize) -> PathBuf {
        let mut out = PathBuf::from("/");
        for c in p.components() {
            match c {
                Component::Prefix(_) | Component::RootDir | Component::CurDir => {}
                Component::ParentDir => {
                    out.pop();
                }
                Component::Normal(name) => {
                    let cand = out.join(name);
                    let link =
                        fs::symlink_metadata(&cand).is_ok_and(|m| m.file_type().is_symlink());
                    out = match fs::read_link(&cand) {
                        Ok(t) if link && depth < 40 => resolve(&out.join(t), depth + 1),
                        _ => cand,
                    };
                }
            }
        }
        out
    }
    resolve(&absolute(p), 0)
}

/// `os.path.abspath`: against the working folder, `.` and `..` taken by
/// their names alone.
fn absolute(p: &Path) -> PathBuf {
    let joined = if p.is_absolute() {
        p.to_path_buf()
    } else {
        std::env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("/"))
            .join(p)
    };
    let mut out = PathBuf::from("/");
    for c in joined.components() {
        match c {
            Component::ParentDir => {
                out.pop();
            }
            Component::Normal(n) => out.push(n),
            _ => {}
        }
    }
    out
}

/// Whether `path`, links followed, is `top` or inside it.
pub fn inside(path: &Path, top: &Path) -> bool {
    realpath(path).starts_with(realpath(top))
}

fn is_file(p: &Path) -> bool {
    fs::metadata(p).is_ok_and(|m| m.is_file())
}

/// A folder holding `members.tsv` or `pedigree.tsv`.
pub fn is_family(d: &Path) -> bool {
    d.is_dir() && layout::MARKERS.iter().any(|m| is_file(&d.join(m)))
}

/// A run name: a letter or digit, then letters, digits, `.`, `_` and `-`.
pub fn valid_name(n: &str) -> bool {
    let mut cs = n.chars();
    cs.next().is_some_and(|c| c.is_ascii_alphanumeric())
        && cs.all(|c| c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-'))
}

/// The family folders under `root`, `(name, dir)` sorted by name.
pub fn family_dirs(root: &Path) -> Vec<(String, PathBuf)> {
    let Ok(rd) = fs::read_dir(root) else {
        return Vec::new();
    };
    let mut names: Vec<String> = rd
        .filter_map(|e| e.ok()?.file_name().into_string().ok())
        .collect();
    names.sort();
    names
        .into_iter()
        .filter_map(|n| {
            let d = root.join(&n);
            (!n.starts_with('.') && valid_name(&n) && inside(&d, root) && is_family(&d))
                .then_some((n, d))
        })
        .collect()
}

// --- the family's files -----------------------------------------------------------------

/// A row of a tab-separated file by its header's names.
pub type Row = IndexMap<String, String>;

/// The rows of a tab-separated file with a header line, comments and blank
/// lines skipped, short rows padded; none when it is not there or is a link
/// out of its own folder (the page serves only what the family folder
/// itself holds).
pub fn read_tsv(path: &Path) -> Result<Vec<Row>, Unreadable> {
    let dir = absolute(path)
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_default();
    if !is_file(path) || !inside(path, &dir) {
        return Ok(Vec::new());
    }
    let f = fs::File::open(path).map_err(|e| Unreadable::of(&e))?;
    let mut lines = Lines::new(BufReader::new(f), Utf8::Replace);
    let mut line = String::new();
    let head: Vec<String> = if lines.read_line(&mut line).map_err(|e| Unreadable::of(&e))? {
        line.trim_end_matches(['\r', '\n'])
            .trim_start_matches('#')
            .split('\t')
            .map(str::to_owned)
            .collect()
    } else {
        vec![String::new()]
    };
    let mut rows = Vec::new();
    while lines.read_line(&mut line).map_err(|e| Unreadable::of(&e))? {
        let l = line.trim_end_matches(['\r', '\n']);
        if pycompat::text::strip(l, None).is_empty() || l.starts_with('#') {
            continue;
        }
        let vals: Vec<&str> = l.split('\t').collect();
        let mut r = Row::new();
        for (i, h) in head.iter().enumerate() {
            r.insert(h.clone(), vals.get(i).copied().unwrap_or("").to_owned());
        }
        rows.push(r);
    }
    Ok(rows)
}

fn cell<'a>(r: &'a Row, k: &str) -> Option<&'a str> {
    r.get(k).map(String::as_str)
}

fn s(v: &str) -> Value {
    Value::Str(v.to_owned())
}

fn opt(v: Option<&str>) -> Value {
    v.map_or(Value::Null, s)
}

fn obj<const N: usize>(pairs: [(&str, Value); N]) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

fn keep(r: &Row, keys: &[&str]) -> Value {
    Value::Object(
        keys.iter()
            .map(|k| ((*k).to_owned(), opt(cell(r, k))))
            .collect(),
    )
}

fn as_obj(v: &mut Value) -> &mut IndexMap<String, Value> {
    match v {
        Value::Object(o) => o,
        _ => unreachable!("built as an object"),
    }
}

/// The people of the family folder, or why the pedigree cannot be read.
pub fn members_of(d: &Path) -> Result<People, String> {
    let (people, _) = layout::load(d)?;
    if people.is_empty() {
        return Err("the pedigree names no one".into());
    }
    Ok(people)
}

/// `(a, b, couple_<a>_<b>)` of each pair of partners, `a` and `b` sorted by
/// name as step 39 names them.
pub fn couple_names(people: &People) -> Vec<(String, String, String)> {
    layout::couples(people)
        .into_iter()
        .filter(|(a, b)| {
            people[a].partner.as_deref() == Some(b.as_str())
                || people[b].partner.as_deref() == Some(a.as_str())
        })
        .map(|(a, b)| {
            let (x, y) = if a <= b { (a, b) } else { (b, a) };
            let base = format!("couple_{x}_{y}");
            (x, y, base)
        })
        .collect()
}

/// The file names of the family folder the page may open: those that are
/// there, inside the folder.
pub fn served_files(d: &Path, people: &People) -> Vec<String> {
    let mut names: Vec<String> = FIXED.iter().map(|f| (*f).to_owned()).collect();
    names.extend(people.keys().map(|s| format!("{s}_segments.svg")));
    for (_, _, base) in couple_names(people) {
        names.push(format!("{base}.md"));
        names.push(format!("{base}.tsv"));
    }
    let mut out: Vec<String> = names
        .into_iter()
        .collect::<HashSet<_>>()
        .into_iter()
        .filter(|n| is_file(&d.join(n)) && inside(&d.join(n), d))
        .collect();
    out.sort();
    out
}

fn people_rows(people: &People) -> Vec<Value> {
    people
        .iter()
        .map(|(s_, m)| {
            obj([
                ("sample", s(s_)),
                ("sex", s(&m.sex)),
                ("father", opt(m.father.as_deref())),
                ("mother", opt(m.mother.as_deref())),
                ("partner", opt(m.partner.as_deref())),
                ("sequenced", Value::Bool(m.sequenced)),
            ])
        })
        .collect()
}

/// state.json's `plugins.family` on the family's page.
pub fn family_state(name: &str, d: &Path) -> Result<Value, Unreadable> {
    let people = match members_of(d) {
        Ok(p) => p,
        Err(e) => {
            return Ok(obj([
                ("kind", s("family")),
                (
                    "error",
                    Value::Str(format!("the pedigree could not be read: {e}")),
                ),
            ]));
        }
    };
    let files = served_files(d, &people);
    let url = format!("/family/{name}/");
    let link = |f: &str| {
        if files.iter().any(|x| x == f) {
            Value::Str(format!("{url}{f}"))
        } else {
            Value::Null
        }
    };
    let children: Vec<Value> = people
        .iter()
        .filter(|(_, m)| m.father.is_some() || m.mother.is_some())
        .map(|(c, m)| {
            obj([
                ("child", s(c)),
                ("father", opt(m.father.as_deref())),
                ("mother", opt(m.mother.as_deref())),
                ("segments", link(&format!("{c}_segments.svg"))),
            ])
        })
        .collect();
    let rel = read_tsv(&d.join("relationships.tsv"))?
        .iter()
        .map(|r| {
            obj([
                ("a", opt(cell(r, "sample_a"))),
                ("b", opt(cell(r, "sample_b"))),
                ("declared", opt(cell(r, "declared"))),
                ("verdict", opt(cell(r, "verdict"))),
                ("reason", opt(cell(r, "reason"))),
            ])
        })
        .collect();
    let ibd = read_tsv(&d.join("ibd_pairs.tsv"))?
        .iter()
        .map(|r| {
            obj([
                ("a", opt(cell(r, "sample_a"))),
                ("b", opt(cell(r, "sample_b"))),
                ("declared", opt(cell(r, "declared"))),
                ("shared_cm", opt(cell(r, "shared_cm"))),
                ("degree", opt(cell(r, "degree"))),
            ])
        })
        .collect();
    let mut couples = Vec::new();
    for (a, b, base) in couple_names(&people) {
        let rows = read_tsv(&d.join(format!("{base}.tsv")))?;
        let count = |w: &str| rows.iter().filter(|r| cell(r, "result") == Some(w)).count() as i128;
        let md = format!("{base}.md");
        couples.push(obj([
            ("a", s(&a)),
            ("b", s(&b)),
            ("report", link(&md)),
            ("table", link(&format!("{base}.tsv"))),
            (
                "screened",
                Value::Bool(!rows.is_empty() || files.contains(&md)),
            ),
            ("shared", Value::Int(count("shared"))),
            ("x_linked", Value::Int(count("x_linked"))),
        ]));
    }
    let plan = read_tsv(&d.join("family_plan.tsv"))?
        .iter()
        .map(|r| {
            obj([
                ("analysis", opt(cell(r, "analysis"))),
                ("members", opt(cell(r, "members"))),
                ("decision", opt(cell(r, "decision"))),
                ("reason", opt(cell(r, "reason"))),
            ])
        })
        .collect();
    let reports = files
        .iter()
        .filter(|f| f.ends_with(".md"))
        .map(|f| obj([("name", s(f)), ("url", Value::Str(format!("{url}{f}")))]))
        .collect();
    Ok(obj([
        ("kind", s("family")),
        ("name", s(name)),
        ("pedigree", Value::Str(format!("{url}pedigree.svg"))),
        ("members", Value::Array(people_rows(&people))),
        ("children", Value::Array(children)),
        ("relationships", Value::Array(rel)),
        ("ibd_pairs", Value::Array(ibd)),
        ("couples", Value::Array(couples)),
        ("plan", Value::Array(plan)),
        ("reports", Value::Array(reports)),
        ("files", Value::Array(files.iter().map(|f| s(f)).collect())),
    ]))
}

/// The pedigree drawing's geometry as data, so the page draws each person
/// as a link.
pub fn layout_value(people: &People) -> Value {
    let l = layout::page_layout(people);
    obj([
        ("width", Value::Int(l.width.into())),
        ("height", Value::Int(l.height.into())),
        ("size", Value::Int(l.size.into())),
        (
            "people",
            Value::Array(
                l.people
                    .iter()
                    .map(|p| {
                        obj([
                            ("sample", s(&p.sample)),
                            ("x", Value::Float(p.x)),
                            ("y", Value::Float(p.y)),
                            ("gen", Value::Int(p.generation.into())),
                            ("sex", s(&p.sex)),
                            ("sequenced", Value::Bool(p.sequenced)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "lines",
            Value::Array(
                l.lines
                    .iter()
                    .map(|ln| Value::Array(ln.iter().map(|v| Value::Float(*v)).collect()))
                    .collect(),
            ),
        ),
    ])
}

/// Step 38 writes the verdict "ok" for a pair that matches its declared
/// relation. The page says that in words: no word on the page reads as an
/// all-clear (state.json keeps the file's own value).
pub fn verdict_words(v: &Value) -> Value {
    match v {
        Value::Str(x) if pycompat::text::lower(pycompat::text::strip(x, None)) == "ok" => {
            s("matches the declared relation")
        }
        _ => v.clone(),
    }
}

fn with_verdict_words(rows: &mut Value) {
    if let Value::Array(rows) = rows {
        for r in rows {
            let o = as_obj(r);
            let w = verdict_words(o.get("verdict").unwrap_or(&Value::Null));
            o.insert("verdict".into(), w);
        }
    }
}

/// family.json without what only the site knows (the member runs): the
/// state, the layout, the full pair rows and couple rows, each child's
/// stretches and the trio tables.
pub fn family_page(name: &str, d: &Path) -> Result<Value, Unreadable> {
    let st = family_state(name, d)?;
    let mut out = obj([
        ("state", st.clone()),
        ("layout", Value::Null),
        ("relationships", Value::Array(Vec::new())),
        ("ibd_pairs", Value::Array(Vec::new())),
        ("couples", Value::Object(IndexMap::new())),
        ("couple_note", s(og_render::notes::COUPLE_NOTE)),
        ("segments", Value::Object(IndexMap::new())),
        ("de_novo", Value::Array(Vec::new())),
        ("compound_het", Value::Array(Vec::new())),
        ("trio_files", Value::Object(IndexMap::new())),
        ("side_by_side", Value::Null),
    ]);
    if get(&st, "error").ok().flatten().is_some() {
        return Ok(out);
    }
    // read again, as the state just did: a pedigree that stopped reading in between is unreadable
    let people = members_of(d).map_err(Unreadable)?;
    let o = as_obj(&mut out);
    o.insert("layout".into(), layout_value(&people));
    let rows = |f: &str, keys: &[&str]| -> Result<Value, Unreadable> {
        Ok(Value::Array(
            read_tsv(&d.join(f))?
                .iter()
                .map(|r| keep(r, keys))
                .collect(),
        ))
    };
    o.insert(
        "relationships".into(),
        rows("relationships.tsv", &REL_KEYS)?,
    );
    o.insert("ibd_pairs".into(), rows("ibd_pairs.tsv", &IBD_KEYS)?);
    let mut couples = IndexMap::new();
    for (a, b, base) in couple_names(&people) {
        let mut rs: Vec<Row> = read_tsv(&d.join(format!("{base}.tsv")))?;
        rs.sort_by_key(|r| match cell(r, "result") {
            Some("shared") => 0,
            Some("x_linked") => 1,
            Some("one_partner") => 2,
            _ => 3,
        });
        couples.insert(
            format!("{a} {b}"),
            Value::Array(rs.iter().map(|r| keep(r, &COUPLE_KEYS)).collect()),
        );
    }
    o.insert("couples".into(), Value::Object(couples));
    for (k, v) in child_depth(d, &people)? {
        o.insert(k, v);
    }
    o.insert(
        "chrom_lengths".into(),
        Value::Array(
            GRCH38
                .iter()
                .map(|(c, n)| Value::Array(vec![s(c), Value::Int((*n).into())]))
                .collect(),
        ),
    );
    let mut state = st;
    if let Some(r) = as_obj(&mut state).get_mut("relationships") {
        with_verdict_words(r);
    }
    o.insert("state".into(), state);
    if let Some(r) = o.get_mut("relationships") {
        with_verdict_words(r);
    }
    Ok(out)
}

/// `here`, `rejected` (a link out of the family folder: it is not read, and
/// the page never shows where it points) or `missing`.
pub fn file_state(d: &Path, f: &str) -> &'static str {
    let p = d.join(f);
    if fs::symlink_metadata(&p).is_err() {
        return "missing";
    }
    if !inside(&p, &absolute(d)) {
        return "rejected";
    }
    if is_file(&p) { "here" } else { "missing" }
}

/// Each child's stretches (step 38's `<child>_segments.tsv`) and the trio
/// tables, as the steps wrote them: the compound-het pairs in trans first.
fn child_depth(d: &Path, people: &People) -> Result<Vec<(String, Value)>, Unreadable> {
    let mut segs = IndexMap::new();
    for (c, m) in people {
        if m.father.is_none() && m.mother.is_none() {
            continue;
        }
        let rows = read_tsv(&d.join(format!("{c}_segments.tsv")))?;
        if !rows.is_empty() {
            segs.insert(
                c.clone(),
                Value::Array(rows.iter().map(|r| keep(r, &SEGMENT_KEYS)).collect()),
            );
        }
    }
    let files: IndexMap<&str, &str> = TRIO_FILES.iter().map(|f| (*f, file_state(d, f))).collect();
    let mut comp = read_tsv(&d.join("compound_het.tsv"))?;
    comp.sort_by_key(|r| match cell(r, "phase") {
        Some("in trans") => 0,
        Some("in cis") => 2,
        _ => 1,
    });
    let de_novo = read_tsv(&d.join("de_novo.tsv"))?;
    let status = trio_status(d, people, &files)?;
    Ok(vec![
        ("segments".into(), Value::Object(segs)),
        (
            "de_novo".into(),
            Value::Array(de_novo.iter().map(|r| keep(r, &DE_NOVO_KEYS)).collect()),
        ),
        (
            "compound_het".into(),
            Value::Array(comp.iter().map(|r| keep(r, &COMP_KEYS)).collect()),
        ),
        (
            "trio_files".into(),
            Value::Object(
                files
                    .iter()
                    .map(|(f, st)| ((*f).to_owned(), Value::Bool(*st == "here")))
                    .collect(),
            ),
        ),
        (
            "trio_file_state".into(),
            Value::Object(
                files
                    .iter()
                    .map(|(f, st)| ((*f).to_owned(), s(st)))
                    .collect(),
            ),
        ),
        ("trio_status".into(), status),
    ])
}

fn line(assessed: Option<bool>, why: String) -> Value {
    obj([
        ("assessed", assessed.map_or(Value::Null, Value::Bool)),
        ("why", Value::Str(why)),
    ])
}

/// Per child, whether step 38 assessed it for de novo candidates and
/// compound-het pairs, and why not. The rule is step 38's (`og family
/// report`): de novo candidates need both parents sequenced, compound-het
/// pairs a sequenced parent; a file that is not here, or is a link out of
/// the folder, means no child was assessed.
fn trio_status(
    d: &Path,
    people: &People,
    files: &IndexMap<&str, &str>,
) -> Result<Value, Unreadable> {
    let mut plan: HashMap<Option<String>, Row> = HashMap::new();
    for r in read_tsv(&d.join("family_plan.tsv"))? {
        plan.insert(cell(&r, "analysis").map(str::to_owned), r);
    }
    let pl = |k: &str| plan.get(&Some(k.to_owned()));
    let filled = |r: Option<&Row>, k: &str| {
        r.and_then(|r| cell(r, k))
            .filter(|v| !v.is_empty())
            .map(str::to_owned)
    };
    let trio = pl("trio_phasing");
    let absent = |f: &str| {
        if files[f] == "rejected" {
            return format!(
                "{f} is a link out of the family folder, so this page does not read it"
            );
        }
        let mut why = format!("{f} is not here");
        if let Some(dec) = filled(trio, "decision") {
            let reason = filled(trio, "reason")
                .map(|r| format!(" ({r})"))
                .unwrap_or_default();
            why.push_str(&format!(": trio phasing {dec}{reason}"));
        }
        why
    };
    let joint = pl("joint_call");
    let called: Option<HashSet<String>> = (joint.and_then(|r| cell(r, "decision")) == Some("runs"))
        .then(|| {
            joint
                .and_then(|r| cell(r, "members"))
                .unwrap_or("")
                .split(',')
                .map(str::to_owned)
                .collect()
        });
    let seq = |p: Option<&str>| {
        p.filter(|p| !p.is_empty())
            .and_then(|p| people.get(p).map(|m| (p, m)))
            .is_some_and(|(p, m)| m.sequenced && called.as_ref().is_none_or(|c| c.contains(p)))
    };
    let need = if called.is_none() {
        "sequenced"
    } else {
        "sequenced and in the joint call"
    };
    let not_called = |p: &str| {
        if people[p].sequenced {
            format!("{p} is not in the joint call (family_plan.tsv)")
        } else {
            format!("{p} is not sequenced")
        }
    };
    let gaps = |m: &layout::Person| {
        [
            ("father", m.father.as_deref()),
            ("mother", m.mother.as_deref()),
        ]
        .into_iter()
        .filter(|(_, p)| !seq(*p))
        .map(|(w, p)| match p {
            None => format!("no {w} is named"),
            Some(p) => format!("{p} (the {w}) is not"),
        })
        .collect::<Vec<_>>()
        .join("; ")
    };
    let mut out = IndexMap::new();
    for (c, m) in people {
        if m.father.is_none() && m.mother.is_none() {
            continue;
        }
        let de_novo = if files["de_novo.tsv"] != "here" {
            line(Some(false), absent("de_novo.tsv"))
        } else if !seq(Some(c)) {
            line(Some(false), not_called(c))
        } else if seq(m.father.as_deref()) && seq(m.mother.as_deref()) {
            line(Some(true), "both parents are sequenced".into())
        } else {
            line(
                Some(false),
                format!("de novo candidates need both parents {need}: {}", gaps(m)),
            )
        };
        let comp = if files["compound_het.tsv"] != "here" {
            line(Some(false), absent("compound_het.tsv"))
        } else if !seq(Some(c)) {
            line(Some(false), not_called(c))
        } else if seq(m.father.as_deref()) || seq(m.mother.as_deref()) {
            line(
                Some(true),
                "a parent is sequenced; pairs come only from the child's clinical-filter candidates (step 23), so \
                 with none given no pair is listed"
                    .into(),
            )
        } else {
            line(
                Some(false),
                format!("compound-het pairs need a parent {need}"),
            )
        };
        out.insert(
            c.clone(),
            obj([("de_novo", de_novo), ("compound_het", comp)]),
        );
    }
    Ok(Value::Object(out))
}

/// Whether the clinical filter (step 23) wrote the child's candidates, the
/// only variants step 38 pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Given {
    /// There, and not newer than `compound_het.tsv`.
    Here,
    /// The run is here without them (a file linking out of the run counts
    /// as absent).
    Absent,
    /// Written after `compound_het.tsv`, so step 38 may not have had them.
    Newer,
}

/// [`Given`] of `sample`'s run folder `d`, read only inside `root`; `None`
/// when the run is not on this site, so the page cannot tell.
pub fn candidates_given(
    root: &Path,
    d: Option<&Path>,
    sample: &str,
    since: Option<SystemTime>,
) -> Option<Given> {
    let d = d?;
    let p = d.join("clinical").join(format!("{sample}_clinical.vcf.gz"));
    if !(is_file(&p) && inside(&p, d) && inside(&p, root)) {
        return Some(Given::Absent);
    }
    let newer = since.is_some_and(|t| {
        fs::metadata(&p)
            .and_then(|m| m.modified())
            .is_ok_and(|m| m > t)
    });
    Some(if newer { Given::Newer } else { Given::Here })
}

/// The trio status with each child's compound-het line told by its step 23
/// candidates. `assessed` is null where the page cannot tell.
pub fn with_candidates(status: &mut Value, given: &HashMap<String, Option<Given>>) {
    let Value::Object(status) = status else {
        return;
    };
    for (child, st) in status.iter_mut() {
        let st = as_obj(st);
        let assessed = st
            .get("compound_het")
            .and_then(|c| get(c, "assessed").ok().flatten())
            .is_some_and(truthy);
        if !assessed {
            continue;
        }
        let new = match given.get(child).copied().flatten() {
            Some(Given::Here) => line(
                Some(true),
                "a parent is sequenced and the child's clinical-filter candidates (step 23) are here; pairs come \
                 only from them"
                    .into(),
            ),
            Some(Given::Absent) => line(
                Some(false),
                format!(
                    "the clinical filter (step 23) wrote no candidates for {child} (clinical/{child}_clinical.vcf.gz), \
                     so step 38 had no variant to pair"
                ),
            ),
            Some(Given::Newer) => line(
                None,
                format!(
                    "{child}'s clinical-filter candidates (step 23) are newer than compound_het.tsv, so step 38 may \
                     not have paired them"
                ),
            ),
            None => line(
                None,
                format!(
                    "pairs come only from the child's clinical-filter candidates (step 23), and {child}'s run is not \
                     on this site to tell whether step 23 wrote them"
                ),
            ),
        };
        st.insert("compound_het".into(), new);
    }
}

/// A member run's summary.json, or `None`: none yet, unreadable, not an
/// object, or a link out of the run folder or out of `root`.
pub fn member_summary(root: &Path, d: Option<&Path>) -> Option<Value> {
    let d = d?;
    let p = d.join("summary.json");
    if !is_file(&p) || !inside(&p, d) || !inside(&p, root) {
        return None;
    }
    let mut text = String::new();
    fs::File::open(&p).ok()?.read_to_string(&mut text).ok()?;
    match loads(&text) {
        Ok(v @ Value::Object(_)) => Some(v),
        _ => None,
    }
}

// --- the members side by side -----------------------------------------------------------

/// A summary of a shape the report never writes: that member's column says
/// so, and only that column.
#[derive(Debug)]
struct Malformed;

type Got<T> = Result<T, Malformed>;

fn empty_obj() -> &'static Value {
    static E: std::sync::OnceLock<Value> = std::sync::OnceLock::new();
    E.get_or_init(|| Value::Object(IndexMap::new()))
}

/// `d.get(k)`, `None` for a missing key; a `d` that is not an object is
/// malformed.
fn at<'a>(d: &'a Value, k: &str) -> Got<&'a Value> {
    static NULL: Value = Value::Null;
    Ok(get(d, k).map_err(|_| Malformed)?.unwrap_or(&NULL))
}

/// `v or {}`.
fn or_obj(v: &Value) -> &Value {
    if truthy(v) { v } else { empty_obj() }
}

/// What `for x in v or []` walks: a list's items, a string's characters, an
/// object's keys; anything else that is true cannot be walked.
fn walk(v: &Value) -> Got<Vec<Value>> {
    if !truthy(v) {
        return Ok(Vec::new());
    }
    match v {
        Value::Array(a) => Ok(a.clone()),
        Value::Str(x) => Ok(x.chars().map(|c| Value::Str(c.to_string())).collect()),
        Value::Object(o) => Ok(o.keys().map(|k| s(k)).collect()),
        _ => Err(Malformed),
    }
}

fn text(v: &Value) -> Got<String> {
    py_str(v).map_err(|_| Malformed)
}

/// A section's data while it is complete or stale (`{}` when it has none),
/// else `None`.
fn section(summary: &Value, key: &str) -> Got<Option<Value>> {
    let sec = or_obj(at(or_obj(at(summary, "sections")?), key)?);
    let state = at(sec, "state")?;
    let data = at(sec, "data")?;
    Ok(matches!(state, Value::Str(x) if x == "ok" || x == "stale").then(|| or_obj(data).clone()))
}

fn stale(summary: &Value, key: &str) -> Got<bool> {
    let sec = or_obj(at(or_obj(at(summary, "sections")?), key)?);
    Ok(matches!(at(sec, "state")?, Value::Str(x) if x == "stale"))
}

/// One cell of the side-by-side table: group, label, value, stale.
type Cell = (String, Value, Value, bool);

/// The cells of one member's summary.json, as its report prints them;
/// nothing computed. A value is stale when its section is older than its
/// input (state `stale`).
fn side_values(summary: Option<&Value>) -> Got<Vec<Cell>> {
    let summary = summary.map_or(empty_obj(), or_obj);
    let mut out: Vec<Cell> = Vec::new();
    let c = or_obj(at(summary, "cyp2d6")?);
    let calls_v = or_obj(at(c, "calls")?);
    let Value::Object(calls_o) = calls_v else {
        return Err(Malformed);
    };
    let mut calls = Vec::new();
    for (k, v) in calls_o {
        if truthy(v) {
            calls.push(format!("{k} {}", text(v)?));
        }
    }
    let calls = calls.join(", ");
    if !calls.is_empty() {
        let agree = at(c, "agree")?;
        let words = if truthy(agree) {
            " (callers agree)"
        } else if *agree == Value::Bool(false) {
            " (callers differ)"
        } else {
            ""
        };
        // summary.json's cyp2d6 block reads pypgx and Cyrius only when complete; its PharmCAT call is cpic's
        let old = truthy(at(calls_v, "PharmCAT")?) && stale(summary, "cpic")?;
        out.push((
            "Medicines".into(),
            s("CYP2D6 (callers)"),
            Value::Str(calls + words),
            old,
        ));
    }
    let cpic = section(summary, "cpic")?;
    for g in walk(at(cpic.as_ref().map_or(empty_obj(), or_obj), "genes")?)? {
        let ph = at(&g, "phenotype")?;
        let dip = at(&g, "diplotype")?;
        let dip = if truthy(dip) { text(dip)? } else { ".".into() };
        let ph = if truthy(ph) {
            format!(" \"{}\"", text(ph)?)
        } else {
            String::new()
        };
        let label = at(&g, "gene")?.clone();
        out.push((
            "Medicines".into(),
            label,
            Value::Str(dip + &ph),
            stale(summary, "cpic")?,
        ));
    }
    let hla = section(summary, "hla")?;
    for l in walk(at(hla.as_ref().map_or(empty_obj(), or_obj), "loci")?)? {
        let alleles = at(&l, "alleles")?;
        let mut names = Vec::new();
        for a in walk(alleles)? {
            match a {
                Value::Str(x) => names.push(x),
                _ => return Err(Malformed),
            }
        }
        let typed = names.join(" ");
        let value = if typed.is_empty() {
            s("not typed")
        } else {
            Value::Str(typed)
        };
        out.push((
            "HLA".into(),
            at(&l, "gene")?.clone(),
            value,
            stale(summary, "hla")?,
        ));
    }
    let dot = s(".");
    let or_dot = |v: &Value| if truthy(v) { v.clone() } else { dot.clone() };
    if let Some(mt) = section(summary, "haplogroup")? {
        out.push((
            "Haplogroups".into(),
            s("Mitochondrial"),
            or_dot(at(&mt, "haplogroup")?),
            stale(summary, "haplogroup")?,
        ));
    }
    let tree = section(summary, "y_tree")?;
    let ykey = if tree.as_ref().is_some_and(truthy) {
        "y_tree"
    } else {
        "y_haplogroup"
    };
    let y = match tree {
        Some(t) if truthy(&t) => Some(t),
        _ => section(summary, "y_haplogroup")?,
    };
    if let Some(y) = y {
        out.push((
            "Haplogroups".into(),
            s("Y chromosome"),
            or_dot(at(&y, "haplogroup")?),
            stale(summary, ykey)?,
        ));
    }
    if let Some(bt) = section(summary, "bloodtype")? {
        let old = stale(summary, "bloodtype")?;
        let mut systems = Vec::new();
        for x in walk(at(&bt, "systems")?)? {
            if truthy(at(&x, "system")?) {
                systems.push(x);
            }
        }
        // {system: phenotype}, a later row of one system winning
        let mut by: Vec<(Value, Value)> = Vec::new();
        for x in &systems {
            let k = at(x, "system")?.clone();
            if matches!(k, Value::Array(_) | Value::Object(_)) {
                return Err(Malformed);
            }
            let v = or_dot(at(x, "phenotype")?);
            match by.iter_mut().find(|(o, _)| same_key(o, &k)) {
                Some(e) => e.1 = v,
                None => by.push((k, v)),
            }
        }
        let look = |name: &str| {
            by.iter()
                .find(|(k, _)| same_key(k, &s(name)))
                .map_or(dot.clone(), |e| e.1.clone())
        };
        let named = |k: &str, name: &str| -> Got<Value> {
            let v = at(&bt, k)?;
            Ok(if truthy(v) { v.clone() } else { look(name) })
        };
        out.push(("Blood groups".into(), s("ABO"), named("abo", "ABO")?, old));
        out.push(("Blood groups".into(), s("RhD"), named("rhd", "RhD")?, old));
        // every other system step 41 typed, in its order (Kell, Duffy, Kidd, MNS...)
        for x in &systems {
            let k = at(x, "system")?;
            if !matches!(k, Value::Str(v) if v == "ABO" || v == "RhD") {
                let v = by
                    .iter()
                    .find(|(o, _)| same_key(o, k))
                    .map_or(dot.clone(), |e| e.1.clone());
                out.push(("Blood groups".into(), k.clone(), v, old));
            }
        }
    }
    Ok(out)
}

/// Whether two labels are the same key, as Python compares them: numbers
/// by value (`True` is 1), strings by text.
fn same_key(a: &Value, b: &Value) -> bool {
    let num = |v: &Value| match v {
        Value::Bool(x) => Some(f64::from(u8::from(*x))),
        Value::Int(i) => Some(*i as f64),
        Value::Float(x) => Some(*x),
        _ => None,
    };
    match (a, b) {
        (Value::Null, Value::Null) => true,
        (Value::Str(x), Value::Str(y)) => x == y,
        _ => matches!((num(a), num(b)), (Some(x), Some(y)) if x == y),
    }
}

/// Step 41's own order of the blood group systems (`og bloodtype systems`).
pub fn step41_order() -> Vec<String> {
    og_bloodtype::rules::SYSTEMS
        .iter()
        .map(|(name, _)| name.clone())
        .collect()
}

/// The members' results in one table: `{members, rows: [{group, label,
/// values, stale}]}`. A member with no summary.json yet says so in its
/// column; a row a member's summary does not have is `.`. `stale` holds,
/// per member, whether that value comes from a section older than its
/// input. A summary whose fields are not of the shape the report writes
/// costs its own column only: it says so there, and the page still draws.
pub fn side_by_side(names: &[String], summaries: &HashMap<String, Option<Value>>) -> Value {
    side_by_side_in(names, summaries, &step41_order())
}

/// [`side_by_side`] with the blood group systems in `step41`'s order.
pub fn side_by_side_in(
    names: &[String],
    summaries: &HashMap<String, Option<Value>>,
    step41: &[String],
) -> Value {
    let mut per: HashMap<&str, Vec<Cell>> = HashMap::new();
    let mut bad: HashSet<&str> = HashSet::new();
    for n in names {
        let sum = summaries.get(n).and_then(Option::as_ref);
        match side_values(sum) {
            // a label that is a list or an object cannot key the table's rows
            Ok(got)
                if !got
                    .iter()
                    .any(|c| matches!(c.1, Value::Array(_) | Value::Object(_))) =>
            {
                per.insert(n, got);
            }
            _ => {
                per.insert(n, Vec::new());
                bad.insert(n);
            }
        }
    }
    let same = |a: &(String, Value), g: &str, l: &Value| a.0 == g && same_key(&a.1, l);
    // Each member's rows keep their own order (step 41's for the blood groups): a row first seen in a later
    // member goes right after the row before it in that member's list, not at the end.
    let mut keys: Vec<(String, Value)> = Vec::new();
    for n in names {
        let mut prev: Option<(String, Value)> = None;
        for (g, label, _, _) in &per[n.as_str()] {
            if !keys.iter().any(|k| same(k, g, label)) {
                let at = match &prev {
                    Some(p) => keys
                        .iter()
                        .position(|k| same(k, &p.0, &p.1))
                        .map_or(keys.len(), |i| i + 1),
                    None => keys.len(),
                };
                keys.insert(at, (g.clone(), label.clone()));
            }
            prev = Some((g.clone(), label.clone()));
        }
    }
    let rank41 = |l: &Value| match l {
        Value::Str(x) => step41.iter().position(|s_| s_ == x).unwrap_or(step41.len()),
        _ => step41.len(),
    };
    let mut ranked: Vec<(usize, usize, usize, (String, Value))> = keys
        .into_iter()
        .enumerate()
        .map(|(i, k)| {
            let g = GROUPS.iter().position(|x| *x == k.0).unwrap_or(9);
            let b = if k.0 == "Blood groups" {
                rank41(&k.1)
            } else {
                0
            };
            (g, b, i, k)
        })
        .collect();
    ranked.sort_by_key(|r| (r.0, r.1, r.2));
    let rows: Vec<Value> = ranked
        .into_iter()
        .map(|(_, _, _, (g, label))| {
            let mut values = Vec::new();
            let mut stale = Vec::new();
            for n in names {
                let (v, st) = if bad.contains(n.as_str()) {
                    (s("summary not readable"), false)
                } else if summaries.get(n).and_then(Option::as_ref).is_none() {
                    (s("no summary yet"), false)
                } else {
                    match per[n.as_str()]
                        .iter()
                        .find(|c| c.0 == g && same_key(&c.1, &label))
                    {
                        Some(c) => (c.2.clone(), c.3),
                        None => (s("."), false),
                    }
                };
                values.push(v);
                stale.push(Value::Bool(st));
            }
            obj([
                ("group", Value::Str(g)),
                ("label", label),
                ("values", Value::Array(values)),
                ("stale", Value::Array(stale)),
            ])
        })
        .collect();
    obj([
        (
            "members",
            Value::Array(names.iter().map(|n| s(n)).collect()),
        ),
        ("rows", Value::Array(rows)),
    ])
}

// --- with what the site knows -----------------------------------------------------------

/// A run folder by member name, from the input's `candidates` or
/// `summaries`; a member left out (or null) has no run on this site.
fn dirs(input: &Value, key: &str) -> Option<HashMap<String, Option<PathBuf>>> {
    let Value::Object(m) = get(input, key).ok().flatten()? else {
        return None;
    };
    Some(
        m.iter()
            .map(|(k, v)| {
                (
                    k.clone(),
                    if let Value::Str(d) = v {
                        Some(PathBuf::from(d))
                    } else {
                        None
                    },
                )
            })
            .collect(),
    )
}

/// family.json as the site serves it: [`family_page`], then the member runs
/// of `input` (`runs`, copied as they are), each child's compound-het line
/// told by its run's step 23 candidates (`candidates`), and the members
/// side by side from each sequenced member's summary.json (`summaries`).
/// A run folder is read only inside `root`.
pub fn site_page(name: &str, d: &Path, root: &Path, input: &Value) -> Result<Value, Unreadable> {
    let mut out = family_page(name, d)?;
    let root = realpath(root);
    let o = as_obj(&mut out);
    if let Some(runs) = get(input, "runs").ok().flatten() {
        o.insert("runs".into(), runs.clone());
    }
    if let Some(given_dirs) = dirs(input, "candidates") {
        let comp = d.join("compound_het.tsv");
        let since = if is_file(&comp) {
            fs::metadata(&comp).and_then(|m| m.modified()).ok()
        } else {
            None
        };
        let mut given = HashMap::new();
        if let Some(Value::Object(st)) = o.get("trio_status") {
            for c in st.keys() {
                let dir = given_dirs.get(c).cloned().flatten();
                given.insert(
                    c.clone(),
                    dir.and_then(|dir| candidates_given(&root, Some(&dir), c, since)),
                );
            }
        }
        if let Some(st) = o.get_mut("trio_status") {
            with_candidates(st, &given);
        }
    }
    if let Some(sum_dirs) = dirs(input, "summaries") {
        let members = o
            .get("state")
            .and_then(|st| get(st, "members").ok().flatten());
        let seq: Vec<String> = match members {
            Some(Value::Array(ms)) => ms
                .iter()
                .filter(|m| get(m, "sequenced").ok().flatten().is_some_and(truthy))
                .filter_map(|m| match get(m, "sample").ok().flatten() {
                    Some(Value::Str(x)) => Some(x.clone()),
                    _ => None,
                })
                .collect(),
            _ => Vec::new(),
        };
        let sums: HashMap<String, Option<Value>> = seq
            .iter()
            .map(|n| {
                let dir = sum_dirs.get(n).cloned().flatten();
                (
                    n.clone(),
                    dir.and_then(|dir| member_summary(&root, Some(&dir))),
                )
            })
            .collect();
        o.insert(
            "side_by_side".into(),
            if seq.is_empty() {
                Value::Null
            } else {
                side_by_side(&seq, &sums)
            },
        );
    }
    Ok(out)
}

/// Each family folder under `root` as /runs.json lists it. `members` is []
/// while the pedigree cannot be read; `unreadable` while family.json cannot
/// be built from the folder.
pub fn families(root: &Path) -> Vec<Value> {
    family_dirs(root)
        .into_iter()
        .map(|(n, d)| {
            let people = members_of(&d).map(|p| people_rows(&p)).unwrap_or_default();
            let readable = !people.is_empty()
                && fs::read_dir(&d).is_ok()
                && family_page(&n, &d).is_ok_and(|p| {
                    let st = get(&p, "state").ok().flatten().unwrap_or(&Value::Null);
                    get(st, "error").ok().flatten().is_none()
                });
            obj([
                ("name", s(&n)),
                ("dir", Value::Str(d.to_string_lossy().into_owned())),
                ("kind", s("family")),
                ("members", Value::Array(people)),
                ("unreadable", Value::Bool(!readable)),
            ])
        })
        .collect()
}

// --- the command line -------------------------------------------------------------------

struct Args {
    pos: Vec<String>,
    opts: HashMap<String, String>,
}

fn parse(args: &[String], known: &[&str]) -> Result<Args, String> {
    let mut a = Args {
        pos: Vec::new(),
        opts: HashMap::new(),
    };
    let mut it = args.iter();
    while let Some(x) = it.next() {
        if let Some(name) = x.strip_prefix("--") {
            if !known.contains(&name) {
                return Err(format!("unknown option --{name}"));
            }
            let v = it.next().ok_or_else(|| format!("--{name} needs a value"))?;
            a.opts.insert(name.to_owned(), v.clone());
        } else {
            a.pos.push(x.clone());
        }
    }
    Ok(a)
}

fn json_line(v: &Value) -> String {
    dumps(v, &Dump::default()).unwrap_or_else(|_| "null".into()) + "\n"
}

fn folder_name(d: &Path) -> String {
    absolute(d)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// `og live family ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let Some(verb) = args.first() else {
        let _ = write!(err, "{USAGE}");
        return 2;
    };
    let known: &[&str] = match verb.as_str() {
        "state" => &["name"],
        "page" => &["name", "root", "input"],
        "families" => &[],
        "svg" => &["title"],
        "files" => &[],
        "-h" | "--help" => {
            let _ = write!(out, "{USAGE}");
            return 0;
        }
        other => {
            let _ = writeln!(err, "{USAGE}og live family: unknown verb '{other}'");
            return 2;
        }
    };
    let a = match parse(&args[1..], known) {
        Ok(a) if a.pos.len() == 1 => a,
        Ok(_) => {
            let _ = writeln!(err, "{USAGE}og live family {verb}: one folder, please");
            return 2;
        }
        Err(e) => {
            let _ = writeln!(err, "{USAGE}og live family {verb}: {e}");
            return 2;
        }
    };
    let d = PathBuf::from(&a.pos[0]);
    let name = a
        .opts
        .get("name")
        .cloned()
        .unwrap_or_else(|| folder_name(&d));
    let answer = match verb.as_str() {
        "state" => family_state(&name, &d).map(|v| json_line(&v)),
        "families" => Ok(json_line(&Value::Array(families(&d)))),
        "files" => match members_of(&d) {
            Ok(p) => Ok(json_line(&Value::Array(
                served_files(&d, &p).iter().map(|f| s(f)).collect(),
            ))),
            Err(e) => {
                let _ = writeln!(err, "ERROR: {e}");
                return 2;
            }
        },
        "svg" => {
            let title = a
                .opts
                .get("title")
                .cloned()
                .unwrap_or_else(|| folder_name(&d));
            match members_of(&d).and_then(|p| layout::draw(&p, &title)) {
                Ok(svg) => Ok(svg),
                Err(e) => {
                    let _ = writeln!(err, "ERROR: {e}");
                    return 2;
                }
            }
        }
        _ => {
            let input = match a.opts.get("input").map(String::as_str) {
                None => Value::Object(IndexMap::new()),
                Some(f) => {
                    let mut text = String::new();
                    let read = if f == "-" {
                        io::stdin().read_to_string(&mut text)
                    } else {
                        fs::File::open(f).and_then(|mut h| h.read_to_string(&mut text))
                    };
                    match read
                        .map_err(|e| e.to_string())
                        .and_then(|_| loads(&text).map_err(|e| e.to_string()))
                    {
                        Ok(v @ Value::Object(_)) => v,
                        Ok(_) | Err(_) => {
                            let _ =
                                writeln!(err, "og live family page: --input is not a JSON object");
                            return 2;
                        }
                    }
                }
            };
            let root = a.opts.get("root").map_or_else(
                || d.parent().map(Path::to_path_buf).unwrap_or_default(),
                PathBuf::from,
            );
            site_page(&name, &d, &root, &input).map(|v| json_line(&v))
        }
    };
    match answer {
        Ok(text) => match out.write_all(text.as_bytes()) {
            Ok(()) => 0,
            Err(e) => {
                let _ = writeln!(err, "og live family: cannot write: {e}");
                1
            }
        },
        Err(Unreadable(kind)) => {
            let _ = writeln!(
                err,
                "og live family: a file of the family cannot be read ({kind})"
            );
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names() {
        assert!(valid_name("family_a") && valid_name("A.b-1"));
        assert!(!valid_name("") && !valid_name(".x") && !valid_name("-x") && !valid_name("a b"));
    }

    #[test]
    fn keys_compare_as_python_does() {
        assert!(
            same_key(&Value::Int(1), &Value::Bool(true))
                && same_key(&Value::Float(1.0), &Value::Int(1))
        );
        assert!(!same_key(&s("1"), &Value::Int(1)) && same_key(&Value::Null, &Value::Null));
    }

    #[test]
    fn walking_a_value() {
        assert_eq!(walk(&s("ab")).unwrap(), vec![s("a"), s("b")]);
        assert!(walk(&Value::Int(3)).is_err() && walk(&Value::Int(0)).unwrap().is_empty());
    }

    #[test]
    fn verdict() {
        assert_eq!(
            verdict_words(&s(" OK ")),
            s("matches the declared relation")
        );
        assert_eq!(verdict_words(&s("mismatch")), s("mismatch"));
        assert_eq!(verdict_words(&Value::Null), Value::Null);
    }
}
