//! What changed since the previous run of the same sample.
//!
//! It compares this run's `summary.json` and `run_manifest.tsv` with the
//! previous run's copy under `history/<stamp>/`:
//!
//! - inputs: versions.env values, each image and the digest it resolved
//!   to, the data rows, the pipeline commit and the declared sex (when and
//!   by which script the manifest was written are left out);
//! - sections: a report section that ran in one run and not the other;
//! - ClinVar hits added, dropped, reclassified, re-starred or with another
//!   genotype; CPSR counts per class; CPIC calls; pypgx's CYP2D6 call and
//!   the CYP2D6 consensus; PRS scores added or dropped, and a percentile
//!   that moved by 5 points or more, or appeared or went away.
//!
//! Each result change names the changed inputs that feed its section, or
//! says that the inputs it tracks are the same. The last database re-check
//! (`recheck/state.json`) is added when there is one.

use crate::files::read_manifest;
use crate::util::{Obj, abspath, basename, normpath, py_eq, s};
use indexmap::IndexMap;
use pycompat::json::{Dump, Value, dumps};
use pycompat::text::html_escape;
use pycompat::value::{py_str, truthy};
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

/// Percentile points a PRS must move to count as a change.
pub const PRS_BAND: f64 = 5.0;

/// The changed inputs that can explain a change in a section: an input
/// whose name contains one of these (case-insensitive). The pipeline
/// commit explains any.
const WHY_TOKENS: [(&str, &[&str]); 6] = [
    (
        "clinvar",
        &[
            "clinvar",
            "DEEPVARIANT",
            "BCFTOOLS",
            "MINIMAP2",
            "BWA",
            "GATK",
            "FASTP",
        ],
    ),
    ("cpsr", &["PCGR", "VEP", "DEEPVARIANT"]),
    (
        "cpic",
        &["PHARMCAT", "PYPGX", "CYRIUS", "DEEPVARIANT", "BCFTOOLS"],
    ),
    ("pypgx", &["PYPGX", "pypgx_bundle"]),
    ("cyp2d6", &["PYPGX", "pypgx_bundle", "CYRIUS", "PHARMCAT"]),
    (
        "prs",
        &["pgs:", "PGSC", "PLINK2", "ancestry_panel", "DEEPVARIANT"],
    ),
];
const IGNORED_RUN_KEYS: [&str; 2] = ["written_utc", "written_by"];

/// The summary section of each result area.
fn section_of_area(area: &str) -> Option<&'static str> {
    Some(match area {
        "ClinVar" => "clinvar",
        "CPSR" => "cpsr",
        "PGx (CPIC)" => "cpic",
        "PGx (pypgx)" | "PGx (CYP2D6 consensus)" => "pypgx",
        "PRS" => "prs",
        _ => return None,
    })
}

/// A summary this comparison cannot read: a value where an object or a
/// list must be.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Malformed(pub String);

impl std::fmt::Display for Malformed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

type R<T> = Result<T, Malformed>;

static NULL: Value = Value::Null;

/// `d.get(key)`, refusing a `d` that is not an object.
fn get<'a>(d: &'a Value, key: &str) -> R<Option<&'a Value>> {
    pycompat::value::get(d, key)
        .map_err(|_| Malformed(format!("not an object where '{key}' was looked up")))
}

/// `d.get(key, default)`.
fn get_or<'a>(d: &'a Value, key: &str, default: &'a Value) -> R<&'a Value> {
    Ok(get(d, key)?.unwrap_or(default))
}

/// `d.get(key)` with `None` as `null`.
fn get_v<'a>(d: &'a Value, key: &str) -> R<&'a Value> {
    get_or(d, key, &NULL)
}

/// The items of a list value.
fn items<'a>(v: &'a Value, what: &str) -> R<&'a [Value]> {
    match v {
        Value::Array(a) => Ok(a),
        _ => Err(Malformed(format!("{what} is not a list"))),
    }
}

/// `str(v)`.
fn st(v: &Value) -> String {
    py_str(v).unwrap_or_default()
}

fn change(
    area: Value,
    kind: &str,
    item: Value,
    before: Value,
    after: Value,
    section: Option<&str>,
) -> Obj {
    let mut c = Obj::new();
    c.insert("area".into(), area);
    c.insert("kind".into(), s(kind));
    c.insert("item".into(), item);
    c.insert("before".into(), before);
    c.insert("after".into(), after);
    if let Some(sec) = section {
        c.insert("section".into(), s(sec));
    }
    c
}

type Manifest = IndexMap<(String, String), String>;

/// The manifest as `{item: value}`, one entry per thing that can change.
/// An image is keyed by its `*_IMAGE` name with `ref (digest)` as its
/// value, so a new tag and a new digest of one tag are each a change.
fn manifest_items(rows: &Manifest) -> IndexMap<String, String> {
    let mut out = IndexMap::new();
    for ((sec, key), val) in rows {
        match sec.as_str() {
            "run" if !IGNORED_RUN_KEYS.contains(&key.as_str()) => {
                out.insert(key.clone(), val.clone());
            }
            "version" if key.ends_with("_IMAGE") => {
                let digest = rows
                    .get(&("image".to_owned(), val.clone()))
                    .map_or("not in the manifest", String::as_str);
                out.insert(key.clone(), format!("{val} ({digest})"));
            }
            "version" | "data" => {
                out.insert(key.clone(), val.clone());
            }
            _ => {}
        }
    }
    out
}

fn diff_inputs(prev: &Manifest, cur: &Manifest) -> Vec<Obj> {
    let (a, b) = (manifest_items(prev), manifest_items(cur));
    let mut keys: Vec<&String> = a.keys().chain(b.keys()).collect();
    keys.sort();
    keys.dedup();
    keys.into_iter()
        .filter(|k| a.get(*k) != b.get(*k))
        .map(|k| {
            let v = |m: &IndexMap<String, String>| s(m.get(k).map_or("absent", String::as_str));
            change(s("inputs"), "input", s(k.as_str()), v(&a), v(&b), None)
        })
        .collect()
}

fn why_for(area: &str, inputs: &[Obj]) -> String {
    let mut tokens: Vec<String> = WHY_TOKENS
        .iter()
        .find(|(a, _)| *a == area)
        .map(|(_, t)| t.iter().map(|x| pycompat::text::lower(x)).collect())
        .unwrap_or_default();
    tokens.push("git_commit".into());
    let hit: Vec<&Obj> = inputs
        .iter()
        .filter(|c| {
            let item = pycompat::text::lower(&st(&c["item"]));
            tokens.iter().any(|t| item.contains(t.as_str()))
        })
        .collect();
    if hit.is_empty() {
        return "the inputs this delta tracks for it are the same; the call itself differs".into();
    }
    let parts: Vec<String> = hit
        .iter()
        .map(|c| {
            format!(
                "{} {} -> {}",
                st(&c["item"]),
                st(&c["before"]),
                st(&c["after"])
            )
        })
        .collect();
    format!("changed inputs: {}", parts.join("; "))
}

/// An empty object, for a missing section.
fn empty_obj() -> &'static Value {
    static EMPTY: LazyLock<Value> = LazyLock::new(|| Value::Object(Obj::new()));
    &EMPTY
}

/// `(summary.get("sections") or {}).get(key) or {}`, borrowed: the summary
/// is never copied.
fn section<'a>(sm: &'a Value, key: &str) -> R<&'a Value> {
    let secs = get_v(sm, "sections")?;
    if !truthy(secs) {
        return Ok(empty_obj());
    }
    let sec = get_v(secs, key)?;
    Ok(if truthy(sec) { sec } else { empty_obj() })
}

fn has_data(sec: &Value) -> R<bool> {
    Ok(matches!(get(sec, "state")?, Some(Value::Str(x)) if x == "ok" || x == "stale"))
}

fn section_keys(sm: &Value) -> R<Vec<String>> {
    match get_v(sm, "sections")? {
        Value::Object(o) => Ok(o.keys().cloned().collect()),
        v if !truthy(v) => Ok(Vec::new()),
        _ => Err(Malformed("sections is not an object".into())),
    }
}

fn diff_states(p: &Value, c: &Value) -> R<Vec<Obj>> {
    let mut keys = section_keys(p)?;
    keys.extend(section_keys(c)?);
    keys.sort();
    keys.dedup();
    let missing = s("missing");
    let mut out = Vec::new();
    for k in keys {
        let (a, b) = (section(p, &k)?, section(c, &k)?);
        // A section only one summary knows (a newer version added it) and
        // that did not run reads as missing on both sides.
        let sa = get_or(a, "state", &missing)?.clone();
        let sb = get_or(b, "state", &missing)?.clone();
        if !py_eq(&sa, &sb) {
            let title = [get_v(b, "title")?, get_v(a, "title")?]
                .into_iter()
                .find(|v| truthy(v))
                .cloned()
                .unwrap_or_else(|| s(k.as_str()));
            out.push(change(title, "section", s(k.as_str()), sa, sb, Some(&k)));
        }
    }
    Ok(out)
}

/// A list of objects looked up by a key value, the later of two equal keys
/// winning in the earlier's place, as a dict comprehension does. String
/// keys, the usual case, are hashed; Python never finds a `str` equal to
/// anything but a `str`, so any other key is compared with [`py_eq`] against
/// the other non-string keys only. That keeps a summary of tens of
/// thousands of ClinVar hits linear.
struct Index<'a> {
    entries: Vec<(Value, &'a Value)>,
    by_str: HashMap<String, usize>,
}

impl<'a> Index<'a> {
    fn new() -> Self {
        Index {
            entries: Vec::new(),
            by_str: HashMap::new(),
        }
    }

    fn position(&self, k: &Value) -> Option<usize> {
        match k {
            Value::Str(t) => self.by_str.get(t).copied(),
            _ => self
                .entries
                .iter()
                .position(|(x, _)| !matches!(x, Value::Str(_)) && py_eq(x, k)),
        }
    }

    /// Adds `k`, or replaces the value of an equal key in its place.
    fn insert(&mut self, k: Value, v: &'a Value) {
        match self.position(&k) {
            Some(i) => self.entries[i].1 = v,
            None => {
                if let Value::Str(t) = &k {
                    self.by_str.insert(t.clone(), self.entries.len());
                }
                self.entries.push((k, v));
            }
        }
    }

    fn get(&self, k: &Value) -> Option<&'a Value> {
        self.position(k).map(|i| self.entries[i].1)
    }
}

fn index_by<'a>(list: &'a [Value], key: impl Fn(&Value) -> R<Value>) -> R<Index<'a>> {
    let mut out = Index::new();
    for v in list {
        out.insert(key(v)?, v);
    }
    Ok(out)
}

/// The union of two indexes' keys, sorted by `str()` (a stable sort).
fn union_keys(a: &Index<'_>, b: &Index<'_>) -> Vec<Value> {
    let mut seen = Index::new();
    let none = Value::Null;
    for (k, _) in a.entries.iter().chain(&b.entries) {
        if seen.position(k).is_none() {
            seen.insert(k.clone(), &none);
        }
    }
    let mut keys: Vec<Value> = seen.entries.into_iter().map(|(k, _)| k).collect();
    keys.sort_by_cached_key(st);
    keys
}

fn lookup<'a>(idx: &Index<'a>, k: &Value) -> Option<&'a Value> {
    idx.get(k)
}

fn data_list<'a>(
    sec: &'a Value,
    key: &str,
    empty_obj: &'a Value,
    empty_list: &'a Value,
) -> R<&'a [Value]> {
    let data = get_or(sec, "data", empty_obj)?;
    items(get_or(data, key, empty_list)?, key)
}

fn diff_clinvar(p: &Value, c: &Value) -> R<Vec<Obj>> {
    let (a, b) = (section(p, "clinvar")?, section(c, "clinvar")?);
    if !(has_data(a)? && has_data(b)?) {
        return Ok(Vec::new());
    }
    let (eo, el) = (Value::Object(Obj::new()), Value::Array(Vec::new()));
    let key = |h: &Value| -> R<Value> {
        Ok(s(format!(
            "{}:{} {}>{}",
            st(get_v(h, "chrom")?),
            st(get_v(h, "pos")?),
            st(get_v(h, "ref")?),
            st(get_v(h, "alt")?)
        )))
    };
    let ha = index_by(data_list(a, "hits", &eo, &el)?, key)?;
    let hb = index_by(data_list(b, "hits", &eo, &el)?, key)?;
    let mut out = Vec::new();
    let dot = s(".");
    for k in union_keys(&ha, &hb) {
        let (x, y) = (lookup(&ha, &k), lookup(&hb, &k));
        let either = [y, x]
            .into_iter()
            .flatten()
            .find(|v| truthy(v))
            .or(x)
            .or(y)
            .expect("in one index");
        let item = s(format!(
            "{} ({})",
            st(&k),
            st(get_or(either, "gene", &dot)?)
        ));
        let called = |h: &Value| -> R<Value> {
            Ok(s(format!(
                "{}, {} stars",
                st(get_v(h, "significance")?),
                st(get_v(h, "stars")?)
            )))
        };
        match (x, y) {
            (None, Some(y)) => out.push(change(
                s("ClinVar"),
                "added",
                item,
                s("absent"),
                called(y)?,
                None,
            )),
            (Some(x), None) => out.push(change(
                s("ClinVar"),
                "dropped",
                item,
                called(x)?,
                s("absent"),
                None,
            )),
            (Some(x), Some(y)) => {
                let (sx, sy) = (get_v(x, "significance")?, get_v(y, "significance")?);
                if !py_eq(sx, sy) {
                    out.push(change(
                        s("ClinVar"),
                        "reclassified",
                        item.clone(),
                        sx.clone(),
                        sy.clone(),
                        None,
                    ));
                }
                if !py_eq(get_v(x, "stars")?, get_v(y, "stars")?)
                    || !py_eq(get_v(x, "review_status")?, get_v(y, "review_status")?)
                {
                    let shown = |h: &Value| -> R<Value> {
                        Ok(s(format!(
                            "{} stars ({})",
                            st(get_v(h, "stars")?),
                            st(get_v(h, "review_status")?)
                        )))
                    };
                    out.push(change(
                        s("ClinVar"),
                        "re-starred",
                        item.clone(),
                        shown(x)?,
                        shown(y)?,
                        None,
                    ));
                }
                let (gx, gy) = (get_v(x, "genotype")?, get_v(y, "genotype")?);
                if !py_eq(gx, gy) {
                    out.push(change(
                        s("ClinVar"),
                        "genotype",
                        item,
                        gx.clone(),
                        gy.clone(),
                        None,
                    ));
                }
            }
            (None, None) => {}
        }
    }
    Ok(out)
}

fn diff_cpsr(p: &Value, c: &Value) -> R<Vec<Obj>> {
    let (a, b) = (section(p, "cpsr")?, section(c, "cpsr")?);
    if !(has_data(a)? && has_data(b)?) {
        return Ok(Vec::new());
    }
    let eo = Value::Object(Obj::new());
    let classes = |sec: &Value| -> R<Obj> {
        match get_v(get_or(sec, "data", &eo)?, "classification")? {
            Value::Object(o) => Ok(o.clone()),
            v if !truthy(v) => Ok(Obj::new()),
            _ => Err(Malformed("classification is not an object".into())),
        }
    };
    let (ca, cb) = (classes(a)?, classes(b)?);
    let mut keys: Vec<&String> = ca.keys().chain(cb.keys()).collect();
    keys.sort();
    keys.dedup();
    let zero = Value::Int(0);
    Ok(keys
        .into_iter()
        .filter_map(|k| {
            let (x, y) = (ca.get(k).unwrap_or(&zero), cb.get(k).unwrap_or(&zero));
            (!py_eq(x, y)).then(|| {
                change(
                    s("CPSR"),
                    "classification count",
                    s(k.as_str()),
                    x.clone(),
                    y.clone(),
                    None,
                )
            })
        })
        .collect())
}

fn diff_cpic(p: &Value, c: &Value) -> R<Vec<Obj>> {
    let (a, b) = (section(p, "cpic")?, section(c, "cpic")?);
    if !(has_data(a)? && has_data(b)?) {
        return Ok(Vec::new());
    }
    let (eo, el) = (Value::Object(Obj::new()), Value::Array(Vec::new()));
    let ga = index_by(data_list(a, "genes", &eo, &el)?, |g| {
        Ok(get_v(g, "gene")?.clone())
    })?;
    let gb = index_by(data_list(b, "genes", &eo, &el)?, |g| {
        Ok(get_v(g, "gene")?.clone())
    })?;
    let call = |g: Option<&Value>| -> R<String> {
        Ok(match g.filter(|g| truthy(g)) {
            Some(g) => format!(
                "{} ({})",
                st(get_v(g, "diplotype")?),
                st(get_v(g, "phenotype")?)
            ),
            None => "absent".into(),
        })
    };
    let mut out = Vec::new();
    for k in union_keys(&ga, &gb) {
        let (fx, fy) = (call(lookup(&ga, &k))?, call(lookup(&gb, &k))?);
        if fx != fy {
            out.push(change(s("PGx (CPIC)"), "call", k, s(fx), s(fy), None));
        }
    }
    Ok(out)
}

fn diff_cyp2d6(p: &Value, c: &Value) -> R<Vec<Obj>> {
    let mut out = Vec::new();
    let (a, b) = (section(p, "pypgx")?, section(c, "pypgx")?);
    let eo = Value::Object(Obj::new());
    if has_data(a)? && has_data(b)? {
        let x = get_v(get_or(a, "data", &eo)?, "cyp2d6")?;
        let y = get_v(get_or(b, "data", &eo)?, "cyp2d6")?;
        if !py_eq(x, y) {
            out.push(change(
                s("PGx (pypgx)"),
                "call",
                s("CYP2D6"),
                x.clone(),
                y.clone(),
                None,
            ));
        }
    }
    let consensus = |sm: &Value| -> R<Value> {
        let c = get_v(sm, "cyp2d6")?;
        let c = if truthy(c) { c } else { &eo };
        let k = get_v(c, "consensus")?;
        let k = if truthy(k) { k } else { &eo };
        Ok(get_v(k, "result")?.clone())
    };
    let (x, y) = (consensus(p)?, consensus(c)?);
    if !py_eq(&x, &y) {
        let shown = |v: Value| if truthy(&v) { v } else { s("none") };
        out.push(change(
            s("PGx (CYP2D6 consensus)"),
            "call",
            s("CYP2D6"),
            shown(x),
            shown(y),
            None,
        ));
    }
    Ok(out)
}

/// `float(v)`, or `None` when it is not a number.
fn pct(v: &Value) -> Option<f64> {
    match v {
        Value::Int(i) => Some(*i as f64),
        Value::Float(x) => Some(*x),
        Value::Bool(b) => Some(f64::from(u8::from(*b))),
        Value::Str(t) => pycompat::num::parse_float(t).ok(),
        _ => None,
    }
}

fn diff_prs(p: &Value, c: &Value) -> R<Vec<Obj>> {
    let (a, b) = (section(p, "prs")?, section(c, "prs")?);
    if !(has_data(a)? && has_data(b)?) {
        return Ok(Vec::new());
    }
    let (eo, el) = (Value::Object(Obj::new()), Value::Array(Vec::new()));
    let sa = index_by(data_list(a, "scores", &eo, &el)?, |r| {
        Ok(get_v(r, "pgs_id")?.clone())
    })?;
    let sb = index_by(data_list(b, "scores", &eo, &el)?, |r| {
        Ok(get_v(r, "pgs_id")?.clone())
    })?;
    let empty = s("");
    let mut out = Vec::new();
    for k in union_keys(&sa, &sb) {
        let (x, y) = (lookup(&sa, &k), lookup(&sb, &k));
        let either = [y, x]
            .into_iter()
            .flatten()
            .find(|v| truthy(v))
            .or(x)
            .or(y)
            .expect("in one index");
        let item = s(format!(
            "{} ({})",
            st(&k),
            st(get_or(either, "condition", &empty)?)
        ));
        let (Some(x), Some(y)) = (x, y) else {
            let added = x.is_none();
            out.push(change(
                s("PRS"),
                if added { "added" } else { "dropped" },
                item,
                s(if added { "absent" } else { "scored" }),
                s(if added { "scored" } else { "absent" }),
                None,
            ));
            continue;
        };
        let (vx, vy) = (get_v(x, "percentile")?, get_v(y, "percentile")?);
        let (px, py) = (pct(vx), pct(vy));
        let moved = match (px, py) {
            (Some(a), Some(b)) => (a - b).abs() >= PRS_BAND,
            (None, None) => false,
            _ => true,
        };
        if moved {
            // A percentile of 0 is a value: only a missing one reads "none".
            let shown = |v: &Value| {
                if *v == Value::Null {
                    s("none")
                } else {
                    v.clone()
                }
            };
            out.push(change(
                s("PRS"),
                "percentile",
                item,
                shown(vx),
                shown(vy),
                None,
            ));
        }
    }
    Ok(out)
}

/// The changes between two runs: inputs, then sections, then results.
pub fn compare(
    prev: &Value,
    prev_manifest: &Manifest,
    cur: &Value,
    manifest: &Manifest,
) -> R<Vec<Obj>> {
    let inputs = diff_inputs(prev_manifest, manifest);
    type Diff = fn(&Value, &Value) -> R<Vec<Obj>>;
    let areas: [(&str, Diff); 5] = [
        ("clinvar", diff_clinvar),
        ("cpsr", diff_cpsr),
        ("cpic", diff_cpic),
        ("cyp2d6", diff_cyp2d6),
        ("prs", diff_prs),
    ];
    let mut results = Vec::new();
    for (area, f) in areas {
        for mut ch in f(prev, cur)? {
            ch.insert("why".into(), s(why_for(area, &inputs)));
            let sec = section_of_area(&st(&ch["area"])).unwrap_or(area);
            ch.insert("section".into(), s(sec));
            results.push(ch);
        }
    }
    let mut out = inputs;
    out.extend(diff_states(prev, cur)?);
    out.extend(results);
    Ok(out)
}

/// The summary sections a delta changes, in order: what the report marks
/// "New since".
pub fn changed_sections(changes: &[Obj]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for c in changes {
        if let Some(sec) = c.get("section").filter(|v| truthy(v)) {
            let sec = st(sec);
            if !out.contains(&sec) {
                out.push(sec);
            }
        }
    }
    out
}

/// The previous run of a sample: the newest `history/<stamp>/` folder
/// (by name, which sorts in run order) that holds a `summary.json`. A
/// half-copied folder without one is passed over. The report step and the
/// summary's "New since" both pick the previous run by this one rule.
pub fn newest_history(sample_dir: &Path) -> Option<PathBuf> {
    let h = sample_dir.join("history");
    let mut names: Vec<std::ffi::OsString> = std::fs::read_dir(&h)
        .ok()?
        .filter_map(|e| e.ok().map(|e| e.file_name()))
        .filter(|n| h.join(n).join("summary.json").is_file())
        .collect();
    pycompat::mtime::sort(&mut names);
    names.pop().map(|n| h.join(n))
}

/// `'2026-09-02 12:00:00 UTC'` as `'20260902120000'`: its first 14 digits.
fn utc14(v: Option<&Value>) -> String {
    let text = v.filter(|v| truthy(v)).map(st).unwrap_or_default();
    text.chars()
        .filter(|c| pycompat::text::is_digit(*c))
        .take(14)
        .collect()
}

/// `<sample>/recheck/state.json` as a short block; `None` without one, or
/// when it checked an earlier run (the sample's `run_manifest.tsv` was
/// written after it, so its "newer than the run used" may be what this run
/// used).
pub fn read_recheck(path: &Path) -> Option<Obj> {
    let st_ = pycompat::json::loads(&crate::util::read_text(path).ok()?).ok()?;
    let Value::Object(o) = &st_ else { return None };
    let Some(Value::Array(list)) = o.get("items") else {
        return None;
    };
    let dir = abspath(path);
    let man_path = dir
        .parent()
        .and_then(Path::parent)
        .map(|d| d.join("run_manifest.tsv"));
    let man = read_manifest(man_path.as_deref()).unwrap_or_default();
    let written = man
        .get(&("run".to_owned(), "written_utc".to_owned()))
        .map(|v| s(v.as_str()));
    if utc14(written.as_ref()) > utc14(o.get("checked_utc")) {
        return None;
    }
    let mut out = Obj::new();
    out.insert(
        "checked_utc".into(),
        o.get("checked_utc").cloned().unwrap_or(Value::Null),
    );
    let items = list
        .iter()
        .filter_map(|i| match i {
            Value::Object(i) => Some(Value::Object(
                ["item", "label", "class", "ran_with", "upstream", "action"]
                    .iter()
                    .map(|k| ((*k).to_owned(), i.get(*k).cloned().unwrap_or(Value::Null)))
                    .collect(),
            )),
            _ => None,
        })
        .collect();
    out.insert("items".into(), Value::Array(items));
    out.insert(
        "affected_steps".into(),
        o.get("affected_steps")
            .filter(|v| truthy(v))
            .cloned()
            .unwrap_or(Value::Array(Vec::new())),
    );
    out.insert(
        "rerun".into(),
        match o.get("rerun") {
            Some(r @ Value::Object(_)) => r.clone(),
            _ => Value::Null,
        },
    );
    Some(out)
}

/// The re-check lines of the text delta; none without a check.
fn recheck_text(r: Option<&Value>) -> Vec<String> {
    let Some(r) = r.filter(|r| truthy(r)) else {
        return Vec::new();
    };
    let f = |v: &Value, k: &str| match v {
        Value::Object(o) => o.get(k).cloned().unwrap_or(Value::Null),
        _ => Value::Null,
    };
    let items = match f(r, "items") {
        Value::Array(a) => a,
        _ => Vec::new(),
    };
    let newer: Vec<&Value> = items
        .iter()
        .filter(|i| {
            let c = f(i, "class");
            !(c == Value::Null || py_eq(&c, &s("no-newer-release")))
        })
        .collect();
    let checked = f(r, "checked_utc");
    let mut out = vec![format!(
        "Database re-check of {}: {}",
        if truthy(&checked) {
            st(&checked)
        } else {
            "unknown date".into()
        },
        if newer.is_empty() {
            "nothing newer than this run".to_owned()
        } else {
            format!("{} item(s) to look at", newer.len())
        }
    )];
    for i in newer {
        let label = f(i, "label");
        let name = if truthy(&label) { label } else { f(i, "item") };
        let action = f(i, "action");
        out.push(format!(
            "  {}: {}{}",
            st(&name),
            st(&f(i, "class")),
            if truthy(&action) {
                format!(" ({})", st(&action))
            } else {
                String::new()
            }
        ));
    }
    let rr = f(r, "rerun");
    if truthy(&rr) {
        let fin = f(&rr, "finished_utc");
        out.push(format!(
            "  last rerun: {} (exit {})",
            if truthy(&fin) {
                st(&fin)
            } else {
                "unfinished".into()
            },
            st(&f(&rr, "rc"))
        ));
    }
    out
}

/// `v or default`, as text.
fn or_text(v: &Value, default: &str) -> String {
    if truthy(v) { st(v) } else { default.to_owned() }
}

fn field<'a>(d: &'a Obj, k: &str) -> &'a Value {
    d.get(k).unwrap_or(&NULL)
}

/// The text delta.
pub fn render_text(d: &Obj) -> String {
    let sample = st(field(d, "sample"));
    let mut w = vec![
        format!("What changed since the previous run: {sample}"),
        String::new(),
    ];
    let rc = recheck_text(d.get("recheck"));
    if *field(d, "previous") == Value::Null {
        w.push("No previous run of this sample on record (its history/ folder is empty).".into());
        w.push("This run is kept as the baseline the next run is compared with.".into());
        w.extend(rc);
        return w.join("\n") + "\n";
    }
    w.push(format!(
        "Previous run: history/{} (summary of {})",
        st(field(d, "previous")),
        or_text(field(d, "previous_generated"), "unknown date")
    ));
    w.push(format!(
        "This run:     summary of {}",
        or_text(field(d, "current_generated"), "unknown date")
    ));
    w.extend(rc);
    w.push(String::new());
    let changes = match field(d, "changes") {
        Value::Array(a) => a.as_slice(),
        _ => &[],
    };
    if changes.is_empty() {
        w.push("No change: the same inputs and the same results as the previous run.".into());
        return w.join("\n") + "\n";
    }
    w.push(format!(
        "New since {}: {} change(s):",
        or_text(field(d, "previous_utc"), "the previous run"),
        changes.len()
    ));
    let mut area: Option<&Value> = None;
    for ch in changes {
        let Value::Object(c) = ch else { continue };
        let a = field(c, "area");
        if area.is_none_or(|x| !py_eq(x, a)) {
            area = Some(a);
            w.push(String::new());
            w.push(format!("## {}", st(a)));
        }
        w.push(format!(
            "  {}: {}: {} -> {}",
            st(field(c, "kind")),
            st(field(c, "item")),
            st(field(c, "before")),
            st(field(c, "after"))
        ));
        if truthy(field(c, "why")) {
            w.push(format!("    why: {}", st(field(c, "why"))));
        }
    }
    w.join("\n") + "\n"
}

/// The HTML delta (dark page, one table).
pub fn render_html(d: &Obj) -> String {
    let e = |t: &str| html_escape(t).into_owned();
    let sample = st(field(d, "sample"));
    let mut body = vec![format!(
        "<h1>What changed since the previous run: {}</h1>",
        e(&sample)
    )];
    let has_prev = *field(d, "previous") != Value::Null;
    if !has_prev {
        body.push(
            "<p>No previous run of this sample on record (its history/ folder is empty). This run is kept as the baseline the next run is compared with.</p>"
                .into(),
        );
    } else {
        body.push(format!(
            "<p>Previous run: history/{} (summary of {}). This run: summary of {}.</p>",
            e(&st(field(d, "previous"))),
            e(&or_text(field(d, "previous_generated"), "unknown date")),
            e(&or_text(field(d, "current_generated"), "unknown date"))
        ));
    }
    let rc = recheck_text(d.get("recheck"));
    if !rc.is_empty() {
        let lines: Vec<String> = rc
            .iter()
            .map(|x| e(pycompat::text::strip(x, None)))
            .collect();
        body.push(format!("<p>{}</p>", lines.join("<br>")));
    }
    if has_prev {
        let changes = match field(d, "changes") {
            Value::Array(a) => a.as_slice(),
            _ => &[],
        };
        if changes.is_empty() {
            body.push(
                "<p>No change: the same inputs and the same results as the previous run.</p>"
                    .into(),
            );
        } else {
            body.push(format!(
                "<p>New since {}: {} change(s).</p>",
                e(&or_text(field(d, "previous_utc"), "the previous run")),
                changes.len()
            ));
            body.push(
                "<table><tr><th>Area</th><th>Change</th><th>Item</th><th>Before</th><th>Now</th><th>Why</th></tr>".into(),
            );
            for ch in changes {
                let Value::Object(c) = ch else { continue };
                let cells: String = ["area", "kind", "item", "before", "after", "why"]
                    .iter()
                    .map(|k| format!("<td>{}</td>", e(&c.get(*k).map_or_else(String::new, st))))
                    .collect();
                body.push(format!("<tr>{cells}</tr>"));
            }
            body.push("</table>".into());
        }
    }
    let style = "body{font-family:sans-serif;margin:2em;background:#111;color:#ddd}table{border-collapse:collapse}td,th{border:1px solid #444;padding:4px 8px;text-align:left;vertical-align:top}";
    format!(
        "<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><title>{}: delta</title><style>{style}</style></head><body>\n{}\n</body></html>\n",
        e(&sample),
        body.join("\n")
    )
}

/// What `og run delta` was given.
pub struct DeltaArgs<'a> {
    pub sample: &'a str,
    pub summary: &'a Path,
    pub manifest: Option<&'a Path>,
    pub previous: Option<&'a Path>,
    pub out: &'a str,
    pub recheck: Option<&'a Path>,
}

/// Why the delta was not written.
#[derive(Debug)]
pub enum DeltaError {
    /// A summary that cannot be read or parsed (exit 2).
    Summary(String),
    /// An output that cannot be written.
    Write(String, std::io::Error),
}

fn load_json(p: &Path) -> Result<Value, String> {
    let text = crate::util::read_text(p).map_err(|e| format!("{}: {e}", p.display()))?;
    pycompat::json::loads(&text).map_err(|e| format!("{}: {e}", p.display()))
}

/// Writes `OUT.txt`, `OUT.html` and `OUT.json`, and returns the text.
pub fn run(a: &DeltaArgs<'_>) -> Result<String, DeltaError> {
    let bad = |e: String| DeltaError::Summary(e);
    let cur = load_json(a.summary).map_err(bad)?;
    let prev = match a.previous {
        Some(p) => Some(load_json(&p.join("summary.json")).map_err(bad)?),
        None => None,
    };
    let malformed = |e: Malformed| DeltaError::Summary(e.0);
    let mut d = Obj::new();
    d.insert("sample".into(), s(a.sample));
    d.insert(
        "previous".into(),
        a.previous
            .map_or(Value::Null, |p| s(basename(&normpath(p)))),
    );
    let gen_of = |v: &Value| get_v(v, "generated_utc").cloned();
    d.insert(
        "previous_generated".into(),
        prev.as_ref()
            .map_or(Ok(Value::Null), gen_of)
            .map_err(malformed)?,
    );
    d.insert("current_generated".into(), gen_of(&cur).map_err(malformed)?);
    d.insert("prs_band".into(), Value::Float(PRS_BAND));
    d.insert("changes".into(), Value::Array(Vec::new()));
    if let (Some(prev), Some(pdir)) = (&prev, a.previous) {
        d.insert("previous_utc".into(), gen_of(prev).map_err(malformed)?);
        let pm =
            read_manifest(Some(&pdir.join("run_manifest.tsv"))).map_err(|e| bad(e.to_string()))?;
        let cm = read_manifest(a.manifest).map_err(|e| bad(e.to_string()))?;
        let changes = compare(prev, &pm, &cur, &cm).map_err(malformed)?;
        let secs = changed_sections(&changes);
        d.insert(
            "changes".into(),
            Value::Array(changes.into_iter().map(Value::Object).collect()),
        );
        d.insert(
            "new_since_sections".into(),
            Value::Array(secs.into_iter().map(Value::Str).collect()),
        );
    }
    let recheck_path = a.recheck.map_or_else(
        || {
            abspath(a.summary)
                .parent()
                .unwrap_or(Path::new("/"))
                .join("recheck")
                .join("state.json")
        },
        Path::to_path_buf,
    );
    if let Some(r) = read_recheck(&recheck_path) {
        d.insert("recheck".into(), Value::Object(r));
    }
    let text = render_text(&d);
    let html = render_html(&d);
    let json = dumps(
        &Value::Object(d),
        &Dump {
            indent: Some(1),
            ..Dump::default()
        },
    )
    .map_err(|e| bad(e.to_string()))?
        + "\n";
    for (ext, body) in [("txt", text.clone()), ("html", html), ("json", json)] {
        let p = format!("{}.{ext}", a.out);
        let write = || -> std::io::Result<()> {
            let mut f = og_io::AtomicFile::create(&p)?;
            f.write_all(body.as_bytes())?;
            f.commit()
        };
        write().map_err(|e| DeltaError::Write(p.clone(), e))?;
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pycompat::json::loads;

    fn summary(prs: &str) -> Value {
        loads(&format!(
            r#"{{"generated_utc":"2026-01-01","sections":{{"prs":{{"state":"ok","title":"PRS","data":{{"scores":[{prs}]}}}}}}}}"#
        ))
        .unwrap()
    }

    #[test]
    fn a_zero_percentile_is_a_value() {
        let a = summary(r#"{"pgs_id":"PGS1","condition":"x","percentile":0}"#);
        let b = summary(r#"{"pgs_id":"PGS1","condition":"x","percentile":7.5}"#);
        let ch = compare(&a, &Manifest::new(), &b, &Manifest::new()).unwrap();
        assert_eq!(ch.len(), 1);
        assert_eq!(ch[0]["before"], Value::Int(0));
        assert_eq!(ch[0]["after"], Value::Float(7.5));
        let none = summary(r#"{"pgs_id":"PGS1","condition":"x"}"#);
        let ch = compare(&none, &Manifest::new(), &b, &Manifest::new()).unwrap();
        assert_eq!(ch[0]["before"], s("none"));
    }

    #[test]
    fn the_newest_complete_history_folder() {
        let d = std::env::temp_dir().join(format!("og-delta-hist-{}", std::process::id()));
        let h = d.join("history");
        for (name, full) in [
            ("2026-01-01T000000Z", true),
            ("2026-02-01T000000Z", true),
            ("2026-03-01T000000Z", false),
        ] {
            std::fs::create_dir_all(h.join(name)).unwrap();
            if full {
                std::fs::write(h.join(name).join("summary.json"), "{}").unwrap();
            }
        }
        assert_eq!(newest_history(&d), Some(h.join("2026-02-01T000000Z")));
        assert_eq!(newest_history(&d.join("none")), None);
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn malformed_summaries_are_refused() {
        let a = loads(r#"{"sections":{"clinvar":{"state":"ok","data":{"hits":[1]}}}}"#).unwrap();
        assert!(compare(&a, &Manifest::new(), &a, &Manifest::new()).is_err());
    }
}
