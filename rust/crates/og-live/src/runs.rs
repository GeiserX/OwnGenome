//! The runs the live site knows, and the files of a run folder a browser
//! may fetch.
//!
//! A run is a folder under the root that holds `logs/run_status.tsv`,
//! `nextflow/` or `run_manifest.tsv` (a person's run), a run known only
//! from its stored messages, or a family folder (`og live family`). Who is a
//! child ([`child_of`]) is computed from the families' files on each
//! request, failing closed while a family cannot be read.
//!
//! A run folder serves only what its report links to by a relative href,
//! the files a plugin names ([`crate::plugins::RUN_GLOBS`]), and the one
//! HTML file of a folder the report links (CPSR's report in `cpsr/`).
//! Nothing that resolves out of the run folder or the root is served.

use crate::family::{inside, realpath};
use crate::static_files::{normpath, unescape};
use indexmap::IndexMap;
use pycompat::json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

/// A run name, as `validate_sample` in `scripts/lib/common.sh` takes it: a
/// letter or digit, then letters, digits, `.`, `_` and `-`.
pub fn valid_run(name: &str) -> bool {
    crate::family::valid_name(name)
}

/// One run the site knows.
#[derive(Debug, Clone, PartialEq)]
pub struct Run {
    pub name: String,
    /// Its folder under the root; `None` for a run known from its messages
    /// alone.
    pub dir: Option<PathBuf>,
    /// `sample` or `family`.
    pub kind: String,
    /// A family's members (`[{sample, father, mother, ...}]`), as
    /// `og live family families` lists them.
    pub members: Vec<Value>,
    /// A family whose page cannot be built from its folder.
    pub unreadable: bool,
}

impl Run {
    /// The folder the run's files are read from: its own, else the one it
    /// would have under the root.
    pub fn folder(&self, root: &Path) -> PathBuf {
        self.dir.clone().unwrap_or_else(|| root.join(&self.name))
    }

    pub fn is_family(&self) -> bool {
        self.kind == "family"
    }

    /// The run as `/runs.json` and `child_of` read a family:
    /// `{name, members, unreadable}`.
    pub fn family_value(&self) -> Value {
        let mut o = IndexMap::new();
        o.insert("name".into(), Value::Str(self.name.clone()));
        o.insert("members".into(), Value::Array(self.members.clone()));
        o.insert("unreadable".into(), Value::Bool(self.unreadable));
        Value::Object(o)
    }
}

fn field<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(key),
        _ => None,
    }
}

fn truthy(v: Option<&Value>) -> bool {
    v.is_some_and(pycompat::value::truthy)
}

/// Whether run `name`'s pages treat it as a child: a run that any readable
/// family lists with a father or a mother is a child, whatever the other
/// families say. While a family cannot be read, a run it might list is
/// treated as a child too (fail closed): every run while its members cannot
/// be read (an empty member list), and the members it lists while another
/// of its files cannot be read (`unreadable`).
pub fn child_of(families: &[Value], name: &str) -> bool {
    let listed = |f: &Value| -> Vec<Value> {
        match field(f, "members") {
            Some(Value::Array(ms)) => ms
                .iter()
                .filter(|m| matches!(field(m, "sample"), Some(Value::Str(s)) if s == name))
                .cloned()
                .collect(),
            _ => Vec::new(),
        }
    };
    let members_empty = |f: &Value| !truthy(field(f, "members"));
    if families.iter().any(|f| {
        listed(f)
            .iter()
            .any(|m| truthy(field(m, "father")) || truthy(field(m, "mother")))
    }) {
        return true;
    }
    families
        .iter()
        .any(|f| members_empty(f) || (truthy(field(f, "unreadable")) && !listed(f).is_empty()))
}

/// The runs under `root`, by name: the person runs' folders, then the runs
/// known from `known` (their stored messages), then the family folders.
pub fn list(root: &Path, known: &[String]) -> BTreeMap<String, Run> {
    let mut out = BTreeMap::new();
    let mut names: Vec<String> = fs::read_dir(root)
        .map(|rd| {
            rd.filter_map(|e| e.ok()?.file_name().into_string().ok())
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    for n in names {
        let d = root.join(&n);
        if n.starts_with('.') || !valid_run(&n) || !inside(&d, root) || !d.is_dir() {
            continue;
        }
        if ["logs/run_status.tsv", "nextflow", "run_manifest.tsv"]
            .iter()
            .any(|x| d.join(x).exists())
        {
            out.insert(
                n.clone(),
                Run {
                    name: n,
                    dir: Some(d),
                    kind: "sample".into(),
                    members: Vec::new(),
                    unreadable: false,
                },
            );
        }
    }
    for n in known {
        out.entry(n.clone()).or_insert_with(|| Run {
            name: n.clone(),
            dir: None,
            kind: "sample".into(),
            members: Vec::new(),
            unreadable: false,
        });
    }
    for f in crate::family::families(root) {
        let Some(Value::Str(name)) = field(&f, "name") else {
            continue;
        };
        if !valid_run(name) || out.contains_key(name) {
            continue;
        }
        let dir = match field(&f, "dir") {
            Some(Value::Str(d)) => Some(PathBuf::from(d)),
            _ => None,
        };
        let members = match field(&f, "members") {
            Some(Value::Array(m)) => m.clone(),
            _ => Vec::new(),
        };
        out.insert(
            name.clone(),
            Run {
                name: name.clone(),
                dir,
                kind: match field(&f, "kind") {
                    Some(Value::Str(k)) => k.clone(),
                    _ => "other".into(),
                },
                members,
                unreadable: truthy(field(&f, "unreadable")),
            },
        );
    }
    out
}

/// The families among `runs`, as [`child_of`] reads them.
pub fn families_of(runs: &BTreeMap<String, Run>) -> Vec<Value> {
    runs.values()
        .filter(|r| r.is_family())
        .map(Run::family_value)
        .collect()
}

/// The relative links of `<run>_report.html`: every `href="..."` up to a
/// `#` or `?`, unescaped and unquoted, that names no scheme, is not
/// absolute and does not leave the run folder.
pub fn report_links(html: &str) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    let mut rest = html;
    while let Some(i) = rest.find("href=\"") {
        rest = &rest[i + 6..];
        let end = rest.find(['"', '#', '?']).unwrap_or(rest.len());
        if end == 0 {
            continue;
        }
        let h = unescape(&rest[..end]);
        rest = &rest[end..];
        if has_scheme(&h) || h.starts_with('/') {
            continue;
        }
        let rel = normpath(&crate::viewer::unquote(&h));
        if !rel.starts_with("..") && rel != "." {
            out.insert(rel);
        }
    }
    out
}

/// `^[a-z][a-z0-9+.-]*:`, any case.
fn has_scheme(h: &str) -> bool {
    let b = h.as_bytes();
    if !b.first().is_some_and(u8::is_ascii_alphabetic) {
        return false;
    }
    for &c in &b[1..] {
        if c == b':' {
            return true;
        }
        if !(c.is_ascii_alphanumeric() || matches!(c, b'+' | b'.' | b'-')) {
            return false;
        }
    }
    false
}

/// The relative paths of run `name`'s folder `d` that its report links to.
pub fn run_files(name: &str, d: Option<&Path>) -> BTreeSet<String> {
    let Some(d) = d else {
        return BTreeSet::new();
    };
    let rep = d.join(format!("{name}_report.html"));
    match fs::read(&rep) {
        Ok(b) if rep.is_file() => report_links(&String::from_utf8_lossy(&b)),
        _ => BTreeSet::new(),
    }
}

/// Whether the browser may fetch `rel` of the run: a report link, a file a
/// plugin names, or the one HTML file of a folder the report links.
pub fn file_allowed(run: &Run, root: &Path, rel: &str) -> bool {
    let files = run_files(&run.name, run.dir.as_deref());
    if files.contains(rel)
        || crate::plugins::RUN_GLOBS
            .iter()
            .any(|g| pycompat::mtime::fnmatch(rel, g))
    {
        return true;
    }
    let (folder, base) = match rel.rsplit_once('/') {
        Some((f, b)) => (f, b),
        None => ("", rel),
    };
    !folder.is_empty() && sole_html(run, root, folder, &files).as_deref() == Some(base)
}

/// The name of the only HTML file directly in a run folder the report
/// links, inside the run and the root; `None` when the folder is not
/// linked or holds none or several.
pub fn sole_html(run: &Run, root: &Path, folder: &str, files: &BTreeSet<String>) -> Option<String> {
    let d = run.dir.as_deref()?;
    if !files.contains(folder) {
        return None;
    }
    let p = d.join(folder);
    if !p.is_dir() {
        return None;
    }
    let top = realpath(d);
    let pages: Vec<String> = fs::read_dir(&p)
        .ok()?
        .filter_map(|e| e.ok()?.file_name().into_string().ok())
        .filter(|f| {
            let fp = p.join(f);
            f.to_lowercase().ends_with(".html")
                && fp.is_file()
                && inside(&fp, &top)
                && inside(&fp, root)
        })
        .collect();
    match pages.as_slice() {
        [one] => Some(one.clone()),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fam(members: &[(&str, &str, &str)], unreadable: bool) -> Value {
        let ms = members
            .iter()
            .map(|(s, f, m)| {
                let mut o = IndexMap::new();
                o.insert("sample".into(), Value::Str((*s).into()));
                o.insert("father".into(), Value::Str((*f).into()));
                o.insert("mother".into(), Value::Str((*m).into()));
                Value::Object(o)
            })
            .collect();
        let mut o = IndexMap::new();
        o.insert("members".into(), Value::Array(ms));
        o.insert("unreadable".into(), Value::Bool(unreadable));
        Value::Object(o)
    }

    #[test]
    fn who_is_a_child() {
        let a = fam(&[("kid", "dad", "mum"), ("dad", "", "")], false);
        let b = fam(&[("kid", "", "")], false);
        assert!(child_of(&[a.clone(), b.clone()], "kid"));
        assert!(child_of(&[b.clone(), a.clone()], "kid"));
        assert!(!child_of(std::slice::from_ref(&a), "dad"));
        assert!(!child_of(std::slice::from_ref(&a), "loner"));
        // A family whose members cannot be read closes everyone.
        assert!(child_of(&[a.clone(), fam(&[], false)], "loner"));
        // One that lists the run but cannot be read closes it alone.
        let c = fam(&[("dad", "", "")], true);
        assert!(child_of(std::slice::from_ref(&c), "dad") && !child_of(&[c], "loner"));
    }

    #[test]
    fn report_links_relative_only() {
        let html = r##"<a href="clinvar/hits.tsv">x</a><a href="#top">t</a><a href="https://x.org/a">o</a>
<a href="/abs">a</a><a href="cpsr/">c</a><a href="a%20b.txt?x=1">q</a><a href="../out">u</a>
<a href="x&amp;y.tsv">e</a><a href="mailto:a@b">m</a><a href=".">d</a>"##;
        let got: Vec<String> = report_links(html).into_iter().collect();
        assert_eq!(got, ["a b.txt", "clinvar/hits.tsv", "cpsr", "x&y.tsv"]);
    }
}
