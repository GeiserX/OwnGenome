//! The database re-check on the live page.
//!
//! `og recheck` writes `<run>/recheck/state.json`: which databases the run
//! used have a newer release, which steps that would change, and the last
//! rerun. The browser may fetch that file (`/runs/<run>/recheck/state.json`,
//! [`crate::plugins::RUN_GLOBS`]), and `state.json` carries a short form of
//! it under `plugins.recheck` ([`state`]): a banner line and, per step, the
//! badges the step map shows ("newer ClinVar available", "redone
//! 2026-10-05"). `assets/live/plugins/recheck.js` draws them. Release names
//! and step numbers only: the check holds no result.

use indexmap::IndexMap;
use pycompat::json::{self, Value};
use pycompat::value::{py_str, truthy};
use std::fs;
use std::path::Path;

/// What each class of item says on the page.
fn act(class: &str, name: &str) -> Option<String> {
    match class {
        "rerun-ready" => Some(format!("a rerun picks up the newer {name}")),
        "refresh-available" => Some(format!("newer {name} available")),
        _ => None,
    }
}

fn text(v: &Value) -> String {
    match v {
        Value::Str(s) => s.clone(),
        other => py_str(other).unwrap_or_default(),
    }
}

fn get<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(key),
        _ => None,
    }
}

/// `v or ""` as text.
fn or_blank(v: Option<&Value>) -> String {
    v.filter(|v| truthy(v)).map(text).unwrap_or_default()
}

/// `'ClinVar (step 06's file)'` as `'ClinVar'`: a last part in brackets
/// cut, `database` for nothing left.
pub fn short(label: &str) -> String {
    let body = label.strip_suffix('\n').unwrap_or(label);
    let cut = label
        .char_indices()
        .filter(|(_, c)| *c == '(')
        .map(|(i, _)| i)
        .find(|&i| {
            let rest = &body.get(i..).unwrap_or("");
            rest.ends_with(')') && !rest.contains('\n')
        })
        .map_or(label, |i| &label[..i]);
    let t = cut.trim();
    if t.is_empty() {
        "database".into()
    } else {
        t.into()
    }
}

/// `'2026-09-02 12:00:00 UTC'` or `'2026-09-02T12:00:00Z'` as
/// `'20260902120000'`: the first 14 digits; `""` without any.
pub fn utc14(text: &str) -> String {
    text.chars().filter(char::is_ascii_digit).take(14).collect()
}

/// The `written_utc` row of `<run_dir>/run_manifest.tsv`, or `""`.
pub fn manifest_written(run_dir: &Path) -> String {
    let Ok(b) = fs::read(run_dir.join("run_manifest.tsv")) else {
        return String::new();
    };
    let t = String::from_utf8_lossy(&b)
        .replace("\r\n", "\n")
        .replace('\r', "\n");
    for line in t.split('\n') {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() >= 3 && c[0] == "run" && c[1] == "written_utc" {
            return c[2].to_owned();
        }
    }
    String::new()
}

/// A step number as a JSON key: text as it is, a number as Python writes
/// a key.
fn key(v: &Value) -> String {
    match v {
        Value::Bool(b) => (if *b { "true" } else { "false" }).into(),
        Value::Null => "null".into(),
        other => text(other),
    }
}

/// `plugins.recheck` of `state.json`: `{file, mtime (ns), checked_utc,
/// banner, steps: {step: [badge]}}`, or `{}` while there is no check, or
/// while the check is older than the run's `run_manifest.tsv` (then it
/// describes an earlier run, as the reports read it: `og run delta`).
pub fn state(run_dir: &Path) -> Value {
    let empty = Value::Object(IndexMap::new());
    let p = run_dir.join("recheck").join("state.json");
    let Ok(mtime) = pycompat::mtime::mtime_ns(&p) else {
        return empty;
    };
    let Some(st) = fs::read(&p)
        .ok()
        .and_then(|b| json::loads(&String::from_utf8_lossy(&b)).ok())
    else {
        return empty;
    };
    if !matches!(st, Value::Object(_)) {
        return empty;
    }
    let items: Vec<&Value> = match get(&st, "items") {
        Some(Value::Array(a)) => a.iter().filter(|i| matches!(i, Value::Object(_))).collect(),
        _ => Vec::new(),
    };
    let checked = get(&st, "checked_utc");
    if utc14(&manifest_written(run_dir)) > utc14(&or_blank(checked)) {
        return empty;
    }
    let mut steps: IndexMap<String, Vec<Value>> = IndexMap::new();
    let mut newer = Vec::new();
    let label = |i: &Value| {
        short(&or_blank(
            get(i, "label").filter(|v| truthy(v)).or(get(i, "item")),
        ))
    };
    for i in items {
        let class = get(i, "class").map(text).unwrap_or_default();
        let is_str = matches!(get(i, "class"), Some(Value::Str(_)));
        if let Some(t) = act(&class, &label(i)).filter(|_| is_str) {
            newer.push(t.clone());
            if let Some(Value::Array(ns)) = get(i, "steps").filter(|v| truthy(v)) {
                for n in ns {
                    steps.entry(key(n)).or_default().push(Value::Str(t.clone()));
                }
            }
        } else if is_str
            && ["in-pipeline-release", "upstream-only", "error"].contains(&class.as_str())
        {
            newer.push(format!("{}: {class}", label(i)));
        }
    }
    let rr = match get(&st, "rerun") {
        Some(v @ Value::Object(_)) => Some(v),
        _ => None,
    };
    let done: String = or_blank(rr.and_then(|r| get(r, "finished_utc")))
        .chars()
        .take(10)
        .collect();
    if !done.is_empty()
        && let Some(Value::Array(ns)) = rr
            .and_then(|r| get(r, "predicted_vs_redone"))
            .and_then(|p| get(p, "redone"))
            .filter(|v| truthy(v))
    {
        for n in ns {
            steps
                .entry(key(n))
                .or_default()
                .push(Value::Str(format!("redone {done}")));
        }
    }
    let when = checked
        .filter(|v| truthy(v))
        .map_or_else(|| "unknown date".to_owned(), text);
    let banner = format!(
        "Database re-check of {when}: {}",
        if newer.is_empty() {
            "nothing newer than this run".to_owned()
        } else {
            newer.join("; ")
        }
    );
    let mut o = IndexMap::new();
    o.insert("file".into(), Value::Str("recheck/state.json".into()));
    o.insert("mtime".into(), Value::Int(mtime));
    o.insert(
        "checked_utc".into(),
        checked.cloned().unwrap_or(Value::Null),
    );
    o.insert("banner".into(), Value::Str(banner));
    o.insert(
        "steps".into(),
        Value::Object(
            steps
                .into_iter()
                .map(|(k, v)| (k, Value::Array(v)))
                .collect(),
        ),
    );
    Value::Object(o)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_names() {
        assert_eq!(short("ClinVar (step 06's file)"), "ClinVar");
        assert_eq!(short("gnomAD"), "gnomAD");
        assert_eq!(short(""), "database");
        assert_eq!(short("(x)"), "database");
        assert_eq!(short("A (b) c"), "A (b) c");
        assert_eq!(short("A (b) (c)"), "A");
    }

    #[test]
    fn utc_digits() {
        assert_eq!(utc14("2026-09-02 12:00:00 UTC"), "20260902120000");
        assert_eq!(utc14("2026-09-02T12:00:00.123Z"), "20260902120000");
        assert_eq!(utc14(""), "");
    }
}
