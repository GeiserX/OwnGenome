//! The files a check writes: `recheck/state.json` (rewritten whole) and
//! `recheck/log.tsv` (only ever appended).

use crate::Check;
use pycompat::json::{Dump, Value, dumps, loads};
use std::fs::OpenOptions;
use std::io::{self, Write};
use std::path::Path;

pub const LOG_HEADER: &str = "utc\tevent\titem\tclass\tdetail\n";

/// `obj` as JSON with indent 1, ASCII only, and a final newline, written
/// under a temporary name and renamed into place (mode 0666 less the umask).
pub fn write_json(path: &Path, obj: &Value) -> io::Result<()> {
    if let Some(d) = path.parent() {
        std::fs::create_dir_all(d)?;
    }
    let text = dumps(
        obj,
        &Dump {
            indent: Some(1),
            ..Dump::default()
        },
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    let mut f = og_io::AtomicFile::create(path)?;
    f.write_all(text.as_bytes())?;
    f.write_all(b"\n")?;
    f.commit()
}

/// Appends rows to the log, one write per row, opened for append only, so
/// the log never shrinks; a new log gets its header first. A tab or line
/// break inside a cell becomes a space.
pub fn append_log(path: &Path, rows: &[[&str; 5]]) -> io::Result<()> {
    if let Some(d) = path.parent() {
        std::fs::create_dir_all(d)?;
    }
    let new = !path.exists();
    let mut f = OpenOptions::new().append(true).create(true).open(path)?;
    if new {
        f.write_all(LOG_HEADER.as_bytes())?;
    }
    for r in rows {
        let cells: Vec<String> = r.iter().map(|c| c.replace(['\t', '\n'], " ")).collect();
        f.write_all(format!("{}\n", cells.join("\t")).as_bytes())?;
    }
    Ok(())
}

/// A value as Python's f-string prints an optional one (`None` when absent).
fn py(v: &Option<String>) -> &str {
    v.as_deref().unwrap_or("None")
}

/// The log rows of a check: one per item.
pub fn check_rows(c: &Check) -> Vec<[String; 5]> {
    c.items
        .iter()
        .map(|i| {
            let detail = match &i.error {
                Some(e) if !e.is_empty() => e.clone(),
                _ => format!(
                    "{} {} {}",
                    py(&i.ran_with),
                    py(&i.installed),
                    py(&i.upstream.clone().flatten())
                ),
            };
            [
                c.checked_utc.clone(),
                "check".into(),
                i.item.clone(),
                i.class.as_str().into(),
                detail,
            ]
        })
        .collect()
}

/// Appends the rows of a check.
pub fn log_check(path: &Path, c: &Check) -> io::Result<()> {
    let rows = check_rows(c);
    let refs: Vec<[&str; 5]> = rows
        .iter()
        .map(|r| r.each_ref().map(String::as_str))
        .collect();
    append_log(path, &refs)
}

/// The `rerun` block of the last `state.json`, if it has one.
pub fn last_rerun(path: &Path) -> Option<Value> {
    let text = std::fs::read_to_string(path).ok()?;
    match loads(&text).ok()? {
        Value::Object(m) => m
            .get("rerun")
            .filter(|v| pycompat::value::truthy(v))
            .cloned(),
        _ => None,
    }
}
