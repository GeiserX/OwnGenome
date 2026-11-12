//! The files a check writes: `recheck/state.json` (rewritten whole) and
//! `recheck/log.tsv` (only ever appended).

use crate::Check;
use pycompat::json::{Dump, Value, dumps, loads};
use std::fs::OpenOptions;
use std::io::{self, Write};
use std::path::Path;

pub const LOG_HEADER: &str = "utc\tevent\titem\tclass\tdetail\n";

fn invalid(e: impl std::fmt::Display) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, e.to_string())
}

fn dump() -> Dump<'static> {
    Dump {
        indent: Some(1),
        ..Dump::default()
    }
}

/// `obj` as JSON with indent 1, ASCII only, and a final newline, written
/// under a temporary name and renamed into place (mode 0666 less the umask).
pub fn write_json(path: &Path, obj: &Value) -> io::Result<()> {
    if let Some(d) = path.parent() {
        std::fs::create_dir_all(d)?;
    }
    let text = dumps(obj, &dump()).map_err(invalid)?;
    let mut f = og_io::AtomicFile::create(path)?;
    f.write_all(text.as_bytes())?;
    f.write_all(b"\n")?;
    f.commit()
}

/// A check as `write_json(path, &c.to_json())` writes it, byte for byte,
/// with the items written one at a time: the whole document never sits in
/// memory, however many PGS scores the run had. Inside the document an item
/// is two levels deep, so each of its lines gets two more spaces (a JSON
/// string holds no raw line break, so a line is always a whole line).
pub fn write_check(path: &Path, c: &Check) -> io::Result<()> {
    if let Some(d) = path.parent() {
        std::fs::create_dir_all(d)?;
    }
    let head = c.json_with(Value::Array(vec![]));
    let text = dumps(&head, &dump()).map_err(invalid)?;
    let mark = "\n \"items\": []";
    let at = text.find(mark).ok_or_else(|| invalid("no items key"))? + mark.len() - 2;
    let mut f = og_io::AtomicFile::create(path)?;
    let mut w = io::BufWriter::new(&mut f);
    w.write_all(&text.as_bytes()[..at])?;
    if c.items.is_empty() {
        w.write_all(b"[]")?;
    } else {
        w.write_all(b"[")?;
        for (n, item) in c.items.iter().enumerate() {
            w.write_all(if n == 0 { b"\n" } else { b",\n" })?;
            let one = dumps(&item.to_json(), &dump()).map_err(invalid)?;
            for (k, line) in one.split('\n').enumerate() {
                if k > 0 {
                    w.write_all(b"\n")?;
                }
                w.write_all(b"  ")?;
                w.write_all(line.as_bytes())?;
            }
        }
        w.write_all(b"\n ]")?;
    }
    w.write_all(&text.as_bytes()[at + 2..])?;
    w.write_all(b"\n")?;
    w.flush()?;
    drop(w);
    f.commit()
}

/// Appends rows to the log, one write per row, opened for append only, so
/// the log never shrinks; a new log gets its header first. A tab or line
/// break inside a cell becomes a space.
pub fn append_log(path: &Path, rows: &[[&str; 5]]) -> io::Result<()> {
    append_rows(path, rows.iter().map(|r| r.map(str::to_owned)))
}

fn append_rows(path: &Path, rows: impl Iterator<Item = [String; 5]>) -> io::Result<()> {
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

/// The log row of one item of a check.
fn check_row(c: &Check, i: &crate::Item) -> [String; 5] {
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
}

/// Appends the rows of a check, one per item, as each is made.
pub fn log_check(path: &Path, c: &Check) -> io::Result<()> {
    append_rows(path, c.items.iter().map(|i| check_row(c, i)))
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
