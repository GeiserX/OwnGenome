//! The small readers the clinical sections (group b) share: a VCF's data
//! lines as cells, INFO and FORMAT as maps, a table whose read errors keep
//! their kind, and Python's way of reading a digit cell.

use crate::{basename, os_error, text_lines};
use pycompat::csv::{DictReader, DictRow, Quoting};
use pycompat::text::{isdigit, strip};
use std::collections::HashMap;

/// Why a file could not be read: an `OSError` (its `str()`), or a fault of
/// another kind, already written as the section's note.
#[derive(Debug, Clone)]
pub enum ReadError {
    Os(String),
    Other(String),
}

impl ReadError {
    /// The note when the Python wrapped the error in `Unreadable(f"{name}: {e}")`.
    pub fn unreadable(&self, path: &str) -> String {
        match self {
            ReadError::Os(e) => format!("{}: {e}", basename(path)),
            ReadError::Other(n) => n.clone(),
        }
    }

    /// The note of the collector's catch-all: `type(e).__name__: e`.
    pub fn note(&self) -> String {
        match self {
            ReadError::Os(e) => format!("{}: {e}", os_kind(e)),
            ReadError::Other(n) => n.clone(),
        }
    }

    /// `str(e)`, as a message that embeds the error prints it.
    pub fn text(&self) -> String {
        match self {
            ReadError::Os(e) => e.clone(),
            ReadError::Other(n) => n.split_once(": ").map_or(n.clone(), |(_, m)| m.to_owned()),
        }
    }
}

fn os_kind(e: &str) -> &'static str {
    if e == GZ_EOF {
        "EOFError"
    } else if e.contains("No such file") {
        "FileNotFoundError"
    } else if e.contains("Permission denied") {
        "PermissionError"
    } else if e.contains("Is a directory") {
        "IsADirectoryError"
    } else {
        "OSError"
    }
}

/// The lines of a text file (gzip by its magic, universal newlines, bad
/// UTF-8 as U+FFFD), each without its line ending.
pub fn lines(path: &str) -> Result<Vec<String>, ReadError> {
    let it = og_io::Lines::new(
        og_io::open(path).map_err(|e| ReadError::Os(os_error(&e, path)))?,
        og_io::Utf8::Replace,
    );
    it.map(|l| l.map_err(|e| ReadError::Os(gz_text(&e))))
        .collect()
}

/// Python's gzip on a stream cut short (an `EOFError`).
const GZ_EOF: &str = "Compressed file ended before the end-of-stream marker was reached";

/// gzip's own words for a stream cut short; the reader's otherwise.
fn gz_text(e: &std::io::Error) -> String {
    if e.kind() == std::io::ErrorKind::UnexpectedEof {
        GZ_EOF.to_owned()
    } else {
        e.to_string()
    }
}

/// Calls `f` with each line of a text file in turn (as [`lines`] reads
/// them), never holding the whole file.
pub fn each_line(path: &str, mut f: impl FnMut(&str)) -> Result<(), ReadError> {
    let it = og_io::Lines::new(
        og_io::open(path).map_err(|e| ReadError::Os(os_error(&e, path)))?,
        og_io::Utf8::Replace,
    );
    for l in it {
        f(&l.map_err(|e| ReadError::Os(gz_text(&e)))?);
    }
    Ok(())
}

/// The rows of a tab-separated table with a header (`csv.DictReader`).
pub fn rows(path: &str) -> Result<Vec<DictRow>, ReadError> {
    let mut lines = text_lines(path).map_err(ReadError::Os)?;
    let iter = std::iter::from_fn(move || {
        let mut s = String::new();
        match lines.read_line(&mut s) {
            Ok(true) => Some(Ok(s)),
            Ok(false) => None,
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                Some(Err(std::io::Error::new(e.kind(), gz_text(&e))))
            }
            Err(e) => Some(Err(e)),
        }
    });
    DictReader::new(iter, '\t', Quoting::Minimal)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| match e.exception() {
            "OSError" => ReadError::Os(e.to_string()),
            kind => ReadError::Other(format!("{kind}: {e}")),
        })
}

/// The header of a table as `csv.reader`'s first row reads it.
pub fn columns(path: &str) -> Result<Vec<String>, ReadError> {
    let first = lines(path)?.into_iter().next().unwrap_or_default();
    if first.is_empty() {
        return Ok(Vec::new());
    }
    Ok(first.split('\t').map(str::to_owned).collect())
}

/// `row.get(key)` where a missing key and a short row both read as `None`.
pub fn cell<'a>(row: &'a DictRow, key: &str) -> Option<&'a str> {
    row.get(key).flatten()
}

/// `row.get(key) or ""`.
pub fn cell_or<'a>(row: &'a DictRow, key: &str) -> &'a str {
    cell(row, key).unwrap_or("")
}

/// `(row.get(key) or "").strip()`.
pub fn stripped<'a>(row: &'a DictRow, key: &str) -> &'a str {
    strip(cell_or(row, key), None)
}

/// `row[key]`: a `KeyError` note when the table has no such column. A
/// short row gives `None`.
pub fn need<'a>(row: &'a DictRow, key: &str) -> Result<Option<&'a str>, String> {
    row.get(key)
        .ok_or_else(|| format!("KeyError: {}", crate::repr(key)))
}

/// `int(x) if x.isdigit() else None`.
pub fn digits(x: Option<&str>) -> Option<i128> {
    x.filter(|v| isdigit(v)).and_then(|v| crate::int(v).ok())
}

/// `info_map`: a VCF INFO field as a map, a later key over an earlier one.
pub fn info_map(field: &str) -> HashMap<&str, &str> {
    field
        .split(';')
        .map(|kv| kv.split_once('=').unwrap_or((kv, "")))
        .collect()
}

/// `dict(zip(FORMAT.split(":"), SAMPLE.split(":")))`.
pub fn format_map<'a>(fmt: &'a str, sample: &'a str) -> HashMap<&'a str, &'a str> {
    fmt.split(':').zip(sample.split(':')).collect()
}

/// One FORMAT field of the first sample, `None` when absent, empty or `.`.
pub fn gt_field(r: &[String], key: &str) -> Option<String> {
    if r.len() < 10 {
        return None;
    }
    format_map(&r[8], &r[9])
        .get(key)
        .filter(|v| !v.is_empty() && **v != ".")
        .map(|v| (*v).to_owned())
}

/// The GT of the first sample, `""` without one.
pub fn gt_of(r: &[String]) -> &str {
    r.get(9).map_or("", |s| s.split(':').next().unwrap_or(""))
}

/// `int(r[1])`, as the note of the collector's catch-all when it is not one.
pub fn pos(r: &[String]) -> Result<i128, String> {
    crate::int(&r[1]).map_err(|e| format!("ValueError: {e}"))
}

/// Python's `max()` over floats: the first of the largest, NaN kept only
/// when it comes first.
pub fn py_max(xs: impl IntoIterator<Item = f64>) -> Option<f64> {
    xs.into_iter().fold(None, |best, x| match best {
        Some(b) if x.partial_cmp(&b) != Some(std::cmp::Ordering::Greater) => Some(b),
        _ => Some(x),
    })
}

/// `os.path.relpath(path, base)` for a path found under `base`.
pub fn relpath(base: &str, path: &str) -> String {
    let b = base.trim_end_matches('/');
    path.strip_prefix(b)
        .and_then(|r| r.strip_prefix('/'))
        .unwrap_or(path)
        .to_owned()
}

/// `os.path.dirname(p)`.
pub fn dirname(p: &str) -> &str {
    p.rfind('/').map_or("", |i| &p[..i])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    #[test]
    fn a_table_cut_short_is_python_s_eof_error() {
        let d = Dir::new("text-cut");
        let body = "A\tB\n".to_owned() + &"1\t2\n".repeat(2000);
        let full = std::fs::read(d.put_gz("t.tsv.gz", &body)).unwrap();
        for n in [30, full.len() / 2, full.len() - 4] {
            let p = d.put("cut.tsv.gz", &full[..n]);
            let e = rows(&p).unwrap_err();
            assert_eq!(e.unreadable(&p), format!("cut.tsv.gz: {GZ_EOF}"), "{n}");
            assert_eq!(e.note(), format!("EOFError: {GZ_EOF}"), "{n}");
        }
    }
}
