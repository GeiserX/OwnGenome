//! What every reader of group c shares: the table reader with its header,
//! the cell lookups with the Python collector's three answers (no such
//! column, a row cut short, a cell), and an ordered JSON object builder.

use crate::collect::sections::{num, py_error};
use crate::{basename, repr, text_lines};
use indexmap::IndexMap;
use pycompat::csv::{DictReader, DictRow, Quoting};
use pycompat::json::Value;

/// A table read as `csv.DictReader` reads it, with its header cells.
pub struct Table {
    /// The header cells as written; empty for an empty file.
    pub head: Vec<String>,
    pub rows: Vec<DictRow>,
}

impl Table {
    /// Whether the header names `col` (`col in reader.fieldnames`).
    pub fn has(&self, col: &str) -> bool {
        self.head.iter().any(|h| h == col)
    }

    /// Whether the first row has every column (`need <= set(rows[0])`):
    /// false for a table without rows.
    pub fn first_has(&self, need: &[&str]) -> bool {
        !self.rows.is_empty() && need.iter().all(|k| self.has(k))
    }
}

/// Reads a tab-separated table with a header. `quote_none` reads `"` as a
/// plain character (`quoting=csv.QUOTE_NONE`).
pub fn table(path: &str, quote_none: bool) -> Result<Table, String> {
    let mut lines = text_lines(path).map_err(|e| os_note(&e))?;
    let iter = std::iter::from_fn(move || {
        let mut s = String::new();
        match lines.read_line(&mut s) {
            Ok(true) => Some(Ok(s)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    });
    let q = if quote_none {
        Quoting::None
    } else {
        Quoting::Minimal
    };
    let mut rd = DictReader::new(iter, '\t', q);
    let csv_err = |e: pycompat::csv::Error| py_error(e.exception(), &e.to_string());
    let head = rd
        .header()
        .map_err(csv_err)?
        .map(|h| h.names().to_vec())
        .unwrap_or_default();
    let rows = rd.collect::<Result<Vec<_>, _>>().map_err(csv_err)?;
    Ok(Table { head, rows })
}

/// The lines of a text file, line ends removed (`open(p, errors="replace")`).
pub fn lines(path: &str) -> Result<Vec<String>, String> {
    let mut rd = text_lines(path).map_err(|e| os_note(&e))?;
    let mut out = Vec::new();
    let mut s = String::new();
    while rd.read_line(&mut s).map_err(|e| {
        py_error(
            "OSError",
            &crate::collect::sections::gzip_message(&e),
        )
    })? {
        out.push(s.strip_suffix('\n').unwrap_or(&s).to_owned());
        s.clear();
    }
    Ok(out)
}

/// An `OSError` from opening a file, as the collector's catch-all wrote it.
pub fn os_note(e: &str) -> String {
    let kind = if e.contains("No such file") {
        "FileNotFoundError"
    } else if e.contains("Permission denied") {
        "PermissionError"
    } else if e.contains("Is a directory") {
        "IsADirectoryError"
    } else {
        "OSError"
    };
    py_error(kind, e)
}

/// `KeyError: 'k'`, a column a row lacks.
pub fn key_error(k: &str) -> String {
    py_error("KeyError", &repr(k))
}

/// `row[k]`: an error when the header has no such column, `None` for a
/// row cut short.
pub fn at<'a>(r: &'a DictRow, k: &str) -> Result<Option<&'a str>, String> {
    r.get(k).ok_or_else(|| key_error(k))
}

/// `row[k]` that must be a text: a short row is an error too.
pub fn text<'a>(r: &'a DictRow, k: &str) -> Result<&'a str, String> {
    at(r, k)?.ok_or_else(|| {
        py_error(
            "TypeError",
            &format!("the cell {} of a short row is None", repr(k)),
        )
    })
}

/// `row.get(k) or "."`.
pub fn or_dot(r: &DictRow, k: &str) -> String {
    or(r, k, ".")
}

/// `row.get(k) or default`.
pub fn or(r: &DictRow, k: &str, default: &str) -> String {
    match r.get(k).flatten() {
        Some(v) if !v.is_empty() => v.to_owned(),
        _ => default.to_owned(),
    }
}

/// `row.get(k, default)`: the default without the column, `null` for a
/// row cut short.
pub fn get_default(r: &DictRow, k: &str, default: &str) -> Value {
    match r.get(k) {
        None => s(default),
        Some(None) => Value::Null,
        Some(Some(v)) => s(v),
    }
}

/// Whether a row has a cell past the header (Python's `None` key) or is
/// cut short of a column of `need`.
pub fn ragged(r: &DictRow, need: &[&str]) -> bool {
    r.extra().is_some() || need.iter().any(|k| !matches!(r.get(k), Some(Some(_))))
}

/// `str.isdigit()` of an ASCII cell.
pub fn isdigit(x: &str) -> bool {
    pycompat::text::isdigit(x)
}

/// `float(x)` that is a number and not NaN.
pub fn number(x: &str) -> Option<f64> {
    num(Some(x))
}

/// `os.path.basename`.
pub fn base(p: &str) -> &str {
    basename(p)
}

pub fn s(x: &str) -> Value {
    Value::Str(x.to_owned())
}

pub fn opt_s(x: Option<&str>) -> Value {
    x.map_or(Value::Null, s)
}

pub fn strs<I, S>(xs: I) -> Value
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    Value::Array(xs.into_iter().map(|x| s(x.as_ref())).collect())
}

pub fn int(n: impl Into<i128>) -> Value {
    Value::Int(n.into())
}

pub fn count(n: usize) -> Value {
    Value::Int(n as i128)
}

/// An object built in insertion order.
#[derive(Debug, Default, Clone)]
pub struct Obj(pub IndexMap<String, Value>);

impl Obj {
    pub fn new() -> Obj {
        Obj::default()
    }

    pub fn set(mut self, k: &str, v: Value) -> Obj {
        self.0.insert(k.to_owned(), v);
        self
    }

    pub fn put(&mut self, k: &str, v: Value) {
        self.0.insert(k.to_owned(), v);
    }

    pub fn val(self) -> Value {
        Value::Object(self.0)
    }
}

/// `{k: row.get(k) or "." for k in keys}`.
pub fn dots(r: &DictRow, keys: &[&str]) -> Value {
    keys.iter()
        .fold(Obj::new(), |o, k| o.set(k, s(&or_dot(r, k))))
        .val()
}

/// `{k: row[k] or "." for k in keys}`: a column the row lacks is an error.
pub fn dots_strict(r: &DictRow, keys: &[&str]) -> Result<Value, String> {
    let mut o = Obj::new();
    for k in keys {
        let v = at(r, k)?.filter(|v| !v.is_empty()).unwrap_or(".");
        o.put(k, s(v));
    }
    Ok(o.val())
}

/// `{r[key]: r.get(value) for r in rows}` of a `metric\tvalue` table: the
/// last row of a key wins, in the place of its first.
pub fn metric_map(t: &Table, key: &str, value: &str) -> IndexMap<Option<String>, Option<String>> {
    let mut out = IndexMap::new();
    for r in &t.rows {
        let k = r.get(key).flatten().map(str::to_owned);
        let v = r.get(value).flatten().map(str::to_owned);
        out.insert(k, v);
    }
    out
}

/// `kv.get(k) or "."` over a [`metric_map`].
pub fn kv_dot(kv: &IndexMap<Option<String>, Option<String>>, k: &str) -> String {
    match kv.get(&Some(k.to_owned())) {
        Some(Some(v)) if !v.is_empty() => v.clone(),
        _ => ".".to_owned(),
    }
}

/// `kv.get(k)` over a [`metric_map`].
pub fn kv_get<'a>(kv: &'a IndexMap<Option<String>, Option<String>>, k: &str) -> Option<&'a str> {
    kv.get(&Some(k.to_owned())).and_then(|v| v.as_deref())
}

/// `[(metric, value)]` of a two-column table with a header line
/// (`read_tsv_pairs`): blank lines skipped, short lines left out.
pub fn pairs(path: &str) -> Result<Vec<(String, String)>, String> {
    let ls = lines(path)?;
    let ls: Vec<&String> = ls
        .iter()
        .filter(|l| !pycompat::num::py_strip(l).is_empty())
        .collect();
    Ok(ls
        .iter()
        .skip(1)
        .filter_map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c.len() >= 2).then(|| (c[0].to_owned(), c[1].to_owned()))
        })
        .collect())
}

/// `round(x, n)` as a JSON float.
pub fn round(x: f64, n: i32) -> Value {
    Value::Float(pycompat::fmt::round(x, n).unwrap_or(x))
}

/// The rows of `assets/<name>`, comment lines left out (`asset_rows`).
pub fn asset_rows(assets: &str, name: &str) -> Result<Vec<DictRow>, String> {
    let p = crate::join(assets, name);
    if !crate::is_file(&p) {
        return Err(format!(
            "{name} is not in {assets}: the report needs the repository's assets/ folder (render_report.py --assets)"
        ));
    }
    let text = std::fs::read(&p).map_err(|e| py_error("OSError", &e.to_string()))?;
    let text =
        String::from_utf8(text).map_err(|e| py_error("UnicodeDecodeError", &e.to_string()))?;
    let lines = text
        .split_inclusive('\n')
        .filter(|l| !l.starts_with('#'))
        .map(|l| Ok::<_, std::io::Error>(l.to_owned()));
    DictReader::new(lines, '\t', Quoting::None)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| py_error(e.exception(), &e.to_string()))
}

/// The header of a table as `csv.reader` reads its first row (empty for an
/// empty file).
pub fn columns(path: &str) -> Result<Vec<String>, String> {
    Ok(table(path, false)?.head)
}

/// The path of a table the Python collector read beside itself in `bin/`
/// (`hla_alleles.tsv`, `str_intermediate_ranges.tsv`): the `bin/` of the
/// checkout the binary sits in (`bin/.og/og-*` and `rust/target/*/og` both
/// do). `OG_BIN_TABLES` names another folder.
pub fn bin_table(name: &str) -> String {
    if let Some(d) = std::env::var_os("OG_BIN_TABLES") {
        return std::path::Path::new(&d).join(name).display().to_string();
    }
    let exe = std::env::current_exe().unwrap_or_default();
    exe.ancestors()
        .skip(1)
        .flat_map(|a| [a.join(name), a.join("bin").join(name)])
        .find(|p| p.is_file())
        .unwrap_or_else(|| std::path::Path::new("bin").join(name))
        .display()
        .to_string()
}

/// The data rows of a `bin/` table: `#` and blank lines left out, the first
/// row the header, each row zipped with it (a short row lacks the keys).
pub fn bin_rows(path: &str) -> Result<Vec<indexmap::IndexMap<String, String>>, String> {
    let text = std::fs::read(path).map_err(|e| {
        py_error(
            if e.kind() == std::io::ErrorKind::NotFound {
                "FileNotFoundError"
            } else {
                "OSError"
            },
            &crate::os_error(&e, path),
        )
    })?;
    let text =
        String::from_utf8(text).map_err(|e| py_error("UnicodeDecodeError", &e.to_string()))?;
    let mut head: Option<Vec<&str>> = None;
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || pycompat::num::py_strip(line).is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        match &head {
            None => head = Some(c),
            Some(h) => out.push(
                h.iter()
                    .zip(c)
                    .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                    .collect(),
            ),
        }
    }
    Ok(out)
}

/// Each data line of a VCF, plain or gzip (BGZF inflated on many cores),
/// without its line end; `#` lines are left out. Universal newlines, as
/// Python reads text. An error is `<file>: <reason>`.
pub fn vcf_lines(path: &str, mut each: impl FnMut(&[u8])) -> Result<(), String> {
    use std::io::BufRead;
    let fail = |e: std::io::Error| {
        format!(
            "{}: {}",
            base(path),
            crate::collect::sections::gzip_message(&e)
        )
    };
    let mut r = og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(fail)?;
    let mut buf = Vec::with_capacity(1 << 16);
    let mut one = |l: &[u8]| {
        if l.first() != Some(&b'#') {
            each(l);
        }
    };
    loop {
        buf.clear();
        if r.read_until(b'\n', &mut buf).map_err(fail)? == 0 {
            break;
        }
        let line = buf.strip_suffix(b"\n").unwrap_or(&buf);
        if memchr::memchr(b'\r', line).is_none() {
            one(line);
        } else {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            line.split(|&c| c == b'\r').for_each(&mut one);
        }
    }
    Ok(())
}

/// A VCF line's cells, as text (`errors="replace"`).
pub fn cells(line: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(line)
        .split('\t')
        .map(str::to_owned)
        .collect()
}
