//! TSV tables as OwnGenome's Python tools read and write them.
//!
//! - [`Options`] and [`Table`]: `csv.DictReader(f, delimiter="\t")` over a
//!   path or `-`, plain or gzip, through [`pycompat::csv::DictReader`], so
//!   short rows, extra cells, a header named twice and quoting behave as in
//!   Python. [`Options::skip`] drops lines before the reader sees them, as
//!   the tools' `(l for l in f if not l.startswith("#"))` does.
//! - [`Row`]: one row. [`Row::get`] is Python's `row.get(col)`; the typed
//!   accessors [`Row::text`], [`Row::float`] and [`Row::int`] parse as
//!   Python's `float()` and `int()` do, and their errors name the file, the
//!   line and the column (`sample.tsv:12: column depth is not a number:
//!   "n/a"`), where Python would stop with a traceback.
//! - [`KvTable`]: the two-column `metric\tvalue` table many tools write,
//!   written atomically.

use og_io::{AtomicFile, Lines};
use pycompat::csv::{self, DictReader, DictRow, FIELD_SIZE_LIMIT};
use pycompat::num;
use std::cell::Cell;
use std::fmt;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::rc::Rc;
use std::sync::Arc;

pub use og_io::Utf8;
pub use pycompat::csv::Quoting;

/// Why a table could not be read, or a cell not parsed. Each message starts
/// with the file, and with the line where there is one.
#[derive(Debug)]
pub enum Error {
    /// The file could not be opened, or a line read (invalid UTF-8
    /// included, when reading strictly): `line` is that line.
    Io {
        path: Arc<str>,
        line: Option<u64>,
        source: io::Error,
    },
    /// `csv.Error`: a row Python's reader refuses.
    Csv {
        path: Arc<str>,
        line: u64,
        source: csv::Error,
    },
    /// The file has no header line.
    NoHeader { path: Arc<str> },
    /// [`Table::require`]: the header lacks these columns.
    MissingColumns {
        path: Arc<str>,
        columns: Vec<String>,
    },
    /// A typed accessor asked for a column the header lacks.
    NoColumn {
        path: Arc<str>,
        line: u64,
        column: String,
    },
    /// The row ends before the column (Python's `None` cell).
    ShortRow {
        path: Arc<str>,
        line: u64,
        column: String,
    },
    /// The cell is not what the accessor parses.
    BadCell {
        path: Arc<str>,
        line: u64,
        column: String,
        value: String,
        want: &'static str,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io {
                path,
                line: None,
                source,
            } => write!(f, "{path}: {source}"),
            Error::Io {
                path,
                line: Some(line),
                source,
            } => write!(f, "{path}:{line}: {source}"),
            Error::Csv { path, line, source } => match source {
                csv::Error::Io(e) => write!(f, "{path}:{line}: {e}"),
                csv::Error::NewlineInUnquotedField { .. } => {
                    write!(
                        f,
                        "{path}:{line}: new-line character seen in unquoted field"
                    )
                }
                csv::Error::FieldTooLarge { .. } => write!(
                    f,
                    "{path}:{line}: field larger than field limit ({FIELD_SIZE_LIMIT})"
                ),
            },
            Error::NoHeader { path } => write!(f, "{path}: no header line"),
            Error::MissingColumns { path, columns } => {
                write!(f, "{path}: no column {} in the header", columns.join(", "))
            }
            Error::NoColumn { path, line, column } => {
                write!(f, "{path}:{line}: no column {column} in the header")
            }
            Error::ShortRow { path, line, column } => {
                write!(f, "{path}:{line}: the row ends before column {column}")
            }
            Error::BadCell {
                path,
                line,
                column,
                value,
                want,
            } => write!(f, "{path}:{line}: column {column} is not {want}: {value:?}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io { source, .. } => Some(source),
            Error::Csv { source, .. } => Some(source),
            _ => None,
        }
    }
}

/// How a [`Table`] reads: `quoting` and UTF-8 as Python's `open()` and
/// `DictReader` take them, and which lines to drop first.
#[derive(Debug, Clone, Copy)]
pub struct Options {
    quoting: Quoting,
    utf8: Utf8,
    skip: Option<fn(&str) -> bool>,
}

impl Default for Options {
    /// `open(path)` and `csv.DictReader(f, delimiter="\t")`: strict UTF-8,
    /// `QUOTE_MINIMAL`, every line read.
    fn default() -> Self {
        Options {
            quoting: Quoting::Minimal,
            utf8: Utf8::Strict,
            skip: None,
        }
    }
}

impl Options {
    /// `quoting=csv.QUOTE_NONE` or the default `QUOTE_MINIMAL`.
    pub fn quoting(self, quoting: Quoting) -> Self {
        Options { quoting, ..self }
    }

    /// `errors="replace"` or the default strict decoding.
    pub fn utf8(self, utf8: Utf8) -> Self {
        Options { utf8, ..self }
    }

    /// Drops each line for which `skip` is true, before the reader sees it.
    /// The line comes with its `\n`, as iterating a Python file gives it, so
    /// `|l| l.starts_with('#')` is `not l.startswith("#")` and
    /// `|l| pycompat::text::strip(l, None).is_empty()` is `x.strip()`. Line
    /// numbers in errors still count every line of the file.
    pub fn skip(self, skip: fn(&str) -> bool) -> Self {
        Options {
            skip: Some(skip),
            ..self
        }
    }

    /// Opens `path` (`-` is stdin), plain or gzip.
    pub fn open(self, path: impl AsRef<Path>) -> Result<Table, Error> {
        let path = path.as_ref();
        let name: Arc<str> = path.display().to_string().into();
        match og_io::open(path) {
            Ok(input) => Ok(self.read(name, input)),
            Err(source) => Err(Error::Io {
                path: name,
                line: None,
                source,
            }),
        }
    }

    /// Reads `input`, naming it `name` in errors.
    pub fn read(self, name: impl Into<Arc<str>>, input: impl BufRead + 'static) -> Table {
        let at = Rc::new(Cell::new(0));
        let lines = Fed {
            lines: Lines::new(Box::new(input), self.utf8),
            skip: self.skip,
            read: 0,
            at: Rc::clone(&at),
        };
        Table {
            path: name.into(),
            rows: DictReader::new(lines, '\t', self.quoting),
            at,
        }
    }
}

/// The lines a [`DictReader`] is fed: the file's, minus the skipped ones,
/// noting the file line number of the last one fed.
struct Fed {
    lines: Lines<Box<dyn BufRead>>,
    skip: Option<fn(&str) -> bool>,
    read: u64,
    at: Rc<Cell<u64>>,
}

impl Iterator for Fed {
    type Item = io::Result<String>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        loop {
            match self.lines.read_line(&mut line) {
                Ok(false) => return None,
                Ok(true) => {
                    self.read += 1;
                    if self.skip.is_some_and(|skip| skip(&line)) {
                        continue;
                    }
                    self.at.set(self.read);
                    return Some(Ok(line));
                }
                Err(e) => {
                    // The reader has moved past the bad line: count it, so
                    // a caller that keeps going gets the right line numbers.
                    self.read += 1;
                    self.at.set(self.read);
                    return Some(Err(e));
                }
            }
        }
    }
}

/// A reader's error as this crate's, at file line `line`.
fn csv_error(path: &Arc<str>, line: u64, source: csv::Error) -> Error {
    match source {
        csv::Error::Io(source) => Error::Io {
            path: Arc::clone(path),
            line: Some(line),
            source,
        },
        source => Error::Csv {
            path: Arc::clone(path),
            line,
            source,
        },
    }
}

/// A TSV table read row by row, `csv.DictReader(f, delimiter="\t")`.
pub struct Table {
    path: Arc<str>,
    rows: DictReader<Fed>,
    at: Rc<Cell<u64>>,
}

impl Table {
    /// `csv.DictReader(open(path), delimiter="\t")`.
    pub fn open(path: impl AsRef<Path>) -> Result<Table, Error> {
        Options::default().open(path)
    }

    /// The name errors give the file.
    pub fn path(&self) -> &str {
        &self.path
    }

    /// `fieldnames`: the header cells, duplicates included; `None` for an
    /// empty file.
    pub fn header(&mut self) -> Result<Option<&[String]>, Error> {
        let (path, at) = (Arc::clone(&self.path), Rc::clone(&self.at));
        match self.rows.header() {
            Ok(h) => Ok(h.map(csv::Header::names)),
            Err(e) => Err(csv_error(&path, at.get(), e)),
        }
    }

    /// Checks that the header has every one of `columns`.
    pub fn require(&mut self, columns: &[&str]) -> Result<(), Error> {
        let path = Arc::clone(&self.path);
        let Some(names) = self.header()? else {
            return Err(Error::NoHeader { path });
        };
        let missing: Vec<String> = columns
            .iter()
            .filter(|c| !names.iter().any(|n| n == *c))
            .map(|c| (*c).to_owned())
            .collect();
        if missing.is_empty() {
            Ok(())
        } else {
            Err(Error::MissingColumns {
                path,
                columns: missing,
            })
        }
    }
}

impl Iterator for Table {
    type Item = Result<Row, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        Some(match self.rows.next()? {
            Ok(dict) => Ok(Row {
                path: Arc::clone(&self.path),
                line: self.at.get(),
                dict,
            }),
            Err(e) => Err(csv_error(&self.path, self.at.get(), e)),
        })
    }
}

/// One row of a [`Table`].
#[derive(Debug, Clone)]
pub struct Row {
    path: Arc<str>,
    line: u64,
    dict: DictRow,
}

impl Row {
    /// The file line the row ends on (a quoted field can span lines).
    pub fn line(&self) -> u64 {
        self.line
    }

    /// The row as Python's dict sees it.
    pub fn dict(&self) -> &DictRow {
        &self.dict
    }

    /// `row.get(column)`: `None` when the header lacks the column or the
    /// row ends before it.
    pub fn get(&self, column: &str) -> Option<&str> {
        self.dict.get(column).flatten()
    }

    /// `row[column]`, an error naming the line when there is no such cell.
    pub fn text(&self, column: &str) -> Result<&str, Error> {
        match self.dict.get(column) {
            Some(Some(cell)) => Ok(cell),
            Some(None) => Err(Error::ShortRow {
                path: Arc::clone(&self.path),
                line: self.line,
                column: column.to_owned(),
            }),
            None => Err(Error::NoColumn {
                path: Arc::clone(&self.path),
                line: self.line,
                column: column.to_owned(),
            }),
        }
    }

    fn bad(&self, column: &str, value: &str, want: &'static str) -> Error {
        Error::BadCell {
            path: Arc::clone(&self.path),
            line: self.line,
            column: column.to_owned(),
            value: value.to_owned(),
            want,
        }
    }

    /// `float(row[column])`: Python's spelling, so `" 1e3 "`, `"nan"` and
    /// `"-inf"` are numbers.
    pub fn float(&self, column: &str) -> Result<f64, Error> {
        let cell = self.text(column)?;
        num::parse_float(cell).map_err(|_| self.bad(column, cell, "a number"))
    }

    /// `int(row[column])`, refused beyond the 64-bit range.
    pub fn int(&self, column: &str) -> Result<i64, Error> {
        let cell = self.text(column)?;
        num::parse_int(cell)
            .ok()
            .and_then(|v| i64::try_from(v).ok())
            .ok_or_else(|| self.bad(column, cell, "an integer"))
    }
}

/// A two-column table: a header line, then one `key\tvalue` line per row,
/// in the order pushed. A tab, `\n` or `\r` inside a cell becomes a space,
/// so no value can break the table.
#[derive(Debug, Clone)]
pub struct KvTable {
    header: (String, String),
    rows: Vec<(String, String)>,
}

fn cell(s: &str) -> String {
    s.replace(['\t', '\n', '\r'], " ")
}

impl KvTable {
    /// A table whose header line is `key\tvalue`.
    pub fn new(key: &str, value: &str) -> Self {
        KvTable {
            header: (cell(key), cell(value)),
            rows: Vec::new(),
        }
    }

    /// The `metric\tvalue` table.
    pub fn metrics() -> Self {
        KvTable::new("metric", "value")
    }

    /// Adds a row.
    pub fn push(&mut self, key: &str, value: &str) -> &mut Self {
        self.rows.push((cell(key), cell(value)));
        self
    }

    /// Writes the table to `out`.
    pub fn write_to(&self, mut out: impl Write) -> io::Result<()> {
        let (k, v) = &self.header;
        writeln!(out, "{k}\t{v}")?;
        for (k, v) in &self.rows {
            writeln!(out, "{k}\t{v}")?;
        }
        out.flush()
    }

    /// Writes the table to `path` under a temporary name and renames it
    /// into place, so a reader never sees half of it.
    pub fn write(&self, path: impl AsRef<Path>) -> io::Result<()> {
        let mut file = AtomicFile::create(path)?;
        self.write_to(&mut file)?;
        file.commit()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(text: &str, opts: Options) -> Table {
        opts.read("t.tsv", io::Cursor::new(text.as_bytes().to_vec()))
    }

    fn rows(t: Table) -> Vec<Row> {
        t.collect::<Result<_, _>>().unwrap()
    }

    #[test]
    fn rows_are_keyed_by_the_header() {
        let mut t = table("a\tb\n1\t2\n\n3\t4\textra\n5\n", Options::default());
        assert_eq!(t.header().unwrap().unwrap(), ["a", "b"]);
        let r = rows(t);
        assert_eq!(r.len(), 3);
        assert_eq!((r[0].get("a"), r[0].get("b")), (Some("1"), Some("2")));
        assert_eq!(r[1].dict().extra(), Some(&["extra".to_owned()][..]));
        // A short row: Python's None, which get() folds with a missing column.
        assert_eq!((r[2].get("b"), r[2].get("zz")), (None, None));
        // The blank line counts: the rows sit on lines 2, 4 and 5.
        assert_eq!(r.iter().map(Row::line).collect::<Vec<_>>(), [2, 4, 5]);
    }

    #[test]
    fn typed_cells_name_file_line_and_column() {
        let r = rows(table(
            "depth\tn\tname\n 12.5 \t7\tx\nn/a\t1e3\n",
            Options::default(),
        ));
        assert_eq!(r[0].float("depth").unwrap(), 12.5);
        assert_eq!(r[0].int("n").unwrap(), 7);
        assert_eq!(r[0].text("name").unwrap(), "x");
        let e = r[1].float("depth").unwrap_err();
        assert_eq!(
            e.to_string(),
            r#"t.tsv:3: column depth is not a number: "n/a""#
        );
        let e = r[1].int("n").unwrap_err();
        assert_eq!(
            e.to_string(),
            r#"t.tsv:3: column n is not an integer: "1e3""#
        );
        let e = r[1].text("name").unwrap_err();
        assert_eq!(e.to_string(), "t.tsv:3: the row ends before column name");
        let e = r[0].float("zz").unwrap_err();
        assert_eq!(e.to_string(), "t.tsv:2: no column zz in the header");
        assert!(r[0].float("n").is_ok());
        let big = rows(table("n\n9223372036854775808\n", Options::default()));
        assert!(big[0].int("n").is_err());
    }

    #[test]
    fn skipped_lines_still_count() {
        let r = rows(table(
            "## banner\n#\na\tb\n# note\n1\t2\n",
            Options::default().skip(|l| l.starts_with('#')),
        ));
        assert_eq!(r.len(), 1);
        assert_eq!((r[0].get("a"), r[0].line()), (Some("1"), 5));
    }

    #[test]
    fn quoting_follows_python() {
        let text = "a\tb\n\"x\ty\"\t2\n";
        let r = rows(table(text, Options::default()));
        assert_eq!(r[0].get("a"), Some("x\ty"));
        let r = rows(table(text, Options::default().quoting(Quoting::None)));
        assert_eq!((r[0].get("a"), r[0].get("b")), (Some("\"x"), Some("y\"")));
        // A quoted field over two lines: the row ends on line 3.
        let r = rows(table("a\n\"x\ny\"\n", Options::default()));
        assert_eq!((r[0].get("a"), r[0].line()), (Some("x\ny"), 3));
    }

    #[test]
    fn header_checks() {
        let mut t = table("", Options::default());
        assert!(t.header().unwrap().is_none());
        assert_eq!(
            t.require(&["a"]).unwrap_err().to_string(),
            "t.tsv: no header line"
        );
        let mut t = table("a\tb\n", Options::default());
        assert!(t.require(&["b", "a"]).is_ok());
        let e = t.require(&["a", "c", "d"]).unwrap_err();
        assert_eq!(e.to_string(), "t.tsv: no column c, d in the header");
    }

    #[test]
    fn bad_input_is_an_error_with_its_line() {
        let mut t =
            Options::default().read("t.tsv", io::Cursor::new(b"a\nok\n\xff\nnext\n".to_vec()));
        let _ = t.next().unwrap().unwrap();
        let e = t.next().unwrap().unwrap_err();
        assert!(matches!(e, Error::Io { line: Some(3), .. }), "{e}");
        // A caller that keeps going past the error keeps the right lines.
        assert_eq!(t.next().unwrap().unwrap().line(), 4);
        let mut t = Options::default()
            .utf8(Utf8::Replace)
            .read("t.tsv", io::Cursor::new(b"a\n\xff\n".to_vec()));
        assert_eq!(t.next().unwrap().unwrap().get("a"), Some("\u{fffd}"));
        let long = format!("a\n{}\n", "x".repeat(FIELD_SIZE_LIMIT + 1));
        let e = table(&long, Options::default())
            .next()
            .unwrap()
            .unwrap_err();
        assert_eq!(
            e.to_string(),
            "t.tsv:2: field larger than field limit (131072)"
        );
        let e = Table::open("/nonexistent/x.tsv").err().unwrap();
        assert!(e.to_string().starts_with("/nonexistent/x.tsv: "));
    }

    #[test]
    fn kv_table_writes_metric_value() {
        let mut t = KvTable::metrics();
        t.push("status", "called").push("note", "a\tb\nc\rd");
        let mut out = Vec::new();
        t.write_to(&mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "metric\tvalue\nstatus\tcalled\nnote\ta b c d\n"
        );
        let dir = std::env::temp_dir().join(format!("og-tsv-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("kv.tsv");
        KvTable::new("key", "value")
            .push("k", "v")
            .write(&path)
            .unwrap();
        // What the table writes, the reader reads back.
        let r = rows(Table::open(&path).unwrap());
        assert_eq!((r[0].get("key"), r[0].get("value")), (Some("k"), Some("v")));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
