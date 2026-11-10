//! Python's `csv` module, as OwnGenome's tools use it.
//!
//! - [`Writer`]: `csv.writer(f, delimiter=D, lineterminator="\n")` with the
//!   default `QUOTE_MINIMAL`.
//! - [`Reader`]: `csv.reader(lines, delimiter=D)` with `QUOTE_MINIMAL` or
//!   `QUOTE_NONE`, a port of CPython's `_csv.c` state machine.
//! - [`DictReader`]: `csv.DictReader` over a [`Reader`]: a short row fills
//!   the missing keys with `None`, extra cells go under the key `None`, and a
//!   header named twice keeps one key whose value is its last column.
//!
//! The other dialect settings keep their defaults: `quotechar='"'`,
//! `doublequote=True`, no `escapechar`, `skipinitialspace=False`,
//! `strict=False`.
//!
//! A reader takes lines *with* their endings, as iterating a Python file
//! yields them (`og_io::Lines::read_line`): a quoted field that spans lines
//! keeps the newline, and only the ending says where it was.

use std::collections::HashMap;
use std::fmt;
use std::io::{self, Write};
use std::sync::Arc;

const QUOTE: char = '"';

/// `csv.field_size_limit()`'s default: a longer field is an error.
pub const FIELD_SIZE_LIMIT: usize = 131_072;

/// The `quoting` of a reader.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quoting {
    /// `csv.QUOTE_MINIMAL` (the default): `"` opens a quoted field.
    Minimal,
    /// `csv.QUOTE_NONE`: `"` is a character like any other.
    None,
}

/// Python's `csv.writer(f, delimiter=D, lineterminator="\n")`.
pub struct Writer<W> {
    inner: W,
    delimiter: char,
    buf: String,
}

impl<W: Write> Writer<W> {
    pub fn new(inner: W, delimiter: char) -> Self {
        Writer {
            inner,
            delimiter,
            buf: String::new(),
        }
    }

    /// `writerow(fields)`. A field is quoted when it holds the delimiter,
    /// `"`, `\n` or `\r`, and a `"` inside is doubled. A row of one empty
    /// field is written `""`, so it does not read back as a blank line.
    /// Python writes `None` as an empty field: pass `""` for it.
    pub fn write_row<I>(&mut self, fields: I) -> io::Result<()>
    where
        I: IntoIterator,
        I::Item: AsRef<str>,
    {
        self.buf.clear();
        let mut count = 0;
        for field in fields {
            if count > 0 {
                self.buf.push(self.delimiter);
            }
            count += 1;
            let field = field.as_ref();
            let d = self.delimiter;
            if field.contains([d, QUOTE, '\n']) {
                self.buf.push(QUOTE);
                for c in field.chars() {
                    if c == QUOTE {
                        self.buf.push(QUOTE);
                    }
                    self.buf.push(c);
                }
                self.buf.push(QUOTE);
            } else {
                self.buf.push_str(field);
            }
        }
        if count == 1 && self.buf.is_empty() {
            self.buf.push_str("\"\"");
        }
        self.buf.push('\n');
        self.inner.write_all(self.buf.as_bytes())
    }

    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    pub fn into_inner(self) -> W {
        self.inner
    }
}

/// Why a [`Reader`] stopped. After an error the reader goes on with the next
/// line, as Python's does.
#[derive(Debug)]
pub enum Error {
    /// `csv.Error`: a `\r` or `\n` with more text after it on the same line
    /// outside quotes. Python's text mode (`newline=None`) never yields one.
    NewlineInUnquotedField { line: u64 },
    /// `csv.Error`: a field longer than [`FIELD_SIZE_LIMIT`] characters.
    FieldTooLarge { line: u64 },
    /// Reading a line failed.
    Io(io::Error),
}

impl Error {
    /// The Python exception class.
    pub fn exception(&self) -> &'static str {
        match self {
            Error::NewlineInUnquotedField { .. } | Error::FieldTooLarge { .. } => "Error",
            Error::Io(_) => "OSError",
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::NewlineInUnquotedField { line } => write!(
                f,
                "line {line}: new-line character seen in unquoted field - do you need to open the file with newline=''?"
            ),
            Error::FieldTooLarge { line } => {
                write!(
                    f,
                    "line {line}: field larger than field limit ({FIELD_SIZE_LIMIT})"
                )
            }
            Error::Io(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Error::Io(e)
    }
}

/// The states of CPython's `_csv.c` parser that the default dialect reaches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    StartRecord,
    StartField,
    InField,
    InQuotedField,
    QuoteInQuotedField,
    EatCrnl,
}

/// One character of a line, or the end of the line (`EOL` in `_csv.c`).
#[derive(Clone, Copy)]
enum Input {
    Char(char),
    Eol,
}

/// Python's `csv.reader(lines, delimiter=D, quoting=Q)`: an iterator of
/// rows over an iterator of lines that keep their endings.
pub struct Reader<I> {
    lines: I,
    delimiter: char,
    quoting: Quoting,
    state: State,
    field: String,
    field_chars: usize,
    fields: Vec<String>,
    line_num: u64,
}

impl<I, L> Reader<I>
where
    I: Iterator<Item = io::Result<L>>,
    L: AsRef<str>,
{
    pub fn new(lines: I, delimiter: char, quoting: Quoting) -> Self {
        Reader {
            lines,
            delimiter,
            quoting,
            state: State::StartRecord,
            field: String::new(),
            field_chars: 0,
            fields: Vec::new(),
            line_num: 0,
        }
    }

    /// `reader.line_num`: the lines read so far (a quoted field spanning
    /// lines counts each).
    pub fn line_num(&self) -> u64 {
        self.line_num
    }

    fn save_field(&mut self) {
        self.fields.push(std::mem::take(&mut self.field));
        self.field_chars = 0;
    }

    fn add_char(&mut self, c: char) -> Result<(), Error> {
        if self.field_chars >= FIELD_SIZE_LIMIT {
            return Err(Error::FieldTooLarge {
                line: self.line_num,
            });
        }
        self.field.push(c);
        self.field_chars += 1;
        Ok(())
    }

    /// `parse_process_char` of `_csv.c`, for the default dialect.
    fn process(&mut self, input: Input) -> Result<(), Error> {
        let quotes = self.quoting == Quoting::Minimal;
        let newline = matches!(input, Input::Char('\n' | '\r'));
        let eol = matches!(input, Input::Eol);
        // The state after a line ends inside a record.
        let ended = if eol {
            State::StartRecord
        } else {
            State::EatCrnl
        };
        match self.state {
            State::StartRecord if eol => {}
            State::StartRecord if newline => self.state = State::EatCrnl,
            State::StartRecord | State::StartField => {
                if newline || eol {
                    self.save_field();
                    self.state = ended;
                } else if let Input::Char(c) = input {
                    if c == QUOTE && quotes {
                        self.state = State::InQuotedField;
                    } else if c == self.delimiter {
                        self.save_field();
                        self.state = State::StartField;
                    } else {
                        self.add_char(c)?;
                        self.state = State::InField;
                    }
                }
            }
            State::InField => {
                if newline || eol {
                    self.save_field();
                    self.state = ended;
                } else if let Input::Char(c) = input {
                    if c == self.delimiter {
                        self.save_field();
                        self.state = State::StartField;
                    } else {
                        self.add_char(c)?;
                    }
                }
            }
            State::InQuotedField => {
                if let Input::Char(c) = input {
                    if c == QUOTE && quotes {
                        self.state = State::QuoteInQuotedField;
                    } else {
                        self.add_char(c)?;
                    }
                }
            }
            State::QuoteInQuotedField => {
                if newline || eol {
                    self.save_field();
                    self.state = ended;
                } else if let Input::Char(c) = input {
                    if c == QUOTE && quotes {
                        self.add_char(c)?;
                        self.state = State::InQuotedField;
                    } else if c == self.delimiter {
                        self.save_field();
                        self.state = State::StartField;
                    } else {
                        self.add_char(c)?;
                        self.state = State::InField;
                    }
                }
            }
            State::EatCrnl => {
                if eol {
                    self.state = State::StartRecord;
                } else if !newline {
                    return Err(Error::NewlineInUnquotedField {
                        line: self.line_num,
                    });
                }
            }
        }
        Ok(())
    }

    fn reset(&mut self) {
        self.state = State::StartRecord;
        self.field.clear();
        self.field_chars = 0;
        self.fields.clear();
    }
}

impl<I, L> Iterator for Reader<I>
where
    I: Iterator<Item = io::Result<L>>,
    L: AsRef<str>,
{
    type Item = Result<Vec<String>, Error>;

    /// The next row: `[]` for a blank line, as Python's reader returns.
    fn next(&mut self) -> Option<Self::Item> {
        self.reset();
        loop {
            let line = match self.lines.next() {
                Some(Ok(line)) => line,
                Some(Err(e)) => return Some(Err(Error::Io(e))),
                None => {
                    // The data ended inside a quoted field: strict=False
                    // keeps what was read.
                    if self.field_chars != 0 || self.state == State::InQuotedField {
                        self.save_field();
                        return Some(Ok(std::mem::take(&mut self.fields)));
                    }
                    return None;
                }
            };
            self.line_num += 1;
            for c in line.as_ref().chars() {
                if let Err(e) = self.process(Input::Char(c)) {
                    return Some(Err(e));
                }
            }
            if let Err(e) = self.process(Input::Eol) {
                return Some(Err(e));
            }
            if self.state == State::StartRecord {
                return Some(Ok(std::mem::take(&mut self.fields)));
            }
        }
    }
}

/// A header row as a dict sees it: each name once, in first-seen order,
/// with the column of its last occurrence.
#[derive(Debug)]
pub struct Header {
    names: Vec<String>,
    keys: Vec<(String, usize)>,
    index: HashMap<String, usize>,
}

impl Header {
    pub fn new(names: Vec<String>) -> Self {
        let mut keys: Vec<(String, usize)> = Vec::new();
        let mut index: HashMap<String, usize> = HashMap::new();
        for (col, name) in names.iter().enumerate() {
            match index.get(name.as_str()) {
                Some(&k) => keys[k].1 = col,
                None => {
                    index.insert(name.clone(), keys.len());
                    keys.push((name.clone(), col));
                }
            }
        }
        Header { names, keys, index }
    }

    /// The header cells as read, duplicates included (`fieldnames`).
    pub fn names(&self) -> &[String] {
        &self.names
    }
}

/// One row of a [`DictReader`]: the dict Python builds, without copying
/// the cells.
#[derive(Debug, Clone)]
pub struct DictRow {
    header: Arc<Header>,
    cells: Vec<String>,
}

impl DictRow {
    /// `row.get(key)` and `key in row` at once: `None` when the header has
    /// no such column; `Some(None)` when the row is too short to reach it
    /// (Python's `restval`); `Some(Some(cell))` otherwise.
    pub fn get(&self, key: &str) -> Option<Option<&str>> {
        let k = *self.header.index.get(key)?;
        Some(self.cells.get(self.header.keys[k].1).map(String::as_str))
    }

    /// The cells past the header's end, which Python keeps under the key
    /// `None` (`restkey`); `None` when there are none.
    pub fn extra(&self) -> Option<&[String]> {
        let n = self.header.names.len();
        (self.cells.len() > n).then(|| &self.cells[n..])
    }

    /// The dict's items in Python's order, without the `None` key.
    pub fn iter(&self) -> impl Iterator<Item = (&str, Option<&str>)> {
        self.header
            .keys
            .iter()
            .map(|(name, col)| (name.as_str(), self.cells.get(*col).map(String::as_str)))
    }

    /// The row's cells as read.
    pub fn cells(&self) -> &[String] {
        &self.cells
    }
}

/// Python's `csv.DictReader(lines, delimiter=D, quoting=Q)`: the first row
/// is the header; blank lines after it are skipped.
pub struct DictReader<I> {
    reader: Reader<I>,
    header: Option<Arc<Header>>,
}

impl<I, L> DictReader<I>
where
    I: Iterator<Item = io::Result<L>>,
    L: AsRef<str>,
{
    pub fn new(lines: I, delimiter: char, quoting: Quoting) -> Self {
        DictReader {
            reader: Reader::new(lines, delimiter, quoting),
            header: None,
        }
    }

    /// `fieldnames`: reads the header while there is none, as Python does,
    /// so after an error the next row is tried as the header. `None` for an
    /// empty input; empty for a blank first line.
    pub fn header(&mut self) -> Result<Option<&Header>, Error> {
        if self.header.is_none()
            && let Some(row) = self.reader.next()
        {
            self.header = Some(Arc::new(Header::new(row?)));
        }
        Ok(self.header.as_deref())
    }

    /// `line_num`: the lines read so far.
    pub fn line_num(&self) -> u64 {
        self.reader.line_num()
    }
}

impl<I, L> Iterator for DictReader<I>
where
    I: Iterator<Item = io::Result<L>>,
    L: AsRef<str>,
{
    type Item = Result<DictRow, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Err(e) = self.header() {
            return Some(Err(e));
        }
        let header = self.header.clone()?;
        loop {
            match self.reader.next()? {
                Ok(cells) if cells.is_empty() => continue,
                Ok(cells) => return Some(Ok(DictRow { header, cells })),
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lines(text: &[&str]) -> impl Iterator<Item = io::Result<String>> {
        text.iter()
            .map(|s| Ok(s.to_string()))
            .collect::<Vec<_>>()
            .into_iter()
    }

    fn rows(text: &[&str], quoting: Quoting) -> Vec<Result<Vec<String>, String>> {
        Reader::new(lines(text), '\t', quoting)
            .map(|r| r.map_err(|e| e.exception().to_owned()))
            .collect()
    }

    fn ok(cells: &[&str]) -> Result<Vec<String>, String> {
        Ok(cells.iter().map(|s| s.to_string()).collect())
    }

    fn written(delimiter: char, fields: &[&str]) -> String {
        let mut w = Writer::new(Vec::new(), delimiter);
        w.write_row(fields).unwrap();
        String::from_utf8(w.into_inner()).unwrap()
    }

    #[test]
    fn writer_quotes_minimally() {
        assert_eq!(written('\t', &["a", "b c", "d,e"]), "a\tb c\td,e\n");
        assert_eq!(written('\t', &["a\tb"]), "\"a\tb\"\n");
        assert_eq!(written('\t', &["x\"y"]), "\"x\"\"y\"\n");
        assert_eq!(written('\t', &["a\rb", "c\nd"]), "\"a\rb\"\t\"c\nd\"\n");
        assert_eq!(written(',', &["a,b", "a\tb"]), "\"a,b\",a\tb\n");
        assert_eq!(written('\t', &[""]), "\"\"\n");
        assert_eq!(written('\t', &["", ""]), "\t\n");
        assert_eq!(written('\t', &[] as &[&str]), "\n");
    }

    #[test]
    fn reader_follows_cpython() {
        let m = Quoting::Minimal;
        assert_eq!(rows(&["a\tb\n", "c\n"], m), [ok(&["a", "b"]), ok(&["c"])]);
        assert_eq!(rows(&["a\"b\"c\n"], m), [ok(&["a\"b\"c"])]);
        assert_eq!(rows(&["\"a\"b\n"], m), [ok(&["ab"])]);
        assert_eq!(rows(&["\"a\n", "b\"\n"], m), [ok(&["a\nb"])]);
        assert_eq!(rows(&["\"a"], m), [ok(&["a"])]);
        assert_eq!(rows(&["\"\"\n"], m), [ok(&[""])]);
        assert_eq!(rows(&["\"x\"\"y\"\tz\n"], m), [ok(&["x\"y", "z"])]);
        assert_eq!(rows(&["\n", "a\n\n"], m), [ok(&[]), ok(&["a"])]);
        assert_eq!(rows(&["a\r\n"], m), [ok(&["a"])]);
        assert_eq!(rows(&["a\tb\t\n"], m), [ok(&["a", "b", ""])]);
        assert_eq!(
            rows(&["a\rb\n", "c\n"], m),
            [Err("Error".to_owned()), ok(&["c"])]
        );
        assert_eq!(rows(&["\"a\"\tb\n"], Quoting::None), [ok(&["\"a\"", "b"])]);
    }

    #[test]
    fn reader_limits_field_size() {
        let long = "x".repeat(FIELD_SIZE_LIMIT);
        let first = rows(&[&format!("{long}\n")], Quoting::Minimal);
        assert_eq!(first, [Ok(vec![long.clone()])]);
        let over = rows(&[&format!("{long}x\n")], Quoting::Minimal);
        assert_eq!(over, [Err("Error".to_owned())]);
    }

    #[test]
    fn dict_reader_keeps_python_dict_rules() {
        let mut d = DictReader::new(
            lines(&["a\tb\ta\n", "1\t2\n", "1\t2\t3\t4\n", "\n", "5\n"]),
            '\t',
            Quoting::Minimal,
        );
        let names = d.header().unwrap().unwrap().names().to_vec();
        assert_eq!(names, ["a", "b", "a"]);
        let r: Vec<DictRow> = d.map(Result::unwrap).collect();
        assert_eq!(r.len(), 3);
        // The last "a" column decides, and a short row leaves it None.
        assert_eq!(
            r[0].iter().collect::<Vec<_>>(),
            [("a", None), ("b", Some("2"))]
        );
        assert_eq!(r[0].extra(), None);
        assert_eq!(r[1].get("a"), Some(Some("3")));
        assert_eq!(r[1].extra(), Some(&["4".to_owned()][..]));
        assert_eq!(r[2].get("b"), Some(None));
        assert_eq!(r[2].get("c"), None);
    }

    #[test]
    fn dict_reader_tries_the_next_row_after_a_bad_header() {
        let mut d = DictReader::new(
            lines(&["a\rb\n", "x\ty\n", "1\t2\n"]),
            '\t',
            Quoting::Minimal,
        );
        assert!(d.next().unwrap().is_err());
        let r: Vec<DictRow> = d.map(Result::unwrap).collect();
        assert_eq!(
            r[0].iter().collect::<Vec<_>>(),
            [("x", Some("1")), ("y", Some("2"))]
        );
    }

    #[test]
    fn dict_reader_blank_and_empty_inputs() {
        let mut empty = DictReader::new(lines(&[]), '\t', Quoting::Minimal);
        assert!(empty.header().unwrap().is_none());
        assert!(empty.next().is_none());
        let blank: Vec<DictRow> = DictReader::new(lines(&["\n", "1\t2\n"]), '\t', Quoting::Minimal)
            .map(Result::unwrap)
            .collect();
        assert_eq!(blank[0].iter().count(), 0);
        assert_eq!(blank[0].extra().map(<[String]>::len), Some(2));
    }
}
