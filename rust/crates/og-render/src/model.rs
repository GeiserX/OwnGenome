//! The two models the report is made of: a typed view of summary.json, read
//! the way `bin/render_report.py` reads it (Python's truthiness, `str()` and
//! `==`), and the document the section modules build from it, which the
//! text and HTML writers lay out.
//!
//! A section module never writes markup. It returns cards of typed lines;
//! every text in them is escaped once, by the HTML writer.

use std::fmt;

use indexmap::IndexMap;
use pycompat::json::Value;
use pycompat::value;

/// A summary that does not have the shape the report reads: the path of the
/// first value that is missing or of the wrong kind.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Malformed(pub String);

impl fmt::Display for Malformed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Malformed {}

/// `str(v)`, which is also `f"{v}"`.
pub fn py(v: &Value) -> Result<String, Malformed> {
    value::py_str(v).map_err(|e| Malformed(e.to_string()))
}

/// `bool(v)`.
pub fn truthy(v: &Value) -> bool {
    value::truthy(v)
}

/// The number of `v` when `isinstance(v, (int, float))`, which holds for a
/// bool too, as in Python.
pub fn num(v: Option<&Value>) -> Option<f64> {
    match v? {
        Value::Int(i) => Some(*i as f64),
        Value::Float(x) => Some(*x),
        Value::Bool(b) => Some(f64::from(u8::from(*b))),
        _ => None,
    }
}

/// `a == b` as Python compares two JSON values: numbers by value across int,
/// float and bool, everything else by kind and content.
pub fn py_eq(a: &Value, b: &Value) -> bool {
    match (num(Some(a)), num(Some(b))) {
        (Some(x), Some(y)) => x == y,
        _ => a == b,
    }
}

/// `v == s` for a string `s`.
pub fn is_str(v: Option<&Value>, s: &str) -> bool {
    matches!(v, Some(Value::Str(x)) if x == s)
}

/// `key in v`: an item of a list, a key of a dict, a substring of a string.
pub fn contains(v: &Value, key: &str) -> bool {
    match v {
        Value::Array(a) => a.iter().any(|x| is_str(Some(x), key)),
        Value::Object(o) => o.contains_key(key),
        Value::Str(s) => s.contains(key),
        _ => false,
    }
}

/// `d.get(key)` on a value that must be a dict; `path` names it in the
/// error.
pub fn get<'a>(d: &'a Value, key: &str, path: &str) -> Result<Option<&'a Value>, Malformed> {
    value::get(d, key).map_err(|_| Malformed(format!("{path} is not an object")))
}

/// `d[key]`: the value, or an error naming `path.key`.
pub fn need<'a>(d: &'a Value, key: &str, path: &str) -> Result<&'a Value, Malformed> {
    get(d, key, path)?.ok_or_else(|| Malformed(format!("{path} has no {key}")))
}

/// A section's state in summary.json.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    Ok,
    Stale,
    Failed,
    Unreadable,
    Missing,
    /// A state no collector writes; the report treats it as Python does, as
    /// none of the five.
    Other,
}

impl State {
    fn parse(v: Option<&Value>) -> State {
        match v {
            Some(Value::Str(s)) => match s.as_str() {
                "ok" => State::Ok,
                "stale" => State::Stale,
                "failed" => State::Failed,
                "unreadable" => State::Unreadable,
                "missing" => State::Missing,
                _ => State::Other,
            },
            _ => State::Missing,
        }
    }

    /// The step ran and its result was read: this run's (`ok`) or an
    /// earlier run's (`stale`).
    pub fn ran(self) -> bool {
        matches!(self, State::Ok | State::Stale)
    }
}

/// One section of the summary, or the empty stand-in of a key the summary
/// does not have (Python's `MISSING`: no title, state missing, no data).
#[derive(Debug, Clone, Copy)]
pub struct Section<'a> {
    pub key: &'a str,
    pub title: &'a Value,
    pub state: State,
    pub note: &'a Value,
    pub data: &'a Value,
    pub file_date: Option<&'a Value>,
    pub source: Option<&'a Value>,
}

impl<'a> Section<'a> {
    /// `data.get(key)`.
    pub fn get(&self, key: &str) -> Result<Option<&'a Value>, Malformed> {
        get(self.data, key, &format!("sections.{}.data", self.key))
    }

    /// `data[key]`.
    pub fn need(&self, key: &str) -> Result<&'a Value, Malformed> {
        need(self.data, key, &format!("sections.{}.data", self.key))
    }
}

/// summary.json, with the values Python's `MISSING` stand-in points at.
#[derive(Debug)]
pub struct Summary {
    root: Value,
    empty: Value,
    null: Value,
    blank: Value,
}

impl Summary {
    /// Reads summary.json's text. The top level must be an object holding
    /// `sections`; the frame's other keys are checked where they are read.
    pub fn parse(text: &str) -> Result<Summary, Malformed> {
        let root = pycompat::json::loads(text).map_err(|e| Malformed(e.to_string()))?;
        need(&root, "sections", "summary")?;
        Ok(Summary {
            root,
            empty: Value::Object(IndexMap::new()),
            null: Value::Null,
            blank: Value::Str(String::new()),
        })
    }

    /// `s[key]`.
    pub fn need(&self, key: &str) -> Result<&Value, Malformed> {
        need(&self.root, key, "summary")
    }

    /// `s.get(key)`.
    pub fn get(&self, key: &str) -> Option<&Value> {
        value::get(&self.root, key).ok().flatten()
    }

    /// `s[a][b]`.
    pub fn need2(&self, a: &str, b: &str) -> Result<&Value, Malformed> {
        need(self.need(a)?, b, a)
    }

    /// `S.get(key, MISSING)`; a section object without a title, note or
    /// data reads as empty there.
    pub fn section<'a>(&'a self, key: &'a str) -> Section<'a> {
        let sec = self
            .get("sections")
            .and_then(|s| value::get(s, key).ok().flatten())
            .filter(|s| matches!(s, Value::Object(_)));
        let field = |k: &str| sec.and_then(|s| value::get(s, k).ok().flatten());
        Section {
            key,
            title: field("title").unwrap_or(&self.blank),
            state: sec.map_or(State::Missing, |_| State::parse(field("state"))),
            note: field("note").unwrap_or(&self.null),
            data: field("data").unwrap_or(&self.empty),
            file_date: field("file_date"),
            source: field("source"),
        }
    }

    /// Every section of the summary, in its order.
    pub fn sections(&self) -> Vec<Section<'_>> {
        match self.get("sections") {
            Some(Value::Object(o)) => o.keys().map(|k| self.section(k)).collect(),
            _ => Vec::new(),
        }
    }

    /// `(s.get("new_since") or {})`: when the summary compares this run with
    /// a previous one.
    pub fn new_since(&self) -> NewSince<'_> {
        let ns = self.get("new_since").filter(|v| truthy(v));
        let field = |k: &str| ns.and_then(|n| value::get(n, k).ok().flatten());
        NewSince {
            previous: field("previous_utc").filter(|v| truthy(v)),
            sections: field("sections").filter(|v| truthy(v)),
        }
    }

    /// The rows of `run_manifest.tsv`, each `[kind, key, value]`; empty
    /// when the summary has none.
    pub fn manifest(&self) -> Result<Vec<[&Value; 3]>, Malformed> {
        let m = self.need("manifest")?;
        if !truthy(m) {
            return Ok(Vec::new());
        }
        let Value::Array(rows) = m else {
            return Err(Malformed("summary.manifest is not a list".into()));
        };
        rows.iter()
            .enumerate()
            .map(|(i, r)| match r {
                Value::Array(c) if c.len() >= 3 => Ok([&c[0], &c[1], &c[2]]),
                _ => Err(Malformed(format!(
                    "summary.manifest[{i}] is not a row of three"
                ))),
            })
            .collect()
    }

    /// A list of texts at the top level (`not_run`, `not_assessed`).
    pub fn texts(&self, key: &str) -> Result<Vec<String>, Malformed> {
        match self.need(key)? {
            Value::Array(a) => a.iter().map(py).collect(),
            Value::Null => Ok(Vec::new()),
            _ => Err(Malformed(format!("summary.{key} is not a list"))),
        }
    }
}

/// What changed since the previous run, as the summary records it.
#[derive(Debug, Clone, Copy)]
pub struct NewSince<'a> {
    previous: Option<&'a Value>,
    sections: Option<&'a Value>,
}

impl NewSince<'_> {
    /// `ns.get("previous_utc") or "the previous run"`.
    pub fn previous(&self) -> Result<String, Malformed> {
        self.previous.map_or(Ok("the previous run".to_owned()), py)
    }

    /// `key in (ns.get("sections") or [])`, the text report's test.
    pub fn lists(&self, key: &str) -> bool {
        self.sections.is_some_and(|s| contains(s, key))
    }

    /// `set(ns.get("sections") or []) & set(keys)` is not empty, the HTML
    /// report's test: a string counts by its characters.
    pub fn any_of(&self, keys: &[&str]) -> bool {
        match self.sections {
            Some(Value::Str(s)) => keys
                .iter()
                .any(|k| k.chars().count() == 1 && s.contains(*k)),
            Some(v) => keys.iter().any(|k| contains(v, k)),
            None => false,
        }
    }
}

// --- the document ------------------------------------------------------------------

/// A badge's colour class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Colour {
    Green,
    Yellow,
    Red,
    Gray,
    State,
}

impl Colour {
    pub fn class(self) -> &'static str {
        match self {
            Colour::Green => "green",
            Colour::Yellow => "yellow",
            Colour::Red => "red",
            Colour::Gray => "gray",
            Colour::State => "state",
        }
    }
}

/// The value side of a stat line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Inline {
    /// Plain text, in the stat's bold weight.
    Text(String),
    /// A coloured badge.
    Badge(String, Colour),
    /// Text in normal weight and a smaller size, for a reason or a remark.
    Quiet(String),
}

/// One line of a card's body.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Line {
    /// A label and its value.
    Stat(String, Inline),
    /// A sentence under the numbers that is not a caveat.
    Note(String),
    /// A paragraph.
    Para(String),
    /// A plain bulleted list.
    List(Vec<String>),
    /// A table of plain cells.
    Table(Table),
    /// A warning in the stale box's style, among the card's lines.
    Warn(String),
    /// The sentence that says how the Viewer column opens the reads.
    ViewerNote,
}

/// A table: its column titles and rows of cells, each cell escaped once.
/// With `viewer`, a last column opens each row's locus in the genome viewer
/// (`chrom:start` or `chrom:start-end`); a row with `None` gets an empty
/// cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Table {
    pub head: Vec<String>,
    pub rows: Vec<Vec<String>>,
    pub viewer: Option<Vec<Option<String>>>,
    pub layout: Layout,
}

/// How a table is laid out in the page's source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layout {
    /// The whole table on one line.
    OneLine,
    /// The table tags and the head row on lines of their own, the rows on
    /// one line each with no indent.
    RowLines,
}

/// One card: a title, the caveats shown under it, the notes saying its
/// result is not this run's own, and its lines.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Card {
    pub title: String,
    /// Spans the whole width of the chapter's grid.
    pub full: bool,
    /// Caveat sentences (`*_NOTE`), in the order the card shows them.
    pub caveats: Vec<&'static str>,
    /// Stale, failed or unreadable notes, in a box under the caveats.
    pub stale: Vec<String>,
    pub body: Vec<Line>,
}

/// A section's part of the HTML report: its cards and the links to the
/// step's own files beside the report.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct HtmlPart {
    pub cards: Vec<Card>,
    pub files: Vec<(String, String)>,
}
