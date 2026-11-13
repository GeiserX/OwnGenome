//! Which pinned databases, bundles and images have fallen behind upstream.
//!
//! Each lookup asks one upstream and returns a typed answer ([`Upstream`])
//! beside the row of the Markdown report ([`DbRow`]). `og dev freshness`
//! renders the report from the rows; `og recheck` reads the typed answers,
//! never the report's wording. Every lookup fails closed: an HTTP error, a
//! timeout, an empty answer or an answer without the expected field is a
//! [`LookupError`], never "current".

pub mod cli;
pub mod http;
pub mod issue;
pub mod lookups;
pub mod pins;
pub mod registry;
pub mod report;
pub mod selftest;
pub mod version;

use pycompat::json::Value;
use std::path::{Path, PathBuf};

/// A lookup that gave no usable answer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LookupError {
    /// The answer failed or said nothing: an HTTP error, no answer, an
    /// empty value.
    Failed(String),
    /// The answer came in a shape the lookup does not know (a list where an
    /// object was expected, a missing key).
    Shape(String),
}

impl std::fmt::Display for LookupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LookupError::Failed(m) => f.write_str(m),
            LookupError::Shape(m) => write!(f, "unexpected answer ({m})"),
        }
    }
}

impl std::error::Error for LookupError {}

pub type Lookup<T> = Result<T, LookupError>;

pub fn failed<T>(msg: impl Into<String>) -> Lookup<T> {
    Err(LookupError::Failed(msg.into()))
}

/// `obj[key]` of a JSON object; anything else is a shape error.
pub fn field<'a>(v: &'a Value, key: &str) -> Lookup<&'a Value> {
    match v {
        Value::Object(m) => m
            .get(key)
            .ok_or_else(|| LookupError::Shape(format!("no key {key:?}"))),
        _ => Err(LookupError::Shape(format!(
            "looked up {key:?} in something that is not an object"
        ))),
    }
}

/// `obj.get(key)` of a JSON object: `None` when it is missing or null.
pub fn opt_field<'a>(v: &'a Value, key: &str) -> Lookup<Option<&'a Value>> {
    match v {
        Value::Object(m) => Ok(m.get(key).filter(|x| !matches!(x, Value::Null))),
        _ => Err(LookupError::Shape(format!(
            "looked up {key:?} in something that is not an object"
        ))),
    }
}

/// A non-empty string value at `key`, or a failed lookup naming `what`
/// (an empty value is an error, never "up to date").
pub fn need_str(v: &Value, key: &str, what: &str) -> Lookup<String> {
    match opt_field(v, key)? {
        None => failed(format!("{what}: empty")),
        Some(Value::Str(s)) if s.is_empty() => failed(format!("{what}: empty")),
        Some(Value::Str(s)) => Ok(s.clone()),
        Some(_) => Err(LookupError::Shape(format!("{what} is not a string"))),
    }
}

/// The elements of a JSON array, or a shape error.
pub fn array<'a>(v: &'a Value, what: &str) -> Lookup<&'a [Value]> {
    match v {
        Value::Array(a) => Ok(a),
        _ => Err(LookupError::Shape(format!("{what} is not a list"))),
    }
}

/// A string value, or a shape error.
pub fn string<'a>(v: &'a Value, what: &str) -> Lookup<&'a str> {
    match v {
        Value::Str(s) => Ok(s),
        _ => Err(LookupError::Shape(format!("{what} is not a string"))),
    }
}

/// One row of the databases table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DbRow {
    pub item: String,
    pub pinned: String,
    pub upstream: String,
    pub state: String,
    pub todo: String,
}

/// What upstream has, as a value a program can compare.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Upstream {
    /// A file's release date, `YYYYMMDD` (ClinVar).
    Date(String),
    /// The newest release, as upstream names it (`117`, `20260620`, `3.4.0`).
    Release(String),
    /// pypgx: the newest version with both an image and a bundle tag that
    /// the held pin could move to, if any.
    CanMoveTo(Option<String>),
    /// A row that reports and compares nothing (gnomAD, held).
    ReportOnly,
}

/// A lookup's answer: its report row and its typed value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Found {
    pub row: DbRow,
    pub upstream: Upstream,
}

/// The checkout an `og` binary belongs to: the first folder above the
/// binary, else above the working directory, that holds `versions.env`
/// and `assets/steps.tsv`.
pub fn checkout_root() -> Option<PathBuf> {
    let is_root =
        |d: &Path| d.join("versions.env").is_file() && d.join("assets/steps.tsv").is_file();
    let from = |p: PathBuf| p.ancestors().find(|d| is_root(d)).map(Path::to_path_buf);
    std::env::current_exe()
        .ok()
        .and_then(from)
        .or_else(|| std::env::current_dir().ok().and_then(from))
}
