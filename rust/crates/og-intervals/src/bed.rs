//! BED rows: `sequence start end [...]`, 0-based and end exclusive.

use og_io::{Lines, Utf8};
use pycompat::num::{parse_int, py_isspace};
use std::collections::HashMap;
use std::fmt;
use std::io::{self, BufRead};

/// The rows of a BED file, kept per sequence.
#[derive(Debug, Default)]
pub struct Rows {
    seqs: HashMap<String, Seq>,
    total: usize,
}

/// One sequence's rows: as read (for the counts), and merged.
#[derive(Debug, Default)]
pub(crate) struct Seq {
    /// `(start, end)` of every row, in file order.
    pub(crate) rows: Vec<(u64, u64)>,
    /// The rows sorted and joined where they overlap or touch: disjoint,
    /// ascending, none empty.
    pub(crate) merged: Vec<(u64, u64)>,
}

/// Why a BED file could not be read.
#[derive(Debug)]
pub enum BedError {
    Io(io::Error),
    /// A row the masking cannot apply: its line (1-based) and what is wrong.
    Row {
        line: usize,
        why: String,
    },
}

impl fmt::Display for BedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BedError::Io(e) => write!(f, "{e}"),
            BedError::Row { line, why } => write!(f, "line {line}: {why}"),
        }
    }
}

impl std::error::Error for BedError {}

impl From<io::Error> for BedError {
    fn from(e: io::Error) -> Self {
        BedError::Io(e)
    }
}

/// A coordinate field: Python's `int()` (signs, leading zeros and `_`
/// between digits read), and never negative.
fn coord(field: &str, what: &str, line: usize) -> Result<u64, BedError> {
    let bad = |why: String| BedError::Row { line, why };
    let v = parse_int(field).map_err(|_| bad(format!("{what} '{field}' is not a whole number")))?;
    u64::try_from(v).map_err(|_| bad(format!("{what} {v} is negative")))
}

impl Rows {
    /// Reads BED text: universal newlines, strict UTF-8. A line is a row
    /// when it has three or more whitespace-separated fields and does not
    /// start with `#`, `track` or `browser`; other lines are skipped. A row
    /// whose start or end is not a whole number, is negative, or whose end
    /// is before its start is an error. A row with end equal to start is
    /// kept and covers no base.
    pub fn read(input: impl BufRead) -> Result<Rows, BedError> {
        let mut lines = Lines::new(input, Utf8::Strict);
        let mut text = String::new();
        let mut rows = Rows::default();
        let mut n = 0;
        while lines.read_line(&mut text)? {
            n += 1;
            let mut fields = text.split(py_isspace).filter(|f| !f.is_empty());
            let (Some(seq), Some(start), Some(end)) = (fields.next(), fields.next(), fields.next())
            else {
                continue;
            };
            if ["#", "track", "browser"]
                .iter()
                .any(|p| text.starts_with(p))
            {
                continue;
            }
            let start = coord(start, "start", n)?;
            let end = coord(end, "end", n)?;
            if end < start {
                return Err(BedError::Row {
                    line: n,
                    why: format!("end {end} is before start {start}"),
                });
            }
            rows.seqs
                .entry(seq.to_owned())
                .or_default()
                .rows
                .push((start, end));
            rows.total += 1;
        }
        for seq in rows.seqs.values_mut() {
            seq.merged = merge(&seq.rows);
        }
        Ok(rows)
    }

    /// How many rows were read, on every sequence.
    pub fn total(&self) -> usize {
        self.total
    }

    /// The disjoint, ascending intervals of `seq`'s rows.
    pub fn merged(&self, seq: &str) -> &[(u64, u64)] {
        self.seqs.get(seq).map_or(&[], |s| &s.merged)
    }

    pub(crate) fn seq(&self, name: &str) -> Option<&Seq> {
        self.seqs.get(name)
    }

    /// `name` as this map's own key, when it has rows.
    pub(crate) fn key(&self, name: &str) -> Option<&str> {
        self.seqs.get_key_value(name).map(|(k, _)| k.as_str())
    }

    pub(crate) fn names(&self) -> impl Iterator<Item = &str> {
        self.seqs.keys().map(String::as_str)
    }
}

/// Sorts `rows` and joins the ones that overlap or touch; empty rows drop out.
pub fn merge(rows: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut sorted: Vec<(u64, u64)> = rows.iter().copied().filter(|(s, e)| s < e).collect();
    sorted.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(sorted.len());
    for (s, e) in sorted {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(text: &str) -> Result<Rows, BedError> {
        Rows::read(text.as_bytes())
    }

    #[test]
    fn skips_headers_and_short_lines() {
        let r = read("browser position s1:1-3\ntrack name=x\n#c\t1\t2\ns1\t5\n\ns1 5 9 name\n")
            .unwrap();
        assert_eq!(r.total(), 1);
        assert_eq!(r.merged("s1"), &[(5, 9)]);
    }

    #[test]
    fn universal_newlines_and_python_whitespace() {
        // A lone \r ends a line in Python's text mode; \x1f splits a field.
        let r = read("s1\t1\t3\rs2\x1f4\x1f6\r\ns1 +2_0 0030\n").unwrap();
        assert_eq!(r.total(), 3);
        assert_eq!(r.merged("s1"), &[(1, 3), (20, 30)]);
        assert_eq!(r.merged("s2"), &[(4, 6)]);
    }

    #[test]
    fn merges_overlapping_and_touching_rows() {
        assert_eq!(
            merge(&[(10, 20), (0, 5), (5, 7), (15, 30), (40, 40), (35, 36)]),
            vec![(0, 7), (10, 30), (35, 36)]
        );
    }

    #[test]
    fn refuses_bad_coordinates() {
        for (bed, want) in [
            ("s1\tx\t3\n", "line 1: start 'x' is not a whole number"),
            ("\ns1\t-1\t3\n", "line 2: start -1 is negative"),
            ("s1\t5\t3\n", "line 1: end 3 is before start 5"),
            ("s1\t5\t3.0\n", "line 1: end '3.0' is not a whole number"),
        ] {
            assert_eq!(read(bed).unwrap_err().to_string(), want, "{bed:?}");
        }
    }

    #[test]
    fn empty_bed() {
        let r = read("").unwrap();
        assert_eq!(r.total(), 0);
        assert!(r.merged("s1").is_empty());
    }
}
