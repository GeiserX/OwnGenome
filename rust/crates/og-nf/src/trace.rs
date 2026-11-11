//! A Nextflow trace file: a header line of column names, then one row per
//! task, tab-separated. Read as text with universal newlines and invalid
//! UTF-8 replaced, one row in memory at a time.

use og_io::{Lines, Utf8};
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

/// One row, its cells keyed by the header. A row shorter than the header
/// lacks the missing columns; cells past the header are dropped; a column
/// named twice gives its last cell.
pub struct Row<'a> {
    header: &'a [String],
    cells: Vec<&'a str>,
}

impl<'a> Row<'a> {
    /// The cell of column `name`, or `None` when the row does not reach it.
    pub fn get(&self, name: &str) -> Option<&'a str> {
        let n = self.header.len().min(self.cells.len());
        (0..n)
            .rev()
            .find(|&i| self.header[i] == name)
            .map(|i| self.cells[i])
    }

    /// True when the row holds only spaces.
    pub fn is_blank(&self) -> bool {
        self.cells
            .iter()
            .all(|c| c.chars().all(pycompat::text::is_space))
    }
}

/// Calls `f` on every row of the trace at `path`, in file order. An error
/// opening or reading the file is returned; a malformed row is still a row.
pub fn for_each_row(path: &Path, mut f: impl FnMut(&Row<'_>)) -> io::Result<()> {
    let mut lines = Lines::new(BufReader::new(File::open(path)?), Utf8::Replace);
    let mut line = String::new();
    if !lines.read_line(&mut line)? {
        return Ok(());
    }
    let header: Vec<String> = line
        .strip_suffix('\n')
        .unwrap_or(&line)
        .split('\t')
        .map(str::to_owned)
        .collect();
    while lines.read_line(&mut line)? {
        let body = line.strip_suffix('\n').unwrap_or(&line);
        let row = Row {
            header: &header,
            cells: body.split('\t').collect(),
        };
        f(&row);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_by_header() {
        let dir = std::env::temp_dir().join(format!("og-nf-trace-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("trace_x.txt");
        std::fs::write(&p, "a\tb\ta\r\n1\t2\t3\t4\r\n5\n\n").unwrap();
        let mut got = Vec::new();
        for_each_row(&p, |r| {
            got.push((
                r.get("a").map(str::to_owned),
                r.get("b").map(str::to_owned),
                r.is_blank(),
            ))
        })
        .unwrap();
        assert_eq!(
            got,
            vec![
                (Some("3".into()), Some("2".into()), false),
                (Some("5".into()), None, false),
                (Some(String::new()), None, true)
            ]
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
