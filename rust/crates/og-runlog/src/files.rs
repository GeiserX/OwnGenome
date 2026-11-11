//! The run files: what run-all.sh, main.nf and the weblog write while a
//! run goes, and the one key of `summary.json` the run state reads.

use indexmap::IndexMap;
use og_io::{Lines, Utf8};
use pycompat::json::{Value, loads};
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

/// Each line of a text file, as Python's `for line in open(path,
/// errors="replace")` gives it, with its `\n` removed; `None` when the
/// file is not a regular file.
fn lines(path: &Path) -> io::Result<Option<Vec<String>>> {
    if !path.is_file() {
        return Ok(None);
    }
    let f = File::open(path)?;
    Lines::new(BufReader::new(f), Utf8::Replace)
        .collect::<io::Result<Vec<_>>>()
        .map(Some)
}

/// `logs/run_status.tsv`: `meta` rows (`started_epoch`, `finished_utc`,
/// `declared_sex`, `command`, ...) and the last row of each step, so
/// `planned` gives way to the step's result.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RunStatus {
    pub meta: IndexMap<String, String>,
    pub steps: IndexMap<String, String>,
}

impl RunStatus {
    /// `meta started_epoch` as a number; `None` when absent or not one.
    pub fn started(&self) -> Option<f64> {
        self.meta
            .get("started_epoch")
            .and_then(|v| pycompat::num::parse_float(v).ok())
    }
}

/// Reads `logs/run_status.tsv`; `None` when there is none.
pub fn read_run_status(path: &Path) -> io::Result<Option<RunStatus>> {
    let Some(lines) = lines(path)? else {
        return Ok(None);
    };
    let mut st = RunStatus::default();
    for line in &lines {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 3 || c[0].starts_with('#') {
            continue;
        }
        match c[0] {
            "meta" => {
                st.meta.insert(c[1].to_owned(), c[2].to_owned());
            }
            "step" => {
                st.steps.insert(c[1].to_owned(), c[2].to_owned());
            }
            _ => {}
        }
    }
    Ok(Some(st))
}

/// One process's row of main.nf's `nextflow/process_status.tsv`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ProcessCounts {
    pub succeeded: i128,
    pub cached: i128,
    pub failed: i128,
    pub ignored: i128,
    pub unfinished: i128,
}

impl ProcessCounts {
    pub fn any(&self) -> bool {
        [
            self.succeeded,
            self.cached,
            self.failed,
            self.ignored,
            self.unfinished,
        ]
        .iter()
        .any(|&n| n != 0)
    }
}

/// `{process name without the workflow prefix: counts}`; `None` without
/// the file. A row whose cells do not match the header, or whose counts are
/// not whole numbers, is left out.
pub fn read_process_status(path: &Path) -> io::Result<Option<IndexMap<String, ProcessCounts>>> {
    let Some(lines) = lines(path)? else {
        return Ok(None);
    };
    let mut it = lines.iter();
    let header: Vec<&str> = it.next().map_or(vec![""], |h| h.split('\t').collect());
    let mut out = IndexMap::new();
    for line in it {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() != header.len() {
            continue;
        }
        let mut r: IndexMap<&str, &str> = IndexMap::new();
        for (k, v) in header.iter().zip(&c) {
            r.insert(k, v);
        }
        let int = |k: &str| r.get(k).and_then(|v| pycompat::num::parse_int(v).ok());
        let (Some(succeeded), Some(cached), Some(failed), Some(ignored), Some(unfinished)) = (
            int("succeeded"),
            int("cached"),
            int("failed"),
            int("ignored"),
            int("unfinished"),
        ) else {
            continue;
        };
        let Some(process) = r.get("process") else {
            continue;
        };
        let name = process.rsplit(':').next().unwrap_or(process);
        out.insert(
            name.to_owned(),
            ProcessCounts {
                succeeded,
                cached,
                failed,
                ignored,
                unfinished,
            },
        );
    }
    Ok(Some(out))
}

/// The weblog messages of an events file: one JSON object per line; a
/// line that is not a JSON object is skipped. No file is no message.
pub fn read_events(path: Option<&Path>) -> io::Result<Vec<Value>> {
    let Some(lines) = (match path {
        Some(p) => lines(p)?,
        None => None,
    }) else {
        return Ok(Vec::new());
    };
    Ok(lines
        .iter()
        .map(|l| pycompat::text::strip(l, None))
        .filter(|l| !l.is_empty())
        .filter_map(|l| loads(l).ok())
        .filter(|v| matches!(v, Value::Object(_)))
        .collect())
}

/// `{section: (source path as summary.json gives it, its file time now)}`
/// from `SAMPLE_DIR/summary.json`; only each section's `source` is read.
/// No summary, or one that cannot be read, gives nothing.
pub fn summary_sources(sample_dir: &Path) -> IndexMap<String, (String, f64)> {
    let mut out = IndexMap::new();
    let Ok(text) = crate::util::read_text(&sample_dir.join("summary.json")) else {
        return out;
    };
    let Ok(Value::Object(top)) = loads(&text) else {
        return out;
    };
    let Some(Value::Object(secs)) = top.get("sections") else {
        return out;
    };
    for (k, v) in secs {
        let Value::Object(sec) = v else { continue };
        let Some(Value::Str(src)) = sec.get("source") else {
            continue;
        };
        if src.is_empty() {
            continue;
        }
        if let Ok(t) = pycompat::mtime::mtime(sample_dir.join(src)) {
            out.insert(k.clone(), (src.clone(), t));
        }
    }
    out
}

/// `run_manifest.tsv` as `{(section, key): value}`; nothing when there is
/// none. A value with tabs keeps them.
pub fn read_manifest(path: Option<&Path>) -> io::Result<IndexMap<(String, String), String>> {
    let mut rows = IndexMap::new();
    let Some(lines) = (match path {
        Some(p) => lines(p)?,
        None => None,
    }) else {
        return Ok(rows);
    };
    for line in &lines {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() >= 3 && !c[0].starts_with('#') {
            rows.insert((c[0].to_owned(), c[1].to_owned()), c[2..].join("\t"));
        }
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("og-files-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn run_status_last_row_wins() {
        let d = tmp("rs");
        let p = d.join("run_status.tsv");
        std::fs::write(
            &p,
            "# x\nmeta\tstarted_epoch\t1767344400\nstep\t11\tplanned\nstep\t11\tok\nstep\t14\tskipped (no BAM)\nshort\tx\n",
        )
        .unwrap();
        let st = read_run_status(&p).unwrap().unwrap();
        assert_eq!(st.started(), Some(1_767_344_400.0));
        assert_eq!(st.steps["11"], "ok");
        assert_eq!(st.steps["14"], "skipped (no BAM)");
        assert!(read_run_status(&d.join("none")).unwrap().is_none());
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn process_status_rows() {
        let d = tmp("ps");
        let p = d.join("ps.tsv");
        std::fs::write(
            &p,
            "process\tsucceeded\tcached\tfailed\tignored\tunfinished\nW:A\t1\t0\t0\t0\t0\nB\tx\t0\t0\t0\t0\nC\t1\t0\n",
        )
        .unwrap();
        let ps = read_process_status(&p).unwrap().unwrap();
        assert_eq!(ps.len(), 1);
        assert_eq!(ps["A"].succeeded, 1);
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn events_skip_bad_lines() {
        let d = tmp("ev");
        let p = d.join("e.jsonl");
        std::fs::write(
            &p,
            "{\"event\":\"started\"}\nnot json\n[1]\n\n  {\"event\":\"completed\"}  \n",
        )
        .unwrap();
        assert_eq!(read_events(Some(&p)).unwrap().len(), 2);
        assert!(read_events(None).unwrap().is_empty());
        std::fs::remove_dir_all(&d).unwrap();
    }
}
