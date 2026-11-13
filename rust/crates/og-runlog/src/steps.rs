//! `assets/steps.tsv`: one row per step of a run.

use indexmap::IndexMap;
use og_io::{Lines, Utf8};
use pycompat::num::parse_float;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

/// One row of the catalogue.
#[derive(Debug, Clone, PartialEq)]
pub struct Step {
    pub step: String,
    pub label: String,
    /// The main.nf `--tools` name; empty when the step has none.
    pub tool: String,
    /// Nextflow processes, by name without the workflow prefix.
    pub processes: Vec<String>,
    /// The summary.json sections it feeds.
    pub sections: Vec<String>,
    /// The steps whose output it reads.
    pub upstream: Vec<String>,
    /// The upper end of its usual wall time, in hours.
    pub hours: Option<f64>,
    /// `measured`, `estimate` or empty.
    pub hours_from: String,
    /// The measured lower end, as the file spells it (empty when unknown).
    pub hours_low: String,
    /// The `run_manifest.tsv` keys of the databases and pins it reads
    /// itself: `KEY`, a prefix `PREFIX*`, either with `@self` when only a
    /// table of its own reads it (empty when the file has no `reads` column).
    pub reads: Vec<String>,
}

/// Why the catalogue could not be read.
#[derive(Debug)]
pub enum StepsError {
    Io(PathBuf, std::io::Error),
    /// A column the catalogue must have is not in its header.
    NoColumn(PathBuf, &'static str),
    /// `hours` is not a number: the file, the line and the cell.
    BadHours(PathBuf, u64, String),
}

impl std::fmt::Display for StepsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StepsError::Io(p, e) => write!(f, "{}: {e}", p.display()),
            StepsError::NoColumn(p, c) => write!(f, "{}: no column '{c}'", p.display()),
            StepsError::BadHours(p, n, v) => {
                write!(f, "{}:{n}: hours is not a number: {v:?}", p.display())
            }
        }
    }
}

impl std::error::Error for StepsError {}

/// The catalogue next to the running binary: the first folder above it that
/// holds `assets/steps.tsv` (`bin/.og/og-linux-*` and
/// `rust/target/*/og` both sit inside a checkout).
pub fn default_path() -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    exe.ancestors()
        .skip(1)
        .map(|d| d.join("assets").join("steps.tsv"))
        .find(|p| p.is_file())
}

/// Reads the catalogue, in file order. Lines starting with `#` and blank
/// lines are skipped; the first other line is the header; a short row has
/// empty cells for the missing columns.
pub fn read_steps(path: &Path) -> Result<Vec<Step>, StepsError> {
    let io = |e| StepsError::Io(path.to_path_buf(), e);
    let file = File::open(path).map_err(io)?;
    let mut header: Option<Vec<String>> = None;
    let mut out = Vec::new();
    for (i, line) in Lines::new(BufReader::new(file), Utf8::Strict).enumerate() {
        let line = line.map_err(io)?;
        let n = i as u64 + 1;
        if line.starts_with('#') || line.chars().all(pycompat::text::is_space) {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let Some(h) = &header else {
            header = Some(cols.iter().map(|c| (*c).to_owned()).collect());
            continue;
        };
        let mut row: IndexMap<&str, &str> = IndexMap::new();
        for (i, name) in h.iter().enumerate() {
            row.insert(name.as_str(), cols.get(i).copied().unwrap_or(""));
        }
        let col = |c: &'static str| {
            row.get(c)
                .copied()
                .ok_or_else(|| StepsError::NoColumn(path.to_path_buf(), c))
        };
        let split = |c: &'static str| -> Result<Vec<String>, StepsError> {
            Ok(col(c)?
                .split(',')
                .filter(|x| !x.is_empty())
                .map(str::to_owned)
                .collect())
        };
        let hours = match col("hours")? {
            "" => None,
            v => Some(
                parse_float(v)
                    .map_err(|_| StepsError::BadHours(path.to_path_buf(), n, v.to_owned()))?,
            ),
        };
        out.push(Step {
            step: col("step")?.to_owned(),
            label: col("label")?.to_owned(),
            tool: col("tool")?.to_owned(),
            processes: split("processes")?,
            sections: split("sections")?,
            upstream: split("upstream")?,
            hours,
            hours_from: col("hours_from")?.to_owned(),
            hours_low: row.get("hours_low").copied().unwrap_or("").to_owned(),
            reads: row
                .get("reads")
                .copied()
                .unwrap_or("")
                .split(',')
                .filter(|x| !x.is_empty())
                .map(str::to_owned)
                .collect(),
        });
    }
    Ok(out)
}

/// `{process: step}`; a process two steps list belongs to the first.
pub fn process_steps(steps: &[Step]) -> IndexMap<&str, &str> {
    let mut out = IndexMap::new();
    for s in steps {
        for p in &s.processes {
            out.entry(p.as_str()).or_insert(s.step.as_str());
        }
    }
    out
}

/// Seconds of silence after which a run with no end record reads
/// `stopped`: twice the longest step's usual hours, never under 6 h. A run
/// goes quiet while one long task runs, so less would call a healthy run
/// stopped.
pub fn stale_after(steps: &[Step]) -> f64 {
    const FLOOR: f64 = 6.0 * 3600.0;
    crate::util::py_max(
        std::iter::once(FLOOR).chain(
            steps
                .iter()
                .filter_map(|s| s.hours.filter(|h| *h != 0.0).map(|h| 2.0 * h * 3600.0)),
        ),
    )
    .unwrap_or(FLOOR)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn repo_steps() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets/steps.tsv")
    }

    #[test]
    fn reads_the_catalogue() {
        let steps = read_steps(&repo_steps()).unwrap();
        assert!(steps.len() > 50);
        let s02 = steps.iter().find(|s| s.step == "02").unwrap();
        assert_eq!(
            s02.processes,
            ["MINIMAP2_INDEX", "ALIGN_MINIMAP2", "CRAM_TO_BAM"]
        );
        assert_eq!(s02.upstream, ["01b"]);
        assert_eq!(s02.hours, Some(5.4));
        assert_eq!(s02.hours_from, "measured");
        // PRS is listed by steps 25 and 26: the first owns it.
        assert_eq!(process_steps(&steps)["PRS"], "25");
        // DeepVariant's 13.5 h makes it 27 h.
        assert_eq!(stale_after(&steps), 27.0 * 3600.0);
    }

    #[test]
    fn refuses_a_bad_hours_cell() {
        let dir = std::env::temp_dir().join(format!("og-steps-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("steps.tsv");
        std::fs::write(&p, "# c\nstep\tlabel\ttool\tprocesses\tsections\tupstream\thours\thours_from\n01\tx\t\tA\t\t\tlong\t\n").unwrap();
        let e = read_steps(&p).unwrap_err().to_string();
        assert!(e.ends_with(":3: hours is not a number: \"long\""), "{e}");
        std::fs::write(&p, "step\tlabel\n01\tx\n").unwrap();
        assert!(
            read_steps(&p)
                .unwrap_err()
                .to_string()
                .contains("no column")
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
