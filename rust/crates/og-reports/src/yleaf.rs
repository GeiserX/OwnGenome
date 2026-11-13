//! `og y yleaf-positions` (step 37): Yleaf's own GRCh38 marker positions,
//! for the pileup step 37 makes beside Yleaf.
//!
//! ```text
//! og y yleaf-positions --data DATA OUT
//! ```
//!
//! Writes `chrY<TAB>pos` for each distinct position in column 4 of
//! `DATA/hg38/new_positions.txt`, in file order, for `samtools mpileup -l`,
//! and checks that `DATA/hg_prediction_tables/tree.json` is there: the two
//! files Yleaf 4.1's `yleaf_constants` names (`HG38`, `NEW_POSITION_FILE`,
//! `TREE_FILE`). Running Yleaf on that pileup stays in `bin/yleaf_run.py
//! predict`, which is glue inside Yleaf's own Python.
//!
//! What changes on purpose: the output is written whole or not at all, and
//! a position is ASCII digits only (the Python's `isdigit` also took other
//! Unicode digits).

use crate::y::{is_file, path_join};
use crate::{Error, Outcome};
use og_filters::args::Spec;
use og_io::Utf8;
use std::collections::HashSet;
use std::fmt::Write as _;
use std::io::Write;

/// Yleaf 4.1's `HG38 / NEW_POSITION_FILE`.
pub const POSITIONS: &str = "hg38/new_positions.txt";
/// Yleaf 4.1's `HG_PREDICTION_FOLDER / TREE_FILE`, relative to its data folder.
pub const TREE: &str = "hg_prediction_tables/tree.json";

/// `Path(data).absolute()`: the working directory joined to a relative path.
fn absolute(data: &str) -> Result<String, Error> {
    if data.starts_with('/') {
        return Ok(data.to_owned());
    }
    let cwd = std::env::current_dir()
        .map_err(|e| Error(format!("cannot read the working directory: {e}")))?;
    let cwd = cwd.to_string_lossy();
    Ok(if data.is_empty() {
        cwd.into_owned()
    } else {
        format!("{cwd}/{data}")
    })
}

fn positions(data: &str, out: &str) -> Outcome {
    let fail = |e| (e, 1);
    let data = absolute(data).map_err(fail)?;
    let src = path_join(&data, POSITIONS);
    for need in [src.clone(), path_join(&data, TREE)] {
        if !is_file(&need) {
            return Err(fail(Error(format!(
                "{need} not found: install Yleaf's data with scripts/setup.sh --yleaf-data <genome_dir>"
            ))));
        }
    }
    let mut seen: HashSet<String> = HashSet::new();
    let mut body = String::new();
    for line in crate::lines(&src, Utf8::Strict).map_err(fail)? {
        let line = line.map_err(fail)?;
        let mut cols = line.strip_suffix('\n').unwrap_or(&line).split('\t');
        // Python checked len(cols) > 3, then the digits and the set.
        if let Some(p) = cols.nth(3)
            && crate::digits(p)
            && !seen.contains(p)
        {
            let _ = writeln!(body, "chrY\t{p}");
            seen.insert(p.to_owned());
        }
    }
    if seen.is_empty() {
        return Err(fail(Error(format!("no marker positions in {src}"))));
    }
    crate::write_atomic(out, body.as_bytes()).map_err(fail)?;
    Ok(format!("{} Y marker positions from {src}\n", seen.len()))
}

const SPEC: Spec = Spec {
    prog: "og y yleaf-positions",
    usage: "usage: og y yleaf-positions [-h] --data DATA out\n",
    options: &["--data"],
    positionals: &["out"],
    required: &["--data"],
};

pub(crate) fn cli(rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    crate::run(&SPEC, rest, out, err, |a| {
        Ok(positions(crate::opt(a, "--data"), &a.pos[0]))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    #[test]
    fn distinct_positions_in_file_order() {
        let d = Dir::new("yleaf-pos");
        d.put(
            POSITIONS,
            "chry\ta\tR\t300\tC->T\tC\tT\nchry\tb\tR\t100\tC->T\tC\tT\nchry\tc\tR\t300\tA->G\tA\tG\n\
             chry\td\tR\t0100\tA->G\tA\tG\nshort\tline\nchry\te\tR\tx1\tA\tA\tG\n",
        );
        d.put(TREE, "{}");
        let data = d.root().display().to_string();
        let line = positions(&data, &d.path("out.txt")).unwrap();
        assert_eq!(d.read("out.txt"), "chrY\t300\nchrY\t100\nchrY\t0100\n");
        assert_eq!(
            line,
            format!("3 Y marker positions from {data}/{POSITIONS}\n")
        );
        std::fs::remove_file(d.path(TREE)).unwrap();
        let e = positions(&data, &d.path("out2.txt")).unwrap_err().0;
        assert_eq!(
            e.0,
            format!(
                "{data}/{TREE} not found: install Yleaf's data with scripts/setup.sh --yleaf-data <genome_dir>"
            )
        );
    }
}
