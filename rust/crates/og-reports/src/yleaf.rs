//! Step 37's two steps around Yleaf's own command line.
//!
//! ```text
//! og y yleaf-positions --data DATA OUT
//! og y yleaf-no-reads --bam BAM --idxstats IDXSTATS OUT
//! ```
//!
//! `yleaf-positions` writes `chrY<TAB>pos` for each distinct position in
//! column 4 of `DATA/hg38/new_positions.txt`, in file order, for
//! `samtools mpileup -l`, and checks that
//! `DATA/hg_prediction_tables/tree.json` is there: the two files Yleaf 4.1's
//! `yleaf_constants` names (`HG38`, `NEW_POSITION_FILE`, `TREE_FILE`).
//! Step 37 then hands that pileup to Yleaf (`Yleaf -rp`).
//!
//! `yleaf-no-reads` writes Yleaf's prediction table for a pileup with no
//! line: Yleaf stops on an empty pileup ("No columns to parse") and writes
//! a header with no row, so the step writes the row Yleaf gives when too
//! few markers had reads: `Hg` NA, `Valid_markers` 0, and `Total_reads`,
//! the mapped reads in `samtools idxstats` (column 3, summed). The sample
//! name is the BAM's file name without its last extension, as Yleaf names
//! it.
//!
//! What changes on purpose from the Python these replace: the output is
//! written whole or not at all, and a position or a read count is ASCII
//! digits only (the Python's `isdigit` also took other Unicode digits).

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

/// The header of Yleaf 4's `hg_prediction.hg`.
pub const HG_HEADER: &str =
    "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n";

/// `Path(bam).name.rsplit('.', 1)[0]`: Yleaf's name for the sample.
fn sample_name(bam: &str) -> &str {
    let name = bam.trim_end_matches('/').rsplit('/').next().unwrap_or("");
    name.rsplit_once('.').map_or(name, |(stem, _)| stem)
}

fn no_reads(bam: &str, idxstats: &str, out: &str) -> Outcome {
    let fail = |e| (e, 1);
    let mut reads: u128 = 0;
    for line in crate::lines(idxstats, Utf8::Strict).map_err(fail)? {
        let line = line.map_err(fail)?;
        if let Some(c) = line.split('\t').nth(2) {
            let c = crate::strip(c);
            if crate::digits(c) {
                let n = crate::int(c).and_then(|n| u128::try_from(n).ok());
                reads = n.and_then(|n| reads.checked_add(n)).ok_or_else(|| {
                    fail(Error(format!("read count {c} too large in {idxstats}")))
                })?;
            }
        }
    }
    let body = format!(
        "{HG_HEADER}{}\tNA\t\t{reads}\t0\tNA\tNA\tNA\tNA\n",
        sample_name(bam)
    );
    crate::write_atomic(out, body.as_bytes()).map_err(fail)?;
    Ok("No read at any Y marker: Hg NA (insufficient markers)\n".to_owned())
}

const NO_READS_SPEC: Spec = Spec {
    prog: "og y yleaf-no-reads",
    usage: "usage: og y yleaf-no-reads [-h] --bam BAM --idxstats IDXSTATS out\n",
    options: &["--bam", "--idxstats"],
    positionals: &["out"],
    required: &["--bam", "--idxstats"],
};

pub(crate) fn no_reads_cli(rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    crate::run(&NO_READS_SPEC, rest, out, err, |a| {
        Ok(no_reads(
            crate::opt(a, "--bam"),
            crate::opt(a, "--idxstats"),
            &a.pos[0],
        ))
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

    #[test]
    fn no_reads_writes_the_na_row_with_the_mapped_reads() {
        let d = Dir::new("yleaf-noreads");
        // idxstats: name, length, mapped, unmapped; the last line is `*`.
        let idx = d.put(
            "idx.txt",
            "chr1\t248956422\t35933\t331\nchrY\t57227415\t1200\t7\nshort\t1\n*\t0\t0\t2048\n",
        );
        let out = d.path("hg_prediction.hg");
        let line = no_reads("/x/aligned/S1_sorted.bam", &idx, &out).unwrap();
        assert_eq!(
            line,
            "No read at any Y marker: Hg NA (insufficient markers)\n"
        );
        assert_eq!(
            d.read("hg_prediction.hg"),
            format!("{HG_HEADER}S1_sorted\tNA\t\t37133\t0\tNA\tNA\tNA\tNA\n")
        );
        assert_eq!(sample_name("a.b.bam"), "a.b");
        assert_eq!(sample_name("noext"), "noext");
        assert_eq!(sample_name("dir/x.bam"), "x");
        let e = no_reads("x.bam", &d.path("missing.txt"), &d.path("o")).unwrap_err();
        assert_eq!(e.1, 1);
        assert!(std::fs::metadata(d.path("o")).is_err());
    }

    #[test]
    fn no_reads_cli_routes_through_og_y() {
        let d = Dir::new("yleaf-noreads-cli");
        let idx = d.put("idx.txt", "chrY\t57227415\t5\t0\n");
        let out = d.path("t.hg");
        let args: Vec<String> = ["yleaf-no-reads", "--bam", "s.bam", "--idxstats", &idx, &out]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let (mut o, mut e) = (Vec::new(), Vec::new());
        let code = crate::y::cli(Some("yleaf-no-reads"), &args[1..], &mut o, &mut e);
        assert_eq!(code, 0, "{}", String::from_utf8_lossy(&e));
        assert_eq!(
            d.read("t.hg"),
            format!("{HG_HEADER}s\tNA\t\t5\t0\tNA\tNA\tNA\tNA\n")
        );
    }
}
