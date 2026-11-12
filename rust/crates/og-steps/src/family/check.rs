//! `og family check`: for every pair of sequenced members, the relatedness
//! the pedigree implies against the one somalier relate measured. Exits 1
//! when a pair disagrees: a sample swap, or a wrong line in the pedigree.

use std::collections::{BTreeSet, HashMap};
use std::path::Path;

use og_pedigree::{Die, Kinship, cell, declared_relationship, degree, read_members, read_table};
use pycompat::fmt::fixed;

use super::args::{Args, Kind, Spec};
use super::{Result, print_out, write_file};

pub const SPEC: Spec = Spec {
    prog: "og family check",
    usage: "usage: og family check [-h] --members MEMBERS --pairs PAIRS --out OUT\n",
    options: &[
        ("--members", Kind::One),
        ("--pairs", Kind::One),
        ("--out", Kind::One),
    ],
    required: &["--members", "--pairs", "--out"],
};

/// somalier relate needs this many sites both samples have a genotype at
/// before its relatedness is read as a verdict.
const MIN_SITES: i128 = 1000;
/// A parent and child share an allele at every site: the share of sites
/// where they are opposite homozygotes (IBS0) stays near 0.
const PARENT_CHILD_MAX_IBS0: f64 = 0.01;

/// (verdict, reason) for one pair.
fn verdict(
    declared: &str,
    expected: f64,
    rel: f64,
    ibs0_rate: f64,
    n: i128,
) -> (&'static str, String) {
    if n < MIN_SITES {
        return (
            "not_checked",
            format!(
                "somalier compared the two at {n} sites, fewer than {MIN_SITES}: too few to tell \
                 (a sliced or targeted sample, or the wrong sites file)"
            ),
        );
    }
    let (want, got) = (degree(expected), degree(rel));
    if expected >= 0.177 {
        if got != want {
            return (
                "mismatch",
                format!("declared {declared} ({want}), but somalier measures {got}"),
            );
        }
        if declared == "parent-child" && ibs0_rate > PARENT_CHILD_MAX_IBS0 {
            return (
                "mismatch",
                format!(
                    "first degree, but opposite homozygotes at {}% of sites: a parent and \
                     child share an allele at every site (siblings do not)",
                    fixed(100.0 * ibs0_rate, 1)
                ),
            );
        }
        return ("ok", format!("{got}, as declared"));
    }
    if rel >= 0.354 {
        return (
            "mismatch",
            format!("declared {declared}, but somalier measures {got}"),
        );
    }
    if got != want && rel > expected {
        return (
            "warn",
            format!("declared {declared}, somalier measures {got}: closer than the pedigree says"),
        );
    }
    ("ok", format!("{got}, as declared"))
}

pub fn run(a: &Args) -> Result<()> {
    let members_path = a.one("--members").unwrap_or_default();
    let pairs_path = a.one("--pairs").unwrap_or_default();
    let members = read_members(members_path)?;
    let mut pairs: HashMap<BTreeSet<String>, indexmap::IndexMap<String, String>> = HashMap::new();
    for r in read_table(pairs_path, true)? {
        let key: BTreeSet<String> = [
            cell(&r, "sample_a", pairs_path)?.to_owned(),
            cell(&r, "sample_b", pairs_path)?.to_owned(),
        ]
        .into_iter()
        .collect();
        pairs.insert(key, r);
    }
    let mut kin = Kinship::new(&members);
    let mut seq: Vec<&String> = members
        .iter()
        .filter(|(_, m)| m.sequenced)
        .map(|(s, _)| s)
        .collect();
    seq.sort();
    let (mut out, mut bad, mut errors) = (Vec::new(), Vec::new(), Vec::new());
    let mut warn = String::new();
    let num_err = |v: &str, k: &str| Die::one(format!("{pairs_path}: {k} '{v}' is not a number"));
    for (i, x) in seq.iter().enumerate() {
        for y in &seq[i + 1..] {
            let key: BTreeSet<String> = [(*x).clone(), (*y).clone()].into_iter().collect();
            let Some(r) = pairs.get(&key) else {
                let all: BTreeSet<&String> = pairs.keys().flatten().collect();
                let all: Vec<&str> = all.into_iter().map(String::as_str).collect();
                errors.push(format!(
                    "somalier relate has no row for {x} and {y}; its samples are: {}",
                    all.join(", ")
                ));
                continue;
            };
            let e = 2.0 * kin.phi(x, y);
            let d = declared_relationship(x, y, &members, e);
            let float = |k: &str| -> std::result::Result<f64, Die> {
                let v = cell(r, k, pairs_path)?;
                pycompat::num::parse_float(v).map_err(|_| num_err(v, k))
            };
            let int = |k: &str| -> std::result::Result<i128, Die> {
                let v = cell(r, k, pairs_path)?;
                pycompat::num::parse_int(v).map_err(|_| num_err(v, k))
            };
            let (rel, ibs0, ibs2, n) =
                (float("relatedness")?, int("ibs0")?, int("ibs2")?, int("n")?);
            let rate = if n != 0 { ibs0 as f64 / n as f64 } else { 0.0 };
            let (v, why) = verdict(&d, e, rel, rate, n);
            if v == "mismatch" {
                bad.push(format!("{x} and {y}: {why}"));
            } else if v != "ok" {
                warn.push_str(&format!("WARNING: {x} and {y}: {why}\n"));
            }
            out.push([
                (*x).clone(),
                (*y).clone(),
                d,
                fixed(e, 3),
                fixed(rel, 3),
                ibs0.to_string(),
                fixed(rate, 4),
                ibs2.to_string(),
                n.to_string(),
                v.to_owned(),
                why,
            ]);
        }
    }
    if !warn.is_empty() {
        eprint!("{warn}");
    }
    if !errors.is_empty() {
        return Err(Die::new(errors).into());
    }
    let mut text = String::from(
        "sample_a\tsample_b\tdeclared\texpected_relatedness\trelatedness\tibs0\tibs0_rate\tibs2\tsites\tverdict\treason\n",
    );
    for row in &out {
        text.push_str(&row.join("\t"));
        text.push('\n');
    }
    write_file(
        Path::new(a.one("--out").unwrap_or_default()),
        text.as_bytes(),
    )?;
    let mut say = String::new();
    for row in &out {
        say.push_str(&format!(
            "{} and {}: {}; relatedness {} (expected {}), {} sites: {}\n",
            row[0], row[1], row[2], row[4], row[3], row[8], row[9]
        ));
    }
    print_out(&say);
    if !bad.is_empty() {
        let mut lines = vec!["the data disagree with the pedigree:".to_owned()];
        lines.extend(bad.iter().map(|b| format!("  {b}")));
        lines.push(
            "A sample swap (at the lab or in the samplesheet), or a wrong line in the pedigree. \
             Check the pedigree and the samplesheet; somalier's own tables are beside relationships.tsv."
                .to_owned(),
        );
        return Err(Die { lines, code: 1 }.into());
    }
    Ok(())
}
