//! The carrier deletion alleles and a male partner's AZF regions, from step
//! 35's interval table (`<id>_cn_intervals.tsv`), with step 18's CNVpytor
//! calls as their cross-check.

use crate::files::{each_line, each_record, info_get, int, is_digits};
use crate::findings::Finding;
use crate::list::CarrierList;
use og_io::Utf8;
use og_pedigree::Die;
use std::collections::HashMap;

/// The rows of `bin/cn_loci.tsv` the screen reads, with their names here.
const CN_ALLELES: [(&str, &str); 5] = [
    ("GJB6_D13S1830", "del(GJB6-D13S1830)"),
    ("CFTR_dele2_3", "CFTRdele2,3"),
    ("HBB_exon3_del", "HBB exon 3 deletion"),
    ("HBB_619bp_del", "HBB 619 bp deletion"),
    ("HBD_HBB_Sicilian_del", "Sicilian (delta-beta)0 deletion"),
];
/// The genes each allele is a change of, as the list names them
/// (del(GJB6-D13S1830) is a DFNB1 allele beside GJB2 changes).
const CN_ALLELE_GENES: [(&str, &[&str], &[&str]); 3] = [
    ("GJB2/GJB6 (DFNB1)", &["GJB2", "GJB6"], &["GJB6_D13S1830"]),
    ("CFTR", &["CFTR"], &["CFTR_dele2_3"]),
    (
        "HBB",
        &["HBB"],
        &["HBB_exon3_del", "HBB_619bp_del", "HBD_HBB_Sicilian_del"],
    ),
];
/// Named as deletions, so present reads as "the deletion is there".
const CN_AZF: [(&str, &str); 3] = [
    ("AZFa", "AZFa deletion"),
    ("AZFb", "AZFb deletion"),
    (
        "AZFc_DAZ",
        "AZFc deletion (DAZ copies; gr/gr and b2/b3 leave two of four)",
    ),
];
/// Step 18 calls on 1 kb bins, so a deletion shorter than two bins is one
/// it cannot see: its silence there is no evidence. A CNVpytor deletion
/// counts for a row when it covers half of the row's span or more.
const CNV_MIN_SPAN: i64 = 2000;

/// One row of the interval table, by column name.
pub type Row = HashMap<String, String>;
/// The interval table: rows by locus, a later row over an earlier one.
pub type Intervals = HashMap<String, Row>;
/// Step 18's deletions: (chrom, 0-based start, end).
pub type Deletions = Vec<(String, i64, i64)>;

/// Reads step 35's interval table.
pub fn read_intervals(path: &str) -> Result<Intervals, Die> {
    let mut head: Option<Vec<String>> = None;
    let mut out = Intervals::new();
    each_line(path, Utf8::Replace, |line| {
        if pycompat::text::strip(line, None).is_empty() {
            return Ok(());
        }
        let c: Vec<&str> = line.split('\t').collect();
        let Some(h) = &head else {
            let want = ["locus", "name", "chrom", "start", "end", "status"];
            let has = |k: &str| c.contains(&k);
            if c.len() < want.len()
                || c[..want.len()] != want
                || !["copies", "expected", "change", "reason"]
                    .iter()
                    .all(|k| has(k))
            {
                return Err(Die::one(format!(
                    "{path}: not step 35's interval table (og cn intervals)"
                )));
            }
            head = Some(c.iter().map(|s| (*s).to_owned()).collect());
            return Ok(());
        };
        if c.len() == h.len() {
            let row: Row = h
                .iter()
                .cloned()
                .zip(c.iter().map(|s| (*s).to_owned()))
                .collect();
            out.insert(c[0].to_owned(), row);
        }
        Ok(())
    })?;
    if head.is_none() {
        return Err(Die::one(format!(
            "{path}: not step 35's interval table (og cn intervals)"
        )));
    }
    Ok(out)
}

/// The deletions in step 18's CNVpytor VCF.
pub fn cnvpytor_deletions(path: &str) -> Result<Deletions, Die> {
    let mut out = Vec::new();
    each_record(path, |r| {
        if r.len() < 8 {
            return Ok(());
        }
        let kind = info_get(r[7], "SVTYPE")
            .filter(|k| !k.is_empty())
            .unwrap_or_else(|| r[4].trim_matches(['<', '>']));
        let end = info_get(r[7], "END").unwrap_or("");
        if kind != "DEL" || !is_digits(end) {
            return Ok(());
        }
        out.push((r[0].to_owned(), int(r[1], path)? - 1, int(end, path)?));
        Ok(())
    })?;
    Ok(out)
}

/// What a row says about its deletion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    Present,
    Absent,
    NotAssessable,
}

impl State {
    fn word(self) -> &'static str {
        match self {
            State::Present => "present",
            State::Absent => "absent",
            State::NotAssessable => "not assessable",
        }
    }
}

fn cell<'a>(row: &'a Row, k: &str) -> &'a str {
    row.get(k).map_or("", String::as_str)
}

/// The state of a row and its depth call in words.
fn state(row: Option<&Row>) -> (State, String) {
    let Some(row) = row else {
        return (
            State::NotAssessable,
            "no row in step 35's interval table".to_owned(),
        );
    };
    if cell(row, "status") != "called" {
        return (
            State::NotAssessable,
            format!("{} ({})", cell(row, "status"), cell(row, "reason")),
        );
    }
    let text = format!(
        "{} of {} copies",
        cell(row, "copies"),
        cell(row, "expected")
    );
    match cell(row, "change") {
        "deletion" => (State::Present, text),
        "no change" => (State::Absent, text),
        other => (State::NotAssessable, format!("{text}, {other}")),
    }
}

/// The CNVpytor cross-check of a row when it disagrees with the depth call.
fn cnv_note(row: Option<&Row>, dels: Option<&Deletions>) -> Result<&'static str, Die> {
    let (Some(row), Some(dels)) = (row, dels) else {
        return Ok("");
    };
    let what = "step 35's interval table";
    let (s, e) = (int(cell(row, "start"), what)?, int(cell(row, "end"), what)?);
    let chrom = cell(row, "chrom");
    let hit = dels
        .iter()
        .any(|(c, ds, de)| c == chrom && 2 * (e.min(*de) - s.max(*ds)) >= e - s);
    Ok(match state(Some(row)).0 {
        State::Present if !hit && e - s >= CNV_MIN_SPAN => {
            "CNVpytor (step 18) calls no deletion here"
        }
        State::Absent if hit => "CNVpytor (step 18) calls a deletion here",
        _ => "",
    })
}

/// The report lines of a partner's carrier deletion alleles and AZF regions.
pub fn lines(
    s: &str,
    sex: Option<&str>,
    intervals: &Intervals,
    dels: Option<&Deletions>,
) -> Result<Vec<String>, Die> {
    let mut out = Vec::new();
    let azf: &[(&str, &str)] = if sex == Some("male") { &CN_AZF } else { &[] };
    for (locus, name) in CN_ALLELES.iter().chain(azf) {
        let r = intervals.get(*locus);
        let (st, text) = state(r);
        let note = cnv_note(r, dels)?;
        let note = if note.is_empty() {
            String::new()
        } else {
            format!("; {note}")
        };
        out.push(format!(
            "- {s}, {name}: **{}** (read depth, step 35: {text}){note}",
            st.word()
        ));
    }
    if sex != Some("male") {
        out.push(format!(
            "- {s}, AZF regions: not applicable ({})",
            if sex == Some("female") {
                "no Y chromosome"
            } else {
                "sex unknown"
            }
        ));
    }
    Ok(out)
}

/// The genes where both partners carry a change once the read-depth alleles
/// count: one has a present allele and the other a present allele or a
/// finding in the same gene. The gene table leaves these alleles out, so
/// without this a couple who both carry CFTRdele2,3 would read as sharing
/// nothing. A gene the list leaves out is no shared risk here either.
pub fn shared(
    a: &str,
    b: &str,
    intervals: &HashMap<String, Intervals>,
    found: &HashMap<String, Vec<Finding>>,
    list: &CarrierList,
) -> Vec<String> {
    let names: HashMap<&str, &str> = CN_ALLELES.iter().copied().collect();
    let mut out = Vec::new();
    let empty = Intervals::new();
    for (label, genes, loci) in CN_ALLELE_GENES {
        if !genes.iter().any(|g| list.contains(g)) {
            continue;
        }
        let mut have: HashMap<&str, String> = HashMap::new();
        let mut by_depth = false;
        for s in [a, b] {
            let iv = intervals.get(s).unwrap_or(&empty);
            let mut parts: Vec<String> = loci
                .iter()
                .filter_map(|x| {
                    let (st, text) = state(iv.get(*x));
                    (st == State::Present)
                        .then(|| format!("{} present (read depth, step 35: {text})", names[x]))
                })
                .collect();
            by_depth |= !parts.is_empty();
            parts.extend(
                found
                    .get(s)
                    .into_iter()
                    .flatten()
                    .filter(|f| genes.contains(&f.gene.as_str()))
                    .map(|f| f.text.clone()),
            );
            if !parts.is_empty() {
                have.insert(s, parts.join("; "));
            }
        }
        if have.len() == 2 && by_depth {
            out.push(format!(
                "**{label}**: {a}: {}; {b}: {}. Both partners may carry a change in {label}: the gene table does not \
                 count read-depth alleles, so a clinical test has to confirm them before the chance for a child can \
                 be given.",
                have[a], have[b]
            ));
        }
    }
    out
}
