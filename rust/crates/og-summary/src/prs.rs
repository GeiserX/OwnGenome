//! `og prs table`: step 25's PRS table and step 26's ancestry table, from
//! pgsc_calc's output folder (Lambert 2024).
//!
//! pgsc_calc writes a match summary (`<set>/match/<set>_summary.csv`), the
//! scores (`<set>/score/<set>_pgs.txt.gz` with an ancestry panel, else
//! `aggregated_scores.txt.gz`), the panel placement
//! (`<set>_popsimilarity.txt.gz`) and match logs (`*_log.csv.gz`). The
//! scoring files `og prs format` wrote give each score's total, label, chrX
//! rows, sex, strand-ambiguous rows and PGS Catalog evaluations.
//!
//! - A score for the other sex keeps its sum and loses its percentile: a
//!   percentile for an organ the person does not have reads as a risk.
//! - The effect is shown only for the group the sample was placed in, and
//!   only when the score was evaluated there.
//! - With a panel, chrX rows the panel run did not match are reported as
//!   left out.
//! - pgsc_calc's zero-match and below-threshold exits become rows, not
//!   errors.
//!
//! The table's columns are [`COLUMNS`], with [`EFFECT_COLUMNS`] before
//! `Input` when any score has evaluations. Stdout has one line per score,
//! which the e2e cases read.

use crate::args::{self, Kind, Opt};
use crate::{basename, finish, float, is_file, join, text_lines, trunc, write_file};
use og_steps::prs::{CHRX_FILE, EVALS_FILE, Eval, Labels, SCORE_NOTES_FILE};
use pycompat::csv::{DictReader, DictRow, Quoting, Reader};
use pycompat::fmt::{fixed, general};
use pycompat::text::{is_space, isdigit};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::io::{self, Write};
use std::path::Path;

/// The table's columns; Input stays last, readers take it as the last one.
pub const COLUMNS: [&str; 13] = [
    "Condition",
    "PGS_ID",
    "Score_SUM",
    "Variants_Matched",
    "Variants_Total",
    "Matched_Pct",
    "Percentile",
    "Ancestry_Group",
    "ChrX",
    "Ambiguous",
    "Applies_To",
    "Sex_Check",
    "Input",
];
/// Written before Input when the score list names PGS Catalog evaluations.
pub const EFFECT_COLUMNS: [&str; 2] = ["Evaluated_In", "Effect"];

/// Sex_Check of a sex-specific score.
pub const SEX_APPLIES: &str = "applies";
pub const SEX_NOT_APPLICABLE: &str = "not applicable to the declared sex";
pub const SEX_NOT_DECLARED: &str = "no sex declared";

const USAGE: &str = "usage: og prs table --sample S --results DIR --sampleset NAME --scores DIR --input-kind {gvcf,vcf}
                    [--panel NAME] [--zero-matches] [--below-threshold LOG] [--sex {,male,female}]
                    --out F [--ancestry-out F]
";

const OPTS: &[Opt] = &[
    Opt("--sample", Kind::Value, true),
    Opt("--results", Kind::Value, true),
    Opt("--sampleset", Kind::Value, true),
    Opt("--scores", Kind::Value, true),
    Opt("--input-kind", Kind::Choice(&["gvcf", "vcf"]), true),
    Opt("--panel", Kind::Value, false),
    Opt("--zero-matches", Kind::Flag, false),
    Opt("--below-threshold", Kind::Value, false),
    Opt("--sex", Kind::Choice(&["", "male", "female"]), false),
    Opt("--out", Kind::Value, true),
    Opt("--ancestry-out", Kind::Value, false),
];

/// How pgsc_calc ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome<'a> {
    /// It scored: read its match summary and scores.
    Scored,
    /// It stopped because no score variant is in the sample's genotypes.
    ZeroMatches,
    /// It stopped because every score matched under its minimum overlap;
    /// the console log at this path names each rate.
    BelowThreshold(&'a str),
}

/// What the table is made from.
#[derive(Debug)]
pub struct Inputs<'a> {
    pub results: &'a str,
    pub sampleset: &'a str,
    pub scores: &'a str,
    pub input_kind: &'a str,
    pub outcome: Outcome<'a>,
    /// The declared sex, `""` for none.
    pub sex: &'a str,
    /// pgsc_calc ran with an ancestry panel.
    pub panel: bool,
}

/// One row of the table: the 15 cells of [`COLUMNS`] with
/// [`EFFECT_COLUMNS`] before Input.
pub type TableRow = [String; 15];

/// The ancestry rows of step 26: (key, value).
pub type Ancestry = Vec<(String, String)>;

/// `PGS000018` for `PGS000018.txt.gz`, `.../PGS000018_hmPOS_GRCh38.txt.gz`.
fn pgs_id_of(path: &str) -> String {
    og_steps::prs::pgs_id_of(Path::new(basename(path)))
}

/// `d.get(k)` of a csv row as an f-string prints it: a short row's missing
/// cell is `None`.
fn shown(v: Option<Option<&str>>, default: &str) -> String {
    match v {
        None => default.to_owned(),
        Some(None) => "None".to_owned(),
        Some(Some(s)) => s.to_owned(),
    }
}

/// `r.get(k, "")` used as text: a short row's missing cell reads as empty.
fn get<'a>(r: &'a DictRow, k: &str) -> &'a str {
    r.get(k).flatten().unwrap_or("")
}

/// The rows of a tab- or comma-separated file with a header, as
/// `csv.DictReader` reads them, one at a time.
fn each_row(
    path: &str,
    delimiter: char,
    mut f: impl FnMut(DictRow) -> Result<(), String>,
) -> Result<(), String> {
    let mut lines = text_lines(path)?;
    let iter = std::iter::from_fn(move || {
        let mut s = String::new();
        match lines.read_line(&mut s) {
            Ok(true) => Some(Ok(s)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    });
    for r in DictReader::new(iter, delimiter, Quoting::Minimal) {
        f(r.map_err(|e| e.to_string())?)?;
    }
    Ok(())
}

/// `f"{float(x):.10g}"`, the number as pgsc_calc wrote it without a
/// trailing `.0`; `NA` for no number.
fn fmt_num(x: Option<&str>) -> String {
    match x.map(float) {
        Some(Ok(v)) => general(v, 10),
        _ => "NA".to_owned(),
    }
}

/// The formatted scoring files: PGS id to (rows, label), sorted by id.
fn totals(scores: &str) -> Result<BTreeMap<String, (i64, Option<String>)>, String> {
    let mut files: Vec<String> = fs::read_dir(scores)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| !n.starts_with('.') && n.ends_with(".txt.gz"))
        .collect();
    files.sort();
    let mut out = BTreeMap::new();
    for name in files {
        let p = join(scores, &name);
        let (mut n, mut label) = (0i64, None);
        let mut lines = text_lines(&p)?;
        let mut line = String::new();
        while lines
            .read_line(&mut line)
            .map_err(|e| format!("{p}: {e}"))?
        {
            if let Some(rest) = line.strip_prefix("#trait_reported=") {
                label = Some(rest.strip_suffix('\n').unwrap_or(rest).to_owned());
            } else if !line.starts_with('#') {
                n += 1;
            }
        }
        // minus the column header
        // A second file of one id keeps the first one's label when it has none.
        let pid = pgs_id_of(&name);
        let label = label.or_else(|| {
            out.get(&pid)
                .and_then(|(_, l): &(i64, Option<String>)| l.clone())
        });
        out.insert(pid, (n - 1, label));
    }
    if out.is_empty() {
        return Err(format!("no formatted scoring file in {scores}"));
    }
    Ok(out)
}

/// `Score X fails minimum matching threshold (N% variants match)` of
/// pgscatalog-match's log: PGS id to the rate with one decimal.
fn below_rates(text: &str) -> Result<HashMap<String, String>, String> {
    const HEAD: &str = "Score ";
    const MID: &str = " fails minimum matching threshold (";
    const TAIL: &str = "% variants match)";
    let mut out = HashMap::new();
    let mut from = 0;
    while let Some(k) = text[from..].find(HEAD) {
        let start = from + k + HEAD.len();
        let id_len: usize = text[start..]
            .chars()
            .take_while(|&c| !is_space(c))
            .map(char::len_utf8)
            .sum();
        let after = start + id_len;
        let matched = (id_len > 0)
            .then(|| text[after..].strip_prefix(MID))
            .flatten()
            .and_then(|rest| {
                let num_len = rest
                    .bytes()
                    .take_while(|b| b.is_ascii_digit() || *b == b'.')
                    .count();
                (num_len > 0)
                    .then(|| rest[num_len..].strip_prefix(TAIL))
                    .flatten()
                    .map(|_| {
                        (
                            &rest[..num_len],
                            text.len() - rest.len() + num_len + TAIL.len(),
                        )
                    })
            });
        match matched {
            Some((rate, end)) => {
                out.insert(pgs_id_of(&text[start..after]), fixed(float(rate)?, 1));
                from = end;
            }
            None => from += k + 1,
        }
    }
    Ok(out)
}

/// PGS id to the chrX rows pgsc_calc's match logs matched; nothing when
/// there is no log or a log lacks the columns.
fn chrx_matched(results: &str, sampleset: &str) -> Result<Option<HashMap<String, i64>>, String> {
    let dir = join(&join(results, sampleset), "match");
    let mut logs: Vec<String> = fs::read_dir(&dir)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| !n.starts_with('.') && n.ends_with("_log.csv.gz"))
        .map(|n| format!("{dir}/{n}"))
        .collect();
    logs.sort();
    if logs.is_empty() {
        return Ok(None);
    }
    let one_row = |line: &str| -> Result<Vec<String>, String> {
        Reader::new(
            std::iter::once(Ok::<_, io::Error>(line)),
            ',',
            Quoting::Minimal,
        )
        .next()
        .transpose()
        .map(Option::unwrap_or_default)
        .map_err(|e| e.to_string())
    };
    let mut out = HashMap::new();
    for p in logs {
        let mut lines = text_lines(&p)?;
        let mut line = String::new();
        lines
            .read_line(&mut line)
            .map_err(|e| format!("{p}: {e}"))?;
        let head = one_row(&line)?;
        let col = |k: &str| head.iter().rposition(|h| h == k);
        let (Some(acc), Some(chr), Some(status)) =
            (col("accession"), col("chr_name"), col("match_status"))
        else {
            return Ok(None);
        };
        // Python's `len(r) > max(c.values())`: every column of the header.
        let need = head.len() - 1;
        while lines
            .read_line(&mut line)
            .map_err(|e| format!("{p}: {e}"))?
        {
            // Most lines are not on chrX; split only those that may be.
            if !line.contains(",X,") && !line.contains(",23,") {
                continue;
            }
            let r = one_row(&line)?;
            if r.len() > need && (r[chr] == "X" || r[chr] == "23") && r[status] == "matched" {
                *out.entry(pgs_id_of(&r[acc])).or_insert(0) += 1;
            }
        }
    }
    Ok(Some(out))
}

/// The table rows and the ancestry rows (key, value).
pub fn table(i: &Inputs) -> Result<(Vec<TableRow>, Ancestry), String> {
    let totals = totals(i.scores)?;
    let side = |name: &str| join(i.scores, name);
    // "8 scored", "8 left out" or "0"; none for scores formatted before the file existed
    let (mut chrx, mut chrx_scored) = (HashMap::new(), HashSet::new());
    if is_file(side(CHRX_FILE)) {
        each_row(&side(CHRX_FILE), '\t', |r| {
            let n = r.get("chrx_rows").unwrap_or(Some(""));
            let pid = get(&r, "pgs_id").to_owned();
            let counted = n != Some("") && n != Some("0");
            let v = if counted {
                format!("{} {}", shown(Some(n), ""), shown(r.get("chrx"), ""))
            } else {
                n.unwrap_or("").to_owned()
            };
            if counted && r.get("chrx").flatten() == Some("scored") {
                chrx_scored.insert(pid.clone());
            }
            chrx.insert(pid, v);
            Ok(())
        })?;
    }
    let (mut applies, mut ambig) = (HashMap::new(), HashMap::new());
    if is_file(side(SCORE_NOTES_FILE)) {
        each_row(&side(SCORE_NOTES_FILE), '\t', |r| {
            let pid = get(&r, "pgs_id").to_owned();
            let n = r.get("ambiguous_rows").unwrap_or(Some(""));
            applies.insert(pid.clone(), get(&r, "applies_to").to_owned());
            let v = if n != Some("") && n != Some("0") {
                format!("{} {}", shown(Some(n), ""), shown(r.get("ambiguous"), ""))
            } else {
                n.unwrap_or("").to_owned()
            };
            ambig.insert(pid, v);
            Ok(())
        })?;
    }
    // The evaluations og prs format copied from assets/pgs_scores.tsv.
    let evals = if is_file(side(EVALS_FILE)) {
        Labels::read(Path::new(&side(EVALS_FILE)))?
    } else {
        Labels::default()
    };
    let sex = pycompat::text::lower(i.sex);
    let (mut matched, mut sums, mut pct, mut group) = (
        HashMap::new(),
        HashMap::new(),
        HashMap::new(),
        HashMap::new(),
    );
    let mut pops: Option<DictRow> = None;
    let mut rate = HashMap::new();
    match i.outcome {
        Outcome::BelowThreshold(log) => rate = below_rates(&crate::read_text(log)?)?,
        Outcome::ZeroMatches => {}
        Outcome::Scored => {
            let (res, ss) = (i.results, i.sampleset);
            let summ = join(&join(&join(res, ss), "match"), &format!("{ss}_summary.csv"));
            if !is_file(&summ) {
                return Err(format!("pgsc_calc wrote no match summary ({summ})"));
            }
            each_row(&summ, ',', |r| {
                let pid = pgs_id_of(get(&r, "accession"));
                if r.get("match_status").flatten() == Some("matched") {
                    let c = match r.get("count").flatten().filter(|c| !c.is_empty()) {
                        Some(c) => trunc(float(c)?)?,
                        None => 0,
                    };
                    *matched.entry(pid).or_insert(0i128) += c;
                }
                Ok(())
            })?;
            let score_dir = join(&join(res, ss), "score");
            let adjusted = format!("{score_dir}/{ss}_pgs.txt.gz");
            let plain = format!("{score_dir}/aggregated_scores.txt.gz");
            let src = if is_file(&adjusted) {
                adjusted
            } else {
                plain.clone()
            };
            if !is_file(&src) {
                return Err(format!("pgsc_calc wrote no scores ({plain})"));
            }
            each_row(&src, '\t', |r| {
                if r.get("sampleset").flatten() != Some(ss) {
                    return Ok(()); // the reference panel's own samples
                }
                let pid = pgs_id_of(get(&r, "PGS"));
                sums.insert(pid.clone(), fmt_num(r.get("SUM").flatten()));
                if let Some(p) = r
                    .get("percentile_MostSimilarPop")
                    .flatten()
                    .filter(|p| !p.is_empty())
                {
                    let v = fmt_num(Some(p));
                    let v = if v == "NA" { v } else { fixed(float(&v)?, 1) };
                    pct.insert(pid, v);
                }
                Ok(())
            })?;
            let popsim = format!("{score_dir}/{ss}_popsimilarity.txt.gz");
            if is_file(&popsim) {
                let mut target = Vec::new();
                each_row(&popsim, '\t', |r| {
                    if r.get("sampleset").flatten() == Some(ss)
                        && r.get("REFERENCE").unwrap_or(Some("False")) != Some("True")
                    {
                        target.push(r);
                    }
                    Ok(())
                })?;
                if target.len() != 1 {
                    return Err(format!(
                        "{} has {} rows for {ss}, expected 1",
                        basename(&popsim),
                        target.len()
                    ));
                }
                let p = target.pop().unwrap_or_else(|| unreachable!());
                let g = p
                    .get("MostSimilarPop")
                    .flatten()
                    .filter(|g| !g.is_empty())
                    .unwrap_or("NA");
                for pid in sums.keys() {
                    group.insert(pid.clone(), g.to_owned());
                }
                pops = Some(p);
            }
            if i.panel && !chrx_scored.is_empty() {
                // chrX rows given to pgsc_calc that its matches did not keep.
                if let Some(xm) = chrx_matched(res, ss)? {
                    for pid in &chrx_scored {
                        let text: &String = &chrx[pid];
                        let first = text.split_whitespace().next().unwrap_or("");
                        let n_x = crate::int(first)?;
                        let got = i128::from(xm.get(pid).copied().unwrap_or(0));
                        if got == 0 {
                            chrx.insert(
                                pid.clone(),
                                format!("{n_x} left out (not matched in the panel run)"),
                            );
                        } else if got < n_x {
                            chrx.insert(
                                pid.clone(),
                                format!(
                                    "{got} scored, {} left out (not matched in the panel run)",
                                    n_x - got
                                ),
                            );
                        }
                    }
                }
            }
        }
    }
    let below = matches!(i.outcome, Outcome::BelowThreshold(_));
    let mut rows = Vec::with_capacity(totals.len());
    for (pid, (t, label)) in &totals {
        let m = matched.get(pid).copied().unwrap_or(0);
        let a = applies.get(pid).map_or("", String::as_str);
        let check = if a.is_empty() {
            ""
        } else if sex.is_empty() {
            SEX_NOT_DECLARED
        } else if sex == a {
            SEX_APPLIES
        } else {
            SEX_NOT_APPLICABLE
        };
        let ev: &[Eval] = evals.evals(pid);
        // The effect is printed only for the group pgsc_calc placed the
        // sample in, and only when that group is one the score was evaluated in.
        let g = group.get(pid).map_or("NA", String::as_str);
        let p = pct.get(pid).map_or("NA", String::as_str);
        let eff = if check == SEX_NOT_APPLICABLE || p == "NA" || below {
            String::new()
        } else {
            ev.iter()
                .find(|e| e.ancestry == g)
                .map(|e| format!("{} ({}, {})", e.effect, e.ancestry, e.ppm))
                .unwrap_or_default()
        };
        let evaluated_in: Vec<&str> = ev.iter().map(|e| e.ancestry.as_str()).collect();
        let label = label.clone().unwrap_or_else(|| pid.clone());
        let extra = [
            chrx.get(pid).cloned().unwrap_or_default(),
            ambig.get(pid).cloned().unwrap_or_default(),
            a.to_owned(),
            check.to_owned(),
            evaluated_in.join(";"),
            eff,
            i.input_kind.to_owned(),
        ];
        let head: [String; 8] = if below {
            [
                label,
                pid.clone(),
                "NA".into(),
                "NA".into(),
                t.to_string(),
                rate.get(pid).cloned().unwrap_or_else(|| "NA".into()),
                "NA".into(),
                "NA".into(),
            ]
        } else {
            let pct_shown = if check == SEX_NOT_APPLICABLE { "NA" } else { p };
            [
                label,
                pid.clone(),
                sums.get(pid).cloned().unwrap_or_else(|| "NA".into()),
                m.to_string(),
                t.to_string(),
                if *t != 0 {
                    fixed(100.0 * m as f64 / *t as f64, 1)
                } else {
                    "0.0".into()
                },
                pct_shown.to_owned(),
                g.to_owned(),
            ]
        };
        let mut row: TableRow = Default::default();
        for (k, v) in head.into_iter().chain(extra).enumerate() {
            row[k] = v;
        }
        rows.push(row);
    }
    let mut ancestry = Vec::new();
    if let Some(p) = pops {
        let val = |k: &str| shown(p.get(k), "");
        ancestry.push(("population".to_owned(), val("MostSimilarPop")));
        ancestry.push((
            "population_low_confidence".to_owned(),
            val("MostSimilarPop_LowConfidence"),
        ));
        let keys: Vec<&str> = p.iter().map(|(k, _)| k).collect();
        let mut rf: Vec<&str> = keys
            .iter()
            .copied()
            .filter(|k| k.starts_with("RF_P_"))
            .collect();
        rf.sort_unstable();
        for k in rf {
            ancestry.push((
                format!("probability_{}", &k[5..]),
                fmt_num(p.get(k).flatten()),
            ));
        }
        let mut pcs: Vec<(i128, &str)> = keys
            .iter()
            .filter(|k| k.starts_with("PC") && isdigit(&k[2..]))
            .filter_map(|k| crate::int(&k[2..]).ok().map(|n| (n, *k)))
            .collect();
        pcs.sort_by_key(|(n, _)| *n);
        for (_, k) in pcs {
            ancestry.push((k.to_owned(), fmt_num(p.get(k).flatten())));
        }
    }
    Ok((rows, ancestry))
}

/// The stdout line of one row.
fn line(r: &TableRow) -> String {
    let mut extra = if r[6] == "NA" {
        String::new()
    } else {
        format!(", percentile {} among {}", r[6], r[7])
    };
    if r[11] == SEX_NOT_APPLICABLE {
        extra = format!(", {}-specific: {SEX_NOT_APPLICABLE}, no percentile", r[10]);
    } else if r[11] == SEX_NOT_DECLARED {
        extra.push_str(&format!(", {}-specific score: {SEX_NOT_DECLARED}", r[10]));
    }
    if !r[13].is_empty() {
        extra.push_str(&format!(", effect in that group {}", r[13]));
    }
    if r[3] == "NA" {
        format!(
            "  {} ({}): no sum, under pgsc_calc's minimum overlap ({}% of {} variants matched)",
            r[0], r[1], r[5], r[4]
        )
    } else {
        format!(
            "  {} ({}): sum {}, {} of {} variants matched ({}%){extra}",
            r[0], r[1], r[2], r[3], r[4], r[5]
        )
    }
}

/// `og prs table ARGS`: the exit code.
pub fn cli(argv: &[String], out: &mut impl Write, err: &mut impl Write) -> u8 {
    let a = match args::parse(argv, OPTS) {
        Ok(a) => a,
        Err(e) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og prs table: error: {e}");
            return 2;
        }
    };
    let outcome = match a.get("--below-threshold").filter(|l| !l.is_empty()) {
        Some(log) => Outcome::BelowThreshold(log),
        None if a.flag("--zero-matches") => Outcome::ZeroMatches,
        None => Outcome::Scored,
    };
    let inputs = Inputs {
        results: a.text("--results"),
        sampleset: a.text("--sampleset"),
        scores: a.text("--scores"),
        input_kind: a.text("--input-kind"),
        outcome,
        sex: a.text("--sex"),
        panel: !a.text("--panel").is_empty(),
    };
    let mut run = || -> Result<(), String> {
        let (rows, ancestry) = table(&inputs)?;
        // A score list without evaluations gives the table its older columns.
        let evaluated = rows.iter().any(|r| !r[12].is_empty());
        let mut body = String::new();
        let cols: Vec<&str> = if evaluated {
            COLUMNS[..12]
                .iter()
                .chain(&EFFECT_COLUMNS)
                .chain(&COLUMNS[12..])
                .copied()
                .collect()
        } else {
            COLUMNS.to_vec()
        };
        body.push_str(&cols.join("\t"));
        body.push('\n');
        for r in &rows {
            let cells: Vec<&str> = if evaluated {
                r.iter().map(String::as_str).collect()
            } else {
                r[..12].iter().chain(&r[14..]).map(String::as_str).collect()
            };
            body.push_str(&cells.join("\t"));
            body.push('\n');
        }
        write_file(a.text("--out"), body.as_bytes())?;
        if let Some(anc) = a
            .get("--ancestry-out")
            .filter(|p| !p.is_empty() && !ancestry.is_empty())
        {
            let mut b = String::from("key\tvalue\n");
            let head = [
                ("sample", a.text("--sample")),
                ("reference_panel", a.text("--panel")),
            ];
            for (k, v) in head
                .iter()
                .map(|(k, v)| (*k, *v))
                .chain(ancestry.iter().map(|(k, v)| (k.as_str(), v.as_str())))
            {
                b.push_str(&format!("{k}\t{v}\n"));
            }
            write_file(anc, b.as_bytes())?;
        }
        for r in &rows {
            writeln!(out, "{}", line(r)).map_err(|e| e.to_string())?;
        }
        Ok(())
    };
    finish(run(), err)
}

#[cfg(test)]
mod tests;
