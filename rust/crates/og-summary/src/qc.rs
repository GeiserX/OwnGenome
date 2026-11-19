//! `og qc sample`: step 33's verdict for one sample, `qc/<id>_sample_qc.tsv`.
//!
//! Is this sample the person it is said to be, and is anyone else's DNA in
//! it (docs/33-sample-qc.md)? It reads somalier relate's `samples.tsv` (the
//! sex it infers from chrX and chrY) and `pairs.tsv` (relatedness to the
//! other samples of the run), and VerifyBamID2's `.selfSM` (FREEMIX, the
//! share of reads that look like another person's).
//!
//! - somalier's own chrX rule is applied again ([`x_sex`]): more than 10
//!   chrX sites, heterozygous over homozygous-ALT under 0.05 is male, over
//!   0.4 female, and under 6% of sites with a middling allele balance.
//! - With the declared sex as its pedigree, somalier starts its sex column
//!   from that sex; a pedigree sex the reads did not confirm is not a call.
//! - `-2` (chrX looks female, chrY has reads) is reported even with no
//!   declared sex.
//! - FREEMIX above `--freemix-warn` (default 0.03) warns.
//! - A relatedness of 0.9 or more is the same person (a duplicate or a swap).
//! - FREEMIX is `uncertain` when VerifyBamID2 skipped its marker check or
//!   kept fewer than [`MIN_CONTAMINATION_MARKERS`] markers.
//! - Identity is `uncertain` below [`MIN_IDENTITY_SITES`] informative sites,
//!   and `not_checked` when no other sample was there to compare with.
//!
//! The table is `key\tvalue` with 25 rows in a fixed order; `og qc
//! sex-loss` appends its `mosaic_*` rows afterwards. Each row is echoed on
//! stdout as `  key: value`.

use crate::args::{self, Kind, Opt};
use crate::{basename, finish, float, text_lines, trunc, write_file};
use pycompat::fmt::{fixed, general};
use pycompat::text::lower;
use std::collections::HashMap;
use std::io::Write;

/// VerifyBamID2's FREEMIX above this is reported as possible contamination.
pub const FREEMIX_WARN: f64 = 0.03;
/// somalier relatedness at or above this: the same person.
pub const SAME_PERSON: f64 = 0.9;
/// Informative sites (genotyped by somalier, or shared by a pair) below
/// which identity is `uncertain`: about one in 18 of somalier's 17,766
/// sites. A whole genome at 30x genotypes most of them; a sliced, targeted
/// or very shallow BAM falls below the floor.
pub const MIN_IDENTITY_SITES: i128 = 1000;
/// Markers VerifyBamID2 kept below which FREEMIX is `uncertain`: five times
/// the 1,000 its own marker check asks for, and half of the 10% of the
/// 100,000-marker panel that check also asks for.
pub const MIN_CONTAMINATION_MARKERS: u64 = 5000;

const USAGE: &str =
    "usage: og qc sample --sample S --somalier-samples F [--somalier-pairs F] [--somalier-id ID]
                    [--selfsm F] [--declared-sex {male,female,}] [--freemix-warn X]
                    [--marker-check {passed,skipped,}] [--verifybamid2-markers N]
                    --out F
";

const OPTS: &[Opt] = &[
    Opt("--sample", Kind::Value, true),
    Opt("--somalier-samples", Kind::Value, true),
    Opt("--somalier-pairs", Kind::Value, false),
    Opt("--somalier-id", Kind::Value, false),
    Opt("--selfsm", Kind::Value, false),
    Opt(
        "--declared-sex",
        Kind::Choice(&["male", "female", ""]),
        false,
    ),
    Opt("--freemix-warn", Kind::Value, false),
    Opt(
        "--marker-check",
        Kind::Choice(&["passed", "skipped", ""]),
        false,
    ),
    Opt("--verifybamid2-markers", Kind::Value, false),
    Opt("--out", Kind::Value, true),
];

/// A sex somalier or the declaration names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
}

impl Sex {
    fn word(self) -> &'static str {
        match self {
            Sex::Male => "male",
            Sex::Female => "female",
        }
    }

    fn from_word(s: &str) -> Option<Sex> {
        match s {
            "male" => Some(Sex::Male),
            "female" => Some(Sex::Female),
            _ => None,
        }
    }
}

/// One row of a somalier or VerifyBamID2 table: column to cell. A header
/// that starts with `#` loses it on its first column only.
pub type Row = HashMap<String, String>;

/// The rows of a tab-separated file whose header may start with `#`. Blank
/// lines are skipped; a row shorter than the header lacks the last columns.
pub fn read_rows(path: &str) -> Result<Vec<Row>, String> {
    let mut lines = text_lines(path)?;
    let mut head: Option<Vec<String>> = None;
    let mut rows = Vec::new();
    for line in &mut lines {
        let line = line.map_err(|e| format!("{path}: {e}"))?;
        if pycompat::text::strip(&line, None).is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        match &head {
            None => {
                let mut h: Vec<String> = cells.iter().map(|c| (*c).to_owned()).collect();
                h[0] = h[0].trim_start_matches('#').to_owned();
                head = Some(h);
            }
            Some(h) => rows.push(
                h.iter()
                    .zip(&cells)
                    .map(|(k, v)| (k.clone(), (*v).to_owned()))
                    .collect(),
            ),
        }
    }
    Ok(rows)
}

/// `num(x)`: a number, or nothing for a missing cell, text that is not a
/// number, and NaN.
fn num(x: Option<&str>) -> Option<f64> {
    x.and_then(|s| float(s).ok()).filter(|v| !v.is_nan())
}

fn cell<'a>(r: &'a Row, k: &str) -> Option<&'a str> {
    r.get(k).map(String::as_str)
}

/// The sex somalier's own rule calls from chrX (relate --infer,
/// `add_parents_and_check_sex`), or nothing when it cannot tell.
pub fn x_sex(row: &Row) -> Option<Sex> {
    let (n, het, hom) = (
        num(cell(row, "X_n"))?,
        num(cell(row, "X_het"))?,
        num(cell(row, "X_hom_alt"))?,
    );
    // A missing column passes (older somalier); a value that is not a
    // number (nan from 0/0) fails, as somalier's own `< 0.06` does.
    if let Some(raw) = cell(row, "p_middling_ab") {
        match num(Some(raw)) {
            Some(m) if m < 0.06 => {}
            _ => return None,
        }
    }
    if n <= 10.0 {
        return None;
    }
    if hom == 0.0 {
        return (het > 0.0).then_some(Sex::Female);
    }
    let ratio = het / hom;
    if ratio < 0.05 {
        Some(Sex::Male)
    } else if ratio > 0.4 {
        Some(Sex::Female)
    } else {
        None
    }
}

/// The inputs of the verdict.
#[derive(Debug, Default)]
pub struct Inputs<'a> {
    pub sample: &'a str,
    pub samples_tsv: &'a str,
    pub selfsm: Option<&'a str>,
    pub pairs_tsv: Option<&'a str>,
    pub declared: Option<Sex>,
    pub freemix_warn: f64,
    pub somalier_id: Option<&'a str>,
    pub marker_check: &'a str,
    /// The markers VerifyBamID2 kept (its log's "N SNP markers remained
    /// after sanity check"), when known.
    pub markers: Option<u64>,
}

/// The sex check: its word in the table and the reason.
#[derive(Debug, PartialEq, Eq)]
enum SexCheck {
    Ok,
    Mismatch,
    NotChecked,
}

impl SexCheck {
    fn word(&self) -> &'static str {
        match self {
            SexCheck::Ok => "ok",
            SexCheck::Mismatch => "mismatch",
            SexCheck::NotChecked => "not_checked",
        }
    }
}

/// The verdict for one sample, as ordered (key, value) rows.
pub fn table(i: &Inputs) -> Result<Vec<(&'static str, String)>, String> {
    let sid = i.somalier_id.filter(|s| !s.is_empty()).unwrap_or(i.sample);
    let rows = read_rows(i.samples_tsv)?;
    let Some(row) = rows.iter().find(|r| cell(r, "sample_id") == Some(sid)) else {
        let names: Vec<&str> = rows
            .iter()
            .map(|r| cell(r, "sample_id").unwrap_or("?"))
            .collect();
        let names = if names.is_empty() {
            "none".to_owned()
        } else {
            names.join(", ")
        };
        return Err(format!(
            "{} has no row for sample '{sid}' (it has: {names})",
            basename(i.samples_tsv)
        ));
    };
    let sex_col = cell(row, "sex").unwrap_or("");
    let mut inferred = match sex_col {
        "1" => Some(Sex::Male),
        "2" => Some(Sex::Female),
        _ => None,
    };
    // With the declared sex as its pedigree (--ped), somalier starts its sex
    // column from that sex and overwrites it only when the reads tell, so
    // the pedigree's sex left standing is not a call.
    let ped = Sex::from_word(&lower(cell(row, "original_pedigree_sex").unwrap_or("")));
    let x_says = ped.and_then(|_| x_sex(row));
    if ped.is_some() && inferred == ped && x_says != ped {
        inferred = None;
    }
    let x_and_y = sex_col == "-2" && (ped.is_none() || x_says == Some(Sex::Female));
    let inferred_word = inferred.map_or("unknown", Sex::word);
    let (check, why) = if x_and_y {
        (
            SexCheck::NotChecked,
            "chrX is heterozygous like a female sample but chrY has reads: a sex-chromosome aneuploidy or a mixed sample".to_owned(),
        )
    } else if let Some(declared) = i.declared {
        match inferred {
            None => (
                SexCheck::NotChecked,
                format!(
                    "somalier could not tell the sex from {} chrX sites (it needs more than 10 with reads, and allele balances that look like one person)",
                    cell(row, "X_n").unwrap_or("0")
                ),
            ),
            Some(s) => (
                if s == declared {
                    SexCheck::Ok
                } else {
                    SexCheck::Mismatch
                },
                format!("declared {}, somalier infers {}", declared.word(), s.word()),
            ),
        }
    } else {
        (SexCheck::NotChecked, "no declared sex".to_owned())
    };
    let gt = num(cell(row, "gt_depth_mean"));
    let ratio = |k: &str| match (gt, num(cell(row, k))) {
        (Some(g), Some(d)) if g != 0.0 => fixed(d / g, 2),
        _ => String::new(),
    };
    let mut sites: i128 = 0;
    for k in ["n_hom_ref", "n_het", "n_hom_alt"] {
        if let Some(v) = num(cell(row, k)) {
            sites += trunc(v).map_err(|e| {
                format!(
                    "{}: {k} of '{sid}' is not a count: {e}",
                    basename(i.samples_tsv)
                )
            })?;
        }
    }
    let raw = |k: &str| cell(row, k).unwrap_or("").to_owned();
    let mut out = vec![
        ("sample", i.sample.to_owned()),
        ("somalier_id", sid.to_owned()),
        ("declared_sex", i.declared.map_or("", Sex::word).to_owned()),
        ("inferred_sex", inferred_word.to_owned()),
        ("sex_check", check.word().to_owned()),
        ("sex_check_reason", why),
        ("sites_genotyped", sites.to_string()),
        ("depth_at_sites", raw("gt_depth_mean")),
        ("x_sites", raw("X_n")),
        ("x_het", raw("X_het")),
        ("x_hom_alt", raw("X_hom_alt")),
        ("x_depth_ratio", ratio("X_depth_mean")),
        ("y_sites", raw("Y_n")),
        ("y_depth_ratio", ratio("Y_depth_mean")),
    ];
    let (mut freemix, mut status, mut markers) = (None, "not_run", String::new());
    if let Some(selfsm) = i.selfsm.filter(|s| !s.is_empty()) {
        let sm = read_rows(selfsm)?;
        freemix = sm.first().and_then(|r| num(cell(r, "FREEMIX")));
        let Some(f) = freemix else {
            return Err(format!("{} has no FREEMIX value", basename(selfsm)));
        };
        status = if f > i.freemix_warn { "warn" } else { "ok" };
        // #SNPS of the selfSM is the size of the panel, not the markers with reads.
        markers = sm
            .first()
            .and_then(|r| cell(r, "#SNPS"))
            .unwrap_or("")
            .to_owned();
    }
    out.extend([
        ("freemix", freemix.map_or_else(String::new, |f| fixed(f, 4))),
        ("freemix_warn_above", general(i.freemix_warn, 6)),
        ("contamination", status.to_owned()),
        ("panel_markers", markers),
        ("verifybamid2_marker_check", i.marker_check.to_owned()),
    ]);
    let mut same = Vec::new();
    // Other samples compared on enough shared sites, and on too few.
    let (mut compared, mut thin) = (0usize, 0usize);
    if let Some(pairs) = i.pairs_tsv.filter(|s| !s.is_empty()) {
        for r in read_rows(pairs)? {
            let (a, b) = (cell(&r, "sample_a"), cell(&r, "sample_b"));
            let rel = num(cell(&r, "relatedness"));
            if (a == Some(sid) || b == Some(sid)) && a != b {
                // `n`: the sites both samples have a genotype at. A pair
                // with fewer, with no `n` column, or with no relatedness
                // (nan: they share no site) is no comparison.
                match (rel, num(cell(&r, "n"))) {
                    (Some(_), Some(n)) if n >= MIN_IDENTITY_SITES as f64 => compared += 1,
                    _ => thin += 1,
                }
                if let Some(rel) = rel.filter(|&x| x >= SAME_PERSON) {
                    let other = if a == Some(sid) { b } else { a };
                    same.push(format!("{} ({})", other.unwrap_or("None"), fixed(rel, 2)));
                }
            }
        }
    }
    out.push(("same_person_as", same.join(", ")));
    let (confidence, confidence_why) = contamination_confidence(freemix, i.marker_check, i.markers);
    let (identity, identity_why) = identity(sites, compared, thin);
    out.extend([
        (
            "verifybamid2_markers",
            i.markers.map_or_else(String::new, |m| m.to_string()),
        ),
        ("contamination_confidence", confidence.to_owned()),
        ("contamination_confidence_reason", confidence_why),
        ("identity", identity.to_owned()),
        ("identity_reason", identity_why),
    ]);
    Ok(out)
}

/// How far FREEMIX can be trusted: `ok`, `uncertain`, or nothing when
/// VerifyBamID2 did not run; and why.
fn contamination_confidence(
    freemix: Option<f64>,
    marker_check: &str,
    markers: Option<u64>,
) -> (&'static str, String) {
    let floor = MIN_CONTAMINATION_MARKERS;
    if freemix.is_none() {
        return ("", String::new());
    }
    match (marker_check, markers) {
        // VerifyBamID2 keeps its check only for more than 1,000 markers
        // that are also more than a tenth of the panel, so on the
        // 100,000-marker panel a skipped check can follow a count above
        // the floor: the reason names the check's rule, not the floor.
        ("skipped", Some(m)) => (
            "uncertain",
            format!(
                "VerifyBamID2 skipped its marker check: its first run kept {m} markers, too few for that check (more than 1,000 and more than a tenth of the panel), so FREEMIX is a rough estimate"
            ),
        ),
        ("skipped", None) => (
            "uncertain",
            "VerifyBamID2 skipped its marker check: too few markers had reads (it needs more than 1,000 and more than a tenth of the panel)".to_owned(),
        ),
        (_, Some(m)) if m < floor => (
            "uncertain",
            format!(
                "FREEMIX rests on {m} markers; a firm estimate needs {}",
                thousands(floor)
            ),
        ),
        (_, Some(m)) => ("ok", format!("FREEMIX rests on {m} markers")),
        ("passed", None) => ("ok", "VerifyBamID2's marker check passed".to_owned()),
        (_, None) => (
            "uncertain",
            "the number of markers FREEMIX rests on is not known".to_owned(),
        ),
    }
}

/// Whether the reads were checked against another sample: `checked`,
/// `uncertain` (too few informative sites) or `not_checked` (nothing to
/// compare with); and why. `compared` counts the other samples of the
/// run with at least [`MIN_IDENTITY_SITES`] shared sites, `thin` the
/// others (fewer, unknown, or no relatedness).
fn identity(sites: i128, compared: usize, thin: usize) -> (&'static str, String) {
    let floor = thousands(MIN_IDENTITY_SITES as u64);
    if sites < MIN_IDENTITY_SITES {
        return (
            "uncertain",
            format!(
                "somalier genotyped {sites} sites; telling one person from another needs {floor}"
            ),
        );
    }
    if compared > 0 {
        let s = if compared == 1 { "" } else { "s" };
        return (
            "checked",
            format!(
                "compared with {compared} other sample{s} of the run on {floor} or more shared sites"
            ),
        );
    }
    if thin > 0 {
        return (
            "uncertain",
            format!(
                "no other sample of the run is known to share {floor} or more genotyped sites with this one"
            ),
        );
    }
    (
        "not_checked",
        "no other sample in the somalier run to compare with: the sex check is the only identity check".to_owned(),
    )
}

/// `n` with a comma between thousands: 5000 is `5,000`.
fn thousands(n: u64) -> String {
    let d = n.to_string();
    let mut out = String::new();
    for (k, c) in d.chars().enumerate() {
        if k > 0 && (d.len() - k).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// `og qc sample ARGS`: the exit code.
pub fn cli(argv: &[String], out: &mut impl Write, err: &mut impl Write) -> u8 {
    let a = match args::parse(argv, OPTS) {
        Ok(a) => a,
        Err(e) => {
            let _ = write!(err, "{USAGE}");
            let _ = writeln!(err, "og qc sample: error: {e}");
            return 2;
        }
    };
    let freemix_warn = match a.get("--freemix-warn") {
        None => FREEMIX_WARN,
        Some(s) => match float(s) {
            Ok(v) => v,
            Err(_) => {
                let _ = write!(err, "{USAGE}");
                let _ = writeln!(
                    err,
                    "og qc sample: error: argument --freemix-warn: invalid float value: {}",
                    crate::repr(s)
                );
                return 2;
            }
        },
    };
    let markers = match a.get("--verifybamid2-markers").filter(|s| !s.is_empty()) {
        None => None,
        Some(s) => match s.parse::<u64>() {
            Ok(v) => Some(v),
            Err(_) => {
                let _ = write!(err, "{USAGE}");
                let _ = writeln!(
                    err,
                    "og qc sample: error: argument --verifybamid2-markers: invalid int value: {}",
                    crate::repr(s)
                );
                return 2;
            }
        },
    };
    let inputs = Inputs {
        sample: a.text("--sample"),
        samples_tsv: a.text("--somalier-samples"),
        selfsm: a.get("--selfsm"),
        pairs_tsv: a.get("--somalier-pairs"),
        declared: Sex::from_word(a.text("--declared-sex")),
        freemix_warn,
        somalier_id: a.get("--somalier-id"),
        marker_check: a.text("--marker-check"),
        markers,
    };
    let mut run = || -> Result<(), String> {
        let rows = table(&inputs)?;
        let mut body = String::from("key\tvalue\n");
        for (k, v) in &rows {
            body.push_str(&format!("{k}\t{v}\n"));
        }
        write_file(a.text("--out"), body.as_bytes())?;
        for (k, v) in &rows {
            writeln!(out, "  {k}: {v}").map_err(|e| e.to_string())?;
        }
        Ok(())
    };
    finish(run(), err)
}

#[cfg(test)]
mod tests;
