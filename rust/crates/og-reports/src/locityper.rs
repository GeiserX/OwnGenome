//! `og locityper report` (step 53): one row per locus from what Locityper
//! genotype wrote.
//!
//! ```text
//! og locityper report --gt-dir DIR --loci BED --exit-code N --out TSV
//! ```
//!
//! DIR is the output folder of `locityper genotype` (`DIR/loci/<locus>/
//! res.json.gz`, with a `success` file beside each finished locus). BED is
//! the locus list the database was built from (`assets/locityper_loci.bed`:
//! chrom, start, end, name; `#` lines are comments). One row per BED locus,
//! in the BED's order: locus, region, status (called, or not called with
//! the reason), genotype and its two haplotypes, Locityper's quality, the
//! reads used and unexplained, the weighted distance, Locityper's warnings,
//! and the filter of the authors' HLA and KIR benchmark (Prodanov et al.
//! 2025, Zenodo 15879519 README): no warning, a weighted distance of 30 or
//! less, and fewer than 1000 or fewer than 20% unexplained reads.
//!
//! A locus Locityper did not finish is "not called", never a genotype. Exit
//! 1 only when the BED cannot be read.
//!
//! What changes on purpose: a `res.json.gz` whose JSON is not an object,
//! or whose genotype is not text, is "not called" with that reason (the
//! Python stopped with a traceback), and a BED coordinate or a read count
//! must be ASCII digits (the Python's `isdigit` also took `²`, then
//! crashed on it).

use crate::{Error, Outcome};
use og_filters::args::Spec;
use og_io::Utf8;
use pycompat::json::Value;
use std::io::{Read, Write};
use std::path::Path;

pub const COLUMNS: [&str; 13] = [
    "locus",
    "region",
    "status",
    "genotype",
    "haplotype_1",
    "haplotype_2",
    "quality",
    "total_reads",
    "unexplained_reads",
    "weighted_distance",
    "warnings",
    "filter",
    "reason",
];

const MAX_WEIGHT_DIST: f64 = 30.0;
const MAX_UNEXPLAINED: u128 = 1000;
const MAX_UNEXPLAINED_SHARE: f64 = 0.2;

/// The loci of the BED: (name, `chrom:start+1-end`).
pub fn read_loci(path: &str) -> Result<Vec<(String, String)>, Error> {
    let mut loci = Vec::new();
    for (n, line) in crate::read_lines(path, Utf8::Strict)?.iter().enumerate() {
        if crate::strip(line).is_empty() || line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(line)
            .split('\t')
            .collect();
        let start = (c.len() >= 4 && crate::digits(c[1]) && crate::digits(c[2]))
            .then(|| crate::int(c[1]))
            .flatten();
        let Some(start) = start else {
            return Err(Error(format!(
                "{path}:{}: want chrom, start, end, name (tab separated)",
                n + 1
            )));
        };
        loci.push((c[3].to_owned(), format!("{}:{}-{}", c[0], start + 1, c[2])));
    }
    if loci.is_empty() {
        return Err(Error(format!("{path} lists no locus")));
    }
    Ok(loci)
}

/// `float(v)` of a JSON value, `None` where Python raises.
fn as_float(v: Option<&Value>) -> Option<f64> {
    match v? {
        Value::Bool(b) => Some(f64::from(u8::from(*b))),
        Value::Int(i) => Some(*i as f64),
        Value::Float(x) => Some(*x),
        Value::Str(s) => crate::float(s),
        _ => None,
    }
}

/// A number of `res.json.gz` with `digits` decimals; "." when it is absent,
/// no number, or NaN.
fn fmt(v: Option<&Value>, digits: usize) -> String {
    match as_float(v) {
        Some(x) if !x.is_nan() => crate::fixed(x, digits),
        _ => ".".into(),
    }
}

/// A read count as a whole number, `None` when it is absent or not one.
fn count(res: &Value, key: &str) -> Option<String> {
    match pycompat::value::get(res, key).ok()?? {
        Value::Int(i) if *i >= 0 => Some(i.to_string()),
        Value::Float(x) if x.is_finite() && x.fract() == 0.0 && *x >= 0.0 => {
            Some(pycompat::fmt::int_str(*x))
        }
        Value::Str(s) if crate::digits(crate::strip(s)) => crate::int(s).map(|v| v.to_string()),
        _ => None,
    }
}

/// Locityper's warnings: none for null, "", "*" or []; a string split at
/// `;` and `,`; each item of a list.
fn warnings_of(res: &Value) -> Result<Vec<String>, String> {
    match pycompat::value::get(res, "warnings").ok().flatten() {
        None | Some(Value::Null) => Ok(Vec::new()),
        Some(Value::Str(s)) if s.is_empty() || s == "*" => Ok(Vec::new()),
        Some(Value::Str(s)) => Ok(s
            .replace(';', ",")
            .split(',')
            .filter(|x| !crate::strip(x).is_empty())
            .map(str::to_owned)
            .collect()),
        Some(Value::Array(a)) => a
            .iter()
            .map(|x| pycompat::value::py_str(x).map_err(|e| e.to_string()))
            .collect(),
        Some(Value::Object(o)) => Ok(o.keys().cloned().collect()),
        Some(_) => Err("its warnings are neither text nor a list".into()),
    }
}

/// The authors' filter: pass, or fail with the reasons.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Filter {
    Pass,
    Fail(Vec<String>),
}

fn judge(res: &Value, warns: &[String]) -> Filter {
    let mut why = Vec::new();
    if !warns.is_empty() {
        why.push(format!("Locityper warned: {}", warns.join(", ")));
    }
    let wd = as_float(pycompat::value::get(res, "weight_dist").ok().flatten()).unwrap_or(f64::NAN);
    if wd.is_nan() {
        why.push("no weighted distance (the database has no haplotype alignments: setup.sh runs locityper augment)".into());
    } else if wd > MAX_WEIGHT_DIST {
        why.push(format!(
            "weighted distance {} is above {}",
            crate::fixed(wd, 2),
            pycompat::fmt::general(MAX_WEIGHT_DIST, 6)
        ));
    }
    let n = |k| count(res, k).and_then(|s| s.parse::<u128>().ok());
    match (n("total_reads"), n("unexpl_reads")) {
        (Some(total), Some(unexpl)) => {
            let share = MAX_UNEXPLAINED_SHARE * total as f64;
            if !(unexpl < MAX_UNEXPLAINED || (total > 0 && (unexpl as f64) < share)) {
                why.push(format!(
                    "{unexpl} of {total} reads unexplained (1000 or more, and 20% or more)"
                ));
            }
        }
        _ => why.push("Locityper gave no read counts (total_reads, unexpl_reads), so the unexplained-read rule cannot be checked".into()),
    }
    if why.is_empty() {
        Filter::Pass
    } else {
        Filter::Fail(why)
    }
}

/// One row of the table, by [`COLUMNS`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LocusRow {
    NotCalled {
        locus: String,
        region: String,
        reason: String,
    },
    Called {
        locus: String,
        region: String,
        haps: [String; 2],
        quality: String,
        total: String,
        unexplained: String,
        weight_dist: String,
        warnings: Vec<String>,
        filter: Filter,
    },
}

impl LocusRow {
    fn cells(&self) -> Vec<String> {
        let dot = |s: Option<String>| s.unwrap_or_else(|| ".".into());
        match self {
            LocusRow::NotCalled {
                locus,
                region,
                reason,
            } => {
                let mut c = vec![locus.clone(), region.clone(), "not called".into()];
                c.extend(std::iter::repeat_n(".".to_owned(), 9));
                c.push(reason.clone());
                c
            }
            LocusRow::Called {
                locus,
                region,
                haps,
                quality,
                total,
                unexplained,
                weight_dist,
                warnings,
                filter,
            } => {
                let (f, why) = match filter {
                    Filter::Pass => ("pass", None),
                    Filter::Fail(w) => ("fail", Some(w.join("; "))),
                };
                vec![
                    locus.clone(),
                    region.clone(),
                    "called".into(),
                    haps.join(","),
                    haps[0].clone(),
                    haps[1].clone(),
                    quality.clone(),
                    total.clone(),
                    unexplained.clone(),
                    weight_dist.clone(),
                    dot((!warnings.is_empty()).then(|| warnings.join(","))),
                    f.into(),
                    dot(why),
                ]
            }
        }
    }
}

/// `res.json.gz` as JSON, or why it cannot be read.
fn read_res(path: &Path) -> Result<Value, String> {
    let raw = std::fs::read(path).map_err(|e| e.to_string())?;
    if raw.len() < 2 || raw[..2] != [0x1f, 0x8b] {
        return Err(format!(
            "Not a gzipped file ({})",
            bytes_repr(&raw[..raw.len().min(2)])
        ));
    }
    let mut text = Vec::new();
    og_io::reader(std::io::Cursor::new(raw))
        .and_then(|mut r| r.read_to_end(&mut text))
        .map_err(|e| e.to_string())?;
    let text = String::from_utf8(text).map_err(|e| e.to_string())?;
    pycompat::json::loads(&text).map_err(|e| e.to_string())
}

/// Python's `repr()` of a short bytes value.
fn bytes_repr(b: &[u8]) -> String {
    let mut s = String::from("b'");
    for &c in b {
        match c {
            b'\'' => s.push_str("\\'"),
            b'\\' => s.push_str("\\\\"),
            b'\n' => s.push_str("\\n"),
            b'\r' => s.push_str("\\r"),
            b'\t' => s.push_str("\\t"),
            0x20..=0x7e => s.push(c as char),
            _ => s.push_str(&format!("\\x{c:02x}")),
        }
    }
    s.push('\'');
    s
}

/// The row of one locus.
pub fn row_for(name: &str, region: &str, gt_dir: &str, exit_code: i64) -> LocusRow {
    let not_called = |reason: String| LocusRow::NotCalled {
        locus: name.into(),
        region: region.into(),
        reason,
    };
    let d = Path::new(gt_dir).join("loci").join(name);
    let res_path = d.join("res.json.gz");
    if !d.join("success").is_file() || !res_path.is_file() {
        return not_called(if exit_code != 0 {
            format!("Locityper stopped before this locus (exit code {exit_code})")
        } else {
            "Locityper wrote no result for this locus".into()
        });
    }
    let res = match read_res(&res_path) {
        Ok(v @ Value::Object(_)) => v,
        Ok(_) => return not_called("res.json.gz is not a JSON object".into()),
        Err(e) => return not_called(format!("res.json.gz cannot be read: {e}")),
    };
    let gt = match pycompat::value::get(&res, "genotype").ok().flatten() {
        Some(Value::Str(s)) => crate::strip(s).to_owned(),
        Some(v) if pycompat::value::truthy(v) => {
            return not_called("Locityper gave a genotype that is not text".into());
        }
        _ => String::new(),
    };
    let haps: Vec<&str> = gt.split(',').filter(|h| !h.is_empty()).collect();
    let [h1, h2] = haps[..] else {
        let shown = if gt.is_empty() { "empty" } else { &gt };
        return not_called(format!("Locityper gave no diploid genotype ({shown})"));
    };
    let warnings = match warnings_of(&res) {
        Ok(w) => w,
        Err(e) => return not_called(format!("res.json.gz cannot be read: {e}")),
    };
    let filter = judge(&res, &warnings);
    let get = |k| pycompat::value::get(&res, k).ok().flatten();
    LocusRow::Called {
        locus: name.into(),
        region: region.into(),
        haps: [h1.into(), h2.into()],
        quality: fmt(get("quality"), 1),
        total: count(&res, "total_reads").unwrap_or_else(|| ".".into()),
        unexplained: count(&res, "unexpl_reads").unwrap_or_else(|| ".".into()),
        weight_dist: fmt(get("weight_dist"), 3),
        warnings,
        filter,
    }
}

const SPEC: Spec = Spec {
    prog: "og locityper report",
    usage: "usage: og locityper report [-h] --gt-dir GT_DIR --loci LOCI [--exit-code EXIT_CODE] --out OUT\n",
    options: &["--gt-dir", "--loci", "--exit-code", "--out"],
    positionals: &[],
    required: &["--gt-dir", "--loci", "--out"],
};

/// `og locityper VERB ARGS`: the exit code.
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    if verb != Some("report") {
        return crate::unknown_verb("locityper", "report", verb, err);
    }
    crate::run(&SPEC, rest, out, err, |a| {
        let exit_code = crate::int_opt(a, "--exit-code", 0)?;
        Ok(report(
            crate::opt(a, "--gt-dir"),
            crate::opt(a, "--loci"),
            exit_code,
            crate::opt(a, "--out"),
        ))
    })
}

fn report(gt_dir: &str, loci: &str, exit_code: i64, out: &str) -> Outcome {
    let loci = read_loci(loci).map_err(|e| (e, 1))?;
    let rows: Vec<LocusRow> = loci
        .iter()
        .map(|(n, r)| row_for(n, r, gt_dir, exit_code))
        .collect();
    let mut w = pycompat::csv::Writer::new(Vec::new(), '\t');
    let io = |e: std::io::Error| (Error::io(&e, out), 1);
    w.write_row(COLUMNS).map_err(io)?;
    for r in &rows {
        w.write_row(r.cells()).map_err(io)?;
    }
    crate::write_atomic(out, &w.into_inner()).map_err(|e| (e, 1))?;
    let called = rows
        .iter()
        .filter(|r| matches!(r, LocusRow::Called { .. }))
        .count();
    let pass = rows
        .iter()
        .filter(|r| {
            matches!(
                r,
                LocusRow::Called {
                    filter: Filter::Pass,
                    ..
                }
            )
        })
        .count();
    Ok(format!(
        "{called} of {} loci called, {pass} pass the authors' filter\n",
        rows.len()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use pycompat::json::loads;

    fn j(s: &str) -> Value {
        loads(s).unwrap()
    }

    #[test]
    fn counts_as_python_reads_them() {
        let r = j(
            r#"{"a": 3, "b": 3.0, "c": " 12 ", "d": true, "e": -1, "f": 2.5, "g": "²", "h": 1e20}"#,
        );
        assert_eq!(count(&r, "a").as_deref(), Some("3"));
        assert_eq!(count(&r, "b").as_deref(), Some("3"));
        assert_eq!(count(&r, "c").as_deref(), Some("12"));
        assert_eq!(count(&r, "d"), None);
        assert_eq!(count(&r, "e"), None);
        assert_eq!(count(&r, "f"), None);
        assert_eq!(count(&r, "g"), None);
        assert_eq!(count(&r, "h").as_deref(), Some("100000000000000000000"));
        assert_eq!(count(&r, "zz"), None);
    }

    #[test]
    fn filter_rules() {
        let r = j(r#"{"weight_dist": 30.0, "total_reads": 4000, "unexpl_reads": 999}"#);
        assert_eq!(judge(&r, &[]), Filter::Pass);
        let r = j(r#"{"weight_dist": 30.5, "total_reads": 10000, "unexpl_reads": 1999}"#);
        assert_eq!(
            judge(&r, &["W".into()]),
            Filter::Fail(vec![
                "Locityper warned: W".into(),
                "weighted distance 30.50 is above 30".into()
            ])
        );
        let r = j(r#"{"weight_dist": NaN, "total_reads": 4000, "unexpl_reads": 1000}"#);
        let Filter::Fail(w) = judge(&r, &[]) else {
            panic!()
        };
        assert!(w[0].starts_with("no weighted distance"));
        assert_eq!(
            w[1],
            "1000 of 4000 reads unexplained (1000 or more, and 20% or more)"
        );
    }

    /// Red first: a weighted distance of exactly 30 passes.
    #[test]
    fn red_first_a_weighted_distance_of_30_passes() {
        let r = j(r#"{"weight_dist": "30", "total_reads": 10, "unexpl_reads": 0}"#);
        assert_eq!(judge(&r, &[]), Filter::Pass);
    }

    #[test]
    fn warnings() {
        assert_eq!(
            warnings_of(&j(r#"{"warnings": "*"}"#)).unwrap(),
            Vec::<String>::new()
        );
        assert_eq!(
            warnings_of(&j(r#"{"warnings": "a; b,,c"}"#)).unwrap(),
            ["a", " b", "c"]
        );
        assert_eq!(
            warnings_of(&j(r#"{"warnings": ["x", 1, null]}"#)).unwrap(),
            ["x", "1", "None"]
        );
        assert_eq!(fmt(Some(&j("1e400")), 1), "inf");
        assert_eq!(fmt(Some(&j("true")), 1), "1.0");
        assert_eq!(fmt(Some(&j("\"x\"")), 1), ".");
    }
}
