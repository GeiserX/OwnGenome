//! The run facts both reports carry: the sample, when the summary was
//! written and when the latest run started, the OwnGenome release and the
//! pipeline commit from `run_manifest.tsv`, the database releases, and the
//! last database re-check.

use pycompat::json::Value;

use crate::escape::e;
use crate::model::{Malformed, Summary, get, is_str, py, truthy};

/// The value of the manifest's `run KEY` row, if it has one.
fn run_row<'a>(rows: &[[&'a Value; 3]], key: &str) -> Option<&'a Value> {
    rows.iter()
        .find(|r| is_str(Some(r[0]), "run") && is_str(Some(r[1]), key))
        .map(|r| r[2])
}

fn field<'a>(i: &'a Value, k: &str) -> Result<Option<&'a Value>, Malformed> {
    get(i, k, "summary.recheck.items[]")
}

/// One line of run facts from the last database re-check: what is newer
/// than this run, or that nothing is.
pub fn recheck_line(r: &Value) -> Result<String, Malformed> {
    let items = get(r, "items", "summary.recheck")?.filter(|v| truthy(v));
    let mut newer = Vec::new();
    if let Some(items) = items {
        let Value::Array(items) = items else {
            return Err(Malformed("summary.recheck.items is not a list".into()));
        };
        for i in items {
            match field(i, "class")? {
                None | Some(Value::Null) => {}
                Some(c) if is_str(Some(c), "no-newer-release") => {}
                Some(c) => {
                    let label = field(i, "label")?.filter(|v| truthy(v));
                    let name = label.or(field(i, "item")?);
                    newer.push(format!(
                        "{} {}",
                        name.map_or(Ok("None".into()), py)?,
                        py(c)?
                    ));
                }
            }
        }
    }
    let date = get(r, "checked_utc", "summary.recheck")?
        .filter(|v| truthy(v))
        .map_or(Ok("unknown date".into()), py)?;
    let head = format!("Database re-check: {date}: ");
    Ok(if newer.is_empty() {
        head + "nothing newer than this run"
    } else {
        head + &newer.join("; ")
    })
}

/// The re-check line, when the summary records a re-check.
fn recheck(s: &Summary) -> Result<Option<String>, Malformed> {
    s.get("recheck")
        .filter(|r| truthy(r))
        .map(recheck_line)
        .transpose()
}

/// When the latest run-all.sh run started, if the run status says.
fn started(s: &Summary) -> Result<Option<&Value>, Malformed> {
    Ok(Some(s.need2("run", "started_utc")?).filter(|v| truthy(v)))
}

/// The text report's opening lines, up to the blank line after the rule.
pub fn text_head(s: &Summary) -> Result<Vec<String>, Malformed> {
    let rule = "=".repeat(80);
    let mut w = vec![
        rule.clone(),
        "  GENOMICS ANALYSIS REPORT".to_owned(),
        format!("  Sample: {}", py(s.need("sample")?)?),
        format!("  Generated: {}", py(s.need("generated_utc")?)?),
    ];
    if let Some(t) = started(s)? {
        w.push(format!("  Latest run-all.sh run: {}", py(t)?));
    }
    if let Some(line) = recheck(s)? {
        w.push(format!("  {line}"));
    }
    w.push(rule);
    w.push(String::new());
    Ok(w)
}

/// The text report's manifest line, before the disclaimer; none without a
/// manifest.
pub fn text_manifest(s: &Summary) -> Result<Option<String>, Malformed> {
    let rows = s.manifest()?;
    if rows.is_empty() {
        return Ok(None);
    }
    let commit = run_row(&rows, "git_commit").map_or(Ok("unknown".into()), py)?;
    let release = match run_row(&rows, "pipeline_version").filter(|v| truthy(v)) {
        Some(v) => format!("OwnGenome {}, ", py(v)?),
        None => String::new(),
    };
    let images = rows.iter().filter(|r| is_str(Some(r[0]), "image")).count();
    Ok(Some(format!(
        "Run manifest: run_manifest.tsv ({release}pipeline commit {commit}, {images} images)"
    )))
}

/// The HTML header's lines, from `<header>` to `</header>`.
pub fn html_header(s: &Summary) -> Result<Vec<String>, Malformed> {
    let rows = s.manifest()?;
    let commit = run_row(&rows, "git_commit").filter(|v| truthy(v));
    let release = run_row(&rows, "pipeline_version").filter(|v| truthy(v));
    let db = |k: &str| -> Result<String, Malformed> {
        let v = s.need2("databases", k)?;
        if truthy(v) {
            py(v)
        } else {
            Ok("unknown".into())
        }
    };
    let mut first = format!(
        "    Sample: <strong>{}</strong> &nbsp;|&nbsp; Generated: {}",
        e(&py(s.need("sample")?)?),
        e(&py(s.need("generated_utc")?)?)
    );
    if let Some(t) = started(s)? {
        first += &format!(" &nbsp;|&nbsp; Latest run: {}", e(&py(t)?));
    }
    let mut second = match release {
        Some(r) => format!("    OwnGenome release: {} &nbsp;|&nbsp;", e(&py(r)?)),
        None => String::new(),
    };
    second += &format!(
        "    Pipeline commit: {} &nbsp;|&nbsp; ClinVar file: {} &nbsp;|&nbsp; ClinVar in CPSR's bundle: {} \
         &nbsp;|&nbsp; HLA database: {}",
        e(&commit.map_or(Ok("unknown".into()), py)?),
        e(&db("clinvar_release")?),
        e(&db("cpsr_clinvar_release")?),
        e(&db("hla_database")?),
    );
    let mut out = vec![
        "<header class=\"header\">".to_owned(),
        "  <h1>Genomic Analysis Report</h1>".to_owned(),
        "  <div class=\"meta\">".to_owned(),
        first,
        "  </div>".to_owned(),
        "  <div class=\"meta\">".to_owned(),
        second,
        "  </div>".to_owned(),
    ];
    if let Some(line) = recheck(s)? {
        out.push(format!("  <div class=\"meta\">{}</div>", e(&line)));
    }
    out.push("</header>".to_owned());
    Ok(out)
}

/// The footer's manifest: every row in a folded table, or a line saying
/// there is none.
pub fn html_manifest(s: &Summary) -> Result<Vec<String>, Malformed> {
    let rows = s.manifest()?;
    if rows.is_empty() {
        return Ok(vec![
            "  <p>No run_manifest.tsv in the sample folder.</p>".to_owned(),
        ]);
    }
    let mut out = vec![
        "  <details><summary>Run manifest (run_manifest.tsv)</summary>".to_owned(),
        "  <table>".to_owned(),
    ];
    for r in &rows {
        out.push(format!(
            "    <tr><td>{}</td><td>{}</td><td>{}</td></tr>",
            e(&py(r[0])?),
            e(&py(r[1])?),
            e(&py(r[2])?)
        ));
    }
    out.push("  </table></details>".to_owned());
    Ok(out)
}
