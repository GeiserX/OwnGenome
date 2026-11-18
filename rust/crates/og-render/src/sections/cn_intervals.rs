//! Step 35's copy-number intervals from read depth: ClinGen's recurrent CNV
//! regions, carrier deletion alleles and the AZF regions. One summary line,
//! then the rows with a deletion or a duplication.

use pycompat::json::Value;

use crate::escape::e;
use crate::model::{Inline, Line, Malformed, get, need, py, truthy};

/// The label of the summary line.
pub const LABEL: &str =
    "Copy-number intervals (ClinGen recurrent regions, carrier deletions, AZF; from depth)";

const PATH: &str = "sections.parascopy.data.cn_intervals";

/// `(iv.get(key) or default)` for a dict or a list.
fn or_empty<'a>(iv: &'a Value, key: &str) -> Result<Option<&'a Value>, Malformed> {
    Ok(get(iv, key, PATH)?.filter(|v| truthy(v)))
}

/// The rows with a change.
fn changes(iv: Option<&Value>) -> Result<&[Value], Malformed> {
    let Some(iv) = iv.filter(|v| truthy(v)) else {
        return Ok(&[]);
    };
    match or_empty(iv, "changes")? {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!("{PATH}.changes is not a list"))),
    }
}

/// The interval table in one line: rows by status and how many changed.
pub fn summary(iv: Option<&Value>) -> Result<String, Malformed> {
    let Some(iv) = iv.filter(|v| !matches!(v, Value::Null)) else {
        return Ok(format!(
            "{LABEL}: not assessed (step 35's <id>_cn_intervals.tsv is not in this run's output)"
        ));
    };
    let by = or_empty(iv, "by_status")?;
    let mut counts: Vec<(String, String)> = Vec::new();
    if let Some(by) = by {
        let Value::Object(o) = by else {
            return Err(Malformed(format!("{PATH}.by_status is not an object")));
        };
        for (k, n) in o {
            counts.push((k.clone(), py(n)?));
        }
    }
    counts.sort();
    let mut counts: String = counts
        .iter()
        .map(|(k, n)| format!("{n} {k}"))
        .collect::<Vec<_>>()
        .join(", ");
    if counts.is_empty() {
        counts = "no row".to_owned();
    }
    let n = changes(Some(iv))?.len();
    let called = match by {
        Some(by) => get(by, "called", PATH)?,
        None => None,
    };
    let end = if n > 0 {
        format!("{n} with a deletion or a duplication")
    } else if let Some(k) = called.filter(|k| truthy(k)) {
        format!(
            "no deletion or duplication among the {} called rows",
            py(k)?
        )
    } else {
        "none called".to_owned()
    };
    let rows = get(iv, "rows", PATH)?.map_or(Ok("0".to_owned()), py)?;
    Ok(format!("{LABEL}: {rows} rows ({counts}); {end}"))
}

/// `chrom:start-end` (1-based) of an interval row of the table (BED).
pub fn place(r: &Value) -> Result<String, Malformed> {
    let at = |k: &str| get(r, k, "an interval row");
    let dot = |k: &str| at(k).and_then(|v| v.map_or(Ok(".".to_owned()), py));
    let b = dot("start")?;
    let start = if pycompat::text::isdigit(&b) {
        let n = pycompat::num::parse_int(&b)
            .map_err(|_| Malformed(format!("an interval row's start {b} is not a number")))?;
        (n + 1).to_string()
    } else {
        b
    };
    Ok(format!("{}:{start}-{}", dot("chrom")?, dot("end")?))
}

/// `str.capitalize()`: the first character upper case, the rest lower.
fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        Some(first) => {
            first.to_uppercase().collect::<String>() + &pycompat::text::lower(c.as_str())
        }
        None => String::new(),
    }
}

/// The text report's lines: the summary, then one line per row with a
/// change.
pub fn lines(iv: Option<&Value>) -> Result<Vec<String>, Malformed> {
    let mut out = vec![summary(iv)?];
    for r in changes(iv)? {
        let f = |k: &str| need(r, k, "an interval row").and_then(py);
        out.push(format!(
            "  {}: {} ({}, {}), {} of {} copies",
            capitalize(&f("change")?),
            f("name")?,
            f("locus")?,
            place(r)?,
            f("copies")?,
            f("expected")?
        ));
    }
    Ok(out)
}

/// The HTML card's lines: the summary as a stat, then a table of the rows
/// with a change.
pub fn html(iv: Option<&Value>) -> Result<Vec<Line>, Malformed> {
    let s = summary(iv)?;
    let (label, value) = s.split_once(": ").unwrap_or((&s, ""));
    let mut out = vec![Line::Stat(label.to_owned(), Inline::Text(value.to_owned()))];
    let rows = changes(iv)?;
    if !rows.is_empty() {
        let mut t = "    <table><tr><th>Region</th><th>Row</th><th>Place</th><th>Copies</th>\
                     <th>Change</th></tr>"
            .to_owned();
        for r in rows {
            let f = |k: &str| need(r, k, "an interval row").and_then(py);
            t += &format!(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{} of {}</td><td>{}</td></tr>",
                e(&f("name")?),
                e(&f("locus")?),
                e(&place(r)?),
                e(&f("copies")?),
                e(&f("expected")?),
                e(&f("change")?)
            );
        }
        t += "</table>";
        out.push(Line::Raw(t));
    }
    Ok(out)
}
