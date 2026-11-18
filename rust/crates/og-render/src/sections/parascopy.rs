//! Step 35 (Parascopy, short reads): SMN1 and SMN2 with the SMN1 marker,
//! the other loci of its catalog, and the depth counts it writes beside them
//! (HBA, CYP21A2, the dosage loci, the copy-number intervals).

use pycompat::json::Value;

use crate::escape::e;
use crate::model::{Line, Malformed, Section, get, is_str, need, py, truthy};
use crate::notes::PARASCOPY_LOCI_NOTE;
use crate::sections::smn1::{SMN1_MARKER_CAVEAT, marker_state};
use crate::sections::{cn_intervals, cyp21a2::cyp21a2_line, dosage::dosage_lines, hba::hba_line};

const PATH: &str = "sections.parascopy.data";

/// `(d.get(key) or {})`.
fn dict_or_empty<'a>(d: &Section<'a>, key: &str, empty: &'a Value) -> Result<&'a Value, Malformed> {
    Ok(d.get(key)?.filter(|v| truthy(v)).unwrap_or(empty))
}

/// `d.get(key, ".")` as text.
fn dot(v: &Value, key: &str, path: &str) -> Result<String, Malformed> {
    get(v, key, path)?.map_or(Ok(".".to_owned()), py)
}

/// Step 35's lines, the same in both reports (the HTML card splits each at
/// its first `": "` into a label and a value).
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let empty = Value::Object(Default::default());
    let smn1 = dict_or_empty(d, "smn1", &empty)?;
    let smn2 = dict_or_empty(d, "smn2", &empty)?;
    let mut out = vec![format!(
        "SMN1 copies: {}",
        dot(smn1, "text", &format!("{PATH}.smn1"))?
    )];
    if let Some(mk) = d.get("smn1_marker")?.filter(|v| truthy(v)) {
        let path = format!("{PATH}.smn1_marker");
        let rsid = match get(mk, "rsid", &path)?.filter(|v| truthy(v)) {
            Some(v) => py(v)?,
            None => ".".to_owned(),
        };
        out.push(format!(
            "SMN1 g.27134T>G silent-carrier marker ({rsid}): {}",
            marker_state(mk)?
        ));
        out.push(format!(
            "SMN1 marker reading: {}",
            dot(mk, "message", &path)?
        ));
        out.push(format!("SMN1 marker caveat: {SMN1_MARKER_CAVEAT}"));
    }
    out.push(format!(
        "SMN2 copies: {}",
        dot(smn2, "text", &format!("{PATH}.smn2"))?
    ));
    let f = |k: &str| d.get(k).and_then(|v| v.map_or(Ok(".".to_owned()), py));
    out.push(format!(
        "SMN1 and SMN2 together (agCN) at exon 7: {} (quality {}, filter {})",
        f("agcn")?,
        f("agcn_qual")?,
        f("agcn_filter")?
    ));
    out.push(format!("Reading: {}", f("reading")?));
    if let Some(h) = d.get("hba")?.filter(|v| truthy(v)) {
        out.push(hba_line(h)?);
    }
    out.push(cyp21a2_line(d.get("cyp21a2")?)?);
    out.extend(dosage_lines(d.get("dosage")?)?);
    Ok(out)
}

/// The other loci of the catalog (`d.get("loci") or []`).
fn loci<'a>(d: &Section<'a>) -> Result<&'a [Value], Malformed> {
    match d.get("loci")?.filter(|v| truthy(v)) {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!("{PATH}.loci is not a list"))),
    }
}

/// `genes` with a space after each comma.
fn genes(x: &Value) -> Result<String, Malformed> {
    Ok(py(need(x, "genes", &format!("{PATH}.loci[]"))?)?.replace(',', ", "))
}

/// One line per catalog locus other than SMN1, the same in both reports.
pub fn locus_lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let path = format!("{PATH}.loci[]");
    let mut out = Vec::new();
    for x in loci(d)? {
        let mut line = format!("{} (agCN): {}", genes(x)?, py(need(x, "text", &path)?)?);
        if is_str(get(x, "status", &path)?, "called")
            && let Some(r) = get(x, "reading", &path)?.filter(|r| !is_str(Some(r), "."))
        {
            line += &format!(". Reading: {}", py(r)?);
        }
        out.push(line);
    }
    Ok(out)
}

/// The text report's block, after its head.
pub fn text_lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let mut out = lines(d)?;
    out.extend(locus_lines(d)?);
    out.extend(cn_intervals::lines(d.get("cn_intervals")?)?);
    Ok(out)
}

/// The HTML card's lines after its status: the step's lines as label and
/// value, the table of the other loci with its caveat, the intervals.
pub fn html_lines(d: &Section<'_>) -> Result<(Vec<Line>, Vec<&'static str>), Malformed> {
    let mut body = Vec::new();
    let mut caveats = Vec::new();
    for line in lines(d)? {
        let (label, value) = line.split_once(": ").unwrap_or((&line, ""));
        body.push(Line::Stat(
            label.to_owned(),
            crate::model::Inline::Text(value.to_owned()),
        ));
    }
    let loci = loci(d)?;
    if !loci.is_empty() {
        let path = format!("{PATH}.loci[]");
        let mut t =
            "    <table><tr><th>Genes</th><th>Copies together (agCN)</th><th>Reading</th></tr>"
                .to_owned();
        for x in loci {
            let reading = if is_str(get(x, "status", &path)?, "called") {
                py(need(x, "reading", &path)?)?
            } else {
                ".".to_owned()
            };
            t += &format!(
                "<tr><td>{}</td><td>{}</td><td>{}</td></tr>",
                e(&genes(x)?),
                e(&py(need(x, "text", &path)?)?),
                e(&reading)
            );
        }
        t += "</table>";
        body.push(Line::Raw(t));
        caveats.push(PARASCOPY_LOCI_NOTE);
    }
    body.extend(cn_intervals::html(d.get("cn_intervals")?)?);
    Ok((body, caveats))
}
