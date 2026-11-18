//! APOE, the Alzheimer's risk factor: the report keeps its lines apart and
//! last, behind a warning, because many people choose not to learn them.
//! The text report prints them after the warning; the HTML report folds
//! them away until opened.

use pycompat::json::Value;

use crate::model::{Malformed, need, py};

/// The warning before the APOE lines.
pub const APOE_WARNING: &str = "APOE and Alzheimer's disease: a risk factor, not a diagnosis. Many \
    people choose not to learn their APOE result; docs/42-traits.md explains it.";

/// The text report's APOE lines, the warning first; none without APOE rows.
pub fn text_lines(rows: &[Value]) -> Result<Vec<String>, Malformed> {
    if rows.is_empty() {
        return Ok(Vec::new());
    }
    let mut out = vec![format!("  {APOE_WARNING}")];
    for r in rows {
        let f = |k: &str| need(r, k, "sections.traits.data.apoe[]").and_then(py);
        out.push(format!(
            "    {} ({}): {}. {}",
            f("trait")?,
            f("rsid")?,
            f("genotype")?,
            f("reading")?
        ));
    }
    Ok(out)
}

/// The HTML report's folded APOE table: `head` is the traits table's
/// header row, `row` writes one row as the traits table does.
pub fn html_lines(
    rows: &[Value],
    head: &str,
    row: impl Fn(&Value) -> Result<String, Malformed>,
) -> Result<Vec<String>, Malformed> {
    if rows.is_empty() {
        return Ok(Vec::new());
    }
    let mut out = vec![
        format!(
            "    <details><summary>{} Open to show it.</summary>",
            crate::escape::e(APOE_WARNING)
        ),
        "    <table>".to_owned(),
        head.to_owned(),
    ];
    for r in rows {
        out.push(format!("      {}", row(r)?));
    }
    out.push("    </table></details>".to_owned());
    Ok(out)
}
