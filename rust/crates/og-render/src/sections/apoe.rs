//! APOE, the Alzheimer's risk factor: the report keeps its lines apart and
//! last, behind a warning, because many people choose not to learn them.
//! The text report prints them after the warning; the HTML report folds
//! them away until opened.

use pycompat::json::Value;

use crate::model::{Cell, Layout, Line, Malformed, Table, need, py};

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

/// The HTML report's APOE table, folded away behind the warning: `head` is
/// the traits table's column titles, `row` builds one row as the traits
/// table does. None without APOE rows.
pub fn html_line(
    rows: &[Value],
    head: Vec<String>,
    row: impl Fn(&Value) -> Result<Vec<Cell>, Malformed>,
) -> Result<Option<Line>, Malformed> {
    if rows.is_empty() {
        return Ok(None);
    }
    Ok(Some(Line::Folded(
        format!("{APOE_WARNING} Open to show it."),
        Table {
            head,
            rows: rows.iter().map(row).collect::<Result<_, _>>()?,
            layout: Layout::Indented,
            ..Table::default()
        },
    )))
}
