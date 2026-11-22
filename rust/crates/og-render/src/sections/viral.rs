//! Step 47's viral DNA in the unmapped reads: the HHV-6 integration call,
//! then each viral genome with reads.
//!
//! One change from `bin/render_report.py`: when the BAM has no unplaced
//! unmapped read, step 47 assessed no genome, and the report says so
//! instead of "No read on any of the N viral genomes".

use pycompat::json::Value;

use crate::model::{
    Card, Cell, HtmlPart, Inline, Layout, Line, Malformed, Section, Summary, Table, is_str, need,
    py,
};
use crate::notes::VIRAL_NOTE;
use crate::{html, text};

const PATH: &str = "sections.viral.data.found[]";

/// The summary's HHV-6 line when the BAM had nothing to look in.
pub const NOT_ASSESSED: &str = "not assessed (the BAM has no unplaced unmapped reads)";

fn found<'a>(d: &Section<'a>) -> Result<&'a [Value], Malformed> {
    match d.get("found")? {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed("sections.viral.data.found is not a list".into())),
    }
}

/// Step 47 looked in no read: the BAM had no unplaced unmapped read.
fn none_assessed(d: &Section<'_>) -> Result<bool, Malformed> {
    Ok(is_str(d.get("hhv6")?, NOT_ASSESSED))
}

const BELOW_THRESHOLD: &str = "below threshold";

/// The integration result, scoped to the reads when threshold data is available.
pub fn hhv6_reading<'a>(reading: &'a str, found: &[Value]) -> &'a str {
    if reading == "no integrated HHV-6"
        && found.iter().any(
            |row| matches!(row, Value::Object(fields) if is_str(fields.get("call"), BELOW_THRESHOLD)),
        )
    {
        "no evidence of inherited integrated HHV-6 in these reads"
    } else {
        reading
    }
}

fn hhv6(d: &Section<'_>, found: &[Value]) -> Result<String, Malformed> {
    let reading = d.get("hhv6")?.map_or(Ok(".".to_owned()), py)?;
    Ok(hhv6_reading(&reading, found).to_owned())
}

/// Explanations for rows with reads below the reporting threshold.
pub fn threshold_notes(found: &[Value]) -> Vec<String> {
    let mut notes = Vec::new();
    for row in found {
        if let Value::Object(fields) = row
            && is_str(fields.get("call"), BELOW_THRESHOLD)
            && let Some(Value::Str(virus)) = fields.get("virus")
            && let Some(Value::Str(note)) = fields.get("note")
            && !matches!(note.trim(), "" | ".")
        {
            notes.push(format!("{virus}: {note}"));
        }
    }
    notes
}

/// `f"{float(x):.0%}"` of the breadth.
fn percent(v: &Value) -> Result<String, Malformed> {
    let x = match v {
        Value::Str(s) => pycompat::num::parse_float(s)
            .map_err(|_| Malformed(format!("{PATH}.breadth {s:?} is not a number")))?,
        v => crate::model::num(Some(v))
            .ok_or_else(|| Malformed(format!("{PATH}.breadth is not a number")))?,
    };
    Ok(pycompat::fmt::percent(x, 0))
}

/// Step 47's lines, the same in both reports: the HHV-6 call, then each
/// genome with reads, then the caveat.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let found = found(d)?;
    let mut out = vec![format!("HHV-6: {}", hhv6(d, found)?)];
    for x in found {
        let copies = need(x, "copies_per_cell", PATH)?;
        let copies = if is_str(Some(copies), ".") {
            String::new()
        } else {
            format!(", {} copies per cell", py(copies)?)
        };
        out.push(format!(
            "{}: {} ({} reads, {} of the genome{copies})",
            py(need(x, "virus", PATH)?)?,
            py(need(x, "call", PATH)?)?,
            py(need(x, "reads", PATH)?)?,
            percent(need(x, "breadth", PATH)?)?
        ));
    }
    if found.is_empty() {
        let n = d.get("genomes")?.map_or(Ok("0".to_owned()), py)?;
        out.push(if none_assessed(d)? {
            format!(
                "Not assessed on any of the {n} viral genomes: the BAM has no unplaced unmapped reads to look in."
            )
        } else {
            format!("No read on any of the {n} viral genomes.")
        });
    }
    out.extend(threshold_notes(found));
    out.push(VIRAL_NOTE.to_owned());
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "viral",
        Some("Viral DNA in the Unmapped Reads (step 47)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        w.extend(
            lines(&s.section("viral"))?
                .into_iter()
                .map(|l| format!("  {l}")),
        );
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or was meant to.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("viral");
    if d.state == crate::model::State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut card = Card {
        title: "Viral DNA in the Unmapped Reads".to_owned(),
        full: true,
        caveats: vec![VIRAL_NOTE],
        ..Card::default()
    };
    if d.state.ran() || d.state == crate::model::State::Other {
        let found = found(&d)?;
        card.body.push(Line::Stat(
            "HHV-6".to_owned(),
            Inline::Text(hhv6(&d, found)?),
        ));
        let mut rows = Vec::new();
        for x in found {
            let mut row = Vec::new();
            for k in [
                "virus",
                "call",
                "reads",
                "breadth",
                "relative_depth",
                "copies_per_cell",
            ] {
                row.push(Cell::Text(py(need(x, k, PATH)?)?));
            }
            rows.push(row);
        }
        card.body.push(if rows.is_empty() {
            Line::Note(if none_assessed(&d)? {
                "Not assessed on any viral genome: the BAM has no unplaced unmapped reads to look in."
                    .to_owned()
            } else {
                "No read on any viral genome.".to_owned()
            })
        } else {
            Line::Table(Table {
                head: [
                    "Virus",
                    "Call",
                    "Reads",
                    "Breadth",
                    "Depth / human",
                    "Copies per cell",
                ]
                .map(str::to_owned)
                .to_vec(),
                rows,
                viewer: false,
                layout: Layout::OneLine,
            })
        });
        card.body
            .extend(threshold_notes(found).into_iter().map(Line::Note));
    }
    Ok(HtmlPart {
        cards: vec![html::card(&d, card)?],
        files: Vec::new(),
    })
}
