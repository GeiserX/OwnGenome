//! Step 10's telomere content (TelomereHunter): the relative content, the
//! intratelomeric reads it rests on, and a warning when they are too few for
//! a stable estimate.

use pycompat::json::Value;

use crate::model::{Card, Colour, HtmlPart, Inline, Line, Malformed, Section, Summary, py};
use crate::notes::TELOMERE_NOTE;
use crate::{html, text};

/// The unit of the telomere content.
pub const TELOMERE_UNIT: &str = "intratelomeric reads per million reads with 48-52% GC";

const WARNING: &str = "WARNING: ";

/// `f"{n:,}"` of a whole number.
fn grouped(v: &Value, path: &str) -> Result<String, Malformed> {
    match v {
        Value::Int(n) => Ok(pycompat::fmt::int_grouped(*n)),
        Value::Bool(b) => Ok(u8::from(*b).to_string()),
        _ => Err(Malformed(format!("{path} is not a whole number"))),
    }
}

/// The read count under the telomere value, and the warning when it is too
/// thin, the same in both reports.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let path = "sections.telomere.data";
    let check = d.get("read_check")?;
    let Some(reads) = d.get("intratelomeric_reads")? else {
        return match check.filter(|v| crate::model::truthy(v)) {
            Some(c) => Ok(vec![format!(
                "Intratelomeric reads: not checked ({})",
                py(c)?
            )]),
            None => Ok(Vec::new()),
        };
    };
    let n = grouped(reads, &format!("{path}.intratelomeric_reads"))?;
    let Some(floor) = d.get("read_floor")? else {
        let check = check.map_or(Ok("not checked".to_owned()), py)?;
        return Ok(vec![format!("Intratelomeric reads: {n}; {check}")]);
    };
    let mut out = vec![format!(
        "Intratelomeric reads: {n} (floor {} at {}x)",
        grouped(floor, &format!("{path}.read_floor"))?,
        py(d.need("depth")?)?
    )];
    if d.get("low_reads")?.is_some_and(crate::model::truthy) {
        out.push(format!(
            "{WARNING}too few telomeric reads for a stable estimate; do not compare this value."
        ));
    }
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "telomere",
        Some("Telomere content (relative, TelomereHunter)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = s.section("telomere");
        w.push(format!(
            "  Telomere content: {} {TELOMERE_UNIT}",
            py(d.need("tel_content")?)?
        ));
        w.extend(lines(&d)?.into_iter().map(|l| format!("  {l}")));
        w.push(format!("  {TELOMERE_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, shown whether or not the step ran.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("telomere");
    let mut body = vec![Line::Stat(
        "Telomere content".to_owned(),
        Inline::Text(d.get("tel_content")?.map_or(Ok(".".to_owned()), py)?),
    )];
    let mut caveats = Vec::new();
    if d.state.ran() || d.state == crate::model::State::Other {
        for l in lines(&d)? {
            body.push(match l.strip_prefix(WARNING) {
                Some(w) => Line::Stat(
                    "Warning".to_owned(),
                    Inline::Badge(w.to_owned(), Colour::Yellow),
                ),
                None => Line::Note(l),
            });
        }
        caveats.push(TELOMERE_NOTE);
    }
    let card = html::card(
        &d,
        Card {
            title: "Telomere content (relative)".to_owned(),
            caveats,
            body,
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
