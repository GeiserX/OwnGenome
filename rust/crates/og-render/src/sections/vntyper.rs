//! Step 52's MUC1 VNTR frameshifts (VNtyper): the result, the call, the
//! depth with its sensitivity, and the advice.

use crate::model::{
    Card, Colour, HtmlPart, Inline, Line, Malformed, Section, State, Summary, is_str, py,
};
use crate::notes::VNTYPER_NOTE;
use crate::{html, text};

/// `d.get(key, ".")` as text.
fn dot(d: &Section<'_>, key: &str) -> Result<String, Malformed> {
    d.get(key)?.map_or(Ok(".".to_owned()), py)
}

/// Step 52's reading, the lines of the text section and of the card.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let mut out = vec![format!("Result: {}", dot(d, "status")?), dot(d, "message")?];
    let variant = d.get("variant")?;
    if variant.is_some() && !is_str(variant, ".") {
        out.push(format!(
            "Variant: {} (Kestrel confidence {}, {} k-mer paths for the variant)",
            dot(d, "variant")?,
            dot(d, "confidence")?,
            dot(d, "alt_depth")?
        ));
    }
    out.push(format!(
        "VNtyper's coverage check: {}",
        dot(d, "coverage_qc")?
    ));
    out.push(dot(d, "sensitivity")?);
    out.push(dot(d, "advice")?);
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "vntyper",
        Some("MUC1 VNTR Frameshifts (ADTKD-MUC1, VNtyper)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        w.extend(
            lines(&s.section("vntyper"))?
                .into_iter()
                .map(|l| format!("  {l}")),
        );
        w.push(format!("  {VNTYPER_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or was meant to.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("vntyper");
    if d.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut card = Card {
        title: "MUC1 VNTR Frameshifts (VNtyper)".to_owned(),
        caveats: vec![VNTYPER_NOTE],
        ..Card::default()
    };
    if d.state.ran() || d.state == State::Other {
        let status = d.get("status")?;
        let colour = if is_str(status, "positive") {
            Colour::Red
        } else if is_str(status, "possible") {
            Colour::Yellow
        } else if is_str(status, "negative") {
            Colour::Green
        } else {
            Colour::Gray
        };
        card.body.push(Line::Stat(
            "Result".to_owned(),
            Inline::Badge(dot(&d, "status")?, colour),
        ));
        card.body.extend(lines(&d)?.into_iter().map(Line::Para));
    }
    Ok(HtmlPart {
        cards: vec![html::card(&d, card)?],
        files: Vec::new(),
    })
}
