//! Step 49's clonal haematopoiesis (CHIP): the calls, the detection limit
//! at this depth, and the genes too shallow to read.

use pycompat::json::Value;

use crate::model::{
    Card, HtmlPart, Inline, Line, Malformed, Section, State, Summary, need, num, py, truthy,
};
use crate::notes::CHIP_NOTE;
use crate::{html, text};

const PATH: &str = "sections.chip.data";

fn calls<'a>(d: &Section<'a>) -> Result<&'a [Value], Malformed> {
    match d.get("calls")? {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!("{PATH}.calls is not a list"))),
    }
}

/// `f"{x:.0f}"` of a number of the summary.
fn whole(v: &Value, key: &str) -> Result<String, Malformed> {
    num(Some(v))
        .map(|x| pycompat::fmt::fixed(x, 0))
        .ok_or_else(|| Malformed(format!("{PATH}.{key} is not a number")))
}

/// The detection limit sentence: the VAF found 8 times in 10 at this depth,
/// or that no gene has the depth to find a clone.
pub fn limit_line(d: &Section<'_>) -> Result<String, Malformed> {
    match d.get("median_limit")?.filter(|v| !matches!(v, Value::Null)) {
        Some(l) => {
            let pct = num(Some(l))
                .ok_or_else(|| Malformed(format!("{PATH}.median_limit is not a number")))?;
            Ok(format!(
                "Detection limit: a clone of VAF about {}% or more is found 8 times in 10 at this depth \
                 (median {}x over the genes' coding exons). CHIP starts at a VAF of 2%: smaller clones are \
                 likely missed, so no call does not rule one out.",
                pycompat::fmt::fixed(pct * 100.0, 0),
                whole(d.need("median_depth")?, "median_depth")?
            ))
        }
        None => Ok("No CHIP gene has the depth to find a clone: no result either way.".to_owned()),
    }
}

/// The CHIP section's lines, the same in both reports, without the caveat:
/// the calls, the detection limit, the genes too shallow to read and the
/// table of the other changes.
pub fn lines(d: &Section<'_>, sample: &str) -> Result<Vec<String>, Malformed> {
    let mut out = Vec::new();
    let calls = calls(d)?;
    for c in calls {
        let p = format!("{PATH}.calls[]");
        let f = |k: &str| need(c, k, &p).and_then(py);
        out.push(format!(
            "{} {}: VAF {} ({} of {} reads), {}",
            f("gene")?,
            f("protein_change")?,
            f("vaf")?,
            f("alt_reads")?,
            f("depth")?,
            f("note")?
        ));
    }
    if calls.is_empty() {
        out.push("No CHIP call.".to_owned());
    }
    out.push(limit_line(d)?);
    if let Some(nc) = d.get("not_covered")?.filter(|v| truthy(v)) {
        let Value::Array(genes) = nc else {
            return Err(Malformed(format!("{PATH}.not_covered is not a list")));
        };
        let genes = genes.iter().map(py).collect::<Result<Vec<_>, _>>()?;
        out.push(format!(
            "Too shallow or not on the reference: {}.",
            genes.join(", ")
        ));
    }
    out.push(format!(
        "{} other protein changes listed, not counted: chip/{sample}_chip.tsv.",
        d.get("other")?.map_or(Ok("0".to_owned()), py)?
    ));
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "chip", Some("Clonal Haematopoiesis (CHIP, step 49)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let sample = py(s.need("sample")?)?;
        w.extend(
            lines(&s.section("chip"), &sample)?
                .into_iter()
                .map(|l| format!("  {l}")),
        );
        w.push(format!("  {CHIP_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or was meant to.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("chip");
    if d.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut card = Card {
        title: "Clonal Haematopoiesis (CHIP)".to_owned(),
        full: true,
        caveats: vec![CHIP_NOTE],
        ..Card::default()
    };
    if d.state.ran() || d.state == State::Other {
        let sample = py(s.need("sample")?)?;
        card.body.push(Line::Stat(
            "CHIP calls".to_owned(),
            Inline::Text(calls(&d)?.len().to_string()),
        ));
        card.body
            .extend(lines(&d, &sample)?.into_iter().map(Line::Para));
    }
    Ok(HtmlPart {
        cards: vec![html::card(&d, card)?],
        files: Vec::new(),
    })
}
