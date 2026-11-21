//! Step 33's mosaic loss of chromosome Y or X: the call, why, and the
//! estimates from depth and from chrX allele fractions.

use pycompat::fmt::fixed;
use pycompat::json::Value;

use crate::model::{Card, Colour, HtmlPart, Inline, Line, Malformed, Section, Summary, num, py};
use crate::notes::MOSAIC_NOTE;
use crate::{html, text};

const KEY: &str = "mosaic_sex_loss";

/// The words and badge colour of each call; a call not listed reads as it
/// is, in grey.
fn call_label(call: &str) -> Option<(&'static str, Colour)> {
    match call {
        "detected" => Some(("detected", Colour::Yellow)),
        "not_detected" => Some(("not detected above 10%", Colour::Gray)),
        "gain" => Some(("an extra copy, not a loss", Colour::Yellow)),
        "not_estimated" => Some(("not estimated", Colour::Gray)),
        _ => None,
    }
}

/// `str(d.get(key))`.
fn get_str(d: &Section<'_>, key: &str) -> Result<String, Malformed> {
    d.get(key)?.map_or(Ok("None".to_owned()), py)
}

/// A share as a percentage with one decimal, or `?` when it is no number.
fn pct(v: Option<&Value>) -> String {
    num(v).map_or_else(|| "?".to_owned(), |x| format!("{}%", fixed(x * 100.0, 1)))
}

/// The estimate from depth and the one from chrX allele fractions, one line
/// each, as both reports print them.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let chrom = d.get("chrom")?;
    let what = if crate::model::is_str(chrom, "chrY") {
        "the Y"
    } else {
        "one X"
    };
    let mut out = Vec::new();
    if num(d.get("loss_depth")?).is_some() {
        let ratio = num(d.get("depth_ratio")?).map_or_else(|| "?".to_owned(), |x| fixed(x, 3));
        out.push(format!(
            "From depth: estimated share of the cells that lost {what} {} (95% range {} to {}); {} to \
             autosome depth ratio {ratio}, {} expected",
            pct(d.get("loss_depth")?),
            pct(d.get("loss_depth_low")?),
            pct(d.get("loss_depth_high")?),
            get_str(d, "chrom")?,
            get_str(d, "depth_ratio_expected")?
        ));
    }
    if num(d.get("loss_af")?).is_some() {
        out.push(format!(
            "From chrX allele fractions at {} heterozygous sites: {} (95% range {} to {})",
            get_str(d, "x_het_sites")?,
            pct(d.get("loss_af")?),
            pct(d.get("loss_af_low")?),
            pct(d.get("loss_af_high")?)
        ));
    }
    Ok(out)
}

/// `q.get(key) or default` as text.
fn or(d: &Section<'_>, key: &str, default: &str) -> Result<String, Malformed> {
    match d.get(key)? {
        Some(v) if crate::model::truthy(v) => py(v),
        _ => Ok(default.to_owned()),
    }
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, KEY, None)? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = s.section(KEY);
        let call = d.get("call")?;
        let label = match call {
            Some(Value::Str(c)) => call_label(c).map(|(l, _)| l.to_owned()),
            _ => None,
        };
        let label = match label {
            Some(l) => l,
            None => call.map_or(Ok("None".to_owned()), py)?,
        };
        w.push(format!(
            "  {}: {label} ({})",
            or(&d, "chrom", "chrY or chrX")?,
            get_str(&d, "reason")?
        ));
        w.extend(lines(&d)?.into_iter().map(|l| format!("  {l}")));
        w.push(format!("  {MOSAIC_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, shown whether or not the step ran.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section(KEY);
    let call = d.get("call")?;
    let known = match call {
        Some(Value::Str(c)) => call_label(c),
        _ => None,
    };
    let (label, colour) = match known {
        Some((l, c)) => (l.to_owned(), c),
        None => (or(&d, "call", "?")?, Colour::Gray),
    };
    let mut body = vec![
        Line::Stat(
            or(&d, "chrom", "chrY or chrX")?,
            Inline::Badge(label, colour),
        ),
        Line::Stat("Why".to_owned(), Inline::Quiet(or(&d, "reason", "")?)),
    ];
    body.extend(lines(&d)?.into_iter().map(Line::Note));
    let card = html::card(
        &d,
        Card {
            title: "Mosaic Loss of Chromosome Y or X".to_owned(),
            caveats: vec![MOSAIC_NOTE],
            body,
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
