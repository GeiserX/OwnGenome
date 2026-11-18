//! Step 31's variant prioritization (slivar): the prioritized and
//! compound-het counts, and the phase of each compound-het pair read from
//! the reads (trans, cis or unknown), with the pairs in trans.

use pycompat::fmt::ljust;
use pycompat::json::Value;

use crate::escape::e;
use crate::model::{Card, HtmlPart, Inline, Line, Malformed, Summary};
use crate::notes::{COMPHET_NOTE, VIEWER_NOTE};
use crate::sections::clinvar::{D, VIEWER_TH, list_or_empty, locus_text, viewer_td_of};
use crate::{html, text};

const TITLE: &str = "Variant Prioritization (slivar)";

/// `ph.get(key, 0)` as text.
fn count(ph: &D<'_>, key: &str) -> Result<String, Malformed> {
    ph.s_or(key, "0")
}

/// `ph.get("trans", 0) - len(trans)`, an int or a float as Python gives it.
fn more_in_trans(ph: &D<'_>, shown: usize) -> Result<Option<String>, Malformed> {
    let Some(t) = ph.get("trans")? else {
        return Ok(None);
    };
    let n = t.num()?;
    if n <= shown as f64 {
        return Ok(None);
    }
    Ok(Some(match t.v {
        Value::Int(i) => (i - shown as i128).to_string(),
        Value::Bool(_) => (n as i128 - shown as i128).to_string(),
        _ => pycompat::fmt::repr(n - shown as f64),
    }))
}

/// The phase lines of the text block.
pub fn phase_text(d: &D<'_>) -> Result<Vec<String>, Malformed> {
    let ph = d.at("compound_het_phase")?;
    let mut out = vec![format!(
        "  Phase from the reads: {} trans, {} cis, {} unknown of {} pairs",
        count(&ph, "trans")?,
        count(&ph, "cis")?,
        count(&ph, "unknown")?,
        d.s_or("compound_het_pairs", "0")?
    )];
    let trans = list_or_empty(d, "compound_het_trans")?;
    for x in &trans {
        out.push(format!(
            "    {} {}  in trans on {} reads or read pairs",
            ljust(&x.s("gene")?, 8),
            x.s("variants")?,
            x.s("reads")?
        ));
    }
    if let Some(more) = more_in_trans(&ph, trans.len())? {
        // The Python takes len() of the list here, which a missing or
        // null list cannot give.
        d.at("compound_het_trans")?.len()?;
        out.push(format!("    ... {more} more in the pairs table"));
    }
    out.push(format!("  ({COMPHET_NOTE})"));
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "slivar", Some(TITLE))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = D::of(&s.section("slivar"));
    w.push(format!("  Prioritized variants: {}", d.s("prioritized")?));
    if d.has("frequency_unknown") {
        w.push(format!(
            "  Frequency unknown (prioritized, not counted as rare): {}",
            d.s("frequency_unknown")?
        ));
    }
    if d.has("compound_het_variants") {
        w.push(format!(
            "  Compound het candidates: {} variants across {} genes",
            d.s("compound_het_variants")?,
            d.s("compound_het_genes")?
        ));
    }
    if d.has("compound_het_phase") {
        w.extend(phase_text(&d)?);
    } else if d.has("compound_het_variants") {
        w.push(
            "  (Candidates only: unphased single-sample data. See docs/31-slivar.md)".to_owned(),
        );
    }
    w.push(String::new());
    Ok(w)
}

/// The phase lines of the card: the counts, then the pairs in trans.
fn phase_html(d: &D<'_>) -> Result<Vec<Line>, Malformed> {
    let ph = d.at("compound_het_phase")?;
    let mut out = vec![Line::Stat(
        "Pairs in trans / cis / unknown".to_owned(),
        Inline::Text(format!(
            "{} / {} / {}",
            count(&ph, "trans")?,
            count(&ph, "cis")?,
            count(&ph, "unknown")?
        )),
    )];
    let trans = list_or_empty(d, "compound_het_trans")?;
    if !trans.is_empty() {
        out.push(Line::Raw("    <table>".to_owned()));
        out.push(Line::Raw(format!(
            "      <tr><th>Gene</th><th>Pair in trans</th><th>Reads</th>{VIEWER_TH}</tr>"
        )));
        for x in &trans {
            let variants = x.at("variants")?;
            out.push(Line::Raw(format!(
                "      <tr><td>{}</td><td>{}</td><td>{}</td>{}</tr>",
                e(&x.s("gene")?),
                e(&variants.text()?),
                e(&x.s("reads")?),
                viewer_td_of(&locus_text(&variants)?)?
            )));
        }
        out.push(Line::Raw("    </table>".to_owned()));
        out.push(Line::Raw(VIEWER_NOTE.to_owned()));
    }
    Ok(out)
}

/// The HTML report's card, always there.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("slivar");
    let d = D::of(&sec);
    let mut body = vec![Line::Stat(
        "Prioritized variants".to_owned(),
        Inline::Text(d.s_or_none("prioritized")?),
    )];
    if d.has("frequency_unknown") {
        body.push(Line::Stat(
            "Frequency unknown (not rare)".to_owned(),
            Inline::Text(d.s("frequency_unknown")?),
        ));
    }
    if d.has("compound_het_variants") {
        body.push(Line::Stat(
            "Compound het candidates".to_owned(),
            Inline::Text(format!(
                "{} variants, {} genes",
                d.s("compound_het_variants")?,
                d.s("compound_het_genes")?
            )),
        ));
    }
    let mut caveats = Vec::new();
    if d.has("compound_het_phase") {
        body.extend(phase_html(&d)?);
        caveats.push(COMPHET_NOTE);
    }
    let card = html::card(
        &sec,
        Card {
            title: TITLE.to_owned(),
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
