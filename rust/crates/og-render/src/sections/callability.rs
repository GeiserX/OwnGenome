//! Step 16b's gene callability: how many of the hard-to-read genes of
//! `assets/callability_genes.tsv` the reads covered, and which were not
//! read, apart from the ones a step of their own calls.

use pycompat::fmt::fixed;

use crate::model::{Card, Colour, HtmlPart, Inline, Line, Malformed, State, Summary};
use crate::notes::CALLABILITY_NOTE;
use crate::sections::clinvar::D;
use crate::{html, text};

const TITLE: &str = "Genes This Genome Could Not Read Reliably";

/// `GENE (37%)`, or `GENE (not on this reference)`.
fn gene(g: &D<'_>) -> Result<String, Malformed> {
    if g.at("status")?.is("not on this reference") {
        return Ok(format!("{} (not on this reference)", g.s("gene")?));
    }
    Ok(format!(
        "{} ({}%)",
        g.s("gene")?,
        fixed(100.0 * g.at("share")?.num()?, 0)
    ))
}

/// The lines of the text block and the card.
pub fn lines(d: &D<'_>) -> Result<Vec<String>, Malformed> {
    let mut out = vec![format!(
        "{} of the {} hard-to-read genes of assets/callability_genes.tsv were read: {}% or more of their \
         bases had {} or more reads of mapping quality {} or more.",
        d.s("read")?,
        d.s("catalog_genes")?,
        fixed(100.0 * d.at("min_share")?.num()?, 0),
        d.s("min_depth")?,
        d.s("min_mapq")?
    )];
    let nr = d.at("not_read")?.items()?;
    let owned = d.at("not_read_owned")?.items()?;
    let listed = if nr.is_empty() {
        "none".to_owned()
    } else {
        nr.iter()
            .map(gene)
            .collect::<Result<Vec<_>, _>>()?
            .join(", ")
    };
    out.push(format!(
        "Not read, and no other step calls them ({}): {listed}",
        nr.len()
    ));
    if !owned.is_empty() {
        let mut parts = Vec::new();
        for g in &owned {
            parts.push(format!("{} by {}", gene(g)?, g.s("owner")?));
        }
        out.push(format!(
            "Not read here, called by a step of their own ({}): {}",
            owned.len(),
            parts.join(", ")
        ));
    }
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "callability", Some(TITLE))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = D::of(&s.section("callability"));
        w.extend(lines(&d)?.iter().map(|l| format!("  {l}")));
        w.push(format!("  {CALLABILITY_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or failed.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("callability");
    if sec.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let d = D::of(&sec);
    let mut body = Vec::new();
    if sec.state.ran() {
        let not_read = d.at("not_read")?.len()?;
        body.push(Line::Stat(
            "Genes read".to_owned(),
            Inline::Text(format!("{} of {}", d.s("read")?, d.s("catalog_genes")?)),
        ));
        body.push(Line::Stat(
            "Not read, no other caller".to_owned(),
            Inline::Badge(
                not_read.to_string(),
                if not_read > 0 {
                    Colour::Red
                } else {
                    Colour::Green
                },
            ),
        ));
        body.extend(lines(&d)?.into_iter().map(Line::Note));
    }
    let card = html::card(
        &sec,
        Card {
            title: TITLE.to_owned(),
            caveats: vec![CALLABILITY_NOTE],
            body,
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
