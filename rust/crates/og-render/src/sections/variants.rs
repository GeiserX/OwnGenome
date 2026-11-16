//! Step 03's variant calls: how many DeepVariant made against the
//! reference, and how many passed its filters.

use crate::model::{Card, HtmlPart, Inline, Line, Malformed, Summary, py};
use crate::notes::VARIANTS_NOTE;
use crate::{html, text};

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "variants", Some("Variant Calling (DeepVariant)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let sec = s.section("variants");
        let n = |k: &str| sec.need(k).and_then(py);
        w.push(format!("  Total variants: {}", n("total")?));
        w.push(format!("  PASS variants:  {}", n("pass")?));
        w.push(format!("  SNPs: {}  Indels: {}", n("snps")?, n("indels")?));
        w.push(format!("  {VARIANTS_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card. A count the summary lacks prints `None`, as the
/// Python did.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("variants");
    let mut body = Vec::new();
    for (label, k) in [
        ("Total variants", "total"),
        ("PASS variants", "pass"),
        ("SNPs", "snps"),
        ("Indels", "indels"),
    ] {
        let v = sec.get(k)?.map_or(Ok("None".to_owned()), py)?;
        body.push(Line::Stat(label.to_owned(), Inline::Text(v)));
    }
    let card = html::card(
        &sec,
        Card {
            title: "Variant Calling".to_owned(),
            caveats: vec![VARIANTS_NOTE],
            body,
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
