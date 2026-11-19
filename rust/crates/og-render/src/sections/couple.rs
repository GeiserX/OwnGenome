//! Step 39's couple carrier screen: the counts, then each shared and
//! X-linked gene with the chance for a child and both partners' findings.

use crate::model::{Card, HtmlPart, Line, Malformed, State, Summary};
use crate::notes::COUPLE_NOTE;
use crate::sections::clinvar::D;
use crate::{html, text};

/// The lines, the same in both reports; the last is [`COUPLE_NOTE`].
pub fn lines(d: &D<'_>) -> Result<Vec<String>, Malformed> {
    let partner = match d.opt("partner")? {
        Some(p) => p.text()?,
        None => "the partner".to_owned(),
    };
    let count = |k: &str| -> Result<usize, Malformed> { d.get(k)?.map_or(Ok(0), |v| v.len()) };
    let mut out = vec![format!(
        "Partner: {partner}. Genes where both carry a change: {}; X-linked findings: {}; genes where only \
         one partner carries a change: {}",
        count("shared")?,
        count("x_linked")?,
        d.s_or("one_partner", "0")?
    )];
    for (label, key) in [("Shared", "shared"), ("X-linked", "x_linked")] {
        let Some(list) = d.get(key)? else { continue };
        for x in list.items()? {
            out.push(format!(
                "{label}: {} ({}): {}. This sample: {}. {partner}: {}",
                x.s("gene")?,
                x.s("condition")?,
                x.s("risk")?,
                x.s("own")?,
                x.s("partner_findings")?
            ));
        }
    }
    out.push(COUPLE_NOTE.to_owned());
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "couple", Some("Couple Carrier Screen (step 39)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = D::of(&s.section("couple"));
        w.extend(lines(&d)?.iter().map(|l| format!("  {l}")));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the screen ran or its step failed.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("couple");
    if sec.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut body = lines(&D::of(&sec))?;
    body.pop();
    let card = html::card(
        &sec,
        Card {
            title: "Couple Carrier Screen".to_owned(),
            full: true,
            caveats: vec![COUPLE_NOTE],
            body: body.into_iter().map(Line::Para).collect(),
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
