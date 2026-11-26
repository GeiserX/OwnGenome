//! Archaic ancestry (step 43, IBDmix): the share of the analysed autosomes
//! in Neanderthal and Denisovan segments, and the longest segments.

use pycompat::fmt::fixed;
use pycompat::json::Value;

use crate::model::{Card, Inline, Line, Malformed, State, Summary, need, py};
use crate::notes::ARCHAIC_NOTE;
use crate::sections::ancestry::{get_text, has, items, need_text, number};
use crate::{html, text};

const PATH: &str = "sections.archaic.data";

/// A length in base pairs as megabases (`50.0 Mb`).
fn mb(d: &Value, key: &str) -> Result<String, Malformed> {
    let bp = number(need(d, key, PATH)?, &format!("{PATH}.{key}"))?;
    Ok(format!("{} Mb", fixed(bp / 1e6, 1)))
}

/// Step 43's lines, the same in both reports.
pub fn lines(d: &Value) -> Result<Vec<String>, Malformed> {
    let n = |k: &str| need_text(d, k, PATH);
    let mut out = Vec::new();
    for (name, key) in [("Neanderthal", "neanderthal"), ("Denisovan", "denisovan")] {
        out.push(format!(
            "{name}: {}% ({} in {} segments)",
            n(&format!("{key}_percent"))?,
            mb(d, &format!("{key}_bp"))?,
            n(&format!("{key}_segments"))?
        ));
    }
    out.push(format!(
        "Analysed: {} of the autosomes; modern frequencies: 1000 Genomes {}",
        mb(d, "analysed_bp")?,
        get_text(d, "population", ".", PATH)?
    ));
    if has(d, "longest", PATH)? {
        out.push("Longest segments:".to_owned());
        for x in items(need(d, "longest", PATH)?, &format!("{PATH}.longest"))? {
            out.push(format!("  {}", py(&x)?));
        }
    }
    Ok(out)
}

/// The text report's block. It ends without a blank line, as the Python
/// report's did.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "archaic", Some("Archaic Ancestry (IBDmix, step 43)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = s.section("archaic").data;
        w.extend(lines(d)?.iter().map(|l| format!("  {l}")));
        w.push(format!("  {ARCHAIC_NOTE}"));
    }
    Ok(w)
}

/// The HTML report's card, when step 43 has a result.
pub fn html(s: &Summary) -> Result<crate::model::HtmlPart, Malformed> {
    let sec = s.section("archaic");
    if sec.state == State::Missing {
        return Ok(crate::model::HtmlPart::default());
    }
    let mut c = Card {
        title: "Archaic Ancestry (Neanderthal, Denisovan)".to_owned(),
        ..Card::default()
    };
    if sec.state.ran() {
        let d = sec.data;
        for (label, key) in [
            ("Neanderthal", "neanderthal_percent"),
            ("Denisovan", "denisovan_percent"),
        ] {
            c.body.push(Line::Stat(
                label.to_owned(),
                Inline::Text(format!("{}%", need_text(d, key, PATH)?)),
            ));
        }
        c.body.extend(lines(d)?.into_iter().map(Line::Para));
        c.caveats.push(ARCHAIC_NOTE);
    }
    Ok(crate::model::HtmlPart {
        cards: vec![html::card(&sec, c)?],
        files: Vec::new(),
    })
}
