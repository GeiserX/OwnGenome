//! The haplogroups: the mitochondrial one with haplocheck's contamination
//! verdict (step 12), and the Y chromosome section, which shows Yleaf's
//! haplogroup (step 37), the second placement on the current public tree
//! ([`crate::sections::y_tree`]) and the Y-STR markers
//! ([`crate::sections::ystr`]). Each haplogroup carries the links where it
//! can be explored further.

use pycompat::json::Value;

use crate::model::{Card, Colour, HtmlPart, Inline, Line, Malformed, Section, State, Summary, py};
use crate::notes::EXPLORE_NOTE;
use crate::sections::ancestry::{explore_html, explore_lines, get_text, need_text};
use crate::sections::{y_tree, ystr};
use crate::{html, text};

/// The Status badge of a section that ran: Stale for an earlier run's
/// result, else Complete.
pub fn done_badge(sec: &Section<'_>) -> Inline {
    if sec.state == State::Stale {
        Inline::Badge("Stale".to_owned(), Colour::Yellow)
    } else {
        Inline::Badge("Complete".to_owned(), Colour::Green)
    }
}

/// Whether a card's lines are shown: the card says Not run, Failed or
/// Unreadable in their place otherwise.
pub fn shows_body(sec: &Section<'_>) -> bool {
    !matches!(
        sec.state,
        State::Missing | State::Failed | State::Unreadable
    )
}

fn path(sec: &Section<'_>) -> String {
    format!("sections.{}.data", sec.key)
}

/// haplocheck's verdict: whether the reads hold two mtDNA haplogroups.
pub fn contamination_text(sec: &Section<'_>) -> Result<String, Malformed> {
    let Some(status) = sec.get("contamination_status")? else {
        return Ok("not checked (needs step 20's Mutect2 chrM calls)".to_owned());
    };
    let word = match status {
        Value::Str(s) => match s.as_str() {
            "YES" => "YES, two mtDNA haplogroups in the reads (another person's DNA?)".to_owned(),
            "NO" => "no".to_owned(),
            "ND" => "not determined".to_owned(),
            _ => s.clone(),
        },
        Value::Array(_) | Value::Object(_) => {
            return Err(Malformed(format!(
                "{}.contamination_status is not a value",
                path(sec)
            )));
        }
        v => py(v)?,
    };
    let level = get_text(sec.data, "contamination_level", ".", &path(sec))?;
    Ok(format!("{word} (level {level})"))
}

/// The text report's mitochondrial block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "haplogroup", Some("Mitochondrial Haplogroup"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let sec = s.section("haplogroup");
        let p = path(&sec);
        w.push(format!(
            "  Haplogroup: {}",
            need_text(sec.data, "haplogroup", &p)?
        ));
        w.push(format!(
            "  Contamination (haplocheck): {}",
            contamination_text(&sec)?
        ));
        w.extend(
            explore_lines(sec.data, &p)?
                .iter()
                .map(|l| format!("  {l}")),
        );
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's mitochondrial card. It is always there, so a run
/// without step 12 says so.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("haplogroup");
    let mut c = Card {
        title: "Mitochondrial Haplogroup".to_owned(),
        ..Card::default()
    };
    if shows_body(&sec) {
        let p = path(&sec);
        c.body = vec![
            Line::Stat("Status".to_owned(), done_badge(&sec)),
            Line::Stat(
                "Haplogroup".to_owned(),
                Inline::Text(get_text(sec.data, "haplogroup", ".", &p)?),
            ),
            Line::Stat(
                "Contamination (haplocheck)".to_owned(),
                Inline::Text(contamination_text(&sec)?),
            ),
        ];
        let (links, note) = explore_html(sec.data, &p)?;
        c.body.extend(links);
        if note {
            c.caveats.push(EXPLORE_NOTE);
        }
    }
    Ok(HtmlPart {
        cards: vec![html::card(&sec, c)?],
        files: Vec::new(),
    })
}

/// The text report's block of Yleaf's haplogroup.
fn y_text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "y_haplogroup", Some("Y-Chromosome Haplogroup (Yleaf)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let sec = s.section("y_haplogroup");
        let p = path(&sec);
        let n = |k: &str| need_text(sec.data, k, &p);
        w.push(format!(
            "  Haplogroup: {} ({} markers, QC-score {})",
            n("haplogroup")?,
            n("valid_markers")?,
            n("qc_score")?
        ));
        w.extend(
            explore_lines(sec.data, &p)?
                .iter()
                .map(|l| format!("  {l}")),
        );
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card of Yleaf's haplogroup, when step 37 has a result.
fn y_card(s: &Summary) -> Result<Option<Card>, Malformed> {
    let sec = s.section("y_haplogroup");
    if sec.state == State::Missing {
        return Ok(None);
    }
    let mut c = Card {
        title: "Y-Chromosome Haplogroup".to_owned(),
        ..Card::default()
    };
    if shows_body(&sec) {
        let p = path(&sec);
        for (label, key) in [
            ("Haplogroup", "haplogroup"),
            ("Markers", "valid_markers"),
            ("QC-score", "qc_score"),
        ] {
            c.body.push(Line::Stat(
                label.to_owned(),
                Inline::Text(get_text(sec.data, key, ".", &p)?),
            ));
        }
        let (links, note) = explore_html(sec.data, &p)?;
        c.body.extend(links);
        if note {
            c.caveats.push(EXPLORE_NOTE);
        }
    }
    Ok(Some(html::card(&sec, c)?))
}

/// The text report's Y chromosome blocks: Yleaf, the second placement,
/// the Y-STR markers.
pub fn y_chromosome_text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let mut w = y_text(s)?;
    w.extend(y_tree::text(s)?);
    w.extend(ystr::text(s)?);
    Ok(w)
}

/// The HTML report's Y chromosome cards, each when its step has a result.
pub fn y_chromosome_html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let cards = [y_card(s)?, y_tree::card(s)?, ystr::card(s)?]
        .into_iter()
        .flatten()
        .collect();
    Ok(HtmlPart {
        cards,
        files: Vec::new(),
    })
}
