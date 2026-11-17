//! The section "Paralog genes and SMN1": Paraphase (step 35b, PacBio HiFi)
//! beside Parascopy and the depth counts of step 35 (short reads), and a
//! "not checked" card when neither left a result, so a report without them
//! never reads as normal.

use pycompat::json::Value;

use crate::escape::e;
use crate::html;
use crate::model::{
    Card, Colour, HtmlPart, Inline, Line, Malformed, Section, State, Summary, get, is_str, need,
    py, truthy,
};
use crate::notes::{PARALOG_NOTE, PARASCOPY_NOTE};
use crate::sections::{parascopy, smn1::smn1_text};
use crate::text;

const PATH: &str = "sections.paralogs.data";

/// The genes Paraphase called (`d.get("genes", [])`).
fn genes<'a>(d: &Section<'a>) -> Result<&'a [Value], Malformed> {
    match d.get("genes")? {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!("{PATH}.genes is not a list"))),
    }
}

/// `sum(d.get("variants", {}).values())`, as Python prints the sum.
fn variants(d: &Section<'_>) -> Result<String, Malformed> {
    let Some(v) = d.get("variants")? else {
        return Ok("0".to_owned());
    };
    let Value::Object(o) = v else {
        return Err(Malformed(format!("{PATH}.variants is not an object")));
    };
    let (mut whole, mut float, mut is_float) = (0_i128, 0.0_f64, false);
    for x in o.values() {
        match x {
            Value::Int(i) => {
                whole += *i;
                float += *i as f64;
            }
            Value::Bool(b) => {
                whole += i128::from(*b);
                float += f64::from(u8::from(*b));
            }
            Value::Float(f) => {
                is_float = true;
                float += f;
            }
            _ => {
                return Err(Malformed(format!(
                    "{PATH}.variants holds a value that is not a number"
                )));
            }
        }
    }
    Ok(if is_float {
        pycompat::fmt::repr(float)
    } else {
        whole.to_string()
    })
}

/// The text report's lines of Paraphase's result.
fn paraphase_lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let path = format!("{PATH}.genes[]");
    let mut out = Vec::new();
    for g in genes(d)? {
        let f = |k: &str| need(g, k, &path).and_then(py);
        let detail = need(g, "detail", &path)?;
        let tail = if is_str(Some(detail), "") || is_str(Some(detail), ".") {
            String::new()
        } else {
            format!("; {}", py(detail)?)
        };
        out.push(format!(
            "{}: {} ({}){tail}",
            f("gene")?,
            f("copy_number")?,
            f("status")?
        ));
    }
    if let Some(smn) = smn1(d)? {
        out.push(format!("SMN1 by both tools: {smn}"));
    }
    if let Some(n) = d.get("other_regions")?.filter(|v| truthy(v)) {
        out.push(format!(
            "Other Paraphase regions: {} ({} with a call)",
            py(n)?,
            py(d.need("other_regions_called")?)?
        ));
    }
    out.push(format!("Variants on phased haplotypes: {}", variants(d)?));
    Ok(out)
}

/// `smn1_text(d.get("smn1") or {})`.
fn smn1(d: &Section<'_>) -> Result<Option<String>, Malformed> {
    match d.get("smn1")?.filter(|v| truthy(v)) {
        Some(m) => smn1_text(m),
        None => Ok(None),
    }
}

/// A key's state as the summary writes it.
fn raw_state(s: &Summary, key: &str) -> String {
    let sec = s.section(key);
    match sec.state {
        State::Ok => "ok".to_owned(),
        State::Stale => "stale".to_owned(),
        State::Failed => "failed".to_owned(),
        State::Unreadable => "unreadable".to_owned(),
        State::Missing => "missing".to_owned(),
        State::Other => s
            .get("sections")
            .and_then(|x| get(x, key, "sections").ok().flatten())
            .and_then(|x| get(x, "state", "sections").ok().flatten())
            .and_then(|v| py(v).ok())
            .unwrap_or_default(),
    }
}

/// `Paralog copy number: not checked (...)` when neither step 35 nor step
/// 35b left a result.
pub fn unchecked(s: &Summary) -> Option<String> {
    let states = [
        ("parascopy", "35", s.section("parascopy").state),
        ("paralogs", "35b", s.section("paralogs").state),
    ];
    if states.iter().any(|(_, _, st)| st.ran()) {
        return None;
    }
    let why = if states.iter().all(|(_, _, st)| *st == State::Missing) {
        "neither step 35 nor step 35b ran".to_owned()
    } else {
        states
            .iter()
            .map(|(k, step, st)| {
                let word = match st {
                    State::Missing => "not run".to_owned(),
                    State::Failed => "failed".to_owned(),
                    State::Unreadable => "unreadable".to_owned(),
                    _ => raw_state(s, k),
                };
                format!("step {step} {word}")
            })
            .collect::<Vec<_>>()
            .join(", ")
    };
    Some(format!("Paralog copy number: not checked ({why})"))
}

/// The text report's block: Paraphase, Parascopy, then the "not checked"
/// block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let mut w = Vec::new();
    if let Some(head) = text::head(
        s,
        "paralogs",
        Some("Paralog Genes (Paraphase, PacBio HiFi)"),
    )? {
        w.extend(head.lines);
        if head.ran {
            let d = s.section("paralogs");
            w.extend(paraphase_lines(&d)?.iter().map(|l| format!("  {l}")));
            w.push(format!("  {PARALOG_NOTE}"));
            w.push(String::new());
        }
    }
    if let Some(head) = text::head(
        s,
        "parascopy",
        Some("SMN1/SMN2 Copy Number (Parascopy, short reads)"),
    )? {
        w.extend(head.lines);
        if head.ran {
            let d = s.section("parascopy");
            w.extend(parascopy::text_lines(&d)?.iter().map(|l| format!("  {l}")));
            w.push(format!("  {PARASCOPY_NOTE}"));
            w.push(String::new());
        }
    }
    if let Some(why) = unchecked(s) {
        w.push("## Paralog Copy Number".to_owned());
        w.push("---".to_owned());
        w.push(format!("  {why}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The Status badge of a section that ran: Complete, or Stale for a result
/// from an earlier run.
pub fn done_badge(sec: &Section<'_>) -> Inline {
    if sec.state == State::Stale {
        Inline::Badge("Stale".to_owned(), Colour::Yellow)
    } else {
        Inline::Badge("Complete".to_owned(), Colour::Green)
    }
}

/// Whether a card shows its own lines: the step ran, or its state is one no
/// collector writes.
pub fn shows_body(sec: &Section<'_>) -> bool {
    !matches!(
        sec.state,
        State::Missing | State::Failed | State::Unreadable
    )
}

/// Paraphase's card.
fn paraphase_card(d: &Section<'_>) -> Result<Card, Malformed> {
    let mut body = vec![Line::Stat("Status".to_owned(), done_badge(d))];
    let mut caveats = Vec::new();
    if shows_body(d) {
        let path = format!("{PATH}.genes[]");
        let mut rows = String::new();
        for g in genes(d)? {
            let f = |k: &str| need(g, k, &path).and_then(py);
            rows += &format!(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                e(&f("gene")?),
                e(&f("copy_number")?),
                e(&f("status")?),
                e(&f("detail")?)
            );
        }
        if !rows.is_empty() {
            body.push(Line::Markup(format!(
                "    <table><tr><th>Gene</th><th>Copy number</th><th>Call</th><th>Detail</th></tr>{rows}</table>"
            )));
        }
        if let Some(smn) = smn1(d)? {
            body.push(Line::Stat(
                "SMN1 by both tools".to_owned(),
                Inline::Text(smn),
            ));
        }
        if let Some(n) = d.get("other_regions")?.filter(|v| truthy(v)) {
            body.push(Line::Stat(
                "Other Paraphase regions".to_owned(),
                Inline::Text(format!(
                    "{} ({} with a call)",
                    py(n)?,
                    py(d.need("other_regions_called")?)?
                )),
            ));
        }
        body.push(Line::Stat(
            "Variants on phased haplotypes".to_owned(),
            Inline::Text(variants(d)?),
        ));
        caveats.push(PARALOG_NOTE);
    }
    html::card(
        d,
        Card {
            title: "Paralog Genes (Paraphase)".to_owned(),
            full: true,
            caveats,
            body,
            ..Card::default()
        },
    )
}

/// Parascopy's card.
fn parascopy_card(d: &Section<'_>) -> Result<Card, Malformed> {
    let mut body = vec![Line::Stat("Status".to_owned(), done_badge(d))];
    let mut caveats = Vec::new();
    if shows_body(d) {
        let (lines, loci_caveats) = parascopy::html_lines(d)?;
        body.extend(lines);
        caveats.extend(loci_caveats);
        caveats.push(PARASCOPY_NOTE);
    }
    html::card(
        d,
        Card {
            title: "SMN1/SMN2 Copy Number (Parascopy)".to_owned(),
            caveats,
            body,
            ..Card::default()
        },
    )
}

/// The HTML report's cards.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let mut cards = Vec::new();
    let pg = s.section("paralogs");
    if pg.state != State::Missing {
        cards.push(paraphase_card(&pg)?);
    }
    let pc = s.section("parascopy");
    if pc.state != State::Missing {
        cards.push(parascopy_card(&pc)?);
    }
    if let Some(why) = unchecked(s) {
        cards.push(Card {
            title: "Paralog Copy Number".to_owned(),
            body: vec![
                Line::Stat(
                    "Status".to_owned(),
                    Inline::Badge("Not checked".to_owned(), Colour::Gray),
                ),
                Line::Note(why),
            ],
            ..Card::default()
        });
    }
    Ok(HtmlPart {
        cards,
        files: Vec::new(),
    })
}
