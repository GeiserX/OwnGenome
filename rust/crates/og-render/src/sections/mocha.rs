//! Step 51's mosaic chromosomal alterations (MoChA): the cell-fraction
//! floor at this depth, then each reported call, or why there is none.

use pycompat::json::Value;

use crate::model::{
    Card, Cell, Colour, HtmlPart, Inline, Layout, Line, Malformed, Section, State, Summary, Table,
    is_str, need, py,
};
use crate::notes::MCA_NOTE;
use crate::{html, text};

const PATH: &str = "sections.mosaic.data.calls[]";

fn calls<'a>(d: &Section<'a>) -> Result<&'a [Value], Malformed> {
    match d.get("calls")? {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed("sections.mosaic.data.calls is not a list".into())),
    }
}

/// `d.get(key, default)` as text.
fn or(d: &Section<'_>, key: &str, default: &str) -> Result<String, Malformed> {
    d.get(key)?.map_or(Ok(default.to_owned()), py)
}

fn assessed(d: &Section<'_>) -> Result<bool, Malformed> {
    Ok(is_str(d.get("status")?, "assessed"))
}

/// Step 51's lines, the same in both reports: the floor, then each call or
/// the reason none is shown, then the caveat.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    if !assessed(d)? {
        return Ok(vec![
            format!("Not assessed: {}", or(d, "message", ".")?),
            MCA_NOTE.to_owned(),
        ]);
    }
    let mut out = vec![format!(
        "Reported at a cell fraction of {} or more at {}x (read-sampling floor {}, minimum {}; {} heterozygous \
         sites on {})",
        or(d, "floor", ".")?,
        or(d, "depth", ".")?,
        or(d, "floor_read", ".")?,
        or(d, "floor_min", ".")?,
        or(d, "het_sites", ".")?,
        or(d, "chromosomes", ".")?
    )];
    let calls = calls(d)?;
    for x in calls {
        let f = |k: &str| need(x, k, PATH).and_then(py);
        out.push(format!(
            "{} on {}, {}:{}-{} ({} Mb): about {} of the cells",
            f("type")?,
            f("arms")?,
            f("chrom")?,
            f("start")?,
            f("end")?,
            f("length_mb")?,
            f("cell_fraction")?
        ));
    }
    if calls.is_empty() {
        out.push(format!(
            "No mosaic alteration reported ({} MoChA calls, none kept).",
            or(d, "calls_mocha", "0")?
        ));
    }
    out.push(MCA_NOTE.to_owned());
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "mosaic",
        Some("Mosaic Chromosomal Alterations (MoChA, step 51)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        w.extend(
            lines(&s.section("mosaic"))?
                .into_iter()
                .map(|l| format!("  {l}")),
        );
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or was meant to.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("mosaic");
    if d.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut card = Card {
        title: "Mosaic Chromosomal Alterations (MoChA)".to_owned(),
        full: true,
        ..Card::default()
    };
    if d.state.ran() {
        card.caveats.push(MCA_NOTE);
        if !assessed(&d)? {
            card.body = vec![
                Line::Stat(
                    "Status".to_owned(),
                    Inline::Badge("Not assessed".to_owned(), Colour::Gray),
                ),
                Line::Note(or(&d, "message", ".")?),
            ];
        } else {
            let mut rows = Vec::new();
            for x in calls(&d)? {
                let f = |k: &str| need(x, k, PATH).and_then(py);
                rows.push(vec![
                    Cell::Text(f("type")?),
                    Cell::Text(f("arms")?),
                    Cell::Text(format!("{}:{}-{}", f("chrom")?, f("start")?, f("end")?)),
                    Cell::Text(f("length_mb")?),
                    Cell::Text(f("cell_fraction")?),
                    Cell::viewer(
                        need(x, "chrom", PATH)?,
                        need(x, "start", PATH)?,
                        Some(need(x, "end", PATH)?),
                    )?,
                ]);
            }
            card.body = vec![
                Line::Stat(
                    "Calls reported".to_owned(),
                    Inline::Text(or(&d, "calls_reported", "0")?),
                ),
                Line::Stat(
                    "Cell-fraction floor".to_owned(),
                    Inline::Text(or(&d, "floor", ".")?),
                ),
            ];
            if rows.is_empty() {
                card.body
                    .push(Line::Note("No mosaic alteration reported.".to_owned()));
            } else {
                card.body.push(Line::Table(Table {
                    head: ["Type", "Arm", "Region", "Mb", "Cell fraction"]
                        .map(str::to_owned)
                        .to_vec(),
                    rows,
                    viewer: true,
                    layout: Layout::OneLine,
                }));
                card.body.push(Line::ViewerNote);
            }
            card.body.push(Line::Note(lines(&d)?.swap_remove(0)));
        }
    }
    Ok(HtmlPart {
        cards: vec![html::card(&d, card)?],
        files: Vec::new(),
    })
}
