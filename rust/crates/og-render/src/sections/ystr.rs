//! The Y-STR markers (step 45, HipSTR): the repeat count of each marker
//! called, in FamilyTreeDNA's Y-111 order and panels, so a count lines up
//! with FamilyTreeDNA results and surname-project tables.

use pycompat::json::Value;

use crate::escape::e;
use crate::model::{Card, Line, Malformed, State, Summary, need, truthy};
use crate::notes::YSTR_NOTE;
use crate::sections::ancestry::{Align, items, need_align, need_text, py_len};
use crate::{html, text};

const PATH: &str = "sections.y_str.data";
const ROW: &str = "sections.y_str.data.markers[]";

/// How many markers were called, and with what confidence.
pub fn counts(d: &Value) -> Result<String, Malformed> {
    let n = |k: &str| need_text(d, k, PATH);
    let high = n("high")?;
    let low = n("low")?;
    let all = py_len(need(d, "markers", PATH)?, &format!("{PATH}.markers"))?;
    Ok(format!(
        "Markers called: {high} high and {low} low confidence of {all}; {} no call, {} not \
         callable from short reads",
        n("no_call")?,
        n("not_callable")?
    ))
}

/// The markers with a count, in the order the summary lists them (the
/// panel's order).
fn called(markers: &Value) -> Result<Vec<Value>, Malformed> {
    let mut out = Vec::new();
    for x in items(markers, &format!("{PATH}.markers"))? {
        if truthy(need(&x, "count", ROW)?) {
            out.push(x);
        }
    }
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "y_str", Some("Y-STR Markers (HipSTR)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = s.section("y_str").data;
    w.push(format!("  {}", counts(d)?));
    for x in called(need(d, "markers", PATH)?)? {
        w.push(format!(
            "  {} {} {} ({})",
            need_align(&x, "panel", 6, Align::Left, ROW)?,
            need_align(&x, "marker", 13, Align::Left, ROW)?,
            need_align(&x, "count", 6, Align::Left, ROW)?,
            need_text(&x, "confidence", ROW)?
        ));
    }
    w.push(format!("  {YSTR_NOTE}"));
    w.push(String::new());
    Ok(w)
}

/// The HTML report's card, when step 45 has a result.
pub fn card(s: &Summary) -> Result<Option<Card>, Malformed> {
    let sec = s.section("y_str");
    if sec.state == State::Missing {
        return Ok(None);
    }
    let mut c = Card {
        title: "Y-STR Markers (HipSTR)".to_owned(),
        full: true,
        ..Card::default()
    };
    if sec.state.ran() {
        let d = sec.data;
        let markers = crate::model::get(d, "markers", PATH)?;
        let mut rows = String::new();
        if let Some(m) = markers {
            for x in called(m)? {
                rows += "<tr>";
                for k in ["panel", "marker", "count", "confidence"] {
                    rows += &format!("<td>{}</td>", e(&need_text(&x, k, ROW)?));
                }
                rows += "</tr>";
            }
        }
        c.body = vec![
            Line::Para(counts(d)?),
            Line::Raw(format!(
                "    <table><tr><th>Panel</th><th>Marker</th><th>Repeats</th><th>Confidence</th></tr>\
                 {rows}</table>"
            )),
        ];
        c.caveats.push(YSTR_NOTE);
    }
    Ok(Some(html::card(&sec, c)?))
}
