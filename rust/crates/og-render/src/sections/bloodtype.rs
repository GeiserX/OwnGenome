//! Step 41's blood groups, predicted from the genome: one row per system
//! (ABO, RhD and the others the marker table defines), with its genotype and
//! confidence. The markers' sources and licences are in the header of
//! `assets/bloodtype_markers.tsv`.

use pycompat::fmt::ljust;
use pycompat::json::Value;

use crate::html;
use crate::model::{
    Card, Cell, HtmlPart, Inline, Layout, Line, Malformed, Summary, Table, need, py,
};
use crate::notes::BLOODTYPE_NOTE;
use crate::sections::paralogs::shows_body;
use crate::text;

const PATH: &str = "sections.bloodtype.data.systems[]";

/// The systems, or an error when the summary holds something else.
fn systems(v: &Value) -> Result<&[Value], Malformed> {
    match v {
        Value::Array(a) => Ok(a),
        _ => Err(Malformed(
            "sections.bloodtype.data.systems is not a list".into(),
        )),
    }
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "bloodtype",
        Some("Blood Groups (predicted from the genome)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = s.section("bloodtype");
        for x in systems(d.need("systems")?)? {
            let f = |k: &str| need(x, k, PATH).and_then(py);
            w.push(format!(
                "  {} {} {} ({})",
                ljust(&f("system")?, 9),
                ljust(&f("phenotype")?, 34),
                ljust(&f("genotype")?, 28),
                f("confidence")?
            ));
        }
        w.push(format!("  {BLOODTYPE_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when step 41 is in the summary.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("bloodtype");
    if d.state == crate::model::State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut body = Vec::new();
    let mut caveats = Vec::new();
    if shows_body(&d) {
        let dot = |k: &str| d.get(k).and_then(|v| v.map_or(Ok(".".to_owned()), py));
        body.push(Line::Stat("ABO".to_owned(), Inline::Text(dot("abo")?)));
        body.push(Line::Stat("RhD".to_owned(), Inline::Text(dot("rhd")?)));
        let mut rows = Vec::new();
        for x in d.get("systems")?.map_or(Ok(&[][..]), systems)? {
            let f = |k: &str| need(x, k, PATH).and_then(py).map(Cell::Text);
            rows.push(vec![
                f("system")?,
                f("phenotype")?,
                f("genotype")?,
                f("confidence")?,
            ]);
        }
        body.push(Line::Table(Table {
            head: ["System", "Predicted phenotype", "Genotype", "Confidence"]
                .map(str::to_owned)
                .to_vec(),
            rows,
            layout: Layout::OneLine,
            ..Table::default()
        }));
        caveats.push(BLOODTYPE_NOTE);
    }
    let card = html::card(
        &d,
        Card {
            title: "Blood Groups (predicted)".to_owned(),
            full: true,
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
