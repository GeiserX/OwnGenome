//! Ancient populations (step 44): where the genome falls on a PCA of
//! present-day West Eurasians, and the ancient groups of the AADR nearest
//! to it there, with the plot.

use pycompat::json::Value;

use crate::escape::e;
use crate::model::{
    Card, Cell, HtmlPart, Layout, Line, Malformed, State, Summary, Table, get, py, truthy,
};
use crate::notes::ANCIENT_NOTE;
use crate::sections::ancestry::{Align, get_text, need_align, need_text};
use crate::{html, text};

const PATH: &str = "sections.ancient.data";
const ROW: &str = "sections.ancient.data.nearest[]";

/// `(d.get("nearest") or [])[:10]`: the ten nearest groups.
fn nearest(d: &Value) -> Result<Vec<&Value>, Malformed> {
    match get(d, "nearest", PATH)?.filter(|v| truthy(v)) {
        None => Ok(Vec::new()),
        Some(Value::Array(a)) => Ok(a.iter().take(10).collect()),
        Some(_) => Err(Malformed(format!("{PATH}.nearest is not a list"))),
    }
}

/// Step 44's lines, the same in both reports: what the PCA was built from,
/// then the nearest groups.
pub fn lines(d: &Value) -> Result<Vec<String>, Malformed> {
    let g = |k: &str| get_text(d, k, ".", PATH);
    let mut out = vec![
        format!(
            "AADR {}: {} of {} sites used; axes from {} present-day West Eurasians, {} ancient \
             people in {} groups projected.",
            g("release")?,
            g("sites_used")?,
            g("panel_sites")?,
            g("present_day_people")?,
            g("ancient_people")?,
            g("ancient_groups")?
        ),
        "Nearest ancient groups:".to_owned(),
    ];
    for r in nearest(d)? {
        let n = |k: &str| need_text(r, k, ROW);
        out.push(format!(
            "  {}. {} ({}, {}; {} people)",
            need_align(r, "rank", 2, Align::Right, ROW)?,
            n("group")?,
            n("dates")?,
            n("country")?,
            n("n")?
        ));
    }
    Ok(out)
}

/// The text report's block, with the plot's file name.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "ancient",
        Some("Ancient Populations (AADR PCA, step 44)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        w.extend(
            lines(s.section("ancient").data)?
                .iter()
                .map(|l| format!("  {l}")),
        );
        w.push(format!(
            "  {ANCIENT_NOTE} Plot: ancient/{}_ancient_pca.svg",
            py(s.need("sample")?)?
        ));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card: the nearest groups as a table, and the plot as
/// an image (shown, never run as part of the page).
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("ancient");
    if sec.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut c = Card {
        title: "Ancient Populations".to_owned(),
        full: true,
        ..Card::default()
    };
    if crate::sections::haplogroups::shows_body(&sec) {
        let d = sec.data;
        let first = lines(d)?.swap_remove(0);
        let mut rows = Vec::new();
        for r in nearest(d)? {
            rows.push(
                ["rank", "group", "dates", "country", "locality", "n"]
                    .iter()
                    .map(|k| need_text(r, k, ROW).map(Cell::Text))
                    .collect::<Result<Vec<_>, Malformed>>()?,
            );
        }
        c.body = vec![
            Line::Para(first),
            Line::Table(Table {
                head: ["#", "Ancient group", "Dates", "Country", "Place", "People"]
                    .map(str::to_owned)
                    .to_vec(),
                rows,
                layout: Layout::RowLines,
                ..Table::default()
            }),
        ];
        if let Some(svg) = get(d, "plot_svg", PATH)?.filter(|v| truthy(v)) {
            c.body.push(Line::Raw(format!(
                "    <img alt=\"PC1 against PC2: present-day West Eurasians in grey, ancient people \
                 coloured by age, the sample marked\" style=\"max-width:100%\" \
                 src=\"data:image/svg+xml;base64,{}\">",
                e(&py(svg)?)
            )));
        }
        c.caveats.push(ANCIENT_NOTE);
    }
    Ok(HtmlPart {
        cards: vec![html::card(&sec, c)?],
        files: Vec::new(),
    })
}
