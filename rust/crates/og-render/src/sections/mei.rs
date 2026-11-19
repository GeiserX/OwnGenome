//! Step 50's mobile element insertions (MEGAnE): the calls per element,
//! then each known disease insertion of the catalog, seen or not.

use pycompat::json::Value;

use crate::model::{
    Card, Cell, HtmlPart, Inline, Layout, Line, Malformed, Section, State, Summary, Table, is_str,
    need, py, truthy,
};
use crate::notes::MEI_NOTE;
use crate::{html, text};

const PATH: &str = "sections.mei.data.sites[]";

fn sites<'a>(d: &Section<'a>) -> Result<&'a [Value], Malformed> {
    match d.get("sites")? {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed("sections.mei.data.sites is not a list".into())),
    }
}

/// `"2 (Alu 1, SVA 1)"`: the calls and their count per element.
fn counts(d: &Section<'_>) -> Result<String, Malformed> {
    let by = match d.get("by_element")?.filter(|v| truthy(v)) {
        None => String::new(),
        Some(Value::Object(o)) => o
            .iter()
            .map(|(k, v)| py(v).map(|v| format!("{k} {v}")))
            .collect::<Result<Vec<_>, _>>()?
            .join(", "),
        Some(_) => {
            return Err(Malformed(
                "sections.mei.data.by_element is not an object".into(),
            ));
        }
    };
    let by = if by.is_empty() { "none".to_owned() } else { by };
    Ok(format!(
        "{} ({by})",
        d.get("calls")?.map_or(Ok("0".to_owned()), py)?
    ))
}

/// Step 50's lines, the same in both reports: the calls per element, then
/// each catalog row, then the caveat.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let mut out = vec![format!("Insertions called: {}", counts(d)?)];
    for x in sites(d)? {
        let f = |k: &str| need(x, k, PATH).and_then(py);
        let call = need(x, "call", PATH)?;
        let call = if is_str(Some(call), ".") {
            String::new()
        } else {
            format!(", {}", py(call)?)
        };
        out.push(format!(
            "{} {} ({}; {}): {} (depth {}{call})",
            f("gene")?,
            f("element")?,
            f("condition")?,
            f("inheritance")?,
            f("status")?,
            f("mean_depth")?
        ));
    }
    out.push(MEI_NOTE.to_owned());
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "mei", Some("Mobile Element Insertions (step 50)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        w.extend(
            lines(&s.section("mei"))?
                .into_iter()
                .map(|l| format!("  {l}")),
        );
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or was meant to.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("mei");
    if d.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut card = Card {
        title: "Mobile Element Insertions".to_owned(),
        full: true,
        caveats: vec![MEI_NOTE],
        ..Card::default()
    };
    if d.state.ran() || d.state == State::Other {
        let mut rows = Vec::new();
        for x in sites(&d)? {
            let mut row = Vec::new();
            for k in [
                "gene",
                "element",
                "condition",
                "inheritance",
                "status",
                "mean_depth",
                "call",
            ] {
                row.push(Cell::Text(py(need(x, k, PATH)?)?));
            }
            rows.push(row);
        }
        card.body = vec![
            Line::Stat("Insertions called".to_owned(), Inline::Text(counts(&d)?)),
            Line::Table(Table {
                head: [
                    "Gene",
                    "Element",
                    "Condition",
                    "Inheritance",
                    "Status",
                    "Depth",
                    "Call",
                ]
                .map(str::to_owned)
                .to_vec(),
                rows,
                viewer: false,
                layout: Layout::OneLine,
            }),
        ];
    }
    Ok(HtmlPart {
        cards: vec![html::card(&d, card)?],
        files: Vec::new(),
    })
}
