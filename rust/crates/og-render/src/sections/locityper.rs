//! Step 53's Locityper second opinion at hard loci: per locus Locityper's
//! two haplotypes with the authors' filter, beside the primary caller's
//! call.

use pycompat::json::Value;

use crate::model::{
    Card, Cell, Colour, HtmlPart, Inline, Layout, Line, Malformed, Section, State, Summary, Table,
    get, is_str, need, py,
};
use crate::notes::LOCITYPER_NOTE;
use crate::{html, text};

const PATH: &str = "sections.locityper.data.loci[]";

fn loci<'a>(d: &Section<'a>) -> Result<&'a [Value], Malformed> {
    match d.get("loci")? {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(
            "sections.locityper.data.loci is not a list".into(),
        )),
    }
}

/// One locus: the two haplotypes with quality and filter, or why there is
/// no call.
pub fn locus_text(x: &Value) -> Result<String, Malformed> {
    let f = |k: &str| need(x, k, PATH).and_then(py);
    if !is_str(Some(need(x, "status", PATH)?), "called") {
        return Ok(format!("not called: {}", f("reason")?));
    }
    let why = if is_str(Some(need(x, "filter", PATH)?), "fail")
        && !is_str(Some(need(x, "reason", PATH)?), ".")
    {
        format!(": {}", f("reason")?)
    } else {
        String::new()
    };
    Ok(format!(
        "{} + {} (quality {}, {} of {} reads unexplained, weighted distance {}; filter {}{why})",
        f("haplotype_1")?,
        f("haplotype_2")?,
        f("quality")?,
        f("unexplained_reads")?,
        f("total_reads")?,
        f("weighted_distance")?,
        f("filter")?
    ))
}

/// `x.get("primary", ".")`.
fn primary(x: &Value) -> Result<String, Malformed> {
    get(x, "primary", PATH)?.map_or(Ok(".".to_owned()), py)
}

/// Step 53's lines: the counts, then per locus Locityper's call and the
/// primary caller's.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let loci = loci(d)?;
    let mut out = vec![format!(
        "Loci called: {} of {}; passing the authors' filter: {}",
        d.get("called")?.map_or(Ok("0".to_owned()), py)?,
        loci.len(),
        d.get("passed")?.map_or(Ok("0".to_owned()), py)?
    )];
    for x in loci {
        let locus = py(need(x, "locus", PATH)?)?;
        out.push(format!(
            "{} Locityper: {}",
            pycompat::fmt::ljust(&locus, 9),
            locus_text(x)?
        ));
        out.push(format!("{} primary:   {}", " ".repeat(9), primary(x)?));
    }
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "locityper",
        Some("Locityper Second Opinion at Hard Loci (step 53)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        w.extend(
            lines(&s.section("locityper"))?
                .into_iter()
                .map(|l| format!("  {l}")),
        );
        w.push(format!("  {LOCITYPER_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or was meant to.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("locityper");
    if d.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut card = Card {
        title: "Locityper Second Opinion at Hard Loci".to_owned(),
        full: true,
        ..Card::default()
    };
    if d.state.ran() {
        let loci = loci(&d)?;
        let mut rows = Vec::new();
        for x in loci {
            let filter = need(x, "filter", PATH)?;
            let colour = if is_str(Some(filter), "pass") {
                Colour::Green
            } else if is_str(Some(filter), "fail") {
                Colour::Yellow
            } else {
                Colour::Gray
            };
            rows.push(vec![
                Cell::Text(py(need(x, "locus", PATH)?)?),
                Cell::Text(locus_text(x)?),
                Cell::Badge(py(filter)?, colour),
                Cell::Text(primary(x)?),
            ]);
        }
        card.caveats.push(LOCITYPER_NOTE);
        card.body = vec![
            Line::Stat(
                "Loci called".to_owned(),
                Inline::Text(format!(
                    "{} of {} ({} pass the filter)",
                    py(d.need("called")?)?,
                    loci.len(),
                    py(d.need("passed")?)?
                )),
            ),
            Line::Table(Table {
                head: ["Locus", "Locityper", "Filter", "Primary caller"]
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
