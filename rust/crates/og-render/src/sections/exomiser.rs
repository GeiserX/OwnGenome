//! Step 48's symptom-driven ranking (Exomiser): the terms searched, then
//! each ranked gene with the variants that put it there.

use pycompat::json::Value;

use crate::model::{
    Card, Cell, HtmlPart, Inline, Line, Malformed, Section, State, Summary, Table, get, is_str,
    need, py, truthy,
};
use crate::notes::EXOMISER_NOTE;
use crate::{html, text};

const PATH: &str = "sections.exomiser.data";

fn list<'a>(v: Option<&'a Value>, path: &str) -> Result<&'a [Value], Malformed> {
    match v {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!("{path} is not a list"))),
    }
}

/// The terms searched, each with its label when it has one; `.` for none.
pub fn terms(d: &Section<'_>) -> Result<String, Malformed> {
    let p = format!("{PATH}.terms[]");
    let mut out = Vec::new();
    for t in list(d.get("terms")?, &format!("{PATH}.terms"))? {
        let id = py(need(t, "id", &p)?)?;
        out.push(match get(t, "label", &p)?.filter(|v| truthy(v)) {
            Some(l) => format!("{} ({id})", py(l)?),
            None => id,
        });
    }
    Ok(if out.is_empty() {
        ".".to_owned()
    } else {
        out.join(", ")
    })
}

/// The loci written in `text` (`chr1:100`, `chr1:5-90`): the span of
/// every position on the first one's chromosome, as the Viewer cell opens
/// it; `None` when the text names no locus.
pub fn locus_span(text: &str) -> Option<(String, u128, u128)> {
    let b = text.as_bytes();
    let word = |c: u8| c.is_ascii_alphanumeric() || c == b'_' || c == b'.';
    let digits = |i: usize| b[i..].iter().take_while(|c| c.is_ascii_digit()).count();
    let mut found: Vec<(&str, u128, Option<u128>)> = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let w = b[i..].iter().take_while(|c| word(**c)).count();
        let colon = i + w;
        let n = if w > 0 && colon < b.len() && b[colon] == b':' {
            digits(colon + 1)
        } else {
            0
        };
        if n == 0 {
            i += 1;
            continue;
        }
        let a = colon + 1;
        let mut end = a + n;
        let mut second = None;
        if end < b.len() && b[end] == b'-' {
            let m = digits(end + 1);
            if m > 0 {
                second = text[end + 1..end + 1 + m].parse().ok();
                end += 1 + m;
            }
        }
        found.push((&text[i..colon], text[a..a + n].parse().ok()?, second));
        i = end;
    }
    let chrom = found.first()?.0;
    let pos: Vec<u128> = found
        .iter()
        .filter(|f| f.0 == chrom)
        .flat_map(|f| std::iter::once(f.1).chain(f.2))
        .collect();
    Some((chrom.to_owned(), *pos.iter().min()?, *pos.iter().max()?))
}

/// The Viewer cell of a row whose locus is written in its variant text.
fn viewer_of(text: &str) -> Cell {
    match locus_span(text) {
        Some((chrom, start, end)) if start != 0 && !chrom.is_empty() => {
            Cell::Viewer(Some(if end != start {
                format!("{chrom}:{start}-{end}")
            } else {
                format!("{chrom}:{start}")
            }))
        }
        _ => Cell::Viewer(None),
    }
}

fn genes<'a>(d: &Section<'a>) -> Result<&'a [Value], Malformed> {
    list(Some(d.need("genes")?), &format!("{PATH}.genes"))
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "exomiser",
        Some("Symptom-Driven Ranking (Exomiser, step 48)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = s.section("exomiser");
        let sample = py(s.need("sample")?)?;
        w.push(format!("  Symptoms searched: {}", terms(&d)?));
        let genes = genes(&d)?;
        if genes.is_empty() {
            w.push(
                "  No gene passed Exomiser's filters with a variant that fits a mode of inheritance."
                    .to_owned(),
            );
        }
        let gp = format!("{PATH}.genes[]");
        let vp = format!("{gp}.variants[]");
        for g in genes {
            let f = |k: &str| need(g, k, &gp).and_then(py);
            w.push(format!(
                "  {}. {} ({}) score {} (symptom match {})",
                f("rank")?,
                f("gene")?,
                f("moi")?,
                f("exomiser_score")?,
                f("phenotype_score")?
            ));
            for v in list(Some(need(g, "variants", &gp)?), &vp)? {
                let f = |k: &str| need(v, k, &vp).and_then(py);
                let disease = need(v, "disease", &vp)?;
                let disease = if is_str(Some(disease), ".") {
                    String::new()
                } else {
                    format!("; {}", py(disease)?)
                };
                w.push(format!(
                    "       {} {} {} {} ACMG {}{disease}",
                    f("variant")?,
                    f("genotype")?,
                    f("consequence")?,
                    f("hgvs")?,
                    f("acmg")?
                ));
            }
        }
        w.push(format!(
            "  {EXOMISER_NOTE} In words: exomiser/{sample}_exomiser.txt"
        ));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or was meant to, with a link
/// to the ranking table beside the report.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("exomiser");
    if d.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut part = HtmlPart::default();
    let mut card = Card {
        title: "Symptom-Driven Ranking (Exomiser)".to_owned(),
        full: true,
        ..Card::default()
    };
    if d.state.ran() {
        let sample = py(s.need("sample")?)?;
        let gp = format!("{PATH}.genes[]");
        let vp = format!("{gp}.variants[]");
        let mut rows = Vec::new();
        for g in genes(&d)? {
            for (i, v) in list(Some(need(g, "variants", &gp)?), &vp)?
                .iter()
                .enumerate()
            {
                let f = |k: &str| need(v, k, &vp).and_then(py);
                let mut row = if i == 0 {
                    ["rank", "gene", "moi", "exomiser_score", "phenotype_score"]
                        .iter()
                        .map(|k| need(g, k, &gp).and_then(py).map(Cell::Text))
                        .collect::<Result<Vec<_>, _>>()?
                } else {
                    vec![Cell::Text(String::new()); 5]
                };
                let variant = f("variant")?;
                row.push(Cell::Sub(
                    format!("{variant} {}", f("genotype")?),
                    format!("{} {}", f("consequence")?, f("hgvs")?),
                ));
                row.push(Cell::Text(f("acmg")?));
                row.push(Cell::Text(f("disease")?));
                row.push(viewer_of(&variant));
                rows.push(row);
            }
        }
        card.body.push(Line::Stat(
            "Symptoms searched".to_owned(),
            Inline::Text(terms(&d)?),
        ));
        if rows.is_empty() {
            card.body.push(Line::Stat(
                "Ranked genes".to_owned(),
                Inline::Fixed("none passed Exomiser's filters"),
            ));
        } else {
            card.body.push(Line::Table(Table {
                head: [
                    "Rank",
                    "Gene",
                    "Inheritance",
                    "Score",
                    "Symptom match",
                    "Variant",
                    "ACMG (Exomiser)",
                    "Linked disease",
                ]
                .map(str::to_owned)
                .to_vec(),
                rows,
                viewer: true,
                one_line: false,
            }));
        }
        card.caveats.push(EXOMISER_NOTE);
        card.body.push(Line::Note(format!(
            "Exomiser's own report: exomiser/{sample}_exomiser.html"
        )));
        if let Some(src) = d.source.filter(|v| truthy(v)) {
            part.files
                .push((py(src)?, "Exomiser ranking table".to_owned()));
        }
    }
    part.cards.push(html::card(&d, card)?);
    Ok(part)
}

#[cfg(test)]
mod tests {
    use super::locus_span;

    #[test]
    fn the_locus_span_reads_every_position_on_the_first_chromosome() {
        assert_eq!(
            locus_span("chr10:123256215 T>G"),
            Some(("chr10".to_owned(), 123_256_215, 123_256_215))
        );
        assert_eq!(
            locus_span("chr1:50 A>G + chr1:90 C>T + chr2:5"),
            Some(("chr1".to_owned(), 50, 90))
        );
        assert_eq!(
            locus_span("DEL chr1:5-90"),
            Some(("chr1".to_owned(), 5, 90))
        );
        assert_eq!(locus_span("chr1-1-A-G"), None);
        assert_eq!(locus_span("no locus: here"), None);
    }
}
