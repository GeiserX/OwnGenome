//! Step 42's traits panel: one row per catalog site with its genotype and
//! the cited reading, the genotype rules that read several sites together
//! (APOL1 with its G2, SERPINA1, HFE and the others), and the APOE lines
//! kept apart behind a warning.

use pycompat::json::Value;

use crate::escape::e;
use crate::html;
use crate::model::{
    Card, Colour, HtmlPart, Inline, Line, Malformed, Section, Summary, is_str, need, py, truthy,
};
use crate::notes::TRAITS_NOTE;
use crate::sections::{apoe, g2};
use crate::text;

const ROWS: &str = "sections.traits.data.rows[]";
const RULES: &str = "sections.traits.data.rules[]";

/// The small grey line under a cell's main text.
const SMALL: &str = "<br><span style=\"font-size:12px;color:#666\">";

/// The header row of the traits table and of the folded APOE table.
const HEAD: &str =
    "      <tr><th>Trait</th><th>Gene</th><th>Genotype</th><th>Reading</th><th>Source</th></tr>";

/// A list in the section's data, `[]` when it is absent or empty.
fn list<'a>(d: &Section<'a>, key: &str, need_it: bool) -> Result<&'a [Value], Malformed> {
    let v = if need_it {
        Some(d.need(key)?)
    } else {
        d.get(key)?.filter(|v| truthy(v))
    };
    match v {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!(
            "sections.traits.data.{key} is not a list"
        ))),
    }
}

/// ` [ancestry]` for a row whose reading depends on ancestry.
fn ancestry(r: &Value, path: &str) -> Result<bool, Malformed> {
    Ok(is_str(Some(need(r, "ancestry", path)?), "yes"))
}

/// The text report's lines (before the closing blank line).
pub fn text_lines(d: &Section<'_>, sample: &str) -> Result<Vec<String>, Malformed> {
    let n = |k: &str| d.need(k).and_then(py);
    let mut out = vec![format!(
        "  Read: {}  No genotype at the position: {}  Other allele: {}",
        n("called")?,
        n("not_covered")?,
        n("other_allele")?
    )];
    let rows = list(d, "rows", true)?;
    let mut missing = Vec::new();
    for r in rows {
        let f = |k: &str| need(r, k, ROWS).and_then(py);
        if is_str(Some(need(r, "status", ROWS)?), "not_covered") {
            missing.push(f("rsid")?);
            continue;
        }
        let mark = if ancestry(r, ROWS)? {
            " [ancestry]"
        } else {
            ""
        };
        out.push(format!(
            "  {} ({} {}): {}{mark}. {} (PMID {})",
            f("trait")?,
            f("gene")?,
            f("rsid")?,
            f("genotype")?,
            f("reading")?,
            f("pmids")?
        ));
    }
    if !missing.is_empty() {
        out.push(format!(
            "  Not read (no genotype in the gVCF): {}",
            missing.join(", ")
        ));
    }
    out.push(format!(
        "  {TRAITS_NOTE} Notes per row: traits/{sample}_traits.tsv"
    ));
    let rules = list(d, "rules", false)?;
    if !rules.is_empty() {
        out.push(format!(
            "  Genotype rules (several sites read together; traits/{sample}_genotype_rules.tsv):"
        ));
        for r in rules {
            let f = |k: &str| need(r, k, RULES).and_then(py);
            out.push(format!(
                "    {} ({} {}): {}. {} (PMID {})",
                f("trait")?,
                f("gene")?,
                f("rsids")?.replace(',', ", "),
                g2::rule_result(r)?,
                f("reading")?,
                f("pmids")?
            ));
            if let Some(line) = g2::g2_text(r)? {
                out.push(line);
            }
        }
    }
    out.extend(apoe::text_lines(list(d, "apoe", false)?)?);
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "traits", Some("Traits Panel (step 42)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let sample = py(s.need("sample")?)?;
        w.extend(text_lines(&s.section("traits"), &sample)?);
        w.push(String::new());
    }
    Ok(w)
}

/// The ancestry badge after a genotype or a result, when the row has one.
fn ancestry_badge(r: &Value, path: &str) -> Result<String, Malformed> {
    Ok(if ancestry(r, path)? {
        format!(" {}", html::badge("ancestry", Colour::Yellow))
    } else {
        String::new()
    })
}

/// One row of the traits table (and of the folded APOE table).
fn row(r: &Value) -> Result<String, Malformed> {
    let f = |k: &str| need(r, k, ROWS).and_then(py);
    Ok(format!(
        "<tr><td>{}</td><td>{}{SMALL}{}</span></td><td>{}{}</td><td>{}{SMALL}{}</span></td>\
         <td>PMID {}</td></tr>",
        e(&f("trait")?),
        e(&f("gene")?),
        e(&f("rsid")?),
        e(&f("genotype")?),
        ancestry_badge(r, ROWS)?,
        e(&f("reading")?),
        e(&f("note")?),
        e(&f("pmids")?)
    ))
}

/// The genotype rules table, none without rules.
fn rules_html(d: &Section<'_>, sample: &str) -> Result<Vec<Line>, Malformed> {
    let rules = list(d, "rules", false)?;
    if rules.is_empty() {
        return Ok(Vec::new());
    }
    let mut out = vec![
        Line::Markup("    <h3>Genotype rules (several sites read together)</h3>".to_owned()),
        Line::Markup("    <table>".to_owned()),
        Line::Markup(
            "      <tr><th>Trait</th><th>Gene and sites</th><th>Genotypes</th><th>Result</th>\
             <th>Reading</th><th>Source</th></tr>"
                .to_owned(),
        ),
    ];
    for r in rules {
        let f = |k: &str| need(r, k, RULES).and_then(py);
        out.push(Line::Markup(format!(
            "      <tr><td>{}</td><td>{}{SMALL}{}</span></td><td>{}</td><td>{}{}{}</td>\
             <td>{}{}{SMALL}{}</span></td><td>PMID {}</td></tr>",
            e(&f("trait")?),
            e(&f("gene")?),
            e(&f("rsids")?.replace(',', ", ")),
            e(&f("genotypes")?.replace(',', ", ")),
            e(&g2::rule_result(r)?),
            ancestry_badge(r, RULES)?,
            g2::g2_result_html(r)?,
            e(&f("reading")?),
            g2::g2_reading_html(r)?,
            e(&f("note")?),
            e(&f("pmids")?)
        )));
    }
    out.push(Line::Markup("    </table>".to_owned()));
    out.push(Line::Note(format!(
        "Every rule line: traits/{sample}_genotype_rules.tsv"
    )));
    Ok(out)
}

/// The traits card's lines: the counts, one table row per site read or not
/// read, the rules, and the APOE lines folded away behind a warning.
pub fn html_lines(d: &Section<'_>, sample: &str) -> Result<Vec<Line>, Malformed> {
    let n = |k: &str| d.need(k).and_then(py);
    let mut out = vec![
        Line::Stat(
            "Read / no genotype / other allele".to_owned(),
            Inline::Text(format!(
                "{} / {} / {}",
                n("called")?,
                n("not_covered")?,
                n("other_allele")?
            )),
        ),
        Line::Markup("    <table>".to_owned()),
        Line::Markup(HEAD.to_owned()),
    ];
    for r in list(d, "rows", true)? {
        out.push(Line::Markup(format!("      {}", row(r)?)));
    }
    out.push(Line::Markup("    </table>".to_owned()));
    out.extend(rules_html(d, sample)?);
    out.extend(
        apoe::html_lines(list(d, "apoe", false)?, HEAD, row)?
            .into_iter()
            .map(Line::Markup),
    );
    out.push(Line::Note(format!(
        "Every row with its note: traits/{sample}_traits.tsv"
    )));
    Ok(out)
}

/// The HTML report's card, shown on every report: "Not run" when step 42
/// is not in the summary.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("traits");
    let (body, caveats) = if sec.state.ran() {
        let sample = py(s.need("sample")?)?;
        (html_lines(&sec, &sample)?, vec![TRAITS_NOTE])
    } else {
        (Vec::new(), Vec::new())
    };
    let card = html::card(
        &sec,
        Card {
            title: "Traits".to_owned(),
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
