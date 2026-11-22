//! Step 42's traits panel: one row per catalog site with its genotype and
//! the cited reading, the genotype rules that read several sites together
//! (APOL1 with its G2, SERPINA1, HFE and the others), and the APOE lines
//! kept apart behind a warning.

use pycompat::json::Value;

use crate::html;
use crate::model::{
    Card, Cell, Colour, HtmlPart, Inline, Line, Malformed, Section, Summary, Table, get, is_str,
    need, py, truthy,
};
use crate::notes::TRAITS_NOTE;
use crate::sections::{apoe, g2};
use crate::text;

const ROWS: &str = "sections.traits.data.rows[]";
const RULES: &str = "sections.traits.data.rules[]";
pub const LOW_GQ_LABEL: &str = "Low-quality genotype calls";
pub const LOW_GQ_READING: &str = "not read (low genotype quality)";

/// A count to display when weak calls were included in the summary.
pub fn low_gq_count(data: &Value) -> Option<i128> {
    let Value::Object(fields) = data else {
        return None;
    };
    match fields.get("low_gq") {
        Some(Value::Int(count)) if *count > 0 => Some(*count),
        _ => None,
    }
}

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
    if let Some(count) = low_gq_count(d.data) {
        out.push(format!("  {LOW_GQ_LABEL}: {count}"));
    }
    let rows = list(d, "rows", true)?;
    let mut missing = Vec::new();
    for r in rows {
        let f = |k: &str| need(r, k, ROWS).and_then(py);
        let status = need(r, "status", ROWS)?;
        if is_str(Some(status), "not_covered") {
            missing.push(f("rsid")?);
            continue;
        }
        if is_str(Some(status), "low_gq") {
            out.push(format!(
                "  {} ({} {}): {LOW_GQ_READING}.",
                f("trait")?,
                f("gene")?,
                f("rsid")?
            ));
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
fn ancestry_badge(r: &Value, path: &str) -> Result<Vec<Inline>, Malformed> {
    Ok(if ancestry(r, path)? {
        vec![
            Inline::Text(" ".to_owned()),
            Inline::Badge("ancestry".to_owned(), Colour::Yellow),
        ]
    } else {
        Vec::new()
    })
}

/// The column titles of the traits table and of the folded APOE table.
fn head() -> Vec<String> {
    ["Trait", "Gene", "Genotype", "Reading", "Source"]
        .map(str::to_owned)
        .to_vec()
}

/// One row of the traits table (and of the folded APOE table).
fn row(r: &Value) -> Result<Vec<Cell>, Malformed> {
    let f = |k: &str| need(r, k, ROWS).and_then(py);
    let trait_ = f("trait")?;
    let gene = Cell::Sub(f("gene")?, f("rsid")?);
    if is_str(get(r, "status", ROWS)?, "low_gq") {
        return Ok(vec![
            Cell::Text(trait_),
            gene,
            Cell::Text(".".to_owned()),
            Cell::Sub(LOW_GQ_READING.to_owned(), f("note")?),
            Cell::Text(format!("PMID {}", f("pmids")?)),
        ]);
    }
    let mut genotype = vec![Inline::Text(f("genotype")?)];
    genotype.extend(ancestry_badge(r, ROWS)?);
    Ok(vec![
        Cell::Text(trait_),
        gene,
        Cell::Parts(genotype),
        Cell::Sub(f("reading")?, f("note")?),
        Cell::Text(format!("PMID {}", f("pmids")?)),
    ])
}

/// The genotype rules table, none without rules.
fn rules_html(d: &Section<'_>, sample: &str) -> Result<Vec<Line>, Malformed> {
    let rules = list(d, "rules", false)?;
    if rules.is_empty() {
        return Ok(Vec::new());
    }
    let mut rows = Vec::new();
    for r in rules {
        let f = |k: &str| need(r, k, RULES).and_then(py);
        let trait_ = f("trait")?;
        let gene = Cell::Sub(f("gene")?, f("rsids")?.replace(',', ", "));
        let genotypes = f("genotypes")?.replace(',', ", ");
        let mut result = vec![Inline::Text(g2::rule_result(r)?)];
        result.extend(ancestry_badge(r, RULES)?);
        result.extend(g2::g2_result_html(r)?);
        let mut reading = vec![Inline::Text(f("reading")?)];
        reading.extend(g2::g2_reading_html(r)?);
        reading.push(Inline::Small(f("note")?));
        rows.push(vec![
            Cell::Text(trait_),
            gene,
            Cell::Text(genotypes),
            Cell::Parts(result),
            Cell::Parts(reading),
            Cell::Text(format!("PMID {}", f("pmids")?)),
        ]);
    }
    Ok(vec![
        Line::Raw("    <h3>Genotype rules (several sites read together)</h3>".to_owned()),
        Line::Table(Table {
            head: [
                "Trait",
                "Gene and sites",
                "Genotypes",
                "Result",
                "Reading",
                "Source",
            ]
            .map(str::to_owned)
            .to_vec(),
            rows,
            ..Table::default()
        }),
        Line::Note(format!(
            "Every rule line: traits/{sample}_genotype_rules.tsv"
        )),
    ])
}

/// The traits card's lines: the counts, one table row per site read or not
/// read, the rules, and the APOE lines folded away behind a warning.
pub fn html_lines(d: &Section<'_>, sample: &str) -> Result<Vec<Line>, Malformed> {
    let n = |k: &str| d.need(k).and_then(py);
    let mut out = vec![Line::Stat(
        "Read / no genotype / other allele".to_owned(),
        Inline::Text(format!(
            "{} / {} / {}",
            n("called")?,
            n("not_covered")?,
            n("other_allele")?
        )),
    )];
    if let Some(count) = low_gq_count(d.data) {
        out.push(Line::Stat(
            LOW_GQ_LABEL.to_owned(),
            Inline::Text(count.to_string()),
        ));
    }
    out.push(Line::Table(Table {
        head: head(),
        rows: list(d, "rows", true)?
            .iter()
            .map(row)
            .collect::<Result<_, _>>()?,
        ..Table::default()
    }));
    out.extend(rules_html(d, sample)?);
    out.extend(apoe::html_line(list(d, "apoe", false)?, head(), row)?);
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
