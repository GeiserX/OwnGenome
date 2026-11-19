//! Step 17's cancer predisposition screen (CPSR): its classification
//! breakdown, the two ClinVar releases, CPSR's own benign calls among the
//! final P/LP rows, the rows where CPSR and ClinVar disagree, the rows below
//! the depth floor and the cDNA read check.

use crate::model::{Card, Colour, HtmlPart, Inline, Line, Malformed, Summary};
use crate::notes::CPSR_NOTE;
use crate::sections::cdna::cdna_text;
use crate::sections::clinvar::{D, call_text};
use crate::{html, text};

/// One P/LP or disagreement row: variant, gene, the classes and codes.
pub fn row_text(x: &D<'_>) -> Result<String, Malformed> {
    let mut first = format!("final {}", x.s("final")?);
    if let Some(a) = x.opt("authority")? {
        first += &format!(" ({})", a.text()?);
    }
    let mut parts = vec![first];
    if let Some(c) = x.opt("clinvar")? {
        parts.push(format!("ClinVar {}", c.text()?));
    }
    if let Some(c) = x.opt("cpsr")? {
        let mut p = format!("CPSR's own {}", c.text()?);
        if let Some(codes) = x.opt("acmg_codes")? {
            p += &format!(" ({})", codes.text()?);
        }
        parts.push(p);
    }
    Ok(format!(
        "{} {}: {}",
        x.s("gene")?,
        x.s("variant")?,
        parts.join(", ")
    ))
}

/// `dbs.get(key) or 'unknown'`.
fn release(dbs: &D<'_>, key: &str) -> Result<String, Malformed> {
    match dbs.opt(key)? {
        Some(v) => v.text(),
        None => Ok("unknown".to_owned()),
    }
}

/// The lines both reports share beyond the breakdown.
pub fn lines(d: &D<'_>, dbs: &D<'_>) -> Result<Vec<String>, Malformed> {
    let mut out = Vec::new();
    if !d.has("plp") {
        return Ok(out);
    }
    out.push(format!(
        "ClinVar release in CPSR's bundle: {}; step 6 used the ClinVar file of {}",
        release(dbs, "cpsr_clinvar_release")?,
        release(dbs, "clinvar_release")?
    ));
    if d.flag("own_class_column")? {
        out.push(format!(
            "Final P/LP rows that CPSR's own rules class benign or likely benign: {}",
            d.s("final_plp_own_benign")?
        ));
        let dis = d.at("disagreements")?;
        out.push(format!(
            "Rows where ClinVar and CPSR's own class disagree on P/LP: {}",
            dis.len()?
        ));
        for x in dis.items()? {
            out.push(format!("  {}", row_text(&x)?));
        }
    } else {
        out.push(
            "The classification table has no CPSR_CLASSIFICATION column, so no disagreement can be shown."
                .to_owned(),
        );
    }
    out.push(format!(
        "P/LP rows below the depth floor (DP below {}): {}",
        d.s("call_floor_dp")?,
        d.s("low_confidence")?
    ));
    let plp = d.at("plp")?.items()?;
    for x in &plp {
        if x.flag("low_confidence")? {
            out.push(format!("  {}; {}", row_text(x)?, call_text(x)?));
        }
    }
    if d.has("cdna_flagged") {
        out.push(format!(
            "cDNA read check of the P/LP rows (alt reads clipped at exon edges): {} flagged",
            d.s("cdna_flagged")?
        ));
        for x in &plp {
            if x.at("cdna")?.at("check")?.is("flagged") {
                out.push(format!("  {}; reads: {}", row_text(x)?, cdna_text(x)?));
            }
        }
    }
    Ok(out)
}

fn databases(s: &Summary) -> Result<D<'_>, Malformed> {
    Ok(D::new(s.need("databases")?, "summary.databases"))
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "cpsr", Some("Cancer Predisposition Screening (CPSR)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = D::of(&s.section("cpsr"));
    if d.flag("html_report")? {
        w.push("  HTML report in cpsr/ (open in a browser)".to_owned());
    }
    if let Some(cls) = d.opt("classification")? {
        w.push(format!(
            "  Classification breakdown ({}):",
            d.s("classification_column")?
        ));
        for (k, v) in cls.pairs()? {
            w.push(format!("    {k}: {}", v.text()?));
        }
    } else if d.has("classification_column") {
        w.push("  Classification table has no CLASSIFICATION column.".to_owned());
    }
    for line in lines(&d, &databases(s)?)? {
        w.push(format!("  {line}"));
    }
    w.push(format!("  {CPSR_NOTE}"));
    w.push(String::new());
    Ok(w)
}

/// Python's `str.strip()`: whitespace as `str.isspace` counts it.
fn strip(s: &str) -> &str {
    pycompat::text::strip(s, None)
}

/// `os.path.dirname(p)`.
fn dirname(p: &str) -> &str {
    let Some(i) = p.rfind('/') else {
        return "";
    };
    let head = &p[..=i];
    if head.chars().all(|c| c == '/') {
        head
    } else {
        head.trim_end_matches('/')
    }
}

/// The HTML report's card, always there, with a link to CPSR's own report
/// beside this one.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("cpsr");
    let d = D::of(&sec);
    let report = if d.flag("html_report")? {
        Inline::Badge("Complete".to_owned(), Colour::State)
    } else {
        Inline::Badge("No HTML".to_owned(), Colour::Gray)
    };
    let mut body = vec![Line::Stat("CPSR report".to_owned(), report)];
    if let Some(cls) = d.opt("classification")? {
        for (k, v) in cls.pairs()? {
            body.push(Line::Stat(k.to_owned(), Inline::Text(v.text()?)));
        }
    }
    for line in lines(&d, &databases(s)?)? {
        body.push(Line::Note(strip(&line).to_owned()));
    }
    let card = html::card(
        &sec,
        Card {
            title: "Cancer Predisposition".to_owned(),
            caveats: vec![CPSR_NOTE],
            body,
            ..Card::default()
        },
    )?;
    let mut files = Vec::new();
    let src = sec.source.filter(|v| crate::model::truthy(v));
    if sec.state.ran()
        && let Some(src) = src
        && d.flag("html_report")?
    {
        let src = D::new(src, "sections.cpsr.source");
        let src = src.string()?;
        files.push(if src.ends_with(".html") {
            (src.to_owned(), "CPSR HTML report".to_owned())
        } else {
            (
                format!("{}/", dirname(src)),
                "the cpsr folder, which holds the CPSR HTML report".to_owned(),
            )
        });
    }
    Ok(HtmlPart {
        cards: vec![card],
        files,
    })
}
