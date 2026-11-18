//! Step 05's AnnotSV class 4 and 5 SVs, with how many callers saw each
//! (step 22's consensus), the exon-edge and benign-set checks, and the
//! calls step 05's rules moved to class 3 with the reason.
//!
//! It is one block of the structural-variants section, beside the SV
//! callers' counts (`sv_callers`): [`text`] and [`html`] write the AnnotSV
//! block and card, and the section's writer puts them after the callers'.

use pycompat::fmt::{general, ljust};

use crate::escape::e;
use crate::model::{Card, Colour, HtmlPart, Inline, Line, Malformed, State, Summary};
use crate::notes::{ANNOTSV_NOTE, VIEWER_NOTE};
use crate::sections::clinvar::{D, VIEWER_TH, locus_text, viewer_td_of};
use crate::{html, text};

const TITLE: &str = "Structural Variants, Class 4 and 5 (AnnotSV)";

/// How many callers saw the SV.
pub fn support_text(x: &D<'_>) -> Result<String, Malformed> {
    let s = x.at("support")?;
    if s.is_none() {
        return Ok("not known (no step 22 consensus)".to_owned());
    }
    if s.is_num(1) {
        return Ok(format!("single caller ({})", x.s("caller")?));
    }
    let callers = x.at("support_callers")?.items()?;
    let names: Vec<&str> = callers.iter().map(D::string).collect::<Result<_, _>>()?;
    Ok(format!("{} callers ({})", s.text()?, names.join(", ")))
}

/// The first three genes, and how many more.
pub fn genes_text(x: &D<'_>) -> Result<String, Malformed> {
    let g = x.at("genes")?;
    if !g.truthy() {
        return Ok("no gene".to_owned());
    }
    let items = g.items()?;
    let first: Vec<&str> = items
        .iter()
        .take(3)
        .map(D::string)
        .collect::<Result<_, _>>()?;
    let mut out = first.join(", ");
    if items.len() > 3 {
        out += &format!(" and {} more", items.len() - 3);
    }
    Ok(out)
}

/// The exon-edge, benign-set and rule flags of one class 4-5 SV.
pub fn flags(x: &D<'_>) -> Result<Vec<String>, Malformed> {
    let mut out = Vec::new();
    let edge = x.at("exon_edge")?;
    if edge.is("yes") {
        out.push(format!(
            "on exon edges: {}",
            x.at("exon_edge_detail")?.string()?
        ));
    } else if edge.is("not checked") {
        out.push("exon edges not checked".to_owned());
    }
    let af = x.at("population_af")?;
    if !af.is_none() {
        let mut f = format!("benign-set AF {}", general(af.num()?, 6));
        let src = x.at("population_source")?;
        if src.truthy() {
            f += &format!(" ({})", src.text()?);
        }
        if x.at("common")?.truthy() {
            f += ", common";
        }
        out.push(f);
    }
    if let Some(detail) = x.opt("rule_detail")? {
        let class = x.s("acmg_class")?;
        if class == "4" || class == "5" {
            out.push(format!("rule {}: {}", x.s("rule")?, detail.text()?));
        }
    }
    Ok(out)
}

/// `v or default`, as text.
fn or_text(x: &D<'_>, key: &str, default: &str) -> Result<String, Malformed> {
    let v = x.at(key)?;
    if v.truthy() {
        v.text()
    } else {
        Ok(default.to_owned())
    }
}

/// One SV step 05's rules moved to class 3: where, from which class, and
/// why.
pub fn moved_text(x: &D<'_>) -> Result<String, Malformed> {
    Ok(format!(
        "{} {}  {}  AnnotSV class {} -> {}  ({}: {})",
        x.s("type")?,
        x.s("sv")?,
        genes_text(x)?,
        x.s("annotsv_class")?,
        x.s("acmg_class")?,
        or_text(x, "rule", "rule not named")?,
        or_text(x, "rule_detail", "no detail")?
    ))
}

/// The class 4-5 lines, the same in both reports.
pub fn lines(d: &D<'_>) -> Result<Vec<String>, Malformed> {
    let count = d.at("count")?;
    let mut head = format!("Class 4-5 SVs: {}", count.text()?);
    if count.truthy() {
        let consensus = d.at("consensus")?.truthy();
        let single = if consensus {
            d.s("single_caller")?
        } else {
            "not known".to_owned()
        };
        let edge = if d.at("exon_edge_checked")?.truthy() {
            d.s("exon_edge")?
        } else {
            "not checked".to_owned()
        };
        head += &format!(
            " (single caller: {single}; deletions on exon edges: {edge}; common: {})",
            d.s("common")?
        );
    }
    let mut out = vec![head];
    if let Some(down) = d.opt("downgraded")? {
        out.push(format!(
            "Moved to class 3 by step 05's rules, not counted above: {}",
            down.len()?
        ));
    }
    if count.truthy() && !d.at("consensus")?.truthy() {
        out.push(match d.opt("consensus_note")? {
            Some(n) => n.text()?,
            None => "Caller support: step 22's consensus did not run, so how many callers saw each SV is \
                     not known."
                .to_owned(),
        });
    }
    Ok(out)
}

/// The AnnotSV block of the text report.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "annotsv", Some(TITLE))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = D::of(&s.section("annotsv"));
    w.extend(lines(&d)?.iter().map(|l| format!("  {l}")));
    for x in d.at("rows")?.items()? {
        w.push(format!(
            "    {} {}  class {}  {}  support: {}",
            ljust(&x.s("type")?, 4),
            x.s("sv")?,
            x.s("acmg_class")?,
            genes_text(&x)?,
            support_text(&x)?
        ));
        for f in flags(&x)? {
            w.push(format!("         {f}"));
        }
    }
    let down = d.opt("downgraded")?;
    if let Some(down) = &down {
        w.push(
            "  Moved to class 3 by step 05's rules (from AnnotSV's class, with the reason):"
                .to_owned(),
        );
        for x in down.items()? {
            w.push(format!("    {}", moved_text(&x)?));
        }
    }
    if d.at("count")?.truthy() || down.is_some() {
        w.push(format!("  {ANNOTSV_NOTE}"));
    }
    w.push(String::new());
    Ok(w)
}

/// `type sv` of a row, as the tables show it.
fn sv_cell(x: &D<'_>) -> Result<String, Malformed> {
    Ok(format!(
        "{} {}",
        x.at("type")?.string()?,
        x.at("sv")?.string()?
    ))
}

/// The AnnotSV card of the HTML report, when the step ran or failed.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("annotsv");
    if sec.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let d = D::of(&sec);
    let done = if sec.state == State::Stale {
        Inline::Badge("Stale".to_owned(), Colour::Yellow)
    } else {
        Inline::Badge("Complete".to_owned(), Colour::Green)
    };
    let mut body = vec![Line::Stat("Status".to_owned(), done)];
    if d.truthy() {
        body.extend(lines(&d)?.into_iter().map(Line::Note));
    }
    let rows = d.opt("rows")?;
    if let Some(rows) = &rows {
        let mut t = format!(
            "    <table><tr><th>SV</th><th>Class</th><th>Genes</th><th>Callers</th><th>Checks</th>{VIEWER_TH}</tr>"
        );
        for x in rows.items()? {
            let checks = flags(&x)?.join("; ");
            t += &format!(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>{}</tr>",
                e(&sv_cell(&x)?),
                e(&x.s("acmg_class")?),
                e(&genes_text(&x)?),
                e(&support_text(&x)?),
                e(if checks.is_empty() { "none" } else { &checks }),
                viewer_td_of(&locus_text(&x.at("sv")?)?)?
            );
        }
        t += "</table>";
        body.push(Line::Raw(t));
    }
    let down = d.opt("downgraded")?;
    if let Some(down) = &down {
        body.push(Line::Raw(
            "    <h3>Moved to class 3 by step 05's rules</h3>".to_owned(),
        ));
        let mut t = format!(
            "    <table><tr><th>SV</th><th>AnnotSV class</th><th>Class</th><th>Genes</th><th>Rule</th>\
             <th>Why</th>{VIEWER_TH}</tr>"
        );
        for x in down.items()? {
            t += &format!(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>{}</tr>",
                e(&sv_cell(&x)?),
                e(&x.s("annotsv_class")?),
                e(&x.s("acmg_class")?),
                e(&genes_text(&x)?),
                e(&or_text(&x, "rule", "not named")?),
                e(&or_text(&x, "rule_detail", "no detail")?),
                viewer_td_of(&locus_text(&x.at("sv")?)?)?
            );
        }
        t += "</table>";
        body.push(Line::Raw(t));
    }
    let mut caveats = Vec::new();
    if rows.is_some() || down.is_some() {
        body.push(Line::Raw(VIEWER_NOTE.to_owned()));
        caveats.push(ANNOTSV_NOTE);
    }
    let card = html::card(
        &sec,
        Card {
            title: TITLE.to_owned(),
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
