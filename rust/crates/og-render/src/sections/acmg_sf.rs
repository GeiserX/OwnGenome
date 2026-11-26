//! Step 23's clinical filter and the ACMG secondary-findings tier: the
//! counts of the filtered variants, then each SF gene's verdict by its
//! reporting rule (`assets/acmg_sf.tsv`): the ClinVar P/LP hits, the rare
//! HIGH-impact records and the class 4/5 SVs it reports, and the variants
//! it does not report, each with the reason. Every phrase comes from the
//! summary as `og summary collect` wrote it.

use pycompat::fmt::ljust;

use crate::model::{Card, Cell, HtmlPart, Inline, Layout, Line, Malformed, Summary, Table};
use crate::notes::ACMG_NOTE;
use crate::sections::annotsv::unknown_reason;
use crate::sections::clinvar::{D, viewer_of};
use crate::{html, text};

/// One row of the SF table: gene, variant, genotype, evidence.
pub type Row = (String, String, String, String);

/// The HIGH-impact records counted as rare: neither a weak call nor of
/// unknown frequency (a summary written before those marks counts every
/// record).
pub fn high_rare(sf: &D<'_>) -> Result<String, Malformed> {
    let all = sf.at("high_impact")?.len()?;
    match sf.get("high_impact_rare")? {
        Some(n) => n.text(),
        None => Ok(all.to_string()),
    }
}

/// The HIGH-impact records listed but not counted as rare.
pub fn high_apart(sf: &D<'_>) -> Result<String, Malformed> {
    let mut parts = Vec::new();
    if let Some(n) = sf.opt("high_impact_weak")? {
        parts.push(format!("{} weak call(s)", n.text()?));
    }
    if let Some(n) = sf.opt("high_impact_frequency_unknown")? {
        parts.push(format!("{} of unknown frequency", n.text()?));
    }
    Ok(if parts.is_empty() {
        String::new()
    } else {
        format!(" (listed apart: {})", parts.join(", "))
    })
}

/// The marks of one HIGH-impact record (unknown frequency, weak call), each
/// after `lead`; `""` for a rare record with a sound call.
pub fn high_marks(x: &D<'_>, lead: &str) -> Result<String, Malformed> {
    let mut marks = Vec::new();
    if x.flag("frequency_unknown")? {
        marks.push("frequency unknown, not counted as rare".to_owned());
    }
    if x.flag("weak_call")? {
        marks.push(format!(
            "weak call ({}): check the reads",
            x.s_or_none("call_note")?
        ));
    }
    if marks.is_empty() {
        return Ok(String::new());
    }
    let close = if lead.ends_with('(') { ")" } else { "" };
    Ok(format!("{lead}{}{close}", marks.join("; ")))
}

/// The severity, rule and two-hit notes of a row.
fn extra(x: &D<'_>) -> Result<String, Malformed> {
    let mut out = String::new();
    for k in ["severity", "rule_note", "two_hit"] {
        if let Some(v) = x.opt(k)? {
            out += &format!("; {}", v.text()?);
        }
    }
    Ok(out)
}

/// Whether step 22's consensus ran, as the AnnotSV section says: then an
/// SV whose support is not known is one the consensus cannot tell.
pub fn consensus_ran(s: &Summary) -> Result<bool, Malformed> {
    let sec = s.section("annotsv");
    let d = D::of(&sec);
    if d.is_none() {
        return Ok(false);
    }
    d.flag("consensus")
}

/// The secondary findings as rows, the same in both reports: ClinVar hits,
/// rare HIGH-impact records, then class 4 and 5 SVs. `consensus` is
/// [`consensus_ran`].
pub fn rows(sf: &D<'_>, consensus: bool) -> Result<Vec<Row>, Malformed> {
    let mut out = Vec::new();
    for x in sf.at("clinvar_hits")?.items()? {
        let low = if x.flag("low_confidence")? {
            "; low-confidence call: check the reads"
        } else {
            ""
        };
        out.push((
            x.s("gene")?,
            x.s("variant")?,
            x.s("genotype")?,
            format!(
                "ClinVar {}{low}{}",
                x.at("significance")?.string()?,
                extra(&x)?
            ),
        ));
    }
    for x in sf.at("high_impact")?.items()? {
        out.push((
            x.s("gene")?,
            x.s("variant")?,
            x.s("genotype")?,
            format!(
                "{}{}{}",
                x.at("consequence")?.string()?,
                high_marks(&x, "; ")?,
                extra(&x)?
            ),
        ));
    }
    if let Some(sv) = sf.get("structural")? {
        for x in sv.items()? {
            let support = x.at("support")?;
            let support = if support.is_none() {
                format!("caller support not known ({})", unknown_reason(consensus))
            } else if support.is_num(1) {
                format!("single caller ({})", x.s("caller")?)
            } else {
                format!("{} callers", support.text()?)
            };
            let via = match x.opt("via")? {
                Some(v) => format!("; through {}", v.text()?),
                None => String::new(),
            };
            out.push((
                x.s("gene")?,
                x.s("variant")?,
                ".".to_owned(),
                format!(
                    "AnnotSV class {}, {support}{via}{}",
                    x.s("acmg_class")?,
                    extra(&x)?
                ),
            ));
        }
    }
    Ok(out)
}

/// The variants in SF genes the gene's rule does not report, or why no
/// rule could be applied.
pub fn apart_lines(sf: &D<'_>) -> Result<Vec<String>, Malformed> {
    if sf.flag("unavailable")? {
        return Ok(vec![format!(
            "Secondary findings: not assessed, {}",
            sf.s("unavailable")?
        )]);
    }
    let mut out = Vec::new();
    if let Some(nr) = sf.get("not_reported")? {
        for x in nr.items()? {
            let gt = match x.opt("genotype")? {
                Some(g) => g.text()?,
                None => ".".to_owned(),
            };
            out.push(format!(
                "{} {} {gt} ({}): not reported, {}",
                x.s("gene")?,
                x.s("variant")?,
                x.s("source")?,
                x.s("reason")?
            ));
        }
    }
    if let Some(n) = sf.opt("two_hit_note")? {
        out.push(n.text()?);
    }
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "clinical", Some("Clinical Variant Filter"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = D::of(&s.section("clinical"));
    let genes = if d.has("genes") {
        format!(" in {} genes", d.s("genes")?)
    } else {
        String::new()
    };
    w.push(format!("  Clinical variants: {}{genes}", d.s("variants")?));
    if let Some(by) = d.opt("by_impact")? {
        let mut pairs = by.pairs()?;
        pairs.sort_by(|a, b| a.0.cmp(b.0));
        let parts: Vec<String> = pairs
            .iter()
            .map(|(k, v)| Ok(format!("{k} {}", v.text()?)))
            .collect::<Result<_, Malformed>>()?;
        w.push(format!("  By impact: {}", parts.join(", ")));
    }
    if d.has("frequency_unknown") {
        w.push(format!(
            "  Frequency unknown (not counted as rare): {}",
            d.s("frequency_unknown")?
        ));
    }
    if d.has("weak_calls") {
        w.push(format!(
            "  Weak genotype calls (check the reads): {}",
            d.s("weak_calls")?
        ));
    }
    if let Some(sf) = d.opt("acmg_sf")? {
        if sf.flag("unavailable")? {
            w.push(format!("  {}", apart_lines(&sf)?[0]));
        } else {
            let mut line = format!(
                "  Secondary-findings genes ({}, {} genes): {} ClinVar P/LP, {} rare HIGH impact{}",
                sf.s("version")?,
                sf.s_or("genes_on_list", "?")?,
                sf.at("clinvar_hits")?.len()?,
                high_rare(&sf)?,
                high_apart(&sf)?
            );
            if sf.flag("structural_checked")? {
                line += &format!(", {} class 4/5 SV(s)", sf.at("structural")?.len()?);
            }
            if let Some(nr) = sf.opt("not_reported")? {
                line += &format!(", {} not reported by the gene's rule", nr.len()?);
            }
            w.push(line);
            for (g, v, z, ev) in rows(&sf, consensus_ran(s)?)? {
                w.push(format!("    {} {v}  {z}  {ev}", ljust(&g, 8)));
            }
            for line in apart_lines(&sf)? {
                w.push(format!("    {line}"));
            }
            w.push(format!("    {ACMG_NOTE}"));
        }
    }
    w.push(String::new());
    Ok(w)
}

/// The counts card, always there; a stat's missing count prints `None`.
fn counts_card(d: &D<'_>) -> Result<Card, Malformed> {
    let stat = |label: String, v: String| Line::Stat(label, Inline::Text(v));
    let mut body = vec![stat(
        "Total interesting variants".to_owned(),
        d.s_or_none("variants")?,
    )];
    if d.has("genes") {
        body.push(stat("Genes".to_owned(), d.s("genes")?));
    }
    if let Some(by) = d.opt("by_impact")? {
        let mut pairs = by.pairs()?;
        pairs.sort_by(|a, b| a.0.cmp(b.0));
        for (k, v) in pairs {
            body.push(stat(format!("{k} impact"), v.text()?));
        }
    }
    if d.has("frequency_unknown") {
        body.push(stat(
            "Frequency unknown (not rare)".to_owned(),
            d.s("frequency_unknown")?,
        ));
    }
    if d.has("weak_calls") {
        body.push(stat("Weak genotype calls".to_owned(), d.s("weak_calls")?));
    }
    Ok(Card {
        title: "Clinical Variant Filter".to_owned(),
        body,
        ..Card::default()
    })
}

/// The SF card, when the tier lists anything or says why it could not.
fn sf_card(sf: &D<'_>, consensus: bool) -> Result<Option<Card>, Malformed> {
    let shown = sf.at("clinvar_hits")?.truthy()
        || sf.at("high_impact")?.truthy()
        || sf.flag("structural")?
        || sf.flag("not_reported")?
        || sf.flag("unavailable")?;
    if !shown {
        return Ok(None);
    }
    let mut cells = Vec::new();
    for (g, v, z, ev) in rows(sf, consensus)? {
        let at = viewer_of(&v)?;
        cells.push(vec![
            Cell::Text(g),
            Cell::Text(v),
            Cell::Text(z),
            Cell::Text(ev),
            at,
        ]);
    }
    let listed = !cells.is_empty();
    let mut body = vec![Line::Table(Table {
        head: ["Gene", "Variant", "Genotype", "Evidence"]
            .map(str::to_owned)
            .to_vec(),
        rows: cells,
        viewer: true,
        layout: Layout::RowLines,
    })];
    body.extend(apart_lines(sf)?.into_iter().map(Line::Para));
    if listed {
        body.push(Line::ViewerNote);
    }
    Ok(Some(Card {
        title: format!("Secondary-Findings Genes ({})", sf.s("version")?),
        full: true,
        caveats: vec![ACMG_NOTE],
        body,
        ..Card::default()
    }))
}

/// The HTML report's cards: the counts, then the SF table.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("clinical");
    let d = D::of(&sec);
    let mut cards = vec![html::card(&sec, counts_card(&d)?)?];
    if let Some(sf) = d.opt("acmg_sf")?
        && let Some(card) = sf_card(&sf, consensus_ran(s)?)?
    {
        cards.push(html::card(&sec, card)?);
    }
    Ok(HtmlPart {
        cards,
        files: Vec::new(),
    })
}
