//! Step 23's clinical filter and the ACMG secondary-findings tier: the
//! counts of the filtered variants, then each SF gene's verdict by its
//! reporting rule (`assets/acmg_sf.tsv`): the ClinVar P/LP hits, the rare
//! HIGH-impact records and the class 4/5 SVs it reports, and the variants
//! it does not report, each with the reason. Every phrase comes from the
//! summary as `og summary collect` wrote it.
//!
//! Each row also carries the marks the collector writes when they apply:
//! the condition ClinVar names for the hit, and whether it is another
//! disease than the one ACMG lists the gene for; the loss-of-function
//! mechanism note; the TTN and LMNA rule condition; a HIGH-impact call
//! with nothing to judge it by; and the pair word of a recessive gene with
//! its phase and note. A summary without those keys prints the rows as
//! before. A run whose step 23 wrote no readable clinical VCF has no
//! `clinical.acmg_sf`; its card is then `sections.acmg_sf`, built from
//! step 6 alone, and says the HIGH-impact source was not checked.

use pycompat::fmt::ljust;

use crate::model::{
    Card, Cell, HtmlPart, Inline, Layout, Line, Malformed, Section, Summary, Table,
};
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
    if let Some(n) = sf.opt("high_impact_unverified")? {
        parts.push(format!(
            "{} with no GQ, DP or AD to judge the call",
            n.text()?
        ));
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
    if x.flag("call_unverified")? {
        marks.push("no GQ, DP or AD to judge the call, not counted as rare".to_owned());
    }
    if marks.is_empty() {
        return Ok(String::new());
    }
    let close = if lead.ends_with('(') { ")" } else { "" };
    Ok(format!("{lead}{}{close}", marks.join("; ")))
}

/// The severity, rule, rule condition, mechanism and two-hit notes of a
/// row, then the pair word of a recessive gene with its phase and note.
fn extra(x: &D<'_>) -> Result<String, Malformed> {
    let mut out = String::new();
    for k in [
        "severity",
        "rule_note",
        "rule_condition",
        "mechanism_note",
        "two_hit",
    ] {
        if let Some(v) = x.opt(k)? {
            out += &format!("; {}", v.text()?);
        }
    }
    if let Some(pair) = x.opt("pair")? {
        let pair = pair.text()?;
        let phase = match x.opt("phase")? {
            Some(p) => p.text()?,
            None => String::new(),
        };
        // The row's own two-hit note can already say it, word and phase
        // alike (a deletion in trans by the same phase set): say it once.
        let said = match x.opt("two_hit")? {
            Some(t) => {
                let t = t.text()?;
                t.starts_with(&pair) && t.ends_with(&phase)
            }
            None => false,
        };
        if !said {
            out += &format!("; {pair}");
            if !phase.is_empty() {
                out += &format!(": {phase}");
            }
        }
    }
    if let Some(n) = x.opt("pair_note")? {
        out += &format!("; {}", n.text()?);
    }
    Ok(out)
}

/// The condition ClinVar names for a hit, and the mark when it is another
/// disease than the one ACMG lists the gene for; `""` for a hit with no
/// named condition.
fn condition(x: &D<'_>) -> Result<String, Malformed> {
    let Some(c) = x.opt("clinvar_condition")? else {
        return Ok(String::new());
    };
    let mut out = format!("; ClinVar condition: {}", c.text()?);
    if x.flag("condition_not_acmg")? {
        let listed = match x.opt("acmg_phenotype")? {
            Some(p) => format!(" ({})", p.text()?),
            None => String::new(),
        };
        out += &format!(
            ", not the disease ACMG lists {} for{listed}: check which disease the variant was \
             reported for",
            x.s("gene")?
        );
    }
    Ok(out)
}

/// The disease ACMG lists `gene` for, from the first row of the gene that
/// carries it (`acmg_phenotype`, on a ClinVar hit whose condition step 6
/// copied).
pub fn listed_condition(sf: &D<'_>, gene: &str) -> Result<Option<String>, Malformed> {
    for k in ["clinvar_hits", "high_impact", "structural"] {
        let Some(list) = sf.get(k)? else { continue };
        if !matches!(list.v, pycompat::json::Value::Array(_)) {
            continue;
        }
        for x in list.items()? {
            if x.get("gene")?.is_some_and(|g| g.is(gene))
                && let Some(p) = x.opt("acmg_phenotype")?
            {
                return Ok(Some(p.text()?));
            }
        }
    }
    Ok(None)
}

/// The opening sentence of a secondary finding, in plain words: what the
/// gene is, why it is reported, and what the result is not. `counted` is
/// false for a row the card lists apart (see [`counted`]): its call or its
/// frequency needs checking before the rule counts it. The list's name,
/// its version and the ClinVar class stay in the evidence rows.
pub fn lede(gene: &str, condition: Option<&str>, counted: bool) -> String {
    let health = match condition {
        Some(c) => format!("a change there can matter for health ({c})"),
        None => "a change there can matter for health".to_owned(),
    };
    let found = if counted {
        "This genome has a change there that the gene's reporting rule counts."
    } else {
        "This genome may have a change there that counts: its call or its frequency needs \
         checking first."
    };
    format!(
        "{gene} is one of the genes doctors report even when nobody looked for it, because \
         {health}. {found} This is a research result, not a diagnosis: a certified laboratory \
         confirms it before anyone acts on it."
    )
}

/// Whether the gene's rule counts each row of [`rows`], in the same order:
/// false for a HIGH-impact record the card lists apart (a weak call, a call
/// with no GQ, DP or AD, or unknown frequency), true for every other row.
pub fn counted(sf: &D<'_>) -> Result<Vec<bool>, Malformed> {
    let mut out = vec![true; sf.at("clinvar_hits")?.len()?];
    for x in sf.at("high_impact")?.items()? {
        out.push(high_marks(&x, "")?.is_empty());
    }
    if let Some(sv) = sf.get("structural")? {
        out.extend(std::iter::repeat_n(true, sv.len()?));
    }
    Ok(out)
}

/// The card's label when step 23's HIGH-impact source was not read.
pub const CLINVAR_ONLY: &str = "ClinVar only, HIGH-impact not checked";

/// Whether the card holds the ClinVar hits and SVs only: step 23 left no
/// readable clinical VCF (`high_impact_checked: false`). A summary written
/// before the key reads as checked.
pub fn clinvar_only(sf: &D<'_>) -> Result<bool, Malformed> {
    Ok(sf
        .get("high_impact_checked")?
        .is_some_and(|c| matches!(c.v, pycompat::json::Value::Bool(false))))
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
                "ClinVar {}{low}{}{}",
                x.at("significance")?.string()?,
                condition(&x)?,
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

/// The card of `sections.acmg_sf`, built from step 6 alone, when the
/// clinical section holds none (step 23 did not run, failed, or wrote no
/// readable clinical VCF). `None` when the clinical section holds the card
/// or the section did not run.
pub fn standalone(s: &Summary) -> Result<Option<Section<'_>>, Malformed> {
    let clinical = s.section("clinical");
    if clinical.state.ran() && D::of(&clinical).opt("acmg_sf")?.is_some() {
        return Ok(None);
    }
    let sec = s.section("acmg_sf");
    Ok((sec.state.ran() && D::of(&sec).truthy()).then_some(sec))
}

/// The text lines of the SF tier, indented as the clinical block indents
/// them.
fn sf_text(sf: &D<'_>, s: &Summary) -> Result<Vec<String>, Malformed> {
    if sf.flag("unavailable")? {
        return Ok(vec![format!("  {}", apart_lines(sf)?[0])]);
    }
    let only = clinvar_only(sf)?;
    let label = if only {
        format!("; {CLINVAR_ONLY}")
    } else {
        String::new()
    };
    let mut line = format!(
        "  Secondary-findings genes ({}, {} genes{label}): {} ClinVar P/LP",
        sf.s("version")?,
        sf.s_or("genes_on_list", "?")?,
        sf.at("clinvar_hits")?.len()?,
    );
    if !only {
        line += &format!(", {} rare HIGH impact{}", high_rare(sf)?, high_apart(sf)?);
    }
    if sf.flag("structural_checked")? {
        line += &format!(", {} class 4/5 SV(s)", sf.at("structural")?.len()?);
    }
    if let Some(nr) = sf.opt("not_reported")? {
        line += &format!(", {} not reported by the gene's rule", nr.len()?);
    }
    let mut w = vec![line];
    for (g, v, z, ev) in rows(sf, consensus_ran(s)?)? {
        w.push(format!("    {} {v}  {z}  {ev}", ljust(&g, 8)));
    }
    for line in apart_lines(sf)? {
        w.push(format!("    {line}"));
    }
    if let Some(n) = sf.opt("high_impact_note")? {
        w.push(format!("    {}", n.text()?));
    }
    w.push(format!("    {ACMG_NOTE}"));
    Ok(w)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "clinical", Some("Clinical Variant Filter"))? else {
        // No clinical section at all: the card from step 6 alone, under its
        // own heading.
        let Some(sec) = standalone(s)? else {
            return Ok(Vec::new());
        };
        let Some(head) = text::head(s, "acmg_sf", Some("Secondary-Findings Genes"))? else {
            return Ok(Vec::new());
        };
        let mut w = head.lines;
        w.extend(sf_text(&D::of(&sec), s)?);
        w.push(String::new());
        return Ok(w);
    };
    let mut w = head.lines;
    if !head.ran {
        if let Some(sec) = standalone(s)? {
            w.extend(sf_text(&D::of(&sec), s)?);
            w.push(String::new());
        }
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
        w.extend(sf_text(&sf, s)?);
    } else if let Some(sec) = standalone(s)? {
        w.extend(sf_text(&D::of(&sec), s)?);
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
        || sf.flag("unavailable")?
        || sf.flag("high_impact_note")?
        || clinvar_only(sf)?;
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
    if let Some(n) = sf.opt("high_impact_note")? {
        body.push(Line::Para(n.text()?));
    }
    if listed {
        body.push(Line::ViewerNote);
    }
    Ok(Some(Card {
        title: if clinvar_only(sf)? {
            format!(
                "Secondary-Findings Genes ({}; {CLINVAR_ONLY})",
                sf.s("version")?
            )
        } else {
            format!("Secondary-Findings Genes ({})", sf.s("version")?)
        },
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
    // The same card the text report shows: the clinical section's while it
    // ran and holds one, else the one built from step 6 alone.
    if let Some(alone) = standalone(s)? {
        if let Some(card) = sf_card(&D::of(&alone), consensus_ran(s)?)? {
            cards.push(html::card(&alone, card)?);
        }
    } else if let Some(sf) = d.opt("acmg_sf")?
        && let Some(card) = sf_card(&sf, consensus_ran(s)?)?
    {
        cards.push(html::card(&sec, card)?);
    }
    Ok(HtmlPart {
        cards,
        files: Vec::new(),
    })
}
