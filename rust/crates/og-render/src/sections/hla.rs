//! HLA typing (step 08): T1K's alleles per locus with their quality, the
//! lines derived from them (the coeliac DQ types, the DQB1*02 dose and the
//! allele table of bin/hla_alleles.tsv), OptiType's second opinion on the
//! class I genes, and the KIR genes.

use pycompat::json::Value;

use crate::html;
use crate::model::{Card, Inline, Line, Malformed, Summary, get, is_str, need, num, py, truthy};
use crate::notes::{HLA_DERIVED_NOTE, HLA_NOTE, KIR_NOTE};
use crate::sections::pgx::{items, left, or_text, text_of};
use crate::text;

/// A quality-0 allele with more candidates than this reads "no confident
/// call" (og summary collect keeps this many of them).
pub const MAX_CANDIDATES: f64 = 3.0;

/// T1K's quality as a number; one that is not a number counts as -1, as in
/// step 36.
fn quality(q: &Value) -> Result<f64, Malformed> {
    match q {
        Value::Str(s) => Ok(pycompat::num::parse_float(s).unwrap_or(-1.0)),
        _ => num(Some(q)).ok_or_else(|| Malformed("an HLA quality is not a number".into())),
    }
}

/// One HLA locus: each allele with its T1K quality, and the low-confidence
/// mark. A quality-0 allele with more than [`MAX_CANDIDATES`] candidates
/// reads "no confident call"; the HFE row carries its label.
pub fn hla_text(locus: &Value) -> Result<String, Malformed> {
    let p = "sections.hla.data.loci[]";
    let q = items(get(locus, "quality", p)?.filter(|v| truthy(v)), p)?;
    let n = items(get(locus, "candidates", p)?.filter(|v| truthy(v)), p)?;
    let mut parts = Vec::new();
    for (i, a) in items(Some(need(locus, "alleles", p)?), p)?
        .iter()
        .enumerate()
    {
        let blank = Value::Str(String::new());
        let qi = q.get(i).unwrap_or(&blank);
        let has_q = !is_str(Some(qi), "");
        let many = match n.get(i) {
            Some(c) => {
                num(Some(c))
                    .ok_or_else(|| Malformed("an HLA candidate count is not a number".into()))?
                    > MAX_CANDIDATES
            }
            None => false,
        };
        if many && has_q && quality(qi)? <= 0.0 {
            parts.push(format!(
                "no confident call ({} candidates, quality {})",
                py(&n[i])?,
                py(qi)?
            ));
        } else {
            let mut t = text_of(a, "an HLA allele")?;
            if has_q {
                t += &format!(" (quality {})", py(qi)?);
            }
            parts.push(t);
        }
    }
    let mut txt = parts.join(" / ");
    if txt.is_empty() {
        txt = "no call".to_owned();
    }
    if get(locus, "low_confidence", p)?.is_some_and(truthy) {
        txt += "; low confidence (T1K quality 0 or below)";
        if get(locus, "withheld_from_pharmcat", p)?.is_some_and(truthy) {
            txt += ", not passed to PharmCAT";
        }
    }
    if let Some(l) = get(locus, "label", p)?.filter(|v| truthy(v)) {
        txt += &format!(" [{}]", py(l)?);
    }
    Ok(txt)
}

/// The order of an allele's state in the table: present ones first.
fn rank(state: &Value) -> u8 {
    match state {
        Value::Str(s) => match s.as_str() {
            "present" => 0,
            "uncertain" => 1,
            "absent" => 2,
            "not_determined" => 3,
            _ => 9,
        },
        _ => 9,
    }
}

/// The derived HLA lines both reports print: the coeliac DQ line, the
/// DQB1*02 dose and the allele table, present ones first.
pub fn hla_derived_lines(d: Option<&Value>) -> Result<Vec<String>, Malformed> {
    let Some(d) = d.filter(|v| truthy(v)) else {
        return Ok(Vec::new());
    };
    let p = "sections.hla.data.derived";
    let empty = Value::Object(Default::default());
    let c = get(d, "coeliac", p)?
        .filter(|v| truthy(v))
        .unwrap_or(&empty);
    let cp = "sections.hla.data.derived.coeliac";
    let mut out = Vec::new();
    if is_str(get(c, "status", cp)?, "derived") {
        let mut dq25 = text_of(need(c, "dq2_5", cp)?, "the DQ2.5 state")?;
        if let Some(ph) = get(c, "dq2_5_phase", cp)?.filter(|v| truthy(v)) {
            dq25 += &format!(", {}", py(ph)?);
        }
        out.push(format!(
            "Coeliac DQ: DQ2.5 {dq25}; DQ2.2 {}; DQ8 {}",
            py(need(c, "dq2_2", cp)?)?,
            py(need(c, "dq8", cp)?)?
        ));
    } else {
        out.push(format!(
            "Coeliac DQ: not determined ({})",
            or_text(get(c, "reason", cp)?, "no DQA1 and DQB1 call")?
        ));
    }
    let dose = match get(c, "dqb1_02_dose", cp)? {
        None | Some(Value::Null) => "not determined".to_owned(),
        Some(v) => py(v)?,
    };
    out.push(format!("DQB1*02 copies: {dose}"));
    let mut alleles: Vec<&Value> = items(
        get(d, "alleles", p)?.filter(|v| truthy(v)),
        "sections.hla.data.derived.alleles",
    )?
    .iter()
    .collect();
    let ap = "sections.hla.data.derived.alleles[]";
    for a in &alleles {
        need(a, "state", ap)?;
    }
    alleles.sort_by_key(|a| rank(need(a, "state", ap).expect("checked above")));
    for a in alleles {
        let mut state = text_of(need(a, "state", ap)?, "an allele state")?.replace('_', " ");
        if let Some(r) = get(a, "reason", ap)?.filter(|v| truthy(v)) {
            state += &format!(" ({})", py(r)?);
        }
        out.push(format!(
            "{} {state}: {} ({})",
            py(need(a, "allele", ap)?)?,
            py(need(a, "association", ap)?)?,
            py(need(a, "citation", ap)?)?
        ));
    }
    Ok(out)
}

/// Step 08's OptiType check: (gene, line) per class I gene, or one
/// (`""`, line) with the reason the callers were not compared.
pub fn hla_optitype_lines(ot: &Value) -> Result<Vec<(String, String)>, Malformed> {
    let p = "sections.hla.data.optitype";
    let genes = need(ot, "genes", p)?;
    if !truthy(genes) {
        return Ok(vec![(
            String::new(),
            format!("not compared: {}", py(need(ot, "reason", p)?)?),
        )]);
    }
    let gp = "sections.hla.data.optitype.genes[]";
    let mut out = Vec::new();
    for g in items(Some(genes), p)? {
        let agree = need(g, "agree", gp)?;
        let word = match agree {
            Value::Str(s) if s == "yes" => "agree".to_owned(),
            Value::Str(s) if s == "no" => "disagree".to_owned(),
            _ => or_text(Some(agree), "not compared")?,
        };
        let mut line = format!(
            "{word}: T1K {}, OptiType {}",
            py(need(g, "t1k", gp)?)?,
            py(need(g, "optitype", gp)?)?
        );
        let note = need(g, "note", gp)?;
        if !(is_str(Some(note), "") || is_str(Some(note), "-")) {
            line += &format!(" ({})", py(note)?);
        }
        out.push((py(need(g, "gene", gp)?)?, line));
    }
    Ok(out)
}

/// Step 08's KIR genes: (gene, alleles) per typed gene, then (`""`, line)
/// for the genes with no allele, or the reason nothing was typed.
pub fn kir_lines(kir: &Value) -> Result<Vec<(String, String)>, Malformed> {
    let p = "sections.hla.data.kir";
    let reason = need(kir, "reason", p)?;
    if truthy(reason) {
        return Ok(vec![(String::new(), py(reason)?)]);
    }
    let gp = "sections.hla.data.kir.typed[]";
    let mut out = Vec::new();
    for g in items(Some(need(kir, "typed", p)?), p)? {
        let alleles = items(Some(need(g, "alleles", gp)?), gp)?;
        let quals = items(Some(need(g, "quality", gp)?), gp)?;
        let mut parts = Vec::new();
        for (a, q) in alleles.iter().zip(quals) {
            let mut t = text_of(a, "a KIR allele")?;
            if !is_str(Some(q), "") {
                t += &format!(" (quality {})", py(q)?);
            }
            parts.push(t);
        }
        out.push((py(need(g, "gene", gp)?)?, parts.join(" / ")));
    }
    let not_called = need(kir, "not_called", p)?;
    if truthy(not_called) {
        let names: Result<Vec<String>, Malformed> = items(Some(not_called), p)?
            .iter()
            .map(|x| text_of(x, "a KIR gene"))
            .collect();
        out.push((
            String::new(),
            format!(
                "No allele called (absent, or too few reads): {}",
                names?.join(", ")
            ),
        ));
    }
    Ok(out)
}

/// `    GENE     line`, or the line alone without a gene.
fn gene_line(g: &str, line: &str) -> String {
    if g.is_empty() {
        format!("    {line}")
    } else {
        format!("    {} {line}", pycompat::fmt::ljust(g, 8))
    }
}

fn database(s: &Summary) -> Result<String, Malformed> {
    or_text(Some(s.need2("databases", "hla_database")?), "unknown")
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "hla", Some("HLA Typing (T1K)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let sec = s.section("hla");
        for l in items(Some(sec.need("loci")?), "sections.hla.data.loci")? {
            w.push(format!(
                "    {} {}",
                left(
                    need(l, "gene", "sections.hla.data.loci[]")?,
                    8,
                    "an HLA gene"
                )?,
                hla_text(l)?
            ));
        }
        let derived = sec.get("derived")?;
        w.extend(
            hla_derived_lines(derived)?
                .into_iter()
                .map(|l| format!("  {l}")),
        );
        if derived.is_some_and(truthy) {
            w.push(format!("  {HLA_DERIVED_NOTE}"));
        }
        w.push(format!("  HLA database: {}", database(s)?));
        w.push(format!("  {HLA_NOTE}"));
        if let Some(ot) = sec.get("optitype")?.filter(|v| **v != Value::Null) {
            w.push("  OptiType check (HLA_QC=true, two fields):".to_owned());
            for (g, line) in hla_optitype_lines(ot)? {
                w.push(gene_line(&g, &line));
            }
        }
        if let Some(kir) = sec.get("kir")?.filter(|v| **v != Value::Null) {
            w.push("  KIR genes (T1K, KIR=true):".to_owned());
            for (g, line) in kir_lines(kir)? {
                w.push(gene_line(&g, &line));
            }
            w.push(format!("  {KIR_NOTE}"));
        }
        w.push(String::new());
    }
    Ok(w)
}

fn stat(label: &str, value: String) -> Line {
    Line::Stat(label.to_owned(), Inline::Text(value))
}

/// The HTML report's card.
pub fn html(s: &Summary) -> Result<crate::model::HtmlPart, Malformed> {
    let sec = s.section("hla");
    let mut body = Vec::new();
    let mut caveats = Vec::new();
    for l in items(sec.get("loci")?, "sections.hla.data.loci")? {
        body.push(stat(
            &py(need(l, "gene", "sections.hla.data.loci[]")?)?,
            hla_text(l)?,
        ));
    }
    let derived = sec.get("derived")?.filter(|v| truthy(v));
    if derived.is_some() {
        body.extend(hla_derived_lines(derived)?.into_iter().map(Line::Note));
        caveats.push(HLA_DERIVED_NOTE);
    }
    body.push(stat("HLA database", database(s)?));
    caveats.push(HLA_NOTE);
    if let Some(ot) = sec.get("optitype")?.filter(|v| **v != Value::Null) {
        for (g, line) in hla_optitype_lines(ot)? {
            body.push(if g.is_empty() {
                stat("OptiType check", line)
            } else {
                stat(&format!("{g}: T1K vs OptiType"), line)
            });
        }
    }
    if let Some(kir) = sec.get("kir")?.filter(|v| **v != Value::Null) {
        for (g, line) in kir_lines(kir)? {
            body.push(if g.is_empty() {
                Line::Note(line)
            } else {
                stat(&g, line)
            });
        }
        caveats.push(KIR_NOTE);
    }
    let card = html::card(
        &sec,
        Card {
            title: "HLA Typing".into(),
            caveats,
            body,
            ..Card::default()
        },
    )?;
    Ok(crate::model::HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
