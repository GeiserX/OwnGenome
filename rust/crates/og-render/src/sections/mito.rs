//! Mitochondrial variants (step 20): the Mutect2 chrM calls with the D-loop
//! poly-C tracts, ClinVar and gnomAD on each call, the confirmed disease
//! variants of assets/mito_confirmed.tsv, eKLIPse's large deletions and the
//! mtDNA copies per cell.

use pycompat::json::Value;

use crate::html;
use crate::model::{
    Card, Colour, HtmlPart, Inline, Layout, Line, Malformed, Section, Summary, Table, get, is_str,
    need, num, py, truthy,
};
use crate::notes::{MITO_COPIES_NOTE, MITO_DELETION_NOTE, MITO_NOTE, MITO_POLYC_NOTE};
use crate::sections::pgx::{items, or_none, or_text, text_of};
use crate::text;

/// `"."` for a value that is missing, empty or `.`, else its text.
fn val(x: Option<&Value>) -> Result<String, Malformed> {
    match x {
        None | Some(Value::Null) => Ok(".".into()),
        Some(Value::Str(s)) if s.is_empty() || s == "." => Ok(".".into()),
        Some(v) => py(v),
    }
}

/// One annotated chrM call as a line.
pub fn mito_variant_text(v: &Value) -> Result<String, Malformed> {
    let p = "sections.mito.data.annotation.variants[]";
    let g = |k: &str| get(v, k, p);
    let cv = val(g("clinvar_significance")?)?;
    let cv = if cv == "." {
        "not in ClinVar".to_owned()
    } else {
        format!(
            "ClinVar {cv} ({}; {})",
            val(g("clinvar_review")?)?,
            val(g("clinvar_condition")?)?
        )
    };
    let gn = if val(g("gnomad_af_hom")?)? == "." {
        "not in gnomAD".to_owned()
    } else {
        format!(
            "gnomAD AF hom {}, het {}",
            py(need(v, "gnomad_af_hom", p)?)?,
            val(g("gnomad_af_het")?)?
        )
    };
    let polyc = if g("polyc")?.is_some_and(truthy) {
        "  poly-C tract"
    } else {
        ""
    };
    Ok(format!(
        "chrM:{} {}>{}  AF {}  {cv}  {gn}{polyc}",
        py(need(v, "pos", p)?)?,
        py(need(v, "ref", p)?)?,
        py(need(v, "alt", p)?)?,
        val(g("af")?)?
    ))
}

/// The annotation block, `{}` when the summary has none.
fn annotation<'a>(d: &Section<'a>) -> Result<Option<&'a Value>, Malformed> {
    Ok(d.get("annotation")?.filter(|v| truthy(v)))
}

/// Step 20's lines under the counts, the same in both reports.
pub fn mito_lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let mut lines = Vec::new();
    let pc = items(
        d.get("polyc_heteroplasmic")?.filter(|v| truthy(v)),
        "sections.mito.data.polyc_heteroplasmic",
    )?;
    if d.get("polyc_heteroplasmic")?.is_some() {
        let tracts: Result<Vec<String>, Malformed> = items(
            d.get("polyc_tracts")?.filter(|v| truthy(v)),
            "sections.mito.data.polyc_tracts",
        )?
        .iter()
        .map(|t| text_of(t, "a poly-C tract"))
        .collect();
        lines.push(format!(
            "Length heteroplasmy in the D-loop poly-C tracts ({}; common, usually technical): {}",
            tracts?.join(", "),
            pc.len()
        ));
        for x in pc {
            lines.push(format!("  {}", polyc_text(x)?));
        }
    }
    if let Some(a) = annotation(d)? {
        let ap = "sections.mito.data.annotation";
        if let Some(src) = get(a, "source", ap)?.filter(|v| truthy(v)) {
            lines.push(format!(
                "Annotation ({}): ClinVar {}; gnomAD {}",
                py(src)?,
                or_text(get(a, "clinvar", ap)?, "unknown")?,
                or_text(get(a, "gnomad", ap)?, "unknown")?
            ));
            lines.push(format!(
                "ClinVar pathogenic or likely pathogenic: {}",
                get(a, "clinvar_plp", ap)?.map_or(Ok("0".into()), py)?
            ));
            let vs = get(a, "variants", ap)?.filter(|v| truthy(v));
            if vs.is_some() {
                lines.push(
                    "Heteroplasmic, or in ClinVar with a class other than benign:".to_owned(),
                );
                for v in items(vs, ap)? {
                    lines.push(format!("  {}", mito_variant_text(v)?));
                }
            }
        } else if let Some(n) = get(a, "note", ap)?.filter(|v| truthy(v)) {
            lines.push(format!("ClinVar and gnomAD: not annotated ({})", py(n)?));
        }
    }
    lines.extend(mito_extended_lines(d)?);
    Ok(lines)
}

/// One poly-C length heteroplasmy: `variant  AF af`.
fn polyc_text(x: &Value) -> Result<String, Malformed> {
    let p = "sections.mito.data.polyc_heteroplasmic[]";
    Ok(format!(
        "{}  AF {}",
        py(need(x, "variant", p)?)?,
        py(need(x, "af", p)?)?
    ))
}

/// One confirmed disease variant.
pub fn mito_confirmed_text(v: &Value) -> Result<String, Malformed> {
    let p = "sections.mito.data.confirmed[]";
    let pmids = text_of(need(v, "pmids", p)?, "a confirmed variant's PubMed ids")?;
    Ok(format!(
        "{} ({})  AF {}  filter {}  {}; ClinVar {}; PMID {}. {}",
        py(need(v, "hgvs", p)?)?,
        py(need(v, "gene", p)?)?,
        or_text(get(v, "af", p)?, ".")?,
        or_text(get(v, "filter", p)?, ".")?,
        py(need(v, "disease", p)?)?,
        py(need(v, "clinvar_id", p)?)?,
        pmids.replace(';', ", "),
        py(need(v, "note", p)?)?
    ))
}

/// `x + d` as Python adds a number to an int or a float.
fn plus(x: &Value, d: i128) -> Result<String, Malformed> {
    match x {
        Value::Int(i) => Ok((i + d).to_string()),
        Value::Bool(b) => Ok((i128::from(*b) + d).to_string()),
        Value::Float(f) => py(&Value::Float(f + d as f64)),
        _ => Err(Malformed("a deletion breakpoint is not a number".into())),
    }
}

/// One large deletion: its span (across the origin when the 3' breakpoint
/// comes first), its level, and the confirmed deletion it matches.
pub fn mito_deletion_text(r: &Value) -> Result<String, Malformed> {
    let p = "sections.mito.data.deletions.listed[]";
    let (b5, b3) = (need(r, "bkp5", p)?, need(r, "bkp3", p)?);
    let (n5, n3) = match (num(Some(b5)), num(Some(b3))) {
        (Some(a), Some(b)) => (a, b),
        _ => return Err(Malformed("a deletion breakpoint is not a number".into())),
    };
    let span = if n3 > n5 {
        format!(
            "m.{}_{}del ({} bp)",
            plus(b5, 1)?,
            plus(b3, -1)?,
            py(need(r, "length", p)?)?
        )
    } else {
        format!(
            "m.{} to m.{} across the origin ({} bp)",
            plus(b5, 1)?,
            plus(b3, -1)?,
            or_none(get(r, "length", p)?)?
        )
    };
    let mut out = format!("{span}  level {}%", py(need(r, "level_pct", p)?)?);
    if let Some(k) = get(r, "known", p)?.filter(|v| truthy(v)) {
        let kp = "sections.mito.data.deletions.listed[].known";
        let pmids = text_of(need(k, "pmids", kp)?, "a deletion's PubMed ids")?;
        out += &format!(
            "  {}; ClinVar {} ({}); PMID {}",
            py(need(k, "disease", kp)?)?,
            py(need(k, "clinvar_id", kp)?)?,
            py(need(k, "hgvs", kp)?)?,
            pmids.replace(';', ", ")
        );
    }
    Ok(out)
}

/// The copy-number line: the copies per cell, or why there are none.
pub fn mito_copies_text(c: &Value) -> Result<String, Malformed> {
    let p = "sections.mito.data.copies";
    if is_str(get(c, "status", p)?, "called") {
        return Ok(format!(
            "mtDNA copy number: {} per cell (chrM depth {} times the autosomal depth)",
            py(need(c, "mtdna_copies", p)?)?,
            py(need(c, "depth_ratio", p)?)?
        ));
    }
    let reason = get(c, "reason", p)?;
    let why = if is_str(reason, "low_autosomal_depth") {
        "too few autosomal reads: a sliced BAM, a panel or a shallow genome".to_owned()
    } else if is_str(reason, "no_chrM_reads") {
        "no chrM reads in the BAM".to_owned()
    } else {
        or_text(reason, "unknown")?
    };
    Ok(format!("mtDNA copy number: not assessed ({why})"))
}

/// The listed deletions and how many of them the confirmed list names.
fn listed(dl: &Value) -> Result<(&[Value], usize), Malformed> {
    let rows = items(
        Some(need(dl, "listed", "sections.mito.data.deletions")?),
        "sections.mito.data.deletions.listed",
    )?;
    let mut named = 0;
    for r in rows {
        if get(r, "known", "sections.mito.data.deletions.listed[]")?.is_some_and(truthy) {
            named += 1;
        }
    }
    Ok((rows, named))
}

/// Step 20's confirmed variants, large deletions and copy number, the same
/// in both reports.
pub fn mito_extended_lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let mut lines = Vec::new();
    if let Some(conf) = d.get("confirmed")? {
        let conf = items(Some(conf), "sections.mito.data.confirmed")?;
        lines.push(format!(
            "Confirmed disease variants (assets/mito_confirmed.tsv): {}",
            conf.len()
        ));
        for v in conf {
            lines.push(format!("  {}", mito_confirmed_text(v)?));
        }
    }
    if let Some(n) = d.get("confirmed_note")?.filter(|v| truthy(v)) {
        lines.push(format!(
            "Confirmed disease variants: not checked ({})",
            py(n)?
        ));
    }
    if let Some(dl) = d.get("deletions")?.filter(|v| truthy(v)) {
        let (rows, named) = listed(dl)?;
        let dp = "sections.mito.data.deletions";
        lines.push(format!(
            "Large deletions (eKLIPse): {}, of them {named} named in assets/mito_confirmed.tsv; listed: named \
             ones and any at {}% or more",
            py(need(dl, "total", dp)?)?,
            py(need(dl, "list_floor_pct", dp)?)?
        ));
        for r in rows {
            lines.push(format!("  {}", mito_deletion_text(r)?));
        }
        lines.push(MITO_DELETION_NOTE.to_owned());
    } else if let Some(n) = d.get("deletions_note")?.filter(|v| truthy(v)) {
        lines.push(format!("Large deletions: not looked for ({})", py(n)?));
    }
    if let Some(c) = d.get("copies")?.filter(|v| truthy(v)) {
        lines.push(mito_copies_text(c)?);
        if is_str(get(c, "status", "sections.mito.data.copies")?, "called") {
            lines.push(MITO_COPIES_NOTE.to_owned());
        }
    } else if let Some(n) = d.get("copies_note")?.filter(|v| truthy(v)) {
        lines.push(format!("mtDNA copy number: not counted ({})", py(n)?));
    }
    Ok(lines)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "mito", Some("Mitochondrial Variants (Mutect2)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = s.section("mito");
        w.push(format!("  PASS variants: {}", py(d.need("pass")?)?));
        w.push(format!(
            "  Heteroplasmic (AF 0.05 to 0.95, outside the D-loop poly-C tracts): {}",
            py(d.need("heteroplasmic")?)?
        ));
        w.extend(mito_lines(&d)?.into_iter().map(|l| format!("  {l}")));
        if d.get("polyc_heteroplasmic")?.is_some_and(truthy) {
            w.push(format!("  {MITO_POLYC_NOTE}"));
        }
        w.push(format!("  {MITO_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

fn stat(label: &str, value: Inline) -> Line {
    Line::Stat(label.to_owned(), value)
}

fn count(n: String, any: bool) -> Inline {
    Inline::Badge(n, if any { Colour::Yellow } else { Colour::Gray })
}

/// The HTML report's card: full width when it holds the table of annotated
/// calls.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let mi = s.section("mito");
    let empty = Value::Object(Default::default());
    let ma = annotation(&mi)?.unwrap_or(&empty);
    let ap = "sections.mito.data.annotation";
    let pc = items(
        mi.get("polyc_heteroplasmic")?.filter(|v| truthy(v)),
        "sections.mito.data.polyc_heteroplasmic",
    )?;
    let mut caveats = Vec::new();
    let mut body = vec![
        stat(
            "chrM variants (PASS)",
            Inline::Text(or_none(mi.get("pass")?)?),
        ),
        stat(
            "Heteroplasmic (AF 0.05 to 0.95, outside the D-loop poly-C tracts)",
            Inline::Text(or_none(mi.get("heteroplasmic")?)?),
        ),
        stat(
            "Length heteroplasmy in the D-loop poly-C tracts",
            Inline::Text(pc.len().to_string()),
        ),
    ];
    for x in pc {
        body.push(Line::Note(polyc_text(x)?));
    }
    if !pc.is_empty() {
        caveats.push(MITO_POLYC_NOTE);
    }
    let variants = get(ma, "variants", ap)?.filter(|v| truthy(v));
    if get(ma, "source", ap)?.is_some_and(truthy) {
        let plp = get(ma, "clinvar_plp", ap)?;
        body.push(stat(
            "ClinVar pathogenic or likely pathogenic",
            count(plp.map_or(Ok("0".into()), py)?, plp.is_some_and(truthy)),
        ));
        if let Some(vs) = variants {
            let vp = "sections.mito.data.annotation.variants[]";
            let mut rows = Vec::new();
            let mut at = Vec::new();
            for v in items(Some(vs), ap)? {
                let g = |k: &str| -> Result<String, Malformed> { or_none(get(v, k, vp)?) };
                let pos = need(v, "pos", vp)?;
                rows.push(vec![
                    format!(
                        "chrM:{} {}>{}",
                        py(pos)?,
                        py(need(v, "ref", vp)?)?,
                        py(need(v, "alt", vp)?)?
                    ),
                    g("af")?,
                    g("clinvar_significance")?,
                    g("clinvar_review")?,
                    g("clinvar_condition")?,
                    g("gnomad_af_hom")?,
                    g("gnomad_af_het")?,
                    if get(v, "polyc", vp)?.is_some_and(truthy) {
                        "yes".into()
                    } else {
                        String::new()
                    },
                ]);
                at.push(if truthy(pos) {
                    Some(format!("chrM:{}", py(pos)?))
                } else {
                    None
                });
            }
            body.push(Line::Table(Table {
                head: [
                    "Variant",
                    "AF",
                    "ClinVar",
                    "Review",
                    "Condition",
                    "gnomAD AF hom",
                    "gnomAD AF het",
                    "Poly-C tract",
                ]
                .map(str::to_owned)
                .to_vec(),
                rows,
                viewer: Some(at),
                layout: Layout::OneLine,
            }));
            body.push(Line::ViewerNote);
        }
        body.push(Line::Note(format!(
            "ClinVar {}; gnomAD {}",
            or_text(get(ma, "clinvar", ap)?, "unknown")?,
            or_text(get(ma, "gnomad", ap)?, "unknown")?
        )));
    } else if let Some(n) = get(ma, "note", ap)?.filter(|v| truthy(v)) {
        body.push(Line::Note(format!(
            "ClinVar and gnomAD: not annotated ({})",
            py(n)?
        )));
    }
    if let Some(conf) = mi.get("confirmed")? {
        let conf = items(Some(conf), "sections.mito.data.confirmed")?;
        body.push(stat(
            "Confirmed disease variants",
            count(conf.len().to_string(), !conf.is_empty()),
        ));
        for v in conf {
            body.push(Line::Note(mito_confirmed_text(v)?));
        }
    }
    if let Some(dl) = mi.get("deletions")?.filter(|v| truthy(v)) {
        let (rows, named) = listed(dl)?;
        body.push(stat(
            "Large deletions (eKLIPse)",
            Inline::Text(py(need(dl, "total", "sections.mito.data.deletions")?)?),
        ));
        body.push(stat(
            "Named in the confirmed list",
            count(named.to_string(), named > 0),
        ));
        for r in rows {
            body.push(Line::Note(mito_deletion_text(r)?));
        }
        caveats.push(MITO_DELETION_NOTE);
    } else if let Some(n) = mi.get("deletions_note")?.filter(|v| truthy(v)) {
        body.push(Line::Note(format!(
            "Large deletions: not looked for ({})",
            py(n)?
        )));
    }
    if let Some(mc) = mi.get("copies")?.filter(|v| truthy(v)) {
        let called = is_str(get(mc, "status", "sections.mito.data.copies")?, "called");
        body.push(stat(
            "mtDNA copies per cell",
            Inline::Text(if called {
                py(need(mc, "mtdna_copies", "sections.mito.data.copies")?)?
            } else {
                "not assessed".into()
            }),
        ));
        body.push(Line::Note(mito_copies_text(mc)?));
        if called {
            caveats.push(MITO_COPIES_NOTE);
        }
    } else if let Some(n) = mi.get("copies_note")?.filter(|v| truthy(v)) {
        body.push(Line::Note(format!(
            "mtDNA copy number: not counted ({})",
            py(n)?
        )));
    }
    caveats.push(MITO_NOTE);
    let card = html::card(
        &mi,
        Card {
            title: "Mitochondrial Analysis".into(),
            full: variants.is_some(),
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
