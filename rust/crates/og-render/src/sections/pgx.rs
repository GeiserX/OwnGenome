//! The medicines sections: PharmCAT (step 07), the CPIC recommendations
//! (step 27), and pypgx (step 32) with Cyrius (step 21) and step 36's
//! CYP2D6 verdict, which both reports show in one section (`cyp2d6`).

use pycompat::json::Value;
use pycompat::text::{lower, slice, upper};

use crate::model::{
    Card, Cell, Colour, HtmlPart, Inline, Layout, Line, Malformed, Section, State, Summary, Table,
    contains, get, is_str, need, py, truthy,
};
use crate::notes::{CPIC_NOTE, PGX_CALLERS_NOTE, PHARMCAT_NOTE, stale_note};
use crate::{html, text};

/// `f"{v:<width}"` for the values a summary holds there: a text, a number
/// (a bool formats as its int), anything else is malformed, as Python's
/// format would refuse it.
pub fn left(v: &Value, width: usize, path: &str) -> Result<String, Malformed> {
    let s = match v {
        Value::Str(_) | Value::Int(_) | Value::Float(_) => py(v)?,
        Value::Bool(b) => u8::from(*b).to_string(),
        _ => return Err(Malformed(format!("{path} cannot be formatted"))),
    };
    Ok(pycompat::fmt::ljust(&s, width))
}

/// `v` as a text, or the path of a value that is not one (where Python
/// joins or concatenates strings).
pub fn text_of(v: &Value, path: &str) -> Result<String, Malformed> {
    match v {
        Value::Str(s) => Ok(s.clone()),
        _ => Err(Malformed(format!("{path} is not a text"))),
    }
}

/// `v or default`, as text.
pub fn or_text(v: Option<&Value>, default: &str) -> Result<String, Malformed> {
    match v.filter(|v| truthy(v)) {
        Some(v) => py(v),
        None => Ok(default.to_owned()),
    }
}

/// `d.get(key)` as Python's f-string prints it: `None` for a missing key.
pub fn or_none(v: Option<&Value>) -> Result<String, Malformed> {
    v.map_or(Ok("None".to_owned()), py)
}

/// The items of a list value (`[]` for a missing key or `null`).
pub fn items<'a>(v: Option<&'a Value>, path: &str) -> Result<&'a [Value], Malformed> {
    match v {
        None | Some(Value::Null) => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!("{path} is not a list"))),
    }
}

/// `str.capitalize()`: the first character upper case, the rest lower.
fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        Some(f) => upper(&f.to_string()) + &lower(c.as_str()),
        None => String::new(),
    }
}

/// `os.path.join(os.path.dirname(path), name)`, on `/` paths.
fn beside(path: &str, name: &str) -> String {
    let head = path.rfind('/').map_or("", |i| &path[..=i]);
    let dir = if !head.is_empty() && head.chars().any(|c| c != '/') {
        head.trim_end_matches('/')
    } else {
        head
    };
    if name.starts_with('/') {
        name.to_owned()
    } else if dir.is_empty() || dir.ends_with('/') {
        format!("{dir}{name}")
    } else {
        format!("{dir}/{name}")
    }
}

fn stat(label: &str, value: Inline) -> Line {
    Line::Stat(label.to_owned(), value)
}

fn sample(s: &Summary) -> Result<String, Malformed> {
    text_of(s.need("sample")?, "summary.sample")
}

/// `<sample>_cpic_recommendations.txt` beside the CPIC report.
fn cpic_file(s: &Summary, sec: &Section<'_>) -> Result<String, Malformed> {
    let src = sec
        .source
        .ok_or_else(|| Malformed("sections.cpic has no source".into()))?;
    let src = text_of(src, "sections.cpic.source")?;
    Ok(beside(
        &src,
        &format!("{}_cpic_recommendations.txt", sample(s)?),
    ))
}

/// The step ran and the summary names its file.
fn ran_with_source(sec: &Section<'_>) -> bool {
    sec.state.ran() && sec.source.is_some_and(truthy)
}

// --- PharmCAT ----------------------------------------------------------------------

/// The text report's PharmCAT block.
pub fn pharmcat_text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "pharmcat", Some("Pharmacogenomics (PharmCAT)"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let sec = s.section("pharmcat");
        w.push(format!(
            "  Report: {} (PharmCAT {})",
            or_none(sec.source)?,
            or_text(sec.get("version")?, "version unknown")?
        ));
        if sec.get("html_report")?.is_some_and(truthy) {
            w.push("  HTML report beside it (open in a browser for full details)".to_owned());
        }
        w.push(format!("  {PHARMCAT_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's PharmCAT card, and the link to PharmCAT's own HTML
/// report beside it.
pub fn pharmcat_html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("pharmcat");
    // Never green: a report that was written is not a finding.
    let mut body = vec![stat(
        "PharmCAT report",
        if sec.state.ran() {
            Inline::Badge("Complete".into(), Colour::State)
        } else {
            Inline::Badge("Not run".into(), Colour::Gray)
        },
    )];
    if let Some(v) = sec.get("version")?.filter(|v| truthy(v)) {
        body.push(stat(
            "PharmCAT version",
            Inline::Text(text_of(v, "sections.pharmcat.data.version")?),
        ));
    }
    let pg = s.section("pypgx");
    if pg.state.ran() && contains(pg.data, "comparison") {
        let c = need(
            pg.need("comparison")?,
            "conflicts",
            "sections.pypgx.data.comparison",
        )?;
        body.push(stat("PharmCAT vs pypgx conflicts", Inline::Text(py(c)?)));
    }
    let card = html::card(
        &sec,
        Card {
            title: "Pharmacogenomics".into(),
            caveats: vec![PHARMCAT_NOTE],
            body,
            ..Card::default()
        },
    )?;
    let mut files = Vec::new();
    if ran_with_source(&sec) && sec.get("html_report")?.is_some_and(truthy) {
        let src = text_of(sec.source.expect("checked"), "sections.pharmcat.source")?;
        files.push((
            format!("{}.html", slice(&src, None, Some(-5))),
            "PharmCAT HTML report".to_owned(),
        ));
    }
    Ok(HtmlPart {
        cards: vec![card],
        files,
    })
}

// --- CPIC --------------------------------------------------------------------------

/// The genes CPIC lists as non-normal or with more than one possible
/// result: (gene, diplotype, phenotype).
fn unusual<'a>(sec: &Section<'a>) -> Result<Vec<[&'a Value; 3]>, Malformed> {
    let mut out = Vec::new();
    for (i, g) in items(Some(sec.need("genes")?), "sections.cpic.data.genes")?
        .iter()
        .enumerate()
    {
        let p = format!("sections.cpic.data.genes[{i}]");
        let st = need(g, "status", &p)?;
        if is_str(Some(st), "non-normal") || is_str(Some(st), "ambiguous") {
            out.push([
                need(g, "gene", &p)?,
                need(g, "diplotype", &p)?,
                need(g, "phenotype", &p)?,
            ]);
        }
    }
    Ok(out)
}

/// The text report's CPIC block.
pub fn cpic_text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "cpic", Some("CPIC Drug-Gene Recommendations"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let sec = s.section("cpic");
        if truthy(sec.need("parse_failed")?) {
            w.push("  The PharmCAT report could not be parsed: see the CPIC report.".to_owned());
        }
        let or0 = |k: &str| -> Result<String, Malformed> { sec.get(k)?.map_or(Ok("0".into()), py) };
        w.push(format!(
            "  Genes with a non-normal phenotype: {}  More than one possible result: {}  Not called: {}  \
             Called without a function phenotype (no drug guidance): {}",
            py(sec.need("non_normal")?)?,
            or0("ambiguous")?,
            py(sec.need("not_called")?)?,
            or0("unclassified")?
        ));
        for [g, d, p] in unusual(&sec)? {
            w.push(format!(
                "    {} {} {}",
                left(g, 10, "a CPIC gene")?,
                left(d, 28, "a CPIC diplotype")?,
                py(p)?
            ));
        }
        for warn in items(Some(sec.need("warnings")?), "sections.cpic.data.warnings")? {
            w.push(format!("  WARNING: {}", py(warn)?));
        }
        w.push(format!("  Drugs per gene: {}", cpic_file(s, &sec)?));
        w.push(format!("  {CPIC_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's CPIC cards: the counts, then the genes with a
/// non-normal or unresolved phenotype; and the link to the recommendations
/// file.
pub fn cpic_html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("cpic");
    let smp = sample(s)?;
    let mut body = vec![stat(
        "Status",
        if sec.state == State::Stale {
            Inline::Badge("Stale".into(), Colour::Yellow)
        } else {
            Inline::Badge("Complete".into(), Colour::Green)
        },
    )];
    if sec.get("parse_failed")?.is_some_and(truthy) {
        body.push(stat(
            "PharmCAT report",
            Inline::Badge("could not be parsed".into(), Colour::Red),
        ));
    } else {
        for (label, k) in [
            ("Genes with a non-normal phenotype", "non_normal"),
            ("Genes with more than one possible result", "ambiguous"),
            ("Genes not called", "not_called"),
            (
                "Genes called without a function phenotype (no drug guidance)",
                "unclassified",
            ),
        ] {
            let v = sec.get(k)?.map_or(Ok("0".into()), py)?;
            body.push(stat(label, Inline::Text(v)));
        }
    }
    body.push(stat(
        "Tip",
        Inline::Quiet(format!(
            "Recommendations per gene: cpic/{smp}_cpic_recommendations.txt"
        )),
    ));
    let mut cards = vec![html::card(
        &sec,
        Card {
            title: "CPIC Drug Recommendations".into(),
            caveats: vec![CPIC_NOTE],
            body,
            ..Card::default()
        },
    )?];
    if sec.state.ran()
        && (truthy(sec.need("non_normal")?) || sec.get("ambiguous")?.is_some_and(truthy))
    {
        let mut rows = Vec::new();
        for [g, d, p] in unusual(&sec)? {
            rows.push(vec![
                Cell::Text(text_of(g, "a CPIC gene")?),
                Cell::Text(text_of(d, "a CPIC diplotype")?),
                Cell::Text(text_of(p, "a CPIC phenotype")?),
            ]);
        }
        let mut body = vec![Line::Table(Table {
            head: vec!["Gene".into(), "Diplotype".into(), "Phenotype".into()],
            rows,
            viewer: false,
            layout: Layout::RowLines,
        })];
        for x in items(Some(sec.need("warnings")?), "sections.cpic.data.warnings")? {
            body.push(Line::Warn(format!(
                "WARNING: {}",
                text_of(x, "a CPIC warning")?
            )));
        }
        body.push(Line::Para(format!(
            "Drugs for each gene: cpic/{smp}_cpic_recommendations.txt"
        )));
        cards.push(html::card(
            &sec,
            Card {
                title: "Genes With a Non-Normal or Unresolved Phenotype (CPIC)".into(),
                full: true,
                body,
                ..Card::default()
            },
        )?);
    }
    let mut files = Vec::new();
    if ran_with_source(&sec) {
        files.push((
            cpic_file(s, &sec)?,
            "CPIC recommendations per gene (text)".to_owned(),
        ));
    }
    Ok(HtmlPart { cards, files })
}

// --- pypgx, Cyrius and CYP2D6 ------------------------------------------------------

/// Why step 32 could not call a gene on this reference.
fn not_callable_why(g: &Value) -> &'static str {
    if is_str(Some(g), "GSTT1") {
        "GRCh38 has it only on the ALT contig chr22_KI270879v1_alt, which a reference without ALT \
         contigs lacks"
    } else {
        "see docs/32-pypgx.md"
    }
}

/// Step 32's lines below its counts, the same in both reports: the genes
/// this reference cannot hold, the calls resting on a low-GQ variant, and
/// the CYP2D6 depth check with its MAPQ shares.
pub fn pypgx_lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let mut lines = Vec::new();
    for g in items(d.get("not_callable")?, "sections.pypgx.data.not_callable")? {
        lines.push(format!(
            "Not assessable on this reference: {} ({})",
            py(g)?,
            not_callable_why(g)
        ));
    }
    for (i, x) in items(d.get("low_gq")?, "sections.pypgx.data.low_gq")?
        .iter()
        .enumerate()
    {
        let p = format!("sections.pypgx.data.low_gq[{i}]");
        lines.push(format!(
            "Low-confidence: {} {} rests on {}: check the reads",
            py(need(x, "gene", &p)?)?,
            py(need(x, "diplotype", &p)?)?,
            or_text(Some(need(x, "weakest", &p)?), "a low-GQ variant")?
        ));
    }
    if let Some(dc) = d.get("depth_check")?.filter(|v| truthy(v)) {
        let p = "sections.pypgx.data.depth_check";
        let share = |k: &str| -> Result<String, Malformed> {
            get(dc, k, p)?.map_or(Ok("not measured".into()), py)
        };
        lines.push(format!(
            "CYP2D6 depth check: {}; share of the depth at MAPQ >= 1 / >= 20: Cyrius's region with REP6 {} / {}, \
             gene body without REP6 {} / {}",
            or_text(get(dc, "status", p)?, "unreadable")?,
            share("gene_mapq1_share")?,
            share("gene_mapq20_share")?,
            share("body_mapq1_share")?,
            share("body_mapq20_share")?
        ));
    }
    Ok(lines)
}

/// (gene, copy number) rows, the same in both reports: pypgx's call from
/// read depth for each BAM-mode gene of assets/pgx_genes.tsv.
pub fn pgx_copy_number(s: &Summary) -> Result<Vec<(String, String)>, Malformed> {
    let pg = s.section("pypgx");
    match pg.get("cnv_call")?.filter(|v| truthy(v)) {
        None => Ok(Vec::new()),
        Some(Value::Object(o)) => o
            .iter()
            .map(|(g, v)| Ok((g.clone(), format!("{} (pypgx, read depth)", py(v)?))))
            .collect(),
        Some(_) => Err(Malformed(
            "sections.pypgx.data.cnv_call is not an object".into(),
        )),
    }
}

/// The summary's CYP2D6 block (calls, agreement, step 36's consensus, the
/// probable call), or `None` for a summary without one (a part collected
/// with `--only`), which shows no CYP2D6 across the callers.
fn block(s: &Summary) -> Result<Option<&Value>, Malformed> {
    match s.get("cyp2d6") {
        None | Some(Value::Null) => Ok(None),
        Some(c @ Value::Object(_)) => Ok(Some(c)),
        Some(_) => Err(Malformed("summary.cyp2d6 is not an object".into())),
    }
}

/// The callers' CYP2D6 calls, in the summary's order.
fn calls(c: &Value) -> Result<Vec<(&String, &Value)>, Malformed> {
    match need(c, "calls", "summary.cyp2d6")? {
        Value::Object(o) => Ok(o.iter().collect()),
        _ => Err(Malformed("summary.cyp2d6.calls is not an object".into())),
    }
}

/// The column label of a caller: PharmCAT calls no CYP2D6 from a VCF, it
/// shows the call step 36 passed to it.
fn caller_label(k: &str) -> &str {
    if k == "PharmCAT" {
        "PharmCAT (step 36's call)"
    } else {
        k
    }
}

/// Step 36's CYP2D6 verdict as one line, or `""` without one.
pub fn cyp2d6_verdict(c: &Value) -> Result<String, Malformed> {
    let Some(v) = get(c, "consensus", "summary.cyp2d6")?.filter(|v| truthy(v)) else {
        return Ok(String::new());
    };
    let p = "summary.cyp2d6.consensus";
    let passed = get(v, "passed_to_pharmcat", p)?.is_some_and(truthy);
    let where_ = if passed {
        "passed to PharmCAT"
    } else {
        "not passed to PharmCAT"
    };
    Ok(format!(
        "{}, {where_} ({})",
        or_text(get(v, "result", p)?, "?")?,
        or_text(get(v, "reason", p)?, "no reason given")?
    ))
}

/// The probable CYP2D6 call as one line.
fn probable_text(p: &Value) -> Result<String, Malformed> {
    let path = "summary.cyp2d6.probable";
    let pheno = match get(p, "phenotype", path)?.filter(|v| truthy(v)) {
        Some(x) => format!(", {}", py(x)?),
        None => String::new(),
    };
    Ok(format!(
        "{}{pheno}; not passed to PharmCAT, so no drug guidance rests on it",
        or_none(get(p, "diplotype", path)?)?
    ))
}

/// The probable call's label and line, when the summary has one.
fn probable(c: &Value) -> Result<Option<(String, String)>, Malformed> {
    let Some(p) = get(c, "probable", "summary.cyp2d6")?.filter(|v| truthy(v)) else {
        return Ok(None);
    };
    let label = text_of(
        need(p, "label", "summary.cyp2d6.probable")?,
        "summary.cyp2d6.probable.label",
    )?;
    Ok(Some((capitalize(&label), probable_text(p)?)))
}

/// The text report's block of the section: pypgx, its copy numbers, then
/// CYP2D6 across the callers.
pub fn cyp2d6_text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let mut w = Vec::new();
    if let Some(head) = text::head(s, "pypgx", Some("pypgx Pharmacogenomics"))? {
        w.extend(head.lines);
        if head.ran {
            let d = s.section("pypgx");
            w.push(format!(
                "  Genes called: {}/{}",
                py(d.need("genes_called")?)?,
                py(d.need("genes_total")?)?
            ));
            w.push(format!("  CYP2D6 diplotype: {}", py(d.need("cyp2d6")?)?));
            if contains(d.data, "comparison") {
                let c = d.need("comparison")?;
                let p = "sections.pypgx.data.comparison";
                w.push(format!(
                    "  PharmCAT vs pypgx conflicts: {}",
                    py(need(c, "conflicts", p)?)?
                ));
                let one = need(c, "one_tool_only", p)?;
                if truthy(one) {
                    w.push(format!(
                        "  Scope differences (only one tool called): {}",
                        py(one)?
                    ));
                }
            }
            w.extend(pypgx_lines(&d)?.into_iter().map(|l| format!("  {l}")));
            w.push(String::new());
        }
    }
    let cn = pgx_copy_number(s)?;
    if !cn.is_empty() {
        w.push("## PGx copy number".to_owned());
        w.push("---".to_owned());
        for (g, v) in cn {
            w.push(format!("  {} {v}", pycompat::fmt::ljust(&g, 8)));
        }
        w.push(String::new());
    }
    let Some(c) = block(s)? else {
        return Ok(w);
    };
    let calls = calls(c)?;
    let verdict = cyp2d6_verdict(c)?;
    if calls.iter().any(|(_, v)| truthy(v)) || !verdict.is_empty() {
        w.push("## CYP2D6 across callers".to_owned());
        w.push("---".to_owned());
        // The consensus change is filed under pypgx.
        if let Some(line) = text::new_since(s, &["pypgx", "cyrius"])? {
            w.push(line);
        }
        for (k, v) in &calls {
            w.push(format!(
                "  {} {}",
                pycompat::fmt::ljust(caller_label(k), 26),
                or_text(Some(v), "not run")?
            ));
        }
        let agree = need(c, "agree", "summary.cyp2d6")?;
        w.push(format!(
            "  pypgx and Cyrius agree: {}",
            if truthy(agree) {
                "yes"
            } else if *agree == Value::Bool(false) {
                "NO, review before acting on CYP2D6"
            } else {
                "fewer than two usable calls"
            }
        ));
        if !verdict.is_empty() {
            w.push(format!("  Step 36: {verdict}"));
        }
        if let Some((label, line)) = probable(c)? {
            w.push(format!("  {label}: {line}"));
        }
        w.push(format!("  {PGX_CALLERS_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's cards of the section: pypgx (when its step ran or
/// failed), its copy numbers, and CYP2D6 across the callers, which is
/// always there.
pub fn cyp2d6_html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let mut cards = Vec::new();
    let pd = s.section("pypgx");
    if pd.state != State::Missing {
        let mut body = vec![
            stat(
                "Genes called",
                Inline::Text(format!(
                    "{}/{}",
                    or_none(pd.get("genes_called")?)?,
                    or_none(pd.get("genes_total")?)?
                )),
            ),
            stat(
                "CYP2D6 diplotype",
                Inline::Text(or_text(pd.get("cyp2d6")?, "not called")?),
            ),
        ];
        body.extend(pypgx_lines(&pd)?.into_iter().map(Line::Para));
        cards.push(html::card(
            &pd,
            Card {
                title: "pypgx Pharmacogenomics".into(),
                body,
                ..Card::default()
            },
        )?);
    }
    let cn = pgx_copy_number(s)?;
    if !cn.is_empty() {
        cards.push(Card {
            title: "PGx Copy Number".into(),
            body: cn
                .into_iter()
                .map(|(g, v)| stat(&g, Inline::Text(v)))
                .collect(),
            ..Card::default()
        });
    }
    let Some(c) = block(s)? else {
        return Ok(HtmlPart {
            cards,
            files: Vec::new(),
        });
    };
    let mut body: Vec<Line> = Vec::new();
    for (k, v) in calls(c)? {
        body.push(stat(
            caller_label(k),
            Inline::Text(or_text(Some(v), "not run")?),
        ));
    }
    let agree = need(c, "agree", "summary.cyp2d6")?;
    body.push(stat(
        "pypgx and Cyrius agree",
        if truthy(agree) {
            Inline::Badge("yes".into(), Colour::Green)
        } else if *agree == Value::Bool(false) {
            Inline::Badge("no".into(), Colour::Red)
        } else {
            Inline::Badge("fewer than two usable calls".into(), Colour::Gray)
        },
    ));
    if let Some(cons) = get(c, "consensus", "summary.cyp2d6")?.filter(|v| truthy(v)) {
        let p = "summary.cyp2d6.consensus";
        let passed = get(cons, "passed_to_pharmcat", p)?.is_some_and(truthy);
        body.push(stat(
            "Step 36",
            Inline::Badge(
                or_text(get(cons, "result", p)?, "?")?,
                if passed {
                    Colour::State
                } else {
                    Colour::Yellow
                },
            ),
        ));
        body.push(Line::Note(format!(
            "{}{}",
            if passed {
                "Passed to PharmCAT: "
            } else {
                "Not passed to PharmCAT: "
            },
            or_text(get(cons, "reason", p)?, "no reason given")?
        )));
    }
    if let Some((label, line)) = probable(c)? {
        body.push(stat(&label, Inline::Text(line)));
    }
    let mut stale = Vec::new();
    for k in ["cpic", "pypgx", "cyrius"] {
        let sec = s.section(k);
        if matches!(sec.state, State::Stale | State::Failed) {
            stale.push(format!(
                "{}: {}",
                text_of(sec.title, &format!("sections.{k}.title"))?,
                stale_note(&sec)?
            ));
        }
    }
    cards.push(Card {
        title: "CYP2D6 Across Callers".into(),
        caveats: vec![PGX_CALLERS_NOTE],
        stale,
        body,
        ..Card::default()
    });
    Ok(HtmlPart {
        cards,
        files: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paths_join_as_os_path_does() {
        assert_eq!(beside("cpic/S_cpic.json", "S_x.txt"), "cpic/S_x.txt");
        assert_eq!(beside("S_cpic.json", "S_x.txt"), "S_x.txt");
        assert_eq!(beside("/r.json", "x"), "/x");
        assert_eq!(beside("a//b", "x"), "a/x");
        assert_eq!(beside("a/b", "/abs"), "/abs");
    }

    #[test]
    fn capitalize_lowers_the_rest() {
        assert_eq!(capitalize("probable CYP2D6 call"), "Probable cyp2d6 call");
        assert_eq!(capitalize(""), "");
    }

    #[test]
    fn left_pads_as_format_does() {
        assert_eq!(left(&Value::Str("ab".into()), 4, "x").unwrap(), "ab  ");
        assert_eq!(left(&Value::Bool(true), 3, "x").unwrap(), "1  ");
        assert!(left(&Value::Null, 3, "x").is_err());
    }
}
