//! Step 46's LPA KIV-2 repeat count with its 95% interval, the KIV-2
//! 4925G>A null allele, step 42's two LPA risk variants, the flag and the
//! advice: the same lines in both reports.

use pycompat::json::Value;

use crate::html;
use crate::model::{
    Card, Colour, HtmlPart, Inline, Line, Malformed, Section, State, Summary, get, is_str, need,
    py, py_eq, truthy,
};
use crate::notes::LPA_NOTE;
use crate::sections::paralogs::shows_body;
use crate::sections::smn1::marker_state;
use crate::text;

/// What a null allele means for the level, and how to measure it.
pub const LPA_NULL_CAVEAT: &str = "A null allele lowers the Lp(a) made from the copies it sits on; \
    the Lp(a) blood test is the measure.";

/// The lines of the text section and of the card.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let dot = |k: &str| d.get(k).and_then(|v| v.map_or(Ok(".".to_owned()), py));
    let need = |k: &str| d.need(k).and_then(py);
    let mut out = vec![if is_str(d.get("status")?, "estimated") {
        format!(
            "KIV-2 repeats (both chromosomes): {} (95% interval {} to {}), {}",
            need("kiv2_copies")?,
            need("ci95_low")?,
            need("ci95_high")?,
            need("class")?
        )
    } else {
        format!("KIV-2 repeats: {}", dot("status")?)
    }];
    out.push(format!("  {}", dot("message")?));
    if let Some(na) = d.get("null_allele")?.filter(|v| truthy(v)) {
        let path = "sections.lpa.data.null_allele";
        let rsid = match get(na, "rsid", path)?.filter(|v| truthy(v)) {
            Some(v) => py(v)?,
            None => ".".to_owned(),
        };
        let rsid = rsid.split(' ').next().unwrap_or_default();
        out.push(format!(
            "KIV-2 4925G>A null allele ({rsid}): {}",
            marker_state(na)?
        ));
        out.push(format!(
            "  {}",
            get(na, "message", path)?.map_or(Ok(".".to_owned()), py)?
        ));
        out.push(format!("  {LPA_NULL_CAVEAT}"));
    }
    let variants: &[Value] = match d.get("variants")? {
        None => &[],
        Some(Value::Array(a)) => a,
        Some(_) => return Err(Malformed("sections.lpa.data.variants is not a list".into())),
    };
    for v in variants {
        let path = "sections.lpa.data.variants[]";
        let n = need_v(v, "risk_alleles", path)?;
        let count = if matches!(n, Value::Null) {
            "not read".to_owned()
        } else {
            let s = if py_eq(n, &Value::Int(1)) { "" } else { "s" };
            format!("{} risk allele{s}", py(n)?)
        };
        out.push(format!(
            "{}: {} ({count})",
            py(need_v(v, "rsid", path)?)?,
            py(need_v(v, "genotype", path)?)?
        ));
    }
    if !d.get("variants")?.is_some_and(truthy) {
        out.push(
            "rs10455872 and rs3798220: not read (step 42, the traits panel, did not run)"
                .to_owned(),
        );
    }
    out.push(dot("flag_text")?);
    out.push(dot("advice")?);
    Ok(out)
}

fn need_v<'a>(v: &'a Value, k: &str, path: &str) -> Result<&'a Value, Malformed> {
    need(v, k, path)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "lpa",
        Some("Lp(a): LPA KIV-2 Repeat Count and Risk Variants"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        w.extend(lines(&s.section("lpa"))?.iter().map(|l| format!("  {l}")));
        w.push(format!("  {LPA_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when step 46 is in the summary.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("lpa");
    if d.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut body = Vec::new();
    let mut caveats = Vec::new();
    if shows_body(&d) {
        if is_str(d.get("flag")?, "high") {
            body.push(Line::Stat(
                "Genetic sign of high Lp(a)".to_owned(),
                Inline::Badge("Yes".to_owned(), Colour::Yellow),
            ));
        }
        body.extend(lines(&d)?.into_iter().map(Line::Para));
        caveats.push(LPA_NOTE);
    }
    let card = html::card(
        &d,
        Card {
            title: "Lp(a): LPA KIV-2 Repeat Count".to_owned(),
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
