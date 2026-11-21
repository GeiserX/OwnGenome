//! Step 40's raw-data export: how many chip sites got a genotype, why the
//! others read `--`, the two files and what to weigh before uploading one.

use pycompat::fmt::fixed;
use pycompat::json::Value;

use crate::model::{Card, HtmlPart, Line, Malformed, Section, State, Summary, num, py, truthy};
use crate::notes::CONSUMER_NOTE;
use crate::{html, text};

const KEY: &str = "consumer_export";
const TITLE: &str = "Raw-Data Export (23andMe v5 layout)";

/// `d.get(key, default)`, and its text.
fn or<'a>(d: &Section<'a>, key: &str, default: &'a Value) -> Result<&'a Value, Malformed> {
    Ok(d.get(key)?.unwrap_or(default))
}

/// Step 40's lines, the same in both reports, without the closing caveat.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    let zero = Value::Int(0);
    let dot = Value::Str(".".to_owned());
    let sites = or(d, "sites", &zero)?;
    let called = or(d, "called", &zero)?;
    let share = if truthy(sites) {
        let (Some(n), Some(of)) = (num(Some(called)), num(Some(sites))) else {
            return Err(Malformed(format!(
                "sections.{KEY}.data: called and sites are not numbers"
            )));
        };
        format!(" ({}%)", fixed(100.0 * n / of, 1))
    } else {
        String::new()
    };
    let mut out = vec![
        format!(
            "Sites with a genotype: {} of {}{share}; '--' at {}",
            py(called)?,
            py(sites)?,
            py(or(d, "no_call", &zero)?)?
        ),
        format!(
            "Homozygous reference from the gVCF's reference blocks: {}",
            py(or(d, "from_reference_blocks", &zero)?)?
        ),
        format!(
            "Files: consumer_export/{} and consumer_export/{}",
            py(or(d, "file_23andme", &dot)?)?,
            py(or(d, "file_ancestrydna", &dot)?)?
        ),
    ];
    if let Some(reasons) = d.get("no_call_reasons")?.filter(|v| truthy(v)) {
        let Value::Object(m) = reasons else {
            return Err(Malformed(format!(
                "sections.{KEY}.data.no_call_reasons is not a dict"
            )));
        };
        let mut pairs = Vec::new();
        for (k, v) in m {
            let n = num(Some(v)).ok_or_else(|| {
                Malformed(format!(
                    "sections.{KEY}.data.no_call_reasons.{k} is not a number"
                ))
            })?;
            pairs.push((k.as_str(), n, v));
        }
        // Most first; a tie keeps the summary's order, as Python's stable sort does.
        pairs.sort_by(|a, b| b.1.total_cmp(&a.1));
        let mut parts = Vec::new();
        for (k, _, v) in pairs {
            parts.push(format!("{k} {}", py(v)?));
        }
        out.push(format!("'--' because: {}", parts.join("; ")));
    }
    Ok(out)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, KEY, Some(TITLE))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = s.section(KEY);
        w.extend(lines(&d)?.into_iter().map(|l| format!("  {l}")));
        w.push(format!("  {CONSUMER_NOTE}"));
        w.push(String::new());
    }
    Ok(w)
}

/// The HTML report's card, when the step ran or was meant to.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section(KEY);
    if d.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let card = html::card(
        &d,
        Card {
            title: TITLE.to_owned(),
            caveats: vec![CONSUMER_NOTE],
            body: lines(&d)?.into_iter().map(Line::Para).collect(),
            ..Card::default()
        },
    )?;
    Ok(HtmlPart {
        cards: vec![card],
        files: Vec::new(),
    })
}
