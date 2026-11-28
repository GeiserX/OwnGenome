//! Polygenic risk scores (step 25, pgsc_calc): each score with the variants
//! it matched, its percentile among the most similar reference group when
//! the ancestry panel is installed, the effect a published evaluation
//! measured in that group, and the line that says what the numbers can be
//! compared with: the ancestry caveat, the chrX rows left out and the
//! sex-specific scores when no sex was declared.
//!
//! The scores are shown as the summary lists them.

use pycompat::json::Value;

use crate::model::{
    Card, Cell, HtmlPart, Layout, Line, Malformed, Summary, Table, get, is_str, need, py, truthy,
};
use crate::notes::PRS_ANCESTRY_NOTE;
use crate::sections::ancestry::{Align, items, need_align, need_text, py_or};
use crate::{html, text};

const PATH: &str = "sections.prs.data";
const ROW: &str = "sections.prs.data.scores[]";

/// What a score of the other sex says in place of a percentile.
pub const NOT_APPLICABLE: &str = "not applicable to the declared sex";

/// Part of the line under the scores, so not a caveat of its own.
pub const EFFECT_TEXT: &str = "Effect: the odds ratio (OR) or hazard ratio (HR) a published \
    evaluation of the score measured in people of the same reference group, in that \
    evaluation's unit (its PGS Catalog PPM page says which). It describes the score in that \
    group, not your own risk, and is shown only for a group the score was evaluated in.";

/// The start of the line under the scores when every percentile was
/// withheld for an input that lists variants only (og summary collect's
/// `percentile_note`): the panel ran, so the no-panel line would be false.
pub const WITHHELD_TEXT: &str = "Percentile withheld: the input lists variants only, so each \
    scoring site it leaves out adds nothing to the sum, and the raw sum cannot be compared with \
    the reference panel's samples, whose genotypes are complete. Genotypes with a reference call \
    (0/0) at one genotyped position in five or more keep the percentile; a gVCF from step 3 or a \
    chip VCF usually has them.";

/// What a score, matched or total cell says when the summary has none: og
/// summary collect writes null for a short row of pgsc_calc's table.
pub const NOT_AVAILABLE: &str = "not available";

/// `str(r[key])`, or [`NOT_AVAILABLE`] for null.
fn count(r: &Value, key: &str) -> Result<String, Malformed> {
    match need(r, key, ROW)? {
        Value::Null => Ok(NOT_AVAILABLE.to_owned()),
        v => py(v),
    }
}

fn field<'a>(r: &'a Value, key: &str) -> Result<Option<&'a Value>, Malformed> {
    get(r, key, ROW)
}

/// `r.get(key)` is true.
fn on<'a>(r: &'a Value, key: &str) -> Result<Option<&'a Value>, Malformed> {
    Ok(field(r, key)?.filter(|v| truthy(v)))
}

fn other_sex(r: &Value) -> Result<bool, Malformed> {
    Ok(is_str(field(r, "sex_check")?, NOT_APPLICABLE))
}

/// The scores of the summary.
fn scores(d: &Value) -> Result<Vec<Value>, Malformed> {
    items(need(d, "scores", PATH)?, &format!("{PATH}.scores"))
}

/// The percentile cell: the percentile with its group, the note for a score
/// of the other sex, or "" when there is none.
pub fn percentile(r: &Value) -> Result<String, Malformed> {
    if other_sex(r)? {
        let sex = py_or(on(r, "applies_to")?, "sex")?;
        return Ok(format!("{NOT_APPLICABLE} ({sex}-specific score)"));
    }
    match on(r, "percentile")? {
        Some(p) => Ok(format!(
            "{} ({})",
            py(p)?,
            py_or(on(r, "group")?, "group unknown")?
        )),
        None => Ok(String::new()),
    }
}

/// `", ".join(v)`: the items of a list of texts, or a text's characters.
fn join(v: &Value, at: &str) -> Result<String, Malformed> {
    let parts = items(v, at)?
        .into_iter()
        .map(|x| match x {
            Value::Str(s) => Ok(s),
            _ => Err(Malformed(format!("{at} holds an item that is not a text"))),
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(parts.join(", "))
}

/// The effect cell: the odds or hazard ratio the PGS Catalog records for
/// the score in the sample's group, a line saying the score has none for
/// that group, or "" when the sample has no group, the score does not apply
/// to the declared sex, or the table names no evaluations. The bool says
/// whether the cell is the ratio itself.
pub fn effect(r: &Value) -> Result<(String, bool), Malformed> {
    let evaluated = match field(r, "evaluated_in")? {
        None | Some(Value::Null) => return Ok((String::new(), false)),
        Some(v) => v,
    };
    if on(r, "percentile")?.is_some() && on(r, "group")?.is_some() && !other_sex(r)? {
        if let Some(x) = on(r, "effect")? {
            return Ok((py(x)?, true));
        }
        let group = py(need(r, "group", ROW)?)?;
        let place = if truthy(evaluated) {
            format!(
                "; evaluated in {}",
                join(evaluated, &format!("{ROW}.evaluated_in"))?
            )
        } else {
            String::new()
        };
        return Ok((format!("no published OR or HR for {group}{place}"), false));
    }
    Ok((String::new(), false))
}

/// The chrX rows left out of the scores because no sex was given.
fn chrx_note(rs: &[Value]) -> Result<String, Malformed> {
    let mut left = Vec::new();
    for r in rs {
        let Some(c) = on(r, "chrx")? else { continue };
        let Value::Str(c) = c else {
            return Err(Malformed(format!("{ROW}.chrx is not a text")));
        };
        if c.ends_with("left out") {
            let first = pycompat::text::split(c.as_str(), None, None)
                .first()
                .copied()
                .unwrap_or_default()
                .to_owned();
            left.push(format!("{} ({first})", need_text(r, "pgs_id", ROW)?));
        }
    }
    if left.is_empty() {
        return Ok(String::new());
    }
    Ok(format!(
        " Rows on chrX were left out, because the sample's sex was not given: {}. Give the sex \
         to score them (docs/25-prs.md).",
        left.join(", ")
    ))
}

/// The sex-specific scores shown for any sex, when no sex was declared.
fn sex_note(rs: &[Value]) -> Result<String, Malformed> {
    let mut named = Vec::new();
    for r in rs {
        if is_str(field(r, "sex_check")?, "no sex declared")
            && let Some(a) = on(r, "applies_to")?
        {
            named.push(format!("{} ({})", need_text(r, "pgs_id", ROW)?, py(a)?));
        }
    }
    if named.is_empty() {
        return Ok(String::new());
    }
    Ok(format!(
        " No sex was declared, so these sex-specific scores are shown for any sex: {}. Give the \
         sex to mark a score of the other sex not applicable (docs/25-prs.md).",
        named.join(", ")
    ))
}

/// Whether a score's percentile was withheld (its `percentile_note`).
fn withheld(rs: &[Value]) -> Result<bool, Malformed> {
    for r in rs {
        if on(r, "percentile_note")?.is_some() {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether the summary's table names evaluations (its Evaluated_In column).
fn names_evaluations(rs: &[Value]) -> Result<bool, Malformed> {
    for r in rs {
        if !matches!(field(r, "evaluated_in")?, None | Some(Value::Null)) {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The line under the scores, in both reports: what the numbers can be
/// compared with.
pub fn note(d: &Value) -> Result<String, Malformed> {
    let rs = scores(d)?;
    let mut tail = chrx_note(&rs)? + &sex_note(&rs)?;
    if let Some(skipped) = get(d, "skipped", PATH)? {
        const SKIPPED: &str = "sections.prs.data.skipped";
        for row in items(skipped, SKIPPED)? {
            let condition = need_text(&row, "condition", SKIPPED)?;
            let id = need_text(&row, "pgs_id", SKIPPED)?;
            let reason = need_text(&row, "reason", SKIPPED)?;
            tail.push_str(&format!(" {condition} ({id}) skipped: {reason}."));
        }
    }
    if !get(d, "adjusted", PATH)?.is_some_and(truthy) {
        if withheld(&rs)? {
            return Ok(format!(
                "{WITHHELD_TEXT} {PRS_ANCESTRY_NOTE} See docs/25-prs.md.{tail}"
            ));
        }
        return Ok(format!(
            "Raw score only: no ancestry reference panel was installed, so there is no percentile \
             and the sum cannot be compared with anyone (scripts/setup.sh --ancestry-panel adds \
             it). {PRS_ANCESTRY_NOTE} See docs/25-prs.md.{tail}"
        ));
    }
    let empty = Value::Object(indexmap::IndexMap::new());
    let anc = get(d, "ancestry", PATH)?
        .filter(|v| truthy(v))
        .unwrap_or(&empty);
    let at = format!("{PATH}.ancestry");
    let mut group = match get(anc, "population", &at)?.filter(|v| truthy(v)) {
        Some(p) => py(p)?,
        None => {
            let mut g = String::new();
            for r in &rs {
                if let Some(x) = on(r, "group")? {
                    g = py(x)?;
                    break;
                }
            }
            g
        }
    };
    if group.is_empty() {
        group = "most similar".to_owned();
    }
    let panel = match get(anc, "panel", &at)?.filter(|v| truthy(v)) {
        Some(p) => format!(" of the {} reference panel", py(p)?),
        None => " of the reference panel".to_owned(),
    };
    // A sample between two groups (`mixed`, from og summary collect) is
    // ranked among the first only: name both, as step 25's note does.
    let second = match get(anc, "second_population", &at)?.filter(|v| truthy(v)) {
        Some(p) if is_str(get(anc, "mixed", &at)?, "True") => Some(py(p)?),
        _ => None,
    };
    let low = match second {
        Some(second) => format!(
            " The ancestry match is low-confidence: the genetic ancestry sits between {group} \
             and {second}, and the percentile ranks it among {group} samples only. Read the \
             percentile with care."
        ),
        None if is_str(get(anc, "low_confidence", &at)?, "True") => {
            " The ancestry match is low-confidence: read the percentile with care.".to_owned()
        }
        None => String::new(),
    };
    let eff = if names_evaluations(&rs)? {
        format!(" {EFFECT_TEXT}")
    } else {
        String::new()
    };
    Ok(format!(
        "Percentile: where the score falls among the {group} samples{panel}, the group whose \
         genetic ancestry is most similar to this sample's (pgsc_calc).{low} A percentile is not \
         a risk.{eff} {PRS_ANCESTRY_NOTE} See docs/25-prs.md.{tail}"
    ))
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "prs", Some("Polygenic Risk Scores"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = s.section("prs").data;
    for r in scores(d)? {
        let condition = need_align(&r, "condition", 35, Align::Left, ROW)?;
        let n = |k: &str| count(&r, k);
        let line = format!(
            "  {condition} {} ({}/{} variants matched) {}",
            n("score")?,
            n("matched")?,
            n("total")?,
            need_text(&r, "pgs_id", ROW)?
        );
        let pct = percentile(&r)?;
        let pct = if other_sex(&r)? {
            format!("  {pct}")
        } else if pct.is_empty() {
            String::new()
        } else {
            format!("  percentile {pct}")
        };
        let (e, is_ratio) = effect(&r)?;
        let eff = if e.is_empty() {
            String::new()
        } else if is_ratio {
            format!("  effect {e}")
        } else {
            format!("  {e}")
        };
        w.push(format!("{line}{pct}{eff}"));
    }
    w.push(format!("  {}", note(d)?));
    w.push(String::new());
    Ok(w)
}

/// The HTML report's card, when step 25 has a result: one row per score,
/// with the effect column when the table names evaluations.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("prs");
    if !sec.state.ran() {
        return Ok(HtmlPart::default());
    }
    let d = sec.data;
    let rs = scores(d)?;
    let adjusted = get(d, "adjusted", PATH)?.is_some_and(truthy);
    let show_effect = adjusted && names_evaluations(&rs)?;
    let mut rows = Vec::new();
    for x in &rs {
        let pct = percentile(x)?;
        let mut row = vec![
            Cell::Text(need_text(x, "condition", ROW)?),
            Cell::Text(need_text(x, "pgs_id", ROW)?),
            Cell::Text(count(x, "score")?),
            Cell::Text(format!("{}/{}", count(x, "matched")?, count(x, "total")?)),
            Cell::Text(if !pct.is_empty() {
                pct
            } else if matches!(need(x, "score", ROW)?, Value::Null) {
                NOT_AVAILABLE.to_owned()
            } else {
                "raw score only".to_owned()
            }),
        ];
        if show_effect {
            row.push(Cell::Text(effect(x)?.0));
        }
        rows.push(row);
    }
    let mut head = [
        "Condition",
        "PGS",
        "Score",
        "Variants matched",
        "Percentile",
    ]
    .map(str::to_owned)
    .to_vec();
    if show_effect {
        head.push("Effect in your group".to_owned());
    }
    let card = Card {
        title: format!(
            "Polygenic Risk Scores{}",
            if adjusted {
                ""
            } else {
                " (raw, not percentiles)"
            }
        ),
        full: true,
        body: vec![
            Line::Table(Table {
                head,
                rows,
                layout: Layout::RowLines,
                ..Table::default()
            }),
            Line::Para(note(d)?),
        ],
        ..Card::default()
    };
    Ok(HtmlPart {
        cards: vec![html::card(&sec, card)?],
        files: Vec::new(),
    })
}
