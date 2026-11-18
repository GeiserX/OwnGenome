//! Fine-grained ancestry (step 26b): the continental and regional groups
//! with their bootstrap ranges, the groups the model leaves out, the share
//! of the genome each group paints, and the chromosome painting with its
//! legend.

use pycompat::json::Value;

use crate::escape::e;
use crate::model::{
    Card, HtmlPart, Inline, Line, Malformed, State, Summary, get, need, py, truthy,
};
use crate::notes::ANCESTRY_FINE_NOTE;
use crate::sections::ancestry::{get_text, has, list, need_text, number, share};
use crate::sections::haplogroups::shows_body;
use crate::{html, text};

const PATH: &str = "sections.ancestry_fine.data";

/// The two levels the model reports, with their headings.
pub const LEVELS: [(&str, &str); 2] = [
    ("continent", "Continental groups"),
    ("region", "Regional groups"),
];

/// The groups of one level worth a line: every continental group of 1% or
/// more, a regional one of 1% or more or whose bootstrap range reaches 5%.
pub fn rows(d: &Value, level: &str) -> Result<Vec<Value>, Malformed> {
    let at = format!("{PATH}.{level}[]");
    let mut out = Vec::new();
    for r in list(d, level, PATH)? {
        let f = need(&r, "fraction", &at)?;
        if matches!(f, Value::Null) {
            continue;
        }
        let keep = number(f, &format!("{at}.fraction"))? >= 0.01
            || (level == "region"
                && match get(&r, "high", &at)?.filter(|v| truthy(v)) {
                    Some(h) => number(h, &format!("{at}.high"))? >= 0.05,
                    None => false,
                });
        if keep {
            out.push(r);
        }
    }
    Ok(out)
}

/// ` (93% to 97%)`: the bootstrap range of a row, when it has one.
pub fn spread(r: &Value, at: &str) -> Result<String, Malformed> {
    match get(r, "low", at)? {
        None | Some(Value::Null) => Ok(String::new()),
        Some(low) => Ok(format!(
            " ({} to {})",
            share(Some(low), at)?,
            share(Some(need(r, "high", at)?), at)?
        )),
    }
}

/// `(.*) \((.*) group\)` matched against the whole name: the group and its
/// kind (`Oceania (continental group)`), or the name alone.
fn name_and_kind(group: &str) -> (&str, &str) {
    let Some(inner) = group.strip_suffix(" group)") else {
        return (group, "");
    };
    match inner.rfind(" (") {
        Some(i) if !inner.contains('\n') => (&inner[..i], &inner[i + 2..]),
        _ => (group, ""),
    }
}

/// The groups the model leaves out, as one sentence, or "".
pub fn left_out(d: &Value) -> Result<String, Malformed> {
    if !has(d, "left_out", PATH)? {
        return Ok(String::new());
    }
    let at = format!("{PATH}.left_out[]");
    // One row per level ("Oceania (continental group)"); a group both
    // levels leave out reads as one name.
    let mut kinds: Vec<(String, Vec<String>)> = Vec::new();
    for r in list(d, "left_out", PATH)? {
        let Value::Str(group) = need(&r, "group", &at)? else {
            return Err(Malformed(format!("{at}.group is not a text")));
        };
        let (name, kind) = name_and_kind(group);
        let i = match kinds.iter().position(|(n, _)| n == name) {
            Some(i) => i,
            None => {
                kinds.push((name.to_owned(), Vec::new()));
                kinds.len() - 1
            }
        };
        if !kind.is_empty() && !kinds[i].1.iter().any(|k| k == kind) {
            kinds[i].1.push(kind.to_owned());
        }
    }
    let names: Vec<String> = kinds
        .iter()
        .map(|(n, k)| {
            if k.is_empty() {
                n.clone()
            } else {
                format!("{n} ({} group)", k.join(" and "))
            }
        })
        .collect();
    Ok(format!(
        "Left out of the model: {}. The panel has fewer than {} unrelated reference people for \
         it, so its ancestry is spread over the groups shown.",
        names.join(", "),
        get_text(d, "min_references", "10", PATH)?
    ))
}

/// The painted shares, or why there are none, and the chromosomes the
/// painting failed on.
pub fn painted(d: &Value) -> Result<String, Malformed> {
    let codes = list(d, "not_painted", PATH)?
        .iter()
        .map(
            |r| match need(r, "code", &format!("{PATH}.not_painted[]"))? {
                Value::Str(c) => Ok(c.clone()),
                _ => Err(Malformed(format!(
                    "{PATH}.not_painted[].code is not a text"
                ))),
            },
        )
        .collect::<Result<Vec<_>, _>>()?;
    let failed = codes.join(", ");
    let mut text = if has(d, "painting", PATH)? {
        let at = format!("{PATH}.painting[]");
        list(d, "painting", PATH)?
            .iter()
            .map(|r| {
                Ok(format!(
                    "{} {}",
                    need_text(r, "group", &at)?,
                    share(Some(need(r, "fraction", &at)?), &at)?
                ))
            })
            .collect::<Result<Vec<_>, Malformed>>()?
            .join(", ")
    } else if failed.is_empty() {
        "none (no chromosome long enough)".to_owned()
    } else {
        "none".to_owned()
    };
    if !failed.is_empty() {
        text += &format!("; not painted: RFMix failed on {failed}");
    }
    Ok(text)
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(
        s,
        "ancestry_fine",
        Some("Fine-Grained Ancestry (ADMIXTURE, RFMix)"),
    )?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = s.section("ancestry_fine").data;
    for (level, title) in LEVELS {
        let rs = rows(d, level)?;
        if rs.is_empty() {
            continue;
        }
        w.push(format!("  {title} (fraction, bootstrap range):"));
        let at = format!("{PATH}.{level}[]");
        for r in &rs {
            let group = crate::sections::ancestry::need_align(
                r,
                "group",
                50,
                crate::sections::ancestry::Align::Left,
                &at,
            )?;
            let frac = share(Some(need(r, "fraction", &at)?), &at)?;
            w.push(format!(
                "    {group} {}{}",
                pycompat::fmt::rjust(&frac, 5),
                spread(r, &at)?
            ));
        }
    }
    let lo = left_out(d)?;
    if !lo.is_empty() {
        w.push(format!("  {lo}"));
    }
    let lead = if has(d, "painting", PATH)? {
        "  Chromosome painting, share of the painted genome: "
    } else {
        "  Chromosome painting: "
    };
    w.push(format!("{lead}{}", painted(d)?));
    w.push(format!("  {ANCESTRY_FINE_NOTE}"));
    w.push(String::new());
    Ok(w)
}

/// The HTML report's card, when step 26b has a result: one table per
/// level, the painted shares, and the painting with its legend.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("ancestry_fine");
    if sec.state == State::Missing {
        return Ok(HtmlPart::default());
    }
    let mut c = Card {
        title: "Fine-Grained Ancestry".to_owned(),
        full: true,
        ..Card::default()
    };
    if shows_body(&sec) {
        let d = sec.data;
        for (level, title) in LEVELS {
            let at = format!("{PATH}.{level}[]");
            let mut rows_html = String::new();
            for r in rows(d, level)? {
                let range = spread(&r, &at)?;
                let range = range.trim_matches([' ', '(', ')']);
                rows_html += &format!(
                    "<tr><td>{}</td><td>{}</td><td>{}</td></tr>\n",
                    e(&need_text(&r, "group", &at)?),
                    e(&share(Some(need(&r, "fraction", &at)?), &at)?),
                    e(if range.is_empty() { "." } else { range })
                );
            }
            if !rows_html.is_empty() {
                c.body.extend([
                    Line::Raw(format!("    <h3>{}</h3>", e(title))),
                    Line::Raw("    <table>".to_owned()),
                    Line::Raw(
                        "      <tr><th>Group</th><th>Fraction</th><th>Bootstrap range</th></tr>"
                            .to_owned(),
                    ),
                    Line::Raw(rows_html.trim_end_matches('\n').to_owned()),
                    Line::Raw("    </table>".to_owned()),
                ]);
            }
        }
        let lo = left_out(d)?;
        if !lo.is_empty() {
            c.body.push(Line::Para(lo));
        }
        c.body.push(Line::Stat(
            "Painted genome".to_owned(),
            Inline::Text(painted(d)?),
        ));
        if let Some(png) = get(d, "painting_png", PATH)?.filter(|v| truthy(v)) {
            c.body.push(Line::Raw(format!(
                "    <img alt=\"Chromosome painting: chromosomes 1 to 22, two haplotypes each, \
                 coloured by continental group\" style=\"max-width:100%\" \
                 src=\"data:image/png;base64,{}\">",
                e(&py(png)?)
            )));
            let at = format!("{PATH}.legend[]");
            let keys = list(d, "legend", PATH)?
                .iter()
                .map(|x| {
                    Ok(format!(
                        "<span style=\"color:{}\">&#9632;</span> {}",
                        css_colour(&need_text(x, "colour", &at)?),
                        e(&need_text(x, "group", &at)?)
                    ))
                })
                .collect::<Result<Vec<_>, Malformed>>()?;
            c.body.push(Line::Raw(format!(
                "    <p>{} <span style=\"color:#e1e1e1\">&#9632;</span> not painted</p>",
                keys.join(" ")
            )));
        }
        c.caveats.push(ANCESTRY_FINE_NOTE);
    }
    Ok(HtmlPart {
        cards: vec![html::card(&sec, c)?],
        files: Vec::new(),
    })
}

/// A legend colour as CSS: the `#rrggbb` the collector writes, else a fixed
/// grey. Escaping keeps it inside the attribute but would still let a `;`
/// add a declaration, such as one that hides the swatch.
fn css_colour(c: &str) -> &str {
    let hex = c.len() == 7 && c.starts_with('#') && c[1..].bytes().all(|b| b.is_ascii_hexdigit());
    if hex { c } else { "#808080" }
}
