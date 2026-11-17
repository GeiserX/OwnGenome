//! Step 09's repeat expansions (ExpansionHunter) with Stranger's status
//! (step 9b): the loci outside the normal range, the loci it does not
//! grade, the readings of the repeats STRchive leaves out, and every locus
//! with its count, status and notes, the key loci first.

use pycompat::fmt::ljust;
use pycompat::json::Value;

use crate::escape::e;
use crate::html;
use crate::model::{
    Card, Colour, HtmlPart, Inline, Line, Malformed, Summary, get, is_str, need, py, truthy,
};
use crate::notes::{EH_NOTE, NOT_GRADED_NOTE, RFC1_NOTE};
use crate::sections::paralogs::{done_badge, shows_body};
use crate::sections::rfc1::{motif_text, rfc1_text};
use crate::text;

const PATH: &str = "sections.expansions.data";

/// What the report says of an allele that was not sized from spanning
/// reads.
fn sizing_text(sizing: &str) -> Option<&'static str> {
    Some(match sizing {
        "flanking_only" => "size not determined (flanking reads only)",
        "in_repeat" => {
            "longer than a read (in-repeat reads only): shown as at least the interval's lower end"
        }
        "low_confidence" => {
            "low confidence (fewer than 2 spanning reads, or the interval reaches the disease range)"
        }
        "no_call" => "no call",
        _ => return None,
    })
}

/// `l.get(key)`.
fn at<'a>(l: &'a Value, key: &str) -> Result<Option<&'a Value>, Malformed> {
    get(l, key, &format!("{PATH}.loci[]"))
}

/// `(l.get(key) or "")` as text.
fn text_or_empty(l: &Value, key: &str) -> Result<String, Malformed> {
    match at(l, key)?.filter(|v| truthy(v)) {
        Some(v) => py(v),
        None => Ok(String::new()),
    }
}

/// `d.get(key) or []`.
fn list<'a>(d: &'a Value, key: &str) -> Result<&'a [Value], Malformed> {
    match get(d, key, PATH)?.filter(|v| truthy(v)) {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!("{PATH}.{key} is not a list"))),
    }
}

/// The repeat count of a locus: an in-repeat allele as at least the lower
/// end of its interval, or `not sized` without one.
pub fn count(l: &Value) -> Result<String, Malformed> {
    let counts = py(need(l, "repeat_count", &format!("{PATH}.loci[]"))?)?;
    let counts: Vec<&str> = counts.split('/').collect();
    let sizing = text_or_empty(l, "allele_sizing")?;
    let sizing: Vec<&str> = sizing.split('/').collect();
    let cis = text_or_empty(l, "ci")?;
    let cis: Vec<&str> = cis.split('/').collect();
    let mut out = Vec::new();
    for (i, c) in counts.iter().enumerate() {
        let z = sizing.get(i).copied().unwrap_or("");
        let ci = if cis.len() == counts.len() {
            cis[i]
        } else if cis.len() == 1 {
            cis[0]
        } else {
            ""
        };
        let lo = ci.split('-').next().unwrap_or("");
        out.push(if z == "in_repeat" {
            if pycompat::text::isdigit(lo) {
                format!(">={lo}")
            } else {
                "not sized".to_owned()
            }
        } else {
            (*c).to_owned()
        });
    }
    Ok(out.join("/"))
}

/// Stranger's status, never `normal` for a locus that was not sized.
pub fn status(l: &Value) -> Result<String, Malformed> {
    let st = at(l, "status")?.filter(|v| truthy(v));
    if st.is_none() && at(l, "not_graded")?.is_some_and(truthy) {
        return Ok("not graded".to_owned());
    }
    if is_str(st, "normal")
        && !is_str(at(l, "sizing")?, "sized")
        && !is_str(at(l, "sizing")?, "unknown")
    {
        return Ok("not sized".to_owned());
    }
    st.map_or(Ok(String::new()), py)
}

/// What the report says beside a locus: sizing, the catalog's normal
/// maximum, an intermediate range, no thresholds, the motif count.
pub fn notes(l: &Value) -> Result<String, Malformed> {
    let mut out = Vec::new();
    if let Some(Value::Str(z)) = at(l, "sizing")?
        && let Some(t) = sizing_text(z)
    {
        let ci = at(l, "ci")?;
        let shown = match ci {
            None | Some(Value::Null) => false,
            Some(v) => !is_str(Some(v), "") && !is_str(Some(v), "."),
        };
        out.push(if shown {
            format!("{t}, interval {}", py(ci.expect("checked above"))?)
        } else {
            t.to_owned()
        });
    }
    if let Some(m) = at(l, "normal_max")?.filter(|v| truthy(v)) {
        out.push(format!("normal max {}", py(m)?));
    }
    if let Some(z) = at(l, "zone")?.filter(|v| truthy(v)) {
        let path = format!("{PATH}.loci[]");
        out.push(format!(
            "{}: {}; {}",
            py(z)?,
            py(need(l, "zone_label", &path)?)?,
            py(need(l, "zone_citation", &path)?)?
        ));
    }
    if let Some(n) = at(l, "not_graded")?.filter(|v| truthy(v)) {
        out.push(format!("no thresholds in the catalog ({})", py(n)?));
    }
    if let Some(m) = at(l, "motif")?.filter(|v| truthy(v)) {
        let locus = py(need(l, "locus", &format!("{PATH}.loci[]"))?)?;
        out.push(motif_text(&locus, m)?);
    }
    Ok(out.join("; "))
}

/// `(locus, count, status, notes)` for every locus, the key loci first; a
/// key locus missing from the output says so.
pub fn rows(d: &Value) -> Result<Vec<[String; 4]>, Malformed> {
    let mut rows = Vec::new();
    for l in list(d, "loci")? {
        rows.push([
            py(need(l, "locus", &format!("{PATH}.loci[]"))?)?,
            count(l)?,
            status(l)?,
            notes(l)?,
        ]);
    }
    let mut out = Vec::new();
    for k in list(d, "key_loci")? {
        let locus = need(k, "locus", &format!("{PATH}.key_loci[]"))?;
        let name = py(locus)?;
        // `k["locus"] not in have`: compared by value, as Python compares.
        let have = list(d, "loci")?.iter().any(|l| {
            get(l, "locus", PATH)
                .ok()
                .flatten()
                .is_some_and(|x| crate::model::py_eq(x, locus))
        });
        if !have {
            out.push([
                name,
                "not in output".to_owned(),
                String::new(),
                String::new(),
            ]);
        }
    }
    out.extend(rows);
    Ok(out)
}

/// The first line of one reading of the repeats STRchive leaves out (CFTR,
/// UGT1A1, 5-HTTLPR).
fn reading_head(r: &Value) -> Result<String, Malformed> {
    let path = format!("{PATH}.readings[]");
    let f = |k: &str| need(r, k, &path).and_then(py);
    Ok(format!(
        "{} {}: {}",
        f("locus")?,
        f("repeat")?,
        reading_what(r)?
    ))
}

/// A reading's genotype when it was read, else its status in words.
fn reading_what(r: &Value) -> Result<String, Malformed> {
    let path = format!("{PATH}.readings[]");
    let st = need(r, "status", &path)?;
    if is_str(Some(st), "read") {
        py(need(r, "genotype", &path)?)
    } else {
        Ok(py(st)?.replace('_', " "))
    }
}

/// The joined items of a list or the keys of a dict, as `", ".join(x)`.
fn joined(v: &Value, key: &str) -> Result<String, Malformed> {
    match v {
        Value::Array(a) => Ok(a.iter().map(py).collect::<Result<Vec<_>, _>>()?.join(", ")),
        Value::Object(o) => Ok(o.keys().cloned().collect::<Vec<_>>().join(", ")),
        Value::Str(s) => Ok(s.chars().map(String::from).collect::<Vec<_>>().join(", ")),
        _ => Err(Malformed(format!("{PATH}.{key} is not a list"))),
    }
}

/// `len(x)` of a list, a dict or a string.
fn len(v: &Value, key: &str) -> Result<usize, Malformed> {
    match v {
        Value::Array(a) => Ok(a.len()),
        Value::Object(o) => Ok(o.len()),
        Value::Str(s) => Ok(s.chars().count()),
        _ => Err(Malformed(format!("{PATH}.{key} has no length"))),
    }
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "expansions", Some("Repeat Expansions (ExpansionHunter)"))?
    else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let sec = s.section("expansions");
    let d = sec.data;
    let g = |k: &str| get(d, k, PATH);
    w.push(format!("  Loci in output: {}", py(sec.need("records")?)?));
    let flagged = list(d, "flagged")?;
    if g("stranger_outdated")?.is_some_and(truthy) {
        w.push(
            "  Stranger's result is older than ExpansionHunter's: rerun step 9b. No locus is \
             marked as normal or expanded."
                .to_owned(),
        );
    } else if !g("stranger")?.is_some_and(truthy) {
        w.push(
            "  Stranger (step 9b) did not run: no locus is marked as normal or expanded."
                .to_owned(),
        );
    } else if !flagged.is_empty() {
        w.push(format!(
            "  Outside the normal range (Stranger): {}",
            len(sec.need("flagged")?, "flagged")?
        ));
        let loci = list(d, "loci")?;
        for l in flagged {
            let path = format!("{PATH}.flagged[]");
            let locus = need(l, "locus", &path)?;
            // `by.get(l["locus"], l)`: the full locus row when there is one;
            // the last row of a locus listed twice, as a dict keeps it.
            let full = loci
                .iter()
                .rev()
                .find(|x| {
                    get(x, "locus", PATH)
                        .ok()
                        .flatten()
                        .is_some_and(|y| crate::model::py_eq(y, locus))
                })
                .unwrap_or(l);
            w.push(format!(
                "    {} {}  {}",
                ljust(&py(locus)?, 8),
                count(full)?,
                py(need(l, "status", &path)?)?
            ));
            if let Some(m) = get(l, "rfc1_motif", &path)?.filter(|v| truthy(v)) {
                w.push(format!("      {}", rfc1_text(m)?));
            }
            if let Some(m) = get(l, "motif", &path)?.filter(|v| truthy(v)) {
                w.push(format!("      {}", motif_text(&py(locus)?, m)?));
            }
        }
    } else {
        w.push("  Outside the normal range (Stranger): none".to_owned());
    }
    if let Some(ng) = g("not_graded")?.filter(|v| truthy(v)) {
        let Value::Object(o) = ng else {
            return Err(Malformed(format!("{PATH}.not_graded is not an object")));
        };
        let items: Vec<String> = o
            .iter()
            .map(|(k, v)| py(v).map(|v| format!("{k}: {v}")))
            .collect::<Result<_, _>>()?;
        w.push(format!(
            "  Not graded (no thresholds in the catalog): {} ({})",
            o.len(),
            items.join("; ")
        ));
        w.push(format!("  {NOT_GRADED_NOTE}"));
    }
    for r in list(d, "readings")? {
        let path = format!("{PATH}.readings[]");
        w.push(format!("  {}", reading_head(r)?));
        w.push(format!(
            "    {} ({})",
            py(need(r, "reading", &path)?)?,
            py(need(r, "source", &path)?)?
        ));
    }
    if g("stranger")?.is_some_and(truthy)
        && let Some(n) = g("no_status")?.filter(|v| truthy(v))
    {
        w.push(format!(
            "  Loci without a Stranger status (not in its catalog): {}",
            py(n)?
        ));
    }
    if let Some(u) = g("unsized")?.filter(|v| truthy(v)) {
        w.push(format!(
            "  Not sized from short reads, so not counted as normal: {} ({})",
            len(u, "unsized")?,
            joined(u, "unsized")?
        ));
    }
    if let Some(i) = g("intermediate")?.filter(|v| truthy(v)) {
        w.push(format!(
            "  In an intermediate range (bin/str_intermediate_ranges.tsv): {}",
            joined(i, "intermediate")?
        ));
    }
    let rows = rows(d)?;
    if !rows.is_empty() {
        w.push("  Every locus (key loci first): repeat count, Stranger status, notes".to_owned());
        for [locus, count, status, notes] in &rows {
            let status = if status.is_empty() { "-" } else { status };
            let tail = if notes.is_empty() {
                String::new()
            } else {
                format!("  {notes}")
            };
            w.push(format!(
                "    {} {} {status}{tail}",
                ljust(locus, 8),
                ljust(count, 10)
            ));
        }
    }
    if flagged.iter().any(|l| {
        get(l, "rfc1_motif", PATH)
            .ok()
            .flatten()
            .is_some_and(truthy)
    }) {
        w.push(format!("  {RFC1_NOTE}"));
    }
    w.push(format!("  {EH_NOTE}"));
    w.push(String::new());
    Ok(w)
}

/// The HTML report's card, shown on every report: "Not run" when step 09
/// is not in the summary.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("expansions");
    let mut body = vec![Line::Stat("Status".to_owned(), done_badge(&sec))];
    let mut caveats = Vec::new();
    if shows_body(&sec) {
        let d = sec.data;
        let g = |k: &str| get(d, k, PATH);
        let label = "Outside the normal range (Stranger)".to_owned();
        let flagged = list(d, "flagged")?;
        body.push(if g("stranger_outdated")?.is_some_and(truthy) {
            Line::Stat(
                label,
                Inline::Badge(
                    "older than ExpansionHunter's result: rerun step 9b".to_owned(),
                    Colour::Yellow,
                ),
            )
        } else if !g("stranger")?.is_some_and(truthy) {
            Line::Stat(label, Inline::Badge("not run".to_owned(), Colour::Gray))
        } else {
            let colour = if flagged.is_empty() {
                Colour::Gray
            } else {
                Colour::Yellow
            };
            Line::Stat(label, Inline::Badge(flagged.len().to_string(), colour))
        });
        if let Some(u) = g("unsized")?.filter(|v| truthy(v)) {
            body.push(Line::Stat(
                "Not sized from short reads (not counted as normal)".to_owned(),
                Inline::Badge(len(u, "unsized")?.to_string(), Colour::Yellow),
            ));
        }
        if let Some(i) = g("intermediate")?.filter(|v| truthy(v)) {
            body.push(Line::Stat(
                "In an intermediate range".to_owned(),
                Inline::Badge(joined(i, "intermediate")?, Colour::Yellow),
            ));
        }
        if let Some(n) = g("not_graded")?.filter(|v| truthy(v)) {
            body.push(Line::Stat(
                "Not graded (no thresholds in the catalog)".to_owned(),
                Inline::Badge(joined(n, "not_graded")?, Colour::Gray),
            ));
            caveats.push(NOT_GRADED_NOTE);
        }
        for r in list(d, "readings")? {
            let path = format!("{PATH}.readings[]");
            let colour = if get(r, "five_t", &path)?.is_some_and(truthy)
                || is_str(get(r, "agreement", &path)?, "differs")
            {
                Colour::Yellow
            } else {
                Colour::Gray
            };
            body.push(Line::Stat(
                format!(
                    "{} {}",
                    py(need(r, "locus", &path)?)?,
                    py(need(r, "repeat", &path)?)?
                ),
                Inline::Badge(reading_what(r)?, colour),
            ));
            body.push(Line::Note(format!(
                "{} ({})",
                py(need(r, "reading", &path)?)?,
                py(need(r, "source", &path)?)?
            )));
        }
        for l in flagged {
            let path = format!("{PATH}.flagged[]");
            if let Some(m) = get(l, "rfc1_motif", &path)?.filter(|v| truthy(v)) {
                body.push(Line::Note(rfc1_text(m)?));
                caveats.push(RFC1_NOTE);
            }
            if let Some(m) = get(l, "motif", &path)?.filter(|v| truthy(v)) {
                body.push(Line::Note(motif_text(&py(need(l, "locus", &path)?)?, m)?));
            }
        }
        let rows = rows(d)?;
        if !rows.is_empty() {
            let mut t = "    <table><tr><th>Locus</th><th>Repeat Count</th><th>Stranger</th>\
                         <th>Notes</th></tr>"
                .to_owned();
            for r in &rows {
                t += "<tr>";
                for x in r {
                    t += &format!("<td>{}</td>", e(x));
                }
                t += "</tr>";
            }
            t += "</table>";
            body.push(Line::Markup(t));
        }
        caveats.push(EH_NOTE);
    }
    let card = html::card(
        &sec,
        Card {
            title: "Repeat Expansions".to_owned(),
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
