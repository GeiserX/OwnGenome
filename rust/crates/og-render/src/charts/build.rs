//! The five chart specs from a summary, read with Python's rules (an int
//! stays an int, `str()`, `int()` and `float()` as Python applies them), so
//! the specs are what `bin/report_charts.py` built.

use indexmap::IndexMap;
use pycompat::json::Value;
use pycompat::num::{parse_float, parse_int};
use pycompat::value::py_repr;

use super::{Kind, Mark, MarkKind, NotDrawn, ROH_MIN_MB, Spec, Tone, chrom_len, draw};
use crate::model::{Malformed, Summary, get, need, py, py_eq, truthy};
use crate::notes;
use crate::sections::clinvar::{HITS_ANCHOR, HITS_SHOWN};
use crate::sections::{ancestry_fine, prs};

/// The specs of the five charts, in report order (ClinVar, genome, mito,
/// fine-grained ancestry, PRS), for the sections that ran, each drawn.
pub fn build(s: &Summary) -> Result<Vec<Spec>, Malformed> {
    let all = [
        clinvar(s)?,
        genome(s)?,
        mito(s)?,
        ancestry(s)?,
        prs_chart(s)?,
    ];
    let mut out = Vec::new();
    for mut sp in all.into_iter().flatten() {
        sp.svg = draw(&sp)?;
        out.push(sp);
    }
    Ok(out)
}

// --- Python's rules ---------------------------------------------------------------------

fn bad(path: &str, why: &str) -> Malformed {
    Malformed(format!("{path}: {why}"))
}

/// `int(v)`: an int, a float cut toward zero, a bool, or a text of digits.
fn py_int(v: &Value, path: &str) -> Result<i128, Malformed> {
    match v {
        Value::Int(i) => Ok(*i),
        Value::Bool(b) => Ok(i128::from(*b)),
        Value::Float(x) if x.is_finite() => Ok(x.trunc() as i128),
        Value::Str(t) => parse_int(t).map_err(|_| bad(path, "not an int")),
        _ => Err(bad(path, "not an int")),
    }
}

/// `float(v)`; `Ok(None)` for a text that is not a number (Python's
/// ValueError), an error for a value of another kind (TypeError).
fn py_float(v: &Value, path: &str) -> Result<Option<f64>, Malformed> {
    match v {
        Value::Int(i) => Ok(Some(*i as f64)),
        Value::Float(x) => Ok(Some(*x)),
        Value::Bool(b) => Ok(Some(f64::from(u8::from(*b)))),
        Value::Str(t) => Ok(parse_float(t).ok()),
        _ => Err(bad(path, "not a number")),
    }
}

/// `float(v)` that must be a number.
fn need_float(v: &Value, path: &str) -> Result<f64, Malformed> {
    py_float(v, path)?.ok_or_else(|| bad(path, "not a number"))
}

/// `100 * v`, an int for an int as in Python.
fn times_100(v: &Value, path: &str) -> Result<Value, Malformed> {
    match v {
        Value::Int(i) => Ok(Value::Int(i * 100)),
        Value::Bool(b) => Ok(Value::Int(i128::from(*b) * 100)),
        Value::Float(x) => Ok(Value::Float(100.0 * x)),
        _ => Err(bad(path, "not a number")),
    }
}

/// `str(d.get(key, default))`.
fn s_or(d: &Value, key: &str, default: &str, path: &str) -> Result<String, Malformed> {
    match get(d, key, path)? {
        Some(v) => py(v),
        None => Ok(default.to_owned()),
    }
}

/// `str(d.get(key) or default)`.
fn s_truthy(d: &Value, key: &str, default: &str, path: &str) -> Result<String, Malformed> {
    match get(d, key, path)?.filter(|v| truthy(v)) {
        Some(v) => py(v),
        None => Ok(default.to_owned()),
    }
}

/// `d.get(key)`, `null` when it is not there.
fn raw(d: &Value, key: &str, path: &str) -> Result<Value, Malformed> {
    Ok(get(d, key, path)?.cloned().unwrap_or(Value::Null))
}

/// `d.get(key) or []` as a list.
fn list<'a>(d: &'a Value, key: &str, path: &str) -> Result<&'a [Value], Malformed> {
    match get(d, key, path)?.filter(|v| truthy(v)) {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(bad(&format!("{path}.{key}"), "not a list")),
    }
}

/// The section's data, when the section ran (`ok` or `stale`).
fn ran<'a>(s: &'a Summary, key: &'a str) -> Option<&'a Value> {
    let sec = s.section(key);
    sec.state.ran().then_some(sec.data)
}

fn caveat(key: &'static str) -> Result<&'static str, Malformed> {
    notes::text(key).ok_or_else(|| Malformed(format!("no caveat called {key}")))
}

fn spec(
    id: &'static str,
    section: &'static str,
    kind: Kind,
    title: &'static str,
    caveat_key: &'static str,
    marks: Vec<Mark>,
) -> Result<Spec, Malformed> {
    Ok(Spec {
        id,
        section,
        kind,
        title,
        caveat_key,
        caveat: caveat(caveat_key)?,
        axis: IndexMap::new(),
        marks,
        not_drawn: Vec::new(),
        text: None,
        svg: String::new(),
    })
}

fn axis(items: &[(&str, Value)]) -> IndexMap<String, Value> {
    items
        .iter()
        .map(|(k, v)| ((*k).to_owned(), v.clone()))
        .collect()
}

fn plural(n: i128) -> &'static str {
    if n == 1 { "" } else { "s" }
}

/// The report's hits table for a hit it lists, the ClinVar section for one
/// past it.
fn hit_anchor(i: usize) -> &'static str {
    if i < HITS_SHOWN {
        HITS_ANCHOR
    } else {
        "clinvar"
    }
}

/// `GENE chrom:pos REF>ALT`.
fn hit_label(h: &Value, path: &str) -> Result<String, Malformed> {
    Ok(format!(
        "{} {}:{} {}>{}",
        s_truthy(h, "gene", ".", path)?,
        py(need(h, "chrom", path)?)?,
        py(need(h, "pos", path)?)?,
        s_or(h, "ref", "", path)?,
        s_or(h, "alt", "", path)?
    ))
}

/// `^(chr[0-9XYM]+):([0-9]+)-([0-9]+)$`, where `$` also matches before a
/// final newline as in Python.
fn region(text: &str) -> Option<(String, i128, i128)> {
    let t = text.strip_suffix('\n').unwrap_or(text);
    let (chrom, rest) = t.split_once(':')?;
    let name = chrom.strip_prefix("chr")?;
    if name.is_empty()
        || !name
            .bytes()
            .all(|b| b.is_ascii_digit() || b"XYM".contains(&b))
    {
        return None;
    }
    let (a, b) = rest.split_once('-')?;
    let digits = |x: &str| !x.is_empty() && x.bytes().all(|c| c.is_ascii_digit());
    if !digits(a) || !digits(b) {
        return None;
    }
    Some((chrom.to_owned(), a.parse().ok()?, b.parse().ok()?))
}

/// Counts by what is not drawn, in the order first seen.
#[derive(Default)]
struct Tally(IndexMap<String, i128>);

impl Tally {
    fn add(&mut self, what: &str) {
        *self.0.entry(what.to_owned()).or_insert(0) += 1;
    }

    fn list(self) -> Vec<NotDrawn> {
        self.0
            .into_iter()
            .map(|(what, count)| NotDrawn { what, count })
            .collect()
    }
}

// --- the five charts --------------------------------------------------------------------

fn clinvar(s: &Summary) -> Result<Option<Spec>, Malformed> {
    let Some(c) = ran(s, "clinvar") else {
        return Ok(None);
    };
    let path = "sections.clinvar.data";
    let empty = Value::Object(IndexMap::new());
    let by = get(c, "by_stars", path)?
        .filter(|v| truthy(v))
        .unwrap_or(&empty);
    let hits = list(c, "hits", path)?;
    let hp = format!("{path}.hits[]");
    let low_conf = |h: &Value| -> Result<bool, Malformed> {
        Ok(get(h, "low_confidence", &hp)?.is_some_and(truthy))
    };
    let mut marks = Vec::new();
    for k in ["4", "3", "2", "1", "0"] {
        let n = match get(by, k, &format!("{path}.by_stars"))? {
            Some(v) => py_int(v, &format!("{path}.by_stars.{k}"))?,
            None => 0,
        };
        let mut idx = Vec::new();
        for (i, h) in hits.iter().enumerate() {
            let stars = match get(h, "stars", &hp)? {
                Some(v) => py(v)?,
                None => "None".to_owned(),
            };
            if stars == k && !low_conf(h)? {
                idx.push(i);
            }
        }
        let shown: Vec<String> = idx
            .iter()
            .take(10)
            .map(|&i| hit_label(&hits[i], &hp))
            .collect::<Result<_, _>>()?;
        let names = shown.join(", ") + if idx.len() > 10 { " and more" } else { "" };
        let one = if k == "1" { "" } else { "s" };
        let mut m = Mark::new(
            format!("stars-{k}"),
            format!("{k} star{one}"),
            Value::Int(n),
        );
        m.tone = if n != 0 { Tone::Attention } else { Tone::Muted };
        m.anchor = Some(idx.first().map_or("clinvar", |&i| hit_anchor(i)).to_owned());
        m.detail = format!("{n} hit{} with {k} review star{one}", plural(n))
            + &if names.is_empty() {
                String::new()
            } else {
                format!(": {names}")
            };
        marks.push(m);
    }
    if let Some(v) = get(c, "low_confidence", path)? {
        let n = py_int(v, &format!("{path}.low_confidence"))?;
        let mut first = None;
        for (i, h) in hits.iter().enumerate() {
            if low_conf(h)? {
                first = Some(i);
                break;
            }
        }
        let f = get(c, "call_floor", path)?
            .filter(|v| truthy(v))
            .unwrap_or(&empty);
        let fp = format!("{path}.call_floor");
        let floor = |k: &str| -> Result<String, Malformed> {
            get(f, k, &fp)?.map_or(Ok("None".to_owned()), py)
        };
        let mut m = Mark::new(
            "low-confidence".into(),
            "Low-confidence, not counted".into(),
            Value::Int(n),
        );
        m.tone = Tone::Muted;
        m.anchor = Some(first.map_or("clinvar", hit_anchor).to_owned());
        m.detail = format!(
            "{n} call{} with GQ below {} or DP below {}: listed in the hits, not counted in the stars",
            plural(n),
            floor("GQ")?,
            floor("DP")?
        );
        marks.push(m);
    }
    let mut sp = spec(
        "clinvar-stars",
        "clinvar",
        Kind::Bars,
        "ClinVar hits by review stars",
        "CLINVAR_FREQ_NOTE",
        marks,
    )?;
    sp.axis = axis(&[("unit", Value::Str("hits".into()))]);
    Ok(Some(sp))
}

fn ancestry(s: &Summary) -> Result<Option<Spec>, Malformed> {
    let Some(d) = ran(s, "ancestry_fine") else {
        return Ok(None);
    };
    let mut marks: Vec<Mark> = Vec::new();
    for (level, title) in ancestry_fine::LEVELS {
        let at = format!("sections.ancestry_fine.data.{level}[]");
        for r in ancestry_fine::rows(d, level)? {
            let v = times_100(need(&r, "fraction", &at)?, &format!("{at}.fraction"))?;
            let bound = |k: &str| -> Result<Value, Malformed> {
                match get(&r, k, &at)? {
                    None | Some(Value::Null) => Ok(Value::Null),
                    Some(x) => times_100(x, &format!("{at}.{k}")),
                }
            };
            let (low, high) = (bound("low")?, bound("high")?);
            let pct =
                |x: &Value| pycompat::fmt::fixed(crate::model::num(Some(x)).unwrap_or(0.0), 0);
            let code = match get(&r, "code", &at)?.filter(|c| truthy(c)) {
                Some(c) => py(c)?,
                None => marks.len().to_string(),
            };
            let mut m = Mark::new(
                format!("{level}-{code}"),
                py(need(&r, "group", &at)?)?,
                v.clone(),
            );
            m.detail = format!("{title}: {}%", pct(&v));
            if !matches!(low, Value::Null) && !matches!(high, Value::Null) {
                m.detail += &format!(" (bootstrap range {}% to {}%)", pct(&low), pct(&high));
            }
            m.low = low;
            m.high = high;
            m.anchor = Some("ancestry-fine".into());
            m.group = Some(title.to_owned());
            marks.push(m);
        }
    }
    if marks.is_empty() {
        return Ok(None);
    }
    let mut sp = spec(
        "ancestry-fine",
        "ancestry_fine",
        Kind::Ranges,
        "Fine-grained ancestry, with bootstrap ranges",
        "ANCESTRY_FINE_NOTE",
        marks,
    )?;
    sp.axis = axis(&[
        ("min", Value::Int(0)),
        ("max", Value::Int(100)),
        ("unit", Value::Str("%".into())),
    ]);
    Ok(Some(sp))
}

fn prs_chart(s: &Summary) -> Result<Option<Spec>, Malformed> {
    let Some(d) = ran(s, "prs") else {
        return Ok(None);
    };
    let path = "sections.prs.data";
    if !get(d, "adjusted", path)?.is_some_and(truthy) {
        let mut sp = spec(
            "prs",
            "prs",
            Kind::Text,
            "Polygenic risk scores",
            "PRS_ANCESTRY_NOTE",
            Vec::new(),
        )?;
        sp.text = Some(prs::note(d)?);
        return Ok(Some(sp));
    }
    let at = format!("{path}.scores[]");
    let mut marks = Vec::new();
    for (i, r) in list(d, "scores", path)?.iter().enumerate() {
        let other =
            matches!(get(r, "sex_check", &at)?, Some(Value::Str(x)) if x == prs::NOT_APPLICABLE);
        let mut value = None;
        if !other && let Some(p) = get(r, "percentile", &at)?.filter(|v| truthy(v)) {
            value = py_float(p, &format!("{at}.percentile"))?.filter(|v| (0.0..=100.0).contains(v));
        }
        let id = s_truthy(r, "pgs_id", &i.to_string(), &at)?;
        let label = format!(
            "{} ({})",
            s_truthy(r, "condition", ".", &at)?,
            s_truthy(r, "pgs_id", ".", &at)?
        );
        let text = prs::percentile(r)?;
        let mut m = Mark::new(
            format!("pgs-{id}"),
            label,
            value.map_or(Value::Null, Value::Float),
        );
        m.tone = if value.is_some() {
            Tone::Info
        } else {
            Tone::Muted
        };
        m.anchor = Some("prs".into());
        m.detail = match (value, text.is_empty()) {
            (Some(_), false) => text,
            (Some(_), true) => "no percentile".into(),
            (None, false) => format!("no percentile: {text}"),
            (None, true) => "no percentile".into(),
        };
        marks.push(m);
    }
    let mut sp = spec(
        "prs",
        "prs",
        Kind::Strip,
        "Polygenic risk score percentiles",
        "PRS_ANCESTRY_NOTE",
        marks,
    )?;
    sp.axis = axis(&[
        ("min", Value::Int(0)),
        ("max", Value::Int(100)),
        ("unit", Value::Str("percentile".into())),
    ]);
    Ok(Some(sp))
}

/// The genome chart sits in the ROH section, so it is drawn only when ROH
/// ran; the ClinVar ticks, AnnotSV spans and MoChA bars come with it when
/// their steps ran too.
fn genome(s: &Summary) -> Result<Option<Spec>, Malformed> {
    let Some(roh) = ran(s, "roh") else {
        return Ok(None);
    };
    let mut marks = Vec::new();
    let mut not_drawn = Tally::default();
    let at = "sections.roh.data.autosomal_over_5mb[]";
    let mut short = 0;
    for (i, r) in list(roh, "autosomal_over_5mb", "sections.roh.data")?
        .iter()
        .enumerate()
    {
        let reg = s_or(r, "region", "", at)?;
        let shown = || -> Result<String, Malformed> {
            match get(r, "region", at)? {
                Some(v) => py_repr(v).map_err(|e| Malformed(e.to_string())),
                None => Ok("None".into()),
            }
        };
        let Some((chrom, start, end)) = region(&reg) else {
            return Err(Malformed(format!(
                "run of homozygosity {}: not a chrom:start-end region",
                shown()?
            )));
        };
        if end < start {
            return Err(Malformed(format!(
                "run of homozygosity {}: the end is before the start",
                shown()?
            )));
        }
        let mb = match get(r, "mb", at)? {
            Some(v) => need_float(v, &format!("{at}.mb"))?,
            None => 0.0,
        };
        if mb < ROH_MIN_MB {
            short += 1;
            continue;
        }
        if chrom_len(&chrom).is_none() {
            not_drawn.add(&format!("runs of homozygosity on {chrom}"));
            continue;
        }
        let mb_text = py(need(r, "mb", at)?)?;
        let mut m = Mark::new(
            format!("roh-{i}"),
            format!("Run of homozygosity {reg}"),
            Value::Float(mb),
        );
        m.locus = Some(format!("{chrom}:{start}-{end}"));
        m.chrom = Some(chrom);
        m.pos = Some(start);
        m.end = Some(end);
        m.anchor = Some("roh".into());
        m.kind = Some(MarkKind::Roh);
        m.detail = format!("{mb_text} MB run of homozygosity at {reg}");
        marks.push(m);
    }
    if short > 0 {
        not_drawn
            .0
            .insert("runs of homozygosity under 5 MB".into(), short);
    }
    if let Some(c) = ran(s, "clinvar") {
        let hp = "sections.clinvar.data.hits[]";
        for (i, h) in list(c, "hits", "sections.clinvar.data")?.iter().enumerate() {
            let chrom = match get(h, "chrom", hp)? {
                Some(Value::Str(x)) => Some(x.as_str()),
                _ => None,
            };
            if chrom == Some("chrM") {
                not_drawn.add("ClinVar hits on chrM");
                continue;
            }
            let Some(chrom) = chrom.filter(|c| chrom_len(c).is_some()) else {
                not_drawn.add("ClinVar hits on other contigs (alt, unplaced)");
                continue;
            };
            let low = get(h, "low_confidence", hp)?.is_some_and(truthy);
            let pos = need(h, "pos", hp)?;
            let stars = raw(h, "stars", hp)?;
            let one = if py_eq(&stars, &Value::Int(1)) {
                ""
            } else {
                "s"
            };
            let mut m = Mark::new(format!("hit-{i}"), hit_label(h, hp)?, stars.clone());
            m.chrom = Some(chrom.to_owned());
            m.pos = Some(py_int(pos, &format!("{hp}.pos"))?);
            m.tone = if low { Tone::Muted } else { Tone::Attention };
            m.anchor = Some(hit_anchor(i).to_owned());
            m.locus = Some(format!("{chrom}:{}", py(pos)?));
            m.kind = Some(MarkKind::Hit);
            m.caveat = Some(caveat("CLINVAR_FREQ_NOTE")?.to_owned());
            m.detail = format!(
                "{}, {} star{one}, {}{}",
                s_truthy(h, "significance", ".", hp)?,
                py(&stars)?,
                s_truthy(h, "genotype", ".", hp)?,
                if low { ", low-confidence call" } else { "" }
            );
            marks.push(m);
        }
    }
    if let Some(a) = ran(s, "annotsv") {
        let rp = "sections.annotsv.data.rows[]";
        for (i, r) in list(a, "rows", "sections.annotsv.data")?.iter().enumerate() {
            let sv = s_or(r, "sv", "", rp)?;
            // AnnotSV writes SV_chrom without the chr prefix: put it back, as the collector does.
            let full = if sv.starts_with("chr") {
                sv.clone()
            } else {
                format!("chr{sv}")
            };
            let Some((chrom, start, end)) =
                region(&full).filter(|(c, _, _)| chrom_len(c).is_some())
            else {
                not_drawn.add("AnnotSV class 4 and 5 SVs on other contigs");
                continue;
            };
            let genes = match get(r, "genes", rp)?.filter(|v| truthy(v)) {
                None => String::new(),
                Some(Value::Array(g)) => g
                    .iter()
                    .take(5)
                    .map(py)
                    .collect::<Result<Vec<_>, _>>()?
                    .join(", "),
                Some(Value::Str(g)) => g
                    .chars()
                    .take(5)
                    .map(String::from)
                    .collect::<Vec<_>>()
                    .join(", "),
                Some(_) => return Err(bad(&format!("{rp}.genes"), "not a list")),
            };
            let genes = if genes.is_empty() {
                "no gene".to_owned()
            } else {
                genes
            };
            let class = raw(r, "acmg_class", rp)?;
            let kind = s_truthy(r, "type", "SV", rp)?;
            let mut m = Mark::new(
                format!("sv-{i}"),
                format!(
                    "{kind} {} (AnnotSV class {})",
                    py(need(r, "sv", rp)?)?,
                    py(&class)?
                ),
                class.clone(),
            );
            m.locus = Some(format!("{chrom}:{start}-{end}"));
            m.chrom = Some(chrom);
            m.pos = Some(start);
            m.end = Some(end.max(start));
            m.tone = Tone::Attention;
            m.anchor = Some("structural-variants".into());
            m.kind = Some(MarkKind::Sv);
            m.caveat = Some(caveat("ANNOTSV_NOTE")?.to_owned());
            m.detail = format!(
                "AnnotSV class {} {kind}, {} bp, genotype {}, genes: {genes}",
                py(&class)?,
                s_truthy(r, "length", ".", rp)?,
                s_truthy(r, "gt", ".", rp)?
            );
            marks.push(m);
        }
    }
    if let Some(mo) = ran(s, "mosaic") {
        let cp = "sections.mosaic.data.calls[]";
        for (i, c) in list(mo, "calls", "sections.mosaic.data")?
            .iter()
            .enumerate()
        {
            let mut chrom = s_truthy(c, "chrom", "", cp)?;
            if !chrom.is_empty() && !chrom.starts_with("chr") {
                chrom = format!("chr{chrom}");
            }
            // int() of start and end; a value that is not one leaves the call unplaced.
            let as_int = |k: &str| -> Option<i128> {
                let v = get(c, k, cp).ok()??;
                py_int(v, k).ok()
            };
            let span = as_int("start").zip(as_int("end"));
            let Some((start, end)) = span.filter(|(a, b)| chrom_len(&chrom).is_some() && b >= a)
            else {
                not_drawn.add("MoChA calls on other contigs or with no position");
                continue;
            };
            let fraction = raw(c, "cell_fraction", cp)?;
            let mut m = Mark::new(
                format!("mca-{i}"),
                format!(
                    "MoChA {} {chrom}:{start}-{end}",
                    s_truthy(c, "type", "call", cp)?
                ),
                fraction,
            );
            m.locus = Some(format!("{chrom}:{start}-{end}"));
            m.chrom = Some(chrom);
            m.pos = Some(start);
            m.end = Some(end);
            m.tone = Tone::Attention;
            m.anchor = Some("mosaic-mca".into());
            m.kind = Some(MarkKind::Mca);
            m.caveat = Some(caveat("MCA_NOTE")?.to_owned());
            m.detail = format!(
                "mosaic {} over {} Mb, cell fraction {}",
                s_truthy(c, "type", "alteration", cp)?,
                s_truthy(c, "length_mb", ".", cp)?,
                s_truthy(c, "cell_fraction", ".", cp)?
            );
            marks.push(m);
        }
    }
    // A position past the end of its chromosome (a file of another assembly)
    // is not drawn; an end past it is cut at the end, so every mark stays
    // inside its chromosome.
    let mut kept = Vec::new();
    for mut m in marks {
        let n = m.chrom.as_deref().and_then(chrom_len).unwrap_or(0);
        if m.pos.is_some_and(|p| p > n) {
            not_drawn.add("marks past the end of their GRCh38 chromosome");
            continue;
        }
        if m.end.is_some_and(|e| e > n) {
            m.end = Some(n);
        }
        kept.push(m);
    }
    let mut sp = spec(
        "genome",
        "roh",
        Kind::Ideogram,
        "Runs of homozygosity, ClinVar hits, structural variants and mosaic alterations on the \
         chromosomes",
        "ROH_NOTE",
        kept,
    )?;
    sp.axis = axis(&[("assembly", Value::Str("GRCh38".into()))]);
    sp.not_drawn = not_drawn.list();
    Ok(Some(sp))
}

/// `^chrM:([0-9]+) `: the position of a poly-C call.
fn polyc_pos(variant: &str) -> Option<i128> {
    let rest = variant.strip_prefix("chrM:")?;
    let n = rest.bytes().take_while(u8::is_ascii_digit).count();
    (n > 0 && rest.as_bytes().get(n) == Some(&b' ')).then(|| rest[..n].parse().ok())?
}

fn mito(s: &Summary) -> Result<Option<Spec>, Malformed> {
    let Some(d) = ran(s, "mito") else {
        return Ok(None);
    };
    let path = "sections.mito.data";
    let sp_at = format!("{path}.heteroplasmic_sites[]");
    let mut marks = Vec::new();
    let mut not_drawn = Vec::new();
    let sites = get(d, "heteroplasmic_sites", path)?.filter(|v| !matches!(v, Value::Null));
    let site_list: &[Value] = match sites {
        None => &[],
        Some(v) if !truthy(v) => &[],
        Some(Value::Array(a)) => a,
        Some(_) => return Err(bad(&format!("{path}.heteroplasmic_sites"), "not a list")),
    };
    for x in site_list {
        let pos = need(x, "pos", &sp_at)?;
        let (p, alt, af) = (
            py(pos)?,
            py(need(x, "alt", &sp_at)?)?,
            need(x, "af", &sp_at)?,
        );
        let mut m = Mark::new(
            format!("m{p}-{alt}"),
            format!("m.{p}{}>{alt}", py(need(x, "ref", &sp_at)?)?),
            Value::Float(need_float(af, &format!("{sp_at}.af"))?),
        );
        m.chrom = Some("chrM".into());
        m.pos = Some(py_int(pos, &format!("{sp_at}.pos"))?);
        m.tone = Tone::Attention;
        m.anchor = Some("mito".into());
        m.locus = Some(format!("chrM:{p}"));
        m.kind = Some(MarkKind::Site);
        m.detail = format!("allele fraction {}", py(af)?);
        marks.push(m);
    }
    let missing = match (sites, get(d, "heteroplasmic", path)?.filter(|v| truthy(v))) {
        (None, Some(h)) => py_int(h, &format!("{path}.heteroplasmic"))?,
        _ => 0,
    };
    if missing != 0 {
        not_drawn.push(NotDrawn {
            what: "heteroplasmic sites this summary does not list (made before it listed them)"
                .into(),
            count: missing,
        });
    }
    let pc_at = format!("{path}.polyc_heteroplasmic[]");
    for x in list(d, "polyc_heteroplasmic", path)? {
        let variant = match get(x, "variant", &pc_at)? {
            None => String::new(),
            Some(Value::Str(v)) => v.clone(),
            Some(_) => return Err(bad(&format!("{pc_at}.variant"), "not a text")),
        };
        let Some(pos) = polyc_pos(&variant) else {
            let shown = match get(x, "variant", &pc_at)? {
                Some(v) => py_repr(v).map_err(|e| Malformed(e.to_string()))?,
                None => "None".into(),
            };
            not_drawn.push(NotDrawn {
                what: format!("poly-C call {shown} with no position"),
                count: 1,
            });
            continue;
        };
        let af = need(x, "af", &pc_at)?;
        let mut m = Mark::new(
            format!("polyc-{pos}"),
            format!("{variant} (D-loop poly-C tract)"),
            Value::Float(need_float(af, &format!("{pc_at}.af"))?),
        );
        m.chrom = Some("chrM".into());
        m.pos = Some(pos);
        m.tone = Tone::Muted;
        m.anchor = Some("mito".into());
        m.locus = Some(format!("chrM:{pos}"));
        m.kind = Some(MarkKind::Polyc);
        m.detail = format!(
            "length heteroplasmy in a D-loop poly-C tract, allele fraction {}; not counted as heteroplasmic",
            py(af)?
        );
        marks.push(m);
    }
    let floor = get(d, "heteroplasmy_floor", path)?
        .cloned()
        .unwrap_or(Value::Float(0.05));
    let mut sp = spec(
        "mito",
        "mito",
        Kind::Scatter,
        "Heteroplasmic chrM sites",
        "MITO_NOTE",
        marks,
    )?;
    sp.axis = axis(&[
        ("min", Value::Int(0)),
        ("max", Value::Int(1)),
        ("unit", Value::Str("allele fraction".into())),
        ("floor", floor),
        ("sites_missing", Value::Int(missing)),
    ]);
    sp.not_drawn = not_drawn;
    Ok(Some(sp))
}
