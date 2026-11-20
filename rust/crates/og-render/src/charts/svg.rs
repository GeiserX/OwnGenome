//! Each chart as inline SVG, the same bytes `bin/report_charts.py` drew: a
//! coordinate with one decimal at most, the arithmetic in Python's order, and
//! `min` and `max` with Python's rules (a NaN passed second never wins).
//!
//! Text uses `currentColor`, so it follows the page's colour on screen and
//! in print.

use pycompat::fmt::{fixed, general};
use pycompat::json::Value;
use pycompat::text::html_escape;

use super::{GRCH38, Kind, MITO_LEN, Mark, MarkKind, Spec, Tone, chrom_len};
use crate::model::{Malformed, num, py};

/// The SVG's width; its height follows the content.
const W: f64 = 640.0;
const FRAME: &str = "#8a909b";
const LINE: &str = "#4a505b";

fn e(s: &str) -> String {
    html_escape(s).into_owned()
}

/// Python's `min(a, b)`: `b` only when `b < a`.
fn pmin(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// Python's `max(a, b)`: `b` only when `b > a`.
fn pmax(a: f64, b: f64) -> f64 {
    if b > a { b } else { a }
}

/// A coordinate: one decimal, trailing zeros and point dropped.
fn f(x: f64) -> Result<String, Malformed> {
    if !x.is_finite() {
        return Err(Malformed(format!(
            "a chart coordinate is {}, not a number",
            pycompat::fmt::repr(x)
        )));
    }
    let s = fixed(x, 1);
    Ok(s.trim_end_matches('0').trim_end_matches('.').to_owned())
}

struct Text<'a> {
    size: u32,
    anchor: &'a str,
    extra: &'a str,
}

const PLAIN: Text<'static> = Text {
    size: 12,
    anchor: "start",
    extra: "",
};

fn text(x: f64, y: f64, s: &str, t: &Text<'_>) -> Result<String, Malformed> {
    Ok(format!(
        "<text x=\"{}\" y=\"{}\" font-size=\"{}\" text-anchor=\"{}\" fill=\"currentColor\"{}>{}</text>",
        f(x)?,
        f(y)?,
        t.size,
        t.anchor,
        t.extra,
        e(s)
    ))
}

fn sized(size: u32, anchor: &str) -> Text<'_> {
    Text {
        size,
        anchor,
        extra: "",
    }
}

/// A mark: a link to its anchor in the report, with its title and id.
fn a(m: &Mark, inner: &str) -> String {
    let title = if m.detail.is_empty() {
        m.label.clone()
    } else {
        format!("{}: {}", m.label, m.detail)
    };
    format!(
        "<a href=\"#{}\" data-mark=\"{}\"><title>{}</title>{inner}</a>",
        e(m.anchor.as_deref().unwrap_or("")),
        e(&m.id),
        e(&title)
    )
}

fn open(sp: &Spec, h: f64) -> Result<String, Malformed> {
    Ok(format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" class=\"chart\" data-chart=\"{}\" viewBox=\"0 0 {W} {}\" \
         width=\"100%\" role=\"img\" aria-label=\"{}\"><title>{}</title>",
        e(sp.id),
        f(h)?,
        e(sp.title),
        e(sp.title)
    ))
}

/// Words of `s` in lines of at most `n` characters (a longer word alone).
fn wrap(s: &str, n: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut cur = String::new();
    for w in pycompat::text::split(s, None, None) {
        if !cur.is_empty() && cur.chars().count() + 1 + w.chars().count() > n {
            lines.push(std::mem::replace(&mut cur, w.to_owned()));
        } else {
            cur = if cur.is_empty() {
                w.to_owned()
            } else {
                format!("{cur} {w}")
            };
        }
    }
    if !cur.is_empty() {
        lines.push(cur);
    }
    lines
}

/// The "Not drawn:" lines under a chart, from `y`, and the y after them.
fn footer(sp: &Spec, y: f64) -> Result<(String, f64), Malformed> {
    if sp.not_drawn.is_empty() {
        return Ok((String::new(), y));
    }
    let parts: Vec<String> = sp
        .not_drawn
        .iter()
        .map(|x| format!("{} {}", x.count, x.what))
        .collect();
    let lines = wrap(&format!("Not drawn: {}.", parts.join("; ")), 100);
    let t = Text {
        size: 11,
        anchor: "start",
        extra: " opacity=\"0.75\"",
    };
    let mut out = String::new();
    for (i, s) in lines.iter().enumerate() {
        out += &text(8.0, y + 14.0 * i as f64, s, &t)?;
    }
    Ok((out, y + 14.0 * lines.len() as f64))
}

/// A mark's value as a number, 0 for `null` (Python's `value or 0`).
fn value(m: &Mark) -> Result<f64, Malformed> {
    match &m.value {
        Value::Null => Ok(0.0),
        v => num(Some(v))
            .ok_or_else(|| Malformed(format!("mark {}: the value is not a number", m.id))),
    }
}

/// A label cut to 40 characters.
fn short(label: &str) -> String {
    if label.chars().count() <= 40 {
        label.to_owned()
    } else {
        label.chars().take(38).collect::<String>() + "..."
    }
}

fn bars(sp: &Spec) -> Result<String, Malformed> {
    let vals: Vec<f64> = sp.marks.iter().map(value).collect::<Result<_, _>>()?;
    if vals.iter().sum::<f64>() == 0.0 {
        let h = 96.0;
        return Ok(open(sp, h)?
            + &format!(
                "<rect x=\"1\" y=\"1\" width=\"{}\" height=\"{}\" rx=\"6\" fill=\"none\" stroke=\"{FRAME}\" \
                 stroke-dasharray=\"6 4\"/>",
                W - 2.0,
                h - 2.0
            )
            + &text(
                W / 2.0,
                40.0,
                "No ClinVar hit in this sample.",
                &sized(14, "middle"),
            )?
            + &text(
                W / 2.0,
                64.0,
                "No hit is not an all-clear: the screen reads only the variants ClinVar lists.",
                &sized(12, "middle"),
            )?
            + "</svg>");
    }
    let top = vals.iter().copied().fold(vals[0], pmax);
    if top == 0.0 {
        return Err(Malformed("ClinVar chart: no bar above 0".into()));
    }
    let (left, right, row) = (190.0, 60.0, 26.0);
    let h = 16.0 + row * sp.marks.len() as f64 + 8.0;
    let mut body = String::new();
    for (i, (m, v)) in sp.marks.iter().zip(&vals).enumerate() {
        let y = 12.0 + i as f64 * row;
        let wbar = (W - left - right) * v / top;
        let bar = if *v != 0.0 {
            format!(
                "<rect x=\"{left}\" y=\"{}\" width=\"{}\" height=\"18\" rx=\"3\" fill=\"{}\"/>",
                f(y)?,
                f(pmax(wbar, 0.0))?,
                m.tone.colour()
            )
        } else {
            format!(
                "<rect x=\"{left}\" y=\"{}\" width=\"2\" height=\"18\" fill=\"{FRAME}\"/>",
                f(y)?
            )
        };
        let inner = text(left - 8.0, y + 13.0, &m.label, &sized(12, "end"))?
            + &bar
            + &text(
                left + pmax(wbar, 2.0) + 6.0,
                y + 13.0,
                &py(&m.value)?,
                &PLAIN,
            )?;
        body += &a(m, &inner);
    }
    let (foot, h2) = footer(sp, h + 4.0)?;
    Ok(open(sp, h2 + 4.0)? + &body + &foot + "</svg>")
}

/// The ticks of a 0 to 100 axis, at `x(t)`.
fn ticks(x: impl Fn(f64) -> f64, label: impl Fn(u32) -> String) -> Result<String, Malformed> {
    let mut out = String::new();
    for t in [0, 25, 50, 75, 100] {
        out += &text(x(f64::from(t)), 10.0, &label(t), &sized(10, "middle"))?;
    }
    Ok(out)
}

fn ranges(sp: &Spec) -> Result<String, Malformed> {
    let (left, right, row) = (260.0, 50.0, 22.0);
    let (x0, x1) = (left, W - right);
    let sx = |v: f64| x0 + (x1 - x0) * pmax(0.0, pmin(100.0, v)) / 100.0;
    let info = Tone::Info.colour();
    let mut y = 22.0;
    let mut body = String::new();
    let mut group: Option<&str> = None;
    let axis = format!("<line x1=\"{x0}\" y1=\"14\" x2=\"{x1}\" y2=\"14\" stroke=\"{LINE}\"/>")
        + &ticks(sx, |t| format!("{t}%"))?;
    for m in &sp.marks {
        if m.group.as_deref() != group {
            group = m.group.as_deref();
            y += 6.0;
            let bold = Text {
                size: 12,
                anchor: "start",
                extra: " font-weight=\"600\"",
            };
            body += &text(8.0, y + 10.0, group.unwrap_or("None"), &bold)?;
            y += row;
        }
        let cy = y + 6.0;
        let v = value(m)?;
        let mut inner = text(left - 8.0, cy + 4.0, &short(&m.label), &sized(12, "end"))?;
        if let (Some(lo), Some(hi)) = (num(Some(&m.low)), num(Some(&m.high))) {
            inner += &format!(
                "<line class=\"whisker\" x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{info}\" stroke-width=\"2\"/>",
                f(sx(lo))?,
                f(cy)?,
                f(sx(hi))?,
                f(cy)?
            );
        }
        inner += &format!(
            "<circle cx=\"{}\" cy=\"{}\" r=\"4.5\" fill=\"{info}\"/>",
            f(sx(v))?,
            f(cy)?
        );
        inner += &text(
            x1 + 6.0,
            cy + 4.0,
            &format!("{}%", fixed(v, 0)),
            &sized(11, "start"),
        )?;
        body += &a(m, &inner);
        y += row;
    }
    Ok(open(sp, y + 8.0)? + &axis + &body + "</svg>")
}

fn strip(sp: &Spec) -> Result<String, Malformed> {
    let (left, right, row) = (260.0, 120.0, 24.0);
    let (x0, x1) = (left, W - right);
    let sx = |v: f64| x0 + (x1 - x0) * v / 100.0;
    let mut body = format!("<line x1=\"{x0}\" y1=\"14\" x2=\"{x1}\" y2=\"14\" stroke=\"{LINE}\"/>")
        + &ticks(sx, |t| t.to_string())?;
    let mut y = 26.0;
    for m in &sp.marks {
        let cy = y + 6.0;
        let mut inner = text(left - 8.0, cy + 4.0, &short(&m.label), &sized(12, "end"))?;
        inner += &format!(
            "<line x1=\"{x0}\" y1=\"{}\" x2=\"{x1}\" y2=\"{}\" stroke=\"{LINE}\"/>",
            f(cy)?,
            f(cy)?
        );
        if matches!(m.value, Value::Null) {
            let faint = Text {
                size: 11,
                anchor: "start",
                extra: " opacity=\"0.8\"",
            };
            inner += &text(x0 + 4.0, cy + 4.0, "no percentile", &faint)?;
        } else {
            let v = value(m)?;
            inner += &format!(
                "<circle cx=\"{}\" cy=\"{}\" r=\"5\" fill=\"{}\"/>",
                f(sx(v))?,
                f(cy)?,
                m.tone.colour()
            );
            inner += &text(x1 + 8.0, cy + 4.0, &general(v, 6), &sized(11, "start"))?;
        }
        body += &a(m, &inner);
        y += row;
    }
    Ok(open(sp, y + 8.0)? + &body + "</svg>")
}

fn text_chart(sp: &Spec) -> Result<String, Malformed> {
    let lines = wrap(sp.text.as_deref().unwrap_or(""), 95);
    let h = 20.0 + 16.0 * lines.len() as f64;
    let mut out = open(sp, h)?
        + &format!(
            "<rect x=\"1\" y=\"1\" width=\"{}\" height=\"{}\" rx=\"6\" fill=\"none\" stroke=\"{FRAME}\"/>",
            W - 2.0,
            f(h - 2.0)?
        );
    for (i, s) in lines.iter().enumerate() {
        out += &text(12.0, 20.0 + 16.0 * i as f64, s, &PLAIN)?;
    }
    Ok(out + "</svg>")
}

fn ideogram(sp: &Spec) -> Result<String, Malformed> {
    let (left, right, row) = (46.0, 16.0, 18.0);
    let (x0, x1) = (left, W - right);
    let top = GRCH38.iter().map(|(_, n)| *n).max().unwrap_or(1) as f64;
    let sx = |bp: f64| x0 + (x1 - x0) * bp / top;
    let row_y = |c: &str| {
        GRCH38
            .iter()
            .position(|(n, _)| *n == c)
            .map(|i| 10 + i as i64 * 18)
    };
    let (attention, info) = (Tone::Attention.colour(), Tone::Info.colour());
    let mut body = String::new();
    for (i, (c, n)) in GRCH38.iter().enumerate() {
        let y = 10 + i as i64 * row as i64;
        body += &text(left - 6.0, (y + 9) as f64, &c[3..], &sized(10, "end"))?;
        body += &format!(
            "<rect x=\"{x0}\" y=\"{}\" width=\"{}\" height=\"8\" rx=\"4\" fill=\"none\" stroke=\"{LINE}\"/>",
            y + 2,
            f(sx(*n as f64) - x0)?
        );
    }
    for m in &sp.marks {
        let chrom = m.chrom.as_deref().unwrap_or("");
        let (Some(y), Some(len)) = (row_y(chrom), chrom_len(chrom)) else {
            return Err(Malformed(format!(
                "mark {}: {chrom} is not a GRCh38 chromosome",
                m.id
            )));
        };
        let (lo, hi) = (x0, sx(len as f64));
        // The x of a mark w px wide, moved inside its chromosome's outline.
        let inside = |a: f64, w: f64| pmax(lo, pmin(a, hi - w));
        let pos = m.pos.unwrap_or(0) as f64;
        let span = || {
            let a = sx(pos);
            let b = sx(m.end.unwrap_or(0) as f64);
            (a, b)
        };
        let inner = match m.kind {
            Some(MarkKind::Roh) => {
                let (a, b) = span();
                let w = pmax(b - a, 1.5);
                format!(
                    "<rect class=\"roh\" x=\"{}\" y=\"{}\" width=\"{}\" height=\"8\" fill=\"{info}\"/>",
                    f(inside(a, w))?,
                    y + 2,
                    f(w)?
                )
            }
            Some(MarkKind::Sv) => {
                // An outlined span over the chromosome, at least 5 px so a small deletion is seen.
                let (a, b) = span();
                let w = pmax(b - a, 5.0);
                format!(
                    "<rect class=\"sv\" x=\"{}\" y=\"{}\" width=\"{}\" height=\"14\" fill=\"none\" \
                     stroke=\"{attention}\" stroke-width=\"2\"/>",
                    f(inside(a - (w - (b - a)) / 2.0, w))?,
                    y - 1,
                    f(w)?
                )
            }
            Some(MarkKind::Mca) => {
                // A bar under the chromosome: a share of the cells, not the whole genome.
                let (a, b) = span();
                let w = pmax(b - a, 3.0);
                format!(
                    "<rect class=\"mca\" x=\"{}\" y=\"{}\" width=\"{}\" height=\"3\" fill=\"{attention}\"/>",
                    f(inside(a, w))?,
                    y + 12,
                    f(w)?
                )
            }
            _ => format!(
                "<rect class=\"tick\" x=\"{}\" y=\"{}\" width=\"3\" height=\"14\" fill=\"{}\"/>",
                f(inside(sx(pos) - 1.5, 3.0))?,
                y - 1,
                m.tone.colour()
            ),
        };
        body += &a(m, &inner);
    }
    let mut y = 10 + GRCH38.len() as i64 * row as i64 + 6;
    let small = sized(11, "start");
    let mut legend =
        format!("<rect x=\"{x0}\" y=\"{y}\" width=\"14\" height=\"8\" fill=\"{info}\"/>")
            + &text(
                x0 + 20.0,
                (y + 8) as f64,
                "run of homozygosity, 5 MB or more",
                &small,
            )?
            + &format!(
                "<rect x=\"{}\" y=\"{}\" width=\"3\" height=\"12\" fill=\"{attention}\"/>",
                x0 + 250.0,
                y - 2
            )
            + &text(x0 + 258.0, (y + 8) as f64, "ClinVar hit", &small)?
            + &format!(
                "<rect x=\"{}\" y=\"{}\" width=\"3\" height=\"12\" fill=\"{}\"/>",
                x0 + 340.0,
                y - 2,
                Tone::Muted.colour()
            )
            + &text(
                x0 + 348.0,
                (y + 8) as f64,
                "low-confidence ClinVar call",
                &small,
            )?;
    let has = |k: MarkKind| sp.marks.iter().any(|m| m.kind == Some(k));
    if has(MarkKind::Sv) || has(MarkKind::Mca) {
        // A second legend row, drawn only when the chart has these marks.
        y += 18;
        if has(MarkKind::Sv) {
            legend += &format!(
                "<rect x=\"{x0}\" y=\"{}\" width=\"14\" height=\"12\" fill=\"none\" stroke=\"{attention}\" \
                 stroke-width=\"2\"/>",
                y - 2
            );
            legend += &text(x0 + 20.0, (y + 8) as f64, "AnnotSV class 4 or 5 SV", &small)?;
        }
        if has(MarkKind::Mca) {
            legend += &format!(
                "<rect x=\"{}\" y=\"{}\" width=\"14\" height=\"3\" fill=\"{attention}\"/>",
                x0 + 250.0,
                y + 4
            );
            legend += &text(
                x0 + 270.0,
                (y + 8) as f64,
                "mosaic alteration (MoChA), a share of cells",
                &small,
            )?;
        }
    }
    let (foot, h) = footer(sp, (y + 26) as f64)?;
    Ok(open(sp, h + 4.0)? + &body + &legend + &foot + "</svg>")
}

fn scatter(sp: &Spec) -> Result<String, Malformed> {
    let (left, right, tp, bottom) = (44.0, 16.0, 12.0, 30.0);
    let h = 200.0;
    let (x0, x1, y0, y1) = (left, W - right, tp, h - bottom);
    let sx = |p: f64| x0 + (x1 - x0) * (p - 1.0) / (MITO_LEN - 1) as f64;
    let sy = |af: f64| y1 - (y1 - y0) * pmax(0.0, pmin(1.0, af));
    let muted = Tone::Muted.colour();
    let mut body = format!(
        "<rect x=\"{x0}\" y=\"{y0}\" width=\"{}\" height=\"{}\" fill=\"none\" stroke=\"{LINE}\"/>",
        x1 - x0,
        y1 - y0
    );
    for t in [0.0, 0.5, 1.0] {
        body += &text(x0 - 6.0, sy(t) + 4.0, &general(t, 6), &sized(10, "end"))?;
    }
    for (t, anchor) in [
        (1, "start"),
        (4000, "middle"),
        (8000, "middle"),
        (12000, "middle"),
        (MITO_LEN, "end"),
    ] {
        body += &text(
            sx(t as f64),
            y1 + 14.0,
            &format!("m.{t}"),
            &sized(10, anchor),
        )?;
    }
    let faint = Text {
        size: 10,
        anchor: "middle",
        extra: " opacity=\"0.8\"",
    };
    body += &text(
        (x0 + x1) / 2.0,
        h - 2.0,
        "position on chrM; height: allele fraction",
        &faint,
    )?;
    for (lo, hi) in [(303.0, 315.0), (16184.0, 16193.0)] {
        body += &format!(
            "<rect x=\"{}\" y=\"{y0}\" width=\"{}\" height=\"{}\" fill=\"{muted}\" opacity=\"0.25\"/>",
            f(sx(lo) - 1.0)?,
            f(pmax(sx(hi) - sx(lo), 2.0))?,
            y1 - y0
        );
    }
    let missing = sp
        .axis
        .get("sites_missing")
        .is_some_and(crate::model::truthy);
    let mid = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    if missing {
        body += &text(
            mid.0,
            mid.1,
            "The positions of the heteroplasmic sites are not in this summary",
            &sized(12, "middle"),
        )?;
    } else if sp.marks.is_empty() {
        body += &text(
            mid.0,
            mid.1,
            "No heteroplasmic site between allele fraction 0.05 and 0.95",
            &sized(12, "middle"),
        )?;
    }
    for m in &sp.marks {
        let (x, y) = (sx(m.pos.unwrap_or(0) as f64), sy(value(m)?));
        let inner = if m.kind == Some(MarkKind::Polyc) {
            format!(
                "<rect class=\"polyc\" x=\"{}\" y=\"{}\" width=\"8\" height=\"8\" fill=\"none\" \
                 stroke=\"{muted}\" stroke-width=\"2\"/>",
                f(x - 4.0)?,
                f(y - 4.0)?
            )
        } else {
            format!(
                "<circle class=\"site\" cx=\"{}\" cy=\"{}\" r=\"5\" fill=\"{}\"/>",
                f(x)?,
                f(y)?,
                Tone::Attention.colour()
            )
        };
        body += &a(m, &inner);
    }
    let (foot, h2) = footer(sp, h + 10.0)?;
    Ok(open(sp, h2 + 4.0)? + &body + &foot + "</svg>")
}

/// The chart drawn as SVG, by its kind.
pub fn draw(sp: &Spec) -> Result<String, Malformed> {
    match sp.kind {
        Kind::Bars => bars(sp),
        Kind::Ranges => ranges(sp),
        Kind::Strip => strip(sp),
        Kind::Text => text_chart(sp),
        Kind::Ideogram => ideogram(sp),
        Kind::Scatter => scatter(sp),
    }
}
