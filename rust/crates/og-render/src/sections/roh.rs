//! Step 11's runs of homozygosity: the segments, the long runs, F_ROH with
//! the nearest relationship band, and a run concentrated on one chromosome
//! (a possible uniparental isodisomy) with the imprinted regions it covers.

use pycompat::json::Value;

use crate::model::{
    Card, Colour, HtmlPart, Inline, Line, Malformed, Section, State, Summary, get, need, num, py,
    truthy,
};
use crate::notes::{self, ROH_BAND_NOTE, ROH_F_NOTE, ROH_NOTE, ROH_UPD_NOTE};
use crate::{html, text};

const PATH: &str = "sections.roh.data";

/// One long run: its region and length, and the part left out in a band.
pub fn region(r: &Value, path: &str) -> Result<String, Malformed> {
    let mut out = format!(
        "{}  {}MB",
        py(need(r, "region", path)?)?,
        py(need(r, "mb", path)?)?
    );
    if let Some(m) = get(r, "masked_mb", path)?.filter(|v| truthy(v)) {
        out += &format!(" ({} MB in a band left out)", py(m)?);
    }
    Ok(out)
}

fn list<'a>(v: Option<&'a Value>, path: &str) -> Result<&'a [Value], Malformed> {
    match v {
        None => Ok(&[]),
        Some(Value::Array(a)) => Ok(a),
        Some(_) => Err(Malformed(format!("{path} is not a list"))),
    }
}

/// One line of [`lines`]: a sentence, or caveats joined by spaces, which
/// the HTML card shows under its title.
enum Item {
    Text(String),
    Caveats(Vec<&'static str>),
}

fn items(d: &Section<'_>) -> Result<Vec<Item>, Malformed> {
    let mut out = Vec::new();
    let f_roh = d.get("f_roh")?.filter(|v| num(Some(v)).is_some());
    if let Some(f) = f_roh {
        let band = match d.get("relationship")?.filter(|v| truthy(v)) {
            Some(b) => py(b)?,
            None => "no band (one chromosome holds most of the runs)".to_owned(),
        };
        let long = d.get("long_mb")?.map_or(Ok("None".to_owned()), py)?;
        out.push(Item::Text(format!(
            "F_ROH: {} ({long} MB in runs of 5 MB or more); nearest: {band}",
            pycompat::fmt::fixed(num(Some(f)).unwrap_or_default(), 4)
        )));
    }
    if let Some(u) = d.get("upd_flag")?.filter(|v| truthy(v)) {
        let p = format!("{PATH}.upd_flag");
        let chrom = py(need(u, "chrom", &p)?)?;
        out.push(Item::Text(format!(
            "** {chrom} carries {} MB, more than 50% of the length in runs of 5 MB or more: both copies may \
             come from one parent (uniparental isodisomy) **",
            py(need(u, "mb", &p)?)?
        )));
        let imprinted = get(u, "imprinted", &p)?;
        for r in list(imprinted, &format!("{p}.imprinted"))? {
            let q = format!("{p}.imprinted[]");
            out.push(Item::Text(format!(
                "Imprinted {}: from the mother: {}; from the father: {} (PMID {})",
                py(need(r, "region", &q)?)?,
                py(need(r, "maternal_upd", &q)?)?,
                py(need(r, "paternal_upd", &q)?)?,
                py(need(r, "pmid", &q)?)?
            )));
        }
        if !imprinted.is_some_and(truthy) {
            out.push(Item::Text(format!(
                "No listed imprinting-disorder region overlaps the flagged runs on {chrom} \
                 (assets/imprinted_regions.tsv); a recessive variant in the run, if there is one, would be \
                 homozygous"
            )));
        }
        out.push(Item::Caveats(vec![ROH_UPD_NOTE]));
    }
    if f_roh.is_some() {
        let mut c = vec![ROH_F_NOTE];
        if d.get("relationship")?.is_some_and(truthy) {
            c.push(ROH_BAND_NOTE);
        }
        out.push(Item::Caveats(c));
    }
    if let Some(m) = d.get("masked_mb")?.filter(|v| truthy(v)) {
        out.push(Item::Text(format!(
            "Left out in centromere and heterochromatin bands: {} MB",
            py(m)?
        )));
    }
    if let Some(n) = d.get("sex_chromosome_segments")?.filter(|v| truthy(v)) {
        out.push(Item::Text(format!(
            "Segments on chrX/chrY, not counted: {}",
            py(n)?
        )));
    }
    let key = notes::roh_note_key(d.get("af_used")?);
    out.push(Item::Caveats(vec![notes::text(key).unwrap_or(ROH_NOTE)]));
    Ok(out)
}

/// The lines under the headline, the same in both reports.
pub fn lines(d: &Section<'_>) -> Result<Vec<String>, Malformed> {
    Ok(items(d)?
        .into_iter()
        .map(|i| match i {
            Item::Text(t) => t,
            Item::Caveats(c) => c.join(" "),
        })
        .collect())
}

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "roh", Some("Runs of Homozygosity"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if head.ran {
        let d = s.section("roh");
        let long = list(
            Some(d.need("autosomal_over_5mb")?),
            &format!("{PATH}.autosomal_over_5mb"),
        )?;
        w.push(format!(
            "  Segments: {}  Total: {} MB  Largest: {} MB",
            py(d.need("segments")?)?,
            py(d.need("total_mb")?)?,
            py(d.need("largest_mb")?)?
        ));
        w.push(format!("  Autosomal ROH > 5MB: {}", long.len()));
        for r in long {
            w.push(format!(
                "    {}",
                region(r, &format!("{PATH}.autosomal_over_5mb[]"))?
            ));
        }
        w.extend(lines(&d)?.into_iter().map(|l| format!("  {l}")));
        w.push(String::new());
    }
    Ok(w)
}

/// `f"{d.get(key, 0):.1f} MB"`.
fn mb(d: &Section<'_>, key: &str) -> Result<String, Malformed> {
    let x = match d.get(key)? {
        None => 0.0,
        v => num(v).ok_or_else(|| Malformed(format!("{PATH}.{key} is not a number")))?,
    };
    Ok(format!("{} MB", pycompat::fmt::fixed(x, 1)))
}

/// The HTML report's card, shown whether or not the step ran.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let d = s.section("roh");
    let mut card = Card {
        title: "Runs of Homozygosity".to_owned(),
        ..Card::default()
    };
    if d.state.ran() || d.state == State::Other {
        let done = if d.state == State::Stale {
            ("Stale", Colour::Yellow)
        } else {
            ("Complete", Colour::Green)
        };
        let long = list(
            d.get("autosomal_over_5mb")?,
            &format!("{PATH}.autosomal_over_5mb"),
        )?;
        let stat = |l: &str, v: String| Line::Stat(l.to_owned(), Inline::Text(v));
        card.body = vec![
            Line::Stat(
                "Status".to_owned(),
                Inline::Badge(done.0.to_owned(), done.1),
            ),
            stat("ROH total", mb(&d, "total_mb")?),
            stat("ROH largest segment", mb(&d, "largest_mb")?),
            stat(
                "Segments",
                d.get("segments")?.map_or(Ok("0".to_owned()), py)?,
            ),
            stat("Autosomal ROH > 5 MB", long.len().to_string()),
        ];
        for r in long {
            card.body.push(Line::Note(region(
                r,
                &format!("{PATH}.autosomal_over_5mb[]"),
            )?));
        }
        for i in items(&d)? {
            match i {
                Item::Text(t) => card.body.push(Line::Note(t)),
                Item::Caveats(c) => card.caveats.extend(c),
            }
        }
    }
    Ok(HtmlPart {
        cards: vec![html::card(&d, card)?],
        files: Vec::new(),
    })
}
