//! Step 06's ClinVar screen: the count of P/LP hits with the call floor and
//! the gnomAD frequency, the hits best-reviewed first (up to 50, then any
//! hit at a known lab-construct site), the cDNA read check and the cDNA
//! scan ([`super::cdna`]) and the known deep intronic sites
//! ([`super::deep_intronic`]). The conflicting records some labs call P/LP
//! (`conflicting`) get a review list of their own, apart from every count.
//!
//! It also holds what the clinical sections share: [`D`], a dict of
//! summary.json read with Python's rules, and the Viewer cell of a findings
//! table.

use pycompat::fmt::{general, ljust, percent};
use pycompat::json::Value;

use crate::model::{
    Card, Cell, Colour, HtmlPart, Inline, Layout, Line, Malformed, Section, Summary, Table,
    contains, num, py, py_eq, truthy,
};
use crate::notes::{CLINVAR_FREQ_NOTE, DEEP_INTRONIC_NOTE};
use crate::sections::{cdna, cdna_scan, deep_intronic};
use crate::{html, text};

/// The hits the report lists, and the anchor of their table (the charts
/// link a listed hit there).
pub const HITS_SHOWN: usize = 50;
pub const HITS_ANCHOR: &str = "clinvar-hits";

// --- a dict of summary.json, read with Python's rules ---------------------------------

/// A value of summary.json and the path that names it in an error.
#[derive(Debug, Clone)]
pub struct D<'a> {
    pub v: &'a Value,
    pub path: String,
}

impl<'a> D<'a> {
    pub fn new(v: &'a Value, path: impl Into<String>) -> D<'a> {
        D {
            v,
            path: path.into(),
        }
    }

    /// The data of section `sec`.
    pub fn of(sec: &Section<'a>) -> D<'a> {
        D::new(sec.data, format!("sections.{}.data", sec.key))
    }

    fn child(&self, v: &'a Value, key: &str) -> D<'a> {
        D::new(v, format!("{}.{key}", self.path))
    }

    /// `d[key]`.
    pub fn at(&self, key: &str) -> Result<D<'a>, Malformed> {
        crate::model::need(self.v, key, &self.path).map(|v| self.child(v, key))
    }

    /// `d.get(key)`: `None` for a missing key.
    pub fn get(&self, key: &str) -> Result<Option<D<'a>>, Malformed> {
        Ok(crate::model::get(self.v, key, &self.path)?.map(|v| self.child(v, key)))
    }

    /// `d.get(key)` when it is true, else `None`.
    pub fn opt(&self, key: &str) -> Result<Option<D<'a>>, Malformed> {
        Ok(self.get(key)?.filter(|x| truthy(x.v)))
    }

    /// `bool(d.get(key))`.
    pub fn flag(&self, key: &str) -> Result<bool, Malformed> {
        Ok(self.opt(key)?.is_some())
    }

    /// `key in d`.
    pub fn has(&self, key: &str) -> bool {
        contains(self.v, key)
    }

    /// `str(d[key])`.
    pub fn s(&self, key: &str) -> Result<String, Malformed> {
        self.at(key)?.text()
    }

    /// `str(d.get(key))`: `None` for a missing key.
    pub fn s_or_none(&self, key: &str) -> Result<String, Malformed> {
        self.get(key)?.map_or(Ok("None".to_owned()), |x| x.text())
    }

    /// `str(d.get(key, default))`.
    pub fn s_or(&self, key: &str, default: &str) -> Result<String, Malformed> {
        self.get(key)?.map_or(Ok(default.to_owned()), |x| x.text())
    }

    /// `str(v)`.
    pub fn text(&self) -> Result<String, Malformed> {
        py(self.v)
    }

    /// `v` where Python needs a string (`html.escape(v)`, `s + v`).
    pub fn string(&self) -> Result<&'a str, Malformed> {
        match self.v {
            Value::Str(s) => Ok(s),
            _ => Err(Malformed(format!("{} is not a string", self.path))),
        }
    }

    /// `bool(v)`.
    pub fn truthy(&self) -> bool {
        truthy(self.v)
    }

    /// `v is None`.
    pub fn is_none(&self) -> bool {
        matches!(self.v, Value::Null)
    }

    /// `v == text`.
    pub fn is(&self, text: &str) -> bool {
        matches!(self.v, Value::Str(s) if s == text)
    }

    /// The items of a list.
    pub fn items(&self) -> Result<Vec<D<'a>>, Malformed> {
        match self.v {
            Value::Array(a) => Ok(a
                .iter()
                .enumerate()
                .map(|(i, v)| D::new(v, format!("{}[{i}]", self.path)))
                .collect()),
            _ => Err(Malformed(format!("{} is not a list", self.path))),
        }
    }

    /// The items of a dict, in its order.
    pub fn pairs(&self) -> Result<Vec<(&'a str, D<'a>)>, Malformed> {
        match self.v {
            Value::Object(o) => Ok(o
                .iter()
                .map(|(k, v)| (k.as_str(), self.child(v, k)))
                .collect()),
            _ => Err(Malformed(format!("{} is not a dict", self.path))),
        }
    }

    /// `len(v)` of a list, dict or string.
    pub fn len(&self) -> Result<usize, Malformed> {
        match self.v {
            Value::Array(a) => Ok(a.len()),
            Value::Object(o) => Ok(o.len()),
            Value::Str(s) => Ok(s.chars().count()),
            _ => Err(Malformed(format!("{} has no length", self.path))),
        }
    }

    /// `not v` of a list, dict or string.
    pub fn is_empty(&self) -> Result<bool, Malformed> {
        Ok(self.len()? == 0)
    }

    /// The number `v` is.
    pub fn num(&self) -> Result<f64, Malformed> {
        num(Some(self.v)).ok_or_else(|| Malformed(format!("{} is not a number", self.path)))
    }

    /// `v == n` for a number `n`.
    pub fn is_num(&self, n: i128) -> bool {
        py_eq(self.v, &Value::Int(n))
    }
}

/// `d.get(key) or []` as a list.
pub fn list_or_empty<'a>(d: &D<'a>, key: &str) -> Result<Vec<D<'a>>, Malformed> {
    d.opt(key)?.map_or(Ok(Vec::new()), |x| x.items())
}

/// `f"{v:<width}"`: a string, a number or a bool padded to `width`.
pub fn lj(v: &D<'_>, width: usize) -> Result<String, Malformed> {
    let t = match v.v {
        Value::Str(s) => s.clone(),
        Value::Bool(b) => u8::from(*b).to_string(),
        Value::Int(_) | Value::Float(_) => v.text()?,
        _ => return Err(Malformed(format!("{} cannot be padded", v.path))),
    };
    Ok(ljust(&t, width))
}

// --- the Viewer column ----------------------------------------------------------------

fn name_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || c == b'_' || c == b'.'
}

/// Every `name:start` or `name:start-end` in `text`, left to right, as the
/// Python's `([A-Za-z0-9_.]+):([0-9]+)(?:-([0-9]+))?(?![0-9])` finds them.
pub fn loci_in(text: &str) -> Vec<(&str, &str, Option<&str>)> {
    let b = text.as_bytes();
    let digits = |from: usize| {
        let mut j = from;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        j
    };
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if !name_char(b[i]) {
            i += 1;
            continue;
        }
        let mut j = i;
        while j < b.len() && name_char(b[j]) {
            j += 1;
        }
        // Every start inside the run ends the name where the run ends.
        if j < b.len() && b[j] == b':' && j + 1 < b.len() && b[j + 1].is_ascii_digit() {
            let k = digits(j + 1);
            let mut end = None;
            let mut next = k;
            if k + 1 < b.len() && b[k] == b'-' && b[k + 1].is_ascii_digit() {
                let m = digits(k + 1);
                end = Some(&text[k + 1..m]);
                next = m;
            }
            out.push((&text[i..j], &text[j + 1..k], end));
            i = next;
        } else {
            i = j.max(i + 1);
        }
    }
    out
}

/// The Viewer cell of a row whose locus is written in its text: the span of
/// every position on the first locus's chromosome.
pub fn viewer_of(text: &str) -> Result<Cell, Malformed> {
    let found = loci_in(text);
    let Some((chrom, _, _)) = found.first().copied() else {
        return Ok(Cell::Viewer(None));
    };
    let pos: Vec<i128> = found
        .iter()
        .filter(|(c, _, _)| *c == chrom)
        .flat_map(|(_, a, b)| [Some(*a), *b])
        .flatten()
        .map(|x| x.parse::<i128>().unwrap_or(i128::MAX))
        .collect();
    let lo = pos.iter().min().copied().unwrap_or_default();
    let hi = pos.iter().max().copied().unwrap_or_default();
    Cell::viewer(
        &Value::Str(chrom.to_owned()),
        &Value::Int(lo),
        Some(&Value::Int(hi)),
    )
}

/// `str(text or "")` for [`viewer_of`].
pub fn locus_text(v: &D<'_>) -> Result<String, Malformed> {
    if v.truthy() {
        v.text()
    } else {
        Ok(String::new())
    }
}

// --- the hits -------------------------------------------------------------------------

/// A hit's combined-allele note (step 06), `""` when it has none.
pub fn combined_cell(h: &D<'_>) -> Result<String, Malformed> {
    match h.get("combined_allele")? {
        None => Ok(String::new()),
        Some(n) if n.is(".") || n.is("") => Ok(String::new()),
        Some(n) => n.text(),
    }
}

/// A hit's gnomAD allele frequency, `""` for a summary written before the
/// field.
pub fn af_text(h: &D<'_>) -> Result<String, Malformed> {
    let Some(af) = h.get("af")? else {
        return Ok(String::new());
    };
    if af.is_none() {
        return Ok("gnomAD AF unknown".to_owned());
    }
    let mut t = format!("gnomAD AF {}", general(af.num()?, 3));
    if h.flag("common")? {
        t += "; common: unlikely to be a high-penetrance finding";
    }
    Ok(t)
}

/// The low-confidence mark of a ClinVar or CPSR row, `""` for a confident
/// one.
pub fn call_text(h: &D<'_>) -> Result<String, Malformed> {
    if !h.flag("low_confidence")? {
        return Ok(String::new());
    }
    let mut vals = Vec::new();
    for (k, key) in [("GQ", "gq"), ("DP", "dp")] {
        if let Some(v) = h.get(key)?.filter(|v| !v.is_none()) {
            vals.push(format!("{k} {}", v.text()?));
        }
    }
    Ok(format!(
        "low-confidence call ({}): check the reads",
        vals.join(", ")
    ))
}

/// How many hits fall below step 06's call floor, or `None` for a summary
/// written before the floor.
pub fn floor_line(d: &D<'_>) -> Result<Option<String>, Malformed> {
    if !d.has("low_confidence") {
        return Ok(None);
    }
    let (gq, dp) = floor(d)?;
    Ok(Some(format!(
        "Low-confidence calls (GQ below {gq} or DP below {dp}), listed but not counted above: {}",
        d.s("low_confidence")?
    )))
}

/// The call floor's GQ and DP, `None` when the summary has none.
fn floor(d: &D<'_>) -> Result<(String, String), Malformed> {
    match d.opt("call_floor")? {
        Some(f) => Ok((f.s_or_none("GQ")?, f.s_or_none("DP")?)),
        None => Ok(("None".to_owned(), "None".to_owned())),
    }
}

/// Where the hits' gnomAD frequencies came from, or why there are none.
pub fn freq_line(d: &D<'_>) -> Result<Option<String>, Malformed> {
    let Some(fr) = d.opt("frequency")? else {
        return Ok(None);
    };
    if let Some(note) = fr.opt("note")? {
        return Ok(Some(format!(
            "gnomAD AF: unknown for every hit ({})",
            note.text()?
        )));
    }
    let mut common = 0;
    if let Some(hits) = d.get("hits")? {
        for h in hits.items()? {
            if h.flag("common")? {
                common += 1;
            }
        }
    }
    Ok(Some(format!(
        "gnomAD AF from {}: {common} hit(s) above {}",
        fr.s("source")?,
        percent(fr.at("common_af")?.num()?, 0)
    )))
}

// --- the conflicting records: a review list, never counted ----------------------------

/// The title of the conflicting records' list in both reports.
pub const CONFLICTING_TITLE: &str =
    "Conflicting classifications with a P/LP submission: review, not counted";

/// What the list is, under it in both reports.
pub const CONFLICTING_NOTE: &str = "ClinVar calls these records conflicting: at least one lab \
    classed the variant pathogenic or likely pathogenic, and another lab did not. They are not in \
    the P/LP count, the secondary findings, the carrier list or the couple screen. A clinician or \
    genetic counsellor reviews each one (docs/06-clinvar-screen.md).";

/// The list's size, or the note that it was not checked.
pub fn conflicting_count(c: &D<'_>) -> Result<String, Malformed> {
    if c.flag("checked")? {
        return c.s("count");
    }
    match c.opt("note")? {
        Some(n) => n.text(),
        None => Ok("not checked".to_owned()),
    }
}

/// `1 star`, `2 stars`.
fn stars_text(h: &D<'_>) -> Result<String, Malformed> {
    let n = h.s_or("stars", "0")?;
    Ok(format!("{n} star{}", if n == "1" { "" } else { "s" }))
}

/// One record as a line: gene, variant, genotype, the labs' split, stars and
/// condition.
pub fn conflicting_line(h: &D<'_>) -> Result<String, Malformed> {
    let mut out = format!(
        "{} {}:{} {}>{} {}: {}; {}",
        h.s("gene")?,
        h.s("chrom")?,
        h.s("pos")?,
        h.s("ref")?,
        h.s("alt")?,
        h.s("genotype")?,
        h.s("submissions")?,
        stars_text(h)?
    );
    if let Some(c) = h.opt("clinvar_condition")? {
        out += &format!("; {}", c.text()?);
    }
    let call = call_text(h)?;
    if !call.is_empty() {
        out += &format!("; {call}");
    }
    Ok(out)
}

/// The text report's lines, when step 06 wrote the list.
fn conflicting_text(d: &D<'_>) -> Result<Vec<String>, Malformed> {
    let Some(c) = d.get("conflicting")? else {
        return Ok(Vec::new());
    };
    let mut w = vec![format!("  {CONFLICTING_TITLE}: {}", conflicting_count(&c)?)];
    let hits = c.at("hits")?.items()?;
    for h in hits.iter().take(HITS_SHOWN) {
        w.push(format!("    {}", conflicting_line(h)?));
    }
    if hits.len() > HITS_SHOWN {
        w.push(format!(
            "    ... and {} more in clinvar/<sample>_clinvar_conflicting.tsv",
            hits.len() - HITS_SHOWN
        ));
    }
    w.push(format!("  {CONFLICTING_NOTE}"));
    Ok(w)
}

/// The HTML card of the list.
fn conflicting_card(c: &D<'_>) -> Result<Card, Malformed> {
    let hits = c.at("hits")?.items()?;
    let checked = c.flag("checked")?;
    let mut body = vec![Line::Stat(
        "Records to review".to_owned(),
        Inline::Badge(
            conflicting_count(c)?,
            if checked && hits.is_empty() {
                Colour::Gray
            } else {
                Colour::Yellow
            },
        ),
    )];
    let mut rows = Vec::new();
    for h in hits.iter().take(HITS_SHOWN) {
        let st =
            |k: &str| -> Result<Cell, Malformed> { Ok(Cell::Text(h.at(k)?.string()?.to_owned())) };
        rows.push(vec![
            st("chrom")?,
            Cell::Text(h.s("pos")?),
            st("ref")?,
            st("alt")?,
            st("genotype")?,
            st("gene")?,
            st("submissions")?,
            st("review_status")?,
            Cell::Text(h.s_or("stars", "0")?),
            Cell::Text(match h.opt("clinvar_condition")? {
                Some(x) => x.text()?,
                None => String::new(),
            }),
            Cell::Text(call_text(h)?),
            Cell::viewer(h.at("chrom")?.v, h.at("pos")?.v, None)?,
        ]);
    }
    if !rows.is_empty() {
        body.push(Line::Table(Table {
            head: [
                "Chr",
                "Position",
                "Ref",
                "Alt",
                "Genotype",
                "Gene",
                "Lab classifications",
                "Review status",
                "Stars",
                "Condition",
                "Call quality",
            ]
            .map(str::to_owned)
            .to_vec(),
            rows,
            viewer: true,
            layout: Layout::Indented,
        }));
        if hits.len() > HITS_SHOWN {
            body.push(Line::Note(format!(
                "{} more in clinvar/<sample>_clinvar_conflicting.tsv",
                hits.len() - HITS_SHOWN
            )));
        }
        body.push(Line::ViewerNote);
    }
    // An unchecked list has no size: its badge says why, and the title
    // never reads "(0 of 0)".
    let title = if checked {
        format!(
            "{CONFLICTING_TITLE} ({} of {})",
            hits.len().min(HITS_SHOWN),
            hits.len()
        )
    } else {
        CONFLICTING_TITLE.to_owned()
    };
    Ok(Card {
        title,
        full: true,
        caveats: vec![CONFLICTING_NOTE],
        body,
        ..Card::default()
    })
}

// --- the text report ------------------------------------------------------------------

/// The text report's block.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let Some(head) = text::head(s, "clinvar", Some("ClinVar Pathogenic Screen"))? else {
        return Ok(Vec::new());
    };
    let mut w = head.lines;
    if !head.ran {
        return Ok(w);
    }
    let d = D::of(&s.section("clinvar"));
    w.push(format!(
        "  Pathogenic/Likely Pathogenic hits: {}",
        d.s("count")?
    ));
    if let Some(low) = floor_line(&d)? {
        w.push(format!("  {low}"));
    }
    w.push(format!("  ClinVar file date: {}", clinvar_release(s)?));
    let hits = d.at("hits")?;
    if hits.truthy() {
        let hits = hits.items()?;
        let by = d.at("by_stars")?;
        let mut parts = Vec::new();
        for k in ["4", "3", "2", "1", "0"] {
            let plural = if k == "1" { "" } else { "s" };
            parts.push(format!("{k} star{plural} {}", by.s(k)?));
        }
        w.push(format!("  By review status: {}", parts.join(", ")));
        if d.has("cdna_flagged") {
            w.push(format!(
                "  cDNA read check (alt reads clipped at exon edges): {} flagged",
                d.s("cdna_flagged")?
            ));
        }
        if let Some(freq) = freq_line(&d)? {
            w.push(format!("  {freq}"));
        }
        w.push(String::new());
        w.push(
            "  Hits (best-reviewed first, up to 50): gene, position, genotype, significance [review status] stars"
                .to_owned(),
        );
        for h in hits.iter().take(HITS_SHOWN) {
            let note = combined_cell(h)?;
            let mut line = format!(
                "    {} {}:{} {}>{}  {}  {} [{}] {}*",
                lj(&h.at("gene")?, 10)?,
                h.s("chrom")?,
                h.s("pos")?,
                h.s("ref")?,
                h.s("alt")?,
                h.s("genotype")?,
                h.s("significance")?,
                h.s("review_status")?,
                h.s("stars")?
            );
            if !note.is_empty() {
                line += &format!("  ({note})");
            }
            w.push(line);
            let af = af_text(h)?;
            if !af.is_empty() {
                w.push(format!("               frequency: {af}"));
            }
            let call = call_text(h)?;
            if !call.is_empty() {
                w.push(format!("               call: {call}"));
            }
            if h.has("cdna") {
                w.push(format!("               reads: {}", cdna::cdna_text(h)?));
            }
        }
        let extra = cdna::construct_past_cap(&hits, HITS_SHOWN)?;
        if !extra.is_empty() {
            w.push("  Hits at a known lab-construct site past the first 50:".to_owned());
            for h in &extra {
                w.push(format!("    {}", cdna::construct_past_cap_text(h)?));
            }
        }
    }
    w.extend(conflicting_text(&d)?);
    w.extend(cdna_scan::text_lines(&s.section("clinvar"))?);
    if d.has("deep_intronic") {
        w.push("  Known deep intronic sites (always reported): gene, variant, call".to_owned());
        for x in d.at("deep_intronic")?.items()? {
            w.push(format!(
                "    {} {} {}:{} {}>{}  {}",
                lj(&x.at("gene")?, 10)?,
                x.s("hgvs")?,
                x.s("chrom")?,
                x.s("pos")?,
                x.s("ref")?,
                x.s("alt")?,
                deep_intronic::call(&x)?
            ));
        }
        w.push(format!("  {DEEP_INTRONIC_NOTE}"));
    }
    w.push(format!("  {CLINVAR_FREQ_NOTE}"));
    w.push(String::new());
    Ok(w)
}

/// `s['databases'][key] or 'unknown'`.
pub fn database(s: &Summary, key: &str) -> Result<String, Malformed> {
    let v = crate::model::need(s.need("databases")?, key, "summary.databases")?;
    if truthy(v) {
        py(v)
    } else {
        Ok("unknown".to_owned())
    }
}

fn clinvar_release(s: &Summary) -> Result<String, Malformed> {
    database(s, "clinvar_release")
}

// --- the HTML report ------------------------------------------------------------------

fn badge(n: &D<'_>, on: Colour) -> Result<Inline, Malformed> {
    let colour = if n.truthy() { on } else { Colour::Gray };
    Ok(Inline::Badge(n.text()?, colour))
}

/// The screen's card: the counts, never green (no hit is not an
/// all-clear).
fn screen_card(s: &Summary, d: &D<'_>) -> Result<Card, Malformed> {
    let count = d.get("count")?;
    let n = count.clone().map_or(Ok("0".to_owned()), |x| x.text())?;
    let n_on = count.as_ref().is_some_and(D::truthy);
    let mut stars = Vec::new();
    let by = d.get("by_stars")?;
    for k in ["4", "3", "2", "1", "0"] {
        stars.push(match &by {
            Some(by) => by.s_or(k, "0")?,
            None => "0".to_owned(),
        });
    }
    let mut body = vec![
        Line::Stat(
            "ClinVar matches".to_owned(),
            Inline::Badge(n, if n_on { Colour::Yellow } else { Colour::Gray }),
        ),
        Line::Stat(
            "By stars (4/3/2/1/0)".to_owned(),
            Inline::Text(stars.join(" / ")),
        ),
        Line::Stat(
            "ClinVar file date".to_owned(),
            Inline::Text(clinvar_release(s)?),
        ),
    ];
    if d.has("cdna_flagged") {
        body.insert(
            2,
            Line::Stat(
                "Flagged by the cDNA read check".to_owned(),
                badge(&d.at("cdna_flagged")?, Colour::Red)?,
            ),
        );
    }
    if let Some(sites) = d.opt("construct_sites")? {
        body.insert(
            2,
            Line::Stat(
                "Hits at a known lab-construct site (check the reads)".to_owned(),
                Inline::Badge(sites.text()?, Colour::Yellow),
            ),
        );
    }
    if let Some(stat) = cdna_scan::stat(&s.section("clinvar"))? {
        body.insert(2, stat);
    }
    if d.has("low_confidence") {
        let (gq, dp) = floor(d)?;
        body.insert(
            1,
            Line::Stat(
                format!("Low-confidence calls (GQ below {gq} or DP below {dp}), not counted"),
                badge(&d.at("low_confidence")?, Colour::Yellow)?,
            ),
        );
    }
    if let Some(freq) = freq_line(d)? {
        body.push(Line::Note(freq));
    }
    Ok(Card {
        title: "ClinVar Screening".to_owned(),
        caveats: vec![CLINVAR_FREQ_NOTE],
        body,
        ..Card::default()
    })
}

/// The table of the first 50 hits, with the columns the summary carries.
fn hits_card(d: &D<'_>) -> Result<Card, Malformed> {
    let reads = d.has("cdna_flagged");
    let freq = d.has("frequency");
    let call = d.has("low_confidence");
    let hits = d.at("hits")?.items()?;
    let nh = hits.len();
    let mut head: Vec<String> = [
        "Chr",
        "Position",
        "Ref",
        "Alt",
        "Genotype",
        "Gene",
        "Significance",
        "Review status",
        "Stars",
        "Combined allele",
    ]
    .map(str::to_owned)
    .to_vec();
    if freq {
        head.push("gnomAD frequency".to_owned());
    }
    if call {
        head.push("Call quality".to_owned());
    }
    if reads {
        head.push("cDNA read check".to_owned());
    }
    let mut rows = Vec::new();
    for h in hits.iter().take(HITS_SHOWN) {
        let st =
            |k: &str| -> Result<Cell, Malformed> { Ok(Cell::Text(h.at(k)?.string()?.to_owned())) };
        let mut row = vec![
            st("chrom")?,
            Cell::Text(h.s("pos")?),
            st("ref")?,
            st("alt")?,
            st("genotype")?,
            st("gene")?,
            st("significance")?,
            st("review_status")?,
            Cell::Text(h.s("stars")?),
            Cell::Text(combined_cell(h)?),
        ];
        if freq {
            row.push(Cell::Text(af_text(h)?));
        }
        if call {
            row.push(Cell::Text(call_text(h)?));
        }
        if reads {
            row.push(Cell::Text(cdna::cdna_text(h)?));
        }
        row.push(Cell::viewer(h.at("chrom")?.v, h.at("pos")?.v, None)?);
        rows.push(row);
    }
    let mut body = vec![
        Line::Raw(format!("    <a id=\"{HITS_ANCHOR}\"></a>")),
        Line::Table(Table {
            head,
            rows,
            viewer: true,
            layout: Layout::Indented,
        }),
    ];
    for h in cdna::construct_past_cap(&hits, HITS_SHOWN)? {
        body.push(Line::Note(format!(
            "Hit at a known lab-construct site past the first 50: {}",
            cdna::construct_past_cap_text(&h)?
        )));
    }
    if nh > 0 {
        body.push(Line::ViewerNote);
    }
    Ok(Card {
        title: format!(
            "ClinVar Hits (best-reviewed first, {} of {nh})",
            nh.min(HITS_SHOWN)
        ),
        full: true,
        body,
        ..Card::default()
    })
}

/// The known deep intronic sites and the cDNA scan's loci.
fn sites_card(sec: &Section<'_>, d: &D<'_>) -> Result<Card, Malformed> {
    let mut body = Vec::new();
    let mut caveats = Vec::new();
    if d.has("deep_intronic") {
        body.extend(deep_intronic::table(&d.at("deep_intronic")?)?);
        caveats.push(DEEP_INTRONIC_NOTE);
    }
    let (scan, caveat) = cdna_scan::html_lines(sec)?;
    body.extend(scan);
    caveats.extend(caveat);
    Ok(Card {
        title: "Known Deep Intronic Sites and the cDNA Scan".to_owned(),
        full: true,
        caveats,
        body,
        ..Card::default()
    })
}

/// The HTML report's cards: the screen always, the hits table when there
/// are hits, the conflicting records to review when step 06 listed them,
/// and the deep intronic sites and the cDNA scan when step 06 wrote them.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let sec = s.section("clinvar");
    let d = D::of(&sec);
    let mut cards = vec![html::card(&sec, screen_card(s, &d)?)?];
    if sec.state.ran() && d.opt("hits")?.is_some() {
        cards.push(html::card(&sec, hits_card(&d)?)?);
    }
    if sec.state.ran()
        && let Some(c) = d.get("conflicting")?
    {
        cards.push(html::card(&sec, conflicting_card(&c)?)?);
    }
    if sec.state.ran() && (d.has("deep_intronic") || d.opt("cdna_scan")?.is_some()) {
        cards.push(html::card(&sec, sites_card(&sec, &d)?)?);
    }
    Ok(HtmlPart {
        cards,
        files: Vec::new(),
    })
}
