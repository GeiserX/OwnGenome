//! The structural-variant callers of the structural-variants section: Manta
//! (step 04), Delly (step 19), CNVpytor (step 18) and the consensus merge
//! (step 22), with the caveat on what short reads find.
//!
//! The section also shows AnnotSV's class 4 and 5 calls, which
//! `sections::annotsv` writes. Its writers in `chapters::writer` put the
//! two together: [`text_blocks`] then AnnotSV's block in the text report,
//! and [`card`] then AnnotSV's card in the HTML report.

use crate::model::{Card, Inline, Line, Malformed, Section, Summary, get, py, truthy};
use crate::notes::{CNVPYTOR_NOTE, SV_CALLERS_NOTE, stale_note};
use crate::text;

/// The four summary keys of the callers, in the order both reports show
/// them.
pub const KEYS: [&str; 4] = ["manta", "delly", "cnvpytor", "sv_consensus"];

/// CNVpytor's significant calls left after its filters, with what each
/// filter left out; `None` for a summary from before the filters.
pub fn cnvpytor_filtered_text(d: &Section<'_>) -> Result<Option<String>, Malformed> {
    let Some(kept) = d.get("significant_filtered")? else {
        return Ok(None);
    };
    let path = "sections.cnvpytor.data.dropped";
    let dropped = d.get("dropped")?.filter(|v| truthy(v));
    let mut parts = Vec::new();
    for (k, label) in [
        ("q0", "q0 0.5 or more"),
        ("gap", "over N gaps (pN 0.5 or more)"),
        ("band", "over a centromere or heterochromatin band"),
        ("chry_female", "on chrY of a declared female"),
    ] {
        if let Some(n) = dropped
            .map(|dr| get(dr, k, path))
            .transpose()?
            .flatten()
            .filter(|v| truthy(v))
        {
            parts.push(format!("{} {label}", py(n)?));
        }
    }
    let mut out = py(kept)?;
    if !parts.is_empty() {
        out += &format!(" (left out: {})", parts.join(", "));
    }
    if !d.get("declared_sex")?.is_some_and(truthy) {
        out += "; no declared sex, so chrY calls are kept";
    }
    Ok(Some(out))
}

/// The text report's blocks of the four callers. The caveat is printed
/// once, under the first of Manta, Delly and CNVpytor with results.
pub fn text_blocks(s: &Summary) -> Result<Vec<String>, Malformed> {
    let mut w = Vec::new();
    let mut note = Some(format!("  {SV_CALLERS_NOTE}"));
    let n = |sec: &Section<'_>, k: &str| sec.need(k).and_then(py);
    for (key, title) in [
        ("manta", "Structural Variants (Manta)"),
        ("delly", "Structural Variants (Delly)"),
    ] {
        if let Some(head) = text::head(s, key, Some(title))? {
            w.extend(head.lines);
            if head.ran {
                let d = s.section(key);
                w.push(format!(
                    "  Total SVs: {}  PASS: {}",
                    n(&d, "total")?,
                    n(&d, "pass")?
                ));
                w.extend(note.take());
                w.push(String::new());
            }
        }
    }
    if let Some(head) = text::head(s, "cnvpytor", Some("Copy Number Variants (CNVpytor)"))? {
        w.extend(head.lines);
        if head.ran {
            let d = s.section("cnvpytor");
            w.push(format!(
                "  Total CNVs: {} ({} deletions, {} duplications)",
                n(&d, "total")?,
                n(&d, "deletions")?,
                n(&d, "duplications")?
            ));
            w.push(format!(
                "  Significant (e-val < 0.01): {}",
                n(&d, "significant")?
            ));
            w.extend(note.take());
            if let Some(filtered) = cnvpytor_filtered_text(&d)?.filter(|t| !t.is_empty()) {
                w.push(format!("  Significant after filters: {filtered}"));
                w.push(format!("  {CNVPYTOR_NOTE}"));
            }
            w.push(String::new());
        }
    }
    if let Some(head) = text::head(s, "sv_consensus", Some("SV Consensus Merge"))? {
        w.extend(head.lines);
        if head.ran {
            let d = s.section("sv_consensus");
            w.push(format!(
                "  Consensus SVs (2+ callers): {}",
                n(&d, "consensus")?
            ));
            w.push(String::new());
        }
    }
    Ok(w)
}

/// `d.get(key, "N/A")` as text.
fn or_na(d: &Section<'_>, key: &str) -> Result<String, Malformed> {
    d.get(key)?.map_or(Ok("N/A".to_owned()), py)
}

/// The HTML report's card of the four callers, shown whether or not they
/// ran: each one's counts, or N/A, and the note of each result that is not
/// this run's own.
pub fn card(s: &Summary) -> Result<Card, Malformed> {
    let [m, dl, cn, sv] = KEYS.map(|k| s.section(k));
    let stat = |l: &str, v: String| Line::Stat(l.to_owned(), Inline::Text(v));
    let mut c = Card {
        title: "Structural Variants".to_owned(),
        ..Card::default()
    };
    for x in [&m, &dl, &cn, &sv] {
        let note = stale_note(x)?;
        if !note.is_empty() {
            c.stale.push(format!("{}: {note}", py(x.title)?));
        }
    }
    c.body = vec![
        stat(
            "Manta SVs (total / PASS)",
            format!("{} / {}", or_na(&m, "total")?, or_na(&m, "pass")?),
        ),
        stat(
            "Delly SVs (total / PASS)",
            format!("{} / {}", or_na(&dl, "total")?, or_na(&dl, "pass")?),
        ),
        stat("CNVpytor CNVs", or_na(&cn, "total")?),
        stat(
            "CNVpytor significant (e-val < 0.01)",
            or_na(&cn, "significant")?,
        ),
    ];
    if let Some(filtered) = cnvpytor_filtered_text(&cn)?.filter(|t| !t.is_empty()) {
        c.body
            .push(stat("CNVpytor significant after filters", filtered));
        c.caveats.push(CNVPYTOR_NOTE);
    }
    c.body
        .push(stat("Consensus SVs (2+ callers)", or_na(&sv, "consensus")?));
    c.caveats.push(SV_CALLERS_NOTE);
    Ok(c)
}
