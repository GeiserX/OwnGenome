//! The quality-control section: mean depth (step 16b), reads kept by the
//! alignment (step 02), the sex measured from the reads by indexcov (step
//! 16) and somalier (step 33) against the declared sex, contamination
//! (VerifyBamID2, step 33), a duplicate or swapped sample, and the
//! sex-chromosome karyotype when a step names one.

use pycompat::json::Value;

use crate::model::{
    Card, Colour, HtmlPart, Inline, Line, Malformed, Section, State, Summary, is_str, num, py,
    py_eq, truthy,
};
use crate::notes::{COVERAGE_NOTE, KARYOTYPE_NOTE, stale_note};
use crate::sections::coverage::mean_depth;
use crate::text;

/// The four summary keys of the section, in the order their notes show.
pub const KEYS: [&str; 4] = ["coverage", "sex_check", "sample_qc", "read_count"];

/// What somalier says when the declared and the measured sex agree.
pub const SEX_CHECK_AGREE: &str = "declared and measured sex agree";

fn or_none(v: Option<&Value>) -> Result<String, Malformed> {
    v.map_or(Ok("None".to_owned()), py)
}

/// `v or default`, as text.
fn or_text(v: Option<&Value>, default: &str) -> Result<String, Malformed> {
    match v.filter(|v| truthy(v)) {
        Some(v) => py(v),
        None => Ok(default.to_owned()),
    }
}

/// Step 02's two counts as one line, in both reports.
pub fn read_count_text(d: &Section<'_>) -> Result<String, Malformed> {
    let mut txt = format!(
        "{} reads in, {} primary records in the BAM",
        py(d.need("input_reads")?)?,
        py(d.need("primary_records")?)?
    );
    if is_str(d.get("read_count_check")?, "mismatch") {
        let m = d.need("missing_reads")?;
        let lost = num(Some(m)).ok_or_else(|| {
            Malformed("sections.read_count.data.missing_reads is not a number".into())
        })?;
        txt += &if lost > 0.0 {
            format!(": {} reads missing", py(m)?)
        } else {
            ": more records than reads in".to_owned()
        };
        txt += " (the BAM was kept by the warn setting, see docs/troubleshooting.md)";
    }
    Ok(txt)
}

/// FREEMIX with its verdict, as both reports print it.
pub fn freemix_text(q: &Section<'_>) -> Result<String, Malformed> {
    let freemix = q.get("freemix")?;
    if is_str(q.get("contamination")?, "not_run") || !freemix.is_some_and(truthy) {
        return Ok("not run".to_owned());
    }
    let freemix = py(freemix.expect("checked above"))?;
    let warn = q.get("freemix_warn_above")?.map_or(Ok(String::new()), py)?;
    let few = if is_str(q.get("verifybamid2_marker_check")?, "skipped") {
        "; fewer than 1,000 panel markers had reads, so VerifyBamID2 ran without its marker check"
    } else {
        ""
    };
    Ok(if is_str(q.get("contamination")?, "warn") {
        format!("{freemix} (above {warn}: possible contamination, see docs/33-sample-qc.md{few})")
    } else {
        format!("{freemix} (warning above {warn}{few})")
    })
}

/// One line per step that names a sex-chromosome karyotype (indexcov, step
/// 16; somalier, step 33). The caller adds [`KARYOTYPE_NOTE`] when there is
/// one.
pub fn karyotype_lines(s: &Summary) -> Result<Vec<String>, Malformed> {
    let mut out = Vec::new();
    let sc = s.section("sex_check");
    if sc.state.ran()
        && let Some(k) = sc.get("karyotype")?.filter(|v| truthy(v))
    {
        out.push(format!(
            "Sex chromosomes: {} (indexcov copy number chrX {}, chrY {})",
            py(k)?,
            or_none(sc.get("cn_chrX")?)?,
            or_none(sc.get("cn_chrY")?)?
        ));
    }
    let sq = s.section("sample_qc");
    if sq.state.ran()
        && let Some(k) = sq.get("karyotype")?.filter(|v| truthy(v))
    {
        out.push(format!(
            "Sex chromosomes: {} (somalier depth over the autosomes chrX {}, chrY {}; one copy reads 0.5)",
            py(k)?,
            or_none(sq.get("x_depth_ratio")?)?,
            or_none(sq.get("y_depth_ratio")?)?
        ));
    }
    Ok(out)
}

/// The text report's block, when any of its four keys is in the summary.
pub fn text(s: &Summary) -> Result<Vec<String>, Malformed> {
    let secs = KEYS.map(|k| s.section(k));
    let [_, sexc, sqc, rdc] = secs;
    if secs.iter().all(|x| x.state == State::Missing) {
        return Ok(Vec::new());
    }
    let mut w = vec!["## Quality control".to_owned(), "---".to_owned()];
    // The fault-new-since-text-qc control drops the keys the badge looks at.
    let since_keys: &[&str] = if cfg!(feature = "fault-new-since-text-qc") {
        &[]
    } else {
        &KEYS
    };
    if let Some(line) = text::new_since(s, since_keys)? {
        w.push(line);
    }
    for sec in &secs {
        let note = stale_note(sec)?;
        if !note.is_empty() {
            w.push(format!("  [{}: {note}]", py(sec.title)?));
        }
    }
    w.push(match mean_depth(s)? {
        Some(md) => format!("  Mean depth: {md}"),
        None => "  Mean depth: not available".to_owned(),
    });
    if rdc.state.ran() {
        let flag = if is_str(rdc.get("read_count_check")?, "mismatch") {
            "  ** READS LOST IN ALIGNMENT **"
        } else {
            ""
        };
        w.push(format!("  Alignment: {}{flag}", read_count_text(&rdc)?));
    }
    let inferred = sexc.get("inferred_sex")?.filter(|v| truthy(v));
    let declared = Some(s.need2("run", "declared_sex")?).filter(|v| truthy(v));
    w.push(format!(
        "  Inferred sex: {}",
        or_text(inferred, "not available")?
    ));
    let clash = match (inferred, declared) {
        (Some(i), Some(d)) if !py_eq(i, d) => "  ** DOES NOT MATCH the inferred sex **",
        _ => "",
    };
    w.push(format!(
        "  Declared sex: {}{clash}",
        or_text(declared, "not declared")?
    ));
    if sqc.state.ran() {
        let check = sqc.get("sex_check")?;
        let tail = if is_str(check, "mismatch") {
            "  ** DOES NOT MATCH the declared sex **".to_owned()
        } else if is_str(check, "ok") {
            format!(" ({SEX_CHECK_AGREE})")
        } else {
            String::new()
        };
        w.push(format!(
            "  Sex from the reads (somalier): {}{tail}",
            or_none(sqc.get("inferred_sex")?)?
        ));
        if is_str(check, "not_checked") {
            w.push(format!(
                "    not checked: {}",
                or_none(sqc.get("sex_check_reason")?)?
            ));
        }
        w.push(format!(
            "  Contamination (VerifyBamID2 FREEMIX): {}",
            freemix_text(&sqc)?
        ));
        if let Some(p) = sqc.get("same_person_as")?.filter(|v| truthy(v)) {
            w.push(format!(
                "  ** The same person as: {} (a duplicate or a sample swap) **",
                py(p)?
            ));
        }
    }
    let kl = karyotype_lines(s)?;
    let karyotype = !kl.is_empty();
    w.extend(kl.into_iter().map(|l| format!("  {l}")));
    if karyotype {
        w.push(format!("  {KARYOTYPE_NOTE}"));
    }
    w.push(format!("  {COVERAGE_NOTE}"));
    w.push(String::new());
    Ok(w)
}

fn stat(label: &str, value: Inline) -> Line {
    Line::Stat(label.to_owned(), value)
}

fn badge(text: String, colour: Colour) -> Inline {
    Inline::Badge(text, colour)
}

/// The HTML report's card. It is always there, so a run with none of the
/// four steps still says what it could not measure.
pub fn html(s: &Summary) -> Result<HtmlPart, Malformed> {
    let [cov, sexc, sqc, rdc] = KEYS.map(|k| s.section(k));
    let inferred = sexc.get("inferred_sex")?.filter(|v| truthy(v));
    let declared = Some(s.need2("run", "declared_sex")?).filter(|v| truthy(v));
    let sex = match (inferred, declared) {
        (Some(i), Some(d)) => badge(
            py(i)?,
            if py_eq(i, d) {
                Colour::Green
            } else {
                Colour::Red
            },
        ),
        _ => Inline::Text(or_text(inferred, "not available")?),
    };
    let mut body = vec![
        stat(
            "Mean depth",
            Inline::Text(mean_depth(s)?.unwrap_or_else(|| "not available".to_owned())),
        ),
        stat("Inferred sex (indexcov)", sex),
        stat(
            "Declared sex",
            Inline::Text(or_text(declared, "not declared")?),
        ),
    ];
    if rdc.state.ran() {
        let colour = if is_str(rdc.get("read_count_check")?, "ok") {
            Colour::Green
        } else {
            Colour::Red
        };
        body.insert(
            1,
            stat(
                "Reads kept by alignment",
                badge(read_count_text(&rdc)?, colour),
            ),
        );
    }
    if sqc.state.ran() {
        let check = sqc.get("sex_check")?;
        let colour = if is_str(check, "ok") {
            Colour::Green
        } else if is_str(check, "mismatch") {
            Colour::Red
        } else {
            Colour::Gray
        };
        let said = sqc.get("inferred_sex")?.map_or(Ok("unknown".into()), py)?;
        body.push(stat("Sex from the reads (somalier)", badge(said, colour)));
        if is_str(check, "not_checked") {
            let why = sqc.get("sex_check_reason")?.map_or(Ok(String::new()), py)?;
            body.push(stat("Not checked", Inline::Quiet(why)));
        }
        let contamination = sqc.get("contamination")?;
        let colour = if is_str(contamination, "ok") {
            Colour::Green
        } else if is_str(contamination, "warn") {
            Colour::Yellow
        } else {
            Colour::Gray
        };
        body.push(stat(
            "Contamination (FREEMIX)",
            badge(freemix_text(&sqc)?, colour),
        ));
        if let Some(p) = sqc.get("same_person_as")?.filter(|v| truthy(v)) {
            body.push(stat("The same person as", badge(py(p)?, Colour::Red)));
        }
    } else {
        body.push(stat(
            "Identity and contamination (step 33)",
            badge("Not run".to_owned(), Colour::Gray),
        ));
    }
    let mut caveats = Vec::new();
    for line in karyotype_lines(s)? {
        let named = line.split_once(": ").map_or("", |(_, b)| b).to_owned();
        body.push(stat("Sex chromosomes", badge(named, Colour::Red)));
        if caveats.is_empty() {
            caveats.push(KARYOTYPE_NOTE);
        }
    }
    caveats.push(COVERAGE_NOTE);
    let mut stale = Vec::new();
    for sec in [&cov, &sexc, &sqc, &rdc] {
        let note = stale_note(sec)?;
        if !note.is_empty() {
            stale.push(format!("{}: {note}", py(sec.title)?));
        }
    }
    Ok(HtmlPart {
        cards: vec![Card {
            title: "Quality Control".to_owned(),
            caveats,
            stale,
            body,
            ..Card::default()
        }],
        files: Vec::new(),
    })
}
