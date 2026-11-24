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

/// The label for a similarity signal that did not settle identity.
pub const POSSIBLE_MATCH: &str = "Possible match to";

/// Why a same-person match is only possible when identity was checked on
/// another pair: the two share fewer than og-summary's
/// `MIN_IDENTITY_SITES` genotyped sites (a test keeps the two in step).
pub const THIN_MATCH: &str = "on fewer than 1,000 shared sites";

/// The other samples at the same-person relatedness (`same_person_as`),
/// split into the firm matches and the possible ones, each a list joined
/// with ", " or `None`. With identity uncertain, every match is possible
/// and the identity line says why. Otherwise the matches step 33 names in
/// `possible_same_person_as` (their pair shares too few sites) are
/// possible, shown with [`THIN_MATCH`] after them, and the rest are firm.
pub fn same_person(q: &Section<'_>) -> Result<(Option<String>, Option<String>), Malformed> {
    let listed = |key: &str| -> Result<Vec<String>, Malformed> {
        Ok(match q.get(key)?.filter(|v| truthy(v)) {
            Some(v) => py(v)?
                .split(", ")
                .filter(|p| !p.is_empty())
                .map(str::to_owned)
                .collect(),
            None => Vec::new(),
        })
    };
    let all = listed("same_person_as")?;
    let joined = |v: Vec<String>| (!v.is_empty()).then(|| v.join(", "));
    if identity_uncertainty(q)?.is_some() {
        return Ok((None, joined(all)));
    }
    let thin = listed("possible_same_person_as")?;
    let (possible, firm): (Vec<String>, Vec<String>) =
        all.into_iter().partition(|p| thin.contains(p));
    Ok((
        joined(firm),
        joined(possible).map(|p| format!("{p} ({THIN_MATCH})")),
    ))
}

/// The identity verdict and the producer's reason when identity is uncertain.
pub fn identity_uncertainty(q: &Section<'_>) -> Result<Option<String>, Malformed> {
    if !is_str(q.get("identity")?, "uncertain") {
        return Ok(None);
    }
    Ok(Some(match q.get("identity_reason")? {
        Some(Value::Str(reason)) if !reason.is_empty() => format!("uncertain ({reason})"),
        _ => "uncertain".to_owned(),
    }))
}

/// Whether VerifyBamID2's estimate relies on uncertain evidence.
pub fn contamination_is_uncertain(q: &Section<'_>) -> Result<bool, Malformed> {
    Ok(is_str(q.get("contamination_confidence")?, "uncertain"))
}

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
    let contamination = q.get("contamination")?;
    let warned = is_str(contamination, "warn");
    let uncertainty = if contamination_is_uncertain(q)? {
        Some(match q.get("contamination_confidence_reason")? {
            Some(Value::Str(reason)) if !reason.is_empty() => format!(": {reason}"),
            _ => String::new(),
        })
    } else {
        None
    };
    if is_str(contamination, "ok")
        && let Some(reason) = &uncertainty
    {
        let mut result = "uncertain".to_owned();
        if let Some(value) = freemix.filter(|v| truthy(v)) {
            result += &format!(" (FREEMIX {})", py(value)?);
        }
        result += reason;
        return Ok(result);
    }
    if is_str(contamination, "not_run") || !freemix.is_some_and(truthy) {
        return Ok("not run".to_owned());
    }
    let freemix = py(freemix.expect("checked above"))?;
    let warn = q.get("freemix_warn_above")?.map_or(Ok(String::new()), py)?;
    let few = if is_str(q.get("verifybamid2_marker_check")?, "skipped") {
        "; fewer than 1,000 panel markers had reads, so VerifyBamID2 ran without its marker check"
    } else {
        ""
    };
    let mut result = if warned {
        format!("{freemix} (above {warn}: possible contamination, see docs/33-sample-qc.md{few})")
    } else {
        format!("{freemix} (warning above {warn}{few})")
    };
    if warned && let Some(reason) = uncertainty {
        result += ". The estimate is uncertain";
        result += &reason;
    }
    Ok(result)
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
        if let Some(result) = identity_uncertainty(&sqc)? {
            w.push(format!("  Sample identity: {result}"));
        }
        let (firm, possible) = same_person(&sqc)?;
        if let Some(p) = firm {
            w.push(format!(
                "  ** The same person as: {p} (a duplicate or a sample swap) **"
            ));
        }
        if let Some(p) = possible {
            w.push(format!("  {POSSIBLE_MATCH}: {p}"));
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
        let colour = if contamination_is_uncertain(&sqc)? || is_str(contamination, "warn") {
            Colour::Yellow
        } else if is_str(contamination, "ok") {
            Colour::Green
        } else {
            Colour::Gray
        };
        body.push(stat(
            "Contamination (FREEMIX)",
            badge(freemix_text(&sqc)?, colour),
        ));
        if let Some(result) = identity_uncertainty(&sqc)? {
            body.push(stat("Sample identity", badge(result, Colour::Yellow)));
        }
        let (firm, possible) = same_person(&sqc)?;
        if let Some(p) = firm {
            body.push(stat("The same person as", badge(p, Colour::Red)));
        }
        if let Some(p) = possible {
            body.push(stat(POSSIBLE_MATCH, badge(p, Colour::Yellow)));
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

#[cfg(test)]
mod tests {
    //! From step 33's own table to both reports: a same-person relatedness
    //! on too few shared sites reads as a possible match, in yellow, even
    //! when another pair of the run makes identity `checked`; a pair
    //! compared on enough sites stays the same person, in red.
    use crate::model::Summary;
    use crate::{html, text};
    use og_summary::qc::{FREEMIX_WARN, Inputs, table};

    const SAMPLES: &str = "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\t\
original_pedigree_sex\tgt_depth_mean\tn_hom_ref\tn_het\tn_hom_alt\tp_middling_ab\t\
X_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\tY_depth_mean\tY_n\n\
S\tS\t-9\t-9\t1\t-9\t-9\t30.0\t5000\t6000\t4000\t0.010\t15.0\t39\t0\t0\t39\t14.0\t5\n";

    /// The keys of the table the QC card reads for identity.
    const KEYS: [&str; 4] = [
        "identity",
        "identity_reason",
        "same_person_as",
        "possible_same_person_as",
    ];

    /// Step 33's table for S with these pairs, as summary.json's
    /// `sample_qc` data, and both reports of that summary.
    fn reports(tag: &str, pairs: &str) -> [String; 2] {
        let d = std::env::temp_dir().join(format!("og-render-thin-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        let (sp, pp) = (d.join("s.samples.tsv"), d.join("s.pairs.tsv"));
        std::fs::write(&sp, SAMPLES).unwrap();
        std::fs::write(&pp, pairs).unwrap();
        let rows = table(&Inputs {
            sample: "S",
            samples_tsv: sp.to_str().unwrap(),
            pairs_tsv: pp.to_str(),
            freemix_warn: FREEMIX_WARN,
            ..Inputs::default()
        })
        .unwrap();
        let _ = std::fs::remove_dir_all(&d);
        let data = rows
            .iter()
            .filter(|(k, _)| KEYS.contains(k))
            .map(|(k, v)| format!("\"{k}\": \"{v}\""))
            .collect::<Vec<_>>()
            .join(", ");
        let summary = Summary::parse(&format!(
            r#"{{"sample": "S", "generated_utc": "2026-01-02 03:04 UTC",
                "run": {{"started_utc": null, "declared_sex": "male"}},
                "databases": {{"clinvar_release": null, "cpsr_clinvar_release": null, "hla_database": null}},
                "cyp2d6": null, "manifest": [], "not_run": [], "not_assessed": [],
                "sections": {{"sample_qc": {{"state": "ok", "data": {{{data}}}}}}}}}"#
        ))
        .unwrap();
        [
            text::report(&summary).unwrap(),
            html::report(&summary).unwrap(),
        ]
    }

    fn stat(label: &str, colour: &str, value: &str) -> String {
        format!(
            r#"<span class="label">{label}</span><span class="value"><span class="badge badge-{colour}">{value}<"#
        )
    }

    #[test]
    fn a_thin_pair_is_a_possible_match_beside_a_well_compared_one() {
        // T: the same person on 12,000 shared sites; U: 0.97 on 300; V:
        // unrelated on 15,000. T and V make identity `checked`.
        let [txt, page] = reports(
            "mixed",
            "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\n\
             S\tT\t0.99\t0\t300\t12000\nU\tS\t0.97\t0\t300\t300\nS\tV\t0.1\t9\t300\t15000\n",
        );
        assert!(
            txt.contains("  ** The same person as: T (0.99) (a duplicate or a sample swap) **\n"),
            "{txt}"
        );
        assert!(
            txt.contains("  Possible match to: U (0.97) (on fewer than 1,000 shared sites)\n"),
            "{txt}"
        );
        assert!(!txt.contains("Sample identity"), "{txt}");
        assert!(
            page.contains(&stat("The same person as", "red", "T (0.99)")),
            "{page}"
        );
        assert!(
            page.contains(&stat(
                "Possible match to",
                "yellow",
                "U (0.97) (on fewer than 1,000 shared sites)"
            )),
            "{page}"
        );
        for report in [&txt, &page] {
            assert!(!report.contains("U (0.97), "), "U is never firm: {report}");
            assert!(!report.contains("T (0.99), U"), "U is never firm: {report}");
        }
    }

    #[test]
    fn well_compared_pairs_stay_the_same_person_in_red() {
        let [txt, page] = reports(
            "firm",
            "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\n\
             S\tT\t0.99\t0\t300\t12000\nU\tS\t0.95\t0\t300\t1000\n",
        );
        assert!(
            txt.contains(
                "  ** The same person as: T (0.99), U (0.95) (a duplicate or a sample swap) **\n"
            ),
            "{txt}"
        );
        assert!(!txt.contains("Possible match to"), "{txt}");
        assert!(
            page.contains(&stat("The same person as", "red", "T (0.99), U (0.95)")),
            "{page}"
        );
        assert!(!page.contains("Possible match to"), "{page}");
    }

    #[test]
    fn thin_pairs_alone_are_possible_matches_with_identity_uncertain() {
        // No pair shares enough: identity is uncertain and its line says
        // why, so the match carries no second reason.
        let [txt, page] = reports(
            "uncertain",
            "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\nU\tS\t0.97\t0\t300\t300\n",
        );
        assert!(txt.contains("  Sample identity: uncertain ("), "{txt}");
        assert!(txt.contains("  Possible match to: U (0.97)\n"), "{txt}");
        let identity = stat("Sample identity", "yellow", "uncertain (");
        assert!(page.contains(identity.trim_end_matches('<')), "{page}");
        assert!(
            page.contains(&stat("Possible match to", "yellow", "U (0.97)")),
            "{page}"
        );
        for report in [&txt, &page] {
            assert!(!report.contains("The same person as"), "{report}");
        }
    }

    #[test]
    fn the_thin_match_note_names_the_identity_floor() {
        assert_eq!(og_summary::qc::MIN_IDENTITY_SITES, 1000);
        assert_eq!(super::THIN_MATCH, "on fewer than 1,000 shared sites");
    }
}
