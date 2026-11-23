use super::*;
use og_render::sections::qc;
use pycompat::fmt::fixed;

pub(super) fn card(id: &str, summary: &Summary) -> Read<Adapter> {
    let mut out = Adapter::default();
    match id {
        "qc" => {
            let coverage = summary.section("coverage");
            let sex = summary.section("sex_check");
            let sample = summary.section("sample_qc");
            let reads = summary.section("read_count");
            if let Some(n) = num(coverage.get("mean_depth")?) {
                out.fact("Mean depth", format!("{}x", fixed(n, 1)));
            }
            if reads.state.ran() {
                let tail = if string(reads.data, "read_count_check") == "mismatch" {
                    " (reads lost: see the report)"
                } else {
                    ""
                };
                out.fact(
                    "Reads kept by alignment",
                    format!("{}{tail}", qc::read_count_text(&reads)?),
                );
            }
            let inferred = val(sex.data, "inferred_sex");
            if truthy(inferred) {
                out.fact("Inferred sex (indexcov)", cell(inferred));
            }
            let run = summary.get("run").unwrap_or(&Value::Null);
            out.fact(
                "Declared sex",
                text_or(val(run, "declared_sex"), "not declared"),
            );
            if sample.state.ran() {
                let d = sample.data;
                let sex_check = string(d, "sex_check");
                let result = match sex_check {
                    "ok" => "Declared and measured sex agree".into(),
                    "mismatch" => "Declared and measured sex differ".into(),
                    "not_checked" => format!(
                        "Not checked: {}",
                        text_or(val(d, "sex_check_reason"), "no reason given")
                    ),
                    _ => text_or(val(d, "sex_check"), "."),
                };
                out.fact(
                    "Sex from the reads (somalier)",
                    format!("{}. {result}", text_or(val(d, "inferred_sex"), ".")),
                );
                let contam = string(d, "contamination");
                let contam_uncertain = qc::contamination_is_uncertain(&sample)?;
                let contamination = if contam == "ok" && !contam_uncertain {
                    format!(
                        "Below the contamination flag (FREEMIX {}, flag above {})",
                        text(val(d, "freemix")),
                        text(val(d, "freemix_warn_above"))
                    )
                } else {
                    qc::freemix_text(&sample)?
                };
                out.fact("Contamination (VerifyBamID2)", contamination);
                let identity = qc::identity_uncertainty(&sample)?;
                if let Some(ref result) = identity {
                    out.fact("Sample identity", result);
                }
                let same = val(d, "same_person_as");
                if truthy(same) {
                    let label = if identity.is_some() {
                        qc::POSSIBLE_MATCH
                    } else {
                        "The same person as"
                    };
                    out.fact(label, cell(same));
                }
                out.attention = sex_check == "mismatch"
                    || contam == "warn"
                    || contam_uncertain
                    || identity.is_some()
                    || truthy(same);
            }
            out.lines = qc::karyotype_lines(summary)?;
            out.attention |= !out.lines.is_empty();
            out.links = vec!["multiqc", "fastp", "_fastp"];
        }
        "callability" => {
            let d = summary.section("callability").data;
            let genes = array(d, "not_read")?;
            out.fact(
                "Hard-to-read genes read",
                format!(
                    "{} of {}",
                    text(val(d, "read")),
                    text(val(d, "catalog_genes"))
                ),
            );
            out.fact(
                "Not read, no other step calls them",
                genes.len().to_string(),
            );
            out.fact(
                "Not read here, called by a step of their own",
                array(d, "not_read_owned")?.len().to_string(),
            );
            let rows = genes
                .iter()
                .map(|g| {
                    object_row(g)?;
                    let share = num(get(g, "share"))
                        .map_or_else(|| ".".into(), |n| format!("{}%", fixed(n * 100.0, 0)));
                    let reasons = array(g, "reasons")?
                        .iter()
                        .map(text)
                        .collect::<Vec<_>>()
                        .join(", ");
                    Ok(vec![
                        cell(val(g, "gene")),
                        cell(val(g, "status")),
                        share,
                        cell(&s(&reasons)),
                    ])
                })
                .collect::<Read<_>>()?;
            out.table(&["Gene", "Status", "Bases read", "Why"], rows);
        }
        "variants" => {
            out.facts(
                summary.section("variants").data,
                &[
                    ("Total variants", "total"),
                    ("PASS variants", "pass"),
                    ("SNPs", "snps"),
                    ("Indels", "indels"),
                ],
            );
            out.links = vec!["vep"];
        }
        _ => unreachable!(),
    }
    Ok(out)
}

pub(super) fn findings(summary: &Summary, found: &mut Vec<Finding>) {
    let d = data(summary, "sample_qc");
    if string(d, "sex_check") == "mismatch" {
        found.push(Finding::new(Level::Trust, "qc", &["sex"], "Declared and measured sex differ",
            "The sex the reads show differs from the one declared: a swapped or mislabelled sample, or a real sex-chromosome difference.")
            .evidence_fields(d, &[("Declared", "declared_sex"), ("From the reads", "inferred_sex")]));
    }
    if string(d, "contamination") == "warn" {
        add_row(found, "qc", Level::Review, || {
            Ok(Some(
                Finding::new(
                    Level::Trust,
                    "qc",
                    &["contamination"],
                    "Another person's DNA may be mixed in",
                    qc::freemix_text(&summary.section("sample_qc"))?,
                )
                .evidence("FREEMIX", val(d, "freemix")),
            ))
        });
    }
    let same = val(d, "same_person_as");
    if truthy(same) {
        add_row(found, "qc", Level::Trust, || {
            let identity = qc::identity_uncertainty(&summary.section("sample_qc"))?;
            let (label, detail) = match identity {
                Some(detail) => (qc::POSSIBLE_MATCH, detail),
                None => (
                    "The same person as",
                    "Two runs of this server hold the same genome.".to_owned(),
                ),
            };
            Ok(Some(
                Finding::new(
                    Level::Trust,
                    "qc",
                    &["same-person"],
                    format!("{label} {}", text(same)),
                    detail,
                )
                .evidence("Matches", same),
            ))
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn live_qc_facts_qualify_uncertain_evidence() {
        let summary = Summary::parse(
            r#"{"sections":{"sample_qc":{"state":"ok","data":{
                "contamination":"ok","freemix":"0.0010",
                "contamination_confidence":"uncertain",
                "contamination_confidence_reason":"FREEMIX rests on 300 markers",
                "identity":"uncertain","identity_reason":"only 300 shared sites",
                "same_person_as":"S2 (0.99)"
            }}}}"#,
        )
        .unwrap();
        let shown = card("qc", &summary).unwrap();
        assert_eq!(
            shown.facts[2..],
            [
                strings([
                    "Contamination (VerifyBamID2)",
                    "uncertain (FREEMIX 0.0010): FREEMIX rests on 300 markers",
                ]),
                strings(["Sample identity", "uncertain (only 300 shared sites)"]),
                strings(["Possible match to", "S2 (0.99)"]),
            ]
        );
        assert!(shown.attention);
    }
}
