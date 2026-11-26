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
                let (firm, possible) = qc::same_person(&sample)?;
                if let Some(ref p) = firm {
                    out.fact("The same person as", p);
                }
                if let Some(ref p) = possible {
                    out.fact(qc::POSSIBLE_MATCH, p);
                }
                out.attention = sex_check == "mismatch"
                    || contam == "warn"
                    || contam_uncertain
                    || identity.is_some()
                    || firm.is_some()
                    || possible.is_some();
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
    if truthy(val(d, "same_person_as")) {
        add_row(found, "qc", Level::Trust, || {
            let sample = summary.section("sample_qc");
            let (firm, _) = qc::same_person(&sample)?;
            let Some(firm) = firm else {
                return Ok(None);
            };
            Ok(Some(
                Finding::new(
                    Level::Trust,
                    "qc",
                    &["same-person"],
                    format!("The same person as {firm}"),
                    "Two runs of this server hold the same genome.",
                )
                .evidence("Matches", &Value::Str(firm)),
            ))
        });
        add_row(found, "qc", Level::Trust, || {
            let sample = summary.section("sample_qc");
            let (_, possible) = qc::same_person(&sample)?;
            let Some(possible) = possible else {
                return Ok(None);
            };
            // With identity uncertain, its reason; else the pair shares
            // too few sites.
            let detail = qc::identity_uncertainty(&sample)?.unwrap_or_else(|| {
                "The two samples share too few genotyped sites to call them the same person."
                    .to_owned()
            });
            Ok(Some(
                Finding::new(
                    Level::Trust,
                    "qc",
                    &["possible-same-person"],
                    format!("{} {possible}", qc::POSSIBLE_MATCH),
                    detail,
                )
                .evidence("Matches", &Value::Str(possible)),
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

        let summary = Summary::parse(
            r#"{"sections":{"sample_qc":{"state":"ok","data":{
                "contamination":"ok","freemix":"0.0010",
                "contamination_confidence":"uncertain",
                "contamination_confidence_reason":"FREEMIX rests on 300 markers",
                "identity":"checked"
            }}}}"#,
        )
        .unwrap();
        assert!(card("qc", &summary).unwrap().attention);
    }

    /// Step 33's own table for S with these pairs, as a summary with only
    /// the `sample_qc` identity rows.
    fn from_step_33(tag: &str, pairs: &str) -> Summary {
        use og_summary::qc::{FREEMIX_WARN, Inputs, table};
        let d = std::env::temp_dir().join(format!("og-live-thin-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        let (sp, pp) = (d.join("s.samples.tsv"), d.join("s.pairs.tsv"));
        std::fs::write(
            &sp,
            "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\t\
original_pedigree_sex\tgt_depth_mean\tn_hom_ref\tn_het\tn_hom_alt\tp_middling_ab\t\
X_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\tY_depth_mean\tY_n\n\
S\tS\t-9\t-9\t1\t-9\t-9\t30.0\t5000\t6000\t4000\t0.010\t15.0\t39\t0\t0\t39\t14.0\t5\n",
        )
        .unwrap();
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
            .filter(|(k, _)| k.contains("identity") || k.contains("same_person"))
            .map(|(k, v)| format!("\"{k}\": \"{v}\""))
            .collect::<Vec<_>>()
            .join(", ");
        Summary::parse(&format!(
            r#"{{"sections":{{"sample_qc":{{"state":"ok","data":{{{data}}}}}}}}}"#
        ))
        .unwrap()
    }

    fn titles(summary: &Summary) -> Vec<(String, String)> {
        let mut found = Vec::new();
        findings(summary, &mut found);
        found.into_iter().map(|f| (f.id, f.title)).collect()
    }

    #[test]
    fn live_qc_shows_a_thin_pair_as_a_possible_match_beside_a_firm_one() {
        // T: the same person on 12,000 shared sites; U: 0.97 on 300; V:
        // unrelated on 15,000. Identity is `checked`.
        let summary = from_step_33(
            "mixed",
            "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\n\
             S\tT\t0.99\t0\t300\t12000\nU\tS\t0.97\t0\t300\t300\nS\tV\t0.1\t9\t300\t15000\n",
        );
        let shown = card("qc", &summary).unwrap();
        let tail = &shown.facts[shown.facts.len() - 2..];
        assert_eq!(
            tail,
            [
                strings(["The same person as", "T (0.99)"]),
                strings([
                    "Possible match to",
                    "U (0.97) (on fewer than 1,000 shared sites)"
                ]),
            ]
        );
        assert!(shown.attention);
        let t = titles(&summary);
        assert_eq!(
            t.iter()
                .map(|(_, title)| title.as_str())
                .collect::<Vec<_>>(),
            [
                "The same person as T (0.99)",
                "Possible match to U (0.97) (on fewer than 1,000 shared sites)"
            ],
            "{t:?}"
        );
        assert_ne!(t[0].0, t[1].0, "two findings, two ids: {t:?}");

        // Both pairs on enough sites: one firm finding, no possible one.
        let summary = from_step_33(
            "firm",
            "#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\n\
             S\tT\t0.99\t0\t300\t12000\nU\tS\t0.95\t0\t300\t1000\n",
        );
        let shown = card("qc", &summary).unwrap();
        assert_eq!(
            shown.facts.last(),
            Some(&strings(["The same person as", "T (0.99), U (0.95)"]))
        );
        let t = titles(&summary);
        assert_eq!(t.len(), 1, "{t:?}");
        assert_eq!(t[0].1, "The same person as T (0.99), U (0.95)");
    }
}
