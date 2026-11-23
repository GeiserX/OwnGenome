use super::clinvar::text_or_null;
use super::*;
use og_render::sections::{hla, pgx};

pub(super) fn card(id: &str, summary: &Summary) -> Read<Adapter> {
    let mut out = Adapter::default();
    match id {
        "pharmcat" => {
            let d = summary.section("pharmcat").data;
            out.fact("PharmCAT report", "Ran");
            out.fact("PharmCAT version", text_or(val(d, "version"), "unknown"));
            out.fact(
                "Its own HTML report",
                if truthy(val(d, "html_report")) {
                    "beside it"
                } else {
                    "not written"
                },
            );
            out.links = vec!["pharmcat"];
        }
        "cpic" => {
            let d = summary.section("cpic").data;
            if truthy(val(d, "parse_failed")) {
                out.lines.push(
                    "The PharmCAT report could not be parsed: see the CPIC report.".to_owned(),
                );
                out.attention = true;
            } else {
                for (label, key) in [
                    ("Genes with a changed-function phenotype", "non_normal"),
                    ("Genes with more than one possible result", "ambiguous"),
                    ("Genes not called", "not_called"),
                ] {
                    out.fact(label, get(d, key).map_or_else(|| "0".into(), cell));
                }
                if let Some(incomplete) =
                    get(d, "incomplete").filter(|v| matches!(v, Value::Int(n) if *n > 0))
                {
                    out.fact(pgx::CPIC_INCOMPLETE_LABEL, cell(incomplete));
                    out.attention = true;
                }
                let rows = array(d, "genes")?
                    .iter()
                    .map(|g| {
                        let phenotype = val(g, "phenotype");
                        vec![
                            cell(val(g, "gene")),
                            cell(val(g, "diplotype")),
                            if truthy(phenotype) {
                                format!("\"{}\"", text(phenotype))
                            } else {
                                ".".into()
                            },
                            cell(val(g, "status")),
                        ]
                    })
                    .collect();
                out.table(
                    &["Gene", "Diplotype", "Phenotype (CPIC's words)", "Status"],
                    rows,
                );
            }
            out.lines.extend(
                array(d, "warnings")?
                    .iter()
                    .map(|w| format!("WARNING: {}", text(w))),
            );
        }
        "cyp2d6" => {
            let d = summary.get("cyp2d6").unwrap_or(&Value::Null);
            let rows = fields(val(d, "calls"))
                .map(|(k, v)| {
                    let name = match k {
                        "PharmCAT" => "PharmCAT (step 36's call)",
                        other => other,
                    };
                    vec![name.to_owned(), text_or(v, "not run")]
                })
                .collect();
            out.table(&["Caller", "CYP2D6"], rows);
            let agree = val(d, "agree");
            out.fact(
                "pypgx and Cyrius agree",
                if truthy(agree) {
                    "yes"
                } else if *agree == Value::Bool(false) {
                    "no"
                } else {
                    "fewer than two usable calls"
                },
            );
            if truthy(d) {
                let verdict = pgx::cyp2d6_verdict(d)?;
                if !verdict.is_empty() {
                    out.fact("Step 36", verdict);
                }
            }
            let pypgx = summary.section("pypgx");
            if pypgx.state.ran() {
                out.lines = pgx::pypgx_lines(&pypgx)?;
            }
            out.attention = *agree == Value::Bool(false);
            if let Some(warning) = pgx::cyp2d6_weak_evidence(&pypgx)? {
                out.lines.push(warning);
                out.attention = true;
            }
        }
        "hla" => {
            let d = summary.section("hla").data;
            let dbs = summary.get("databases").unwrap_or(&Value::Null);
            out.fact("HLA database", text_or(val(dbs, "hla_database"), "unknown"));
            let rows = array(d, "loci")?
                .iter()
                .map(|l| Ok(vec![cell(val(l, "gene")), cell(&s(&hla::hla_text(l)?))]))
                .collect::<Read<_>>()?;
            out.table(&["Gene", "Alleles (T1K quality)"], rows);
            out.lines = hla::hla_derived_lines(get(d, "derived"))?;
        }
        "locityper" => {
            let d = summary.section("locityper").data;
            let loci = array(d, "loci")?;
            out.fact(
                "Loci called",
                format!("{} of {}", default(d, "called", "0"), loci.len()),
            );
            out.fact(
                "Passing the authors' filter",
                get(d, "passed").map_or_else(|| "0".into(), cell),
            );
            out.fields_table(
                &[
                    "Locus",
                    "Panel haplotypes",
                    "Quality",
                    "Filter",
                    "Primary caller",
                ],
                loci,
                &["locus", "genotype", "quality", "filter", "primary"],
            );
        }
        _ => unreachable!(),
    }
    Ok(out)
}

pub(super) fn findings(summary: &Summary, found: &mut Vec<Finding>) {
    for g in rows(val(data(summary, "cpic"), "genes")) {
        add_row(found, "cpic", Level::Review, || {
            object_row(g)?;
            let status = string(g, "status");
            if !matches!(status, "non-normal" | "ambiguous") {
                return Ok(None);
            }
            let gene = val(g, "gene");
            let name = text(gene);
            let diplotype = text(val(g, "diplotype"));
            let phenotype = text(val(g, "phenotype"));
            let title = if status == "ambiguous" {
                format!("{name}: more than one possible result")
            } else {
                format!("{name}: {phenotype}")
            };
            Ok(Some(Finding::new(Level::Medicine, "cpic", &["pgx", &text_or_null(gene)], title,
                format!("{name} {diplotype} reads as \"{phenotype}\" in CPIC's words. The drugs it bears on are in the CPIC report; tell whoever prescribes."))
                .gene(gene).detail(format!("{name} {diplotype}"))
                .evidence_fields(g, &[("Diplotype", "diplotype"), ("Phenotype", "phenotype"), ("Status", "status")])))
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn incomplete_cpic_calls_set_attention_and_are_counted_and_named() {
        let summary = Summary::parse(
            r#"{"sections":{"cpic":{"state":"ok","data":{"incomplete":1,"warnings":[],
            "genes":[{"gene":"CYP2C19","diplotype":"*1/*1","status":"incomplete",
            "phenotype":"Normal Metabolizer [incomplete: positions not read: rs4244285]"}]}}}}"#,
        )
        .expect("summary");
        let out = card("cpic", &summary).expect("CPIC card");
        assert!(out.attention);
        assert!(out.facts.contains(&strings([
            "Genes called with positions PharmCAT never read (no drug guidance)",
            "1",
        ])));
        let table = out.table.expect("CPIC rows");
        let rows = array(&table, "rows").expect("CPIC rows");
        assert_eq!(
            rows,
            &[strings([
                "CYP2C19",
                "*1/*1",
                "\"Normal Metabolizer [incomplete: positions not read: rs4244285]\"",
                "incomplete",
            ])]
        );
    }
}
