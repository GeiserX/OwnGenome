use super::clinvar::text_or_null;
use super::*;
use crate::widgets;
use og_render::sections::{acmg_sf, carrier, clinvar::D, comphet, cpsr, parascopy, smn1, two_hits};

fn sf_rows(summary: &Summary, sf: &Value) -> Read<Vec<acmg_sf::Row>> {
    if !truthy(sf) || truthy(val(sf, "unavailable")) {
        return Ok(Vec::new());
    }
    let consensus = truthy(val(data(summary, "annotsv"), "consensus"));
    acmg_sf::rows(&D::new(sf, "clinical.acmg_sf"), consensus)
}

fn carrier_rows(carrier: &Value) -> impl Iterator<Item = &Value> {
    rows(val(carrier, "carriers"))
        .iter()
        .chain(rows(val(carrier, "hom_or_hemizygous")))
}

fn two_hit_genes(carrier: &Value) -> impl Iterator<Item = &Value> {
    rows(val(val(carrier, "two_hits"), "genes"))
        .iter()
        .map(|g| {
            if matches!(g, Value::Str(_)) {
                g
            } else {
                val(g, "gene")
            }
        })
        .filter(|g| matches!(g, Value::Str(_)))
}

fn smn1_pair(summary: &Summary) -> &Value {
    val(data(summary, "paralogs"), "smn1")
}

pub(super) fn card(id: &str, summary: &Summary) -> Read<Adapter> {
    let mut out = Adapter::default();
    match id {
        "clinical-filter" => {
            let d = summary.section("clinical").data;
            let sf = val(d, "acmg_sf");
            let bi = val(d, "by_impact");
            out.fact(
                "Variants after the clinical filter",
                cell(val(d, "variants")),
            );
            out.fact(
                "By impact (high / moderate / low)",
                ["HIGH", "MODERATE", "LOW"]
                    .map(|k| default(bi, k, "0"))
                    .join(" / "),
            );
            out.fact("Genes", cell(val(d, "genes")));
            let record = D::new(sf, "clinical.acmg_sf");
            if truthy(val(sf, "unavailable")) {
                out.lines = acmg_sf::apart_lines(&record)?;
                out.attention = true;
            } else {
                let rows = sf_rows(summary, sf)?
                    .into_iter()
                    .map(|(gene, variant, kind, rule)| {
                        [gene, variant, kind, rule]
                            .into_iter()
                            .map(|t| cell(&s(&t)))
                            .collect()
                    })
                    .collect();
                out.table(&["Gene", "Variant", "Kind", "Rule"], rows);
                let high = if truthy(val(sf, "high_impact_checked")) {
                    format!(
                        "{}{}",
                        acmg_sf::high_rare(&record)?,
                        acmg_sf::high_apart(&record)?
                    )
                } else {
                    "not checked".into()
                };
                let version = text_or(val(sf, "version"), "ACMG SF");
                out.fact(
                    format!("{version}: ClinVar P/LP in its genes"),
                    list(val(sf, "clinvar_hits")).len().to_string(),
                );
                out.fact("High-impact rare variants in its genes", high);
                out.attention = ["clinvar_hits", "high_impact", "structural"]
                    .iter()
                    .any(|k| truthy(val(sf, k)));
            }
        }
        "cpsr" => {
            let sec = summary.section("cpsr");
            let d = sec.data;
            out.fact(
                "CPSR report",
                if truthy(val(d, "html_report")) {
                    "Ran"
                } else {
                    "No HTML"
                },
            );
            for (k, v) in fields(val(d, "classification")) {
                out.fact(k.replace('_', " "), cell(v));
            }
            let rows = array(d, "plp")?
                .iter()
                .map(|x| {
                    let mut r = row(x, &["gene", "variant", "final"]);
                    r.push(text_or(val(x, "clinvar"), "."));
                    r.push(
                        if truthy(val(x, "low_confidence")) {
                            "low confidence"
                        } else {
                            "."
                        }
                        .into(),
                    );
                    r
                })
                .collect();
            out.table(&["Gene", "Variant", "Final class", "ClinVar", "Call"], rows);
            let dbs = summary.get("databases").unwrap_or(&Value::Null);
            out.lines = cpsr::lines(&D::of(&sec), &D::new(dbs, "databases"))?;
            out.links = vec!["cpsr"];
        }
        "slivar" => {
            let sec = summary.section("slivar");
            let d = sec.data;
            out.facts(
                d,
                &[
                    ("Prioritized variants", "prioritized"),
                    ("Compound-het candidate variants", "compound_het_variants"),
                    ("Compound-het candidate genes", "compound_het_genes"),
                ],
            );
            if get(d, "compound_het_phase").is_some() {
                out.lines = comphet::phase_text(&D::of(&sec))?;
            }
        }
        "exomiser" => {
            let d = summary.section("exomiser").data;
            let terms = array(d, "terms")?
                .iter()
                .map(|t| format!("{} ({})", text(val(t, "label")), text(val(t, "id"))))
                .collect::<Vec<_>>()
                .join(", ");
            out.fact("Symptoms searched", text_or(&s(&terms), "none"));
            let mut hits = Vec::new();
            for g in array(d, "genes")? {
                object_row(g)?;
                let variants = array(g, "variants")?;
                let entries = if variants.is_empty() {
                    &[Value::Null]
                } else {
                    variants
                };
                for v in entries {
                    let mut r = row(g, &["rank", "gene", "moi", "exomiser_score"]);
                    r.push(if truthy(val(v, "hgvs")) {
                        cell(val(v, "hgvs"))
                    } else {
                        cell(val(v, "variant"))
                    });
                    r.extend(row(v, &["genotype", "acmg", "disease"]));
                    hits.push(r);
                }
            }
            out.table(
                &[
                    "Rank",
                    "Gene",
                    "Inheritance",
                    "Score",
                    "Variant",
                    "Genotype",
                    "ACMG",
                    "Disease",
                ],
                hits,
            );
        }
        "carrier-list" => {
            let sec = summary.section("carrier");
            let d = sec.data;
            if truthy(val(d, "unavailable")) {
                out.lines = carrier::lines(&D::of(&sec))?;
                out.attention = true;
            } else {
                let calls = array(d, "carriers")?
                    .iter()
                    .chain(array(d, "hom_or_hemizygous")?);
                let rows: Vec<_> = calls
                    .map(|c| {
                        row(
                            c,
                            &[
                                "gene",
                                "variant",
                                "genotype",
                                "significance",
                                "stars",
                                "inheritance",
                                "condition",
                            ],
                        )
                    })
                    .collect();
                out.fact("Genes on the carrier list", cell(val(d, "genes_on_list")));
                out.fact("Carrier findings", rows.len().to_string());
                let two = val(d, "two_hits");
                let genes = array(two, "genes")?;
                let genes_text = genes
                    .iter()
                    .map(|g| match g {
                        Value::Str(t) => t.clone(),
                        _ => default(g, "gene", "."),
                    })
                    .collect::<Vec<_>>()
                    .join(", ");
                out.fact("Genes with two hits", text_or(&s(&genes_text), "none"));
                out.table(
                    &[
                        "Gene",
                        "Variant",
                        "Genotype",
                        "Significance",
                        "Stars",
                        "Inheritance",
                        "Condition",
                    ],
                    rows,
                );
                let two = D::new(two, "carrier.two_hits");
                out.lines = two_hits::lines(Some(&two))?;
                out.attention = !genes.is_empty();
            }
        }
        "couple" => {
            let d = summary.section("couple").data;
            out.fact("Partner", cell(val(d, "partner")));
            let shared = array(d, "shared")?;
            let x_linked = array(d, "x_linked")?;
            out.fact("Genes where both carry a change", shared.len().to_string());
            out.fact("X-linked findings", x_linked.len().to_string());
            out.fact(
                "Genes where only one partner carries a change",
                format!(
                    "{}, unless the other partner carries something the screen does not see",
                    default(d, "one_partner", "0")
                ),
            );
            let rows = [("Both carry", shared), ("X-linked", x_linked)]
                .into_iter()
                .flat_map(|(label, items)| {
                    items.iter().map(move |x| {
                        let mut r = vec![label.to_owned()];
                        r.extend(row(
                            x,
                            &["gene", "condition", "risk", "own", "partner_findings"],
                        ));
                        r
                    })
                })
                .collect();
            out.table(
                &[
                    "Result",
                    "Gene",
                    "Condition",
                    "Chance for a child",
                    "This person",
                    "Partner",
                ],
                rows,
            );
        }
        "paralogs" => {
            let parascopy = summary.section("parascopy");
            let pc = data(summary, "parascopy");
            if parascopy.state.ran() {
                out.lines = parascopy::text_lines(&parascopy)?;
                if truthy(val(pc, "smn1_marker")) {
                    out.fact("SMN1 copies", default(val(pc, "smn1"), "text", "."));
                    out.fact(
                        "SMN1 g.27134T>G silent-carrier marker",
                        smn1::marker_state(val(pc, "smn1_marker"))?,
                    );
                }
            }
            let p = data(summary, "paralogs");
            if string(smn1_pair(summary), "agree") == "disagree"
                && let Some(verdict) = smn1::smn1_text(smn1_pair(summary))?
            {
                out.lines.push(format!("SMN1 by both tools: {verdict}"));
            }
            let genes = array(p, "genes")?;
            out.fields_table(
                &["Gene", "Copies", "Status", "Detail"],
                genes,
                &["gene", "copy_number", "status", "detail"],
            );
            out.extra("widget", widgets::copies(genes, pc));
        }
        _ => unreachable!(),
    }
    Ok(out)
}

pub(super) fn clinical_findings<'a>(
    summary: &'a Summary,
    found: &mut Vec<Finding>,
) -> Vec<&'a Value> {
    let sf = val(data(summary, "clinical"), "acmg_sf");
    let car = data(summary, "carrier");
    let sf_excluded = rows(val(sf, "clinvar_hits"))
        .iter()
        .filter(|x| matches!(x, Value::Object(_)))
        .map(|x| val(x, "gene"));
    let carrier_excluded = carrier_rows(car)
        .filter(|x| matches!(x, Value::Object(_)))
        .map(|x| val(x, "gene"));
    let two_hit: BTreeSet<_> = two_hit_genes(car)
        .map(str_val)
        .filter(|g| !g.is_empty())
        .collect();
    let sf_record = D::new(sf, "clinical.acmg_sf");
    let sf_list = text_or(val(sf, "version"), "ACMG SF");
    match sf_rows(summary, sf).and_then(|rows| {
        rows.into_iter()
            .map(|r| {
                let condition = acmg_sf::listed_condition(&sf_record, &r.0)?;
                Ok((acmg_sf::lede(&r.0, condition.as_deref()), r))
            })
            .collect::<Read<Vec<_>>>()
    }) {
        Ok(rows) => {
            for (lede, (gene, variant, genotype, why)) in rows {
                found.push(
                    Finding::new(
                        Level::Act,
                        "clinical-filter",
                        &["sf", &gene, &variant],
                        format!("Secondary finding in {gene}"),
                        lede,
                    )
                    .gene(&s(&gene))
                    .detail(format!(
                        "{gene} {variant} {}",
                        if genotype == "." { "" } else { &genotype }
                    ))
                    .evidence("Variant", &s(&variant))
                    .evidence("Genotype", &s(&genotype))
                    .evidence("Reporting list", &s(&sf_list))
                    .evidence("Why it is reported", &s(&why)),
                );
            }
        }
        Err(err) => found.push(unreadable(
            "clinical-filter",
            Level::Act,
            found.len(),
            error_kind(&err),
        )),
    }
    for x in rows(val(data(summary, "cpsr"), "plp")) {
        add_row(found, "cpsr", Level::Act, || {
            object_row(x)?;
            let gene = val(x, "gene");
            let variant = val(x, "variant");
            Ok(Some(Finding::new(Level::Act, "cpsr", &["cpsr", &text_or_null(gene), &text_or_null(variant)],
                format!("Cancer predisposition ({})", text(gene)),
                format!("CPSR grades this change in {} as {}: it can raise the risk of some cancers. It is not a diagnosis of cancer.",
                    text(gene), text(val(x, "final")))).gene(gene)
                .detail(format!("{} {}", text(gene), text(variant)))
                .evidence_fields(x, &[("Class", "final"), ("ClinVar", "clinvar")])
                .evidence("Call", &if truthy(val(x, "low_confidence")) { s("low confidence") } else { Value::Null })))
        });
    }
    for gene in two_hit {
        found.push(Finding::new(Level::Act, "carrier-list", &["two-hits", gene], format!("Two changes in {gene}"),
            format!("Two changes were found in {gene}, a gene where it takes two to cause the condition. Whether they sit on different copies decides whether it matters."))
            .gene(&s(gene)));
    }
    sf_excluded
        .chain(carrier_excluded)
        .chain(two_hit_genes(car))
        .collect()
}

pub(super) fn family_findings(summary: &Summary, found: &mut Vec<Finding>) {
    let car = data(summary, "carrier");
    let paired: Vec<_> = two_hit_genes(car).collect();
    for c in carrier_rows(car) {
        add_row(found, "carrier-list", Level::Act, || {
            object_row(c)?;
            let gene = val(c, "gene");
            if paired.contains(&gene) {
                return Ok(None);
            }
            let g = text(gene);
            let variant = val(c, "variant");
            let gt = text_or(val(c, "genotype"), "");
            let condition = val(c, "condition");
            let key = ["carrier", &text_or_null(gene), &text_or_null(variant)];
            let evidence = [
                ("Significance", "significance"),
                ("Stars", "stars"),
                ("Inheritance", "inheritance"),
            ];
            let finding = if matches!(gt.as_str(), "hom" | "hemizygous") {
                let (title, copies) = if gt == "hom" {
                    ("two changed copies", "two changed copies")
                } else {
                    ("the only copy changed", "no unchanged copy")
                };
                Finding::new(Level::Act, "carrier-list", &key, format!("{}: {title}", text_or(condition, &g)),
                    format!("This genome has {copies} of {g}, a gene of {} ({}). That is not the pattern of a carrier: a doctor or genetic counsellor reads it.",
                        text_or(condition, "a condition on the carrier list"), text_or(val(c, "inheritance"), "inheritance not given")))
                    .detail(format!("{g} {} {gt}", text_or(variant, ""))).evidence("Genotype", &s(&gt))
            } else {
                let title = if truthy(condition) {
                    format!("{} carrier", text(condition))
                } else {
                    format!("Carrier in {g}")
                };
                Finding::new(Level::Carrier, "carrier-list", &key, title,
                    format!("A carrier has one changed copy of a gene where it takes two to cause the condition. This genome has one changed copy of {g}: a carrier, not affected."))
                    .detail(format!("{g} {}", text_or(variant, "")))
            };
            // The call's mark against the hits file's floor, and the pair word
            // of two heterozygous hits in the gene with its phase.
            let call = if truthy(val(c, "low_confidence")) {
                s("low confidence: check the reads")
            } else {
                Value::Null
            };
            let pair = match (text_or(val(c, "pair"), ""), text_or(val(c, "phase"), "")) {
                (p, ph) if p.is_empty() => ph,
                (p, ph) if ph.is_empty() => p,
                (p, ph) => format!("{p}: {ph}"),
            };
            Ok(Some(
                finding
                    .gene(gene)
                    .evidence_fields(c, &evidence)
                    .evidence("Call", &call)
                    .evidence("Second hit in the gene", &s(&pair)),
            ))
        });
    }
    let cp = data(summary, "couple");
    for key in ["shared", "x_linked"] {
        for x in rows(val(cp, key)) {
            add_row(found, "couple", Level::Review, || {
                object_row(x)?;
                let gene = val(x, "gene");
                let g = text(gene);
                let label = if key == "shared" {
                    format!("Both partners carry {g}")
                } else {
                    format!("X-linked finding in {g}")
                };
                Ok(Some(
                    Finding::new(
                        Level::Carrier,
                        "couple",
                        &["couple", key, &text_or_null(gene)],
                        label,
                        format!(
                            "{}: {} for each pregnancy, with {}.",
                            text_or(val(x, "condition"), &g),
                            text_or(val(x, "risk"), "see the couple screen"),
                            text(val(cp, "partner"))
                        ),
                    )
                    .gene(gene)
                    .evidence("Chance for a child", val(x, "risk"))
                    .evidence("Partner", val(cp, "partner")),
                ))
            });
        }
    }
    add_row(found, "paralogs", Level::Review, || {
        let pair = smn1_pair(summary);
        if string(pair, "agree") == "disagree" && smn1::smn1_text(pair)?.is_some() {
            return Ok(None);
        }
        let mut called = None;
        for g in array(data(summary, "paralogs"), "genes")? {
            object_row(g)?;
            if string(g, "gene") == "SMN1" && string(g, "status") == "called" && called.is_none() {
                called = Some(val(g, "copy_number"));
            }
        }
        let (one, by) = if let Some(cn) = called {
            (text(cn) == "1", "Paraphase")
        } else {
            (
                num(Some(val(val(data(summary, "parascopy"), "smn1"), "cn"))) == Some(1.0),
                "Parascopy",
            )
        };
        Ok(one.then(|| Finding::new(Level::Carrier, "paralogs", &["smn1"], "Spinal muscular atrophy carrier (SMN1)",
            "One copy of SMN1 where most people have two: a carrier of spinal muscular atrophy, not affected.")
            .gene(&s("SMN1")).evidence("SMN1 copies", &s("1")).evidence("Called by", &s(by))))
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(sections: &str) -> Summary {
        Summary::parse(&format!(r#"{{"sections": {{{sections}}}}}"#)).unwrap()
    }

    /// The evidence of a finding as "label: value" lines.
    fn evidence_of(f: &Finding) -> Vec<String> {
        f.evidence
            .iter()
            .map(|e| list(e).iter().map(text).collect::<Vec<_>>().join(": "))
            .collect()
    }

    /// A ClinVar hit in BRCA1 on the SF card: one ClinVar P/LP hit whose
    /// condition step 6 copied, so it carries the disease ACMG lists.
    const BRCA1: &str = r#""clinical": {"state": "ok", "data": {"variants": 1, "acmg_sf": {
        "version": "ACMG SF v3.3", "genes_on_list": 84, "high_impact": [], "high_impact_checked": true,
        "clinvar_hits": [{"gene": "BRCA1", "variant": "chr17:43045712 G>A", "genotype": "het",
          "significance": "Pathogenic", "stars": 3, "low_confidence": false,
          "clinvar_condition": "Hereditary breast ovarian cancer syndrome", "clinvar_condition_ids": [],
          "acmg_phenotype": "hereditary breast and ovarian cancer", "condition_not_acmg": false,
          "rule": "All P and LP", "rule_note": null, "severity": null}]}}}"#;

    /// The SF finding opens in plain words: the gene is one doctors report
    /// even unasked, the disease it is listed for, and that this is no
    /// diagnosis. The list's name and the ClinVar class stay in the
    /// evidence rows.
    #[test]
    fn the_sf_finding_says_why_in_plain_words() {
        let mut found = Vec::new();
        clinical_findings(&parse(BRCA1), &mut found);
        let f = found
            .iter()
            .find(|f| f.title == "Secondary finding in BRCA1")
            .expect("the SF finding");
        assert!(
            f.text.contains("hereditary breast and ovarian cancer"),
            "{}",
            f.text
        );
        assert!(f.text.contains("not a diagnosis"), "{}", f.text);
        assert!(f.text.contains("certified laboratory"), "{}", f.text);
        assert!(!f.text.contains("ACMG"), "{}", f.text);
        let ev = evidence_of(f);
        assert!(
            ev.contains(&"Reporting list: ACMG SF v3.3".to_owned()),
            "{ev:?}"
        );
        assert!(
            ev.iter().any(|e| e.starts_with(
                "Why it is reported: ClinVar Pathogenic; \
                ClinVar condition: Hereditary breast ovarian cancer syndrome"
            )),
            "{ev:?}"
        );
    }

    /// A carrier row below the call floor, and two heterozygous hits of one
    /// recessive gene, show the mark and the pair on the carrier finding; a
    /// row with neither gains no evidence line.
    #[test]
    fn a_carrier_finding_shows_its_call_and_its_pair() {
        let summary = parse(
            r#""carrier": {"state": "ok", "data": {"genes_on_list": 3, "set_aside": [], "hom_or_hemizygous": [],
              "carriers": [
                {"gene": "GENE1", "variant": "chr1:100 C>T", "genotype": "het", "significance": "Pathogenic",
                 "stars": 2, "inheritance": "AR", "condition": "Condition one", "low_confidence": true},
                {"gene": "GENE2", "variant": "chr2:100 C>T", "genotype": "het", "significance": "Pathogenic",
                 "stars": 2, "inheritance": "AR", "condition": "Condition two", "low_confidence": false,
                 "pair": "possible two hits", "phase": "phase unknown"},
                {"gene": "GENE3", "variant": "chr3:100 C>T", "genotype": "het", "significance": "Pathogenic",
                 "stars": 2, "inheritance": "AR", "condition": "Condition three", "low_confidence": false}]}}"#,
        );
        let mut found = Vec::new();
        family_findings(&summary, &mut found);
        let ev = |gene: &str| {
            evidence_of(
                found
                    .iter()
                    .find(|f| f.gene == s(gene))
                    .expect("the carrier finding"),
            )
        };
        assert!(ev("GENE1").contains(&"Call: low confidence: check the reads".to_owned()));
        assert!(
            ev("GENE2")
                .contains(&"Second hit in the gene: possible two hits: phase unknown".to_owned())
        );
        assert!(
            !ev("GENE3")
                .iter()
                .any(|e| e.starts_with("Call") || e.starts_with("Second hit")),
            "{:?}",
            ev("GENE3")
        );
    }
}
