use super::*;
use og_render::sections::traits::{LOW_GQ_LABEL, low_gq_count};
use og_render::sections::{apoe, bloodtype, consumer, g2, prs, telomere, viral};

pub(super) fn card(id: &str, summary: &Summary) -> Read<Adapter> {
    let mut out = Adapter::default();
    match id {
        "traits-panel" => {
            let d = summary.section("traits").data;
            let rows = array(d, "rows")?
                .iter()
                .filter(|r| {
                    matches!(string(r, "status"), "called" | "low_gq")
                        && !string(r, "gene").eq_ignore_ascii_case("APOE")
                })
                .map(|r| {
                    if string(r, "status") == "low_gq" {
                        return vec![
                            cell(val(r, "gene")),
                            cell(val(r, "trait")),
                            ".".to_owned(),
                            text_or(val(r, "reading"), "."),
                        ];
                    }
                    let genotype = text_or(val(r, "genotype"), ".");
                    let ancestry = if string(r, "ancestry") == "yes" {
                        " (ancestry)"
                    } else {
                        ""
                    };
                    let reading = text_or(val(r, "reading"), ".");
                    let note = val(r, "note");
                    let reading = if truthy(note) {
                        format!("{reading}. {}", text(note))
                    } else {
                        reading
                    };
                    vec![
                        cell(val(r, "gene")),
                        cell(val(r, "trait")),
                        format!("{genotype}{ancestry}"),
                        reading,
                    ]
                })
                .collect();
            out.table(&["Gene", "Trait", "Genotype", "Reading"], rows);
            out.facts(
                d,
                &[
                    ("Rows read", "called"),
                    ("No genotype at the position", "not_covered"),
                ],
            );
            if let Some(count) = low_gq_count(d) {
                out.fact(LOW_GQ_LABEL, count.to_string());
            }
            for r in array(d, "rules")? {
                out.lines.push(format!(
                    "{} ({} {}): {}. {}",
                    text(val(r, "trait")),
                    text(val(r, "gene")),
                    text_or(val(r, "rsids"), "").replace(',', ", "),
                    g2::rule_result(r)?,
                    text(val(r, "reading"))
                ));
                if let Some(g2) = g2::g2_text(r)? {
                    out.lines.push(g2);
                }
            }
            let apoe = val(d, "apoe");
            if truthy(apoe) {
                let lines = super::rows(apoe)
                    .iter()
                    .map(|r| {
                        if matches!(r, Value::Object(_)) {
                            og_render::sections::apoe::text_lines(std::slice::from_ref(r))
                                .map(|mut lines| lines.pop().unwrap_or_default().trim().to_owned())
                        } else {
                            Ok(cell(r))
                        }
                    })
                    .collect::<Read<Vec<_>>>()?;
                out.extra(
                    "hidden",
                    obj([
                        ("summary", s(apoe::APOE_WARNING)),
                        ("text", s(&plain(&lines.join("; ")))),
                    ]),
                );
            }
        }
        "prs" => {
            let d = summary.section("prs").data;
            let scores = array(d, "scores")?;
            if !truthy(val(d, "adjusted")) {
                out.fact("Scores", scores.len().to_string());
                out.lines.push("Raw scores only: without the ancestry panel no percentile is shown (docs/25-prs.md).".into());
            } else {
                let rows = scores
                    .iter()
                    .map(|r| {
                        let effect = prs::effect(r)?.0;
                        Ok(vec![
                            cell(val(r, "condition")),
                            cell(val(r, "pgs_id")),
                            cell(&s(&prs::percentile(r)?)),
                            format!("{}/{}", text(val(r, "matched")), text(val(r, "total"))),
                            cell(&s(&effect)),
                        ])
                    })
                    .collect::<Read<_>>()?;
                out.table(
                    &[
                        "Condition",
                        "Score",
                        "Percentile",
                        "Variants matched",
                        "Effect",
                    ],
                    rows,
                );
                out.lines.push(prs::note(d)?);
            }
        }
        "bloodtype" => {
            let d = summary.section("bloodtype").data;
            out.facts(d, &[("ABO", "abo"), ("RhD", "rhd")]);
            let systems = array(d, "systems")?;
            out.fields_table(
                &["System", "Predicted phenotype", "Genotype", "Confidence"],
                systems,
                &["system", "phenotype", "genotype", "confidence"],
            );
            for row in systems {
                if let Some((system, note)) = bloodtype::system_note(row)? {
                    out.lines.push(format!("{system}: {note}"));
                }
            }
        }
        "telomere" => {
            let sec = summary.section("telomere");
            out.fact(
                "Telomere content",
                format!(
                    "{} ({})",
                    text(val(sec.data, "tel_content")),
                    telomere::TELOMERE_UNIT
                ),
            );
            out.lines = telomere::lines(&sec)?;
            out.attention = truthy(val(sec.data, "low_reads"));
        }
        "viral" => {
            let d = summary.section("viral").data;
            let found = array(d, "found")?;
            let hhv6 = cell(val(d, "hhv6"));
            out.fact("HHV-6", viral::hhv6_reading(&hhv6, found));
            out.fact("Genomes searched", cell(val(d, "genomes")));
            out.fields_table(
                &["Virus", "Reads", "Breadth", "Copies per cell", "Call"],
                found,
                &["virus", "reads", "breadth", "copies_per_cell", "call"],
            );
            out.lines = viral::threshold_notes(found);
        }
        "consumer-export" => {
            let sec = summary.section("consumer_export");
            out.facts(
                sec.data,
                &[
                    ("Sites", "sites"),
                    ("Called", "called"),
                    ("No call", "no_call"),
                ],
            );
            out.lines = consumer::lines(&sec)?;
        }
        _ => unreachable!(),
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const GQ_READING: &str = "A low-quality call at this position (genotype quality 12, below 20): not read; check the reads.";
    const FILTERED_READING: &str = "A low-quality call at this position (the caller filtered the call: LowQual): not read; check the reads.";
    const THRESHOLD_NOTE: &str = "7 counted reads, under the 10 a call needs: too few to tell the virus from contamination at this depth.";

    fn summary() -> Summary {
        Summary::parse(&format!(
            r#"{{"sections": {{
                "traits": {{"data": {{"called": 0, "not_covered": 0, "low_gq": 2, "rows": [
                    {{"gene": "GENE1", "trait": "Example GQ trait", "genotype": "A/G",
                     "status": "low_gq", "reading": "{GQ_READING}", "note": "Study note"}},
                    {{"gene": "GENE2", "trait": "Example filtered trait", "genotype": "C/T",
                     "status": "low_gq", "reading": "{FILTERED_READING}", "note": "Study note"}}
                ]}}}},
                "bloodtype": {{"data": {{"systems": [
                    {{"system": "ABO", "note": "Rare ABO alleles are outside this panel"}},
                    {{"system": "RhD", "note": "RhD needs read depth"}},
                    {{"system": "Duffy", "note": "Duffy phase was unknown"}},
                    {{"system": "Kell", "note": "Kell study note"}}
                ]}}}},
                "viral": {{"data": {{"hhv6": "no integrated HHV-6", "genomes": 24, "found": [
                    {{"virus": "HHV-6B", "call": "below threshold", "reads": "7",
                     "breadth": "0.0041", "copies_per_cell": ".", "note": "{THRESHOLD_NOTE}"}}
                ]}}}}
            }}}}"#
        ))
        .unwrap()
    }

    #[test]
    fn counts_low_quality_trait_calls() {
        let traits = card("traits-panel", &summary()).unwrap();
        assert!(
            traits
                .facts
                .contains(&strings(["Low-quality genotype calls", "2"]))
        );
    }

    #[test]
    fn shows_unread_trait_rows_with_the_producers_reasons() {
        let traits = card("traits-panel", &summary()).unwrap();
        let rows = val(traits.table.as_ref().unwrap(), "rows");
        assert_eq!(
            *rows,
            Value::Array(vec![
                strings(["GENE1", "Example GQ trait", ".", GQ_READING]),
                strings(["GENE2", "Example filtered trait", ".", FILTERED_READING]),
            ])
        );
    }

    #[test]
    fn shows_only_the_abo_rhd_and_duffy_notes() {
        let bloodtype = card("bloodtype", &summary()).unwrap();
        assert_eq!(
            bloodtype.lines,
            [
                "ABO: Rare ABO alleles are outside this panel",
                "RhD: RhD needs read depth",
                "Duffy: Duffy phase was unknown",
            ]
        );
    }

    #[test]
    fn uses_the_hhv6_reading_for_reads_below_threshold() {
        let viral = card("viral", &summary()).unwrap();
        assert!(
            viral
                .facts
                .contains(&strings(["HHV-6", "no integration detected"]))
        );
    }

    #[test]
    fn explains_why_viral_reads_are_below_threshold() {
        let viral = card("viral", &summary()).unwrap();
        assert_eq!(viral.lines, [format!("HHV-6B: {THRESHOLD_NOTE}")]);
    }
}
