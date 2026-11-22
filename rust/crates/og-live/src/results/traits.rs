use super::*;
use og_render::sections::traits::{LOW_GQ_LABEL, LOW_GQ_READING};
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
                        let reason = string(r, "note");
                        let reading = if matches!(reason.trim(), "" | ".") {
                            LOW_GQ_READING.to_owned()
                        } else {
                            format!("{LOW_GQ_READING}. {reason}")
                        };
                        return vec![
                            cell(val(r, "gene")),
                            cell(val(r, "trait")),
                            ".".to_owned(),
                            reading,
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
            if !matches!(val(d, "low_gq"), Value::Null) {
                out.fact(LOW_GQ_LABEL, cell(val(d, "low_gq")));
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
            let sec = summary.section("viral");
            let d = sec.data;
            let hhv6_value = val(d, "hhv6");
            out.fact(
                "HHV-6",
                if matches!(hhv6_value, Value::Null) {
                    cell(hhv6_value)
                } else {
                    viral::hhv6(&sec)?
                },
            );
            out.fact("Genomes searched", cell(val(d, "genomes")));
            out.fields_table(
                &["Virus", "Reads", "Breadth", "Copies per cell", "Call"],
                array(d, "found")?,
                &["virus", "reads", "breadth", "copies_per_cell", "call"],
            );
            out.lines = viral::threshold_notes(&sec)?;
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
