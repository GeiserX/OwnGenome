use super::clinvar::text_or_null;
use super::*;
use crate::widgets;
use og_render::sections::{
    annotsv, chip, clinvar::D, lpa, mito, mocha, sex_loss, smn1, str as repeats, vntyper,
};
use pycompat::fmt::fixed;

fn annotsv_support(summary: &Summary, x: &D<'_>) -> Read<String> {
    annotsv::support_text(x, truthy(val(data(summary, "annotsv"), "consensus")))
}

fn repeat_status(locus: &Value) -> Read<&str> {
    let status = val(locus, "status");
    if !truthy(status) {
        return Ok("");
    }
    match status {
        Value::Str(s) => Ok(s),
        _ => Err(Malformed("repeat status is not a string".into())),
    }
}

fn repeat_flags(status: &str) -> impl Iterator<Item = &str> {
    status
        .split(',')
        .filter(|x| !x.is_empty() && *x != "normal")
}

pub(super) fn card(id: &str, summary: &Summary) -> Read<Adapter> {
    let mut out = Adapter::default();
    match id {
        "structural-variants" => {
            let m = summary.section("manta").data;
            let dl = summary.section("delly").data;
            let cn = summary.section("cnvpytor").data;
            let sv = summary.section("sv_consensus").data;
            let a = summary.section("annotsv").data;
            out.fact(
                "Manta SVs (total / PASS)",
                format!("{} / {}", default(m, "total", "."), default(m, "pass", ".")),
            );
            out.fact(
                "Delly SVs (total / PASS)",
                format!(
                    "{} / {}",
                    default(dl, "total", "."),
                    default(dl, "pass", ".")
                ),
            );
            out.fact(
                "CNVpytor CNVs (significant)",
                format!(
                    "{} ({})",
                    default(cn, "total", "."),
                    default(cn, "significant", ".")
                ),
            );
            out.fact("Consensus SVs (2+ callers)", cell(val(sv, "consensus")));
            out.fact("AnnotSV class 4 and 5", cell(val(a, "count")));
            let mut rows = Vec::new();
            for x in array(a, "rows")? {
                object_row(x)?;
                let record = D::new(x, "AnnotSV row");
                let mut r = row(x, &["type", "sv", "gt", "length", "acmg_class"]);
                let genes = if matches!(val(x, "genes"), Value::Null) {
                    ".".into()
                } else {
                    annotsv::genes_text(&record)?
                };
                let flags = if record.has("exon_edge") {
                    annotsv::flags(&record)?.join("; ")
                } else {
                    ".".into()
                };
                r.extend([
                    cell(&s(&genes)),
                    cell(&s(&annotsv::support_text(
                        &record,
                        truthy(val(a, "consensus")),
                    )?)),
                    cell(&s(&flags)),
                ]);
                rows.push(r);
            }
            out.table(
                &[
                    "Type", "Where", "Genotype", "Length", "Class", "Genes", "Support", "Checks",
                ],
                rows,
            );
            out.attention = truthy(val(a, "count"));
        }
        "mei" => {
            let d = summary.section("mei").data;
            let sites = array(d, "sites")?;
            out.fact("Insertions called", cell(val(d, "calls")));
            let elements = fields(val(d, "by_element"))
                .map(|(k, v)| format!("{k} {}", text(v)))
                .collect::<Vec<_>>()
                .join(", ");
            out.fact("By element", text_or(&s(&elements), "none"));
            out.fields_table(
                &[
                    "Gene",
                    "Element",
                    "Condition",
                    "Inheritance",
                    "Known site",
                    "Depth",
                ],
                sites,
                &[
                    "gene",
                    "element",
                    "condition",
                    "inheritance",
                    "status",
                    "mean_depth",
                ],
            );
            out.attention = sites
                .iter()
                .any(|x| !matches!(string(x, "status"), "" | "not seen"));
        }
        "expansions" => {
            let d = summary.section("expansions").data;
            let mut outside = 0;
            let mut rows = Vec::new();
            let mut gauge = Vec::new();
            for l in array(d, "loci")? {
                object_row(l)?;
                let flagged = repeat_flags(repeat_status(l)?).next().is_some();
                outside += usize::from(flagged);
                let zone = val(l, "zone");
                out.attention |= flagged || truthy(zone);
                let status = repeats::status(l)?;
                let shown = match status.as_str() {
                    "" => "not graded".to_owned(),
                    "not graded" | "not sized" => status.clone(),
                    _ => format!("\"{status}\""),
                };
                let count = if matches!(val(l, "repeat_count"), Value::Str(_)) {
                    repeats::count(l)?
                } else {
                    cell(val(l, "repeat_count"))
                };
                let mut r = vec![cell(val(l, "locus")), cell(&s(&count)), shown];
                r.extend(row(l, &["normal_max", "pathologic_min", "ci", "sizing"]));
                r.push(if truthy(zone) {
                    format!("{} ({})", text(val(l, "zone_label")), text(zone))
                } else {
                    ".".into()
                });
                rows.push(r);
                gauge.push((l.clone(), status));
            }
            out.fact("Loci tested", cell(val(d, "records")));
            out.fact(
                "Outside Stranger's \"normal\"",
                if truthy(val(d, "stranger")) {
                    outside.to_string()
                } else {
                    "Stranger not run".into()
                },
            );
            for (label, key) in [
                ("In an intermediate zone", "intermediate"),
                ("Not sized from short reads", "unsized"),
            ] {
                let value = array(d, key)?
                    .iter()
                    .map(text)
                    .collect::<Vec<_>>()
                    .join(", ");
                out.fact(label, text_or(&s(&value), "none"));
            }
            out.table(
                &[
                    "Locus",
                    "Repeats",
                    "Stranger status",
                    "Normal up to",
                    "Pathologic from",
                    "95% interval",
                    "Sizing",
                    "Intermediate zone",
                ],
                rows,
            );
            out.extra("widget", widgets::str_gauge(gauge));
        }
        "vntyper" => {
            let sec = summary.section("vntyper");
            out.facts(
                sec.data,
                &[
                    ("Result", "status"),
                    ("Confidence", "confidence"),
                    ("Mean depth over the repeat", "vntr_mean_depth"),
                ],
            );
            out.lines = vntyper::lines(&sec)?;
            out.attention = matches!(string(sec.data, "status"), "positive" | "possible");
        }
        "mito" => {
            let sec = summary.section("mito");
            let d = sec.data;
            out.facts(
                d,
                &[
                    ("PASS variants on chrM", "pass"),
                    ("Heteroplasmic sites", "heteroplasmic"),
                ],
            );
            let a = val(d, "annotation");
            if let Some(n) = get(a, "clinvar_plp") {
                out.fact("ClinVar pathogenic or likely pathogenic", cell(n));
            }
            let rows = array(d, "heteroplasmic_sites")?
                .iter()
                .map(|x| {
                    vec![
                        format!(
                            "m.{}{}>{}",
                            text(val(x, "pos")),
                            text(val(x, "ref")),
                            text(val(x, "alt"))
                        ),
                        num(get(x, "af"))
                            .map_or_else(|| ".".into(), |n| format!("{}%", fixed(n * 100.0, 1))),
                    ]
                })
                .collect();
            out.table(&["Site", "Allele fraction"], rows);
            out.lines = mito::mito_lines(&sec)?;
            out.attention =
                num(get(a, "clinvar_plp")).is_some_and(|n| n > 0.0) || truthy(val(d, "confirmed"));
        }
        "lpa" => {
            let sec = summary.section("lpa");
            let d = sec.data;
            let repeats = if string(d, "status") == "estimated" {
                format!(
                    "{} (95% interval {} to {})",
                    text(val(d, "kiv2_copies")),
                    text(val(d, "ci95_low")),
                    text(val(d, "ci95_high"))
                )
            } else {
                cell(val(d, "status"))
            };
            out.fact("KIV-2 repeats", repeats);
            out.fact("Class", cell(val(d, "class")));
            if truthy(val(d, "null_allele")) {
                out.fact(
                    "KIV-2 4925G>A null allele",
                    smn1::marker_state(val(d, "null_allele"))?,
                );
            }
            out.lines = lpa::lines(&sec)?;
            out.attention = string(d, "flag") == "high";
            out.extra("widget", widgets::lpa(d));
        }
        "chip" => {
            let sec = summary.section("chip");
            let d = sec.data;
            let calls = array(d, "calls")?;
            out.fact("CHIP calls", calls.len().to_string());
            out.fact(
                "Typical detection limit at this depth (VAF)",
                cell(val(d, "median_limit")),
            );
            out.fields_table(
                &["Gene", "Change", "VAF", "Alt reads", "Depth", "Status"],
                calls,
                &[
                    "gene",
                    "protein_change",
                    "vaf",
                    "alt_reads",
                    "depth",
                    "status",
                ],
            );
            out.lines = chip::lines(&sec, &text(val(summary.root(), "sample")))?;
            out.attention = !calls.is_empty();
            out.extra("widget", widgets::chip(array(d, "genes")?, calls));
        }
        "mosaic" => {
            let sec = summary.section("mosaic_sex_loss");
            let d = sec.data;
            let call = string(d, "call");
            let label = match call {
                "not_detected" => "not detected above 10%",
                "gain" => "an extra copy, not a loss",
                "not_estimated" => "not estimated",
                _ => call,
            };
            out.fact(
                text_or(val(d, "chrom"), "chrY or chrX"),
                if label.is_empty() { "." } else { label },
            );
            out.lines.push(text_or(val(d, "reason"), ""));
            out.lines.extend(sex_loss::lines(&sec)?);
            out.attention = matches!(call, "detected" | "gain");
        }
        "mosaic-mca" => {
            let sec = summary.section("mosaic");
            let d = sec.data;
            let calls = array(d, "calls")?;
            out.fact("Status", cell(val(d, "status")));
            out.fact(
                "Calls reported",
                get(d, "calls_reported").map_or_else(|| calls.len().to_string(), cell),
            );
            out.fact("Smallest cell fraction it can see", cell(val(d, "floor")));
            out.fields_table(
                &[
                    "Chromosome",
                    "Start",
                    "End",
                    "Length (Mb)",
                    "Type",
                    "Cell fraction",
                ],
                calls,
                &[
                    "chrom",
                    "start",
                    "end",
                    "length_mb",
                    "type",
                    "cell_fraction",
                ],
            );
            out.lines = mocha::lines(&sec)?;
            out.attention = !calls.is_empty();
        }
        _ => unreachable!(),
    }
    Ok(out)
}

pub(super) fn findings(summary: &Summary, found: &mut Vec<Finding>) {
    for x in rows(val(data(summary, "annotsv"), "rows")) {
        add_row(found, "structural-variants", Level::Act, || {
            object_row(x)?;
            let genes = array(x, "genes")?;
            let gene = genes.first().unwrap_or(&Value::Null);
            let genes = genes
                .iter()
                .map(|g| match g {
                    Value::Str(t) => Ok(t.as_str()),
                    _ => Err(Malformed("gene is not a string".into())),
                })
                .collect::<Read<Vec<_>>>()?
                .join(", ");
            let kind = text_or(val(x, "type"), "structural change");
            let class = val(x, "acmg_class");
            let sv = val(x, "sv");
            let length = num(get(x, "length"))
                .or_else(|| string(x, "length").trim().parse().ok())
                .filter(|n| n.is_finite())
                .unwrap_or(0.0)
                .abs() as i128;
            let length = pycompat::fmt::int_grouped(length);
            let summary_genes = if genes.is_empty() { "no gene" } else { &genes };
            Ok(Some(Finding::new(if text(class) == "5" { Level::Act } else { Level::Review },
                "structural-variants", &["sv", &text_or_null(sv)], format!("Large {kind} over {summary_genes}"),
                format!("A {kind} of {length} bases that AnnotSV ranks class {}. Short reads can misplace large changes: confirm it.", text(class)))
                .gene(gene).detail(text_or(sv, ""))
                .evidence("Class", class).evidence("Callers", &s(&annotsv_support(summary, &D::new(x, "AnnotSV row"))?))
                .evidence("Genotype", val(x, "gt")).viewer(text_or(sv, "").replace(' ', ""))))
        });
    }
    for l in rows(val(data(summary, "expansions"), "loci")) {
        add_row(found, "expansions", Level::Act, || {
            object_row(l)?;
            let mut flags = repeat_flags(repeat_status(l)?);
            let first = flags.next();
            if first.is_none() && !truthy(val(l, "zone")) {
                return Ok(None);
            }
            let full = first
                .filter(|s| s.contains("full"))
                .or_else(|| flags.find(|s| s.contains("full")));
            let what = match full.or(first) {
                Some("pre_mutation") => "premutation range".to_owned(),
                Some("full_mutation") => "full-mutation range".to_owned(),
                Some(s) => s.replace('_', " "),
                None => format!("{} range", text_or(val(l, "zone_label"), "intermediate")),
            };
            let count = if truthy(val(l, "repeat_count")) {
                repeats::count(l)?
            } else {
                "not given".into()
            };
            let locus = val(l, "locus");
            let name = text(locus);
            let next = if full.is_some() {
                "See a clinical geneticist."
            } else {
                "It may grow when passed on; worth a word with a genetic counsellor."
            };
            let sizing = val(l, "allele_sizing");
            Ok(Some(
                Finding::new(
                    if full.is_some() {
                        Level::Act
                    } else {
                        Level::Review
                    },
                    "expansions",
                    &["str", &text_or_null(locus)],
                    format!("{name} repeat in the {what}"),
                    format!("{name} has {count} repeats: Stranger reads it as \"{what}\". {next}"),
                )
                .gene(locus)
                .detail(format!("{name} {count} repeats"))
                .evidence("Repeats", &s(&count))
                .evidence_fields(
                    l,
                    &[
                        ("95% interval", "ci"),
                        ("Pathologic from", "pathologic_min"),
                    ],
                )
                .evidence(
                    "Sizing",
                    if truthy(sizing) {
                        sizing
                    } else {
                        val(l, "sizing")
                    },
                ),
            ))
        });
    }
    let vn = data(summary, "vntyper");
    let status = string(vn, "status");
    if matches!(status, "positive" | "possible") {
        found.push(
            Finding::new(
                if status == "positive" {
                    Level::Act
                } else {
                    Level::Review
                },
                "vntyper",
                &["muc1", status],
                format!("MUC1 kidney gene: {status}"),
                format!(
                    "VNtyper reads the MUC1 repeat as {status}; it needs a clinical confirmation."
                ),
            )
            .gene(&s("MUC1"))
            .evidence("Confidence", val(vn, "confidence")),
        );
    }
    let lp = data(summary, "lpa");
    if string(lp, "flag") == "high" {
        found.push(
            Finding::new(
                Level::Review,
                "lpa",
                &["lpa"],
                "Lp(a) likely high",
                "The KIV-2 count suggests a high Lp(a); one blood measurement settles it.",
            )
            .gene(&s("LPA"))
            .evidence_fields(lp, &[("KIV-2 repeats", "kiv2_copies"), ("Class", "class")]),
        );
    }
    for c in rows(val(data(summary, "chip"), "calls")) {
        add_row(found, "chip", Level::Review, || {
            object_row(c)?;
            let gene = val(c, "gene");
            let change = val(c, "protein_change");
            Ok(Some(Finding::new(Level::Review, "chip", &["chip", &text_or_null(gene), &text_or_null(change)],
                format!("A blood clone in {} (CHIP)", text(gene)),
                format!("{} {} is carried by a share of the blood cells (VAF {}). It becomes common with age; most people with it never get a blood cancer.",
                    text(gene), text_or(change, ""), text(val(c, "vaf"))))
                .gene(gene).evidence("VAF", val(c, "vaf"))
                .evidence("Reads", &s(&format!("{} of {}", text(val(c, "alt_reads")), text(val(c, "depth")))))))
        });
    }
    for c in rows(val(data(summary, "mosaic"), "calls")) {
        add_row(found, "mosaic-mca", Level::Review, || {
            object_row(c)?;
            let kind = text_or(val(c, "type"), "change");
            let chrom = val(c, "chrom");
            let start = val(c, "start");
            let fraction = val(c, "cell_fraction");
            Ok(Some(
                Finding::new(
                    Level::Review,
                    "mosaic-mca",
                    &["mca", &text_or_null(chrom), &text_or_null(start)],
                    format!("A {kind} on {} in some blood cells", text(chrom)),
                    format!(
                        "MoChA sees a {kind} of {} Mb on {} in a share of the blood cells ({}).",
                        text(val(c, "length_mb")),
                        text(chrom),
                        text(fraction)
                    ),
                )
                .evidence(
                    "Where",
                    &s(&format!(
                        "{}:{}-{}",
                        text(chrom),
                        text(start),
                        text(val(c, "end"))
                    )),
                )
                .evidence("Cell fraction", fraction),
            ))
        });
    }
    for x in rows(val(data(summary, "exomiser"), "genes")).iter().take(3) {
        add_row(found, "exomiser", Level::Review, || {
            object_row(x)?;
            let gene = val(x, "gene");
            Ok(Some(
                Finding::new(
                    Level::Review,
                    "exomiser",
                    &["exomiser", &text_or_null(gene)],
                    format!(
                        "Symptom ranking: {} ranks {}",
                        text(gene),
                        text(val(x, "rank"))
                    ),
                    "A ranking to bring to a clinical geneticist, not a diagnosis.",
                )
                .gene(gene)
                .evidence_fields(x, &[("Score", "exomiser_score"), ("Inheritance", "moi")]),
            ))
        });
    }
}
