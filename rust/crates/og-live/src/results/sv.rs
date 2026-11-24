use super::clinvar::text_or_null;
use super::*;
use crate::widgets;
use og_render::sections::{
    annotsv, chip, clinvar::D, lpa, mito, mocha, sex_loss, smn1, str as repeats, vntyper,
};
use pycompat::fmt::fixed;

/// A class 4-5 row's benign-set AF for its finding: "no frequency data"
/// when step 05 says AnnotSV gave none, never blank, which would read as
/// rare. Null (no evidence line) for a row from before step 05 wrote its
/// per-record evidence, whatever its AF, so such a finding reads as before.
fn frequency(x: &Value) -> Read<Value> {
    let row = D::new(x, "AnnotSV row");
    if annotsv::no_frequency(&row)? {
        return Ok(s("no frequency data (unknown, not the same as rare)"));
    }
    if !annotsv::has_evidence(&row)? {
        return Ok(Value::Null);
    }
    Ok(match val(x, "population_af") {
        Value::Null => Value::Null,
        af => s(&format!(
            "{} (AnnotSV's aggregate over the benign records listed)",
            text(af)
        )),
    })
}

/// Step 05's `Benign_match` of a row in words, as the card's checks say it
/// ("only larger SVs around it"); null without the key.
fn benign_records(x: &Value) -> Read<Value> {
    Ok(match D::new(x, "AnnotSV row").opt("benign_match")? {
        Some(m) => s(annotsv::benign_match_text(m.string()?)),
        None => Value::Null,
    })
}

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
            if truthy(a) {
                // Step 05's per-record evidence and unknown frequencies, said
                // only when a row carries them, as the report says them.
                let d = D::new(a, "sections.annotsv.data");
                if let Some(n) = annotsv::same_size_common(&d)? {
                    out.fact(annotsv::SAME_SIZE_LABEL, n);
                }
                if annotsv::shows_evidence_note(&d)? {
                    out.lines.push(annotsv::BENIGN_EVIDENCE_NOTE.to_owned());
                }
            }
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
            // The call MEGAnE made at the site, so "call filtered out (LC)"
            // names the flagged call as the report does.
            out.fields_table(
                &[
                    "Gene",
                    "Element",
                    "Condition",
                    "Inheritance",
                    "Known site",
                    "Depth",
                    "Call",
                ],
                sites,
                &[
                    "gene",
                    "element",
                    "condition",
                    "inheritance",
                    "status",
                    "mean_depth",
                    "call",
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
            let mut penetrance = Vec::new();
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
                // Stranger reads a count of reduced penetrance as a plain full
                // mutation: the range says what it is.
                let shown = match repeats::penetrance_text(l)? {
                    Some(p) => format!("{shown}; {p}"),
                    None => shown,
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
                penetrance.push(val(l, "penetrance").clone());
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
            let mut widget = widgets::str_gauge(gauge);
            // person.js names a locus's range of reduced penetrance under its
            // gauge; a locus without one carries no key.
            if let Value::Object(w) = &mut widget
                && let Some(Value::Array(loci)) = w.get_mut("loci")
            {
                for (g, p) in loci.iter_mut().zip(penetrance) {
                    if let Value::Object(g) = g
                        && truthy(&p)
                    {
                        g.insert("penetrance".into(), p);
                    }
                }
            }
            out.extra("widget", widget);
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
            let message = string(d, "message");
            let repeats = if string(d, "status") == "estimated" {
                format!(
                    "{} ({})",
                    text(val(d, "kiv2_copies")),
                    lpa::interval_text(
                        message,
                        &text(val(d, "ci95_low")),
                        &text(val(d, "ci95_high"))
                    )
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
            let mut widget = widgets::lpa(d);
            // The whisker's label says what the interval covers (step 46's
            // message); person.js reads it, a widget without it keeps "95%
            // interval".
            if let (Value::Object(w), Some(label)) = (&mut widget, lpa::interval_label(message)) {
                w.insert("interval".into(), s(label));
            }
            out.extra("widget", widget);
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
                .evidence("Genotype", val(x, "gt"))
                .evidence("Benign records", &benign_records(x)?)
                .evidence_fields(x, &[("Same-size benign record at 1% or more", "clinical_common")])
                .evidence("Benign-set AF", &frequency(x)?)
                .viewer(text_or(sv, "").replace(' ', ""))))
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
            // A count Stranger reads as a full mutation can sit in a range of
            // reduced penetrance (HTT 36-39). Each allele is read on its own:
            // the title follows the most severe one, so a pair with an allele
            // above the range (or one not sized) keeps Stranger's reading,
            // and the sentence names the allele of reduced penetrance.
            let alleles = repeats::penetrance_alleles(l)?;
            let (title, said) = match repeats::penetrance_text(l)? {
                Some(p) if matches!(alleles, Some((_, false))) => (
                    format!("{name} repeat in a reduced-penetrance range"),
                    format!(
                        "{name} has {count} repeats: Stranger reads it as \"{what}\", but {p}. {next}"
                    ),
                ),
                Some(p) => (
                    format!("{name} repeat in the {what}"),
                    format!(
                        "{name} has {count} repeats: Stranger reads it as \"{what}\"; {p}. {next}"
                    ),
                ),
                None => (
                    format!("{name} repeat in the {what}"),
                    format!("{name} has {count} repeats: Stranger reads it as \"{what}\". {next}"),
                ),
            };
            Ok(Some(
                Finding::new(
                    if full.is_some() {
                        Level::Act
                    } else {
                        Level::Review
                    },
                    "expansions",
                    &["str", &text_or_null(locus)],
                    title,
                    said,
                )
                .gene(locus)
                .detail(format!("{name} {count} repeats"))
                .evidence("Repeats", &s(&count))
                .evidence_fields(
                    l,
                    &[
                        ("95% interval", "ci"),
                        ("Pathologic from", "pathologic_min"),
                        ("Reduced penetrance", "penetrance"),
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

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(sections: &str) -> Summary {
        Summary::parse(&format!(r#"{{"sections": {{{sections}}}}}"#)).unwrap()
    }

    fn findings_of(summary: &Summary) -> Vec<Finding> {
        let mut found = Vec::new();
        findings(summary, &mut found);
        found
    }

    /// The live card and finding say "no frequency data", tell the
    /// aggregate AF from the same-size record, and count those records; a
    /// row from before the evidence columns reads as before.
    #[test]
    fn live_sv_shows_unknown_frequency_and_the_per_record_evidence() {
        let summary = parse(
            r#""annotsv": {"state": "ok", "data": {"count": 2, "consensus": null, "rows": [
                {"type": "DEL", "sv": "chr1:1000-5000", "gt": "0/1", "length": "-4000", "acmg_class": "4",
                 "genes": ["G1"], "caller": "manta", "support": 1, "support_callers": ["manta"],
                 "exon_edge": "no", "population_af": null, "population_source": "no frequency data",
                 "common": false, "benign_match": "none", "clinical_common": "no"},
                {"type": "DEL", "sv": "chr3:1000-9000", "gt": "0/1", "length": "-8000", "acmg_class": "5",
                 "genes": ["G3"], "caller": "manta", "support": 1, "support_callers": ["manta"],
                 "exon_edge": "no", "population_af": 0.05, "population_source": "gnomAD-SV_v4.1",
                 "common": true, "benign_match": "same", "clinical_common": "yes"}]}}"#,
        );
        let shown = card("structural-variants", &summary).unwrap();
        let table = shown.table.unwrap();
        let rows = val(&table, "rows");
        let checks = |i: usize| match rows {
            Value::Array(r) => text(&list(&r[i]).last().cloned().unwrap_or(Value::Null)),
            _ => String::new(),
        };
        assert_eq!(
            checks(0),
            "benign-set AF: no frequency data (unknown, not the same as rare); benign records: \
             none listed, same-size record at 1% or more: no"
        );
        assert_eq!(
            checks(1),
            "benign-set AF 0.05 (gnomAD-SV_v4.1; the highest AF of the benign records listed, \
             larger SVs around it included), common; benign records: one of about the same size, \
             same-size record at 1% or more: yes"
        );
        assert!(
            shown.facts.contains(&strings([
                "Calls with a confirmed same-size benign record at 1% or more",
                "1"
            ])),
            "{:?}",
            shown.facts
        );
        assert_eq!(shown.lines, [annotsv::BENIGN_EVIDENCE_NOTE]);

        let found = findings_of(&summary);
        let evidence = |i: usize| {
            found[i]
                .evidence
                .iter()
                .map(|e| list(e).iter().map(text).collect::<Vec<_>>().join(": "))
                .collect::<Vec<_>>()
        };
        assert!(
            evidence(0).contains(
                &"Benign-set AF: no frequency data (unknown, not the same as rare)".into()
            ),
            "{:?}",
            evidence(0)
        );
        assert!(
            evidence(1).contains(&"Same-size benign record at 1% or more: yes".into()),
            "{:?}",
            evidence(1)
        );

        // Control: a row from before step 05 wrote the evidence.
        let old = parse(
            r#""annotsv": {"state": "ok", "data": {"count": 1, "consensus": null, "rows": [
                {"type": "DEL", "sv": "chr1:1000-5000", "acmg_class": "4", "genes": ["G1"],
                 "caller": "manta", "support": 1, "exon_edge": "no", "population_af": null,
                 "population_source": ""}]}}"#,
        );
        let shown = card("structural-variants", &old).unwrap();
        assert!(shown.lines.is_empty(), "{:?}", shown.lines);
        assert_eq!(shown.facts.len(), 5, "{:?}", shown.facts);
        let found = findings_of(&old);
        assert!(
            !found[0]
                .evidence
                .iter()
                .any(|e| text(&list(e)[0]).starts_with("Benign-set")),
            "{:?}",
            found[0].evidence
        );
    }

    /// The evidence of a finding as "label: value" lines.
    fn evidence_of(f: &Finding) -> Vec<String> {
        f.evidence
            .iter()
            .map(|e| list(e).iter().map(text).collect::<Vec<_>>().join(": "))
            .collect()
    }

    /// A legacy row (no per-record evidence, no "no frequency data" marker)
    /// with a numeric AF gains no new evidence line: its finding reads as
    /// before step 05 wrote the columns.
    #[test]
    fn live_sv_legacy_numeric_af_adds_no_evidence_line() {
        let old = parse(
            r#""annotsv": {"state": "ok", "data": {"count": 1, "consensus": null, "rows": [
                {"type": "DEL", "sv": "chr2:2000-2500", "acmg_class": "4", "genes": ["G1"],
                 "caller": "manta", "support": 1, "exon_edge": "no", "population_af": 0.05,
                 "population_source": "gnomAD-SV_v4.1", "common": true}]}}"#,
        );
        let found = findings_of(&old);
        let evidence = evidence_of(&found[0]);
        assert!(
            !evidence
                .iter()
                .any(|e| e.starts_with("Benign") || e.starts_with("Same-size")),
            "{evidence:?}"
        );
    }

    /// A finding names step 05's Benign_match in words, as the card does:
    /// "only larger SVs around it" is what tells a contained match with no
    /// same-size common record from a same-size one.
    #[test]
    fn live_sv_finding_names_the_benign_match() {
        let summary = parse(
            r#""annotsv": {"state": "ok", "data": {"count": 1, "consensus": null, "rows": [
                {"type": "DEL", "sv": "chr2:2000-2500", "gt": "0/1", "length": "-500", "acmg_class": "4",
                 "genes": ["G1"], "caller": "manta", "support": 1, "support_callers": ["manta"],
                 "exon_edge": "no", "population_af": 0.05, "population_source": "gnomAD-SV_v4.1",
                 "common": true, "benign_match": "contained", "clinical_common": "no"}]}}"#,
        );
        let evidence = evidence_of(&findings_of(&summary)[0]);
        assert!(
            evidence.contains(&"Benign records: only larger SVs around it".into()),
            "{evidence:?}"
        );
        assert!(
            evidence.contains(&"Same-size benign record at 1% or more: no".into()),
            "{evidence:?}"
        );
    }

    /// The card counts confirmed same-size matches and says how many the
    /// evidence could not assess, so every row unknown never reads as a bare
    /// 0; its caveat names no count the card does not show.
    #[test]
    fn live_sv_same_size_count_says_confirmed_and_not_assessed() {
        let row = |sv: &str| {
            format!(
                r#"{{"type": "DEL", "sv": "{sv}", "acmg_class": "4", "genes": ["G1"],
                    "caller": "manta", "support": 1, "exon_edge": "no", "population_af": 0.05,
                    "population_source": "gnomAD-SV_v4.1", "common": true,
                    "benign_match": "same", "clinical_common": "unknown"}}"#
            )
        };
        let summary = parse(&format!(
            r#""annotsv": {{"state": "ok", "data": {{"count": 2, "consensus": null, "rows": [{}, {}]}}}}"#,
            row("chr1:1000-5000"),
            row("chr3:1000-9000")
        ));
        let shown = card("structural-variants", &summary).unwrap();
        assert!(
            shown.facts.contains(&strings([
                "Calls with a confirmed same-size benign record at 1% or more",
                "0, not assessed: 2"
            ])),
            "{:?}",
            shown.facts
        );
        assert_eq!(shown.lines.len(), 1, "{:?}", shown.lines);
        assert!(
            !shown.lines[0].contains("common count"),
            "{:?}",
            shown.lines
        );
    }

    /// The MEI table names the call MEGAnE flagged, beside its status, and
    /// shows a site with too few reads as not callable.
    #[test]
    fn live_mei_names_the_filtered_call_and_the_site_not_assessed() {
        let summary = parse(
            r#""mei": {"state": "ok", "data": {"calls": 1, "by_element": {"SVA": 1}, "sites": [
                {"gene": "FKTN", "element": "SVA", "condition": "FCMD", "inheritance": "AR",
                 "status": "call filtered out (LC)", "mean_depth": "30.0",
                 "call": "chr9:105639660 SVA_F 1 copy LC"},
                {"gene": "TAF1", "element": "SVA", "condition": "XDP", "inheritance": "XL",
                 "status": "not callable", "mean_depth": "4.0", "call": "."}]}}"#,
        );
        let shown = card("mei", &summary).unwrap();
        let table = shown.table.unwrap();
        assert_eq!(text(list(val(&table, "columns")).last().unwrap()), "Call");
        let Value::Array(rows) = val(&table, "rows") else {
            panic!("{table:?}")
        };
        let row = |i: usize| list(&rows[i]).iter().map(text).collect::<Vec<_>>();
        assert_eq!(
            row(0)[4..],
            [
                "call filtered out (LC)",
                "30.0",
                "chr9:105639660 SVA_F 1 copy LC"
            ]
        );
        assert_eq!(row(1)[4], "not callable");
        assert!(shown.attention);
    }

    /// The range of reduced penetrance `og summary collect` puts on an HTT
    /// locus with a sized count of 36 to 39.
    const PEN: &str = r#", "penetrance": "36-39", "penetrance_label": "reduced penetrance (may or may not cause Huntington disease; full penetrance from 40)", "penetrance_citation": "GeneReviews""#;

    /// The invented counts of a mixed HTT pair: one allele in the range of
    /// reduced penetrance, the other above it, a full mutation.
    const MIXED: &str = "38/42";

    /// Invented controls: one count of reduced penetrance; one full
    /// mutation (of 42, and of 40, the first count of full penetrance).
    const REDUCED: &str = "17/38";
    const FULL: &str = "17/42";
    const FULL_40: &str = "17/40";

    /// A summary whose only expansions locus is HTT with `counts` (sized,
    /// flagged by Stranger as a full mutation) and `extra` keys.
    fn htt(counts: &str, extra: &str) -> Summary {
        let ci = counts
            .split('/')
            .map(|c| format!("{c}-{c}"))
            .collect::<Vec<_>>()
            .join("/");
        parse(&format!(
            r#""expansions": {{"state": "ok", "data": {{"records": 1, "stranger": true,
                "loci": [{{"locus": "HTT", "repeat_count": "{counts}", "status": "full_mutation",
                "normal_max": "26", "pathologic_min": "36", "ci": "{ci}", "sizing": "sized",
                "allele_sizing": "sized/sized", "zone": ""{extra}}}], "unsized": [],
                "intermediate": []}}}}"#
        ))
    }

    /// An HTT count of reduced penetrance says so in the table, the gauge
    /// and the finding; a count of 40 reads as before.
    #[test]
    fn live_htt_of_reduced_penetrance_says_so() {
        let summary = htt(REDUCED, PEN);
        let shown = card("expansions", &summary).unwrap();
        let table = shown.table.unwrap();
        let Value::Array(rows) = val(&table, "rows") else {
            panic!("{table:?}")
        };
        assert_eq!(
            text(&list(&rows[0])[2]),
            "\"full_mutation\"; 38 repeats in 36-39: reduced penetrance (may or may not cause \
             Huntington disease; full penetrance from 40)"
        );
        let widget = shown.extra.get("widget").unwrap();
        assert_eq!(
            text(val(&list(val(widget, "loci"))[0], "penetrance")),
            "36-39"
        );
        let found = findings_of(&summary);
        assert_eq!(found[0].title, "HTT repeat in a reduced-penetrance range");
        assert!(
            found[0]
                .text
                .contains("but 38 repeats in 36-39: reduced penetrance"),
            "{}",
            found[0].text
        );

        let summary = htt(FULL_40, "");
        let shown = card("expansions", &summary).unwrap();
        let widget = shown.extra.get("widget").unwrap();
        assert!(get(&list(val(widget, "loci"))[0], "penetrance").is_none());
        let found = findings_of(&summary);
        assert_eq!(found[0].title, "HTT repeat in the full-mutation range");
    }

    /// Each sized HTT allele is read on its own: in a mixed pair the
    /// finding keeps the full-mutation reading of the allele above the
    /// range and names the other as reduced penetrance. Controls: a pair
    /// with only a count of reduced penetrance, and one with only a full
    /// mutation.
    #[test]
    fn live_htt_mixed_pair_keeps_the_full_mutation_reading() {
        let found = findings_of(&htt(MIXED, PEN));
        assert_eq!(found[0].title, "HTT repeat in the full-mutation range");
        assert!(
            found[0].text.contains(
                "Stranger reads it as \"full-mutation range\"; 38 repeats in 36-39: reduced \
                 penetrance (may or may not cause Huntington disease; full penetrance from 40). \
                 See a clinical geneticist."
            ),
            "{}",
            found[0].text
        );
        assert!(!found[0].text.contains(", but "), "{}", found[0].text);

        let found = findings_of(&htt(REDUCED, PEN));
        assert_eq!(found[0].title, "HTT repeat in a reduced-penetrance range");
        let found = findings_of(&htt(FULL, ""));
        assert_eq!(found[0].title, "HTT repeat in the full-mutation range");
        assert!(!found[0].text.contains("reduced"), "{}", found[0].text);
    }

    /// The KIV-2 fact and the whisker say the interval covers read sampling
    /// only; a message from before keeps "95% interval".
    #[test]
    fn live_lpa_labels_the_interval_read_sampling() {
        let at = |message: &str| {
            parse(&format!(
                r#""lpa": {{"state": "ok", "data": {{"status": "estimated", "kiv2_copies": "29.4",
                    "ci95_low": "23.8", "ci95_high": "35.8", "class": "typical", "message": "{message}",
                    "advice": ".", "flag": "", "flag_text": "."}}}}"#
            ))
        };
        let shown = card(
            "lpa",
            &at(
                "about 29 KIV-2 repeats in total over both copies (95% interval for read \
                 sampling, 23.8 to 35.8; aligner or library bias not included): typical",
            ),
        )
        .unwrap();
        assert_eq!(
            shown.facts[0],
            strings([
                "KIV-2 repeats",
                "29.4 (95% interval for read sampling only, 23.8 to 35.8)"
            ])
        );
        let widget = shown.extra.get("widget").unwrap();
        assert_eq!(
            text(val(widget, "interval")),
            "95% interval for read sampling only"
        );

        let shown = card("lpa", &at("ok")).unwrap();
        assert_eq!(
            shown.facts[0],
            strings(["KIV-2 repeats", "29.4 (95% interval 23.8 to 35.8)"])
        );
        assert!(get(shown.extra.get("widget").unwrap(), "interval").is_none());
    }
}
