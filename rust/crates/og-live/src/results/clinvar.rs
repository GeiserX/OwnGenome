use super::*;
use og_render::sections::{
    cdna, cdna_scan,
    clinvar::{self, D, HITS_SHOWN},
    deep_intronic,
};

fn variant(h: &Value) -> String {
    format!(
        "{}:{} {}>{}",
        text(val(h, "chrom")),
        text(val(h, "pos")),
        text(val(h, "ref")),
        text(val(h, "alt"))
    )
}

fn locus(h: &Value) -> String {
    format!("{}:{}", text(val(h, "chrom")), text(val(h, "pos")))
}

fn stars(h: &Value, empty: &str) -> Read<String> {
    let stars = match val(h, "stars") {
        v if !truthy(v) => 0,
        Value::Int(n) => *n,
        Value::Float(n) if n.is_finite() => *n as i128,
        Value::Bool(n) => i128::from(*n),
        Value::Str(s) => s
            .trim()
            .parse()
            .map_err(|_| Malformed("stars is not an integer".into()))?,
        _ => return Err(Malformed("stars is not an integer".into())),
    };
    Ok(if stars <= 0 {
        empty.into()
    } else {
        "★".repeat(stars.min(1000) as usize)
    })
}

fn call(h: &D<'_>) -> Read<String> {
    let mut parts = Vec::new();
    let call = clinvar::call_text(h)?;
    if !call.is_empty() {
        parts.push(call);
    }
    if h.has("cdna") {
        parts.push(format!("reads: {}", cdna::cdna_text(h)?));
    }
    let combined = clinvar::combined_cell(h)?;
    if !combined.is_empty() {
        parts.push(combined);
    }
    Ok(if parts.is_empty() {
        ".".into()
    } else {
        parts.join("; ")
    })
}

fn frequency(h: &D<'_>) -> Read<String> {
    let frequency = clinvar::af_text(h)?;
    Ok(if frequency.is_empty() {
        "frequency unknown".into()
    } else {
        frequency
    })
}

pub(super) fn card(summary: &Summary) -> Read<Adapter> {
    let sec = summary.section("clinvar");
    let d = sec.data;
    let mut out = Adapter::default();
    out.fact("P/LP matches", default(d, "count", "0"));
    let by = val(d, "by_stars");
    out.fact(
        "By review stars (4/3/2/1/0)",
        ["4", "3", "2", "1", "0"]
            .map(|n| default(by, n, "0"))
            .join(" / "),
    );
    let dbs = summary.get("databases").unwrap_or(&Value::Null);
    out.fact(
        "ClinVar file date",
        text_or(val(dbs, "clinvar_release"), "unknown"),
    );
    if let Some(n) = get(d, "low_confidence") {
        let f = val(d, "call_floor");
        out.fact(
            format!(
                "Low-confidence calls (GQ below {} or DP below {}), not counted",
                text(val(f, "GQ")),
                text(val(f, "DP"))
            ),
            cell(n),
        );
    }
    let hits = array(d, "hits")?;
    let mut table = Vec::new();
    let mut viewer = Vec::new();
    for h in hits.iter().take(HITS_SHOWN) {
        object_row(h)?;
        let hit = D::new(h, "ClinVar hit");
        table.push(vec![
            cell(val(h, "gene")),
            variant(h),
            cell(val(h, "genotype")),
            cell(val(h, "significance")),
            stars(h, "0")?,
            frequency(&hit)?,
            call(&hit)?,
        ]);
        viewer.push(locus(h));
    }
    out.table(
        &[
            "Gene",
            "Variant",
            "Genotype",
            "Significance",
            "Stars",
            "Frequency",
            "Call",
        ],
        table,
    );
    out.extra("viewer", strings(viewer));
    if let Some(line) = clinvar::freq_line(&D::of(&sec))? {
        out.lines.push(line);
    }
    out.lines.extend(cdna_scan::lines(&sec)?);
    for x in array(d, "deep_intronic")? {
        out.lines.push(format!(
            "Known deep intronic site (always reported): {} {} {}: {}",
            text(val(x, "gene")),
            text(val(x, "hgvs")),
            variant(x),
            deep_intronic::call(&D::new(x, "deep intronic site"))?
        ));
    }
    Ok(out)
}

pub(super) fn findings(summary: &Summary, excluded: &[&Value], found: &mut Vec<Finding>) {
    for h in rows(val(data(summary, "clinvar"), "hits"))
        .iter()
        .take(HITS_SHOWN)
    {
        add_row(found, "clinvar", Level::Review, || {
            object_row(h)?;
            let gene = val(h, "gene");
            if excluded.contains(&gene) {
                return Ok(None);
            }
            let d = D::new(h, "ClinVar hit");
            let significance = val(h, "significance");
            let genotype = val(h, "genotype");
            let copies = match str_val(genotype) {
                "hom" => "two copies".into(),
                "het" => "one copy".into(),
                _ => text_or(genotype, "genotype not given"),
            };
            let key = ["gene", "chrom", "pos", "ref", "alt"].map(|k| text_or_null(val(h, k)));
            let mut finding = Finding::new(Level::Review, "clinvar",
                &["clinvar", &key[0], &key[1], &key[2], &key[3], &key[4]],
                format!("{} change in {}", text_or(significance, "Pathogenic"), text(gene)),
                format!("ClinVar lists this change in {} as {}; this genome has {copies}. Whether it matters depends on how the gene is inherited: a genetic counsellor reads it with the rest of this report.",
                    text(gene), text_or(significance, "").to_lowercase()))
                .gene(gene).detail(format!("{} {} {}", text(gene), variant(h), text_or(genotype, "")))
                .evidence("Significance", significance).evidence("ClinVar review", &s(&stars(h, "0 stars")?))
                .evidence("Genotype", genotype).evidence("Frequency", &s(&frequency(&d)?))
                .evidence("Call", &s(&call(&d)?)).viewer(locus(h));
            finding.af = val(h, "af").clone();
            Ok(Some(finding))
        });
    }
}

pub(super) fn text_or_null(v: &Value) -> String {
    if matches!(v, Value::Null) {
        String::new()
    } else {
        text(v)
    }
}
