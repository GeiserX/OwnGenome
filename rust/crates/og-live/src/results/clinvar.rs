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
    conflicting(&D::of(&sec), &mut out)?;
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

/// The conflicting records some labs call P/LP: a fact and one line each,
/// never a finding and never in the counts above.
fn conflicting(d: &D<'_>, out: &mut Adapter) -> Read<()> {
    let Some(c) = d.get("conflicting")? else {
        return Ok(());
    };
    out.fact(clinvar::CONFLICTING_TITLE, clinvar::conflicting_count(&c)?);
    let hits = c.at("hits")?.items()?;
    for h in hits.iter().take(HITS_SHOWN) {
        out.lines
            .push(format!("To review: {}", clinvar::conflicting_line(h)?));
    }
    if hits.len() > HITS_SHOWN {
        out.lines.push(format!(
            "To review: {} more in clinvar/<sample>_clinvar_conflicting.tsv",
            hits.len() - HITS_SHOWN
        ));
    }
    if !hits.is_empty() || !c.flag("checked")? {
        out.lines.push(clinvar::CONFLICTING_NOTE.to_owned());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn summary(conflicting: &str) -> Summary {
        Summary::parse(&format!(
            r#"{{"databases": {{"clinvar_release": "2026-09-28"}}, "sections": {{"clinvar": {{"state": "ok", "data": {{
                "count": 1, "by_stars": {{"4": 0, "3": 0, "2": 1, "1": 0, "0": 0}},
                "hits": [{{"chrom": "chr2", "pos": 47403191, "ref": "C", "alt": "T", "genotype": "het",
                  "gene": "DEMOA", "significance": "Pathogenic", "stars": 2}}]{conflicting}}}}}}}}}"#
        ))
        .unwrap()
    }

    const LIST: &str = r#", "conflicting": {"checked": true, "note": null, "count": 1, "hits": [
        {"chrom": "chr5", "pos": 31400123, "ref": "G", "alt": "A", "genotype": "het", "gene": "DEMOD",
         "stars": 1, "submissions": "Pathogenic(2), Uncertain_significance(1)",
         "clinvar_condition": "Demo condition D", "low_confidence": false}]}"#;

    #[test]
    fn the_review_list_is_on_the_card_and_never_a_finding() {
        let with = card(&summary(LIST)).unwrap();
        let plain = card(&summary("")).unwrap();
        assert_eq!(with.facts[..plain.facts.len()], plain.facts[..]);
        assert_eq!(
            with.facts[plain.facts.len()],
            strings([clinvar::CONFLICTING_TITLE, "1"])
        );
        assert!(with.lines.contains(
            &"To review: DEMOD chr5:31400123 G>A het: Pathogenic(2), Uncertain_significance(1); 1 star; Demo condition D".to_owned()
        ));
        assert!(with.lines.contains(&clinvar::CONFLICTING_NOTE.to_owned()));
        assert_eq!(with.table, plain.table);
        let (mut a, mut b) = (Vec::new(), Vec::new());
        findings(&summary(LIST), &[], &mut a);
        findings(&summary(""), &[], &mut b);
        assert_eq!(a.len(), 1);
        assert_eq!(
            a.iter().map(|f| f.id.clone()).collect::<Vec<_>>(),
            b.iter().map(|f| f.id.clone()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn an_unchecked_list_says_why() {
        let shown = card(&summary(
            r#", "conflicting": {"checked": false, "note": "not checked: run setup.sh --refresh clinvar", "count": null, "hits": []}"#,
        ))
        .unwrap();
        assert!(shown.facts.contains(&strings([
            clinvar::CONFLICTING_TITLE,
            "not checked: run setup.sh --refresh clinvar"
        ])));
    }
}
