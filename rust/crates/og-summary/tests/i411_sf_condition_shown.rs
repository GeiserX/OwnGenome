//! A secondary finding keeps the condition ClinVar names for its variant.
//!
//! The ACMG tier, the carrier list and the ClinVar section map a step 06
//! hit to a gene through GENEINFO, which names the gene, not the disease.
//! So a RET variant that ClinVar asserts only for Hirschsprung disease
//! filled the MEN2 slot of the card like a MEN2 variant, and nothing in
//! summary.json said which disease the variant was reported for. These
//! cases run the collectors on synthetic hits with the shipped `assets/`.

use indexmap::IndexMap;
use og_summary::collect::acmg::{Inputs, acmg_sf};
use og_summary::collect::carrier::Carrier;
use og_summary::collect::clinvar::{ClinVar, Floor, HitRecords, HitsVcf};
use og_summary::collect::two_hits::DeletionInputs;
use pycompat::json::Value;

const ASSETS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets");

/// ClinVar's condition of a RET variant asserted for Hirschsprung disease
/// only, as step 06 copies it.
const HIRSCHSPRUNG: &str = ";CLNDN=Hirschsprung_disease,_susceptibility_to,_1\
                            ;CLNDISDB=MONDO:MONDO:0018309,MedGen:C1851851,OMIM:142623";
/// A MEN2 variant: one of its conditions is the ACMG disease.
const MEN2: &str = ";CLNDN=Multiple_endocrine_neoplasia,_type_2|not_provided\
                    ;CLNDISDB=MONDO:MONDO:0012000,MedGen:C0025269,OMIM:171400|MedGen:CN517202";

/// One step 06 hit: chr10 (RET's chromosome) unless `chrom` says, `gene`
/// in GENEINFO, Pathogenic with 2 stars, a confident heterozygous call,
/// and `info` added to INFO.
fn hit(chrom: &str, pos: i128, gene: &str, info: &str) -> Vec<String> {
    let info = format!(
        "GENEINFO={gene}:1;CLNSIG=Pathogenic;\
         CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts{info}"
    );
    [
        chrom,
        &pos.to_string(),
        ".",
        "C",
        "T",
        ".",
        "PASS",
        &info,
        "GT:GQ:DP",
        "0/1:60:35",
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect()
}

fn no_deletions() -> DeletionInputs {
    DeletionInputs {
        annotsv: Ok(None),
        step35: Ok(None),
    }
}

/// The ACMG tier over `hits`, with no clinical VCF.
fn tier(assets: &str, hits: Vec<Vec<String>>) -> Value {
    let h: Result<HitRecords, String> = Ok((Floor::default(), hits));
    let d = no_deletions();
    acmg_sf(&Inputs {
        assets,
        hits: Some(&h),
        clinical: None,
        deletions: &d,
    })
    .unwrap()
}

type Row = IndexMap<String, Value>;

fn rows<'a>(v: &'a Value, k: &str) -> Vec<&'a Row> {
    let Value::Object(m) = v else { panic!("{v:?}") };
    let Some(Value::Array(a)) = m.get(k) else {
        panic!("{k}: {v:?}")
    };
    a.iter()
        .map(|x| match x {
            Value::Object(o) => o,
            _ => panic!("{x:?}"),
        })
        .collect()
}

/// The reported ClinVar row of the variant at `pos`.
fn sf_row(v: &Value, pos: i128) -> &Row {
    let at = format!(":{pos} ");
    let found: Vec<&Row> = rows(v, "clinvar_hits")
        .into_iter()
        .filter(|m| matches!(m.get("variant"), Some(Value::Str(s)) if s.contains(&at)))
        .collect();
    assert_eq!(found.len(), 1, "one reported row at {pos}: {v:?}");
    found[0]
}

fn s<'a>(m: &'a Row, k: &str) -> &'a str {
    match m.get(k) {
        Some(Value::Str(x)) => x,
        x => panic!("{k}: {x:?} in {m:?}"),
    }
}

fn mark(m: &Row) -> Option<&Value> {
    m.get("condition_not_acmg")
}

#[test]
fn a_ret_hirschsprung_variant_is_reported_with_the_mark() {
    let v = tier(ASSETS, vec![hit("chr10", 43_100_000, "RET", HIRSCHSPRUNG)]);
    let m = sf_row(&v, 43_100_000);
    assert_eq!(s(m, "gene"), "RET");
    assert_eq!(
        s(m, "clinvar_condition"),
        "Hirschsprung disease, susceptibility to, 1"
    );
    assert_eq!(
        m.get("clinvar_condition_ids"),
        Some(&Value::Array(vec![
            Value::Str("MONDO:MONDO:0018309".into()),
            Value::Str("MedGen:C1851851".into()),
            Value::Str("OMIM:142623".into()),
        ]))
    );
    assert!(
        s(m, "acmg_phenotype").contains("multiple endocrine neoplasia type 2"),
        "{m:?}"
    );
    assert_eq!(mark(m), Some(&Value::Bool(true)), "{m:?}");
    // Shown and marked, never dropped: nothing of RET is set apart.
    assert!(
        rows(&v, "not_reported")
            .iter()
            .all(|r| s(r, "gene") != "RET"),
        "{v:?}"
    );
}

#[test]
fn a_ret_men2_variant_has_no_mark() {
    let v = tier(ASSETS, vec![hit("chr10", 43_114_500, "RET", MEN2)]);
    let m = sf_row(&v, 43_114_500);
    // "not provided" is left out of the name.
    assert_eq!(
        s(m, "clinvar_condition"),
        "Multiple endocrine neoplasia, type 2"
    );
    assert_eq!(mark(m), Some(&Value::Bool(false)), "{m:?}");
}

#[test]
fn the_two_ret_variants_side_by_side() {
    let v = tier(
        ASSETS,
        vec![
            hit("chr10", 43_100_000, "RET", HIRSCHSPRUNG),
            hit("chr10", 43_114_500, "RET", MEN2),
        ],
    );
    assert_eq!(mark(sf_row(&v, 43_100_000)), Some(&Value::Bool(true)));
    assert_eq!(mark(sf_row(&v, 43_114_500)), Some(&Value::Bool(false)));
}

#[test]
fn a_name_without_ids_matches_by_its_words() {
    // A hand-built ClinVar file can carry CLNDN alone.
    let v = tier(
        ASSETS,
        vec![hit("chr10", 43_114_500, "RET", ";CLNDN=MEN2A")],
    );
    assert_eq!(mark(sf_row(&v, 43_114_500)), Some(&Value::Bool(false)));
    let v = tier(
        ASSETS,
        vec![hit(
            "chr10",
            43_100_000,
            "RET",
            ";CLNDN=Hirschsprung_disease_1",
        )],
    );
    assert_eq!(mark(sf_row(&v, 43_100_000)), Some(&Value::Bool(true)));
}

#[test]
fn placeholders_and_panel_names_never_mark() {
    for info in [
        ";CLNDN=not_provided;CLNDISDB=MedGen:CN517202",
        ";CLNDN=Hereditary_cancer-predisposing_syndrome;CLNDISDB=MedGen:C0027672",
        ";CLNDN=not_specified|Inborn_genetic_diseases",
    ] {
        let v = tier(ASSETS, vec![hit("chr10", 43_114_500, "RET", info)]);
        let m = sf_row(&v, 43_114_500);
        assert_eq!(mark(m), Some(&Value::Bool(false)), "{info}: {m:?}");
    }
    // Only placeholders: the condition is there, with no name.
    let v = tier(
        ASSETS,
        vec![hit("chr10", 43_114_500, "RET", ";CLNDN=not_provided")],
    );
    assert_eq!(
        sf_row(&v, 43_114_500).get("clinvar_condition"),
        Some(&Value::Null)
    );
}

#[test]
fn a_hit_without_the_condition_gets_no_new_field() {
    // A hits VCF from before step 06 kept CLNDN: the row is as it was.
    let v = tier(ASSETS, vec![hit("chr10", 43_114_500, "RET", "")]);
    let m = sf_row(&v, 43_114_500);
    for k in [
        "clinvar_condition",
        "clinvar_condition_ids",
        "acmg_phenotype",
        "condition_not_acmg",
    ] {
        assert!(!m.contains_key(k), "{k} in {m:?}");
    }
}

#[test]
fn a_list_without_phenotype_columns_keeps_the_condition_unjudged() {
    // An acmg_sf.tsv of its own, without the two columns: the condition is
    // shown, and no mark claims a match or a mismatch.
    let shipped = std::fs::read_to_string(format!("{ASSETS}/acmg_sf.tsv")).unwrap();
    let mut cut = String::new();
    let mut drop: Option<usize> = None;
    for l in shipped.lines() {
        if l.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = l.split('\t').collect();
        let n = *drop.get_or_insert_with(|| {
            cells
                .iter()
                .position(|c| *c == "phenotype")
                .unwrap_or(cells.len())
        });
        cut.push_str(&cells[..n.min(cells.len())].join("\t"));
        cut.push('\n');
    }
    let dir = std::env::temp_dir().join(format!("og-sf-condition-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("acmg_sf.tsv"), cut).unwrap();
    let v = tier(
        &dir.display().to_string(),
        vec![hit("chr10", 43_100_000, "RET", HIRSCHSPRUNG)],
    );
    std::fs::remove_dir_all(&dir).unwrap();
    let m = sf_row(&v, 43_100_000);
    assert_eq!(
        s(m, "clinvar_condition"),
        "Hirschsprung disease, susceptibility to, 1"
    );
    assert!(!m.contains_key("acmg_phenotype"), "{m:?}");
    assert!(!m.contains_key("condition_not_acmg"), "{m:?}");
}

#[test]
fn every_acmg_gene_has_a_phenotype_and_keys() {
    let shipped = std::fs::read_to_string(format!("{ASSETS}/acmg_sf.tsv")).unwrap();
    let mut lines = shipped.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let at = |k: &str| head.iter().position(|c| *c == k).unwrap();
    let (p, k) = (at("phenotype"), at("phenotype_keys"));
    let mut n = 0;
    for l in lines {
        let c: Vec<&str> = l.split('\t').collect();
        n += 1;
        assert!(!matches!(c[p], "" | "."), "{}: no phenotype", c[0]);
        assert!(
            c[k].split(';').any(|x| x.trim().starts_with("OMIM:")),
            "{}: no MIM disorder number in {}",
            c[0],
            c[k]
        );
    }
    assert_eq!(n, 84);
}

#[test]
fn a_carrier_row_keeps_clinvar_condition() {
    let h: Result<HitRecords, String> = Ok((
        Floor::default(),
        vec![
            hit(
                "chr7",
                117_559_590,
                "CFTR",
                ";CLNDN=Cystic_fibrosis|not_provided;CLNDISDB=MONDO:MONDO:0009061,OMIM:219700|MedGen:CN517202",
            ),
            hit("chr7", 117_559_600, "CFTR", ""),
        ],
    ));
    let c = Carrier::read(ASSETS, &h, &no_deletions())
        .unwrap()
        .to_json();
    let carriers = rows(&c, "carriers");
    let at = |pos: &str| {
        *carriers
            .iter()
            .find(|m| s(m, "variant").contains(pos))
            .unwrap_or_else(|| panic!("{pos}: {c:?}"))
    };
    let with = at(":117559590 ");
    assert_eq!(s(with, "clinvar_condition"), "Cystic fibrosis");
    assert_eq!(
        with.get("clinvar_condition_ids"),
        Some(&Value::Array(vec![
            Value::Str("MONDO:MONDO:0009061".into()),
            Value::Str("OMIM:219700".into()),
        ]))
    );
    // The list's own condition for the gene stays.
    assert!(with.contains_key("condition"));
    let without = at(":117559600 ");
    assert!(!without.contains_key("clinvar_condition"), "{without:?}");
}

#[test]
fn the_clinvar_section_hit_keeps_clinvar_condition() {
    let hits = HitsVcf {
        path: "hits.vcf".into(),
        rel: "clinvar/S_clinvar_hits.vcf".into(),
        read: Ok((
            Floor::default(),
            vec![
                hit("chr10", 43_100_000, "RET", HIRSCHSPRUNG),
                hit("chr10", 43_114_500, "RET", ""),
            ],
        )),
    };
    let none = |_: &str| -> Option<(String, String)> { None };
    let cv = ClinVar::read(&hits, "", &none, "S").unwrap().to_json();
    let all = rows(&cv, "hits");
    let at = |pos: i128| {
        *all.iter()
            .find(|m| m.get("pos") == Some(&Value::Int(pos)))
            .unwrap_or_else(|| panic!("{pos}: {cv:?}"))
    };
    assert_eq!(
        s(at(43_100_000), "clinvar_condition"),
        "Hirschsprung disease, susceptibility to, 1"
    );
    assert!(!at(43_114_500).contains_key("clinvar_condition"));
    // The ACMG judgement is the card's, not the ClinVar table's.
    assert!(!at(43_100_000).contains_key("condition_not_acmg"));
}
