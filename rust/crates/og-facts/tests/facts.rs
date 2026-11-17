//! `og facts` on an invented sample folder: one test per rule OwnHealth's
//! importer and genomics.md hold the file to. Nothing here is a real person:
//! the genome is a made-up repeating sequence, the calls are invented, and
//! the PharmCAT report is PharmCAT's own example
//! (tests/fixtures/pharmcat/pharmcat-docs-example.json).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use og_facts::build::{Built, Genome, build};
use og_facts::panels::Panels;
use og_facts::rules::{fact_id, sha256_hex, term_ok};
use pycompat::json::{self, Value};

const SAMPLE: &str = "SYN-01";

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn panels() -> Panels {
    Panels::load(&repo().join("assets")).expect("assets load")
}

/// A made-up genome: `ACGT` repeated, with a few bases set.
struct Synth(HashMap<(String, u64), u8>);

fn base(p0: u64) -> u8 {
    b"ACGT"[(p0 % 4) as usize]
}

impl Genome for Synth {
    fn get(&mut self, chrom: &str, beg: u64, end: u64) -> Result<Vec<u8>, String> {
        Ok((beg..end)
            .map(|p| *self.0.get(&(chrom.to_owned(), p)).unwrap_or(&base(p)))
            .collect())
    }
}

/// `pos` (1-based): its reference base and another one.
fn snv(pos: u64) -> (String, String) {
    let r = base(pos - 1);
    let a = base(pos);
    (char::from(r).to_string(), char::from(a).to_string())
}

fn hit(
    chrom: &str,
    pos: u64,
    gene: &str,
    zyg: &str,
    sig: &str,
    stars: u8,
    ref_: Option<&str>,
) -> String {
    let (r, a) = snv(pos);
    format!(
        r#"{{"stars":{stars},"chrom":"{chrom}","pos":{pos},"ref":"{}","alt":"{a}","genotype":"{zyg}","gene":"{gene}",
        "significance":"{sig}","dp":35,"gq":60,"low_confidence":false,"rsids":["rs999{pos}"]}}"#,
        ref_.unwrap_or(&r)
    )
}

/// step 36's CYP2D6 block: PharmCAT made no call, pypgx a probable one.
const CYP2D6_PROBABLE: &str = r#"{"calls":{"PharmCAT":"not called","pypgx":"*4/*5"},"agree":null,"consensus":null,
    "probable":{"diplotype":"*4/*5","phenotype":"Poor Metabolizer","label":"probable (one caller: pypgx)"}}"#;

/// Writes the invented sample folder and returns it.
fn folder(name: &str, sections: &str) -> PathBuf {
    folder_with(name, sections, CYP2D6_PROBABLE)
}

fn folder_with(name: &str, sections: &str, cyp2d6: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("og-facts-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    for d in ["vcf", "pharmcat", "clinvar"] {
        std::fs::create_dir_all(dir.join(d)).unwrap();
    }
    std::fs::copy(
        repo().join("tests/fixtures/pharmcat/pharmcat-docs-example.json"),
        dir.join(format!("pharmcat/{SAMPLE}.report.json")),
    )
    .unwrap();
    std::fs::write(dir.join(format!("vcf/{SAMPLE}.vcf.gz")), b"calls").unwrap();
    std::fs::write(dir.join(format!("vcf/{SAMPLE}.vcf.gz.tbi")), b"calls index").unwrap();
    std::fs::write(dir.join(format!("vcf/{SAMPLE}.g.vcf.gz")), b"gvcf").unwrap();
    let (r, a) = snv(117_300_000);
    std::fs::write(
        dir.join(format!("clinvar/{SAMPLE}_clinvar_hits.tsv")),
        format!("chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat\nchr7\t117300000\t{r}\t{a}\t0/1\t7105\tCFTR:1080\tPathogenic\t.\n"),
    )
    .unwrap();
    let summary = format!(
        r#"{{"schema_version":2,"sample":"{SAMPLE}","generated_utc":"2026-10-01 09:15 UTC",
        "run":{{"declared_sex":"female"}},
        "manifest":[["run","pipeline_version","1.4.2"],["run","git_commit","abc1234"],
                    ["data","clinvar_file_date","2026-09-28 (clinvar/clinvar.vcf.gz)"],["image","x","y"]],
        "cyp2d6":{cyp2d6},
        "sections":{{{sections}}}}}"#
    );
    std::fs::write(dir.join("summary.json"), summary).unwrap();
    dir
}

fn full_sections() -> String {
    let p = panels();
    let lpa: Vec<String> = p
        .lpa
        .iter()
        .map(|s| {
            format!(
                r#"{{"rsid":"{}","genotype":"{}/{}","status":"called"}}"#,
                s.rsid, s.ref_, s.ref_
            )
        })
        .collect();
    format!(
        r#"
        "read_count":{{"state":"ok","data":{{}}}},
        "coverage":{{"state":"ok","data":{{"mean_depth":31.5}}}},
        "clinvar":{{"state":"ok","source":"clinvar/{SAMPLE}_clinvar_hits.vcf","data":{{"call_floor":{{"GQ":30,"DP":10}},"hits":[
            {},{},{},{},{},{},{}]}}}},
        "clinical":{{"state":"ok","data":{{"acmg_sf":{{"clinvar_hits":[{{"gene":"BRCA2","variant":"chr13:32320000 {}>{}"}}]}}}}}},
        "pharmcat":{{"state":"ok","source":"pharmcat/{SAMPLE}.report.json","data":{{"version":"3.4.0"}}}},
        "cpic":{{"state":"ok","source":"cpic/{SAMPLE}_phenotypes.tsv","data":{{"genes":[
            {{"gene":"CYP2C19","diplotype":"*38/*38","phenotype":"Normal Metabolizer","status":"normal"}},
            {{"gene":"ABCG2","diplotype":"rs2231142 reference (G)/rs2231142 reference (G)","phenotype":"Normal Function","status":"normal"}},
            {{"gene":"VKORC1","diplotype":"rs9923231 reference (C)/rs9923231 reference (C)","phenotype":"-1639 GG","status":"normal"}},
            {{"gene":"DPYD","diplotype":"Unknown/Unknown","phenotype":"No Result","status":"not called"}},
            {{"gene":"CACNA1S","diplotype":"Reference/Reference","phenotype":"Uncertain Susceptibility","status":"normal"}},
            {{"gene":"MT-RNR1","diplotype":"Reference","phenotype":"normal risk of aminoglycoside-induced hearing loss","status":"normal"}},
            {{"gene":"CYP2D6","diplotype":"Unknown/Unknown","phenotype":"No Result","status":"not called"}},
            {{"gene":"IFNL3","diplotype":"rs12979860 reference (C)/rs12979860 reference (C)","phenotype":"n/a","status":"normal"}},
            {{"gene":"CYP3A4","diplotype":"*1/*1","phenotype":"Normal Metabolizer","status":"normal"}},
            {{"gene":"HLA-B","diplotype":"*57:01/*58:01","phenotype":"HLA-B*58:01 positive","status":"non-normal"}}]}}}},
        "hla":{{"state":"ok","source":"hla_t1k/{SAMPLE}_hla_genotype.tsv","data":{{"loci":[
            {{"gene":"HLA-A","alleles":["HLA-A*02:01:01","HLA-A*24:02:01"],"low_confidence":false}},
            {{"gene":"HLA-B","alleles":["HLA-B*07:02:01","HLA-B*58:01:01"],"low_confidence":false}}],
            "derived":{{"alleles":[
            {{"allele":"B*57:01","kind":"drug","association":"abacavir hypersensitivity","pmids":"18256392","citation":"Mallal 2008","state":"absent"}},
            {{"allele":"B*58:01","kind":"drug","association":"allopurinol severe skin reactions","pmids":"15743917","citation":"Hung 2005","state":"present"}},
            {{"allele":"A*31:01","kind":"drug","association":"carbamazepine","pmids":"21428769","citation":"McCormack 2011","state":"not_determined"}},
            {{"allele":"DRB1*15:01","kind":"disease","association":"multiple sclerosis","pmids":"12557126","citation":"x","state":"present"}}]}}}}}},
        "prs":{{"state":"ok","source":"prs/{SAMPLE}_prs_summary.tsv","data":{{"scores":[
            {{"condition":"Coronary artery disease","pgs_id":"PGS000018","matched":"1900","total":"2000","percentile":"86.0","group":"EUR","evaluated_in":["EUR"]}},
            {{"condition":"Alzheimer's disease","pgs_id":"PGS000334","matched":"80","total":"100","percentile":"40.2","group":"EUR","evaluated_in":["EUR"]}},
            {{"condition":"Height","pgs_id":"PGS000999","matched":"1","total":"2","percentile":"NA","group":"EUR","evaluated_in":[]}}]}}}},
        "traits":{{"state":"ok","source":"traits/{SAMPLE}_traits.tsv","data":{{"apoe":[
            {{"rsid":"rs429358","genotype":"T/C","status":"called"}},
            {{"rsid":"rs7412","genotype":"C/C","status":"called"}},
            {{"rsid":"rs429358+rs7412","genotype":"e3/e4","status":"called","reading":"One e4."}}]}}}},
        "lpa":{{"state":"ok","source":"lpa/{SAMPLE}_lpa_kiv2.tsv","data":{{"status":"estimated","kiv2_copies":"37.0","ci95_low":"36.3","ci95_high":"37.8",
            "risk_alleles":0,"message":"about 37 KIV-2 repeats","flag_text":"No genetic sign of a high Lp(a).","variants":[{}]}}}},
        "callability":{{"state":"ok","data":{{"not_read":[
            {{"gene":"PMS2","reasons":["segdup"],"share":0.81,"status":"not read"}},
            {{"gene":"OR4F5","reasons":["segdup"],"share":0.0,"status":"not read"}},
            {{"gene":"GSTT1","reasons":["ALT-only"],"share":null,"status":"not on this reference"}}],"not_read_owned":[]}}}}
        "#,
        hit("chr7", 117_300_000, "CFTR", "het", "Pathogenic", 3, None),
        hit(
            "chr7",
            117_300_004,
            "CFTR",
            "hom",
            "Likely pathogenic",
            1,
            None
        ),
        hit("chr13", 32_320_000, "BRCA2", "het", "Pathogenic", 2, None),
        hit(
            "chr7",
            117_300_008,
            "CFTR",
            "het",
            "Pathogenic",
            2,
            Some("NN")
        ),
        hit(
            "chr7",
            117_300_012,
            "CFTR",
            "het",
            "Conflicting classifications of pathogenicity",
            1,
            None
        ),
        hit("chr13", 1_000_000, "BRCA2", "het", "Pathogenic", 2, None),
        hit("chr2", 1_000, "", "het", "Pathogenic", 2, None),
        snv(32_320_000).0,
        snv(32_320_000).1,
        lpa.join(","),
    )
}

fn genome() -> Synth {
    let p = panels();
    let mut m = HashMap::new();
    for s in p
        .apoe
        .iter()
        .chain(&p.lpa)
        .chain(p.pgx_sites.iter().map(|(_, s)| s))
    {
        m.insert((s.chrom.clone(), s.pos - 1), s.ref_.as_bytes()[0]);
    }
    Synth(m)
}

fn run(name: &str, sections: &str) -> (Built, PathBuf) {
    let dir = folder(name, sections);
    let built = build(SAMPLE, &dir, &panels(), &mut genome()).expect("builds");
    (built, dir)
}

fn get<'a>(v: &'a Value, path: &[&str]) -> Option<&'a Value> {
    path.iter().try_fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k),
        _ => None,
    })
}

fn st<'a>(v: &'a Value, path: &[&str]) -> &'a str {
    match get(v, path) {
        Some(Value::Str(t)) => t,
        _ => "",
    }
}

fn arr<'a>(v: &'a Value, key: &str) -> &'a [Value] {
    match get(v, &[key]) {
        Some(Value::Array(a)) => a,
        _ => &[],
    }
}

fn facts(b: &Built) -> &[Value] {
    arr(&b.doc, "facts")
}

fn fact<'a>(b: &'a Built, category: &str, gene: &str) -> Vec<&'a Value> {
    facts(b)
        .iter()
        .filter(|f| st(f, &["category"]) == category && st(f, &["gene"]) == gene)
        .collect()
}

fn coverage<'a>(b: &'a Built, gene: &str) -> Vec<&'a Value> {
    arr(&b.doc, "coverage")
        .iter()
        .filter(|c| st(c, &["id"]) == gene)
        .collect()
}

#[test]
fn ids_are_the_contract_hash_of_category_gene_and_key() {
    let (b, _) = run("ids", &full_sections());
    assert!(facts(&b).len() >= 10, "{}", facts(&b).len());
    for f in facts(&b) {
        let (cat, gene) = (st(f, &["category"]), st(f, &["gene"]));
        let key = match cat {
            "pgx" => gene.to_owned(),
            "hla" => st(f, &["value", "marker"]).to_owned(),
            "prs" => st(f, &["value", "pgs_id"]).to_owned(),
            "risk_marker" => gene.to_owned(),
            _ => st(&arr(f, "variants")[0], &["vrs_id"]).to_owned(),
        };
        assert_eq!(st(f, &["id"]), fact_id(cat, gene, &key), "{cat} {gene}");
    }
}

#[test]
fn every_variant_has_a_vrs_id_and_a_call_state() {
    let (b, _) = run("variants", &full_sections());
    let mut n = 0;
    for f in facts(&b) {
        for v in arr(f, "variants") {
            n += 1;
            assert!(st(v, &["vrs_id"]).starts_with("ga4gh:VA."), "{v:?}");
            assert!(!st(v, &["genotype", "call_state"]).is_empty(), "{v:?}");
        }
    }
    assert_eq!(
        n,
        3 + 2 + 2 + 3,
        "three ClinVar hits, APOE's two sites, LPA's two, the sites of ABCG2, VKORC1 and IFNL3"
    );
    let apoe = fact(&b, "risk_marker", "APOE");
    let states: Vec<&str> = arr(apoe[0], "variants")
        .iter()
        .map(|v| st(v, &["genotype", "call_state"]))
        .collect();
    assert_eq!(states, ["alt", "ref"]);
}

#[test]
fn clinvar_hits_are_carrier_acmg_sf_or_monogenic() {
    let (b, _) = run("clinvar", &full_sections());
    let cftr = |cat| fact(&b, cat, "CFTR");
    assert_eq!(
        cftr("carrier").len(),
        1,
        "het in a recessive carrier-list gene"
    );
    assert_eq!(
        cftr("monogenic_pathogenic").len(),
        1,
        "hom is not a carrier"
    );
    assert_eq!(st(cftr("carrier")[0], &["condition", "inheritance"]), "AR");
    assert_eq!(
        get(
            &arr(cftr("carrier")[0], "variants")[0],
            &["clinvar_variation_id"]
        ),
        Some(&Value::Int(7105))
    );
    assert_eq!(
        st(cftr("carrier")[0], &["evidence", "level"]),
        "E1",
        "3 stars"
    );
    assert_eq!(
        st(cftr("monogenic_pathogenic")[0], &["evidence", "level"]),
        "E2",
        "1 star"
    );
    assert_eq!(
        fact(&b, "acmg_sf", "BRCA2").len(),
        1,
        "on step 23's ACMG SF list"
    );
    assert_eq!(st(cftr("carrier")[0], &["confidence", "level"]), "high");
    // The wrong-reference hit and the conflicting one are left out, the first with a reason.
    assert_eq!(
        fact(&b, "carrier", "CFTR").len() + fact(&b, "monogenic_pathogenic", "CFTR").len(),
        2
    );
    assert!(
        b.warnings
            .iter()
            .any(|w| w.contains("117300008") && w.contains("no VRS identifier"))
    );
}

#[test]
fn sensitive_and_adult_onset_follow_ownhealth_and_stay_in_its_spans() {
    let (b, _) = run("sensitivity", &full_sections());
    let p = panels();
    let brca2 = fact(&b, "acmg_sf", "BRCA2")[0];
    assert_eq!(st(brca2, &["sensitivity", "level"]), "sensitive");
    assert_eq!(st(brca2, &["sensitivity", "reason"]), "acmg_sf_cancer");
    assert_eq!(
        st(
            fact(&b, "risk_marker", "APOE")[0],
            &["sensitivity", "level"]
        ),
        "adult_onset"
    );
    let alz = facts(&b)
        .iter()
        .find(|f| st(f, &["value", "pgs_id"]) == "PGS000334")
        .unwrap();
    assert_eq!(st(alz, &["sensitivity", "level"]), "adult_onset");
    assert_eq!(
        st(fact(&b, "pgx", "CYP2C19")[0], &["sensitivity", "level"]),
        "normal"
    );
    // Importer rule 9: what a sensitive or carrier fact cites lies in OwnHealth's spans.
    for f in facts(&b) {
        let sensitive = st(f, &["sensitivity", "level"]) != "normal";
        for v in arr(f, "variants") {
            let pos = match get(v, &["pos"]) {
                Some(Value::Int(n)) => *n as u64,
                _ => 0,
            };
            if sensitive {
                assert!(p.in_sensitive_span(st(v, &["chrom"]), pos), "{f:?}");
            }
            if st(f, &["category"]) == "carrier" {
                assert!(p.in_carrier_span(st(v, &["chrom"]), pos), "{f:?}");
            }
        }
    }
}

#[test]
fn pgx_calls_carry_terms_and_guidance_and_cyp2d6_its_caller() {
    let (b, _) = run("pgx", &full_sections());
    let c19 = fact(&b, "pgx", "CYP2C19")[0];
    assert_eq!(st(c19, &["value", "diplotype"]), "*38/*38");
    let drugs: Vec<(&str, &str)> = arr(c19, "drug_guidance")
        .iter()
        .map(|g| (st(g, &["drug"]), st(g, &["source"])))
        .collect();
    assert!(drugs.contains(&("clopidogrel", "CPIC")), "{drugs:?}");
    // One-site genes: PharmCAT's label becomes a term; a phenotype that is not one is left out.
    assert_eq!(
        st(fact(&b, "pgx", "ABCG2")[0], &["value", "diplotype"]),
        "rs2231142 G/G"
    );
    let vk = fact(&b, "pgx", "VKORC1")[0];
    assert_eq!(st(vk, &["value", "diplotype"]), "rs9923231 C/C");
    assert!(get(vk, &["value", "phenotype"]).is_none());
    let d6 = fact(&b, "pgx", "CYP2D6")[0];
    assert_eq!(st(d6, &["value", "diplotype"]), "*4/*5");
    assert_eq!(st(d6, &["confidence", "level"]), "medium", "one caller");
    assert!(fact(&b, "pgx", "DPYD").is_empty(), "not called: no fact");
    assert!(fact(&b, "pgx", "HLA-B").is_empty(), "HLA comes from T1K");
    for f in facts(&b).iter().filter(|f| st(f, &["category"]) == "pgx") {
        let gene = st(f, &["gene"]);
        for k in ["diplotype", "phenotype"] {
            if let Some(Value::Str(t)) = get(f, &["value", k]) {
                assert!(term_ok(t, Some(gene)), "{gene} {k} {t}");
            }
        }
    }
}

#[test]
fn hla_prs_and_lpa() {
    let (b, _) = run("other", &full_sections());
    let hla = fact(&b, "hla", "HLA-B");
    let present: Vec<(&str, Option<&Value>)> = hla
        .iter()
        .map(|f| (st(f, &["value", "marker"]), get(f, &["value", "present"])))
        .collect();
    assert_eq!(
        present,
        [
            ("HLA-B*57:01", Some(&Value::Bool(false))),
            ("HLA-B*58:01", Some(&Value::Bool(true)))
        ],
        "drug alleles typed; not determined and disease alleles left out"
    );
    assert!(
        fact(&b, "hla", "HLA-A").is_empty(),
        "A*31:01 was not determined"
    );
    let prs: Vec<&str> = facts(&b)
        .iter()
        .filter(|f| st(f, &["category"]) == "prs")
        .map(|f| st(f, &["value", "pgs_id"]))
        .collect();
    assert_eq!(
        prs,
        ["PGS000018", "PGS000334"],
        "a score without a percentile is left out"
    );
    let cad = facts(&b)
        .iter()
        .find(|f| st(f, &["value", "pgs_id"]) == "PGS000018")
        .unwrap();
    assert_eq!(
        get(cad, &["value", "matched_pct"]),
        Some(&Value::Float(95.0))
    );
    assert_eq!(st(cad, &["confidence", "level"]), "medium");
    assert_eq!(
        st(fact(&b, "risk_marker", "LPA")[0], &["evidence", "level"]),
        "E2"
    );
}

#[test]
fn coverage_names_every_panel_gene_once_with_its_state() {
    let (b, _) = run("coverage", &full_sections());
    let p = panels();
    let ids: Vec<&str> = arr(&b.doc, "coverage")
        .iter()
        .map(|c| st(c, &["id"]))
        .collect();
    let mut want: Vec<&str> = p.acmg.iter().map(|(g, _)| g.as_str()).collect();
    want.extend(p.cpic.iter().map(String::as_str));
    want.extend(p.carrier.iter().map(|c| c.span.gene.as_str()));
    for g in &want {
        assert_eq!(ids.iter().filter(|i| *i == g).count(), 1, "{g}");
    }
    let one = |g: &str| coverage(&b, g)[0];
    let pms2 = one("PMS2");
    assert_eq!(
        (
            st(pms2, &["panel"]),
            st(pms2, &["state"]),
            st(pms2, &["reason"])
        ),
        ("acmg_sf_v3.3", "partial", "paralog")
    );
    assert_eq!(get(pms2, &["fraction_callable"]), Some(&Value::Float(0.81)));
    assert_eq!(st(one("DPYD"), &["reason"]), "required_positions_missing");
    assert_eq!(st(one("CYP2D6"), &["state"]), "partial");
    assert_eq!(
        st(one("COMT"), &["reason"]),
        "not_run",
        "PharmCAT does not call COMT"
    );
    // A call the file leaves out makes its gene partial, never "assessed".
    let hla_a = one("HLA-A");
    assert_eq!(
        (st(hla_a, &["state"]), st(hla_a, &["reason"])),
        ("partial", "other")
    );
    assert!(st(hla_a, &["note"]).contains("HLA-A*31:01"), "{hla_a:?}");
    let brca2 = one("BRCA2");
    assert_eq!(st(brca2, &["state"]), "partial");
    assert!(st(brca2, &["note"]).contains("chr13:1000000"), "{brca2:?}");
    assert!(st(one("CFTR"), &["note"]).contains("117300008"));
    let nameless = one("chr2:1000");
    assert_eq!(
        (st(nameless, &["scope"]), st(nameless, &["state"])),
        ("locus", "partial")
    );
    assert_eq!(st(one("BRCA1"), &["state"]), "assessed");
    assert_eq!(st(one("MT-RNR1"), &["state"]), "assessed");
    assert_eq!(
        (st(one("OR4F5"), &["panel"]), st(one("OR4F5"), &["state"])),
        ("callability", "not_assessed")
    );
    assert_eq!(st(one("GSTT1"), &["reason"]), "other");
    // Every gene with a PGx or HLA call has a row, never "not assessed"
    // beside its call: CYP3A4 is on no panel.
    for f in facts(&b) {
        if matches!(st(f, &["category"]), "pgx" | "hla") {
            let g = st(f, &["gene"]);
            assert_eq!(coverage(&b, g).len(), 1, "{g}");
            assert_ne!(st(one(g), &["state"]), "not_assessed", "{g}");
        }
    }
    let cyp3a4 = one("CYP3A4");
    assert_eq!(
        (st(cyp3a4, &["panel"]), st(cyp3a4, &["state"])),
        ("other", "assessed")
    );
}

#[test]
fn a_section_the_run_did_not_do_is_not_assessed_never_negative() {
    let (b, _) = run(
        "missing",
        r#""read_count":{"state":"ok","data":{}},"clinvar":{"state":"missing","data":{}},"prs":{"state":"missing","data":{}}"#,
    );
    assert!(facts(&b).is_empty(), "{:?}", facts(&b));
    let brca1 = coverage(&b, "BRCA1")[0];
    assert_eq!(
        (st(brca1, &["state"]), st(brca1, &["reason"])),
        ("not_assessed", "not_run")
    );
    assert_eq!(st(coverage(&b, "CYP2C19")[0], &["reason"]), "not_run");
    for id in ["APOE", "LPA", "prs", "callability"] {
        assert_eq!(st(coverage(&b, id)[0], &["state"]), "not_assessed", "{id}");
    }
}

#[test]
fn files_and_provenance() {
    let (b, _) = run("files", &full_sections());
    let files = arr(&b.doc, "files");
    assert_eq!(
        files.len(),
        1,
        "the gVCF has no index here, so it is not listed"
    );
    assert_eq!(st(&files[0], &["role"]), "calls_vcf");
    assert_eq!(st(&files[0], &["sha256"]), sha256_hex(b"calls"));
    assert_eq!(st(&files[0], &["index_sha256"]), sha256_hex(b"calls index"));
    assert_eq!(
        st(&b.doc, &["provenance", "generated_utc"]),
        "2026-10-01T09:15:00Z"
    );
    assert_eq!(st(&b.doc, &["provenance", "version"]), "1.4.2");
    assert_eq!(
        st(&b.doc, &["provenance", "db_versions", "clinvar"]),
        "2026-09-28"
    );
    assert_eq!(st(&b.doc, &["subject", "declared_sex"]), "female");
    assert_eq!(st(&b.doc, &["assay", "type"]), "wgs");
}

#[test]
fn the_file_validates_against_ownhealths_schema() {
    let (b, _) = run("schema", &full_sections());
    let schema_text =
        std::fs::read_to_string(repo().join("assets/ownhealth/genome-facts.v1.json")).unwrap();
    let schema = json::loads(&schema_text).unwrap();
    let errors = og_dev::schema::validate(&schema, &b.doc).unwrap();
    assert!(errors.is_empty(), "{errors:?}");
    // A positive control: the same document with one fact id broken.
    let mut bad = b.doc.clone();
    if let Value::Object(m) = &mut bad
        && let Some(Value::Array(f)) = m.get_mut("facts")
        && let Value::Object(first) = &mut f[0]
    {
        first.insert("id".into(), Value::Str("f_nothex".into()));
    }
    assert_eq!(og_dev::schema::validate(&schema, &bad).unwrap().len(), 1);
}

#[test]
fn the_vendored_files_are_ownhealths_byte_for_byte() {
    let dir = repo().join("assets/ownhealth");
    let table = std::fs::read_to_string(dir.join("VENDORED.tsv")).unwrap();
    let mut n = 0;
    for line in table.lines().filter(|l| !l.starts_with('#')).skip(1) {
        let c: Vec<&str> = line.split('\t').collect();
        let bytes = std::fs::read(dir.join(c[0])).unwrap();
        assert_eq!(sha256_hex(&bytes), c[3], "{}", c[0]);
        n += 1;
    }
    assert_eq!(n, 2);
}

#[test]
fn pgx_and_hla_facts_carry_their_call_state() {
    let (b, _) = run("call-state", &full_sections());
    let state = |gene: &str| -> String {
        let f = fact(&b, "pgx", gene);
        assert_eq!(f.len(), 1, "{gene}");
        st(f[0], &["value", "call_state"]).to_owned()
    };
    // PharmCAT's reference alleles: Reference, *38 for CYP2C19, rsN reference (X).
    assert_eq!(state("CACNA1S"), "ref");
    assert_eq!(state("CYP2C19"), "ref");
    assert_eq!(state("ABCG2"), "ref");
    assert_eq!(state("VKORC1"), "ref");
    assert_eq!(state("CYP2D6"), "alt", "*4/*5");
    // A call whose words OwnHealth does not read is kept with its call_state.
    assert_eq!(state("MT-RNR1"), "ref");
    let mt = fact(&b, "pgx", "MT-RNR1")[0];
    assert!(
        get(mt, &["value", "diplotype"]).is_none() && get(mt, &["value", "phenotype"]).is_none()
    );
    // A one-site gene cites its site, which OwnHealth reads the call from.
    let vk = arr(fact(&b, "pgx", "VKORC1")[0], "variants");
    assert_eq!(vk.len(), 1);
    assert_eq!(
        (
            st(&vk[0], &["rsid"]),
            st(&vk[0], &["genotype", "call_state"]),
            st(&vk[0], &["genotype", "gt"])
        ),
        ("rs9923231", "ref", "0/0")
    );
    assert!(st(&vk[0], &["vrs_id"]).starts_with("ga4gh:VA."));
    assert!(arr(fact(&b, "pgx", "CYP2C19")[0], "variants").is_empty());
    // PharmCAT's placeholder phenotype is no phenotype: the call keeps its
    // call_state, and the summary does not show the word.
    let ifnl3 = fact(&b, "pgx", "IFNL3")[0];
    assert_eq!(state("IFNL3"), "ref");
    assert!(get(ifnl3, &["value", "phenotype"]).is_none());
    assert!(
        !st(ifnl3, &["summary"]).contains("n/a"),
        "{}",
        st(ifnl3, &["summary"])
    );
    let hla: Vec<(&str, &str)> = fact(&b, "hla", "HLA-B")
        .iter()
        .map(|f| (st(f, &["value", "marker"]), st(f, &["value", "call_state"])))
        .collect();
    assert_eq!(hla, [("HLA-B*57:01", "ref"), ("HLA-B*58:01", "alt")]);
    assert_eq!(og_facts::build::pgx_call_state("NAT2", "*4/*4"), "ref");
    assert_eq!(og_facts::build::pgx_call_state("TPMT", "*1/*3A"), "alt");
    assert_eq!(
        og_facts::build::pgx_call_state("TPMT", "*1/Unknown"),
        "no_call"
    );
}

/// The sections of `full_sections` with `name`'s replaced by `with`, or
/// dropped when `with` is empty.
fn swap(name: &str, with: &str) -> String {
    let all = full_sections();
    let key = format!("\"{name}\":{{\"state\"");
    let at = all.find(&key).unwrap_or_else(|| panic!("{name}"));
    // The section runs to the comma before the next top-level key, or the end.
    let rest = &all[at..];
    let mut depth = 0i32;
    let mut end = rest.len();
    for (i, c) in rest.char_indices() {
        match c {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    end = i + 1;
                    break;
                }
            }
            _ => {}
        }
    }
    let tail = rest[end..]
        .trim_start()
        .strip_prefix(',')
        .unwrap_or(&rest[end..]);
    let head = all[..at].trim_end();
    match (with.is_empty(), tail.trim().is_empty()) {
        (true, true) => head.trim_end_matches(',').to_owned(),
        (true, false) => format!("{head}{tail}"),
        (false, true) => format!("{head}{with}"),
        (false, false) => format!("{head}{with},{tail}"),
    }
}

#[test]
fn a_gene_is_assessed_only_with_evidence_that_it_was_read() {
    let low = swap(
        "coverage",
        r#""coverage":{"state":"ok","data":{"mean_depth":0.11}}"#,
    );
    let (b, _) = run("low-depth", &low);
    let brca1 = coverage(&b, "BRCA1")[0];
    assert_eq!(
        (st(brca1, &["state"]), st(brca1, &["reason"])),
        ("not_assessed", "low_depth")
    );
    assert!(st(brca1, &["note"]).contains("0.11"), "{brca1:?}");
    let (b, _) = run("no-callability", &swap("callability", ""));
    let brca1 = coverage(&b, "BRCA1")[0];
    assert_eq!(
        (st(brca1, &["state"]), st(brca1, &["reason"])),
        ("not_assessed", "not_run")
    );
    // Callability ran but the coverage step gave no mean depth: not proven either.
    let (b, _) = run(
        "no-depth",
        &swap("coverage", r#""coverage":{"state":"ok","data":{}}"#),
    );
    let brca1 = coverage(&b, "BRCA1")[0];
    assert_eq!(
        (st(brca1, &["state"]), st(brca1, &["reason"])),
        ("not_assessed", "not_run")
    );
    assert!(st(brca1, &["note"]).contains("mean depth"), "{brca1:?}");
    let (b, _) = run("read", &full_sections());
    assert_eq!(st(coverage(&b, "BRCA1")[0], &["state"]), "assessed");
}

#[test]
fn a_call_never_sits_beside_not_assessed_nor_assessed_with_an_allele_missing() {
    // Without step 16b, CACNA1S and RYR1 (on both the ACMG SF and CPIC
    // lists) take the screen's "not proven read" row; their PGx call makes
    // them partial, and the note keeps both parts.
    let (b, _) = run("call-no-callability", &swap("callability", ""));
    for f in facts(&b) {
        if matches!(st(f, &["category"]), "pgx" | "hla") {
            let g = st(f, &["gene"]);
            let row = coverage(&b, g)[0];
            assert!(
                matches!(st(row, &["state"]), "assessed" | "partial"),
                "{g}: {row:?}"
            );
        }
    }
    let cacna1s = coverage(&b, "CACNA1S")[0];
    assert_eq!(
        (st(cacna1s, &["state"]), st(cacna1s, &["reason"])),
        ("partial", "other")
    );
    let note = st(cacna1s, &["note"]);
    assert!(
        note.contains("PGx call") && note.contains("step 16b"),
        "{note}"
    );
    // PharmCAT could not call one allele: partial, on a panel gene and on a
    // gene no panel names.
    let one_unknown = full_sections()
        .replace(
            r#""gene":"CYP2C19","diplotype":"*38/*38""#,
            r#""gene":"CYP2C19","diplotype":"*38/Unknown""#,
        )
        .replace(
            r#""gene":"CYP3A4","diplotype":"*1/*1""#,
            r#""gene":"CYP3A4","diplotype":"*1/Unknown""#,
        );
    let (b, _) = run("call-one-unknown", &one_unknown);
    for g in ["CYP2C19", "CYP3A4"] {
        assert_eq!(
            st(fact(&b, "pgx", g)[0], &["value", "call_state"]),
            "no_call"
        );
        let row = coverage(&b, g)[0];
        assert_eq!(
            (st(row, &["state"]), st(row, &["reason"])),
            ("partial", "other"),
            "{g}"
        );
        assert!(st(row, &["note"]).contains("one allele"), "{row:?}");
    }
}

#[test]
fn a_run_from_a_bam_has_reads_and_a_vcf_only_run_is_capped() {
    // No read count (alignment was skipped), but mosdepth and callability ran.
    let (b, _) = run("from-bam", &swap("read_count", ""));
    assert_eq!(st(&b.doc, &["assay", "type"]), "wgs");
    assert_eq!(
        st(fact(&b, "carrier", "CFTR")[0], &["confidence", "level"]),
        "high"
    );
    let vcf_only = swap("callability", "");
    let vcf_only = vcf_only
        .replace(r#""read_count":{"state":"ok","data":{}},"#, "")
        .replace(
            r#""coverage":{"state":"ok","data":{"mean_depth":31.5}},"#,
            "",
        );
    let (b, _) = run("vcf-only", &vcf_only);
    assert_eq!(st(&b.doc, &["assay", "type"]), "vcf_only");
    let cftr = fact(&b, "carrier", "CFTR")[0];
    assert_eq!(st(cftr, &["confidence", "level"]), "medium");
    assert!(st(cftr, &["confidence", "basis"]).contains("VCF given as input"));
    assert_eq!(
        st(fact(&b, "pgx", "CYP2D6")[0], &["confidence", "level"]),
        "low"
    );
    assert_eq!(st(coverage(&b, "BRCA1")[0], &["state"]), "not_assessed");
}

#[test]
fn lpa_states_only_what_was_read() {
    let (b, _) = run("lpa", &full_sections());
    let lpa = fact(&b, "risk_marker", "LPA")[0];
    assert_eq!(get(lpa, &["value", "risk_alleles"]), Some(&Value::Int(0)));
    assert!(!st(lpa, &["summary"]).contains("Measure"), "{lpa:?}");
    assert!(
        st(lpa, &["summary"]).contains("KIV-2 about 37 copies (95% interval 36.3 to 37.8)"),
        "{lpa:?}"
    );
    // Step 46 ran but read neither the KIV-2 count nor the sites (no traits table).
    let unread = swap(
        "lpa",
        r#""lpa":{"state":"ok","source":"lpa/x.tsv","data":{"status":"not assessed","kiv2_copies":".","risk_alleles":0,"variants":[]}}"#,
    );
    let (b, _) = run("lpa-unread", &unread);
    assert!(fact(&b, "risk_marker", "LPA").is_empty());
    let row = coverage(&b, "LPA")[0];
    assert_eq!(
        (st(row, &["state"]), st(row, &["reason"])),
        ("not_assessed", "other")
    );
    // The KIV-2 count without the traits table: no sites, no risk allele count.
    let kiv2 = swap(
        "lpa",
        r#""lpa":{"state":"ok","source":"lpa/x.tsv","data":{"status":"estimated","kiv2_copies":"37.0","risk_alleles":0,"variants":[]}}"#,
    );
    let (b, _) = run("lpa-kiv2", &kiv2);
    let lpa = fact(&b, "risk_marker", "LPA")[0];
    assert!(get(lpa, &["value", "risk_alleles"]).is_none(), "{lpa:?}");
    assert!(get(lpa, &["variants"]).is_none(), "{lpa:?}");
    // A site with no VRS identifier is not cited, and the coverage map says so.
    let mut p = panels();
    let rsid = p.lpa[0].rsid.clone();
    p.lpa[0].chrom = "chrUn".to_owned();
    let dir = folder("lpa-no-vrs", &full_sections());
    let b = build(SAMPLE, &dir, &p, &mut genome()).unwrap();
    let lpa = fact(&b, "risk_marker", "LPA")[0];
    assert!(!format!("{:?}", get(lpa, &["variants"])).contains(&rsid));
    let row = coverage(&b, "LPA")[0];
    assert_eq!(st(row, &["state"]), "partial", "{row:?}");
    assert!(st(row, &["note"]).contains(&rsid), "{row:?}");
}

#[test]
fn cyp2d6_is_high_only_when_two_callers_agree() {
    let called = full_sections().replace(
        r#"{"gene":"CYP2D6","diplotype":"Unknown/Unknown","phenotype":"No Result","status":"not called"}"#,
        r#"{"gene":"CYP2D6","diplotype":"*1/*4","phenotype":"Intermediate Metabolizer","status":"non-normal"}"#,
    );
    let disagree = r#"{"calls":{"PharmCAT":"*1/*4","pypgx":"*1/*4","Cyrius":"*1/*1"},"agree":false,
        "consensus":{"result":"*1/*4","passed_to_pharmcat":true,"reason":"x","source":null},"probable":null}"#;
    let dir = folder_with("d6-disagree", &called, disagree);
    let b = build(SAMPLE, &dir, &panels(), &mut genome()).unwrap();
    assert_eq!(
        st(fact(&b, "pgx", "CYP2D6")[0], &["confidence", "level"]),
        "medium"
    );
    // Two callers agree, but step 36 withheld the call (its depth check failed): no fact.
    let withheld = r#"{"calls":{"PharmCAT":"not called","pypgx":"*1/*1","Cyrius":"*1/*1"},"agree":true,
        "consensus":{"result":"indeterminate","passed_to_pharmcat":false,"reason":"the depth check failed","source":null},"probable":null}"#;
    let dir = folder_with("d6-withheld", &full_sections(), withheld);
    let b = build(SAMPLE, &dir, &panels(), &mut genome()).unwrap();
    assert!(fact(&b, "pgx", "CYP2D6").is_empty());
    // The two read-level callers agree on another call than PharmCAT's, and
    // step 36 did not hand theirs to PharmCAT: the fact shows PharmCAT's call, unconfirmed.
    let other = r#"{"calls":{"PharmCAT":"*1/*4","pypgx":"*1/*1","Cyrius":"*1/*1"},"agree":true,
        "consensus":{"result":"*1/*1","passed_to_pharmcat":false,"reason":"x","source":null},"probable":null}"#;
    let dir = folder_with("d6-other", &called, other);
    let b = build(SAMPLE, &dir, &panels(), &mut genome()).unwrap();
    let d6 = fact(&b, "pgx", "CYP2D6")[0];
    assert_eq!(st(d6, &["confidence", "level"]), "medium");
    assert_eq!(st(d6, &["value", "diplotype"]), "*1/*4");
    // Not handed over, but PharmCAT called the same pair (in the other order): confirmed.
    let same = r#"{"calls":{"PharmCAT":"*1/*4","pypgx":"*4/*1","Cyrius":"*4/*1"},"agree":true,
        "consensus":{"result":"*4/*1","passed_to_pharmcat":false,"reason":"pypgx and Cyrius agree","source":null},"probable":null}"#;
    let dir = folder_with("d6-same", &called, same);
    let b = build(SAMPLE, &dir, &panels(), &mut genome()).unwrap();
    assert_eq!(
        st(fact(&b, "pgx", "CYP2D6")[0], &["confidence", "level"]),
        "high"
    );
    let agree = r#"{"calls":{"PharmCAT":"*1/*4","pypgx":"*1/*4","Cyrius":"*1/*4"},"agree":true,
        "consensus":{"result":"*1/*4","passed_to_pharmcat":true,"reason":"pypgx and Cyrius agree","source":null},"probable":null}"#;
    let dir = folder_with("d6-agree", &called, agree);
    let b = build(SAMPLE, &dir, &panels(), &mut genome()).unwrap();
    assert_eq!(
        st(fact(&b, "pgx", "CYP2D6")[0], &["confidence", "level"]),
        "high"
    );
    // PharmCAT's table says "not called", but step 36 handed it the agreed
    // call: the fact is the agreed call, and its gene is assessed, never
    // "not assessed" beside the call.
    let handed = r#"{"calls":{"PharmCAT":"not called","pypgx":"*1/*4","Cyrius":"*1/*4"},"agree":true,
        "consensus":{"result":"*1/*4","passed_to_pharmcat":true,"reason":"pypgx and Cyrius agree","source":null},"probable":null}"#;
    let dir = folder_with("d6-handed", &full_sections(), handed);
    let b = build(SAMPLE, &dir, &panels(), &mut genome()).unwrap();
    assert_eq!(
        st(fact(&b, "pgx", "CYP2D6")[0], &["value", "diplotype"]),
        "*1/*4"
    );
    assert_eq!(st(coverage(&b, "CYP2D6")[0], &["state"]), "assessed");
    // A probable call from one caller stays partial.
    let (b, _) = run("d6-probable", &full_sections());
    assert_eq!(st(coverage(&b, "CYP2D6")[0], &["state"]), "partial");
}

#[test]
fn a_reference_without_grch38_names_is_refused() {
    let dir = std::env::temp_dir().join(format!("og-facts-{}-b37", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fa = dir.join("ref.fa");
    std::fs::write(&fa, b">1\nACGT\n").unwrap();
    std::fs::write(dir.join("ref.fa.fai"), b"1\t4\t3\t4\t5\n").unwrap();
    let idx = og_intervals::fasta::Indexed::open(&fa).unwrap();
    let e = og_facts::check_build(&idx, &panels()).unwrap_err();
    assert!(e.contains("none of GRCh38's chromosomes"), "{e}");
    std::fs::write(&fa, b">chrM\nACGT\n").unwrap();
    std::fs::write(dir.join("ref.fa.fai"), b"chrM\t4\t6\t4\t5\n").unwrap();
    let idx = og_intervals::fasta::Indexed::open(&fa).unwrap();
    assert!(
        og_facts::check_build(&idx, &panels())
            .unwrap_err()
            .contains("16569")
    );
}
