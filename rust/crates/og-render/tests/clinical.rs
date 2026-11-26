//! The clinical sections, end to end: synthetic sample folders read by
//! `og summary collect` (in the process) and written by `og report render`.
//! Nothing here describes a person: the calls are made up at real gene
//! positions.
//!
//! This was three Python files; each test names the checks it carries.
//!
//! - `tests/test_acmg_sf_rules.py`: the ACMG SF v3.3 catalog and its rules
//!   (MUTYH, HFE, the X-linked genes, SVs, BTD's partial deficiency), the
//!   one-person carrier list with its set-aside alleles, both reports, and
//!   a run without the assets folder. Its couple screen checks are in
//!   `rust/crates/og-couple/tests/port_report_clinical.rs`.
//! - `tests/test_cross_engine_biallelic.py`: the two-hit rule across
//!   callers (a deletion plus a small variant), the hemizygous genotype, the
//!   SMN1 and HBA pairs, the gap line when no deletion caller was read, and
//!   both reports. Its couple screen checks are in og-couple too.
//! - `tests/test_callability.py`: the callability catalog, step 16b's
//!   reader with its planted-bug controls, and both reports. Its ALT contig
//!   refusal checks are in `rust/crates/og-cn/tests/port_report_clinical.rs`.

use std::fs;
use std::io::Write;
use std::path::PathBuf;

use og_render::model::Summary;
use og_render::sections::{acmg_sf, carrier, clinvar::D};
use og_summary::collect::{Group, Input};
use pycompat::json::{Dump, Value, dumps};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");
const NULL: Value = Value::Null;

// --- scratch folders and sample folders ----------------------------------------------

/// A scratch folder, removed when dropped.
struct Dir(PathBuf);

impl Dir {
    fn new(name: &str) -> Dir {
        let p =
            std::env::temp_dir().join(format!("og-render-clinical-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&p);
        fs::create_dir_all(&p).unwrap();
        Dir(p)
    }

    fn path(&self, rel: &str) -> String {
        self.0.join(rel).to_string_lossy().into_owned()
    }

    fn put(&self, rel: &str, text: &str) -> String {
        let p = self.0.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }

    /// Writes `text` BGZF-compressed.
    fn put_gz(&self, rel: &str, text: &str) -> String {
        let p = self.0.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        let f = fs::File::create(&p).unwrap();
        let mut w = og_io::bgzf::Writer::new(f, og_io::bgzf::workers()).unwrap();
        w.write_all(text.as_bytes()).unwrap();
        w.finish().unwrap();
        p.to_string_lossy().into_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// summary.json of sample `S` in `dir`, as `og summary collect` writes it.
fn collect_with(dir: &str, assets: &str) -> Value {
    og_summary::collect::post::collect(&Input {
        sample: "S".to_owned(),
        sample_dir: dir.to_owned(),
        declared_sex: None,
        assets: assets.to_owned(),
        group: Group::All,
        tables: format!("{REPO}/bin"),
    })
    .unwrap()
}

fn collect(dir: &str) -> Value {
    collect_with(dir, &format!("{REPO}/assets"))
}

/// Both reports of a summary: (text, HTML).
fn reports(v: &Value) -> (String, String) {
    let s = Summary::parse(&dumps(v, &Dump::default()).unwrap()).unwrap();
    (
        og_render::text::report(&s).unwrap(),
        og_render::html::report(&s).unwrap(),
    )
}

/// `v[a][b]...`, null where a key is missing.
fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    let mut v = v;
    for k in path {
        v = match v {
            Value::Object(o) => o.get(*k).unwrap_or(&NULL),
            _ => &NULL,
        };
    }
    v
}

fn list(v: &Value) -> &[Value] {
    match v {
        Value::Array(a) => a,
        _ => &[],
    }
}

fn s(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => "",
    }
}

fn set(v: &mut Value, path: &[&str], to: Value) {
    let mut v = v;
    for k in &path[..path.len() - 1] {
        let Value::Object(o) = v else { panic!("{k}") };
        v = o.get_mut(*k).unwrap();
    }
    let Value::Object(o) = v else { panic!() };
    o.insert(path[path.len() - 1].to_owned(), to);
}

const VCF_HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

/// A ClinVar P hit, reviewed by several submitters; `ps`: a phase set.
fn hit(v: (&str, i64, &str, &str), gene: &str, gt: &str, ps: Option<i64>) -> String {
    let (chrom, pos, r, a) = v;
    let (fmt, val) = match ps {
        Some(ps) => ("GT:GQ:DP:PS", format!("{gt}:50:30:{ps}")),
        None => ("GT:GQ:DP", format!("{gt}:50:30")),
    };
    format!(
        "{chrom}\t{pos}\t.\t{r}\t{a}\t50\tPASS\tGENEINFO={gene}:1;CLNSIG=Pathogenic;\
         CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\t{fmt}\t{val}\n"
    )
}

fn h(chrom: &str, pos: i64, r: &str, a: &str, gene: &str, gt: &str) -> String {
    hit((chrom, pos, r, a), gene, gt, None)
}

const ANNOTSV_COLS: &str = "SV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tGene_name\tACMG_class\tCaller\t\
                            Exon_edge\tExon_edge_detail\tPopulation_AF\tPopulation_source\tCommon";

/// AnnotSV's class 4-5 table as the SF test wrote it: no genotype columns.
fn tab13(rows: &str) -> String {
    format!("{ANNOTSV_COLS}\n{rows}")
}

/// The same table with the genotype of each SV, as the two-hit test wrote it.
fn tab15(rows: &str) -> String {
    format!("{ANNOTSV_COLS}\tFORMAT\tS\n{rows}")
}

/// A class `cls` deletion AnnotSV reports, with its genotype.
fn sv(
    span: (&str, i64, i64),
    genes: &str,
    cls: &str,
    edge: &str,
    gt: &str,
    ps: Option<i64>,
) -> String {
    let (chrom, a, b) = span;
    let (fmt, val) = match ps {
        Some(ps) => ("GT:PS", format!("{gt}:{ps}")),
        None => ("GT", gt.to_owned()),
    };
    format!(
        "{chrom}\t{a}\t{b}\t-{}\tDEL\t{genes}\t{cls}\tmanta\t{edge}\t.\t\t\tno\t{fmt}\t{val}\n",
        b - a
    )
}

/// A sample folder of ClinVar hits and, when given, AnnotSV's class 4-5
/// table. Paths are absolute: `Dir::put` takes them as they are.
fn sample(t: &Dir, name: &str, hits: &[String], annotsv: Option<&str>) -> String {
    let d = t.path(&format!("{name}/S"));
    t.put(
        &format!("{d}/clinvar/S_clinvar_hits.vcf"),
        &format!("{VCF_HEAD}{}", hits.concat()),
    );
    if let Some(table) = annotsv {
        t.put(&format!("{d}/annotsv/S_sv_class4_5_checked.tsv"), table);
    }
    d
}

/// The SF tier of a folder whose clinical VCF has no record (the Nextflow
/// layout: the tier comes from the hits alone).
fn sf(t: &Dir, d: &str) -> Value {
    let vcf = format!("{d}/clinical/S_clinical.vcf.gz");
    if !std::path::Path::new(&vcf).exists() {
        t.put_gz(&vcf, VCF_HEAD);
    }
    at(&collect(d), &["sections", "clinical", "data", "acmg_sf"]).clone()
}

fn genes_of<'a>(sf: &'a Value, key: &str, gene: &str) -> Vec<&'a Value> {
    list(at(sf, &[key]))
        .iter()
        .filter(|x| s(at(x, &["gene"])) == gene)
        .collect()
}

fn reported<'a>(sf: &'a Value, gene: &str) -> Vec<&'a Value> {
    genes_of(sf, "clinvar_hits", gene)
}

fn apart<'a>(sf: &'a Value, gene: &str) -> Vec<&'a Value> {
    genes_of(sf, "not_reported", gene)
}

/// The rows of an asset table, by column name.
fn asset(name: &str) -> Vec<std::collections::HashMap<String, String>> {
    let text = fs::read_to_string(format!("{REPO}/assets/{name}")).unwrap();
    table(&text)
}

fn table(text: &str) -> Vec<std::collections::HashMap<String, String>> {
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let head: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_owned)
        .collect();
    lines
        .map(|l| {
            head.iter()
                .cloned()
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect()
}

// --- the ACMG SF v3.3 catalog ----------------------------------------------------------

/// ACMG SF v3.3 Table 1 (Lee et al., Genet Med 2025;27(8):101454), gene by
/// gene.
const SF_V33: &str = "APC RET BRCA1 BRCA2 PALB2 SDHD SDHAF2 SDHC SDHB MAX TMEM127 BMPR1A SMAD4 TP53 MLH1 \
    MSH2 MSH6 PMS2 MEN1 MUTYH NF2 STK11 PTEN RB1 TSC1 TSC2 VHL WT1 FBN1 TGFBR1 TGFBR2 SMAD3 ACTA2 MYH11 PKP2 \
    DSP DSC2 TMEM43 DSG2 RYR2 CASQ2 TRDN TNNT2 LMNA FLNC TTN BAG3 DES RBM20 TNNC1 PLN COL3A1 LDLR APOB PCSK9 \
    MYH7 MYBPC3 TNNI3 TPM1 MYL3 ACTC1 PRKAG2 MYL2 KCNQ1 KCNH2 SCN5A CALM1 CALM2 CALM3 BTD CYP27A1 GLA HFE OTC \
    GAA ABCD1 ACVRL1 ENG TTR RYR1 CACNA1S HNF1A RPE65 ATP7B";

/// test_acmg_sf_rules.py check 1.
#[test]
fn the_sf_catalog_is_table_1_of_sf_v3_3_with_its_rules() {
    let rows = asset("acmg_sf.tsv");
    let genes: Vec<&str> = rows.iter().map(|r| r["gene"].as_str()).collect();
    let mut uniq = genes.clone();
    uniq.sort_unstable();
    uniq.dedup();
    assert_eq!((genes.len(), uniq.len()), (84, 84), "84 rows, one per gene");
    let mut want: Vec<&str> = SF_V33.split_whitespace().collect();
    want.sort_unstable();
    assert_eq!(want.len(), 84);
    assert_eq!(uniq, want, "the genes are SF v3.3 Table 1's");
    assert!(rows.iter().all(|r| r["hgnc_id"].starts_with("HGNC:")));
    let mut xl: Vec<&str> = rows
        .iter()
        .filter(|r| r["inheritance"] == "XL")
        .map(|r| r["gene"].as_str())
        .collect();
    xl.sort_unstable();
    assert_eq!(xl, ["ABCD1", "GLA", "OTC"]);
    let rules = og_summary::collect::assets::acmg_sf_rules(&format!("{REPO}/assets")).unwrap();
    for (g, text, rule) in [
        (
            "GLA",
            "All hemi, het, homozygous P and LP",
            "het_female_too",
        ),
        (
            "OTC",
            "All hemi, het, homozygous P and LP",
            "het_female_too",
        ),
        (
            "ABCD1",
            "All hemi, homozygous or two het. P and/or LP",
            "hemizygous",
        ),
    ] {
        assert_eq!(
            (rules[g].table_text.as_str(), rules[g].report_if.as_str()),
            (text, rule),
            "{g}"
        );
    }
    let mut two: Vec<&str> = rules
        .iter()
        .filter(|(_, r)| r.report_if == "biallelic")
        .map(|(g, _)| g.as_str())
        .collect();
    two.sort_unstable();
    assert_eq!(
        two,
        [
            "ATP7B", "BTD", "CASQ2", "CYP27A1", "GAA", "MUTYH", "RPE65", "TRDN"
        ]
    );
    let hfe = &rules["HFE"];
    assert_eq!(hfe.report_if, "hom_variant");
    let keys: Vec<_> = hfe.variants.keys().cloned().collect();
    assert_eq!(
        keys,
        [("chr6".to_owned(), 26092913, "G".to_owned(), "A".to_owned())],
        "HFE counts p.C282Y homozygotes only"
    );
    let mut svs: Vec<&str> = rules
        .iter()
        .filter(|(_, r)| r.sv_cn)
        .map(|(g, _)| g.as_str())
        .collect();
    svs.sort_unstable();
    assert_eq!(svs, ["BRCA1", "BRCA2", "LDLR", "MLH1", "MSH2"]);
    // Control: a copy without TTN has 83 genes, so the count can fail.
    let t = Dir::new("sf-dropped");
    let text = fs::read_to_string(format!("{REPO}/assets/acmg_sf.tsv")).unwrap();
    let dropped: String = text
        .split_inclusive('\n')
        .filter(|l| !l.starts_with("TTN\t"))
        .collect();
    t.put("acmg_sf.tsv", &dropped);
    let n = og_summary::collect::assets::acmg_sf_rules(&t.path(""))
        .unwrap()
        .len();
    assert_eq!(n, 83);
    // A severity entry with a kind and no variant is refused.
    assert!(og_summary::collect::assets::parse_alleles(Some("hypomorph"), true).is_err());
}

// --- the SF rules ----------------------------------------------------------------------

// MUTYH p.Gly396Asp and p.Tyr179Cys, HFE p.C282Y (GRCh38).
const G396: (&str, i64, &str, &str) = ("chr1", 45332803, "C", "T");
const Y179: (&str, i64, &str, &str) = ("chr1", 45334483, "T", "C");
const C282Y: (&str, i64, &str, &str) = ("chr6", 26092913, "G", "A");

/// test_acmg_sf_rules.py checks 2 and 3.
#[test]
fn mutyh_needs_two_variants_and_hfe_only_c282y_homozygous() {
    let t = Dir::new("sf-mutyh");
    let sfv = sf(
        &t,
        &sample(&t, "het", &[hit(G396, "MUTYH", "0/1", None)], None),
    );
    assert!(reported(&sfv, "MUTYH").is_empty(), "one het is not an SF");
    let why: Vec<&str> = apart(&sfv, "MUTYH")
        .iter()
        .map(|x| s(at(x, &["reason"])))
        .collect();
    assert_eq!(
        why,
        ["one heterozygous variant; MUTYH is reported with P and LP (2 variants)"]
    );
    let sfv = sf(
        &t,
        &sample(
            &t,
            "two",
            &[
                hit(G396, "MUTYH", "0/1", None),
                hit(Y179, "MUTYH", "0/1", None),
            ],
            None,
        ),
    );
    assert_eq!(
        reported(&sfv, "MUTYH").len(),
        2,
        "two hets are (phase unknown)"
    );
    assert!(apart(&sfv, "MUTYH").is_empty());
    let sfv = sf(
        &t,
        &sample(&t, "hom", &[hit(G396, "MUTYH", "1/1", None)], None),
    );
    assert_eq!(reported(&sfv, "MUTYH").len(), 1, "a homozygote is");

    let sfv = sf(
        &t,
        &sample(&t, "hfe-het", &[hit(C282Y, "HFE", "0/1", None)], None),
    );
    assert!(reported(&sfv, "HFE").is_empty());
    assert_eq!(apart(&sfv, "HFE").len(), 1);
    let sfv = sf(
        &t,
        &sample(&t, "hfe-hom", &[hit(C282Y, "HFE", "1/1", None)], None),
    );
    assert_eq!(reported(&sfv, "HFE").len(), 1);
    let sfv = sf(
        &t,
        &sample(
            &t,
            "h63d",
            &[h("chr6", 26090951, "C", "G", "HFE", "1/1")],
            None,
        ),
    );
    assert!(
        reported(&sfv, "HFE").is_empty(),
        "H63D homozygous is not an SF"
    );
}

/// test_acmg_sf_rules.py check 4.
#[test]
fn x_linked_genes_follow_their_table_1_rule() {
    let t = Dir::new("sf-xl");
    let sfv = sf(
        &t,
        &sample(
            &t,
            "female",
            &[
                h("chrX", 101398000, "C", "T", "GLA", "0/1"),
                h("chrX", 38400000, "G", "A", "OTC", "0/1"),
                h("chrX", 153725000, "C", "T", "ABCD1", "0/1"),
            ],
            None,
        ),
    );
    assert_eq!(reported(&sfv, "GLA").len(), 1, "GLA: het female reported");
    assert_eq!(reported(&sfv, "OTC").len(), 1, "OTC: het female reported");
    assert!(
        reported(&sfv, "ABCD1").is_empty(),
        "ABCD1: het female is not"
    );
    assert_eq!(apart(&sfv, "ABCD1").len(), 1);
    let sfv = sf(
        &t,
        &sample(
            &t,
            "male",
            &[h("chrX", 153725000, "C", "T", "ABCD1", "1")],
            None,
        ),
    );
    assert_eq!(
        reported(&sfv, "ABCD1").len(),
        1,
        "a hemizygous ABCD1 call is"
    );
    let sfv = sf(
        &t,
        &sample(
            &t,
            "two",
            &[
                h("chrX", 153725000, "C", "T", "ABCD1", "0/1"),
                h("chrX", 153726000, "G", "A", "ABCD1", "0/1"),
            ],
            None,
        ),
    );
    assert_eq!(
        reported(&sfv, "ABCD1").len(),
        2,
        "two het ABCD1 variants are"
    );
}

/// AnnotSV's rows as the SF test wrote them: no genotype columns.
fn sv13(chrom: &str, a: i64, b: i64, genes: &str) -> String {
    format!(
        "{chrom}\t{a}\t{b}\t-{}\tDEL\t{genes}\t5\tmanta\tno\t.\t\t\tno\n",
        b - a
    )
}

/// test_acmg_sf_rules.py checks 5 and 6.
#[test]
fn sf_deletions_and_btd_severity() {
    let t = Dir::new("sf-sv");
    let rows = sv13("17", 43070000, 43080000, "BRCA1")
        + &sv13("2", 47380000, 47390000, "EPCAM")
        + &sv13("1", 45330000, 45336000, "MUTYH");
    let sfv = sf(&t, &sample(&t, "sv", &[], Some(&tab13(&rows))));
    let st: Vec<(String, String)> = list(at(&sfv, &["structural"]))
        .iter()
        .map(|x| {
            (
                s(at(x, &["gene"])).to_owned(),
                s(at(x, &["via"])).to_owned(),
            )
        })
        .collect();
    assert!(st.contains(&("BRCA1".into(), String::new())), "{st:?}");
    assert!(
        st.iter().any(|(g, v)| g == "MSH2" && v.contains("EPCAM")),
        "an EPCAM 3' deletion counts for MSH2: {st:?}"
    );
    assert!(
        !st.iter().any(|(g, _)| g == "MUTYH"),
        "MUTYH has no SV rule"
    );

    let sfv = sf(
        &t,
        &sample(
            &t,
            "d444h",
            &[h("chr3", 15645186, "G", "C", "BTD", "1/1")],
            None,
        ),
    );
    let btd = reported(&sfv, "BTD");
    assert_eq!(btd.len(), 1);
    assert!(s(at(btd[0], &["severity"])).starts_with("partial, not profound"));
    let sfv = sf(
        &t,
        &sample(
            &t,
            "btd",
            &[h("chr3", 15647000, "G", "A", "BTD", "1/1")],
            None,
        ),
    );
    let btd = reported(&sfv, "BTD");
    assert_eq!(btd.len(), 1);
    assert!(!og_render::model::truthy(at(btd[0], &["severity"])));
}

/// test_acmg_sf_rules.py check 7 (the one-sample part).
#[test]
fn the_carrier_list_of_one_sample_sets_aside_a_pseudodeficiency_allele() {
    let t = Dir::new("carrier");
    let d = sample(
        &t,
        "c",
        &[
            h("chr7", 117559590, "ATCT", "A", "CFTR", "0/1"),
            h("chr15", 72350584, "G", "A", "HEXA", "0/1"),
            h("chr15", 72346579, "G", "GGATA", "HEXA", "0/1"),
        ],
        None,
    );
    let summ = collect(&d);
    let car = at(&summ, &["sections", "carrier", "data"]);
    let lines: Vec<(&str, &str)> = list(at(car, &["carriers"]))
        .iter()
        .map(|x| (s(at(x, &["gene"])), s(at(x, &["variant"]))))
        .collect();
    assert!(
        lines.contains(&("CFTR", "chr7:117559590 ATCT>A")),
        "{lines:?}"
    );
    assert!(!lines.contains(&("HEXA", "chr15:72350584 G>A")));
    assert!(lines.contains(&("HEXA", "chr15:72346579 G>GGATA")));
    let kinds: Vec<&str> = list(at(car, &["set_aside"]))
        .iter()
        .map(|x| s(at(x, &["kind"])))
        .collect();
    assert_eq!(kinds, ["pseudodeficiency"]);
}

/// test_acmg_sf_rules.py checks 8 and 9.
#[test]
fn both_reports_print_the_sf_verdicts_and_the_carrier_list() {
    let t = Dir::new("sf-report");
    let d = sample(
        &t,
        "r",
        &[
            hit(G396, "MUTYH", "0/1", None),
            hit(C282Y, "HFE", "0/1", None),
            h("chr3", 15645186, "G", "C", "BTD", "1/1"),
            h("chr7", 117559590, "ATCT", "A", "CFTR", "0/1"),
        ],
        Some(&tab13(&sv13("17", 43070000, 43080000, "BRCA1"))),
    );
    t.put_gz("r/S/clinical/S_clinical.vcf.gz", VCF_HEAD);
    let mut summ = collect(&d);
    let (txt, html) = reports(&summ);
    let html_plain = html.replace("&gt;", ">");
    for needle in [
        "MUTYH chr1:45332803 C>T het (ClinVar): not reported, one heterozygous variant",
        "HFE chr6:26092913 G>A het (ClinVar): not reported",
        "partial, not profound",
        "DEL chr17:43070000-43080000",
        "Carrier list (",
    ] {
        assert!(txt.contains(needle), "text: {needle}");
        assert!(html_plain.contains(needle), "HTML: {needle}");
    }
    assert!(
        txt.contains("caller support not known") && !txt.contains("None caller"),
        "without step 22 the SF deletion's caller support is not known"
    );
    set(
        &mut summ,
        &["sections", "carrier", "state"],
        Value::Str("failed".into()),
    );
    set(
        &mut summ,
        &["sections", "carrier", "data"],
        Value::Object(Default::default()),
    );
    let (_, html) = reports(&summ);
    assert!(
        html.contains("Carrier List (this sample)"),
        "a failed carrier section still shows its card"
    );

    // Without the assets folder the hits stay and both lists say why.
    let summ = collect_with(&d, &t.path("no-assets"));
    let cv = at(&summ, &["sections", "clinvar"]);
    assert_eq!(s(at(cv, &["state"])), "ok");
    assert_eq!(at(cv, &["data", "count"]), &Value::Int(4));
    let sfv = at(
        &summ,
        &["sections", "clinical", "data", "acmg_sf", "unavailable"],
    );
    assert!(s(sfv).contains("acmg_sf.tsv"), "{sfv:?}");
    let car = at(&summ, &["sections", "carrier", "data", "unavailable"]);
    assert!(s(car).contains("carrier_genes.tsv"), "{car:?}");
    let (txt, _) = reports(&summ);
    assert!(txt.contains("Secondary findings: not assessed"));
    assert!(txt.contains("Carrier list: not made"));
}

// --- two hits across callers -----------------------------------------------------------

// MUTYH spans chr1:45329163-45340893 (assets/carrier_genes.tsv).
const SNV: (&str, i64, &str, &str) = ("chr1", 45332000, "C", "T");
const DEL: (&str, i64, i64) = ("chr1", 45335000, 45337000);
const INSIDE: (&str, i64, &str, &str) = ("chr1", 45336000, "G", "A");
// RBM8A c.-21G>A (rs139428292), on the plus strand C>T.
const TAR_ALLELE: (&str, i64, &str, &str) = ("chr1", 145927447, "C", "T");

fn del(gt: &str, ps: Option<i64>) -> String {
    tab15(&sv(DEL, "MUTYH", "4", "yes", gt, ps))
}

fn carrier_of(d: &str) -> Value {
    at(&collect(d), &["sections", "carrier", "data"]).clone()
}

fn two(d: &str, gene: &str) -> Vec<Value> {
    list(at(&carrier_of(d), &["two_hits", "genes"]))
        .iter()
        .filter(|x| s(at(x, &["gene"])) == gene)
        .cloned()
        .collect()
}

fn field<'a>(v: &'a [Value], k: &str) -> &'a str {
    v.first().map_or("", |x| s(at(x, &[k])))
}

/// Step 35's tables: Parascopy's SMN table, the interval table and the HBA
/// count.
fn step35(t: &Dir, d: &str, smn1: &str, intervals: &str, hba: Option<&str>) {
    let rel = d;
    t.put(
        &format!("{rel}/paralogs/S_smn_copy_number.tsv"),
        &format!(
            "chrom\tstart\tend\tagCN\tpsCN\tpsCN_qual\tpsCN_filter\tagCN_qual\tagCN_filter\thomologous_regions\n\
             chr5\t70951000\t70952000\t4\t{smn1},2\t40,40\tPASS\t40\tPASS\tchr5:70076000-70077000\n"
        ),
    );
    if !intervals.is_empty() {
        t.put(
            &format!("{rel}/paralogs/S_cn_intervals.tsv"),
            &format!(
                "locus\tname\tchrom\tstart\tend\tstatus\tcopies\texpected\tchange\testimate\treason\n{intervals}"
            ),
        );
    }
    if let Some(g) = hba {
        t.put(
            &format!("{rel}/paralogs/S_hba_copy_number.tsv"),
            &format!(
                "metric\tvalue\nstatus\tcalled\nalpha_globin_genes\t{g}\nestimate\t{g}.02\n\
                 message\tone gene lost on one chromosome\n"
            ),
        );
    }
}

/// test_cross_engine_biallelic.py checks 1 to 3.
#[test]
fn a_deletion_and_a_hit_in_one_gene_are_two_hits_unless_in_cis() {
    let t = Dir::new("two-mutyh");
    let d = sample(
        &t,
        "m",
        &[hit(SNV, "MUTYH", "0/1", None)],
        Some(&del("0/1", None)),
    );
    let tw = two(&d, "MUTYH");
    assert_eq!(tw.len(), 1);
    assert_eq!(field(&tw, "verdict"), "possibly biallelic");
    assert!(field(&tw, "deletion").contains("DEL chr1:45335000-45337000"));
    assert_eq!(field(&tw, "variant"), "chr1:45332000 C>T");
    assert!(field(&tw, "why").contains("possibly in trans"));
    let sfv = sf(&t, &d);
    let rep = reported(&sfv, "MUTYH");
    assert_eq!(rep.len(), 1, "the SF tier reports the hit");
    assert!(s(at(rep[0], &["two_hit"])).contains("possibly biallelic with DEL"));
    assert!(apart(&sfv, "MUTYH").is_empty());

    let d = sample(
        &t,
        "trans",
        &[hit(SNV, "MUTYH", "0|1", Some(45330000))],
        Some(&del("1|0", Some(45330000))),
    );
    let tw = two(&d, "MUTYH");
    assert_eq!(tw.len(), 1);
    assert_eq!(field(&tw, "verdict"), "biallelic");
    assert!(field(&tw, "why").contains("in trans"));
    let d = sample(
        &t,
        "otherps",
        &[hit(SNV, "MUTYH", "0|1", Some(45330000))],
        Some(&del("1|0", Some(45000000))),
    );
    assert_eq!(field(&two(&d, "MUTYH"), "verdict"), "possibly biallelic");

    let d = sample(
        &t,
        "cis",
        &[hit(SNV, "MUTYH", "1|0", Some(45330000))],
        Some(&del("1|0", Some(45330000))),
    );
    assert!(two(&d, "MUTYH").is_empty(), "in cis: no two-hit line");
    let sfv = sf(&t, &d);
    assert!(reported(&sfv, "MUTYH").is_empty());
    assert!(
        !apart(&sfv, "MUTYH").is_empty(),
        "in cis: the one hit stays apart"
    );
}

/// test_cross_engine_biallelic.py check 4.
#[test]
fn the_tar_deletion_and_the_rbm8a_allele_give_the_tar_line() {
    let t = Dir::new("two-tar");
    let row = "1q21.1_TAR_syndrome_proximal_BP2-BP3\t1q21.1_TAR_syndrome_proximal_BP2-BP3_t1\tchr1\t145686999\t\
               146048495\tcalled\t1\t2\tdeletion\t1.04\t.\n";
    let d = sample(&t, "tar", &[hit(TAR_ALLELE, "RBM8A", "1/1", None)], None);
    step35(&t, &d, "2", row, None);
    let c = carrier_of(&d);
    let tw: Vec<&Value> = list(at(&c, &["two_hits", "genes"]))
        .iter()
        .filter(|x| s(at(x, &["gene"])) == "RBM8A")
        .collect();
    assert_eq!(tw.len(), 1);
    assert_eq!(s(at(tw[0], &["verdict"])), "possibly biallelic");
    let lines = carrier::lines(&D::new(&c, "carrier")).unwrap();
    let line: Vec<&String> = lines
        .iter()
        .filter(|x| x.contains("RBM8A") && x.contains("both copies"))
        .collect();
    assert_eq!(line.len(), 1, "{lines:?}");
    for part in [
        "Thrombocytopenia-absent radius syndrome",
        "1q21.1_TAR_syndrome",
        "rs139428292",
    ] {
        assert!(line[0].contains(part), "{part}: {}", line[0]);
    }
    let gene_in = |k: &str| {
        list(at(&c, &[k]))
            .iter()
            .any(|x| s(at(x, &["gene"])) == "RBM8A")
    };
    assert!(gene_in("set_aside") && !gene_in("carriers"));
    let d = sample(&t, "alone", &[hit(TAR_ALLELE, "RBM8A", "0/1", None)], None);
    step35(
        &t,
        &d,
        "2",
        &row.replace("\t1\t2\tdeletion\t1.04", "\t2\t2\tno change\t2.01"),
        None,
    );
    assert!(two(&d, "RBM8A").is_empty(), "the allele without a deletion");
}

/// test_cross_engine_biallelic.py check 5.
#[test]
fn a_homozygous_call_inside_a_deletion_is_hemizygous() {
    let t = Dir::new("two-hemi");
    let d = sample(
        &t,
        "hemi",
        &[hit(INSIDE, "MUTYH", "1/1", None)],
        Some(&del("0/1", None)),
    );
    let c = carrier_of(&d);
    let item: Vec<&Value> = list(at(&c, &["hom_or_hemizygous"]))
        .iter()
        .filter(|x| s(at(x, &["gene"])) == "MUTYH")
        .collect();
    assert_eq!(item.len(), 1);
    assert_eq!(s(at(item[0], &["genotype"])), "hemizygous");
    assert!(s(at(item[0], &["hemizygous"])).contains("inside DEL"));
    let tw = two(&d, "MUTYH");
    assert_eq!(field(&tw, "verdict"), "biallelic");
    assert_eq!(field(&tw, "genotype"), "hemizygous");
    let rep = sf(&t, &d);
    assert_eq!(
        s(at(reported(&rep, "MUTYH")[0], &["genotype"])),
        "hemizygous"
    );

    let genotype = |name: &str, v: (&str, i64, &str, &str)| {
        let d = sample(
            &t,
            name,
            &[hit(v, "MUTYH", "1/1", None)],
            Some(&del("0/1", None)),
        );
        let c = carrier_of(&d);
        list(at(&c, &["hom_or_hemizygous"]))
            .iter()
            .find(|x| s(at(x, &["gene"])) == "MUTYH")
            .map(|x| s(at(x, &["genotype"])).to_owned())
            .unwrap_or_default()
    };
    assert_eq!(genotype("out", SNV), "hom", "outside the deletion");
    assert_eq!(
        genotype("edge", ("chr1", DEL.1, "G", "A")),
        "hom",
        "AnnotSV's SV_start is the base before the deleted bases"
    );
    assert_eq!(
        genotype("first", ("chr1", DEL.1 + 1, "G", "A")),
        "hemizygous",
        "the first deleted base"
    );
}

/// test_cross_engine_biallelic.py checks 6 and 7.
#[test]
fn smn1_and_hba_counts_pair_with_a_hit_and_no_caller_is_named() {
    let t = Dir::new("two-cn");
    let smn = ("chr5", 70951950, "G", "A");
    let d = sample(&t, "smn1", &[hit(smn, "SMN1", "0/1", None)], None);
    step35(&t, &d, "1", "", None);
    assert!(field(&two(&d, "SMN1"), "deletion").contains("one SMN1 copy"));
    let d = sample(&t, "smn2", &[hit(smn, "SMN1", "0/1", None)], None);
    step35(&t, &d, "2", "", None);
    assert!(two(&d, "SMN1").is_empty());
    let hba2 = ("chr16", 173600, "C", "T");
    for (n, want) in [("3", true), ("4", false), ("0", false)] {
        let d = sample(
            &t,
            &format!("hba{n}"),
            &[hit(hba2, "HBA2", "0/1", None)],
            None,
        );
        step35(&t, &d, "2", "", Some(n));
        let tw = two(&d, "HBA2");
        assert_eq!(!tw.is_empty(), want, "{n} alpha-globin genes");
        if want {
            assert_eq!(field(&tw, "verdict"), "possibly biallelic");
            assert!(field(&tw, "deletion").contains("3 of 4 alpha-globin genes"));
        }
    }

    let d = sample(&t, "none", &[hit(SNV, "MUTYH", "0/1", None)], None);
    let c = carrier_of(&d);
    let th = at(&c, &["two_hits"]);
    assert_eq!(at(th, &["genes"]), &Value::Array(Vec::new()));
    let not_read: Vec<&str> = list(at(th, &["deletions_not_read"]))
        .iter()
        .map(s)
        .collect();
    assert!(
        not_read.iter().any(|x| x.contains("AnnotSV")),
        "{not_read:?}"
    );
    assert!(
        not_read.iter().any(|x| x.contains("step 35")),
        "{not_read:?}"
    );
    let lines = carrier::lines(&D::new(&c, "carrier")).unwrap();
    assert!(
        lines.iter().any(|x| x.contains("not read: AnnotSV")),
        "{lines:?}"
    );
}

/// test_cross_engine_biallelic.py check 8.
#[test]
fn both_reports_print_the_two_hit_line_and_the_hemizygous_genotype() {
    let t = Dir::new("two-report");
    let d = sample(
        &t,
        "r",
        &[
            hit(INSIDE, "MUTYH", "1/1", None),
            h("chr1", 45330500, "A", "G", "MUTYH", "0/1"),
        ],
        Some(&del("0/1", None)),
    );
    t.put_gz("r/S/clinical/S_clinical.vcf.gz", VCF_HEAD);
    let (txt, html) = reports(&collect(&d));
    assert!(txt.contains("Biallelic: MUTYH"));
    assert!(html.contains("Biallelic: MUTYH"));
    // The hit has no unchanged copy: the summary lists it apart from the
    // carriers (hom_or_hemizygous), and both reports still print it in the
    // carrier block, with the hemizygous note, and not as a carrier.
    assert!(html.contains("hemizygous (called homozygous inside DEL"));
    assert!(txt.contains("hemizygous (called homozygous inside DEL"));
    for r in [&txt, &html] {
        assert!(r.contains("1 finding(s), 1 with no unchanged copy"), "{r}");
        assert!(
            r.contains("the only copy changed, not the pattern of a carrier"),
            "{r}"
        );
    }
    assert!(
        html.replace("&gt;", ">")
            .contains("chr1:45336000 G>A hemizygous")
    );
    assert!(txt.contains("chr1:45336000 G>A hemizygous"));
    assert!(txt.contains("biallelic with DEL chr1:45335000-45337000"));
    assert!(html.contains("biallelic with DEL chr1:45335000-45337000"));
}

/// A homozygous P/LP row has no unchanged copy: the carrier block lists it
/// apart from the carriers and does not count it as one.
#[test]
fn a_homozygous_row_is_listed_apart_from_the_carriers() {
    let t = Dir::new("hom-apart");
    let d = sample(&t, "hom", &[hit(SNV, "MUTYH", "1/1", None)], None);
    let c = carrier_of(&d);
    let lines = carrier::lines(&D::new(&c, "carrier")).unwrap();
    assert!(
        lines[0].contains("0 finding(s), 1 with no unchanged copy"),
        "{lines:?}"
    );
    let line: Vec<&String> = lines.iter().filter(|x| x.starts_with("MUTYH (")).collect();
    assert_eq!(line.len(), 1, "{lines:?}");
    assert!(line[0].contains(" hom, ClinVar "), "{}", line[0]);
    assert!(
        line[0].ends_with(": two changed copies, not the pattern of a carrier"),
        "{}",
        line[0]
    );
}

/// test_cross_engine_biallelic.py check 10.
#[test]
fn two_hit_edges() {
    let t = Dir::new("two-edges");
    // MUTYH ends at chr1:45340893: a deletion whose padding base is that one
    // removes none of it.
    let d = sample(
        &t,
        "pad",
        &[hit(SNV, "MUTYH", "0/1", None)],
        Some(&tab15(&sv(
            ("chr1", 45340893, 45342000),
            "MUTYH",
            "4",
            "yes",
            "0/1",
            None,
        ))),
    );
    assert!(two(&d, "MUTYH").is_empty());
    let d = sample(
        &t,
        "cis-inside",
        &[hit(INSIDE, "MUTYH", "1|0", Some(45330000))],
        Some(&del("1|0", Some(45330000))),
    );
    assert!(
        two(&d, "MUTYH").is_empty(),
        "a het hit inside a deletion, in cis"
    );

    // A HIGH-impact record (no ClinVar hit) phased in trans with the deletion.
    let csq = format!(
        "##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. \
         Format: Allele|Consequence|IMPACT|SYMBOL|MAX_AF\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
         chr1\t{}\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|MUTYH|0.0001\tGT:GQ:DP:AD:PS\t0|1:50:30:15,15:45330000\n",
        SNV.1
    );
    let d = sample(&t, "high", &[], Some(&del("1|0", Some(45330000))));
    t.put_gz("high/S/clinical/S_clinical.vcf.gz", &csq);
    let sfv = sf(&t, &d);
    let rep = genes_of(&sfv, "high_impact", "MUTYH");
    assert_eq!(rep.len(), 1);
    assert!(s(at(rep[0], &["two_hit"])).starts_with("biallelic with DEL"));

    // A cis and a trans hit in one gene: only the trans one carries the note.
    let d = sample(
        &t,
        "cis-trans",
        &[
            hit(("chr1", 45330500, "A", "G"), "MUTYH", "1|0", Some(45330000)),
            hit(SNV, "MUTYH", "0|1", Some(45330000)),
        ],
        Some(&del("1|0", Some(45330000))),
    );
    let sfv = sf(&t, &d);
    let note = |variant: &str| {
        reported(&sfv, "MUTYH")
            .iter()
            .find(|x| s(at(x, &["variant"])) == variant)
            .map(|x| at(x, &["two_hit"]).clone())
    };
    assert!(s(&note("chr1:45332000 C>T").unwrap()).starts_with("biallelic with DEL"));
    assert_eq!(note("chr1:45330500 A>G"), Some(Value::Null));

    // An unreadable carrier list: the SF card says what was not read.
    let d = sample(
        &t,
        "no-list",
        &[hit(SNV, "MUTYH", "0/1", None)],
        Some(&del("0/1", None)),
    );
    let broken = t.path("assets-broken");
    copy_tree(
        &PathBuf::from(format!("{REPO}/assets")),
        &PathBuf::from(&broken),
    );
    let list_path = format!("{broken}/carrier_genes.tsv");
    let text = fs::read_to_string(&list_path).unwrap();
    fs::write(
        &list_path,
        text.replacen("gene\tchrom\tstart", "name\tchrom\tstart", 1),
    )
    .unwrap();
    t.put_gz("no-list/S/clinical/S_clinical.vcf.gz", VCF_HEAD);
    let summ = collect_with(&d, &broken);
    let sfv = at(&summ, &["sections", "clinical", "data", "acmg_sf"]);
    assert!(s(at(sfv, &["two_hit_note"])).contains("carrier_genes.tsv"));
    let lines = acmg_sf::apart_lines(&D::new(sfv, "acmg_sf")).unwrap();
    assert!(
        lines.iter().any(|x| x.contains("carrier_genes.tsv")),
        "{lines:?}"
    );
}

/// `shutil.copytree(from, to)`.
fn copy_tree(from: &std::path::Path, to: &std::path::Path) {
    fs::create_dir_all(to).unwrap();
    for e in fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let dest = to.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_tree(&e.path(), &dest);
        } else {
            fs::copy(e.path(), dest).unwrap();
        }
    }
}

// --- callability -------------------------------------------------------------------

const REASONS: [&str; 5] = [
    "CMRG",
    "segdup",
    "false_duplication",
    "collapsed",
    "ALT-only",
];

/// test_callability.py check 1.
#[test]
fn the_callability_catalog() {
    let text = fs::read_to_string(format!("{REPO}/assets/callability_genes.tsv")).unwrap();
    let header = text.lines().find(|l| !l.starts_with('#')).unwrap();
    assert_eq!(header, "gene\thgnc_id\tchrom\tstart\tend\treasons\towner");
    let cat = asset("callability_genes.tsv");
    let mut genes: Vec<&str> = cat.iter().map(|r| r["gene"].as_str()).collect();
    let n = genes.len();
    genes.sort_unstable();
    genes.dedup();
    assert!(genes.len() == n && n > 1000, "one row per gene: {n}");
    assert!(cat.iter().all(|r| r["hgnc_id"].starts_with("HGNC:")));
    assert!(
        cat.iter()
            .all(|r| r["reasons"].split(',').all(|x| REASONS.contains(&x)))
    );
    let main: Vec<String> = (1..=22)
        .map(|c| format!("chr{c}"))
        .chain(["chrX".into(), "chrY".into()])
        .collect();
    let bad: Vec<&str> = cat
        .iter()
        .filter(|r| {
            let alt = r["reasons"].contains("ALT-only");
            let none = r["start"] == ".";
            alt != none
                || (none && !r["chrom"].ends_with("_alt"))
                || (!none
                    && (!main.contains(&r["chrom"])
                        || r["start"].parse::<i64>().unwrap() >= r["end"].parse::<i64>().unwrap()))
        })
        .map(|r| r["gene"].as_str())
        .collect();
    assert!(bad.is_empty(), "{bad:?}");
    assert_eq!(
        cat.iter().filter(|r| r["reasons"].contains("CMRG")).count(),
        273
    );
    let by = |g: &str| {
        cat.iter()
            .find(|r| r["gene"] == g)
            .cloned()
            .unwrap_or_default()
    };
    for g in ["GSTT1", "HLA-DRB3"] {
        assert_eq!(
            (by(g)["reasons"].as_str(), by(g)["owner"].as_str()),
            ("ALT-only", ".")
        );
    }
    for g in ["KCNE1", "CBS", "CRYAA"] {
        assert!(by(g)["reasons"].contains("false_duplication"), "{g}");
    }
    let parascopy: Vec<String> =
        table(&fs::read_to_string(format!("{REPO}/assets/paralog_loci.tsv")).unwrap())
            .iter()
            .flat_map(|r| r["genes"].split(',').map(str::to_owned).collect::<Vec<_>>())
            .collect();
    let pypgx: Vec<String> = asset("pgx_genes.tsv")
        .iter()
        .filter(|r| r["mode"] == "bam" && r["callable_without_alt"] != "no")
        .map(|r| r["gene"].clone())
        .collect();
    let wrong: Vec<&str> = cat
        .iter()
        .filter(|r| {
            let g = &r["gene"];
            let o = r["owner"].as_str();
            (parascopy.contains(g) && o != "parascopy")
                || (!parascopy.contains(g) && pypgx.contains(g) && o != "pypgx")
                || (![".", "parascopy", "pypgx"].contains(&o) && !o.starts_with("cn_loci:"))
        })
        .map(|r| r["gene"].as_str())
        .collect();
    assert!(wrong.is_empty(), "owners: {wrong:?}");
    assert_eq!(by("HBA1")["owner"], "cn_loci:HBA");
    assert_eq!(by("HBA2")["owner"], "cn_loci:HBA");
    assert_eq!(
        by("RSRP1")["owner"],
        ".",
        "RSRP1 only holds RHD's depth target"
    );
    assert_eq!(by("RHD")["owner"], "cn_loci:RHD");
}

/// mosdepth's thresholds table over the catalog's genes: every base at 10x,
/// except the genes in `zero` (none) and in `drop` (no row).
fn thresholds(zero: &[&str], drop: &dyn Fn(&str, &str) -> bool) -> String {
    let mut out = "#chrom\tstart\tend\tregion\t10X\n".to_owned();
    for r in asset("callability_genes.tsv") {
        if r["start"] == "." || drop(&r["gene"], &r["chrom"]) {
            continue;
        }
        let n = r["end"].parse::<i64>().unwrap() - r["start"].parse::<i64>().unwrap();
        let at10 = if zero.contains(&r["gene"].as_str()) {
            0
        } else {
            n
        };
        out += &format!(
            "{}\t{}\t{}\t{}\t{at10}\n",
            r["chrom"], r["start"], r["end"], r["gene"]
        );
    }
    out
}

/// The callability section og summary collect reads from a table at `rel`.
fn callability(t: &Dir, name: &str, rel: &str, table: &str) -> Value {
    t.put_gz(&format!("{name}/S/{rel}"), table);
    let summ = collect(&t.path(&format!("{name}/S")));
    at(&summ, &["sections", "callability"]).clone()
}

/// `{gene: status}` of the genes not read.
fn readings(d: &Value) -> Vec<(String, String)> {
    list(at(d, &["not_read"]))
        .iter()
        .chain(list(at(d, &["not_read_owned"])))
        .map(|g| {
            (
                s(at(g, &["gene"])).to_owned(),
                s(at(g, &["status"])).to_owned(),
            )
        })
        .collect()
}

const MAPQ0: &str = "KCNE1"; // every read MAPQ 0, as GRCh38's false duplication makes them
const UNIQUE: &str = "CBS"; // every base at 10x
const OWNED: &str = "HBA1"; // not read, but step 35's HBA count calls it

/// What is wrong with check 2 on data `d`: the MAPQ 0 gene and the
/// ALT-only genes not read, the unique gene read.
fn mapq0_and_alt(d: &Value) -> Vec<String> {
    let got = readings(d);
    let status = |g: &str| got.iter().find(|(x, _)| x == g).map(|(_, s)| s.as_str());
    let mut out = Vec::new();
    if status(MAPQ0) != Some("not read") {
        out.push(format!("{MAPQ0}: {:?}", status(MAPQ0)));
    }
    if let Some(st) = status(UNIQUE) {
        out.push(format!("{UNIQUE}: {st}"));
    }
    for g in ["GSTT1", "HLA-DRB3"] {
        if status(g) != Some("not on this reference") {
            out.push(format!("{g}: {:?}", status(g)));
        }
    }
    out
}

/// Check 3: the unowned genes not read are the catalog query's.
fn count_matches(d: &Value) -> bool {
    let query: Vec<String> = asset("callability_genes.tsv")
        .iter()
        .filter(|r| r["owner"] == "." && (r["start"] == "." || r["gene"] == MAPQ0))
        .map(|r| r["gene"].clone())
        .collect();
    let got: Vec<String> = list(at(d, &["not_read"]))
        .iter()
        .map(|g| s(at(g, &["gene"])).to_owned())
        .collect();
    got == query
}

/// The planted bug of the controls: an ALT-only gene taken as read.
fn alt_read_bug(d: &Value) -> Value {
    let mut d = d.clone();
    let kept: Vec<Value> = list(at(&d, &["not_read"]))
        .iter()
        .filter(|g| s(at(g, &["status"])) != "not on this reference")
        .cloned()
        .collect();
    let moved = list(at(&d, &["not_on_reference"])).len() as i128;
    if let Value::Int(n) = at(&d, &["read"]) {
        let read = Value::Int(n + moved);
        set(&mut d, &["read"], read);
    }
    set(&mut d, &["not_read"], Value::Array(kept));
    d
}

/// test_callability.py checks 2 to 4.
#[test]
fn step_16b_reads_the_genes_not_read_and_both_reports_list_them() {
    let t = Dir::new("callability");
    let table = thresholds(&[MAPQ0, OWNED], &|_, _| false);
    let sec = callability(&t, "mosdepth", "mosdepth/S_callability.bed.gz", &table);
    assert_eq!(s(at(&sec, &["state"])), "ok", "{sec:?}");
    let d = at(&sec, &["data"]).clone();
    assert!(
        !mapq0_and_alt(&alt_read_bug(&d)).is_empty(),
        "control: a reader that takes an ALT-only gene as read fails check 2"
    );
    let problems = mapq0_and_alt(&d);
    assert!(problems.is_empty(), "{problems:?}");
    assert!(
        !count_matches(&alt_read_bug(&d)),
        "control: the planted bug's list differs from the catalog query"
    );
    assert!(
        count_matches(&d),
        "the genes not read with no owner are the catalog query's"
    );
    let owned = list(at(&d, &["not_read_owned"]));
    assert!(
        owned
            .iter()
            .any(|g| s(at(g, &["gene"])) == OWNED && s(at(g, &["owner"])) == "cn_loci:HBA")
    );
    assert!(owned.iter().all(|g| s(at(g, &["gene"])) != "GSTT1"));
    assert!(
        list(at(&d, &["not_read"]))
            .iter()
            .any(|g| s(at(g, &["gene"])) == "GSTT1")
    );

    // The section and both reports, from the Nextflow layout (coverage/).
    t.put_gz("nf/S/coverage/S_callability.bed.gz", &table);
    let summ = collect(&t.path("nf/S"));
    assert_eq!(s(at(&summ, &["sections", "callability", "state"])), "ok");
    let n = list(at(&d, &["not_read"])).len();
    let (txt, html) = reports(&summ);
    let line = format!("Not read, and no other step calls them ({n}): ");
    assert!(txt.contains(&line), "text: {line}");
    assert!(txt.contains(&format!("{MAPQ0} (0%)")));
    assert!(txt.contains("HLA-DRB3 (not on this reference)"));
    assert!(html.contains(&*og_render::escape::e(&line)), "HTML: {line}");
    assert!(html.contains("Genes This Genome Could Not Read Reliably"));

    // A table from another catalog, in mosdepth/ beside coverage/: unreadable,
    // not a reading.
    let sec = callability(
        &t,
        "nf",
        "mosdepth/S_callability.bed.gz",
        &thresholds(&[], &|g, _| g == UNIQUE),
    );
    assert_eq!(s(at(&sec, &["state"])), "unreadable");
    assert!(s(at(&sec, &["note"])).contains(UNIQUE));

    // A reference without chr3 (as the e2e fixture's): its genes are not on
    // this reference, not unreadable.
    let sec = callability(
        &t,
        "nf",
        "mosdepth/S_callability.bed.gz",
        &thresholds(&[], &|_, c| c == "chr3"),
    );
    assert_eq!(s(at(&sec, &["state"])), "ok");
    let got = readings(at(&sec, &["data"]));
    let chr3: Vec<String> = asset("callability_genes.tsv")
        .iter()
        .filter(|r| r["chrom"] == "chr3")
        .map(|r| r["gene"].clone())
        .collect();
    assert!(!chr3.is_empty());
    for g in &chr3 {
        assert!(
            got.iter()
                .any(|(x, st)| x == g && st == "not on this reference"),
            "{g}"
        );
    }
}
