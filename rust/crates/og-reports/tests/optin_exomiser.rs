//! `og exomiser prepare` and `og exomiser report` (step 48): the
//! symptom-driven ranking. Moved from tests/test_exomiser_report.py (its
//! cases 1 to 4 and the stub table's columns); the summary's reading is in
//! og-summary's `tests/optin_summary.rs` and the report lines in og-render's
//! `tests/optin_report.rs`.
//!
//! Every input is synthetic: a few-term OBO file, one- and several-sample
//! VCF headers, a pedigree, and Exomiser variants TSVs in the column layout
//! of Exomiser 15's TsvVariantResultsWriter.

use pycompat::json::Value;
use std::collections::HashMap;
use std::fs;
use std::io::Write as _;
use std::num::NonZero;
use std::path::PathBuf;

const MODULE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../../modules/local/exomiser/main.nf"
);

const OBO: &str = "format-version: 1.2

[Term]
id: HP:0001156
name: Brachydactyly
alt_id: HP:0009803

[Term]
id: HP:0001363
name: Craniosynostosis

[Term]
id: HP:0011304
name: Broad thumb

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0001250
name: obsolete Seizure
is_obsolete: true
replaced_by: HP:0012345

[Typedef]
id: part_of
name: part of
";

const HEADER: [&str; 41] = [
    "#RANK",
    "ID",
    "GENE_SYMBOL",
    "ENTREZ_GENE_ID",
    "MOI",
    "P-VALUE",
    "EXOMISER_GENE_COMBINED_SCORE",
    "EXOMISER_GENE_PHENO_SCORE",
    "EXOMISER_GENE_VARIANT_SCORE",
    "EXOMISER_VARIANT_SCORE",
    "CONTRIBUTING_VARIANT",
    "WHITELIST_VARIANT",
    "VCF_ID",
    "RS_ID",
    "CONTIG",
    "START",
    "END",
    "REF",
    "ALT",
    "CHANGE_LENGTH",
    "QUAL",
    "FILTER",
    "GENOTYPE",
    "FUNCTIONAL_CLASS",
    "HGVS",
    "EXOMISER_ACMG_CLASSIFICATION",
    "EXOMISER_ACMG_EVIDENCE",
    "EXOMISER_ACMG_DISEASE_ID",
    "EXOMISER_ACMG_DISEASE_NAME",
    "CLINVAR_VARIATION_ID",
    "CLINVAR_PRIMARY_INTERPRETATION",
    "CLINVAR_STAR_RATING",
    "GENE_CONSTRAINT_LOEUF",
    "GENE_CONSTRAINT_LOEUF_LOWER",
    "GENE_CONSTRAINT_LOEUF_UPPER",
    "MAX_FREQ_SOURCE",
    "MAX_FREQ",
    "ALL_FREQ",
    "MAX_PATH_SOURCE",
    "MAX_PATH",
    "ALL_PATH",
];

type Variant = HashMap<&'static str, String>;

/// One variants TSV row; `extra` sets further columns.
#[allow(clippy::too_many_arguments)]
fn variant(
    rank: u32,
    gene: &str,
    moi: &str,
    contig: &str,
    pos: u64,
    refa: &str,
    alt: &str,
    gt: &str,
    effect: &str,
    hgvs: &str,
    acmg: &str,
    extra: &[(&'static str, &str)],
) -> Variant {
    let mut v: Variant = HEADER.iter().map(|k| (*k, String::new())).collect();
    for (k, x) in [
        ("#RANK", rank.to_string()),
        ("ID", format!("{contig}-{pos}-{refa}-{alt}_{moi}")),
        ("GENE_SYMBOL", gene.into()),
        ("ENTREZ_GENE_ID", "1".into()),
        ("MOI", moi.into()),
        ("P-VALUE", "0.0010".into()),
        ("EXOMISER_GENE_COMBINED_SCORE", "0.9000".into()),
        ("EXOMISER_GENE_PHENO_SCORE", "0.8000".into()),
        ("EXOMISER_GENE_VARIANT_SCORE", "0.9500".into()),
        ("EXOMISER_VARIANT_SCORE", "0.9500".into()),
        ("CONTRIBUTING_VARIANT", "1".into()),
        ("WHITELIST_VARIANT", "0".into()),
        ("CONTIG", contig.into()),
        ("START", pos.to_string()),
        ("END", pos.to_string()),
        ("REF", refa.into()),
        ("ALT", alt.into()),
        ("CHANGE_LENGTH", "0".into()),
        ("QUAL", "50.0000".into()),
        ("FILTER", "PASS".into()),
        ("GENOTYPE", gt.into()),
        ("FUNCTIONAL_CLASS", effect.into()),
        ("HGVS", hgvs.into()),
        ("EXOMISER_ACMG_CLASSIFICATION", acmg.into()),
        ("CLINVAR_PRIMARY_INTERPRETATION", "NOT_PROVIDED".into()),
        ("CLINVAR_STAR_RATING", "0".into()),
    ] {
        v.insert(k, x);
    }
    for (k, x) in extra {
        v.insert(k, (*x).to_owned());
    }
    v
}

fn rows() -> Vec<Variant> {
    vec![
        variant(
            1,
            "GENEA",
            "AD",
            "chr10",
            123_256_215,
            "T",
            "G",
            "0/1",
            "missense_variant",
            "GENEA:ENST1:c.1688A>C:p.(Glu563Ala)",
            "LIKELY_PATHOGENIC",
            &[
                ("EXOMISER_ACMG_EVIDENCE", "PM2,PP3_Strong,PP4"),
                ("EXOMISER_ACMG_DISEASE_ID", "OMIM:100001"),
                ("EXOMISER_ACMG_DISEASE_NAME", "Synthetic syndrome one"),
                ("CLINVAR_PRIMARY_INTERPRETATION", "LIKELY_PATHOGENIC"),
                ("CLINVAR_STAR_RATING", "1"),
                ("EXOMISER_GENE_COMBINED_SCORE", "0.9981"),
            ],
        ),
        variant(
            2,
            "GENEB",
            "AR",
            "chr6",
            1000,
            "G",
            "A",
            "0/1",
            "splice_donor_variant",
            "GENEB:ENST2:c.10+1G>A:p.?",
            "UNCERTAIN_SIGNIFICANCE",
            &[
                ("EXOMISER_ACMG_EVIDENCE", "PVS1_Strong"),
                ("MAX_FREQ_SOURCE", "GNOMAD_E_SAS"),
                ("MAX_FREQ", "0.0032"),
            ],
        ),
        variant(
            2,
            "GENEB",
            "AR",
            "chr6",
            2000,
            "C",
            "T",
            "1|0",
            "stop_gained",
            "GENEB:ENST2:c.300C>T:p.(Gln100*)",
            "PATHOGENIC",
            &[("EXOMISER_ACMG_EVIDENCE", "PVS1,PM2")],
        ),
        variant(
            2,
            "GENEB",
            "AR",
            "chr6",
            3000,
            "A",
            "C",
            "0/1",
            "missense_variant",
            "GENEB:ENST2:c.400A>C:p.(Lys134Thr)",
            "NOT_AVAILABLE",
            &[("CONTRIBUTING_VARIANT", "0")],
        ),
        variant(
            3,
            "GENEC",
            "XR",
            "chrX",
            5000,
            "A",
            "G",
            "1",
            "missense_variant",
            "GENEC:ENST3:c.5A>G:p.(Asp2Gly)",
            "UNCERTAIN_SIGNIFICANCE",
            &[],
        ),
        variant(
            4,
            "GENED",
            "AD",
            "chr1",
            7000,
            "A",
            "T",
            "0/1",
            "missense_variant",
            "GENED:ENST4:c.7A>T:p.(Asn3Ile)",
            "VUS",
            &[],
        ),
    ]
}

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d =
            std::env::temp_dir().join(format!("og-optin-exomiser-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn p(&self, n: &str) -> String {
        self.0.join(n).to_str().unwrap().to_owned()
    }
    fn put(&self, n: &str, t: &str) -> String {
        fs::write(self.p(n), t).unwrap();
        self.p(n)
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn og(verb: &str, args: &[String]) -> (u8, String) {
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = og_reports::exomiser::cli(Some(verb), args, &mut o, &mut e);
    (code, String::from_utf8(e).unwrap())
}

fn strings(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| (*s).to_owned()).collect()
}

/// A bgzipped VCF with one record and these samples.
fn vcf(d: &Dir, samples: &[&str]) -> String {
    let mut text =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT".to_owned();
    for s in samples {
        text += &format!("\t{s}");
    }
    text += "\nchr1\t100\t.\tA\tG\t50\tPASS\t.\tGT";
    for _ in samples {
        text += "\t0/1";
    }
    text += "\n";
    let mut w = og_io::bgzf::Writer::new(Vec::new(), NonZero::new(1).unwrap()).unwrap();
    w.write_all(text.as_bytes()).unwrap();
    let p = d.p("in.vcf.gz");
    fs::write(&p, w.finish().unwrap()).unwrap();
    p
}

fn prepare(d: &Dir, hpo: &str, samples: &[&str], extra: &[&str], obo: bool) -> (u8, String) {
    let v = vcf(d, samples);
    let mut args = strings(&[
        "--hpo",
        hpo,
        "--sample",
        "S1",
        "--vcf",
        &v,
        "--assembly",
        "hg38",
        "--out-dir",
        "results",
        "--analysis",
        &d.p("a.yml"),
        "--terms",
        &d.p("terms.tsv"),
    ]);
    if obo {
        args.extend(strings(&["--obo", &d.put("hp.obo", OBO)]));
    }
    args.extend(strings(extra));
    og("prepare", &args)
}

fn read_rows(path: &str) -> Vec<HashMap<String, String>> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.trim_end_matches('\n').split('\n');
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or(&Value::Null),
        _ => &Value::Null,
    })
}

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a
            .iter()
            .map(|x| match x {
                Value::Str(s) => s.clone(),
                Value::Object(m) => m.keys().next().cloned().unwrap_or_default(),
                _ => String::new(),
            })
            .collect(),
        _ => Vec::new(),
    }
}

fn s(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => "",
    }
}

#[test]
fn prepare_checks_the_terms_and_writes_the_analysis() {
    let d = Dir::new("prepare");
    let (code, err) = prepare(
        &d,
        "HP:0001156, HP:0001363,HP:0009803 HP:0011304",
        &["S1"],
        &[],
        true,
    );
    assert_eq!(code, 0, "{err}");
    let a = pycompat::json::loads(&fs::read_to_string(d.p("a.yml")).unwrap()).unwrap();
    let an = at(&a, &["analysis"]);
    assert_eq!(
        texts(at(an, &["hpoIds"])),
        ["HP:0001156", "HP:0001363", "HP:0011304"],
        "the alternative id is its term and the repeat is dropped"
    );
    assert!(s(at(an, &["vcf"])).ends_with("in.vcf.gz") && s(at(an, &["genomeAssembly"])) == "hg38");
    let steps = texts(at(an, &["steps"]));
    assert_eq!(
        steps[steps.len() - 2..],
        ["omimPrioritiser", "hiPhivePrioritiser"]
    );
    assert!(
        steps.iter().any(|x| x == "variantEffectFilter")
            && steps.iter().any(|x| x == "frequencyFilter")
    );
    assert!(matches!(at(an, &["ped"]), Value::Null) && matches!(at(an, &["proband"]), Value::Null));
    let oo = at(&a, &["outputOptions"]);
    let mut formats = texts(at(oo, &["outputFormats"]));
    formats.sort();
    assert_eq!(
        (
            s(at(oo, &["outputFileName"])),
            s(at(oo, &["outputDirectory"])),
            formats
        ),
        (
            "S1_exomiser",
            "results",
            strings(&["HTML", "TSV_GENE", "TSV_VARIANT"])
        )
    );
    let terms: Vec<(String, String)> = read_rows(&d.p("terms.tsv"))
        .into_iter()
        .map(|t| (t["hpo_id"].clone(), t["label"].clone()))
        .collect();
    assert_eq!(
        terms,
        [
            ("HP:0001156".into(), "Brachydactyly".into()),
            ("HP:0001363".into(), "Craniosynostosis".into()),
            ("HP:0011304".into(), "Broad thumb".into())
        ]
    );
    let (code, err) = prepare(&d, "HP:0001156", &["S1"], &[], false);
    assert_eq!(code, 0, "{err}");
    assert_eq!(
        read_rows(&d.p("terms.tsv"))[0]["label"],
        "",
        "no ontology: no label"
    );
}

#[test]
fn prepare_stops() {
    let d = Dir::new("stops");
    for (hpo, want) in [
        ("HP:12345", "not an HPO id"),
        ("Seizure", "not an HPO id"),
        ("HP:0009999", "not a term of the HPO release"),
        ("HP:0001250", "obsolete"),
        (" , ", "no HPO term given"),
    ] {
        let (code, err) = prepare(&d, hpo, &["S1"], &[], true);
        assert!(code == 1 && err.contains(want), "{hpo}: {code} {err}");
    }
    assert!(
        prepare(&d, "HP:0001250", &["S1"], &[], true)
            .1
            .contains("use HP:0012345")
    );
    let (code, err) = prepare(&d, "HP:0001156", &["S1", "S2"], &[], true);
    assert!(code == 1 && err.contains("has 2 samples"), "{err}");
    let (code, err) = prepare(&d, "HP:0001156", &["S1"], &["--affected", "S2"], true);
    assert!(code == 1 && err.contains("need --pedigree"), "{err}");
}

const PEDIGREE: &str = "sample\tsex\tfather\tmother\tpartner
S1\tfemale\tDAD\tMUM\t.
SIB\tmale\tDAD\tMUM\t.
DAD\tmale\t.\t.\tMUM
MUM\tfemale\t.\t.\tDAD
";

#[test]
fn the_family_file() {
    let d = Dir::new("family");
    let ped = d.put("pedigree.tsv", PEDIGREE);
    let out = d.p("fam.ped");
    let fam = |extra: &[&str]| {
        let mut a = vec!["--pedigree", ped.as_str(), "--ped", out.as_str()];
        a.extend(extra);
        a.iter().map(|x| (*x).to_owned()).collect::<Vec<String>>()
    };
    let three = ["S1", "SIB", "DAD"];
    let lines = || -> HashMap<String, Vec<String>> {
        fs::read_to_string(&out)
            .unwrap()
            .lines()
            .map(|l| {
                let x: Vec<String> = l.split('\t').map(str::to_owned).collect();
                (x[1].clone(), x)
            })
            .collect()
    };
    let e = fam(&["--unaffected", "DAD", "--affected", "SIB"]);
    let e: Vec<&str> = e.iter().map(String::as_str).collect();
    let (code, err) = prepare(&d, "HP:0001156", &three, &e, true);
    assert_eq!(code, 0, "{err}");
    let got = lines();
    let mut names: Vec<&String> = got.keys().collect();
    names.sort();
    assert_eq!(names, ["DAD", "S1", "SIB"]);
    assert_eq!(got["S1"], ["family", "S1", "DAD", "0", "2", "2"]);
    assert!(got["SIB"][5] == "2" && got["SIB"][4] == "1");
    assert_eq!(got["DAD"][5], "1");
    let an = pycompat::json::loads(&fs::read_to_string(d.p("a.yml")).unwrap()).unwrap();
    assert_eq!(
        (
            s(at(&an, &["analysis", "ped"])),
            s(at(&an, &["analysis", "proband"]))
        ),
        (out.as_str(), "S1")
    );
    let e = fam(&[]);
    let e: Vec<&str> = e.iter().map(String::as_str).collect();
    prepare(&d, "HP:0001156", &three, &e, true);
    let got = lines();
    assert!(
        got["SIB"][5] == "0" && got["DAD"][5] == "0",
        "unnamed members are unknown"
    );
    for (extra, want) in [
        (&["--affected", "MUM"][..], "not a member of the family VCF"),
        (&["--unaffected", "S1"][..], "cannot be unaffected"),
        (
            &["--affected", "SIB", "--unaffected", "SIB"][..],
            "both affected and unaffected",
        ),
    ] {
        let e = fam(extra);
        let e: Vec<&str> = e.iter().map(String::as_str).collect();
        let (code, err) = prepare(&d, "HP:0001156", &three, &e, true);
        assert!(code == 1 && err.contains(want), "{want}: {code} {err}");
    }
    let e = fam(&[]);
    let e: Vec<&str> = e.iter().map(String::as_str).collect();
    let (code, err) = prepare(&d, "HP:0001156", &["S1", "STRANGER"], &e, true);
    assert!(code == 1 && err.contains("STRANGER"), "{err}");
}

fn tsv(variants: &[Variant], header: &[&str]) -> String {
    let mut t = header.join("\t") + "\n";
    for v in variants {
        t += &header
            .iter()
            .map(|k| v[k].as_str())
            .collect::<Vec<_>>()
            .join("\t");
        t += "\n";
    }
    t
}

fn report(d: &Dir, raw: &str) -> (u8, String) {
    let v = d.put("v.tsv", raw);
    let t = d.put(
        "terms.tsv",
        "hpo_id\tlabel\tgiven_as\nHP:0001156\tBrachydactyly\tHP:0001156\n",
    );
    og(
        "report",
        &strings(&[
            "--sample",
            "S1",
            "--variants",
            &v,
            "--terms",
            &t,
            "--out",
            &d.p("top.tsv"),
            "--text",
            &d.p("words.txt"),
            "--top",
            "3",
        ]),
    )
}

#[test]
fn the_report_keeps_the_top_ranks_and_their_contributing_variants() {
    let d = Dir::new("report");
    let (code, err) = report(&d, &tsv(&rows(), &HEADER));
    assert_eq!(code, 0, "{err}");
    let top = read_rows(&d.p("top.tsv"));
    let got: Vec<[&str; 3]> = top
        .iter()
        .map(|x| [x["rank"].as_str(), &x["gene"], &x["variant"]])
        .collect();
    assert_eq!(
        got,
        [
            ["1", "GENEA", "chr10:123256215 T>G"],
            ["2", "GENEB", "chr6:1000 G>A"],
            ["2", "GENEB", "chr6:2000 C>T"],
            ["3", "GENEC", "chrX:5000 A>G"]
        ]
    );
    let a = &top[0];
    assert_eq!(
        [
            "exomiser_score",
            "phenotype_score",
            "variant_score",
            "moi",
            "acmg",
            "acmg_evidence"
        ]
        .map(|k| a[k].as_str()),
        [
            "0.998",
            "0.800",
            "0.950",
            "AD",
            "LIKELY_PATHOGENIC",
            "PM2,PP3_Strong,PP4"
        ]
    );
    assert_eq!(
        (a["disease"].as_str(), a["clinvar"].as_str()),
        (
            "Synthetic syndrome one (OMIM:100001)",
            "LIKELY_PATHOGENIC (1 star)"
        )
    );
    assert_eq!(
        (
            top[1]["max_frequency"].as_str(),
            a["max_frequency"].as_str()
        ),
        ("GNOMAD_E_SAS 0.0032%", "not seen")
    );
    let words = fs::read_to_string(d.p("words.txt")).unwrap();
    for want in [
        "Symptoms searched: Brachydactyly (HP:0001156).",
        "1. GENEA (autosomal dominant): Exomiser score 0.998 (symptom match 0.800, variant 0.950).",
        "chr10:123256215 T>G heterozygous, missense_variant",
        "ACMG likely pathogenic (PM2,PP3_Strong,PP4)",
        "chr6:1000 G>A heterozygous",
        "chr6:2000 C>T heterozygous",
        "chrX:5000 A>G hemizygous",
        "X-linked recessive",
        "not a diagnosis",
    ] {
        assert!(words.contains(want), "{want} not in {words}");
    }
    assert!(!words.contains("GENED"), "rank 4 is past the top 3");
}

#[test]
fn a_family_genotype_names_each_member() {
    let d = Dir::new("famreport");
    d.put(
        "S1_exomiser.ped",
        "family\tM\t0\t0\t2\t1\nfamily\tF\t0\t0\t1\t1\nfamily\tS1\tF\tM\t1\t2\n",
    );
    let v = variant(
        1,
        "GENEA",
        "AD",
        "chr10",
        123_256_215,
        "T",
        "G",
        "0/0:0/1:0/1",
        "missense_variant",
        "GENEA:ENST1:c.1688A>C:p.(Glu563Ala)",
        "PATHOGENIC",
        &[],
    );
    let (code, err) = report(&d, &tsv(&[v], &HEADER));
    let words = fs::read_to_string(d.p("words.txt")).unwrap_or_default();
    assert!(
        code == 0
            && words.contains(
                "chr10:123256215 T>G [M does not carry it, F heterozygous, S1 heterozygous], missense_variant"
            ),
        "{err}{words}"
    );
    assert_eq!(read_rows(&d.p("top.tsv"))[0]["genotype"], "0/0:0/1:0/1");
}

#[test]
fn the_report_on_no_rank_and_on_damaged_tables() {
    let d = Dir::new("damaged");
    let (code, err) = report(&d, &tsv(&[], &HEADER));
    assert_eq!(code, 0, "{err}");
    assert!(read_rows(&d.p("top.tsv")).is_empty());
    assert!(
        fs::read_to_string(d.p("words.txt"))
            .unwrap()
            .contains("No gene passed")
    );
    let less: Vec<&str> = HEADER
        .iter()
        .copied()
        .filter(|h| *h != "EXOMISER_ACMG_CLASSIFICATION")
        .collect();
    let (code, err) = report(&d, &tsv(&rows(), &less));
    assert!(
        code == 1 && err.contains("EXOMISER_ACMG_CLASSIFICATION"),
        "{err}"
    );
    let (code, err) = report(&d, &format!("{}\n1\tx\tGENEA\n", HEADER.join("\t")));
    assert!(code == 1 && err.contains("one field per column"), "{err}");
    let (code, err) = report(&d, "");
    assert!(
        code == 1 && err.contains("not an Exomiser variants TSV"),
        "{err}"
    );
}

#[test]
fn the_modules_stub_table_has_the_reports_columns() {
    let d = Dir::new("stub");
    report(&d, &tsv(&rows(), &HEADER));
    let out_header = fs::read_to_string(d.p("top.tsv"))
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_owned();
    let nf = fs::read_to_string(MODULE).unwrap();
    let stub = nf.rsplit("stub:").next().unwrap();
    let line = stub
        .lines()
        .find(|l| l.contains("_exomiser_top.tsv"))
        .expect("the stub writes the top table");
    let head = line
        .split('\'')
        .nth(1)
        .unwrap()
        .replace("\\\\t", "\t")
        .replace("\\\\n", "\n");
    assert_eq!(head, format!("{out_header}\n"));
}
