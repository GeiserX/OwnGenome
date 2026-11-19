//! A recessive secondary finding needs two different copies of the gene.
//!
//! The ACMG SF tier used to add up the alt copies of a gene's records. So
//! two MUTYH variants that read-backed phase puts on one copy (cis), a
//! combined allele counted with its own parts, or a record with no alt copy
//! called could complete the two hits MUTYH needs. These cases run the tier
//! on synthetic hits with the shipped `assets/acmg_sf.tsv`.

use indexmap::IndexMap;
use og_summary::collect::acmg::{ClinicalVcf, HighRecord, Inputs, acmg_sf};
use og_summary::collect::clinvar::{Floor, HitRecords};
use og_summary::collect::two_hits::DeletionInputs;
use pycompat::json::Value;

const ASSETS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets");

// Two public MUTYH pathogenic sites (GRCh38).
const A: (i128, &str, &str) = (45_331_556, "C", "T");
const B: (i128, &str, &str) = (45_332_803, "T", "C");

/// One step 06 hit in `gene`, Pathogenic with 2 stars, GT and PS as given,
/// and `info` added to INFO.
fn hit(
    site: (i128, &str, &str),
    gene: &str,
    gt: &str,
    ps: Option<&str>,
    info: &str,
) -> Vec<String> {
    let (fmt, sample) = match ps {
        Some(p) => ("GT:GQ:DP:PS", format!("{gt}:50:30:{p}")),
        None => ("GT:GQ:DP", format!("{gt}:50:30")),
    };
    let info = format!(
        "GENEINFO={gene}:1;CLNSIG=Pathogenic;\
         CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts{info}"
    );
    [
        "chr1",
        &site.0.to_string(),
        ".",
        site.1,
        site.2,
        ".",
        "PASS",
        &info,
        fmt,
        &sample,
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

fn tier(
    hits: Vec<Vec<String>>,
    high: IndexMap<(String, i128, String, String), HighRecord>,
) -> Value {
    let h: Result<HitRecords, String> = Ok((Floor::default(), hits));
    let c: Result<ClinicalVcf, String> = Ok(ClinicalVcf {
        records: high.len(),
        high: Ok(high),
    });
    let d = no_deletions();
    acmg_sf(&Inputs {
        assets: ASSETS,
        hits: Some(&h),
        clinical: Some(&c),
        deletions: &d,
    })
    .unwrap()
}

fn rows<'a>(v: &'a Value, k: &str) -> Vec<&'a IndexMap<String, Value>> {
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

fn s<'a>(m: &'a IndexMap<String, Value>, k: &str) -> &'a str {
    match m.get(k) {
        Some(Value::Str(x)) => x,
        x => panic!("{k}: {x:?} in {m:?}"),
    }
}

type Row<'a> = &'a IndexMap<String, Value>;

/// The MUTYH rows of the reported ClinVar hits and of `not_reported`.
fn mutyh(v: &Value) -> (Vec<Row<'_>>, Vec<Row<'_>>) {
    let rep = rows(v, "clinvar_hits")
        .into_iter()
        .filter(|m| s(m, "gene") == "MUTYH")
        .collect();
    let apart = rows(v, "not_reported")
        .into_iter()
        .filter(|m| s(m, "gene") == "MUTYH")
        .collect();
    (rep, apart)
}

fn none() -> IndexMap<(String, i128, String, String), HighRecord> {
    IndexMap::new()
}

#[test]
fn two_variants_phased_on_one_copy_are_no_finding() {
    let v = tier(
        vec![
            hit(A, "MUTYH", "0|1", Some("7"), ""),
            hit(B, "MUTYH", "0|1", Some("7"), ""),
        ],
        none(),
    );
    let (rep, apart) = mutyh(&v);
    assert!(rep.is_empty(), "a cis pair was reported: {rep:?}");
    assert_eq!(apart.len(), 2, "{apart:?}");
    for m in apart {
        let r = s(m, "reason");
        assert!(
            r.starts_with("two heterozygous variants in cis") && r.contains("phase set 7"),
            "{r}"
        );
    }
}

#[test]
fn two_variants_phased_on_different_copies_are_biallelic() {
    let v = tier(
        vec![
            hit(A, "MUTYH", "0|1", Some("7"), ""),
            hit(B, "MUTYH", "1|0", Some("7"), ""),
        ],
        none(),
    );
    let (rep, apart) = mutyh(&v);
    assert!(apart.is_empty(), "{apart:?}");
    assert_eq!(rep.len(), 2);
    for m in rep {
        assert_eq!(s(m, "pair"), "biallelic");
        assert!(s(m, "phase").starts_with("in trans"), "{}", s(m, "phase"));
    }
}

#[test]
fn unphased_or_unlinked_pairs_are_reported_as_possibly_biallelic() {
    for (ga, pa, gb, pb) in [
        ("0/1", None, "0/1", None),
        // Two phase sets link nothing.
        ("0|1", Some("7"), "0|1", Some("8")),
        // A "|" GT with no phase set is not phased.
        ("0|1", None, "0|1", None),
        // A half call never puts its alt on a known copy.
        ("1|.", Some("7"), "0|1", Some("7")),
    ] {
        let v = tier(
            vec![hit(A, "MUTYH", ga, pa, ""), hit(B, "MUTYH", gb, pb, "")],
            none(),
        );
        let (rep, apart) = mutyh(&v);
        assert!(apart.is_empty(), "{ga} {gb}: {apart:?}");
        assert_eq!(rep.len(), 2, "{ga} {gb}");
        for m in rep {
            assert_eq!(s(m, "pair"), "possibly biallelic", "{ga} {gb}");
            assert!(s(m, "phase").starts_with("phase unknown"), "{ga} {gb}");
        }
    }
}

#[test]
fn the_two_alt_alleles_of_one_split_call_are_biallelic() {
    let v = tier(
        vec![
            hit((A.0, "C", "T"), "MUTYH", "1/.", None, ""),
            hit((A.0, "C", "G"), "MUTYH", "./1", None, ""),
        ],
        none(),
    );
    let (rep, _) = mutyh(&v);
    assert_eq!(rep.len(), 2);
    for m in rep {
        assert_eq!(s(m, "pair"), "biallelic");
        assert!(
            s(m, "phase").starts_with("in trans: the two alt alleles of one call"),
            "{}",
            s(m, "phase")
        );
    }
}

#[test]
fn a_homozygote_is_biallelic_and_one_heterozygote_is_not() {
    let v = tier(vec![hit(A, "MUTYH", "1/1", None, "")], none());
    let (rep, _) = mutyh(&v);
    assert_eq!(rep.len(), 1);
    assert_eq!(s(rep[0], "pair"), "biallelic");
    assert!(s(rep[0], "phase").starts_with("homozygous"));
    let v = tier(vec![hit(A, "MUTYH", "0/1", None, "")], none());
    let (rep, apart) = mutyh(&v);
    assert!(rep.is_empty());
    assert!(s(apart[0], "reason").starts_with("one heterozygous variant;"));
}

#[test]
fn a_record_with_no_alt_copy_called_completes_no_pair() {
    for gt in ["./.", "0/0"] {
        let v = tier(
            vec![
                hit(A, "MUTYH", "0/1", None, ""),
                hit(B, "MUTYH", gt, None, ""),
            ],
            none(),
        );
        let (rep, apart) = mutyh(&v);
        assert!(rep.is_empty(), "GT {gt}: {rep:?}");
        assert_eq!(apart.len(), 2);
        let none_row = apart
            .iter()
            .find(|m| s(m, "variant").contains("45332803"))
            .unwrap();
        assert!(
            s(none_row, "reason").starts_with("no alt copy called"),
            "{gt}"
        );
    }
}

// A combined allele of step 06 and its parts (chr1:100 AGC>A and chr1:105
// T>TGGG joined as chr1:101 GCCAT>CATGGG).
const PART1: (i128, &str, &str) = (100, "AGC", "A");
const PART2: (i128, &str, &str) = (105, "T", "TGGG");
const WHOLE: (i128, &str, &str) = (101, "GCCAT", "CATGGG");
const JOINED: &str = ";PGP_JOINED_FROM=chr1:100:AGC:A|chr1:105:T:TGGG";
const PART_OF: &str = ";PGP_PART_OF=chr1:101:GCCAT:CATGGG";

#[test]
fn a_combined_allele_counts_once_with_its_part() {
    let v = tier(
        vec![
            hit(PART1, "MUTYH", "0/1", None, PART_OF),
            hit(WHOLE, "MUTYH", "0/1", None, JOINED),
        ],
        none(),
    );
    let (rep, apart) = mutyh(&v);
    assert!(rep.is_empty(), "one allele was reported as two: {rep:?}");
    assert_eq!(apart.len(), 2);
    assert!(
        s(apart[0], "reason").contains("a combined allele and its parts count once"),
        "{}",
        s(apart[0], "reason")
    );
}

#[test]
fn a_phased_combined_allele_counts_once_with_both_parts() {
    let v = tier(
        vec![
            hit(PART1, "MUTYH", "0|1", Some("90"), PART_OF),
            hit(WHOLE, "MUTYH", "0|1", Some("90"), JOINED),
            hit(PART2, "MUTYH", "0|1", Some("90"), PART_OF),
        ],
        none(),
    );
    let (rep, apart) = mutyh(&v);
    assert!(rep.is_empty(), "{rep:?}");
    assert_eq!(apart.len(), 3);
}

#[test]
fn an_unphased_combined_allele_with_two_candidate_parts_counts_its_parts() {
    // The two parts may sit on different copies: possibly biallelic.
    let v = tier(
        vec![
            hit(PART1, "MUTYH", "0/1", None, PART_OF),
            hit(WHOLE, "MUTYH", "0/1", None, JOINED),
            hit(PART2, "MUTYH", "0/1", None, PART_OF),
        ],
        none(),
    );
    let (rep, apart) = mutyh(&v);
    assert!(apart.is_empty(), "{apart:?}");
    assert_eq!(rep.len(), 3);
    assert!(rep.iter().all(|m| s(m, "pair") == "possibly biallelic"));
    let whole = rep
        .iter()
        .find(|m| s(m, "variant").contains("GCCAT"))
        .unwrap();
    assert!(
        s(whole, "combined_allele")
            .ends_with("unphased, so the parts may sit on different haplotypes"),
        "{}",
        s(whole, "combined_allele")
    );
}

#[test]
fn a_part_found_as_a_high_impact_record_counts_once_with_the_combined_allele() {
    let mut genes = IndexMap::new();
    genes.insert("MUTYH".to_owned(), "frameshift_variant".to_owned());
    let mut high = IndexMap::new();
    high.insert(
        (
            "chr1".to_owned(),
            PART1.0,
            PART1.1.to_owned(),
            PART1.2.to_owned(),
        ),
        HighRecord {
            genotype: "0/1".into(),
            ps: None,
            genes,
            transcript: "MANE",
            frequency: "rare",
            weak: false,
            why: Vec::new(),
            gq: Some(50),
            dp: Some(30),
            ad: None,
        },
    );
    let v = tier(vec![hit(WHOLE, "MUTYH", "0/1", None, JOINED)], high);
    let (rep, _) = mutyh(&v);
    assert!(rep.is_empty(), "{rep:?}");
    assert!(
        rows(&v, "high_impact")
            .iter()
            .all(|m| s(m, "gene") != "MUTYH"),
        "the part was reported as a second hit"
    );
}

#[test]
fn a_combined_allele_with_an_unrelated_variant_is_still_two() {
    let v = tier(
        vec![
            hit(PART1, "MUTYH", "0/1", None, PART_OF),
            hit(WHOLE, "MUTYH", "0/1", None, JOINED),
            hit(A, "MUTYH", "0/1", None, ""),
        ],
        none(),
    );
    let (rep, apart) = mutyh(&v);
    assert!(apart.is_empty(), "{apart:?}");
    assert_eq!(rep.len(), 3);
}

#[test]
fn a_hemizygous_rule_gene_also_counts_copies_not_records() {
    // ABCD1 (rule hemizygous): two heterozygous variants of a woman in cis
    // are one copy.
    let v = tier(
        vec![
            hit(A, "ABCD1", "0|1", Some("3"), ""),
            hit(B, "ABCD1", "0|1", Some("3"), ""),
        ],
        none(),
    );
    assert!(
        rows(&v, "clinvar_hits")
            .iter()
            .all(|m| s(m, "gene") != "ABCD1")
    );
    let apart: Vec<_> = rows(&v, "not_reported")
        .into_iter()
        .filter(|m| s(m, "gene") == "ABCD1")
        .collect();
    assert_eq!(apart.len(), 2);
    assert!(s(apart[0], "reason").contains("in cis"));
}
