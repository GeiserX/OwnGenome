//! The carrier list keeps each hit's call quality and zygosity, and counts
//! two hits in one gene as copies.
//!
//! The list used to judge every hit at GQ 30 and DP 10, whatever floor the
//! hits file states, and it dropped the low-confidence mark the ClinVar
//! section shows for the same row. It listed a homozygous or hemizygous
//! P/LP row as a carrier beside the heterozygous ones, two different
//! heterozygous hits of one recessive gene as two plain carrier rows, and a
//! deletion plus a low-confidence call as a confident two-hit gene. These
//! cases run the list on synthetic hits with the shipped
//! `assets/carrier_genes.tsv`. The last cases run the ACMG SF card on the
//! same hits (`assets/acmg_sf.tsv`): a pair that rests on a low-confidence
//! call reads the same there.

use indexmap::IndexMap;
use og_summary::collect::acmg::{Inputs, acmg_sf};
use og_summary::collect::carrier::Carrier;
use og_summary::collect::clinvar::{Floor, HitRecords, HitsVcf};
use og_summary::collect::two_hits::{DeletionInputs, SvRow};
use pycompat::json::Value;

const ASSETS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets");

// MUTYH and CFTR spans in assets/carrier_genes.tsv: chr1:45329163-45340893
// and chr7:117287120-117715971. MUTYH is on the ACMG SF list, CFTR is not.
const MUTYH_A: (&str, i128, &str, &str) = ("chr1", 45_331_556, "C", "T");
const MUTYH_B: (&str, i128, &str, &str) = ("chr1", 45_332_803, "T", "C");
const CFTR_A: (&str, i128, &str, &str) = ("chr7", 117_559_590, "ATCT", "A");
const CFTR_B: (&str, i128, &str, &str) = ("chr7", 117_587_806, "G", "T");
const OTC_A: (&str, i128, &str, &str) = ("chrX", 38_400_000, "C", "T");

/// One step 06 hit: Pathogenic with 2 stars, the FORMAT values given (GT
/// first), and `info` added to INFO.
fn hit(
    site: (&str, i128, &str, &str),
    gene: &str,
    fmt: &str,
    sample: &str,
    info: &str,
) -> Vec<String> {
    let info = format!(
        "GENEINFO={gene}:1;CLNSIG=Pathogenic;\
         CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts{info}"
    );
    [
        site.0,
        &site.1.to_string(),
        ".",
        site.2,
        site.3,
        ".",
        "PASS",
        &info,
        fmt,
        sample,
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect()
}

/// A hit with good call quality (GQ 50, DP 30).
fn good(site: (&str, i128, &str, &str), gene: &str, gt: &str, ps: Option<&str>) -> Vec<String> {
    match ps {
        Some(p) => hit(site, gene, "GT:GQ:DP:PS", &format!("{gt}:50:30:{p}"), ""),
        None => hit(site, gene, "GT:GQ:DP", &format!("{gt}:50:30"), ""),
    }
}

fn no_deletions() -> DeletionInputs {
    DeletionInputs {
        annotsv: Ok(None),
        step35: Ok(None),
    }
}

/// A class 5 MUTYH deletion, chr1:45335000-45337000, GT and PS as given.
fn mutyh_deletion(gt: &str, ps: Option<&str>) -> DeletionInputs {
    DeletionInputs {
        annotsv: Ok(Some(vec![SvRow {
            sv: "chr1:45335000-45337000".into(),
            sv_type: Some("DEL".into()),
            gt: Some(gt.into()),
            ps: ps.map(str::to_owned),
            acmg_class: Some("5".into()),
            genes: vec!["MUTYH".into()],
            caller: Some("manta".into()),
            support: None,
            common: false,
        }])),
        step35: Ok(None),
    }
}

fn list(hits: Vec<Vec<String>>, dels: &DeletionInputs) -> IndexMap<String, Value> {
    let h: Result<HitRecords, String> = Ok((Floor::default(), hits));
    obj(Carrier::read(ASSETS, &h, dels).unwrap().to_json())
}

fn obj(v: Value) -> IndexMap<String, Value> {
    match v {
        Value::Object(m) => m,
        x => panic!("{x:?}"),
    }
}

fn rows<'a>(m: &'a IndexMap<String, Value>, k: &str) -> Vec<&'a IndexMap<String, Value>> {
    let Some(Value::Array(a)) = m.get(k) else {
        panic!("{k}: {m:?}")
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

fn of<'a>(m: &'a IndexMap<String, Value>, k: &str, gene: &str) -> Vec<&'a IndexMap<String, Value>> {
    rows(m, k)
        .into_iter()
        .filter(|r| s(r, "gene") == gene)
        .collect()
}

#[test]
fn a_weak_call_keeps_its_low_confidence_mark() {
    let c = list(
        vec![hit(MUTYH_A, "MUTYH", "GT:GQ:DP", "0/1:1:1", "")],
        &no_deletions(),
    );
    let r = of(&c, "carriers", "MUTYH");
    assert_eq!(r.len(), 1, "{c:?}");
    assert_eq!(
        r[0].get("low_confidence"),
        Some(&Value::Bool(true)),
        "{r:?}"
    );
    assert_eq!(r[0].get("gq"), Some(&Value::Int(1)));
    assert_eq!(r[0].get("dp"), Some(&Value::Int(1)));
    assert_eq!(s(r[0], "zygosity"), "het");

    let c = list(vec![good(MUTYH_A, "MUTYH", "0/1", None)], &no_deletions());
    let r = of(&c, "carriers", "MUTYH");
    assert_eq!(
        r[0].get("low_confidence"),
        Some(&Value::Bool(false)),
        "{r:?}"
    );
}

#[test]
fn the_hits_file_floor_judges_the_carrier_rows() {
    let dir = std::env::temp_dir().join(format!("og-i411-carrier-floor-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("S_clinvar_hits.vcf");
    let rec = hit(MUTYH_A, "MUTYH", "GT:GQ:DP", "0/1:35:20", "").join("\t");
    std::fs::write(
        &path,
        format!(
            "##fileformat=VCFv4.2\n##PGP_CALL_FLOOR=GQ=40,DP=15\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n{rec}\n"
        ),
    )
    .unwrap();
    let p = path.to_string_lossy().into_owned();
    let vcf = HitsVcf::read(p.clone(), p);
    std::fs::remove_dir_all(&dir).unwrap();
    assert_eq!(vcf.read.as_ref().unwrap().0, Floor { gq: 40, dp: 15 });
    let c = obj(Carrier::read(ASSETS, &vcf.read, &no_deletions())
        .unwrap()
        .to_json());
    let r = of(&c, "carriers", "MUTYH");
    assert_eq!(r.len(), 1, "{c:?}");
    assert_eq!(
        r[0].get("low_confidence"),
        Some(&Value::Bool(true)),
        "GQ 35 is below the file's floor of GQ 40: {r:?}"
    );
}

#[test]
fn a_homozygous_or_hemizygous_row_is_no_carrier() {
    let c = list(
        vec![
            good(CFTR_A, "CFTR", "1/1", None),
            good(OTC_A, "OTC", "1", None),
            good(MUTYH_A, "MUTYH", "0/1", None),
        ],
        &no_deletions(),
    );
    let carriers: Vec<&str> = rows(&c, "carriers").iter().map(|r| s(r, "gene")).collect();
    assert_eq!(carriers, ["MUTYH"], "{c:?}");
    let both: Vec<(&str, &str)> = rows(&c, "hom_or_hemizygous")
        .iter()
        .map(|r| (s(r, "gene"), s(r, "zygosity")))
        .collect();
    assert_eq!(both, [("CFTR", "hom"), ("OTC", "hemizygous")]);
}

#[test]
fn a_homozygous_call_inside_a_deletion_is_listed_hemizygous_apart() {
    let inside = ("chr1", 45_336_000, "G", "A");
    let c = list(
        vec![good(inside, "MUTYH", "1/1", None)],
        &mutyh_deletion("0/1", None),
    );
    assert!(of(&c, "carriers", "MUTYH").is_empty(), "{c:?}");
    let r = of(&c, "hom_or_hemizygous", "MUTYH");
    assert_eq!(r.len(), 1, "{c:?}");
    assert_eq!(s(r[0], "genotype"), "hemizygous");
    assert_eq!(s(r[0], "zygosity"), "hemizygous");
}

#[test]
fn two_unphased_hets_in_one_recessive_gene_are_possible_two_hits() {
    let c = list(
        vec![
            good(CFTR_A, "CFTR", "0/1", None),
            good(CFTR_B, "CFTR", "0/1", None),
        ],
        &no_deletions(),
    );
    let r = of(&c, "carriers", "CFTR");
    assert_eq!(r.len(), 2, "{c:?}");
    for x in r {
        assert_eq!(s(x, "pair"), "possible two hits", "{x:?}");
        assert!(s(x, "phase").starts_with("phase unknown"), "{x:?}");
    }
}

#[test]
fn read_backed_phase_decides_two_hits_or_one_copy() {
    let trans = list(
        vec![
            good(CFTR_A, "CFTR", "0|1", Some("9")),
            good(CFTR_B, "CFTR", "1|0", Some("9")),
        ],
        &no_deletions(),
    );
    for x in of(&trans, "carriers", "CFTR") {
        assert_eq!(s(x, "pair"), "two hits", "{x:?}");
        assert!(s(x, "phase").contains("phase set 9"), "{x:?}");
    }
    let cis = list(
        vec![
            good(CFTR_A, "CFTR", "0|1", Some("9")),
            good(CFTR_B, "CFTR", "0|1", Some("9")),
        ],
        &no_deletions(),
    );
    for x in of(&cis, "carriers", "CFTR") {
        assert!(x.get("pair").is_none(), "a cis pair is one copy: {x:?}");
        assert!(
            s(x, "phase").starts_with("two heterozygous variants in cis"),
            "{x:?}"
        );
    }
    // A phased GT without a phase set says nothing about the other call.
    let no_ps = list(
        vec![
            good(CFTR_A, "CFTR", "0|1", None),
            good(CFTR_B, "CFTR", "1|0", None),
        ],
        &no_deletions(),
    );
    for x in of(&no_ps, "carriers", "CFTR") {
        assert_eq!(s(x, "pair"), "possible two hits", "{x:?}");
    }
}

#[test]
fn a_low_confidence_call_makes_two_hits_only_possible() {
    let c = list(
        vec![
            good(CFTR_A, "CFTR", "0|1", Some("9")),
            hit(CFTR_B, "CFTR", "GT:GQ:DP:PS", "1|0:5:4:9", ""),
        ],
        &no_deletions(),
    );
    assert_eq!(of(&c, "carriers", "CFTR").len(), 2, "{c:?}");
    for x in of(&c, "carriers", "CFTR") {
        assert_eq!(
            s(x, "pair"),
            "possible two hits (low-confidence call)",
            "{x:?}"
        );
    }
}

#[test]
fn a_combined_allele_and_its_part_are_one_allele() {
    let joined = format!(
        ";PGP_JOINED_FROM={}:{}:{}:{}|chr7:117559600:G:A",
        CFTR_A.0, CFTR_A.1, CFTR_A.2, CFTR_A.3
    );
    let c = list(
        vec![
            good(CFTR_A, "CFTR", "0/1", None),
            hit(CFTR_B, "CFTR", "GT:GQ:DP", "0/1:50:30", &joined),
        ],
        &no_deletions(),
    );
    assert_eq!(of(&c, "carriers", "CFTR").len(), 2, "{c:?}");
    for x in of(&c, "carriers", "CFTR") {
        assert!(x.get("pair").is_none(), "{x:?}");
        assert!(x.get("phase").is_none(), "{x:?}");
    }
}

fn two_hit_verdict(c: &IndexMap<String, Value>, gene: &str) -> String {
    let Some(Value::Object(t)) = c.get("two_hits") else {
        panic!("{c:?}")
    };
    let g = rows(t, "genes");
    let g: Vec<_> = g.into_iter().filter(|r| s(r, "gene") == gene).collect();
    assert_eq!(g.len(), 1, "{t:?}");
    s(g[0], "verdict").to_owned()
}

#[test]
fn a_deletion_plus_a_low_confidence_call_is_only_possibly_biallelic() {
    let dels = mutyh_deletion("1|0", Some("5"));
    let sure = list(vec![good(MUTYH_A, "MUTYH", "0|1", Some("5"))], &dels);
    assert_eq!(two_hit_verdict(&sure, "MUTYH"), "biallelic");
    let weak = list(
        vec![hit(MUTYH_A, "MUTYH", "GT:GQ:DP:PS", "0|1:12:6:5", "")],
        &dels,
    );
    assert_eq!(
        two_hit_verdict(&weak, "MUTYH"),
        "possibly biallelic (low-confidence call)"
    );
}

/// The ACMG SF card on step 06's hits alone (no clinical VCF).
fn sf(hits: Vec<Vec<String>>, dels: &DeletionInputs) -> IndexMap<String, Value> {
    let h: Result<HitRecords, String> = Ok((Floor::default(), hits));
    obj(acmg_sf(&Inputs {
        assets: ASSETS,
        hits: Some(&h),
        clinical: None,
        deletions: dels,
    })
    .unwrap())
}

/// The `pair` word of each reported ClinVar row of `gene` on the SF card.
fn sf_pairs(c: &IndexMap<String, Value>, gene: &str) -> Vec<String> {
    let r = of(c, "clinvar_hits", gene);
    assert!(!r.is_empty(), "{c:?}");
    r.iter().map(|x| s(x, "pair").to_owned()).collect()
}

#[test]
fn the_sf_card_caps_a_deletion_pair_on_a_low_confidence_call() {
    let dels = mutyh_deletion("1|0", Some("9"));
    let weak = sf(
        vec![hit(MUTYH_A, "MUTYH", "GT:GQ:DP:PS", "0|1:5:4:9", "")],
        &dels,
    );
    assert_eq!(
        sf_pairs(&weak, "MUTYH"),
        ["possibly biallelic (low-confidence call)"]
    );
    let sure = sf(vec![good(MUTYH_A, "MUTYH", "0|1", Some("9"))], &dels);
    assert_eq!(sf_pairs(&sure, "MUTYH"), ["biallelic"]);
}

#[test]
fn the_sf_card_caps_a_phased_pair_on_a_low_confidence_call() {
    let weak = sf(
        vec![
            good(MUTYH_A, "MUTYH", "0|1", Some("9")),
            hit(MUTYH_B, "MUTYH", "GT:GQ:DP:PS", "1|0:5:4:9", ""),
        ],
        &no_deletions(),
    );
    let low = "possibly biallelic (low-confidence call)";
    assert_eq!(sf_pairs(&weak, "MUTYH"), [low, low]);
    for x in of(&weak, "clinvar_hits", "MUTYH") {
        assert!(s(x, "phase").contains("low confidence"), "{x:?}");
    }
    // Unphased, one of them weak: the same mark.
    let unphased = sf(
        vec![
            good(MUTYH_A, "MUTYH", "0/1", None),
            hit(MUTYH_B, "MUTYH", "GT:GQ:DP", "0/1:5:4", ""),
        ],
        &no_deletions(),
    );
    assert_eq!(sf_pairs(&unphased, "MUTYH"), [low, low]);
    // Both confident: biallelic, and possibly biallelic unphased.
    let sure = sf(
        vec![
            good(MUTYH_A, "MUTYH", "0|1", Some("9")),
            good(MUTYH_B, "MUTYH", "1|0", Some("9")),
        ],
        &no_deletions(),
    );
    assert_eq!(sf_pairs(&sure, "MUTYH"), ["biallelic", "biallelic"]);
    let open = sf(
        vec![
            good(MUTYH_A, "MUTYH", "0/1", None),
            good(MUTYH_B, "MUTYH", "0/1", None),
        ],
        &no_deletions(),
    );
    assert_eq!(
        sf_pairs(&open, "MUTYH"),
        ["possibly biallelic", "possibly biallelic"]
    );
    // A confident homozygous call proves both copies on its own; a weak
    // second variant does not weaken that.
    let hom = sf(
        vec![
            good(MUTYH_A, "MUTYH", "1/1", None),
            hit(MUTYH_B, "MUTYH", "GT:GQ:DP", "0/1:5:4", ""),
        ],
        &no_deletions(),
    );
    assert_eq!(sf_pairs(&hom, "MUTYH"), ["biallelic", "biallelic"]);
}
