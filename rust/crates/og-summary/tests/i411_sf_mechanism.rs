//! A predicted loss-of-function call counts as a secondary finding only in
//! genes where loss of function causes the disease, and TTN and LMNA keep
//! their Table 1 limits.
//!
//! The ACMG SF tier used to count a rare VEP HIGH-impact record in every
//! gene of the list, so a PCSK9 stop-gain (PCSK9 causes familial
//! hypercholesterolemia by gain of function) read as a finding. A call with
//! no GQ, DP or AD counted as a sound rare call. TTN's "truncating variants
//! only" and LMNA's "missense only with cardiac evidence" were notes, not
//! rules. These cases run the tier on a synthetic clinical VCF, read from a
//! file as the collector reads it, with the shipped `assets/acmg_sf.tsv`.

use indexmap::IndexMap;
use og_summary::collect::acmg::{ClinicalVcf, Inputs, acmg_sf};
use og_summary::collect::clinvar::{Floor, HitRecords};
use og_summary::collect::two_hits::DeletionInputs;
use pycompat::json::Value;
use std::sync::atomic::{AtomicUsize, Ordering};

const ASSETS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets");

const HEADER: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. \
Format: Allele|Consequence|IMPACT|SYMBOL|Feature|CANONICAL|MANE_SELECT|MAX_AF\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

/// A sound heterozygous call.
const OK: (&str, &str) = ("GT:GQ:DP:AD", "0/1:50:30:15,15");

// GRCh38 positions: PCSK9, BRCA2, LMNA, TGFBR2 and MUTYH inside their
// genes; TTN in an A-band exon (exon 300 of ENST00000589042 lies near
// 178.60 Mb) and in the I-band (178.65 Mb).
const PCSK9: (&str, i128) = ("chr1", 55_051_215);
const BRCA2: (&str, i128) = ("chr13", 32_340_300);
const TTN_A: (&str, i128) = ("chr2", 178_600_000);
const TTN_I: (&str, i128) = ("chr2", 178_650_000);
const LMNA: (&str, i128) = ("chr1", 156_130_000);
const TGFBR2: (&str, i128) = ("chr3", 30_650_000);
const MUTYH_A: (&str, i128) = ("chr1", 45_331_556);
const MUTYH_B: (&str, i128) = ("chr1", 45_332_803);

/// One clinical VCF record: `cons|impact|gene` on the MANE transcript with
/// `max_af` (empty: unknown frequency), and the call.
fn record(
    at: (&str, i128),
    cons: &str,
    impact: &str,
    gene: &str,
    max_af: &str,
    call: (&str, &str),
) -> String {
    format!(
        "{}\t{}\t.\tC\tT\t50\tPASS\tCSQ=T|{cons}|{impact}|{gene}|ENST1|YES|NM_1|{max_af}\t{}\t{}\n",
        at.0, at.1, call.0, call.1
    )
}

fn stop(at: (&str, i128), gene: &str) -> String {
    record(at, "stop_gained", "HIGH", gene, "0.0001", OK)
}

static N: AtomicUsize = AtomicUsize::new(0);

/// The clinical VCF of `records`, read from a file.
fn clinical(records: &[String]) -> Result<ClinicalVcf, String> {
    let dir = std::env::temp_dir().join(format!(
        "og-i411-sf-mechanism-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("S_clinical.vcf");
    std::fs::write(&p, format!("{HEADER}{}", records.concat())).unwrap();
    let c = ClinicalVcf::read(&p.display().to_string());
    std::fs::remove_dir_all(&dir).unwrap();
    c
}

/// One step 06 hit: Pathogenic, 2 stars, a sound het call.
fn hit(at: (&str, i128), gene: &str) -> Vec<String> {
    [
        at.0,
        &at.1.to_string(),
        ".",
        "C",
        "T",
        ".",
        "PASS",
        &format!(
            "GENEINFO={gene}:1;CLNSIG=Pathogenic;\
             CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts"
        ),
        "GT:GQ:DP",
        "0/1:50:30",
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect()
}

/// The tier from step 06 `hits` and the clinical VCF (`None`: no VCF).
fn tier(hits: Vec<Vec<String>>, records: Option<&[String]>) -> Value {
    let h: Result<HitRecords, String> = Ok((Floor::default(), hits));
    let c = records.map(clinical);
    let d = DeletionInputs {
        annotsv: Ok(None),
        step35: Ok(None),
    };
    acmg_sf(&Inputs {
        assets: ASSETS,
        hits: Some(&h),
        clinical: c.as_ref(),
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

fn s<'a>(m: &'a Row, k: &str) -> &'a str {
    match m.get(k) {
        Some(Value::Str(x)) => x,
        x => panic!("{k}: {x:?} in {m:?}"),
    }
}

fn int(v: &Value, k: &str) -> Option<i128> {
    let Value::Object(m) = v else { panic!("{v:?}") };
    match m.get(k) {
        Some(Value::Int(n)) => Some(*n),
        None => None,
        x => panic!("{k}: {x:?}"),
    }
}

fn of<'a>(v: &'a Value, list: &str, gene: &str) -> Vec<&'a Row> {
    rows(v, list)
        .into_iter()
        .filter(|m| s(m, "gene") == gene)
        .collect()
}

#[test]
fn a_pcsk9_stop_gain_is_no_secondary_finding() {
    let v = tier(Vec::new(), Some(&[stop(PCSK9, "PCSK9")]));
    assert!(of(&v, "high_impact", "PCSK9").is_empty(), "{v:?}");
    let apart = of(&v, "not_reported", "PCSK9");
    assert_eq!(apart.len(), 1, "{v:?}");
    let r = s(apart[0], "reason");
    assert!(
        r.contains("does not count in PCSK9") && r.contains("ClinVar P/LP"),
        "{r}"
    );
    assert_eq!(s(apart[0], "source"), "HIGH impact");
    assert_eq!(int(&v, "high_impact_rare"), Some(0));
    // A ClinVar P/LP variant of PCSK9 still counts: only the predicted
    // loss of function is set aside.
    let v = tier(vec![hit(PCSK9, "PCSK9")], Some(&[]));
    assert_eq!(of(&v, "clinvar_hits", "PCSK9").len(), 1, "{v:?}");
}

#[test]
fn a_brca2_stop_gain_is_one() {
    let v = tier(Vec::new(), Some(&[stop(BRCA2, "BRCA2")]));
    let rep = of(&v, "high_impact", "BRCA2");
    assert_eq!(rep.len(), 1, "{v:?}");
    assert!(!rep[0].contains_key("call_unverified"));
    assert!(!rep[0].contains_key("mechanism_note"));
    assert_eq!(int(&v, "high_impact_rare"), Some(1));
    assert_eq!(int(&v, "high_impact_unverified"), None);
}

#[test]
fn a_call_with_no_gq_dp_or_ad_is_unverified() {
    for call in [("GT", "0/1"), ("GT:GQ:DP:AD", "0/1:.:.:.")] {
        let v = tier(
            Vec::new(),
            Some(&[record(
                BRCA2,
                "stop_gained",
                "HIGH",
                "BRCA2",
                "0.0001",
                call,
            )]),
        );
        let rep = of(&v, "high_impact", "BRCA2");
        assert_eq!(rep.len(), 1, "{call:?}: {v:?}");
        assert_eq!(
            rep[0].get("call_unverified"),
            Some(&Value::Bool(true)),
            "{call:?}"
        );
        assert_eq!(rep[0].get("weak_call"), Some(&Value::Bool(false)));
        assert_eq!(int(&v, "high_impact_rare"), Some(0), "{call:?}");
        assert_eq!(int(&v, "high_impact_unverified"), Some(1), "{call:?}");
    }
}

#[test]
fn ttn_counts_a_truncating_variant_in_its_a_band_exons_only() {
    let v = tier(Vec::new(), Some(&[stop(TTN_I, "TTN")]));
    assert!(of(&v, "high_impact", "TTN").is_empty(), "{v:?}");
    let r = s(of(&v, "not_reported", "TTN")[0], "reason").to_owned();
    assert!(r.contains("footnote f") && r.contains("outside"), "{r}");

    let v = tier(Vec::new(), Some(&[stop(TTN_A, "TTN")]));
    assert_eq!(of(&v, "high_impact", "TTN").len(), 1, "{v:?}");
    assert!(of(&v, "not_reported", "TTN").is_empty());

    // A ClinVar hit in the A-band counts when its record in the clinical
    // VCF truncates, not when it is a missense one.
    let missense = record(TTN_A, "missense_variant", "MODERATE", "TTN", "0.0001", OK);
    let v = tier(vec![hit(TTN_A, "TTN")], Some(&[missense]));
    assert!(of(&v, "clinvar_hits", "TTN").is_empty(), "{v:?}");
    let r = s(of(&v, "not_reported", "TTN")[0], "reason").to_owned();
    assert!(r.contains("no truncating call"), "{r}");
    let v = tier(vec![hit(TTN_A, "TTN")], Some(&[stop(TTN_A, "TTN")]));
    let rep = of(&v, "clinvar_hits", "TTN");
    assert_eq!(rep.len(), 1, "{v:?}");
    assert!(!rep[0].contains_key("rule_condition"));
    // Outside the A-band a ClinVar hit is no finding either.
    let v = tier(vec![hit(TTN_I, "TTN")], Some(&[stop(TTN_I, "TTN")]));
    assert!(of(&v, "clinvar_hits", "TTN").is_empty(), "{v:?}");
    // With no clinical VCF the consequence is not known: reported, with the
    // condition to check.
    let v = tier(vec![hit(TTN_A, "TTN")], None);
    let rep = of(&v, "clinvar_hits", "TTN");
    assert_eq!(rep.len(), 1, "{v:?}");
    assert!(s(rep[0], "rule_condition").contains("not known here"));
}

#[test]
fn an_lmna_variant_that_may_be_missense_carries_the_case_evidence_condition() {
    let v = tier(vec![hit(LMNA, "LMNA")], Some(&[]));
    let rep = of(&v, "clinvar_hits", "LMNA");
    assert_eq!(rep.len(), 1, "{v:?}");
    let c = s(rep[0], "rule_condition");
    assert!(c.contains("footnote e") && c.contains("cardiac"), "{c}");
    // A truncating one is reported as it is.
    let v = tier(vec![hit(LMNA, "LMNA")], Some(&[stop(LMNA, "LMNA")]));
    for list in ["clinvar_hits", "high_impact"] {
        let rep = of(&v, list, "LMNA");
        assert_eq!(rep.len(), 1, "{list}: {v:?}");
        assert!(!rep[0].contains_key("rule_condition"), "{list}");
    }
}

#[test]
fn tgfbr2_counts_a_truncating_call_with_a_mark() {
    let v = tier(Vec::new(), Some(&[stop(TGFBR2, "TGFBR2")]));
    let rep = of(&v, "high_impact", "TGFBR2");
    assert_eq!(rep.len(), 1, "{v:?}");
    assert!(s(rep[0], "mechanism_note").contains("missense"));
}

#[test]
fn a_pair_that_needs_a_call_of_unknown_frequency_says_so() {
    let unknown = record(MUTYH_B, "frameshift_variant", "HIGH", "MUTYH", "", OK);
    let v = tier(vec![hit(MUTYH_A, "MUTYH")], Some(&[unknown]));
    let cv = of(&v, "clinvar_hits", "MUTYH");
    assert_eq!(cv.len(), 1, "{v:?}");
    assert_eq!(s(cv[0], "pair"), "possibly biallelic");
    assert!(s(cv[0], "pair_note").contains("unknown frequency"));
    let hi = of(&v, "high_impact", "MUTYH");
    assert_eq!(hi[0].get("frequency_unknown"), Some(&Value::Bool(true)));
    assert!(hi[0].contains_key("pair_note"));
    // A sound rare call needs no note.
    let v = tier(vec![hit(MUTYH_A, "MUTYH")], Some(&[stop(MUTYH_B, "MUTYH")]));
    assert!(!of(&v, "clinvar_hits", "MUTYH")[0].contains_key("pair_note"));
}

#[test]
fn a_weak_call_never_completes_a_recessive_pair() {
    let weak = record(
        MUTYH_B,
        "stop_gained",
        "HIGH",
        "MUTYH",
        "0.0001",
        ("GT:GQ:DP:AD", "0/1:8:30:15,15"),
    );
    let v = tier(vec![hit(MUTYH_A, "MUTYH")], Some(&[weak]));
    for m in of(&v, "clinvar_hits", "MUTYH") {
        assert_eq!(
            s(m, "pair"),
            "possibly biallelic (low-confidence call)",
            "{m:?}"
        );
    }
}
