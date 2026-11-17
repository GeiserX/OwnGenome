//! A deletion the sample was never called for pairs with nothing.
//!
//! AnnotSV writes a row for every SV of its input, with the sample's own
//! GT. A merged or multi-sample input can hold a deletion whose GT for this
//! sample is `./.` or missing: the sample was not called for it. Such a
//! row is no copy loss, so a homozygous call inside its span must not read
//! as hemizygous on the one copy the deletion leaves.

use indexmap::IndexMap;
use og_summary::collect::two_hits::{
    DeletionInputs, Spans, SvRow, Variant, Verdict, gene_deletions, two_hits,
};

/// AnnotSV class 5 DEL chr1:1000-5000 over the recessive gene G, with the
/// sample's GT.
fn annotsv_del(gt: Option<&str>) -> DeletionInputs {
    DeletionInputs {
        annotsv: Ok(Some(vec![SvRow {
            sv: "1:1000-5000".into(),
            sv_type: Some("DEL".into()),
            gt: gt.map(str::to_owned),
            ps: None,
            acmg_class: Some("5".into()),
            genes: vec!["G".into()],
            caller: Some("manta".into()),
            support: None,
            common: false,
        }])),
        step35: Ok(None),
    }
}

/// A P/LP SNV at chr1:2000 in G, inside the deletion's span.
fn snv(gt: &str, copies: usize) -> Variant {
    Variant {
        gene: "G".into(),
        key: ("chr1".into(), 2000, "A".into(), "T".into()),
        chrom: "chr1".into(),
        pos: 2000,
        text: "chr1:2000 A>T".into(),
        genotype: if copies >= 2 { "hom" } else { "het" }.into(),
        copies,
        ploidy: 2,
        gt: gt.into(),
        ps: None,
        allele: None,
    }
}

fn spans() -> Spans {
    IndexMap::from([("G".to_owned(), ("chr1".to_owned(), 500, 6000))])
}

fn pair(del_gt: Option<&str>, v: &Variant) -> og_summary::collect::two_hits::TwoHits {
    let (dels, read, _) = gene_deletions(&annotsv_del(del_gt), &spans());
    assert_eq!(read, ["AnnotSV class 4 and 5 deletions (step 05)"]);
    two_hits(std::slice::from_ref(v), &dels)
}

#[test]
fn a_no_call_deletion_makes_no_homozygous_variant_biallelic() {
    for gt in ["./.", ".|.", ".", "0/.", "./0"] {
        let found = pair(Some(gt), &snv("1/1", 2));
        assert!(
            found.entries.is_empty(),
            "GT {gt}: an uncalled deletion paired: {:?}",
            found.entries
        );
        assert!(found.hemi.is_empty(), "GT {gt}: read as hemizygous");
        assert!(found.pairs.is_empty(), "GT {gt}: a pair was written");
    }
}

#[test]
fn a_deletion_with_no_gt_makes_no_homozygous_variant_biallelic() {
    let found = pair(None, &snv("1/1", 2));
    assert!(found.entries.is_empty(), "{:?}", found.entries);
    assert!(found.hemi.is_empty());
}

#[test]
fn a_reference_genotype_still_pairs_with_nothing() {
    for gt in ["0/0", "0|0", "0"] {
        assert!(pair(Some(gt), &snv("1/1", 2)).entries.is_empty(), "GT {gt}");
    }
}

#[test]
fn a_called_deletion_still_makes_a_homozygous_variant_hemizygous() {
    for gt in ["0/1", "1/0", "./1", "1/1", "0|1"] {
        let v = snv("1/1", 2);
        let found = pair(Some(gt), &v);
        let e = &found.entries["G"];
        assert_eq!(e.verdict, Verdict::Biallelic, "GT {gt}");
        assert_eq!(e.genotype, "hemizygous", "GT {gt}");
        assert!(
            e.why
                .starts_with("the variant is called homozygous inside the deletion"),
            "GT {gt}: {}",
            e.why
        );
        assert!(found.hemi.contains_key(&("G".to_owned(), v.key.clone())));
    }
}

#[test]
fn a_called_deletion_still_pairs_with_a_het_variant() {
    let found = pair(Some("0/1"), &snv("0/1", 1));
    assert_eq!(found.entries["G"].verdict, Verdict::PossiblyBiallelic);
}
