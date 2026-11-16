//! The disease-gene drug warnings of `og pgx cpic-report` against the shipped
//! `assets/pgx_contraindications.tsv` (issue #411, pgx-21, pgx-22, pgx-23):
//!
//! - one heterozygous pathogenic variant in a recessive urea-cycle gene
//!   (ASL, ASS1, CPS1, NAGS) makes a carrier: a NOTE, never a WARNING;
//! - both copies changed, or two variants in the gene, give the WARNING;
//! - one POLG variant still warns (a known valproate risk);
//! - SCN1A gives nothing: no label or guideline names it;
//! - the CPIC level is printed as the provisional level it is.

use std::collections::HashSet;

use og_pgx::sections::contraindication_lines;
use og_pgx::tables::{ClinvarHit, Contraindication, read_contraindications};

fn table() -> Vec<Contraindication> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../../assets/pgx_contraindications.tsv"
    );
    read_contraindications(path).unwrap()
}

fn hit(gene: &str, pos: &str, alt: &str, genotype: &str) -> ClinvarHit {
    ClinvarHit {
        chrom: "chr13".to_owned(),
        pos: pos.to_owned(),
        ref_: "G".to_owned(),
        alt: alt.to_owned(),
        genotype: genotype.to_owned(),
        geneinfo: Some(format!("{gene}:1")),
        clnsig: Some("Pathogenic".to_owned()),
    }
}

fn lines(hits: &[ClinvarHit]) -> Vec<String> {
    let meds: HashSet<String> = ["valproate".to_owned()].into();
    contraindication_lines(hits, &table(), &meds)
}

fn warnings(out: &[String]) -> Vec<&str> {
    out.iter()
        .map(|l| l.trim())
        .filter(|l| l.starts_with("WARNING:"))
        .collect()
}

#[test]
fn one_heterozygous_urea_cycle_hit_is_a_note() {
    for gene in ["ASL", "ASS1", "CPS1", "NAGS"] {
        let out = lines(&[hit(gene, "40000000", "A", "0/1")]);
        assert!(warnings(&out).is_empty(), "{gene}: {out:#?}");
        let note = format!(
            "  NOTE: valproic acid -- {gene} chr13:40000000 G>A, genotype 0/1, ClinVar Pathogenic; ON THE MEDICATION LIST."
        );
        assert_eq!(out.first(), Some(&note), "{gene}: {out:#?}");
        assert!(
            out.iter().any(|l| l.contains("a carrier, not a patient")),
            "{gene}: {out:#?}"
        );
    }
    // A phased genotype is still one changed copy.
    for gt in ["0|1", "1|0"] {
        let out = lines(&[hit("ASS1", "40000000", "A", gt)]);
        assert!(warnings(&out).is_empty(), "{gt}: {out:#?}");
    }
}

#[test]
fn both_copies_or_two_variants_in_a_urea_cycle_gene_warn() {
    for gt in ["1/1", "1|1", "1"] {
        let out = lines(&[hit("ASS1", "40000000", "A", gt)]);
        assert_eq!(warnings(&out).len(), 1, "{gt}: {out:#?}");
        assert!(
            warnings(&out)[0].contains("both copies changed"),
            "{gt}: {out:#?}"
        );
    }
    // Step 06 splits a 1/2 record into 1/. and ./1: the missing allele is
    // the other ALT, so both copies are changed.
    for gt in ["1/.", "./1"] {
        let out = lines(&[hit("ASS1", "40000000", "A", gt)]);
        assert_eq!(
            warnings(&out),
            [format!(
                "WARNING: valproic acid -- ASS1 chr13:40000000 G>A, genotype {gt}, ClinVar Pathogenic, both copies changed in this recessive gene (the other copy is not the reference); ON THE MEDICATION LIST."
            )],
            "{gt}: {out:#?}"
        );
        assert!(!out.iter().any(|l| l.contains("NOTE:")), "{gt}: {out:#?}");
    }
    // Two variants in trans, at two positions or as the two alleles of one
    // split multi-allelic site: one WARNING for each.
    for (pos2, alt2) in [("40000500", "A"), ("40000000", "C")] {
        let out = lines(&[
            hit("CPS1", "40000000", "A", "0/1"),
            hit("CPS1", pos2, alt2, "1/0"),
        ]);
        let w = warnings(&out);
        assert_eq!(w.len(), 2, "{pos2}>{alt2}: {out:#?}");
        assert!(w.iter().all(|l| l.contains("2 variants in CPS1")), "{w:#?}");
        assert!(!out.iter().any(|l| l.contains("NOTE:")), "{out:#?}");
    }
    // The same variant listed twice (two ClinVar records) is one variant.
    let out = lines(&[
        hit("NAGS", "40000000", "A", "0/1"),
        hit("NAGS", "40000000", "A", "0/1"),
    ]);
    assert!(warnings(&out).is_empty(), "{out:#?}");
}

#[test]
fn one_polg_or_otc_variant_still_warns() {
    for gene in ["POLG", "OTC"] {
        let out = lines(&[hit(gene, "89300000", "T", "0/1")]);
        assert_eq!(
            warnings(&out),
            [format!(
                "WARNING: valproic acid -- {gene} chr13:89300000 G>T, genotype 0/1, ClinVar Pathogenic; ON THE MEDICATION LIST."
            )],
            "{out:#?}"
        );
    }
}

#[test]
fn an_scn1a_hit_gives_no_drug_line() {
    assert!(!table().iter().any(|r| r.gene == "SCN1A"));
    let out = lines(&[hit("SCN1A", "166000000", "T", "0/1")]);
    assert!(out.is_empty(), "{out:#?}");
}

#[test]
fn every_row_quotes_a_label() {
    let bad: Vec<String> = table()
        .into_iter()
        .filter(|r| !r.label_source.starts_with("https://dailymed.nlm.nih.gov/"))
        .map(|r| format!("{} {}", r.gene, r.drug))
        .collect();
    assert!(bad.is_empty(), "rows no label names: {bad:?}");
}

#[test]
fn the_cpic_level_reads_as_provisional() {
    let out = lines(&[hit("POLG", "89300000", "T", "0/1")]);
    assert!(
        out.iter().any(
            |l| l.ends_with("(provisional CPIC level A/B: no CPIC guideline covers this pair)")
        ),
        "{out:#?}"
    );
    assert!(
        !out.iter().any(|l| l.ends_with("(CPIC level A/B)")),
        "{out:#?}"
    );
}
