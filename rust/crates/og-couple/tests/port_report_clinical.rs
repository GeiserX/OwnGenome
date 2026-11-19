//! The couple screen's set-aside alleles, on synthetic ClinVar hits at real
//! gene positions: a HEXA pseudodeficiency allele and the RBM8A
//! low-expression allele (rs139428292) are never a carrier finding on their
//! own.
//!
//! These were the couple screen checks of two Python tests whose other
//! checks are in `rust/crates/og-render/tests/clinical.rs`:
//! check 7 of `tests/test_acmg_sf_rules.py` and check 9 of
//! `tests/test_cross_engine_biallelic.py`.

mod common;

use common::{Dir, og, read, repo, rows};

const VCF_HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
const PEDIGREE: &str =
    "sample\tsex\tfather\tmother\tpartner\nA\tfemale\t.\t.\tB\nB\tmale\t.\t.\tA\n";

/// A ClinVar P hit, reviewed by several submitters.
fn hit(chrom: &str, pos: i64, r: &str, a: &str, gene: &str, gt: &str) -> String {
    format!(
        "{chrom}\t{pos}\t.\t{r}\t{a}\t50\tPASS\tGENEINFO={gene}:1;CLNSIG=Pathogenic;\
         CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\tGT:GQ:DP\t{gt}:50:30\n"
    )
}

/// test_acmg_sf_rules.py check 7, the couple screen part.
#[test]
fn a_hexa_pseudodeficiency_allele_is_set_aside_and_the_screen_goes_on() {
    let t = Dir::new("clinical-hexa");
    let hits = t.put(
        "S_clinvar_hits.vcf",
        &format!(
            "{VCF_HEAD}{}{}{}",
            hit("chr7", 117559590, "ATCT", "A", "CFTR", "0/1"),
            hit("chr15", 72350584, "G", "A", "HEXA", "0/1"),
            hit("chr15", 72346579, "G", "GGATA", "HEXA", "0/1")
        ),
    );
    let ped = t.put("pedigree.tsv", PEDIGREE);
    let genes = repo("assets/carrier_genes.tsv");
    let out = t.path("out");
    let (a, b) = (format!("A={hits}"), format!("B={hits}"));
    let (rc, log) = og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &genes,
            "--clinvar",
            &a,
            "--clinvar",
            &b,
            "--out-dir",
            &out,
        ],
    );
    let md = read(&format!("{out}/couple_A_B.md"));
    let mut found: Vec<String> = rows(&format!("{out}/couple_A_B.tsv"))
        .iter()
        .map(|r| r["gene"].clone())
        .collect();
    found.sort();
    let aside: Vec<&str> = md
        .lines()
        .filter(|x| x.contains("Set aside for A"))
        .collect();
    assert_eq!(found, ["CFTR", "HEXA"], "{log}");
    assert_eq!(aside.len(), 1, "{md}");
    assert!(aside[0].contains("pseudodeficiency"), "{}", aside[0]);
    assert_eq!(rc, 0, "{log}");
    assert!(md.contains("Set aside for A, not a carrier finding on its own"));
}

// RBM8A c.-21G>A (rs139428292), on the plus strand C>T.
const TAR: (&str, i64, &str, &str) = ("chr1", 145927447, "C", "T");

/// test_cross_engine_biallelic.py check 9.
#[test]
fn the_rbm8a_allele_counts_only_against_a_deletion_in_the_partner() {
    let t = Dir::new("clinical-tar");
    let tar = t.put(
        "tar.vcf",
        &format!(
            "{VCF_HEAD}{}",
            hit(TAR.0, TAR.1, TAR.2, TAR.3, "RBM8A", "0/1")
        ),
    );
    let empty = t.put("empty.vcf", VCF_HEAD);
    let ped = t.put("pedigree.tsv", PEDIGREE);
    // B's null allele: a 300 kb Delly deletion over RBM8A's exons (no RDCN,
    // so read depth is no question).
    let null_sv = t.put(
        "B_sv.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
         chr1\t145700000\tDEL0\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=146000000\tGT\t0/1\n",
    );
    let exons: String = (145921000..145928000)
        .step_by(1200)
        .map(|s| format!("chr1\t{s}\t{}\n", s + 100))
        .collect();
    let exons = t.put("exons.bed", &exons);
    let genes = repo("assets/carrier_genes.tsv");
    let a = format!("A={tar}");
    let couple = |name: &str, b_hits: &str, extra: &[&str]| {
        let out = t.path(name);
        let b = format!("B={b_hits}");
        let mut args = vec![
            "--pedigree",
            &ped,
            "--genes",
            &genes,
            "--clinvar",
            &a,
            "--clinvar",
            &b,
        ];
        args.extend_from_slice(extra);
        args.extend_from_slice(&["--out-dir", &out]);
        let (rc, log) = og("screen", &args);
        let (rows, md) = if rc == 0 {
            (
                rows(&format!("{out}/couple_A_B.tsv")),
                read(&format!("{out}/couple_A_B.md")),
            )
        } else {
            (Vec::new(), String::new())
        };
        (rc, log, rows, md)
    };

    let delly = format!("B={null_sv}");
    let (_, log, got, _) = couple("null", &empty, &["--delly", &delly, "--exons", &exons]);
    let x: Vec<_> = got.iter().filter(|r| r["gene"] == "RBM8A").collect();
    assert_eq!(x.len(), 1, "{log}");
    assert_eq!(x[0]["result"], "shared");
    assert_eq!(x[0]["offspring_risk"], "1 in 4 per pregnancy");
    assert!(x[0]["findings_a"].contains("rs139428292"));
    assert!(x[0]["explanation"].contains("hypomorph"));

    let (rc, log, got, _) = couple("both", &tar, &[]);
    assert_eq!(rc, 0, "{log}");
    assert!(
        got.iter().all(|r| r["gene"] != "RBM8A"),
        "the allele in both partners is set aside"
    );
    let (rc, log, got, md) = couple("one", &empty, &[]);
    assert_eq!(rc, 0, "{log}");
    assert!(got.iter().all(|r| r["gene"] != "RBM8A"), "the allele alone");
    assert!(
        md.contains("| A | 0 / 1 |") && md.contains("Set aside for A"),
        "no finding of the listed genes (0 / 1), and a gap line: {md}"
    );
}
