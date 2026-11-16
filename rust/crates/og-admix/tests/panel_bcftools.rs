//! The panel reader against what bcftools prints for the same files
//! (`tests/data/README.md` holds the commands): every biallelic SNV, in
//! file order, with the chosen people's phased genotypes in the order
//! asked for.

use og_admix::panel::Haplotypes;
use std::path::{Path, PathBuf};

fn here(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

/// The panel read back as bcftools query's lines.
fn query(path: &Path, people: &[&str]) -> String {
    let h = Haplotypes::read(path, people, None).unwrap();
    let mut out = String::new();
    for (m, s) in h.sites().iter().enumerate() {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}",
            h.chrom(),
            s.pos,
            s.ref_ as char,
            s.alt as char
        ));
        for i in 0..people.len() {
            out.push_str(&format!(
                "\t{}|{}",
                h.allele(m, 2 * i),
                h.allele(m, 2 * i + 1)
            ));
        }
        out.push('\n');
    }
    out
}

/// The red control reads the people in another order than bcftools was
/// asked for, so the comparison must fail.
fn asked<'a>(people: &[&'a str]) -> Vec<&'a str> {
    let mut p = people.to_vec();
    if cfg!(feature = "planted-red") {
        p.reverse();
    }
    p
}

#[test]
fn bgzipped_panel_matches_bcftools() {
    let want = std::fs::read_to_string(here("tests/data/panel.P4-P1-P3.bcftools.tsv")).unwrap();
    assert_eq!(want.lines().count(), 6, "the golden itself");
    let got = query(
        &here("tests/data/panel.vcf.gz"),
        &asked(&["P4", "P1", "P3"]),
    );
    assert_eq!(got, want);
}

#[test]
fn parity_family_fixture_matches_bcftools() {
    let want =
        std::fs::read_to_string(here("tests/data/family-phased.mum-dad.bcftools.tsv")).unwrap();
    assert_eq!(want.lines().count(), 3, "the golden itself");
    let got = query(
        &here("../../../tests/parity/inputs/family/phased.vcf"),
        &asked(&["mum", "dad"]),
    );
    assert_eq!(got, want);
}

#[test]
fn skipped_records_are_counted() {
    let path = here("tests/data/panel.vcf.gz");
    let mut r = og_admix::panel::PanelReader::open(&path, &["P2"]).unwrap();
    let mut row = vec![0u64; 1];
    let mut n = 0;
    while r.next_site(&mut row).unwrap().is_some() {
        n += 1;
    }
    // A multi-allelic site, an indel, no ALT and a `*` allele.
    assert_eq!((n, r.skipped()), (6, 4));
}
