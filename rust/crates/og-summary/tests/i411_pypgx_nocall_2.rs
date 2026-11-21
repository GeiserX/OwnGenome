//! summary.json's pypgx block counts a gene as called only when its
//! diplotype is a call by the rule steps 27 and 36 use
//! (`og_steps::pgx_rules::is_call`): a diplotype with one `Indeterminate`
//! allele, such as `*1/Indeterminate`, in any case, is no call, and so is
//! one that names no two alleles or several diplotypes. A real call still
//! counts (issue #411, pgx-06).

use og_summary::collect::pgx::Pypgx;
use pycompat::json::Value;
use std::fs;
use std::path::PathBuf;

const HEADER: &str = "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\tStatus\tDefining_variants\tMin_GQ\tWeakest_variant\tDefining_filter\tQuality\n";

/// The pypgx summary of sample S with one row per (gene, diplotype), read
/// from a folder of its own.
fn read(tag: &str, rows: &[(&str, &str)]) -> Pypgx {
    let d: PathBuf = std::env::temp_dir().join(format!(
        "og-i411-pypgx-nocall-2-{}-{tag}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    let mut body = HEADER.to_owned();
    for (gene, dip) in rows {
        body.push_str(&format!("{gene}\t{dip}\t\t\tVCF\tok\t\t\t\t\tok\n"));
    }
    let path = d.join("S_pypgx_summary.tsv");
    fs::write(&path, body).unwrap();
    let got = Pypgx::read(path.to_str().unwrap(), "S").unwrap();
    fs::remove_dir_all(&d).unwrap();
    got
}

fn called_in_json(p: &Pypgx) -> Value {
    let Value::Object(m) = p.to_json() else {
        panic!("the pypgx block is not an object")
    };
    m["genes_called"].clone()
}

#[test]
fn a_half_indeterminate_diplotype_is_not_a_called_gene() {
    let p = read(
        "half",
        &[
            ("CYP2D6", "*1/Indeterminate"),
            ("CYP2C9", "Indeterminate/*2"),
            ("CYP2C19", "*1/indeterminate"),
            ("UGT1A1", " *1/INDETERMINATE "),
            ("CYP3A5", "Indeterminate/Indeterminate"),
        ],
    );
    assert_eq!(p.genes_total, 5);
    assert_eq!(p.genes_called, 0);
    assert_eq!(called_in_json(&p), Value::Int(0));
}

#[test]
fn a_real_call_still_counts_beside_the_no_calls() {
    let p = read(
        "mixed",
        &[
            ("CYP2D6", "*1/Indeterminate"),
            ("CYP2C19", "*1/*2"),
            ("CYP2B6", " *1/*6 "),
            ("UGT1A1", "Indeterminate"),
            ("CYP2A6", "FAILED"),
            ("CYP2E1", "."),
            ("NUDT15", ""),
        ],
    );
    assert_eq!(p.genes_total, 7);
    assert_eq!(p.genes_called, 2);
    assert_eq!(called_in_json(&p), Value::Int(2));
    // the CYP2D6 diplotype itself is still shown as pypgx wrote it
    assert_eq!(p.cyp2d6, "*1/Indeterminate");
}

#[test]
fn a_diplotype_the_steps_refuse_is_not_counted() {
    // one allele, several diplotypes, and a missing allele are no call in
    // steps 27 and 36, so the summary does not count them either
    let p = read(
        "refused",
        &[
            ("CYP2D6", "*1"),
            ("CYP2C9", "*1/*2;*1/*3"),
            ("CYP2C19", "*1/"),
            ("TPMT", "*1/*3A"),
        ],
    );
    assert_eq!(p.genes_total, 4);
    assert_eq!(p.genes_called, 1);
}
