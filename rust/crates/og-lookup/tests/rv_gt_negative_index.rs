//! A genotype index below zero names no allele. VCF writes GT as allele
//! numbers from 0 (REF) up, or `.`; a reader that indexed `[REF] + ALT` as a
//! Python list read `-1` as the last ALT, so `GT -1/-1` at REF A, ALT T
//! answered `T/T`. The lookup reads it as a no-call, as it already read an
//! index past the last allele.

use std::fs;
use std::path::{Path, PathBuf};

fn input(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/parity/inputs/lookup")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

/// Runs `og lookup query` for each of `queries` over `gvcf`: stdout.
fn lookup(gvcf: &str, queries: &[&str]) -> String {
    let mut args: Vec<String> = [
        "query",
        "--sites",
        &input("sites.tsv.gz"),
        "--genes",
        &input("gencode.v50.basic.genes.gtf"),
        "--table",
        &input("traits.tsv"),
        "--gvcf",
        gvcf,
        "--",
    ]
    .map(String::from)
    .to_vec();
    args.extend(queries.iter().map(|q| (*q).to_owned()));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let rc = og_lookup::cli(&args, &mut out, &mut err);
    let (out, err) = (
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    );
    assert_eq!(rc, 0, "{out}{err}");
    out
}

#[test]
fn a_negative_genotype_index_is_no_call() {
    let dir: PathBuf = Path::new(env!("CARGO_TARGET_TMPDIR")).join("rv_gt_negative_index");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let gvcf = dir.join("g.vcf");
    let body = [
        "chr20\t100\t.\tA\tT\t30\tPASS\t.\tGT:GQ\t-1/-1:30",
        "chr20\t200\t.\tA\tT\t30\tPASS\t.\tGT:GQ\t1|-1:30",
        "chr20\t300\t.\tA\tT,<*>\t30\tPASS\t.\tGT:GQ\t0/1:30",
        "chr20\t400\t.\tA\tT,<*>\t30\tPASS\t.\tGT:GQ\t1/1:30",
    ];
    let text = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n{}\n",
        body.join("\n")
    );
    fs::write(&gvcf, text).unwrap();
    let gvcf = gvcf.to_string_lossy().into_owned();

    let genotype = |q: &str| -> String {
        let out = lookup(&gvcf, &[q]);
        out.lines()
            .find(|l| l.trim_start().starts_with("genotype: "))
            .unwrap_or_else(|| panic!("no genotype line for {q}: {out}"))
            .trim_start()
            .to_owned()
    };

    // -1 is no allele: a no-call, never the ALT T.
    assert_eq!(
        genotype("chr20:100"),
        "genotype: ./. (no call; REF A, ALT T)"
    );
    assert_eq!(
        genotype("chr20:200"),
        "genotype: T|. (half call; REF A, ALT T)"
    );
    // The answers of readable genotypes do not change.
    assert_eq!(
        genotype("chr20:300"),
        "genotype: A/T (heterozygous; REF A, ALT T)"
    );
    assert_eq!(
        genotype("chr20:400"),
        "genotype: T/T (homozygous; REF A, ALT T)"
    );

    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn gene_coverage_requires_an_index_that_names_an_allele() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("rv_gt_coverage_lookup");
    fs::create_dir_all(&dir).unwrap();
    let gvcf = dir.join("g.vcf");
    fs::write(
        &gvcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
         chr7\t100\t.\tA\tT\t30\tPASS\tEND=109\tGT\t-1/-1\n\
         chr7\t200\t.\tA\tT\t30\tPASS\tEND=219\tGT\t2/2\n\
         chr7\t300\t.\tA\t.\t30\tPASS\tEND=302\tGT\t1/1\n\
         chr7\t400\t.\tA\tT,G\t30\tPASS\tEND=404\tGT\t2|-1\n\
         chr7\t500\t.\tA\t<*>\t30\tPASS\tEND=506\tGT\t0/.\n",
    )
    .unwrap();
    let out = lookup(gvcf.to_str().unwrap(), &["GENEA"]);

    assert!(out.contains("12 of 611 bases have a genotype call"), "{out}");
    fs::remove_dir_all(&dir).unwrap();
}
