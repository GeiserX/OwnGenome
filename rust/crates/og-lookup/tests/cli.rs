//! `og lookup` on the synthetic gVCF, annotated VCF, gene lines, site table
//! and summary of `tests/parity/inputs/lookup`: what the old
//! `tests/test_genome_lookup.py` checked, the coverage union, and the same
//! answers through a tabix index as from a whole-file read.

use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../tests/parity/inputs/lookup")
}

fn p(name: &str) -> String {
    inputs().join(name).to_string_lossy().into_owned()
}

/// Runs `og lookup ARGS`: exit code, stdout, stderr.
fn og(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let rc = og_lookup::cli(&args, &mut out, &mut err);
    (
        rc,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

fn resolve(extra: &[&str]) -> (u8, String) {
    let mut a = vec!["resolve"];
    a.extend_from_slice(extra);
    let (rc, out, _) = og(&a);
    (rc, out)
}

fn report(
    kind: &str,
    query: &str,
    chrom: &str,
    start: u64,
    end: u64,
    files: &[&str],
) -> (u8, String) {
    let (s, e) = (start.to_string(), end.to_string());
    let mut a = vec![
        "report", "--query", query, "--kind", kind, "--chrom", chrom, "--start", &s, "--end", &e,
    ];
    a.extend_from_slice(files);
    let (rc, out, err) = og(&a);
    assert!(err.is_empty() || rc != 0, "stderr: {err}");
    (rc, out)
}

#[test]
fn resolve_finds_rsids_and_genes_where_the_python_did() {
    let (sites, ann, gtf, traits) = (
        p("sites.tsv.gz"),
        p("ann.vcf"),
        p("gencode.v50.basic.genes.gtf"),
        p("traits.tsv"),
    );
    let (rc, out) = resolve(&["--query", "rs2002", "--sites", &sites, "--annotated", &ann]);
    assert_eq!(rc, 0);
    // The table names rs2002 twice at the same GRCh38 position: one row.
    assert_eq!(
        out,
        "site\tchr7\t150\t150\tA\tthe consumer site table (dbSNP positions)\n"
    );
    let (rc, out) = resolve(&["--query", "RS6006", "--table", &traits]);
    assert_eq!(
        (rc, out.as_str()),
        (0, "site\tchr7\t180\t180\tG\ttraits.tsv\n")
    );
    // rs1001 is in Existing_variation; rs10011 is another rsID.
    let (rc, out) = resolve(&["--query", "rs1001", "--sites", &sites, "--annotated", &ann]);
    assert_eq!(rc, 0);
    assert_eq!(
        out,
        "site\tchr7\t250\t250\tC\tthe annotated VCF (a variant called in this sample)\n"
    );
    let (rc, out) = resolve(&["--query", "rs999", "--sites", &sites, "--annotated", &ann]);
    assert_eq!((rc, out.as_str()), (3, ""));
    let (rc, out) = resolve(&["--query", "genea", "--genes", &gtf]);
    assert_eq!(
        (rc, out.as_str()),
        (0, "gene\tchr7\t90\t700\t.\tGENCODE 50 (GENEA)\n")
    );
    let (rc, out) = resolve(&["--query", "GENEX", "--genes", &gtf]);
    assert_eq!(rc, 0);
    assert_eq!(out.lines().count(), 2, "both pseudoautosomal copies: {out}");
    assert_eq!(resolve(&["--query", "NOPE1", "--genes", &gtf]).0, 3);
    assert_eq!(resolve(&["--query", "GENEA"]).0, 3);
    assert_eq!(resolve(&["--query", "a b"]).0, 2);
}

#[test]
fn a_site_reads_the_gvcf_as_the_python_did() {
    let (g, a, s) = (p("g.vcf"), p("ann.vcf"), p("summary.json"));
    let files = [
        "--gvcf",
        g.as_str(),
        "--annotated",
        a.as_str(),
        "--summary",
        s.as_str(),
    ];
    let site = |q: &str, pos: u64, extra: &[&str]| {
        let mut f = files.to_vec();
        f.extend_from_slice(extra);
        report("site", q, "chr7", pos, pos, &f).1
    };
    let out = site("rs2002", 150, &["--ref", "A"]);
    assert!(
        out.contains(
            "genotype: A/A (homozygous reference, inside the reference block chr7:100-200)"
        ),
        "{out}"
    );
    assert!(
        out.contains("lowest depth in the block 31") && out.contains("genotype quality 50"),
        "{out}"
    );
    assert!(
        out.contains("no variant of this sample here in the annotated VCF"),
        "{out}"
    );
    assert!(
        out.contains("Traits panel (step 42) (stale: from an earlier run):"),
        "{out}"
    );
    assert!(site("chr7:150", 150, &[]).contains("genotype: ref/ref (homozygous reference"));
    let out = site("chr7:350", 350, &[]);
    assert!(
        out.contains("genotype: no call (no gVCF record covers this position"),
        "{out}"
    );
    assert!(
        site("chr7:405", 405, &[])
            .contains("genotype: no call (reference block chr7:400-410 has no genotype)")
    );
    let out = site("rs1001", 250, &[]);
    assert!(
        out.contains("genotype: C/T (heterozygous; REF C, ALT T)"),
        "{out}"
    );
    assert!(
        out.contains("depth 33") && out.contains("QUAL 42.5"),
        "{out}"
    );
    // The most severe entry: missense in GENEA, not GENEZ upstream.
    assert!(
        out.contains("gene: GENEA") && !out.contains("gene: GENEZ"),
        "{out}"
    );
    assert!(
        out.contains("ClinVar (pipeline's ClinVar file): Pathogenic"),
        "{out}"
    );
    assert!(
        out.contains("ClinVar review status: criteria provided, multiple submitters"),
        "{out}"
    );
    assert!(
        out.contains("PharmCAT (step 07):") && !out.contains("MUST-NOT-PRINT"),
        "{out}"
    );
    assert!(site("chr7:600", 600, &[]).contains("filter: RefCall"));
    let out = site("chr7:660", 660, &[]);
    assert!(
        out.contains("genotype: C|G (heterozygous; REF A, ALT C,G)"),
        "{out}"
    );
}

#[test]
fn a_gene_lists_its_variants_and_its_called_bases() {
    let (g, a, s) = (p("g.vcf"), p("ann.vcf"), p("summary.json"));
    let (rc, out) = report(
        "gene",
        "GENEA",
        "chr7",
        90,
        700,
        &[
            "--gvcf",
            &g,
            "--annotated",
            &a,
            "--summary",
            &s,
            "--source",
            "GENCODE 50 (GENEA)",
        ],
    );
    assert_eq!(rc, 0);
    assert!(
        out.contains("chr7:250 C>T C/T (heterozygous; COSV1,rs1001"),
        "{out}"
    );
    assert!(
        !out.contains("chr7:650"),
        "a variant of another gene: {out}"
    );
    // 100-200, 250, 251-299, 420-430, 500-520 (the deletion 500-503 and the
    // block 501-520 count once), 600, 650, 660, 670: 187 of 611.
    assert!(
        out.contains("187 of 611 bases have a genotype call (30.6%)"),
        "{out}"
    );
    let (_, out) = report("gene", "GENEA", "chr7", 90, 700, &["--gvcf", &g]);
    assert!(
        out.contains("chr7:650 T>G G/G (homozygous") && out.contains("filter RefCall"),
        "{out}"
    );
    assert!(
        !out.contains("chr7:100 "),
        "a block is not a variant: {out}"
    );
    assert!(out.contains("no summary.json for this sample"), "{out}");
}

#[test]
fn coverage_counts_overlapping_records_once() {
    // A 10-base gene, a 5-base deletion record at 100 and a reference block
    // over 101-103: 100-104 and 106-109 have a call, 9 bases. The Python
    // summed the records and said 10 (12 clamped to the span).
    let d1 = p("d1.g.vcf");
    let (rc, out) = report("gene", "GENED", "chr7", 100, 109, &["--gvcf", &d1]);
    assert_eq!(rc, 0);
    assert!(
        out.contains("  gVCF: 9 of 10 bases have a genotype call (90.0%);"),
        "{out}"
    );
}

#[test]
fn the_index_and_a_whole_file_read_give_the_same_answer() {
    for (q, kind, chrom, s, e) in [
        ("rs1001", "site", "chr7", 250, 250),
        ("chr7:502", "site", "chr7", 502, 502),
        ("GENEA", "gene", "chr7", 90, 700),
        ("GENEX", "gene", "chrX", 900, 1100),
        ("MT-RNR1", "gene", "chrM", 70, 1000),
    ] {
        let plain = report(
            kind,
            q,
            chrom,
            s,
            e,
            &["--gvcf", &p("g.vcf"), "--annotated", &p("ann.vcf")],
        );
        let tbi = report(
            kind,
            q,
            chrom,
            s,
            e,
            &["--gvcf", &p("g.vcf.gz"), "--annotated", &p("ann.vcf.gz")],
        );
        let csi = report(
            kind,
            q,
            chrom,
            s,
            e,
            &[
                "--gvcf",
                &p("csi.g.vcf.gz"),
                "--annotated",
                &p("ann.vcf.gz"),
            ],
        );
        assert_eq!(plain, tbi, "{q}");
        assert_eq!(plain, csi, "{q}");
    }
}

#[test]
fn query_answers_every_query_and_says_how_it_read() {
    let base = [
        "query",
        "--sites",
        &p("sites.tsv.gz"),
        "--genes",
        &p("gencode.v50.basic.genes.gtf"),
        "--table",
        &p("traits.tsv"),
    ]
    .map(String::from);
    let run = |gvcf: &str, ann: &str, qs: &[&str]| {
        let mut a: Vec<&str> = base.iter().map(String::as_str).collect();
        a.extend(["--gvcf", gvcf, "--annotated", ann, "--"]);
        a.extend_from_slice(qs);
        og(&a)
    };
    let (rc, out, err) = run(
        &p("g.vcf.gz"),
        &p("ann.vcf.gz"),
        &["rs2002", "rs999", "GENEZ"],
    );
    assert_eq!(rc, 3, "the last failing query's code");
    assert!(
        err.contains("ERROR: rs999 was not found") && !err.contains("NOTE"),
        "{err}"
    );
    let blocks: Vec<&str> = out.split("\n\n").collect();
    // rs2002, a blank line for rs999, GENEZ.
    assert!(blocks[0].starts_with("== rs2002: chr7:150 "), "{out}");
    assert!(out.contains("\n\n\n== GENEZ: chr7:640-660 "), "{out}");
    assert!(out.ends_with("\n\n"), "{out}");
    let (rc, _, err) = run(&p("g.vcf"), &p("ann.vcf"), &["chr7:250"]);
    assert_eq!(rc, 0);
    assert_eq!(
        err.matches("has no .tbi index; reading the whole file")
            .count(),
        2,
        "{err}"
    );
    let (rc, _, err) = og(&["query", "--gvcf", &p("g.vcf")]);
    assert_eq!(rc, 2, "{err}");
}
