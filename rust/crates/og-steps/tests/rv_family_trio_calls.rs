//! `og family report` on small hand-made trios: the calls the trio report
//! missed or trusted too much.
//!
//! 1. Compound heterozygotes: a child call with two different ALTs (1/2,
//!    as GLnexus writes an unsplit multiallelic record) is heterozygous, and
//!    each ALT that is a candidate is paired with its origin, so a 1/2 from a
//!    0/1 father and a 0/2 mother is one pair in trans. The two ALTs of one
//!    1/2 call are in trans even when the parents do not show where each
//!    came from. A 0/1 child gives the same pairs as before.
//! 2. A son's sex chromosomes outside the PARs: his Y comes from his father
//!    and his X from his mother, so a de novo call there needs only that
//!    parent's call, and the depth and quality floors apply only to him and
//!    that parent.
//! 3. Parental reads: a parent called 0/0 with reads on the new allele (a
//!    mosaic parent) keeps a de novo candidate out of `pass`.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

fn par_bed() -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../assets/par_grch38.bed")
        .display()
        .to_string()
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-family-trio-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// A trio: `kid` (sex `kid_sex`), father `dad`, mother `mum`, all three
/// sequenced.
fn members(d: &Path, kid_sex: &str) -> String {
    let p = d.join("members.tsv");
    fs::write(
        &p,
        format!(
            "sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced\n\
             kid\t{kid_sex}\tdad\tmum\t\tyes\tyes\n\
             dad\tmale\t\t\tmum\tyes\tyes\n\
             mum\tfemale\t\t\tdad\tyes\tyes\n"
        ),
    )
    .unwrap();
    p.display().to_string()
}

/// A joint VCF with the columns kid, dad, mum and FORMAT GT:AD:DP:GQ; each
/// record is `(chrom, pos, ref, alt, [kid, dad, mum])`.
fn joint_vcf(d: &Path, records: &[(&str, i64, &str, &str, [&str; 3])]) -> String {
    let mut out = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n\
         ##contig=<ID=chr2,length=242193529>\n##contig=<ID=chrX,length=156040895>\n\
         ##contig=<ID=chrY,length=57227415>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n\
         ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkid\tdad\tmum\n",
    );
    for (chrom, pos, r, a, s) in records {
        out += &format!(
            "{chrom}\t{pos}\t.\t{r}\t{a}\t50\tPASS\t.\tGT:AD:DP:GQ\t{}\t{}\t{}\n",
            s[0], s[1], s[2]
        );
    }
    let p = d.join("family.vcf");
    fs::write(&p, out).unwrap();
    p.display().to_string()
}

/// The child's candidate VCF as step 23 leaves it: split per ALT, the gene
/// in INFO/GENE. Each record is `(chrom, pos, ref, alt, gene)`.
fn candidates(d: &Path, records: &[(&str, i64, &str, &str, &str)]) -> String {
    let mut out =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (chrom, pos, r, a, g) in records {
        out += &format!("{chrom}\t{pos}\t.\t{r}\t{a}\t50\tPASS\tGENE={g}\n");
    }
    let p = d.join("kid.candidates.vcf");
    fs::write(&p, out).unwrap();
    p.display().to_string()
}

/// Runs `og family report` in `d` and returns the output folder.
fn report(d: &Path, members: &str, vcf: &str, cands: Option<&str>) -> PathBuf {
    let rel = d.join("relationships.tsv");
    fs::write(
        &rel,
        "sample_a\tsample_b\tdeclared\texpected_relatedness\trelatedness\tibs0_rate\tsites\tverdict\treason\n",
    )
    .unwrap();
    let out = d.join("out");
    let par = par_bed();
    let cand = cands.map(|c| format!("kid={c}"));
    let mut args = vec![
        "report",
        "--members",
        members,
        "--relationships",
        rel.to_str().unwrap(),
        "--vcf",
        vcf,
        "--par-bed",
        &par,
        "--out-dir",
        out.to_str().unwrap(),
    ];
    if let Some(c) = &cand {
        args.extend(["--candidates", c.as_str()]);
    }
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    assert_eq!(og_steps::family::cli(&args), 0, "og family report exits 0");
    out
}

/// The rows of a TSV as maps from the header.
fn read_tsv(path: &Path) -> Vec<HashMap<String, String>> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .zip(l.split('\t'))
                .map(|(h, v)| ((*h).to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}

fn pairs(out: &Path) -> Vec<Vec<String>> {
    read_tsv(&out.join("compound_het.tsv"))
        .into_iter()
        .map(|r| {
            [
                "child",
                "gene",
                "variant_a",
                "variant_b",
                "a_from",
                "b_from",
                "phase",
            ]
            .iter()
            .map(|k| r[*k].clone())
            .collect()
        })
        .collect()
}

/// The de novo rows at `chrom:pos`, as (child_gt, quality).
fn de_novo_at(out: &Path, chrom: &str, pos: i64) -> Vec<(String, String)> {
    read_tsv(&out.join("de_novo.tsv"))
        .into_iter()
        .filter(|r| r["chrom"] == chrom && r["pos"] == pos.to_string())
        .map(|r| (r["child_gt"].clone(), r["quality"].clone()))
        .collect()
}

const GOOD_REF: &str = "0/0:30,0:30:99";

// ---- 1. compound heterozygotes ------------------------------------------

#[test]
fn a_two_alt_call_from_both_parents_is_one_pair_in_trans() {
    let d = scratch("two-alt");
    let m = members(&d, "female");
    let vcf = joint_vcf(
        &d,
        &[(
            "chr1",
            1_000_000,
            "C",
            "A,G",
            [
                "1/2:0,15,15:30:99",
                "0/1:15,15,0:30:99",
                "0/2:15,0,15:30:99",
            ],
        )],
    );
    let c = candidates(
        &d,
        &[
            ("chr1", 1_000_000, "C", "A", "G1"),
            ("chr1", 1_000_000, "C", "G", "G1"),
        ],
    );
    let out = report(&d, &m, &vcf, Some(&c));
    assert_eq!(
        pairs(&out),
        vec![vec![
            "kid".to_owned(),
            "G1".into(),
            "chr1:1000000:C>A".into(),
            "chr1:1000000:C>G".into(),
            "father".into(),
            "mother".into(),
            "in trans".into(),
        ]],
        "a 1/2 child of a 0/1 father and a 0/2 mother: one pair in trans"
    );
    let md = fs::read_to_string(out.join("family_report.md")).unwrap();
    assert!(
        md.contains(
            "in one gene: 1; in trans (one from each parent, both copies of the gene hit) 1"
        ),
        "the report counts the pair: {md}"
    );
}

#[test]
fn a_two_alt_call_is_in_trans_even_when_the_parents_do_not_tell() {
    // both parents 1/2: the origin of each ALT is unknown, but the two ALTs
    // of one diploid call are on different copies by definition
    let d = scratch("two-alt-unknown-origin");
    let m = members(&d, "female");
    let vcf = joint_vcf(
        &d,
        &[(
            "chr1",
            1_000_000,
            "C",
            "A,G",
            [
                "1/2:0,15,15:30:99",
                "1/2:0,15,15:30:99",
                "1/2:0,15,15:30:99",
            ],
        )],
    );
    let c = candidates(
        &d,
        &[
            ("chr1", 1_000_000, "C", "A", "G1"),
            ("chr1", 1_000_000, "C", "G", "G1"),
        ],
    );
    let out = report(&d, &m, &vcf, Some(&c));
    assert_eq!(
        pairs(&out),
        vec![vec![
            "kid".to_owned(),
            "G1".into(),
            "chr1:1000000:C>A".into(),
            "chr1:1000000:C>G".into(),
            ".".into(),
            ".".into(),
            "in trans".into(),
        ]],
        "a 1/2 child of two 1/2 parents: one pair in trans, origins unknown"
    );
}

#[test]
fn a_two_alt_call_pairs_with_a_ref_alt_call_in_the_gene() {
    let d = scratch("two-alt-plus-one");
    let m = members(&d, "male");
    let vcf = joint_vcf(
        &d,
        &[
            (
                "chr1",
                1_000_000,
                "C",
                "A,G",
                [
                    "1/2:0,15,15:30:99",
                    "0/1:15,15,0:30:99",
                    "0/2:15,0,15:30:99",
                ],
            ),
            (
                "chr1",
                1_000_500,
                "T",
                "C",
                ["0/1:15,15:30:99", "0/1:15,15:30:99", GOOD_REF],
            ),
        ],
    );
    // only the G of the 1/2 is a candidate, with the 0/1 site from the father
    let c = candidates(
        &d,
        &[
            ("chr1", 1_000_000, "C", "G", "G1"),
            ("chr1", 1_000_500, "T", "C", "G1"),
        ],
    );
    let out = report(&d, &m, &vcf, Some(&c));
    assert_eq!(
        pairs(&out),
        vec![vec![
            "kid".to_owned(),
            "G1".into(),
            "chr1:1000000:C>G".into(),
            "chr1:1000500:T>C".into(),
            "mother".into(),
            "father".into(),
            "in trans".into(),
        ]],
        "the G of a 1/2 (from the mother) and a 0/1 from the father: one pair in trans"
    );
}

#[test]
fn ref_alt_calls_pair_as_before() {
    let d = scratch("ref-alt");
    let m = members(&d, "female");
    let vcf = joint_vcf(
        &d,
        &[
            // from the father
            (
                "chr1",
                2_000_000,
                "A",
                "G",
                ["0/1:15,15:30:99", "0/1:15,15:30:99", GOOD_REF],
            ),
            // from the mother, written ALT first
            (
                "chr1",
                2_000_100,
                "C",
                "T",
                ["1/0:15,15:30:99", GOOD_REF, "0/1:15,15:30:99"],
            ),
            // from the father again
            (
                "chr1",
                2_000_200,
                "G",
                "A",
                ["0/1:15,15:30:99", "1/1:0,30:30:99", GOOD_REF],
            ),
            // both parents het: origin unknown
            (
                "chr2",
                3_000_000,
                "T",
                "G",
                ["0/1:15,15:30:99", "0/1:15,15:30:99", "0/1:15,15:30:99"],
            ),
        ],
    );
    let c = candidates(
        &d,
        &[
            ("chr1", 2_000_000, "A", "G", "G2"),
            ("chr1", 2_000_100, "C", "T", "G2"),
            ("chr1", 2_000_200, "G", "A", "G2"),
            ("chr2", 3_000_000, "T", "G", "G2"),
        ],
    );
    let out = report(&d, &m, &vcf, Some(&c));
    let row = |a: &str, b: &str, fa: &str, fb: &str, ph: &str| {
        ["kid", "G2", a, b, fa, fb, ph]
            .iter()
            .map(|s| (*s).to_owned())
            .collect::<Vec<_>>()
    };
    assert_eq!(
        pairs(&out),
        vec![
            row(
                "chr1:2000000:A>G",
                "chr1:2000100:C>T",
                "father",
                "mother",
                "in trans"
            ),
            row(
                "chr1:2000000:A>G",
                "chr1:2000200:G>A",
                "father",
                "father",
                "in cis"
            ),
            row(
                "chr1:2000000:A>G",
                "chr2:3000000:T>G",
                "father",
                ".",
                "unknown"
            ),
            row(
                "chr1:2000100:C>T",
                "chr1:2000200:G>A",
                "mother",
                "father",
                "in trans"
            ),
            row(
                "chr1:2000100:C>T",
                "chr2:3000000:T>G",
                "mother",
                ".",
                "unknown"
            ),
            row(
                "chr1:2000200:G>A",
                "chr2:3000000:T>G",
                "father",
                ".",
                "unknown"
            ),
        ],
        "0/1 children pair as before"
    );
}

// ---- 2. a son's X and Y outside the PARs -----------------------------------

#[test]
fn a_sons_y_de_novo_needs_only_the_father() {
    let d = scratch("son-y");
    let m = members(&d, "male");
    let vcf = joint_vcf(
        &d,
        &[
            // the mother has no call on Y, as always
            (
                "chrY",
                15_000_000,
                "A",
                "T",
                ["1:0,30:30:99", "0:30,0:30:99", "./.:.:.:."],
            ),
            // the father has no call either: nothing tells
            (
                "chrY",
                15_000_100,
                "C",
                "G",
                ["1:0,30:30:99", "./.:.:.:.", "./.:.:.:."],
            ),
            // the father's allele: inherited, no row
            (
                "chrY",
                15_000_200,
                "G",
                "T",
                ["1:0,30:30:99", "1:0,30:30:99", "./.:.:.:."],
            ),
        ],
    );
    let out = report(&d, &m, &vcf, None);
    assert_eq!(
        de_novo_at(&out, "chrY", 15_000_000),
        vec![("1".to_owned(), "pass".to_owned())],
        "a son's Y allele his father does not have, the mother with no call: one pass row"
    );
    assert!(
        de_novo_at(&out, "chrY", 15_000_100).is_empty(),
        "no father call on Y: no row"
    );
    assert!(
        de_novo_at(&out, "chrY", 15_000_200).is_empty(),
        "the father's allele on Y: no row"
    );
}

#[test]
fn a_sons_x_de_novo_needs_only_the_mother() {
    let d = scratch("son-x");
    let m = members(&d, "male");
    let vcf = joint_vcf(
        &d,
        &[
            // the father has no call on X
            (
                "chrX",
                20_000_000,
                "A",
                "T",
                ["1:0,30:30:99", "./.:.:.:.", GOOD_REF],
            ),
            // the father's call is weak: it does not count on a son's X
            (
                "chrX",
                20_000_100,
                "C",
                "G",
                ["1:0,30:30:99", "0:3,0:3:5", GOOD_REF],
            ),
            // the mother's call is weak: low
            (
                "chrX",
                20_000_200,
                "G",
                "T",
                ["1:0,30:30:99", "0:30,0:30:99", "0/0:5,0:5:10"],
            ),
            // the father carries the allele, the mother does not: still new on a son's X
            (
                "chrX",
                20_000_300,
                "T",
                "C",
                ["1:0,30:30:99", "1:0,30:30:99", GOOD_REF],
            ),
        ],
    );
    let out = report(&d, &m, &vcf, None);
    assert_eq!(
        de_novo_at(&out, "chrX", 20_000_000),
        vec![("1".to_owned(), "pass".to_owned())],
        "a son's X allele his mother does not have, the father with no call: one pass row"
    );
    assert_eq!(
        de_novo_at(&out, "chrX", 20_000_100),
        vec![("1".to_owned(), "pass".to_owned())],
        "the father's depth and quality do not count on a son's X"
    );
    assert_eq!(
        de_novo_at(&out, "chrX", 20_000_200),
        vec![("1".to_owned(), "low".to_owned())],
        "the mother's depth and quality still count"
    );
    assert_eq!(
        de_novo_at(&out, "chrX", 20_000_300),
        vec![("1".to_owned(), "pass".to_owned())],
        "a father's X allele is not a son's"
    );
}

// ---- 3. parental reads on the new allele -----------------------------------

#[test]
fn a_parents_reads_on_the_new_allele_keep_it_out_of_pass() {
    let d = scratch("parent-reads");
    let m = members(&d, "female");
    let kid = "0/1:15,15:30:99";
    let vcf = joint_vcf(
        &d,
        &[
            // a mosaic father: 3 of 30 reads on the new allele
            (
                "chr2",
                5_000_000,
                "A",
                "G",
                [kid, "0/0:27,3:30:40", GOOD_REF],
            ),
            // the same trio, no reads on it
            (
                "chr2",
                5_000_100,
                "A",
                "G",
                [kid, "0/0:30,0:30:40", GOOD_REF],
            ),
            // one read of 15 (7%) in the mother
            (
                "chr2",
                5_000_200,
                "C",
                "T",
                [kid, GOOD_REF, "0/0:14,1:15:40"],
            ),
            // one read of 30 (3%) in the mother
            (
                "chr2",
                5_000_300,
                "C",
                "T",
                [kid, GOOD_REF, "0/0:29,1:30:40"],
            ),
            // two reads on a new allele of a 1/2 child, in the mother
            (
                "chr2",
                5_000_400,
                "C",
                "A,T",
                ["1/2:0,15,15:30:99", "0/1:15,15,0:30:99", "0/0:58,0,2:60:40"],
            ),
        ],
    );
    let out = report(&d, &m, &vcf, None);
    let at = |p| de_novo_at(&out, "chr2", p);
    assert_eq!(
        at(5_000_000),
        vec![("0/1".to_owned(), "low".to_owned())],
        "a father 0/0 with AD=27,3: not pass"
    );
    assert_eq!(
        at(5_000_100),
        vec![("0/1".to_owned(), "pass".to_owned())],
        "a father 0/0 with AD=30,0: pass"
    );
    assert_eq!(
        at(5_000_200),
        vec![("0/1".to_owned(), "low".to_owned())],
        "one read of 15 in the mother (5% or more): not pass"
    );
    assert_eq!(
        at(5_000_300),
        vec![("0/1".to_owned(), "pass".to_owned())],
        "one read of 30 in the mother (under 5%): pass"
    );
    assert_eq!(
        at(5_000_400),
        vec![("1/2".to_owned(), "low".to_owned())],
        "two reads on the new ALT of a 1/2 in the mother: not pass"
    );
}
