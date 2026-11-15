//! What a finding means for a child, beyond a recessive gene, and how sure
//! the call is: a gene that also causes a dominant condition, an X-linked
//! dominant gene, a male's call in a PAR when no PAR file is given, a
//! haploid call on an autosome, a call below step 06's floor, and the
//! condition ClinVar names for the allele. Each check reads the
//! repository's `assets/carrier_genes.tsv` at a gene's real place.

mod common;

use common::{Dir, HEAD, by_gene, cell, genes, og, read, repo};
use std::collections::{BTreeSet, HashMap};

const PED: &str = "sample\tsex\tfather\tmother\tpartner\nA\tfemale\t.\t.\tB\nB\tmale\t.\t.\tA\n";
const REVIEWED: &str = "criteria_provided,_multiple_submitters,_no_conflicts";
const FORMAT_HEAD: &str = "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n";
const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

/// A ClinVar pathogenic hit of step 06 in the middle of `gene`, with its
/// FORMAT and sample fields as given, and CLNDN when `clndn` is not empty.
fn hit(gene: &str, format: &str, sample: &str, clndn: &str) -> String {
    let g = &genes(&repo("assets/carrier_genes.tsv"))[gene];
    let dn = if clndn.is_empty() {
        String::new()
    } else {
        format!(";CLNDN={clndn}")
    };
    format!(
        "{}\t{}\t9001\tA\tG\t.\tPASS\tGENEINFO={gene}:1;CLNSIG=Pathogenic;CLNREVSTAT={REVIEWED}{dn}\t{format}\t{sample}\n",
        g.chrom,
        (g.start + g.end) / 2
    )
}

/// `og couple screen` on A's and B's hits (A female, B male); `par`: pass
/// the repository's PAR file. The couple table by gene, and the report.
fn screen(
    t: &Dir,
    name: &str,
    a: &str,
    b: &str,
    par: bool,
) -> (HashMap<String, HashMap<String, String>>, String) {
    let ped = t.put(&format!("{name}/pedigree.tsv"), PED);
    let va = t.put(&format!("{name}/A.vcf"), a);
    let vb = t.put(&format!("{name}/B.vcf"), b);
    let out = t.path(&format!("{name}/out"));
    let genes = repo("assets/carrier_genes.tsv");
    let par_bed = repo("assets/par_grch38.bed");
    let (ca, cb) = (format!("A={va}"), format!("B={vb}"));
    let mut args = vec![
        "--pedigree",
        &ped,
        "--genes",
        &genes,
        "--clinvar",
        &ca,
        "--clinvar",
        &cb,
        "--out-dir",
        &out,
    ];
    if par {
        args.extend(["--par-bed", &par_bed]);
    }
    let (rc, log) = og("screen", &args);
    assert_eq!(rc, 0, "{log}");
    (
        by_gene(&format!("{out}/couple_A_B.tsv")),
        read(&format!("{out}/couple_A_B.md")),
    )
}

fn plain(rows: &[String]) -> String {
    format!("{HEAD}{}", rows.concat())
}

#[test]
fn a_gene_with_a_dominant_condition_reads_that_a_child_can_be_affected() {
    let t = Dir::new("dominant");
    let a = plain(&[
        hit("MLH1", "GT", "0/1", "Lynch_syndrome_1|not_provided"),
        hit("CFTR", "GT", "0/1", ""),
    ]);
    let b = plain(&[hit(
        "BRCA1",
        "GT",
        "0/1",
        "Hereditary_breast_ovarian_cancer_syndrome",
    )]);
    let (rows, md) = screen(&t, "d", &a, &b, true);
    for (gene, who) in [("MLH1", "A"), ("BRCA1", "B")] {
        let why = cell(&rows, gene, "explanation");
        assert!(
            why.contains("a child who inherits it can be affected")
                && !why.contains("not affected"),
            "{gene}: {why}"
        );
        assert_eq!(cell(&rows, gene, "result"), "one_partner");
        assert_eq!(
            cell(&rows, gene, "offspring_risk"),
            "each child 1 in 2, dominant"
        );
        assert!(md.contains(why), "the report says it too");
        assert!(why.contains(&format!("Only {who} carries")), "{why}");
    }
    // ClinVar's condition for the allele is named, without its placeholder.
    let fa = cell(&rows, "MLH1", "findings_a");
    assert!(fa.ends_with(", Lynch syndrome 1"), "{fa}");
    assert!(!fa.contains("not provided"), "{fa}");
    assert!(
        cell(&rows, "BRCA1", "findings_b").ends_with(", Hereditary breast ovarian cancer syndrome")
    );
    // Control: a recessive gene with no dominant condition keeps its words.
    let cf = cell(&rows, "CFTR", "explanation");
    assert!(
        cf.contains("not affected") && !cf.contains("dominant"),
        "{cf}"
    );
    assert_eq!(cell(&rows, "CFTR", "offspring_risk"), "none shared");
    assert!(
        cell(&rows, "CFTR", "findings_a").ends_with(", ClinVar 9001"),
        "no CLNDN, no condition"
    );
}

#[test]
fn the_list_marks_the_dominant_genes_of_the_acmg_secondary_findings() {
    let list = genes(&repo("assets/carrier_genes.tsv"));
    let mut acmg: HashMap<String, String> = HashMap::new();
    for line in read(&repo("assets/acmg_sf.tsv")).lines() {
        let c: Vec<&str> = line.split('\t').collect();
        if !line.starts_with('#') && c[0] != "gene" {
            acmg.insert(c[0].to_owned(), c[2].to_owned());
        }
    }
    let marked: BTreeSet<&str> = list
        .iter()
        .filter(|(_, g)| matches!(g.dominant.as_str(), "AD" | "SD"))
        .map(|(n, _)| n.as_str())
        .collect();
    let want: BTreeSet<&str> = list
        .iter()
        .filter(|(n, g)| {
            g.inheritance == "AR" && acmg.get(*n).is_some_and(|m| m == "AD" || m == "SD")
        })
        .map(|(n, _)| n.as_str())
        .collect();
    assert_eq!(
        marked, want,
        "AD or SD exactly on the ACMG SF dominant genes"
    );
    assert!(want.len() > 40 && want.contains("MLH1") && want.contains("BRCA1"));
    for n in &marked {
        assert_eq!(
            list[*n].dominant, acmg[*n],
            "{n}: the value acmg_sf.tsv gives"
        );
    }
    for (n, g) in &list {
        assert!(
            g.dominant != "XLD" || g.inheritance == "XL",
            "{n}: an XLD mark is on an XL gene only"
        );
    }
}

#[test]
fn x_linked_dominant_genes_read_that_a_daughter_can_be_affected() {
    let list = genes(&repo("assets/carrier_genes.tsv"));
    for g in ["MECP2", "CDKL5", "IKBKG", "OFD1", "PCDH19"] {
        assert_eq!(
            list[g].inheritance, "XLD",
            "{g}: its main disorder is dominant"
        );
    }
    for g in ["GJB1", "GLA", "PLP1", "FMR1"] {
        assert_eq!(list[g].inheritance, "XL", "{g} stays XL");
        assert_eq!(list[g].dominant, "XLD", "{g}: females can show signs");
    }
    assert_eq!(list["OTC"].dominant, ".");
    let t = Dir::new("xld");
    let a = plain(&[
        hit("MECP2", "GT", "0/1", "Rett_syndrome"),
        hit("GLA", "GT", "0/1", ""),
    ]);
    let b = plain(&[hit("PCDH19", "GT", "1", ""), hit("OTC", "GT", "1", "")]);
    let (rows, _) = screen(&t, "x", &a, &b, true);
    let mecp2 = cell(&rows, "MECP2", "explanation");
    assert!(
        mecp2.contains("a daughter or a son who inherits it can be affected")
            && !mecp2.contains("carrier"),
        "{mecp2}"
    );
    assert_eq!(cell(&rows, "MECP2", "inheritance"), "XLD");
    assert_eq!(
        cell(&rows, "MECP2", "offspring_risk"),
        "A: each child 1 in 2, can be affected"
    );
    let pcdh19 = cell(&rows, "PCDH19", "explanation");
    assert!(
        pcdh19.contains("Every daughter inherits it and can be affected")
            && !pcdh19.contains("carrier"),
        "{pcdh19}"
    );
    let gla = cell(&rows, "GLA", "explanation");
    assert!(
        gla.contains("females can show signs") && gla.contains("carrier"),
        "{gla}"
    );
    assert_eq!(cell(&rows, "GLA", "inheritance"), "XL");
    // Control: a recessive X-linked gene with no dominant disorder.
    let otc = cell(&rows, "OTC", "explanation");
    assert!(
        otc.contains("Every daughter inherits it and is a carrier") && !otc.contains("show signs"),
        "{otc}"
    );
    assert_eq!(cell(&rows, "OTC", "result"), "x_linked");
}

#[test]
fn without_a_par_file_a_male_call_in_par1_is_one_copy_and_outside_hemizygous() {
    let list = genes(&repo("assets/carrier_genes.tsv"));
    let par = og_pedigree::Par::read(Some(&repo("assets/par_grch38.bed"))).unwrap();
    let in_par = |g: &str| par.contains(&list[g].chrom, (list[g].start + list[g].end) / 2);
    assert!(in_par("CSF2RA") && !in_par("DMD"));
    let t = Dir::new("par");
    let a = plain(&[hit("CSF2RA", "GT", "0/1", "")]);
    let b = plain(&[hit("CSF2RA", "GT", "0/1", ""), hit("DMD", "GT", "0/1", "")]);
    for with_file in [false, true] {
        let (rows, _) = screen(
            &t,
            if with_file { "file" } else { "none" },
            &a,
            &b,
            with_file,
        );
        assert_eq!(
            cell(&rows, "CSF2RA", "offspring_risk"),
            "1 in 4 per pregnancy",
            "PAR file given: {with_file}"
        );
        assert!(cell(&rows, "CSF2RA", "findings_b").contains(" het, "));
        assert!(cell(&rows, "DMD", "findings_b").contains(" hemizygous, "));
        assert_eq!(
            cell(&rows, "DMD", "offspring_risk"),
            "B: daughters carriers"
        );
    }
}

#[test]
fn a_haploid_call_on_an_autosome_is_one_copy() {
    let t = Dir::new("haploid");
    let a = plain(&[hit("CFTR", "GT", "0/1", "")]);
    let (rows, _) = screen(&t, "one", &a, &plain(&[hit("CFTR", "GT", "1", "")]), true);
    assert_eq!(
        cell(&rows, "CFTR", "offspring_risk"),
        "1 in 4 per pregnancy"
    );
    let fb = cell(&rows, "CFTR", "findings_b");
    assert!(
        fb.contains("one copy, a haploid call (1)") && !fb.contains(" hom"),
        "{fb}"
    );
    // Control: two copies are still two.
    let (rows, _) = screen(&t, "two", &a, &plain(&[hit("CFTR", "GT", "1/1", "")]), true);
    assert_eq!(
        cell(&rows, "CFTR", "offspring_risk"),
        "1 in 2 per pregnancy"
    );
}

#[test]
fn a_call_below_the_floor_is_marked_low_confidence_and_still_counts() {
    let t = Dir::new("low");
    let weak = format!(
        "{FORMAT_HEAD}{COLUMNS}{}",
        hit("CFTR", "GT:GQ:DP", "0/1:1:1", "")
    );
    let strong = format!(
        "{FORMAT_HEAD}{COLUMNS}{}",
        hit("CFTR", "GT:GQ:DP", "0/1:99:40", "")
    );
    let (rows, md) = screen(&t, "low", &weak, &strong, true);
    assert_eq!(cell(&rows, "CFTR", "result"), "shared", "it still counts");
    assert!(
        cell(&rows, "CFTR", "findings_a")
            .ends_with(", low-confidence call (GQ 1, DP 1): check the reads")
    );
    assert!(!cell(&rows, "CFTR", "findings_b").contains("low-confidence"));
    let why = cell(&rows, "CFTR", "explanation");
    assert!(
        why.contains("The finding of A rests on a low-confidence call"),
        "{why}"
    );
    assert!(md.contains("| A | 1 / 1 (1 low-confidence) |"), "{md}");
    // Control: the run's own floor (the file's ##PGP_CALL_FLOOR line) of 0 marks nothing.
    let floor0 = format!(
        "{FORMAT_HEAD}##PGP_CALL_FLOOR=GQ=0,DP=0\n{COLUMNS}{}",
        hit("CFTR", "GT:GQ:DP", "0/1:1:1", "")
    );
    let (rows, _) = screen(&t, "floor0", &floor0, &strong, true);
    assert!(!cell(&rows, "CFTR", "findings_a").contains("low-confidence"));
    assert!(!cell(&rows, "CFTR", "explanation").contains("low-confidence"));
}

#[test]
fn genes_keeps_x_linked_dominant_only_where_the_main_disorder_is() {
    let t = Dir::new("genes-xld");
    let par = repo("assets/par_grch38.bed");
    let gtf = t.put(
        "g.gtf",
        &[
            ("DOMX", "chrX", 5000000, 5000100),
            ("MIXX", "chrX", 6000000, 6000100),
            ("MECP2", "chrX", 7000000, 7000100),
            ("RECX", "chrX", 8000000, 8000100),
            ("AUT1", "chr3", 10, 20),
        ]
        .iter()
        .map(|(g, c, s, e)| {
            format!("{c}\tHAVANA\tgene\t{s}\t{e}\t.\t+\t.\tgene_id \"X\"; gene_name \"{g}\";\n")
        })
        .collect::<String>(),
    );
    let modes = [
        ("1", "X-linked dominant"),
        ("2", "X-linked recessive"),
        ("3", "Autosomal recessive"),
    ];
    let p9 = t.put(
        "p9.xml",
        &format!(
            "<JDBOR><DisorderList>{}</DisorderList></JDBOR>",
            modes
                .iter()
                .map(|(c, m)| format!(
                    "<Disorder><OrphaCode>{c}</OrphaCode><TypeOfInheritanceList><TypeOfInheritance><Name>{m}</Name>\
                     </TypeOfInheritance></TypeOfInheritanceList></Disorder>"
                ))
                .collect::<String>()
        ),
    );
    let dis = |code: &str, name: &str, sym: &str| {
        format!(
            "<Disorder><OrphaCode>{code}</OrphaCode><Name>{name}</Name><DisorderGeneAssociationList>\
             <DisorderGeneAssociation><Gene><Symbol>{sym}</Symbol></Gene><DisorderGeneAssociationType><Name>\
             Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType>\
             <DisorderGeneAssociationStatus><Name>Assessed</Name></DisorderGeneAssociationStatus>\
             </DisorderGeneAssociation></DisorderGeneAssociationList></Disorder>"
        )
    };
    let p6 = t.put(
        "p6.xml",
        &format!(
            "<JDBOR date=\"2026-01-02 03:04:05\"><DisorderList>{}</DisorderList></JDBOR>",
            [
                dis("1", "Dominant only", "DOMX"),
                dis("1", "Z dominant", "MIXX"),
                dis("2", "A recessive", "MIXX"),
                dis("1", "Z rett", "MECP2"),
                dis("2", "A recessive mecp2", "MECP2"),
                dis("2", "Recessive only", "RECX"),
                dis("3", "Autosomal one", "AUT1"),
            ]
            .concat()
        ),
    );
    let seed = t.put(
        "seed.tsv",
        "gene\tchrom\tstart\tend\tinheritance\tlists\tcondition\tseverity\tdominant\n\
         AUT1\t.\t0\t0\tAR\tOrphanet\tx\t.\tAD\n",
    );
    let new = t.path("new.tsv");
    let (rc, log) = og(
        "genes",
        &[
            "--table",
            &seed,
            "--orphanet-genes",
            &p6,
            "--orphanet-inheritance",
            &p9,
            "--gtf",
            &gtf,
            "--par-bed",
            &par,
            "--out",
            &new,
        ],
    );
    assert_eq!(rc, 0, "{log}");
    let g = genes(&new);
    let row = |n: &str| {
        (
            g[n].inheritance.as_str(),
            g[n].dominant.as_str(),
            g[n].condition.as_str(),
        )
    };
    assert_eq!(row("DOMX"), ("XLD", ".", "Dominant only"));
    assert_eq!(
        row("MIXX"),
        ("XL", "XLD", "A recessive; and 1 other"),
        "mixed: XL, females can show signs"
    );
    assert_eq!(
        row("MECP2"),
        ("XLD", ".", "Z rett; and 1 other"),
        "a named mixed gene is XLD, named by its dominant disorder"
    );
    assert_eq!(row("RECX"), ("XL", ".", "Recessive only"));
    assert_eq!(
        row("AUT1"),
        ("AR", "AD", "Autosomal one"),
        "the AD mark is kept"
    );
    assert!(read(&new).contains("\tseverity\tdominant\n"));
}
