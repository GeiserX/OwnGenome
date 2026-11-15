//! `og couple screen` and `og couple genes` (step 39 and the Nextflow
//! COUPLE_CARRIER process) on synthetic step outputs at the real coordinates
//! of assets/carrier_genes.tsv. These are the checks of the Python test the
//! port replaced (tests/test_couple_report.py), section by section:
//!
//!  1. the carrier list: every ACMG tier 3 gene (97 autosomal recessive, 16
//!     X-linked) is in it, the ACOG genes are marked, every row has a GRCh38
//!     span on chr1 to chrX, and a chrX gene outside the PARs is X-linked;
//!  2. one variant in the same recessive gene in each partner is one shared
//!     gene, a variant in another gene is one partner only, a gene outside
//!     the list is left out;
//!  3. Delly deletions: over an exon, not inside an intron, not filtered,
//!     not without the exon list; the read-depth rule for long ones, the
//!     length rule and the common-deletion rule;
//!  4. a homozygous hit gives 1 in 2, two het hits in one partner say the
//!     chance can be higher;
//!  5. X-linked genes and repeats;
//!  6. SMN1 (Parascopy) and CYP21A2 not assessed;
//!  7. the residual-risk lines name what was not assessed; the partner
//!     tables are the couple table;
//!  8. red controls;
//!  9. the same inputs give the same bytes;
//! 10. `genes`: the list rebuilt from a small Orphanet and GENCODE extract;
//! 11. the carrier deletion alleles and AZF regions of step 35's interval
//!     table, and the shared risk they add.
//!
//! Section 12 (a split half call is one copy) is the unit test
//! `zygosity_counts_a_missing_allele_in_the_ploidy` beside the rule, and the
//! half-call wording is checked here.

mod common;

use common::{Dir, HEAD, by_gene, cell, genes, og, read, repo, rows};
use std::collections::{BTreeSet, HashMap};

const ACMG: &str = "HBB XPC TYR CYP21A2 PAH CFTR TNXB HEXA GJB2 DHCR7 ATP7B ASPA ACADM NPHS1 PMM2 FKTN SLC26A4 ERCC2
DYNC2H1 CEP290 GBE1 GAA CHRNE G6PC1 OCA2 COL7A1 ABCC8 ALDOB FANCC GRIP1 BCKDHB ANO10 NAGA SMPD1 USH2A MMUT CPT2
AHI1 DHDDS SLC19A3 GALT CYP11A1 TF MMACHC ABCA3 GBA1 MCOLN1 GNPTAB AGA PCDH15 FAH AIRE BBS2 CYP27A1 CCDC88C FMO3
TMEM216 CNGB3 MCPH1 SLC37A4 PRF1 SCO2 AGXT ACADVL ASL EVC2 ARSA MVK PKHD1 BTD ALPL BBS1 CLCN1 CYP27B1 POLG MCCC2
MLC1 ACAT1 CC2D2A SLC26A2 CBS LRP2 IDUA FKRP RNASEH2B RARS2 HBA1 HBA2 SMN1 HPS1 HPS3 ELP1 FXN DLD NEB CLRN1
BLM";
const ACMG_XL: &str = "ABCD1 AFF2 ARX DMD F8 F9 FMR1 GLA L1CAM MID1 NR0B1 OTC PLP1 RPGR RS1 SLC6A8";
const ACOG: &str = "CFTR SMN1 HBB HBA1 HBA2 FMR1 HEXA ASPA ELP1 BLM ABCC8 FANCC GBA1 G6PC1 TMEM216 BCKDHB MCOLN1 SMPD1 PCDH15 CLRN1";
const REVIEWED: &str = "criteria_provided,_multiple_submitters,_no_conflicts";

fn words(s: &str) -> BTreeSet<String> {
    s.split_whitespace().map(str::to_owned).collect()
}

/// The carrier list and the PAR file of the repository.
struct Fixture {
    list: HashMap<String, common::Gene>,
    genes: String,
    par: String,
}

impl Fixture {
    fn new() -> Fixture {
        let genes = repo("assets/carrier_genes.tsv");
        Fixture {
            list: common::genes(&genes),
            genes,
            par: repo("assets/par_grch38.bed"),
        }
    }

    /// The middle of a gene of the list, moved by `off`.
    fn mid(&self, gene: &str, off: i64) -> (String, i64) {
        let g = &self.list[gene];
        (g.chrom.clone(), (g.start + g.end).div_euclid(2) + off)
    }

    /// A ClinVar pathogenic hit of step 06 in `gene`.
    fn hit(&self, gene: &str, gt: &str, off: i64) -> String {
        let (c, p) = self.mid(gene, off);
        hit_at(gene, gt, off, &c, p)
    }
}

fn hit_at(gene: &str, gt: &str, off: i64, c: &str, p: i64) -> String {
    format!(
        "{c}\t{p}\t9000{}\tA\tG\t.\tPASS\tGENEINFO={gene}:1;CLNSIG=Pathogenic;CLNREVSTAT={REVIEWED}\tGT\t{gt}\n",
        off % 97
    )
}

fn deletion(c: &str, s: i64, e: i64, filt: &str) -> String {
    format!("{c}\t{s}\tDEL0\tN\t<DEL>\t.\t{filt}\tSVTYPE=DEL;END={e}\tGT:FT\t0/1:PASS\n")
}

/// A Delly 2.x deletion with its read-depth copy number (FORMAT RDCN).
fn deldepth(c: &str, s: i64, e: i64, rdcn: &str, gt: &str, precise: &str) -> String {
    format!(
        "{c}\t{s}\tDEL1\tN\t<DEL>\t.\tPASS\t{precise};SVTYPE=DEL;END={e}\tGT:FT:RDCN\t{gt}:PASS:{rdcn}\n"
    )
}

/// An ExpansionHunter record; alleles: the repeat count of each allele.
fn eh(locus: &str, alleles: &[u32]) -> String {
    let gt = vec!["1"; alleles.len()].join("/");
    let sizes: Vec<String> = alleles.iter().map(u32::to_string).collect();
    format!(
        "chrX\t1\t.\tC\t<STR>\t.\tPASS\tVARID={locus};REPID={locus}\tGT:REPCN\t{gt}:{}\n",
        sizes.join("/")
    )
}

fn smn(cn: u32, qual: &str) -> String {
    format!(
        "{}chr5\t70951000\t70953000\tSMN1\tPASS\t{}\t40\tPASS\t{cn},2\t{qual},40\tchr5:70076000-70078000\n",
        common::PARASCOPY_HEAD,
        cn + 2
    )
}

const PED: &str = "sample\tsex\tfather\tmother\tpartner\nA\tfemale\t.\t.\tB\nB\tmale\t.\t.\tA\n";

#[test]
fn the_carrier_list_holds_the_panels() {
    let f = Fixture::new();
    let acmg: BTreeSet<String> = f
        .list
        .iter()
        .filter(|(_, g)| g.lists.contains("ACMG"))
        .map(|(n, _)| n.clone())
        .collect();
    let want: BTreeSet<String> = words(ACMG).union(&words(ACMG_XL)).cloned().collect();
    assert_eq!(acmg, want, "all 113 ACMG tier 3 genes are in the list");
    assert_eq!(
        words(ACMG)
            .iter()
            .filter(|g| f.list[*g].inheritance == "AR")
            .count(),
        97
    );
    assert!(words(ACMG_XL).iter().all(|g| f.list[g].inheritance == "XL"));
    let acog: BTreeSet<String> = f
        .list
        .iter()
        .filter(|(_, g)| g.lists.contains("ACOG"))
        .map(|(n, _)| n.clone())
        .collect();
    assert_eq!(acog, words(ACOG), "the ACOG genes are marked");
    let primary: Vec<String> = (1..=22)
        .map(|i| format!("chr{i}"))
        .chain(["chrX".to_owned()])
        .collect();
    assert!(
        f.list
            .values()
            .all(|g| primary.contains(&g.chrom) && 0 < g.start && g.start <= g.end),
        "every row has a span on chr1 to chrX"
    );
    let par = og_pedigree::Par::read(Some(&f.par)).unwrap();
    for (n, g) in &f.list {
        let outside = g.chrom == "chrX" && !par.contains("chrX", g.start);
        assert_eq!(
            g.inheritance == "XL",
            outside,
            "{n}: a chrX gene outside the PARs is X-linked, inside a PAR autosomal"
        );
    }
    assert!(
        f.list.len() > 3000,
        "more than the panel: the Orphanet genes"
    );
}

#[test]
fn a_shared_gene_deletions_over_exons_and_the_report() {
    let f = Fixture::new();
    let t = Dir::new("shared");
    let ped = t.put("pedigree.tsv", PED);
    let (cx, cy) = f.mid("CYP21A2", 0);
    let exons = t.put(
        "exons.bed",
        &format!(
            "{cx}\t{}\t{}\n{cx}\t{}\t{}\n",
            cy - 200,
            cy + 200,
            cy + 5000,
            cy + 5100
        ),
    );
    let (tx, ty) = f.mid("TNXB", 20000);
    let a_hits = t.put(
        "A.vcf",
        &format!(
            "{HEAD}{}{}",
            f.hit("NAGA", "0/1", 0),
            hit_at("APOE", "0/1", 0, "chr19", 44908684)
        ),
    );
    let b_hits = t.put(
        "B.vcf",
        &format!(
            "{HEAD}{}{}{}{}",
            f.hit("NAGA", "0/1", 7),
            f.hit("UGT1A1", "0/1", 0),
            f.hit("CYP21A2", "0/1", 0),
            f.hit("TNXB", "0/1", 20000)
        ),
    );
    let a_sv = t.put(
        "A_sv.vcf",
        &format!(
            "{HEAD}{}{}{}",
            deletion(&cx, cy - 1000, cy + 1000, "PASS"),
            // inside a TNXB intron: no exon in the list
            deletion(&tx, ty - 300, ty + 300, "PASS"),
            deletion(&cx, cy - 1000, cy + 1000, "LowQual")
        ),
    );
    let ca = format!("A={a_hits}");
    let cb = format!("B={b_hits}");
    let dl = format!("A={a_sv}");
    let base = [
        "--pedigree",
        &ped,
        "--genes",
        &f.genes,
        "--par-bed",
        &f.par,
        "--clinvar",
        &ca,
        "--clinvar",
        &cb,
        "--delly",
        &dl,
    ];
    let out = t.path("out");
    let mut args = base.to_vec();
    args.extend(["--exons", &exons, "--out-dir", &out]);
    let (rc, log) = og("screen", &args);
    assert_eq!(rc, 0, "{log}");
    let tab = rows(&format!("{out}/couple_A_B.tsv"));
    let res = by_gene(&format!("{out}/couple_A_B.tsv"));
    // 2.
    assert_eq!(cell(&res, "NAGA", "result"), "shared");
    assert_eq!(cell(&res, "NAGA", "offspring_risk"), "1 in 4 per pregnancy");
    assert_eq!(cell(&res, "UGT1A1", "result"), "one_partner");
    assert!(!res.contains_key("APOE") && !f.list.contains_key("APOE"));
    // 3.
    let fa = cell(&res, "CYP21A2", "findings_a");
    assert_eq!(cell(&res, "CYP21A2", "result"), "shared", "{fa}");
    assert!(
        fa.contains("deletion") && fa.contains("1 exon of the gene (Delly)"),
        "{fa}"
    );
    assert_eq!(
        fa.matches("deletion").count(),
        1,
        "the filtered copy is not counted: {fa}"
    );
    assert_eq!(cell(&res, "TNXB", "result"), "one_partner");
    assert_eq!(cell(&res, "TNXB", "findings_a"), ".");
    let results: Vec<&str> = tab.iter().map(|r| r["result"].as_str()).collect();
    assert_eq!(
        &results[..2],
        ["shared", "shared"],
        "shared genes come first"
    );
    let table = read(&format!("{out}/couple_A_B.tsv"));
    assert_eq!(read(&format!("{out}/A_couple_carrier.tsv")), table);
    assert_eq!(read(&format!("{out}/B_couple_carrier.tsv")), table);
    let md = read(&format!("{out}/couple_A_B.md"));
    assert!(md.split("## X-linked").next().unwrap().contains("**NAGA**"));
    assert!(md.contains("genetic counsellor"));
    // 7.
    assert!(md.contains("Deletions of B: no Delly calls"), "{md}");
    assert!(md.contains("SMN1 copy number of A: no Parascopy"), "{md}");
    assert!(
        md.contains("FMR1 and FXN repeats of B: no ExpansionHunter"),
        "{md}"
    );
    assert!(log.contains("A and B: 2 shared genes"), "{log}");
    // 9.
    let out2 = t.path("out2");
    let mut args = base.to_vec();
    args.extend(["--exons", &exons, "--out-dir", &out2]);
    og("screen", &args);
    for n in ["couple_A_B.tsv", "couple_A_B.md"] {
        assert_eq!(
            read(&format!("{out}/{n}")),
            read(&format!("{out2}/{n}")),
            "{n}"
        );
    }
    // 3. without the exon list
    let noexons = t.path("noexons");
    let mut args = base.to_vec();
    args.extend(["--out-dir", &noexons]);
    og("screen", &args);
    let r3 = by_gene(&format!("{noexons}/couple_A_B.tsv"));
    assert_eq!(cell(&r3, "CYP21A2", "result"), "one_partner");
    assert!(read(&format!("{noexons}/couple_A_B.md")).contains("without the GENCODE exon list"));
}

#[test]
fn homozygous_and_two_hets() {
    let f = Fixture::new();
    let t = Dir::new("hom");
    let ped = t.put("pedigree.tsv", PED);
    let a = t.put("A_hom.vcf", &format!("{HEAD}{}", f.hit("NAGA", "1/1", 0)));
    let b = t.put(
        "B_two.vcf",
        &format!(
            "{HEAD}{}{}",
            f.hit("NAGA", "0/1", 7),
            f.hit("NAGA", "0/1", 9)
        ),
    );
    let out = t.path("hom");
    let (ca, cb) = (format!("A={a}"), format!("B={b}"));
    og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &f.genes,
            "--clinvar",
            &ca,
            "--clinvar",
            &cb,
            "--out-dir",
            &out,
        ],
    );
    let r4 = &rows(&format!("{out}/couple_A_B.tsv"))[0];
    assert_eq!(r4["offspring_risk"], "1 in 2 per pregnancy");
    assert!(
        r4["explanation"].contains("B has two changes in the gene"),
        "{r4:?}"
    );
}

#[test]
fn x_linked_genes_and_repeats() {
    let f = Fixture::new();
    let t = Dir::new("xl");
    let ped = t.put("pedigree.tsv", PED);
    let a_eh = t.put(
        "A_eh.vcf",
        &format!(
            "{}{}{}",
            HEAD.replace("\tS\n", "\tA\n"),
            eh("FMR1", &[30, 85]),
            eh("FXN", &[9, 120])
        ),
    );
    let b_eh = t.put(
        "B_eh.vcf",
        &format!("{HEAD}{}{}", eh("FMR1", &[29]), eh("FXN", &[12, 700])),
    );
    let b_dmd = t.put("B_dmd.vcf", &format!("{HEAD}{}", f.hit("DMD", "1", 0)));
    let a_none = t.put("A_none.vcf", HEAD);
    let out = t.path("xl");
    let args: Vec<String> = [
        "--pedigree",
        &ped,
        "--genes",
        &f.genes,
        "--par-bed",
        &f.par,
        "--clinvar",
        &format!("A={a_none}"),
        "--clinvar",
        &format!("B={b_dmd}"),
        "--eh",
        &format!("A={a_eh}"),
        "--eh",
        &format!("B={b_eh}"),
        "--out-dir",
        &out,
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    og("screen", &args);
    let r5 = by_gene(&format!("{out}/couple_A_B.tsv"));
    assert_eq!(cell(&r5, "FMR1", "result"), "x_linked");
    assert_eq!(cell(&r5, "FMR1", "offspring_risk"), "A: sons 1 in 2");
    assert!(cell(&r5, "FMR1", "findings_a").contains("premutation (55 to 200)"));
    assert!(cell(&r5, "FMR1", "explanation").contains("grow into a full mutation"));
    assert_eq!(
        cell(&r5, "FMR1", "findings_b"),
        ".",
        "the male's FMR1 allele of 29 is no finding"
    );
    assert_eq!(cell(&r5, "DMD", "offspring_risk"), "B: daughters carriers");
    assert!(cell(&r5, "DMD", "findings_b").contains("hemizygous"));
    assert_eq!(cell(&r5, "FXN", "result"), "shared", "FXN expanded in both");
}

/// Runs a screen of A and B with no ClinVar hit plus `extra`: the table by
/// gene and the report.
fn empty_screen(
    f: &Fixture,
    t: &Dir,
    name: &str,
    extra: &[&str],
) -> (u8, String, HashMap<String, HashMap<String, String>>, String) {
    let ped = t.put("pedigree.tsv", PED);
    let none = t.put("none.vcf", HEAD);
    let (ca, cb) = (format!("A={none}"), format!("B={none}"));
    let out = t.path(name);
    let mut args = vec![
        "--pedigree",
        &ped,
        "--genes",
        &f.genes,
        "--par-bed",
        &f.par,
        "--clinvar",
        &ca,
        "--clinvar",
        &cb,
    ];
    args.extend_from_slice(extra);
    args.extend(["--out-dir", &out]);
    let (rc, log) = og("screen", &args);
    (
        rc,
        log,
        by_gene(&format!("{out}/couple_A_B.tsv")),
        read(&format!("{out}/couple_A_B.md")),
    )
}

#[test]
fn smn1_and_cyp21a2() {
    let f = Fixture::new();
    let t = Dir::new("smn");
    let (sa, sb, su) = (
        t.put("A_smn.tsv", &smn(1, "45")),
        t.put("B_smn.tsv", &smn(1, "45")),
        t.put("U_smn.tsv", &smn(1, "5")),
    );
    let (a, b) = (format!("A={sa}"), format!("B={sb}"));
    let (_, _, r6, _) = empty_screen(&f, &t, "smn", &["--smn", &a, "--smn", &b]);
    assert_eq!(cell(&r6, "SMN1", "result"), "shared", "{r6:?}");
    let u = format!("A={su}");
    let (_, _, r, md) = empty_screen(&f, &t, "smnu", &["--smn", &u, "--smn", &b]);
    assert_eq!(
        cell(&r, "SMN1", "findings_a"),
        ".",
        "a call of quality 5 is no finding: {r:?}"
    );
    assert!(
        md.contains(
            "SMN1 copy number of A: SMN1 copy number 1 (quality 5, filter PASS), not resolved"
        ),
        "{md}"
    );
    let first = common::PARASCOPY_HEAD;
    let sn = t.put(
        "N_smn.tsv",
        &format!(
            "{first}chr5\t70951945\t70951946\tSMN1\tnot assessed: no read over chr5:70951946-70951946\t.\t.\t\
             not assessed: no read over chr5:70951946-70951946\t.\t.\t.\n"
        ),
    );
    let n = format!("A={sn}");
    let (_, _, r, md) = empty_screen(&f, &t, "smnn", &["--smn", &n, "--smn", &b]);
    assert_eq!(cell(&r, "SMN1", "findings_a"), ".", "{r:?}");
    assert!(
        md.contains(
            "SMN1 copy number of A: SMN1 not assessed: no read over chr5:70951946-70951946; copy number not resolved."
        ),
        "{md}"
    );
    // 6b. CYP21A2 not assessed: the table's message, the reason code for an older table
    let why = "not assessed: the BAM has few or no multi-mapping reads here (reads of mapping quality 0 make up \
               0.0% of the depth over CYP21A2 and CYP21A1P, under 2%), the sign of an aligner that leaves out the \
               reads that fit both copies";
    let na = "metric\tvalue\nstatus\tnot assessed\ncyp21a2_copies\t.\nfinding\t.\nreason\tlow_mapq0_share\n";
    let ca = t.put("A_cyp.tsv", &format!("{na}message\t{why}\n"));
    let cb = t.put("B_cyp.tsv", na);
    let (a, b) = (format!("A={ca}"), format!("B={cb}"));
    let (_, _, _, md) = empty_screen(&f, &t, "cyp", &["--cyp21a2", &a, "--cyp21a2", &b]);
    assert!(
        md.contains(&format!(
            "CYP21A2 deletion and conversions of A: CYP21A2 {why}."
        )),
        "{md}"
    );
    assert!(
        md.contains(
            "CYP21A2 deletion and conversions of B: CYP21A2 not assessed (low_mapq0_share)."
        ),
        "{md}"
    );
}

#[test]
fn read_depth_decides_long_deletions() {
    let f = Fixture::new();
    let t = Dir::new("depth");
    let m: HashMap<&str, (String, i64)> = ["GBA1", "PAH", "CFTR", "HBB", "DMD", "OTC"]
        .iter()
        .map(|g| (*g, f.mid(g, 0)))
        .collect();
    let mut spots: Vec<&(String, i64)> = m.values().collect();
    spots.sort();
    let dexons = t.put(
        "depth_exons.bed",
        &spots
            .iter()
            .map(|(c, p)| format!("{c}\t{}\t{}\n", p - 50, p + 50))
            .collect::<String>(),
    );
    let mb = 25_000_000;
    let (gc, gp) = &m["GBA1"];
    let (pc, pp) = &m["PAH"];
    let (fc, fp) = &m["CFTR"];
    let (hc, hp) = &m["HBB"];
    let (dc, dp) = &m["DMD"];
    let (oc, op) = &m["OTC"];
    let a = t.put(
        "A_depth.vcf",
        &[
            HEAD.to_owned(),
            deldepth(gc, gp - mb, gp + mb, "2", "0/1", "IMPRECISE"), // 50 Mb, depth shows 2 copies
            deldepth(pc, pp - mb, pp + mb, "1", "0/1", "IMPRECISE"), // 50 Mb, depth shows 1 copy
            deldepth(fc, fp - 300, fp + 300, "2", "0/1", "PRECISE"), // 600 bp, split reads alone
            deldepth(hc, hp - 1000, hp + 1000, ".", "0/1", "IMPRECISE"), // 2 kb, no estimate
        ]
        .concat(),
    );
    let b = t.put(
        "B_depth.vcf",
        &[
            HEAD.to_owned(),
            deldepth(dc, dp - 1_000_000, dp + 1_000_000, "1", "1", "IMPRECISE"), // male chrX: 1 is normal
            deldepth(oc, op - 2500, op + 2500, "0", "1", "IMPRECISE"), // male chrX: 0, a loss
        ]
        .concat(),
    );
    let (da, db) = (format!("A={a}"), format!("B={b}"));
    let (rc, log, rd, md) = empty_screen(
        &f,
        &t,
        "depth",
        &["--delly", &da, "--delly", &db, "--exons", &dexons],
    );
    assert_eq!(rc, 0, "{log}");
    assert!(
        !rd.contains_key("GBA1"),
        "a 50 Mb deletion with RDCN 2 is set aside"
    );
    assert_eq!(cell(&rd, "PAH", "result"), "one_partner");
    assert!(
        cell(&rd, "PAH", "findings_a").contains("(50000.0 kb) het, 1 exon of the gene (Delly)")
    );
    assert_eq!(cell(&rd, "CFTR", "result"), "one_partner");
    assert!(cell(&rd, "CFTR", "findings_a").contains("(600 bp) het"));
    assert!(
        !rd.contains_key("HBB"),
        "a 2 kb deletion without a depth estimate is set aside"
    );
    assert!(
        !rd.contains_key("DMD"),
        "a male's chrX deletion with RDCN 1 is set aside"
    );
    assert_eq!(cell(&rd, "OTC", "result"), "x_linked");
    assert!(cell(&rd, "OTC", "findings_b").contains("hemizygous"));
    assert!(md.contains(
        "| A | 0 / 0 | 2 over exons of listed genes; 2 longer than 1000 bp set aside, read depth does not confirm the loss |"
    ));
    assert!(md.contains(
        "| B | 0 / 0 | 1 over exons of listed genes; 1 longer than 1000 bp set aside, read depth does not confirm the loss |"
    ));
}

#[test]
fn length_and_common_deletions() {
    let f = Fixture::new();
    let t = Dir::new("common");
    let names = [
        "ASPA", "HEXA", "SMPD1", "MCOLN1", "GJB2", "TPP1", "GALC", "NPC1", "CLN3",
    ];
    let m: HashMap<&str, (String, i64)> = names.iter().map(|g| (*g, f.mid(g, 0))).collect();
    let mut spots: Vec<&(String, i64)> = m.values().collect();
    spots.sort();
    let cexons = t.put(
        "common_exons.bed",
        &spots
            .iter()
            .map(|(c, p)| format!("{c}\t{}\t{}\n", p - 50, p + 50))
            .collect::<String>(),
    );
    let d = |g: &str, l: i64, r: i64, rdcn: &str, precise: &str| {
        let (c, p) = &m[g];
        deldepth(c, p - l, p + r, rdcn, "0/1", precise)
    };
    let a = t.put(
        "A_len.vcf",
        &[
            HEAD.to_owned(),
            d("ASPA", 15, 15, "2", "PRECISE"), // 30 bp over an ASPA exon
            d("HEXA", 2500, 2500, "1", "IMPRECISE"), // 5 kb, a common deletion
            d("SMPD1", 2500, 2500, "1", "IMPRECISE"), // 5 kb, no common match
            d("MCOLN1", 2500, 2500, "1", "IMPRECISE"), // 5 kb inside a 20 kb common one
            d("GJB2", 2500, 2500, "1", "IMPRECISE"), // 5 kb, a rare benign record
            d("TPP1", 25, 24, "2", "PRECISE"), // 49 bp, one short of an SV
            d("GALC", 25, 25, "2", "PRECISE"), // 50 bp, the smallest SV
            d("NPC1", 2500, 2500, "1", "IMPRECISE"), // 5 kb, a record of exactly 1%
            d("CLN3", 2500, 2500, "1", "IMPRECISE"), // 5 kb, half a 10 kb common one
        ]
        .concat(),
    );
    // The benign set as AnnotSV writes it: chromosome without "chr", 0-based
    // start, source, coordinates, frequency.
    let rec = |g: &str, l: i64, r: i64, af: &str| {
        let (c, p) = &m[g];
        (c[3..].to_owned(), p - l, p + r, af.to_owned())
    };
    let mut recs = [
        rec("HEXA", 2400, 2700, "0.3100"), // the HEXA deletion, breakpoints 100-200 bp off
        rec("MCOLN1", 10000, 10000, "0.2000"), // covers MCOLN1's deletion, 25% of itself
        rec("GJB2", 2500, 2500, "0.0040"), // GJB2's exactly, but below 1%
        rec("NPC1", 2500, 2500, "0.0100"), // NPC1's exactly, at 1%
        rec("CLN3", 2500, 7500, "0.2000"), // CLN3's 5 kb is 50% of this 10 kb record
    ];
    recs.sort_by(|x, y| (format!("chr{}", x.0), x.1, x.2).cmp(&(format!("chr{}", y.0), y.1, y.2)));
    let benign = t.put(
        "benign_Loss_SV_GRCh38.sorted.bed",
        &recs
            .iter()
            .map(|(c, s, e, af)| format!("{c}\t{s}\t{e}\tgnomAD-SV_fixture\t{c}:{s}-{e}\t{af}\n"))
            .collect::<String>(),
    );
    let da = format!("A={a}");
    let (rc, log, rc_, mdc) = empty_screen(
        &f,
        &t,
        "common",
        &["--delly", &da, "--exons", &cexons, "--benign-sv", &benign],
    );
    assert_eq!(rc, 0, "{log}");
    assert!(!rc_.contains_key("ASPA"), "30 bp is not reported");
    assert!(!rc_.contains_key("TPP1"), "49 bp is not reported");
    assert_eq!(cell(&rc_, "GALC", "result"), "one_partner");
    assert!(cell(&rc_, "GALC", "findings_a").contains("(50 bp) het"));
    assert!(
        !rc_.contains_key("HEXA"),
        "a common deletion is not reported"
    );
    assert!(
        !rc_.contains_key("NPC1"),
        "a record of exactly 1% sets it aside"
    );
    assert!(
        !rc_.contains_key("CLN3"),
        "exactly 50% of a common record sets it aside"
    );
    assert_eq!(cell(&rc_, "SMPD1", "result"), "one_partner");
    assert!(cell(&rc_, "SMPD1", "findings_a").contains("(5.0 kb) het"));
    assert_eq!(
        cell(&rc_, "MCOLN1", "result"),
        "one_partner",
        "25% of a common one still counts"
    );
    assert_eq!(
        cell(&rc_, "GJB2", "result"),
        "one_partner",
        "a record below 1% sets nothing aside"
    );
    assert!(mdc.contains(
        "| A | 0 / 0 | 4 over exons of listed genes; 2 shorter than 50 bp set aside; 3 set aside as common deletions |"
    ), "{mdc}");
    assert!(mdc.contains(
        "a deletion matching one of 1% or more in AnnotSV's benign set (4 records, 50% reciprocal overlap) is set aside"
    ), "{mdc}");
    let (rc, log, rn, mdn) =
        empty_screen(&f, &t, "nobenign", &["--delly", &da, "--exons", &cexons]);
    assert_eq!(rc, 0, "{log}");
    assert!(!rn.contains_key("ASPA"));
    assert_eq!(
        cell(&rn, "HEXA", "result"),
        "one_partner",
        "red: without a benign set HEXA is reported"
    );
    assert!(mdn.contains("Common deletions: not filtered, no benign SV set was given"));
    assert!(mdn.contains("Common deletions of A: not filtered"));
    let missing = t.path("no_such.bed");
    let (rc, log, _, _) = empty_screen(
        &f,
        &t,
        "nobed",
        &["--delly", &da, "--exons", &cexons, "--benign-sv", &missing],
    );
    assert!(rc == 2 && log.contains("no_such.bed"), "{rc} {log}");
}

#[test]
fn red_controls() {
    let f = Fixture::new();
    let t = Dir::new("red");
    let ped = t.put("pedigree.tsv", PED);
    let a_hits = t.put("A.vcf", &format!("{HEAD}{}", f.hit("NAGA", "0/1", 0)));
    let b_hits = t.put("B.vcf", &format!("{HEAD}{}", f.hit("NAGA", "0/1", 7)));
    let (ca, cb) = (format!("A={a_hits}"), format!("B={b_hits}"));
    let empty = t.put(
        "empty.tsv",
        "# emptied\ngene\tchrom\tstart\tend\tinheritance\tlists\tcondition\n",
    );
    let red = t.path("red");
    let (rc, log) = og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &empty,
            "--clinvar",
            &ca,
            "--clinvar",
            &cb,
            "--out-dir",
            &red,
        ],
    );
    assert!(rc == 2 && log.contains("holds no gene"), "{rc} {log}");
    assert!(!std::path::Path::new(&format!("{red}/couple_A_B.tsv")).exists());
    let lone = t.put(
        "lone.tsv",
        "sample\tsex\tfather\tmother\tpartner\nA\tfemale\t.\t.\t.\nB\tmale\t.\t.\t.\n",
    );
    let o = t.path("lone");
    let (rc, log) = og(
        "screen",
        &[
            "--pedigree",
            &lone,
            "--genes",
            &f.genes,
            "--clinvar",
            &ca,
            "--out-dir",
            &o,
        ],
    );
    assert!(
        rc == 1 && log.contains("no two people as partners"),
        "{rc} {log}"
    );
    let half = t.path("half");
    let (rc, log) = og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &f.genes,
            "--clinvar",
            &ca,
            "--out-dir",
            &half,
        ],
    );
    assert!(
        rc == 1 && log.contains("B has no ClinVar screen"),
        "{rc} {log}"
    );
    assert!(!std::path::Path::new(&format!("{half}/couple_A_B.tsv")).exists());
    let c = format!("C={a_hits}");
    let o = t.path("c");
    let (rc, log) = og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &f.genes,
            "--clinvar",
            &c,
            "--out-dir",
            &o,
        ],
    );
    assert!(
        rc == 2 && log.contains("'C' is not in the pedigree"),
        "{rc} {log}"
    );
}

const IVH: &str =
    "locus\tname\tchrom\tstart\tend\tstatus\tcopies\texpected\tchange\testimate\treason\n";

/// One row of step 35's interval table.
fn ivrow(locus: &str, chrom: &str, s: i64, e: i64, kv: &[(&str, &str)]) -> String {
    let get = |k: &str, d: &'static str| {
        kv.iter()
            .find(|(x, _)| *x == k)
            .map_or(d.to_owned(), |(_, v)| (*v).to_owned())
    };
    format!(
        "{locus}\t{locus} name\t{chrom}\t{s}\t{e}\t{}\t{}\t{}\t{}\t.\t{}\n",
        get("status", "called"),
        get("copies", "2"),
        get("expected", "2"),
        get("change", "no change"),
        get("reason", ".")
    )
}

const DEL1: &[(&str, &str)] = &[("copies", "1"), ("change", "deletion")];

#[test]
fn carrier_deletion_alleles_and_azf() {
    let f = Fixture::new();
    let t = Dir::new("cn");
    let ped = t.put("pedigree.tsv", PED);
    let a_hits = t.put("A.vcf", &format!("{HEAD}{}", f.hit("NAGA", "0/1", 0)));
    let b_hits = t.put("B.vcf", &format!("{HEAD}{}", f.hit("NAGA", "0/1", 7)));
    let a_iv = t.put(
        "A_cn_intervals.tsv",
        &[
            IVH.to_owned(),
            ivrow("GJB6_D13S1830", "chr13", 20223037, 20531806, DEL1),
            ivrow("CFTR_dele2_3", "chr7", 117498312, 117519392, &[]),
            ivrow(
                "HBB_exon3_del",
                "chr11",
                5225254,
                5225875,
                &[
                    ("status", "not assessed"),
                    ("copies", "."),
                    ("change", "."),
                    ("reason", "between_counts"),
                ],
            ),
            ivrow("HBB_619bp_del", "chr11", 5225388, 5226007, &[]),
            ivrow("HBD_HBB_Sicilian_del", "chr11", 5219877, 5233280, &[]),
            ivrow(
                "AZFa",
                "chrY",
                12537649,
                12920478,
                &[
                    ("status", "not applicable"),
                    ("copies", "."),
                    ("expected", "."),
                    ("change", "."),
                    ("reason", "no_y"),
                ],
            ),
        ]
        .concat(),
    );
    let b_iv = t.put(
        "B_cn_intervals.tsv",
        &[
            IVH.to_owned(),
            ivrow("GJB6_D13S1830", "chr13", 20223037, 20531806, &[]),
            ivrow("CFTR_dele2_3", "chr7", 117498312, 117519392, &[]),
            ivrow("HBB_exon3_del", "chr11", 5225254, 5225875, DEL1),
            ivrow("HBB_619bp_del", "chr11", 5225388, 5226007, DEL1),
            ivrow("HBD_HBB_Sicilian_del", "chr11", 5219877, 5233280, DEL1),
            ivrow(
                "AZFa",
                "chrY",
                12537649,
                12920478,
                &[("copies", "1"), ("expected", "1")],
            ),
            ivrow(
                "AZFb",
                "chrY",
                19702297,
                20593154,
                &[("copies", "1"), ("expected", "1")],
            ),
            ivrow(
                "AZFc_DAZ",
                "chrY",
                23129354,
                24907040,
                &[("copies", "2"), ("expected", "4"), ("change", "deletion")],
            ),
        ]
        .concat(),
    );
    // CNVpytor: no call over A's GJB6 deletion, a deletion over B's CFTRdele2,3 (absent by depth)
    let cnv_head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let b_cnv = t.put(
        "B_cnvs.vcf",
        &format!("{cnv_head}chr7\t117498001\tCNV1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=117520000\n"),
    );
    let a_cnv = t.put("A_cnvs.vcf", cnv_head);
    let out = t.path("cn");
    let (ca, cb) = (format!("A={a_hits}"), format!("B={b_hits}"));
    let (ia, ib) = (format!("A={a_iv}"), format!("B={b_iv}"));
    let (na, nb) = (format!("A={a_cnv}"), format!("B={b_cnv}"));
    let (rc, log) = og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &f.genes,
            "--clinvar",
            &ca,
            "--clinvar",
            &cb,
            "--cn-intervals",
            &ia,
            "--cn-intervals",
            &ib,
            "--cnvpytor",
            &na,
            "--cnvpytor",
            &nb,
            "--out-dir",
            &out,
        ],
    );
    assert_eq!(rc, 0, "{log}");
    let md = read(&format!("{out}/couple_A_B.md"));
    for p in ["A", "B"] {
        for n in [
            "del(GJB6-D13S1830)",
            "CFTRdele2,3",
            "HBB exon 3 deletion",
            "HBB 619 bp deletion",
            "Sicilian (delta-beta)0 deletion",
        ] {
            let states = ["present", "absent", "not assessable"]
                .iter()
                .filter(|st| md.contains(&format!("- {p}, {n}: **{st}**")))
                .count();
            assert_eq!(states, 1, "{p} {n}: one of the three states\n{md}");
        }
    }
    for want in [
        "- A, del(GJB6-D13S1830): **present** (read depth, step 35: 1 of 2 copies)",
        "- A, HBB exon 3 deletion: **not assessable** (read depth, step 35: not assessed (between_counts))",
        "- B, CFTRdele2,3: **absent**",
        "- B, Sicilian (delta-beta)0 deletion: **present** (read depth, step 35: 1 of 2 copies)",
        "- B, HBB 619 bp deletion: **present** (read depth, step 35: 1 of 2 copies)",
        "- A, Sicilian (delta-beta)0 deletion: **absent** (read depth, step 35: 2 of 2 copies)",
        // where CNVpytor and the depth call disagree, both are printed
        "- A, del(GJB6-D13S1830): **present** (read depth, step 35: 1 of 2 copies); CNVpytor (step 18) calls no deletion here",
        "- B, CFTRdele2,3: **absent** (read depth, step 35: 2 of 2 copies); CNVpytor (step 18) calls a deletion here",
        // no CNVpytor note under two of its 1 kb bins; the 13.4 kb Sicilian gets it
        "- B, HBB 619 bp deletion: **present** (read depth, step 35: 1 of 2 copies)\n",
        "- B, HBB exon 3 deletion: **present** (read depth, step 35: 1 of 2 copies)\n",
        "- B, Sicilian (delta-beta)0 deletion: **present** (read depth, step 35: 1 of 2 copies); CNVpytor (step 18) calls no deletion here",
        // the male partner's AZF regions; the female's not applicable
        "- B, AZFc deletion (DAZ copies; gr/gr and b2/b3 leave two of four): **present** (read depth, step 35: 2 of 4 copies)",
        "- B, AZFa deletion: **absent** (read depth, step 35: 1 of 1 copies)",
        "- A, AZF regions: not applicable (no Y chromosome)",
    ] {
        assert!(md.contains(want), "{want}\n{md}");
    }
    assert!(!md.contains("- B, GJB6_D13S1830"));
    assert!(
        format!("{md}\n")
            .contains("- B, del(GJB6-D13S1830): **absent** (read depth, step 35: 2 of 2 copies)\n")
    );
    assert!(!md.contains("- A, AZFa"));
    let nocn = t.path("nocn");
    let (rc, _) = og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &f.genes,
            "--clinvar",
            &ca,
            "--clinvar",
            &cb,
            "--out-dir",
            &nocn,
        ],
    );
    assert_eq!(rc, 0);
    assert!(!read(&format!("{nocn}/couple_A_B.md")).contains("Carrier deletion alleles"));
}

/// The shared-risk section of a couple report.
fn shared_head(md: &str) -> String {
    let after = md.split_once("## Shared risk").map_or(md, |(_, b)| b);
    after.split("## X-linked").next().unwrap().to_owned()
}

#[test]
fn read_depth_alleles_add_a_shared_risk() {
    let f = Fixture::new();
    let t = Dir::new("cn2");
    let ped2 = t.put(
        "ped2.tsv",
        "sample\tsex\tfather\tmother\tpartner\nC\tfemale\t.\t.\tD\nD\tmale\t.\t.\tC\n",
    );
    let c_hits = t.put("C.vcf", &format!("{HEAD}{}", f.hit("UGT1A1", "0/1", 0)));
    let d_hits = t.put("D.vcf", &format!("{HEAD}{}", f.hit("HBB", "0/1", 0)));
    let c_iv = t.put(
        "C_cn_intervals.tsv",
        &[
            IVH.to_owned(),
            ivrow("CFTR_dele2_3", "chr7", 117498312, 117519392, DEL1),
            ivrow("HBB_exon3_del", "chr11", 5225254, 5225875, DEL1),
        ]
        .concat(),
    );
    let d_iv = t.put(
        "D_cn_intervals.tsv",
        &[
            IVH.to_owned(),
            ivrow("CFTR_dele2_3", "chr7", 117498312, 117519392, DEL1),
        ]
        .concat(),
    );
    let screen2 = |genes: &str, d_vcf: &str, civ: &str, div: Option<&str>, name: &str| {
        let out = t.path(name);
        let (cc, cd, ic) = (
            format!("C={c_hits}"),
            format!("D={d_vcf}"),
            format!("C={civ}"),
        );
        let id = div.map(|d| format!("D={d}"));
        let mut args = vec![
            "--pedigree",
            &ped2,
            "--genes",
            genes,
            "--clinvar",
            &cc,
            "--clinvar",
            &cd,
            "--cn-intervals",
            &ic,
        ];
        if let Some(id) = &id {
            args.extend(["--cn-intervals", id]);
        }
        args.extend(["--out-dir", &out]);
        let (rc, log) = og("screen", &args);
        (
            rc,
            log,
            read(&format!("{out}/couple_C_D.md")),
            rows(&format!("{out}/couple_C_D.tsv")),
        )
    };
    let (rc, log, md, tab) = screen2(&f.genes, &d_hits, &c_iv, Some(&d_iv), "cn2");
    let head = shared_head(&md);
    assert_eq!(rc, 0, "{log}");
    assert!(!tab.iter().any(|r| r["result"] == "shared"));
    assert!(head.contains(
        "No gene of the list where both partners carry a change in the gene table; read depth (step 35) adds a possible \
         shared risk, below."
    ), "{head}");
    assert!(head.contains(
        "- **CFTR**: C: CFTRdele2,3 present (read depth, step 35: 1 of 2 copies); D: CFTRdele2,3 present (read depth, \
         step 35: 1 of 2 copies). Both partners may carry a change in CFTR"
    ), "{head}");
    assert!(head.contains(
        "- **HBB**: C: HBB exon 3 deletion present (read depth, step 35: 1 of 2 copies); D: "
    ));
    assert!(
        log.contains("possibly shared, from read depth: CFTR")
            && log.contains("possibly shared, from read depth: HBB")
    );
    let none_iv = t.put(
        "N_cn_intervals.tsv",
        &[
            IVH.to_owned(),
            ivrow("CFTR_dele2_3", "chr7", 117498312, 117519392, &[]),
        ]
        .concat(),
    );
    let (rc, log, md, _) = screen2(&f.genes, &d_hits, &none_iv, Some(&d_iv), "cn3");
    let head = shared_head(&md);
    assert_eq!(rc, 0);
    assert!(
        head.contains("No gene of the list where both partners carry a change.\n"),
        "red: {head}"
    );
    assert!(!head.contains("Possible shared risk") && !log.contains("possibly shared"));
    // DFNB1: one partner's del(GJB6-D13S1830) beside the other's GJB2 hit
    let g_hits = t.put("G.vcf", &format!("{HEAD}{}", f.hit("GJB2", "0/1", 0)));
    let g_iv = t.put(
        "G_cn_intervals.tsv",
        &[
            IVH.to_owned(),
            ivrow("GJB6_D13S1830", "chr13", 20223037, 20531806, DEL1),
        ]
        .concat(),
    );
    let (rc, log, md, _) = screen2(&f.genes, &g_hits, &g_iv, None, "cn4");
    let head = shared_head(&md);
    assert_eq!(rc, 0, "{log}");
    let dfnb1 = "- **GJB2/GJB6 (DFNB1)**: C: del(GJB6-D13S1830) present (read depth, step 35: 1 of 2 copies); D: ";
    assert!(head.contains(dfnb1), "{head}");
    assert!(head.split_once("; D: ").unwrap().1.contains("GJB2"));
    // A carrier list without CFTR: CFTRdele2,3 in both is no shared risk; HBB still is
    let no_cftr = t.put(
        "no_cftr.tsv",
        &read(&f.genes)
            .lines()
            .filter(|l| !l.starts_with("CFTR\t"))
            .map(|l| format!("{l}\n"))
            .collect::<String>(),
    );
    let (rc, log, md, _) = screen2(&no_cftr, &d_hits, &c_iv, Some(&d_iv), "cn5");
    let head = shared_head(&md);
    assert_eq!(rc, 0, "{log}");
    assert!(!head.contains("**CFTR**"));
    assert!(head.contains("- **HBB**: C: HBB exon 3 deletion present"));
    assert!(!log.contains("possibly shared, from read depth: CFTR"));
}

#[test]
fn a_half_call_is_one_copy() {
    let f = Fixture::new();
    let t = Dir::new("half");
    let v = t.put("half.vcf", &format!("{HEAD}{}", f.hit("CFTR", "./1", 0)));
    let none = t.put("none.vcf", HEAD);
    let ped = t.put("pedigree.tsv", PED);
    let (ca, cb) = (format!("A={v}"), format!("B={none}"));
    let out = t.path("out");
    let (rc, log) = og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &f.genes,
            "--par-bed",
            &f.par,
            "--clinvar",
            &ca,
            "--clinvar",
            &cb,
            "--out-dir",
            &out,
        ],
    );
    assert_eq!(rc, 0, "{log}");
    let r = by_gene(&format!("{out}/couple_A_B.tsv"));
    let fa = cell(&r, "CFTR", "findings_a");
    assert!(
        fa.contains(" one copy, other copy not called ref (./1), "),
        "{fa}"
    );
    // both partners with the half call: one copy each, so 1 in 4, never 1 in 2
    let both = format!("B={v}");
    let out2 = t.path("out2");
    og(
        "screen",
        &[
            "--pedigree",
            &ped,
            "--genes",
            &f.genes,
            "--clinvar",
            &ca,
            "--clinvar",
            &both,
            "--out-dir",
            &out2,
        ],
    );
    let r = by_gene(&format!("{out2}/couple_A_B.tsv"));
    assert_eq!(cell(&r, "CFTR", "offspring_risk"), "1 in 4 per pregnancy");
}

#[test]
fn genes_rebuilds_the_list() {
    let t = Dir::new("genes");
    let par = repo("assets/par_grch38.bed");
    let gtf = t.put(
        "g.gtf",
        &format!(
            "##description: evidence-based annotation of the human genome (GRCh38), version 99\n{}",
            [
                ("NAGA", "chr22", 100, 900),
                ("AAA1", "chr3", 10, 20),
                ("XXL1", "chrX", 5000000, 5000100),
                ("PARG1", "chrX", 20000, 30000),
                ("PARG1", "chrY", 20000, 30000),
                ("YYY1", "chrY", 9000000, 9000100),
                ("C9orf72", "chr9", 1, 2),
            ]
            .iter()
            .map(|(g, c, s, e)| format!(
                "{c}\tHAVANA\tgene\t{s}\t{e}\t.\t+\t.\tgene_id \"X\"; gene_name \"{g}\";\n"
            ))
            .collect::<String>()
        ),
    );
    let p9 = t.put(
        "p9.xml",
        &format!(
            "<JDBOR><DisorderList>{}</DisorderList></JDBOR>",
            [("1", "Autosomal recessive"), ("2", "X-linked recessive"), ("3", "Autosomal dominant")]
                .iter()
                .map(|(c, m)| format!(
                    "<Disorder><OrphaCode>{c}</OrphaCode><TypeOfInheritanceList><TypeOfInheritance><Name>{m}</Name>\
                     </TypeOfInheritance></TypeOfInheritanceList></Disorder>"
                ))
                .collect::<String>()
        ),
    );
    let dis = |code: &str, name: &str, sym: &str, kind: &str| {
        format!(
            "<Disorder><OrphaCode>{code}</OrphaCode><Name>{name}</Name><DisorderGeneAssociationList>\
             <DisorderGeneAssociation><Gene><Symbol>{sym}</Symbol></Gene><DisorderGeneAssociationType><Name>{kind}\
             </Name></DisorderGeneAssociationType><DisorderGeneAssociationStatus><Name>Assessed</Name>\
             </DisorderGeneAssociationStatus></DisorderGeneAssociation></DisorderGeneAssociationList></Disorder>"
        )
    };
    let dc = "Disease-causing germline mutation(s) in";
    let p6 = t.put(
        "p6.xml",
        &format!(
            "<JDBOR date=\"2026-01-02 03:04:05\"><DisorderList>{}{}{}{}{}{}{}</DisorderList></JDBOR>",
            dis("1", "Recessive one", "AAA1", dc),
            dis("2", "X one", "XXL1", dc),
            dis("1", "Pseudo", "PARG1", dc),
            dis("1", "On Y", "YYY1", dc),
            dis("3", "Dominant", "DOM1", dc),
            dis("1", "Candidate", "AAA1", "Candidate gene tested in"),
            dis("1", "Upper case orf", "C9ORF72", dc)
        ),
    );
    let seed = t.put(
        "seed.tsv",
        "gene\tchrom\tstart\tend\tinheritance\tlists\tcondition\nNAGA\t.\t0\t0\tAR\tACMG\tSchindler disease\n\
         OLD1\tchr1\t1\t2\tAR\tOrphanet\tdropped\n",
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
    assert_eq!(g["NAGA"].start, 100, "the ACMG row keeps its GENCODE span");
    assert_eq!(g["NAGA"].condition, "Schindler disease");
    assert_eq!(g["AAA1"].inheritance, "AR");
    assert_eq!(g["XXL1"].inheritance, "XL");
    assert_eq!(g["PARG1"].inheritance, "AR", "a PAR gene is autosomal");
    assert_eq!(g["PARG1"].chrom, "chrX");
    assert!(
        !g.contains_key("YYY1") && !g.contains_key("DOM1"),
        "chrY and dominant left out"
    );
    assert!(
        !g.contains_key("OLD1"),
        "an old Orphanet row is made again, not kept"
    );
    assert!(
        g.contains_key("C9orf72"),
        "an upper-case orf symbol is matched"
    );
    assert_eq!(g["AAA1"].condition, "Recessive one");
    let text = read(&new);
    assert!(
        text.contains("release 99") && text.contains("2026-01-02"),
        "the header names the releases"
    );
    let bad = t.put(
        "bad.tsv",
        "gene\tchrom\tstart\tend\tinheritance\tlists\tcondition\nNOSUCH\t.\t0\t0\tAR\tACMG\tx\n",
    );
    let bad_out = t.path("bad_out.tsv");
    let (rc, log) = og(
        "genes",
        &[
            "--table",
            &bad,
            "--orphanet-genes",
            &p6,
            "--orphanet-inheritance",
            &p9,
            "--gtf",
            &gtf,
            "--par-bed",
            &par,
            "--out",
            &bad_out,
        ],
    );
    assert!(
        rc == 2 && log.contains("NOSUCH"),
        "a panel gene GENCODE does not have stops it: {rc} {log}"
    );
}
