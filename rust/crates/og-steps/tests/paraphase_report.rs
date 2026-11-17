//! og paraphase: step 35b's read-group check and its two tables (the checks
//! of the Python test the port replaced, tests/test_paraphase_report.py;
//! its section 4, SMN1 from step 35's Parascopy table beside Paraphase, is
//! og-couple's `paralog_readers` test).
//!
//!  1. check-bam: a BAM whose read groups all say PL:PACBIO passes (exit 0);
//!     one with PL:ILLUMINA, one with no @RG line, one with a PacBio and an
//!     Oxford Nanopore read group, and a file that is not a BAM stop (exit
//!     3, or 1 for the file that cannot be read), each saying why;
//!  2. table, on a JSON shaped like Paraphase 4.1.0's: SMN1, SMN2 and the
//!     SMN exon 7-8 deletion, PMS2, STRC and NCF1 from gene_cn, HBA1 and
//!     HBA2 from the haplotype names, RCCX, GBA1 and OTOA from their gene1
//!     haplotypes, a region without gene labels as the whole family;
//!  3. the variant table: one row per record a haplotype carries;
//!  5. the region catalog, assets/paraphase_genes.tsv: its rows, the old
//!     eight regions byte for byte, every new entry, a region without its
//!     field not assessed (except HBA's genotype and RCCX's alleles), a row
//!     reading either field on another region refused;
//!  6. the region count: docs/35-paralogs.md and scripts/35b-paraphase.sh
//!     state PARAPHASE_REGIONS, the regions of the pinned Paraphase 4.1.0's
//!     data/38/config.yaml (the e2e case paralog-rows-1 reads the constant
//!     below and counts them in the image);
//!  7. og's own fixes: a header announcing a negative or huge text length,
//!     a JSON object where a copy number belongs, and no copy-number table
//!     left behind when the variant table cannot be written.

use pycompat::json::{Dump, Value, dumps, loads};
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

/// The regions of the pinned Paraphase 4.1.0's data/38/config.yaml (its
/// top-level keys); the e2e case paralog-rows-1 counts them in PARAPHASE_IMAGE.
const PARAPHASE_REGIONS: usize = 165;

fn repo(rel: &str) -> String {
    format!("{}/../../../{rel}", env!("CARGO_MANIFEST_DIR"))
}

struct Dir(PathBuf);

impl Dir {
    fn new(name: &str) -> Dir {
        let p = std::env::temp_dir().join(format!("og-paraphase-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&p);
        std::fs::create_dir_all(&p).unwrap();
        Dir(p)
    }

    fn path(&self, rel: &str) -> String {
        self.0.join(rel).to_string_lossy().into_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// `og paraphase ARGS`: the exit code and stdout followed by stderr.
fn run(args: &[&str]) -> (u8, String) {
    let argv: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let rc = og_steps::paraphase::cli(&argv, &mut out, &mut err);
    out.extend(err);
    (rc, String::from_utf8_lossy(&out).into_owned())
}

/// A BAM with this header and no reference or read (one gzip member, as
/// BGZF); `l_text`: the text length it announces, when not the real one.
fn bam(path: &str, header: &str, l_text: Option<i32>) {
    let text = header.as_bytes();
    let n = l_text.unwrap_or(i32::try_from(text.len()).unwrap());
    let mut raw = b"BAM\x01".to_vec();
    raw.extend(n.to_le_bytes());
    raw.extend(text);
    raw.extend(0i32.to_le_bytes());
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(&raw).unwrap();
    std::fs::write(path, gz.finish().unwrap()).unwrap();
}

const PARAPHASE_JSON: &str = r#"{
 "smn1": {"total_cn": 4, "genome_depth": 31.5, "region_depth": {"median": 120, "percentile80": 130},
          "region_specific_info": {"smn1_cn": 2, "smn2_cn": null, "smn_del78_cn": 0}},
 "hba": {"total_cn": 3, "genome_depth": 31.5, "region_depth": {"median": 90},
         "final_haplotypes": {"1121": "hba_hba1hap1", "1211": "hba_hba1hap2", "2111": "hba_hba2hap1",
                              "2221": "hba_4p2delhap1", "1112": "hba_homology_hap1"},
         "two_copy_haplotypes": [],
         "region_specific_info": {"genotype": "-a/aa",
                                  "sv_called": {"hba_4p2delhap1": "4p2del,chr16:169742-174000"}}},
 "rccx": {"total_cn": 4, "region_specific_info": {"annotated_alleles": ["WT", "WT"], "phasing_success": true}},
 "pms2": {"total_cn": 4, "region_specific_info": {"gene_cn": 2}},
 "strc": {"total_cn": 4, "region_specific_info": {"gene_cn": 2, "intergenic_depth": 30}},
 "ncf1": {"total_cn": 6, "region_specific_info": {"gene_cn": 2, "gene_reads": 20, "pseudo_reads": 40}},
 "GBA": {"total_cn": 4, "genes_in_region": "GBA1",
         "final_haplotypes": {"1": "GBA_gene1hap1", "2": "GBA_gene1hap2", "3": "GBA_gene2hap1", "4": "GBA_gene2hap2"},
         "two_copy_haplotypes": [], "region_specific_info": {"fusions_called": {}}},
 "OTOA": {"total_cn": 4, "genes_in_region": "OTOA",
          "final_haplotypes": {"1": "OTOA_gene1hap1", "2": "OTOA_gene2hap1", "3": "OTOA_unknownhap1"},
          "two_copy_haplotypes": ["OTOA_gene1hap1"], "region_specific_info": {}},
 "CFC1": {"total_cn": 3, "genes_in_region": "CFC1,CFC1B",
          "final_haplotypes": {"1": "CFC1_hap1", "2": "CFC1_hap2", "3": "CFC1_hap3"}, "region_specific_info": {}},
 "TNXB": {"total_cn": null, "genes_in_region": ["TNXB"], "region_specific_info": {}}
}"#;

const VCF: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpms2_pms2hap1\tpms2_pms2hap2\tpms2_pms2clhap1\n\
chr7\t5970000\t.\tC\tT\t.\tPASS\tHPBOUND=1-2\tGT:DP:AD\t1:20:0,19\t0:22:22,0\t1:18:1,17\n\
chr7\t5980000\t.\tG\tA\t.\tPASS\tHPBOUND=1-2\tGT:DP:AD\t0:20:20,0\t0:22:22,0\t.:3:2,1\n";

/// The table the hard-coded regions wrote for PARAPHASE_JSON before
/// assets/paraphase_genes.tsv (bin/paraphase_report.py at main 71d5d58d).
const GOLDEN_OLD8: &str = "region\tgene\tcopy_number\tstatus\tdetail\ttotal_cn\tregion_depth\tgenome_depth\n\
hba\tHBA1\t2\tcalled\tgenotype -a/aa (3 alpha-globin genes); 4.2 kb deletion (-a4.2): HBA2 lost (hba_4p2delhap1)\t3\t90\t31.5\n\
hba\tHBA2\t1\tcalled\tgenotype -a/aa (3 alpha-globin genes); 4.2 kb deletion (-a4.2): HBA2 lost (hba_4p2delhap1)\t3\t90\t31.5\n\
ncf1\tNCF1\t2\tcalled\t.\t6\t.\t.\n\
pms2\tPMS2\t2\tcalled\t.\t4\t.\t.\n\
rccx\tCYP21A2 (RCCX)\t4\tcalled\tRCCX module copies (each carries CYP21A2 or its pseudogene CYP21A1P); CYP21A2 alleles: WT / WT\t4\t.\t.\n\
smn1\tSMN1\t2\tcalled\t.\t4\t120\t31.5\n\
smn1\tSMN2\t.\tno call\tone haplotype seen and the depth does not settle one copy or two\t4\t120\t31.5\n\
smn1\tSMN exon 7-8 deletion\t0\tcalled\t.\t4\t120\t31.5\n\
strc\tSTRC\t2\tcalled\t.\t4\t.\t.\n\
CFC1\tCFC1,CFC1B\t3\tcalled\tcopies of the whole gene family (gene and paralogs together)\t3\t.\t.\n\
GBA\tGBA1\t2\tcalled\tcopies of GBA1 alone (Paraphase's gene1 haplotypes); pseudogene copies 2 (gene2 haplotypes)\t4\t.\t.\n\
OTOA\tOTOA\t2\tcalled\tcopies of OTOA alone (Paraphase's gene1 haplotypes); pseudogene copies 1 (gene2 haplotypes); \
1 haplotype(s) not placed on either\t4\t.\t.\n\
TNXB\tTNXB\t.\tno call\tcopies of the whole gene family (gene and paralogs together)\t.\t.\t.\n";

/// The rows assets/paraphase_genes.tsv must hold, in its order.
const EXPECTED_ROWS: [(&str, &str); 20] = [
    ("smn1", "SMN1"),
    ("smn1", "SMN2"),
    ("smn1", "SMN exon 7-8 deletion"),
    ("hba", "HBA1"),
    ("hba", "HBA2"),
    ("rccx", "CYP21A2 (RCCX)"),
    ("GBA", "GBA1"),
    ("pms2", "PMS2"),
    ("strc", "STRC"),
    ("ncf1", "NCF1"),
    ("OTOA", "OTOA"),
    ("f8", "F8 intron 22"),
    ("opn1lw", "OPN1LW"),
    ("opn1lw", "OPN1MW"),
    ("ikbkg", "IKBKG"),
    ("neb", "NEB triplicate repeat"),
    ("CYP11B1", "CYP11B1,CYP11B2"),
    ("TPSAB1", "TPSAB1,TPSB2"),
    ("AMY1A", "AMY1A,AMY1B,AMY1C"),
    ("CFHclust", "CFH,CFHR1,CFHR2,CFHR3,CFHR4"),
];

type Want = (&'static str, &'static str, &'static str, &'static str);

/// Synthetic regions shaped as Paraphase 4.1.0 writes them (its
/// src/genes/*.rs, src/toolkit/update_calls.rs and docs/*.md), and the rows
/// each must give: (gene, copy_number, status, a piece of the detail).
const NEW_REGIONS: [(&str, &str, &[Want]); 8] = [
    (
        "f8",
        r#"{"total_cn": 3, "region_specific_info": {"sv_called": {"f8_int22invhap1": "inversion"}, "exon1_to_exon22_depth": 30}}"#,
        &[("F8 intron 22", "3", "called", "inversion (f8_int22invhap1)")],
    ),
    (
        "opn1lw",
        r#"{"total_cn": 2, "region_specific_info": {"opn1lw_cn": 1, "opn1mw_cn": 1, "phasing_success": true}}"#,
        &[("OPN1LW", "1", "called", ""), ("OPN1MW", "1", "called", "")],
    ),
    (
        "ikbkg",
        r#"{"total_cn": 4, "region_specific_info": {"gene_cn": 2, "deletion_haplotypes": ["ikbkg_pseudohap1"], "del_read_number": 12}}"#,
        &[("IKBKG", "2", "called", "11.7 kb deletion: ikbkg_pseudohap1")],
    ),
    (
        "neb",
        r#"{"total_cn": 6, "region_specific_info": {"repeat_name": {"tri1": ["neb_hap1", "neb_hap2"], "tri2": ["neb_hap3", "neb_hap4"], "tri3": ["neb_hap5", "neb_hap6"]}}}"#,
        &[(
            "NEB triplicate repeat",
            "6",
            "called",
            "TRI1: 2 haplotype(s)",
        )],
    ),
    (
        "CYP11B1",
        r#"{"total_cn": 3, "region_specific_info": {"fusions_called": {"CYP11B1_hap3": {"type": "deletion", "breakpoint": [[1, 2], [3, 4]], "sequence": "12"}}}}"#,
        &[("CYP11B1,CYP11B2", "3", "called", "deletion (CYP11B1_hap3)")],
    ),
    (
        "TPSAB1",
        r#"{"total_cn": 5, "genes_in_region": "TPSAB1,TPSB2", "region_specific_info": {}}"#,
        &[("TPSAB1,TPSB2", "5", "called", "whole gene family")],
    ),
    (
        "AMY1A",
        r#"{"total_cn": 6, "genes_in_region": "AMY1A,AMY1B,AMY1C", "region_specific_info": {}}"#,
        &[("AMY1A,AMY1B,AMY1C", "6", "called", "whole gene family")],
    ),
    (
        "CFHclust",
        r#"{"total_cn": 4, "genes_in_region": "CFH,CFHR1,CFHR2,CFHR3,CFHR4", "region_specific_info": {"fusions_called": {}}}"#,
        &[(
            "CFH,CFHR1,CFHR2,CFHR3,CFHR4",
            "4",
            "called",
            "no deletion or duplication",
        )],
    ),
];

fn json(text: &str) -> Value {
    loads(text).unwrap()
}

fn obj(v: &mut Value) -> &mut indexmap::IndexMap<String, Value> {
    match v {
        Value::Object(m) => m,
        _ => panic!("{v:?}"),
    }
}

/// Sets `v[path...] = x`.
fn set(v: &mut Value, path: &[&str], x: Value) {
    let (last, head) = path.split_last().unwrap();
    let mut cur = v;
    for k in head {
        cur = obj(cur).get_mut(*k).unwrap();
    }
    obj(cur).insert((*last).to_owned(), x);
}

/// The copy-number rows og paraphase table writes for this JSON, or one
/// ["error", stderr] row.
fn cn_table(data: &Value, genes: &str) -> Vec<Vec<String>> {
    let t = Dir::new(&format!("cn-{:x}", data_hash(data)));
    let (j, cn, var) = (t.path("p.json"), t.path("cn.tsv"), t.path("var.tsv"));
    std::fs::write(&j, dumps(data, &Dump::default()).unwrap()).unwrap();
    let (rc, out) = run(&[
        "table",
        "--sample",
        "S",
        "--json",
        &j,
        "--out-cn",
        &cn,
        "--out-variants",
        &var,
        "--genes",
        genes,
    ]);
    if rc != 0 {
        return vec![vec!["error".into(), out]];
    }
    std::fs::read_to_string(&cn)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

/// A name for each JSON's scratch folder, so parallel tests never share one.
fn data_hash(v: &Value) -> u64 {
    use std::hash::{DefaultHasher, Hash, Hasher};
    let mut h = DefaultHasher::new();
    dumps(v, &Dump::default()).unwrap().hash(&mut h);
    std::thread::current().name().unwrap_or("").hash(&mut h);
    h.finish()
}

fn rows(path: &str) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .map(|h| (*h).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect()
}

/// (region, gene) of the catalog, in the file's order.
fn catalog_rows(text: &str) -> Vec<(String, String)> {
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].to_owned(), c[1].to_owned())
        })
        .collect()
}

#[test]
fn check_bam_reads_the_platform() {
    let t = Dir::new("bam");
    let b = t.path("x.bam");
    bam(
        &b,
        "@HD\tVN:1.6\n@RG\tID:a\tPL:PACBIO\tSM:S\n@RG\tID:b\tPL:PacBio\tSM:S\n",
        None,
    );
    let (rc, out) = run(&["check-bam", &b]);
    assert!(
        rc == 0 && out.contains("2 read group(s), all PL:PACBIO"),
        "{out}"
    );
    bam(&b, "@HD\tVN:1.6\n@RG\tID:a\tPL:ILLUMINA\tSM:S\n", None);
    let (rc, out) = run(&["check-bam", &b]);
    assert!(
        rc == 3
            && out.contains("PL ILLUMINA, not PACBIO")
            && out.contains("licence")
            && out.contains("step 35"),
        "{rc} {out}"
    );
    bam(&b, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n", None);
    let (rc, out) = run(&["check-bam", &b]);
    assert!(rc == 3 && out.contains("has no @RG line"), "{rc} {out}");
    bam(
        &b,
        "@RG\tID:a\tPL:PACBIO\n@RG\tID:b\tPL:ONT\n@RG\tID:c\n",
        None,
    );
    let (rc, out) = run(&["check-bam", &b]);
    assert!(
        rc == 3 && out.contains("PL ONT, none, not PACBIO"),
        "{rc} {out}"
    );
    std::fs::write(&b, b"not a bam at all").unwrap();
    let (rc, out) = run(&["check-bam", &b]);
    assert!(rc == 1 && out.contains("ERROR"), "{rc} {out}");
    // 7. a header announcing a negative or huge text length is refused
    for n in [-1, (1 << 30) + 1] {
        bam(&b, "@RG\tID:a\tPL:PACBIO\n", Some(n));
        let (rc, out) = run(&["check-bam", &b]);
        assert!(
            rc == 1 && out.contains("not a BAM header"),
            "{n}: {rc} {out}"
        );
    }
}

#[test]
fn table_reads_each_region() {
    let genes = repo("assets/paraphase_genes.tsv");
    let t = Dir::new("table");
    let j = t.path("S.paraphase.json");
    let vd = t.path("S_paraphase_vcfs");
    std::fs::create_dir_all(&vd).unwrap();
    std::fs::write(&j, PARAPHASE_JSON).unwrap();
    std::fs::write(format!("{vd}/S_pms2.vcf"), VCF).unwrap();
    let (cn, var) = (t.path("cn.tsv"), t.path("var.tsv"));
    let (rc, out) = run(&[
        "table",
        "--sample",
        "S",
        "--json",
        &j,
        "--vcf-dir",
        &vd,
        "--out-cn",
        &cn,
        "--out-variants",
        &var,
        "--genes",
        &genes,
    ]);
    assert_eq!(rc, 0, "{out}");
    let all = rows(&cn);
    let r: HashMap<(String, String), &HashMap<String, String>> = all
        .iter()
        .map(|x| ((x["region"].clone(), x["gene"].clone()), x))
        .collect();
    let g = |reg: &str, gene: &str, k: &str| -> String {
        r.get(&(reg.to_owned(), gene.to_owned()))
            .and_then(|x| x.get(k))
            .cloned()
            .unwrap_or_default()
    };
    let cnv = |reg: &str, gene: &str| g(reg, gene, "copy_number");
    assert_eq!(
        (
            cnv("smn1", "SMN1"),
            cnv("smn1", "SMN2"),
            g("smn1", "SMN2", "status"),
            cnv("smn1", "SMN exon 7-8 deletion")
        ),
        ("2".into(), ".".into(), "no call".into(), "0".into())
    );
    assert!(g("smn1", "SMN2", "detail").contains("depth does not settle one copy or two"));
    assert_eq!(
        (
            cnv("pms2", "PMS2"),
            cnv("strc", "STRC"),
            cnv("ncf1", "NCF1")
        ),
        ("2".into(), "2".into(), "2".into())
    );
    assert_eq!(
        (cnv("hba", "HBA1"), cnv("hba", "HBA2")),
        ("2".into(), "1".into())
    );
    let hd = g("hba", "HBA2", "detail");
    assert!(
        hd.contains("genotype -a/aa (3 alpha-globin genes)")
            && hd.contains("4.2 kb deletion (-a4.2): HBA2 lost")
    );
    assert_eq!(g("smn1", "SMN1", "region_depth"), "120");
    assert_eq!(g("smn1", "SMN1", "genome_depth"), "31.5");
    assert_eq!(cnv("rccx", "CYP21A2 (RCCX)"), "4");
    assert!(g("rccx", "CYP21A2 (RCCX)", "detail").contains("WT / WT"));
    assert_eq!(
        cnv("GBA", "GBA1"),
        "2",
        "from its gene1 haplotypes, not the family's 4"
    );
    assert_eq!(g("GBA", "GBA1", "total_cn"), "4");
    let gd = g("GBA", "GBA1", "detail");
    assert!(
        gd.contains("copies of GBA1 alone") && gd.contains("pseudogene copies 2"),
        "{gd}"
    );
    assert_eq!(cnv("OTOA", "OTOA"), "2", "a two-copy gene1 haplotype twice");
    assert!(g("OTOA", "OTOA", "detail").contains("1 haplotype(s) not placed"));
    assert_eq!(cnv("CFC1", "CFC1,CFC1B"), "3");
    assert!(g("CFC1", "CFC1,CFC1B", "detail").contains("whole gene family"));
    assert_eq!(
        (cnv("TNXB", "TNXB"), g("TNXB", "TNXB", "status")),
        (".".into(), "no call".into())
    );
    let order: Vec<&str> = all.iter().map(|x| x["region"].as_str()).collect();
    let pos = |k: &str| order.iter().position(|x| *x == k).unwrap();
    assert!(
        pos("TNXB") > pos("hba").max(pos("smn1")).max(pos("OTOA")),
        "the named regions come first"
    );
    // 3. the variant table
    let v = rows(&var);
    assert_eq!(v.len(), 1, "{v:?}");
    assert_eq!(v[0]["pos"], "5970000");
    assert_eq!(v[0]["haplotypes"], "pms2_pms2hap1,pms2_pms2clhap1");
    assert_eq!(v[0]["support"], "pms2_pms2hap1:19/20,pms2_pms2clhap1:17/18");
    assert_eq!(v[0]["region"], "pms2");
    // HBA "--/aa" and "aa/aa"
    let mut data = json(PARAPHASE_JSON);
    set(
        &mut data,
        &["hba", "region_specific_info"],
        json(r#"{"genotype": "--/aa", "sv_called": {}}"#),
    );
    set(
        &mut data,
        &["hba", "final_haplotypes"],
        json(r#"{"1": "hba_hba1hap1", "2": "hba_hba2hap1"}"#),
    );
    let tb = cn_table(&data, &genes);
    let row = |tb: &[Vec<String>], gene: &str| {
        tb.iter()
            .find(|x| x[0] == "hba" && x[1] == gene)
            .cloned()
            .unwrap()
    };
    let h2 = row(&tb, "HBA2");
    assert!(
        h2[4].contains("large deletion such as --SEA") && h2[4].contains("(2 alpha-globin genes)"),
        "{h2:?}"
    );
    assert_eq!((row(&tb, "HBA1")[2].as_str(), h2[2].as_str()), ("1", "1"));
    set(
        &mut data,
        &["hba", "region_specific_info"],
        json(r#"{"genotype": "aa/aa", "sv_called": {}}"#),
    );
    set(
        &mut data,
        &["hba", "two_copy_haplotypes"],
        json(r#"["hba_hba1hap1", "hba_hba2hap1"]"#),
    );
    let tb = cn_table(&data, &genes);
    let h2 = row(&tb, "HBA2");
    assert_eq!((row(&tb, "HBA1")[2].as_str(), h2[2].as_str()), ("2", "2"));
    assert!(
        h2[4].contains("no 3.7 kb or 4.2 kb deletion or duplication found"),
        "{h2:?}"
    );
}

#[test]
fn the_region_catalog() {
    let genes = repo("assets/paraphase_genes.tsv");
    let text = std::fs::read_to_string(&genes).unwrap();
    let got = catalog_rows(&text);
    let want: Vec<(String, String)> = EXPECTED_ROWS
        .iter()
        .map(|(a, b)| ((*a).into(), (*b).into()))
        .collect();
    assert_eq!(got, want, "the catalog holds the expected rows, in order");
    let no_neb: String = text
        .lines()
        .filter(|l| !l.starts_with("neb\t"))
        .map(|l| format!("{l}\n"))
        .collect();
    let dropped = catalog_rows(&no_neb);
    assert!(
        dropped != want && dropped.len() == want.len() - 1,
        "control: a dropped NEB row"
    );
    let t = Dir::new("catalog");
    let j = t.path("S.paraphase.json");
    std::fs::write(&j, PARAPHASE_JSON).unwrap();
    let (cn, var) = (t.path("cn.tsv"), t.path("var.tsv"));
    let table = |g: &str| {
        run(&[
            "table",
            "--sample",
            "S",
            "--json",
            &j,
            "--out-cn",
            &cn,
            "--out-variants",
            &var,
            "--genes",
            g,
        ])
    };
    let (rc, out) = table(&genes);
    assert_eq!(rc, 0, "{out}");
    assert_eq!(
        std::fs::read_to_string(&cn).unwrap(),
        GOLDEN_OLD8,
        "the old eight regions, byte for byte"
    );
    let mut moved = json(PARAPHASE_JSON);
    set(
        &mut moved,
        &["pms2", "region_specific_info", "gene_cn"],
        Value::Int(3),
    );
    let joined: String = cn_table(&moved, &genes)
        .iter()
        .map(|x| format!("{}\n", x.join("\t")))
        .collect();
    assert!(
        !GOLDEN_OLD8.contains(&joined),
        "control: one changed field is not the old table"
    );
    let mut new = Value::Object(indexmap::IndexMap::new());
    for (name, region, _) in NEW_REGIONS {
        let mut r = json(region);
        set(&mut r, &["genome_depth"], Value::Int(30));
        set(&mut r, &["region_depth"], json(r#"{"median": 60}"#));
        obj(&mut new).insert(name.to_owned(), r);
    }
    let tb = cn_table(&new, &genes);
    for (name, _, want) in NEW_REGIONS {
        let rs: Vec<&Vec<String>> = tb.iter().filter(|x| x[0] == name).collect();
        assert_eq!(rs.len(), want.len(), "{name}: {rs:?}");
        for (x, w) in rs.iter().zip(want) {
            assert_eq!(
                (x[1].as_str(), x[2].as_str(), x[3].as_str()),
                (w.0, w.1, w.2),
                "{name}"
            );
            assert!(x[4].contains(w.3), "{name}: {x:?}");
        }
    }
    let mut both = Value::Object(indexmap::IndexMap::new());
    let mut tps = json(NEW_REGIONS[5].1);
    set(
        &mut tps,
        &["region_specific_info"],
        json(r#"{"gene_cn": 2}"#),
    );
    obj(&mut both).insert("TPSAB1".into(), tps);
    let rs = cn_table(&both, &genes);
    assert!(
        rs.len() == 1 && rs[0][2] == "5" && rs[0][4].contains("whole gene family"),
        "{rs:?}"
    );
    let first: Vec<&str> = tb.iter().take(4).map(|x| x[0].as_str()).collect();
    assert_eq!(
        first,
        ["CFHclust", "CYP11B1", "f8", "ikbkg"],
        "the new regions with a field of their own first"
    );
    for name in ["f8", "opn1lw", "ikbkg", "neb", "CYP11B1", "CFHclust"] {
        let (_, region, _) = NEW_REGIONS.iter().find(|r| r.0 == name).unwrap();
        let mut r = json(region);
        set(&mut r, &["region_specific_info"], json("{}"));
        let mut bare = Value::Object(indexmap::IndexMap::new());
        obj(&mut bare).insert(name.into(), r);
        let rs = cn_table(&bare, &genes);
        assert!(
            !rs.is_empty()
                && rs.iter().all(|x| x[2] == "."
                    && x[3] == "not assessed"
                    && x[4].contains("Paraphase wrote no")),
            "{name}: {rs:?}"
        );
    }
    let bad = t.path("genes.tsv");
    std::fs::write(&bad, "region\tgene\tfield\n").unwrap();
    let (rc, out) = table(&bad);
    assert!(
        rc == 1 && out.contains("header"),
        "a catalog with the wrong columns: {rc} {out}"
    );
    // genotype is read by HBA's code and annotated_alleles by RCCX's: on
    // another region the catalog is refused.
    let head = "region\tgene\tfield\treading\tsource\n";
    for (region, field) in [("smn1", "genotype"), ("pms2", "annotated_alleles")] {
        let p = t.path(&format!("genes-{field}.tsv"));
        std::fs::write(&p, format!("{head}{region}\tX\t{field}\tr\ts\n")).unwrap();
        let (rc, out) = table(&p);
        assert!(
            rc == 1 && out.contains(&format!("field {field} on region {region}")),
            "{rc} {out}"
        );
    }
    let own = t.path("genes-own.tsv");
    std::fs::write(
        &own,
        format!("{head}hba\tHBA1\tgenotype\tr\ts\nrccx\tCYP21A2 (RCCX)\tannotated_alleles\tr\ts\n"),
    )
    .unwrap();
    let (rc, out) = table(&own);
    assert_eq!(
        rc, 0,
        "control: the same fields on their own regions load: {out}"
    );
    // HBA's genotype and RCCX's alleles: the copy number does not come from
    // the field, so the rows keep it.
    let base = json(PARAPHASE_JSON);
    let mut bare = Value::Object(indexmap::IndexMap::new());
    for k in ["hba", "rccx"] {
        let Value::Object(m) = &base else {
            unreachable!()
        };
        let mut r = m[k].clone();
        set(&mut r, &["region_specific_info"], json("{}"));
        obj(&mut bare).insert(k.into(), r);
    }
    let got: HashMap<String, Vec<String>> = cn_table(&bare, &genes)
        .into_iter()
        .map(|x| (x[1].clone(), x))
        .collect();
    assert_eq!(
        (got["HBA1"][2].as_str(), got["HBA1"][3].as_str()),
        ("2", "no call"),
        "{got:?}"
    );
    assert_eq!(
        (got["HBA2"][2].as_str(), got["HBA2"][3].as_str()),
        ("1", "no call")
    );
    assert!(got["HBA1"][4].contains("could not resolve the alleles"));
    let rc_ = &got["CYP21A2 (RCCX)"];
    assert_eq!((rc_[2].as_str(), rc_[3].as_str()), ("4", "called"));
    assert!(!rc_[4].contains("alleles:"), "{rc_:?}");
}

/// Every count of Paraphase regions a text states ("162 GRCh38 regions",
/// "all 162", "Paraphase's 162").
fn region_counts(text: &str) -> Vec<usize> {
    let mut out = Vec::new();
    let num_before = |s: &str| -> Option<usize> {
        let d: String = s.chars().rev().take_while(char::is_ascii_digit).collect();
        d.chars().rev().collect::<String>().parse().ok()
    };
    let num_after = |s: &str| -> Option<usize> {
        let d: String = s.chars().take_while(char::is_ascii_digit).collect();
        let next = s[d.len()..].chars().next();
        (!d.is_empty() && !next.is_some_and(|c| c.is_alphanumeric() || c == '_'))
            .then(|| d.parse().ok())?
    };
    let mut i = 0;
    while i < text.len() {
        let rest = &text[i..];
        if rest.starts_with(" GRCh38 regions")
            && let Some(n) = num_before(&text[..i])
        {
            out.push(n);
        } else if let Some(r) = rest.strip_prefix("all ")
            && let Some(n) = num_after(r)
        {
            out.push(n);
        } else if let Some(r) = rest.strip_prefix("Paraphase's ")
            && let Some(n) = r.chars().next().filter(char::is_ascii_digit).and_then(|_| {
                let d: String = r.chars().take_while(char::is_ascii_digit).collect();
                d.parse().ok()
            })
        {
            out.push(n);
        }
        i += rest.chars().next().map_or(1, char::len_utf8);
    }
    out
}

#[test]
fn the_region_count_in_the_docs_and_the_step() {
    for rel in ["docs/35-paralogs.md", "scripts/35b-paraphase.sh"] {
        let said = region_counts(&std::fs::read_to_string(repo(rel)).unwrap());
        assert!(
            !said.is_empty() && said.iter().all(|n| *n == PARAPHASE_REGIONS),
            "{rel}: {said:?}"
        );
    }
    assert!(
        region_counts("It covers 162 GRCh38 regions.").contains(&162),
        "control: a planted 162"
    );
}

#[test]
fn og_s_own_fixes_in_the_table() {
    let genes = repo("assets/paraphase_genes.tsv");
    let t = Dir::new("fixes");
    let bad = t.path("bad.json");
    std::fs::write(
        &bad,
        r#"{"pms2": {"total_cn": 4, "region_specific_info": {"gene_cn": {"a": 1}}}}"#,
    )
    .unwrap();
    let (cn, var) = (t.path("cn.tsv"), t.path("var.tsv"));
    let (rc, out) = run(&[
        "table",
        "--sample",
        "S",
        "--json",
        &bad,
        "--out-cn",
        &cn,
        "--out-variants",
        &var,
        "--genes",
        &genes,
    ]);
    assert!(
        rc == 1 && out.contains("pms2.gene_cn is not a number"),
        "{rc} {out}"
    );
    let j = t.path("S.paraphase.json");
    std::fs::write(&j, PARAPHASE_JSON).unwrap();
    let _ = std::fs::remove_file(&cn);
    let nowhere = t.path("no-such-dir/var.tsv");
    let (rc, out) = run(&[
        "table",
        "--sample",
        "S",
        "--json",
        &j,
        "--out-cn",
        &cn,
        "--out-variants",
        &nowhere,
        "--genes",
        &genes,
    ]);
    assert!(rc == 1 && !std::path::Path::new(&cn).exists(), "{rc} {out}");
}
