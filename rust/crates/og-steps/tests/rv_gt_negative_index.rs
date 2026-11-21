//! A genotype index below zero names no allele. VCF writes GT as allele
//! numbers from 0 (REF) up, or `.`; the consumer export indexed
//! `[REF] + ALT` as a Python list, so `GT -1/-1` at REF A, ALT T exported
//! `TT`, and a record with a negative index counted as non-reference, so a
//! longer REF with `GT -1/-1` hid the reference block after it as if it
//! were a called deletion. A negative index is now an unreadable genotype.

use std::num::NonZero;
use std::path::{Path, PathBuf};

use og_steps::consumer::export::{ExportArgs, export};
use og_steps::consumer::zip::{DosDate, read_all};

fn input(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/parity/inputs/consumer")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

/// Exports the plain gVCF `text` over the unit site table (rs101 chr20:100
/// A, rs102 :150 C, rs103 :160 G, rs104 :170 T, rs105 :180 A); the 23andMe
/// calls by rsID and the summary table.
fn run(name: &str, text: &str) -> (Vec<(String, String)>, String) {
    let dir: PathBuf = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let gvcf = dir.join("in.g.vcf");
    std::fs::write(&gvcf, text).unwrap();
    let out = dir.join("out");
    let (gvcf, out_dir) = (gvcf.to_string_lossy(), out.to_string_lossy());
    let (sites, par) = (input("unit-sites.tsv.gz"), input("par_grch38.bed"));
    let a = ExportArgs {
        gvcf: &gvcf,
        sites: &sites,
        sample: "S",
        sex: "male",
        par: Some(&par),
        out_dir: &out_dir,
        date: DosDate::parse("2026-01-02").unwrap(),
    };
    export(&a, NonZero::new(2).unwrap()).expect("export");
    let z = std::fs::read(out.join("S_23andme_v5.zip")).unwrap();
    let members = read_all(&z).unwrap();
    let calls = String::from_utf8(members[0].1.clone())
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].to_owned(), c[3].to_owned())
        })
        .collect();
    let tsv = std::fs::read_to_string(out.join("S_consumer_export.tsv")).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    (calls, tsv)
}

#[test]
fn a_negative_genotype_index_is_unreadable() {
    let body = [
        // -1/-1 is no allele: unreadable, never the ALT T.
        "chr20\t100\t.\tA\tT\t30\tPASS\t.\tGT:GQ\t-1/-1:30",
        // Readable genotypes keep their bases.
        "chr20\t150\t.\tC\tT,<*>\t30\tPASS\t.\tGT:GQ\t0/1:30",
        "chr20\t160\t.\tG\tA,<*>\t30\tPASS\t.\tGT:GQ\t1/1:30",
        // A longer REF whose genotype is unreadable is not a called
        // deletion: the reference block after it still answers.
        "chr20\t168\t.\tGTTT\tG,<*>\t30\tPASS\t.\tGT:GQ\t-1/-1:30",
        "chr20\t169\t.\tT\t<*>\t0\t.\tEND=179\tGT:GQ\t0/0:30",
        // One negative part is enough.
        "chr20\t180\t.\tA\tT\t30\tPASS\t.\tGT:GQ\t1/-1:30",
    ];
    let text = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n{}\n",
        body.join("\n")
    );
    let (calls, tsv) = run("rv_gt_negative_index_consumer", &text);
    let call = |id: &str| -> &str {
        &calls
            .iter()
            .find(|(r, _)| r == id)
            .unwrap_or_else(|| panic!("{id} in the export"))
            .1
    };
    // Ids and genotypes apart, so no line pairs an rsID with letters.
    let want = ["--", "CT", "AA", "TT", "--"];
    for (n, g) in (101..=105).zip(want) {
        let id = format!("rs{n}");
        assert_eq!(call(&id), g, "{id}");
    }
    let row = "no_call: unreadable genotype\t2";
    assert!(tsv.lines().any(|l| l == row), "{row} in {tsv}");
}

#[test]
fn a_past_the_end_index_does_not_hide_the_reference_block() {
    let export_call = |label, gt| {
        let text = format!(
            "##fileformat=VCFv4.2\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
             chr20\t168\t.\tGTTT\tG,<*>\t30\tPASS\t.\tGT:GQ\t{gt}:30\n\
             chr20\t169\t.\tT\t<*>\t0\t.\tEND=179\tGT:GQ\t0/0:30\n"
        );
        let (calls, _) = run(label, &text);
        calls
            .into_iter()
            .find(|(id, _)| id == &format!("rs{}", 104))
            .expect("site in the export")
            .1
    };

    assert_eq!(export_call("rv_gt_past_end_consumer", "5/5"), "TT");
    assert_eq!(export_call("rv_gt_one_alt_consumer", "1/5"), "--");
}
