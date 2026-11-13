//! `og bloodtype report` end to end on the repository's marker table and
//! catalog, with synthetic VCFs and depth tables.

use og_bloodtype::genotypes::{Access, read_vcf};
use og_bloodtype::markers::{read_markers, rhd_regions, rhd_row};
use og_bloodtype::{Files, check_markers, report};
use og_cn::catalog::Catalog;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

fn markers() -> String {
    format!("{REPO}/assets/bloodtype_markers.tsv")
}

fn cn_loci() -> String {
    format!("{REPO}/bin/cn_loci.tsv")
}

/// A fresh folder for one test.
fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-bloodtype-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

const HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

/// A one-sample VCF record at the marker `name`.
fn rec(name: &str, gt: &str) -> String {
    let m = &read_markers(&markers()).unwrap()[name];
    format!(
        "{}\t{}\t.\t{}\t{}\t50\tPASS\t.\tGT\t{gt}\n",
        m.chrom, m.pos, m.ref_, m.alt
    )
}

/// A bedcov table: every site at 30x, RHD at `copies` with exon `e`
/// at `exons[e]` copies.
fn depth(d: &Path, copies: f64, exons: &HashMap<i64, f64>) -> String {
    let cat = Catalog::load(&cn_loci()).unwrap();
    let mut text = String::new();
    for m in read_markers(&markers()).unwrap().values() {
        let len = m.pos - m.start;
        text.push_str(&format!(
            "{}\t{}\t{}\tsite\t{}\t{}\t{}\n",
            m.chrom,
            m.start,
            m.pos,
            m.system,
            m.name,
            30 * len
        ));
    }
    for (name, r) in rhd_regions(rhd_row(&cat).unwrap()) {
        let c = match name.as_str() {
            "RHD_FLANK_UP" | "RHD_FLANK_DOWN" => 2.0,
            "RHD_GENE" => copies,
            n => exons
                .get(&n["RHD_EX".len()..].parse::<i64>().unwrap())
                .copied()
                .unwrap_or(copies),
        };
        let mean = 15.0 * c;
        text.push_str(&format!(
            "{}\t{}\t{}\t{name}\t{}\n",
            r.chrom,
            r.start,
            r.end,
            mean * (r.end - r.start) as f64
        ));
    }
    let p = d.join("d.tsv");
    std::fs::write(&p, text).unwrap();
    p.to_string_lossy().into_owned()
}

/// {system: (phenotype, genotype, confidence, note)} of one report.
fn run(d: &Path, records: &[String], depth: Option<String>) -> HashMap<String, Vec<String>> {
    let vcf = d.join("s.vcf");
    std::fs::write(&vcf, format!("{HEAD}{}", records.concat())).unwrap();
    let files = Files {
        vcf: vcf.to_string_lossy().into_owned(),
        markers: markers(),
        depth,
        cn_loci: Some(cn_loci()),
        out: d.join("t.tsv").to_string_lossy().into_owned(),
        evidence: d.join("e.tsv").to_string_lossy().into_owned(),
    };
    report(&files).unwrap();
    std::fs::read_to_string(&files.out)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<String> = l.split('\t').map(str::to_owned).collect();
            (c[0].clone(), c[1..].to_vec())
        })
        .collect()
}

#[test]
fn the_shipped_table_agrees_with_the_rules() {
    assert_eq!(
        check_markers(&markers(), Some(&cn_loci())).unwrap(),
        Vec::<String>::new()
    );
}

#[test]
fn a_reference_genome_with_reads() {
    let d = dir("ref");
    let t = run(&d, &[], Some(depth(&d, 2.0, &HashMap::new())));
    for (system, want) in [
        ("ABO", "O"),
        ("RhD", "D+"),
        ("RhCE", "C-c+ E-e+"),
        ("Kell", "K-k+ Kp(a-b+) Js(a-b+)"),
        ("Duffy", "Fy(a+b-)"),
        ("Kidd", "Jk(a+b-)"),
        ("HPA-15", "HPA-15(a-b+)"),
        ("H", "H+"),
    ] {
        assert_eq!(t[system][0], want, "{system}");
    }
    assert_eq!(t["ABO"][2], "high");
    assert_eq!(t["RhCE"][2], "low");
}

#[test]
fn abo_a1b_and_two_kidd_nulls() {
    let d = dir("abo");
    let recs = [
        rec("ABO_261", "1/1"),
        rec("ABO_803", "0/1"),
        rec("ABO_796", "0/1"),
        rec("JK_01N06", "1/1"),
    ];
    let t = run(&d, &recs, Some(depth(&d, 2.0, &HashMap::new())));
    assert_eq!(t["ABO"][..2], ["A1B", "A1/B"]);
    assert_eq!(t["Kidd"][..2], ["Jk(a-b-)", "null/null"]);
}

#[test]
fn rhd_hybrid_needs_an_exon_under_half_a_copy() {
    let d = dir("hybrid");
    let exons = |c: f64| (3..=7).map(|e| (e, c)).collect::<HashMap<i64, f64>>();
    let t = run(&d, &[], Some(depth(&d, 1.6, &exons(0.45))));
    assert_eq!(t["RhD"][0], "not determined");
    assert!(t["RhD"][3].contains("exons 3, 4, 5, 6, 7 without RHD reads"));
    let t = run(&d, &[], Some(depth(&d, 1.6, &exons(0.55))));
    assert_eq!(t["RhD"][0], "D+");
}

#[test]
fn rhd_deleted_on_both_chromosomes() {
    let d = dir("rhd0");
    let t = run(&d, &[], Some(depth(&d, 0.0, &HashMap::new())));
    assert_eq!(
        t["RhD"][..3],
        [
            "D-",
            "RHD deleted/RHD deleted (RHD copies by read depth 0.00)",
            "high"
        ]
    );
}

/// Copies `from` to `to` as a new file (a fresh mtime, unlike `fs::copy`).
fn put(from: PathBuf, to: impl AsRef<Path>) {
    std::fs::write(to, std::fs::read(from).unwrap()).unwrap();
}

/// The parity inputs: one VCF bgzipped with its .tbi, and the same VCF plain.
fn parity(name: &str) -> PathBuf {
    PathBuf::from(format!("{REPO}/tests/parity/inputs/bloodtype/{name}"))
}

#[test]
fn the_index_and_a_full_pass_read_the_same_sites() {
    for stem in ["abo", "blocks", "gvcf", "vcf-quirks", "nochr", "bombay"] {
        let d = dir(&format!("idx-{stem}"));
        let gz = d.join(format!("{stem}.vcf.gz"));
        put(parity(&format!("{stem}.vcf.gz")), &gz);
        put(
            parity(&format!("{stem}.vcf.gz.tbi")),
            d.join(format!("{stem}.vcf.gz.tbi")),
        );
        let (indexed, how) =
            read_vcf(&gz.to_string_lossy(), read_markers(&markers()).unwrap()).unwrap();
        assert_eq!(how, Access::Index, "{stem}");
        let plain = parity(&format!("{stem}.vcf"));
        let (whole, how) =
            read_vcf(&plain.to_string_lossy(), read_markers(&markers()).unwrap()).unwrap();
        assert_eq!(how, Access::Stream, "{stem}");
        assert_eq!(indexed, whole, "{stem}");
        let called = whole
            .values()
            .filter(|s| s.status != og_bloodtype::genotypes::Status::Absent)
            .count();
        assert!(called > 0, "{stem}: no site read");
    }
}

#[test]
fn an_index_older_than_its_vcf_is_not_used() {
    let d = dir("stale");
    let gz = d.join("abo.vcf.gz");
    put(parity("abo.vcf.gz.tbi"), d.join("abo.vcf.gz.tbi"));
    std::thread::sleep(std::time::Duration::from_millis(20));
    put(parity("abo.vcf.gz"), &gz);
    let (_, how) = read_vcf(&gz.to_string_lossy(), read_markers(&markers()).unwrap()).unwrap();
    assert_eq!(how, Access::Stream);
}
