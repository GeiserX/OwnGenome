//! Step 42's readings say no more than their studies found. These cases pin:
//!
//! 1. the Ethiopian lactase persistence allele -14009G (rs869051967, MCM6,
//!    Jones et al. 2013) has its row, at dbSNP's GRCh38 position with the
//!    plus-strand bases, and `og traits table` reads it;
//! 2. no copy of -13910T never reads as proof of lactose intolerance, the
//!    lactase notes say the rows hold the main known alleles and not every
//!    one, and the weaker -22018 marker does not count as a persistence
//!    allele;
//! 3. the readings that said more than their studies are gone: ALDH2
//!    "cleared normally", MC1R tanning from one copy, "defines e4/e2", the
//!    Gilbert "carrier", CYP2A6 nicotine clearance from one site, and the
//!    unsupported GC comparison;
//! 4. the GC rule's protein levels are each quoted as that study's group
//!    mean, never as a level of your own;
//! 5. the docs of the genotype-rules table point a script at the APOL1
//!    result with G2 in `summary.json`, since the table's `result` column
//!    never shows G2.

use std::collections::HashMap;
use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

fn repo(rel: &str) -> String {
    format!("{REPO}/{rel}")
}

fn read(rel: &str) -> String {
    std::fs::read_to_string(repo(rel)).unwrap_or_else(|e| panic!("{rel}: {e}"))
}

/// The traits table as maps, header names as keys, one per data row.
fn traits() -> Vec<HashMap<String, String>> {
    let text = read("assets/traits.tsv");
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            header
                .iter()
                .map(|h| (*h).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect()
}

fn row(rows: &[HashMap<String, String>], rsid: &str) -> HashMap<String, String> {
    rows.iter()
        .find(|r| r["rsid"] == rsid)
        .unwrap_or_else(|| panic!("assets/traits.tsv has no {rsid} row"))
        .clone()
}

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!(
            "og-i411-traits-readings-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn path(&self, rel: &str) -> String {
        self.0.join(rel).to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

const LACTASE: [&str; 6] = [
    "rs4988235",
    "rs182549",
    "rs145946881",
    "rs41380347",
    "rs41525747",
    "rs869051967",
];

#[test]
fn the_ethiopian_lactase_allele_has_its_row() {
    let r = row(&traits(), "rs869051967");
    // dbSNP's GRCh38 record: NC_000002.12:g.135851175A>C, -14009T>G on the
    // LCT strand (the minus strand).
    assert_eq!(
        [&r["chrom"], &r["pos"], &r["gene"], &r["ref"], &r["alt"]],
        ["chr2", "135851175", "MCM6", "A", "C"]
    );
    assert_eq!(r["ancestry"], "yes");
    assert!(r["pmids"].split(',').any(|p| p == "23993196"), "{r:?}");
    assert!(r["trait"].starts_with("Lactase persistence"), "{r:?}");
    assert!(r["het"].contains("-14009G"), "{r:?}");
    assert!(r["note"].contains("Ethiopia"), "{r:?}");
}

#[test]
fn og_traits_table_reads_the_ethiopian_allele() {
    let w = Scratch::new("table");
    std::fs::create_dir_all(w.path("traits")).unwrap();
    let g = w.path("g.tsv");
    std::fs::write(&g, "chr2\t135851175\tA\tC\tA/C\t40\n").unwrap();
    let (out, rout) = (
        w.path("traits/S_traits.tsv"),
        w.path("traits/S_genotype_rules.tsv"),
    );
    let args: Vec<String> = [
        "table",
        "--table",
        &repo("assets/traits.tsv"),
        "--rules",
        &repo("assets/genotype_rules.tsv"),
        "--genotypes",
        &g,
        "--sample",
        "S",
        "--out",
        &out,
        "--rules-out",
        &rout,
    ]
    .iter()
    .map(|a| (*a).to_owned())
    .collect();
    assert_eq!(og_steps::traits::cli(&args), 0);
    let text = std::fs::read_to_string(&out).unwrap();
    let line = text
        .lines()
        .find(|l| l.starts_with("rs869051967\t"))
        .unwrap_or_else(|| panic!("no rs869051967 line in {text}"));
    let c: Vec<&str> = line.split('\t').collect();
    let want = row(&traits(), "rs869051967");
    assert_eq!(
        (c[5], c[6], c[7]),
        ("A/C", "called", want["het"].as_str()),
        "{line}"
    );
}

#[test]
fn no_persistence_allele_here_is_no_proof_of_intolerance() {
    let rows = traits();
    let lct = row(&rows, "rs4988235");
    let hom_ref = &lct["hom_ref"];
    assert!(hom_ref.contains("does not prove"), "{hom_ref}");
    assert!(!hom_ref.contains("unless another row"), "{hom_ref}");
    for rsid in LACTASE {
        let note = &row(&rows, rsid)["note"];
        assert!(!note.contains("five"), "{rsid}: {note}");
        assert!(note.contains("not every one"), "{rsid}: {note}");
        assert!(note.contains("-14009"), "{rsid}: {note}");
    }
    let tag = row(&rows, "rs182549");
    for k in ["het", "hom_alt"] {
        assert!(tag[k].contains("not a persistence allele"), "{k}: {tag:?}");
    }
}

#[test]
fn readings_that_said_more_than_their_studies_are_gone() {
    let rows = traits();
    let gone: [(&str, &str, &str); 11] = [
        ("rs671", "hom_ref", "cleared normally"),
        ("rs1805007", "het", "tans poorly is more likely"),
        ("rs1805008", "het", "tans poorly is more likely"),
        ("rs1805009", "het", "tans poorly is more likely"),
        ("rs1805007", "note", "necessary but not always sufficient"),
        ("rs429358", "het", "defines e4"),
        ("rs7412", "het", "defines e2"),
        ("rs887829", "het", "a carrier"),
        ("rs887829", "hom_alt", "the genotype of Gilbert syndrome ("),
        ("rs1801272", "hom_alt", "nicotine is cleared slowly"),
        ("rs2282679", "note", "matter more"),
    ];
    for (rsid, k, phrase) in gone {
        let text = &row(&rows, rsid)[k];
        assert!(!text.contains(phrase), "{rsid} {k}: {text}");
    }
    assert!(row(&rows, "rs671")["hom_ref"].contains("not read here"));
    let fut2 = row(&rows, "rs601338");
    for k in ["hom_ref", "het"] {
        assert!(fut2[k].contains("not read here"), "{k}: {fut2:?}");
    }
    assert!(row(&rows, "rs2282679")["hom_alt"].contains("Framingham"));
}

#[test]
fn gc_protein_levels_are_the_study_means_not_your_own() {
    let text = read("assets/genotype_rules.tsv");
    let mut n = 0;
    for l in text.lines() {
        let c: Vec<&str> = l.split('\t').collect();
        if c.first() == Some(&"combination") && c.get(1) == Some(&"GC") {
            let reading = c[4];
            if reading.contains("micrograms per millilitre") {
                n += 1;
                assert!(
                    reading.contains("not an estimate of your own level"),
                    "{reading}"
                );
            }
        }
    }
    assert_eq!(n, 4, "the four GC pairs with a level of Powe 2013");
}

#[test]
fn the_genotype_rules_docs_point_scripts_at_the_g2_result() {
    for doc in ["docs/output-contract.md", "docs/42-traits.md"] {
        let text = read(doc);
        let line = text
            .lines()
            .find(|l| l.contains("sections.traits.data.rules"))
            .unwrap_or_else(|| panic!("{doc} does not name sections.traits.data.rules"));
        for want in ["_genotype_rules.tsv", "APOL1", "`g2`"] {
            assert!(line.contains(want), "{doc}: {want} not in: {line}");
        }
    }
}
