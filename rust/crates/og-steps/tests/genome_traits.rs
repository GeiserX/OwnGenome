//! Step 42 through `og traits table`, and the rule alleles against dbSNP:
//! the checks of the old `tests/test_traits_report.py` and part 6 of
//! `tests/test_lpa_kiv2.py`, by meaning, in cargo test.
//!
//! 1. `scripts/ci/check-rule-alleles.sh` passes on `assets/genotype_rules.tsv`
//!    and names each planted mistake: ref and alt swapped, the minus-strand
//!    bases, an rsID the pinned dbSNP table lacks, a wrong position, a rules
//!    file with no site, and a dbSNP table of another build than versions.env
//!    pins.
//! 2. The APOE row of the traits table for eleven genotype pairs is the
//!    golden row (its sha256), so any change to the APOE reading shows.
//! 3. Each multi-SNP rule (APOL1, SERPINA1, HFE, GC, CAD9p21, MAPT) reads its
//!    combinations; on GRCh38 reference blocks only it reads the label written
//!    here, independently of the rules table; one allele flipped is not the
//!    same result; APOL1 G0/G1 says G2 is not read and the other copy could
//!    carry it.
//! 4. A table with CSV quoting or an embedded tab is refused, and the
//!    outputs of the run before stay as they were.
//! 5. `assets/traits.tsv` holds the two LPA risk variants with the alleles
//!    the Lp(a) flag counts, and the blood-test advice.
//!
//! The reference-block labels are the guard against a rules table that
//! changes a label: the `fault-wrong-reference` feature plants WRONG_REFERENCE
//! as the label of APOL1 G0/G0 and CAD9p21 nonrisk/nonrisk, and the two tests
//! of those rules must then fail (the control in .github/workflows/guard.yml).

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

/// The columns of the traits table (`og traits table --out`).
const OUT_COLUMNS: [&str; 11] = [
    "rsid", "gene", "chrom", "pos", "trait", "genotype", "status", "reading", "ancestry", "note",
    "pmids",
];

/// The sha256 of the APOE row (its OUT_COLUMNS joined by tabs) for each pair
/// of genotypes at rs429358 and rs7412 ("." leaves the site out).
const APOE_GOLDEN: [(&str, &str, &str); 11] = [
    (
        "T/T",
        "C/C",
        "e3d7512b649f1d23f5acbcaa08da19bd6a009d5b39ac57cb6856591c50e7cfda",
    ),
    (
        "T/T",
        "C/T",
        "62e910d8ec7ac45ab15e9c9ed85c93b128ac41703d1149a332082dc0d7b5f389",
    ),
    (
        "T/T",
        "T/T",
        "41b4ad51c67028bb19af1b7ab7e4f414f5c0354e4c4a7c3c19309186b429f960",
    ),
    (
        "T/C",
        "C/C",
        "08c5ed9e2898f2ae29cbd49d5ac97394758fa852f9d35e3c7173b2aa8d3111ad",
    ),
    (
        "T/C",
        "C/T",
        "555c9cc0e8b59b10bdc5bdc2e33d93fb4ccbee1e6ea12916b17e73ffd6fdb8d8",
    ),
    (
        "T/C",
        "T/T",
        "043f51b41486f217cac3a71e705e0e29c1a1ced9745e60abde2230df2fe901db",
    ),
    (
        "C/C",
        "C/C",
        "ca7f36540cf1269b84649b4ad6ace19e0abc03eab7984d495de584540bb11dcc",
    ),
    (
        "C/C",
        "C/T",
        "db499d1bfebedf77ed34a9d0525030d34be6336e685758a8efd53125197e3d00",
    ),
    (
        "C/C",
        "T/T",
        "16e5214f90c5f2c401e5849dffd76a08dc3fc5ba2deb324d7fb02c76a92113be",
    ),
    (
        "T/C",
        ".",
        "34624db16a92981c0882edb3d8ce8d291162aa78b7c5208eba2360141e6ccbd8",
    ),
    (
        "T/G",
        "C/C",
        "1d0aee0b8d241561ebd7477d6f4a8be274a48157a297abecd444c6c6a0a6abac",
    ),
];

/// One reading of a rule: the genotypes at its sites and the result.
type Reading = (&'static [&'static str], &'static str);

/// The multi-SNP rules: genotypes at the rule's sites, in its order, and the
/// result they read.
const MULTI_SNP: [(&str, &[Reading]); 6] = [
    (
        "APOL1",
        &[(&["G/G", "G/G"], "G1/G1"), (&["A/G", "T/G"], "G0/G1")],
    ),
    (
        "SERPINA1",
        &[
            (&["T/T", "T/T"], "PiZZ"),
            (&["T/A", "C/T"], "PiSZ"),
            (&["A/A", "C/C"], "PiSS"),
        ],
    ),
    (
        "HFE",
        &[
            (&["A/A", "C/C"], "C282Y/C282Y"),
            (&["G/A", "C/G"], "C282Y/H63D"),
            (&["G/G", "G/G"], "H63D/H63D"),
        ],
    ),
    (
        "GC",
        &[
            (&["C/C", "G/G"], "Gc1S/Gc1S"),
            (&["A/C", "G/T"], "Gc1S/Gc2"),
            (&["A/A", "T/T"], "Gc2/Gc2"),
        ],
    ),
    (
        "CAD9p21",
        &[
            (&["G/G", "C/C"], "two risk copies"),
            (&["A/G", "G/C"], "one risk copy"),
        ],
    ),
    (
        "MAPT",
        &[(&["G/G", "G/G"], "H2/H2"), (&["A/G", "A/G"], "H1/H2")],
    ),
];

/// The label of homozygous GRCh38 reference blocks at every site of a rule,
/// written here and not read from the rules table.
const MULTI_SNP_REF: [(&str, &str); 6] = [
    ("APOL1", "G0/G0"),
    ("SERPINA1", "PiMM"),
    ("HFE", "wt/wt"),
    ("GC", "Gc1F/Gc1F"),
    ("CAD9p21", "no risk copy"),
    ("MAPT", "H1/H1"),
];

fn repo(rel: &str) -> String {
    format!("{REPO}/{rel}")
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_default()
}

/// A scratch folder of its own, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d =
            std::env::temp_dir().join(format!("og-genome-traits-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn path(&self, rel: &str) -> String {
        self.0.join(rel).to_string_lossy().into_owned()
    }

    fn put(&self, rel: &str, text: &str) -> String {
        let p = self.0.join(rel);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The rules table the tests read: the checkout's, or with the
/// `fault-wrong-reference` feature a copy whose APOL1 G0/G0 and CAD9p21
/// nonrisk/nonrisk labels are WRONG_REFERENCE.
fn rules_path(w: &Scratch) -> String {
    let text = read(&repo("assets/genotype_rules.tsv"));
    if !cfg!(feature = "fault-wrong-reference") {
        return repo("assets/genotype_rules.tsv");
    }
    let mut planted = text;
    for (rid, key, label) in [
        ("APOL1", "G0/G0", "G0/G0"),
        ("CAD9p21", "nonrisk/nonrisk", "no risk copy"),
    ] {
        let old = format!("combination\t{rid}\t{key}\t{label}\t");
        assert!(planted.contains(&old), "{old}");
        planted = planted.replace(
            &old,
            &format!("combination\t{rid}\t{key}\tWRONG_REFERENCE\t"),
        );
    }
    w.put("planted_rules.tsv", &planted)
}

/// The site lines of each rule: rsid, chrom, pos, ref, alt.
fn sites() -> HashMap<String, Vec<[String; 5]>> {
    let mut out: HashMap<String, Vec<[String; 5]>> = HashMap::new();
    for line in read(&repo("assets/genotype_rules.tsv")).lines() {
        let c: Vec<&str> = line.split('\t').collect();
        if c.first() == Some(&"site") && c.len() >= 7 {
            out.entry(c[1].to_owned())
                .or_default()
                .push([c[2], c[3], c[4], c[5], c[6]].map(str::to_owned));
        }
    }
    out
}

/// A TSV with a header, as rows by column name.
fn rows(text: &str) -> Vec<HashMap<String, String>> {
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .map(|h| (*h).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect()
}

type Rows = Vec<HashMap<String, String>>;

/// `og traits table` on genotype records (each a tab-joined line): the
/// traits table and the rules table, read back.
fn table(w: &Scratch, rules: &str, records: &[Vec<String>]) -> (Rows, Rows) {
    let g = w.put(
        "g.tsv",
        &records
            .iter()
            .map(|r| r.join("\t") + "\n")
            .collect::<String>(),
    );
    let (out, rout) = (
        w.path("traits/S_traits.tsv"),
        w.path("traits/S_genotype_rules.tsv"),
    );
    std::fs::create_dir_all(w.path("traits")).unwrap();
    let args: Vec<String> = [
        "table",
        "--table",
        &repo("assets/traits.tsv"),
        "--rules",
        rules,
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
    assert_eq!(og_steps::traits::cli(&args), 0, "og traits table exits 0");
    (rows(&read(&out)), rows(&read(&rout)))
}

fn rule<'a>(lines: &'a Rows, rid: &str) -> &'a HashMap<String, String> {
    lines
        .iter()
        .find(|r| r["rule"] == rid)
        .unwrap_or_else(|| panic!("no {rid} line in {lines:?}"))
}

/// Records of a rule's sites with these genotypes.
fn records(sites: &[[String; 5]], gts: &[&str]) -> Vec<Vec<String>> {
    sites
        .iter()
        .zip(gts)
        .map(|(s, gt)| {
            vec![
                s[1].clone(),
                s[2].clone(),
                s[3].clone(),
                s[4].clone(),
                (*gt).to_owned(),
            ]
        })
        .collect()
}

fn sha256(text: &str) -> String {
    for (cmd, args) in [("sha256sum", &[][..]), ("shasum", &["-a", "256"][..])] {
        let Ok(mut child) = Command::new(cmd)
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
        else {
            continue;
        };
        child
            .stdin
            .take()
            .unwrap()
            .write_all(text.as_bytes())
            .unwrap();
        let o = child.wait_with_output().unwrap();
        if o.status.success() {
            return String::from_utf8(o.stdout)
                .unwrap()
                .split_whitespace()
                .next()
                .unwrap()
                .to_owned();
        }
    }
    panic!("neither sha256sum nor shasum ran");
}

// 1. check-rule-alleles.sh

/// `bash scripts/ci/check-rule-alleles.sh ARGS`: exit code and stdout + stderr.
fn dbsnp(script: &str, args: &[&str]) -> (i32, String) {
    let o = Command::new("bash")
        .arg(script)
        .args(args)
        .output()
        .unwrap();
    (
        o.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&o.stdout).into_owned() + &String::from_utf8_lossy(&o.stderr),
    )
}

#[test]
fn every_rule_site_matches_dbsnp() {
    let script = repo("scripts/ci/check-rule-alleles.sh");
    let n: usize = sites().values().map(Vec::len).sum();
    let (rc, log) = dbsnp(&script, &[&repo("assets/genotype_rules.tsv")]);
    assert!(
        rc == 0 && log.contains(&format!("{n} rule sites checked")),
        "{rc} {log}"
    );
}

#[test]
fn check_rule_alleles_names_each_planted_mistake() {
    let w = Scratch::new("dbsnp");
    let script = repo("scripts/ci/check-rule-alleles.sh");
    let text = read(&repo("assets/genotype_rules.tsv"));
    let site = "site\tTAS2R38\trs1726866\tchr7\t141972905\t";
    assert!(text.contains(&format!("{site}G\tA")));
    let plant = |from: &str, to: &str| {
        let p = w.put("planted.tsv", &text.replace(from, to));
        dbsnp(&script, &[&p])
    };
    for (new, want) in [
        ("A\tG", "rs1726866: ref and alt are swapped"),
        ("C\tT", "rs1726866: the alleles are the minus-strand bases"),
    ] {
        let (rc, log) = plant(&format!("{site}G\tA"), &format!("{site}{new}"));
        assert!(rc == 1 && log.contains(want), "{want}: {rc} {log}");
    }
    let (rc, log) = plant(
        "\trs1726866\tchr7\t141972905\t",
        "\trs1726867\tchr7\t141972905\t",
    );
    assert!(
        rc == 1 && log.contains("rs1726867: not in the pinned dbSNP table"),
        "{rc} {log}"
    );
    let (rc, log) = plant(
        "\trs1726866\tchr7\t141972905\t",
        "\trs1726866\tchr7\t141972906\t",
    );
    assert!(
        rc == 1 && log.contains("dbSNP at chr7:141972905"),
        "{rc} {log}"
    );
    let p = w.put("planted.tsv", "# no sites\n");
    let (rc, log) = dbsnp(&script, &[&p]);
    assert!(rc == 1 && log.contains("no site line"), "{rc} {log}");
    // A copy of the script whose table names another build than versions.env.
    let root = w.path("root");
    std::fs::create_dir_all(format!("{root}/scripts/ci")).unwrap();
    std::fs::copy(repo("versions.env"), format!("{root}/versions.env")).unwrap();
    let copy = format!("{root}/scripts/ci/check-rule-alleles.sh");
    let src = read(&script);
    assert!(src.contains("\n#: # build "));
    std::fs::write(&copy, src.replace("\n#: # build ", "\n#: # build x")).unwrap();
    let o = Command::new("bash")
        .arg(&copy)
        .arg(repo("assets/genotype_rules.tsv"))
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&o.stderr).contains("versions.env pins"));
}

// 2. APOE

#[test]
fn the_apoe_row_is_the_golden_row_for_each_genotype_pair() {
    let w = Scratch::new("apoe");
    let rules = rules_path(&w);
    let mut bad = Vec::new();
    for (a, b, want) in APOE_GOLDEN {
        let recs: Vec<Vec<String>> = [("44908684", "T", "C", a), ("44908822", "C", "T", b)]
            .iter()
            .filter(|(.., gt)| *gt != ".")
            .map(|(pos, r, alt, gt)| {
                ["chr19", pos, r, alt, gt]
                    .iter()
                    .map(|s| (*s).to_owned())
                    .collect()
            })
            .collect();
        let (traits, _) = table(&w, &rules, &recs);
        let apoe = traits
            .iter()
            .find(|r| r["rsid"].contains('+'))
            .expect("the APOE row");
        let got = sha256(&OUT_COLUMNS.map(|k| apoe[k].as_str()).join("\t"));
        if got != want {
            bad.push(format!("{a}, {b}: {got}"));
        }
    }
    assert!(bad.is_empty(), "APOE goldens differ: {bad:?}");
}

// 3. multi-SNP rules

#[test]
fn each_multi_snp_rule_reads_its_combinations() {
    let w = Scratch::new("multi");
    let rules = rules_path(&w);
    let sites = sites();
    for (rid, cases) in MULTI_SNP {
        let s = sites
            .get(rid)
            .unwrap_or_else(|| panic!("{rid} is not in the rules table"));
        for (gts, want) in cases {
            let (_, lines) = table(&w, &rules, &records(s, gts));
            let got = rule(&lines, rid);
            assert_eq!(
                (got["status"].as_str(), got["result"].as_str()),
                ("called", *want),
                "{rid} {gts:?}: {got:?}"
            );
        }
        // One allele of the first site flipped is not the first result.
        let (gts, want) = cases[0];
        let (a, b) = gts[0].split_once('/').unwrap();
        let flip = if b == s[0][4] { &s[0][3] } else { &s[0][4] };
        let mut flipped: Vec<String> = gts.iter().map(|g| (*g).to_owned()).collect();
        flipped[0] = format!("{a}/{flip}");
        let flipped: Vec<&str> = flipped.iter().map(String::as_str).collect();
        let (_, lines) = table(&w, &rules, &records(s, &flipped));
        assert_ne!(
            rule(&lines, rid)["result"],
            want,
            "control: {rid} with one allele flipped"
        );
    }
}

/// Reference blocks only (`chrom pos ref . ref/ref 40` at every site) read
/// the label MULTI_SNP_REF holds for the rule.
fn reference_blocks_read(rid: &str) {
    let w = Scratch::new(&format!("ref-{rid}"));
    let rules = rules_path(&w);
    let want = MULTI_SNP_REF.iter().find(|(r, _)| *r == rid).unwrap().1;
    let recs: Vec<Vec<String>> = sites()[rid]
        .iter()
        .map(|s| {
            vec![
                s[1].clone(),
                s[2].clone(),
                s[3].clone(),
                ".".into(),
                format!("{0}/{0}", s[3]),
                "40".into(),
            ]
        })
        .collect();
    let (_, lines) = table(&w, &rules, &recs);
    let got = rule(&lines, rid);
    assert_eq!(
        (got["status"].as_str(), got["result"].as_str()),
        ("called", want),
        "{rid}, reference blocks only: {got:?}"
    );
}

#[test]
fn reference_blocks_apol1() {
    reference_blocks_read("APOL1");
}

#[test]
fn reference_blocks_serpina1() {
    reference_blocks_read("SERPINA1");
}

#[test]
fn reference_blocks_hfe() {
    reference_blocks_read("HFE");
}

#[test]
fn reference_blocks_gc() {
    reference_blocks_read("GC");
}

#[test]
fn reference_blocks_cad9p21() {
    reference_blocks_read("CAD9p21");
}

#[test]
fn reference_blocks_mapt() {
    reference_blocks_read("MAPT");
}

#[test]
fn every_multi_snp_rule_is_in_the_rules_table() {
    let s = sites();
    for (rid, _) in MULTI_SNP {
        assert!(s.contains_key(rid), "{rid}");
    }
}

#[test]
fn apol1_g0_g1_says_g2_is_not_read_and_the_other_copy_could_carry_it() {
    let w = Scratch::new("apol1");
    let (_, lines) = table(
        &w,
        &rules_path(&w),
        &records(&sites()["APOL1"], &["A/G", "T/G"]),
    );
    let reading = &rule(&lines, "APOL1")["reading"];
    assert!(
        reading.contains("G2 is not read") && reading.contains("G1/G2"),
        "{reading}"
    );
}

// 4. plain TSV

#[test]
fn quoting_or_an_embedded_tab_is_refused_and_the_outputs_before_stay() {
    let w = Scratch::new("plain");
    let g = w.put("g.tsv", "");
    std::fs::create_dir_all(w.path("traits")).unwrap();
    let (out, rout) = (
        w.path("traits/S_traits.tsv"),
        w.path("traits/S_genotype_rules.tsv"),
    );
    let source = read(&repo("assets/traits.tsv"));
    let rules = read(&repo("assets/genotype_rules.tsv"));
    let run = |t: &str, r: &str| {
        let args: Vec<String> = [
            "table",
            "--table",
            t,
            "--rules",
            r,
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
        og_steps::traits::cli(&args)
    };
    for (name, text, rules_text) in [
        (
            "quoted header",
            source.replacen("chrom\t", "\"chrom\"\t", 1),
            rules.clone(),
        ),
        (
            "quoted catalog field",
            source.replacen("\tMCM6\t", "\t\"MCM6\"\t", 1),
            rules.clone(),
        ),
        (
            "quoted rule text",
            source.clone(),
            rules.replacen("\tAPOE\tAPOE\t", "\tAPOE\t\"APOE\"\t", 1),
        ),
        (
            "embedded tab",
            source.replacen("\tMCM6\t", "\t\"MCM6\textra\"\t", 1),
            rules.clone(),
        ),
    ] {
        assert!(
            text != source || rules_text != rules,
            "{name}: nothing planted"
        );
        assert_eq!(
            run(
                &repo("assets/traits.tsv"),
                &repo("assets/genotype_rules.tsv")
            ),
            0
        );
        let before = (std::fs::read(&out).unwrap(), std::fs::read(&rout).unwrap());
        assert!(!before.0.is_empty() && !before.1.is_empty());
        let (bad, badrules) = (w.put("bad.tsv", &text), w.put("badrules.tsv", &rules_text));
        assert_eq!(run(&bad, &badrules), 1, "{name} is refused");
        assert_eq!(
            (std::fs::read(&out).unwrap(), std::fs::read(&rout).unwrap()),
            before,
            "{name}: both outputs stay as they were"
        );
    }
}

// 5. the LPA rows of traits.tsv (part 6 of tests/test_lpa_kiv2.py)

#[test]
fn traits_tsv_holds_the_two_lpa_risk_variants_and_the_blood_test_advice() {
    let text = read(&repo("assets/traits.tsv"));
    let body: String = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| l.to_owned() + "\n")
        .collect();
    let table = rows(&body);
    for (rs, alt) in [("rs10455872", "G"), ("rs3798220", "C")] {
        let r = table
            .iter()
            .find(|r| r["rsid"] == rs)
            .unwrap_or_else(|| panic!("no {rs} row"));
        assert_eq!(
            (r["gene"].as_str(), r["alt"].as_str()),
            ("LPA", alt),
            "{rs}"
        );
        let note = &r["note"];
        assert!(
            note.to_lowercase().contains("measure") && note.contains("blood"),
            "{rs}: {note}"
        );
    }
}

#[test]
fn the_test_reads_the_checkout() {
    assert!(Path::new(&repo("assets/genotype_rules.tsv")).is_file());
}
