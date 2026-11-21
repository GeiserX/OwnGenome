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
//!    Gilbert "carrier", CYP2A6 two *2 copies or nicotine clearance from one
//!    site, and the unsupported vitamin D comparisons; what replaced them
//!    names its limits (the *2/*4 ambiguity, the rare e1, the *80 tag, the
//!    population limits of FUT2 and UGT1A1, the MC1R cohort);
//! 4. the GC rule gives no protein level and no ranking: the monoclonal
//!    assay of Powe 2013 reads the types differently (Denburg 2016);
//! 5. the APOE rule says its odds ratios are of clinic and autopsy studies,
//!    with their age and sex limits, and the HFE rule cites the cohorts its
//!    penetrance comes from;
//! 6. `og traits table` on synthetic genotypes: the -22018 tag reads its
//!    "only when -13910 has no genotype" hint whether or not -13910 is
//!    covered, and the e1/e4 and e1/e2 pairs read e1 on both site rows;
//! 7. the docs of the genotype-rules table say its `result` never shows G2
//!    and point a script at `summary.json`, and that path holds a G2 result
//!    the table's result lacks (the summary-collect parity golden, og's own
//!    output).

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
    let gone: [(&str, &str, &str); 14] = [
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
        (
            "rs1801272",
            "hom_alt",
            "both copies make an inactive enzyme",
        ),
        ("rs2282679", "note", "matter more"),
        ("rs2282679", "note", "about as large"),
        ("rs182549", "note", "read on its own"),
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

    // CYP2A6: T/T can be *2/*2 or *2 with a deleted copy (*4).
    let cyp = row(&rows, "rs1801272");
    for want in ["*4", "deleted", "cannot tell", "step 32"] {
        assert!(cyp["hom_alt"].contains(want), "{want}: {}", cyp["hom_alt"]);
    }
    // APOE: one site alone never names the allele; the rare e1 is said.
    for rsid in ["rs429358", "rs7412"] {
        let r = row(&rows, rsid);
        for k in ["het", "hom_alt", "note"] {
            assert!(r[k].contains("e1"), "{rsid} {k}: {}", r[k]);
        }
    }
    // The -22018 tag: when to use it.
    assert!(
        row(&rows, "rs182549")["note"].contains("only when row rs4988235 has no genotype"),
        "{:?}",
        row(&rows, "rs182549")
    );
    // MC1R: Flanagan 2000 cited on every row, its findings as cohort ones.
    for rsid in ["rs1805007", "rs1805008", "rs1805009"] {
        let r = row(&rows, rsid);
        assert!(r["pmids"].split(',').any(|p| p == "11030758"), "{r:?}");
        for want in ["Flanagan 2000", "only one variant copy", "not rules"] {
            assert!(r["note"].contains(want), "{rsid} {want}: {}", r["note"]);
        }
        assert!(r["het"].contains("not a prediction for you"), "{r:?}");
    }
    // FUT2: population limits, the opt-in blood-group step.
    let fut2 = row(&rows, "rs601338");
    assert_eq!(fut2["ancestry"], "yes");
    for want in ["differs between populations", "weak secretor", "Step 41"] {
        assert!(fut2["note"].contains(want), "{want}: {}", fut2["note"]);
    }
    // Vitamin D: the Framingham means against supplements and season, and
    // the explained variance of the discovery cohorts.
    let vd = &row(&rows, "rs2282679")["note"];
    for want in ["83.4", "74.7", "winter", "comparable", "discovery cohorts"] {
        assert!(vd.contains(want), "{want}: {vd}");
    }
    // UGT1A1: the *80 tag of *28, European ancestry.
    let ugt = row(&rows, "rs887829");
    assert_eq!(ugt["ancestry"], "yes");
    for want in ["*80", "*28", "European ancestry"] {
        assert!(ugt["note"].contains(want), "{want}: {}", ugt["note"]);
    }
    assert!(ugt["het"].contains("bilirubin"), "{ugt:?}");
}

/// The lines of a rule in the rules table, split at tabs.
fn rule_lines(text: &str, rule: &str) -> Vec<Vec<String>> {
    text.lines()
        .map(|l| l.split('\t').map(str::to_owned).collect::<Vec<_>>())
        .filter(|c| c.len() > 2 && c[1] == rule)
        .collect()
}

#[test]
fn gc_readings_give_no_protein_level_and_no_ranking() {
    let lines = rule_lines(&read("assets/genotype_rules.tsv"), "GC");
    let mut n = 0;
    for c in lines.iter().filter(|c| c[0] == "combination") {
        n += 1;
        let reading = &c[4];
        for bad in ["micrograms", "mean", "lowest", "highest", "intermediate"] {
            assert!(!reading.contains(bad), "{bad}: {reading}");
        }
        assert!(reading.contains("not a level"), "{reading}");
    }
    assert_eq!(n, 6, "the six GC pairs");
    let rule = lines.iter().find(|c| c[0] == "rule").unwrap();
    let (note, pmids) = (&rule[6], &rule[7]);
    assert!(!note.contains("79.4"), "{note}");
    assert!(note.contains("Denburg 2016"), "{note}");
    assert!(note.contains("no protein level"), "{note}");
    assert!(pmids.split(',').any(|p| p == "27250744"), "{pmids}");
}

#[test]
fn apoe_and_hfe_rules_name_the_limits_of_their_numbers() {
    let text = read("assets/genotype_rules.tsv");
    let apoe = rule_lines(&text, "APOE");
    let rule = apoe.iter().find(|c| c[0] == "rule").unwrap();
    assert!(rule[3].contains("e1"), "{}", rule[3]);
    for want in [
        "clinic and autopsy",
        "not a lifetime risk",
        "after 70",
        "men and women",
        "rare e1",
    ] {
        assert!(rule[6].contains(want), "{want}: {}", rule[6]);
    }
    let hfe = rule_lines(&text, "HFE");
    let rule = hfe.iter().find(|c| c[0] == "rule").unwrap();
    for want in ["Allen 2008", "Gurrin 2009", "Beutler 2002", "28% of male"] {
        assert!(rule[6].contains(want), "{want}: {}", rule[6]);
    }
    for p in ["8696333", "18199861", "19554541", "11812557"] {
        assert!(rule[7].split(',').any(|x| x == p), "{p}: {}", rule[7]);
    }
}

/// `og traits table` over the shipped tables and these genotypes; the
/// traits output as one map per line, keyed by rsid.
fn table_of(name: &str, genotypes: &str) -> HashMap<String, Vec<String>> {
    let w = Scratch::new(name);
    std::fs::create_dir_all(w.path("traits")).unwrap();
    let g = w.path("g.tsv");
    std::fs::write(&g, genotypes).unwrap();
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
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<String> = l.split('\t').map(str::to_owned).collect();
            (c[0].clone(), c)
        })
        .collect()
}

#[test]
fn the_lactase_tag_reads_its_hint_with_or_without_the_primary_marker() {
    let rows = traits();
    let tag = row(&rows, "rs182549");
    // Only the tag covered: -13910 is not read, the tag gives its hint.
    let only_tag = table_of("tag-only", "chr2\t135859184\tC\tT\tC/T\t40\n");
    assert_eq!(only_tag["rs4988235"][6], "not_covered");
    let t = &only_tag["rs182549"];
    assert_eq!((t[5].as_str(), t[6].as_str()), ("C/T", "called"));
    assert_eq!(t[7], tag["het"]);
    assert!(t[7].contains("only when that row has no genotype"), "{t:?}");
    // Both covered: -13910 reads persistence itself; the tag still says it
    // is not a persistence allele.
    let both = table_of(
        "tag-both",
        "chr2\t135851076\tG\tA\tG/A\t40\nchr2\t135859184\tC\tT\tC/T\t40\n",
    );
    let lct = &both["rs4988235"];
    assert_eq!(
        (lct[6].as_str(), lct[7].as_str()),
        ("called", row(&rows, "rs4988235")["het"].as_str())
    );
    assert!(both["rs182549"][7].contains("not a persistence allele"));
}

#[test]
fn e1_pairs_read_e1_on_both_site_rows() {
    let rows = traits();
    // (rs429358, rs7412, the column of each site row's reading)
    for (a, b, ka, kb) in [
        ("C/C", "C/T", "hom_alt", "het"),
        ("T/C", "T/T", "het", "hom_alt"),
    ] {
        let g = format!("chr19\t44908684\tT\tC\t{a}\t40\nchr19\t44908822\tC\tT\t{b}\t40\n");
        let out = table_of(&format!("e1-{}", ka), &g);
        let (ra, rb) = (&out["rs429358"], &out["rs7412"]);
        assert_eq!(ra[7], row(&rows, "rs429358")[ka], "{a}/{b}");
        assert_eq!(rb[7], row(&rows, "rs7412")[kb], "{a}/{b}");
        assert!(
            ra[7].contains("e1") && rb[7].contains("e1"),
            "{ra:?} {rb:?}"
        );
        let apoe = out
            .iter()
            .find(|(k, _)| k.contains('+'))
            .map(|(_, v)| v)
            .expect("the APOE line");
        assert!(apoe[4].contains("e1"), "{apoe:?}");
        assert!(apoe[7].contains("rare e1 allele: not read"), "{apoe:?}");
    }
}

#[test]
fn the_genotype_rules_docs_point_scripts_at_the_g2_result() {
    for doc in ["docs/output-contract.md", "docs/42-traits.md"] {
        let text = read(doc);
        let line = text
            .lines()
            .find(|l| l.contains("sections.traits.data.rules"))
            .unwrap_or_else(|| panic!("{doc} does not name sections.traits.data.rules"));
        for want in [
            "_genotype_rules.tsv",
            "never shows G2",
            "`summary.json`",
            "sections.traits.data.rules",
            "APOL1",
            "`g2`",
            "`result`",
            "`reading`",
        ] {
            assert!(line.contains(want), "{doc}: {want} not in: {line}");
        }
    }
    // The documented path in og's own output: the summary-collect parity
    // golden of a sample with a G1 copy and a G2 copy. The APOL1 line's
    // `result` (the table's) has no G2; its `g2` has the G1 plus G2 answer.
    use pycompat::json::Value;
    let golden = read("tests/parity/golden/summary-collect/c-full/full.json");
    let v = pycompat::json::loads(&golden).unwrap();
    let get = |v: &Value, k: &str| -> Value {
        match v {
            Value::Object(m) => m.get(k).cloned().unwrap_or(Value::Null),
            _ => Value::Null,
        }
    };
    let text = |v: Value| match v {
        Value::Str(s) => s,
        other => panic!("not a string: {other:?}"),
    };
    let rules = get(&get(&get(&get(&v, "sections"), "traits"), "data"), "rules");
    let Value::Array(rules) = rules else {
        panic!("sections.traits.data.rules is not a list")
    };
    let apol1 = rules
        .iter()
        .find(|r| get(r, "rule") == Value::Str("APOL1".into()))
        .expect("an APOL1 line");
    let g2 = get(apol1, "g2");
    let (table_result, result, reading) = (
        text(get(apol1, "result")),
        text(get(&g2, "result")),
        text(get(&g2, "reading")),
    );
    assert!(!table_result.contains("G2"), "{table_result}");
    assert!(result.contains("G2"), "{result}");
    assert!(reading.contains(&result), "{reading}");
}
