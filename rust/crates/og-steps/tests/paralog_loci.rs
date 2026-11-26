//! assets/paralog_loci.tsv against the models of
//! assets/parascopy_models_v1.7.txt (the checks of the Python test the port
//! replaced, tests/test_paralog_loci.py; its readers' checks are og-couple's
//! `paralog_readers` test).
//!
//! Step 35 runs `parascopy cn-using` once per row of the catalog, with the
//! model file the row names. The models are the genes of the model list
//! (every region of each model of Parascopy's pinned data,
//! models_GRCh38_1KGP_v1.7.tar.gz, Zenodo record 15019940), each as
//! `<gene>.gz`, the file name list.txt gives it.
//!
//!  1. the model list: 167 genes, each with exactly one main region, each
//!     name `<gene>.gz`;
//!  2. the catalog: its columns in order, and every row valid: a model on
//!     the list, a GRCh38 region inside the model's main region,
//!     pscn_reliable yes or no, a known tag and cross-check, a table name
//!     used once, a source; the SMN1 row writes the table `smn`;
//!  3. controls: each planted fault is refused;
//!  4. one owner per locus: no catalog model overlaps a cn-depth row of
//!     bin/cn_loci.tsv but RHCE's cross-check of RHD; every catalog gene
//!     assets/callability_genes.tsv lists has the owner parascopy, and the
//!     catalog genes it does not list are the 17 of NOT_IN_CALLABILITY.

use std::collections::{BTreeSet, HashMap, HashSet};

fn repo(rel: &str) -> String {
    format!("{}/../../../{rel}", env!("CARGO_MANIFEST_DIR"))
}

const COLUMNS: [&str; 9] = [
    "model",
    "genes",
    "region",
    "reading",
    "pscn_reliable",
    "tag",
    "cross_check",
    "table",
    "source",
];
const TAGS: [&str; 5] = ["carrier", "SF", "disease", "PGx", "trait"];
const CROSS_CHECKS: [&str; 3] = ["paraphase", "cn-depth", "none"];
/// The one catalog model allowed to overlap a cn-depth row of
/// bin/cn_loci.tsv: RHCE's model covers most of RHD, and the catalog reads
/// it as a cross-check.
const CN_DEPTH_EXEMPT: [(&str, &str); 1] = [("RHCE.gz", "RHD")];
/// The catalog genes assets/callability_genes.tsv does not list.
/// build-callability.sh takes its genes from GIAB's CMRG list and genome
/// stratifications only, never from this catalog: seven are pseudogenes,
/// which it leaves out as not protein-coding, and none of its reasons
/// applies to the other ten. A new catalog gene missing there, or one of
/// these listed, changes this set and goes red: rebuild, or extend the set.
const NOT_IN_CALLABILITY: [&str; 17] = [
    "CELP", "GBAP1", "NCF1B", "NCF1C", "OTOAP1", "PMS2CL", "STRCP1", "ABCC6", "CFH", "CFHR4",
    "CHRNA7", "DEFA3", "KRT86", "NEB", "PIK3CA", "PLG", "VWF",
];

/// Data lines (no comment, no blank), split on tabs.
fn data_lines(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap()
}

/// One `<gene>.gz` per main line of the region file.
fn read_models(text: &str) -> Vec<String> {
    data_lines(text)
        .into_iter()
        .filter(|r| r.len() == 3 && r[2] == "main")
        .map(|r| format!("{}.gz", r[0]))
        .collect()
}

/// `^\d+$`, as the Python's `\d+`.
fn digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// `chr[0-9XY]+:\d+-\d+` -> (chrom, start, end), 1-based and inclusive.
fn span(text: &str) -> Option<(String, i64, i64)> {
    let (c, rest) = text.split_once(':')?;
    let (s, e) = rest.split_once('-')?;
    let ok = c.len() > 3
        && c.starts_with("chr")
        && c[3..]
            .bytes()
            .all(|b| b.is_ascii_digit() || b == b'X' || b == b'Y')
        && digits(s)
        && digits(e);
    if !ok {
        return None;
    }
    Some((c.to_owned(), s.parse().ok()?, e.parse().ok()?))
}

type Regions = HashMap<String, Vec<(String, i64, i64, String)>>;

/// `<gene>.gz` -> every region of the model, 1-based.
fn model_regions(text: &str) -> Regions {
    let mut out: Regions = HashMap::new();
    for r in data_lines(text) {
        let (c, rest) = r[1].split_once(':').unwrap();
        let (s, e) = rest.split_once('-').unwrap();
        out.entry(format!("{}.gz", r[0])).or_default().push((
            c.to_owned(),
            s.parse().unwrap(),
            e.parse().unwrap(),
            r[2].clone(),
        ));
    }
    out
}

fn table_rows(text: &str) -> Vec<HashMap<String, String>> {
    let lines = data_lines(text);
    lines[1..]
        .iter()
        .map(|r| lines[0].iter().cloned().zip(r.iter().cloned()).collect())
        .collect()
}

/// (locus, chrom, start, end) of every target and control of the cn_loci
/// rows the depth engine owns, 1-based inclusive (the file is BED).
fn cn_depth_parts(text: &str) -> Vec<(String, String, i64, i64)> {
    let mut out = Vec::new();
    for r in table_rows(text) {
        if r["owner"] != "cn-depth" {
            continue;
        }
        for part in format!("{};{}", r["targets"], r["controls"]).split(';') {
            let c: Vec<&str> = part.split(':').collect();
            if c.len() >= 2
                && let Some((s, e)) = c[1].split_once('-')
            {
                out.push((
                    r["locus"].clone(),
                    r["chrom"].clone(),
                    s.parse::<i64>().unwrap() + 1,
                    e.parse().unwrap(),
                ));
            }
        }
    }
    out
}

/// Every overlap of a catalog model's region with a cn-depth row's target
/// or control that the exemptions do not name.
fn owner_problems(
    models: &[String],
    regions: &Regions,
    parts: &[(String, String, i64, i64)],
    exempt: &[(&str, &str)],
) -> Vec<String> {
    let mut out = Vec::new();
    for m in models {
        for (c, s, e, kind) in regions.get(m).into_iter().flatten() {
            for (locus, cc, cs, ce) in parts {
                if c == cc && s <= ce && cs <= e && !exempt.contains(&(m.as_str(), locus.as_str()))
                {
                    out.push(format!(
                        "{m} ({kind} {c}:{s}-{e}) overlaps cn_loci row {locus} ({cc}:{cs}-{ce})"
                    ));
                }
            }
        }
    }
    out
}

/// Every catalog gene the callability catalog lists without the owner parascopy.
fn callability_problems(
    catalog: &HashSet<String>,
    rows: &[HashMap<String, String>],
) -> Vec<String> {
    rows.iter()
        .filter(|r| catalog.contains(&r["gene"]) && r["owner"] != "parascopy")
        .map(|r| format!("{} has owner {}", r["gene"], r["owner"]))
        .collect()
}

/// A word of the shapes scripts/35-paralogs.sh accepts.
fn model_word(s: &str) -> bool {
    s.len() > 3
        && s.ends_with(".gz")
        && s[..s.len() - 3]
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
}

fn region_word(s: &str) -> bool {
    let nz = |x: &str| digits(x) && !x.starts_with('0');
    span(s).is_some_and(|_| {
        let (_, rest) = s.split_once(':').unwrap();
        let (a, b) = rest.split_once('-').unwrap();
        nz(a) && nz(b)
    })
}

fn table_word(s: &str) -> bool {
    let mut b = s.bytes();
    b.next().is_some_and(|c| c.is_ascii_alphanumeric())
        && b.all(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'-')
}

/// Every problem of a catalog (its text, comments included).
fn problems(text: &str, models: &HashSet<String>, regions: Option<&Regions>) -> Vec<String> {
    let rows = data_lines(text);
    if rows.is_empty() || rows[0] != COLUMNS {
        return vec![format!("header is {:?}, want {COLUMNS:?}", rows.first())];
    }
    let mut out = Vec::new();
    let mut tables = HashSet::new();
    for (i, r) in rows[1..].iter().enumerate() {
        let n = i + 2;
        if r.len() != COLUMNS.len() {
            out.push(format!("row {n}: {} columns", r.len()));
            continue;
        }
        let c: HashMap<&str, &str> = COLUMNS
            .iter()
            .copied()
            .zip(r.iter().map(String::as_str))
            .collect();
        if !models.contains(c["model"]) {
            out.push(format!(
                "row {n}: model {} is not in the pinned model list",
                c["model"]
            ));
        }
        let contigs: Vec<String> = (1..=22)
            .map(|i| format!("chr{i}"))
            .chain(["chrX".into(), "chrY".into()])
            .collect();
        match span(c["region"]) {
            Some((ch, s, e)) if contigs.contains(&ch) && s >= 1 && s <= e => {
                if let Some(m) = regions.and_then(|g| g.get(c["model"]))
                    && !m
                        .iter()
                        .any(|(mc, ms, me, k)| k == "main" && *mc == ch && *ms <= s && e <= *me)
                {
                    out.push(format!(
                        "row {n}: region {} is not inside the main region of {}",
                        c["region"], c["model"]
                    ));
                }
            }
            _ => out.push(format!(
                "row {n}: region {} is not chrom:start-end on GRCh38",
                c["region"]
            )),
        }
        if !CROSS_CHECKS.contains(&c["cross_check"]) {
            out.push(format!("row {n}: cross_check {}", c["cross_check"]));
        }
        if !["yes", "no"].contains(&c["pscn_reliable"]) {
            out.push(format!("row {n}: pscn_reliable {}", c["pscn_reliable"]));
        }
        if !TAGS.contains(&c["tag"]) {
            out.push(format!("row {n}: tag {}", c["tag"]));
        }
        // model, region and table reach a shell loop and file names: the
        // patterns of scripts/35-paralogs.sh (a table is never "." or "..")
        for (k, ok) in [
            ("model", model_word(c["model"])),
            ("region", region_word(c["region"])),
            ("table", table_word(c["table"])),
        ] {
            if !ok {
                out.push(format!(
                    "row {n}: {k} {:?} is not what the step accepts",
                    c[k]
                ));
            }
        }
        if tables.contains(c["table"]) || ["hba", "paraphase"].contains(&c["table"]) {
            out.push(format!("row {n}: table {} is used twice", c["table"]));
        }
        tables.insert(c["table"].to_owned());
        for k in ["genes", "reading", "source"] {
            if c[k].trim().is_empty() {
                out.push(format!("row {n}: empty {k}"));
            }
        }
    }
    if !rows[1..].iter().any(|r| {
        r.first().map(String::as_str) == Some("SMN1.gz")
            && r.len() == COLUMNS.len()
            && r[7] == "smn"
    }) {
        out.push("no SMN1.gz row writing the table smn, which the reports read".into());
    }
    out
}

struct Fix {
    models: HashSet<String>,
    regions: Regions,
    catalog: String,
}

fn fix() -> Fix {
    let text = read(&repo("assets/parascopy_models_v1.7.txt"));
    Fix {
        models: read_models(&text).into_iter().collect(),
        regions: model_regions(&text),
        catalog: read(&repo("assets/paralog_loci.tsv")),
    }
}

#[test]
fn the_model_list() {
    let text = read(&repo("assets/parascopy_models_v1.7.txt"));
    let models = read_models(&text);
    let genes: HashSet<String> = data_lines(&text)
        .into_iter()
        .map(|r| r[0].clone())
        .collect();
    assert_eq!(models.len(), 167, "the model list has 167 names");
    let unique: HashSet<&String> = models.iter().collect();
    assert!(
        unique.len() == models.len() && models.len() == genes.len(),
        "each model has one main region"
    );
    assert!(
        models.iter().all(|m| model_word(m)),
        "each name is <model>.gz"
    );
    assert!(unique.contains(&"SMN1.gz".to_owned()) && unique.contains(&"GBA.gz".to_owned()));
}

#[test]
fn the_catalog_and_its_controls() {
    let f = fix();
    let found = problems(&f.catalog, &f.models, Some(&f.regions));
    assert!(found.is_empty(), "every row is valid: {}", found.join("; "));
    assert_eq!(table_rows(&f.catalog).len(), 37, "the SMN1 row and 36 more");
    let head = format!("{}\n", COLUMNS.join("\t"));
    let smn = "SMN1.gz\tSMN1,SMN2\tchr5:70951946-70951946\tx\tyes\tcarrier\tparaphase\tsmn\ts\n";
    let alone = problems(&format!("{head}{smn}"), &f.models, None);
    assert!(
        alone.is_empty(),
        "control: the SMN1 row alone is valid: {alone:?}"
    );
    for (desc, row) in [
        (
            "a model not on the list (GBA1.gz)",
            "GBA1.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba\ts\n",
        ),
        (
            "a model the data does not have (CYP21A2.gz)",
            "CYP21A2.gz\tCYP21A2\tchr6:32038265-32041670\tx\tno\tcarrier\tnone\tcyp21a2\ts\n",
        ),
        (
            "a region outside its model's main region",
            "GBA.gz\tGBA1\tchr5:70951946-70951946\tx\tno\tcarrier\tnone\tgba\ts\n",
        ),
        (
            "an unknown cross_check",
            "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tgauchian\tgba\ts\n",
        ),
        (
            "a region on no GRCh38 contig",
            "GBA.gz\tGBA1\tchr1_alt:1-2\tx\tno\tcarrier\tnone\tgba\ts\n",
        ),
        (
            "a region that ends before it starts",
            "GBA.gz\tGBA1\tchr1:20-10\tx\tno\tcarrier\tnone\tgba\ts\n",
        ),
        (
            "a repeated table name",
            "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tsmn\ts\n",
        ),
        (
            "a table name with a space",
            "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tg a\ts\n",
        ),
        (
            "a table name the step refuses (gba:2)",
            "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba:2\ts\n",
        ),
        (
            "a table name that is a path (..)",
            "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\t..\ts\n",
        ),
        (
            "an unknown tag",
            "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcancer\tnone\tgba\ts\n",
        ),
    ] {
        assert!(
            !problems(&format!("{head}{smn}{row}"), &f.models, Some(&f.regions)).is_empty(),
            "control: {desc} is refused"
        );
    }
    let no_smn = "GBA.gz\tGBA1\tchr1:155234452-155244699\tx\tno\tcarrier\tnone\tgba\ts\n";
    assert!(
        !problems(&format!("{head}{no_smn}"), &f.models, None).is_empty(),
        "a catalog without the SMN1 row"
    );
}

#[test]
fn one_owner_per_locus() {
    let f = fix();
    let rows = table_rows(&f.catalog);
    let parts = cn_depth_parts(&read(&repo("bin/cn_loci.tsv")));
    let loci: BTreeSet<&str> = parts.iter().map(|p| p.0.as_str()).collect();
    assert!(
        loci.contains("HBA") && loci.contains("RHD"),
        "cn-depth rows to compare with: {loci:?}"
    );
    let models: Vec<String> = rows.iter().map(|r| r["model"].clone()).collect();
    let found = owner_problems(&models, &f.regions, &parts, &CN_DEPTH_EXEMPT);
    assert!(
        found.is_empty(),
        "no overlap but the RHCE cross-check: {:?}",
        &found[..found.len().min(5)]
    );
    assert!(
        owner_problems(&models, &f.regions, &parts, &[])
            .iter()
            .any(|x| x.contains("RHCE.gz") && x.contains("RHD")),
        "control: without the exemption RHCE's model is refused"
    );
    let mut planted = f.regions.clone();
    planted.insert(
        "HBAX.gz".into(),
        vec![("chr16".into(), 172000, 178000, "main".into())],
    );
    assert!(
        !owner_problems(&["HBAX.gz".into()], &planted, &parts, &CN_DEPTH_EXEMPT).is_empty(),
        "control: a model over HBA1 and HBA2 is refused"
    );
    let rhce = rows
        .iter()
        .find(|r| r["model"] == "RHCE.gz")
        .expect("an RHCE row");
    assert_eq!(rhce["cross_check"], "cn-depth");
    assert!(!rhce["genes"].split(',').any(|g| g == "RHD"));
    let catalog: HashSet<String> = rows
        .iter()
        .flat_map(|r| r["genes"].split(',').map(str::to_owned))
        .collect();
    let cal = table_rows(&read(&repo("assets/callability_genes.tsv")));
    let found = callability_problems(&catalog, &cal);
    assert!(found.is_empty(), "{}", found.join("; "));
    let listed: HashSet<&str> = cal
        .iter()
        .map(|r| r["gene"].as_str())
        .filter(|g| catalog.contains(*g))
        .collect();
    assert!(
        ["GBA1", "C4A", "AMY1C", "SMN1"]
            .iter()
            .all(|g| listed.contains(g))
    );
    let unlisted = |cal: &[HashMap<String, String>]| -> BTreeSet<String> {
        let have: HashSet<&str> = cal.iter().map(|r| r["gene"].as_str()).collect();
        catalog
            .iter()
            .filter(|g| !have.contains(g.as_str()))
            .cloned()
            .collect()
    };
    let want: BTreeSet<String> = NOT_IN_CALLABILITY.iter().map(|s| (*s).to_owned()).collect();
    assert_eq!(
        unlisted(&cal),
        want,
        "the 17 known catalog genes the callability catalog does not list"
    );
    let mut planted_neb = cal.clone();
    let mut neb = cal[0].clone();
    neb.insert("gene".into(), "NEB".into());
    neb.insert("owner".into(), "parascopy".into());
    planted_neb.push(neb);
    assert_ne!(
        unlisted(&planted_neb),
        want,
        "control: a row for NEB makes the set differ"
    );
    let planted_cal: Vec<HashMap<String, String>> = cal
        .iter()
        .map(|r| {
            let mut r = r.clone();
            if r["gene"] == "C4A" {
                r.insert("owner".into(), ".".into());
            }
            r
        })
        .collect();
    assert!(
        !callability_problems(&catalog, &planted_cal).is_empty(),
        "control: owner '.' is refused"
    );
}
