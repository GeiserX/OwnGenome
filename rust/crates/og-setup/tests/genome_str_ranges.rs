//! The repeat ranges: the `str_ranges.tsv` og str catalog writes, the subset
//! the report cites (`bin/str_intermediate_ranges.tsv`) and the docs. These
//! are the checks of the old `tests/test_str_ranges.py`, by meaning.
//!
//! 1. `str_ranges.tsv`, beside the catalog built from the pinned STRchive
//!    files: one row per ExpansionHunter locus (82) traced to its STRchive
//!    id, with the catalog's NormalMax and PathologicMin, empty where the
//!    catalog has none and the reason; the intermediate range is STRchive's.
//!    Controls: one row dropped, or one row twice, is caught.
//! 2. `bin/str_intermediate_ranges.tsv`: the columns og summary collect
//!    reads; HTT, FMR1 and DMPK; one row per locus; whole numbers with min <=
//!    max ending right below the next class; a label and a GeneReviews
//!    citation; `docs/results/structural-and-repeats.md` prints the same range
//!    in the locus's row; and the range lies above the catalog's normal
//!    maximum and below its pathogenic minimum (control: HTT 20-35 is
//!    caught). That og reads every row is og-summary's own test.
//! 3. The "Key Disease Thresholds" table of `docs/09-str-expansions.md`
//!    gives the catalog's thresholds (control: HTT's minimum changed to 37 is
//!    caught).
//!
//! Part 2's own checks run everywhere; what compares with the catalog
//! downloads STRchive, so it runs where `CI` is set or with
//! `OG_NET_TESTS=1`, as `tests/strchive_catalog.rs` does.

use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use pycompat::json::{self, Value};

const COLUMNS: [&str; 7] = [
    "locus",
    "intermediate_min",
    "intermediate_max",
    "next_min",
    "next_class",
    "label",
    "citation",
];
const RANGE_COLUMNS: [&str; 7] = [
    "locus",
    "source_id",
    "normal_max",
    "intermediate_min",
    "intermediate_max",
    "pathogenic_min",
    "no_thresholds",
];

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn read(rel: &str) -> String {
    std::fs::read_to_string(repo().join(rel)).unwrap()
}

fn online(test: &str) -> bool {
    let on = std::env::var_os("CI").is_some_and(|v| !v.is_empty())
        || std::env::var("OG_NET_TESTS").as_deref() == Ok("1");
    if !on {
        eprintln!("skip: {test} downloads STRchive; set OG_NET_TESTS=1 (CI always runs it)");
    }
    on
}

fn env(name: &str) -> String {
    let text = read("versions.env");
    let prefix = format!("{name}=");
    let line = text
        .lines()
        .find(|l| l.starts_with(&prefix))
        .unwrap_or_else(|| panic!("no {name} in versions.env"));
    line[prefix.len()..]
        .split('#')
        .next()
        .unwrap()
        .trim()
        .trim_matches('"')
        .to_owned()
}

fn sha256(path: &Path) -> String {
    for (cmd, args) in [("sha256sum", &[][..]), ("shasum", &["-a", "256"][..])] {
        if let Ok(o) = Command::new(cmd).args(args).arg(path).output()
            && o.status.success()
        {
            return String::from_utf8(o.stdout)
                .unwrap()
                .split_whitespace()
                .next()
                .unwrap_or_default()
                .to_owned();
        }
    }
    panic!("neither sha256sum nor shasum could hash {}", path.display());
}

/// The table rows of a TSV: no `#` line, no blank line; (header, rows).
fn tsv(text: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let head = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .map(str::to_owned)
        .collect();
    (
        head,
        lines
            .map(|l| l.split('\t').map(str::to_owned).collect())
            .collect(),
    )
}

/// What the converter writes on the pinned files: (catalog entries, the
/// ranges table, the STRchive loci records by id).
struct Built {
    catalog: Vec<Value>,
    ranges: String,
    loci: HashMap<String, Value>,
}

fn built() -> &'static Built {
    static B: OnceLock<Built> = OnceLock::new();
    B.get_or_init(|| {
        let rel = env("STRCHIVE_VERSION");
        let dir = std::env::var_os("STRCHIVE_DIR").map_or_else(
            || std::env::temp_dir().join(format!("pgp-strchive-{rel}")),
            PathBuf::from,
        );
        std::fs::create_dir_all(&dir).unwrap();
        let mut paths = Vec::new();
        for (name, var) in [
            (
                format!("STRchive-disease-loci-{rel}.hg38.expansionhunter.json"),
                "STRCHIVE_EH_SHA256",
            ),
            (format!("STRchive-loci-{rel}.json"), "STRCHIVE_LOCI_SHA256"),
        ] {
            let path = dir.join(&name);
            if !(path.exists() && sha256(&path) == env(var)) {
                let part = dir.join(format!("{name}.part"));
                let url = format!(
                    "https://github.com/dashnowlab/STRchive/releases/download/{rel}/{name}"
                );
                let o = Command::new("curl")
                    .args([
                        "-sSLf",
                        "--retry",
                        "3",
                        "--retry-delay",
                        "5",
                        "--max-time",
                        "120",
                        "-o",
                    ])
                    .arg(&part)
                    .arg(&url)
                    .output()
                    .unwrap();
                assert!(
                    o.status.success(),
                    "curl {url}: {}",
                    String::from_utf8_lossy(&o.stderr)
                );
                std::fs::rename(&part, &path).unwrap();
            }
            assert_eq!(
                sha256(&path),
                env(var),
                "{name} does not have the sha256 {var} pins"
            );
            paths.push(path.to_string_lossy().into_owned());
        }
        let out = std::env::temp_dir().join(format!("og-genome-str-ranges-{}", std::process::id()));
        std::fs::create_dir_all(&out).unwrap();
        let (cat, rng) = (out.join("catalog.json"), out.join("str_ranges.tsv"));
        let args: Vec<String> = [
            "--eh",
            &paths[0],
            "--loci",
            &paths[1],
            "--release",
            &rel,
            "--out",
            &cat.to_string_lossy(),
            "--ranges",
            &rng.to_string_lossy(),
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
        assert_eq!(
            og_setup::strchive::cli(&args),
            0,
            "the converter on the pinned files"
        );
        let catalog = match json::loads(&std::fs::read_to_string(&cat).unwrap()).unwrap() {
            Value::Array(a) => a,
            v => panic!("{v:?}"),
        };
        let ranges = std::fs::read_to_string(&rng).unwrap();
        let loci = match json::loads(&std::fs::read_to_string(&paths[1]).unwrap()).unwrap() {
            Value::Array(a) => a
                .into_iter()
                .map(|x| (s(f(&x, "id")).to_owned(), x))
                .collect(),
            v => panic!("{v:?}"),
        };
        let _ = std::fs::remove_dir_all(&out);
        Built {
            catalog,
            ranges,
            loci,
        }
    })
}

fn f<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(k),
        _ => None,
    }
}

fn s(v: Option<&Value>) -> &str {
    match v {
        Some(Value::Str(s)) => s,
        _ => "",
    }
}

fn int(v: Option<&Value>) -> Option<i128> {
    match v {
        Some(Value::Int(i)) => Some(*i),
        _ => None,
    }
}

fn num(x: &str) -> Option<i128> {
    if x.is_empty() {
        None
    } else {
        Some(x.parse().unwrap())
    }
}

/// What is wrong with ranges rows against the catalog (empty: nothing).
fn range_problems(rows: &[Vec<String>], cat: &[Value]) -> Vec<String> {
    let by: HashMap<&str, &Value> = cat.iter().map(|e| (s(f(e, "LocusId")), e)).collect();
    let mut bad = Vec::new();
    let loci: Vec<&str> = rows.iter().map(|r| r[0].as_str()).collect();
    let (have, want): (BTreeSet<&str>, BTreeSet<&str>) =
        (loci.iter().copied().collect(), by.keys().copied().collect());
    if have != want {
        bad.push(format!(
            "loci differ: missing {:?}, extra {:?}",
            &want - &have,
            &have - &want
        ));
    }
    if loci.len() != have.len() {
        bad.push("a locus twice".into());
    }
    for r in rows {
        let Some(e) = by
            .get(r[0].as_str())
            .filter(|_| r.len() == RANGE_COLUMNS.len())
        else {
            bad.push(format!("bad row {r:?}"));
            continue;
        };
        if r[1] != s(f(e, "SourceId")) {
            bad.push(format!("{}: source_id", r[0]));
        }
        if (num(&r[2]), num(&r[5])) != (int(f(e, "NormalMax")), int(f(e, "PathologicMin"))) {
            bad.push(format!(
                "{}: {}/{} are not the catalog's thresholds",
                r[0], r[2], r[5]
            ));
        }
        if r[6] != s(f(e, "NoThresholds")) {
            bad.push(format!("{}: no_thresholds", r[0]));
        }
    }
    bad
}

fn generated() -> HashMap<String, HashMap<&'static str, String>> {
    let (_, rows) = tsv(&built().ranges);
    rows.into_iter()
        .map(|r| (r[0].clone(), RANGE_COLUMNS.iter().copied().zip(r).collect()))
        .collect()
}

#[test]
fn str_ranges_tsv_has_one_row_per_catalog_locus_with_its_thresholds() {
    if !online("str_ranges_tsv_has_one_row_per_catalog_locus_with_its_thresholds") {
        return;
    }
    let b = built();
    let (head, rows) = tsv(&b.ranges);
    assert_eq!(head, RANGE_COLUMNS);
    assert_eq!(rows.len(), 82);
    assert_eq!(range_problems(&rows, &b.catalog), Vec::<String>::new());
    assert!(
        !range_problems(&rows[1..], &b.catalog).is_empty(),
        "control: one row dropped"
    );
    let twice: Vec<Vec<String>> = rows.iter().chain(&rows[..1]).cloned().collect();
    assert!(
        !range_problems(&twice, &b.catalog).is_empty(),
        "control: one row twice"
    );
    let ranges = generated();
    for (locus, d) in &ranges {
        let rec = &b.loci[&d["source_id"]];
        assert_eq!(
            (num(&d["intermediate_min"]), num(&d["intermediate_max"])),
            (
                int(f(rec, "intermediate_min")),
                int(f(rec, "intermediate_max"))
            ),
            "{locus}: STRchive's intermediate range"
        );
    }
    let atxn2 = &ranges["ATXN2"];
    assert_eq!(
        RANGE_COLUMNS[2..6]
            .iter()
            .map(|k| atxn2[k].as_str())
            .collect::<Vec<_>>(),
        ["28", "29", "34", "35"],
        "ATXN2: normal to 28, intermediate 29-34, pathogenic from 35"
    );
    let muc1 = &ranges["MUC1"];
    assert!(
        muc1["normal_max"].is_empty()
            && muc1["no_thresholds"].contains("short reads cannot size it"),
        "{muc1:?}"
    );
}

/// `GeneReviews: <title>, https://www.ncbi.nlm.nih.gov/books/NBK<digits>/`.
fn cited(c: &str) -> bool {
    let Some(rest) = c.strip_prefix("GeneReviews: ") else {
        return false;
    };
    let Some((title, url)) = rest.rsplit_once(", ") else {
        return false;
    };
    let id = url
        .strip_prefix("https://www.ncbi.nlm.nih.gov/books/NBK")
        .and_then(|x| x.strip_suffix('/'));
    !title.is_empty() && id.is_some_and(|d| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()))
}

/// Locus -> its row of the repeat table in docs/results/structural-and-repeats.md.
fn doc_rows() -> HashMap<String, String> {
    let mut out = HashMap::new();
    for line in read("docs/results/structural-and-repeats.md").lines() {
        let c: Vec<&str> = line
            .trim()
            .trim_matches('|')
            .split('|')
            .map(str::trim)
            .collect();
        if line.starts_with('|') && c.len() >= 6 {
            out.entry(c[0].to_uppercase())
                .or_insert_with(|| line.to_owned());
        }
    }
    out
}

/// `lo-hi` in `row` with no digit or dash right before or after it.
fn prints_range(row: &str, lo: i64, hi: i64) -> bool {
    let needle = format!("{lo}-{hi}");
    row.match_indices(&needle).any(|(i, _)| {
        let before = row[..i].chars().last();
        let after = row[i + needle.len()..].chars().next();
        let free = |c: Option<char>| !c.is_some_and(|c| c.is_ascii_digit() || c == '-');
        free(before) && free(after)
    })
}

/// The subset's rows, checked on their own.
fn subset() -> Vec<HashMap<&'static str, String>> {
    let (head, body) = tsv(&read("bin/str_intermediate_ranges.tsv"));
    assert_eq!(head, COLUMNS, "the columns og summary collect reads");
    body.into_iter()
        .map(|r| {
            assert_eq!(r.len(), COLUMNS.len(), "{r:?}");
            COLUMNS.iter().copied().zip(r).collect()
        })
        .collect()
}

#[test]
fn the_cited_subset_is_well_formed_and_the_docs_print_its_ranges() {
    let rows = subset();
    let loci: Vec<&str> = rows.iter().map(|r| r["locus"].as_str()).collect();
    for l in ["HTT", "FMR1", "DMPK"] {
        assert!(loci.contains(&l), "{l}");
    }
    assert_eq!(
        loci.iter().collect::<BTreeSet<_>>().len(),
        loci.len(),
        "one row per locus"
    );
    let docs = doc_rows();
    for d in &rows {
        let l = &d["locus"];
        let n: Vec<i64> = ["intermediate_min", "intermediate_max", "next_min"]
            .iter()
            .map(|k| {
                d[k].parse()
                    .unwrap_or_else(|_| panic!("{l}: {k} is not a whole number"))
            })
            .collect();
        let (lo, hi, next) = (n[0], n[1], n[2]);
        assert!(lo <= hi, "{l}: {lo}-{hi}");
        assert_eq!(
            hi + 1,
            next,
            "{l}: the range ends right below {}",
            d["next_class"]
        );
        assert!(
            !d["label"].trim().is_empty() && !d["next_class"].trim().is_empty(),
            "{l}"
        );
        assert!(cited(&d["citation"]), "{l}: {}", d["citation"]);
        let row = docs.get(&l.to_uppercase()).map_or("", String::as_str);
        assert!(
            prints_range(row, lo, hi),
            "{l}: the results doc prints {lo}-{hi} in its row: {row}"
        );
    }
    assert!(
        !cited("GeneReviews: x, https://example.org/NBK1/"),
        "control: a citation elsewhere"
    );
    assert!(
        !prints_range("| HTT | 127-35 |", 27, 35),
        "control: a longer number"
    );
}

#[test]
fn each_cited_range_lies_between_the_catalogs_normal_and_pathogenic() {
    if !online("each_cited_range_lies_between_the_catalogs_normal_and_pathogenic") {
        return;
    }
    let ranges = generated();
    let inside = |locus: &str, lo: i64, hi: i64| {
        ranges.get(locus).is_some_and(|g| {
            !g["normal_max"].is_empty()
                && g["normal_max"].parse::<i64>().unwrap() < lo
                && hi < g["pathogenic_min"].parse::<i64>().unwrap()
        })
    };
    for d in subset() {
        let (lo, hi) = (
            d["intermediate_min"].parse().unwrap(),
            d["intermediate_max"].parse().unwrap(),
        );
        assert!(
            inside(&d["locus"], lo, hi),
            "{}: {lo}-{hi} against {:?}",
            d["locus"],
            ranges.get(&d["locus"])
        );
    }
    assert!(
        !inside("HTT", 20, 35),
        "control: a planted HTT range of 20-35"
    );
}

/// Gene -> (normal max, pre-mutation band, full-mutation minimum) of the
/// Key Disease Thresholds table.
fn threshold_rows(text: &str) -> HashMap<String, (String, String, String)> {
    let (mut out, mut on) = (HashMap::new(), false);
    for line in text.lines() {
        if line.starts_with("## ") {
            on = line.trim() == "## Key Disease Thresholds";
            continue;
        }
        let c: Vec<&str> = line
            .trim()
            .trim_matches('|')
            .split('|')
            .map(str::trim)
            .collect();
        if on && line.starts_with('|') && c.len() == 6 && c[1] != "Gene" && c[1] != "---" {
            out.insert(
                c[1].to_owned(),
                (c[3].to_owned(), c[4].to_owned(), c[5].to_owned()),
            );
        }
    }
    out
}

fn threshold_problems(
    rows: &HashMap<String, (String, String, String)>,
    cat: &[Value],
) -> Vec<String> {
    let by: HashMap<&str, &Value> = cat.iter().map(|e| (s(f(e, "LocusId")), e)).collect();
    let mut bad = if rows.is_empty() {
        vec!["no rows".to_owned()]
    } else {
        Vec::new()
    };
    for (gene, got) in rows {
        let (Some(n), Some(p)) = (
            by.get(gene.as_str()).and_then(|e| int(f(e, "NormalMax"))),
            by.get(gene.as_str())
                .and_then(|e| int(f(e, "PathologicMin"))),
        ) else {
            bad.push(format!("{gene}: no thresholds in the catalog"));
            continue;
        };
        let want = (
            format!("≤{n}"),
            format!("{}-{}", n + 1, p - 1),
            format!("≥{p}"),
        );
        if *got != want {
            bad.push(format!(
                "{gene}: the docs say {got:?}, the catalog {want:?}"
            ));
        }
    }
    bad
}

#[test]
fn the_key_disease_thresholds_table_lists_the_report_loci() {
    let rows = threshold_rows(&read("docs/09-str-expansions.md"));
    for g in ["HTT", "FMR1", "FXN", "C9ORF72", "DMPK", "ATXN1", "ATXN2"] {
        assert!(rows.contains_key(g), "{g}: {rows:?}");
    }
}

#[test]
fn the_key_disease_thresholds_table_gives_the_catalogs_thresholds() {
    if !online("the_key_disease_thresholds_table_gives_the_catalogs_thresholds") {
        return;
    }
    let text = read("docs/09-str-expansions.md");
    let rows = threshold_rows(&text);
    assert_eq!(
        threshold_problems(&rows, &built().catalog),
        Vec::<String>::new()
    );
    let from = "| HTT | CAG | ≤26 | 27-35 | ≥36 |";
    assert!(text.contains(from));
    let planted = threshold_rows(&text.replace(from, "| HTT | CAG | ≤26 | 27-35 | ≥37 |"));
    assert!(
        planted != rows && !threshold_problems(&planted, &built().catalog).is_empty(),
        "control: HTT at 37"
    );
}
