//! `og cn` and `bin/cn_loci.tsv`: one depth-ratio engine for every
//! copy-number locus. These are the library and command-line checks of the
//! old `tests/test_copy_number.py`, by meaning; its checks of step 35's
//! shell blocks, of the Nextflow wiring, of `og dev build cn-intervals`, of
//! the summary and of the couple screen are `rust/og/tests/genome_steps.rs`.
//!
//!  1. the catalog: the gene, dosage and CYP21A2 rows first, then the
//!     interval rows; every row is well formed; a broken row is refused, and
//!     `og cn bed --catalog` refuses it with the same message;
//!  2. one owner per locus: no row's targets or controls overlap a region of
//!     a model the pinned Parascopy data holds, except the pairs EXEMPT
//!     lists; planted AMY1C, C4B, GYPA and DEFB4 rows are caught;
//!  3. parity: made-up HBA and LPA depths give, through `og cn call`, the
//!     tables the code before the engine wrote (`tests/data/genome_gold.json`,
//!     frozen from it);
//!  4. synthetic depths at 0 to 4 copies of a made-up gene are called right;
//!     the MAPQ 0 ceiling; a kind=interval row with pooled depths;
//!  5. a control with no read beside covered controls is BROKEN_CONTROL
//!     (exit 3 through `call`); no read anywhere is not assessed;
//!  6. the command line: `bed`, `call --table`, the bedcov reader's refusals;
//!  8. the gene-dosage rows (HP, LCE3BC, CCL3L1): 0 to 4 copies; CCL3L1's
//!     readings name no disease; a row shifted by 1 kb refuses its fixture;
//!  9. the interval rows: their content, their calls by sex, and
//!     `og cn intervals` over one bedcov table;
//! 10. the CYP21A2 row (call psv): its PSV sites, the refused variants, the
//!     committed fixtures as the generator here makes them, and the calls.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::PathBuf;

use indexmap::IndexMap;
use og_cn::catalog::{COLUMNS, CallRule, Catalog, Kind, Locus, PSV_FINDINGS, parse_row};
use og_cn::engine::{AssessError, Block, Means, Reason, Sex, Status, assess};
use pycompat::json::{self, Value};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

fn repo(rel: &str) -> String {
    format!("{REPO}/{rel}")
}

fn catalog_path() -> String {
    repo("bin/cn_loci.tsv")
}

fn cat() -> Catalog {
    Catalog::load(&catalog_path()).expect("bin/cn_loci.tsv loads")
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_default()
}

/// A scratch folder of its own, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!("og-genome-cn-{name}-{}", std::process::id()));
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

/// `og cn ARGS`: (exit code, stdout, stderr).
fn og_cn(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_cn::cli::main(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

/// The catalog's rows as written, by locus: column -> cell.
fn raw_rows() -> IndexMap<String, HashMap<String, String>> {
    read(&catalog_path())
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .skip(1)
        .map(|l| {
            let m: HashMap<String, String> = COLUMNS
                .iter()
                .map(|c| (*c).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect();
            (m["locus"].clone(), m)
        })
        .collect()
}

/// A made-up bedcov row: one 1 kb gene between two 10 kb flanks, called by
/// rounding, with a MAPQ 0 ceiling of 0.5; `kw` replaces columns.
fn row(kw: &[(&str, &str)]) -> HashMap<String, String> {
    let mut r: HashMap<String, String> = [
        ("locus", "TESTG"),
        ("kind", "gene"),
        ("chrom", "chr3"),
        ("targets", "TESTG:50000-51000"),
        ("controls", "fl_a:30000-40000;fl_b:60000-70000"),
        ("details", "."),
        ("target_depth", "sum"),
        ("control_depth", "flanks"),
        ("scale", "2"),
        ("offset", "0"),
        ("bedcov", ".,-Q 1"),
        ("mapq0_min", "."),
        ("mapq0_max", "0.5"),
        ("min_control_depth", "8"),
        ("max_control_ratio", "1.5"),
        ("call", "round:0.3"),
        ("interval", "."),
        ("readings", "0=none|1=one|2=two|more={n} copies"),
        ("input", "bedcov"),
        ("owner", "cn-depth"),
        ("cross_check", "none"),
        ("pmids", "."),
    ]
    .iter()
    .map(|(k, v)| ((*k).to_owned(), (*v).to_owned()))
    .collect();
    for (k, v) in kw {
        r.insert((*k).to_owned(), (*v).to_owned());
    }
    r
}

fn line(r: &HashMap<String, String>) -> String {
    COLUMNS
        .iter()
        .map(|c| r[*c].as_str())
        .collect::<Vec<_>>()
        .join("\t")
}

fn write_catalog(path: &str, rows: &[HashMap<String, String>]) {
    let mut text = COLUMNS.join("\t") + "\n";
    for r in rows {
        text += &(line(r) + "\n");
    }
    std::fs::write(path, text).unwrap();
}

fn parsed(r: &HashMap<String, String>) -> Locus {
    let l = line(r);
    let cells: Vec<&str> = l.split('\t').collect();
    parse_row(&cells).unwrap_or_else(|e| panic!("{e}"))
}

/// `Catalog::load` on these rows: the loaded catalog's size, or its error.
fn load_rows(w: &Scratch, rows: &[HashMap<String, String>]) -> Result<usize, String> {
    let p = w.path("c.tsv");
    write_catalog(&p, rows);
    Catalog::load(&p).map(|c| c.loci.len()).map_err(|e| e.0)
}

fn means(pairs: &[(&str, f64)]) -> Means {
    pairs.iter().map(|(k, v)| ((*k).to_owned(), *v)).collect()
}

/// Python's `round()` of a float to a whole number (ties to even).
fn round(v: f64) -> i64 {
    v.round_ties_even() as i64
}

/// samtools bedcov over a row's segments at these mean depths: whole sums,
/// or with `exact` the sums themselves, so each mean comes back as given.
fn bedcov(path: &str, row: &Locus, depth: &Means, exact: bool) {
    let mut text = String::new();
    for (c, s, e, n) in row.regions() {
        let v = depth[n] * (e - s) as f64;
        if exact {
            text += &format!("{c}\t{s}\t{e}\t{n}\t{v:?}\n");
        } else {
            text += &format!("{c}\t{s}\t{e}\t{n}\t{}\n", round(v));
        }
    }
    std::fs::write(path, text).unwrap();
}

/// A `metric<TAB>value` table, in its order; empty when it is not there.
fn table(path: &str) -> Vec<(String, String)> {
    read(path)
        .lines()
        .skip(1)
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_owned(), v.to_owned()))
        .collect()
}

fn get<'a>(t: &'a [(String, String)], k: &str) -> &'a str {
    t.iter()
        .find(|(m, _)| m == k)
        .map_or("", |(_, v)| v.as_str())
}

fn status(b: &Block) -> String {
    b.status.to_string()
}

/// The change column of an interval row for a block.
fn change(b: &Block) -> &'static str {
    match (b.copies, b.expected) {
        (Some(c), Some(e)) if c < e => "deletion",
        (Some(c), Some(e)) if c > e => "duplication",
        (Some(_), Some(_)) => "no change",
        _ => ".",
    }
}

// 1. the catalog

const DOSAGE: [&str; 3] = ["HP", "LCE3BC", "CCL3L1"];

/// `-X` or `-X N` options, one or more, space separated.
fn samtools_options(p: &str) -> bool {
    let mut words = p.split(' ').peekable();
    let mut any = false;
    while let Some(w) = words.next() {
        let b = w.as_bytes();
        if !(b.len() == 2 && b[0] == b'-' && b[1].is_ascii_alphabetic()) {
            return false;
        }
        any = true;
        if words
            .peek()
            .is_some_and(|n| !n.is_empty() && n.bytes().all(|c| c.is_ascii_digit()))
        {
            words.next();
        }
    }
    any
}

/// `none` or `tool:Name[,tool:Name...]`.
fn cross_check_ok(x: &str) -> bool {
    x == "none"
        || x.split(',').all(|p| {
            p.split_once(':').is_some_and(|(a, b)| {
                !a.is_empty()
                    && a.bytes().all(|c| c.is_ascii_lowercase())
                    && !b.is_empty()
                    && b.bytes().all(|c| c.is_ascii_alphanumeric() || c == b'_')
            })
        })
}

/// `.` or comma-separated numbers.
fn pmids_ok(x: &str, dot: bool) -> bool {
    (dot && x == ".")
        || x.split(',')
            .all(|n| !n.is_empty() && n.bytes().all(|c| c.is_ascii_digit()))
}

#[test]
fn the_gene_dosage_and_cyp21a2_rows_come_first_then_the_interval_rows() {
    let c = cat();
    let ids: Vec<&str> = c.loci.keys().map(String::as_str).collect();
    assert_eq!(
        &ids[..7],
        [
            "HBA", "LPA_KIV2", "RHD", "HP", "LCE3BC", "CCL3L1", "CYP21A2"
        ]
    );
    for r in c.loci.values().skip(7) {
        assert!(
            r.kind == Kind::Interval && r.depth().passes == ["."],
            "{}",
            r.id
        );
    }
}

#[test]
fn every_catalog_row_is_well_formed() {
    let c = cat();
    let raw = raw_rows();
    for (id, r) in &c.loci {
        let w = &raw[id];
        assert_eq!(w["owner"], "cn-depth", "{id}");
        assert!(
            cross_check_ok(&w["cross_check"]),
            "{id}: {}",
            w["cross_check"]
        );
        assert!(pmids_ok(&w["pmids"], true), "{id}: {}", w["pmids"]);
        if w["input"] != "bedcov" {
            continue;
        }
        let passes = &r.depth().passes;
        assert!(
            passes.iter().all(|p| p == "." || samtools_options(p)),
            "{id}: {passes:?}"
        );
        if passes.len() > 1 {
            let strip = |p: &str| p.replace('.', "").replace("-Q 1", "").trim().to_owned();
            assert_eq!(
                strip(&passes[1]),
                passes[0].replace('.', "").trim(),
                "{id}: a second pass adds -Q 1"
            );
        }
        for s in r.targets.iter().chain(&r.controls).chain(&r.details) {
            assert!(
                0 < s.end - s.start && s.end - s.start < 1_000_000,
                "{id}: {}",
                s.name
            );
        }
        for ctl in &r.controls {
            assert!(
                !r.targets
                    .iter()
                    .any(|t| ctl.start < t.end && t.start < ctl.end),
                "{id}: control {} overlaps a target",
                ctl.name
            );
        }
        let want: Option<Vec<String>> = match &r.call {
            CallRule::Class { .. } => Some(["low", "typical", "high"].map(str::to_owned).to_vec()),
            CallRule::Cut(cuts) => Some((0..=cuts.len()).map(|i| i.to_string()).collect()),
            _ => None,
        };
        if let Some(want) = want {
            for k in want {
                assert!(r.readings.contains_key(&k), "{id}: no reading for {k}");
            }
        }
    }
}

/// A broken row: what, the changed columns, a part of the error.
type Broken = (
    &'static str,
    Vec<(&'static str, &'static str)>,
    &'static str,
);

/// The broken rows.
fn broken() -> Vec<Broken> {
    vec![
        ("a bad kind", vec![("kind", "exon")], "kind"),
        (
            "a segment with no span",
            vec![("targets", "TESTG:51000-50000")],
            "empty span",
        ),
        (
            "a class call without an interval",
            vec![("call", "class:30,45")],
            "needs interval cv",
        ),
        (
            "a MAPQ 0 guard with one bedcov pass",
            vec![("bedcov", ".")],
            "two bedcov passes",
        ),
        (
            "a segment name used twice",
            vec![("details", "fl_a:1-2")],
            "used twice",
        ),
        (
            "a NaN call threshold",
            vec![("call", "round:nan")],
            "is not round",
        ),
        (
            "an infinite cut threshold",
            vec![("call", "cut:0.5,inf")],
            "is not round",
        ),
        (
            "a NaN number column",
            vec![("mapq0_max", "nan")],
            "not a number",
        ),
    ]
}

#[test]
fn broken_rows_are_refused_by_the_reader_and_by_og_cn_bed() {
    let w = Scratch::new("broken");
    assert_eq!(
        load_rows(&w, &[row(&[])]),
        Ok(1),
        "the plain made-up row loads"
    );
    for (desc, kw, err) in broken() {
        let got = load_rows(&w, &[row(&kw)]);
        let msg = got.clone().err().unwrap_or_else(|| panic!("{desc} loaded"));
        assert!(msg.contains(err), "{desc}: {msg}");
        let (rc, _, stderr) = og_cn(&["bed", "--catalog", &w.path("c.tsv")]);
        assert_eq!(
            (rc, stderr),
            (2, format!("ERROR: {msg}\n")),
            "{desc} through og cn bed"
        );
    }
    let two = row(&[
        ("locus", "TESTH"),
        ("targets", "TESTG:80000-81000"),
        ("controls", "fl_c:1-2;fl_d:3-4"),
    ]);
    let got = load_rows(&w, &[row(&[]), two]);
    assert!(
        got.clone().is_err_and(|e| e.contains("used by two loci")),
        "{got:?}"
    );
}

// 2. one owner per locus

/// The only rows allowed to overlap a Parascopy model, each for a stated
/// reason: step 41 runs by default and Parascopy is opt-in, so the RhD
/// deletion stays a depth row; RHCE's model holds RHD as a paralog copy and
/// is its cross-check. A cross_check in the catalog does not exempt a row.
const EXEMPT: [(&str, &str); 1] = [("RHD", "RHCE")];

type Models = BTreeMap<String, Vec<(String, i64, i64, String)>>;

fn models() -> Models {
    let mut out: Models = BTreeMap::new();
    for l in read(&repo("assets/parascopy_models_v1.7.txt")).lines() {
        if l.starts_with('#') || l.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = l.split_whitespace().collect();
        let (c, span) = f[1].split_once(':').unwrap();
        let (s, e) = span.split_once('-').unwrap();
        out.entry(f[0].to_owned()).or_default().push((
            c.to_owned(),
            s.parse().unwrap(),
            e.parse().unwrap(),
            f[2].to_owned(),
        ));
    }
    out
}

/// [(model, segment)] of a row's targets and controls that overlap any
/// region of a model (1-based, both ends included, against 0-based
/// segments).
fn overlaps(r: &Locus, models: &Models) -> Vec<(String, String)> {
    let mut hits = Vec::new();
    for s in r.targets.iter().chain(&r.controls) {
        for (name, regs) in models {
            if regs
                .iter()
                .any(|(c, ms, me, _)| *c == r.chrom && s.start < *me && ms - 1 < s.end)
            {
                hits.push((name.clone(), s.name.clone()));
            }
        }
    }
    hits
}

fn owner_violations(rows: &[&Locus], models: &Models) -> Vec<(String, String, String)> {
    rows.iter()
        .flat_map(|r| {
            overlaps(r, models)
                .into_iter()
                .map(|(m, s)| (r.id.clone(), m, s))
        })
        .filter(|(l, m, _)| !EXEMPT.contains(&(l.as_str(), m.as_str())))
        .collect()
}

#[test]
fn the_parascopy_model_list_is_the_pinned_one() {
    let m = models();
    assert_eq!(m.len(), 167);
    for g in [
        "AMY1C", "C4A", "CFH", "CFHR4", "FCGR3A", "CHRNA7", "APOBEC3A", "GYPA", "RHCE",
    ] {
        assert!(m.contains_key(g), "{g}");
    }
    assert!(
        m.values()
            .all(|r| r.iter().filter(|x| x.3 == "main").count() == 1
                && r.iter().all(|x| x.3 == "main" || x.3 == "copy"))
    );
    assert!(m.values().map(Vec::len).sum::<usize>() > 1000);
}

#[test]
fn no_row_counts_a_gene_a_pinned_parascopy_model_owns() {
    let (c, m) = (cat(), models());
    let raw = raw_rows();
    for (id, w) in &raw {
        for x in w["cross_check"].split(',') {
            if let Some(model) = x.strip_prefix("parascopy:") {
                assert!(
                    m.contains_key(model),
                    "{id}: cross_check {x} is not a pinned model"
                );
            }
        }
    }
    let rows: Vec<&Locus> = c.loci.values().collect();
    assert_eq!(owner_violations(&rows, &m), Vec::new());
    for (l, model) in EXEMPT {
        let r = c
            .loci
            .get(l)
            .unwrap_or_else(|| panic!("the exempt row {l} is gone"));
        assert!(
            overlaps(r, &m).iter().any(|(x, _)| x == model),
            "the exemption {l}/{model} is stale"
        );
    }
    let rhd: HashSet<String> = overlaps(&c.loci["RHD"], &m)
        .into_iter()
        .map(|(x, _)| x)
        .collect();
    assert_eq!(
        rhd,
        HashSet::from(["RHCE".to_owned()]),
        "RHD overlaps RHCE's model only"
    );
}

#[test]
fn planted_rows_on_a_parascopy_model_are_caught() {
    let (c, m) = (cat(), models());
    for (kw, want) in [
        (
            vec![
                ("locus", "AMY1"),
                ("chrom", "chr1"),
                ("targets", "AMY1C:103745500-103760000"),
                (
                    "controls",
                    "amy_l:103000000-103010000;amy_r:104000000-104010000",
                ),
            ],
            ("AMY1", "AMY1C", "AMY1C"),
        ),
        (
            vec![
                ("locus", "C4B"),
                ("chrom", "chr6"),
                ("targets", "C4B:32014800-32035400"),
                ("controls", "c4_l:31000000-31010000;c4_r:33500000-33510000"),
            ],
            ("C4B", "C4A", "C4B"),
        ),
        (
            vec![
                ("locus", "GYPA"),
                ("chrom", "chr4"),
                ("cross_check", "parascopy:GYPA"),
                ("targets", "GYPA_T:144140000-144150000"),
                (
                    "controls",
                    "gy_l:140000000-140010000;gy_r:150000000-150010000",
                ),
            ],
            ("GYPA", "GYPA", "GYPA_T"),
        ),
        // DEFB4 has no row: the beta-defensin repeat lies in SPAG11A's copies.
        (
            vec![
                ("locus", "DEFB4"),
                ("chrom", "chr8"),
                ("targets", "DEFB_UNIT_B:7414854-7430348"),
                ("controls", "defb_l:7622000-7662000;defb_r:8250000-8290000"),
            ],
            ("DEFB4", "SPAG11A", "DEFB_UNIT_B"),
        ),
    ] {
        let planted = parsed(&row(&kw));
        let mut rows: Vec<&Locus> = c.loci.values().collect();
        rows.push(&planted);
        let v = owner_violations(&rows, &m);
        let want = (want.0.to_owned(), want.1.to_owned(), want.2.to_owned());
        assert!(v.contains(&want), "{want:?} not caught: {v:?}");
        if want.0 != "GYPA" {
            assert_eq!(v, vec![want]);
        }
    }
}

// 3. parity with the code before the engine

fn gold() -> Value {
    json::loads(&read(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/genome_gold.json"
    )))
    .unwrap()
}

fn num(v: &Value) -> f64 {
    match v {
        Value::Int(i) => *i as f64,
        Value::Float(f) => *f,
        _ => panic!("{v:?}"),
    }
}

fn arr(v: &Value) -> &[Value] {
    match v {
        Value::Array(a) => a,
        _ => panic!("{v:?}"),
    }
}

/// Every field of the frozen row: numbers as numbers, text as it is.
fn same(got: &[(String, String)], want: &Value) -> Result<(), String> {
    let Value::Object(want) = want else { panic!() };
    for (k, v) in want {
        let g = get(got, k);
        let ok = match v {
            Value::Str(s) => g == s,
            _ => g.parse::<f64>().is_ok_and(|x| (x - num(v)).abs() < 1e-9),
        };
        if !ok {
            return Err(format!("{k}: got {g:?}, want {v:?}"));
        }
    }
    Ok(())
}

#[test]
fn parity_hba_tables_are_the_ones_before_the_engine() {
    let (c, w) = (cat(), Scratch::new("gold-hba"));
    let (a, q, tab) = (w.path("a.tsv"), w.path("q.tsv"), w.path("t.tsv"));
    let hrow = &c.loci["HBA"];
    let g = gold();
    let Value::Object(g) = &g else { panic!() };
    let cases = arr(&g["hba"]);
    assert_eq!(cases.len(), 12);
    for case in cases {
        let (inp, want) = (arr(&arr(case)[0]), &arr(case)[1]);
        let [h2, h1, left, right, share] = [0, 1, 2, 3, 4].map(|i| num(&inp[i]));
        let m = means(&[
            ("flank_left", left),
            ("HBA2", h2),
            ("HBA1", h1),
            ("flank_right", right),
        ]);
        let mq: Means = m
            .iter()
            .map(|(n, v)| {
                (
                    n.clone(),
                    if n.starts_with("HBA") {
                        v * (1.0 - share)
                    } else {
                        *v
                    },
                )
            })
            .collect();
        bedcov(&a, hrow, &m, true);
        bedcov(&q, hrow, &mq, true);
        let (rc, _, err) = og_cn(&[
            "call",
            "--catalog",
            &catalog_path(),
            "--locus",
            "HBA",
            "--bedcov",
            &a,
            "--bedcov-mapq1",
            &q,
            "--table",
            &tab,
        ]);
        assert_eq!(rc, 0, "{err}");
        if let Err(e) = same(&table(&tab), want) {
            panic!("HBA {h2}/{h1} flanks {left}/{right} MAPQ0 {share}: {e}");
        }
    }
}

#[test]
fn parity_lpa_tables_are_the_ones_before_the_engine() {
    let (c, w) = (cat(), Scratch::new("gold-lpa"));
    let (a, q, tab) = (w.path("a.tsv"), w.path("q.tsv"), w.path("t.tsv"));
    // The code before the engine had no systematic term: the parity runs on
    // the catalog with LPA_KIV2's interval as plain cv (the term itself is
    // checked in genome_lpa.rs and og-cn's engine tests).
    let plain_cv = w.path("cn_loci.tsv");
    let col = COLUMNS.iter().position(|c| *c == "interval").unwrap();
    let text: String = read(&catalog_path())
        .lines()
        .map(|l| {
            let mut f: Vec<&str> = l.split('\t').collect();
            if f[0] == "LPA_KIV2" {
                assert!(f[col].starts_with("cv:"), "LPA_KIV2 has cv:SYS");
                f[col] = "cv";
            }
            f.join("\t") + "\n"
        })
        .collect();
    std::fs::write(&plain_cv, text).unwrap();
    let lrow = &c.loci["LPA_KIV2"];
    let g = gold();
    let Value::Object(g) = &g else { panic!() };
    let cases = arr(&g["lpa"]);
    assert_eq!(cases.len(), 11);
    for case in cases {
        let (inp, want) = (arr(&arr(case)[0]), &arr(case)[1]);
        let [n, flank, ks, fs, share] = [0, 1, 2, 3, 4].map(|i| num(&inp[i]));
        let k = flank * (n + 2.0) / 14.0;
        let mut m = Means::new();
        for (i, s) in lrow.targets.iter().enumerate() {
            m.insert(
                s.name.clone(),
                k * (1.0 + if i % 2 == 0 { ks } else { -ks }),
            );
        }
        for (i, s) in lrow.controls.iter().enumerate() {
            m.insert(
                s.name.clone(),
                flank * (1.0 + if i % 2 == 0 { fs } else { -fs }),
            );
        }
        let mq: Means = m
            .iter()
            .map(|(x, v)| {
                (
                    x.clone(),
                    if x.starts_with("kiv2_") {
                        v * (1.0 - share)
                    } else {
                        *v
                    },
                )
            })
            .collect();
        bedcov(&a, lrow, &m, true);
        bedcov(&q, lrow, &mq, true);
        let (rc, _, err) = og_cn(&[
            "call",
            "--catalog",
            &plain_cv,
            "--locus",
            "LPA_KIV2",
            "--bedcov",
            &a,
            "--bedcov-mapq1",
            &q,
            "--table",
            &tab,
        ]);
        assert_eq!(rc, 0, "{err}");
        if let Err(e) = same(&table(&tab), want) {
            panic!("LPA {n} copies at {flank}x: {e}");
        }
    }
}

// 4. synthetic copies, the MAPQ 0 ceiling, a kind=interval row

#[test]
fn synthetic_depths_at_0_to_4_copies_are_called_right() {
    let g = parsed(&row(&[]));
    for copies in 0..=4 {
        let m = means(&[
            ("fl_a", 30.0),
            ("fl_b", 30.0),
            ("TESTG", 15.0 * copies as f64),
        ]);
        let q: Means = m.iter().map(|(n, v)| (n.clone(), v * 0.9)).collect();
        let b = assess(&g, &m, Some(&q), Sex::Unknown).unwrap();
        assert_eq!(
            (b.status, b.copies),
            (Status::Called, Some(copies)),
            "{copies} copies"
        );
    }
}

#[test]
fn the_mapq0_ceiling_between_counts_and_no_control_depth() {
    let g = parsed(&row(&[]));
    let m = means(&[("fl_a", 30.0), ("fl_b", 30.0), ("TESTG", 30.0)]);
    let b = assess(
        &g,
        &m,
        Some(&means(&[("fl_a", 30.0), ("fl_b", 30.0), ("TESTG", 12.0)])),
        Sex::Unknown,
    )
    .unwrap();
    assert_eq!(
        (b.status, b.reason),
        (Status::NotAssessed, Some(Reason::HighMapq0Share)),
        "60% from MAPQ 0"
    );
    let b = assess(
        &g,
        &m,
        Some(&means(&[("fl_a", 30.0), ("fl_b", 30.0), ("TESTG", 16.0)])),
        Sex::Unknown,
    )
    .unwrap();
    assert_eq!(
        (b.status, b.copies),
        (Status::Called, Some(2)),
        "47% from MAPQ 0"
    );
    let m25 = means(&[("fl_a", 30.0), ("fl_b", 30.0), ("TESTG", 37.5)]);
    let b = assess(&g, &m25, Some(&m25), Sex::Unknown).unwrap();
    assert_eq!(b.reason, Some(Reason::BetweenCounts), "2.5 copies");
    // No minimum control depth set: no depth, or one empty flank, is a reason, not a crash.
    let z = parsed(&row(&[("min_control_depth", ".")]));
    let m0 = means(&[("fl_a", 0.0), ("fl_b", 0.0), ("TESTG", 30.0)]);
    let b = assess(&z, &m0, Some(&m0), Sex::Unknown).unwrap();
    assert_eq!(b.reason, Some(Reason::LowControlDepth));
    let m1 = means(&[("fl_a", 0.0), ("fl_b", 30.0), ("TESTG", 30.0)]);
    let b = assess(&z, &m1, Some(&m1), Sex::Unknown).unwrap();
    assert_eq!(b.reason, Some(Reason::ControlsDisagree));
}

#[test]
fn an_interval_row_with_pooled_depths_and_the_mapq1_depths_rule() {
    let iv = parsed(&row(&[
        ("kind", "interval"),
        ("target_depth", "pooled"),
        ("control_depth", "pooled"),
        ("mapq0_max", "."),
        ("bedcov", "."),
        ("targets", "seg1:50000-50500;seg2:50500-51000"),
    ]));
    let m = means(&[
        ("fl_a", 20.0),
        ("fl_b", 40.0),
        ("seg1", 45.0),
        ("seg2", 45.0),
    ]);
    let b = assess(&iv, &m, None, Sex::Unknown).unwrap();
    assert_eq!(
        (b.status, b.copies),
        (Status::Called, Some(3)),
        "pooled 45x over pooled 30x"
    );
    let c = cat();
    let rhd = &c.loci["RHD"];
    let m: Means = rhd
        .regions()
        .iter()
        .map(|r| (r.3.to_owned(), 30.0))
        .collect();
    let b = assess(rhd, &m, None, Sex::Unknown).unwrap();
    assert!(
        b.status == Status::Called && b.mapq0_share.is_none(),
        "{b:?}"
    );
    let g = parsed(&row(&[]));
    let m = means(&[("fl_a", 30.0), ("fl_b", 30.0), ("TESTG", 37.5)]);
    match assess(&g, &m, None, Sex::Unknown) {
        Err(AssessError::NeedsMapq1(e)) => {
            assert!(e.contains("MAPQ >= 1 depths are needed"), "{e}")
        }
        other => panic!("{other:?}"),
    }
}

// 5. a broken control

#[test]
fn a_control_with_no_read_is_broken_and_no_read_anywhere_is_not_assessed() {
    let g = parsed(&row(&[]));
    let bad = means(&[("fl_a", 0.0), ("fl_b", 30.0), ("TESTG", 30.0)]);
    match assess(&g, &bad, Some(&bad), Sex::Unknown) {
        Err(AssessError::BrokenControl(e)) => {
            assert!(
                e.starts_with("BROKEN_CONTROL: control segment fl_a of TESTG"),
                "{e}"
            );
        }
        other => panic!("{other:?}"),
    }
    let b = assess(
        &g,
        &means(&[("fl_a", 0.0), ("fl_b", 0.0), ("TESTG", 0.0)]),
        Some(&Means::new()),
        Sex::Unknown,
    )
    .unwrap();
    assert_eq!(
        (b.status, b.reason),
        (Status::NotAssessed, Some(Reason::LowControlDepth))
    );
    let b = assess(
        &g,
        &means(&[("fl_a", 0.0), ("fl_b", 5.0), ("TESTG", 5.0)]),
        Some(&Means::new()),
        Sex::Unknown,
    )
    .unwrap();
    assert_eq!(
        b.reason,
        Some(Reason::LowControlDepth),
        "one flank empty beside a thin one"
    );
}

#[test]
fn call_exits_3_on_a_broken_control_and_writes_nothing() {
    let (c, w) = (cat(), Scratch::new("broken-call"));
    let lrow = &c.loci["LPA_KIV2"];
    let mut d: Means = lrow
        .regions()
        .iter()
        .map(|r| (r.3.to_owned(), 30.0))
        .collect();
    d.insert("flank_e19".into(), 0.0);
    let (a, q, tab) = (w.path("a.tsv"), w.path("q.tsv"), w.path("t.tsv"));
    bedcov(&a, lrow, &d, false);
    bedcov(&q, lrow, &d, false);
    let (rc, _, err) = og_cn(&[
        "call",
        "--catalog",
        &catalog_path(),
        "--locus",
        "LPA_KIV2",
        "--bedcov",
        &a,
        "--bedcov-mapq1",
        &q,
        "--table",
        &tab,
    ]);
    assert_eq!(rc, 3, "{err}");
    assert!(
        err.contains("BROKEN_CONTROL: control segment flank_e19 of LPA_KIV2"),
        "{err}"
    );
    assert!(!std::path::Path::new(&tab).exists());
}

// 6. the command line

#[test]
fn the_bedcov_reader_refuses_a_region_twice_or_an_empty_region() {
    let w = Scratch::new("reader");
    let p = w.path("c.tsv");
    write_catalog(&p, &[row(&[])]);
    let g = parsed(&row(&[]));
    let b = w.path("dup.bedcov");
    bedcov(
        &b,
        &g,
        &means(&[("fl_a", 30.0), ("fl_b", 30.0), ("TESTG", 30.0)]),
        false,
    );
    let first = read(&b).lines().next().unwrap().to_owned();
    let mut text = read(&b);
    text += &(first.rsplit_once('\t').unwrap().0.to_owned() + "\t0\n");
    std::fs::write(&b, text).unwrap();
    let (rc, _, err) = og_cn(&[
        "call",
        "--catalog",
        &p,
        "--locus",
        "TESTG",
        "--bedcov",
        &b,
        "--bedcov-mapq1",
        &b,
    ]);
    assert!(rc == 2 && err.contains("appears twice"), "{rc} {err}");
    std::fs::write(
        &b,
        format!(
            "{}\t50000\t50000\tTESTG\t0\n",
            first.split('\t').next().unwrap()
        ),
    )
    .unwrap();
    let (rc, _, err) = og_cn(&[
        "call",
        "--catalog",
        &p,
        "--locus",
        "TESTG",
        "--bedcov",
        &b,
        "--bedcov-mapq1",
        &b,
    ]);
    assert!(
        rc == 2 && err.contains("is empty or reversed"),
        "{rc} {err}"
    );
}

#[test]
fn bed_prints_every_segment_of_the_bedcov_loci_once() {
    let c = cat();
    let (rc, out, err) = og_cn(&["bed", "--catalog", &catalog_path()]);
    assert_eq!(rc, 0, "{err}");
    let names: Vec<&str> = out.lines().map(|l| l.split('\t').nth(3).unwrap()).collect();
    let uniq: HashSet<&&str> = names.iter().collect();
    let want: usize = c
        .loci
        .values()
        .filter(|r| r.input == og_cn::catalog::Input::Bedcov)
        .map(|r| r.regions().len())
        .sum();
    assert_eq!((names.len(), uniq.len()), (want, want));
}

/// Step 46's table, in its order.
const LPA_ORDER: [&str; 10] = [
    "status",
    "kiv2_copies",
    "ci95_low",
    "ci95_high",
    "class",
    "kiv2_depth",
    "flank_depth",
    "mapq0_share",
    "message",
    "advice",
];

#[test]
fn call_writes_step_46s_table_needs_mapq1_and_reads_empty_tables() {
    let (c, w) = (cat(), Scratch::new("cli"));
    let lrow = &c.loci["LPA_KIV2"];
    let d: Means = lrow
        .regions()
        .iter()
        .map(|r| {
            (
                r.3.to_owned(),
                if r.3.starts_with("kiv2_") {
                    30.0 * 22.0 / 14.0
                } else {
                    30.0
                },
            )
        })
        .collect();
    let q: Means = d
        .iter()
        .map(|(n, v)| (n.clone(), if n.starts_with("kiv2_") { v * 0.1 } else { *v }))
        .collect();
    let (a, b, tab) = (w.path("a.tsv"), w.path("q.tsv"), w.path("t.tsv"));
    bedcov(&a, lrow, &d, false);
    bedcov(&b, lrow, &q, false);
    let cp = catalog_path();
    let (rc, _, err) = og_cn(&[
        "call",
        "--catalog",
        &cp,
        "--locus",
        "LPA_KIV2",
        "--bedcov",
        &a,
        "--bedcov-mapq1",
        &b,
        "--table",
        &tab,
    ]);
    let t = table(&tab);
    assert_eq!(rc, 0, "{err}");
    assert_eq!(
        (get(&t, "status"), get(&t, "kiv2_copies")),
        ("estimated", "20.0")
    );
    assert_eq!(
        t.iter().map(|(k, _)| k.as_str()).collect::<Vec<_>>(),
        LPA_ORDER
    );
    let (rc, _, err) = og_cn(&[
        "call",
        "--catalog",
        &cp,
        "--locus",
        "LPA_KIV2",
        "--bedcov",
        &a,
    ]);
    assert!(rc == 2 && err.contains("--bedcov-mapq1"), "{rc} {err}");
    let empty = w.put("empty.tsv", "");
    for (loc, key) in [("HBA", "alpha_globin_genes"), ("LPA_KIV2", "kiv2_copies")] {
        let _ = std::fs::remove_file(&tab);
        let (rc, _, err) = og_cn(&[
            "call",
            "--catalog",
            &cp,
            "--locus",
            loc,
            "--bedcov",
            &empty,
            "--bedcov-mapq1",
            &empty,
            "--table",
            &tab,
        ]);
        let t = table(&tab);
        assert_eq!(rc, 0, "{loc}: {err}");
        assert_eq!(
            (get(&t, "status"), get(&t, key)),
            ("not assessed", "."),
            "{loc}"
        );
        assert!(get(&t, "message").contains("no chr"), "{loc}: {t:?}");
    }
}

// 8. the gene-dosage rows

/// Mean depths of a dosage row at a copy number: controls at 30x, the
/// target copies spread evenly over the row's targets, and the MAPQ >= 1
/// depths with half of the targets'.
fn dosage_depths(r: &Locus, copies: f64) -> (Means, Means) {
    let tg: HashSet<&str> = r.targets.iter().map(|s| s.name.as_str()).collect();
    let m: Means = r
        .regions()
        .iter()
        .map(|x| {
            (
                x.3.to_owned(),
                if tg.contains(x.3) {
                    15.0 * copies / tg.len() as f64
                } else {
                    30.0
                },
            )
        })
        .collect();
    let q = m
        .iter()
        .map(|(n, v)| {
            (
                n.clone(),
                if tg.contains(n.as_str()) { v * 0.5 } else { *v },
            )
        })
        .collect();
    (m, q)
}

/// Words a reading of CCL3L1 must not carry: its disease links are disputed.
const DENY: [&str; 16] = [
    "hiv",
    "aids",
    "infect",
    "disease",
    "risk",
    "suscept",
    "progress",
    "malaria",
    "hepatitis",
    "arthritis",
    "lupus",
    "kawasaki",
    "psoria",
    "diabet",
    "crohn",
    "vaccin",
];

fn denied(text: &str) -> Vec<&'static str> {
    let low = text.to_lowercase();
    DENY.iter().copied().filter(|w| low.contains(w)).collect()
}

#[test]
fn the_dosage_rows_count_0_to_4_copies_with_their_readings() {
    let (c, raw) = (cat(), raw_rows());
    for loc in DOSAGE {
        let r = &c.loci[loc];
        assert_eq!(
            r.depth().passes,
            [".", "-Q 1"],
            "{loc}: all reads, then MAPQ >= 1"
        );
        assert!(matches!(r.call, CallRule::Round { .. }), "{loc}");
        for k in ["0", "1", "more"] {
            assert!(r.readings.contains_key(k), "{loc}: no reading for {k}");
        }
        assert!(pmids_ok(&raw[loc]["pmids"], false), "{loc}: a PMID");
        for copies in 0..=4 {
            let (m, q) = dosage_depths(r, copies as f64);
            let b = assess(r, &m, Some(&q), Sex::Unknown).unwrap();
            assert_eq!(
                (b.status, b.copies),
                (Status::Called, Some(copies)),
                "{loc} {copies}"
            );
            assert_eq!(b.reading, r.reading_count(copies), "{loc} {copies}");
            assert!(!b.reading.contains("{n}"), "{loc}: {}", b.reading);
        }
        let (m, q) = dosage_depths(r, 2.5);
        assert_eq!(
            assess(r, &m, Some(&q), Sex::Unknown).unwrap().reason,
            Some(Reason::BetweenCounts),
            "{loc}"
        );
    }
    let hp = &c.loci["HP"];
    let (mut m, q) = dosage_depths(hp, 4.0);
    m.insert("HP_UNIT1".into(), 60.0);
    m.insert("HP_UNIT2".into(), 0.0);
    assert_eq!(
        assess(hp, &m, Some(&q), Sex::Unknown).unwrap().copies,
        Some(4),
        "both copies on one reference copy"
    );
}

#[test]
fn no_ccl3l1_reading_names_a_disease() {
    let c = cat();
    let readings: Vec<&str> = c.loci["CCL3L1"]
        .readings
        .values()
        .map(String::as_str)
        .collect();
    let all = readings.join(" ");
    assert_eq!(denied(&all), Vec::<&str>::new());
    assert!(
        !denied(&(all + " more copies go with slower HIV progression")).is_empty(),
        "control: a planted HIV"
    );
}

#[test]
fn a_dosage_row_shifted_by_1_kb_refuses_its_fixture() {
    let w = Scratch::new("shift");
    let (c, raw) = (cat(), raw_rows());
    let p = w.path("c.tsv");
    for loc in DOSAGE {
        let mut r = raw[loc].clone();
        r.insert(
            "targets".into(),
            c.loci[loc]
                .targets
                .iter()
                .map(|s| format!("{}:{}-{}", s.name, s.start + 1000, s.end + 1000))
                .collect::<Vec<_>>()
                .join(";"),
        );
        write_catalog(&p, &[r]);
        let f = repo(&format!("tests/fixtures/cn/{loc}.bedcov.tsv"));
        let (rc, _, err) = og_cn(&[
            "call",
            "--catalog",
            &p,
            "--locus",
            loc,
            "--bedcov",
            &f,
            "--bedcov-mapq1",
            &f,
        ]);
        assert!(
            rc == 2 && err.contains("regions differ"),
            "{loc}: {rc} {err}"
        );
        // and the unshifted row reads it
        let (rc, _, err) = og_cn(&[
            "call",
            "--catalog",
            &catalog_path(),
            "--locus",
            loc,
            "--bedcov",
            &f,
            "--bedcov-mapq1",
            &f,
        ]);
        assert_eq!(rc, 0, "{loc}: {err}");
    }
}

// 9. the interval rows

const HAND: [&str; 10] = [
    "1p36_terminal",
    "Xq28_MECP2",
    "GJB6_D13S1830",
    "CFTR_dele2_3",
    "HBB_exon3_del",
    "HBB_619bp_del",
    "HBD_HBB_Sicilian_del",
    "AZFa",
    "AZFb",
    "AZFc_DAZ",
];

fn interval_rows(c: &Catalog) -> Vec<&Locus> {
    c.loci
        .values()
        .filter(|r| r.kind == Kind::Interval && r.depth().passes == ["."])
        .collect()
}

/// Targets at `target`, every other segment at `control`.
fn depths(r: &Locus, target: f64, control: f64) -> Means {
    let tg: HashSet<&str> = r.targets.iter().map(|s| s.name.as_str()).collect();
    r.regions()
        .iter()
        .map(|x| {
            (
                x.3.to_owned(),
                if tg.contains(x.3) { target } else { control },
            )
        })
        .collect()
}

#[test]
fn the_interval_rows_are_the_clingen_regions_then_the_hand_rows() {
    let (c, raw) = (cat(), raw_rows());
    let iv = interval_rows(&c);
    let names: Vec<&str> = iv.iter().map(|r| r.id.as_str()).collect();
    assert!(
        iv.len() > 80 && names[names.len() - 10..] == HAND,
        "{} {:?}",
        iv.len(),
        &names[names.len() - 10..]
    );
    let spans = |l: &str| {
        c.loci[l]
            .targets
            .iter()
            .map(|s| (s.start, s.end))
            .collect::<Vec<_>>()
    };
    assert_eq!(spans("HBB_619bp_del"), [(5_225_388, 5_226_007)]);
    assert_eq!(spans("HBD_HBB_Sicilian_del"), [(5_219_877, 5_233_280)]);
    assert!(
        iv.iter()
            .all(|r| raw[&r.id]["cross_check"] == "cnvpytor:cnvs")
    );
    let q15: Vec<&&Locus> = iv.iter().filter(|r| r.id.starts_with("15q13.3")).collect();
    assert!(!q15.is_empty());
    assert!(
        q15.iter().all(|r| r
            .targets
            .iter()
            .all(|s| !(32_154_858 < s.end && s.start < 32_178_018))),
        "CHRNA7"
    );
    let mut hits = Vec::new();
    for r in &iv {
        for ctl in &r.controls {
            for o in &iv {
                if o.id != r.id
                    && o.chrom == r.chrom
                    && o.targets
                        .iter()
                        .any(|t| t.start < ctl.end && ctl.start < t.end)
                {
                    hits.push((r.id.clone(), o.id.clone()));
                }
            }
        }
    }
    assert_eq!(
        hits,
        Vec::<(String, String)>::new(),
        "a control on another row's target"
    );
}

#[test]
fn interval_rows_call_deletions_duplications_and_follow_the_sex() {
    let c = cat();
    let (a_d, p17) = (
        &c.loci["22q11.2_DGS_proximal_A-D"],
        &c.loci["17p12_HNPP_CMT1A"],
    );
    for (r, target, want, ch) in [
        (a_d, 15.0, 1, "deletion"),
        (p17, 45.0, 3, "duplication"),
        (a_d, 30.0, 2, "no change"),
        (p17, 30.0, 2, "no change"),
    ] {
        let b = assess(r, &depths(r, target, 30.0), None, Sex::Unknown).unwrap();
        assert_eq!(
            (status(&b).as_str(), b.copies, change(&b)),
            ("called", Some(want), ch),
            "{} at {target}x",
            r.id
        );
    }
    let az = &c.loci["AZFc_DAZ"];
    assert_eq!(
        assess(az, &depths(az, 15.0, 15.0), None, Sex::Male)
            .unwrap()
            .copies,
        Some(4)
    );
    let b = assess(az, &depths(az, 7.5, 15.0), None, Sex::Male).unwrap();
    assert_eq!(
        (b.copies, change(&b)),
        (Some(2), "deletion"),
        "AZFc at half its DAZ depth"
    );
    let a = &c.loci["AZFa"];
    let b = assess(a, &depths(a, 15.0, 15.0), None, Sex::Male).unwrap();
    assert_eq!((b.copies, b.expected), (Some(1), Some(1)));
    let b = assess(a, &depths(a, 0.0, 15.0), None, Sex::Male).unwrap();
    assert_eq!((b.copies, change(&b)), (Some(0), "deletion"));
    for loc in ["AZFa", "AZFb", "AZFc_DAZ"] {
        let r = &c.loci[loc];
        let b = assess(r, &depths(r, 15.0, 15.0), None, Sex::Female).unwrap();
        assert_eq!(
            (b.status, b.reason),
            (Status::NotApplicable, Some(Reason::NoY)),
            "{loc} female"
        );
        let b = assess(r, &depths(r, 15.0, 15.0), None, Sex::Unknown).unwrap();
        assert_eq!(
            (b.status, b.reason),
            (Status::NotAssessed, Some(Reason::SexUnknown)),
            "{loc} unknown"
        );
    }
    let xp = &c.loci["Xp22.31"];
    let b = assess(xp, &depths(xp, 15.0, 15.0), None, Sex::Male).unwrap();
    assert_eq!((b.copies, b.expected), (Some(1), Some(1)));
    let b = assess(xp, &depths(xp, 0.0, 15.0), None, Sex::Male).unwrap();
    assert_eq!(b.copies, Some(0), "the STS deletion");
}

/// samtools bedcov over the lines of a BED at 30x, the segments in `half`
/// at 15x and those in `zero` at 0.
fn fakecov(bed: &str, half: &HashSet<String>, zero: &HashSet<String>) -> String {
    bed.lines()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let d = if zero.contains(c[3]) {
                0.0
            } else if half.contains(c[3]) {
                15.0
            } else {
                30.0
            };
            let len: i64 = c[2].parse::<i64>().unwrap() - c[1].parse::<i64>().unwrap();
            format!("{l}\t{}\n", (d * len as f64) as i64)
        })
        .collect()
}

pub const INTERVAL_COLUMNS: [&str; 11] = [
    "locus", "name", "chrom", "start", "end", "status", "copies", "expected", "change", "estimate",
    "reason",
];

#[test]
fn og_cn_intervals_writes_one_line_per_interval_row() {
    let (c, w) = (cat(), Scratch::new("intervals"));
    let iv = interval_rows(&c);
    let names: Vec<&str> = iv.iter().map(|r| r.id.as_str()).collect();
    let cp = catalog_path();
    let (rc, bed, err) = og_cn(&["bed", "--catalog", &cp, "--kind", "interval"]);
    assert_eq!(rc, 0, "{err}");
    assert_eq!(
        bed.lines().count(),
        iv.iter().map(|r| r.regions().len()).sum::<usize>()
    );
    let (a_d, p17) = (
        &c.loci["22q11.2_DGS_proximal_A-D"],
        &c.loci["17p12_HNPP_CMT1A"],
    );
    let half: HashSet<String> = a_d.targets.iter().map(|s| s.name.clone()).collect();
    let zero = HashSet::from([p17.controls[0].name.clone()]);
    let cov = w.put("cov.tsv", &fakecov(&bed, &half, &zero));
    let tab = w.path("iv.tsv");
    let (rc, _, err) = og_cn(&[
        "intervals",
        "--catalog",
        &cp,
        "--bedcov",
        &cov,
        "--sex",
        "female",
        "--out",
        &tab,
    ]);
    assert_eq!(rc, 0, "{err}");
    let text = read(&tab);
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some(INTERVAL_COLUMNS.join("\t").as_str()));
    let rows: IndexMap<&str, Vec<&str>> = lines
        .map(|l| (l.split('\t').next().unwrap(), l.split('\t').collect()))
        .collect();
    assert_eq!(rows.keys().copied().collect::<Vec<_>>(), names);
    assert_eq!(
        rows[a_d.id.as_str()][5..9],
        ["called", "1", "2", "deletion"]
    );
    assert_eq!(
        rows[p17.id.as_str()][10],
        "broken_control",
        "a control with no read"
    );
    assert_eq!(
        rows["3q29"][5..9],
        ["called", "2", "2", "no change"],
        "the rest are called"
    );
    for loc in ["AZFa", "AZFb", "AZFc_DAZ"] {
        assert_eq!(rows[loc][5], "not applicable", "{loc} on a female");
    }
    let (rc, _, err) = og_cn(&["intervals", "--catalog", &cp, "--failed", "--out", &tab]);
    assert_eq!(rc, 0, "{err}");
    let text = read(&tab);
    let failed: Vec<&str> = text.lines().skip(1).collect();
    assert_eq!(failed.len(), iv.len());
    assert!(
        failed
            .iter()
            .all(|l| l.split('\t').nth(10) == Some("check_failed"))
    );
}

// 10. the CYP21A2 row (call psv)

const CYP_SITES: [&str; 6] = [
    "PSV_In2G",
    "PSV_I172N",
    "PSV_V237E",
    "PSV_V281L",
    "PSV_Q318X",
    "PSV_R356W",
];

#[test]
fn the_cyp21a2_row_and_its_psv_sites() {
    let c = cat();
    let r = &c.loci["CYP21A2"];
    let names: Vec<&str> = r.psv.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(names, CYP_SITES);
    assert!(matches!(r.call, CallRule::Psv { .. }));
    for f in PSV_FINDINGS {
        assert!(r.readings.contains_key(f), "a reading for {f}");
    }
    let tg: Vec<&str> = r.targets.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(tg, ["CYP21A2_BODY", "CYP21A1P_BODY"]);
    let (g, p) = (&r.targets[0], &r.targets[1]);
    for s in &r.psv {
        assert!(g.start <= s.gene_pos && s.gene_pos < g.end, "{}", s.name);
        assert!(
            p.start <= s.pseudo_pos && s.pseudo_pos < p.end,
            "{}",
            s.name
        );
        assert!(
            (32_700..32_800).contains(&(s.gene_pos - s.pseudo_pos)),
            "{}",
            s.name
        );
    }
    let mut gb = r.psv[0].gene.clone();
    gb.sort_unstable();
    assert_eq!(
        (gb.as_slice(), r.psv[0].pseudo.as_slice()),
        (&b"AC"[..], &b"G"[..]),
        "In2G c.293-13A/C>G"
    );
    let d = r.depth();
    assert_eq!(d.passes, [".", "-Q 1"]);
    assert_eq!((d.mapq0_min, d.mapq0_max), (Some(0.02), None));
}

#[test]
fn broken_cyp21a2_rows_are_refused() {
    let w = Scratch::new("cyp-broken");
    let base = raw_rows()["CYP21A2"].clone();
    let det = base["details"].clone();
    let readings: String = base["readings"]
        .split('|')
        .filter(|x| !x.starts_with("conversion="))
        .collect::<Vec<_>>()
        .join("|");
    for (desc, col, val, err) in [
        (
            "a PSV site of two bases",
            "details",
            det.replace("PSV_In2G:32039080-32039081", "PSV_In2G:32039080-32039082"),
            "must be one base",
        ),
        (
            "a site with the same base on both genes",
            "details",
            det.replace("32006352 AC>G", "32006352 AC>C"),
            "distinct gene and pseudogene",
        ),
        (
            "a site label without its pseudogene position",
            "details",
            det.replace("32006352 AC>G", "here AC>G"),
            "is not 'PSEUDO_POS",
        ),
        (
            "a psv row without a conversion reading",
            "readings",
            readings,
            "a reading for each",
        ),
    ] {
        let mut r = base.clone();
        assert_ne!(r[col], val, "{desc}: nothing planted");
        r.insert(col.into(), val);
        let got = load_rows(&w, &[r]);
        assert!(
            got.clone().is_err_and(|e| e.contains(err)),
            "{desc}: {got:?}"
        );
    }
}

#[test]
fn bed_psv_prints_twelve_one_base_positions_and_refuses_a_row_without_sites() {
    let cp = catalog_path();
    let (rc, out, err) = og_cn(&["bed", "--catalog", &cp, "--psv", "--locus", "CYP21A2"]);
    assert_eq!(rc, 0, "{err}");
    let starts: Vec<i64> = out
        .lines()
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect();
    assert_eq!(starts.len(), 12);
    assert!(out.lines().all(|l| {
        let c: Vec<i64> = l
            .split('\t')
            .skip(1)
            .take(2)
            .map(|x| x.parse().unwrap())
            .collect();
        c[1] - c[0] == 1
    }));
    let mut sorted = starts.clone();
    sorted.sort_unstable();
    assert_eq!(starts, sorted);
    let (rc, _, err) = og_cn(&["bed", "--catalog", &cp, "--psv", "--locus", "HBA"]);
    assert!(rc == 2 && err.contains("no psv locus"), "{rc} {err}");
}

/// (bedcov text, mpileup text) of a CYP21A2 genome: 15x per copy, flanks at
/// 30x, two CYP21A2 copies on CYP21A2's body and the rest of the module
/// total on CYP21A1P's. At each site the reads carrying CYP21A2's base
/// (`a2` copies, 2 by default) and CYP21A1P's (the rest) are split 4 to 1
/// between the CYP21A2 and the CYP21A1P position, with read starts, ends and
/// an indel in the bases column.
fn cyp_fixture(r: &Locus, total: f64, a2: &[(&str, f64)]) -> (String, String) {
    let mut depth: Means = r.regions().iter().map(|x| (x.3.to_owned(), 30.0)).collect();
    depth.insert("CYP21A2_BODY".into(), 30.0);
    depth.insert("CYP21A1P_BODY".into(), 15.0 * (total - 2.0));
    let bed: String = r
        .regions()
        .iter()
        .map(|(c, s, e, n)| {
            format!(
                "{c}\t{s}\t{e}\t{n}\t{}\n",
                round(depth[*n] * (e - s) as f64)
            )
        })
        .collect();
    let mut lines: Vec<(i64, String)> = Vec::new();
    for site in &r.psv {
        let copies = a2
            .iter()
            .find(|(n, _)| *n == site.name)
            .map_or(2.0, |x| x.1);
        let (gene, pseudo) = (15.0 * copies, 15.0 * (total - copies));
        let gbase = *site.gene.iter().max().unwrap() as char;
        let pbase = *site.pseudo.iter().min().unwrap() as char;
        for (pos, share) in [(site.gene_pos, 0.8), (site.pseudo_pos, 0.2)] {
            let (ng, np) = (round(gene * share), round(pseudo * share));
            let mut reads: Vec<String> = std::iter::repeat_n(gbase, ng as usize)
                .chain(std::iter::repeat_n(pbase, np as usize))
                .enumerate()
                .map(|(i, x)| {
                    if i % 2 == 1 {
                        x.to_ascii_lowercase().to_string()
                    } else {
                        x.to_string()
                    }
                })
                .collect();
            if !reads.is_empty() {
                reads[0] = format!("^]{}", reads[0]);
                let last = reads.len() - 1;
                reads[last] += "$";
            }
            if reads.len() > 2 {
                reads[1] += "+2AG";
            }
            let n = reads.len();
            lines.push((
                pos + 1,
                format!(
                    "{}\t{}\tN\t{n}\t{}\t{}\n",
                    r.chrom,
                    pos + 1,
                    reads.concat(),
                    "I".repeat(n)
                ),
            ));
        }
    }
    lines.sort();
    (bed, lines.into_iter().map(|x| x.1).collect())
}

/// What `og cn call --locus CYP21A2` wrote for a genome, read back.
#[derive(Debug)]
struct Cyp {
    table: Vec<(String, String)>,
    line: String,
    sites: HashMap<String, (f64, f64)>,
}

impl Cyp {
    fn v(&self, k: &str) -> Option<&str> {
        Some(get(&self.table, k)).filter(|v| !v.is_empty() && *v != ".")
    }

    fn copies(&self) -> Option<i64> {
        self.v("cyp21a2_copies").map(|v| v.parse().unwrap())
    }

    fn total(&self) -> Option<i64> {
        self.v("module_total").map(|v| v.parse().unwrap())
    }

    fn conversion_sites(&self) -> Vec<&str> {
        self.v("conversion_sites")
            .map_or(Vec::new(), |v| v.split(',').collect())
    }
}

/// `copies X (EST); ratio RATIO` in a psv_ row: (EST, RATIO).
fn site_numbers(v: &str) -> Option<(f64, f64)> {
    let est = v
        .split_once("copies ")?
        .1
        .split_once(" (")?
        .1
        .split_once(')')?
        .0
        .parse()
        .ok()?;
    let ratio = v
        .split_once("ratio ")?
        .1
        .split(|c: char| !(c.is_ascii_digit() || c == '.'))
        .next()?
        .parse()
        .ok()?;
    Some((est, ratio))
}

/// og cn call on a CYP21A2 genome; `mapq0` is the share of every segment's
/// depth from MAPQ 0 reads (0.1 as HG002's reads give over the two targets).
fn cyp_call(bed: &str, pileup: &str, mapq0: f64) -> Cyp {
    static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let w = Scratch::new(&format!(
        "cyp-{}",
        N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
    ));
    let b = w.put("b.tsv", bed);
    let q = w.put(
        "q.tsv",
        &bed.lines()
            .map(|l| {
                let c: Vec<&str> = l.split('\t').collect();
                format!(
                    "{}\t{:?}\n",
                    c[..4].join("\t"),
                    c[4].parse::<f64>().unwrap() * (1.0 - mapq0)
                )
            })
            .collect::<String>(),
    );
    let p = w.put("p.tsv", pileup);
    let tab = w.path("t.tsv");
    let (_, out, _) = og_cn(&[
        "call",
        "--catalog",
        &catalog_path(),
        "--locus",
        "CYP21A2",
        "--bedcov",
        &b,
        "--bedcov-mapq1",
        &q,
        "--pileup",
        &p,
        "--table",
        &tab,
    ]);
    let table = table(&tab);
    let sites = table
        .iter()
        .filter_map(|(k, v)| Some((k.strip_prefix("psv_")?.to_owned(), site_numbers(v)?)))
        .collect();
    Cyp {
        table,
        line: out.trim().to_owned(),
        sites,
    }
}

fn cyp_row() -> Locus {
    cat().loci["CYP21A2"].clone()
}

/// A CYP21A2 case: its name, the module total and the site readings.
type CypCase = (&'static str, f64, &'static [(&'static str, f64)]);

const CYP_CASES: [CypCase; 3] = [
    ("normal", 4.0, &[]),
    // the 30 kb chimeric deletion: one module fewer, CYP21A1P's sequence at the 5' sites
    ("deletion", 3.0, &[("PSV_In2G", 1.0), ("PSV_I172N", 1.0)]),
    // In2G converted to the pseudogene base on one of the two CYP21A2 copies
    ("conversion", 4.0, &[("PSV_In2G", 1.0)]),
];

#[test]
fn the_committed_cyp21a2_fixtures_are_the_generators_and_call_right() {
    let r = cyp_row();
    let mut got = HashMap::new();
    for (case, total, a2) in CYP_CASES {
        let (bed, pile) = cyp_fixture(&r, total, a2);
        let fb = read(&repo(&format!(
            "tests/fixtures/cn/CYP21A2_{case}.bedcov.tsv"
        )));
        let fp = read(&repo(&format!(
            "tests/fixtures/cn/CYP21A2_{case}.pileup.tsv"
        )));
        assert!(
            fb == bed && fp == pile,
            "tests/fixtures/cn/CYP21A2_{case} differs from the generator"
        );
        got.insert(case, cyp_call(&bed, &pile, 0.1));
    }
    let b = &got["normal"];
    assert_eq!(
        (b.v("status"), b.copies(), b.total(), b.v("finding")),
        (Some("called"), Some(2), Some(4), Some("none"))
    );
    assert_eq!(b.v("reading").unwrap_or_default(), r.reading_word("none"));
    let b = &got["deletion"];
    assert_eq!(
        (b.v("status"), b.copies(), b.total(), b.v("finding")),
        (Some("called"), Some(1), Some(3), Some("deletion"))
    );
    let b = &got["conversion"];
    let ratio = b.sites.get("PSV_In2G").map(|x| x.1).unwrap_or(f64::NAN);
    assert!((ratio - 0.5).abs() <= 0.1, "{b:?}");
    assert_eq!(
        (b.v("finding"), b.conversion_sites(), b.total()),
        (Some("conversion"), vec!["PSV_In2G"], Some(4))
    );
    let reading = b.v("reading").unwrap_or_default();
    assert!(
        reading.contains("c.293-13A/C>G (In2G") && reading.contains("possible conversion"),
        "{reading}"
    );
}

#[test]
fn the_cyp21a2_call_rests_on_the_psv_reads() {
    let r = cyp_row();
    let (bed, _) = cyp_fixture(&r, 4.0, &[("PSV_In2G", 1.0)]);
    let (_, even) = cyp_fixture(&r, 4.0, &[]);
    assert_eq!(
        cyp_call(&bed, &even, 0.1).v("finding"),
        Some("none"),
        "the conversion's depths, every site even"
    );
    let (bed, pile) = cyp_fixture(&r, 3.0, &[]);
    let b = cyp_call(&bed, &pile, 0.1);
    assert_eq!(
        (b.copies(), b.v("finding")),
        (Some(2), Some("pseudogene_deletion"))
    );
    // One noisy 3' site (V281L at 1.47, rounds to 1) inside the 0.6 margin.
    for (total, want) in [(3.0, "pseudogene_deletion"), (4.0, "none")] {
        let (bed, pile) = cyp_fixture(&r, total, &[("PSV_V281L", 1.47)]);
        let b = cyp_call(&bed, &pile, 0.1);
        let est = b.sites.get("PSV_V281L").map_or(f64::NAN, |x| x.0);
        assert!(1.4 < est && est < 1.5, "{b:?}");
        assert_eq!(
            (b.copies(), b.v("finding")),
            (Some(2), Some(want)),
            "module total {total}"
        );
    }
    let all: Vec<(&str, f64)> = CYP_SITES.iter().map(|n| (*n, 1.0)).collect();
    let (bed, pile) = cyp_fixture(&r, 4.0, &all);
    let b = cyp_call(&bed, &pile, 0.1);
    assert_eq!(
        (b.v("finding"), b.conversion_sites(), b.copies()),
        (Some("conversion"), CYP_SITES.to_vec(), Some(1))
    );
    let (bed, _) = cyp_fixture(&r, 4.0, &[]);
    let thin: String = even
        .lines()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            format!("{}\t{}\tN\t1\tC\tI\n", c[0], c[1])
        })
        .collect();
    let b = cyp_call(&bed, &thin, 0.1);
    assert_eq!(
        (b.v("status"), b.v("reason")),
        (Some("not assessed"), Some("psv_low_depth"))
    );
    let (bed, pile) = cyp_fixture(&r, 4.0, &[]);
    let from = format!("\t{}\n", round(30.0 * 2689.0));
    assert!(bed.contains(&from));
    let b = cyp_call(
        &bed.replace(&from, &format!("\t{}\n", round(37.5 * 2689.0))),
        &pile,
        0.1,
    );
    assert_eq!(
        (b.v("status"), b.v("reason"), b.v("finding")),
        (Some("not assessed"), Some("between_counts"), None)
    );
}

#[test]
fn the_cyp21a2_mapq0_guard() {
    let r = cyp_row();
    let (bed, pile) = cyp_fixture(&r, 3.0, &[("PSV_In2G", 1.0), ("PSV_I172N", 1.0)]);
    for (share, want) in [
        (
            0.0,
            (Some("not assessed"), Some("low_mapq0_share"), None, None),
        ),
        (
            0.019,
            (Some("not assessed"), Some("low_mapq0_share"), None, None),
        ),
        (0.021, (Some("called"), None, Some("deletion"), Some(1))),
        (0.1, (Some("called"), None, Some("deletion"), Some(1))),
    ] {
        let b = cyp_call(&bed, &pile, share);
        assert_eq!(
            (b.v("status"), b.v("reason"), b.v("finding"), b.copies()),
            want,
            "share {share}"
        );
    }
    let b = cyp_call(&bed, &pile, 0.0);
    let msg = get(&b.table, "message");
    assert!(
        msg.starts_with("not assessed: the BAM has few or no multi-mapping reads here"),
        "{msg}"
    );
    assert!(msg.contains("CYP21A2 and CYP21A1P"), "{msg}");
    assert_eq!(get(&b.table, "mapq0_share"), "0.000");
    assert_eq!(
        (
            get(&b.table, "cyp21a2_copies"),
            get(&b.table, "module_total"),
            get(&b.table, "finding")
        ),
        (".", ".", ".")
    );
    assert!(b.line.contains("low_mapq0_share"), "{}", b.line);
    let (bed, pile) = cyp_fixture(&r, 4.0, &[]);
    let t = cyp_call(&bed, &pile, 0.1).table;
    assert_eq!(
        (
            get(&t, "status"),
            get(&t, "finding"),
            get(&t, "mapq0_share"),
            get(&t, "message")
        ),
        ("called", "none", "0.100", ".")
    );
}
