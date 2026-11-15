//! `og run delta` on two synthetic runs of the same sample. It replaces
//! `tests/test_run_delta.py`.
//!
//! 1. no previous run: the delta says so, lists no change, and still
//!    writes the text, HTML and JSON files;
//! 2. two identical runs (only the time each was written differs): no
//!    change;
//! 3. a planted ClinVar reclassification and a changed image digest: the
//!    delta names exactly those two changes, the ClinVar one with the
//!    digest change as a possible cause only when that image feeds ClinVar,
//!    and with "the call itself differs" when it does not;
//! 4. a hit added, one dropped, a re-starred one, a CPIC call, pypgx's
//!    CYP2D6, a PRS percentile that moved past the band and one that stayed
//!    inside it, a CPSR count, a section that stopped running, a new
//!    ClinVar file date and a new versions.env tag: each is listed once,
//!    and the reclassified hit names the ClinVar file date as the reason; a
//!    new section that did not run is no change, one that ran is;
//! 5. the HTML escapes what it prints;
//! 6. with a previous run the delta says "New since <the previous run's
//!    date>", names the sections whose results changed, and adds the last
//!    database re-check (recheck/state.json beside the summary); with no
//!    previous run there is no previous_utc, and the re-check line is still
//!    printed in the text and the HTML;
//! 7. a re-check older than the run's manifest is left out, and a rerun
//!    block that is not an object is dropped;
//! 8. `--newest-history SAMPLE_DIR` compares with the newest
//!    history/<stamp>/ that holds a summary.json, passing over a newer one
//!    without it; a sample with no history/ is a first run; `--previous`
//!    and `--newest-history` together exit 2.

use og_runlog::cli::{Env, run};
use pycompat::json::{Dump, Value, dumps, loads};
use std::fs;
use std::path::{Path, PathBuf};

const SUMMARY: &str = r#"{"schema_version": 1, "sample": "S", "generated_utc": "GENERATED",
 "sections": {
  "clinvar": {"title": "ClinVar screen", "step": "x", "state": "ok", "source": null, "file_date": null, "note": null,
   "data": {"count": 2, "hits": [
    {"stars": 2, "chrom": "chr1", "pos": 100, "ref": "A", "alt": "G", "genotype": "het", "gene": "GENEA",
     "significance": "Pathogenic", "review_status": "criteria provided, multiple submitters, no conflicts"},
    {"stars": 2, "chrom": "chr1", "pos": 200, "ref": "A", "alt": "G", "genotype": "het", "gene": "GENEA",
     "significance": "Likely pathogenic", "review_status": "criteria provided, multiple submitters, no conflicts"}]}},
  "cpsr": {"title": "Cancer predisposition (CPSR)", "step": "x", "state": "ok", "source": null, "file_date": null,
   "note": null, "data": {"classification": {"Benign": 10, "VUS": 3}}},
  "cpic": {"title": "CPIC drug-gene recommendations", "step": "x", "state": "ok", "source": null, "file_date": null,
   "note": null, "data": {"genes": [
    {"gene": "CYP2C19", "diplotype": "*1/*1", "phenotype": "Normal Metabolizer", "status": "normal"}]}},
  "pypgx": {"title": "pypgx", "step": "x", "state": "ok", "source": null, "file_date": null, "note": null,
   "data": {"cyp2d6": "*1/*2"}},
  "prs": {"title": "Polygenic risk scores", "step": "x", "state": "ok", "source": null, "file_date": null,
   "note": null, "data": {"scores": [
    {"pgs_id": "PGS000001", "condition": "Trait one", "percentile": "40.0"},
    {"pgs_id": "PGS000002", "condition": "Trait two", "percentile": "50.0"}]}},
  "roh": {"title": "Runs of homozygosity", "step": "x", "state": "ok", "source": null, "file_date": null,
   "note": null, "data": {}}},
 "cyp2d6": {"consensus": null}}"#;

const MANIFEST: [(&str, &str, &str); 12] = [
    ("run", "written_utc", "2026-01-01 00:00:00 UTC"),
    ("run", "written_by", "run-all.sh"),
    ("run", "git_commit", "abc"),
    ("version", "PYTHON_IMAGE", "python-ref"),
    ("version", "DEEPVARIANT_IMAGE", "deepvariant-ref"),
    ("version", "PHARMCAT_IMAGE", "pharmcat-ref"),
    ("version", "VEP_CACHE_RELEASE", "115"),
    ("image", "python-ref", "digest-aaa"),
    ("image", "deepvariant-ref", "digest-bbb"),
    ("image", "pharmcat-ref", "digest-ccc"),
    (
        "data",
        "clinvar_file_date",
        "20260101 (clinvar/clinvar.vcf.gz)",
    ),
    (
        "data",
        "ancestry_panel",
        "not installed (setup.sh --ancestry-panel)",
    ),
];

type Manifest = Vec<(String, String, String)>;

fn manifest() -> Manifest {
    MANIFEST
        .iter()
        .map(|(a, b, c)| ((*a).to_owned(), (*b).to_owned(), (*c).to_owned()))
        .collect()
}

fn set_row(m: &Manifest, sec: &str, key: &str, value: &str) -> Manifest {
    m.iter()
        .map(|(s, k, v)| {
            let v = if s == sec && k == key { value } else { v };
            (s.clone(), k.clone(), v.to_owned())
        })
        .collect()
}

fn summary(generated: &str) -> Value {
    loads(&SUMMARY.replace("GENERATED", generated)).unwrap()
}

fn s(x: &str) -> Value {
    Value::Str(x.to_owned())
}

fn json(text: &str) -> Value {
    loads(text).unwrap()
}

/// The value at a dotted path (object keys, array indexes), for editing.
fn at<'a>(v: &'a mut Value, path: &str) -> &'a mut Value {
    path.split('.').fold(v, |v, k| match v {
        Value::Object(o) => o.entry(k.to_owned()).or_insert(Value::Null),
        Value::Array(a) => &mut a[k.parse::<usize>().unwrap()],
        _ => panic!("{path}: {k} is not inside an object or an array"),
    })
}

fn hits(v: &mut Value) -> &mut Vec<Value> {
    match at(v, "sections.clinvar.data.hits") {
        Value::Array(a) => a,
        _ => unreachable!(),
    }
}

fn hit(pos: i128, sig: &str, gene: &str) -> Value {
    let mut h = summary("x");
    let mut one = hits(&mut h)[0].clone();
    *at(&mut one, "pos") = Value::Int(pos);
    *at(&mut one, "significance") = s(sig);
    *at(&mut one, "gene") = s(gene);
    one
}

fn get<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(o) => o.get(k),
        _ => None,
    }
}

fn text(v: Option<&Value>) -> &str {
    match v {
        Some(Value::Str(x)) => x,
        _ => "",
    }
}

fn write_run(d: &Path, summ: &Value, m: &Manifest) {
    fs::create_dir_all(d).unwrap();
    fs::write(
        d.join("summary.json"),
        dumps(summ, &Dump::default()).unwrap(),
    )
    .unwrap();
    let rows: String = m
        .iter()
        .map(|(a, b, c)| format!("{a}\t{b}\t{c}\n"))
        .collect();
    fs::write(d.join("run_manifest.tsv"), rows).unwrap();
}

fn og_run(args: &[String]) -> (u8, String) {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = run(args, &Env::system(), &mut out, &mut err);
    (code, String::from_utf8(err).unwrap())
}

struct W(PathBuf);

impl W {
    fn args(&self, prev: Option<&Path>, cur: &Path, name: &str, extra: &[&str]) -> Vec<String> {
        let mut a: Vec<String> = [
            "delta",
            "--sample",
            "S",
            "--summary",
            cur.join("summary.json").to_str().unwrap(),
            "--manifest",
            cur.join("run_manifest.tsv").to_str().unwrap(),
            "--out",
            self.0.join(name).to_str().unwrap(),
        ]
        .map(str::to_owned)
        .to_vec();
        a.extend(extra.iter().map(|x| (*x).to_owned()));
        if let Some(p) = prev {
            a.extend(["--previous".to_owned(), p.display().to_string()]);
        }
        a
    }

    /// The written (JSON, text, HTML) of one delta, which must exit 0.
    fn delta(
        &self,
        prev: Option<&Path>,
        cur: &Path,
        name: &str,
        extra: &[&str],
    ) -> (Value, String, String) {
        let (code, err) = og_run(&self.args(prev, cur, name, extra));
        assert_eq!(code, 0, "{name}: exits 0: {err}");
        let read = |ext: &str| fs::read_to_string(self.0.join(format!("{name}.{ext}"))).unwrap();
        (json(&read("json")), read("txt"), read("html"))
    }
}

fn changes(d: &Value) -> Vec<Value> {
    match get(d, "changes") {
        Some(Value::Array(a)) => a.clone(),
        other => panic!("changes is not a list: {other:?}"),
    }
}

/// (area, kind, item) of each change, sorted.
fn kinds(d: &Value) -> Vec<(String, String, String)> {
    let mut v: Vec<_> = changes(d)
        .iter()
        .map(|c| {
            (
                text(get(c, "area")).to_owned(),
                text(get(c, "kind")).to_owned(),
                text(get(c, "item")).to_owned(),
            )
        })
        .collect();
    v.sort();
    v
}

fn find(d: &Value, key: &str, value: &str) -> Value {
    changes(d)
        .into_iter()
        .find(|c| text(get(c, key)) == value)
        .unwrap_or(Value::Null)
}

fn triples(v: &[(&str, &str, &str)]) -> Vec<(String, String, String)> {
    let mut v: Vec<_> = v
        .iter()
        .map(|(a, b, c)| ((*a).to_owned(), (*b).to_owned(), (*c).to_owned()))
        .collect();
    v.sort();
    v
}

const FEB: &str = "2026-02-01T00:00:00Z";

#[test]
fn deltas_of_two_runs() {
    let root = std::env::temp_dir().join(format!("og-run-delta-{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let w = W(root.clone());
    let prev = root.join("history/2026-01-01T000000Z");
    let cur = root.join("cur");
    write_run(&prev, &summary("2026-01-01T00:00:00Z"), &manifest());

    // 1. no previous run
    write_run(
        &cur,
        &summary(FEB),
        &set_row(&manifest(), "run", "written_utc", "later"),
    );
    let (d, txt, page) = w.delta(None, &cur, "first", &[]);
    assert!(
        get(&d, "previous") == Some(&Value::Null) && changes(&d).is_empty(),
        "first run: previous is null and no change: {d:?}"
    );
    assert!(
        txt.contains("No previous run of this sample"),
        "first run: the text says so"
    );
    assert!(
        page.contains("No previous run of this sample"),
        "first run: the HTML says so"
    );

    // 2. identical runs
    let (d, txt, _) = w.delta(Some(&prev), &cur, "same", &[]);
    assert!(changes(&d).is_empty(), "identical runs: no change: {d:?}");
    assert!(
        txt.contains("No change:"),
        "identical runs: the text says no change"
    );
    assert_eq!(
        text(get(&d, "previous")),
        "2026-01-01T000000Z",
        "previous is the history folder name"
    );

    // 3. one reclassification and one new digest
    let mut sm = summary(FEB);
    *at(&mut sm, "sections.clinvar.data.hits.1.significance") = s("Benign");
    write_run(
        &cur,
        &sm,
        &set_row(&manifest(), "image", "deepvariant-ref", "digest-ddd"),
    );
    let (d, txt, _) = w.delta(Some(&prev), &cur, "planted", &[]);
    assert_eq!(
        kinds(&d),
        triples(&[
            ("ClinVar", "reclassified", "chr1:200 A>G (GENEA)"),
            ("inputs", "input", "DEEPVARIANT_IMAGE")
        ]),
        "planted: exactly the two planted changes"
    );
    let cv = find(&d, "kind", "reclassified");
    assert_eq!(
        (text(get(&cv, "before")), text(get(&cv, "after"))),
        ("Likely pathogenic", "Benign"),
        "planted: the reclassification reads Likely pathogenic -> Benign"
    );
    let why = text(get(&cv, "why"));
    assert!(
        why.contains("DEEPVARIANT_IMAGE") && why.contains("digest-ddd"),
        "planted: the new DeepVariant digest is named as a possible cause: {why}"
    );
    let img = find(&d, "item", "DEEPVARIANT_IMAGE");
    assert!(
        text(get(&img, "before")).contains("digest-bbb")
            && text(get(&img, "after")).contains("digest-ddd"),
        "planted: the digest change shows both digests: {img:?}"
    );
    assert!(
        txt.contains("2 change(s):"),
        "planted: the text lists 2 changes: {txt}"
    );

    // 3b. the same reclassification with only the PharmCAT digest changed
    write_run(
        &cur,
        &sm,
        &set_row(&manifest(), "image", "pharmcat-ref", "digest-eee"),
    );
    let (d, _, _) = w.delta(Some(&prev), &cur, "unrelated", &[]);
    assert_eq!(
        text(get(&find(&d, "kind", "reclassified"), "why")),
        "the inputs this delta tracks for it are the same; the call itself differs",
        "a PharmCAT digest does not explain a ClinVar change"
    );

    // 4. one of each
    let mut sm = summary(FEB);
    {
        let h = hits(&mut sm);
        *at(&mut h[0], "stars") = Value::Int(3);
        *at(&mut h[0], "review_status") = s("reviewed by expert panel");
        *at(&mut h[1], "significance") = s("Pathogenic");
        h.push(hit(300, "Pathogenic", "GENEB"));
        let first = h[0].clone();
        h.retain(|x| get(x, "pos") != Some(&Value::Int(100)));
        h.push(hit(400, "Pathogenic", "GENEA"));
        let mut again = first;
        *at(&mut again, "pos") = Value::Int(500);
        h.push(again);
    }
    *at(&mut sm, "sections.cpic.data.genes.0.diplotype") = s("*1/*2");
    *at(&mut sm, "sections.cpic.data.genes.0.phenotype") = s("Intermediate Metabolizer");
    *at(&mut sm, "sections.pypgx.data.cyp2d6") = s("*1/*4");
    *at(&mut sm, "sections.prs.data.scores.0.percentile") = s("52.0");
    *at(&mut sm, "sections.prs.data.scores.1.percentile") = s("53.0");
    *at(&mut sm, "sections.cpsr.data.classification") = json(r#"{"Benign": 10, "VUS": 4}"#);
    *at(&mut sm, "sections.roh.state") = s("missing");
    let m = set_row(
        &manifest(),
        "data",
        "clinvar_file_date",
        "20260201 (clinvar/clinvar.vcf.gz)",
    );
    write_run(
        &cur,
        &sm,
        &set_row(&m, "version", "VEP_CACHE_RELEASE", "116"),
    );
    let (d, _, _) = w.delta(Some(&prev), &cur, "many", &[]);
    assert_eq!(
        kinds(&d),
        triples(&[
            ("inputs", "input", "VEP_CACHE_RELEASE"),
            ("inputs", "input", "clinvar_file_date"),
            ("Runs of homozygosity", "section", "roh"),
            ("ClinVar", "dropped", "chr1:100 A>G (GENEA)"),
            ("ClinVar", "reclassified", "chr1:200 A>G (GENEA)"),
            ("ClinVar", "added", "chr1:300 A>G (GENEB)"),
            ("ClinVar", "added", "chr1:400 A>G (GENEA)"),
            ("ClinVar", "added", "chr1:500 A>G (GENEA)"),
            ("CPSR", "classification count", "VUS"),
            ("PGx (CPIC)", "call", "CYP2C19"),
            ("PGx (pypgx)", "call", "CYP2D6"),
            ("PRS", "percentile", "PGS000001 (Trait one)"),
        ]),
        "each change once, the 3-point PRS move left out"
    );
    let why = text(get(&find(&d, "kind", "reclassified"), "why")).to_owned();
    assert!(
        why.contains("clinvar_file_date 20260101") && why.contains("20260201"),
        "the reclassification names the new ClinVar file date: {why}"
    );
    let why = text(get(&find(&d, "area", "CPSR"), "why")).to_owned();
    assert!(
        why.contains("VEP_CACHE_RELEASE 115 -> 116"),
        "the CPSR count names the new VEP cache release: {why}"
    );

    // 4b. re-starred, on its own
    let mut sm = summary(FEB);
    *at(&mut sm, "sections.clinvar.data.hits.0.stars") = Value::Int(3);
    *at(&mut sm, "sections.clinvar.data.hits.0.review_status") = s("reviewed by expert panel");
    write_run(&cur, &sm, &manifest());
    let (d, _, _) = w.delta(Some(&prev), &cur, "restarred", &[]);
    let got: Vec<[String; 3]> = changes(&d)
        .iter()
        .map(|c| ["kind", "before", "after"].map(|k| text(get(c, k)).to_owned()))
        .collect();
    assert_eq!(
        got,
        [[
            "re-starred",
            "2 stars (criteria provided, multiple submitters, no conflicts)",
            "3 stars (reviewed by expert panel)"
        ]
        .map(str::to_owned)],
        "a new review status is one re-starred change"
    );

    // 4c. a section only the newer summary knows (a pipeline upgrade added it)
    let mut sm = summary(FEB);
    let new_section = |state: &str| {
        json(&format!(
            r#"{{"title": "New step", "step": "x", "state": "{state}", "source": null, "file_date": null, "note": null, "data": {{}}}}"#
        ))
    };
    *at(&mut sm, "sections.newstep") = new_section("missing");
    write_run(&cur, &sm, &manifest());
    let (d, _, _) = w.delta(Some(&prev), &cur, "newsection", &[]);
    assert!(
        changes(&d).is_empty(),
        "a new section that did not run is no change: {d:?}"
    );
    *at(&mut sm, "sections.newstep") = new_section("ok");
    write_run(&cur, &sm, &manifest());
    let (d, _, _) = w.delta(Some(&prev), &cur, "newsection-ran", &[]);
    let got: Vec<[String; 4]> = changes(&d)
        .iter()
        .map(|c| ["kind", "item", "before", "after"].map(|k| text(get(c, k)).to_owned()))
        .collect();
    assert_eq!(
        got,
        [["section", "newstep", "missing", "ok"].map(str::to_owned)],
        "a new section that ran is one section change"
    );

    // 5. HTML escaping
    let mut sm = summary(FEB);
    *at(&mut sm, "sections.clinvar.data.hits.1.significance") = s("<b>x</b>");
    write_run(&cur, &sm, &manifest());
    let (_, _, page) = w.delta(Some(&prev), &cur, "escape", &[]);
    assert!(
        page.contains("&lt;b&gt;x&lt;/b&gt;") && !page.contains("<b>x</b>"),
        "the HTML escapes the values"
    );

    // 6. new since, the changed sections and the database re-check
    let mut sm = summary(FEB);
    *at(&mut sm, "sections.clinvar.data.hits.1.significance") = s("Benign");
    *at(&mut sm, "sections.roh.state") = s("missing");
    write_run(&cur, &sm, &manifest());
    let state = r#"{"checked_utc": "2026-02-02T00:00:00Z", "affected_steps": ["06", "24"], "items": [
        {"item": "clinvar_file_date", "label": "ClinVar", "class": "refresh-available", "action": "refresh it"},
        {"item": "VEP_CACHE_RELEASE", "label": "VEP", "class": "no-newer-release"}]}"#;
    fs::create_dir_all(cur.join("recheck")).unwrap();
    fs::write(cur.join("recheck/state.json"), state).unwrap();
    let (d, txt, _) = w.delta(Some(&prev), &cur, "newsince", &[]);
    assert_eq!(
        text(get(&d, "previous_utc")),
        "2026-01-01T00:00:00Z",
        "new since: previous_utc is the previous run's date"
    );
    assert!(
        txt.contains("New since 2026-01-01T00:00:00Z: 2 change(s):"),
        "new since: the text heads the changes with it: {txt}"
    );
    assert_eq!(
        get(&d, "new_since_sections"),
        Some(&json(r#"["roh", "clinvar"]"#)),
        "new since: the changed sections are roh and clinvar"
    );
    let steps = get(&d, "recheck").and_then(|r| get(r, "affected_steps"));
    assert!(
        steps == Some(&json(r#"["06", "24"]"#))
            && txt.contains("Database re-check of 2026-02-02T00:00:00Z: 1 item(s) to look at")
            && txt.contains("ClinVar: refresh-available"),
        "new since: the re-check is in the delta: {txt}"
    );
    let (d, txt, _) = w.delta(None, &cur, "newsince-first", &[]);
    assert!(
        get(&d, "previous_utc").is_none()
            && get(&d, "new_since_sections").is_none()
            && !txt.contains("New since"),
        "new since: no previous run, no previous_utc and no changed sections: {d:?}"
    );
    let (_, txt, page) = w.delta(None, &cur, "recheck-first", &[]);
    assert!(
        txt.contains("Database re-check of 2026-02-02T00:00:00Z"),
        "first run: the text shows the re-check"
    );
    assert!(
        page.contains("Database re-check of 2026-02-02T00:00:00Z"),
        "first run: the HTML shows the re-check"
    );

    // 7. a check older than the run's manifest is left out; a rerun block
    // that is not an object is dropped and the delta is still written.
    let late = root.join("late");
    let odd = root.join("odd");
    write_run(
        &late,
        &sm,
        &set_row(&manifest(), "run", "written_utc", "2026-03-01 00:00:00 UTC"),
    );
    write_run(&odd, &sm, &manifest());
    for (dir, rerun) in [(&late, Value::Null), (&odd, json("[1]"))] {
        let mut st = json(state);
        *at(&mut st, "rerun") = rerun;
        fs::create_dir_all(dir.join("recheck")).unwrap();
        fs::write(
            dir.join("recheck/state.json"),
            dumps(&st, &Dump::default()).unwrap(),
        )
        .unwrap();
    }
    let quiet = |dir: &Path, name: &str| -> (u8, Value, String) {
        let (code, _) = og_run(&w.args(Some(&prev), dir, name, &[]));
        let d = fs::read_to_string(root.join(format!("{name}.json")))
            .ok()
            .and_then(|t| loads(&t).ok())
            .unwrap_or(Value::Null);
        let t = fs::read_to_string(root.join(format!("{name}.txt"))).unwrap_or_default();
        (code, d, t)
    };
    let (code, d, t) = quiet(&late, "rc-late");
    assert!(
        code == 0 && get(&d, "recheck").is_none() && !t.contains("Database re-check"),
        "re-check: one older than the run's manifest is left out: {code} {d:?}"
    );
    let (code, d, _) = quiet(&odd, "rc-odd");
    let r = get(&d, "recheck");
    assert!(
        code == 0
            && r.and_then(|r| get(r, "affected_steps")) == Some(&json(r#"["06", "24"]"#))
            && matches!(r.and_then(|r| get(r, "rerun")), None | Some(Value::Null)),
        "re-check: a rerun block that is not an object is dropped, the delta still written: {code} {r:?}"
    );

    // 8. --newest-history: og picks the previous run, as generate-report.sh asks it to.
    let sd = root.join("sample");
    let h = sd.join("history");
    write_run(
        &h.join("2026-01-01T000000Z"),
        &summary("2026-01-01T00:00:00Z"),
        &manifest(),
    );
    write_run(
        &h.join("2026-01-01T000000Z-02"),
        &summary("2026-01-15T00:00:00Z"),
        &manifest(),
    );
    fs::create_dir_all(h.join("2026-03-01T000000Z")).unwrap(); // stopped before its summary.json
    write_run(&sd, &summary(FEB), &manifest());
    let sds = sd.to_str().unwrap();
    let (d, txt, _) = w.delta(None, &sd, "newest", &["--newest-history", sds]);
    assert_eq!(
        (text(get(&d, "previous")), text(get(&d, "previous_utc"))),
        ("2026-01-01T000000Z-02", "2026-01-15T00:00:00Z"),
        "--newest-history: the newest folder with a summary.json, by name"
    );
    assert!(
        txt.contains("Previous run: history/2026-01-01T000000Z-02 "),
        "--newest-history: the text names that folder: {txt}"
    );
    fs::remove_dir_all(&h).unwrap();
    let (d, txt, _) = w.delta(None, &sd, "newest-none", &["--newest-history", sds]);
    assert!(
        get(&d, "previous") == Some(&Value::Null)
            && txt.contains("No previous run of this sample on record"),
        "--newest-history without history/: a first run: {d:?}"
    );
    let both: Vec<String> = [
        "delta",
        "--sample",
        "S",
        "--summary",
        sd.join("summary.json").to_str().unwrap(),
        "--out",
        root.join("both").to_str().unwrap(),
        "--previous",
        prev.to_str().unwrap(),
        "--newest-history",
        sds,
    ]
    .map(str::to_owned)
    .to_vec();
    let (code, err) = og_run(&both);
    assert!(
        code == 2 && err.contains("not allowed with argument --previous"),
        "--previous and --newest-history together exit 2: {code} {err}"
    );
    fs::remove_dir_all(&root).unwrap();
}
