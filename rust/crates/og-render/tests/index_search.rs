//! `og report index` on the synthetic sample of
//! `tests/fixtures/results_index/SYNTH` (invented): ClinVar hits in CFTR and
//! BRCA2, a CPIC table (CYP2C19, CYP2D6), and the tables of steps 23 and 31
//! in both forms, the summary TSVs the bash steps write and the VCFs the
//! Nextflow report stages. The summary is collected by `og summary
//! collect`'s own code.
//!
//! Checks: the terms and their items; that `og lookup` finds every term the
//! index finds in summary.json (one names rule); a prose gene field and a
//! `locus` key; the TSV and the VCF form give the same index; stale, failed
//! and unreadable sections; a table that cannot be read; the query rule; the
//! command line; the veil; and search.js's `match()`, run under node,
//! answering each query as [`query`] does (node is required on CI and the
//! check is skipped with a note without it).

use std::io::Write;
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::Command;

use og_render::index::{Closable, Kind, build, cli, query, query_lines};
use og_summary::collect::{Group, Input};
use pycompat::json::{Dump, Value, dumps, loads};

/// The repository's root.
fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// A scratch folder, removed when dropped.
struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let p = std::env::temp_dir().join(format!(
            "og-render-index-{tag}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = std::fs::remove_dir_all(&p);
        std::fs::create_dir_all(&p).unwrap();
        Dir(p)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).display().to_string()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn copy_tree(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for e in std::fs::read_dir(from).unwrap() {
        let e = e.unwrap();
        let dest = to.join(e.file_name());
        if e.file_type().unwrap().is_dir() {
            copy_tree(&e.path(), &dest);
        } else {
            std::fs::copy(e.path(), dest).unwrap();
        }
    }
}

fn gzip(path: &Path) {
    let body = std::fs::read(path).unwrap();
    let mut gz = PathBuf::from(path).into_os_string();
    gz.push(".gz");
    let f = std::fs::File::create(&gz).unwrap();
    let mut w = og_io::bgzf::Writer::new(f, NonZero::new(1).unwrap()).unwrap();
    w.write_all(&body).unwrap();
    w.finish().unwrap();
    std::fs::remove_file(path).unwrap();
}

/// A copy of the fixture under `work/form/SYNTH` with the step 23 and 31
/// tables in one form: "tsv" (bash: the summary TSVs, the VCFs too) or
/// "vcf" (Nextflow: the VCFs only).
fn sample(work: &Dir, form: &str) -> String {
    let d = work.0.join(form).join("SYNTH");
    copy_tree(&repo().join("tests/fixtures/results_index/SYNTH"), &d);
    gzip(&d.join("clinical/SYNTH_clinical.vcf"));
    gzip(&d.join("slivar/SYNTH_prioritized.vcf"));
    if form == "vcf" {
        std::fs::remove_file(d.join("clinical/SYNTH_clinical_summary.tsv")).unwrap();
        std::fs::remove_file(d.join("slivar/SYNTH_slivar_summary.tsv")).unwrap();
    }
    d.display().to_string()
}

/// summary.json of a sample folder, as `og summary collect` writes it.
fn collect(dir: &str) -> Value {
    let input = Input {
        sample: "SYNTH".to_owned(),
        sample_dir: dir.to_owned(),
        declared_sex: None,
        assets: repo().join("assets").display().to_string(),
        group: Group::All,
        tables: repo().join("bin").display().to_string(),
    };
    og_summary::collect::post::collect(&input).expect("the summary")
}

fn index(summary: &Value, dir: Option<&str>) -> Value {
    build(summary, dir, &[], &mut std::io::sink()).expect("the index")
}

fn get<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(o) => o.get(*k).unwrap_or(&Value::Null),
        _ => &Value::Null,
    })
}

fn terms(idx: &Value) -> Vec<String> {
    match get(idx, &["terms"]) {
        Value::Object(o) => o.keys().cloned().collect(),
        _ => Vec::new(),
    }
}

fn items<'a>(idx: &'a Value, term: &str) -> &'a [Value] {
    match get(idx, &["terms", term]) {
        Value::Array(a) => a,
        _ => &[],
    }
}

fn s(v: &Value) -> String {
    match v {
        Value::Str(x) => x.clone(),
        Value::Int(n) => n.to_string(),
        Value::Null => "null".to_owned(),
        other => format!("{other:?}"),
    }
}

/// (section, anchor, count) of each item of `term`.
fn sections(idx: &Value, term: &str) -> Vec<(String, String, String)> {
    items(idx, term)
        .iter()
        .map(|it| {
            (
                s(get(it, &["section"])),
                s(get(it, &["anchor"])),
                s(get(it, &["count"])),
            )
        })
        .collect()
}

/// {term: [mention text]} of the items of section `key`.
fn table_entries(idx: &Value, key: &str) -> Vec<(String, Vec<String>)> {
    terms(idx)
        .into_iter()
        .filter_map(|t| {
            let texts: Vec<String> = items(idx, &t)
                .iter()
                .filter(|it| s(get(it, &["section"])) == key)
                .flat_map(|it| match get(it, &["mentions"]) {
                    Value::Array(m) => m.iter().map(|x| s(get(x, &["text"]))).collect(),
                    _ => Vec::new(),
                })
                .collect();
            (!texts.is_empty()).then_some((t, texts))
        })
        .collect()
}

fn set(v: &mut Value, path: &[&str], to: Value) {
    let mut cur = v;
    for k in &path[..path.len() - 1] {
        let Value::Object(o) = cur else { panic!("{k}") };
        cur = o.get_mut(*k).expect("the key");
    }
    let Value::Object(o) = cur else {
        panic!("not an object")
    };
    o.insert(path[path.len() - 1].to_owned(), to);
}

fn json(t: &str) -> Value {
    loads(t).unwrap()
}

fn write_json(path: &str, v: &Value) {
    std::fs::write(path, dumps(v, &Dump::default()).unwrap()).unwrap();
}

/// What `og lookup report` prints under "in the reports" for one gene or
/// rsID, on an empty gVCF.
fn lookup_says(work: &Dir, summary: &Value, term: &str) -> String {
    let spath = work.path("lookup-summary.json");
    write_json(&spath, summary);
    let gvcf = work.path("empty.g.vcf");
    std::fs::write(
        &gvcf,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
    )
    .unwrap();
    let kind = if term.starts_with("rs") {
        "site"
    } else {
        "gene"
    };
    let args: Vec<String> = [
        "report",
        "--query",
        term,
        "--kind",
        kind,
        "--chrom",
        "chr1",
        "--start",
        "1",
        "--end",
        "1",
        "--gvcf",
        &gvcf,
        "--summary",
        &spath,
    ]
    .iter()
    .map(|x| (*x).to_owned())
    .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_lookup::cli(&args, &mut out, &mut err);
    let out = String::from_utf8(out).unwrap();
    match out.split_once("  in the reports:\n") {
        Some((_, rest)) if code == 0 => rest.to_owned(),
        _ => format!(
            "og lookup report failed: {code} {out} {}",
            String::from_utf8_lossy(&err)
        ),
    }
}

fn run(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|x| (*x).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

#[test]
fn the_index_holds_every_gene_and_rsid_with_its_sections() {
    let work = Dir::new("terms");
    let tsv = sample(&work, "tsv");
    let sum = collect(&tsv);
    assert_eq!(
        ["clinical", "slivar"].map(|k| s(get(&sum, &["sections", k, "source"]))),
        [
            "clinical/SYNTH_clinical_summary.tsv",
            "slivar/SYNTH_slivar_summary.tsv"
        ]
    );
    let idx = index(&sum, Some(&tsv));
    assert_eq!(
        terms(&idx),
        [
            "BRCA1",
            "BRCA2",
            "CFTR",
            "CYP2C19",
            "CYP2D6",
            "rs1000002",
            "rs123"
        ]
    );
    assert_eq!(get(&idx, &["schema_version"]), &Value::Int(1));
    assert_eq!(s(get(&idx, &["sample"])), "SYNTH");
    // One item per section with its report anchor, one mention per variant:
    // the secondary findings and the table both name chr13:2000.
    let want = [
        ("clinvar", "clinvar"),
        ("carrier", "carrier-list"),
        ("clinical", "clinical-filter"),
        ("slivar", "slivar"),
    ]
    .map(|(a, b)| (a.to_owned(), b.to_owned(), "1".to_owned()));
    assert_eq!(sections(&idx, "BRCA2"), want);
    // Existing_variation 'rs123&COSV1' gives rs123 only.
    assert!(!terms(&idx).iter().any(|t| t.eq_ignore_ascii_case("cosv1")));
}

#[test]
fn og_lookup_finds_every_term_the_index_finds_in_the_summary() {
    let work = Dir::new("lookup");
    let tsv = sample(&work, "tsv");
    let sum = collect(&tsv);
    let idx = index(&sum, Some(&tsv));
    let from_summary: Vec<String> = terms(&idx)
        .into_iter()
        .filter(|t| {
            items(&idx, t)
                .iter()
                .any(|it| !["clinical", "slivar"].contains(&s(get(it, &["section"])).as_str()))
        })
        .collect();
    assert_eq!(from_summary, ["BRCA2", "CFTR", "CYP2C19", "CYP2D6"]);
    for t in &from_summary {
        let said = lookup_says(&work, &sum, t);
        assert!(
            !said.contains("none of the reports mentions"),
            "{t}: {said}"
        );
        for it in items(&idx, t) {
            let key = s(get(it, &["section"]));
            if key == "clinical" || key == "slivar" {
                continue;
            }
            let title = s(get(it, &["title"]));
            assert!(said.contains(&title), "{t}: {title} missing from {said}");
        }
    }
    // A ClinVar hit's rsids: the index and the lookup both find it there,
    // and without them neither does (so the check above can fail).
    let mut rs = sum.clone();
    let Value::Array(hits) = get(&sum, &["sections", "clinvar", "data", "hits"]).clone() else {
        panic!("clinvar hits")
    };
    let mut hits = hits;
    let Value::Object(h0) = &mut hits[0] else {
        panic!("a hit")
    };
    h0.insert("rsids".to_owned(), json(r#"["rs7000001"]"#));
    set(
        &mut rs,
        &["sections", "clinvar", "data", "hits"],
        Value::Array(hits),
    );
    let ri = index(&rs, Some(&tsv));
    assert_eq!(sections(&ri, "rs7000001")[0].0, "clinvar");
    assert!(lookup_says(&work, &rs, "rs7000001").contains("ClinVar screen (step 06)"));
    assert!(
        lookup_says(&work, &sum, "rs7000001").contains("none of the reports mentions rs7000001")
    );
    assert!(items(&index(&sum, Some(&tsv)), "rs7000001").is_empty());
}

#[test]
fn a_prose_gene_field_and_a_locus_key_give_their_symbols() {
    let work = Dir::new("prose");
    let tsv = sample(&work, "tsv");
    let sum = collect(&tsv);
    let mut prose = sum.clone();
    set(
        &mut prose,
        &["sections", "mito", "state"],
        Value::Str("ok".into()),
    );
    set(
        &mut prose,
        &["sections", "mito", "data"],
        json(r#"{"variants": [{"gene": "SMN exon 7-8 deletion"}, {"gene": "VDBP (GC)"}]}"#),
    );
    let pt = terms(&index(&prose, Some(&tsv)));
    for t in ["SMN", "VDBP", "GC"] {
        assert!(pt.contains(&t.to_owned()), "{t} not in {pt:?}");
    }
    for t in ["EXON", "DELETION", "7-8", "(GC)", "(GC"] {
        assert!(!pt.contains(&t.to_owned()), "{t} in {pt:?}");
    }
    assert!(!lookup_says(&work, &prose, "GC").contains("none of the reports mentions"));

    let mut loc = sum.clone();
    for (k, data) in [
        (
            "expansions",
            r#"{"flagged": [{"locus": "HTT", "status": "pathogenic"}], "key_loci": [{"locus": "FMR1", "repeat_count": "."}],
                "readings": [{"locus": "UGT1A1", "repeat": "promoter TA repeat"}]}"#,
        ),
        (
            "locityper",
            r#"{"loci": [{"locus": "FCGR", "region": "chr1:161505430-161678654", "status": "called"},
                         {"locus": "HLA-DRB1", "region": "chr6:32484348-32593538", "status": "called"}]}"#,
        ),
    ] {
        set(&mut loc, &["sections", k, "state"], Value::Str("ok".into()));
        set(&mut loc, &["sections", k, "data"], json(data));
    }
    let li = index(&loc, Some(&tsv));
    for (t, sec) in [
        ("HTT", "expansions"),
        ("FMR1", "expansions"),
        ("UGT1A1", "expansions"),
        ("FCGR", "locityper"),
        ("HLA-DRB1", "locityper"),
    ] {
        assert!(
            sections(&li, t).iter().any(|x| x.0 == sec),
            "{t}: {:?}",
            sections(&li, t)
        );
    }
    assert!(
        !terms(&li).iter().any(|t| t.starts_with("CHR")),
        "a region is no term"
    );
    assert!(!lookup_says(&work, &loc, "HTT").contains("none of the reports mentions"));
}

#[test]
fn the_tsv_and_the_vcf_form_give_the_same_index() {
    let work = Dir::new("forms");
    let (tsv, vcf) = (sample(&work, "tsv"), sample(&work, "vcf"));
    let (st, sv) = (collect(&tsv), collect(&vcf));
    assert_eq!(
        ["clinical", "slivar"].map(|k| s(get(&sv, &["sections", k, "source"]))),
        [
            "clinical/SYNTH_clinical.vcf.gz",
            "slivar/SYNTH_prioritized.vcf.gz"
        ]
    );
    let (a, b) = (index(&st, Some(&tsv)), index(&sv, Some(&vcf)));
    for key in ["clinical", "slivar"] {
        assert!(!table_entries(&a, key).is_empty(), "{key}");
        assert_eq!(table_entries(&a, key), table_entries(&b, key), "{key}");
    }
    let names = |k| {
        table_entries(&a, k)
            .into_iter()
            .map(|x| x.0)
            .collect::<Vec<_>>()
    };
    assert_eq!(names("clinical"), ["BRCA1", "BRCA2"]);
    assert_eq!(names("slivar"), ["BRCA2", "CFTR", "rs1000002", "rs123"]);
    // Without the sample folder no table is read (red for the check above).
    assert!(table_entries(&index(&sv, None), "slivar").is_empty());
    // The whole index is the same, the file names and the time aside.
    let strip = |mut v: Value| {
        set(&mut v, &["generated_utc"], Value::Null);
        if let Value::Object(o) = &mut v
            && let Some(Value::Object(t)) = o.get_mut("terms")
        {
            for list in t.values_mut() {
                if let Value::Array(items) = list {
                    for it in items {
                        set(it, &["source"], Value::Null);
                    }
                }
            }
        }
        v
    };
    assert_eq!(strip(a), strip(b));
}

#[test]
fn only_ok_and_stale_sections_count_and_a_bad_table_adds_nothing() {
    let work = Dir::new("states");
    let tsv = sample(&work, "tsv");
    let sum = collect(&tsv);
    let mut st = sum.clone();
    set(
        &mut st,
        &["sections", "cpic", "state"],
        Value::Str("stale".into()),
    );
    for (k, state, data) in [
        (
            "mito",
            "failed",
            r#"{"variants": [{"gene": "MT-ND4", "variant": "chrM:11778 G>A"}]}"#,
        ),
        (
            "lpa",
            "missing",
            r#"{"variants": [{"rsid": "rs10455872", "gene": "LPA"}]}"#,
        ),
        ("haplogroup", "unreadable", r#"{"gene": "MT-CYB"}"#),
    ] {
        set(&mut st, &["sections", k, "state"], Value::Str(state.into()));
        set(&mut st, &["sections", k, "data"], json(data));
    }
    let si = index(&st, Some(&tsv));
    let states: Vec<String> = items(&si, "CYP2D6")
        .iter()
        .map(|it| s(get(it, &["state"])))
        .collect();
    assert_eq!(states, ["stale"]);
    for t in ["MT-ND4", "LPA", "rs10455872", "MT-CYB"] {
        assert!(!terms(&si).contains(&t.to_owned()), "{t}");
    }
    let mut ok = st.clone();
    set(
        &mut ok,
        &["sections", "mito", "state"],
        Value::Str("ok".into()),
    );
    let oi = index(&ok, Some(&tsv));
    assert_eq!(
        sections(&oi, "MT-ND4")[0].0,
        "mito",
        "red: the same data in an ok section"
    );
    let it = &items(&oi, "MT-ND4")[0];
    let anchor = og_render::chapters::section_of("mito").map_or("null", |x| x.id);
    assert_eq!(s(get(it, &["anchor"])), anchor);
    assert_eq!(s(get(it, &["count"])), "1");
    assert_eq!(
        get(it, &["mentions"]),
        &json(r#"[{"text": "gene=MT-ND4, variant=chrM:11778 G>A", "locus": "chrM:11778"}]"#)
    );

    // The summary of a readable run, the tables of a damaged copy: the
    // slivar VCF cut short after step 31 ran.
    let vcf = sample(&work, "vcf");
    let sv = collect(&vcf);
    let gz = Path::new(&vcf).join("slivar/SYNTH_prioritized.vcf.gz");
    let head = std::fs::read(&gz).unwrap()[..40].to_vec();
    std::fs::write(&gz, head).unwrap();
    let mut warn = Vec::new();
    let bi = build(&sv, Some(&vcf), &[], &mut warn).expect("the index is still written");
    let warn = String::from_utf8(warn).unwrap();
    assert!(
        warn.starts_with(
            "[WARN] slivar: cannot read slivar/SYNTH_prioritized.vcf.gz for the search index: "
        ),
        "{warn}"
    );
    assert!(table_entries(&bi, "slivar").is_empty());
    assert!(!table_entries(&bi, "clinical").is_empty());
}

#[test]
fn a_query_is_an_rsid_a_gene_or_its_prefix_or_a_position() {
    let work = Dir::new("query");
    let tsv = sample(&work, "tsv");
    let idx = index(&collect(&tsv), Some(&tsv));
    let a = query(&idx, "RS123");
    assert_eq!((a.kind, a.term.as_deref()), (Kind::Rsid, Some("rs123")));
    assert_eq!(a.items.len(), 1);
    assert_eq!(s(get(&a.items[0], &["section"])), "slivar");
    let Value::Array(m) = get(&a.items[0], &["mentions"]) else {
        panic!("mentions")
    };
    assert_eq!(s(get(&m[0], &["locus"])), "chr7:3000");
    let a = query(&idx, "brca");
    assert_eq!((a.kind, a.items.len()), (Kind::Gene, 0));
    assert_eq!(a.suggestions, ["BRCA1", "BRCA2"]);
    let a = query(&idx, "brca2");
    let secs: Vec<String> = a.items.iter().map(|it| s(get(it, &["section"]))).collect();
    assert_eq!(secs, ["clinvar", "carrier", "clinical", "slivar"]);
    assert!(a.suggestions.is_empty());
    let a = query(&idx, "chr7:1,000");
    assert_eq!(a.kind, Kind::Site);
    assert!(
        a.items.is_empty()
            && a.hint
                .unwrap()
                .contains("owngenome-lookup.sh <sample> chr7:1,000")
    );
    let a = query(&idx, "TP53");
    assert_eq!((a.kind, a.items.len()), (Kind::Gene, 0));
    assert_eq!(query_lines(&a), ["none of these results mentions TP53"]);
    let a = query(&idx, "two words");
    assert_eq!(a.kind, Kind::Error);
    assert_eq!(
        a.hint.as_deref(),
        Some(
            "not an rsID (rs123), a position (chr7:117559590) or a gene symbol (CFTR): 'two words'"
        )
    );
    assert_eq!(query(&idx, "chr1:0").kind, Kind::Error);
}

#[test]
fn the_command_line_builds_and_queries_the_index() {
    let work = Dir::new("cli");
    let tsv = sample(&work, "tsv");
    let sum = collect(&tsv);
    let spath = work.path("summary.json");
    write_json(&spath, &sum);
    let out = work.path("idx.json");
    let (code, stdout, stderr) = run(&[
        "build",
        "--summary",
        &spath,
        "--sample-dir",
        &tsv,
        "-o",
        &out,
    ]);
    assert_eq!(
        (code, stdout.as_str(), stderr.as_str()),
        (0, "7 terms\n", "")
    );
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(
        text.starts_with("{\n \"schema_version\": 1,\n \"sample\": \"SYNTH\""),
        "{text}"
    );
    assert!(text.ends_with("}\n"));
    assert_eq!(loads(&text).unwrap(), index(&sum, Some(&tsv)));

    let (code, stdout, _) = run(&["query", "--index", &out, "rs123"]);
    assert_eq!(code, 0);
    assert!(
        stdout.starts_with("rs123 is in 1 section(s):\n  Variant prioritization"),
        "{stdout}"
    );
    assert!(stdout.contains("    - chr7:3000 T>C\n"), "{stdout}");
    let (code, stdout, _) = run(&["query", "--index", &out, "brc"]);
    assert_eq!(code, 0);
    assert!(
        stdout.ends_with("terms that start with it: BRCA1, BRCA2\n"),
        "{stdout}"
    );
    let (code, stdout, _) = run(&["query", "--index", &out, "7:117559590"]);
    assert_eq!(code, 2);
    assert!(stdout.contains("owngenome-lookup.sh"), "{stdout}");
    let (code, _, stderr) = run(&["build", "--summary", &spath]);
    assert_eq!(code, 2, "{stderr}");
    let (code, _, stderr) = run(&["build", "--summary", &work.path("absent.json"), "-o", &out]);
    assert_eq!(code, 1, "{stderr}");
    let (code, _, stderr) = run(&["serve"]);
    assert!(
        code == 2 && stderr.contains("invalid choice: 'serve'"),
        "{stderr}"
    );
}

#[test]
fn the_veil_tags_the_items_that_name_a_closed_result() {
    let work = Dir::new("veil");
    let tsv = sample(&work, "tsv");
    let sum = collect(&tsv);
    let spans = [
        Closable {
            id: "sf-BRCA2".to_owned(),
            gene: Some("BRCA2".to_owned()),
            locus: None,
        },
        Closable {
            id: "clinvar-chr7-1000".to_owned(),
            gene: None,
            locus: Some(("chr7".to_owned(), 1000, 1000)),
        },
    ];
    let idx = build(&sum, Some(&tsv), &spans, &mut std::io::sink()).unwrap();
    let veil = |t: &str| -> Vec<String> {
        items(&idx, t)
            .iter()
            .map(|it| s(get(it, &["veil"])))
            .collect()
    };
    assert_eq!(veil("BRCA2"), ["sf-BRCA2"; 4]);
    // Each item by its own mentions: CFTR's ClinVar hit and its carrier row
    // name chr7:1000 itself, while its slivar row at chr7:3000 names no
    // closed result, so the locus the other sections name does not veil it.
    let cftr: Vec<(String, String)> = items(&idx, "CFTR")
        .iter()
        .map(|it| (s(get(it, &["section"])), s(get(it, &["veil"]))))
        .collect();
    assert_eq!(
        cftr,
        [
            ("clinvar".to_owned(), "clinvar-chr7-1000".to_owned()),
            ("carrier".to_owned(), "clinvar-chr7-1000".to_owned()),
            ("slivar".to_owned(), "null".to_owned()),
        ]
    );
    assert_eq!(veil("CYP2D6"), ["null"]);
    // No closed result, no field.
    assert!(
        !dumps(&index(&sum, Some(&tsv)), &Dump::default())
            .unwrap()
            .contains("\"veil\"")
    );
}

#[test]
fn search_js_answers_each_query_as_the_rust_rule_does() {
    let js_path = repo().join("assets/live/plugins/search.js");
    let js = std::fs::read_to_string(&js_path).unwrap();
    // The page asks for nothing but the index.
    assert_eq!(js.matches("fetch(").count(), 1);
    for bad in [
        "XMLHttpRequest",
        "sendBeacon",
        "WebSocket",
        "eval(",
        "innerHTML",
    ] {
        assert!(!js.contains(bad), "search.js uses {bad}");
    }
    let node = Command::new("node").arg("--version").output();
    if node.is_err() {
        assert!(
            std::env::var_os("CI").is_none(),
            "node is not on the PATH, and CI needs it for this check"
        );
        eprintln!("note: node is not on the PATH, so search.js is not run against the Rust rule");
        return;
    }
    let work = Dir::new("node");
    let tsv = sample(&work, "tsv");
    let idx = index(&collect(&tsv), Some(&tsv));
    let ipath = work.path("idx.json");
    write_json(&ipath, &idx);
    let queries = [
        "RS123",
        "rs1000002",
        "brca",
        "BRCA2",
        "cyp2",
        "C",
        "TP53",
        "chr7:1,000",
        "13-2000",
        "chr1:0",
        "chrMT_5",
        "two words",
        "",
        "  cftr ",
        "HLA-DRB1",
        "rs",
        "rs12x",
        "-x",
    ];
    let prog = "const m=require(process.argv[1]).match;\
        const idx=JSON.parse(require('fs').readFileSync(process.argv[2]));\
        const qs=JSON.parse(process.argv[3]);\
        console.log(JSON.stringify(qs.map(q=>{const r=m(idx,q);\
        return [r.kind,r.term,r.items.map(i=>i.section),r.suggestions]})))";
    let qs = Value::Array(
        queries
            .iter()
            .map(|q| Value::Str((*q).to_owned()))
            .collect(),
    );
    let out = Command::new("node")
        .args(["-e", prog, &js_path.display().to_string(), &ipath])
        .arg(dumps(&qs, &Dump::default()).unwrap())
        .output()
        .expect("node runs");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let got = loads(String::from_utf8(out.stdout).unwrap().trim()).unwrap();
    let want = Value::Array(
        queries
            .iter()
            .map(|q| {
                let a = query(&idx, q);
                Value::Array(vec![
                    Value::Str(a.kind.word().to_owned()),
                    a.term.map_or(Value::Null, Value::Str),
                    Value::Array(
                        a.items
                            .iter()
                            .map(|it| get(it, &["section"]).clone())
                            .collect(),
                    ),
                    Value::Array(a.suggestions.into_iter().map(Value::Str).collect()),
                ])
            })
            .collect(),
    );
    let Value::Array(g) = &got else {
        panic!("{got:?}")
    };
    let Value::Array(w) = &want else {
        unreachable!()
    };
    for ((q, g), w) in queries.iter().zip(g).zip(w) {
        assert_eq!(
            g, w,
            "query {q:?}: search.js (left) and the Rust rule (right) differ"
        );
    }
    assert_eq!(g.len(), queries.len());
}

#[test]
fn a_step_6_card_alone_is_indexed_and_a_clinical_card_once() {
    // Step 6's card with no step 23: the clinical section holds no card, so
    // the reports show sections.acmg_sf, and its genes must reach the index.
    let alone = json(
        r#"{"sections": {"acmg_sf": {"state": "ok", "step": "06", "title": "Secondary findings",
            "data": {"rows": [{"gene": "MUTYH", "variant": "chr1:45332803 C>T"}]}}}}"#,
    );
    let ai = index(&alone, None);
    let anchor = og_render::chapters::section_of("acmg_sf").map_or("null", |x| x.id);
    assert_eq!(
        sections(&ai, "MUTYH"),
        [("acmg_sf".to_owned(), anchor.to_owned(), "1".to_owned())],
        "red on main: the card of step 6 alone was left out of the index"
    );
    // With step 23's card in the clinical section the same gene is listed
    // once, under the clinical section, as the reports show it there.
    let both = json(
        r#"{"sections": {"acmg_sf": {"state": "ok", "step": "06", "title": "Secondary findings",
            "data": {"rows": [{"gene": "MUTYH", "variant": "chr1:45332803 C>T"}]}},
            "clinical": {"state": "ok", "step": "23", "title": "Clinical filter",
            "data": {"acmg_sf": {"rows": [{"gene": "MUTYH", "variant": "chr1:45332803 C>T"}]}}}}}"#,
    );
    let bi = index(&both, None);
    let keys: Vec<String> = sections(&bi, "MUTYH").into_iter().map(|x| x.0).collect();
    assert_eq!(keys, ["clinical"]);
    // A clinical section that did not run (or failed) is not read, so the
    // card of step 6 is.
    let mut failed = both.clone();
    set(
        &mut failed,
        &["sections", "clinical", "state"],
        Value::Str("failed".into()),
    );
    let fi = index(&failed, None);
    let keys: Vec<String> = sections(&fi, "MUTYH").into_iter().map(|x| x.0).collect();
    assert_eq!(keys, ["acmg_sf"]);
    // A clinical card that is empty or false is no card, as in the reports
    // and on the live site (both test it for truth): step 6's card is read.
    for empty in [Value::Object(Default::default()), Value::Bool(false)] {
        let mut e = both.clone();
        set(&mut e, &["sections", "clinical", "data", "acmg_sf"], empty.clone());
        let ei = index(&e, None);
        let keys: Vec<String> = sections(&ei, "MUTYH").into_iter().map(|x| x.0).collect();
        assert_eq!(keys, ["acmg_sf"], "clinical.data.acmg_sf = {empty:?}");
    }
}
