//! A step the run selected that wrote nothing for the sample reads
//! "selected, no result", not "not run". The Nextflow report-only steps
//! retry once and are then dropped, so their section has no file, the same
//! as a step nobody asked for. HTML_REPORT writes the `--tools` list as the
//! `run tools` row of run_manifest.tsv; og summary collect maps each tool
//! to its sections through assets/steps.tsv, lists every selected section
//! that is `missing` under `selected_no_result`, labels its `not_run` line
//! and gives it a note. A run that records no selection (the bash runner)
//! reads as before.

use og_summary::collect::{Group, Input, post};
use pycompat::json::Value;
use std::fs;
use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!(
            "og-selected-no-result-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn put(&self, rel: &str, body: &str) {
        let p = self.0.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(&p, body).unwrap();
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn collect_with(d: &Dir, assets: &str) -> Value {
    post::collect(&Input {
        sample: "S".into(),
        sample_dir: d.0.to_str().unwrap().to_owned(),
        declared_sex: None,
        assets: assets.to_owned(),
        group: Group::All,
        tables: format!("{REPO}/bin"),
    })
    .unwrap()
}

fn collect(d: &Dir) -> Value {
    collect_with(d, &format!("{REPO}/assets"))
}

fn field<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(key),
        _ => panic!("not an object: {v:?}"),
    }
}

fn get<'a>(v: &'a Value, key: &str) -> &'a Value {
    field(v, key).unwrap_or_else(|| panic!("no {key} in {v:?}"))
}

fn texts(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a
            .iter()
            .map(|x| match x {
                Value::Str(s) => s.clone(),
                _ => panic!("not a string: {x:?}"),
            })
            .collect(),
        _ => panic!("not a list: {v:?}"),
    }
}

fn st(v: &Value) -> Option<&str> {
    match v {
        Value::Str(s) => Some(s),
        Value::Null => None,
        _ => panic!("not a string or null: {v:?}"),
    }
}

fn section<'a>(v: &'a Value, key: &str) -> &'a Value {
    get(get(v, "sections"), key)
}

fn title(v: &Value, key: &str) -> String {
    st(get(section(v, key), "title")).unwrap().to_owned()
}

/// A Nextflow run that selected ClinVar (step 06), VNtyper (step 52),
/// mosdepth (step 16b) and the report: only mosdepth's depth summary is
/// there, so step 06 and step 52 wrote nothing.
fn nextflow_folder(tag: &str) -> Dir {
    let d = Dir::new(tag);
    d.put(
        "run_manifest.tsv",
        "run\tpipeline_version\t1.0.0\nrun\ttools\tclinvar vntyper mosdepth html_report\n",
    );
    d.put(
        "coverage/S.mosdepth.summary.txt",
        "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t1\t1\t30.5\t0\t9\n",
    );
    d
}

#[test]
fn a_selected_step_with_no_folder_reads_selected_no_result() {
    let d = nextflow_folder("nf");
    let v = collect(&d);
    // Step 16b's other section, callability, has no table here: mosdepth
    // was selected and left it out, so it is listed too.
    assert_eq!(
        texts(get(&v, "selected_no_result")),
        ["callability", "clinvar", "vntyper"],
    );
    for (key, tool, step) in [("clinvar", "clinvar", "06"), ("vntyper", "vntyper", "52")] {
        let s = section(&v, key);
        assert_eq!(st(get(s, "state")), Some("missing"), "{key}");
        assert_eq!(
            st(get(s, "note")),
            Some(
                format!(
                    "selected, no result: --tools named {tool}, but step {step} wrote nothing for this sample"
                )
                .as_str()
            ),
            "{key}"
        );
    }
    let not_run = texts(get(&v, "not_run"));
    for key in ["clinvar", "vntyper"] {
        let line = format!("{} (selected, no result)", title(&v, key));
        assert!(not_run.contains(&line), "{line:?} not in {not_run:?}");
        assert!(!not_run.contains(&title(&v, key)), "{key}: plain line left");
    }
    // A step nobody selected still reads not run, with no note.
    assert!(not_run.contains(&title(&v, "telomere")), "{not_run:?}");
    assert_eq!(st(get(section(&v, "telomere"), "note")), None);
    // The section whose file is there is not listed.
    assert_ne!(st(get(section(&v, "coverage"), "state")), Some("missing"));
}

#[test]
fn a_section_that_is_no_step_of_its_own_is_not_listed() {
    // Step 06 feeds the carrier list too; its own section names the step,
    // as not_run does.
    let v = collect(&nextflow_folder("carrier"));
    assert!(!texts(get(&v, "selected_no_result")).contains(&"carrier".to_owned()));
    assert_eq!(st(get(section(&v, "carrier"), "note")), None);
}

#[test]
fn a_failed_step_keeps_its_failed_line() {
    let d = nextflow_folder("failed");
    d.put(
        "logs/run_status.tsv",
        "meta\tstarted_epoch\t1577836800\nstep\t06\tfailed (exit 1)\n",
    );
    let v = collect(&d);
    assert_eq!(st(get(section(&v, "clinvar"), "state")), Some("failed"));
    assert!(!texts(get(&v, "selected_no_result")).contains(&"clinvar".to_owned()));
    let line = format!("{} (step 06 failed)", title(&v, "clinvar"));
    assert!(texts(get(&v, "not_run")).contains(&line));
}

#[test]
fn a_run_that_records_no_selection_reads_as_before() {
    let d = Dir::new("bash");
    d.put("run_manifest.tsv", "run\tpipeline_version\t1.0.0\n");
    let v = collect(&d);
    assert!(field(&v, "selected_no_result").is_none());
    assert_eq!(st(get(section(&v, "clinvar"), "note")), None);
    assert!(texts(get(&v, "not_run")).contains(&title(&v, "clinvar")));
}

#[test]
fn a_selected_step_whose_files_are_there_is_not_listed() {
    let d = Dir::new("all-there");
    d.put("run_manifest.tsv", "run\ttools\tmosdepth\n");
    d.put(
        "coverage/S.mosdepth.summary.txt",
        "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t1\t1\t30.5\t0\t9\n",
    );
    d.put("coverage/S_callability.bed.gz", "");
    let v = collect(&d);
    assert_eq!(texts(get(&v, "selected_no_result")), Vec::<String>::new());
}

#[test]
fn the_catalogue_beside_the_binary_serves_an_assets_folder_without_one() {
    // An assets folder without steps.tsv (an older one): the catalogue of
    // the checkout the binary sits in is read instead, when there is one.
    let assets = Dir::new("assets-only");
    fs::copy(
        format!("{REPO}/assets/callability_genes.tsv"),
        assets.0.join("callability_genes.tsv"),
    )
    .unwrap();
    let v = collect_with(&nextflow_folder("old-assets"), assets.0.to_str().unwrap());
    if og_runlog::steps::default_path().is_some() {
        assert!(texts(get(&v, "selected_no_result")).contains(&"clinvar".to_owned()));
    } else {
        assert!(field(&v, "selected_no_result").is_none());
    }
}

#[test]
fn every_section_the_catalogue_names_is_a_summary_section() {
    // A section name in assets/steps.tsv that summary.json does not have
    // would never be marked.
    let steps =
        og_runlog::steps::read_steps(&PathBuf::from(format!("{REPO}/assets/steps.tsv"))).unwrap();
    for s in steps.iter().filter(|s| !s.tool.is_empty()) {
        for k in &s.sections {
            assert!(
                post::ORDER.contains(&k.as_str()),
                "step {} names section {k}, which summary.json does not have",
                s.step
            );
        }
    }
}
