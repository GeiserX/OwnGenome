//! The journal's wording, caveats, levels, lead and widgets on synthetic runs.

#[path = "port_live_server_util.rs"]
mod util;

use flate2::read::GzDecoder;
use og_render::{chapters::CHAPTERS, notes};
use pycompat::json::{Value, dumps, loads};
use regex::Regex;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use util::{Work, repo};

const JOURNAL: &str = r#"{"sample":"JOURNAL","sections":{
    "cpsr":{"state":"ok","data":{"plp":[{"gene":"GENEA","variant":"c.100del",
        "final":"Pathogenic","clinvar":"Pathogenic"}]}},
    "clinvar":{"state":"ok","data":{"hits":[{"gene":"GENEB","chrom":"chr2","pos":1000,
        "ref":"C","alt":"T","significance":"Pathogenic","stars":2,"genotype":"het","af":0.0001}]}},
    "carrier":{"state":"ok","data":{"carriers":[
        {"gene":"GENEC","variant":"c.5A>G","condition":"Invented condition","significance":"Pathogenic",
         "stars":2,"inheritance":"AR","genotype":"het"},
        {"gene":"GENED","variant":"c.9C>T","condition":"Another invented condition","significance":"Pathogenic",
         "stars":2,"inheritance":"AR","genotype":"hom"}]}},
    "expansions":{"state":"ok","data":{"loci":[
        {"locus":"LOCUS1","repeat_count":"30/250","status":"pre_mutation,full_mutation",
         "ci":"29-31/240-260","allele_sizing":"spanning/flanking"}]}},
    "cpic":{"state":"ok","data":{"genes":[
        {"gene":"GENEE","diplotype":"*2/*2","phenotype":"Poor Metabolizer","status":"non-normal"},
        {"gene":"GENEF","diplotype":"*1/*1","phenotype":"Normal Metabolizer","status":"normal"}]}},
    "sample_qc":{"state":"ok","data":{"contamination":"warn","freemix":0.05}},
    "coverage":{"state":"ok","data":{"mean_depth":31.2}},
    "lpa":{"state":"failed","note":"invented failure","data":{}}
}}"#;

const CLINVAR_EMPTY: &str =
    r#"{"sections":{"clinvar":{"state":"ok","data":{"hits":[],"count":0}}}}"#;
const NO_SCREENS: &str = r#"{"sections":{"clinvar":{"state":"missing","data":{}}}}"#;
const LPA: &str = r#"{"sections":{"lpa":{"state":"ok","data":{"status":"estimated",
    "kiv2_copies":"36.9","ci95_low":"36.2","ci95_high":"37.6","class":"typical","flag":"."}}}}"#;

fn results(path: &Path, flags: &[&str]) -> Value {
    let mut args = vec![path.to_str().unwrap().to_owned()];
    args.extend(flags.iter().map(|arg| (*arg).to_owned()));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(
        og_live::results::cli(&args, &mut out, &mut err),
        0,
        "{}",
        String::from_utf8_lossy(&err)
    );
    loads(&String::from_utf8(out).unwrap()).unwrap()
}

fn summary(work: &Work, text: &str) -> PathBuf {
    let path = work.dir.join("summary.json");
    fs::write(&path, text).unwrap();
    path
}

fn field<'a>(doc: &'a Value, key: &str) -> &'a Value {
    let Value::Object(fields) = doc else {
        panic!("expected an object with {key}");
    };
    fields
        .get(key)
        .unwrap_or_else(|| panic!("missing {key} in {doc:?}"))
}

fn list(doc: &Value) -> &[Value] {
    let Value::Array(items) = doc else {
        panic!("expected an array, got {doc:?}");
    };
    items
}

fn text(doc: &Value) -> &str {
    let Value::Str(s) = doc else {
        panic!("expected text, got {doc:?}");
    };
    s
}

fn cards(doc: &Value) -> impl Iterator<Item = &Value> {
    list(field(doc, "chapters"))
        .iter()
        .flat_map(|chapter| list(field(chapter, "sections")))
}

fn card<'a>(doc: &'a Value, id: &str) -> &'a Value {
    cards(doc)
        .find(|card| text(field(card, "id")) == id)
        .unwrap_or_else(|| panic!("missing card {id}"))
}

/// Exempt exact caveats, links, raw summary data and quoted tool terms.
const DOCUMENT_EXEMPTIONS: &[&str] = &["src", "url", "caveat", "summary"];

struct Wording {
    forbidden: Regex,
    quoted: Regex,
}

impl Wording {
    fn new() -> Self {
        Self {
            forbidden: Regex::new(
                r"(?i)\b(ok|fine|healthy|all.clear|nothing found|normal range)\b",
            )
            .unwrap(),
            quoted: Regex::new(r#""[^"]*""#).unwrap(),
        }
    }

    fn check(&self, at: &str, value: &str) {
        let unquoted = self.quoted.replace_all(value, "");
        for word in self.forbidden.find_iter(&unquoted) {
            // The old check exempted compound terms such as fine-grained.
            let compound = word.as_str().eq_ignore_ascii_case("fine")
                && (unquoted[..word.start()].ends_with(['-', '/'])
                    || unquoted[word.end()..]
                        .to_ascii_lowercase()
                        .starts_with("-grained"));
            assert!(compound, "{at} contains all-clear wording: {value:?}");
        }
    }

    fn document(&self, at: &str, value: &Value, exemptions: &[&str]) {
        match value {
            Value::Object(fields) => {
                for (key, value) in fields {
                    if !exemptions.contains(&key.as_str()) {
                        self.document(&format!("{at}.{key}"), value, exemptions);
                    }
                }
            }
            Value::Array(items) => {
                for (i, value) in items.iter().enumerate() {
                    self.document(&format!("{at}[{i}]"), value, exemptions);
                }
            }
            Value::Str(value) => self.check(at, value),
            _ => {}
        }
    }
}

#[test]
fn results_overviews_and_widgets_use_no_all_clear_wording() {
    let work = Work::new("wording-results");
    let wording = Wording::new();
    for (name, input) in [("journal", JOURNAL), ("widgets", LPA)] {
        let path = summary(&work, input);
        for flags in [&[][..], &["--overview"][..]] {
            wording.document(name, &results(&path, flags), DOCUMENT_EXEMPTIONS);
        }
    }
    let doc = results(&summary(&work, LPA), &[]);
    wording.document("Lp(a) widget", field(card(&doc, "lpa"), "widget"), &[]);
}

#[test]
fn family_page_uses_no_all_clear_wording() {
    let dir = repo().join("tests/fixtures/families_live/family_a");
    let family = og_live::family::family_page("family_a", &dir).unwrap();
    assert_ne!(field(&family, "layout"), &Value::Null);
    Wording::new().document("family.json", &family, DOCUMENT_EXEMPTIONS);
}

fn asset_files(dir: &Path) -> Vec<PathBuf> {
    fs::read_dir(dir)
        .unwrap()
        .flat_map(|entry| {
            let entry = entry.unwrap();
            if entry.file_type().unwrap().is_dir() {
                asset_files(&entry.path())
            } else {
                vec![entry.path()]
            }
        })
        .collect()
}

#[test]
fn live_asset_text_uses_no_all_clear_wording() {
    let html_code = Regex::new(r"(?is)<(?:script|style)\b.*?</(?:script|style)>").unwrap();
    let tags = Regex::new(r"<[^>]*>").unwrap();
    // Consume comments as whole tokens so their quotes are ignored.
    let js = Regex::new(
        r#"(?ms)/\*.*?\*/|^[\t ]*//[^\n]*|'((?:\\.|[^'\\\n])*)'|"((?:\\.|[^"\\\n])*)"|`((?:\\.|[^`\\])*)`"#,
    )
    .unwrap();
    let wording = Wording::new();
    let mut checked = 0;
    for path in asset_files(&repo().join("assets/live")) {
        let at = path.to_str().unwrap();
        match path.extension().and_then(|ext| ext.to_str()) {
            Some("html") => {
                let source = fs::read_to_string(&path).unwrap();
                let visible = html_code.replace_all(&source, " ");
                for segment in tags.split(&visible) {
                    wording.check(at, segment);
                }
            }
            Some("js") => {
                let source = fs::read_to_string(&path).unwrap();
                for token in js.captures_iter(&source) {
                    if let Some(literal) = (1..=3).find_map(|group| token.get(group)) {
                        wording.check(at, literal.as_str());
                    }
                }
            }
            _ => continue,
        }
        checked += 1;
    }
    assert!(checked >= 10);
}

fn workspace_og() -> Option<PathBuf> {
    std::env::var_os("OG")
        .filter(|path| !path.is_empty())
        .map(PathBuf::from)
        .or_else(|| {
            let exe = std::env::current_exe().unwrap();
            let og = exe.parent()?.parent()?.join("og");
            og.is_file().then_some(og)
        })
}

fn checked(output: Output) -> Vec<u8> {
    assert!(
        output.status.success(),
        "{}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

#[test]
fn demo_results_overviews_and_map_documents_use_no_all_clear_wording() {
    let Some(og) = workspace_og() else {
        assert_ne!(std::env::var("CI").as_deref(), Ok("true"), "build og in CI");
        eprintln!("[SKIP] build og to check the demo results and map wording");
        return;
    };
    let work = Work::new("wording-demo");
    let root = work.root();
    let repo = repo();
    let reference = root.join("reference");
    fs::create_dir_all(&reference).unwrap();
    for file in ["gencode.v50.basic.genes.gtf", "gencode.v50.basic.exons.bed"] {
        fs::copy(
            repo.join("tests/fixtures/map/gencode-slice/reference")
                .join(file),
            reference.join(file),
        )
        .unwrap();
    }

    let wording = Wording::new();
    for (sample, seed) in [("DEMO-001", "1"), ("DEMO-003", "3")] {
        checked(
            Command::new(&og)
                .args(["dev", "demo-sample", "--genome-dir"])
                .arg(&root)
                .args([
                    "--records",
                    "2000",
                    "--sample",
                    sample,
                    "--seed",
                    seed,
                    "--repo",
                ])
                .arg(&repo)
                .output()
                .unwrap(),
        );
        let path = root.join(sample).join("summary.json");
        checked(
            Command::new(&og)
                .args(["summary", "collect", "--sample", sample, "--sample-dir"])
                .arg(root.join(sample))
                .arg("--assets")
                .arg(repo.join("assets"))
                .output()
                .unwrap(),
        );
        let doc = results(&path, &[]);
        assert_eq!(text(field(&doc, "sample")), sample);
        wording.document(sample, &doc, DOCUMENT_EXEMPTIONS);
        wording.document(
            sample,
            &results(&path, &["--overview"]),
            DOCUMENT_EXEMPTIONS,
        );
        for card in cards(&doc) {
            if let Value::Object(fields) = card
                && let Some(widget) = fields.get("widget")
            {
                wording.document(text(field(card, "id")), widget, &[]);
            }
        }

        let input = format!(
            r#"{{"results":{}}}"#,
            dumps(&doc, &Default::default()).unwrap()
        );
        let mut map = Command::new(&og)
            .args(["map", "build", "--root"])
            .arg(&root)
            .args(["--run", sample, "--assets"])
            .arg(repo.join("assets"))
            .arg("--cache-dir")
            .arg(work.dir.join("map-cache"))
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        map.stdin
            .take()
            .unwrap()
            .write_all(input.as_bytes())
            .unwrap();
        let gzip = checked(map.wait_with_output().unwrap());
        let mut decoded = String::new();
        GzDecoder::new(gzip.as_slice())
            .read_to_string(&mut decoded)
            .unwrap();
        let doc = loads(&decoded).unwrap();
        // Layer state is a protocol code.
        wording.document(&format!("{sample} map"), &doc, &["state"]);
    }
}

#[test]
fn every_report_section_has_a_named_caveat() {
    for chapter in CHAPTERS {
        for section in chapter.sections {
            assert!(
                !notes::section_notes(section.id).is_empty(),
                "{}",
                section.id
            );
        }
    }
    for (id, names) in notes::SECTION_NOTES {
        assert!(
            CHAPTERS
                .iter()
                .any(|ch| ch.sections.iter().any(|s| s.id == id)),
            "{id}"
        );
        for name in names {
            assert!(
                notes::text(name).is_some_and(|text| !text.trim().is_empty()),
                "{id}: {name}"
            );
        }
    }
}

#[test]
fn cards_keep_the_report_order_and_caveats_word_for_word() {
    let work = Work::new("journal-caveats");
    let doc = results(&summary(&work, JOURNAL), &[]);
    let chapters = list(field(&doc, "chapters"));
    assert_eq!(chapters.len(), CHAPTERS.len());
    for (chapter, expected) in chapters.iter().zip(CHAPTERS) {
        assert_eq!(text(field(chapter, "key")), expected.id);
        let sections = list(field(chapter, "sections"));
        assert_eq!(sections.len(), expected.sections.len(), "{}", expected.id);
        for (card, section) in sections.iter().zip(expected.sections) {
            assert_eq!(text(field(card, "id")), section.id);
            let want: Vec<_> = notes::section_notes(section.id)
                .iter()
                .map(|name| notes::text(name).unwrap())
                .collect();
            let got: Vec<_> = list(field(card, "caveat")).iter().map(text).collect();
            assert!(!got.is_empty(), "{}", section.id);
            assert_eq!(got, want, "{}", section.id);
        }
    }
}

#[test]
fn clinvar_has_a_table_for_hits_and_a_caveat_for_no_hits() {
    let work = Work::new("journal-clinvar");
    let hits = results(&summary(&work, JOURNAL), &[]);
    let table = field(card(&hits, "clinvar"), "table");
    let rows = list(field(table, "rows"));
    assert_eq!(rows.len(), 1);
    assert_eq!(text(&list(&rows[0])[0]), "GENEB");

    let empty = results(&summary(&work, CLINVAR_EMPTY), &[]);
    let clinvar = card(&empty, "clinvar");
    assert_eq!(field(clinvar, "table"), &Value::Null);
    assert_eq!(
        text(field(clinvar, "empty_text")),
        "No P/LP match in this ClinVar file at the call floor. ClinVar lists only part of what can cause disease, so this is information, and the caveats below still apply."
    );
}

#[test]
fn journal_rules_assign_levels_and_use_safe_finding_ids() {
    let work = Work::new("journal-levels");
    let doc = results(&summary(&work, JOURNAL), &[]);
    let findings = list(field(&doc, "findings"));
    let want = [
        ("cpsr-", "act", 1),
        ("clinvar-", "review", 4),
        ("carrier-genec", "carrier", 3),
        ("carrier-gened", "act", 1),
        ("str-locus1", "act", 1),
        ("pgx-genee", "medicine", 2),
        ("contamination", "trust", 5),
    ];
    assert_eq!(findings.len(), want.len());
    for (prefix, level, section) in want {
        let found: Vec<_> = findings
            .iter()
            .filter(|f| text(field(f, "id")).starts_with(prefix))
            .collect();
        assert_eq!(found.len(), 1, "{prefix}");
        let id = text(field(found[0], "id"));
        assert!(og_live::plugins::valid_finding_id(id), "{id}");
        assert_eq!(text(field(found[0], "level")), level, "{id}");
        assert_eq!(field(found[0], "section"), &Value::Int(section), "{id}");
        if id == "str-locus1" {
            assert!(text(field(found[0], "title")).contains("full-mutation"));
        }
    }
    for (id, level) in [
        ("cpsr", Some("act")),
        ("cpic", Some("medicine")),
        ("carrier-list", Some("act")),
        ("clinvar", Some("review")),
        ("qc", Some("trust")),
        ("lpa", None),
    ] {
        let want = level.map_or(Value::Null, |level| Value::Str(level.to_owned()));
        assert_eq!(field(card(&doc, id), "level"), &want, "{id}");
    }
}

#[test]
fn lead_puts_data_quality_first_and_names_the_findings_for_a_doctor() {
    let work = Work::new("journal-lead-quality");
    let path = summary(&work, JOURNAL);
    let doc = results(&path, &[]);
    let lead = results(&path, &["--overview"]);
    assert_eq!(field(&lead, "screened"), &Value::Bool(true));
    let open = text(field(&lead, "lead"));
    assert_eq!(
        open,
        concat!(
            "Read the data-quality section first: Another person's DNA may be mixed in. ",
            "Of everything this genome was screened for, three findings should go to a doctor: ",
            "Cancer predisposition (GENEA), LOCUS1 repeat in the full-mutation range and ",
            "Another invented condition: two changed copies. ",
            "One medicine result bears on how drugs work (GENEE): tell whoever prescribes. ",
            "A carrier of Invented condition carrier: carriers are not affected; ",
            "it matters when a partner carries the same gene. ",
            "One result is worth a look: Pathogenic change in GENEB."
        )
    );
    let closed = text(field(&lead, "lead_closed"));
    for f in list(field(&doc, "findings")) {
        if text(field(f, "level")) == "act" {
            let title = text(field(f, "title"));
            assert!(open.contains(title), "{open}");
            assert!(!closed.contains(title), "{closed}");
        }
    }
}

#[test]
fn lead_distinguishes_a_completed_screen_from_a_screen_with_no_result() {
    let work = Work::new("journal-lead-screen");
    for (input, screened, want) in [
        (
            CLINVAR_EMPTY,
            true,
            "No finding in what these steps assessed needs a doctor; sections 8 and 9 say what they do not see.",
        ),
        (
            NO_SCREENS,
            false,
            "The clinical screens have no result yet; section 9 lists the steps that did not run.",
        ),
    ] {
        let path = summary(&work, input);
        let overview = results(&path, &["--overview"]);
        assert_eq!(field(&overview, "screened"), &Value::Bool(screened));
        assert_eq!(text(field(&overview, "lead")), want);
        assert_eq!(text(field(&overview, "lead_closed")), want);
    }
}

#[test]
fn lpa_draws_an_interval_only_for_an_estimated_count() {
    let work = Work::new("journal-lpa");
    let doc = results(&summary(&work, LPA), &[]);
    let widget = field(card(&doc, "lpa"), "widget");
    assert_eq!(text(field(widget, "kind")), "ci");
    let (low, count, high) = (36.2, 36.9, 37.6);
    assert_eq!(field(widget, "low"), &Value::Float(low));
    assert_eq!(field(widget, "value"), &Value::Float(count));
    assert_eq!(field(widget, "high"), &Value::Float(high));
    assert!(matches!(field(widget, "max"), Value::Int(max) if *max >= 38));

    // Keep the numeric inputs so only the status can suppress the widget.
    let no_estimate = LPA.replace("\"estimated\"", "\"too little depth\"");
    let doc = results(&summary(&work, &no_estimate), &[]);
    let Value::Object(fields) = card(&doc, "lpa") else {
        panic!("expected the Lp(a) card");
    };
    assert!(
        fields
            .get("widget")
            .is_none_or(|value| value == &Value::Null)
    );
}

#[test]
fn journal_adapters_keep_the_report_details_and_unavailable_states() {
    let work = Work::new("journal-adapters");
    let input = r#"{"sections":{
        "ancient":{"state":"ok","data":{"nearest":[],"sites_used":1,"panel_sites":2,"plot_svg":"PHN2Zz48L3N2Zz4="}},
        "slivar":{"state":"ok","data":{"prioritized":2,"compound_het_variants":2,"compound_het_genes":1,
            "compound_het_pairs":1,"compound_het_phase":{"trans":1},
            "compound_het_trans":[{"gene":"GENE2","variants":"v1 v2","reads":3}]}},
        "exomiser":{"state":"ok","data":{"genes":[{"rank":1,"gene":"GENE3",
            "variants":[{"hgvs":"c.1A>G"},{"hgvs":"c.2C>T"}]}]}},
        "clinical":{"state":"ok","data":{"acmg_sf":{"unavailable":"no VEP output"}}},
        "carrier":{"state":"ok","data":{"unavailable":"no ClinVar screen"}},
        "cpic":{"state":"ok","data":{"parse_failed":true,"genes":[],"warnings":[]}},
        "annotsv":{"state":"ok","data":{"count":1,"rows":[{"type":"DEL","sv":"chr1:100-5000","gt":"0/1",
            "length":4900,"acmg_class":5,"genes":["GENE1"],"exon_edge":"yes","exon_edge_detail":"GENE1 exon 2",
            "population_af":null,"population_source":"","common":false,"support":null}]}},
        "traits":{"state":"ok","data":{"rows":[],"rules":[{"trait":"Rule trait","gene":"GENE4","rsids":"rs1,rs2",
            "status":"called","result":"AB","phase":"read_backed","reading":"rule reading"}]}},
        "clinvar":{"state":"ok","data":{"hits":[],"deep_intronic":[
            {"gene":"GENE5","hgvs":"c.100+500A>G","chrom":"chr1","pos":1,"ref":"A","alt":"G",
             "call":"het","detail":"GQ 40"}]}}
    }}"#;
    let doc = results(&summary(&work, input), &[]);
    let image = field(card(&doc, "ancient"), "image");
    assert_eq!(
        text(field(image, "src")),
        "data:image/svg+xml;base64,PHN2Zz48L3N2Zz4="
    );
    let exomiser = field(card(&doc, "exomiser"), "table");
    assert_eq!(list(field(exomiser, "rows")).len(), 2);
    for (id, phrase) in [
        ("slivar", "in trans on 3 reads"),
        ("clinical-filter", "not assessed"),
        ("carrier-list", "not made"),
        ("cpic", "could not be parsed"),
        ("structural-variants", "on exon edges"),
        ("traits-panel", "Rule trait"),
        ("clinvar", "GENE5"),
    ] {
        let rendered = dumps(card(&doc, id), &Default::default()).unwrap();
        assert!(
            rendered.contains(phrase),
            "{id} lacked {phrase:?}: {rendered}"
        );
    }
    let sf = dumps(
        field(card(&doc, "clinical-filter"), "facts"),
        &Default::default(),
    )
    .unwrap();
    assert!(!sf.contains("ClinVar P/LP in its genes"));
    for id in ["carrier-list", "cpic"] {
        assert!(list(field(card(&doc, id), "facts")).is_empty(), "{id}");
    }
}

#[test]
fn apoe_text_stays_behind_its_warning() {
    let work = Work::new("journal-apoe");
    let input = r#"{"sections":{"traits":{"state":"ok","data":{"rows":[],"apoe":[
        {"rsid":"rs1","gene":"APOE","trait":"Invented marker","genotype":"A/T","reading":"invented reading"}
    ]}}}}"#;
    let doc = results(&summary(&work, input), &[]);
    let traits = card(&doc, "traits-panel");
    let hidden = field(traits, "hidden");
    assert_eq!(
        text(field(hidden, "summary")),
        og_render::sections::apoe::APOE_WARNING
    );
    assert_eq!(
        text(field(hidden, "text")),
        "Invented marker (rs1): A/T. invented reading"
    );
    for key in ["facts", "lines", "table"] {
        let visible = dumps(field(traits, key), &Default::default()).unwrap();
        assert!(!visible.contains("Invented marker"), "{key}");
    }
}

#[test]
fn widgets_keep_their_source_values_and_have_browser_drawers() {
    let work = Work::new("journal-widget-values");
    let input = r#"{"sections":{
        "expansions":{"state":"ok","data":{"stranger":true,"loci":[
            {"locus":"LOCUS1","repeat_count":"17/42","status":"normal,full_mutation","normal_max":"26",
             "pathologic_min":"40","ci":"17-17/40-45","sizing":"sized","allele_sizing":"sized/sized"},
            {"locus":"LOCUS2","repeat_count":"30","status":"normal","sizing":"flanking_only","allele_sizing":"flanking_only"}
        ]}},
        "chip":{"state":"ok","data":{"calls":[
            {"gene":"GENE1","protein_change":"p.A1T","vaf":"0.081","alt_reads":"4","depth":"49","note":"invented"},
            {"gene":"GENE2","protein_change":"p.A2T","vaf":".","alt_reads":"3","depth":"40","note":"invented"}
        ],"genes":[
            {"gene":"GENE1","detection_limit_vaf":"0.133","status":"covered"},
            {"gene":"GENE2","detection_limit_vaf":"0.142","status":"covered"}
        ]}},
        "paralogs":{"state":"ok","data":{"genes":[
            {"gene":"SMN1","copy_number":"2","status":"called"},
            {"gene":"SMN2","copy_number":"1","status":"called"},
            {"gene":"HBA1","copy_number":"2.4","status":"no_call"},
            {"gene":"HBA2","copy_number":"2","status":"called"},
            {"gene":"CYP21A2 (RCCX)","copy_number":"4","status":"called"},
            {"gene":"GENE3","copy_number":"2","status":"called"}
        ]}},
        "parascopy":{"state":"ok","data":{"rows":1,"smn1":{"cn":1,"text":"resolved"},
            "smn2":{"cn":null,"text":"not resolved"}}},
        "y_tree":{"state":"ok","data":{"haplogroup":"TEST-TIP","tree":"Invented tree","path":[
            {"Level":"1","Branch":"ROOT","Derived":"11","Ancestral":"0","No_call":"1"},
            {"Level":"2","Branch":"TEST-TIP","Derived":"9","Ancestral":"1","No_call":"0"}
        ]}}
    }}"#;
    let doc = results(&summary(&work, input), &[]);
    let lpa = results(&summary(&work, LPA), &[]);
    let source = fs::read_to_string(repo().join("assets/live/person.js")).unwrap();
    let drawers = source
        .split_once("var WIDGETS = {")
        .unwrap()
        .1
        .split_once('}')
        .unwrap()
        .0;
    let wording = Wording::new();
    for (doc, id, kind, drawer) in [
        (&doc, "expansions", "str", "strGauge"),
        (&lpa, "lpa", "ci", "ciBar"),
        (&doc, "chip", "vaf", "vafChart"),
        (&doc, "paralogs", "copies", "copyDots"),
        (&doc, "y-chromosome", "ypath", "yPath"),
    ] {
        let widget = field(card(doc, id), "widget");
        assert_eq!(text(field(widget, "kind")), kind, "{id}");
        assert!(drawers.contains(&format!("{kind}: {drawer}")), "{id}");
        assert!(source.contains(&format!("function {drawer}(")), "{id}");
        wording.document(id, widget, &[]);
    }

    let expansions = card(&doc, "expansions");
    let loci = list(field(field(expansions, "widget"), "loci"));
    assert_eq!(loci.len(), 2);
    let sized = &loci[0];
    for (key, want) in [
        ("locus", "LOCUS1"),
        ("repeat_count", "17/42"),
        ("normal_max", "26"),
        ("pathologic_min", "40"),
        ("ci", "17-17/40-45"),
    ] {
        assert_eq!(text(field(sized, key)), want, "{key}");
    }
    assert_ne!(text(field(&loci[1], "status")), "normal");
    let rows = list(field(field(expansions, "table"), "rows"));
    assert!(list(&rows[1]).contains(&Value::Str("not sized".into())));

    let chip = field(card(&doc, "chip"), "widget");
    let genes = list(field(chip, "genes"));
    assert_eq!(genes.len(), 2);
    for (gene, name, limit) in [(0, "GENE1", 0.133), (1, "GENE2", 0.142)] {
        assert_eq!(text(field(&genes[gene], "gene")), name);
        assert_eq!(field(&genes[gene], "limit"), &Value::Float(limit));
    }
    let calls = list(field(chip, "calls"));
    assert_eq!(calls.len(), 1);
    assert_eq!(text(field(&calls[0], "gene")), "GENE1");
    assert_eq!(field(&calls[0], "vaf"), &Value::Float(0.081));
    assert_eq!(field(&calls[0], "limit"), &Value::Float(0.133));

    let copies = field(field(card(&doc, "paralogs"), "widget"), "genes");
    let want = loads(
        r#"[
        {"gene":"SMN1","copies":2,"status":"called"},
        {"gene":"SMN2","copies":1,"status":"called"},
        {"gene":"HBA1","copies":null,"status":"no_call"},
        {"gene":"HBA2","copies":2,"status":"called"},
        {"gene":"CYP21A2 (RCCX)","copies":4,"status":"called"},
        {"gene":"SMN1 (Parascopy)","copies":1,"status":"called"},
        {"gene":"SMN2 (Parascopy)","copies":null,"status":"not resolved"}
    ]"#,
    )
    .unwrap();
    assert_eq!(copies, &want);

    let y = field(card(&doc, "y-chromosome"), "widget");
    assert_eq!(text(field(y, "tip")), "TEST-TIP");
    let want = loads(
        r#"[
        {"level":"1","branch":"ROOT","derived":11.0,"ancestral":0.0,"no_call":1.0},
        {"level":"2","branch":"TEST-TIP","derived":9.0,"ancestral":1.0,"no_call":0.0}
    ]"#,
    )
    .unwrap();
    assert_eq!(field(y, "path"), &want);
}
