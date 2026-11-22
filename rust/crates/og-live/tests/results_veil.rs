use pycompat::json::{Value, loads};
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};

fn input(summary: &str) -> PathBuf {
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    let path = std::env::temp_dir().join(format!(
        "og-results-test-{}-{}.json",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::Relaxed)
    ));
    fs::write(&path, summary).unwrap();
    path
}

fn results(summary: &str, mode: &[&str]) -> Value {
    let path = input(summary);
    let mut args = vec![path.to_string_lossy().into_owned()];
    args.extend(mode.iter().map(|arg| (*arg).to_owned()));
    let mut out = Vec::new();
    let mut err = Vec::new();
    let code = og_live::results::cli(&args, &mut out, &mut err);
    fs::remove_file(path).unwrap();
    assert_eq!(code, 0, "{}", String::from_utf8_lossy(&err));
    loads(&String::from_utf8(out).unwrap()).unwrap()
}

fn field<'a>(doc: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(doc, |node, key| {
        let Value::Object(object) = node else {
            panic!("expected object at {key}")
        };
        object
            .get(*key)
            .unwrap_or_else(|| panic!("missing key {key}"))
    })
}

fn text_list(values: &[Option<&str>]) -> Value {
    Value::Array(
        values
            .iter()
            .map(|item| item.map_or(Value::Null, |text| Value::Str(text.into())))
            .collect(),
    )
}

const CLOSED: &str = r#"{
  "sample":"TEST",
  "sections":{
    "clinical":{"state":"ok","data":{"acmg_sf":{"clinvar_hits":[
      {"gene":"GENE1","variant":"chr2:110 A>G","genotype":"het","significance":"Pathogenic"}
    ],"high_impact":[]}}},
    "clinvar":{"state":"ok","data":{"count":3,"hits":[
      {"gene":"GENE1","chrom":"chr2","pos":110,"ref":"A","alt":"G","stars":2},
      {"gene":"GENE1","chrom":"chr2","pos":130,"ref":"C","alt":"T","stars":1},
      {"gene":"GENE2","chrom":"chr3","pos":7,"ref":"T","alt":"C","stars":1}
    ],"deep_intronic":[
      {"gene":"GENE1","chrom":"chr2","pos":123,"ref":"A","alt":"T","call":"het"},
      {"gene":"GENE2","chrom":"chr3","pos":8,"ref":"C","alt":"A","call":"het"}
    ]}},
    "annotsv":{"state":"ok","data":{"consensus":true,"count":1,"rows":[
      {"type":"DEL","sv":"chr2:100-200","length":100,"acmg_class":"5","genes":["GENE1"],"support":null}
    ]}},
    "exomiser":{"state":"ok","data":{"genes":[{"gene":"GENE1","rank":1},{"gene":"GENE2","rank":2}]}},
    "chip":{"state":"ok","data":{"genes":[{"gene":"GENE1","mean_depth":75,"detection_limit_vaf":0.06}],
      "calls":[{"gene":"GENE2","vaf":0.08},{"gene":"GENE1","vaf":0.12}]}},
    "carrier":{"state":"ok","data":{"two_hits":{"genes":["SMN1"]}}},
    "paralogs":{"state":"ok","data":{"genes":[{"gene":"SMN1","copy_number":1,"status":"called"},
      {"gene":"SMN2","copy_number":2,"status":"called"}]}},
    "parascopy":{"state":"ok","data":{"smn1":{"text":"1","cn":1},"smn1_marker":{"status":"not assessed"}}}
  }
}"#;

#[test]
fn veil_distinguishes_exact_loci_spans_and_unrelated_rows() {
    let doc = results(CLOSED, &[]);
    let veil = field(&doc, &["veil"]);
    let sf = Some("sf-gene1-chr2-110-a-g");
    let sv = Some("sv-chr2-100-200");
    assert_eq!(
        field(veil, &["rows", "clinvar"]),
        &text_list(&[sf, sf, None])
    );
    assert_eq!(
        field(veil, &["viewer", "clinvar"]),
        &text_list(&[sf, sf, None])
    );
    assert_eq!(
        field(veil, &["rows", "structural-variants"]),
        &text_list(&[sv])
    );
    assert_eq!(field(veil, &["rows", "exomiser"]), &text_list(&[sf, None]));
    assert_eq!(field(veil, &["lines", "clinvar"]), &text_list(&[sf, None]));
    assert_eq!(
        field(veil, &["findings", "exomiser-gene1"]),
        &Value::Str(sf.unwrap().into())
    );
    assert!(matches!(field(veil, &["marks", "genome"]), Value::Object(m) if !m.is_empty()));
    assert!(matches!(field(veil, &["spans"]), Value::Array(a) if a.len() >= 3));
}

#[test]
fn veil_tags_widget_items_facts_and_closed_lead() {
    let doc = results(CLOSED, &[]);
    let veil = field(&doc, &["veil"]);
    assert_eq!(
        field(veil, &["widget", "chip"]),
        &text_list(&[None, Some("sf-gene1-chr2-110-a-g")])
    );
    assert_eq!(
        field(veil, &["widget", "paralogs"]),
        &text_list(&[Some("two-hits-smn1"), None])
    );
    assert_eq!(
        field(veil, &["facts", "paralogs"]),
        &text_list(&[Some("two-hits-smn1"), Some("two-hits-smn1")])
    );
    let overview = results(CLOSED, &["--overview"]);
    assert!(matches!(field(&overview, &["lead_closed"]), Value::Str(t)
        if t.contains("closed until opened on purpose") && !t.contains("GENE1")));
    let Value::Array(top) = field(&overview, &["top"]) else {
        panic!("missing top")
    };
    assert!(
        !top.iter()
            .any(|f| field(f, &["id"]) == &Value::Str("exomiser-gene1".into()))
    );
}

#[test]
fn no_act_has_no_tags_and_gene_only_closes_the_viewer_link() {
    let base = r#"{"sections":{"clinvar":{"state":"ok","data":{"hits":[
        {"gene":"GENE1","chrom":"chr1","pos":7}
    ]}}}}"#;
    assert_eq!(field(&results(base, &[]), &["veil"]), &loads("{}").unwrap());
    let with_gene_only = base.replace(
        r#""sections":{"#,
        r#""sections":{
        "carrier":{"state":"ok","data":{"two_hits":{"genes":["GENE1"]}}},"#,
    );
    let doc = results(&with_gene_only, &[]);
    assert_eq!(
        field(&doc, &["veil", "viewer", "clinvar"]),
        &text_list(&[Some("two-hits-gene1")])
    );
}

#[test]
fn discordant_smn1_counts_use_the_report_verdict_without_carrier_finding() {
    let doc = results(
        r#"{"sections":{
      "paralogs":{"state":"ok","data":{"genes":[{"gene":"SMN1","copy_number":1,"status":"called","detail":"."}],
        "smn1":{"agree":"disagree","paraphase":1,"parascopy":{"text":"2","cn":2}}}},
      "parascopy":{"state":"ok","data":{"smn1":{"text":"2","cn":2}}}
    }}"#,
        &[],
    );
    let Value::Array(findings) = field(&doc, &["findings"]) else {
        panic!("missing findings")
    };
    assert!(
        !findings
            .iter()
            .any(|f| field(f, &["id"]) == &Value::Str("smn1".into()))
    );
    assert!(
        pycompat::json::dumps(&doc, &Default::default())
            .unwrap()
            .contains("they DISAGREE: check both before using either")
    );
}
