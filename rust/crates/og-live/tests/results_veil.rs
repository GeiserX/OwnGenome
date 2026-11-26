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
      {"gene":"GENE1","chrom":"chr2","pos":123,"ref":"A","alt":"T","call":"het","detail":"GQ 40"},
      {"gene":"GENE2","chrom":"chr3","pos":8,"ref":"C","alt":"A","call":"het","detail":"GQ 40"}
    ]}},
    "annotsv":{"state":"ok","data":{"consensus":true,"count":1,"rows":[
      {"type":"DEL","sv":"chr2:100-200","length":100,"acmg_class":"5","genes":["GENE1"],"support":null}
    ]}},
    "roh":{"state":"ok","data":{"autosomal_over_5mb":[{"region":"chr2:100-9000000","mb":9}]}},"exomiser":{"state":"ok","data":{"genes":[{"gene":"GENE1","rank":1},{"gene":"GENE2","rank":2}]}},
    "chip":{"state":"ok","data":{"genes":[{"gene":"GENE1","mean_depth":75,"detection_limit_vaf":0.06}],
      "calls":[{"gene":"GENE2","vaf":0.08,"protein_change":"p.A1T","depth":75,"alt_reads":6,"note":"clone"},{"gene":"GENE1","vaf":0.12,"protein_change":"p.G2A","depth":75,"alt_reads":9,"note":"clone"}]}},
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
        &Value::Str("sf-gene1-chr2-110-a-g".into())
    );
    assert_eq!(
        field(veil, &["marks", "genome"]),
        &loads(r#"{"hit-0":"sf-gene1-chr2-110-a-g","hit-1":"sf-gene1-chr2-110-a-g","sv-0":"sv-chr2-100-200"}"#).unwrap()
    );
    assert_eq!(
        field(veil, &["marks", "clinvar-stars"]),
        &loads(r#"{"stars-2":"sf-gene1-chr2-110-a-g","stars-1":"sf-gene1-chr2-110-a-g"}"#).unwrap()
    );
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

#[test]
fn search_index_tags_only_the_items_that_name_a_closed_result() {
    let request = format!(
        r#"{{"summary":{CLOSED},"terms":{{
          "GENE1":[{{"section":"clinical","mentions":[]}},{{"section":"exomiser","mentions":[]}}],
          "GENE10":[{{"section":"traits","mentions":[{{"text":"GENE10","locus":"chr2:201"}}]}}],
          "rs123":[
            {{"section":"clinvar","mentions":[{{"text":"At the exact site","locus":"2:110"}}]}},
            {{"section":"slivar","mentions":[{{"text":"Inside the deletion","locus":"chr2:155"}}]}},
            {{"section":"traits","mentions":[{{"text":"Outside the deletion","locus":"chr2:201"}}]}}
          ],
          "rs456":[{{"section":"cpsr","mentions":[{{"text":"A change in GENE1.","locus":null}}]}}]
        }}}}"#,
    );
    let got = results(&request, &["--tag-search"]);
    let groups = [
        ("GENE1", vec![Some("sf-gene1-chr2-110-a-g"); 2]),
        ("GENE10", vec![None]),
        (
            "rs123",
            vec![Some("sf-gene1-chr2-110-a-g"), Some("sv-chr2-100-200"), None],
        ),
        ("rs456", vec![Some("sf-gene1-chr2-110-a-g")]),
    ];
    for (term, expected) in groups {
        let Value::Array(items) = field(&got, &[term]) else {
            panic!("expected entries for {term}");
        };
        assert_eq!(items.len(), expected.len());
        for (item, expected) in items.iter().zip(expected) {
            let Value::Object(item) = item else {
                panic!("expected index item");
            };
            assert_eq!(
                item.get("veil").cloned(),
                expected.map(|id| Value::Str(id.into())),
                "{term}"
            );
        }
    }
}

#[test]
fn unreadable_actionable_rows_keep_a_finding_at_the_card() {
    let doc = results(
        r#"{"sections":{
          "annotsv":{"state":"ok","data":{"rows":[{"sv":"chr1:1000-9000","type":"deletion",
            "length":-8000,"acmg_class":"5","genes":["GENEX","GENEY","GENEZ",null],"gt":"0/1"}]}},
          "cpsr":{"state":"ok","data":{"plp":["GENEY c.1del"]}},
          "carrier":{"state":"ok","data":{"carriers":[["GENEX"]]}},
          "clinvar":{"state":"ok","data":{"hits":[["GENEZ"]]}}
        }}"#,
        &[],
    );
    let Value::Array(found) = field(&doc, &["findings"]) else {
        panic!("expected findings");
    };
    for (card, level) in [
        ("structural-variants", "act"),
        ("cpsr", "act"),
        ("carrier-list", "act"),
        ("clinvar", "review"),
    ] {
        let f = found
            .iter()
            .find(|f| field(f, &["card"]) == &Value::Str(card.into()))
            .unwrap();
        assert_eq!(field(f, &["level"]), &Value::Str(level.into()), "{card}");
        assert_eq!(field(f, &["gene"]), &Value::Null, "{card}");
        for key in ["title", "text"] {
            assert!(
                !matches!(field(f, &[key]), Value::Str(text) if text.contains("GENE")),
                "{card} {key}"
            );
        }
    }
}

#[test]
fn smn1_marker_card_keeps_the_report_facts_and_caveat() {
    let doc = results(
        r#"{"sections":{"parascopy":{"state":"ok","data":{
          "smn1":{"text":"2 (quality 33, filter PASS)","cn":2},
          "smn1_marker":{"status":"present","alt_reads":5,"depth":40,"copies":0.5}
        }}}}"#,
        &[],
    );
    let Value::Array(chapters) = field(&doc, &["chapters"]) else {
        panic!("expected chapters");
    };
    let card = chapters
        .iter()
        .filter_map(|ch| match field(ch, &["sections"]) {
            Value::Array(cards) => Some(cards),
            _ => None,
        })
        .flatten()
        .find(|c| field(c, &["id"]) == &Value::Str("paralogs".into()))
        .expect("paralogs card");
    assert_eq!(
        field(card, &["facts"]),
        &loads(r#"[["SMN1 copies","2 (quality 33, filter PASS)"],["SMN1 g.27134T>G silent-carrier marker","present, 5 of 40 reads carry it, about 0.5 copies"]]"#).unwrap()
    );
    let Value::Array(lines) = field(card, &["lines"]) else {
        panic!("expected card lines");
    };
    let caveat = og_render::sections::smn1::SMN1_MARKER_CAVEAT;
    assert!(lines.contains(&Value::Str(format!("SMN1 marker caveat: {caveat}"))));
}
