use pycompat::json::{Dump, Value, dumps, loads};
use pycompat::text::html_escape;
use std::fs;
use std::path::Path;

const EXAMPLE: &str =
    include_str!("../../../../tests/fixtures/pharmcat/pharmcat-docs-example.json");
const CREDIT: &str = "Recommendation text from CPIC (https://cpicpgx.org/, CC0 1.0) via PharmCAT. \
    PharmCAT also carries PharmGKB (ClinPGx) content under CC BY-SA 4.0";
const CREDIT_LICENSE: &str = "https://creativecommons.org/licenses/by-sa/4.0/";

fn cpic_report(case: &str, json: &str) -> String {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("pgx-attribution-{case}"));
    fs::create_dir_all(&dir).unwrap();
    let report = dir.join("report.json");
    fs::write(&report, json).unwrap();
    let args = [
        "cpic-report",
        "--sample",
        "DEMO-001",
        "--report",
        report.to_str().unwrap(),
        "--outdir",
        dir.to_str().unwrap(),
    ];
    assert_eq!(og_pgx::main(&args.map(str::to_owned)), 0);
    fs::read_to_string(dir.join("DEMO-001_cpic_recommendations.txt")).unwrap()
}

#[test]
fn credits_recommendation_sources_when_text_is_used() {
    let text = cpic_report("credit", EXAMPLE);
    assert!(
        text.lines()
            .any(|line| line.contains(CREDIT) && line.contains(CREDIT_LICENSE))
    );
}

#[test]
fn a_report_with_no_recommendation_text_has_no_pharmgkb_credit() {
    let Value::Object(mut data) = loads(EXAMPLE).unwrap() else {
        panic!("expected a PharmCAT report object");
    };
    data.shift_remove("drugs");
    let json = dumps(&Value::Object(data), &Dump::default()).unwrap();
    let text = cpic_report("no-text", &json);
    assert!(text.contains("CYP2D6 -- Intermediate Metabolizer"));
    assert!(!text.contains("PharmGKB"));
}

#[test]
fn cpic_drugs_without_recommendation_text_have_no_credit() {
    fn clear_recommendations(value: &mut Value) {
        match value {
            Value::Object(fields) => {
                for (key, field) in fields {
                    if key == "drugRecommendation" {
                        *field = Value::Null;
                    } else {
                        clear_recommendations(field);
                    }
                }
            }
            Value::Array(items) => items.iter_mut().for_each(clear_recommendations),
            _ => {}
        }
    }

    let mut data = loads(EXAMPLE).unwrap();
    clear_recommendations(&mut data);
    let json = dumps(&data, &Dump::default()).unwrap();
    let text = cpic_report("null-recommendations", &json);
    assert!(
        text.lines()
            .any(|line| line == "      - atomoxetine [CPIC, Moderate]")
    );
    assert!(!text.contains("PharmGKB"));
}

#[test]
fn text_report_decodes_the_pharmcat_threshold() {
    let text = cpic_report("plain-text", EXAMPLE);
    assert!(text.contains("If <200 ng/mL"));
    assert!(!text.contains("&lt;"));
}

#[test]
fn text_report_survives_one_html_escape_layer() {
    let text = cpic_report("html-escaped", EXAMPLE);
    let html = html_escape(&text);
    assert!(html.contains("If &lt;200 ng/mL"));
    assert!(!html.contains("If <200 ng/mL"));
    assert!(!html.contains("If &amp;lt;200 ng/mL"));
}

#[test]
fn decodes_the_five_entities_once_after_removing_literal_markup() {
    let json = EXAMPLE.replace(
        "If &lt;200 ng/mL",
        "If &lt;200 ng/mL and &gt;0 <b>&quot;A&amp;B&#39;s&quot;</b> (&amp;lt; and &nbsp;)",
    );
    let text = cpic_report("single-decode", &json);
    assert!(text.contains(r#"If <200 ng/mL and >0 "A&B's" (&lt; and &nbsp;)"#));
}
