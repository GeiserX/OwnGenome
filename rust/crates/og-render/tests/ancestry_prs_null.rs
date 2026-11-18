//! A short row of pgsc_calc's table: og summary collect writes its score,
//! matched and total as null (the golden summary
//! tests/parity/golden/summary-collect/all/all.json has one, Short/PGS1).
//! Both reports render it, and the row reads "not available" in those
//! cells and, in the HTML report, in the result cell: with no score there
//! is no raw score either. bin/render_report.py stops the whole HTML report on that row
//! (`html.escape(None)`) and writes "None" in the text report, so this is a
//! known difference from the Python, outside the parity corpus.

use std::path::Path;

use og_render::model::Summary;
use og_render::{html, text};

fn golden() -> Summary {
    let p = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/parity/golden/summary-collect/all/all.json");
    Summary::parse(&std::fs::read_to_string(p).unwrap()).unwrap()
}

#[test]
fn a_null_score_reads_not_available_in_both_reports() {
    let s = golden();
    let page = html::report(&s).unwrap();
    assert!(
        page.contains(
            "<tr><td>Short</td><td>PGS1</td><td>not available</td><td>not available/not \
             available</td><td>not available</td>"
        ),
        "the HTML row"
    );
    assert!(
        page.contains("<tr><td>CAD</td><td>PGS000018</td><td>0.5</td><td>100/120</td>"),
        "the other rows are unchanged"
    );
    let txt = text::report(&s).unwrap();
    let row = txt
        .lines()
        .find(|l| l.starts_with("  Short "))
        .expect("the text row");
    assert_eq!(
        row.trim_end(),
        format!(
            "  {:<35} not available (not available/not available variants matched) PGS1",
            "Short"
        )
    );
    assert!(!txt.contains("None (None/None"), "{txt}");
}
