//! The report's charts (og report charts), on the invented sections of
//! tests/fixtures/report_charts/sections.json and on changed copies of them.
//! Moved here from tests/test_report_charts.py, check for check; each check
//! runs on a broken input too and must catch it.
//!
//! 1. ClinVar: one bar per star 4 to 0 with the by_stars count, then a grey
//!    bar with the low-confidence count; no green; with no hit a grey frame
//!    says that no hit is not an all-clear. A hit past the 50-row table links
//!    to the section.
//! 2. PRS: no strip without a panel, only the report's note; with a panel a
//!    score with no percentile, or one of the other sex, says "no
//!    percentile" and gets no dot.
//! 3. Fine-grained ancestry: the groups are the ones the report lists, and
//!    each dot lies inside its range line.
//! 4. Genome: an unreadable run of homozygosity stops the chart; a run under
//!    5 MB is not drawn; hits on chrM and on alt contigs, SVs and MoChA calls
//!    with no place are counted in not_drawn; every mark is drawn inside its
//!    chromosome's outline.
//! 5. Mito: one dot per heteroplasmic site, the poly-C call marked apart.
//! 6. Every caveat is the report's own, and no caveat sentence is written in
//!    the charts' source.
//! 7. The charts JSON is what bin/report_charts.py wrote, byte for byte, on
//!    the fixture and its changed copies (frozen as hashes of the oracle's
//!    files), the same on every run, and the report paths write it.

use std::path::{Path, PathBuf};

use og_render::charts::{self, GRCH38, Kind, MarkKind, Spec, build, chrom_len, draw, figure, json};
use og_render::model::Summary;
use og_render::notes;
use og_render::sections::{ancestry_fine, prs};
use pycompat::json::{Dump, Value, dumps, loads};

const GREENS: [&str; 5] = ["#d5f5e3", "#1d3b2a", "#8ee0ae", "#196f3d", "green"];

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// The fixture as a summary value: sample S, a fixed clock and its sections.
fn fixture() -> Value {
    let text =
        std::fs::read_to_string(repo().join("tests/fixtures/report_charts/sections.json")).unwrap();
    let f = loads(&text).unwrap();
    let Value::Object(o) = f else { panic!() };
    let mut s = indexmap::IndexMap::new();
    s.insert("sample".to_owned(), Value::Str("S".into()));
    s.insert(
        "generated_utc".to_owned(),
        Value::Str("2026-01-01 00:00 UTC".into()),
    );
    s.insert("sections".to_owned(), o["sections"].clone());
    Value::Object(s)
}

/// `v[a][b]...` for changing a copy.
fn at<'a>(v: &'a mut Value, path: &[&str]) -> &'a mut Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(o) => o.get_mut(*k).unwrap_or_else(|| panic!("no {k}")),
        Value::Array(a) => &mut a[k.parse::<usize>().unwrap()],
        _ => panic!("{k}: not a container"),
    })
}

fn set(v: &mut Value, path: &[&str], x: &str) {
    let (last, head) = path.split_last().unwrap();
    let Value::Object(o) = at(v, head) else {
        panic!()
    };
    o.insert((*last).to_owned(), loads(x).unwrap());
}

fn summary(v: &Value) -> Summary {
    Summary::parse(&dumps(v, &Dump::default()).unwrap()).unwrap()
}

fn charts_of(v: &Value) -> Vec<Spec> {
    build(&summary(v)).expect("the charts")
}

fn chart(v: &Value, id: &str) -> Spec {
    charts_of(v)
        .into_iter()
        .find(|c| c.id == id)
        .unwrap_or_else(|| panic!("no chart {id}"))
}

fn data<'a>(v: &'a Value, sec: &str) -> &'a Value {
    let Value::Object(o) = v else { panic!() };
    let Value::Object(s) = &o["sections"] else {
        panic!()
    };
    let Value::Object(x) = &s[sec] else { panic!() };
    &x["data"]
}

fn int(v: &Value) -> f64 {
    match v {
        Value::Int(i) => *i as f64,
        Value::Float(x) => *x,
        Value::Str(s) => s.parse().unwrap(),
        _ => panic!("{v:?}"),
    }
}

fn num(v: &Value) -> Option<f64> {
    match v {
        Value::Null => None,
        v => Some(int(v)),
    }
}

/// The number in attribute `name` of an SVG element's text `piece`.
fn attr(piece: &str, name: &str) -> Option<f64> {
    let k = format!(" {name}=\"");
    let piece = format!(" {piece}");
    let i = piece.find(&k)? + k.len();
    piece[i..].split('"').next()?.parse().ok()
}

fn redraw(sp: &mut Spec) {
    sp.svg = draw(sp).unwrap();
}

// --- the detectors, each also run on a broken input ---------------------------------------

fn clinvar_problems(sp: &Spec, d: &Value) -> Vec<String> {
    let mut out = Vec::new();
    let Value::Object(o) = d else { panic!() };
    let Value::Object(by) = &o["by_stars"] else {
        panic!()
    };
    let mut want: Vec<f64> = ["4", "3", "2", "1", "0"]
        .iter()
        .map(|k| int(&by[*k]))
        .collect();
    want.push(int(&o["low_confidence"]));
    let got: Vec<f64> = sp.marks.iter().map(|m| int(&m.value)).collect();
    if got != want {
        out.push(format!(
            "bars {got:?} != by_stars and low_confidence {want:?}"
        ));
    }
    if sp.marks.last().unwrap().tone != charts::Tone::Muted {
        out.push("the low-confidence bar is not grey".into());
    }
    let low = sp.svg.to_lowercase();
    for g in GREENS {
        if low.contains(g) {
            out.push(format!("green {g} in the chart"));
        }
    }
    out
}

/// A dot per score with a percentile of the declared sex, "no percentile"
/// for the others.
fn prs_problems(sp: &Spec, scores: &[Value]) -> Vec<String> {
    let mut out = Vec::new();
    for (m, r) in sp.marks.iter().zip(scores) {
        let Value::Object(r) = r else { panic!() };
        let pct = &r["percentile"];
        let placeable = !matches!(pct, Value::Str(s) if s.is_empty())
            && !matches!(&r["sex_check"], Value::Str(s) if s == prs::NOT_APPLICABLE);
        if placeable && num(&m.value) != Some(int(pct)) {
            out.push(format!(
                "{}: value {:?} is not the percentile",
                m.id, m.value
            ));
        }
        if !placeable && num(&m.value).is_some() {
            out.push(format!(
                "{}: placed with no percentile of the declared sex",
                m.id
            ));
        }
    }
    let dots = sp.svg.matches("<circle ").count();
    let want = sp.marks.iter().filter(|m| num(&m.value).is_some()).count();
    if dots != want {
        out.push(format!("{dots} dots for {want} percentiles"));
    }
    if sp.svg.matches(">no percentile</text>").count() != sp.marks.len() - want {
        out.push("not every score without a percentile says 'no percentile'".into());
    }
    out
}

/// Each mark with a range: its dot between the two ends of its line.
fn whisker_problems(svg: &str) -> Vec<String> {
    let mut out = Vec::new();
    for g in svg.split("<a ").skip(1) {
        let g = g.split("</a>").next().unwrap();
        let (Some(line), Some(dot)) = (g.find("class=\"whisker\""), g.find("<circle ")) else {
            continue;
        };
        let (x1, x2) = (
            attr(&g[line..], "x1").unwrap(),
            attr(&g[line..], "x2").unwrap(),
        );
        let cx = attr(&g[dot..], "cx").unwrap();
        if !(x1 <= cx && cx <= x2) {
            out.push(format!("dot {cx} outside {x1}..{x2}"));
        }
    }
    out
}

fn with_chr(c: &str) -> String {
    if c.starts_with("chr") {
        c.to_owned()
    } else {
        format!("chr{c}")
    }
}

/// What the genome chart misses of the AnnotSV rows and MoChA calls on the
/// primary chromosomes.
fn sv_mca_problems(sp: &Spec, svs: &[(&str, bool)], mcas: &[bool]) -> Vec<String> {
    let mut out = Vec::new();
    let mark = |id: &str| sp.marks.iter().find(|m| m.id == id);
    for (i, (sv, _)) in svs.iter().enumerate() {
        if chrom_len(&with_chr(sv.split(':').next().unwrap())).is_none() {
            continue;
        }
        let id = format!("sv-{i}");
        if mark(&id).is_none_or(|m| m.anchor.as_deref() != Some("structural-variants"))
            || !sp.svg.contains(&format!("data-mark=\"{id}\""))
        {
            out.push(format!("AnnotSV row {sv} not drawn"));
        }
    }
    for (i, placeable) in mcas.iter().enumerate() {
        let id = format!("mca-{i}");
        if !placeable {
            if mark(&id).is_some() {
                out.push(format!(
                    "MoChA call {i} drawn, but it has no place on the chart"
                ));
            }
            continue;
        }
        if mark(&id).is_none_or(|m| m.anchor.as_deref() != Some("mosaic-mca"))
            || !sp.svg.contains(&format!("data-mark=\"{id}\""))
        {
            out.push(format!("MoChA call {i} not drawn"));
        }
    }
    if sp.svg.matches("class=\"sv\"").count() != svs.iter().filter(|(_, p)| *p).count() {
        out.push("the outlined spans do not match the rows".into());
    }
    if sp.svg.matches("class=\"mca\"").count() != mcas.iter().filter(|p| **p).count() {
        out.push("the bars do not match the calls".into());
    }
    out
}

/// The marks of an ideogram drawn past either end of their chromosome's
/// outline: each mark's y is its row's plus a fixed offset per kind.
fn outside_marks(svg: &str) -> Vec<String> {
    let mut rows = Vec::new();
    for p in svg.split("<rect ").skip(1) {
        let p = p.split("/>").next().unwrap();
        if p.contains("height=\"8\" rx=\"4\" fill=\"none\"") {
            let (x, y, w) = (
                attr(p, "x").unwrap(),
                attr(p, "y").unwrap(),
                attr(p, "width").unwrap(),
            );
            rows.push((y - 2.0, x, x + w));
        }
    }
    let mut out = Vec::new();
    for p in svg.split("<rect class=\"").skip(1) {
        let kind = p.split('"').next().unwrap();
        let off = match kind {
            "roh" => 2.0,
            "sv" | "tick" => -1.0,
            "mca" => 12.0,
            _ => continue,
        };
        let p = p.split("/>").next().unwrap();
        let (x, y, w) = (
            attr(p, "x").unwrap(),
            attr(p, "y").unwrap(),
            attr(p, "width").unwrap(),
        );
        match rows.iter().find(|(ry, _, _)| (ry - (y - off)).abs() < 0.01) {
            Some((_, lo, hi)) if x >= lo - 0.01 && x + w <= hi + 0.01 => {}
            _ => out.push(format!("{kind} at x={x} y={y} width={w}")),
        }
    }
    out
}

/// Caveat sentences written into a source file instead of taken from notes.
fn note_literals(src: &str) -> Vec<String> {
    let mut found: Vec<String> = notes::NOTES
        .iter()
        .filter(|(_, t)| t.len() > 30 && src.contains(&t[..30]))
        .map(|(k, _)| (*k).to_owned())
        .collect();
    if src.lines().any(|l| l.contains("_NOTE: &str =")) {
        found.push("a *_NOTE constant".into());
    }
    found
}

// --- the checks ---------------------------------------------------------------------------

#[test]
fn clinvar_bars_count_the_stars_in_yellow_and_grey() {
    let s = fixture();
    let cv = chart(&s, "clinvar-stars");
    let d = data(&s, "clinvar");
    assert_eq!(clinvar_problems(&cv, d), Vec::<String>::new());
    assert!(
        cv.marks
            .iter()
            .all(|m| matches!(m.tone, charts::Tone::Attention | charts::Tone::Muted))
    );
    // Red controls: a green bar, and a bar off by one.
    let mut green = cv.clone();
    green.svg = green
        .svg
        .replace(charts::Tone::Attention.colour(), "#d5f5e3");
    assert!(!clinvar_problems(&green, d).is_empty());
    let mut off = cv.clone();
    off.marks[1].value = Value::Int(int(&off.marks[1].value) as i128 + 1);
    assert!(!clinvar_problems(&off, d).is_empty());
}

#[test]
fn no_clinvar_hit_is_a_grey_frame_that_is_not_an_all_clear() {
    let mut z = fixture();
    set(
        &mut z,
        &["sections", "clinvar", "data", "by_stars"],
        r#"{"4":0,"3":0,"2":0,"1":0,"0":0}"#,
    );
    set(&mut z, &["sections", "clinvar", "data", "hits"], "[]");
    set(&mut z, &["sections", "clinvar", "data", "count"], "0");
    set(
        &mut z,
        &["sections", "clinvar", "data", "low_confidence"],
        "0",
    );
    let svg = chart(&z, "clinvar-stars").svg;
    assert!(svg.contains("No hit is not an all-clear"), "{svg}");
    assert!(!svg.contains("<a "));
    assert!(!GREENS.iter().any(|g| svg.to_lowercase().contains(g)));
}

#[test]
fn a_hit_past_the_table_links_to_the_section() {
    let s = fixture();
    let g = chart(&s, "genome");
    let Value::Object(d) = data(&s, "clinvar") else {
        panic!()
    };
    let Value::Array(hits) = &d["hits"] else {
        panic!()
    };
    let row = |locus: &str| {
        hits.iter()
            .position(|h| {
                let Value::Object(h) = h else { panic!() };
                format!(
                    "{}:{}",
                    pycompat::value::py_str(&h["chrom"]).unwrap(),
                    pycompat::value::py_str(&h["pos"]).unwrap()
                ) == locus
            })
            .unwrap()
    };
    let ticks: Vec<_> = g
        .marks
        .iter()
        .filter(|m| m.kind == Some(MarkKind::Hit))
        .collect();
    let want: Vec<&str> = ticks
        .iter()
        .map(|m| {
            if row(m.locus.as_deref().unwrap()) < 50 {
                "clinvar-hits"
            } else {
                "clinvar"
            }
        })
        .collect();
    let got: Vec<&str> = ticks.iter().map(|m| m.anchor.as_deref().unwrap()).collect();
    assert_eq!(got, want);
    assert!(want.contains(&"clinvar"));
}

#[test]
fn prs_strip_places_only_percentiles_of_the_declared_sex() {
    let s = fixture();
    let pr = chart(&s, "prs");
    let Value::Object(d) = data(&s, "prs") else {
        panic!()
    };
    let Value::Array(scores) = &d["scores"] else {
        panic!()
    };
    assert_eq!(pr.kind, Kind::Strip);
    assert_eq!(prs_problems(&pr, scores), Vec::<String>::new());
    // Red controls: the raw Score_SUM drawn as a dot, and the other sex's percentile.
    let field = |i: usize, k: &str| {
        let Value::Object(r) = &scores[i] else {
            panic!()
        };
        int(&r[k])
    };
    let mut raw = pr.clone();
    raw.marks[1].value = Value::Float(field(1, "score"));
    redraw(&mut raw);
    assert!(!prs_problems(&raw, scores).is_empty());
    let mut other = pr.clone();
    other.marks[2].value = Value::Float(field(2, "percentile"));
    redraw(&mut other);
    assert!(!prs_problems(&other, scores).is_empty());
}

#[test]
fn prs_without_a_panel_is_the_reports_note() {
    let mut na = fixture();
    set(&mut na, &["sections", "prs", "data", "adjusted"], "false");
    let nc = chart(&na, "prs");
    assert_eq!(nc.kind, Kind::Text);
    assert!(nc.marks.is_empty() && !nc.svg.contains("<circle"));
    assert_eq!(
        nc.text.as_deref(),
        Some(prs::note(data(&na, "prs")).unwrap().as_str())
    );
    assert_eq!(figure(&nc), None, "the HTML report leaves a text chart out");
}

#[test]
fn ancestry_groups_are_the_reports_with_each_dot_inside_its_range() {
    let s = fixture();
    let an = chart(&s, "ancestry-fine");
    let d = data(&s, "ancestry_fine");
    let mut want = Vec::new();
    for (level, _) in ancestry_fine::LEVELS {
        for r in ancestry_fine::rows(d, level).unwrap() {
            let Value::Object(r) = r else { panic!() };
            want.push(pycompat::value::py_str(&r["group"]).unwrap());
        }
    }
    let got: Vec<String> = an.marks.iter().map(|m| m.label.clone()).collect();
    assert_eq!(got, want);
    assert!(want.iter().any(|g| g == "Invented region B"));
    assert!(!want.iter().any(|g| g == "Invented group C"));
    assert!(an.svg.contains("whisker"));
    assert_eq!(whisker_problems(&an.svg), Vec::<String>::new());
    // Red control: a dot moved outside its range.
    let mut off = an.clone();
    off.marks[0].value = Value::Float(99.0);
    redraw(&mut off);
    assert!(!whisker_problems(&off.svg).is_empty());
}

#[test]
fn an_ancestry_range_with_one_bound_is_left_out_of_the_detail() {
    let mut half = fixture();
    set(
        &mut half,
        &["sections", "ancestry_fine", "data"],
        r#"{"continent": [{"group": "Invented group D", "fraction": 0.3, "low": 0.2}]}"#,
    );
    let m = &chart(&half, "ancestry-fine").marks[0];
    assert!(!m.detail.contains("bootstrap range"), "{}", m.detail);
    assert_eq!(m.high, Value::Null);
}

#[test]
fn genome_draws_runs_of_5_mb_and_counts_what_has_no_place() {
    let s = fixture();
    let g = chart(&s, "genome");
    let roh: Vec<&str> = g
        .marks
        .iter()
        .filter(|m| m.kind == Some(MarkKind::Roh))
        .map(|m| m.locus.as_deref().unwrap())
        .collect();
    assert_eq!(roh, ["chr3:10000000-22300000", "chr15:40000000-46100000"]);
    let nd = |what: &str| g.not_drawn.iter().find(|x| x.what == what).map(|x| x.count);
    assert_eq!(nd("runs of homozygosity under 5 MB"), Some(1));
    assert_eq!(nd("ClinVar hits on chrM"), Some(1));
    assert_eq!(nd("ClinVar hits on other contigs (alt, unplaced)"), Some(1));
    assert!(g.svg.contains("Not drawn:"));
    assert!(
        !g.marks
            .iter()
            .any(|m| matches!(m.chrom.as_deref(), Some("chrM" | "chr22_KI270879v1_alt")))
    );
    assert!(!g.svg.contains("AnnotSV class 4 or 5 SV") && !g.svg.contains("MoChA"));
    // The floor is the roh section's, the 5 MB the docs promise.
    assert_eq!(charts::ROH_MIN_MB, 5.0);
    // No genome chart when ROH did not run, even with ClinVar hits.
    let mut nr = s.clone();
    set(&mut nr, &["sections", "roh", "state"], "\"missing\"");
    assert!(!charts_of(&nr).iter().any(|c| c.id == "genome"));
}

#[test]
fn an_unreadable_run_of_homozygosity_stops_the_chart() {
    for region in ["chr3:10000000_22300000", "3:1-2", "chr3:20-10"] {
        let mut bad = fixture();
        set(
            &mut bad,
            &[
                "sections",
                "roh",
                "data",
                "autosomal_over_5mb",
                "0",
                "region",
            ],
            &format!("\"{region}\""),
        );
        let e = build(&summary(&bad)).expect_err(region);
        assert!(
            e.0.starts_with(&format!("run of homozygosity '{region}'")),
            "{e}"
        );
    }
    // Control: the readable regions do not stop it.
    assert!(build(&summary(&fixture())).is_ok());
}

const SV_ROWS: &str = r#"[
  {"sv": "chr7:117480000-117530000", "type": "DEL", "gt": "0/1", "length": "50000", "acmg_class": "4", "genes": ["CFTR"]},
  {"sv": "chr22_KI270879v1_alt:100-900", "type": "DUP", "gt": "0/1", "length": "800", "acmg_class": "5", "genes": []},
  {"sv": "17:7660000-7690000", "type": "DEL", "gt": "0/1", "length": "30000", "acmg_class": "5", "genes": ["TP53"]}]"#;
const MCA_CALLS: &str = r#"[
  {"chrom": "chr20", "start": "30000000", "end": "42000000", "length_mb": "12.0", "type": "Del", "cell_fraction": "0.12"},
  {"chrom": "chrM", "start": "100", "end": "900", "length_mb": ".", "type": "Del", "cell_fraction": "."},
  {"chrom": "chr3", "start": ".", "end": ".", "length_mb": ".", "type": "Dup", "cell_fraction": "."},
  {"chrom": "20", "start": "50000000", "end": "52000000", "length_mb": "2.0", "type": "Dup", "cell_fraction": "0.05"}]"#;

fn with_sv_mca(rows: &str, calls: &str) -> Value {
    let mut s = fixture();
    set(
        &mut s,
        &["sections", "annotsv"],
        &format!(r#"{{"state": "ok", "data": {{"count": 3, "rows": {rows}}}}}"#),
    );
    set(
        &mut s,
        &["sections", "mosaic"],
        &format!(r#"{{"state": "ok", "data": {{"status": "assessed", "calls": {calls}}}}}"#),
    );
    s
}

#[test]
fn genome_draws_annotsv_spans_and_mocha_bars() {
    let gm = chart(&with_sv_mca(SV_ROWS, MCA_CALLS), "genome");
    let svs = [
        ("chr7:117480000-117530000", true),
        ("chr22_KI270879v1_alt:100-900", false),
        ("17:7660000-7690000", true),
    ];
    let mcas = [true, false, false, true];
    assert_eq!(sv_mca_problems(&gm, &svs, &mcas), Vec::<String>::new());
    let nd = |what: &str| {
        gm.not_drawn
            .iter()
            .find(|x| x.what == what)
            .map(|x| x.count)
    };
    assert_eq!(nd("AnnotSV class 4 and 5 SVs on other contigs"), Some(1));
    assert_eq!(
        nd("MoChA calls on other contigs or with no position"),
        Some(2)
    );
    assert!(
        gm.svg.contains("AnnotSV class 4 or 5 SV") && gm.svg.contains("mosaic alteration (MoChA)")
    );
    // Red control: a genome chart that leaves out the MoChA calls.
    let mut lost = gm.clone();
    lost.marks.retain(|m| m.kind != Some(MarkKind::Mca));
    redraw(&mut lost);
    assert!(!sv_mca_problems(&lost, &svs, &mcas).is_empty());
}

#[test]
fn every_genome_mark_stays_inside_its_chromosome() {
    let edge = with_sv_mca(
        // past the end of chr21 (46,709,983 bp), and a 10 bp one at the very start of chr1
        r#"[{"sv": "chr21:46700000-46800000", "type": "DUP", "gt": "0/1", "length": "100000", "acmg_class": "4", "genes": []},
            {"sv": "chr1:1-10", "type": "DEL", "gt": "0/1", "length": "10", "acmg_class": "4", "genes": []}]"#,
        r#"[{"chrom": "chr22", "start": "50000000", "end": "51000000", "length_mb": "1.0", "type": "Dup", "cell_fraction": "0.1"},
            {"chrom": "chr22", "start": "52000000", "end": "53000000", "length_mb": "1.0", "type": "Del", "cell_fraction": "0.1"}]"#,
    );
    let ge = chart(&edge, "genome");
    assert_eq!(outside_marks(&ge.svg), Vec::<String>::new());
    let nd = ge
        .not_drawn
        .iter()
        .find(|x| x.what == "marks past the end of their GRCh38 chromosome");
    assert_eq!(nd.map(|x| x.count), Some(1));
    assert!(ge.marks.iter().all(|m| {
        m.end
            .is_none_or(|e| e <= chrom_len(m.chrom.as_deref().unwrap()).unwrap())
    }));
    assert_eq!(ge.svg.matches("class=\"sv\"").count(), 2);
    assert_eq!(ge.svg.matches("class=\"mca\"").count(), 1);
    let cut: Vec<(&str, i128)> = ge
        .marks
        .iter()
        .filter(|m| {
            matches!(m.kind, Some(MarkKind::Sv | MarkKind::Mca)) && m.pos.unwrap() > 40_000_000
        })
        .map(|m| (m.chrom.as_deref().unwrap(), m.end.unwrap()))
        .collect();
    assert_eq!(cut, [("chr21", 46_709_983), ("chr22", 50_818_468)]);
    // Red control: a span drawn past its chromosome's end.
    let i = ge.svg.find("<rect class=\"sv\" x=\"").unwrap() + "<rect class=\"sv\" x=\"".len();
    let j = i + ge.svg[i..].find('"').unwrap();
    let shifted = format!("{}2000{}", &ge.svg[..i], &ge.svg[j..]);
    assert!(!outside_marks(&shifted).is_empty());
}

#[test]
fn mito_dots_the_sites_and_marks_the_poly_c_call_apart() {
    let s = fixture();
    let mi = chart(&s, "mito");
    assert_eq!(mi.svg.matches("class=\"site\"").count(), 2);
    assert_eq!(mi.svg.matches("class=\"polyc\"").count(), 1);
    // A summary made before it listed the sites: none drawn, how many said.
    let mut old = s.clone();
    let Value::Object(d) = at(&mut old, &["sections", "mito", "data"]) else {
        panic!()
    };
    d.shift_remove("heteroplasmic_sites");
    let om = chart(&old, "mito");
    assert!(!om.svg.contains("class=\"site\""));
    assert_eq!(om.not_drawn[0].count, 2);
    assert!(
        om.svg
            .contains("positions of the heteroplasmic sites are not in this summary")
    );
    assert!(!om.svg.contains("No heteroplasmic site"));
}

#[test]
fn every_caveat_is_the_reports_own() {
    for sp in charts_of(&fixture()) {
        assert_eq!(Some(sp.caveat), notes::text(sp.caveat_key), "{}", sp.id);
    }
    let dir = repo().join("rust/crates/og-render/src/charts");
    let mut src = String::new();
    for e in std::fs::read_dir(&dir).unwrap() {
        src += &std::fs::read_to_string(e.unwrap().path()).unwrap();
    }
    assert_eq!(note_literals(&src), Vec::<String>::new());
    // Red control: a copied caveat is caught.
    let copied = format!("{src}\nconst X: &str = \"{}\";\n", notes::ROH_NOTE);
    assert_eq!(note_literals(&copied), ["ROH_NOTE"]);
}

#[test]
fn the_grch38_lengths_are_ucscs() {
    assert_eq!(GRCH38.len(), 24);
    assert_eq!(chrom_len("chr1"), Some(248_956_422));
    assert_eq!(chrom_len("chrY"), Some(57_227_415));
    assert_eq!(chrom_len("chrM"), None);
}

/// FNV-1a of a file's bytes: a frozen vector of what the oracle wrote.
fn fnv(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, b| {
        (h ^ u64::from(*b)).wrapping_mul(0x0100_0000_01b3)
    })
}

fn charts_file(v: &Value) -> String {
    let opts = Dump {
        indent: Some(1),
        ..Dump::default()
    };
    dumps(&json(&summary(v)).unwrap(), &opts).unwrap() + "\n"
}

/// The summaries whose charts JSON is frozen below. With OG_CHARTS_DUMP set
/// to a folder, each is written there, to render it with the oracle
/// (`bin/render_report.py --summary F --charts-json OUT`).
fn frozen() -> Vec<(&'static str, Value)> {
    let mut z = fixture();
    set(
        &mut z,
        &["sections", "clinvar", "data", "by_stars"],
        r#"{"4":0,"3":0,"2":0,"1":0,"0":0}"#,
    );
    set(&mut z, &["sections", "clinvar", "data", "hits"], "[]");
    set(
        &mut z,
        &["sections", "clinvar", "data", "low_confidence"],
        "0",
    );
    let mut na = fixture();
    set(&mut na, &["sections", "prs", "data", "adjusted"], "false");
    let mut old = fixture();
    if let Value::Object(d) = at(&mut old, &["sections", "mito", "data"]) {
        d.shift_remove("heteroplasmic_sites");
    }
    let edge = with_sv_mca(
        r#"[{"sv": "chr21:46700000-46800000", "type": "DUP", "gt": "0/1", "length": "100000", "acmg_class": "4", "genes": []},
            {"sv": "chr1:1-10", "type": "DEL", "gt": "0/1", "length": "10", "acmg_class": "4", "genes": []}]"#,
        r#"[{"chrom": "chr22", "start": "50000000", "end": "51000000", "length_mb": "1.0", "type": "Dup", "cell_fraction": "0.1"},
            {"chrom": "chr22", "start": "52000000", "end": "53000000", "length_mb": "1.0", "type": "Del", "cell_fraction": "0.1"}]"#,
    );
    let all = vec![
        ("fixture", fixture()),
        ("no-hit", z),
        ("no-panel", na),
        ("no-sites", old),
        ("sv-mca", with_sv_mca(SV_ROWS, MCA_CALLS)),
        ("edge", edge),
    ];
    if let Some(dir) = std::env::var_os("OG_CHARTS_DUMP") {
        for (name, v) in &all {
            let p = Path::new(&dir).join(format!("{name}.json"));
            std::fs::write(p, dumps(v, &Dump::default()).unwrap()).unwrap();
        }
    }
    all
}

/// The FNV-1a of the charts JSON `bin/report_charts.py` wrote at commit
/// 65095e97 for each frozen summary (render_report.py --charts-json). A
/// deliberate change to a chart changes its vector: say so in the PR.
const ORACLE: [(&str, u64); 6] = [
    ("fixture", 0x0fb8_e455_9f1d_5688),
    ("no-hit", 0x97e4_dae4_149c_3c39),
    ("no-panel", 0xd855_b5d4_785b_2c3a),
    ("no-sites", 0xa570_9eae_3ca9_1850),
    ("sv-mca", 0x003a_8f9c_2441_ba65),
    ("edge", 0xca4a_fec8_1f14_13d8),
];

#[test]
fn the_charts_json_is_what_the_python_wrote() {
    let got: Vec<(&str, u64)> = frozen()
        .iter()
        .map(|(n, v)| (*n, fnv(charts_file(v).as_bytes())))
        .collect();
    let shown: Vec<String> = got
        .iter()
        .map(|(n, h)| format!("(\"{n}\", {h:#018x})"))
        .collect();
    assert_eq!(got, ORACLE, "now: {}", shown.join(", "));
}

#[test]
fn the_charts_json_is_the_same_on_every_run() {
    let s = fixture();
    let (a, b) = (charts_file(&s), charts_file(&s));
    assert_eq!(a, b);
    let v = loads(&a).unwrap();
    let Value::Object(o) = &v else { panic!() };
    assert_eq!(o["schema_version"], Value::Int(1));
    let Value::Array(c) = &o["charts"] else {
        panic!()
    };
    let ids: Vec<String> = c
        .iter()
        .map(|x| {
            let Value::Object(x) = x else { panic!() };
            pycompat::value::py_str(&x["id"]).unwrap()
        })
        .collect();
    assert_eq!(
        ids,
        ["clinvar-stars", "genome", "mito", "ancestry-fine", "prs"]
    );
    // Red control: one changed value is a different file.
    let mut moved = s.clone();
    set(
        &mut moved,
        &["sections", "mito", "data", "heteroplasmic_sites", "0", "af"],
        "0.31",
    );
    assert_ne!(charts_file(&moved), a);
}

#[test]
fn the_verb_writes_the_file_and_counts_the_sections() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("charts-verb");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let inp = dir.join("s.json");
    std::fs::write(&inp, dumps(&fixture(), &Dump::default()).unwrap()).unwrap();
    let outp = dir.join("S_charts.json");
    let args: Vec<String> = [
        "--summary",
        inp.to_str().unwrap(),
        "--charts-json",
        outp.to_str().unwrap(),
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(
        charts::cli(&args, &mut out, &mut err),
        0,
        "{}",
        String::from_utf8_lossy(&err)
    );
    assert_eq!(String::from_utf8(out).unwrap(), "Sections: ok 5\n");
    assert_eq!(
        std::fs::read_to_string(&outp).unwrap(),
        charts_file(&fixture())
    );
    // No summary is a usage error; an unreadable one fails.
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(
        charts::cli(
            &["--charts-json".into(), "x.json".into()],
            &mut out,
            &mut err
        ),
        2
    );
    let missing = dir.join("none.json").to_str().unwrap().to_owned();
    assert_eq!(
        charts::cli(&["--summary".into(), missing], &mut out, &mut err),
        1
    );
}

#[test]
fn both_report_paths_write_the_charts_beside_the_report() {
    let nf = std::fs::read_to_string(repo().join("modules/local/html_report/main.nf")).unwrap();
    let sh = std::fs::read_to_string(repo().join("scripts/24-html-report.sh")).unwrap();
    assert!(nf.contains("--charts-json ${id}_charts.json"));
    assert!(nf.contains("path(\"${meta.id}_charts.json\")"));
    assert!(
        nf.split("stub:")
            .nth(1)
            .unwrap()
            .contains("${meta.id}_charts.json")
    );
    assert!(sh.contains("--charts-json \"/genome/${SAMPLE}/${SAMPLE}_charts.json\""));
}
