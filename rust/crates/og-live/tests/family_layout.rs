//! The pedigree drawing (`og_pedigree::layout`, `og live family svg`) on
//! tests/fixtures/families_live, read back from its data- attributes and
//! coordinates; moved here from tests/test_pedigree_svg.py, check for check:
//!
//! 1. family_a (members.tsv, three generations): every person is drawn
//!    once, with the sex and sequenced state of the file; each generation is
//!    one row, parents above their children;
//! 2. each line of descent starts at its two parents (the middle of their
//!    partner line) and ends on top of each child it names, and the parents
//!    it names are the ones this test expects (written here, not read from
//!    the fixture);
//! 3. family_b (pedigree.tsv only): the couple is drawn joined, both filled;
//! 4. a person whose sex is not known is a diamond; the command line writes
//!    the same SVG and exits 2 on a folder without a pedigree.
//!
//! Then the family page's geometry (`page_layout`, family.json's `layout`)
//! is the drawing's, number for number.
//!
//! Checks 1 and 2 run again on a copy of the fixture with two parents
//! swapped and must fail there: a check that cannot fail is no check.

use og_pedigree::layout::{self, SIZE, load, page_layout};
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn fix(name: &str) -> PathBuf {
    repo().join("tests/fixtures/families_live").join(name)
}

fn tmp(name: &str) -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "family-layout-{name}-{}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// What family_a declares, written out by hand: child, father, mother.
const EXPECTED: [(&str, &str, &str); 4] = [
    ("child_a", "father_a", "mother_a"),
    ("child_b", "father_a", "mother_a"),
    ("father_a", "gf_a", "gm_a"),
    ("mother_a", ".", "gm_m"),
];
const SEX: [(&str, &str); 7] = [
    ("child_a", "female"),
    ("child_b", "male"),
    ("father_a", "male"),
    ("mother_a", "female"),
    ("gf_a", "male"),
    ("gm_a", "female"),
    ("gm_m", "female"),
];

#[derive(Debug, Clone)]
struct Drawn {
    sex: String,
    generation: i64,
    seq: String,
    shape: String,
    x: f64,
    y: f64,
}

#[derive(Debug, Clone)]
struct Descent {
    father: String,
    mother: String,
    children: Vec<String>,
    lines: Vec<[f64; 4]>,
}

struct Svg {
    people: BTreeMap<String, Vec<Drawn>>,
    descents: Vec<Descent>,
    couples: Vec<(String, String, [f64; 4])>,
}

fn num(s: &str) -> f64 {
    s.parse().unwrap()
}

fn parse(svg: &str) -> Svg {
    let person = Regex::new(
        r#"<g class="person" data-sample="([^"]*)" data-sex="([^"]*)" data-gen="(\d+)" data-sequenced="([a-z]+)">(?:<rect x="([-\d.]+)" y="([-\d.]+)"|<circle cx="([-\d.]+)" cy="([-\d.]+)"|<polygon points="([^"]*)")"#,
    )
    .unwrap();
    let mut people: BTreeMap<String, Vec<Drawn>> = BTreeMap::new();
    for c in person.captures_iter(svg) {
        let h = SIZE as f64 / 2.0;
        let (shape, x, y) = if let Some(x) = c.get(5) {
            ("rect", num(x.as_str()) + h, num(&c[6]) + h)
        } else if let Some(x) = c.get(7) {
            ("circle", num(x.as_str()), num(&c[8]))
        } else {
            let pts: Vec<(f64, f64)> = c[9]
                .split(' ')
                .map(|p| {
                    let (a, b) = p.split_once(',').unwrap();
                    (num(a), num(b))
                })
                .collect();
            ("polygon", pts[0].0, (pts[0].1 + pts[2].1) / 2.0)
        };
        people.entry(c[1].to_owned()).or_default().push(Drawn {
            sex: c[2].to_owned(),
            generation: c[3].parse().unwrap(),
            seq: c[4].to_owned(),
            shape: shape.to_owned(),
            x,
            y,
        });
    }
    let line = Regex::new(r#"<line x1="([-\d.]+)" y1="([-\d.]+)" x2="([-\d.]+)" y2="([-\d.]+)"/>"#)
        .unwrap();
    let descent =
        Regex::new(r#"<g class="descent" data-father="([^"]*)" data-mother="([^"]*)" data-children="([^"]*)"[^>]*>(.*?)</g>"#)
            .unwrap();
    let descents = descent
        .captures_iter(svg)
        .map(|c| Descent {
            father: c[1].to_owned(),
            mother: c[2].to_owned(),
            children: c[3].split(' ').map(str::to_owned).collect(),
            lines: line
                .captures_iter(&c[4])
                .map(|l| [num(&l[1]), num(&l[2]), num(&l[3]), num(&l[4])])
                .collect(),
        })
        .collect();
    let couple = Regex::new(
        r#"<line class="couple" data-a="([^"]*)" data-b="([^"]*)" x1="([-\d.]+)" y1="([-\d.]+)" x2="([-\d.]+)" y2="([-\d.]+)""#,
    )
    .unwrap();
    let couples = couple
        .captures_iter(svg)
        .map(|c| {
            (
                c[1].to_owned(),
                c[2].to_owned(),
                [num(&c[3]), num(&c[4]), num(&c[5]), num(&c[6])],
            )
        })
        .collect();
    Svg {
        people,
        descents,
        couples,
    }
}

fn people_problems(svg: &Svg) -> Vec<String> {
    let mut out = Vec::new();
    for (s, sex) in SEX {
        let want_seq = if s == "gm_m" { "no" } else { "yes" };
        let want_shape = if sex == "male" { "rect" } else { "circle" };
        match svg.people.get(s).map(Vec::as_slice) {
            Some([d]) if d.sex == sex && d.seq == want_seq && d.shape == want_shape => {}
            got => out.push(format!("{s}: {got:?}")),
        }
    }
    let names: BTreeSet<&str> = svg.people.keys().map(String::as_str).collect();
    if names != SEX.iter().map(|(s, _)| *s).collect() {
        out.push(format!("people {names:?}"));
    }
    out
}

fn row_problems(svg: &Svg) -> Vec<String> {
    let mut out = Vec::new();
    let y = |s: &str| svg.people[s][0].y;
    for (child, f, m) in EXPECTED {
        for p in [f, m] {
            if p != "." && y(p) >= y(child) {
                out.push(format!("{p} not above {child}"));
            }
        }
    }
    let mut ys: BTreeMap<i64, BTreeSet<String>> = BTreeMap::new();
    for v in svg.people.values() {
        ys.entry(v[0].generation)
            .or_default()
            .insert(format!("{}", v[0].y));
    }
    if ys.values().any(|y| y.len() != 1) {
        out.push(format!("rows {ys:?}"));
    }
    out
}

fn descent_problems(svg: &Svg) -> Vec<String> {
    let mut out = Vec::new();
    let mut seen: BTreeMap<String, (String, String)> = BTreeMap::new();
    for d in &svg.descents {
        for c in &d.children {
            seen.insert(c.clone(), (d.father.clone(), d.mother.clone()));
        }
        let ps: Vec<&Drawn> = [&d.father, &d.mother]
            .into_iter()
            .filter(|p| *p != ".")
            .map(|p| &svg.people[p][0])
            .collect();
        let ox = ps.iter().map(|p| p.x).sum::<f64>() / ps.len() as f64;
        let first = d.lines[0];
        if (first[0] - ox).abs() > 0.5 || (first[2] - ox).abs() > 0.5 {
            out.push(format!("start of {} {}", d.father, d.mother));
        }
        for c in &d.children {
            let p = &svg.people[c][0];
            let h = SIZE as f64 / 2.0;
            if !d
                .lines
                .iter()
                .any(|l| (l[0] - p.x).abs() < 0.5 && (l[3] - (p.y - h)).abs() < 0.5)
            {
                out.push(format!("end at {c}"));
            }
        }
    }
    let want: BTreeMap<String, (String, String)> = EXPECTED
        .iter()
        .map(|(c, f, m)| ((*c).into(), ((*f).into(), (*m).into())))
        .collect();
    if seen != want {
        out.push(format!("parents {seen:?}"));
    }
    out
}

fn drawing(d: &Path, title: &str) -> String {
    let (people, _) = load(d).unwrap();
    layout::draw(&people, title).unwrap()
}

fn cli(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_live::family::cli(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

#[test]
fn family_a_is_drawn_as_declared() {
    let svg = parse(&drawing(&fix("family_a"), "family_a"));
    assert_eq!(
        people_problems(&svg),
        Vec::<String>::new(),
        "every person drawn once with the file's sex and sequenced state"
    );
    assert_eq!(
        row_problems(&svg),
        Vec::<String>::new(),
        "one row per generation, parents above their children"
    );
    assert_eq!(
        descent_problems(&svg),
        Vec::<String>::new(),
        "each line of descent from the expected parents to each child"
    );
    let couples: BTreeSet<BTreeSet<String>> = svg
        .couples
        .iter()
        .map(|(a, b, _)| [a.clone(), b.clone()].into_iter().collect())
        .collect();
    let want: BTreeSet<BTreeSet<String>> = [["father_a", "mother_a"], ["gf_a", "gm_a"]]
        .iter()
        .map(|p| p.iter().map(|s| (*s).to_owned()).collect())
        .collect();
    assert_eq!(
        couples, want,
        "partners joined: father_a with mother_a, gf_a with gm_a"
    );
}

#[test]
fn control_child_drawn_under_its_grandparents() {
    // child_a's parents swapped for father_a's, and father_a left without parents
    let bad = tmp("swapped");
    let text = fs::read_to_string(fix("family_a").join("members.tsv")).unwrap();
    let lines: Vec<String> = text
        .split('\n')
        .map(|ln| {
            let mut c: Vec<&str> = ln.split('\t').collect();
            match c.first().copied() {
                Some("child_a") => {
                    c[2] = "gf_a";
                    c[3] = "gm_a";
                }
                Some("father_a") => {
                    c[2] = "";
                    c[3] = "";
                }
                _ => {}
            }
            c.join("\t")
        })
        .collect();
    fs::write(bad.join("members.tsv"), lines.join("\n")).unwrap();
    let svg = parse(&drawing(&bad, "family_a"));
    assert!(
        !descent_problems(&svg).is_empty() || !row_problems(&svg).is_empty(),
        "control: child_a drawn under its grandparents was NOT caught"
    );
}

#[test]
fn family_b_from_its_pedigree_is_a_joined_filled_couple() {
    let svg = parse(&drawing(&fix("family_b"), "family_b"));
    let couples: Vec<(&str, &str)> = svg
        .couples
        .iter()
        .map(|(a, b, _)| (a.as_str(), b.as_str()))
        .collect();
    assert_eq!(couples, [("partner_x", "partner_y")]);
    assert!(svg.descents.is_empty());
    assert!(
        ["partner_x", "partner_y"]
            .iter()
            .all(|s| svg.people[*s][0].seq == "yes"),
        "{:?}",
        svg.people
    );
}

#[test]
fn unknown_sex_is_a_diamond_and_named_parents_are_not_sequenced() {
    let d = tmp("unknown");
    fs::write(
        d.join("pedigree.tsv"),
        "sample\tsex\tfather\tmother\tpartner\nkid_u\tunknown\tdad_u\tmum_u\t.\n",
    )
    .unwrap();
    let svg = parse(&drawing(&d, "Pedigree"));
    assert_eq!(svg.people["kid_u"][0].shape, "polygon");
    assert_eq!(svg.people["dad_u"][0].seq, "no");
    assert_eq!(svg.descents.len(), 1);
    assert_eq!(svg.descents[0].children, ["kid_u"]);
}

#[test]
fn the_command_line_draws_the_same_and_refuses_an_empty_folder() {
    let a = fix("family_a");
    let (code, out, err) = cli(&["svg", a.to_str().unwrap()]);
    assert_eq!(
        (code, out.as_str()),
        (0, drawing(&a, "family_a").as_str()),
        "{err}"
    );
    let (code, _, err) = cli(&["svg", a.to_str().unwrap(), "--title", "x<y"]);
    assert!(code == 0 && err.is_empty());
    let empty = tmp("empty");
    let (code, out, err) = cli(&["svg", empty.to_str().unwrap()]);
    assert!(
        code == 2 && out.is_empty() && err.contains("members.tsv"),
        "{code} {err}"
    );
    // a header and no one: nothing to draw
    fs::write(
        empty.join("members.tsv"),
        "sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced\n",
    )
    .unwrap();
    let (code, _, err) = cli(&["svg", empty.to_str().unwrap()]);
    assert!(code == 2 && err.contains("names no one"), "{code} {err}");
}

#[test]
fn the_title_is_escaped() {
    let svg = drawing(&fix("family_b"), "a<b & \"c\"");
    assert!(svg.contains(">a&lt;b &amp; &quot;c&quot;</text>"), "{svg}");
}

/// family.json's layout and the SVG agree number for number: each
/// person's centre, every partner line and every line of descent.
fn geometry_problems(d: &Path) -> Vec<String> {
    let (people, _) = load(d).unwrap();
    let svg = parse(&layout::draw(&people, "t").unwrap());
    let page = page_layout(&people);
    let mut out = Vec::new();
    for p in &page.people {
        let drawn = &svg.people[&p.sample][0];
        if (drawn.x, drawn.y, drawn.generation) != (p.x, p.y, p.generation) || drawn.sex != p.sex {
            out.push(format!(
                "{}: page ({}, {}) svg ({}, {})",
                p.sample, p.x, p.y, drawn.x, drawn.y
            ));
        }
    }
    let mut lines: Vec<[f64; 4]> = svg.couples.iter().map(|(_, _, l)| *l).collect();
    lines.extend(svg.descents.iter().flat_map(|d| d.lines.clone()));
    if lines != page.lines {
        out.push(format!("lines: page {:?} svg {:?}", page.lines, lines));
    }
    if page.size != SIZE || page.people.len() != svg.people.len() {
        out.push(format!("size {} people {}", page.size, page.people.len()));
    }
    out
}

#[test]
fn the_page_geometry_is_the_drawings() {
    for d in [
        fix("family_a"),
        fix("family_b"),
        repo().join("tests/fixtures/live/demo/DEMO-FAMILY"),
    ] {
        assert_eq!(
            geometry_problems(&d),
            Vec::<String>::new(),
            "{}",
            d.display()
        );
    }
    let (people, _) = load(&fix("family_a")).unwrap();
    let page = page_layout(&people);
    // three rows of 2, 2 and 2 or 3 people at most: 2 * 40 + 36 + 2 * 120 wide, 50 + 36 + 2 * 130 + 30 high
    assert_eq!((page.width, page.height), (356, 376));
}

#[test]
fn control_a_moved_symbol_breaks_the_geometry_check() {
    // the same folder, one partner line drawn from the other partner: the page and the drawing then differ
    let d = tmp("moved");
    fs::copy(fix("family_a").join("members.tsv"), d.join("members.tsv")).unwrap();
    let (people, _) = load(&d).unwrap();
    let svg = layout::draw(&people, "t")
        .unwrap()
        .replacen(r#"x1="136.0""#, r#"x1="137.0""#, 1);
    let parsed = parse(&svg);
    let page = page_layout(&people);
    let first = parsed.couples[0].2;
    assert_ne!(
        first, page.lines[0],
        "control: a moved partner line was NOT caught"
    );
}
