//! The drawing of a pedigree: where each person goes, the lines that join
//! partners and drop to their children, and the SVG of it all. The family
//! page draws the same geometry itself, so each person is a link
//! ([`page_layout`]); `og live family svg` writes the picture ([`draw`]).
//!
//! It only draws what the file says. It reads no genotype and checks no
//! relationship: step 38's relationship check does that.
//!
//! The file is `members.tsv` (step 38's table of every person the pedigree
//! names, with `in_pedigree` and `sequenced`) or, when there is none,
//! `pedigree.tsv` (the rows a person wrote; parents and partners named
//! without a row are added, as step 38 adds them, by [`crate::build_members`]).
//!
//! The usual pedigree symbols: a square for a male, a circle for a female,
//! a diamond when the sex is not known. A filled symbol is a sequenced
//! person (`members.tsv`'s `sequenced` column; from `pedigree.tsv`, a person
//! with a row of their own). A horizontal line joins two partners, or the
//! two parents of a child; a line drops from it to a bar over their
//! children. Each generation is one row, parents above their children.
//!
//! Every symbol is a `<g class="person" data-sample data-sex data-gen
//! data-sequenced>`, every partner line a `<line class="couple" data-a
//! data-b>`, and every line of descent a `<g class="descent" data-father
//! data-mother data-children>`, so a test can read the drawing back.

use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use indexmap::IndexMap;
use pycompat::text::{lower, strip};

/// A symbol's width and height.
pub const SIZE: i64 = 36;
/// From one person to the next in a row.
pub const STEP_X: i64 = 120;
/// From one generation to the next.
pub const STEP_Y: i64 = 130;
pub const MARGIN: i64 = 40;
/// Room for the title.
pub const TOP: i64 = 50;
const BG: &str = "#16181c";
const FG: &str = "#e6e6e6";
const LINE: &str = "#c0c4cc";
const FILL: &str = "#7aa2f7";
/// The files of a family folder that declare it, in the order they are read.
pub const MARKERS: [&str; 2] = ["members.tsv", "pedigree.tsv"];

/// One person as the drawing needs them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Person {
    /// `male`, `female` or `unknown`.
    pub sex: String,
    pub father: Option<String>,
    pub mother: Option<String>,
    pub partner: Option<String>,
    pub sequenced: bool,
}

impl Person {
    fn rel(&self, k: &str) -> Option<&str> {
        match k {
            "father" => self.father.as_deref(),
            "mother" => self.mother.as_deref(),
            _ => self.partner.as_deref(),
        }
    }
}

/// Everyone the file names, in its order.
pub type People = IndexMap<String, Person>;

/// A value the file leaves empty.
fn blank(v: Option<&str>) -> bool {
    v.is_none_or(|v| matches!(strip(v, None), "" | "." | "-9" | "0" | "NA" | "na"))
}

fn sex_of(v: &str) -> String {
    let s = lower(strip(v, None));
    if s == "male" || s == "female" {
        s
    } else {
        "unknown".into()
    }
}

/// The header of a tab-separated file: its first line, a leading `#`
/// dropped.
fn header(line: &str) -> Vec<String> {
    line.trim_start_matches('#')
        .split('\t')
        .map(str::to_owned)
        .collect()
}

/// A line's values by column, short lines padded with empty values.
fn record<'a>(head: &'a [String], line: &'a str) -> HashMap<&'a str, &'a str> {
    let vals: Vec<&str> = line.split('\t').collect();
    let mut out = HashMap::new();
    for (i, h) in head.iter().enumerate() {
        out.insert(h.as_str(), vals.get(i).copied().unwrap_or(""));
    }
    out
}

fn lines(path: &Path) -> Result<Vec<String>, String> {
    let p = path.to_string_lossy();
    crate::text_lines(&p).map_err(|d| d.lines.join("; "))
}

/// The people of `members.tsv`: a row with an empty name or a name already
/// seen is skipped; a sex other than male or female is `unknown`.
fn read_members_tsv(path: &Path, text: &[String]) -> Result<People, String> {
    let head = header(text.first().map_or("", String::as_str));
    if !head.iter().any(|h| h == "sample") {
        return Err(format!(
            "{}: the first line must be a header with a 'sample' column",
            path.display()
        ));
    }
    let mut out = People::new();
    for line in text.iter().skip(1) {
        if strip(line, None).is_empty() || line.starts_with('#') {
            continue;
        }
        let r = record(&head, line);
        let s = strip(r["sample"], None);
        if blank(Some(s)) || out.contains_key(s) {
            continue;
        }
        let name = |k: &str| {
            let v = r.get(k).copied();
            if blank(v) {
                None
            } else {
                v.map(|v| strip(v, None).to_owned())
            }
        };
        out.insert(
            s.to_owned(),
            Person {
                sex: sex_of(r.get("sex").copied().unwrap_or("")),
                father: name("father"),
                mother: name("mother"),
                partner: name("partner"),
                sequenced: strip(r.get("sequenced").copied().unwrap_or(""), None) == "yes",
            },
        );
    }
    Ok(out)
}

/// The people of a pedigree file, read by step 38's own rules; any problem
/// with the file is an error naming them all.
fn read_pedigree_tsv(path: &Path) -> Result<People, String> {
    let p = path.to_string_lossy();
    let (rows, errors) = crate::read_pedigree(&p).map_err(|d| d.lines.join("; "))?;
    if !errors.is_empty() {
        return Err(errors.join("; "));
    }
    let (members, more) = crate::build_members(&rows, &HashMap::new(), &p, "the samplesheet");
    if !more.is_empty() {
        return Err(more.join("; "));
    }
    Ok(members
        .into_iter()
        .map(|(s, m)| {
            let sex = match m.sex.as_deref() {
                Some(x @ ("male" | "female")) => x.to_owned(),
                _ => "unknown".into(),
            };
            (
                s,
                Person {
                    sex,
                    father: m.father,
                    mother: m.mother,
                    partner: m.partner,
                    sequenced: m.row,
                },
            )
        })
        .collect())
}

/// The people of a family folder (its `members.tsv`, else its
/// `pedigree.tsv`) or of one such file, and the file read. A file whose
/// header has a `sequenced` column is read as `members.tsv`.
pub fn load(path: &Path) -> Result<(People, PathBuf), String> {
    if path.is_dir() {
        for name in MARKERS {
            let p = path.join(name);
            if p.is_file() {
                return load(&p);
            }
        }
        return Err(format!(
            "{} holds no members.tsv or pedigree.tsv",
            path.display()
        ));
    }
    let text = lines(path)?;
    let head = header(text.first().map_or("", String::as_str));
    let people = if head.iter().any(|h| h == "sequenced") {
        read_members_tsv(path, &text)?
    } else {
        read_pedigree_tsv(path)?
    };
    Ok((people, path.to_path_buf()))
}

/// The pairs drawn joined: partners, and the two parents of a child, in
/// file order, each pair once.
pub fn couples(people: &People) -> Vec<(String, String)> {
    let mut seen: HashSet<(String, String)> = HashSet::new();
    let mut out = Vec::new();
    for (s, m) in people {
        for (a, b) in [
            (Some(s.as_str()), m.rel("partner")),
            (m.rel("father"), m.rel("mother")),
        ] {
            let (Some(a), Some(b)) = (a, b) else { continue };
            if !people.contains_key(a) || !people.contains_key(b) {
                continue;
            }
            let key = if a <= b {
                (a.to_owned(), b.to_owned())
            } else {
                (b.to_owned(), a.to_owned())
            };
            if seen.insert(key) {
                out.push((a.to_owned(), b.to_owned()));
            }
        }
    }
    out
}

/// The generation of each person: 0 with no parent named, else one below
/// the lower parent. A partner with no parent named moves down to the
/// other partner's row; two partners who both have parents keep their own
/// rows (a line then joins them across rows).
pub fn generations(people: &People, pairs: &[(String, String)]) -> HashMap<String, i64> {
    let mut generation: HashMap<String, i64> = people.keys().map(|s| (s.clone(), 0)).collect();
    let parents = |s: &str| -> Vec<String> {
        let m = &people[s];
        [&m.father, &m.mother]
            .into_iter()
            .flatten()
            .filter(|p| people.contains_key(p.as_str()))
            .cloned()
            .collect()
    };
    for _ in 0..people.len() + 2 {
        let mut changed = false;
        for s in people.keys() {
            let g = parents(s)
                .iter()
                .map(|p| generation[p] + 1)
                .chain([generation[s]])
                .max()
                .unwrap_or(0);
            if g != generation[s] {
                generation.insert(s.clone(), g);
                changed = true;
            }
        }
        for (a, b) in pairs {
            for (x, y) in [(a, b), (b, a)] {
                if parents(x).is_empty() && generation[x] < generation[y] {
                    let g = generation[y];
                    generation.insert(x.clone(), g);
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    generation
}

/// Where everyone goes.
#[derive(Debug, Clone, PartialEq)]
pub struct Layout {
    /// The centre of each person's symbol.
    pub pos: HashMap<String, (f64, f64)>,
    pub generation: HashMap<String, i64>,
    /// The pairs drawn joined ([`couples`]).
    pub pairs: Vec<(String, String)>,
    /// People in the longest row.
    pub widest: i64,
    /// Rows: the last generation plus one.
    pub ngen: i64,
}

/// The centre of each symbol. People with a parent drawn go under their
/// parents (by the mean place of the parents in their row); one without (a
/// partner who married in, or a whole first generation) goes beside the
/// partner already placed, else in file order. Each row is centred on the
/// widest.
pub fn layout(people: &People) -> Layout {
    let pairs = couples(people);
    let generation = generations(people, &pairs);
    let mut partners: HashMap<&str, Vec<&str>> = HashMap::new();
    for (a, b) in &pairs {
        partners.entry(a).or_default().push(b);
        partners.entry(b).or_default().push(a);
    }
    let rank: HashMap<&str, usize> = people
        .keys()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut gens: Vec<i64> = generation
        .values()
        .copied()
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    gens.sort_unstable();
    let mut rows: Vec<(i64, Vec<String>)> = Vec::new();
    let mut x_of: HashMap<String, usize> = HashMap::new();
    for &g in &gens {
        let row_people: Vec<&str> = people
            .keys()
            .filter(|s| generation[*s] == g)
            .map(String::as_str)
            .collect();
        let key = |s: &str, x_of: &HashMap<String, usize>| -> Option<f64> {
            let m = &people[s];
            let ps: Vec<f64> = [&m.father, &m.mother]
                .into_iter()
                .flatten()
                .filter_map(|p| x_of.get(p).map(|&x| x as f64))
                .collect();
            (!ps.is_empty()).then(|| ps.iter().sum::<f64>() / ps.len() as f64)
        };
        let mut anchored: Vec<(f64, usize, &str)> = row_people
            .iter()
            .filter_map(|&s| key(s, &x_of).map(|k| (k, rank[s], s)))
            .collect();
        anchored.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut row: Vec<String> = Vec::new();
        let mut placed: HashSet<&str> = HashSet::new();
        fn put<'a>(
            s: &'a str,
            g: i64,
            row: &mut Vec<String>,
            placed: &mut HashSet<&'a str>,
            partners: &HashMap<&'a str, Vec<&'a str>>,
            unanchored: &dyn Fn(&str) -> bool,
            generation: &HashMap<String, i64>,
        ) {
            if !placed.insert(s) {
                return;
            }
            row.push(s.to_owned());
            for &p in partners.get(s).map_or(&[][..], Vec::as_slice) {
                if generation.get(p) == Some(&g) && unanchored(p) {
                    put(p, g, row, placed, partners, unanchored, generation);
                }
            }
        }
        let unanchored = |p: &str| key(p, &x_of).is_none();
        let order: Vec<&str> = anchored
            .iter()
            .map(|a| a.2)
            .chain(
                row_people
                    .iter()
                    .copied()
                    .filter(|&s| key(s, &x_of).is_none()),
            )
            .collect();
        for s in order {
            put(
                s,
                g,
                &mut row,
                &mut placed,
                &partners,
                &unanchored,
                &generation,
            );
        }
        for (i, s) in row.iter().enumerate() {
            x_of.insert(s.clone(), i);
        }
        rows.push((g, row));
    }
    let widest = rows.iter().map(|(_, r)| r.len()).max().unwrap_or(0) as i64;
    let mut pos = HashMap::new();
    for (g, row) in &rows {
        let shift = (widest - row.len() as i64) as f64 * STEP_X as f64 / 2.0;
        for (i, s) in row.iter().enumerate() {
            pos.insert(
                s.clone(),
                (
                    MARGIN as f64 + SIZE as f64 / 2.0 + shift + (i as i64 * STEP_X) as f64,
                    TOP as f64 + SIZE as f64 / 2.0 + (g * STEP_Y) as f64,
                ),
            );
        }
    }
    let ngen = rows.last().map_or(0, |(g, _)| g + 1);
    Layout {
        pos,
        generation,
        pairs,
        widest,
        ngen,
    }
}

/// The children of each pair of parents drawn, in file order: (father,
/// mother, children), a parent who is not drawn left out.
pub fn sibships(
    people: &People,
    pos: &HashMap<String, (f64, f64)>,
) -> Vec<(Option<String>, Option<String>, Vec<String>)> {
    let mut out: Vec<(Option<String>, Option<String>, Vec<String>)> = Vec::new();
    for (s, m) in people {
        let drawn = |p: &Option<String>| p.clone().filter(|p| pos.contains_key(p));
        let (f, mo) = (drawn(&m.father), drawn(&m.mother));
        if f.is_none() && mo.is_none() {
            continue;
        }
        match out.iter_mut().find(|(a, b, _)| *a == f && *b == mo) {
            Some(e) => e.2.push(s.clone()),
            None => out.push((f, mo, vec![s.clone()])),
        }
    }
    out
}

/// One sibship's lines: from the parents (the middle of their partner
/// line, or under a lone parent) down to a bar over the children, the bar,
/// and one line down to each child. Each line is `[x1, y1, x2, y2]`.
fn descent_lines(
    f: Option<&str>,
    mo: Option<&str>,
    kids: &[String],
    pos: &HashMap<String, (f64, f64)>,
) -> Vec<[f64; 4]> {
    let h = SIZE as f64 / 2.0;
    let (ox, oy) = match (f, mo) {
        (Some(f), Some(mo)) => ((pos[f].0 + pos[mo].0) / 2.0, (pos[f].1 + pos[mo].1) / 2.0),
        _ => {
            let p = pos[f.or(mo).unwrap_or_default()];
            (p.0, p.1 + h)
        }
    };
    let bar = kids.iter().map(|k| pos[k].1).fold(f64::INFINITY, f64::min) - h - 22.0;
    let xs: Vec<f64> = kids.iter().map(|k| pos[k].0).chain([ox]).collect();
    let (lo, hi) = (
        xs.iter().copied().fold(f64::INFINITY, f64::min),
        xs.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    );
    let mut out = vec![[ox, oy, ox, bar], [lo, bar, hi, bar]];
    out.extend(kids.iter().map(|k| [pos[k].0, bar, pos[k].0, pos[k].1 - h]));
    out
}

/// One person of [`PageLayout`].
#[derive(Debug, Clone, PartialEq)]
pub struct PagePerson {
    pub sample: String,
    pub x: f64,
    pub y: f64,
    pub generation: i64,
    pub sex: String,
    pub sequenced: bool,
}

/// The drawing's geometry as data, for the family page, which draws each
/// person as a link: the people with their centre, then every line
/// (partner lines first, then the descent lines of each sibship), each
/// number rounded to one decimal.
#[derive(Debug, Clone, PartialEq)]
pub struct PageLayout {
    pub width: i64,
    pub height: i64,
    pub size: i64,
    pub people: Vec<PagePerson>,
    pub lines: Vec<[f64; 4]>,
}

fn round1(v: f64) -> f64 {
    pycompat::fmt::round(v, 1).unwrap_or(v)
}

/// The geometry [`draw`] draws, as the family page draws it.
pub fn page_layout(people: &People) -> PageLayout {
    let l = layout(people);
    let h = SIZE as f64 / 2.0;
    let persons = people
        .iter()
        .map(|(s, m)| PagePerson {
            sample: s.clone(),
            x: round1(l.pos[s].0),
            y: round1(l.pos[s].1),
            generation: l.generation[s],
            sex: m.sex.clone(),
            sequenced: m.sequenced,
        })
        .collect();
    let mut lines = Vec::new();
    for (a, b) in &l.pairs {
        let (pa, pb) = (l.pos[a], l.pos[b]);
        // the two ends in order of x, then y
        let ((xa, ya), (xb, yb)) = if (pa.0, pa.1) <= (pb.0, pb.1) {
            (pa, pb)
        } else {
            (pb, pa)
        };
        lines.push([xa + h, ya, xb - h, yb]);
    }
    for (f, mo, kids) in sibships(people, &l.pos) {
        lines.extend(descent_lines(f.as_deref(), mo.as_deref(), &kids, &l.pos));
    }
    PageLayout {
        width: 2 * MARGIN + SIZE + (l.widest - 1) * STEP_X,
        height: TOP + SIZE + (l.ngen - 1) * STEP_Y + 30,
        size: SIZE,
        people: persons,
        lines: lines.into_iter().map(|ln| ln.map(round1)).collect(),
    }
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn symbol(s: &str, m: &Person, x: f64, y: f64, g: i64) -> String {
    let h = SIZE as f64 / 2.0;
    let fill = if m.sequenced { FILL } else { BG };
    let shape = match m.sex.as_str() {
        "male" => format!(
            r#"<rect x="{:.1}" y="{:.1}" width="{SIZE}" height="{SIZE}""#,
            x - h,
            y - h
        ),
        "female" => format!(r#"<circle cx="{x:.1}" cy="{y:.1}" r="{h:.1}""#),
        _ => format!(
            r#"<polygon points="{:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1}""#,
            x,
            y - h,
            x + h,
            y,
            x,
            y + h,
            x - h,
            y
        ),
    };
    let e = esc(s);
    format!(
        r#"<g class="person" data-sample="{e}" data-sex="{sex}" data-gen="{g}" data-sequenced="{seq}">{shape} fill="{fill}" stroke="{LINE}" stroke-width="2"/><text x="{x:.1}" y="{ty:.1}" fill="{FG}" text-anchor="middle">{e}</text><title>{e}: {sex}, {said}</title></g>"#,
        sex = m.sex,
        seq = if m.sequenced { "yes" } else { "no" },
        ty = y + h + 16.0,
        said = if m.sequenced {
            "sequenced"
        } else {
            "not sequenced"
        },
    )
}

/// The pedigree as an SVG, titled `title`; an error when the file names
/// no one.
pub fn draw(people: &People, title: &str) -> Result<String, String> {
    if people.is_empty() {
        return Err("no person to draw".into());
    }
    let l = layout(people);
    let width = 2 * MARGIN + SIZE + (l.widest - 1) * STEP_X;
    let legend_y = TOP + SIZE + (l.ngen - 1) * STEP_Y + 50;
    let height = legend_y + 30;
    let mut out: Vec<String> = vec![
        format!(
            r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="sans-serif" font-size="12">"#
        ),
        format!(r#"<rect width="100%" height="100%" fill="{BG}"/>"#),
        format!(
            r#"<text x="10" y="22" fill="{FG}" font-size="14">{}</text>"#,
            esc(title)
        ),
    ];
    let h = SIZE as f64 / 2.0;
    for (a, b) in &l.pairs {
        let (mut pa, mut pb) = (l.pos[a], l.pos[b]);
        if pa.0 > pb.0 {
            std::mem::swap(&mut pa, &mut pb);
        }
        out.push(format!(
            r#"<line class="couple" data-a="{}" data-b="{}" x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{LINE}" stroke-width="2"/>"#,
            esc(a),
            esc(b),
            pa.0 + h,
            pa.1,
            pb.0 - h,
            pb.1
        ));
    }
    for (f, mo, kids) in sibships(people, &l.pos) {
        let mut parts = String::new();
        for [x1, y1, x2, y2] in descent_lines(f.as_deref(), mo.as_deref(), &kids, &l.pos) {
            let _ = write!(
                parts,
                r#"<line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}"/>"#
            );
        }
        out.push(format!(
            r#"<g class="descent" data-father="{}" data-mother="{}" data-children="{}" stroke="{LINE}" stroke-width="2">{parts}</g>"#,
            esc(f.as_deref().unwrap_or(".")),
            esc(mo.as_deref().unwrap_or(".")),
            esc(&kids.join(" "))
        ));
    }
    for (s, m) in people {
        let (x, y) = l.pos[s];
        out.push(symbol(s, m, x, y, l.generation[s]));
    }
    out.push(format!(
        r#"<text x="10" y="{legend_y}" fill="{FG}">Square: male. Circle: female. Diamond: sex not known. Filled: sequenced. A line joins partners; children hang from it.</text>"#
    ));
    out.push("</svg>".into());
    Ok(out.join("\n") + "\n")
}
