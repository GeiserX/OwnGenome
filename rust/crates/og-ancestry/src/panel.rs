//! `panel-map` and `check-map`: the panel's sample table.
//!
//! `samples.tsv` has one row per panel person: sample, population,
//! continental and regional group (`-` when in none), whether the person
//! is a reference (unrelated, in a group, and the regional group has
//! [`MIN_REFERENCES`] unrelated people), and whether it is one of the
//! painting references (at most N per continental group, taken in turn from
//! each population).

use crate::{Args, Fail, Level, MIN_REFERENCES, Pending, Result, Spec, groups_of, label};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{self, Write};
use std::path::{Path, PathBuf};

pub(crate) const PANEL_MAP: Spec = Spec {
    verb: "panel-map",
    usage: "usage: og ancestry-fine panel-map --meta TSV --bcf-samples LIST [--exclude IDS] [--paint-per-group N] --out samples.tsv\n",
    options: &[
        "--meta",
        "--bcf-samples",
        "--exclude",
        "--paint-per-group",
        "--out",
    ],
    required: &["--meta", "--bcf-samples", "--out"],
};

pub(crate) const CHECK_MAP: Spec = Spec {
    verb: "check-map",
    usage: "usage: og ancestry-fine check-map --map samples.tsv\n",
    options: &["--map"],
    required: &["--map"],
};

/// The header of `samples.tsv`.
pub const HEADER: [&str; 6] = [
    "sample",
    "population",
    "continent",
    "region",
    "reference",
    "paint",
];

/// The meta TSV columns `panel-map` reads.
const META_SAMPLE: &str = "s";
const META_POP: &str = "population";
const META_RELATED: &str = "sample_filters.all_samples_related";

/// One row of `samples.tsv`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleRow {
    pub sample: String,
    pub population: String,
    pub continent: String,
    pub region: String,
    pub reference: bool,
    pub paint: bool,
}

impl SampleRow {
    fn group(&self, level: Level) -> &str {
        match level {
            Level::Continent => &self.continent,
            Level::Region => &self.region,
        }
    }
}

/// One metadata row: the three columns `panel-map` reads, `None` where the
/// row is too short to have one.
#[derive(Debug, Clone)]
pub struct Meta {
    pub sample: Option<String>,
    pub population: Option<String>,
    pub related: Option<String>,
}

/// The rows of `samples.tsv` for each panel sample, in the panel's order,
/// leaving out `exclude`. A person is a reference when unrelated, in a
/// group, and its regional group has `min_references` unrelated people;
/// a continental group is the union of its regional groups, so it keeps
/// that many too or leaves the model with them.
pub fn panel_map(
    meta: &[Meta],
    bcf_samples: &[String],
    exclude: &HashSet<String>,
    paint_per_group: i64,
    min_references: usize,
) -> std::result::Result<Vec<SampleRow>, String> {
    // The painting references of a group are its RFMix reference, so the
    // cap must not undo the minimum.
    let least = min_references.max(2);
    if paint_per_group < least as i64 {
        return Err(format!(
            "--paint-per-group is {paint_per_group}; each group needs at least {least} painting references"
        ));
    }
    let mut by_sample: HashMap<&str, &Meta> = HashMap::new();
    for m in meta {
        if let Some(s) = &m.sample {
            by_sample.insert(s, m);
        }
    }
    let missing: Vec<&String> = bcf_samples
        .iter()
        .filter(|s| !by_sample.contains_key(s.as_str()))
        .collect();
    if let Some(first) = missing.first() {
        return Err(format!(
            "{} panel samples have no metadata row, e.g. {first}",
            missing.len()
        ));
    }
    let mut rows: Vec<SampleRow> = bcf_samples
        .iter()
        .filter(|s| !exclude.contains(*s))
        .map(|s| {
            let m = by_sample[s.as_str()];
            // A short row has no population: Python wrote it as "None".
            let population = m.population.clone().unwrap_or_else(|| "None".to_owned());
            let groups = m.population.as_deref().and_then(groups_of);
            let (continent, region) = groups.unwrap_or(("-", "-"));
            SampleRow {
                sample: s.clone(),
                population,
                continent: continent.to_owned(),
                region: region.to_owned(),
                reference: groups.is_some() && m.related.as_deref() == Some("false"),
                paint: false,
            }
        })
        .collect();
    let mut per_region: HashMap<String, usize> = HashMap::new();
    for r in rows.iter().filter(|r| r.reference) {
        *per_region.entry(r.region.clone()).or_default() += 1;
    }
    for r in &mut rows {
        if r.reference && per_region[&r.region] < min_references {
            r.reference = false;
        }
    }
    // Painting references: in turn from each population of a continental
    // group (populations by name, people by sample id), up to the cap.
    let mut by_group: BTreeMap<&str, BTreeMap<&str, Vec<usize>>> = BTreeMap::new();
    for (i, r) in rows.iter().enumerate().filter(|(_, r)| r.reference) {
        by_group
            .entry(r.continent.as_str())
            .or_default()
            .entry(r.population.as_str())
            .or_default()
            .push(i);
    }
    let mut chosen = Vec::new();
    for pops in by_group.values() {
        let mut queues: Vec<std::collections::VecDeque<usize>> = pops
            .values()
            .map(|v| {
                let mut v = v.clone();
                v.sort_by(|&a, &b| rows[a].sample.cmp(&rows[b].sample));
                v.into()
            })
            .collect();
        let mut taken = 0i64;
        while taken < paint_per_group && queues.iter().any(|q| !q.is_empty()) {
            for q in &mut queues {
                if taken < paint_per_group
                    && let Some(i) = q.pop_front()
                {
                    chosen.push(i);
                    taken += 1;
                }
            }
        }
    }
    for i in chosen {
        rows[i].paint = true;
    }
    Ok(rows)
}

/// What makes a sample table unfit for the model: a regional group with
/// some but fewer than `min_references` reference people, a continental
/// group with reference people but fewer than 2 painting references, a
/// painting reference who is no reference, a reference person without both
/// groups.
pub fn map_problems(rows: &[SampleRow], min_references: usize) -> Vec<String> {
    let mut refs: BTreeMap<&str, usize> = BTreeMap::new();
    let mut conts: BTreeSet<&str> = BTreeSet::new();
    let mut paint: HashMap<&str, usize> = HashMap::new();
    let (mut stray, mut nogroup) = (Vec::new(), Vec::new());
    let blank = |g: &str| g == "-" || g.is_empty();
    for r in rows {
        if r.reference && (blank(&r.continent) || blank(&r.region)) {
            nogroup.push(r.sample.as_str());
        } else if r.reference {
            *refs.entry(&r.region).or_default() += 1;
            conts.insert(&r.continent);
            if r.paint {
                *paint.entry(&r.continent).or_default() += 1;
            }
        } else if r.paint {
            stray.push(r.sample.as_str());
        }
    }
    let mut out: Vec<String> = refs
        .iter()
        .filter(|&(_, &n)| n < min_references)
        .map(|(g, n)| {
            format!("regional group {g} has {n} of the {min_references} reference people it needs")
        })
        .collect();
    for g in &conts {
        let n = paint.get(g).copied().unwrap_or(0);
        if n < 2 {
            out.push(format!(
                "continental group {g} has {n} of the 2 painting references RFMix needs"
            ));
        }
    }
    if let Some(first) = stray.first() {
        out.push(format!(
            "painting references who are no reference people: {}, e.g. {first}",
            stray.len()
        ));
    }
    if let Some(first) = nogroup.first() {
        out.push(format!(
            "reference people without a continental and a regional group: {}, e.g. {first}",
            nogroup.len()
        ));
    }
    out
}

/// The groups `samples.tsv` rows have people of but the model leaves out,
/// as (level, code), for each level the model has (`codes`).
pub fn left_out_groups(
    rows: &[Vec<String>],
    codes: &BTreeMap<Level, HashSet<String>>,
) -> Vec<(Level, String)> {
    let mut out = Vec::new();
    for level in Level::ALL {
        let Some(have) = codes.get(&level) else {
            continue;
        };
        let col = level.column();
        let present: BTreeSet<&str> = rows
            .iter()
            .filter_map(|r| r.get(col).map(String::as_str))
            .filter(|c| *c != "-" && !c.is_empty())
            .collect();
        out.extend(
            present
                .into_iter()
                .filter(|c| !have.contains(*c))
                .map(|c| (level, c.to_owned())),
        );
    }
    out
}

/// Each left-out group as one name: "Oceania (continental and regional
/// group)" when both levels leave out a group of one name.
pub fn left_out_text(gone: &[(Level, String)]) -> String {
    let mut kinds: Vec<(&str, Vec<&str>)> = Vec::new();
    for (level, code) in gone {
        let name = label(*level, code);
        let kind = match level {
            Level::Continent => "continental",
            Level::Region => "regional",
        };
        match kinds.iter_mut().find(|(n, _)| *n == name) {
            Some((_, k)) => k.push(kind),
            None => kinds.push((name, vec![kind])),
        }
    }
    kinds
        .iter()
        .map(|(name, k)| format!("{name} ({} group)", k.join(" and ")))
        .collect::<Vec<_>>()
        .join(", ")
}

fn flag(b: bool) -> &'static str {
    if b { "1" } else { "0" }
}

pub(crate) fn panel_map_main(a: &Args, out: &mut dyn Write) -> Result<()> {
    let meta_path = a.req("--meta");
    let mut table = og_tsv::Table::open(meta_path).map_err(|e| Fail(e.to_string()))?;
    let mut meta = Vec::new();
    for row in table.by_ref() {
        let row = row.map_err(|e| Fail(e.to_string()))?;
        let get = |c: &str| row.get(c).map(str::to_owned);
        meta.push(Meta {
            sample: get(META_SAMPLE),
            population: get(META_POP),
            related: get(META_RELATED),
        });
    }
    if !meta.is_empty() {
        let names = table
            .header()
            .map_err(|e| Fail(e.to_string()))?
            .unwrap_or_default()
            .to_vec();
        for col in [META_SAMPLE, META_POP, META_RELATED] {
            if !names.iter().any(|n| n == col) {
                return Err(Fail(format!("{meta_path} has no column {col}")));
            }
        }
    }
    let samples: Vec<String> = crate::read_lines(Path::new(a.req("--bcf-samples")))?
        .iter()
        .map(|l| pycompat::num::py_strip(l))
        .filter(|l| !l.is_empty())
        .map(str::to_owned)
        .collect();
    let exclude_arg = a.get("--exclude").unwrap_or("").replace(',', " ");
    let exclude: HashSet<String> = crate::split_ws(&exclude_arg).map(str::to_owned).collect();
    let known: HashSet<&str> = samples.iter().map(String::as_str).collect();
    let unknown: BTreeSet<&str> = exclude
        .iter()
        .map(String::as_str)
        .filter(|x| !known.contains(x))
        .collect();
    if !unknown.is_empty() {
        return Err(Fail(format!(
            "--exclude names samples the panel does not have: {}",
            unknown.into_iter().collect::<Vec<_>>().join(", ")
        )));
    }
    let paint_per_group = a.int("--paint-per-group", 100).unwrap_or(100);
    let rows = panel_map(&meta, &samples, &exclude, paint_per_group, MIN_REFERENCES)?;
    let mut text = String::new();
    text.push_str(&HEADER.join("\t"));
    text.push('\n');
    for r in &rows {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            r.sample,
            r.population,
            r.continent,
            r.region,
            flag(r.reference),
            flag(r.paint)
        ));
    }
    let mut pending = Pending::new();
    pending.put(PathBuf::from(a.req("--out")), text.as_bytes())?;
    pending.commit_all()?;
    let n_ref = rows.iter().filter(|r| r.reference).count();
    let n_paint = rows.iter().filter(|r| r.paint).count();
    let say = |out: &mut dyn Write, s: String| -> Result<()> {
        writeln!(out, "{s}").map_err(|e: io::Error| Fail(format!("cannot write stdout: {e}")))
    };
    say(
        out,
        format!(
            "{} panel samples, {n_ref} reference individuals, {n_paint} painting references",
            rows.len()
        ),
    )?;
    let mut codes: BTreeMap<Level, HashSet<String>> = BTreeMap::new();
    for level in Level::ALL {
        codes.insert(
            level,
            rows.iter()
                .filter(|r| r.reference)
                .map(|r| r.group(level).to_owned())
                .collect(),
        );
    }
    let as_cells: Vec<Vec<String>> = rows
        .iter()
        .map(|r| {
            vec![
                r.sample.clone(),
                r.population.clone(),
                r.continent.clone(),
                r.region.clone(),
            ]
        })
        .collect();
    let gone = left_out_groups(&as_cells, &codes);
    if !gone.is_empty() {
        say(
            out,
            format!(
                "Left out of the model, fewer than {MIN_REFERENCES} unrelated reference people: {}",
                left_out_text(&gone)
            ),
        )?;
    }
    Ok(())
}

pub(crate) fn check_map_main(a: &Args, _out: &mut dyn Write) -> Result<()> {
    let path = a.req("--map");
    let lines = crate::read_lines(Path::new(path))?;
    let cells: Vec<Vec<&str>> = lines
        .iter()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    if let Some(bad) = cells.iter().find(|r| r.len() != 6) {
        return Err(Fail(format!(
            "{path}: a row without six columns, e.g. {}",
            crate::repr_list(bad)
        )));
    }
    let is_flag = |s: &str| s == "0" || s == "1";
    if let Some(bad) = cells.iter().find(|r| !is_flag(r[4]) || !is_flag(r[5])) {
        return Err(Fail(format!(
            "{path}: reference and paint must be 0 or 1, e.g. {}",
            crate::repr_list(bad)
        )));
    }
    let rows: Vec<SampleRow> = cells
        .iter()
        .map(|r| SampleRow {
            sample: r[0].to_owned(),
            population: r[1].to_owned(),
            continent: r[2].to_owned(),
            region: r[3].to_owned(),
            reference: r[4] == "1",
            paint: r[5] == "1",
        })
        .collect();
    if !rows.iter().any(|r| r.reference) {
        return Err(Fail(format!("{path}: no reference person")));
    }
    let problems = map_problems(&rows, MIN_REFERENCES);
    if !problems.is_empty() {
        return Err(Fail(format!("{path}: {}", problems.join("; "))));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{Dir, run};

    fn meta(s: &str, pop: &str, related: &str) -> Meta {
        Meta {
            sample: Some(s.to_owned()),
            population: Some(pop.to_owned()),
            related: Some(related.to_owned()),
        }
    }

    fn strings(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| (*s).to_owned()).collect()
    }

    fn small() -> Vec<Meta> {
        // 10 unrelated Yoruba, 3 unrelated San, one unrelated Papuan and two
        // related ones: the shape of gnomAD's Oceanians, 29 of 30 related.
        let mut m: Vec<Meta> = (0..10)
            .map(|i| meta(&format!("y{i}"), "YRI", "false"))
            .collect();
        m.extend((0..3).map(|i| meta(&format!("s{i}"), "San", "false")));
        m.push(meta("p0", "PapuanHighlands", "false"));
        m.push(meta("p1", "Bougainville", "true"));
        m.push(meta("p2", "PapuanSepik", "true"));
        m
    }

    fn ids(m: &[Meta]) -> Vec<String> {
        m.iter().map(|x| x.sample.clone().unwrap()).collect()
    }

    #[test]
    fn groups_references_and_painting_in_turn() {
        let m = vec![
            meta("y1", "YRI", "false"),
            meta("y2", "YRI", "false"),
            meta("y3", "YRI", "true"),
            meta("m1", "Mandenka", "false"),
            meta("i1", "IBS", "false"),
            meta("b1", "Basque", "false"),
            meta("t1", "TSI", "false"),
            meta("a1", "ASW", "false"),
            meta("x1", "CHB", "false"),
        ];
        let ex: HashSet<String> = ["x1".to_owned()].into();
        let rows = panel_map(&m, &ids(&m), &ex, 2, 1).unwrap();
        let by = |s: &str| rows.iter().find(|r| r.sample == s).cloned();
        assert!(by("x1").is_none(), "an excluded person is not in the table");
        let i1 = by("i1").unwrap();
        assert_eq!((i1.continent.as_str(), i1.region.as_str()), ("WEA", "IBE"));
        assert_eq!(by("b1").unwrap().region, "IBE");
        assert!(
            !by("y3").unwrap().reference,
            "a related person is no reference"
        );
        let a1 = by("a1").unwrap();
        assert_eq!((a1.continent.as_str(), a1.reference), ("-", false));
        let mut afr: Vec<&str> = rows
            .iter()
            .filter(|r| r.continent == "AFR" && r.paint)
            .map(|r| r.sample.as_str())
            .collect();
        afr.sort_unstable();
        assert_eq!(afr, ["m1", "y1"], "in turn: Mandenka before a second YRI");
    }

    #[test]
    fn refusals() {
        let m = vec![meta("y1", "YRI", "false")];
        let e = panel_map(&m, &strings(&["y1", "zz"]), &HashSet::new(), 100, 10).unwrap_err();
        assert_eq!(e, "1 panel samples have no metadata row, e.g. zz");
        let e = panel_map(&m, &strings(&["y1"]), &HashSet::new(), 9, 10).unwrap_err();
        assert_eq!(
            e,
            "--paint-per-group is 9; each group needs at least 10 painting references"
        );
        assert!(panel_map(&m, &strings(&["y1"]), &HashSet::new(), 1, 1).is_err());
    }

    #[test]
    fn a_small_group_is_left_out() {
        let m = small();
        let rows = panel_map(&m, &ids(&m), &HashSet::new(), 50, MIN_REFERENCES).unwrap();
        let by = |s: &str| rows.iter().find(|r| r.sample == s).unwrap();
        assert_eq!((by("p0").reference, by("p0").paint), (false, false));
        assert_eq!(by("p0").continent, "OCE");
        assert!((0..3).all(|i| !by(&format!("s{i}")).reference));
        assert!((0..10).all(|i| by(&format!("y{i}")).reference && by(&format!("y{i}")).paint));
        let codes: BTreeMap<Level, HashSet<String>> = [
            (Level::Continent, ["AFR".to_owned()].into()),
            (Level::Region, ["WAF".to_owned()].into()),
        ]
        .into();
        let cells: Vec<Vec<String>> = rows
            .iter()
            .map(|r| {
                vec![
                    r.sample.clone(),
                    r.population.clone(),
                    r.continent.clone(),
                    r.region.clone(),
                ]
            })
            .collect();
        let gone = left_out_groups(&cells, &codes);
        assert_eq!(
            gone,
            [
                (Level::Continent, "OCE".to_owned()),
                (Level::Region, "OCE".to_owned()),
                (Level::Region, "SAF".to_owned())
            ]
        );
        assert_eq!(
            left_out_text(&gone),
            "Oceania (continental and regional group), Southern Africa (regional group)"
        );
        assert!(map_problems(&rows, MIN_REFERENCES).is_empty());
        let mut one = rows.clone();
        for r in &mut one {
            if r.sample == "p0" {
                r.reference = true;
                r.paint = true;
            }
        }
        assert_eq!(
            map_problems(&one, MIN_REFERENCES),
            [
                "regional group OCE has 1 of the 10 reference people it needs",
                "continental group OCE has 1 of the 2 painting references RFMix needs"
            ]
        );
        let nopaint: Vec<SampleRow> = rows
            .iter()
            .cloned()
            .map(|r| SampleRow { paint: false, ..r })
            .collect();
        assert_eq!(
            map_problems(&nopaint, MIN_REFERENCES),
            ["continental group AFR has 0 of the 2 painting references RFMix needs"]
        );
        let stray: Vec<SampleRow> = rows
            .iter()
            .cloned()
            .map(|r| SampleRow {
                paint: r.paint || r.sample == "p1",
                ..r
            })
            .collect();
        assert_eq!(
            map_problems(&stray, MIN_REFERENCES),
            ["painting references who are no reference people: 1, e.g. p1"]
        );
        let nogrp: Vec<SampleRow> = rows
            .iter()
            .cloned()
            .map(|r| {
                if r.sample == "y0" {
                    SampleRow {
                        region: "-".into(),
                        ..r
                    }
                } else {
                    r
                }
            })
            .collect();
        assert_eq!(
            map_problems(&nogrp, MIN_REFERENCES),
            [
                "regional group WAF has 9 of the 10 reference people it needs",
                "reference people without a continental and a regional group: 1, e.g. y0"
            ]
        );
    }

    #[test]
    fn the_commands_write_and_refuse() {
        let d = Dir::new("panel");
        let m = small();
        let mut text = format!("{META_SAMPLE}\t{META_POP}\t{META_RELATED}\n");
        for x in &m {
            text.push_str(&format!(
                "{}\t{}\t{}\n",
                x.sample.as_deref().unwrap(),
                x.population.as_deref().unwrap(),
                x.related.as_deref().unwrap()
            ));
        }
        let meta_p = d.put("meta.tsv", &text);
        let bcf = d.put(
            "bcf.txt",
            &ids(&m).iter().map(|s| format!("{s}\n")).collect::<String>(),
        );
        let out = d.path("samples.tsv");
        let (rc, stdout, err) = run(&[
            "panel-map",
            "--meta",
            &meta_p,
            "--bcf-samples",
            &bcf,
            "--out",
            &out,
        ]);
        assert_eq!(rc, 0, "{err}");
        assert!(stdout.contains(
            "Left out of the model, fewer than 10 unrelated reference people: Oceania (continental and regional group), Southern Africa (regional group)"
        ));
        assert!(
            d.read("samples.tsv")
                .contains("p0\tPapuanHighlands\tOCE\tOCE\t0\t0\n")
        );
        assert_eq!(run(&["check-map", "--map", &out]).0, 0);
        let (rc, _, err) = run(&[
            "panel-map",
            "--meta",
            &meta_p,
            "--bcf-samples",
            &bcf,
            "--exclude",
            "y1,zz",
            "--out",
            &out,
        ]);
        assert_eq!(
            (rc, err.as_str()),
            (
                1,
                "ERROR: --exclude names samples the panel does not have: zz\n"
            )
        );
        let short = d.put(
            "short.tsv",
            "sample\tpopulation\tcontinent\tregion\treference\tpaint\ny0\tYRI\tAFR\n",
        );
        let (rc, _, err) = run(&["check-map", "--map", &short]);
        assert_eq!(rc, 1);
        assert!(
            err.ends_with("a row without six columns, e.g. ['y0', 'YRI', 'AFR']\n"),
            "{err}"
        );
        let flag = d.put("flag.tsv", "h\ny10\tYRI\tAFR\tWAF\ttrue\t0\n");
        assert!(
            run(&["check-map", "--map", &flag])
                .2
                .contains("reference and paint must be 0 or 1")
        );
        let empty = d.put(
            "empty.tsv",
            "sample\tpopulation\tcontinent\tregion\treference\tpaint\n",
        );
        assert!(
            run(&["check-map", "--map", &empty])
                .2
                .ends_with(": no reference person\n")
        );
        let (rc, _, err) = run(&["check-map", "--map", &d.path("nope.tsv")]);
        assert_eq!(rc, 1);
        assert!(err.starts_with("ERROR: cannot read "), "{err}");
    }
}
