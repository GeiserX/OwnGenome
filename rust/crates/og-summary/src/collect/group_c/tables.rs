//! The sections that read one step table each: the couple screen, the
//! polygenic scores, fine ancestry, telomeres, the haplogroups, the Y-STR
//! panel, the raw-data export, blood groups, archaic and ancient ancestry,
//! Lp(a), VNtyper, Locityper, viral DNA, CHIP, Exomiser, mobile elements
//! and MoChA.

use super::super::model::State;
use super::util::{
    Obj, Table, at, base, count, dots, dots_strict, get_default, isdigit, kv_dot, kv_get,
    metric_map, opt_s, or, or_dot, ragged, round, s, strs, table,
};
use super::{Ctx, Found, Others, Out};
use crate::collect::sections::{num, py_error};
use indexmap::IndexMap;
use pycompat::json::Value;

fn unreadable<T>(msg: String) -> Result<T, String> {
    Err(msg)
}

/// Python's `sorted(need)` of a set, as its list repr.
fn sorted_repr(need: &[&str]) -> String {
    let mut v: Vec<&str> = need.to_vec();
    v.sort_unstable();
    format!(
        "[{}]",
        v.iter()
            .map(|k| crate::repr(k))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn float_err(x: &str) -> String {
    py_error("ValueError", &format!("could not convert string to float: {}", crate::repr(x)))
}

fn int_of(x: &str) -> Result<i128, String> {
    crate::int(x).map_err(|e| py_error("ValueError", &e))
}

fn opt_num(x: Option<&str>) -> Value {
    num(x).map_or(Value::Null, Value::Float)
}

/// `{k: r[k] or "."}` over the header's own first line, read as text.
fn first_line(path: &str) -> Result<Vec<String>, String> {
    let l = super::util::lines(path)?;
    let h = l.first().map_or("", String::as_str);
    Ok(h.trim_end_matches('\r').split('\t').map(str::to_owned).collect())
}

/// Base64 of a file's bytes (`base64.b64encode`).
fn b64(bytes: &[u8]) -> String {
    const T: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for ch in bytes.chunks(3) {
        let b = [ch[0], *ch.get(1).unwrap_or(&0), *ch.get(2).unwrap_or(&0)];
        let n = (u32::from(b[0]) << 16) | (u32::from(b[1]) << 8) | u32::from(b[2]);
        for i in 0..4 {
            if i <= ch.len() {
                out.push(T[((n >> (18 - 6 * i)) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

/// The base64 of a picture under `limit` bytes, else `null`.
fn picture(f: Option<Found>, limit: u64) -> Result<Value, String> {
    let Some(f) = f else { return Ok(Value::Null) };
    let size = std::fs::metadata(&f.path)
        .map_err(|e| super::util::os_note(&crate::os_error(&e, &f.path)))?
        .len();
    if size >= limit {
        return Ok(Value::Null);
    }
    let bytes =
        std::fs::read(&f.path).map_err(|e| super::util::os_note(&crate::os_error(&e, &f.path)))?;
    Ok(s(&b64(&bytes)))
}

// --- step 39 --------------------------------------------------------------------

pub fn couple(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("couple/{}_couple_carrier.tsv", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, true)?;
    let need = [
        "gene",
        "condition",
        "result",
        "offspring_risk",
        "partner_a",
        "findings_a",
        "partner_b",
        "findings_b",
        "explanation",
    ];
    if need.iter().any(|k| !t.has(k)) {
        return unreadable(format!(
            "{} lacks a column of {}",
            base(&p.path),
            need.join(", ")
        ));
    }
    let sample = c.sample;
    let cell = |r: &pycompat::csv::DictRow, k: &str| r.get(k).flatten().unwrap_or("").to_owned();
    for (i, r) in t.rows.iter().enumerate() {
        let bad = ragged(r, &need)
            || (cell(r, "partner_a") != sample && cell(r, "partner_b") != sample)
            || !matches!(
                cell(r, "result").as_str(),
                "shared" | "x_linked" | "one_partner"
            );
        if bad {
            return unreadable(format!(
                "{} line {} is not a row of {sample}'s couple table",
                base(&p.path),
                i + 2
            ));
        }
    }
    let side = |r: &pycompat::csv::DictRow| {
        let (own, other) = if cell(r, "partner_a") == sample {
            ("a", "b")
        } else {
            ("b", "a")
        };
        Obj::new()
            .set("gene", s(&cell(r, "gene")))
            .set("condition", s(&cell(r, "condition")))
            .set("risk", s(&cell(r, "offspring_risk")))
            .set("own", s(&cell(r, &format!("findings_{own}"))))
            .set("partner_findings", s(&cell(r, &format!("findings_{other}"))))
            .set("explanation", s(&cell(r, "explanation")))
            .val()
    };
    let partner = t.rows.first().map(|r| {
        if cell(r, "partner_a") == sample {
            cell(r, "partner_b")
        } else {
            cell(r, "partner_a")
        }
    });
    let of = |kind: &str| {
        Value::Array(
            t.rows
                .iter()
                .filter(|r| cell(r, "result") == kind)
                .map(side)
                .collect(),
        )
    };
    let v = Obj::new()
        .set("partner", opt_s(partner.as_deref()))
        .set("shared", of("shared"))
        .set("x_linked", of("x_linked"))
        .set(
            "one_partner",
            count(t.rows.iter().filter(|r| cell(r, "result") == "one_partner").count()),
        )
        .val();
    Ok(Some((p, v)))
}

// --- step 25 --------------------------------------------------------------------

pub fn prs(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find1(format!("prs/{s_}_prs_summary.tsv")) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    let mut rows = Vec::new();
    let mut adjusted = false;
    for r in &t.rows {
        let pct = or(r, "Percentile", "");
        let pct = if pct == "NA" { String::new() } else { pct };
        adjusted |= !pct.is_empty();
        let group = match r.get("Ancestry_Group") {
            None => s(""),
            Some(Some("NA")) => s(""),
            Some(v) => opt_s(v),
        };
        let evaluated = match r.get("Evaluated_In") {
            Some(Some(x)) => strs(x.split(';').filter(|x| !x.is_empty())),
            _ => Value::Null,
        };
        rows.push(
            Obj::new()
                .set("condition", get_default(r, "Condition", ""))
                .set("pgs_id", get_default(r, "PGS_ID", ""))
                .set("score", get_default(r, "Score_SUM", ""))
                .set("matched", get_default(r, "Variants_Matched", ""))
                .set("total", get_default(r, "Variants_Total", ""))
                .set("chrx", get_default(r, "ChrX", ""))
                .set("ambiguous", get_default(r, "Ambiguous", ""))
                .set("applies_to", get_default(r, "Applies_To", ""))
                .set("sex_check", get_default(r, "Sex_Check", ""))
                .set("percentile", s(&pct))
                .set("group", group)
                .set("evaluated_in", evaluated)
                .set("effect", s(&or(r, "Effect", "")))
                .val(),
        );
    }
    let mut out = Obj::new()
        .set("scores", Value::Array(rows))
        .set("adjusted", Value::Bool(adjusted));
    if let Some(a) = c.find(&[
        format!("ancestry/{s_}_ancestry.tsv"),
        format!("prs/{s_}_ancestry.tsv"),
    ]) {
        let mut kv: IndexMap<String, Value> = IndexMap::new();
        for r in table(&a.path, false)?.rows {
            if let Some(Some(k)) = r.get("key")
                && !k.is_empty()
            {
                kv.insert(k.to_owned(), get_default(&r, "value", ""));
            }
        }
        let g = |k: &str| kv.get(k).cloned().unwrap_or_else(|| s(""));
        out.put(
            "ancestry",
            Obj::new()
                .set("population", g("population"))
                .set("panel", g("reference_panel"))
                .set("low_confidence", g("population_low_confidence"))
                .val(),
        );
    }
    Ok(Some((p, out.val())))
}

// --- step 26b -------------------------------------------------------------------

pub fn ancestry_fine(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find1(format!("ancestry_fine/{s_}_ancestry_fine.tsv")) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    let need = ["level", "code", "group", "fraction", "se", "low", "high"];
    if !t.first_has(&need) {
        return unreadable(format!(
            "{} has no rows or lacks a column of {}",
            base(&p.path),
            sorted_repr(&need)
        ));
    }
    let levels = ["continent", "region", "painting", "left_out", "not_painted"];
    let mut by: IndexMap<&str, Vec<Value>> = levels.iter().map(|l| (*l, Vec::new())).collect();
    for r in &t.rows {
        let Some(Some(level)) = r.get("level") else {
            continue;
        };
        if let Some(list) = by.get_mut(level) {
            list.push(
                Obj::new()
                    .set("code", opt_s(at(r, "code")?))
                    .set("group", opt_s(at(r, "group")?))
                    .set("fraction", opt_num(at(r, "fraction")?))
                    .set("low", opt_num(at(r, "low")?))
                    .set("high", opt_num(at(r, "high")?))
                    .val(),
            );
        }
    }
    let mut out = Obj::new();
    for (k, v) in by {
        out.put(k, Value::Array(v));
    }
    out.put("min_references", count(og_ancestry::MIN_REFERENCES));
    let png = c.find1(format!("ancestry_fine/{s_}_painting.png"));
    out.put("painting_png", Value::Null);
    out.put(
        "legend",
        Value::Array(
            og_ancestry::COLOURS
                .iter()
                .map(|(code, rgb)| {
                    Obj::new()
                        .set("code", s(code))
                        .set(
                            "group",
                            s(og_ancestry::label(og_ancestry::Level::Continent, code)),
                        )
                        .set(
                            "colour",
                            s(&format!("#{:02x}{:02x}{:02x}", rgb[0], rgb[1], rgb[2])),
                        )
                        .val()
                })
                .collect(),
        ),
    );
    out.put("painting_png", picture(png, 2_000_000)?);
    Ok(Some((p, out.val())))
}

// --- step 44 --------------------------------------------------------------------

pub fn ancient(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find1(format!("ancient/{s_}_ancient.tsv")) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    let need = ["rank", "group", "n", "distance", "dates", "country", "locality"];
    if !t.first_has(&need) {
        return unreadable(format!(
            "{} has no rows or lacks a column of {}",
            base(&p.path),
            sorted_repr(&need)
        ));
    }
    let info = match c.find1(format!("ancient/{s_}_ancient_info.tsv")) {
        Some(ip) => metric_map(&table(&ip.path, false)?, "key", "value"),
        None => IndexMap::new(),
    };
    let nearest = t
        .rows
        .iter()
        .map(|r| {
            need.iter()
                .fold(Obj::new(), |o, k| o.set(k, opt_s(r.get(k).flatten())))
                .val()
        })
        .collect();
    let mut out = Obj::new().set("nearest", Value::Array(nearest));
    for k in [
        "release",
        "sites_used",
        "panel_sites",
        "ancient_people",
        "ancient_groups",
        "present_day_people",
    ] {
        out.put(k, s(&kv_dot(&info, k)));
    }
    let svg = c.find1(format!("ancient/{s_}_ancient_pca.svg"));
    out.put("plot_svg", picture(svg, 3_000_000)?);
    Ok(Some((p, out.val())))
}

// --- step 10 --------------------------------------------------------------------

const GENOME_BASES: f64 = 3.1e9;
const TELOMERE_FLOOR_PER_30X: f64 = 2000.0;

pub fn telomere(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find(&[
        format!("telomere/{s_}/{s_}/{s_}_summary.tsv"),
        format!("telomere/{s_}/{s_}_summary.tsv"),
        "telomere/**/*_summary.tsv".to_owned(),
    ]) else {
        return Ok(None);
    };
    let lines: Vec<Vec<String>> = super::util::lines(&p.path)?
        .into_iter()
        .filter(|l| !pycompat::num::py_strip(l).is_empty())
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    if lines.len() < 2 {
        return Ok(Some((p, Obj::new().set("tel_content", s(".")).val())));
    }
    let (head, row) = (&lines[0], &lines[1]);
    let i = head.iter().position(|h| h == "tel_content").unwrap_or(10);
    let mut out = Obj::new().set("tel_content", s(row.get(i).map_or(".", String::as_str)));
    let col = |names: &[&str]| -> Option<&str> {
        names.iter().find_map(|n| {
            head.iter()
                .position(|h| h == n)
                .and_then(|i| row.get(i))
                .map(String::as_str)
        })
    };
    let intratel = num(col(&["intratel_reads", "intratelomeric_reads"]));
    let total = num(col(&["total_reads"]));
    let lengths: Vec<f64> = col(&["read_lengths", "read_length"])
        .unwrap_or("")
        .split([',', ';', '|', ' '])
        .filter_map(|x| num(Some(x)))
        .filter(|x| *x != 0.0)
        .collect();
    let Some(intratel) = intratel else {
        out.put(
            "read_check",
            s("no intratelomeric read count in the summary"),
        );
        return Ok(Some((p, out.val())));
    };
    out.put(
        "intratelomeric_reads",
        Value::Int(crate::trunc(intratel).map_err(|e| py_error("OverflowError", &e))?),
    );
    let total = total.filter(|t| *t != 0.0);
    let (Some(total), false) = (total, lengths.is_empty()) else {
        out.put(
            "read_check",
            s("depth unknown (no total_reads or read_lengths in the summary): not checked"),
        );
        return Ok(Some((p, out.val())));
    };
    let longest = lengths.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let depth = total * longest / GENOME_BASES;
    let floor = TELOMERE_FLOOR_PER_30X * depth / 30.0;
    out.put("depth", round(depth, 1));
    let fl = pycompat::fmt::round_int(floor)
        .map_err(|e| py_error(e.exception(), "cannot round"))
        .and_then(|x| crate::trunc(x).map_err(|e| py_error("OverflowError", &e)))?;
    out.put("read_floor", Value::Int(fl));
    out.put("low_reads", Value::Bool(intratel < floor));
    Ok(Some((p, out.val())))
}

// --- steps 12 and 37 ------------------------------------------------------------

const EXPLORE_DOCS: &str = "https://github.com/GeiserX/OwnGenome/blob/main/docs/";

/// `urllib.parse.quote(x, safe="")`.
fn quote(x: &str) -> String {
    let mut out = String::new();
    for b in x.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-' | b'~') {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

fn explore(kind: &str, hg: &str) -> Value {
    let mut hg = pycompat::num::py_strip(hg);
    if kind == "y" {
        hg = hg.split('*').next().unwrap_or("");
    }
    if matches!(hg, "" | "." | "NA" | "insufficient markers") {
        return Value::Array(Vec::new());
    }
    let q = quote(hg);
    let link = |name: &str, url: String| Obj::new().set("name", s(name)).set("url", s(&url)).val();
    let (ftdna, yfull, doc) = if kind == "y" {
        (
            format!("https://discover.familytreedna.com/y-dna/{q}"),
            ("YFull YTree", format!("https://www.yfull.com/tree/{q}/")),
            "37-y-haplogroup.md#explore-your-branch-further",
        )
    } else {
        (
            format!("https://discover.familytreedna.com/mtdna/{q}"),
            ("YFull MTree", format!("https://www.yfull.com/mtree/{q}/")),
            "12-mito-haplogroup.md#explore-your-branch-further",
        )
    };
    Value::Array(vec![
        link("FamilyTreeDNA Discover (free, no upload)", ftdna),
        link(yfull.0, yfull.1),
        link(
            "More routes, and the files to share",
            format!("{EXPLORE_DOCS}{doc}"),
        ),
    ])
}

pub fn haplogroup(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find1(format!("mito/{s_}_haplogroup.txt")) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    let hg = match t.rows.first() {
        None => String::new(),
        Some(r) => {
            let named = r.get("Haplogroup").flatten().filter(|x| !x.is_empty());
            let second = || {
                r.iter()
                    .nth(1)
                    .and_then(|(_, v)| v)
                    .filter(|x| !x.is_empty())
            };
            named.or_else(second).unwrap_or("").trim_matches('"').to_owned()
        }
    };
    let mut out = Obj::new()
        .set("haplogroup", s(if hg.is_empty() { "." } else { &hg }))
        .set("explore", explore("mt", &hg));
    if let Some(hc) = c.find1(format!("mito/{s_}_haplocheck.txt")) {
        let t = table(&hc.path, false)?;
        let first: Option<IndexMap<String, String>> = t.rows.first().map(|r| {
            r.iter()
                .filter(|(k, _)| !k.is_empty())
                .map(|(k, v)| {
                    (
                        k.trim_matches('"').to_owned(),
                        v.unwrap_or("").trim_matches('"').to_owned(),
                    )
                })
                .collect()
        });
        let Some(first) = first.filter(|m| m.contains_key("Contamination Status")) else {
            return unreadable(format!(
                "{} has no Contamination Status column",
                base(&hc.path)
            ));
        };
        let dot = |k: &str| {
            first
                .get(k)
                .filter(|v| !v.is_empty())
                .map_or(".", String::as_str)
                .to_owned()
        };
        out.put("contamination_status", s(&dot("Contamination Status")));
        out.put("contamination_level", s(&dot("Contamination Level")));
    }
    Ok(Some((p, out.val())))
}

pub fn y_haplogroup(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("y_haplogroup/{}_y_haplogroup.txt", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    let Some(r) = t.rows.first().filter(|_| t.has("Hg")) else {
        return unreadable(format!(
            "{} has no Hg column or no sample row",
            base(&p.path)
        ));
    };
    let hg = pycompat::num::py_strip(r.get("Hg").flatten().unwrap_or("")).to_owned();
    let v = Obj::new()
        .set(
            "haplogroup",
            s(if hg.is_empty() || hg == "NA" {
                "insufficient markers"
            } else {
                &hg
            }),
        )
        .set("valid_markers", s(&or_dot(r, "Valid_markers")))
        .set("qc_score", s(&or_dot(r, "QC-score")))
        .set("explore", explore("y", &hg))
        .val();
    Ok(Some((p, v)))
}

pub fn y_tree(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find1(format!("y_haplogroup/{s_}_y_tree.txt")) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    let need = ["Hg", "Level", "Tree"];
    let Some(r) = t.rows.first().filter(|_| need.iter().all(|k| t.has(k))) else {
        return unreadable(format!(
            "{} lacks a column of {}",
            base(&p.path),
            need.join(", ")
        ));
    };
    let hg = pycompat::num::py_strip(r.get("Hg").flatten().unwrap_or("")).to_owned();
    let mut out = Obj::new().set(
        "haplogroup",
        s(if hg.is_empty() || hg == "NA" {
            "no branch has a derived call"
        } else {
            &hg
        }),
    );
    for (k, col) in [
        ("level", "Level"),
        ("tree", "Tree"),
        ("path_derived", "Path_derived"),
        ("path_ancestral", "Path_ancestral"),
        ("path_no_call", "Path_no_call"),
        ("ancestral_in_terminal", "Ancestral_in_terminal"),
        ("outvoted_ancestral", "Outvoted_ancestral"),
        ("snps_called", "SNPs_called"),
        ("min_depth", "Min_depth"),
    ] {
        out.put(k, s(&or_dot(r, col)));
    }
    out.put("explore", explore("y", &hg));
    if let Some(lvl) = c.find1(format!("y_haplogroup/{s_}_y_tree_path.tsv")) {
        let pt = table(&lvl.path, false)?;
        if !pt.has("Branch") {
            return unreadable(format!("{} has no Branch column", base(&lvl.path)));
        }
        let keys = [
            "Level",
            "Branch",
            "SNPs",
            "On_tables",
            "Derived",
            "Ancestral",
            "No_call",
            "Ancestral_SNPs",
        ];
        out.put(
            "path",
            Value::Array(
                pt.rows
                    .iter()
                    .map(|row| {
                        keys.iter()
                            .fold(Obj::new(), |o, k| o.set(k, s(&or(row, k, ""))))
                            .val()
                    })
                    .collect(),
            ),
        );
    }
    Ok(Some((p, out.val())))
}

// --- step 45 --------------------------------------------------------------------

pub fn y_str(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("y_str/{}_y_str.tsv", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    let need = ["order", "panel", "marker", "count", "confidence"];
    if !t.first_has(&need) {
        return unreadable(format!(
            "{} has no rows or lacks a column of {}",
            base(&p.path),
            need.join(", ")
        ));
    }
    let keep: Vec<Vec<String>> = t
        .rows
        .iter()
        .map(|r| need.iter().map(|k| or(r, k, "")).collect())
        .collect();
    let kinds = ["high", "low", "no call", "not callable"];
    let counts: Vec<usize> = kinds
        .iter()
        .map(|c| keep.iter().filter(|r| r[4] == *c).count())
        .collect();
    if counts.iter().sum::<usize>() != keep.len() {
        return unreadable(format!(
            "{} has a confidence other than high, low, no call or not callable",
            base(&p.path)
        ));
    }
    let mut out = Obj::new().set(
        "markers",
        Value::Array(
            keep.iter()
                .map(|r| {
                    need.iter()
                        .zip(r)
                        .fold(Obj::new(), |o, (k, v)| o.set(k, s(v)))
                        .val()
                })
                .collect(),
        ),
    );
    for (k, n) in kinds.iter().zip(counts) {
        out.put(&k.replace(' ', "_"), count(n));
    }
    Ok(Some((p, out.val())))
}

// --- step 40 --------------------------------------------------------------------

pub fn consumer_export(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("consumer_export/{}_consumer_export.tsv", c.s())) else {
        return Ok(None);
    };
    let mut data: IndexMap<String, Option<String>> = IndexMap::new();
    for r in table(&p.path, false)?.rows {
        if let Some(Some(k)) = r.get("key")
            && !k.is_empty()
        {
            let v = match r.get("value") {
                None => Some(String::new()),
                Some(v) => v.map(str::to_owned),
            };
            data.insert(k.to_owned(), v);
        }
    }
    let text = |k: &str| data.get(k).cloned().flatten().unwrap_or_default();
    for k in ["sites", "called", "no_call"] {
        if !isdigit(&text(k)) {
            return unreadable(format!(
                "{} has no whole-number '{k}' row",
                base(&p.path)
            ));
        }
    }
    let mut reasons = Obj::new();
    for (k, v) in &data {
        if let Some(why) = k.strip_prefix("no_call: ") {
            let v = v.clone().unwrap_or_default();
            if !isdigit(&v) {
                return unreadable(format!(
                    "{} has a '{k}' row without a whole-number count",
                    base(&p.path)
                ));
            }
            reasons.put(why, Value::Int(int_of(&v)?));
        }
    }
    let blocks = text("from_reference_blocks");
    let blocks = if blocks.is_empty() { 0 } else { int_of(&blocks)? };
    let dot = |k: &str| {
        let v = text(k);
        s(if v.is_empty() { "." } else { &v })
    };
    let v = Obj::new()
        .set("sites", Value::Int(int_of(&text("sites"))?))
        .set("called", Value::Int(int_of(&text("called"))?))
        .set("no_call", Value::Int(int_of(&text("no_call"))?))
        .set("from_reference_blocks", Value::Int(blocks))
        .set("template_sites", dot("template_sites"))
        .set("file_23andme", dot("file_23andme"))
        .set("file_ancestrydna", dot("file_ancestrydna"))
        .set("no_call_reasons", reasons.val())
        .val();
    Ok(Some((p, v)))
}

// --- step 41 --------------------------------------------------------------------

pub fn bloodtype(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("bloodtype/{}_bloodtype.tsv", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    if !t.first_has(&["system", "phenotype", "confidence"]) {
        return unreadable(format!(
            "{} has no system, phenotype and confidence columns, or no row",
            base(&p.path)
        ));
    }
    let mut by: IndexMap<Option<String>, String> = IndexMap::new();
    let mut systems = Vec::new();
    for r in &t.rows {
        let system = r.get("system").flatten();
        let pheno = or_dot(r, "phenotype");
        by.insert(system.map(str::to_owned), pheno.clone());
        systems.push(
            Obj::new()
                .set("system", opt_s(system))
                .set("phenotype", s(&pheno))
                .set("genotype", s(&or_dot(r, "genotype")))
                .set("confidence", s(&or_dot(r, "confidence")))
                .val(),
        );
    }
    let get = |k: &str| s(by.get(&Some(k.to_owned())).map_or(".", String::as_str));
    let v = Obj::new()
        .set("systems", Value::Array(systems))
        .set("abo", get("ABO"))
        .set("rhd", get("RhD"))
        .val();
    Ok(Some((p, v)))
}

// --- step 43 --------------------------------------------------------------------

fn all_digits(x: &str) -> bool {
    !x.is_empty() && x.bytes().all(|b| b.is_ascii_digit())
}

pub fn archaic(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("archaic/{}_archaic.tsv", c.s())) else {
        return Ok(None);
    };
    let mut data: IndexMap<String, String> = IndexMap::new();
    for r in table(&p.path, false)?.rows {
        if let Some(Some(k)) = r.get("key")
            && !k.is_empty()
        {
            data.insert(k.to_owned(), or(&r, "value", ""));
        }
    }
    let mut longest: Vec<(i128, &String)> = data
        .iter()
        .filter_map(|(k, v)| {
            let n = k.strip_prefix("longest_").filter(|n| all_digits(n))?;
            Some((crate::int(n).ok()?, v))
        })
        .collect();
    longest.sort_by_key(|x| x.0);
    let mut out = Obj::new()
        .set(
            "population",
            s(data.get("population").map_or(".", String::as_str)),
        )
        .set("longest", strs(longest.iter().map(|x| x.1.as_str())));
    for k in [
        "neanderthal_segments",
        "neanderthal_bp",
        "denisovan_segments",
        "denisovan_bp",
        "analysed_bp",
    ] {
        let v = data.get(k).map_or("", String::as_str);
        if !all_digits(v) {
            return unreadable(format!(
                "{} has no whole number for {k}",
                base(&p.path)
            ));
        }
        out.put(k, Value::Int(int_of(v)?));
    }
    for k in ["neanderthal_percent", "denisovan_percent"] {
        let v = data.get(k).map_or("", String::as_str);
        let ok = v
            .split_once('.')
            .is_some_and(|(a, b)| all_digits(a) && all_digits(b));
        if !ok {
            return unreadable(format!(
                "{} has no percentage for {k}",
                base(&p.path)
            ));
        }
        out.put(k, s(v));
    }
    Ok(Some((p, out.val())))
}

// --- step 46 --------------------------------------------------------------------

const LPA_RISK_ALLELES: [(&str, &str); 2] = [("rs10455872", "G"), ("rs3798220", "C")];

pub fn lpa(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find1(format!("lpa/{s_}_lpa_kiv2.tsv")) else {
        return Ok(None);
    };
    let kv = metric_map(&table(&p.path, false)?, "metric", "value");
    if !matches!(
        kv_get(&kv, "status"),
        Some("estimated" | "too little depth" | "not assessed")
    ) {
        return unreadable(format!(
            "{} has no status row (estimated, too little depth or not assessed)",
            base(&p.path)
        ));
    }
    let mut data = Obj::new();
    for k in [
        "status",
        "kiv2_copies",
        "ci95_low",
        "ci95_high",
        "class",
        "kiv2_depth",
        "flank_depth",
        "mapq0_share",
        "message",
        "advice",
    ] {
        data.put(k, s(&kv_dot(&kv, k)));
    }
    let dir = p.path.rsplit_once('/').map_or("", |x| x.0);
    data.put(
        "null_allele",
        super::cn::read_marker_table(
            &crate::join(dir, &format!("{}_lpa_kiv2_4925.tsv", c.sample)),
            "46",
            "LPA_KIV2_4925",
        ),
    );
    let mut variants = Vec::new();
    let mut risk = 0i128;
    let mut any_read = false;
    if let Some(t) = c.find1(format!("traits/{s_}_traits.tsv")) {
        for r in table(&t.path, false)?.rows {
            let rsid = r.get("rsid").flatten();
            let Some(alt) = LPA_RISK_ALLELES
                .iter()
                .find(|(id, _)| Some(*id) == rsid)
                .map(|x| x.1)
            else {
                continue;
            };
            let called = r.get("status").flatten() == Some("called");
            let gt = r.get("genotype").flatten().unwrap_or("");
            let n = called.then(|| gt.split(['/', '|']).filter(|a| *a == alt).count() as i128);
            if let Some(n) = n {
                risk += n;
                any_read = true;
            }
            variants.push(
                Obj::new()
                    .set("rsid", opt_s(rsid))
                    .set("genotype", s(&or_dot(&r, "genotype")))
                    .set("status", s(&or_dot(&r, "status")))
                    .set("risk_alleles", n.map_or(Value::Null, Value::Int))
                    .val(),
            );
        }
    }
    let class = kv_dot(&kv, "class");
    let mut reasons = Vec::new();
    if risk != 0 {
        reasons.push(format!(
            "{risk} risk allele{} at rs10455872 and rs3798220",
            if risk != 1 { "s" } else { "" }
        ));
    }
    if class == "low" {
        reasons.push("a low KIV-2 repeat count".to_owned());
    }
    let read = kv_dot(&kv, "status") == "estimated" || any_read;
    let flag_text = if !reasons.is_empty() {
        format!(
            "Genetic signs of a high Lp(a): {}. Measure Lp(a) in blood.",
            reasons.join(" and ")
        )
    } else if read {
        "No genetic sign of a high Lp(a) from these markers. They explain only part of it: the KIV-2 count explained 21% to 27% of the variation in a Danish population study.".to_owned()
    } else {
        "Not determined: neither the KIV-2 count nor the two risk variants could be read. Measure Lp(a) in blood.".to_owned()
    };
    data.put("variants", Value::Array(variants));
    data.put("risk_alleles", Value::Int(risk));
    data.put("flag", s(if reasons.is_empty() { "." } else { "high" }));
    data.put("flag_text", s(&flag_text));
    Ok(Some((p, data.val())))
}

// --- step 52 --------------------------------------------------------------------

const VNTYPER_STATUSES: [&str; 4] = ["positive", "possible", "negative", "not assessed"];

pub fn vntyper(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("vntyper/{}_vntyper.tsv", c.s())) else {
        return Ok(None);
    };
    let kv = metric_map(&table(&p.path, false)?, "metric", "value");
    if !kv_get(&kv, "status").is_some_and(|st| VNTYPER_STATUSES.contains(&st)) {
        return unreadable(format!(
            "{} has no status row ({})",
            base(&p.path),
            VNTYPER_STATUSES.join(", ")
        ));
    }
    let v = [
        "status",
        "confidence",
        "variant",
        "flag",
        "alt_depth",
        "vntr_mean_depth",
        "coverage_qc",
        "sensitivity",
        "message",
        "advice",
    ]
    .iter()
    .fold(Obj::new(), |o, k| o.set(k, s(&kv_dot(&kv, k))))
    .val();
    Ok(Some((p, v)))
}

// --- step 53 --------------------------------------------------------------------

const LOCITYPER_COLUMNS: [&str; 13] = [
    "locus",
    "region",
    "status",
    "genotype",
    "haplotype_1",
    "haplotype_2",
    "quality",
    "total_reads",
    "unexplained_reads",
    "weighted_distance",
    "warnings",
    "filter",
    "reason",
];

pub fn locityper(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("locityper/{}_locityper.tsv", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    if !t.first_has(&LOCITYPER_COLUMNS) {
        return unreadable(format!(
            "{} lacks a row or one of the columns {}",
            base(&p.path),
            LOCITYPER_COLUMNS.join(", ")
        ));
    }
    let loci: Vec<Value> = t.rows.iter().map(|r| dots(r, &LOCITYPER_COLUMNS)).collect();
    let n = |k: &str, want: &str| {
        t.rows
            .iter()
            .filter(|r| or_dot(r, k) == want)
            .count()
    };
    let v = Obj::new()
        .set("called", count(n("status", "called")))
        .set("passed", count(n("filter", "pass")))
        .val();
    let mut out = Obj::new().set("loci", Value::Array(loci));
    if let Value::Object(m) = v {
        out.0.extend(m);
    }
    Ok(Some((p, out.val())))
}

fn text_at<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(k),
        _ => None,
    }
}

/// Each Locityper locus's primary call (`locityper_primary`), as text.
pub fn locityper_post(v: &mut Value, hla_state: State, hla: Option<&Value>, others: &Others) {
    let old = |st: State| {
        if st == State::Stale {
            " (stale: from an earlier run)"
        } else {
            ""
        }
    };
    let hla = hla.filter(|_| matches!(hla_state, State::Ok | State::Stale));
    let mut by_gene: IndexMap<String, &Value> = IndexMap::new();
    if let Some(Value::Array(loci)) = hla.and_then(|h| text_at(h, "loci")) {
        for x in loci {
            if let Some(Value::Str(g)) = text_at(x, "gene") {
                by_gene.insert(g.clone(), x);
            }
        }
    }
    let mut primary: IndexMap<&str, String> = IndexMap::new();
    for g in ["HLA-DRB1", "HLA-DQB1"] {
        let x = by_gene
            .get(g)
            .or_else(|| by_gene.get(&g.replace("HLA-", "")));
        primary.insert(
            g,
            match x {
                None => "T1K (step 08) did not run".to_owned(),
                Some(x) => {
                    let alleles = match text_at(x, "alleles") {
                        Some(Value::Array(a)) => a
                            .iter()
                            .map(|a| pycompat::value::py_str(a).unwrap_or_default())
                            .collect::<Vec<_>>()
                            .join("/"),
                        _ => String::new(),
                    };
                    let low = text_at(x, "low_confidence").is_some_and(pycompat::value::truthy);
                    format!(
                        "T1K (step 08): {}{}{}",
                        if alleles.is_empty() {
                            "no allele"
                        } else {
                            &alleles
                        },
                        if low { " (low confidence)" } else { "" },
                        old(hla_state)
                    )
                }
            },
        );
    }
    let typed = hla
        .and_then(|h| text_at(h, "kir"))
        .and_then(|k| text_at(k, "typed"))
        .filter(|t| pycompat::value::truthy(t));
    primary.insert(
        "KIR",
        match typed {
            Some(Value::Array(t)) => format!(
                "T1K KIR pass (step 08): {} genes typed{}",
                t.len(),
                old(hla_state)
            ),
            Some(Value::Object(t)) => format!(
                "T1K KIR pass (step 08): {} genes typed{}",
                t.len(),
                old(hla_state)
            ),
            _ => "T1K's KIR pass (step 08, KIR=true) gave no call".to_owned(),
        },
    );
    let mut cyp = Vec::new();
    if let Some((st, call)) = &others.pypgx {
        cyp.push(format!("pypgx (step 32): {call}{}", old(*st)));
    }
    if let Some((st, genotype, filter)) = &others.cyrius {
        let usable = if filter.as_deref() == Some("PASS") {
            String::new()
        } else {
            format!(
                " (not usable: Filter {})",
                filter.as_deref().filter(|f| !f.is_empty()).unwrap_or("missing")
            )
        };
        cyp.push(format!(
            "Cyrius (step 21): {}{usable}{}",
            genotype.as_deref().unwrap_or("none"),
            old(*st)
        ));
    }
    primary.insert(
        "CYP2D6",
        if cyp.is_empty() {
            "pypgx (step 32) and Cyrius (step 21) did not run".to_owned()
        } else {
            cyp.join("; ")
        },
    );
    primary.insert("FCGR", "no other step calls FCGR".to_owned());
    let Value::Object(m) = v else { return };
    let Some(Value::Array(loci)) = m.get_mut("loci") else {
        return;
    };
    for x in loci {
        let Value::Object(x) = x else { continue };
        let locus = match x.get("locus") {
            Some(Value::Str(l)) => l.clone(),
            _ => String::new(),
        };
        let p = primary
            .get(locus.as_str())
            .cloned()
            .unwrap_or_else(|| "no other step calls this locus".to_owned());
        x.insert("primary".into(), s(&p));
    }
}

// --- step 47 --------------------------------------------------------------------

pub fn viral(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("viral/{}_viral.tsv", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    let need = [
        "virus",
        "accession",
        "reads",
        "breadth",
        "relative_depth",
        "copies_per_cell",
        "call",
    ];
    if !t.first_has(&need) {
        return unreadable(format!(
            "{} has no rows or lacks a column of {}",
            base(&p.path),
            need.join(", ")
        ));
    }
    let all: Vec<&str> = need.iter().copied().chain(["human_depth", "note"]).collect();
    let keep = t
        .rows
        .iter()
        .map(|r| dots_strict(r, &all))
        .collect::<Result<Vec<_>, _>>()?;
    let get = |v: &Value, k: &str| match text_at(v, k) {
        Some(Value::Str(x)) => x.clone(),
        _ => String::new(),
    };
    let ici = keep
        .iter()
        .find(|r| get(r, "call") == "inherited integrated HHV-6 likely")
        .map(|r| get(r, "virus"));
    let human = get(&keep[0], "human_depth");
    let hhv6 = match ici {
        Some(v) => format!("inherited integrated {v} likely"),
        None if human == "." => "not assessed (no human depth)".to_owned(),
        None => "no integrated HHV-6".to_owned(),
    };
    let found: Vec<Value> = keep
        .iter()
        .filter(|r| {
            let n = get(r, "reads");
            isdigit(&n) && crate::int(&n).is_ok_and(|n| n > 0)
        })
        .cloned()
        .collect();
    let v = Obj::new()
        .set("hhv6", s(&hhv6))
        .set("human_depth", s(&human))
        .set("found", Value::Array(found))
        .set("genomes", count(keep.len()))
        .val();
    Ok(Some((p, v)))
}

// --- step 49 --------------------------------------------------------------------

pub fn chip(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let (Some(p), Some(dp)) = (
        c.find1(format!("chip/{s_}_chip.tsv")),
        c.find1(format!("chip/{s_}_chip_depth.tsv")),
    ) else {
        return Ok(None);
    };
    let calls = table(&p.path, false)?;
    let depth = table(&dp.path, false)?;
    let need = [
        "gene",
        "protein_change",
        "vaf",
        "alt_reads",
        "depth",
        "status",
        "note",
    ];
    let head = first_line(&p.path)?;
    if need.iter().any(|k| !head.iter().any(|h| h == k)) {
        return unreadable(format!(
            "{} lacks a column of {}",
            base(&p.path),
            need.join(", ")
        ));
    }
    let dneed = ["gene", "mean_depth", "detection_limit_vaf", "status"];
    if !depth.first_has(&dneed) {
        return unreadable(format!(
            "{} has no gene, mean_depth, detection_limit_vaf and status columns, or no row",
            base(&dp.path)
        ));
    }
    let chip: Vec<Value> = calls
        .rows
        .iter()
        .filter(|r| r.get("status").flatten() == Some("chip"))
        .map(|r| dots(r, &need))
        .collect();
    let mut covered = Vec::new();
    for r in &depth.rows {
        if r.get("status").flatten() == Some("covered") {
            let x = r.get("mean_depth").flatten().ok_or_else(|| {
                py_error(
                    "TypeError",
                    "float() argument must be a string or a real number, not 'NoneType'",
                )
            })?;
            covered.push(pycompat::num::parse_float(x).map_err(|_| float_err(x))?);
        }
    }
    covered.sort_by(f64::total_cmp);
    let med = match covered.len() {
        0 => None,
        n if n % 2 == 1 => Some(covered[n / 2]),
        n => Some((covered[n / 2 - 1] + covered[n / 2]) / 2.0),
    };
    let limit = match med {
        None => Value::Null,
        Some(m) => og_chip::binom::detection_limit(m)
            .map_err(|e| py_error("ValueError", &e.to_string()))?
            .map_or(Value::Null, Value::Float),
    };
    let v = Obj::new()
        .set("calls", Value::Array(chip.clone()))
        .set("other", count(calls.rows.len() - chip.len()))
        .set(
            "genes",
            Value::Array(depth.rows.iter().map(|r| dots(r, &dneed)).collect()),
        )
        .set("median_depth", med.map_or(Value::Null, Value::Float))
        .set("median_limit", limit)
        .set(
            "not_covered",
            Value::Array(
                depth
                    .rows
                    .iter()
                    .filter(|r| r.get("status").flatten() != Some("covered"))
                    .map(|r| opt_s(r.get("gene").flatten()))
                    .collect(),
            ),
        )
        .val();
    Ok(Some((p, v)))
}

// --- step 50 --------------------------------------------------------------------

pub fn mei(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("mei/{}_mei_sites.tsv", c.s())) else {
        return Ok(None);
    };
    let sites = table(&p.path, false)?;
    let need = [
        "gene",
        "element",
        "condition",
        "inheritance",
        "status",
        "mean_depth",
        "call",
    ];
    if !sites.first_has(&need) {
        return unreadable(format!(
            "{} has no rows or lacks a column of {}",
            base(&p.path),
            need.join(", ")
        ));
    }
    let dir = p.path.rsplit_once('/').map_or("", |x| x.0);
    let cp = crate::join(dir, &format!("{}_mei.tsv", c.sample));
    if !crate::is_file(&cp) {
        return unreadable(format!(
            "{} is missing beside {}",
            base(&cp),
            base(&p.path)
        ));
    }
    let calls = table(&cp, false)?;
    if ["chrom", "pos", "element", "copies"]
        .iter()
        .any(|k| !calls.has(k))
    {
        return unreadable(format!(
            "{} lacks a column of chrom, pos, element, copies",
            base(&cp)
        ));
    }
    let mut counts: std::collections::BTreeMap<String, usize> = Default::default();
    for r in &calls.rows {
        let e = r.get("element").flatten().unwrap_or("null").to_owned();
        *counts.entry(e).or_default() += 1;
    }
    let v = Obj::new()
        .set("calls", count(calls.rows.len()))
        .set(
            "by_element",
            Value::Object(counts.into_iter().map(|(k, n)| (k, count(n))).collect()),
        )
        .set(
            "sites",
            Value::Array(sites.rows.iter().map(|r| dots(r, &need)).collect()),
        )
        .val();
    Ok(Some((p, v)))
}

// --- step 51 --------------------------------------------------------------------

const MOCHA_SUMMARY_KEYS: [&str; 12] = [
    "status",
    "message",
    "depth",
    "het_sites",
    "het_per_10mb",
    "chromosomes",
    "floor_read",
    "floor_min",
    "floor",
    "calls_mocha",
    "calls_reported",
    "phasing",
];
const MOCHA_CALL_KEYS: [&str; 11] = [
    "chrom",
    "start",
    "end",
    "length_mb",
    "arms",
    "type",
    "cell_fraction",
    "n_hets",
    "lod_baf_phase",
    "reported",
    "reason",
];

pub fn mosaic(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find1(format!("mosaic/{s_}_mocha_summary.tsv")) else {
        return Ok(None);
    };
    let kv = metric_map(&table(&p.path, false)?, "metric", "value");
    if !matches!(kv_get(&kv, "status"), Some("assessed" | "not assessed")) {
        return unreadable(format!(
            "{} has no status row (assessed or not assessed)",
            base(&p.path)
        ));
    }
    let mut data = MOCHA_SUMMARY_KEYS
        .iter()
        .fold(Obj::new(), |o, k| o.set(k, s(&kv_dot(&kv, k))));
    let Some(cf) = c.find1(format!("mosaic/{s_}_mocha.tsv")) else {
        return unreadable(format!(
            "{} has no {}_mocha.tsv beside it",
            base(&p.path),
            c.sample
        ));
    };
    let head = first_line(&cf.path)?;
    if MOCHA_CALL_KEYS.iter().any(|k| !head.iter().any(|h| h == k)) {
        return unreadable(format!(
            "{} lacks a column of {}",
            base(&cf.path),
            MOCHA_CALL_KEYS.join(", ")
        ));
    }
    let rows: Table = table(&cf.path, false)?;
    let reported = |r: &pycompat::csv::DictRow| r.get("reported").flatten() == Some("yes");
    let calls: Vec<Value> = rows
        .rows
        .iter()
        .filter(|r| reported(r))
        .map(|r| dots(r, &MOCHA_CALL_KEYS))
        .collect();
    let n = calls.len();
    data.put("calls", Value::Array(calls));
    data.put(
        "not_reported",
        count(rows.rows.iter().filter(|r| !reported(r)).count()),
    );
    let said = kv_dot(&kv, "calls_reported");
    if said != "." && said != n.to_string() {
        return unreadable(format!(
            "{} has {n} reported calls, {} says {said}",
            base(&cf.path),
            base(&p.path)
        ));
    }
    Ok(Some((p, data.val())))
}

// --- step 48 --------------------------------------------------------------------

pub fn exomiser(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("exomiser/{}_exomiser_top.tsv", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, true)?;
    let header = first_line(&p.path)?;
    let need = [
        "rank",
        "gene",
        "moi",
        "exomiser_score",
        "phenotype_score",
        "variant",
        "genotype",
        "consequence",
        "hgvs",
        "acmg",
        "disease",
    ];
    if need.iter().any(|k| !header.iter().any(|h| h == k)) {
        return unreadable(format!(
            "{} lacks a column of {}",
            base(&p.path),
            need.join(", ")
        ));
    }
    if t.rows.iter().any(|r| ragged(r, &need)) {
        return unreadable(format!(
            "{} has a line that does not have one field per column",
            base(&p.path)
        ));
    }
    let dir = p.path.rsplit_once('/').map_or("", |x| x.0);
    let tp = crate::join(dir, &format!("{}_exomiser_terms.tsv", c.sample));
    let mut terms = Vec::new();
    if crate::is_file(&tp) {
        for r in table(&tp, true)?.rows {
            terms.push(
                Obj::new()
                    .set("id", get_default(&r, "hpo_id", ""))
                    .set("label", get_default(&r, "label", ""))
                    .val(),
            );
        }
    }
    let cell = |r: &pycompat::csv::DictRow, k: &str| s(r.get(k).flatten().unwrap_or(""));
    let mut genes: Vec<(String, Obj, Vec<Value>)> = Vec::new();
    for r in &t.rows {
        let rank = r.get("rank").flatten().unwrap_or("").to_owned();
        if genes.last().is_none_or(|g| g.0 != rank) {
            let o = ["rank", "gene", "moi", "exomiser_score", "phenotype_score"]
                .iter()
                .fold(Obj::new(), |o, k| o.set(k, cell(r, k)));
            genes.push((rank, o, Vec::new()));
        }
        let v = ["variant", "genotype", "consequence", "hgvs", "acmg", "disease"]
            .iter()
            .fold(Obj::new(), |o, k| o.set(k, cell(r, k)))
            .val();
        if let Some(g) = genes.last_mut() {
            g.2.push(v);
        }
    }
    let v = Obj::new()
        .set("terms", Value::Array(terms))
        .set(
            "genes",
            Value::Array(
                genes
                    .into_iter()
                    .map(|(_, o, vs)| o.set("variants", Value::Array(vs)).val())
                    .collect(),
            ),
        )
        .val();
    Ok(Some((p, v)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base64_matches_python() {
        assert_eq!(b64(b""), "");
        assert_eq!(b64(b"f"), "Zg==");
        assert_eq!(b64(b"fo"), "Zm8=");
        assert_eq!(b64(b"foo"), "Zm9v");
        assert_eq!(b64(b"foob"), "Zm9vYg==");
    }

    #[test]
    fn quote_keeps_unreserved_only() {
        assert_eq!(quote("R-M269"), "R-M269");
        assert_eq!(quote("H5a7 x/y"), "H5a7%20x%2Fy");
        assert_eq!(quote("é"), "%C3%A9");
    }
}
