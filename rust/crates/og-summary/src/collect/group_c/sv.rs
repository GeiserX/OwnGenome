//! Structural variants: Manta's and Delly's record counts, CNVpytor's
//! calls with the significance and quality filters, step 22's consensus
//! count, and step 05's class 4 and 5 calls with their support from the
//! consensus.

use super::roh::masked_overlap;
use super::util::{Obj, at, base, cells, count, int, lines, number, s, strs, table, vcf_lines};
use super::{Ctx, Found, Out};
use pycompat::csv::DictRow;
use pycompat::json::Value;
use std::collections::HashMap;

/// A record count and the records that PASS.
fn count_vcf(path: &str) -> Result<(usize, usize), String> {
    let (mut total, mut pass) = (0, 0);
    vcf_lines(path, |l| {
        total += 1;
        if l.split(|&c| c == b'\t').nth(6) == Some(b"PASS") {
            pass += 1;
        }
    })?;
    Ok((total, pass))
}

fn totals(p: Found) -> Out {
    let (t, ps) = count_vcf(&p.path)?;
    Ok(Some((
        p,
        Obj::new()
            .set("total", count(t))
            .set("pass", count(ps))
            .val(),
    )))
}

/// Step 04 (Manta).
pub fn manta(c: &Ctx<'_>) -> Out {
    match c.find(&[
        "manta/results/variants/diploidSV.vcf.gz".to_owned(),
        "manta/*diploidSV.vcf.gz".to_owned(),
        "manta/**/diploidSV.vcf.gz".to_owned(),
    ]) {
        Some(p) => totals(p),
        None => Ok(None),
    }
}

/// Step 19 (Delly).
pub fn delly(c: &Ctx<'_>) -> Out {
    match c.find1(format!("delly/{}_sv.vcf.gz", c.s())) {
        Some(p) => totals(p),
        None => Ok(None),
    }
}

/// Step 22's consensus count.
pub fn sv_consensus(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("sv_merged/{}_sv_consensus.vcf.gz", c.s())) else {
        return Ok(None);
    };
    let (t, _) = count_vcf(&p.path)?;
    Ok(Some((p, Obj::new().set("consensus", count(t)).val())))
}

/// A CNVpytor call is kept in the filtered count when fewer than half of
/// its reads have mapping quality 0 and fewer than half of its bases are N,
/// and it touches no masked band.
const CNV_MAX_Q0: f64 = 0.5;
const CNV_MAX_PN: f64 = 0.5;

/// `(\S+):(\d+)-(\d+)` over a whole cell.
fn span(x: &str) -> Option<(&str, i128, i128)> {
    let (c, rest) = x.rsplit_once(':')?;
    let (a, b) = rest.split_once('-')?;
    let digits = |v: &str| !v.is_empty() && v.bytes().all(|d| d.is_ascii_digit());
    if c.is_empty() || c.chars().any(char::is_whitespace) || !digits(a) || !digits(b) {
        return None;
    }
    Some((c, a.parse().ok()?, b.parse().ok()?))
}

/// Step 18 (CNVpytor): the calls, the significant ones (e-value under
/// 0.01), and those kept after the quality and band filters.
pub fn cnvpytor(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("cnvpytor/{}_cnvs.txt", c.s())) else {
        return Ok(None);
    };
    let (mut total, mut dels, mut dups, mut sig, mut kept, mut kept_y) = (0, 0, 0, 0, 0, 0);
    let (mut q0s, mut gaps, mut bands) = (0, 0, 0);
    for line in lines(&p.path)? {
        let c = pycompat::text::split(&line, None, None);
        if c.is_empty() {
            continue;
        }
        total += 1;
        dels += usize::from(c[0] == "deletion");
        dups += usize::from(c[0] == "duplication");
        let m = if c.len() > 4 { span(c[1]) } else { None };
        match c.get(4).map(|v| pycompat::num::parse_float(v)) {
            Some(Ok(e)) if e.partial_cmp(&0.01).is_none_or(|o| o.is_lt()) => {}
            _ => continue,
        }
        sig += 1;
        let (q0, pn) = if c.len() > 9 {
            (number(c[8]), number(c[9]))
        } else {
            (None, None)
        };
        if q0.is_some_and(|q| q >= CNV_MAX_Q0) {
            q0s += 1;
        } else if pn.is_some_and(|q| q >= CNV_MAX_PN) {
            gaps += 1;
        } else if let Some((ch, a, b)) = m
            && masked_overlap(ch, a, b) > 0
        {
            bands += 1;
        } else {
            kept += 1;
            kept_y += usize::from(m.is_some_and(|(ch, ..)| ch == "chrY" || ch == "Y"));
        }
    }
    let data = Obj::new()
        .set("total", count(total))
        .set("deletions", count(dels))
        .set("duplications", count(dups))
        .set("significant", count(sig))
        .set("significant_filtered", count(kept))
        .set("filtered_chry", count(kept_y))
        .set(
            "dropped",
            Obj::new()
                .set("q0", count(q0s))
                .set("gap", count(gaps))
                .set("band", count(bands))
                .val(),
        );
    Ok(Some((p, data.val())))
}

/// Takes the chrY calls out of the filtered count of a declared female;
/// with no declared sex they stay, and the data says so.
pub fn cnvpytor_sex_filter(v: &mut Value, declared_sex: Option<&str>) {
    let Value::Object(m) = v else { return };
    let (Some(Value::Int(kept)), Some(Value::Int(y))) = (
        m.get("significant_filtered").cloned(),
        m.get("filtered_chry").cloned(),
    ) else {
        return;
    };
    let declared = declared_sex.filter(|x| !x.is_empty());
    m.insert("declared_sex".into(), declared.map_or(Value::Null, s));
    if declared == Some("female") {
        if let Some(Value::Object(d)) = m.get_mut("dropped") {
            d.insert("chry_female".into(), Value::Int(y));
        }
        m.insert("significant_filtered".into(), Value::Int(kept - y));
    }
}

// --- step 05, AnnotSV class 4 and 5 ---------------------------------------------

/// Two breakpoints within this many bases of a consensus record's are the
/// same SV.
const CONSENSUS_BP: i128 = 1000;

/// AnnotSV writes chromosomes without `chr`; step 22 writes them with it.
fn chrom_key(c: &str) -> &str {
    if c.get(..3).is_some_and(|p| p.eq_ignore_ascii_case("chr")) {
        &c[3..]
    } else {
        c
    }
}

fn info_value<'a>(info: &'a str, key: &str) -> &'a str {
    info.split(';')
        .find_map(|kv| {
            let (k, v) = kv.split_once('=').unwrap_or((kv, ""));
            (k == key).then_some(v)
        })
        .unwrap_or("")
}

fn kind(t: &str) -> String {
    let t = pycompat::text::upper(t);
    if t.starts_with("DUP") {
        "DUP".into()
    } else {
        t
    }
}

struct Rec {
    pos: i128,
    end: i128,
    supp: i128,
    vec: String,
}

/// Step 22's consensus: caller names, and the records by (chrom, type).
/// The header and the records of the consensus VCF by (chrom, id).
type Consensus = (Vec<String>, HashMap<(String, String), Vec<Rec>>);

fn read_consensus(path: &str) -> Result<Consensus, String> {
    let mut names = Vec::new();
    let mut by: HashMap<(String, String), Vec<Rec>> = HashMap::new();
    for line in lines(path)? {
        if line.starts_with("##") {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if line.starts_with("#CHROM") {
            names = c.iter().skip(9).map(|x| (*x).to_owned()).collect();
            continue;
        }
        if c.len() < 8 {
            continue;
        }
        let Ok(pos) = crate::int(c[1]) else { continue };
        let end = match info_value(c[7], "END") {
            "" => Ok(pos),
            v => crate::int(v),
        };
        let supp = match info_value(c[7], "SUPP") {
            "" => Ok(0),
            v => crate::int(v),
        };
        let (Ok(end), Ok(supp)) = (end, supp) else {
            continue;
        };
        by.entry((chrom_key(c[0]).to_owned(), kind(info_value(c[7], "SVTYPE"))))
            .or_default()
            .push(Rec {
                pos,
                end,
                supp,
                vec: info_value(c[7], "SUPP_VEC").to_owned(),
            });
    }
    Ok((names, by))
}

/// Support and the supporting callers of each row from the consensus: the
/// record of the same type whose two breakpoints lie within 1,000 bp (the
/// nearest, the first of a tie), else 1 and the row's own caller; `None`
/// for both without a consensus file.
fn consensus_support(rows: &[DictRow], path: Option<&str>) -> Result<Vec<(Value, Value)>, String> {
    let Some(path) = path.filter(|p| crate::is_file(p)) else {
        return Ok(rows.iter().map(|_| (Value::Null, Value::Null)).collect());
    };
    let (names, by) = read_consensus(path)?;
    let mut out = Vec::new();
    for r in rows {
        let cell = |k: &str| r.get(k).flatten().unwrap_or_default();
        let mut best: Option<(i128, &Rec)> = None;
        if let (Ok(start), Ok(end)) = (crate::int(cell("SV_start")), crate::int(cell("SV_end"))) {
            let key = (
                chrom_key(cell("SV_chrom")).to_owned(),
                kind(cell("SV_type")),
            );
            for rec in by.get(&key).map_or(&[][..], Vec::as_slice) {
                let (dp, de) = ((rec.pos - start).abs(), (rec.end - end).abs());
                if dp <= CONSENSUS_BP && de <= CONSENSUS_BP && best.is_none_or(|(d, _)| dp + de < d)
                {
                    best = Some((dp + de, rec));
                }
            }
        }
        out.push(match best {
            Some((_, rec)) => {
                let callers: Vec<String> = rec
                    .vec
                    .chars()
                    .enumerate()
                    .filter(|&(_, x)| x == '1')
                    .map(|(i, _)| {
                        names
                            .get(i)
                            .cloned()
                            .unwrap_or_else(|| format!("caller {}", i + 1))
                    })
                    .collect();
                (int(rec.supp), s(&callers.join(",")))
            }
            None => (int(1), s(&super::util::or(r, "Caller", "."))),
        });
    }
    Ok(out)
}

/// `{FORMAT key: value}` of the sample column (the column after FORMAT,
/// unless it is `Annotation_mode`).
fn sample_format(r: &DictRow, head: &[String]) -> HashMap<String, String> {
    let keys: Vec<&str> = r.iter().map(|(k, _)| k).collect();
    let _ = head;
    let Some(i) = keys.iter().position(|k| *k == "FORMAT").map(|i| i + 1) else {
        return HashMap::new();
    };
    let Some(col) = keys.get(i).filter(|k| **k != "Annotation_mode") else {
        return HashMap::new();
    };
    let fmt = r.get("FORMAT").flatten().unwrap_or_default();
    let val = r.get(col).flatten().unwrap_or_default();
    fmt.split(':')
        .zip(val.split(':'))
        .map(|(k, v)| (k.to_owned(), v.to_owned()))
        .collect()
}

/// Step 05's class 4 and 5 SVs, each with the callers that support it from
/// step 22's consensus; a call step 05's rules moved to class 3 listed
/// apart.
pub fn annotsv(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("annotsv/{}_sv_class4_5_checked.tsv", c.s())) else {
        return Ok(None);
    };
    let size = std::fs::metadata(&p.path).map_or(0, |m| m.len());
    let t = if size > 0 {
        table(&p.path, false)?
    } else {
        super::util::Table {
            head: Vec::new(),
            rows: Vec::new(),
        }
    };
    let need = [
        "SV_chrom",
        "SV_start",
        "SV_end",
        "SV_type",
        "ACMG_class",
        "Caller",
        "Exon_edge",
        "Population_AF",
        "Common",
    ];
    if !t.rows.is_empty() && !need.iter().all(|k| t.has(k)) {
        let mut sorted = need.to_vec();
        sorted.sort_unstable();
        return Err(format!(
            "{} lacks a column of [{}]",
            base(&p.path),
            sorted
                .iter()
                .map(|k| crate::repr(k))
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    let mut cons = c.find1(format!("sv_merged/{}_sv_consensus.vcf.gz", c.s()));
    let mut cons_note = Value::Null;
    let support = match consensus_support(&t.rows, cons.as_ref().map(|f| f.path.as_str())) {
        Ok(x) => x,
        Err(e) => {
            let f = cons.take().map(|f| f.path).unwrap_or_default();
            let why = e.split_once(": ").map_or(e.as_str(), |x| x.1);
            cons_note = s(&format!("{} could not be read: {why}", base(&f)));
            consensus_support(&t.rows, None)?
        }
    };
    let mut keep = Vec::new();
    let mut moved = Vec::new();
    for (r, (sp, callers)) in t.rows.iter().zip(support) {
        // A row cut short reads None there, as Python's DictReader gives it.
        let cell =
            |k: &str| -> Result<String, String> { Ok(at(r, k)?.unwrap_or("None").to_owned()) };
        let cellv = |k: &str| -> Result<Value, String> { Ok(super::util::opt_s(at(r, k)?)) };
        let af = match at(r, "Population_AF")? {
            None | Some("") => Value::Null,
            Some(v) => pycompat::num::parse_float(v).map_or(Value::Null, Value::Float),
        };
        let genes: Vec<String> = super::util::or(r, "Gene_name", "")
            .split(';')
            .filter(|g| !g.is_empty())
            .map(str::to_owned)
            .collect();
        let fmt = sample_format(r, &t.head);
        let fmt_get = |k: &str| {
            fmt.get(k)
                .filter(|v| !v.is_empty())
                .map_or(Value::Null, |v| s(v))
        };
        let acmg = cell("ACMG_class")?;
        let support_callers = match &callers {
            Value::Str(x) if !x.is_empty() => strs(x.split(',')),
            _ => Value::Array(Vec::new()),
        };
        let row = Obj::new()
            .set(
                "sv",
                s(&format!(
                    "{}:{}-{}",
                    cell("SV_chrom")?,
                    cell("SV_start")?,
                    cell("SV_end")?
                )),
            )
            .set("type", cellv("SV_type")?)
            .set("gt", fmt_get("GT"))
            .set("ps", fmt_get("PS"))
            .set("length", s(&super::util::or(r, "SV_length", "")))
            .set("acmg_class", cellv("ACMG_class")?)
            .set("genes", strs(&genes))
            .set("caller", cellv("Caller")?)
            .set("support", sp)
            .set("support_callers", support_callers)
            .set("exon_edge", cellv("Exon_edge")?)
            .set(
                "exon_edge_detail",
                s(&super::util::or(r, "Exon_edge_detail", "")),
            )
            .set("population_af", af)
            .set(
                "population_source",
                s(&super::util::or(r, "Population_source", "")),
            )
            .set("common", Value::Bool(at(r, "Common")? == Some("yes")))
            .set(
                "annotsv_class",
                match r.get("AnnotSV_class").flatten().filter(|v| !v.is_empty()) {
                    Some(v) => s(v),
                    None => cellv("ACMG_class")?,
                },
            )
            .set("rule", s(&super::util::or(r, "Rule", "")))
            .set("rule_detail", s(&super::util::or(r, "Rule_detail", "")));
        if matches!(acmg.as_str(), "4" | "5") {
            keep.push(row.val());
        } else {
            moved.push(row.val());
        }
    }
    let field = |v: &Value, k: &str| -> Value {
        match v {
            Value::Object(m) => m.get(k).cloned().unwrap_or(Value::Null),
            _ => Value::Null,
        }
    };
    let single = keep
        .iter()
        .filter(|x| field(x, "support") == Value::Int(1))
        .count();
    let edge = keep
        .iter()
        .filter(|x| field(x, "exon_edge") == s("yes"))
        .count();
    let edge_checked = !keep
        .iter()
        .any(|x| field(x, "exon_edge") == s("not checked"));
    let common = keep
        .iter()
        .filter(|x| field(x, "common") == Value::Bool(true))
        .count();
    let data = Obj::new()
        .set("count", count(keep.len()))
        .set("rows", Value::Array(keep))
        .set("downgraded", Value::Array(moved))
        .set("consensus", cons.map_or(Value::Null, |f| s(&f.rel)))
        .set("consensus_note", cons_note)
        .set("single_caller", count(single))
        .set("exon_edge", count(edge))
        .set("exon_edge_checked", Value::Bool(edge_checked))
        .set("common", count(common));
    let _ = cells;
    Ok(Some((p, data.val())))
}
