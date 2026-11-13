//! Step 20's chrM section: the PASS calls and the heteroplasmic ones, the
//! poly-C tracts counted apart, the annotation table with its ClinVar
//! classes, the confirmed variants of `assets/mito_confirmed.tsv`, the large
//! deletions and the mtDNA copy number.

use super::util::{Obj, asset_rows, cells, lines, opt_s, pairs, round, s, strs, vcf_lines};
use super::{Ctx, Found, Out};
use crate::collect::sections::py_error;
use indexmap::IndexMap;
use pycompat::csv::{DictReader, DictRow, Quoting};
use pycompat::json::Value;

/// The D-loop poly-C tracts, rCRS m.303-315 and m.16184-16193.
const POLYC_TRACTS: [(i128, i128); 2] = [(303, 315), (16184, 16193)];
const MITO_ANNOT_COLS: [&str; 13] = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "filter",
    "af",
    "clinvar_id",
    "clinvar_significance",
    "clinvar_review",
    "clinvar_condition",
    "gnomad_af_hom",
    "gnomad_af_het",
    "gnomad_max_hl",
];
const CLINVAR_PLP: [&str; 3] = [
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
];
const CLINVAR_BENIGN: [&str; 4] = [
    "Benign",
    "Likely_benign",
    "Benign/Likely_benign",
    "not_provided",
];
const MITO_DEL_TOLERANCE: i128 = 50;
const MITO_DEL_LIST_PCT: f64 = 1.0;
const MITO_CN_COLS: [&str; 8] = [
    "status",
    "reason",
    "chrm_reads",
    "chrm_length",
    "autosomal_reads",
    "autosomal_length",
    "depth_ratio",
    "mtdna_copies",
];
const HET_LOW: f64 = 0.05;
const HET_HIGH: f64 = 0.95;

fn n_chars(x: &str) -> i128 {
    x.chars().count() as i128
}

/// The rCRS positions a chrM allele changes.
fn mito_span(pos: i128, r: &str, a: &str) -> (i128, i128) {
    if n_chars(r) == n_chars(a) {
        (pos, pos + n_chars(r) - 1)
    } else {
        (pos + 1, pos + (n_chars(r) - 1).max(1))
    }
}

pub fn in_polyc(pos: i128, r: &str, a: &str) -> bool {
    let (lo, hi) = mito_span(pos, r, a);
    POLYC_TRACTS.iter().any(|(a, b)| lo <= *b && hi >= *a)
}

/// `(pos, ref, alt)` without the bases both share at the end, then at the
/// start, one base always kept.
pub fn min_allele(mut pos: i128, r: &str, a: &str) -> (i128, String, String) {
    let mut r: Vec<char> = r.chars().collect();
    let mut a: Vec<char> = a.chars().collect();
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    let mut i = 0;
    while r.len() - i > 1 && a.len() - i > 1 && r[i] == a[i] {
        i += 1;
        pos += 1;
    }
    (pos, r[i..].iter().collect(), a[i..].iter().collect())
}

fn py_int(x: Option<&str>) -> Result<i128, String> {
    match x {
        None => Err(py_error(
            "TypeError",
            "int() argument must be a string, a bytes-like object or a real number, not 'NoneType'",
        )),
        Some(x) => crate::int(x).map_err(|e| py_error("ValueError", &e)),
    }
}

fn py_float(x: Option<&str>) -> Result<f64, String> {
    match x {
        None => Err(py_error(
            "TypeError",
            "float() argument must be a string or a real number, not 'NoneType'",
        )),
        Some(x) => crate::float(x).map_err(|e| py_error("ValueError", &e)),
    }
}

/// `row[k]`, the column required.
fn col<'a>(r: &'a DictRow, k: &str) -> Result<Option<&'a str>, String> {
    super::util::at(r, k)
}

fn confirmed_fields(k: &DictRow) -> Result<Obj, String> {
    let mut o = Obj::new();
    for (key, c) in [
        ("hgvs", "hgvs"),
        ("gene", "gene"),
        ("disease", "disease"),
        ("clinvar_id", "clinvar_id"),
        ("pmids", "pmids"),
        ("note", "heteroplasmy_note"),
    ] {
        o.put(key, opt_s(col(k, c)?));
    }
    Ok(o)
}

/// Python's `max` of floats: the first, replaced by each later one that is
/// greater.
fn py_max(xs: &[f64]) -> Option<f64> {
    let mut it = xs.iter().copied();
    let mut cur = it.next()?;
    for x in it {
        if x > cur {
            cur = x;
        }
    }
    Some(cur)
}

type Confirmed = (IndexMap<(i128, String, String), DictRow>, Vec<DictRow>);

fn mito_confirmed(assets: &str) -> Result<Result<Confirmed, String>, String> {
    let rows = match asset_rows(assets, "mito_confirmed.tsv") {
        Ok(r) => r,
        // A missing assets folder is a note, not an unreadable section.
        Err(e) if e.starts_with("mito_confirmed.tsv is not in ") => return Ok(Err(e)),
        Err(e) => return Err(e),
    };
    let mut snv = IndexMap::new();
    let mut dels = Vec::new();
    for r in rows {
        match col(&r, "kind")? {
            Some("snv") => {
                let key = (
                    py_int(col(&r, "pos")?)?,
                    col(&r, "ref")?.unwrap_or("").to_owned(),
                    col(&r, "alt")?.unwrap_or("").to_owned(),
                );
                snv.insert(key, r);
            }
            Some("deletion") => dels.push(r),
            _ => {}
        }
    }
    Ok(Ok((snv, dels)))
}

fn mito_extended(c: &Ctx<'_>, vcf: &Found) -> Result<Obj, String> {
    let mut out = Obj::new();
    let (snv, dels) = match mito_confirmed(c.assets)? {
        Ok((snv, dels)) => (Some(snv), dels),
        Err(note) => {
            out.put("confirmed_note", s(&note));
            (None, Vec::new())
        }
    };
    if let Some(snv) = &snv {
        let mut records = Vec::new();
        vcf_lines(&vcf.path, |l| records.push(cells(l)))?;
        let mut confirmed = Vec::new();
        for r in records {
            if r.len() < 8 {
                continue;
            }
            let (fmt, val): (Vec<&str>, Vec<&str>) = if r.len() > 9 {
                (r[8].split(':').collect(), r[9].split(':').collect())
            } else {
                (Vec::new(), Vec::new())
            };
            let afs: Vec<&str> = match fmt.iter().position(|f| *f == "AF") {
                Some(i) if i < val.len() => val[i].split(',').collect(),
                _ => Vec::new(),
            };
            let pos = py_int(Some(&r[1]))?;
            for (i, alt) in r[4].split(',').enumerate() {
                if let Some(k) = snv.get(&min_allele(pos, &r[3], alt)) {
                    confirmed.push(
                        confirmed_fields(k)?
                            .set("pos", Value::Int(py_int(col(k, "pos")?)?))
                            .set("filter", s(&r[6]))
                            .set("af", s(afs.get(i).copied().unwrap_or(".")))
                            .val(),
                    );
                }
            }
        }
        out.put("confirmed", Value::Array(confirmed));
    }
    let s_ = c.s();
    match c.find1(format!("mito/{s_}_chrM_deletions.tsv")) {
        None => out.put(
            "deletions_note",
            s("no deletion table (mito/<id>_chrM_deletions.tsv): step 20 ran before it looked for large deletions, eKLIPse wrote no table (see the step's log), or the Nextflow MITO_DELETIONS task failed (see the Nextflow log)"),
        ),
        Some(p) => {
            let rows = super::util::table(&p.path, false)?.rows;
            let mut listed = Vec::new();
            for r in &rows {
                let bkp5 = py_int(col(r, "bkp5")?)?;
                let bkp3 = py_int(col(r, "bkp3")?)?;
                let length = match col(r, "length")? {
                    Some(".") => Value::Null,
                    x => Value::Int(py_int(x)?),
                };
                let level = py_float(col(r, "level_pct")?)?;
                let mut known = None;
                for k in &dels {
                    let kp = py_int(col(k, "pos")?)?;
                    if (bkp5 + 1 - kp).abs() <= MITO_DEL_TOLERANCE {
                        let ke = py_int(col(k, "end")?)?;
                        if (bkp3 - 1 - ke).abs() <= MITO_DEL_TOLERANCE {
                            known = Some(k);
                            break;
                        }
                    }
                }
                if known.is_some() || level >= MITO_DEL_LIST_PCT {
                    let mut row = Obj::new()
                        .set("bkp5", Value::Int(bkp5))
                        .set("bkp3", Value::Int(bkp3))
                        .set("length", length)
                        .set("level_pct", Value::Float(level));
                    if let Some(k) = known {
                        row.put("known", confirmed_fields(k)?.val());
                    }
                    listed.push(row.val());
                }
            }
            out.put(
                "deletions",
                Obj::new()
                    .set("source", s(&p.rel))
                    .set("total", super::util::count(rows.len()))
                    .set("listed", Value::Array(listed))
                    .set("list_floor_pct", Value::Float(MITO_DEL_LIST_PCT))
                    .val(),
            );
        }
    }
    match c.find1(format!("mito/{s_}_mtdna_copies.tsv")) {
        Some(p) => {
            let m: IndexMap<String, String> = pairs(&p.path)?.into_iter().collect();
            out.put(
                "copies",
                MITO_CN_COLS
                    .iter()
                    .fold(Obj::new(), |o, k| {
                        o.set(k, opt_s(m.get(*k).map(String::as_str).filter(|v| *v != ".")))
                    })
                    .val(),
            );
        }
        None => out.put(
            "copies_note",
            s("no copy-number table (mito/<id>_mtdna_copies.tsv): step 20 ran before it counted mtDNA copies, or the Nextflow MITO_DELETIONS task failed (see the Nextflow log)"),
        ),
    }
    Ok(out)
}

fn mtime(p: &str) -> Result<f64, String> {
    pycompat::mtime::mtime(p).map_err(|e| super::util::os_note(&crate::os_error(&e, p)))
}

/// The annotation table's sources and its PASS rows, or a note.
fn mito_annotation(c: &Ctx<'_>, vcf: &Found) -> Result<(Obj, Vec<DictRow>), String> {
    let Some(p) = c.find1(format!("mito/{}_chrM_annotated.tsv", c.s())) else {
        return Ok((
            Obj::new().set("note", s("no annotation table (mito/<id>_chrM_annotated.tsv): step 20 ran before it annotated, or the Nextflow MITO_VARIANTS_ANNOTATE task failed (see the Nextflow log)")),
            Vec::new(),
        ));
    };
    if mtime(&p.path)? < mtime(&vcf.path)? {
        return Ok((
            Obj::new().set(
                "note",
                s(&format!(
                    "{} is older than the chrM calls: rerun step 20",
                    p.rel
                )),
            ),
            Vec::new(),
        ));
    }
    let mut sources: IndexMap<String, String> = IndexMap::new();
    let mut body = Vec::new();
    for line in lines(&p.path)? {
        if let Some(rest) = line.strip_prefix("# ") {
            let c: Vec<&str> = rest.split('\t').collect();
            if c.len() >= 3 {
                let v = if c[1] == "not annotated" {
                    format!("not annotated: {}", c[2])
                } else if c[2] != "-" {
                    format!("{} ({})", c[1], c[2])
                } else {
                    c[1].to_owned()
                };
                sources.insert(c[0].to_owned(), v);
            }
        } else if !pycompat::num::py_strip(&line).is_empty() {
            body.push(Ok::<_, std::io::Error>(line + "\n"));
        }
    }
    let rows = DictReader::new(body.into_iter(), '\t', Quoting::Minimal)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| py_error(e.exception(), &e.to_string()))?
        .into_iter()
        .filter(|r| r.get("filter").flatten() == Some("PASS"))
        .collect();
    let src = |k: &str| opt_s(sources.get(k).map(String::as_str));
    Ok((
        Obj::new()
            .set("source", s(&p.rel))
            .set("clinvar", src("clinvar"))
            .set("gnomad", src("gnomad")),
        rows,
    ))
}

fn clinvar_terms(sig: Option<&str>) -> Vec<&str> {
    sig.unwrap_or("")
        .split(['|', ','])
        .filter(|t| !t.is_empty() && *t != ".")
        .collect()
}

pub fn mito(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("mito/{}_chrM_filtered.vcf.gz", c.s())) else {
        return Ok(None);
    };
    let mut records = Vec::new();
    vcf_lines(&p.path, |l| records.push(cells(l)))?;
    let (mut passed, mut het) = (0usize, 0usize);
    let (mut polyc, mut sites) = (Vec::new(), Vec::new());
    for r in &records {
        if r.len() < 10 || r[6] != "PASS" {
            continue;
        }
        passed += 1;
        let fmt: Vec<&str> = r[8].split(':').collect();
        let val: Vec<&str> = r[9].split(':').collect();
        let Some(i) = fmt
            .iter()
            .position(|f| *f == "AF")
            .filter(|i| *i < val.len())
        else {
            continue;
        };
        let afs: Result<Vec<f64>, _> = val[i]
            .split(',')
            .filter(|x| !matches!(*x, "" | "."))
            .map(pycompat::num::parse_float)
            .collect();
        let Some(af) = afs.ok().and_then(|a| py_max(&a)) else {
            continue;
        };
        if !(HET_LOW..HET_HIGH).contains(&af) {
            continue;
        }
        let pos = py_int(Some(&r[1]))?;
        if r[4].split(',').any(|a| in_polyc(pos, &r[3], a)) {
            polyc.push(
                Obj::new()
                    .set(
                        "variant",
                        s(&format!("{}:{} {}>{}", r[0], r[1], r[3], r[4])),
                    )
                    .set("af", round(af, 3))
                    .val(),
            );
        } else {
            het += 1;
            sites.push(
                Obj::new()
                    .set("pos", Value::Int(pos))
                    .set("ref", s(&r[3]))
                    .set("alt", s(&r[4]))
                    .set("af", round(af, 3))
                    .val(),
            );
        }
    }
    let (mut ann, rows) = mito_annotation(c, &p)?;
    let mut listed = Vec::new();
    let mut plp = 0i128;
    for r in &rows {
        let terms = clinvar_terms(r.get("clinvar_significance").flatten());
        let is_plp = terms.iter().any(|t| CLINVAR_PLP.contains(t));
        plp += i128::from(is_plp);
        let het_row = r
            .get("af")
            .flatten()
            .and_then(|x| pycompat::num::parse_float(x).ok())
            .is_some_and(|af| (HET_LOW..HET_HIGH).contains(&af));
        if is_plp
            || het_row
            || (!terms.is_empty() && !terms.iter().all(|t| CLINVAR_BENIGN.contains(t)))
        {
            let mut o = Obj::new();
            for k in MITO_ANNOT_COLS {
                if k != "chrom" && k != "filter" {
                    o.put(k, opt_s(r.get(k).flatten()));
                }
            }
            let pos = py_int(col(r, "pos")?)?;
            let rf = col(r, "ref")?.unwrap_or("");
            let al = col(r, "alt")?.unwrap_or("");
            o.put("polyc", Value::Bool(in_polyc(pos, rf, al)));
            listed.push(o.val());
        }
    }
    if ann.0.contains_key("source") {
        ann.put("variants", Value::Array(listed));
        ann.put("clinvar_plp", Value::Int(plp));
    }
    let mut out = mito_extended(c, &p)?;
    out.put("pass", super::util::count(passed));
    out.put("heteroplasmic", super::util::count(het));
    out.put("heteroplasmic_sites", Value::Array(sites));
    out.put("heteroplasmy_floor", Value::Float(HET_LOW));
    out.put("polyc_heteroplasmic", Value::Array(polyc));
    out.put(
        "polyc_tracts",
        strs(POLYC_TRACTS.iter().map(|(a, b)| format!("m.{a}-{b}"))),
    );
    out.put("annotation", ann.val());
    Ok(Some((p, out.val())))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_allele_trims_end_then_start() {
        assert_eq!(min_allele(302, "AC", "ACC"), (302, "A".into(), "AC".into()));
        assert_eq!(min_allele(10, "ACT", "AGT"), (11, "C".into(), "G".into()));
        assert_eq!(min_allele(5, "A", "G"), (5, "A".into(), "G".into()));
    }

    #[test]
    fn polyc_span() {
        assert!(in_polyc(302, "A", "AC"));
        assert!(!in_polyc(300, "A", "G"));
        assert!(in_polyc(16189, "T", "C"));
    }

    #[test]
    fn python_max_keeps_a_leading_nan() {
        assert!(py_max(&[f64::NAN, 0.5]).is_some_and(f64::is_nan));
        assert_eq!(py_max(&[0.5, f64::NAN]), Some(0.5));
    }
}
