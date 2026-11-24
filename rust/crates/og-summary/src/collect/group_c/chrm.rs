//! Step 20's chrM section: the PASS calls and the heteroplasmic ones, the
//! poly-C tracts counted apart, the annotation table with its ClinVar
//! classes, the confirmed variants of `assets/mito_confirmed.tsv`, the large
//! deletions and the mtDNA copy number. A heteroplasmic call near the
//! origin of the circular reference, or under 10% AF, stays counted and
//! carries a flag that says why it is less reliable. When step 20 called the
//! control region again on the rotated chrM (the `##og_shifted_pass=` line
//! of its calls), no call sits next to a reference end, so none carries the
//! origin flag.

use super::util::{Obj, asset_rows, base, cells, lines, opt_s, pairs, round, s, strs, vcf_lines};
use super::{Ctx, Found, Out};
use crate::collect::sections::py_error;
use indexmap::IndexMap;
use pycompat::csv::{DictReader, DictRow, Quoting};
use pycompat::json::Value;

/// The D-loop poly-C tracts, rCRS m.303-315 and m.16184-16193.
const POLYC_TRACTS: [(i128, i128); 2] = [(303, 315), (16184, 16193)];
/// The two ends of the control region, next to the origin of the circular
/// reference: reads that span the origin align poorly at both edges, unless
/// the shifted-reference pass called them again ([`SHIFTED_MARK`]).
const ORIGIN_EDGE: [(i128, i128); 2] = [(1, 500), (16000, 16569)];
/// The header line `og vcf shift-back` (og-vcf's `shift::MARK`) writes into
/// the chrM calls it merges: the control region came from Mutect2 on chrM
/// rotated by 8000 bp, where the origin lies mid-sequence.
const SHIFTED_MARK: &str = "##og_shifted_pass=";
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
/// Under this AF a heteroplasmic call is low level: often somatic or
/// technical, and the range where NuMT reads pass as heteroplasmy.
const LOW_AF: f64 = 0.10;

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

fn overlaps(regions: &[(i128, i128)], pos: i128, r: &str, a: &str) -> bool {
    let (lo, hi) = mito_span(pos, r, a);
    regions.iter().any(|(a, b)| lo <= *b && hi >= *a)
}

pub fn in_polyc(pos: i128, r: &str, a: &str) -> bool {
    overlaps(&POLYC_TRACTS, pos, r, a)
}

/// Whether the bases an allele changes touch either end of the control
/// region (`ORIGIN_EDGE`). The allele is trimmed first (`min_allele`), so a
/// padded multi-allelic record reads the same as its minimal form.
pub fn at_origin_edge(pos: i128, r: &str, a: &str) -> bool {
    let (pos, r, a) = min_allele(pos, r, a);
    overlaps(&ORIGIN_EDGE, pos, &r, &a)
}

fn low_af(af: f64) -> bool {
    af < LOW_AF
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

/// The value of the `##og_shifted_pass=` line of the chrM calls, if the
/// header has one.
fn shifted_pass(path: &str) -> Result<Option<String>, String> {
    use std::io::BufRead;
    let fail = |e: std::io::Error| {
        format!(
            "{}: {}",
            base(path),
            crate::collect::sections::gzip_message(&e)
        )
    };
    let r = og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(fail)?;
    for line in r.split(b'\n') {
        let line = line.map_err(fail)?;
        if line.first() != Some(&b'#') {
            break;
        }
        if let Some(v) = line.strip_prefix(SHIFTED_MARK.as_bytes()) {
            let v = v.strip_suffix(b"\r").unwrap_or(v);
            return Ok(Some(String::from_utf8_lossy(v).into_owned()));
        }
    }
    Ok(None)
}

pub fn mito(c: &Ctx<'_>) -> Out {
    let Some(p) = c.find1(format!("mito/{}_chrM_filtered.vcf.gz", c.s())) else {
        return Ok(None);
    };
    let mut records = Vec::new();
    vcf_lines(&p.path, |l| records.push(cells(l)))?;
    let shifted = shifted_pass(&p.path)?;
    // After the shifted pass no call sits next to a reference end.
    let edge = |pos: i128, r: &str, a: &str| shifted.is_none() && at_origin_edge(pos, r, a);
    let (mut passed, mut het, mut het_edge, mut het_low) = (0usize, 0usize, 0usize, 0usize);
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
            // Flagged calls stay counted: the flags say why one is less
            // reliable, they never hide it.
            let at_edge = r[4].split(',').any(|a| edge(pos, &r[3], a));
            let low = low_af(af);
            het += 1;
            het_edge += usize::from(at_edge);
            het_low += usize::from(low);
            sites.push(
                Obj::new()
                    .set("pos", Value::Int(pos))
                    .set("ref", s(&r[3]))
                    .set("alt", s(&r[4]))
                    .set("af", round(af, 3))
                    .set("origin_edge", Value::Bool(at_edge))
                    .set("low_af", Value::Bool(low))
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
        let af = r
            .get("af")
            .flatten()
            .and_then(|x| pycompat::num::parse_float(x).ok());
        let het_row = af.is_some_and(|af| (HET_LOW..HET_HIGH).contains(&af));
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
            o.put("origin_edge", Value::Bool(edge(pos, rf, al)));
            // Low level only when heteroplasmic: a homoplasmic or absent AF
            // is not a low-level call.
            o.put("low_af", Value::Bool(het_row && af.is_some_and(low_af)));
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
    out.put("heteroplasmic_origin_edge", super::util::count(het_edge));
    out.put("heteroplasmic_low_af", super::util::count(het_low));
    out.put(
        "origin_edge_regions",
        strs(ORIGIN_EDGE.iter().map(|(a, b)| format!("m.{a}-{b}"))),
    );
    out.put("low_af_ceiling", Value::Float(LOW_AF));
    // Only when the pass ran, so a summary without it reads as before.
    if let Some(v) = shifted {
        out.put("shifted_pass", s(&v));
    }
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
        // The tract's last base is in it, the next base is not; an SNV just
        // before a tract is not, an insertion anchored there is.
        assert!(in_polyc(315, "C", "T"));
        assert!(!in_polyc(316, "C", "T"));
        assert!(!in_polyc(302, "A", "G"));
    }

    #[test]
    fn python_max_keeps_a_leading_nan() {
        assert!(py_max(&[f64::NAN, 0.5]).is_some_and(f64::is_nan));
        assert_eq!(py_max(&[0.5, f64::NAN]), Some(0.5));
    }

    const HEAD: &str =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
    const ANNOT_HEAD: &str = "# clinvar\t2026-10-01\t-\n# gnomad\tv3.1\t-\n\
        chrom\tpos\tref\talt\tfilter\taf\tclinvar_id\tclinvar_significance\tclinvar_review\t\
        clinvar_condition\tgnomad_af_hom\tgnomad_af_het\tgnomad_max_hl\n";

    /// Step 20's section of a sample folder holding these chrM calls and
    /// annotation rows, read with the repository's `assets/`.
    fn section(vcf_body: &str, annot_body: &str) -> Value {
        section_with(HEAD, vcf_body, annot_body)
    }

    /// [`section`] with the VCF header `head`.
    fn section_with(head: &str, vcf_body: &str, annot_body: &str) -> Value {
        let d = crate::testutil::Dir::new("chrm-flags");
        d.put_gz("mito/S_chrM_filtered.vcf.gz", &format!("{head}{vcf_body}"));
        d.put(
            "mito/S_chrM_annotated.tsv",
            format!("{ANNOT_HEAD}{annot_body}"),
        );
        let index = crate::collect::scan::FolderIndex::walk(&d.0);
        let dir = d.0.display().to_string();
        let assets = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets");
        let c = Ctx {
            sample: "S",
            dir: &dir,
            assets,
            tables: "",
            index: &index,
        };
        mito(&c).unwrap().expect("a mito section").1
    }

    fn at<'v>(v: &'v Value, path: &[&str]) -> &'v Value {
        path.iter().fold(v, |v, k| match v {
            Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no key {k} in {v:?}")),
            _ => panic!("{k}: not an object: {v:?}"),
        })
    }

    fn items(v: &Value) -> &[Value] {
        match v {
            Value::Array(a) => a,
            _ => panic!("not an array: {v:?}"),
        }
    }

    #[test]
    fn origin_edge_call_is_flagged_and_stays_counted() {
        let d = section(
            "chrM\t16\t.\tA\tT\t50\tPASS\t.\tGT:AF\t0/1:0.12\n\
             chrM\t16050\t.\tC\tT\t50\tPASS\t.\tGT:AF\t0/1:0.30\n\
             chrM\t8000\t.\tC\tT\t50\tPASS\t.\tGT:AF\t0/1:0.40\n",
            "chrM\t16\tA\tT\tPASS\t0.12\t5\tPathogenic\t.\t.\t.\t.\t.\n",
        );
        assert_eq!(at(&d, &["heteroplasmic"]), &Value::Int(3));
        assert_eq!(at(&d, &["heteroplasmic_origin_edge"]), &Value::Int(2));
        assert_eq!(at(&d, &["heteroplasmic_low_af"]), &Value::Int(0));
        let sites = items(at(&d, &["heteroplasmic_sites"]));
        let flags: Vec<(&Value, &Value, &Value)> = sites
            .iter()
            .map(|x| (at(x, &["pos"]), at(x, &["origin_edge"]), at(x, &["low_af"])))
            .collect();
        let (t, f) = (Value::Bool(true), Value::Bool(false));
        assert_eq!(
            flags,
            [
                (&Value::Int(16), &t, &f),
                (&Value::Int(16050), &t, &f),
                (&Value::Int(8000), &f, &f)
            ]
        );
        // The ClinVar class still counts; the listed call carries the flag.
        assert_eq!(at(&d, &["annotation", "clinvar_plp"]), &Value::Int(1));
        let v = &items(at(&d, &["annotation", "variants"]))[0];
        assert_eq!(at(v, &["origin_edge"]), &t);
        assert_eq!(at(v, &["low_af"]), &f);
        assert_eq!(
            at(&d, &["origin_edge_regions"]),
            &strs(["m.1-500", "m.16000-16569"])
        );
        assert_eq!(at(&d, &["low_af_ceiling"]), &Value::Float(0.10));
    }

    #[test]
    fn after_the_shifted_pass_no_call_is_flagged_at_the_origin() {
        let body = "chrM\t16\t.\tA\tT\t50\tPASS\t.\tGT:AF\t0/1:0.12\n\
             chrM\t16050\t.\tC\tT\t50\tPASS\t.\tGT:AF\t0/1:0.30\n";
        let annot = "chrM\t16\tA\tT\tPASS\t0.12\t5\tPathogenic\t.\t.\t.\t.\t.\n";
        let mark = "chrM:16025-16569 and chrM:1-575 from the rotated chrM";
        let head = format!(
            "##fileformat=VCFv4.2\n{SHIFTED_MARK}{mark}\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        );
        let f = Value::Bool(false);
        let d = section_with(&head, body, annot);
        assert_eq!(at(&d, &["heteroplasmic"]), &Value::Int(2));
        assert_eq!(at(&d, &["heteroplasmic_origin_edge"]), &Value::Int(0));
        for site in items(at(&d, &["heteroplasmic_sites"])) {
            assert_eq!(at(site, &["origin_edge"]), &f, "{site:?}");
        }
        let v = &items(at(&d, &["annotation", "variants"]))[0];
        assert_eq!(at(v, &["origin_edge"]), &f);
        assert_eq!(at(&d, &["shifted_pass"]), &s(mark));
        // The same calls without the line: both flagged, and no new key.
        let d = section(body, annot);
        assert_eq!(at(&d, &["heteroplasmic_origin_edge"]), &Value::Int(2));
        match &d {
            Value::Object(m) => assert!(!m.contains_key("shifted_pass")),
            _ => panic!("not an object"),
        }
    }

    #[test]
    fn edge_bounds_are_inclusive() {
        assert!(at_origin_edge(500, "C", "T"));
        assert!(!at_origin_edge(501, "C", "T"));
        assert!(!at_origin_edge(15999, "C", "T"));
        assert!(at_origin_edge(16000, "C", "T"));
        // A deletion anchored before the edge reaches into it.
        assert!(at_origin_edge(15998, "CAT", "C"));
        // A padded insertion reads as its minimal form: AC>ACT at m.499 is
        // C>CT at m.500, an insertion after the edge.
        assert!(!at_origin_edge(500, "C", "CT"));
        assert!(!at_origin_edge(499, "AC", "ACT"));
    }

    #[test]
    fn low_level_m3243_is_flagged_and_keeps_attention() {
        let d = section(
            "chrM\t3243\t.\tA\tG\t50\tPASS\t.\tGT:AF\t0/1:0.07\n\
             chrM\t9000\t.\tC\tT\t50\tPASS\t.\tGT:AF\t1:0.99\n",
            "chrM\t3243\tA\tG\tPASS\t0.07\t9589\tPathogenic\tcriteria\tMELAS\t0\t0\t0\n\
             chrM\t9000\tC\tT\tPASS\t0.99\t7\tUncertain_significance\t.\t.\t.\t.\t.\n",
        );
        let (t, f) = (Value::Bool(true), Value::Bool(false));
        assert_eq!(at(&d, &["heteroplasmic"]), &Value::Int(1));
        assert_eq!(at(&d, &["heteroplasmic_low_af"]), &Value::Int(1));
        assert_eq!(at(&d, &["heteroplasmic_origin_edge"]), &Value::Int(0));
        let site = &items(at(&d, &["heteroplasmic_sites"]))[0];
        assert_eq!(
            (at(site, &["low_af"]), at(site, &["origin_edge"])),
            (&t, &f)
        );
        // What draws attention stays: the ClinVar P/LP count and the
        // confirmed variant, with its PASS filter and AF.
        assert_eq!(at(&d, &["annotation", "clinvar_plp"]), &Value::Int(1));
        let conf = items(at(&d, &["confirmed"]));
        assert_eq!(conf.len(), 1);
        assert_eq!(at(&conf[0], &["hgvs"]), &s("m.3243A>G"));
        assert_eq!(at(&conf[0], &["filter"]), &s("PASS"));
        assert_eq!(at(&conf[0], &["af"]), &s("0.07"));
        let listed = items(at(&d, &["annotation", "variants"]));
        assert_eq!(listed.len(), 2);
        assert_eq!(at(&listed[0], &["low_af"]), &t);
        // A homoplasmic call is never low level.
        assert_eq!(at(&listed[1], &["low_af"]), &f);
        assert_eq!(at(&listed[1], &["origin_edge"]), &f);
    }
}
