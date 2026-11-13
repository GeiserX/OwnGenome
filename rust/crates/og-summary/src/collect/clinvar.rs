//! The ClinVar screen (step 06): each P/LP hit with its review stars, its
//! call quality against step 06's floor, gnomAD AF and rsIDs from step 13's
//! VEP output, and the cDNA read check, the cDNA scan and the always-report
//! deep intronic sites when step 06 wrote them.
//!
//! The hits VCF is small and three sections read it (this one, the ACMG
//! tier of the clinical section and the carrier list), so it is read once
//! ([`HitsVcf`]) and each section parses the rows it needs.

use super::text::{self, ReadError, digits, format_map, gt_of, info_map, py_max};
use indexmap::IndexMap;
use og_clinvar::{ReviewStatus, combined_note, half_call, zygosity};
use pycompat::json::Value;
use std::collections::HashMap;

/// gnomAD AF above which a P/LP hit is marked common.
pub const COMMON_AF: f64 = 0.01;

/// Step 06's call floor: a hit whose GQ or DP is below it is low confidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Floor {
    pub gq: i128,
    pub dp: i128,
}

impl Default for Floor {
    fn default() -> Self {
        Floor { gq: 30, dp: 10 }
    }
}

impl Floor {
    fn to_json(self) -> Value {
        obj([("GQ", Value::Int(self.gq)), ("DP", Value::Int(self.dp))])
    }
}

/// `(GQ|DP)=(\d+)` over a header line: every match, the last one wins.
pub fn floor_values(line: &str, word_start: bool, floor: &mut Floor) {
    let b = line.as_bytes();
    for key in ["GQ=", "DP="] {
        let mut from = 0;
        while let Some(i) = line[from..].find(key).map(|i| i + from) {
            from = i + 1;
            if word_start && i > 0 && (b[i - 1].is_ascii_alphanumeric() || b[i - 1] == b'_') {
                continue;
            }
            let rest = &line[i + 3..];
            let n = rest.bytes().take_while(u8::is_ascii_digit).count();
            if n > 0
                && let Ok(v) = crate::int(&rest[..n])
            {
                // Python takes the matches in text order; with one key per
                // pass, the last of each key is the one that stays.
                if key == "GQ=" {
                    floor.gq = v;
                } else {
                    floor.dp = v;
                }
            }
        }
    }
}

/// A variant as `(chrom, pos, ref, alt)`.
pub type VarKey = (String, i128, String, String);
/// The hits VCF's call floor and data lines.
pub type HitRecords = (Floor, Vec<Vec<String>>);

/// The hits VCF, read once: its call floor and its data lines.
#[derive(Debug, Clone)]
pub struct HitsVcf {
    /// The path (full) and the path relative to the sample folder.
    pub path: String,
    pub rel: String,
    /// The floor and the records, or the note of a file that cannot be read.
    pub read: Result<HitRecords, String>,
}

impl HitsVcf {
    pub fn read(path: String, rel: String) -> HitsVcf {
        let read = text::lines(&path)
            .map(|lines| {
                let mut floor = Floor::default();
                for l in lines.iter().take_while(|l| l.starts_with("##")) {
                    if l.starts_with("##PGP_CALL_FLOOR=") {
                        floor_values(l, false, &mut floor);
                    }
                }
                let recs = lines
                    .iter()
                    .filter(|l| !l.starts_with('#'))
                    .map(|l| l.split('\t').map(str::to_owned).collect())
                    .collect();
                (floor, recs)
            })
            .map_err(|e: ReadError| e.unreadable(&path));
        HitsVcf { path, rel, read }
    }
}

/// One hit as `clinvar_row` reads it.
#[derive(Debug, Clone, PartialEq)]
pub struct Hit {
    pub stars: u8,
    pub chrom: String,
    pub pos: i128,
    pub reference: String,
    pub alt: String,
    pub genotype: String,
    pub gene: String,
    pub significance: String,
    pub review_status: String,
    pub combined_allele: String,
    pub dp: Option<i128>,
    pub gq: Option<i128>,
    pub low_confidence: bool,
    pub frequency: Option<Frequency>,
    pub cdna: Option<Value>,
}

/// gnomAD and rsIDs of one hit from VEP.
#[derive(Debug, Clone, PartialEq)]
pub struct Frequency {
    pub exome: Option<f64>,
    pub genome: Option<f64>,
    pub rsids: Vec<String>,
}

impl Hit {
    /// `clinvar_row`: an error is the collector's note (a POS that is no
    /// integer).
    pub fn parse(r: &[String], floor: Floor) -> Result<Hit, String> {
        let info = if r.len() > 7 {
            info_map(&r[7])
        } else {
            HashMap::new()
        };
        let get = |k: &str| info.get(k).copied().unwrap_or("");
        let rev = get("CLNREVSTAT");
        let genes: Vec<&str> = get("GENEINFO")
            .split('|')
            .filter(|g| !g.is_empty())
            .map(|g| g.split(':').next().unwrap_or(""))
            .collect();
        let gt = gt_of(r);
        let fmt = if r.len() > 9 {
            format_map(&r[8], &r[9])
        } else {
            HashMap::new()
        };
        let dp = digits(fmt.get("DP").copied());
        let gq = digits(fmt.get("GQ").copied());
        let pos = text::pos(r)?;
        let s = |v: &str| v.replace('_', " ");
        let or_dot = |v: String| if v.is_empty() { ".".to_owned() } else { v };
        Ok(Hit {
            stars: ReviewStatus::parse(rev.as_bytes()).stars(),
            chrom: r[0].clone(),
            pos,
            reference: r[3].clone(),
            alt: r[4].clone(),
            genotype: String::from_utf8_lossy(&zygosity(gt.as_bytes())).into_owned(),
            gene: or_dot(genes.join(",")),
            significance: or_dot(s(get("CLNSIG"))),
            review_status: or_dot(s(rev)),
            combined_allele: String::from_utf8_lossy(&combined_note(
                get("PGP_JOINED_FROM").as_bytes(),
                get("PGP_PART_OF").as_bytes(),
                get("PGP_ONE_ALLELE").as_bytes(),
                gt.as_bytes(),
            ))
            .into_owned(),
            dp,
            gq,
            low_confidence: gq.is_some_and(|v| v < floor.gq) || dp.is_some_and(|v| v < floor.dp),
            frequency: None,
            cdna: None,
        })
    }

    pub fn key(&self) -> (String, i128, String, String) {
        (
            self.chrom.clone(),
            self.pos,
            self.reference.clone(),
            self.alt.clone(),
        )
    }

    /// `chrom:pos ref>alt`.
    pub fn variant(&self) -> String {
        format!(
            "{}:{} {}>{}",
            self.chrom, self.pos, self.reference, self.alt
        )
    }

    fn to_json(&self) -> Value {
        let mut m = IndexMap::new();
        let put = |m: &mut IndexMap<String, Value>, k: &str, v: Value| {
            m.insert(k.to_owned(), v);
        };
        put(&mut m, "stars", Value::Int(self.stars.into()));
        put(&mut m, "chrom", Value::Str(self.chrom.clone()));
        put(&mut m, "pos", Value::Int(self.pos));
        put(&mut m, "ref", Value::Str(self.reference.clone()));
        put(&mut m, "alt", Value::Str(self.alt.clone()));
        put(&mut m, "genotype", Value::Str(self.genotype.clone()));
        put(&mut m, "gene", Value::Str(self.gene.clone()));
        put(
            &mut m,
            "significance",
            Value::Str(self.significance.clone()),
        );
        put(
            &mut m,
            "review_status",
            Value::Str(self.review_status.clone()),
        );
        put(
            &mut m,
            "combined_allele",
            Value::Str(self.combined_allele.clone()),
        );
        put(&mut m, "dp", opt_int(self.dp));
        put(&mut m, "gq", opt_int(self.gq));
        put(&mut m, "low_confidence", Value::Bool(self.low_confidence));
        if let Some(f) = &self.frequency {
            let af = py_max(f.exome.into_iter().chain(f.genome));
            put(&mut m, "gnomad_exome_af", opt_f(f.exome));
            put(&mut m, "gnomad_genome_af", opt_f(f.genome));
            put(
                &mut m,
                "rsids",
                Value::Array(f.rsids.iter().cloned().map(Value::Str).collect()),
            );
            put(&mut m, "af", opt_f(af));
            put(
                &mut m,
                "common",
                Value::Bool(af.is_some_and(|a| a > COMMON_AF)),
            );
        }
        if let Some(c) = &self.cdna {
            put(&mut m, "cdna", c.clone());
        }
        Value::Object(m)
    }
}

/// A hit's genotype as the carrier list words it.
pub fn carrier_genotype(gt: &str, hit: &Hit) -> String {
    if half_call(gt.as_bytes()) {
        hit.genotype.clone()
    } else {
        let (copies, ploidy) = allele_copies(gt);
        zygosity_word(copies, ploidy).to_owned()
    }
}

/// `(alt copies, ploidy)` of a GT: a missing allele is a copy not called,
/// not a missing copy, so the half of a split call is never hemizygous.
pub fn allele_copies(gt: &str) -> (usize, usize) {
    let a: Vec<&str> = gt.split(['/', '|']).filter(|x| !x.is_empty()).collect();
    (
        a.iter().filter(|x| !matches!(**x, "0" | ".")).count(),
        a.len(),
    )
}

pub fn zygosity_word(copies: usize, ploidy: usize) -> &'static str {
    if copies == 0 {
        "unknown"
    } else if ploidy == 1 {
        "hemizygous"
    } else if copies >= 2 {
        "hom"
    } else {
        "het"
    }
}

pub fn obj<'a>(pairs: impl IntoIterator<Item = (&'a str, Value)>) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

pub fn opt_int(x: Option<i128>) -> Value {
    x.map_or(Value::Null, Value::Int)
}

pub fn opt_f(x: Option<f64>) -> Value {
    x.map_or(Value::Null, Value::Float)
}

pub fn opt_s(x: Option<&str>) -> Value {
    x.map_or(Value::Null, |s| Value::Str(s.to_owned()))
}

// --- the cDNA read check, the scan, the deep intronic sites ----------------------

/// A row key of a cDNA check table: chrom and the alleles may be missing.
pub type CdnaKey = (Option<String>, i128, Option<String>, Option<String>);

/// `clinvar_cdna_checks`: the read check of each row, keyed by variant.
pub fn cdna_checks(path: &str) -> Result<HashMap<CdnaKey, Value>, String> {
    let rows = text::rows(path).map_err(|e| match e {
        ReadError::Os(_) => e.unreadable(path),
        ReadError::Other(n) => n,
    })?;
    let mut out = HashMap::new();
    for r in &rows {
        let Some(pos) = digits(text::cell(r, "pos")) else {
            continue;
        };
        let c = |k: &str| text::cell(r, k).map(str::to_owned);
        let or = |k: &str, d: &str| {
            Value::Str(
                text::cell(r, k)
                    .filter(|v| !v.is_empty())
                    .unwrap_or(d)
                    .to_owned(),
            )
        };
        let n = |k: &str| opt_int(digits(text::cell(r, k)));
        out.insert(
            (c("chrom"), pos, c("ref"), c("alt")),
            obj([
                ("check", or("cdna_check", "not_checked")),
                ("note", or("cdna_note", "")),
                ("alt_reads", n("alt_reads")),
                ("alt_clipped", n("alt_clipped")),
                ("ref_reads", n("ref_reads")),
                ("ref_clipped", n("ref_clipped")),
                (
                    "construct",
                    opt_s(text::cell(r, "construct_site").filter(|v| !v.is_empty())),
                ),
            ]),
        );
    }
    Ok(out)
}

/// The check of a variant the table does not hold.
pub fn not_in_table() -> Value {
    obj([
        ("check", Value::Str("not_checked".into())),
        ("note", Value::Str("not in the read-check table".into())),
        ("alt_reads", Value::Null),
        ("alt_clipped", Value::Null),
        ("ref_reads", Value::Null),
        ("ref_clipped", Value::Null),
    ])
}

pub fn cdna_key(chrom: &str, pos: i128, r: &str, a: &str) -> CdnaKey {
    (
        Some(chrom.to_owned()),
        pos,
        Some(r.to_owned()),
        Some(a.to_owned()),
    )
}

pub fn check_word(c: &Value) -> Option<&str> {
    match c {
        Value::Object(m) => match m.get("check") {
            Some(Value::Str(s)) => Some(s),
            _ => None,
        },
        _ => None,
    }
}

fn has_construct(c: &Value) -> bool {
    matches!(c, Value::Object(m) if m.get("construct").is_some_and(pycompat::value::truthy))
}

/// `clinvar_cdna_scan`: the scan's note and loci.
fn cdna_scan(path: &str) -> Result<Value, String> {
    let lines = text::lines(path).map_err(|e| e.unreadable(path))?;
    let note = lines
        .first()
        .and_then(|l| l.strip_prefix("# cdna_scan: "))
        .unwrap_or("")
        .to_owned();
    let body = lines
        .iter()
        .filter(|l| !l.starts_with('#'))
        .map(|l| Ok::<_, std::io::Error>(l.clone()));
    let mut loci = Vec::new();
    for r in pycompat::csv::DictReader::new(body, '\t', pycompat::csv::Quoting::Minimal) {
        let r = r.map_err(|e| format!("{}: {e}", e.exception()))?;
        let Some(start) = digits(text::cell(&r, "start")) else {
            continue;
        };
        let int0 = |k: &str| -> Result<Value, String> {
            match text::cell(&r, k).filter(|v| !v.is_empty()) {
                None => Ok(Value::Int(0)),
                Some(v) => crate::int(v)
                    .map(Value::Int)
                    .map_err(|e| format!("ValueError: {e}")),
            }
        };
        loci.push(obj([
            ("chrom", opt_s(text::cell(&r, "chrom"))),
            ("start", Value::Int(start)),
            ("end", int0("end")?),
            ("junctions", int0("junctions")?),
            ("reads", int0("reads")?),
            (
                "junction_list",
                Value::Str(text::cell_or(&r, "junction_list").to_owned()),
            ),
        ]));
    }
    Ok(obj([
        ("note", Value::Str(note)),
        ("loci", Value::Array(loci)),
    ]))
}

/// `clinvar_deep_intronic`: one row per always-report site.
fn deep_intronic(path: &str) -> Result<Value, String> {
    let rows = text::rows(path).map_err(|e| match e {
        ReadError::Os(_) => e.unreadable(path),
        ReadError::Other(n) => n,
    })?;
    let mut out = Vec::new();
    for r in &rows {
        let Some(pos) = digits(text::cell(r, "pos")) else {
            continue;
        };
        let or = |k: &str, d: &str| {
            Value::Str(
                text::cell(r, k)
                    .filter(|v| !v.is_empty())
                    .unwrap_or(d)
                    .to_owned(),
            )
        };
        out.push(obj([
            ("chrom", opt_s(text::cell(r, "chrom"))),
            ("pos", Value::Int(pos)),
            ("ref", opt_s(text::cell(r, "ref"))),
            ("alt", opt_s(text::cell(r, "alt"))),
            ("gene", or("gene", ".")),
            ("hgvs", or("hgvs", "")),
            ("clinvar_id", or("clinvar_id", "")),
            ("call", or("call", "no call")),
            ("detail", or("detail", "")),
        ]));
    }
    Ok(Value::Array(out))
}

// --- gnomAD AF and rsIDs from step 13's VEP output -----------------------------------

/// Where the frequencies came from, or why there are none.
#[derive(Debug, Clone, PartialEq)]
pub struct FrequencySource {
    pub source: Option<String>,
    pub note: Option<String>,
}

/// The VEP output's frequencies of every hit, in one pass per contig over
/// the hits sorted by position: one indexed query per hit position, served
/// by the reader's block cache when hits share a block.
pub fn vep_frequencies(vep: &str, keys: &[VarKey]) -> Result<HashMap<VarKey, Frequency>, String> {
    let mut fmt: Option<Vec<String>> = None;
    {
        let mut lines = og_io::Lines::new(
            og_io::open(vep).map_err(|e| crate::os_error(&e, vep))?,
            og_io::Utf8::Replace,
        );
        let mut line = String::new();
        while lines.read_line(&mut line).map_err(|e| e.to_string())? {
            if !line.starts_with('#') {
                break;
            }
            if line.starts_with("##INFO=<ID=CSQ") {
                fmt = line.split_once("Format: ").map(|(_, f)| {
                    f.split('"')
                        .next()
                        .unwrap_or("")
                        .split('|')
                        .map(str::to_owned)
                        .collect()
                });
                // `([^"]+)` needs one character: an empty format is none.
                if fmt
                    .as_ref()
                    .is_some_and(|f| f.len() == 1 && f[0].is_empty())
                {
                    fmt = None;
                }
            }
        }
    }
    let Some(fmt) = fmt.filter(|f| f.iter().any(|c| c == "Allele")) else {
        return Ok(HashMap::new());
    };
    let index = |c: &str| fmt.iter().position(|x| x == c);
    let cols = [index("gnomADe_AF"), index("gnomADg_AF")];
    let ev = index("Existing_variation");
    let allele_ix = index("Allele").unwrap_or(0);
    let mut sorted: Vec<&VarKey> = keys.iter().collect();
    sorted.sort();
    sorted.dedup();
    let mut out = HashMap::new();
    if sorted.is_empty() {
        return Ok(out);
    }
    let mut reader =
        og_io::tabix::Reader::with_index(vep, format!("{vep}.tbi")).map_err(|e| e.to_string())?;
    for key in sorted {
        let (chrom, pos, reference, alt) = key;
        let Ok(p) = u64::try_from(*pos) else { continue };
        if p == 0 {
            continue;
        }
        let lines: Vec<Vec<u8>> = reader
            .query(chrom, p - 1, p)
            .collect::<Result<_, _>>()
            .map_err(|e| e.to_string())?;
        for raw in lines {
            let line = String::from_utf8_lossy(&raw);
            let r: Vec<&str> = line.split('\t').collect();
            if r.len() < 3 || r[1].parse::<i128>().ok() != Some(*pos) {
                continue;
            }
            let alts: Vec<&str> = if r.len() > 7 {
                r[4].split(',').collect()
            } else {
                Vec::new()
            };
            if r.get(3) != Some(&reference.as_str()) || !alts.contains(&alt.as_str()) {
                continue;
            }
            // VEP's Allele drops the base an indel's alleles share.
            let first = |s: &str| s.chars().next();
            let snv = reference.chars().count() == 1 && alts.iter().all(|a| a.chars().count() == 1);
            let allele = if snv || alts.iter().any(|a| first(a) != first(reference)) {
                alt.clone()
            } else {
                let rest: String = alt.chars().skip(1).collect();
                if rest.is_empty() {
                    "-".to_owned()
                } else {
                    rest
                }
            };
            let mut got = [None, None];
            let mut ids: Vec<String> = Vec::new();
            for tr in info_map(r[7]).get("CSQ").copied().unwrap_or("").split(',') {
                let v: Vec<&str> = tr.split('|').collect();
                if v.len() != fmt.len() || (alts.len() > 1 && v[allele_ix] != allele) {
                    continue;
                }
                for (i, c) in cols.iter().enumerate() {
                    if let Some(c) = c
                        && got[i].is_none()
                    {
                        got[i] = py_max(
                            v[*c]
                                .split('&')
                                .filter_map(|x| pycompat::num::parse_float(x).ok()),
                        );
                    }
                }
                if let Some(e) = ev {
                    for x in v[e].split('&') {
                        if is_rsid(x) {
                            let x = x.to_lowercase();
                            if !ids.contains(&x) {
                                ids.push(x);
                            }
                        }
                    }
                }
            }
            ids.sort_by_key(|x| x[2..].parse::<u128>().unwrap_or(u128::MAX));
            out.insert(
                key.clone(),
                Frequency {
                    exome: got[0],
                    genome: got[1],
                    rsids: ids,
                },
            );
        }
    }
    Ok(out)
}

/// `(?i)rs\d+`, the whole field.
fn is_rsid(x: &str) -> bool {
    let b = x.as_bytes();
    b.len() > 2 && b[..2].eq_ignore_ascii_case(b"rs") && b[2..].iter().all(u8::is_ascii_digit)
}

// --- the section --------------------------------------------------------------------

/// The data of the ClinVar section.
#[derive(Debug, Clone, PartialEq)]
pub struct ClinVar {
    pub hits: Vec<Hit>,
    pub floor: Floor,
    pub frequency: FrequencySource,
    pub cdna: bool,
    pub cdna_scan: Option<Value>,
    pub deep_intronic: Option<Value>,
}

impl ClinVar {
    /// `sec_clinvar` over the hits read once. `side(rel)` is the file of
    /// another of step 06's or 13's outputs, found in the folder.
    pub fn read(
        hits: &HitsVcf,
        dir: &str,
        side: &dyn Fn(&str) -> Option<(String, String)>,
        s: &str,
    ) -> Result<ClinVar, String> {
        let (floor, recs) = hits.read.clone()?;
        let mut rows = Vec::new();
        for r in recs.iter().filter(|r| r.len() >= 8) {
            rows.push(Hit::parse(r, floor)?);
        }
        rows.sort_by(|a, b| {
            b.stars
                .cmp(&a.stars)
                .then_with(|| a.chrom.cmp(&b.chrom))
                .then_with(|| a.pos.cmp(&b.pos))
        });
        let frequency = frequencies(dir, side, s, &mut rows);
        let mut cdna = false;
        if let Some((p, _)) = side(&format!("clinvar/{s}_clinvar_cdna_check.tsv")) {
            let table = cdna_checks(&p)?;
            for x in &mut rows {
                x.cdna = Some(
                    table
                        .get(&cdna_key(&x.chrom, x.pos, &x.reference, &x.alt))
                        .cloned()
                        .unwrap_or_else(not_in_table),
                );
            }
            cdna = true;
        }
        let cdna_scan = match side(&format!("clinvar/{s}_clinvar_cdna_scan.tsv")) {
            Some((p, _)) => Some(cdna_scan(&p)?),
            None => None,
        };
        let deep_intronic = match side(&format!("clinvar/{s}_clinvar_deep_intronic.tsv")) {
            Some((p, _)) => Some(deep_intronic(&p)?),
            None => None,
        };
        Ok(ClinVar {
            hits: rows,
            floor,
            frequency,
            cdna,
            cdna_scan,
            deep_intronic,
        })
    }

    pub fn to_json(&self) -> Value {
        let sure: Vec<&Hit> = self.hits.iter().filter(|x| !x.low_confidence).collect();
        let by = (0..=4u8)
            .rev()
            .map(|k| {
                (
                    k.to_string(),
                    Value::Int(sure.iter().filter(|x| x.stars == k).count() as i128),
                )
            })
            .collect();
        let mut m: IndexMap<String, Value> = IndexMap::new();
        m.insert("count".into(), Value::Int(sure.len() as i128));
        m.insert("by_stars".into(), Value::Object(by));
        m.insert(
            "hits".into(),
            Value::Array(self.hits.iter().map(Hit::to_json).collect()),
        );
        m.insert(
            "low_confidence".into(),
            Value::Int((self.hits.len() - sure.len()) as i128),
        );
        m.insert("call_floor".into(), self.floor.to_json());
        m.insert(
            "frequency".into(),
            obj([
                ("source", opt_s(self.frequency.source.as_deref())),
                ("note", opt_s(self.frequency.note.as_deref())),
                ("common_af", Value::Float(COMMON_AF)),
            ]),
        );
        if self.cdna {
            let checks: Vec<&Value> = self.hits.iter().filter_map(|x| x.cdna.as_ref()).collect();
            m.insert(
                "cdna_flagged".into(),
                Value::Int(
                    checks
                        .iter()
                        .filter(|c| check_word(c) == Some("flagged"))
                        .count() as i128,
                ),
            );
            m.insert(
                "construct_sites".into(),
                Value::Int(checks.iter().filter(|c| has_construct(c)).count() as i128),
            );
        }
        if let Some(v) = &self.cdna_scan {
            m.insert("cdna_scan".into(), v.clone());
        }
        if let Some(v) = &self.deep_intronic {
            m.insert("deep_intronic".into(), v.clone());
        }
        Value::Object(m)
    }
}

/// `clinvar_frequencies`: puts the frequencies on each hit and says where
/// they came from.
fn frequencies(
    dir: &str,
    side: &dyn Fn(&str) -> Option<(String, String)>,
    s: &str,
    rows: &mut [Hit],
) -> FrequencySource {
    let vep = side(&format!("vep/{s}_vep.vcf.gz"));
    let _ = dir;
    let mut note = None;
    let mut found = HashMap::new();
    match &vep {
        None => note = Some("no VEP output (step 13) in this folder".to_owned()),
        Some((p, rel)) if !std::path::Path::new(&format!("{p}.tbi")).exists() => {
            note = Some(format!("{rel} has no .tbi index"));
        }
        Some((p, rel)) => {
            let keys: Vec<_> = rows.iter().map(Hit::key).collect();
            match vep_frequencies(p, &keys) {
                Ok(f) => found = f,
                Err(e) => note = Some(format!("{rel} could not be read: {e}")),
            }
        }
    }
    for x in rows.iter_mut() {
        x.frequency = Some(found.get(&x.key()).cloned().unwrap_or(Frequency {
            exome: None,
            genome: None,
            rsids: Vec::new(),
        }));
    }
    FrequencySource {
        source: if note.is_none() {
            vep.map(|(_, rel)| rel)
        } else {
            None
        },
        note,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floor_takes_the_last_value_of_each_key() {
        let mut f = Floor::default();
        floor_values("##PGP_CALL_FLOOR=GQ=20,DP=8,GQ=25", false, &mut f);
        assert_eq!(f, Floor { gq: 25, dp: 8 });
        let mut f = Floor::default();
        floor_values("##PGP_CALL_FLOOR=XGQ=5,DP=7", true, &mut f);
        assert_eq!(f, Floor { gq: 30, dp: 7 });
    }

    #[test]
    fn copies_and_words() {
        assert_eq!(allele_copies("0/1"), (1, 2));
        assert_eq!(allele_copies("1|1"), (2, 2));
        assert_eq!(allele_copies("1"), (1, 1));
        assert_eq!(allele_copies("1/."), (1, 2));
        assert_eq!(allele_copies(""), (0, 0));
        assert_eq!(zygosity_word(1, 1), "hemizygous");
        assert_eq!(zygosity_word(0, 2), "unknown");
    }

    #[test]
    fn rsid_fields() {
        assert!(is_rsid("rs123"));
        assert!(is_rsid("RS9"));
        assert!(!is_rsid("rs"));
        assert!(!is_rsid("COSV1"));
    }
}
