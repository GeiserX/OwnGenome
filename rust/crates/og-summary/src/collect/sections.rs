//! The readers of group a: depth, gene callability, sex and karyotype,
//! sample identity, mosaic sex-chromosome loss, read count and variant
//! counts. Each returns the file it read and its typed data, `None` when
//! the step wrote nothing, or the note of an unreadable section.

use super::model::opt_str;
use crate::{basename, float, int, text_lines};
use indexmap::IndexMap;
use pycompat::csv::{DictReader, DictRow, Quoting};
use pycompat::json::Value;
use std::io::BufRead;

/// The data of one section.
#[derive(Debug, Clone, PartialEq)]
pub enum Data {
    Coverage(Coverage),
    Callability(Callability),
    SexCheck(SexCheck),
    SampleQc(SampleQc),
    MosaicSexLoss(Mosaic),
    ReadCount(ReadCount),
    Variants(Variants),
}

impl Data {
    pub fn to_json(&self) -> Value {
        match self {
            Data::Coverage(d) => d.to_json(),
            Data::Callability(d) => d.to_json(),
            Data::SexCheck(d) => d.to_json(),
            Data::SampleQc(d) => d.to_json(),
            Data::MosaicSexLoss(d) => d.to_json(),
            Data::ReadCount(d) => d.to_json(),
            Data::Variants(d) => d.to_json(),
        }
    }
}

fn obj(pairs: impl IntoIterator<Item = (&'static str, Value)>) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

fn opt_float(x: Option<f64>) -> Value {
    x.map_or(Value::Null, Value::Float)
}

/// `Exception` text as the Python collector's catch-all wrote it.
pub fn py_error(kind: &str, msg: &str) -> String {
    format!("{kind}: {msg}")
}

/// `float(x)` that is a number and not NaN, else `None` (`num` in Python).
pub fn num(x: Option<&str>) -> Option<f64> {
    pycompat::num::parse_float(x?).ok().filter(|v| !v.is_nan())
}

// --- tables ------------------------------------------------------------------------

/// The rows of a tab-separated table with a header, as `csv.DictReader`
/// reads it. An error is the note of the section.
pub fn read_tsv(path: &str) -> Result<Vec<DictRow>, String> {
    let mut lines = text_lines(path).map_err(|e| os_note(&e))?;
    let iter = std::iter::from_fn(move || {
        let mut s = String::new();
        match lines.read_line(&mut s) {
            Ok(true) => Some(Ok(s)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    });
    DictReader::new(iter, '\t', Quoting::Minimal)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| py_error(e.exception(), &e.to_string()))
}

/// An `OSError` from opening a file, as `type(e).__name__: e`.
fn os_note(e: &str) -> String {
    let kind = if e.contains("No such file") {
        "FileNotFoundError"
    } else if e.contains("Permission denied") {
        "PermissionError"
    } else if e.contains("Is a directory") {
        "IsADirectoryError"
    } else {
        "OSError"
    };
    py_error(kind, e)
}

/// The key and value rows of a step's `key\tvalue` table: a row with an
/// empty key is left out, a later key replaces an earlier one in its place.
/// The value is `None` when the row has no value cell, `""` when the table
/// has no value column.
fn key_values(rows: &[DictRow]) -> IndexMap<String, Option<String>> {
    let mut out = IndexMap::new();
    for r in rows {
        if let Some(Some(k)) = r.get("key")
            && !k.is_empty()
        {
            let v = match r.get("value") {
                None => Some(String::new()),
                Some(v) => v.map(str::to_owned),
            };
            out.insert(k.to_owned(), v);
        }
    }
    out
}

fn kv_json(rows: &IndexMap<String, Option<String>>) -> IndexMap<String, Value> {
    rows.iter()
        .map(|(k, v)| (k.clone(), opt_str(v.as_deref())))
        .collect()
}

// --- coverage (step 16b, mosdepth) -----------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct Coverage {
    /// The `mean` of the `total` row, the last one when there are several.
    pub mean_depth: Option<f64>,
}

impl Coverage {
    pub fn read(path: &str) -> Result<Coverage, String> {
        let mut mean = None;
        for r in read_tsv(path)? {
            if r.get("chrom") == Some(Some("total")) {
                mean = Some(r.get("mean").flatten().map(str::to_owned));
            }
        }
        let mean_depth = match mean.flatten() {
            Some(m) if !m.is_empty() => Some(float(&m).map_err(|e| py_error("ValueError", &e))?),
            _ => None,
        };
        Ok(Coverage { mean_depth })
    }

    fn to_json(&self) -> Value {
        obj([("mean_depth", opt_float(self.mean_depth))])
    }
}

// --- gene callability (step 16b) -----------------------------------------------------

/// A gene counts as read when this share of its bases reached the depth.
pub const CALLABLE_SHARE: f64 = 0.9;

/// The share of a gene's bases at the depth, and what it means.
#[derive(Debug, Clone, PartialEq)]
pub enum Reading {
    Read(f64),
    NotRead(f64),
    NotOnReference,
}

impl Reading {
    fn word(&self) -> &'static str {
        match self {
            Reading::Read(_) => "read",
            Reading::NotRead(_) => "not read",
            Reading::NotOnReference => "not on this reference",
        }
    }

    fn share(&self) -> Value {
        match self {
            Reading::Read(s) | Reading::NotRead(s) => Value::Float(*s),
            Reading::NotOnReference => Value::Null,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Gene {
    pub gene: String,
    pub hgnc_id: String,
    pub reasons: Vec<String>,
    pub owner: Option<String>,
    pub reading: Reading,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Callability {
    pub min_depth: i128,
    pub genes: Vec<Gene>,
}

/// One row of `assets/callability_genes.tsv`.
pub struct CatalogRow {
    pub gene: String,
    pub hgnc_id: String,
    pub chrom: String,
    pub start: String,
    pub reasons: String,
    pub owner: String,
}

/// Reads the catalog: its comment lines left out, cells as written.
pub fn read_catalog(assets: &str) -> Result<Vec<CatalogRow>, String> {
    const NAME: &str = "callability_genes.tsv";
    let p = crate::join(assets, NAME);
    if !crate::is_file(&p) {
        return Err(format!(
            "{NAME} is not in {assets}: the report needs the repository's assets/ folder (render_report.py --assets)"
        ));
    }
    let text = std::fs::read(&p).map_err(|e| py_error("OSError", &e.to_string()))?;
    let text =
        String::from_utf8(text).map_err(|e| py_error("UnicodeDecodeError", &e.to_string()))?;
    let lines = text
        .split_inclusive('\n')
        .filter(|l| !l.starts_with('#'))
        .map(|l| Ok::<_, std::io::Error>(l.to_owned()));
    let mut out = Vec::new();
    for r in DictReader::new(lines, '\t', Quoting::None) {
        let r = r.map_err(|e| py_error(e.exception(), &e.to_string()))?;
        let cell = |k: &str| -> Result<String, String> {
            match r.get(k) {
                None => Err(py_error("KeyError", &crate::repr(k))),
                Some(v) => Ok(v.unwrap_or_default().to_owned()),
            }
        };
        out.push(CatalogRow {
            gene: cell("gene")?,
            hgnc_id: cell("hgnc_id")?,
            chrom: cell("chrom")?,
            start: cell("start")?,
            reasons: cell("reasons")?,
            owner: cell("owner")?,
        });
    }
    Ok(out)
}

impl Callability {
    /// Step 16b's per-gene table (mosdepth `--thresholds` of one depth over
    /// the catalog's genes: chrom, start, end, gene, bases at the depth)
    /// read against the catalog.
    pub fn read(catalog: &[CatalogRow], table: &str) -> Result<Callability, String> {
        let name = basename(table);
        let mut rows: Vec<Vec<String>> = Vec::new();
        let mut lines = text_lines(table).map_err(|e| os_note(&e))?;
        let mut line = String::new();
        while lines
            .read_line(&mut line)
            .map_err(|e| py_error("OSError", &e.to_string()))?
        {
            if pycompat::num::py_strip(&line).is_empty() {
                continue;
            }
            let body = line.trim_end_matches('\n');
            rows.push(body.split('\t').map(str::to_owned).collect());
        }
        let depth_col = rows.first().filter(|h| {
            h[0].starts_with('#')
                && h.len() == 5
                && h[4]
                    .strip_suffix('X')
                    .is_some_and(|d| !d.is_empty() && d.bytes().all(|c| c.is_ascii_digit()))
        });
        let Some(head) = depth_col else {
            return Err(format!(
                "{name} is not mosdepth's thresholds table of one depth"
            ));
        };
        let min_depth =
            int(&head[4][..head[4].len() - 1]).map_err(|e| py_error("ValueError", &e))?;
        let mut bases: std::collections::HashMap<&str, (i128, i128)> =
            std::collections::HashMap::new();
        for c in &rows[1..] {
            if c.len() != 5 {
                return Err(format!("{name}: a row has {} columns, not 5", c.len()));
            }
            let v = |s: &str| int(s).map_err(|e| py_error("ValueError", &e));
            let (start, end, ok) = (v(&c[1])?, v(&c[2])?, v(&c[4])?);
            bases.insert(&c[3], (end - start, ok));
        }
        if bases.is_empty() {
            return Err(format!("{name} has no gene rows"));
        }
        // mosdepth writes no row on a contig the BAM lacks: a chromosome with
        // no row at all is not on the sample's reference; a gene missing from
        // a chromosome that has rows means another catalog wrote the table.
        let contigs: std::collections::HashSet<&str> =
            rows[1..].iter().map(|c| c[0].as_str()).collect();
        let mut genes = Vec::with_capacity(catalog.len());
        for r in catalog {
            let reading = if r.start == "." || !contigs.contains(r.chrom.as_str()) {
                Reading::NotOnReference
            } else {
                let Some(&(n, ok)) = bases.get(r.gene.as_str()) else {
                    return Err(format!(
                        "{name} has no row for {}: it was written from another assets/callability_genes.tsv; run step 16b again",
                        r.gene
                    ));
                };
                let share = if n == 0 {
                    0.0
                } else {
                    pycompat::fmt::round(ok as f64 / n as f64, 4)
                        .map_err(|e| py_error(e.exception(), &e.to_string()))?
                };
                if share >= CALLABLE_SHARE {
                    Reading::Read(share)
                } else {
                    Reading::NotRead(share)
                }
            };
            genes.push(Gene {
                gene: r.gene.clone(),
                hgnc_id: r.hgnc_id.clone(),
                reasons: r.reasons.split(',').map(str::to_owned).collect(),
                owner: (r.owner != ".").then(|| r.owner.clone()),
                reading,
            });
        }
        Ok(Callability { min_depth, genes })
    }

    fn to_json(&self) -> Value {
        let not_read = |owned: bool| {
            Value::Array(
                self.genes
                    .iter()
                    .filter(|g| {
                        !matches!(g.reading, Reading::Read(_)) && g.owner.is_some() == owned
                    })
                    .map(|g| {
                        let who = if owned {
                            ("owner", opt_str(g.owner.as_deref()))
                        } else {
                            ("hgnc_id", Value::Str(g.hgnc_id.clone()))
                        };
                        let mut m = vec![("gene", Value::Str(g.gene.clone())), who];
                        if !owned {
                            m.push((
                                "reasons",
                                Value::Array(
                                    g.reasons.iter().map(|r| Value::Str(r.clone())).collect(),
                                ),
                            ));
                        }
                        m.push(("share", g.reading.share()));
                        m.push(("status", Value::Str(g.reading.word().into())));
                        obj(m)
                    })
                    .collect(),
            )
        };
        let read = self
            .genes
            .iter()
            .filter(|g| matches!(g.reading, Reading::Read(_)))
            .count();
        obj([
            ("min_depth", Value::Int(self.min_depth)),
            ("min_mapq", Value::Int(20)),
            ("min_share", Value::Float(CALLABLE_SHARE)),
            ("catalog_genes", Value::Int(self.genes.len() as i128)),
            ("read", Value::Int(read as i128)),
            ("not_read", not_read(false)),
            ("not_read_owned", not_read(true)),
            (
                "not_on_reference",
                Value::Array(
                    self.genes
                        .iter()
                        .filter(|g| g.reading == Reading::NotOnReference)
                        .map(|g| Value::Str(g.gene.clone()))
                        .collect(),
                ),
            ),
        ])
    }
}

// --- karyotype ---------------------------------------------------------------------

/// How far a copy number may sit from a whole number and still be read as it.
pub const KARYOTYPE_TOL: f64 = 0.3;

/// The sex-chromosome aneuploidies named, by (chrX, chrY) copies. Any other
/// pair of whole numbers names nothing: in a run it is far more often depth
/// noise or a broken input than a karyotype.
const KARYOTYPES: [((i64, i64), &str); 9] = [
    ((1, 0), "45,X"),
    ((2, 1), "47,XXY"),
    ((1, 2), "47,XYY"),
    ((3, 0), "47,XXX"),
    ((2, 2), "48,XXYY"),
    ((3, 1), "48,XXXY"),
    ((4, 0), "48,XXXX"),
    ((4, 1), "49,XXXXY"),
    ((5, 0), "49,XXXXX"),
];

/// The karyotype two copy numbers name; `None` for 46,XX and 46,XY, for a
/// missing value, a value between whole numbers, and any other pair.
pub fn karyotype(x: Option<f64>, y: Option<f64>) -> Option<&'static str> {
    let (x, y) = (x?, y?);
    if !(-1.0 < x && x < 10.0 && -1.0 < y && y < 10.0) {
        return None;
    }
    let (rx, ry) = (x.round_ties_even(), y.round_ties_even());
    if (x - rx).abs() > KARYOTYPE_TOL || (y - ry).abs() > KARYOTYPE_TOL {
        return None;
    }
    KARYOTYPES
        .iter()
        .find(|(k, _)| *k == (rx as i64, ry as i64))
        .map(|(_, n)| *n)
}

// --- sex check (step 16, indexcov) ----------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
pub struct SexCheck {
    /// `male`, `female` or `unknown`, from the ped's sex code.
    pub inferred_sex: &'static str,
    pub cn_chrx: String,
    pub cn_chry: String,
    pub karyotype: Option<&'static str>,
}

impl SexCheck {
    /// The last sample row of indexcov's ped file.
    pub fn read(path: &str) -> Result<SexCheck, String> {
        let text = std::fs::read(path).map_err(|e| os_note(&crate::os_error(&e, path)))?;
        // Plain text whatever the name, as Python's open() reads it.
        let mut lines = og_io::Lines::new(&text[..], og_io::Utf8::Replace);
        let mut rows: Vec<Vec<String>> = Vec::new();
        let mut l = String::new();
        while lines
            .read_line(&mut l)
            .map_err(|e| py_error("OSError", &e.to_string()))?
        {
            if !pycompat::num::py_strip(&l).is_empty() {
                rows.push(
                    l.trim_end_matches('\n')
                        .split('\t')
                        .map(str::to_owned)
                        .collect(),
                );
            }
        }
        let head: Vec<&str> = rows.first().map_or(Vec::new(), |h| {
            h.iter().map(|c| c.trim_start_matches('#')).collect()
        });
        if !head.contains(&"sex") || rows.len() < 2 {
            return Err(format!(
                "{} has no sex column or no sample row",
                basename(path)
            ));
        }
        let last = rows.last().expect("two rows");
        let mut row: IndexMap<&str, &str> = IndexMap::new();
        for (k, v) in head.iter().zip(last) {
            row.insert(k, v);
        }
        let cell = |k: &str| row.get(k).copied();
        Ok(SexCheck {
            inferred_sex: match cell("sex") {
                Some("1") => "male",
                Some("2") => "female",
                _ => "unknown",
            },
            cn_chrx: cell("CNchrX").unwrap_or("").to_owned(),
            cn_chry: cell("CNchrY").unwrap_or("").to_owned(),
            karyotype: karyotype(num(cell("CNchrX")), num(cell("CNchrY"))),
        })
    }

    fn to_json(&self) -> Value {
        obj([
            ("inferred_sex", Value::Str(self.inferred_sex.into())),
            ("cn_chrX", Value::Str(self.cn_chrx.clone())),
            ("cn_chrY", Value::Str(self.cn_chry.clone())),
            ("karyotype", opt_str(self.karyotype)),
        ])
    }
}

// --- sample identity and contamination (step 33) ----------------------------------------

/// Step 33's verdict table, row by row, and the karyotype somalier's depth
/// ratios name.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleQc {
    pub rows: IndexMap<String, Option<String>>,
    pub karyotype: Option<&'static str>,
}

impl SampleQc {
    pub fn read(path: &str, rows: &[DictRow]) -> Result<SampleQc, String> {
        let rows = key_values(rows);
        if !rows.contains_key("inferred_sex") {
            return Err(format!("{} has no inferred_sex row", basename(path)));
        }
        let get = |k: &str| rows.get(k).and_then(|v| v.as_deref());
        // somalier's chrX and chrY depth over the autosomal depth: one copy
        // reads 0.5. Named only from more than 10 chrX sites, the floor
        // somalier's own sex call needs; a thinner ratio is noise.
        let (xr, yr, xs) = (
            num(get("x_depth_ratio")),
            num(get("y_depth_ratio")),
            num(get("x_sites")),
        );
        let karyotype = match (xr, yr, xs) {
            (Some(xr), Some(yr), Some(xs)) if xs.is_finite() && xs > 10.0 => {
                karyotype(Some(2.0 * xr), Some(2.0 * yr))
            }
            _ => None,
        };
        Ok(SampleQc { rows, karyotype })
    }

    /// The table's `declared_sex` row.
    pub fn declared_sex(&self) -> Option<&str> {
        self.rows.get("declared_sex").and_then(|v| v.as_deref())
    }

    fn to_json(&self) -> Value {
        let mut m = kv_json(&self.rows);
        m.insert("karyotype".into(), opt_str(self.karyotype));
        Value::Object(m)
    }
}

// --- mosaic loss of chrY or of one chrX (step 33, og qc sex-loss) -------------------------

pub const MOSAIC_CALLS: [&str; 4] = ["detected", "not_detected", "gain", "not_estimated"];
pub const MOSAIC_NUMBERS: [&str; 9] = [
    "depth_ratio",
    "depth_ratio_expected",
    "loss_depth",
    "loss_depth_low",
    "loss_depth_high",
    "loss_af",
    "loss_af_low",
    "loss_af_high",
    "floor",
];

/// The `mosaic_*` rows of the step 33 table, without the prefix.
#[derive(Debug, Clone, PartialEq)]
pub struct Mosaic {
    /// Every row as written, in table order.
    pub rows: IndexMap<String, String>,
    /// The numbers, finite or `None`, in [`MOSAIC_NUMBERS`] order.
    pub numbers: [Option<f64>; 9],
    pub x_het_sites: Option<i128>,
}

impl Mosaic {
    /// `None` for a table written before the estimate existed.
    pub fn read(path: &str, rows: &[DictRow]) -> Result<Option<Mosaic>, String> {
        let mut m = IndexMap::new();
        for r in rows {
            if let Some(Some(k)) = r.get("key")
                && let Some(rest) = k.strip_prefix("mosaic_")
            {
                m.insert(
                    rest.to_owned(),
                    r.get("value").flatten().unwrap_or("").to_owned(),
                );
            }
        }
        if m.is_empty() {
            return Ok(None);
        }
        if !m
            .get("call")
            .is_some_and(|c| MOSAIC_CALLS.contains(&c.as_str()))
        {
            return Err(format!(
                "{} has no mosaic_call row ({})",
                basename(path),
                MOSAIC_CALLS.join(", ")
            ));
        }
        let numbers =
            MOSAIC_NUMBERS.map(|k| num(m.get(k).map(String::as_str)).filter(|v| v.is_finite()));
        let x_het_sites = match m.get("x_het_sites") {
            Some(v) if pycompat::text::isdigit(v) => {
                Some(int(v).map_err(|e| py_error("ValueError", &e))?)
            }
            _ => None,
        };
        Ok(Some(Mosaic {
            rows: m,
            numbers,
            x_het_sites,
        }))
    }

    fn to_json(&self) -> Value {
        let mut m: IndexMap<String, Value> = self
            .rows
            .iter()
            .map(|(k, v)| (k.clone(), Value::Str(v.clone())))
            .collect();
        for (k, v) in MOSAIC_NUMBERS.iter().zip(self.numbers) {
            m.insert((*k).to_owned(), opt_float(v));
        }
        m.insert(
            "x_het_sites".into(),
            self.x_het_sites.map_or(Value::Null, Value::Int),
        );
        Value::Object(m)
    }
}

// --- read count (step 02) ----------------------------------------------------------------

pub const READ_COUNTS: [&str; 3] = ["input_reads", "primary_records", "missing_reads"];

/// Step 02's read count: the reads that went into the aligner and the
/// primary records of the BAM it wrote.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadCount {
    pub rows: IndexMap<String, Option<String>>,
    /// [`READ_COUNTS`] in order.
    pub counts: [i128; 3],
}

impl ReadCount {
    pub fn read(path: &str) -> Result<ReadCount, String> {
        let rows = key_values(&read_tsv(path)?);
        let mut counts = [0; 3];
        for (i, k) in READ_COUNTS.iter().enumerate() {
            let v = rows.get(*k).and_then(|v| v.as_deref()).unwrap_or("");
            let digits = v.strip_prefix('-').unwrap_or(v);
            if digits.is_empty() || !digits.bytes().all(|c| c.is_ascii_digit()) {
                return Err(format!("{} has no whole number for {k}", basename(path)));
            }
            counts[i] = int(v).map_err(|e| py_error("ValueError", &e))?;
        }
        let check = rows.get("read_count_check").and_then(|v| v.as_deref());
        if !matches!(check, Some("ok" | "mismatch")) {
            return Err(format!(
                "{} has no read_count_check row (ok or mismatch)",
                basename(path)
            ));
        }
        Ok(ReadCount { rows, counts })
    }

    fn to_json(&self) -> Value {
        let mut m = kv_json(&self.rows);
        for (k, v) in READ_COUNTS.iter().zip(self.counts) {
            m.insert((*k).to_owned(), Value::Int(v));
        }
        Value::Object(m)
    }
}

// --- variant counts (step 03) --------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Variants {
    pub total: u64,
    pub pass: u64,
    pub snps: u64,
    pub indels: u64,
}

impl Variants {
    /// One pass over the sample VCF, BGZF inflated on many cores. A record
    /// counts with seven cells or more; an ALT of `.`, `*` or `<NON_REF>`
    /// is no allele; a record is a SNP when one ALT and the REF are one
    /// base each, an indel when one ALT that is no symbolic allele differs
    /// in length from the REF. Lengths are in characters, as Python's are.
    pub fn read(path: &str) -> Result<Variants, String> {
        let fail = |e: std::io::Error| format!("{}: {}", basename(path), gzip_message(&e));
        let mut r = og_io::bgzf::open(path, og_io::bgzf::workers()).map_err(fail)?;
        let mut v = Variants::default();
        let mut buf = Vec::with_capacity(1 << 16);
        loop {
            buf.clear();
            if r.read_until(b'\n', &mut buf).map_err(fail)? == 0 {
                break;
            }
            let line = buf.strip_suffix(b"\n").unwrap_or(&buf);
            if memchr::memchr(b'\r', line).is_none() {
                v.record(line);
            } else {
                // Universal newlines: \r\n and a lone \r end a line too.
                let line = line.strip_suffix(b"\r").unwrap_or(line);
                line.split(|&c| c == b'\r').for_each(|l| v.record(l));
            }
        }
        Ok(v)
    }

    fn record(&mut self, line: &[u8]) {
        if line.first() == Some(&b'#') {
            return;
        }
        let mut cells = line.splitn(8, |&c| c == b'\t');
        let (Some(_), Some(_), Some(_), Some(reference), Some(alt), Some(_), Some(filter)) = (
            cells.next(),
            cells.next(),
            cells.next(),
            cells.next(),
            cells.next(),
            cells.next(),
            cells.next(),
        ) else {
            return;
        };
        // The seventh cell ends at the next tab, which splitn left in the rest.
        let filter = filter.split(|&c| c == b'\t').next().unwrap_or(filter);
        self.total += 1;
        if filter == b"PASS" {
            self.pass += 1;
        }
        let ref_len = chars(reference);
        let (mut snp, mut indel) = (false, false);
        for a in alt.split(|&c| c == b',') {
            if matches!(a, b"." | b"*" | b"<NON_REF>") {
                continue;
            }
            let n = chars(a);
            snp |= n == 1 && ref_len == 1;
            indel |= n != ref_len && a.first() != Some(&b'<');
        }
        self.snps += u64::from(snp);
        self.indels += u64::from(indel);
    }

    fn to_json(self) -> Value {
        obj([
            ("total", Value::Int(self.total.into())),
            ("pass", Value::Int(self.pass.into())),
            ("snps", Value::Int(self.snps.into())),
            ("indels", Value::Int(self.indels.into())),
        ])
    }
}

/// The length Python's `len()` gives a cell read with `errors="replace"`.
fn chars(b: &[u8]) -> usize {
    if b.is_ascii() {
        b.len()
    } else {
        String::from_utf8_lossy(b).chars().count()
    }
}

/// What Python's gzip says for the same fault, where the two are the same
/// fault; otherwise the reader's own words.
fn gzip_message(e: &std::io::Error) -> String {
    if e.kind() == std::io::ErrorKind::UnexpectedEof {
        "Compressed file ended before the end-of-stream marker was reached".to_owned()
    } else {
        e.to_string()
    }
}
