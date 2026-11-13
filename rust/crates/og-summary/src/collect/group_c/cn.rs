//! Step 35 (Parascopy and the depth counts, short reads) and step 35b
//! (Paraphase, PacBio HiFi): SMN1 and SMN2 at exon 7 with the SMA reading,
//! the HBA, CYP21A2, interval and gene-dosage tables, every other catalog
//! locus, and the Paraphase genes with SMN1 from both tools.

use super::util::{
    Obj, Table, asset_rows, at, base, columns, count, isdigit, lines, opt_s, or, or_dot, pairs, s,
    strs, table,
};
use super::{Ctx, Found, Out};
use og_cn::markers::{MARKERS, Status};
use og_cn::smn::SMN1_EXON7;
use pycompat::csv::DictRow;
use pycompat::json::Value;

/// SMN2 exon 7, c.840 (the base that tells SMN2 from SMN1), GRCh38 1-based.
pub const SMN2_EXON7: (&str, i64) = ("chr5", 70_076_526);

/// Tables in paralogs/ that other parts of steps 35 and 35b write, not
/// catalog rows.
const NOT_CATALOG_TABLES: [&str; 2] = ["hba", "paraphase"];
/// The columns of a Parascopy copy-number table.
const PARASCOPY_COLUMNS: [&str; 9] = [
    "chrom",
    "start",
    "end",
    "agCN",
    "agCN_filter",
    "agCN_qual",
    "psCN",
    "psCN_qual",
    "psCN_filter",
];
const NO_SMN_TABLE: &str = "no SMN1 table (step 35 wrote no SMN copy-number table)";
const MARKER_FIELDS: [&str; 8] = [
    "marker",
    "rsid",
    "status",
    "alt_reads",
    "depth",
    "fraction",
    "copies",
    "message",
];
const DOSAGE_COLUMNS: [&str; 6] = ["locus", "status", "copies", "estimate", "reason", "reading"];
const INTERVAL_COLUMNS: [&str; 9] = [
    "locus", "name", "chrom", "start", "end", "status", "copies", "expected", "change",
];

/// A copy number read at one position from a Parascopy table.
#[derive(Debug, Clone, PartialEq)]
pub struct CnAt {
    pub cn: Option<i128>,
    pub text: String,
    pub agcn: String,
    pub agcn_qual: String,
    pub agcn_filter: String,
}

impl CnAt {
    fn none(text: &str) -> CnAt {
        CnAt {
            cn: None,
            text: text.to_owned(),
            agcn: ".".into(),
            agcn_qual: ".".into(),
            agcn_filter: ".".into(),
        }
    }

    fn json(&self) -> Value {
        Obj::new()
            .set("cn", cn(self.cn))
            .set("text", s(&self.text))
            .set("agcn", s(&self.agcn))
            .set("agcn_qual", s(&self.agcn_qual))
            .set("agcn_filter", s(&self.agcn_filter))
            .val()
    }
}

fn cn(x: Option<i128>) -> Value {
    x.map_or(Value::Null, Value::Int)
}

/// The region of a row and its homologous regions, 1-based closed; `None`
/// for a row without a whole start and end.
fn regions(get: &dyn Fn(&str) -> Option<String>) -> Option<Vec<(String, i128, i128)>> {
    let chrom = get("chrom")?;
    let start = crate::int(&get("start")?).ok()?;
    let end = crate::int(&get("end")?).ok()?;
    let mut out = vec![(chrom, start + 1, end)];
    let homologous = get("homologous_regions").filter(|h| !h.is_empty());
    for h in homologous.as_deref().unwrap_or("*").split(',') {
        if let Some(r) = region(pycompat::num::py_strip(h)) {
            out.push(r);
        }
    }
    Some(out)
}

/// `^([^:]+):(\d+)-(\d+)` at the start of a text.
fn region(h: &str) -> Option<(String, i128, i128)> {
    let (c, rest) = h.split_once(':')?;
    if c.is_empty() {
        return None;
    }
    let digits = |x: &str| -> Option<(i128, usize)> {
        let n = x.bytes().take_while(u8::is_ascii_digit).count();
        (n > 0)
            .then(|| (x[..n].parse().ok(), n))
            .and_then(|(v, n)| Some((v?, n)))
    };
    let (b, n) = digits(rest)?;
    let rest = rest[n..].strip_prefix('-')?;
    let (e, _) = digits(rest)?;
    Some((c.to_owned(), b, e))
}

/// The paralog-specific copy number (psCN) of the row, or of its
/// homologous region, that holds a position, with the row's aggregate
/// value (agCN). `cn` is `None` unless the filter is PASS, the quality 20
/// or more and the value a whole number (Parascopy's threshold).
fn cn_at_rows<'r>(
    rows: impl Iterator<Item = &'r dyn Fn(&str) -> Option<String>>,
    chrom: &str,
    pos: i64,
    nothing: &str,
) -> CnAt {
    let pos = i128::from(pos);
    for get in rows {
        let Some(regs) = regions(get) else { continue };
        let Some(idx) = regs
            .iter()
            .position(|(c, b, e)| c == chrom && *b <= pos && pos <= *e)
        else {
            continue;
        };
        let val = |k: &str| get(k).filter(|v| !v.is_empty());
        let cns = val("psCN").unwrap_or_default();
        let quals = val("psCN_qual").unwrap_or_default();
        let cns: Vec<&str> = cns.split(',').collect();
        let quals: Vec<&str> = quals.split(',').collect();
        let cn_txt = cns.get(idx).copied().unwrap_or("?");
        let q_txt = quals.get(idx).copied().unwrap_or("*");
        let filt = val("psCN_filter").unwrap_or_else(|| ".".into());
        let q = pycompat::num::parse_float(q_txt).ok();
        let ok = filt == "PASS" && q.is_some_and(|q| q >= 20.0) && isdigit(cn_txt);
        return CnAt {
            cn: if ok { crate::int(cn_txt).ok() } else { None },
            text: format!(
                "{cn_txt} (quality {q_txt}, filter {filt}){}",
                if ok { "" } else { ", not resolved" }
            ),
            agcn: val("agCN").unwrap_or_else(|| ".".into()),
            agcn_qual: val("agCN_qual").unwrap_or_else(|| ".".into()),
            agcn_filter: val("agCN_filter").unwrap_or_else(|| ".".into()),
        };
    }
    CnAt::none(nothing)
}

/// `parascopy_cn_at`: over the rows of step 35's table.
/// A lookup of a row's cell by column.
type Getter<'a> = Box<dyn Fn(&str) -> Option<String> + 'a>;

pub fn parascopy_cn_at(rows: &[DictRow], chrom: &str, pos: i64) -> CnAt {
    let getters: Vec<Getter<'_>> = rows
        .iter()
        .map(|r| {
            Box::new(move |k: &str| r.get(k).flatten().map(str::to_owned))
                as Box<dyn Fn(&str) -> Option<String>>
        })
        .collect();
    cn_at_rows(
        getters
            .iter()
            .map(|g| g.as_ref() as &dyn Fn(&str) -> Option<String>),
        chrom,
        pos,
        "no row covers exon 7",
    )
}

/// `paraphase_report.parascopy_smn1`: SMN1 at exon 7 from step 35's SMN
/// table read line by line, each line zipped with the header.
pub fn parascopy_smn1(path: &str) -> Result<(Option<i128>, String), String> {
    let ls = lines(path)?;
    let mut it = ls.iter();
    let header: Vec<&str> = it.next().map_or_else(Vec::new, |h| h.split('\t').collect());
    let rows: Vec<std::collections::HashMap<&str, &str>> = it
        .map(|l| header.iter().copied().zip(l.split('\t')).collect())
        .collect();
    let getters: Vec<Getter<'_>> = rows
        .iter()
        .map(|r| {
            Box::new(move |k: &str| r.get(k).map(|v| (*v).to_owned()))
                as Box<dyn Fn(&str) -> Option<String>>
        })
        .collect();
    let r = cn_at_rows(
        getters
            .iter()
            .map(|g| g.as_ref() as &dyn Fn(&str) -> Option<String>),
        SMN1_EXON7.0,
        SMN1_EXON7.1,
        "no row covers SMN1 exon 7",
    );
    Ok((r.cn, r.text))
}

/// What an SMN1 copy number means for SMA, in one line.
pub fn smn_reading(smn1: Option<i128>) -> String {
    match smn1 {
        None => {
            "SMN1 not resolved (quality under 20 or filter not PASS): no carrier reading".into()
        }
        Some(0) => "no SMN1 copy: the pattern of spinal muscular atrophy (SMA) itself".into(),
        Some(1) => "one SMN1 copy: the pattern of an SMA carrier".into(),
        Some(2) => "two SMN1 copies: not a carrier by copy number; a 2+0 carrier (both copies on one chromosome) cannot be ruled out".into(),
        Some(n) => format!("{n} SMN1 copies: not a carrier by copy number"),
    }
}

/// A copy-aware marker table (metric, value) beside a count; a table that
/// is not there or has no known status reads as not assessed, with the
/// reason in its message.
pub fn read_marker_table(path: &str, step: &str, key: &str) -> Value {
    let m = MARKERS
        .iter()
        .find(|m| m.key == key)
        .expect("a marker of og-cn");
    let statuses = [
        Status::Present,
        Status::NotSeen,
        Status::Uncertain,
        Status::TooLittleDepth,
        Status::NotAssessed,
    ]
    .map(Status::as_str);
    let mut out: indexmap::IndexMap<&str, String> =
        MARKER_FIELDS.iter().map(|k| (*k, ".".to_owned())).collect();
    out.insert("marker", m.name.to_owned());
    out.insert("rsid", m.dbsnp.0.to_owned());
    out.insert("status", "not assessed".to_owned());
    let name = base(path);
    let done = |out: indexmap::IndexMap<&str, String>| {
        Value::Object(
            out.into_iter()
                .map(|(k, v)| (k.to_owned(), s(&v)))
                .collect(),
        )
    };
    if !crate::is_file(path) {
        out.insert(
            "message",
            format!(
                "{}: not assessed, step {step} wrote no {name} (a run from before this reader, or the check did not run).",
                m.name
            ),
        );
        return done(out);
    }
    let kv: std::collections::HashMap<String, String> =
        pairs(path).unwrap_or_default().into_iter().collect();
    if !kv
        .get("status")
        .is_some_and(|v| statuses.contains(&v.as_str()))
    {
        out.insert(
            "message",
            format!(
                "{}: not assessed, {name} has no status row ({}).",
                m.name,
                statuses.join(", ")
            ),
        );
        return done(out);
    }
    for k in MARKER_FIELDS {
        if let Some(v) = kv.get(k).filter(|v| !v.is_empty()) {
            out.insert(k, v.clone());
        }
    }
    done(out)
}

/// `{table: path}` of step 35's catalog tables (`<s>_<table>_copy_number.tsv`).
fn parascopy_tables(c: &Ctx<'_>) -> Vec<(String, Found)> {
    let pat = format!("paralogs/{}_*_copy_number.tsv", c.s());
    let skip = c.sample.chars().count() + 1;
    c.index
        .all(&pat)
        .into_iter()
        .filter_map(|rel| {
            let b = base(rel);
            let n = b.chars().count();
            let name: String = b
                .chars()
                .skip(skip)
                .take(n.saturating_sub(skip + "_copy_number.tsv".len()))
                .collect();
            (!NOT_CATALOG_TABLES.contains(&name.as_str())).then(|| {
                (
                    name,
                    Found {
                        path: c.path(rel),
                        rel: rel.to_owned(),
                    },
                )
            })
        })
        .collect()
}

/// A table whose header must hold `need`.
fn has_all(cols: &[String], need: &[&str]) -> bool {
    need.iter().all(|k| cols.iter().any(|c| c == k))
}

/// The first `agCN_filter` that starts with "not assessed".
fn not_assessed(rows: &[DictRow]) -> Result<Option<String>, String> {
    for r in rows {
        if r.get("agCN_filter")
            .flatten()
            .is_some_and(|f| f.starts_with("not assessed"))
        {
            return Ok(at(r, "agCN_filter")?.map(str::to_owned));
        }
    }
    Ok(None)
}

/// Step 35 (Parascopy, short reads): SMN1 and SMN2 at exon 7, their
/// aggregate copy number and the SMA reading; the HBA, interval, CYP21A2
/// and dosage tables; and every other catalog table, one locus each.
pub fn parascopy(c: &Ctx<'_>) -> Out {
    let sm = c.sample;
    let tables = parascopy_tables(c);
    let direct = |name: &str| {
        let rel = format!("paralogs/{sm}_{name}");
        Found {
            path: c.path(&rel),
            rel,
        }
    };
    let hp = direct("hba_copy_number.tsv");
    let ip = direct("cn_intervals.tsv");
    let cp = direct("cyp21a2_psv.tsv");
    let dp = direct("cn_dosage.tsv");
    let smn = tables.iter().find(|(n, _)| n == "smn").map(|x| x.1.clone());
    let Some(p) = smn
        .clone()
        .or_else(|| tables.first().map(|x| x.1.clone()))
        .or_else(|| {
            [&hp, &ip, &cp, &dp]
                .into_iter()
                .find(|f| crate::is_file(&f.path))
                .cloned()
        })
    else {
        return Ok(None);
    };
    let (smn1, smn2, reading, nrows) = match &smn {
        Some(f) => {
            let t = table(&f.path, false)?;
            if !has_all(&columns(&f.path)?, &PARASCOPY_COLUMNS) {
                return Err(format!(
                    "{} lacks the chrom, start, end, agCN, agCN_filter, agCN_qual, psCN, psCN_qual or psCN_filter column",
                    base(&f.path)
                ));
            }
            let smn1 = parascopy_cn_at(&t.rows, SMN1_EXON7.0, SMN1_EXON7.1);
            let smn2 = parascopy_cn_at(&t.rows, SMN2_EXON7.0, SMN2_EXON7.1);
            let reading = match not_assessed(&t.rows)? {
                Some(na) => format!("SMN1 and SMN2 {na}: no carrier reading"),
                None if smn1.text == "no row covers exon 7" => {
                    "no row of step 35's table covers SMN1 exon 7: no carrier reading".into()
                }
                None => smn_reading(smn1.cn),
            };
            (smn1, smn2, reading, t.rows.len())
        }
        None => (
            CnAt::none(NO_SMN_TABLE),
            CnAt::none(NO_SMN_TABLE),
            format!("{NO_SMN_TABLE}: no carrier reading"),
            0,
        ),
    };
    let others: Vec<&(String, Found)> = tables.iter().filter(|(n, _)| n != "smn").collect();
    let mut catalog: indexmap::IndexMap<String, DictRow> = indexmap::IndexMap::new();
    if !others.is_empty() {
        for r in asset_rows(c.assets, "paralog_loci.tsv")? {
            let k = at(&r, "table")?.unwrap_or("None").to_owned();
            catalog.insert(k, r);
        }
    }
    let mut order: Vec<&(String, Found)> = others;
    order.sort_by(|a, b| {
        let key = |n: &str| catalog.get_index_of(n).unwrap_or(catalog.len());
        key(&a.0).cmp(&key(&b.0)).then_with(|| a.0.cmp(&b.0))
    });
    let mut loci = Vec::new();
    for (n, f) in order {
        loci.push(parascopy_locus(n, &f.path, catalog.get(n))?);
    }
    let opt = |f: &Found, read: fn(&str) -> Result<Value, String>| -> Result<Value, String> {
        if crate::is_file(&f.path) {
            read(&f.path)
        } else {
            Ok(Value::Null)
        }
    };
    let data = Obj::new()
        .set("smn1", smn1.json())
        .set(
            "smn2",
            Obj::new()
                .set("cn", cn(smn2.cn))
                .set("text", s(&smn2.text))
                .val(),
        )
        .set("agcn", s(&smn1.agcn))
        .set("agcn_qual", s(&smn1.agcn_qual))
        .set("agcn_filter", s(&smn1.agcn_filter))
        .set(
            "smn1_marker",
            read_marker_table(
                &c.path(&format!("paralogs/{sm}_smn1_marker.tsv")),
                "35",
                "SMN1_27134",
            ),
        )
        .set("rows", count(nrows))
        .set("hba", opt(&hp, read_hba)?)
        .set("cn_intervals", opt(&ip, read_cn_intervals)?)
        .set("cyp21a2", opt(&cp, read_hba)?)
        .set("dosage", opt(&dp, read_cn_dosage)?)
        .set("reading", s(&reading))
        .set("loci", Value::Array(loci));
    Ok(Some((p, data.val())))
}

/// One catalog table that is not SMN1's: the locus's aggregate copy number
/// per region of its profile, the catalog's genes, tag and reading.
fn parascopy_locus(name: &str, path: &str, row: Option<&DictRow>) -> Result<Value, String> {
    let cat =
        |k: &str, default: &str| row.map_or_else(|| default.to_owned(), |r| or(r, k, default));
    let out = Obj::new()
        .set("table", s(name))
        .set("genes", s(&cat("genes", name)))
        .set("tag", s(&cat("tag", ".")))
        .set("reading", s(&cat("reading", ".")));
    let t = table(path, false)?;
    let with = |o: Obj, status: &str, cnv: Value, text: &str| {
        o.set("status", s(status))
            .set("cn", cnv)
            .set("text", s(text))
            .val()
    };
    if !has_all(&t.head, &PARASCOPY_COLUMNS) {
        let text = format!("unreadable: {} lacks the copy-number columns", base(path));
        return Ok(with(out, "unreadable", Value::Null, &text));
    }
    if let Some(na) = not_assessed(&t.rows)? {
        return Ok(with(out, "not assessed", Value::Null, &na));
    }
    struct Reg {
        chrom: Option<String>,
        start: Option<String>,
        end: Option<String>,
        cn: String,
        qual: String,
        filter: String,
        ok: bool,
    }
    let regs: Vec<Reg> = t
        .rows
        .iter()
        .map(|r| {
            let get = |k: &str| r.get(k).flatten().map(str::to_owned);
            let q = pycompat::num::parse_float(&or(r, "agCN_qual", "")).ok();
            let cnv = or_dot(r, "agCN");
            let ok = isdigit(&cnv)
                && r.get("agCN_filter").flatten() == Some("PASS")
                && q.is_some_and(|q| q >= 20.0);
            Reg {
                chrom: get("chrom"),
                start: get("start"),
                end: get("end"),
                cn: cnv,
                qual: or_dot(r, "agCN_qual"),
                filter: or_dot(r, "agCN_filter"),
                ok,
            }
        })
        .collect();
    if regs.is_empty() {
        return Ok(with(
            out,
            "not assessed",
            Value::Null,
            "not assessed: the table has no region row",
        ));
    }
    let values: std::collections::HashSet<&str> = regs
        .iter()
        .filter(|r| r.ok)
        .map(|r| r.cn.as_str())
        .collect();
    if values.len() == 1 && regs.iter().all(|r| r.ok) {
        let n = crate::int(values.into_iter().next().unwrap_or_default())
            .map_err(|e| super::super::sections::py_error("ValueError", &e))?;
        let mut min = f64::INFINITY;
        for r in &regs {
            min = min.min(
                crate::float(&r.qual)
                    .map_err(|e| super::super::sections::py_error("ValueError", &e))?,
            );
        }
        let text = format!(
            "{n} ({} region(s), quality {} or more, filter PASS)",
            regs.len(),
            pycompat::fmt::general(min, 6)
        );
        return Ok(with(out, "called", Value::Int(n), &text));
    }
    let none = |x: &Option<String>| x.clone().unwrap_or_else(|| "None".into());
    let text = regs
        .iter()
        .map(|r| {
            let b = match &r.start {
                Some(b) if isdigit(b) => {
                    crate::int(b).map_or_else(|_| b.clone(), |v| (v + 1).to_string())
                }
                b => none(b),
            };
            format!(
                "{}:{b}-{} {} (quality {}, filter {})",
                none(&r.chrom),
                none(&r.end),
                r.cn,
                r.qual,
                r.filter
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    let called = regs.iter().any(|r| r.ok);
    let text = format!(
        "{text}{}",
        if called {
            ""
        } else {
            ": no region with a whole copy number, filter PASS and quality 20 or more"
        }
    );
    let out = out.set("reading", s("."));
    Ok(with(
        out,
        if called { "called" } else { "not resolved" },
        Value::Null,
        &text,
    ))
}

/// Step 35's interval table: the rows with a deletion or a duplication, and
/// a count of the rows by status.
fn read_cn_intervals(path: &str) -> Result<Value, String> {
    let t = table(path, false)?;
    if !has_all(&t.head, &INTERVAL_COLUMNS) {
        return Err(format!(
            "{} lacks the locus, name, status, copies, expected or change column",
            base(path)
        ));
    }
    short_row(path, &t, &INTERVAL_COLUMNS)?;
    let mut by: indexmap::IndexMap<String, i128> = indexmap::IndexMap::new();
    let mut changes = Vec::new();
    for r in &t.rows {
        let cell = |k: &str| r.get(k).flatten().unwrap_or_default().to_owned();
        *by.entry(cell("status")).or_default() += 1;
        if matches!(cell("change").as_str(), "deletion" | "duplication") {
            let keep = [
                "locus", "name", "chrom", "start", "end", "copies", "expected", "change",
            ];
            changes.push(
                keep.iter()
                    .fold(Obj::new(), |o, k| o.set(k, s(&cell(k))))
                    .val(),
            );
        }
    }
    Ok(Obj::new()
        .set("rows", count(t.rows.len()))
        .set(
            "by_status",
            Value::Object(by.into_iter().map(|(k, v)| (k, Value::Int(v))).collect()),
        )
        .set("changes", Value::Array(changes))
        .val())
}

/// Step 35's gene-dosage table: one entry per locus, in the table's order.
fn read_cn_dosage(path: &str) -> Result<Value, String> {
    if !has_all(&columns(path)?, &DOSAGE_COLUMNS) {
        return Err(format!(
            "{} lacks the {} or reading column",
            base(path),
            DOSAGE_COLUMNS[..5].join(", ")
        ));
    }
    let t = table(path, false)?;
    short_row(path, &t, &DOSAGE_COLUMNS)?;
    Ok(Value::Array(
        t.rows
            .iter()
            .map(|r| super::util::dots(r, &DOSAGE_COLUMNS))
            .collect(),
    ))
}

/// A row cut short of a column makes the table unreadable.
fn short_row(path: &str, t: &Table, cols: &[&str]) -> Result<(), String> {
    match t
        .rows
        .iter()
        .position(|r| cols.iter().any(|k| !matches!(r.get(k), Some(Some(_)))))
    {
        Some(i) => Err(format!("{} line {} lacks a field", base(path), i + 2)),
        None => Ok(()),
    }
}

/// Step 35's HBA table (metric, value), and its CYP21A2 table the same way;
/// a table without a status row reads as not assessed.
fn read_hba(path: &str) -> Result<Value, String> {
    let mut out: indexmap::IndexMap<String, Value> = indexmap::IndexMap::new();
    for l in lines(path)? {
        if let Some((k, v)) = l.split_once('\t')
            && k != "metric"
        {
            out.insert(k.to_owned(), s(v));
        }
    }
    let ok =
        matches!(out.get("status"), Some(Value::Str(v)) if v == "called" || v == "not assessed");
    if !ok {
        return Ok(Obj::new()
            .set("status", s("not assessed"))
            .set(
                "message",
                s(&format!("not assessed: {} has no status row", base(path))),
            )
            .val());
    }
    Ok(Value::Object(out))
}

/// The Paraphase regions shown gene by gene, in the order of
/// assets/paraphase_genes.tsv.
fn paralog_regions(assets: &str) -> Result<Vec<String>, String> {
    let mut out: Vec<String> = Vec::new();
    for r in asset_rows(assets, "paraphase_genes.tsv")? {
        let g = at(&r, "region")?.unwrap_or("None").to_owned();
        if !out.contains(&g) {
            out.push(g);
        }
    }
    Ok(out)
}

/// Step 35b (Paraphase, PacBio HiFi): the genes of the catalog regions, a
/// count of the other regions, the variant rows per region, and SMN1
/// beside step 35's Parascopy call.
pub fn paralogs(c: &Ctx<'_>) -> Out {
    let s_ = c.sample;
    let Some(p) = c.find1(format!("paralogs/{}_paraphase_copy_number.tsv", c.s())) else {
        return Ok(None);
    };
    let t = table(&p.path, false)?;
    if !t.rows.is_empty() && !has_all(&t.head, &["region", "gene", "copy_number", "status"]) {
        return Err(format!(
            "{} lacks the region, gene, copy_number or status column",
            base(&p.path)
        ));
    }
    let regions = paralog_regions(c.assets)?;
    let region_of = |r: &DictRow| r.get("region").flatten().map(str::to_owned);
    let mut genes: Vec<(usize, Value, Option<String>, Option<String>)> = Vec::new();
    let mut others = 0;
    let mut others_called = 0;
    for r in &t.rows {
        let reg = region_of(r);
        match reg
            .as_ref()
            .and_then(|g| regions.iter().position(|x| x == g))
        {
            Some(i) => {
                let cell = |k: &str| opt_s(r.get(k).flatten());
                genes.push((
                    i,
                    Obj::new()
                        .set("region", cell("region"))
                        .set("gene", cell("gene"))
                        .set("copy_number", cell("copy_number"))
                        .set("status", cell("status"))
                        .set("detail", s(&or_dot(r, "detail")))
                        .val(),
                    r.get("gene").flatten().map(str::to_owned),
                    r.get("copy_number").flatten().map(str::to_owned),
                ));
            }
            None => {
                others += 1;
                if r.get("status").flatten() == Some("called") {
                    others_called += 1;
                }
            }
        }
    }
    genes.sort_by_key(|g| g.0);
    let mut variants: indexmap::IndexMap<String, i128> = indexmap::IndexMap::new();
    if let Some(vp) = c.beside(&p, &format!("{s_}_paraphase_variants.tsv")) {
        for r in table(&vp.path, false)?.rows {
            let k = match r.get("region") {
                None => String::new(),
                Some(None) => "null".into(),
                Some(Some(v)) => v.to_owned(),
            };
            *variants.entry(k).or_default() += 1;
        }
    }
    let smn1 = genes
        .iter()
        .find(|g| g.2.as_deref() == Some("SMN1"))
        .and_then(|g| g.3.clone())
        .filter(|x| isdigit(x))
        .and_then(|x| crate::int(&x).ok());
    let parascopy = match c.beside(&p, &format!("{s_}_smn_copy_number.tsv")) {
        Some(pp) => Some(parascopy_smn1(&pp.path)?),
        None => None,
    };
    let agree = match (smn1, &parascopy) {
        (Some(a), Some((Some(b), _))) => s(if a == *b { "agree" } else { "disagree" }),
        _ => Value::Null,
    };
    let parascopy_json = parascopy.map_or(Value::Null, |(n, text)| {
        Obj::new().set("cn", cn(n)).set("text", s(&text)).val()
    });
    let data = Obj::new()
        .set(
            "genes",
            Value::Array(genes.into_iter().map(|g| g.1).collect()),
        )
        .set("other_regions", count(others))
        .set("other_regions_called", count(others_called))
        .set(
            "variants",
            Value::Object(
                variants
                    .into_iter()
                    .map(|(k, v)| (k, Value::Int(v)))
                    .collect(),
            ),
        )
        .set(
            "smn1",
            Obj::new()
                .set("paraphase", cn(smn1))
                .set("parascopy", parascopy_json)
                .set("agree", agree)
                .val(),
        );
    let _ = strs::<[&str; 0], &str>;
    Ok(Some((p, data.val())))
}
