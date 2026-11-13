//! Step 35b's two ends around Paraphase (`docs/35-paralogs.md`).
//!
//! ```text
//! og paraphase check-bam BAM
//! og paraphase table --sample S --json S.paraphase.json [--vcf-dir DIR]
//!                    --out-cn CN.tsv --out-variants VARIANTS.tsv [--genes FILE]
//! ```
//!
//! `check-bam` reads the BAM's header and exits 0 only when every read
//! group says `PL:PACBIO`. Paraphase is written for PacBio HiFi reads, and
//! its licence (the PacBio Software License Agreement, from Paraphase 3.3.0
//! on) allows it only on data from a PacBio instrument. Another BAM stops
//! here with exit 3 and the reason; an unreadable one with exit 1.
//!
//! `table` turns Paraphase's JSON into one row per gene (copy number, call
//! status, what the region-specific fields say) and its per-region VCFs
//! into one row per variant a haplotype carries, with the reads behind it.
//! The regions of `assets/paraphase_genes.tsv` (or `--genes`) get the rows
//! that catalog names, each from the JSON field it names; every other
//! region gets one row. Both tables are written together or not at all.
//!
//! This is `bin/paraphase_report.py`'s `check-bam` and `table` rewritten.
//! On purpose: the header read is bounded (a corrupt `l_text` cannot make
//! it read the whole BAM), and a JSON object, list or boolean where a
//! number belongs stops the table with `<region>.<field> is not a number`
//! instead of printing Python syntax into a cell.

use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use indexmap::IndexMap;
use og_io::{Lines, Utf8};
use pycompat::json::Value;
use pycompat::value::{py_str, truthy};

use crate::pgx_rules::{self, Error, TextLines, strip};

/// The copy-number table's columns.
pub const CN_HEADER: [&str; 8] = [
    "region",
    "gene",
    "copy_number",
    "status",
    "detail",
    "total_cn",
    "region_depth",
    "genome_depth",
];
/// The variant table's columns.
pub const VARIANT_HEADER: [&str; 8] = [
    "region",
    "chrom",
    "pos",
    "ref",
    "alt",
    "filter",
    "haplotypes",
    "support",
];
/// The region catalog's columns.
pub const GENES_COLUMNS: [&str; 5] = ["region", "gene", "field", "reading", "source"];
/// Fields read by code written for one region's JSON, and that region.
const REGION_FIELDS: [(&str, &str); 2] = [("genotype", "hba"), ("annotated_alleles", "rccx")];
/// The largest SAM header text a BAM may announce: 1 GiB.
const MAX_HEADER: i32 = 1 << 30;

const WHOLE_FAMILY: &str = "copies of the whole gene family (gene and paralogs together)";
const HBA_LARGE: &str = "an allele with both genes lost (--) is a large deletion such as --SEA, --MED or --FIL; its breakpoints lie outside the region Paraphase reads, so which one it is cannot be told";

/// The HBA structural variant Paraphase names from a hybrid haplotype.
fn hba_sv(kind: &str) -> &str {
    match kind {
        "3p7del" => {
            "3.7 kb deletion (-a3.7): HBA2 and HBA1 fused into one gene, one alpha-globin gene lost"
        }
        "4p2del" => "4.2 kb deletion (-a4.2): HBA2 lost",
        "3p7dup" => "3.7 kb duplication (aaa anti-3.7): one extra alpha-globin gene",
        "4p2dup" => "4.2 kb duplication (aaa anti-4.2): one extra HBA2",
        other => other,
    }
}

// --- check-bam ---------------------------------------------------------------

/// The SAM header text of a BAM, non-ASCII bytes read as U+FFFD.
pub fn bam_header(path: &str) -> Result<String, Error> {
    let fail = |e: io::Error| Error(format!("{path}: {e}"));
    let mut r = og_io::open(path).map_err(|e| Error(og_filters::os_error(&e, path)))?;
    let mut magic = Vec::with_capacity(4);
    (&mut r).take(4).read_to_end(&mut magic).map_err(fail)?;
    if magic != b"BAM\x01" {
        return Err(Error(format!("{path} is not a BAM file")));
    }
    let mut len = [0u8; 4];
    r.read_exact(&mut len).map_err(fail)?;
    let l_text = i32::from_le_bytes(len);
    if !(0..=MAX_HEADER).contains(&l_text) {
        return Err(Error(format!(
            "{path}: not a BAM header (it announces {l_text} bytes of text)"
        )));
    }
    let mut text = Vec::new();
    r.take(u64::from(l_text.unsigned_abs()))
        .read_to_end(&mut text)
        .map_err(fail)?;
    Ok(text
        .iter()
        .map(|&b| if b.is_ascii() { b as char } else { '\u{fffd}' })
        .collect())
}

/// `(read group ID, PL)` of the header's `@RG` lines; PL empty when absent.
pub fn platforms(header: &str) -> Vec<(String, String)> {
    header
        .split('\n')
        .filter(|l| l.starts_with("@RG"))
        .map(|line| {
            let (mut id, mut pl) = ("?".to_owned(), String::new());
            for tag in line.split('\t').skip(1) {
                match tag.split_once(':') {
                    Some(("ID", v)) => v.clone_into(&mut id),
                    Some(("PL", v)) => v.clone_into(&mut pl),
                    _ => {}
                }
            }
            (id, pl)
        })
        .collect()
}

/// `Ok(Ok(message))` when the BAM has read groups and all are `PL:PACBIO`,
/// `Ok(Err(reason))` when it has none or another platform; `Err` when the
/// BAM cannot be read.
pub fn check_bam(path: &str) -> Result<Result<String, String>, Error> {
    let rgs = platforms(&bam_header(path)?);
    if rgs.is_empty() {
        return Ok(Err(format!(
            "{path} has no @RG line, so nothing says the reads are PacBio HiFi. Paraphase needs a BAM of PacBio HiFi reads with PL:PACBIO in each read group (step 02b with PLATFORM=hifi writes one, as pbmm2 does)."
        )));
    }
    let mut other: Vec<String> = rgs
        .iter()
        .filter(|(_, pl)| pycompat::text::upper(pl) != "PACBIO")
        .map(|(_, pl)| {
            if pl.is_empty() {
                "none".to_owned()
            } else {
                pl.clone()
            }
        })
        .collect();
    other.sort();
    other.dedup();
    if !other.is_empty() {
        return Ok(Err(format!(
            "{path} has read groups with PL {}, not PACBIO. Paraphase is written for PacBio HiFi reads, and its licence (PacBio Software License Agreement) allows it only on data from a PacBio instrument. For short reads, step 35 (Parascopy) estimates SMN1 and SMN2 copy number.",
            other.join(", ")
        )));
    }
    Ok(Ok(format!(
        "{path}: {} read group(s), all PL:PACBIO",
        rgs.len()
    )))
}

// --- the region catalog --------------------------------------------------------

/// The region catalog: per region, its `(gene, field)` rows in file order.
pub type Genes = Vec<(String, Vec<(String, String)>)>;

fn list_repr(items: &[&str]) -> String {
    let inner: Vec<String> = items.iter().map(|s| og_filters::py_repr(s)).collect();
    format!("[{}]", inner.join(", "))
}

/// Reads `assets/paraphase_genes.tsv`.
pub fn load_genes(path: &str) -> Result<Genes, Error> {
    let lines = TextLines::open(path, Utf8::Strict)
        .map_err(|e| Error(format!("the region catalog {path}: {e}")))?;
    let mut rows: Vec<String> = Vec::new();
    for line in lines {
        let line = line.map_err(|e| {
            Error(format!(
                "the region catalog {path}: {}",
                Error::io(&e, path)
            ))
        })?;
        if !strip(&line).is_empty() && !line.starts_with('#') {
            rows.push(line.strip_suffix('\n').unwrap_or(&line).to_owned());
        }
    }
    let cells: Vec<Vec<&str>> = rows.iter().map(|r| r.split('\t').collect()).collect();
    match cells.first() {
        Some(h) if h.as_slice() == GENES_COLUMNS => {}
        h => {
            let got = h.map_or_else(|| "None".to_owned(), |h| list_repr(h));
            return Err(Error(format!(
                "{path}: header {got}, want {}",
                list_repr(&GENES_COLUMNS)
            )));
        }
    }
    let mut out: Genes = Vec::new();
    for (n, c) in cells.iter().enumerate().skip(1) {
        let n = n + 1;
        if c.len() != GENES_COLUMNS.len() || c.iter().any(|x| strip(x).is_empty()) {
            return Err(Error(format!(
                "{path}: row {n} does not have the five columns filled"
            )));
        }
        if let Some((field, only)) = REGION_FIELDS.iter().find(|(f, _)| *f == c[2])
            && *only != c[0]
        {
            return Err(Error(format!(
                "{path}: row {n} reads field {field} on region {}; that field is read only on region {only}",
                c[0]
            )));
        }
        let entry = (c[1].to_owned(), c[2].to_owned());
        match out.iter_mut().find(|(r, _)| r == c[0]) {
            Some((_, e)) => e.push(entry),
            None => out.push((c[0].to_owned(), vec![entry])),
        }
    }
    Ok(out)
}

/// `assets/paraphase_genes.tsv` beside the running og: the first folder
/// above it that holds it (`bin/.og/og-linux-*` in a checkout,
/// `/pgp-bin/.og/...` beside `/pgp-assets` does not count: the step passes
/// `--genes`).
fn default_genes() -> Option<PathBuf> {
    let exe = std::env::current_exe().ok()?;
    exe.ancestors()
        .skip(1)
        .map(|d| d.join("assets").join("paraphase_genes.tsv"))
        .find(|p| p.is_file())
}

// --- reading Paraphase's JSON ----------------------------------------------------

type Object = IndexMap<String, Value>;

/// One region of the JSON, with its name for messages.
struct Region<'a> {
    name: &'a str,
    v: &'a Object,
}

/// A copy number as its cell prints it, and whether there is one.
struct Cn {
    text: String,
    some: bool,
}

impl Cn {
    fn count(n: usize) -> Cn {
        Cn {
            text: n.to_string(),
            some: true,
        }
    }

    fn none() -> Cn {
        Cn {
            text: ".".to_owned(),
            some: false,
        }
    }

    fn status(&self) -> &'static str {
        if self.some { "called" } else { "no call" }
    }
}

/// A number cell: `.` for nothing or an empty string, an int for an
/// integral float. An object, a list or a boolean is refused.
fn number(v: Option<&Value>, region: &str, field: &str) -> Result<Cn, Error> {
    let text = match v {
        None | Some(Value::Null) => return Ok(Cn::none()),
        Some(Value::Str(s)) if s.is_empty() => ".".to_owned(),
        Some(Value::Str(s)) => s.clone(),
        Some(Value::Int(i)) => i.to_string(),
        Some(Value::Float(x)) if x.is_finite() && x.fract() == 0.0 => pycompat::fmt::int_str(*x),
        Some(Value::Float(x)) => pycompat::fmt::repr(*x),
        Some(_) => return Err(Error(format!("{region}.{field} is not a number"))),
    };
    Ok(Cn { text, some: true })
}

/// `str(v)`, as an error when Python's would be one this port refuses.
fn text_of(v: &Value, region: &str, field: &str) -> Result<String, Error> {
    py_str(v).map_err(|e| Error(format!("{region}.{field}: {e}")))
}

impl<'a> Region<'a> {
    fn get(&self, key: &str) -> Option<&'a Value> {
        self.v.get(key)
    }

    /// `region.get(key) or default`, as present and true.
    fn truthy(&self, key: &str) -> Option<&'a Value> {
        self.get(key).filter(|v| truthy(v))
    }

    fn not_a(&self, field: &str, what: &str) -> Error {
        Error(format!("{}.{field} is not {what}", self.name))
    }

    /// `region.get("region_specific_info") or {}`.
    fn rsi(&self) -> Result<Option<&'a Object>, Error> {
        match self.truthy("region_specific_info") {
            None => Ok(None),
            Some(Value::Object(o)) => Ok(Some(o)),
            Some(_) => Err(self.not_a("region_specific_info", "an object")),
        }
    }

    fn rsi_get(&self, key: &str) -> Result<Option<&'a Value>, Error> {
        Ok(self.rsi()?.and_then(|o| o.get(key)))
    }

    fn rsi_has(&self, key: &str) -> Result<bool, Error> {
        Ok(self.rsi()?.is_some_and(|o| o.contains_key(key)))
    }

    /// `rsi.get(key) or {}` that must be an object, items sorted by key.
    fn rsi_items(&self, key: &str) -> Result<Vec<(&'a String, &'a Value)>, Error> {
        match self.rsi_get(key)?.filter(|v| truthy(v)) {
            None => Ok(Vec::new()),
            Some(Value::Object(o)) => {
                let mut items: Vec<_> = o.iter().collect();
                items.sort_by(|a, b| a.0.cmp(b.0));
                Ok(items)
            }
            Some(_) => Err(self.not_a(key, "an object")),
        }
    }

    /// `str()` of each item Python iterates: a list's items, an object's
    /// keys, a string's characters.
    fn iter_texts(&self, v: &Value, field: &str) -> Result<Vec<String>, Error> {
        match v {
            Value::Array(a) => a.iter().map(|x| text_of(x, self.name, field)).collect(),
            Value::Object(o) => Ok(o.keys().cloned().collect()),
            Value::Str(s) => Ok(s.chars().map(String::from).collect()),
            _ => Err(self.not_a(field, "a list")),
        }
    }

    fn cn(&self, key: &str) -> Result<Cn, Error> {
        number(self.get(key), self.name, key)
    }

    /// Names of the final haplotypes, each two-copy haplotype counted twice.
    fn hap_names(&self) -> Result<Vec<&'a str>, Error> {
        let field = "final_haplotypes";
        let values: Vec<&Value> = match self.truthy(field) {
            None => Vec::new(),
            Some(Value::Object(o)) => o.values().collect(),
            Some(Value::Array(a)) => a.iter().collect(),
            Some(_) => return Err(self.not_a(field, "a list or an object")),
        };
        let mut names: Vec<&str> = values
            .into_iter()
            .map(|v| match v {
                Value::Str(s) => Ok(s.as_str()),
                _ => Err(self.not_a(field, "a list of haplotype names")),
            })
            .collect::<Result<_, _>>()?;
        let twice: Vec<String> = match self.truthy("two_copy_haplotypes") {
            None => Vec::new(),
            Some(Value::Array(a)) => a
                .iter()
                .filter_map(|v| match v {
                    Value::Str(s) => Some(s.clone()),
                    _ => None,
                })
                .collect(),
            Some(Value::Object(o)) => o.keys().cloned().collect(),
            Some(Value::Str(s)) => s.chars().map(String::from).collect(),
            Some(_) => return Err(self.not_a("two_copy_haplotypes", "a list")),
        };
        let extra: Vec<&str> = twice
            .iter()
            .filter_map(|h| names.iter().find(|n| **n == h.as_str()).copied())
            .collect();
        names.extend(extra);
        Ok(names)
    }
}

/// One row of the copy-number table before its region columns.
struct GeneRow {
    gene: String,
    cn: Cn,
    status: &'static str,
    detail: String,
}

impl GeneRow {
    fn new(gene: &str, cn: Cn, detail: String) -> GeneRow {
        let status = cn.status();
        GeneRow {
            gene: gene.to_owned(),
            cn,
            status,
            detail,
        }
    }

    fn not_assessed(gene: &str, field: &str) -> GeneRow {
        GeneRow {
            gene: gene.to_owned(),
            cn: Cn::none(),
            status: "not assessed",
            detail: format!("not assessed: Paraphase wrote no {field} for this region"),
        }
    }
}

fn hba_rows(r: &Region) -> Result<Vec<GeneRow>, Error> {
    let names = r.hap_names()?;
    let count = |tag: &str| names.iter().filter(|h| h.contains(tag)).count();
    let (hba1, hba2) = (count("hba1hap"), count("hba2hap"));
    let genotype = match r.rsi_get("genotype")?.filter(|v| truthy(v)) {
        None => None,
        Some(Value::Str(s)) => Some(s.as_str()),
        Some(_) => return Err(r.not_a("genotype", "a string")),
    };
    let mut named = Vec::new();
    for (hap, val) in r.rsi_items("sv_called")? {
        let kind = text_of(val, r.name, "sv_called")?;
        let kind = kind.split(',').next().unwrap_or_default();
        named.push(format!("{} ({hap})", hba_sv(kind)));
    }
    let mut detail = Vec::new();
    match genotype {
        Some(g) => {
            let alpha = g.matches('a').count();
            detail.push(format!("genotype {g} ({alpha} alpha-globin genes)"));
            if g.contains("--") {
                detail.push(HBA_LARGE.to_owned());
            }
        }
        None => detail.push("no genotype: Paraphase could not resolve the alleles".to_owned()),
    }
    let none_named = named.is_empty();
    detail.extend(named);
    if none_named && genotype.is_some() {
        detail.push("no 3.7 kb or 4.2 kb deletion or duplication found".to_owned());
    }
    let status = if genotype.is_some() {
        "called"
    } else {
        "no call"
    };
    let detail = detail.join("; ");
    let row = |gene: &str, n: usize| GeneRow {
        gene: gene.to_owned(),
        cn: if names.is_empty() {
            Cn::none()
        } else {
            Cn::count(n)
        },
        status,
        detail: detail.clone(),
    };
    Ok(vec![row("HBA1", hba1), row("HBA2", hba2)])
}

fn rccx_row(r: &Region, gene: &str) -> Result<GeneRow, Error> {
    let field = "annotated_alleles";
    let alleles = match r.rsi_get(field)?.filter(|v| truthy(v)) {
        None => Vec::new(),
        Some(Value::Object(o)) => o
            .values()
            .map(|v| text_of(v, r.name, field))
            .collect::<Result<_, _>>()?,
        Some(v) => r.iter_texts(v, field)?,
    };
    let mut detail =
        vec!["RCCX module copies (each carries CYP21A2 or its pseudogene CYP21A1P)".to_owned()];
    if !alleles.is_empty() {
        detail.push(format!("CYP21A2 alleles: {}", alleles.join(" / ")));
    }
    if matches!(r.rsi_get("phasing_success")?, Some(Value::Bool(false))) {
        detail.push("the copies could not be phased into two alleles".to_owned());
    }
    Ok(GeneRow::new(gene, r.cn("total_cn")?, detail.join("; ")))
}

fn generic_row(r: &Region) -> Result<GeneRow, Error> {
    let genes: Vec<String> = match r.truthy("genes_in_region") {
        None => vec![r.name.to_owned()],
        Some(Value::Str(s)) => s.split(',').map(str::to_owned).collect(),
        Some(Value::Array(a)) => a
            .iter()
            .map(|v| match v {
                Value::Str(s) => Ok(s.clone()),
                _ => Err(r.not_a("genes_in_region", "a list of gene names")),
            })
            .collect::<Result<_, _>>()?,
        Some(_) => return Err(r.not_a("genes_in_region", "a list of gene names")),
    };
    let mut gene = genes.join(",");
    let mut detail = Vec::new();
    if let Some(f) = r.rsi_get("fusions_called")?.filter(|v| truthy(v)) {
        let opts = pycompat::json::Dump {
            sort_keys: true,
            ..Default::default()
        };
        let text = pycompat::json::dumps(f, &opts)
            .map_err(|e| Error(format!("{}.fusions_called: {e}", r.name)))?;
        detail.push(format!("fusions: {text}"));
    }
    let names = r.hap_names()?;
    let cn = if r.rsi_has("gene_cn")? {
        detail.insert(0, "copies of the gene, not of its pseudogenes".to_owned());
        number(r.rsi_get("gene_cn")?, r.name, "gene_cn")?
    } else if names
        .iter()
        .any(|h| h.contains("_gene1hap") || h.contains("_gene2hap"))
    {
        // Paraphase labels each haplotype as the region's first gene
        // (gene1) or its second gene or pseudogene (gene2): GBA1/GBAP1,
        // OTOA/OTOAP1. One it cannot place is "unknownhap".
        gene.clone_from(&genes[0]);
        let other = genes.get(1).map_or("pseudogene", String::as_str);
        let count = |tag: &str| names.iter().filter(|h| h.contains(tag)).count();
        let unknown = count("_unknownhap");
        let mut d = format!(
            "copies of {gene} alone (Paraphase's gene1 haplotypes); {other} copies {} (gene2 haplotypes)",
            count("_gene2hap")
        );
        if unknown > 0 {
            let _ = write!(d, "; {unknown} haplotype(s) not placed on either");
        }
        detail.insert(0, d);
        Cn::count(count("_gene1hap"))
    } else {
        detail.insert(0, WHOLE_FAMILY.to_owned());
        r.cn("total_cn")?
    };
    Ok(GeneRow::new(&gene, cn, detail.join("; ")))
}

fn number_row(r: &Region, gene: &str, field: &str) -> Result<GeneRow, Error> {
    if !r.rsi_has(field)? {
        return Ok(GeneRow::not_assessed(gene, field));
    }
    let cn = number(r.rsi_get(field)?, r.name, field)?;
    let mut detail = if cn.some {
        String::new()
    } else {
        "one haplotype seen and the depth does not settle one copy or two".to_owned()
    };
    if let Some(dels) = r.rsi_get("deletion_haplotypes")?.filter(|v| truthy(v)) {
        let dels = r.iter_texts(dels, "deletion_haplotypes")?.join(", ");
        let with = format!("haplotypes with the 11.7 kb deletion: {dels}");
        detail = if detail.is_empty() {
            with
        } else {
            format!("{detail}; {with}")
        };
    }
    Ok(GeneRow::new(gene, cn, detail))
}

fn sv_row(r: &Region, gene: &str) -> Result<GeneRow, Error> {
    if !r.rsi_has("sv_called")? {
        return Ok(GeneRow::not_assessed(gene, "sv_called"));
    }
    let svs = r
        .rsi_items("sv_called")?
        .into_iter()
        .map(|(hap, kind)| Ok(format!("{} ({hap})", text_of(kind, r.name, "sv_called")?)))
        .collect::<Result<Vec<_>, Error>>()?;
    let detail = if svs.is_empty() {
        "no intron 22 inversion or h1/h2 deletion found".to_owned()
    } else {
        svs.join("; ")
    };
    Ok(GeneRow::new(gene, r.cn("total_cn")?, detail))
}

fn repeat_row(r: &Region, gene: &str) -> Result<GeneRow, Error> {
    if !r.rsi_has("repeat_name")? {
        return Ok(GeneRow::not_assessed(gene, "repeat_name"));
    }
    let mut parts = Vec::new();
    for (k, v) in r.rsi_items("repeat_name")? {
        let n = match v {
            Value::Array(a) => a.len().to_string(),
            v => text_of(v, r.name, "repeat_name")?,
        };
        parts.push(format!("{}: {n} haplotype(s)", pycompat::text::upper(k)));
    }
    let detail = if parts.is_empty() {
        "no haplotype assigned to a repeat copy".to_owned()
    } else {
        parts.join("; ")
    };
    Ok(GeneRow::new(gene, r.cn("total_cn")?, detail))
}

fn fusion_row(r: &Region, gene: &str) -> Result<GeneRow, Error> {
    let field = "fusions_called";
    if !r.rsi_has(field)? {
        return Ok(GeneRow::not_assessed(gene, field));
    }
    let mut parts = Vec::new();
    for (hap, v) in r.rsi_items(field)? {
        let kind = match v {
            Value::Object(o) => match o.get("type").filter(|t| truthy(t)) {
                Some(t) => text_of(t, r.name, field)?,
                None => "fusion".to_owned(),
            },
            v => text_of(v, r.name, field)?,
        };
        parts.push(format!("{kind} ({hap})"));
    }
    let mut detail = vec![WHOLE_FAMILY.to_owned()];
    detail.push(if parts.is_empty() {
        "no deletion or duplication between the copies found".to_owned()
    } else {
        parts.join("; ")
    });
    Ok(GeneRow::new(gene, r.cn("total_cn")?, detail.join("; ")))
}

/// The rows of one region: its catalog rows, or one generic row.
fn gene_rows(r: &Region, genes: &Genes) -> Result<Vec<GeneRow>, Error> {
    let Some((_, entries)) = genes
        .iter()
        .find(|(name, e)| name == r.name && !e.is_empty())
    else {
        return Ok(vec![generic_row(r)?]);
    };
    let mut rows = Vec::new();
    for (gene, field) in entries {
        match field.as_str() {
            // HBA counts its genes from the haplotype names and RCCX reads
            // the module's copies from total_cn, so neither is "not
            // assessed" without its field.
            "genotype" => rows.extend(hba_rows(r)?.into_iter().filter(|x| x.gene == *gene)),
            "annotated_alleles" => rows.push(rccx_row(r, gene)?),
            "final_haplotypes" => {
                let g = generic_row(r)?;
                rows.push(GeneRow {
                    gene: gene.clone(),
                    ..g
                });
            }
            // The whole family, even when Paraphase also writes a gene_cn.
            "total_cn" => rows.push(GeneRow::new(
                gene,
                r.cn("total_cn")?,
                WHOLE_FAMILY.to_owned(),
            )),
            "sv_called" => rows.push(sv_row(r, gene)?),
            "repeat_name" => rows.push(repeat_row(r, gene)?),
            "fusions_called" => rows.push(fusion_row(r, gene)?),
            _ => rows.push(number_row(r, gene, field)?),
        }
    }
    Ok(rows)
}

/// The copy-number table of Paraphase's JSON: regions whose rows come
/// from a field of their own first, then by name, case folded.
pub fn cn_table(data: &Object, genes: &Genes) -> Result<Vec<[String; 8]>, Error> {
    let first = |name: &str| {
        genes.iter().any(|(r, e)| {
            r == name
                && e.iter()
                    .any(|(_, f)| f != "final_haplotypes" && f != "total_cn")
        })
    };
    let mut names: Vec<(bool, String, &String)> = data
        .keys()
        .map(|k| (!first(k), pycompat::text::lower(k), k))
        .collect();
    names.sort_by(|a, b| (a.0, &a.1).cmp(&(b.0, &b.1)));
    let mut out = Vec::new();
    for (_, _, name) in names {
        let Value::Object(v) = &data[name] else {
            continue;
        };
        let r = Region { name, v };
        let total = r.cn("total_cn")?.text;
        let depth = match r.get("region_depth") {
            Some(Value::Object(o)) => number(o.get("median"), name, "region_depth")?,
            Some(Value::Array(a)) if !a.is_empty() => number(a.first(), name, "region_depth")?,
            d => number(d, name, "region_depth")?,
        }
        .text;
        let genome = r.cn("genome_depth")?.text;
        for g in gene_rows(&r, genes)? {
            let detail = if g.detail.is_empty() {
                ".".to_owned()
            } else {
                g.detail
            };
            out.push([
                name.clone(),
                g.gene,
                g.cn.text,
                g.status.to_owned(),
                detail,
                total.clone(),
                depth.clone(),
                genome.clone(),
            ]);
        }
    }
    Ok(out)
}

/// Reads Paraphase's JSON: an object of regions.
pub fn load_json(path: &str) -> Result<Object, Error> {
    let base = Path::new(path)
        .file_name()
        .map_or_else(|| path.to_owned(), |b| b.to_string_lossy().into_owned());
    let bytes = std::fs::read(path)
        .map_err(|e| Error(format!("{base}: {}", og_filters::os_error(&e, path))))?;
    let text = String::from_utf8(bytes).map_err(|e| Error(format!("{base}: {e}")))?;
    match pycompat::json::loads(&text) {
        Ok(Value::Object(o)) => Ok(o),
        Ok(_) => Err(Error(format!("{base} is not a JSON object of regions"))),
        Err(e) => Err(Error(format!("{base}: {e}"))),
    }
}

// --- the variant table ---------------------------------------------------------

/// The rows of one Paraphase VCF: every record with a haplotype that
/// carries the ALT allele, with the reads behind each.
fn variant_rows(region: &str, path: &Path, out: &mut Vec<[String; 8]>) -> Result<(), Error> {
    let shown = path.display().to_string();
    let base = path
        .file_name()
        .map_or_else(|| shown.clone(), |b| b.to_string_lossy().into_owned());
    let fail = |e: &io::Error| Error(format!("{base}: {}", og_filters::os_error(e, &shown)));
    let f = File::open(path).map_err(|e| fail(&e))?;
    let mut lines = Lines::new(BufReader::new(f), Utf8::Replace);
    let mut samples: Vec<String> = Vec::new();
    let mut line = String::new();
    while lines.read_line(&mut line).map_err(|e| fail(&e))? {
        if line.starts_with("##") {
            continue;
        }
        let c: Vec<&str> = line
            .strip_suffix('\n')
            .unwrap_or(&line)
            .split('\t')
            .collect();
        if line.starts_with("#CHROM") {
            samples = c.iter().skip(9).map(|s| (*s).to_owned()).collect();
            continue;
        }
        if c.len() < 10 {
            continue;
        }
        let keys: Vec<&str> = c[8].split(':').collect();
        let (mut carriers, mut support) = (Vec::new(), Vec::new());
        for (hap, val) in samples.iter().zip(&c[9..]) {
            let fmt = |key: &str| {
                keys.iter()
                    .zip(val.split(':'))
                    .filter(|(k, _)| **k == key)
                    .map(|(_, v)| v)
                    .last()
            };
            if !matches!(fmt("GT"), Some("1" | "1/1" | "1|1")) {
                continue;
            }
            carriers.push(hap.as_str());
            let ad: Vec<&str> = fmt("AD").unwrap_or("").split(',').collect();
            let alt = ad.get(1).copied().unwrap_or(".");
            support.push(format!("{hap}:{alt}/{}", fmt("DP").unwrap_or(".")));
        }
        if !carriers.is_empty() {
            out.push([
                region.to_owned(),
                c[0].to_owned(),
                c[1].to_owned(),
                c[3].to_owned(),
                c[4].to_owned(),
                c[6].to_owned(),
                carriers.join(","),
                support.join(","),
            ]);
        }
    }
    Ok(())
}

/// The variant table of every `*.vcf` in `dir`, in name order; none
/// without the folder.
pub fn variants_table(dir: Option<&str>, sample: &str) -> Result<Vec<[String; 8]>, Error> {
    let mut out = Vec::new();
    let Some(dir) = dir.filter(|d| !d.is_empty() && Path::new(d).is_dir()) else {
        return Ok(out);
    };
    let mut names = Vec::new();
    // An entry that cannot be read stops the run: a partial listing would
    // publish a variant table missing regions.
    for entry in std::fs::read_dir(dir).map_err(|e| Error(og_filters::os_error(&e, dir)))? {
        let name = entry
            .map_err(|e| Error(og_filters::os_error(&e, dir)))?
            .file_name()
            .to_string_lossy()
            .into_owned();
        if name.ends_with(".vcf") {
            names.push(name);
        }
    }
    names.sort();
    let prefix = format!("{sample}_");
    for name in names {
        let stem = &name[..name.len() - ".vcf".len()];
        let region = stem.strip_prefix(&prefix).unwrap_or(stem);
        variant_rows(region, &Path::new(dir).join(&name), &mut out)?;
    }
    Ok(out)
}

fn tsv(header: &[&str], rows: &[[String; 8]]) -> Vec<u8> {
    let mut s = header.join("\t");
    s.push('\n');
    for r in rows {
        s.push_str(&r.join("\t"));
        s.push('\n');
    }
    s.into_bytes()
}

/// What `og paraphase table` reads and writes.
pub struct TableArgs<'a> {
    pub sample: &'a str,
    pub json: &'a str,
    pub vcf_dir: Option<&'a str>,
    pub out_cn: &'a str,
    pub out_variants: &'a str,
    pub genes: Option<&'a str>,
}

/// `og paraphase table`: writes both tables and returns the log.
pub fn table(a: &TableArgs) -> Result<String, Error> {
    let genes_path = match a.genes {
        Some(g) => g.to_owned(),
        None => default_genes()
            .ok_or_else(|| {
                Error("no assets/paraphase_genes.tsv above og; give --genes".to_owned())
            })?
            .display()
            .to_string(),
    };
    let genes = load_genes(&genes_path)?;
    let rows = cn_table(&load_json(a.json)?, &genes)?;
    let variants = variants_table(a.vcf_dir, a.sample)?;
    if rows.is_empty() {
        return Err(Error(format!("{} has no region", a.json)));
    }
    pgx_rules::write_both(
        (a.out_cn, &tsv(&CN_HEADER, &rows)),
        (a.out_variants, &tsv(&VARIANT_HEADER, &variants)),
    )?;
    let mut log = String::new();
    for r in &rows {
        if genes.iter().any(|(g, _)| *g == r[0]) {
            let _ = writeln!(log, "  {}: {} ({})", r[1], r[2], r[3]);
        }
    }
    let _ = writeln!(
        log,
        "  {} gene rows, {} variant rows",
        rows.len(),
        variants.len()
    );
    Ok(log)
}

// --- the command line ------------------------------------------------------------

const USAGE: &str = "usage: og paraphase check-bam BAM
       og paraphase table --sample SAMPLE --json JSON [--vcf-dir VCF_DIR] --out-cn OUT_CN
                          --out-variants OUT_VARIANTS [--genes GENES]
";

/// `og paraphase ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let refuse = |err: &mut dyn Write, msg: &str| {
        let _ = write!(err, "{USAGE}");
        let _ = writeln!(err, "og paraphase: error: {msg}");
        2
    };
    match args.first().map(String::as_str) {
        Some("-h" | "--help") => {
            let _ = out.write_all(USAGE.as_bytes());
            0
        }
        Some("check-bam") => {
            let rest = &args[1..];
            if rest.iter().any(|a| a == "-h" || a == "--help") {
                let _ = out.write_all(USAGE.as_bytes());
                return 0;
            }
            let [bam] = rest else {
                return refuse(err, "check-bam takes one BAM");
            };
            match check_bam(bam) {
                Ok(Ok(msg)) => pgx_rules::finish(Ok(format!("{msg}\n")), out, err),
                Ok(Err(why)) => {
                    let _ = writeln!(err, "ERROR: {why}");
                    3
                }
                Err(e) => pgx_rules::finish(Err(e), out, err),
            }
        }
        Some("table") => {
            let req = ["--sample", "--json", "--out-cn", "--out-variants"];
            let o = match pgx_rules::parse(
                &args[1..],
                USAGE,
                "og paraphase table",
                &[
                    "--sample",
                    "--json",
                    "--vcf-dir",
                    "--out-cn",
                    "--out-variants",
                    "--genes",
                ],
                &[],
                &req,
                err,
            ) {
                Ok(o) => o,
                Err(code) => return code,
            };
            let get = |n| o.get(n).unwrap_or_default();
            let a = TableArgs {
                sample: get("--sample"),
                json: get("--json"),
                vcf_dir: o.get("--vcf-dir"),
                out_cn: get("--out-cn"),
                out_variants: get("--out-variants"),
                genes: o.get("--genes"),
            };
            pgx_rules::finish(table(&a), out, err)
        }
        Some(other) => refuse(err, &format!("unknown paraphase tool '{other}'")),
        None => refuse(err, "the following arguments are required: paraphase tool"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(tag: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!("og-paraphase-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&p);
        std::fs::create_dir_all(&p).unwrap();
        p
    }

    /// A BGZF-free gzip of a BAM start: magic, `l_text`, the text.
    fn bam(dir: &Path, l_text: i32, text: &str) -> String {
        use flate2::{Compression, write::GzEncoder};
        let mut raw = b"BAM\x01".to_vec();
        raw.extend(l_text.to_le_bytes());
        raw.extend(text.as_bytes());
        raw.extend([0u8; 64]);
        let p = dir.join(format!("x{l_text}.bam"));
        let mut gz = GzEncoder::new(File::create(&p).unwrap(), Compression::fast());
        gz.write_all(&raw).unwrap();
        gz.finish().unwrap();
        p.display().to_string()
    }

    #[test]
    fn only_pacbio_read_groups_pass() {
        let d = tmp("rg");
        let hifi = "@HD\tVN:1.6\n@RG\tID:a\tPL:PACBIO\n@RG\tID:b\tPL:pacbio\n";
        let path = bam(&d, i32::try_from(hifi.len()).unwrap(), hifi);
        assert!(check_bam(&path).unwrap().is_ok());
        // The red-first check of the plan: a non-PacBio read group must be refused.
        let mixed = "@RG\tID:a\tPL:PACBIO\n@RG\tID:b\tPL:ILLUMINA\n@RG\tID:c\n";
        let path = bam(&d, i32::try_from(mixed.len()).unwrap(), mixed);
        let why = check_bam(&path).unwrap().unwrap_err();
        assert!(why.contains("PL ILLUMINA, none, not PACBIO"), "{why}");
        let path = bam(&d, 5, "@HD\tV");
        assert!(
            check_bam(&path)
                .unwrap()
                .unwrap_err()
                .contains("no @RG line")
        );
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn a_corrupt_header_length_is_refused_not_read() {
        let d = tmp("ltext");
        for n in [-1, MAX_HEADER + 1] {
            let e = bam_header(&bam(&d, n, "")).unwrap_err();
            assert!(e.0.contains("not a BAM header"), "{e}");
        }
        let plain = d.join("plain.bam");
        std::fs::write(&plain, "not a bam").unwrap();
        let e = bam_header(&plain.display().to_string()).unwrap_err();
        assert!(e.0.ends_with("is not a BAM file"), "{e}");
        let _ = std::fs::remove_dir_all(d);
    }

    fn obj(json: &str) -> Object {
        match pycompat::json::loads(json).unwrap() {
            Value::Object(o) => o,
            _ => panic!("not an object"),
        }
    }

    fn genes() -> Genes {
        vec![
            (
                "smn1".into(),
                vec![
                    ("SMN1".into(), "smn1_cn".into()),
                    ("SMN2".into(), "smn2_cn".into()),
                ],
            ),
            (
                "hba".into(),
                vec![
                    ("HBA1".into(), "genotype".into()),
                    ("HBA2".into(), "genotype".into()),
                ],
            ),
        ]
    }

    #[test]
    fn numbers_print_as_python_and_objects_are_refused() {
        let data = obj(
            r#"{"GBA": {"total_cn": 4.0, "region_depth": [31.5], "genome_depth": 30,
                        "genes_in_region": ["GBA1", "GBAP1"],
                        "final_haplotypes": {"1": "GBA_gene1hap1", "2": "GBA_gene2hap1"},
                        "two_copy_haplotypes": ["GBA_gene1hap1"]},
                "smn1": {"total_cn": 4, "region_depth": {"median": 60.2},
                         "region_specific_info": {"smn1_cn": 2, "smn2_cn": null}}}"#,
        );
        let rows = cn_table(&data, &genes()).unwrap();
        assert_eq!(rows[0][..5], ["smn1", "SMN1", "2", "called", "."]);
        assert_eq!(rows[1][2..4], [".", "no call"]);
        assert_eq!(rows[2][1..3], ["GBA1", "2"]);
        assert_eq!(rows[2][5..], ["4", "31.5", "30"]);
        let bad = obj(r#"{"x": {"total_cn": {"a": 1}}}"#);
        let e = cn_table(&bad, &genes()).unwrap_err();
        assert_eq!(e.0, "x.total_cn is not a number");
        let bad = obj(r#"{"smn1": {"region_specific_info": {"smn1_cn": true}}}"#);
        assert_eq!(
            cn_table(&bad, &genes()).unwrap_err().0,
            "smn1.smn1_cn is not a number"
        );
    }

    #[test]
    fn hba_counts_haplotypes_and_names_svs() {
        let data = obj(r#"{"hba": {"total_cn": 3,
                        "final_haplotypes": ["hba_hba2hap1", "hba_hba1hap1", "hba_hba1hap2"],
                        "region_specific_info": {"genotype": "-a/aa", "sv_called": {"h1": "3p7del,x"}}}}"#);
        let rows = cn_table(&data, &genes()).unwrap();
        assert_eq!(rows[0][1..4], ["HBA1", "2", "called"]);
        assert_eq!(
            rows[0][4],
            "genotype -a/aa (3 alpha-globin genes); 3.7 kb deletion (-a3.7): HBA2 and HBA1 fused into one gene, one alpha-globin gene lost (h1)"
        );
        assert_eq!(rows[1][1..3], ["HBA2", "1"]);
    }

    #[test]
    fn both_tables_or_neither() {
        let d = tmp("both");
        std::fs::write(d.join("p.json"), r#"{"smn1": {"total_cn": 4}}"#).unwrap();
        std::fs::write(d.join("g.tsv"), "region\tgene\tfield\treading\tsource\n").unwrap();
        let s = |p: &str| d.join(p).display().to_string();
        let (json, g, cn) = (s("p.json"), s("g.tsv"), s("cn.tsv"));
        let a = TableArgs {
            sample: "S",
            json: &json,
            vcf_dir: None,
            out_cn: &cn,
            out_variants: &s("missing/v.tsv"),
            genes: Some(&g),
        };
        assert!(table(&a).is_err());
        assert!(
            !d.join("cn.tsv").exists(),
            "the copy-number table stands alone"
        );
        let _ = std::fs::remove_dir_all(d);
    }
}
