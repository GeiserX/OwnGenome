//! Step 32's summary of pypgx's results (`docs/32-pypgx.md`).
//!
//! ```text
//! og pgx pypgx-summary --outbase DIR --bam-genes "G ..." --vcf-genes "G ..."
//!                      --vcf VCF --out SUMMARY.tsv [--not-callable GENE ...]
//!                      [--gq-floor 20] [--pypgx-data DIR]
//! ```
//!
//! Per gene: the diplotype, phenotype and copy-number call pypgx wrote in
//! `DIR/<gene>/results.zip`, and for each star allele called the variants
//! that define it (pypgx's core variants on GRCh38), looked up in the VCF
//! pypgx read: how many it holds, the lowest GQ among them, the weakest one
//! with its AD, their FILTER, and `low-confidence: check the reads` when one
//! has a GQ below the floor. A gene in `--not-callable` (GSTT1 on a
//! reference without ALT contigs) is not callable on this reference, not
//! FAILED.
//!
//! It reads pypgx's own files instead of running pypgx's Python: the
//! results archive is a zip whose `data.tsv` is the sample table `pypgx
//! print-data` prints, and the star-allele definitions and phenotype rules
//! are pypgx's packaged `allele-table.csv` and `gene-table.csv`, found
//! beside the `pypgx` on `PATH` (`--pypgx-data` names the folder). The VCF
//! is read through its tabix index. pandas reads those tables, so a cell
//! pandas reads as missing (`N/A`, `NA`, empty, ...) is missing here too.
//!
//! This is `bin/pgx_outside_calls.py pypgx-summary` rewritten. On purpose:
//! a VCF without an index, or pypgx tables that cannot be found, give one
//! WARNING and `not looked up`, where the Python stopped with a traceback;
//! the summary is written whole or not at all.

use std::collections::HashSet;
use std::fmt::Write as _;
use std::io::{self, BufRead, Read, Write};
use std::path::{Path, PathBuf};

use og_io::Utf8;
use pycompat::csv::{Quoting, Reader};

use crate::pgx_rules::{self, Error, NO_CALL, TextLines, is_call, strip};

/// The summary's columns.
pub const SUMMARY_HEADER: [&str; 11] = [
    "Gene",
    "Diplotype",
    "Phenotype",
    "CNV_call",
    "Source",
    "Status",
    "Defining_variants",
    "Min_GQ",
    "Weakest_variant",
    "Defining_filter",
    "Quality",
];
/// The status of a gene the reference cannot hold.
pub const NOT_CALLABLE: &str = "not callable on this reference";
/// The quality flag of a call resting on a variant below the GQ floor.
pub const LOW_GQ: &str = "low-confidence: check the reads";
/// The quality flag of a call resting on a variant the VCF gives no GQ.
pub const NO_GQ: &str = "GQ not in the VCF: check the reads";
/// The GQ below which a defining variant flags the call.
pub const DEFAULT_GQ_FLOOR: i64 = 20;

/// The cells pandas' `read_csv` reads as missing (pandas 2.3, the pypgx
/// 0.26.0 image).
const PANDAS_NA: [&str; 19] = [
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN",
    "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null",
];

fn pandas_na(s: &str) -> bool {
    PANDAS_NA.contains(&s)
}

// --- star alleles and their evidence -----------------------------------------

/// The star alleles of a pypgx diplotype, each once: `*2x2/*4` gives `*2`,
/// `*4` (a copy count is not an allele), `*36+*10/*4` gives `*36`, `*10`,
/// `*4`; a no-call gives none.
pub fn star_alleles(diplotype: &str) -> Vec<String> {
    let d = strip(diplotype);
    if NO_CALL.contains(&pycompat::text::lower(d).as_str()) || !d.contains('/') {
        return Vec::new();
    }
    let mut out: Vec<String> = Vec::new();
    for half in d.split('/') {
        for a in half.split('+') {
            let a = without_copy_count(strip(a));
            if !a.is_empty() && !out.iter().any(|x| x == a) {
                out.push(a.to_owned());
            }
        }
    }
    out
}

/// `re.sub(r"x(\d+|N)$", "", a)`: a trailing copy count cut.
fn without_copy_count(a: &str) -> &str {
    if let Some(head) = a.strip_suffix("xN") {
        return head;
    }
    let digits = a.len() - a.trim_end_matches(|c: char| c.is_ascii_digit()).len();
    if digits > 0 && a[..a.len() - digits].ends_with('x') {
        return &a[..a.len() - digits - 1];
    }
    a
}

/// What the VCF says at one defining variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    /// The FILTER column (`PASS`, `.`, or the filters).
    pub filter: String,
    /// The first sample's GQ, when it has one.
    pub gq: Option<i64>,
    /// The first sample's AD, `.` without one.
    pub ad: String,
}

/// `chr22:42128945 C>T` for pypgx's `22-42128945-C-T`.
fn variant_name(v: &str) -> String {
    let v = v.replacen("chr", "", 1);
    let p: Vec<&str> = v.splitn(4, '-').collect();
    match p.as_slice() {
        [c, pos, r, a] => format!("chr{c}:{pos} {r}>{a}"),
        _ => format!("chr{v}"),
    }
}

/// The five evidence cells (`Defining_variants Min_GQ Weakest_variant
/// Defining_filter Quality`) of one call. `variants_of(allele)`: the
/// variants pypgx defines it by, empty for one defined by none (`*1`, a
/// deletion), `None` for one pypgx's table lacks. `lookup(variant)`: the
/// VCF's record there, `None` when it has none.
pub fn defining_evidence(
    diplotype: &str,
    mut variants_of: impl FnMut(&str) -> Option<Vec<String>>,
    mut lookup: impl FnMut(&str) -> Result<Option<Site>, Error>,
    floor: i64,
) -> Result<[String; 5], Error> {
    let na = || "N/A".to_owned();
    let alleles = star_alleles(diplotype);
    if alleles.is_empty() {
        return Ok([na(), na(), na(), na(), na()]);
    }
    let (mut variants, mut unknown): (Vec<String>, Vec<String>) = (Vec::new(), Vec::new());
    for a in &alleles {
        match variants_of(a) {
            None => unknown.push(a.clone()),
            Some(vs) => {
                for v in vs {
                    if !variants.contains(&v) {
                        variants.push(v);
                    }
                }
            }
        }
    }
    let not_in_table =
        (!unknown.is_empty()).then(|| format!("not in pypgx's table: {}", unknown.join(", ")));
    if variants.is_empty() {
        return Ok(match not_in_table {
            Some(t) => ["unknown".to_owned(), na(), na(), t, na()],
            None => ["none".to_owned(), na(), na(), "none".to_owned(), na()],
        });
    }
    let (mut found, mut missing, mut not_pass) = (Vec::new(), Vec::new(), Vec::new());
    for v in &variants {
        let name = variant_name(v);
        let Some(site) = lookup(v)? else {
            missing.push(name);
            continue;
        };
        if site.filter != "PASS" {
            let f = if site.filter == "." {
                "no filter"
            } else {
                site.filter.as_str()
            };
            not_pass.push(format!("{name} {f}"));
        }
        found.push((site.gq, name, site.ad));
    }
    let with_gq: Vec<(i64, &String, &String)> = found
        .iter()
        .filter_map(|(gq, name, ad)| gq.map(|g| (g, name, ad)))
        .collect();
    let weakest = with_gq.iter().min().copied();
    let filt = if missing.is_empty() && not_pass.is_empty() && unknown.is_empty() {
        "PASS".to_owned()
    } else {
        let mut parts = Vec::new();
        if !not_pass.is_empty() {
            parts.push(format!("not PASS: {}", not_pass.join(", ")));
        }
        if !missing.is_empty() {
            parts.push(format!("not in the VCF: {}", missing.join(", ")));
        }
        parts.extend(not_in_table);
        parts.join("; ")
    };
    let quality = if weakest.is_some_and(|w| w.0 < floor) {
        LOW_GQ
    } else if with_gq.len() < found.len() {
        NO_GQ
    } else {
        "ok"
    };
    Ok([
        format!("{}/{} in the VCF", found.len(), variants.len()),
        weakest.map_or_else(na, |w| w.0.to_string()),
        weakest.map_or_else(na, |(gq, name, ad)| format!("{name} GQ {gq} AD {ad}")),
        filt,
        quality.to_owned(),
    ])
}

// --- pypgx's tables ------------------------------------------------------------

/// pypgx's gene and allele tables: what `has_phenotype` and `list_variants`
/// read.
pub struct Tables {
    /// Per gene: a target gene, and one with a phenotype rule.
    genes: Vec<(String, bool, bool)>,
    /// `(gene, star allele, GRCh38 core variants)`, in the table's order.
    alleles: Vec<(String, String, Option<String>)>,
}

/// The rows of a comma-separated table with a header, as `pd.read_csv`
/// reads them for the columns asked for.
fn csv_columns(path: &Path, cols: &[&str]) -> Result<Vec<Vec<String>>, Error> {
    let shown = path.display().to_string();
    let lines = TextLines::open(&shown, Utf8::Strict)?;
    let mut rows = Reader::new(lines, ',', Quoting::Minimal);
    let header = match rows.next() {
        Some(h) => h.map_err(|e| Error(format!("{shown}: {e}")))?,
        None => return Err(Error(format!("{shown} is empty"))),
    };
    let at: Vec<usize> = cols
        .iter()
        .map(|c| {
            header
                .iter()
                .position(|h| h == c)
                .ok_or_else(|| Error(format!("{shown} has no {c} column")))
        })
        .collect::<Result<_, _>>()?;
    let mut out = Vec::new();
    for row in rows {
        let row = row.map_err(|e| Error(format!("{shown}: {e}")))?;
        if row.iter().all(String::is_empty) {
            continue;
        }
        out.push(
            at.iter()
                .map(|&i| row.get(i).cloned().unwrap_or_default())
                .collect(),
        );
    }
    Ok(out)
}

impl Tables {
    /// Reads `gene-table.csv` and `allele-table.csv` from `dir`.
    pub fn read(dir: &Path) -> Result<Tables, Error> {
        let genes = csv_columns(
            &dir.join("gene-table.csv"),
            &["Gene", "Target", "PhenotypeMethod"],
        )?
        .into_iter()
        .map(|r| {
            let [gene, target, method] = <[String; 3]>::try_from(r).expect("three columns");
            (
                gene,
                target.eq_ignore_ascii_case("true"),
                !pandas_na(&method),
            )
        })
        .collect();
        let alleles = csv_columns(
            &dir.join("allele-table.csv"),
            &["Gene", "StarAllele", "GRCh38Core"],
        )?
        .into_iter()
        .map(|r| {
            let [gene, allele, core] = <[String; 3]>::try_from(r).expect("three columns");
            (gene, allele, (!pandas_na(&core)).then_some(core))
        })
        .collect();
        Ok(Tables { genes, alleles })
    }

    fn target(&self, gene: &str) -> bool {
        self.genes.iter().any(|(g, t, _)| g == gene && *t)
    }

    /// `pypgx.has_phenotype(gene)`; `Err` names a gene that is not one of
    /// pypgx's target genes, where pypgx raises.
    pub fn has_phenotype(&self, gene: &str) -> Result<bool, String> {
        if !self.target(gene) {
            return Err(format!("{gene} is not one of pypgx's target genes"));
        }
        Ok(self.genes.iter().any(|(g, _, p)| g == gene && *p))
    }

    /// `pypgx.list_variants(gene, alleles=[allele], mode="core",
    /// assembly="GRCh38")`, or `None` where pypgx raises (not a target gene,
    /// an allele its table lacks).
    pub fn core_variants(&self, gene: &str, allele: &str) -> Option<Vec<String>> {
        if !self.target(gene) {
            return None;
        }
        let (_, _, core) = self
            .alleles
            .iter()
            .find(|(g, a, _)| g == gene && a == allele)?;
        let mut vs: Vec<String> = core
            .as_deref()
            .map(|c| c.split(',').map(str::to_owned).collect())
            .unwrap_or_default();
        vs.sort_by_key(|v| sort_key(v));
        vs.dedup();
        Some(vs)
    }
}

/// fuc's `sort_variants` order: chromosome (1 to 22, X, Y, M), position,
/// REF, ALT.
fn sort_key(v: &str) -> (usize, String, u64, String, String) {
    let p: Vec<&str> = v.splitn(4, '-').collect();
    let chrom = p.first().copied().unwrap_or_default();
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let rank = match bare {
        "X" => 23,
        "Y" => 24,
        "M" | "MT" => 25,
        n => n
            .parse::<usize>()
            .ok()
            .filter(|n| (1..=22).contains(n))
            .unwrap_or(26),
    };
    let pos = p.get(1).and_then(|x| x.parse().ok()).unwrap_or(0);
    let s = |i: usize| p.get(i).copied().unwrap_or_default().to_owned();
    (rank, chrom.to_owned(), pos, s(2), s(3))
}

/// pypgx's packaged tables beside the `pypgx` on `PATH`:
/// `<prefix>/bin/pypgx` and `<prefix>/lib/pythonX.Y/site-packages/pypgx/api/data`,
/// for the Python its `#!` line names (an image may hold a `python3.1`
/// link beside `python3.12`), else the newest Python that has them.
fn find_data() -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    let exe = std::env::split_paths(&path)
        .map(|d| d.join("pypgx"))
        .find(|p| p.is_file())?;
    let lib = exe.parent()?.parent()?.join("lib");
    let data = |python: &str| lib.join(python).join("site-packages/pypgx/api/data");
    let has = |p: &Path| p.join("allele-table.csv").is_file();
    let mut first = String::new();
    if let Ok(f) = std::fs::File::open(&exe) {
        let _ = io::BufReader::new(f).take(256).read_line(&mut first);
    }
    let named = first
        .strip_prefix("#!")
        .and_then(|l| l.split_whitespace().next())
        .and_then(|i| i.rsplit('/').next())
        .map(data)
        .filter(|p| has(p));
    named.or_else(|| {
        let version = |name: &str| -> Option<(u32, u32)> {
            let (a, b) = name.strip_prefix("python")?.split_once('.')?;
            Some((a.parse().ok()?, b.parse().ok()?))
        };
        std::fs::read_dir(&lib)
            .ok()?
            .filter_map(Result::ok)
            .filter_map(|e| {
                let name = e.file_name().to_string_lossy().into_owned();
                Some((version(&name)?, data(&name)))
            })
            .filter(|(_, p)| has(p))
            .max_by_key(|(v, _)| *v)
            .map(|(_, p)| p)
    })
}

// --- results.zip ---------------------------------------------------------------

/// The first sample's cells of a pypgx results archive, by column name:
/// what `pypgx print-data` prints, a cell pandas reads as missing printed
/// empty.
fn read_results(path: &Path) -> Result<Vec<(String, String)>, String> {
    let mut zip = og_setup::zip::Archive::open(path).map_err(|e| e.to_string())?;
    let parent = zip
        .entries()
        .first()
        .map(|e| e.name.split('/').next().unwrap_or_default().to_owned())
        .ok_or("the archive is empty")?;
    let mut member = |name: &str| -> Result<String, String> {
        let full = format!("{parent}/{name}");
        let e = zip
            .find(&full)
            .cloned()
            .ok_or_else(|| format!("no {full} in the archive"))?;
        let mut body = Vec::new();
        zip.copy_to(&e, &mut body).map_err(|e| e.to_string())?;
        String::from_utf8(body).map_err(|e| format!("{full}: {e}"))
    };
    let meta = member("metadata.txt")?;
    let kind = meta
        .lines()
        .filter_map(|l| strip(l).split_once('='))
        .rev()
        .find(|(k, _)| *k == "SemanticType")
        .map(|(_, v)| v.to_owned())
        .unwrap_or_default();
    if !kind.contains("SampleTable") {
        return Err(format!("the archive holds {kind}, not a sample table"));
    }
    let data = member("data.tsv")?;
    let lines = data
        .split_inclusive('\n')
        .map(|l| Ok::<_, io::Error>(l.to_owned()));
    let mut rows = Reader::new(lines, '\t', Quoting::Minimal);
    let mut next = || -> Result<Vec<String>, String> {
        rows.next()
            .transpose()
            .map_err(|e| e.to_string())
            .map(Option::unwrap_or_default)
    };
    let (header, first) = (next()?, next()?);
    Ok(header
        .into_iter()
        .zip(first)
        .map(|(h, v)| (h, if pandas_na(&v) { String::new() } else { v }))
        .collect())
}

// --- the VCF ---------------------------------------------------------------------

/// The VCF pypgx read, looked up by variant through its tabix index.
pub struct Vcf {
    reader: og_io::tabix::Reader,
    contigs: HashSet<String>,
}

impl Vcf {
    /// Opens a bgzipped VCF with its `.tbi` or `.csi`, and reads the
    /// contigs its index and its header name (a header may have no
    /// `##contig` lines).
    pub fn open(path: &str) -> io::Result<Vcf> {
        let reader = og_io::tabix::Reader::open(path)?;
        let mut contigs: HashSet<String> = reader.index().names().iter().cloned().collect();
        let mut head = og_io::open(path)?;
        let mut line = String::new();
        loop {
            line.clear();
            if head.read_line(&mut line)? == 0 || !line.starts_with('#') {
                break;
            }
            if let Some(rest) = line.strip_prefix("##contig=<") {
                let id = rest
                    .trim_end()
                    .trim_end_matches('>')
                    .split(',')
                    .find_map(|kv| kv.strip_prefix("ID="));
                if let Some(id) = id {
                    contigs.insert(id.to_owned());
                }
            }
        }
        Ok(Vcf { reader, contigs })
    }

    /// The first record at pypgx's `22-42128945-C-T` with that REF and that
    /// ALT among its ALTs: its FILTER and the first sample's GQ and AD.
    pub fn lookup(&mut self, v: &str) -> io::Result<Option<Site>> {
        let p: Vec<&str> = v.splitn(4, '-').collect();
        let [chrom, pos, r, a] = p.as_slice() else {
            return Ok(None);
        };
        let chrom = chrom.replacen("chr", "", 1);
        let with = format!("chr{chrom}");
        let c = if self.contigs.contains(&with) {
            with
        } else if self.contigs.contains(&chrom) {
            chrom
        } else {
            return Ok(None);
        };
        let Ok(pos) = pos.parse::<u64>() else {
            return Ok(None);
        };
        for line in self.reader.query(&c, pos.saturating_sub(1), pos) {
            let line = line?;
            let line = String::from_utf8_lossy(&line);
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8
                || f[1].parse::<u64>().ok() != Some(pos)
                || f[3] != *r
                || !f[4].split(',').any(|x| x == *a)
            {
                continue;
            }
            let filter = if f[6].is_empty() { "." } else { f[6] };
            let sample = |key: &str| -> Option<&str> {
                let at = f.get(8)?.split(':').position(|k| k == key)?;
                f.get(9)?.split(':').nth(at)
            };
            let gq = sample("GQ")
                .and_then(|g| pycompat::num::parse_float(g).ok())
                .filter(|g| g.is_finite())
                .map(|g| g.trunc() as i64);
            let ad = match sample("AD") {
                None | Some("" | ".") => ".".to_owned(),
                Some(ad) => ad
                    .split(',')
                    .map(|x| match x.parse::<i64>() {
                        Ok(n) => n.to_string(),
                        Err(_) => x.to_owned(),
                    })
                    .collect::<Vec<_>>()
                    .join(","),
            };
            return Ok(Some(Site {
                filter: filter.to_owned(),
                gq,
                ad,
            }));
        }
        Ok(None)
    }
}

// --- the summary ---------------------------------------------------------------

/// What `og pgx pypgx-summary` reads.
#[derive(Debug, Clone)]
pub struct Args {
    pub outbase: String,
    pub bam_genes: Vec<String>,
    pub vcf_genes: Vec<String>,
    pub vcf: String,
    pub out: String,
    pub not_callable: Vec<String>,
    pub gq_floor: i64,
    pub pypgx_data: Option<String>,
}

/// Writes the summary; the stdout text, with its WARNING lines on `err`.
pub fn summary(a: &Args, err: &mut dyn Write) -> Result<String, Error> {
    let genes: Vec<&String> = a.bam_genes.iter().chain(&a.vcf_genes).collect();
    let data = a.pypgx_data.as_ref().map(PathBuf::from).or_else(find_data);
    let tables = match data.as_deref().map(Tables::read) {
        Some(Ok(t)) => Some(t),
        Some(Err(e)) => {
            let _ = writeln!(
                err,
                "WARNING: could not read pypgx gene table ({e}); phenotypes left as pypgx wrote them and the defining variants are not looked up"
            );
            None
        }
        None => {
            let _ = writeln!(
                err,
                "WARNING: no pypgx tables beside the pypgx on PATH (give --pypgx-data); phenotypes left as pypgx wrote them and the defining variants are not looked up"
            );
            None
        }
    };
    // pypgx writes Indeterminate for a gene it has no phenotype rule for:
    // say that instead, so the label is not read as an undecided call.
    let mut no_rule: HashSet<&str> = HashSet::new();
    if let Some(t) = &tables {
        for g in &genes {
            match t.has_phenotype(g) {
                Ok(true) => {}
                Ok(false) => {
                    no_rule.insert(g.as_str());
                }
                Err(e) => {
                    let _ = writeln!(
                        err,
                        "WARNING: could not read pypgx gene table ({e}); phenotypes left as pypgx wrote them"
                    );
                    no_rule.clear();
                    break;
                }
            }
        }
    }
    let mut vcf = match &tables {
        None => None,
        Some(_) => match Vcf::open(&a.vcf) {
            Ok(v) => Some(v),
            Err(e) => {
                let _ = writeln!(
                    err,
                    "WARNING: cannot read {} ({e}); the defining variants are not looked up",
                    a.vcf
                );
                None
            }
        },
    };
    let mut rows: Vec<Vec<String>> = Vec::new();
    let na = || "N/A".to_owned();
    for gene in &genes {
        let gene = gene.as_str();
        let source = if a.bam_genes.iter().any(|g| g == gene) {
            "BAM"
        } else {
            "VCF"
        };
        let head = |d: &str, p: &str, c: &str, status: &str| {
            [gene, d, p, c, source, status].map(str::to_owned).to_vec()
        };
        if a.not_callable.iter().any(|g| g == gene) {
            let mut r = head("N/A", "N/A", "N/A", NOT_CALLABLE);
            r.extend(std::iter::repeat_with(na).take(5));
            rows.push(r);
            continue;
        }
        let results = Path::new(&a.outbase).join(gene).join("results.zip");
        if !results.is_file() {
            let mut r = head("FAILED", "N/A", "N/A", "failed");
            r.extend(std::iter::repeat_with(na).take(5));
            rows.push(r);
            continue;
        }
        let (mut diplotype, mut phenotype, mut cnv) = (na(), na(), na());
        match read_results(&results) {
            Err(e) => {
                let _ = writeln!(err, "WARNING: cannot read {gene}'s results.zip: {e}");
            }
            Ok(cells) => {
                let col = |k: &str| cells.iter().find(|(h, _)| h == k).map(|(_, v)| v.clone());
                diplotype = col("Genotype").filter(|v| !v.is_empty()).unwrap_or_else(na);
                phenotype = col("Phenotype")
                    .filter(|v| !v.is_empty())
                    .unwrap_or_else(na);
                if no_rule.contains(gene) {
                    phenotype = "no pypgx phenotype rule".to_owned();
                }
                // Copy number is called from read depth: meaningful only for
                // the BAM-based genes, pypgx's own value, never a guess.
                if source == "BAM" {
                    cnv = col("CNV").unwrap_or_else(|| {
                        let _ = writeln!(err, "WARNING: {gene} results have no CNV column");
                        na()
                    });
                }
            }
        }
        let evidence = match (&tables, &mut vcf) {
            (Some(t), Some(v)) => defining_evidence(
                &diplotype,
                |allele| t.core_variants(gene, allele),
                |variant| v.lookup(variant).map_err(|e| Error::io(&e, &a.vcf)),
                a.gq_floor,
            )?
            .to_vec(),
            _ => {
                let mut e = vec!["not looked up".to_owned()];
                e.extend(std::iter::repeat_with(na).take(4));
                e
            }
        };
        let status = if is_call(&diplotype) {
            "called"
        } else {
            "no call"
        };
        let mut r = head(&diplotype, &phenotype, &cnv, status);
        r.extend(evidence);
        rows.push(r);
    }
    let mut w = pycompat::csv::Writer::new(Vec::new(), '\t');
    let _ = w.write_row(SUMMARY_HEADER);
    for r in &rows {
        let _ = w.write_row(r);
    }
    let fail = |e: io::Error| Error::io(&e, &a.out);
    let mut f = og_io::AtomicFile::create(&a.out).map_err(fail)?;
    f.write_all(&w.into_inner()).map_err(fail)?;
    f.commit().map_err(fail)?;
    let called = rows.iter().filter(|r| r[5] == "called").count();
    let mut log = format!("Genes called: {called}/{}\n", rows.len());
    let low: Vec<String> = rows
        .iter()
        .filter(|r| r[10] == LOW_GQ)
        .map(|r| format!("{} {} ({})", r[0], r[1], r[8]))
        .collect();
    if !low.is_empty() {
        let _ = writeln!(
            log,
            "WARNING: star alleles resting on a variant below GQ {}, check the reads: {}",
            a.gq_floor,
            low.join("; ")
        );
    }
    Ok(log)
}

const USAGE: &str =
    "usage: og pgx pypgx-summary --outbase OUTBASE --bam-genes BAM_GENES --vcf-genes VCF_GENES
                            --vcf VCF --out OUT [--not-callable NOT_CALLABLE] [--gq-floor GQ_FLOOR]
                            [--pypgx-data DIR]
";

/// `og pgx pypgx-summary ARGS`.
pub fn cli(args: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let req = ["--outbase", "--bam-genes", "--vcf-genes", "--vcf", "--out"];
    let o = match pgx_rules::parse(
        args,
        USAGE,
        "og pgx pypgx-summary",
        &[
            "--outbase",
            "--bam-genes",
            "--vcf-genes",
            "--vcf",
            "--out",
            "--not-callable",
            "--gq-floor",
            "--pypgx-data",
        ],
        &[],
        &req,
        err,
    ) {
        Ok(o) => o,
        Err(code) => return code,
    };
    let gq_floor = match o.get("--gq-floor") {
        None => DEFAULT_GQ_FLOOR,
        Some(v) => match pycompat::num::parse_int(v)
            .ok()
            .and_then(|n| i64::try_from(n).ok())
        {
            Some(n) => n,
            None => {
                let _ = write!(err, "{USAGE}");
                let _ = writeln!(
                    err,
                    "og pgx pypgx-summary: error: argument --gq-floor: invalid int value: {}",
                    og_filters::py_repr(v)
                );
                return 2;
            }
        },
    };
    let words = |n: &str| -> Vec<String> {
        pycompat::text::split(o.get(n).unwrap_or_default(), None, None)
            .into_iter()
            .map(str::to_owned)
            .collect()
    };
    let a = Args {
        outbase: o.get("--outbase").unwrap_or_default().to_owned(),
        bam_genes: words("--bam-genes"),
        vcf_genes: words("--vcf-genes"),
        vcf: o.get("--vcf").unwrap_or_default().to_owned(),
        out: o.get("--out").unwrap_or_default().to_owned(),
        not_callable: o.all("--not-callable"),
        gq_floor,
        pypgx_data: o.get("--pypgx-data").map(str::to_owned),
    };
    let r = summary(&a, err);
    pgx_rules::finish(r, out, err)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn star_alleles_of_a_diplotype() {
        assert_eq!(star_alleles("*2x2/*4"), ["*2", "*4"]);
        assert_eq!(star_alleles("*36+*10/*10"), ["*36", "*10"]);
        assert_eq!(star_alleles("*1xN/*4"), ["*1", "*4"]);
        assert!(star_alleles("Indeterminate").is_empty());
        assert!(star_alleles("*1").is_empty());
    }

    fn defs(a: &str) -> Option<Vec<String>> {
        let v = |xs: &[&str]| Some(xs.iter().map(|x| (*x).to_owned()).collect());
        match a {
            "*1" => v(&[]),
            "*2" => v(&["10-94775367-A-G", "10-94781859-G-A"]),
            "*17" => v(&["10-94761900-C-T"]),
            _ => None,
        }
    }

    fn vcf(
        gq: Option<i64>,
        filt: &str,
        drop: &str,
    ) -> impl FnMut(&str) -> Result<Option<Site>, Error> {
        let filt = filt.to_owned();
        let drop = drop.to_owned();
        move |v: &str| {
            let site = |f: &str, gq, ad: &str| Site {
                filter: f.to_owned(),
                gq,
                ad: ad.to_owned(),
            };
            Ok(match v {
                _ if v == drop => None,
                "10-94775367-A-G" => Some(site("PASS", Some(48), "20,19")),
                "10-94781859-G-A" => Some(site(&filt, gq, "26,3")),
                "10-94761900-C-T" => Some(site("PASS", Some(60), "0,31")),
                _ => None,
            })
        }
    }

    fn ev(d: &str, gq: Option<i64>, filt: &str, drop: &str) -> [String; 5] {
        defining_evidence(d, defs, vcf(gq, filt, drop), DEFAULT_GQ_FLOOR).unwrap()
    }

    #[test]
    fn a_low_gq_defining_variant_flags_the_call() {
        let e = ev("*2/*17", Some(12), "PASS", "");
        assert_eq!(e[4], LOW_GQ);
        assert_eq!(e[1..3], ["12", "chr10:94781859 G>A GQ 12 AD 26,3"]);
        assert_eq!((e[0].as_str(), e[3].as_str()), ("3/3 in the VCF", "PASS"));
        assert_eq!(ev("*2/*17", Some(50), "PASS", "")[4], "ok");
        assert_eq!(ev("*2/*17", Some(19), "PASS", "")[4], LOW_GQ);
        assert_eq!(ev("*2/*17", Some(20), "PASS", "")[4], "ok");
        assert_eq!(ev("*2/*17", None, "PASS", "")[4], NO_GQ);
    }

    #[test]
    fn missing_unfiltered_and_unknown_variants_are_named() {
        assert_eq!(
            ev("*2/*17", Some(50), "PASS", "10-94775367-A-G")[3],
            "not in the VCF: chr10:94775367 A>G"
        );
        assert_eq!(
            ev("*2/*17", Some(50), "RefCall", "")[3],
            "not PASS: chr10:94781859 G>A RefCall"
        );
        assert_eq!(
            ev("*2/*17", Some(50), ".", "")[3],
            "not PASS: chr10:94781859 G>A no filter"
        );
        assert_eq!(
            ev("*1/*1", Some(50), "PASS", ""),
            ["none", "N/A", "N/A", "none", "N/A"]
        );
        assert_eq!(
            ev("*2/*99", Some(50), "PASS", "")[3],
            "not in pypgx's table: *99"
        );
        assert_eq!(
            ev("*98/*99", Some(50), "PASS", ""),
            [
                "unknown",
                "N/A",
                "N/A",
                "not in pypgx's table: *98, *99",
                "N/A"
            ]
        );
        assert_eq!(
            ev("Indeterminate", Some(50), "PASS", ""),
            ["N/A", "N/A", "N/A", "N/A", "N/A"]
        );
    }

    #[test]
    fn core_variants_sorted_and_unknowns_refused() {
        let t = Tables {
            genes: vec![
                ("CYP2C19".into(), true, true),
                ("GSTM1".into(), true, false),
                ("CYP2D7".into(), false, false),
            ],
            alleles: vec![
                (
                    "CYP2C19".into(),
                    "*2".into(),
                    Some("10-94781859-G-A,10-94775367-A-G".into()),
                ),
                ("CYP2C19".into(), "*1".into(), None),
            ],
        };
        assert_eq!(
            t.core_variants("CYP2C19", "*2").unwrap(),
            ["10-94775367-A-G", "10-94781859-G-A"]
        );
        assert_eq!(
            t.core_variants("CYP2C19", "*1").unwrap(),
            Vec::<String>::new()
        );
        assert_eq!(t.core_variants("CYP2C19", "*99"), None);
        assert_eq!(t.core_variants("CYP2D7", "*1"), None);
        assert_eq!(t.has_phenotype("GSTM1"), Ok(false));
        assert!(t.has_phenotype("CYP2D7").is_err());
    }

    #[test]
    fn a_vcf_without_contig_lines_is_still_looked_up() {
        let path = format!(
            "{}/../../../tests/fixtures/pgx/summary/no-contig.vcf.gz",
            env!("CARGO_MANIFEST_DIR")
        );
        let mut v = Vcf::open(&path).unwrap();
        let site = v.lookup("10-94781859-G-A").unwrap().unwrap();
        assert_eq!(site.gq, Some(12));
        assert_eq!(site.filter, "PASS");
    }
}
