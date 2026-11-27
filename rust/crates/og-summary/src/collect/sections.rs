//! The readers of group a: depth, gene callability, sex and karyotype,
//! sample identity, mosaic sex-chromosome loss, read count and variant
//! counts. Each returns the file it read and its typed data, `None` when
//! the step wrote nothing, or the note of an unreadable section.

use super::model::{Manifest, RUN_ALL_STEPS, RunStatus, Section, State, is_stale, opt_str};
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
    /// A section of group c, kept as the ordered JSON it writes.
    Json(Value),
    ClinVar(super::clinvar::ClinVar),
    Carrier(super::carrier::Carrier),
    AcmgSf(super::acmg_card::AcmgCard),
    PharmCat(super::pgx::PharmCat),
    Cpic(super::pgx::Cpic),
    Pypgx(super::pgx::Pypgx),
    Cyrius(super::pgx::Cyrius),
    Cpsr(super::cpsr::Cpsr),
    Clinical(super::clinical::Clinical),
    Slivar(super::clinical::Slivar),
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
            Data::Json(v) => v.clone(),
            Data::ClinVar(d) => d.to_json(),
            Data::Carrier(d) => d.to_json(),
            Data::AcmgSf(d) => d.to_json(),
            Data::PharmCat(d) => d.to_json(),
            Data::Cpic(d) => d.to_json(),
            Data::Pypgx(d) => d.to_json(),
            Data::Cyrius(d) => d.to_json(),
            Data::Cpsr(d) => d.to_json(),
            Data::Clinical(d) => d.to_json(),
            Data::Slivar(d) => d.to_json(),
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
    /// The `mean` of the `total` row, the last one when there are several;
    /// for an exome, of the `total_region` row.
    pub mean_depth: Option<f64>,
    /// The run is an exome (by the rule `read` describes), so the mean is
    /// over the capture regions. summary.json carries `on_target: true`
    /// only then, so a whole genome's section keeps its one key; og facts
    /// reads it to write the assay type `wes`.
    pub on_target: bool,
}

impl Coverage {
    /// Reads mosdepth's summary. A whole genome's mean depth is the `total`
    /// row. An exome's is the `total_region` row, the mean over the capture
    /// regions mosdepth was given: `total` spreads the reads over the whole
    /// genome, a few X for a 100x exome, and every reader of the mean (the
    /// QC lines, the live headline, the genome facts' depth floor) would
    /// read the exome as too thin. An exome summary without `total_region`
    /// has no mean depth.
    ///
    /// The run is an exome when `exome` says so (the run manifest's
    /// `run model_type` is `WES`), or when the summary shows a capture: a
    /// `total_region` under half as long as `total`. 500 bp windows tile
    /// every contig mosdepth reads, so their `total_region` is exactly as
    /// long as `total`; a capture kit covers 1-3% of the genome. The second
    /// rule serves every path that writes no manifest row: run-all.sh with
    /// `--model_type WES`, a re-render with 24-html-report.sh, og live's
    /// re-read of a sample folder, and step 16b with `CAPTURE_BED`.
    pub fn read(path: &str, exome: bool) -> Result<Coverage, String> {
        // The `length` and `mean` of the last `total` and `total_region` rows.
        let (mut total, mut region) = (None, None);
        for r in read_tsv(path)? {
            let row = match r.get("chrom") {
                Some(Some("total")) => &mut total,
                Some(Some("total_region")) => &mut region,
                _ => continue,
            };
            let cell = |k: &str| r.get(k).flatten().map(str::to_owned);
            *row = Some((cell("length"), cell("mean")));
        }
        let length = |row: &Option<(Option<String>, Option<String>)>| {
            row.as_ref()
                .and_then(|(l, _)| l.as_deref())
                .and_then(|l| l.parse::<u64>().ok())
        };
        let capture = matches!((length(&region), length(&total)), (Some(r), Some(t)) if r.saturating_mul(2) < t);
        let on_target = exome || capture;
        let row = if on_target { region } else { total };
        let mean_depth = match row.and_then(|(_, m)| m) {
            Some(m) if !m.is_empty() => Some(float(&m).map_err(|e| py_error("ValueError", &e))?),
            _ => None,
        };
        Ok(Coverage {
            mean_depth,
            on_target,
        })
    }

    fn to_json(&self) -> Value {
        let mut pairs = vec![("mean_depth", opt_float(self.mean_depth))];
        if self.on_target {
            pairs.push(("on_target", Value::Bool(true)));
        }
        obj(pairs)
    }
}

// --- selected steps that left no result ------------------------------------------

/// What a section's note and its `not_run` line say of a step the run
/// selected that wrote nothing for the sample.
pub const SELECTED_NO_RESULT: &str = "selected, no result";

/// The sections whose step the run selected and that have no file, in the
/// order of `sections`, each with a note that says so.
///
/// The selection is the `run tools` row of `run_manifest.tsv`: the
/// `--tools` names, separated by spaces (HTML_REPORT writes it from
/// main.nf's tool list). `assets/steps.tsv` names the step of each tool and
/// the sections each step feeds; the one in `assets` is read, else the one
/// beside the binary. A section is marked when it is `missing` (a `failed`
/// one already says its step failed), a selected step feeds it, and it is
/// not in `skip` (the sections that are no step of their own).
/// Nextflow runs the Y haplogroup and Y-STR steps only for declared males.
///
/// `None` when the run did not say what it selected (no `run tools` row:
/// the bash runner writes none), or the catalogue cannot be read: the
/// sections then read as before, `missing` with no note.
pub fn selected_no_result(
    sections: &mut [Section],
    manifest: Option<&Manifest>,
    assets: &str,
    skip: &[&str],
    declared_sex: Option<&str>,
) -> Option<Vec<&'static str>> {
    let tools: Vec<&str> = manifest?
        .value("run", "tools")?
        .split_whitespace()
        .collect();
    let here = std::path::Path::new(assets).join("steps.tsv");
    let catalogue = if here.is_file() {
        here
    } else {
        og_runlog::steps::default_path()?
    };
    let steps = og_runlog::steps::read_steps(&catalogue).ok()?;
    let mut out = Vec::new();
    const MALE_ONLY: [&str; 3] = ["y_haplogroup", "y_tree", "y_str"];
    for s in sections.iter_mut().filter(|s| {
        let key = s.spec.key;
        s.state == State::Missing
            && !skip.contains(&key)
            && (declared_sex == Some("male") || !MALE_ONLY.contains(&key))
    }) {
        let Some(step) = steps.iter().find(|st| {
            !st.tool.is_empty()
                && tools.contains(&st.tool.as_str())
                && st.sections.iter().any(|k| k == s.spec.key)
        }) else {
            continue;
        };
        s.note = Some(format!(
            "{SELECTED_NO_RESULT}: --tools named {}, but step {} wrote nothing for this sample",
            step.tool, step.step
        ));
        out.push(s.spec.key);
    }
    Some(out)
}

// --- gene callability (step 16b) -----------------------------------------------------

/// A gene counts as read when this share of its bases reached the depth.
pub const CALLABLE_SHARE: f64 = 0.9;

/// The share of a gene's bases at the depth, and what it means.
#[derive(Debug, Clone, PartialEq)]
pub enum Reading {
    Read(f64),
    NotRead(f64),
    /// The body reached the share, a coding exon did not (step 16b's exon
    /// table): the gene is not read.
    ExonNotRead(f64),
    NotOnReference,
}

impl Reading {
    fn word(&self) -> &'static str {
        match self {
            Reading::Read(_) => "read",
            Reading::NotRead(_) => "not read",
            Reading::ExonNotRead(_) => "exon not read",
            Reading::NotOnReference => "not on this reference",
        }
    }

    fn share(&self) -> Value {
        match self {
            Reading::Read(s) | Reading::NotRead(s) | Reading::ExonNotRead(s) => Value::Float(*s),
            Reading::NotOnReference => Value::Null,
        }
    }
}

/// A gene's coding exon with the lowest share at the depth.
#[derive(Debug, Clone, PartialEq)]
pub struct WorstExon {
    /// `GENE:exonN`, as `assets/callability_exons.tsv` names it.
    pub exon: String,
    pub share: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Gene {
    pub gene: String,
    pub hgnc_id: String,
    pub reasons: Vec<String>,
    pub owner: Option<String>,
    /// The owner's step and whether it ran, for an owner [`owner_step`] knows.
    pub owner_run: Option<OwnerRun>,
    pub reading: Reading,
    /// Its worst coding exon, when step 16b measured exons and one was
    /// under the share.
    pub worst_exon: Option<WorstExon>,
}

/// The step that calls an owned gene, and what it left for this sample.
#[derive(Debug, Clone, PartialEq)]
pub struct OwnerRun {
    pub step: &'static str,
    pub status: OwnerStatus,
}

/// What the owner's step left in the sample folder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OwnerStatus {
    /// Its file is there, from this run.
    Ran,
    /// Its file is not there.
    NotRun,
    /// Its file is older than the last run and the step was not ok in that
    /// run: the rule that makes a section stale.
    Earlier,
}

impl OwnerStatus {
    fn word(self) -> &'static str {
        match self {
            OwnerStatus::Ran => "ran",
            OwnerStatus::NotRun => "not run",
            OwnerStatus::Earlier => "from an earlier run",
        }
    }
}

/// The step behind an owner of `assets/callability_genes.tsv` and the file
/// that step writes in the sample folder (`{s}` the sample, `{t}` the
/// gene's table in `assets/paralog_loci.tsv`) when it ran. `None` for an
/// owner this list does not know.
fn owner_step(owner: &str) -> Option<(&'static str, &'static str)> {
    Some(match owner {
        "parascopy" => ("35", "paralogs/{s}_{t}_copy_number.tsv"),
        "pypgx" => ("32", "pypgx/{s}_pypgx_summary.tsv"),
        "cn_loci:HBA" => ("35", "paralogs/{s}_hba_copy_number.tsv"),
        "cn_loci:CYP21A2" => ("35", "paralogs/{s}_cyp21a2_psv.tsv"),
        "cn_loci:HP" | "cn_loci:CCL3L1" | "cn_loci:LCE3BC" => ("35", "paralogs/{s}_cn_dosage.tsv"),
        "cn_loci:LPA_KIV2" => ("46", "lpa/{s}_lpa_kiv2.tsv"),
        "cn_loci:RHD" => ("41", "bloodtype/{s}_bloodtype.tsv"),
        _ => return None,
    })
}

/// The file the owner's step writes for this gene in `dir`, the sample
/// folder, and its step. A Parascopy gene without a table has no file:
/// step 35 never answers for it.
fn owner_file(
    r: &CatalogRow,
    dir: &std::path::Path,
    s: &str,
) -> Option<(&'static str, Option<std::path::PathBuf>)> {
    let (step, file) = owner_step(&r.owner)?;
    let file = if r.owner == "parascopy" {
        r.table.as_deref().map(|t| file.replace("{t}", t))
    } else {
        Some(file.to_owned())
    };
    Some((step, file.map(|f| dir.join(f.replace("{s}", s)))))
}

/// What the owner's step left at `file`, judged as a section is
/// ([`super::model::Section::new`]): a file older than the run's start, from
/// a step that was not ok in that run, is from an earlier run.
fn owner_status(
    step: &str,
    file: Option<&std::path::Path>,
    run: Option<&RunStatus>,
) -> OwnerStatus {
    let Some(mtime) = file
        .filter(|f| f.is_file())
        .and_then(|f| pycompat::mtime::mtime(f).ok())
    else {
        return OwnerStatus::NotRun;
    };
    if let Some(run) = run
        && let Some(started) = run.started_epoch
        && RUN_ALL_STEPS.contains(&step)
        && run.steps.get(step).map(String::as_str) != Some("ok")
        && is_stale(mtime, started)
    {
        return OwnerStatus::Earlier;
    }
    OwnerStatus::Ran
}

#[derive(Debug, Clone, PartialEq)]
pub struct Callability {
    pub min_depth: i128,
    pub genes: Vec<Gene>,
    /// The coding exons of catalog genes step 16b's exon table measured;
    /// `None` when there was no exon table (a folder from before it).
    pub exons_checked: Option<usize>,
}

/// The data rows of a mosdepth thresholds table of one depth: the depth,
/// and (chrom, region, bases, bases at the depth) of each row.
type Thresholds = (i128, Vec<(String, String, i128, i128)>);

fn thresholds(table: &str) -> Result<Thresholds, String> {
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
    let min_depth = int(&head[4][..head[4].len() - 1]).map_err(|e| py_error("ValueError", &e))?;
    let mut out = Vec::with_capacity(rows.len() - 1);
    for c in &rows[1..] {
        if c.len() != 5 {
            return Err(format!("{name}: a row has {} columns, not 5", c.len()));
        }
        let v = |s: &str| int(s).map_err(|e| py_error("ValueError", &e));
        let (start, end, ok) = (v(&c[1])?, v(&c[2])?, v(&c[4])?);
        out.push((c[0].clone(), c[3].clone(), end - start, ok));
    }
    Ok((min_depth, out))
}

/// The share of `n` bases that `ok` are, to four decimals; 0 for no bases.
fn share(n: i128, ok: i128) -> Result<f64, String> {
    if n == 0 {
        return Ok(0.0);
    }
    pycompat::fmt::round(ok as f64 / n as f64, 4)
        .map_err(|e| py_error(e.exception(), &e.to_string()))
}

/// Step 16b's exon table (`<s>_callability_exons.bed.gz`, rows named
/// `GENE:exonN`) read against the catalog: the worst coding exon of each
/// catalog gene it has rows for (the lower exon number on a tie), and how
/// many exon rows of catalog genes it holds. Rows of other genes are left
/// out: the table follows the catalog it was measured with.
fn worst_exons(
    catalog: &[CatalogRow],
    table: &str,
) -> Result<(std::collections::HashMap<String, WorstExon>, usize), String> {
    let name = basename(table);
    let (_, rows) = thresholds(table)?;
    let genes: std::collections::HashSet<&str> = catalog.iter().map(|r| r.gene.as_str()).collect();
    let number = |exon: &str| {
        exon.rsplit_once(":exon")
            .and_then(|(_, n)| n.parse::<u64>().ok())
    };
    let mut worst: std::collections::HashMap<String, WorstExon> = std::collections::HashMap::new();
    let mut checked = 0;
    for (_, region, n, ok) in rows {
        let Some((gene, _)) = region
            .rsplit_once(":exon")
            .filter(|_| number(&region).is_some())
        else {
            return Err(format!("{name}: {region:?} is not named GENE:exonN"));
        };
        if !genes.contains(gene) {
            continue;
        }
        checked += 1;
        let s = share(n, ok)?;
        let gene = gene.to_owned();
        let better = worst
            .get(&gene)
            .is_none_or(|w| (s, number(&region)) < (w.share, number(&w.exon)));
        if better {
            worst.insert(
                gene,
                WorstExon {
                    exon: region,
                    share: s,
                },
            );
        }
    }
    Ok((worst, checked))
}

/// One row of `assets/callability_genes.tsv`.
pub struct CatalogRow {
    pub gene: String,
    pub hgnc_id: String,
    pub chrom: String,
    pub start: String,
    pub reasons: String,
    pub owner: String,
    /// A Parascopy gene's table in `assets/paralog_loci.tsv`.
    pub table: Option<String>,
}

/// Reads the catalog: its comment lines left out, cells as written.
pub fn read_catalog(assets: &str) -> Result<Vec<CatalogRow>, String> {
    const NAME: &str = "callability_genes.tsv";
    let p = crate::join(assets, NAME);
    if !crate::is_file(&p) {
        return Err(format!(
            "{NAME} is not in {assets}: the report needs the repository's assets/ folder (og report render --assets)"
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
            table: None,
        });
    }
    // Step 35 writes one Parascopy table per row of paralog_loci.tsv, for
    // the genes the row names.
    if out.iter().any(|r| r.owner == "parascopy")
        && crate::is_file(crate::join(assets, "paralog_loci.tsv"))
    {
        let mut tables = std::collections::HashMap::new();
        for r in super::assets::asset_rows(assets, "paralog_loci.tsv")? {
            let cell = |k: &str| r.get(k).flatten().unwrap_or_default().to_owned();
            for g in cell("genes").split(',') {
                tables.entry(g.to_owned()).or_insert_with(|| cell("table"));
            }
        }
        for r in out.iter_mut().filter(|r| r.owner == "parascopy") {
            r.table = tables.get(&r.gene).cloned();
        }
    }
    Ok(out)
}

impl Callability {
    /// Step 16b's per-gene table (mosdepth `--thresholds` of one depth over
    /// the catalog's genes: chrom, start, end, gene, bases at the depth)
    /// read against the catalog.
    pub fn read(catalog: &[CatalogRow], table: &str) -> Result<Callability, String> {
        let name = basename(table);
        let (min_depth, rows) = thresholds(table)?;
        let mut bases: std::collections::HashMap<&str, (i128, i128)> =
            std::collections::HashMap::new();
        for (_, region, n, ok) in &rows {
            bases.insert(region, (*n, *ok));
        }
        if bases.is_empty() {
            return Err(format!("{name} has no gene rows"));
        }
        // mosdepth writes no row on a contig the BAM lacks: a chromosome with
        // no row at all is not on the sample's reference; a gene missing from
        // a chromosome that has rows means another catalog wrote the table.
        let contigs: std::collections::HashSet<&str> = rows.iter().map(|c| c.0.as_str()).collect();
        // The table is <sample folder>/mosdepth/ (bash) or coverage/ (Nextflow)
        // <s>_callability.bed.gz: the owner steps write into the same folder.
        let folder = std::path::Path::new(table)
            .parent()
            .and_then(|p| p.parent());
        let sample = name.strip_suffix("_callability.bed.gz");
        // run-all.sh's run status, which tells a file of this run from one an
        // earlier run left (Folder::open has read it already; a run status it
        // could not read stopped the summary before this reader).
        let run = folder
            .and_then(|d| std::fs::read_to_string(d.join("logs").join("run_status.tsv")).ok())
            .map(|t| RunStatus::parse(&t));
        // Step 16b's exon table beside it, when that step measured exons:
        // the worst coding exon of each catalog gene.
        let exon_table = sample.and_then(|s| {
            let p =
                std::path::Path::new(table).with_file_name(format!("{s}_callability_exons.bed.gz"));
            p.is_file().then(|| p.to_string_lossy().into_owned())
        });
        let (worst, exons_checked) = match &exon_table {
            Some(p) => {
                let (w, n) = worst_exons(catalog, p)?;
                (w, Some(n))
            }
            None => (std::collections::HashMap::new(), None),
        };
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
                let share = share(n, ok)?;
                if share < CALLABLE_SHARE {
                    Reading::NotRead(share)
                } else if worst
                    .get(r.gene.as_str())
                    .is_some_and(|w| w.share < CALLABLE_SHARE)
                {
                    Reading::ExonNotRead(share)
                } else {
                    Reading::Read(share)
                }
            };
            // The worst coding exon, named only when it was not read.
            let worst_exon = worst
                .get(r.gene.as_str())
                .filter(|w| w.share < CALLABLE_SHARE && reading != Reading::NotOnReference)
                .cloned();
            let owner = (r.owner != ".").then(|| r.owner.clone());
            let owner_run = match (folder, sample) {
                (Some(dir), Some(s)) => owner_file(r, dir, s).map(|(step, file)| OwnerRun {
                    step,
                    status: owner_status(step, file.as_deref(), run.as_ref()),
                }),
                _ => None,
            };
            genes.push(Gene {
                gene: r.gene.clone(),
                hgnc_id: r.hgnc_id.clone(),
                reasons: r.reasons.split(',').map(str::to_owned).collect(),
                owner,
                owner_run,
                reading,
                worst_exon,
            });
        }
        Ok(Callability {
            min_depth,
            genes,
            exons_checked,
        })
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
                        if owned {
                            let run = g.owner_run.as_ref();
                            m.push(("owner_step", opt_str(run.map(|r| r.step))));
                            m.push(("owner_status", opt_str(run.map(|r| r.status.word()))));
                        }
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
                        if let Some(w) = &g.worst_exon {
                            m.push((
                                "worst_exon",
                                obj([
                                    ("exon", Value::Str(w.exon.clone())),
                                    ("share", Value::Float(w.share)),
                                ]),
                            ));
                        }
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
        let mut m = vec![
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
        ];
        // Only when step 16b measured exons, so a folder from before keeps its summary.
        if let Some(n) = self.exons_checked {
            m.push(("exons_checked", Value::Int(n as i128)));
        }
        obj(m)
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
pub(crate) fn gzip_message(e: &std::io::Error) -> String {
    if e.kind() == std::io::ErrorKind::UnexpectedEof {
        "Compressed file ended before the end-of-stream marker was reached".to_owned()
    } else {
        e.to_string()
    }
}
