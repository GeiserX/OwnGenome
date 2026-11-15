//! `og couple genes`: rebuilds the carrier list, `assets/carrier_genes.tsv`.
//!
//! ```text
//! og couple genes --table assets/carrier_genes.tsv --orphanet-genes en_product6.xml
//!                 --orphanet-inheritance en_product9_ages.xml
//!                 --gtf gencode.v<N>.basic.annotation.gtf.gz --par-bed assets/par_grch38.bed --out NEW.tsv
//! ```
//!
//! The ACMG and ACOG rows of `--table` are kept as they are written, every
//! other row is made again from Orphanet (the genes of disorders with
//! autosomal recessive or X-linked inheritance), and every row gets its
//! GRCh38 span from the GENCODE gene lines. Each input is read once, as a
//! stream.
//!
//! An X-linked gene is `XLD` when its main disorder is X-linked dominant:
//! every X-linked disorder Orphanet gives the gene is dominant, or the gene
//! is one of [`MAIN_DOMINANT`]. Every other X-linked gene is `XL`; when one
//! of its disorders is X-linked dominant, its `dominant` cell is `XLD`
//! (females who carry a change can show signs). The `severity` cell and an
//! autosomal gene's `dominant` cell (`AD` or `SD`, from the ACMG secondary
//! findings list) are written by hand and kept.

use crate::args::{self, Kind, Spec};
use crate::files::{each_line, int};
use crate::list::{DOMINANT_COLUMN, GENE_COLUMNS, SEVERITY_COLUMN, header_ok};
use crate::xml::{self, Event};
use indexmap::IndexMap;
use og_io::Utf8;
use og_pedigree::{Die, Par};
use std::collections::{BTreeSet, HashMap};
use std::io::Write;

const USAGE: &str =
    "usage: og couple genes --table TABLE --orphanet-genes XML --orphanet-inheritance XML
                       --gtf GTF --par-bed BED --out OUT
";

const SPEC: Spec = Spec {
    prog: "og couple genes",
    usage: USAGE,
    options: &[
        ("--table", Kind::One),
        ("--orphanet-genes", Kind::One),
        ("--orphanet-inheritance", Kind::One),
        ("--gtf", Kind::One),
        ("--par-bed", Kind::One),
        ("--out", Kind::One),
    ],
    required: &[
        "--table",
        "--orphanet-genes",
        "--orphanet-inheritance",
        "--gtf",
        "--par-bed",
        "--out",
    ],
};

/// The Orphanet gene associations kept: a disease-causing germline change.
const KEEP: [&str; 3] = [
    "Disease-causing germline mutation(s) in",
    "Disease-causing germline mutation(s) (loss of function) in",
    "Disease-causing germline mutation(s) (gain of function) in",
];

/// The chromosomes a span may be on, in the list's order.
const PRIMARY: [&str; 24] = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22", "chrX", "chrY",
];

fn primary(chrom: &str) -> Option<usize> {
    PRIMARY.iter().position(|c| *c == chrom)
}

/// The genes with both X-linked dominant and X-linked recessive disorders
/// whose main disorder is the dominant one, with that disorder. A gene whose
/// X-linked disorders are all dominant needs no entry.
pub const MAIN_DOMINANT: [(&str, &str); 4] = [
    ("CDKL5", "CDKL5-deficiency disorder"),
    ("IKBKG", "Incontinentia pigmenti"),
    ("MECP2", "Rett syndrome"),
    ("OFD1", "Orofaciodigital syndrome type 1"),
];

/// An Orphanet mode of inheritance the list keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Mode {
    Ar,
    Xlr,
    Xld,
}

fn mode(name: &str) -> Option<Mode> {
    match name {
        "Autosomal recessive" => Some(Mode::Ar),
        "X-linked recessive" => Some(Mode::Xlr),
        "X-linked dominant" => Some(Mode::Xld),
        _ => None,
    }
}

/// `og couple genes ARGS`.
pub fn cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match args::parse(&SPEC, argv, out, err) {
        Ok(a) => a,
        Err(code) => return code,
    };
    let one = |k: &str| a.one(k).unwrap_or_default();
    match run(
        one("--table"),
        one("--orphanet-genes"),
        one("--orphanet-inheritance"),
        one("--gtf"),
        one("--par-bed"),
        one("--out"),
    ) {
        Ok(line) => {
            let _ = writeln!(out, "{line}");
            0
        }
        Err(die) => {
            for line in &die.lines {
                let _ = writeln!(err, "ERROR: {line}");
            }
            die.code
        }
    }
}

/// A GENCODE gene's place: chromosome, 1-based start and end.
type Place = (String, i64, i64);

/// The first `gene_name "..."` of a GTF attribute field.
fn gene_name(attrs: &str) -> Option<&str> {
    let mut rest = attrs;
    while let Some(i) = rest.find("gene_name \"") {
        let v = &rest[i + 11..];
        if let Some(j) = v.find('"')
            && j > 0
        {
            return Some(&v[..j]);
        }
        rest = v;
    }
    None
}

/// GENCODE's gene lines: each name's place (a name on several
/// chromosomes keeps the first of chr1 to chrY, so a PAR gene is on chrX),
/// in the order the file first names them, and the release.
fn gencode(path: &str) -> Result<(IndexMap<String, Place>, String), Die> {
    let mut found: IndexMap<String, Vec<(usize, i64, i64)>> = IndexMap::new();
    let mut release = "?".to_owned();
    each_line(path, Utf8::Replace, |line| {
        if line.starts_with('#') {
            // re.search(r"version (\d+)", line)
            let found = line.match_indices("version ").find_map(|(i, _)| {
                let v = &line[i + 8..];
                let n = v.bytes().take_while(u8::is_ascii_digit).count();
                (n > 0).then(|| v[..n].to_owned())
            });
            if let Some(v) = found {
                release = v;
            }
            return Ok(());
        }
        let c: Vec<&str> = line.splitn(10, '\t').collect();
        if c.len() < 9 || c[2] != "gene" {
            return Ok(());
        }
        let Some(chrom) = primary(c[0]) else {
            return Ok(());
        };
        if let Some(name) = gene_name(c[8]) {
            found.entry(name.to_owned()).or_default().push((
                chrom,
                int(c[3], path)?,
                int(c[4], path)?,
            ));
        }
        Ok(())
    })?;
    let out = found
        .into_iter()
        .map(|(g, mut locs)| {
            locs.sort_by_key(|&(c, s, _)| (c, s));
            let (c, s, e) = locs[0];
            (g, (PRIMARY[c].to_owned(), s, e))
        })
        .collect();
    Ok((out, release))
}

/// The kept gene associations of a disorder: (type, status, symbol).
type Raw = (String, String, String);

/// One element being read: its name, its own text before its first child,
/// the first texts of the child paths the reader asks for, and, for a
/// Disorder, the modes and associations found inside it.
#[derive(Default)]
struct Frame {
    name: String,
    text: String,
    child: bool,
    first: HashMap<String, String>,
    modes: BTreeSet<Mode>,
    assoc: Vec<Raw>,
}

impl Frame {
    /// `findtext(path)`: the first text of a kept child path, "" if none.
    fn get(&self, path: &str) -> String {
        self.first.get(path).cloned().unwrap_or_default()
    }
}

/// The child paths whose first text is kept, as `findtext` reads them.
const PATHS: [&str; 6] = [
    "OrphaCode",
    "Name",
    "Symbol",
    "DisorderGeneAssociationType/Name",
    "DisorderGeneAssociationStatus/Name",
    "Gene/Symbol",
];

/// Reads an Orphadata product. Each `TypeOfInheritance` adds its mode, and
/// each `DisorderGeneAssociation` its association, to every Disorder around
/// it; `disorder` is called with each Disorder as it closes. The root's
/// `date` attribute is returned.
fn walk(path: &str, mut disorder: impl FnMut(&Frame)) -> Result<Option<String>, Die> {
    let mut stack: Vec<Frame> = Vec::new();
    let mut date: Option<Option<String>> = None;
    xml::each_event(path, |e| {
        match e {
            Event::Start(name, attrs) => {
                if date.is_none() {
                    date = Some(xml::attr(attrs, "date"));
                }
                if let Some(top) = stack.last_mut() {
                    top.child = true;
                }
                stack.push(Frame {
                    name: name.to_owned(),
                    ..Frame::default()
                });
            }
            Event::Text(t) => {
                if let Some(top) = stack.last_mut()
                    && !top.child
                {
                    top.text.push_str(&t);
                }
            }
            Event::End(_) => {
                let Some(done) = stack.pop() else {
                    return Err(Die::one(format!(
                        "{path}: not well-formed XML (an end tag too many)"
                    )));
                };
                let around = stack.iter_mut().filter(|f| f.name == "Disorder");
                match done.name.as_str() {
                    "TypeOfInheritance" => {
                        if let Some(m) = mode(&done.get("Name")) {
                            around.for_each(|d| {
                                d.modes.insert(m);
                            });
                        }
                    }
                    "DisorderGeneAssociation" => {
                        let a = (
                            done.get("DisorderGeneAssociationType/Name"),
                            done.get("DisorderGeneAssociationStatus/Name"),
                            done.get("Gene/Symbol"),
                        );
                        around.for_each(|d| d.assoc.push(a.clone()));
                    }
                    "Disorder" => disorder(&done),
                    _ => {}
                }
                if let Some(parent) = stack.last_mut() {
                    if PATHS.contains(&done.name.as_str()) {
                        parent
                            .first
                            .entry(done.name.clone())
                            .or_insert_with(|| done.text.clone());
                    }
                    for (k, v) in &done.first {
                        let key = format!("{}/{k}", done.name);
                        if PATHS.contains(&key.as_str()) {
                            parent.first.entry(key).or_insert_with(|| v.clone());
                        }
                    }
                }
            }
        }
        Ok(())
    })?;
    Ok(date.flatten())
}

/// One kept association of a gene: the disorder's name and its modes.
type Association = (String, BTreeSet<Mode>);

/// The Orphanet modes of inheritance of each disorder (product 9), by
/// OrphaCode.
fn inheritance(path: &str) -> Result<HashMap<String, BTreeSet<Mode>>, Die> {
    let mut out = HashMap::new();
    walk(path, |d| {
        out.insert(d.get("OrphaCode"), d.modes.clone());
    })?;
    Ok(out)
}

/// Gene associations by gene symbol, in first-named order.
type ByGene = IndexMap<String, Vec<Association>>;

/// The kept gene associations of every disorder with a kept mode (product
/// 6), by gene symbol in the order the product first names them, and the
/// product's date.
fn associations(
    path: &str,
    inh: &HashMap<String, BTreeSet<Mode>>,
) -> Result<(ByGene, Option<String>), Die> {
    let mut out = ByGene::new();
    let date = walk(path, |d| {
        let Some(modes) = inh.get(&d.get("OrphaCode")).filter(|m| !m.is_empty()) else {
            return;
        };
        for (kind, status, symbol) in &d.assoc {
            if KEEP.contains(&kind.as_str()) && status == "Assessed" {
                out.entry(symbol.clone())
                    .or_default()
                    .push((d.get("Name"), modes.clone()));
            }
        }
    })?;
    Ok((out, date))
}

fn run(
    table: &str,
    p6: &str,
    p9: &str,
    gtf: &str,
    par_bed: &str,
    out_path: &str,
) -> Result<String, Die> {
    let mut lines: Vec<Vec<String>> = Vec::new();
    each_line(table, Utf8::Strict, |line| {
        if !line.starts_with('#') {
            lines.push(
                line.trim_end_matches('\r')
                    .split('\t')
                    .map(str::to_owned)
                    .collect(),
            );
        }
        Ok(())
    })?;
    let cols = lines.first().cloned().unwrap_or_default();
    if !header_ok(&cols) {
        return Err(Die::one(format!(
            "{table}: the header must be {}, optionally followed by {SEVERITY_COLUMN} and {DOMINANT_COLUMN}",
            GENE_COLUMNS.join("\t")
        )));
    }
    let rows_in: Vec<HashMap<&str, &str>> = lines[1..]
        .iter()
        .filter(|r| !r[0].is_empty())
        .map(|r| {
            cols.iter()
                .map(String::as_str)
                .zip(r.iter().map(String::as_str))
                .collect()
        })
        .collect();
    let cell = |r: &HashMap<&str, &str>, k: &str| r.get(k).copied().unwrap_or_default().to_owned();
    // The severity cell, and the dominant cell of an autosomal gene, are
    // written by hand: every gene keeps its own.
    let mut severity: HashMap<String, String> = HashMap::new();
    let mut dominant: HashMap<String, String> = HashMap::new();
    for r in &rows_in {
        let s = cell(r, SEVERITY_COLUMN);
        severity.insert(
            cell(r, "gene"),
            if s.is_empty() { ".".to_owned() } else { s },
        );
        let d = cell(r, DOMINANT_COLUMN);
        if matches!(d.as_str(), "AD" | "SD") {
            dominant.insert(cell(r, "gene"), d);
        }
    }
    let panel: Vec<&HashMap<&str, &str>> = rows_in
        .iter()
        .filter(|r| {
            let l = cell(r, "lists");
            l.contains("ACMG") || l.contains("ACOG")
        })
        .collect();
    let (genes, release) = gencode(gtf)?;
    let mut upper: HashMap<String, &str> = HashMap::new();
    for g in genes.keys() {
        upper.insert(pycompat::text::upper(g), g);
    }
    let par = Par::read(Some(par_bed))?;
    let canonical = |g: &str| -> Option<&str> {
        if genes.contains_key(g) {
            genes.get_key_value(g).map(|(k, _)| k.as_str())
        } else {
            upper.get(&pycompat::text::upper(g)).copied()
        }
    };
    let place = |g: &str| canonical(g).and_then(|k| genes.get(k));
    let missing: Vec<String> = panel
        .iter()
        .map(|r| cell(r, "gene"))
        .filter(|g| place(g).is_none())
        .collect();
    if !missing.is_empty() {
        return Err(Die::one(format!(
            "GENCODE has no gene line for {}; give the current HGNC symbol",
            missing.join(", ")
        )));
    }
    let inh = inheritance(p9)?;
    let (disorders, date) = associations(p6, &inh)?;

    // The X-linked disorders of each listed gene: (name, dominant, recessive).
    let mut x_of: HashMap<&str, Vec<(&str, bool, bool)>> = HashMap::new();
    for (sym, ds) in &disorders {
        let Some(gene) = canonical(sym) else { continue };
        if x_of.contains_key(gene) {
            continue;
        }
        let x: Vec<(&str, bool, bool)> = ds
            .iter()
            .filter(|(_, m)| m.contains(&Mode::Xld) || m.contains(&Mode::Xlr))
            .map(|(n, m)| (n.as_str(), m.contains(&Mode::Xld), m.contains(&Mode::Xlr)))
            .collect();
        x_of.insert(gene, x);
    }
    let some_dominant = |g: &str| x_of.get(g).is_some_and(|x| x.iter().any(|d| d.1));
    // gene -> its seven cells, in order.
    let mut rows: IndexMap<String, [String; 7]> = IndexMap::new();
    for r in &panel {
        let (c, s, e) = place(&cell(r, "gene")).expect("checked above");
        rows.insert(
            cell(r, "gene"),
            [
                cell(r, "gene"),
                c.clone(),
                s.to_string(),
                e.to_string(),
                cell(r, "inheritance"),
                cell(r, "lists"),
                cell(r, "condition"),
            ],
        );
    }
    for (sym, ds) in &disorders {
        let Some(gene) = canonical(sym) else { continue };
        if rows.contains_key(gene) {
            continue;
        }
        let (c, s, e) = &genes[gene];
        // The chromosome decides how a gene is passed on: X outside the
        // PARs is X-linked, an autosome or a PAR is autosomal.
        let x = c == "chrX" && !par.contains(c, *s);
        let names: BTreeSet<&str> = ds
            .iter()
            .filter(|(_, modes)| {
                if x {
                    modes.contains(&Mode::Xlr) || modes.contains(&Mode::Xld)
                } else {
                    modes.contains(&Mode::Ar)
                }
            })
            .map(|(n, _)| n.as_str())
            .collect();
        if names.is_empty() || c == "chrY" {
            continue;
        }
        // Rule (y): X-linked dominant only when the main disorder is.
        let xd = &x_of[gene];
        let xld = x
            && xd.iter().any(|d| d.1)
            && (xd.iter().all(|d| !d.2) || MAIN_DOMINANT.iter().any(|(g, _)| *g == gene));
        // An X-linked dominant gene is named by a dominant disorder first.
        let mut names: Vec<&str> = names.into_iter().collect();
        if xld {
            names.sort_by_key(|n| !xd.iter().any(|d| d.0 == *n && d.1));
        }
        let first = names[0];
        let cond = match names.len() {
            1 => first.to_owned(),
            n => format!(
                "{first}; and {} other{}",
                n - 1,
                if n == 2 { "" } else { "s" }
            ),
        };
        rows.insert(
            gene.to_owned(),
            [
                gene.to_owned(),
                c.clone(),
                s.to_string(),
                e.to_string(),
                if xld {
                    "XLD"
                } else if x {
                    "XL"
                } else {
                    "AR"
                }
                .to_owned(),
                "Orphanet".to_owned(),
                cond.replace('\t', " "),
            ],
        );
    }
    let date6 = date
        .as_deref()
        .and_then(|d| d.split_whitespace().next())
        .unwrap_or("?")
        .to_owned();
    let mut text = header(&release, &date6);
    text.push_str(&GENE_COLUMNS.join("\t"));
    text.push('\t');
    text.push_str(SEVERITY_COLUMN);
    text.push('\t');
    text.push_str(DOMINANT_COLUMN);
    text.push('\n');
    let mut order: Vec<&String> = rows.keys().collect();
    order.sort_by_key(|g| {
        let r = &rows[*g];
        (primary(&r[1]), r[2].parse::<i64>().unwrap_or(0), *g)
    });
    for g in order {
        text.push_str(&rows[g].join("\t"));
        text.push('\t');
        text.push_str(severity.get(g).map_or(".", String::as_str));
        text.push('\t');
        let d = match rows[g][4].as_str() {
            "XL" if some_dominant(g) => "XLD",
            "AR" => dominant.get(g.as_str()).map_or(".", String::as_str),
            _ => ".",
        };
        text.push_str(d);
        text.push('\n');
    }
    let mut f = og_io::AtomicFile::create(out_path).map_err(|e| Die::io(out_path, &e))?;
    f.write_all(text.as_bytes())
        .and_then(|()| f.commit())
        .map_err(|e| Die::io(out_path, &e))?;
    Ok(format!(
        "{} genes: {} ACMG or ACOG, {} Orphanet",
        rows.len(),
        panel.len(),
        rows.len() - panel.len()
    ))
}

/// The comment lines at the top of the list.
fn header(release: &str, date6: &str) -> String {
    [
        "# carrier_genes.tsv: the genes of the couple carrier screen (step 39, og couple screen).".to_owned(),
        "# One row per gene. chrom, start and end: the gene's GRCh38 span (1-based, inclusive) from the".to_owned(),
        format!("# gene lines of GENCODE release {release} (basic annotation). inheritance: AR (autosomal, or in a"),
        "# pseudoautosomal region), XL (chrX outside the PARs) or XLD (X-linked with a dominant main".to_owned(),
        "# disorder: every X-linked disorder Orphanet gives the gene is dominant, or the gene is CDKL5,".to_owned(),
        "# IKBKG, MECP2 or OFD1). lists: where the gene comes from:".to_owned(),
        "#   ACMG      Gregg et al., Genet Med 2021;23(10):1793-1806 (doi:10.1038/s41436-021-01203-z),".to_owned(),
        "#             tier 3: the 97 autosomal recessive genes of its tables 1 to 5 and the 16 X-linked".to_owned(),
        "#             genes of table 6; condition names from those tables".to_owned(),
        "#   ACOG      the conditions of ACOG Committee Opinions 690 and 691 (2017): cystic fibrosis, spinal".to_owned(),
        "#             muscular atrophy, the hemoglobinopathies, fragile X and the Ashkenazi Jewish panel".to_owned(),
        "#   Orphanet  every other gene of an Orphanet disorder with autosomal recessive or X-linked".to_owned(),
        format!("#             inheritance (Orphadata en_product6 and en_product9_ages of {date6}, CC BY 4.0,"),
        "#             https://www.orphadata.com), with a disease-causing germline association that".to_owned(),
        "#             Orphanet assessed; condition: the first such disorder by name (for XLD, the first".to_owned(),
        "#             X-linked dominant one)".to_owned(),
        "# severity: alleles of the gene that are not a full loss of function, separated by \";\": each is".to_owned(),
        "#   hypomorph or pseudodeficiency, then chrom:pos:ref:alt (GRCh38), then a name; \".\" for none. Such".to_owned(),
        "#   an allele never makes a carrier finding on its own (the couple screen, step 39, and the".to_owned(),
        "#   carrier list of the sample report). With a called deletion of the gene it is a second hit (the".to_owned(),
        "#   two-hit lines of the sample report, og summary collect). In the couple screen a".to_owned(),
        "#   hypomorph beside a full loss of function in the other partner's same gene is a shared gene.".to_owned(),
        "# dominant: a dominant condition of the gene besides the one it is listed for. AD or SD: the".to_owned(),
        "#   gene's condition on the ACMG secondary findings list (assets/acmg_sf.tsv) is autosomal".to_owned(),
        "#   dominant or semidominant, written by hand; XLD: an XL gene with an Orphanet X-linked dominant".to_owned(),
        "#   disorder, so females who carry a change can show signs; \".\" for none.".to_owned(),
        "# Rebuilt with: og couple genes (docs/39-couple.md says how).".to_owned(),
    ]
    .iter()
    .map(|l| format!("{l}\n"))
    .collect()
}
