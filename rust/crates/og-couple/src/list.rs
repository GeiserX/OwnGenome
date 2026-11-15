//! The carrier list, `assets/carrier_genes.tsv`: one row per gene with its
//! GRCh38 span, how it is passed on, where it comes from, its condition,
//! the alleles that are not a full loss of function, and a dominant
//! condition of the gene besides the one it is listed for.

use crate::files::{each_line, is_digits};
use indexmap::IndexMap;
use og_io::Utf8;
use og_pedigree::Die;
use og_summary::collect::assets::{Allele, Key, parse_alleles};
use std::collections::HashMap;

/// The columns every list has, in order.
pub const GENE_COLUMNS: [&str; 7] = [
    "gene",
    "chrom",
    "start",
    "end",
    "inheritance",
    "lists",
    "condition",
];
/// The optional column after them: alleles that are not a full loss of
/// function.
pub const SEVERITY_COLUMN: &str = "severity";
/// The optional column after `severity`: a dominant condition of the gene
/// besides the one its inheritance gives ([`Dominant`]).
pub const DOMINANT_COLUMN: &str = "dominant";

/// Whether `cols` is a list's header: the gene columns, then optionally
/// `severity`, then optionally `dominant`.
pub fn header_ok<S: AsRef<str>>(cols: &[S]) -> bool {
    let n = GENE_COLUMNS.len();
    let more = [SEVERITY_COLUMN, DOMINANT_COLUMN];
    cols.len() >= n
        && cols.len() <= n + more.len()
        && cols.iter().map(AsRef::as_ref).eq(GENE_COLUMNS
            .iter()
            .chain(more.iter())
            .take(cols.len())
            .copied())
}

/// How a gene of the list is passed on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Inheritance {
    /// Autosomal recessive, or in a pseudoautosomal region.
    Ar,
    /// X-linked: chrX outside the PARs, the main disorder recessive.
    Xl,
    /// X-linked dominant: chrX outside the PARs, the main disorder
    /// dominant, so a daughter or a son who inherits a change can be
    /// affected.
    Xld,
}

impl Inheritance {
    pub fn parse(s: &str) -> Option<Inheritance> {
        match s {
            "AR" => Some(Inheritance::Ar),
            "XL" => Some(Inheritance::Xl),
            "XLD" => Some(Inheritance::Xld),
            _ => None,
        }
    }

    /// The value the list and the couple table write.
    pub fn code(self) -> &'static str {
        match self {
            Inheritance::Ar => "AR",
            Inheritance::Xl => "XL",
            Inheritance::Xld => "XLD",
        }
    }

    /// X-linked, recessive or dominant.
    pub fn x_linked(self) -> bool {
        self != Inheritance::Ar
    }
}

/// A dominant condition of a gene besides the one its inheritance gives:
/// the `dominant` column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dominant {
    /// `.`: none.
    No,
    /// `AD` or `SD`: the gene's condition on the ACMG secondary findings
    /// list (`assets/acmg_sf.tsv`) is autosomal dominant or semidominant,
    /// so one change can affect a child who inherits it.
    Autosomal,
    /// `XLD` on an X-linked gene: an Orphanet disorder of the gene is X-linked
    /// dominant, so females who carry a change can show signs.
    XLinked,
}

impl Dominant {
    pub fn parse(s: &str) -> Option<Dominant> {
        match s {
            "." => Some(Dominant::No),
            "AD" | "SD" => Some(Dominant::Autosomal),
            "XLD" => Some(Dominant::XLinked),
            _ => None,
        }
    }
}

/// One gene of the list.
#[derive(Debug, Clone)]
pub struct Gene {
    pub name: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub inheritance: Inheritance,
    pub lists: String,
    pub condition: String,
    /// The severity column: hypomorph and pseudodeficiency alleles.
    pub alleles: IndexMap<Key, Allele>,
    /// The dominant column.
    pub dominant: Dominant,
}

/// The list: every gene by name, a gene on two rows keeping its last.
#[derive(Debug, Clone, Default)]
pub struct CarrierList {
    pub genes: HashMap<String, Gene>,
}

impl CarrierList {
    /// Reads the list at `path`. A list without rows stops the screen: an
    /// empty list would report no shared gene for every couple.
    pub fn read(path: &str) -> Result<CarrierList, Die> {
        let mut genes = HashMap::new();
        let mut header: Option<Vec<String>> = None;
        let mut n = 0usize;
        each_line(path, Utf8::Strict, |raw| {
            n += 1;
            let line = raw.trim_end_matches(['\r', '\n']);
            if line.is_empty() || line.starts_with('#') {
                return Ok(());
            }
            let vals: Vec<&str> = line.split('\t').collect();
            let Some(head) = &header else {
                if !header_ok(&vals) {
                    return Err(Die::one(format!(
                        "{path}: the header must be {}, optionally followed by {SEVERITY_COLUMN} and \
                         {DOMINANT_COLUMN}",
                        GENE_COLUMNS.join("\t")
                    )));
                }
                header = Some(vals.iter().map(|v| (*v).to_owned()).collect());
                return Ok(());
            };
            if vals.len() != head.len() {
                return Err(Die::one(format!(
                    "{path} line {n}: {} fields, the header has {}",
                    vals.len(),
                    head.len()
                )));
            }
            let inheritance = Inheritance::parse(vals[4]);
            let dominant = Dominant::parse(vals.get(8).copied().unwrap_or("."));
            let (Some(inheritance), Some(dominant), true, true) = (
                inheritance,
                dominant,
                is_digits(vals[2]),
                is_digits(vals[3]),
            ) else {
                return Err(Die::one(format!(
                    "{path} line {n}: inheritance must be AR, XL or XLD, dominant ., AD, SD or XLD, start and \
                     end whole numbers"
                )));
            };
            let num = |s: &str| crate::files::int(s, &format!("{path} line {n}"));
            let alleles = parse_alleles(vals.get(7).copied().or(Some(".")), true)
                .map_err(|e| Die::one(format!("{path} line {n}: {e}")))?;
            let gene = Gene {
                name: vals[0].to_owned(),
                chrom: vals[1].to_owned(),
                start: num(vals[2])?,
                end: num(vals[3])?,
                inheritance,
                lists: vals[5].to_owned(),
                condition: vals[6].to_owned(),
                alleles,
                dominant,
            };
            genes.insert(gene.name.clone(), gene);
            Ok(())
        })?;
        if genes.is_empty() {
            return Err(Die::one(format!(
                "{path} holds no gene: the screen would find nothing. Restore assets/carrier_genes.tsv."
            )));
        }
        Ok(CarrierList { genes })
    }

    pub fn get(&self, gene: &str) -> Option<&Gene> {
        self.genes.get(gene)
    }

    pub fn contains(&self, gene: &str) -> bool {
        self.genes.contains_key(gene)
    }

    pub fn len(&self) -> usize {
        self.genes.len()
    }
}

/// One chromosome's spans sorted by start, end, name, and its longest span.
type Spans = (Vec<(i64, i64, String)>, i64);

/// The genes of the list by chromosome, for the genes a deletion overlaps.
pub struct GeneIndex {
    by: HashMap<String, Spans>,
}

impl GeneIndex {
    pub fn new(list: &CarrierList) -> GeneIndex {
        let mut by: HashMap<String, Spans> = HashMap::new();
        for g in list.genes.values() {
            let e = by.entry(g.chrom.clone()).or_default();
            e.0.push((g.start, g.end, g.name.clone()));
            e.1 = e.1.max(g.end - g.start);
        }
        for (rows, _) in by.values_mut() {
            rows.sort();
        }
        GeneIndex { by }
    }

    /// The genes whose 1-based span overlaps the 1-based `[s, e]`, in
    /// order of start.
    pub fn overlapping(&self, chrom: &str, s: i64, e: i64) -> Vec<(i64, i64, &str)> {
        let Some((rows, longest)) = self.by.get(chrom) else {
            return Vec::new();
        };
        let lo = rows.partition_point(|r| r.0 < s - longest - 1);
        let hi = rows.partition_point(|r| r.0 <= e);
        rows[lo..hi.max(lo)]
            .iter()
            .filter(|r| r.0 <= e && r.1 >= s)
            .map(|r| (r.0, r.1, r.2.as_str()))
            .collect()
    }
}
