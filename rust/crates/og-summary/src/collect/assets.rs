//! The two gene lists the clinical sections read from `assets/`:
//! `acmg_sf.tsv` (the ACMG SF v3.3 genes and the reporting rule of each,
//! Lee et al., Genet Med 2025, PMID 40568962) and `carrier_genes.tsv` (the
//! genes of the carrier list and their spans).
//!
//! A missing or broken list is an error with the words the reports print;
//! the section that needs it says the list is unavailable and keeps its
//! other values.

use crate::{basename, join};
use indexmap::IndexMap;
use pycompat::csv::{DictReader, DictRow, Quoting};

/// The list version both reports print.
pub const ACMG_SF_VERSION: &str = "ACMG SF v3.3";
/// The rules a gene of the ACMG list can have.
pub const SF_REPORT_IF: [&str; 5] = [
    "any",
    "biallelic",
    "hom_variant",
    "hemizygous",
    "het_female_too",
];
/// The kinds of a severity entry.
pub const SEVERITY_KINDS: [&str; 2] = ["hypomorph", "pseudodeficiency"];

/// A variant as `(chrom, pos, ref, alt)`.
pub type Key = (String, i128, String, String);

/// What a severity entry or a variants entry names.
#[derive(Debug, Clone, PartialEq)]
pub struct Allele {
    /// `hypomorph` or `pseudodeficiency`; `None` in a variants cell.
    pub kind: Option<String>,
    pub name: String,
}

/// The reporting rule of one ACMG gene.
#[derive(Debug, Clone)]
pub struct Rule {
    pub report_if: String,
    pub variants: IndexMap<Key, Allele>,
    /// AnnotSV's class 4 and 5 SVs count in this gene.
    pub sv_cn: bool,
    /// Regions where a deletion counts for the gene: chrom, start, end, name.
    pub sv_regions: Vec<(String, i128, i128, String)>,
    pub severity: IndexMap<Key, Allele>,
    pub table_text: String,
    pub note: Option<String>,
}

impl Rule {
    pub fn severity_kind(&self, key: &Key) -> Option<&str> {
        self.severity.get(key).and_then(|a| a.kind.as_deref())
    }
}

/// One gene of the carrier list.
#[derive(Debug, Clone)]
pub struct CarrierGene {
    pub inheritance: String,
    pub condition: String,
    pub severity: IndexMap<Key, Allele>,
    /// GRCh38, 1-based, inclusive.
    pub span: (String, i128, i128),
}

/// The rows of `assets/<name>` without its `#` lines, cells as written
/// (no quoting).
pub fn asset_rows(assets: &str, name: &str) -> Result<Vec<DictRow>, String> {
    let p = join(assets, name);
    let missing = || {
        format!(
            "{name} is not in {assets}: the report needs the repository's assets/ folder \
             (render_report.py --assets)"
        )
    };
    if !crate::is_file(&p) {
        return Err(missing());
    }
    let bytes = std::fs::read(&p).map_err(|e| crate::os_error(&e, &p))?;
    let text = String::from_utf8(bytes).map_err(|e| {
        format!(
            "UnicodeDecodeError: 'utf-8' codec can't decode {}: {e}",
            basename(&p)
        )
    })?;
    let lines = text
        .split_inclusive('\n')
        .filter(|l| !l.starts_with('#'))
        .map(|l| Ok::<_, std::io::Error>(l.to_owned()));
    DictReader::new(lines, '\t', Quoting::None)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("{}: {e}", e.exception()))
}

/// `row[key]`, its error as `str(KeyError)`.
fn col<'a>(r: &'a DictRow, key: &str) -> Result<Option<&'a str>, String> {
    r.get(key).ok_or_else(|| crate::repr(key))
}

fn int(s: Option<&str>) -> Result<i128, String> {
    match s {
        None => Err(
            "int() argument must be a string, a bytes-like object or a real number, not 'NoneType'"
                .into(),
        ),
        Some(s) => crate::int(s),
    }
}

/// `'chr6:26092913:G:A'` as a key.
pub fn allele_key(locus: &str) -> Result<Key, String> {
    let f: Vec<&str> = locus.split(':').collect();
    match f.as_slice() {
        [c, pos, r, a] => Ok((
            (*c).to_owned(),
            crate::int(pos)?,
            (*r).to_owned(),
            (*a).to_owned(),
        )),
        _ if f.len() < 4 => Err(format!(
            "not enough values to unpack (expected 4, got {})",
            f.len()
        )),
        _ => Err(format!(
            "too many values to unpack (expected 4, got {})",
            f.len()
        )),
    }
}

/// The variants or severity cell: `.` or entries separated by `;`, each
/// `[kind ]chrom:pos:ref:alt name`.
pub fn parse_alleles(text: Option<&str>, kinds: bool) -> Result<IndexMap<Key, Allele>, String> {
    let mut out = IndexMap::new();
    for e in text.filter(|t| !t.is_empty()).unwrap_or(".").split(';') {
        let e = pycompat::text::strip(e, None);
        if e.is_empty() || e == "." {
            continue;
        }
        let mut parts: Vec<&str> = e.splitn(if kinds { 3 } else { 2 }, ' ').collect();
        if kinds && parts.len() < 2 {
            return Err(format!(
                "severity entry needs a kind and a chrom:pos:ref:alt: {e}"
            ));
        }
        let kind = if kinds {
            let k = parts.remove(0);
            if !SEVERITY_KINDS.contains(&k) {
                return Err(format!(
                    "severity '{k}' is not one of {}: {e}",
                    SEVERITY_KINDS.join(", ")
                ));
            }
            Some(k.to_owned())
        } else {
            None
        };
        let name = parts.get(1).copied().unwrap_or(parts[0]).to_owned();
        out.insert(allele_key(parts[0])?, Allele { kind, name });
    }
    Ok(out)
}

/// `{gene: rule}` of `acmg_sf.tsv`, in the file's order.
pub fn acmg_sf_rules(assets: &str) -> Result<IndexMap<String, Rule>, String> {
    let mut rules: IndexMap<String, Rule> = IndexMap::new();
    for (n, r) in asset_rows(assets, "acmg_sf.tsv")?.iter().enumerate() {
        let one = || -> Result<(String, Rule), String> {
            let g = col(r, "gene")?.unwrap_or("None").to_owned();
            let report_if = col(r, "report_if")?;
            let sv_cn = col(r, "sv_cn")?;
            let ri_ok = report_if.is_some_and(|x| SF_REPORT_IF.contains(&x));
            if !ri_ok || !matches!(sv_cn, Some("yes" | "no")) {
                return Err(format!(
                    "report_if must be one of {}, sv_cn yes or no",
                    SF_REPORT_IF.join(", ")
                ));
            }
            let report_if = report_if.unwrap_or_default();
            if sv_cn == Some("yes") && !matches!(report_if, "any" | "het_female_too") {
                return Err(
                    "sv_cn yes needs a rule that one variant meets (any, het_female_too)".into(),
                );
            }
            let mut regions = Vec::new();
            for e in col(r, "sv_regions")?
                .filter(|x| !x.is_empty())
                .unwrap_or(".")
                .split(';')
            {
                let e = pycompat::text::strip(e, None);
                if e.is_empty() || e == "." {
                    continue;
                }
                let (loc, name) = e.split_once(' ').unwrap_or((e, ""));
                let (c, span) = loc.split_once(':').unwrap_or((loc, ""));
                let (a, b) = span.split_once('-').unwrap_or((span, ""));
                regions.push((
                    c.to_owned(),
                    crate::int(a)?,
                    crate::int(b)?,
                    if name.is_empty() { loc } else { name }.to_owned(),
                ));
            }
            let variants = parse_alleles(col(r, "variants")?, false)?;
            if report_if == "hom_variant" && variants.is_empty() {
                return Err("report_if hom_variant needs its variants".into());
            }
            if rules.contains_key(&g) {
                return Err("the gene is on two rows".into());
            }
            let inheritance = col(r, "inheritance")?;
            let _ = inheritance;
            let severity = parse_alleles(col(r, "severity")?, true)?;
            let table_text = col(r, "table_text")?.unwrap_or("None").to_owned();
            let note = col(r, "note")?
                .filter(|x| !matches!(*x, "" | "."))
                .map(str::to_owned);
            Ok((
                g,
                Rule {
                    report_if: report_if.to_owned(),
                    variants,
                    sv_cn: sv_cn == Some("yes"),
                    sv_regions: regions,
                    severity,
                    table_text,
                    note,
                },
            ))
        };
        let (g, rule) = one().map_err(|e| format!("acmg_sf.tsv row {}: {e}", n + 1))?;
        rules.insert(g, rule);
    }
    if rules.is_empty() {
        return Err("acmg_sf.tsv has no gene".into());
    }
    Ok(rules)
}

/// `{gene: row}` of `carrier_genes.tsv`; a gene on two rows keeps its
/// first place and its last row.
pub fn carrier_genes(assets: &str) -> Result<IndexMap<String, CarrierGene>, String> {
    let mut out = IndexMap::new();
    for (n, r) in asset_rows(assets, "carrier_genes.tsv")?.iter().enumerate() {
        let one = || -> Result<(String, CarrierGene), String> {
            let g = col(r, "gene")?.unwrap_or("None").to_owned();
            let inheritance = col(r, "inheritance")?.unwrap_or("None").to_owned();
            col(r, "lists")?;
            let condition = col(r, "condition")?.unwrap_or("None").to_owned();
            let severity = parse_alleles(r.get("severity").flatten(), true)?;
            let chrom = col(r, "chrom")?.unwrap_or("None").to_owned();
            let span = (chrom, int(col(r, "start")?)?, int(col(r, "end")?)?);
            Ok((
                g,
                CarrierGene {
                    inheritance,
                    condition,
                    severity,
                    span,
                },
            ))
        };
        let (g, x) = one().map_err(|e| format!("carrier_genes.tsv row {}: {e}", n + 1))?;
        out.insert(g, x);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alleles_read_with_and_without_kinds() {
        let v = parse_alleles(Some("chr6:10:G:A C282Y; chr6:20:C:G"), false).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(
            v[&("chr6".into(), 10, "G".into(), "A".into())].name,
            "C282Y"
        );
        assert_eq!(
            v[&("chr6".into(), 20, "C".into(), "G".into())].name,
            "chr6:20:C:G"
        );
        let s = parse_alleles(Some("hypomorph chr1:5:A:T p.X"), true).unwrap();
        assert_eq!(
            s.values().next().unwrap().kind.as_deref(),
            Some("hypomorph")
        );
        assert!(parse_alleles(Some("odd chr1:5:A:T"), true).is_err());
        assert!(parse_alleles(Some("chr1:5:A"), false).is_err());
        assert!(parse_alleles(None, true).unwrap().is_empty());
    }
}
