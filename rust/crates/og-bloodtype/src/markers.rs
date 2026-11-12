//! The inputs that define the systems: the marker table
//! (`assets/bloodtype_markers.tsv`, one row per defining variant) and the
//! RhD depth regions (the RHD row of `bin/cn_loci.tsv`).

use crate::Error;
use crate::rules::{NULL_HOSTS, NULL_SYSTEMS};
use indexmap::IndexMap;
use og_cn::catalog::{Catalog, Locus};
use pycompat::csv::{DictReader, Quoting};
use std::io;

/// One `site` row of the marker table. A column the row lacks reads empty.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Marker {
    pub chrom: String,
    /// 0-based start, as the BED column holds it.
    pub start: i64,
    /// 1-based end: the VCF position of the variant.
    pub pos: i64,
    pub system: String,
    pub name: String,
    pub ref_: String,
    /// The ALT alleles, `,`-joined (a null row may list several).
    pub alt: String,
    pub ref_means: String,
    pub alt_means: String,
    pub hgvs: String,
    pub rsid: String,
    pub caveat: String,
    pub meaning: String,
    pub high_prevalence: String,
    pub null_on: String,
}

impl Marker {
    /// A row whose ALT makes no antigen.
    pub fn is_null(&self) -> bool {
        self.alt_means == "null"
    }
}

/// One RhD depth region: a segment of the RHD catalog row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    /// The segment's label (`RHD exon 3`).
    pub label: String,
}

/// The marker table's sites by name, in table order.
pub type Sites = IndexMap<String, Marker>;
/// The RhD depth regions by name: flanks, gene, exons, as the catalog row
/// lists its controls, targets and details.
pub type Regions = IndexMap<String, Region>;

/// The RHD row of the catalog.
pub fn rhd_row(cat: &Catalog) -> Result<&Locus, Error> {
    cat.loci
        .get("RHD")
        .ok_or_else(|| Error(format!("{}: no RHD row", cat.path)))
}

/// The RhD depth regions of the catalog's RHD row.
pub fn rhd_regions(row: &Locus) -> Regions {
    row.controls
        .iter()
        .chain(&row.targets)
        .chain(&row.details)
        .map(|s| {
            (
                s.name.clone(),
                Region {
                    chrom: row.chrom.clone(),
                    start: s.start,
                    end: s.end,
                    label: s.label.clone(),
                },
            )
        })
        .collect()
}

fn whole(v: &str, col: &str, path: &str) -> Result<i64, Error> {
    pycompat::num::parse_int(v)
        .ok()
        .and_then(|n| i64::try_from(n).ok())
        .ok_or_else(|| {
            Error(format!(
                "{path}: {col} {} is not a whole number",
                og_filters::py_repr(v)
            ))
        })
}

/// The site rows of the marker table at `path`, by name, in table order.
/// `#` lines and blank lines are skipped; the first other line is the
/// header, which must have a `kind` column. A row of another kind, a row
/// with an empty REF or ALT, or a null row in a system that reads no null
/// rows, is refused.
pub fn read_markers(path: &str) -> Result<Sites, Error> {
    let text = std::fs::read(path).map_err(|e| Error::io(&e, path))?;
    let mut lines = og_io::Lines::new(&text[..], og_io::Utf8::Strict);
    let mut kept: Vec<String> = Vec::new();
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        if !line.starts_with('#') && !pycompat::text::strip(&line, None).is_empty() {
            kept.push(line.clone());
        }
    }
    let mut rows = DictReader::new(kept.into_iter().map(io::Result::Ok), '\t', Quoting::Minimal);
    let no_rows = || {
        Error(format!(
            "{path}: no marker rows (header with a 'kind' column expected)"
        ))
    };
    let has_kind = rows
        .header()
        .map_err(|e| Error(format!("{path}: {e}")))?
        .is_some_and(|h| h.names().iter().any(|n| n == "kind"));
    let mut sites = Sites::new();
    let mut any = false;
    for row in rows.by_ref() {
        let row = row.map_err(|e| Error(format!("{path}: {e}")))?;
        if !has_kind {
            return Err(no_rows());
        }
        any = true;
        let get = |k: &str| row.get(k).flatten().unwrap_or_default().to_owned();
        let start = whole(&get("start"), "start", path)?;
        let end = whole(&get("end"), "end", path)?;
        let (name, kind) = (get("name"), get("kind"));
        if kind != "site" {
            return Err(Error(format!(
                "{path}: row {name} has kind '{kind}' (site; the RhD depth regions are the RHD row of bin/cn_loci.tsv)"
            )));
        }
        let (ref_, alt) = (get("ref"), get("alt"));
        if ref_.is_empty() || alt.split(',').any(str::is_empty) {
            return Err(Error(format!(
                "{path}: row {name} has an empty ref or alt: every site needs both"
            )));
        }
        let m = Marker {
            chrom: get("chrom"),
            start,
            pos: end,
            system: get("system"),
            name: name.clone(),
            ref_,
            alt,
            ref_means: get("ref_means"),
            alt_means: get("alt_means"),
            hgvs: get("hgvs"),
            rsid: get("rsid"),
            caveat: get("caveat"),
            meaning: get("meaning"),
            high_prevalence: get("high_prevalence"),
            null_on: get("null_on"),
        };
        if m.is_null()
            && !NULL_SYSTEMS.iter().any(|s| s.system == m.system)
            && !NULL_HOSTS.contains(&m.system.as_str())
        {
            return Err(Error(format!(
                "{path}: null row {name} is in {}, which reads no null rows",
                m.system
            )));
        }
        sites.insert(name, m);
    }
    if !any {
        return Err(no_rows());
    }
    Ok(sites)
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEAD: &str = "chrom\tstart\tend\tkind\tsystem\tname\tref\talt\tref_means\talt_means\n";

    fn table(name: &str, body: &str) -> String {
        let p = std::env::temp_dir().join(format!(
            "og-bloodtype-markers-{}-{name}",
            std::process::id()
        ));
        std::fs::write(&p, format!("# a comment\n{HEAD}{body}")).unwrap();
        p.to_string_lossy().into_owned()
    }

    #[test]
    fn a_site_row_reads() {
        let p = table("ok", "chr7\t10\t11\tsite\tKell\tKEL_578\tG\tA\tk\tK\n");
        let s = read_markers(&p).unwrap();
        assert_eq!((s["KEL_578"].pos, s["KEL_578"].alt.as_str()), (11, "A"));
    }

    #[test]
    fn a_row_with_an_empty_ref_or_alt_is_refused() {
        for (name, body) in [
            ("noalt", "chr7\t10\t11\tsite\tKell\tKEL_578\tG\t\tk\tK\n"),
            ("noref", "chr7\t10\t11\tsite\tKell\tKEL_578\t\tA\tk\tK\n"),
            (
                "emptyalt",
                "chr7\t10\t11\tsite\tKell\tKEL_578\tG\tA,\tk\tK\n",
            ),
        ] {
            let e = read_markers(&table(name, body)).unwrap_err();
            assert!(
                e.0.ends_with("row KEL_578 has an empty ref or alt: every site needs both"),
                "{e}"
            );
        }
    }
}
