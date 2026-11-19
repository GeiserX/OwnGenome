//! The fixed lists `og facts` reads from `assets/`: the panels the coverage
//! map names (ACMG SF v3.3, the CPIC genes, the carrier genes), OwnHealth's
//! sensitive gene spans, the GRCh38 refget identifiers, the APOE sites, the
//! two LPA sites of the traits panel and the site of each one-site PGx gene.

use std::collections::HashMap;
use std::path::Path;

/// OwnHealth pads every gene span by this many bases on each side
/// (`SPAN_PAD_BP` of its `sensitive.rs`).
pub const SPAN_PAD_BP: u64 = 2_000;

/// A gene's GRCh38 span, 1-based and closed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    pub gene: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Span {
    /// Whether `pos` is inside the span padded as OwnHealth pads it.
    pub fn holds(&self, chrom: &str, pos: u64) -> bool {
        self.chrom == chrom
            && self.start.saturating_sub(SPAN_PAD_BP).max(1) <= pos
            && pos <= self.end + SPAN_PAD_BP
    }
}

/// A gene of the carrier list.
#[derive(Debug, Clone)]
pub struct Carrier {
    pub span: Span,
    pub inheritance: String,
    pub condition: String,
    /// The alleles the severity column sets aside (hypomorph or
    /// pseudodeficiency), as `chrom:pos:ref:alt` with their kind: such an
    /// allele never makes a carrier finding on its own.
    pub set_aside: Vec<(String, String)>,
}

/// One site of the APOE rule.
#[derive(Debug, Clone)]
pub struct Site {
    pub rsid: String,
    pub chrom: String,
    pub pos: u64,
    pub ref_: String,
    pub alt: String,
}

/// Every list, loaded.
#[derive(Debug, Clone, Default)]
pub struct Panels {
    /// ACMG SF v3.3 genes and their inheritance.
    pub acmg: Vec<(String, String)>,
    pub cpic: Vec<String>,
    pub carrier: Vec<Carrier>,
    /// OwnHealth's sensitive genes, with their group.
    pub sensitive: Vec<(Span, String)>,
    /// chromosome -> (refget identifier, length)
    pub refget: HashMap<String, (String, u64)>,
    /// The two APOE sites: rs429358 (e4) and rs7412 (e2).
    pub apoe: Vec<Site>,
    /// The LPA sites of the traits panel (rs10455872, rs3798220).
    pub lpa: Vec<Site>,
    /// The site of each one-site PGx gene, which a PGx fact cites.
    pub pgx_sites: Vec<(String, Site)>,
}

/// The PGx genes PharmCAT calls from one site, GRCh38 (PharmCAT's allele
/// definitions): gene, rsid, chromosome, position, reference, alternate.
const PGX_SITES: [(&str, &str, &str, u64, &str, &str); 3] = [
    ("ABCG2", "rs2231142", "chr4", 88_131_171, "G", "T"),
    ("IFNL3", "rs12979860", "chr19", 39_248_147, "C", "T"),
    ("VKORC1", "rs9923231", "chr16", 31_096_368, "C", "T"),
];

/// The rows of a TSV as maps from header names, `#` lines skipped.
fn rows(path: &Path) -> Result<Vec<HashMap<String, String>>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty());
    let header: Vec<&str> = lines
        .next()
        .ok_or_else(|| format!("{}: no header", path.display()))?
        .split('\t')
        .collect();
    Ok(lines
        .map(|l| {
            header
                .iter()
                .zip(l.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect())
}

fn col<'a>(row: &'a HashMap<String, String>, key: &str, path: &Path) -> Result<&'a str, String> {
    row.get(key)
        .map(String::as_str)
        .ok_or_else(|| format!("{}: a row has no {key}", path.display()))
}

fn span(row: &HashMap<String, String>, path: &Path) -> Result<Span, String> {
    let num = |k: &str| -> Result<u64, String> {
        col(row, k, path)?
            .parse()
            .map_err(|_| format!("{}: {k} is not a number", path.display()))
    };
    Ok(Span {
        gene: col(row, "gene", path)?.to_owned(),
        chrom: col(row, "chrom", path)?.to_owned(),
        start: num("start")?,
        end: num("end")?,
    })
}

impl Panels {
    /// Reads the lists from the assets folder.
    pub fn load(assets: &Path) -> Result<Panels, String> {
        let oh = assets.join("ownhealth");
        let mut p = Panels::default();
        let f = assets.join("acmg_sf.tsv");
        for r in rows(&f)? {
            p.acmg.push((
                col(&r, "gene", &f)?.to_owned(),
                col(&r, "inheritance", &f)?.to_owned(),
            ));
        }
        let f = oh.join("cpic_genes.tsv");
        for r in rows(&f)? {
            p.cpic.push(col(&r, "gene", &f)?.to_owned());
        }
        let f = assets.join("carrier_genes.tsv");
        for r in rows(&f)? {
            // severity: "kind chrom:pos:ref:alt name" entries split by ";", or ".".
            let set_aside = col(&r, "severity", &f)?
                .split(';')
                .filter_map(|e| {
                    let mut w = e.split_whitespace();
                    Some((w.next()?, w.next()?))
                })
                .map(|(kind, key)| (key.to_owned(), kind.to_owned()))
                .collect();
            p.carrier.push(Carrier {
                span: span(&r, &f)?,
                inheritance: col(&r, "inheritance", &f)?.to_owned(),
                condition: col(&r, "condition", &f)?.to_owned(),
                set_aside,
            });
        }
        let f = oh.join("sensitive_genes.tsv");
        for r in rows(&f)? {
            p.sensitive
                .push((span(&r, &f)?, col(&r, "group", &f)?.to_owned()));
        }
        let f = oh.join("grch38_refget.tsv");
        for r in rows(&f)? {
            let len = col(&r, "length", &f)?
                .parse()
                .map_err(|_| format!("{}: length is not a number", f.display()))?;
            p.refget.insert(
                col(&r, "chrom", &f)?.to_owned(),
                (col(&r, "refget", &f)?.to_owned(), len),
            );
        }
        let f = assets.join("genotype_rules.tsv");
        let text = std::fs::read_to_string(&f).map_err(|e| format!("{}: {e}", f.display()))?;
        for l in text.lines() {
            let c: Vec<&str> = l.split('\t').collect();
            if c.len() >= 7 && c[0] == "site" && c[1] == "APOE" {
                p.apoe.push(Site {
                    rsid: c[2].to_owned(),
                    chrom: c[3].to_owned(),
                    pos: c[4]
                        .parse()
                        .map_err(|_| format!("{}: APOE position {}", f.display(), c[4]))?,
                    ref_: c[5].to_owned(),
                    alt: c[6].to_owned(),
                });
            }
        }
        if p.apoe.len() != 2 {
            return Err(format!("{}: want the two APOE sites", f.display()));
        }
        let f = assets.join("traits.tsv");
        for r in rows(&f)? {
            if col(&r, "gene", &f)? == "LPA" {
                p.lpa.push(Site {
                    rsid: col(&r, "rsid", &f)?.to_owned(),
                    chrom: col(&r, "chrom", &f)?.to_owned(),
                    pos: col(&r, "pos", &f)?
                        .parse()
                        .map_err(|_| format!("{}: an LPA position", f.display()))?,
                    ref_: col(&r, "ref", &f)?.to_owned(),
                    alt: col(&r, "alt", &f)?.to_owned(),
                });
            }
        }
        p.pgx_sites = PGX_SITES
            .iter()
            .map(|(gene, rsid, chrom, pos, r, a)| {
                (
                    (*gene).to_owned(),
                    Site {
                        rsid: (*rsid).to_owned(),
                        chrom: (*chrom).to_owned(),
                        pos: *pos,
                        ref_: (*r).to_owned(),
                        alt: (*a).to_owned(),
                    },
                )
            })
            .collect();
        Ok(p)
    }

    /// The site of a one-site PGx gene.
    pub fn pgx_site(&self, gene: &str) -> Option<&Site> {
        self.pgx_sites
            .iter()
            .find(|(g, _)| g == gene)
            .map(|(_, s)| s)
    }

    /// The sensitive group of `gene`, when OwnHealth lists it.
    pub fn sensitive_group(&self, gene: &str) -> Option<&str> {
        self.sensitive
            .iter()
            .find(|(s, _)| s.gene == gene)
            .map(|(_, g)| g.as_str())
    }

    /// Whether a position is inside OwnHealth's sensitive gene spans.
    pub fn in_sensitive_span(&self, chrom: &str, pos: u64) -> bool {
        self.sensitive.iter().any(|(s, _)| s.holds(chrom, pos))
    }

    /// Whether a position is inside the carrier-panel gene spans.
    pub fn in_carrier_span(&self, chrom: &str, pos: u64) -> bool {
        self.carrier.iter().any(|c| c.span.holds(chrom, pos))
    }

    /// The carrier-list entry of `gene`.
    pub fn carrier_gene(&self, gene: &str) -> Option<&Carrier> {
        self.carrier.iter().find(|c| c.span.gene == gene)
    }

    /// Whether `gene` is on ACMG SF v3.3.
    pub fn is_acmg(&self, gene: &str) -> bool {
        self.acmg.iter().any(|(g, _)| g == gene)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spans_are_padded_as_ownhealth_pads_them() {
        let s = Span {
            gene: "G".into(),
            chrom: "chr1".into(),
            start: 10_000,
            end: 20_000,
        };
        assert!(s.holds("chr1", 8_000) && s.holds("chr1", 22_000));
        assert!(!s.holds("chr1", 7_999) && !s.holds("chr1", 22_001));
        assert!(!s.holds("chr2", 15_000));
    }
}
