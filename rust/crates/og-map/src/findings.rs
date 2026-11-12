//! The findings layer: each finding of the results body with a place on
//! the map.
//!
//! - A finding's `viewer` (`chrom:pos` or `chrom:start-end`, 1-based, the
//!   page's reads link) places it: a position is `[pos - 1, pos)`, a span
//!   `[start - 1, end)`. The `chr` prefix may be left out.
//! - Without a usable `viewer`, a finding that names a gene is placed on
//!   that gene's span, read from the GENCODE genes file the gene layer
//!   reads (`genes::parse_gtf`): medicine genes, HLA, LPA, a gene-level
//!   result. A name on several chromosomes takes the first in map order.
//! - A finding with neither is not on the map. The layer note counts those
//!   that no closed result hides; a closable one, and one the results'
//!   veil ties to a closable result (`Finding::veiled`), is left out of
//!   every count.
//! - The layer's count and totals take the rows with no veil tag only: a
//!   closable finding, and an open one over a closable one's span, are in
//!   no number.
//! - `zyg` comes from the finding's own `Genotype` evidence, else `unknown`:
//!   no VCF is read for a label.

use crate::facts::{self, Finding};
use crate::genes;
use crate::reference::{CHROMS, Reference, chrom_index};
use crate::veil::{self, Veil};
use crate::{Layer, State, text};
use std::collections::{HashMap, HashSet};
use std::path::Path;

/// One finding on the map.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    /// Its index in the page's order of the findings.
    pub order: usize,
    pub finding: Finding,
    /// Index into [`CHROMS`].
    pub chrom: usize,
    /// The 1-based position of a finding placed at one base, else `None`.
    pub pos: Option<u64>,
    /// 0-based start.
    pub s: u64,
    /// Exclusive end.
    pub e: u64,
    /// `het`, `hom`, `hemi` or `unknown`.
    pub zyg: &'static str,
}

impl Row {
    /// The place as the reads viewer takes it: 1-based and inclusive.
    pub fn reads(&self) -> String {
        match self.pos {
            Some(p) => format!("{}:{p}", CHROMS[self.chrom]),
            None => format!("{}:{}-{}", CHROMS[self.chrom], self.s + 1, self.e),
        }
    }

    /// The finding's page, relative to the person's page.
    pub fn href(&self) -> String {
        format!("finding/{}/", self.finding.id)
    }

    /// The sentence the map shows for it.
    pub fn phrase(&self) -> String {
        let f = &self.finding;
        match (&f.gate, &f.gene) {
            (Some(_), Some(g)) => text::gated_phrase(g, CHROMS[self.chrom]),
            _ => text::finding_phrase(&f.level, &f.title, CHROMS[self.chrom]),
        }
    }
}

impl veil::Span for Row {
    fn chrom(&self) -> usize {
        self.chrom
    }
    fn span(&self) -> (u64, u64) {
        (self.s, self.e)
    }
}

/// The layer, its rows sorted by chromosome and start, and the veil of
/// its closable findings.
#[derive(Debug, Clone, PartialEq)]
pub struct Placed {
    pub layer: Layer,
    pub rows: Vec<Row>,
    pub veil: Veil,
}

/// `chrom:pos` or `chrom:start-end` (1-based, commas and spaces ignored)
/// as (chromosome, 0-based start, end, the position of a single base).
pub fn locus(v: &str) -> Option<(usize, u64, u64, Option<u64>)> {
    let v: String = v.chars().filter(|c| !matches!(c, ',' | ' ')).collect();
    let (chrom, at) = v.split_once(':')?;
    let chrom = if chrom.starts_with("chr") {
        chrom_index(chrom)
    } else {
        chrom_index(&format!("chr{chrom}"))
    }?;
    let num = |x: &str| x.parse::<u64>().ok().filter(|&n| n >= 1);
    match at.split_once('-') {
        Some((a, b)) => {
            let (a, b) = (num(a)?, num(b)?);
            (b >= a).then_some((chrom, a - 1, b, None))
        }
        None => {
            let p = num(at)?;
            Some((chrom, p - 1, p, Some(p)))
        }
    }
}

/// The zygosity word of a genotype as the results body writes it.
pub fn zygosity(genotype: Option<&str>) -> &'static str {
    let g = genotype.unwrap_or("").trim().to_ascii_lowercase();
    match g.as_str() {
        "het" | "heterozygous" | "0/1" | "0|1" | "1|0" | "1/0" => "het",
        "hom" | "homozygous" | "1/1" | "1|1" => "hom",
        "hemi" | "hemizygous" => "hemi",
        _ => "unknown",
    }
}

/// The span of each gene in `names`, from the highest installed GENCODE
/// release under `root`. Empty when it is not installed or cannot be read:
/// those findings then stay off the map, counted in the note.
fn gene_spans(root: &Path, names: &HashSet<&str>) -> HashMap<String, (usize, u64, u64)> {
    let mut out = HashMap::new();
    if names.is_empty() {
        return out;
    }
    let Ok(inputs) = genes::locate(root) else {
        return out;
    };
    let name = inputs
        .gtf
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let Ok(input) = og_io::open(&inputs.gtf) else {
        return out;
    };
    let Ok(parsed) = genes::parse_gtf(input, &name) else {
        return out;
    };
    // Sorted by chromosome then start: the first one kept is first in map order.
    for g in parsed {
        if names.contains(g.name.as_str()) {
            out.entry(g.name).or_insert((g.chrom, g.s, g.e));
        }
    }
    out
}

/// Places the findings of the results body on the map.
pub fn read(findings: &facts::Findings, root: &Path, run: &str, reference: &Reference) -> Placed {
    let mut layer = Layer::new("24", "Findings");
    layer.caveat = text::FINDINGS_CAVEAT.to_owned();
    layer.source = "{sample}/summary.json".replace("{sample}", run);
    let list = match findings {
        facts::Findings::Missing => {
            layer.note = Some(text::FINDINGS_MISSING.to_owned());
            return Placed {
                layer,
                rows: Vec::new(),
                veil: Veil::default(),
            };
        }
        facts::Findings::Unreadable(why) => {
            layer.state = State::Unreadable;
            layer.note = Some(why.clone());
            return Placed {
                layer,
                rows: Vec::new(),
                veil: Veil::default(),
            };
        }
        facts::Findings::Read(list) => list,
    };
    let within = |(c, s, e, p): (usize, u64, u64, Option<u64>)| {
        let len = reference.chroms[c].len;
        (s < len).then_some((c, s, e.min(len), p))
    };
    let at_viewer: Vec<_> = list
        .iter()
        .map(|f| f.viewer.as_deref().and_then(locus).and_then(within))
        .collect();
    let wanted: HashSet<&str> = list
        .iter()
        .zip(&at_viewer)
        .filter(|(_, at)| at.is_none())
        .filter_map(|(f, _)| f.gene.as_deref())
        .collect();
    let spans = gene_spans(root, &wanted);
    let mut rows = Vec::with_capacity(list.len());
    let mut open_unplaced = 0usize;
    for (order, (f, at)) in list.iter().zip(at_viewer).enumerate() {
        let place = at.or_else(|| {
            let &(c, s, e) = spans.get(f.gene.as_deref()?)?;
            within((c, s, e, None))
        });
        let Some((chrom, s, e, pos)) = place else {
            if !veil::closable(f) && !f.veiled {
                open_unplaced += 1;
            }
            continue;
        };
        rows.push(Row {
            order,
            finding: f.clone(),
            chrom,
            pos,
            s,
            e,
            zyg: zygosity(f.genotype.as_deref()),
        });
    }
    rows.sort_by_key(|r| (r.chrom, r.s, r.e, r.order));
    // Gated findings first, then the page's order: a row over both a gated
    // finding and an act finding carries the gated one, which every page
    // keeps closed.
    let veil = Veil::new(rows.iter().filter(|r| veil::closable(&r.finding)).map(|r| {
        let key = (r.finding.gate.is_none(), r.order);
        (key, r.finding.id.as_str(), r.chrom, (r.s, r.e))
    }));
    // A row of an open finding over a closable one's span is tagged too, so
    // the counts take the untagged rows only.
    let kept = veil::untagged(&veil.tags(&rows));
    layer.state = State::Ok;
    layer.count = kept;
    layer.totals = Some((kept + open_unplaced, kept));
    layer.empty = Some(text::FINDINGS_EMPTY.to_owned());
    if open_unplaced > 0 {
        layer.note = Some(text::findings_unplaced(open_unplaced));
    }
    Placed { layer, rows, veil }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loci() {
        let c = |n: &str| chrom_index(n).unwrap();
        assert_eq!(
            locus("chr13:32,340,000"),
            Some((c("chr13"), 32_339_999, 32_340_000, Some(32_340_000)))
        );
        assert_eq!(locus("17:100-200"), Some((c("chr17"), 99, 200, None)));
        assert_eq!(locus("chrX:5-5"), Some((c("chrX"), 4, 5, None)));
        for bad in [
            "chrM:5",
            "chr1:0",
            "chr1:9-3",
            "chr1",
            "chr1:x",
            "chr22_KI270879v1_alt:100-900",
        ] {
            assert_eq!(locus(bad), None, "{bad}");
        }
    }

    #[test]
    fn zygosity_words() {
        assert_eq!(zygosity(Some("0/1")), "het");
        assert_eq!(zygosity(Some("HOM")), "hom");
        assert_eq!(zygosity(Some("hemizygous")), "hemi");
        assert_eq!(zygosity(Some("./.")), "unknown");
        assert_eq!(zygosity(None), "unknown");
    }
}
