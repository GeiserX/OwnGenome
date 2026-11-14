//! What the map draws one by one at each level: the rule `drawsAt(layer, row, level)` of
//! `assets/live/map.js`, written once more here so a test can hold the two to the same
//! answers. The page's canvas, its whole-genome strip and its "In view" table all ask
//! that one JS function; `og map build --predicate-table` writes this module's answers
//! for every row of a run's document, and the browser check compares them row by row.
//!
//! A row that is not drawn one by one at a level is either drawn as a density (the
//! shorter runs of homozygosity at the genome level, as 1 Mb bins) or not drawn at all
//! (genes above the region level). The veil is a separate rule: the page drops a closed
//! row whatever this module says.

use crate::genes::{self, Genes};
use crate::reference::CHROMS;
use crate::{BuildError, Model, Request};
use std::fmt::Write as _;

/// The map's zoom levels, from the farthest out. `assets/live/map.js` `levelOf()` picks
/// one from the camera: the genome while no chromosome is open, then by bp per CSS pixel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Level {
    Genome,
    Chromosome,
    Region,
    Gene,
    Reads,
}

impl Level {
    /// Every level, in order.
    pub const ALL: [Level; 5] = [
        Level::Genome,
        Level::Chromosome,
        Level::Region,
        Level::Gene,
        Level::Reads,
    ];

    /// The level's name, as `map.js` writes it.
    pub fn word(self) -> &'static str {
        match self {
            Level::Genome => "genome",
            Level::Chromosome => "chromosome",
            Level::Region => "region",
            Level::Gene => "gene",
            Level::Reads => "reads",
        }
    }
}

/// One row of a layer, with only what the rule reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Feature {
    /// A finding with a place.
    Finding,
    /// A run of homozygosity; `long` is the document's `long` column (5 Mb or more,
    /// not an artefact).
    Roh { long: bool },
    /// A gene of the GENCODE files.
    Gene,
}

impl Feature {
    /// The layer's name, as `map.js` writes it.
    pub fn layer(self) -> &'static str {
        match self {
            Feature::Finding => "findings",
            Feature::Roh { .. } => "roh",
            Feature::Gene => "genes",
        }
    }
}

/// Whether the map draws this row one by one at this level.
///
/// - A finding: at every level (a glyph on its row at the genome level, a mark in its
///   lane once its chromosome is open).
/// - A run of homozygosity: a long one at every level; a shorter one only once its
///   chromosome is open (at the genome level it is part of the 1 Mb coverage bins).
/// - A gene: from the region level down.
pub fn draws_at(f: Feature, level: Level) -> bool {
    match f {
        Feature::Finding => true,
        Feature::Roh { long } => long || level > Level::Genome,
        Feature::Gene => level >= Level::Region,
    }
}

/// The header of [`table`].
pub const HEADER: &str = "layer\tchrom\ts\te\tid\tgenome\tchromosome\tregion\tgene\treads\n";

/// One line of [`table`]: where the row is and its answer at each level.
fn line(out: &mut String, f: Feature, chrom: usize, s: u64, e: u64, id: &str) {
    let _ = write!(out, "{}\t{}\t{s}\t{e}\t{id}", f.layer(), CHROMS[chrom]);
    for lv in Level::ALL {
        out.push_str(if draws_at(f, lv) { "\t1" } else { "\t0" });
    }
    out.push('\n');
}

/// The answers for every row of the model's layers (findings by id, runs of
/// homozygosity) and, when given, every gene: a header, then one tab-separated line
/// per row in the document's order, `1` where the row draws one by one and `0` where
/// it does not. Every row is listed, a closable finding's too: the veil is not this
/// rule.
pub fn table(m: &Model, genes: Option<&Genes>) -> String {
    let mut out = String::from(HEADER);
    for r in &m.findings.rows {
        line(&mut out, Feature::Finding, r.chrom, r.s, r.e, &r.finding.id);
    }
    for r in &m.roh.rows {
        line(
            &mut out,
            Feature::Roh { long: r.long },
            r.chrom,
            r.s,
            r.e,
            "",
        );
    }
    for g in genes.map_or(&[][..], |g| &g.genes) {
        line(&mut out, Feature::Gene, g.chrom, g.s, g.e, &g.name);
    }
    out
}

/// `og map build --predicate-table`: [`table`] for the request's run, with the genes of
/// the highest GENCODE release installed under the root (none while none is). Built
/// fresh each time: the document cache is never read or written.
pub fn build_table(req: &Request<'_>) -> Result<String, BuildError> {
    let model = crate::build_person(req)?;
    let genes = match genes::locate(req.root) {
        Ok(inputs) => Some(genes::load(&inputs)?),
        Err(_) => None,
    };
    Ok(table(&model, genes.as_ref()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn row(f: Feature) -> Vec<bool> {
        Level::ALL.iter().map(|&lv| draws_at(f, lv)).collect()
    }

    #[test]
    fn each_layer_at_each_level() {
        assert_eq!(row(Feature::Finding), [true; 5]);
        assert_eq!(row(Feature::Roh { long: true }), [true; 5]);
        assert_eq!(
            row(Feature::Roh { long: false }),
            [false, true, true, true, true]
        );
        assert_eq!(row(Feature::Gene), [false, false, true, true, true]);
    }

    #[test]
    fn levels_in_order_with_their_words() {
        let words: Vec<_> = Level::ALL.iter().map(|l| l.word()).collect();
        assert_eq!(words, ["genome", "chromosome", "region", "gene", "reads"]);
        assert!(Level::ALL.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(
            HEADER.trim_end().split('\t').skip(5).collect::<Vec<_>>(),
            words
        );
    }

    #[test]
    fn the_table_of_an_invented_run() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR"));
        let req = Request {
            root: &dir.join("tests/data/genome"),
            run: "DEMO-001",
            assets: &dir.join("../../../assets"),
            made_by: "test",
            stdin: b"",
        };
        let t = build_table(&req).unwrap();
        let lines: Vec<&str> = t.lines().collect();
        assert_eq!(format!("{}\n", lines[0]), HEADER);
        let roh: Vec<Vec<&str>> = lines[1..]
            .iter()
            .map(|l| l.split('\t').collect())
            .filter(|c: &Vec<&str>| c[0] == "roh")
            .collect();
        assert_eq!(
            roh.len() + 1,
            lines.len(),
            "only runs: no findings, no genes"
        );
        let m = crate::build_person(&req).unwrap();
        assert_eq!(roh.len(), m.roh.rows.len());
        for (c, r) in roh.iter().zip(&m.roh.rows) {
            assert_eq!(c[1], CHROMS[r.chrom]);
            assert_eq!((c[2], c[3]), (&*r.s.to_string(), &*r.e.to_string()));
            assert_eq!(c[5], if r.long { "1" } else { "0" }, "{c:?}");
            assert_eq!(&c[6..], ["1"; 4]);
        }
        assert!(roh.iter().any(|c| c[5] == "1") && roh.iter().any(|c| c[5] == "0"));
    }

    #[test]
    fn a_line_per_row() {
        let mut out = String::new();
        line(&mut out, Feature::Roh { long: false }, 10, 4, 9, "");
        line(&mut out, Feature::Gene, 23, 0, 1, "SRY");
        assert_eq!(
            out,
            "roh\tchr11\t4\t9\t\t0\t1\t1\t1\t1\ngenes\tchrY\t0\t1\tSRY\t0\t0\t1\t1\t1\n"
        );
    }
}
