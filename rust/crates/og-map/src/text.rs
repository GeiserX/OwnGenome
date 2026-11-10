//! Every sentence the map document carries. The page draws these words and
//! writes none that make a claim of its own.

use pycompat::fmt::fixed;

pub const ROH_CAVEAT: &str = "Runs of homozygosity from bcftools roh on this person's variant calls. \
A run counts as long at 5 Mb or more. A run lying mostly in centromere or heterochromatin bands, \
or on a male's single X or Y, is marked as an artefact and never counts as long.";

pub const ROH_EMPTY: &str = "No run of homozygosity in what step 11 assessed.";

/// The note of a missing ROH layer.
pub fn roh_missing(run: &str) -> String {
    format!("Step 11 has not written {run}_roh.txt for this run.")
}

/// One piece of a sentence; `b` marks the words drawn in bold.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Piece {
    pub t: String,
    pub b: bool,
}

fn runs(n: usize) -> String {
    if n == 1 {
        "1 run".to_owned()
    } else {
        format!("{n} runs")
    }
}

/// The genome-level sentence of the long runs: their count and total length.
pub fn roh_long(n: usize, bp: u64) -> Vec<Piece> {
    if n == 0 {
        return vec![Piece {
            t: "No run of homozygosity of 5 Mb or more.".to_owned(),
            b: false,
        }];
    }
    vec![
        Piece {
            t: runs(n),
            b: true,
        },
        Piece {
            t: format!(
                " of homozygosity of 5 Mb or more, {} Mb in all.",
                fixed(bp as f64 / 1e6, 1)
            ),
            b: false,
        },
    ]
}

/// The sentence of one chromosome with long runs.
pub fn roh_long_chrom(n: usize, bp: u64) -> Vec<Piece> {
    vec![Piece {
        t: format!(
            "{} of homozygosity of 5 Mb or more here, {} Mb in all.",
            runs(n),
            fixed(bp as f64 / 1e6, 1)
        ),
        b: false,
    }]
}

pub const FINDINGS_CAVEAT: &str = "Each finding sits at its own position, or on its gene's span \
when it names a gene but no position. A finding with neither stays on the page only.";

pub const FINDINGS_EMPTY: &str = "No finding has a place on the map.";

/// The note of a findings layer with no results body on stdin.
pub const FINDINGS_MISSING: &str =
    "No results for this run yet: the findings come from its summary (step 24).";

/// The note of the findings with no place, counted over those no closed
/// result hides.
pub fn findings_unplaced(n: usize) -> String {
    if n == 1 {
        "1 finding with no position stays on the page.".to_owned()
    } else {
        format!("{n} findings with no position stay on the page.")
    }
}

/// The count of findings on the map, in its three forms. The page counts
/// the untagged findings plus the tagged ones that are open, picks the form
/// and puts the number in place of `{n}`.
pub fn findings_count() -> [(&'static str, Vec<Piece>); 3] {
    let p = |t: &str, b: bool| Piece { t: t.to_owned(), b };
    [
        ("none", vec![p(FINDINGS_EMPTY, false)]),
        (
            "one",
            vec![p("1 finding", true), p(" has a place on the map.", false)],
        ),
        (
            "many",
            vec![
                p("{n} findings", true),
                p(" have a place on the map.", false),
            ],
        ),
    ]
}

/// The page's words for a level (`LEVEL_WORDS` in bin/live_results.py).
pub fn level_words(level: &str) -> &'static str {
    match level {
        "act" => "Take to a doctor",
        "medicine" => "Tell your prescriber",
        "review" => "Worth a look",
        "carrier" => "For family planning",
        "trust" => "Data quality",
        _ => "For information",
    }
}

/// `chr7` as the words "chromosome 7".
fn chromosome(name: &str) -> String {
    format!("chromosome {}", name.strip_prefix("chr").unwrap_or(name))
}

/// The sentence of one finding on the map.
pub fn finding_phrase(level: &str, title: &str, chrom: &str) -> String {
    let title = title.trim().trim_end_matches('.');
    if title.is_empty() {
        format!("{}, on {}.", level_words(level), chromosome(chrom))
    } else {
        format!("{}: {title}, on {}.", level_words(level), chromosome(chrom))
    }
}

/// The sentence of a gated finding: it names the gene, never the result.
pub fn gated_phrase(gene: &str, chrom: &str) -> String {
    format!(
        "{gene}, on {}: kept behind its warning on every page until opened on purpose.",
        chromosome(chrom)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sentences() {
        assert_eq!(
            roh_long(0, 0)[0].t,
            "No run of homozygosity of 5 Mb or more."
        );
        let s = roh_long(2, 12_345_678);
        assert_eq!((s[0].t.as_str(), s[0].b), ("2 runs", true));
        assert_eq!(s[1].t, " of homozygosity of 5 Mb or more, 12.3 Mb in all.");
        assert_eq!(
            roh_long_chrom(1, 5_050_000)[0].t,
            "1 run of homozygosity of 5 Mb or more here, 5.0 Mb in all."
        );
        assert_eq!(
            finding_phrase("act", "Secondary finding in BRCA2.", "chr13"),
            "Take to a doctor: Secondary finding in BRCA2, on chromosome 13."
        );
        assert_eq!(
            finding_phrase("whatever", "", "chrX"),
            "For information, on chromosome X."
        );
        assert_eq!(
            gated_phrase("APOE", "chr19"),
            "APOE, on chromosome 19: kept behind its warning on every page until opened on purpose."
        );
        assert_eq!(
            findings_unplaced(1),
            "1 finding with no position stays on the page."
        );
        assert_eq!(
            findings_unplaced(3),
            "3 findings with no position stay on the page."
        );
    }
}
