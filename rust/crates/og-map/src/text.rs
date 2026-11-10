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
    }
}
