//! Step 06's ClinVar screen: what OwnGenome itself decides about a hit.
//!
//! - [`ReviewStatus`]: ClinVar's review status and the one table that
//!   turns it into stars (4 practice guideline, 3 expert panel, 2 multiple
//!   submitters with no conflict, 1 single submitter or conflicting
//!   classifications, 0 anything else).
//! - [`zygosity`], [`combined_note`] and [`CallFloor`]: how a hit's genotype,
//!   its combined-allele tags and its call quality are worded.
//! - [`hap`]: the haplotype signature that says two records can be one
//!   allele.
//! - [`join`]: the combined allele of two neighbouring records
//!   (`og clinvar join-adjacent`).
//! - [`hits`]: one row per hit (`og clinvar hits`).
//! - [`sites`]: the always-report list (`og clinvar deep-intronic`).

pub mod cli;
pub mod hap;
pub mod hits;
pub mod join;
pub mod sites;
pub mod text;

/// ClinVar's review status of a record (INFO/CLNREVSTAT).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReviewStatus {
    PracticeGuideline,
    ExpertPanel,
    MultipleNoConflict,
    SingleSubmitter,
    /// Conflicting classifications (ClinVar's name since 2024) or
    /// interpretations (the older name).
    Conflicting,
    /// No assertion criteria, no classification, or a status this table
    /// does not know.
    Other,
}

impl ReviewStatus {
    /// Reads CLNREVSTAT as ClinVar writes it (`_` for spaces).
    pub fn parse(raw: &[u8]) -> ReviewStatus {
        match raw {
            b"practice_guideline" => ReviewStatus::PracticeGuideline,
            b"reviewed_by_expert_panel" => ReviewStatus::ExpertPanel,
            b"criteria_provided,_multiple_submitters,_no_conflicts" => {
                ReviewStatus::MultipleNoConflict
            }
            b"criteria_provided,_single_submitter" => ReviewStatus::SingleSubmitter,
            b"criteria_provided,_conflicting_classifications"
            | b"criteria_provided,_conflicting_interpretations" => ReviewStatus::Conflicting,
            _ => ReviewStatus::Other,
        }
    }

    /// ClinVar's review stars.
    pub fn stars(self) -> u8 {
        match self {
            ReviewStatus::PracticeGuideline => 4,
            ReviewStatus::ExpertPanel => 3,
            ReviewStatus::MultipleNoConflict => 2,
            ReviewStatus::SingleSubmitter | ReviewStatus::Conflicting => 1,
            ReviewStatus::Other => 0,
        }
    }
}

/// A hit's GT in words: `het`, `hom`, `one copy, other copy not called ref
/// (GT)` for one alt allele beside a missing one (the half of a split
/// multiallelic call, never het, which would say the other copy is the
/// reference), or the GT as written.
pub fn zygosity(gt: &[u8]) -> Vec<u8> {
    match gt {
        b"0/1" | b"1/0" | b"0|1" | b"1|0" => b"het".to_vec(),
        b"1/1" | b"1|1" => b"hom".to_vec(),
        _ if half_call(gt) => [b"one copy, other copy not called ref (", gt, b")"].concat(),
        _ => gt.to_vec(),
    }
}

/// One alt allele (1 to 9) and one missing allele: `1/.`, `./1`, `1|.`.
pub fn half_call(gt: &[u8]) -> bool {
    let alt = |c: u8| (b'1'..=b'9').contains(&c);
    let sep = |c: u8| c == b'/' || c == b'|';
    matches!(gt, &[a, s, b'.'] if alt(a) && sep(s))
        || matches!(gt, &[b'.', s, b] if alt(b) && sep(s))
}

/// The note on a hit from step 06's combined-allele tags: a combined
/// record, a part of one, or one that may be one allele with another hit;
/// `.` when none. On an unphased het (`0/1`) the first two say the parts
/// may sit on different haplotypes.
pub fn combined_note(joined: &[u8], part_of: &[u8], one: &[u8], gt: &[u8]) -> Vec<u8> {
    let phase: &[u8] = match gt.iter().position(|&c| c == b'/') {
        Some(i)
            if text::is_digits(&gt[..i])
                && text::is_digits(&gt[i + 1..])
                && text::cmp_digits(&gt[..i], &gt[i + 1..]).is_ne() =>
        {
            b"; unphased, so the parts may sit on different haplotypes"
        }
        _ => b"",
    };
    if !joined.is_empty() {
        [b"combined allele of ", &variants(joined, b" + ")[..], phase].concat()
    } else if !part_of.is_empty() {
        [
            b"part of a combined allele: ",
            &variants(part_of, b", ")[..],
            phase,
        ]
        .concat()
    } else if !one.is_empty() {
        [
            b"possibly one allele with ",
            &variants(one, b", ")[..],
            b": check ClinVar for a combined record",
        ]
        .concat()
    } else {
        b".".to_vec()
    }
}

/// `chr1:10:C:T|chr1:15:G:A` as `chr1:10 C>T<sep>chr1:15 G>A`.
fn variants(v: &[u8], sep: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, x) in text::split(v, b'|').into_iter().enumerate() {
        if i > 0 {
            out.extend_from_slice(sep);
        }
        let f = text::split(x, b':');
        let g = |k| text::nth(&f, k);
        out.extend([g(1), b":", g(2), b" ", g(3), b">", g(4)].concat());
    }
    out
}

/// Step 06's call floor: a call with a GQ or DP below it is low confidence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CallFloor {
    pub gq: f64,
    pub dp: f64,
}

impl Default for CallFloor {
    fn default() -> Self {
        CallFloor { gq: 30.0, dp: 10.0 }
    }
}

impl CallFloor {
    /// Takes the `GQ=` and `DP=` of a `##PGP_CALL_FLOOR=` header line (the
    /// first of each followed by digits); a value the line lacks stays.
    pub fn read_line(&mut self, line: &[u8]) {
        if let Some(v) = first_number(line, b"GQ=") {
            self.gq = v;
        }
        if let Some(v) = first_number(line, b"DP=") {
            self.dp = v;
        }
    }

    /// True when the GQ or the DP (digits, when the call has one) is below
    /// the floor.
    pub fn low(&self, gq: Option<&[u8]>, dp: Option<&[u8]>) -> bool {
        let below =
            |v: Option<&[u8]>, floor: f64| v.and_then(text::number).is_some_and(|x| x < floor);
        below(gq, self.gq) || below(dp, self.dp)
    }
}

/// The digits after the first `key` that has any.
fn first_number(line: &[u8], key: &[u8]) -> Option<f64> {
    memchr::memmem::find_iter(line, key).find_map(|at| {
        let rest = &line[at + key.len()..];
        let n = rest.iter().take_while(|c| c.is_ascii_digit()).count();
        (n > 0).then(|| text::number(&rest[..n])).flatten()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_review_status_has_its_stars() {
        let table: [(&[u8], u8); 9] = [
            (b"practice_guideline", 4),
            (b"reviewed_by_expert_panel", 3),
            (b"criteria_provided,_multiple_submitters,_no_conflicts", 2),
            (b"criteria_provided,_single_submitter", 1),
            (b"criteria_provided,_conflicting_classifications", 1),
            (b"criteria_provided,_conflicting_interpretations", 1),
            (b"no_assertion_criteria_provided", 0),
            (b"no_classification_provided", 0),
            (b"", 0),
        ];
        for (raw, stars) in table {
            assert_eq!(
                ReviewStatus::parse(raw).stars(),
                stars,
                "{}",
                String::from_utf8_lossy(raw)
            );
        }
    }

    #[test]
    fn genotypes_in_words() {
        assert_eq!(zygosity(b"1|0"), b"het");
        assert_eq!(zygosity(b"1/1"), b"hom");
        assert_eq!(
            zygosity(b"./1"),
            b"one copy, other copy not called ref (./1)"
        );
        assert_eq!(
            zygosity(b"1|."),
            b"one copy, other copy not called ref (1|.)"
        );
        assert_eq!(zygosity(b"0/."), b"0/.");
        assert_eq!(zygosity(b"1/2"), b"1/2");
    }

    #[test]
    fn notes_say_when_an_unphased_het_may_be_two_haplotypes() {
        let j = b"chr1:10:C:T|chr1:15:G:A";
        assert_eq!(
            combined_note(j, b"", b"", b"0/1"),
            b"combined allele of chr1:10 C>T + chr1:15 G>A; unphased, so the parts may sit on different haplotypes"
        );
        assert_eq!(
            combined_note(j, b"", b"", b"0|1"),
            b"combined allele of chr1:10 C>T + chr1:15 G>A"
        );
        assert_eq!(
            combined_note(b"", j, b"", b"01/1"),
            b"part of a combined allele: chr1:10 C>T, chr1:15 G>A"
        );
        assert_eq!(
            combined_note(b"", b"", b"x", b"0/1"),
            b"possibly one allele with x: >: check ClinVar for a combined record"
        );
        assert_eq!(combined_note(b"", b"", b"", b"0/1"), b".");
    }

    #[test]
    fn the_floor_line_and_low_confidence() {
        let mut f = CallFloor::default();
        f.read_line(b"##PGP_CALL_FLOOR=GQ=x,GQ=20,DP=8,HET_AB=0.2-0.8");
        assert_eq!(f, CallFloor { gq: 20.0, dp: 8.0 });
        assert!(f.low(Some(b"19"), None));
        assert!(!f.low(Some(b"20"), Some(b"8")));
        assert!(f.low(None, Some(b"7")));
        assert!(!f.low(None, None));
    }
}
