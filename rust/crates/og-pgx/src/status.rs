//! How a gene's phenotype reads: normal, a changed function, or neither.

use pycompat::text::{lower, split, strip};
use std::fmt;

/// What a gene's call means for the report, checked in this order: not
/// called, ambiguous, normal, non-normal (listed with its drugs), and
/// unclassified (called, with no phenotype that says normal or changed
/// function: not counted, no drug guidance).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    NotCalled,
    Ambiguous,
    Normal,
    NonNormal,
    Unclassified,
}

impl Status {
    /// The word the phenotypes table writes.
    pub fn as_str(self) -> &'static str {
        match self {
            Status::NotCalled => "not called",
            Status::Ambiguous => "ambiguous",
            Status::Normal => "normal",
            Status::NonNormal => "non-normal",
            Status::Unclassified => "unclassified",
        }
    }
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

const NORMAL_WORDS: [&str; 3] = ["normal", "typical", "extensive"];

/// Normal function, or an HLA result whose every allele test is negative.
pub fn is_normal(phenotype: &str) -> bool {
    let parts: Vec<String> = phenotype
        .split(';')
        .filter(|p| !strip(p, None).is_empty())
        .map(|p| lower(strip(p, None)))
        .collect();
    if parts.is_empty() {
        return false;
    }
    if parts.iter().all(|p| p.ends_with(" negative")) {
        return true;
    }
    parts.iter().all(|p| {
        split(&p.replace('-', " "), None, None)
            .iter()
            .any(|w| NORMAL_WORDS.contains(w))
    })
}

/// `s` without `prefix`, when it starts with one of `prefixes`.
fn after<'a>(s: &'a str, prefixes: &[&str]) -> Option<&'a str> {
    prefixes.iter().find_map(|p| s.strip_prefix(p))
}

/// An optional `likely ` or `possible `.
fn hedge(s: &str) -> &str {
    after(s, &["likely ", "possible "]).unwrap_or(s)
}

/// One `;`-part of a phenotype, lower-cased, hyphens as spaces and its
/// whitespace collapsed, names a changed metabolism or function: a poor,
/// intermediate, rapid or ultrarapid metabolizer; decreased, increased, poor
/// or no function; a positive HLA test; a G6PD deficiency; malignant
/// hyperthermia susceptibility; an increased risk (MT-RNR1); an
/// ivacaftor-responsive CFTR; a VKORC1 -1639 genotype with the A allele.
fn changed_function(p: &str) -> bool {
    let metabolizer = hedge(p).strip_suffix(" metabolizer").is_some_and(|m| {
        matches!(
            m,
            "poor" | "intermediate" | "rapid" | "ultrarapid" | "ultra rapid"
        )
    });
    let function = hedge(p)
        .strip_suffix(" function")
        .is_some_and(|m| matches!(m, "decreased" | "increased" | "poor" | "no"));
    metabolizer
        || function
        || p.strip_suffix(" positive").is_some_and(|s| !s.is_empty())
        || matches!(
            p,
            "deficient"
                | "deficient with cnsha"
                | "variable"
                | "malignant hyperthermia susceptible"
                | "malignant hyperthermia susceptibility"
                | "ivacaftor responsive"
                | "ivacaftor responsive in cf patients"
                | "1639 ga"
                | "1639 aa"
        )
        || p.strip_prefix("increased risk of ")
            .is_some_and(|s| !s.is_empty())
}

/// A phenotype with at least one part that names a changed function
/// (`changed_function`). Anything else PharmCAT writes for a called gene
/// (n/a, `no phenotype assigned`, VKORC1's reference `-1639 GG`, `Uncertain
/// Susceptibility`, `Indeterminate`) names none.
pub fn is_non_normal(phenotype: &str) -> bool {
    phenotype.split(';').any(|p| {
        let p = split(&lower(p).replace('-', " "), None, None).join(" ");
        !p.is_empty() && changed_function(&p)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_normal_phenotypes() {
        for (phen, want) in [
            ("Likely Poor Metabolizer", true),
            ("Possible Intermediate Metabolizer", true),
            ("Ultrarapid Metabolizer", true),
            ("Ultra-rapid Metabolizer", true),
            ("Rapid  Metabolizer", true),
            ("Possible Decreased Function", true),
            ("Increased Function", true),
            ("No Function", true),
            ("Deficient", true),
            ("Deficient with CNSHA", true),
            ("Variable", true),
            ("Malignant Hyperthermia Susceptibility", true),
            ("Malignant Hyperthermia Susceptible", true),
            ("*57:01 positive", true),
            (" positive", false),
            (
                "increased risk of aminoglycoside-induced hearing loss",
                true,
            ),
            ("increased risk of", false),
            ("ivacaftor responsive in CF patients", true),
            ("ivacaftor non-responsive in CF patients", false),
            ("Indeterminate", false),
            ("n/a", false),
            ("no phenotype assigned", false),
            ("-1639 GA", true),
            ("-1639 AA", true),
            ("-1639 GG", false),
            ("Uncertain Susceptibility", false),
            ("*58:01 negative", false),
            ("*58:01 negative; *57:01 positive", true),
            ("Probably Poor Metabolizer", false),
            ("Normal Metabolizer", false),
            ("", false),
        ] {
            assert_eq!(is_non_normal(phen), want, "{phen:?}");
        }
    }

    #[test]
    fn normal_phenotypes() {
        for (phen, want) in [
            ("Normal Metabolizer", true),
            ("Normal Function", true),
            ("Extensive Metabolizer", true),
            ("Typical", true),
            ("*57:01 negative; *58:01 negative", true),
            ("*57:01 negative; *58:01 positive", false),
            ("Normal Metabolizer; Poor Function", false),
            ("normal-function", true),
            ("abnormal", false),
            ("", false),
            (" ; ", false),
        ] {
            assert_eq!(is_normal(phen), want, "{phen:?}");
        }
    }
}
