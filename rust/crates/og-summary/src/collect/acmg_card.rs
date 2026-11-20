//! The ACMG secondary-findings card as its own summary section,
//! `acmg_sf`, whose state follows step 06: ok when step 06's ClinVar hits
//! VCF is there, not run (or failed) with step 06's reason otherwise.
//!
//! The card's main input needs no VEP: step 06's P/LP hits, step 05's
//! class 4 and 5 SVs and step 35's deletions, with each gene's rule from
//! `assets/acmg_sf.tsv`. Only the rare HIGH-impact source needs step 23's
//! clinical VCF, so a run without the VEP cache still gets the card. The
//! verdicts are [`acmg::acmg_sf`]'s, the one function the clinical section
//! also calls: with the clinical VCF when it is there and readable, without
//! it otherwise. Without it the card says `high_impact_checked: false` and
//! `high_impact_note` names why and what that can miss. With a clinical
//! VCF left by an earlier run (step 23 not ok in this run, the VCF older
//! than the run's start) the card reads it and its `high_impact_note` says
//! the HIGH-impact source is from that earlier run.

use super::acmg::{self, ClinicalVcf};
use super::clinvar::HitsVcf;
use super::two_hits::DeletionInputs;
use pycompat::json::Value;

/// Why step 23's rare HIGH-impact records are not in the card.
#[derive(Debug, Clone, PartialEq)]
pub enum NoHighImpact {
    /// The run status gives step 23's result (`failed (exit 1)`,
    /// `skipped (needs VEP, ...)`, or `ok` with no clinical VCF left).
    Step23(String),
    /// No run status names step 23 (a Nextflow report, a folder without
    /// `logs/run_status.tsv`).
    NotRun,
    /// The clinical VCF is there and cannot be read.
    Unreadable(String),
}

impl NoHighImpact {
    /// `high_impact_note`: why, then what the card can miss.
    pub fn note(&self) -> String {
        let why = match self {
            NoHighImpact::Step23(r) if r.starts_with("failed") => {
                format!("step 23 (clinical filter) {r}")
            }
            NoHighImpact::Step23(r) if r.starts_with("skipped") => {
                format!("step 23 (clinical filter) was {r}")
            }
            NoHighImpact::Step23(r) => {
                format!("step 23 (clinical filter) was {r} but wrote no clinical VCF")
            }
            NoHighImpact::NotRun => {
                "step 23 (clinical filter) did not run; it needs VEP's annotation (step 13) and \
                 the VEP cache"
                    .to_owned()
            }
            NoHighImpact::Unreadable(e) => format!("step 23's clinical VCF cannot be read: {e}"),
        };
        format!(
            "Rare HIGH-impact variants were not screened: {why}. This card holds the ClinVar P/LP \
             hits and the class 4 and 5 structural variants only, so it can miss a predicted \
             loss-of-function variant that is not in ClinVar, alone or as the second hit of a \
             recessive gene."
        )
    }
}

/// What the card reads: the step 06 hits, the clinical VCF when step 23
/// left one, the deletions of the two-hit rule, step 23's result word
/// from the run status, and the clinical VCF's file date when it is from
/// an earlier run.
pub struct Inputs<'a> {
    pub assets: &'a str,
    pub hits: &'a HitsVcf,
    pub clinical: Option<&'a Result<ClinicalVcf, String>>,
    pub deletions: &'a DeletionInputs,
    pub step23: Option<&'a str>,
    /// The clinical VCF's file date (UTC minute) when step 23 was not ok in
    /// this run and the VCF was written before the run started.
    pub earlier_vcf: Option<String>,
}

/// `high_impact_note` for a card that read a clinical VCF of an earlier run.
pub fn earlier_note(date: &str, step23: Option<&str>) -> String {
    let result = step23.filter(|r| !r.is_empty()).unwrap_or("not run");
    format!(
        "Rare HIGH-impact variants come from the clinical VCF of an earlier run ({date}); step \
         23 (clinical filter) was {result} in this run."
    )
}

/// The data of the `acmg_sf` section.
#[derive(Debug, Clone, PartialEq)]
pub struct AcmgCard(pub Value);

impl AcmgCard {
    /// The card. An `Err` is the section's note (a hits record whose POS
    /// is no integer); a rules list or a hits VCF that cannot be read gives
    /// the card's `unavailable` form, as in the clinical section.
    pub fn read(i: &Inputs) -> Result<AcmgCard, String> {
        // A clinical VCF that cannot be read costs the HIGH-impact source
        // only, never the ClinVar rows.
        let (clinical, missing) = match i.clinical {
            Some(Ok(c)) => match &c.high {
                Ok(_) => (i.clinical, None),
                Err(e) => (None, Some(NoHighImpact::Unreadable(e.clone()))),
            },
            Some(Err(e)) => (None, Some(NoHighImpact::Unreadable(e.clone()))),
            None => (
                None,
                Some(
                    i.step23
                        .map_or(NoHighImpact::NotRun, |r| NoHighImpact::Step23(r.to_owned())),
                ),
            ),
        };
        let mut v = acmg::acmg_sf(&acmg::Inputs {
            assets: i.assets,
            hits: Some(&i.hits.read),
            clinical,
            deletions: i.deletions,
        })?;
        if let Value::Object(m) = &mut v {
            let checked = m.get("high_impact_checked").cloned();
            match (missing, &i.earlier_vcf) {
                (Some(why), _) if checked == Some(Value::Bool(false)) => {
                    m.insert("high_impact_note".into(), Value::Str(why.note()));
                }
                (None, Some(date)) if checked == Some(Value::Bool(true)) => {
                    m.insert(
                        "high_impact_note".into(),
                        Value::Str(earlier_note(date, i.step23)),
                    );
                }
                _ => {}
            }
        }
        Ok(AcmgCard(v))
    }

    pub fn to_json(&self) -> Value {
        self.0.clone()
    }
}
