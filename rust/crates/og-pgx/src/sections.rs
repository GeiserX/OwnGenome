//! The sections of the recommendations text beside the gene results, each
//! a list of lines.

use crate::guidance::{GuidanceMap, gene_drugs};
use crate::pyio::PyErr;
use crate::report::{GeneCall, STEP36_GENES, position_read, text, truthy_field};
use crate::status::{Status, is_normal};
use crate::tables::{
    ClinvarHit, ConsensusRow, Contraindication, Inhibitors, Medication, PypgxCall, norm_drug,
    pypgx_called,
};
use pycompat::fmt::general;
use pycompat::json::Value;
use pycompat::text::{lower, strip};
use std::collections::{BTreeSet, HashMap, HashSet};

type Related = HashMap<String, Vec<String>>;

// --- single-position fallback --------------------------------------------

/// The CPIC warfarin guideline, which the single-position fallback cites.
pub const WARFARIN_GUIDELINE: &str =
    "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/";
const FALLBACK_LABEL: &str = "single-position fallback, not a PharmCAT call";

/// A gene whose CPIC warfarin guidance rests on one position. When PharmCAT
/// gives no single result for it (extra unphased heterozygous calls
/// elsewhere in CYP4F2 can leave no diplotype that matches), the genotype
/// PharmCAT itself read there, in the gene's `variants` list, still says
/// what the rule needs.
pub struct SinglePosition {
    pub gene: &'static str,
    pub rsid: &'static str,
    pub reference: &'static str,
    pub variant: &'static str,
    /// The genotype in words, from the number of variant alleles.
    pub spell: fn(usize) -> String,
    pub rule: &'static str,
}

/// CYP4F2 *3 and VKORC1 -1639. rs9923231 is VKORC1 -1639G>A; the gene is on
/// the minus strand, so the genomic C is -1639 G and T is -1639 A.
pub const SINGLE_POSITION: [SinglePosition; 2] = [
    SinglePosition {
        gene: "CYP4F2",
        rsid: "rs2108622",
        reference: "C",
        variant: "T",
        spell: |n| {
            let copies = ["no copy", "one copy", "two copies"][n];
            format!("{copies} of the *3 variant (T, p.V433M)")
        },
        rule: "the CPIC warfarin dosing algorithm for patients not of African ancestry adds 5-10% to the dose of a CYP4F2 *3 carrier.",
    },
    SinglePosition {
        gene: "VKORC1",
        rsid: "rs9923231",
        reference: "C",
        variant: "T",
        spell: |n| format!("-1639 {}{}", "G".repeat(2 - n), "A".repeat(n)),
        rule: "the CPIC warfarin guideline doses by a pharmacogenetic algorithm that takes VKORC1 -1639G>A with CYP2C9; each A allele lowers the expected dose.",
    },
];

/// The single-position fallback for each [`SINGLE_POSITION`] gene PharmCAT
/// left without a single result. A variant carrier is a WARNING line, so the
/// summary and the reports carry it; every line says it is not PharmCAT's
/// call. PharmCAT's own result for the gene is left as it is.
pub fn fallback_lines(calls: &[GeneCall], data: &Value) -> Vec<String> {
    let mut lines = Vec::new();
    for c in calls {
        let Some(sp) = SINGLE_POSITION.iter().find(|s| s.gene == c.gene) else {
            continue;
        };
        let status = c.status();
        if !matches!(status, Status::NotCalled | Status::Ambiguous) || c.outside {
            continue;
        }
        let what = if status == Status::NotCalled {
            "no result"
        } else {
            "no single result"
        };
        let Some(v) = position_read(data, &c.gene, sp.rsid) else {
            lines.push(format!(
                "  {} -- PharmCAT has {what}, and its report lists no genotype at {}: not callable.",
                c.gene, sp.rsid
            ));
            continue;
        };
        let call = truthy_field(v, "call").map(text).unwrap_or_default();
        let alleles: Vec<&str> = if call.is_empty() {
            Vec::new()
        } else {
            call.split(['/', '|']).collect()
        };
        if alleles.len() != 2
            || alleles
                .iter()
                .any(|a| *a != sp.reference && *a != sp.variant)
        {
            let got = if call.is_empty() {
                "has no genotype in the VCF PharmCAT read".to_owned()
            } else {
                format!(
                    "reads {call}, not two of {} and {}",
                    sp.reference, sp.variant
                )
            };
            lines.push(format!(
                "  {} -- PharmCAT has {what}; {} {got}: not callable.",
                c.gene, sp.rsid
            ));
            continue;
        }
        let n = alleles.iter().filter(|a| **a == sp.variant).count();
        let head = format!(
            "PharmCAT has {what} for {}; {FALLBACK_LABEL}: {} reads {call}, {}.",
            c.gene,
            sp.rsid,
            (sp.spell)(n)
        );
        if n > 0 {
            lines.push(format!("  WARNING: {head}"));
            lines.push(format!("    Warfarin (CPIC): {}", sp.rule));
            lines.push(format!("    {WARFARIN_GUIDELINE}"));
        } else {
            lines.push(format!(
                "  {head} The CPIC warfarin rule for this position does not apply."
            ));
        }
    }
    lines
}

// --- pypgx ----------------------------------------------------------------

/// One allele name in a form both tools share. PharmCAT names the alleles
/// of a single-position gene (VKORC1, IFNL3) `rs9923231 reference (C)` and
/// `rs9923231 variant (T)`, where pypgx writes `Reference` and `rs9923231`.
fn norm_allele(allele: &str) -> String {
    let a = strip(allele, None);
    let parts: Vec<&str> = a.split(' ').collect();
    if let [id, kind, bases] = parts[..]
        && !id.is_empty()
        && !id.chars().any(pycompat::text::is_space)
        && matches!(kind, "variant" | "reference")
        && bases
            .strip_prefix('(')
            .and_then(|b| b.strip_suffix(')'))
            .is_some_and(|b| !b.is_empty() && b.chars().all(|c| "ACGT".contains(c)))
    {
        return if kind == "reference" {
            "Reference".to_owned()
        } else {
            id.to_owned()
        };
    }
    if lower(a) == "reference" {
        "Reference".to_owned()
    } else {
        a.to_owned()
    }
}

/// A diplotype's alleles in a shared form and order: either tool may put
/// them in either order.
pub fn norm_diplotype(diplotype: &str) -> Vec<String> {
    let mut alleles: Vec<String> = diplotype.split('/').map(norm_allele).collect();
    alleles.sort();
    alleles
}

/// One row of the PharmCAT/pypgx comparison.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Comparison {
    pub gene: String,
    pub pharmcat: String,
    pub pypgx: String,
    pub matched: &'static str,
    pub called_by: &'static str,
    pub pharmcat_phenotype: String,
    pub pypgx_phenotype: String,
}

impl Comparison {
    pub fn cells(&self) -> [&str; 7] {
        [
            &self.gene,
            &self.pharmcat,
            &self.pypgx,
            self.matched,
            self.called_by,
            &self.pharmcat_phenotype,
            &self.pypgx_phenotype,
        ]
    }
}

/// Rows for the genes either tool called, by gene name. The phenotypes sit
/// beside the diplotypes so a gene pypgx has no phenotype rule for (step 32
/// writes `no pypgx phenotype rule`) still shows the phenotype PharmCAT gave.
pub fn compare(calls: &[GeneCall], pypgx: &HashMap<String, PypgxCall>) -> Vec<Comparison> {
    let genes: BTreeSet<&str> = calls
        .iter()
        .map(|c| c.gene.as_str())
        .chain(pypgx.keys().map(String::as_str))
        .collect();
    let mut rows = Vec::new();
    for gene in genes {
        let c = calls.iter().find(|c| c.gene == gene);
        let pc_ok = c.is_some_and(|c| c.called && !c.ambiguous());
        let pc_dip = match c {
            Some(c) if c.ambiguous() => {
                format!("ambiguous ({} possible diplotypes)", c.labels.len())
            }
            Some(c) => c.diplotype.clone(),
            None => "Not called".to_owned(),
        };
        let (pg_dip, pg_phen) = match pypgx.get(gene) {
            Some(p) if !p.diplotype.is_empty() => (p.diplotype.clone(), p.phenotype.clone()),
            Some(p) => ("Not called".to_owned(), p.phenotype.clone()),
            None => ("Not called".to_owned(), String::new()),
        };
        let pg_ok = pypgx_called(&pg_dip);
        let (matched, called_by) = match (pc_ok, pg_ok) {
            (false, false) => continue,
            (false, true) => ("pypgx only", "pypgx only"),
            (true, false) => ("PharmCAT only", "PharmCAT only"),
            (true, true) if norm_diplotype(&pc_dip) == norm_diplotype(&pg_dip) => ("Yes", "both"),
            (true, true) => ("No", "both"),
        };
        rows.push(Comparison {
            gene: gene.to_owned(),
            pharmcat: pc_dip,
            pypgx: pg_dip,
            matched,
            called_by,
            pharmcat_phenotype: c.map(GeneCall::shown).unwrap_or_default(),
            pypgx_phenotype: pg_phen,
        });
    }
    rows
}

/// Lines for genes PharmCAT reports on but could not call while pypgx did:
/// the medications section is silent on them, so say so. Genes step 36 held
/// back are left to its section.
pub fn pypgx_warnings(
    calls: &[GeneCall],
    pypgx: &HashMap<String, PypgxCall>,
    guidance: &GuidanceMap,
    related: &Related,
    held_back: &HashSet<String>,
) -> Vec<String> {
    let mut lines = Vec::new();
    for c in calls {
        if (c.called && !c.ambiguous()) || held_back.contains(&c.gene) {
            continue;
        }
        let Some(p) = pypgx.get(&c.gene) else {
            continue;
        };
        if !pypgx_called(&p.diplotype) {
            continue;
        }
        let what = if c.called {
            format!("no single result ({} possible diplotypes)", c.labels.len())
        } else {
            "no result".to_owned()
        };
        let (dip, phen) = (&p.diplotype, &p.phenotype);
        if is_normal(phen) {
            // pypgx's call is normal: worth knowing, nothing to review.
            lines.push(format!(
                "  NOTE: PharmCAT has {what} for {}; pypgx (step 32) called {dip} ({phen}).",
                c.gene
            ));
            continue;
        }
        let phen = if phen.is_empty() {
            "no phenotype"
        } else {
            phen
        };
        lines.push(format!(
            "  WARNING: PharmCAT has {what} for {gene}, but pypgx (step 32) called {dip} ({phen}). The medications above do not cover {gene}.",
            gene = c.gene
        ));
        // Every drug PharmCAT links to the gene, not only the ones whose
        // guidance happened to match another gene of this sample.
        let names = gene_drugs(&c.gene, guidance, related);
        if !names.is_empty() {
            lines.push(format!(
                "    Drugs affected by {}: {}",
                c.gene,
                names.join(", ")
            ));
        }
        lines.push("    Review them with the pypgx call and docs/32-pypgx.md.".to_owned());
    }
    lines
}

// --- outside calls (step 36) ---------------------------------------------

/// Marks every outside call of a step 36 gene that the consensus did not
/// pass on, with the same diplotype, as not called, and returns their
/// indexes. Such a call comes from an earlier step 36 run that step 07 has
/// not replaced, so it must give no drug guidance.
pub fn unconfirm_outside(calls: &mut [GeneCall], rows: &[ConsensusRow]) -> Vec<usize> {
    let mut passed: HashMap<&str, &str> = HashMap::new();
    for r in rows.iter().filter(|r| r.passed) {
        passed.insert(&r.gene, &r.result);
    }
    let mut stale = Vec::new();
    for (i, c) in calls.iter_mut().enumerate() {
        let confirmed = passed
            .get(c.gene.as_str())
            .is_some_and(|p| norm_diplotype(p) == norm_diplotype(&c.diplotype));
        if c.outside && c.called && STEP36_GENES.contains(&c.gene.as_str()) && !confirmed {
            c.called = false;
            c.phenotype = "outside call not confirmed by step 36".to_owned();
            stale.push(i);
        }
    }
    stale
}

/// `" Drugs affected by GENE: a, b."`, or nothing when there are none.
fn affected(gene: &str, names: &[String]) -> String {
    if names.is_empty() {
        String::new()
    } else {
        format!(" Drugs affected by {gene}: {}.", names.join(", "))
    }
}

/// The warning for each outside call the consensus does not pass on.
pub fn stale_lines(stale: &[&GeneCall], guidance: &GuidanceMap, related: &Related) -> Vec<String> {
    let mut lines = Vec::new();
    for c in stale {
        let names = gene_drugs(&c.gene, guidance, related);
        lines.push(format!(
            "  WARNING: PharmCAT's report has the outside call {} {}, which the consensus does not pass on: the report is older than step 36's last run.",
            c.gene, c.diplotype
        ));
        lines.push(format!(
            "           Rerun step 07, then this step. No drug guidance is given for {} here.{}",
            c.gene,
            affected(&c.gene, &names)
        ));
    }
    lines
}

/// What other tools gave PharmCAT and what was held back.
pub fn consensus_lines(
    rows: &[ConsensusRow],
    calls: &[GeneCall],
    guidance: &GuidanceMap,
    related: &Related,
) -> Vec<String> {
    let mut lines = Vec::new();
    for r in rows {
        let (gene, result) = (&r.gene, &r.result);
        if r.passed {
            let seen = if calls
                .iter()
                .any(|c| c.gene == *gene && c.outside && c.called)
            {
                "PharmCAT reports it as an outside call"
            } else {
                "PharmCAT's report does not list it as an outside call: rerun step 07"
            };
            lines.push(format!(
                "  {gene:<8} {result:<24} passed to PharmCAT from {}; {seen}.",
                r.reason
            ));
            lines.push(format!(
                "           {}",
                r.evidence.as_deref().unwrap_or("")
            ));
            continue;
        }
        lines.push(format!(
            "  {gene:<8} {result:<24} not passed to PharmCAT: {}.",
            r.reason
        ));
        if let Some(ev) = r.evidence.as_deref().filter(|e| *e != "-") {
            lines.push(format!("           {ev}"));
        }
        if result == "indeterminate" {
            let names = gene_drugs(gene, guidance, related);
            lines.push(format!(
                "           No drug guidance is given for {gene} here.{}",
                affected(gene, &names)
            ));
        }
    }
    lines
}

// --- phenoconversion -----------------------------------------------------

/// The genes a medication list is checked for, with the inhibitor table:
/// the FDA's inhibitors of each, and for CYP2D6 the factor CPIC's opioid
/// guideline multiplies the activity score by (0 strong, 0.5 moderate).
pub const PHENOCONVERSION_GENES: [&str; 5] = ["CYP2D6", "CYP2C19", "CYP2C9", "CYP3A4", "CYP3A5"];

/// What the phenotypes table writes for a gene the list leaves alone.
pub const NO_CONVERSION: &str = "-";

/// CPIC's wording of its CYP2D6 rule for a strength.
fn cpic_cyp2d6_rule(strength: &str) -> Option<&'static str> {
    match strength {
        "strong" => Some(
            "for patients on strong CYP2D6 inhibitors, the CYP2D6 activity score is adjusted to 0 and the predicted phenotype is a poor metabolizer.",
        ),
        "moderate" => Some(
            "for patients who are on moderate CYP2D6 inhibitors, the activity score is multiplied by 0.5 and then converted to the predicted phenotype.",
        ),
        _ => None,
    }
}

/// CPIC's CYP2D6 phenotype for an activity score.
pub fn cyp2d6_band(score: f64) -> &'static str {
    if score == 0.0 {
        "Poor Metabolizer"
    } else if score < 1.25 {
        "Intermediate Metabolizer"
    } else if score <= 2.25 {
        "Normal Metabolizer"
    } else {
        "Ultrarapid Metabolizer"
    }
}

/// The Phenoconverted_phenotype column of each gene the list touches, and
/// the section's lines.
pub struct Phenoconversion {
    pub column: HashMap<&'static str, String>,
    pub lines: Vec<String>,
}

/// Phenoconversion by an opt-in medication list. A line naming no drug of
/// the inhibitor table is listed as not checked. The genotype phenotype is
/// never changed: the converted one sits beside it. With several inhibitors
/// of one gene the strongest counts (CPIC gives no rule for combining them).
///
/// # Errors
///
/// A CYP2D6 row whose strength CPIC has no rule for: `KeyError`, as Python
/// raised, which the report quotes.
pub fn phenoconversion(
    calls: &[GeneCall],
    meds: &[Medication],
    inhibitors: &Inhibitors,
    guidance: &GuidanceMap,
    related: &Related,
) -> Result<Phenoconversion, PyErr> {
    let mut column = HashMap::new();
    let mut lines = Vec::new();
    let mut unknown: Vec<&str> = Vec::new();
    let mut names: Vec<&str> = Vec::new();
    for m in meds {
        if !inhibitors.keys().any(|(d, _)| *d == m.name) && !unknown.contains(&m.text.as_str()) {
            unknown.push(&m.text);
        }
        if !names.contains(&m.name.as_str()) {
            names.push(&m.name);
        }
    }
    for gene in PHENOCONVERSION_GENES {
        let hits: Vec<(&str, &str, Option<f64>)> = names
            .iter()
            .filter_map(|n| {
                inhibitors
                    .get(&((*n).to_owned(), gene.to_owned()))
                    .map(|i| (*n, i.strength.as_str(), i.factor))
            })
            .collect();
        if hits.is_empty() {
            continue;
        }
        let on = hits
            .iter()
            .map(|(n, st, _)| format!("{n} ({st} {gene} inhibitor)"))
            .collect::<Vec<_>>()
            .join(", ");
        let c = calls.iter().find(|c| c.gene == gene);
        let status = c.map(GeneCall::status);
        let Some(c) = c.filter(|_| matches!(status, Some(Status::Normal | Status::NonNormal)))
        else {
            let why = match status {
                None => "the report has no result for the gene".to_owned(),
                Some(s) => format!("the gene is {s} here"),
            };
            column.insert(gene, format!("not converted: {why}"));
            lines.push(format!(
                "  {gene} -- on {on}: not converted, {why}, so there is no genotype phenotype to convert."
            ));
            continue;
        };
        let drugs = gene_drugs(gene, guidance, related);
        if gene == "CYP2D6" {
            // CYP2D6 rows always carry a factor (read_inhibitors checks it).
            let factor = hits
                .iter()
                .filter_map(|h| h.2)
                .fold(f64::INFINITY, f64::min);
            let strongest: Vec<&str> = hits
                .iter()
                .filter(|h| h.2 == Some(factor))
                .map(|h| h.0)
                .collect();
            let strength = hits
                .iter()
                .find(|h| h.2 == Some(factor))
                .map_or("", |h| h.1);
            let scores: Vec<Option<f64>> = c.scores.clone();
            let all: Option<Vec<f64>> = scores.iter().copied().collect();
            let (converted, how) = if factor == 0.0 {
                (
                    "Poor Metabolizer",
                    "activity score adjusted to 0".to_owned(),
                )
            } else if let Some(all) = all.filter(|a| {
                !a.is_empty()
                    && a.iter()
                        .map(|s| cyp2d6_band(s * factor))
                        .collect::<HashSet<_>>()
                        .len()
                        == 1
            }) {
                let mut parts: Vec<String> = Vec::new();
                for s in &all {
                    let p = format!(
                        "{} x {} = {}",
                        general(*s, 6),
                        general(factor, 6),
                        general(s * factor, 6)
                    );
                    if !parts.contains(&p) {
                        parts.push(p);
                    }
                }
                (
                    cyp2d6_band(all[0] * factor),
                    format!("activity score {}", parts.join(", ")),
                )
            } else {
                let why = if scores.is_empty() || scores.contains(&None) {
                    "the report gives no activity score"
                } else {
                    "the possible diplotypes convert to different phenotypes"
                };
                column.insert(gene, format!("not converted: {why}"));
                lines.push(format!("  {gene} -- on {on}: not converted, {why}."));
                continue;
            };
            let label = format!("{converted} (phenoconverted by {})", strongest.join(", "));
            lines.push(format!(
                "  WARNING: CYP2D6 genotype phenotype {}; phenoconverted phenotype {label}, {how}.",
                c.shown()
            ));
            column.insert(gene, label);
            let rule = cpic_cyp2d6_rule(strength)
                .ok_or_else(|| PyErr::new("KeyError", crate::pyio::repr(strength)))?;
            lines.push(format!("    CPIC: {rule}"));
        } else {
            let why = if gene == "CYP2C19" {
                "CPIC has no CYP2C19 activity score".to_owned()
            } else {
                format!("CPIC gives no rule that converts the {gene} phenotype")
            };
            let inhibiting = hits.iter().map(|h| h.0).collect::<Vec<_>>().join(", ");
            column.insert(
                gene,
                format!("not computed: {inhibiting} inhibits {gene}; {why}"),
            );
            let mut first = why.chars();
            let why_cap: String = first
                .next()
                .map(|f| f.to_uppercase().chain(first).collect())
                .unwrap_or_default();
            lines.push(format!(
                "  WARNING: {gene} genotype phenotype {}; on {on}. {why_cap}, so no converted phenotype is computed.",
                c.shown()
            ));
            // A poor metabolizer (CYP3A5 *3/*3, the usual genotype) already
            // has no or little activity for an inhibitor to lower.
            if !lower(&c.phenotype).contains("poor metabolizer") {
                lines.push(format!(
                    "    The inhibitor can lower {gene} activity below what the genotype predicts."
                ));
            }
        }
        if !drugs.is_empty() {
            lines.push(format!(
                "    Drugs affected by {gene}: {}.",
                drugs.join(", ")
            ));
        }
    }
    if !unknown.is_empty() {
        lines.push(format!(
            "  Not in the inhibitor table: {}.",
            unknown.join(", ")
        ));
        lines.push("  These were not checked. A drug missing from the table (a brand name, a weak inhibitor)".to_owned());
        lines.push(format!(
            "  can still change the activity of {}.",
            PHENOCONVERSION_GENES.join(", ")
        ));
    }
    Ok(Phenoconversion { column, lines })
}

// --- disease-variant contraindications -----------------------------------

/// One part of CLNSIG (split at `|`, `/` and `,`) is Pathogenic or
/// Likely_pathogenic. `Conflicting_classifications_of_pathogenicity` and
/// `risk_factor` alone do not count.
pub fn is_plp(clnsig: Option<&str>) -> bool {
    clnsig.unwrap_or("").split(['|', '/', ',']).any(|p| {
        matches!(
            lower(strip(p, Some(" _"))).as_str(),
            "pathogenic" | "likely_pathogenic"
        )
    })
}

/// The gene symbols of ClinVar's GENEINFO (`POLG:5428|FOO:1`).
pub fn hit_genes(geneinfo: Option<&str>) -> BTreeSet<String> {
    geneinfo
        .unwrap_or("")
        .split('|')
        .map(|g| strip(g.split(':').next().unwrap_or(""), None))
        .filter(|g| !g.is_empty())
        .map(str::to_owned)
        .collect()
}

/// WARNING lines for each pathogenic or likely pathogenic hit in a gene of
/// the table, one per drug the table names for the gene. A drug whose names
/// hold a name of the medication list is marked as on it.
pub fn contraindication_lines(
    hits: &[ClinvarHit],
    table: &[Contraindication],
    med_names: &HashSet<String>,
) -> Vec<String> {
    let mut lines = Vec::new();
    for h in hits {
        if !is_plp(h.clnsig.as_deref()) {
            continue;
        }
        for gene in hit_genes(h.geneinfo.as_deref()) {
            let where_ = format!("{gene} {}:{} {}>{}", h.chrom, h.pos, h.ref_, h.alt);
            for r in table.iter().filter(|r| r.gene == gene) {
                let on_list = r
                    .names
                    .split(',')
                    .filter(|n| !strip(n, None).is_empty())
                    .any(|n| med_names.contains(&norm_drug(n)));
                let on = if on_list {
                    "; ON THE MEDICATION LIST"
                } else {
                    ""
                };
                lines.push(format!(
                    "  WARNING: {} -- {where_}, genotype {}, ClinVar {}{on}.",
                    r.drug,
                    h.genotype,
                    h.clnsig.as_deref().unwrap_or("None")
                ));
                lines.push(format!("    {} (CPIC level {})", r.note, r.cpic_level));
                let source = if r.label_source == "-" {
                    &r.pair_source
                } else {
                    &r.label_source
                };
                lines.push(format!("    {source}"));
            }
        }
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_position_alleles_match_pypgx() {
        assert_eq!(
            norm_diplotype("rs9923231 reference (C)/rs9923231 variant (T)"),
            norm_diplotype("Reference/rs9923231")
        );
        assert_ne!(
            norm_diplotype("rs12979860 variant (T)/rs12979860 variant (T)"),
            norm_diplotype("Reference/rs12979860")
        );
        assert_eq!(norm_diplotype("*2/*1"), norm_diplotype("*1/*2"));
        assert_eq!(norm_allele(" reference "), "Reference");
        assert_eq!(norm_allele("rs1 variant (TN)"), "rs1 variant (TN)");
    }

    #[test]
    fn cpic_bands() {
        for (s, want) in [
            (0.0, "Poor Metabolizer"),
            (0.25, "Intermediate Metabolizer"),
            (1.0, "Intermediate Metabolizer"),
            (1.25, "Normal Metabolizer"),
            (2.25, "Normal Metabolizer"),
            (2.5, "Ultrarapid Metabolizer"),
        ] {
            assert_eq!(cyp2d6_band(s), want, "{s}");
        }
    }

    #[test]
    fn clinvar_classes_and_genes() {
        assert!(is_plp(Some("Pathogenic/Likely_pathogenic")));
        assert!(is_plp(Some("Likely_pathogenic|risk_factor")));
        assert!(is_plp(Some("_Pathogenic_")));
        assert!(!is_plp(Some(
            "Conflicting_classifications_of_pathogenicity"
        )));
        assert!(!is_plp(Some("risk_factor")));
        assert!(!is_plp(None));
        let g: Vec<String> = hit_genes(Some("POLG:5428| FOO:1|:3")).into_iter().collect();
        assert_eq!(g, ["FOO", "POLG"]);
    }
}
