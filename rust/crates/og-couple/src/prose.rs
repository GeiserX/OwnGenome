//! The couple table and report: what two partners' findings in one gene
//! mean for a child, gene by gene, then in words.

use crate::findings::{DEPTH_MIN, Finding, SV_MIN_SIZE};
use crate::list::{CarrierList, Dominant, Gene, Inheritance};
use std::collections::{BTreeSet, HashMap};
use std::fmt::Write as _;

/// The columns of the couple table, in order.
pub const TABLE_COLUMNS: [&str; 11] = [
    "gene",
    "inheritance",
    "lists",
    "condition",
    "result",
    "offspring_risk",
    "partner_a",
    "findings_a",
    "partner_b",
    "findings_b",
    "explanation",
];

/// What a gene's findings mean for the couple, in table order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Outcome {
    /// Both partners carry a change in one recessive gene.
    Shared,
    /// A change in an X-linked gene.
    XLinked,
    /// Only one partner carries a change.
    OnePartner,
}

impl Outcome {
    pub fn code(self) -> &'static str {
        match self {
            Outcome::Shared => "shared",
            Outcome::XLinked => "x_linked",
            Outcome::OnePartner => "one_partner",
        }
    }
}

/// One row of the couple table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub gene: String,
    pub inheritance: &'static str,
    pub lists: String,
    pub condition: String,
    pub outcome: Outcome,
    pub risk: String,
    pub partner_a: String,
    pub findings_a: String,
    pub partner_b: String,
    pub findings_b: String,
    pub explanation: String,
}

impl Row {
    fn cells(&self) -> [&str; 11] {
        [
            &self.gene,
            self.inheritance,
            &self.lists,
            &self.condition,
            self.outcome.code(),
            &self.risk,
            &self.partner_a,
            &self.findings_a,
            &self.partner_b,
            &self.findings_b,
            &self.explanation,
        ]
    }
}

/// `1 in 4`, `1 in 2`, `every`.
fn odds(quarters: u8) -> &'static str {
    match quarters {
        1 => "1 in 4",
        2 => "1 in 2",
        _ => "every",
    }
}

fn two_copies(f: &[&Finding]) -> bool {
    f.iter().any(|x| x.copies == 2)
}

/// Both partners carry a change in one autosomal recessive gene.
fn ar_row(
    g: &Gene,
    a: &str,
    b: &str,
    fa: &[&Finding],
    fb: &[&Finding],
) -> (Outcome, String, String) {
    // The chance a child inherits each parent's change, in halves.
    let half = |f: &[&Finding]| if two_copies(f) { 2u8 } else { 1 };
    let quarters = half(fa) * half(fb);
    let risk = if quarters == 4 {
        format!("{} pregnancy", odds(quarters))
    } else {
        format!("{} per pregnancy", odds(quarters))
    };
    let mut why = format!(
        "Both partners carry a change in {}, linked to {} (autosomal recessive). A child who inherits the change \
         of each parent can have the condition: {risk}.",
        g.name, g.condition
    );
    let low: Vec<&str> = [(a, fa), (b, fb)]
        .iter()
        .filter(|(_, f)| !f.is_empty() && f.iter().all(|x| x.hypomorph.is_some()))
        .map(|(n, _)| *n)
        .collect();
    if !low.is_empty() {
        let _ = write!(
            why,
            " The change of {} is a low-function (hypomorph) allele, not a full loss of function: with the other \
             parent's loss-of-function change a child can have the condition, for some genes in a milder form; two \
             such alleles alone are set aside.",
            low.join(" and ")
        );
    }
    let more: Vec<&str> = [(a, fa), (b, fb)]
        .iter()
        .filter(|(_, f)| f.len() > 1 && f.iter().all(|x| x.copies == 1))
        .map(|(n, _)| *n)
        .collect();
    if !more.is_empty() {
        let _ = write!(
            why,
            " {} {} two changes in the gene; when they sit on different copies, every child gets one, and the \
             chance is higher.",
            more.join(" and "),
            if more.len() == 1 { "has" } else { "have" }
        );
    }
    (Outcome::Shared, risk, why)
}

/// One partner's lines for an X-linked dominant gene: the risk and the
/// explanation. A child who inherits the change can be affected, a
/// daughter as well as a son.
fn xld_lines(g: &Gene, who: &str, sex: Option<&str>, f: &[&Finding]) -> (String, String) {
    match sex {
        Some("female") if two_copies(f) => (
            format!("{who}: every child, can be affected"),
            format!(
                "{who} carries the change in {} on both X chromosomes (X-linked dominant: {}): every son and \
                 every daughter inherits it and can be affected.",
                g.name, g.condition
            ),
        ),
        Some("female") => (
            format!("{who}: each child 1 in 2, can be affected"),
            format!(
                "{who} carries a change in {} on one X chromosome (X-linked dominant: {}). Each child has a 1 in \
                 2 chance of inheriting it, and a daughter or a son who inherits it can be affected; {who} can \
                 have signs herself.",
                g.name, g.condition
            ),
        ),
        Some("male") => (
            format!("{who}: every daughter, can be affected"),
            format!(
                "{who} has a change in {} on his X chromosome (X-linked dominant: {}). Every daughter inherits it \
                 and can be affected; a son gets his Y chromosome, not this X.",
                g.name, g.condition
            ),
        ),
        _ => (
            format!("{who}: sex unknown"),
            format!(
                "{who} has a change in {} (X-linked dominant: {}). The pedigree gives no sex for {who}, so the \
                 chance for sons and daughters cannot be said.",
                g.name, g.condition
            ),
        ),
    }
}

/// One partner's lines for an X-linked gene: the risk and the explanation.
fn xl_lines(g: &Gene, who: &str, sex: Option<&str>, f: &[&Finding]) -> (String, String) {
    if g.inheritance == Inheritance::Xld {
        return xld_lines(g, who, sex, f);
    }
    let (risk, mut why) = xlr_lines(g, who, sex, f);
    if g.dominant == Dominant::XLinked {
        let _ = write!(
            why,
            " {} also has an X-linked dominant disorder: females can show signs.",
            g.name
        );
    }
    (risk, why)
}

/// One partner's lines for an X-linked recessive gene.
fn xlr_lines(g: &Gene, who: &str, sex: Option<&str>, f: &[&Finding]) -> (String, String) {
    let fmr1 = g.name == "FMR1" && f.iter().any(|x| x.text.contains("premutation"));
    match sex {
        Some("female") => {
            let (risk, mut why) = if two_copies(f) {
                (
                    format!("{who}: every son, every daughter a carrier"),
                    format!(
                        "{who} carries the change on both X chromosomes: every son inherits it and can have {}; \
                         every daughter is a carrier.",
                        g.condition
                    ),
                )
            } else {
                (
                    format!("{who}: sons 1 in 2"),
                    format!(
                        "{who} carries a change in {} on one X chromosome (X-linked: {}). Each son has a 1 in 2 \
                         chance of inheriting it and can have the condition; each daughter has a 1 in 2 chance of \
                         being a carrier.",
                        g.name, g.condition
                    ),
                )
            };
            if fmr1 {
                why.push_str(
                    " A premutation can grow into a full mutation when a mother passes it on.",
                );
            }
            (risk, why)
        }
        Some("male") => (
            format!("{who}: daughters carriers"),
            format!(
                "{who} has a change in {} on his X chromosome (X-linked: {}). Every daughter inherits it and is a \
                 carrier; a son gets his Y chromosome, not this X.",
                g.name, g.condition
            ),
        ),
        _ => (
            format!("{who}: sex unknown"),
            format!(
                "{who} has a change in {} (X-linked: {}). The pedigree gives no sex for {who}, so the chance for \
                 sons and daughters cannot be said.",
                g.name, g.condition
            ),
        ),
    }
}

fn joined(f: &[&Finding]) -> String {
    if f.is_empty() {
        ".".to_owned()
    } else {
        f.iter()
            .map(|x| x.text.as_str())
            .collect::<Vec<_>>()
            .join("; ")
    }
}

/// The rows of the couple table: shared genes first, then X-linked, then
/// the rest, each by gene. A hypomorph counts only when the other partner
/// has a full loss of function in the same gene.
pub fn screen_couple(
    a: &str,
    b: &str,
    sexes: &HashMap<String, Option<String>>,
    found: &HashMap<String, Vec<Finding>>,
    list: &CarrierList,
) -> Vec<Row> {
    let none = Vec::new();
    let of = |s: &str| found.get(s).unwrap_or(&none);
    let full = |s: &str| -> BTreeSet<&str> {
        of(s)
            .iter()
            .filter(|x| x.hypomorph.is_none())
            .map(|x| x.gene.as_str())
            .collect()
    };
    let (full_a, full_b) = (full(a), full(b));
    let sex = |s: &str| sexes.get(s).and_then(|x| x.as_deref());
    let mut rows = Vec::new();
    for name in full_a.union(&full_b) {
        let Some(g) = list.get(name) else { continue };
        let pick = |s: &str, other: &BTreeSet<&str>| -> Vec<&Finding> {
            of(s)
                .iter()
                .filter(|x| x.gene == *name && (x.hypomorph.is_none() || other.contains(name)))
                .collect()
        };
        let (fa, fb) = (pick(a, &full_b), pick(b, &full_a));
        let (outcome, risk, mut why) = if g.inheritance.x_linked() {
            let parts: Vec<(String, String)> = [(a, &fa), (b, &fb)]
                .iter()
                .filter(|(_, f)| !f.is_empty())
                .map(|(w, f)| xl_lines(g, w, sex(w), f))
                .collect();
            (
                Outcome::XLinked,
                parts
                    .iter()
                    .map(|p| p.0.as_str())
                    .collect::<Vec<_>>()
                    .join("; "),
                parts
                    .iter()
                    .map(|p| p.1.as_str())
                    .collect::<Vec<_>>()
                    .join(" "),
            )
        } else if !fa.is_empty() && !fb.is_empty() {
            let (o, r, mut w) = ar_row(g, a, b, &fa, &fb);
            if g.dominant == Dominant::Autosomal {
                let _ = write!(
                    w,
                    " {name} also causes a dominant condition (it is on the ACMG secondary findings list): a child \
                     who inherits either parent's change can be affected by it."
                );
            }
            (o, r, w)
        } else if g.dominant == Dominant::Autosomal {
            let (who, other, f) = if fa.is_empty() {
                (b, a, &fb)
            } else {
                (a, b, &fa)
            };
            if two_copies(f) {
                // Both copies changed: every child inherits one of them.
                (
                    Outcome::OnePartner,
                    "every child, dominant".to_owned(),
                    format!(
                        "Only {who} carries a change in {name}, on both copies. {name} also causes a dominant \
                         condition (it is on the ACMG secondary findings list): every child inherits it and can be \
                         affected, and {who} can be too. For {} (autosomal recessive), a child can also have the \
                         condition if {other} carries a change this screen does not see.",
                        g.condition
                    ),
                )
            } else {
                (
                    Outcome::OnePartner,
                    "each child 1 in 2, dominant".to_owned(),
                    format!(
                        "Only {who} carries a change in {name}. {name} also causes a dominant condition (it is on \
                         the ACMG secondary findings list): a child who inherits it can be affected, 1 in 2 per \
                         child, and {who} can be too. For {} (autosomal recessive), a child can be a carrier like \
                         {who}, unless {other} carries a change this screen does not see.",
                        g.condition
                    ),
                )
            }
        } else {
            let (who, other) = if fa.is_empty() { (b, a) } else { (a, b) };
            (
                Outcome::OnePartner,
                "none shared".to_owned(),
                format!(
                    "Only {who} carries a change in {name} ({}, autosomal recessive). A child can be a carrier \
                     like {who}, not affected, unless {other} carries a change this screen does not see.",
                    g.condition
                ),
            )
        };
        let weak: Vec<&str> = [(a, &fa), (b, &fb)]
            .iter()
            .filter(|(_, f)| f.iter().any(|x| x.low))
            .map(|(w, _)| *w)
            .collect();
        if !weak.is_empty() {
            let _ = write!(
                why,
                " The finding of {} rests on a low-confidence call (GQ or DP below step 06's floor): confirm it in \
                 the reads, or with a clinical test, before reading the chance.",
                weak.join(" and ")
            );
        }
        rows.push(Row {
            gene: (*name).to_owned(),
            inheritance: g.inheritance.code(),
            lists: g.lists.clone(),
            condition: g.condition.clone(),
            outcome,
            risk,
            partner_a: a.to_owned(),
            findings_a: joined(&fa),
            partner_b: b.to_owned(),
            findings_b: joined(&fb),
            explanation: why,
        });
    }
    rows.sort_by(|x, y| (x.outcome, &x.gene).cmp(&(y.outcome, &y.gene)));
    rows
}

/// The couple table as text: a header, then a line per row. A tab or a
/// newline inside a cell reads as a space.
pub fn table(rows: &[Row]) -> String {
    let mut out = TABLE_COLUMNS.join("\t");
    out.push('\n');
    for r in rows {
        let cells: Vec<String> = r
            .cells()
            .iter()
            .map(|c| c.replace(['\t', '\n'], " "))
            .collect();
        out.push_str(&cells.join("\t"));
        out.push('\n');
    }
    out
}

/// What the screen does not see, whatever the inputs.
fn residual() -> [String; 9] {
    [
        "Changes ClinVar does not list as pathogenic or likely pathogenic: a new loss-of-function change, or one \
         of uncertain significance, is not counted."
            .to_owned(),
        "Genes outside the carrier list (assets/carrier_genes.tsv).".to_owned(),
        "SMN1 silent carriers: two copies on one chromosome and none on the other read as two copies.".to_owned(),
        "Changes in genes with a near-identical copy (HBA1 and HBA2, CYP21A2, GBA1, SMN1): short reads place them \
         poorly, and a gene conversion is not seen. A deletion there counts only when Delly calls it; for CYP21A2, \
         step 35's table (--cyp21a2) adds the 30 kb chimeric deletion and conversions at six of its pseudogene \
         variant sites, nothing else."
            .to_owned(),
        "Repeat expansions other than FMR1 and FXN, and the exact size of a long FMR1 or FXN allele: short reads \
         do not span it."
            .to_owned(),
        format!(
            "Deletions Delly does not call, deletions longer than {DEPTH_MIN} bp that read depth (Delly's RDCN) does \
             not confirm, duplications, inversions and insertions."
        ),
        format!(
            "Delly deletions shorter than {SV_MIN_SIZE} bp: they count only when the small-variant caller also calls \
             them and ClinVar lists them, so one it misses is missed here."
        ),
        format!(
            "With AnnotSV's benign set, a deletion that matches a common one there (1% or more, {}% reciprocal \
             overlap) is not counted: a rare pathogenic deletion at the same place as a common one is missed.",
            crate::findings::COMMON_OVERLAP_PERCENT
        ),
        "Low-level mosaicism, and whether two changes of one person in one gene sit on the same copy.".to_owned(),
    ]
}

const COUNSEL: &str = "This is a screen, not a diagnosis, and not medical advice. Before acting on a result, or on \
its absence, see a genetic counsellor: a clinical carrier test confirms a finding, and the counsellor can say what a \
negative screen leaves out for your ancestry.";

/// What was read for one partner, for the report's table.
#[derive(Debug, Clone, Default)]
pub struct Assessed {
    pub clinvar: String,
    pub delly: String,
    pub smn: String,
    pub eh: String,
}

/// Everything the report of one couple says.
pub struct Report<'a> {
    pub a: &'a str,
    pub b: &'a str,
    pub sexes: &'a HashMap<String, Option<String>>,
    pub rows: &'a [Row],
    pub assessed: &'a HashMap<String, Assessed>,
    pub gaps: &'a [String],
    pub genes: usize,
    pub common: &'a str,
    pub cn: &'a [String],
    pub shared_cn: &'a [String],
}

impl Report<'_> {
    fn section(&self, out: &mut Vec<String>, title: &str, keep: Outcome, none: &str) {
        out.push(format!("## {title}"));
        out.push(String::new());
        let sel: Vec<&Row> = self.rows.iter().filter(|r| r.outcome == keep).collect();
        if sel.is_empty() {
            out.push(none.to_owned());
            out.push(String::new());
            return;
        }
        for r in sel {
            out.push(format!(
                "- **{}** ({}; {}): {}",
                r.gene, r.condition, r.lists, r.risk
            ));
            for (p, f) in [(&r.partner_a, &r.findings_a), (&r.partner_b, &r.findings_b)] {
                if f != "." {
                    out.push(format!("  - {p}: {f}"));
                }
            }
            out.push(format!("  - {}", r.explanation));
        }
        out.push(String::new());
    }

    /// The report as Markdown.
    pub fn markdown(&self) -> String {
        let sx = |s: &str| {
            self.sexes
                .get(s)
                .and_then(|x| x.as_deref())
                .unwrap_or("sex unknown")
                .to_owned()
        };
        let (a, b) = (self.a, self.b);
        let mut out = vec![
            format!("# Couple carrier screen: {a} and {b}"),
            String::new(),
            format!(
                "{a} ({}) and {b} ({}), partners in the pedigree. The screen covers {} genes of autosomal \
                 recessive and X-linked conditions (assets/carrier_genes.tsv: the ACMG and ACOG carrier screening \
                 genes and the recessive and X-linked genes of Orphanet).",
                sx(a),
                sx(b),
                self.genes
            ),
            String::new(),
            COUNSEL.to_owned(),
            String::new(),
        ];
        let shared_none = if self.shared_cn.is_empty() {
            "No gene of the list where both partners carry a change."
        } else {
            "No gene of the list where both partners carry a change in the gene table; read depth (step 35) adds a \
             possible shared risk, below."
        };
        self.section(
            &mut out,
            "Shared risk: both partners carry a change in the same gene",
            Outcome::Shared,
            shared_none,
        );
        if !self.shared_cn.is_empty() {
            out.push(
                "Possible shared risk from read depth (step 35's carrier deletion alleles, not in the gene table):"
                    .to_owned(),
            );
            out.push(String::new());
            out.extend(self.shared_cn.iter().map(|x| format!("- {x}")));
            out.push(String::new());
        }
        self.section(
            &mut out,
            "X-linked findings",
            Outcome::XLinked,
            "No change in an X-linked gene of the list.",
        );
        self.section(
            &mut out,
            "Genes where only one partner carries a change",
            Outcome::OnePartner,
            "None.",
        );
        if !self.cn.is_empty() {
            out.push(
                "## Carrier deletion alleles and AZF regions (read depth, step 35)".to_owned(),
            );
            out.push(String::new());
            out.push(
                "Each allele is present, absent or not assessable from the read depth over it. A present deletion \
                 needs a clinical test. The AZF regions are read on a male partner only. The three beta-globin rows \
                 overlap: the Sicilian deletion covers both HBB rows, and the 619 bp and exon 3 rows share 487 bp, \
                 so one deletion can read as present on two or three of them; only a clinical test tells which \
                 allele it is."
                    .to_owned(),
            );
            out.push(String::new());
            out.extend(self.cn.iter().cloned());
            out.push(String::new());
        }
        out.push("## What was assessed".to_owned());
        out.push(String::new());
        out.push(
            "| Partner | ClinVar hits (in listed genes / all) | Deletions | SMN1 | FMR1 and FXN repeats |".to_owned(),
        );
        out.push("|---|---|---|---|---|".to_owned());
        for s in [a, b] {
            let x = &self.assessed[s];
            out.push(format!(
                "| {s} | {} | {} | {} | {} |",
                x.clinvar, x.delly, x.smn, x.eh
            ));
        }
        out.push(String::new());
        out.push(format!(
            "Deletions shorter than {SV_MIN_SIZE} bp are left to the small-variant caller. Common deletions: {}.",
            self.common
        ));
        out.push(String::new());
        out.push("## What this screen does not see".to_owned());
        out.push(String::new());
        out.extend(
            self.gaps
                .iter()
                .cloned()
                .chain(residual())
                .map(|g| format!("- {g}")),
        );
        out.push(String::new());
        out.push(
            "A screen that finds no shared gene lowers the chance of a child with one of these conditions; it does \
             not remove it."
                .to_owned(),
        );
        out.push(String::new());
        out.join("\n")
    }
}
