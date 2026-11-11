//! The drug guidance PharmCAT's report matched to the sample, and the
//! drug lists of a gene.

use crate::report::{GeneCall, field, list, text, truthy_field};
use pycompat::json::Value;
use pycompat::text::{lower, split, strip};
use std::collections::{BTreeSet, HashMap};

/// Fallback only: the drugs of a gene when the report itself names none
/// (no `drugs` section or `relatedDrugs`). Source:
/// <https://cpicpgx.org/guidelines/>. `docs/27-cpic-lookup.md` lists them.
pub const STATIC_DRUGS: &[(&str, &str)] = &[
    ("ABCG2", "rosuvastatin"),
    (
        "CACNA1S",
        "desflurane,enflurane,halothane,isoflurane,methoxyflurane,sevoflurane,succinylcholine",
    ),
    ("CYP2B6", "efavirenz,sertraline"),
    (
        "CYP2C19",
        "clopidogrel,voriconazole,escitalopram,citalopram,sertraline,amitriptyline,clomipramine,doxepin,imipramine,trimipramine,lansoprazole,omeprazole,pantoprazole,dexlansoprazole",
    ),
    (
        "CYP2C9",
        "warfarin,phenytoin,flurbiprofen,celecoxib,ibuprofen,lornoxicam,meloxicam,piroxicam,tenoxicam,siponimod",
    ),
    (
        "CYP2D6",
        "codeine,tramadol,hydrocodone,oxycodone,amitriptyline,clomipramine,desipramine,doxepin,imipramine,nortriptyline,trimipramine,fluvoxamine,paroxetine,atomoxetine,ondansetron,tropisetron,tamoxifen,eliglustat",
    ),
    ("CYP3A5", "tacrolimus"),
    ("CYP4F2", "warfarin"),
    ("DPYD", "fluorouracil,capecitabine,tegafur"),
    (
        "G6PD",
        "rasburicase,dapsone,chloroquine,primaquine,nitrofurantoin,methylene-blue",
    ),
    ("HLA-A", "carbamazepine,oxcarbazepine"),
    (
        "HLA-B",
        "abacavir,carbamazepine,oxcarbazepine,phenytoin,allopurinol",
    ),
    ("IFNL3", "peginterferon-alfa-2a,peginterferon-alfa-2b"),
    ("MT-RNR1", "aminoglycosides"),
    ("NUDT15", "azathioprine,mercaptopurine,thioguanine"),
    (
        "RYR1",
        "desflurane,enflurane,halothane,isoflurane,methoxyflurane,sevoflurane,succinylcholine",
    ),
    (
        "SLCO1B1",
        "simvastatin,atorvastatin,rosuvastatin,pravastatin,pitavastatin,fluvastatin,lovastatin",
    ),
    ("TPMT", "azathioprine,mercaptopurine,thioguanine"),
    ("UGT1A1", "atazanavir,belinostat,irinotecan"),
    ("VKORC1", "warfarin"),
];

/// The fallback drugs of `gene`, comma separated.
pub fn static_drugs(gene: &str) -> Option<&'static str> {
    STATIC_DRUGS
        .iter()
        .find(|(g, _)| *g == gene)
        .map(|(_, d)| *d)
}

/// A guidance source, in the order they are listed.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Source {
    Cpic,
    Dpwg,
    FdaLabel,
    FdaAssociation,
    /// Any other source, by its key in the report.
    Other(String),
}

impl Source {
    /// The source of a key of the report's `drugs` section.
    pub fn of(key: &str) -> Source {
        if key.starts_with("CPIC") {
            Source::Cpic
        } else if key.starts_with("DPWG") {
            Source::Dpwg
        } else if key.starts_with("FDA Label") {
            Source::FdaLabel
        } else if key.starts_with("FDA PGx") || key == "FDA association" {
            Source::FdaAssociation
        } else if key == "FDA label" {
            // A key spelled as a short name is that source.
            Source::FdaLabel
        } else {
            Source::Other(key.to_owned())
        }
    }

    /// The short name the text uses.
    pub fn name(&self) -> &str {
        match self {
            Source::Cpic => "CPIC",
            Source::Dpwg => "DPWG",
            Source::FdaLabel => "FDA label",
            Source::FdaAssociation => "FDA association",
            Source::Other(k) => k,
        }
    }

    /// The listing order: the four named sources, then every other one.
    fn rank(&self) -> u8 {
        match self {
            Source::Cpic => 0,
            Source::Dpwg => 1,
            Source::FdaLabel => 2,
            Source::FdaAssociation => 3,
            Source::Other(_) => 4,
        }
    }
}

/// One drug annotation of the report, for one gene. PharmCAT lists an
/// annotation for every diplotype the sample may have; `labels` and
/// `phenotype` say which diplotypes of the gene it is for.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Guidance {
    pub drug: String,
    pub source: Source,
    pub class: String,
    /// The recommendation, its whitespace collapsed.
    pub rec: String,
    /// The gene's diplotype labels the annotation names.
    pub labels: BTreeSet<String>,
    /// The phenotype the annotation gives the gene, or "".
    pub phenotype: String,
}

impl Guidance {
    /// What a drug list shows of it: drug, source, class and text.
    fn shown(&self) -> (&str, &Source, &str, &str) {
        (&self.drug, &self.source, &self.class, &self.rec)
    }
}

/// The guidance of each gene, best source first, then by drug name.
pub type GuidanceMap = HashMap<String, Vec<Guidance>>;

/// `" ".join(s.split())`.
pub fn collapse(s: &str) -> String {
    split(s, None, None).join(" ")
}

/// The report's own `drugs` section, per gene.
pub fn drug_guidance(data: &Value) -> GuidanceMap {
    let mut out: GuidanceMap = HashMap::new();
    let Some(Value::Object(drugs)) = field(data, "drugs") else {
        return out;
    };
    for (src_key, by_drug) in drugs {
        let Value::Object(by_drug) = by_drug else {
            continue;
        };
        let source = Source::of(src_key);
        for (drug_name, rep) in by_drug {
            if !matches!(rep, Value::Object(_)) {
                continue;
            }
            let drug = truthy_field(rep, "name").map_or_else(|| drug_name.clone(), text);
            for gl in list(rep, "guidelines") {
                for ann in list(gl, "annotations") {
                    if !matches!(ann, Value::Object(_)) {
                        continue;
                    }
                    add_annotation(&mut out, ann, &drug, &source);
                }
            }
        }
    }
    for entries in out.values_mut() {
        // A stable sort: equal keys keep the report's order.
        entries.sort_by_cached_key(|e| (e.source.rank(), lower(&e.drug)));
    }
    out
}

fn add_annotation(out: &mut GuidanceMap, ann: &Value, drug: &str, source: &Source) {
    let phen = match truthy_field(ann, "phenotypes") {
        Some(Value::Object(p)) => Some(p),
        _ => None,
    };
    let mut labels: Vec<(String, BTreeSet<String>)> = Vec::new();
    for gt in list(ann, "genotypes") {
        for dp in list(gt, "diplotypes") {
            let Some(gene) = truthy_field(dp, "gene") else {
                continue;
            };
            let gene = text(gene);
            let label = truthy_field(dp, "label").map(text).unwrap_or_default();
            match labels.iter_mut().find(|(g, _)| *g == gene) {
                Some((_, set)) => {
                    set.insert(label);
                }
                None => labels.push((gene, BTreeSet::from([label]))),
            }
        }
    }
    let rec = collapse(
        &truthy_field(ann, "drugRecommendation")
            .map(text)
            .unwrap_or_default(),
    );
    let class = truthy_field(ann, "classification")
        .map(text)
        .unwrap_or_default();
    let mut genes: Vec<&str> = phen
        .map(|p| p.keys().map(String::as_str).collect())
        .unwrap_or_default();
    for (g, _) in &labels {
        if !genes.contains(&g.as_str()) {
            genes.push(g);
        }
    }
    for gene in genes {
        let entry = Guidance {
            drug: drug.to_owned(),
            source: source.clone(),
            class: class.clone(),
            rec: rec.clone(),
            labels: labels
                .iter()
                .find(|(g, _)| g == gene)
                .map(|(_, s)| s.clone())
                .unwrap_or_default(),
            phenotype: phen
                .and_then(|p| p.get(gene))
                .filter(|v| pycompat::value::truthy(v))
                .map(text)
                .unwrap_or_default(),
        };
        let entries = out.entry(gene.to_owned()).or_default();
        if !entries.contains(&entry) {
            entries.push(entry);
        }
    }
}

/// The guidance entries for the diplotype PharmCAT called: those naming the
/// called diplotype, else those for the called phenotype exactly.
pub fn matching<'a>(call: &GeneCall, entries: &'a [Guidance]) -> Vec<&'a Guidance> {
    let by_label: Vec<&Guidance> = entries
        .iter()
        .filter(|e| e.labels.contains(&call.labels[0]))
        .collect();
    if !by_label.is_empty() {
        return by_label;
    }
    let called: Vec<&str> = call
        .phenotype
        .split(';')
        .map(|p| strip(p, None))
        .filter(|p| !p.is_empty())
        .collect();
    entries
        .iter()
        .filter(|e| {
            e.labels.is_empty()
                && !e.phenotype.is_empty()
                && called.contains(&strip(&e.phenotype, None))
        })
        .collect()
}

/// `re.sub(r"<[^>]+>", " ", s)`: each HTML tag becomes a space.
fn strip_tags(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut rest = s;
    while let Some(i) = rest.find('<') {
        let after = &rest[i + 1..];
        match after.find('>') {
            Some(j) if j > 0 => {
                out.push_str(&rest[..i]);
                out.push(' ');
                rest = &after[j + 1..];
            }
            _ => {
                out.push_str(&rest[..=i]);
                rest = after;
            }
        }
    }
    out.push_str(rest);
    out
}

/// The drugs of a gene: PharmCAT's links, else the fallback table, else
/// those its guidance names; each once.
pub fn gene_drugs(
    gene: &str,
    guidance: &GuidanceMap,
    related: &HashMap<String, Vec<String>>,
) -> Vec<String> {
    let names: Vec<String> = match related.get(gene) {
        Some(r) => r.clone(),
        None => match static_drugs(gene) {
            Some(s) => s.split(',').map(str::to_owned).collect(),
            None => guidance
                .get(gene)
                .map(|g| g.iter().map(|e| e.drug.clone()).collect())
                .unwrap_or_default(),
        },
    };
    let mut out: Vec<String> = Vec::with_capacity(names.len());
    for n in names {
        if !out.contains(&n) {
            out.push(n);
        }
    }
    out
}

/// The lines that list a gene's medications, best source first: CPIC's
/// matched recommendation per drug, then the drugs other sources name.
pub fn drug_lines(
    call: &GeneCall,
    guidance: &GuidanceMap,
    related: &HashMap<String, Vec<String>>,
) -> Vec<String> {
    let gene = &call.gene;
    let entries = matching(call, guidance.get(gene).map_or(&[], Vec::as_slice));
    if !entries.is_empty() {
        let mut lines = vec!["    Drugs with guidance in this PharmCAT report:".to_owned()];
        let mut shown: Vec<&str> = Vec::new();
        for e in &entries {
            if e.source != Source::Cpic || shown.contains(&e.drug.as_str()) {
                continue;
            }
            shown.push(&e.drug);
            let rec = collapse(&strip_tags(&e.rec));
            let class = if !e.class.is_empty() && e.class != "Unspecified" {
                format!(", {}", e.class)
            } else {
                String::new()
            };
            let mut head = format!("      - {} [CPIC{class}]", e.drug);
            if !rec.is_empty() {
                // whole: a cut loses dose qualifiers at the end
                head.push_str(": ");
                head.push_str(&rec);
            }
            lines.push(head);
        }
        let others: BTreeSet<&str> = entries
            .iter()
            .filter(|e| e.source != Source::Cpic && !shown.contains(&e.drug.as_str()))
            .map(|e| e.drug.as_str())
            .collect();
        if !others.is_empty() {
            let srcs: BTreeSet<(u8, &str)> = entries
                .iter()
                .filter(|e| e.source != Source::Cpic)
                .map(|e| (e.source.rank(), e.source.name()))
                .collect();
            let srcs: Vec<&str> = srcs.into_iter().map(|(_, n)| n).collect();
            lines.push(format!(
                "      also named by {}: {}",
                srcs.join(", "),
                others.into_iter().collect::<Vec<_>>().join(", ")
            ));
        }
        return lines;
    }
    if let Some(r) = related.get(gene) {
        return vec![format!(
            "    Drugs PharmCAT links to {gene} (no matched guidance text in the report): {}",
            r.join(", ")
        )];
    }
    if let Some(s) = static_drugs(gene) {
        return vec![format!(
            "    Drugs (pipeline fallback table; the report names none): {}",
            s.replace(',', ", ")
        )];
    }
    vec![format!(
        "    NOTE: {gene} has a non-normal phenotype but is not in the drug table; see the PharmCAT HTML report."
    )]
}

/// Where the phenotype of a gene came from when [`drug_level_phenotype`]
/// gave it.
pub const INDETERMINATE_NOTE: &str =
    "gene level Indeterminate; phenotype from PharmCAT's drug annotations";

/// The phenotype PharmCAT's drug annotations give a called gene whose only
/// gene-level phenotype is Indeterminate, when they all give the same one.
///
/// PharmCAT can leave a gene Indeterminate at gene level while the
/// annotations of its drugs assign it a phenotype; without this the two
/// levels of one report disagree. Annotations for any of the gene's
/// possible diplotypes are used when there are any, else every annotation
/// of the gene, so diplotypes PharmCAT cannot tell apart must agree. An
/// annotation that says Indeterminate, No Result or n/a disagrees; one with
/// no phenotype for the gene says nothing. Several possible diplotypes must
/// also get the same guidance, because [`drug_lines`] shows the first one's.
pub fn drug_level_phenotype(call: &GeneCall, guidance: &GuidanceMap) -> Option<String> {
    if !call.called || call.ambiguous() || lower(strip(&call.phenotype, None)) != "indeterminate" {
        return None;
    }
    let entries = guidance.get(&call.gene).map_or(&[][..], Vec::as_slice);
    let per_label: Vec<BTreeSet<(&str, &Source, &str, &str)>> = call
        .labels
        .iter()
        .map(|label| {
            entries
                .iter()
                .filter(|e| e.labels.contains(label))
                .map(Guidance::shown)
                .collect()
        })
        .collect();
    if per_label.iter().any(|g| *g != per_label[0]) {
        return None;
    }
    let for_labels: Vec<&Guidance> = entries
        .iter()
        .filter(|e| call.labels.iter().any(|l| e.labels.contains(l)))
        .collect();
    let used: Vec<&Guidance> = if for_labels.is_empty() {
        entries.iter().collect()
    } else {
        for_labels
    };
    let names: BTreeSet<&str> = used
        .iter()
        .map(|e| strip(&e.phenotype, None))
        .filter(|p| !p.is_empty())
        .collect();
    let [name] = names.into_iter().collect::<Vec<_>>()[..] else {
        return None;
    };
    (!matches!(lower(name).as_str(), "indeterminate" | "no result" | "n/a"))
        .then(|| name.to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tags_become_spaces() {
        assert_eq!(strip_tags("a<b>c</b>d"), "a c d");
        assert_eq!(strip_tags("a <> b"), "a <> b");
        assert_eq!(strip_tags("x<a<b>y"), "x y");
        assert_eq!(strip_tags("no close <here"), "no close <here");
    }

    #[test]
    fn sources_rank_cpic_first() {
        let mut s = [
            Source::of("FDA PGx Association"),
            Source::of("CPIC Guideline Annotation"),
            Source::of("DPWG Guideline Annotation"),
            Source::of("FDA Label Annotation"),
        ];
        s.sort_by_key(Source::rank);
        let names: Vec<&str> = s.iter().map(Source::name).collect();
        assert_eq!(names, ["CPIC", "DPWG", "FDA label", "FDA association"]);
        assert_eq!(Source::of("Other").name(), "Other");
    }

    #[test]
    fn the_static_table_has_every_gene_once() {
        let mut genes: Vec<&str> = STATIC_DRUGS.iter().map(|(g, _)| *g).collect();
        let n = genes.len();
        genes.dedup();
        assert_eq!(genes.len(), n);
        assert!(genes.is_sorted());
    }
}
