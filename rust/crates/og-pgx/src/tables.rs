//! The tables `cpic-report` reads beside the report: pypgx's summary, step
//! 36's consensus, the medication list, the inhibitor and contraindication
//! tables, and step 06's ClinVar hits. Each is read as Python's
//! `csv.DictReader(open(path), delimiter="\t")` reads it; a failure is the
//! exception Python would have raised, which the report quotes.

use crate::pyio::{PyErr, csv_error, lines, read_text, repr};
use pycompat::csv::{DictReader, DictRow, Quoting};
use pycompat::text::{lower, split, strip};
use std::collections::{HashMap, HashSet};

/// The header (`None` for an empty input) and rows of a TSV text.
fn dict_rows(text: &str) -> Result<(Option<Vec<String>>, Vec<DictRow>), PyErr> {
    let mut reader = DictReader::new(lines(text), '\t', Quoting::Minimal);
    let header = reader
        .header()
        .map_err(|e| csv_error(&e))?
        .map(|h| h.names().to_vec());
    let rows = reader
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| csv_error(&e))?;
    Ok((header, rows))
}

/// The text of `path` without the lines that start with `#`.
fn uncommented(path: &str) -> Result<String, PyErr> {
    let text = read_text(path)?;
    Ok(text
        .split_inclusive('\n')
        .filter(|l| !l.starts_with('#'))
        .collect())
}

/// `row.get(key)`: `None` for a missing column or a short row.
pub fn cell<'a>(row: &'a DictRow, key: &str) -> Option<&'a str> {
    row.get(key).flatten()
}

/// `str(row.get(key, default))`: `default` for a missing column, `None` for
/// a short row.
pub fn cell_or(row: &DictRow, key: &str, default: &str) -> String {
    match row.get(key) {
        None => default.to_owned(),
        Some(None) => "None".to_owned(),
        Some(Some(c)) => c.to_owned(),
    }
}

/// `repr()` of a list of strings.
fn repr_list<S: AsRef<str>>(items: &[S]) -> String {
    let parts: Vec<String> = items.iter().map(|s| repr(s.as_ref())).collect();
    format!("[{}]", parts.join(", "))
}

/// `repr()` of a row's dict, extra cells under the key `None`.
fn repr_row(row: &DictRow) -> String {
    let mut parts: Vec<String> = row
        .iter()
        .map(|(k, v)| format!("{}: {}", repr(k), v.map_or_else(|| "None".to_owned(), repr)))
        .collect();
    if let Some(extra) = row.extra() {
        parts.push(format!("None: {}", repr_list(extra)));
    }
    format!("{{{}}}", parts.join(", "))
}

// --- pypgx ---------------------------------------------------------------

/// pypgx's call of a gene: diplotype and phenotype.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PypgxCall {
    pub diplotype: String,
    pub phenotype: String,
}

/// `{gene: call}` from step 32's / PYPGX's summary TSV.
pub fn read_pypgx(path: &str) -> Result<HashMap<String, PypgxCall>, PyErr> {
    let (_, rows) = dict_rows(&read_text(path)?)?;
    let mut out = HashMap::new();
    for r in &rows {
        let gene = strip(cell(r, "Gene").unwrap_or(""), None);
        if !gene.is_empty() {
            let get = |k| strip(cell(r, k).unwrap_or(""), None).to_owned();
            out.insert(
                gene.to_owned(),
                PypgxCall {
                    diplotype: get("Diplotype"),
                    phenotype: get("Phenotype"),
                },
            );
        }
    }
    Ok(out)
}

/// A diplotype pypgx actually called.
pub fn pypgx_called(diplotype: &str) -> bool {
    !matches!(diplotype, "" | "FAILED" | "N/A" | "Not called")
}

// --- step 36 -------------------------------------------------------------

/// One row of the consensus table `bin/pgx_outside_calls.py` writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConsensusRow {
    pub gene: String,
    pub result: String,
    /// Passed to PharmCAT as an outside call.
    pub passed: bool,
    pub reason: String,
    /// `None` when the table has no Evidence column.
    pub evidence: Option<String>,
}

/// The rows of step 36's consensus table that name a gene.
pub fn read_consensus(path: &str) -> Result<Vec<ConsensusRow>, PyErr> {
    let (_, rows) = dict_rows(&read_text(path)?)?;
    Ok(rows
        .iter()
        .filter(|r| !strip(cell(r, "Gene").unwrap_or(""), None).is_empty())
        .map(|r| ConsensusRow {
            gene: cell(r, "Gene").unwrap_or_default().to_owned(),
            result: cell_or(r, "Result", ""),
            passed: cell(r, "Outside_call") == Some("yes"),
            reason: cell_or(r, "Reason", ""),
            evidence: r.get("Evidence").map(|_| cell_or(r, "Evidence", "")),
        })
        .collect())
}

// --- medications and inhibitors ------------------------------------------

/// A strength and CPIC's activity-score factor for it, when it has one.
#[derive(Debug, Clone, PartialEq)]
pub struct Inhibition {
    pub strength: String,
    pub factor: Option<f64>,
}

/// `{(drug, gene): inhibition}` of the inhibitor table.
pub type Inhibitors = HashMap<(String, String), Inhibition>;

/// A drug name as the tables compare it: lower case, whitespace collapsed.
pub fn norm_drug(name: &str) -> String {
    split(&lower(name), None, None).join(" ")
}

/// A number from a cell, when it is a finite one.
fn finite(s: &str) -> Option<f64> {
    pycompat::num::parse_float(s).ok().filter(|x| x.is_finite())
}

/// The inhibitor table (`assets/pgx_inhibitors.tsv`); lines starting with
/// `#` are comments. A row without a drug, a gene or a strength, a factor
/// that is not a number, and a CYP2D6 row without a factor from 0 to 1 are
/// errors: CYP2D6 rows feed CPIC's activity-score rule.
pub fn read_inhibitors(path: &str) -> Result<Inhibitors, PyErr> {
    let (_, rows) = dict_rows(&uncommented(path)?)?;
    let mut out = HashMap::new();
    for r in &rows {
        let text = |k| strip(cell(r, k).unwrap_or(""), None);
        if ["Drug", "Gene", "Strength"]
            .iter()
            .any(|k| text(k).is_empty())
        {
            return Err(PyErr::value(format!(
                "{path}: a row has an empty or missing Drug, Gene or Strength cell"
            )));
        }
        let (drug, gene) = (text("Drug"), text("Gene"));
        let fac = text("Activity_score_factor");
        let factor = if matches!(fac, "" | "-") {
            None
        } else {
            finite(fac)
        };
        if !matches!(fac, "" | "-") && factor.is_none() {
            return Err(PyErr::value(format!(
                "{path}: {} has a factor that is not a number: {}",
                repr(drug),
                repr(fac)
            )));
        }
        if gene == "CYP2D6" && !factor.is_some_and(|f| (0.0..=1.0).contains(&f)) {
            return Err(PyErr::value(format!(
                "{path}: CYP2D6 row {} needs a factor from 0 to 1",
                repr(drug)
            )));
        }
        out.insert(
            (norm_drug(cell(r, "Drug").unwrap_or("")), gene.to_owned()),
            Inhibition {
                strength: text("Strength").to_owned(),
                factor,
            },
        );
    }
    Ok(out)
}

/// One line of a medication list: as written, and the drug it names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Medication {
    pub text: String,
    pub name: String,
}

const DRUG_SEPARATORS: &str = " \t,;/+";

/// The lines of a medication list, one drug per line; `#` starts a comment
/// and blank lines are skipped. The name is the first word, lower case, so
/// a dose after it (`20 mg`) does not hide the drug. A slash or a plus also
/// ends a word, so `lopinavir/ritonavir` names lopinavir: the first drug of
/// a combination is the one checked.
pub fn read_medications(path: &str) -> Result<Vec<Medication>, PyErr> {
    let text = read_text(path)?;
    let mut out = Vec::new();
    for line in text.split_inclusive('\n') {
        let t = strip(line.split('#').next().unwrap_or(""), None);
        if t.is_empty() {
            continue;
        }
        let first = strip(t, Some(DRUG_SEPARATORS))
            .split(|c: char| pycompat::text::is_space(c) || ",;/+".contains(c))
            .next()
            .unwrap_or("");
        let name = norm_drug(first);
        if !name.is_empty() {
            out.push(Medication {
                text: t.to_owned(),
                name,
            });
        }
    }
    Ok(out)
}

// --- contraindications ---------------------------------------------------

/// The columns of `assets/pgx_contraindications.tsv`, in order.
pub const CONTRA_COLUMNS: [&str; 7] = [
    "Gene",
    "Drug",
    "Names",
    "CPIC_level",
    "Note",
    "Pair_source",
    "Label_source",
];

/// One rule of the contraindication table: a disease gene PharmCAT does not
/// call, and a drug a pathogenic variant in it makes risky.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contraindication {
    pub gene: String,
    pub drug: String,
    /// The drug's names, comma separated, lower case.
    pub names: String,
    pub cpic_level: String,
    pub note: String,
    pub pair_source: String,
    pub label_source: String,
}

/// The rows of the contraindication table; lines starting with `#` are
/// comments. A missing column or an empty cell is an error, never a row
/// read as half a rule.
pub fn read_contraindications(path: &str) -> Result<Vec<Contraindication>, PyErr> {
    let (header, rows) = dict_rows(&uncommented(path)?)?;
    if header.as_deref() != Some(&CONTRA_COLUMNS.map(str::to_owned)[..]) {
        let got = header
            .as_deref()
            .map_or_else(|| "None".to_owned(), repr_list);
        return Err(PyErr::value(format!(
            "{path}: the header is {got}, not {}",
            repr_list(&CONTRA_COLUMNS)
        )));
    }
    rows.iter()
        .map(|r| {
            let cells: Vec<&str> = CONTRA_COLUMNS
                .iter()
                .map(|k| strip(cell(r, k).unwrap_or(""), None))
                .collect();
            if cells.iter().any(|c| c.is_empty()) {
                return Err(PyErr::value(format!(
                    "{path}: a row has an empty or missing cell: {}",
                    repr_row(r)
                )));
            }
            let own = |i: usize| cells[i].to_owned();
            Ok(Contraindication {
                gene: own(0),
                drug: own(1),
                names: own(2),
                cpic_level: own(3),
                note: own(4),
                pair_source: own(5),
                label_source: own(6),
            })
        })
        .collect()
}

/// The columns of step 06's hits table this step reads.
const HIT_COLUMNS: [&str; 7] = [
    "chrom", "pos", "ref", "alt", "genotype", "geneinfo", "clnsig",
];

/// One row of step 06's / CLINVAR_SCREEN's hits table. A cell of a short
/// row reads as `None`, as Python formats it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClinvarHit {
    pub chrom: String,
    pub pos: String,
    pub ref_: String,
    pub alt: String,
    pub genotype: String,
    pub geneinfo: Option<String>,
    pub clnsig: Option<String>,
}

/// The rows of step 06's hits table.
pub fn read_clinvar_hits(path: &str) -> Result<Vec<ClinvarHit>, PyErr> {
    let (header, rows) = dict_rows(&read_text(path)?)?;
    let names: HashSet<&str> = header.iter().flatten().map(String::as_str).collect();
    let mut missing: Vec<&str> = HIT_COLUMNS
        .iter()
        .copied()
        .filter(|c| !names.contains(c))
        .collect();
    if !missing.is_empty() {
        missing.sort_unstable();
        return Err(PyErr::value(format!(
            "{path}: no column {}",
            missing.join(", ")
        )));
    }
    Ok(rows
        .iter()
        .map(|r| ClinvarHit {
            chrom: cell_or(r, "chrom", ""),
            pos: cell_or(r, "pos", ""),
            ref_: cell_or(r, "ref", ""),
            alt: cell_or(r, "alt", ""),
            genotype: cell_or(r, "genotype", ""),
            geneinfo: cell(r, "geneinfo").map(str::to_owned),
            clnsig: cell(r, "clnsig").map(str::to_owned),
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str, body: &str) -> String {
        let dir = std::env::temp_dir().join(format!("og-pgx-tables-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join(name);
        std::fs::write(&p, body).unwrap();
        p.to_str().unwrap().to_owned()
    }

    #[test]
    fn medication_lines_name_their_first_drug() {
        let p = tmp(
            "meds.txt",
            "Paroxetine 20 mg daily\n, paroxetine\n;\n# a comment\n\nlopinavir/ritonavir 400/100 mg\r\nValproate  # evening\n",
        );
        let got: Vec<(String, String)> = read_medications(&p)
            .unwrap()
            .into_iter()
            .map(|m| (m.text, m.name))
            .collect();
        let want = [
            ("Paroxetine 20 mg daily", "paroxetine"),
            (", paroxetine", "paroxetine"),
            ("lopinavir/ritonavir 400/100 mg", "lopinavir"),
            ("Valproate", "valproate"),
        ];
        assert_eq!(got, want.map(|(a, b)| (a.to_owned(), b.to_owned())));
    }

    #[test]
    fn the_shipped_tables_read() {
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets/");
        let inh = read_inhibitors(&format!("{root}pgx_inhibitors.tsv")).unwrap();
        assert_eq!(
            inh[&("paroxetine".to_owned(), "CYP2D6".to_owned())].factor,
            Some(0.0)
        );
        let contra = read_contraindications(&format!("{root}pgx_contraindications.tsv")).unwrap();
        assert!(
            contra
                .iter()
                .any(|c| c.gene == "POLG" && c.drug == "valproic acid")
        );
    }

    #[test]
    fn bad_tables_say_what_is_wrong() {
        let p = tmp(
            "short.tsv",
            "Drug\tGene\tStrength\tActivity_score_factor\nparoxetine\n",
        );
        assert_eq!(
            read_inhibitors(&p).unwrap_err().to_string(),
            format!("ValueError: {p}: a row has an empty or missing Drug, Gene or Strength cell")
        );
        let p = tmp(
            "nofactor.tsv",
            "Drug\tGene\tStrength\tActivity_score_factor\nparoxetine\tCYP2D6\tstrong\t\n",
        );
        assert_eq!(
            read_inhibitors(&p).unwrap_err().to_string(),
            format!("ValueError: {p}: CYP2D6 row 'paroxetine' needs a factor from 0 to 1")
        );
        let p = tmp(
            "nan.tsv",
            "Drug\tGene\tStrength\tActivity_score_factor\nx\tCYP2C9\tstrong\tlots\n",
        );
        assert_eq!(
            read_inhibitors(&p).unwrap_err().to_string(),
            format!("ValueError: {p}: 'x' has a factor that is not a number: 'lots'")
        );
        let p = tmp("contra.tsv", "Gene\tDrug\n");
        assert_eq!(
            read_contraindications(&p).unwrap_err().to_string(),
            format!(
                "ValueError: {p}: the header is ['Gene', 'Drug'], not ['Gene', 'Drug', 'Names', 'CPIC_level', 'Note', 'Pair_source', 'Label_source']"
            )
        );
        let p = tmp(
            "contra2.tsv",
            "Gene\tDrug\tNames\tCPIC_level\tNote\tPair_source\tLabel_source\nPOLG\tx\ty\n",
        );
        assert_eq!(
            read_contraindications(&p).unwrap_err().to_string(),
            format!(
                "ValueError: {p}: a row has an empty or missing cell: {{'Gene': 'POLG', 'Drug': 'x', 'Names': 'y', 'CPIC_level': None, 'Note': None, 'Pair_source': None, 'Label_source': None}}"
            )
        );
        let p = tmp("hits.tsv", "nothing\there\n");
        assert_eq!(
            read_clinvar_hits(&p).unwrap_err().to_string(),
            format!("ValueError: {p}: no column alt, chrom, clnsig, geneinfo, genotype, pos, ref")
        );
    }
}
