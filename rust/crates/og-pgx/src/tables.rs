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

/// A diplotype pypgx actually called. `Indeterminate` is no call, in any
/// case and in any form: pypgx writes it for a gene it could not call,
/// step 32 writes it for CYP2D6 when the depth check fails, and a
/// diplotype with an Indeterminate allele (`*1/Indeterminate`) names no
/// diplotype either. The rule is `og_steps::pgx_rules::is_call`, the one
/// summary.json counts with (og-pgx cannot depend on og-steps): two
/// alleles, neither a no-call text, and one diplotype, not pypgx's
/// `;`-joined list of possible ones.
pub fn pypgx_called(diplotype: &str) -> bool {
    const NO_CALL: [&str; 9] = [
        "",
        "none",
        "n/a",
        "na",
        "failed",
        "indeterminate",
        "not called",
        "unknown",
        ".",
    ];
    let no_call = |a: &str| NO_CALL.contains(&lower(strip(a, None)).as_str());
    let d = strip(diplotype, None);
    !no_call(d)
        && d.contains('/')
        && !d.contains(';')
        && !d.split('/').any(no_call)
        && !lower(d).contains("indeterminate")
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

/// One line of a medication list: as written, and every drug name in it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Medication {
    pub text: String,
    /// The words of the line that can name a drug, lower case, in order,
    /// each once. Never empty.
    pub names: Vec<String>,
}

/// Words of a medication line that name no drug: units, dose forms and how
/// often it is taken.
const DOSE_WORDS: &[&str] = &[
    "mg", "mcg", "µg", "ug", "g", "kg", "ml", "l", "iu", "unit", "units", "%", "tablet", "tablets",
    "tab", "tabs", "capsule", "capsules", "cap", "caps", "er", "xr", "sr", "xl", "daily", "once",
    "twice", "a", "day", "per", "weekly", "morning", "evening", "night",
];

/// A word that is a dose, a unit or a frequency, not a drug: it starts with
/// a digit or a point (`20`, `400mg`, `.5`), or it is one of the dose words.
fn dose_word(w: &str) -> bool {
    w.starts_with(|c: char| c.is_ascii_digit() || c == '.') || DOSE_WORDS.contains(&w)
}

/// The parts of a hyphenated word that can name a drug, in order: each part
/// that starts with a letter, has three or more characters and is no dose
/// word. `olanzapine-fluoxetine` gives olanzapine and fluoxetine; the `co`
/// of `co-trimoxazole`, the `l` of `l-dopa` and the `5` of `5-htp` are
/// never drugs, so they give nothing. A word without a hyphen has no parts.
fn hyphen_parts(w: &str) -> impl Iterator<Item = &str> {
    w.split('-').filter(move |p| {
        *p != w && p.chars().count() >= 3 && p.starts_with(char::is_alphabetic) && !dose_word(p)
    })
}

/// The lines of a medication list, one product per line; `#` starts a
/// comment and blank lines are skipped. Every word of a line can name a
/// drug, so each drug of a combination product is checked: whitespace,
/// commas, semicolons, `/`, `+` and brackets split the line, punctuation at
/// a word's ends is dropped (`paroxetine.` is paroxetine), and the word
/// `and` is dropped (`olanzapine/fluoxetine`, `naltrexone/bupropion 8/90 mg`
/// and `lopinavir and ritonavir` each name two drugs; `Contrave
/// (naltrexone/bupropion)` names three words). Doses, units and frequencies
/// (`20`, `mg`, `daily`) are dropped too. A hyphenated word names itself
/// first, so a hyphenated name of a table still matches whole, then its
/// parts (`hyphen_parts`), so `olanzapine-fluoxetine` is checked for
/// fluoxetine. A part is only a name: it hits only when it is itself a name
/// of a table, and no table name is a part of another one (a test checks
/// the shipped tables). A line of nothing but those keeps
/// its first word, so it is still listed as not checked; a line of only
/// separators and punctuation is skipped.
pub fn read_medications(path: &str) -> Result<Vec<Medication>, PyErr> {
    let text = read_text(path)?;
    let mut out = Vec::new();
    for line in text.split_inclusive('\n') {
        let t = strip(line.split('#').next().unwrap_or(""), None);
        if t.is_empty() {
            continue;
        }
        let words: Vec<String> = t
            .split(|c: char| pycompat::text::is_space(c) || ",;/+()[]".contains(c))
            .map(|w| w.trim_matches(|c: char| !c.is_alphanumeric()))
            .filter(|w| !w.is_empty())
            .map(lower)
            .collect();
        let mut names: Vec<String> = Vec::new();
        for w in &words {
            if w == "and" || dose_word(w) {
                continue;
            }
            for n in std::iter::once(w.as_str()).chain(hyphen_parts(w)) {
                if !names.iter().any(|m| m == n) {
                    names.push(n.to_owned());
                }
            }
        }
        if names.is_empty()
            && let Some(first) = words.first()
        {
            names.push(first.clone());
        }
        if !names.is_empty() {
            out.push(Medication {
                text: t.to_owned(),
                names,
            });
        }
    }
    Ok(out)
}

// --- contraindications ---------------------------------------------------

/// The columns of `assets/pgx_contraindications.tsv`, in order.
pub const CONTRA_COLUMNS: [&str; 8] = [
    "Gene",
    "Drug",
    "Names",
    "CPIC_level",
    "Note",
    "Pair_source",
    "Label_source",
    "Report_if",
];

/// When a pathogenic hit in a row's gene gives a WARNING (the table's
/// `Report_if` column, in the words of `assets/acmg_sf.tsv`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportIf {
    /// One variant is enough (`any`).
    Any,
    /// The disorder is recessive: both copies changed, or two variants in
    /// the gene (`biallelic`). One heterozygous variant is a carrier.
    Biallelic,
}

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
    pub report_if: ReportIf,
}

/// The rows of the contraindication table; lines starting with `#` are
/// comments. A missing column, an empty cell or a `Report_if` other than
/// `any` or `biallelic` is an error, never a row read as half a rule.
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
            let report_if = match cells[7] {
                "any" => ReportIf::Any,
                "biallelic" => ReportIf::Biallelic,
                other => {
                    return Err(PyErr::value(format!(
                        "{path}: Report_if is {}, not any or biallelic: {}",
                        repr(other),
                        repr_row(r)
                    )));
                }
            };
            let own = |i: usize| cells[i].to_owned();
            Ok(Contraindication {
                gene: own(0),
                drug: own(1),
                names: own(2),
                cpic_level: own(3),
                note: own(4),
                pair_source: own(5),
                label_source: own(6),
                report_if,
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
    fn medication_lines_name_every_drug() {
        let p = tmp(
            "meds.txt",
            "Paroxetine 20 mg daily\n, paroxetine\n;\n# a comment\n\nlopinavir/ritonavir 400/100 mg\r\nValproate  # evening\nnaltrexone/bupropion 8/90 mg\nOlanzapine + fluoxetine; olanzapine\nlopinavir and ritonavir\n20 mg\nContrave (naltrexone/bupropion) 8/90 mg\nEpilim (sodium valproate) 500 mg\nDepakote(divalproex sodium) [delayed-release]\nparoxetine.\nolanzapine-fluoxetine 6/25 mg\nCo-trimoxazole 960 mg\nL-dopa 100 mg\n5-HTP 100 mg\n( )\n",
        );
        let got: Vec<(String, Vec<String>)> = read_medications(&p)
            .unwrap()
            .into_iter()
            .map(|m| (m.text, m.names))
            .collect();
        let want: [(&str, &[&str]); 16] = [
            ("Paroxetine 20 mg daily", &["paroxetine"]),
            (", paroxetine", &["paroxetine"]),
            (
                "lopinavir/ritonavir 400/100 mg",
                &["lopinavir", "ritonavir"],
            ),
            ("Valproate", &["valproate"]),
            ("naltrexone/bupropion 8/90 mg", &["naltrexone", "bupropion"]),
            (
                "Olanzapine + fluoxetine; olanzapine",
                &["olanzapine", "fluoxetine"],
            ),
            ("lopinavir and ritonavir", &["lopinavir", "ritonavir"]),
            ("20 mg", &["20"]),
            // Brackets split a line, and punctuation at a word's ends is
            // dropped, so a brand name with its generic in brackets names
            // both, and a line of only brackets names nothing.
            (
                "Contrave (naltrexone/bupropion) 8/90 mg",
                &["contrave", "naltrexone", "bupropion"],
            ),
            (
                "Epilim (sodium valproate) 500 mg",
                &["epilim", "sodium", "valproate"],
            ),
            (
                "Depakote(divalproex sodium) [delayed-release]",
                &[
                    "depakote",
                    "divalproex",
                    "sodium",
                    "delayed-release",
                    "delayed",
                    "release",
                ],
            ),
            ("paroxetine.", &["paroxetine"]),
            // A hyphenated word names itself, then each part that can name
            // a drug: never a prefix, a lone letter or a number.
            (
                "olanzapine-fluoxetine 6/25 mg",
                &["olanzapine-fluoxetine", "olanzapine", "fluoxetine"],
            ),
            ("Co-trimoxazole 960 mg", &["co-trimoxazole", "trimoxazole"]),
            ("L-dopa 100 mg", &["l-dopa", "dopa"]),
            ("5-HTP 100 mg", &["5-htp"]),
        ];
        let want: Vec<(String, Vec<String>)> = want
            .iter()
            .map(|(t, n)| ((*t).to_owned(), n.iter().map(|s| (*s).to_owned()).collect()))
            .collect();
        assert_eq!(got, want);
    }

    /// The names of `names` that are a hyphen part of another of them.
    fn part_clashes(names: &[String]) -> Vec<String> {
        names
            .iter()
            .flat_map(|n| hyphen_parts(n).filter(|p| names.iter().any(|m| m == p)))
            .map(str::to_owned)
            .collect()
    }

    #[test]
    fn no_table_name_is_a_hyphen_part_of_another() {
        // A hyphen part of a medication word is checked even when the whole
        // word is a table name; that adds no hit only while no name of the
        // tables is a part of another one.
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets/");
        let inh = read_inhibitors(&format!("{root}pgx_inhibitors.tsv")).unwrap();
        let contra = read_contraindications(&format!("{root}pgx_contraindications.tsv")).unwrap();
        let mut names: Vec<String> = inh.keys().map(|(d, _)| d.clone()).collect();
        for r in &contra {
            names.extend(r.names.split(',').map(norm_drug).filter(|n| !n.is_empty()));
        }
        assert!(names.len() > 40, "{names:?}");
        assert_eq!(part_clashes(&names), Vec::<String>::new());
        // Control: a planted clash is found.
        let planted = ["co-trimoxazole".to_owned(), "trimoxazole".to_owned()];
        assert_eq!(part_clashes(&planted), ["trimoxazole"]);
    }

    #[test]
    fn indeterminate_is_no_pypgx_call() {
        for d in ["*1/*4", " *1/*1 ", "Reference/rs9923231"] {
            assert!(pypgx_called(d), "{d}");
        }
        for d in [
            "",
            "FAILED",
            "N/A",
            "Not called",
            "Indeterminate",
            "indeterminate",
            " INDETERMINATE ",
            "*1/Indeterminate",
            "Indeterminate/Indeterminate",
            // The shapes summary.json's count refuses, refused here too.
            "failed",
            "none",
            "unknown",
            ".",
            "*1",
            "*1/",
            "*1/.",
            "*1/*2;*1/*3",
        ] {
            assert!(!pypgx_called(d), "{d}");
        }
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
        assert!(contra.iter().any(|c| c.gene == "POLG"
            && c.drug == "valproic acid"
            && c.report_if == ReportIf::Any));
        assert!(
            contra
                .iter()
                .any(|c| c.gene == "ASS1" && c.report_if == ReportIf::Biallelic)
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
                "ValueError: {p}: the header is ['Gene', 'Drug'], not ['Gene', 'Drug', 'Names', 'CPIC_level', 'Note', 'Pair_source', 'Label_source', 'Report_if']"
            )
        );
        let p = tmp(
            "contra2.tsv",
            "Gene\tDrug\tNames\tCPIC_level\tNote\tPair_source\tLabel_source\tReport_if\nPOLG\tx\ty\n",
        );
        assert_eq!(
            read_contraindications(&p).unwrap_err().to_string(),
            format!(
                "ValueError: {p}: a row has an empty or missing cell: {{'Gene': 'POLG', 'Drug': 'x', 'Names': 'y', 'CPIC_level': None, 'Note': None, 'Pair_source': None, 'Label_source': None, 'Report_if': None}}"
            )
        );
        let p = tmp(
            "contra3.tsv",
            "Gene\tDrug\tNames\tCPIC_level\tNote\tPair_source\tLabel_source\tReport_if\nASL\tv\tv\tB\tn\tp\tl\trecessive\n",
        );
        assert_eq!(
            read_contraindications(&p).unwrap_err().to_string(),
            format!(
                "ValueError: {p}: Report_if is 'recessive', not any or biallelic: {{'Gene': 'ASL', 'Drug': 'v', 'Names': 'v', 'CPIC_level': 'B', 'Note': 'n', 'Pair_source': 'p', 'Label_source': 'l', 'Report_if': 'recessive'}}"
            )
        );
        let p = tmp("hits.tsv", "nothing\there\n");
        assert_eq!(
            read_clinvar_hits(&p).unwrap_err().to_string(),
            format!("ValueError: {p}: no column alt, chrom, clnsig, geneinfo, genotype, pos, ref")
        );
    }
}
