//! `og pgx cpic-report` and `og pgx step36-outside`.

use crate::guidance::{INDETERMINATE_NOTE, drug_guidance, drug_level_phenotype, drug_lines};
use crate::pyio::{PyErr, basename, load_json};
use crate::report::{GeneCall, Parsed, parse_genes, related_drugs, text};
use crate::sections::{
    NO_CONVERSION, PHENOCONVERSION_GENES, Phenoconversion, compare, consensus_lines,
    contraindication_lines, fallback_lines, phenoconversion, pypgx_warnings, stale_lines,
    unconfirm_outside,
};
use crate::status::Status;
use crate::tables::{
    Medication, read_clinvar_hits, read_consensus, read_contraindications, read_inhibitors,
    read_medications, read_pypgx,
};
use og_io::AtomicFile;
use pycompat::csv::Writer;
use pycompat::json::Value;
use std::collections::HashSet;
use std::io::{self, Write};
use std::path::Path;

const USAGE: &str = "usage: og pgx cpic-report --sample S --report R.json --outdir DIR
                          [--pypgx S_pypgx_summary.tsv --comparison OUT.tsv]
                          [--consensus S_pgx_consensus.tsv]
                          [--medications LIST --inhibitors pgx_inhibitors.tsv]
                          [--contraindications pgx_contraindications.tsv
                           [--clinvar-hits S_clinvar_hits.tsv]]
       og pgx step36-outside --report R.json
";

/// The options of `cpic-report`. An option given as an empty string counts
/// as not given, as the Python's truth tests read it.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CpicArgs {
    pub sample: String,
    pub report: String,
    pub outdir: String,
    pub pypgx: Option<String>,
    pub comparison: Option<String>,
    pub consensus: Option<String>,
    pub medications: Option<String>,
    pub inhibitors: Option<String>,
    pub contraindications: Option<String>,
    pub clinvar_hits: Option<String>,
    /// The `Date:` line (`YYYY-MM-DD`); today in UTC when `None`. Hidden:
    /// for fixed test outputs.
    pub date: Option<String>,
}

/// argparse's refusal: usage and the error on stderr, exit 2.
fn usage_error(msg: &str) -> u8 {
    eprint!("{USAGE}");
    eprintln!("og pgx: error: {msg}");
    2
}

/// `--name VALUE` and `--name=VALUE`, the last one given winning, as
/// argparse reads them.
fn options(args: &[String], names: &[&str]) -> Result<Vec<Option<String>>, String> {
    let mut values = vec![None; names.len()];
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let (key, inline) = match a.split_once('=') {
            Some((k, v)) if k.starts_with("--") => (k, Some(v.to_owned())),
            _ => (a.as_str(), None),
        };
        let Some(k) = names.iter().position(|n| *n == key) else {
            return Err(format!("unrecognized arguments: {a}"));
        };
        let v = match inline {
            Some(v) => v,
            None => match it.next() {
                Some(v) if !v.starts_with("--") || v == "--" => v.clone(),
                _ => return Err(format!("argument {key}: expected one argument")),
            },
        };
        values[k] = Some(v);
    }
    Ok(values)
}

/// `og pgx <verb> [args]`: the exit code.
pub fn main(args: &[String]) -> u8 {
    match args.first().map(String::as_str) {
        Some("cpic-report") => {
            const NAMES: [&str; 11] = [
                "--sample",
                "--report",
                "--outdir",
                "--pypgx",
                "--comparison",
                "--consensus",
                "--medications",
                "--inhibitors",
                "--contraindications",
                "--clinvar-hits",
                "--date",
            ];
            let v = match options(&args[1..], &NAMES) {
                Ok(v) => v,
                Err(e) => return usage_error(&e),
            };
            let missing: Vec<&str> = NAMES[..3]
                .iter()
                .zip(&v)
                .filter(|(_, v)| v.is_none())
                .map(|(n, _)| *n)
                .collect();
            if !missing.is_empty() {
                return usage_error(&format!(
                    "the following arguments are required: {}",
                    missing.join(", ")
                ));
            }
            let given = |i: usize| v[i].clone().filter(|s| !s.is_empty());
            let a = CpicArgs {
                sample: v[0].clone().unwrap_or_default(),
                report: v[1].clone().unwrap_or_default(),
                outdir: v[2].clone().unwrap_or_default(),
                pypgx: given(3),
                comparison: given(4),
                consensus: given(5),
                medications: given(6),
                inhibitors: given(7),
                contraindications: given(8),
                clinvar_hits: given(9),
                date: given(10),
            };
            if a.medications.is_some() && a.inhibitors.is_none() {
                return usage_error("--medications needs --inhibitors (assets/pgx_inhibitors.tsv)");
            }
            if a.clinvar_hits.is_some() && a.contraindications.is_none() {
                return usage_error(
                    "--clinvar-hits needs --contraindications (assets/pgx_contraindications.tsv)",
                );
            }
            if let Some(d) = &a.date
                && !is_iso_date(d)
            {
                return usage_error(&format!("argument --date: not a YYYY-MM-DD date: {d}"));
            }
            cpic_report(&a)
        }
        Some("step36-outside") => match options(&args[1..], &["--report"]) {
            Ok(v) => match &v[0] {
                Some(report) => step36_outside(report),
                None => usage_error("the following arguments are required: --report"),
            },
            Err(e) => usage_error(&e),
        },
        Some(other) => usage_error(&format!("unknown pgx tool '{other}'")),
        None => usage_error("the following arguments are required: pgx tool"),
    }
}

fn is_iso_date(d: &str) -> bool {
    let b = d.as_bytes();
    b.len() == 10
        && b[4] == b'-'
        && b[7] == b'-'
        && b.iter()
            .enumerate()
            .all(|(i, c)| i == 4 || i == 7 || c.is_ascii_digit())
}

/// Today's date in UTC, `YYYY-MM-DD` (the containers the step runs in keep
/// UTC, which Python's `date.today()` read).
fn today_utc() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    let days = i64::try_from(secs / 86_400).unwrap_or(0);
    // Howard Hinnant's civil_from_days.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!("{year:04}-{month:02}-{day:02}")
}

/// Writes `body` to `path` atomically, so a killed task leaves no half file
/// a later `[ -s ]` check would trust.
fn write_atomic(path: &Path, body: &[u8]) -> io::Result<()> {
    let mut f = AtomicFile::create(path)?;
    f.write_all(body)?;
    f.commit()
}

/// The rows as Python's `csv.writer(f, delimiter="\t", lineterminator="\n")`.
fn tsv<R, I>(rows: R) -> Vec<u8>
where
    R: IntoIterator<Item = I>,
    I: IntoIterator,
    I::Item: AsRef<str>,
{
    let mut w = Writer::new(Vec::new(), '\t');
    for r in rows {
        // Writing to a Vec cannot fail.
        let _ = w.write_row(r);
    }
    w.into_inner()
}

/// Exit 1 with `what` and the error, for a write that failed.
fn write_failed(path: &Path, e: &io::Error) -> u8 {
    eprintln!("ERROR: could not write {}: {e}", path.display());
    1
}

/// The report that says the parse failed, never an all-clear, and an empty
/// phenotypes table.
fn write_failure(a: &CpicArgs, date: &str, why: &str) -> io::Result<()> {
    let dir = Path::new(&a.outdir);
    let body = format!(
        "Pharmacogenomic Drug Recommendations\n=====================================\nSample: {}\nDate: {date}\n\n\
         !! PARSING FAILED -- RESULTS CANNOT BE TRUSTED !!\n{why}\n\
         This is NOT a clean result. Read the PharmCAT HTML report directly,\n\
         and report this as a pipeline bug (likely a PharmCAT output-format change).\n",
        a.sample
    );
    write_atomic(
        &dir.join(format!("{}_cpic_recommendations.txt", a.sample)),
        body.as_bytes(),
    )?;
    write_atomic(
        &dir.join(format!("{}_phenotypes.tsv", a.sample)),
        b"Gene\tDiplotype\tPhenotype\tStatus\n",
    )
}

/// `-` repeated 72 times, the rule under each heading.
fn rule() -> String {
    "-".repeat(72)
}

/// The medication list and its phenoconversion, or why there is none.
struct Medications {
    meds: Vec<Medication>,
    error: Option<PyErr>,
    conversion: Phenoconversion,
}

fn medications(
    a: &CpicArgs,
    calls: &[GeneCall],
    guidance: &crate::guidance::GuidanceMap,
    related: &std::collections::HashMap<String, Vec<String>>,
) -> Medications {
    let empty = || Phenoconversion {
        column: std::collections::HashMap::new(),
        lines: Vec::new(),
    };
    let (Some(list), Some(inhibitors)) = (&a.medications, &a.inhibitors) else {
        return Medications {
            meds: Vec::new(),
            error: None,
            conversion: empty(),
        };
    };
    let run = || -> Result<(Vec<Medication>, Phenoconversion), PyErr> {
        let meds = read_medications(list)?;
        let inhibitors = read_inhibitors(inhibitors)?;
        let conversion = phenoconversion(calls, &meds, &inhibitors, guidance, related)?;
        Ok((meds, conversion))
    };
    match run() {
        Ok((meds, conversion)) => Medications {
            meds,
            error: None,
            conversion,
        },
        Err(e) => {
            // Opt-in and beside the report: say why it is missing, never fail the step.
            eprintln!(
                "WARNING: could not read the medication list or the inhibitor table ({e}); no phenotype is converted"
            );
            Medications {
                meds: Vec::new(),
                error: Some(e),
                conversion: empty(),
            }
        }
    }
}

/// The `Drug Contraindications From Disease Variants` section: the
/// contraindication table against step 06's ClinVar hits. Never fails the
/// step: an unreadable file is said in the section.
fn contraindication_section(a: &CpicArgs, table_path: &str) -> Vec<String> {
    let mut out = vec![
        "Drug Contraindications From Disease Variants (ClinVar hits, step 06):".to_owned(),
        rule(),
        String::new(),
    ];
    let table = match read_contraindications(table_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("WARNING: could not read the contraindication table ({e})");
            out.extend([
                format!("  Could not read the contraindication table ({e})."),
                String::new(),
            ]);
            return out;
        }
    };
    let mut genes: Vec<&str> = Vec::new();
    for r in &table {
        if !genes.contains(&r.gene.as_str()) {
            genes.push(&r.gene);
        }
    }
    let genes = genes.join(", ");
    let Some(hits_path) = &a.clinvar_hits else {
        out.extend([
            "  Not checked: no ClinVar hits table (step 06, or clinvar in --tools, did not run)."
                .to_owned(),
            format!("  Genes it would check: {genes}."),
            String::new(),
        ]);
        return out;
    };
    let hits = match read_clinvar_hits(hits_path) {
        Ok(h) => h,
        Err(e) => {
            eprintln!("WARNING: could not read the ClinVar hits table ({e})");
            out.extend([
                format!("  Could not read the ClinVar hits table ({e})."),
                format!("  Not checked: {genes}. Run step 06 again, then this step."),
                String::new(),
            ]);
            return out;
        }
    };
    // The phenoconversion section already says when the list is unreadable.
    let med_names: HashSet<String> = a
        .medications
        .as_deref()
        .and_then(|m| read_medications(m).ok())
        .map(|m| m.into_iter().map(|m| m.name).collect())
        .unwrap_or_default();
    let found = contraindication_lines(&hits, &table, &med_names);
    if found.is_empty() {
        out.extend([
            format!("  No pathogenic or likely pathogenic ClinVar hit in {genes}."),
            String::new(),
        ]);
        return out;
    }
    out.extend(found);
    out.extend(
        [
            "  Whether a warning applies depends on the disorder the gene causes and on how many",
            "  copies carry the variant: one copy in a recessive gene usually makes a carrier, not a",
            "  patient. Review it with your clinician before taking or stopping the drug.",
            "",
        ]
        .map(str::to_owned),
    );
    out
}

/// `og pgx cpic-report`: `DIR/S_phenotypes.tsv`, `DIR/S_cpic_recommendations.txt`
/// and, with `--pypgx` and `--comparison`, the PharmCAT/pypgx comparison.
/// Exit 1 when the report cannot be read or yields no gene, after writing a
/// report that says so: a report that parses to nothing is a format change,
/// never an all-clear.
pub fn cpic_report(a: &CpicArgs) -> u8 {
    let date = a.date.clone().unwrap_or_else(today_utc);
    let dir = Path::new(&a.outdir);
    if let Err(e) = std::fs::create_dir_all(dir) {
        eprintln!("ERROR: could not create {}: {e}", a.outdir);
        return 1;
    }
    let data = match load_json(&a.report) {
        Ok(d) => d,
        Err(e) => {
            let why = format!("Could not read {}: {}", a.report, e.msg);
            if let Err(w) = write_failure(a, &date, &why) {
                return write_failed(dir, &w);
            }
            eprintln!(
                "ERROR: could not read PharmCAT report {}: {}",
                a.report, e.msg
            );
            return 1;
        }
    };
    let mut calls = match parse_genes(&data) {
        Parsed::Ok(c) => c,
        failed => {
            let (status, why) = if failed == Parsed::UnknownFormat {
                (
                    "UNKNOWN_FORMAT",
                    "The report has no gene results in any known layout.",
                )
            } else {
                (
                    "PARSE_EMPTY",
                    "The report's gene section was recognised but no gene diplotype could be read.",
                )
            };
            if let Err(w) = write_failure(a, &date, why) {
                return write_failed(dir, &w);
            }
            eprintln!(
                "ERROR: {status}: no gene parsed from {}; refusing to write an all-clear report",
                a.report
            );
            return 1;
        }
    };
    let guidance = drug_guidance(&data);
    let related = related_drugs(&data);

    let mut consensus = Vec::new();
    let mut consensus_error = None;
    if let Some(path) = &a.consensus {
        match read_consensus(path) {
            Ok(c) => consensus = c,
            Err(e) => {
                eprintln!("WARNING: could not read the consensus table {path} ({e})");
                consensus_error = Some(e);
            }
        }
    }
    let held_back: HashSet<String> = consensus
        .iter()
        .filter(|r| !r.passed)
        .map(|r| r.gene.clone())
        .collect();
    let stale = if a.consensus.is_some() {
        unconfirm_outside(&mut calls, &consensus)
    } else {
        Vec::new()
    };
    for c in &mut calls {
        if let Some(p) = drug_level_phenotype(c, &guidance) {
            c.phenotype = p;
            c.note = Some(INDETERMINATE_NOTE);
        }
    }

    let mut pypgx = std::collections::HashMap::new();
    let mut pypgx_error = None;
    if let Some(path) = &a.pypgx {
        match read_pypgx(path) {
            Ok(p) => pypgx = p,
            Err(e) => {
                // The CPIC report does not depend on pypgx: write it, say
                // why the comparison is missing, and write no comparison.
                eprintln!(
                    "WARNING: could not read the pypgx summary {path} ({e}); no PharmCAT/pypgx comparison is made"
                );
                pypgx_error = Some(e);
            }
        }
    }
    let meds = medications(a, &calls, &guidance, &related);

    // The phenotypes table: the converted phenotype beside the genotype
    // one, never in its place.
    let table = if meds.meds.is_empty() {
        let head = vec!["Gene", "Diplotype", "Phenotype", "Status"]
            .into_iter()
            .map(str::to_owned)
            .collect::<Vec<_>>();
        tsv(std::iter::once(head).chain(calls.iter().map(|c| {
            vec![
                c.gene.clone(),
                c.diplotype.clone(),
                c.shown(),
                c.status().to_string(),
            ]
        })))
    } else {
        let head = [
            "Gene",
            "Diplotype",
            "Phenotype",
            "Status",
            "Phenoconverted_phenotype",
        ]
        .map(str::to_owned)
        .to_vec();
        tsv(std::iter::once(head).chain(calls.iter().map(|c| {
            let conv = meds
                .conversion
                .column
                .get(c.gene.as_str())
                .map_or(NO_CONVERSION, String::as_str);
            vec![
                c.gene.clone(),
                c.diplotype.clone(),
                c.shown(),
                c.status().to_string(),
                conv.to_owned(),
            ]
        })))
    };
    let pheno_path = dir.join(format!("{}_phenotypes.tsv", a.sample));
    if let Err(e) = write_atomic(&pheno_path, &table) {
        return write_failed(&pheno_path, &e);
    }

    let version = match crate::report::field(&data, "pharmcatVersion") {
        Some(v) => text(v),
        None => String::new(),
    };
    let source = format!(
        "Source: PharmCAT {version} report ({})",
        basename(&a.report)
    );
    let mut lines: Vec<String> = vec![
        "Pharmacogenomic Drug Recommendations".to_owned(),
        "=====================================".to_owned(),
        format!("Sample: {}", a.sample),
        format!("Date: {date}"),
        pycompat::text::rstrip(&source, None).to_owned(),
        String::new(),
        "Gene Results:".to_owned(),
        rule(),
        format!("{:<12} {:<30} Phenotype", "Gene", "Diplotype"),
    ];
    for c in &calls {
        let outside = if c.outside { "  [outside call]" } else { "" };
        lines.push(format!(
            "{:<12} {:<30} {}{outside}",
            c.gene,
            c.diplotype,
            c.shown()
        ));
    }
    let genes = PHENOCONVERSION_GENES.join(", ");
    if let Some(list) = &a.medications {
        lines.extend([
            String::new(),
            format!("Medications and Phenoconversion ({genes}):"),
            rule(),
            String::new(),
        ]);
        if let Some(e) = &meds.error {
            lines.push(format!(
                "  Could not read the medication list or the inhibitor table ({e})."
            ));
            lines.push("  No phenotype is converted.".to_owned());
        } else if meds.meds.is_empty() {
            lines
                .push("  The medication list names no drug: no phenotype is converted.".to_owned());
        } else {
            let n = meds.meds.len();
            lines.push(format!(
                "  Medication list: {} ({n} line{}).",
                basename(list),
                if n == 1 { "" } else { "s" }
            ));
            if meds.conversion.lines.is_empty() {
                lines.push(format!(
                    "  No drug on the list is an inhibitor of {genes} in the inhibitor table."
                ));
            } else {
                lines.extend(meds.conversion.lines.iter().cloned());
            }
            lines.push("  The genotype phenotype in the table above is unchanged. A phenoconverted phenotype".to_owned());
            lines.push("  holds while the inhibitor is taken. This is an aid for your clinician, not advice.".to_owned());
        }
    }
    lines.extend([
        String::new(),
        "Affected Medications:".to_owned(),
        rule(),
        String::new(),
    ]);
    let by_status = |s: Status| calls.iter().filter(move |c| c.status() == s);
    let contra = a
        .contraindications
        .as_deref()
        .map(|t| contraindication_section(a, t))
        .unwrap_or_default();
    let listed: Vec<&GeneCall> = by_status(Status::NonNormal).collect();
    for c in &listed {
        lines.push(format!("  {} -- {}:", c.gene, c.shown()));
        lines.push(format!("    Diplotype: {}", c.diplotype));
        lines.extend(drug_lines(c, &guidance, &related));
        lines.push(
            "    Action: Consult CPIC guidelines at https://cpicpgx.org/guidelines/".to_owned(),
        );
        lines.push(String::new());
    }
    if listed.is_empty() {
        lines.extend([
            "  No gene with a non-normal phenotype.".to_owned(),
            String::new(),
        ]);
    }
    lines.extend(contra);

    let ambiguous: Vec<&GeneCall> = by_status(Status::Ambiguous).collect();
    if !ambiguous.is_empty() {
        lines.extend([
            "Genes With More Than One Possible Result:".to_owned(),
            rule(),
            String::new(),
        ]);
        for c in &ambiguous {
            let more = if c.phenotypes.len() > 5 { ", ..." } else { "" };
            lines.push(format!(
                "  {} -- {} possible diplotypes with different phenotypes ({}{more}).",
                c.gene,
                c.labels.len(),
                c.phenotypes[..c.phenotypes.len().min(5)].join(", ")
            ));
        }
        lines.extend(
            [
                "  The data cannot tell them apart (often positions missing from the VCF), so no",
                "  drug guidance is given for them here. See the PharmCAT HTML report.",
                "",
            ]
            .map(str::to_owned),
        );
    }

    let unclassified: Vec<&GeneCall> = by_status(Status::Unclassified).collect();
    if !unclassified.is_empty() {
        lines.extend([
            "Called Genes Without a Function Phenotype:".to_owned(),
            rule(),
            String::new(),
        ]);
        for c in &unclassified {
            lines.push(format!("  {} -- {} ({})", c.gene, c.shown(), c.diplotype));
        }
        lines.extend(
            [
                "  PharmCAT called these genes but its phenotype (n/a, a reference genotype,",
                "  'no phenotype assigned', an uncertain result) names no changed function, so",
                "  they are not counted as non-normal and get no drug list here. A variant allele",
                "  among them can still matter: see the PharmCAT HTML report.",
                "",
            ]
            .map(str::to_owned),
        );
    }

    if a.consensus.is_some() {
        lines.extend([
            "Calls From Other Tools (outside calls, step 36):".to_owned(),
            rule(),
            String::new(),
        ]);
        if let Some(e) = &consensus_error {
            lines.push(format!("  Could not read the consensus table ({e})."));
            lines.push("  Run step 36 again, then steps 07 and 27.".to_owned());
        } else {
            let cl = consensus_lines(&consensus, &calls, &guidance, &related);
            if cl.is_empty() {
                lines.push("  The table lists no gene.".to_owned());
            } else {
                lines.extend(cl);
            }
        }
        let stale: Vec<&GeneCall> = stale.iter().map(|&i| &calls[i]).collect();
        lines.extend(stale_lines(&stale, &guidance, &related));
        lines.push(String::new());
    }

    lines.extend(["Uncallable Genes:".to_owned(), rule(), String::new()]);
    let uncalled: Vec<&GeneCall> = by_status(Status::NotCalled).collect();
    for c in &uncalled {
        lines.push(format!(
            "  {} -- {} (not callable from available data)",
            c.gene, c.phenotype
        ));
    }
    if uncalled.is_empty() {
        lines.push("  None -- all genes were successfully called.".to_owned());
    }
    lines.push(String::new());

    let fallback = fallback_lines(&calls, &data);
    if !fallback.is_empty() {
        lines.extend([
            "Single-Position Fallback (not a PharmCAT call):".to_owned(),
            rule(),
            String::new(),
        ]);
        lines.extend(fallback);
        lines.extend(
            [
                "  Read from the genotype PharmCAT lists at the position in its report.json. Confirm",
                "  it before acting on it: PharmCAT itself made no call for the gene.",
                "",
            ]
            .map(str::to_owned),
        );
    }

    if a.pypgx.is_some() {
        lines.extend(["PharmCAT and pypgx:".to_owned(), rule(), String::new()]);
        if let Some(e) = &pypgx_error {
            lines.push(format!(
                "  Could not read the pypgx summary ({e}); no comparison was made."
            ));
            lines.push("  Run step 32 again, then this step.".to_owned());
        } else {
            let warn = pypgx_warnings(&calls, &pypgx, &guidance, &related, &held_back);
            if warn.is_empty() {
                let tail = if held_back.iter().any(|g| pypgx.contains_key(g)) {
                    let mut held: Vec<&str> = held_back.iter().map(String::as_str).collect();
                    held.sort_unstable();
                    format!(
                        " outside the ones step 36 held back ({}: see above).",
                        held.join(", ")
                    )
                } else {
                    ".".to_owned()
                };
                lines.push(format!(
                    "  No gene that PharmCAT could not call has a pypgx call{tail}"
                ));
            } else {
                lines.extend(warn);
            }
        }
        lines.push(String::new());
    }

    lines.extend(
        [
            "NOTE: Uncallable genes may lack coverage, have complex structural",
            "variants, or require data not present in your VCF. Their absence",
            "from the recommendations section does NOT mean normal function.",
            "",
        ]
        .map(str::to_owned),
    );
    lines.push(rule());
    lines.extend(
        [
            "DISCLAIMER: These recommendations are based on CPIC clinical",
            "guidelines. Always consult your healthcare provider before",
            "making any medication changes. This is NOT medical advice.",
        ]
        .map(str::to_owned),
    );
    lines.push(rule());
    let mut body = lines.join("\n");
    body.push('\n');
    let rec_path = dir.join(format!("{}_cpic_recommendations.txt", a.sample));
    if let Err(e) = write_atomic(&rec_path, body.as_bytes()) {
        return write_failed(&rec_path, &e);
    }

    if let (Some(_), Some(comparison)) = (&a.pypgx, &a.comparison) {
        let path = Path::new(comparison);
        if pypgx_error.is_some() {
            // An earlier run's comparison must not stand in for this one.
            if path.exists()
                && let Err(e) = std::fs::remove_file(path)
            {
                return write_failed(path, &e);
            }
        } else {
            let rows = compare(&calls, &pypgx);
            let head = [
                "Gene",
                "PharmCAT_diplotype",
                "pypgx_diplotype",
                "Match",
                "Called_by",
                "PharmCAT_phenotype",
                "pypgx_phenotype",
            ];
            let body = tsv(std::iter::once(head).chain(rows.iter().map(|r| r.cells())));
            if let Err(e) = write_atomic(path, &body) {
                return write_failed(path, &e);
            }
            let same = rows.iter().filter(|r| r.matched == "Yes").count();
            println!(
                "Comparison written: {comparison} (concordant {same}, other {})",
                rows.len() - same
            );
        }
    }
    println!(
        "Genes parsed: {}; non-normal: {}; ambiguous: {}; unclassified: {}; not called: {}",
        calls.len(),
        listed.len(),
        ambiguous.len(),
        unclassified.len(),
        uncalled.len()
    );
    0
}

/// `og pgx step36-outside`: prints, one per line, the step 36 genes
/// (CYP2D6, HLA-A, HLA-B) PharmCAT's report calls from an outside call: the
/// ones a missing consensus table leaves unconfirmed. Step 27 runs it to
/// decide whether to ask for step 36 again.
pub fn step36_outside(report: &str) -> u8 {
    let data: Value = match load_json(report) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("ERROR: could not read PharmCAT report {report}: {}", e.msg);
            return 1;
        }
    };
    if let Parsed::Ok(calls) = parse_genes(&data) {
        for c in calls.iter().filter(|c| c.step36_outside()) {
            println!("{}", c.gene);
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn options_read_as_argparse_reads_them() {
        let args: Vec<String> = ["--b=2", "--a", "1", "--a", "3"]
            .map(str::to_owned)
            .to_vec();
        assert_eq!(
            options(&args, &["--a", "--b"]).unwrap(),
            [Some("3".to_owned()), Some("2".to_owned())]
        );
        let args: Vec<String> = ["--a"].map(str::to_owned).to_vec();
        assert!(options(&args, &["--a"]).is_err());
        let args: Vec<String> = ["--c", "x"].map(str::to_owned).to_vec();
        assert!(options(&args, &["--a"]).is_err());
    }

    #[test]
    fn dates() {
        assert!(is_iso_date("2026-11-10"));
        assert!(!is_iso_date("2026-1-10"));
        assert!(is_iso_date(&today_utc()), "{}", today_utc());
    }
}
