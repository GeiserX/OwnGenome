//! `og exomiser prepare` and `og exomiser report` (step 48): the
//! symptom-driven ranking of Exomiser.
//!
//! ```text
//! og exomiser prepare --hpo HP:0001156,HP:0001363 --sample S --vcf S.vcf.gz
//!                     --assembly hg38 --out-dir DIR --analysis S_analysis.yml
//!                     --terms S_exomiser_terms.tsv [--obo hp.obo]
//!                     [--pedigree pedigree.tsv --ped S_exomiser.ped
//!                      --affected a,b --unaffected c,d]
//! og exomiser report  --sample S --variants S_exomiser.variants.tsv
//!                     --terms S_exomiser_terms.tsv
//!                     --out S_exomiser_top.tsv --text S_exomiser.txt [--top 10]
//! ```
//!
//! `prepare` checks the HPO terms (HP: and seven digits; with `--obo`, each
//! must be a current term of the ontology Exomiser's data holds: an
//! alternative id is read as its term, an obsolete one stops the step with
//! its replacement) and writes the analysis file (JSON, which YAML reads):
//! the VCF, the terms, Exomiser's exome preset and the output options. With
//! `--pedigree` (step 38's pedigree.tsv) the VCF is the family's joint call
//! and `--ped` gets the PED file of its members: the proband and the
//! `--affected` are affected, the `--unaffected` are not, everyone else is
//! unknown.
//!
//! `report` writes the top `--top` ranks (a gene under one mode of
//! inheritance is one rank) with each variant that contributes to the
//! gene's score, in a table and in words. In a family run the PED beside
//! the variants names each member's genotype.
//!
//! What changes on purpose: the summary `report` prints on stdout is the
//! "Symptoms searched" line and one line per rank kept; the Python counted
//! variant rows, so a rank with several variants also printed the closing
//! note (stdout only; the text file is the same).

use crate::{Error, Outcome};
use indexmap::IndexMap;
use og_filters::args::Spec;
use og_io::Utf8;
use og_tsv::{Options, Quoting};
use pycompat::json::{Dump, Value};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;

pub const TERMS_HEADER: [&str; 3] = ["hpo_id", "label", "given_as"];
pub const OUT_HEADER: [&str; 16] = [
    "rank",
    "gene",
    "moi",
    "exomiser_score",
    "phenotype_score",
    "variant_score",
    "p_value",
    "variant",
    "genotype",
    "consequence",
    "hgvs",
    "acmg",
    "acmg_evidence",
    "disease",
    "clinvar",
    "max_frequency",
];
/// Columns of Exomiser's variants TSV (TsvVariantResultsWriter, 15.x).
pub const NEED: [&str; 23] = [
    "#RANK",
    "GENE_SYMBOL",
    "MOI",
    "P-VALUE",
    "EXOMISER_GENE_COMBINED_SCORE",
    "EXOMISER_GENE_PHENO_SCORE",
    "EXOMISER_GENE_VARIANT_SCORE",
    "CONTRIBUTING_VARIANT",
    "CONTIG",
    "START",
    "REF",
    "ALT",
    "GENOTYPE",
    "FUNCTIONAL_CLASS",
    "HGVS",
    "EXOMISER_ACMG_CLASSIFICATION",
    "EXOMISER_ACMG_EVIDENCE",
    "EXOMISER_ACMG_DISEASE_ID",
    "EXOMISER_ACMG_DISEASE_NAME",
    "CLINVAR_PRIMARY_INTERPRETATION",
    "CLINVAR_STAR_RATING",
    "MAX_FREQ_SOURCE",
    "MAX_FREQ",
];

fn moi_words(m: &str) -> &str {
    match m {
        "AD" => "autosomal dominant",
        "AR" => "autosomal recessive",
        "XD" => "X-linked dominant",
        "XR" => "X-linked recessive",
        "MT" => "mitochondrial",
        "ANY" => "any mode of inheritance",
        o => o,
    }
}

/// Exomiser's exome preset (examples/preset-exome-analysis.yml of 15.1.0).
const INHERITANCE: [(&str, f64); 7] = [
    ("AUTOSOMAL_DOMINANT", 0.1),
    ("AUTOSOMAL_RECESSIVE_HOM_ALT", 0.1),
    ("AUTOSOMAL_RECESSIVE_COMP_HET", 2.0),
    ("X_DOMINANT", 0.1),
    ("X_RECESSIVE_HOM_ALT", 0.1),
    ("X_RECESSIVE_COMP_HET", 2.0),
    ("MITOCHONDRIAL", 0.2),
];
const FREQUENCY_SOURCES: [&str; 11] = [
    "UK10K",
    "GNOMAD_E_AFR",
    "GNOMAD_E_AMR",
    "GNOMAD_E_EAS",
    "GNOMAD_E_NFE",
    "GNOMAD_E_SAS",
    "GNOMAD_G_AFR",
    "GNOMAD_G_AMR",
    "GNOMAD_G_EAS",
    "GNOMAD_G_NFE",
    "GNOMAD_G_SAS",
];
const PATHOGENICITY_SOURCES: [&str; 4] = ["REVEL", "MVP", "ALPHA_MISSENSE", "SPLICE_AI"];
const REMOVED_EFFECTS: [&str; 11] = [
    "FIVE_PRIME_UTR_EXON_VARIANT",
    "FIVE_PRIME_UTR_INTRON_VARIANT",
    "THREE_PRIME_UTR_EXON_VARIANT",
    "THREE_PRIME_UTR_INTRON_VARIANT",
    "NON_CODING_TRANSCRIPT_EXON_VARIANT",
    "NON_CODING_TRANSCRIPT_INTRON_VARIANT",
    "CODING_TRANSCRIPT_INTRON_VARIANT",
    "UPSTREAM_GENE_VARIANT",
    "DOWNSTREAM_GENE_VARIANT",
    "INTERGENIC_VARIANT",
    "REGULATORY_REGION_VARIANT",
];
/// Genes written to Exomiser's own files: the report shows `--top` of them.
const OUTPUT_GENES: i128 = 50;

pub const NOTE: &str = "A ranking to bring to a clinical geneticist, not a diagnosis. A high score means a rare variant, predicted to be damaging, in a gene linked to symptoms like the ones given; it does not show that the variant causes them. A low score or an empty list does not rule a genetic cause out: this search reads coding and splice variants against the terms given, and not structural variants, repeat expansions, mitochondrial heteroplasmy or variants outside genes.";

/// `re.split(r"[,\s]+", text)` without the empty items.
fn split_list(text: &str) -> Vec<&str> {
    text.split(|c: char| c == ',' || pycompat::text::is_space(c))
        .filter(|x| !x.is_empty())
        .collect()
}

/// `HP:` and seven ASCII digits.
fn is_term(t: &str) -> bool {
    t.strip_prefix("HP:")
        .is_some_and(|d| d.len() == 7 && d.bytes().all(|b| b.is_ascii_digit()))
}

/// An HPO term of the OBO file.
#[derive(Debug, Clone, Default)]
struct Term {
    name: String,
    obsolete: bool,
    replaced_by: String,
    consider: Vec<String>,
}

/// The terms of an OBO file by id, and its alternative ids (alternative id to term id).
type Obo = (HashMap<String, Term>, HashMap<String, String>);

/// Reads an OBO file's `[Term]` stanzas.
fn read_obo(path: &str) -> Result<Obo, Error> {
    #[derive(Default)]
    struct Open {
        id: Option<String>,
        term: Term,
        alt: Vec<String>,
    }
    let (mut terms, mut alt) = (HashMap::new(), HashMap::new());
    let mut close = |cur: Option<Open>| {
        if let Some(c) = cur
            && let Some(id) = c.id.filter(|i| i.starts_with("HP:"))
        {
            for a in c.alt {
                alt.insert(a, id.clone());
            }
            terms.insert(id, c.term);
        }
    };
    let mut cur: Option<Open> = None;
    for raw in crate::read_lines(path, Utf8::Replace)? {
        let line = raw.strip_suffix('\n').unwrap_or(&raw);
        if line.starts_with('[') {
            close(cur.take());
            cur = (line == "[Term]").then(Open::default);
            continue;
        }
        let Some(c) = cur.as_mut() else { continue };
        let Some((k, v)) = line.split_once(": ") else {
            continue;
        };
        let v = crate::strip(v.split(" ! ").next().unwrap_or(v)).to_owned();
        match k {
            "id" => c.id = Some(v),
            "name" => c.term.name = v,
            "replaced_by" => c.term.replaced_by = v,
            "is_obsolete" => c.term.obsolete = v == "true",
            "alt_id" => c.alt.push(v),
            "consider" => c.term.consider.push(v),
            _ => {}
        }
    }
    close(cur);
    Ok((terms, alt))
}

/// `(id, label, given_as)` of the given terms, in order, once each.
fn check_terms(given: &[&str], obo: Option<&str>) -> Result<Vec<[String; 3]>, Error> {
    if given.is_empty() {
        return Err(Error("no HPO term given: name the symptoms as HPO ids, e.g. HPO=HP:0001250,HP:0001263 (find them at https://hpo.jax.org)".into()));
    }
    let bad: Vec<&str> = given.iter().copied().filter(|t| !is_term(t)).collect();
    if !bad.is_empty() {
        return Err(Error(format!(
            "not an HPO id: {}. An HPO id is HP: and seven digits, e.g. HP:0001250 (Seizure)",
            bad.join(", ")
        )));
    }
    let (terms, alt) = match obo {
        Some(p) => read_obo(p)?,
        None => Default::default(),
    };
    let (mut out, mut seen, mut problems) = (Vec::new(), HashSet::new(), Vec::new());
    for &t in given {
        let tid = alt.get(t).map_or(t, String::as_str);
        let mut name = String::new();
        if let Some(path) = obo {
            let Some(term) = terms.get(tid) else {
                problems.push(format!("{t} is not a term of the HPO release in {path}"));
                continue;
            };
            if term.obsolete {
                let hint = if !term.replaced_by.is_empty() {
                    format!(" use {} instead", term.replaced_by)
                } else if !term.consider.is_empty() {
                    format!(" consider {}", term.consider.join(", "))
                } else {
                    " pick a current term".into()
                };
                problems.push(format!(
                    "{t} ({}) is obsolete in this HPO release;{hint}",
                    term.name
                ));
                continue;
            }
            name.clone_from(&term.name);
        }
        if seen.insert(tid.to_owned()) {
            out.push([tid.to_owned(), name, t.to_owned()]);
        }
    }
    if !problems.is_empty() {
        return Err(Error(
            problems.join("; ") + ". Look the terms up at https://hpo.jax.org",
        ));
    }
    Ok(out)
}

/// The sample names of a VCF's `#CHROM` line.
fn vcf_samples(path: &str) -> Result<Vec<String>, Error> {
    let file = og_io::open(path).map_err(|e| Error::io(&e, path))?;
    let mut lines = og_io::Lines::new(file, Utf8::Replace);
    let mut line = String::new();
    while lines
        .read_line(&mut line)
        .map_err(|e| Error::io(&e, path))?
    {
        if line.starts_with("#CHROM") {
            let l = line.strip_suffix('\n').unwrap_or(&line);
            return Ok(l.split('\t').skip(9).map(str::to_owned).collect());
        }
        if !line.starts_with('#') {
            break;
        }
    }
    Err(Error(format!("{path} has no #CHROM header line")))
}

/// Rows of step 38's pedigree TSV by sample, each a column-to-value map.
fn read_pedigree(path: &str) -> Result<HashMap<String, HashMap<String, String>>, Error> {
    let mut rows = HashMap::new();
    let mut header: Option<Vec<String>> = None;
    for raw in crate::lines(path, Utf8::Replace)? {
        let raw = raw?;
        let line = raw.trim_end_matches(['\r', '\n']);
        if crate::strip(line).is_empty() {
            continue;
        }
        let Some(h) = &header else {
            let h: Vec<String> = line
                .trim_start_matches('#')
                .split('\t')
                .map(|c| pycompat::text::lower(crate::strip(c)))
                .collect();
            if !h.iter().any(|c| c == "sample") {
                return Err(Error(format!(
                    "{path}: the first line must be a header with a 'sample' column"
                )));
            }
            header = Some(h);
            continue;
        };
        if line.starts_with('#') {
            continue;
        }
        let mut vals: Vec<String> = line
            .split('\t')
            .map(|v| crate::strip(v).to_owned())
            .collect();
        if vals.len() < h.len() {
            vals.resize(h.len(), String::new());
        }
        let r: HashMap<String, String> = h.iter().cloned().zip(vals).collect();
        rows.insert(r.get("sample").cloned().unwrap_or_default(), r);
    }
    if header.is_none() {
        return Err(Error(format!("{path} is empty")));
    }
    Ok(rows)
}

/// The PED file Exomiser reads: one line per member of the VCF.
fn ped_text(
    pedigree: &str,
    members: &[String],
    proband: &str,
    affected: &[&str],
    unaffected: &[&str],
) -> Result<String, Error> {
    let rows = read_pedigree(pedigree)?;
    let has = |n: &str| members.iter().any(|m| m == n);
    let missing: Vec<&str> = members
        .iter()
        .map(String::as_str)
        .filter(|m| !rows.contains_key(*m))
        .collect();
    if !missing.is_empty() {
        return Err(Error(format!(
            "the family VCF holds {}, which {pedigree} does not name",
            missing.join(", ")
        )));
    }
    if !has(proband) {
        return Err(Error(format!(
            "{proband} is not a sample of the family VCF ({})",
            members.join(", ")
        )));
    }
    for (label, names) in [("--affected", affected), ("--unaffected", unaffected)] {
        let stray: Vec<&str> = names.iter().copied().filter(|n| !has(n)).collect();
        if !stray.is_empty() {
            return Err(Error(format!(
                "{label} names {}, not a member of the family VCF ({})",
                stray.join(", "),
                members.join(", ")
            )));
        }
    }
    let mut both: Vec<&str> = affected
        .iter()
        .copied()
        .filter(|a| unaffected.contains(a))
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    both.sort_unstable();
    if !both.is_empty() {
        return Err(Error(format!(
            "{} is named both affected and unaffected",
            both.join(", ")
        )));
    }
    if unaffected.contains(&proband) {
        return Err(Error(format!(
            "{proband} is the proband, the person the symptoms are of: it cannot be unaffected"
        )));
    }
    let mut out = String::new();
    for m in members {
        let r = &rows[m];
        let get = |k: &str| r.get(k).map_or("", String::as_str);
        let parent = |k: &str| if has(get(k)) { get(k) } else { "0" };
        let sex = match pycompat::text::lower(get("sex")).as_str() {
            "male" => "1",
            "female" => "2",
            _ => "0",
        };
        let status = if m == proband || affected.contains(&m.as_str()) {
            "2"
        } else if unaffected.contains(&m.as_str()) {
            "1"
        } else {
            "0"
        };
        out.push_str(&["family", m, parent("father"), parent("mother"), sex, status].join("\t"));
        out.push('\n');
    }
    Ok(out)
}

fn s(v: &str) -> Value {
    Value::Str(v.into())
}

fn obj(items: Vec<(&str, Value)>) -> Value {
    Value::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

fn strs(v: &[&str]) -> Value {
    Value::Array(v.iter().map(|x| s(x)).collect())
}

/// The arguments of `prepare`.
pub struct Prepare<'a> {
    pub hpo: &'a str,
    pub obo: Option<&'a str>,
    pub sample: &'a str,
    pub vcf: &'a str,
    pub assembly: &'a str,
    pub out_dir: &'a str,
    pub analysis: &'a str,
    pub terms: &'a str,
    pub pedigree: Option<&'a str>,
    pub ped: Option<&'a str>,
    pub affected: &'a str,
    pub unaffected: &'a str,
}

/// `og exomiser prepare`: writes the files; the stdout line.
pub fn prepare(a: &Prepare<'_>) -> Result<String, Error> {
    let terms = check_terms(&split_list(a.hpo), a.obo)?;
    let empty = || IndexMap::new();
    let steps = Value::Array(
        [
            ("geneBlacklistFilter", Value::Object(empty())),
            ("failedVariantFilter", Value::Object(empty())),
            ("alleleBalanceFilter", Value::Object(empty())),
            (
                "variantEffectFilter",
                obj(vec![("remove", strs(&REMOVED_EFFECTS))]),
            ),
            (
                "frequencyFilter",
                obj(vec![("maxFrequency", Value::Float(2.0))]),
            ),
            (
                "pathogenicityFilter",
                obj(vec![("keepNonPathogenic", Value::Bool(true))]),
            ),
            ("inheritanceFilter", Value::Object(empty())),
            ("omimPrioritiser", Value::Object(empty())),
            ("hiPhivePrioritiser", Value::Object(empty())),
        ]
        .into_iter()
        .map(|(k, v)| obj(vec![(k, v)]))
        .collect(),
    );
    let mut analysis = vec![
        ("genomeAssembly", s(a.assembly)),
        ("vcf", s(a.vcf)),
        (
            "hpoIds",
            Value::Array(terms.iter().map(|t| s(&t[0])).collect()),
        ),
        (
            "inheritanceModes",
            obj(INHERITANCE
                .iter()
                .map(|(k, v)| (*k, Value::Float(*v)))
                .collect()),
        ),
        ("analysisMode", s("PASS_ONLY")),
        ("frequencySources", strs(&FREQUENCY_SOURCES)),
        ("pathogenicitySources", strs(&PATHOGENICITY_SOURCES)),
        ("steps", steps),
    ];
    let (affected, unaffected) = (split_list(a.affected), split_list(a.unaffected));
    let mut ped = None;
    if let Some(pedigree) = a.pedigree {
        let Some(out) = a.ped else {
            return Err(Error(
                "--pedigree needs --ped, the PED file to write".into(),
            ));
        };
        let text = ped_text(
            pedigree,
            &vcf_samples(a.vcf)?,
            a.sample,
            &affected,
            &unaffected,
        )?;
        ped = Some((out, text));
        analysis.push(("ped", s(out)));
        analysis.push(("proband", s(a.sample)));
    } else if !affected.is_empty() || !unaffected.is_empty() {
        return Err(Error("--affected and --unaffected name family members: they need --pedigree and the family VCF".into()));
    } else {
        let n = vcf_samples(a.vcf)?.len();
        if n != 1 {
            return Err(Error(format!(
                "{} has {n} samples; one person's VCF has one (a family VCF goes with --pedigree)",
                a.vcf
            )));
        }
    }
    let doc = obj(vec![
        ("analysis", obj(analysis)),
        (
            "outputOptions",
            obj(vec![
                ("outputContributingVariantsOnly", Value::Bool(false)),
                ("numGenes", Value::Int(OUTPUT_GENES)),
                ("outputDirectory", s(a.out_dir)),
                ("outputFileName", s(&format!("{}_exomiser", a.sample))),
                ("outputFormats", strs(&["TSV_VARIANT", "TSV_GENE", "HTML"])),
            ]),
        ),
    ]);
    let dump = Dump {
        indent: Some(2),
        ..Dump::default()
    };
    let json = pycompat::json::dumps(&doc, &dump).map_err(|e| Error(e.to_string()))? + "\n";
    let mut tsv = TERMS_HEADER.join("\t") + "\n";
    for t in &terms {
        tsv.push_str(&t.join("\t"));
        tsv.push('\n');
    }
    if let Some((path, text)) = &ped {
        crate::write_atomic(path, text.as_bytes())?;
    }
    crate::write_atomic(a.analysis, json.as_bytes())?;
    crate::write_atomic(a.terms, tsv.as_bytes())?;
    let shown: Vec<String> = terms
        .iter()
        .map(|t| {
            if t[1].is_empty() {
                t[0].clone()
            } else {
                format!("{} ({})", t[0], t[1])
            }
        })
        .collect();
    Ok(format!("HPO terms: {}\n", shown.join(", ")))
}

/// `(hpo_id, label)` of the terms TSV `prepare` wrote; the label is `None`
/// on a row that ends before it.
fn read_terms(path: &str) -> Result<Vec<(String, Option<String>)>, Error> {
    let (header, rows) = crate::read_table(path, Options::default().quoting(Quoting::None))?;
    if rows.is_empty() || TERMS_HEADER.iter().any(|k| !header.iter().any(|h| h == k)) {
        return Err(Error(format!(
            "{path} has no term rows or lacks a column of {}",
            TERMS_HEADER.join(", ")
        )));
    }
    Ok(rows
        .iter()
        .map(|r| {
            (
                r.get("hpo_id").unwrap_or("None").to_owned(),
                r.get("label").map(str::to_owned),
            )
        })
        .collect())
}

/// `"%.3f" % float(v)`, "." when it is no number.
fn number(v: &str) -> String {
    crate::float(v).map_or(".".into(), |x| crate::fixed(x, 3))
}

fn genotype_word(gt: &str) -> String {
    match gt.replace('|', "/").as_str() {
        "0/1" | "1/0" => "heterozygous".into(),
        "1/1" => "homozygous".into(),
        "1" => "hemizygous".into(),
        "0/0" | "0" => "does not carry it".into(),
        _ => gt.into(),
    }
}

/// One person's genotype in a word; a family's as each member's, named.
fn genotype_text(gt: &str, members: &[String]) -> String {
    let gts: Vec<&str> = gt.split(':').collect();
    if members.len() > 1 && gts.len() == members.len() {
        let each: Vec<String> = members
            .iter()
            .zip(&gts)
            .map(|(m, g)| format!("{m} {}", genotype_word(g)))
            .collect();
        return format!("[{}]", each.join(", "));
    }
    genotype_word(gt)
}

/// The members of the PED `prepare` wrote, in its line order; none when it
/// cannot be read.
fn ped_members(path: &Path) -> Vec<String> {
    crate::read_lines(&path.display().to_string(), Utf8::Replace)
        .map(|lines| {
            lines
                .iter()
                .filter_map(|l| {
                    let mut c = l.strip_suffix('\n').unwrap_or(l).split('\t');
                    c.nth(1).map(str::to_owned)
                })
                .collect()
        })
        .unwrap_or_default()
}

/// Exomiser's variants TSV: its rows, each checked to have one cell per
/// column.
fn read_variants(path: &str) -> Result<Vec<og_tsv::Row>, Error> {
    let opts = Options::default()
        .quoting(Quoting::None)
        .utf8(Utf8::Replace);
    let (header, rows) = crate::read_table(path, opts)?;
    if rows.is_empty() && header.first().map(String::as_str) != Some("#RANK") {
        return Err(Error(format!(
            "{path} is not an Exomiser variants TSV (no #RANK header)"
        )));
    }
    let missing: Vec<&str> = NEED
        .iter()
        .copied()
        .filter(|k| !header.iter().any(|h| h == k))
        .collect();
    if !missing.is_empty() {
        return Err(Error(format!(
            "{path} lacks the column{} {}: is it Exomiser 15's TSV_VARIANT output?",
            if missing.len() > 1 { "s" } else { "" },
            missing.join(", ")
        )));
    }
    if let Some(i) = rows
        .iter()
        .position(|r| r.dict().extra().is_some() || NEED.iter().any(|k| r.get(k).is_none()))
    {
        return Err(Error(format!(
            "{path} line {} does not have one field per column",
            i + 2
        )));
    }
    Ok(rows)
}

/// One row of the top table, by [`OUT_HEADER`].
type TopRow = [String; 16];

/// The contributing variants of the first `top` ranks, in rank order.
fn top_rows(rows: &[og_tsv::Row], top: i64) -> Result<Vec<TopRow>, Error> {
    let (mut out, mut ranks) = (Vec::new(), Vec::<i128>::new());
    for r in rows {
        let g = |k: &str| r.get(k).unwrap_or("");
        let or_dot = |k: &str| {
            if g(k).is_empty() {
                ".".to_owned()
            } else {
                g(k).to_owned()
            }
        };
        let rank = crate::int(g("#RANK"))
            .ok_or_else(|| Error(format!("rank '{}' is not a number", g("#RANK"))))?;
        if g("CONTRIBUTING_VARIANT") != "1" {
            continue;
        }
        if !ranks.contains(&rank) {
            if ranks.len() as i64 >= top {
                continue;
            }
            ranks.push(rank);
        }
        let (name, id) = (
            g("EXOMISER_ACMG_DISEASE_NAME"),
            g("EXOMISER_ACMG_DISEASE_ID"),
        );
        let disease = match (name.is_empty(), id.is_empty()) {
            (_, true) => name.to_owned(),
            (false, false) => format!("{name} ({id})"),
            (true, false) => id.to_owned(),
        };
        let cv = g("CLINVAR_PRIMARY_INTERPRETATION");
        let stars = g("CLINVAR_STAR_RATING");
        let clinvar = if cv.is_empty() || cv == "NOT_PROVIDED" {
            ".".to_owned()
        } else {
            format!(
                "{cv} ({} star{})",
                if stars.is_empty() { "0" } else { stars },
                if stars == "1" { "" } else { "s" }
            )
        };
        let freq = if g("MAX_FREQ").is_empty() {
            "not seen".to_owned()
        } else {
            format!("{} {}%", g("MAX_FREQ_SOURCE"), g("MAX_FREQ"))
        };
        out.push([
            rank.to_string(),
            g("GENE_SYMBOL").into(),
            g("MOI").into(),
            number(g("EXOMISER_GENE_COMBINED_SCORE")),
            number(g("EXOMISER_GENE_PHENO_SCORE")),
            number(g("EXOMISER_GENE_VARIANT_SCORE")),
            or_dot("P-VALUE"),
            format!("{}:{} {}>{}", g("CONTIG"), g("START"), g("REF"), g("ALT")),
            or_dot("GENOTYPE"),
            or_dot("FUNCTIONAL_CLASS"),
            or_dot("HGVS"),
            or_dot("EXOMISER_ACMG_CLASSIFICATION"),
            or_dot("EXOMISER_ACMG_EVIDENCE"),
            if disease.is_empty() {
                ".".into()
            } else {
                disease
            },
            clinvar,
            freq,
        ]);
    }
    Ok(out)
}

/// The text file's lines: what was searched, each rank in a sentence, then
/// the note. Also the number of rank lines.
fn words(
    rows: &[TopRow],
    terms: &[(String, Option<String>)],
    members: &[String],
) -> (Vec<String>, usize) {
    let named: Vec<String> = terms
        .iter()
        .map(|(id, label)| match label.as_deref() {
            Some(l) if !l.is_empty() => format!("{l} ({id})"),
            _ => id.clone(),
        })
        .collect();
    let mut lines = vec![format!("Symptoms searched: {}.", named.join(", "))];
    if rows.is_empty() {
        lines.push(
            "No gene passed Exomiser's filters with a variant that fits a mode of inheritance."
                .into(),
        );
    }
    let mut by_rank: IndexMap<&str, Vec<&TopRow>> = IndexMap::new();
    for r in rows {
        by_rank.entry(r[0].as_str()).or_default().push(r);
    }
    for (rank, vs) in &by_rank {
        let g = vs[0];
        let head = format!(
            "{rank}. {} ({}): Exomiser score {} (symptom match {}, variant {}).",
            g[1],
            moi_words(&g[2]),
            g[3],
            g[4],
            g[5]
        );
        let parts: Vec<String> = vs
            .iter()
            .map(|v| {
                let mut p = format!(
                    "{} {}, {} {}, ACMG {}",
                    v[7],
                    genotype_text(&v[8], members),
                    v[9],
                    v[10],
                    pycompat::text::lower(&v[11].replace('_', " "))
                );
                if v[12] != "." {
                    p.push_str(&format!(" ({})", v[12]));
                }
                if v[13] != "." {
                    p.push_str(&format!(", linked disease {}", v[13]));
                }
                if v[14] != "." {
                    p.push_str(&format!(
                        ", ClinVar {}",
                        pycompat::text::lower(&v[14].replace('_', " "))
                    ));
                }
                p.push_str(&format!(", highest population frequency {}", v[15]));
                p
            })
            .collect();
        lines.push(format!("{head} {}.", parts.join("; ")));
    }
    lines.push(NOTE.into());
    (lines, by_rank.len())
}

/// The arguments of `report`.
pub struct Report<'a> {
    pub sample: &'a str,
    pub variants: &'a str,
    pub terms: &'a str,
    pub out: &'a str,
    pub text: &'a str,
    pub top: i64,
}

/// `og exomiser report`: writes the table and the text; the stdout lines.
pub fn report(a: &Report<'_>) -> Result<String, Error> {
    let terms = read_terms(a.terms)?;
    let rows = top_rows(&read_variants(a.variants)?, a.top)?;
    let mut tsv = OUT_HEADER.join("\t") + "\n";
    for r in &rows {
        tsv.push_str(
            &r.iter()
                .map(|c| c.replace('\t', " "))
                .collect::<Vec<_>>()
                .join("\t"),
        );
        tsv.push('\n');
    }
    crate::write_atomic(a.out, tsv.as_bytes())?;
    let dir = std::path::absolute(a.variants)
        .ok()
        .and_then(|p| p.parent().map(Path::to_path_buf))
        .unwrap_or_default();
    let members = ped_members(&dir.join(format!("{}_exomiser.ped", a.sample)));
    let (text, ranks) = words(&rows, &terms, &members);
    crate::write_atomic(a.text, (text.join("\n") + "\n").as_bytes())?;
    // The searched line and one line per rank kept.
    let shown = if rows.is_empty() { 1 } else { 1 + ranks };
    Ok(text[..shown].join("\n") + "\n")
}

const PREPARE: Spec = Spec {
    prog: "og exomiser prepare",
    usage: "usage: og exomiser prepare [-h] --hpo HPO [--obo OBO] --sample SAMPLE --vcf VCF [--assembly {hg38,hg19}] --out-dir OUT_DIR --analysis ANALYSIS --terms TERMS [--pedigree PEDIGREE] [--ped PED] [--affected AFFECTED] [--unaffected UNAFFECTED]\n",
    options: &[
        "--hpo",
        "--obo",
        "--sample",
        "--vcf",
        "--assembly",
        "--out-dir",
        "--analysis",
        "--terms",
        "--pedigree",
        "--ped",
        "--affected",
        "--unaffected",
    ],
    positionals: &[],
    required: &[
        "--hpo",
        "--sample",
        "--vcf",
        "--out-dir",
        "--analysis",
        "--terms",
    ],
};

const REPORT: Spec = Spec {
    prog: "og exomiser report",
    usage: "usage: og exomiser report [-h] --sample SAMPLE --variants VARIANTS --terms TERMS --out OUT --text TEXT [--top TOP]\n",
    options: &[
        "--sample",
        "--variants",
        "--terms",
        "--out",
        "--text",
        "--top",
    ],
    positionals: &[],
    required: &["--sample", "--variants", "--terms", "--out", "--text"],
};

/// `og exomiser VERB ARGS`: the exit code (1 for a problem, 2 for a
/// refused command line).
pub fn cli(verb: Option<&str>, rest: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let fail = |r: Result<String, Error>| -> Outcome { r.map_err(|e| (e, 1)) };
    match verb {
        Some("prepare") => crate::run(&PREPARE, rest, out, err, |a| {
            let assembly = a.get("--assembly").unwrap_or("hg38");
            if assembly != "hg38" && assembly != "hg19" {
                return Err(og_filters::args::ArgError::Usage(format!(
                    "argument --assembly: invalid choice: {} (choose from 'hg38', 'hg19')",
                    og_filters::py_repr(assembly)
                )));
            }
            Ok(fail(prepare(&Prepare {
                hpo: crate::opt(a, "--hpo"),
                obo: a.get("--obo"),
                sample: crate::opt(a, "--sample"),
                vcf: crate::opt(a, "--vcf"),
                assembly,
                out_dir: crate::opt(a, "--out-dir"),
                analysis: crate::opt(a, "--analysis"),
                terms: crate::opt(a, "--terms"),
                pedigree: a.get("--pedigree").filter(|p| !p.is_empty()),
                ped: a.get("--ped").filter(|p| !p.is_empty()),
                affected: a.get("--affected").unwrap_or(""),
                unaffected: a.get("--unaffected").unwrap_or(""),
            })))
        }),
        Some("report") => crate::run(&REPORT, rest, out, err, |a| {
            let top = crate::int_opt(a, "--top", 10)?;
            Ok(fail(report(&Report {
                sample: crate::opt(a, "--sample"),
                variants: crate::opt(a, "--variants"),
                terms: crate::opt(a, "--terms"),
                out: crate::opt(a, "--out"),
                text: crate::opt(a, "--text"),
                top,
            })))
        }),
        other => crate::unknown_verb("exomiser", "prepare|report", other, err),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Dir;

    const OBO: &str = "format-version: 1.2\n\n[Term]\nid: HP:0001250\nname: Seizure\nalt_id: HP:0002279\n\n[Term]\nid: HP:0000001\nname: Old\nis_obsolete: true\nreplaced_by: HP:0001250 ! Seizure\n\n[Typedef]\nid: part_of\n";

    #[test]
    fn terms_against_the_ontology() {
        let d = Dir::new("exo-obo");
        let obo = d.put("hp.obo", OBO);
        let t = check_terms(&["HP:0002279", "HP:0001250"], Some(&obo)).unwrap();
        assert_eq!(
            t,
            [[
                "HP:0001250".to_owned(),
                "Seizure".into(),
                "HP:0002279".into()
            ]]
        );
        let e = check_terms(&["HP:0000001"], Some(&obo)).unwrap_err();
        assert_eq!(
            e.0,
            "HP:0000001 (Old) is obsolete in this HPO release; use HP:0001250 instead. Look the terms up at https://hpo.jax.org"
        );
        assert!(
            check_terms(&["HP:123"], None)
                .unwrap_err()
                .0
                .starts_with("not an HPO id: HP:123.")
        );
        assert_eq!(split_list(" HP:1,HP:2\tHP:3 ,, "), ["HP:1", "HP:2", "HP:3"]);
    }

    const VHEAD: &str = "#RANK\tID\tGENE_SYMBOL\tENTREZ_GENE_ID\tMOI\tP-VALUE\tEXOMISER_GENE_COMBINED_SCORE\tEXOMISER_GENE_PHENO_SCORE\tEXOMISER_GENE_VARIANT_SCORE\tEXOMISER_VARIANT_SCORE\tCONTRIBUTING_VARIANT\tWHITELIST_VARIANT\tFUNCTIONAL_CLASS\tHGVS\tEXOMISER_ACMG_CLASSIFICATION\tEXOMISER_ACMG_EVIDENCE\tEXOMISER_ACMG_DISEASE_ID\tEXOMISER_ACMG_DISEASE_NAME\tCLINVAR_VARIANT_ID\tCLINVAR_PRIMARY_INTERPRETATION\tCLINVAR_STAR_RATING\tGENE_CONSTRAINT_LOEUF\tGENE_CONSTRAINT_LOEUF_LOWER\tGENE_CONSTRAINT_LOEUF_UPPER\tMAX_FREQ_SOURCE\tMAX_FREQ\tALL_FREQ\tMAX_PATH_SOURCE\tMAX_PATH\tALL_PATH\tCONTIG\tSTART\tEND\tREF\tALT\tGENOTYPE\n";

    fn v(rank: u32, gene: &str, contributing: &str, start: u32) -> String {
        format!(
            "{rank}\tid\t{gene}\t1\tAD\t0.01\t0.9\t0.8\t0.7\t0.7\t{contributing}\t0\tmissense_variant\t{gene}:c.1A>G\tLIKELY_PATHOGENIC\tPM2\tOMIM:1\tA disease\t1\tPATHOGENIC\t2\t.\t.\t.\tGNOMAD_E_NFE\t0.01\t.\t.\t.\t.\t1\t{start}\t{start}\tA\tG\t0/1\n"
        )
    }

    /// Red first: a variant Exomiser marks 0 does not contribute.
    #[test]
    fn red_first_a_non_contributing_variant_is_left_out() {
        let d = Dir::new("exo-top");
        let p = d.put(
            "v.tsv",
            &format!(
                "{VHEAD}{}{}{}{}",
                v(1, "A", "1", 10),
                v(1, "A", "1", 20),
                v(1, "A", "0", 30),
                v(2, "B", "1", 40)
            ),
        );
        let rows = top_rows(&read_variants(&p).unwrap(), 10).unwrap();
        assert_eq!(
            rows.iter().map(|r| r[7].as_str()).collect::<Vec<_>>(),
            ["1:10 A>G", "1:20 A>G", "1:40 A>G"]
        );
        assert_eq!(rows[0][14], "PATHOGENIC (2 stars)");
        assert_eq!(rows[0][15], "GNOMAD_E_NFE 0.01%");
        let rows1 = top_rows(&read_variants(&p).unwrap(), 1).unwrap();
        assert_eq!(rows1.len(), 2);
        let (text, ranks) = words(&rows, &[("HP:0001250".into(), Some("Seizure".into()))], &[]);
        assert_eq!(ranks, 2);
        assert_eq!(text[0], "Symptoms searched: Seizure (HP:0001250).");
        assert!(text[1].starts_with("1. A (autosomal dominant): Exomiser score 0.900 (symptom match 0.800, variant 0.700). 1:10 A>G heterozygous, missense_variant A:c.1A>G, ACMG likely pathogenic (PM2), linked disease A disease (OMIM:1), ClinVar pathogenic (2 stars), highest population frequency GNOMAD_E_NFE 0.01%; 1:20"));
        assert_eq!(text.len(), 4);
    }

    #[test]
    fn family_genotypes_are_named() {
        let m = ["kid".to_owned(), "mum".into()];
        assert_eq!(
            genotype_text("0/1:0/0", &m),
            "[kid heterozygous, mum does not carry it]"
        );
        assert_eq!(genotype_text("1|1", &[]), "homozygous");
    }
}
