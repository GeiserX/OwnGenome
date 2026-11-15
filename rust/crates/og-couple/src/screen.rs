//! `og couple screen`: step 39's couple carrier screen (`scripts/39-couple.sh`
//! and the Nextflow `COUPLE_CARRIER` process run the same command).
//!
//! ```text
//! og couple screen --pedigree P.tsv --genes assets/carrier_genes.tsv
//!                  --clinvar ID=HITS.vcf ... [--delly ID=SV.vcf.gz ...]
//!                  [--eh ID=EH.vcf ...] [--smn ID=SMN.tsv ...]
//!                  [--exons EXONS.bed] [--benign-sv BENIGN_LOSS.bed]
//!                  [--cn-intervals ID=CN_INTERVALS.tsv ...] [--cnvpytor ID=CNVS.vcf.gz ...]
//!                  [--cyp21a2 ID=CYP21A2_PSV.tsv ...] [--par-bed PAR.bed]
//!                  [--declared-sex ID=SEX ...] --out-dir D
//! ```
//!
//! For every two samples the pedigree marks as partners, and that both have
//! a ClinVar screen (step 06's `<id>_clinvar_hits.vcf`), it writes
//! `D/couple_A_B.tsv` (one row per gene where either partner has a finding:
//! shared, x_linked or one_partner), `D/couple_A_B.md` (the same in words,
//! what was not assessed, and the residual risk) and a copy of the table for
//! each partner's own report, `D/A_couple_carrier.tsv` and
//! `D/B_couple_carrier.tsv`.
//!
//! Without `--par-bed` the PARs are GRCh38's, the rows of `assets/par_grch38.bed`.
//!
//! Every couple is computed before anything is written: an input problem in
//! any couple stops the screen with no table at all, never with the tables
//! of the couples before it. Exits 1 when no couple could be screened, 2 on
//! an input problem.

use crate::args::{self, Kind, Spec};
use crate::cn::{self, Intervals};
use crate::findings::{self, Benign, DEPTH_MIN, Exons, Finding, Pars, SV_MIN_SIZE};
use crate::list::{CarrierList, GeneIndex};
use crate::prose::{self, Assessed, Outcome, Report};
use og_pedigree::Die;
use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::path::Path;

const USAGE: &str =
    "usage: og couple screen --pedigree PEDIGREE --genes GENES [--clinvar ID=VCF] [--delly ID=VCF]
                        [--eh ID=VCF] [--smn ID=TSV] [--cn-intervals ID=TSV] [--cyp21a2 ID=TSV]
                        [--cnvpytor ID=VCF] [--exons BED] [--benign-sv BED] [--par-bed BED]
                        [--declared-sex ID=SEX] --out-dir OUT_DIR
";

const SPEC: Spec = Spec {
    prog: "og couple screen",
    usage: USAGE,
    options: &[
        ("--pedigree", Kind::One),
        ("--genes", Kind::One),
        ("--clinvar", Kind::Append),
        ("--delly", Kind::Append),
        ("--eh", Kind::Append),
        ("--smn", Kind::Append),
        ("--cn-intervals", Kind::Append),
        ("--cyp21a2", Kind::Append),
        ("--cnvpytor", Kind::Append),
        ("--exons", Kind::One),
        ("--benign-sv", Kind::One),
        ("--par-bed", Kind::One),
        ("--declared-sex", Kind::Append),
        ("--out-dir", Kind::One),
    ],
    required: &["--pedigree", "--genes", "--out-dir"],
};

/// The per-sample inputs, in the order they are checked.
const INPUTS: [&str; 7] = [
    "clinvar",
    "delly",
    "eh",
    "smn",
    "cn-intervals",
    "cnvpytor",
    "cyp21a2",
];

/// `og couple screen ARGS`.
pub fn cli(argv: &[String], out: &mut dyn Write, err: &mut dyn Write) -> u8 {
    let a = match args::parse(&SPEC, argv, out, err) {
        Ok(a) => a,
        Err(code) => return code,
    };
    let mut said = Vec::new();
    let result = run(&a, &mut said);
    // An empty line holds the place of a couple that was computed but, the
    // screen having stopped, never written: it says nothing.
    for line in said.iter().filter(|l| !l.is_empty()) {
        let _ = writeln!(out, "{line}");
    }
    match result {
        Ok(()) => 0,
        Err(die) => {
            let _ = out.flush();
            for line in &die.lines {
                let _ = writeln!(err, "ERROR: {line}");
            }
            die.code
        }
    }
}

/// `--KIND ID=PATH ...` as `[(ID, PATH)]` in order; a malformed or repeated
/// ID stops the screen.
fn kv(values: &[String], what: &str) -> Result<Vec<(String, String)>, Die> {
    let mut out: Vec<(String, String)> = Vec::new();
    for v in values {
        let (k, p) = v.split_once('=').unwrap_or((v, ""));
        if k.is_empty() || p.is_empty() || !v.contains('=') {
            return Err(Die::one(format!("--{what} takes ID=PATH, got '{v}'")));
        }
        if out.iter().any(|(o, _)| o == k) {
            return Err(Die::one(format!("--{what} names '{k}' twice")));
        }
        out.push((k.to_owned(), p.to_owned()));
    }
    Ok(out)
}

/// The files of one couple, ready to write.
struct Written {
    base: String,
    table: String,
    partners: [String; 2],
    markdown: String,
}

fn run(a: &args::Args, said: &mut Vec<String>) -> Result<(), Die> {
    let req = |k: &str| a.one(k).unwrap_or_default();
    let genes_path = req("--genes");
    let list = CarrierList::read(genes_path)?;
    let index = GeneIndex::new(&list);
    let pedigree = req("--pedigree");
    let (rows, errors) = og_pedigree::read_pedigree(pedigree)?;
    let mut declared = HashMap::new();
    for v in a.many("--declared-sex") {
        let (k, s) = v.split_once('=').unwrap_or((v, ""));
        let word = pycompat::text::lower(pycompat::text::strip(s, None));
        if let Some(Some(sex)) = og_pedigree::sex_word(&word) {
            declared.insert(k.to_owned(), sex.to_owned());
        }
    }
    let (members, more) = og_pedigree::build_members(&rows, &declared, pedigree, "the samplesheet");
    if !errors.is_empty() || !more.is_empty() {
        return Err(Die::new(errors.into_iter().chain(more).collect()));
    }
    let mut inputs: HashMap<&str, Vec<(String, String)>> = HashMap::new();
    for what in INPUTS {
        inputs.insert(what, kv(a.many(&format!("--{what}")), what)?);
    }
    for what in INPUTS {
        for (k, p) in &inputs[what] {
            if !members.contains_key(k) {
                return Err(Die::one(format!(
                    "--{what} {k}: '{k}' is not in the pedigree"
                )));
            }
            if !Path::new(p).is_file() {
                return Err(Die::one(format!("--{what} {k}: no file {p}")));
            }
        }
    }
    let of = |what: &str, s: &str| -> Option<&str> {
        inputs[what]
            .iter()
            .find(|(k, _)| k == s)
            .map(|(_, p)| p.as_str())
    };
    let exons = a.one("--exons").map(Exons::read).transpose()?;
    let mut benign = None;
    let mut common = "not filtered, no benign SV set was given (AnnotSV's annotations, step 05): a deletion many \
                      people carry can show as a finding"
        .to_owned();
    if let Some(p) = a.one("--benign-sv") {
        if !Path::new(p).is_file() {
            return Err(Die::one(format!("--benign-sv: no file {p}")));
        }
        let b = Benign::read(p)?;
        if b.count > 0 {
            common = format!(
                "a deletion matching one of 1% or more in AnnotSV's benign set ({} record{}, {}% reciprocal \
                 overlap) is set aside",
                b.count,
                if b.count == 1 { "" } else { "s" },
                findings::COMMON_OVERLAP_PERCENT
            );
            benign = Some(b);
        } else {
            common = "not filtered, the benign SV set holds no deletion of 1% or more".to_owned();
        }
    }
    let par = Pars::read(a.one("--par-bed"))?;
    let sexes: HashMap<String, Option<String>> = members
        .iter()
        .map(|(s, m)| (s.clone(), m.sex.clone()))
        .collect();
    let couples: BTreeSet<(String, String)> = members
        .iter()
        .filter_map(|(s, m)| {
            let p = m.partner.as_ref()?;
            Some(if s <= p {
                (s.clone(), p.clone())
            } else {
                (p.clone(), s.clone())
            })
        })
        .collect();
    if couples.is_empty() {
        return Err(Die {
            lines: vec![format!(
                "{pedigree} marks no two people as partners (the partner column): there is no couple to screen."
            )],
            code: 1,
        });
    }
    let out_dir = req("--out-dir");
    std::fs::create_dir_all(out_dir)
        .map_err(|e| Die::one(format!("cannot create {out_dir}: {e}")))?;
    let mut done: Vec<Written> = Vec::new();
    // What stdout says, line by line, for each couple once it is written.
    let mut lines_of: Vec<Vec<String>> = Vec::new();
    for (x, y) in &couples {
        let lack: Vec<&str> = [x, y]
            .into_iter()
            .filter(|s| of("clinvar", s).is_none())
            .map(String::as_str)
            .collect();
        if !lack.is_empty() {
            said.push(format!(
                "{x} and {y}: not screened, {} {} no ClinVar screen (step 06's <sample>_clinvar_hits.vcf)",
                lack.join(" and "),
                if lack.len() == 1 { "has" } else { "have" }
            ));
            continue;
        }
        let mut found: HashMap<String, Vec<Finding>> = HashMap::new();
        let mut assessed: HashMap<String, Assessed> = HashMap::new();
        let mut gaps: Vec<String> = Vec::new();
        let mut cn_lines: Vec<String> = Vec::new();
        let mut ivs: HashMap<String, Intervals> = HashMap::new();
        for s in [x, y] {
            let sex = sexes.get(s).and_then(|v| v.as_deref());
            if let Some(p) = of("cn-intervals", s) {
                let iv = cn::read_intervals(p)?;
                let dels = of("cnvpytor", s).map(cn::cnvpytor_deletions).transpose()?;
                cn_lines.extend(cn::lines(s, sex, &iv, dels.as_ref())?);
                ivs.insert(s.clone(), iv);
            } else if !inputs["cn-intervals"].is_empty() {
                gaps.push(format!(
                    "Carrier deletion alleles of {s}: no interval table (step 35)."
                ));
            }
            let cv = findings::clinvar(of("clinvar", s).unwrap_or_default(), &list, sex, &par)?;
            let mut f = cv.findings;
            let low = f.iter().filter(|x| x.hypomorph.is_none() && x.low).count();
            let mut info = Assessed {
                clinvar: format!(
                    "{} / {}{}",
                    f.iter().filter(|x| x.hypomorph.is_none()).count(),
                    cv.total,
                    if low > 0 {
                        format!(" ({low} low-confidence)")
                    } else {
                        String::new()
                    }
                ),
                ..Assessed::default()
            };
            for note in &cv.aside {
                gaps.push(format!(
                    "Set aside for {s}, not a carrier finding on its own: {note}."
                ));
            }
            match (of("delly", s), &exons) {
                (None, _) => {
                    info.delly = "not assessed (no Delly calls, step 19)".to_owned();
                    gaps.push(format!("Deletions of {s}: no Delly calls (step 19)."));
                }
                (Some(_), None) => {
                    info.delly = "not assessed (no exon list)".to_owned();
                    gaps.push(format!(
                        "Deletions of {s}: Delly ran, but without the GENCODE exon list a deletion cannot be \
                         placed on exons."
                    ));
                }
                (Some(p), Some(exons)) => {
                    let (d, skip) = findings::delly(p, &index, exons, sex, &par, benign.as_ref())?;
                    info.delly = format!("{} over exons of listed genes", d.len());
                    if skip.short > 0 {
                        info.delly +=
                            &format!("; {} shorter than {SV_MIN_SIZE} bp set aside", skip.short);
                    }
                    if skip.no_loss > 0 {
                        info.delly += &format!(
                            "; {} longer than {DEPTH_MIN} bp set aside, read depth does not confirm the loss",
                            skip.no_loss
                        );
                    }
                    if skip.common > 0 {
                        info.delly += &format!("; {} set aside as common deletions", skip.common);
                    }
                    if benign.is_none() {
                        gaps.push(format!("Common deletions of {s}: {common}."));
                    }
                    f.extend(d);
                }
            }
            if let Some(p) = of("smn", s) {
                let (m, txt) = findings::smn(p, &list)?;
                if m.is_empty()
                    && (txt.contains("not resolved")
                        || txt.contains("not assessed")
                        || txt.contains("no row"))
                {
                    gaps.push(format!("SMN1 copy number of {s}: {txt}."));
                }
                info.smn = txt;
                f.extend(m);
            } else {
                info.smn = "not assessed (no Parascopy, step 35)".to_owned();
                gaps.push(format!("SMN1 copy number of {s}: no Parascopy (step 35)."));
            }
            if let Some(p) = of("cyp21a2", s) {
                let c = findings::cyp21a2(p, &list)?;
                f.extend(c.findings);
                if let Some(gap) = c.gap {
                    gaps.push(format!("CYP21A2 deletion and conversions of {s}: {gap}."));
                }
            } else {
                gaps.push(format!(
                    "CYP21A2 deletion and conversions of {s}: step 35's CYP21A2 table was not passed to the \
                     couple screen (--cyp21a2)."
                ));
            }
            if let Some(p) = of("eh", s) {
                let r = findings::repeats(p, &list)?;
                info.eh = format!("{} expanded", r.len());
                f.extend(r);
            } else {
                info.eh = "not assessed (no ExpansionHunter, step 09)".to_owned();
                gaps.push(format!(
                    "FMR1 and FXN repeats of {s}: no ExpansionHunter (step 09)."
                ));
            }
            found.insert(s.clone(), f);
            assessed.insert(s.clone(), info);
        }
        let rows = prose::screen_couple(x, y, &sexes, &found, &list);
        let shared_cn = cn::shared(x, y, &ivs, &found, &list);
        let base = Path::new(out_dir)
            .join(format!("couple_{x}_{y}"))
            .to_string_lossy()
            .into_owned();
        let markdown = Report {
            a: x,
            b: y,
            sexes: &sexes,
            rows: &rows,
            assessed: &assessed,
            gaps: &gaps,
            genes: list.len(),
            common: &common,
            cn: &cn_lines,
            shared_cn: &shared_cn,
        }
        .markdown();
        let count = |k: Outcome| rows.iter().filter(|r| r.outcome == k).count();
        let shared = count(Outcome::Shared);
        let mut lines = vec![format!(
            "{x} and {y}: {shared} shared gene{}, {} X-linked, {} with one partner only ({base}.tsv)",
            if shared == 1 { "" } else { "s" },
            count(Outcome::XLinked),
            count(Outcome::OnePartner)
        )];
        lines.extend(
            rows.iter()
                .filter(|r| r.outcome == Outcome::Shared)
                .map(|r| format!("  shared: {} ({}): {}", r.gene, r.condition, r.risk)),
        );
        lines.extend(shared_cn.iter().map(|l| {
            let head = l.split(':').next().unwrap_or("");
            format!(
                "  possibly shared, from read depth: {}",
                head.trim_matches('*')
            )
        }));
        // A placeholder keeps "not screened" lines in their place.
        said.push(String::new());
        lines_of.push(lines);
        done.push(Written {
            base,
            table: prose::table(&rows),
            partners: [x.clone(), y.clone()],
            markdown,
        });
    }
    if done.is_empty() {
        return Err(Die {
            lines: vec![
                "no couple of the pedigree has a ClinVar screen for both partners: run step 06 for each".to_owned(),
            ],
            code: 1,
        });
    }
    // Every couple is computed: now write them all, then say so.
    for w in &done {
        write(&format!("{}.tsv", w.base), &w.table)?;
        for s in &w.partners {
            let p = Path::new(out_dir).join(format!("{s}_couple_carrier.tsv"));
            write(&p.to_string_lossy(), &w.table)?;
        }
        write(&format!("{}.md", w.base), &w.markdown)?;
    }
    let mut each = lines_of.into_iter();
    let mut full = Vec::new();
    for line in said.drain(..) {
        if line.is_empty() {
            full.extend(each.next().unwrap_or_default());
        } else {
            full.push(line);
        }
    }
    *said = full;
    Ok(())
}

/// Writes `text` to `path` through a temporary file renamed into place.
fn write(path: &str, text: &str) -> Result<(), Die> {
    let mut f = og_io::AtomicFile::create(path)
        .map_err(|e| Die::one(format!("cannot write {path}: {e}")))?;
    f.write_all(text.as_bytes())
        .and_then(|()| f.commit())
        .map_err(|e| Die::one(format!("cannot write {path}: {e}")))
}
