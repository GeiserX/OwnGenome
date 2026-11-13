//! Step 08 (T1K): the HLA alleles per gene with their qualities, the
//! candidate cap, the HFE label, the coeliac DQ haplotypes and the table
//! alleles marked present or absent; OptiType's concordance and the KIR
//! pass when they ran. The rules are docs/08-hla-typing.md's "What the
//! report derives".

use super::util::{Obj, at, bin_rows, bin_table, count, key_error, lines, s, strs, table};
use super::{Ctx, Out};
use og_steps::pgx_rules::HLA_GENES;
use pycompat::json::Value;

/// A quality-0 allele with more candidates than this is "no confident call".
pub const MAX_CANDIDATES: usize = 3;
const HFE_LABEL: &str = "T1K sequence allele, not the HFE C282Y/H63D genotype: that genotype comes from the variant calls (ClinVar section)";

/// One gene's row of T1K's genotype table.
#[derive(Debug, Clone)]
pub struct Locus {
    pub gene: String,
    pub alleles: Vec<String>,
    pub quality: Vec<String>,
    pub low_confidence: bool,
    pub candidates: Vec<usize>,
    pub one_allele_type: bool,
}

impl Locus {
    fn json(&self) -> Value {
        let mut o = Obj::new()
            .set("gene", s(&self.gene))
            .set("alleles", strs(&self.alleles))
            .set("quality", strs(&self.quality))
            .set("low_confidence", Value::Bool(self.low_confidence))
            .set(
                "withheld_from_pharmcat",
                Value::Bool(self.low_confidence && HLA_GENES.contains(&self.gene.as_str())),
            )
            .set(
                "candidates",
                Value::Array(self.candidates.iter().map(|n| count(*n)).collect()),
            )
            .set("one_allele_type", Value::Bool(self.one_allele_type));
        if name(&self.gene) == "HFE" {
            o.put("label", s(HFE_LABEL));
        }
        o.val()
    }
}

/// T1K's quality as a number; one that is not a number counts as -1, as in
/// step 36.
pub fn quality(q: &str) -> f64 {
    pycompat::num::parse_float(q).unwrap_or(-1.0)
}

/// `HLA-DQB1*02:01:01` or `DQB1*02:01:01` is `DQB1*02:01:01`.
fn name(x: &str) -> &str {
    let x = pycompat::num::py_strip(x);
    x.strip_prefix("HLA-").unwrap_or(x)
}

/// T1K's allele field: one allele, or equally likely candidates separated
/// by commas.
fn candidates(allele: &str) -> Vec<&str> {
    allele.split(',').map(name).filter(|a| !a.is_empty()).collect()
}

/// The alleles and their qualities of columns 3 and 6 (`.` is no allele).
fn alleles(c: &[&str], skip_empty: bool) -> (Vec<String>, Vec<String>) {
    let (mut a, mut q) = (Vec::new(), Vec::new());
    for i in [2, 5] {
        let x = c.get(i).copied().unwrap_or("");
        if (skip_empty && x.is_empty()) || x == "." || x.is_empty() {
            continue;
        }
        a.push(x.to_owned());
        q.push(c.get(i + 2).copied().unwrap_or("").to_owned());
    }
    (a, q)
}

/// One locus with the candidate cap, the candidate counts and whether T1K
/// typed one allele.
fn locus(gene: &str, alleles: Vec<String>, quals: Vec<String>, n_diff: &str) -> Locus {
    // Step 36's rule: an allele quality of 0 or below withholds the gene.
    let low = quals.iter().any(|q| quality(q) <= 0.0);
    let mut counts = Vec::new();
    let mut kept = Vec::new();
    for (i, a) in alleles.iter().enumerate() {
        let q = quals.get(i).map_or("", String::as_str);
        let cands = candidates(a);
        counts.push(cands.len());
        kept.push(if quality(q) <= 0.0 && cands.len() > MAX_CANDIDATES {
            cands[..MAX_CANDIDATES].join(",")
        } else {
            a.clone()
        });
    }
    let one = pycompat::num::py_strip(n_diff) == "1" && kept.len() == 1;
    Locus {
        gene: gene.to_owned(),
        alleles: kept,
        quality: quals,
        low_confidence: low,
        candidates: counts,
        one_allele_type: one,
    }
}

/// `yes` when every candidate of the allele is in the group, `no` when
/// none is, `mixed` when only some are.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Match {
    Yes,
    No,
    Mixed,
}

fn matches(allele: &str, group: &str) -> Match {
    let prefix = format!("{group}:");
    let hits: Vec<bool> = candidates(allele)
        .iter()
        .map(|c| *c == group || c.starts_with(&prefix))
        .collect();
    if !hits.is_empty() && hits.iter().all(|h| *h) {
        Match::Yes
    } else if hits.iter().any(|h| *h) {
        Match::Mixed
    } else {
        Match::No
    }
}

/// Over both alleles of a locus.
fn carries(l: &Locus, group: &str) -> Match {
    let m: Vec<Match> = l.alleles.iter().map(|a| matches(a, group)).collect();
    if m.contains(&Match::Yes) {
        Match::Yes
    } else if m.contains(&Match::Mixed) {
        Match::Mixed
    } else {
        Match::No
    }
}

fn both(a: Match, b: Match) -> &'static str {
    if a == Match::Yes && b == Match::Yes {
        "present"
    } else if a == Match::No || b == Match::No {
        "absent"
    } else {
        "uncertain"
    }
}

/// The DQ haplotypes, the DQB1*02 dose and each table allele present or
/// absent.
fn derived(loci: &[Locus], table: &[indexmap::IndexMap<String, String>]) -> Result<Value, String> {
    let mut by: indexmap::IndexMap<&str, &Locus> = indexmap::IndexMap::new();
    for l in loci {
        by.entry(name(&l.gene)).or_insert(l);
    }
    let usable = |gene: &str| -> (Option<&Locus>, String) {
        match by.get(gene) {
            Some(l) if !l.alleles.is_empty() => {
                if l.low_confidence {
                    (None, format!("HLA-{gene} low confidence"))
                } else {
                    (Some(*l), String::new())
                }
            }
            _ => (None, format!("HLA-{gene} not typed")),
        }
    };
    let (dqa, ra) = usable("DQA1");
    let (dqb, rb) = usable("DQB1");
    let reason = [ra, rb]
        .into_iter()
        .filter(|r| !r.is_empty())
        .collect::<Vec<_>>()
        .join("; ");
    let mut co = Obj::new()
        .set("status", s("not_determined"))
        .set("reason", s(&reason))
        .set("dq2_5", Value::Null)
        .set("dq2_5_phase", Value::Null)
        .set("dq2_2", Value::Null)
        .set("dq8", Value::Null)
        .set("dqb1_02_dose", Value::Null);
    if let Some(dqb) = dqb {
        let m: Vec<Match> = dqb.alleles.iter().map(|a| matches(a, "DQB1*02")).collect();
        if !m.contains(&Match::Mixed) {
            let n = m.iter().filter(|x| **x == Match::Yes).count();
            co.put(
                "dqb1_02_dose",
                count(if dqb.one_allele_type { 2 * n } else { n }),
            );
        }
    }
    if let (Some(dqa), Some(dqb)) = (dqa, dqb) {
        co.put("status", s("derived"));
        co.put("reason", s(""));
        let dq25 = both(carries(dqa, "DQA1*05"), carries(dqb, "DQB1*02"));
        co.put("dq2_5", s(dq25));
        co.put(
            "dq2_2",
            s(both(carries(dqa, "DQA1*02:01"), carries(dqb, "DQB1*02:02"))),
        );
        co.put(
            "dq8",
            s(both(carries(dqa, "DQA1*03"), carries(dqb, "DQB1*03:02"))),
        );
        if dq25 == "present" {
            let (drb1, why) = usable("DRB1");
            let phase = match drb1.map(|l| carries(l, "DRB1*03")) {
                Some(Match::Yes) => "likely cis (DRB1*03 present)".to_owned(),
                Some(Match::No) => "possible trans (no DRB1*03)".to_owned(),
                _ => format!(
                    "phase unknown ({})",
                    if why.is_empty() {
                        "DRB1*03 uncertain"
                    } else {
                        &why
                    }
                ),
            };
            co.put("dq2_5_phase", s(&phase));
        }
    }
    let mut alleles = Vec::new();
    for r in table {
        let cell = |k: &str| r.get(k).cloned().ok_or_else(|| key_error(k));
        let allele = cell("allele")?;
        let gene = allele.split('*').next().unwrap_or_default();
        let (l, why) = usable(gene);
        let state = match l.map(|l| carries(l, &allele)) {
            Some(Match::Yes) => "present",
            Some(Match::No) => "absent",
            Some(Match::Mixed) => "uncertain",
            None => "not_determined",
        };
        alleles.push(
            Obj::new()
                .set("allele", s(&allele))
                .set("kind", s(&cell("kind")?))
                .set("association", s(&cell("association")?))
                .set("pmids", s(&cell("pmids")?))
                .set("citation", s(&cell("citation")?))
                .set("state", s(state))
                .set("reason", s(&why))
                .val(),
        );
    }
    Ok(Obj::new()
        .set("coeliac", co.val())
        .set("alleles", Value::Array(alleles))
        .val())
}

/// The first `#` line of a file step 08 wrote instead of a result.
fn first_comment(path: &str) -> Result<String, String> {
    Ok(lines(path)?
        .iter()
        .find(|l| l.starts_with('#'))
        .map(|l| pycompat::num::py_strip(l.trim_start_matches('#')).to_owned())
        .unwrap_or_default())
}

/// OptiType's concordance with T1K per class I gene, or why it did not run.
fn optitype(c: &Ctx<'_>) -> Result<Option<Value>, String> {
    let conc = c.find1(format!("hla_optitype/{}_hla_concordance.tsv", c.s()));
    let res = c.find1(format!("hla_optitype/{}_optitype_result.tsv", c.s()));
    let Some(conc) = conc else {
        let Some(res) = res else { return Ok(None) };
        let why = first_comment(&res.path)?;
        let why = if why.is_empty() {
            "no concordance table: step 08 did not compare the callers".to_owned()
        } else {
            why
        };
        return Ok(Some(
            Obj::new()
                .set("genes", Value::Array(Vec::new()))
                .set("reason", s(&why))
                .val(),
        ));
    };
    let mut genes = Vec::new();
    for r in table(&conc.path, false)?.rows {
        let cell = |k: &str| pycompat::num::py_strip(r.get(k).flatten().unwrap_or("")).to_owned();
        if !cell("Gene").is_empty() {
            genes.push(
                Obj::new()
                    .set("gene", s(&cell("Gene")))
                    .set("t1k", s(&cell("T1K")))
                    .set("optitype", s(&cell("OptiType")))
                    .set("agree", s(&cell("Agree")))
                    .set("note", s(&cell("Note")))
                    .val(),
            );
        }
    }
    let reason = if genes.is_empty() {
        "the concordance table has no rows"
    } else {
        ""
    };
    Ok(Some(
        Obj::new()
            .set("genes", Value::Array(genes))
            .set("reason", s(reason))
            .val(),
    ))
}

/// The KIR pass: the genes T1K gave an allele, with each quality, and the
/// genes it gave none.
fn kir(c: &Ctx<'_>) -> Result<Option<Value>, String> {
    let s_ = c.s();
    let Some(p) = c.find(&[
        format!("kir_t1k/{s_}_kir_genotype.tsv"),
        format!("kir/{s_}_kir_genotype.tsv"),
    ]) else {
        return Ok(None);
    };
    let (mut typed, mut not_called) = (Vec::new(), Vec::new());
    for l in lines(&p.path)? {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 3 || c[0].starts_with('#') {
            continue;
        }
        let (a, q) = alleles(&c, true);
        if a.is_empty() {
            not_called.push(c[0].to_owned());
        } else {
            typed.push(
                Obj::new()
                    .set("gene", s(c[0]))
                    .set("alleles", strs(&a))
                    .set("quality", strs(&q))
                    .val(),
            );
        }
    }
    let reason = if !typed.is_empty() || !not_called.is_empty() {
        String::new()
    } else {
        let why = first_comment(&p.path)?;
        if why.is_empty() {
            "the KIR genotype file has no gene rows".to_owned()
        } else {
            why
        }
    };
    Ok(Some(
        Obj::new()
            .set("typed", Value::Array(typed))
            .set("not_called", strs(&not_called))
            .set("reason", s(&reason))
            .val(),
    ))
}

/// Step 08's T1K genotype table, the derived calls, and OptiType and KIR
/// when they ran.
pub fn hla(c: &Ctx<'_>) -> Out {
    let s_ = c.s();
    let Some(p) = c.find(&[
        format!("hla_t1k/{s_}_hla_genotype.tsv"),
        format!("hla/{s_}_hla_genotype.tsv"),
        "hla/*_hla_genotype.tsv".to_owned(),
    ]) else {
        return Ok(None);
    };
    let mut loci = Vec::new();
    for l in lines(&p.path)? {
        let c: Vec<&str> = l.split('\t').collect();
        if c.len() < 3 || c[0].starts_with('#') {
            continue;
        }
        let (a, q) = alleles(&c, false);
        loci.push(locus(c[0], a, q, c.get(1).copied().unwrap_or("")));
    }
    let table = bin_rows(&bin_table("hla_alleles.tsv"))?;
    let mut out = Obj::new()
        .set("loci", Value::Array(loci.iter().map(Locus::json).collect()))
        .set("derived", derived(&loci, &table)?);
    if let Some(o) = optitype(c)? {
        out.put("optitype", o);
    }
    if let Some(k) = kir(c)? {
        out.put("kir", k);
    }
    let _ = at;
    Ok(Some((p, out.val())))
}
