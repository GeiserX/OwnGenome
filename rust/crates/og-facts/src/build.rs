//! The facts document of one finished sample folder, from its stable
//! contract files (docs/output-contract.md): `summary.json`, the ClinVar hits
//! table, PharmCAT's report, and the calls VCF and gVCF it lists by hash.
//!
//! What goes in, and how each field is decided, is in docs/ownhealth.md. In
//! short: the ClinVar P/LP hits (as `acmg_sf`, `carrier` or
//! `monogenic_pathogenic`), the PGx calls with their CPIC and DPWG guidance,
//! the drug-relevant HLA alleles, the polygenic scores, APOE and Lp(a). A
//! section the run did not produce adds no fact; the coverage map says
//! `not_assessed` for what it would have covered.

use std::collections::{BTreeMap, HashSet};
use std::io::Read;
use std::path::Path;

use pycompat::json::{self, Value};

use crate::panels::{Panels, Site};
use crate::rules::{clean_text, drug_ok, fact_id, gene_ok, hex, term_ok};
use crate::vrs::{self, Bases};

/// The genome the calls were made on, read for VRS normalisation.
pub trait Genome {
    /// The bases `[beg, end)` of `chrom`, 0-based.
    fn get(&mut self, chrom: &str, beg: u64, end: u64) -> Result<Vec<u8>, String>;
}

/// One chromosome of a [`Genome`] as VRS reads it.
struct On<'a> {
    genome: &'a mut dyn Genome,
    chrom: &'a str,
    len: u64,
}

impl Bases for On<'_> {
    fn chrom_len(&self) -> u64 {
        self.len
    }
    fn get(&mut self, beg: u64, end: u64) -> Result<Vec<u8>, String> {
        let mut b = self.genome.get(self.chrom, beg, end)?;
        b.make_ascii_uppercase();
        Ok(b)
    }
}

/// The document and what was left out of it, one line each.
pub struct Built {
    pub doc: Value,
    pub warnings: Vec<String>,
}

/// The ACMG SF version the rules of `assets/acmg_sf.tsv` follow.
pub const ACMG_SF: &str = "acmg_sf_v3.3";

fn s(text: &str) -> Value {
    Value::Str(text.to_owned())
}

fn obj<const N: usize>(pairs: [(&str, Value); N]) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

/// Sets `key` on an object.
fn put(o: &mut Value, key: &str, v: Value) {
    if let Value::Object(m) = o {
        m.insert(key.to_owned(), v);
    }
}

fn get<'a>(v: &'a Value, path: &[&str]) -> Option<&'a Value> {
    path.iter().try_fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k),
        _ => None,
    })
}

fn text<'a>(v: &'a Value, path: &[&str]) -> Option<&'a str> {
    match get(v, path) {
        Some(Value::Str(t)) => Some(t),
        _ => None,
    }
}

/// A number, also from a string ("76.9" in the PRS table).
fn number(v: Option<&Value>) -> Option<f64> {
    match v {
        Some(Value::Int(n)) => Some(*n as f64),
        Some(Value::Float(x)) if x.is_finite() => Some(*x),
        Some(Value::Str(t)) => t.trim().parse().ok().filter(|x: &f64| x.is_finite()),
        _ => None,
    }
}

fn int(v: Option<&Value>) -> Option<i128> {
    match v {
        Some(Value::Int(n)) => Some(*n),
        Some(Value::Str(t)) => t.trim().parse().ok(),
        _ => None,
    }
}

fn list<'a>(v: &'a Value, path: &[&str]) -> &'a [Value] {
    match get(v, path) {
        Some(Value::Array(a)) => a,
        _ => &[],
    }
}

/// The data of a section of summary.json when its state is `ok`.
fn section<'a>(summary: &'a Value, name: &str) -> Option<&'a Value> {
    let sec = get(summary, &["sections", name])?;
    (text(sec, &["state"]) == Some("ok")).then(|| get(sec, &["data"]))?
}

/// `2026-11-16 04:19 UTC` (summary.json) as `2026-11-16T04:19:00Z`.
pub fn iso_utc(t: &str) -> Option<String> {
    let t = t.trim().strip_suffix("UTC")?.trim();
    let (date, time) = t.split_once(' ')?;
    let ok = |p: &str, n: usize| p.len() == n && p.bytes().all(|b| b.is_ascii_digit());
    let d: Vec<&str> = date.split('-').collect();
    let h: Vec<&str> = time.split(':').collect();
    let date_ok = d.len() == 3 && ok(d[0], 4) && ok(d[1], 2) && ok(d[2], 2);
    let time_ok = (2..=3).contains(&h.len()) && h.iter().all(|p| ok(p, 2));
    (date_ok && time_ok).then(|| {
        let sec = h.get(2).copied().unwrap_or("00");
        format!("{date}T{}:{}:{sec}Z", h[0], h[1])
    })
}

/// sha256 and size of a file, read in 1 MiB pieces.
fn file_sha256(path: &Path) -> Result<(String, u64), String> {
    use sha2::{Digest, Sha256};
    let mut f = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut n = 0u64;
    loop {
        let k = f
            .read(&mut buf)
            .map_err(|e| format!("{}: {e}", path.display()))?;
        if k == 0 {
            break;
        }
        h.update(&buf[..k]);
        n += k as u64;
    }
    Ok((hex(&h.finalize()), n))
}

/// What the build needs to know as it goes.
struct Ctx<'a> {
    panels: &'a Panels,
    genome: &'a mut dyn Genome,
    facts: Vec<Value>,
    ids: HashSet<String>,
    warnings: Vec<String>,
    /// vcf_only: the run started from a VCF, so no read-level checks.
    vcf_only: bool,
    /// What a finding the build left out does to the coverage map: the
    /// gene (or `chrom:pos` locus) and a note naming the call.
    left_out: Vec<(String, String)>,
}

impl Ctx<'_> {
    fn push(&mut self, fact: Value) {
        let id = text(&fact, &["id"]).unwrap_or_default().to_owned();
        if self.ids.insert(id) {
            self.facts.push(fact);
        }
    }

    fn warn(&mut self, w: String) {
        self.warnings.push(w);
    }

    /// A finding the file cannot hold: a warning, and its gene (or locus)
    /// is `partial` in the coverage map, so "no finding" never stands for
    /// a call that was made and left out.
    fn leave_out(&mut self, id: &str, call: &str, why: &str) {
        self.warn(format!("{call}: left out, {why}"));
        self.left_out.push((
            id.to_owned(),
            format!("{call} was left out of this file: {why}"),
        ));
    }

    /// The VRS identifier of `chrom:pos ref>alt` on GRCh38.
    fn vrs_id(&mut self, chrom: &str, pos: u64, ref_: &str, alt: &str) -> Result<String, String> {
        let Some((refget, len)) = self.panels.refget.get(chrom).cloned() else {
            return Err(format!("{chrom} is not a GRCh38 chromosome"));
        };
        let mut on = On {
            genome: &mut *self.genome,
            chrom,
            len,
        };
        let allele = vrs::normalize(&mut on, pos, ref_, alt)?;
        Ok(vrs::allele_id(&refget, &allele))
    }
}

/// The schema's `inheritance` for a table's value, when it has one.
fn inheritance(t: &str) -> Option<&'static str> {
    match t {
        "AD" => Some("AD"),
        "AR" => Some("AR"),
        "XL" | "XLR" | "XLD" => Some("XL"),
        "MT" => Some("MT"),
        _ => None,
    }
}

fn gt_of(zygosity: &str) -> Option<&'static str> {
    match zygosity {
        "het" => Some("0/1"),
        "hom" => Some("1/1"),
        "hemi" => Some("1"),
        _ => None,
    }
}

/// The ClinVar P/LP hits of step 06.
fn clinvar_facts(ctx: &mut Ctx<'_>, sample: &str, dir: &Path, summary: &Value, female: bool) {
    let Some(data) = section(summary, "clinvar") else {
        return;
    };
    // ClinVar VariationIDs, from the hits table (its clinvar_id column).
    let mut variation: BTreeMap<String, i128> = BTreeMap::new();
    let table = dir.join(format!("clinvar/{sample}_clinvar_hits.tsv"));
    if let Ok(t) = std::fs::read_to_string(&table) {
        for line in t.lines().skip(1) {
            let c: Vec<&str> = line.split('\t').collect();
            if let (Some(k), Some(Ok(id))) = (c.get(..4), c.get(5).map(|x| x.parse::<i128>())) {
                variation.insert(k.join(":"), id);
            }
        }
    }
    // The hits step 23 reports under ACMG SF v3.3 (its report_if rules met).
    let sf: HashSet<&str> = list(
        summary,
        &["sections", "clinical", "data", "acmg_sf", "clinvar_hits"],
    )
    .iter()
    .filter_map(|h| text(h, &["variant"]))
    .collect();
    let floor_gq = int(get(data, &["call_floor", "GQ"])).unwrap_or(30);
    let floor_dp = int(get(data, &["call_floor", "DP"])).unwrap_or(10);
    for h in list(data, &["hits"]) {
        let (Some(chrom), Some(pos), Some(r), Some(a)) = (
            text(h, &["chrom"]),
            int(get(h, &["pos"])).and_then(|p| u64::try_from(p).ok()),
            text(h, &["ref"]),
            text(h, &["alt"]),
        ) else {
            continue;
        };
        let at = format!("{chrom}:{pos} {r}>{a}");
        let sig = text(h, &["significance"]).unwrap_or_default();
        let lower = sig.to_ascii_lowercase();
        if !lower.contains("pathogenic") || lower.contains("conflicting") {
            continue;
        }
        let gene = text(h, &["gene"])
            .unwrap_or_default()
            .split(['|', ',', '&', ';'])
            .next()
            .unwrap_or_default()
            .to_owned();
        if !gene_ok(&gene) {
            let locus = format!("{chrom}:{pos}");
            ctx.leave_out(
                &locus,
                &format!("{sig} {at}"),
                &format!("no HGNC gene symbol ({gene:?})"),
            );
            continue;
        }
        let zygosity = text(h, &["genotype"]).unwrap_or("unknown");
        let vrs_id = match ctx.vrs_id(chrom, pos, r, a) {
            Ok(id) => id,
            Err(e) => {
                ctx.leave_out(
                    &gene,
                    &format!("{sig} {gene} {at}"),
                    &format!("no VRS identifier ({e})"),
                );
                continue;
            }
        };
        let carrier = ctx.panels.carrier_gene(&gene).cloned();
        let recessive = carrier.as_ref().is_some_and(|c| {
            c.inheritance == "AR" || (female && inheritance(&c.inheritance) == Some("XL"))
        });
        let category = if sf.contains(at.as_str()) {
            "acmg_sf"
        } else if recessive && zygosity == "het" && ctx.panels.in_carrier_span(chrom, pos) {
            "carrier"
        } else {
            "monogenic_pathogenic"
        };
        let group = ctx.panels.sensitive_group(&gene).map(str::to_owned);
        if group.is_some() && !ctx.panels.in_sensitive_span(chrom, pos) {
            ctx.leave_out(
                &gene,
                &format!("{sig} {gene} {at}"),
                &format!("outside OwnHealth's span of {gene}, which gates its lookups"),
            );
            continue;
        }
        let stars = int(get(h, &["stars"])).unwrap_or(0).clamp(0, 4);
        let level = match stars {
            2.. => "E1",
            1 => "E2",
            _ => "E3",
        };
        let (dp, gq) = (int(get(h, &["dp"])), int(get(h, &["gq"])));
        let low = get(h, &["low_confidence"]) == Some(&Value::Bool(true));
        let (conf, basis) = match (dp, gq) {
            _ if ctx.vcf_only => (
                "medium",
                "a call from a VCF given as input; the reads were not checked here".to_owned(),
            ),
            (Some(d), Some(g)) if !low && g >= floor_gq && d >= floor_dp => {
                ("high", format!("GQ {g}, DP {d}"))
            }
            (Some(d), Some(g)) => (
                "medium",
                format!("GQ {g}, DP {d}: below the call floor of GQ {floor_gq} and DP {floor_dp}"),
            ),
            _ => (
                "medium",
                "no depth or genotype quality on the call".to_owned(),
            ),
        };
        let mut genotype = obj([
            ("call_state", s("alt")),
            (
                "zygosity",
                s(if gt_of(zygosity).is_some() {
                    zygosity
                } else {
                    "unknown"
                }),
            ),
        ]);
        if let Some(g) = gt_of(zygosity) {
            put(&mut genotype, "gt", s(g));
        }
        if let Some(d) = dp {
            put(&mut genotype, "dp", Value::Int(d));
        }
        if let Some(g) = gq {
            put(&mut genotype, "gq", Value::Int(g));
        }
        let mut variant = obj([
            ("chrom", s(chrom)),
            ("pos", Value::Int(i128::from(pos))),
            ("ref", s(r)),
            ("alt", s(a)),
        ]);
        if let Some(rs) = list(h, &["rsids"]).iter().find_map(|x| match x {
            Value::Str(t)
                if t.starts_with("rs")
                    && t[2..].bytes().all(|b| b.is_ascii_digit())
                    && t.len() > 2 =>
            {
                Some(t.clone())
            }
            _ => None,
        }) {
            put(&mut variant, "rsid", Value::Str(rs));
        }
        put(&mut variant, "vrs_id", s(&vrs_id));
        if let Some(id) = variation.get(&format!("{chrom}:{pos}:{r}:{a}")) {
            put(&mut variant, "clinvar_variation_id", Value::Int(*id));
        }
        put(&mut variant, "genotype", genotype);
        let star_text = if stars == 1 {
            "1 star".to_owned()
        } else {
            format!("{stars} stars")
        };
        let summary_text = match category {
            "carrier" => format!(
                "Carrier of one {} {gene} variant ({}). Carriers of a recessive condition are usually not affected; \
                 it matters when a partner is tested.",
                lower,
                carrier.as_ref().map_or("", |c| c.condition.as_str())
            ),
            "acmg_sf" => format!(
                "{sig} {gene} variant ({zygosity}), a gene on the ACMG secondary-findings list v3.3. ClinVar \
                 review: {star_text}. A result from research-grade sequencing: a clinical test confirms it."
            ),
            _ => format!(
                "{sig} {gene} variant ({zygosity}). ClinVar review: {star_text}. A result from research-grade \
                 sequencing: a clinical test confirms it."
            ),
        };
        let mut fact = obj([
            ("id", s(&fact_id(category, &gene, &vrs_id))),
            ("category", s(category)),
            ("gene", s(&gene)),
            ("variants", Value::Array(vec![variant])),
            ("summary", s(&clean_text(&summary_text))),
        ]);
        let cond = match (&carrier, category) {
            (Some(c), "carrier" | "monogenic_pathogenic") => {
                let mut o = obj([("name", s(&clean_text(&c.condition)))]);
                if let Some(i) = inheritance(&c.inheritance) {
                    put(&mut o, "inheritance", s(i));
                }
                Some(o)
            }
            (_, "acmg_sf") => ctx
                .panels
                .acmg
                .iter()
                .find(|(g, _)| *g == gene)
                .and_then(|(_, i)| inheritance(i))
                .map(|i| obj([("inheritance", s(i))])),
            _ => None,
        };
        if let Some(c) = cond {
            put(&mut fact, "condition", c);
        }
        put(
            &mut fact,
            "evidence",
            obj([
                ("level", s(level)),
                (
                    "native",
                    s(&clean_text(&format!("ClinVar {sig}, {star_text}"))),
                ),
                ("clinvar_stars", Value::Int(stars)),
            ]),
        );
        put(
            &mut fact,
            "confidence",
            obj([("level", s(conf)), ("basis", s(&basis))]),
        );
        let sens = match &group {
            Some(g) => obj([("level", s("sensitive")), ("reason", s(g))]),
            None => obj([("level", s("normal"))]),
        };
        put(&mut fact, "sensitivity", sens);
        put(&mut fact, "confirm_clinically", Value::Bool(true));
        put(&mut fact, "origin", origin(summary, "clinvar", "06"));
        ctx.push(fact);
    }
}

/// The PGx calls of steps 07, 27, 32 and 36 (summary.json's `cpic` and
/// `cyp2d6`), with PharmCAT's CPIC and DPWG guidance for the called
/// diplotype.
fn pgx_facts(ctx: &mut Ctx<'_>, dir: &Path, summary: &Value) {
    let Some(cpic) = section(summary, "cpic") else {
        return;
    };
    let report_path = get(summary, &["sections", "pharmcat"])
        .and_then(|p| text(p, &["source"]))
        .map(|p| dir.join(p));
    let report = report_path
        .as_ref()
        .and_then(|p| std::fs::read_to_string(p).ok())
        .and_then(|t| json::loads(&t).ok());
    let (calls, guidance) = match &report {
        Some(r) => (
            match og_pgx::report::parse_genes(r) {
                og_pgx::report::Parsed::Ok(c) => c,
                _ => Vec::new(),
            },
            og_pgx::guidance::drug_guidance(r),
        ),
        None => (Vec::new(), Default::default()),
    };
    let pharmcat = section(summary, "pharmcat")
        .and_then(|d| text(d, &["version"]))
        .unwrap_or("unknown version")
        .to_owned();
    let low_gq: Vec<(String, String)> = section(summary, "pypgx")
        .map(|d| {
            list(d, &["low_gq"])
                .iter()
                .filter_map(|x| {
                    Some((
                        text(x, &["gene"])?.to_owned(),
                        text(x, &["weakest"])?.to_owned(),
                    ))
                })
                .collect()
        })
        .unwrap_or_default();
    for row in list(cpic, &["genes"]) {
        let Some(gene) = text(row, &["gene"]) else {
            continue;
        };
        if gene.starts_with("HLA-") || !gene_ok(gene) {
            continue; // HLA comes from T1K (hla_facts)
        }
        let status = text(row, &["status"]).unwrap_or_default();
        let mut diplotype = text(row, &["diplotype"]).unwrap_or_default().to_owned();
        let mut phenotype = text(row, &["phenotype"]).unwrap_or_default().to_owned();
        let mut tools = vec![s("PharmCAT")];
        let (mut conf, mut basis) = if gene == "CYP2D6" {
            // High only when two callers agree (below).
            (
                if ctx.vcf_only { "low" } else { "medium" },
                if ctx.vcf_only {
                    "CYP2D6 from a VCF: copy number and hybrids are not seen".to_owned()
                } else {
                    format!(
                        "PharmCAT {pharmcat}; not confirmed by two read-level callers that agree"
                    )
                },
            )
        } else {
            ("high", format!("PharmCAT {pharmcat} on the genome's calls"))
        };
        if gene == "CYP2D6" {
            let cyp = get(summary, &["cyp2d6"]);
            let consensus = cyp.and_then(|c| text(c, &["consensus", "result"]));
            let passed = cyp.and_then(|c| get(c, &["consensus", "passed_to_pharmcat"]))
                == Some(&Value::Bool(true));
            let probable = cyp
                .and_then(|c| get(c, &["probable"]))
                .filter(|p| matches!(p, Value::Object(_)));
            // Step 36 handed the agreed call to PharmCAT, or PharmCAT's own
            // call is the one the two read-level callers agree on.
            if let Some(c) = consensus
                && cyp.and_then(|c| get(c, &["agree"])) == Some(&Value::Bool(true))
                && (passed || (status != "not called" && same_diplotype(&diplotype, c)))
            {
                let reason = cyp
                    .and_then(|c| text(c, &["consensus", "reason"]))
                    .unwrap_or("two callers agree");
                if !ctx.vcf_only {
                    (conf, basis) = ("high", reason.to_owned());
                }
                tools.extend([s("pypgx"), s("Cyrius")]);
                if status == "not called" {
                    diplotype = c.to_owned();
                }
            } else if status == "not called"
                && let Some(p) = probable
            {
                diplotype = text(p, &["diplotype"]).unwrap_or_default().to_owned();
                phenotype = text(p, &["phenotype"]).unwrap_or_default().to_owned();
                let label = text(p, &["label"]).unwrap_or("probable, one caller");
                (conf, basis) = (
                    if ctx.vcf_only { "low" } else { "medium" },
                    label.to_owned(),
                );
                tools = vec![s("pypgx")];
            } else if status == "not called" {
                continue;
            }
        } else if status == "not called" {
            continue;
        }
        if let Some((_, weakest)) = low_gq.iter().find(|(g, _)| g == gene)
            && conf == "high"
        {
            (conf, basis) = (
                "medium",
                format!("a defining site has a low genotype quality: {weakest}"),
            );
        }
        // genomics.md 2.2 rule 1: a call without variants[] carries its own
        // call_state; a one-site gene also cites its site.
        let call_state = pgx_call_state(gene, &diplotype);
        let site = ctx.panels.pgx_site(gene).cloned();
        let shown = single_site(&diplotype);
        let mut variants = Vec::new();
        if let Some(site) = &site
            && let Some(bases) = shown.strip_prefix(&format!("{} ", site.rsid))
            && bases
                .split('/')
                .all(|b| !b.is_empty() && b.bytes().all(|x| b"ACGT".contains(&x)))
        {
            match site_variant(ctx, site, "called", bases) {
                Ok(v) => variants.push(v),
                Err(e) => ctx.warn(format!(
                    "{gene}: {} not cited, no VRS identifier ({e})",
                    site.rsid
                )),
            }
        }
        let mut value = obj([]);
        if term_ok(&shown, Some(gene)) {
            put(&mut value, "diplotype", s(&shown));
        } else if !shown.is_empty() {
            ctx.warn(format!(
                "{gene}: diplotype {shown:?} left out, not a term OwnHealth reads"
            ));
        }
        // PharmCAT's placeholders say it named no phenotype; OwnHealth then
        // reads the call_state instead of counting the word as a finding.
        if PHENOTYPE_PLACEHOLDERS.contains(&phenotype.trim().to_lowercase().as_str()) {
            phenotype.clear();
        }
        let pheno_ok = term_ok(&phenotype, Some(gene));
        if pheno_ok {
            put(&mut value, "phenotype", s(&phenotype));
        } else if !phenotype.is_empty() {
            ctx.warn(format!(
                "{gene}: phenotype {phenotype:?} left out, not a term OwnHealth reads"
            ));
        }
        put(&mut value, "call_state", s(call_state));
        // PharmCAT's guidance for the diplotype it called (none for a
        // CYP2D6 call PharmCAT did not make).
        let mut dg = Vec::new();
        let mut seen = HashSet::new();
        if let Some(call) = calls.iter().find(|c| c.gene == gene && c.called)
            && let Some(entries) = guidance.get(gene)
        {
            for g in og_pgx::guidance::matching(call, entries) {
                let src = match g.source {
                    og_pgx::guidance::Source::Cpic => "CPIC",
                    og_pgx::guidance::Source::Dpwg => "DPWG",
                    _ => continue,
                };
                let drug = g.drug.to_lowercase();
                if !drug_ok(&drug) || !seen.insert((drug.clone(), src)) {
                    continue;
                }
                let class = if g.class.trim().is_empty() {
                    "not rated"
                } else {
                    g.class.trim()
                };
                let mut o = obj([("drug", s(&drug)), ("source", s(src)), ("level", s(class))]);
                if pheno_ok {
                    put(&mut o, "phenotype", s(&phenotype));
                }
                if !g.rec.trim().is_empty() {
                    put(&mut o, "recommendation", s(&clean_text(&g.rec)));
                }
                dg.push(o);
            }
        }
        let on_cpic = ctx.panels.cpic.iter().any(|g| g == gene);
        let (level, native) = if on_cpic {
            ("E1", "CPIC guideline gene")
        } else {
            ("E2", "PharmCAT call; no CPIC guideline for this gene")
        };
        let what = [shown.as_str(), phenotype.as_str()]
            .iter()
            .filter(|x| !x.is_empty())
            .copied()
            .collect::<Vec<_>>()
            .join(": ");
        let mut fact = obj([
            ("id", s(&fact_id("pgx", gene, gene))),
            ("category", s("pgx")),
            ("gene", s(gene)),
        ]);
        if !variants.is_empty() {
            put(&mut fact, "variants", Value::Array(variants));
        }
        put(&mut fact, "value", value);
        put(
            &mut fact,
            "summary",
            s(&clean_text(&format!("{gene} {what}."))),
        );
        if !dg.is_empty() {
            put(&mut fact, "drug_guidance", Value::Array(dg));
        }
        put(
            &mut fact,
            "evidence",
            obj([
                ("level", s(level)),
                ("native", s(native)),
                ("tools", Value::Array(tools)),
            ]),
        );
        put(
            &mut fact,
            "confidence",
            obj([("level", s(conf)), ("basis", s(&clean_text(&basis)))]),
        );
        put(&mut fact, "sensitivity", obj([("level", s("normal"))]));
        put(&mut fact, "origin", origin(summary, "cpic", "27"));
        ctx.push(fact);
    }
}

/// The phenotypes PharmCAT writes when it names none.
const PHENOTYPE_PLACEHOLDERS: &[&str] = &[
    "n/a",
    "no result",
    "indeterminate",
    "unclassified",
    "no phenotype assigned",
];

/// The call_state of a star-allele call: `ref` when every allele is the
/// gene's reference allele as PharmCAT names it (`*1`, `Reference`, `B
/// (reference)`, `rsN reference (X)`, CYP2C19 `*38`, NAT2 `*4`), `no_call`
/// when an allele is unknown, else `alt`.
pub fn pgx_call_state(gene: &str, diplotype: &str) -> &'static str {
    let alleles: Vec<&str> = diplotype.split('/').map(str::trim).collect();
    let unknown = |a: &str| {
        a.is_empty() || a == "?" || matches!(a.to_lowercase().as_str(), "unknown" | "indeterminate")
    };
    let reference = |a: &str| {
        a.to_lowercase().contains("reference")
            || matches!(a, "*1" | "*1A" | "*1.001")
            || (gene == "CYP2C19" && a == "*38")
            || (gene == "NAT2" && a == "*4")
    };
    if alleles.iter().any(|a| unknown(a)) {
        "no_call"
    } else if alleles.iter().all(|a| reference(a)) {
        "ref"
    } else {
        "alt"
    }
}

/// PharmCAT's diplotype of a one-site gene, `rsN reference (X)/rsN variant
/// (Y)`, as `rsN X/Y`, a term OwnHealth reads; any other diplotype unchanged.
fn single_site(diplotype: &str) -> String {
    let half = |h: &str| -> Option<(String, String)> {
        let (rs, rest) = h.trim().split_once(' ')?;
        let base = rest.trim().strip_suffix(')')?.rsplit_once('(')?.1;
        let ok = rs.len() > 2
            && rs.starts_with("rs")
            && rs[2..].bytes().all(|b| b.is_ascii_digit())
            && !base.is_empty()
            && base.bytes().all(|b| b"ACGT".contains(&b));
        ok.then(|| (rs.to_owned(), base.to_owned()))
    };
    match diplotype.split_once('/').map(|(a, b)| (half(a), half(b))) {
        Some((Some((rs, a)), Some((rs2, b)))) if rs == rs2 => format!("{rs} {a}/{b}"),
        _ => diplotype.to_owned(),
    }
}

/// The drug-relevant HLA alleles T1K typed (step 08).
fn hla_facts(ctx: &mut Ctx<'_>, summary: &Value) {
    let Some(hla) = section(summary, "hla") else {
        return;
    };
    for a in list(hla, &["derived", "alleles"]) {
        if text(a, &["kind"]) != Some("drug") {
            continue;
        }
        let (Some(allele), Some(state)) = (text(a, &["allele"]), text(a, &["state"])) else {
            continue;
        };
        let gene = format!("HLA-{}", allele.split('*').next().unwrap_or_default());
        if !gene_ok(&gene) {
            continue;
        }
        let present = match state {
            "present" => true,
            "absent" => false,
            _ => {
                ctx.leave_out(
                    &gene,
                    &format!("HLA-{allele}"),
                    "T1K's typing could not say whether it is present",
                );
                continue;
            }
        };
        let marker = format!("HLA-{allele}");
        let locus = list(hla, &["loci"])
            .iter()
            .find(|l| text(l, &["gene"]) == Some(gene.as_str()));
        let alleles: Vec<Value> = locus
            .map(|l| list(l, &["alleles"]).to_vec())
            .unwrap_or_default();
        let low = locus.and_then(|l| get(l, &["low_confidence"])) == Some(&Value::Bool(true));
        let association = text(a, &["association"]).unwrap_or_default();
        let citation = text(a, &["citation"]).unwrap_or_default();
        let found = if present { "detected" } else { "not detected" };
        let pmids: Vec<Value> = text(a, &["pmids"])
            .unwrap_or_default()
            .split(',')
            .filter_map(|p| p.trim().parse::<i128>().ok().map(Value::Int))
            .collect();
        let mut evidence = obj([("level", s("E1")), ("native", s("CPIC guideline allele"))]);
        if !pmids.is_empty() {
            put(&mut evidence, "pmids", Value::Array(pmids));
        }
        put(&mut evidence, "tools", Value::Array(vec![s("T1K")]));
        let fact = obj([
            ("id", s(&fact_id("hla", &gene, &marker))),
            ("category", s("hla")),
            ("gene", s(&gene)),
            (
                "value",
                obj([
                    ("alleles", Value::Array(alleles)),
                    ("marker", s(&marker)),
                    ("present", Value::Bool(present)),
                    ("call_state", s(if present { "alt" } else { "ref" })),
                ]),
            ),
            (
                "summary",
                s(&clean_text(&format!(
                    "{marker} {found}. Linked to {association} ({citation})."
                ))),
            ),
            ("evidence", evidence),
            (
                "confidence",
                obj([
                    ("level", s(if low { "low" } else { "medium" })),
                    (
                        "basis",
                        s(if low {
                            "T1K typing from reads, flagged low confidence"
                        } else {
                            "T1K typing from reads"
                        }),
                    ),
                ]),
            ),
            ("sensitivity", obj([("level", s("normal"))])),
            ("origin", origin(summary, "hla", "08")),
        ]);
        ctx.push(fact);
    }
}

/// `{step, file}`: the step and the file a section of summary.json read.
fn origin(summary: &Value, name: &str, step: &str) -> Value {
    let file = get(summary, &["sections", name])
        .and_then(|sec| text(sec, &["source"]))
        .unwrap_or_default();
    obj([("step", s(step)), ("file", s(file))])
}

/// The polygenic scores of step 25.
fn prs_facts(ctx: &mut Ctx<'_>, summary: &Value) {
    let Some(prs) = section(summary, "prs") else {
        return;
    };
    for row in list(prs, &["scores"]) {
        let (Some(pgs), Some(pct)) = (text(row, &["pgs_id"]), number(get(row, &["percentile"])))
        else {
            continue;
        };
        if !pgs.starts_with("PGS") {
            continue;
        }
        let condition = text(row, &["condition"]).unwrap_or_default();
        let group = text(row, &["group"]).unwrap_or("reference");
        let evaluated: Vec<Value> = list(row, &["evaluated_in"]).to_vec();
        let inside = evaluated.iter().any(|e| e == &s(group));
        let matched = match (number(get(row, &["matched"])), number(get(row, &["total"]))) {
            (Some(m), Some(t)) if t > 0.0 => Some((m / t * 1000.0).round() / 10.0),
            _ => None,
        };
        let conf = if matched.is_some_and(|m| m >= 90.0) && inside {
            "medium"
        } else {
            "low"
        };
        let basis = format!(
            "matched {} of the score's variants; {group} ancestry {} the groups it was evaluated in",
            matched.map_or("an unknown share".to_owned(), |m| format!("{m}%")),
            if inside { "inside" } else { "outside" }
        );
        let mut value = obj([
            ("condition", s(&clean_text(condition))),
            ("pgs_id", s(pgs)),
            ("percentile", Value::Float(pct)),
            ("most_similar_pop", s(group)),
            ("evaluated_in", Value::Array(evaluated)),
        ]);
        if let Some(m) = matched {
            put(&mut value, "matched_pct", Value::Float(m));
        }
        let adult = condition.to_lowercase().contains("alzheimer");
        let pct_text = pycompat::fmt::repr(pct);
        let mut fact = obj([
            ("id", s(&fact_id("prs", "", pgs))),
            ("category", s("prs")),
            ("value", value),
            (
                "summary",
                s(&clean_text(&format!(
                    "{condition}: higher than {pct_text}% of the {group} reference group for this score."
                ))),
            ),
            (
                "evidence",
                obj([
                    ("level", s("E3")),
                    ("native", s(&format!("PGS Catalog score {pgs}"))),
                    ("tools", Value::Array(vec![s("pgsc_calc")])),
                ]),
            ),
            (
                "confidence",
                obj([("level", s(conf)), ("basis", s(&basis))]),
            ),
            (
                "limitations",
                Value::Array(vec![s(&format!(
                    "A percentile within the {group} reference group, not an absolute risk."
                ))]),
            ),
        ]);
        put(
            &mut fact,
            "sensitivity",
            if adult {
                obj([
                    ("level", s("adult_onset")),
                    ("reason", s("Alzheimer's disease score")),
                ])
            } else {
                obj([("level", s("normal"))])
            },
        );
        put(&mut fact, "origin", origin(summary, "prs", "25"));
        ctx.push(fact);
    }
}

/// `called` with genotype `G/A` at a site of `ref>alt` as a VCF genotype.
fn site_genotype(site: &Site, status: &str, genotype: &str) -> Value {
    let state = match status {
        "called" => {
            let copies = genotype
                .split(['/', '|'])
                .filter(|b| *b == site.alt)
                .count();
            if copies == 0 { "ref" } else { "alt" }
        }
        "not_covered" | "not covered" => "not_covered",
        _ => "no_call",
    };
    let mut g = obj([("call_state", s(state))]);
    if status == "called" {
        let alleles: Vec<&str> = genotype.split(['/', '|']).collect();
        let copies = alleles.iter().filter(|b| **b == site.alt).count();
        let gt: Vec<&str> = alleles
            .iter()
            .map(|b| if *b == site.alt { "1" } else { "0" })
            .collect();
        put(&mut g, "gt", s(&gt.join("/")));
        let zyg = match (alleles.len(), copies) {
            (1, 0) | (2, 0) => "ref",
            (1, _) => "hemi",
            (_, 1) => "het",
            _ => "hom",
        };
        put(&mut g, "zygosity", s(zyg));
    }
    g
}

/// Two diplotypes name the same pair of alleles, in either order.
fn same_diplotype(a: &str, b: &str) -> bool {
    fn pair(d: &str) -> Vec<&str> {
        let mut v: Vec<&str> = d.split('/').map(str::trim).collect();
        v.sort_unstable();
        v
    }
    pair(a) == pair(b)
}

fn site_variant(
    ctx: &mut Ctx<'_>,
    site: &Site,
    status: &str,
    genotype: &str,
) -> Result<Value, String> {
    let vrs_id = ctx.vrs_id(&site.chrom, site.pos, &site.ref_, &site.alt)?;
    Ok(obj([
        ("chrom", s(&site.chrom)),
        ("pos", Value::Int(i128::from(site.pos))),
        ("ref", s(&site.ref_)),
        ("alt", s(&site.alt)),
        ("rsid", s(&site.rsid)),
        ("vrs_id", s(&vrs_id)),
        ("genotype", site_genotype(site, status, genotype)),
    ]))
}

/// APOE from step 42's rule (traits section): the e2/e3/e4 pair.
fn apoe_fact(ctx: &mut Ctx<'_>, summary: &Value) -> bool {
    let Some(traits) = section(summary, "traits") else {
        return false;
    };
    let rows = list(traits, &["apoe"]);
    let rule = rows
        .iter()
        .find(|r| text(r, &["rsid"]).is_some_and(|x| x.contains('+')));
    let Some(pair) = rule
        .filter(|r| text(r, &["status"]) == Some("called"))
        .and_then(|r| text(r, &["genotype"]))
        .filter(|g| term_ok(g, Some("APOE")))
        .map(str::to_owned)
    else {
        return false;
    };
    let mut variants = Vec::new();
    for site in ctx.panels.apoe.clone() {
        let row = rows
            .iter()
            .find(|r| text(r, &["rsid"]) == Some(site.rsid.as_str()));
        let status = row
            .and_then(|r| text(r, &["status"]))
            .unwrap_or("not_covered");
        let genotype = row.and_then(|r| text(r, &["genotype"])).unwrap_or_default();
        match site_variant(ctx, &site, status, genotype) {
            Ok(v) => variants.push(v),
            Err(e) => {
                ctx.warn(format!(
                    "APOE: left out, no VRS identifier for {} ({e})",
                    site.rsid
                ));
                return false;
            }
        }
    }
    let e4 = pair.contains("e4");
    let reading = rule.and_then(|r| text(r, &["reading"])).unwrap_or_default();
    let fact = obj([
        ("id", s(&fact_id("risk_marker", "APOE", "APOE"))),
        ("category", s("risk_marker")),
        ("gene", s("APOE")),
        ("variants", Value::Array(variants)),
        ("value", obj([("haplotype", s(&pair))])),
        (
            "summary",
            s(&clean_text(&format!(
                "APOE {pair}. {reading} It is not a diagnosis."
            ))),
        ),
        (
            "evidence",
            obj([
                ("level", s("E2")),
                ("native", s("replicated association with effect size")),
                (
                    "pmids",
                    Value::Array(vec![Value::Int(8346443), Value::Int(9343467)]),
                ),
            ]),
        ),
        (
            "confidence",
            obj([
                ("level", s(if ctx.vcf_only { "medium" } else { "high" })),
                ("basis", s("both sites read from the gVCF")),
            ]),
        ),
        // Every APOE result is adult-onset, e4 or not: a hidden row only
        // when e4 is there would itself say e4 is there.
        (
            "sensitivity",
            obj([
                ("level", s("adult_onset")),
                ("reason", s(if e4 { "APOE e4" } else { "APOE result" })),
            ]),
        ),
        ("origin", origin(summary, "traits", "42")),
    ]);
    ctx.push(fact);
    true
}

/// Lp(a) from step 46 (the KIV-2 count) and step 42's two LPA sites. Only
/// what was read goes in: a site with no step 42 row is not cited, the risk
/// allele count needs a called site, and with neither a KIV-2 count nor a
/// called site there is no fact. `None` when written, else why not.
fn lpa_fact(ctx: &mut Ctx<'_>, summary: &Value) -> Option<&'static str> {
    let Some(lpa) = section(summary, "lpa") else {
        return Some("the Lp(a) step (46) did not run");
    };
    let mut variants = Vec::new();
    let mut called = Vec::new();
    for site in ctx.panels.lpa.clone() {
        let Some(row) = list(lpa, &["variants"])
            .iter()
            .find(|r| text(r, &["rsid"]) == Some(site.rsid.as_str()))
        else {
            continue;
        };
        let status = text(row, &["status"]).unwrap_or_default();
        let genotype = text(row, &["genotype"]).unwrap_or_default();
        if status == "called" {
            called.push(site.rsid.clone());
        }
        match site_variant(ctx, &site, status, genotype) {
            Ok(v) => variants.push(v),
            Err(e) => ctx.leave_out(
                "LPA",
                &format!("LPA {}", site.rsid),
                &format!("no VRS identifier ({e})"),
            ),
        }
    }
    let copies = (text(lpa, &["status"]) == Some("estimated"))
        .then(|| number(get(lpa, &["kiv2_copies"])))
        .flatten();
    if copies.is_none() && called.is_empty() {
        return Some(
            "neither the KIV-2 count (step 46) nor the two LPA sites (step 42) could be read",
        );
    }
    let mut value = obj([]);
    let mut parts = Vec::new();
    if !called.is_empty() {
        let n = int(get(lpa, &["risk_alleles"])).unwrap_or(0);
        put(&mut value, "risk_alleles", Value::Int(n));
        parts.push(format!(
            "{n} risk allele{} at {}",
            if n == 1 { "" } else { "s" },
            called.join(" and ")
        ));
    }
    if let Some(c) = copies {
        put(&mut value, "kiv2_copies", Value::Float(c));
        let (lo, hi) = (
            number(get(lpa, &["ci95_low"])),
            number(get(lpa, &["ci95_high"])),
        );
        let mut kiv2 = format!("KIV-2 about {c} copies");
        if let (Some(lo), Some(hi)) = (lo, hi) {
            put(&mut value, "kiv2_ci95_low", Value::Float(lo));
            put(&mut value, "kiv2_ci95_high", Value::Float(hi));
            kiv2.push_str(&format!(" (95% interval {lo} to {hi})"));
        }
        if text(lpa, &["class"]) == Some("low") {
            kiv2.push_str(", a low count; low counts go with a high Lp(a)");
        }
        parts.push(kiv2);
    } else {
        parts.push("KIV-2 count not estimated".to_owned());
    }
    let mut fact = obj([
        ("id", s(&fact_id("risk_marker", "LPA", "LPA"))),
        ("category", s("risk_marker")),
        ("gene", s("LPA")),
    ]);
    if !variants.is_empty() {
        put(&mut fact, "variants", Value::Array(variants));
    }
    put(&mut fact, "value", value);
    put(
        &mut fact,
        "summary",
        s(&clean_text(&format!("LPA: {}.", parts.join("; ")))),
    );
    put(
        &mut fact,
        "evidence",
        obj([
            ("level", s("E2")),
            ("native", s("replicated association with effect size")),
            ("pmids", Value::Array(vec![Value::Int(20032323)])),
        ]),
    );
    let basis = if copies.is_some() {
        "KIV-2 count is a read-depth estimate with an interval"
    } else {
        "the two LPA sites only; the KIV-2 count was not estimated"
    };
    put(
        &mut fact,
        "confidence",
        obj([("level", s("medium")), ("basis", s(basis))]),
    );
    put(
        &mut fact,
        "limitations",
        Value::Array(vec![s(
            "A research estimate from read depth; Lp(a) itself is measured in blood.",
        )]),
    );
    put(&mut fact, "sensitivity", obj([("level", s("normal"))]));
    put(&mut fact, "origin", origin(summary, "lpa", "46"));
    ctx.push(fact);
    None
}

/// One coverage row.
fn cov(
    id: &str,
    panel: &str,
    state: &str,
    reason: Option<&str>,
    fraction: Option<f64>,
    note: Option<String>,
) -> Value {
    let mut o = obj([
        ("scope", s("gene")),
        ("id", s(id)),
        ("panel", s(panel)),
        ("state", s(state)),
    ]);
    if let Some(r) = reason {
        put(&mut o, "reason", s(r));
    }
    if let Some(f) = fraction {
        put(&mut o, "fraction_callable", Value::Float(f));
    }
    if let Some(n) = note {
        put(&mut o, "note", s(&clean_text(&n)));
    }
    o
}

/// The callability catalogue's verdict on a gene short reads did not read.
fn callability_row(id: &str, panel: &str, e: &Value) -> Value {
    let share = number(get(e, &["share"]));
    let reasons: Vec<&str> = list(e, &["reasons"])
        .iter()
        .filter_map(|r| match r {
            Value::Str(t) => Some(t.as_str()),
            _ => None,
        })
        .collect();
    let owner = text(e, &["owner"]);
    let Some(share) = share else {
        return cov(
            id,
            panel,
            "not_assessed",
            Some("other"),
            None,
            Some("not on the reference this run used".to_owned()),
        );
    };
    let reason = if reasons
        .iter()
        .any(|r| matches!(*r, "segdup" | "collapsed" | "false_duplication"))
    {
        "paralog"
    } else if reasons.contains(&"ALT-only") {
        "other"
    } else {
        "low_depth"
    };
    let mut note = format!(
        "{:.0}% of the gene had 10 or more reads of mapping quality 20 or more",
        share * 100.0
    );
    if !reasons.is_empty() {
        note.push_str(&format!(" ({})", reasons.join(", ")));
    }
    if let Some(o) = owner {
        note.push_str(&format!("; called by its own method ({o})"));
    }
    let state = if share > 0.0 || owner.is_some() {
        "partial"
    } else {
        "not_assessed"
    };
    cov(id, panel, state, Some(reason), Some(share), Some(note))
}

/// Every gene of the fixed panels once (ACMG SF, then CPIC, then carrier),
/// then every catalogue gene the run could not read.
fn coverage(ctx: &Ctx<'_>, summary: &Value) -> Vec<Value> {
    let call = section(summary, "callability");
    let mut unread: BTreeMap<String, &Value> = BTreeMap::new();
    if let Some(c) = call {
        for e in list(c, &["not_read"])
            .iter()
            .chain(list(c, &["not_read_owned"]))
        {
            if let Some(g) = text(e, &["gene"]) {
                unread.insert(g.to_owned(), e);
            }
        }
    }
    let clinvar = section(summary, "clinvar").is_some();
    let cpic = section(summary, "cpic");
    let hla = section(summary, "hla");
    let cyp2d6_probable =
        get(summary, &["cyp2d6", "probable"]).is_some_and(|p| matches!(p, Value::Object(_)));
    // A gene counts as assessed by the ClinVar screen only with evidence
    // that its reads were there: the callability catalogue ran, and the
    // genome's mean depth reached its depth floor.
    let min_depth = call
        .and_then(|c| number(get(c, &["min_depth"])))
        .unwrap_or(10.0);
    let depth = section(summary, "coverage").and_then(|c| number(get(c, &["mean_depth"])));
    let unproven: Option<(&str, String)> = match (call, depth) {
        (None, _) => Some((
            "not_run",
            "gene callability (step 16b) did not run, so whether this gene was read is not known"
                .to_owned(),
        )),
        (_, None) => Some((
            "not_run",
            "the coverage step (16) did not report a mean depth, so whether this gene was read is not known"
                .to_owned(),
        )),
        (_, Some(d)) if d < min_depth => Some((
            "low_depth",
            format!("the genome's mean depth is {d}, below the {min_depth} reads a gene needs"),
        )),
        _ => None,
    };
    let mut rows = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let screen = |gene: &str, panel: &str| -> Value {
        if !clinvar {
            return cov(
                gene,
                panel,
                "not_assessed",
                Some("not_run"),
                None,
                Some("the ClinVar screen (step 06) did not run".to_owned()),
            );
        }
        match (unread.get(gene), &unproven) {
            (Some(e), _) => callability_row(gene, panel, e),
            (None, Some((reason, note))) => cov(
                gene,
                panel,
                "not_assessed",
                Some(reason),
                None,
                Some(note.clone()),
            ),
            (None, None) => cov(gene, panel, "assessed", None, None, None),
        }
    };
    for (gene, _) in &ctx.panels.acmg {
        if seen.insert(gene.clone()) {
            rows.push(screen(gene, ACMG_SF));
        }
    }
    for gene in &ctx.panels.cpic {
        if !seen.insert(gene.clone()) {
            continue;
        }
        let row = if gene.starts_with("HLA-") {
            let typed = hla.is_some_and(|h| {
                list(h, &["loci"]).iter().any(|l| {
                    text(l, &["gene"]) == Some(gene.as_str()) && !list(l, &["alleles"]).is_empty()
                })
            });
            if typed {
                cov(gene, "cpic", "assessed", None, None, None)
            } else {
                cov(
                    gene,
                    "cpic",
                    "not_assessed",
                    Some("not_run"),
                    None,
                    Some("HLA typing (step 08) did not type this gene".to_owned()),
                )
            }
        } else if let Some(c) = cpic {
            let r = list(c, &["genes"])
                .iter()
                .find(|r| text(r, &["gene"]) == Some(gene.as_str()));
            match r.and_then(|r| text(r, &["status"])) {
                None => cov(
                    gene,
                    "cpic",
                    "not_assessed",
                    Some("not_run"),
                    None,
                    Some("PharmCAT does not call this gene".to_owned()),
                ),
                Some("not called") if gene == "CYP2D6" && cyp2d6_probable => cov(
                    gene,
                    "cpic",
                    "partial",
                    Some("other"),
                    None,
                    Some("a probable call from one caller; PharmCAT made none".to_owned()),
                ),
                Some("not called") => cov(
                    gene,
                    "cpic",
                    "not_assessed",
                    Some("required_positions_missing"),
                    None,
                    Some("PharmCAT made no call".to_owned()),
                ),
                Some(_) => cov(gene, "cpic", "assessed", None, None, None),
            }
        } else {
            cov(
                gene,
                "cpic",
                "not_assessed",
                Some("not_run"),
                None,
                Some("the PGx steps (07, 27) did not run".to_owned()),
            )
        };
        rows.push(row);
    }
    for c in &ctx.panels.carrier {
        if seen.insert(c.span.gene.clone()) {
            rows.push(screen(&c.span.gene, "carrier"));
        }
    }
    if call.is_none() {
        rows.push(obj([
            ("scope", s("panel")),
            ("id", s("callability")),
            ("panel", s("callability")),
            ("state", s("not_assessed")),
            ("reason", s("not_run")),
            ("note", s("gene callability (step 16b) did not run, so hard genes the run could not read are not listed")),
        ]));
    }
    for (gene, e) in &unread {
        if seen.insert(gene.clone()) {
            rows.push(callability_row(gene, "callability", e));
        }
    }
    // A gene with a PGx or HLA call that no panel names (CYP3A4 has DPWG
    // guidance only) gets a row, so it never reads "not assessed" beside
    // its call; a low-confidence call is partial.
    for f in &ctx.facts {
        let (Some(gene), Some("pgx" | "hla")) = (text(f, &["gene"]), text(f, &["category"])) else {
            continue;
        };
        if seen.insert(gene.to_owned()) {
            rows.push(match text(f, &["confidence", "level"]) {
                Some("low") => cov(
                    gene,
                    "other",
                    "partial",
                    Some("other"),
                    None,
                    text(f, &["confidence", "basis"]).map(str::to_owned),
                ),
                _ => cov(gene, "other", "assessed", None, None, None),
            });
        }
    }
    // A call the file leaves out makes its gene partial, never "assessed,
    // no finding"; a gene or locus no panel names gets a row of its own.
    for (id, note) in &ctx.left_out {
        let row = rows
            .iter_mut()
            .find(|r| text(r, &["id"]) == Some(id.as_str()));
        match row {
            Some(r) => {
                if text(r, &["state"]) == Some("assessed") {
                    put(r, "state", s("partial"));
                    put(r, "reason", s("other"));
                }
                let n = match text(r, &["note"]) {
                    Some(old) => format!("{old}; {note}"),
                    None => note.clone(),
                };
                put(r, "note", s(&clean_text(&n)));
            }
            None => {
                let mut r = cov(
                    id,
                    "other",
                    "partial",
                    Some("other"),
                    None,
                    Some(note.clone()),
                );
                if id.contains(':') {
                    put(&mut r, "scope", s("locus"));
                }
                rows.push(r);
            }
        }
    }
    rows
}

/// A locus the run did not assess, as a coverage row.
fn locus(id: &str, why: &str) -> Value {
    let reason = if why.contains("did not run") {
        "not_run"
    } else {
        "other"
    };
    obj([
        ("scope", s("locus")),
        ("id", s(id)),
        ("panel", s("other")),
        ("state", s("not_assessed")),
        ("reason", s(reason)),
        ("note", s(why)),
    ])
}

/// The database versions of the run: the manifest's `data` rows, PharmCAT's
/// version and data version, and the ACMG SF list.
fn db_versions(summary: &Value, dir: &Path) -> Value {
    let mut m = obj([]);
    for row in list(summary, &["manifest"]) {
        if let Value::Array(c) = row
            && let [Value::Str(kind), Value::Str(key), Value::Str(val)] = c.as_slice()
            && kind == "data"
        {
            put(&mut m, key, s(val));
        }
    }
    if let Some(Value::Str(c)) = get(&m, &["clinvar_file_date"]).cloned() {
        put(
            &mut m,
            "clinvar",
            s(c.split_whitespace().next().unwrap_or_default()),
        );
    }
    if let Some(v) = section(summary, "pharmcat").and_then(|d| text(d, &["version"])) {
        put(&mut m, "pharmcat", s(v));
    }
    let report = get(summary, &["sections", "pharmcat", "source"])
        .and_then(|p| match p {
            Value::Str(p) => std::fs::read_to_string(dir.join(p)).ok(),
            _ => None,
        })
        .and_then(|t| json::loads(&t).ok());
    if let Some(dv) = report.as_ref().and_then(|r| text(r, &["dataVersion"])) {
        put(&mut m, "pharmcat_data", s(dv));
        // PharmCAT's data release carries the CPIC tables it quotes.
        put(&mut m, "cpic", s(dv.get(..10).unwrap_or(dv)));
    }
    put(&mut m, "acmg_sf", s("v3.3"));
    m
}

/// The calls VCF and gVCF with their tabix indexes, by hash.
fn files(sample: &str, dir: &Path) -> Result<Vec<Value>, String> {
    let mut out = Vec::new();
    for (role, rel) in [
        ("calls_vcf", format!("vcf/{sample}.vcf.gz")),
        ("gvcf", format!("vcf/{sample}.g.vcf.gz")),
    ] {
        let (path, index) = (dir.join(&rel), dir.join(format!("{rel}.tbi")));
        if !path.is_file() || !index.is_file() {
            continue;
        }
        let (sha, bytes) = file_sha256(&path)?;
        let (isha, _) = file_sha256(&index)?;
        out.push(obj([
            ("role", s(role)),
            ("path", s(&rel)),
            ("sha256", s(&sha)),
            ("index_path", s(&format!("{rel}.tbi"))),
            ("index_sha256", s(&isha)),
            ("bytes", Value::Int(i128::from(bytes))),
        ]));
    }
    Ok(out)
}

/// The facts document of the sample folder `dir`.
pub fn build(
    sample: &str,
    dir: &Path,
    panels: &Panels,
    genome: &mut dyn Genome,
) -> Result<Built, String> {
    let path = dir.join("summary.json");
    let raw = std::fs::read_to_string(&path)
        .map_err(|e| format!("{}: {e} (step 24 writes it)", path.display()))?;
    let summary = json::loads(&raw).map_err(|e| format!("{}: not JSON ({e:?})", path.display()))?;
    let declared = match text(&summary, &["run", "declared_sex"]) {
        Some("female") => "female",
        Some("male") => "male",
        _ => "unknown",
    };
    // Reads were there when a section that needs the BAM is: the read
    // count (alignment), the depth (mosdepth) or the callability catalogue.
    // A run that starts from a BAM has no read count but has the other two.
    let reads = ["read_count", "coverage", "callability"]
        .iter()
        .any(|k| section(&summary, k).is_some());
    let mut ctx = Ctx {
        panels,
        genome,
        facts: Vec::new(),
        ids: HashSet::new(),
        warnings: Vec::new(),
        vcf_only: !reads,
        left_out: Vec::new(),
    };
    clinvar_facts(&mut ctx, sample, dir, &summary, declared == "female");
    pgx_facts(&mut ctx, dir, &summary);
    hla_facts(&mut ctx, &summary);
    prs_facts(&mut ctx, &summary);
    let apoe = apoe_fact(&mut ctx, &summary);
    let lpa = lpa_fact(&mut ctx, &summary);
    let mut cover = coverage(&ctx, &summary);
    if !apoe {
        cover.push(locus(
            "APOE",
            "the traits panel (step 42) did not read both APOE sites",
        ));
    }
    if let Some(why) = lpa {
        cover.push(locus("LPA", why));
    }
    if section(&summary, "prs").is_none() {
        cover.push(obj([
            ("scope", s("panel")),
            ("id", s("prs")),
            ("panel", s("other")),
            ("state", s("not_assessed")),
            ("reason", s("not_run")),
            ("note", s("polygenic scores (step 25) did not run")),
        ]));
    }
    let generated = text(&summary, &["generated_utc"])
        .and_then(iso_utc)
        .ok_or_else(|| {
            format!(
                "{}: generated_utc is not 'YYYY-MM-DD HH:MM UTC'",
                path.display()
            )
        })?;
    let manifest = |key: &str| -> Option<String> {
        list(&summary, &["manifest"]).iter().find_map(|r| match r {
            Value::Array(c) => match c.as_slice() {
                [Value::Str(k), Value::Str(n), Value::Str(v)] if k == "run" && n == key => {
                    Some(v.clone())
                }
                _ => None,
            },
            _ => None,
        })
    };
    let mut provenance = obj([
        ("pipeline", s("OwnGenome")),
        (
            "version",
            s(&manifest("pipeline_version").unwrap_or_else(|| "unknown".to_owned())),
        ),
    ]);
    if let Some(c) = manifest("git_commit") {
        put(&mut provenance, "commit", s(&c));
    }
    put(&mut provenance, "generated_utc", s(&generated));
    put(&mut provenance, "db_versions", db_versions(&summary, dir));
    let mut assay = obj([
        ("type", s(if reads { "wgs" } else { "vcf_only" })),
        ("build", s("GRCh38")),
    ]);
    if let Some(d) = section(&summary, "coverage").and_then(|c| number(get(c, &["mean_depth"]))) {
        put(&mut assay, "mean_depth", Value::Float(d));
    }
    let doc = obj([
        ("schema", s("ownhealth.genome-facts/v1")),
        (
            "subject",
            obj([("sample_id", s(sample)), ("declared_sex", s(declared))]),
        ),
        ("assay", assay),
        ("provenance", provenance),
        ("files", Value::Array(files(sample, dir)?)),
        ("coverage", Value::Array(cover)),
        ("facts", Value::Array(std::mem::take(&mut ctx.facts))),
    ]);
    Ok(Built {
        doc,
        warnings: ctx.warnings,
    })
}
