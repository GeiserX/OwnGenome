//! The demo sample's medicines chapter: T1K, pypgx, Cyrius and the CYP2D6
//! depth check (steps 08, 32, 21), the PharmCAT report step 07 writes once
//! step 36 passed their consensus on, step 27's CPIC tables, and the
//! Locityper second opinion (step 53).

use super::{Ctx, Rng, Scratch, arg, copy, put, read, versions_env};
use pycompat::json::{self, Dump, Value};
use std::collections::HashMap;

/// `v[k]` of a JSON object.
fn get<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(k),
        _ => None,
    }
}

/// `v[k]` of a JSON object, to change.
fn get_mut<'a>(v: &'a mut Value, k: &str) -> Option<&'a mut Value> {
    match v {
        Value::Object(m) => m.get_mut(k),
        _ => None,
    }
}

fn s(x: &str) -> Value {
    Value::Str(x.to_owned())
}

fn obj(pairs: Vec<(&str, Value)>) -> Value {
    Value::Object(pairs.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}

/// The release of `HLA_DB_RELEASE`'s default in `scripts/lib/common.sh`.
fn hla_release(ctx: &Ctx) -> Result<String, String> {
    let text = read(&ctx.asset("scripts/lib/common.sh"))?;
    text.lines()
        .find_map(|l| {
            let rest = l.trim().strip_prefix("HLA_DB_RELEASE=${HLA_DB_RELEASE:-")?;
            let (release, _) = rest.split_once('}')?;
            (!release.is_empty() && release.bytes().all(|b| b.is_ascii_digit() || b == b'.'))
                .then(|| release.to_owned())
        })
        .ok_or_else(|| "HLA_DB_RELEASE is not in scripts/lib/common.sh".to_owned())
}

/// The tag of an image line of `versions.env` (PHARMCAT_IMAGE -> 3.4.0).
fn image_tag(ctx: &Ctx, var: &str) -> Result<String, String> {
    let image = versions_env(&ctx.repo, var)?;
    let name = image.split('@').next().unwrap_or("");
    Ok(name.rsplit(':').next().unwrap_or("").to_owned())
}

/// T1K, pypgx, Cyrius and the depth check, then the PharmCAT report that
/// step 07 would write once step 36 passed their consensus on, then step
/// 27's tables from it with `og pgx cpic-report`.
///
/// The report starts from PharmCAT's own example report
/// (`tests/fixtures/pharmcat/pharmcat-docs-example.json`, PharmCAT's
/// invented example sample); HLA-A, HLA-B and CYP2D6 become the outside
/// calls step 36 passes on, worked out by `og pgx outside-calls` itself.
/// T1K's table and the depth check are the invented fixtures of
/// `tests/fixtures/pgx`.
pub fn pgx(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let sample = ctx.sample.as_str();
    let example_path = ctx.asset("tests/fixtures/pharmcat/pharmcat-docs-example.json");
    let example = json::loads(&read(&example_path)?)
        .map_err(|e| format!("{}: {e}", example_path.display()))?;
    let genes = get(&example, "genes").ok_or("the PharmCAT example has no genes")?;
    // The example's first source diplotype of a gene: its label and its phenotypes, "; "-joined.
    let called = |gene: &str| -> Result<(String, String), String> {
        let first = get(genes, gene)
            .and_then(|g| get(g, "sourceDiplotypes"))
            .and_then(|d| match d {
                Value::Array(a) => a.first(),
                _ => None,
            })
            .ok_or_else(|| format!("the PharmCAT example has no diplotype of {gene}"))?;
        let label = match get(first, "label") {
            Some(Value::Str(x)) => x.clone(),
            _ => String::new(),
        };
        let phen = match get(first, "phenotypes") {
            Some(Value::Array(a)) => a
                .iter()
                .filter_map(|p| match p {
                    Value::Str(x) => Some(x.as_str()),
                    _ => None,
                })
                .collect::<Vec<_>>()
                .join("; "),
            _ => String::new(),
        };
        Ok((label, phen))
    };

    // T1K (step 08) and the CYP2D6 depth check: the invented fixtures.
    let t1k = ctx.path(&format!("hla_t1k/{sample}_hla_genotype.tsv"));
    copy(&ctx.asset("tests/fixtures/pgx/t1k_genotype.tsv"), &t1k)?;
    let depth = ctx.path(&format!("pypgx/{sample}_cyp2d6_depth_check.tsv"));
    copy(&ctx.asset("tests/fixtures/pgx/depth_check_ok.tsv"), &depth)?;
    // The HLA database release the manifest prints: one header line of
    // hla.dat, naming the release setup.sh installs.
    put(
        &ctx.genome_dir.join("t1k_idx/hlaidx/hla.dat"),
        &format!(
            "ID   HLA00001; standard; DNA; HUM; 3503 BP.\nDT   15/07/2025 (Rel. {}, Created)\n",
            hla_release(ctx)?
        ),
    )?;

    // pypgx (step 32) and Cyrius (step 21) agree with PharmCAT's example
    // calls, CYP2D6 included, so the three CYP2D6 callers agree.
    let (cyp, cyp_phen) = called("CYP2D6")?;
    let mut rows = String::from("Gene\tDiplotype\tPhenotype\tCNV_call\tSource\n");
    for gene in [
        "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "NUDT15", "SLCO1B1", "TPMT", "UGT1A1",
    ] {
        let (dip, phen) = called(gene)?;
        let (cnv, source) = if gene == "CYP2D6" {
            ("Normal", "BAM")
        } else {
            ("N/A", "VCF")
        };
        rows.push_str(&format!("{gene}\t{dip}\t{phen}\t{cnv}\t{source}\n"));
    }
    let pypgx = ctx.path(&format!("pypgx/{sample}_pypgx_summary.tsv"));
    put(&pypgx, &rows)?;
    let cyrius = ctx.path(&format!("cyrius/{sample}_cyp2d6.tsv"));
    put(
        &cyrius,
        &format!("Sample\tGenotype\tFilter\n{sample}\t{cyp}\tPASS\n"),
    )?;

    // What step 36 will pass to PharmCAT, from its own code.
    let tmp = Scratch::new("pgx")?;
    let calls = tmp.join("calls.tsv");
    let consensus = tmp.join("consensus.tsv");
    ctx.og(&[
        "pgx",
        "outside-calls",
        "--calls",
        arg(&calls),
        "--consensus",
        arg(&consensus),
        "--hla",
        arg(&t1k),
        "--pypgx",
        arg(&pypgx),
        "--cyrius",
        arg(&cyrius),
        "--depth-check",
        arg(&depth),
    ])?;
    let mut outside: Vec<(String, String)> = Vec::new();
    for line in read(&calls)?.lines() {
        let cells: Vec<&str> = line.split('\t').collect();
        let [gene, dip] = cells[..] else {
            return Err(format!("og pgx outside-calls wrote {line:?}"));
        };
        outside.retain(|(g, _)| g != gene);
        outside.push((gene.to_owned(), dip.to_owned()));
    }

    // The PharmCAT report (step 07, run after step 36): PharmCAT's example,
    // with this sample's outside calls. HLA phenotypes as PharmCAT writes them.
    let mut report = example.clone();
    if let Value::Object(m) = &mut report {
        m.shift_remove("_fixture_note");
        m.insert("title".to_owned(), s(sample));
        m.insert("timestamp".to_owned(), s("2026-10-01T09:00:00.000Z"));
        m.insert(
            "pharmcatVersion".to_owned(),
            s(&image_tag(ctx, "PHARMCAT_IMAGE")?),
        );
    }
    let hla_tests: HashMap<&str, &[&str]> = HashMap::from([
        ("HLA-A", &["*31:01"][..]),
        ("HLA-B", &["*15:02", "*57:01", "*58:01"][..]),
    ]);
    let mut changed: HashMap<String, String> = HashMap::new();
    for (gene, dip) in &outside {
        let alleles: Vec<&str> = dip.split('/').collect();
        let [a1, a2] = alleles[..] else {
            return Err(format!(
                "og pgx outside-calls: {gene} {dip} is not two alleles"
            ));
        };
        let phen: Vec<Value> = match hla_tests.get(gene.as_str()) {
            Some(tests) => tests
                .iter()
                .map(|t| {
                    let sign = if alleles.contains(t) {
                        "positive"
                    } else {
                        "negative"
                    };
                    s(&format!("{t} {sign}"))
                })
                .collect(),
            None if *dip == cyp => vec![s(&cyp_phen)],
            None => vec![s("Indeterminate")],
        };
        if *dip != called(gene)?.0 {
            changed.insert(gene.clone(), dip.clone());
        }
        let entry = get_mut(&mut report, "genes")
            .and_then(|g| get_mut(g, gene))
            .ok_or_else(|| format!("the PharmCAT example has no gene {gene}"))?;
        if let Value::Object(m) = entry {
            m.insert("callSource".to_owned(), s("OUTSIDE"));
            let allele = |name: &str| {
                obj(vec![
                    ("gene", s(gene)),
                    ("name", s(name)),
                    ("function", Value::Null),
                ])
            };
            m.insert(
                "sourceDiplotypes".to_owned(),
                Value::Array(vec![obj(vec![
                    ("gene", s(gene)),
                    ("label", s(dip)),
                    ("phenotypes", Value::Array(phen)),
                    ("activityScore", Value::Null),
                    ("outsidePhenotype", Value::Bool(false)),
                    ("allele1", allele(a1)),
                    ("allele2", allele(a2)),
                ])]),
            );
        }
    }
    // MT-RNR1 is an outside call in PharmCAT's example; this pipeline passes
    // none for it, so it is read from the VCF here, as the reference allele.
    let mt = get_mut(&mut report, "genes")
        .and_then(|g| get_mut(g, "MT-RNR1"))
        .ok_or("the PharmCAT example has no MT-RNR1")?;
    if let Value::Object(m) = mt {
        m.insert("callSource".to_owned(), s("MATCHER"));
        m.insert(
            "sourceDiplotypes".to_owned(),
            Value::Array(vec![obj(vec![
                ("gene", s("MT-RNR1")),
                ("label", s("Reference")),
                (
                    "phenotypes",
                    Value::Array(vec![s(
                        "normal risk of aminoglycoside-induced hearing loss",
                    )]),
                ),
                ("activityScore", Value::Null),
                ("outsidePhenotype", Value::Bool(false)),
                (
                    "allele1",
                    obj(vec![
                        ("gene", s("MT-RNR1")),
                        ("name", s("Reference")),
                        ("function", Value::Null),
                    ]),
                ),
                ("allele2", Value::Null),
            ])]),
        );
    }
    changed.insert("MT-RNR1".to_owned(), "Reference".to_owned());
    prune_annotations(&mut report, &changed);
    let text = json::dumps(
        &report,
        &Dump {
            indent: Some(1),
            ..Dump::default()
        },
    )
    .map_err(|e| e.to_string())?;
    let report_path = ctx.path(&format!("pharmcat/{sample}.report.json"));
    put(&report_path, &(text + "\n"))?;
    // Step 27's tables from that report, by the parser step 27 runs.
    let cpic = ctx.path("cpic");
    ctx.og(&[
        "pgx",
        "cpic-report",
        "--sample",
        sample,
        "--report",
        arg(&report_path),
        "--outdir",
        arg(&cpic),
    ])
}

/// Drops the example's drug annotations written for a diplotype of a gene
/// whose call changed: PharmCAT lists only those that fit the sample. A
/// guideline left with no annotation goes, then a drug with no guideline,
/// then a source with no drug.
pub fn prune_annotations(report: &mut Value, changed: &HashMap<String, String>) {
    let Some(Value::Object(drugs)) = get_mut(report, "drugs") else {
        return;
    };
    for by_drug in drugs.values_mut() {
        let Value::Object(by_drug) = by_drug else {
            continue;
        };
        for rep in by_drug.values_mut() {
            let Value::Object(rep) = rep else {
                continue;
            };
            let mut guidelines = match rep.get("guidelines") {
                Some(Value::Array(g)) => g.clone(),
                _ => Vec::new(),
            };
            for gl in &mut guidelines {
                let anns = match get(gl, "annotations") {
                    Some(Value::Array(a)) => a.clone(),
                    _ => Vec::new(),
                };
                let keep: Vec<Value> = anns.into_iter().filter(|a| !stale(a, changed)).collect();
                if let Value::Object(m) = gl {
                    m.insert("annotations".to_owned(), Value::Array(keep));
                }
            }
            guidelines.retain(
                |gl| matches!(get(gl, "annotations"), Some(Value::Array(a)) if !a.is_empty()),
            );
            rep.insert("guidelines".to_owned(), Value::Array(guidelines));
        }
        by_drug.retain(
            |_, rep| matches!(get(rep, "guidelines"), Some(Value::Array(g)) if !g.is_empty()),
        );
    }
    drugs.retain(|_, by_drug| !matches!(by_drug, Value::Object(m) if m.is_empty()));
}

/// Was this annotation written for another call of a changed gene: a
/// diplotype of it with another label, or a phenotype of it with no
/// diplotype of it?
fn stale(ann: &Value, changed: &HashMap<String, String>) -> bool {
    let mut labels: Vec<(Option<String>, Option<String>)> = Vec::new();
    if let Some(Value::Array(gts)) = get(ann, "genotypes") {
        for gt in gts {
            if let Some(Value::Array(dps)) = get(gt, "diplotypes") {
                for dp in dps {
                    if let Value::Object(m) = dp {
                        let text = |k: &str| match m.get(k) {
                            Some(Value::Str(x)) => Some(x.clone()),
                            _ => None,
                        };
                        labels.push((text("gene"), text("label")));
                    }
                }
            }
        }
    }
    let other_label = labels.iter().any(|(g, lab)| {
        g.as_ref()
            .and_then(|g| changed.get(g))
            .is_some_and(|want| lab.as_deref() != Some(want.as_str()))
    });
    let phenotype_genes: Vec<String> = match get(ann, "phenotypes") {
        Some(Value::Object(m)) => m.keys().cloned().collect(),
        Some(Value::Array(a)) => a
            .iter()
            .filter_map(|x| match x {
                Value::Str(t) => Some(t.clone()),
                _ => None,
            })
            .collect(),
        _ => Vec::new(),
    };
    other_label
        || phenotype_genes.iter().any(|g| {
            changed.contains_key(g) && !labels.iter().any(|(x, _)| x.as_deref() == Some(g.as_str()))
        })
}

/// Step 53: what `og locityper report` writes for a made-up Locityper run
/// that called every locus of `assets/locityper_loci.bed` with two HPRC
/// panel haplotypes, made-up qualities and read counts, and passed the
/// authors' filter (no warning, a weighted distance under 30, few
/// unexplained reads).
pub fn locityper(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let panel = [
        "HG00438.1",
        "HG00733.2",
        "HG01891.1",
        "HG02257.2",
        "NA18906.1",
        "NA19240.2",
    ];
    let mut out = String::from(
        "locus\tregion\tstatus\tgenotype\thaplotype_1\thaplotype_2\tquality\ttotal_reads\tunexplained_reads\tweighted_distance\twarnings\tfilter\treason\n",
    );
    let bed = read(&ctx.asset("assets/locityper_loci.bed"))?;
    for line in bed.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        let [chrom, start, end, name, ..] = c[..] else {
            return Err(format!("assets/locityper_loci.bed: short row {line:?}"));
        };
        let start: i64 = start
            .parse()
            .map_err(|_| format!("assets/locityper_loci.bed: start {start:?}"))?;
        let haps = rng.sample_idx(panel.len(), 2);
        let (h1, h2) = (panel[haps[0]], panel[haps[1]]);
        let quality = rng.uniform(25.0, 60.0);
        let total = rng.int(800, 4000);
        let unexplained = rng.int(0, 30);
        let dist = rng.uniform(0.1, 3.0);
        out.push_str(&format!(
            "{name}\t{chrom}:{}-{end}\tcalled\t{h1},{h2}\t{h1}\t{h2}\t{quality:.1}\t{total}\t{unexplained}\t{dist:.3}\t.\tpass\t.\n",
            start + 1
        ));
    }
    put(
        &ctx.path(&format!("locityper/{}_locityper.tsv", ctx.sample)),
        &out,
    )
}
