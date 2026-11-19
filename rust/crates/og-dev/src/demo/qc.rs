//! The demo sample's quality chapter: the run status, the variant calls
//! (step 03), coverage (steps 16 and 16b), the sample checks and mosaic
//! rows (step 33), the reads (step 01b and samtools), the read count and
//! the genes the genome could not read (step 16b).

use super::{
    BASES, Ctx, MEAN_DEPTH, RUN_STARTED, Rng, arg, autosomes, cell, erfc, finish_gz, gz_writer,
    put, put_gz, read, table, utc,
};
use indexmap::IndexMap;
use pycompat::json::{self, Dump, Value};
use std::io::Write;

/// `logs/run_status.tsv` as run-all.sh writes it, every step it records ok
/// (the `recorded=yes` rows of `assets/steps.tsv`, in sorted order).
pub fn run_status(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let steps_text = read(&ctx.asset("assets/steps.tsv"))?;
    let mut steps: Vec<String> = table(&steps_text)
        .iter()
        .filter(|r| cell(r, "recorded") == "yes")
        .map(|r| cell(r, "step").to_owned())
        .collect();
    if steps.is_empty() {
        return Err("assets/steps.tsv has no recorded step".to_owned());
    }
    steps.sort();
    let mut text = format!(
        "# run-all.sh: when this run started and how each step ended\nmeta\tstarted_epoch\t{RUN_STARTED}\nmeta\tstarted_utc\t{}\nmeta\tdeclared_sex\t{}\n",
        utc(RUN_STARTED),
        ctx.sex()
    );
    for s in steps {
        text.push_str(&format!("step\t{s}\tok\n"));
    }
    put(&ctx.path("logs/run_status.tsv"), &text)
}

/// Variants per base of a sex chromosome, against an autosome (one X, one Y).
fn density(chrom: &str) -> f64 {
    match chrom {
        "chrX" => 0.45,
        "chrY" => 0.04,
        _ => 1.0,
    }
}

/// `vcf/S.vcf.gz`: `--records` records at random positions, about 94% PASS,
/// 85% SNVs.
pub fn variants(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let chroms = ctx.own_chroms();
    let weights: Vec<f64> = chroms.iter().map(|(c, n)| *n as f64 * density(c)).collect();
    let total: f64 = weights.iter().sum();
    let path = ctx.path(&format!("vcf/{}.vcf.gz", ctx.sample));
    let mut gz = gz_writer(&path, 1)?;
    let head = ctx.vcf_header(&[
        "##FILTER=<ID=RefCall,Description=\"Genotyping model thinks this site is reference.\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Conditional genotype quality\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
    ]);
    let io = |e: std::io::Error| format!("{}: {e}", path.display());
    gz.write_all(head.as_bytes()).map_err(io)?;
    let indel = |rng: &mut Rng| -> String {
        (0..1 + rng.below(6))
            .map(|_| BASES[rng.index(4)] as char)
            .collect()
    };
    for ((chrom, length), w) in chroms.iter().zip(&weights) {
        let k = (ctx.records as f64 * w / total).round() as usize;
        let hap = ctx.haploid(chrom);
        let mut out = String::with_capacity(k * 48);
        for pos in rng.distinct_sorted(10001, length - 10000, k) {
            let b = rng.index(4);
            let mut reference = (BASES[b] as char).to_string();
            let r = rng.random();
            let alt = if r < 0.855 {
                (BASES[(b + 1 + rng.index(3)) % 4] as char).to_string()
            } else if r < 0.93 {
                format!("{reference}{}", indel(rng))
            } else {
                let alt = reference.clone();
                reference.push_str(&indel(rng));
                alt
            };
            let dp = (rng.gauss(MEAN_DEPTH / if hap { 2.0 } else { 1.0 }, 6.0) as i64).max(3);
            if rng.random() < 0.942 {
                let gt = if hap || rng.random() < 0.38 {
                    "1/1"
                } else {
                    "0/1"
                };
                out.push_str(&format!(
                    "{chrom}\t{pos}\t.\t{reference}\t{alt}\t{}.{}\tPASS\t.\tGT:GQ:DP\t{gt}:{}:{dp}\n",
                    rng.int(20, 66),
                    rng.below(10),
                    rng.int(20, 60)
                ));
            } else {
                out.push_str(&format!(
                    "{chrom}\t{pos}\t.\t{reference}\t{alt}\t{}.{}\tRefCall\t.\tGT:GQ:DP\t0/0:{}:{dp}\n",
                    rng.int(0, 9),
                    rng.below(10),
                    rng.int(5, 30)
                ));
            }
        }
        gz.write_all(out.as_bytes()).map_err(io)?;
    }
    finish_gz(gz, &path)
}

/// The share of bases with at least each depth: a normal curve, plus the
/// gaps of the reference (no reads) at depth 0.
fn cumulative(mean: f64, sd: f64, gap: f64) -> Vec<(i64, f64)> {
    let mut out = Vec::new();
    for dep in 0..=(mean * 2.6) as i64 {
        let p = if dep == 0 {
            1.0
        } else {
            (1.0 - gap) * 0.5 * erfc((dep as f64 - 0.5 - mean) / (sd * std::f64::consts::SQRT_2))
        };
        if dep > 0 && p < 0.005 {
            break;
        }
        out.push((dep, p));
    }
    out
}

/// mosdepth (step 16b) and indexcov (step 16).
pub fn coverage(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let mut summary = vec!["chrom\tlength\tbases\tmean\tmin\tmax".to_owned()];
    let (mut tot_len, mut tot_bases) = (0i64, 0i64);
    let mut curves = Vec::new();
    for (chrom, length) in ctx.own_chroms() {
        let mean =
            MEAN_DEPTH / if ctx.haploid(chrom) { 2.0 } else { 1.0 } * rng.uniform(0.97, 1.03);
        let gap = if chrom == "chrY" {
            0.55
        } else {
            rng.uniform(0.01, 0.03)
        };
        let bases = (length as f64 * mean * (1.0 - gap)) as i64;
        summary.push(format!(
            "{chrom}\t{length}\t{bases}\t{:.2}\t0\t{}",
            bases as f64 / length as f64,
            (mean * 9.0) as i64
        ));
        tot_len += length;
        tot_bases += bases;
        curves.push((chrom, cumulative(mean, mean * 0.24, gap)));
    }
    summary.push(format!(
        "total\t{tot_len}\t{tot_bases}\t{:.2}\t0\t{}",
        tot_bases as f64 / tot_len as f64,
        (MEAN_DEPTH * 9.0) as i64
    ));
    let mut dist = Vec::new();
    for (chrom, curve) in &curves {
        dist.extend(curve.iter().map(|(d, p)| format!("{chrom}\t{d}\t{p:.2}")));
    }
    dist.extend(
        cumulative(MEAN_DEPTH, MEAN_DEPTH * 0.25, 0.02)
            .iter()
            .map(|(d, p)| format!("total\t{d}\t{p:.2}")),
    );
    put(
        &ctx.path(&format!("mosdepth/{s}.mosdepth.summary.txt")),
        &(summary.join("\n") + "\n"),
    )?;
    put(
        &ctx.path(&format!("mosdepth/{s}.mosdepth.global.dist.txt")),
        &(dist.join("\n") + "\n"),
    )?;

    let (sex, cnx, cny) = if ctx.male {
        (1, "1.01", "0.98")
    } else {
        (2, "1.99", "0.01")
    };
    put(
        &ctx.path("indexcov/indexcov-indexcov.ped"),
        &format!(
            "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\tbins.out\tbins.lo\tbins.hi\tbins.in\tslope\tp.out\tPC1\tPC2\tPC3\tPC4\tPC5\nunknown\t{s}\t-9\t-9\t{sex}\t-9\t{cnx}\t{cny}\t312\t205\t107\t174820\t0.993\t0.0018\t0\t0\t0\t0\t0\n"
        ),
    )?;
    let mut roc = vec![format!("#chrom\tcov\t{s}")];
    for (chrom, _) in ctx.own_chroms() {
        let centre = if ctx.haploid(chrom) { 0.5 } else { 1.0 };
        for k in 0..41 {
            let cov = f64::from(k) * 0.05;
            let frac = 0.5 * erfc((cov - centre) / (0.12 * std::f64::consts::SQRT_2));
            roc.push(format!("{chrom}\t{cov:.2}\t{:.3}", frac.min(1.0)));
        }
    }
    put(
        &ctx.path("indexcov/indexcov-indexcov.roc"),
        &(roc.join("\n") + "\n"),
    )
}

/// somalier and VerifyBamID2 files, then step 33's own table writer,
/// `og qc sample`. `sex` is the sex somalier infers from the reads;
/// `original_pedigree_sex` stays -9 (unknown), as in a real run: step 33
/// gives somalier no pedigree file.
pub fn sample_qc(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let sex = if ctx.male { 1 } else { 2 };
    let x = if ctx.male {
        "15.6\t612\t352\t4\t256\t14.9\t41"
    } else {
        "31.2\t612\t300\t210\t102\t0.1\t0"
    };
    let samples = ctx.path(&format!("qc/somalier/{s}.samples.tsv"));
    put(
        &samples,
        &format!(
            "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\toriginal_pedigree_sex\tgt_depth_mean\tgt_depth_sd\tdepth_mean\tdepth_sd\tab_mean\tab_std\tn_hom_ref\tn_het\tn_hom_alt\tn_unknown\tp_middling_ab\tX_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\tY_depth_mean\tY_n\n{s}\t{s}\t-9\t-9\t{sex}\t-9\t-9\t31.4\t7.2\t31.4\t7.2\t0.41\t0.08\t9817\t6904\t4383\t196\t0.01\t{x}\n"
        ),
    )?;
    let pairs = ctx.path(&format!("qc/somalier/{s}.pairs.tsv"));
    put(
        &pairs,
        "#sample_a\tsample_b\trelatedness\tibs0\tibs2\thom_concordance\thets_a\thets_b\thets_ab\tshared_hets\thom_alts_a\thom_alts_b\tshared_hom_alts\tn\tx_ibs0\tx_ibs2\texpected_relatedness\n",
    )?;
    let selfsm = ctx.path(&format!("qc/verifybamid2/{s}.selfSM"));
    put(
        &selfsm,
        &format!(
            "#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\tFREELK1\tFREELK0\tFREE_RH\tFREE_RA\tCHIPMIX\tCHIPLK1\tCHIPLK0\tCHIP_RH\tCHIP_RA\tDPREF\tRDPHET\tRDPALT\n{s}\tNA\tNA\t99846\t3087112\t30.92\t0.00412\t-1418823.2\t-1419102.7\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
        ),
    )?;
    let out = ctx.path(&format!("qc/{s}_sample_qc.tsv"));
    ctx.og(&[
        "qc",
        "sample",
        "--sample",
        s,
        "--somalier-samples",
        arg(&samples),
        "--somalier-pairs",
        arg(&pairs),
        "--selfsm",
        arg(&selfsm),
        "--declared-sex",
        ctx.sex(),
        "--marker-check",
        "passed",
        "--out",
        arg(&out),
    ])
}

fn float(x: f64) -> Value {
    Value::Float(x)
}

fn int(x: i64) -> Value {
    Value::Int(i128::from(x))
}

fn obj(pairs: Vec<(&str, Value)>) -> Value {
    Value::Object(
        pairs
            .into_iter()
            .map(|(k, v)| (k.to_owned(), v))
            .collect::<IndexMap<_, _>>(),
    )
}

/// Python's `round(x, n)`.
fn round(x: f64, n: i32) -> f64 {
    pycompat::fmt::round(x, n).unwrap_or(x)
}

/// fastp's JSON (step 01b) and samtools flagstat (step 28 writes it from the BAM).
pub fn reads(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let pairs: i64 = 412_906_118;
    let before = pairs * 2;
    let share = |x: i64, f: f64| (x as f64 * f) as i64;
    let passed = share(before, 0.9891);
    let low_q = share(before, 0.0081);
    let too_short = share(before, 0.0026);
    let n_reads = before - passed - low_q - too_short;
    let length: i64 = 151;
    let curves = |after: bool, rng: &mut Rng| -> Value {
        let q: Vec<Value> = (0..length)
            .map(|i| {
                let i = i as f64;
                let top = (31.5 + 6.0 * (1.0 - (-i / 4.0).exp())).min(37.4);
                let fall = (i / length as f64).powi(3) * if after { 2.2 } else { 3.5 };
                float(round(top - fall + rng.uniform(-0.15, 0.15), 2))
            })
            .collect();
        let mut base: Vec<(&str, Vec<f64>)> = Vec::new();
        for (b, p) in [("A", 0.295), ("T", 0.296), ("C", 0.204), ("G", 0.205)] {
            base.push((
                b,
                (0..length)
                    .map(|_| round(p + rng.uniform(-0.004, 0.004), 4))
                    .collect(),
            ));
        }
        let gc: Vec<Value> = (0..length as usize)
            .map(|i| float(round(base[2].1[i] + base[3].1[i], 4)))
            .collect();
        let mut content: Vec<(&str, Value)> = base
            .iter()
            .map(|(b, v)| (*b, Value::Array(v.iter().map(|x| float(*x)).collect())))
            .collect();
        content.push(("N", Value::Array(vec![float(0.0002); length as usize])));
        content.push(("GC", Value::Array(gc)));
        let n = if after { passed } else { before } / 2;
        let mut quality: Vec<(&str, Value)> = ["A", "T", "C", "G"]
            .iter()
            .map(|b| (*b, Value::Array(q.clone())))
            .collect();
        quality.push(("mean", Value::Array(q)));
        obj(vec![
            ("total_reads", int(n)),
            ("total_bases", int(n * length)),
            ("q20_bases", int((n as f64 * length as f64 * 0.968) as i64)),
            ("q30_bases", int((n as f64 * length as f64 * 0.921) as i64)),
            ("total_cycles", int(length)),
            ("quality_curves", obj(quality)),
            ("content_curves", obj(content)),
        ])
    };
    let stats = |n: i64, gc: f64| {
        obj(vec![
            ("total_reads", int(n)),
            ("total_bases", int(n * length)),
            ("q20_bases", int((n as f64 * length as f64 * 0.968) as i64)),
            ("q30_bases", int((n as f64 * length as f64 * 0.921) as i64)),
            ("q20_rate", float(0.968)),
            ("q30_rate", float(0.921)),
            ("read1_mean_length", int(length)),
            ("read2_mean_length", int(length)),
            ("gc_content", float(gc)),
        ])
    };
    let hist: Vec<Value> = (0..1000)
        .map(|i| {
            let z = (f64::from(i) - 382.0) / 118.0;
            int((4.2e6 * (-z * z / 2.0).exp()) as i64)
        })
        .collect();
    let part = format!("/genome/{s}/fastq_trimmed.part/{s}");
    let fastp = obj(vec![
        (
            "summary",
            obj(vec![
                ("fastp_version", Value::Str("1.0.1".to_owned())),
                (
                    "sequencing",
                    Value::Str("paired end (151 cycles + 151 cycles)".to_owned()),
                ),
                ("before_filtering", stats(before, 0.4105)),
                ("after_filtering", stats(passed, 0.4098)),
            ]),
        ),
        (
            "filtering_result",
            obj(vec![
                ("passed_filter_reads", int(passed)),
                ("low_quality_reads", int(low_q)),
                ("too_many_N_reads", int(n_reads)),
                ("too_short_reads", int(too_short)),
                ("too_long_reads", int(0)),
            ]),
        ),
        ("duplication", obj(vec![("rate", float(0.0712))])),
        (
            "insert_size",
            obj(vec![
                ("peak", int(382)),
                ("unknown", int(5_102_117)),
                ("histogram", Value::Array(hist)),
            ]),
        ),
        (
            "adapter_cutting",
            obj(vec![
                ("adapter_trimmed_reads", int(share(before, 0.074))),
                ("adapter_trimmed_bases", int(share(before, 1.9))),
                (
                    "read1_adapter_sequence",
                    Value::Str("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_owned()),
                ),
                (
                    "read2_adapter_sequence",
                    Value::Str("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_owned()),
                ),
            ]),
        ),
        ("read1_before_filtering", curves(false, rng)),
        ("read2_before_filtering", curves(false, rng)),
        ("read1_after_filtering", curves(true, rng)),
        ("read2_after_filtering", curves(true, rng)),
        (
            "command",
            Value::Str(format!(
                "fastp -i /genome/{s}/fastq/{s}_R1.fastq.gz -I /genome/{s}/fastq/{s}_R2.fastq.gz -o {part}_R1.fastq.gz -O {part}_R2.fastq.gz -j {part}_fastp.json -h {part}_fastp.html"
            )),
        ),
    ]);
    let text = json::dumps(
        &fastp,
        &Dump {
            indent: Some(1),
            ..Dump::default()
        },
    )
    .map_err(|e| e.to_string())?;
    put(
        &ctx.path(&format!("fastq_trimmed/{s}_fastp.json")),
        &(text + "\n"),
    )?;

    let primary = passed;
    let supp = share(primary, 0.0034);
    let total = primary + supp;
    let mapped = share(primary, 0.9962) + supp;
    let dups = share(primary, 0.0712);
    let paired = primary;
    let proper = share(paired, 0.981);
    let single = share(paired, 0.0019);
    let both = share(paired, 0.9943);
    let diff = share(paired, 0.0061);
    let pct = |a: i64, b: i64| 100.0 * a as f64 / b as f64;
    let flag = [
        format!("{total} + 0 in total (QC-passed reads + QC-failed reads)"),
        format!("{primary} + 0 primary"),
        "0 + 0 secondary".to_owned(),
        format!("{supp} + 0 supplementary"),
        format!("{dups} + 0 duplicates"),
        format!("{dups} + 0 primary duplicates"),
        format!("{mapped} + 0 mapped ({:.2}% : N/A)", pct(mapped, total)),
        format!(
            "{} + 0 primary mapped ({:.2}% : N/A)",
            mapped - supp,
            pct(mapped - supp, primary)
        ),
        format!("{paired} + 0 paired in sequencing"),
        format!("{} + 0 read1", paired / 2),
        format!("{} + 0 read2", paired / 2),
        format!(
            "{proper} + 0 properly paired ({:.2}% : N/A)",
            pct(proper, paired)
        ),
        format!("{both} + 0 with itself and mate mapped"),
        format!(
            "{single} + 0 singletons ({:.2}% : N/A)",
            pct(single, paired)
        ),
        format!("{diff} + 0 with mate mapped to a different chr"),
        format!(
            "{} + 0 with mate mapped to a different chr (mapQ>=5)",
            share(diff, 0.62)
        ),
    ];
    put(
        &ctx.path(&format!("aligned/{s}_flagstat.txt")),
        &(flag.join("\n") + "\n"),
    )
}

/// `qc/S_read_count.tsv` as `read_count_check` (scripts/lib/common.sh)
/// writes it: the BAM holds every read fastp kept, so the check is ok.
pub fn read_count(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let path = ctx.path(&format!("fastq_trimmed/{s}_fastp.json"));
    let doc = json::loads(&read(&path)?).map_err(|e| format!("{}: {e}", path.display()))?;
    let kept = match &doc {
        Value::Object(m) => match m.get("filtering_result") {
            Some(Value::Object(f)) => match f.get("passed_filter_reads") {
                Some(Value::Int(n)) => *n,
                _ => 0,
            },
            _ => 0,
        },
        _ => 0,
    };
    if kept == 0 {
        return Err(format!("{} has no passed_filter_reads", path.display()));
    }
    let half = kept / 2;
    put(
        &ctx.path(&format!("qc/{s}_read_count.tsv")),
        &format!(
            "key\tvalue\ninput_reads\t{}\ninput_source\tfastq_trimmed/{s}_R1.fastq.gz {half} records + fastq_trimmed/{s}_R2.fastq.gz {half} records\nprimary_records\t{}\nmissing_reads\t0\nread_count_check\tok\n",
            half * 2,
            half * 2
        ),
    )
}

/// Genes the demo leaves below step 16b's share: four that another step
/// owns and four with no other reading, so both lines of the card have
/// names in them.
const DEMO_NOT_READ: [&str; 8] = [
    "SMN1", "SMN2", "NCF1", "PMS2", "NBPF1", "ANKRD36C", "DUX4", "MUC12",
];

/// Step 16b's gene table (mosdepth `--thresholds 10` over the genes of
/// `assets/callability_genes.tsv`: chrom, start, end, gene, bases at 10x),
/// written from the repository's catalog, so it always matches it. The
/// ALT-only genes (start ".") get no row, as in the real step.
pub fn callability(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let text = read(&ctx.asset("assets/callability_genes.tsv"))?;
    let mut rows = vec!["#chrom\tstart\tend\tregion\t10X".to_owned()];
    for r in table(&text) {
        if cell(&r, "start") == "." {
            continue;
        }
        let num = |k: &str| {
            cell(&r, k)
                .parse::<i64>()
                .map_err(|_| format!("assets/callability_genes.tsv: {k} '{}'", cell(&r, k)))
        };
        let (start, end) = (num("start")?, num("end")?);
        let gene = cell(&r, "gene");
        let share = if DEMO_NOT_READ.contains(&gene) {
            rng.uniform(0.35, 0.8)
        } else {
            rng.uniform(0.93, 1.0)
        };
        rows.push(format!(
            "{}\t{start}\t{end}\t{gene}\t{}",
            cell(&r, "chrom"),
            ((end - start) as f64 * share) as i64
        ));
    }
    put_gz(
        &ctx.path(&format!("mosdepth/{}_callability.bed.gz", ctx.sample)),
        &(rows.join("\n") + "\n"),
    )
}

/// The first ten single-copy pieces of chrY and chrX that `og qc sex-loss`
/// reads (og-qc `SINGLE_COPY_Y` and `SINGLE_COPY_X`).
const SINGLE_COPY_Y: [(i64, i64); 10] = [
    (2781479, 3050044),
    (6748713, 7574481),
    (11037032, 11592902),
    (11674123, 13984473),
    (14058733, 15875093),
    (15905214, 16159795),
    (16425965, 17456266),
    (18870334, 20207793),
    (20257793, 21305953),
    (26415435, 26673214),
];
const SINGLE_COPY_X: [(i64, i64); 10] = [
    (2781479, 37099262),
    (37285837, 48343310),
    (48434604, 49157254),
    (49528394, 50228964),
    (50278964, 52032464),
    (52223402, 52489200),
    (52520537, 55437684),
    (55547739, 58100000),
    (61000000, 62412542),
    (62462542, 63185573),
];

/// mosdepth's 500 bp windows (step 16b's `regions.bed.gz`), cut down to
/// three 200 kb stretches of each autosome and the single-copy part of chrY
/// (a female's chrX), then `og qc sex-loss` adds its rows to step 33's
/// table. One Y in every cell: the ratio sits near 0.5 and no loss is seen.
pub fn mosaic(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let path = ctx.path(&format!("mosdepth/{}.regions.bed.gz", ctx.sample));
    let mut lines = Vec::new();
    for (chrom, length) in autosomes() {
        let mean = MEAN_DEPTH * rng.uniform(0.98, 1.02);
        for k in 0..3 {
            let start = (length / 4) * (k + 1) / 500 * 500;
            for s in (start..start + 200_000).step_by(500) {
                lines.push(format!(
                    "{chrom}\t{s}\t{}\t{:.2}",
                    s + 500,
                    rng.gauss(mean, mean * 0.12).max(0.0)
                ));
            }
        }
    }
    let (sx, pieces) = if ctx.male {
        ("Y", &SINGLE_COPY_Y)
    } else {
        ("X", &SINGLE_COPY_X)
    };
    let ymean = MEAN_DEPTH / if ctx.male { 2.0 } else { 1.0 } * rng.uniform(0.98, 1.0);
    for (start, end) in pieces {
        let s0 = (start / 500 + 1) * 500;
        for s in (s0..(*end).min(s0 + 1_000_000) - 500).step_by(500) {
            lines.push(format!(
                "chr{sx}\t{s}\t{}\t{:.2}",
                s + 500,
                rng.gauss(ymean, ymean * 0.14).max(0.0)
            ));
        }
    }
    put_gz(&path, &(lines.join("\n") + "\n"))?;
    let qc_table = ctx.path(&format!("qc/{}_sample_qc.tsv", ctx.sample));
    ctx.og(&[
        "qc",
        "sex-loss",
        "--table",
        arg(&qc_table),
        "--regions",
        arg(&path),
    ])
}
