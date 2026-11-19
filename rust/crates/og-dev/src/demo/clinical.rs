//! The demo sample's clinical chapter, but for ClinVar
//! ([`super::clinical_clinvar`]) and the structural variants
//! ([`super::sv`]): the clinical filter, slivar and CPSR (steps 17, 23,
//! 31), repeat expansions (09, 09b), mitochondrial variants (20), Lp(a)
//! (46), CHIP (49), the Exomiser ranking (48), MoChA (51) and VNtyper (52).

use super::clinical_clinvar::PLANTED;
use super::{Ctx, Rng, autosomes, put, put_gz, read, table};

/// A gene with an invented name, so no real gene is tied to an invented person.
fn gene(rng: &mut Rng) -> String {
    format!("G{:04}", rng.int(1, 900))
}

/// Steps 17, 23 and 31: the clinical filter's table, slivar's tables and
/// CPSR's classification, with invented genes.
pub fn clinical(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let impacts = std::iter::repeat_n("HIGH", 41)
        .chain(std::iter::repeat_n("MODERATE", 233))
        .chain(std::iter::repeat_n("LOW", 44));
    let mut rows = String::from(
        "CHROM\tPOS\tREF\tALT\tGT\tIMPACT\tGENE\tConsequence\tMAX_AF\tCADD_PHRED\tREVEL\tAM_CLASS\n",
    );
    for imp in impacts {
        let (c, length) = autosomes()[rng.index(22)];
        let pos = rng.int(1, length);
        rows.push_str(&format!(
            "{c}\t{pos}\tC\tT\t0/1\t{imp}\t{}\tmissense_variant\t0.0004\t24.1\t0.41\t.\n",
            gene(rng)
        ));
    }
    put(
        &ctx.path(&format!("clinical/{s}_clinical_summary.tsv")),
        &rows,
    )?;

    let head = "CHROM\tPOS\tREF\tALT\tIMPACT\tGENE\n";
    let mut sl = String::from(head);
    for _ in 0..57 {
        let (c, p) = (rng.int(1, 22), rng.int(1, 90_000_000));
        let imp = *rng.pick(&["HIGH", "MODERATE"]);
        sl.push_str(&format!("chr{c}\t{p}\tA\tG\t{imp}\t{}\n", gene(rng)));
    }
    put(&ctx.path(&format!("slivar/{s}_slivar_summary.tsv")), &sl)?;
    let mut ch = String::from(head);
    for _ in 0..3 {
        let g = gene(rng);
        for _ in 0..2 {
            let (c, p) = (rng.int(1, 22), rng.int(1, 90_000_000));
            ch.push_str(&format!("chr{c}\t{p}\tA\tG\tMODERATE\t{g}\n"));
        }
    }
    put(&ctx.path(&format!("slivar/{s}_compound_hets.tsv")), &ch)?;

    let mut cls = String::from("SAMPLE_ID\tGENOMIC_CHANGE\tCPSR_CLASSIFICATION\tCLASSIFICATION\n");
    for (k, n) in [("Benign", 289), ("Likely_Benign", 112), ("VUS", 38)] {
        for _ in 0..n {
            let (c, p) = (rng.int(1, 22), rng.int(1, 90_000_000));
            cls.push_str(&format!("{s}\tchr{c}:g.{p}C>T\t{k}\t{k}\n"));
        }
    }
    put_gz(
        &ctx.path(&format!("cpsr/{s}.cpsr.grch38.classification.tsv.gz")),
        &cls,
    )?;
    put(
        &ctx.path(&format!("cpsr/{s}.cpsr.grch38.html")),
        &format!(
            "<!DOCTYPE html><html><head><title>CPSR {s}</title></head><body>Invented sample.</body></html>\n"
        ),
    )
}

/// ExpansionHunter's key loci with the normal range each repeat count is
/// drawn from: locus, chrom, pos, unit, low, high.
const EH_RANGES: [(&str, &str, i64, &str, i64, i64); 8] = [
    ("HTT", "chr4", 3074877, "CAG", 10, 26),
    ("FMR1", "chrX", 147912051, "CGG", 20, 40),
    ("C9ORF72", "chr9", 27573529, "GGCCCC", 2, 10),
    ("ATXN1", "chr6", 16327636, "TGC", 25, 35),
    ("DMPK", "chr19", 45770205, "CAG", 5, 30),
    ("ATXN3", "chr14", 92071011, "GCT", 14, 30),
    ("AR", "chrX", 67545317, "GCA", 17, 26),
    ("FXN", "chr9", 69037287, "GAA", 6, 30),
];

/// Stranger's thresholds per locus (normal up to, pathologic from),
/// invented for the demo's catalog.
fn str_thresholds(locus: &str) -> (i64, i64) {
    match locus {
        "HTT" => (26, 40),
        "FMR1" => (54, 200),
        "C9ORF72" => (20, 30),
        "ATXN1" => (35, 39),
        "DMPK" => (34, 50),
        "ATXN3" => (44, 60),
        "AR" => (34, 38),
        _ => (33, 66), // FXN
    }
}

/// Steps 09 and 09b: ExpansionHunter's VCF, and Stranger's copy of it with
/// the catalog's thresholds, the worst status of the alleles, and
/// ExpansionHunter's 95% interval of each (REPCI), a few repeats wide past
/// 50. DEMO-002 carries an FMR1 allele in the premutation range.
pub fn expansions(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let mut eh = String::new();
    let mut stranger = String::new();
    for (rep, c, pos, unit, lo, hi) in EH_RANGES {
        let hap = ctx.haploid(c);
        let mut a: Vec<i64> = (0..if hap { 1 } else { 2 })
            .map(|_| rng.int(lo, hi))
            .collect();
        a.sort_unstable();
        if rep == "FMR1" && ctx.has("fmr1") && !hap {
            a = vec![a[0], 62]; // above 54, under 200: Stranger's premutation range
        }
        let reference = (lo + hi) / 2;
        let gt = if hap {
            "1"
        } else if a[0] != a[1] {
            "1/2"
        } else {
            "1/1"
        };
        let cn = a.iter().map(i64::to_string).collect::<Vec<_>>().join("/");
        let mut distinct = a.clone();
        distinct.dedup();
        let alts = distinct
            .iter()
            .map(|x| format!("<STR{x}>"))
            .collect::<Vec<_>>()
            .join(",");
        let ulen = unit.len() as i64;
        let info = format!(
            "END={};REF={reference};RL={};RU={unit};VARID={rep};REPID={rep}",
            pos + reference * ulen,
            reference * ulen
        );
        let lc = format!("{:.6}", rng.uniform(28.0, 34.0));
        let first = &unit[..1];
        eh.push_str(&format!(
            "{c}\t{pos}\t.\t{first}\t{alts}\t.\tPASS\t{info}\tGT:SO:REPCN:LC\t{gt}:SPANNING:{cn}:{lc}\n"
        ));
        let (nmax, pmin) = str_thresholds(rep);
        let worst = if a.iter().any(|x| *x >= pmin) {
            "full_mutation"
        } else if a.iter().any(|x| *x > nmax) {
            "pre_mutation"
        } else {
            "normal"
        };
        let ci = a
            .iter()
            .map(|x| {
                if *x < 50 {
                    format!("{x}-{x}")
                } else {
                    format!("{}-{}", x - 4, x + 5)
                }
            })
            .collect::<Vec<_>>()
            .join("/");
        stranger.push_str(&format!(
            "{c}\t{pos}\t.\t{first}\t{alts}\t.\tPASS\t{info};STR_STATUS={worst};STR_NORMAL_MAX={nmax};STR_PATHOLOGIC_MIN={pmin}\tGT:SO:REPCN:REPCI:LC\t{gt}:SPANNING:{cn}:{ci}:{lc}\n"
        ));
    }
    let head = ctx.vcf_header(&["##ALT=<ID=STR,Description=\"Short tandem repeat\">"]);
    let s = &ctx.sample;
    put(
        &ctx.path(&format!("expansion_hunter/{s}_eh.vcf")),
        &format!("{head}{eh}"),
    )?;
    put(
        &ctx.path(&format!("expansion_hunter/{s}_eh_stranger.vcf")),
        &format!("{head}{stranger}"),
    )
}

/// Step 20's filtered chrM calls: 37 homoplasmic PASS records, two of them
/// heteroplasmic, and four weak ones.
pub fn mito(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let mut text = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{s}\n"
    );
    for (i, pos) in rng.distinct_sorted(60, 16500, 41).into_iter().enumerate() {
        let b = rng.index(4);
        let (r, a) = (super::BASES[b] as char, super::BASES[(b + 2) % 4] as char);
        let af = if i == 11 || i == 23 {
            rng.uniform(0.11, 0.38)
        } else if i < 37 {
            rng.uniform(0.985, 1.0)
        } else {
            rng.uniform(0.01, 0.03)
        };
        let filt = if i < 37 { "PASS" } else { "weak_evidence" };
        text.push_str(&format!(
            "chrM\t{pos}\t.\t{r}\t{a}\t.\t{filt}\t.\tGT:AF:DP\t0/1:{af:.3}:{}\n",
            rng.int(2400, 3300)
        ));
    }
    put_gz(&ctx.path(&format!("mito/{s}_chrM_filtered.vcf.gz")), &text)
}

/// The table `og cn call --locus LPA_KIV2` writes for the demo's step 46
/// depths: about 37 KIV-2 repeats at a flank depth near 31, half the KIV-2
/// depth at MAPQ 0.
pub const LPA_TABLE: &str = concat!(
    "metric\tvalue\n",
    "status\testimated\n",
    "kiv2_copies\t37.0\n",
    "ci95_low\t36.3\n",
    "ci95_high\t37.8\n",
    "class\ttypical\n",
    "kiv2_depth\t86.61\n",
    "flank_depth\t31.09\n",
    "mapq0_share\t0.500\n",
    "message\tabout 37 KIV-2 repeats in total (95% interval 36.3 to 37.8): within or overlapping the range most people have (30 to 45 for about 70% of the 1000 Genomes samples)\n",
    "advice\tA research estimate from read depth, not a diagnostic test. Lp(a) is measured in blood; a 2022 European Atherosclerosis Society consensus advises measuring it at least once in every adult.\n",
);

/// Step 46: [`LPA_TABLE`].
pub fn lpa(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    put(
        &ctx.path(&format!("lpa/{}_lpa_kiv2.tsv", ctx.sample)),
        LPA_TABLE,
    )
}

const CHIP_COLUMNS: [&str; 15] = [
    "gene",
    "chrom",
    "pos",
    "ref",
    "alt",
    "consequence",
    "protein_change",
    "exon",
    "vaf",
    "alt_reads",
    "depth",
    "filter",
    "status",
    "class",
    "note",
];

/// og chip's detection limit: the smallest VAF that gives `k` or more reads
/// with the change with probability `power` at `depth` (binomial, by
/// bisection).
pub fn chip_limit(depth: f64, k: u32, power: f64) -> f64 {
    let n = depth as u32;
    let p_at_least = |f: f64| {
        let mut below = 0.0;
        let mut comb = 1.0; // C(n, i)
        for i in 0..k {
            below += comb * f.powi(i as i32) * (1.0 - f).powi((n - i) as i32);
            comb = comb * f64::from(n - i) / f64::from(i + 1);
        }
        1.0 - below
    };
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if p_at_least(mid) >= power {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

/// Step 49: no CHIP call (DEMO-002: one small DNMT3A clone), and one depth
/// row per gene of `assets/chip_genes.tsv` at a made-up mean depth near 31,
/// with og chip's detection limit.
pub fn chip(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let genes = table(&read(&ctx.asset("assets/chip_genes.tsv"))?);
    let mut depth = String::from(
        "gene\ttranscript\tchrom\tcoding_bases\tmean_depth\tdetection_limit_vaf\tstatus\n",
    );
    for g in &genes {
        let mean = rng.uniform(29.0, 33.0);
        depth.push_str(&format!(
            "{}\t\t\t{}\t{mean:.1}\t{:.3}\tcovered\n",
            super::cell(g, "gene"),
            rng.int(1500, 7000),
            chip_limit(mean, 3, 0.8)
        ));
    }
    let mut calls = CHIP_COLUMNS.join("\t") + "\n";
    if ctx.has("chip") {
        // One small DNMT3A clone, above the detection limit at this depth.
        calls.push_str("DNMT3A\tchr2\t25234373\tC\tT\tmissense_variant\tp.R882H\t23\t0.161\t8\t50\tPASS\tchip\thotspot\ta published CHIP hotspot\n");
    }
    let s = &ctx.sample;
    put(&ctx.path(&format!("chip/{s}_chip.tsv")), &calls)?;
    put(&ctx.path(&format!("chip/{s}_chip_depth.tsv")), &depth)
}

/// One ranked gene: rank, gene, the four scores, variant, HGVS, disease.
type Ranked = (u32, String, [&'static str; 4], String, String, &'static str);

/// Step 48: what `og exomiser report` writes for made-up Exomiser rows for
/// two searched symptoms: two weak candidates of uncertain significance,
/// the shape of a search that finds nothing to act on. DEMO-003's search
/// also ranks the planted gene third, at the planted hit: a symptom search
/// can rank the gene of a result his page keeps closed.
pub fn exomiser(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let terms = [("HP:0001250", "Seizure"), ("HP:0002315", "Headache")];
    let mut ranks: Vec<Ranked> = vec![
        (
            1,
            "SCN1A".to_owned(),
            ["0.412", "0.788", "0.534", "0.021"],
            "2:165991230 C>T".to_owned(),
            "SCN1A:ENST00000674923.1:c.3199G>A:p.(Ala1067Thr)".to_owned(),
            "Dravet syndrome (OMIM:607208)",
        ),
        (
            2,
            "CACNA1A".to_owned(),
            ["0.297", "0.701", "0.421", "0.048"],
            "19:13299371 G>A".to_owned(),
            "CACNA1A:ENST00000360228.11:c.2134C>T:p.(Arg712Cys)".to_owned(),
            "Migraine, familial hemiplegic, 1 (OMIM:141500)",
        ),
    ];
    if ctx.has("act") {
        let (c, pos, r, a) = PLANTED.hit;
        ranks.push((
            3,
            PLANTED.gene.to_owned(),
            ["0.205", "0.512", "0.399", "0.061"],
            format!("{}:{pos} {r}>{a}", c.trim_start_matches("chr")),
            format!("{}:ENST00000262186.10:c.1C>T:p.(Met1?)", PLANTED.gene),
            "Long QT syndrome 2 (OMIM:613688)",
        ));
    }
    let freq = "GNOMAD_E_NFE 0.0123%";
    let mut top = String::from(
        "rank\tgene\tmoi\texomiser_score\tphenotype_score\tvariant_score\tp_value\tvariant\tgenotype\tconsequence\thgvs\tacmg\tacmg_evidence\tdisease\tclinvar\tmax_frequency\n",
    );
    let searched: Vec<String> = terms.iter().map(|(t, l)| format!("{l} ({t})")).collect();
    let mut words = vec![format!("Symptoms searched: {}.", searched.join(", "))];
    for (rank, g, sc, var, hgvs, disease) in &ranks {
        top.push_str(&format!(
            "{rank}\t{g}\tAD\t{}\t{var}\t0/1\tmissense_variant\t{hgvs}\tUNCERTAIN_SIGNIFICANCE\t.\t{disease}\t.\t{freq}\n",
            sc.join("\t")
        ));
        words.push(format!(
            "{rank}. {g} (autosomal dominant): Exomiser score {} (symptom match {}, variant {}). {var} heterozygous, missense_variant {hgvs}, ACMG uncertain significance, linked disease {disease}, highest population frequency {freq}.",
            sc[0], sc[1], sc[2]
        ));
    }
    words.push("A ranking to bring to a clinical geneticist, not a diagnosis. A high score means a rare variant, predicted to be damaging, in a gene linked to symptoms like the ones given; it does not show that the variant causes them. A low score or an empty list does not rule a genetic cause out: this search reads coding and splice variants against the terms given, and not structural variants, repeat expansions, mitochondrial heteroplasmy or variants outside genes.".to_owned());
    let s = &ctx.sample;
    let mut given = String::from("hpo_id\tlabel\tgiven_as\n");
    for (t, l) in terms {
        given.push_str(&format!("{t}\t{l}\t{t}\n"));
    }
    put(
        &ctx.path(&format!("exomiser/{s}_exomiser_terms.tsv")),
        &given,
    )?;
    put(&ctx.path(&format!("exomiser/{s}_exomiser_top.tsv")), &top)?;
    put(
        &ctx.path(&format!("exomiser/{s}_exomiser.txt")),
        &(words.join("\n") + "\n"),
    )
}

/// Step 51: what `og mocha report` writes for made-up MoChA tables of a 30x
/// genome (chr1 to chr22 read, 95,000 hets per chromosome of 100 Mb, the
/// sample filter met), with no mosaic call; DEMO-002's carries a loss of
/// 20q in about 12% of the blood cells, which MoChA's filters keep and the
/// floor lets through.
pub fn mocha(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let mut calls = String::from(
        "chrom\tstart\tend\tlength_mb\tarms\ttype\tcell_fraction\tbdev\trel_cov\tn_hets\tlod_baf_phase\tlod_lrr_baf\treported\treason\n",
    );
    if ctx.has("mocha") {
        calls.push_str("chr20\t30100000\t42300000\t12.20\t20q\tloss\t0.120\t0.0310\t1.8600\t4210\t48.2\t35.1\tyes\t.\n");
    }
    let n = if ctx.has("mocha") { "1" } else { "0" };
    let chroms: Vec<String> = (1..=22).map(|i| format!("chr{i}")).collect();
    let chroms = chroms.join(",");
    let summary = [
        ("status", "assessed"),
        ("message", "."),
        ("depth", "31.0"),
        ("het_sites", "2090000"),
        ("phased_het_sites", "2068000"),
        ("het_per_10mb", "9500"),
        ("analysed_mb", "2200.0"),
        ("chromosomes", chroms.as_str()),
        ("call_rate", "0.9981"),
        ("baf_auto", "0.0124"),
        ("floor_read", "0.011"),
        ("floor_min", "0.05"),
        ("floor", "0.050"),
        ("calls_mocha", n),
        ("calls_reported", n),
        ("phasing", "beagle"),
    ];
    let mut text = String::from("metric\tvalue\n");
    for (k, v) in summary {
        text.push_str(&format!("{k}\t{v}\n"));
    }
    let s = &ctx.sample;
    put(&ctx.path(&format!("mosaic/{s}_mocha.tsv")), &calls)?;
    put(&ctx.path(&format!("mosaic/{s}_mocha_summary.tsv")), &text)
}

/// Step 52: what `og vntyper report` writes for a made-up VNtyper run with
/// no MUC1 frameshift, at a VNTR depth of 212x that its coverage check
/// passes.
pub fn vntyper(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let rows = [
        ("status", "negative"),
        ("confidence", "Negative"),
        ("variant", "."),
        ("flag", "."),
        ("alt_depth", "."),
        ("region_depth", "."),
        ("depth_score", "."),
        ("vntr_mean_depth", "212.40"),
        ("vntr_uncovered_pct", "0.00"),
        ("coverage_qc", "PASS"),
        (
            "sensitivity",
            "Mean depth over the MUC1 VNTR: 212x. At that depth (at 100x or more) VNtyper found the known variant in 160 of 161 runs (99%), in VNtyper 2.0.47's downsampling of 42 confirmed positive GRCh38 exomes (repeated runs, not independent tests).",
        ),
        ("message", "VNtyper calls no MUC1 VNTR frameshift."),
        (
            "advice",
            "A negative here lowers the chance of the frameshifts VNtyper looks for; it does not rule ADTKD-MUC1 out. A family history of kidney failure with gout or a bland urine sediment still deserves a clinical test (SNaPshot for dupC).",
        ),
    ];
    let mut text = String::from("metric\tvalue\n");
    for (k, v) in rows {
        text.push_str(&format!("{k}\t{v}\n"));
    }
    put(
        &ctx.path(&format!("vntyper/{}_vntyper.tsv", ctx.sample)),
        &text,
    )
}
