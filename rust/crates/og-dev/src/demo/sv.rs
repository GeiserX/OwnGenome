//! The demo sample's structural variants and copy numbers: Manta, Delly,
//! SURVIVOR and CNVpytor (steps 04, 18, 19, 22), AnnotSV's class 4 and 5
//! table (05), mobile element insertions (50), and the paralog genes:
//! Paraphase (35b), Parascopy's SMN1 and the HBA1/HBA2 copy number (35).

use super::clinical_clinvar::PLANTED;
use super::{CHROMS, Ctx, MEAN_DEPTH, Rng, Scratch, arg, chrom_index, put, put_gz, read};

/// One SV record: chrom, pos, id, type, end, svlen, filter, genotype, qual.
type Sv = (
    &'static str,
    i64,
    String,
    &'static str,
    i64,
    i64,
    &'static str,
    &'static str,
    i64,
);

fn sv_records(rng: &mut Rng, n: i64, caller: &str, pass_rate: f64) -> Vec<Sv> {
    let types = ["DEL", "DUP", "INV", "INS"];
    let mut out = Vec::new();
    for i in 0..n {
        let (c, length) = CHROMS[rng.index(CHROMS.len() - 1)];
        let t = types[rng.weighted(&[62.0, 14.0, 6.0, 18.0])];
        let pos = rng.int(20_000, length - 2_000_000);
        let size = rng.uniform(60f64.ln(), 400_000f64.ln()).exp() as i64;
        let end = pos + if t == "INS" { 1 } else { size };
        let svlen = if t == "DEL" { -size } else { size };
        let filt = if rng.random() < pass_rate {
            "PASS"
        } else {
            *rng.pick(&["MinQUAL", "MinGQ", "LowQual"])
        };
        let gt = if rng.random() < 0.2 { "1/1" } else { "0/1" };
        out.push((
            c,
            pos,
            format!("{caller}{t}:{i}"),
            t,
            end,
            svlen,
            filt,
            gt,
            rng.int(20, 999),
        ));
    }
    out.sort_by_key(|r| (chrom_index(r.0), r.1));
    out
}

/// Manta's, Delly's and SURVIVOR's VCFs, CNVpytor's calls and step 05's
/// checked AnnotSV table: none of class 4 or 5, as in most genomes, but
/// DEMO-003's class 4 deletion over an invented gene and his planted class
/// 5 deletion over the planted ClinVar hit.
pub fn structural(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let head = ctx.vcf_header(&[
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of the SV\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    ]);
    let write = |rel: String, recs: &[Sv], extra: &[String]| {
        let mut body = head.clone();
        for (k, (c, p, id, t, e, l, f, g, q)) in recs.iter().enumerate() {
            let x = extra.get(k).map_or("", String::as_str);
            body.push_str(&format!(
                "{c}\t{p}\t{id}\tN\t<{t}>\t{q}\t{f}\tEND={e};SVTYPE={t};SVLEN={l}{x}\tGT\t{g}\n"
            ));
        }
        put_gz(&ctx.path(&rel), &body)
    };
    let n = rng.int(7600, 8000);
    write(
        "manta/results/variants/diploidSV.vcf.gz".to_owned(),
        &sv_records(rng, n, "Manta", 0.82),
        &[],
    )?;
    let n = rng.int(5900, 6300);
    write(
        format!("delly/{s}_sv.vcf.gz"),
        &sv_records(rng, n, "Delly", 0.86),
        &[],
    )?;
    // SURVIVOR's merge of Manta, Delly and CNVpytor (step 22): SUPP and SUPP_VEC.
    let n = rng.int(2700, 3100);
    let merged = sv_records(rng, n, "SURVIVOR", 1.0);
    let vecs = ["110", "101", "011", "111"];
    let extra: Vec<String> = merged
        .iter()
        .map(|_| {
            let v = vecs[rng.weighted(&[70.0, 9.0, 6.0, 15.0])];
            format!(";SUPP={};SUPP_VEC={v}", v.matches('1').count())
        })
        .collect();
    write(
        format!("sv_merged/{s}_sv_consensus.vcf.gz"),
        &merged,
        &extra,
    )?;

    let mut lines = Vec::new();
    for _ in 0..rng.int(1100, 1250) {
        let (c, length) = CHROMS[rng.index(CHROMS.len() - 1)];
        let t = if rng.random() < 0.83 {
            "deletion"
        } else {
            "duplication"
        };
        let size = 1000 + 100 * rng.below(3990) as i64;
        let start = 10001 + 100 * rng.below(((length - size - 10001) as u64).div_ceil(100)) as i64;
        let cn = if t == "deletion" {
            rng.uniform(0.02, 0.6)
        } else {
            rng.uniform(1.3, 2.1)
        };
        let p = 10f64.powf(rng.uniform(-40.0, -0.5));
        let sci = |x: f64| pycompat::fmt::sci(x, 4);
        lines.push(format!(
            "{t}\t{c}:{start}-{}\t{size}\t{cn:.4}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}",
            start + size - 1,
            sci(p),
            sci(p * 3.0),
            sci((p * 5.0).min(1.0)),
            sci((p * 9.0).min(1.0)),
            rng.uniform(0.0, 0.5),
            rng.uniform(0.0, 0.1),
            rng.int(0, 9)
        ));
    }
    put(
        &ctx.path(&format!("cnvpytor/{s}_cnvs.txt")),
        &(lines.join("\n") + "\n"),
    )?;

    // Step 05's class 4 and 5 SVs, checked (the columns og sv postfilter writes).
    let mut text = String::from(
        "AnnotSV_ID\tID\tSV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tGene_name\tACMG_class\tAnnotSV_ranking_score\tCaller\tExon_edge\tExon_edge_detail\tPopulation_AF\tPopulation_source\tCommon\tAnnotSV_class\tDepth_support\tRule\tRule_detail\n",
    );
    let row = |id: String,
               mid: &str,
               c: &str,
               a: i64,
               b: i64,
               gene: &str,
               class: &str,
               score: &str| {
        format!(
            "{id}\t{mid}\t{c}\t{a}\t{b}\t{}\tDEL\t{gene}\t{class}\t{score}\tManta\tno\t.\t\t\tno\t{class}\t.\t\t\n",
            a - b
        )
    };
    if ctx.has("annotsv") {
        // One class 4 deletion over an invented gene, from Manta, kept by step 05's rules.
        text.push_str(&row(
            "chr16_2084000_2096000_DEL_1".to_owned(),
            "MantaDEL:1",
            "chr16",
            2_084_000,
            2_096_000,
            "G0533",
            "4",
            "0.92",
        ));
    }
    if ctx.has("act") {
        // The class 5 deletion over the planted ClinVar hit: a result to act on.
        let (c, a, b) = PLANTED.sv;
        text.push_str(&row(
            format!("{c}_{a}_{b}_DEL_1"),
            "MantaDEL:2",
            c,
            a,
            b,
            PLANTED.gene,
            "5",
            "0.99",
        ));
    }
    put(
        &ctx.path(&format!("annotsv/{s}_sv_class4_5_checked.tsv")),
        &text,
    )
}

/// Step 35's HBA1/HBA2 depths: mean depth per region of the HBA row of
/// `bin/cn_loci.tsv`, every read, and the share of each left at MAPQ 1 or
/// more. og-cn's `port_tests_cn_traits` reads these three literals from
/// this file and checks that `og cn call --locus HBA` writes [`HBA_TABLE`]
/// for them.
pub const HBA_DEPTH: [(&str, f64); 4] = [
    ("flank_left", 31.2),
    ("HBA2", 30.9),
    ("HBA1", 30.1),
    ("flank_right", 30.4),
];
pub const HBA_Q1: [(&str, f64); 4] = [
    ("flank_left", 0.98),
    ("HBA2", 0.76),
    ("HBA1", 0.76),
    ("flank_right", 0.98),
];
/// The table `og cn call --locus HBA` writes for those depths: four genes.
pub const HBA_TABLE: &str = concat!(
    "metric\tvalue\n",
    "status\tcalled\n",
    "alpha_globin_genes\t4\n",
    "estimate\t3.96\n",
    "hba2_depth\t30.9\n",
    "hba1_depth\t30.1\n",
    "flank_left_depth\t31.2\n",
    "flank_right_depth\t30.4\n",
    "mapq0_share\t0.240\n",
    "message\tfour alpha-globin genes (aa/aa): no deletion seen\n",
);

/// Paraphase's JSON for one region.
fn region(total: Option<i64>, rsi: &str, extra: &str) -> String {
    let median = total.unwrap_or(4) as f64 * 15.6;
    let total = total.map_or("null".to_owned(), |t| t.to_string());
    format!(
        "{{\"total_cn\": {total}, \"genome_depth\": {}, \"region_depth\": {{\"median\": {}}}, \"region_specific_info\": {rsi}{extra}}}",
        pycompat::fmt::repr(MEAN_DEPTH),
        pycompat::fmt::repr(median)
    )
}

/// Step 35b: a Paraphase JSON with the usual copy numbers (two of each
/// gene, four in each gene family), as `og paraphase table` tables it.
/// Paraphase reads PacBio HiFi; DEMO-001 carries the section so its card is
/// in the pictures. Then step 35's tables from short reads: Parascopy's two
/// SMN1 and two SMN2 copies at exon 7, and [`HBA_TABLE`].
pub fn paralogs(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let haps = |g: &str, a: &str, b: &str| {
        format!(
            ", \"final_haplotypes\": {{\"1\": \"{g}_{a}hap1\", \"2\": \"{g}_{a}hap2\", \"3\": \"{g}_{b}hap1\", \"4\": \"{g}_{b}hap2\"}}"
        )
    };
    let regions = [
        (
            "smn1",
            region(
                Some(4),
                "{\"smn1_cn\": 2, \"smn2_cn\": 2, \"smn_del78_cn\": 0}",
                "",
            ),
        ),
        (
            "hba",
            region(
                Some(4),
                "{\"genotype\": \"aa/aa\", \"sv_called\": {}}",
                &haps("hba", "hba2", "hba1"),
            ),
        ),
        (
            "rccx",
            region(
                Some(4),
                "{\"annotated_alleles\": [\"WT\", \"WT\"], \"phasing_success\": true}",
                "",
            ),
        ),
        (
            "GBA",
            region(
                Some(4),
                "{}",
                &format!(
                    ", \"genes_in_region\": [\"GBA1\", \"GBAP1\"]{}",
                    haps("GBA", "gene1", "gene2")
                ),
            ),
        ),
        ("pms2", region(Some(4), "{\"gene_cn\": 2}", "")),
        ("strc", region(Some(4), "{\"gene_cn\": 2}", "")),
        ("ncf1", region(Some(6), "{\"gene_cn\": 2}", "")),
        (
            "OTOA",
            region(
                Some(4),
                "{}",
                &format!(
                    ", \"genes_in_region\": [\"OTOA\", \"OTOAP1\"]{}",
                    haps("OTOA", "gene1", "gene2")
                ),
            ),
        ),
        (
            "CFC1",
            region(
                Some(4),
                "{}",
                ", \"genes_in_region\": [\"CFC1\", \"CFC1B\"]",
            ),
        ),
        (
            "NEB",
            region(Some(6), "{}", ", \"genes_in_region\": [\"NEB\"]"),
        ),
        (
            "TNXB",
            region(None, "{}", ", \"genes_in_region\": [\"TNXB\"]"),
        ),
    ];
    let body: Vec<String> = regions
        .iter()
        .map(|(k, v)| format!("\"{k}\": {v}"))
        .collect();
    let tmp = Scratch::new("paraphase")?;
    let json = tmp.join(&format!("{s}.paraphase.json"));
    put(&json, &format!("{{{}}}", body.join(", ")))?;
    let genes = ctx.asset("assets/paraphase_genes.tsv");
    let cn = ctx.path(&format!("paralogs/{s}_paraphase_copy_number.tsv"));
    let variants = ctx.path(&format!("paralogs/{s}_paraphase_variants.tsv"));
    std::fs::create_dir_all(ctx.path("paralogs")).map_err(|e| e.to_string())?;
    ctx.og(&[
        "paraphase",
        "table",
        "--sample",
        s,
        "--json",
        arg(&json),
        "--genes",
        arg(&genes),
        "--out-cn",
        arg(&cn),
        "--out-variants",
        arg(&variants),
    ])?;
    put(
        &ctx.path(&format!("paralogs/{s}_smn_copy_number.tsv")),
        "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\nchr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tPASS\t2,2\t33,30\tchr5:70049000-70077000:-\n",
    )?;
    put(
        &ctx.path(&format!("paralogs/{s}_hba_copy_number.tsv")),
        HBA_TABLE,
    )
}

/// Step 50: what `og mei report` writes for a made-up MEGAnE VCF (two Alu
/// and one L1 insertion, away from the catalog) and a made-up depth near
/// 31x over every window of `assets/mei_sites.tsv`, so every catalog row
/// reads "not seen".
pub fn mei(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let mut calls = String::from(
        "chrom\tpos\tend\telement\tsubfamily\tcopies\tfilter\tconfidence\tgenes\tlocation\tbenign_af\tacmg_class\n",
    );
    for (c, p, e, sub) in [
        ("chr3", 41258120, "Alu", "AluYa5"),
        ("chr7", 102931455, "Alu", "AluYb8"),
        ("chr11", 63044871, "L1", "L1HS"),
    ] {
        calls.push_str(&format!(
            "{c}\t{p}\t{p}\t{e}\t{sub}\t1\tPASS\tH\t.\t.\t.\t.\n"
        ));
    }
    let text = read(&ctx.asset("assets/mei_sites.tsv"))?;
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let head = lines.next().ok_or("assets/mei_sites.tsv is empty")?;
    let mut sites = format!("{head}\tstatus\tmean_depth\tcall\n");
    for line in lines {
        let c: Vec<&str> = line.split('\t').collect();
        let num = |i: usize| {
            c.get(i)
                .and_then(|x| x.parse::<i64>().ok())
                .ok_or_else(|| format!("assets/mei_sites.tsv: {line:?}"))
        };
        let width = (num(3)? - num(2)?) as f64;
        let mean = (width * rng.uniform(28.0, 34.0)).round_ties_even() / width;
        sites.push_str(&format!("{line}\tnot seen\t{mean:.1}\t.\n"));
    }
    let s = &ctx.sample;
    put(&ctx.path(&format!("mei/{s}_mei.tsv")), &calls)?;
    put(&ctx.path(&format!("mei/{s}_mei_sites.tsv")), &sites)
}
