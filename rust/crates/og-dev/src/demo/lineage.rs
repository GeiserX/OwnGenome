//! The demo sample's ancestry chapter: polygenic scores and the ancestry
//! placement (step 25), fine-grained ancestry (26b), the mitochondrial and
//! Y haplogroups (12, 37), the Y-STR markers (45), archaic (43) and ancient
//! (44) ancestry, and runs of homozygosity (11).

use super::{
    Ctx, Rng, Scratch, arg, autosomes, cell, chrom_index, copy, put, put_gz, read, table,
    versions_env,
};
use std::collections::HashMap;

/// `rng.sample("ACGT", 2)`, as a string pair.
fn allele_pair(rng: &mut Rng) -> (String, String) {
    let (a, b) = rng.two_bases();
    (a.to_string(), b.to_string())
}

/// `;`-joined entries, each stripped; empty when the cell is blank.
fn semis(x: &str) -> String {
    if x.trim().is_empty() {
        String::new()
    } else {
        x.split(';').map(str::trim).collect::<Vec<_>>().join(";")
    }
}

/// Step 25 through pgsc_calc's files: the harmonised score files (their
/// headers feed the run manifest), the score folder `og prs format` would
/// write from them (without the allele list, which nothing here reads:
/// every row is autosomal, so it is written, an A/T or C/G pair counts as
/// strand-ambiguous, and no chrX row is left out), pgsc_calc's results for
/// sampleset "sample" with the 1000 Genomes panel, then `og prs table`.
pub fn prs(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let sample = ctx.sample.as_str();
    let listed = table(&read(&ctx.asset("assets/pgs_scores.tsv"))?);
    // pgs_id -> (trait, applies_to, evaluations), in file order.
    let mut ids: Vec<String> = Vec::new();
    let mut info: HashMap<String, (String, String, Option<String>)> = HashMap::new();
    for r in &listed {
        let pid = cell(r, "pgs_id").trim();
        if pid.is_empty() {
            continue;
        }
        let keys = ["ppm_id", "ancestry", "effect"];
        let evals = keys.iter().any(|k| !cell(r, k).trim().is_empty()).then(|| {
            keys.iter()
                .map(|k| semis(cell(r, k)))
                .collect::<Vec<_>>()
                .join("\t")
        });
        if !info.contains_key(pid) {
            ids.push(pid.to_owned());
        }
        info.insert(
            pid.to_owned(),
            (
                cell(r, "trait_reported").trim().to_owned(),
                cell(r, "applies_to").trim().to_lowercase(),
                evals,
            ),
        );
    }
    if ids.is_empty() {
        return Err("assets/pgs_scores.tsv lists no score".to_owned());
    }
    let harmonised = ctx.genome_dir.join("prs_scores");
    let work = ctx.path("prs/pgsc_calc");
    let scores = work.join("scores");
    let mut counts: HashMap<String, (usize, usize)> = HashMap::new();
    for pid in &ids {
        let trait_ = &info[pid].0;
        let n = rng.int(250, 4000);
        let mut rows = String::new();
        let mut written = String::new();
        let mut ambiguous = 0;
        for _ in 0..n {
            let (c, length) = autosomes()[rng.index(22)];
            let (ea, oa) = allele_pair(rng);
            let pos = rng.int(1, length);
            let w = format!("{:.5}", rng.gauss(0.0, 0.05));
            let hm_pos = rng.int(1, length);
            let c = &c[3..];
            rows.push_str(&format!("{c}\t{pos}\t{ea}\t{oa}\t{w}\t{c}\t{hm_pos}\n"));
            written.push_str(&format!("{c}\t{hm_pos}\t{ea}\t{oa}\t{w}\n"));
            if ["AT", "TA", "CG", "GC"].contains(&format!("{ea}{oa}").as_str()) {
                ambiguous += 1;
            }
        }
        counts.insert(pid.clone(), (n as usize, ambiguous));
        put_gz(
            &harmonised.join(format!("{pid}.txt.gz")),
            &format!(
                "#format_version=2.0\n#pgs_id={pid}\n#pgs_name={pid}\n#trait_reported={trait_}\n#variants_number={n}\n#HmPOS_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\n{rows}"
            ),
        )?;
        let shown = if trait_.is_empty() { pid } else { trait_ }.replace('=', "-");
        put_gz(
            &scores.join(format!("{pid}.txt.gz")),
            &format!(
                "#pgs_id={pid}\n#pgs_name={pid}\n#trait_reported={shown}\n#genome_build=GRCh38\nchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n{written}"
            ),
        )?;
    }
    let mut sorted = ids.clone();
    sorted.sort();
    let mut chrx = String::from("pgs_id\tchrx_rows\tchrx\n");
    let mut notes = String::from("pgs_id\tapplies_to\tambiguous_rows\tambiguous\n");
    let mut evals = String::from("pgs_id\tppm_id\tancestry\teffect\n");
    for p in &sorted {
        chrx.push_str(&format!("{p}\t0\tleft out\n"));
        notes.push_str(&format!("{p}\t{}\t{}\tdropped\n", info[p].1, counts[p].1));
        if let Some(e) = &info[p].2 {
            evals.push_str(&format!("{p}\t{e}\n"));
        }
    }
    put(&scores.join("chrx_rows.tsv"), &chrx)?;
    put(&scores.join("score_notes.tsv"), &notes)?;
    put(&scores.join("score_evaluations.tsv"), &evals)?;

    let res = work.join("results/sample");
    let mut summ = vec![
        "dataset,accession,ambiguous,is_multiallelic,match_flipped,duplicate_best_match,duplicate_ID,match_IDs,match_status,count,score_pass,match_rate".to_owned(),
    ];
    let mut score = Vec::new();
    for pid in &sorted {
        let total = counts[pid].0 as f64;
        let m = (total * rng.uniform(0.86, 0.99)).round_ties_even() as i64;
        let rate = m as f64 / total;
        summ.push(format!(
            "sample,{pid},false,false,false,false,false,NA,matched,{m},true,{rate:.2}"
        ));
        summ.push(format!(
            "sample,{pid},false,false,false,false,false,NA,unmatched,{},true,{rate:.2}",
            total as i64 - m
        ));
        score.push(format!(
            "sample\t{sample}\t{sample}\t{pid}\t{:.4}\t{:.3}\t{:.2}",
            rng.gauss(0.0, 0.4),
            rng.gauss(0.0, 1.0),
            rng.uniform(18.0, 82.0)
        ));
    }
    put(
        &res.join("match/sample_summary.csv"),
        &(summ.join("\n") + "\n"),
    )?;
    put_gz(
        &res.join("score/sample_pgs.txt.gz"),
        &format!(
            "sampleset\tFID\tIID\tPGS\tSUM\tZ_MostSimilarPop\tpercentile_MostSimilarPop\n{}\n",
            score.join("\n")
        ),
    )?;
    let pcs: Vec<String> = (0..10)
        .map(|_| format!("{:.4}", rng.gauss(0.0, 0.01)))
        .collect();
    let pc_names: Vec<String> = (1..=10).map(|i| format!("PC{i}")).collect();
    put_gz(
        &res.join("score/sample_popsimilarity.txt.gz"),
        &format!(
            "sampleset\tFID\tIID\t{}\tRF_P_AFR\tRF_P_AMR\tRF_P_EAS\tRF_P_EUR\tRF_P_SAS\tMostSimilarPop\tMostSimilarPop_LowConfidence\tREFERENCE\nsample\t{sample}\t{sample}\t{}\t0.00\t0.02\t0.00\t0.97\t0.01\tEUR\tFalse\tFalse\n",
            pc_names.join("\t"),
            pcs.join("\t")
        ),
    )?;
    std::fs::create_dir_all(ctx.path("ancestry")).map_err(|e| e.to_string())?;
    let panel =
        versions_env(&ctx.repo, "PGSC_PANEL").unwrap_or_else(|_| "pgsc_1000G_v1".to_owned());
    let results = work.join("results");
    let ancestry = ctx.path(&format!("ancestry/{sample}_ancestry.tsv"));
    let out = ctx.path(&format!("prs/{sample}_prs_summary.tsv"));
    ctx.og(&[
        "prs",
        "table",
        "--sample",
        sample,
        "--results",
        arg(&results),
        "--sampleset",
        "sample",
        "--scores",
        arg(&scores),
        "--input-kind",
        "gvcf",
        "--panel",
        &panel,
        "--ancestry-out",
        arg(&ancestry),
        "--out",
        arg(&out),
    ])
}

/// Three fixed runs on chr11, the chromosome the live demo's map checks
/// open: with its ClinVar finding they make the four features in view there,
/// and one of them is the run whose row flies the map.
const ROH_CHR11: [(i64, i64); 3] = [
    (60_960_259, 64_330_971),
    (70_291_155, 71_284_241),
    (74_599_677, 76_877_148),
];

/// Step 11: about 44 Mb of runs of homozygosity, bcftools roh's RG rows:
/// the fixed chr11 runs, then random runs on the other autosomes.
pub fn roh(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let chr11 = chrom_index("chr11");
    let mut segs: Vec<(usize, &str, i64, i64, i64)> = ROH_CHR11
        .iter()
        .map(|&(s, e)| (chr11, "chr11", s, e, e - s))
        .collect();
    let mut total: i64 = segs.iter().map(|x| x.4).sum();
    while total < 44_000_000 {
        let i = rng.index(22);
        if i == chr11 {
            continue;
        }
        let (c, length) = autosomes()[i];
        let size = rng.uniform(0.5e6, 3.9e6) as i64;
        let start = rng.int(1_000_000, length - size - 1_000_000);
        segs.push((i, c, start, start + size, size));
        total += size;
    }
    segs.sort();
    let mut text = String::from(
        "# RG, Regions [2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality\n",
    );
    for (_, c, s, e, n) in segs {
        text.push_str(&format!(
            "RG\t{}\t{c}\t{s}\t{e}\t{n}\t{}\t{:.1}\n",
            ctx.sample,
            n / 1400,
            rng.uniform(40.0, 90.0)
        ));
    }
    put(&ctx.path(&format!("vcf/{}_roh.txt", ctx.sample)), &text)
}

/// Steps 12 and 37: haplogrep's mitochondrial haplogroup and haplocheck's
/// contamination line, and for a male Yleaf's Y haplogroup with its
/// placement on the current tree (`og y tree-place`'s columns: "." for an
/// empty SNP list, and its default depth of 10).
pub fn haplogroups(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    put(
        &ctx.path(&format!("mito/{s}_haplogroup.txt")),
        &format!(
            "\"SampleID\"\t\"Haplogroup\"\t\"Rank\"\t\"Quality\"\t\"Range\"\n\"{s}\"\t\"H1c3\"\t\"1\"\t\"0.9312\"\t\"1-16569;\"\n"
        ),
    )?;
    put(
        &ctx.path(&format!("mito/{s}_haplocheck.txt")),
        &format!(
            "\"Sample\"\t\"Contamination Status\"\t\"Contamination Level\"\t\"Distance\"\t\"Sample Coverage\"\n\"{s}\"\t\"NO\"\t\"ND\"\t\"0\"\t\"2861\"\n"
        ),
    )?;
    if !ctx.male {
        return Ok(()); // no Y chromosome: step 37 has nothing to place
    }
    put(
        &ctx.path(&format!("y_haplogroup/{s}_y_haplogroup.txt")),
        &format!(
            "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\n{s}_sorted\tI-M253\tM253\t4021876\t1187\t0.982\t1.0\t0.991\t0.982\n"
        ),
    )?;
    put(
        &ctx.path(&format!("y_haplogroup/{s}_y_tree.txt")),
        &format!(
            "Sample_name\tHg\tLevel\tPath_derived\tPath_ancestral\tPath_no_call\tAncestral_in_terminal\tSNPs_called\tMin_depth\tTree\tPositions\tOutvoted_ancestral\n{s}\tI-M253\t5\t41\t0\t3\t.\t1187\t10\tYFull YTree v14.06.0\tYleaf hg38 position tables\t.\n"
        ),
    )?;
    put(
        &ctx.path(&format!("y_haplogroup/{s}_y_tree_path.tsv")),
        "Level\tBranch\tSNPs\tOn_tables\tDerived\tAncestral\tNo_call\tAncestral_SNPs\n1\tCT\t12\t12\t11\t0\t1\t.\n2\tF\t9\t9\t9\t0\t0\t.\n3\tIJ\t8\t8\t7\t0\t1\t.\n4\tI\t7\t7\t7\t0\t0\t.\n5\tI-M253\t8\t8\t7\t0\t1\t.\n",
    )
}

/// Step 26b's table (ADMIXTURE fractions with their bootstrap range, then
/// the painted shares) and its painting, written by `og ancestry-fine
/// report` from made-up ADMIXTURE and RFMix output. Mostly West Eurasian,
/// as step 25's EUR placement.
pub fn ancestry_fine(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let s = &ctx.sample;
    let continent: &[(&str, f64)] = &[
        ("WEA", 0.942),
        ("CSA", 0.031),
        ("AFR", 0.012),
        ("EAS", 0.008),
        ("AMR", 0.005),
        ("OCE", 0.002),
    ];
    let region: &[(&str, f64)] = &[
        ("NWE", 0.463),
        ("IBE", 0.274),
        ("ITA", 0.118),
        ("NEE", 0.052),
        ("SAR", 0.009),
        ("CAU", 0.011),
        ("LEV", 0.008),
        ("BDN", 0.002),
        ("NAF", 0.005),
        ("CSA", 0.031),
        ("WAF", 0.007),
        ("EAF", 0.003),
        ("CAF", 0.001),
        ("SAF", 0.001),
        ("EAS", 0.006),
        ("SIB", 0.002),
        ("AMR", 0.005),
        ("OCE", 0.002),
    ];
    let work = Scratch::new("ancestry-fine")?;
    let (model, adm, paint, out) = (
        work.join("model"),
        work.join("admixture"),
        work.join("paint"),
        work.join("out"),
    );
    let mut levels = String::new();
    for (level, fr) in [("continent", continent), ("region", region)] {
        let codes: Vec<&str> = fr.iter().map(|(c, _)| *c).collect();
        levels.push_str(&format!("{level}\t{}\t{}\n", codes.len(), codes.join(",")));
        let q: Vec<String> = fr.iter().map(|(_, f)| format!("{f:.6}")).collect();
        put(
            &adm.join(format!("{s}_{level}.{}.Q", codes.len())),
            &(q.join(" ") + "\n"),
        )?;
        // 20 bootstrap replicates, each fraction moved by a few percent of its spread.
        let mut boot = String::new();
        for _ in 0..20 {
            let row: Vec<f64> = fr
                .iter()
                .map(|(_, f)| {
                    let sd = ((f * (1.0 - f)).sqrt() * 0.06).max(0.002);
                    (f + rng.gauss(0.0, sd)).max(0.0)
                })
                .collect();
            let sum: f64 = row.iter().sum();
            let cells: Vec<String> = row.iter().map(|v| format!("{:.6}", v / sum)).collect();
            boot.push_str(&(cells.join(" ") + "\n"));
        }
        put(&adm.join(format!("{s}_{level}.boot.Q")), &boot)?;
    }
    put(&model.join("levels.tsv"), &levels)?;
    // The painting: each chromosome in windows of 3 to 15 Mb, each
    // haplotype West Eurasian but for a few windows.
    let groups = ["AFR", "CSA", "WEA"];
    let order: Vec<String> = groups
        .iter()
        .enumerate()
        .map(|(k, g)| format!("{g}={k}"))
        .collect();
    for (i, (_, n)) in autosomes().iter().enumerate() {
        let c = format!("chr{}", i + 1);
        let mut rows = vec![
            format!("#Subpopulation order/codes: {}", order.join("\t")),
            format!("#chm\tspos\tepos\tsgpos\tegpos\tn snps\t{s}.0\t{s}.1"),
        ];
        let mut pos = rng.int(200_000, 1_500_000);
        let end_of = n - rng.int(200_000, 1_500_000);
        while pos < end_of - 2_000_000 {
            let end = end_of.min(pos + rng.int(3_000_000, 15_000_000));
            let mut hap = || {
                if rng.random() < 0.85 {
                    2
                } else {
                    *rng.pick(&[1, 1, 0])
                }
            };
            let (h0, h1) = (hap(), hap());
            rows.push(format!(
                "{c}\t{pos}\t{end}\t{:.3}\t{:.3}\t100\t{h0}\t{h1}",
                pos as f64 / 1e6,
                end as f64 / 1e6
            ));
            pos = end;
        }
        put(
            &paint.join(format!("{c}.msp.tsv")),
            &(rows.join("\n") + "\n"),
        )?;
    }
    ctx.og(&[
        "ancestry-fine",
        "report",
        "--sample",
        s,
        "--model",
        arg(&model),
        "--admixture",
        arg(&adm),
        "--paint",
        arg(&paint),
        "--panel",
        "demo",
        "--out",
        arg(&out),
    ])?;
    for name in [
        format!("{s}_ancestry_fine.tsv"),
        format!("{s}_painting.png"),
    ] {
        copy(
            &out.join(&name),
            &ctx.path(&format!("ancestry_fine/{name}")),
        )?;
    }
    Ok(())
}

/// Step 45 (a male only): one row per marker of `assets/ystr_panel.tsv` in
/// the columns of `og y str-report`. The multi-copy, absent and too-long
/// markers (a repeat with its 10-base flanks longer than 150 bp reads) read
/// not callable, with og's notes, as in a real run; every other marker gets
/// a made-up count at confidence high.
pub fn y_str(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    if !ctx.male {
        return Ok(());
    }
    let panel = table(&read(&ctx.asset("assets/ystr_panel.tsv"))?);
    let cols = [
        "order",
        "panel",
        "marker",
        "count",
        "confidence",
        "reads",
        "stutter_reads",
        "flank_indel_reads",
        "q",
        "note",
    ];
    let mut out = cols.join("\t") + "\n";
    for r in &panel {
        let regions = cell(r, "regions");
        let mut spans: Vec<(i64, i64)> = Vec::new();
        if regions != "." {
            for part in regions.split(';') {
                let range = part
                    .split('/')
                    .next()
                    .and_then(|x| x.split_once(':'))
                    .and_then(|(_, r)| r.split_once('-'))
                    .and_then(|(b, e)| Some((b.parse().ok()?, e.parse().ok()?)))
                    .ok_or_else(|| format!("assets/ystr_panel.tsv: region {part:?}"))?;
                spans.push(range);
            }
        }
        let kind = cell(r, "kind");
        if kind == "span" && !spans.is_empty() {
            spans = vec![(spans[0].0, spans[spans.len() - 1].1)];
        }
        let longest = spans.iter().map(|(b, e)| e - b + 1).max().unwrap_or(0);
        let mut row: HashMap<&str, String> = HashMap::new();
        row.insert("order", cell(r, "order").to_owned());
        row.insert("panel", cell(r, "panel").to_owned());
        row.insert("marker", cell(r, "marker").to_owned());
        row.insert("confidence", "not callable".to_owned());
        let note = if kind == "multicopy" {
            "not callable from short reads: a multi-copy marker in duplicated Y sequence".to_owned()
        } else if kind == "absent" {
            "not callable: no region for it in the HipSTR marker set".to_owned()
        } else if longest + 2 * 10 > 150 {
            format!(
                "not callable from short reads: a {longest} bp repeat is longer than 150 bp reads can span with 10 bases on each side"
            )
        } else {
            row.insert("count", rng.int(10, 30).to_string());
            row.insert("confidence", "high".to_owned());
            row.insert("reads", rng.int(18, 34).to_string());
            row.insert("stutter_reads", "0".to_owned());
            row.insert("flank_indel_reads", "0".to_owned());
            row.insert("q", "0.99".to_owned());
            ".".to_owned()
        };
        row.insert("note", note);
        let cells: Vec<&str> = cols
            .iter()
            .map(|k| row.get(k).map_or("", String::as_str))
            .collect();
        out.push_str(&(cells.join("\t") + "\n"));
    }
    put(&ctx.path(&format!("y_str/{}_y_str.tsv", ctx.sample)), &out)
}

/// Step 43: the archaic report's key/value summary, with made-up shares
/// near the published non-African range and five made-up longest segments.
pub fn archaic(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let analysed: i64 = 2_870_000_000;
    let (nea_n, den_n) = (rng.int(420, 520), rng.int(2, 6));
    let nea_bp = (analysed as f64 * rng.uniform(0.017, 0.021)) as i64;
    let den_bp = (analysed as f64 * rng.uniform(0.0002, 0.0006)) as i64;
    let pct = |bp: i64| format!("{:.2}", 100.0 * bp as f64 / analysed as f64);
    let mut rows: Vec<(String, String)> = [
        ("sample", ctx.sample.clone()),
        ("method", "IBDmix".to_owned()),
        (
            "archaic_genomes",
            "Altai Neanderthal, Vindija 33.19 Neanderthal, Altai Denisovan".to_owned(),
        ),
        ("population", "EUR".to_owned()),
        ("min_lod", "4".to_owned()),
        ("min_length_bp", "50000".to_owned()),
        ("analysed_bp", analysed.to_string()),
        ("analysed_chromosomes", "22".to_owned()),
        ("analysed_sites", rng.int(5_500_000, 6_500_000).to_string()),
        ("genes", "GENCODE protein-coding".to_owned()),
        ("neanderthal_segments", nea_n.to_string()),
        ("neanderthal_bp", nea_bp.to_string()),
        ("neanderthal_percent", pct(nea_bp)),
        ("denisovan_segments", den_n.to_string()),
        ("denisovan_bp", den_bp.to_string()),
        ("denisovan_percent", pct(den_bp)),
    ]
    .into_iter()
    .map(|(k, v)| (k.to_owned(), v))
    .collect();
    let segs = [
        ("chr1", "Neanderthal"),
        ("chr3", "Neanderthal"),
        ("chr6", "Neanderthal"),
        ("chr9", "Neanderthal"),
        ("chr12", "Denisovan"),
    ];
    for (i, (c, anc)) in segs.iter().enumerate() {
        let i = i + 1;
        let start = rng.int(10_000_000, 90_000_000);
        let n = (rng.uniform(0.9, 1.3 - 0.08 * i as f64) * 1e6) as i64;
        let who = if *anc == "Denisovan" {
            "Altai Denisovan"
        } else {
            "Altai/Vindija"
        };
        rows.push((
            format!("longest_{i}"),
            format!(
                "{c}:{start}-{} {anc} {:.2} Mb ({who}); genes: none",
                start + n,
                n as f64 / 1e6
            ),
        ));
    }
    let mut text = String::from("key\tvalue\n");
    for (k, v) in rows {
        text.push_str(&format!("{k}\t{v}\n"));
    }
    put(
        &ctx.path(&format!("archaic/{}_archaic.tsv", ctx.sample)),
        &text,
    )
}

/// Step 44: what `og ancient report` writes (the nearest-group table, the
/// info table, the text and the PC plot) for a made-up panel: present-day
/// people around the origin, ancient people of a few AADR-style group
/// labels each around its own point, and the person among them.
pub fn ancient(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    let groups = [
        (
            "England_BellBeaker",
            4200,
            "United Kingdom",
            "Amesbury Down",
        ),
        ("Germany_CordedWare", 4600, "Germany", "Esperstedt"),
        ("Spain_BellBeaker", 4200, "Spain", "Camino de las Yeseras"),
        ("Czechia_EBA_Unetice", 3900, "Czech Republic", "Moravia"),
        ("Hungary_LBA", 3250, "Hungary", "Tiszafüred"),
        ("Italy_Imperial_Rome", 1800, "Italy", "Rome"),
        ("Sweden_Viking", 1030, "Sweden", "Sigtuna"),
        ("France_Neolithic", 5700, "France", "Gurgy"),
        ("Turkey_Neolithic", 8200, "Turkey", "Barcın"),
        ("Russia_Yamnaya_Samara", 4950, "Russia", "Samara"),
    ];
    let person = (0.01234_f64, -0.00456_f64);
    let evals = "#eigvals: 9.000 4.000 2.000 1.000\n";
    let mut moderns = String::new();
    for i in 0..300 {
        let pcs: Vec<String> = [0.03, 0.025, 0.02, 0.015]
            .iter()
            .map(|sd| format!("{:.6}", rng.gauss(0.0, *sd)))
            .collect();
        moderns.push_str(&format!(
            "m{i} {} {}\n",
            pcs.join(" "),
            ["French", "Sardinian", "Basque"][i % 3]
        ));
    }
    let mut ancients = String::new();
    let mut anc_rows =
        String::from("id\tgroup\tdate_bp\tdate_sd\tfull_date\tlocality\tcountry\tlat\tlon\tsnps\n");
    for (i, (g, bp, country, loc)) in groups.iter().enumerate() {
        let spread = 0.004 + 0.006 * i as f64;
        let cx = person.0 + rng.uniform(-spread, spread);
        let cy = person.1 + rng.uniform(-spread, spread);
        for j in 0..rng.int(3, 12) {
            let aid = format!("a{i}_{j}");
            ancients.push_str(&format!(
                "{aid} {:.6} {:.6} 0.000000 0.000000 Ancient\n",
                cx + rng.gauss(0.0, 0.004),
                cy + rng.gauss(0.0, 0.004)
            ));
            let date = bp + rng.int(-150, 150);
            anc_rows.push_str(&format!(
                "{aid}\t{g}\t{date}\t50\t.\t{loc}\t{country}\t50\t10\t500000\n"
            ));
        }
    }
    let work = Scratch::new("ancient")?;
    let (panel, run, anc, stats) = (
        work.join("panel.evec"),
        work.join("run.evec"),
        work.join("ancients.tsv"),
        work.join("stats.tsv"),
    );
    put(&panel, &format!("{evals}{moderns}{ancients}"))?;
    put(
        &run,
        &format!(
            "{evals}{moderns}{} {:.6} {:.6} 0.000000 0.000000 Person\n",
            ctx.sample, person.0, person.1
        ),
    )?;
    put(&anc, &anc_rows)?;
    put(
        &stats,
        &format!(
            "panel_sites\t597573\nsites_genotyped\t560000\nother_allele\t120\nsites_used\t{}\n",
            rng.int(520_000, 560_000)
        ),
    )?;
    let s = &ctx.sample;
    let out = |x: &str| ctx.path(&format!("ancient/{s}_{x}"));
    let (tsv, info, txt, svg) = (
        out("ancient.tsv"),
        out("ancient_info.tsv"),
        out("ancient.txt"),
        out("ancient_pca.svg"),
    );
    std::fs::create_dir_all(ctx.path("ancient")).map_err(|e| e.to_string())?;
    ctx.og(&[
        "ancient",
        "report",
        "--panel-evec",
        arg(&panel),
        "--run-evec",
        arg(&run),
        "--ancients",
        arg(&anc),
        "--sample",
        s,
        "--release",
        "v62.0",
        "--stats",
        arg(&stats),
        "--out-tsv",
        arg(&tsv),
        "--out-info",
        arg(&info),
        "--out-txt",
        arg(&txt),
        "--out-svg",
        arg(&svg),
    ])
}
