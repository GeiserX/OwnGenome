//! The demo sample's ClinVar results: the ClinVar screen (step 06), with
//! the son's planted result to act on, and the couple carrier screen
//! (step 39) the carrier chapter shows.

use super::{Ctx, Rng, Scratch, arg, cell, chrom_index, copy, put, put_gz, read, table};

/// The son's results to act on (DEMO-003, finding "act"): a pathogenic
/// ClinVar hit in KCNH2, an ACMG secondary-findings gene that is on no
/// carrier list, and an AnnotSV class 5 deletion over it. Public gene name,
/// invented positions and alleles. The live-site tests check that his page
/// names neither until each is opened on purpose, and that the family page
/// keeps his places closed with them, the compound-het pairs of
/// `tests/fixtures/live/demo/DEMO-FAMILY` too.
pub struct Planted {
    pub gene: &'static str,
    /// chrom, pos, ref, alt.
    pub hit: (&'static str, i64, &'static str, &'static str),
    /// chrom, start, end.
    pub sv: (&'static str, i64, i64),
}

pub const PLANTED: Planted = Planted {
    gene: "KCNH2",
    hit: ("chr7", 150_960_000, "C", "T"),
    sv: ("chr7", 150_955_000, 150_965_000),
};

/// Invented ClinVar genes and the chromosome of each hit.
pub const INVENTED_GENES: [(&str, &str); 3] =
    [("DEMOA", "chr2"), ("DEMOB", "chr11"), ("DEMOC", "chr16")];

/// Three heterozygous hits in genes with invented names (the picture stops
/// above the hits table, and no real gene is tied to an invented person).
/// DEMO-003 gets a fourth, the planted result to act on: a real gene,
/// needed for the secondary-findings rule, at an invented position. Also
/// the ClinVar file itself, cut down to the header lines the run manifest
/// reads.
pub fn clinvar(ctx: &Ctx, _rng: &mut Rng) -> Result<(), String> {
    let [(a, chr_a), (b, chr_b), (c, chr_c)] = INVENTED_GENES;
    let multi = "criteria_provided,_multiple_submitters,_no_conflicts";
    let mut hits: Vec<(&str, i64, &str, &str, String, &str, &str)> = vec![
        (
            chr_a,
            47403191,
            "C",
            "T",
            format!("{a}:90001"),
            "Pathogenic",
            multi,
        ),
        (
            chr_b,
            108259012,
            "G",
            "A",
            format!("{b}:90002"),
            "Likely_pathogenic",
            "criteria_provided,_single_submitter",
        ),
        (
            chr_c,
            2087914,
            "CT",
            "C",
            format!("{c}:90003"),
            "Pathogenic",
            multi,
        ),
    ];
    if ctx.has("act") {
        let (c, p, r, a) = PLANTED.hit;
        hits.push((
            c,
            p,
            r,
            a,
            format!("{}:3757", PLANTED.gene),
            "Pathogenic",
            multi,
        ));
        hits.sort_by_key(|h| (chrom_index(h.0), h.1));
    }
    let mut body = ctx.vcf_header(&[]);
    for (i, (c, p, r, a, g, s, rev)) in hits.iter().enumerate() {
        body.push_str(&format!(
            "{c}\t{p}\t{}\t{r}\t{a}\t50\tPASS\tGENEINFO={g};CLNSIG={s};CLNREVSTAT={rev}\tGT\t0/1\n",
            90000 + i
        ));
    }
    put(
        &ctx.path(&format!("clinvar/{}_clinvar_hits.vcf", ctx.sample)),
        &body,
    )?;
    put_gz(
        &ctx.genome_dir.join("clinvar/clinvar.vcf.gz"),
        "##fileformat=VCFv4.1\n##fileDate=2026-09-28\n##source=ClinVar\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
}

/// Step 39 with an invented partner, DEMO-002: each carries one change in a
/// different recessive gene of the ACMG carrier list, so no gene is shared
/// and the card names none. The rows are `og couple screen`'s, fed a
/// made-up ClinVar screen of each partner. Only DEMO-001 carries the table:
/// the family page reads one copy of it.
pub fn couple(ctx: &Ctx, rng: &mut Rng) -> Result<(), String> {
    if ctx.sample != super::SAMPLE {
        return Ok(());
    }
    let genes_tsv = ctx.asset("assets/carrier_genes.tsv");
    let rows = table(&read(&genes_tsv)?);
    let mut ar: Vec<&Vec<(String, String)>> = rows
        .iter()
        .filter(|r| cell(r, "inheritance") == "AR" && cell(r, "lists").contains("ACMG"))
        .collect();
    ar.sort_by(|a, b| cell(a, "gene").cmp(cell(b, "gene")));
    if ar.len() < 2 {
        return Err(
            "assets/carrier_genes.tsv has fewer than two AR genes of the ACMG list".to_owned(),
        );
    }
    let work = Scratch::new("couple")?;
    let mut args: Vec<String> = Vec::new();
    for (who, i) in [ctx.sample.as_str(), "DEMO-002"]
        .into_iter()
        .zip(rng.sample_idx(ar.len(), 2))
    {
        let g = ar[i];
        let num = |k: &str| {
            cell(g, k)
                .parse::<i64>()
                .map_err(|_| format!("assets/carrier_genes.tsv: {} {k}", cell(g, "gene")))
        };
        let pos = rng.int(num("start")?, num("end")?);
        let (r, a) = rng.two_bases();
        let vcf = work.join(&format!("{who}.vcf"));
        put(
            &vcf,
            &format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n{}\t{pos}\t.\t{r}\t{a}\t.\tPASS\tGENEINFO={}:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\tGT\t0/1\n",
                cell(g, "chrom"),
                cell(g, "gene")
            ),
        )?;
        args.push("--clinvar".to_owned());
        args.push(format!("{who}={}", arg(&vcf)));
    }
    let pedigree = work.join("pedigree.tsv");
    put(
        &pedigree,
        &format!(
            "sample\tsex\tfather\tmother\tpartner\n{s}\t{}\t.\t.\tDEMO-002\nDEMO-002\tfemale\t.\t.\t{s}\n",
            ctx.sex(),
            s = ctx.sample
        ),
    )?;
    let par = ctx.asset("assets/par_grch38.bed");
    let mut call: Vec<&str> = vec![
        "couple",
        "screen",
        "--pedigree",
        arg(&pedigree),
        "--genes",
        arg(&genes_tsv),
        "--par-bed",
        arg(&par),
    ];
    call.extend(args.iter().map(String::as_str));
    call.extend(["--out-dir", arg(&work.0)]);
    ctx.og(&call)?;
    let name = format!("{}_couple_carrier.tsv", ctx.sample);
    copy(&work.join(&name), &ctx.path(&format!("couple/{name}")))
}
