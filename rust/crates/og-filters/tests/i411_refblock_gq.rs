//! A low-quality reference call never reaches PharmCAT as a confident 0/0
//! (issue #411, pgx-07; floor D10).
//!
//! Step 07 expands the gVCF's reference blocks into one record per base and
//! `og pgx refblock-fill` passes that stream on to PharmCAT. Before this fix
//! the GQ 20 floor held only for the 0/0 records the filter writes itself at
//! positions longer than one base: a one-base block, or a DeepVariant
//! `RefCall`, with GQ 5 at CYP2C19 rs4244285 (the *2 position) passed through
//! as 0/0, and PharmCAT read the sample as confidently not *2. Now every
//! reference call (GT all 0) below the floor, or without a GQ, is left out
//! of the stream, so PharmCAT reads the position as missing. A variant call
//! passes whatever its GQ: the floor is for reference calls only.

use og_filters::refblock::{MIN_GQ, Outcome, run};

const HEADER: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP\">\n\
##contig=<ID=chr10,length=133797422>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

/// PharmCAT's positions: CYP2C19 *2 (one base) and *17 (one base), and a
/// CYP2C9 deletion (two bases), as in its positions file.
const POSITIONS: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n\
chr10\t94761900\trs12248560\tC\tT\t.\tPASS\tPX=CYP2C19\tGT\t0/0\n\
chr10\t94781859\trs4244285\tG\tA\t.\tPASS\tPX=CYP2C19\tGT\t0/0\n\
chr10\t94942205\trs1\tCT\tC\t.\tPASS\tPX=CYP2C9\tGT\t0/0\n";

const RS4244285: &str = "chr10\t94781859";

struct Run {
    out: String,
    outcome: Outcome,
}

/// The verb on `records` (the gVCF's records over the region, as step 07
/// writes them with `bcftools view -R`) and `stream` (the expanded records
/// on stdin), with `min_gq`.
fn fill(records: &str, stream: &str, min_gq: i128) -> Run {
    let dir = std::env::temp_dir().join(format!(
        "og-i411-refblock-gq-{}-{}",
        std::process::id(),
        std::thread::current()
            .name()
            .unwrap_or("t")
            .replace("::", "-")
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let positions = dir.join("positions.vcf");
    let gvcf = dir.join("gvcf.vcf");
    std::fs::write(&positions, POSITIONS).unwrap();
    std::fs::write(&gvcf, format!("{HEADER}{records}")).unwrap();
    let mut out = Vec::new();
    let (_, outcome) = run(
        positions.to_str().unwrap(),
        gvcf.to_str().unwrap(),
        min_gq,
        format!("{HEADER}{stream}").as_bytes(),
        &mut out,
    )
    .unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    Run {
        out: String::from_utf8(out).unwrap(),
        outcome,
    }
}

/// The one record the output holds at rs4244285, if any.
fn at_rs4244285(out: &str) -> Option<String> {
    let found: Vec<&str> = out
        .lines()
        .filter(|l| l.starts_with(&format!("{RS4244285}\t")))
        .collect();
    assert!(found.len() <= 1, "one record at most: {found:?}");
    found.first().map(|l| (*l).to_owned())
}

/// A one-base reference block over rs4244285 with `gq`: the gVCF record and
/// what `bcftools convert --gvcf2vcf | bcftools view --trim-alt-alleles`
/// makes of it.
fn one_base_block(gq: u32) -> (String, String) {
    (
        format!("{RS4244285}\t.\tG\t<*>\t0\t.\tEND=94781859\tGT:GQ:MIN_DP\t0/0:{gq}:3\n"),
        format!("{RS4244285}\t.\tG\t.\t0\t.\t.\tGT:GQ:MIN_DP\t0/0:{gq}:3\n"),
    )
}

#[test]
fn red_first_a_gq_5_block_at_rs4244285_is_a_no_call() {
    let (records, stream) = one_base_block(5);
    let r = fill(&records, &stream, MIN_GQ);
    assert_eq!(
        at_rs4244285(&r.out),
        None,
        "a one-base block with GQ 5 at CYP2C19 rs4244285 must not reach PharmCAT as 0/0"
    );
}

#[test]
fn a_gq_25_block_at_rs4244285_stays_0_0() {
    let (records, stream) = one_base_block(25);
    let r = fill(&records, &stream, MIN_GQ);
    assert_eq!(
        at_rs4244285(&r.out).as_deref(),
        Some(stream.trim_end()),
        "a block with GQ 25 is confident reference and passes unchanged"
    );
    assert_eq!(r.outcome.weak, 0);
    assert!(r.outcome.weak_at.is_empty());
}

#[test]
fn the_floor_is_d10s_20_and_follows_min_gq() {
    for (gq, min_gq, kept) in [
        (20, MIN_GQ, true),
        (19, MIN_GQ, false),
        (25, 30, false),
        (30, 30, true),
    ] {
        let (records, stream) = one_base_block(gq);
        let r = fill(&records, &stream, min_gq);
        assert_eq!(
            at_rs4244285(&r.out).is_some(),
            kept,
            "GQ {gq} with --min-gq {min_gq}: kept {kept}"
        );
    }
    assert_eq!(MIN_GQ, 20, "the shared floor D10");
}

#[test]
fn refcall_records_take_the_same_floor() {
    // DeepVariant's RefCall: a candidate it called reference, GT 0/0, its ALT
    // trimmed to . by --trim-alt-alleles.
    for (gq, kept) in [(8, false), (30, true)] {
        let records =
            format!("{RS4244285}\t.\tG\tA,<*>\t0\tRefCall\t.\tGT:GQ:MIN_DP\t0/0:{gq}:3\n");
        let stream = format!("{RS4244285}\t.\tG\t.\t0\tRefCall\t.\tGT:GQ:MIN_DP\t0/0:{gq}:3\n");
        let r = fill(&records, &stream, MIN_GQ);
        assert_eq!(
            at_rs4244285(&r.out).is_some(),
            kept,
            "a RefCall with GQ {gq}: kept {kept}"
        );
    }
}

#[test]
fn other_reference_calls_without_confidence_are_left_out() {
    for (what, sample) in [
        ("haploid, GQ 5", "0:5:3"),
        ("phased, GQ 5", "0|0:5:3"),
        ("no GQ value", "0/0:.:3"),
    ] {
        let stream = format!("{RS4244285}\t.\tG\t.\t0\t.\t.\tGT:GQ:MIN_DP\t{sample}\n");
        let r = fill("", &stream, MIN_GQ);
        assert_eq!(at_rs4244285(&r.out), None, "{what}: left out");
    }
    let stream = format!("{RS4244285}\t.\tG\t.\t0\t.\t.\tGT:MIN_DP\t0/0:3\n");
    let r = fill("", &stream, MIN_GQ);
    assert_eq!(at_rs4244285(&r.out), None, "no GQ field: left out");
}

#[test]
fn variant_calls_and_no_sample_lines_pass_whatever_the_gq() {
    for sample in ["0/1:5:3", "1/1:3:3", "1:2:3"] {
        let stream = format!("{RS4244285}\t.\tG\tA\t10\tPASS\t.\tGT:GQ:MIN_DP\t{sample}\n");
        let r = fill("", &stream, MIN_GQ);
        assert_eq!(
            at_rs4244285(&r.out).as_deref(),
            Some(stream.trim_end()),
            "a variant call {sample} passes unchanged"
        );
    }
}

#[test]
fn the_outcome_counts_them_and_names_pharmcat_positions() {
    // A GQ 5 block over three bases: rs4244285 and the two after it.
    let records = format!("{RS4244285}\t.\tG\t<*>\t0\t.\tEND=94781861\tGT:GQ:MIN_DP\t0/0:5:3\n");
    let stream: String = (94_781_859..=94_781_861)
        .map(|p| format!("chr10\t{p}\t.\tN\t.\t0\t.\t.\tGT:GQ:MIN_DP\t0/0:5:3\n"))
        .collect();
    let r = fill(&records, &stream, MIN_GQ);
    assert_eq!(
        r.out.lines().filter(|l| !l.starts_with('#')).count(),
        0,
        "all three left out"
    );
    assert_eq!(r.outcome.weak, 3);
    assert_eq!(r.outcome.weak_at, vec!["CYP2C19 chr10:94781859".to_owned()]);
}

#[test]
fn a_low_gq_block_inside_a_multi_base_span_still_fills_nothing() {
    // The CYP2C9 deletion at 94942205 (CT): a GQ 5 block over both bases. No
    // 0/0 is written for it, and its one-base records are left out too.
    let records = "chr10\t94942200\t.\tA\t<*>\t0\t.\tEND=94942210\tGT:GQ:MIN_DP\t0/0:5:3\n";
    let stream: String = (94_942_200..=94_942_210)
        .map(|p| format!("chr10\t{p}\t.\tN\t.\t0\t.\t.\tGT:GQ:MIN_DP\t0/0:5:3\n"))
        .collect();
    let r = fill(records, &stream, MIN_GQ);
    assert!(r.outcome.filled.is_empty(), "nothing filled");
    assert_eq!(
        r.out.lines().filter(|l| !l.starts_with('#')).count(),
        0,
        "no record reaches PharmCAT"
    );
    assert_eq!(r.outcome.weak_at, vec!["CYP2C9 chr10:94942205".to_owned()]);
}
