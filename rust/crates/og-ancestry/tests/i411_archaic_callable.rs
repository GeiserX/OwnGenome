//! `og archaic genotypes` reads a table site with no VCF record as the
//! reference only where a confident reference block of the gVCF covers it.
//! Before this, every such site was written as 0 (REF/REF), so a stretch
//! with no reads at all could score as a long archaic segment.

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use og_ancestry::archaic::cli;
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

const TABLE_HEAD: &str = "chrom\tpos\tref\talt\tAltai\tVindija\tDenisova\tAF_ALL\tAF_AFR\tAF_AMR\tAF_EAS\tAF_EUR\tAF_SAS\n";
const VCF_HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

/// A fresh folder for one test, with an empty `table/`.
fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-i411-archaic-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("table")).unwrap();
    d
}

/// One chromosome of the table: an A>G row at each position, the three
/// archaic genomes 0/0 and the ALT at 0.95 in every population, the shape
/// where a site read as the reference matches the archaic genomes and
/// rarely the modern people.
fn table(d: &Path, chrom: &str, positions: impl IntoIterator<Item = i64>) {
    let mut text = String::from(TABLE_HEAD);
    for p in positions {
        text += &format!("{chrom}\t{p}\tA\tG\t0\t0\t0\t0.95\t0.95\t0.95\t0.95\t0.95\t0.95\n");
    }
    let path = d.join(format!("table/{chrom}.tsv.gz"));
    let mut gz = GzEncoder::new(std::fs::File::create(path).unwrap(), Compression::new(1));
    gz.write_all(text.as_bytes()).unwrap();
    gz.finish().unwrap();
}

/// A reference block over `from`-`to` with these FILTER and FORMAT values.
fn block(chrom: &str, from: i64, to: i64, filter: &str, format: &str, values: &str) -> String {
    format!("{chrom}\t{from}\t.\tA\t<*>\t0\t{filter}\tEND={to}\t{format}\t{values}\n")
}

/// `og archaic genotypes` with these VCF records and, when given, these gVCF
/// records: the exit code and stderr.
fn run(d: &Path, vcf: &str, gvcf: Option<&str>) -> (u8, String) {
    std::fs::write(d.join("s.vcf"), format!("{VCF_HEAD}{vcf}")).unwrap();
    let mut args: Vec<String> = [
        "genotypes",
        "--table-dir",
        d.join("table").to_str().unwrap(),
        "--vcf",
        d.join("s.vcf").to_str().unwrap(),
        "--pseudo",
        "2",
        "--out-dir",
        d.join("gt").to_str().unwrap(),
    ]
    .map(str::to_owned)
    .to_vec();
    if let Some(g) = gvcf {
        std::fs::write(d.join("s.g.vcf"), format!("{VCF_HEAD}{g}")).unwrap();
        args.extend([
            "--gvcf".to_owned(),
            d.join("s.g.vcf").to_str().unwrap().to_owned(),
        ]);
    }
    let mut err = Vec::new();
    let code = cli(&args, &mut Vec::new(), &mut err);
    (code, String::from_utf8(err).unwrap())
}

/// Your code (the SAMPLE column) at each position of a chromosome's
/// genotype file.
fn codes(d: &Path, chrom: &str) -> BTreeMap<i64, String> {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(d.join(format!("gt/{chrom}.gt.gz"))).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[1].parse().unwrap(), c[7].to_owned())
        })
        .collect()
}

fn spans(d: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(d.join("gt/spans.tsv"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

/// The issue's case: 201 table sites from 100,000 to 300,000 with no
/// record and no reference block over them, and 51 more covered by one.
/// The uncovered sites go to IBDmix as unknown (9), so they cannot add up
/// to a segment, and the span runs over the covered ones only.
#[test]
fn a_gap_in_coverage_is_unknown_not_the_reference() {
    let d = dir("gap");
    table(
        &d,
        "chr22",
        (100_000..=300_000)
            .step_by(1000)
            .chain((301_000..=351_000).step_by(1000)),
    );
    let vcf = "chr22\t50000\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n";
    let gvcf = block("chr22", 1, 99_999, ".", "GT:GQ:MIN_DP", "0/0:50:30")
        + &block("chr22", 300_500, 360_000, ".", "GT:GQ:MIN_DP", "0/0:50:30");
    let (code, err) = run(&d, vcf, Some(&gvcf));
    assert_eq!(code, 0, "{err}");
    let c = codes(&d, "chr22");
    let gap: Vec<&String> = c.range(100_000..=300_000).map(|(_, v)| v).collect();
    assert_eq!(gap.len(), 201);
    assert!(
        gap.iter().all(|v| *v == "9"),
        "the gap must be unknown: {gap:?}"
    );
    assert!(
        c.range(301_000..).all(|(_, v)| v == "0"),
        "a covered site is the reference"
    );
    assert_eq!(spans(&d), [["chr22", "301000", "351000", "51", "0", "0"]]);
    assert!(err.contains("51 sites read, 201 not read"), "{err}");
    let _ = std::fs::remove_dir_all(&d);
}

/// Without a gVCF nothing says a site with no record was read: it is
/// unknown, and a VCF with no call at any table site has nothing to analyse.
#[test]
fn without_a_gvcf_a_site_with_no_record_is_unknown() {
    let d = dir("no-gvcf");
    table(&d, "chr22", (100_000..=300_000).step_by(1000));
    let (code, err) = run(&d, "chr22\t50000\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n", None);
    assert_eq!(code, 1, "{err}");
    assert!(err.contains("no site read"), "{err}");
    assert!(!d.join("gt/chr22.gt.gz").exists());
    // A VCF with a record at every site (an all-sites VCF) needs no gVCF.
    let all: String = (100_000..=300_000)
        .step_by(1000)
        .map(|p| format!("chr22\t{p}\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\n"))
        .collect();
    let (code, err) = run(&d, &all, None);
    assert_eq!(code, 0, "{err}");
    assert!(codes(&d, "chr22").values().all(|v| v == "0"));
    let _ = std::fs::remove_dir_all(&d);
}

/// Which gVCF records count as a confident reference call, and that a VCF
/// record always wins over the gVCF.
#[test]
fn only_a_confident_reference_block_reads_as_the_reference() {
    let d = dir("confident");
    table(&d, "chr22", (1..=12).map(|i| i * 1000));
    let gvcf = [
        block("chr22", 1000, 1000, ".", "GT:GQ:MIN_DP", "0/0:20:10"), // reference
        block("chr22", 2000, 2000, "PASS", "GT", "0|0"),              // reference, phased, no GQ
        block("chr22", 3000, 3000, ".", "GT:GQ:MIN_DP", "0/0:0:10"),  // GQ 0
        block("chr22", 4000, 4000, ".", "GT:GQ:MIN_DP", "0/0:30:0"),  // no reads
        block("chr22", 5000, 5000, ".", "GT:GQ:DP", "0/0:30:0"),      // no reads (DP)
        block("chr22", 6000, 6000, "RefCall", "GT:GQ", "0/0:30"),     // filtered
        block("chr22", 7000, 7000, ".", "GT:GQ", "./.:0"),            // no call
        "chr22\t8000\t.\tA\tG,<*>\t30\tPASS\t.\tGT:GQ\t0/1:30\n".to_owned(), // a variant
        block("chr22", 9000, 9000, ".", "GT:GQ:MIN_DP", "0:20:10"),   // haploid reference
        // 10000 is in no record; 11000 is a VCF call outside every block.
        block("chr22", 12000, 12000, ".", "GT:GQ:MIN_DP", "0/0:.:."), // missing values
    ]
    .concat();
    let vcf = "chr22\t11000\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\n";
    let (code, err) = run(&d, vcf, Some(&gvcf));
    assert_eq!(code, 0, "{err}");
    let got: Vec<(i64, String)> = codes(&d, "chr22").into_iter().collect();
    let want: Vec<(i64, String)> = [
        (1000, "0"),
        (2000, "0"),
        (3000, "9"),
        (4000, "9"),
        (5000, "9"),
        (6000, "9"),
        (7000, "9"),
        (8000, "9"),
        (9000, "0"),
        (10000, "9"),
        (11000, "2"),
        (12000, "0"),
    ]
    .map(|(p, c)| (p, c.to_owned()))
    .to_vec();
    assert_eq!(got, want);
    assert_eq!(spans(&d), [["chr22", "1000", "12000", "5", "1", "0"]]);
    let _ = std::fs::remove_dir_all(&d);
}

/// A long block covers the sites up to its END, also past a shorter record
/// that starts inside it.
#[test]
fn a_site_is_covered_up_to_the_furthest_end() {
    let d = dir("overlap");
    table(&d, "chr22", [1000, 1500, 2500, 3500]);
    let gvcf = block("chr22", 900, 3000, ".", "GT:GQ", "0/0:40")
        + &block("chr22", 1200, 1300, ".", "GT:GQ", "./.:0");
    let (code, err) = run(
        &d,
        "chr22\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
        Some(&gvcf),
    );
    assert_eq!(code, 0, "{err}");
    let c = codes(&d, "chr22");
    assert_eq!(c[&1500], "0");
    assert_eq!(c[&2500], "0");
    assert_eq!(c[&3500], "9");
    let _ = std::fs::remove_dir_all(&d);
}

/// The gVCF is read alongside the VCF: the chromosomes the table lacks are
/// read past, and a gVCF in another chromosome order, or unsorted, stops
/// the step (exit 2) instead of reading as no reference calls.
#[test]
fn the_gvcf_is_read_in_the_vcf_order() {
    let d = dir("order");
    table(&d, "chr21", [1000, 2000]);
    table(&d, "chr22", [1000, 2000]);
    let vcf =
        "chr21\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\nchr22\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
    let chr21 = block("chr21", 1, 5000, ".", "GT:GQ", "0/0:40");
    let chr22 = block("chr22", 1, 5000, ".", "GT:GQ", "0/0:40");
    let other = block("chrM", 1, 100, ".", "GT:GQ", "0/0:40");
    let (code, err) = run(
        &d,
        vcf,
        Some(&format!(
            "{other}{chr21}{other_x}{chr22}",
            other_x = block("chrX", 1, 9, ".", "GT", "0/0")
        )),
    );
    assert_eq!(code, 0, "{err}");
    assert_eq!(codes(&d, "chr21")[&2000], "0");
    assert_eq!(codes(&d, "chr22")[&2000], "0");

    let (code, err) = run(&d, vcf, Some(&format!("{chr22}{chr21}")));
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("one chromosome order"), "{err}");

    let back = format!(
        "{chr21}{chr22}{}",
        block("chr21", 9000, 9100, ".", "GT", "0/0")
    );
    let (code, err) = run(&d, vcf, Some(&back));
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("comes back"), "{err}");

    let unsorted = format!(
        "{}{chr21}{chr22}",
        block("chr21", 3000, 3100, ".", "GT", "0/0")
    );
    let (code, err) = run(&d, vcf, Some(&unsorted));
    assert_eq!(code, 2, "{err}");
    assert!(err.contains("is not sorted"), "{err}");
    let _ = std::fs::remove_dir_all(&d);
}

/// A gVCF of two samples, or one without a header line, is refused.
#[test]
fn a_gvcf_of_one_sample_only() {
    let d = dir("samples");
    table(&d, "chr22", [1000]);
    let vcf = "chr22\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
    std::fs::write(d.join("s.vcf"), format!("{VCF_HEAD}{vcf}")).unwrap();
    let two = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n";
    for (gvcf, want) in [
        (two, "holds 2 samples"),
        (
            "chr22\t1\t.\tA\t<*>\t0\t.\tEND=9\tGT\t0/0\n",
            "no #CHROM header",
        ),
    ] {
        std::fs::write(d.join("x.g.vcf"), gvcf).unwrap();
        let args: Vec<String> = [
            "genotypes",
            "--table-dir",
            d.join("table").to_str().unwrap(),
            "--vcf",
            d.join("s.vcf").to_str().unwrap(),
            "--gvcf",
            d.join("x.g.vcf").to_str().unwrap(),
            "--out-dir",
            d.join("gt").to_str().unwrap(),
        ]
        .map(str::to_owned)
        .to_vec();
        let mut err = Vec::new();
        assert_eq!(cli(&args, &mut Vec::new(), &mut err), 1);
        let err = String::from_utf8(err).unwrap();
        assert!(err.contains(want), "{err}");
    }
    let _ = std::fs::remove_dir_all(&d);
}
