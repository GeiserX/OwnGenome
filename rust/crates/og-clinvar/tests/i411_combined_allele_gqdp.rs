//! The call quality of a combined allele (`og clinvar join-adjacent`), as
//! `og clinvar hits` then judges it against step 06's call floor.
//!
//! A combined record is only as sure as its weakest part. Before, it carried
//! the GT (and PS) only, so `og clinvar hits` and the summary found no GQ or
//! DP on it and never marked it low confidence, even when both parts were
//! GQ 3, DP 2 calls, and it counted in the headline. Now it carries, for
//! each of GQ and DP, the lower of its parts' values, and `0` when a part
//! lacks a value the VCF measures (a `##FORMAT` line with `Type=Integer`, or
//! the other part has it while no line defines it). A metric the VCF does
//! not measure at all stays off the record, as it is off every single call
//! of that VCF, and so does one the header defines with another type (the
//! `Type=Float` GQ FreeBayes writes) when a part has no whole number for it.
//!
//! Each check builds the combined record of two neighbouring calls on a
//! small reference, puts ClinVar's columns on it as step 06 does, and reads
//! the row `og clinvar hits` writes: its DP, GQ and LOW_CONFIDENCE columns.

use std::path::PathBuf;

use og_clinvar::hits;
use og_clinvar::join::{self, Mode, WINDOW};
use og_intervals::fasta::Indexed;

/// GRCh38 chr22 from 30000200, as contig t from position 1.
const SEQ: &str = "GCCTGTAATCCCAGCACTTTGGGAGGCTGAG";

const MEASURED: &str = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n";

/// The FORMAT lines FreeBayes writes: GQ is a Float, and with no `-=` no
/// record carries it.
const FREEBAYES: &str = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n";

const GT_ONLY: &str = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

const CV: &str = "GENEINFO=SYNTH:1;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts";

/// A scratch reference of its own, removed when dropped.
struct Ref(PathBuf);

impl Ref {
    fn new(tag: &str) -> Ref {
        let dir =
            std::env::temp_dir().join(format!("og-clinvar-i411-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("ref.fa"), format!(">t\n{SEQ}\n")).unwrap();
        std::fs::write(
            dir.join("ref.fa.fai"),
            format!("t\t{}\t3\t{}\t{}\n", SEQ.len(), SEQ.len(), SEQ.len() + 1),
        )
        .unwrap();
        Ref(dir)
    }
}

impl Drop for Ref {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The two parts: T>A at 4 and T>C at 9, each with its FORMAT and sample.
fn vcf(header: &str, first: (&str, &str), second: (&str, &str)) -> String {
    format!(
        "##fileformat=VCFv4.2\n{header}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
         t\t4\t.\tT\tA\t50\tPASS\t.\t{}\t{}\n\
         t\t9\t.\tT\tC\t50\tPASS\t.\t{}\t{}\n",
        first.0, first.1, second.0, second.1
    )
}

/// The one combined record build writes, split into its columns.
fn combined(tag: &str, text: &str) -> Vec<String> {
    let r = Ref::new(tag);
    let mut fa = Indexed::open(r.0.join("ref.fa")).unwrap();
    let mut out = Vec::new();
    join::join(Mode::Build(&mut fa), WINDOW, text.as_bytes(), &mut out).unwrap();
    let built = String::from_utf8(out).unwrap();
    let recs: Vec<&str> = built.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(recs.len(), 1, "one combined record: {built}");
    let mut prep = Vec::new();
    join::join(Mode::Prepare, WINDOW, text.as_bytes(), &mut prep).unwrap();
    let prep = String::from_utf8(prep).unwrap();
    let cols: Vec<String> = recs[0].split('\t').map(str::to_owned).collect();
    let pcols: Vec<&str> = prep
        .lines()
        .find(|l| !l.starts_with('#'))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(
        (pcols[8], pcols[9]),
        (cols[8].as_str(), cols[9].as_str()),
        "prepare writes the FORMAT and sample columns build writes"
    );
    cols
}

/// The DP, GQ and LOW_CONFIDENCE columns of the hits row of the combined
/// record, with ClinVar's columns on it, under the floor line `floor` (none
/// for the default GQ 30, DP 10).
fn judged(cols: &[String], floor: &str) -> [String; 3] {
    let mut c = cols.to_vec();
    c[7] = format!("{CV};{}", c[7]);
    let file = format!("{floor}#CHROM\n{}\n", c.join("\t"));
    let mut out = Vec::new();
    hits::hits(file.as_bytes(), &mut out).unwrap();
    let out = String::from_utf8(out).unwrap();
    let row: Vec<&str> = out.trim_end().split('\t').collect();
    assert_eq!(row.len(), 13, "{out}");
    assert!(
        row[9].starts_with("combined allele of t:4 T>A + t:9 T>C"),
        "{out}"
    );
    [row[10].to_owned(), row[11].to_owned(), row[12].to_owned()]
}

fn s(v: [&str; 3]) -> [String; 3] {
    v.map(str::to_owned)
}

#[test]
fn two_weak_parts_give_a_low_confidence_combined_allele() {
    let rec = combined(
        "weak",
        &vcf(MEASURED, ("GT:GQ:DP", "0/1:3:2"), ("GT:GQ:DP", "0/1:3:2")),
    );
    assert_eq!((rec[8].as_str(), rec[9].as_str()), ("GT:GQ:DP", "0/1:3:2"));
    assert_eq!(judged(&rec, ""), s(["2", "3", "yes"]));
}

#[test]
fn the_combined_allele_takes_the_weaker_part_for_each_metric() {
    // The first part has the lower GQ, the second the lower DP; the keys can
    // come in any order.
    let rec = combined(
        "mixed",
        &vcf(
            MEASURED,
            ("GT:GQ:DP", "1/0:12:40"),
            ("DP:GQ:GT", "9:50:0/1"),
        ),
    );
    assert_eq!((rec[8].as_str(), rec[9].as_str()), ("GT:GQ:DP", "0/1:12:9"));
    assert_eq!(judged(&rec, ""), s(["9", "12", "yes"]));
    // One weak part beside a strong one is enough.
    let rec = combined(
        "one-weak",
        &vcf(MEASURED, ("GT:GQ:DP", "0/1:60:45"), ("GT:GQ:DP", "0/1:3:2")),
    );
    assert_eq!(judged(&rec, ""), s(["2", "3", "yes"]));
}

#[test]
fn a_part_without_gq_or_dp_gives_low_confidence() {
    // The second part has no GQ or DP although the VCF measures both.
    let rec = combined(
        "unmeasured",
        &vcf(MEASURED, ("GT:GQ:DP", "0/1:60:45"), ("GT", "0/1")),
    );
    assert_eq!((rec[8].as_str(), rec[9].as_str()), ("GT:GQ:DP", "0/1:0:0"));
    assert_eq!(judged(&rec, ""), s(["0", "0", "yes"]));
    // A value that is no whole number is no value either.
    let rec = combined(
        "dot",
        &vcf(
            MEASURED,
            ("GT:GQ:DP", "0/1:.:45"),
            ("GT:GQ:DP", "0/1:60:45"),
        ),
    );
    assert_eq!(judged(&rec, ""), s(["45", "0", "yes"]));
    // Neither part has them, and the header defines them.
    let rec = combined("neither", &vcf(MEASURED, ("GT", "0/1"), ("GT", "0/1")));
    assert_eq!(judged(&rec, ""), s(["0", "0", "yes"]));
    // One part has a GQ the header does not define: GQ is measured.
    let rec = combined(
        "undeclared",
        &vcf(GT_ONLY, ("GT:GQ", "0/1:60"), ("GT", "0/1")),
    );
    assert_eq!((rec[8].as_str(), rec[9].as_str()), ("GT:GQ", "0/1:0"));
    assert_eq!(judged(&rec, ""), s([".", "0", "yes"]));
}

#[test]
fn two_strong_parts_stay_confident() {
    let rec = combined(
        "strong",
        &vcf(
            MEASURED,
            ("GT:GQ:DP", "0/1:50:40"),
            ("GT:GQ:DP", "0/1:45:35"),
        ),
    );
    assert_eq!(
        (rec[8].as_str(), rec[9].as_str()),
        ("GT:GQ:DP", "0/1:45:35")
    );
    assert_eq!(judged(&rec, ""), s(["35", "45", "no"]));
    // A run's own floor applies: GQ 3, DP 2 parts pass a floor of GQ 3, DP 2.
    let rec = combined(
        "floor",
        &vcf(MEASURED, ("GT:GQ:DP", "0/1:3:2"), ("GT:GQ:DP", "0/1:4:2")),
    );
    assert_eq!(
        judged(&rec, "##PGP_CALL_FLOOR=GQ=3,DP=2\n"),
        s(["2", "3", "no"])
    );
}

#[test]
fn a_phased_pair_keeps_its_phase_set_beside_gq_and_dp() {
    let rec = combined(
        "phased",
        &vcf(
            MEASURED,
            ("GT:PS:GQ:DP", "0|1:7:40:30"),
            ("GT:GQ:DP:PS", "0|1:20:31:7"),
        ),
    );
    assert_eq!(
        (rec[8].as_str(), rec[9].as_str()),
        ("GT:PS:GQ:DP", "0|1:7:20:30")
    );
    assert_eq!(judged(&rec, ""), s(["30", "20", "yes"]));
}

#[test]
fn a_vcf_that_measures_neither_writes_neither() {
    // A VCF with GT only (a converted chip file, say): no single call of it
    // is ever marked, and the combined record stays as it was.
    let rec = combined("gt-only", &vcf(GT_ONLY, ("GT", "0/1"), ("GT", "1/0")));
    assert_eq!((rec[8].as_str(), rec[9].as_str()), ("GT", "0/1"));
    assert_eq!(judged(&rec, ""), s([".", ".", "no"]));
}

#[test]
fn a_float_gq_is_never_read_as_zero() {
    // FreeBayes defines GQ as a Float and, run without `-=`, writes it on no
    // record. Two strong parts stay confident; DP still takes the weaker.
    let rec = combined(
        "freebayes",
        &vcf(FREEBAYES, ("GT:DP", "0/1:60"), ("GT:DP", "0/1:58")),
    );
    assert_eq!((rec[8].as_str(), rec[9].as_str()), ("GT:DP", "0/1:58"));
    assert_eq!(judged(&rec, ""), s(["58", ".", "no"]));
    // With `-=` the values are not whole numbers, which hits never reads.
    let rec = combined(
        "freebayes-gq",
        &vcf(
            FREEBAYES,
            ("GT:GQ:DP", "0/1:160.046:60"),
            ("GT:GQ:DP", "0/1:2.5:8"),
        ),
    );
    assert_eq!((rec[8].as_str(), rec[9].as_str()), ("GT:DP", "0/1:8"));
    assert_eq!(judged(&rec, ""), s(["8", ".", "yes"]));
    // A Description that mentions Type=Integer does not make GQ whole.
    let said = FREEBAYES.replace("\"Genotype Quality\"", "\"Previously Type=Integer\"");
    let rec = combined(
        "freebayes-said",
        &vcf(&said, ("GT:DP", "0/1:60"), ("GT:DP", "0/1:58")),
    );
    assert_eq!((rec[8].as_str(), rec[9].as_str()), ("GT:DP", "0/1:58"));
}
