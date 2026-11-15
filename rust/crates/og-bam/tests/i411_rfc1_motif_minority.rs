//! Synthetic BAM regressions for pathogenic and unknown motifs under a benign majority.
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use rust_htslib::bam::{self, header::HeaderRecord, record::{Cigar, CigarString}};

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn count(counts: &[(&str, usize)], dominant: bool) -> (String, String) {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let dir = Scratch(std::env::temp_dir().join(format!(
        "og-i411-motif-{}-{}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)
    )));
    fs::create_dir_all(&dir.0).unwrap();
    let bam_path = dir.0.join("reads.bam");
    let vcf = dir.0.join("calls.vcf");
    let out = dir.0.join("motif.tsv");
    let catalog = dir.0.join("catalog.json");
    let locus = if dominant { "FAMEX" } else { "RFC1" };
    fs::write(&vcf, format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr4\t101\t.\tA\t<STR>\t.\tPASS\tEND=155;REPID={locus};STR_STATUS=full_mutation\n"
    )).unwrap();
    fs::write(&catalog, r#"[{"LocusId":"FAMEX","ReferenceRegion":"chr4:100-155","Disease":"FAME","InheritanceMode":"AD","BenignMotifs":["TTTTA"],"PathogenicMotifs":["TTTCA"]}]"#).unwrap();
    let mut header = bam::Header::new();
    header.push_record(HeaderRecord::new(b"HD").push_tag(b"SO", "coordinate"));
    header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", "chr4").push_tag(b"LN", 1000));
    {
        let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        for (motif, n) in counts {
            for i in 0..*n {
                let mut rec = bam::Record::new();
                // Both strands must keep the same evidence. 150 bases exceed the reference run.
                let reverse = i % 2 == 1;
                let seq: Vec<u8> = if reverse {
                    motif.bytes().rev().map(|c| match c {
                        b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', _ => unreachable!()
                    }).collect::<Vec<_>>().repeat(30)
                } else { motif.as_bytes().repeat(30) };
                rec.set(format!("{motif}-{i}").as_bytes(), Some(&CigarString(vec![Cigar::Match(150)])), &seq, &[30; 150]);
                rec.set_tid(0);
                rec.set_pos(75);
                rec.set_flags(if reverse { 16 } else { 0 });
                rec.set_mapq(60);
                writer.write(&rec).unwrap();
            }
        }
    }
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();
    let mut args = vec!["--bam".into(), bam_path.display().to_string(), "--vcf".into(), vcf.display().to_string(), "--out".into(), out.display().to_string()];
    if dominant {
        args.extend(["--catalog".into(), catalog.display().to_string(), "--flagged".into()]);
    }
    let (mut stdout, mut stderr) = (Vec::new(), Vec::new());
    assert_eq!(og_bam::motif::cli(&args, &mut stdout, &mut stderr), 0, "{}", String::from_utf8_lossy(&stderr));
    (fs::read_to_string(out).unwrap(), String::from_utf8(stdout).unwrap())
}

fn row(table: &str) -> Vec<&str> {
    table.lines().nth(1).unwrap().split('\t').collect()
}

fn uncertain(counts: &[(&str, usize)], dominant: bool, evidence: &str) {
    let (table, _) = count(counts, dominant);
    let r = row(&table);
    assert_eq!(r[9], "uncertain_motif", "{table}");
    assert!(r[10].contains(evidence), "{table}");
    assert!(r[10].contains("possible"), "{table}");
    assert!(!r[10].contains("the flag is not supported"), "{table}");
}

#[test]
fn rfc1_keeps_fifteen_percent_pathogenic_reads() {
    uncertain(&[("AAAAG", 17), ("AAGGG", 3)], false, "AAGGG in 3 of 20");
}

#[test]
fn rfc1_keeps_one_or_two_pathogenic_reads_without_calling_them_pathogenic() {
    uncertain(&[("AAAAG", 4), ("AAGGG", 2)], false, "AAGGG in 2 of 6");
    uncertain(&[("AAAAG", 19), ("AAGGG", 1)], false, "AAGGG in 1 of 20");
}

#[test]
fn dominant_insert_can_be_a_minority_and_still_needs_confirmation() {
    uncertain(&[("TTTTA", 17), ("TTTCA", 3)], true, "TTTCA in 3 of 20");
    uncertain(&[("TTTTA", 999), ("TTTCA", 1)], true, "TTTCA in 1 of 1000");
}

#[test]
fn unknown_minority_and_tie_are_not_dismissed() {
    // AGAGG is named AAGAG after strand/rotation normalization.
    uncertain(&[("AAAAG", 6), ("AGAGG", 4)], false, "AAGAG in 4 of 10");
    uncertain(&[("AAAAG", 5), ("AGAGG", 5)], false, "AAGAG in 5 of 10");
    uncertain(&[("AAAAG", 19), ("AGAGG", 1)], false, "AAGAG in 1 of 20");
}

#[test]
fn twenty_and_twenty_five_percent_are_possible_without_resolving_both_copies() {
    for counts in [[("AAAAG", 16), ("AAGGG", 4)], [("AAAAG", 15), ("AAGGG", 5)]] {
        let (table, _) = count(&counts, false);
        let r = row(&table);
        assert_eq!(r[9], "pathogenic_motif");
        assert!(r[10].contains("confirm with repeat-primed PCR"));
        assert!(r[10].contains("needs a pathogenic expansion on both copies"));
        assert!(!r[10].contains("carrier"));
    }
}

#[test]
fn pure_and_mixed_known_benign_motifs_remain_benign() {
    for counts in [vec![("AAAAG", 10)], vec![("AAAAG", 6), ("AAAGG", 4)]] {
        let (table, _) = count(&counts, false);
        assert_eq!(row(&table)[9], "benign_motif");
    }
    let (table, _) = count(&[("TTTTA", 10)], true);
    assert_eq!(row(&table)[9], "benign_motif");
}

#[test]
fn fewer_than_three_total_reads_remain_unresolved() {
    let (table, _) = count(&[("AAGGG", 2)], false);
    assert_eq!(row(&table)[9], "not_resolved");
}

#[test]
fn minority_outputs_match_the_deliberate_parity_expectations() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../tests/parity/want/rfc1");
    for (case, counts, dominant) in [
        ("pathogenic-minority", [("AAAAG", 17), ("AAGGG", 3)], false),
        ("dominant-minority", [("TTTTA", 17), ("TTTCA", 3)], true),
        ("unknown-tie", [("AAAAG", 5), ("AGAGG", 5)], false),
    ] {
        let (table, stdout) = count(&counts, dominant);
        let want = fs::read_to_string(root.join(format!("{case}.tsv"))).unwrap();
        assert_eq!(table, want);
        let r = row(&table);
        assert_eq!(stdout, format!("{} motif count: {} ({})\n", r[0], r[9], r[10]));
    }
}
