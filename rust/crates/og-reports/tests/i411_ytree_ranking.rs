//! The Y placement ranks lineages by their evidence, says what it set aside,
//! and the private-SNP scan reads one man and checks the allele (step 37).
//!
//! The tree is the one of `y_tree_place.rs`:
//!
//! ```text
//! (root) -> A (a1, a2/a2alt) -> B (b1, b2, b3) -> C (c1, c2)
//!                            -> X (x1)        -> Y (y1) -> Z (z1) -> W (w1)
//!                            -> P (p1)        -> Q (q1, q1b; one site)
//! ```
//!
//! 1. compete: x1 ancestral and outvoted by y1 and z1, against a clean P-Q.
//!    Q has more support from A (3 against 1), so Q wins; the old walk took
//!    Z because it is deeper. Z is reported as the conflict.
//! 2. the skeptic's three cases: s1, a clean B (3) against the outvoted Z
//!    (1): B; s2, the outvoted lineage with more support (W, 2) against a
//!    single derived P (1): W stays; s3, equal support on two lineages: the
//!    deeper one wins.
//! 3. on one lineage the outvote still takes the walk to the deepest branch.
//! 4. `Informative_callable` and `Path_breadth` count the readable sites.
//! 5. the current tree's pileup flags are written out (`-A -Q20 -q20`), in
//!    the step script, the Nextflow module and the e2e mirror: a fused
//!    `-AQ20q1` makes samtools read `20q1` as the `-Q` value, so the
//!    mapping-quality filter is lost and MQ0 reads count.
//! 6. private SNPs: another base at a known position stays private, the
//!    known change (either orientation) does not; a VCF with two samples is
//!    refused unless `--sample` names one, and then that column is read.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

const TREE: &str = r#"{"id": "", "snps": "", "children": [
 {"id": "A", "snps": "a1, a2/a2alt", "children": [
  {"id": "B", "snps": "b1, b2, b3", "children": [
   {"id": "C", "snps": "c1, c2", "children": []}]},
  {"id": "X", "snps": "x1", "children": [
   {"id": "Y", "snps": "y1", "children": [
    {"id": "Z", "snps": "z1", "children": [
     {"id": "W", "snps": "w1", "children": []}]}]}]},
  {"id": "P", "snps": "p1", "children": [
   {"id": "Q", "snps": "q1, q1b", "children": []}]}]}]}"#;

const SNPS: [(&str, u32, char, char); 14] = [
    ("a1", 3000100, 'C', 'T'),
    ("a2alt", 3000200, 'G', 'A'),
    ("b1", 3000300, 'A', 'G'),
    ("b2", 3000400, 'T', 'C'),
    ("b3", 3000500, 'C', 'G'),
    ("c1", 3000600, 'G', 'T'),
    ("c2", 3000700, 'A', 'T'),
    ("x1", 3000800, 'A', 'C'),
    ("y1", 3000900, 'T', 'A'),
    ("z1", 3001000, 'G', 'C'),
    ("w1", 3001100, 'C', 'A'),
    ("p1", 3001200, 'T', 'G'),
    ("q1", 3001300, 'A', 'T'),
    ("q1b", 3001300, 'A', 'T'),
];

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-i411-ytree-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn write_tree(d: &Path) {
    let rows: String = SNPS
        .iter()
        .map(|(n, p, a, de)| format!("chry\t{n}\tR-{n}\t{p}\t{a}->{de}\t{a}\t{de}\n"))
        .collect();
    fs::write(d.join("current_tree.json"), TREE).unwrap();
    fs::write(d.join("positions_ftdna_hg38.txt"), rows).unwrap();
    fs::write(d.join("positions_yfull_hg38.txt"), "").unwrap();
    fs::write(d.join("version.txt"), "tree\tt\npositions\tp\n").unwrap();
}

fn y(verb: &str, args: &[String]) -> (u8, String, String) {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_reports::y::cli(Some(verb), args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

fn strs(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| (*s).to_owned()).collect()
}

/// The result row by column name, and the step line, for these SNPs read
/// derived ('D') or ancestral ('A') at 12 reads each.
fn place(d: &Path, tag: &str, states: &[(&str, char)]) -> (HashMap<String, String>, String) {
    let mut text = String::new();
    for &(name, st) in states {
        let &(_, pos, anc, der) = SNPS.iter().find(|r| r.0 == name).unwrap();
        let base = if st == 'D' { der } else { anc };
        text += &format!(
            "chrY\t{pos}\tN\t12\t{}\t{}\n",
            base.to_string().repeat(12),
            "I".repeat(12)
        );
    }
    let pileup = d.join(format!("{tag}.pileup"));
    fs::write(&pileup, text).unwrap();
    let prefix = d.join(tag);
    let (code, out, err) = y(
        "tree-place",
        &strs(&[
            "--tree-dir",
            d.to_str().unwrap(),
            "--pileup",
            pileup.to_str().unwrap(),
            "--sample",
            "S",
            "--out-prefix",
            prefix.to_str().unwrap(),
        ]),
    );
    assert_eq!(code, 0, "{tag}: {err}");
    let t = fs::read_to_string(format!("{}_y_tree.txt", prefix.display())).unwrap();
    let mut lines = t.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let row = lines.next().unwrap().split('\t');
    let row = head
        .iter()
        .zip(row)
        .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
        .collect();
    (row, out)
}

fn cols<'a>(r: &'a HashMap<String, String>, keys: &[&str]) -> Vec<&'a str> {
    keys.iter().map(|k| r[*k].as_str()).collect()
}

const CONFLICT: [&str; 4] = ["Hg", "Conflict", "Conflict_support", "Terminal_support"];

/// Red first: the depth-first walk returned Z here.
#[test]
fn red_first_a_clean_lineage_beats_a_deeper_outvoted_one() {
    let d = scratch("compete");
    write_tree(&d);
    let (r, line) = place(
        &d,
        "compete",
        &[
            ("a1", 'D'),
            ("a2alt", 'D'),
            ("x1", 'A'),
            ("y1", 'D'),
            ("z1", 'D'),
            ("p1", 'D'),
            ("q1", 'D'),
        ],
    );
    assert_eq!(r["Hg"], "Q", "the clean lineage with more support wins");
    assert_eq!(cols(&r, &CONFLICT), ["Q", "Z", "1", "3"], "{r:?}");
    assert_eq!(r["Outvoted_ancestral"], ".", "x1 is not on Q's path");
    assert!(
        line.contains("another lineage has support too: Z (1 against 3"),
        "the step line names the conflict: {line}"
    );
    fs::remove_dir_all(d).unwrap();
}

#[test]
fn the_skeptic_cases() {
    let d = scratch("skeptic");
    write_tree(&d);
    let a = [("a1", 'D'), ("a2alt", 'D')];
    let with = |more: &[(&'static str, char)]| -> Vec<(&'static str, char)> {
        a.iter().chain(more).copied().collect()
    };
    // s1: the old outvote case with B clean: B (3) beats Z (1).
    let (r, _) = place(
        &d,
        "s1",
        &with(&[
            ("b1", 'D'),
            ("b2", 'D'),
            ("b3", 'D'),
            ("x1", 'A'),
            ("y1", 'D'),
            ("z1", 'D'),
        ]),
    );
    assert_eq!(r["Hg"], "B", "s1");
    assert_eq!(cols(&r, &CONFLICT), ["B", "Z", "1", "3"], "s1 {r:?}");
    // s2: the outvoted lineage has more support (-1 + 3) than one derived p1.
    let (r, line) = place(
        &d,
        "s2",
        &with(&[
            ("x1", 'A'),
            ("y1", 'D'),
            ("z1", 'D'),
            ("w1", 'D'),
            ("p1", 'D'),
        ]),
    );
    assert_eq!(cols(&r, &CONFLICT), ["W", "P", "1", "2"], "s2 {r:?}");
    assert_eq!(r["Outvoted_ancestral"], "x1", "s2 {r:?}");
    assert!(
        line.contains("outvoted") && line.contains("P (1 against 2"),
        "{line}"
    );
    // s3: equal support on two lineages: the deeper one.
    let (r, _) = place(
        &d,
        "s3",
        &with(&[("b1", 'D'), ("b2", 'D'), ("x1", 'D'), ("y1", 'D')]),
    );
    assert_eq!(cols(&r, &CONFLICT), ["Y", "B", "2", "2"], "s3 {r:?}");
    // A lineage the outvote does not carry (one site below p1) is no conflict.
    let (r, line) = place(
        &d,
        "lone",
        &with(&[("b1", 'D'), ("b2", 'D'), ("p1", 'A'), ("q1", 'D')]),
    );
    assert_eq!(cols(&r, &CONFLICT), ["B", ".", ".", "."], "lone {r:?}");
    assert!(!line.contains("another lineage"), "{line}");
    fs::remove_dir_all(d).unwrap();
}

#[test]
fn one_lineage_still_goes_deepest() {
    let d = scratch("linear");
    write_tree(&d);
    let (r, _) = place(
        &d,
        "linear",
        &[
            ("a1", 'D'),
            ("a2alt", 'D'),
            ("x1", 'A'),
            ("y1", 'D'),
            ("z1", 'D'),
            ("w1", 'D'),
        ],
    );
    assert_eq!(
        cols(&r, &["Hg", "Level", "Outvoted_ancestral", "Conflict"]),
        ["W", "5", "x1", "."]
    );
    fs::remove_dir_all(d).unwrap();
}

#[test]
fn callable_and_breadth() {
    let d = scratch("breadth");
    write_tree(&d);
    // b3 has no read: B's path has 5 SNPs on the tables, 4 of them readable.
    let (r, _) = place(
        &d,
        "breadth",
        &[
            ("a1", 'D'),
            ("a2alt", 'D'),
            ("b1", 'D'),
            ("b2", 'D'),
            ("x1", 'A'),
        ],
    );
    assert_eq!(
        cols(
            &r,
            &["Hg", "Informative_callable", "Path_breadth", "SNPs_called"]
        ),
        ["B", "5", "0.800", "5"]
    );
    // Nothing derived: no path, no breadth.
    let (r, _) = place(&d, "none", &[("a1", 'A')]);
    assert_eq!(
        cols(
            &r,
            &["Hg", "Informative_callable", "Path_breadth", "Conflict"]
        ),
        ["NA", "1", ".", "."]
    );
    fs::remove_dir_all(d).unwrap();
}

/// A flag cluster that runs letters on after a `-Q` value (`-AQ20q1`):
/// samtools reads the rest as part of the value and drops those flags.
fn fused_after_q(line: &str) -> Option<&str> {
    line.split_whitespace().find(|t| {
        let Some(flags) = t.strip_prefix('-').filter(|f| !f.starts_with('-')) else {
            return false;
        };
        let Some(q) = flags.find('Q') else {
            return false;
        };
        let rest = &flags[q + 1..];
        let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
        digits > 0 && rest[digits..].starts_with(|c: char| c.is_ascii_alphabetic())
    })
}

#[test]
fn the_tree_pileup_flags_are_written_out() {
    assert_eq!(
        fused_after_q("samtools mpileup -l p -AQ20q1 b"),
        Some("-AQ20q1")
    );
    assert_eq!(fused_after_q("samtools mpileup -l p -A -Q20 -q20 b"), None);
    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let mut tree_lines = 0;
    for f in [
        "scripts/37-y-haplogroup.sh",
        "modules/local/y_haplogroup/main.nf",
        "tests/e2e/y-tree-current-1-placement.sh",
    ] {
        let text = fs::read_to_string(repo.join(f)).unwrap();
        for (n, line) in text.lines().enumerate() {
            let code = line.trim_start();
            if !code.contains("mpileup") || code.starts_with('#') || code.starts_with("//") {
                continue;
            }
            if code.contains("tree") {
                tree_lines += 1;
                assert_eq!(fused_after_q(code), None, "{f}:{}: {line}", n + 1);
                assert!(code.contains(" -A -Q20 -q20 "), "{f}:{}: {line}", n + 1);
            } else {
                // Yleaf's own pileup keeps Yleaf's flags byte for byte (the
                // YLEAF_IMAGE smoke row checks them against Yleaf's source).
                assert!(code.contains(" -ABQ20q1 "), "{f}:{}: {line}", n + 1);
            }
        }
    }
    assert_eq!(tree_lines, 3, "one tree pileup in each of the three files");
}

/// The private-SNP scan on a plain mask and one table row per SNP.
struct Private {
    d: PathBuf,
}

impl Private {
    fn new(name: &str) -> Private {
        let d = scratch(name);
        fs::write(
            d.join("positions_ftdna_hg38.txt"),
            "chrY\tk1\tR\t7000100\tG->A\tG\tA\nchrY\tk2\tR\t7000200\tC->T\tC\tT\n",
        )
        .unwrap();
        fs::write(
            d.join("positions_yfull_hg38.txt"),
            "chrY\tk3\tR\t7000300\tG->A\tG\tA\nchrY\tk4\tR\t7000400\t.\t.\t.\n",
        )
        .unwrap();
        let seq: String = "P".repeat(7_001_000);
        let mut mask = String::from(">chrY\n");
        for chunk in seq.as_bytes().chunks(60) {
            mask += std::str::from_utf8(chunk).unwrap();
            mask.push('\n');
        }
        fs::write(d.join("m.fa"), mask).unwrap();
        Private { d }
    }

    fn run(&self, vcf: &str, sample: &str) -> (u8, String, Vec<String>) {
        let v = self.d.join("v.vcf");
        fs::write(&v, vcf).unwrap();
        let prefix = self.d.join("o");
        let (code, _, err) = y(
            "private-snps",
            &strs(&[
                "--vcf",
                v.to_str().unwrap(),
                "--tree-dir",
                self.d.to_str().unwrap(),
                "--callable-mask",
                self.d.join("m.fa").to_str().unwrap(),
                "--sample",
                sample,
                "--out-prefix",
                prefix.to_str().unwrap(),
            ]),
        );
        let rows = fs::read_to_string(format!("{}_y_private_snps.tsv", prefix.display()))
            .map(|t| {
                t.lines()
                    .filter(|l| l.starts_with("70"))
                    .map(|l| l.split('\t').take(3).collect::<Vec<_>>().join(" "))
                    .collect()
            })
            .unwrap_or_default();
        let _ = fs::remove_file(format!("{}_y_private_snps.tsv", prefix.display()));
        (code, err, rows)
    }
}

impl Drop for Private {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.d);
    }
}

const HEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

fn rec(pos: u32, r: &str, a: &str, samples: &[&str]) -> String {
    format!(
        "chrY\t{pos}\t.\t{r}\t{a}\t50\tPASS\t.\tGT:AD:DP:GQ\t{}\n",
        samples.join("\t")
    )
}

const CALL: &str = "1:0,30:30:50";
const REF: &str = "0:30,0:30:50";

/// Red first: the position alone dropped 7000100 G>T and 7000200 C>G.
#[test]
fn red_first_another_base_at_a_known_position_stays_private() {
    let p = Private::new("allele");
    let vcf = format!(
        "{HEAD}\tS\n{}{}{}{}{}{}",
        rec(7000100, "G", "A", &[CALL]), // the known change: dropped
        rec(7000100, "G", "T", &[CALL]), // another base there: private
        rec(7000200, "T", "C", &[CALL]), // the known change, reference derived: dropped
        rec(7000200, "C", "G", &[CALL]), // private
        rec(7000300, "g", "a", &[CALL]), // the second table, lower case: dropped
        rec(7000400, "C", "T", &[CALL]), // a row with no usable bases: dropped
    );
    let (code, err, rows) = p.run(&vcf, "S");
    assert_eq!(code, 0, "{err}");
    assert_eq!(rows, ["7000100 G T", "7000200 C G"]);
}

#[test]
fn red_first_a_two_sample_vcf_needs_a_name() {
    let p = Private::new("samples");
    let body = format!(
        "{}{}",
        rec(7000500, "C", "T", &[REF, CALL]),
        rec(7000600, "C", "T", &[CALL, REF])
    );
    let (code, err, _) = p.run(&format!("{HEAD}\tdad\tson\n{body}"), "S");
    assert_eq!(code, 1, "refused: {err}");
    assert!(
        err.contains("2 samples (dad, son)") && err.contains("--sample"),
        "{err}"
    );
    let (code, err, rows) = p.run(&format!("{HEAD}\tdad\tson\n{body}"), "son");
    assert_eq!((code, rows), (0, vec!["7000500 C T".to_owned()]), "{err}");
    let (code, _, rows) = p.run(&format!("{HEAD}\tdad\tson\n{body}"), "dad");
    assert_eq!((code, rows), (0, vec!["7000600 C T".to_owned()]));
    // No #CHROM line to name the columns: refused.
    let (code, err, _) = p.run(&format!("##fileformat=VCFv4.2\n{body}"), "son");
    assert!(
        code == 1 && err.contains("more than one sample column"),
        "{err}"
    );
    // One sample: read whatever its name.
    let one = rec(7000700, "C", "T", &[CALL]);
    let (code, _, rows) = p.run(&format!("{HEAD}\tsomeone\n{one}"), "S");
    assert_eq!((code, rows), (0, vec!["7000700 C T".to_owned()]));
}
