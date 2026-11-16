//! `og clinvar join-adjacent`, `og clinvar hits` and
//! `og clinvar deep-intronic` on hand-made records at real GRCh38 chr22
//! positions: what the old `tests/test_clinvar_join_adjacent.py` checked,
//! by meaning, in cargo test.
//!
//! 1. Two records within 10 bp, both het with the same unphased GT, become
//!    one REF/ALT with the reference bases between them; 0/1 and 1/0 are
//!    the same het call.
//! 2. Phased records join only in the same phase set with the same GT.
//! 3. A het next to a hom, an unphased half call, overlapping records and
//!    records 11 bp apart are not joined; 10 bp apart are; a phased half
//!    call joins a record with its ALT on the same haplotype.
//! 4. A single SNV gives no combined record.
//! 5. finish drops a combined record whose part disagrees with the
//!    reference.
//! 6. The combined records come out in POS order.
//! 7. The tag pass marks the parts of a matched combined allele, flags two
//!    hits in one gene that can be one allele, and leaves the rest alone.
//! 8. The fixture records of `tests/fixtures/vcf/combined_allele_*.vcf`,
//!    read in place: the one combined record, left-trimmed as step 06's
//!    `bcftools norm` writes it, is ClinVar's 2-star record; its parts are
//!    ClinVar's two 1-star records; the tag pass marks both and leaves the
//!    single SNV hit unchanged, four hits in all. The real step 06 with
//!    bcftools on the same files is the e2e case `clinvar-combined-allele-1`.
//! 9. A 1/2 record at a poly-A tract split as step 06 splits it (1/. for
//!    ClinVar's allele): the hit is worded "one copy, other copy not called
//!    ref", never het; a plain 0/1 hit stays het; and the always-report
//!    list gives every site a call. `tests/test_isec_columns.sh` runs the
//!    real split through step 06.
//!
//! Prepare, finish and build must agree: every join check runs prepare,
//! a fill of `PGP_JOIN_SEQ` as `bcftools +fill-from-fasta` does it, then
//! finish, and build on a FASTA, and both must write the same records.
//! Set `WINDOW` to 0 in `src/join.rs` and checks 1, 2, 3, 6, 7 and 8 fail.

use std::io::Cursor;
use std::path::{Path, PathBuf};

use og_clinvar::join::{self, Mode, WINDOW};
use og_clinvar::{hits, sites};
use og_intervals::fasta::Indexed;

/// GRCh38 chr22 (the same in the NCBI no-ALT analysis set), two short
/// stretches.
const REF: [(u64, &str); 2] = [
    (
        30000000,
        "AGACTCTGTCGCAAAAAAAAAATTAATGTAACCATTTAGCCATTACCCAGTTCCTTAAAAT",
    ),
    (30000200, "GCCTGTAATCCCAGCACTTTGGGAGGCTGAG"),
];

const HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr22,length=50818468>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// A scratch folder of its own for one test, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!("og-clinvar-join-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn ref_at(pos: u64, n: usize) -> String {
    for (start, seq) in REF {
        if start <= pos && pos as usize + n <= start as usize + seq.len() {
            let i = (pos - start) as usize;
            return seq[i..i + n].to_owned();
        }
    }
    panic!("chr22:{pos} (+{n}) is outside the test's reference stretches");
}

/// chr22 up to the last stretch, N everywhere else, 60 bases a line, with
/// its .fai.
fn write_reference(dir: &Path) -> String {
    let length = REF
        .iter()
        .map(|(s, q)| *s as usize + q.len())
        .max()
        .unwrap()
        + 50;
    let mut seq = vec![b'N'; length];
    for (start, q) in REF {
        let s = start as usize - 1;
        seq[s..s + q.len()].copy_from_slice(q.as_bytes());
    }
    let mut fa = b">chr22\n".to_vec();
    for chunk in seq.chunks(60) {
        fa.extend_from_slice(chunk);
        fa.push(b'\n');
    }
    let path = dir.join("ref.fa");
    std::fs::write(&path, fa).unwrap();
    std::fs::write(
        dir.join("ref.fa.fai"),
        format!("chr22\t{length}\t7\t60\t61\n"),
    )
    .unwrap();
    path.to_string_lossy().into_owned()
}

fn rec(pos: u64, r: &str, alt: &str, gt: &str, ps: Option<u32>) -> String {
    let (fmt, smp) = match ps {
        None => ("GT".to_owned(), gt.to_owned()),
        Some(p) => ("GT:PS".to_owned(), format!("{gt}:{p}")),
    };
    format!("chr22\t{pos}\t.\t{r}\t{alt}\t50\tPASS\t.\t{fmt}\t{smp}\n")
}

/// What `bcftools +fill-from-fasta -- -c INFO/PGP_JOIN_SEQ` does: the
/// reference bases from POS for the length of REF.
fn fill(vcf: &str) -> String {
    vcf.lines()
        .map(|line| {
            if line.starts_with('#') {
                return line.to_owned() + "\n";
            }
            let mut f: Vec<String> = line.split('\t').map(str::to_owned).collect();
            let seq = ref_at(f[1].parse().unwrap(), f[3].len());
            f[7] = f[7].replacen("PGP_JOIN_SEQ=.", &format!("PGP_JOIN_SEQ={seq}"), 1);
            f.join("\t") + "\n"
        })
        .collect()
}

fn body(text: &str) -> Vec<String> {
    text.lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// prepare | fill | finish, and build, which must write the same records:
/// the records as fields, and the whole output of finish.
fn joined(fasta: &str, records: &[String]) -> (Vec<Vec<String>>, String) {
    let input = HEAD.to_owned() + &records.concat();
    let mut prep = Vec::new();
    join::join(
        Mode::Prepare,
        WINDOW,
        Cursor::new(input.as_bytes()),
        &mut prep,
    )
    .expect("prepare");
    let mut done = Vec::new();
    join::finish(
        Cursor::new(fill(&String::from_utf8(prep).unwrap()).into_bytes()),
        &mut done,
    )
    .expect("finish");
    let done = String::from_utf8(done).unwrap();
    let mut ix = Indexed::open(fasta).expect("the test reference");
    let mut built = Vec::new();
    join::join(
        Mode::Build(&mut ix),
        WINDOW,
        Cursor::new(input.as_bytes()),
        &mut built,
    )
    .expect("build");
    let built = String::from_utf8(built).unwrap();
    assert_eq!(
        body(&built),
        body(&done),
        "build must write the records of prepare and finish"
    );
    let recs = body(&done)
        .iter()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    (recs, done)
}

/// POS, REF, ALT, FORMAT and the sample of each record.
fn core(recs: &[Vec<String>]) -> Vec<[&str; 5]> {
    recs.iter()
        .map(|r| [&*r[1], &*r[3], &*r[4], &*r[8], &*r[9]])
        .collect()
}

#[test]
fn neighbouring_records_join_into_one_allele() {
    let w = Scratch::new("join");
    let fasta = write_reference(&w.0);
    let j = |records: &[String]| joined(&fasta, records);

    // 1. unphased het pair, 3 bp apart: AGC>A (a 2 bp deletion) and T>TGGG
    let (out, raw) = j(&[
        rec(30000037, "AGC", "A", "0/1", None),
        rec(30000042, "T", "TGGG", "0/1", None),
    ]);
    assert_eq!(
        core(&out),
        [["30000037", "AGCCAT", "ACATGGG", "GT", "0/1"]],
        "two het records 3 bp apart become one REF/ALT with the reference bases between them"
    );
    assert_eq!(
        out[0][7], "PGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG",
        "the combined record names its two parts"
    );
    assert!(
        raw.contains("ID=PGP_JOINED_FROM") && !raw.contains("ID=PGP_JOIN_SEQ"),
        "the header keeps PGP_JOINED_FROM and drops the temporary PGP_JOIN_SEQ: {raw}"
    );
    let (out, _) = j(&[
        rec(30000203, "T", "A", "1/0", None),
        rec(30000208, "T", "C", "0/1", None),
    ]);
    assert_eq!(
        core(&out),
        [["30000203", "TGTAAT", "AGTAAC", "GT", "0/1"]],
        "1/0 and 0/1 (unphased) are the same het call and are joined"
    );

    // 2. phase sets
    let (out, _) = j(&[
        rec(30000203, "T", "A", "0|1", Some(100)),
        rec(30000208, "T", "C", "0|1", Some(100)),
    ]);
    assert_eq!(
        core(&out),
        [["30000203", "TGTAAT", "AGTAAC", "GT:PS", "0|1:100"]],
        "phased, same phase set and GT: joined, with GT and PS"
    );
    let none: [[&str; 5]; 0] = [];
    for (desc, a, b) in [
        (
            "phased on different phase sets: not joined",
            rec(30000203, "T", "A", "0|1", Some(100)),
            rec(30000208, "T", "C", "0|1", Some(200)),
        ),
        (
            "phased on the two haplotypes of one phase set: not joined",
            rec(30000203, "T", "A", "0|1", Some(100)),
            rec(30000208, "T", "C", "1|0", Some(100)),
        ),
        (
            "one phased and one unphased record: not joined",
            rec(30000203, "T", "A", "0|1", Some(100)),
            rec(30000208, "T", "C", "0/1", None),
        ),
        // 3. what is not one allele, and the window's edge
        (
            "a het next to a hom: not joined",
            rec(30000203, "T", "A", "0/1", None),
            rec(30000208, "T", "C", "1/1", None),
        ),
        (
            "unphased half calls (./1): not joined",
            rec(30000203, "T", "A", "./1", None),
            rec(30000208, "T", "C", "./1", None),
        ),
        (
            "a phased half call and an ALT on the other haplotype: not joined",
            rec(30000203, "T", "A", "1|.", Some(100)),
            rec(30000208, "T", "C", "0|1", Some(100)),
        ),
        (
            "overlapping records: not joined",
            rec(30000203, "TGT", "T", "0/1", None),
            rec(30000204, "G", "A", "0/1", None),
        ),
        (
            "11 bp apart: not joined",
            rec(30000203, "T", "A", "0/1", None),
            rec(30000214, "C", "A", "0/1", None),
        ),
        // 5. a part whose REF is not the reference's base (30000203 is T)
        (
            "a part that disagrees with the reference: dropped by finish",
            rec(30000203, "C", "A", "0/1", None),
            rec(30000208, "T", "C", "0/1", None),
        ),
    ] {
        let (out, raw) = j(&[a, b]);
        assert_eq!(core(&out), none, "{desc}: {raw}");
    }
    for (a, b, gt) in [
        ("1|.", "1|0", "1|."),
        ("1|0", "1|.", "1|."),
        (".|1", ".|1", ".|1"),
    ] {
        let (out, _) = j(&[
            rec(30000203, "T", "A", a, Some(100)),
            rec(30000208, "T", "C", b, Some(100)),
        ]);
        let sample = format!("{gt}:100");
        assert_eq!(
            core(&out),
            [["30000203", "TGTAAT", "AGTAAC", "GT:PS", sample.as_str()]],
            "phased half calls {a} and {b} in one phase set: joined, '.' kept on the other haplotype"
        );
    }
    let (out, _) = j(&[
        rec(30000203, "T", "A", "0/1", None),
        rec(30000213, "G", "T", "0/1", None),
    ]);
    let three: Vec<[&str; 3]> = core(&out).iter().map(|c| [c[0], c[1], c[2]]).collect();
    assert_eq!(
        three,
        [["30000203", "TGTAATCCCAG", "AGTAATCCCAT"]],
        "10 bp apart: joined"
    );

    // 4. a single SNV
    let (out, _) = j(&[rec(30000210, "C", "T", "0/1", None)]);
    assert_eq!(core(&out), none, "a single SNV gives no combined record");

    // 6. POS order: a 12 bp deletion that ends after the two records inside it
    let (out, raw) = j(&[
        rec(30000201, "CCTGTAATCCCA", "C", "0/1", None),
        rec(30000205, "T", "G", "0/1", None),
        rec(30000208, "T", "C", "0/1", None),
        rec(30000214, "C", "A", "0/1", None),
    ]);
    let pos: Vec<&str> = out.iter().map(|r| r[1].as_str()).collect();
    assert_eq!(
        pos,
        ["30000201", "30000205", "30000205", "30000208"],
        "four combined records, in POS order: {raw}"
    );

    // mode is required
    let args = ["join-adjacent".to_owned()];
    assert_eq!(
        og_clinvar::cli::main(&args),
        2,
        "no mode: exit 2 with a message"
    );
}

const CV: &str = "GENEINFO={g};CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter";

fn cv(g: &str) -> String {
    CV.replace("{g}", g)
}

/// The tag pass: its output and the INFO of each record by POS.
fn tag(hits: &str) -> (String, std::collections::HashMap<String, String>) {
    let mut out = Vec::new();
    join::tag(Cursor::new(hits.as_bytes()), WINDOW, &mut out).expect("tag");
    let out = String::from_utf8(out).unwrap();
    let info = body(&out)
        .iter()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[1].to_owned(), f[7].to_owned())
        })
        .collect();
    (out, info)
}

#[test]
fn the_tag_pass_marks_parts_and_one_allele_pairs() {
    let lines = [
        format!(
            "chr22\t30000037\t9900001\tAGC\tA\t50\tPASS\t{}\tGT\t0/1",
            cv("SYNTH1:1")
        ),
        format!(
            "chr22\t30000038\t9900003\tGCCAT\tCATGGG\t.\tPASS\tPGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG;{}\tGT\t0/1",
            cv("SYNTH1:1").replace("single_submitter", "multiple_submitters,_no_conflicts")
        ),
        format!(
            "chr22\t30000042\t9900002\tT\tTGGG\t50\tPASS\t{}\tGT\t0/1",
            cv("SYNTH1:1")
        ),
        format!(
            "chr22\t30000210\t9900004\tC\tT\t50\tPASS\t{}\tGT\t0/1",
            cv("SYNTH2:2")
        ),
        format!(
            "chr22\t30000300\t9900005\tA\tG\t50\tPASS\t{}\tGT\t0/1",
            cv("SYNTH3:3")
        ),
        format!(
            "chr22\t30000305\t9900006\tC\tT\t50\tPASS\t{}\tGT\t0/1",
            cv("SYNTH3:3")
        ),
        format!(
            "chr22\t30000400\t9900007\tA\tG\t50\tPASS\t{}\tGT\t0/1",
            cv("SYNTH4:4")
        ),
        format!(
            "chr22\t30000404\t9900008\tC\tT\t50\tPASS\t{}\tGT\t0/1",
            cv("SYNTH5:5")
        ),
        format!(
            "chr22\t30000500\t9900009\tA\tG\t50\tPASS\t{}\tGT:PS\t0|1:7",
            cv("SYNTH6:6")
        ),
        format!(
            "chr22\t30000504\t9900010\tC\tT\t50\tPASS\t{}\tGT:PS\t1|0:7",
            cv("SYNTH6:6")
        ),
        format!(
            "chr22\t30000600\t9900011\tA\tG\t50\tPASS\t{}\tGT:PS\t1|.:8",
            cv("SYNTH7:7")
        ),
        format!(
            "chr22\t30000605\t9900012\tC\tT\t50\tPASS\t{}\tGT:PS\t1|0:8",
            cv("SYNTH7:7")
        ),
    ];
    let (out, info) = tag(&(HEAD.to_owned() + &lines.join("\n") + "\n"));
    let got = body(&out);
    assert!(
        got.len() == lines.len() && out.starts_with("##fileformat"),
        "tag: header and every record kept: {out}"
    );
    let i = |p: &str| info.get(p).map_or("", String::as_str);
    let part = "PGP_PART_OF=chr22:30000038:GCCAT:CATGGG";
    assert!(
        i("30000037").ends_with(part) && i("30000042").ends_with(part),
        "tag: both parts are marked part of the combined allele: {info:?}"
    );
    assert!(
        !i("30000038").contains("PGP_PART_OF") && !i("30000038").contains("PGP_ONE_ALLELE"),
        "tag: the combined hit itself gets no part mark: {info:?}"
    );
    assert!(
        got.contains(&lines[3]),
        "tag: the single SNV hit is unchanged"
    );
    assert!(
        i("30000300").ends_with("PGP_ONE_ALLELE=chr22:30000305:C:T")
            && i("30000305").ends_with("PGP_ONE_ALLELE=chr22:30000300:A:G"),
        "tag: two hits in one gene, 5 bp apart, no combined record: flagged both ways: {info:?}"
    );
    assert!(
        lines[6..10].iter().all(|l| got.contains(l)),
        "tag: hits in two genes, or on the two haplotypes of a phase set, are unchanged"
    );
    assert!(
        i("30000600").ends_with("PGP_ONE_ALLELE=chr22:30000605:C:T")
            && i("30000605").ends_with("PGP_ONE_ALLELE=chr22:30000600:A:G"),
        "tag: a phased half call (1|.) and 1|0 in one gene and phase set: flagged both ways: {info:?}"
    );
}

fn fixture(name: &str) -> String {
    let p = repo().join("tests/fixtures/vcf").join(name);
    std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

/// VCF left-trimming as `bcftools norm` writes a combined record: drop a
/// first base REF and ALT share while both keep one.
fn left_trim(pos: u64, r: &str, a: &str) -> (u64, String, String) {
    let (mut pos, mut r, mut a) = (pos, r.to_owned(), a.to_owned());
    while r.len() > 1 && a.len() > 1 && r.as_bytes()[0] == a.as_bytes()[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }
    (pos, r, a)
}

#[test]
fn the_fixture_delins_matches_clinvar_once_as_one_allele() {
    let w = Scratch::new("fixture");
    let fasta = write_reference(&w.0);
    let sample = fixture("combined_allele_sample.vcf");
    let clinvar = fixture("combined_allele_clinvar.vcf");
    let sample_recs: Vec<Vec<String>> = body(&sample)
        .iter()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    for r in &sample_recs {
        assert_eq!(
            ref_at(r[1].parse().unwrap(), r[3].len()),
            r[3],
            "the fixture's REF bases are the reference's at chr22:{}",
            r[1]
        );
    }

    let mut ix = Indexed::open(&fasta).expect("the test reference");
    let mut built = Vec::new();
    join::join(
        Mode::Build(&mut ix),
        WINDOW,
        Cursor::new(sample.as_bytes()),
        &mut built,
    )
    .expect("build");
    let combined: Vec<Vec<String>> = body(&String::from_utf8(built).unwrap())
        .iter()
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect();
    assert_eq!(
        combined.len(),
        1,
        "the delins written as two records gives one combined record: {combined:?}"
    );
    let c = &combined[0];
    let (pos, r, a) = left_trim(c[1].parse().unwrap(), &c[3], &c[4]);

    // ClinVar's records by allele.
    let cv: std::collections::HashMap<(u64, String, String), Vec<String>> = body(&clinvar)
        .iter()
        .map(|l| {
            let f: Vec<String> = l.split('\t').map(str::to_owned).collect();
            ((f[1].parse().unwrap(), f[3].clone(), f[4].clone()), f)
        })
        .collect();
    let two_star = cv.get(&(pos, r.clone(), a.clone())).unwrap_or_else(|| {
        panic!("the combined allele chr22:{pos} {r}>{a} is no ClinVar record of the fixture")
    });
    assert!(
        two_star[2] == "9900003"
            && two_star[7]
                .contains("CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts"),
        "the combined allele is ClinVar's 2-star record 9900003: {two_star:?}"
    );
    assert_eq!(
        c[7], "PGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG",
        "the combined record names its two parts"
    );
    for (p, r, a, id) in [
        (30000037, "AGC", "A", "9900001"),
        (30000042, "T", "TGGG", "9900002"),
    ] {
        assert_eq!(
            cv.get(&(p, r.to_owned(), a.to_owned()))
                .map(|f| f[2].as_str()),
            Some(id),
            "the part {p} {r}>{a} is ClinVar's 1-star record {id}"
        );
    }

    // The hits as step 06 writes them before the tag pass: each candidate
    // (the sample's records and the combined one) that is a ClinVar allele,
    // with ClinVar's ID and INFO, in POS order.
    let mut candidates = sample_recs.clone();
    let mut trimmed = c.clone();
    (trimmed[1], trimmed[3], trimmed[4]) = (pos.to_string(), r, a);
    candidates.push(trimmed);
    candidates.sort_by_key(|f| f[1].parse::<u64>().unwrap());
    let hits_lines: Vec<String> = candidates
        .iter()
        .filter_map(|f| {
            let m = cv.get(&(f[1].parse().unwrap(), f[3].clone(), f[4].clone()))?;
            let info = if f[7] == "." {
                m[7].clone()
            } else {
                format!("{};{}", f[7], m[7])
            };
            let mut h = f.clone();
            (h[2], h[7]) = (m[2].clone(), info);
            Some(h.join("\t"))
        })
        .collect();
    let (out, info) = tag(&(HEAD.to_owned() + &hits_lines.join("\n") + "\n"));
    let got = body(&out);
    assert_eq!(
        got.len(),
        4,
        "four hits in all: the combined allele is matched once: {out}"
    );
    let i = |p: &str| info.get(p).map_or("", String::as_str);
    assert!(
        i("30000038").contains("PGP_JOINED_FROM=chr22:30000037:AGC:A|chr22:30000042:T:TGGG"),
        "the 2-star hit is the combined allele: {info:?}"
    );
    assert!(
        ["30000037", "30000042"]
            .iter()
            .all(|p| i(p).contains("PGP_PART_OF=chr22:30000038:GCCAT:CATGGG")),
        "both 1-star parts are hits, marked part of the combined allele: {info:?}"
    );
    let snv = got
        .iter()
        .find(|l| l.split('\t').nth(1) == Some("30000210"))
        .expect("the SNV hit");
    assert!(
        snv.split('\t').nth(2) == Some("9900004") && !snv.contains("PGP_"),
        "the single SNV hit is unchanged (its ClinVar fields, no combined-allele tag): {snv}"
    );
}

#[test]
fn a_split_multiallelic_call_is_not_het_and_every_listed_site_gets_a_call() {
    let tract = ref_at(30000011, 5);
    assert_eq!(
        tract, "CAAAA",
        "poly-A: the reference holds C and a run of A at chr22:30000011"
    );

    // Step 06's hits after `bcftools norm -m -any --multi-overlaps .`: the
    // ClinVar allele of the 1/2 record (the 4-base deletion) is 1/., the
    // other copy holds CAA.
    let hits_vcf = format!(
        "{HEAD}chr22\t30000011\t9900005\t{tract}\tC\t50\tPASS\tGENEINFO=SYNTH4:990004;CLNSIG=Pathogenic;\
CLNREVSTAT=criteria_provided,_single_submitter\tGT:GQ:DP\t1/.:45:30\n\
chr22\t30000210\t9900004\tC\tT\t50\tPASS\tGENEINFO=SYNTH2:990002;CLNSIG=Pathogenic;\
CLNREVSTAT=criteria_provided,_single_submitter\tGT:GQ:DP\t0/1:40:30\n"
    );
    let mut rows = Vec::new();
    hits::hits(Cursor::new(hits_vcf.as_bytes()), &mut rows).expect("og clinvar hits");
    let rows = String::from_utf8(rows).unwrap();
    let zyg: std::collections::HashMap<&str, &str> = rows
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[2], f[5])
        })
        .collect();
    assert_eq!(
        zyg.get("30000011").copied(),
        Some("one copy, other copy not called ref (1/.)"),
        "poly-A: the hit is not worded het: {rows}"
    );
    assert_eq!(
        zyg.get("30000210").copied(),
        Some("het"),
        "poly-A: the plain 0/1 SNV hit beside it is still het"
    );

    // The always-report list: four sites, and the calls as step 06's
    // `bcftools query '%CHROM\t%POS\t%END\t%REF\t%ALT\t[%GT]\t[%GQ]\t[%DP]\n'`
    // writes them from the split PASS VCF and the gVCF.
    let list = format!(
        "# four sites of the test\nchrom\tpos\tref\talt\tgene\thgvs\tclinvar_id\n\
chr1\t100\tA\tG\tSYNTHA\tc.1A>G\t1\n\
chr22\t30000011\t{tract}\tC\tSYNTH4\tc.2del\t9900005\n\
chr22\t30000030\tA\tG\tSYNTHB\tc.3A>G\t2\n\
chr22\t30000210\tC\tT\tSYNTH2\tc.4C>T\t9900004\n"
    );
    let vcf_calls = format!(
        "chr22\t30000011\t30000015\t{tract}\tC\t1/.\t45\t30\n\
chr22\t30000011\t30000015\t{tract}\tCAA\t./1\t45\t30\n\
chr22\t30000210\t30000210\tC\tT\t0/1\t40\t30\n"
    );
    let gvcf_calls = format!(
        "chr22\t30000011\t30000015\t{tract}\tC,CAA,<*>\t1/2\t45\t30\n\
chr22\t30000016\t30000209\tA\t<*>\t0/0\t50\t30\n\
chr22\t30000210\t30000210\tC\tT,<*>\t0/1\t40\t30\n"
    );
    let o = sites::Options {
        min_gq: "30".to_owned(),
        min_gq_value: 30,
        gvcf: true,
    };
    let mut out = Vec::new();
    sites::deep_intronic(
        list.as_bytes(),
        vcf_calls.as_bytes(),
        gvcf_calls.as_bytes(),
        &o,
        &mut out,
    )
    .expect("og clinvar deep-intronic");
    let out = String::from_utf8(out).unwrap();
    let got: std::collections::BTreeMap<&str, (&str, &str)> = out
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[1], (f[7], f[8]))
        })
        .collect();
    assert_eq!(
        got.keys().copied().collect::<Vec<_>>(),
        ["100", "30000011", "30000030", "30000210"],
        "always-report: one row per site, all four: {out}"
    );
    assert_eq!(
        got["30000210"],
        ("het", "GQ 40, DP 30"),
        "always-report: the carried SNV gives het with its GQ and DP"
    );
    assert_eq!(
        got["30000030"],
        ("ref", "GQ 50 (gVCF block)"),
        "always-report: a site inside a gVCF reference block gives ref with the block's GQ"
    );
    assert_eq!(
        got["30000011"],
        ("one copy, other copy not called ref (1/.)", "GQ 45, DP 30"),
        "always-report: the poly-A site gives the 1/. call, not het"
    );
    assert_eq!(
        got["100"],
        (
            "no call",
            "no record in the VCF or the gVCF covers the site"
        ),
        "always-report: a site no record covers gives no call, and why"
    );
}
