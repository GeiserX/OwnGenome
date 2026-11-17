//! `og dev demo-sample`, the made-up sample of the docs pictures and the
//! live-site demo (it replaced `tests/demo/make_demo_sample.py` as the
//! writer the retake script and the live demo run):
//!
//!   1. DEMO-001 is written, step 24 reads every section of it as ok, and
//!      `og dev demo-sample check` says so; with one section's file taken
//!      away the check fails and names that section;
//!   2. the same seed writes the same bytes, and another seed other values;
//!   3. DEMO-003 carries the planted results to act on (the KCNH2 ClinVar
//!      hit, the class 5 deletion over it, the gene in the symptom ranking)
//!      and DEMO-001 none of them; DEMO-002 is female, with no chrY and no Y
//!      files, and carries the FMR1 premutation, the CHIP call and the
//!      MoChA loss;
//!   4. the command line: a name not of the form DEMO-001, an unknown
//!      option and a folder that exists are refused.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

/// A fresh folder for one case.
fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-port-demo-sample-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// `og ARGS`, run from the repository (step 24 reads `bin/` tables from there).
fn og(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_og"))
        .args(args)
        .current_dir(REPO)
        .output()
        .expect("og runs")
}

/// Writes one person into `root`; panics with og's stderr when it fails.
fn demo(root: &Path, extra: &[&str]) {
    let mut args = vec![
        "dev",
        "demo-sample",
        "--genome-dir",
        root.to_str().unwrap(),
        "--records",
        "3000",
        "--repo",
        REPO,
    ];
    args.extend(extra);
    let out = og(&args);
    assert!(
        out.status.success(),
        "og dev demo-sample {extra:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// Step 24's summary.json of `sample` in `root`, then the check form:
/// (exit code, stderr).
fn collect_and_check(root: &Path, sample: &str) -> (i32, String) {
    let d = root.join(sample);
    let out = og(&[
        "summary",
        "collect",
        "--sample",
        sample,
        "--sample-dir",
        d.to_str().unwrap(),
        "--assets",
        &format!("{REPO}/assets"),
        "--out",
        d.join("summary.json").to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let out = og(&[
        "dev",
        "demo-sample",
        "check",
        "--genome-dir",
        root.to_str().unwrap(),
        "--sample",
        sample,
    ]);
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// The text of a file, gunzipped when it ends in `.gz`.
fn text(path: &Path) -> String {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    if path.extension().is_some_and(|e| e == "gz") {
        let out = Command::new("gzip").arg("-dc").arg(path).output().unwrap();
        String::from_utf8(out.stdout).unwrap()
    } else {
        String::from_utf8(bytes).unwrap()
    }
}

/// Every file under `root`, relative, sorted.
fn files(root: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut todo = vec![root.to_path_buf()];
    while let Some(d) = todo.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                todo.push(p);
            } else {
                out.push(p.strip_prefix(root).unwrap().to_path_buf());
            }
        }
    }
    out.sort();
    out
}

#[test]
fn demo_001_feeds_every_section_and_the_check_catches_a_missing_one() {
    let root = dir("check");
    demo(&root, &[]);
    let (code, err) = collect_and_check(&root, "DEMO-001");
    assert_eq!(code, 0, "{err}");
    // The VCF holds about the records asked for, sorted, on every chromosome.
    let vcf = text(&root.join("DEMO-001/vcf/DEMO-001.vcf.gz"));
    let recs: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(
        (2980..=3020).contains(&recs.len()),
        "{} records",
        recs.len()
    );
    assert!(recs.iter().any(|l| l.starts_with("chrY\t")));
    // Red control: without step 52's table the check fails and names it.
    std::fs::remove_file(root.join("DEMO-001/vntyper/DEMO-001_vntyper.tsv")).unwrap();
    let (code, err) = collect_and_check(&root, "DEMO-001");
    assert_eq!(code, 1, "{err}");
    assert!(err.contains("FAIL: vntyper is missing"), "{err}");
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn the_same_seed_writes_the_same_bytes() {
    let (a, b, c) = (dir("seed-a"), dir("seed-b"), dir("seed-c"));
    demo(&a, &[]);
    demo(&b, &[]);
    demo(&c, &["--seed", "2"]);
    let (fa, fb) = (files(&a), files(&b));
    assert_eq!(fa, fb);
    assert!(fa.len() > 80, "{} files", fa.len());
    for f in &fa {
        assert_eq!(
            std::fs::read(a.join(f)).unwrap(),
            std::fs::read(b.join(f)).unwrap(),
            "{} differs between two runs of seed 1",
            f.display()
        );
    }
    let vcf = "DEMO-001/vcf/DEMO-001.vcf.gz";
    assert_ne!(
        text(&a.join(vcf)),
        text(&c.join(vcf)),
        "seed 2 wrote seed 1's VCF"
    );
    for d in [a, b, c] {
        let _ = std::fs::remove_dir_all(&d);
    }
}

#[test]
fn the_planted_findings_land_on_their_people() {
    let root = dir("planted");
    demo(&root, &[]);
    demo(
        &root,
        &["--sample", "DEMO-002", "--sex", "female", "--seed", "2"],
    );
    demo(&root, &["--sample", "DEMO-003", "--seed", "3"]);
    let read = |rel: &str| text(&root.join(rel));

    let hits = read("DEMO-003/clinvar/DEMO-003_clinvar_hits.vcf");
    assert!(
        hits.lines().any(|l| l.starts_with("chr7\t150960000\t")
            && l.contains("\tC\tT\t")
            && l.contains("GENEINFO=KCNH2:3757;CLNSIG=Pathogenic;")),
        "{hits}"
    );
    let annotsv = read("DEMO-003/annotsv/DEMO-003_sv_class4_5_checked.tsv");
    assert!(
        annotsv.contains("chr7_150955000_150965000_DEL_1\tMantaDEL:2\tchr7\t150955000\t150965000\t-10000\tDEL\tKCNH2\t5\t"),
        "{annotsv}"
    );
    assert!(annotsv.contains("\tG0533\t4\t"), "{annotsv}");
    let ranked = read("DEMO-003/exomiser/DEMO-003_exomiser_top.tsv");
    assert!(
        ranked
            .lines()
            .any(|l| l.starts_with("3\tKCNH2\t") && l.contains("7:150960000 C>T")),
        "{ranked}"
    );
    for rel in [
        "DEMO-001/clinvar/DEMO-001_clinvar_hits.vcf",
        "DEMO-001/annotsv/DEMO-001_sv_class4_5_checked.tsv",
        "DEMO-001/exomiser/DEMO-001_exomiser_top.tsv",
    ] {
        assert!(!read(rel).contains("KCNH2"), "{rel} names the planted gene");
    }

    assert!(!read("DEMO-002/vcf/DEMO-002.vcf.gz").contains("\nchrY\t"));
    assert!(!root.join("DEMO-002/y_haplogroup").exists());
    assert!(!root.join("DEMO-002/y_str").exists());
    let stranger = read("DEMO-002/expansion_hunter/DEMO-002_eh_stranger.vcf");
    let fmr1 = stranger
        .lines()
        .find(|l| l.contains("VARID=FMR1;"))
        .unwrap();
    assert!(
        fmr1.contains("STR_STATUS=pre_mutation") && fmr1.contains("/62:"),
        "{fmr1}"
    );
    assert!(read("DEMO-002/chip/DEMO-002_chip.tsv").contains("DNMT3A\tchr2\t25234373\t"));
    assert!(read("DEMO-002/mosaic/DEMO-002_mocha.tsv").contains("chr20\t30100000\t42300000\t"));
    assert!(read("DEMO-002/indexcov/indexcov-indexcov.ped").contains("\t2\t-9\t1.99\t0.01\t"));
    // Only DEMO-001 carries the couple table: the family page reads one copy.
    assert!(
        root.join("DEMO-001/couple/DEMO-001_couple_carrier.tsv")
            .is_file()
    );
    assert!(!root.join("DEMO-002/couple").exists());
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn the_command_line_refuses_bad_input() {
    let root = dir("cli");
    let r = root.to_str().unwrap();
    let out = og(&["dev", "demo-sample", "--genome-dir", r, "--sample", "ALICE"]);
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("an invented name of the form DEMO-001"));
    let out = og(&[
        "dev",
        "demo-sample",
        "--genome-dir",
        r,
        "--og-scaffold-probe",
    ]);
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("unrecognized argument"));
    std::fs::create_dir_all(root.join("DEMO-001")).unwrap();
    let out = og(&["dev", "demo-sample", "--genome-dir", r]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("exists: give an empty --genome-dir"));
    let _ = std::fs::remove_dir_all(&root);
}
