//! og setup check reads the component table og setup install installs from.
//! These tests hold the two together: every installed path of the table is
//! read by a check or listed in NOT_CHECKED with its reason, and the lists
//! the check keeps beside the table (common.sh's DATA_FILES, the academic
//! scores of assets/pgs_scores.tsv) say what their sources say. The output
//! itself is compared with validate-setup.sh's by tests/setup/check-ab.sh.

use og_setup::check::{
    self, Check, Colors, DATA_FILES, NOT_CHECKED, Report, days_since, fixed, pgs_rows, uncovered,
    validate_sample, vcf_header_samples,
};
use og_setup::components::TABLE;
use og_setup::install::Vars;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Command;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .canonicalize()
        .unwrap()
}

/// A scratch folder of this test process, removed first.
fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-setup-check-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// The (component, template) pairs one whole check asks the table for, run
/// over GENOME with no container engine (so no image is inspected).
fn seen_over(genome: &Path) -> BTreeSet<(String, String)> {
    let mut env = BTreeMap::new();
    for k in ["PATH", "HOME"] {
        if let Ok(v) = std::env::var(k) {
            env.insert(k.to_owned(), v);
        }
    }
    env.insert("GENOME_DIR".into(), genome.to_string_lossy().into_owned());
    env.insert("CONTAINER_ENGINE".into(), "og-check-test-no-engine".into());
    let c = Check::new(&root(), Vars::from_map(env), "").unwrap();
    let mut out = Vec::new();
    let mut r = Report::new(&mut out, Colors::default());
    c.run(&mut r);
    c.seen.into_inner()
}

fn rows() -> Vec<(&'static str, &'static [&'static str])> {
    TABLE.iter().map(|c| (c.key, c.installed)).collect()
}

#[test]
fn every_installed_path_is_checked_or_listed() {
    let g = scratch("cover");
    let seen = seen_over(&g);
    let bad = uncovered(rows(), &seen);
    assert!(bad.is_empty(), "{}", bad.join("\n"));
    // A filled tree asks for the same paths as an empty one: every path is
    // asked for before the branch that reads it.
    std::fs::create_dir_all(g.join("reference")).unwrap();
    std::fs::write(g.join("reference/Homo_sapiens_assembly38.fasta"), "").unwrap();
    assert_eq!(seen_over(&g), seen);
    let _ = std::fs::remove_dir_all(&g);
}

#[test]
fn a_new_row_whose_path_no_check_reads_fails() {
    let g = scratch("red");
    let seen = seen_over(&g);
    let mut planted = rows();
    planted.push(("planted-data", &["reference/planted/DONE"]));
    let bad = uncovered(planted, &seen);
    assert_eq!(
        bad,
        ["planted-data: reference/planted/DONE is read by no check and not in NOT_CHECKED"]
    );
    // And a path moved inside a row: the old one is no longer in any row.
    let moved: Vec<_> = rows()
        .into_iter()
        .map(|(k, t)| {
            if k == "kir-data" {
                (k, &["kir/IPD-KIR_${KIR_DB_RELEASE}/moved.dat"][..])
            } else {
                (k, t)
            }
        })
        .collect();
    let bad = uncovered(moved, &seen);
    assert!(
        bad.iter().any(|l| l.contains("moved.dat is read by no check"))
            && bad
                .iter()
                .any(|l| l.contains("kir.dat is read by a check but in no row")),
        "{bad:?}"
    );
    let _ = std::fs::remove_dir_all(&g);
}

#[test]
fn not_checked_gives_a_reason_for_each_path() {
    for (k, t, why) in NOT_CHECKED {
        assert!(!why.trim().is_empty(), "{k}: {t} has no reason");
    }
}

#[test]
#[should_panic(expected = "the kir-data row has no installed path")]
fn a_check_of_a_path_the_row_does_not_have_panics() {
    let g = scratch("panic");
    let env = BTreeMap::from([(
        "GENOME_DIR".to_owned(),
        g.to_string_lossy().into_owned(),
    )]);
    let c = Check::new(&root(), Vars::from_map(env), "").unwrap();
    c.table("kir-data", "kir/elsewhere.dat");
}

#[test]
fn academic_scores_are_the_prs_academic_row() {
    let tsv: BTreeSet<String> = pgs_rows(&root().join("assets/pgs_scores.tsv"))
        .into_iter()
        .filter(|(_, lic)| lic == "academic-research")
        .map(|(id, _)| format!("prs_scores/{id}.txt.gz"))
        .collect();
    let row: BTreeSet<String> = og_setup::components::by_key("prs-academic")
        .unwrap()
        .installed
        .iter()
        .map(|t| (*t).to_owned())
        .collect();
    assert!(!tsv.is_empty());
    assert_eq!(tsv, row);
}

#[test]
fn data_files_are_common_sh_s() {
    let common = std::fs::read_to_string(root().join("scripts/lib/common.sh")).unwrap();
    let line = common
        .lines()
        .find_map(|l| l.strip_prefix("DATA_FILES=\""))
        .and_then(|l| l.strip_suffix('"'))
        .expect("DATA_FILES in common.sh");
    assert_eq!(line.split(' ').collect::<Vec<_>>(), DATA_FILES);
}

/// Runs bash SCRIPT with common.sh sourced; its stdout and status.
fn bash(script: &str, args: &[&str]) -> (String, bool) {
    let common = root().join("scripts/lib/common.sh");
    let o = Command::new("bash")
        .arg("-c")
        .arg(format!(". \"$0\" >/dev/null 2>&1; {script}"))
        .arg(&common)
        .args(args)
        .env("GENOME_DIR", std::env::temp_dir())
        .output()
        .unwrap();
    (
        String::from_utf8_lossy(&o.stdout).into_owned(),
        o.status.success(),
    )
}

#[test]
fn sample_names_are_refused_as_common_sh_refuses_them() {
    for name in [
        "HG002", "sample_1", "a.b-c", "-x", "../up", "a/b", "a b", ".", "..", "x;y", "Ñ", "",
        "a\nb", ".hidden",
    ] {
        let (_, ok) = bash("( validate_sample \"$1\" ) 2>/dev/null", &[name]);
        assert_eq!(
            validate_sample(name).is_ok(),
            ok,
            "validate_sample {name:?}"
        );
    }
}

#[test]
fn vcf_samples_are_common_sh_s() {
    for header in [
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone two\tx",
        "",
    ] {
        let (sh, _) = bash("printf '%s' \"$1\" | vcf_header_samples", &[header]);
        assert_eq!(
            vcf_header_samples(header),
            sh.trim_end_matches('\n'),
            "{header:?}"
        );
    }
}

#[test]
fn rounding_is_awk_s() {
    for x in [
        0.0, 0.05, 0.15, 0.25, 0.35, 0.45, 1.25, 2.5, 3.5, 15.95, 16.0, 7.449_999, 0.5, 1.5,
        1023.75, 3.149_999_999,
    ] {
        for n in [0usize, 1] {
            let o = Command::new("awk")
                .arg(format!("BEGIN {{ printf \"%.{n}f\", {x:?} }}"))
                .output()
                .unwrap();
            assert_eq!(fixed(x, n), String::from_utf8_lossy(&o.stdout), "{x} .{n}");
        }
    }
}

#[test]
fn clinvar_age_is_date_s() {
    assert_eq!(days_since("2024-02-30"), None);
    assert_eq!(days_since("2023-13-01"), None);
    let o = Command::new("date").args(["-u", "+%Y-%m-%d"]).output().unwrap();
    let today = String::from_utf8_lossy(&o.stdout).trim().to_owned();
    assert_eq!(days_since(&today), Some(0));
    assert!(days_since("2024-02-29").unwrap() > 0);
}

#[test]
fn help_is_on_stdout_and_exits_0() {
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(check::cli(&["--help".into()], &mut out, &mut err), 0);
    assert!(String::from_utf8_lossy(&out).starts_with("usage: og setup check [SAMPLE]"));
    assert!(err.is_empty());
}

/// A Check over GENOME_DIR=/genome with REF_FASTA, as the scripts set them.
fn check_with_ref(reference: &str) -> Check {
    let env = BTreeMap::from([
        ("GENOME_DIR".to_owned(), "/genome".to_owned()),
        ("REF_FASTA".to_owned(), reference.to_owned()),
    ]);
    Check::new(&root(), Vars::from_map(env), "").unwrap()
}

/// check_bam_masked_reference read out of validate-setup.sh, run by bash
/// with REF_FASTA on HEADER and validate-setup.sh's own line formats.
fn bash_masked(reference: &str, header: &str) -> String {
    let fun = Command::new("awk")
        .arg(r"/^check_bam_masked_reference\(\) \{$/ {on = 1} on {print} on && /^}$/ {exit}")
        .arg(root().join("scripts/validate-setup.sh"))
        .output()
        .unwrap();
    let fun = String::from_utf8(fun.stdout).unwrap();
    assert!(fun.contains("_maskedGRC"), "check_bam_masked_reference not found");
    let o = Command::new("bash")
        .arg("-c")
        .arg(format!(
            "set -euo pipefail
pass() {{ echo \"  [OK]    $1\"; }}
warn() {{ echo \"  [WARN]  $1\"; }}
fail() {{ echo \"  [FAIL]  $1\"; }}
{fun}
REF_FASTA=$0 check_bam_masked_reference \"$1\""
        ))
        .arg(reference)
        .arg(header)
        .output()
        .unwrap();
    String::from_utf8_lossy(&o.stdout).into_owned()
}

/// The hand-made headers of tests/test_validate_setup_masked.sh: the @PG
/// lines of minimap2 and bwa-mem2 on either reference, none, one without a
/// command line, merged ids, and _maskedGRC only in read names.
#[test]
fn masked_reference_lines_are_validate_setup_s() {
    let stock = "/genome/reference/GRCh38_no_alt_analysis_set";
    let masked = format!("{stock}_maskedGRC");
    let sq = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr21\tLN:46709983";
    let mm2 = |r: &str| format!("{sq}\n@PG\tID:minimap2\tPN:minimap2\tVN:2.28\tCL:minimap2 -t 16 -a -x sr {r}.sr.mmi R1.fq.gz R2.fq.gz\n@PG\tID:samtools\tPN:samtools\tPP:minimap2\tCL:samtools sort -o x.bam");
    let bwa = |r: &str| format!("{sq}\n@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:2.2.1\tCL:bwa-mem2 mem -t 16 {r}.fasta R1.fq.gz R2.fq.gz");
    let none = format!("{sq}\n@PG\tID:DRAGEN\tPN:DRAGEN\tCL:dragen -r /ref");
    let nocl = format!("{sq}\n@PG\tID:minimap2\tPN:minimap2\tVN:2.28");
    let rgname = |r: &str| format!("{sq}\n@PG\tID:minimap2\tPN:minimap2\tCL:minimap2 -a -x sr -R @RG\\tID:s_maskedGRC\\tSM:s_maskedGRC {r}.sr.mmi s_maskedGRC_R1.fastq.gz s_maskedGRC_R2.fq.gz");
    let merged = |r: &str, id: &str| {
        mm2(r)
            .replacen("ID:minimap2", &format!("ID:{id}"), 1)
            .replacen("PN:minimap2", "PN:mm", 1)
    };
    let both = format!(
        "{}\n{}",
        mm2(stock),
        bwa(&masked).lines().filter(|l| l.starts_with("@PG")).collect::<Vec<_>>().join("\n")
    );
    let quoted = format!("{sq}\n@PG\tID:bwa\tPN:bwa\tCL:bwa mem  \"{masked}.fasta\" 'R1.fq.gz'");
    let headers = [
        mm2(&masked),
        mm2(stock),
        bwa(&masked),
        bwa(stock),
        none.clone(),
        nocl.clone(),
        rgname(stock),
        rgname(&masked),
        merged(stock, "minimap2.1"),
        merged(&masked, "minimap2.1.1"),
        both,
        quoted,
        sq.to_owned(),
    ];
    let mut n = 0;
    for reference in [format!("{masked}.fasta"), format!("{stock}.fasta")] {
        let c = check_with_ref(&reference);
        for h in &headers {
            let mut out = Vec::new();
            let mut r = Report::new(&mut out, Colors::default());
            c.check_bam_masked_reference(&mut r, h, check::M5::Unknown);
            let sh = bash_masked(&reference, h);
            assert_eq!(String::from_utf8(out).unwrap(), sh, "REF_FASTA={reference}\n{h}");
            n += usize::from(!sh.is_empty());
        }
    }
    // Most cases print a line, so an empty compare cannot pass them all.
    assert!(n >= 10, "only {n} cases printed a line");
}

/// read_count_check of common.sh, run by bash, on the counts of
/// tests/test_read_count_gate.sh (with no table, as validate-setup.sh calls it).
#[test]
fn read_count_lines_are_common_sh_s() {
    let cases: &[(&str, &[&str], Option<&str>)] = &[
        ("10", &["5", "5"], None),
        ("6", &["5", "5"], None),
        ("6", &["5", "5"], Some("warn")),
        ("11", &["5", "5"], None),
        ("", &["5", "5"], None),
        ("10", &["5", "x"], None),
        ("10", &[], None),
        ("10", &["10"], Some("maybe")),
        ("0", &["0", "0"], None),
        ("3000000000", &["1500000000", "1500000000"], None),
        ("2999999999", &["1500000000", "1500000000"], None),
        ("007", &["7"], None),
        ("1", &["3"], Some("fail")),
    ];
    for (primary, counts, mode) in cases {
        let common = root().join("scripts/lib/common.sh");
        let mut c2 = Command::new("bash");
        c2.arg("-c")
            .arg(". \"$0\" >/dev/null 2>&1; p=$1; shift; read_count_check \"$p\" \"R1 + R2\" \"\" \"$@\" 2>&1; echo \"rc=$?\"")
            .arg(&common)
            .arg(primary)
            .args(counts.iter())
            .env("GENOME_DIR", std::env::temp_dir());
        match mode {
            Some(m) => c2.env("READ_COUNT_CHECK", m),
            None => c2.env_remove("READ_COUNT_CHECK"),
        };
        let o = c2.output().unwrap();
        let sh = String::from_utf8_lossy(&o.stdout).into_owned();
        let (rc, msg) = check::read_count_check(primary, "R1 + R2", counts, *mode);
        let ours = if msg.is_empty() {
            format!("rc={rc}\n")
        } else {
            format!("{msg}\nrc={rc}\n")
        };
        assert_eq!(ours, sh, "read_count_check {primary} {counts:?} {mode:?}");
    }
}
