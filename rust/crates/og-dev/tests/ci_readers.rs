//! Black-box controls for the CI readers. Every clean fixture is paired with
//! a changed or unreadable input that must fail.
use og_dev::{nf_trace, pharmcat, smoke};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "og-ci-readers-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
    fn file(&self, name: &str, text: &str) -> String {
        let p = self.0.join(name);
        std::fs::write(&p, text).unwrap();
        p.to_str().unwrap().into()
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().into()
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
type Cli = fn(&[String], &mut dyn std::io::Write, &mut dyn std::io::Write) -> u8;
fn run(cli: Cli, args: &[&str]) -> (u8, String, String) {
    let mut out = Vec::new();
    let mut err = Vec::new();
    let code = cli(
        &args.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
        &mut out,
        &mut err,
    );
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}
#[test]
fn unreadable_pharmcat_is_unassessed() {
    let t = Temp::new();
    for p in [
        t.path("missing"),
        t.file("broken", "{\"genes\":"),
        t.file("wrong", "[]"),
        t.file("bad-shape", "{\"genes\": []}"),
    ] {
        for flag in ["--called", "--called-lines"] {
            let (rc, out, err) = run(pharmcat::cli, &["diplotypes", &p, flag]);
            assert_eq!(rc, 2, "{err}");
            assert!(out.is_empty());
            assert!(err.starts_with("unreadable report:"));
        }
    }
    let p = t.file("empty", "{\"genes\": {}}");
    assert_eq!(
        run(pharmcat::cli, &["diplotypes", &p, "--called"]),
        (0, "0\n".into(), "\n".into())
    );
}
#[test]
fn pharmcat_shapes_fallback_and_deduplication() {
    let t = Temp::new();
    let p=t.file("r",r#"{"genes":{"CPIC":{"B":{"sourceDiplotypes":[],"recommendationDiplotypes":[{"allele1":{"name":"*1"}}]},"A":{"sourceDiplotypes":[{"allele1":{"name":"Unknown"},"allele2":{"name":"?"}}]}},"DPWG":{"B":{"sourceDiplotypes":[{"allele1":{"name":"*1"}}]}},"C":{"sourceDiplotypes":[{"allele1":{"name":"NONE"}}]}}}"#);
    assert_eq!(
        run(pharmcat::cli, &["diplotypes", &p]),
        (0, "A\tUnknown/?\nB\t*1/\nC\tNONE/\n".into(), "".into())
    );
    assert_eq!(
        run(pharmcat::cli, &["diplotypes", &p, "--called"]),
        (0, "1\n".into(), "B\n".into())
    );
    let flat = t.file(
        "flat",
        r#"{"genes":{"B":{"sourceDiplotypes":[{"allele1":{"name":"*1"}}]}}}"#,
    );
    assert_eq!(
        run(pharmcat::cli, &["diplotypes", &flat, "--called-lines"]),
        (0, "B\n".into(), "".into())
    );
}
const TRACE: &str = "\"task_id\"\tname\tattempt\texit\tmemory\tcpus\ttime\n1\tX:ALIGN_MINIMAP2 (s)\t1\t137\t34359738368\t8\t86400000\n2\tX:ALIGN_MINIMAP2 (s)\t2\t0\t68719476736\t8\t86400000\n";
#[test]
fn trace_attempts_and_duplicate_controls() {
    let t = Temp::new();
    let p = t.file("trace", TRACE);
    let args = |p: &str| {
        run(
            nf_trace::cli,
            &[
                "attempts",
                p,
                "--name-glob",
                "*ALIGN_MINIMAP2*",
                "--want",
                "1:137:32,2:0:64",
                "--nextflow-exit",
                "0",
            ],
        )
    };
    assert_eq!(args(&p).0, 0);
    let duplicate = t.file(
        "dup",
        &format!("{TRACE}1\tX:ALIGN_MINIMAP2 (s)\t1\t-\t34359738368\t8\t86400000\n"),
    );
    assert_eq!(args(&duplicate).0, 0);
    for (name, text) in [
        ("memory", TRACE.replace("68719476736", "34359738368")),
        ("exit", TRACE.replace("\t137\t", "\t1\t")),
        (
            "missing",
            TRACE.lines().take(2).collect::<Vec<_>>().join("\n"),
        ),
        ("empty", String::new()),
        ("dash", TRACE.replace("68719476736", "-")),
        (
            "conflict",
            format!("{TRACE}1\tX:ALIGN_MINIMAP2 (s)\t1\t0\t34359738368\t8\t86400000\n"),
        ),
    ] {
        assert_ne!(args(&t.file(name, &text)).0, 0, "{name}");
    }
    assert_eq!(args("trace_*.txt").0, 2);
    assert_eq!(run(nf_trace::cli, &["attempts", &p, &p]).0, 2);
}
#[test]
fn trace_requests_and_budget_controls() {
    let t = Temp::new();
    let p = t.file("trace", TRACE);
    let want = t.file(
        "want",
        "process\tcpus\tmemory_gb\ttime_h\nALIGN_MINIMAP2\t8\t32\t24\n",
    );
    let args = |p: &str, cpus: &str, mem: &str| {
        run(
            nf_trace::cli,
            &[
                "requests",
                p,
                "--want",
                &want,
                "--executor-cpus",
                cpus,
                "--executor-mem-gb",
                mem,
            ],
        )
    };
    assert_eq!(args(&p, "16", "96").0, 0);
    assert_eq!(args(&p, "15", "96").0, 1);
    assert_eq!(args(&p, "16", "95").0, 1);
    for (name, text) in [
        ("cpus", TRACE.replace("\t8\t", "\t4\t")),
        ("time", TRACE.replace("86400000", "3600000")),
        ("third", TRACE.replace("\t2\t0\t", "\t3\t0\t")),
        (
            "duplicate-name",
            format!("{TRACE}3\tX:ALIGN_MINIMAP2 (s)\t1\t137\t34359738368\t8\t86400000\n"),
        ),
    ] {
        assert_eq!(args(&t.file(name, &text), "512", "2048").0, 1, "{name}");
    }
}
#[test]
fn quoted_csv_cell_and_missing_column() {
    let t = Temp::new();
    let p = t.file(
        "csv",
        "Type,Filter,Value\r\n\"S,N\",PASS,\"a,b\nnext\"\r\nSNP,PASS,0.9\r\n",
    );
    assert_eq!(
        run(
            smoke::cli,
            &["csv-cell", &p, "Type=SNP,Filter=PASS", "Value"]
        ),
        (0, "0.9\n".into(), "".into())
    );
    assert_eq!(
        run(smoke::cli, &["csv-cell", &p, "Filter=PASS", "Value"]).1,
        "a,b\nnext\n"
    );
    assert_eq!(
        run(smoke::cli, &["csv-cell", &p, "Filter=PASS", "missing"]).0,
        2
    );
}
#[test]
fn smoke_tables_reject_removed_properties() {
    let t = Temp::new();
    for (kind, good, bad) in [
        (
            "hap-ibd",
            "A 1 B 1 20 1 10000000 10\n",
            "A 1 C 1 20 1 10000000 10\n",
        ),
        (
            "ibdmix",
            "h\nS 22 2000000 2600000 4\n",
            "h\nS 22 2000000 2300000 4\n",
        ),
        (
            "rfmix",
            "#codes: POPA=0 POPB=1\n#h\n20 0 10000000 0 0 0 0 0\n20 10000000 20000000 0 0 0 1 1\n",
            "#codes: POPA=0 POPB=1\n20 0 20000000 0 0 0 0 0\n",
        ),
        (
            "eigen",
            "#h\na -2 0 PopA\nb 2 0 PopB\np -1 0 Proj\n",
            "#h\na -2 0 PopA\nb 2 0 PopB\np 1 0 Proj\n",
        ),
        (
            "parascopy",
            "chr5\t1\t2\ta\tb\tc\t4\t99\n",
            "chr5\t1\t2\ta\tb\tc\t*\t99\n",
        ),
    ] {
        let p = t.file("table", good);
        let clean = run(smoke::cli, &["check", kind, &p]);
        assert_eq!(clean.0, 0, "{kind}");
        if kind == "eigen" {
            assert_eq!(clean.1, "eigen: 3 rows; OK\n");
        }
        let p = t.file("table", bad);
        assert_eq!(run(smoke::cli, &["check", kind, &p]).0, 1, "{kind}");
    }
}
#[test]
fn smoke_generators_keep_coordinates_and_reverse_complements() {
    let t = Temp::new();
    let seq = format!(
        "{}{}{}",
        "A".repeat(21000),
        "CA".repeat(12),
        "G".repeat(478976)
    );
    let src = t.file("src", &format!(">slice\n{seq}\n"));
    let orig = t.path("orig");
    let mini = t.path("mini");
    assert_eq!(
        run(
            smoke::cli,
            &["mini-ref", &src, &orig, &mini, "450000", "3000"]
        )
        .0,
        0
    );
    let fa = std::fs::read_to_string(&mini).unwrap();
    let s = fa.lines().skip(1).collect::<String>();
    assert_eq!(s.len(), 503000);
    assert_eq!(&s[..450000], &seq[..450000]);
    assert_eq!(&s[453000..], &seq[450000..]);
    let cat = t.path("cat");
    assert_eq!(
        run(smoke::cli, &["eh-catalog", &mini, &cat, "449000"]).1,
        "repeat 21000 21024 12 units\n"
    );
    assert!(
        std::fs::read_to_string(&cat)
            .unwrap()
            .contains("chr20:21000-21024")
    );
    let small = t.file("small", &format!(">chr20\n{}\n", "A".repeat(11000)));
    let fq = t.path("reads.gz");
    assert_eq!(
        run(smoke::cli, &["long-reads", &small, &small, &fq]).1,
        "4 reads\n"
    );
    use std::io::Read;
    let mut text = String::new();
    og_io::open(&fq).unwrap().read_to_string(&mut text).unwrap();
    let ls = text.lines().collect::<Vec<_>>();
    assert_eq!(ls.len(), 16);
    assert_eq!(ls[1], "A".repeat(10000));
    assert_eq!(ls[5], "T".repeat(10000));
    assert_eq!(ls[3], "?".repeat(10000));
    let bad = t.file("bad", ">x\nA\n");
    assert_eq!(
        run(
            smoke::cli,
            &["mini-ref", &bad, &orig, &mini, "450000", "3000"]
        )
        .0,
        2
    );
    assert_eq!(
        run(smoke::cli, &["eh-catalog", &small, &cat, "449000"]).0,
        2
    );
}

#[test]
fn pca_nonfinite_values_follow_python_comparisons() {
    let t = Temp::new();
    for (name, extra) in [
        ("nan", "n NaN 0 PopA\n"),
        ("inf", "n inf 0 PopA\n"),
        ("overflow", "n 1e400 0 PopA\n"),
    ] {
        let p = t.file(
            name,
            &format!("#h\n{extra}a -2 0 PopA\nb 2 0 PopB\np -1 0 Proj\n"),
        );
        assert_eq!(run(smoke::cli, &["check", "eigen", &p]).0, 1, "{name}");
    }
    let p = t.file("short.csv", "Type,Value\nSNP\n");
    assert_eq!(
        run(smoke::cli, &["csv-cell", &p, "Type=SNP", "Value"]).1,
        "None\n"
    );
    let p = t.file("quoted.tsv", "# a comment\n\"chr5\"\t1\t2\ta\tb\tc\t\t99\n");
    assert_eq!(run(smoke::cli, &["check", "parascopy", &p]).0, 1);
}

#[test]
fn smoke_coordinates_require_integers() {
    let t = Temp::new();
    for bad in ["-inf", "inf", "nan", "2000000.0"] {
        for text in [
            format!("h\nS 22 {bad} 2600000 4\n"),
            format!("h\nS 22 2000000 {bad} 4\n"),
        ] {
            let p = t.file("ibdmix", &text);
            assert_ne!(run(smoke::cli, &["check", "ibdmix", &p]).0, 0, "{text}");
        }
        for text in [
            format!(
                "#codes: POPA=0 POPB=1\n20 {bad} 10000000 0 0 0 0 0\n20 10000000 20000000 0 0 0 1 1\n"
            ),
            format!(
                "#codes: POPA=0 POPB=1\n20 0 {bad} 0 0 0 0 0\n20 10000000 20000000 0 0 0 1 1\n"
            ),
        ] {
            let p = t.file("rfmix", &text);
            assert_ne!(run(smoke::cli, &["check", "rfmix", &p]).0, 0, "{text}");
        }
    }
    let p = t.file("ints", "h\nS 22 +2_000_000 2_600_000 4.0\n");
    assert_eq!(run(smoke::cli, &["check", "ibdmix", &p]).0, 0);
}

#[test]
fn trace_retains_concrete_duplicate_exit() {
    let t = Temp::new();
    let first = "task_id\tname\tattempt\texit\tmemory\n1\tPRS (s)\t1\t-\t25769803776\n";
    let p = t.file(
        "trace",
        &format!("{first}1\tPRS (s)\t1\t137\t25769803776\n"),
    );
    let args = [
        "attempts",
        &p,
        "--name-glob",
        "PRS (s)",
        "--want",
        "1:137:24",
        "--nextflow-exit",
        "0",
    ];
    assert_eq!(run(nf_trace::cli, &args).0, 0);
    t.file(
        "trace",
        &format!("{first}1\tPRS (s)\t1\t137\t25769803776\n1\tPRS (s)\t1\t0\t25769803776\n"),
    );
    assert_eq!(run(nf_trace::cli, &args).0, 2);
}

#[test]
fn trace_integer_policy() {
    let t = Temp::new();
    let p = t.file("spaces", &TRACE.replace("\t8\t", "\t +8 \t"));
    let want = t.file(
        "want",
        "process\tcpus\tmemory_gb\ttime_h\nALIGN_MINIMAP2\t8\t32\t24\n",
    );
    let args = [
        "requests",
        &p,
        "--want",
        &want,
        "--executor-cpus",
        "16",
        "--executor-mem-gb",
        "96",
    ];
    assert_eq!(run(nf_trace::cli, &args).0, 0);
    t.file("spaces", &TRACE.replace("\t8\t", "\t٨\t"));
    assert_eq!(run(nf_trace::cli, &args).0, 2);
}

/// Run explicitly with a full git history and Python available. The oracle
/// is the actual replaced assertion, not a second maintained implementation.
#[test]
#[ignore = "requires archived source and a local Python interpreter"]
fn differential_coordinate_controls() {
    use std::process::Command;
    let t = Temp::new();
    let repo = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let archived = Command::new("git")
        .arg("-C")
        .arg(&repo)
        .args([
            "show",
            "669a8196cf4c377200f9b4af9c235ef85fb774a1:tests/smoke/commands.tsv",
        ])
        .output()
        .unwrap();
    assert!(archived.status.success());
    let table = String::from_utf8(archived.stdout).unwrap();
    for (kind, var, name) in [
        ("ibdmix", "IBDMIX_IMAGE", "o.txt"),
        ("rfmix", "RFMIX_IMAGE", "o.msp.tsv"),
    ] {
        let row = table
            .lines()
            .find(|l| l.starts_with(&format!("{var}\t")))
            .unwrap();
        let expect = row.split('\t').nth(4).unwrap();
        let code = expect
            .split("\" \"")
            .nth(1)
            .unwrap()
            .strip_suffix('"')
            .unwrap();
        for value in [
            "2000000",
            "+2_000_000",
            "٢٠٠٠٠٠٠",
            "𝟚𝟘𝟘𝟘𝟘𝟘𝟘",
            "-inf",
            "inf",
            "nan",
            "2000000.0",
        ] {
            for column in [0, 1] {
                let mut bounds = ["2000000", "2600000"];
                bounds[column] = value;
                let text = if kind == "ibdmix" {
                    format!("h\nS 22 {} {} 4\n", bounds[0], bounds[1])
                } else {
                    let a = if column == 0 { value } else { "0" };
                    let b = if column == 1 { value } else { "10000000" };
                    format!(
                        "#codes: POPA=0 POPB=1\n20 {a} {b} 0 0 0 0 0\n20 10000000 20000000 0 0 0 1 1\n"
                    )
                };
                let p = t.file(name, &text);
                let old = Command::new("python3")
                    .args(["-I", "-c", code])
                    .current_dir(&t.0)
                    .output()
                    .unwrap();
                let new = run(smoke::cli, &["check", kind, &p]);
                assert_eq!(
                    old.status.success(),
                    new.0 == 0,
                    "{kind} col={column} value={value}: {}",
                    new.2
                );
                println!(
                    "{kind} col={column} value={value}: Python={} Rust={}",
                    old.status.success(),
                    new.0
                );
            }
        }
    }
}

#[test]
#[ignore = "requires archived source, jq and a local Python interpreter"]
fn differential_json_controls() {
    use std::process::Command;
    let t = Temp::new();
    let repo = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let archived = Command::new("git")
        .arg("-C")
        .arg(&repo)
        .args([
            "show",
            "669a8196cf4c377200f9b4af9c235ef85fb774a1:tests/smoke/commands.tsv",
        ])
        .output()
        .unwrap();
    assert!(archived.status.success());
    let old_table = String::from_utf8(archived.stdout).unwrap();
    let new_table = std::fs::read_to_string(repo.join("tests/smoke/commands.tsv")).unwrap();
    let script = std::fs::read_to_string(repo.join("scripts/ci/image-smoke.sh")).unwrap();
    let helper = script
        .split("jqcheck() {")
        .nth(1)
        .unwrap()
        .split("\n}\n")
        .next()
        .unwrap();
    let helper = format!(
        "CHECK_FAILS=0; pass() {{ :; }}; fail() {{ CHECK_FAILS=1; }}; jqcheck() {{{helper}\n}}\njqcheck control \"$1\" \"$2\"; exit \"$CHECK_FAILS\""
    );
    for (var, index, file, template, good) in [
        (
            "FASTP_IMAGE",
            0,
            "fastp.json",
            "{\"summary\":{\"before_filtering\":{\"total_reads\":100},\"after_filtering\":{\"total_reads\":VALUE}}}",
            "90",
        ),
        (
            "EXPANSIONHUNTER_IMAGE",
            0,
            "bundled_catalog.json",
            "VALUE",
            "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]",
        ),
        (
            "EXPANSIONHUNTER_IMAGE",
            1,
            "eh.json",
            "{\"LocusResults\":{\"SMOKE1\":{\"Variants\":{\"SMOKE1\":{\"Genotype\":VALUE}}}}}",
            "\"12/12\"",
        ),
        (
            "PARAPHASE_IMAGE",
            0,
            "smn1.json",
            "{\"smn1\":{\"region_specific_info\":{\"smn1_cn\":2,\"smn2_cn\":VALUE,\"smn_del78_cn\":1}}}",
            "1",
        ),
        (
            "PYTHON_IMAGE",
            0,
            "HG002_cyp2d6.json",
            "{\"HG002_cyrius\":{\"Median_depth\":VALUE,\"Coverage_MAD\":0,\"Genotype\":\"test\"}}",
            "20",
        ),
    ] {
        let row = |table: &str| {
            table
                .lines()
                .find(|l| l.starts_with(&format!("{var}\t")))
                .unwrap()
                .split('\t')
                .nth(4)
                .unwrap()
                .to_owned()
        };
        let old = row(&old_table);
        let code = old
            .split("py \"")
            .nth(index + 1)
            .unwrap()
            .split("\" \"")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap();
        let new = row(&new_table);
        let filter = new
            .split("jqcheck \"")
            .nth(index + 1)
            .unwrap()
            .split('\'')
            .nth(1)
            .unwrap();
        for value in [
            good,
            "\"broken\"",
            "null",
            "[]",
            "{}",
            "true",
            "false",
            "20",
            "\"aaaaaaaaaaaaaaaaaaaaaaaaaaaa\"",
        ] {
            let json = template.replace("VALUE", value);
            let p = t.file(file, &json);
            let python = Command::new("python3")
                .args(["-I", "-c", code])
                .current_dir(&t.0)
                .output()
                .unwrap();
            let rust = Command::new("bash")
                .args(["-c", &helper, "_", &p, filter])
                .output()
                .unwrap();
            let deliberate = (var == "EXPANSIONHUNTER_IMAGE"
                && ((index == 0 && value.starts_with('"'))
                    || (index == 1 && ["true", "20"].contains(&value))))
                || (["PARAPHASE_IMAGE", "PYTHON_IMAGE"].contains(&var) && value == "true");
            if deliberate {
                assert!(!rust.status.success(), "typed policy: {var} {value}");
            } else {
                assert_eq!(
                    python.status.success(),
                    rust.status.success(),
                    "{var} {index} {value}: {}",
                    String::from_utf8_lossy(&python.stderr)
                );
            }
            println!(
                "{var}/{index} {value}: Python={} jq={} deliberate={deliberate}",
                python.status.success(),
                rust.status.success()
            );
        }
        let valid = template.replace("VALUE", good);
        let p = t.file(file, &format!("{valid}\n{valid}\n"));
        assert!(
            !Command::new("python3")
                .args(["-I", "-c", code])
                .current_dir(&t.0)
                .output()
                .unwrap()
                .status
                .success()
        );
        assert!(
            !Command::new("bash")
                .args(["-c", &helper, "_", &p, filter])
                .output()
                .unwrap()
                .status
                .success()
        );
    }
}
