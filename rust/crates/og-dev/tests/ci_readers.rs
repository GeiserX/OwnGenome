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
        assert_eq!(run(smoke::cli, &["check", kind, &p]).0, 0, "{kind}");
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
