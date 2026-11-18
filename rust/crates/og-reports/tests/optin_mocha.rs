//! `og mocha report` (step 51): mosaic chromosomal alterations from MoChA's
//! tables. Moved from tests/test_mocha_report.py (its cases 1 to 5 and 7);
//! the summary's reading is in og-summary's `tests/optin_summary.rs` and the
//! report lines in og-render's `tests/optin_report.rs`.
//!
//! Each case writes MoChA's `--stats` and `--calls` tables and the
//! per-chromosome site counts with chosen numbers, runs the verb and reads
//! what it writes.

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

const VERSIONS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../versions.env");
const STATS_HEAD: &str = "sample_id\tcomputed_gender\tcall_rate\tcov_median\tcov_sd\tcov_auto\tbaf_corr\tbaf_conc\tbaf_auto\tn_sites\tn_hets\n";
const CALLS_HEAD: &str = "sample_id\tcomputed_gender\tchrom\tbeg_GRCh38\tend_GRCh38\tlength\tp_arm\tq_arm\tn_sites\tn_hets\tn50_hets\tbdev\tbdev_se\trel_cov\trel_cov_se\tlod_lrr_baf\tlod_baf_phase\tn_flips\tbaf_conc\tlod_baf_conc\ttype\tcf\n";

/// One MoChA call; `call()` is a 20 Mb CN-LOH at a cell fraction of 0.15.
#[derive(Clone)]
struct Call {
    chrom: &'static str,
    beg: u64,
    end: u64,
    p: &'static str,
    q: &'static str,
    n50: u64,
    bdev: f64,
    rel: f64,
    typ: &'static str,
    cf: &'static str,
}

fn call() -> Call {
    Call {
        chrom: "chr20",
        beg: 60_000,
        end: 20_000_000,
        p: "T",
        q: "N",
        n50: 3000,
        bdev: 0.075,
        rel: 2.0,
        typ: "CN-LOH",
        cf: "0.15",
    }
}

impl Call {
    fn line(&self) -> String {
        format!(
            "S\tU\t{}\t{}\t{}\t{}\t{}\t{}\t20000\t6600\t{}\t{}\t0.001\t{:.1}\t0.01\tnan\t900.0\t0\tnan\tnan\t{}\t{}\n",
            self.chrom,
            self.beg,
            self.end,
            self.end - self.beg,
            self.p,
            self.q,
            self.n50,
            self.bdev,
            self.rel,
            self.typ,
            self.cf
        )
    }
}

/// The sample's numbers.
struct Sample {
    depth: f64,
    hets_per_10mb: u64,
    call_rate: &'static str,
    baf_auto: &'static str,
    phased: bool,
    stats: Option<&'static str>,
    sites: Option<&'static str>,
}

fn sample() -> Sample {
    Sample {
        depth: 30.0,
        hets_per_10mb: 3000,
        call_rate: "0.999",
        baf_auto: "0.001",
        phased: true,
        stats: None,
        sites: None,
    }
}

struct Run {
    code: u8,
    rows: Vec<HashMap<String, String>>,
    summary: HashMap<String, String>,
}

fn run(tag: &str, calls: &[Call], s: &Sample) -> Run {
    let d: PathBuf =
        std::env::temp_dir().join(format!("og-optin-mocha-{}-{tag}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    let p = |n: &str| d.join(n).to_str().unwrap().to_owned();
    let span: u64 = 60_000_000;
    let h = s.hets_per_10mb * span / 10_000_000;
    let stats = s.stats.map_or_else(
        || {
            format!(
                "{STATS_HEAD}S\tU\t{}\t{:?}\t3\t0.01\t0.001\t0.5\t{}\t100000\t{h}\n",
                s.call_rate, s.depth, s.baf_auto
            )
        },
        str::to_owned,
    );
    fs::write(p("stats.tsv"), stats).unwrap();
    fs::write(
        p("calls.tsv"),
        format!(
            "{CALLS_HEAD}{}",
            calls.iter().map(Call::line).collect::<String>()
        ),
    )
    .unwrap();
    let sites = s.sites.map_or_else(
        || {
            format!(
                "chr20\t1\t{span}\t{}\t{h}\t{}\n",
                3 * h,
                if s.phased { h } else { 0 }
            )
        },
        str::to_owned,
    );
    fs::write(p("sites.tsv"), sites).unwrap();
    let args: Vec<String> = [
        "--sample",
        "S",
        "--calls",
        &p("calls.tsv"),
        "--stats",
        &p("stats.tsv"),
        "--sites",
        &p("sites.tsv"),
        "--phasing",
        "input",
        "--out-calls",
        &p("S_mocha.tsv"),
        "--out-summary",
        &p("S_mocha_summary.tsv"),
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = og_reports::mocha::cli(Some("report"), &args, &mut o, &mut e);
    let mut r = Run {
        code,
        rows: Vec::new(),
        summary: HashMap::new(),
    };
    if code == 0 {
        let text = fs::read_to_string(p("S_mocha.tsv")).unwrap();
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
        r.rows = lines
            .map(|l| {
                head.iter()
                    .zip(l.split('\t'))
                    .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                    .collect()
            })
            .collect();
        r.summary = fs::read_to_string(p("S_mocha_summary.tsv"))
            .unwrap()
            .lines()
            .map(|l| {
                let (k, v) = l.split_once('\t').unwrap();
                (k.to_owned(), v.to_owned())
            })
            .collect();
    }
    let _ = fs::remove_dir_all(&d);
    r
}

#[test]
fn the_floor() {
    let s = run("floor30", &[call()], &sample()).summary;
    assert_eq!(
        s["floor_read"],
        format!("{:.3}", 6.0 / (30.0f64 * 3000.0).sqrt())
    );
    assert_eq!(
        [s["floor"].as_str(), &s["floor_min"], &s["floor_read"]],
        ["0.050", "0.05", "0.020"]
    );
    let thin = Sample {
        depth: 5.0,
        hets_per_10mb: 300,
        ..sample()
    };
    let s = run("floor5", &[call()], &thin).summary;
    assert_eq!(s["floor"], format!("{:.3}", 6.0 / (5.0f64 * 300.0).sqrt()));
    assert_eq!(s["het_per_10mb"], "300");
}

#[test]
fn reported_and_under_the_floor() {
    let r = run("reported", &[call()], &sample());
    let got: Vec<[&str; 3]> = r
        .rows
        .iter()
        .map(|x| [x["type"].as_str(), &x["cell_fraction"], &x["reported"]])
        .collect();
    assert_eq!(got, [["CN-LOH", "0.150", "yes"]]);
    assert_eq!(r.rows[0]["arms"], "20p, to the p telomere");
    assert_eq!(
        [
            r.summary["status"].as_str(),
            &r.summary["calls_reported"],
            &r.summary["calls_mocha"]
        ],
        ["assessed", "1", "1"]
    );
    let low = run(
        "low",
        &[Call {
            bdev: 0.02,
            cf: "0.04",
            ..call()
        }],
        &sample(),
    );
    assert_eq!(
        [low.rows[0]["reported"].as_str(), &low.rows[0]["reason"]],
        ["no", "cell fraction 0.040 under the floor 0.050"]
    );
    let thin = Sample {
        depth: 5.0,
        hets_per_10mb: 300,
        ..sample()
    };
    assert_eq!(run("thin", &[call()], &thin).rows[0]["reported"], "no");
}

#[test]
fn mochas_own_filter() {
    let r = run(
        "filter",
        &[
            Call {
                typ: "CNP_Loss",
                ..call()
            },
            Call {
                bdev: 0.2,
                n50: 500_000,
                cf: "0.4",
                ..call()
            },
            Call {
                rel: 3.0,
                bdev: 0.16,
                typ: "Gain",
                cf: "0.6",
                ..call()
            },
            Call {
                chrom: "chrX",
                ..call()
            },
            Call {
                cf: "nan",
                ..call()
            },
            Call {
                beg: 30_000_000,
                end: 60_000_000,
                p: "N",
                q: "T",
                typ: "Loss",
                bdev: 0.03,
                cf: "0.12",
                ..call()
            },
        ],
        &sample(),
    );
    let why: Vec<&str> = r.rows.iter().map(|x| x["reason"].as_str()).collect();
    assert!(
        why[0].starts_with("germline copy-number polymorphism"),
        "{why:?}"
    );
    assert!(why[1].starts_with("likely constitutional"), "{why:?}");
    assert!(why[2].starts_with("likely germline duplication"), "{why:?}");
    assert!(why[3].starts_with("not an autosome"), "{why:?}");
    assert_eq!(why[4], "cell fraction not estimated");
    let kept: Vec<(&str, &str)> = r
        .rows
        .iter()
        .filter(|x| x["reported"] == "yes")
        .map(|x| (x["type"].as_str(), x["arms"].as_str()))
        .collect();
    assert_eq!(kept, [("loss", "20q, to the q telomere")]);
    assert_eq!(
        [
            r.summary["calls_reported"].as_str(),
            &r.summary["calls_mocha"]
        ],
        ["1", "6"]
    );
}

#[test]
fn the_sample_filter() {
    for (tag, s, why) in [
        (
            "rate",
            Sample {
                call_rate: "0.95",
                ..sample()
            },
            "call rate",
        ),
        (
            "baf",
            Sample {
                baf_auto: "0.05",
                ..sample()
            },
            "phased BAF autocorrelation",
        ),
        (
            "unphased",
            Sample {
                phased: false,
                ..sample()
            },
            "no heterozygous site is phased",
        ),
        (
            "nanrate",
            Sample {
                call_rate: "nan",
                ..sample()
            },
            "MoChA gave no call rate or phased BAF autocorrelation",
        ),
        (
            "nanbaf",
            Sample {
                baf_auto: "nan",
                ..sample()
            },
            "MoChA gave no call rate or phased BAF autocorrelation",
        ),
    ] {
        let r = run(tag, &[call()], &s);
        let msg = &r.summary["message"];
        assert!(
            r.summary["status"] == "not assessed" && (msg.starts_with(why) || msg.contains(why)),
            "{why}: {:?}",
            r.summary
        );
        let rep: Vec<&str> = r.rows.iter().map(|x| x["reported"].as_str()).collect();
        assert_eq!(rep, ["no"], "{why}");
    }
}

#[test]
fn the_cli_refuses_tables_it_cannot_read() {
    let no_cov = Sample {
        stats: Some("sample_id\tcall_rate\nS\t0.99\n"),
        ..sample()
    };
    assert_eq!(run("nocov", &[call()], &no_cov).code, 2);
    let short = Sample {
        sites: Some("chr20\t1\t100\n"),
        ..sample()
    };
    assert_eq!(run("short", &[call()], &short).code, 2);
}

/// The value of `NAME="..."` in versions.env.
fn pin(env: &str, name: &str) -> String {
    env.lines()
        .find_map(|l| l.strip_prefix(&format!("{name}=\"")))
        .and_then(|v| v.split('"').next())
        .unwrap_or_default()
        .to_owned()
}

#[test]
fn versions_env_pins_mocha_as_the_e2e_case_does() {
    let env = fs::read_to_string(VERSIONS).unwrap();
    let image = pin(&env, "MOCHA_IMAGE");
    let (name, digest) = image.split_once('@').unwrap_or((&image, ""));
    assert_eq!(
        (
            name.split('/').next().unwrap(),
            name.rsplit(':').next().unwrap()
        ),
        ("us.gcr.io", "1.22-20250819")
    );
    assert_eq!(
        digest,
        "sha256:c9398385cfe3a968918ddf470e2769b7ee42cbc0bb8f564b4744c460343550a5"
    );
    assert_eq!(
        pin(&env, "MOCHA_COMMIT"),
        "95686b7b65f53a490513be76bb120e5fc20a8bcf"
    );
    assert_eq!(pin(&env, "MOCHA_VERSION"), "2025-08-19");
}
