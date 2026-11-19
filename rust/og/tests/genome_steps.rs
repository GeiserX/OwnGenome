//! The genome's copy-number, repeat and traits results from the step that
//! writes them to the report: the checks of the old `tests/test_copy_number.py`,
//! `tests/test_lpa_kiv2.py` and `tests/test_traits_report.py` that need more
//! than one tool (og cn, og traits, og dev, og couple, og summary collect, og
//! report render), so they live in the crate that holds every verb.
//!
//! - The wiring: steps 35, 41 and 46 and their Nextflow modules call `og cn`
//!   with the bedcov passes of their catalog row, and never the old Python.
//! - Step 35's dosage, CYP21A2 and interval blocks, cut out of
//!   `scripts/35-paralogs.sh` and run with `run_in` stubbed (og runs for
//!   real; samtools prints the fixtures of `tests/fixtures/cn/`), write the
//!   known tables, and a failed pass writes not assessed.
//! - `og dev build cn-intervals` on a small ClinGen BED.
//! - The summary reads step 35's interval and CYP21A2 tables, the couple
//!   screen reads the CYP21A2 tables, and the reports print them.
//! - Step 46's Lp(a) count and KIV-2 4925G>A in both reports.
//! - Step 42's rules in the HTML report, and a refused traits table leaves
//!   the outputs the collector reads.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use og_cn::catalog::{COLUMNS, Catalog, Locus};
use og_summary::collect::{Group, Input};
use pycompat::json::{self, Dump, Value};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

fn repo(rel: &str) -> String {
    format!("{REPO}/{rel}")
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_default()
}

fn cat() -> Catalog {
    Catalog::load(&repo("bin/cn_loci.tsv")).unwrap()
}

/// A scratch folder of its own, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!("og-genome-steps-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn path(&self, rel: &str) -> String {
        self.0.join(rel).to_string_lossy().into_owned()
    }

    fn put(&self, rel: &str, text: &str) -> String {
        let p = self.0.join(rel);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn dir(&self, rel: &str) -> String {
        let p = self.0.join(rel);
        std::fs::create_dir_all(&p).unwrap();
        p.to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A copy of og beside the catalog, as `bin/.og/og-*` sits in a checkout:
/// `og cn` finds `cn_loci.tsv` beside its own binary.
fn og() -> &'static str {
    static OG: OnceLock<String> = OnceLock::new();
    OG.get_or_init(|| {
        let d = std::env::temp_dir().join(format!("og-genome-steps-bin-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        std::fs::copy(env!("CARGO_BIN_EXE_og"), d.join("og")).unwrap();
        std::fs::copy(repo("bin/cn_loci.tsv"), d.join("cn_loci.tsv")).unwrap();
        d.join("og").to_string_lossy().into_owned()
    })
}

fn og_cn(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_cn::cli::main(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

/// A `metric<TAB>value` table, in its order.
fn table(path: &str) -> Vec<(String, String)> {
    read(path)
        .lines()
        .skip(1)
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_owned(), v.to_owned()))
        .collect()
}

fn get<'a>(t: &'a [(String, String)], k: &str) -> &'a str {
    t.iter()
        .find(|(m, _)| m == k)
        .map_or("", |(_, v)| v.as_str())
}

/// The whole summary.json og summary collect writes for sample S in `d`.
fn summary(d: &str) -> Value {
    og_summary::collect::post::collect(&Input {
        sample: "S".into(),
        sample_dir: d.into(),
        declared_sex: None,
        assets: repo("assets"),
        group: Group::All,
        tables: repo("bin"),
    })
    .unwrap()
}

fn at<'a>(v: &'a Value, path: &[&str]) -> Option<&'a Value> {
    path.iter().try_fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k),
        _ => None,
    })
}

fn s(v: Option<&Value>) -> &str {
    match v {
        Some(Value::Str(s)) => s,
        _ => "",
    }
}

/// Both reports of a summary: (text, HTML).
fn reports(v: &Value) -> (String, String) {
    let text = json::dumps(v, &Dump::default()).unwrap();
    let sm = og_render::model::Summary::parse(&text).unwrap();
    (
        og_render::text::report(&sm).unwrap(),
        og_render::html::report(&sm).unwrap(),
    )
}

fn e(t: &str) -> String {
    og_render::escape::e(t).into_owned()
}

// The wiring

/// The option strings of each `samtools bedcov` command in `text` (options,
/// then a path, a quoted argument or a variable).
fn bedcov_calls(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    for (i, _) in text.match_indices("samtools bedcov") {
        let rest = &text[i + "samtools bedcov".len()..];
        let mut opts = String::new();
        let mut r = rest;
        loop {
            let b = r.as_bytes();
            if b.len() >= 3 && b[0] == b' ' && b[1] == b'-' && b[2].is_ascii_alphabetic() {
                opts += &r[..3];
                r = &r[3..];
                let digits: usize = r
                    .strip_prefix(' ')
                    .map_or(0, |x| x.bytes().take_while(u8::is_ascii_digit).count());
                if digits > 0 && !r[1 + digits..].starts_with(|c: char| !c.is_whitespace()) {
                    opts += &r[..1 + digits];
                    r = &r[1 + digits..];
                }
                continue;
            }
            break;
        }
        let Some(arg) = r.strip_prefix(' ') else {
            continue;
        };
        let arg: &str = arg.split_whitespace().next().unwrap_or("");
        let is_bed = arg.find(".bed").is_some_and(|j| {
            j > 0 && !arg[j + 4..].starts_with(|c: char| c.is_alphanumeric() || c == '_')
        });
        if arg.starts_with(['"', '$', '/']) || is_bed {
            out.push(opts.trim().to_owned());
        }
    }
    out
}

#[test]
fn steps_35_41_46_and_their_modules_call_og_cn_with_the_catalogs_passes() {
    let c = cat();
    for (path, loc) in [
        ("scripts/35-paralogs.sh", "HBA"),
        ("modules/local/parascopy/main.nf", "HBA"),
        ("scripts/41-bloodtype.sh", "RHD"),
        ("modules/local/bloodtype/main.nf", "RHD"),
        ("scripts/46-lpa.sh", "LPA_KIV2"),
        ("modules/local/lpa/main.nf", "LPA_KIV2"),
    ] {
        let text = read(&repo(path));
        assert!(
            text.contains(&format!("og cn bed --locus {loc}")),
            "{path}: og cn bed --locus {loc}"
        );
        if loc != "RHD" {
            assert!(
                text.contains(&format!("og cn call --locus {loc}")),
                "{path}: og cn call --locus {loc}"
            );
        }
        let mut passes: Vec<String> = c.loci[loc]
            .depth()
            .passes
            .iter()
            .map(|p| if p == "." { String::new() } else { p.clone() })
            .collect();
        if loc == "HBA" {
            // and the dosage loop's two passes, the interval pass, the CYP21A2 passes
            passes.extend(["", "-Q 1", "", "", "-Q 1"].map(str::to_owned));
        }
        let mut calls = bedcov_calls(&text);
        calls.sort();
        passes.sort();
        assert_eq!(calls, passes, "{path}: its bedcov passes are the catalog's");
        for old in [
            "copy_number.py ",
            "hba_copy_number.py ",
            "lpa_kiv2.py ",
            "traits_report.py --marker",
        ] {
            assert!(!text.contains(old), "{path} calls the old Python: {old}");
        }
    }
}

const DOSAGE: [&str; 3] = ["HP", "LCE3BC", "CCL3L1"];

#[test]
fn step_35_and_the_parascopy_module_run_the_dosage_interval_and_cyp21a2_checks() {
    let (step, module) = (
        read(&repo("scripts/35-paralogs.sh")),
        read(&repo("modules/local/parascopy/main.nf")),
    );
    let list = |text: &str, key: &str| -> Option<Vec<String>> {
        text.lines()
            .find_map(|l| l.trim_start().strip_prefix(key))
            .and_then(|v| v.strip_prefix('"')?.strip_suffix('"'))
            .map(|v| v.split_whitespace().map(str::to_owned).collect())
    };
    assert_eq!(
        list(&step, "DOSAGE_LOCI=").unwrap_or_default(),
        DOSAGE,
        "step 35: DOSAGE_LOCI"
    );
    assert!(step.lines().any(|l| l.starts_with("DOSAGE_LOCI=")));
    assert_eq!(
        list(&module, "dosage_loci=").unwrap_or_default(),
        DOSAGE,
        "the module's dosage_loci"
    );
    for (path, text, q) in [("step 35", &step, "\""), ("the module", &module, "\"\\")] {
        for cmd in ["bed", "call"] {
            assert!(
                text.contains(&format!("og cn {cmd} --locus {q}$loc")),
                "{path}: og cn {cmd} on each dosage row"
            );
        }
        for cmd in [
            "og cn bed --kind interval",
            "og cn intervals --bedcov",
            "og cn intervals --failed",
        ] {
            assert!(text.contains(cmd), "{path}: {cmd}");
        }
    }
    assert!(module.contains("def sex = meta.sex in ['male', 'female'] ? meta.sex : 'unknown'"));
    assert!(module.contains("--sex ${sex}"));
    for (path, text, q1) in [
        ("step 35", &step, "bedcov_q1.tsv\")\""),
        ("the module", &module, "cyp21a2_bedcov_q1.tsv"),
    ] {
        for cmd in [
            "og cn bed --locus CYP21A2",
            "og cn bed --psv --locus CYP21A2",
            "--bedcov-mapq1",
            q1,
            "samtools mpileup -A -B -q 0 -Q 13 -d 100000",
            "og cn call --locus CYP21A2 --bedcov",
            "--pileup",
            "_cyp21a2_psv.tsv",
        ] {
            assert!(text.contains(cmd), "{path}: {cmd}");
        }
    }
    assert!(module.contains("path(\"${meta.id}_cyp21a2_psv.tsv\"),     emit: cyp21a2"));
}

/// The lines of step 35 from `start` to the line `end`, both included.
fn block(start: &str, end: &str) -> String {
    let step = read(&repo("scripts/35-paralogs.sh"));
    let i = step
        .find(&format!("\n{start}"))
        .unwrap_or_else(|| panic!("step 35 holds the {start} block"))
        + 1;
    let j = step[i..]
        .find(&format!("\n{end}\n"))
        .unwrap_or_else(|| panic!("the end of the {start} block"))
        + i;
    step[i..j + end.len() + 2].to_owned()
}

/// `run_in` and `cpath` for a block of step 35: og runs for real, the
/// `samtools` lines are the case's.
fn stub(samtools: &str) -> String {
    format!(
        r#"run_in() {{
  while [ $# -gt 0 ]; do case $1 in -v|--cpus|--memory) shift 2 ;; *) break ;; esac; done
  shift
  case "$1 $2" in
    "/pgp-bin/og cn") shift; "${{OG}}" "$@" ;;
{samtools}
    *) return 1 ;;
  esac
}}
cpath() {{ printf '%s' "$1"; }}
"#
    )
}

fn run_block(script: &str, env: &[(&str, &str)]) -> (i32, String) {
    let o = Command::new("bash")
        .args(["-euo", "pipefail", "-c", script])
        .env("OG", og())
        .env("REPO", REPO)
        .env("PGP_ROOT", REPO)
        .env("SAMPLE", "S")
        .env("PARASCOPY_IMAGE", "img")
        .env("FIX", repo("tests/fixtures/cn"))
        .envs(env.iter().copied())
        .output()
        .unwrap();
    (
        o.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

#[test]
fn step_35s_dosage_block_counts_the_fixtures() {
    let w = Scratch::new("dosage");
    let script = stub(
        r#"    "samtools bedcov") shift 2; q=""; if [ "$1" = -Q ]; then q=.q1; shift 2; fi
      d=$(dirname "$1"); cat "${FIX}/${d##*/}${q}.bedcov.tsv" ;;"#,
    ) + &block("DOSAGE_LOCI=", "rm -rf \"$DOSAGE_WORK\"");
    let out = w.path("");
    let (rc, err) = run_block(&script, &[("OUTDIR", &out), ("BAM", &w.path("S.bam"))]);
    let text = read(&w.path("S_cn_dosage.tsv"));
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some("locus\tstatus\tcopies\testimate\treason\treading")
    );
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    assert_eq!(rc, 0, "{err}");
    assert_eq!(rows.iter().map(|r| r[0]).collect::<Vec<_>>(), DOSAGE);
    let c = cat();
    for (loc, n) in [("HP", 3), ("LCE3BC", 1), ("CCL3L1", 5)] {
        let r = rows.iter().find(|r| r[0] == loc).unwrap();
        assert_eq!(r[1..3], ["called", &n.to_string()], "{loc}");
        assert_eq!(
            r[5..],
            [c.loci[loc].reading_count(n).as_str()],
            "{loc}: its reading"
        );
    }
}

/// samtools bedcov for the interval block: 30x, the segments named in
/// FAKE_HALF at 15x, those in FAKE_ZERO at 0; FAKE_FAIL fails the pass.
const FAKECOV: &str = r#"    "samtools bedcov") shift 2; [ -z "${FAKE_FAIL:-}" ] || return 1
      awk -F'\t' -v OFS='\t' -v half="${FAKE_HALF:-}" -v zero="${FAKE_ZERO:-}" '
        BEGIN { n = split(half, h, ","); for (i = 1; i <= n; i++) H[h[i]] = 1
                n = split(zero, z, ","); for (i = 1; i <= n; i++) Z[z[i]] = 1 }
        { d = ($4 in Z) ? 0 : ($4 in H) ? 15 : 30; print $0, int(d * ($3 - $2)) }' "$1" ;;"#;

fn interval_rows(c: &Catalog) -> Vec<&Locus> {
    c.loci.values().filter(|r| r.is_interval_row()).collect()
}

fn smn_header() -> &'static str {
    "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n"
}

#[test]
fn step_35s_interval_block_and_the_summary_of_its_table() {
    let (c, w) = (cat(), Scratch::new("intervals"));
    let names: Vec<&str> = interval_rows(&c).iter().map(|r| r.id.as_str()).collect();
    let half: Vec<&str> = c.loci["AZFc_DAZ"]
        .targets
        .iter()
        .map(|s| s.name.as_str())
        .collect();
    let half = half.join(",");
    let script = stub(FAKECOV) + &block("INTERVAL_OUT=", "rm -rf \"$INTERVAL_WORK\"");
    let out = w.dir("out");
    let env = [
        ("OUTDIR", out.as_str()),
        ("BAM", "S.bam"),
        ("SEX", "male"),
        ("FAKE_HALF", half.as_str()),
    ];
    let (rc, err) = run_block(&script, &env);
    let tab = format!("{out}/S_cn_intervals.tsv");
    let text = read(&tab);
    let rows: HashMap<&str, Vec<&str>> = text
        .lines()
        .skip(1)
        .map(|l| (l.split('\t').next().unwrap(), l.split('\t').collect()))
        .collect();
    let order: Vec<&str> = text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(rc, 0, "{err}");
    assert_eq!(order, names);
    // AZFc's DAZ at 15x against controls at 30x on a male: 4 x 15 / 30 = 2 of 4.
    assert_eq!(rows["AZFc_DAZ"][6..9], ["2", "4", "deletion"]);
    assert_eq!(rows["AZFa"][6..9], ["1", "1", "no change"]);
    let good = w.put("good.tsv", &text);
    let mut env_fail = env.to_vec();
    env_fail.push(("FAKE_FAIL", "1"));
    let (rc, err) = run_block(&script, &env_fail);
    let failed: Vec<String> = read(&tab).lines().skip(1).map(str::to_owned).collect();
    assert_eq!(rc, 0, "{err}");
    assert_eq!(failed.len(), names.len());
    assert!(
        failed
            .iter()
            .all(|l| l.split('\t').nth(10) == Some("check_failed")),
        "a failed bedcov"
    );
    // og summary collect: step 35's section carries the rows with a change.
    let sd = w.dir("sample");
    w.put("sample/paralogs/S_smn_copy_number.tsv", smn_header());
    std::fs::copy(&good, format!("{sd}/paralogs/S_cn_intervals.tsv")).unwrap();
    for pass in ["with the SMN table", "without another table"] {
        let v = summary(&sd);
        let ci = at(&v, &["sections", "parascopy", "data", "cn_intervals"])
            .unwrap_or_else(|| panic!("{pass}: {v:?}"));
        let changes: Vec<&str> = match at(ci, &["changes"]) {
            Some(Value::Array(a)) => a.iter().map(|x| s(at(x, &["locus"]))).collect(),
            _ => Vec::new(),
        };
        assert_eq!(changes, ["AZFc_DAZ"], "{pass}");
        assert!(
            matches!(at(ci, &["rows"]), Some(Value::Int(n)) if *n as usize == names.len()),
            "{pass}"
        );
        std::fs::remove_file(format!("{sd}/paralogs/S_smn_copy_number.tsv")).unwrap_or_default();
    }
}

/// samtools for the CYP21A2 block: the fixtures of CYP_CASE; the MAPQ >= 1
/// pass at 90% of each depth, or all of it with CYP_NO_MAPQ0; mpileup fails
/// with MPILEUP_FAILS.
const STUB_CYP: &str = r#"    "samtools bedcov")
      if [ "${3:-}" = "-Q" ] && [ -z "${CYP_NO_MAPQ0:-}" ]; then
        awk -F'\t' -v OFS='\t' '{$5 = int($5 * 0.9); print}' "${FIX}/CYP21A2_${CYP_CASE}.bedcov.tsv"
      else
        cat "${FIX}/CYP21A2_${CYP_CASE}.bedcov.tsv"
      fi ;;
    "samtools mpileup") [ -z "${MPILEUP_FAILS:-}" ] && cat "${FIX}/CYP21A2_${CYP_CASE}.pileup.tsv" ;;"#;

/// Step 35's CYP21A2 block on each case; the tables by case.
fn cyp_tables(w: &Scratch) -> HashMap<&'static str, String> {
    let script = stub(STUB_CYP) + &block("CYP21_OUT=", "rm -rf \"$CYP21_WORK\"");
    let out = w.dir("cyp");
    let mut tabs = HashMap::new();
    for (case, fixture, extra) in [
        ("normal", "normal", None),
        ("deletion", "deletion", None),
        ("conversion", "conversion", None),
        ("fail", "normal", Some(("MPILEUP_FAILS", "1"))),
        ("no_mapq0", "deletion", Some(("CYP_NO_MAPQ0", "1"))),
    ] {
        let mut env = vec![
            ("OUTDIR", out.as_str()),
            ("BAM", "S.bam"),
            ("CYP_CASE", fixture),
        ];
        env.extend(extra);
        let (rc, err) = run_block(&script, &env);
        let tab = format!("{out}/S_cyp21a2_psv.tsv");
        assert!(rc == 0 && !table(&tab).is_empty(), "{case}: {rc} {err}");
        let keep = w.path(&format!("{case}.tsv"));
        std::fs::rename(&tab, &keep).unwrap();
        tabs.insert(case, keep);
    }
    tabs
}

#[test]
fn step_35s_cyp21a2_block_the_summary_the_report_and_the_couple_screen() {
    let w = Scratch::new("cyp21a2");
    let tabs = cyp_tables(&w);
    let row = |case: &str| {
        let t = table(&tabs[case]);
        ["status", "cyp21a2_copies", "module_total", "finding"].map(|k| get(&t, k).to_owned())
    };
    assert_eq!(row("normal"), ["called", "2", "4", "none"]);
    assert_eq!(row("deletion"), ["called", "1", "3", "deletion"]);
    assert_eq!(row("conversion"), ["called", "1", "4", "conversion"]);
    assert_eq!(row("fail"), ["not assessed", ".", ".", "."]);
    assert_eq!(row("no_mapq0"), ["not assessed", ".", ".", "."]);
    assert_eq!(
        read(&tabs["fail"]).matches("reason\tcheck_failed").count(),
        1
    );
    let (nm, nl, fl) = (
        table(&tabs["no_mapq0"]),
        table(&tabs["normal"]),
        table(&tabs["fail"]),
    );
    assert_eq!(
        (get(&nm, "reason"), get(&nm, "mapq0_share")),
        ("low_mapq0_share", "0.000")
    );
    assert!(
        get(&nm, "message")
            .starts_with("not assessed: the BAM has few or no multi-mapping reads here")
    );
    assert_eq!(
        (get(&nl, "mapq0_share"), get(&nl, "message")),
        ("0.100", "."),
        "the MAPQ >= 1 pass ran"
    );
    assert_eq!(
        fl.iter().map(|x| x.0.as_str()).collect::<Vec<_>>(),
        nl.iter()
            .map(|x| x.0.as_str())
            .filter(|k| !k.starts_with("psv_"))
            .collect::<Vec<_>>(),
        "a failed mpileup writes the rows of a called table but the sites'"
    );

    // the couple screen reads each table
    let ped = w.put(
        "ped.tsv",
        "sample\tsex\tfather\tmother\tpartner\nA\tfemale\t.\t.\tB\nB\tmale\t.\t.\tA\n",
    );
    let vcf = w.put(
        "none.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
    );
    let carrier = repo("assets/carrier_genes.tsv");
    let couple = |name: &str,
                  extra: &[String],
                  genes: &str|
     -> (u8, Option<HashMap<String, String>>, String) {
        let od = w.path(&format!("couple-{name}"));
        let mut argv: Vec<String> = [
            "--pedigree",
            &ped,
            "--genes",
            genes,
            "--clinvar",
            &format!("A={vcf}"),
            "--clinvar",
            &format!("B={vcf}"),
            "--out-dir",
            &od,
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
        argv.extend(extra.iter().cloned());
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let rc = og_couple::screen::cli(&argv, &mut out, &mut err);
        if rc != 0 {
            return (rc, None, String::from_utf8_lossy(&err).into_owned());
        }
        let text = read(&format!("{od}/couple_A_B.tsv"));
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let row = lines.find(|l| l.starts_with("CYP21A2\t")).map(|l| {
            head.iter()
                .map(|h| (*h).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        });
        (0, row, read(&format!("{od}/couple_A_B.md")))
    };
    let arg = |s: &str, case: &str| vec!["--cyp21a2".to_owned(), format!("{s}={}", tabs[case])];
    let (rc, r, md) = couple(
        "deletion",
        &[arg("A", "deletion"), arg("B", "deletion")].concat(),
        &carrier,
    );
    let r = r.unwrap_or_else(|| panic!("{rc} {md}"));
    assert_eq!(
        (r["result"].as_str(), r["offspring_risk"].as_str()),
        ("shared", "1 in 4 per pregnancy")
    );
    assert!(r["findings_a"].contains("chimeric deletion"), "{r:?}");
    let (_, r, _) = couple("conversion", &arg("A", "conversion"), &carrier);
    assert!(r.unwrap()["findings_a"].contains("possible CYP21A2 conversion at In2G"));
    let gap = "CYP21A2 deletion and conversions of A: ";
    let (rc, r, md) = couple("normal", &arg("A", "normal"), &carrier);
    assert!(
        rc == 0 && r.is_none() && !md.contains(gap),
        "the unmodified table: {md}"
    );
    let (rc, r, md) = couple("no_mapq0", &arg("A", "no_mapq0"), &carrier);
    assert!(rc == 0 && r.is_none(), "{md}");
    assert!(
        md.contains(&format!(
            "{gap}CYP21A2 not assessed: the BAM has few or no multi-mapping reads"
        )),
        "{md}"
    );
    let (rc, r, md) = couple("fail", &arg("A", "fail"), &carrier);
    assert!(
        rc == 0 && r.is_none() && md.contains(&format!("{gap}CYP21A2 not assessed")),
        "{md}"
    );
    let bad = w.put(
        "bad.tsv",
        &read(&tabs["normal"]).replace("finding\tnone", "finding\tlost"),
    );
    let (rc, _, md) = couple(
        "bad",
        &["--cyp21a2".to_owned(), format!("A={bad}")],
        &carrier,
    );
    assert_eq!(rc, 2, "an unknown finding stops the screen: {md}");
    let no_cyp: String = read(&carrier)
        .lines()
        .filter(|l| !l.starts_with("CYP21A2\t"))
        .map(|l| l.to_owned() + "\n")
        .collect();
    let no_cyp = w.put("no_cyp21a2.tsv", &no_cyp);
    let (rc, r, _) = couple("no-list", &arg("A", "deletion"), &no_cyp);
    assert!(
        rc == 0 && r.is_none(),
        "a CYP21A2 finding needs CYP21A2 in the carrier list"
    );
    for (desc, extra, want) in [
        ("no partner has one", Vec::new(), HashSet::from(["A", "B"])),
        (
            "one partner has one",
            arg("A", "normal"),
            HashSet::from(["B"]),
        ),
    ] {
        let (_, _, md) = couple(&format!("gap-{}", extra.len()), &extra, &carrier);
        let got: HashSet<&str> = ["A", "B"]
            .into_iter()
            .filter(|p| {
                md.contains(&format!(
                    "CYP21A2 deletion and conversions of {p}: step 35's CYP21A2 table was not passed to the couple screen (--cyp21a2)."
                ))
            })
            .collect();
        assert_eq!(got, want, "{desc}");
    }

    // og summary collect carries the table, and the reports say why one is not assessed
    let sd = w.dir("sample");
    w.put("sample/paralogs/S_smn_copy_number.tsv", smn_header());
    std::fs::copy(
        &tabs["deletion"],
        format!("{sd}/paralogs/S_cyp21a2_psv.tsv"),
    )
    .unwrap();
    let v = summary(&sd);
    let cy = at(&v, &["sections", "parascopy", "data", "cyp21a2"]);
    assert_eq!(
        (
            s(cy.and_then(|c| at(c, &["finding"]))),
            s(cy.and_then(|c| at(c, &["cyp21a2_copies"])))
        ),
        ("deletion", "1")
    );
    std::fs::copy(
        &tabs["no_mapq0"],
        format!("{sd}/paralogs/S_cyp21a2_psv.tsv"),
    )
    .unwrap();
    let v = summary(&sd);
    let line = og_render::sections::cyp21a2::cyp21a2_line(at(
        &v,
        &["sections", "parascopy", "data", "cyp21a2"],
    ))
    .unwrap();
    assert!(
        line.contains(": not assessed (the BAM has few or no multi-mapping reads here")
            && !line.contains("deletion,"),
        "{line}"
    );
    let (txt, _) = reports(&v);
    assert!(txt.contains(&line), "the text report prints the line");
}

// og dev build cn-intervals

const SMALL_BED: &str = "track name=\"test\"
#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb
chr22\t18159724\t18924717\tBP-A\t0\t.\t18159724\t18924717\t0,0,0
chr22\t18924718\t21111383\t22q11.2 recurrent (DGS) region (proximal, A-D) (includes TBX1)\t0\t.\t18924718\t21111383\t255,111,0
chr22\t21111384\t21562827\tBP-D\t0\t.\t21111384\t21562827\t0,0,0
chrx\t6527561\t6537770\tXp22.31 distal breakpoint\t0\t.\t6527561\t6537770\t0,0,0
chrx\t6537771\t8156913\tXp22.31 recurrent region (includes STS)\t0\t.\t6537771\t8156913\t255,111,0
";

/// og dev build cn-intervals --stdout on a ClinGen BED and no masks: (exit
/// code, rows by column, stderr).
fn build_rows(w: &Scratch, bed: &str) -> (i32, Vec<HashMap<String, String>>, String) {
    let p = w.put("clingen.bed", bed);
    let o = Command::new(og())
        .args(["dev", "build", "cn-intervals", "--clingen", &p])
        .args([
            "--segdups",
            "/dev/null",
            "--gaps",
            "/dev/null",
            "--satellites",
            "/dev/null",
            "--parascopy",
            "/dev/null",
        ])
        .arg("--stdout")
        .output()
        .unwrap();
    let rows = String::from_utf8_lossy(&o.stdout)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            COLUMNS
                .iter()
                .map(|c| (*c).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect();
    (
        o.status.code().unwrap_or(-1),
        rows,
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

fn load(w: &Scratch, rows: &[HashMap<String, String>]) -> Result<usize, String> {
    let mut text = COLUMNS.join("\t") + "\n";
    for r in rows {
        text += &(COLUMNS
            .iter()
            .map(|c| r[*c].as_str())
            .collect::<Vec<_>>()
            .join("\t")
            + "\n");
    }
    let p = w.put("c.tsv", &text);
    Catalog::load(&p).map(|c| c.loci.len()).map_err(|e| e.0)
}

#[test]
fn og_dev_build_cn_intervals_on_a_small_clingen_bed() {
    let w = Scratch::new("build");
    let (code, rows, err) = build_rows(&w, SMALL_BED);
    let firsts: Vec<(&str, &str)> = rows
        .iter()
        .take(2)
        .map(|r| (r["locus"].as_str(), r["chrom"].as_str()))
        .collect();
    assert_eq!(code, 0, "{err}");
    assert_eq!(
        firsts,
        [("22q11.2_DGS_proximal_A-D", "chr22"), ("Xp22.31", "chrX")]
    );
    assert_eq!(load(&w, &rows), Ok(rows.len()), "its rows load");
    let lower: Vec<HashMap<String, String>> = rows
        .iter()
        .map(|r| {
            let mut r = r.clone();
            if r["chrom"] == "chrX" {
                r.insert("chrom".into(), "chrx".into());
            }
            r
        })
        .collect();
    assert!(
        load(&w, &lower).is_err_and(|e| e.contains("chrom 'chrx'")),
        "the rows left as chrx are refused"
    );
    let (code, _, err) = build_rows(&w, &SMALL_BED.replace("chrx", "chrQ"));
    assert!(
        code == 2 && err.contains("'chrQ' is not a GRCh38 chromosome"),
        "{code} {err}"
    );
    let r = &rows[0];
    let line: Vec<&str> = COLUMNS.iter().map(|c| r[*c].as_str()).collect();
    let loc = og_cn::catalog::parse_row(&line).unwrap();
    for t in &loc.targets {
        assert!(
            !(t.start < 18_924_717 && 18_159_724 < t.end
                || t.start < 21_562_827 && 21_111_384 < t.end),
            "{t:?}"
        );
    }
}

// Step 46's Lp(a) section in both reports

const TRAITS_COLUMNS: [&str; 11] = [
    "rsid", "gene", "chrom", "pos", "trait", "genotype", "status", "reading", "ancestry", "note",
    "pmids",
];
const R1: &str = "rs10455872";
const R2: &str = "rs3798220";
const ADVICE: &str = "A research estimate from read depth, not a diagnostic test. Lp(a) is measured in blood; a 2022 \
                      European Atherosclerosis Society consensus advises measuring it at least once in every adult.";
const COUNT38: &str = "KIV-2 repeats (both chromosomes): 38.0";

/// Writes step 46's table for `n` copies at a flank depth into `d/lpa/`.
fn lpa_table(d: &str, n: f64, flank: f64) {
    let r = cat().loci["LPA_KIV2"].clone();
    let k = flank * (n + 2.0) / 14.0;
    let mut m: HashMap<&str, f64> = HashMap::new();
    for (i, s) in r.targets.iter().enumerate() {
        m.insert(&s.name, k * if i % 2 == 0 { 1.03 } else { 0.97 });
    }
    for (i, s) in r.controls.iter().enumerate() {
        m.insert(&s.name, flank * if i % 2 == 0 { 1.03 } else { 0.97 });
    }
    let bed = |share: f64| -> String {
        r.regions()
            .iter()
            .map(|(c, s, e, nm)| {
                let v = m[nm]
                    * if nm.starts_with("kiv2_") {
                        1.0 - share
                    } else {
                        1.0
                    };
                format!(
                    "{c}\t{s}\t{e}\t{nm}\t{}\n",
                    (v * (e - s) as f64).round_ties_even() as i64
                )
            })
            .collect()
    };
    std::fs::create_dir_all(format!("{d}/lpa")).unwrap();
    let (a, b) = (format!("{d}/a.tsv"), format!("{d}/b.tsv"));
    std::fs::write(&a, bed(0.0)).unwrap();
    std::fs::write(&b, bed(0.9)).unwrap();
    let cp = repo("bin/cn_loci.tsv");
    let out = format!("{d}/lpa/S_lpa_kiv2.tsv");
    let (rc, _, err) = og_cn(&[
        "call",
        "--catalog",
        &cp,
        "--locus",
        "LPA_KIV2",
        "--bedcov",
        &a,
        "--bedcov-mapq1",
        &b,
        "--table",
        &out,
    ]);
    assert_eq!(rc, 0, "{err}");
    std::fs::remove_file(a).unwrap();
    std::fs::remove_file(b).unwrap();
}

/// A sample with step 46's table and, with `genotypes`, step 42's rows for
/// the two LPA variants: (summary, text, HTML).
fn lpa_sample(
    w: &Scratch,
    name: &str,
    n: f64,
    flank: f64,
    genotypes: &[(&str, &str)],
) -> (Value, String, String) {
    let d = w.dir(&format!("{name}/S"));
    lpa_table(&d, n, flank);
    if !genotypes.is_empty() {
        let mut text = TRAITS_COLUMNS.join("\t") + "\n";
        for (rs, gt) in genotypes {
            let status = if *gt == "." { "not_covered" } else { "called" };
            text += &TRAITS_COLUMNS
                .map(|k| match k {
                    "rsid" => *rs,
                    "gene" => "LPA",
                    "genotype" => gt,
                    "status" => status,
                    _ => ".",
                })
                .join("\t");
            text += "\n";
        }
        w.put(&format!("{name}/S/traits/S_traits.tsv"), &text);
    }
    let v = summary(&d);
    let (t, h) = reports(&v);
    (v, t, h)
}

#[test]
fn the_lpa_section_in_both_reports() {
    let w = Scratch::new("lpa");
    let (v, txt, html) = lpa_sample(&w, "low", 20.0, 30.0, &[(R1, "A/G"), (R2, "T/T")]);
    let sec = at(&v, &["sections", "lpa"]).expect("the section");
    assert!(
        !["missing", "unreadable", ""].contains(&s(at(sec, &["state"]))),
        "{sec:?}"
    );
    assert_eq!(s(at(sec, &["data", "flag"])), "high");
    assert!(
        s(at(sec, &["data", "flag_text"]))
            .contains("1 risk allele at rs10455872 and rs3798220 and a low KIV-2 repeat count"),
        "{sec:?}"
    );
    for want in [
        "## Lp(a): LPA KIV-2 Repeat Count and Risk Variants",
        "KIV-2 repeats (both chromosomes): 20.0",
        &format!("{R1}: A/G (1 risk allele)"),
        &format!("{R2}: T/T (0 risk alleles)"),
        ADVICE,
    ] {
        assert!(txt.contains(want), "the text section: {want}");
    }
    assert!(
        html.contains("<h2>Lp(a): LPA KIV-2 Repeat Count</h2>")
            && html.contains("Genetic sign of high Lp(a)")
    );
    assert!(html.contains(&e(ADVICE)));
    let (v, txt, html) = lpa_sample(&w, "typ", 38.0, 30.0, &[(R1, "A/A"), (R2, "T|T")]);
    assert_eq!(s(at(&v, &["sections", "lpa", "data", "flag"])), ".");
    assert!(
        txt.contains("No genetic sign of a high Lp(a)")
            && !html.contains("Genetic sign of high Lp(a)")
    );
    let (_, txt, _) = lpa_sample(&w, "nov", 38.0, 30.0, &[]);
    assert!(
        txt.contains("rs10455872 and rs3798220: not read (step 42, the traits panel, did not run)"),
        "{txt}"
    );
    let (_, txt, _) = lpa_sample(&w, "thin", 38.0, 3.0, &[(R1, "."), (R2, "C/C")]);
    for want in [
        "KIV-2 repeats: too little depth",
        &format!("{R1}: . (not read)"),
        &format!("{R2}: C/C (2 risk alleles)"),
    ] {
        assert!(txt.contains(want), "{want}");
    }
    let (v, txt, _) = lpa_sample(&w, "thin0", 38.0, 3.0, &[]);
    assert_eq!(s(at(&v, &["sections", "lpa", "data", "flag"])), ".");
    assert!(
        txt.contains("Not determined") && !txt.contains("No genetic sign"),
        "never a negative result"
    );
    let none = w.dir("none/S");
    let (t, h) = reports(&summary(&none));
    assert!(
        !h.contains("<h2>Lp(a)") && !t.contains("## Lp(a)"),
        "no card and no text section without step 46"
    );
}

#[test]
fn kiv2_4925_prints_between_the_count_and_the_risk_variants() {
    let w = Scratch::new("lpa-null");
    let cav = og_render::sections::lpa::LPA_NULL_CAVEAT;
    let (v, txt, _) = lpa_sample(&w, "nul0", 38.0, 30.0, &[]);
    assert_eq!(
        s(at(
            &v,
            &["sections", "lpa", "data", "null_allele", "status"]
        )),
        "not assessed"
    );
    let beside = |t: &str| {
        t.split_once(COUNT38)
            .map_or(String::new(), |x| x.1.split(R1).next().unwrap().to_owned())
    };
    assert!(
        beside(&txt).contains("KIV-2 4925G>A null allele (rs1484919544): not assessed\n"),
        "{txt}"
    );
    assert!(txt.contains("step 46 wrote no S_lpa_kiv2_4925.tsv") && txt.contains(cav));
    let d = w.path("nul0/S/lpa");
    let p = w.put(
        "nul0.pileup",
        &[
            160_646_214,
            160_640_667,
            160_635_123,
            160_629_574,
            160_624_030,
            160_618_484,
            160_612_932,
        ]
        .iter()
        .map(|p| {
            format!(
                "chr6\t{p}\tN\t100\t^]{}{}$\t{}\n",
                "C".repeat(95),
                "t".repeat(5),
                "I".repeat(100)
            )
        })
        .collect::<String>(),
    );
    let (rc, _, err) = og_cn(&[
        "marker-read",
        "LPA_KIV2_4925",
        "--catalog",
        &repo("bin/cn_loci.tsv"),
        "--pileup",
        &p,
        "--copy-table",
        &format!("{d}/S_lpa_kiv2.tsv"),
        "--out",
        &format!("{d}/S_lpa_kiv2_4925.tsv"),
    ]);
    assert_eq!(rc, 0, "{err}");
    let (txt, html) = reports(&summary(&w.path("nul0/S")));
    let b = beside(&txt);
    assert!(b.contains("KIV-2 4925G>A null allele (rs1484919544): present, 35 of 700 reads carry it, about 2.0 copies\n"), "{b}");
    assert!(
        b.contains("cannot say which chromosome or which repeat") && b.contains(cav),
        "{b}"
    );
    let card = html
        .split_once("<h2>Lp(a): LPA KIV-2 Repeat Count</h2>")
        .map_or("", |x| x.1.split("<h2>").next().unwrap());
    assert!(
        card.contains(&e(
            "KIV-2 4925G>A null allele (rs1484919544): present, 35 of 700 reads carry it"
        )),
        "{card}"
    );
    assert!(card.contains(&e(cav)));
    std::fs::write(
        format!("{d}/S_lpa_kiv2_4925.tsv"),
        "metric\tvalue\nmessage\tno status\n",
    )
    .unwrap();
    let (txt, _) = reports(&summary(&w.path("nul0/S")));
    assert!(txt.contains("KIV-2 4925G>A null allele (rs1484919544): not assessed\n"));
    assert!(
        txt.contains("S_lpa_kiv2_4925.tsv has no status row"),
        "never a result"
    );
}

// Step 42's rules in the report

#[test]
fn each_multi_snp_rule_prints_its_trait_and_result_in_the_html_report() {
    let rules = read(&repo("assets/genotype_rules.tsv"));
    let sites = |rid: &str| -> Vec<Vec<String>> {
        rules
            .lines()
            .map(|l| l.split('\t').collect::<Vec<_>>())
            .filter(|c| c.len() >= 7 && c[0] == "site" && c[1] == rid)
            .map(|c| c[3..7].iter().map(|x| (*x).to_owned()).collect())
            .collect()
    };
    for (rid, gts, want) in [
        ("APOL1", ["G/G", "G/G"], "G1/G1"),
        ("SERPINA1", ["T/T", "T/T"], "PiZZ"),
        ("HFE", ["A/A", "C/C"], "C282Y/C282Y"),
        ("GC", ["C/C", "G/G"], "Gc1S/Gc1S"),
        ("CAD9p21", ["G/G", "C/C"], "two risk copies"),
        ("MAPT", ["G/G", "G/G"], "H2/H2"),
    ] {
        let w = Scratch::new(&format!("rules-{rid}"));
        let g: String = sites(rid)
            .iter()
            .zip(gts)
            .map(|(s, gt)| format!("{}\t{gt}\n", s.join("\t")))
            .collect();
        let g = w.put("g.tsv", &g);
        let d = w.dir("S");
        let (out, rout) = (
            w.put("S/traits/S_traits.tsv", ""),
            w.path("S/traits/S_genotype_rules.tsv"),
        );
        let args: Vec<String> = [
            "table",
            "--table",
            &repo("assets/traits.tsv"),
            "--rules",
            &repo("assets/genotype_rules.tsv"),
            "--genotypes",
            &g,
            "--sample",
            "S",
            "--out",
            &out,
            "--rules-out",
            &rout,
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
        assert_eq!(og_steps::traits::cli(&args), 0);
        let text = read(&rout);
        let mut lines = text.lines();
        let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let line: HashMap<&str, &str> = lines
            .map(|l| {
                head.iter()
                    .copied()
                    .zip(l.split('\t'))
                    .collect::<HashMap<_, _>>()
            })
            .find(|r| r["rule"] == rid)
            .unwrap();
        assert_eq!(line["result"], want, "{rid}");
        let (_, html) = reports(&summary(&d));
        assert!(
            html.contains(&e(line["trait"])) && html.contains(&e(want)),
            "{rid}: the HTML report"
        );
    }
}

#[test]
fn a_refused_traits_table_leaves_outputs_the_collector_reads() {
    let w = Scratch::new("plain");
    let g = w.put("g.tsv", "");
    let d = w.dir("S");
    let (out, rout) = (
        w.path("S/traits/S_traits.tsv"),
        w.path("S/traits/S_genotype_rules.tsv"),
    );
    std::fs::create_dir_all(Path::new(&out).parent().unwrap()).unwrap();
    let run = |t: &str| {
        let args: Vec<String> = [
            "table",
            "--table",
            t,
            "--rules",
            &repo("assets/genotype_rules.tsv"),
            "--genotypes",
            &g,
            "--sample",
            "S",
            "--out",
            &out,
            "--rules-out",
            &rout,
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
        og_steps::traits::cli(&args)
    };
    assert_eq!(run(&repo("assets/traits.tsv")), 0);
    let bad = w.put(
        "bad.tsv",
        &read(&repo("assets/traits.tsv")).replacen("\tMCM6\t", "\t\"MCM6\"\t", 1),
    );
    assert_eq!(run(&bad), 1);
    let v = summary(&d);
    let sec = at(&v, &["sections", "traits"]).expect("the traits section");
    assert_eq!(s(at(sec, &["state"])), "ok");
    assert!(matches!(at(sec, &["data", "rows"]), Some(Value::Array(a)) if !a.is_empty()));
}
