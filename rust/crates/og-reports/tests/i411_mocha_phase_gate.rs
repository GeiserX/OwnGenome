//! Issue #411 (trait-65): `og mocha report` marks a sample "not assessed"
//! when an autosome has too few of its heterozygous sites phased for MoChA
//! to follow a haplotype, not only when no site at all is phased.
//!
//! Before the fix one phased site in the whole genome was enough for
//! "assessed", and a real 6% CN-LOH above the floor then went unreported
//! under a summary that read as a negative result. The minimum is a share
//! per autosome; `docs/51-mocha.md` gives the shares measured on two public
//! genomes and how the minimum follows from them.

use std::path::{Path, PathBuf};

/// The minimum, in percent: `PHASED_MIN_PERCENT` in `src/mocha.rs`.
const MIN_PERCENT: u64 = 80;

const STATS: &str =
    "sample_id\tcall_rate\tcov_median\tbaf_auto\tn_hets\nS\t1.0000\t31\t0.0124\t2090000\n";
const CALLS: &str = "sample_id\tchrom\tbeg_GRCh38\tend_GRCh38\tlength\tp_arm\tq_arm\tn_hets\tn50_hets\tbdev\tbdev_se\trel_cov\tlod_baf_phase\tlod_lrr_baf\tlod_baf_conc\ttype\tcf\n\
S\tchr20\t30100000\t42300000\t12200000\tN\tY\t4210\t350000\t0.031\t0.004\t1.86\t48.2\t35.1\tnan\tLoss\t0.120\n";

struct Run {
    dir: PathBuf,
}

impl Run {
    fn new(name: &str) -> Run {
        let dir = std::env::temp_dir().join(format!("og-i411-mocha-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Run { dir }
    }

    fn path(&self, f: &str) -> String {
        self.dir.join(f).to_string_lossy().into_owned()
    }

    /// Runs `og mocha report` on these sites rows (chrom, hets, phased);
    /// the summary as (metric, value) pairs and the one call row.
    fn report(&self, sites: &[(&str, u64, u64)]) -> (Vec<(String, String)>, Vec<String>) {
        let mut text = String::new();
        for (c, hets, phased) in sites {
            text.push_str(&format!(
                "{c}\t10001\t100010001\t{}\t{hets}\t{phased}\n",
                hets * 3
            ));
        }
        std::fs::write(self.path("sites.tsv"), text).unwrap();
        std::fs::write(self.path("stats.tsv"), STATS).unwrap();
        std::fs::write(self.path("calls.tsv"), CALLS).unwrap();
        let args: Vec<String> = [
            "--sample",
            "S",
            "--calls",
            &self.path("calls.tsv"),
            "--stats",
            &self.path("stats.tsv"),
            "--sites",
            &self.path("sites.tsv"),
            "--phasing",
            "beagle",
            "--out-calls",
            &self.path("o.tsv"),
            "--out-summary",
            &self.path("s.tsv"),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let code = og_reports::mocha::cli(Some("report"), &args, &mut out, &mut err);
        assert_eq!(code, 0, "{}", String::from_utf8_lossy(&err));
        let summary = read(&self.dir.join("s.tsv"))
            .lines()
            .skip(1)
            .map(|l| {
                let (k, v) = l.split_once('\t').unwrap();
                (k.to_owned(), v.to_owned())
            })
            .collect();
        let call = read(&self.dir.join("o.tsv"))
            .lines()
            .nth(1)
            .unwrap()
            .split('\t')
            .map(str::to_owned)
            .collect();
        (summary, call)
    }
}

impl Drop for Run {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn read(p: &Path) -> String {
    std::fs::read_to_string(p).unwrap()
}

fn get<'a>(summary: &'a [(String, String)], k: &str) -> &'a str {
    summary
        .iter()
        .find(|(n, _)| n == k)
        .map(|(_, v)| v.as_str())
        .unwrap()
}

fn autosomes(hets: u64, phased: impl Fn(usize) -> u64) -> Vec<(String, u64, u64)> {
    (1..=22)
        .map(|i| (format!("chr{i}"), hets, phased(i)))
        .collect()
}

fn rows(v: &[(String, u64, u64)]) -> Vec<(&str, u64, u64)> {
    v.iter().map(|(c, h, p)| (c.as_str(), *h, *p)).collect()
}

/// Red first: the issue's case, 22 autosomes of 95,000 heterozygous sites
/// with one site phased in the whole genome. Before the fix: "assessed".
#[test]
fn red_first_one_phased_site_in_the_genome_is_not_assessed() {
    let run = Run::new("one-site");
    let sites = autosomes(95_000, |i| u64::from(i == 1));
    let (summary, call) = run.report(&rows(&sites));
    assert_eq!(get(&summary, "status"), "not assessed", "{summary:?}");
    let msg = get(&summary, "message");
    assert!(msg.starts_with("under 0.80 of the heterozygous sites are phased on chr1 (1 of 95000), chr2 (0 of 95000)"), "{msg}");
    assert!(
        msg.ends_with("chr22 (0 of 95000), so MoChA cannot follow a haplotype there"),
        "{msg}"
    );
    assert_eq!(get(&summary, "phased_het_sites"), "1");
    assert_eq!(call[12], "no");
    assert_eq!(call[13], "sample not assessed");
}

/// Red first: one autosome left unphased (a panel cut to regions) among 21
/// phased ones. The genome-wide share is still 0.95, so only a check per
/// autosome sees it. Before the fix: "assessed".
#[test]
fn red_first_one_unphased_autosome_is_not_assessed() {
    let run = Run::new("one-chrom");
    let sites = autosomes(95_000, |i| if i == 7 { 12 } else { 93_000 });
    let (summary, call) = run.report(&rows(&sites));
    assert_eq!(get(&summary, "status"), "not assessed", "{summary:?}");
    assert_eq!(
        get(&summary, "message"),
        "under 0.80 of the heterozygous sites are phased on chr7 (12 of 95000), so MoChA cannot follow a haplotype there"
    );
    assert_eq!(call[13], "sample not assessed");
}

/// Red first at the minimum minus one site: 95,000 heterozygous sites on
/// chr5 need 76,000 phased; 75,999 is not assessed, 76,000 is.
#[test]
fn red_first_the_minimum_minus_one_site_is_not_assessed() {
    let need = 95_000 * MIN_PERCENT / 100;
    assert_eq!(need, 76_000);
    for (phased5, status) in [(need - 1, "not assessed"), (need, "assessed")] {
        let run = Run::new(&format!("min-{phased5}"));
        let sites = autosomes(95_000, |i| if i == 5 { phased5 } else { 93_000 });
        let (summary, call) = run.report(&rows(&sites));
        assert_eq!(
            get(&summary, "status"),
            status,
            "chr5 {phased5}: {summary:?}"
        );
        let reported = if status == "assessed" { "yes" } else { "no" };
        assert_eq!(call[12], reported, "chr5 {phased5}: {call:?}");
    }
}

/// Control: a whole genome with 0.92 of every autosome phased, under the
/// lowest share measured on one autosome of the two public genomes (0.933,
/// docs/51-mocha.md), is assessed and its call reported; chrX, which the
/// report never reads, is not gated even when none of it is phased.
#[test]
fn the_measured_shares_are_assessed() {
    let run = Run::new("measured");
    let mut sites = autosomes(100_000, |_| 92_000);
    sites.push(("chrX".into(), 40_000, 0));
    let (summary, call) = run.report(&rows(&sites));
    assert_eq!(get(&summary, "status"), "assessed", "{summary:?}");
    assert_eq!(get(&summary, "message"), ".");
    assert_eq!(call[12], "yes", "{call:?}");
}
