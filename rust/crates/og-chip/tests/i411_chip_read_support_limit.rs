//! Issue #411, trait-59: the CHIP read-support limit comes from each coding
//! base's depth, not from a gene's mean depth.
//!
//! Before the fix, `og chip report` took one mean over a gene's coding bases
//! and gave the binomial limit at that mean. A gene with a fifth of its
//! bases unread then read `covered` with a limit near the one of a gene read
//! evenly, and the summary line stated the limit at the median depth.

use og_chip::binom::{detection_limit, p_at_least};
use pycompat::fmt::fixed;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

const GENES: &str = "gene\ttruncating\tmissense\tnote\tpmids\n\
DNMT3A\tall\tR882\tn\t25426837\n\
TET2\tall\tnone\tn\t25426837\n";

/// One `og chip report` run: its stdout and its depth table by gene.
struct Run {
    stdout: String,
    genes: HashMap<String, HashMap<String, String>>,
}

impl Run {
    fn cell(&self, gene: &str, col: &str) -> &str {
        self.genes
            .get(gene)
            .and_then(|r| r.get(col))
            .map_or("<missing>", String::as_str)
    }
}

/// Runs `og chip report` on a BED of `(contig, start, end, gene)` rows and
/// a depth file of `(contig, pos, depth)` lines, in a folder of its own.
fn report(name: &str, bed: &[(&str, u64, u64, &str)], depth: &[(&str, u64, u64)]) -> Run {
    let dir: PathBuf =
        std::env::temp_dir().join(format!("og-chip-i411-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let put = |file: &str, body: String| {
        let p = dir.join(file);
        fs::write(&p, body).unwrap();
        p.to_str().unwrap().to_owned()
    };
    let genes = put("genes.tsv", GENES.to_owned());
    let gff = put("chip.gff3", "##gff-version 3\n".to_owned());
    let calls = put("calls.tsv", String::new());
    let bed = put(
        "chip.bed",
        bed.iter()
            .map(|(c, s, e, g)| format!("{c}\t{s}\t{e}\t{g}\n"))
            .collect(),
    );
    let depth = put(
        "depth.tsv",
        depth
            .iter()
            .map(|(c, p, d)| format!("{c}\t{p}\t{d}\n"))
            .collect(),
    );
    let out = dir.join("S_chip.tsv").to_str().unwrap().to_owned();
    let dout = dir.join("S_chip_depth.tsv").to_str().unwrap().to_owned();
    let args: Vec<String> = [
        "--genes",
        &genes,
        "--gff",
        &gff,
        "--bed",
        &bed,
        "--calls",
        &calls,
        "--depth",
        &depth,
        "--sample",
        "S",
        "--out",
        &out,
        "--depth-out",
        &dout,
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let (mut so, mut se) = (Vec::new(), Vec::new());
    let code = og_chip::cli("report", &args, &mut so, &mut se);
    assert_eq!(code, 0, "{}", String::from_utf8_lossy(&se));
    let text = fs::read_to_string(&dout).unwrap();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let genes = lines
        .map(|ln| {
            let row: HashMap<String, String> = head
                .iter()
                .zip(ln.split('\t'))
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect();
            (row["gene"].clone(), row)
        })
        .collect();
    let _ = fs::remove_dir_all(&dir);
    Run {
        stdout: String::from_utf8(so).unwrap(),
        genes,
    }
}

/// The test's own oracle: the smallest VAF at which the mean over bases of
/// `P(at least 3 reads with the change)` reaches 0.8, by bisection.
fn per_base_limit(depths: &[(u64, u64)]) -> f64 {
    let all: u64 = depths.iter().map(|d| d.1).sum();
    let power = |f: f64| {
        depths
            .iter()
            .map(|&(d, n)| n as f64 * p_at_least(3, d, f))
            .sum::<f64>()
            / all as f64
    };
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if power(mid) >= 0.8 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

/// Every base of `start..end` (0-based, half-open) at `depth`.
fn even(contig: &'static str, start: u64, end: u64, depth: u64) -> Vec<(&'static str, u64, u64)> {
    (start + 1..=end).map(|p| (contig, p, depth)).collect()
}

/// DNMT3A as ten 50-base coding exons on chr2, eight at 33x, exon 3 with no
/// pileup line and exon 7 with lines at depth 0 (reads under mapping
/// quality 20 look so): a fifth of its 500 bases unread. TET2 is 500 bases
/// on chr4 at 33x everywhere.
fn dark_fifth() -> Run {
    let exons: Vec<(u64, u64)> = (0..10).map(|i| (1000 + i * 200, 1050 + i * 200)).collect();
    let mut bed: Vec<(&str, u64, u64, &str)> = exons
        .iter()
        .map(|&(s, e)| ("chr2", s, e, "DNMT3A"))
        .collect();
    bed.push(("chr4", 5000, 5500, "TET2"));
    let mut depth = Vec::new();
    for (i, &(s, e)) in exons.iter().enumerate() {
        match i {
            2 => {}
            6 => depth.extend(even("chr2", s, e, 0)),
            _ => depth.extend(even("chr2", s, e, 33)),
        }
    }
    depth.extend(even("chr4", 5000, 5500, 33));
    report("dark", &bed, &depth)
}

#[test]
fn red_first_a_gene_with_a_fifth_of_its_bases_unread_has_no_limit() {
    let r = dark_fifth();
    // Before: "26.4  0.157  covered", next to 0.152 for 27x everywhere.
    assert_eq!(r.cell("DNMT3A", "mean_depth"), "26.4");
    assert_eq!(r.cell("DNMT3A", "status"), "partly_covered");
    assert_eq!(r.cell("DNMT3A", "detection_limit_vaf"), "");
    assert_eq!(r.cell("DNMT3A", "callable_fraction"), "0.800");
    // The stated limit is TET2's, at 33x: DNMT3A shows a clone of that VAF
    // at 80% of its coding bases.
    assert_eq!(r.cell("DNMT3A", "assessable_fraction"), "0.800");
    let lim33 = detection_limit(33.0).unwrap().unwrap();
    assert_eq!(r.cell("TET2", "status"), "covered");
    assert_eq!(r.cell("TET2", "detection_limit_vaf"), fixed(lim33, 3));
    assert_eq!(r.cell("TET2", "callable_fraction"), "1.000");
    assert_eq!(r.cell("TET2", "assessable_fraction"), "1.000");
    assert!(
        r.stdout.contains("Too few reads for any limit in 1 gene;")
            && r.stdout.contains(": DNMT3A 80%."),
        "{}",
        r.stdout
    );
}

#[test]
fn red_first_a_30x_run_with_uneven_depth_states_the_higher_vaf() {
    // Half of DNMT3A's bases at 20x and half at 40x: a mean of 30x, where
    // one mean gave 13.7%. Base by base the limit is higher.
    let mut depth = even("chr2", 1000, 1250, 20);
    depth.extend(even("chr2", 1250, 1500, 40));
    let r = report("uneven", &[("chr2", 1000, 1500, "DNMT3A")], &depth);
    let want = per_base_limit(&[(20, 250), (40, 250)]);
    let at_mean = detection_limit(30.0).unwrap().unwrap();
    assert!(want > at_mean + 0.01, "{want} vs {at_mean}");
    assert_eq!(r.cell("DNMT3A", "mean_depth"), "30.0");
    assert_eq!(r.cell("DNMT3A", "status"), "covered");
    assert_eq!(r.cell("DNMT3A", "detection_limit_vaf"), fixed(want, 3));
    assert!(
        r.stdout.contains(&format!(
            "a clone with a VAF of about {}% or more",
            fixed(want * 100.0, 0)
        )),
        "{}",
        r.stdout
    );
}

#[test]
fn uniform_30x_states_the_binomial_limit_and_says_what_it_is() {
    let r = report(
        "even30",
        &[("chr2", 1000, 1500, "DNMT3A")],
        &even("chr2", 1000, 1500, 30),
    );
    assert_eq!(r.cell("DNMT3A", "status"), "covered");
    assert_eq!(r.cell("DNMT3A", "detection_limit_vaf"), "0.137");
    assert_eq!(r.cell("DNMT3A", "assessable_fraction"), "1.000");
    assert_eq!(r.cell("TET2", "status"), "not_on_reference");
    assert_eq!(r.cell("TET2", "callable_fraction"), "");
    for want in [
        "Read-support limit (median 30x",
        "a clone with a VAF of about 14% or more shows at least 3 reads with the change 8 times in 10",
        "Mutect2 and its filters call fewer than that.",
        "does not rule one out. Not a diagnosis.",
    ] {
        assert!(r.stdout.contains(want), "{want:?} not in {}", r.stdout);
    }
    assert!(!r.stdout.contains("Too few reads"), "{}", r.stdout);
}

#[test]
fn shallow_everywhere_and_shallow_in_part_are_told_apart() {
    // DNMT3A at 2x everywhere: no base can show 3 reads. TET2 at 60x on
    // half its bases: some can, too few for 8 in 10.
    let mut depth = even("chr2", 1000, 1500, 2);
    depth.extend(even("chr4", 5000, 5250, 60));
    let r = report(
        "shallow",
        &[("chr2", 1000, 1500, "DNMT3A"), ("chr4", 5000, 5500, "TET2")],
        &depth,
    );
    assert_eq!(r.cell("DNMT3A", "status"), "too_shallow");
    assert_eq!(r.cell("DNMT3A", "callable_fraction"), "0.000");
    assert_eq!(r.cell("TET2", "status"), "partly_covered");
    assert_eq!(r.cell("TET2", "callable_fraction"), "0.500");
    // No covered gene: no stated limit, no share at it.
    assert_eq!(r.cell("TET2", "assessable_fraction"), "");
    assert!(r.stdout.contains("no result either way"), "{}", r.stdout);
}

#[test]
fn red_first_a_base_two_genes_share_is_read_for_both() {
    // DNMT3A 1000..1500 and TET2 1400..1600 on chr2, all at 33x. The last
    // 100 bases of DNMT3A are TET2's too. Before, they counted only for
    // TET2, the last gene of the BED, so DNMT3A read them as depth 0:
    // "partly_covered", callable 0.800. Its mean stays the Python's.
    let r = report(
        "shared",
        &[("chr2", 1000, 1500, "DNMT3A"), ("chr2", 1400, 1600, "TET2")],
        &even("chr2", 1000, 1600, 33),
    );
    let lim33 = fixed(detection_limit(33.0).unwrap().unwrap(), 3);
    assert_eq!(r.cell("DNMT3A", "mean_depth"), "26.4");
    assert_eq!(r.cell("DNMT3A", "status"), "covered");
    assert_eq!(r.cell("DNMT3A", "detection_limit_vaf"), lim33);
    assert_eq!(r.cell("DNMT3A", "callable_fraction"), "1.000");
    assert_eq!(r.cell("DNMT3A", "assessable_fraction"), "1.000");
    assert_eq!(r.cell("TET2", "mean_depth"), "33.0");
    assert_eq!(r.cell("TET2", "detection_limit_vaf"), lim33);
    assert_eq!(r.cell("TET2", "callable_fraction"), "1.000");
}
