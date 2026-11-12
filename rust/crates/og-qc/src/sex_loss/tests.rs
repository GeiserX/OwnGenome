//! The estimate on synthetic mosdepth windows and VCFs: what each case must
//! read, by meaning, and two controls that show the method's corrections
//! matter. The byte-for-byte check against the Python this replaced is the
//! parity harness (`tests/parity`, tool `qc-sex-loss`); the full-precision
//! values below come from CPython 3.15 on the same inputs.

use super::*;
use std::fmt::Write as _;
use std::path::PathBuf;

/// A small deterministic generator (SplitMix64) with a normal draw.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn gauss(&mut self, mu: f64, sigma: f64) -> f64 {
        let u1 = 1.0 - self.uniform();
        let u2 = self.uniform();
        mu + sigma * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo + 1)
    }
}

struct Dir(PathBuf);

impl Dir {
    fn new(name: &str) -> Dir {
        let p = std::env::temp_dir().join(format!("og-sex-loss-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&p);
        std::fs::create_dir_all(&p).unwrap();
        Dir(p)
    }
    fn put(&self, name: &str, body: &str) -> String {
        let p = self.0.join(name);
        std::fs::write(&p, body).unwrap();
        p.to_str().unwrap().to_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// mosdepth `--by 500` windows over `[start, end)` at `depth` (sd 10%).
fn windows(out: &mut String, chrom: &str, start: u64, end: u64, depth: f64, rng: &mut Rng) {
    let mut s = start - start % 500;
    while s < end {
        let d = if depth > 0.0 {
            rng.gauss(depth, depth * 0.1).max(0.0)
        } else {
            0.0
        };
        let _ = writeln!(out, "{chrom}\t{s}\t{}\t{d:.2}", s + 500);
        s += 500;
    }
}

#[derive(Clone, Copy)]
struct Genome {
    male: bool,
    x_depth: Option<f64>,
    y_depth: f64,
    prefix: &'static str,
    slice_only: bool,
    seed: u64,
}

const MAN: Genome = Genome {
    male: true,
    x_depth: None,
    y_depth: 15.0,
    prefix: "chr",
    slice_only: false,
    seed: 1,
};
const WOMAN: Genome = Genome { male: false, ..MAN };

/// A 30x genome's windows: chr1 and chr2 slices with a gap and a
/// satellite, chrX (PAR1, single-copy, XTR) and for a man every kind of
/// chrY piece.
fn regions(g: Genome) -> String {
    let mut rng = Rng(g.seed);
    let mut o = String::new();
    let c = |n: &str| format!("{}{n}", g.prefix);
    windows(&mut o, &c("1"), 1_000_000, 4_000_000, 30.0, &mut rng);
    windows(&mut o, &c("1"), 4_000_000, 4_100_000, 0.0, &mut rng); // a gap
    windows(&mut o, &c("2"), 5_000_000, 6_000_000, 30.0, &mut rng);
    windows(&mut o, &c("2"), 6_000_000, 6_050_000, 400.0, &mut rng); // a satellite
    if g.slice_only {
        windows(&mut o, &c("Y"), 2_790_000, 2_800_000, 15.0, &mut rng);
        return o;
    }
    let xd = g.x_depth.unwrap_or(if g.male { 15.0 } else { 30.0 });
    windows(&mut o, &c("X"), 1_000_000, 2_000_000, 30.0, &mut rng); // PAR1
    windows(&mut o, &c("X"), 73_200_000, 74_500_000, xd, &mut rng);
    let xtr = if g.male { 30.0 } else { 0.0 };
    windows(&mut o, &c("X"), 89_500_000, 90_500_000, xtr, &mut rng); // XTR
    if g.male {
        let yd = g.y_depth;
        windows(&mut o, &c("Y"), 0, 2_781_479, 0.0, &mut rng); // PAR, masked
        windows(&mut o, &c("Y"), 2_781_479, 3_050_044, yd, &mut rng); // single copy
        windows(&mut o, &c("Y"), 3_050_044, 6_235_111, 15.0, &mut rng); // XTR
        windows(&mut o, &c("Y"), 6_748_713, 7_574_481, yd, &mut rng); // single copy
        windows(&mut o, &c("Y"), 7_574_481, 9_046_914, 60.0, &mut rng); // ampliconic
        windows(&mut o, &c("Y"), 10_316_944, 10_544_039, 900.0, &mut rng); // centromere
        windows(&mut o, &c("Y"), 26_673_214, 56_673_214, 0.0, &mut rng); // heterochromatin
    }
    o
}

/// Depth that drifts: each 200 kb stretch of chr1-22 at its own level (sd
/// 8%); chrY only on 2.78-3.00 Mb, at `y_level` times one copy.
fn drift_regions(y_level: f64) -> String {
    let mut rng = Rng(5);
    let mut o = String::new();
    for (chrom, start) in [("chr1", 1_000_000), ("chr2", 5_000_000)] {
        for s in (start..start + 10_000_000).step_by(200_000) {
            let level = 30.0 * rng.gauss(1.0, 0.08);
            windows(&mut o, chrom, s, s + 200_000, level, &mut rng);
        }
    }
    windows(
        &mut o,
        "chrY",
        2_781_479,
        3_000_000,
        15.0 * y_level,
        &mut rng,
    );
    o
}

/// Heterozygous SNVs: chr1 at 0.5, chrX at 0.5 ± d; `spread` adds a
/// per-site spread of the true fraction on every chromosome (a noisy
/// library). Then records the estimate must skip.
fn vcf(d: f64, spread: f64, seed: u64) -> String {
    let mut rng = Rng(seed);
    let mut o = String::from(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
    );
    let mut site = |o: &mut String, chrom: &str, pos: u64, p: f64| {
        let n = rng.range(20, 40);
        let p = if spread > 0.0 {
            (p + rng.gauss(0.0, spread)).clamp(0.05, 0.95)
        } else {
            p
        };
        let a = (0..n).filter(|_| rng.uniform() < p).count() as u64;
        let _ = writeln!(
            o,
            "{chrom}\t{pos}\t.\tA\tG\t40\tPASS\t.\tGT:AD:DP\t0/1:{},{a}:{n}",
            n - a
        );
    };
    for i in 0..3000 {
        site(&mut o, "chr1", 1_000_000 + 1000 * i, 0.5);
    }
    let mut side = Rng(seed + 100);
    for i in 0..600 {
        let p = if side.uniform() < 0.5 {
            0.5 + d
        } else {
            0.5 - d
        };
        site(&mut o, "chrX", 73_200_000 + 2000 * i, p);
    }
    // hom-ALT, not PASS, an indel, PAR1, depth below 10
    o.push_str("chrX\t73100000\t.\tA\tG\t40\tPASS\t.\tGT:AD\t1/1:0,30\n");
    o.push_str("chrX\t73100100\t.\tA\tG\t40\tRefCall\t.\tGT:AD\t0/1:5,25\n");
    o.push_str("chrX\t73100200\t.\tAT\tA\t40\tPASS\t.\tGT:AD\t0/1:5,25\n");
    o.push_str("chrX\t1500000\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:2,28\n");
    o.push_str("chrX\t73100300\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:1,5\n");
    o
}

fn table(inferred: &str, declared: &str) -> Vec<(String, String)> {
    vec![
        ("inferred_sex".to_owned(), inferred.to_owned()),
        ("declared_sex".to_owned(), declared.to_owned()),
    ]
}

struct Row(Vec<(&'static str, String)>);

impl Row {
    fn get(&self, k: &str) -> &str {
        &self.0.iter().find(|(key, _)| *key == k).unwrap().1
    }
    fn num(&self, k: &str) -> f64 {
        self.get(k).parse().unwrap_or(f64::NAN)
    }
    fn near(&self, k: &str, want: f64, tol: f64) -> bool {
        (self.num(k) - want).abs() <= tol
    }
    fn zero_inside(&self) -> bool {
        self.num("mosaic_loss_depth_low") <= 0.0 && 0.0 <= self.num("mosaic_loss_depth_high")
    }
}

impl std::fmt::Debug for Row {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_map()
            .entries(self.0.iter().map(|(k, v)| (k, v)))
            .finish()
    }
}

fn run(t: &[(String, String)], regions: Option<&str>, vcf: Option<&str>) -> Row {
    Row(estimate(t, regions, vcf))
}

#[test]
fn a_man_with_one_y_reads_half_and_no_loss() {
    let d = Dir::new("man");
    let r = run(
        &table("male", "male"),
        Some(&d.put("m.bed", &regions(MAN))),
        None,
    );
    assert!(r.near("mosaic_depth_ratio", 0.5, 0.02), "{r:?}");
    assert!(r.near("mosaic_loss_depth", 0.0, 0.03), "{r:?}");
    assert!(r.zero_inside(), "{r:?}");
    assert_eq!(r.get("mosaic_call"), "not_detected");
    assert!(
        r.get("mosaic_reason")
            .contains("below that nothing can be told at 30x")
    );
    assert_eq!(r.get("mosaic_chrom"), "chrY");
    assert_eq!(r.get("mosaic_floor"), "0.100");
    assert_eq!(r.get("mosaic_depth_ratio_expected"), "0.5");
    assert!(r.num("mosaic_depth_blocks") > 100.0, "{r:?}");
    assert_eq!(r.get("mosaic_loss_af"), "");
    assert_eq!(r.get("mosaic_x_het_sites"), "");
}

#[test]
fn drift_along_the_genome_widens_the_range() {
    let d = Dir::new("drift");
    let path = d.put("drift.bed", &drift_regions(1.13));
    let r = run(&table("male", "male"), Some(&path), None);
    assert!(r.zero_inside(), "{r:?}");
    assert_eq!(r.get("mosaic_call"), "not_detected");
    // Control: from the blocks alone, without the drift, the range would
    // leave zero out and read a gain.
    let b = depth_blocks(&path, Sex::Male).unwrap();
    let ((ma, sea), (ms, ses)) = (median_se(&b.auto), median_se(&b.sex));
    let f = 1.0 - 2.0 * ms / ma;
    let se_f = 2.0 * ms / ma * ((ses / ms).powi(2) + (sea / ma).powi(2)).sqrt();
    assert!(f + 1.96 * se_f < 0.0, "{f} {se_f}");
}

#[test]
fn a_lost_y_is_detected_and_an_extra_y_is_a_gain() {
    let d = Dir::new("loy");
    let g = Genome {
        y_depth: 10.5,
        ..MAN
    };
    let r = run(
        &table("male", "male"),
        Some(&d.put("m30.bed", &regions(g))),
        None,
    );
    assert!(r.near("mosaic_loss_depth", 0.30, 0.03), "{r:?}");
    assert_eq!(r.get("mosaic_call"), "detected");
    assert!(r.num("mosaic_loss_depth_low") > 0.10, "{r:?}");
    let g = Genome {
        y_depth: 30.0,
        ..MAN
    };
    let r = run(
        &table("male", "male"),
        Some(&d.put("xyy.bed", &regions(g))),
        None,
    );
    assert_eq!(r.get("mosaic_call"), "gain");
    assert!(r.near("mosaic_loss_depth", -1.0, 0.05), "{r:?}");
    assert!(r.get("mosaic_reason").contains("47,XYY"));
}

#[test]
fn a_woman_reads_from_depth_and_from_allele_fractions() {
    let d = Dir::new("woman");
    let reg = d.put("f.bed", &regions(WOMAN));
    let r = run(
        &table("female", "female"),
        Some(&reg),
        Some(&d.put("f.vcf", &vcf(0.0, 0.0, 2))),
    );
    assert!(r.near("mosaic_depth_ratio", 1.0, 0.03), "{r:?}");
    assert!(r.near("mosaic_loss_depth", 0.0, 0.05), "{r:?}");
    assert!(r.num("mosaic_loss_af") < 0.10, "{r:?}");
    assert_eq!(r.num("mosaic_loss_af_low"), 0.0);
    assert_eq!(r.get("mosaic_x_het_sites"), "600");
    assert_eq!(r.get("mosaic_call"), "not_detected");

    let d30 = 0.3 / (2.0 * (2.0 - 0.3));
    let g = Genome {
        x_depth: Some(30.0 * 0.85),
        ..WOMAN
    };
    let r = run(
        &table("female", "female"),
        Some(&d.put("f30.bed", &regions(g))),
        Some(&d.put("f30.vcf", &vcf(d30, 0.0, 3))),
    );
    assert!(r.near("mosaic_loss_depth", 0.30, 0.04), "{r:?}");
    assert!(r.near("mosaic_loss_af", 0.30, 0.08), "{r:?}");
    assert_eq!(r.get("mosaic_call"), "detected");
    assert!(r.get("mosaic_reason").contains("lost one X"));

    let r = run(
        &table("female", "female"),
        Some(&reg),
        Some(&d.put("fcn.vcf", &vcf(d30, 0.0, 4))),
    );
    assert_eq!(r.get("mosaic_call"), "detected");
    assert!(r.get("mosaic_reason").contains("copy-neutral"));
    assert!(r.near("mosaic_loss_depth", 0.0, 0.05), "{r:?}");
}

#[test]
fn a_noisy_library_is_not_read_as_a_shift() {
    let d = Dir::new("noisy");
    let path = d.put("noisy.vcf", &vcf(0.0, 0.06, 6));
    let r = run(
        &table("female", "female"),
        Some(&d.put("f.bed", &regions(WOMAN))),
        Some(&path),
    );
    // 600 noisy sites leave a wide range: the point estimate may sit above
    // the floor, its low end must not.
    assert_eq!(r.get("mosaic_call"), "not_detected");
    assert!(r.num("mosaic_loss_af_low") < 0.10, "{r:?}");
    // Control: against binomial spread alone the same sites would read a loss.
    let h = het_sites(&path).unwrap();
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let xn: Vec<f64> = h.x_n.iter().map(|&n| n as f64).collect();
    let naive = (mean(&h.x_u) - 1.0) / (4.0 * (mean(&xn) - 1.0));
    assert!(naive > 0.0 && af_loss(naive) > 0.10, "{naive}");
}

#[test]
fn nothing_to_read_is_not_estimated_with_the_reason() {
    let d = Dir::new("none");
    let reg = d.put("f.bed", &regions(WOMAN));
    let r = run(&table("unknown", ""), Some(&reg), None);
    assert_eq!(r.get("mosaic_call"), "not_estimated");
    assert!(r.get("mosaic_reason").contains("no sex"));
    let r = run(&table("male", "male"), None, None);
    assert_eq!(r.get("mosaic_call"), "not_estimated");
    assert!(r.get("mosaic_reason").contains("step 16b"));
    let g = Genome {
        slice_only: true,
        ..MAN
    };
    let r = run(
        &table("male", "male"),
        Some(&d.put("slice.bed", &regions(g))),
        None,
    );
    assert_eq!(r.get("mosaic_call"), "not_estimated");
    assert!(r.get("mosaic_reason").contains("too few 10 kb blocks"));
    let r = run(
        &table("female", "female"),
        None,
        Some(&d.put("p.vcf.gz", "placeholder\n")),
    );
    assert_eq!(r.get("mosaic_call"), "not_estimated");
    assert!(
        r.get("mosaic_reason")
            .contains("allele fractions not estimated")
    );
    let r = run(&table("female", "female"), Some(&reg), None);
    assert_eq!(r.get("mosaic_call"), "not_detected");
    assert!(r.get("mosaic_reason").contains("no VCF"));
}

#[test]
fn the_sex_from_the_reads_wins_and_bare_contig_names_work() {
    let d = Dir::new("sex");
    let reg = d.put("m.bed", &regions(MAN));
    let r = run(&table("male", "female"), Some(&reg), None);
    assert_eq!(
        (r.get("mosaic_sex"), r.get("mosaic_sex_from")),
        ("male", "reads")
    );
    let r = run(&table("unknown", "male"), Some(&reg), None);
    assert_eq!(
        (r.get("mosaic_sex"), r.get("mosaic_sex_from")),
        ("male", "declared")
    );
    let g = Genome {
        y_depth: 10.5,
        prefix: "",
        ..MAN
    };
    let r = run(
        &table("male", "male"),
        Some(&d.put("nochr.bed", &regions(g))),
        None,
    );
    assert!(r.near("mosaic_loss_depth", 0.30, 0.03), "{r:?}");
}

#[test]
fn unreadable_inputs_name_the_python_error() {
    let d = Dir::new("errors");
    let t = table("male", "male");
    let r = run(&t, Some("/nonexistent/S.regions.bed.gz"), None);
    assert_eq!(
        r.get("mosaic_reason"),
        "depth not estimated: could not read S.regions.bed.gz: FileNotFoundError: \
         [Errno 2] No such file or directory: '/nonexistent/S.regions.bed.gz'"
    );
    let r = run(
        &t,
        Some(&d.put("short.bed", "chr1\t0\t500\t12\nchr1\t500\n")),
        None,
    );
    assert!(r.get("mosaic_reason").ends_with(
        "could not read short.bed: ValueError: not enough values to unpack (expected 4, got 2)"
    ));
    let r = run(&t, Some(&d.put("int.bed", "chrY\tx\t2790500\t15\n")), None);
    assert!(
        r.get("mosaic_reason")
            .ends_with("ValueError: invalid literal for int() with base 10: 'x'")
    );
    let r = run(&t, Some(&d.put("float.bed", "chr1\t0\t500\tdeep\n")), None);
    assert!(
        r.get("mosaic_reason")
            .ends_with("ValueError: could not convert string to float: 'deep\\n'")
    );
    let r = run(
        &t,
        Some(&d.put("nan.bed", "chr1\t0\t500\t12\nchr1\t500\t1000\tnan\n")),
        None,
    );
    assert!(
        r.get("mosaic_reason")
            .ends_with("ValueError: line 2: a depth of nan, not a finite number")
    );
}

#[test]
fn ratio_is_correctly_rounded() {
    // num / den in CPython, as float.hex(), on random big ints.
    let cases: [(u128, u128, f64); 10] = [
        ((1 << 60) + 1, 3, 3.843_071_682_022_823e17),
        (100_000_000_000_000_000_000, 7, 1.428_571_428_571_428_7e19),
        (463951277996445748708, 2795742298, 0.0),
        (576413965836357413709867874642195916, 1796035749, 0.0),
        (699997782514302555833621942, 531725357, 0.0),
        (723381256811754296, 2503952635, 0.0),
        (14794762486492979, 2301595701, 0.0),
        (10548175869641935, 1703729694, 0.0),
        // A tie: 2**53 + 1 rounds to even.
        ((1 << 54) + 2, 2, 9_007_199_254_740_992.0),
        (25, 15, 25.0 / 15.0),
    ];
    let hex = [
        "",
        "",
        "0x1.351ac234de39fp+37",
        "0x1.097902c916c53p+88",
        "0x1.24505dff67902p+60",
        "0x1.13832fef4688bp+28",
        "0x1.88563709fe4eap+22",
        "0x1.79e1e8af01e66p+22",
        "",
        "",
    ];
    for (i, &(num, den, want)) in cases.iter().enumerate() {
        let want = if hex[i].is_empty() {
            want
        } else {
            from_hex(hex[i])
        };
        assert_eq!(ratio(num, den), want, "{num} / {den}");
    }
}

/// A positive normal float from Python's `float.hex()`.
fn from_hex(s: &str) -> f64 {
    let (mant, exp) = s.trim_start_matches("0x1.").split_once("p+").unwrap();
    let frac = u64::from_str_radix(mant, 16).unwrap();
    let exp: u64 = exp.parse().unwrap();
    f64::from_bits(((exp + 1023) << 52) | frac)
}

#[test]
fn a_read_count_past_32_bits_skips_the_site() {
    assert_eq!(het("GT:AD", "0/1:0,4294967296"), None);
    assert_eq!(
        het("GT:AD", "0/1:1,4294967294"),
        Some((4294967295, ratio(4294967293u128.pow(2), 4294967295)))
    );
    assert_eq!(het("GT:AD", "0|1:5,15"), Some((20, 5.0)));
    assert_eq!(het("GT:AD", "1/1:5,15"), None);
    assert_eq!(het("GT:AD", "0/1:5"), None);
    assert_eq!(het("GT", "0/1"), None);
}

#[test]
fn autosomes_are_one_to_twenty_two_as_python_writes_them() {
    for n in 1..=22 {
        assert!(is_autosome(&n.to_string()));
    }
    for s in ["0", "23", "01", "+1", "X", "", "1 ", "220"] {
        assert!(!is_autosome(s), "{s}");
    }
}

#[test]
fn single_copy_needs_the_whole_window_in_one_piece() {
    assert!(single_copy(SINGLE_COPY_Y, 2_781_479, 2_781_979));
    assert!(!single_copy(SINGLE_COPY_Y, 2_781_478, 2_781_978));
    assert!(!single_copy(SINGLE_COPY_Y, 3_049_800, 3_050_300));
    assert!(!single_copy(SINGLE_COPY_Y, -5, 10));
}

/// CPython 3.15 (`PYTHON_IMAGE`) on the same two files, through the Python
/// this replaced: `af_estimate` and `depth_estimate`, to the last bit. The
/// VCF's seed is one where a naive two-pass float variance moves the last
/// digit of both ends of the range (0.2700004351269259 for 0.27000043512692595,
/// 0.32504770224985213 for 0.3250477022498521), so it fails here.
#[test]
fn the_estimates_match_cpython_to_the_last_bit() {
    let d = Dir::new("bits");
    let d30 = 0.3 / (2.0 * (2.0 - 0.3));
    let af = af_estimate(&d.put("af.vcf", &vcf(d30, 0.0, 6)))
        .unwrap()
        .unwrap();
    assert_eq!(
        af,
        Af {
            sites: 600,
            loss: 0.2994563206373466,
            low: 0.27000043512692595,
            high: 0.3250477022498521,
        }
    );
    let dep = depth_estimate(&d.put("depth.bed", &drift_regions(1.13)), Sex::Male)
        .unwrap()
        .unwrap();
    assert_eq!(
        dep,
        Depth {
            ratio: 0.556056015068381,
            blocks: 22,
            loss: -0.11211203013676196,
            low: -0.24953092365800184,
            high: 0.025306863384477907,
        }
    );
}

#[test]
fn a_table_that_is_not_utf8_names_its_own_line() {
    let d = Dir::new("table-utf8");
    let p = d.0.join("S_sample_qc.tsv");
    std::fs::write(&p, b"key\tvalue\n\n\nsex\t\xff\n").unwrap();
    assert_eq!(
        read_table(p.to_str().unwrap()).unwrap_err(),
        "line 4 is not UTF-8"
    );
}
