//! `og cn call --locus LPA_KIV2`: step 46's LPA KIV-2 repeat count from
//! short-read depth, and KIV-2 4925G>A beside it (`og cn marker-read`, `og cn
//! call --marker`). These are the tool checks of the old
//! `tests/test_lpa_kiv2.py`, by meaning; its report checks are
//! `rust/og/tests/genome_steps.rs` and its traits-table check is
//! og-steps' `tests/genome_traits.rs`.
//!
//! On `samtools bedcov` tables made up for each case (summed depth = mean
//! depth x exon length), 30x over the single-copy LPA exons, with 90% of the
//! KIV-2 depth from reads of MAPQ 0 unless a case says otherwise:
//!  1. the arithmetic: a block depth of 30 x (N + 2) / 14 gives N, for N =
//!     12, 20 and 40;
//!  2. the interval: equal depths give a zero-width interval; planted spread
//!     gives sqrt(cv_kiv2^2 / 14 + cv_flank^2 / 15) on N + 2, computed here
//!     by hand; twice the spread gives a wider interval;
//!  3. the class: low, typical and high by the whole interval against 30
//!     and 45;
//!  4. too little depth, no MAPQ 0 read over KIV-2, no chr6 and a depth under
//!     the floor of any count are not counts;
//!  5. the command line: 29 regions, the table, and the refused inputs;
//!  8. KIV-2 4925G>A at its seven GRCh38 copies: present, not seen,
//!     uncertain, too little depth, not assessed; other positions and indel
//!     bases are not counted; step 46's table is left byte for byte;
//!     `call --marker` writes what `marker-read` writes;
//!  9. the demo sample's step 46 table (`og dev demo-sample`'s `LPA_TABLE`) is
//!     the one `og cn call` writes for its depths.

use std::collections::HashMap;
use std::path::PathBuf;

use og_cn::catalog::{Catalog, Locus};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

/// The seven copies of KIV-2 4925G>A on GRCh38.
const COPIES: [i64; 7] = [
    160_646_214,
    160_640_667,
    160_635_123,
    160_629_574,
    160_624_030,
    160_618_484,
    160_612_932,
];

const ADVICE: &str = "A research estimate from read depth, not a diagnostic test. Lp(a) is measured in blood; a 2022 \
                      European Atherosclerosis Society consensus advises measuring it at least once in every adult.";

fn catalog_path() -> String {
    format!("{REPO}/bin/cn_loci.tsv")
}

fn lpa_row() -> Locus {
    Catalog::load(&catalog_path()).unwrap().loci["LPA_KIV2"].clone()
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_default()
}

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!("og-genome-lpa-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn path(&self, rel: &str) -> String {
        let p = self.0.join(rel);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        p.to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
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

type Depth = HashMap<String, f64>;

/// Mean depth per region: the flanks at `flank`, the KIV-2 block at the
/// depth `n` copies give; a spread alternates +x and -x around each mean.
fn depths(r: &Locus, n: f64, flank: f64, ks: f64, fs: f64) -> Depth {
    let k = flank * (n + 2.0) / 14.0;
    let mut m = Depth::new();
    for (i, s) in r.targets.iter().enumerate() {
        m.insert(
            s.name.clone(),
            k * (1.0 + if i % 2 == 0 { ks } else { -ks }),
        );
    }
    for (i, s) in r.controls.iter().enumerate() {
        m.insert(
            s.name.clone(),
            flank * (1.0 + if i % 2 == 0 { fs } else { -fs }),
        );
    }
    m
}

/// The MAPQ >= 1 depths: `share` of the KIV-2 depth is from MAPQ 0.
fn q1(m: &Depth, share: f64) -> Depth {
    m.iter()
        .map(|(n, v)| {
            (
                n.clone(),
                if n.starts_with("kiv2_") {
                    v * (1.0 - share)
                } else {
                    *v
                },
            )
        })
        .collect()
}

fn bedcov(path: &str, r: &Locus, m: &Depth, exact: bool) {
    let text: String = r
        .regions()
        .iter()
        .map(|(c, s, e, n)| {
            let v = m[*n] * (e - s) as f64;
            if exact {
                format!("{c}\t{s}\t{e}\t{n}\t{v:?}\n")
            } else {
                format!("{c}\t{s}\t{e}\t{n}\t{}\n", v.round_ties_even() as i64)
            }
        })
        .collect();
    std::fs::write(path, text).unwrap();
}

fn table(path: &str) -> HashMap<String, String> {
    read(path)
        .lines()
        .skip(1)
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_owned(), v.to_owned()))
        .collect()
}

fn g<'a>(t: &'a HashMap<String, String>, k: &str) -> &'a str {
    t.get(k).map_or("", String::as_str)
}

fn estimate(a: &str, b: &str, out: &str) -> (u8, String) {
    let (rc, _, err) = og_cn(&[
        "call",
        "--catalog",
        &catalog_path(),
        "--locus",
        "LPA_KIV2",
        "--bedcov",
        a,
        "--bedcov-mapq1",
        b,
        "--table",
        out,
    ]);
    (rc, err)
}

/// Step 46's table for mean depths m (every read) and q (MAPQ >= 1); empty
/// tables when `m` is None.
fn assess(name: &str, mq: Option<(&Depth, &Depth)>) -> HashMap<String, String> {
    let w = Scratch::new(name);
    let r = lpa_row();
    let (a, b, o) = (w.path("a.tsv"), w.path("b.tsv"), w.path("o.tsv"));
    match mq {
        Some((m, q)) => {
            bedcov(&a, &r, m, true);
            bedcov(&b, &r, q, true);
        }
        None => {
            std::fs::write(&a, "").unwrap();
            std::fs::write(&b, "").unwrap();
        }
    }
    let (rc, err) = estimate(&a, &b, &o);
    assert_eq!(rc, 0, "{err}");
    table(&o)
}

fn plain(name: &str, n: f64, flank: f64, share: f64) -> HashMap<String, String> {
    let m = depths(&lpa_row(), n, flank, 0.0, 0.0);
    assess(name, Some((&m, &q1(&m, share))))
}

fn f(t: &HashMap<String, String>, k: &str) -> f64 {
    g(t, k).parse().unwrap_or_else(|_| panic!("{k}: {t:?}"))
}

fn cv(v: &[f64]) -> f64 {
    let mu = v.iter().sum::<f64>() / v.len() as f64;
    (v.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt() / mu
}

#[test]
fn the_arithmetic_a_block_at_30_n_plus_2_over_14_is_n() {
    for n in [12.0, 20.0, 40.0] {
        let t = plain(&format!("n{n}"), n, 30.0, 0.9);
        assert_eq!(
            (g(&t, "status"), g(&t, "kiv2_copies")),
            ("estimated", format!("{n:.1}").as_str()),
            "{t:?}"
        );
    }
}

#[test]
fn the_interval_by_hand() {
    let r = lpa_row();
    let t = plain("zero", 20.0, 30.0, 0.9);
    assert_eq!(
        (g(&t, "ci95_low"), g(&t, "ci95_high")),
        ("20.0", "20.0"),
        "equal depths: zero width"
    );
    let m = depths(&r, 20.0, 30.0, 0.10, 0.15);
    let t = assess("spread", Some((&m, &q1(&m, 0.9))));
    let kv: Vec<f64> = r.targets.iter().map(|s| m[&s.name]).collect();
    let fv: Vec<f64> = r.controls.iter().map(|s| m[&s.name]).collect();
    let rel = (cv(&kv).powi(2) / 14.0 + cv(&fv).powi(2) / 15.0).sqrt();
    let pooled = |segs: &[og_cn::catalog::Segment]| {
        segs.iter()
            .map(|s| m[&s.name] * (s.end - s.start) as f64)
            .sum::<f64>()
            / segs.iter().map(|s| (s.end - s.start) as f64).sum::<f64>()
    };
    let total = 14.0 * pooled(&r.targets) / pooled(&r.controls);
    let want = (
        format!("{:.1}", total - 2.0),
        format!("{:.1}", total * (-1.96 * rel).exp() - 2.0),
        format!("{:.1}", total * (1.96 * rel).exp() - 2.0),
    );
    let got = (
        g(&t, "kiv2_copies").to_owned(),
        g(&t, "ci95_low").to_owned(),
        g(&t, "ci95_high").to_owned(),
    );
    assert_eq!(got, want);
    assert!(f(&t, "ci95_low") < f(&t, "kiv2_copies") && f(&t, "kiv2_copies") < f(&t, "ci95_high"));
    let m2 = depths(&r, 20.0, 30.0, 0.20, 0.30);
    let t2 = assess("spread2", Some((&m2, &q1(&m2, 0.9))));
    assert!(
        f(&t2, "ci95_high") - f(&t2, "ci95_low") > f(&t, "ci95_high") - f(&t, "ci95_low"),
        "twice the spread"
    );
}

#[test]
fn the_class_by_the_whole_interval() {
    let r = lpa_row();
    for (n, sp, cls) in [
        (20.0, 0.03, "low"),
        (29.0, 0.03, "low"),
        (29.0, 0.25, "typical"),
        (38.0, 0.03, "typical"),
        (47.0, 0.25, "typical"),
        (60.0, 0.03, "high"),
    ] {
        let m = depths(&r, n, 30.0, sp, sp);
        let t = assess(&format!("class-{n}-{sp}"), Some((&m, &q1(&m, 0.9))));
        assert_eq!(g(&t, "class"), cls, "{n} copies, spread {sp}: {t:?}");
    }
    let t = plain("lowmsg", 20.0, 30.0, 0.9);
    assert!(
        g(&t, "message")
            .contains("a low count goes with small apo(a) isoforms and a higher Lp(a) level")
    );
}

#[test]
fn what_is_not_a_count() {
    let t = plain("8x", 20.0, 8.0, 0.9);
    assert!(
        g(&t, "status") == "too little depth"
            && g(&t, "kiv2_copies") == "."
            && g(&t, "message").contains("8.0x"),
        "{t:?}"
    );
    assert_eq!(g(&plain("10x", 20.0, 10.0, 0.9), "status"), "estimated");
    let t = plain("nomq0", 20.0, 30.0, 0.0);
    assert!(
        g(&t, "status") == "not assessed"
            && g(&t, "kiv2_copies") == "."
            && g(&t, "message").contains("novoalign"),
        "{t:?}"
    );
    assert_eq!(
        g(&plain("quarter", 20.0, 30.0, 0.25), "status"),
        "estimated"
    );
    let t = assess("nochr6", None);
    assert!(
        g(&t, "status") == "not assessed" && g(&t, "message").contains("no chr6"),
        "{t:?}"
    );
    let t = plain("floor", -1.5, 30.0, 0.9);
    assert!(
        g(&t, "status") == "not assessed"
            && g(&t, "kiv2_copies") == "."
            && g(&t, "message").contains("too little for any repeat count"),
        "{t:?}"
    );
}

#[test]
fn the_command_line() {
    let (rc, out, _) = og_cn(&["bed", "--catalog", &catalog_path(), "--locus", "LPA_KIV2"]);
    let rows: Vec<Vec<&str>> = out.lines().map(|l| l.split('\t').collect()).collect();
    assert!(rc == 0 && rows.len() == 29);
    let kiv: Vec<&Vec<&str>> = rows.iter().filter(|r| r[3].starts_with("kiv2_")).collect();
    assert_eq!(kiv.len(), 14);
    assert_eq!(
        kiv.iter()
            .map(|r| r[2].parse::<i64>().unwrap() - r[1].parse::<i64>().unwrap())
            .sum::<i64>(),
        2394
    );
    assert_eq!(
        rows.iter().filter(|r| r[3].starts_with("flank_")).count(),
        15
    );
    let (r, w) = (lpa_row(), Scratch::new("cli"));
    let (a, b, o) = (w.path("a.tsv"), w.path("b.tsv"), w.path("o.tsv"));
    let m = depths(&r, 20.0, 30.0, 0.0, 0.0);
    bedcov(&a, &r, &m, false);
    bedcov(&b, &r, &q1(&m, 0.9), false);
    let (rc, err) = estimate(&a, &b, &o);
    let t = table(&o);
    assert_eq!(rc, 0, "{err}");
    assert_eq!(
        (g(&t, "kiv2_copies"), g(&t, "status")),
        ("20.0", "estimated")
    );
    assert!(g(&t, "advice").contains("blood"));
    std::fs::write(&a, "").unwrap();
    std::fs::write(&b, "").unwrap();
    let (rc, _) = estimate(&a, &b, &o);
    assert_eq!(
        (rc, g(&table(&o), "status")),
        (0, "not assessed"),
        "no chr6"
    );
    bedcov(&a, &r, &m, false);
    let (rc, err) = estimate(&a, &b, &o);
    assert!(
        rc == 2 && err.contains("is empty but the other table is not"),
        "{rc} {err}"
    );
    bedcov(&b, &r, &q1(&m, 0.9), false);
    let mut text = read(&a);
    text += "chr6\t1\t100\tother\t500\n";
    std::fs::write(&a, text).unwrap();
    let (rc, err) = estimate(&a, &b, &o);
    assert!(
        rc == 2 && err.contains("differ"),
        "other regions: {rc} {err}"
    );
    std::fs::write(&a, "chr6\t160646213\t160646395\n").unwrap();
    let (rc, err) = estimate(&a, &b, &o);
    assert!(
        rc == 2 && err.contains("expected chrom"),
        "a short row: {rc} {err}"
    );
}

// 8. KIV-2 4925G>A

/// mpileup lines at the seven copies; `(c, t)` reads with C and T at each.
fn pileup(path: &str, c: usize, t: usize, extra: &str) {
    let mut text: String = COPIES
        .iter()
        .map(|p| {
            format!(
                "chr6\t{p}\tN\t{}\t^]{}{}$\t{}\n",
                c + t,
                "C".repeat(c),
                "t".repeat(t),
                "I".repeat(c + t)
            )
        })
        .collect();
    text += extra;
    std::fs::write(path, text).unwrap();
}

/// Step 46's table at 38 copies, a pileup, then marker-read: (exit, the
/// 4925G>A table, paths of the LPA table, the pileup and the 4925G>A table).
fn after(
    w: &Scratch,
    name: &str,
    c: usize,
    t: usize,
    extra: &str,
) -> (u8, HashMap<String, String>, [String; 3]) {
    let r = lpa_row();
    let (a, b, o, p, n) = (
        w.path(&format!("{name}/a.tsv")),
        w.path(&format!("{name}/b.tsv")),
        w.path(&format!("{name}/S_lpa_kiv2.tsv")),
        w.path(&format!("{name}/p.tsv")),
        w.path(&format!("{name}/S_lpa_kiv2_4925.tsv")),
    );
    let m = depths(&r, 38.0, 30.0, 0.0, 0.0);
    bedcov(&a, &r, &m, false);
    bedcov(&b, &r, &q1(&m, 0.9), false);
    assert_eq!(estimate(&a, &b, &o).0, 0);
    pileup(&p, c, t, extra);
    let (rc, _, _) = og_cn(&[
        "marker-read",
        "LPA_KIV2_4925",
        "--catalog",
        &catalog_path(),
        "--pileup",
        &p,
        "--copy-table",
        &o,
        "--out",
        &n,
    ]);
    (rc, table(&n), [o, p, n])
}

#[test]
fn kiv2_4925_present_reads_copies_and_the_same_table_through_call() {
    let w = Scratch::new("4925");
    let (rc, t, [o, p, n]) = after(&w, "present", 95, 5, "");
    assert_eq!(
        (rc, g(&t, "status"), g(&t, "alt_reads"), g(&t, "depth")),
        (0, "present", "35", "700")
    );
    assert_eq!(
        (g(&t, "fraction"), g(&t, "copies"), g(&t, "kiv2_copies")),
        ("0.0500", "2.0", "38.0"),
        "0.05 x 40"
    );
    assert!(g(&t, "message").contains("cannot say which chromosome or which repeat"));
    assert_eq!(
        read(&n).lines().take(3).collect::<Vec<_>>(),
        [
            "metric\tvalue",
            "marker\tLPA KIV-2 4925G>A",
            "rsid\trs1484919544 (first KIV-2 copy)"
        ]
    );
    // Step 46's table is left byte for byte.
    let r = lpa_row();
    let (a0, b0, o0) = (
        w.path("plain/a.tsv"),
        w.path("plain/b.tsv"),
        w.path("plain/S_lpa_kiv2.tsv"),
    );
    let m = depths(&r, 38.0, 30.0, 0.0, 0.0);
    bedcov(&a0, &r, &m, false);
    bedcov(&b0, &r, &q1(&m, 0.9), false);
    assert_eq!(estimate(&a0, &b0, &o0).0, 0);
    assert_eq!(read(&o), read(&o0));
    let (n3, o3) = (w.path("viacall.tsv"), w.path("viacall_lpa.tsv"));
    let (rc, out, err) = og_cn(&[
        "call",
        "--catalog",
        &catalog_path(),
        "--locus",
        "LPA_KIV2",
        "--bedcov",
        &a0,
        "--bedcov-mapq1",
        &b0,
        "--table",
        &o3,
        "--marker",
        "LPA_KIV2_4925",
        "--marker-pileup",
        &p,
        "--marker-out",
        &n3,
    ]);
    assert_eq!(rc, 0, "{err}");
    assert_eq!(
        read(&o3),
        read(&o0),
        "call --marker writes the same count table"
    );
    assert_eq!(
        read(&n3),
        read(&n),
        "and the same 4925G>A table as marker-read"
    );
    assert!(
        out.lines()
            .last()
            .unwrap_or_default()
            .starts_with("LPA LPA KIV-2 4925G>A"),
        "{out}"
    );
    // Without step 46's table the reading stands and the copies are '.'.
    let (rc, _, _) = og_cn(&[
        "marker-read",
        "LPA_KIV2_4925",
        "--catalog",
        &catalog_path(),
        "--pileup",
        &p,
        "--copy-table",
        &w.path("none.tsv"),
        "--out",
        &n3,
    ]);
    let t = table(&n3);
    assert_eq!(
        (rc, g(&t, "status"), g(&t, "copies"), g(&t, "kiv2_copies")),
        (0, "present", ".", ".")
    );
    for (share, want) in [("0.050", "not assessed"), ("0.900", "present")] {
        let lt = w.path(&format!("share{share}.tsv"));
        std::fs::write(
            &lt,
            format!("metric\tvalue\nstatus\tnot assessed\nkiv2_copies\t.\nmapq0_share\t{share}\n"),
        )
        .unwrap();
        let (rc, _, _) = og_cn(&[
            "marker-read",
            "LPA_KIV2_4925",
            "--catalog",
            &catalog_path(),
            "--pileup",
            &p,
            "--copy-table",
            &lt,
            "--out",
            &n3,
        ]);
        let t = table(&n3);
        assert_eq!((rc, g(&t, "status")), (0, want), "MAPQ 0 share {share}");
        if want != "present" {
            assert!(
                g(&t, "message").contains("left out the reads that fit several copies"),
                "{t:?}"
            );
        }
    }
}

#[test]
fn kiv2_4925_statuses_and_what_is_not_counted() {
    let w = Scratch::new("4925s");
    let (_, t, _) = after(&w, "ref", 100, 0, "");
    assert_eq!((g(&t, "status"), g(&t, "copies")), ("not seen", "."));
    assert!(g(&t, "message").contains("not proof of absence"));
    let (_, t, _) = after(&w, "unc", 100, 0, "chr6\t160640667\tN\t3\tTTT\tIII\n");
    assert_eq!((g(&t, "status"), g(&t, "alt_reads")), ("uncertain", "3"));
    let (_, t, _) = after(&w, "thin", 10, 1, "");
    assert_eq!(g(&t, "status"), "too little depth");
    let other = format!(
        "chr6\t160640668\tN\t50\t{t}\t.\nchr7\t160640667\tN\t50\t{t}\t.\nchr6\t160635123\tN\t2\tC+2TTC-1T\tII\n",
        t = "T".repeat(50)
    );
    let (_, t, _) = after(&w, "other", 100, 0, &other);
    assert_eq!(
        (g(&t, "status"), g(&t, "alt_reads"), g(&t, "depth")),
        ("not seen", "0", "702")
    );
    let (_, _, [o, p, n]) = after(&w, "bad", 95, 5, "");
    std::fs::write(&p, "not a pileup\n").unwrap();
    let (rc, _, _) = og_cn(&[
        "marker-read",
        "LPA_KIV2_4925",
        "--catalog",
        &catalog_path(),
        "--pileup",
        &p,
        "--copy-table",
        &o,
        "--out",
        &n,
    ]);
    assert_eq!(
        (rc, g(&table(&n), "status")),
        (0, "not assessed"),
        "an unreadable pileup"
    );
    let (rc, _, _) = og_cn(&[
        "marker-read",
        "LPA_KIV2_4925",
        "--catalog",
        &catalog_path(),
        "--pileup",
        &p,
        "--copy-table",
        &o,
        "--out",
        &format!("{}/no/such/dir.tsv", w.0.display()),
    ]);
    assert_eq!(rc, 2, "an output that cannot be written");
}

// 9. the demo's step 46 table

/// The string the `concat!` of string literals assigned to `const NAME` in
/// `src` (Rust source) joins to.
fn rs_str(src: &str, name: &str) -> String {
    let start = src
        .find(&format!("const {name}: &str = concat!("))
        .unwrap_or_else(|| panic!("no {name} literal"));
    let rest = &src[start..];
    let body = &rest[rest.find('(').unwrap() + 1..rest.find("\n);").unwrap()];
    body.split('"')
        .skip(1)
        .step_by(2)
        .map(|p| p.replace("\\t", "\t").replace("\\n", "\n"))
        .collect()
}

/// The demo writer `og dev demo-sample` (it replaced
/// `tests/demo/make_demo_sample.py`) writes step 46's table as a constant:
/// it must be what `og cn call` writes for the demo's depths, about 37
/// KIV-2 repeats at a flank depth near 31 (KIV-2 exons 86.5 x 1.02 and
/// x 0.98 in turn, single-copy exons 31.0 x 1.03 and x 0.97), half the
/// KIV-2 depth at MAPQ 0 (MAPQ >= 1: KIV-2 x 0.5, the rest x 0.98).
#[test]
fn the_demo_step_46_table_is_what_og_cn_call_writes_for_its_depths() {
    let src = read(&format!("{REPO}/rust/crates/og-dev/src/demo/clinical.rs"));
    let want = rs_str(&src, "LPA_TABLE");
    assert!(
        want.starts_with("metric\tvalue\nstatus\testimated\n"),
        "{want:?}"
    );
    let r = lpa_row();
    let mut m = Depth::new();
    for (i, s) in r.targets.iter().enumerate() {
        m.insert(s.name.clone(), 86.5 * if i % 2 == 0 { 1.02 } else { 0.98 });
    }
    for (i, s) in r.controls.iter().enumerate() {
        m.insert(s.name.clone(), 31.0 * if i % 2 == 0 { 1.03 } else { 0.97 });
    }
    let q: Depth = m
        .iter()
        .map(|(n, v)| {
            (
                n.clone(),
                v * if n.starts_with("kiv2_") { 0.5 } else { 0.98 },
            )
        })
        .collect();
    let w = Scratch::new("demo");
    let (a, b, o) = (w.path("a.tsv"), w.path("b.tsv"), w.path("o.tsv"));
    bedcov(&a, &r, &m, true);
    bedcov(&b, &r, &q, true);
    assert_eq!(estimate(&a, &b, &o).0, 0);
    assert_eq!(read(&o), want);
}

#[test]
fn the_advice_row() {
    let t = plain("advice", 20.0, 30.0, 0.9);
    assert_eq!(g(&t, "advice"), ADVICE);
}
