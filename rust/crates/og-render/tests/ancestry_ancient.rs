//! `og ancient` (step 44) on synthetic files, run as the pipeline runs it,
//! then the summary and both reports. It replaces
//! `tests/test_ancient_report.py`.
//!
//! 1. ranking: on a synthetic .evec the nearest ancient group is the one
//!    whose mean position is closest on the weighted PCs, a group far on PC1
//!    but with a small eigenvalue on PC3 is ranked by the weighted distance,
//!    and each row carries its dates (oldest first, AADR style BCE), country
//!    and member count;
//! 2. axes: a run whose PC2 has the opposite sign is turned back before the
//!    person is placed, and a run whose present-day coordinates differ
//!    stops with the reason instead of placing the person on other axes;
//! 3. liftover: lift keeps a site inside a block, a site on a minus-strand
//!    chain with its alleles complemented, and drops a site in a gap, one
//!    whose alleles miss the GRCh38 base, chrX and a monomorphic marker;
//! 4. genotypes: person-geno appends one column, 9 where the person has no
//!    usable genotype, and stops on a genotype file shorter than the .snp
//!    file;
//! 5. select: the present-day reference people come from MODERN_GROUPS, the
//!    ancient ones need a date, a passing assessment, enough sites and a
//!    place in West Eurasia, one genetic id per individual;
//! 6. the report command writes the table, the summary, the text and an SVG
//!    with the person marked;
//! 7. both reports show them: an HTML card with the table and the plot, and
//!    a section of the text report.
//!
//! The unit tests in og-ancestry and og-intervals check the date labels,
//! the dosages and the single-position lift in detail.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use og_ancestry::ancient::cli;
use og_render::model::Summary;
use og_render::{html, text};
use pycompat::json::{Dump, dumps};

const CAVEAT: &str = "Closeness on a PCA is genetic similarity, not descent";
const EVALS: [f64; 4] = [9.0, 4.0, 0.01, 0.01];

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// A scratch folder that is removed when the test ends, also when it fails.
struct Tmp(PathBuf);

impl std::ops::Deref for Tmp {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<Path> for Tmp {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<std::ffi::OsStr> for Tmp {
    fn as_ref(&self) -> &std::ffi::OsStr {
        self.0.as_os_str()
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn tmp(tag: &str) -> Tmp {
    let d = std::env::temp_dir().join(format!("og-ancestry-ancient-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    Tmp(d)
}

fn write(d: &Path, name: &str, text: &str) -> String {
    let p = d.join(name);
    std::fs::write(&p, text).unwrap();
    p.display().to_string()
}

fn read(p: &str) -> String {
    std::fs::read_to_string(p).unwrap_or_default()
}

fn og(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = cli(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

/// One person of an .evec: id, label, four PCs.
type Pc = (String, &'static str, [f64; 4]);

fn evec(rows: &[Pc], evals: [f64; 4]) -> String {
    let mut out = vec![format!(
        "#eigvals: {}",
        evals
            .iter()
            .map(|v| format!("{v:.3}"))
            .collect::<Vec<_>>()
            .join(" ")
    )];
    for (i, label, pcs) in rows {
        out.push(format!(
            "{i} {} {label}",
            pcs.iter()
                .map(|v| format!("{v:.6}"))
                .collect::<Vec<_>>()
                .join(" ")
        ));
    }
    out.join("\n") + "\n"
}

/// Present-day people M0..M19: they build the axes.
fn moderns() -> Vec<Pc> {
    (0..20)
        .map(|i| {
            let f = f64::from(i);
            (
                format!("M{i}"),
                if i % 2 == 1 { "French" } else { "Sardinian" },
                [
                    f * 0.01 - 0.1,
                    f64::from(i % 5) * 0.02 - 0.04,
                    f64::from(i % 3) * 0.01 - 0.01,
                    f64::from(i % 4) * 0.01 - 0.015,
                ],
            )
        })
        .collect()
}

/// The panel: the present-day people and three ancient groups.
fn panel() -> Vec<Pc> {
    let mut p = moderns();
    p.extend([
        ("A1".to_owned(), "Ancient", [0.10, 0.10, 0.0, 0.0]), // Steppe
        ("A2".to_owned(), "Ancient", [0.12, 0.10, 0.0, 0.0]),
        ("B1".to_owned(), "Ancient", [-0.20, 0.00, 0.0, 0.0]), // Farmer
        // Close on PC1 and PC2, far on PC3 (eigenvalue 0.01).
        ("C1".to_owned(), "Ancient", [0.00, 0.00, 0.9, 0.0]),
    ]);
    p
}

/// The person's run: the present-day people with PC2 flipped, as smartpca
/// may write it, and the person.
fn run_rows() -> Vec<Pc> {
    let mut r: Vec<Pc> = moderns()
        .into_iter()
        .map(|(k, l, p)| (k, l, [p[0], -p[1], p[2], p[3]]))
        .collect();
    r.push(("S".to_owned(), "Person", [0.08, -0.09, 0.0, 0.0]));
    r
}

/// The files of one report run.
struct Report {
    tsv: String,
    info: String,
    txt: String,
    svg: String,
    ancients: String,
    stats: String,
}

impl Report {
    fn new(d: &Path) -> Report {
        let mut anc = "id\tgroup\tdate_bp\tdate_sd\tfull_date\tlocality\tcountry\tlat\tlon\tsnps\n"
            .to_owned();
        for (id, group, bp, country, locality) in [
            ("A1", "Steppe_EBA", "4800", "Russia", "Samara"),
            ("A2", "Steppe_EBA", "4600", "Russia", "Samara"),
            ("B1", "Anatolia_N", "8200", "Turkey", "Barcin"),
            ("C1", "Odd_C", "3000", "Italy", "Somewhere"),
        ] {
            anc += &format!("{id}\t{group}\t{bp}\t50\tx\t{locality}\t{country}\t50\t40\t500000\n");
        }
        let o = |k: &str| d.join(format!("S_{k}")).display().to_string();
        Report {
            tsv: o("ancient.tsv"),
            info: o("info.tsv"),
            txt: o("ancient.txt"),
            svg: o("pca.svg"),
            ancients: write(d, "anc.tsv", &anc),
            stats: write(
                d,
                "stats.tsv",
                "panel_sites\t4\nsites_genotyped\t3\nother_allele\t1\nsites_used\t2\n",
            ),
        }
    }

    fn run(&self, panel_evec: &str, run_evec: &str, with_stats: bool) -> (u8, String, String) {
        let mut a = vec![
            "report",
            "--panel-evec",
            panel_evec,
            "--run-evec",
            run_evec,
            "--ancients",
            &self.ancients,
            "--sample",
            "S",
            "--release",
            "test",
        ];
        if with_stats {
            a.extend(["--stats", &self.stats]);
        }
        a.extend([
            "--out-tsv",
            &self.tsv,
            "--out-info",
            &self.info,
            "--out-txt",
            &self.txt,
            "--out-svg",
            &self.svg,
        ]);
        og(&a)
    }

    fn table(&self) -> Vec<BTreeMap<String, String>> {
        let t = read(&self.tsv);
        let mut lines = t.lines();
        let head: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
        lines
            .map(|l| {
                head.iter()
                    .zip(l.split('\t'))
                    .map(|(h, v)| ((*h).to_owned(), v.to_owned()))
                    .collect()
            })
            .collect()
    }

    fn groups(&self) -> Vec<String> {
        self.table()
            .into_iter()
            .map(|r| r["group"].clone())
            .collect()
    }
}

#[test]
fn the_nearest_group_on_the_weighted_pcs_comes_first() {
    let d = tmp("rank");
    let r = Report::new(&d);
    let pe = write(&d, "panel.evec", &evec(&panel(), EVALS));
    let re = write(&d, "run.evec", &evec(&run_rows(), EVALS));
    let (code, _, err) = r.run(&pe, &re, true);
    assert_eq!(code, 0, "{err}");
    assert_eq!(
        r.groups(),
        ["Steppe_EBA", "Odd_C", "Anatolia_N"],
        "Steppe_EBA first from the flipped run turned back; Odd_C, far only on a PC with a tiny \
         eigenvalue, second"
    );
    let rows = r.table();
    let st = &rows[0];
    assert_eq!(
        [&st["n"], &st["dates"], &st["country"], &st["locality"]],
        ["2", "2850-2650 BCE", "Russia", "Samara"]
    );
    let ranks: Vec<&str> = rows.iter().map(|x| x["rank"].as_str()).collect();
    assert_eq!(ranks, ["1", "2", "3"]);

    // Control: without the eigenvalue weights Odd_C falls behind Anatolia_N.
    let fe = write(&d, "flat_panel.evec", &evec(&panel(), [1.0; 4]));
    let fr = write(&d, "flat_run.evec", &evec(&run_rows(), [1.0; 4]));
    let (code, _, err) = r.run(&fe, &fr, true);
    assert_eq!(code, 0, "{err}");
    assert_eq!(r.groups(), ["Steppe_EBA", "Anatolia_N", "Odd_C"]);
}

#[test]
fn a_flipped_axis_is_turned_back_and_other_axes_stop() {
    let d = tmp("align");
    let pe = write(&d, "panel.evec", &evec(&panel(), EVALS));
    let re = write(&d, "run.evec", &evec(&run_rows(), EVALS));
    let (code, out, err) = og(&[
        "align",
        "--panel-evec",
        &pe,
        "--run-evec",
        &re,
        "--pcs",
        "2",
    ]);
    assert_eq!(code, 0, "{err}");
    assert!(out.contains("PC1+ PC2-"), "{out}");
    let bad: Vec<Pc> = moderns()
        .into_iter()
        .map(|(k, l, p)| {
            let n: u32 = k[1..].parse().unwrap();
            let s = if n.is_multiple_of(3) { -1.0 } else { 1.0 };
            (k, l, [p[0] * s, p[1], p[2], p[3]])
        })
        .collect();
    let be = write(&d, "bad.evec", &evec(&bad, EVALS));
    let (code, _, err) = og(&[
        "align",
        "--panel-evec",
        &pe,
        "--run-evec",
        &be,
        "--pcs",
        "2",
    ]);
    assert_eq!(code, 1, "{err}");
    assert!(err.contains("PC1"), "{err}");
}

#[test]
fn lift_keeps_sites_in_blocks_and_complements_the_minus_strand() {
    let d = tmp("lift");
    let chain = write(
        &d,
        "test.chain",
        &([
            // chr1 source 1000-1100 -> target 2000-2100 (+), gap of 50, then
            // 1150-1200 -> 2150-2200.
            "chain 100 chr1 100000 + 1000 1200 chr1 100000 + 2000 2200 1",
            "100 50 50",
            "50",
            "",
            // chr2 source 500-600 -> target on the minus strand of a 10,000 bp
            // chr2, reverse start 300.
            "chain 90 chr2 50000 + 500 600 chr2 10000 - 300 400 2",
            "100",
            "",
            // chr3: two chains over 100-200.
            "chain 80 chr3 50000 + 100 200 chr3 50000 + 100 200 3",
            "100",
            "",
            "chain 70 chr3 50000 + 150 250 chr3 50000 + 5150 5250 4",
            "100",
            "",
        ]
        .join("\n")
            + "\n"),
    );
    let mut seq1 = vec![b'N'; 3000];
    seq1[2050] = b'G'; // 1-based 2051
    seq1[2159] = b'T'; // 1-based 2160
    let mut seq2 = vec![b'N'; 10000];
    seq2[9689] = b'C'; // 1-based 9690: the complement of the source G
    let (mut fa, mut fai, mut off) = (String::new(), String::new(), 0usize);
    for (name, seq) in [("chr1", &seq1), ("chr2", &seq2)] {
        let head = format!(">{name}\n");
        fa += &head;
        off += head.len();
        for c in seq.chunks(60) {
            fa += std::str::from_utf8(c).unwrap();
            fa.push('\n');
        }
        fai += &format!("{name}\t{}\t{off}\t60\t61\n", seq.len());
        off += seq.len() + seq.len().div_ceil(60);
    }
    let fa = write(&d, "ref.fa", &fa);
    write(&d, "ref.fa.fai", &fai);
    let snp = write(
        &d,
        "src.snp",
        &([
            "rs1 1 0.0 1051 G A",  // ref G at target: kept
            "rs2 1 0.0 1160 C T",  // second allele is the target base: kept
            "rs3 1 0.0 1120 A G",  // gap: dropped
            "rs4 2 0.0 511 G T",   // minus strand: alleles C/A, target C: kept
            "rs5 1 0.0 1051 A C",  // neither allele is G: dropped
            "rs6 23 0.0 1051 G A", // X: dropped
            "rs7 1 0.0 1051 G X",  // monomorphic marker: dropped
        ]
        .join("\n")
            + "\n"),
    );
    let out = d.join("sites.tsv").display().to_string();
    let (code, _, err) = og(&[
        "lift", "--snp", &snp, "--chain", &chain, "--fasta", &fa, "--out", &out,
    ]);
    assert_eq!(code, 0, "{err}");
    assert_eq!(
        read(&out).lines().collect::<Vec<_>>(),
        [
            "snp\tchrom\tpos\tallele1\tallele2",
            "rs1\tchr1\t2051\tG\tA",
            "rs2\tchr1\t2160\tC\tT",
            "rs4\tchr2\t9690\tC\tA"
        ]
    );

    // person-geno on those sites.
    let msnp = write(
        &d,
        "moderns.snp",
        "rs1 1 0.0 1051 G A\nrs2 1 0.0 1160 C T\nrsX 23 0.0 5 A G\nrs4 2 0.0 511 G T\n",
    );
    let mgeno = write(&d, "moderns.geno", "012\n220\n999\n101\n");
    let gts = write(
        &d,
        "gts.tsv",
        "chr1\t2051\tG\tA\tG/A\nchr1\t2160\tT\t.\tT/T\nchr2\t9690\tC\tG\tC/G\n",
    );
    let merged = d.join("merged.geno").display().to_string();
    let pstats = d.join("pstats.tsv").display().to_string();
    let geno = |m: &str| {
        og(&[
            "person-geno",
            "--sites",
            &out,
            "--genotypes",
            &gts,
            "--snp",
            &msnp,
            "--moderns-geno",
            m,
            "--out",
            &merged,
            "--stats",
            &pstats,
        ])
    };
    let (code, _, err) = geno(&mgeno);
    assert_eq!(code, 0, "{err}");
    assert_eq!(
        read(&merged).lines().collect::<Vec<_>>(),
        ["0121", "2200", "9999", "1019"],
        "one column: 1, 0, 9 for X, 9 for an allele outside the pair"
    );
    let st: BTreeMap<String, String> = read(&pstats)
        .lines()
        .map(|l| {
            let (k, v) = l.split_once('\t').unwrap();
            (k.to_owned(), v.to_owned())
        })
        .collect();
    let want: BTreeMap<String, String> = [
        ("panel_sites", "4"),
        ("sites_genotyped", "3"),
        ("other_allele", "1"),
        ("sites_used", "2"),
    ]
    .map(|(k, v)| (k.to_owned(), v.to_owned()))
    .into();
    assert_eq!(st, want);
    let short = write(&d, "short.geno", "012\n220\n");
    let (code, _, err) = geno(&short);
    assert_eq!(
        code, 1,
        "a genotype file shorter than the .snp file stops: {err}"
    );
}

#[test]
fn select_takes_reference_and_ancient_people_by_the_rules() {
    let d = tmp("select");
    let head = [
        "Genetic ID (long text)",
        "Persistent Genetic ID",
        "Individual ID",
        "\"Date mean in BP in years before 1950\"",
        "Date standard deviation in BP [x]",
        "Full Date One of two formats",
        "Group ID",
        "Locality",
        "Political Entity",
        "Latitude",
        "Longitude",
        "Data type",
        "SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)",
        "SNPs hit on autosomal targets (Computed using easystats on HO snpset)",
        "ASSESSMENT",
        "ASSESSMENT WARNINGS (long)",
    ]
    .map(str::to_owned)
    .to_vec();
    #[allow(clippy::too_many_arguments)]
    fn arow(
        gid: &str,
        ind: &str,
        bp: u32,
        group: &str,
        lat: &str,
        lon: &str,
        dtype: &str,
        snps: u32,
        ass: &str,
    ) -> Vec<String> {
        [
            gid,
            gid,
            ind,
            &bp.to_string(),
            "50",
            "x",
            group,
            "Loc",
            "Land",
            lat,
            lon,
            dtype,
            &snps.to_string(),
            "0",
            ass,
            "",
        ]
        .map(str::to_owned)
        .to_vec()
    }
    let mut rows: Vec<Vec<String>> = (0..60)
        .map(|i| {
            arow(
                &format!("F{i}.DG"),
                &format!("F{i}"),
                0,
                "French",
                "46",
                "2",
                "Shotgun.diploid",
                1_100_000,
                "PROVISIONAL_PASS",
            )
        })
        .collect();
    rows.extend([
        // The same individual with fewer sites.
        arow(
            "F0.SG",
            "F0",
            0,
            "French",
            "46",
            "2",
            "Shotgun.diploid",
            900_000,
            "PASS",
        ),
        // Not a reference group.
        arow(
            "Y1.DG",
            "Y1",
            0,
            "Yakut",
            "62",
            "129",
            "Shotgun.diploid",
            1_100_000,
            "PASS",
        ),
        // Not shotgun diploid.
        arow(
            "H1.HO",
            "H1",
            0,
            "French",
            "46",
            "2",
            "HumanOrigins",
            600_000,
            "PASS",
        ),
        arow(
            "A1.AG",
            "A1",
            4800,
            "Russia_Steppe",
            "53",
            "50",
            "1240k",
            500_000,
            "PASS",
        ),
        // Too few sites.
        arow(
            "A2.AG",
            "A2",
            4800,
            "Russia_Steppe",
            "53",
            "50",
            "1240k",
            50_000,
            "PASS",
        ),
        arow(
            "A3.AG",
            "A3",
            4800,
            "Russia_Steppe",
            "53",
            "50",
            "1240k",
            500_000,
            "Questionable",
        ),
        // Outside West Eurasia.
        arow(
            "A4.AG", "A4", 4800, "China_N", "40", "116", "1240k", 500_000, "PASS",
        ),
        // No place.
        arow(
            "A5.AG", "A5", 4800, "Unknown", "..", "..", "1240k", 500_000, "PASS",
        ),
    ]);
    let anno_text: String = std::iter::once(&head)
        .chain(rows.iter())
        .map(|r| r.join("\t") + "\n")
        .collect();
    let anno = write(&d, "x.anno", &anno_text);
    let ind_text: String = rows
        .iter()
        .map(|r| format!("{} U {}\n", r[0], r[6]))
        .collect();
    let ind = write(&d, "x.ind", &ind_text);
    let o = |n: &str| d.join(n).display().to_string();
    let (oind, plist, anc) = (o("panel.ind"), o("poplist.txt"), o("ancients.tsv"));
    let (code, _, err) = og(&[
        "select",
        "--anno",
        &anno,
        "--ind",
        &ind,
        "--out-ind",
        &oind,
        "--poplist",
        &plist,
        "--ancients",
        &anc,
    ]);
    assert_eq!(code, 0, "{err}");
    let labels: Vec<(String, String)> = read(&oind)
        .lines()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].to_owned(), c[2].to_owned())
        })
        .collect();
    let label = |k: &str| labels.iter().find(|(i, _)| i == k).map(|(_, v)| v.as_str());
    assert_eq!(
        labels.iter().filter(|(_, v)| v == "French").count(),
        60,
        "the 60 French shotgun genomes are the reference"
    );
    assert_eq!(label("F0.SG"), Some("Ignore"), "one id per person (F0.DG)");
    assert_eq!(label("Y1.DG"), Some("Ignore"));
    assert_eq!(label("H1.HO"), Some("Ignore"));
    let projected: Vec<&str> = labels
        .iter()
        .filter(|(_, v)| v == "Ancient")
        .map(|(k, _)| k.as_str())
        .collect();
    assert_eq!(projected, ["A1.AG"]);
    assert_eq!(read(&plist), "French\n");
    assert!(
        read(&anc)
            .lines()
            .nth(1)
            .unwrap()
            .starts_with("A1.AG\tRussia_Steppe\t4800")
    );
}

/// The summary og summary collect writes for S, as text.
fn collect(sd: &Path) -> String {
    let input = og_summary::collect::Input {
        sample: "S".to_owned(),
        sample_dir: sd.display().to_string(),
        declared_sex: None,
        assets: repo().join("assets").display().to_string(),
        group: og_summary::collect::Group::All,
        tables: repo().join("bin").display().to_string(),
    };
    dumps(
        &og_summary::collect::post::collect(&input).unwrap(),
        &Dump::default(),
    )
    .unwrap()
}

#[test]
fn the_report_writes_its_files_and_both_reports_show_them() {
    let d = tmp("report");
    let r = Report::new(&d);
    let pe = write(&d, "panel.evec", &evec(&panel(), EVALS));
    let re = write(&d, "run.evec", &evec(&run_rows(), EVALS));
    let (code, _, err) = r.run(&pe, &re, true);
    assert_eq!(code, 0, "{err}");
    let info: BTreeMap<String, String> = read(&r.info)
        .lines()
        .skip(1)
        .map(|l| {
            let (k, v) = l.split_once('\t').unwrap();
            (k.to_owned(), v.to_owned())
        })
        .collect();
    assert_eq!(
        info["pc2"], "0.09000",
        "the person's PC2 in the panel's sign"
    );
    assert_eq!(info["sites_used"], "2");
    let svg = read(&r.svg);
    assert!(svg.starts_with("<svg"), "{svg}");
    assert!(svg.contains("#15171c"), "a dark background");
    assert!(svg.contains(">S</text>"), "the person's label");
    assert!(svg.trim_end().ends_with("</svg>"));
    let txt = read(&r.txt);
    assert!(txt.contains(CAVEAT) && txt.contains("qpAdm"), "{txt}");

    // Both reports read the outputs.
    let sd = d.join("S");
    std::fs::create_dir_all(sd.join("ancient")).unwrap();
    for (from, to) in [
        (&r.tsv, "S_ancient.tsv"),
        (&r.info, "S_ancient_info.tsv"),
        (&r.txt, "S_ancient.txt"),
        (&r.svg, "S_ancient_pca.svg"),
    ] {
        std::fs::copy(from, sd.join("ancient").join(to)).unwrap();
    }
    let s = Summary::parse(&collect(&sd)).unwrap();
    let page = html::report(&s).unwrap();
    assert!(page.contains("<h2>Ancient Populations</h2>"));
    assert!(page.contains("<td>1</td><td>Steppe_EBA</td><td>2850-2650 BCE</td><td>Russia</td>"));
    assert!(page.contains("data:image/svg+xml;base64,"), "the plot");
    let report = text::report(&s).unwrap();
    assert!(report.contains("Ancient Populations (AADR PCA, step 44)"));
    assert!(report.contains("1. Steppe_EBA ("), "{report}");
    assert!(report.contains("not descent"));
    assert!(report.contains("Plot: ancient/S_ancient_pca.svg"));

    let (code, _, err) = r.run(&pe, &pe, false);
    assert_eq!(code, 1, "a run without the person stops");
    assert!(err.contains("no row for S"), "{err}");
}
