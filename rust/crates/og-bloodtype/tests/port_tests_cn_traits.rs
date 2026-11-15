//! `og bloodtype report`: the blood group, platelet (HPA) and neutrophil
//! (HNA) antigens a genome predicts (step 41), through the command line.
//!
//! Each case writes a synthetic one-sample VCF (only the records it plants)
//! and, when the case has a BAM, a `samtools bedcov` table over
//! `assets/bloodtype_markers.tsv` and the RhD regions with chosen depths,
//! runs `og bloodtype report`, and checks rows of the table it writes. The
//! marker table is read here on its own, not through the crate's reader, so
//! a reader bug cannot hide a table bug.
//!
//! Red cells (this was `tests/test_bloodtype_report.py`):
//!   1. the table and the code agree (`og bloodtype check-markers`); a table
//!      missing a rule's site fails the check
//!   2. nothing planted, reads everywhere: the reference alleles at
//!      confidence high; without the BAM the same calls at medium, and RhD
//!      not determined
//!   3. ABO: A1, A2, B, AB, O.02, and markers that do not add up
//!   4. RhD from depth: two copies, one, none, the RHD pseudogene, a hybrid
//!      gene, low flank depth
//!   5. Duffy and Kell flip with planted variants
//!   6. the VCF: too few reads, a FILTER-failed call, a RefCall 0/0, contigs
//!      without chr, a multi-allelic record, a REF that is not the marker's,
//!      a second sample
//!   7. a BAM without some contigs
//!   8. each null and HPA-7 to -35 row's REF and ALT read its c. change on
//!      the gene's strand
//!   9. every such row, positive and negative
//!  10. gVCF records and reference blocks
//!  11. null alleles in detail
//!  12. common variants are not null alleles
//!  13. another ALT at a null row's position
//!  14. two null copies need no pair site
//!
//! Platelets and neutrophils (this was `tests/test_bloodtype_hpa_hna.py`):
//!  15. every HPA and HNA system has its site(s), at the GRCh38 position and
//!      REF base Ensembl gave for its dbSNP id
//!  16. each two-allele antigen for 0/0, 0/1 and 1/1
//!  17. the meaning a person who lacks an antigen gets in the note
//!  18. HNA-1: FCGR3B*01, *02 and *03, a c.233 no call, markers that
//!      disagree, c.194 that disagrees, and confidence always low

use indexmap::IndexMap;
use std::collections::{BTreeSet, HashMap};
use std::fmt::Debug;
use std::io::Write;
use std::path::PathBuf;
use std::sync::LazyLock;
use std::sync::atomic::{AtomicU64, Ordering};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

fn markers_path() -> String {
    format!("{REPO}/assets/bloodtype_markers.tsv")
}

fn cn_loci() -> String {
    format!("{REPO}/bin/cn_loci.tsv")
}

/// One marker row: its columns by name.
type Site = IndexMap<String, String>;

struct Table {
    /// The marker table's sites by name, in table order.
    sites: IndexMap<String, Site>,
    /// The RhD depth regions `og cn bed --locus RHD` prints, by name:
    /// (chrom, start, end).
    regions: IndexMap<String, (String, i64, i64)>,
}

static TABLE: LazyLock<Table> = LazyLock::new(|| {
    let text = std::fs::read_to_string(markers_path()).unwrap();
    let mut lines = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let mut sites = IndexMap::new();
    for l in lines {
        let cols: Vec<&str> = l.split('\t').collect();
        let row: Site = header
            .iter()
            .enumerate()
            .map(|(i, h)| {
                (
                    (*h).to_owned(),
                    cols.get(i).copied().unwrap_or("").to_owned(),
                )
            })
            .collect();
        sites.insert(row["name"].clone(), row);
    }
    let args: Vec<String> = ["bed", "--catalog", &cn_loci(), "--locus", "RHD"]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    assert_eq!(og_cn::cli::main(&args, &mut out, &mut err), 0);
    let regions = String::from_utf8(out)
        .unwrap()
        .lines()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                c[3].to_owned(),
                (
                    c[0].to_owned(),
                    c[1].parse().unwrap(),
                    c[2].parse().unwrap(),
                ),
            )
        })
        .collect();
    Table { sites, regions }
});

fn site(name: &str) -> &'static Site {
    TABLE
        .sites
        .get(name)
        .unwrap_or_else(|| panic!("no marker {name}"))
}

fn start(s: &Site) -> i64 {
    s["start"].parse().unwrap()
}

/// The VCF position: the BED end.
fn pos(s: &Site) -> i64 {
    s["end"].parse().unwrap()
}

fn nulls() -> BTreeSet<String> {
    TABLE
        .sites
        .iter()
        .filter(|(_, r)| r["alt_means"] == "null")
        .map(|(n, _)| n.clone())
        .collect()
}

/// Collects failed checks and reports them all at once.
#[derive(Default)]
struct Checks(Vec<String>);

impl Checks {
    fn eq<T: PartialEq + Debug>(&mut self, desc: &str, got: T, want: T) {
        if got != want {
            self.0
                .push(format!("FAIL {desc}: got {got:?}, want {want:?}"));
        }
    }
    fn done(self) {
        assert!(self.0.is_empty(), "\n{}", self.0.join("\n"));
    }
}

/// A VCF record at the marker `name`, with the options a case changes.
#[derive(Clone, Copy)]
struct Line<'a> {
    filt: &'a str,
    chrom: Option<&'a str>,
    alt: Option<&'a str>,
    ref_: Option<&'a str>,
}

const PLAIN: Line<'static> = Line {
    filt: "PASS",
    chrom: None,
    alt: None,
    ref_: None,
};

fn vcf_with(name: &str, gt: &str, o: Line<'_>) -> String {
    let r = site(name);
    [
        o.chrom.unwrap_or(&r["chrom"]),
        &pos(r).to_string(),
        ".",
        o.ref_.unwrap_or(&r["ref"]),
        o.alt.unwrap_or(&r["alt"]),
        "50",
        o.filt,
        ".",
        "GT",
        gt,
    ]
    .join("\t")
        + "\n"
}

fn vl(name: &str, gt: &str) -> String {
    vcf_with(name, gt, PLAIN)
}

fn raw_line(chrom: &str, pos: i64, ref_: &str, alt: &str, gt: &str) -> String {
    format!("{chrom}\t{pos}\t.\t{ref_}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n")
}

/// Mean depth per region name, `*` for every other one.
type Depth = HashMap<String, f64>;

fn depth_of(pairs: &[(&str, f64)]) -> Depth {
    pairs.iter().map(|(k, v)| ((*k).to_owned(), *v)).collect()
}

/// Depth means for an RHD copy number (2 = the flanks' depth), with some
/// exons at other copy numbers.
fn rhd_ex(copies: f64, flank: f64, exons: &[(i64, f64)]) -> Depth {
    let mut d = depth_of(&[
        ("*", 30.0),
        ("RHD_FLANK_UP", flank),
        ("RHD_FLANK_DOWN", flank),
        ("RHD_GENE", flank * copies / 2.0),
    ]);
    for n in TABLE.regions.keys() {
        if let Some(e) = n.strip_prefix("RHD_EX") {
            let e: i64 = e.parse().unwrap();
            let c = exons
                .iter()
                .find(|(x, _)| *x == e)
                .map_or(copies, |(_, c)| *c);
            d.insert(n.clone(), flank * c / 2.0);
        }
    }
    d
}

fn rhd(copies: f64) -> Depth {
    rhd_ex(copies, 30.0, &[])
}

fn with(mut d: Depth, name: &str, mean: f64) -> Depth {
    d.insert(name.to_owned(), mean);
    d
}

/// The shape of the marker rows of a depth table.
#[derive(Clone, Copy)]
enum Cols {
    /// The table's 16 columns, as step 41 gives samtools bedcov.
    Wide,
    /// The first six columns.
    Narrow,
}

fn depth_table(depth: &Depth, drop: &[&str], cols: Cols) -> String {
    let mean = |name: &str| {
        depth
            .get(name)
            .or_else(|| depth.get("*"))
            .copied()
            .unwrap_or(30.0)
    };
    let mut text = String::new();
    for (name, r) in &TABLE.sites {
        if drop.contains(&r["chrom"].as_str()) {
            continue; // a contig the BAM has no @SQ line for: step 41 leaves its rows out
        }
        let (s, e) = (start(r), pos(r));
        let total = (mean(name) * (e - s) as f64) as i64;
        let mut c: Vec<String> = ["chrom", "start", "end", "kind", "system", "name"]
            .iter()
            .map(|k| r[*k].clone())
            .collect();
        if let Cols::Wide = cols {
            c.extend(std::iter::repeat_n(String::new(), 10));
        }
        c.push(total.to_string());
        text += &(c.join("\t") + "\n");
    }
    for (name, (chrom, s, e)) in &TABLE.regions {
        if drop.contains(&chrom.as_str()) {
            continue;
        }
        let total = (mean(name) * (e - s) as f64) as i64;
        text += &format!("{chrom}\t{s}\t{e}\t{name}\t{total}\n");
    }
    text
}

static SEQ: AtomicU64 = AtomicU64::new(0);

/// A fresh folder for one run.
fn scratch() -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "og-bloodtype-port-{}-{}",
        std::process::id(),
        SEQ.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// A run's options beyond its records.
#[derive(Clone, Copy)]
struct Opts<'a> {
    depth: Option<&'a Depth>,
    cols: Cols,
    gz: bool,
    samples: &'a [&'a str],
    drop: &'a [&'a str],
}

const NO_BAM: Opts<'static> = Opts {
    depth: None,
    cols: Cols::Wide,
    gz: false,
    samples: &["S"],
    drop: &[],
};

/// The report's rows in table order, or the error text of a refused run.
type Rows = Result<Vec<Site>, String>;

fn report(records: &[String], o: Opts<'_>) -> Rows {
    let d = scratch();
    let text = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n{}",
        o.samples.join("\t"),
        records.concat()
    );
    let vcf = d.join(if o.gz { "s.vcf.gz" } else { "s.vcf" });
    if o.gz {
        let mut g = flate2::write::GzEncoder::new(
            std::fs::File::create(&vcf).unwrap(),
            flate2::Compression::default(),
        );
        g.write_all(text.as_bytes()).unwrap();
        g.finish().unwrap();
    } else {
        std::fs::write(&vcf, text).unwrap();
    }
    let p = |n: &str| d.join(n).to_string_lossy().into_owned();
    let mut args = vec![
        "--cn-loci".to_owned(),
        cn_loci(),
        "--sample".to_owned(),
        "S".to_owned(),
        "--vcf".to_owned(),
        vcf.to_string_lossy().into_owned(),
        "--markers".to_owned(),
        markers_path(),
        "--out".to_owned(),
        p("t.tsv"),
        "--evidence".to_owned(),
        p("e.tsv"),
    ];
    if let Some(depth) = o.depth {
        std::fs::write(d.join("d.tsv"), depth_table(depth, o.drop, o.cols)).unwrap();
        args.extend(["--depth".to_owned(), p("d.tsv")]);
    }
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_bloodtype::cli("report", &args, &mut out, &mut err);
    let result = if code == 0 {
        let t = std::fs::read_to_string(p("t.tsv")).unwrap();
        let mut lines = t.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        Ok(lines
            .map(|l| {
                header
                    .iter()
                    .zip(l.split('\t'))
                    .map(|(h, v)| ((*h).to_owned(), v.to_owned()))
                    .collect()
            })
            .collect())
    } else {
        let e = String::from_utf8_lossy(&err).trim().to_owned();
        Err(if e.is_empty() {
            format!("exit {code}")
        } else {
            e
        })
    };
    let _ = std::fs::remove_dir_all(&d);
    result
}

/// `{system: row}` of a run.
struct T(Rows);

impl T {
    fn get(&self, system: &str, col: &str) -> String {
        match &self.0 {
            Err(e) => e.clone(),
            Ok(rows) => rows
                .iter()
                .find(|r| r["system"] == system)
                .map_or_else(|| format!("no {system} row"), |r| r[col].clone()),
        }
    }
    fn ph(&self, system: &str) -> String {
        self.get(system, "phenotype")
    }
    fn cf(&self, system: &str) -> String {
        self.get(system, "confidence")
    }
    fn gt(&self, system: &str) -> String {
        self.get(system, "genotype")
    }
    fn note(&self, system: &str) -> String {
        self.get(system, "note")
    }
    fn order(&self) -> Vec<String> {
        match &self.0 {
            Ok(rows) => rows.iter().map(|r| r["system"].clone()).collect(),
            Err(e) => vec![e.clone()],
        }
    }
}

/// A run with reads (the table's 16 columns in the depth rows).
fn run(records: &[String], depth: &Depth) -> T {
    T(report(
        records,
        Opts {
            depth: Some(depth),
            ..NO_BAM
        },
    ))
}

/// A run with 30 reads everywhere and six-column depth rows.
fn run30(records: &[String]) -> T {
    static FLAT: LazyLock<Depth> = LazyLock::new(|| depth_of(&[("*", 30.0)]));
    T(report(
        records,
        Opts {
            depth: Some(&FLAT),
            cols: Cols::Narrow,
            ..NO_BAM
        },
    ))
}

/// The genotype's first part, before " (RHD copies ...)".
fn head(s: &str) -> String {
    s.split(" (").next().unwrap_or("").to_owned()
}

fn has(s: &str, part: &str) -> bool {
    s.contains(part)
}

fn check_markers(path: &str) -> (u8, String) {
    let args = [path.to_owned(), "--cn-loci".to_owned(), cn_loci()];
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_bloodtype::cli("check-markers", &args, &mut out, &mut err);
    (code, String::from_utf8_lossy(&err).into_owned())
}

#[test]
fn s01_the_marker_table_and_the_rules_agree() {
    let mut c = Checks::default();
    c.eq(
        "the marker table and the rules agree",
        check_markers(&markers_path()),
        (0, String::new()),
    );
    let d = scratch();
    let broken = d.join("m.tsv");
    let text: String = std::fs::read_to_string(markers_path())
        .unwrap()
        .split_inclusive('\n')
        .filter(|l| !l.contains("\tKEL_841\t"))
        .collect();
    std::fs::write(&broken, text).unwrap();
    let (code, err) = check_markers(&broken.to_string_lossy());
    c.eq(
        "control: a table without KEL_841 fails the check",
        (code, has(&err, "rule sites the table lacks: KEL_841")),
        (1, true),
    );
    let _ = std::fs::remove_dir_all(&d);
    c.done();
}

#[test]
fn s02_nothing_planted() {
    let mut c = Checks::default();
    let r = run(&[], &rhd(2.0));
    for (system, want) in [
        ("ABO", "O"),
        ("RhD", "D+"),
        ("RhCE", "C-c+ E-e+"),
        ("Kell", "K-k+ Kp(a-b+) Js(a-b+)"),
        ("Duffy", "Fy(a+b-)"),
        ("Kidd", "Jk(a+b-)"),
        ("MNS", "M-N+ S-s+"),
        ("Lutheran", "Lu(a-b+)"),
        ("Diego", "Di(a-b+)"),
        ("Dombrock", "Do(a-b+)"),
        ("Colton", "Co(a+b-)"),
        ("HPA-1", "HPA-1(a+b-)"),
        ("HPA-2", "HPA-2(a+b-)"),
        ("HPA-3", "HPA-3(a+b-)"),
        ("HPA-4", "HPA-4(a+b-)"),
        ("HPA-5", "HPA-5(a+b-)"),
        ("HPA-15", "HPA-15(a-b+)"),
    ] {
        c.eq(
            &format!("reference alleles, reads everywhere: {system}"),
            r.ph(system).as_str(),
            want,
        );
    }
    c.eq("ABO genotype", r.gt("ABO").as_str(), "O.01/O.01");
    c.eq(
        "ABO confidence high with reads",
        r.cf("ABO").as_str(),
        "high",
    );
    c.eq(
        "RhCE low: c.307 sits in RHCE's RHD-like exon 2",
        r.cf("RhCE").as_str(),
        "low",
    );
    c.eq("Kell high", r.cf("Kell").as_str(), "high");
    let nobam = T(report(&[], NO_BAM));
    c.eq(
        "without the BAM, ABO is still O",
        nobam.ph("ABO").as_str(),
        "O",
    );
    c.eq(
        "without the BAM, confidence medium",
        nobam.cf("ABO").as_str(),
        "medium",
    );
    c.eq(
        "without the BAM, RhD is not determined",
        nobam.ph("RhD").as_str(),
        "not determined",
    );
    c.eq(
        "without the BAM, the note says why",
        has(&nobam.note("RhD"), "deletes the whole RHD gene"),
        true,
    );
    c.done();
}

#[test]
fn s03_abo() {
    let mut c = Checks::default();
    let d2 = rhd(2.0);
    c.eq(
        "A1/O: one c.261 insertion",
        run(&[vl("ABO_261", "0/1")], &d2).ph("ABO").as_str(),
        "A1",
    );
    let a2 = run(&[vl("ABO_261", "1/1"), vl("ABO_1061", "1/1")], &d2);
    c.eq(
        "A2/A2: c.1061delC on both",
        (a2.ph("ABO"), a2.gt("ABO")),
        ("A2".into(), "A2/A2".into()),
    );
    c.eq(
        "A1/A2 reads A1",
        run(&[vl("ABO_261", "1/1"), vl("ABO_1061", "0/1")], &d2)
            .ph("ABO")
            .as_str(),
        "A1",
    );
    let b = run(
        &[
            vl("ABO_261", "0/1"),
            vl("ABO_803", "0/1"),
            vl("ABO_796", "0/1"),
        ],
        &d2,
    );
    c.eq(
        "B/O",
        (b.ph("ABO"), b.gt("ABO")),
        ("B".into(), "B/O.01".into()),
    );
    let ab = run(
        &[
            vl("ABO_261", "1/1"),
            vl("ABO_803", "0/1"),
            vl("ABO_796", "0/1"),
        ],
        &d2,
    );
    c.eq(
        "A1B",
        (ab.ph("ABO"), ab.cf("ABO")),
        ("A1B".into(), "high".into()),
    );
    c.eq(
        "the A1 or A2 caveat is in the note",
        has(&ab.note("ABO"), "c.1061delC, a deletion in a run of C"),
        true,
    );
    let o2 = run(&[vl("ABO_261", "0/1"), vl("ABO_802", "0/1")], &d2);
    c.eq(
        "O.01/O.02 is O",
        (o2.ph("ABO"), o2.gt("ABO")),
        ("O".into(), "O.01/O.02".into()),
    );
    c.eq(
        "two B markers on one functional allele: not determined",
        run(&[vl("ABO_261", "0/1"), vl("ABO_803", "1/1")], &d2)
            .ph("ABO")
            .as_str(),
        "not determined",
    );
    let odd = run(&[vl("ABO_261", "0/1"), vl("ABO_803", "0/1")], &d2);
    c.eq(
        "c.803 without c.796: B, confidence low",
        (odd.ph("ABO"), odd.cf("ABO")),
        ("B".into(), "low".into()),
    );
    c.done();
}

#[test]
fn s04_rhd_from_depth() {
    let mut c = Checks::default();
    c.eq(
        "two copies: D+",
        run(&[], &rhd(2.0)).ph("RhD").as_str(),
        "D+",
    );
    let het = run(&[], &rhd(1.0));
    c.eq(
        "one copy (one RHD deleted): D+",
        (het.ph("RhD"), head(&het.gt("RhD"))),
        ("D+".into(), "RHD/RHD deleted".into()),
    );
    let gone = run(&[], &rhd(0.0));
    c.eq(
        "no copy (RHD deleted on both): D-",
        (gone.ph("RhD"), gone.cf("RhD")),
        ("D-".into(), "high".into()),
    );
    c.eq(
        "the deleted genotype",
        head(&gone.gt("RhD")).as_str(),
        "RHD deleted/RHD deleted",
    );
    let psi = [vl("RHD_PSI_DUP", "0/1"), vl("RHD_PSI_609", "0/1")];
    c.eq(
        "one RHD copy and it is the pseudogene: D-",
        run(&psi, &rhd(1.0)).ph("RhD").as_str(),
        "D-",
    );
    c.eq(
        "two copies, one of them the pseudogene: D+",
        run(&psi, &rhd(2.0)).ph("RhD").as_str(),
        "D+",
    );
    let exons: Vec<(i64, f64)> = (3..=7).map(|e| (e, 0.0)).collect();
    let hyb = run(&[], &rhd_ex(1.6, 30.0, &exons));
    c.eq(
        "exons 3 to 7 without reads: a hybrid, RhD not determined",
        hyb.ph("RhD").as_str(),
        "not determined",
    );
    c.eq(
        "the hybrid note names the exons",
        has(&hyb.note("RhD"), "exons 3, 4, 5, 6, 7"),
        true,
    );
    c.eq(
        "flanks under 10x: not determined",
        run(&[], &rhd_ex(2.0, 6.0, &[])).ph("RhD").as_str(),
        "not determined",
    );
    c.eq(
        "near a boundary (1.45 copies): confidence low",
        run(&[], &rhd(1.45)).cf("RhD").as_str(),
        "low",
    );
    c.done();
}

#[test]
fn s05_duffy_and_kell_flip() {
    let mut c = Checks::default();
    let d2 = rhd(2.0);
    c.eq(
        "Fyb het: Fy(a+b+)",
        run(&[vl("FY_125", "0/1")], &d2).ph("Duffy").as_str(),
        "Fy(a+b+)",
    );
    let gata = run(&[vl("FY_125", "0/1"), vl("FY_GATA", "0/1")], &d2);
    c.eq(
        "Fyb het with the GATA box: Fy(a+b-)",
        (gata.ph("Duffy"), gata.gt("Duffy")),
        ("Fy(a+b-)".into(), "FY*A/FY*B GATA".into()),
    );
    c.eq(
        "Fyb hom, GATA hom: Fy(a-b-)",
        run(&[vl("FY_125", "1/1"), vl("FY_GATA", "1/1")], &d2)
            .ph("Duffy")
            .as_str(),
        "Fy(a-b-)",
    );
    c.eq(
        "FY*X: Fy(a-b+w)",
        run(&[vl("FY_125", "1/1"), vl("FY_265", "1/1")], &d2)
            .ph("Duffy")
            .as_str(),
        "Fy(a-b+w)",
    );
    let k = run(&[vl("KEL_578", "0/1")], &d2);
    c.eq(
        "K het: K+k+",
        (k.ph("Kell"), k.gt("Kell")),
        (
            "K+k+ Kp(a-b+) Js(a-b+)".into(),
            "K/k Kpb/Kpb Jsb/Jsb".into(),
        ),
    );
    c.eq(
        "K hom: K+k-",
        run(&[vl("KEL_578", "1/1")], &d2)
            .ph("Kell")
            .split(' ')
            .next()
            .unwrap_or(""),
        "K+k-",
    );
    c.eq(
        "HPA-15: the reference is 15b, the alt 15a",
        run(&[vl("HPA15_2108", "1/1")], &d2).ph("HPA-15").as_str(),
        "HPA-15(a+b-)",
    );
    c.done();
}

#[test]
fn s06_the_vcf() {
    let mut c = Checks::default();
    let d2 = rhd(2.0);
    let low = run(&[], &with(rhd(2.0), "KEL_578", 3.0));
    c.eq(
        "a missing site with 3 reads: that pair is not determined",
        (low.ph("Kell"), low.cf("Kell")),
        ("K/k not determined Kp(a-b+) Js(a-b+)".into(), "low".into()),
    );
    let lowqual = Line {
        filt: "LowQual",
        ..PLAIN
    };
    c.eq(
        "a FILTER-failed call is no call",
        run(&[vcf_with("JK_838", "1/1", lowqual)], &d2)
            .ph("Kidd")
            .as_str(),
        "not determined",
    );
    let refcall = Line {
        filt: "RefCall",
        ..PLAIN
    };
    c.eq(
        "a RefCall 0/0 is a reference call",
        run(&[vcf_with("JK_838", "0/0", refcall)], &d2)
            .ph("Kidd")
            .as_str(),
        "Jk(a+b-)",
    );
    let nochr = Line {
        chrom: Some("18"),
        ..PLAIN
    };
    let gz = T(report(
        &[vcf_with("JK_838", "1/1", nochr)],
        Opts {
            depth: Some(&d2),
            gz: true,
            ..NO_BAM
        },
    ));
    c.eq(
        "contigs without chr, gzip",
        gz.ph("Kidd").as_str(),
        "Jk(a-b+)",
    );
    let multi = run(
        &[vcf_with(
            "KEL_578",
            "1/2",
            Line {
                alt: Some("C,A"),
                ..PLAIN
            },
        )],
        &d2,
    );
    c.eq(
        "a multi-allelic record: the second ALT is K",
        multi.ph("Kell").split(' ').next().unwrap_or(""),
        "K+k+",
    );
    c.eq(
        "another allele at the site lowers the confidence",
        multi.cf("Kell").as_str(),
        "low",
    );
    let mism = run(
        &[vcf_with(
            "KEL_578",
            "1/1",
            Line {
                ref_: Some("C"),
                alt: Some("T"),
                ..PLAIN
            },
        )],
        &d2,
    );
    let words: Vec<String> = mism
        .ph("Kell")
        .split(' ')
        .take(3)
        .map(str::to_owned)
        .collect();
    c.eq(
        "a record whose REF is not the marker's REF: that pair is not determined",
        (words, mism.cf("Kell")),
        (
            vec!["K/k".into(), "not".into(), "determined".into()],
            "low".into(),
        ),
    );
    c.eq(
        "and the note names both REFs",
        has(
            &mism.note("Kell"),
            "REF C at this position, the marker's REF is G",
        ),
        true,
    );
    let span = run(
        &[vcf_with(
            "KEL_578",
            "0/0",
            Line {
                ref_: Some("GT"),
                alt: Some("G"),
                ..PLAIN
            },
        )],
        &d2,
    );
    c.eq(
        "a reference call over a longer REF that starts with the marker's REF is reference",
        span.ph("Kell").split(' ').next().unwrap_or(""),
        "K-k+",
    );
    c.eq(
        "a no-call genotype",
        run(&[vl("JK_838", "./.")], &d2).ph("Kidd").as_str(),
        "not determined",
    );
    let two = report(
        &[],
        Opts {
            samples: &["S", "T"],
            ..NO_BAM
        },
    );
    c.eq(
        "two samples stop the report",
        matches!(&two, Err(e) if e.contains("2 samples")),
        true,
    );
    c.done();
}

#[test]
fn s07_a_bam_without_some_contigs() {
    let mut c = Checks::default();
    let d2 = rhd(2.0);
    // A BAM without chr1 and chr9: no depth rows for RH, Duffy and ABO.
    let part = T(report(
        &[vl("JK_838", "1/1")],
        Opts {
            depth: Some(&d2),
            drop: &["chr1", "chr9"],
            ..NO_BAM
        },
    ));
    c.eq(
        "no depth over RHD: RhD not determined",
        part.ph("RhD").as_str(),
        "not determined",
    );
    c.eq(
        "the RhD note names the missing contig",
        has(&part.note("RhD"), "the BAM has no chr1 contig"),
        true,
    );
    c.eq(
        "ABO with no record and no depth row: not determined",
        part.ph("ABO").as_str(),
        "not determined",
    );
    c.eq(
        "the ABO note says why",
        has(&part.note("ABO"), "no depth row"),
        true,
    );
    c.eq(
        "Duffy (chr1) not determined",
        part.ph("Duffy").as_str(),
        "not determined",
    );
    c.eq(
        "Kidd (chr18, called) is read as usual",
        (part.ph("Kidd"), part.cf("Kidd")),
        ("Jk(a-b+)".into(), "high".into()),
    );
    c.eq(
        "Kell (chr7, covered) is read as usual",
        (part.ph("Kell"), part.cf("Kell")),
        ("K-k+ Kp(a-b+) Js(a-b+)".into(), "high".into()),
    );
    let contigs: Vec<&str> = TABLE
        .sites
        .values()
        .map(|r| r["chrom"].as_str())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let empty_depth = Depth::new();
    let empty = T(report(
        &[],
        Opts {
            depth: Some(&empty_depth),
            drop: &contigs,
            ..NO_BAM
        },
    ));
    c.eq(
        "an empty depth file (no marker contig in the BAM): RhD not determined",
        empty.ph("RhD").as_str(),
        "not determined",
    );
    c.eq(
        "and a system with no record is not determined",
        empty.ph("Kidd").as_str(),
        "not determined",
    );
    c.done();
}

// 8. The orientation of the null and HPA-7 to -35 rows: each row's REF and
// ALT are its c. change on the gene's strand (Ensembl GRCh38, 2026-11-03).
// A row with REF and ALT swapped, or on the wrong strand, fails here.

fn strand(gene: &str) -> i8 {
    match gene {
        "SMIM1" | "FUT2" | "SLC14A1" | "ITGB3" | "GP1BB" | "ITGA2" | "GP9" => 1,
        "ABCB6" | "ABCG2" | "FUT1" | "KEL" | "ITGA2B" => -1,
        other => panic!("no strand for gene {other}: add it to strand()"),
    }
}

fn is_bases(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b"ACGT".contains(&b))
}

/// A c. position: digits with an optional leading `-` or `*` and `+`
/// offsets (`[-*0-9+]+`).
fn is_cpos(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b"-*+0123456789".contains(&b))
}

/// `P` or `P_P`: the one or two positions of a change.
fn positions(s: &str) -> Option<(&str, Option<&str>)> {
    match s.split_once('_') {
        Some((a, b)) if is_cpos(a) && is_cpos(b) => Some((a, Some(b))),
        None if is_cpos(s) => Some((s, None)),
        _ => None,
    }
}

fn span(a: &str, b: Option<&str>) -> usize {
    let at = |x: &str| -> i64 {
        let (base, off) = x.split_once('+').unwrap_or((x, ""));
        base.parse::<i64>().unwrap()
            + if off.is_empty() {
                0
            } else {
                off.parse::<i64>().unwrap()
            }
    };
    b.map_or(1, |b| usize::try_from(at(b) - at(a) + 1).unwrap())
}

/// The VCF REF/ALT agree with one c. change on the gene's strand: the bases
/// of a substitution or a delins, the length of a deletion, duplication or
/// insertion.
fn oriented(gene: &str, ref_: &str, alt: &str, change: &str) -> bool {
    let flip = |s: &str| -> String {
        if strand(gene) == -1 {
            s.bytes()
                .rev()
                .map(|b| match b {
                    b'A' => 'T',
                    b'C' => 'G',
                    b'G' => 'C',
                    b'T' => 'A',
                    o => o as char,
                })
                .collect()
        } else {
            s.to_owned()
        }
    };
    let Some(c) = change.strip_prefix("c.") else {
        return false;
    };
    // A substitution: c.[-*0-9+_]+X>Y
    if let Some((left, y)) = c.split_once('>') {
        let Some(x) = left.chars().last().map(|x| x.to_string()) else {
            return false;
        };
        let p = &left[..left.len() - 1];
        let ok = !p.is_empty() && p.bytes().all(|b| b"-*+_0123456789".contains(&b));
        return ok
            && is_bases(&x)
            && y.len() == 1
            && is_bases(y)
            && (flip(ref_), flip(alt)) == (x, y.to_owned());
    }
    if let Some((p, bases)) = c.split_once("delins") {
        return match positions(p) {
            Some((a, b)) if is_bases(bases) => ref_.len() == span(a, b) && flip(alt) == bases,
            _ => false,
        };
    }
    for kind in ["del", "dup"] {
        if let Some((p, bases)) = c.split_once(kind) {
            let Some((a, b)) = positions(p) else {
                return false;
            };
            if !(bases.is_empty() || is_bases(bases)) {
                return false;
            }
            let n = if bases.is_empty() {
                span(a, b)
            } else {
                bases.len()
            };
            let d = if kind == "del" {
                ref_.len() as i64 - alt.len() as i64
            } else {
                alt.len() as i64 - ref_.len() as i64
            };
            if d != n as i64 || ref_.as_bytes().first() != alt.as_bytes().first() {
                return false;
            }
            // Bases given: the VCF's deleted or added bases on the gene's
            // strand are those bases, or a rotation of them (VCF left-aligns
            // in a repeat, HGVS shifts 3').
            let moved = flip(&(if kind == "del" { ref_ } else { alt })[1..]);
            return bases.is_empty() || (moved.len() == n && bases.repeat(2).contains(&moved));
        }
    }
    if let Some((p, bases)) = c.split_once("ins") {
        let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
        return match p.split_once('_') {
            Some((a, b)) if digits(a) && digits(b) && is_bases(bases) => {
                alt.len() as i64 - ref_.len() as i64 == bases.len() as i64
                    && ref_.as_bytes().first() == alt.as_bytes().first()
            }
            _ => false,
        };
    }
    false
}

/// The null rows and HPA-7 to -35 but HPA-15: one site each, the a antigen
/// on the reference allele.
fn new_rows() -> Vec<String> {
    let rare = TABLE.sites.keys().filter(|n| {
        n.strip_prefix("HPA").is_some_and(|rest| {
            rest.split('_').next().unwrap().parse::<u32>().unwrap() >= 7 && !n.starts_with("HPA15_")
        })
    });
    let mut all: BTreeSet<String> = nulls();
    all.extend(rare.cloned());
    all.into_iter().collect()
}

#[test]
fn s08_null_and_rare_hpa_rows_read_their_c_change_on_the_gene_strand() {
    let mut c = Checks::default();
    let rows = new_rows();
    assert!(rows.len() > 20, "too few null and rare HPA rows: {rows:?}");
    for n in &rows {
        let r = site(n);
        let (gene, changes) = r["hgvs"].split_once(' ').unwrap_or((&r["hgvs"], ""));
        let changes: Vec<&str> = changes.split(", ").collect();
        let alts: Vec<&str> = r["alt"].split(',').collect();
        c.eq(
            &format!("{n}: one c. change per ALT"),
            changes.len(),
            alts.len(),
        );
        let wrong: Vec<&str> = alts
            .iter()
            .zip(&changes)
            .filter(|(a, ch)| !oriented(gene, &r["ref"], a, ch))
            .map(|(_, ch)| *ch)
            .collect();
        c.eq(
            &format!("{n}: REF and ALT read the c. change on the {gene} strand"),
            wrong,
            vec![],
        );
    }
    c.eq(
        "control: a swapped REF and ALT fails the orientation",
        oriented("KEL", "A", "G", "c.578C>T"),
        false,
    );
    c.eq(
        "control: the right one passes",
        oriented("KEL", "G", "A", "c.578C>T"),
        true,
    );
    c.eq(
        "control: a deletion one base short fails",
        oriented("FUT2", "CG", "C", "c.719_720delGT"),
        false,
    );
    c.eq(
        "control: a deletion of the right length but other bases fails",
        oriented("FUT2", "CAC", "C", "c.719_720delGT"),
        false,
    );
    c.eq(
        "control: a left-aligned deletion in a repeat passes",
        oriented("FUT2", "CTG", "C", "c.719_720delGT"),
        true,
    );
    c.done();
}

#[test]
fn s09_every_null_and_rare_hpa_row_positive_and_negative() {
    let mut c = Checks::default();
    let neg = |s: &str| match s {
        "Vel" => "Vel-".to_owned(),
        "Langereis" => "Lan-".to_owned(),
        "Junior" => "Jr(a-)".to_owned(),
        "H" => "H- (Bombay, Oh)".to_owned(),
        "Secretor" => "non-secretor".to_owned(),
        "Kell" => "K-k- Kp(a-b-) Js(a-b-)".to_owned(),
        "Kidd" => "Jk(a-b-)".to_owned(),
        s => format!("{s}(a-b+)"),
    };
    let pos_ = |s: &str| match s {
        "Vel" => "Vel+".to_owned(),
        "Langereis" => "Lan+".to_owned(),
        "Junior" => "Jr(a+)".to_owned(),
        "H" => "H+".to_owned(),
        "Secretor" => "secretor".to_owned(),
        "Kell" => "K-k+ Kp(a-b+) Js(a-b+)".to_owned(),
        "Kidd" => "Jk(a+b-)".to_owned(),
        s => format!("{s}(a+b-)"),
    };
    let base = run30(&[]);
    for n in new_rows() {
        let r = site(&n);
        let s = r["system"].as_str();
        c.eq(
            &format!("{n}: the reference genome is {}", pos_(s)),
            base.ph(s),
            pos_(s),
        );
        for a in r["alt"].split(',') {
            let mut recs = vec![vcf_with(
                &n,
                "1/1",
                Line {
                    alt: Some(a),
                    ..PLAIN
                },
            )];
            if s == "H" {
                recs.push(vl("FUT2_01N02", "1/1")); // FUT1 null with a non-secretor: Bombay
            }
            c.eq(
                &format!("{n} {}>{a} on both chromosomes: {}", r["ref"], neg(s)),
                run30(&recs).ph(s),
                neg(s),
            );
        }
    }
    c.eq(
        "a reference genome lacks no high-prevalence antigen: ABO is the first row",
        base.order().into_iter().take(3).collect::<Vec<_>>(),
        vec!["ABO".to_owned(), "RhD".into(), "RhCE".into()],
    );
    c.done();
}

#[test]
fn s10_gvcf_records_and_reference_blocks() {
    let mut c = Checks::default();
    // Jk(a-b-) from SLC14A1 c.342-1G>A and HPA-21b/21b, each with the <*>
    // allele DeepVariant's gVCF adds.
    let (jk, hpa) = (site("JK_01N06"), site("HPA21_1960"));
    let gv = vec![
        format!(
            "{}\t{}\t.\tG\tA,<*>\t40\tPASS\t.\tGT:GQ\t1/1:40\n",
            jk["chrom"],
            pos(jk)
        ),
        format!(
            "{}\t{}\t.\t{}\t{},<*>\t40\tPASS\t.\tGT:GQ\t1/1:40\n",
            hpa["chrom"],
            pos(hpa),
            hpa["ref"],
            hpa["alt"]
        ),
    ];
    let g = run30(&gv);
    c.eq(
        "gVCF: c.342-1G>A 1/1 is Jk(a-b-), lacking Jk3",
        (
            g.ph("Kidd"),
            g.gt("Kidd"),
            has(&g.note("Kidd"), "lacks Jk3, carried by most people"),
        ),
        ("Jk(a-b-)".into(), "null/null".into(), true),
    );
    c.eq(
        "gVCF: HPA-21b/21b",
        (g.ph("HPA-21"), g.gt("HPA-21")),
        ("HPA-21(a-b+)".into(), "21b/21b".into()),
    );
    c.eq(
        "the systems that lack a high-prevalence antigen come first",
        g.order().into_iter().take(2).collect::<Vec<_>>(),
        vec!["Kidd".to_owned(), "HPA-21".into()],
    );
    // A reference block at the sites gives the common phenotypes.
    let blocks: Vec<String> = [jk, hpa]
        .iter()
        .map(|r| {
            format!(
                "{}\t{}\t.\t{}\t<*>\t0\t.\tEND={}\tGT:GQ\t0/0:50\n",
                r["chrom"],
                pos(r),
                &r["ref"][..1],
                pos(r) + 5
            )
        })
        .collect();
    let gb = run30(&blocks);
    c.eq(
        "gVCF reference blocks: Jk(a+b-) and HPA-21(a+b-), confidence high",
        (gb.ph("Kidd"), gb.ph("HPA-21"), gb.cf("Kidd")),
        ("Jk(a+b-)".into(), "HPA-21(a+b-)".into(), "high".into()),
    );
    c.done();
}

#[test]
fn s11_null_alleles_in_detail() {
    let mut c = Checks::default();
    let one = run30(&[vl("VEL_01N01", "0/1")]);
    c.eq(
        "one SMIM1 null copy: Vel+, Vel/null, the note names it",
        (
            one.ph("Vel"),
            one.gt("Vel"),
            has(&one.note("Vel"), "SMIM1 c.64_80del"),
        ),
        ("Vel+".into(), "Vel/null".into(), true),
    );
    let cpd = run30(&[vl("LAN_01N02", "0/1"), vl("LAN_01N13", "0/1")]);
    c.eq(
        "two different ABCB6 nulls: Lan-, confidence low (no phase), the note says so",
        (
            cpd.ph("Langereis"),
            cpd.cf("Langereis"),
            has(&cpd.note("Langereis"), "phase is not read"),
        ),
        ("Lan-".into(), "low".into(), true),
    );
    let kk = run30(&[vl("KEL_578", "0/1"), vl("KEL_02N03", "0/1")]);
    c.eq(
        "K/k with one KEL null: the null goes on the k chromosome, K+k-",
        (kk.ph("Kell"), kk.gt("Kell")),
        (
            "K+k- Kp(a-b+) Js(a-b+)".into(),
            "K/null Kpb/null Jsb/null".into(),
        ),
    );
    c.eq(
        "and that lacks k, a high-prevalence antigen, so Kell comes first",
        (
            kk.order().first().cloned().unwrap_or_default(),
            has(
                &kk.note("Kell"),
                "lacks k, carried by most people (high prevalence)",
            ),
        ),
        ("Kell".into(), true),
    );
    let jb = run30(&[vl("JK_838", "0/1"), vl("JK_02N06", "0/1")]);
    c.eq(
        "Jka/Jkb with JK*02N.06 (c.871T>C): the null goes on Jkb, Jk(a+b-)",
        (jb.ph("Kidd"), jb.gt("Kidd")),
        ("Jk(a+b-)".into(), "Jka/null".into()),
    );
    let para = run30(&[vl("FUT1_01N08", "1/1"), vl("ABO_261", "0/1")]);
    c.eq(
        "FUT1 null with a working FUT2: para-Bombay",
        para.ph("H").as_str(),
        "H- (para-Bombay)",
    );
    c.eq(
        "and the ABO row says the red cells group as O",
        has(&para.note("ABO"), "the cells usually group as O"),
        true,
    );
    let se = run30(&[vl("FUT2_01N02", "1/1")]);
    let kept: Vec<(String, String)> = se
        .order()
        .into_iter()
        .filter(|s| ["H", "ABO", "Secretor"].contains(&s.as_str()))
        .map(|s| {
            let p = se.ph(&s);
            (s, p)
        })
        .collect();
    c.eq(
        "a non-secretor alone leaves H and ABO as they are",
        kept,
        [("ABO", "O"), ("H", "H+"), ("Secretor", "non-secretor")]
            .map(|(s, p)| (s.to_owned(), p.to_owned()))
            .to_vec(),
    );
    c.eq(
        "the only Vel null row with 2 reads: not determined",
        run(&[], &with(rhd(2.0), "VEL_01N01", 2.0))
            .ph("Vel")
            .as_str(),
        "not determined",
    );
    let lowk = run(&[], &with(rhd(2.0), "KEL_02N03", 2.0));
    c.eq(
        "one KEL null row with 2 reads: Kell still read, confidence medium",
        (lowk.ph("Kell"), lowk.cf("Kell")),
        ("K-k+ Kp(a-b+) Js(a-b+)".into(), "medium".into()),
    );
    let multi_null = run30(&[vcf_with(
        "KEL_02N01",
        "1/1",
        Line {
            alt: Some("A"),
            ..PLAIN
        },
    )]);
    c.eq(
        "a null row with several ALTs counts each (c.223+1G>T, the third)",
        multi_null.ph("Kell").as_str(),
        "K-k- Kp(a-b-) Js(a-b-)",
    );
    c.done();
}

#[test]
fn s12_common_variants_are_not_null_alleles() {
    // FUT2 c.772G>A (rs602662, about half of all chromosomes in gnomAD) and
    // c.335C>T (rs200157007, 0.14 in South Asians) ride on the c.461G>A
    // haplotype, and SLC14A1 c.810G>A (rs17675299, 0.08 in Finns) is far
    // too common for a Jk null. As null rows they turned Se/se carriers into
    // non-secretors and Jk(a+b+) into Jk(a+b-) at confidence high.
    let mut c = Checks::default();
    for (desc, recs, system, want) in [
        (
            "FUT2 c.461G>A 0/1 with c.772G>A 0/1: secretor",
            vec![
                vl("FUT2_01N02", "0/1"),
                raw_line("chr19", 48_703_728, "G", "A", "0/1"),
            ],
            "Secretor",
            "secretor",
        ),
        (
            "FUT2 c.772G>A 1/1 alone: secretor",
            vec![raw_line("chr19", 48_703_728, "G", "A", "1/1")],
            "Secretor",
            "secretor",
        ),
        (
            "FUT2 c.461G>A 0/1 with c.335C>T 0/1: secretor",
            vec![
                vl("FUT2_01N02", "0/1"),
                raw_line("chr19", 48_703_291, "C", "T", "0/1"),
            ],
            "Secretor",
            "secretor",
        ),
        (
            "JK_838 0/1 with SLC14A1 c.810G>A 0/1: Jk(a+b+)",
            vec![
                vl("JK_838", "0/1"),
                raw_line("chr18", 45_739_309, "G", "A", "0/1"),
            ],
            "Kidd",
            "Jk(a+b+)",
        ),
    ] {
        c.eq(desc, run30(&recs).ph(system).as_str(), want);
    }
    let common = [
        ("chr19", 48_703_728),
        ("chr19", 48_703_291),
        ("chr18", 45_739_309),
    ];
    c.eq(
        "those three positions are no null row",
        nulls()
            .into_iter()
            .filter(|n| {
                let r = site(n);
                common.contains(&(r["chrom"].as_str(), pos(r)))
            })
            .collect::<Vec<_>>(),
        vec![],
    );
    c.done();
}

#[test]
fn s13_another_alt_at_a_null_row_is_not_that_null_allele() {
    let mut c = Checks::default();
    // KEL c.1546C>T's position, G>T: the antigens stay fully read at
    // confidence high, and the note names it.
    let kel = site("KEL_02N17");
    let oth = run30(&[raw_line(&kel["chrom"], pos(kel), "G", "T", "0/1")]);
    c.eq(
        "another ALT at a KEL null row: K-k+, confidence high, the note names it",
        (
            oth.ph("Kell"),
            oth.cf("Kell"),
            has(&oth.note("Kell"), "another allele at this position"),
        ),
        ("K-k+ Kp(a-b+) Js(a-b+)".into(), "high".into(), true),
    );
    c.done();
}

#[test]
fn s14_two_null_copies_need_no_pair_site() {
    let mut c = Checks::default();
    let nc = run30(&[vl("JK_01N06", "1/1"), vl("JK_838", "./.")]);
    c.eq(
        "Jk null 1/1 with c.838 a no call: Jk(a-b-), null/null",
        (nc.ph("Kidd"), nc.gt("Kidd")),
        ("Jk(a-b-)".into(), "null/null".into()),
    );
    let kn = run30(&[vl("KEL_02N03", "1/1"), vl("KEL_578", "./.")]);
    c.eq(
        "KEL null 1/1 with c.578 a no call: K-k- Kp(a-b-) Js(a-b-), null/null",
        (kn.ph("Kell"), kn.gt("Kell")),
        ("K-k- Kp(a-b-) Js(a-b-)".into(), "null/null".into()),
    );
    let two = run30(&[
        vl("JK_838", "0/1"),
        vl("JK_02N21", "0/1"),
        vl("JK_01N32", "0/1"),
    ]);
    c.eq(
        "two Kidd nulls with different null_on: the note names each allele's placement",
        (
            has(&two.note("Kidd"), "SLC14A1 c.118G>A: Jkb"),
            has(&two.note("Kidd"), "SLC14A1 c.151+1G>C: Jka"),
        ),
        (true, true),
    );
    c.done();
}

// Platelets (HPA) and neutrophils (HNA).

/// (system, chrom, GRCh38 position, REF, ALT, rsid)
type Want = (
    &'static str,
    &'static str,
    i64,
    &'static str,
    &'static str,
    &'static str,
);

/// name: [`Want`]
const HPA_HNA_SITES: [(&str, Want); 16] = [
    (
        "HPA1_176",
        ("HPA-1", "chr17", 47_283_364, "T", "C", "rs5918"),
    ),
    (
        "HPA2_482",
        ("HPA-2", "chr17", 4_933_086, "C", "T", "rs6065"),
    ),
    (
        "HPA3_2621",
        ("HPA-3", "chr17", 44_375_697, "A", "C", "rs5911"),
    ),
    (
        "HPA4_506",
        ("HPA-4", "chr17", 47_284_587, "G", "A", "rs5917"),
    ),
    (
        "HPA5_1600",
        ("HPA-5", "chr5", 53_062_927, "G", "A", "rs1801106"),
    ),
    (
        "HPA6_1544",
        ("HPA-6", "chr17", 47_292_422, "G", "A", "rs13306487"),
    ),
    (
        "HPA15_2108",
        ("HPA-15", "chr6", 73_783_709, "A", "C", "rs10455097"),
    ),
    (
        "FCGR3B_108",
        ("HNA-1", "chr1", 161_629_989, "G", "C", "rs200688856"),
    ),
    (
        "FCGR3B_114",
        ("HNA-1", "chr1", 161_629_983, "A", "G", "rs527909462"),
    ),
    (
        "FCGR3B_194",
        ("HNA-1", "chr1", 161_629_903, "T", "C", "rs448740"),
    ),
    (
        "FCGR3B_233",
        ("HNA-1", "chr1", 161_629_864, "G", "T", "rs5030738"),
    ),
    (
        "FCGR3B_244",
        ("HNA-1", "chr1", 161_629_853, "T", "C", "rs147574249"),
    ),
    (
        "FCGR3B_316",
        ("HNA-1", "chr1", 161_629_781, "T", "C", "rs2290834"),
    ),
    (
        "HNA3_461",
        ("HNA-3", "chr19", 10_631_494, "A", "G", "rs2288904"),
    ),
    (
        "HNA4_230",
        ("HNA-4", "chr16", 31_265_490, "G", "A", "rs1143679"),
    ),
    (
        "HNA5_2372",
        ("HNA-5", "chr16", 30_506_720, "G", "C", "rs2230433"),
    ),
];

fn is_hpa_hna(system: &str) -> bool {
    system.starts_with("HPA-") || system.starts_with("HNA-")
}

#[test]
fn s15_every_hpa_and_hna_site_is_at_its_grch38_position() {
    let mut c = Checks::default();
    for (n, w) in HPA_HNA_SITES {
        let got = TABLE.sites.get(n).map(|r| {
            (
                r["system"].clone(),
                r["chrom"].clone(),
                pos(r),
                r["ref"].clone(),
                r["alt"].clone(),
                r["rsid"].clone(),
            )
        });
        let want = (
            w.0.to_owned(),
            w.1.to_owned(),
            w.2,
            w.3.to_owned(),
            w.4.to_owned(),
            w.5.to_owned(),
        );
        c.eq(&format!("{n} at its GRCh38 position"), got, Some(want));
    }
    let (mut out, mut err) = (Vec::new(), Vec::new());
    c.eq(
        "og bloodtype systems exits 0",
        og_bloodtype::cli("systems", &[], &mut out, &mut err),
        0,
    );
    let out = String::from_utf8(out).unwrap();
    let systems: Vec<&str> = out.split('\n').collect();
    for s in ["HPA-6", "HNA-1", "HNA-3", "HNA-4", "HNA-5"] {
        c.eq(
            &format!("{s} is a system of the table"),
            systems.contains(&s),
            true,
        );
    }
    let names = |want_meaning: bool| -> Vec<&String> {
        TABLE
            .sites
            .iter()
            .filter(|(_, r)| is_hpa_hna(&r["system"]) == want_meaning)
            .filter(|(_, r)| r["meaning"].is_empty() == want_meaning)
            .map(|(n, _)| n)
            .collect()
    };
    c.eq("every HPA and HNA site has a meaning", names(true), vec![]);
    c.eq("no red cell site has a meaning", names(false), vec![]);
    c.eq(
        "every FCGR3B site carries the FCGR3A and copy-number caveat",
        HPA_HNA_SITES
            .iter()
            .map(|(n, _)| *n)
            .filter(|n| {
                n.starts_with("FCGR3B")
                    && !TABLE
                        .sites
                        .get(*n)
                        .is_some_and(|r| r["caveat"].contains("copy number"))
            })
            .collect::<Vec<_>>(),
        vec![],
    );
    c.done();
}

#[test]
fn s16_each_two_allele_antigen_at_0_0_0_1_and_1_1() {
    let mut c = Checks::default();
    // name: (system, phenotype at 0/0, 0/1, 1/1, genotype at 0/1)
    let pairs = [
        (
            "HPA1_176",
            (
                "HPA-1",
                "HPA-1(a+b-)",
                "HPA-1(a+b+)",
                "HPA-1(a-b+)",
                "1a/1b",
            ),
        ),
        (
            "HPA2_482",
            (
                "HPA-2",
                "HPA-2(a+b-)",
                "HPA-2(a+b+)",
                "HPA-2(a-b+)",
                "2a/2b",
            ),
        ),
        (
            "HPA3_2621",
            (
                "HPA-3",
                "HPA-3(a+b-)",
                "HPA-3(a+b+)",
                "HPA-3(a-b+)",
                "3a/3b",
            ),
        ),
        (
            "HPA4_506",
            (
                "HPA-4",
                "HPA-4(a+b-)",
                "HPA-4(a+b+)",
                "HPA-4(a-b+)",
                "4a/4b",
            ),
        ),
        (
            "HPA5_1600",
            (
                "HPA-5",
                "HPA-5(a+b-)",
                "HPA-5(a+b+)",
                "HPA-5(a-b+)",
                "5a/5b",
            ),
        ),
        (
            "HPA6_1544",
            (
                "HPA-6",
                "HPA-6(a+b-)",
                "HPA-6(a+b+)",
                "HPA-6(a-b+)",
                "6a/6b",
            ),
        ),
        (
            "HPA15_2108",
            (
                "HPA-15",
                "HPA-15(a-b+)",
                "HPA-15(a+b+)",
                "HPA-15(a+b-)",
                "15a/15b",
            ),
        ),
        (
            "HNA3_461",
            (
                "HNA-3",
                "HNA-3(a-b+)",
                "HNA-3(a+b+)",
                "HNA-3(a+b-)",
                "3a/3b",
            ),
        ),
        (
            "HNA4_230",
            (
                "HNA-4",
                "HNA-4(a+b-)",
                "HNA-4(a+b+)",
                "HNA-4(a-b+)",
                "4a/4b",
            ),
        ),
        (
            "HNA5_2372",
            (
                "HNA-5",
                "HNA-5(a+b-)",
                "HNA-5(a+b+)",
                "HNA-5(a-b+)",
                "5a/5b",
            ),
        ),
    ];
    let none = run30(&[]);
    for (n, (system, p0, p1, p2, g1)) in pairs {
        let het = run30(&[vl(n, "0/1")]);
        let hom = run30(&[vl(n, "1/1")]);
        c.eq(
            &format!("{system}: 0/0, 0/1, 1/1"),
            (none.ph(system), het.ph(system), hom.ph(system)),
            (p0.into(), p1.into(), p2.into()),
        );
        c.eq(
            &format!("{system}: the 0/1 genotype at confidence high"),
            (het.gt(system), het.cf(system)),
            (g1.into(), "high".into()),
        );
        c.eq(
            &format!("{system}: both antigens, no meaning in the note"),
            has(&het.note(system), "lacks "),
            false,
        );
    }
    c.done();
}

#[test]
fn s17_a_person_who_lacks_an_antigen_gets_its_meaning() {
    let mut c = Checks::default();
    let none = run30(&[]);
    let mum = run30(&[vl("HPA1_176", "1/1")]);
    let note = mum.note("HPA-1");
    c.eq(
        "HPA-1b/1b: the note says the person lacks HPA-1a",
        note.starts_with("lacks HPA-1a: "),
        true,
    );
    c.eq(
        "and names FNAIT and anti-HPA-1a",
        has(&note, "FNAIT") && has(&note, "anti-HPA-1a"),
        true,
    );
    c.eq(
        "HPA-1a/1a lacks HPA-1b",
        none.note("HPA-1").starts_with("lacks HPA-1b: "),
        true,
    );
    let hna3 = none.note("HNA-3");
    c.eq(
        "GRCh38 carries HNA-3b: no record lacks HNA-3a, with the neutrophil meaning",
        (
            hna3.starts_with("lacks HNA-3a: "),
            has(&hna3, "neonatal alloimmune neutropenia"),
        ),
        (true, true),
    );
    c.eq(
        "a red cell system never gets a meaning",
        has(&mum.note("Kidd"), "lacks "),
        false,
    );
    c.done();
}

#[test]
fn s18_hna_1() {
    let mut c = Checks::default();
    let a4 = ["FCGR3B_108", "FCGR3B_114", "FCGR3B_244", "FCGR3B_316"];
    let all = |gt: &str, extra: &[(&str, &str)]| -> Vec<String> {
        a4.iter()
            .map(|n| vl(n, gt))
            .chain(extra.iter().map(|(n, g)| vl(n, g)))
            .collect()
    };
    let b = run30(&[]);
    c.eq(
        "the GRCh38 bases (FCGR3B*02 at four markers): HNA-1(a-b+c-), FCGR3B*02 twice",
        (b.ph("HNA-1"), b.gt("HNA-1")),
        ("HNA-1(a-b+c-)".into(), "FCGR3B*02/FCGR3B*02".into()),
    );
    c.eq(
        "confidence is always low (FCGR3A, copy number)",
        b.cf("HNA-1").as_str(),
        "low",
    );
    c.eq(
        "c.194 that does not match is named in the note",
        has(&b.note("HNA-1"), "FCGR3B c.194A>G (x0) does not match"),
        true,
    );
    c.eq(
        "the note lacks HNA-1a and HNA-1c",
        has(&b.note("HNA-1"), "lacks HNA-1a, HNA-1c: "),
        true,
    );
    let ab = run30(&all("0/1", &[("FCGR3B_194", "0/1")]));
    c.eq(
        "HNA-1a/1b",
        (ab.ph("HNA-1"), ab.gt("HNA-1")),
        ("HNA-1(a+b+c-)".into(), "FCGR3B*01/FCGR3B*02".into()),
    );
    c.eq(
        "c.194 that matches: no note about it",
        has(&ab.note("HNA-1"), "does not match"),
        false,
    );
    let aa = run30(&all("1/1", &[("FCGR3B_194", "1/1")]));
    c.eq(
        "HNA-1a/1a",
        (aa.ph("HNA-1"), aa.gt("HNA-1")),
        ("HNA-1(a+b-c-)".into(), "FCGR3B*01/FCGR3B*01".into()),
    );
    c.eq(
        "HNA-1a/1a lacks HNA-1b",
        has(&aa.note("HNA-1"), "lacks HNA-1b, HNA-1c: "),
        true,
    );
    let cc = run30(&[vl("FCGR3B_233", "0/1"), vl("FCGR3B_194", "1/1")]);
    c.eq(
        "HNA-1c on an FCGR3B*02: HNA-1(a-b+c+), FCGR3B*02/FCGR3B*03",
        (cc.ph("HNA-1"), cc.gt("HNA-1")),
        ("HNA-1(a-b+c+)".into(), "FCGR3B*02/FCGR3B*03".into()),
    );
    let ac = run30(&all("1/1", &[("FCGR3B_194", "1/1"), ("FCGR3B_233", "0/1")]));
    c.eq(
        "HNA-1c with two FCGR3B*01: b+ (FCGR3B*03 carries HNA-1b) and the extra-copy note",
        (ac.ph("HNA-1"), has(&ac.note("HNA-1"), "extra FCGR3B copy")),
        ("HNA-1(a+b+c+)".into(), true),
    );
    let nc = run30(&[vl("FCGR3B_233", "./."), vl("FCGR3B_194", "1/1")]);
    c.eq(
        "c.233 with no call: the HNA-1b alleles read *02 or *03, HNA-1c not determined",
        (
            nc.ph("HNA-1"),
            nc.gt("HNA-1"),
            has(&nc.note("HNA-1"), "HNA-1c not determined"),
        ),
        (
            "HNA-1(a-b+)".into(),
            "FCGR3B*02 or *03/FCGR3B*02 or *03".into(),
            true,
        ),
    );
    let mixed = run30(&[vl("FCGR3B_108", "1/1"), vl("FCGR3B_114", "0/1")]);
    c.eq(
        "markers that disagree: not determined",
        mixed.ph("HNA-1").as_str(),
        "not determined",
    );
    c.eq(
        "and the note says why",
        has(&mixed.note("HNA-1"), "do not agree on HNA-1a or HNA-1b"),
        true,
    );
    c.done();
}
