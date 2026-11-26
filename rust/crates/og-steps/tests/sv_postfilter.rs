//! `og sv postfilter` on a made-up class 4 and 5 table: what the old
//! `tests/test_annotsv_postfilter.py` checked, by meaning, in cargo test.
//!
//! 1. Exon edges, on a made-up AnnotSV genes file: a deletion from the last
//!    base of exon 1 to the base before exon 3 of one transcript is flagged
//!    and named; one that starts 10 bp inside the intron is not; on a minus
//!    strand transcript the exons are numbered from its own start; a
//!    duplication is not a deletion ("."); without the genes file the
//!    column says "not checked".
//! 2. The population frequency: B_loss_AFmax for a deletion, B_gain_AFmax
//!    for a duplication, AF 0.05 is common, an empty one is not.
//! 3. The caller from INFO SVCALLER.
//! 4. The support from step 22's consensus, and its red control: with an
//!    empty consensus the two-caller deletion has support 1. Manta's ID
//!    `MantaDUP:2` matches SURVIVOR's `MantaDUP_2`.
//! 5. The command line writes the table, and `--consensus` adds the support.
//! 6. The rule table and the depth check (`--rules`
//!    assets/sv_rule_regions.tsv, `--sex`, `--par-bed`), and the red
//!    controls: without `--rules` the DSPP and IGH deletions stay class 4;
//!    the rule table itself has a known rule on every row, IGH and the MHC
//!    where GRCh38 has them, and at least 2,000 coding tandem repeats.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use og_steps::sv::Sex;
use og_steps::sv::postfilter::{self, Opts};

type Row = HashMap<String, String>;

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

fn asset(name: &str) -> String {
    repo()
        .join("assets")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

/// A scratch folder of its own for one test, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d =
            std::env::temp_dir().join(format!("og-sv-postfilter-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn put(&self, name: &str, text: &str) -> String {
        let p = self.0.join(name);
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// The options of `og sv postfilter --class45 C45 --out OUT` with nothing
/// else given.
fn opts(class45: &str, out: &str) -> Opts {
    Opts {
        class45: class45.to_owned(),
        genes: None,
        rules: None,
        consensus: None,
        sex: Sex::Unknown,
        par_bed: None,
        out: out.to_owned(),
    }
}

/// The written table: its header and its rows by column name.
fn table(path: &str) -> (Vec<String>, Vec<Row>) {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    let mut lines = text.lines();
    let head: Vec<String> = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .map(str::to_owned)
        .collect();
    let rows = lines
        .map(|l| {
            head.iter()
                .cloned()
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect();
    (head, rows)
}

/// Runs the step; the written rows keyed by `key`, and the summary line
/// (or the error, when it stopped).
fn run(o: &Opts, key: &str) -> (HashMap<String, Row>, String) {
    match postfilter::run(o) {
        Ok(summary) => {
            let (_, rows) = table(&o.out);
            let by = rows
                .into_iter()
                .map(|r| (r.get(key).cloned().unwrap_or_default(), r))
                .collect();
            (by, summary)
        }
        Err(e) => (HashMap::new(), format!("ERROR: {e}")),
    }
}

fn col<'a>(by: &'a HashMap<String, Row>, id: &str, c: &str) -> &'a str {
    by.get(id)
        .and_then(|r| r.get(c))
        .map_or("<missing>", String::as_str)
}

// chrom txStart txEnd strand gene tx cdsStart cdsEnd exonStarts exonEnds
const GENES: &str = "20\t1000000\t1010000\t+\tFIXG1\tNM_900001\t1000100\t1009900\t1000000,1006000,1009800,\t1000200,1006100,1010000,\n\
20\t3000000\t3005100\t-\tMING1\tNM_900005\t3000000\t3005100\t3000000,3005000,\t3000100,3005100,\n";

const HEAD: [&str; 16] = [
    "AnnotSV_ID",
    "SV_chrom",
    "SV_start",
    "SV_end",
    "SV_length",
    "SV_type",
    "ID",
    "INFO",
    "Annotation_mode",
    "Gene_name",
    "B_gain_source",
    "B_gain_AFmax",
    "B_loss_source",
    "B_loss_AFmax",
    "ACMG_class",
    "AnnotSV_ranking_score",
];

/// One class 4 or 5 row of the first table: chr20, AnnotSV's columns.
#[allow(clippy::too_many_arguments)]
fn row(
    i: u32,
    start: i64,
    end: i64,
    t: &str,
    cls: u8,
    caller: &str,
    loss_af: &str,
    loss_src: &str,
) -> String {
    let len = if t == "DEL" { start - end } else { end - start };
    [
        format!("20_{start}_{end}_{t}_1"),
        "20".into(),
        start.to_string(),
        end.to_string(),
        len.to_string(),
        t.into(),
        format!("SV{i}"),
        format!("END={end};SVTYPE={t};SVCALLER={caller}"),
        "full".into(),
        format!("GENE{i}"),
        String::new(),
        String::new(),
        loss_src.into(),
        loss_af.into(),
        cls.to_string(),
        "0.95".into(),
    ]
    .join("\t")
        + "\n"
}

fn class45() -> String {
    HEAD.join("\t")
        + "\n"
        + &row(
            1,
            1000200,
            1009800,
            "DEL",
            4,
            "manta",
            "0.05",
            "gnomAD-SV_v4.1_DEL_fixture",
        )
        + &row(2, 10100999, 10102999, "DEL", 5, "manta", "", "")
        + &row(3, 20000000, 20450000, "DEL", 4, "manta", "", "")
        + &row(4, 30000000, 30500000, "DUP", 4, "delly", "", "")
        + &row(5, 1000190, 1009800, "DEL", 4, "cnvpytor", "", "")
        + &row(6, 3000100, 3005000, "DEL", 4, "manta", "", "")
}

const CONS_HEAD: &str = "##fileformat=VCFv4.1\n##source=SURVIVOR\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\n";

fn consensus() -> String {
    CONS_HEAD.to_owned()
        + "chr20\t10101000\tMantaDEL:1\tN\t<DEL>\t.\tPASS\tSUPP=2;SUPP_VEC=11;SVLEN=-2000;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;\
CHR2=chr20;END=10103000;STRANDS=+-\tGT:TY:ID\t0/1:DEL:SV2\t0/1:DEL:DELLY2\n\
chr20\t20000010\tMantaDUP:2\tN\t<DUP>\t.\tPASS\tSUPP=2;SUPP_VEC=11;SVLEN=450000;SVTYPE=DUP;END=20450010\tGT:TY:ID\t0/1:DUP:MantaDUP_2\t0/1:DUP:DELLY3\n"
}

/// The columns of the checked table: AnnotSV's kept ones, the checks', and
/// with `--consensus` the support.
const KEEP: [&str; 10] = [
    "AnnotSV_ID",
    "ID",
    "SV_chrom",
    "SV_start",
    "SV_end",
    "SV_length",
    "SV_type",
    "Gene_name",
    "ACMG_class",
    "AnnotSV_ranking_score",
];
const ADDED: [&str; 10] = [
    "Caller",
    "Exon_edge",
    "Exon_edge_detail",
    "Population_AF",
    "Population_source",
    "Common",
    "AnnotSV_class",
    "Depth_support",
    "Rule",
    "Rule_detail",
];
const SUPPORT: [&str; 2] = ["Support", "Support_callers"];

#[test]
fn exon_edges_population_caller_and_support() {
    let w = Scratch::new("main");
    let genes = w.put("genes.RefSeq.sorted.bed", GENES);
    let c45 = w.put("S_sv_class4_5.tsv", &class45());
    let dst = w.path("checked.tsv");
    let checked = |genes: Option<&str>, consensus: Option<&str>| {
        let o = Opts {
            genes: genes.map(str::to_owned),
            consensus: consensus.map(str::to_owned),
            ..opts(&c45, &dst)
        };
        let (by, summary) = run(&o, "Gene_name");
        assert!(
            !summary.starts_with("ERROR"),
            "og sv postfilter stopped: {summary}"
        );
        by
    };

    // 1. exon edges
    let out = checked(Some(&genes), None);
    assert_eq!(
        col(&out, "GENE1", "Exon_edge"),
        "yes",
        "an intron-exact deletion (exon 1 end to exon 3 start) is on exon edges"
    );
    assert_eq!(
        col(&out, "GENE1", "Exon_edge_detail"),
        "NM_900001 (FIXG1): exon 1 and exon 3 edges",
        "it names the transcript and the exons"
    );
    assert_eq!(
        col(&out, "GENE5", "Exon_edge"),
        "no",
        "a deletion starting 10 bp inside the intron is not"
    );
    assert_eq!(
        col(&out, "GENE3", "Exon_edge"),
        "no",
        "a deletion in no transcript is not"
    );
    assert_eq!(
        col(&out, "GENE6", "Exon_edge_detail"),
        "NM_900005 (MING1): exon 1 and exon 2 edges",
        "minus strand: exons numbered from the transcript's start"
    );
    assert_eq!(
        col(&out, "GENE4", "Exon_edge"),
        ".",
        "a duplication is not checked as a deletion"
    );
    let no_genes = checked(None, None);
    assert_eq!(
        col(&no_genes, "GENE1", "Exon_edge"),
        "not checked",
        "without the genes file: not checked"
    );

    // 2. population frequency
    assert_eq!(
        (
            col(&out, "GENE1", "Population_AF"),
            col(&out, "GENE1", "Population_source")
        ),
        ("0.05", "gnomAD-SV_v4.1_DEL_fixture"),
        "the deletion's B_loss_AFmax and its source"
    );
    assert_eq!(col(&out, "GENE1", "Common"), "yes", "AF 0.05 is common");
    assert_eq!(
        (
            col(&out, "GENE2", "Population_AF"),
            col(&out, "GENE2", "Common")
        ),
        ("", "no"),
        "no benign record: no AF, not common"
    );
    assert_eq!(
        col(&out, "GENE4", "Population_AF"),
        "",
        "a duplication reads B_gain_AFmax (empty here)"
    );

    // 3. caller
    assert_eq!(
        (col(&out, "GENE2", "Caller"), col(&out, "GENE4", "Caller")),
        ("manta", "delly"),
        "the caller from INFO SVCALLER"
    );

    // 4. consensus support
    let cons = w.put("S_sv_consensus.vcf", &consensus());
    let sup = checked(Some(&genes), Some(&cons));
    assert_eq!(
        (
            col(&sup, "GENE2", "Support"),
            col(&sup, "GENE2", "Support_callers")
        ),
        ("2", "manta,delly"),
        "a deletion in the consensus: support 2, manta and delly"
    );
    assert_eq!(
        (
            col(&sup, "GENE3", "Support"),
            col(&sup, "GENE3", "Support_callers")
        ),
        ("1", "manta"),
        "a deletion only Manta made: support 1 (the consensus DUP at its place does not count)"
    );
    assert_eq!(
        (
            col(&sup, "GENE5", "Support"),
            col(&sup, "GENE5", "Support_callers")
        ),
        ("", ""),
        "CNVpytor is no column of this consensus, so its call's support is not known"
    );
    // SURVIVOR writes Manta's `MantaDUP:2` as `MantaDUP_2` in FORMAT/ID.
    let c45_dup = w.put(
        "S_sv_class4_5_dup.tsv",
        &(class45()
            + &row(7, 20000000, 20450000, "DUP", 4, "manta", "", "").replacen(
                "\tSV7\t",
                "\tMantaDUP:2\t",
                1,
            )),
    );
    let (dup, _) = run(
        &Opts {
            consensus: Some(cons.clone()),
            ..opts(&c45_dup, &w.path("dup_checked.tsv"))
        },
        "Gene_name",
    );
    assert_eq!(
        (
            col(&dup, "GENE7", "Support"),
            col(&dup, "GENE7", "Support_callers")
        ),
        ("2", "manta,delly"),
        "Manta's MantaDUP:2 is the consensus member MantaDUP_2: support 2, manta and delly"
    );
    let gone = checked(None, Some(&w.path("none.vcf")));
    assert!(
        !gone.is_empty() && gone.values().all(|r| r["Support"].is_empty()),
        "without a consensus file the support is unknown, not 1: {gone:?}"
    );

    // Red control: an empty consensus (header only) gives the two-caller
    // deletion support 1, so the support-2 check above can fail.
    let empty = w.put("empty.vcf", CONS_HEAD);
    let sup_e = checked(None, Some(&empty));
    assert_eq!(
        col(&sup_e, "GENE2", "Support"),
        "1",
        "red control: with an empty consensus the two-caller deletion has support 1"
    );
    assert!(
        sup_e
            .values()
            .all(|r| r["Support"] == if r["Caller"] == "cnvpytor" { "" } else { "1" }),
        "red control: and no row has support 2 (each Manta and Delly row 1): {sup_e:?}"
    );
}

#[test]
fn the_command_line_writes_the_table() {
    let w = Scratch::new("cli");
    let genes = w.put("genes.bed", GENES);
    let c45 = w.put("S_sv_class4_5.tsv", &class45());
    let cons = w.put("S_sv_consensus.vcf", &consensus());
    let dst = w.path("checked.tsv");
    let argv = |extra: &[&str]| -> Vec<String> {
        ["--class45", &c45, "--out", &dst]
            .iter()
            .chain(extra)
            .map(|s| (*s).to_owned())
            .collect()
    };

    assert_eq!(
        postfilter::cli(&argv(&["--genes", &genes, "--consensus", &cons])),
        0,
        "the command exits 0"
    );
    let (head, got) = table(&dst);
    assert_eq!(got.len(), 6, "one row per SV");
    let want: Vec<&str> = KEEP.iter().chain(&ADDED).chain(&SUPPORT).copied().collect();
    assert_eq!(head, want, "the columns");
    let g2 = got
        .iter()
        .find(|r| r["Gene_name"] == "GENE2")
        .expect("GENE2");
    assert_eq!(
        (g2["Support"].as_str(), g2["Support_callers"].as_str()),
        ("2", "manta,delly"),
        "the written support"
    );
    let summary = postfilter::run(&Opts {
        genes: Some(genes.clone()),
        consensus: Some(cons.clone()),
        ..opts(&c45, &dst)
    })
    .expect("og sv postfilter");
    assert!(
        summary.contains("deletions on exon edges: 2; common (benign-set AF >= 1%): 1"),
        "the summary line: {summary}"
    );

    assert_eq!(postfilter::cli(&argv(&[])), 0);
    let (head, _) = table(&dst);
    let summary = postfilter::run(&opts(&c45, &dst)).expect("og sv postfilter");
    assert!(
        !head.iter().any(|h| h == "Support")
            && summary.contains("deletions on exon edges: not checked"),
        "without --consensus no support columns, without --genes not checked: {head:?} {summary}"
    );

    std::fs::write(&c45, "").unwrap();
    assert_eq!(postfilter::cli(&argv(&["--genes", &genes])), 0);
    let text = std::fs::read_to_string(&dst).unwrap();
    assert_eq!(
        text.lines().count(),
        1,
        "an empty class 4 and 5 table gives a header and no row: {text:?}"
    );
}

// --- the rule table and the depth check --------------------------------------------------

const HEAD7: [&str; 16] = [
    "AnnotSV_ID",
    "SV_chrom",
    "SV_start",
    "SV_end",
    "SV_length",
    "SV_type",
    "Samples_ID",
    "ID",
    "INFO",
    "FORMAT",
    "S1",
    "Annotation_mode",
    "Gene_name",
    "B_loss_AFmax",
    "ACMG_class",
    "AnnotSV_ranking_score",
];

/// One row of the rule and depth table: (ID, chrom, start, end, type,
/// class, caller, gene, FORMAT, sample).
type R7<'a> = (
    &'a str,
    &'a str,
    i64,
    i64,
    &'a str,
    u8,
    &'a str,
    &'a str,
    &'a str,
    &'a str,
);

fn row7(r: &R7) -> String {
    let (id, chrom, start, end, t, cls, caller, gene, fmt, sample) = *r;
    let len = if t == "DEL" { start - end } else { end - start };
    [
        format!("{chrom}_{start}_{end}_{t}_1"),
        chrom.into(),
        start.to_string(),
        end.to_string(),
        len.to_string(),
        t.into(),
        "S1".into(),
        id.into(),
        format!("END={end};SVTYPE={t};SVCALLER={caller}"),
        fmt.into(),
        sample.into(),
        "full".into(),
        gene.into(),
        String::new(),
        cls.to_string(),
        "0.99".into(),
    ]
    .join("\t")
        + "\n"
}

const ROWS7: [R7; 19] = [
    (
        "TR_IN", "4", 87614600, 87615500, "DEL", 4, "manta", "DSPP", "GT:DHFFC", "0/1:0.45",
    ),
    (
        "TR_OUT", "4", 87614600, 87616700, "DEL", 4, "manta", "DSPP", "GT:DHFFC", "0/1:0.45",
    ),
    (
        "IGH_IN", "14", 106000000, 106100000, "DEL", 4, "manta", "IGHV3-23", "GT:DHFFC", "0/1:0.5",
    ),
    (
        "IGH_EDGE", "14", 105400000, 105700000, "DEL", 4, "manta", "TMEM121", "GT:DHFFC", "0/1:0.5",
    ),
    (
        "IGH_MOST", "14", 105550000, 105650000, "DEL", 4, "manta", "IGHA2", "GT:DHFFC", "0/1:0.5",
    ),
    (
        "MHC_IN", "6", 31000000, 31100000, "DEL", 4, "manta", "HLA-C", "GT:DHFFC", "0/1:0.5",
    ),
    (
        "ALIAS", "10", 46634800, 46635600, "DEL", 5, "manta", "RHEB", "GT:DHFFC", "0/1:0.3",
    ),
    (
        "REAL", "7", 151470000, 151480000, "DEL", 5, "manta", "RHEB", "GT:DHFFC", "0/1:0.3",
    ),
    (
        "D_LOW", "20", 40000000, 40300000, "DEL", 4, "manta", "G7", "GT:DHFFC", "0/1:0.95",
    ),
    (
        "D_OK", "20", 40000000, 40300000, "DEL", 4, "manta", "G8", "GT:DHFFC", "0/1:0.5",
    ),
    (
        "DUP_LOW", "20", 41000000, 41300000, "DUP", 4, "manta", "G9", "GT:DHBFC", "0/1:1.1",
    ),
    (
        "DELLY2", "20", 42000000, 42300000, "DEL", 4, "delly", "G10", "GT:RDCN", "0/1:2",
    ),
    (
        "DELLY1", "20", 43000000, 43300000, "DEL", 4, "delly", "G11", "GT:RDCN", "0/1:1",
    ),
    (
        "XDUP", "X", 50000000, 50300000, "DUP", 4, "delly", "G12", "GT:RDCN", "0/1:2",
    ),
    (
        "NODEPTH", "20", 44000000, 44300000, "DEL", 4, "manta", "G13", "GT", "0/1",
    ),
    (
        "CNV", "20", 45000000, 45300000, "DEL", 4, "cnvpytor", "G14", "GT", "0/1",
    ),
    (
        "DELLY_GT", "20", 46000000, 46300000, "DEL", 4, "delly", "G15", "GT:RDCN", "0/1:1",
    ),
    (
        "DELLY_NORDCN",
        "20",
        46500000,
        46800000,
        "DEL",
        4,
        "delly",
        "G17",
        "GT",
        "0/1",
    ),
    (
        "DELLY_SHORT",
        "20",
        47000000,
        47000800,
        "DEL",
        4,
        "delly",
        "G16",
        "GT",
        "0/1",
    ),
];

/// RHEB as AnnotSV's genes file would name it: once on the pseudogene
/// RHEBP1 (chr10, the region of the rule table) and once on RHEB itself
/// (chr7).
const GENES7: &str = "10\t46634910\t46635466\t+\tRHEB\tNM_990001\t46634910\t46635466\t46634910,\t46635466,\n\
7\t151466000\t151520000\t-\tRHEB\tNM_005614\t151466100\t151519900\t151466000,151519000,\t151467000,151520000,\n";

const SEX_ROWS: [R7; 10] = [
    (
        "XDEL1", "X", 50000000, 50050000, "DEL", 4, "delly", "G20", "GT:RDCN", "1/1:1",
    ),
    (
        "XDEL0", "X", 51000000, 51050000, "DEL", 4, "delly", "G21", "GT:RDCN", "1/1:0",
    ),
    (
        "XDUP2", "X", 52000000, 52050000, "DUP", 4, "delly", "G22", "GT:RDCN", "0/1:2",
    ),
    (
        "PARDUP2", "X", 1000000, 1050000, "DUP", 4, "delly", "G23", "GT:RDCN", "0/1:2",
    ),
    (
        "PARDEL1", "X", 1000000, 1050000, "DEL", 4, "delly", "G24", "GT:RDCN", "0/1:1",
    ),
    (
        "YDEL1", "Y", 20000000, 20050000, "DEL", 4, "delly", "G25", "GT:RDCN", "1/1:1",
    ),
    (
        "YDUP2", "Y", 21000000, 21050000, "DUP", 4, "delly", "G26", "GT:RDCN", "0/1:2",
    ),
    // PAR1 ends at chrX:2,781,479 (assets/par_grch38.bed): the last base in
    // it, then the first base past it.
    (
        "PAREDGE_IN",
        "X",
        2781479,
        2831479,
        "DUP",
        4,
        "delly",
        "G27",
        "GT:RDCN",
        "0/1:2",
    ),
    (
        "PAREDGE_OUT",
        "X",
        2781480,
        2831480,
        "DUP",
        4,
        "delly",
        "G28",
        "GT:RDCN",
        "0/1:2",
    ),
    (
        "BIGDUP",
        "20",
        48000000,
        48050000,
        "DUP",
        4,
        "delly",
        "G29",
        "GT:RDCN",
        "0/1:1000000",
    ),
];

fn table7(header: &[&str], rows: &[R7]) -> String {
    header.join("\t") + "\n" + &rows.iter().map(row7).collect::<String>()
}

/// The step on ROWS7 with the genes file, and with or without the rules.
fn run7(w: &Scratch, name: &str, with_rules: bool) -> (HashMap<String, Row>, String) {
    let c45 = w.put(&format!("{name}.tsv"), &table7(&HEAD7, &ROWS7));
    let genes = w.put(&format!("{name}_genes.bed"), GENES7);
    let o = Opts {
        genes: Some(genes),
        rules: with_rules.then(|| asset("sv_rule_regions.tsv")),
        ..opts(&c45, &w.path(&format!("{name}_checked.tsv")))
    };
    run(&o, "ID")
}

#[test]
fn region_rules_and_the_depth_check() {
    let w = Scratch::new("rules");
    let (got, summary) = run7(&w, "rules", true);
    let c = |id: &str, k: &str| col(&got, id, k).to_owned();

    // the region rules
    assert_eq!(
        (
            c("TR_IN", "ACMG_class"),
            c("TR_IN", "AnnotSV_class"),
            c("TR_IN", "Rule")
        ),
        ("3".into(), "4".into(), "tandem_repeat".into()),
        "an in-frame deletion inside DSPP's coding tandem repeat is class 3, not likely pathogenic"
    );
    assert_eq!(
        (c("TR_OUT", "ACMG_class"), c("TR_OUT", "Rule")),
        ("4".into(), String::new()),
        "the same deletion reaching 254 bp past the repeat stays class 4"
    );
    assert_eq!(
        (c("IGH_IN", "ACMG_class"), c("IGH_IN", "Rule")),
        ("3".into(), "ig_tr_locus".into()),
        "a deletion inside IGH is class 3 (V(D)J)"
    );
    assert_eq!(
        c("IGH_EDGE", "ACMG_class"),
        "4",
        "a deletion 39% inside IGH stays class 4"
    );
    assert_eq!(
        c("IGH_MOST", "ACMG_class"),
        "3",
        "a deletion 66% inside IGH is class 3 (gnomAD-SV's rule: over half)"
    );
    assert_eq!(
        (c("MHC_IN", "ACMG_class"), c("MHC_IN", "Rule")),
        ("3".into(), "mhc".into()),
        "a deletion inside the MHC is class 3"
    );
    assert_eq!(
        (
            c("ALIAS", "Gene_name"),
            c("ALIAS", "ACMG_class"),
            c("ALIAS", "Rule")
        ),
        ("RHEBP1".into(), "3".into(), "gene_alias".into()),
        "a class 5 deletion on the pseudogene RHEBP1 is not mapped to RHEB"
    );
    assert_eq!(
        (
            c("REAL", "Gene_name"),
            c("REAL", "ACMG_class"),
            c("REAL", "Rule")
        ),
        ("RHEB".into(), "5".into(), String::new()),
        "a class 5 deletion on RHEB itself stays RHEB and class 5"
    );

    // the depth check
    assert_eq!(
        (
            c("D_LOW", "Depth_support"),
            c("D_LOW", "ACMG_class"),
            c("D_LOW", "Rule")
        ),
        ("no".into(), "3".into(), "no_depth_support".into()),
        "a Manta deletion with DHFFC 0.95 has no depth support and is class 3"
    );
    assert_eq!(
        (c("D_OK", "Depth_support"), c("D_OK", "ACMG_class")),
        ("yes".into(), "4".into()),
        "the same deletion with DHFFC 0.5 stays class 4"
    );
    assert_eq!(
        (c("DUP_LOW", "Depth_support"), c("DUP_LOW", "ACMG_class")),
        ("no".into(), "3".into()),
        "a duplication with DHBFC 1.1 is class 3"
    );
    assert_eq!(
        (c("DELLY2", "Depth_support"), c("DELLY1", "Depth_support")),
        ("no".into(), "yes".into()),
        "Delly RDCN 2 on a deletion: no support; RDCN 1: support"
    );
    assert_eq!(
        c("XDUP", "Depth_support"),
        "yes",
        "a chrX duplication with RDCN 2 counts when the sex is not given"
    );
    assert_eq!(
        (c("NODEPTH", "Depth_support"), c("NODEPTH", "ACMG_class")),
        ("not checked".into(), "4".into()),
        "no depth value: not checked, class kept"
    );
    assert_eq!(
        c("CNV", "Depth_support"),
        "yes",
        "a CNVpytor call has depth support"
    );
    assert_eq!(
        (c("DELLY_GT", "Depth_support"), c("DELLY_GT", "ACMG_class")),
        ("yes".into(), "4".into()),
        "a Delly call step 05 added carries the RDCN that passed there (FORMAT GT:RDCN): support"
    );
    assert_eq!(
        (
            c("DELLY_NORDCN", "Depth_support"),
            c("DELLY_NORDCN", "ACMG_class")
        ),
        ("unknown".into(), "4".into()),
        "a long Delly call with no RDCN (a VCF not written by Delly 2.x) keeps its class with unknown support"
    );
    assert_eq!(
        (
            c("DELLY_SHORT", "Depth_support"),
            c("DELLY_SHORT", "ACMG_class")
        ),
        ("not checked".into(), "4".into()),
        "the same call of 800 bp had no RDCN gate: not checked, class kept"
    );

    let mut moved: Vec<&str> = got
        .iter()
        .filter(|(_, r)| r["ACMG_class"] != r["AnnotSV_class"])
        .map(|(i, _)| i.as_str())
        .collect();
    moved.sort_unstable();
    let mut want = vec![
        "TR_IN", "IGH_IN", "IGH_MOST", "MHC_IN", "ALIAS", "D_LOW", "DUP_LOW", "DELLY2",
    ];
    want.sort_unstable();
    assert_eq!(
        moved, want,
        "exactly the eight rows above are moved to class 3"
    );
    assert!(
        summary.contains("moved to class 3: 8; depth support no: 3, not checked: 2, unknown: 1"),
        "the summary line counts them: {summary}"
    );
    assert_eq!(
        (c("D_OK", "FORMAT"), c("D_OK", "S1"), c("NODEPTH", "S1")),
        ("GT:DHFFC".into(), "0/1:0.5".into(), "0/1".into()),
        "every row keeps AnnotSV's FORMAT and sample column (GT and the depth fields)"
    );

    // A sample named like a checked column stops the step.
    let clash_head: Vec<&str> = HEAD7
        .iter()
        .map(|&h| if h == "S1" { "Caller" } else { h })
        .collect();
    let clash = w.put("clash.tsv", &table7(&clash_head, &ROWS7));
    let err = postfilter::run(&opts(&clash, &w.path("clash_checked.tsv")))
        .expect_err("a sample named Caller must stop the step");
    assert!(
        err.to_string().contains("rename the sample"),
        "a sample named Caller stops the step instead of losing its GT and PS: {err}"
    );

    // Red control: without --rules the DSPP and IGH deletions stay class 4.
    let (plain, _) = run7(&w, "norules", false);
    assert_eq!(
        (
            col(&plain, "TR_IN", "ACMG_class"),
            col(&plain, "IGH_IN", "ACMG_class")
        ),
        ("4", "4"),
        "red control: without --rules the DSPP and IGH deletions stay class 4"
    );
}

#[test]
fn delly_depth_by_sex_and_par() {
    let w = Scratch::new("sex");
    let par = asset("par_grch38.bed");
    let c45 = w.put("sex.tsv", &table7(&HEAD7, &SEX_ROWS));
    let (yes, no) = (("yes", "4"), ("no", "3"));
    let ids = [
        "XDEL1",
        "XDUP2",
        "PARDUP2",
        "XDEL0",
        "PARDEL1",
        "YDEL1",
        "YDUP2",
        "PAREDGE_IN",
        "PAREDGE_OUT",
        "BIGDUP",
    ];
    for (sex, want) in [
        (Sex::Male, [no, yes, no, yes, yes, no, yes, no, yes, yes]),
        (Sex::Female, [yes, no, no, yes, yes, yes, no, no, no, yes]),
        (
            Sex::Unknown,
            [yes, yes, no, yes, yes, yes, yes, no, yes, yes],
        ),
    ] {
        let o = Opts {
            sex,
            par_bed: Some(par.clone()),
            ..opts(&c45, &w.path(&format!("sex_{sex:?}_checked.tsv")))
        };
        let (got, summary) = run(&o, "ID");
        assert!(!summary.starts_with("ERROR"), "{summary}");
        let got: Vec<(&str, &str)> = ids
            .iter()
            .map(|i| (col(&got, i, "Depth_support"), col(&got, i, "ACMG_class")))
            .collect();
        assert_eq!(
            got, want,
            "sex {sex:?}: chrX and chrY RDCN against 1 copy (male), 2 (female), either (not given); \
             inside the PAR, up to its last base, 2 copies; RDCN 1000000 is a gain ({ids:?})"
        );
    }
}

#[test]
fn the_rule_table() {
    let text = std::fs::read_to_string(asset("sv_rule_regions.tsv")).unwrap();
    let table: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    let rules: std::collections::BTreeSet<&str> = table
        .iter()
        .map(|r| r.get(3).copied().unwrap_or(""))
        .collect();
    assert!(
        table.iter().all(|r| r.len() == 7)
            && rules
                == ["gene_alias", "ig_tr_locus", "mhc", "tandem_repeat"]
                    .into_iter()
                    .collect(),
        "the rule table: every row has 7 columns and a known rule: {rules:?}"
    );
    let igh: Vec<&Vec<&str>> = table
        .iter()
        .filter(|r| r[3] == "ig_tr_locus" && r[5] == "IGH")
        .collect();
    let n = |s: &str| s.parse::<i64>().unwrap();
    assert!(
        igh.len() == 1
            && igh[0][0] == "chr14"
            && 105_500_000 < n(igh[0][1])
            && n(igh[0][1]) < n(igh[0][2])
            && n(igh[0][2]) < 107_000_000,
        "the rule table: IGH on chr14 between 105.5 and 107 Mb: {igh:?}"
    );
    let mhc: Vec<&[&str]> = table
        .iter()
        .filter(|r| r[3] == "mhc")
        .map(|r| &r[..3])
        .collect();
    assert_eq!(
        mhc,
        [["chr6", "28510119", "33480577"]],
        "the rule table: the MHC is GRCh38's 28,510,120-33,480,577"
    );
    assert!(
        table.iter().filter(|r| r[3] == "tandem_repeat").count() >= 2000,
        "the rule table: at least 2,000 coding tandem repeats"
    );
}
