//! The summary sections of the opt-in steps: what `og summary collect`
//! reads from the tables of `og chip`, `og viral`, `og mei`, `og mocha`,
//! `og vntyper`, `og locityper` and `og exomiser`. Moved from the summary
//! checks of tests/test_{chip,viral,mei,mocha,vntyper,locityper,exomiser}_report.py.
//!
//! Each case writes the tables as the step's tool writes them (the tools'
//! own tests check those columns) into a sample folder, then collects the
//! whole summary. The CHIP case runs `og chip report` itself.

use og_summary::collect::{Group, Input, post};
use pycompat::json::Value;
use std::fs;
use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!("og-optin-summary-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn p(&self, n: &str) -> String {
        self.0.join(n).to_str().unwrap().to_owned()
    }
    /// Writes `text` to `rel`, making its folder.
    fn put(&self, rel: &str, text: &str) -> String {
        let p = self.0.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(&p, text).unwrap();
        p.to_str().unwrap().to_owned()
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The whole summary of sample `s` in folder `dir`.
fn collect(dir: &str, s: &str) -> Value {
    post::collect(&Input {
        sample: s.into(),
        sample_dir: dir.into(),
        declared_sex: None,
        assets: format!("{REPO}/assets"),
        group: Group::All,
        tables: format!("{REPO}/bin"),
    })
    .unwrap()
}

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no {k} in {v:?}")),
        Value::Array(a) => &a[k.parse::<usize>().unwrap()],
        _ => panic!("{k} of {v:?}"),
    })
}

fn text(v: &Value) -> String {
    match v {
        Value::Str(s) => s.clone(),
        Value::Int(n) => n.to_string(),
        Value::Float(f) => f.to_string(),
        Value::Null => "null".into(),
        _ => format!("{v:?}"),
    }
}

fn list(v: &Value) -> &[Value] {
    match v {
        Value::Array(a) => a,
        _ => panic!("not a list: {v:?}"),
    }
}

/// The texts of a list of objects at `key`.
fn column(v: &Value, key: &str) -> Vec<String> {
    list(v).iter().map(|x| text(at(x, &[key]))).collect()
}

fn texts(v: &Value) -> Vec<String> {
    list(v).iter().map(text).collect()
}

/// The section's state.
fn state(s: &Value, key: &str) -> String {
    text(at(s, &["sections", key, "state"]))
}

fn data<'a>(s: &'a Value, key: &str) -> &'a Value {
    at(s, &["sections", key, "data"])
}

/// The section has no file: no source and no data.
fn missing(s: &Value, key: &str) -> bool {
    matches!(at(s, &["sections", key, "source"]), Value::Null)
        && matches!(data(s, key), Value::Object(m) if m.is_empty())
}

fn not_assessed(s: &Value) -> Vec<String> {
    texts(at(s, &["not_assessed"]))
}

fn now() -> f64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs_f64()
}

// --- step 49, CHIP ----------------------------------------------------------------

const CHIP_COLUMNS: &str = "gene\tchrom\tpos\tref\talt\tconsequence\tprotein_change\texon\tvaf\talt_reads\tdepth\tfilter\tstatus\tclass\tnote\n";
const CHIP_DEPTH: &str = "gene\ttranscript\tchrom\tcoding_bases\tmean_depth\tdetection_limit_vaf\tstatus\tcallable_fraction\tassessable_fraction\n";

#[test]
fn chip_reads_the_call_the_median_limit_and_the_uncovered_genes() {
    let d = Dir::new("chip");
    let l30 = og_chip::binom::detection_limit(30.0).unwrap().unwrap();
    d.put(
        "S/chip/S_chip.tsv",
        &format!(
            "{CHIP_COLUMNS}DNMT3A\tchr2\t1150\tC\tT\tmissense\tR882H\t2\t0.130\t6\t46\tPASS\tchip\tmissense\thotspot R882\n"
        ),
    );
    d.put(
        "S/chip/S_chip_depth.tsv",
        &format!(
            "{CHIP_DEPTH}DNMT3A\tT2\tchr2\t392\t30.0\t{l30:.3}\tcovered\nJAK2\tT3\tchr9\t0\t\t\tnot_on_reference\n"
        ),
    );
    let s = collect(&d.p("S"), "S");
    let c = data(&s, "chip");
    assert_eq!(column(at(c, &["calls"]), "protein_change"), ["R882H"]);
    let Value::Float(m) = at(c, &["median_limit"]) else {
        panic!("{c:?}")
    };
    assert!((m - 0.137).abs() < 0.001, "{m}");
    assert_eq!(texts(at(c, &["not_covered"])), ["JAK2"]);

    // With an even count of covered genes the median depth is the mean of
    // the two middle depths, and the limit is the median of the genes' own
    // limits, not the limit at the median depth (#411).
    let lim = |x: f64| og_chip::binom::detection_limit(x).unwrap().unwrap();
    let rows: String = [20.0, 30.0, 60.0, 80.0]
        .iter()
        .map(|x| format!("G{x}\tT\tchr1\t10\t{x:.1}\t{:.4}\tcovered\n", lim(*x)))
        .collect();
    d.put("E/chip/E_chip.tsv", CHIP_COLUMNS);
    d.put("E/chip/E_chip_depth.tsv", &format!("{CHIP_DEPTH}{rows}"));
    let s = collect(&d.p("E"), "E");
    let c = data(&s, "chip");
    let four = |x: f64| format!("{:.4}", lim(x)).parse::<f64>().unwrap();
    let mid = f64::midpoint(four(30.0), four(60.0));
    assert_eq!(text(at(c, &["median_depth"])), "45");
    let Value::Float(m) = at(c, &["median_limit"]) else {
        panic!("{c:?}")
    };
    assert!((m - mid).abs() < 1e-9, "{m} {mid}");
}

#[test]
fn chip_reads_what_og_chip_report_writes() {
    let d = Dir::new("chip-og");
    let genes = d.put(
        "genes.tsv",
        "gene\ttruncating\tmissense\tnote\tpmids\nDNMT3A\tall\tR882\tn\t25426837\nJAK2\tnone\tV617F\tn\t15858187\n",
    );
    let gff = d.put(
        "o.gff3",
        "##gff-version 3\n\
chr2\tensembl\tgene\t1000\t5000\t.\t-\t.\tID=gene:G1;Name=DNMT3A\n\
chr2\tensembl\tmRNA\t1000\t5000\t.\t-\t.\tID=transcript:T2;Parent=gene:G1;transcript_id=T2\n\
chr2\tensembl\texon\t1000\t1300\t.\t-\t.\tParent=transcript:T2;rank=2\n\
chr2\tensembl\texon\t4000\t5000\t.\t-\t.\tParent=transcript:T2;rank=1\n",
    );
    let bed = d.put(
        "o.bed",
        "chr2\t1039\t1310\tDNMT3A\nchr2\t3989\t4110\tDNMT3A\n",
    );
    let calls = d.put(
        "calls.tsv",
        "chr2\t1150\tC\tT\tPASS\tmissense|DNMT3A|T2|protein_coding|-|882R>882H|1150C>T\t40,6\t0.130\t46\n",
    );
    let depth: String = (1040..1311)
        .chain(3990..4111)
        .map(|p| format!("chr2\t{p}\t30\n"))
        .collect();
    let depth = d.put("depth.tsv", &depth);
    fs::create_dir_all(d.p("S/chip")).unwrap();
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
        &d.p("S/chip/S_chip.tsv"),
        "--depth-out",
        &d.p("S/chip/S_chip_depth.tsv"),
    ]
    .iter()
    .map(|s| (*s).to_owned())
    .collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    assert_eq!(
        og_chip::cli("report", &args, &mut o, &mut e),
        0,
        "{}",
        String::from_utf8_lossy(&e)
    );
    let s = collect(&d.p("S"), "S");
    let c = data(&s, "chip");
    assert_eq!(column(at(c, &["calls"]), "protein_change"), ["R882H"]);
    assert_eq!(texts(at(c, &["not_covered"])), ["JAK2"]);
}

// --- step 47, viral -----------------------------------------------------------------

const VIRAL_HEAD: &str = "virus\taccession\tlength\treads\tbreadth\tmean_depth\tmode_depth\thuman_depth\trelative_depth\tcopies_per_cell\tcall\tnote\n";

/// One row of the viral table.
fn viral_row(virus: &str, acc: &str, reads: &str, human: &str, call: &str) -> String {
    let rel = if human == "." { "." } else { "0.500" };
    format!("{virus}\t{acc}\t2000\t{reads}\t1.0000\t15.000\t15\t{human}\t{rel}\t.\t{call}\t.\n")
}

fn viral_table(rows: &[String]) -> String {
    format!("{VIRAL_HEAD}{}", rows.concat())
}

#[test]
fn viral_reads_the_table() {
    let d = Dir::new("viral");
    d.put(
        "S/viral/S_viral.tsv",
        &viral_table(&[
            viral_row("HHV-6A", "NC_001664.4", "0", "30.00", "not found"),
            viral_row(
                "HHV-6B",
                "NC_000898.1",
                "400",
                "30.00",
                "inherited integrated HHV-6 likely",
            ),
            viral_row(
                "Human gammaherpesvirus 4",
                "NC_007605.1",
                "30",
                "30.00",
                "viral DNA found",
            ),
            viral_row(
                "Escherichia phage phiX174",
                "NC_001422.1",
                "20",
                "30.00",
                "lab control (phiX)",
            ),
        ]),
    );
    let s = collect(&d.p("S"), "S");
    assert!(text(at(&s, &["sections", "viral", "source"])).ends_with("viral/S_viral.tsv"));
    let v = data(&s, "viral");
    assert_eq!(text(at(v, &["hhv6"])), "inherited integrated HHV-6B likely");
    assert_eq!(
        column(at(v, &["found"]), "virus"),
        [
            "HHV-6B",
            "Human gammaherpesvirus 4",
            "Escherichia phage phiX174"
        ]
    );
    fs::create_dir_all(d.p("none")).unwrap();
    assert!(
        missing(&collect(&d.p("none"), "S"), "viral"),
        "no table: no section"
    );
    d.put("S/viral/S_viral.tsv", "virus\tcall\n");
    assert_eq!(state(&collect(&d.p("S"), "S"), "viral"), "unreadable");
}

#[test]
fn viral_hhv6_line_without_a_human_depth_and_without_integration() {
    let d = Dir::new("viral-hhv6");
    for (human, want) in [
        (".", "not assessed (no human depth)"),
        ("30.00", "no integrated HHV-6"),
    ] {
        d.put(
            "S/viral/S_viral.tsv",
            &viral_table(&[
                viral_row("HHV-6B", "NC_000898.1", "0", human, "not found"),
                viral_row(
                    "Human gammaherpesvirus 4",
                    "NC_007605.1",
                    "30",
                    human,
                    "viral DNA found",
                ),
            ]),
        );
        let s = collect(&d.p("S"), "S");
        assert_eq!(text(at(data(&s, "viral"), &["hhv6"])), want);
    }
}

#[test]
fn viral_lists_a_genome_under_ten_reads_as_below_threshold() {
    let d = Dir::new("viral-below");
    d.put(
        "S/viral/S_viral.tsv",
        &viral_table(&[
            viral_row("HHV-6B", "NC_000898.1", "0", "30.00", "not found"),
            viral_row(
                "Hepatitis C virus genotype 1",
                "NC_004102.1",
                "5",
                "30.00",
                "below threshold",
            ),
        ]),
    );
    let s = collect(&d.p("S"), "S");
    let found = at(data(&s, "viral"), &["found"]);
    assert_eq!(column(found, "virus"), ["Hepatitis C virus genotype 1"]);
    assert_eq!(column(found, "call"), ["below threshold"]);
    assert_eq!(column(found, "reads"), ["5"]);
}

// --- step 50, MEI -------------------------------------------------------------------

const MEI_NOT_ASSESSED: &str = "Mobile element insertions (Alu, L1 and SVA, the FKTN, TAF1 and BBS1 SVA insertions among them): not looked for (step 50, MEGAnE, is opt-in)";
const MEI_CALLS: &str = "chrom\tpos\tend\telement\tsubfamily\tcopies\tfilter\tconfidence\tgenes\tlocation\tbenign_af\tacmg_class\n\
chr9\t105639660\t105639668\tSVA\tSVA_F\t1\tPASS\tH\t.\t.\t.\t.\n\
chr20\t5040478\t5040486\tAlu\tSVA_F\t1\tPASS\tH\t.\t.\t.\t.\n";
const MEI_SITES: &str = "gene\tchrom\tstart\tend\tsite\telement\tinheritance\tcondition\tclinvar\tpmid\tstatus\tmean_depth\tcall\n\
FKTN\tchr9\t105639156\t105640156\t105639656\tSVA\tautosomal recessive\tFukuyama congenital muscular dystrophy\t3199\t9690476\tseen\t30.0\tchr9:105639660 SVA_F 1 copy PASS\n\
TAF1\tchrX\t71439993\t71440993\t71440493\tSVA\tX-linked\tX-linked dystonia-parkinsonism\t9862\t29474918\tnot seen\t30.0\t.\n";

#[test]
fn mei_section_and_its_not_assessed_line() {
    let d = Dir::new("mei");
    fs::create_dir_all(d.p("S")).unwrap();
    let s = collect(&d.p("S"), "S");
    assert!(not_assessed(&s).iter().any(|x| x == MEI_NOT_ASSESSED));
    assert!(missing(&s, "mei"));
    // A step 50 that ran and failed was looked for: not run, not opt-in.
    d.put(
        "F/logs/run_status.tsv",
        &format!(
            "meta\tstarted_epoch\t{}\nstep\t50\tfailed (exit 1)\n",
            now() - 60.0
        ),
    );
    let s = collect(&d.p("F"), "F");
    assert_eq!(state(&s, "mei"), "failed");
    assert!(!not_assessed(&s).iter().any(|x| x == MEI_NOT_ASSESSED));
    assert!(
        texts(at(&s, &["not_run"]))
            .iter()
            .any(|x| x.contains("step 50 failed"))
    );

    let calls = d.put("S/mei/S_mei.tsv", MEI_CALLS);
    d.put("S/mei/S_mei_sites.tsv", MEI_SITES);
    let s = collect(&d.p("S"), "S");
    assert!(text(at(&s, &["sections", "mei", "source"])).ends_with("mei/S_mei_sites.tsv"));
    let m = data(&s, "mei");
    assert_eq!(text(at(m, &["calls"])), "2");
    assert_eq!(
        (
            text(at(m, &["by_element", "Alu"])),
            text(at(m, &["by_element", "SVA"]))
        ),
        ("1".into(), "1".into())
    );
    let sites: Vec<(String, String)> = list(at(m, &["sites"]))
        .iter()
        .map(|x| (text(at(x, &["gene"])), text(at(x, &["status"]))))
        .collect();
    assert_eq!(
        sites,
        [
            ("FKTN".into(), "seen".into()),
            ("TAF1".into(), "not seen".into())
        ]
    );
    assert!(!not_assessed(&s).iter().any(|x| x == MEI_NOT_ASSESSED));

    fs::write(&calls, format!("{}\n", MEI_CALLS.lines().next().unwrap())).unwrap();
    assert_eq!(
        text(at(data(&collect(&d.p("S"), "S"), "mei"), &["calls"])),
        "0"
    );
    for bad in ["", "chrom\tpos\n"] {
        fs::write(&calls, bad).unwrap();
        assert_eq!(
            state(&collect(&d.p("S"), "S"), "mei"),
            "unreadable",
            "{bad:?}"
        );
    }
    fs::write(&calls, MEI_CALLS).unwrap();
    d.put("S/mei/S_mei_sites.tsv", "gene\tstatus\nFKTN\tseen\n");
    assert_eq!(state(&collect(&d.p("S"), "S"), "mei"), "unreadable");
}

// --- step 51, MoChA -----------------------------------------------------------------

const MCA_NOT_ASSESSED: &str = "Autosomal mosaic chromosomal alterations (large gains, losses and copy-neutral LOH in a share of the blood cells): step 51 (MoChA) did not assess this sample";
const MOCHA_CALLS: &str = "chrom\tstart\tend\tlength_mb\tarms\ttype\tcell_fraction\tbdev\trel_cov\tn_hets\tlod_baf_phase\tlod_lrr_baf\treported\treason\n";

fn mocha_summary(status: &str, message: &str, reported: usize, all: usize) -> String {
    format!(
        "metric\tvalue\nstatus\t{status}\nmessage\t{message}\ndepth\t30.0\nhet_sites\t18000\nphased_het_sites\t18000\n\
het_per_10mb\t3000\nanalysed_mb\t60.0\nchromosomes\tchr20\ncall_rate\t0.9990\nbaf_auto\t0.0010\nfloor_read\t0.020\n\
floor_min\t0.05\nfloor\t0.050\ncalls_mocha\t{all}\ncalls_reported\t{reported}\nphasing\tinput\n"
    )
}

const CNLOH: &str = "chr20\t60000\t20000000\t19.94\t20p, to the p telomere\tCN-LOH\t0.150\t0.075\t2.0\t6600\t900.0\tnan\tyes\t.\n";
const CHRX: &str = "chrX\t60000\t20000000\t19.94\tXp, to the p telomere\tCN-LOH\t0.150\t0.075\t2.0\t6600\t900.0\tnan\tno\tnot an autosome (sex chromosomes are step 33's)\n";

#[test]
fn mosaic_lists_the_reported_call_and_the_not_assessed_line() {
    let d = Dir::new("mocha");
    let sum = d.put(
        "S/mosaic/S_mocha_summary.tsv",
        &mocha_summary("assessed", ".", 1, 2),
    );
    let calls = d.put(
        "S/mosaic/S_mocha.tsv",
        &format!("{MOCHA_CALLS}{CNLOH}{CHRX}"),
    );
    let s = collect(&d.p("S"), "S");
    let m = data(&s, "mosaic");
    let got: Vec<[String; 3]> = list(at(m, &["calls"]))
        .iter()
        .map(|x| ["type", "chrom", "cell_fraction"].map(|k| text(at(x, &[k]))))
        .collect();
    assert_eq!(got, [["CN-LOH", "chr20", "0.150"].map(str::to_owned)]);
    assert_eq!(text(at(m, &["not_reported"])), "1");
    assert!(!not_assessed(&s).iter().any(|x| x == MCA_NOT_ASSESSED));
    // Large mosaic alterations stay out of the chrM-only line.
    let mosaic: Vec<String> = not_assessed(&s)
        .into_iter()
        .filter(|x| x.starts_with("Mosaic"))
        .collect();
    assert_eq!(
        mosaic,
        [
            "Mosaic and low-fraction small variants (only chrM heteroplasmy is reported, and with step 49 the CHIP genes)"
        ]
    );
    fs::remove_file(&sum).unwrap();
    let s = collect(&d.p("S"), "S");
    assert!(
        not_assessed(&s).iter().any(|x| x == MCA_NOT_ASSESSED),
        "without step 51"
    );

    d.put(
        "S/mosaic/S_mocha_summary.tsv",
        &mocha_summary(
            "not assessed",
            "phased BAF autocorrelation 0.0500 over 0.03",
            0,
            1,
        ),
    );
    fs::write(
        &calls,
        format!(
            "{MOCHA_CALLS}{}",
            CNLOH.replace("\tyes\t.", "\tno\tsample not assessed")
        ),
    )
    .unwrap();
    let s = collect(&d.p("S"), "S");
    assert_eq!(text(at(data(&s, "mosaic"), &["status"])), "not assessed");
    assert!(
        not_assessed(&s).iter().any(|x| x == MCA_NOT_ASSESSED),
        "not assessed keeps the line"
    );
}

#[test]
fn mosaic_is_unreadable_without_a_call_table_that_agrees() {
    let d = Dir::new("mocha-bad");
    d.put(
        "S/mosaic/S_mocha_summary.tsv",
        &mocha_summary("assessed", ".", 1, 1),
    );
    let s = collect(&d.p("S"), "S");
    assert_eq!(state(&s, "mosaic"), "unreadable", "no call table");
    let calls = d.put("S/mosaic/S_mocha.tsv", MOCHA_CALLS);
    assert_eq!(
        state(&collect(&d.p("S"), "S"), "mosaic"),
        "unreadable",
        "fewer reported calls"
    );
    fs::write(&calls, "chrom\tstart\n").unwrap();
    assert_eq!(
        state(&collect(&d.p("S"), "S"), "mosaic"),
        "unreadable",
        "not MoChA's columns"
    );
    fs::write(&calls, format!("{MOCHA_CALLS}{CNLOH}")).unwrap();
    assert_eq!(state(&collect(&d.p("S"), "S"), "mosaic"), "ok", "control");
}

// --- step 52, VNtyper ---------------------------------------------------------------

const MUC1: &str = "MUC1 VNTR frameshifts (ADTKD-MUC1): step 52 (VNtyper, opt-in) gave no reading; no standard caller sees them";

fn vntyper_table(status: &str) -> String {
    format!(
        "metric\tvalue\nstatus\t{status}\nconfidence\tHigh_Precision*\nvariant\t59dupC\nflag\tNot flagged\n\
alt_depth\t416\nregion_depth\t7110\ndepth_score\t0.0585\nvntr_mean_depth\t1075.32\nvntr_uncovered_pct\t0.00\n\
coverage_qc\tPASS\nsensitivity\tAt 1075x the 100x-or-more band: 160 of 161 runs (99%).\nmessage\tm\nadvice\ta\n"
    )
}

#[test]
fn vntyper_section_and_the_muc1_line() {
    let d = Dir::new("vntyper");
    for (status, listed) in [
        ("positive", false),
        ("negative", false),
        ("not assessed", true),
    ] {
        d.put("S/vntyper/S_vntyper.tsv", &vntyper_table(status));
        let s = collect(&d.p("S"), "S");
        assert_eq!(state(&s, "vntyper"), "ok");
        assert_eq!(text(at(data(&s, "vntyper"), &["status"])), status);
        assert_eq!(
            not_assessed(&s).iter().any(|x| x == MUC1),
            listed,
            "{status}"
        );
    }
    d.put("S/vntyper/S_vntyper.tsv", &vntyper_table("positive"));
    d.put(
        "S/logs/run_status.tsv",
        &format!("meta\tstarted_epoch\t{}\nstep\t11\tok\n", now() + 60.0),
    );
    let s = collect(&d.p("S"), "S");
    assert_eq!(
        state(&s, "vntyper"),
        "stale",
        "older than a run that did not run step 52"
    );
    assert!(not_assessed(&s).iter().any(|x| x == MUC1));
    fs::create_dir_all(d.p("none")).unwrap();
    let s = collect(&d.p("none"), "S");
    assert!(not_assessed(&s).iter().any(|x| x == MUC1) && missing(&s, "vntyper"));
}

// --- step 53, Locityper -------------------------------------------------------------

const LOCITYPER: &str = "locus\tregion\tstatus\tgenotype\thaplotype_1\thaplotype_2\tquality\ttotal_reads\tunexplained_reads\tweighted_distance\twarnings\tfilter\treason\n\
HLA-DQB1\tchr6:32652604-32669077\tcalled\tHG00733.1,NA19240.2\tHG00733.1\tNA19240.2\t27.8\t900\t12\t50.000\t.\tfail\tweighted distance 50.00 over 30\n\
CYP2D6\tchr22:42122285-42152960\tcalled\tHG00733.1,NA19240.2\tHG00733.1\tNA19240.2\t27.8\t900\t12\t0.200\t.\tpass\t.\n\
FCGR\tchr1:161506762-161681047\tnot called\t.\t.\t.\t.\t.\t.\t.\t.\t.\tLocityper wrote no result for this locus\n";

fn cyp2d6(s: &Value) -> String {
    list(at(data(s, "locityper"), &["loci"]))
        .iter()
        .find(|x| text(at(x, &["locus"])) == "CYP2D6")
        .map(|x| text(at(x, &["primary"])))
        .unwrap()
}

#[test]
fn locityper_counts_and_the_primary_calls() {
    let d = Dir::new("locityper");
    d.put("S/locityper/S_locityper.tsv", LOCITYPER);
    let s = collect(&d.p("S"), "S");
    let l = data(&s, "locityper");
    assert_eq!(state(&s, "locityper"), "ok");
    assert_eq!(
        (text(at(l, &["called"])), text(at(l, &["passed"]))),
        ("2".into(), "1".into())
    );
    let primary: Vec<String> = column(at(l, &["loci"]), "primary")
        .into_iter()
        .map(|p| p.chars().take(5).collect())
        .collect();
    assert_eq!(primary, ["T1K (", "pypgx", "no ot"]);
    d.put(
        "S/cyrius/S_cyp2d6.tsv",
        "Sample\tGenotype\tFilter\nS\tfixture_call\tLowQ_high_CN\n",
    );
    let s = collect(&d.p("S"), "S");
    assert!(
        cyp2d6(&s).contains("Cyrius (step 21): fixture_call (not usable: Filter LowQ_high_CN)"),
        "{}",
        cyp2d6(&s)
    );
    d.put(
        "S/logs/run_status.tsv",
        &format!("meta\tstarted_epoch\t{}\nstep\t53\tok\n", now() + 60.0),
    );
    let s = collect(&d.p("S"), "S");
    assert!(
        cyp2d6(&s).contains("(stale: from an earlier run)"),
        "{}",
        cyp2d6(&s)
    );
}

// --- step 48, Exomiser --------------------------------------------------------------

const EXOMISER_TOP: &str = "rank\tgene\tmoi\texomiser_score\tphenotype_score\tvariant_score\tp_value\tvariant\tgenotype\tconsequence\thgvs\tacmg\tacmg_evidence\tdisease\tclinvar\tmax_frequency\n";

fn top_row(rank: &str, gene: &str, variant: &str) -> String {
    format!(
        "{rank}\t{gene}\tAD\t0.998\t0.800\t0.950\t0.0010\t{variant}\t0/1\tmissense_variant\tp.(Glu563Ala)\tPATHOGENIC\tPM2\t.\t.\tnot seen\n"
    )
}

#[test]
fn exomiser_section() {
    let d = Dir::new("exomiser");
    let top = d.put(
        "S1/exomiser/S1_exomiser_top.tsv",
        &format!(
            "{EXOMISER_TOP}{}{}{}{}",
            top_row("1", "GENEA", "chr10:123256215 T>G"),
            top_row("2", "GENEB", "chr6:1000 G>A"),
            top_row("2", "GENEB", "chr6:2000 C>T"),
            top_row("3", "GENEC", "chrX:5000 A>G")
        ),
    );
    d.put(
        "S1/exomiser/S1_exomiser_terms.tsv",
        "hpo_id\tlabel\tgiven_as\nHP:0001156\tBrachydactyly\tHP:0001156\n",
    );
    let s = collect(&d.p("S1"), "S1");
    assert!(!matches!(
        at(&s, &["sections", "exomiser", "source"]),
        Value::Null
    ));
    let e = data(&s, "exomiser");
    let genes: Vec<(String, String, usize)> = list(at(e, &["genes"]))
        .iter()
        .map(|g| {
            (
                text(at(g, &["rank"])),
                text(at(g, &["gene"])),
                list(at(g, &["variants"])).len(),
            )
        })
        .collect();
    assert_eq!(
        genes,
        [
            ("1".into(), "GENEA".into(), 1),
            ("2".into(), "GENEB".into(), 2),
            ("3".into(), "GENEC".into(), 1)
        ]
    );
    assert_eq!(column(at(e, &["terms"]), "id"), ["HP:0001156"]);
    assert_eq!(column(at(e, &["terms"]), "label"), ["Brachydactyly"]);
    assert_eq!(text(at(&s, &["sections", "exomiser", "step"])), "48");
    fs::write(&top, "rank\tgene\n1\tGENEA\n").unwrap();
    assert_eq!(state(&collect(&d.p("S1"), "S1"), "exomiser"), "unreadable");
    // The module's -stub table: the header alone reads as no gene.
    fs::write(&top, EXOMISER_TOP).unwrap();
    let s = collect(&d.p("S1"), "S1");
    assert!(list(at(data(&s, "exomiser"), &["genes"])).is_empty());
}
