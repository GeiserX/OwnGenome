//! `og mei report` and `og mei annotsv-input` (step 50): the Alu, L1 and
//! SVA insertions MEGAnE calls, and the known disease insertions. Moved
//! from tests/test_mei_report.py (its cases 1 to 4, 6 and 7); the summary's
//! reading is in og-summary's `tests/optin_summary.rs` and the report lines
//! in og-render's `tests/optin_report.rs`.

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

const CATALOG: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets/mei_sites.tsv");

const HEAD: &str = "##fileformat=VCFv4.1
##source=MEI search version \"v1.0.1 2022/04/12\"
##ALT=<ID=INS:ME,Description=\"Insertion of mobile element\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=MEPRED,Number=1,Type=String,Description=\"ME prediction status\">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
";

const SITES: &str = "# a test catalog
gene\tchrom\tstart\tend\tsite\telement\tinheritance\tcondition\tclinvar\tpmid
FKTN\tchr9\t105639156\t105640156\t105639656\tSVA\tautosomal recessive\tFukuyama congenital muscular dystrophy\t3199\t9690476
TAF1\tchrX\t71439993\t71440993\t71440493\tSVA\tX-linked\tX-linked dystonia-parkinsonism\t9862\t29474918
";

/// chrom, bedcov start, end, gene, summed depth over the 1001 bases, reads.
const DEPTH_30X: &str =
    "chr9\t105639155\t105640156\tFKTN\t30030\t200\nchrX\t71439992\t71440993\tTAF1\t30030\t200\n";

/// One record as MEGAnE writes it.
struct Rec<'a> {
    chrom: &'a str,
    pos: u64,
    svtype: &'a str,
    mei: &'a str,
    filt: &'a str,
    mepred: &'a str,
    cn: &'a str,
    conf: &'a str,
}

fn rec<'a>(chrom: &'a str, pos: u64, svtype: &'a str) -> Rec<'a> {
    Rec {
        chrom,
        pos,
        svtype,
        mei: "MEI=SVA_F,1/1500,+/+",
        filt: "PASS",
        mepred: "PASS",
        cn: "1",
        conf: "H",
    }
}

impl Rec<'_> {
    fn line(&self) -> String {
        format!(
            "{}\t{}\t.\tN\t<INS:ME>\t.\t{}\tSVTYPE={};MEPRED={};HOMLEN=9;MEI_rpos={};{};CN_conf={}\tCN\t{}\n",
            self.chrom,
            self.pos,
            self.filt,
            self.svtype,
            self.mepred,
            self.pos + 8,
            self.mei,
            self.conf,
            self.cn
        )
    }
}

fn fktn_sva<'a>() -> Rec<'a> {
    rec("chr9", 105_639_660, "Retroposon/SVA")
}

fn vcf(recs: &[Rec<'_>]) -> String {
    format!("{HEAD}{}", recs.iter().map(Rec::line).collect::<String>())
}

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!("og-optin-mei-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn p(&self, n: &str) -> String {
        self.0.join(n).to_str().unwrap().to_owned()
    }
    fn put(&self, n: &str, t: &str) -> String {
        fs::write(self.p(n), t).unwrap();
        self.p(n)
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn og(verb: &str, args: &[&str]) -> (u8, String) {
    let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
    let (mut o, mut e) = (Vec::new(), Vec::new());
    let code = og_reports::mei::cli(Some(verb), &args, &mut o, &mut e);
    (code, String::from_utf8(e).unwrap())
}

type Row = HashMap<String, String>;

/// The calls, the sites, and the header of the calls table.
type Tables = (Vec<Row>, Vec<Row>, Vec<String>);

fn rows(text: &str) -> (Vec<String>, Vec<Row>) {
    let mut lines = text.lines();
    let head: Vec<String> = lines
        .next()
        .unwrap_or("")
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

/// (calls, sites) `og mei report` writes, or its stderr.
fn run_with(
    tag: &str,
    vcf_text: &str,
    depth: &str,
    sites: &str,
    annotsv: Option<&str>,
) -> Result<Tables, String> {
    let d = Dir::new(tag);
    let (v, s, dp) = (
        d.put("m.vcf", vcf_text),
        d.put("sites.tsv", sites),
        d.put("depth.tsv", depth),
    );
    let (c, so) = (d.p("c.tsv"), d.p("s.tsv"));
    let mut args = vec![
        "--sample",
        "S",
        "--vcf",
        &v,
        "--sites",
        &s,
        "--depth",
        &dp,
        "--out-calls",
        &c,
        "--out-sites",
        &so,
    ];
    let a;
    if let Some(t) = annotsv {
        a = d.put("a.tsv", t);
        args.extend(["--annotsv", &a]);
    }
    let (code, err) = og("report", &args);
    if code != 0 {
        return Err(err);
    }
    let (head, calls) = rows(&fs::read_to_string(&c).unwrap());
    let (_, sites) = rows(&fs::read_to_string(&so).unwrap());
    Ok((calls, sites, head))
}

fn run(tag: &str, recs: &[Rec<'_>], depth: &str) -> (Vec<Row>, Vec<Row>) {
    let (c, s, _) = run_with(tag, &vcf(recs), depth, SITES, None).unwrap();
    (c, s)
}

fn site<'a>(r: &'a (Vec<Row>, Vec<Row>), gene: &str) -> &'a Row {
    r.1.iter().find(|x| x["gene"] == gene).unwrap()
}

#[test]
fn the_planted_fktn_sva() {
    let r = run(
        "planted",
        &[rec("chr20", 5_040_478, "SINE/Alu"), fktn_sva()],
        DEPTH_30X,
    );
    let s = site(&r, "FKTN");
    assert_eq!(
        [s["status"].as_str(), &s["call"], &s["mean_depth"]],
        ["seen", "chr9:105639660 SVA_F 1 copy PASS", "30.0"]
    );
    let status = |tag, recs: &[Rec<'_>]| site(&run(tag, recs, DEPTH_30X), "FKTN")["status"].clone();
    assert_eq!(
        status("without", &[rec("chr20", 5_040_478, "SINE/Alu")]),
        "not seen"
    );
    assert_eq!(
        status("past", &[rec("chr9", 105_641_200, "Retroposon/SVA")]),
        "not seen",
        "1 kb past the window"
    );
    assert_eq!(
        status(
            "alu",
            &[Rec {
                mei: "MEI=AluYa5,0/288,+/+",
                ..rec("chr9", 105_639_660, "SINE/Alu")
            }]
        ),
        "not seen",
        "an Alu is not the catalog's SVA"
    );
    let two = run(
        "two",
        &[Rec {
            cn: "2",
            ..fktn_sva()
        }],
        DEPTH_30X,
    );
    assert_eq!(
        site(&two, "FKTN")["call"],
        "chr9:105639660 SVA_F 2 copies PASS"
    );
    assert_eq!(
        status(
            "failed",
            &[Rec {
                mepred: "FAILED",
                ..fktn_sva()
            }]
        ),
        "not seen"
    );
}

#[test]
fn too_few_reads() {
    let thin = "chr9\t105639155\t105640156\tFKTN\t9009\t60\n";
    let r = run("thin", &[], thin);
    assert_eq!(
        [
            site(&r, "FKTN")["status"].as_str(),
            &site(&r, "FKTN")["mean_depth"],
            &site(&r, "TAF1")["status"],
            &site(&r, "TAF1")["mean_depth"]
        ],
        ["not callable", "9.0", "not callable", "."]
    );
    let r = run("thin-call", &[fktn_sva()], thin);
    assert_eq!(site(&r, "FKTN")["status"], "seen", "a call at low depth");
    let r = run("nodepth", &[fktn_sva()], "");
    let all: Vec<&str> = r.1.iter().map(|x| x["status"].as_str()).collect();
    assert_eq!(all, ["not callable", "not callable"]);
}

#[test]
fn the_calls_table() {
    let recs = [
        Rec {
            mei: "MEI=AluYb8,0/298,+/+|AluYb9,0/298,+/+",
            filt: "LC",
            conf: "L",
            ..rec("chr20", 5_040_478, "SINE/Alu")
        },
        Rec {
            mei: "MEI_lbp=L1HS,5581,+;MEI_rbp=pA",
            cn: "2",
            ..rec("chr2", 1000, "LINE/L1")
        },
        Rec {
            mei: "MEI_lbp=NA,NA,+;MEI_rbp=pA",
            ..rec("chr2", 2000, "LINE/L1")
        },
        Rec {
            mei: "MEI=SVA_E,1/900,+/+",
            ..rec("chrX", 3000, "Retroposon/SVA")
        },
        Rec {
            mei: "MEI=Charlie5,2546/2611,+/+",
            ..rec("chr20", 6_665_193, "DNA/hAT-Charlie")
        },
        Rec {
            mei: "MEI=MER41A,91/60,-/-",
            ..rec("chr20", 16_134_625, "LTR/ERV1")
        },
        Rec {
            mepred: "FAILED",
            ..rec("chr20", 28_576_054, "LINE/L1")
        },
        Rec {
            cn: "0",
            ..rec("chr20", 30_000_000, "SINE/Alu")
        },
    ];
    let (c, _, head) = run_with("calls", &vcf(&recs), DEPTH_30X, SITES, None).unwrap();
    let kept: Vec<(&str, &str, &str)> = c
        .iter()
        .map(|x| {
            (
                x["chrom"].as_str(),
                x["pos"].as_str(),
                x["element"].as_str(),
            )
        })
        .collect();
    assert_eq!(
        kept,
        [
            ("chr2", "1000", "L1"),
            ("chr2", "2000", "L1"),
            ("chr20", "5040478", "Alu"),
            ("chrX", "3000", "SVA")
        ]
    );
    assert_eq!(c[2]["subfamily"], "AluYb8", "the first MEI name");
    assert_eq!(c[0]["subfamily"], "L1HS", "from a breakpoint");
    assert_eq!(c[1]["subfamily"], ".");
    assert_eq!(
        ["copies", "filter", "confidence", "end"].map(|k| c[2][k].as_str()),
        ["1", "LC", "L", "5040486"]
    );
    assert_eq!(
        head,
        [
            "chrom",
            "pos",
            "end",
            "element",
            "subfamily",
            "copies",
            "filter",
            "confidence",
            "genes",
            "location",
            "benign_af",
            "acmg_class"
        ]
    );
    assert!(
        ["genes", "location", "benign_af", "acmg_class"]
            .iter()
            .all(|k| c[0][*k] == "."),
        "without --annotsv"
    );
    let r = run("empty", &[], DEPTH_30X);
    let all: Vec<&str> = r.1.iter().map(|x| x["status"].as_str()).collect();
    assert!(r.0.is_empty() && all == ["not seen", "not seen"]);
}

#[test]
fn damaged_input_stops_og() {
    let err = |tag, v: &str, depth, sites| run_with(tag, v, depth, sites, None).unwrap_err();
    assert!(err("nochrom", &fktn_sva().line(), DEPTH_30X, SITES).contains("not a VCF"));
    let other = format!("{HEAD}chr9\t1\t.\tN\t<INS:ME>\t.\tPASS\tEND=5\tCN\t1\n");
    assert!(err("other", &other, DEPTH_30X, SITES).contains("not MEGAnE's VCF"));
    assert!(
        err("bedcov", &vcf(&[]), "chr9\t105639155\t105640156\n", SITES).contains("samtools bedcov")
    );
    assert!(err("cat", &vcf(&[]), DEPTH_30X, "gene\tchrom\nFKTN\tchr9\n").contains("column"));
}

#[test]
fn the_real_catalog() {
    let text = fs::read_to_string(CATALOG).unwrap();
    let body: String = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect();
    let (_, rows) = rows(&body);
    let genes: Vec<&str> = rows.iter().map(|x| x["gene"].as_str()).collect();
    assert_eq!(genes, ["FKTN", "TAF1", "BBS1"]);
    assert!(rows.iter().all(|x| x["element"] == "SVA"));
    let n = |x: &Row, k: &str| x[k].parse::<i64>().unwrap();
    assert!(
        rows.iter()
            .all(|x| n(x, "start") == n(x, "site") - 500 && n(x, "end") == n(x, "site") + 500)
    );
    let (_, s, _) = run_with("real", &vcf(&[fktn_sva()]), DEPTH_30X, &text, None).unwrap();
    let all: Vec<&str> = s.iter().map(|x| x["status"].as_str()).collect();
    assert_eq!(all, ["seen", "not seen", "not callable"]);
}

const AHEAD: &str = "AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tID\tAnnotation_mode\tGene_name\tLocation\tB_ins_AFmax\tACMG_class\n";

fn table() -> String {
    format!(
        "{AHEAD}\
20_5040477_5040478_INS_1\t20\t5040477\t5040478\tINS\tchr20_5040478_Alu\tfull\tFIXG1;FIXG2\t\t0.0312\tNA\n\
20_5040477_5040478_INS_1\t20\t5040477\t5040478\tINS\tchr20_5040478_Alu\tsplit\tFIXG1\tintron2-intron2\t\tfull=NA\n\
20_5040477_5040478_INS_1\t20\t5040477\t5040478\tINS\tchr20_5040478_Alu\tsplit\tFIXG2\texon1-exon1\t\tfull=NA\n\
9_105639659_105639660_INS_1\t9\t105639659\t105639660\tINS\tchr9_105639660_SVA\tfull\t\t\t\tNA\n"
    )
}

fn recs<'a>() -> Vec<Rec<'a>> {
    vec![
        Rec {
            mei: "MEI=AluYa5,0/288,+/+",
            ..rec("chr20", 5_040_478, "SINE/Alu")
        },
        fktn_sva(),
        Rec {
            cn: "2",
            ..rec("chr20", 7_000_000, "LINE/L1")
        },
        Rec {
            mepred: "FAILED",
            ..rec("chr20", 8_000_000, "SINE/Alu")
        },
        rec("chr20", 9_000_000, "DNA/hAT-Charlie"),
    ]
}

/// (the AnnotSV input VCF, the calls joined with `table`), or the stderr.
fn annotsv_run(tag: &str, recs: &[Rec<'_>], table: &str) -> Result<(String, Vec<Row>), String> {
    let d = Dir::new(tag);
    let head = HEAD.replacen("##ALT", "##contig=<ID=chr20,length=64444167>\n##ALT", 1);
    let v = d.put(
        "m.vcf",
        &format!("{head}{}", recs.iter().map(Rec::line).collect::<String>()),
    );
    let input = d.p("in.vcf");
    let (code, err) = og(
        "annotsv-input",
        &["--sample", "S", "--vcf", &v, "--out", &input],
    );
    if code != 0 {
        return Err(err);
    }
    let input = fs::read_to_string(&input).unwrap();
    let (c, _, _) = run_with(
        &format!("{tag}-r"),
        &fs::read_to_string(&v).unwrap(),
        DEPTH_30X,
        SITES,
        Some(table),
    )?;
    Ok((input, c))
}

fn by_pos(c: &[Row]) -> HashMap<String, Row> {
    c.iter().map(|x| (x["pos"].clone(), x.clone())).collect()
}

const ANN: [&str; 4] = ["genes", "location", "benign_af", "acmg_class"];

#[test]
fn annotsv_input_and_join() {
    let (input, calls) = annotsv_run("ann", &recs(), &table()).unwrap();
    let body: Vec<Vec<&str>> = input
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    let ids: Vec<&str> = body.iter().map(|x| x[2]).collect();
    assert_eq!(
        ids,
        [
            "chr20_5040478_Alu",
            "chr9_105639660_SVA",
            "chr20_7000000_L1"
        ]
    );
    let recs_: Vec<[&str; 4]> = body.iter().map(|x| [x[4], x[6], x[7], x[9]]).collect();
    assert_eq!(
        recs_,
        [
            ["<INS:ME>", "PASS", "SVTYPE=INS;END=5040478", "0/1"],
            ["<INS:ME>", "PASS", "SVTYPE=INS;END=105639660", "0/1"],
            ["<INS:ME>", "PASS", "SVTYPE=INS;END=7000000", "1/1"]
        ]
    );
    assert!(input.contains("##contig=<ID=chr20,length=64444167>"));
    let c = by_pos(&calls);
    assert_eq!(
        ANN.map(|k| c["5040478"][k].as_str()),
        [
            "FIXG1,FIXG2",
            "FIXG1:intron2-intron2,FIXG2:exon1-exon1",
            "0.0312",
            "NA"
        ]
    );
    assert_eq!(
        ANN.map(|k| c["105639660"][k].as_str()),
        [".", ".", ".", "NA"]
    );
    assert_eq!(
        ANN.map(|k| c["7000000"][k].as_str()),
        ["."; 4],
        "not annotated"
    );
    // Red control: an ID that matches no call leaves the Alu unjoined.
    let (_, red) = annotsv_run(
        "red",
        &recs(),
        &table().replace("chr20_5040478_Alu", "chr20_5040479_Alu"),
    )
    .unwrap();
    assert_eq!(red[1]["genes"], ".");
    let (_, empty) = annotsv_run("emptyann", &recs(), "").unwrap();
    assert!(empty.iter().all(|x| x["acmg_class"] == "."));
    let (_, quoted) = annotsv_run(
        "quoted",
        &recs(),
        &table().replace("\tsplit\tFIXG1\t", "\tsplit\t\"FIXG1\t"),
    )
    .unwrap();
    let q = by_pos(&quoted);
    assert_eq!(
        (
            q["105639660"]["acmg_class"].as_str(),
            q["5040478"]["genes"].as_str()
        ),
        ("NA", "FIXG1,FIXG2"),
        "a field that starts with a quote is read as it is"
    );
    let err = annotsv_run("notann", &recs(), "SV_chrom\tSV_start\n20\t1\n").unwrap_err();
    assert!(err.contains("not AnnotSV's table"), "{err}");
    let (dup, _) = annotsv_run(
        "dup",
        &[rec("chr20", 100, "SINE/Alu"), rec("chr20", 100, "SINE/Alu")],
        AHEAD,
    )
    .unwrap();
    let ids: Vec<&str> = dup
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(2).unwrap())
        .collect();
    assert_eq!(ids, ["chr20_100_Alu", "chr20_100_Alu_2"]);
}

#[test]
fn the_report_options_are_required() {
    let d = Dir::new("required");
    let v = d.put("m.vcf", HEAD);
    // The usage message goes to the process's stderr (og_filters::args).
    let (code, _) = og("report", &["--sample", "S", "--vcf", &v]);
    assert_eq!(code, 2);
}
