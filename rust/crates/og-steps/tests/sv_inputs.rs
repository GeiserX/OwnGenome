//! `og sv inv-filter`, `og sv svlen-fix` and `og sv add-callers` on made-up
//! VCFs: what the old `tests/test_annotsv_sv_inputs.py` checked, by
//! meaning, in cargo test.
//!
//! 1. inv-filter, on Manta's file and a `convertInversion.py`-style
//!    converted file: a reciprocal pair and a pair that shares one Manta
//!    EVENT keep their INV records (INV_JUNCTIONS=2); a single 50 Mb INV3
//!    junction gives back its two breakends (INV_JUNCTIONS=1); a single
//!    100 kb junction stays an INV, a single 400 bp one also gets
//!    INVDUP_CANDIDATE; the deletion passes unchanged; the counts. Control:
//!    the converted file holds that 50 Mb INV.
//! 2. svlen-fix: REF 137 bp / ALT 2 bp gets SVLEN -135 and END POS+134 with
//!    the old values kept (control: the input says -136); a pure deletion
//!    keeps SVLEN and gets END POS+135; a net insertion gets its net SVLEN
//!    and keeps END; symbolic, breakend and multi-ALT records pass
//!    unchanged; a second run changes nothing; two sample columns stop it.
//! 3. add-callers: which Delly and CNVpytor calls are added beside Manta's,
//!    no-call records that must not hide another caller's call, and the
//!    control that Manta's file alone lacks the CNVpytor deletion.
//! 4. add-callers and Delly's read depth (FORMAT RDCN), with and without
//!    the sex and the PAR file; control: every planted call is PASS DEL or
//!    DUP, so without the depth check each would be added.

use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

use og_filters::TextLines;
use og_filters::inv::{self, Bounds};
use og_filters::svlen::{self, SvlenError};
use og_steps::sv::Sex;
use og_steps::sv::add_callers::{self, Opts};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// A scratch folder of its own for one test, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!("og-sv-inputs-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn put(&self, name: &str, text: &str) -> String {
        let p = self.0.join(name);
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn put_gz(&self, name: &str, text: &str) -> String {
        let p = self.0.join(name);
        let mut gz = flate2::write::GzEncoder::new(std::fs::File::create(&p).unwrap(), flate2::Compression::default());
        gz.write_all(text.as_bytes()).unwrap();
        gz.finish().unwrap();
        p.to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn records(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

/// INFO as a map; a flag maps to `None`.
fn info(rec: &[String]) -> HashMap<String, Option<String>> {
    rec[7]
        .split(';')
        .map(|item| match item.split_once('=') {
            Some((k, v)) => (k.to_owned(), Some(v.to_owned())),
            None => (item.to_owned(), None),
        })
        .collect()
}

fn info_of(rec: &[String], key: &str) -> Option<String> {
    info(rec).get(key).cloned().flatten()
}

fn by_id(recs: &[Vec<String>]) -> HashMap<String, Vec<String>> {
    recs.iter().map(|r| (r[2].clone(), r.clone())).collect()
}

// --- 1. og sv inv-filter -----------------------------------------------------------------

const HEAD: &str = "##fileformat=VCFv4.1\n\
##contig=<ID=chr1,length=248956422>\n\
##contig=<ID=chr20,length=64444167>\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";

fn bnd(chrom: &str, pos: u64, rid: &str, alt: &str, mate: &str) -> String {
    format!("{chrom}\t{pos}\t{rid}\tN\t{alt}\t500\tPASS\tSVTYPE=BND;MATEID={mate}\tGT\t0/1\n")
}

fn inv(chrom: &str, pos: u64, rid: &str, end: u64, side: u8, extra: &str) -> String {
    format!(
        "{chrom}\t{pos}\t{rid}\tN\t<INV>\t500\tPASS\tEND={end};SVTYPE=INV;SVLEN={}{extra};INV{side}\tGT\t0/1\n",
        end - pos
    )
}

const DEL: &str = "chr20\t40000000\tMantaDEL:9:0:0:0:0:0\tN\t<DEL>\t500\tPASS\tEND=40010000;SVTYPE=DEL;SVLEN=-10000\tGT\t0/1\n";

/// Manta's own file: each junction a pair of breakends.
fn raw() -> String {
    [
        // reciprocal inversion chr20:1,000,000-2,000,000: INV3 and INV5 junctions
        bnd("chr20", 1000000, "MantaBND:10:0:1:0:0:0:0", "N]chr20:2000000]", "MantaBND:10:0:1:0:0:0:1"),
        bnd("chr20", 1000101, "MantaBND:10:0:2:0:0:0:0", "[chr20:2000101[N", "MantaBND:10:0:2:0:0:0:1"),
        bnd("chr20", 2000000, "MantaBND:10:0:1:0:0:0:1", "N]chr20:1000000]", "MantaBND:10:0:1:0:0:0:0"),
        bnd("chr20", 2000101, "MantaBND:10:0:2:0:0:0:1", "[chr20:1000101[N", "MantaBND:10:0:2:0:0:0:0"),
        // one INV3 junction, 50 Mb
        bnd("chr20", 5000000, "MantaBND:20:0:1:0:0:0:0", "N]chr20:55000000]", "MantaBND:20:0:1:0:0:0:1"),
        bnd("chr20", 55000000, "MantaBND:20:0:1:0:0:0:1", "N]chr20:5000000]", "MantaBND:20:0:1:0:0:0:0"),
        // one INV5 junction, 400 bp
        bnd("chr20", 30000000, "MantaBND:30:0:1:0:0:0:0", "[chr20:30000400[N", "MantaBND:30:0:1:0:0:0:1"),
        bnd("chr20", 30000400, "MantaBND:30:0:1:0:0:0:1", "[chr20:30000000[N", "MantaBND:30:0:1:0:0:0:0"),
        // one INV3 junction, 100 kb
        bnd("chr1", 3000000, "MantaBND:40:0:1:0:0:0:0", "N]chr1:3100000]", "MantaBND:40:0:1:0:0:0:1"),
        bnd("chr1", 3100000, "MantaBND:40:0:1:0:0:0:1", "N]chr1:3000000]", "MantaBND:40:0:1:0:0:0:0"),
        DEL.to_owned(),
    ]
    .iter()
    .fold(HEAD.to_owned(), |a, l| a + l)
}

/// What `convertInversion.py` writes for it (an INV5 record moves one base
/// left).
fn converted() -> String {
    let ev = ";EVENT=MantaINV:50:0:1:0:0:0";
    [
        inv("chr1", 3000000, "MantaINV:40:0:1:0:0:0", 3100000, 3, ""),
        // two junctions 3 kb apart that Manta put in one EVENT
        inv("chr1", 10000000, "MantaINV:50:0:1:0:0:0", 10500000, 3, ev),
        inv("chr1", 10003000, "MantaINV:50:0:2:0:0:0", 10503000, 5, ev),
        inv("chr20", 1000000, "MantaINV:10:0:1:0:0:0", 2000000, 3, ""),
        inv("chr20", 1000100, "MantaINV:10:0:2:0:0:0", 2000100, 5, ""),
        inv("chr20", 5000000, "MantaINV:20:0:1:0:0:0", 55000000, 3, ""),
        inv("chr20", 29999999, "MantaINV:30:0:1:0:0:0", 30000399, 5, ""),
        DEL.to_owned(),
    ]
    .iter()
    .fold(HEAD.to_owned(), |a, l| a + l)
}

fn inv_filter(raw: &str, conv: &str, b: Bounds) -> (String, inv::Counts) {
    let mut out = Vec::new();
    let c = inv::run(raw, conv, b, &mut out).unwrap_or_else(|e| panic!("og sv inv-filter: {e}"));
    (String::from_utf8(out).unwrap(), c)
}

#[test]
fn inv_filter_keeps_paired_inversions_and_returns_long_single_junctions() {
    let w = Scratch::new("inv");
    let raw_path = w.put_gz("raw.vcf.gz", &raw());
    let conv = converted();
    let conv_path = w.put("conv.vcf", &conv);

    let long_inv = records(&conv)
        .into_iter()
        .filter(|r| {
            info_of(r, "SVTYPE").as_deref() == Some("INV")
                && info_of(r, "END").unwrap().parse::<u64>().unwrap() - r[1].parse::<u64>().unwrap() >= 1_000_000
                && r[2] == "MantaINV:20:0:1:0:0:0"
        })
        .count();
    assert_eq!(long_inv, 1, "control: the converted file holds the 50 Mb single-junction INV");

    let (out, counts) = inv_filter(&raw_path, &conv_path, Bounds::default());
    let recs = records(&out);
    let by = by_id(&recs);
    for rid in ["MantaINV:10:0:1:0:0:0", "MantaINV:10:0:2:0:0:0", "MantaINV:50:0:1:0:0:0", "MantaINV:50:0:2:0:0:0"] {
        assert_eq!(
            by.get(rid).and_then(|r| info_of(r, "INV_JUNCTIONS")).as_deref(),
            Some("2"),
            "both junctions present: {rid} keeps its INV, INV_JUNCTIONS=2"
        );
    }
    assert!(!by.contains_key("MantaINV:20:0:1:0:0:0"), "the 50 Mb single junction is no INV record any more");
    let mb_invs: Vec<&String> = recs
        .iter()
        .filter(|r| {
            info_of(r, "SVTYPE").as_deref() == Some("INV")
                && info_of(r, "END").unwrap().parse::<u64>().unwrap() - r[1].parse::<u64>().unwrap() > 1_000_000
        })
        .map(|r| &r[2])
        .collect();
    assert!(mb_invs.is_empty(), "no INV record longer than the largest paired inversion (1 Mb): {mb_invs:?}");
    for rid in ["MantaBND:20:0:1:0:0:0:0", "MantaBND:20:0:1:0:0:0:1"] {
        let r = by.get(rid).unwrap_or_else(|| panic!("its breakend {rid} is back"));
        assert_eq!(
            (info_of(r, "SVTYPE").as_deref(), info_of(r, "INV_JUNCTIONS").as_deref()),
            (Some("BND"), Some("1")),
            "its breakend {rid} is back, INV_JUNCTIONS=1"
        );
    }
    let r = &by["MantaINV:40:0:1:0:0:0"];
    assert!(
        info_of(r, "INV_JUNCTIONS").as_deref() == Some("1") && !info(r).contains_key("INVDUP_CANDIDATE"),
        "the 100 kb single junction stays an INV, INV_JUNCTIONS=1, no INVDUP_CANDIDATE: {r:?}"
    );
    let r = &by["MantaINV:30:0:1:0:0:0"];
    assert!(
        info_of(r, "INV_JUNCTIONS").as_deref() == Some("1") && info(r).get("INVDUP_CANDIDATE") == Some(&None),
        "the 400 bp single junction gets INVDUP_CANDIDATE: {r:?}"
    );
    let del: Vec<String> = DEL.trim_end().split('\t').map(str::to_owned).collect();
    assert!(recs.contains(&del), "the deletion passes unchanged");
    let order: Vec<(usize, u64)> = recs
        .iter()
        .map(|r| (usize::from(r[0] == "chr20"), r[1].parse().unwrap()))
        .collect();
    let mut sorted = order.clone();
    sorted.sort_unstable();
    assert_eq!(order, sorted, "records sorted by the ##contig order, then POS");
    assert!(
        out.contains("##INFO=<ID=INV_JUNCTIONS,") && out.contains("##INFO=<ID=INVDUP_CANDIDATE,"),
        "the header declares INV_JUNCTIONS and INVDUP_CANDIDATE"
    );
    assert!(
        out.find("##INFO=<ID=INV_JUNCTIONS,") < out.find("#CHROM"),
        "the header lines come before #CHROM"
    );
    // The summary line og prints is these four counts.
    assert_eq!(
        (counts.paired, counts.single, counts.short, counts.long),
        (4, 2, 1, 1),
        "4 INV records with both junctions, 2 single-junction INV records kept (1 shorter than 1000 bp \
         marked INVDUP_CANDIDATE), 1 single-junction inversion of 1000000 bp or more"
    );

    let (out, _) = inv_filter(
        &raw_path,
        &conv_path,
        Bounds {
            max_single_bp: 100_000_000,
            ..Bounds::default()
        },
    );
    assert!(
        records(&out)
            .iter()
            .any(|r| r[2] == "MantaINV:20:0:1:0:0:0" && info_of(r, "INV_JUNCTIONS").as_deref() == Some("1")),
        "with a 100 Mb limit the 50 Mb junction stays an INV (INV_JUNCTIONS=1)"
    );
}

// --- 2. og sv svlen-fix ------------------------------------------------------------------

const SV: &str = "##fileformat=VCFv4.1\n##contig=<ID=chr20,length=64444167>\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";

fn svlen_fix(path: &str) -> Result<(String, usize), SvlenError> {
    let mut out = Vec::new();
    let input = og_io::open(path).expect("open");
    let n = svlen::run(TextLines::file(input), &mut out)?;
    Ok((String::from_utf8(out).unwrap(), n))
}

#[test]
fn svlen_fix_sets_the_net_length_change() {
    let w = Scratch::new("svlen");
    let (ref137, alt2) = ("A".repeat(137), "AC".to_owned());
    let ins_alt = format!("A{}", "T".repeat(59));
    let rows = [
        ("1005400", "complex", ref137.as_str(), alt2.as_str(), "END=1005536;SVTYPE=DEL;SVLEN=-136"),
        ("1100000", "pure", ref137.as_str(), "A", "END=1100136;SVTYPE=DEL;SVLEN=-136"),
        ("1200000", "ins", "ACG", ins_alt.as_str(), "END=1200002;SVTYPE=INS;SVLEN=59"),
        ("1300000", "sym", "N", "<DEL>", "END=1310000;SVTYPE=DEL;SVLEN=-10000"),
        ("1400000", "bnd", "N", "N]chr20:1500000]", "SVTYPE=BND"),
        ("1600000", "multi", ref137.as_str(), "A,AC", "END=1600136;SVTYPE=DEL;SVLEN=-136"),
    ];
    let vcf = rows.iter().fold(SV.to_owned(), |a, (pos, id, r, alt, inf)| {
        a + &format!("chr20\t{pos}\t{id}\t{r}\t{alt}\t50\tPASS\t{inf}\tGT\t0/1\n")
    });
    let src = w.put("sv.vcf", &vcf);
    let input = by_id(&records(&vcf));
    assert_eq!(
        info_of(&input["complex"], "SVLEN").as_deref(),
        Some("-136"),
        "control: the input gives the REF 137 / ALT 2 record SVLEN -136"
    );

    let (out, n) = svlen_fix(&src).expect("og sv svlen-fix");
    let by = by_id(&records(&out));
    let c = info(&by["complex"]);
    let get = |m: &HashMap<String, Option<String>>, k: &str| m.get(k).cloned().flatten();
    assert_eq!(
        (get(&c, "SVLEN").as_deref(), get(&c, "END").as_deref()),
        (Some("-135"), Some("1005534")),
        "REF 137 / ALT 2: SVLEN -135, END POS+134 (an interval of 135 bp for AnnotSV)"
    );
    assert_eq!(
        (get(&c, "ORIG_SVLEN").as_deref(), get(&c, "ORIG_END").as_deref()),
        (Some("-136"), Some("1005536")),
        "REF 137 / ALT 2: the old values kept"
    );
    let p = info(&by["pure"]);
    assert!(
        get(&p, "SVLEN").as_deref() == Some("-136")
            && get(&p, "END").as_deref() == Some("1100135")
            && !p.contains_key("ORIG_SVLEN"),
        "pure deletion: SVLEN stays -136, END becomes POS+135: {p:?}"
    );
    let i = info(&by["ins"]);
    assert_eq!(
        (get(&i, "SVLEN").as_deref(), get(&i, "END").as_deref()),
        (Some("57"), Some("1200002")),
        "net insertion: SVLEN 57, END kept"
    );
    for rid in ["sym", "bnd", "multi"] {
        assert_eq!(by[rid], input[rid], "{rid} passes unchanged");
    }
    assert!(
        out.contains("##INFO=<ID=ORIG_SVLEN,") && out.contains("##INFO=<ID=ORIG_END,"),
        "the header declares ORIG_SVLEN and ORIG_END"
    );
    assert_eq!(n, 3, "the summary counts the 3 rewritten sequence-resolved records");
    let again = w.put("sv2.vcf", &out);
    let (out2, _) = svlen_fix(&again).expect("a second run");
    assert_eq!(out2, out, "a second run changes nothing");

    // One sample only: the reports read the GT and PS of the first sample
    // column, so a second person's column (here the carrier, after a 0/0)
    // must stop the step, not be read silently.
    let two = w.put(
        "two.vcf",
        &(SV.replace("\ts1\n", "\ts1\ts2\n")
            + "chr20\t1300000\tsym\tN\t<DEL>\t50\tPASS\tEND=1310000;SVTYPE=DEL\tGT\t0/0\t0/1\n"),
    );
    match svlen_fix(&two) {
        Err(e @ SvlenError::MultiSample(_)) => assert!(
            e.to_string().contains("holds 2 samples (s1, s2)"),
            "a VCF with two sample columns stops and names them: {e}"
        ),
        other => panic!("a VCF with two sample columns must stop (exit 2), not {other:?}"),
    }
}

// --- 3. og sv add-callers ----------------------------------------------------------------

const BASE: &str = "##fileformat=VCFv4.1\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr20,length=64444167>\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
chr20\t1000000\tMantaDEL:1\tN\t<DEL>\t500\tPASS\tEND=1100000;SVTYPE=DEL;SVLEN=-100000\tGT\t0/1\n\
chr20\t3000000\tMantaINV:2\tN\t<INV>\t500\tPASS\tEND=4000000;SVTYPE=INV;SVLEN=1000000;INV3\tGT\t0/1\n";

const DELLY: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
chr20\t1000500\tDEL00000001\tN\t<DEL>\t60\tPASS\tIMPRECISE;SVTYPE=DEL;END=1099000;CIPOS=-50,50;CIEND=-50,50\tGT:GQ\t0/1:40\n\
chr20\t2000000\tDEL00000002\tN\t<DEL>\t60\tPASS\tPRECISE;SVTYPE=DEL;END=2050000\tGT:GQ\t1/1:40\n\
chr20\t2500000\tDEL00000003\tN\t<DEL>\t10\tLowQual\tIMPRECISE;SVTYPE=DEL;END=2600000\tGT\t0/1\n\
chr20\t3000000\tINV00000004\tN\t<INV>\t60\tPASS\tIMPRECISE;SVTYPE=INV;END=4000000;CT=3to3\tGT\t0/1\n";

fn cnv(chrom: &str, pos: u64, rid: &str, alt: &str, end: u64, p1: &str, q0: &str, gt: &str) -> String {
    let t = &alt[1..4];
    format!(
        "chr{chrom}\t{pos}\t{rid}\tN\t{alt}\t.\tPASS\tEND={end};IMPRECISE;SVLEN={};SVTYPE={t};pytorRD=0.5;pytorP1={p1};\
pytorP2=0;pytorP3=0;pytorP4=0;pytorQ0={q0};pytorPN=0;pytorDG=1000;pytorCL=rd_mean_shift\tGT:CN\t{gt}:1\n",
        end - pos
    )
}

fn cnvpytor() -> String {
    "##fileformat=VCFv4.2\n##source=CNVpytor\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n".to_owned()
        + &cnv("20", 2001000, "CNVpytor_del0", "<DEL>", 2049000, "1e-06", "0.01", "1/1") // Delly's added call has it
        + &cnv("20", 3200000, "CNVpytor_del1", "<DEL>", 3700000, "1e-08", "0.02", "0/1") // inside Manta's INV: added
        + &cnv("20", 5000000, "CNVpytor_del2", "<DEL>", 5100000, "0.5", "0.0", "0/1") // e-val too high
        + &cnv("20", 6000000, "CNVpytor_del3", "<DEL>", 6100000, "1e-08", "0.9", "0/1") // q0 too high
        + &cnv("20", 7000000, "CNVpytor_loh4", "<LOH>", 7100000, "1e-08", "0.0", "0/1") // LOH
        + &cnv("1", 500000, "CNVpytor_dup5", "<DUP>", 600000, "1e-05", "0.0", "0/1") // added
}

fn opts(base: &str) -> Opts {
    Opts {
        base: base.to_owned(),
        base_name: "manta".to_owned(),
        delly: None,
        sex: Sex::Unknown,
        par_bed: None,
        cnvpytor: None,
        max_eval: 0.01,
        max_q0: 0.5,
        min_overlap: 0.5,
    }
}

/// `og sv add-callers`: the VCF it writes and the note it prints.
fn add(o: &Opts) -> (String, String) {
    let mut out = Vec::new();
    let note = add_callers::run(o, &mut out).unwrap_or_else(|e| panic!("og sv add-callers: {e}"));
    (String::from_utf8(out).unwrap(), note)
}

fn ids(vcf: &str) -> Vec<String> {
    records(vcf).into_iter().map(|r| r[2].clone()).collect()
}

#[test]
fn add_callers_adds_what_manta_lacks() {
    let w = Scratch::new("callers");
    let base = w.put("base.vcf", BASE);
    let delly = w.put_gz("delly.vcf.gz", DELLY);
    let cnvp = w.put_gz("cnv.vcf.gz", &cnvpytor());
    assert!(!BASE.contains("CNVpytor_del1"), "control: Manta's file alone lacks the CNVpytor deletion");

    let all = Opts {
        delly: Some(delly.clone()),
        cnvpytor: Some(cnvp.clone()),
        ..opts(&base)
    };
    let (out, note) = add(&all);
    assert_eq!(
        ids(&out),
        ["CNVpytor_dup5", "MantaDEL:1", "DEL00000002", "MantaINV:2", "CNVpytor_del1"],
        "the added records are exactly Delly's DEL00000002 and CNVpytor's del1 and dup5, sorted"
    );
    let by = by_id(&records(&out));
    assert!(
        ["MantaDEL:1", "MantaINV:2"]
            .iter()
            .all(|i| info_of(&by[*i], "SVCALLER").as_deref() == Some("manta")),
        "Manta's records name their caller"
    );
    let d = &by["DEL00000002"];
    assert!(
        info_of(d, "SVCALLER").as_deref() == Some("delly")
            && d[4] == "<DEL>"
            && info_of(d, "END").as_deref() == Some("2050000")
            && d[8..] == ["GT", "1/1"],
        "Delly's call: SVCALLER=delly, <DEL>, END kept, its genotype: {d:?}"
    );
    let c = &by["CNVpytor_del1"];
    assert!(
        info_of(c, "SVCALLER").as_deref() == Some("cnvpytor")
            && info_of(c, "END").as_deref() == Some("3700000")
            && !info(c).contains_key("SVLEN")
            && info(c).get("IMPRECISE") == Some(&None)
            && c[8..] == ["GT", "0/1"],
        "CNVpytor's deletion: SVCALLER=cnvpytor, END, no SVLEN (its positive one), IMPRECISE, its genotype: {c:?}"
    );
    let c = &by["CNVpytor_dup5"];
    assert!(
        c[4] == "<DUP>" && info_of(c, "SVTYPE").as_deref() == Some("DUP"),
        "CNVpytor's duplication is a <DUP>: {c:?}"
    );
    assert!(out.contains("##INFO=<ID=SVCALLER,"), "the header declares SVCALLER");
    assert!(
        note.contains("delly: 2 PASS DEL/DUP calls, 1 added")
            && note.contains("cnvpytor: 3 DEL/DUP calls with e-val < 0.01 and q0 < 0.5, 2 added"),
        "the summary counts what each caller gave: {note}"
    );

    // A Manta record that is filtered or genotyped 0/0 is no call: it must
    // not hide the Delly PASS deletion it overlaps.
    let manta9 = |filt: &str, gt: &str| {
        format!("chr20\t2000100\tMantaDEL:9\tN\t<DEL>\t500\t{filt}\tEND=2049900;SVTYPE=DEL;SVLEN=-49800\tGT\t{gt}\n")
    };
    for (label, gt, filt) in [("a 0/0", "0/0", "PASS"), ("a MinGQ", "0/1", "MinGQ"), ("a ./.", "./.", "PASS")] {
        let b = w.put("base_nocall.vcf", &(BASE.to_owned() + &manta9(filt, gt)));
        let got = ids(&add(&Opts { delly: Some(delly.clone()), ..opts(&b) }).0);
        assert!(
            got.iter().any(|i| i == "MantaDEL:9") && got.iter().any(|i| i == "DEL00000002"),
            "{label} Manta deletion stays and does not hide Delly's DEL00000002: {got:?}"
        );
    }
    let b = w.put("base_nocall.vcf", &(BASE.to_owned() + &manta9("PASS", "0|1")));
    let got = ids(&add(&Opts { delly: Some(delly.clone()), ..opts(&b) }).0);
    assert!(
        !got.iter().any(|i| i == "DEL00000002"),
        "control: a PASS 0|1 Manta deletion still covers Delly's DEL00000002: {got:?}"
    );
    let sites: String = (BASE.to_owned()
        + "chr20\t2000100\tMantaDEL:9\tN\t<DEL>\t500\tPASS\tEND=2049900;SVTYPE=DEL\tGT\t0/0\n")
        .lines()
        .map(|l| l.split('\t').take(8).collect::<Vec<_>>().join("\t") + "\n")
        .collect();
    let b = w.put("sites.vcf", &sites);
    let got = ids(&add(&Opts { delly: Some(delly.clone()), ..opts(&b) }).0);
    assert!(
        !got.iter().any(|i| i == "DEL00000002"),
        "a sites-only base deletion covers Delly's DEL00000002 (no genotype to read): {got:?}"
    );

    // A Delly PASS deletion genotyped 0/0 or ./. is no call either: it is
    // not added, and CNVpytor's deletion at the same place is added instead.
    for gt in ["0/0", "./."] {
        let d = w.put_gz("delly_nocall.vcf.gz", &DELLY.replace("\t1/1:40\n", &format!("\t{gt}:40\n")));
        let got = ids(&add(&Opts {
            delly: Some(d),
            cnvpytor: Some(cnvp.clone()),
            ..opts(&base)
        })
        .0);
        assert!(
            !got.iter().any(|i| i == "DEL00000002") && got.iter().any(|i| i == "CNVpytor_del0"),
            "a {gt} Delly deletion is not added and does not hide CNVpytor_del0: {got:?}"
        );
    }

    let (out, note) = add(&Opts {
        base_name: "sv_vcf".to_owned(),
        ..opts(&base)
    });
    let callers: Vec<Option<String>> = records(&out).iter().map(|r| info_of(r, "SVCALLER")).collect();
    assert!(
        note.contains("none given") && callers == [Some("sv_vcf".to_owned()), Some("sv_vcf".to_owned())],
        "no other caller: the base out with SVCALLER=sv_vcf: {note} {callers:?}"
    );
}

// --- 4. Delly's read depth ---------------------------------------------------------------

const DELLY_FMT: &str = "GT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV";

fn delly_rdcn(chrom: &str, pos: u64, rid: &str, kind: &str, end: u64, rdcn: &str, precise: bool) -> String {
    let flag = if precise { "PRECISE" } else { "IMPRECISE" };
    format!(
        "{chrom}\t{pos}\t{rid}\tN\t<{kind}>\t60\tPASS\t{flag};SVTYPE={kind};END={end}\t{DELLY_FMT}\t\
0/1:-10,0,-10:40:PASS:100:100:100:{rdcn}:10:10:0:0\n"
    )
}

fn depth() -> String {
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n".to_owned()
        + &delly_rdcn("chr20", 5000000, "DEL_rd2", "DEL", 55000000, "2", false) // 50 Mb, no loss: set aside
        + &delly_rdcn("chr1", 10000000, "DEL_rd1", "DEL", 60000000, "1", false) // 50 Mb, one copy: added
        + &delly_rdcn("chr1", 100000000, "DUP_rd2", "DUP", 100500000, "2", false) // no gain: set aside
        + &delly_rdcn("chr1", 120000000, "DUP_rd3", "DUP", 120500000, "3", false) // a gain: added
        + &delly_rdcn("chr1", 150000000, "DEL_short", "DEL", 150000500, "2", true) // 500 bp: added
        + &delly_rdcn("chr1", 160000000, "DEL_rdneg", "DEL", 161000000, "-1", false) // no estimate: set aside
        + &delly_rdcn("chrX", 50000000, "DEL_x1", "DEL", 51000000, "1", false) // male: normal, set aside
        + &delly_rdcn("chrX", 60000000, "DUP_x2", "DUP", 60500000, "2", false) // male: a gain; female: none
        + &delly_rdcn("chrX", 1000000, "DEL_par1", "DEL", 1500000, "1", false) // PAR: a loss for anyone
}

#[test]
fn add_callers_checks_delly_read_depth() {
    let w = Scratch::new("depth");
    let base = w.put("base.vcf", BASE);
    let delly = w.put_gz("delly_rdcn.vcf.gz", &depth());
    let par = repo().join("assets/par_grch38.bed").to_string_lossy().into_owned();
    let planted = records(&depth());
    assert!(
        planted.len() == 9 && planted.iter().all(|r| r[6] == "PASS" && (r[4] == "<DEL>" || r[4] == "<DUP>")),
        "control: every planted Delly call is PASS DEL or DUP"
    );

    let added = |sex: Sex, par: Option<&str>| {
        let (out, note) = add(&Opts {
            delly: Some(delly.clone()),
            sex,
            par_bed: par.map(str::to_owned),
            ..opts(&base)
        });
        let ids: BTreeSet<String> = records(&out)
            .iter()
            .filter(|r| info_of(r, "SVCALLER").as_deref() == Some("delly"))
            .map(|r| r[2].clone())
            .collect();
        assert!(!ids.is_empty(), "og sv add-callers adds Delly calls: {note}");
        (ids, out, note)
    };
    let has = |ids: &BTreeSet<String>, i: &str| ids.contains(i);

    // No sex and no PAR file first: the arguments step 05 gave before.
    let (ids, _, _) = added(Sex::Unknown, None);
    assert!(!has(&ids, "DEL_rd2"), "no sex: the 50 Mb deletion with RDCN 2 is set aside: {ids:?}");
    assert!(!has(&ids, "DUP_rd2"), "no sex: the duplication with RDCN 2 is set aside: {ids:?}");
    assert!(!has(&ids, "DEL_rdneg"), "no sex: the long deletion with RDCN -1 is set aside: {ids:?}");
    assert!(
        ["DEL_rd1", "DUP_rd3", "DEL_short"].iter().all(|i| has(&ids, i)),
        "no sex: the deletion with RDCN 1, the duplication with RDCN 3 and the short call are added: {ids:?}"
    );
    assert!(has(&ids, "DEL_x1"), "no sex: the chrX deletion with RDCN 1 is added: {ids:?}");
    assert!(has(&ids, "DUP_x2"), "no sex: the chrX duplication with RDCN 2 is added: {ids:?}");

    let (ids, out, note) = added(Sex::Male, Some(&par));
    let want: BTreeSet<String> = ["DEL_rd1", "DUP_rd3", "DEL_short", "DUP_x2", "DEL_par1"]
        .iter()
        .map(|s| (*s).to_owned())
        .collect();
    assert_eq!(
        ids, want,
        "male: RDCN 2 deletion and duplication, RDCN -1 and the chrX RDCN 1 deletion (his normal) set aside; \
         RDCN 1 deletion, RDCN 3 duplication, the short call, the chrX RDCN 2 duplication and the PAR RDCN 1 \
         deletion added"
    );
    let line = "delly: 5 PASS DEL/DUP calls, 5 added, and 4 PASS calls longer than 1000 bp set aside";
    assert!(note.contains(line), "male: the summary counts the 4 calls set aside: {note}");
    assert!(
        out.contains(&format!("##svcallers_added=\"Callers beside manta: {line}"))
            && out.find("##svcallers_added=") < out.find("#CHROM"),
        "male: the input VCF's header states the same, set-aside count included"
    );

    let (ids, _, _) = added(Sex::Female, Some(&par));
    assert!(!has(&ids, "DUP_x2"), "female: the chrX duplication with RDCN 2 is set aside: {ids:?}");
    assert!(has(&ids, "DEL_x1"), "female: the chrX deletion with RDCN 1 is added: {ids:?}");
}
