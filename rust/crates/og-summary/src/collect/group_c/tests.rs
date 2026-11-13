//! Group c through `collect`: each section's state, and the post pass
//! (CNVpytor's chrY filter, Locityper's primary calls, the repeat readings).
//! The section readers are checked against the Python collector by the
//! summary-collect parity cases (tests/parity/cases.tsv).

use super::super::{Group, Input, collect};
use super::{GROUP_C, Others};
use crate::testutil::Dir;
use pycompat::json::Value;

fn summary(d: &Dir, declared_sex: Option<&str>) -> Value {
    let input = Input {
        sample: "S".into(),
        sample_dir: d.path("S"),
        declared_sex: declared_sex.map(str::to_owned),
        assets: d.path("assets"),
        group: Group::C,
        tables: concat!(env!("CARGO_MANIFEST_DIR"), "/../../../bin").into(),
    };
    collect(&input).unwrap().to_json()
}

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no {k} in {v:?}")),
        Value::Array(a) => &a[k.parse::<usize>().unwrap()],
        _ => panic!("{k} of {v:?}"),
    })
}

fn text(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => panic!("not a text: {v:?}"),
    }
}

const CNVS: &str = "deletion\tchrY:2800001-2900000\t100000\t0.0\t1e-15\t1\t1\t1\t0.05\t0\t9\n\
deletion\tchr2:10000001-10100000\t100000\t0.4\t1e-9\t1\t1\t1\t0.01\t0\t9\n\
deletion\tchr1:122000001-123000000\t1000000\t0.1\t1e-12\t1\t1\t1\t0.2\t0\t9\n";

#[test]
fn every_section_of_group_c_is_written_in_order() {
    let d = Dir::new("c-order");
    let v = summary(&d, None);
    let Value::Object(secs) = at(&v, &["sections"]) else {
        panic!()
    };
    let keys: Vec<&str> = secs.keys().map(String::as_str).collect();
    let want: Vec<&str> = GROUP_C.iter().map(|s| s.key).collect();
    assert_eq!(keys, want);
    assert!(secs.values().all(|s| text(at(s, &["state"])) == "missing"));
}

#[test]
fn cnvpytor_drops_chry_for_a_declared_female_only() {
    let d = Dir::new("c-cnv");
    d.put("S/cnvpytor/S_cnvs.txt", CNVS);
    let female = summary(&d, Some("female"));
    let data = at(&female, &["sections", "cnvpytor", "data"]);
    assert_eq!(at(data, &["significant_filtered"]), &Value::Int(1));
    assert_eq!(at(data, &["dropped", "chry_female"]), &Value::Int(1));
    assert_eq!(at(data, &["dropped", "band"]), &Value::Int(1));
    let male = summary(&d, Some("male"));
    let data = at(&male, &["sections", "cnvpytor", "data"]);
    assert_eq!(at(data, &["significant_filtered"]), &Value::Int(2));
    assert_eq!(at(data, &["declared_sex"]), &Value::Str("male".into()));
}

#[test]
fn locityper_names_the_primary_call_of_each_locus() {
    let d = Dir::new("c-loci");
    d.put(
        "S/hla_t1k/S_hla_genotype.tsv",
        "HLA-DRB1\t2\tDRB1*03:01\t1\t40\tDRB1*07:01\t1\t0\n",
    );
    let head = "locus\tregion\tstatus\tgenotype\thaplotype_1\thaplotype_2\tquality\ttotal_reads\tunexplained_reads\tweighted_distance\twarnings\tfilter\treason\n";
    d.put(
        "S/locityper/S_locityper.tsv",
        format!(
            "{head}HLA-DRB1\tx\tcalled\th\th\th\t1\t1\t1\t1\t\tpass\t\n\
             HLA-DQB1\tx\tcalled\th\th\th\t1\t1\t1\t1\t\tpass\t\n\
             CYP2D6\tx\tcalled\th\th\th\t1\t1\t1\t1\t\tpass\t\n\
             OTHER\tx\tcalled\th\th\th\t1\t1\t1\t1\t\tpass\t\n"
        ),
    );
    let v = summary(&d, None);
    let loci = at(&v, &["sections", "locityper", "data", "loci"]);
    let primary = |i: &str| text(at(loci, &[i, "primary"])).to_owned();
    assert_eq!(
        primary("0"),
        "T1K (step 08): DRB1*03:01/DRB1*07:01 (low confidence)"
    );
    assert_eq!(primary("1"), "T1K (step 08) did not run");
    assert_eq!(
        primary("2"),
        "pypgx (step 32) and Cyrius (step 21) did not run"
    );
    assert_eq!(primary("3"), "no other step calls this locus");
}

#[test]
fn locityper_reads_pypgx_and_cyrius_from_group_b() {
    let mut v = pycompat::json::loads(r#"{"loci": [{"locus": "CYP2D6"}]}"#).unwrap();
    let others = Others {
        pypgx: Some((super::State::Stale, "*1/*4".into())),
        cyrius: Some((super::State::Ok, Some("*1/*2".into()), None)),
        ..Others::default()
    };
    super::tables::locityper_post(&mut v, super::State::Missing, None, &others);
    assert_eq!(
        text(at(&v, &["loci", "0", "primary"])),
        "pypgx (step 32): *1/*4 (stale: from an earlier run); Cyrius (step 21): *1/*2 (not usable: Filter missing)"
    );
}

#[test]
fn a_cftr_5t_with_r117h_reads_the_phase() {
    let d = Dir::new("c-cftr");
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
chr7\t1\t.\tT\t<STR11>\t.\tPASS\tVARID=CFTR_TG;REPID=CFTR_TG\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:11/12:11-11/12-12:8/9\n\
chr7\t2\t.\tT\t<STR5>\t.\tPASS\tVARID=CFTR_T;REPID=CFTR_T\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:5/7:5-5/7-7:8/9\n";
    d.put("S/expansion_hunter/S_eh.vcf", vcf);
    d.put(
        "S/traits/S_traits.tsv",
        "rsid\tgene\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids\n\
rs78655421\tCFTR\tR117H\tG/A\tcalled\tr\ta\tn\tp\n",
    );
    let v = summary(&d, None);
    let cftr = at(&v, &["sections", "expansions", "data", "readings", "0"]);
    assert_eq!(text(at(cftr, &["phase"])), "unknown");
    assert_eq!(text(at(cftr, &["r117h"])), "G/A");
    assert!(text(at(cftr, &["reading"])).contains("phase unknown"));
}

#[test]
fn a_stale_section_says_from_which_run() {
    let d = Dir::new("c-stale");
    d.put(
        "S/logs/run_status.tsv",
        "meta\tstarted_epoch\t4102444800\nstep\t11\tskipped\n",
    );
    d.put(
        "S/vcf/S_roh.txt",
        "RG\tS\tchr2\t10000001\t17000000\t7000000\t9\t9\n",
    );
    let v = summary(&d, None);
    let roh = at(&v, &["sections", "roh"]);
    assert_eq!(text(at(roh, &["state"])), "stale");
    assert!(text(at(roh, &["note"])).contains("step 11 was skipped"));
    assert_eq!(at(roh, &["data", "largest_mb"]), &Value::Float(7.0));
}
