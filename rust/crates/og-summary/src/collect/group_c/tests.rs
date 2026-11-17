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
fn chip_nan_depth_is_unreadable() {
    let d = Dir::new("c-chip-nan");
    d.put(
        "S/chip/S_chip.tsv",
        "gene\tprotein_change\tvaf\talt_reads\tdepth\tstatus\tnote\n",
    );
    d.put("S/chip/S_chip_depth.tsv", "gene\tmean_depth\tdetection_limit_vaf\tstatus\nA\t30\t.\tcovered\nB\tnan\t.\tcovered\nC\t40\t.\tcovered\n");
    let v = summary(&d, None);
    let chip = at(&v, &["sections", "chip"]);
    assert_eq!(text(at(chip, &["state"])), "unreadable");
    assert_eq!(at(chip, &["data"]), &Value::Object(Default::default()));
}

#[test]
fn chip_infinity_only_fails_when_the_median_is_nonfinite() {
    let d = Dir::new("c-chip-inf");
    d.put(
        "S/chip/S_chip.tsv",
        "gene\tprotein_change\tvaf\talt_reads\tdepth\tstatus\tnote\n",
    );
    for (depths, median) in [
        (vec!["30", "inf", "40"], Some(40.0)),
        (vec!["30", "-inf", "40"], Some(30.0)),
        (vec!["inf"], None),
        (vec!["-inf"], None),
        (vec!["-inf", "inf"], None),
    ] {
        let mut body = "gene\tmean_depth\tdetection_limit_vaf\tstatus\n".to_owned();
        for (i, value) in depths.iter().enumerate() {
            body.push_str(&format!("G{i}\t{value}\t.\tcovered\n"));
        }
        d.put("S/chip/S_chip_depth.tsv", body);
        let v = summary(&d, None);
        let chip = at(&v, &["sections", "chip"]);
        if let Some(median) = median {
            assert_eq!(text(at(chip, &["state"])), "ok", "{depths:?}");
            assert_eq!(at(chip, &["data", "median_depth"]), &Value::Float(median));
        } else {
            assert_eq!(text(at(chip, &["state"])), "unreadable", "{depths:?}");
            assert_eq!(at(chip, &["data"]), &Value::Object(Default::default()));
        }
    }
}

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
    // Invented calls. Keep FORMAT separate from the values, as in the parity generator.
    let vcf = concat!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        "chr7\t1\t.\tT\t<STR11>\t.\tPASS\tVARID=CFTR_TG;REPID=CFTR_TG\t",
        "GT:SO:REPCN:REPCI:ADSP\t",
        "0/1:SPANNING/SPANNING:11/12:11-11/12-12:8/9\n",
        "chr7\t2\t.\tT\t<STR5>\t.\tPASS\tVARID=CFTR_T;REPID=CFTR_T\t",
        "GT:SO:REPCN:REPCI:ADSP\t",
        "0/1:SPANNING/SPANNING:5/7:5-5/7-7:8/9\n",
    );
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

/// Every locus of the checkout's `bin/str_intermediate_ranges.tsv` is read
/// (tests/test_str_ranges.py checked the Python reader the same way).
#[test]
fn every_row_of_the_repeat_ranges_table_is_read() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../../bin/str_intermediate_ranges.tsv"
    );
    let rows = super::util::bin_rows(path).unwrap();
    let text = std::fs::read_to_string(path).unwrap();
    let loci: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .skip(1)
        .map(|l| l.split('\t').next().unwrap_or(""))
        .collect();
    let read: Vec<&str> = rows.iter().map(|r| r["locus"].as_str()).collect();
    assert_eq!(read, loci);
    assert!(loci.len() >= 3, "{loci:?}");
}

// i411-roh-scope: ST rows describe assessed chromosomes, including ones
// with no long run. A run list alone cannot establish whole-genome scope.
fn roh_scope_summary(sites: &[u8], runs: &str) -> Value {
    let d = Dir::new("i411-roh-scope");
    let mut body = String::new();
    for c in sites {
        body.push_str(&format!("ST\tS\tchr{c}\t40000000\t0\t60\n"));
    }
    body.push_str(runs);
    d.put("S/roh/S_roh.txt", body);
    at(&summary(&d, None), &["sections", "roh", "data"]).clone()
}

#[test]
fn roh_scope_chr20_uses_its_own_unmasked_denominator() {
    let r = roh_scope_summary(&[20], "RG\tS\tchr20\t1000001\t22000000\t21000000\n");
    assert_eq!(at(&r, &["f_roh"]), &Value::Float(0.3515));
    assert_eq!(at(&r, &["relationship"]), &Value::Null);
    assert_eq!(at(&r, &["upd_flag"]), &Value::Null);
    assert_eq!(at(&r, &["scope", "unmasked_bp"]), &Value::Int(59_744_167));
    assert_eq!(
        at(&r, &["scope", "autosomes"]),
        &Value::Array(vec![Value::Str("chr20".into())])
    );
}

#[test]
fn roh_scope_counts_sites_without_runs_and_keeps_partial_labels_empty() {
    let r = roh_scope_summary(
        &[22, 20, 20],
        "RG\tS\tchr20\t1000001\t9000000\t8000000\nRG\tS\tchr22\t20000001\t28000000\t8000000\n",
    );
    assert_eq!(at(&r, &["f_roh"]), &Value::Float(0.1717));
    assert_eq!(at(&r, &["relationship"]), &Value::Null);
    assert_eq!(at(&r, &["scope", "unmasked_bp"]), &Value::Int(93_162_635));
    let mixed = roh_scope_summary(&[20, 22], "RG\tS\tchr20\t1000001\t9000000\t8000000\n");
    assert_eq!(at(&mixed, &["f_roh"]), &Value::Float(0.0859));
    assert_eq!(
        at(&mixed, &["scope", "unmasked_bp"]),
        &Value::Int(93_162_635)
    );
    let no_runs = roh_scope_summary(&[20], "");
    assert_eq!(at(&no_runs, &["f_roh"]), &Value::Float(0.0));
    assert_eq!(at(&no_runs, &["relationship"]), &Value::Null);
    // Value::Float(-0.0) == Value::Float(0.0), so check the sign apart: a
    // sample with no long run must read 0.0, never -0.0, in both reports.
    for key in ["f_roh", "long_mb"] {
        assert!(
            matches!(at(&no_runs, &[key]), Value::Float(f) if f.is_sign_positive()),
            "{key} must be +0.0 without long runs, got {:?}",
            at(&no_runs, &[key])
        );
    }
}

#[test]
fn roh_scope_sparse_all_autosomes_keep_the_whole_genome_band() {
    let sites: Vec<u8> = (1..=22).collect();
    let r = roh_scope_summary(&sites, "RG\tS\tchr2\t10000001\t18000000\t8000000\n");
    assert_eq!(at(&r, &["f_roh"]), &Value::Float(0.003));
    assert_eq!(
        at(&r, &["relationship"]),
        &Value::Str("no sign of related parents".into())
    );
    assert_eq!(
        at(&r, &["scope", "unmasked_bp"]),
        &Value::Int(2_658_001_522)
    );
    assert_eq!(at(&r, &["scope", "partial"]), &Value::Bool(false));
}

#[test]
fn roh_scope_a_segment_only_file_never_establishes_scope() {
    let r = roh_scope_summary(&[], "RG\tS\tchr20\t1000001\t22000000\t21000000\n");
    assert_eq!(at(&r, &["f_roh"]), &Value::Float(0.3515));
    assert_eq!(at(&r, &["relationship"]), &Value::Null);
    assert_eq!(at(&r, &["upd_flag"]), &Value::Null);
    assert_eq!(
        at(&r, &["scope", "source"]),
        &Value::Str("segments_only".into())
    );
    let empty = roh_scope_summary(&[], "# empty\n");
    assert_eq!(at(&empty, &["f_roh"]), &Value::Null);
}

#[test]
fn roh_scope_upd_is_a_region_hint_without_syndrome_names() {
    let sites: Vec<u8> = (1..=22).collect();
    // Distal chr15 and chr11 runs do not establish the imprinted locus's
    // parental origin, but mixed iso/heterodisomy can involve that locus.
    for (chrom, start, end, length, region) in [
        (
            15,
            40000001,
            101991189,
            61991189,
            "15q11-q13 (SNRPN, UBE3A)",
        ),
        (
            11,
            60000001,
            100000000,
            40000000,
            "11p15 (H19/IGF2, KCNQ1OT1)",
        ),
        (20, 1000001, 22000000, 21000000, "20q13.32 (GNAS)"),
    ] {
        let r = roh_scope_summary(
            &sites,
            &format!("RG\tS\tchr{chrom}\t{start}\t{end}\t{length}\n"),
        );
        let u = at(&r, &["upd_flag"]);
        assert_eq!(
            at(u, &["imprinted", "0", "region"]),
            &Value::Str(region.into())
        );
        assert_eq!(
            at(u, &["imprinted", "0", "maternal_upd"]),
            &Value::Str("-".into())
        );
        assert_eq!(
            at(u, &["imprinted", "0", "paternal_upd"]),
            &Value::Str("-".into())
        );
        let note = text(at(u, &["note"]));
        assert!(note.contains("deletion") && note.contains("acquired") && note.contains("parent"));
        assert_eq!(at(&r, &["relationship"]), &Value::Null);
    }
}
