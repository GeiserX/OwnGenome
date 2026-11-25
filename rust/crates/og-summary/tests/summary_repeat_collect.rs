//! `og summary collect` writes the same summary twice from the same inputs.
//! Only two fields follow time: `generated_utc` (the clock when the
//! summary is written) and each section's `file_date` (the minute its file
//! was written). clinvar-conflicting-1 compared two summaries whose hits
//! VCF step 06 had written on either side of a minute, and only that
//! stamp moved.
//!
//! 1. Thirty collects of one folder are byte-identical once
//!    `generated_utc` is dropped. Each `HashMap` gets its own random seed,
//!    so an array built in hash order would show here.
//! 2. The same hits VCF, written again one minute later, changes the
//!    `file_date` of the sections read from it and nothing else.

use og_summary::collect::{Group, Input, post};
use pycompat::json::{Dump, Value, dumps};
use std::path::{Path, PathBuf};
use std::time::{Duration, UNIX_EPOCH};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

const PL: &str = "CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts";

/// Step 06's hits: an ACMG gene (BRCA2), an ACMG pair (MUTYH) and a
/// carrier gene (CFTR), so the card and the carrier list both have rows.
fn hits() -> String {
    format!(
        "##fileformat=VCFv4.2\n##PGP_CALL_FLOOR=GQ=20,DP=8\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
chr13\t32340300\t1\tC\tT\t.\tPASS\tGENEINFO=BRCA2:675;CLNSIG=Pathogenic;{PL}\tGT:GQ:DP\t0/1:60:35\n\
chr1\t45331556\t2\tC\tT\t.\tPASS\tGENEINFO=MUTYH:4595;CLNSIG=Pathogenic;{PL}\tGT:GQ:DP\t0/1:50:30\n\
chr1\t45332400\t3\tG\tA\t.\tPASS\tGENEINFO=MUTYH:4595;CLNSIG=Likely_pathogenic;{PL}\tGT:GQ:DP\t0/1:55:32\n\
chr7\t117559590\t4\tATCT\tA\t.\tPASS\tGENEINFO=CFTR:1080;CLNSIG=Pathogenic;{PL}\tGT:GQ:DP\t0/1:50:40\n"
    )
}

/// Writes the hits VCF with its modification time at `secs`.
fn write_hits(path: &Path, secs: u64) {
    std::fs::write(path, hits()).unwrap();
    std::fs::File::options()
        .write(true)
        .open(path)
        .unwrap()
        .set_modified(UNIX_EPOCH + Duration::from_secs(secs))
        .unwrap();
}

/// summary.json as `og summary collect` writes it, without `generated_utc`.
fn collect(dir: &Path) -> (Value, String) {
    let input = Input {
        sample: "S".into(),
        sample_dir: dir.to_str().unwrap().to_owned(),
        declared_sex: None,
        assets: repo().join("assets").to_str().unwrap().to_owned(),
        group: Group::All,
        tables: og_summary::collect::group_c::util::default_tables(),
    };
    let mut v = post::collect(&input).unwrap();
    match &mut v {
        Value::Object(m) => assert!(m.shift_remove("generated_utc").is_some()),
        _ => panic!("summary.json is no object"),
    }
    let text = dumps(
        &v,
        &Dump {
            indent: Some(2),
            ..Dump::default()
        },
    )
    .unwrap();
    (v, text)
}

/// The paths where `a` and `b` differ, key order included.
fn diff(a: &Value, b: &Value, at: &str, out: &mut Vec<String>) {
    match (a, b) {
        (Value::Object(x), Value::Object(y)) if x.keys().eq(y.keys()) => {
            for ((k, va), vb) in x.iter().zip(y.values()) {
                diff(va, vb, &format!("{at}.{k}"), out);
            }
        }
        (Value::Array(x), Value::Array(y)) if x.len() == y.len() => {
            for (i, (va, vb)) in x.iter().zip(y).enumerate() {
                diff(va, vb, &format!("{at}[{i}]"), out);
            }
        }
        // Two objects here have other keys or the same keys in another
        // order; `==` on an object ignores key order, so it cannot judge.
        (Value::Object(_), Value::Object(_)) => out.push(at.to_owned()),
        _ if a != b => out.push(at.to_owned()),
        _ => {}
    }
}

fn get<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no {k} in {path:?}")),
        _ => panic!("{path:?} is no object at {k}"),
    })
}

fn rows(v: &Value) -> usize {
    match v {
        Value::Array(a) => a.len(),
        _ => 0,
    }
}

/// Removes the scratch folder when the test ends, also when an assertion fails.
struct Scratch(PathBuf);

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn the_same_folder_gives_the_same_summary_and_only_file_dates_follow_the_file() {
    let d = std::env::temp_dir().join(format!("og-summary-repeat-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    let _scratch = Scratch(d.clone());
    let sd = d.join("S");
    std::fs::create_dir_all(sd.join("clinvar")).unwrap();
    let hits = sd.join("clinvar/S_clinvar_hits.vcf");
    // 2026-11-25 02:48:57 UTC, three seconds before a minute.
    let minute = 1_795_574_880;
    write_hits(&hits, minute + 57);

    let (first, text) = collect(&sd);
    assert!(
        rows(get(
            &first,
            &["sections", "acmg_sf", "data", "clinvar_hits"]
        )) >= 3
    );
    assert_eq!(
        get(&first, &["sections", "carrier", "state"]),
        &Value::Str("ok".into())
    );
    for run in 2..=30 {
        let (v, t) = collect(&sd);
        let mut where_ = Vec::new();
        diff(&first, &v, "", &mut where_);
        assert!(where_.is_empty(), "collect {run} differs at {where_:?}");
        assert_eq!(t, text, "collect {run} writes other bytes (key order)");
    }

    // Step 06 runs again six seconds later: the same hits, the next minute.
    write_hits(&hits, minute + 63);
    let (again, _) = collect(&sd);
    let mut where_ = Vec::new();
    diff(&first, &again, "", &mut where_);
    where_.sort();
    assert_eq!(
        where_,
        [
            ".sections.acmg_sf.file_date",
            ".sections.carrier.file_date",
            ".sections.clinvar.file_date",
        ],
        "only the file stamps of the sections read from the hits VCF move"
    );
    assert_eq!(
        get(&again, &["sections", "acmg_sf", "file_date"]),
        &Value::Str("2026-11-25 02:49 UTC".into())
    );
}
