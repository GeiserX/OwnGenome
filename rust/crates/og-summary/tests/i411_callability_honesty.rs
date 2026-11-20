//! A gene a step of its own calls (Parascopy, pypgx, the copy-number
//! engine) and that step 16b did not read is listed under `not_read_owned`
//! with its owner. Its entry also says which step that is and whether it ran
//! for this sample (`owner_step`, `owner_status`), so a gene whose own step
//! never ran no longer reads as handed to a step that answered for it. The
//! share keeps its four decimals: 0.8988 stays 0.8988.

use og_summary::collect::{Group, Input, post};
use pycompat::json::Value;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!(
            "og-callability-honesty-{}-{tag}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn put(&self, rel: &str, body: &[u8]) {
        let p = self.0.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(&p, body).unwrap();
    }
    fn put_gz(&self, rel: &str, text: &str) {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        e.write_all(text.as_bytes()).unwrap();
        self.put(rel, &e.finish().unwrap());
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The shipped catalog's rows: gene, chrom, start, end, owner.
fn catalog() -> Vec<[String; 5]> {
    let text = fs::read_to_string(format!("{REPO}/assets/callability_genes.tsv")).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            [f[0], f[2], f[3], f[4], f[6]].map(str::to_owned)
        })
        .collect()
}

/// A sample folder with step 16b's table at `table` (under `mosdepth/` or
/// `coverage/`): every catalog gene with a position, `ok` of its bases at
/// 10x, as a share of each gene's length.
fn folder(tag: &str, table: &str, ok: impl Fn(&str, i64) -> i64) -> Dir {
    let d = Dir::new(tag);
    let mut bed = String::from("#chrom\tstart\tend\tregion\t10X\n");
    for [gene, chrom, start, end, _] in catalog() {
        if start == "." {
            continue;
        }
        let (s, e): (i64, i64) = (start.parse().unwrap(), end.parse().unwrap());
        bed.push_str(&format!(
            "{chrom}\t{s}\t{e}\t{gene}\t{}\n",
            ok(&gene, e - s)
        ));
    }
    d.put_gz(&format!("{table}/S_callability.bed.gz"), &bed);
    d
}

fn collect(d: &Dir) -> Value {
    post::collect(&Input {
        sample: "S".into(),
        sample_dir: d.0.to_str().unwrap().to_owned(),
        declared_sex: None,
        assets: format!("{REPO}/assets"),
        group: Group::A,
        tables: format!("{REPO}/bin"),
    })
    .unwrap()
}

fn get<'a>(v: &'a Value, key: &str) -> &'a Value {
    match v {
        Value::Object(m) => m.get(key).unwrap_or_else(|| panic!("no {key} in {v:?}")),
        _ => panic!("not an object: {v:?}"),
    }
}

fn items(v: &Value) -> &[Value] {
    match v {
        Value::Array(a) => a,
        _ => panic!("not a list: {v:?}"),
    }
}

fn st(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => panic!("not a string: {v:?}"),
    }
}

/// `not_read_owned` of the summary: (gene, owner, owner_step, owner_status).
fn owned(v: &Value) -> Vec<(String, String, String, String)> {
    let data = get(get(get(v, "sections"), "callability"), "data");
    items(get(data, "not_read_owned"))
        .iter()
        .map(|g| {
            (
                st(get(g, "gene")).to_owned(),
                st(get(g, "owner")).to_owned(),
                st(get(g, "owner_step")).to_owned(),
                st(get(g, "owner_status")).to_owned(),
            )
        })
        .collect()
}

/// The owned genes of the catalog that have a position.
fn owned_in_catalog() -> usize {
    catalog()
        .iter()
        .filter(|r| r[4] != "." && r[2] != ".")
        .count()
}

#[test]
fn an_owned_gene_whose_step_never_ran_says_not_run() {
    let d = folder("never", "mosdepth", |_, _| 0);
    let got = owned(&collect(&d));
    assert_eq!(got.len(), owned_in_catalog());
    assert!(got.len() >= 50, "the catalog owns {} genes", got.len());
    for (gene, owner, step, status) in &got {
        assert_eq!(status, "not run", "{gene} ({owner}, step {step})");
    }
    let step = |g: &str| got.iter().find(|r| r.0 == g).map(|r| r.2.clone());
    assert_eq!(step("SMN1").as_deref(), Some("35"));
    assert_eq!(step("CYP2D6").as_deref(), Some("32"));
    assert_eq!(step("HBA1").as_deref(), Some("35"));
    assert_eq!(step("LPA").as_deref(), Some("46"));
    assert_eq!(step("RHD").as_deref(), Some("41"));
}

#[test]
fn an_owned_gene_whose_step_ran_says_ran() {
    // the Nextflow layout: the table under coverage/
    let d = folder("ran", "coverage", |_, _| 0);
    for f in [
        "paralogs/S_smn_copy_number.tsv",
        "paralogs/S_hba_copy_number.tsv",
        "paralogs/S_cyp21a2_psv.tsv",
        "paralogs/S_cn_dosage.tsv",
        "pypgx/S_pypgx_summary.tsv",
        "lpa/S_lpa_kiv2.tsv",
        "bloodtype/S_bloodtype.tsv",
    ] {
        d.put(f, b"x\n");
    }
    let got = owned(&collect(&d));
    assert_eq!(got.len(), owned_in_catalog());
    for (gene, owner, step, status) in &got {
        assert_eq!(status, "ran", "{gene} ({owner}, step {step})");
    }
}

#[test]
fn each_owner_is_judged_by_its_own_file() {
    // step 35's HBA table and step 35b's Paraphase table share paralogs/, and
    // neither is a Parascopy answer; CPIC's comparison shares pypgx/.
    let d = folder("own-file", "mosdepth", |_, _| 0);
    d.put("paralogs/S_hba_copy_number.tsv", b"x\n");
    d.put("paralogs/S_paraphase_copy_number.tsv", b"x\n");
    d.put("pypgx/S_pharmcat_comparison.tsv", b"x\n");
    let got = owned(&collect(&d));
    let status = |g: &str| got.iter().find(|r| r.0 == g).map(|r| r.3.clone()).unwrap();
    assert_eq!(status("HBA1"), "ran");
    assert_eq!(status("HBA2"), "ran");
    assert_eq!(status("SMN1"), "not run");
    assert_eq!(status("CYP2D6"), "not run");
    assert_eq!(status("LPA"), "not run");
    for (gene, owner, _, status) in &got {
        if owner == "parascopy" || owner == "pypgx" {
            assert_eq!(status, "not run", "{gene} ({owner})");
        }
    }
}

#[test]
fn a_share_under_ninety_percent_keeps_four_decimals() {
    // 8988 of 10000 bases: 0.8988, not 0.9 and not read
    let d = folder("share", "mosdepth", |g, n| {
        if g == "SMN1" { n * 8988 / 10000 } else { 0 }
    });
    let v = collect(&d);
    let data = get(get(get(&v, "sections"), "callability"), "data");
    let smn1 = items(get(data, "not_read_owned"))
        .iter()
        .find(|g| st(get(g, "gene")) == "SMN1")
        .expect("SMN1 is listed");
    let [_, _, start, end, _] = catalog().into_iter().find(|r| r[0] == "SMN1").unwrap();
    let n: i64 = end.parse::<i64>().unwrap() - start.parse::<i64>().unwrap();
    let want = ((n * 8988 / 10000) as f64 / n as f64 * 10000.0).round() / 10000.0;
    assert_eq!(want, 0.8988, "SMN1 is {n} bases long");
    assert_eq!(get(smn1, "share"), &Value::Float(0.8988));
    assert_eq!(st(get(smn1, "status")), "not read");
}
