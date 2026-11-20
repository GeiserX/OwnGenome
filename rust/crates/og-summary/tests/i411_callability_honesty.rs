//! A gene a step of its own calls (Parascopy, pypgx, the copy-number
//! engine) and that step 16b did not read is listed under `not_read_owned`
//! with its owner. Its entry also says which step that is and whether it ran
//! for this sample (`owner_step`, `owner_status`: ran, not run, or from an
//! earlier run, the collector's stale rule), so a gene whose own step never
//! ran no longer reads as handed to a step that answered for it. A Parascopy
//! gene is judged by its own table. The share keeps its four decimals:
//! 0.8988 stays 0.8988.

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

/// The `table` column of step 35's Parascopy catalog: one table per row.
fn parascopy_tables() -> Vec<String> {
    let text = fs::read_to_string(format!("{REPO}/assets/paralog_loci.tsv")).unwrap();
    let mut rows = text.lines().filter(|l| !l.starts_with('#'));
    let head: Vec<&str> = rows.next().unwrap().split('\t').collect();
    let col = head.iter().position(|h| *h == "table").unwrap();
    rows.map(|l| l.split('\t').nth(col).unwrap().to_owned())
        .collect()
}

#[test]
fn an_owned_gene_whose_step_ran_says_ran() {
    // the Nextflow layout: the table under coverage/
    let d = folder("ran", "coverage", |_, _| 0);
    for t in parascopy_tables() {
        d.put(&format!("paralogs/S_{t}_copy_number.tsv"), b"x\n");
    }
    for f in [
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
fn a_parascopy_gene_is_judged_by_its_own_table() {
    // GBA1's table is there, SMN1's is not: step 35 answered GBA1 only.
    let d = folder("own-table", "mosdepth", |_, _| 0);
    d.put("paralogs/S_gba1_copy_number.tsv", b"x\n");
    let got = owned(&collect(&d));
    let status = |g: &str| got.iter().find(|r| r.0 == g).map(|r| r.3.clone()).unwrap();
    assert_eq!(status("GBA1"), "ran");
    assert_eq!(status("SMN1"), "not run");
    assert_eq!(status("SMN2"), "not run");
}

#[test]
fn a_table_from_an_earlier_run_of_a_failed_step_says_so() {
    // step 35 failed in this run; SMN1's table is left from an earlier run,
    // GBA1's was written in this run before the step failed. The collector
    // calls such a section stale; the owned gene says the same.
    let d = folder("earlier", "mosdepth", |_, _| 0);
    let now = std::time::SystemTime::now();
    let started = now.duration_since(std::time::UNIX_EPOCH).unwrap().as_secs() - 60;
    d.put(
        "logs/run_status.tsv",
        format!("meta\tstarted_epoch\t{started}\nstep\t35\tfailed (exit 1)\n").as_bytes(),
    );
    d.put("paralogs/S_smn_copy_number.tsv", b"x\n");
    d.put("paralogs/S_gba1_copy_number.tsv", b"x\n");
    fs::File::options()
        .write(true)
        .open(d.0.join("paralogs/S_smn_copy_number.tsv"))
        .unwrap()
        .set_modified(now - std::time::Duration::from_secs(3600))
        .unwrap();
    let got = owned(&collect(&d));
    let status = |g: &str| got.iter().find(|r| r.0 == g).map(|r| r.3.clone()).unwrap();
    assert_eq!(status("SMN1"), "from an earlier run");
    assert_eq!(status("SMN2"), "from an earlier run");
    assert_eq!(status("GBA1"), "ran");
    assert_eq!(status("CYP2D6"), "not run");
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
