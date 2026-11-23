//! A gene counts as read when 90% of its whole body reaches 10x at MAPQ 20,
//! introns included, so a coding exon can be missed while the gene passes.
//! Step 16b also measures the coding exons of each catalog gene
//! (`assets/callability_exons.tsv`) into `<s>_callability_exons.bed.gz`
//! beside the gene table. A gene whose body is read but one of whose coding
//! exons is under 90% reads "exon not read" and names that exon
//! (`worst_exon`). A folder without the exon table keeps its summary as it
//! was: no gene changes and no new key appears.

use og_summary::collect::{Group, Input, post};
use pycompat::json::Value;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d =
            std::env::temp_dir().join(format!("og-callability-exons-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn put_gz(&self, rel: &str, text: &str) {
        let p = self.0.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        e.write_all(text.as_bytes()).unwrap();
        fs::write(&p, e.finish().unwrap()).unwrap();
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The data rows of a committed asset table, split on tabs, header first.
fn asset(name: &str) -> Vec<Vec<String>> {
    fs::read_to_string(format!("{REPO}/assets/{name}"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

/// The catalog genes with a position: (gene, chrom, start, end, owner).
fn genes() -> Vec<(String, String, i64, i64, String)> {
    asset("callability_genes.tsv")[1..]
        .iter()
        .filter(|r| r[3] != ".")
        .map(|r| {
            (
                r[0].clone(),
                r[2].clone(),
                r[3].parse().unwrap(),
                r[4].parse().unwrap(),
                r[6].clone(),
            )
        })
        .collect()
}

/// The committed exon rows: (exon, gene, chrom, start, end).
fn exons() -> Vec<(String, String, String, i64, i64)> {
    let rows = asset("callability_exons.tsv");
    assert_eq!(
        rows[0],
        ["exon", "gene", "transcript", "chrom", "start", "end"],
        "the exon table's header"
    );
    rows[1..]
        .iter()
        .map(|r| {
            (
                r[0].clone(),
                r[1].clone(),
                r[3].clone(),
                r[4].parse().unwrap(),
                r[5].parse().unwrap(),
            )
        })
        .collect()
}

/// The first catalog gene with `owner` ("." for none) and two or more
/// coding exons.
fn pick(owned: bool) -> String {
    let ex = exons();
    genes()
        .into_iter()
        .find(|g| (g.4 != ".") == owned && ex.iter().filter(|e| e.1 == g.0).count() >= 2)
        .expect("a catalog gene with two coding exons")
        .0
}

/// The share of a gene's body at 10x and of each of its exons.
struct Plan<'a> {
    body: &'a dyn Fn(&str) -> f64,
    exon: &'a dyn Fn(&str) -> f64,
}

const HEAD: &str = "#chrom\tstart\tend\tregion\t10X\n";

/// A sample folder with step 16b's gene table under `layout` (`mosdepth`
/// or `coverage`), and its exon table when `exon_table`.
fn folder(tag: &str, layout: &str, plan: &Plan, exon_table: bool) -> Dir {
    let d = Dir::new(tag);
    let bases = |n: i64, s: f64| ((n as f64 * s).ceil() as i64).min(n);
    let mut bed = String::from(HEAD);
    for (g, chrom, s, e, _) in genes() {
        bed.push_str(&format!(
            "{chrom}\t{s}\t{e}\t{g}\t{}\n",
            bases(e - s, (plan.body)(&g))
        ));
    }
    d.put_gz(&format!("{layout}/S_callability.bed.gz"), &bed);
    if exon_table {
        let mut bed = String::from(HEAD);
        for (x, _, chrom, s, e) in exons() {
            bed.push_str(&format!(
                "{chrom}\t{s}\t{e}\t{x}\t{}\n",
                bases(e - s, (plan.exon)(&x))
            ));
        }
        d.put_gz(&format!("{layout}/S_callability_exons.bed.gz"), &bed);
    }
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

fn get<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(key),
        _ => panic!("not an object: {v:?}"),
    }
}

fn at<'a>(v: &'a Value, key: &str) -> &'a Value {
    get(v, key).unwrap_or_else(|| panic!("no {key} in {v:?}"))
}

fn section(v: &Value) -> &Value {
    at(at(v, "sections"), "callability")
}

fn data(v: &Value) -> &Value {
    at(section(v), "data")
}

fn list<'a>(v: &'a Value, key: &str) -> &'a [Value] {
    match at(data(v), key) {
        Value::Array(a) => a,
        x => panic!("not a list: {x:?}"),
    }
}

fn st(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => panic!("not a string: {v:?}"),
    }
}

fn num(v: &Value) -> f64 {
    match v {
        Value::Float(f) => *f,
        Value::Int(i) => *i as f64,
        _ => panic!("not a number: {v:?}"),
    }
}

/// The genes of `key` (`not_read` or `not_read_owned`) on this reference:
/// the ALT-only genes are listed as not on it whatever the reads say.
fn unread<'a>(v: &'a Value, key: &str) -> Vec<&'a Value> {
    list(v, key)
        .iter()
        .filter(|g| st(at(g, "status")) != "not on this reference")
        .collect()
}

/// The entry of `gene` in `key` (`not_read` or `not_read_owned`).
fn entry<'a>(v: &'a Value, key: &str, gene: &str) -> Option<&'a Value> {
    list(v, key).iter().find(|g| st(at(g, "gene")) == gene)
}

/// The exon of `gene` set to 0%: its second coding exon.
fn second_exon(gene: &str) -> String {
    exons()
        .into_iter()
        .filter(|e| e.1 == gene)
        .nth(1)
        .unwrap()
        .0
}

#[test]
fn a_gene_read_on_its_body_with_an_unread_exon_reads_exon_not_read() {
    let gene = pick(false);
    let bad = second_exon(&gene);
    let plan = Plan {
        body: &|g| if g == gene { 0.95 } else { 1.0 },
        exon: &|x| if x == bad { 0.0 } else { 1.0 },
    };
    let v = collect(&folder("exon-0", "mosdepth", &plan, true));
    assert_eq!(st(at(section(&v), "state")), "ok");
    let g = entry(&v, "not_read", &gene)
        .unwrap_or_else(|| panic!("{gene} is not listed as not read: {:?}", data(&v)));
    assert_eq!(st(at(g, "status")), "exon not read");
    let body = num(at(g, "share"));
    assert!(
        (0.95..0.96).contains(&body),
        "the share stays the body's: {body}"
    );
    let w = at(g, "worst_exon");
    assert_eq!(
        st(at(w, "exon")),
        bad,
        "the exon that was not read is named"
    );
    assert_eq!(num(at(w, "share")), 0.0);
    assert_eq!(unread(&v, "not_read").len(), 1, "every other gene was read");
    let catalog = genes().len() as f64;
    assert_eq!(num(at(data(&v), "read")), catalog - 1.0);
    assert_eq!(
        num(at(data(&v), "exons_checked")),
        exons().len() as f64,
        "every exon row of the committed table belongs to a catalog gene"
    );
}

#[test]
fn control_a_gene_with_every_exon_read_is_read() {
    let gene = pick(false);
    let plan = Plan {
        body: &|g| if g == gene { 0.95 } else { 1.0 },
        exon: &|_| 0.9,
    };
    let v = collect(&folder("exon-90", "mosdepth", &plan, true));
    assert!(entry(&v, "not_read", &gene).is_none(), "{:?}", data(&v));
    assert!(unread(&v, "not_read").is_empty());
}

#[test]
fn a_gene_not_read_on_its_body_names_its_worst_exon_too() {
    let gene = pick(false);
    let bad = second_exon(&gene);
    let plan = Plan {
        body: &|g| if g == gene { 0.5 } else { 1.0 },
        exon: &|x| if x == bad { 0.0 } else { 0.5 },
    };
    let v = collect(&folder("body", "coverage", &plan, true));
    let g = entry(&v, "not_read", &gene).unwrap();
    assert_eq!(st(at(g, "status")), "not read");
    assert_eq!(st(at(at(g, "worst_exon"), "exon")), bad);
    assert_eq!(num(at(at(g, "worst_exon"), "share")), 0.0);
}

#[test]
fn an_owned_gene_with_an_unread_exon_reads_exon_not_read() {
    let gene = pick(true);
    let bad = second_exon(&gene);
    let plan = Plan {
        body: &|_| 1.0,
        exon: &|x| if x == bad { 0.5 } else { 1.0 },
    };
    // the Nextflow layout: both tables under coverage/
    let v = collect(&folder("owned", "coverage", &plan, true));
    let g = entry(&v, "not_read_owned", &gene).unwrap();
    assert_eq!(st(at(g, "status")), "exon not read");
    assert_eq!(st(at(at(g, "worst_exon"), "exon")), bad);
    assert!(unread(&v, "not_read").is_empty());
}

#[test]
fn without_the_exon_table_the_summary_is_as_before() {
    let gene = pick(false);
    let plan = Plan {
        body: &|g| if g == gene { 0.95 } else { 1.0 },
        exon: &|_| 0.0,
    };
    let v = collect(&folder("no-table", "mosdepth", &plan, false));
    assert!(unread(&v, "not_read").is_empty(), "{:?}", data(&v));
    assert!(get(data(&v), "exons_checked").is_none());
}

#[test]
fn a_row_not_named_gene_exon_n_makes_the_section_unreadable() {
    let d = folder(
        "bad-name",
        "mosdepth",
        &Plan {
            body: &|_| 1.0,
            exon: &|_| 1.0,
        },
        false,
    );
    d.put_gz(
        "mosdepth/S_callability_exons.bed.gz",
        &format!("{HEAD}chr1\t10\t20\tBRCA1\t10\n"),
    );
    let v = collect(&d);
    assert_eq!(st(at(section(&v), "state")), "unreadable");
    assert!(
        st(at(section(&v), "note"))
            .contains("S_callability_exons.bed.gz: \"BRCA1\" is not named GENE:exonN"),
        "{:?}",
        section(&v)
    );
}

#[test]
fn exon_rows_of_genes_outside_the_catalog_are_left_out() {
    let d = folder(
        "other",
        "mosdepth",
        &Plan {
            body: &|_| 1.0,
            exon: &|_| 1.0,
        },
        false,
    );
    d.put_gz(
        "mosdepth/S_callability_exons.bed.gz",
        &format!("{HEAD}chr1\t10\t20\tNOT_A_CATALOG_GENE:exon1\t0\n"),
    );
    let v = collect(&d);
    assert!(unread(&v, "not_read").is_empty());
    assert_eq!(num(at(data(&v), "exons_checked")), 0.0);
}
