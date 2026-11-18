//! Step 42's traits table and genotype rules, from `og traits table` through
//! og summary collect to both reports. These are the checks of
//! `tests/test_collect_summary.py` (section 18 and the rules part of section
//! 26) and `tests/test_render_report.py` (section 11) that ran the traits
//! table through the helper of the old `tests/test_traits_report.py`; they
//! moved here with it, and now check the Rust reports.

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use og_render::notes::TRAITS_NOTE;
use og_render::sections::apoe::APOE_WARNING;
use og_summary::collect::{Group, Input};
use pycompat::json::{self, Dump, Value};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");

fn repo(rel: &str) -> String {
    format!("{REPO}/{rel}")
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_default()
}

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!(
            "og-genome-traits-report-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn path(&self, rel: &str) -> String {
        let p = self.0.join(rel);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn put(&self, rel: &str, text: &str) -> String {
        let p = self.path(rel);
        std::fs::write(&p, text).unwrap();
        p
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// `og traits table` on these genotype lines, writing the sample's traits
/// table into `d/traits/` and the rules table to `rules_out`.
fn traits_tables(w: &Scratch, genotypes: &str, d: &str, rules_out: &str) {
    let g = w.put("g.tsv", genotypes);
    std::fs::create_dir_all(format!("{d}/traits")).unwrap();
    let args: Vec<String> = [
        "table",
        "--table",
        &repo("assets/traits.tsv"),
        "--rules",
        &repo("assets/genotype_rules.tsv"),
        "--genotypes",
        &g,
        "--sample",
        "S",
        "--out",
        &format!("{d}/traits/S_traits.tsv"),
        "--rules-out",
        rules_out,
    ]
    .iter()
    .map(|a| (*a).to_owned())
    .collect();
    assert_eq!(og_steps::traits::cli(&args), 0);
}

fn summary(d: &str) -> Value {
    og_summary::collect::post::collect(&Input {
        sample: "S".into(),
        sample_dir: d.into(),
        declared_sex: None,
        assets: repo("assets"),
        group: Group::All,
        tables: repo("bin"),
    })
    .unwrap()
}

/// (summary, text report, HTML report) of the sample in `d`.
fn render(d: &str) -> (Value, String, String) {
    let v = summary(d);
    let text = json::dumps(&v, &Dump::default()).unwrap();
    let sm = og_render::model::Summary::parse(&text).unwrap();
    (
        v.clone(),
        og_render::text::report(&sm).unwrap(),
        og_render::html::report(&sm).unwrap(),
    )
}

fn at<'a>(v: &'a Value, path: &[&str]) -> Option<&'a Value> {
    path.iter().try_fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k),
        _ => None,
    })
}

fn s(v: Option<&Value>) -> &str {
    match v {
        Some(Value::Str(s)) => s,
        _ => "",
    }
}

fn list(v: Option<&Value>) -> &[Value] {
    match v {
        Some(Value::Array(a)) => a,
        _ => &[],
    }
}

fn int(v: Option<&Value>) -> i128 {
    match v {
        Some(Value::Int(i)) => *i,
        other => panic!("{other:?}"),
    }
}

/// The HTML with the escapes the report writes undone.
fn unescape(h: &str) -> String {
    h.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&#x27;", "'")
        .replace("&amp;", "&")
}

/// The Traits card of an HTML report.
fn card(html: &str) -> String {
    html.split_once("<h2>Traits</h2>")
        .map_or(String::new(), |x| {
            x.1.split("<h2>").next().unwrap().to_owned()
        })
}

/// The rows of assets/traits.tsv, as maps.
fn traits_rows() -> Vec<HashMap<String, String>> {
    let text = read(&repo("assets/traits.tsv"));
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            head.iter()
                .map(|h| (*h).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect()
        })
        .collect()
}

/// og dev schema on the summary: true when it validates.
fn validates(w: &Scratch, v: &Value) -> bool {
    let p = w.put("summary.json", &json::dumps(v, &Dump::default()).unwrap());
    Command::new(env!("CARGO_BIN_EXE_og"))
        .args([
            "dev",
            "schema",
            &repo("tests/schema/summary.schema.json"),
            &p,
        ])
        .status()
        .unwrap()
        .success()
}

#[test]
fn the_genotype_rules_reach_both_reports() {
    let w = Scratch::new("rules");
    let t = w.path("S/x").replace("/x", "");
    // TAS2R38 homozygous for the taster haplotype PAV (rs1726866's G is the reference base).
    let g = "chr7\t141972804\tT\tC\tC/C\t40\nchr7\t141972905\tG\t.\tG/G\t40\nchr7\t141973545\tC\tG\tG/G\t40\n\
             chr19\t44908684\tT\t.\tT/T\t40\nchr19\t44908822\tC\t.\tC/C\t40\n";
    let saved = w.path("rules_saved.tsv");
    traits_tables(&w, g, &t, &saved);
    let (v, txt, html) = render(&t);
    assert!(
        at(&v, &["sections", "traits", "data", "rules"]).is_none_or(|r| *r == Value::Null),
        "no rules table"
    );
    assert!(!(txt.clone() + &html).contains("Genotype rules (several sites"));
    std::fs::copy(&saved, format!("{t}/traits/S_genotype_rules.tsv")).unwrap();
    let (v, txt, html) = render(&t);
    let html = unescape(&html);
    let rl: Vec<[&str; 4]> = list(at(&v, &["sections", "traits", "data", "rules"]))
        .iter()
        .filter(|r| ["TAS2R38", "APOE"].contains(&s(at(r, &["rule"]))))
        .map(|r| ["rule", "result", "status", "phase"].map(|k| s(at(r, &[k]))))
        .collect();
    assert_eq!(
        rl,
        [["TAS2R38", "PAV/PAV", "called", "not_needed"]],
        "TAS2R38 read, APOE left out"
    );
    let text = read(&saved);
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let trait_ = lines
        .map(|l| {
            head.iter()
                .copied()
                .zip(l.split('\t'))
                .collect::<HashMap<_, _>>()
        })
        .find(|r| r["rule"] == "TAS2R38")
        .unwrap()["trait"]
        .to_owned();
    assert!(
        txt.contains(&format!(
            "{trait_} (TAS2R38 rs713598, rs1726866, rs10246939): PAV/PAV (no phase needed). "
        )),
        "{txt}"
    );
    assert!(txt.contains("traits/S_genotype_rules.tsv"));
    let c = card(&html);
    assert!(
        c.contains("<h3>Genotype rules (several sites read together)</h3>"),
        "{c}"
    );
    assert!(
        c.contains("<td>G/G, G/G, C/C</td><td>PAV/PAV (no phase needed)</td>"),
        "{c}"
    );
    let folded = c.split_once("<details>").map_or("", |x| x.1);
    let rules_part = c
        .split_once("<h3>Genotype rules")
        .map_or("", |x| x.1)
        .split("<details>")
        .next()
        .unwrap();
    assert!(
        folded.contains("rs429358+rs7412") && !rules_part.contains("APOE"),
        "the APOE lines stay folded away"
    );
    let rows = traits_rows();
    let tas: Vec<&HashMap<String, String>> =
        rows.iter().filter(|r| r["gene"] == "TAS2R38").collect();
    assert!(
        tas.iter().all(
            |r| r["note"].contains("genotype rule") && !r["note"].contains("is not known here")
        )
    );
    assert!(validates(&w, &v), "the summary validates");
    w.put(
        "S/traits/S_genotype_rules.tsv",
        "rule\tgene\tresult\nTAS2R38\tTAS2R38\tPAV/PAV\n",
    );
    assert_eq!(
        s(at(&summary(&t), &["sections", "traits", "state"])),
        "unreadable",
        "no reading column"
    );
}

#[test]
fn the_traits_section_and_both_reports() {
    let w = Scratch::new("panel");
    let t = w.path("S/x").replace("/x", "");
    let rows = traits_rows();
    let by: HashMap<&str, &HashMap<String, String>> =
        rows.iter().map(|r| (r["rsid"].as_str(), r)).collect();
    let (lct, ear) = (by["rs4988235"], by["rs17822931"]);
    let g = format!(
        "{}\t{}\t{}\t{}\t{}/{}\n{}\t{}\t{}\t.\t{}/{}\nchr19\t44908684\tT\t.\tT/T\nchr19\t44908822\tC\t.\tC/C\n",
        lct["chrom"],
        lct["pos"],
        lct["ref"],
        lct["alt"],
        lct["ref"],
        lct["alt"],
        ear["chrom"],
        ear["pos"],
        ear["ref"],
        ear["ref"],
        ear["ref"]
    );
    traits_tables(&w, &g, &t, &w.path("traits_rules.tsv"));
    let (v, txt, html) = render(&t);
    let td = at(&v, &["sections", "traits"]).unwrap();
    assert_eq!(s(at(td, &["state"])), "ok");
    assert_eq!(
        (
            int(at(td, &["data", "called"])),
            int(at(td, &["data", "not_covered"]))
        ),
        (2, rows.len() as i128 - 4),
        "counts without the APOE rows"
    );
    assert_eq!(list(at(td, &["data", "apoe"])).len(), 3);
    assert!(
        !list(at(td, &["data", "rows"]))
            .iter()
            .any(|r| s(at(r, &["gene"])) == "APOE")
    );
    let read_row = format!(
        "(MCM6 rs4988235): {}/{} [ancestry]. {} (PMID {})",
        lct["ref"], lct["alt"], lct["het"], lct["pmids"]
    );
    assert!(txt.contains(&read_row), "{read_row}");
    let not_read = "Not read (no genotype in the gVCF): ";
    assert!(txt.contains(not_read) && txt.split_once(not_read).unwrap().1.contains("rs671"));
    let apoe_at = txt.find(APOE_WARNING).expect("the APOE warning");
    assert!(
        txt[apoe_at..].contains("(rs429358+rs7412): e3/e3.") && !txt[..apoe_at].contains("e3/e3")
    );
    let c = card(&html);
    assert!(
        c.contains(&format!(
            "rs4988235</span></td><td>{}/{} ",
            lct["ref"], lct["alt"]
        )),
        "{c}"
    );
    assert!(c.contains("badge-yellow\">ancestry<"), "{c}");
    let (open, folded) = c.split_once("<details>").expect("a folded block");
    assert!(folded.contains("rs429358+rs7412") && !open.contains("rs429358"));
    assert!(validates(&w, &v));
    let tab = format!("{t}/traits/S_traits.tsv");
    std::fs::write(&tab, "rsid\tgene\tstatus\nrs1\tX\tcalled\n").unwrap();
    assert_eq!(
        s(at(&summary(&t), &["sections", "traits", "state"])),
        "unreadable",
        "no reading column"
    );
    std::fs::write(&tab, "rsid\tgene\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids\nrs1\tX\tT\tA/A\tcalled\n").unwrap();
    assert_eq!(
        s(at(&summary(&t), &["sections", "traits", "state"])),
        "unreadable",
        "a row cut short"
    );
    std::fs::remove_file(&tab).unwrap();
    let (v, _, html) = render(&t);
    assert_eq!(s(at(&v, &["sections", "traits", "state"])), "missing");
    assert!(
        list(at(&v, &["not_run"]))
            .iter()
            .any(|x| s(Some(x)) == "Traits panel")
    );
    assert!(card(&html).contains(">Not run<"));
}

#[test]
fn the_traits_card_note_no_reads_and_a_stale_table() {
    let w = Scratch::new("card");
    let d = w.path("S/x").replace("/x", "");
    traits_tables(&w, "", &d, &w.path("traits_rules.tsv"));
    let (v, txt, html) = render(&d);
    assert!(
        txt.contains(TRAITS_NOTE) && html.contains(TRAITS_NOTE),
        "both reports say a reading is not a diagnosis"
    );
    assert_eq!(int(at(&v, &["sections", "traits", "data", "called"])), 0);
    let apoe = list(at(&v, &["sections", "traits", "data", "apoe"]));
    assert_eq!(
        s(apoe.last().and_then(|r| at(r, &["status"]))),
        "not_covered"
    );
    let tab = format!("{d}/traits/S_traits.tsv");
    let now = SystemTime::now();
    let old = now - Duration::from_secs(86_400);
    std::fs::File::options()
        .write(true)
        .open(&tab)
        .unwrap()
        .set_modified(old)
        .unwrap();
    let started = now.duration_since(UNIX_EPOCH).unwrap().as_secs_f64() - 3600.0;
    w.put(
        "S/logs/run_status.tsv",
        &format!("meta\tstarted_epoch\t{started}\nstep\t42\tskipped (no gVCF: S/vcf/S.g.vcf.gz)\n"),
    );
    let (v, txt, html) = render(&d);
    assert_eq!(s(at(&v, &["sections", "traits", "state"])), "stale");
    assert!(txt.contains("[STALE: from an earlier run"));
    assert!(card(&html).contains("class=\"stale\">STALE: from an earlier run"));
    std::fs::File::options()
        .write(true)
        .open(&tab)
        .unwrap()
        .set_modified(SystemTime::now())
        .unwrap();
    assert_eq!(
        s(at(&summary(&d), &["sections", "traits", "state"])),
        "ok",
        "newer than the run"
    );
}
