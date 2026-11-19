//! Step 09's repeat counts as `og summary collect` reads them.
//!
//! 1. HTT 36 to 39 CAG is the reduced-penetrance range (GeneReviews:
//!    Huntington Disease). Stranger grades it with the catalog's pathogenic
//!    minimum alone, so 37 reads as plain `full_mutation` under STRchive's
//!    thresholds. A sized HTT allele in that range carries the range, its
//!    label and its citation from `assets/str_penetrance_ranges.tsv`.
//! 2. An allele whose ExpansionHunter estimate rests on in-repeat reads only
//!    (ADSP and ADFL 0, ADIR above 0) is longer than a read: it is
//!    `in_repeat`, not sized, whatever its SO field says. Its point estimate
//!    must never earn a range label.

use og_summary::collect::{Group, Input, collect};
use pycompat::json::Value;
use std::fs;
use std::path::PathBuf;

const HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
const FMT: &str = "GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC";

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let p = std::env::temp_dir().join(format!(
            "og-summary-i411-str-{tag}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = fs::remove_dir_all(&p);
        fs::create_dir_all(&p).unwrap();
        Dir(p)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).display().to_string()
    }

    fn put(&self, name: &str, body: &str) {
        let p = self.0.join(name);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(p, body).unwrap();
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// One ExpansionHunter record with Stranger's INFO fields. The counts are
/// invented.
fn rec(locus: &str, so: &str, cn: &str, ci: &str, ad: [&str; 3], st: &str) -> String {
    format!(
        "chr4\t3074877\t.\tC\t<STR40>\t.\tPASS\tEND=3074940;REPID={locus};VARID={locus};\
         STR_STATUS={st};STR_NORMAL_MAX=26;STR_PATHOLOGIC_MIN=36\t{FMT}\t\
         0/1:{so}:{cn}:{ci}:{}:{}:{}:30.0\n",
        ad[0], ad[1], ad[2]
    )
}

/// The expansions section's loci by name, for one ExpansionHunter VCF that
/// Stranger graded.
fn loci(tag: &str, records: &str) -> Vec<(String, Value)> {
    let d = Dir::new(tag);
    let vcf = format!("{HEAD}{records}");
    d.put("S/expansion_hunter/S_eh.vcf", &vcf);
    d.put("S/expansion_hunter/S_eh_stranger.vcf", &vcf);
    let input = Input {
        sample: "S".into(),
        sample_dir: d.path("S"),
        declared_sex: None,
        assets: concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets").into(),
        group: Group::C,
        tables: concat!(env!("CARGO_MANIFEST_DIR"), "/../../../bin").into(),
    };
    let v = collect(&input).unwrap().to_json();
    let sec = at(&v, &["sections", "expansions"]);
    assert_eq!(at(sec, &["state"]), &Value::Str("ok".into()), "{sec:?}");
    match at(sec, &["data", "loci"]) {
        Value::Array(xs) => xs
            .iter()
            .map(|x| (text(at(x, &["locus"])).to_owned(), x.clone()))
            .collect(),
        other => panic!("loci: {other:?}"),
    }
}

fn locus<'a>(ls: &'a [(String, Value)], name: &str) -> &'a Value {
    &ls.iter()
        .find(|(k, _)| k == name)
        .unwrap_or_else(|| panic!("no {name}"))
        .1
}

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no {k} in {v:?}")),
        _ => panic!("{k} of {v:?}"),
    })
}

fn get<'a>(v: &'a Value, k: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(k),
        _ => None,
    }
}

fn text(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => panic!("not a text: {v:?}"),
    }
}

const SPAN2: &str = "SPANNING/SPANNING";

#[test]
fn htt_37_reads_reduced_penetrance() {
    let ls = loci(
        "htt37",
        &rec(
            "HTT",
            SPAN2,
            "17/37",
            "17-17/37-37",
            ["12/9", "0/0", "0/0"],
            "normal,full_mutation",
        ),
    );
    let htt = locus(&ls, "HTT");
    assert_eq!(text(at(htt, &["sizing"])), "sized", "{htt:?}");
    assert_eq!(text(at(htt, &["penetrance"])), "36-39", "{htt:?}");
    assert!(
        text(at(htt, &["penetrance_label"])).starts_with("reduced penetrance"),
        "{htt:?}"
    );
    assert!(
        text(at(htt, &["penetrance_citation"])).contains("NBK1305"),
        "{htt:?}"
    );
    // Stranger's own grade stays as it wrote it.
    assert_eq!(text(at(htt, &["status"])), "normal,full_mutation");
}

#[test]
fn the_flagged_list_carries_the_penetrance_range() {
    let d = Dir::new("flagged");
    let vcf = format!(
        "{HEAD}{}",
        rec(
            "HTT",
            SPAN2,
            "20/38",
            "20-20/38-38",
            ["10/8", "0/0", "0/0"],
            "normal,full_mutation"
        )
    );
    d.put("S/expansion_hunter/S_eh.vcf", &vcf);
    d.put("S/expansion_hunter/S_eh_stranger.vcf", &vcf);
    let input = Input {
        sample: "S".into(),
        sample_dir: d.path("S"),
        declared_sex: None,
        assets: concat!(env!("CARGO_MANIFEST_DIR"), "/../../../assets").into(),
        group: Group::C,
        tables: concat!(env!("CARGO_MANIFEST_DIR"), "/../../../bin").into(),
    };
    let v = collect(&input).unwrap().to_json();
    let flagged = match at(&v, &["sections", "expansions", "data", "flagged"]) {
        Value::Array(xs) => xs.clone(),
        other => panic!("{other:?}"),
    };
    assert_eq!(flagged.len(), 1, "{flagged:?}");
    assert_eq!(text(at(&flagged[0], &["penetrance"])), "36-39");
}

#[test]
fn the_ends_of_the_range_and_their_neighbours() {
    for (cn, want) in [
        ("20/35", None),
        ("20/36", Some("36-39")),
        ("20/39", Some("36-39")),
        ("20/40", None),
    ] {
        let (a, b) = cn.split_once('/').unwrap();
        let ls = loci(
            &format!("ends-{b}"),
            &rec(
                "HTT",
                SPAN2,
                cn,
                &format!("{a}-{a}/{b}-{b}"),
                ["10/8", "0/0", "0/0"],
                "normal",
            ),
        );
        let htt = locus(&ls, "HTT");
        assert_eq!(get(htt, "penetrance").map(text), want, "{cn}: {htt:?}");
        // The keys come together or not at all.
        assert_eq!(get(htt, "penetrance_label").is_some(), want.is_some());
        assert_eq!(get(htt, "penetrance_citation").is_some(), want.is_some());
    }
}

#[test]
fn only_a_sized_allele_earns_the_label_and_only_at_its_locus() {
    // A low-confidence 37 (one spanning read), and 37 at another locus.
    let ls = loci(
        "unsized37",
        &(rec(
            "HTT",
            SPAN2,
            "17/37",
            "17-17/37-37",
            ["10/1", "0/0", "0/0"],
            "normal,full_mutation",
        ) + &rec(
            "ATXN1",
            SPAN2,
            "17/37",
            "17-17/37-37",
            ["10/8", "0/0", "0/0"],
            "normal",
        )),
    );
    let htt = locus(&ls, "HTT");
    assert_eq!(text(at(htt, &["sizing"])), "low_confidence");
    assert!(get(htt, "penetrance").is_none(), "{htt:?}");
    assert!(get(locus(&ls, "ATXN1"), "penetrance").is_none());
}

#[test]
fn an_allele_of_in_repeat_reads_only_is_not_sized() {
    // No SO value, and SO=FLANKING, on an allele that only in-repeat reads
    // support: each is longer than a read. The first allele has spanning
    // reads (sized, or unknown where SO is empty).
    for (so, want) in [
        (".", "unknown/in_repeat"),
        ("SPANNING/.", "sized/in_repeat"),
        ("SPANNING/FLANKING", "sized/in_repeat"),
    ] {
        let ls = loci(
            &format!("irr-{}", so.replace(['/', '.'], "_")),
            &rec(
                "HTT",
                so,
                "17/37",
                "17-17/37-60",
                ["12/0", "0/0", "0/14"],
                "normal,full_mutation",
            ),
        );
        let htt = locus(&ls, "HTT");
        assert_eq!(text(at(htt, &["allele_sizing"])), want, "{so}: {htt:?}");
        assert_eq!(text(at(htt, &["sizing"])), "in_repeat", "{so}");
        // The point estimate 37 of an unsized allele earns no label.
        assert!(get(htt, "penetrance").is_none(), "{so}: {htt:?}");
        assert_eq!(text(at(htt, &["zone"])), "", "{so}");
    }
}

#[test]
fn spanning_or_flanking_reads_keep_the_allele_as_it_was() {
    // Reads beside the in-repeat ones: the IRR rule does not fire.
    let ls = loci(
        "not-irr",
        &(rec(
            "HTT",
            ".",
            "17/37",
            "17-17/37-37",
            ["12/3", "0/0", "0/14"],
            "normal,full_mutation",
        ) + &rec(
            "ATXN1",
            "SPANNING/FLANKING",
            "17/50",
            "17-17/40-60",
            ["12/0", "0/4", "0/14"],
            "normal",
        ) + &rec(
            "DMPK",
            SPAN2,
            "5/12",
            "5-5/12-12",
            ["12/9", "0/0", "0/0"],
            "normal",
        )),
    );
    assert_eq!(
        text(at(locus(&ls, "HTT"), &["allele_sizing"])),
        "unknown/unknown"
    );
    assert_eq!(
        text(at(locus(&ls, "ATXN1"), &["allele_sizing"])),
        "sized/flanking_only"
    );
    assert_eq!(
        text(at(locus(&ls, "DMPK"), &["allele_sizing"])),
        "sized/sized"
    );
}
