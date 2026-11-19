//! The copy-aware markers (`og cn marker-regions`, `og cn marker-read
//! SMN1_27134`): the marker checks of the old `tests/test_traits_report.py`,
//! by meaning.
//!
//! 1. Each marker's pinned dbSNP record is one of its copies with its ref
//!    and alt; the controls are caught: ref and alt swapped, the minus-strand
//!    bases, a dbSNP position one base off the copies (the planted record of
//!    the old guard.yml control), and a context whose middle base is not the
//!    ref. KIV-2 4925G>A has seven copies.
//! 2. `og cn marker-regions` prints one region per copy, as written here.
//! 3. SMN1 g.27134T>G read from the pileup at both copies beside step 35's
//!    SMN table: present, not seen, uncertain, too little depth, not
//!    assessed; reads elsewhere are not counted; the SMN1 psCN at exon 7 (or
//!    through its homologous region) with its quality; and every message,
//!    for every status and psCN, says the marker does not prove 2+0.

use std::collections::HashMap;
use std::path::PathBuf;

use og_cn::markers::{LPA_KIV2_4925, MARKERS, Marker, SMN1_27134, marker_errors};

/// The copies of each marker, written here as the test's own copy.
const COPIES: [(&str, &[i64]); 2] = [
    (
        "LPA_KIV2_4925",
        &[
            160_646_214,
            160_640_667,
            160_635_123,
            160_629_574,
            160_624_030,
            160_618_484,
            160_612_932,
        ],
    ),
    ("SMN1_27134", &[70_952_074, 70_076_654]),
];

const SMN_HEADER: &str = "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n";
const NOT_PROOF: &str = "This marker does not prove 2+0, and its absence does not rule 2+0 out.";

fn og_cn(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_cn::cli::main(&args, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        x => x,
    }
}

/// The context with its middle base replaced.
fn with_middle(m: &Marker, base: u8) -> &'static str {
    let mut c = m.context.as_bytes().to_vec();
    let mid = c.len() / 2;
    c[mid] = base;
    Box::leak(String::from_utf8(c).unwrap().into_boxed_str())
}

#[test]
fn each_marker_agrees_with_its_dbsnp_record_and_the_controls_are_caught() {
    for m in &MARKERS {
        assert_eq!(marker_errors(m), Vec::<String>::new(), "{}", m.key);
        let mut swapped = *m;
        (swapped.reference, swapped.alt) = (m.alt, m.reference);
        assert!(
            marker_errors(&swapped)
                .iter()
                .any(|e| e.contains("swapped")),
            "{}: swap",
            m.key
        );
        let mut minus = *m;
        (minus.reference, minus.alt) = (complement(m.reference), complement(m.alt));
        minus.context = with_middle(m, complement(m.reference));
        assert!(
            marker_errors(&minus)
                .iter()
                .any(|e| e.contains("minus-strand")),
            "{}: minus strand",
            m.key
        );
        let mut off = *m;
        off.dbsnp.2 += 1;
        let e = marker_errors(&off);
        assert!(
            e.iter()
                .any(|e| e.contains("not at one of the marker's copies")),
            "{}: pos + 1: {e:?}",
            m.key
        );
        let mut ctx = *m;
        ctx.context = with_middle(m, m.alt);
        assert!(
            marker_errors(&ctx)
                .iter()
                .any(|e| e.contains("middle base")),
            "{}: context",
            m.key
        );
    }
    assert_eq!(LPA_KIV2_4925.copies.len(), 7, "six KIV-2 units and KIV-1");
}

#[test]
fn marker_regions_prints_one_region_per_copy() {
    for (key, copies) in COPIES {
        let m = MARKERS.iter().find(|m| m.key == key).unwrap();
        let (rc, out, err) = og_cn(&["marker-regions", key]);
        let want: Vec<String> = copies
            .iter()
            .map(|p| format!("{}:{p}-{p}", m.chrom))
            .collect();
        assert_eq!(
            (
                rc,
                out.split_whitespace()
                    .map(str::to_owned)
                    .collect::<Vec<_>>()
            ),
            (0, want),
            "{key}: {err}"
        );
    }
}

struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d =
            std::env::temp_dir().join(format!("og-genome-markers-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn put(&self, rel: &str, text: &str) -> String {
        let p = self.0.join(rel);
        std::fs::write(&p, text).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn path(&self, rel: &str) -> String {
        self.0.join(rel).to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A pileup with `t` T reads and `g` G reads at each SMN1 copy.
fn pileup(w: &Scratch, t: usize, g: usize, extra: &str) -> String {
    let text: String = SMN1_27134
        .copies
        .iter()
        .map(|p| {
            format!(
                "chr5\t{p}\tN\t{}\t^I{}{}$\t.\n",
                t + g,
                "T".repeat(t),
                "g".repeat(g)
            )
        })
        .collect();
    w.put("p.tsv", &(text + extra))
}

/// Step 35's SMN table with one row over SMN1's exon 7.
fn smn(w: &Scratch, pscn: &str, pqual: &str) -> String {
    w.put(
        "smn.tsv",
        &format!(
            "{SMN_HEADER}chr5\t70925000\t70953000\tSMN1\tPASS\t4\t50\tPASS\t{pscn}\t{pqual}\tchr5:70049000-70078000\n"
        ),
    )
}

/// marker-read SMN1_27134 on a pileup and an SMN table: exit and the table.
fn read_marker(w: &Scratch, p: &str, smn_path: &str) -> (u8, HashMap<String, String>) {
    let out = w.path("o.tsv");
    let _ = std::fs::remove_file(&out);
    let (rc, _, _) = og_cn(&[
        "marker-read",
        "SMN1_27134",
        "--pileup",
        p,
        "--copy-table",
        smn_path,
        "--out",
        &out,
    ]);
    let t = std::fs::read_to_string(&out)
        .unwrap_or_default()
        .lines()
        .skip(1)
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, v)| (k.to_owned(), v.to_owned()))
        .collect();
    (rc, t)
}

fn g<'a>(t: &'a HashMap<String, String>, k: &str) -> &'a str {
    t.get(k).map_or("", String::as_str)
}

#[test]
fn smn1_marker_statuses_and_the_smn1_pscn() {
    let w = Scratch::new("smn1");
    let s = smn(&w, "2,2", "30,30");
    let (rc, t) = read_marker(&w, &pileup(&w, 15, 5, ""), &s);
    assert_eq!(
        (
            rc,
            g(&t, "status"),
            g(&t, "alt_reads"),
            g(&t, "depth"),
            g(&t, "fraction")
        ),
        (0, "present", "10", "40", "0.2500")
    );
    assert_eq!(
        (g(&t, "copies"), g(&t, "agcn")),
        ("1.0", "4"),
        "copies = fraction x agCN"
    );
    assert_eq!(g(&t, "smn1_pscn"), "2 (quality 30, filter PASS)");
    assert!(
        g(&t, "message").contains("higher chance")
            && g(&t, "message").contains("does not prove 2+0"),
        "{t:?}"
    );
    let (_, t) = read_marker(&w, &pileup(&w, 20, 0, ""), &s);
    assert!(
        g(&t, "status") == "not seen"
            && g(&t, "copies") == "."
            && g(&t, "message").contains("lower"),
        "{t:?}"
    );
    let (_, t) = read_marker(&w, &pileup(&w, 20, 0, "chr5\t70952074\tN\t2\tGG\t.\n"), &s);
    assert_eq!(g(&t, "status"), "uncertain");
    let (_, t) = read_marker(&w, &pileup(&w, 3, 1, ""), &s);
    assert_eq!(g(&t, "status"), "too little depth");
    let other = format!("chr5\t70952075\tN\t30\t{}\t.\n", "G".repeat(30));
    let (_, t) = read_marker(&w, &pileup(&w, 15, 5, &other), &s);
    assert_eq!(
        g(&t, "alt_reads"),
        "10",
        "reads at another position are not counted"
    );
    let bad = w.put("bad.tsv", "nonsense\n");
    let (rc, t) = read_marker(&w, &bad, &s);
    assert_eq!(
        (rc, g(&t, "status")),
        (0, "not assessed"),
        "an unreadable pileup"
    );
    let (_, t) = read_marker(&w, &pileup(&w, 15, 5, ""), &smn(&w, "1,3", "30,30"));
    assert!(g(&t, "message").contains("SMN1 has 1 copy here"), "{t:?}");
    let (_, t) = read_marker(&w, &pileup(&w, 15, 5, ""), &smn(&w, "2,2", "12,30"));
    assert_eq!(
        g(&t, "smn1_pscn"),
        "2 (quality 12, filter PASS), not resolved"
    );
    assert!(
        g(&t, "message").contains("SMN1 copy number: 2 (quality 12, filter PASS), not resolved;"),
        "{t:?}"
    );
    // SMN2's region first: exon 7 of SMN1 is its homologous region.
    let s2 = w.put(
        "smn2.tsv",
        &format!("{SMN_HEADER}chr5\t70049000\t70078000\tSMN1\tPASS\t4\t50\tPASS\t2,2\t30,40\tchr5:70925000-70953000\n"),
    );
    let (_, t) = read_marker(&w, &pileup(&w, 15, 5, ""), &s2);
    assert_eq!(
        (g(&t, "smn1_pscn"), g(&t, "agcn")),
        ("2 (quality 40, filter PASS)", "4")
    );
    let (_, t) = read_marker(&w, &pileup(&w, 15, 5, ""), &w.path("missing.tsv"));
    assert!(g(&t, "smn1_pscn").contains("no SMN table"), "{t:?}");
    let dir = w.0.to_string_lossy().into_owned();
    let (rc, t) = read_marker(&w, &pileup(&w, 15, 5, ""), &dir);
    assert_eq!((rc, g(&t, "status"), g(&t, "agcn")), (0, "present", "."));
    assert!(g(&t, "smn1_pscn").contains("could not be read"), "{t:?}");
}

#[test]
fn every_smn1_marker_message_says_it_does_not_prove_2_plus_0() {
    let w = Scratch::new("smn1-msg");
    let bad = w.put("bad.tsv", "nonsense\n");
    let mut msgs = Vec::new();
    for counts in [
        Some((30, 10)),
        Some((40, 0)),
        Some((38, 2)),
        Some((3, 1)),
        None,
    ] {
        let p = match counts {
            None => bad.clone(),
            Some((t, g)) => w.put(
                "all.tsv",
                &format!(
                    "chr5\t{}\tN\t{}\t{}{}\t.\n",
                    SMN1_27134.copies[0],
                    t + g,
                    "T".repeat(t),
                    "G".repeat(g)
                ),
            ),
        };
        for path in [Some("2,2"), Some("1,3"), Some("5,5"), None] {
            let s = match path {
                None => w.path("missing.tsv"),
                Some("5,5") => smn(&w, "2,2", "5,5"),
                Some(ps) => smn(&w, ps, "30,30"),
            };
            msgs.push(g(&read_marker(&w, &p, &s).1, "message").to_owned());
        }
    }
    assert_eq!(msgs.len(), 20);
    let without: Vec<&String> = msgs
        .iter()
        .filter(|m| !m.contains("does not prove 2+0"))
        .collect();
    assert!(without.is_empty(), "{without:?}");
    assert!(msgs[0].contains(NOT_PROOF), "{}", msgs[0]);
    assert!(
        !msgs[0]
            .replace(NOT_PROOF, "")
            .contains("does not prove 2+0"),
        "control: without the sentence"
    );
}
