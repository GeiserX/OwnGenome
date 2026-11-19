//! `og mei report` (step 50): what makes a known disease insertion site
//! "seen" (#411 sv-33, sv-34).
//!
//! MEGAnE writes every insertion it identifies, with its own filter result
//! in FILTER: `PASS`, or flags such as `LC` (low confidence). A call that
//! failed MEGAnE's filters must not make a catalog site "seen"; the site
//! names the filtered call instead. Without a call, a window with too few
//! reads (no bedcov line, or a mean under 10x, zero included) reads
//! "not callable", never "not seen".
//!
//! Each case writes a synthetic MEGAnE VCF, a two-row catalog (the FKTN SVA
//! window and the TAF1 SVA window) and `samtools bedcov -c` lines, then
//! reads the sites table `report` writes.

use og_reports::mei::{Files, report};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

const HEAD: &str = "##fileformat=VCFv4.1
##ALT=<ID=INS:ME,Description=\"Insertion of mobile element\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
";

const SITES: &str = "gene\tchrom\tstart\tend\tsite\telement\tinheritance\tcondition\tclinvar\tpmid
FKTN\tchr9\t105639156\t105640156\t105639656\tSVA\tautosomal recessive\tFukuyama congenital muscular dystrophy\t3199\t9690476
TAF1\tchrX\t71439993\t71440993\t71440493\tSVA\tX-linked\tX-linked dystonia-parkinsonism\t9862\t29474918
";

/// 30x over both windows (30030 bases summed over 1001).
const DEPTH_30X: &str =
    "chr9\t105639155\t105640156\tFKTN\t30030\t200\nchrX\t71439992\t71440993\tTAF1\t30030\t200\n";

/// An SVA call MEGAnE identified (MEPRED=PASS), one copy, with this FILTER.
fn sva(chrom: &str, pos: u64, filter: &str) -> String {
    format!(
        "{chrom}\t{pos}\t.\tN\t<INS:ME>\t.\t{filter}\tSVTYPE=Retroposon/SVA;MEPRED=PASS;MEI_rpos={};MEI=SVA_F,1/1500,+/+;CN_conf=H\tCN\t1\n",
        pos + 8
    )
}

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!("og-i411-mei-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn p(&self, n: &str) -> String {
        self.0.join(n).to_str().unwrap().to_owned()
    }
    fn put(&self, n: &str, t: &str) -> String {
        fs::write(self.p(n), t).unwrap();
        self.p(n)
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

type Row = HashMap<String, String>;

/// The sites table, one map per catalog row, keyed by gene.
fn sites(tag: &str, records: &[String], depth: &str) -> HashMap<String, Row> {
    let d = Dir::new(tag);
    let f = Files {
        vcf: &d.put("m.vcf", &format!("{HEAD}{}", records.concat())),
        sites: &d.put("sites.tsv", SITES),
        depth: &d.put("depth.tsv", depth),
        annotsv: None,
        out_calls: &d.p("calls.tsv"),
        out_sites: &d.p("sites_out.tsv"),
    };
    report("S", &f).unwrap();
    let text = fs::read_to_string(d.p("sites_out.tsv")).unwrap();
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|l| {
            let row: Row = head
                .iter()
                .map(|h| (*h).to_owned())
                .zip(l.split('\t').map(str::to_owned))
                .collect();
            (row["gene"].clone(), row)
        })
        .collect()
}

fn status(rows: &HashMap<String, Row>, gene: &str) -> (String, String) {
    (rows[gene]["status"].clone(), rows[gene]["call"].clone())
}

/// Red first: an SVA call at the FKTN site that MEGAnE flagged LC is not
/// "seen". The site says the call was filtered out, and why, and still
/// names the call.
#[test]
fn red_first_a_filtered_call_is_not_seen() {
    let rows = sites("lc", &[sva("chr9", 105_639_660, "LC")], DEPTH_30X);
    assert_eq!(
        status(&rows, "FKTN"),
        (
            "call filtered out (LC)".to_owned(),
            "chr9:105639660 SVA_F 1 copy LC".to_owned()
        )
    );
    assert_eq!(status(&rows, "TAF1").0, "not seen");
}

/// Red first: a filtered call at a thin window still names the call; it
/// does not hide behind "not callable".
#[test]
fn red_first_a_filtered_call_at_low_depth() {
    let thin = "chr9\t105639155\t105640156\tFKTN\t5005\t30\n";
    let rows = sites("lc-thin", &[sva("chr9", 105_639_660, "LC")], thin);
    assert_eq!(status(&rows, "FKTN").0, "call filtered out (LC)");
    assert_eq!(rows["FKTN"]["mean_depth"], "5.0");
}

/// A PASS call wins over a filtered one in the same window, whichever
/// comes first.
#[test]
fn a_pass_call_in_the_window_is_seen() {
    let rows = sites(
        "pass-and-lc",
        &[
            sva("chr9", 105_639_400, "LC"),
            sva("chr9", 105_639_660, "PASS"),
        ],
        DEPTH_30X,
    );
    assert_eq!(
        status(&rows, "FKTN"),
        (
            "seen".to_owned(),
            "chr9:105639660 SVA_F 1 copy PASS".to_owned()
        )
    );
    // A PASS call at low depth is still seen: the call is the evidence.
    let thin = "chr9\t105639155\t105640156\tFKTN\t5005\t30\n";
    let rows = sites("pass-thin", &[sva("chr9", 105_639_660, "PASS")], thin);
    assert_eq!(status(&rows, "FKTN").0, "seen");
}

/// Without a call, zero depth, or no bedcov line at all, is "not callable"
/// with no call, never "not seen".
#[test]
fn zero_depth_is_not_callable() {
    let zero = "chr9\t105639155\t105640156\tFKTN\t0\t0\nchrX\t71439992\t71440993\tTAF1\t0\t0\n";
    let rows = sites("zero", &[], zero);
    for gene in ["FKTN", "TAF1"] {
        assert_eq!(
            [
                rows[gene]["status"].as_str(),
                &rows[gene]["mean_depth"],
                &rows[gene]["call"]
            ],
            ["not callable", "0.0", "."],
            "{gene}"
        );
    }
    let rows = sites("no-line", &[], "");
    assert_eq!(status(&rows, "FKTN").0, "not callable");
    assert_eq!(rows["FKTN"]["mean_depth"], ".");
    // Only a filtered call far from the window: depth decides, as without
    // any call.
    let rows = sites("lc-far", &[sva("chr9", 105_700_000, "LC")], zero);
    assert_eq!(status(&rows, "FKTN").0, "not callable");
}
