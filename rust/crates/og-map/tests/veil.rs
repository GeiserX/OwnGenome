//! The veil tags on a synthetic person (`tests/data/facts`): an act finding
//! with a position, a second one beside it, a gene-only act finding, and
//! the APOE result an adult's page keeps behind its warning, with rows of a
//! source layer over and beside each. Every row over a closable finding
//! carries its id, no row beside one does, and no bin counts a tagged row.

use og_map::facts::{self, Findings};
use og_map::findings::{self, Placed};
use og_map::reference::{Reference, chrom_index};
use og_map::veil::{self, Span};
use std::fs;
use std::path::{Path, PathBuf};

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/facts")
}

fn assets() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets")
}

/// A GENOME_DIR with the synthetic run and the public GENCODE slice.
pub fn genome_dir(name: &str) -> PathBuf {
    let root = std::env::temp_dir().join(format!("og-map-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(root.join("DEMO-005")).unwrap();
    fs::create_dir_all(root.join("reference")).unwrap();
    fs::copy(
        data().join("DEMO-005/DEMO-005_sex_check.tsv"),
        root.join("DEMO-005/DEMO-005_sex_check.tsv"),
    )
    .unwrap();
    let slice = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../tests/fixtures/map/gencode-slice/reference");
    for f in ["gencode.v50.basic.genes.gtf", "gencode.v50.basic.exons.bed"] {
        fs::copy(slice.join(f), root.join("reference").join(f)).unwrap();
    }
    root
}

struct Row {
    chrom: usize,
    s: u64,
    e: u64,
    want: Option<String>,
}

impl Span for Row {
    fn chrom(&self) -> usize {
        self.chrom
    }
    fn span(&self) -> (u64, u64) {
        (self.s, self.e)
    }
}

fn rows() -> Vec<Row> {
    fs::read_to_string(data().join("rows.tsv"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Row {
                chrom: chrom_index(f[0]).unwrap(),
                s: f[1].parse().unwrap(),
                e: f[2].parse().unwrap(),
                want: (f[3] != "-").then(|| f[3].to_owned()),
            }
        })
        .collect()
}

fn placed(root: &Path, reference: &Reference) -> Placed {
    let stdin = fs::read(data().join("DEMO-005.json")).unwrap();
    let f = facts::parse(&stdin).unwrap();
    assert!(matches!(f.findings, Findings::Read(_)));
    findings::read(&f.findings, root, "DEMO-005", reference)
}

#[test]
fn every_row_over_a_closable_finding_carries_its_id() {
    let root = genome_dir("veil");
    let reference = Reference::load(&assets()).unwrap();
    let p = placed(&root, &reference);

    // The three act findings with a place, then APOE, in the page's order.
    // The act finding with no place veils nothing and is not listed.
    assert_eq!(
        p.veil.ids,
        [
            "sf-brca2-c-1a-g",
            "cpsr-brca2-c-2c-t",
            "two-hits-cftr",
            "apoe"
        ]
    );

    let rows = rows();
    assert_eq!(rows.len(), 11);
    let tags = p.veil.tags(&rows);
    for (r, t) in rows.iter().zip(&tags) {
        let got = t.map(|k| p.veil.ids[k as usize].as_str());
        assert_eq!(
            got,
            r.want.as_deref(),
            "row {}:{}-{}",
            og_map::reference::CHROMS[r.chrom],
            r.s,
            r.e
        );
    }

    // No bin counts a tagged row: chr13's bin 32 holds five rows, four of
    // them over a closed result; chr19's bin 44 holds only the APOE row.
    let bins = veil::bins(&rows, &tags, &reference, 1_000_000);
    let bin = |c: &str, i: usize| {
        let k = chrom_index(c).unwrap();
        bins.iter().find(|(x, _)| *x == k).map_or(0, |(_, b)| b[i])
    };
    assert_eq!(bin("chr13", 32), 1);
    assert_eq!(bin("chr7", 117), 1);
    assert_eq!(bin("chr19", 44), 0);
    assert_eq!(bin("chr19", 1), 1);
    let counted: u64 = bins.iter().flat_map(|(_, b)| b).sum();
    assert_eq!(counted as usize, veil::untagged(&tags));
    assert_eq!(veil::untagged(&tags), 5);

    fs::remove_dir_all(&root).unwrap();
}

#[test]
fn apoe_is_closable_on_an_adult_page() {
    // DEMO-005 is an adult with no family: nothing about the page closes
    // its act findings, but the APOE result is still tagged, because the
    // page keeps it behind its warning on every page.
    let root = genome_dir("apoe");
    let reference = Reference::load(&assets()).unwrap();
    let p = placed(&root, &reference);
    let apoe = p
        .rows
        .iter()
        .find(|r| r.finding.id == "apoe")
        .expect("APOE is on the map");
    assert!(veil::closable(&apoe.finding));
    assert_eq!(apoe.finding.gate.as_deref(), Some("apoe"));
    // On the APOE gene of the slice (GTF 44903787-44909396).
    assert_eq!(
        (og_map::reference::CHROMS[apoe.chrom], apoe.s, apoe.e),
        ("chr19", 44_903_786, 44_909_396)
    );
    assert_eq!(apoe.zyg, "unknown");
    assert_eq!(p.veil.tag(apoe.chrom, (44_905_000, 44_905_001)), Some(3));
    fs::remove_dir_all(&root).unwrap();
}
