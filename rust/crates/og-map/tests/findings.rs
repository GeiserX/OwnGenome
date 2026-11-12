//! The findings layer in the person document, on the synthetic person under
//! `tests/data/facts` and the public GENCODE slice: each finding with a
//! place lands on it, a gene-only one on its gene, one with neither is
//! counted and not drawn, the hidden APOE text never enters the document,
//! and the bytes and the cache behave as before.

use flate2::read::GzDecoder;
use og_map::{Request, State};
use pycompat::json::{self, Value};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

const MADE_BY: &str = "0123456789abcdef0123456789abcdef01234567";

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/facts")
}

fn assets() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets")
}

/// A GENOME_DIR with the synthetic run and, when `genes`, the GENCODE slice.
fn genome_dir(name: &str, genes: bool) -> PathBuf {
    let root = std::env::temp_dir().join(format!("og-map-f-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(root.join("DEMO-005")).unwrap();
    fs::copy(
        data().join("DEMO-005/DEMO-005_sex_check.tsv"),
        root.join("DEMO-005/DEMO-005_sex_check.tsv"),
    )
    .unwrap();
    if genes {
        fs::create_dir_all(root.join("reference")).unwrap();
        let slice = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../../tests/fixtures/map/gencode-slice/reference");
        for f in ["gencode.v50.basic.genes.gtf", "gencode.v50.basic.exons.bed"] {
            fs::copy(slice.join(f), root.join("reference").join(f)).unwrap();
        }
    }
    root
}

fn request<'a>(root: &'a Path, assets: &'a Path, stdin: &'a [u8]) -> Request<'a> {
    Request {
        root,
        run: "DEMO-005",
        assets,
        made_by: MADE_BY,
        stdin,
    }
}

fn text(doc: &[u8]) -> String {
    let mut s = String::new();
    GzDecoder::new(doc).read_to_string(&mut s).unwrap();
    s
}

fn get<'v>(v: &'v Value, path: &[&str]) -> &'v Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no key {k}")),
        _ => panic!("not an object at {k}"),
    })
}

fn js(s: &str) -> Value {
    json::loads(s).unwrap()
}

#[test]
fn findings_land_on_their_place_and_the_rest_are_counted() {
    let (root, assets) = (genome_dir("doc", true), assets());
    let stdin = fs::read(data().join("DEMO-005.json")).unwrap();
    let model = og_map::build_person(&request(&root, &assets, &stdin)).unwrap();
    let raw = text(&og_map::doc::write(&model).unwrap());
    let doc = js(&raw);

    let f = get(&doc, &["per", "DEMO-005", "findings"]);
    // Map order: chromosome, then start.
    assert_eq!(
        get(f, &["id"]),
        &js(
            r#"["two-hits-cftr","sf-brca2-c-1a-g","clinvar-brca2-chr13-32340000-a-g","cpsr-brca2-c-2c-t","clinvar-tp53-chr17-7675000-c-t","apoe","pgx-cyp2d6"]"#
        )
    );
    assert_eq!(
        get(f, &["chrom"]),
        &js(r#"["chr7","chr13","chr13","chr13","chr17","chr19","chr22"]"#)
    );
    // A position is [pos - 1, pos); a gene-only finding takes its gene's
    // span from the GTF (1-based inclusive there).
    assert_eq!(
        get(f, &["pos"]),
        &js("[null,32340000,32340000,32390000,7675000,null,null]")
    );
    assert_eq!(
        get(f, &["s"]),
        &js("[117287119,32339999,32339999,32389999,7674999,44903786,42120961]")
    );
    assert_eq!(
        get(f, &["e"]),
        &js("[117715971,32340000,32340000,32390000,7675000,44909396,42131727]")
    );
    assert_eq!(
        get(f, &["reads"]),
        &js(
            r#"["chr7:117287120-117715971","chr13:32340000","chr13:32340000","chr13:32390000","chr17:7675000","chr19:44903787-44909396","chr22:42120962-42131727"]"#
        )
    );
    assert_eq!(
        get(f, &["gate"]),
        &js(r#"[null,null,null,null,null,"apoe",null]"#)
    );
    assert_eq!(
        get(f, &["zyg"]),
        &js(r#"["unknown","het","unknown","unknown","het","unknown","unknown"]"#)
    );
    assert_eq!(
        get(f, &["hgvs"]),
        &js(r#"[null,"c.1A>G",null,null,null,null,null]"#)
    );
    // APOE is gated, so it comes first in dict.veil. The open ClinVar row
    // at the act finding's base carries that finding's tag.
    assert_eq!(get(f, &["veil"]), &js("[3,1,1,2,-1,0,-1]"));
    // `veils` holds every closable finding a row overlaps (one each here;
    // src/veil.rs tests a row over several), its first being `veil`.
    assert_eq!(get(f, &["veils"]), &js("[[3],[1],[1],[2],[],[0],[]]"));
    assert_eq!(
        get(f, &["href"]),
        &js(
            r#"["finding/two-hits-cftr/","finding/sf-brca2-c-1a-g/","finding/clinvar-brca2-chr13-32340000-a-g/","finding/cpsr-brca2-c-2c-t/","finding/clinvar-tp53-chr17-7675000-c-t/","finding/apoe/","finding/pgx-cyp2d6/"]"#
        )
    );
    assert_eq!(
        get(f, &["steps"]),
        &js(r#"[["29"],["24"],["13"],["30"],["13"],["42"],["17"]]"#)
    );
    assert_eq!(
        get(f, &["copy"]),
        &js("[null,null,null,null,null,null,null]")
    );
    let Value::Array(phrases) = get(f, &["phrase"]) else {
        panic!()
    };
    assert_eq!(
        phrases[1],
        Value::Str("Take to a doctor: Secondary finding in BRCA2, on chromosome 13.".into())
    );
    assert_eq!(
        phrases[5],
        Value::Str("APOE, on chromosome 19: kept behind its warning on every page until opened on purpose.".into())
    );
    assert_eq!(
        get(&doc, &["dict", "veil"]),
        &js(r#"["apoe","sf-brca2-c-1a-g","cpsr-brca2-c-2c-t","two-hits-cftr"]"#)
    );

    // Counts cover the rows no closed result hides: TP53 and CYP2D6; the
    // sex mismatch, LPA (not in the slice) and the chrM hit add to `all`
    // and stay on the page. The act finding with no place is in no count,
    // and neither is the open ClinVar row tagged by the act finding it
    // sits on.
    let layer = get(&doc, &["layers", "findings"]);
    assert_eq!(get(layer, &["state"]), &Value::Str("ok".into()));
    assert_eq!(
        get(layer, &["source"]),
        &Value::Str("DEMO-005/summary.json".into())
    );
    assert_eq!(get(layer, &["count"]), &Value::Int(2));
    assert_eq!(get(layer, &["totals"]), &js(r#"{"all":5,"kept":2}"#));
    assert_eq!(
        get(layer, &["note"]),
        &Value::Str("3 findings with no place on the map stay on the page.".into())
    );
    assert_eq!(
        get(&doc, &["text", "findings", "many"]),
        &js(r#"[{"t":"{n} findings","b":true},{"t":" have a place on the map.","b":false}]"#)
    );

    // The sex step 16 measured wins over the declared one; the GENCODE
    // release is the one installed.
    assert_eq!(
        get(&doc, &["people"]),
        &js(r#"[{"id":"DEMO-005","label":"DEMO-005","sex":"female","has_y":false}]"#)
    );
    assert_eq!(get(&doc, &["ref", "gencode"]), &Value::Str("v50".into()));

    // The APOE result itself never enters the document, only that it exists.
    assert!(!raw.contains("HIDDEN-TEXT-MARKER"));
    assert!(!raw.contains("Alzheimer"));

    let keys: Vec<&String> = match get(&doc, &["layers"]) {
        Value::Object(m) => m.keys().collect(),
        _ => panic!(),
    };
    assert_eq!(keys, ["findings", "roh"]);
    fs::remove_dir_all(&root).unwrap();
}

#[test]
fn without_genes_or_results() {
    let (root, assets) = (genome_dir("nogenes", false), assets());
    let stdin = fs::read(data().join("DEMO-005.json")).unwrap();
    let m = og_map::build_person(&request(&root, &assets, &stdin)).unwrap();
    // No GENCODE installed: only the findings with a position are placed.
    let ids: Vec<&str> = m
        .findings
        .rows
        .iter()
        .map(|r| r.finding.id.as_str())
        .collect();
    assert_eq!(
        ids,
        [
            "sf-brca2-c-1a-g",
            "clinvar-brca2-chr13-32340000-a-g",
            "cpsr-brca2-c-2c-t",
            "clinvar-tp53-chr17-7675000-c-t"
        ]
    );
    assert_eq!(m.findings.layer.totals, Some((5, 1)));
    assert_eq!(m.findings.layer.count, 1);
    assert_eq!(m.gencode, None);
    assert_eq!(
        m.findings.veil.ids,
        ["sf-brca2-c-1a-g", "cpsr-brca2-c-2c-t"]
    );

    // No results body: the layer is missing, with its words.
    let m = og_map::build_person(&request(&root, &assets, b"")).unwrap();
    assert_eq!(m.findings.layer.state, State::Missing);
    assert_eq!(
        m.findings.layer.note.as_deref(),
        Some("No results for this run yet: the findings come from its summary (step 24).")
    );
    let m =
        og_map::build_person(&request(&root, &assets, br#"{"results":{"findings":7}}"#)).unwrap();
    assert_eq!(m.findings.layer.state, State::Unreadable);
    fs::remove_dir_all(&root).unwrap();
}

#[test]
fn same_bytes_and_the_cache() {
    let (root, assets) = (genome_dir("cache", true), assets());
    let stdin = fs::read(data().join("DEMO-005.json")).unwrap();
    let req = request(&root, &assets, &stdin);
    let a = og_map::doc::write(&og_map::build_person(&req).unwrap()).unwrap();
    let b = og_map::doc::write(&og_map::build_person(&req).unwrap()).unwrap();
    assert_eq!(a, b);

    let cache = root.join("cache");
    let first = og_map::build_cached(&req, &cache).unwrap();
    assert_eq!(first, a);
    let path = og_map::stamp::cache_path(&cache, "DEMO-005", &req.stamp());
    fs::write(&path, b"cached").unwrap();
    assert_eq!(og_map::build_cached(&req, &cache).unwrap(), b"cached");
    // A changed sex file is a new stamp.
    let tsv = root.join("DEMO-005/DEMO-005_sex_check.tsv");
    fs::write(&tsv, "inferred_sex\tCNchrX\tCNchrY\nmale\t1.0\t1.0\n").unwrap();
    let again = og_map::build_cached(&req, &cache).unwrap();
    assert!(text(&again).contains(r#""sex":"male","has_y":true"#));
    // Removing the genes file is a new stamp too.
    fs::remove_file(root.join("reference/gencode.v50.basic.genes.gtf")).unwrap();
    let no_genes = og_map::build_cached(&req, &cache).unwrap();
    assert!(text(&no_genes).contains(r#""gencode":null"#));
    fs::remove_dir_all(&root).unwrap();
}
