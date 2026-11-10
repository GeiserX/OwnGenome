//! `og map build` end to end on the invented runs under `tests/data/genome`:
//! the document rebuilds to the reader's rows, the same input gives the
//! same bytes, and each layer state is reported.

use flate2::read::GzDecoder;
use og_map::facts::Sex;
use og_map::layers::roh;
use og_map::reference::{CHROMS, Reference};
use og_map::{BuildError, Request, State};
use pycompat::json::{self, Value};
use std::io::Read;
use std::path::{Path, PathBuf};

const MALE: &[u8] = br#"{"results":{"person":{"inferred_sex":"male","declared_sex":"male"}}}"#;

fn data() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/genome")
}

fn assets() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets")
}

fn request<'a>(root: &'a Path, assets: &'a Path, run: &'a str, stdin: &'a [u8]) -> Request<'a> {
    Request {
        root,
        run,
        assets,
        made_by: "0123456789abcdef0123456789abcdef01234567",
        stdin,
    }
}

fn decode(doc: &[u8]) -> Value {
    let mut text = String::new();
    GzDecoder::new(doc).read_to_string(&mut text).unwrap();
    json::loads(&text).unwrap()
}

fn get<'v>(v: &'v Value, path: &[&str]) -> &'v Value {
    path.iter().fold(v, |v, k| match v {
        Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no key {k}")),
        _ => panic!("not an object at {k}"),
    })
}

fn column(v: &Value) -> Vec<i128> {
    match v {
        Value::Array(a) => a
            .iter()
            .map(|x| match x {
                Value::Int(i) => *i,
                _ => panic!("not an int"),
            })
            .collect(),
        _ => panic!("not an array"),
    }
}

/// Rebuilds (s, e) from a gap column and a length column, as the page does.
fn spans(cols: &Value) -> Vec<(i128, i128)> {
    let mut s = 0;
    column(get(cols, &["ds"]))
        .into_iter()
        .zip(column(get(cols, &["l"])))
        .map(|(d, l)| {
            s += d;
            (s, s + l)
        })
        .collect()
}

#[test]
fn roh_rows_round_trip_through_the_document() {
    let (root, assets) = (data(), assets());
    let req = request(&root, &assets, "DEMO-001", MALE);
    let model = og_map::build_person(&req).unwrap();
    let doc = decode(&og_map::doc::write(&model).unwrap());

    let file = std::fs::File::open(root.join("DEMO-001/vcf/DEMO-001_roh.txt")).unwrap();
    let reference = Reference::load(&assets).unwrap();
    let (rows, all) = roh::parse(
        std::io::BufReader::new(file),
        "DEMO-001_roh.txt",
        &reference,
        Sex::Male,
    )
    .unwrap();
    assert_eq!((all, rows.len()), (9, 7));

    let per = get(&doc, &["per", "DEMO-001", "roh"]);
    let Value::Object(by_chrom) = per else {
        panic!("roh is not an object")
    };
    let mut rebuilt = Vec::new();
    for (chrom, cols) in by_chrom {
        let k = CHROMS.iter().position(|c| c == chrom).unwrap();
        let long = column(get(cols, &["long"]));
        let art = column(get(cols, &["art"]));
        let markers = column(get(cols, &["markers"]));
        for (i, (s, e)) in spans(cols).into_iter().enumerate() {
            rebuilt.push((k, s, e, markers[i], long[i] == 1, art[i]));
        }
    }
    let want: Vec<_> = rows
        .iter()
        .map(|r| {
            (
                r.chrom,
                i128::from(r.s),
                i128::from(r.e),
                i128::from(r.markers),
                r.long,
                r.art as i128,
            )
        })
        .collect();
    assert_eq!(rebuilt, want);
    // The first run, RG 1000001-7500000, is [1000000, 7500000) on the map.
    assert_eq!((rebuilt[0].1, rebuilt[0].2), (1_000_000, 7_500_000));

    // The bands rebuild to the asset's rows.
    for c in &reference.chroms {
        let got = spans(get(&doc, &["ref", "bands", c.name]));
        let want: Vec<(i128, i128)> = c
            .bands
            .iter()
            .map(|b| (i128::from(b.s), i128::from(b.e)))
            .collect();
        assert_eq!(got, want, "{}", c.name);
    }

    let layer = get(&doc, &["layers", "roh"]);
    assert_eq!(get(layer, &["state"]), &Value::Str("ok".into()));
    assert_eq!(
        get(layer, &["source"]),
        &Value::Str("DEMO-001/vcf/DEMO-001_roh.txt".into())
    );
    assert_eq!(get(layer, &["count"]), &Value::Int(7));
    assert_eq!(get(layer, &["totals", "all"]), &Value::Int(9));
    // Two long runs (chr1 and chr7); chr13's lies in its short arm and the
    // male's chrX run is flagged.
    assert_eq!(
        get(&doc, &["text", "genome"]),
        &json::loads(r#"[{"t":"2 runs","b":true},{"t":" of homozygosity of 5 Mb or more, 12.5 Mb in all.","b":false}]"#).unwrap()
    );
    assert_eq!(
        get(&doc, &["people"]),
        &json::loads(r#"[{"id":"DEMO-001","label":"DEMO-001","sex":"male","has_y":true}]"#)
            .unwrap()
    );
    let bins = column(get(&doc, &["per", "DEMO-001", "roh_bins", "chr1"]));
    assert_eq!(bins.len(), 249);
    assert_eq!(
        &bins[..10],
        [
            0, 1_000_000, 1_000_000, 1_000_000, 1_000_000, 1_000_000, 1_000_000, 500_000, 0,
            700_000
        ]
    );
}

#[test]
fn same_input_same_bytes() {
    let (root, assets) = (data(), assets());
    let req = request(&root, &assets, "DEMO-001", MALE);
    let a = og_map::doc::write(&og_map::build_person(&req).unwrap()).unwrap();
    let b = og_map::doc::write(&og_map::build_person(&req).unwrap()).unwrap();
    assert_eq!(a, b);
    // One gzip member, mtime 0, no file name (FLG has no FNAME bit).
    assert_eq!(&a[..4], [0x1f, 0x8b, 8, 0]);
    assert_eq!(&a[4..8], [0, 0, 0, 0]);
    let text = {
        let mut s = String::new();
        GzDecoder::new(&a[..]).read_to_string(&mut s).unwrap();
        s
    };
    assert!(
        text.starts_with(r#"{"v":1,"build":"GRCh38","kind":"person","run":"DEMO-001","stamp":""#)
    );
    let keys: Vec<String> = match decode(&a) {
        Value::Object(m) => m.keys().cloned().collect(),
        _ => panic!(),
    };
    assert_eq!(
        keys,
        [
            "v", "build", "kind", "run", "stamp", "made_by", "people", "ref", "layers", "per",
            "dict", "text"
        ]
    );
}

#[test]
fn layer_states() {
    let (root, assets) = (data(), assets());
    let state = |run: &str| {
        let m = og_map::build_person(&request(&root, &assets, run, b"")).unwrap();
        (
            m.roh.layer.state,
            m.roh.layer.source,
            m.roh.layer.note,
            m.roh.layer.count,
        )
    };
    // The pipeline's path, read when the single-step one is absent; no RG row.
    assert_eq!(
        state("DEMO-002"),
        (State::Ok, "DEMO-002/roh/DEMO-002_roh.txt".into(), None, 0)
    );
    assert_eq!(
        state("DEMO-003"),
        (
            State::Unreadable,
            "DEMO-003/vcf/DEMO-003_roh.txt".into(),
            Some("DEMO-003_roh.txt:5: expected 8 columns, found 5".into()),
            0
        )
    );
    let tmp = std::env::temp_dir().join(format!("og-map-states-{}", std::process::id()));
    std::fs::create_dir_all(tmp.join("DEMO-004")).unwrap();
    let m = og_map::build_person(&request(&tmp, &assets, "DEMO-004", b"")).unwrap();
    assert_eq!(m.roh.layer.state, State::Missing);
    assert_eq!(m.roh.layer.source, "DEMO-004/vcf/DEMO-004_roh.txt");
    assert_eq!(
        m.roh.layer.note.as_deref(),
        Some("Step 11 has not written DEMO-004_roh.txt for this run.")
    );
    let doc = decode(&og_map::doc::write(&m).unwrap());
    assert_eq!(get(&doc, &["text", "genome"]), &Value::Array(vec![]));
    assert_eq!(
        get(&doc, &["people"]),
        &json::loads(r#"[{"id":"DEMO-004","label":"DEMO-004","sex":"unknown","has_y":null}]"#)
            .unwrap()
    );

    // No folder for the run, or a name that is not one: no map.
    let e = og_map::build_person(&request(&tmp, &assets, "DEMO-005", b"")).unwrap_err();
    assert!(matches!(e, BuildError::NoMap(_)), "{e}");
    for bad in ["../DEMO-001", ".x", "", "a/b"] {
        let e = og_map::build_person(&request(&root, &assets, bad, b"")).unwrap_err();
        assert!(matches!(e, BuildError::Usage(_)), "{bad}: {e}");
    }
    let e = og_map::build_person(&request(&root, &assets, "DEMO-001", b"[")).unwrap_err();
    assert!(matches!(e, BuildError::Usage(_)), "{e}");
    std::fs::remove_dir_all(&tmp).unwrap();
}

#[test]
fn cache_hit_and_a_changed_input() {
    let (root, assets) = (data(), assets());
    let cache = std::env::temp_dir().join(format!("og-map-cache-{}", std::process::id()));
    let req = request(&root, &assets, "DEMO-001", MALE);
    let first = og_map::build_cached(&req, &cache).unwrap();
    let path = og_map::stamp::cache_path(&cache, "DEMO-001", &req.stamp());
    assert_eq!(std::fs::read(&path).unwrap(), first);
    // A hit returns the stored bytes, whatever they are.
    std::fs::write(&path, b"cached").unwrap();
    assert_eq!(og_map::build_cached(&req, &cache).unwrap(), b"cached");
    // Other facts on stdin are another stamp: built again, the old one removed.
    let req2 = request(&root, &assets, "DEMO-001", b"{}");
    let second = og_map::build_cached(&req2, &cache).unwrap();
    assert_ne!(second, first);
    assert!(!path.exists());
    assert_eq!(
        std::fs::read_dir(cache.join("DEMO-001")).unwrap().count(),
        1
    );
    std::fs::remove_dir_all(&cache).unwrap();
}
