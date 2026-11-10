//! `og map genes` end to end on the public GENCODE 50 slice under
//! `tests/fixtures/map/gencode-slice`: the files rebuild to the input's
//! genes and exons, the same input gives the same bytes, and the cache is
//! reused until an input changes.

use flate2::read::GzDecoder;
use og_map::BuildError;
use og_map::genes::{self, Want};
use og_map::reference::{CHROMS, chrom_index};
use pycompat::json::{self, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

const MADE_BY: &str = "0123456789abcdef0123456789abcdef01234567";

fn slice() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../tests/fixtures/map/gencode-slice")
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("og-map-genes-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// A copy of the slice as a GENOME_DIR the test may change.
fn copy_slice(dir: &Path) -> PathBuf {
    let root = dir.join("genome");
    fs::create_dir_all(root.join("reference")).unwrap();
    for f in ["gencode.v50.basic.genes.gtf", "gencode.v50.basic.exons.bed"] {
        fs::copy(
            slice().join("reference").join(f),
            root.join("reference").join(f),
        )
        .unwrap();
    }
    root
}

fn decode(doc: &[u8]) -> Value {
    let mut text = String::new();
    GzDecoder::new(doc).read_to_string(&mut text).unwrap();
    json::loads(&text).unwrap()
}

fn get<'v>(v: &'v Value, k: &str) -> &'v Value {
    match v {
        Value::Object(m) => m.get(k).unwrap_or_else(|| panic!("no key {k}")),
        _ => panic!("not an object at {k}"),
    }
}

fn ints(v: &Value) -> Vec<u64> {
    match v {
        Value::Array(a) => a
            .iter()
            .map(|x| match x {
                Value::Int(i) => u64::try_from(*i).unwrap(),
                _ => panic!("not an int"),
            })
            .collect(),
        _ => panic!("not an array"),
    }
}

fn strs(v: &Value) -> Vec<String> {
    match v {
        Value::Array(a) => a
            .iter()
            .map(|x| match x {
                Value::Str(s) => s.clone(),
                _ => panic!("not a string"),
            })
            .collect(),
        _ => panic!("not an array"),
    }
}

/// Rebuilds (s, e) from a gap column and a length column, as the page does.
fn spans(cols: &Value) -> Vec<(u64, u64)> {
    let mut s = 0;
    ints(get(cols, "ds"))
        .into_iter()
        .zip(ints(get(cols, "l")))
        .map(|(d, l)| {
            s += d;
            (s, s + l)
        })
        .collect()
}

/// The slice's gene lines, read here on their own: (chrom, s, e, strand,
/// name, type) with s = start - 1, for the map's chromosomes.
fn expected_genes() -> Vec<(String, u64, u64, String, String, String)> {
    let text = fs::read_to_string(slice().join("reference/gencode.v50.basic.genes.gtf")).unwrap();
    let attr = |a: &str, k: &str| {
        let i = a.find(&format!("{k} \"")).unwrap() + k.len() + 2;
        a[i..i + a[i..].find('"').unwrap()].to_owned()
    };
    text.lines()
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| chrom_index(f[0]).is_some())
        .map(|f| {
            (
                f[0].to_owned(),
                f[3].parse::<u64>().unwrap() - 1,
                f[4].parse().unwrap(),
                f[6].to_owned(),
                attr(f[8], "gene_name"),
                attr(f[8], "gene_type"),
            )
        })
        .collect()
}

/// The union of the slice's exons per chromosome, found base by base: a
/// different way from the crate's sort-and-merge.
fn expected_blocks() -> BTreeMap<String, Vec<(u64, u64)>> {
    let text = fs::read_to_string(slice().join("reference/gencode.v50.basic.exons.bed")).unwrap();
    let mut bases: BTreeMap<String, BTreeSet<u64>> = BTreeMap::new();
    for l in text.lines() {
        let f: Vec<&str> = l.split('\t').collect();
        if chrom_index(f[0]).is_none() {
            continue;
        }
        let (s, e): (u64, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        bases.entry(f[0].to_owned()).or_default().extend(s..e);
    }
    bases
        .into_iter()
        .map(|(c, set)| {
            let mut blocks: Vec<(u64, u64)> = Vec::new();
            for b in set {
                match blocks.last_mut() {
                    Some(last) if last.1 == b => last.1 = b + 1,
                    _ => blocks.push((b, b + 1)),
                }
            }
            (c, blocks)
        })
        .collect()
}

fn fetch(root: &Path, cache: &Path, want: Want) -> Result<Vec<u8>, BuildError> {
    genes::genes_cached(root, cache, MADE_BY, want)
}

#[test]
fn genes_and_blocks_round_trip_through_the_chromosome_files() {
    let dir = scratch("round-trip");
    let cache = dir.join("cache");
    let want_genes = expected_genes();
    let want_blocks = expected_blocks();
    assert!(
        want_genes.len() > 250,
        "the slice holds {} genes",
        want_genes.len()
    );
    let mut got_genes = Vec::new();
    for (c, name) in CHROMS.iter().enumerate() {
        let doc = decode(&fetch(&slice(), &cache, Want::Chrom(c)).unwrap());
        assert_eq!(get(&doc, "v"), &Value::Int(1));
        assert_eq!(get(&doc, "gencode"), &Value::Str("v50".into()));
        assert_eq!(get(&doc, "chrom"), &Value::Str((*name).into()));
        let g = get(&doc, "g");
        let rows = spans(g);
        let (strand, names, kinds) = (
            strs(get(g, "strand")),
            strs(get(g, "name")),
            strs(get(g, "type")),
        );
        assert!(
            rows.windows(2).all(|w| w[0] <= w[1]),
            "{name}: sorted by start"
        );
        for (i, (s, e)) in rows.into_iter().enumerate() {
            got_genes.push((
                (*name).to_owned(),
                s,
                e,
                strand[i].clone(),
                names[i].clone(),
                kinds[i].clone(),
            ));
        }
        let blocks = spans(get(&doc, "x"));
        assert_eq!(
            blocks,
            want_blocks.get(*name).cloned().unwrap_or_default(),
            "{name}: exon blocks"
        );
    }
    let mut want_sorted = want_genes.clone();
    want_sorted.sort();
    got_genes.sort();
    assert_eq!(got_genes, want_sorted);
    // The slice has genes off the map (chrM), and the map has a chromosome
    // with none in the slice (chr3): both cases are covered.
    assert!(
        fs::read_to_string(slice().join("reference/gencode.v50.basic.genes.gtf"))
            .unwrap()
            .contains("chrM\t")
    );
    assert!(!want_genes.iter().any(|g| g.0 == "chr3"));
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn names_file_lists_every_gene_by_name() {
    let dir = scratch("names");
    let doc = decode(&fetch(&slice(), &dir.join("cache"), Want::Names).unwrap());
    let names = strs(get(&doc, "name"));
    let chroms = strs(get(&doc, "chrom"));
    let (s, e) = (ints(get(&doc, "s")), ints(get(&doc, "e")));
    assert!(names.windows(2).all(|w| w[0] <= w[1]), "sorted by name");
    let mut got: Vec<_> = (0..names.len())
        .map(|i| (names[i].clone(), chroms[i].clone(), s[i], e[i]))
        .collect();
    let mut want: Vec<_> = expected_genes()
        .into_iter()
        .map(|g| (g.4, g.0, g.1, g.2))
        .collect();
    got.sort();
    want.sort();
    assert_eq!(got, want);
    let cyp = names.iter().position(|n| n == "CYP2D6").unwrap();
    // GENCODE 50 line: chr22 42120962 42131727 (1-based, inclusive).
    assert_eq!(
        (chroms[cyp].as_str(), s[cyp], e[cyp]),
        ("chr22", 42120961, 42131727)
    );
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn same_input_same_bytes() {
    let inputs = genes::locate(&slice()).unwrap();
    let a = genes::write_all(&genes::load(&inputs).unwrap()).unwrap();
    let b = genes::write_all(&genes::load(&inputs).unwrap()).unwrap();
    assert_eq!(a, b);
    assert_eq!(a.len(), 25);
    for (_, bytes) in &a {
        // gzip header: mtime 0 (bytes 4-7) and no file name (FLG bit 3 off).
        assert_eq!(&bytes[4..8], &[0, 0, 0, 0]);
        assert_eq!(bytes[3] & 0x08, 0);
    }
}

#[test]
fn the_cache_is_reused_until_an_input_changes() {
    let dir = scratch("cache");
    let root = copy_slice(&dir);
    let cache = dir.join("cache");
    // A set in progress left by a killed build an hour ago is removed; one a
    // live build is still writing is kept.
    let stale = cache.join("genes/.0000000000000000.tmp-1");
    let live = cache.join("genes/.1111111111111111.tmp-2");
    fs::create_dir_all(&stale).unwrap();
    fs::create_dir_all(&live).unwrap();
    let old = std::time::SystemTime::now() - std::time::Duration::from_secs(7200);
    fs::File::open(&stale).unwrap().set_modified(old).unwrap();
    let first = fetch(&root, &cache, Want::Chrom(6)).unwrap();
    assert!(!stale.exists(), "a stale set in progress is removed");
    assert!(live.is_dir(), "a live set in progress is kept");
    fs::remove_dir(&live).unwrap();
    let sets: Vec<PathBuf> = fs::read_dir(cache.join("genes"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .collect();
    assert_eq!(sets.len(), 1, "one set: {sets:?}");
    assert_eq!(fs::read_dir(&sets[0]).unwrap().count(), 25);

    // A planted byte in the cached file comes back: the set was not rebuilt.
    fs::write(sets[0].join("genes.chr7.json.gz"), b"cached").unwrap();
    assert_eq!(fetch(&root, &cache, Want::Chrom(6)).unwrap(), b"cached");

    // A changed GTF is a new stamp: the files are built again, and the old
    // set is removed.
    let gtf = root.join("reference/gencode.v50.basic.genes.gtf");
    let mut text = fs::read_to_string(&gtf).unwrap();
    text.push_str(
        "chr7\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id \"ENSG0.1\"; gene_type \"lncRNA\"; gene_name \"PLANTED\";\n",
    );
    fs::write(&gtf, text).unwrap();
    let rebuilt = fetch(&root, &cache, Want::Chrom(6)).unwrap();
    assert_ne!(rebuilt, first);
    assert!(strs(get(get(&decode(&rebuilt), "g"), "name")).contains(&"PLANTED".to_owned()));
    let now: Vec<PathBuf> = fs::read_dir(cache.join("genes"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .collect();
    assert_eq!(now.len(), 1);
    assert_ne!(now[0], sets[0]);

    // A cache that cannot be written still serves the file.
    let blocked = dir.join("blocked");
    fs::write(&blocked, b"a file, not a folder").unwrap();
    assert_eq!(fetch(&root, &blocked, Want::Chrom(6)).unwrap(), rebuilt);
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn missing_inputs_say_what_is_needed() {
    let dir = scratch("missing");
    let root = copy_slice(&dir);
    let cache = dir.join("cache");
    fs::remove_file(root.join("reference/gencode.v50.basic.exons.bed")).unwrap();
    let Err(BuildError::NoMap(why)) = fetch(&root, &cache, Want::Names) else {
        panic!("a missing exon file is not a gene set");
    };
    assert_eq!(
        why,
        "Gene names need step 08's GENCODE file: reference/gencode.v50.basic.exons.bed is missing (scripts/setup.sh installs it)"
    );
    fs::remove_file(root.join("reference/gencode.v50.basic.genes.gtf")).unwrap();
    let Err(BuildError::NoMap(why)) = fetch(&root, &cache, Want::Names) else {
        panic!("a missing GTF is not a gene set");
    };
    assert!(
        why.starts_with("Gene names need step 08's GENCODE file: "),
        "{why}"
    );
    assert!(!cache.exists(), "nothing is cached without inputs");
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_broken_line_is_named() {
    let dir = scratch("broken");
    let root = copy_slice(&dir);
    let bed = root.join("reference/gencode.v50.basic.exons.bed");
    let mut text = fs::read_to_string(&bed).unwrap();
    text.push_str("chr1\t500\t400\n");
    let line = text.lines().count();
    fs::write(&bed, text).unwrap();
    let Err(BuildError::NoMap(why)) = fetch(&root, &dir.join("cache"), Want::Names) else {
        panic!("a broken BED is not a gene set");
    };
    assert_eq!(
        why,
        format!(
            "the gene files cannot be built: gencode.v50.basic.exons.bed:{line}: not a BED span: 500-400"
        )
    );
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn the_newest_release_is_read() {
    let dir = scratch("release");
    let root = copy_slice(&dir);
    let r = root.join("reference");
    fs::write(
        r.join("gencode.v9.basic.genes.gtf"),
        "chr1\tH\tgene\t1\t5\t.\t+\t.\tgene_id \"OLD\";\n",
    )
    .unwrap();
    fs::write(r.join("gencode.v9.basic.exons.bed"), "chr1\t0\t5\n").unwrap();
    let inputs = genes::locate(&root).unwrap();
    assert_eq!(inputs.release, "v50");
    assert_eq!(inputs.gtf, r.join("gencode.v50.basic.genes.gtf"));
    fs::remove_dir_all(&dir).unwrap();
}
