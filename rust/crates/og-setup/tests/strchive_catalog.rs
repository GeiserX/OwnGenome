//! `og str catalog` on the two STRchive files `versions.env` pins, and
//! Stranger on the catalog it builds: what the old
//! `tests/test_str_catalog.py` checked, by meaning, in cargo test.
//!
//! 1. The files: downloaded once into a cache folder (`STRCHIVE_DIR`, else
//!    `<tmp>/pgp-strchive-<version>`) and checked against their sha256; a
//!    copy with one byte flipped fails that check.
//! 2. The converter: one entry per STRchive ExpansionHunter locus (82), no
//!    name twice, every ExpansionHunter LocusId traced to one entry by its
//!    STRchive id; the report's key loci keep their names; the thresholds
//!    of HTT, FMR1 and SAMD12; the motif loci; the loci without thresholds
//!    and why; and a planted ExpansionHunter locus with no STRchive record
//!    stops it (exit 2, nothing written).
//! 3. Stranger (`STRANGER_IMAGE`, through docker) on the catalog: HTT 26,
//!    27 and 36 give normal, pre_mutation and full_mutation; SAMD12 1 and 2
//!    normal and pre_mutation; MUC1 no status; control: Stranger's bundled
//!    GRCh38 catalog does not call HTT 27 pre_mutation.
//! 4. The repeats STRchive leaves out (`assets/str_extra_loci.json`,
//!    `--extra`): merged after STRchive's 82 entries, copied as they are,
//!    with ranges rows and no Stranger thresholds; each region is, in
//!    GRCh38 (NCBI's GRCh38.p14 sequence), a whole run of its motif that
//!    the bases beside it do not continue; controls: an extra catalog
//!    without UGT1A1 is caught, four planted bad entries each stop the
//!    converter, and a region shifted by one base fails the sequence check.
//!
//! The tests download (GitHub and NCBI) and run docker, so they run where
//! `CI` is set (every CI job) or with `OG_NET_TESTS=1`; elsewhere each says
//! it skipped. In CI a failed download or a missing docker fails the test.
//! What `og str motif` counts on the catalog's motifs is `tests/parity.rs`
//! of og-bam and the e2e case `str-catalog-1-strchive`.

use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use indexmap::IndexMap;
use pycompat::json::{self, Dump, Value};

fn repo() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

/// Whether the tests that download and run docker run here.
fn online(test: &str) -> bool {
    let on = std::env::var_os("CI").is_some_and(|v| !v.is_empty())
        || std::env::var("OG_NET_TESTS").as_deref() == Ok("1");
    if !on {
        eprintln!("skip: {test} downloads STRchive and runs docker; set OG_NET_TESTS=1 (CI always runs it)");
    }
    on
}

/// A value of versions.env.
fn env(name: &str) -> String {
    let text = std::fs::read_to_string(repo().join("versions.env")).unwrap();
    let prefix = format!("{name}=");
    let line = text
        .lines()
        .find(|l| l.starts_with(&prefix))
        .unwrap_or_else(|| panic!("no {name} line in versions.env"));
    let v = line[prefix.len()..].split('#').next().unwrap().trim();
    v.trim_matches('"').to_owned()
}

fn sha256(path: &Path) -> String {
    for (cmd, args) in [("sha256sum", &[][..]), ("shasum", &["-a", "256"][..])] {
        if let Ok(o) = Command::new(cmd).args(args).arg(path).output()
            && o.status.success()
        {
            let out = String::from_utf8(o.stdout).unwrap();
            return out.split_whitespace().next().unwrap_or_default().to_owned();
        }
    }
    panic!("neither sha256sum nor shasum could hash {}", path.display());
}

/// GET `url` into `dest` with curl, three tries.
fn download(url: &str, dest: &Path) -> Result<(), String> {
    let o = Command::new("curl")
        .args(["-sSLf", "--retry", "3", "--retry-delay", "5", "--max-time", "120", "-o"])
        .arg(dest)
        .arg(url)
        .output()
        .map_err(|e| format!("curl: {e}"))?;
    if o.status.success() {
        Ok(())
    } else {
        Err(format!("curl {url}: {}", String::from_utf8_lossy(&o.stderr)))
    }
}

struct Files {
    release: String,
    eh: String,
    loci: String,
}

/// The release and the two pinned files, downloaded once and checked.
fn files() -> &'static Files {
    static FILES: OnceLock<Files> = OnceLock::new();
    FILES.get_or_init(|| {
        let rel = env("STRCHIVE_VERSION");
        let dir = std::env::var_os("STRCHIVE_DIR")
            .map_or_else(|| std::env::temp_dir().join(format!("pgp-strchive-{rel}")), PathBuf::from);
        std::fs::create_dir_all(&dir).unwrap();
        let mut out = Vec::new();
        for (name, var) in [
            (format!("STRchive-disease-loci-{rel}.hg38.expansionhunter.json"), "STRCHIVE_EH_SHA256"),
            (format!("STRchive-loci-{rel}.json"), "STRCHIVE_LOCI_SHA256"),
        ] {
            let path = dir.join(&name);
            if !(path.exists() && sha256(&path) == env(var)) {
                let part = dir.join(format!("{name}.part"));
                let url = format!("https://github.com/dashnowlab/STRchive/releases/download/{rel}/{name}");
                download(&url, &part).unwrap();
                std::fs::rename(&part, &path).unwrap();
            }
            assert_eq!(sha256(&path), env(var), "{name} does not have the sha256 {var} pins");
            out.push(path.to_string_lossy().into_owned());
        }
        Files {
            release: rel,
            eh: out[0].clone(),
            loci: out[1].clone(),
        }
    })
}

/// A scratch folder of its own for one test, removed when dropped.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Scratch {
        let d = std::env::temp_dir().join(format!("og-strchive-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Scratch(d)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn load(path: &str) -> Value {
    json::loads(&std::fs::read_to_string(path).unwrap()).unwrap()
}

fn dump(path: &str, v: &Value) {
    std::fs::write(path, json::dumps(v, &Dump::default()).unwrap()).unwrap();
}

fn list(v: &Value) -> &[Value] {
    match v {
        Value::Array(a) => a,
        other => panic!("not a list: {other:?}"),
    }
}

fn obj(v: &Value) -> &IndexMap<String, Value> {
    match v {
        Value::Object(o) => o,
        other => panic!("not an object: {other:?}"),
    }
}

fn s(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        other => panic!("not text: {other:?}"),
    }
}

fn strs(v: Option<&Value>) -> Option<Vec<&str>> {
    v.map(|v| list(v).iter().map(s).collect())
}

fn int(v: Option<&Value>) -> Option<i128> {
    match v {
        Some(Value::Int(i)) => Some(*i),
        _ => None,
    }
}

/// `og str catalog ARGS` on the pinned files: the exit code.
fn catalog(extra_args: &[&str]) -> u8 {
    let f = files();
    let mut args: Vec<String> = ["--eh", &f.eh, "--loci", &f.loci, "--release", &f.release]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
    args.extend(extra_args.iter().map(|a| (*a).to_owned()));
    og_setup::strchive::cli(&args)
}

#[test]
fn a_flipped_byte_fails_the_checksum() {
    if !online("a_flipped_byte_fails_the_checksum") {
        return;
    }
    let w = Scratch::new("flip");
    let mut data = std::fs::read(&files().eh).unwrap();
    let mid = data.len() / 2;
    data[mid] ^= 1;
    let bad = w.0.join("flipped.json");
    std::fs::write(&bad, data).unwrap();
    assert_ne!(sha256(&bad), env("STRCHIVE_EH_SHA256"), "control: a copy with one byte flipped fails the sha256 check");
}

#[test]
fn the_converter_on_the_pinned_files() {
    if !online("the_converter_on_the_pinned_files") {
        return;
    }
    let w = Scratch::new("convert");
    let (out, rng) = (w.path("catalog.json"), w.path("str_ranges.tsv"));
    assert_eq!(catalog(&["--out", &out, "--ranges", &rng]), 0, "the converter runs on the pinned files");
    let src_v = load(&files().eh);
    let src = list(&src_v);
    let cat_v = load(&out);
    let cat: Vec<&IndexMap<String, Value>> = list(&cat_v).iter().map(obj).collect();
    assert_eq!(src.len(), 82, "STRchive's ExpansionHunter catalog has 82 loci");
    assert_eq!(cat.len(), src.len(), "the catalog has one entry per ExpansionHunter locus");

    let mut by_src: HashMap<&str, Vec<&str>> = HashMap::new();
    for e in &cat {
        by_src.entry(s(&e["SourceId"])).or_default().push(s(&e["LocusId"]));
    }
    let src_ids: BTreeSet<&str> = src.iter().map(|x| s(&obj(x)["LocusId"])).collect();
    assert!(
        by_src.keys().copied().collect::<BTreeSet<_>>() == src_ids && by_src.values().all(|v| v.len() == 1),
        "every ExpansionHunter LocusId has exactly one entry: {:?}",
        by_src.iter().filter(|(_, v)| v.len() != 1).collect::<Vec<_>>()
    );
    let names: Vec<&str> = cat
        .iter()
        .flat_map(|e| match e.get("VariantId") {
            Some(Value::Array(a)) => a.iter().map(s).collect(),
            _ => vec![s(&e["LocusId"])],
        })
        .collect();
    let unique: BTreeSet<&str> = names.iter().copied().collect();
    assert_eq!(unique.len(), names.len(), "no variant name twice");

    let by: HashMap<&str, &IndexMap<String, Value>> = cat.iter().map(|e| (s(&e["LocusId"]), *e)).collect();
    for k in ["HTT", "FMR1", "C9ORF72", "ATXN1", "DMPK", "RFC1"] {
        assert!(by.contains_key(k), "the report's key loci and RFC1 keep their names: {k}");
    }
    for k in ["EIEE1_ARX", "PRTS_ARX", "HFG_HOXA13-I"] {
        assert!(by.contains_key(k), "one gene with two loci keeps STRchive's ids: {k}");
    }
    let htt = by["HTT"];
    assert_eq!(
        (int(htt.get("NormalMax")), int(htt.get("PathologicMin")), strs(htt.get("VariantId"))),
        (Some(26), Some(36), Some(vec!["HTT", "HTT_CCG"])),
        "HTT: normal to 26, pathogenic from 36 (STRchive), the CAG repeat is REPID HTT"
    );
    let fmr1 = by["FMR1"];
    assert_eq!(
        (int(fmr1.get("NormalMax")), int(fmr1.get("PathologicMin"))),
        (Some(44), Some(201)),
        "FMR1: normal to 44, pathogenic from 201"
    );
    let s12 = by["SAMD12"];
    assert_eq!(
        (
            int(s12.get("NormalMax")),
            int(s12.get("PathologicMin")),
            strs(s12.get("BenignMotifs")),
            strs(s12.get("PathogenicMotifs")),
            strs(s12.get("VariantId")),
        ),
        (Some(1), Some(105), Some(vec!["TAAAA"]), Some(vec!["TGAAA"]), Some(vec!["SAMD12_TAAAA", "SAMD12"])),
        "SAMD12 (FAME1): no benign maximum in STRchive, so normal is the reference count (1 TGAAA), pathogenic \
         from 105; TTTTA benign, TTTCA pathogenic"
    );
    let rfc1 = by["RFC1"];
    assert_eq!(
        (strs(rfc1.get("BenignMotifs")), strs(rfc1.get("PathogenicMotifs"))),
        (Some(vec!["AAAAG", "AAAGG"]), Some(vec!["AAGGG", "ACAGG", "CAGGG"])),
        "RFC1: AAAGG is counted benign, CAGGG pathogenic"
    );
    let motif_loci: BTreeSet<&str> = cat
        .iter()
        .filter(|e| e.contains_key("PathogenicMotifs"))
        .map(|e| s(&e["LocusId"]))
        .collect();
    for k in ["RFC1", "SAMD12", "STARD7", "YEATS2", "MARCHF6", "RAPGEF2", "TNRC6A", "RAI1", "BEAN1", "DAB1"] {
        assert!(motif_loci.contains(k), "motif loci: RFC1, the FAME loci, SCA31 and SCA37: {k}");
    }
    let none: BTreeSet<&str> = cat
        .iter()
        .filter(|e| e.contains_key("NoThresholds"))
        .map(|e| s(&e["LocusId"]))
        .collect();
    assert_eq!(
        none.iter().copied().collect::<Vec<_>>(),
        ["CEL", "MIR7-2", "MUC1", "PLIN4", "POLG", "TYMS", "VWA1"],
        "no thresholds: the four repeats longer than two reads, two without a pathogenic size, one contraction"
    );
    assert!(
        none.iter().all(|k| !by[k].contains_key("NormalMax") && !by[k].contains_key("PathologicMin")),
        "a locus without thresholds has no NormalMax, so Stranger gives it no status"
    );
    assert!(
        cat.iter()
            .filter(|e| e.contains_key("NormalMax"))
            .all(|e| int(e.get("NormalMax")) < int(e.get("PathologicMin"))),
        "every threshold pair is normal below pathogenic"
    );
    assert_eq!(
        htt.get("SourceDisplay").map(s),
        Some(format!("STRchive-{}", files().release).as_str()),
        "the release is in SourceDisplay"
    );

    // Control: a planted ExpansionHunter locus with no STRchive record.
    let mut planted = src.to_vec();
    planted.push(json::loads(r#"{"LocusId": "PLANTED_GENE", "LocusStructure": "(CAG)*", "ReferenceRegion": "chr1:1000-1030", "VariantType": "Repeat"}"#).unwrap());
    let eh = w.path("planted_eh.json");
    dump(&eh, &Value::Array(planted));
    let (o2, r2) = (w.path("out.json"), w.path("ranges.tsv"));
    let f = files();
    let rc = og_setup::strchive::cli(
        &["--eh", &eh, "--loci", &f.loci, "--release", &f.release, "--out", &o2, "--ranges", &r2].map(str::to_owned),
    );
    assert!(
        rc == 2 && !Path::new(&o2).exists() && !Path::new(&r2).exists(),
        "control: a planted ExpansionHunter locus with no STRchive record stops the converter (exit 2, no catalog \
         and no ranges table): exit {rc}"
    );
}

// --- Stranger on the catalog -------------------------------------------------------------

const FMT: &str = "GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC";

fn id(flag: &str) -> String {
    let o = Command::new("id").arg(flag).output().unwrap();
    String::from_utf8(o.stdout).unwrap().trim().to_owned()
}

/// STR_STATUS per record id from STRANGER_IMAGE on a VCF of (id, chrom,
/// pos, end, count) records; `catalog` None is the GRCh38 catalog Stranger
/// bundles, as step 9b passes it.
fn stranger(dir: &Path, catalog: Option<&str>, records: &[(&str, &str, u64, u64, u32)]) -> HashMap<String, String> {
    let smoke = std::fs::read_to_string(repo().join("tests/smoke/stranger_input.vcf.in")).unwrap();
    let mut vcf: String = smoke.lines().filter(|l| l.starts_with("##")).map(|l| l.to_owned() + "\n").collect();
    let alts: BTreeSet<u32> = records.iter().map(|r| r.4).collect();
    for a in alts {
        vcf += &format!("##ALT=<ID=STR{a},Description=\"Allele comprised of {a} repeat units\">\n");
    }
    let chroms: BTreeSet<&str> = records.iter().map(|r| r.1).collect();
    for c in chroms {
        vcf += &format!("##contig=<ID={c},length=250000000>\n");
    }
    vcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
    for (i, (rid, chrom, pos, end, n)) in records.iter().enumerate() {
        vcf += &format!(
            "{chrom}\t{pos}\t{rid}.{i}\tA\t<STR{n}>\t.\tPASS\tSVTYPE=STR;END={end};REF=1;RL=3;RU=CAG;VARID={rid};REPID={rid}\t\
{FMT}\t1/1:SPANNING/SPANNING:{n}/{n}:{n}-{n}/{n}-{n}:20/20:0/0:0/0:30.0\n"
        );
    }
    std::fs::write(dir.join("in.vcf"), vcf).unwrap();
    let mut cmd = Command::new("docker");
    cmd.args(["run", "--rm", "--network", "none", "-u", &format!("{}:{}", id("-u"), id("-g")), "-v"])
        .arg(format!("{}:/w", dir.display()))
        .arg(env("STRANGER_IMAGE"));
    match catalog {
        Some(c) => cmd.args(["stranger", "--repeats-file", &format!("/w/{c}"), "/w/in.vcf"]),
        None => cmd.args([
            "sh",
            "-c",
            "exec stranger --repeats-file \"$(python -c 'import os, stranger; print(os.path.join(\
             os.path.dirname(stranger.__file__), \"resources\", \"variant_catalog_grch38.json\"))')\" /w/in.vcf",
        ]),
    };
    let o = cmd.output().unwrap_or_else(|e| panic!("docker: {e}"));
    assert!(o.status.success(), "Stranger runs: {}", String::from_utf8_lossy(&o.stderr));
    String::from_utf8(o.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let status = f[7]
                .split(';')
                .find_map(|kv| kv.strip_prefix("STR_STATUS="))
                .unwrap_or_default();
            (f[2].to_owned(), status.to_owned())
        })
        .collect()
}

#[test]
fn stranger_grades_with_the_catalog() {
    if !online("stranger_grades_with_the_catalog") {
        return;
    }
    let w = Scratch::new("stranger");
    assert_eq!(catalog(&["--out", &w.path("catalog.json")]), 0, "the converter runs on the pinned files");
    // The container runs as the caller and reads the folder.
    assert!(Command::new("chmod").args(["755"]).arg(&w.0).status().unwrap().success());
    let recs = [
        ("HTT", "chr4", 3074877, 3074933, 26),
        ("HTT", "chr4", 3074877, 3074933, 27),
        ("HTT", "chr4", 3074877, 3074933, 36),
        ("SAMD12", "chr8", 118366913, 118366918, 1),
        ("SAMD12", "chr8", 118366913, 118366918, 2),
        ("MUC1", "chr1", 155188506, 155192239, 60),
    ];
    let got = stranger(&w.0, Some("catalog.json"), &recs);
    let want: HashMap<String, String> = [
        ("HTT.0", "normal"),
        ("HTT.1", "pre_mutation"),
        ("HTT.2", "full_mutation"),
        ("SAMD12.3", "normal"),
        ("SAMD12.4", "pre_mutation"),
        ("MUC1.5", ""),
    ]
    .iter()
    .map(|(k, v)| ((*k).to_owned(), (*v).to_owned()))
    .collect();
    assert_eq!(
        got, want,
        "Stranger on the catalog: HTT 26 normal, 27 pre_mutation, 36 full_mutation; SAMD12 1 normal, 2 \
         pre_mutation; MUC1 no status"
    );
    let bundled = stranger(&w.0, None, &recs[..3]);
    assert_ne!(
        bundled.get("HTT.1").map(String::as_str),
        Some("pre_mutation"),
        "control: Stranger's bundled GRCh38 catalog grades HTT 27 otherwise: {bundled:?}"
    );
}

// --- the extra repeats -------------------------------------------------------------------

/// GRCh38 bases [start, end) (0-based) of chrom, from NCBI (GRCh38.p14,
/// as the no-ALT set's chr2 and chr7).
fn grch38(chrom: &str, start: u64, end: u64) -> String {
    let acc = match chrom {
        "chr2" => "NC_000002.12",
        "chr7" => "NC_000007.14",
        other => panic!("no accession for {other}"),
    };
    let url = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&id={acc}&seq_start={}&seq_stop={end}",
        start + 1
    );
    let w = Scratch::new(&format!("ncbi-{chrom}-{start}"));
    let dest = w.0.join("seq.fa");
    download(&url, &dest).unwrap();
    std::fs::read_to_string(&dest)
        .unwrap()
        .lines()
        .skip(1)
        .collect::<String>()
        .to_uppercase()
}

fn region(text: &str) -> (&str, u64, u64) {
    let (chrom, span) = text.split_once(':').unwrap();
    let (a, b) = span.split_once('-').unwrap();
    (chrom, a.parse().unwrap(), b.parse().unwrap())
}

/// Why each region of an extra entry is not a whole run of its motif in
/// GRCh38 (empty when it is).
fn region_errors(entry: &IndexMap<String, Value>, shift: u64) -> Vec<String> {
    let regions: Vec<&str> = match &entry["ReferenceRegion"] {
        Value::Array(a) => a.iter().map(s).collect(),
        v => vec![s(v)],
    };
    let motifs = og_setup::strchive::motifs_of(s(&entry["LocusStructure"]));
    let mut errs = Vec::new();
    for (m, r) in motifs.iter().zip(regions) {
        let (chrom, a, b) = region(r);
        let (a, b) = (a + shift, b + shift);
        let k = m.len();
        let seq = grch38(chrom, a - k as u64, b + k as u64);
        let (left, inner, right) = (&seq[..k], &seq[k..seq.len() - k], &seq[seq.len() - k..]);
        if inner != m.repeat(((b - a) as usize) / k) {
            errs.push(format!("{r}: GRCh38 has {inner}, not {m} repeated"));
        } else if left == *m || right == *m {
            errs.push(format!("{r}: the {m} run goes on past the region ({left}|...|{right})"));
        }
    }
    errs
}

#[test]
fn the_extra_repeats() {
    if !online("the_extra_repeats") {
        return;
    }
    let w = Scratch::new("extra");
    let extra_path = repo().join("assets/str_extra_loci.json").to_string_lossy().into_owned();
    let extra_v = load(&extra_path);
    let extra = list(&extra_v);
    let shape: Vec<(&str, &str)> = extra
        .iter()
        .map(|e| (s(&obj(e)["LocusId"]), s(&obj(e)["LocusStructure"])))
        .collect();
    assert_eq!(
        shape,
        [("CFTR", "(TG)*(T)*"), ("UGT1A1", "(TA)*")],
        "the extra catalog names CFTR (TG and T) and UGT1A1 (TA)"
    );
    let missing = |cat: &Value| -> Vec<String> {
        let have: BTreeSet<&str> = list(cat).iter().map(|e| s(&obj(e)["LocusId"])).collect();
        extra
            .iter()
            .map(|e| s(&obj(e)["LocusId"]))
            .filter(|i| !have.contains(i))
            .map(str::to_owned)
            .collect()
    };

    let (out, rng) = (w.path("catalog.json"), w.path("ranges.tsv"));
    let rc = catalog(&["--out", &out, "--ranges", &rng, "--extra", &extra_path]);
    assert_eq!(rc, 0, "--extra: the converter runs");
    let cat_v = load(&out);
    let cat = list(&cat_v);
    assert_eq!(cat.len(), 82 + extra.len(), "--extra: STRchive's 82 entries and the extra loci");
    assert!(missing(&cat_v).is_empty(), "every extra locus is in the merged catalog");
    let tail = &cat[cat.len() - extra.len()..];
    assert_eq!(tail, extra, "the extra loci are copied as they are, last");
    assert!(
        tail.iter().map(obj).all(|e| !e.contains_key("NormalMax")
            && !e.contains_key("PathologicMin")
            && e.get("NoThresholds").is_some_and(|v| !s(v).is_empty())),
        "no extra locus has Stranger thresholds, each says why"
    );
    let ranges = std::fs::read_to_string(&rng).unwrap();
    let rows: Vec<Vec<&str>> = ranges.lines().skip(2).map(|l| l.split('\t').collect()).collect();
    let want_tail: Vec<Vec<&str>> = extra
        .iter()
        .map(|e| vec![s(&obj(e)["LocusId"]), s(&obj(e)["SourceId"]), "", "", "", ""])
        .collect();
    let got_tail: Vec<Vec<&str>> = rows[rows.len() - extra.len()..].iter().map(|r| r[..6].to_vec()).collect();
    assert!(
        rows.len() == cat.len() && got_tail == want_tail,
        "the ranges table has a row per entry, the extra ones without STRchive ranges: {got_tail:?}"
    );

    // Control: an extra catalog without UGT1A1 gives a merged catalog the
    // missing-locus check flags.
    let short = w.path("short.json");
    dump(
        &short,
        &Value::Array(extra.iter().filter(|e| s(&obj(e)["LocusId"]) != "UGT1A1").cloned().collect()),
    );
    catalog(&["--out", &out, "--extra", &short]);
    assert_eq!(
        missing(&load(&out)),
        ["UGT1A1"],
        "control: an extra catalog without UGT1A1 gives a merged catalog the check flags"
    );

    // Controls: four planted bad entries each stop the converter.
    let cftr = obj(&extra[0]).clone();
    let with = |k: &str, v: Value| {
        let mut e = cftr.clone();
        e.insert(k.to_owned(), v);
        Value::Object(e)
    };
    let planted = [
        ("no NoThresholds", with("NoThresholds", Value::Str(String::new()))),
        (
            "a region that is not a whole number of motifs",
            with(
                "ReferenceRegion",
                Value::Array(vec![
                    Value::Str("chr7:117548606-117548627".into()),
                    Value::Str("chr7:117548627-117548635".into()),
                ]),
            ),
        ),
        ("a STRchive locus name", with("LocusId", Value::Str("HTT".into()))),
        ("Stranger thresholds", {
            let mut e = cftr.clone();
            e.insert("NormalMax".into(), Value::Int(11));
            e.insert("PathologicMin".into(), Value::Int(13));
            Value::Object(e)
        }),
    ];
    for (desc, bad) in planted {
        let (p, o, g) = (w.path("bad.json"), w.path("bad_out.json"), w.path("bad_ranges.tsv"));
        dump(&p, &Value::Array(vec![bad]));
        let rc = catalog(&["--out", &o, "--ranges", &g, "--extra", &p]);
        assert!(
            rc == 2 && !Path::new(&o).exists() && !Path::new(&g).exists(),
            "control: an extra entry with {desc} stops the converter (exit 2, nothing written): exit {rc}"
        );
    }

    for e in extra.iter().map(obj) {
        let errs = region_errors(e, 0);
        assert!(errs.is_empty(), "{}: every region is a whole run of its motif in GRCh38: {errs:?}", s(&e["LocusId"]));
    }
    assert!(
        !region_errors(obj(&extra[1]), 1).is_empty(),
        "control: the UGT1A1 region shifted by one base fails the sequence check"
    );
}
