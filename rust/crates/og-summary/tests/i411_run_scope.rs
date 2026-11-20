//! What part of the genome a run covered (`scope` in summary.json), and the
//! partial mark on each genome-wide section of a run that covered less.
//!
//! A run on a part of the genome (`--intervals`, a sliced VCF or BAM) gives
//! ROH, archaic, ancestry, PRS, callability and telomere results that read
//! as whole-genome results. The collector reads the contigs of the sample
//! VCF's records, writes them as `scope`, and marks each of those sections
//! "partial: covers N of 22 autosomes" when fewer than 22 are there.

use og_summary::collect::{Group, Input, post};
use pycompat::json::Value;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

/// The sections whose result is a share or a sum over the whole genome.
const GENOME_WIDE: [&str; 6] = [
    "callability",
    "prs",
    "ancestry_fine",
    "telomere",
    "roh",
    "archaic",
];

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!("og-run-scope-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn path(&self) -> String {
        self.0.to_str().unwrap().to_owned()
    }
    fn put(&self, rel: &str, body: &[u8]) {
        let p = self.0.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(&p, body).unwrap();
    }
    fn put_gz(&self, rel: &str, text: &str) {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        e.write_all(text.as_bytes()).unwrap();
        self.put(rel, &e.finish().unwrap());
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn collect(dir: &Dir) -> Value {
    post::collect(&Input {
        sample: "S".into(),
        sample_dir: dir.path(),
        declared_sex: None,
        assets: format!("{REPO}/assets"),
        group: Group::All,
        tables: format!("{REPO}/bin"),
    })
    .unwrap()
}

fn get<'a>(v: &'a Value, key: &str) -> Option<&'a Value> {
    match v {
        Value::Object(m) => m.get(key),
        _ => None,
    }
}

fn at<'a>(v: &'a Value, path: &[&str]) -> &'a Value {
    path.iter().fold(v, |v, k| {
        get(v, k).unwrap_or_else(|| panic!("no {k} in {v:?}"))
    })
}

fn st(v: &Value) -> &str {
    match v {
        Value::Str(s) => s,
        _ => panic!("not a string: {v:?}"),
    }
}

/// A sample folder whose sample VCF has records on `contigs`, with a result
/// in each genome-wide section, as the steps write them.
fn folder(tag: &str, contigs: &[&str]) -> Dir {
    let d = Dir::new(tag);
    let mut vcf = String::from(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
    );
    for (i, c) in contigs.iter().enumerate() {
        for pos in [1_000_000, 2_000_000] {
            vcf.push_str(&format!(
                "{c}\t{}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
                pos + i
            ));
        }
    }
    d.put_gz("vcf/S.vcf.gz", &vcf);
    // Step 11: bcftools roh's sites and one run on each chromosome called.
    let mut roh = String::from("# ST\n# RG\n");
    for c in contigs {
        roh.push_str(&format!("ST\tS\t{c}\t1000000\t1\t50\n"));
    }
    roh.push_str(&format!(
        "RG\tS\t{}\t1000000\t7000000\t6000001\t9\t50\n",
        contigs[0]
    ));
    d.put("vcf/S_roh.txt", roh.as_bytes());
    // Step 16b: the catalog's genes on the chromosomes the BAM has.
    let catalog = fs::read_to_string(format!("{REPO}/assets/callability_genes.tsv")).unwrap();
    let mut bed = String::from("#chrom\tstart\tend\tregion\t10X\n");
    for line in catalog.lines().filter(|l| !l.starts_with('#')).skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        let plain = f[2].trim_start_matches("chr");
        if f[3] != "." && contigs.iter().any(|c| c.trim_start_matches("chr") == plain) {
            let (a, b): (i64, i64) = (f[3].parse().unwrap(), f[4].parse().unwrap());
            bed.push_str(&format!("{}\t{a}\t{b}\t{}\t{}\n", f[2], f[0], b - a));
        }
    }
    d.put_gz("mosdepth/S_callability.bed.gz", &bed);
    d.put(
        "prs/S_prs_summary.tsv",
        b"Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\tPercentile\n\
          Coronary artery disease\tPGS000018\t0.5\t10\t100\tNA\n",
    );
    d.put(
        "ancestry_fine/S_ancestry_fine.tsv",
        b"level\tcode\tgroup\tfraction\tse\tlow\thigh\ncontinent\tEUR\tEurope\t1.0\t0.01\t0.98\t1.0\n",
    );
    d.put(
        "telomere/S/S/S_summary.tsv",
        b"PID\tsample\ttotal_reads\tread_lengths\tintratel_reads\ttel_content\nS\ttumor\t1000\t150\t10\t12.5\n",
    );
    d.put(
        "archaic/S_archaic.tsv",
        b"key\tvalue\npopulation\tEUR\nneanderthal_segments\t3\nneanderthal_bp\t300000\n\
          denisovan_segments\t0\ndenisovan_bp\t0\nanalysed_bp\t1000000\n\
          neanderthal_percent\t1.5\ndenisovan_percent\t0.0\n",
    );
    d
}

fn autosomes() -> Vec<String> {
    (1..=22).map(|i| format!("chr{i}")).collect()
}

#[test]
fn a_run_on_chr20_and_chr22_is_partial_in_every_genome_wide_section() {
    let d = folder("two", &["chr20", "chr22"]);
    let v = collect(&d);
    let scope = at(&v, &["scope"]);
    assert_eq!(at(scope, &["autosomes_covered"]), &Value::Int(2));
    assert_eq!(
        at(scope, &["autosomes"]),
        &Value::Array(vec![Value::Str("chr20".into()), Value::Str("chr22".into())])
    );
    assert_eq!(at(scope, &["chrx"]), &Value::Bool(false));
    assert_eq!(at(scope, &["intervals"]), &Value::Null);
    assert_eq!(at(scope, &["partial"]), &Value::Bool(true));
    assert_eq!(
        st(at(scope, &["mark"])),
        "partial: covers 2 of 22 autosomes"
    );
    for key in GENOME_WIDE {
        let sec = at(&v, &["sections", key]);
        assert_eq!(st(at(sec, &["state"])), "ok", "{key}: {sec:?}");
        assert_eq!(
            st(at(sec, &["partial"])),
            "partial: covers 2 of 22 autosomes",
            "{key}"
        );
    }
    // A section about the calls themselves, not a share of the genome.
    let variants = at(&v, &["sections", "variants"]);
    assert_eq!(st(at(variants, &["state"])), "ok");
    assert!(get(variants, "partial").is_none());
}

#[test]
fn a_whole_genome_run_carries_no_partial_mark() {
    let mut names = autosomes();
    names.push("chrX".into());
    let plain: Vec<String> = names.iter().map(|c| c[3..].to_owned()).collect();
    for (tag, names) in [("chr", &names), ("plain", &plain)] {
        let contigs: Vec<&str> = names.iter().map(String::as_str).collect();
        let d = folder(tag, &contigs);
        let v = collect(&d);
        let scope = at(&v, &["scope"]);
        assert_eq!(at(scope, &["autosomes_covered"]), &Value::Int(22), "{tag}");
        assert_eq!(
            at(scope, &["autosomes"]),
            &Value::Array(autosomes().into_iter().map(Value::Str).collect()),
            "{tag}"
        );
        assert_eq!(at(scope, &["chrx"]), &Value::Bool(true), "{tag}");
        assert_eq!(at(scope, &["partial"]), &Value::Bool(false), "{tag}");
        assert_eq!(at(scope, &["mark"]), &Value::Null, "{tag}");
        for key in GENOME_WIDE {
            let sec = at(&v, &["sections", key]);
            assert_eq!(st(at(sec, &["state"])), "ok", "{tag} {key}: {sec:?}");
            assert!(get(sec, "partial").is_none(), "{tag} {key}");
        }
    }
}

#[test]
fn intervals_the_run_manifest_records_make_a_whole_genome_run_partial() {
    let names = autosomes();
    let contigs: Vec<&str> = names.iter().map(String::as_str).collect();
    let d = folder("intervals", &contigs);
    d.put(
        "run_manifest.tsv",
        b"run\tintervals\tchr1:1-1000000 chr2:1-1000000\n",
    );
    let v = collect(&d);
    let scope = at(&v, &["scope"]);
    assert_eq!(at(scope, &["autosomes_covered"]), &Value::Int(22));
    assert_eq!(
        st(at(scope, &["intervals"])),
        "chr1:1-1000000 chr2:1-1000000"
    );
    assert_eq!(at(scope, &["partial"]), &Value::Bool(true));
    let mark = "partial: limited to the run's intervals (chr1:1-1000000 chr2:1-1000000)";
    assert_eq!(st(at(scope, &["mark"])), mark);
    for key in GENOME_WIDE {
        assert_eq!(st(at(&v, &["sections", key, "partial"])), mark, "{key}");
    }
}

#[test]
fn without_a_readable_sample_vcf_no_scope_is_written_and_nothing_is_marked() {
    let d = folder("novcf", &["chr20"]);
    fs::remove_file(d.0.join("vcf/S.vcf.gz")).unwrap();
    let v = collect(&d);
    assert!(get(&v, "scope").is_none());
    for key in GENOME_WIDE {
        assert!(
            get(at(&v, &["sections", key]), "partial").is_none(),
            "{key}"
        );
    }
    // A VCF cut off mid-stream is no scope either: its contigs are unknown.
    let d = folder("badvcf", &["chr20"]);
    d.put("vcf/S.vcf.gz", b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03");
    let v = collect(&d);
    assert!(get(&v, "scope").is_none());
}

/// Marks the ROH result (step 11) as from an earlier run: the run started
/// after every file was written, and step 11 failed in it.
fn roh_from_an_earlier_run(d: &Dir) {
    let started = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs()
        + 86_400;
    let mut status = format!("meta\tstarted_epoch\t{started}\n");
    for step in og_summary::collect::model::RUN_ALL_STEPS {
        let result = if *step == "11" {
            "failed (exit 1)"
        } else {
            "ok"
        };
        status.push_str(&format!("step\t{step}\t{result}\n"));
    }
    d.put("logs/run_status.tsv", status.as_bytes());
}

#[test]
fn a_stale_genome_wide_result_never_reads_as_complete() {
    // A whole-genome run whose ROH step failed shows the ROH result of an
    // earlier run. That run's scope was not recorded: it may have been a
    // chr22-only run, so the mark says the scope is not known.
    let names = autosomes();
    let contigs: Vec<&str> = names.iter().map(String::as_str).collect();
    for (tag, contigs, mark) in [
        ("stale-whole", contigs.clone(), None),
        (
            "stale-two",
            vec!["chr20", "chr22"],
            Some("partial: covers 2 of 22 autosomes"),
        ),
    ] {
        let d = folder(tag, &contigs);
        roh_from_an_earlier_run(&d);
        let v = collect(&d);
        let roh = at(&v, &["sections", "roh"]);
        assert_eq!(st(at(roh, &["state"])), "stale", "{tag}: {roh:?}");
        assert_eq!(
            st(at(roh, &["partial"])),
            "partial: not known, from an earlier run",
            "{tag}"
        );
        for key in GENOME_WIDE.iter().filter(|k| **k != "roh") {
            let sec = at(&v, &["sections", key]);
            assert_eq!(st(at(sec, &["state"])), "ok", "{tag} {key}: {sec:?}");
            assert_eq!(get(sec, "partial").map(st), mark, "{tag} {key}");
        }
    }
}
