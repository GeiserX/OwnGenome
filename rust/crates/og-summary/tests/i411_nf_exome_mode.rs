//! An exome's mean depth (`--model_type WES` on the Nextflow path).
//!
//! mosdepth's summary has a `total` row, the reads spread over the whole
//! genome, and with `--by <capture BED>` a `total_region` row, the mean over
//! the capture regions. A 100x exome reads about 2x in `total`. The report's
//! mean depth (the QC lines, the live headline, the genome facts' depth
//! floor) must be the on-target one for an exome, and stay `total` for a
//! whole genome, whose 500 bp windows also give a `total_region` row.
//! HTML_REPORT writes `run model_type WES` into the run manifest for an
//! exome; nothing is written for a whole genome.

use og_summary::collect::{Group, Input, post};
use pycompat::json::Value;
use std::fs;
use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

/// A mosdepth summary where the genome-wide and on-target means differ.
const SUMMARY: &str = "chrom\tlength\tbases\tmean\tmin\tmax\n\
chr20\t64444167\t135332751\t2.10\t0\t900\n\
chr20_region\t1250000\t123750000\t99.00\t0\t900\n\
total\t64444167\t135332751\t2.10\t0\t900\n\
total_region\t1250000\t123750000\t99.00\t0\t900\n";

struct Dir(PathBuf);

impl Dir {
    fn new(tag: &str) -> Dir {
        let d = std::env::temp_dir().join(format!("og-nf-exome-{}-{tag}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        Dir(d)
    }
    fn put(&self, rel: &str, body: &str) {
        let p = self.0.join(rel);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(&p, body).unwrap();
    }
}

impl Drop for Dir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// A Nextflow sample folder with `summary` as mosdepth's summary and
/// `manifest` as its run manifest (none when empty).
fn folder(tag: &str, summary: &str, manifest: &str) -> Dir {
    let d = Dir::new(tag);
    d.put("coverage/S.mosdepth.summary.txt", summary);
    if !manifest.is_empty() {
        d.put("run_manifest.tsv", manifest);
    }
    d
}

fn mean_depth(d: &Dir) -> Value {
    let v = post::collect(&Input {
        sample: "S".into(),
        sample_dir: d.0.to_str().unwrap().to_owned(),
        declared_sex: None,
        assets: format!("{REPO}/assets"),
        group: Group::All,
        tables: format!("{REPO}/bin"),
    })
    .unwrap();
    ["sections", "coverage", "data", "mean_depth"]
        .iter()
        .fold(&v, |v, k| match v {
            Value::Object(m) => m.get(*k).unwrap_or_else(|| panic!("no {k} in {v:?}")),
            _ => panic!("no {k} in {v:?}"),
        })
        .clone()
}

#[test]
fn an_exome_reads_the_on_target_mean_depth() {
    let d = folder(
        "wes",
        SUMMARY,
        "run\ttools\tmosdepth\nrun\tmodel_type\tWES\n",
    );
    assert_eq!(mean_depth(&d), Value::Float(99.0));
}

#[test]
fn a_whole_genome_keeps_the_genome_wide_mean_depth() {
    // No model_type row (a whole genome, or the bash runner), and the row
    // a WGS model_type would give if one were written: both read `total`.
    for (tag, manifest) in [
        ("wgs-none", "run\ttools\tmosdepth\n"),
        ("wgs-row", "run\tmodel_type\tWGS\n"),
        ("wgs-no-manifest", ""),
    ] {
        let d = folder(tag, SUMMARY, manifest);
        assert_eq!(mean_depth(&d), Value::Float(2.1), "{tag}");
    }
}

#[test]
fn an_exome_summary_without_capture_regions_has_no_mean_depth() {
    // mosdepth run without --by on an exome: the genome-wide `total` is the
    // only mean, and it is not the exome's depth.
    let summary = "chrom\tlength\tbases\tmean\tmin\tmax\n\
chr20\t64444167\t135332751\t2.10\t0\t900\n\
total\t64444167\t135332751\t2.10\t0\t900\n";
    let d = folder("wes-no-region", summary, "run\tmodel_type\tWES\n");
    assert_eq!(mean_depth(&d), Value::Null);
}
