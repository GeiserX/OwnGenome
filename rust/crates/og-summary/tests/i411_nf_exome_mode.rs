//! An exome's mean depth (`--model_type WES` on the Nextflow path).
//!
//! mosdepth's summary has a `total` row, the reads spread over the whole
//! genome, and with `--by <capture BED>` a `total_region` row, the mean over
//! the capture regions. A 100x exome reads about 2x in `total`. The report's
//! mean depth (the QC lines, the live headline, the genome facts' depth
//! floor) must be the on-target one for an exome, and stay `total` for a
//! whole genome, whose 500 bp windows also give a `total_region` row.
//! The run reads as an exome when the run manifest says so (HTML_REPORT
//! writes `run model_type WES`), or when the summary itself shows a capture:
//! a `total_region` row under half as long as the `total` row. 500 bp
//! windows tile every contig mosdepth reads, so their `total_region` is as
//! long as `total`; a capture kit covers 1-3% of the genome. The second rule
//! covers the paths that write no manifest row: run-all.sh with
//! `--model_type WES`, a re-render with 24-html-report.sh, og live's own
//! re-read of a sample folder, and step 16b with `CAPTURE_BED`.

use og_summary::collect::{Group, Input, post};
use pycompat::json::Value;
use std::fs;
use std::path::PathBuf;

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

/// A mosdepth summary of an exome: `--by <capture BED>`, so `total_region`
/// is the capture's 1,250,000 bp of chr20's 64,444,167 (about 2%).
const SUMMARY: &str = "chrom\tlength\tbases\tmean\tmin\tmax\n\
chr20\t64444167\t135332751\t2.10\t0\t900\n\
chr20_region\t1250000\t123750000\t99.00\t0\t900\n\
total\t64444167\t135332751\t2.10\t0\t900\n\
total_region\t1250000\t123750000\t99.00\t0\t900\n";

/// A mosdepth summary of a whole genome: `--by 500`, so `total_region` is
/// as long as `total` (mosdepth 0.3.14 writes the same length in both).
/// Real windows give the same mean in both rows too; the two means differ
/// here only so the test sees which row is read.
const WINDOWS: &str = "chrom\tlength\tbases\tmean\tmin\tmax\n\
chr20\t64444167\t1933325010\t30.00\t0\t900\n\
chr20_region\t64444167\t1939769427\t30.10\t0\t900\n\
total\t64444167\t1933325010\t30.00\t0\t900\n\
total_region\t64444167\t1939769427\t30.10\t0\t900\n";

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
fn an_exome_without_the_manifest_row_reads_the_on_target_mean_depth() {
    // The paths that write no model_type row: run-all.sh rewrites the run
    // manifest with bin/write_manifest.sh, as does a re-render with
    // 24-html-report.sh; og live re-reads the folder; step 16b with
    // CAPTURE_BED has no manifest of its own. The capture-length
    // `total_region` says exome.
    for (tag, manifest) in [
        ("wes-no-row", "run\ttools\tmosdepth\n"),
        ("wes-run-all", "run\tscript\trun-all.sh\n"),
        ("wes-no-manifest", ""),
    ] {
        let d = folder(tag, SUMMARY, manifest);
        assert_eq!(mean_depth(&d), Value::Float(99.0), "{tag}");
    }
}

#[test]
fn a_whole_genome_keeps_the_genome_wide_mean_depth() {
    // 500 bp windows: `total_region` as long as `total`, with no model_type
    // row, a WGS row, or no manifest: all read `total`.
    for (tag, manifest) in [
        ("wgs-none", "run\ttools\tmosdepth\n"),
        ("wgs-row", "run\tmodel_type\tWGS\n"),
        ("wgs-no-manifest", ""),
    ] {
        let d = folder(tag, WINDOWS, manifest);
        assert_eq!(mean_depth(&d), Value::Float(30.0), "{tag}");
    }
    // Regions half as long as the genome or more are not a capture; and a
    // summary without --by has `total` only.
    let half = "chrom\tlength\tbases\tmean\tmin\tmax\n\
total\t64444167\t1933325010\t30.00\t0\t900\n\
total_region\t32222084\t969864728\t30.10\t0\t900\n";
    let no_by = "chrom\tlength\tbases\tmean\tmin\tmax\n\
total\t64444167\t1933325010\t30.00\t0\t900\n";
    for (tag, summary) in [("wgs-half", half), ("wgs-no-by", no_by)] {
        let d = folder(tag, summary, "");
        assert_eq!(mean_depth(&d), Value::Float(30.0), "{tag}");
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
