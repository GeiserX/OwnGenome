//! Public slices of real files (tests/data/README.md) read, parsed and
//! written back byte for byte, through the sequential and the parallel
//! readers and the parallel writer.

use std::io::{BufRead, Read, Write};
use std::num::NonZero;
use std::path::PathBuf;

use og_io::bgzf;
use og_vcf::{CsqLayout, Gt, Header, Lines, Record, csq_field, write_info};

const SLICES: [&str; 4] = [
    "hg002_vep.vcf.gz",
    "hg002_benchmark.vcf.gz",
    "hg002_smvar_trf.vcf.gz",
    "hg002_chr20.g.vcf.gz",
];

fn path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn workers(n: usize) -> NonZero<usize> {
    NonZero::new(n).unwrap()
}

fn read_all(mut r: impl Read) -> Vec<u8> {
    let mut out = Vec::new();
    r.read_to_end(&mut out).unwrap();
    out
}

/// How a test rebuilds INFO from its items.
type Rebuild = fn(&Record, &mut Vec<u8>);

fn in_order(rec: &Record, out: &mut Vec<u8>) {
    write_info(rec.info_items(), out);
}

/// The planted fault: INFO items written back in another order.
fn reordered(rec: &Record, out: &mut Vec<u8>) {
    let mut items: Vec<_> = rec.info_items().collect();
    items.sort();
    write_info(items, out);
}

/// Reads a whole VCF and writes it back, header and records, rebuilding each
/// record's INFO with `rebuild`. Returns the output and the record count.
fn rewrite(input: impl BufRead, rebuild: Rebuild) -> (Vec<u8>, usize) {
    let mut lines = Lines::new(input);
    let header = Header::read(&mut lines).unwrap();
    let mut out = Vec::new();
    header.write(&mut out, &[]).unwrap();
    let mut info = Vec::new();
    let mut records = 0;
    while let Some(line) = lines.next_line().unwrap() {
        let rec = Record::parse(line).unwrap();
        info.clear();
        rebuild(&rec, &mut info);
        rec.write_with_info(&info, &mut out);
        out.push(b'\n');
        records += 1;
    }
    (out, records)
}

#[test]
fn slices_round_trip_byte_for_byte() {
    for name in SLICES {
        let want = read_all(og_io::open(path(name)).unwrap());
        assert!(want.len() > 5000, "{name}: too small to mean anything");
        for n in [1, 3, 8] {
            let got = read_all(bgzf::open(path(name), workers(n)).unwrap());
            assert!(
                got == want,
                "{name}: the parallel reader with {n} workers differs"
            );
        }
        let (out, records) = rewrite(bgzf::open(path(name), workers(4)).unwrap(), in_order);
        assert!(records >= 12, "{name}: {records} records");
        assert!(out == want, "{name}: written back differently");
    }
}

#[test]
fn a_planted_info_reorder_is_caught() {
    // The benchmark slice has 8 to 10 INFO items per record, never sorted.
    let input = og_io::open(path("hg002_benchmark.vcf.gz")).unwrap();
    let want = read_all(og_io::open(path("hg002_benchmark.vcf.gz")).unwrap());
    let (out, _) = rewrite(input, reordered);
    assert!(
        out != want,
        "reordering INFO went unnoticed: the round trip proves nothing"
    );
}

#[test]
fn the_parallel_writer_round_trips() {
    for name in SLICES {
        let want = read_all(og_io::open(path(name)).unwrap());
        let mut w = bgzf::Writer::new(Vec::new(), workers(4)).unwrap();
        w.write_all(&want).unwrap();
        let file = w.finish().unwrap();
        assert!(file.ends_with(&bgzf::EOF_BLOCK));
        assert!(read_all(og_io::reader(std::io::Cursor::new(file.clone())).unwrap()) == want);
        assert!(read_all(bgzf::reader(std::io::Cursor::new(file), workers(2)).unwrap()) == want);
    }
}

#[test]
fn every_genotype_and_csq_entry_reads() {
    for name in SLICES {
        let mut lines = Lines::new(og_io::open(path(name)).unwrap());
        let header = Header::read(&mut lines).unwrap();
        let samples = header.samples().count();
        assert_eq!(samples, 1, "{name}");
        assert!(header.contigs().any(|c| c == b"chr20"), "{name}");
        let csq = header.csq_layout();
        while let Some(line) = lines.next_line().unwrap() {
            let rec = Record::parse(line).unwrap();
            assert!(rec.pos().unwrap() > 0);
            let gt: Gt = rec
                .gt(0)
                .unwrap_or_else(|| panic!("{name}: no GT in {}", String::from_utf8_lossy(line)));
            assert!(!gt.alleles.is_empty());
            if let Some(layout) = &csq {
                check_csq(layout, &rec);
            }
        }
    }
}

fn check_csq(layout: &CsqLayout, rec: &Record) {
    let csq = rec.info_get(b"CSQ").flatten().expect("a CSQ value");
    for entry in csq.split(|&c| c == b',') {
        assert_eq!(entry.split(|&c| c == b'|').count(), layout.names.len());
        let impact = csq_field(entry, layout.impact.unwrap());
        assert!(
            matches!(impact, b"HIGH" | b"MODERATE" | b"LOW" | b"MODIFIER"),
            "{impact:?}"
        );
        assert!(!csq_field(entry, layout.consequence.unwrap()).is_empty());
    }
}

#[test]
fn the_vep_slice_has_the_fields_ranking_needs() {
    let mut lines = Lines::new(og_io::open(path("hg002_vep.vcf.gz")).unwrap());
    let layout = Header::read(&mut lines).unwrap().csq_layout().unwrap();
    for field in [
        layout.impact,
        layout.consequence,
        layout.symbol,
        layout.canonical,
        layout.mane,
    ] {
        assert!(field.is_some());
    }
    // This VEP run had no frequency fields: the layout must say so.
    let has_freq = layout
        .names
        .iter()
        .any(|n| matches!(n.as_str(), "MAX_AF" | "gnomADe_AF" | "gnomADg_AF"));
    assert_eq!(layout.freq.is_empty(), !has_freq);
    assert!(layout.max_needed <= layout.names.len());
}
