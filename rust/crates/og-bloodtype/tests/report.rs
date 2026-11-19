//! `og bloodtype report` end to end on the repository's marker table and
//! catalog, with synthetic VCFs and depth tables.

use og_bloodtype::genotypes::{Access, read_vcf};
use og_bloodtype::markers::{read_markers, rhd_regions, rhd_row};
use og_bloodtype::{Files, check_markers, report};
use og_cn::catalog::Catalog;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");

fn markers() -> String {
    format!("{REPO}/assets/bloodtype_markers.tsv")
}

fn cn_loci() -> String {
    format!("{REPO}/bin/cn_loci.tsv")
}

/// A fresh folder for one test.
fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-bloodtype-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

const HEAD: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

/// A one-sample VCF record at the marker `name`.
fn rec(name: &str, gt: &str) -> String {
    let m = &read_markers(&markers()).unwrap()[name];
    format!(
        "{}\t{}\t.\t{}\t{}\t50\tPASS\t.\tGT\t{gt}\n",
        m.chrom, m.pos, m.ref_, m.alt
    )
}

/// A bedcov table: every site at 30x, RHD at `copies` with exon `e`
/// at `exons[e]` copies.
fn depth(d: &Path, copies: f64, exons: &HashMap<i64, f64>) -> String {
    let cat = Catalog::load(&cn_loci()).unwrap();
    let mut text = String::new();
    for m in read_markers(&markers()).unwrap().values() {
        let len = m.pos - m.start;
        text.push_str(&format!(
            "{}\t{}\t{}\tsite\t{}\t{}\t{}\n",
            m.chrom,
            m.start,
            m.pos,
            m.system,
            m.name,
            30 * len
        ));
    }
    for (name, r) in rhd_regions(rhd_row(&cat).unwrap()) {
        let c = match name.as_str() {
            "RHD_FLANK_UP" | "RHD_FLANK_DOWN" => 2.0,
            "RHD_GENE" => copies,
            n => exons
                .get(&n["RHD_EX".len()..].parse::<i64>().unwrap())
                .copied()
                .unwrap_or(copies),
        };
        let mean = 15.0 * c;
        text.push_str(&format!(
            "{}\t{}\t{}\t{name}\t{}\n",
            r.chrom,
            r.start,
            r.end,
            mean * (r.end - r.start) as f64
        ));
    }
    let p = d.join("d.tsv");
    std::fs::write(&p, text).unwrap();
    p.to_string_lossy().into_owned()
}

/// {system: (phenotype, genotype, confidence, note)} of one report.
fn run(d: &Path, records: &[String], depth: Option<String>) -> HashMap<String, Vec<String>> {
    let vcf = d.join("s.vcf");
    std::fs::write(&vcf, format!("{HEAD}{}", records.concat())).unwrap();
    let files = Files {
        vcf: vcf.to_string_lossy().into_owned(),
        markers: markers(),
        depth,
        cn_loci: Some(cn_loci()),
        out: d.join("t.tsv").to_string_lossy().into_owned(),
        evidence: d.join("e.tsv").to_string_lossy().into_owned(),
    };
    report(&files).unwrap();
    std::fs::read_to_string(&files.out)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<String> = l.split('\t').map(str::to_owned).collect();
            (c[0].clone(), c[1..].to_vec())
        })
        .collect()
}

#[test]
fn the_shipped_table_agrees_with_the_rules() {
    assert_eq!(
        check_markers(&markers(), Some(&cn_loci())).unwrap(),
        Vec::<String>::new()
    );
}

#[test]
fn a_reference_genome_with_reads() {
    let d = dir("ref");
    let t = run(&d, &[], Some(depth(&d, 2.0, &HashMap::new())));
    for (system, want) in [
        ("ABO", "O"),
        ("RhD", "D+"),
        ("RhCE", "C-c+ E-e+"),
        ("Kell", "K-k+ Kp(a-b+) Js(a-b+)"),
        ("Duffy", "Fy(a+b-)"),
        ("Kidd", "Jk(a+b-)"),
        ("HPA-15", "HPA-15(a-b+)"),
        ("H", "H+"),
    ] {
        assert_eq!(t[system][0], want, "{system}");
    }
    assert_eq!(t["ABO"][2], "high");
    assert_eq!(t["RhCE"][2], "low");
}

#[test]
fn abo_a1b_and_two_kidd_nulls() {
    let d = dir("abo");
    let recs = [
        rec("ABO_261", "1/1"),
        rec("ABO_803", "0/1"),
        rec("ABO_796", "0/1"),
        rec("JK_01N06", "1/1"),
    ];
    let t = run(&d, &recs, Some(depth(&d, 2.0, &HashMap::new())));
    assert_eq!(t["ABO"][..2], ["A1B", "A1/B"]);
    assert_eq!(t["Kidd"][..2], ["Jk(a-b-)", "null/null"]);
}

#[test]
fn rhd_hybrid_needs_an_exon_under_half_a_copy() {
    let d = dir("hybrid");
    let exons = |c: f64| (3..=7).map(|e| (e, c)).collect::<HashMap<i64, f64>>();
    let t = run(&d, &[], Some(depth(&d, 1.6, &exons(0.45))));
    assert_eq!(t["RhD"][0], "not determined");
    assert!(t["RhD"][3].contains("exons 3, 4, 5, 6, 7 without RHD reads"));
    let t = run(&d, &[], Some(depth(&d, 1.6, &exons(0.55))));
    assert_eq!(t["RhD"][0], "D+");
}

#[test]
fn rhd_deleted_on_both_chromosomes() {
    let d = dir("rhd0");
    let t = run(&d, &[], Some(depth(&d, 0.0, &HashMap::new())));
    assert_eq!(
        t["RhD"][..3],
        [
            "D-",
            "RHD deleted/RHD deleted (RHD copies by read depth 0.00)",
            "high"
        ]
    );
}

/// Copies `from` to `to` as a new file (a fresh mtime, unlike `fs::copy`).
fn put(from: PathBuf, to: impl AsRef<Path>) {
    std::fs::write(to, std::fs::read(from).unwrap()).unwrap();
}

/// The parity inputs: one VCF bgzipped with its .tbi, and the same VCF plain.
fn parity(name: &str) -> PathBuf {
    PathBuf::from(format!("{REPO}/tests/parity/inputs/bloodtype/{name}"))
}

#[test]
fn the_index_and_a_full_pass_read_the_same_sites() {
    for stem in ["abo", "blocks", "gvcf", "vcf-quirks", "nochr", "bombay"] {
        let d = dir(&format!("idx-{stem}"));
        let gz = d.join(format!("{stem}.vcf.gz"));
        put(parity(&format!("{stem}.vcf.gz")), &gz);
        put(
            parity(&format!("{stem}.vcf.gz.tbi")),
            d.join(format!("{stem}.vcf.gz.tbi")),
        );
        let (indexed, how) =
            read_vcf(&gz.to_string_lossy(), read_markers(&markers()).unwrap()).unwrap();
        assert_eq!(how, Access::Index, "{stem}");
        let plain = parity(&format!("{stem}.vcf"));
        let (whole, how) =
            read_vcf(&plain.to_string_lossy(), read_markers(&markers()).unwrap()).unwrap();
        assert_eq!(how, Access::Stream, "{stem}");
        assert_eq!(indexed, whole, "{stem}");
        let called = whole
            .values()
            .filter(|s| s.status != og_bloodtype::genotypes::Status::Absent)
            .count();
        assert!(called > 0, "{stem}: no site read");
    }
}

#[test]
fn an_index_older_than_its_vcf_is_not_used() {
    let d = dir("stale");
    let gz = d.join("abo.vcf.gz");
    put(parity("abo.vcf.gz.tbi"), d.join("abo.vcf.gz.tbi"));
    std::thread::sleep(std::time::Duration::from_millis(20));
    put(parity("abo.vcf.gz"), &gz);
    let (_, how) = read_vcf(&gz.to_string_lossy(), read_markers(&markers()).unwrap()).unwrap();
    assert_eq!(how, Access::Stream);
}

/// A one-sample VCF record line.
fn line(chrom: &str, pos: i64, r: &str, alt: &str, filt: &str, gt: &str) -> String {
    format!("{chrom}\t{pos}\t.\t{r}\t{alt}\t50\t{filt}\t.\tGT\t{gt}\n")
}

/// `records` (sorted, under 64 KiB) as a plain VCF and as a bgzipped VCF
/// with a `.tbi`: one bin per sequence whose one chunk spans that
/// sequence's lines, the smallest index the reader takes.
fn both_ways(d: &Path, records: &[String]) -> (PathBuf, PathBuf) {
    use std::io::Write;
    let text = format!("{HEAD}{}", records.concat());
    let plain = d.join("s.vcf");
    std::fs::write(&plain, &text).unwrap();
    let gz = d.join("s.vcf.gz");
    let mut w = og_io::bgzf::Writer::new(
        std::fs::File::create(&gz).unwrap(),
        std::num::NonZero::new(1).unwrap(),
    )
    .unwrap();
    w.write_all(text.as_bytes()).unwrap();
    w.finish().unwrap();
    // Each sequence's (first byte, end byte) in the one block's text.
    let mut seqs: Vec<(String, u64, u64)> = Vec::new();
    let mut at = HEAD.len() as u64;
    for r in records {
        let chrom = r.split('\t').next().unwrap().to_owned();
        let end = at + r.len() as u64;
        match seqs.last_mut() {
            Some(s) if s.0 == chrom => s.2 = end,
            _ => seqs.push((chrom, at, end)),
        }
        at = end;
    }
    let names: Vec<u8> = seqs.iter().flat_map(|s| s.0.bytes().chain([0])).collect();
    let mut tbi: Vec<u8> = b"TBI\x01".to_vec();
    let i32s = |v: &mut Vec<u8>, xs: &[i32]| xs.iter().for_each(|x| v.extend(x.to_le_bytes()));
    // n_ref, VCF preset, columns 1 2 0, meta '#', skip 0, names.
    i32s(&mut tbi, &[seqs.len() as i32, 2, 1, 2, 0, b'#' as i32, 0]);
    i32s(&mut tbi, &[names.len() as i32]);
    tbi.extend(&names);
    for (_, beg, end) in &seqs {
        // One bin (0, the whole sequence), one chunk, no linear index.
        i32s(&mut tbi, &[1, 0, 1]);
        tbi.extend(beg.to_le_bytes());
        tbi.extend(end.to_le_bytes());
        i32s(&mut tbi, &[0]);
    }
    std::fs::write(d.join("s.vcf.gz.tbi"), tbi).unwrap();
    (plain, gz)
}

/// (status, alt copies, note) of the site `name`, read through the index
/// and by a full pass, which must agree.
fn site_both_ways(d: &Path, records: &[String], name: &str) -> (String, Option<i64>, String) {
    let (plain, gz) = both_ways(d, records);
    let (whole, how) =
        read_vcf(&plain.to_string_lossy(), read_markers(&markers()).unwrap()).unwrap();
    assert_eq!(how, Access::Stream);
    let (indexed, how) =
        read_vcf(&gz.to_string_lossy(), read_markers(&markers()).unwrap()).unwrap();
    assert_eq!(how, Access::Index);
    assert_eq!(indexed, whole, "{records:?}");
    let s = &whole[name];
    (s.status.as_str().to_owned(), s.n, s.note.clone())
}

#[test]
fn the_abo_insertion_counts_after_a_reference_or_other_record() {
    let d = dir("abo261");
    let ins = |gt: &str| line("chr9", 133257521, "T", "TC", "PASS", gt);
    for (first, gt, copies) in [
        // An all-sites VCF (bcftools call -m without -v): the reference
        // record comes first, then the insertion.
        (line("chr9", 133257521, "T", ".", ".", "0/0"), "0/1", 1),
        (line("chr9", 133257521, "T", ".", ".", "0/0"), "1/1", 2),
        // A gVCF-style reference record.
        (line("chr9", 133257521, "T", "<*>", ".", "0/0"), "0/1", 1),
    ] {
        for recs in [[first.clone(), ins(gt)], [ins(gt), first.clone()]] {
            let (status, n, note) = site_both_ways(&d, &recs, "ABO_261");
            assert_eq!(
                (status.as_str(), n, note.as_str()),
                ("called", Some(copies), ""),
                "{recs:?}"
            );
        }
    }
    // A split multi-allelic record (bcftools norm -m-any writes the SNV
    // first): the insertion counts and the SNV is named.
    let snv = line("chr9", 133257521, "T", "G", "PASS", "1/0");
    for recs in [[snv.clone(), ins("0/1")], [ins("0/1"), snv.clone()]] {
        let (status, n, note) = site_both_ways(&d, &recs, "ABO_261");
        assert_eq!((status.as_str(), n), ("called", Some(1)), "{recs:?}");
        assert_eq!(note, "another allele at this position (T>G, genotype 1/0)");
    }
    // The whole report: A1/O.01, not O.
    let t = run(
        &d,
        &[line("chr9", 133257521, "T", ".", ".", "0/0"), ins("0/1")],
        Some(depth(&d, 2.0, &HashMap::new())),
    );
    assert_eq!(t["ABO"][..3], ["A1", "A1/O.01", "low"]);
}

#[test]
fn a_later_record_of_the_marker_allele_that_is_no_call_wins_over_a_reference_record() {
    let d = dir("abo261-nocall");
    let r = line("chr9", 133257521, "T", ".", ".", "0/0");
    for (later, note) in [
        (
            line("chr9", 133257521, "T", "TC", "PASS", "./."),
            "genotype ./.",
        ),
        (
            line("chr9", 133257521, "T", "TC", "LowQual", "0/1"),
            "FILTER LowQual",
        ),
    ] {
        for recs in [[r.clone(), later.clone()], [later.clone(), r.clone()]] {
            let (status, n, got) = site_both_ways(&d, &recs, "ABO_261");
            assert_eq!(
                (status.as_str(), n, got.as_str()),
                ("no call", None, note),
                "{recs:?}"
            );
        }
    }
}

#[test]
fn a_shorter_reference_record_leaves_a_deletion_marker_to_the_reads() {
    let d = dir("abo1061");
    let ins = line("chr9", 133257521, "T", "TC", "PASS", "0/1");
    // An all-sites VCF writes `C . 0/0` at the A2 marker (REF CG): it says
    // nothing past its own base, so the reads decide, as with no record.
    for r in [
        line("chr9", 133255669, "C", ".", ".", "0/0"),
        line("chr9", 133255669, "C", "<*>", ".", "0/0"),
    ] {
        let recs = [r.clone(), ins.clone()];
        assert_eq!(
            site_both_ways(&d, &recs, "ABO_1061"),
            ("absent".to_owned(), None, String::new())
        );
        let t = run(&d, &recs, Some(depth(&d, 2.0, &HashMap::new())));
        assert_eq!(t["ABO"][..3], ["A1", "A1/O.01", "low"], "{r}");
        // The deletion after it still counts, in either order.
        let del = line("chr9", 133255669, "CG", "C", "PASS", "0/1");
        for recs in [[r.clone(), del.clone()], [del.clone(), r.clone()]] {
            let (status, n, _) = site_both_ways(&d, &recs, "ABO_1061");
            assert_eq!((status.as_str(), n), ("called", Some(1)), "{recs:?}");
        }
    }
    // A shorter record with a non-reference call stays a no call.
    let (status, _, note) = site_both_ways(
        &d,
        &[line("chr9", 133255669, "C", "T", "PASS", "0/1")],
        "ABO_1061",
    );
    assert_eq!(
        (status.as_str(), note.as_str()),
        ("no call", "REF C at this position, the marker's REF is CG")
    );
}

#[test]
fn a_record_that_spans_a_marker_with_another_base_there_is_no_call() {
    let d = dir("kel578");
    let reads = || Some(depth(&d, 2.0, &HashMap::new()));
    // KEL c.578C>T (K) is chr7:142957921 G>A. A padded record one base
    // upstream carries the same change (CG>CA), or a deletion of that
    // base with it (DeepVariant's CG>C,CA 1/2).
    for (r, alt, gt) in [
        ("CG", "CA", "1/1"),
        ("CG", "C,CA", "1/2"),
        ("CG", "C", "0/1"),
    ] {
        let recs = [line("chr7", 142957920, r, alt, "PASS", gt)];
        let (status, n, note) = site_both_ways(&d, &recs, "KEL_578");
        assert_eq!(
            (status.as_str(), n, note.as_str()),
            (
                "no call",
                None,
                format!("a record at chr7:142957920 ({r}>{alt}, genotype {gt}) spans this site")
                    .as_str()
            )
        );
        let t = run(&d, &recs, reads());
        assert!(t["Kell"][0].starts_with("K/k not determined"), "{recs:?}");
        assert_ne!(t["Kell"][2], "high", "{recs:?}");
    }
    // Split (bcftools norm -m-any): the deletion before the SNV record
    // still makes the SNV's copy count unknown (one copy may be deleted).
    let recs = [
        line("chr7", 142957920, "CG", "C", "PASS", "0/1"),
        line("chr7", 142957921, "G", "A", "PASS", "0/1"),
    ];
    let (status, ..) = site_both_ways(&d, &recs, "KEL_578");
    assert_eq!(status, "no call");
    // A spanning record that keeps the marker's base, or calls reference,
    // leaves the site to its own records and the reads.
    for recs in [
        [line("chr7", 142957920, "CG", "TG", "PASS", "1/1")],
        [line("chr7", 142957920, "CG", "C", "PASS", "0/0")],
        [line("chr7", 142957920, "CG", "C", "PASS", "./.")],
    ] {
        let (status, ..) = site_both_ways(&d, &recs, "KEL_578");
        assert_eq!(status, "absent", "{recs:?}");
    }
}

/// The complete table and evidence must agree through both VCF readers.
fn report_both_ways(d: &Path, records: &[String]) -> (String, String) {
    let (plain, gz) = both_ways(d, records);
    let mut outputs = Vec::new();
    for (vcf, access) in [(plain, Access::Stream), (gz, Access::Index)] {
        let files = Files {
            vcf: vcf.to_string_lossy().into_owned(),
            markers: markers(),
            depth: Some(depth(d, 2.0, &HashMap::new())),
            cn_loci: Some(cn_loci()),
            out: d.join("table.tsv").to_string_lossy().into_owned(),
            evidence: d.join("evidence.tsv").to_string_lossy().into_owned(),
        };
        assert_eq!(report(&files).unwrap().1, access);
        outputs.push((
            std::fs::read_to_string(&files.out).unwrap(),
            std::fs::read_to_string(&files.evidence).unwrap(),
        ));
    }
    assert_eq!(outputs[0], outputs[1]);
    outputs.remove(0)
}

fn split_null_marker(name: &str, system: &str, phenotype: &str) {
    let d = dir(&format!("split-{name}"));
    let m = read_markers(&markers()).unwrap()[name].clone();
    let alts: Vec<&str> = m.alt.split(',').collect();
    let record = |alt: &str, gt: &str| line(&m.chrom, m.pos, &m.ref_, alt, "PASS", gt);
    let unsplit = record(&m.alt, "1/2");
    let expected = report_both_ways(&d, std::slice::from_ref(&unsplit));
    let row: Vec<&str> = expected
        .0
        .lines()
        .find(|l| l.starts_with(&format!("{system}\t")))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(row[..4], [system, phenotype, "null/null", "high"]);
    let a = record(alts[0], "1/0");
    let b = record(alts[1], "0/1");
    for records in [
        vec![a.clone(), b.clone()],
        vec![b.clone(), a.clone()],
        vec![a.clone(), a.clone(), b.clone(), b.clone()],
        vec![unsplit.clone(), a.clone(), b.clone()],
        vec![b.clone(), a.clone(), unsplit],
        vec![record(&format!("{},{}", alts[1], alts[0]), "2/1")],
        vec![record(".", "0/0"), a, b],
    ] {
        assert_eq!(report_both_ways(&d, &records), expected, "{records:?}");
    }
    // A repeated call of one ALT is still one copy, even when another
    // record lists the same ALT at a different allele index.
    let one = record(alts[0], "0/1");
    for gt in ["0/1", "1/1"] {
        let single = record(alts[0], gt);
        let expected = report_both_ways(&d, std::slice::from_ref(&single));
        let copies = if gt == "0/1" { 1 } else { 2 };
        let reordered = record(
            &format!("{},{}", alts[1], alts[0]),
            if copies == 1 { "0/2" } else { "2/2" },
        );
        let records = [single.clone(), reordered, single];
        assert_eq!(report_both_ways(&d, &records), expected);
        assert_eq!(
            site_both_ways(&d, &records, name),
            ("called".to_owned(), Some(copies), String::new())
        );
    }
    assert_eq!(site_both_ways(&d, &[one], name).1, Some(1));
}

#[test]
fn split_kidd_null_alts_match_an_unsplit_call() {
    split_null_marker("JK_01N06", "Kidd", "Jk(a-b-)");
}

#[test]
fn split_kell_null_alts_match_an_unsplit_call() {
    split_null_marker("KEL_02N01", "Kell", "K-k- Kp(a-b-) Js(a-b-)");
}

#[test]
fn incompatible_null_alt_dosages_are_no_call_in_both_readers() {
    for name in ["JK_01N06", "KEL_02N01"] {
        let d = dir(&format!("incompatible-{name}"));
        let m = read_markers(&markers()).unwrap()[name].clone();
        let alts: Vec<&str> = m.alt.split(',').collect();
        let record = |alt: &str, gt: &str| line(&m.chrom, m.pos, &m.ref_, alt, "PASS", gt);
        // Two distinct ALTs cannot have three copies. Different dosages
        // for the same ALT cannot both be true either.
        let hom = record(alts[0], "1/1");
        for other in [record(alts[1], "0/1"), record(alts[0], "0/1")] {
            for records in [
                vec![hom.clone(), other.clone()],
                vec![other.clone(), hom.clone()],
            ] {
                let (status, n, note) = site_both_ways(&d, &records, name);
                assert_eq!((status.as_str(), n), ("no call", None), "{records:?}");
                assert_eq!(note, "incompatible ALT dosages at this position");
                let (table, _) = report_both_ways(&d, &records);
                let row: Vec<&str> = table
                    .lines()
                    .find(|l| l.starts_with(&format!("{}\t", m.system)))
                    .unwrap()
                    .split('\t')
                    .collect();
                assert_ne!(row[3], "high");
            }
        }
    }
}

/// A one-sample VCF record at the marker `name` with these FORMAT keys and values.
fn recf(name: &str, fmt: &str, smp: &str) -> String {
    let m = &read_markers(&markers()).unwrap()[name];
    format!(
        "{}\t{}\t.\t{}\t{}\t50\tPASS\t.\t{fmt}\t{smp}\n",
        m.chrom, m.pos, m.ref_, m.alt
    )
}

#[test]
fn a_call_under_the_gq_and_dp_floor_lowers_confidence() {
    let d = dir("weak");
    let reads = || Some(depth(&d, 2.0, &HashMap::new()));
    // KEL c.578C>T 0/1 with GQ 0 and one read: K+k+ still, at confidence low.
    let t = run(&d, &[recf("KEL_578", "GT:GQ:DP", "0/1:0:1")], reads());
    assert_eq!(t["Kell"][0], "K+k+ Kp(a-b+) Js(a-b+)");
    assert_eq!(t["Kell"][2], "low");
    assert!(
        t["Kell"][3].contains("KEL c.578C>T: weak call (GQ 0 under 20, DP 1 under 8)"),
        "{}",
        t["Kell"][3]
    );
    let ev = std::fs::read_to_string(d.join("e.tsv")).unwrap();
    let row = ev.lines().find(|l| l.contains("\tKEL_578\t")).unwrap();
    assert!(row.ends_with("GQ 0 under 20, DP 1 under 8"), "{row}");
    // The same call with GQ 20 and 8 reads, at the floors, stays high.
    let t = run(&d, &[recf("KEL_578", "GT:GQ:DP", "0/1:20:8")], reads());
    assert_eq!(
        t["Kell"][..3],
        ["K+k+ Kp(a-b+) Js(a-b+)", "K/k Kpb/Kpb Jsb/Jsb", "high"]
    );
    // A gVCF reference block with GQ 3 at a site: low too.
    let t = run(&d, &[recf("JK_838", "GT:GQ:MIN_DP", "0/0:3:20")], reads());
    assert_eq!(t["Kidd"][..3], ["Jk(a+b-)", "Jka/Jka", "low"]);
}

#[test]
fn a_vcf_with_no_chr9_record_reads_abo_not_assessed_without_a_bam() {
    let d = dir("nochr9");
    // A VCF with records elsewhere, none on chr9, and no BAM: ABO is not
    // assessed, not O; Kell, with its record, is read.
    let t = run(&d, &[rec("KEL_578", "0/1")], None);
    assert_eq!(t["ABO"][..3], ["not assessed", ".", "not assessed"]);
    assert!(t["ABO"][3].contains("no VCF record at any of its 5 sites and no BAM"));
    assert!(t["ABO"][3].contains("chr9"));
    assert_eq!(t["Kell"][0], "K+k+ Kp(a-b+) Js(a-b+)");
    // One ABO record is enough to read the others as the reference.
    let t = run(&d, &[rec("ABO_802", "0/0")], None);
    assert_eq!(t["ABO"][..3], ["O", "O.01/O.01", "medium"]);
    // A VCF with no record at all is no error: every system not assessed,
    // and RhD not determined (it needs the BAM).
    let t = run(&d, &[], None);
    assert_eq!(t["ABO"][0], "not assessed");
    assert_eq!(t["Vel"][0], "not assessed");
    assert_eq!(t["RhD"][0], "not determined");
    // With reads the reads decide, as before.
    let t = run(&d, &[], Some(depth(&d, 2.0, &HashMap::new())));
    assert_eq!(t["ABO"][..3], ["O", "O.01/O.01", "high"]);
}

#[test]
fn a_phased_gata_box_goes_on_the_chromosome_the_vcf_gives() {
    let d = dir("duffy-phase");
    let reads = || Some(depth(&d, 2.0, &HashMap::new()));
    let duffy = |fy: &str, gata: &str| {
        let t = run(
            &d,
            &[recf("FY_GATA", "GT:PS", gata), recf("FY_125", "GT:PS", fy)],
            reads(),
        );
        t["Duffy"].clone()
    };
    // On the Fya chromosome: Fya silenced on red cells, Fyb expressed.
    let r = duffy("0|1:7", "1|0:7");
    assert_eq!(r[..3], ["Fy(a-b+)", "FY*A GATA/FY*B", "high"]);
    assert!(
        r[3].contains("the VCF's phase puts it on the Fya chromosome"),
        "{}",
        r[3]
    );
    // On the Fyb chromosome.
    let r = duffy("0|1:7", "0|1:7");
    assert_eq!(r[..2], ["Fy(a+b-)", "FY*A/FY*B GATA"]);
    assert!(r[3].contains("the VCF's phase puts it on the Fyb chromosome"));
    // Two phase sets, or no phase: Fyb first, as before.
    for (fy, gata) in [("0|1:7", "1|0:9"), ("0/1:.", "0/1:.")] {
        let r = duffy(fy, gata);
        assert_eq!(r[..2], ["Fy(a+b-)", "FY*A/FY*B GATA"], "{fy} {gata}");
        assert!(r[3].contains("with no phase it is placed on the Fyb chromosome first"));
    }
}

#[test]
fn an_abo_a_allele_reads_low_and_o_or_b_keeps_its_confidence() {
    let d = dir("abo-a-low");
    let reads = || Some(depth(&d, 2.0, &HashMap::new()));
    let t = run(&d, &[rec("ABO_261", "0/1")], reads());
    assert_eq!(t["ABO"][..3], ["A1", "A1/O.01", "low"]);
    assert!(t["ABO"][3].contains("an A allele is read from a c.261 insertion"));
    let t = run(
        &d,
        &[
            rec("ABO_261", "0/1"),
            rec("ABO_803", "0/1"),
            rec("ABO_796", "0/1"),
        ],
        reads(),
    );
    assert_eq!(t["ABO"][..3], ["B", "B/O.01", "high"]);
    assert!(t["ABO"][3].contains("ABO alleles outside these five variants"));
    assert!(!t["ABO"][3].contains("an A allele is read"));
}

#[test]
fn rows_state_what_dna_cannot_check() {
    let d = dir("limits");
    let t = run(&d, &[], Some(depth(&d, 2.0, &HashMap::new())));
    assert!(t["RhD"][3].contains("RHD-CE-D hybrids that keep exons 3 to 7"));
    for system in [
        "Vel",
        "Langereis",
        "Junior",
        "H",
        "Secretor",
        "Kell",
        "Kidd",
    ] {
        assert!(
            t[system][3].contains("a null allele the table does not hold"),
            "{system}: {}",
            t[system][3]
        );
    }
    assert!(!t["Lutheran"][3].contains("null allele"));
    // Both deleted: D-, and no D+ limit to state.
    let t = run(&d, &[], Some(depth(&d, 0.0, &HashMap::new())));
    assert!(!t["RhD"][3].contains("RHD-CE-D hybrids"));
}
