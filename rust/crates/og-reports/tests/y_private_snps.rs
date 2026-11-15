//! `og y private-snps` on a synthetic VCF, synthetic position tables and a
//! synthetic callable mask (no real sample). It replaces
//! `tests/test_y_private_snps.py`.
//!
//! 1. a chrY call at a position one of the public SNP tables names is never
//!    private; both tables are checked, not only the first;
//! 2. a call at a position no table names is listed, with its depth, ALT
//!    fraction, genotype quality and genotype;
//! 3. the gates: a position in a pseudoautosomal region, a record that did
//!    not pass the filter, a genotype carrying the reference allele, an
//!    indel, a second ALT allele, depth below `--min-depth`, genotype
//!    quality below `--min-gq` and an ALT fraction below `--min-fraction`
//!    are all dropped, each counted in the file's header;
//! 4. depth far above the sample's own chrY median is dropped, and
//!    `--max-depth-factor 0` turns that gate off;
//! 5. a gzipped VCF reads the same as a plain one;
//! 6. a missing position table stops the tool and names `setup.sh
//!    --yleaf-data`;
//! 7. a call outside the callable mask (a base not marked `P`) is dropped
//!    and counted as outside_callable, whatever its quality;
//!    `--callable-mask` points at another mask; a missing mask stops the
//!    tool.

use std::fs;
use std::io::Write;
use std::num::NonZero;
use std::path::{Path, PathBuf};

const HEADER: &str = "##fileformat=VCFv4.2\n##contig=<ID=chrY,length=57227415>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n";

/// pos, ref, alt, filter, GT:AD:DP:GQ
const RECORDS: [(u32, &str, &str, &str, &str); 15] = [
    (7000100, "G", "A", "PASS", "1:0,30:30:50"), // on the first table
    (7000200, "G", "A", "PASS", "1:0,30:30:50"), // on the second table
    (7000300, "C", "T", "PASS", "1:0,30:30:50"), // private
    (7000400, "C", "T", "PASS", "1/1:1,29:30:50"), // private, diploid call
    (2000000, "C", "T", "PASS", "1:0,30:30:50"), // PAR1
    (56900000, "C", "T", "PASS", "1:0,30:30:50"), // PAR2
    (7000500, "C", "T", "RefCall", "1:0,30:30:50"), // not PASS
    (7000600, "C", "T", "PASS", "0/1:15,15:30:50"), // carries the reference allele
    (7000700, "C", "CTT", "PASS", "1:0,30:30:50"), // insertion
    (7000800, "C", "T,G", "PASS", "1:0,20,10:30:50"), // second ALT allele
    (7000900, "C", "T", "PASS", "1:0,9:9:50"),   // depth below the gate
    (7001000, "C", "T", "PASS", "1:0,30:30:5"),  // genotype quality below the gate
    (7001100, "C", "T", "PASS", "1:6,24:30:50"), // ALT fraction 0.8
    (7001200, "C", "T", "PASS", "1:0,300:300:50"), // ten times the median depth
    (7001300, "C", "T", "PASS", "1:0,30:30:50"), // outside the callable mask
];

/// Every base passes (`P`) up to MASK_LEN, except a run of `Z` (too many
/// zero mapping-quality reads) around 7001300.
const MASK_LEN: usize = 7_002_000;
const MASK_Z: (usize, usize) = (7_001_250, 7_001_350);

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("og-y-private-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// `text` gzipped (as BGZF, which gzip readers read).
fn gz(path: &Path, text: &[u8]) {
    let mut w =
        og_io::bgzf::Writer::new(fs::File::create(path).unwrap(), NonZero::<usize>::MIN).unwrap();
    w.write_all(text).unwrap();
    w.finish().unwrap();
}

fn write_mask(path: &Path, z: (usize, usize)) {
    let mut seq = vec![b'P'; MASK_LEN];
    seq[z.0 - 1..z.1].fill(b'Z');
    let mut text = b">chrY\n".to_vec();
    for line in seq.chunks(60) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }
    gz(path, &text);
}

/// The columns of the two Yleaf GRCh38 tables: chromosome, name,
/// haplogroup, position, mutation, ancestral, derived.
fn write_tables(d: &Path) {
    fs::write(
        d.join("positions_ftdna_hg38.txt"),
        "chrY\tpub1\tR\t7000100\tG->A\tG\tA\n",
    )
    .unwrap();
    fs::write(
        d.join("positions_yfull_hg38.txt"),
        "chrY\tpub2\tR\t7000200\tG->A\tG\tA\n",
    )
    .unwrap();
}

fn write_vcf(path: &Path) {
    let mut text = HEADER.to_owned();
    for (pos, r, a, f, s) in RECORDS {
        text += &format!("chrY\t{pos}\t.\t{r}\t{a}\t50\t{f}\t.\tGT:AD:DP:GQ\t{s}\n");
    }
    if path.extension().is_some_and(|e| e == "gz") {
        gz(path, text.as_bytes());
    } else {
        fs::write(path, text).unwrap();
    }
}

struct Run {
    code: u8,
    stdout: String,
    stderr: String,
    head: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Run {
    fn pos(&self) -> Vec<&str> {
        self.rows.iter().map(|r| r[0].as_str()).collect()
    }
    fn head_has(&self, label: &str) -> bool {
        self.head.iter().any(|h| h.contains(label))
    }
}

fn run(d: &Path, vcf: &Path, prefix: &str, extra: &[&str]) -> Run {
    let mut args: Vec<String> = [
        "--vcf",
        vcf.to_str().unwrap(),
        "--tree-dir",
        d.to_str().unwrap(),
        "--sample",
        "S",
        "--out-prefix",
        d.join(prefix).to_str().unwrap(),
    ]
    .map(str::to_owned)
    .to_vec();
    args.extend(extra.iter().map(|s| (*s).to_owned()));
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = og_reports::y::cli(Some("private-snps"), &args, &mut out, &mut err);
    let text =
        fs::read_to_string(d.join(format!("{prefix}_y_private_snps.tsv"))).unwrap_or_default();
    let mut head = Vec::new();
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            head.push(line.to_owned());
        } else if !line.starts_with("Pos\t") {
            rows.push(line.split('\t').map(str::to_owned).collect());
        }
    }
    Run {
        code,
        stdout: String::from_utf8(out).unwrap(),
        stderr: String::from_utf8(err).unwrap(),
        head,
        rows,
    }
}

#[test]
fn private_snps_gates_tables_and_mask() {
    let d = scratch("main");
    write_tables(&d);
    write_mask(&d.join("callable_mask_chrY.fasta.gz"), MASK_Z);
    let vcf = d.join("sample.vcf");
    write_vcf(&vcf);
    let r = run(&d, &vcf, "S", &[]);
    assert_eq!(r.code, 0, "the tool exits 0: {}", r.stderr);
    let pos = r.pos();
    assert!(
        !pos.contains(&"7000100"),
        "1. pub1 (first table) is not private: {pos:?}"
    );
    assert!(
        !pos.contains(&"7000200"),
        "1. pub2 (second table) is not private: {pos:?}"
    );
    assert!(
        r.head_has("on_public_lists=1"),
        "1. the header counts the two public positions: {:?}",
        r.head
    );
    assert!(
        pos.contains(&"7000300"),
        "2. the private haploid call is listed"
    );
    assert!(
        pos.contains(&"7000400"),
        "2. the private diploid ALT/ALT call is listed"
    );
    assert!(
        r.rows.contains(
            &["7000300", "C", "T", "30", "1.000", "50", "1"]
                .map(str::to_owned)
                .to_vec()
        ),
        "2. its row carries the depth, ALT fraction, quality and genotype: {:?}",
        r.rows
    );
    let mut sorted = pos.clone();
    sorted.sort_unstable();
    assert_eq!(
        sorted,
        ["7000300", "7000400"],
        "3. no other position is listed"
    );
    for label in [
        "PAR=2",
        "not_PASS=1",
        "genotype=1",
        "depth=1",
        "quality=1",
        "alt_fraction=1",
        "excess_depth=1",
        "outside_callable=1",
    ] {
        assert!(
            r.head_has(label),
            "3. the header counts {label}: {:?}",
            r.head
        );
    }
    assert!(
        r.stdout
            .lines()
            .any(|l| l.contains("13 chrY SNV records read")),
        "3. the indel and the second-ALT record never reach the gates: {}",
        r.stdout
    );
    assert!(
        r.stdout.contains("unconfirmed until another man"),
        "the step line says a private SNP is unconfirmed: {}",
        r.stdout
    );

    // 4. the depth gate off
    let off = run(&d, &vcf, "off", &["--max-depth-factor", "0"]);
    assert!(
        off.pos().contains(&"7001200"),
        "4. --max-depth-factor 0 keeps the deep call"
    );

    // 5. gzipped
    let gzv = d.join("s.vcf.gz");
    write_vcf(&gzv);
    assert_eq!(
        run(&d, &gzv, "gz", &[]).pos(),
        ["7000300", "7000400"],
        "5. a gzipped VCF gives the same rows"
    );

    // 7. the callable mask
    assert!(
        !pos.contains(&"7001300"),
        "7. the call outside the mask is not listed"
    );
    assert!(
        r.head_has("callable_mask=callable_mask_chrY.fasta.gz"),
        "7. the header names the mask it used: {:?}",
        r.head
    );
    let other = d.join("other_mask.fa.gz");
    write_mask(&other, (7_000_250, 7_000_350));
    assert_eq!(
        run(
            &d,
            &vcf,
            "alt",
            &["--callable-mask", other.to_str().unwrap()]
        )
        .pos(),
        ["7000400", "7001300"],
        "7. --callable-mask reads another mask (7000300 masked there, 7001300 callable)"
    );
    fs::remove_file(d.join("callable_mask_chrY.fasta.gz")).unwrap();
    let r = run(&d, &vcf, "nomask", &[]);
    assert!(
        r.code != 0 && r.stderr.contains("setup.sh --yleaf-data"),
        "7. a missing mask stops the tool and names setup.sh --yleaf-data: {}",
        r.stderr
    );
    write_mask(&d.join("callable_mask_chrY.fasta.gz"), MASK_Z);

    // 6. a missing table
    fs::remove_file(d.join("positions_yfull_hg38.txt")).unwrap();
    let r = run(&d, &vcf, "missing", &[]);
    assert!(
        r.code != 0 && r.stderr.contains("setup.sh --yleaf-data"),
        "6. a missing table stops the tool and names setup.sh --yleaf-data: {}",
        r.stderr
    );
    fs::remove_dir_all(d).unwrap();
}
