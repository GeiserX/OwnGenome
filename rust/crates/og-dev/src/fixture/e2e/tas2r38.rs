//! `tas2r38 REF_DIR GENOME_DIR`: the made-up inputs of the step 42 rule
//! case on TAS2R38 (`tests/e2e/genotype-rules-1-phase.sh`).
//!
//! - `REF_DIR/chr7.fa`: chr7 of GRCh38's length, N everywhere except a 2.6
//!   kb window of random bases (`random.seed(42)`) whose three TAS2R38 sites
//!   carry the GRCh38 bases: rs10246939 T, rs1726866 G, rs713598 C.
//! - GRSPAN, GRSWAP and GRNOSPAN: a gVCF heterozygous at the three sites
//!   between reference blocks of GQ 42, the variant-only VCF, and reads:
//!   pairs spanning the sites from the PAV and AVI copies (GRSPAN), from the
//!   G,A,C and C,G,T copies (GRSWAP), or single 60-base reads that each cover
//!   one site (GRNOSPAN).
//! - GRREF: a gVCF that is one reference block, and no reads.
//!
//! Each goes to `GENOME_DIR/<sample>/vcf/<sample>.g.vcf`, `.vcf` and
//! `GENOME_DIR/<sample>/aligned/<sample>.sam`.

use super::str_extra::Lines;
use super::{Fail, PyRandom, create, exactly};
use std::fmt::Write as _;
use std::io::Write;
use std::path::Path;

/// GRCh38 chr7's length, and the window [W0, W1), 1-based.
const LEN: usize = 159_345_973;
const W0: usize = 141_972_000;
const W1: usize = 141_974_600;
/// rs10246939, rs1726866 and rs713598, in the order of their positions.
const S3: usize = 141_972_804;
const S2: usize = 141_972_905;
const S1: usize = 141_973_545;
const SITES: [(usize, u8, u8); 3] = [(S3, b'T', b'C'), (S2, b'G', b'A'), (S1, b'C', b'G')];

/// The window's bases, with the GRCh38 bases at the three sites.
fn window() -> Vec<u8> {
    let mut win = PyRandom::new(42).bases(W1 - W0);
    for (p, b, _) in SITES {
        win[p - W0] = b;
    }
    win
}

fn base(win: &[u8], p: usize) -> char {
    char::from(win[p - W0])
}

/// The window with `alleles` at S1, S2 and S3 (in that order).
fn hap(win: &[u8], alleles: &[u8; 3]) -> Vec<u8> {
    let mut h = win.to_vec();
    for (p, &b) in [S1, S2, S3].iter().zip(alleles) {
        h[p - W0] = b;
    }
    h
}

fn gvcf(win: &[u8], het: bool) -> String {
    let mut s = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chr7,length={LEN}>\n\
         ##ALT=<ID=*,Description=\"Any other allele\">\n\
         ##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of the block\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
         ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum depth\">\n\
         ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred likelihoods\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    );
    let block = |s: &mut String, start: usize, end: usize| {
        let _ = writeln!(
            s,
            "chr7\t{start}\t.\t{}\t<*>\t0\t.\tEND={end}\tGT:GQ:MIN_DP:PL\t0/0:42:30:0,90,900",
            base(win, start)
        );
    };
    let mut start = W0;
    if het {
        for (p, _, alt) in SITES {
            block(&mut s, start, p - 1);
            let _ = writeln!(
                s,
                "chr7\t{p}\t.\t{}\t{},<*>\t40\tPASS\t.\tGT:GQ:DP:PL\t0/1:35:40:400,0,400,990,990,990",
                base(win, p),
                char::from(alt)
            );
            start = p + 1;
        }
    }
    block(&mut s, start, W1 - 1);
    s
}

fn vcf(win: &[u8]) -> String {
    let mut s = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chr7,length={LEN}>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    );
    for (p, _, alt) in SITES {
        let _ = writeln!(
            s,
            "chr7\t{p}\t.\t{}\t{}\t40\tPASS\t.\tGT\t0/1",
            base(win, p),
            char::from(alt)
        );
    }
    s
}

fn bases(h: &[u8], at: usize, n: usize) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(&h[at - W0..at - W0 + n])
}

/// 20 read pairs: read 1 covers S3 and S2, its mate S1.
fn pairs(s: &mut String, h: &[u8], tag: &str) {
    let q = "I".repeat(150);
    for i in 0..20 {
        let (a, b) = (S3 - 40 + i, S1 - 60 + i);
        let tlen = b + 150 - a;
        let _ = writeln!(
            s,
            "{tag}{i}\t99\tchr7\t{a}\t60\t150M\t=\t{b}\t{tlen}\t{}\t{q}",
            bases(h, a, 150)
        );
        let _ = writeln!(
            s,
            "{tag}{i}\t147\tchr7\t{b}\t60\t150M\t=\t{a}\t-{tlen}\t{}\t{q}",
            bases(h, b, 150)
        );
    }
}

/// 20 single 60-base reads around each site: none covers two.
fn singles(s: &mut String, h: &[u8], tag: &str) {
    let q = "I".repeat(60);
    for (p, _, _) in SITES {
        for i in 0..20 {
            let a = p - 35 + (i % 10);
            let flag = if i % 2 == 1 { 16 } else { 0 };
            let _ = writeln!(
                s,
                "{tag}{p}_{i}\t{flag}\tchr7\t{a}\t60\t60M\t*\t0\t0\t{}\t{q}",
                bases(h, a, 60)
            );
        }
    }
}

fn sam(win: &[u8], sample: &str) -> String {
    let mut s = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr7\tLN:{LEN}\n");
    let (pav, avi) = (hap(win, b"GGC"), hap(win, b"CAT"));
    match sample {
        "GRSPAN" => {
            pairs(&mut s, &pav, "pav");
            pairs(&mut s, &avi, "avi");
        }
        "GRSWAP" => {
            pairs(&mut s, &hap(win, b"GAC"), "gac");
            pairs(&mut s, &hap(win, b"CGT"), "cgt");
        }
        _ => {
            singles(&mut s, &pav, "pav");
            singles(&mut s, &avi, "avi");
        }
    }
    s
}

fn write(path: &Path, text: &str) -> Result<(), Fail> {
    if let Some(d) = path.parent() {
        std::fs::create_dir_all(d).map_err(|e| Fail::Run(format!("{}: {e}", d.display())))?;
    }
    let mut out = create(path)?;
    out.write_all(text.as_bytes())?;
    out.flush()?;
    Ok(())
}

fn reference(w: &mut dyn Write, win: &[u8]) -> std::io::Result<()> {
    w.write_all(b">chr7\n")?;
    let mut lines = Lines::new(w);
    lines.ns(W0 - 1)?;
    lines.bases(win)?;
    lines.ns(LEN - (W1 - 1))?;
    lines.end()
}

pub fn run(args: &[String]) -> Result<(), Fail> {
    let a = exactly(args, 2, "REF_DIR GENOME_DIR")?;
    let (rd, g) = (Path::new(&a[0]), Path::new(&a[1]));
    let win = window();
    std::fs::create_dir_all(rd).map_err(|e| Fail::Run(format!("{}: {e}", rd.display())))?;
    let mut out = create(&rd.join("chr7.fa"))?;
    reference(&mut out, &win)?;
    out.flush()?;
    for s in ["GRSPAN", "GRSWAP", "GRNOSPAN"] {
        write(
            &g.join(s).join("vcf").join(format!("{s}.g.vcf")),
            &gvcf(&win, true),
        )?;
        write(&g.join(s).join("vcf").join(format!("{s}.vcf")), &vcf(&win))?;
        write(
            &g.join(s).join("aligned").join(format!("{s}.sam")),
            &sam(&win, s),
        )?;
    }
    write(
        &g.join("GRREF").join("vcf").join("GRREF.g.vcf"),
        &gvcf(&win, false),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixture::e2e::tests::{Hashed, fnv};

    // Pins of the bytes the inline Python of
    // tests/e2e/genotype-rules-1-phase.sh wrote.
    #[test]
    fn writes_the_pythons_bytes() {
        let win = window();
        let got: Vec<u64> = [
            gvcf(&win, true),
            gvcf(&win, false),
            vcf(&win),
            sam(&win, "GRSPAN"),
            sam(&win, "GRSWAP"),
            sam(&win, "GRNOSPAN"),
        ]
        .iter()
        .map(|t| fnv(t.as_bytes()))
        .collect();
        assert_eq!(got, PINS);
        let mut h = Hashed::default();
        reference(&mut h, &win).unwrap();
        assert_eq!(h.lines, 1 + LEN.div_ceil(60));
        assert_eq!(h.hash, PIN_REF);
    }

    #[test]
    fn the_copies_carry_their_alleles_at_the_three_sites() {
        let win = window();
        let pav = hap(&win, b"GGC");
        assert_eq!([pav[S3 - W0], pav[S2 - W0], pav[S1 - W0]], *b"CGG");
        assert_eq!([win[S3 - W0], win[S2 - W0], win[S1 - W0]], *b"TGC");
        let spans = sam(&win, "GRSPAN");
        assert_eq!(spans.lines().filter(|l| !l.starts_with('@')).count(), 80);
        // Read 1 of the first pair starts 40 bases before S3 and covers S2.
        let r1: Vec<&str> = spans.lines().nth(2).unwrap().split('\t').collect();
        assert_eq!(r1[3], (S3 - 40).to_string());
        assert_eq!(r1[9].as_bytes()[40], b'C');
        assert_eq!(r1[9].as_bytes()[40 + S2 - S3], b'G');
    }

    #[test]
    fn writes_every_file() {
        let d = crate::fixture::e2e::tests::scratch("tas2r38");
        let (code, _, err) = crate::fixture::e2e::tests::run(
            "tas2r38",
            &[
                d.join("ref").to_str().unwrap(),
                d.join("g").to_str().unwrap(),
            ],
        );
        assert_eq!(code, Some(0), "{err}");
        for f in [
            "GRSPAN/vcf/GRSPAN.g.vcf",
            "GRSWAP/aligned/GRSWAP.sam",
            "GRNOSPAN/vcf/GRNOSPAN.vcf",
            "GRREF/vcf/GRREF.g.vcf",
        ] {
            assert!(d.join("g").join(f).is_file(), "{f}");
        }
        assert!(!d.join("g/GRREF/aligned/GRREF.sam").exists());
        assert_eq!(
            std::fs::metadata(d.join("ref/chr7.fa")).unwrap().len() as usize,
            6 + LEN + LEN.div_ceil(60)
        );
        std::fs::remove_dir_all(d).unwrap();
    }

    const PINS: [u64; 6] = [
        15_484_360_573_342_373_983,
        9_756_800_734_539_015_102,
        16_850_025_374_021_543_854,
        12_356_292_411_453_828_175,
        143_055_192_633_579_987,
        6_433_961_028_917_585_143,
    ];
    const PIN_REF: u64 = 12_401_145_397_378_992_347;
}
