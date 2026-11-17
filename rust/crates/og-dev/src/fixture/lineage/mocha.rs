//! Two synthetic phased VCFs for step 51 (MoChA): one with a CN-LOH on 20p.

use super::{Fail, Fasta, PyRandom, bad, create, want};
use std::io::Write;

const CHROMS: [&str; 4] = ["chr16", "chr19", "chr20", "chr22"];
/// 20p from the telomere.
const EVENT: (&str, u64, u64) = ("chr20", 1, 20_000_000);
/// CN-LOH in 15% of the cells: haplotype balance 0.575.
const BDEV: f64 = 0.075;

/// Python's Knuth sampler of the case: multiply uniforms until the product
/// drops to `exp(-lam)`.
fn poisson(rnd: &mut PyRandom, lam: f64) -> u64 {
    let l = (-lam).exp();
    let (mut k, mut p) = (0, 1.0);
    loop {
        p *= rnd.random();
        if p <= l {
            return k;
        }
        k += 1;
    }
}

/// `sum(random() < f for _ in range(n))`.
fn binom(rnd: &mut PyRandom, n: u64, f: f64) -> u64 {
    (0..n).filter(|_| rnd.random() < f).count() as u64
}

/// `og dev fixture mocha-vcf FASTA OUT_SHIFT OUT_FLAT`: phased SNVs about
/// every kilobase of chr16, chr19, chr20 and chr22 with allelic depths
/// (seed 51 for each file): OUT_SHIFT (sample MOSAIC) carries a CN-LOH on
/// chr20 up to 20 Mb in 15% of the cells, OUT_FLAT (MOSAIC0) none. Prints
/// each file with its line count.
pub fn vcf(args: &[String], out: &mut dyn Write, _err: &mut dyn Write) -> Result<(), Fail> {
    want(args, 3)?;
    let mut fa = Fasta::open(&args[0])?;
    for (path, shift) in [(&args[1], true), (&args[2], false)] {
        let mut rnd = PyRandom::new(51);
        let mut lines = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
        ];
        for c in CHROMS {
            let len = fa
                .len(c)
                .ok_or_else(|| bad(format!("{}.fai: no {c}", args[0])))?;
            lines.push(format!("##contig=<ID={c},length={len}>"));
        }
        lines.extend(
            [
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
            ]
            .map(String::from),
        );
        let name = if shift { "MOSAIC" } else { "MOSAIC0" };
        lines.push(format!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}"
        ));
        for c in CHROMS {
            let length = fa.len(c).unwrap_or(0);
            let mut pos: u64 = 100_000;
            while pos + 100_000 < length {
                pos += rnd.randint(800, 1200) as u64;
                let base = fa.base(c, pos)?;
                let p = rnd.uniform(0.05, 0.5);
                let r = rnd.random();
                let depth = poisson(&mut rnd, 30.0);
                let o = rnd.random() < 0.5;
                let bal = if shift && c == EVENT.0 && EVENT.1 <= pos && pos <= EVENT.2 {
                    0.5 + BDEV
                } else {
                    0.5
                };
                let k_het = binom(&mut rnd, depth, bal);
                let k_err = binom(&mut rnd, depth, 0.005);
                let (r_base, alt_base) = match base {
                    Some(b'A') => ('A', 'G'),
                    Some(b'G') => ('G', 'A'),
                    Some(b'C') => ('C', 'T'),
                    Some(b'T') => ('T', 'C'),
                    _ => continue,
                };
                let het = 2.0 * p * (1.0 - p);
                let (gt, alt) = if r < het {
                    if o {
                        ("0|1", depth - k_het)
                    } else {
                        ("1|0", k_het)
                    }
                } else if r < het + p * p {
                    ("1|1", depth - k_err)
                } else {
                    continue;
                };
                lines.push(format!(
                    "{c}\t{pos}\t.\t{r_base}\t{alt_base}\t50\tPASS\t.\tGT:GQ:DP:AD\t{gt}:50:{depth}:{},{alt}",
                    depth - alt
                ));
            }
        }
        let mut o = create(path)?;
        o.write_all(lines.join("\n").as_bytes())?;
        o.write_all(b"\n")?;
        o.flush()?;
        writeln!(out, "{path} {} lines", lines.len())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::testutil::{dir, fasta, fnv, run};

    #[test]
    fn vcfs_are_the_pythons_bytes() {
        let d = dir("mocha");
        fasta(
            &d.join("ref.fa"),
            &[
                ("chr16", 260_000),
                ("chr19", 230_000),
                ("chr20", 300_000),
                ("chr22", 201_000),
            ],
            21,
        );
        let p = |n: &str| d.join(n).to_str().expect("utf8").to_string();
        let (code, out, err) = run("mocha-vcf", &[&p("ref.fa"), &p("s.vcf"), &p("f.vcf")]);
        assert_eq!(code, 0, "{err}");
        assert!(out.starts_with(&format!("{} ", p("s.vcf"))), "{out}");
        // FNV-1a of what the Python heredoc of mocha-1.sh wrote.
        for (f, want) in [
            ("s.vcf", 0xdd2b_0f2a_63b7_1109),
            ("f.vcf", 0x6b70_c06a_6cd2_a61d),
        ] {
            let got = fnv(&std::fs::read(d.join(f)).expect("read"));
            assert_eq!(got, want, "{f}: {got:#x}");
        }
    }
}
