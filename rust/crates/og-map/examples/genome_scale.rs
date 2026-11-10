//! A seeded, invented genome-scale run for the map's budgets:
//!
//! ```text
//! cargo run --release -p og-map --example genome_scale -- OUT_DIR [SEED]
//! ```
//!
//! writes `OUT_DIR/DEMO-009/vcf/DEMO-009_roh.txt`: a `bcftools roh` text of
//! about 200 MB, with an `ST` line every 520 bases on average over chr1-22,
//! X and Y and about 4,000 `RG` rows (a 30x genome's count), a few of them
//! 5 Mb or longer. The same seed writes the same bytes.

use og_map::reference::Reference;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const SAMPLE: &str = "DEMO-009";
const RUNS: u64 = 4_300;

/// splitmix64: small, seeded, the same numbers on every platform.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Uniform in `[lo, hi)`.
    fn range(&mut self, lo: u64, hi: u64) -> u64 {
        lo + self.next() % (hi - lo)
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn main() -> std::io::Result<()> {
    let mut args = std::env::args().skip(1);
    let Some(out) = args.next().map(PathBuf::from) else {
        eprintln!("usage: genome_scale OUT_DIR [SEED]");
        std::process::exit(2);
    };
    let seed = args
        .next()
        .map_or(1, |s| s.parse().expect("SEED is a number"));
    let assets = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../assets");
    let reference = Reference::load(&assets).expect("the band asset reads");
    let total: u64 = reference.chroms.iter().map(|c| c.len).sum();

    let dir = out.join(SAMPLE).join("vcf");
    fs::create_dir_all(&dir)?;
    let mut w = BufWriter::with_capacity(
        1 << 20,
        File::create(dir.join(format!("{SAMPLE}_roh.txt")))?,
    );
    writeln!(
        w,
        "# This file was produced by: bcftools roh (invented: og-map genome_scale, seed {seed})"
    )?;
    writeln!(
        w,
        "# ST\t[2]Sample\t[3]Chromosome\t[4]Position\t[5]State (0:HW, 1:AZ)\t[6]Quality"
    )?;
    writeln!(
        w,
        "# RG\t[2]Sample\t[3]Chromosome\t[4]Start\t[5]End\t[6]Length (bp)\t[7]Number of markers\t[8]Quality"
    )?;

    let mut rng = Rng(seed);
    let (mut sites, mut runs) = (0u64, 0u64);
    for c in &reference.chroms {
        let n = (RUNS * c.len).div_ceil(total);
        // Runs one after another along the chromosome, never overlapping:
        // the gaps are sized so that about `n` runs fit.
        let mean_gap = (c.len / n).saturating_sub(560_000).max(20_000);
        let mut rg = Vec::with_capacity(n as usize);
        let mut cursor = 0;
        for _ in 0..n {
            let len = if rng.range(0, 100) == 0 {
                rng.range(5_000_000, 12_000_000)
            } else {
                // Log-uniform between 100 kb and 1.5 Mb.
                (100_000.0 * 15f64.powf(rng.unit())) as u64
            };
            let start = cursor + 1 + rng.range(0, 2 * mean_gap);
            let end = start + len - 1;
            if end >= c.len {
                break;
            }
            cursor = end;
            let markers = len / 1_400 + 1;
            let q = 30.0 + 60.0 * rng.unit();
            rg.push(format!(
                "RG\t{SAMPLE}\t{}\t{start}\t{end}\t{len}\t{markers}\t{q:.1}\n",
                c.name
            ));
        }
        let mut pos = 10_000;
        while pos < c.len {
            let state = u8::from(rng.range(0, 10) < 3);
            writeln!(w, "ST\t{SAMPLE}\t{}\t{pos}\t{state}\t99.0", c.name)?;
            sites += 1;
            pos += rng.range(1, 1_040);
        }
        for line in &rg {
            w.write_all(line.as_bytes())?;
        }
        runs += rg.len() as u64;
    }
    w.flush()?;
    eprintln!(
        "genome_scale: {sites} ST lines, {runs} RG rows in {}",
        dir.display()
    );
    Ok(())
}
