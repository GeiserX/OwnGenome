use super::{Fail, PyRandom, create, exactly};
use flate2::{Compression, write::GzEncoder};
use std::fs;
use std::io::Write;
use std::path::Path;

const SAMPLES: [&str; 5] = ["DAD", "KID", "SIB", "MUM", "OTHER"];
const PAIRS: [(&str, &str, usize, usize); 5] = [
    ("DAD", "KID", 1, 1),
    ("DAD", "SIB", 2, 1),
    ("KID", "SIB", 1, 1),
    ("KID", "SIB", 2, 2),
    ("MUM", "KID", 1, 2),
];
const GENOTYPES: [&str; 4] = ["0|0", "0|1", "0|1", "1|1"];
const MB: usize = 1_000_000;
const FIRST: usize = MB;
const MARKERS: usize = 50_000;
const STEP: usize = 1000;
const LAST: usize = FIRST + STEP * (MARKERS - 1);

/// The fixed-seed 22-chromosome IBD speed fixture.
pub(super) fn run(args: &[String]) -> Result<(), Fail> {
    let [dir] = exactly(args, 1, "DIRECTORY")? else {
        unreachable!()
    };
    let dir = Path::new(dir);
    for part in ["maps", "phased", "ibd"] {
        fs::create_dir_all(dir.join(part))?;
    }
    fs::write(
        dir.join("pedigree.tsv"),
        "sample\tsex\tfather\tmother\tpartner\nKID\tmale\tDAD\tMUM\t.\nSIB\tfemale\tDAD\tMUM\t.\n",
    )?;

    let mut rng = PyRandom::new(38);
    let mut chroms = create(&dir.join("chroms.txt"))?;
    for number in 1..=22 {
        let chrom = format!("chr{number}");
        writeln!(chroms, "{chrom}")?;
        let mut map = create(&dir.join("maps").join(format!("{chrom}.map")))?;
        for i in 0..=52 {
            writeln!(map, "{chrom}\t.\t{i}.0000\t{}", (i * MB).max(1))?;
        }
        map.flush()?;

        let path = dir.join("ibd").join(format!("{chrom}.ibd.gz"));
        let mut ibd = GzEncoder::new(create(&path)?, Compression::default());
        for (a, b, a_hap, b_hap) in PAIRS {
            let mut start = FIRST;
            while start < FIRST + STEP * MARKERS - 5 * MB {
                let end = (start + (8 + rng.below(23)) * MB).min(LAST);
                let cm = (end - start) as f64 / MB as f64;
                writeln!(
                    ibd,
                    "{a}\t{a_hap}\t{b}\t{b_hap}\t{chrom}\t{start}\t{end}\t{}",
                    pycompat::fmt::repr(cm)
                )?;
                start = end + rng.below(4) * MB;
            }
        }
        ibd.finish()?.flush()?;

        // The old fixture draws all markers of each sample before the next.
        let genotypes = SAMPLES.map(|_| {
            (0..MARKERS)
                .map(|_| rng.below(GENOTYPES.len()) as u8)
                .collect::<Vec<_>>()
        });
        let path = dir.join("phased").join(format!("{chrom}.vcf.gz"));
        let mut vcf = GzEncoder::new(create(&path)?, Compression::default());
        writeln!(
            vcf,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
            SAMPLES.join("\t")
        )?;
        for marker in 0..MARKERS {
            write!(
                vcf,
                "{chrom}\t{}\t.\tA\tG\t.\tPASS\t.\tGT",
                FIRST + marker * STEP
            )?;
            for sample in &genotypes {
                write!(vcf, "\t{}", GENOTYPES[sample[marker] as usize])?;
            }
            writeln!(vcf)?;
        }
        vcf.finish()?.flush()?;
    }
    chroms.flush()?;
    Ok(())
}
