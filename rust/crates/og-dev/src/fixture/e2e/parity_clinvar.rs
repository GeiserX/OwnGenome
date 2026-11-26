use super::{Fail, PyRandom, create, exactly};
use std::io::Write;
use std::path::Path;

const CONTIGS: [&str; 3] = ["s1", "s2", "s3"];
const LENGTH: usize = 300_000;
const WIDTH: usize = 70;
const BASES: &[u8] = b"ACGT";
const GENOTYPES: [&str; 19] = [
    "0/1", "0/1", "0/1", "0/1", "0/1", "0/1", "1/0", "1/1", "./1", "1/.", "0/0", "0|1", "1|0",
    "0|1", "1|0", "1|.", ".|1", "1|1", "2/0",
];

/// The fixed-seed reference and dense variant calls for adjacent allele tests.
pub(super) fn run(args: &[String]) -> Result<(), Fail> {
    let [fasta, index, vcf] = exactly(args, 3, "FASTA FAI VCF")? else {
        unreachable!()
    };
    let mut fasta = create(Path::new(fasta))?;
    let mut index = create(Path::new(index))?;
    let mut vcf = create(Path::new(vcf))?;
    let mut rng = PyRandom::new(20261111);
    let mut sequences = Vec::with_capacity(CONTIGS.len());
    let mut offset = 0;

    for contig in CONTIGS {
        let mut sequence = rng.bases(LENGTH);
        for _ in 0..40 {
            let start = rng.below(LENGTH - 2000);
            let end = start + 50 + rng.below(1450);
            sequence[start..end].make_ascii_lowercase();
        }
        for _ in 0..5 {
            let start = rng.below(LENGTH - 500);
            let end = start + 10 + rng.below(290);
            sequence[start..end].fill(b'N');
        }
        let header = format!(">{contig} synthetic\n");
        fasta.write_all(header.as_bytes())?;
        offset += header.len();
        writeln!(
            index,
            "{contig}\t{LENGTH}\t{offset}\t{WIDTH}\t{}",
            WIDTH + 1
        )?;
        for row in sequence.chunks(WIDTH) {
            fasta.write_all(row)?;
            fasta.write_all(b"\n")?;
            offset += row.len() + 1;
        }
        sequence.make_ascii_uppercase();
        sequences.push(sequence);
    }

    writeln!(
        vcf,
        "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">"
    )?;
    for contig in CONTIGS {
        writeln!(vcf, "##contig=<ID={contig},length={LENGTH}>")?;
    }
    writeln!(
        vcf,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
        ##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG"
    )?;
    let steps = [1, 1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 15, 20];
    for (contig, sequence) in CONTIGS.into_iter().zip(sequences) {
        let mut pos = 100;
        let mut phase_set = 100;
        while pos < LENGTH - 100 {
            pos += rng.choice(&steps);
            if rng.random() < 0.02 {
                phase_set = pos;
            }
            let kind = rng.random();
            let base = sequence[pos - 1];
            let (mut reference, alternate) = if kind < 0.6 {
                let other: Vec<_> = BASES.iter().copied().filter(|b| *b != base).collect();
                (vec![base], vec![rng.choice(&other)])
            } else if kind < 0.8 {
                let length = 2 + rng.below(6);
                (sequence[pos - 1..pos - 1 + length].to_vec(), vec![base])
            } else {
                let length = 1 + rng.below(5);
                let mut alt = vec![base];
                alt.extend(rng.bases(length));
                (vec![base], alt)
            };
            if rng.random() < 0.05 {
                reference[0] = rng.choice(BASES);
            }
            if rng.random() < 0.03 {
                reference.make_ascii_lowercase();
            }
            let genotype = rng.choice(&GENOTYPES);
            let phased = genotype.contains('|') && rng.random() < 0.8;
            let gq = rng.below(60);
            let dp = rng.below(40);
            let quality = 1 + rng.below(59);
            write!(vcf, "{contig}\t{pos}\t.\t")?;
            vcf.write_all(&reference)?;
            vcf.write_all(b"\t")?;
            vcf.write_all(&alternate)?;
            write!(vcf, "\t{quality}\tPASS\t.\tGT:GQ:DP")?;
            if phased {
                writeln!(vcf, ":PS\t{genotype}:{gq}:{dp}:{phase_set}")?;
            } else {
                writeln!(vcf, "\t{genotype}:{gq}:{dp}")?;
            }
        }
    }
    fasta.flush()?;
    index.flush()?;
    vcf.flush()?;
    Ok(())
}
