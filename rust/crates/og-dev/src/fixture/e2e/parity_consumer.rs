use super::{Fail, create, exactly};
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::Path;

const BASES: &[u8; 4] = b"ACGT";
const CONTIGS: [(&str, &str, usize); 5] = [
    ("chr1", "1", 10_800_000),
    ("chr2", "2", 8_400_000),
    ("chrX", "X", 6_000_000),
    ("chrY", "Y", 1_800_000),
    ("chrM", "MT", 16_569),
];

fn base_index(contig: &str, pos: usize) -> usize {
    ((pos as u64 * 2_654_435_761 + contig.len() as u64) % 7 % 4) as usize
}

/// The fixed-seed gVCF and consumer template used by the parity case.
pub(super) fn run(args: &[String]) -> Result<(), Fail> {
    let [vcf, template] = exactly(args, 2, "GVCF.gz SITES.tsv.gz")? else {
        unreachable!()
    };
    let mut vcf = og_io::bgzf::Writer::new(create(Path::new(vcf))?, NonZeroUsize::MIN)?;
    let mut sites = Vec::new();
    let mut state = 12345u64;
    let mut random = |n| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 33) as usize) % n
    };

    writeln!(vcf, "##fileformat=VCFv4.2")?;
    for (contig, _, _) in CONTIGS {
        writeln!(vcf, "##contig=<ID={contig}>")?;
    }
    writeln!(
        vcf,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
    )?;
    let mut site_count = 0;
    for (contig, label, size) in CONTIGS {
        let mut pos = 1 + random(50);
        while pos < size {
            let reference = base_index(contig, pos);
            let base = char::from(BASES[reference]);
            let kind = random(100);
            let haploid = matches!(contig, "chrY" | "chrM") && random(2) != 0;
            match kind {
                0..62 => {
                    let end = size.min(pos + random(40));
                    let genotype = if haploid { "0" } else { "0/0" };
                    writeln!(
                        vcf,
                        "{contig}\t{pos}\t.\t{base}\t<*>\t0\t.\tEND={end}\tGT:GQ:MIN_DP\t{genotype}:30:20"
                    )?;
                    pos = end + 1;
                }
                62..85 => {
                    let alt = char::from(BASES[(reference + 1 + random(3)) % 4]);
                    let genotype = if haploid {
                        "1"
                    } else {
                        ["0/1", "1/1", "0|1"][random(3)]
                    };
                    let filter = if random(20) == 0 { "LowQual" } else { "PASS" };
                    writeln!(
                        vcf,
                        "{contig}\t{pos}\t.\t{base}\t{alt},<*>\t40\t{filter}\t.\tGT:GQ\t{genotype}:40"
                    )?;
                    if random(5) == 0 {
                        writeln!(
                            vcf,
                            "{contig}\t{pos}\t.\t{base}\t<*>\t0\t.\tEND={}\tGT:GQ\t0/0:20",
                            pos + 3
                        )?;
                        pos += 3;
                    }
                    pos += 1;
                }
                85..93 => {
                    let length = 2 + random(4);
                    let reference: String = (pos..pos + length)
                        .map(|p| char::from(BASES[base_index(contig, p)]))
                        .collect();
                    let genotype = ["0/1", "1/1", "0/0"][random(3)];
                    let filter = if genotype == "0/0" { "RefCall" } else { "PASS" };
                    writeln!(
                        vcf,
                        "{contig}\t{pos}\t.\t{reference}\t{base},<*>\t30\t{filter}\t.\tGT:GQ\t{genotype}:30"
                    )?;
                    let next_base = char::from(BASES[base_index(contig, pos + 1)]);
                    writeln!(
                        vcf,
                        "{contig}\t{}\t.\t{next_base}\t<*>\t0\t.\tEND={}\tGT:GQ\t0/0:20",
                        pos + 1,
                        pos + length + 2
                    )?;
                    pos += length + 3;
                }
                93..97 => {
                    let end = size.min(pos + random(10));
                    writeln!(
                        vcf,
                        "{contig}\t{pos}\t.\t{base}\t<*>\t0\t.\tEND={end}\tGT:GQ\t./.:0"
                    )?;
                    pos = end + 1;
                }
                _ => pos += 1 + random(30),
            }
        }
        for pos in (1 + random(5))..size {
            if random(150) != 0 {
                continue;
            }
            let reference = base_index(contig, pos);
            let base = char::from(BASES[(reference + usize::from(random(400) == 0)) % 4]);
            let strand = if random(4) == 0 { '-' } else { '+' };
            site_count += 1;
            writeln!(
                sites,
                "rs{site_count}\t{label}\t{}\t{base}\t{pos}\t{base}\t{strand}",
                pos + 10_000
            )?;
        }
    }
    vcf.finish()?.flush()?;
    let mut template =
        flate2::write::GzEncoder::new(create(Path::new(template))?, flate2::Compression::default());
    writeln!(
        template,
        "## template_sites: {}\n## dropped: 1000\n#rsid\tchrom\tpos37\tref37\tpos38\tref38\tstrand",
        site_count + 1000
    )?;
    template.write_all(&sites)?;
    template.finish()?.flush()?;
    Ok(())
}
