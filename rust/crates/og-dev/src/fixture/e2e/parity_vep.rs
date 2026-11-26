use super::{Fail, PyRandom, create, exactly};
use pycompat::fmt;
use std::io::Write;
use std::path::Path;

const FIELDS: &str = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|SIFT|PolyPhen|MAX_AF|gnomADe_AF|gnomADg_AF|CLIN_SIG";

/// The original fixed-seed 500,000-record fixture, including tied transcripts.
pub(super) fn run(args: &[String]) -> Result<(), Fail> {
    let [path] = exactly(args, 1, "OUT.vcf")? else {
        unreachable!()
    };
    let mut out = create(Path::new(path))?;
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(
        out,
        r#"##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {FIELDS}">"#
    )?;
    writeln!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
    )?;

    let mut rng = PyRandom::new(20261110);
    let impacts = ["MODIFIER", "LOW", "MODERATE", "HIGH"];
    for index in 0..500_000 {
        let mane = rng.below(9);
        write!(
            out,
            "chr1\t{}\t.\tA\tG\t50\tPASS\tAC=1;AN=2;CSQ=",
            1000 + index * 10
        )?;
        for transcript in 0..(1 + rng.below(8)) {
            let impact = rng.choice(&impacts);
            let canonical = if rng.random() < 0.4 { "YES" } else { "" };
            let mane = if transcript + 1 == mane || rng.random() < 0.05 {
                "NM_1.1"
            } else {
                ""
            };
            let af = fmt::fixed(rng.random() / 10.0, 5);
            let exome_af = fmt::general(rng.random() / 10.0, 4);
            let gene = format!("GENE{transcript}");
            let feature = format!("ENST{transcript:011}");
            let mut fields = [""; 32];
            fields[0] = "G";
            fields[1] = "intron_variant";
            fields[2] = impact;
            fields[3] = &gene;
            fields[6] = &feature;
            fields[23] = canonical;
            fields[24] = mane;
            fields[28] = &af;
            fields[29] = &exome_af;
            if transcript > 0 {
                write!(out, ",")?;
            }
            write!(out, "{}", fields.join("|"))?;
        }
        let quality = rng.below(61);
        let depth = rng.below(51);
        let reference_depth = rng.below(31);
        let alternate_depth = rng.below(31);
        let genotype = rng.choice(&["0/1", "1/1", "0|1", "1/2"]);
        writeln!(
            out,
            "\tGT:GQ:DP:AD\t{genotype}:{quality}:{depth}:{reference_depth},{alternate_depth}"
        )?;
    }
    out.flush()?;
    Ok(())
}
