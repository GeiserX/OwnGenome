//! The setup components: one row per flag of `scripts/setup.sh`, the one
//! table `og setup install`, `og setup flags` and `og setup ported` read,
//! and that `og setup check` and `og dev everything-check` read too.
//!
//! A row says what the flag installs (`what`), for which steps, which paths
//! mean it is installed (`installed`, relative to the install folder, with
//! `${NAME}` from versions.env and `{ref}`, `{ref_stem}` and `{ref_name}`
//! from the reference), the setup.sh arguments that retry it, its licence
//! notices and the images it pulls. Its
//! downloads are the rows of `assets/data-sources.tsv` whose installer column
//! names its key ([`Component::sources`]).
//!
//! Whether og runs a flag yet is not a column: the install module that runs
//! it lists it in its `PORTED` ([`crate::install::ported`]), so a later
//! slice of the port edits only its own module. A flag no module lists is
//! run by `scripts/setup.sh` with `OG_SETUP=bash`.

use crate::fetch::{Ctx, Source, Sources};
use std::path::{Path, PathBuf};

/// What a flag is: an installer of data, or a way to run setup.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    /// Installs one component and exits.
    Installer,
    /// Changes how a whole setup runs (`--everything`, `--no-clinical-data`,
    /// `--refresh`, `--og`, `--pull-only`).
    Mode,
}

/// One flag of setup.sh.
#[derive(Debug, Clone, Copy)]
pub struct Component {
    /// The flag without its dashes, as the table's installer column names it.
    pub key: &'static str,
    pub flag: &'static str,
    pub kind: Kind,
    /// The steps that read what it installs ("" for a mode).
    pub steps: &'static str,
    pub what: &'static str,
    /// Paths that mean it is installed, relative to the install folder.
    pub installed: &'static [&'static str],
    /// The setup.sh arguments that retry it, before the install folder.
    pub retry: &'static str,
    /// The licence notices setup prints before it installs.
    pub terms: &'static [&'static str],
    /// The versions.env variables of the images it pulls.
    pub images: &'static [&'static str],
}

pub const CYRIUS_TERMS: &str = "Cyrius: PolyForm Strict 1.0.0, non-commercial use only.";
pub const PRS_ACADEMIC_TERMS: [&str; 4] = [
    "PGS000014: freely available to the academic community for research use.",
    "PGS000016: freely available to the academic community for research use.",
    "PGS000017: freely available to the academic community for research use.",
    "PGS000027: freely available to the academic community for research use.",
];
pub const ANNOTSV_TERMS: &str =
    "AnnotSV's annotation bundle includes third-party data with their own terms.";
pub const EXOMISER_TERMS: &str =
    "Exomiser's hg38 and phenotype bundles include third-party data with their own terms.";
pub const UCSC_CHAIN_TERMS: &str = "UCSC allows free non-commercial use of its hg19ToHg38 chain by independent researchers and nonprofit organisations. Other uses require a commercial licence.";
pub const ANNOTATION_SCORE_TERMS: [&str; 4] = [
    "CADD: free for non-commercial and academic use. Commercial use requires a licence from the University of Washington.",
    "SpliceAI precomputed scores: academic and not-for-profit use only. Other use requires a commercial licence from Illumina.",
    "REVEL: the authors allow non-commercial use. Zenodo lists ODbL 1.0. Contact the maintainers for other uses.",
    "AlphaMissense: CC BY-NC-SA 4.0 (non-commercial, share-alike).",
];
pub const PCGR_TERMS: &str =
    "PCGR/CPSR reference data includes third-party databases with their own terms.";

const fn installer(
    key: &'static str,
    flag: &'static str,
    steps: &'static str,
    what: &'static str,
    installed: &'static [&'static str],
    images: &'static [&'static str],
) -> Component {
    Component {
        key,
        flag,
        kind: Kind::Installer,
        steps,
        what,
        installed,
        retry: flag,
        terms: &[],
        images,
    }
}

const fn mode(key: &'static str, flag: &'static str, what: &'static str) -> Component {
    Component {
        key,
        flag,
        kind: Kind::Mode,
        steps: "",
        what,
        installed: &[],
        retry: flag,
        terms: &[],
        images: &[],
    }
}

const fn terms(c: Component, t: &'static [&'static str]) -> Component {
    Component { terms: t, ..c }
}

/// Every flag of setup.sh, in the order of its usage text.
pub const TABLE: &[Component] = &[
    mode(
        "everything",
        "--everything",
        "a full setup, then every installer of assets/everything.tsv",
    ),
    mode(
        "og",
        "--og",
        "og's binaries in bin/.og (built from rust/, or a release's)",
    ),
    terms(
        installer(
            "cyrius",
            "--cyrius",
            "21",
            "Cyrius",
            &["tools/cyrius-${CYRIUS_VERSION}/INSTALLED"],
            &["PYTHON_IMAGE"],
        ),
        &[CYRIUS_TERMS],
    ),
    terms(
        installer(
            "prs-academic",
            "--prs-academic",
            "25",
            "the academic-only PGS Catalog scores",
            &[
                "prs_scores/PGS000014.txt.gz",
                "prs_scores/PGS000016.txt.gz",
                "prs_scores/PGS000017.txt.gz",
                "prs_scores/PGS000027.txt.gz",
            ],
            &[],
        ),
        &PRS_ACADEMIC_TERMS,
    ),
    installer(
        "parascopy-data",
        "--parascopy-data",
        "35",
        "Parascopy's homology table and models",
        &[
            "reference/parascopy-${PARASCOPY_DATA_VERSION}/homology_table/GRCh38.bed.gz",
            "reference/parascopy-${PARASCOPY_DATA_VERSION}/models_GRCh38_1KGP/EUR/SMN1.gz",
        ],
        &["PARASCOPY_IMAGE"],
    ),
    installer(
        "kir-data",
        "--kir-data",
        "08",
        "the IPD-KIR database",
        &["kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat"],
        &[],
    ),
    installer(
        "ancestry-panel",
        "--ancestry-panel",
        "25,26",
        "pgsc_calc's ancestry reference panel",
        &["reference/pgsc_calc/${PGSC_PANEL}.tar.zst"],
        &[],
    ),
    installer(
        "ancestry-fine-panel",
        "--ancestry-fine-panel",
        "26b",
        "the fine-grained ancestry panel",
        &[
            "reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/panel/DONE",
            "reference/ancestry_fine/${ANCESTRY_FINE_PANEL}/model/DONE",
        ],
        &["OG_IMAGE", "BCFTOOLS_HTTPS_IMAGE", "BEAGLE_IMAGE"],
    ),
    installer(
        "yleaf-data",
        "--yleaf-data",
        "37",
        "Yleaf's Y marker tables and trees",
        &[
            "reference/yleaf-${YLEAF_DATA_VERSION}/data/hg38/new_positions.txt",
            "reference/yleaf-${YLEAF_DATA_VERSION}/data/hg_prediction_tables/tree.json",
        ],
        &["YLEAF_IMAGE"],
    ),
    installer(
        "consumer-sites",
        "--consumer-sites",
        "40",
        "the site table of the raw-data export",
        &["reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"],
        &[],
    ),
    installer(
        "ancient-panel",
        "--ancient-panel",
        "44",
        "the ancient populations panel",
        &["reference/ancient/aadr_${AADR_VERSION}/DONE"],
        &["EIGENSOFT_IMAGE"],
    ),
    installer(
        "viral-data",
        "--viral-data",
        "47",
        "the viral genomes",
        &["reference/viral-${VIRAL_REF_VERSION}/viral_genomes.fa"],
        &[],
    ),
    installer(
        "masked-reference",
        "--masked-reference",
        "02,02a",
        "the reference with GRCh38's false duplications masked",
        &[
            "{ref_stem}_maskedGRC.fasta",
            "{ref_stem}_maskedGRC.fasta.fai",
            "{ref_stem}_maskedGRC.dict",
        ],
        &["OG_IMAGE"],
    ),
    terms(
        installer(
            "exomiser-data",
            "--exomiser-data",
            "48",
            "Exomiser's hg38 and phenotype data",
            &[
                "reference/exomiser-${EXOMISER_DATA_VERSION}/${EXOMISER_DATA_VERSION}_hg38",
                "reference/exomiser-${EXOMISER_DATA_VERSION}/${EXOMISER_DATA_VERSION}_phenotype",
            ],
            &["EXOMISER_IMAGE"],
        ),
        &[EXOMISER_TERMS],
    ),
    installer(
        "chip-data",
        "--chip-data",
        "49",
        "Ensembl's GRCh38 gene annotation",
        &["reference/ensembl/Homo_sapiens.GRCh38.${ENSEMBL_GFF_RELEASE}.chr.gff3.gz"],
        &[],
    ),
    installer(
        "strchive-data",
        "--strchive-data",
        "09,9b",
        "the STRchive repeat catalog",
        &[
            "reference/strchive-${STRCHIVE_VERSION}/catalog.json",
            "reference/strchive-${STRCHIVE_VERSION}/str_ranges.tsv",
            "reference/strchive-${STRCHIVE_VERSION}/catalog.stamp",
        ],
        &["OG_IMAGE"],
    ),
    installer(
        "megane-kmers",
        "--megane-kmers",
        "50",
        "MEGAnE's k-mer set of the reference",
        &[
            "reference/megane/{ref_name}.mk",
            "reference/megane/{ref_name}.mi",
        ],
        &["MEGANE_IMAGE"],
    ),
    installer(
        "mocha-data",
        "--mocha-data",
        "51",
        "MoChA's phasing panel",
        &["reference/mocha-${MOCHA_PANEL_VERSION}/DONE"],
        &["MOCHA_IMAGE", "BEAGLE_IMAGE"],
    ),
    installer(
        "locityper-data",
        "--locityper-data",
        "53",
        "Locityper's database",
        &[
            "reference/locityper-${LOCITYPER_PANEL_VERSION}/DONE",
            "reference/locityper-${LOCITYPER_PANEL_VERSION}/counts.jf",
            "reference/locityper-${LOCITYPER_PANEL_VERSION}/loci.bed",
            "reference/locityper-${LOCITYPER_PANEL_VERSION}/db",
        ],
        &["LOCITYPER_IMAGE"],
    ),
    installer(
        "vep-cache",
        "--vep-cache",
        "13",
        "the VEP cache",
        &["vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt"],
        &[],
    ),
    terms(
        installer(
            "annotsv-data",
            "--annotsv-data",
            "05",
            "AnnotSV's data, gnomAD-SV's common SVs and the MANE Select list",
            &[
                "annotsv_annotations/Annotations_Human/Genes/GRCh38",
                "annotsv_annotations/mane_select_refseq.txt",
                "annotsv_annotations/mane_select_refseq.release",
            ],
            &[],
        ),
        &[ANNOTSV_TERMS],
    ),
    terms(
        installer(
            "annotation-scores",
            "--annotation-scores",
            "13",
            "CADD, masked SpliceAI, REVEL and AlphaMissense",
            &["annotation_scores"],
            &[],
        ),
        &ANNOTATION_SCORE_TERMS,
    ),
    installer(
        "somatic-data",
        "--somatic-data",
        "29,49",
        "GATK's somatic resources",
        &[
            "somatic/af-only-gnomad.hg38.vcf.gz",
            "somatic/af-only-gnomad.hg38.vcf.gz.tbi",
            "somatic/1000g_pon.hg38.vcf.gz",
            "somatic/1000g_pon.hg38.vcf.gz.tbi",
            "somatic/small_exac_common_3.hg38.vcf.gz",
            "somatic/small_exac_common_3.hg38.vcf.gz.tbi",
        ],
        &[],
    ),
    installer(
        "bwa-index",
        "--bwa-index",
        "04b",
        "the classic BWA index of the reference",
        &[
            "{ref}.amb",
            "{ref}.ann",
            "{ref}.bwt",
            "{ref}.pac",
            "{ref}.sa",
        ],
        &["BWA_IMAGE"],
    ),
    installer(
        "cnvpytor-data",
        "--cnvpytor-data",
        "18",
        "CNVpytor's reference files",
        &[
            "reference/cnvpytor/gc_hg38.pytor",
            "reference/cnvpytor/mask_hg38.pytor",
            "reference/cnvpytor/gc_hg19.pytor",
            "reference/cnvpytor/mask_hg19.pytor",
            "reference/cnvpytor/gc_chm13v2.0.pytor",
            "reference/cnvpytor/gc_chm13v1.1.pytor",
            "reference/cnvpytor/gc_kn99.pytor",
        ],
        &[],
    ),
    installer(
        "pypgx-data",
        "--pypgx-data",
        "32",
        "the PyPGx bundle",
        &["reference/pypgx-bundle"],
        &[],
    ),
    terms(
        installer(
            "archaic",
            "--archaic",
            "43",
            "the archaic genomes and 1000 Genomes frequencies",
            &["reference/archaic/${ARCHAIC_DATA_VERSION}/DONE"],
            &["BCFTOOLS_IMAGE", "PICARD_IMAGE", "OG_IMAGE", "IBDMIX_IMAGE"],
        ),
        &[UCSC_CHAIN_TERMS],
    ),
    terms(
        installer(
            "clinical-data",
            "--clinical-data",
            "13,17",
            "VEP 116, the PCGR/CPSR bundle and VEP 115",
            &[
                "vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt",
                "pcgr_data/${PCGR_DATA_BUNDLE}/data",
            ],
            &[],
        ),
        &[PCGR_TERMS],
    ),
    mode(
        "no-clinical-data",
        "--no-clinical-data",
        "a full setup without the VEP caches and the PCGR/CPSR bundle",
    ),
    mode(
        "pull-only",
        "--pull-only",
        "the Docker images of versions.env, nothing else",
    ),
    installer(
        "sample-qc-data",
        "--sample-qc-data",
        "33",
        "somalier's sites and VerifyBamID2's panel",
        &[
            "reference/somalier/sites.hg38.vcf.gz",
            "reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.UD",
            "reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.mu",
            "reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.bed",
        ],
        &[],
    ),
    mode(
        "refresh",
        "--refresh",
        "ClinVar replaced with NCBI's current release (--refresh clinvar)",
    ),
];

/// The row of FLAG (with its dashes).
pub fn by_flag(flag: &str) -> Option<&'static Component> {
    TABLE.iter().find(|c| c.flag == flag)
}

/// The row of KEY (a flag without its dashes).
pub fn by_key(key: &str) -> Option<&'static Component> {
    TABLE.iter().find(|c| c.key == key)
}

/// Every flag, in the table's order.
pub fn flags() -> impl Iterator<Item = &'static str> {
    TABLE.iter().map(|c| c.flag)
}

impl Component {
    /// og runs it (`og setup install`); otherwise setup.sh's bash does.
    pub fn ported(&self) -> bool {
        crate::install::ported(self.key)
    }

    /// The rows of the source table this component downloads.
    pub fn sources<'a>(&self, s: &'a Sources) -> Vec<&'a Source> {
        s.rows
            .values()
            .filter(|r| r.installers.iter().any(|i| i == self.key))
            .collect()
    }

    /// The installed paths, expanded: `${NAME}` from versions.env and
    /// `{ref}`, `{ref_stem}`, `{ref_name}` from REF (the reference FASTA),
    /// each under GENOME unless it names the reference.
    pub fn installed_paths(
        &self,
        ctx: &Ctx,
        genome: &Path,
        reference: &Path,
    ) -> Result<Vec<PathBuf>, String> {
        let r = reference.to_string_lossy().into_owned();
        let stem = match r.rfind('.') {
            Some(i) if !r[i..].contains('/') => r[..i].to_owned(),
            _ => r.clone(),
        };
        let name = {
            let base = reference
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default();
            let base = base.strip_suffix(".gz").unwrap_or(&base).to_owned();
            match base.rfind('.') {
                Some(i) => base[..i].to_owned(),
                None => base,
            }
        };
        self.installed
            .iter()
            .map(|t| {
                let e = ctx.expand(t, None)?;
                Ok(if let Some(rest) = e.strip_prefix("{ref_stem}") {
                    PathBuf::from(format!("{stem}{rest}"))
                } else if let Some(rest) = e.strip_prefix("{ref}") {
                    PathBuf::from(format!("{r}{rest}"))
                } else {
                    genome.join(e.replace("{ref_name}", &name))
                })
            })
            .collect()
    }

    /// The command that retries it: setup.sh (SCRIPT) with its arguments and
    /// the install folder, each shell-quoted as `printf %q` quotes it.
    pub fn retry_command(&self, script: &str, genome: &str) -> String {
        format!(
            "{} {} {}",
            crate::install::shell_quote(script),
            self.retry,
            crate::install::shell_quote(genome)
        )
    }
}
