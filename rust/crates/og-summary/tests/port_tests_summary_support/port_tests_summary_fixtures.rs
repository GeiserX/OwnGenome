//! Synthetic sample folders used by the collector and report integration tests.

use og_summary::collect::{Group, Input, post};
use pycompat::json::Value;
use std::fs::{self, File, FileTimes};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Duration, UNIX_EPOCH};

pub const REPO: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../..");
pub const VCF: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
pub const CLINVAR: &str = "clinvar/S_clinvar_hits.vcf";
pub const COUPLE: &str = "couple/B_couple_carrier.tsv";
pub const COUPLE_HEADER: &str = "gene\tinheritance\tlists\tcondition\tresult\toffspring_risk\tpartner_a\tfindings_a\tpartner_b\tfindings_b\texplanation\n";
pub const PRS: &str = "prs/S_prs_summary.tsv";
pub const ROH: &str = "vcf/S_roh.txt";
pub const SLIVAR: &str = "slivar/S_slivar_summary.tsv";
pub const PAIRS: &str = "slivar/S_compound_het_pairs.tsv";
pub const CPIC: &str = "cpic/S_phenotypes.tsv";
pub const PYP: &str = "pypgx/S_pypgx_summary.tsv";
pub const CYRIUS: &str = "cyrius/S_cyp2d6.tsv";
pub const CONSENSUS: &str = "pgx_consensus/S_pgx_consensus.tsv";
pub const RUN: &str = "logs/run_status.tsv";
pub const EXPORT: &str = "consumer_export/S_consumer_export.tsv";
pub const QC: &str = "qc/S_sample_qc.tsv";
pub const STR: &str = "expansion_hunter/S_eh.vcf";
pub const STRANGER: &str = "expansion_hunter/S_eh_stranger.vcf";
pub const STR_CATALOG: &str = "expansion_hunter/S_eh_catalog.json";
pub const STARTED: u64 = 1_700_000_000;

pub struct Sample {
    dir: PathBuf,
    pub id: &'static str,
}

impl Sample {
    pub fn new(case: &str, id: &'static str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("og-summary-port-{}-{case}", std::process::id()));
        fs::create_dir(&dir).expect("create the synthetic sample folder");
        Self { dir, id }
    }

    pub fn path(&self, name: &str) -> PathBuf {
        self.dir.join(name)
    }

    pub fn put(&self, name: &str, body: impl AsRef<[u8]>) {
        let path = self.path(name);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, body).unwrap();
    }

    pub fn gzip(&self, name: &str, body: &str) {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(body.as_bytes()).unwrap();
        self.put(name, encoder.finish().unwrap());
    }

    pub fn old(&self, name: &str) {
        File::open(self.path(name))
            .unwrap()
            .set_times(
                FileTimes::new().set_modified(UNIX_EPOCH + Duration::from_secs(STARTED - 60)),
            )
            .unwrap();
    }

    pub fn collect(&self) -> Value {
        post::collect(&Input {
            sample: self.id.into(),
            sample_dir: self.dir.display().to_string(),
            declared_sex: None,
            assets: Path::new(REPO).join("assets").display().to_string(),
            group: Group::All,
            tables: Path::new(REPO).join("bin").display().to_string(),
        })
        .expect("collect the synthetic sample")
    }
}

impl Drop for Sample {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.dir);
    }
}

pub fn put_bash(d: &Sample) {
    d.gzip(
        "vcf/S.vcf.gz",
        &format!("{VCF}chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\nchr1\t2\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1\nchr1\t3\t.\tC\tT\t5\tRefCall\t.\tGT\t0/0\n"),
    );
    d.put(CLINVAR, format!(
        "{VCF}chr1\t10\t1\tC\tT\t50\tPASS\tGENEINFO=GENEA:11;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\tGT\t0/1\n\
         chr1\t20\t2\tT\tG\t50\tPASS\tGENEINFO=GENEB:22;CLNSIG=Likely_pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT\t1/1\n\
         chr2\t5\t3\tA\tG\t50\tPASS\tGENEINFO=GENEC:33;CLNSIG=Pathogenic;CLNREVSTAT=no_assertion_criteria_provided\tGT\t0|1\n"
    ));
    let mut mito = VCF.to_owned();
    for (i, af) in ["0.02", "0.30", "0.99", "0.05", "0.949"]
        .into_iter()
        .enumerate()
    {
        use std::fmt::Write;
        writeln!(
            mito,
            "chrM\t{}\t.\tA\tG\t.\tPASS\t.\tGT:AF\t0/1:{af}",
            i + 1
        )
        .unwrap();
    }
    d.gzip("mito/S_chrM_filtered.vcf.gz", &mito);
    d.put(
        ROH,
        "RG\tS\tchr1\t100\t7000100\t7000000\t3\t50\nRG\tS\tchrX\t1\t9000000\t9000000\t3\t50\n",
    );
    d.put("mosdepth/S.mosdepth.summary.txt", "chrom\tlength\tbases\tmean\tmin\tmax\nchr1\t10\t300\t30.0\t0\t50\ntotal\t10\t300\t29.5\t0\t50\n");
    d.put("indexcov/indexcov-indexcov.ped", "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\nS\tS\t-9\t-9\t2\t-9\t2.0\t0.0\n");
    d.put(PRS, "Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\nCoronary artery disease\tPGS000018\t0.12\t100\t200\nType 2 diabetes\tPGS000014\t-0.4\t50\t90\n");
    d.gzip("cpsr/S.cpsr.grch38.classification.tsv.gz", "SAMPLE_ID\tGENOMIC_CHANGE\tCPSR_CLASSIFICATION\tCLASSIFICATION\nS\ta\tVUS\tVUS\nS\tb\tVUS\tLikely_Pathogenic\nS\tc\tB\tVUS\n");
    d.put(
        SLIVAR,
        "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\nchr1\t1\tA\tG\tHIGH\tX\n",
    );
    d.put(
        "hla_t1k/S_hla_genotype.tsv",
        "HLA-A\t2\tA*01:01:01\t30\t60\tA*02:01:01\t20\t55\t\n",
    );
    d.put(CPIC, "Gene\tDiplotype\tPhenotype\tStatus\nCYP2D6\t*1/*4\tIntermediate Metabolizer\tnon-normal\nCYP2C19\t*1/*1\tNormal Metabolizer\tnormal\nVKORC1\trs9923231 reference (C)/rs9923231 reference (C)\t-1639 GG\tunclassified\n");
    d.put(
        PYP,
        "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\nCYP2D6\t*4/*1\tIM\t.\tbam\n",
    );
    put_cyrius(d, "*1/*4", "PASS");
}

pub fn put_cyrius(d: &Sample, diplotype: &str, filter: &str) {
    d.put(
        CYRIUS,
        format!(
            "Sample\tGenotype\tFilter\n{}\t{diplotype}\t{filter}\n",
            d.id
        ),
    );
}

pub fn put_nextflow(d: &Sample) {
    d.put("roh/S_roh.txt", "RG\tS\tchr2\t1\t6000001\t6000000\t3\t50\n");
    d.put("roh/S_roh_summary.txt", "summary\n");
    d.put(
        "coverage/S.mosdepth.summary.txt",
        "chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t1\t1\t31.25\t0\t9\n",
    );
    d.put(
        "hla/S_hla_genotype.tsv",
        "HLA-B\t2\tB*07:02:01\t30\t60\tB*08:01:01\t20\t55\t\n",
    );
    d.put(
        "pharmcat/S.report.json",
        r#"{"pharmcatVersion":"3.2.0","genes":{}}"#,
    );
}

pub fn put_stale(d: &Sample) {
    put_bash(d);
    d.old(SLIVAR);
    d.old(ROH);
    d.put(RUN, format!("meta\tstarted_epoch\t{STARTED}\nmeta\tdeclared_sex\tmale\nstep\t31\tskipped (needs VEP, step 13 skipped)\nstep\t11\tok\n"));
}

pub fn put_failed(d: &Sample) {
    d.put(
        RUN,
        format!(
            "meta\tstarted_epoch\t{STARTED}\nstep\t10\tfailed\nstep\t11\tok\nstep\t32\tfailed\n"
        ),
    );
}

pub fn put_couple(d: &Sample) {
    d.put(COUPLE, format!(
        "{COUPLE_HEADER}NAGA\tAR\tACMG\tNAGA disease\tshared\t1 in 4 per pregnancy\tA\ta-snv\tB\tb-snv\tx\n\
         FMR1\tXL\tACMG\tFMR1 disease\tx_linked\tA: sons 1 in 2\tA\ta-repeat\tB\t.\tx\n\
         UGT1A1\tAR\tACMG\tUGT1A1 disease\tone_partner\tnone shared\tA\t.\tB\tb-only\tx\n"
    ));
}

pub fn put_comphet(d: &Sample, nextflow: bool) {
    let pair_header = "GENE\tCHROM\tPOS_1\tREF_1\tALT_1\tPOS_2\tREF_2\tALT_2\tDISTANCE\tPHASE\tCIS_READS\tTRANS_READS\tPHASE_NOTE\n";
    let trans = "G2\tchr1\t50\tA\tG\t90\tC\tT\t40\ttrans\t1\t19\tx\n";
    if nextflow {
        d.gzip(
            "slivar/S_prioritized.vcf.gz",
            &format!("{VCF}chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"),
        );
        d.put(PAIRS, format!("{pair_header}{trans}"));
    } else {
        d.put(
            SLIVAR,
            "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\nchr1\t1\tA\tG\tHIGH\tX\n",
        );
        d.put("slivar/S_compound_hets.tsv", "GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\tPHASE\nG1\tchr1\t10\tA\tG\tHIGH\tstop_gained\t0/1\t.\nG1\tchr1\t13\tA\tG\tHIGH\tstop_gained\t0/1\t.\nG2\tchr1\t50\tA\tG\tHIGH\tstop_gained\t0/1\t.\nG2\tchr1\t90\tA\tG\tHIGH\tstop_gained\t0/1\t.\n");
        d.put(PAIRS, format!("{pair_header}G1\tchr1\t10\tA\tG\t13\tC\tT\t3\tcis\t22\t0\tx\n{trans}G3\tchr1\t100\tA\tG\t50100\tC\tT\t50000\tunknown\t0\t0\tno read or read pair covers both positions\n"));
    }
}

pub fn put_cyp_unusable(d: &Sample) {
    d.put(PYP, "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\nCYP2D6\tIndeterminate\tIndeterminate\t.\tbam\nCYP2C19\t*1/*2\tIM\t.\tbam\nCYP2C9\tindeterminate\t.\t.\tbam\n");
    put_cyrius(d, "*5/*5", "CYP2D6_depth_unreliable");
    d.put(CONSENSUS, "Gene\tResult\tOutside_call\tReason\tEvidence\nHLA-A\tnot typed\tno\tHLA typing (step 08) did not run\t-\nCYP2D6\tindeterminate\tno\tthe CYP2D6 depth check failed\tpypgx: no call (Indeterminate)\n");
}

pub fn put_prs_evaluation(d: &Sample, group: &str, effect: &str) {
    d.put(PRS, format!("Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\tMatched_Pct\tPercentile\tAncestry_Group\tChrX\tAmbiguous\tApplies_To\tSex_Check\tEvaluated_In\tEffect\tInput\nTest condition\tPGS000018\t0.12\t100\t200\t50\t71\t{group}\t0\t0\tall\tapplies\tEUR\t{effect}\tgenotypes\n"));
}

pub fn put_lineage(d: &Sample) {
    d.put(
        "mito/S_haplogroup.txt",
        "\"SampleID\"\t\"Haplogroup\"\t\"Rank\"\n\"S\"\t\"H1a\"\t\"0.95\"\n",
    );
    d.put("mito/S_haplocheck.txt", "\"Sample\"\t\"Contamination Status\"\t\"Contamination Level\"\t\"Distance\"\n\"S\"\t\"YES\"\t\"0.12\"\t\"5\"\n");
    d.put("y_haplogroup/S_y_haplogroup.txt", "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\nS_sorted\tR-M269\tM269\t100\t42\t0.97\t1\t1\t0.97\n");
    d.put("y_haplogroup/S_y_tree.txt", "Sample_name\tHg\tLevel\tPath_derived\tPath_ancestral\tPath_no_call\tAncestral_in_terminal\tSNPs_called\tMin_depth\tTree\tPositions\tOutvoted_ancestral\nS\tR-FT12345\t24\t61\t2\t310\tFT999\t4210\t10\tYFull YTree 14.06.0\tYleaf 4.1.6 tables\tFT777\n");
    d.put("y_haplogroup/S_y_tree_path.tsv", "Level\tBranch\tSNPs\tOn_tables\tDerived\tAncestral\tNo_call\tAncestral_SNPs\n1\tA0-T\t4\t4\t3\t0\t1\t.\n24\tR-FT12345\t6\t5\t2\t1\t2\tFT999\n");
}

pub fn put_ungraded_repeat(d: &Sample) {
    let genotype = "0/1";
    let records = [
        ("chr1", 155188000, 'C', 9, "MUC1", (20, 25), ""),
        (
            "chr8",
            118366913,
            'A',
            40,
            "SAMD12",
            (1, 40),
            ";STR_STATUS=pre_mutation",
        ),
        ("chr9", 100, 'A', 9, "NOTINCAT", (5, 9), ""),
    ];
    for (name, graded) in [(STR, false), (STRANGER, true)] {
        use std::fmt::Write;
        let mut body = VCF.to_owned();
        for (chrom, pos, reference, alt, locus, (a, b), status) in records {
            let status = if graded { status } else { "" };
            writeln!(body, "{chrom}\t{pos}\t.\t{reference}\t<STR{alt}>\t.\tPASS\tEND=1;REPID={locus};VARID={locus}{status}\tGT:REPCN\t{genotype}:{a}/{b}").unwrap();
        }
        d.put(name, body);
        if !graded {
            d.old(name);
        }
    }
    d.put(STR_CATALOG, r#"[{"LocusId":"MUC1","NoThresholds":"the reference repeat is 1200 bases, longer than 300: short reads cannot size it"},{"LocusId":"SAMD12","NormalMax":1,"PathologicMin":105}]"#);
    d.put("expansion_hunter/S_str_motifs.tsv", "locus\tchrom\tstart\tend\trepeat_reads\tdominant_motif\tdominant_class\tdominant_share\tmotif_shares\tverdict\tnote\nSAMD12\tchr8\t118366812\t118366918\t5\tTGAAA\tpathogenic\t1.0\tTGAAA=1.0\tpathogenic_motif\tmotif TGAAA (pathogenic) in 5 of 5 repeat reads: confirm with repeat-primed PCR\n");
}

pub fn put_export(d: &Sample) {
    d.put(EXPORT, "key\tvalue\nsites\t10\ncalled\t7\nno_call\t3\nfrom_reference_blocks\t5\ntemplate_sites\t12\nfile_23andme\tS_23andme_v5.zip\nfile_ancestrydna\tS_ancestrydna_layout.zip\nno_call: no-call in the gVCF\t3\n");
}

pub fn put_mosaic(d: &Sample) {
    d.put(QC, "key\tvalue\nsample\tS\ninferred_sex\tmale\nsex_check\tok\nmosaic_sex\tmale\nmosaic_sex_from\treads\nmosaic_chrom\tchrY\nmosaic_depth_ratio\t0.351\nmosaic_depth_ratio_expected\t0.5\nmosaic_depth_blocks\t109\nmosaic_loss_depth\t0.298\nmosaic_loss_depth_low\t0.293\nmosaic_loss_depth_high\t0.303\nmosaic_x_het_sites\t\nmosaic_loss_af\t\nmosaic_loss_af_low\t\nmosaic_loss_af_high\t\nmosaic_floor\t0.100\nmosaic_call\tdetected\nmosaic_reason\tthe low end of the 95% range is above 10%: a share of the cells lost the Y\n");
}

pub fn put_cpsr(d: &Sample) {
    d.gzip(
        "cpsr/S.cpsr.grch38.classification.tsv.gz",
        "VAR_ID\tGENOMIC_CHANGE\tSYMBOL\tCLINVAR_CLASSIFICATION\tCPSR_CLASSIFICATION\tACMG_CODE\tASSERTION_AUTHORITY\tCLASSIFICATION\tDP_CONTROL\n\
         13_32340300_G_A\t13:g.32340300G>A\tBRCA2\tPathogenic\tBenign\tBA1\tClinVar\tPathogenic\t40\n\
         17_43045700_C_T\t17:g.43045700C>T\tBRCA1\tLikely Pathogenic\tLikely Benign\tBS1\tClinVar\tLikely Pathogenic\t6\n\
         2_47800000_A_G\t2:g.47800000A>G\tMSH2\t\tLikely Pathogenic\tPVS1,PM2\tCPSR\tLikely Pathogenic\t30\n\
         5_112800000_T_C\t5:g.112800000T>C\tAPC\tBenign\tPathogenic\tPVS1,PM2\tClinVar\tBenign\t30\n\
         3_100_A_G\t3:g.100A>G\tMLH1\tVUS\tVUS\tPM2\tCPSR\tVUS\t30\n",
    );
    d.put("cpsr/S.cpsr.grch38.html", "<html></html>\n");
    d.put(CLINVAR, format!("{VCF}##PGP_CALL_FLOOR=GQ=30,DP=10\n"));
    d.put("cpsr/S_cpsr_cdna_check.tsv", "chrom\tpos\tref\talt\talt_reads\talt_clipped\tref_reads\tref_clipped\tcdna_check\tcdna_note\n13\t32340300\tG\tA\t12\t0\t10\t0\tclear\t\n2\t47800000\tA\tG\t9\t8\t11\t0\tflagged\tpossible cDNA/transgene contamination: check the reads\n");
    d.put(
        "run_manifest.tsv",
        "data\tclinvar_file_date\t20261001\ndata\tcpsr_clinvar_release\t20250309\n",
    );
}

pub const CONSTRUCT_SCAN: &str = "clinvar/S_clinvar_cdna_scan.tsv";
pub const ALWAYS_REPORT: &str = "clinvar/S_clinvar_deep_intronic.tsv";

pub fn put_construct(d: &Sample) {
    d.put(CLINVAR, format!(
        "{VCF}chr7\t117548628\t1\tGTTT\tG\t50\tPASS\tGENEINFO=CFTR:1080;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\tGT:GQ:DP\t1/.:45:30\n\
        chr21\t25897626\t2\tTC\tGA\t50\tPASS\tGENEINFO=APP:351;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\tGT:GQ:DP\t0/1:50:35\n"
    ));
    d.put("clinvar/S_clinvar_cdna_check.tsv", "chrom\tpos\tref\talt\talt_reads\talt_clipped\tref_reads\tref_clipped\tcdna_check\tcdna_note\tconstruct_site\nchr7\t117548628\tGTTT\tG\t\t\t\t\tnot_checked\tnot a SNV, MNV or simple indel\t\nchr21\t25897626\tTC\tGA\t12\t0\t10\t0\tclear\t\tAPP Swedish K670N/M671L\n");
    d.put(CONSTRUCT_SCAN, "# cdna_scan: read 3 exon regions on 1 contigs; 1 loci listed\nchrom\tstart\tend\tjunctions\treads\tjunction_list\nchr21\t25880551\t25891785\t3\t24\t25880551-25881101:10,25881201-25891785:8,25880551-25891785:6\n");
    d.put(ALWAYS_REPORT, "chrom\tpos\tref\talt\tgene\thgvs\tclinvar_id\tcall\tdetail\nchr1\t215891198\tT\tC\tUSH2A\tNM_206933.4:c.7595-2144A>G\t30722\tref\tGQ 50 (gVCF block)\nchr7\t117639961\tC\tT\tCFTR\tNM_000492.4:c.3718-2477C>T\t7166\thet\tGQ 40, DP 30\nchrX\t101399747\tC\tT\tGLA\tNM_000169.3:c.640-801G>A\t10768\tno call\tno record in the VCF or the gVCF covers the site\n");
}

pub const SMN: &str = "paralogs/S_smn_copy_number.tsv";
pub const HBA: &str = "paralogs/S_hba_copy_number.tsv";
pub const MARKER: &str = "paralogs/S_smn1_marker.tsv";
pub const SMN_HEADER: &str = "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n";

pub fn put_parascopy(d: &Sample) {
    d.put(SMN, format!("{SMN_HEADER}chr5\t70925378\t70934928\tSMN1\tPASS\t4\t80\tPASS\t2,2\t40,40\tchr5:70049961-70059504:+\nchr5\t70950920\t70953320\tSMN1\tPASS\t4\t60\tPASS\t2,2\t33,30\tchr5:70075496-70077900:+\n"));
    d.put(HBA, "metric\tvalue\nstatus\tcalled\nalpha_globin_genes\t4\nestimate\t4.02\nhba2_depth\t30.1\nhba1_depth\t30.2\nflank_left_depth\t30.0\nflank_right_depth\t30.0\nmessage\tfour alpha-globin genes (aa/aa): no deletion seen\n");
    d.put(MARKER, "metric\tvalue\nmarker\tSMN1 g.27134T>G\nrsid\trs143838139\nstatus\tpresent\nalt_reads\t5\ndepth\t40\nfraction\t0.1250\ncopies\t0.5\nsmn1_pscn\t2 (quality 33, filter PASS)\nagcn\t4\nmessage\tThis marker does not prove 2+0, and its absence does not rule 2+0 out.\n");
}

pub const CLINICAL_VCF: &str = "clinical/S_clinical.vcf.gz";

pub fn put_clinical_mane(d: &Sample, zero_floor: bool) {
    let floor = if zero_floor {
        "##PGP_CALL_FLOOR=GQ=0,DP=0,HET_AB=0-1,HOM_AB=0\n"
    } else {
        ""
    };
    let head = VCF.replacen("#CHROM", &format!("{floor}##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Feature|CANONICAL|MANE_SELECT|MAX_AF\">\n#CHROM"), 1);
    d.gzip(CLINICAL_VCF, &format!(
        "{head}chr2\t100\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|TTN|ENST2|||0.0001,T|intron_variant|MODIFIER|TTN|ENST1|YES|NM_1|0.0001\tGT:GQ:DP:AD\t0/1:50:30:15,15\n\
        chr13\t200\t.\tA\tAT\t50\tPASS\tCSQ=AT|frameshift_variant|HIGH|BRCA2|ENST3|YES|NM_2|0.0001\tGT:GQ:DP:AD\t0/1:12:.:3,2\n\
        chr19\t300\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|LDLR|ENST4|YES|NM_3|0.0001\tGT:GQ:DP:AD\t1/1:11:17:12,5\n\
        chr11\t400\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|MYBPC3|ENST5|YES|NM_4|0.0001\tGT:GQ:DP:AD\t1/1:50:30:0,30\n\
        chr14\t500\t.\tCTTT\tC\t50\tPASS\tCSQ=-|frameshift_variant|HIGH|MYH7|ENST6|YES|NM_5|\tGT:GQ:DP:AD\t0/1:50:30:15,15\n\
        chr11\t600\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|KCNQ1|ENST7|YES|NM_6|0.3\tGT:GQ:DP:AD\t0/1:50:30:15,15\n\
        chr3\t700\t.\tC\tT\t50\tPASS\tCSQ=T|missense_variant|MODERATE|SCN5A|ENST8|||0.0001,T|stop_gained|HIGH|SCN5A|ENST9|YES||0.0001\tGT:GQ:DP:AD\t0/1:50:30:15,15\n"
    ));
    d.put("clinical/S_clinical_summary.tsv",
        "CHROM\tPOS\tREF\tALT\tGT\tIMPACT\tGENE\tConsequence\tMAX_AF\tCADD_PHRED\tREVEL\tAM_CLASS\tTRANSCRIPT\tWORST_IMPACT\tWORST_Consequence\tFREQ_TIER\tGQ\tDP\tAD\tCALL_CHECK\tCALL_WHY\n\
        chr13\t200\tA\tAT\t0/1\tHIGH\tBRCA2\tframeshift_variant\t0.0001\t.\t.\t.\tMANE\tHIGH\tframeshift_variant\trare\t12\t.\t3,2\tweak\tGQ12\n\
        chr14\t500\tCTTT\tC\t0/1\tHIGH\tMYH7\tframeshift_variant\t.\t.\t.\t.\tMANE\tHIGH\tframeshift_variant\tunknown\t50\t30\t15,15\tok\t.\n");
    d.put(SLIVAR,
        "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\tConsequence\tExisting_variation\tGT\tMAX_AF\tgnomADe_AF\tgnomADg_AF\tFREQ_TIER\tTRANSCRIPT\tMANE_Consequence\n\
        chr14\t500\tCTTT\tC\tHIGH\tMYH7\tframeshift_variant\t.\t0/1\t.\t.\t.\tunknown\tMANE\tframeshift_variant\n\
        chr13\t200\tA\tAT\tHIGH\tBRCA2\tframeshift_variant\t.\t0/1\t0.0001\t0.0001\t.\trare\tMANE\tframeshift_variant\n");
}

pub fn put_hla_qc(d: &Sample, nextflow: bool) {
    let dir = if nextflow { "hla" } else { "hla_t1k" };
    d.put(
        &format!("{dir}/S_hla_genotype.tsv"),
        "HLA-A\t2\tHLA-A*01:01:01\t30\t60\tHLA-A*26:01:01\t20\t55\n\
         HLA-B\t2\tHLA-B*57:01:01\t30\t60\tHLA-B*08:01:01\t20\t40\n\
         HLA-C\t2\tHLA-C*07:01:01\t5\t0\tHLA-C*07:02:01\t5\t30\n",
    );
    if nextflow {
        d.put(
            "hla_optitype/S_optitype_result.tsv",
            "# OptiType not run: no read pairs at HLA-A, HLA-B or HLA-C in this BAM\n",
        );
        d.put(
            "kir/S_kir_genotype.tsv",
            "# KIR not typed: T1K found too few reads at the KIR genes in this BAM\n",
        );
    } else {
        d.put("hla_optitype/S_hla_concordance.tsv", "Gene\tT1K\tOptiType\tAgree\tNote\n\
            HLA-A\t*01:01/*26:01\t*26:01/*01:01\tyes\t\n\
            HLA-B\t*57:01/*08:01\t*58:01/*08:01\tno\tthe callers disagree at two fields\n\
            HLA-C\t*07:01/*07:02\t*07:01/*07:02\tnot compared\tT1K gives an allele quality 0 or below\n");
        d.put(
            "kir_t1k/S_kir_genotype.tsv",
            "KIR2DL1\t2\tKIR2DL1*00302\t40\t60\tKIR2DL1*0010101\t30\t35\n\
            KIR3DL1\t1\tKIR3DL1*0150201\t50\t50\t.\t0\t-1\n\
            KIR2DS1\t0\t.\t0\t-1\t.\t0\t-1\n",
        );
    }
}

pub const DEPTH_CHECK: &str = "pypgx/S_cyp2d6_depth_check.tsv";
pub fn put_pgx_evidence(d: &Sample, depth: &str) {
    d.put(PYP, "Gene\tDiplotype\tPhenotype\tCNV_call\tSource\tStatus\tDefining_variants\tMin_GQ\tWeakest_variant\tDefining_filter\tQuality\n\
        CYP2D6\t*1/*2\tNormal Metabolizer\tNormal\tBAM\tcalled\t1/1 in the VCF\t50\tchr22:42126611 C>G GQ 50 AD 20,22\tPASS\tok\n\
        CYP2A6\t*1/*1\tno pypgx phenotype rule\tNormal\tBAM\tcalled\tnone\tN/A\tN/A\tnone\tN/A\n\
        GSTM1\t*0/*0\tno pypgx phenotype rule\tWholeDel2\tBAM\tcalled\tnone\tN/A\tN/A\tnone\tN/A\n\
        GSTT1\tN/A\tN/A\tN/A\tBAM\tnot callable on this reference\tN/A\tN/A\tN/A\tN/A\tN/A\n\
        CYP2C19\t*1/*2\tIntermediate Metabolizer\tN/A\tVCF\tcalled\t2/2 in the VCF\t12\tchr10:94781859 G>A GQ 12 AD 26,3\tPASS\tlow-confidence: check the reads\n");
    d.put(DEPTH_CHECK, format!("metric\tvalue\nstatus\t{depth}\nmessage\tCYP2D6 depth check\ngene_depth_all\t18.53\ngene_mapq1_share\t0.91\ngene_mapq20_share\t0.84\nbody_depth_all\t30.10\nbody_mapq1_share\t0.99\nbody_mapq20_share\t0.97\n"));
}

pub fn put_paraphase(d: &Sample, smn1: u8) {
    d.put(
        "paralogs/S_paraphase_copy_number.tsv",
        "region\tgene\tcopy_number\tstatus\tdetail\ttotal_cn\tregion_depth\tgenome_depth\n\
        smn1\tSMN1\t2\tcalled\t.\t4\t120\t31\n\
        smn1\tSMN2\t.\tno call\tone haplotype seen\t4\t120\t31\n\
        hba\tHBA1\t2\tcalled\tgenotype -a/aa (3 alpha-globin genes)\t3\t90\t31\n\
        hba\tHBA2\t1\tcalled\t4.2 kb deletion (-a4.2): HBA2 lost\t3\t90\t31\n\
        TNXB\tTNXB\t4\tcalled\tcopies of the whole gene family\t4\t60\t31\n\
        SSX2\tSSX2,SSX2B\t.\tno call\t.\t.\t0\t31\n",
    );
    d.put("paralogs/S_paraphase_variants.tsv", "region\tchrom\tpos\tref\talt\tfilter\thaplotypes\tsupport\nhba\tchr16\t172900\tC\tT\tPASS\thba_hba2hap1\thba_hba2hap1:20/21\n");
    d.put(SMN, format!("{SMN_HEADER}chr5\t70925000\t70953000\tSMN1\tPASS\t{}\t40\tPASS\t{smn1},2\t33,30\tchr5:70049000-70077000:-\n", smn1 + 2));
}

pub fn put_repeat_sizing(d: &Sample) {
    use std::fmt::Write;

    let spanning = "SPANNING/SPANNING";
    let genotype = "1/2";
    let mut rows = VCF.to_owned();
    for (locus, counts, source, limits, reads, status, normal, pathogenic) in [
        (
            "HTT",
            (17, 30),
            spanning,
            "17-17/30-30",
            "12/9",
            "normal",
            26,
            40,
        ),
        (
            "C9ORF72",
            (2, 8),
            spanning,
            "2-2/8-8",
            "10/11",
            "normal",
            24,
            30,
        ),
        (
            "ATXN1",
            (29, 31),
            spanning,
            "29-29/31-31",
            "8/7",
            "normal",
            35,
            39,
        ),
        (
            "DMPK",
            (5, 40),
            spanning,
            "5-5/40-40",
            "15/6",
            "pre_mutation",
            34,
            50,
        ),
        (
            "LOCUS_FL",
            (15, 60),
            "FLANKING/FLANKING",
            "0-119/4-507",
            "0/0",
            "normal",
            20,
            40,
        ),
        (
            "LOCUS_SP",
            (15, 60),
            spanning,
            "15-15/60-60",
            "5/5",
            "normal",
            100,
            200,
        ),
        (
            "LOCUS_LOW",
            (12, 18),
            spanning,
            "12-12/16-45",
            "6/1",
            "normal",
            20,
            40,
        ),
        (
            "RFC1",
            (100, 150),
            "INREPEAT/INREPEAT",
            "80-120/130-170",
            "0/0",
            "full_mutation",
            11,
            12,
        ),
    ] {
        let (a, b) = counts;
        let no_reads = "0/0";
        writeln!(rows, "chr1\t{}\t.\tC\t<STR50>\t.\tPASS\tEND=200;REPID={locus};VARID={locus};STR_STATUS={status};STR_NORMAL_MAX={normal};STR_PATHOLOGIC_MIN={pathogenic}\tGT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC\t{genotype}:{source}:{a}/{b}:{limits}:{reads}:{no_reads}:{no_reads}:30.0", 100 + locus.len()).unwrap();
    }
    d.put(STR, &rows);
    d.old(STR);
    d.put(STRANGER, rows);
    d.put("expansion_hunter/S_rfc1_motif.tsv", "locus\tchrom\tstart\tend\trepeat_reads\tdominant_motif\tdominant_class\tdominant_share\tmotif_shares\tverdict\tnote\nRFC1\tchr4\t39348425\t39348479\t24\tAAAAG\tbenign\t0.92\tAAAAG=0.92,AAAGG=0.08\tbenign_motif\tmotif AAAAG (benign) in 22 of 24 repeat reads (92%): the flag is not supported\n");
}

pub fn put_extra_repeats(d: &Sample, five_t: bool) {
    use std::fmt::Write;

    let genotype = "0/1";
    let mut rows = VCF.to_owned();
    let (tg, t, ta) = if five_t {
        ((11, 12), (7, 5), (7, 8))
    } else {
        ((11, 11), (7, 7), (7, 7))
    };
    for (chrom, pos, reference, end, locus, id, (a, b)) in [
        ("chr7", 117548606, 'T', 117548628, "CFTR", "CFTR_TG", tg),
        ("chr7", 117548628, 'T', 117548635, "CFTR", "CFTR_T", t),
        ("chr2", 233760234, 'A', 233760248, "UGT1A1", "UGT1A1", ta),
    ] {
        writeln!(rows, "{chrom}\t{pos}\t.\t{reference}\t<STR>\t.\tPASS\tEND={end};REPID={locus};VARID={id}\tGT:REPCN\t{genotype}:{a}/{b}").unwrap();
    }
    d.put(STR, &rows);
    d.old(STR);
    d.put(STRANGER, rows);
    d.put(STR_CATALOG, r#"[{"LocusId":"CFTR","NoThresholds":"the CFTR tract"},{"LocusId":"UGT1A1","NoThresholds":"the UGT1A1 repeat"}]"#);
    d.put(
        CPIC,
        "Gene\tDiplotype\tPhenotype\tStatus\nUGT1A1\t*1/*28\tsome function\tnon-normal\n",
    );
    let reference_call = "G/G";
    d.put("traits/S_traits.tsv", format!("rsid\tgene\tchrom\tpos\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids\nrs78655421\tCFTR\tchr7\t117530975\tCFTR R117H\t{reference_call}\tcalled\tr\tno\tn\t32404922\n"));
}
