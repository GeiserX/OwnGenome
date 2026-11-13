#!/usr/bin/env bash
# summary-collect-c.sh: writes the sample folders of the summary-collect
# parity cases of group c in the current directory, the same every time.
# Each folder is one case of og summary collect --only c against
# bin/collect_summary.py: every section of the group with a typical table,
# the Nextflow and other file names, files from an earlier run, unreadable
# tables of each kind, and the APOL1 G2 deletion written left-aligned,
# shifted, twice, under a reference block and not at all.
# py/ holds the modules the oracle imports, as they were at the oracle's
# commit, and the two bin/ tables it reads beside itself; assets/ is a copy
# of the repository's catalogs these sections read.
# scripts/ci/oracle-diff.sh runs it.
set -euo pipefail
export LC_ALL=C TZ=UTC
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "${here}/../../.." && pwd)
ORACLE_SHA=98aa9252ff77825a8a797d388b449931fa9a234e

T0=1577836800        # 2020-01-01 00:00:00 UTC, the start of the run
OLD=201906011230.00  # before the run
NEW=202103040506.07  # after it, the time of every file not set otherwise

at() { touch -t "$1" "${@:2}"; }   # [[CC]YY]MMDDhhmm[.SS], in UTC
put() { mkdir -p "$(dirname "$1")"; printf '%b' "$2" > "$1"; at "$NEW" "$1"; }
putgz() { mkdir -p "$(dirname "$1")"; printf '%b' "$2" | gzip -n -c > "$1"; at "$NEW" "$1"; }

# --- the oracle's modules -------------------------------------------------------
git -C "$root" cat-file -e "${ORACLE_SHA}^{commit}" 2> /dev/null \
  || git -C "$root" fetch -q --depth=1 origin "$ORACLE_SHA"
mkdir -p py
for m in paraphase_report ancestry_fine_report annotsv_postfilter cnvpytor_to_vcf family_report chip_report \
  traits_report genotype_rules_indels pgx_outside_calls; do
  git -C "$root" show "${ORACLE_SHA}:bin/${m}.py" > "py/${m}.py"
done
for t in hla_alleles.tsv str_intermediate_ranges.tsv; do
  git -C "$root" show "${ORACLE_SHA}:bin/${t}" > "py/${t}"
done

mkdir -p assets
for a in paraphase_genes.tsv paralog_loci.tsv mito_confirmed.tsv callability_genes.tsv; do
  cp "${root}/assets/${a}" assets/
done

VCF_HEAD='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'

# full: every section, the bash names -------------------------------------------
f=full
put $f/couple/S_couple_carrier.tsv "gene\tcondition\tresult\toffspring_risk\tpartner_a\tfindings_a\tpartner_b\tfindings_b\texplanation
CFTR\tCystic fibrosis\tshared\t1 in 4\tS\tc.1521_1523del (P)\tP2\tc.350G>A (LP)\tBoth carry a change.
DMD\tDuchenne\tx_linked\t1 in 4 sons\tP2\t.\tS\tc.10G>T (P)\tOn X.
GJB2\tDeafness\tone_partner\t.\tS\tc.35del (P)\tP2\t.\tOne side.
HBB\tBeta thalassaemia\tone_partner\t.\tP2\tc.20A>T\tS\t.\tOne side.\n"
put $f/paralogs/S_smn_copy_number.tsv "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions
chr5\t70950000\t70953000\tSMN1\tPASS\t4\t45.5\tPASS\t1,3\t31,25\tchr5:70075000-70078000
chr5\t70960000\t70970000\tSMN1\tPASS\t4\t40\tLowQual\t2,2\t10,12\t*\n"
put $f/paralogs/S_smn1_marker.tsv "metric\tvalue\nmarker\tSMN1 g.27134T>G\nstatus\tnot seen\nalt_reads\t0\ndepth\t31\nfraction\t0.000\ncopies\t0\nmessage\tno read carries the G\n"
put $f/paralogs/S_gba1_copy_number.tsv "chrom\tstart\tend\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual
chr1\t155231478\t155240000\tPASS\t4\t60\tPASS\t2,2\t30,30
chr1\t155240000\t155244699\tPASS\t4\t35.25\tPASS\t2,2\t30,30\n"
put $f/paralogs/S_pms2_copy_number.tsv "chrom\tstart\tend\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual
chr7\t5970000\t5980000\tPASS\t4\t50\tPASS\t2,2\t30,30
chr7\t5980000\t5990000\tLowQual\t3\t12\tPASS\t2,1\t30,30\n"
put $f/paralogs/S_mytable_copy_number.tsv "chrom\tstart\tend\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual
chr9\t1\t100\tnot assessed: no read over the locus\t.\t.\t.\t.\t.\n"
put $f/paralogs/S_strc_copy_number.tsv "chrom\tstart\tend\tagCN\n"
put $f/paralogs/S_hba_copy_number.tsv "metric\tvalue\nstatus\tcalled\ncopies\t4\nreading\ttwo copies of each gene\n"
put $f/paralogs/S_cn_intervals.tsv "locus\tname\tchrom\tstart\tend\tstatus\tcopies\texpected\tchange
DGCR\t22q11.2\tchr22\t18900000\t21500000\tcalled\t2\t2\tnone
AZFc\tAZFc\tchrY\t23000000\t26000000\tcalled\t0\t1\tdeletion
X1\tdup\tchr1\t1\t2\tcalled\t3\t2\tduplication
X2\tlow\tchr2\t1\t2\tnot assessed\t.\t2\t.\n"
put $f/paralogs/S_cyp21a2_psv.tsv "metric\tvalue\nstatus\tcalled\ncopies\t2\nmodule_total\t4\nfinding\tnone\n"
put $f/paralogs/S_cn_dosage.tsv "locus\tstatus\tcopies\testimate\treason\treading
HP\tcalled\t2\t2.1\t\tusual
LCE3BC\tcalled\t1\t1.0\t.\tone deleted
CCL3L1\tnot assessed\t\t\tlow depth\t\n"
put $f/paralogs/S_paraphase_copy_number.tsv "region\tgene\tcopy_number\tstatus\tdetail
pms2\tPMS2\t2\tcalled\t
smn1\tSMN1\t2\tcalled\tok
smn1\tSMN2\t2\tcalled\t
zz_other\tX\t2\tcalled\t.
zz_other2\tY\t.\tno call\t.
hba\tHBA1\t1\tcalled\t
hba\tHBA2\t2\tcalled\t\n"
put $f/paralogs/S_paraphase_variants.tsv "region\tpos\nsmn1\t1\nsmn1\t2\npms2\t3\n"
put $f/hla_t1k/S_hla_genotype.tsv "#comment\tline
HLA-A\t2\tA*02:01:01\t10\t60\tA*01:01:01\t9\t50
HLA-B\t2\tB*57:01:01:02\t10\t0\tB*08:01:01,B*08:01:02,B*08:01:03,B*08:01:04,B*08:01:05\t9\t0
HLA-C\t1\tC*07:01:01\t10\t40\t.\t0\t-1
HLA-DQA1\t2\tDQA1*05:01:01\t10\t40\tDQA1*02:01\t9\t30
HLA-DQB1\t2\tDQB1*02:01:01\t10\t40\tDQB1*02:02:01\t9\t30
HLA-DRB1\t2\tDRB1*03:01:01\t10\t40\tDRB1*07:01:01\t9\tnan
HFE\t1\tHFE*001\t3\t20\t.\t0\t-1
short\t1\n"
put $f/hla_optitype/S_hla_concordance.tsv "Gene\tT1K\tOptiType\tAgree\tNote\nHLA-A\t A*02:01 \tA*02:01\tyes\t\n\t\t\t\t\nHLA-B\tB*57:01\tB*57:01\tyes\tok\n"
put $f/kir_t1k/S_kir_genotype.tsv "KIR2DL1\t2\tKIR2DL1*001\t5\t30\tKIR2DL1*002\t4\t20\nKIR3DL3\t0\t.\t0\t-1\t.\t0\t-1\n"
put $f/prs/S_prs_summary.tsv "Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\tChrX\tAmbiguous\tApplies_To\tSex_Check\tPercentile\tAncestry_Group\tEvaluated_In\tEffect
CAD\tPGS000018\t0.5\t100\t120\t0\t12 dropped\tall\tok\t75\tEUR\tEUR;;AFR\tHR 1.5 (EUR)
T2D\tPGS000014\t0.1\t10\t20\t0\t0\tall\tok\tNA\tNA\t\t
Short\tPGS1\n"
put $f/ancestry/S_ancestry.tsv "key\tvalue\npopulation\tEUR\nreference_panel\tHGDP+1kGP\n\tskipped\npopulation_low_confidence\n"
put $f/ancestry_fine/S_ancestry_fine.tsv "level\tcode\tgroup\tfraction\tse\tlow\thigh
continent\tWEA\tWest Eurasia\t0.95\t0.01\t0.93\t0.97
continent\tAFR\tAfrica\t0.05\t0.01\tnan\tx
region\tIBE\tIberia\t0.60\t0.05\t0.5\t0.7
left_out\tOCE\tOceania\t.\t.\t.\t.
not_painted\tchr21\t.\t.\t.\t.\t.
other\tZ\tZ\t1\t1\t1\t1\n"
put $f/ancestry_fine/S_painting.png "\x89PNG\r\n\x1a\nnot a real picture\x00\xff"
putgz $f/manta/results/variants/diploidSV.vcf.gz "${VCF_HEAD}chr1\t100\tM1\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL\tGT\t0/1\nchr1\t200\tM2\tA\t<DUP>\t50\tMinQUAL\tSVTYPE=DUP\tGT\t0/1\nchr2\t5\n"
putgz $f/delly/S_sv.vcf.gz "${VCF_HEAD}chr1\t100\tD1\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL\tGT\t0/1\n"
put $f/cnvpytor/S_cnvs.txt "deletion\tchr2:10000001-10100000\t100000\t0.4\t1e-9\t1\t1\t1\t0.01\t0\t9
deletion\tchr1:122000001-123000000\t1000000\t0.1\t1e-12\t1\t1\t1\t0.2\t0\t9
duplication\tchr3:1000001-1200000\t200000\t1.8\t1e-6\t1\t1\t1\t0.8\t0\t9
deletion\tchr5:1-100000\t100000\t0.0\t1e-20\t1\t1\t1\t0\t0.9\t0
deletion\tchrY:2800001-2900000\t100000\t0.0\t1e-15\t1\t1\t1\t0.05\t0\t9
deletion\tchr4:1000001-1050000\t50000\t0.6\t0.2\t1\t1\t1\t0\t0\t9
duplication\tchr6:1-2\t1\t2\tnan\t1\t1\t1\t0\t0\t9
duplication\tweird\t1\t2\t0.001

\n"
putgz $f/sv_merged/S_sv_consensus.vcf.gz "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\ttiddit
chr1\t1000\tc1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;SUPP=2;SUPP_VEC=110\tGT\t0/1\t0/1\t./.
chr2\t100\tc2\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP:TANDEM;END=900;SUPP=3;SUPP_VEC=1111\tGT\t0/1\t0/1\t0/1
chr3\tx\tc3\tN\t<INV>\t.\tPASS\tSVTYPE=INV\tGT\t0/1\t0/1\t0/1\n"
put $f/annotsv/S_sv_class4_5_checked.tsv "AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tGene_name\tACMG_class\tFORMAT\tS\tAnnotation_mode\tCaller\tExon_edge\tExon_edge_detail\tPopulation_AF\tPopulation_source\tCommon\tAnnotSV_class\tRule\tRule_detail
a1\t1\t1100\t5100\t-4000\tDEL\tGENE1;GENE2;\t5\tGT:PS\t0/1:77\tsplit\tmanta\tyes\texon 2\t0.001\tgnomAD\tno\t5\t\t
a2\t2\t150\t880\t730\tDUP\tGENE3\t4\tGT\t1/1\tsplit\tdelly\tno\t\t\t\tno\t\t\t
a3\t7\t1\t100\t99\tINV\t\t3\tGT\t0/1\tsplit\tmanta\tnot checked\t\tx\tgnomAD\tyes\t4\tmhc\tthe MHC
a4\t9\t1\t100\t99\tDEL\tG9\t4\t\t\tfull\ttiddit\tno\t\t0.2\tgnomAD\tyes\t\t\t\n"
put $f/expansion_hunter/S_eh.vcf "${VCF_HEAD}chr4\t3074877\t.\tC\t<STR20>\t.\tPASS\tEND=3074940;REF=21;RL=63;RU=CAG;VARID=HTT;REPID=HTT\tGT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC\t0/1:SPANNING/SPANNING:17/20:17-17/20-20:10/12:5/5:0/0:30\n"
EH_RECS="chr4\t3074877\t.\tC\t<STR20>\t.\tPASS\tEND=3074940;REF=21;RL=63;RU=CAG;VARID=HTT;REPID=HTT;STR_STATUS=normal,normal;STR_NORMAL_MAX=26;STR_PATHOLOGIC_MIN=40\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:17/20:17-17/20-20:10/12
chrX\t147912051\t.\tC\t<STR45>\t.\tPASS\tEND=147912110;VARID=FMR1;REPID=FMR1;STR_STATUS=premutation;STR_NORMAL_MAX=44;STR_PATHOLOGIC_MIN=200\tGT:SO:REPCN:REPCI:ADSP\t1:SPANNING:58:55-60:1
chr4\t39348425\t.\tA\t<STR100>\t.\tPASS\tEND=39348480;VARID=RFC1;REPID=RFC1;STR_STATUS=full_mutation;STR_NORMAL_MAX=.;STR_PATHOLOGIC_MIN=400\tGT:SO:REPCN:REPCI:ADSP\t1/1:INREPEAT/INREPEAT:500/600:400-700/500-800:0/0
chr7\t117548628\t.\tT\t<STR11>\t.\tPASS\tEND=117548650;VARID=CFTR_TG;REPID=CFTR_TG\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:11/12:11-11/12-12:8/9
chr7\t117548651\t.\tT\t<STR5>\t.\tPASS\tEND=117548660;VARID=CFTR_T;REPID=CFTR_T\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:5/7:5-5/7-7:8/9
chr2\t233760233\t.\tA\t<STR8>\t.\tPASS\tEND=233760250;VARID=UGT1A1;REPID=UGT1A1;STR_STATUS=normal\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:7/8:7-7/8-8:5/6
chr6\t16327636\t.\tT\t<STR40>\t.\tPASS\tEND=16327700;VARID=ATXN1;REPID=ATXN1;STR_STATUS=normal;STR_NORMAL_MAX=35;STR_PATHOLOGIC_MIN=39\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/FLANKING:30/37:30-30/35-40:4/0
chr9\t27573528\t.\tC\t<STR8>\t.\tPASS\tEND=27573540;VARID=C9ORF72;REPID=C9ORF72\tGT:SO:REPCN:REPCI:ADSP\t./.:.:.:.:.
chr1\t10\t.\tC\t<STR8>\t.\tPASS\tEND=20;VARID=EXTRA1\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING:8/9:8-8/9-9:1
short\t1\t.\tC"
put $f/expansion_hunter/S_eh_stranger.vcf "${VCF_HEAD}${EH_RECS}\n"
at "$OLD" $f/expansion_hunter/S_eh.vcf
put $f/expansion_hunter/S_rfc1_motif.tsv "verdict\tnote\tdominant_motif\tdominant_class\tdominant_share\tmotif_shares\trepeat_reads\npathogenic\tAAGGG dominant\tAAGGG\tpathogenic\t0.9\tAAGGG:0.9\t40\n"
put $f/expansion_hunter/S_str_motifs.tsv "locus\tverdict\tnote\tdominant_motif\nFMR1\tinterrupted\tAGG\tCGG\n\tx\tx\tx\n"
put $f/expansion_hunter/S_eh_catalog.json '[{"LocusId": "EXTRA1", "NoThresholds": "no published range"}, {"LocusId": "X"}, 3]'
put $f/telomere/S/S/S_summary.tsv "PID\tsample\ttotal_reads\tread_lengths\tintratel_reads\ttel_content\nS\ttumor\t620000000\t100,150\t40\t1.5\n"
put $f/vcf/S_roh.txt "# RG\nRG\tS\tchrY\t2781480\t47781479\t45000000\t900\t50
RG\tS\tchrEBV\t1\t6000000\t6000000\t90\t50
RG\tS\tchr1\t125200001\t143200000\t18000000\t300\t50
RG\tS\tchr9\t42300001\t49300000\t7000000\t200\t50
RG\tS\tchr1\t120000001\t145000000\t25000000\t500\t50
RG\tS\tchr2\t10000001\t17000000\t7000000\t900\t50
RG\tS\tchr15\t23000001\t50000000\t27000000\t900\t50
RG\tS\tchr3\t1\t2\tx\t1\t1\n"
put $f/mito/S_haplogroup.txt "\"SampleID\"\t\"Haplogroup\"\t\"Rank\"\n\"S\"\t\"H5a1\"\t\"0.9\"\n"
put $f/mito/S_haplocheck.txt "\"Sample\"\t\"Contamination Status\"\t\"Contamination Level\"\n\"S\"\t\"NO\"\t\"\"\n"
put $f/y_haplogroup/S_y_haplogroup.txt "Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\nS\tR1b1a1b1a1a2*(xR-L21)\tR-M269\t100\t40\t0.98\n"
put $f/y_haplogroup/S_y_tree.txt "Hg\tLevel\tTree\tPath_derived\tSNPs_called\nR-Z2105\t12\tISOGG 2026\t30\t400\n"
put $f/y_haplogroup/S_y_tree_path.tsv "Level\tBranch\tSNPs\tDerived\n1\tA0-T\tM42\t1\n2\tR\n"
put $f/y_str/S_y_str.tsv "order\tpanel\tmarker\tcount\tconfidence\n1\tY-12\tDYS393\t13\thigh\n2\tY-12\tDYS390\t24\tlow\n3\tY-12\tDYS19\t\tno call\n4\tY-25\tDYS464\t.\tnot callable\n"
put $f/consumer_export/S_consumer_export.tsv "key\tvalue\nsites\t600000\ncalled\t590000\nno_call\t10000\nfrom_reference_blocks\t\ntemplate_sites\t631000\nno_call: low depth\t9000\nno_call: filtered\t1000\nfile_23andme\tS_23andme.txt\n"
put $f/bloodtype/S_bloodtype.tsv "system\tphenotype\tgenotype\tconfidence\nABO\tA\tA1/O1\thigh\nRhD\tRhD positive\t\tmedium\nKell\t\tK-/k+\t\n"
put $f/lpa/S_lpa_kiv2.tsv "metric\tvalue\nstatus\testimated\nkiv2_copies\t35\nci95_low\t30\nci95_high\t40\nclass\tlow\nmessage\tok\n"
put $f/lpa/S_lpa_kiv2_4925.tsv "metric\tvalue\nstatus\tpresent\nalt_reads\t4\ndepth\t200\n"
put $f/traits/S_traits.tsv "rsid\tgene\tchrom\tpos\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids
rs10455872\tLPA\tchr6\t1\tLp(a)\tA/G\tcalled\tone copy\tall\t.\t1
rs3798220\tLPA\tchr6\t2\tLp(a)\tT|C\tcalled\tone copy\tall\t.\t1
rs78655421\tCFTR\tchr7\t3\tR117H\tG/A\tcalled\tone copy\tall\t.\t1
rs429358\tAPOE\tchr19\t4\tAPOE\tT/C\tcalled\te3/e4\tall\t.\t1
rs1\tX\tchr1\t5\tx\t.\tnot_covered\tno reads\tall\t.\t1
rs2\tX\tchr1\t6\tx\tC/C\tother_allele\tother\tall\t.\t1\n"
put $f/traits/S_genotype_rules.tsv "rule\tgene\ttrait\trsids\tchrom\tpos\tgenotypes\tgq\tphase\tresult\tstatus\treading\tancestry\tnote\tpmids
APOE\tAPOE\tAPOE\trs429358,rs7412\tchr19\t1,2\tT/C,C/C\t40\t.\te3/e4\tcalled\tr\tall\t.\t1
TAS2R38\tTAS2R38\tbitter\trs1\tchr7\t1\tC/G\t40\tunphased\tPAV/AVI\tcalled\ttaster\tall\t.\t1
APOL1\tAPOL1\tkidney\trs73885319,rs60910145\tchr22\t1,2\tA/G,T/T\t40\t.\tG0/G1\tcalled\tone G1\tAfrican\t.\t20647424\n"
putgz $f/vcf/S.vcf.gz "${VCF_HEAD}chr22\t36265995\trs71785313\tAATAATT\tA\t50\tPASS\t.\tGT:GQ\t0/1:42\n"
put $f/cpic/S_phenotypes.tsv "Gene\tDiplotype\tPhenotype\tStatus\nUGT1A1\t*1/*28\tIntermediate\tnon-normal\n"
put $f/vntyper/S_vntyper.tsv "metric\tvalue\nstatus\tnegative\nconfidence\t\nvariant\t.\nflag\tnone\nadvice\tnone\n"
put $f/locityper/S_locityper.tsv "locus\tregion\tstatus\tgenotype\thaplotype_1\thaplotype_2\tquality\ttotal_reads\tunexplained_reads\tweighted_distance\twarnings\tfilter\treason
HLA-DRB1\tchr6:1-2\tcalled\th1,h2\th1\th2\t30\t100\t1\t0.1\t\tpass\t
KIR\tchr19:1-2\tcalled\tk1,k2\tk1\tk2\t20\t100\t1\t0.1\t\tpass\t
CYP2D6\tchr22:1-2\tno call\t.\t.\t.\t.\t.\t.\t.\tlow\tfail\tlow depth
FCGR\tchr1:1-2\tcalled\tf\tf\tf\t9\t1\t1\t1\t\tfail\tq
ZZZ\tchr1:1-2\tcalled\tz\tz\tz\t9\t1\t1\t1\t\tpass\t\n"
putgz $f/mito/S_chrM_filtered.vcf.gz "${VCF_HEAD}chrM\t302\t.\tA\tAC,ACC\t50\tPASS\t.\tGT:AF\t0/1:0.3,0.2
chrM\t3243\t.\tA\tG\t50\tPASS\t.\tGT:AF\t0/1:0.40
chrM\t8344\t.\tA\tG\t50\tweak_evidence\t.\tGT:AF\t0/1:0.30
chrM\t11778\t.\tG\tA\t50\tPASS\t.\tGT:AF\t1:0.99
chrM\t16189\t.\tT\tC\t50\tPASS\t.\tGT:AF\t0/1:0.5
chrM\t100\t.\tG\tA\t50\tPASS\t.\tGT:AF\t0/1:.
chrM\t150\t.\tC\tT\t50\tPASS\t.\tGT\t0/1
chrM\t200\t.\tC\tT\t50\tPASS\t.\tGT:AF\t0/1:nan,0.5\n"
put $f/mito/S_chrM_annotated.tsv "# clinvar\t2026-10-01\t-
# gnomad\tnot annotated\tno file
chrom\tpos\tref\talt\tfilter\taf\tclinvar_id\tclinvar_significance\tclinvar_review\tclinvar_condition\tgnomad_af_hom\tgnomad_af_het\tgnomad_max_hl
chrM\t3243\tA\tG\tPASS\t0.40\t9589\tPathogenic\tcriteria\tMELAS\t0\t0\t0
chrM\t302\tA\tAC\tPASS\t0.3\t.\t.\t.\t.\t.\t.\t.
chrM\t16189\tT\tC\tPASS\t0.99\t1\tBenign|not_provided\t.\t.\t.\t.\t.
chrM\t500\tC\tT\tPASS\t0.99\t2\tUncertain_significance\t.\t.\t.\t.\t.
chrM\t8344\tA\tG\tweak_evidence\t0.3\t9579\tPathogenic\t.\t.\t.\t.\t.\n"
put $f/mito/S_chrM_deletions.tsv "bkp5\tbkp3\tlength\tlevel_pct\n8469\t13447\t4977\t0.5\n100\t200\t.\t2.5\n300\t400\t99\t0.1\n"
put $f/mito/S_mtdna_copies.tsv "metric\tvalue\nstatus\testimated\nreason\t.\nchrm_reads\t100000\nmtdna_copies\t250\n"
put $f/archaic/S_archaic.tsv "key\tvalue\npopulation\tCEU\nneanderthal_segments\t120\nneanderthal_bp\t50000000\ndenisovan_segments\t3\ndenisovan_bp\t200000\nanalysed_bp\t2500000000\nneanderthal_percent\t2.00\ndenisovan_percent\t0.01\nlongest_10\tj\nlongest_2\tb\nlongest_1\ta\n"
put $f/ancient/S_ancient.tsv "rank\tgroup\tn\tdistance\tdates\tcountry\tlocality\n1\tSpain_BA\t10\t0.01\t2000 BCE\tSpain\tx\n2\tFrance_N\t5\t0.02\t\tFrance\t\n"
put $f/ancient/S_ancient_info.tsv "key\tvalue\nrelease\tv62\nsites_used\t500000\n"
put $f/ancient/S_ancient_pca.svg '<svg xmlns="http://www.w3.org/2000/svg"></svg>\n'
put $f/viral/S_viral.tsv "virus\taccession\treads\tbreadth\trelative_depth\tcopies_per_cell\tcall\thuman_depth\tnote
phiX\tNC_001422\t500\t1.0\t.\t.\tspike-in\t30\t
HHV-6B\tNC_000898\t2000\t0.99\t1.0\t1.0\tinherited integrated HHV-6 likely\t30\t
EBV\tNC_007605\t0\t0\t0\t0\tnot found\t30\t\n"
put $f/exomiser/S_exomiser_top.tsv "rank\tgene\tmoi\texomiser_score\tphenotype_score\tvariant\tgenotype\tconsequence\thgvs\tacmg\tdisease
1\tGENE1\tAD\t0.9\t0.8\tchr1-1-A-G\t0/1\tmissense\tp.X\tLP\tD1
1\tGENE1\tAD\t0.9\t0.8\tchr1-2-A-G\t0/1\tstop\tp.Y\tP\tD1
2\tGENE2\tAR\t0.5\t0.4\tchr2-1-A-G\t1/1\tsplice\tc.1\tVUS\tD2\n"
put $f/exomiser/S_exomiser_terms.tsv "hpo_id\tlabel\nHP:0001250\tSeizure\nHP:0000001\n"
put $f/chip/S_chip.tsv "gene\tprotein_change\tvaf\talt_reads\tdepth\tstatus\tnote\tchrom
DNMT3A\tR882H\t0.05\t3\t60\tchip\t\tchr2
TET2\tQ1\t0.5\t30\t60\tgermline\tlikely inherited\tchr4\n"
put $f/chip/S_chip_depth.tsv "gene\tmean_depth\tdetection_limit_vaf\tstatus\nDNMT3A\t30.5\t0.1\tcovered\nTET2\t40\t0.08\tcovered\nASXL1\t2\t.\tlow depth\nJAK2\t35\t0.1\tcovered\nTP53\t25\t0.1\tcovered\n"
put $f/mei/S_mei_sites.tsv "gene\telement\tcondition\tinheritance\tstatus\tmean_depth\tcall
FKTN\tSVA\tFukuyama\tAR\tcallable\t30\tnot seen
TAF1\tSVA\tXDP\tXL\tnot callable\t2\t\n"
put $f/mei/S_mei.tsv "chrom\tpos\telement\tcopies\nchr1\t1\tALU\t1\nchr2\t2\tL1\t1\nchr3\t3\tALU\t2\n"
put $f/mosaic/S_mocha_summary.tsv "metric\tvalue\nstatus\tassessed\ndepth\t30\ncalls_reported\t1\nfloor\t0.05\n"
put $f/mosaic/S_mocha.tsv "chrom\tstart\tend\tlength_mb\tarms\ttype\tcell_fraction\tn_hets\tlod_baf_phase\treported\treason
chr20\t1\t2\t10\tq\tLoss\t0.1\t50\t20\tyes\t
chr1\t1\t2\t0.5\tp\tCNN-LOH\t0.02\t5\t3\tno\tbelow floor\n"

# nf: the Nextflow and other names, a run with stale files and a female --------
f=nf
put $f/logs/run_status.tsv "meta\tstarted_epoch\t${T0}\nmeta\tdeclared_sex\tfemale\nstep\t08\tok\nstep\t18\tskipped (opt-in)\nstep\t10\tfailed (exit 1)\nstep\t11\tok\nstep\t12\tskipped\nstep\t22\tfailed (exit 1)\nstep\t25\tfailed (exit 2)\nstep\t37\tok\nstep\t53\tskipped\n"
put $f/hla/S_hla_genotype.tsv "HLA-DQA1\t2\tDQA1*03:01\t10\t0\tDQA1*01:01\t9\t30\nHLA-DQB1\t1\tDQB1*02:01\t10\t40\t.\t0\t-1\nHLA-DRB1\t2\tDRB1*04:01\t10\t40\tDRB1*07:01\t9\t40\n"
put $f/hla_optitype/S_optitype_result.tsv "# OptiType did not run: no class I reads\n"
put $f/kir/S_kir_genotype.tsv "# too few KIR reads\n"
put $f/cnvpytor/S_cnvs.txt "deletion\tchrY:2800001-2900000\t100000\t0.0\t1e-15\t1\t1\t1\t0.05\t0\t9\ndeletion\tchr2:10000001-10100000\t100000\t0.4\t1e-9\n"
at "$OLD" $f/cnvpytor/S_cnvs.txt
put $f/telomere/x/deep/T_summary.tsv "PID\tsample\tintratelomeric_reads\ttel_content\nS\ttumor\t61873\t334\n"
put $f/roh/S_roh.txt "RG\tS\tchr2\t10000001\t17000000\t7000000\t900\t50\n"
at "$OLD" $f/roh/S_roh.txt
put $f/mito/S_haplogroup.txt "SampleID\tRank\nS\tU5a\n"
at "$OLD" $f/mito/S_haplogroup.txt
put $f/y_haplogroup/S_y_haplogroup.txt "Sample_name\tHg\nS\tNA\n"
put $f/y_haplogroup/S_y_tree.txt "Hg\tLevel\tTree\n\t\t\n"
put $f/y_haplogroup/S_y_tree_path.tsv "Level\tBranch\tSNPs\n"
put $f/manta/x/S.diploidSV.vcf.gz ""
putgz $f/manta/x/diploidSV.vcf.gz "${VCF_HEAD}"
putgz $f/sv_merged/S_sv_consensus.vcf.gz "${VCF_HEAD}chr1\t1\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=10;SUPP=1;SUPP_VEC=1\tGT\t0/1\n"
at "$OLD" $f/sv_merged/S_sv_consensus.vcf.gz
put $f/annotsv/S_sv_class4_5_checked.tsv ""
put $f/prs/S_prs_summary.tsv "Condition\tPGS_ID\tPercentile\nX\tPGS1\tNA\n"
put $f/ancestry/S_ancestry.tsv "key\tvalue\n"
at "$OLD" $f/prs/S_prs_summary.tsv
put $f/expansion_hunter/S_eh_stranger.vcf "${VCF_HEAD}chr4\t3074877\t.\tC\t<STR20>\t.\tPASS\tVARID=HTT;REPID=HTT;STR_STATUS=normal\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:17/20:17-17/20-20:10/12\n"
at "$OLD" $f/expansion_hunter/S_eh_stranger.vcf
put $f/expansion_hunter/S_eh.vcf "${VCF_HEAD}chr7\t117548628\t.\tT\t<STR11>\t.\tPASS\tVARID=CFTR_TG;REPID=CFTR_TG\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:12/12:12-12/12-12:8/9
chr7\t117548651\t.\tT\t<STR5>\t.\tPASS\tVARID=CFTR_T;REPID=CFTR_T\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:5/5:5-5/5-5:8/9
chr2\t233760233\t.\tA\t<STR8>\t.\tPASS\tVARID=UGT1A1;REPID=UGT1A1\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:8/10:8-8/10-10:5/6\n"
put $f/locityper/S_locityper.tsv "locus\tregion\tstatus\tgenotype\thaplotype_1\thaplotype_2\tquality\ttotal_reads\tunexplained_reads\tweighted_distance\twarnings\tfilter\treason
HLA-DRB1\tx\tcalled\th\th\th\t1\t1\t1\t1\t\tpass\t
HLA-DQB1\tx\tcalled\th\th\th\t1\t1\t1\t1\t\tpass\t
CYP2D6\tx\tcalled\th\th\th\t1\t1\t1\t1\t\tpass\t\n"
at "$OLD" $f/locityper/S_locityper.tsv
put $f/lpa/S_lpa_kiv2.tsv "metric\tvalue\nstatus\ttoo little depth\n"
put $f/traits/S_traits.tsv "rsid\tgene\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids
rs10455872\tLPA\tLp(a)\t.\tnot_covered\tno reads\tall\t.\t1\n"
putgz $f/mito/S_chrM_filtered.vcf.gz "${VCF_HEAD}chrM\t3243\t.\tA\tG\t50\tPASS\t.\tGT:AF\t0/1:0.40\n"
put $f/mito/S_chrM_annotated.tsv "chrom\tpos\tref\talt\tfilter\taf\n"
at "$OLD" $f/mito/S_chrM_annotated.tsv
put $f/mosaic/S_mocha_summary.tsv "metric\tvalue\nstatus\tnot assessed\nmessage\ttoo shallow\ncalls_reported\t.\n"
put $f/mosaic/S_mocha.tsv "chrom\tstart\tend\tlength_mb\tarms\ttype\tcell_fraction\tn_hets\tlod_baf_phase\treported\treason\n"
put $f/vntyper/S_vntyper.tsv "metric\tvalue\nstatus\tnot assessed\n"
put $f/exomiser/S_exomiser_top.tsv "rank\tgene\tmoi\texomiser_score\tphenotype_score\tvariant\tgenotype\tconsequence\thgvs\tacmg\tdisease\n"
put $f/chip/S_chip.tsv "gene\tprotein_change\tvaf\talt_reads\tdepth\tstatus\tnote\n"
put $f/chip/S_chip_depth.tsv "gene\tmean_depth\tdetection_limit_vaf\tstatus\nDNMT3A\t2\t.\tlow depth\n"
put $f/viral/S_viral.tsv "virus\taccession\treads\tbreadth\trelative_depth\tcopies_per_cell\tcall\thuman_depth\tnote\nphiX\tNC\t0\t0\t.\t.\tnot found\t\t\n"
put $f/mei/S_mei_sites.tsv "gene\telement\tcondition\tinheritance\tstatus\tmean_depth\tcall\nFKTN\tSVA\tF\tAR\tcallable\t30\tseen\n"
put $f/mei/S_mei.tsv "chrom\tpos\telement\tcopies\n"
put $f/consumer_export/S_consumer_export.tsv "key\tvalue\nsites\t10\ncalled\t10\nno_call\t0\n"
put $f/paralogs/S_hba_copy_number.tsv "metric\tvalue\nnote\tno status\n"
put $f/paralogs/S_paraphase_copy_number.tsv "region\tgene\tcopy_number\tstatus\nsmn1\tSMN1\tx\tcalled\n"
put $f/paralogs/S_smn_copy_number.tsv "chrom\tstart\tend\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual
chr5\t1\t2\tnot assessed: no model\t.\t.\t.\t.\t.\n"

# bad1, bad2: an unreadable table of each kind ------------------------------------
f=bad1
put $f/couple/S_couple_carrier.tsv "gene\tcondition\tresult\toffspring_risk\tpartner_a\tfindings_a\tpartner_b\tfindings_b\texplanation\nG\tc\tshared\t1\tA\tx\tB\ty\tz\n"
put $f/ancestry_fine/S_ancestry_fine.tsv "level\tcode\n"
put $f/hla_t1k/S_hla_genotype.tsv "HLA-A\t2\tA*02:01\t1\tx\n"
put $f/mito/S_haplogroup.txt ""
put $f/mito/S_haplocheck.txt "Sample\tStatus\nS\tNO\n"
put $f/y_haplogroup/S_y_haplogroup.txt "Sample\tGroup\nS\tR\n"
put $f/y_haplogroup/S_y_tree.txt "Hg\tLevel\nR\t1\n"
put $f/y_str/S_y_str.tsv "order\tpanel\tmarker\tcount\tconfidence\n1\tY\tD\t1\tmaybe\n"
put $f/consumer_export/S_consumer_export.tsv "key\tvalue\nsites\tmany\n"
put $f/bloodtype/S_bloodtype.tsv "system\tphenotype\n"
put $f/archaic/S_archaic.tsv "key\tvalue\nneanderthal_segments\t1\nneanderthal_bp\t1\ndenisovan_segments\t1\ndenisovan_bp\t1\nanalysed_bp\t1\nneanderthal_percent\t2\n"
put $f/lpa/S_lpa_kiv2.tsv "metric\tvalue\nstatus\tmaybe\n"
put $f/vntyper/S_vntyper.tsv "metric\tvalue\n"
put $f/locityper/S_locityper.tsv "locus\tstatus\nX\tcalled\n"
put $f/viral/S_viral.tsv "virus\treads\nX\t1\n"
put $f/chip/S_chip.tsv "gene\tvaf\n"
put $f/chip/S_chip_depth.tsv "gene\tmean_depth\tdetection_limit_vaf\tstatus\nX\t1\t1\tcovered\n"
put $f/mei/S_mei_sites.tsv "gene\telement\tcondition\tinheritance\tstatus\tmean_depth\tcall\nX\tSVA\tc\tAR\tcallable\t1\tseen\n"
put $f/mosaic/S_mocha_summary.tsv "metric\tvalue\nstatus\tassessed\ncalls_reported\t2\n"
put $f/mosaic/S_mocha.tsv "chrom\tstart\tend\tlength_mb\tarms\ttype\tcell_fraction\tn_hets\tlod_baf_phase\treported\treason\nchr1\t1\t2\t1\tp\tLoss\t0.1\t1\t1\tyes\t\n"
put $f/exomiser/S_exomiser_top.tsv "rank\tgene\n"
put $f/traits/S_traits.tsv "rsid\tgene\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids\nrs1\tX\tx\tA/A\tcalled\tr\n"
put $f/ancient/S_ancient.tsv "rank\tgroup\n1\tx\n"
put $f/annotsv/S_sv_class4_5_checked.tsv "SV_chrom\tSV_start\n1\t2\n"
put $f/paralogs/S_smn_copy_number.tsv "chrom\tstart\tend\n"
put $f/paralogs/S_paraphase_copy_number.tsv "region\tgene\nsmn1\tSMN1\n"
put $f/paralogs/S_cn_dosage.tsv "locus\tstatus\n"
put $f/telomere/S/S_summary.tsv "PID\tsample\tintratel_reads\tread_lengths\ttel_content\nS\ttumor\t40\t150\t1.5\n"
put $f/expansion_hunter/S_eh.vcf "${VCF_HEAD}chr7\t1\t.\tT\t<STR11>\t.\tPASS\tVARID=CFTR_TG\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/FLANKING:11/12:11-11/12-12:8/0\nchr7\t2\t.\tT\t<STR5>\t.\tPASS\tVARID=CFTR_T\tGT:SO:REPCN:REPCI:ADSP\t0/1:SPANNING/SPANNING:7/9:7-7/9-9:8/9\n"
put $f/mito/S_chrM_deletions.tsv "bkp5\tbkp3\n1\t2\n"
putgz $f/mito/S_chrM_filtered.vcf.gz "${VCF_HEAD}"
put $f/vcf/S_roh.txt "RG\tS\tchr1\t1\t2\tx\n"

f=bad2
put $f/couple/S_couple_carrier.tsv "gene\tcondition\n"
put $f/mei/S_mei_sites.tsv "gene\telement\tcondition\tinheritance\tstatus\tmean_depth\tcall\nX\tSVA\tc\tAR\tcallable\t1\tseen\n"
put $f/mei/S_mei.tsv "chrom\tpos\n"
put $f/mosaic/S_mocha_summary.tsv "metric\tvalue\nstatus\tassessed\n"
put $f/traits/S_traits.tsv "rsid\tgene\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids\nrs1\tX\tx\tA/A\tcalled\tr\ta\tn\tp\n"
put $f/traits/S_genotype_rules.tsv "rule\tgene\nAPOL1\tAPOL1\n"
put $f/archaic/S_archaic.tsv "key\tvalue\nneanderthal_segments\tmany\n"
put $f/viral/S_viral.tsv "virus\taccession\treads\tbreadth\trelative_depth\tcopies_per_cell\tcall\nX\tA\t1\t1\t1\t1\tc\n"
put $f/consumer_export/S_consumer_export.tsv "key\tvalue\nsites\t1\ncalled\t1\nno_call\t0\nno_call: other\tsome\n"
put $f/exomiser/S_exomiser_top.tsv "rank\tgene\tmoi\texomiser_score\tphenotype_score\tvariant\tgenotype\tconsequence\thgvs\tacmg\tdisease\n1\tG\n"
put $f/y_haplogroup/S_y_tree.txt "Hg\tLevel\tTree\nR\t1\tT\n"
put $f/y_haplogroup/S_y_tree_path.tsv "Level\tSNPs\n1\tM\n"
put $f/paralogs/S_cn_intervals.tsv "locus\tname\tchrom\tstart\tend\tstatus\tcopies\texpected\tchange\nX\tx\tchr1\t1\t2\tcalled\n"
put $f/mito/S_chrM_filtered.vcf.gz "not gzip\n"
put $f/delly/S_sv.vcf.gz "\x1f\x8b\x08\x00broken"
put $f/expansion_hunter/S_eh_stranger.vcf "${VCF_HEAD}chr4\t39348425\t.\tA\t<STR100>\t.\tPASS\tVARID=RFC1;REPID=RFC1;STR_STATUS=full_mutation\tGT:SO:REPCN:REPCI:ADSP\t1/1:INREPEAT/INREPEAT:500/600:400-700/500-800:0/0\n"

# g2-*: the APOL1 G2 deletion beside step 42's APOL1 line --------------------------
g2() {
  local d=$1 result=$2 status=$3
  put "$d/traits/S_traits.tsv" "rsid\tgene\ttrait\tgenotype\tstatus\treading\tancestry\tnote\tpmids\nrs1\tX\tx\tA/A\tcalled\tr\ta\tn\tp\n"
  put "$d/traits/S_genotype_rules.tsv" "rule\tgene\ttrait\trsids\tchrom\tpos\tgenotypes\tgq\tphase\tresult\tstatus\treading\tancestry\tnote\tpmids\nAPOL1\tAPOL1\tkidney\trs1\tchr22\t1\tA/G\t40\t.\t${result}\t${status}\tr\tAfrican\t.\t1\n"
}
g2 g2-left G0/G1 called
putgz g2-left/vcf/S.vcf.gz "${VCF_HEAD}chr21\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\nchr22\t36265995\trs71785313\tAATAATT\tA\t50\tPASS\t.\tGT:GQ\t0/1:42\nchr22\t36266100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
# The same deletion written one base to the right: the same sequence, so G2.
g2 g2-shifted G0/G0 called
putgz g2-shifted/vcf/S.vcf.gz "${VCF_HEAD}chr22\t36265996\t.\tATAATTA\tA\t50\tPASS\t.\tGT:GQ\t1/1:30\n"
g2 g2-two G1/G1 called
putgz g2-two/vcf/S.vcf.gz "${VCF_HEAD}chr22\t36265995\t.\tAATAATT\tA\t50\tPASS\t.\tGT:GQ\t0/1:42\nchr22\t36265996\t.\tATAATTA\tA\t50\tPASS\t.\tGT:GQ\t0/1:40\n"
g2 g2-gvcf . .
putgz g2-gvcf/vcf/S.g.vcf.gz "${VCF_HEAD}chr22\t36265990\t.\tA\t<NON_REF>\t.\t.\tEND=36266010\tGT:GQ\t0/0:50\n"
putgz g2-gvcf/vcf/S.vcf.gz "${VCF_HEAD}chr22\t36265995\t.\tAATAATT\tA\t50\tPASS\t.\tGT:GQ\t1/1:42\n"
g2 g2-none G0/G0 called
putgz g2-none/vcf/S.vcf.gz "${VCF_HEAD}chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
g2 g2-other G1/G1 called
putgz g2-other/vcf/S.vcf.gz "${VCF_HEAD}chr22\t36265998\t.\tATTATA\tA\t50\tPASS\t.\tGT:GQ\t0/1:42\nchr22\t36266000\t.\tT\tC\t50\tLowQual\t.\tGT\t0/1\n"
g2 g2-nocall G0/G2 called
putgz g2-nocall/vcf/S.g.vcf.gz "${VCF_HEAD}chr22\t36265995\t.\tAATAATT\tA,<NON_REF>\t50\tPASS\t.\tGT:GQ\t./.:0\n"
g2 g2-novcf G1/G1 called
