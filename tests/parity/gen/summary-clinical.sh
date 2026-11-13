#!/usr/bin/env bash
# summary-clinical.sh: writes the sample folders of the summary-collect
# parity cases of group b (og summary collect --only b) in the current
# directory, the same every time. The clinical sections: ClinVar hits with
# each review status, a call floor, combined alleles, half calls and
# phased calls, the cDNA check, scan and deep intronic tables, gnomAD AF
# from an indexed VEP output; the carrier list; PharmCAT, CPIC, pypgx and
# Cyrius and step 36's CYP2D6 row; CPSR with its own class beside
# ClinVar's; the clinical filter's table and VCF with the ACMG SF tier;
# slivar's table and VCF with the phased compound-het pairs; and the
# deletions the two-hit rule pairs with a variant (AnnotSV, step 35's
# interval, HBA and SMN1 tables), in trans, in cis and inside the deletion.
#
# The VEP output is the committed pair cb_vep.vcf.gz and its .tbi, staged
# beside the folders; vep/S_vep.vcf.gz links to it. assets/ is a copy of
# the repository's gene lists. pymods/ holds the helper modules the Python
# collector imports, at a pinned commit, so the oracle reads them instead
# of stubs. scripts/ci/oracle-diff.sh runs it.
set -euo pipefail
export LC_ALL=C TZ=UTC
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "${here}/../../.." && pwd)

put() { mkdir -p "$(dirname "$1")"; printf '%b' "$2" > "$1"; }
putgz() { mkdir -p "$(dirname "$1")"; printf '%b' "$2" | gzip -n -c > "$1"; }
at() { touch -t "$1" "${@:2}"; }   # [[CC]YY]MMDDhhmm[.SS], in UTC

T0=1577836800        # 2020-01-01 00:00:00 UTC, the start of the run
OLD=201906011230.00  # before the run
NEW=202103040506.07  # after it
NEWER=202103040507.00

# The helper modules of bin/collect_summary.py, as they were when the
# collector was last the summary's writer.
HELPERS=7c1470bc50933f48a1a0dde3c0e0b2b87bc48718
git -C "$root" cat-file -e "${HELPERS}^{commit}" 2> /dev/null \
  || git -C "$root" fetch -q --depth=1 origin "$HELPERS"
mkdir -p pymods
for m in pgx_outside_calls clip_check annotsv_postfilter paraphase_report traits_report; do
  git -C "$root" show "${HELPERS}:bin/${m}.py" > "pymods/${m}.py"
done

mkdir -p assets
for a in acmg_sf.tsv carrier_genes.tsv paralog_loci.tsv; do
  cp "${root}/assets/${a}" assets/
done

VCF_HEAD='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
PL='criteria_provided,_multiple_submitters,_no_conflicts'

# --- the ClinVar hits of step 06 ---------------------------------------------------
HITS="##fileformat=VCFv4.2\n##PGP_CALL_FLOOR=GQ=20,DP=8\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
HITS+="chr1\t45331556\t1\tC\tT\t.\tPASS\tGENEINFO=MUTYH:4595;CLNSIG=Pathogenic;CLNREVSTAT=${PL}\tGT:GQ:DP:PS\t0|1:50:30:45331000\n"
HITS+="chr6\t26092913\t2\tG\tA\t.\tPASS\tGENEINFO=HFE:3077;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\tGT:GQ:DP\t0/1:99:40\n"
HITS+="chr7\t117559590\t3\tATCT\tA\t.\tPASS\tGENEINFO=CFTR:1080;CLNSIG=Pathogenic;CLNREVSTAT=practice_guideline\tGT:GQ:DP\t0/1:15:30\n"
HITS+="chr13\t20189473\t4\tC\tT\t.\tPASS\tGENEINFO=GJB2:2706;CLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;PGP_JOINED_FROM=chr13:20189473:C:T|chr13:20189475:G:A\tGT:GQ:DP\t1/.:40:5\n"
HITS+="chr3\t15645186\t5\tG\tC\t.\tPASS\tGENEINFO=BTD:686;CLNSIG=Pathogenic/Likely_pathogenic;CLNREVSTAT=criteria_provided,_conflicting_classifications\tGT\t0/1\n"
HITS+="chr17\t43057100\t6\tC\tT\t.\tPASS\tGENEINFO=BRCA1:672|NBR2:10230;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;PGP_PART_OF=chr17:43057100:C:T|chr17:43057102:A:G\tGT:GQ:DP\t0/1:60:35\n"
HITS+="chrX\t101398000\t7\tC\tT\t.\tPASS\tGENEINFO=GLA:2717;CLNSIG=Pathogenic;CLNREVSTAT=no_assertion_criteria_provided;PGP_ONE_ALLELE=chrX:101398004:G:A\tGT\t1\n"
HITS+="chr16\t176700\t8\tG\tA\t.\tPASS\tGENEINFO=HBA1:3039;CLNSIG=Pathogenic;CLNREVSTAT=${PL}\tGT:GQ:DP\t0/1:45:31\n"
HITS+="chr5\t70951946\t9\tC\tT\t.\tPASS\tGENEINFO=SMN1:6606;CLNSIG=Pathogenic;CLNREVSTAT=${PL}\tGT:GQ:DP\t1/1:45:31\n"
HITS+="chr2\t100\t10\tA\tG\t.\tPASS\tCLNREVSTAT=;CLNSIG=\tGT:DP\t0/1:x\n"
HITS+="chr2\t200\t11\tA\tG\t.\tPASS\n"

CDNA='chrom\tpos\tref\talt\talt_reads\talt_clipped\tref_reads\tref_clipped\tcdna_check\tcdna_note\tconstruct_site\n'
CDNA+='chr1\t45331556\tC\tT\t12\t6\t10\t0\tflagged\tpossible cDNA/transgene contamination: check the reads\tyes\n'
CDNA+='chr6\t26092913\tG\tA\t20\t0\t20\t0\tclean\t\t\n'
CDNA+='chr6\tx\tG\tA\t1\t1\t1\t1\tclean\t\t\n'
SCAN='# cdna_scan: 2 genes with spliced reads\nchrom\tstart\tend\tjunctions\treads\tjunction_list\nchr1\t100\t200\t2\t5\tA:1-2,A:3-4\nchr2\tx\t1\t1\t1\t\n'
DEEP='chrom\tpos\tref\talt\tgene\thgvs\tclinvar_id\tcall\tdetail\nchr7\t117639961\tC\tT\tCFTR\tc.3718-2477C>T\t7129\t\t\nchr7\tNA\tC\tT\tCFTR\t\t\t\t\n'

# --- pharmacogenomics ----------------------------------------------------------------
CPIC='Gene\tDiplotype\tPhenotype\tStatus\nCYP2D6\t*1/*4\tIntermediate Metabolizer\tnon-normal\nCYP2C19\t*1/*17\tRapid\tambiguous\nCYP2C9\t*1/*1\tNormal\tnormal\nDPYD\n'
CYRIUS_HEADER='Sample\tGenotype\tFilter\n'
REC='CPIC lookup\n  WARNING:  CYP2C19 has two possible phenotypes\nnothing\nWARNING:last\n'
PYPGX='Gene\tDiplotype\tPhenotype\tStatus\tSource\tCNV_call\tQuality\tMin_GQ\tWeakest_variant\tDefining_filter\n'
PYPGX+='CYP2D6\t*1/*2\t Normal Metabolizer \tok\tBAM\tNormal\tok\t\t\tPASS\n'
PYPGX+='CYP2A6\tFAILED\t\tfailed\tBAM\t\t\t\t\t\n'
PYPGX+='CYP2B6\t*1/*6\tIM\tdepth check failed\tbam\tN/A\tok\t\t\t\n'
PYPGX+='CYP2E1\t.\t\tnot callable on this reference\tBAM\t.\t\t\t\t\n'
PYPGX+='UGT1A1\tIndeterminate\t\tok\tVCF\t\tlow-confidence: check the reads\t 12 \tchr2:1:A:G\t\n'
PYPGX+='SLCO1B1\t*1/*5\t\tok\tBAM\t.\tlow-confidence: check the reads\tx\t\t\n'
DEPTH='metric\tvalue\nstatus\tok\ngene_mapq1_share\t0.99\nbody_mapq20_share\t0.95\nother\t1\n'
COMP='Gene\tMatch\nCYP2D6\tYes\nCYP2C19\tNo\nCYP3A5\tpypgx only\nTPMT\tPharmCAT only\n'
CONS='Gene\tResult\tOutside_call\tReason\tEvidence\nHLA-A\tx\tno\t\t\nCYP2D6\t*1/*2\tyes\tpypgx and Cyrius agree\tpypgx *1/*2; Cyrius *2/*1\n'

# --- CPSR (step 17) --------------------------------------------------------------------
CPSR='VAR_ID\tGENOMIC_CHANGE\tSYMBOL\tCLASSIFICATION\tCPSR_CLASSIFICATION\tCLINVAR_CLASSIFICATION\tASSERTION_AUTHORITY\tACMG_CODE\tDP_CONTROL\n'
CPSR+='17_43057100_C_T\t17:g.43057100C>T\tBRCA1\tPathogenic\tLikely_Benign\tPathogenic\tClinVar\tPVS1\t6\n'
CPSR+='\t13:g.32338000G>A\tBRCA2\tLikely_Pathogenic\tVUS\t\tCPSR\t\t40\n'
CPSR+='2_47403000_A_G\t\tMSH2\tVUS\tLP\tB\t\t\t\n'
CPSR+='x\tbad\tATM\tBenign\tBenign\t\t\t\t3\n'
CPSR+='5_112839000_A_T\t\tAPC\t\tPathogenic\t\t\t\t30\n'
CPSR_CDNA="${CDNA}17\t43057100\tC\tT\t5\t5\t9\t0\tflagged\tsplice\t\n"

# --- the clinical filter (step 23) ------------------------------------------------------
CSQF='Allele|Consequence|IMPACT|SYMBOL|Feature|CANONICAL|MANE_SELECT|gnomADe_AF|gnomADg_AF'
CLIN="##fileformat=VCFv4.2\n##PGP_CALL_FLOOR=GQ=20,DP=8,HET_AB=0.25-0.75,HOM_AB=0.85\n"
CLIN+="##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: ${CSQF}\">\n"
CLIN+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
CLIN+='chr1\t45332000\t.\tG\tA\t50\tPASS\tCSQ=A|stop_gained|HIGH|MUTYH|T1|YES|NM_1|0.0001|\tGT:GQ:DP:AD:PS\t1|0:40:30:15,15:45331000\n'
CLIN+='chr2\t47414000\t.\tA\tAT\t50\tPASS\tCSQ=AT|frameshift_variant|HIGH|MSH2|T2|YES||0.00002|0.00001\tGT:GQ:DP:AD\t0/1:50:30:28,2\n'
CLIN+='chr10\t100\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|PTEN|T3|YES||0.2|\tGT\t0/1\n'
CLIN+='chr6\t26093000\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|HFE|T4|YES||&|\tGT:GQ:DP:AD\t1/1:12:5:3,4\n'
CLIN+='chr17\t43057100\t.\tC\tT\t50\tPASS\tCSQ=T|missense_variant|MODERATE|BRCA1|T5|YES|NM_2||,T|stop_gained|HIGH|BRCA1|T6|||0.0001|\tGT\t0/1\n'
CLIN+='chr11\t47332000\t.\tG\tA\t50\tPASS\tCSQ=A|splice_donor_variant&intron_variant|HIGH|MYBPC3|T7|||0.001|,A|stop_gained|HIGH|MYBPC3|T8|||0.0001|\tGT:GQ\t0/1:.\n'
CLINT='GENE\tIMPACT\tFREQ_TIER\tCALL_CHECK\nMUTYH\tHIGH\trare\tok\nMSH2\tHIGH\trare\tweak\nHFE\tHIGH\tunknown\tweak\n.\tMODERATE\trare\tok\n\t\t\t\n'

# --- slivar (step 31) -------------------------------------------------------------------
SLIV='CHROM\tPOS\tGENE\tFREQ_TIER\nchr1\t1\tA\trare\nchr1\t2\tA\tunknown\n'
CHETS='GENE\tCHROM\nA\tchr1\nA\tchr1\nB\tchr2\n\tchr3\n'
PAIRS='GENE\tCHROM\tPOS_1\tREF_1\tALT_1\tPOS_2\tREF_2\tALT_2\tPHASE\tTRANS_READS\tCIS_READS\n'
for i in $(seq 1 23); do PAIRS+="G${i}\tchr1\t${i}0\tA\tG\t${i}5\tC\tT\ttrans\t${i}\t1\n"; done
PAIRS+='H\tchr2\t1\tA\tG\t2\tC\tT\tcis\t0\t4\nK\tchr2\t3\tA\tG\nL\tchr2\t4\tA\tG\t5\tC\tT\tunknown\tx\t\n'

# --- the deletions of the two-hit rule --------------------------------------------------
ANN='SV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tGene_name\tACMG_class\tCaller\tExon_edge\tExon_edge_detail\tPopulation_AF\tPopulation_source\tCommon\tAnnotSV_class\tFORMAT\tS\tAnnotation_mode\n'
ANN+='1\t45338000\t45340000\t-2000\tDEL\tMUTYH\t5\tManta\tyes\t\t\t\tno\t5\tGT:PS\t1|0:45331000\tfull\n'
ANN+='chr2\t47386000\t47388000\t-2000\tDEL\tEPCAM\t4\tDelly\tno\t\t\t\tno\t4\tGT\t0/1\tfull\n'
ANN+='17\t43000000\t43100000\t100000\tDUP\tBRCA1;NBR2\t4\tCNVpytor\tno\t\t0.001\tgnomAD\tno\t4\tGT\t0/1\tfull\n'
ANN+='13\t20180000\t20200000\t-20000\tDEL\tGJB2;NOTAGENE\t5\tManta\tno\t\t\t\tno\t5\tGT\t0/0\tfull\n'
ANN+='13\t20180000\t20200000\t-20000\tDEL\tGJB2;XYZ1\t4\tManta\tno\t\t\t\tno\t4\tGT\t0/1\tfull\n'
ANN+='3\t15600000\t15700000\t-100000\tDEL\tBTD\t3\tManta\tno\t\t\t\tno\t5\tGT\t0/1\tfull\n'
ANN+='3\t15600000\t15700000\t-100000\tDEL\tBTD\t5\tManta\tno\t\t0.2\tgnomAD\tyes\t5\tGT\t0/1\tfull\n'
CONSV='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\tcnvpytor\n'
CONSV+='chr1\t45338010\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=45340020;SUPP=2;SUPP_VEC=110\n'
CONSV+='chr1\t45338500\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=45340900;SUPP=3;SUPP_VEC=111\n'
CONSV+='chr17\t43000000\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP:TANDEM;END=43100000;SUPP=1;SUPP_VEC=0011\n'
SMN='chrom\tstart\tend\tagCN\tagCN_filter\tagCN_qual\tpsCN\tpsCN_qual\tpsCN_filter\thomologous_regions\n'
SMN+='chr5\t70925000\t70953000\t4\tPASS\t40\t1,3\t30,30\tPASS\tchr5:70049000-70078000\n'
HBA='metric\tvalue\nstatus\tcalled\nalpha_globin_genes\t3\n'
IV='locus\tname\tchrom\tstart\tend\tstatus\tcopies\texpected\tchange\n'
IV+='CFTR_del\tCFTR region\tchr7\t117500000\t117600000\tcalled\t1\t2\tdeletion\n'
IV+='DUP1\tx\tchr1\t1\t100\tcalled\t3\t2\tduplication\n'
IV+='GJB6\tx\tchr13\tx\t2\tcalled\t1\t2\tdeletion\n'

# full: every input, the bash layout, a run whose steps were all ok.
full() {
  local d=$1
  put "$d/clinvar/S_clinvar_hits.vcf" "$HITS"
  put "$d/clinvar/S_clinvar_cdna_check.tsv" "$CDNA"
  put "$d/clinvar/S_clinvar_cdna_scan.tsv" "$SCAN"
  put "$d/clinvar/S_clinvar_deep_intronic.tsv" "$DEEP"
  mkdir -p "$d/vep"
  ln -s ../../cb_vep.vcf.gz "$d/vep/S_vep.vcf.gz"
  ln -s ../../cb_vep.vcf.gz.tbi "$d/vep/S_vep.vcf.gz.tbi"
  put "$d/pharmcat/S.report.json" '{"pharmcatVersion": "3.4.0", "genes": {}}'
  put "$d/pharmcat/S.report.html" '<html></html>'
  put "$d/cpic/S_phenotypes.tsv" "$CPIC"
  put "$d/cpic/S_cpic_recommendations.txt" "$REC"
  put "$d/pypgx/S_pypgx_summary.tsv" "$PYPGX"
  put "$d/pypgx/S_cyp2d6_depth_check.tsv" "$DEPTH"
  put "$d/pypgx/S_pharmcat_comparison.tsv" "$COMP"
  put "$d/cyrius/S_cyp2d6.tsv" "${CYRIUS_HEADER}S\t*2/*1\tPASS\n"
  put "$d/pgx_consensus/S_pgx_consensus.tsv" "$CONS"
  put "$d/cpsr/S.cpsr.grch38.html" '<html></html>'
  putgz "$d/cpsr/S.cpsr.grch38.classification.tsv.gz" "$CPSR"
  put "$d/cpsr/S_cpsr_cdna_check.tsv" "$CPSR_CDNA"
  put "$d/clinical/S_clinical_summary.tsv" "$CLINT"
  putgz "$d/clinical/S_clinical.vcf.gz" "$CLIN"
  put "$d/slivar/S_slivar_summary.tsv" "$SLIV"
  put "$d/slivar/S_compound_hets.tsv" "$CHETS"
  put "$d/slivar/S_compound_het_pairs.tsv" "$PAIRS"
  put "$d/annotsv/S_sv_class4_5_checked.tsv" "$ANN"
  putgz "$d/sv_merged/S_sv_consensus.vcf.gz" "$CONSV"
  put "$d/paralogs/S_smn_copy_number.tsv" "$SMN"
  put "$d/paralogs/S_hba_copy_number.tsv" "$HBA"
  put "$d/paralogs/S_cn_intervals.tsv" "$IV"
  put "$d/logs/run_status.tsv" "meta\tstarted_epoch\t${T0}\nstep\t06\tok\nstep\t07\tok\nstep\t17\tok\nstep\t21\tok\nstep\t23\tok\nstep\t27\tok\nstep\t31\tok\nstep\t32\tok\n"
}
full cb-full

# Change one header at a time. Quoted TSV fields have the same names as
# the unquoted baseline, including the dosage table's validation columns.
full cb-quoted-smn
put cb-quoted-smn/paralogs/S_smn_copy_number.tsv "\"chrom\"${SMN#chrom}"
full cb-quoted-interval
put cb-quoted-interval/paralogs/S_cn_intervals.tsv "\"locus\"${IV#locus}"
full cb-quoted-dosage
put cb-quoted-dosage/paralogs/S_cn_dosage.tsv '"locus"\tstatus\tcopies\testimate\treason\treading\n'

# A valid UTF-8 contig with a code point crossing byte offset three must
# leave the consensus readable. The last chr17 record still matches.
full cb-unicode-contig
putgz cb-unicode-contig/sv_merged/S_sv_consensus.vcf.gz "${CONSV//chr1\\t/éé\\t}"

# nf: the Nextflow names. The clinical VCF and slivar's VCF without their
# tables, two PharmCAT reports (the newer wins), Cyrius under another name,
# CPSR by its wildcard names, no VEP output, a consensus older than pypgx,
# PharmCAT's step failed, CPIC from an earlier run, pypgx failed.
put cb-nf/clinvar/S_clinvar_hits.vcf "$HITS"
putgz cb-nf/clinical/S_clinical.vcf.gz "$CLIN"
putgz cb-nf/slivar/S_prioritized.vcf.gz "${VCF_HEAD}chr1\t1\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\nchr1\t2\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\n"
put cb-nf/slivar/S_compound_het_pairs.tsv "$PAIRS"
put cb-nf/vcf/A.report.json '{"pharmcatVersion": 3.2}'
put cb-nf/vcf/B.report.json '[1, 2]'
put cb-nf/cpic/S_phenotypes.tsv 'Gene\tDiplotype\tPhenotype\tStatus\nCYP2D6\t*4/*4\tPoor\tnot called\n'
put cb-nf/cyrius/T_cyp2d6.tsv "${CYRIUS_HEADER}T\t*1/*2\tLowQ\n"
put cb-nf/pypgx/S_pypgx_summary.tsv "$PYPGX"
put cb-nf/pgx_consensus/S_pgx_consensus.tsv "$CONS"
put cb-nf/cpsr/T.cpsr.grch37.html '<html></html>'
putgz cb-nf/cpsr/T.cpsr.grch37.classification.tsv.gz 'SYMBOL\tFINAL_CLASSIFICATION\tCPSR_CLASSIFICATION\nA\tPathogenic\tBenign\nB\t\t\n'
put cb-nf/annotsv/S_sv_class4_5_checked.tsv ''
put cb-nf/paralogs/S_cn_intervals.tsv "$IV"
put cb-nf/logs/run_status.tsv "meta\tstarted_epoch\t${T0}\nstep\t07\tfailed (exit 1)\nstep\t27\tskipped\nstep\t32\tok\nstep\t21\tok\n"

# probable: pypgx alone, Cyrius not usable, the depth check ok; no consensus.
put cb-probable/pypgx/S_pypgx_summary.tsv "$PYPGX"
put cb-probable/pypgx/S_cyp2d6_depth_check.tsv "$DEPTH"
put cb-probable/cyrius/S_cyp2d6.tsv "${CYRIUS_HEADER}S\t*1/*3\tLowQ\n"
put cb-probable/cpic/S_phenotypes.tsv 'Gene\tDiplotype\tPhenotype\tStatus\nCYP2C9\t*1/*1\tNormal\tnormal\n'
put cb-probable/cpic/S_cpic_recommendations.txt 'PARSING FAILED for one gene\n'

# failed: CPIC, pypgx and Cyrius failed in the run and wrote nothing.
put cb-failed/logs/run_status.tsv "meta\tstarted_utc\t2020-02-02 02:02 UTC\nstep\t21\tfailed (exit 1)\nstep\t27\tfailed\nstep\t32\tfailed (exit 3)\nstep\t06\tok\n"
put cb-failed/pgx_consensus/S_pgx_consensus.tsv 'Gene\tResult\nCYP2C19\tx\n'

# bad: a hit whose POS is no number, a report that is no JSON, a CPSR
# table cut short, pypgx without a Gene column, an AnnotSV table without
# its columns, a Parascopy table without its columns, a clinical floor that
# is no number, slivar's VCF cut short.
put cb-bad/clinvar/S_clinvar_hits.vcf "${VCF_HEAD}chr1\tx\t.\tA\tG\t.\tPASS\tGENEINFO=MUTYH:1;CLNSIG=Pathogenic\tGT\t0/1\n"
put cb-bad/pharmcat/S.report.json '{"pharmcatVersion": '
putgz cb-bad/cpsr/S.cpsr.grch38.classification.tsv.gz "$CPSR"
head -c 60 cb-bad/cpsr/S.cpsr.grch38.classification.tsv.gz > cb-bad.tmp && mv cb-bad.tmp cb-bad/cpsr/S.cpsr.grch38.classification.tsv.gz
put cb-bad/pypgx/S_pypgx_summary.tsv 'Name\tSource\tCNV_call\nCYP2D6\tBAM\tNormal\n'
put cb-bad/annotsv/S_sv_class4_5_checked.tsv 'SV_chrom\tSV_start\n1\t1\n'
put cb-bad/paralogs/S_smn_copy_number.tsv 'chrom\tstart\tend\nchr5\t1\t2\n'
putgz cb-bad/clinical/S_clinical.vcf.gz "##PGP_CALL_FLOOR=GQ=20,HOM_AB=0..8\n#CHROM\n"
put cb-bad/clinical/S_clinical_summary.tsv 'GENE\tIMPACT\nA\tHIGH\n'
put cb-bad/cyrius/S_cyp2d6.tsv ''
put cb-bad/cpic/S_phenotypes.tsv 'Gene\tDiplotype\n'
putgz cb-bad/slivar/S_prioritized.vcf.gz "$CLIN"
head -c 100 cb-bad/slivar/S_prioritized.vcf.gz > cb-bad.tmp && mv cb-bad.tmp cb-bad/slivar/S_prioritized.vcf.gz

# nf-bad: the Nextflow clinical VCF cut short; Parascopy's catalog tables
# without their catalog; an HBA table without a status.
put cb-nfbad/clinvar/S_clinvar_hits.vcf "$HITS"
putgz cb-nfbad/clinical/S_clinical.vcf.gz "$CLIN"
head -c 100 cb-nfbad/clinical/S_clinical.vcf.gz > cb-nfbad.tmp && mv cb-nfbad.tmp cb-nfbad/clinical/S_clinical.vcf.gz
put cb-nfbad/paralogs/S_hba_copy_number.tsv 'metric\tvalue\nalpha_globin_genes\t2\n'
put cb-nfbad/paralogs/S_rhd_copy_number.tsv "$SMN"

# notbi: the VEP output without its index; the HBA count of no gene left.
put cb-notbi/clinvar/S_clinvar_hits.vcf "$HITS"
mkdir -p cb-notbi/vep
ln -s ../../cb_vep.vcf.gz cb-notbi/vep/S_vep.vcf.gz
put cb-notbi/paralogs/S_hba_copy_number.tsv 'metric\tvalue\nstatus\tcalled\nalpha_globin_genes\t0\n'

# notgz: the clinical and slivar VCFs as plain text under a .gz name. The
# Python tells gzip by the name and finds both unreadable; og tells it by
# the first bytes and reads them (a golden).
put cb-notgz/clinical/S_clinical.vcf.gz "$CLIN"
put cb-notgz/slivar/S_prioritized.vcf.gz "${VCF_HEAD}chr1\t1\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\n"

mkdir -p cb-empty

# Every file after the run started; the CPIC table of nf from before it;
# nf's second PharmCAT report and its consensus' sources newer.
find . -type f -exec touch -t "$NEW" {} +
at "$OLD" cb-nf/cpic/S_phenotypes.tsv cb-nf/pgx_consensus/S_pgx_consensus.tsv
at "$NEWER" cb-nf/vcf/B.report.json
