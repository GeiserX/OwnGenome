#!/usr/bin/env bash
# ClinVar records labelled "Conflicting classifications of pathogenicity" that
# some labs call P/LP, end to end: setup's subset rule (og clinvar
# conflicting-subset in the bcftools image, as setup.sh runs it), step 06's
# conflicting table, the summary's sections.clinvar.data.conflicting, and the
# CLINVAR_SCREEN module.
#   - A small ClinVar at real chr22 positions of the fixture reference holds
#     three records the sample carries. Their INFO tags are copied as ClinVar
#     writes them from three records of NCBI's 2026-11-22 release (public):
#     481516, a conflicting BRCA2 record with Pathogenic_(1)|Likely_pathogenic_(1)|
#     Uncertain_significance_(2) (kept by the subset); 882070, a conflicting BRCA2
#     record with Uncertain_significance_(2)|Likely_benign_(3) only (dropped);
#     and 53622, a plain Pathogenic CFTR record (the pathogenic subset's, never
#     the conflicting one's).
#   - With the subset, step 06 writes one conflicting row with the submitter
#     split, the stars and the condition, and the plain record is the only hit.
#   - The tier adds, it never changes: without the subset the hits VCF and TSV
#     are byte-identical, and so are the summary's clinvar hits, ACMG card,
#     carrier list and couple section; the table then says "not checked", never
#     zero rows. The couple screen reads only the hits VCF.
#   - The module writes the same table as the bash step, and "not checked"
#     without --clinvar_conflicting; the Nextflow summary.json lists the record
#     under sections.clinvar.data.conflicting, as og summary collect does on
#     the bash run.
. "$(dirname "$0")/lib.sh"

S2=ccf
G2="${E2E_WORK}/ccf-genome"   # a GENOME_DIR of its own: its ClinVar files are this case's
rm -rf "${G2:?}"
mkdir -p "${G2}/reference" "${G2}/clinvar" "${G2}/${S2}/vcf"
for f in GRCh38_no_alt_analysis_set.fasta GRCh38_no_alt_analysis_set.fasta.fai GRCh38_no_alt_analysis_set.dict; do
  ln "${GENOME_DIR}/reference/${f}" "${G2}/reference/${f}" 2>/dev/null || cp "${GENOME_DIR}/reference/${f}" "${G2}/reference/${f}"
done
export GENOME_DIR="$G2"   # lib.sh's bcf and sam, and the steps, read this one now
U="$(id -u):$(id -g)"
# in_g2 IMAGE CMD...: a command in IMAGE as the runner's user, with bin/ as /pgp-bin.
in_g2() { local img=$1; shift; docker run --rm --user "$U" -v "${G2}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w /genome "$img" "$@"; }

# Three SNVs at real positions: REF is the reference's base, ALT another one.
POS=(30000300 30000400 30000500)
declare -A REFB ALTB
for p in "${POS[@]}"; do
  b=$(sam faidx reference/GRCh38_no_alt_analysis_set.fasta "chr22:${p}-${p}" | sed 1d | tr -d '\n' | tr '[:lower:]' '[:upper:]')
  REFB[$p]=$b
  case "$b" in A) ALTB[$p]=C ;; [CGT]) ALTB[$p]=A ;; *) ALTB[$p]=N ;; esac
done
check "the reference has a base at each of the three positions (${REFB[30000300]}${REFB[30000400]}${REFB[30000500]})" \
  has '^[ACGT]{3}$' "${REFB[30000300]}${REFB[30000400]}${REFB[30000500]}"

KEPT=30000300 BENIGN=30000400 PLAIN=30000500
{
  printf '##fileformat=VCFv4.1\n##fileDate=20261101\n##contig=<ID=chr22,length=50818468>\n'
  printf '##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s) for the variant reported as gene symbol:gene id">\n'
  printf '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Aggregate germline classification for this single variant">\n'
  printf '##INFO=<ID=CLNSIGCONF,Number=.,Type=String,Description="Conflicting germline classification for this single variant">\n'
  printf '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status of the germline classification">\n'
  printf '##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar preferred disease name">\n'
  printf '##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr22\t%s\t481516\t%s\t%s\t.\t.\tGENEINFO=BRCA2:675;CLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Pathogenic_(1)|Likely_pathogenic_(1)|Uncertain_significance_(2);CLNREVSTAT=criteria_provided,_conflicting_classifications;CLNDN=Hereditary_cancer-predisposing_syndrome|Hereditary_breast_ovarian_cancer_syndrome|Breast-ovarian_cancer,_familial,_susceptibility_to,_2|Familial_cancer_of_breast;CLNDISDB=MONDO:MONDO:0015356,MeSH:D009386,MedGen:C0027672,Orphanet:140162|MONDO:MONDO:0003582,MeSH:D061325,MedGen:C0677776,Orphanet:145|MONDO:MONDO:0012933,MedGen:C2675520,OMIM:612555,Orphanet:145|MONDO:MONDO:0016419,MedGen:C0346153,OMIM:114480,Orphanet:227535\n' \
    "$KEPT" "${REFB[$KEPT]}" "${ALTB[$KEPT]}"
  printf 'chr22\t%s\t882070\t%s\t%s\t.\t.\tGENEINFO=BRCA2:675|LOC106721785:106721785;CLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Uncertain_significance_(2)|Likely_benign_(3);CLNREVSTAT=criteria_provided,_conflicting_classifications\n' \
    "$BENIGN" "${REFB[$BENIGN]}" "${ALTB[$BENIGN]}"
  printf 'chr22\t%s\t53622\t%s\t%s\t.\t.\tGENEINFO=CFTR:1080;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNDN=CFTR-related_disorder|Congenital_bilateral_aplasia_of_vas_deferens_from_CFTR_mutation|Cystic_fibrosis|not_provided|Bronchiectasis_with_or_without_elevated_sweat_chloride_1|Hereditary_pancreatitis;CLNDISDB=MONDO:MONDO:7770004,MedGen:C5924204|MONDO:MONDO:0010178,MedGen:C0403814,OMIM:277180,Orphanet:48|MONDO:MONDO:0009061,MedGen:C0010674,OMIM:219700,Orphanet:586|MedGen:C3661900|MONDO:MONDO:0008887,MedGen:C2749757,OMIM:211400,Orphanet:60033|MONDO:MONDO:0008185,MedGen:C0238339,OMIM:167800,Orphanet:676\n' \
    "$PLAIN" "${REFB[$PLAIN]}" "${ALTB[$PLAIN]}"
} > "${G2}/clinvar/clinvar_chr.vcf"
{
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr22,length=50818468>\n##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
  printf '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\n'
  for p in "${POS[@]}"; do printf 'chr22\t%s\t.\t%s\t%s\t50\tPASS\t.\tGT:GQ:DP\t0/1:50:40\n' "$p" "${REFB[$p]}" "${ALTB[$p]}"; done
} > "${G2}/${S2}/vcf/${S2}.vcf"

# The two subsets, built as setup.sh builds them: the pathogenic one with
# bcftools' CLNSIG rule, the conflicting one with og clinvar conflicting-subset.
in_g2 "$BCFTOOLS_IMAGE" bash -o pipefail -c "set -e
  bcftools view -Oz -o ${S2}/vcf/${S2}.vcf.gz ${S2}/vcf/${S2}.vcf
  bcftools index -f -t ${S2}/vcf/${S2}.vcf.gz
  bcftools view -Oz -o clinvar/clinvar_chr.vcf.gz clinvar/clinvar_chr.vcf
  bcftools index -f -t clinvar/clinvar_chr.vcf.gz
  bcftools view -i 'CLNSIG~\"Pathogenic\" || CLNSIG~\"Likely_pathogenic\"' clinvar/clinvar_chr.vcf.gz -Oz -o clinvar/clinvar_pathogenic_chr.vcf.gz
  bcftools index -f -t clinvar/clinvar_pathogenic_chr.vcf.gz
  /pgp-bin/og clinvar conflicting-subset clinvar/clinvar_chr.vcf.gz | bcftools view --no-version -Oz -o clinvar/clinvar_conflicting_plp_chr.vcf.gz -
  bcftools index -f -t clinvar/clinvar_conflicting_plp_chr.vcf.gz"
check_eq "the conflicting subset holds only the record two labs call Pathogenic" \
  "$(bcf query -f '%ID\n' clinvar/clinvar_conflicting_plp_chr.vcf.gz 2>/dev/null | paste -sd, -)" 481516
check_eq "the pathogenic subset holds only the plain Pathogenic record" \
  "$(bcf query -f '%ID\n' clinvar/clinvar_pathogenic_chr.vcf.gz 2>/dev/null | paste -sd, -)" 53622

OUT="${G2}/${S2}/clinvar"
CONF="${OUT}/${S2}_clinvar_conflicting.tsv"
# summary_part FILE: the parts of summary.json the tier must never change.
summary_part() {
  jq -S '{clinvar: (.sections.clinvar.data | del(.conflicting)), acmg_sf: .sections.acmg_sf, carrier: .sections.carrier, couple: .sections.couple}' "$1"
}
collect() {
  "${REPO}/bin/og" summary collect --sample "$S2" --sample-dir "${G2}/${S2}" --assets "${REPO}/assets" \
    --out "$1" > "${CASE_TMP}/collect.log" 2>&1
}

# --- step 06 with the subset ------------------------------------------------------
run_step 06-clinvar-screen.sh "$S2"
check_step_exit 06-clinvar-screen.sh
check_eq "one conflicting row" "$(awk 'NR > 1' "$CONF" 2>/dev/null | grep -c . || true)" 1
ROW=$(awk -F'\t' -v p="$KEPT" 'NR > 1 && $2 == p' "$CONF" 2>/dev/null)
check "the row is the kept record, with its gene, ClinVar id and 1 star" \
  has "^chr22	${KEPT}	${REFB[$KEPT]}	${ALTB[$KEPT]}	het	BRCA2	481516	1	" "${ROW:-none}"
check "the row has the submitter split (without ClinVar's _ before each count) and two P/LP submissions" \
  has '	Pathogenic\(1\), Likely_pathogenic\(1\), Uncertain_significance\(2\)	2	' "${ROW:-none}"
check_eq "the row names the conditions" "$(cut -f12 <<< "${ROW:-none}")" \
  'Hereditary cancer-predisposing syndrome / Hereditary breast ovarian cancer syndrome / Breast-ovarian cancer, familial, susceptibility to, 2 / Familial cancer of breast'
check_eq "the row lists each condition id once" "$(cut -f13 <<< "${ROW:-none}")" \
  'MONDO:MONDO:0015356,MeSH:D009386,MedGen:C0027672,Orphanet:140162,MONDO:MONDO:0003582,MeSH:D061325,MedGen:C0677776,Orphanet:145,MONDO:MONDO:0012933,MedGen:C2675520,OMIM:612555,MONDO:MONDO:0016419,MedGen:C0346153,OMIM:114480,Orphanet:227535'
check_eq "the row's depth, quality and confidence" "$(cut -f14- <<< "${ROW:-none}")" "$(printf '40\t50\tno')"
check_eq "the table's header" "$(head -1 "$CONF" 2>/dev/null)" \
  "$(printf 'chrom\tpos\tref\talt\tgenotype\tgene\tclinvar_id\tstars\treview_status\tsubmissions\tplp_submissions\tcondition\tcondition_ids\tdp\tgq\tlow_confidence')"
check_eq "the only hit is the plain Pathogenic record" \
  "$(grep -v '^#' "${OUT}/${S2}_clinvar_hits.vcf" 2>/dev/null | cut -f2,3 | paste -sd, -)" "${PLAIN}	53622"
check "the step's log counts one conflicting record apart from the hits" \
  grep -q 'Conflicting records some labs call pathogenic (review with a clinician; not hits): 1 in ' "$STEP_LOG"
check "the step's log still counts one pathogenic hit" grep -q '^Count: 1 pathogenic hits$' "$STEP_LOG"
cp "${OUT}/${S2}_clinvar_hits.vcf" "${CASE_TMP}/with.vcf"
cp "${OUT}/${S2}_clinvar_hits.tsv" "${CASE_TMP}/with.tsv"
cp "$CONF" "${CASE_TMP}/bash_conflicting.tsv"
collect "${CASE_TMP}/with.json"
check_eq "og summary collect exits 0 (with the subset)" "$?" 0
check_eq "summary: sections.clinvar.data.conflicting is checked with one record" \
  "$(jq -c '.sections.clinvar.data.conflicting | [.checked, .count, .hits[0].clinvar_id, .hits[0].plp_submissions, .hits[0].clinvar_condition]' "${CASE_TMP}/with.json" 2>/dev/null)" \
  '[true,1,"481516",2,"Hereditary cancer-predisposing syndrome / Hereditary breast ovarian cancer syndrome / Breast-ovarian cancer, familial, susceptibility to, 2 / Familial cancer of breast"]'
check_eq "summary: the clinvar hit count leaves the conflicting record out" \
  "$(jq -c '.sections.clinvar.data.count' "${CASE_TMP}/with.json" 2>/dev/null)" 1
check "summary: BRCA2 is in neither the ACMG card nor the carrier list" \
  lacks BRCA2 "$(jq -c '[.sections.acmg_sf, .sections.carrier]' "${CASE_TMP}/with.json" 2>/dev/null)"
check "summary.json validates against tests/schema/summary.schema.json" \
  "${REPO}/bin/og" dev schema "${REPO}/tests/schema/summary.schema.json" "${CASE_TMP}/with.json"

# --- step 06 without the subset (an old setup) -----------------------------------
mkdir -p "${CASE_TMP}/aside"
mv "${G2}/clinvar/clinvar_conflicting_plp_chr.vcf.gz" "${G2}/clinvar/clinvar_conflicting_plp_chr.vcf.gz.tbi" "${CASE_TMP}/aside/"
run_step 06-clinvar-screen.sh "$S2"
check_step_exit 06-clinvar-screen.sh
check_eq "without the subset the table says not checked" "$(head -1 "$CONF" 2>/dev/null)" \
  '# not checked: run setup.sh --refresh clinvar'
check_eq "without the subset the table has no row" "$(grep -v '^#' "$CONF" 2>/dev/null | awk 'NR > 1' | grep -c . || true)" 0
check "the step's log says not checked" \
  grep -q 'Conflicting records some labs call pathogenic: not checked: run setup.sh --refresh clinvar' "$STEP_LOG"
# bcftools stamps "; Date=..." into its command header lines, so compare without it.
check "the hits VCF is byte-identical with and without the tier, apart from bcftools date stamps" \
  cmp -s <(sed "s/; Date=.*//" "${CASE_TMP}/with.vcf") <(sed "s/; Date=.*//" "${OUT}/${S2}_clinvar_hits.vcf")
check "the hits TSV is byte-identical with and without the tier" cmp -s "${CASE_TMP}/with.tsv" "${OUT}/${S2}_clinvar_hits.tsv"
collect "${CASE_TMP}/without.json"
check_eq "og summary collect exits 0 (without the subset)" "$?" 0
check_eq "summary: not checked, never zero" \
  "$(jq -c '.sections.clinvar.data.conflicting | [.checked, .count, (.hits | length), .note]' "${CASE_TMP}/without.json" 2>/dev/null)" \
  '[false,null,0,"not checked: run setup.sh --refresh clinvar"]'
summary_part "${CASE_TMP}/with.json" > "${CASE_TMP}/with.part.json" 2>/dev/null
summary_part "${CASE_TMP}/without.json" > "${CASE_TMP}/without.part.json" 2>/dev/null
check_ge "the compared summary parts are not empty" "$(wc -c < "${CASE_TMP}/with.part.json" | tr -d ' ')" 100
check "the summary's clinvar hits, ACMG card, carrier list and couple section are byte-identical with and without the tier" \
  cmp -s "${CASE_TMP}/with.part.json" "${CASE_TMP}/without.part.json"
mv "${CASE_TMP}/aside/"* "${G2}/clinvar/"

# --- the Nextflow module and report ----------------------------------------------
if command -v nextflow >/dev/null; then
  SHEET="${CASE_TMP}/sheet.csv"
  printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$S2" "${G2}/${S2}/vcf/${S2}.vcf.gz" "${G2}/${S2}/vcf/${S2}.vcf.gz.tbi" > "$SHEET"
  # nf_run OUTDIR [ARGS...]: main.nf with --tools clinvar,html_report on this
  # case's sample (HTML_REPORT writes <outdir>/<sample>/summary.json).
  nf_run() {
    local out=$1; shift
    echo "+ nextflow run main.nf --tools clinvar,html_report $*"
    (
      cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
        -work-dir "${E2E_WORK}/nf-work-ccf" --input "$SHEET" \
        --reference "${G2}/reference/GRCh38_no_alt_analysis_set.fasta" --tools clinvar,html_report \
        --clinvar "${G2}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
        --clinvar_index "${G2}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
        --outdir "$out" --max_cpus 4 --max_memory 14.GB "$@"
    ) > "${CASE_TMP}/nf.log" 2>&1
    NF_RC=$?
    sed 's/^/    | /' "${CASE_TMP}/nf.log" | grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists'
  }
  nf_run "${G2}/nf" --clinvar_conflicting "${G2}/clinvar/clinvar_conflicting_plp_chr.vcf.gz"
  check_eq "nextflow run with --clinvar_conflicting exits 0" "$NF_RC" 0
  check "CLINVAR_SCREEN writes the same conflicting table as step 06" \
    cmp -s "${G2}/nf/${S2}/clinvar/${S2}_clinvar_conflicting.tsv" "${CASE_TMP}/bash_conflicting.tsv"
  check_eq "CLINVAR_SCREEN hits (the plain record only)" \
    "$(grep -v '^#' "${G2}/nf/${S2}/clinvar/${S2}_clinvar_hits.vcf" 2>/dev/null | cut -f2,3 | paste -sd, -)" "${PLAIN}	53622"
  # The Nextflow summary.json carries the table as og summary collect does on
  # the bash run: the same section, field for field.
  check_eq "the Nextflow summary.json lists the record under sections.clinvar.data.conflicting, as the bash summary does" \
    "$(jq -cS '.sections.clinvar.data.conflicting' "${G2}/nf/${S2}/summary.json" 2>/dev/null || echo missing)" \
    "$(jq -cS '.sections.clinvar.data.conflicting' "${CASE_TMP}/with.json" 2>/dev/null)"
  check_eq "the Nextflow summary.json counts one conflicting record" \
    "$(jq -c '.sections.clinvar.data.conflicting | [.checked, .count]' "${G2}/nf/${S2}/summary.json" 2>/dev/null)" '[true,1]'
  nf_run "${G2}/nf-none"
  check_eq "nextflow run without --clinvar_conflicting exits 0" "$NF_RC" 0
  check "CLINVAR_SCREEN without --clinvar_conflicting says not checked" \
    grep -q '^# not checked: no --clinvar_conflicting' "${G2}/nf-none/${S2}/clinvar/${S2}_clinvar_conflicting.tsv"
  check "CLINVAR_SCREEN writes the same hits VCF without the tier, apart from bcftools date stamps" \
    cmp -s <(sed "s/; Date=.*//" "${G2}/nf/${S2}/clinvar/${S2}_clinvar_hits.vcf") <(sed "s/; Date=.*//" "${G2}/nf-none/${S2}/clinvar/${S2}_clinvar_hits.vcf")
  check_eq "the Nextflow summary.json without --clinvar_conflicting says not checked, never zero" \
    "$(jq -c '.sections.clinvar.data.conflicting | [.checked, .count, (.hits | length), .note]' "${G2}/nf-none/${S2}/summary.json" 2>/dev/null)" \
    '[false,null,0,"not checked: no --clinvar_conflicting; run setup.sh --refresh clinvar and pass clinvar/clinvar_conflicting_plp_chr.vcf.gz"]'
else
  fail "nextflow is not on PATH"
fi

finish
