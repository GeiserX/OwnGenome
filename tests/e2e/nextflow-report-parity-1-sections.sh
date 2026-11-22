#!/usr/bin/env bash
# The Nextflow report carries the sections the bash report has: the RFC1
# motif count (step 9b), the chrM ClinVar and gnomAD annotation (step 20),
# the cDNA read check (06), the CYP2D6 depth check (32), the phased
# compound-het pairs (31) and AnnotSV's class 4 and 5 SVs (05); and VEP gets
# the same normalised records from VEP_NORM as from step 13.
#   1. VEP_NORM against step 13, on the fixture VCF with one planted 1/2
#      record (so the split has something to do). VEP itself is a docker
#      wrapper that keeps the file step 13 gives it (its cache is 26 GB).
#      The records must match byte for byte, and the header but bcftools'
#      command lines.
#   2. Sample NFPAR: step 20 on the fixture BAM with ClinVar and the pinned
#      gnomAD chrM file, then step 9b on made-up reads at RFC1 (the reads of
#      og dev fixture rfc1 pathogenic, which Stranger flags), plus the tables of
#      earlier cases (cDNA check, depth check, pairs) and a class 4-5 table.
#      STRANGER_RFC1 and MITO_VARIANTS_ANNOTATE, run for real, must write
#      the tables steps 9b and 20 wrote, byte for byte.
#   3. Step 24 on NFPAR and HTML_REPORT on the same tables (the RFC1 and
#      chrM ones from the Nextflow processes; step 20's mtDNA copy number and
#      deletion tables, and the confirmed mtDNA list, as main.nf stages them):
#      the clinvar, pypgx, annotsv,
#      expansions, mito and slivar sections of the two summary.json files
#      must be equal. Control: the comparison reports a changed section.
#   4. main.nf passes every new table on: a stub run of every tool
#      (-profile test_all) stages each one in each HTML_REPORT task that has
#      the step's main output, and every VEP task reads VEP_NORM's file.
#      Control: a copy of the pipeline without the five joins and with VEP
#      reading the raw VCF stages none.
# What this cannot show: real carriers of an RFC1 expansion, and real class 4-5
# SVs; the reads and the SV table are made up.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

S=NFPAR
SD="${GENOME_DIR}/${S}"
REF="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta"
CV="${GENOME_DIR}/clinvar/clinvar_chr.vcf.gz"
GN="${GENOME_DIR}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
HG_VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
HG_BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
C="${SAMPLE}cdna"
for f in "$REF" "${REF}.fai" "$CV" "${CV}.tbi" "$HG_VCF" "$HG_BAM" \
         "${GENOME_DIR}/${C}/clinvar/${C}_clinvar_hits.vcf" "${GENOME_DIR}/${C}/clinvar/${C}_clinvar_cdna_check.tsv" \
         "${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv" "${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_cyp2d6_depth_check.tsv" \
         "${GENOME_DIR}/CHPH/slivar/CHPH_prioritized.vcf.gz" "${GENOME_DIR}/CHPH/slivar/CHPH_compound_het_pairs.tsv"; do
  check "an earlier case or the fixture left ${f#"$GENOME_DIR"/}" test -s "$f"
done

# A two-line pipeline in its own project folder, which gets the repository's bin/
NP="${CASE_TMP}/project"
rm -rf "$NP"
mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
# nf_run NAME [CONFIG LINES]: run NP/main.nf with the real containers, output under CASE_TMP/NAME
nf_run() {
  local name=$1 out="${CASE_TMP}/$1"; shift
  rm -rf "$out"
  {
    echo "params.outdir = '${out}/results'"
    echo "params.publish_dir_mode = 'copy'"
    echo "docker.enabled = true"
    echo "docker.runOptions = '--user $(id -u):$(id -g) --network none'"
    printf '%s\n' "$@"
  } > "${NP}/nextflow.config"
  echo "+ nextflow run (${name})"
  (
    cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false \
      -c "${REPO}/conf/containers.config" -work-dir "${CASE_TMP}/work-${name}"
  ) > "${CASE_TMP}/${name}.log" 2>&1
  local rc=$?
  check_eq "nextflow run (${name}) exits 0" "$rc" 0
  [ "$rc" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/${name}.log" | head -40
}

# --- 1. VEP_NORM against step 13 ---------------------------------------------------
GD="${CASE_TMP}/genome"
rm -rf "$GD"
mkdir -p "${GD}/reference" "${GD}/${S}/vcf" "${GD}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
echo "species homo_sapiens" > "${GD}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt"
for x in "" .fai; do ln -f "${REF}${x}" "${GD}/reference/ref.fa${x}" 2>/dev/null || cp "${REF}${x}" "${GD}/reference/ref.fa${x}"; done
# The planted record: a 1/2 SNV at the first free chr20 position from 10,250,001,
# with the reference's own base as REF (read through the .fai offsets).
POS=$(gzip -dc "$HG_VCF" | awk -F'\t' '$1 == "chr20" {t[$2 + 0] = 1} END {p = 10250001; while (p in t) p++; print p}')
OFF=$(awk -F'\t' -v p="$POS" '$1 == "chr20" {printf "%.0f\n", $3 + int((p - 1) / $4) * $5 + (p - 1) % $4}' "${REF}.fai")
B=$(dd if="$REF" bs=1 skip="$OFF" count=1 2>/dev/null | tr '[:lower:]' '[:upper:]')
ALTS=$(printf '%s\n' A C G T | grep -vxF "$B" | head -n 2 | paste -sd, -)
# The record goes before the first later chr20 record, or before chr22, chrX,
# chrY or chrM, or at the end.
gzip -dc "$HG_VCF" | awk -F'\t' -v pos="$POS" -v rec="chr20\t${POS}\t.\t${B}\t${ALTS}\t50\tPASS\t.\tGT:AD\t1/2:0,15,13" '
  !done && !/^#/ && (($1 == "chr20" && $2 + 0 > pos + 0) || $1 == "chr22" || $1 == "chrX" || $1 == "chrY" || $1 == "chrM") {
    print rec; done = 1
  }
  { print }
  END { if (!done) print rec }' > "${GD}/${S}/vcf/${S}.vcf"
echo "planted chr20:${POS} ${B}>${ALTS}" > "${CASE_TMP}/planted.txt"
cat "${CASE_TMP}/planted.txt"
docker run --rm --user "$(id -u):$(id -g)" -v "${GD}:/genome" -w /genome "$BCFTOOLS_IMAGE" \
  sh -c "bcftools view -Oz -o ${S}/vcf/${S}.vcf.gz ${S}/vcf/${S}.vcf && bcftools index -t ${S}/vcf/${S}.vcf.gz"
check "the VCF with the planted 1/2 record is indexed" test -s "${GD}/${S}/vcf/${S}.vcf.gz.tbi"
check_eq "it has the planted record and AD in its header" \
  "$(gzip -dc "${GD}/${S}/vcf/${S}.vcf.gz" | grep -c -e $'\t1/2:0,15,13$' -e '^##FORMAT=<ID=AD,')" 2

WRAP="${CASE_TMP}/wrap"
mkdir -p "$WRAP"
cat > "${WRAP}/docker" <<'WRAPPER'
#!/usr/bin/env bash
# The VEP image only keeps its input; every other call goes to the docker
# the PATH without this folder finds.
set -euo pipefail
self=$(cd "$(dirname "$0")" && pwd)
rest=""
IFS=: read -r -a dirs <<< "$PATH"
for d in "${dirs[@]}"; do
  [ "$d" = "$self" ] || rest+="${rest:+:}${d}"
done
export PATH="$rest"
args=("$@") in="" out="" vep=0
for i in "${!args[@]}"; do
  case "${args[$i]}" in
    "$WRAP_VEP_IMAGE") vep=1 ;;
    --input_file) in=${args[$((i + 1))]} ;;
    -o) out=${args[$((i + 1))]} ;;
  esac
done
if [ "$vep" = 1 ] && [ -n "$in" ] && [ -n "$out" ]; then
  cp "${WRAP_GENOME}${in#/genome}" "$WRAP_SEEN"
  cp "${WRAP_GENOME}${in#/genome}" "${WRAP_GENOME}${out#/genome}"
  exit 0
fi
exec docker "$@"
WRAPPER
chmod +x "${WRAP}/docker"
SEEN="${CASE_TMP}/step13_vep_input.vcf.gz"
rm -f "$SEEN"
# GENOME_DIR first: scripts/ci/e2e-run.sh reads a line that starts NAME_IMAGE="..." as an image to pull
GENOME_DIR="$GD" WRAP_VEP_IMAGE="$VEP_IMAGE" WRAP_GENOME="$GD" WRAP_SEEN="$SEEN" \
  REF_FASTA="${GD}/reference/ref.fa" PATH="${WRAP}:${PATH}" run_step 13-vep-annotation.sh "$S"
check_step_exit 13-vep-annotation.sh
check "step 13 gave VEP a file" test -s "$SEEN"

cat > "${NP}/main.nf" <<NF
include { VEP_NORM } from '${REPO}/modules/local/vep/main'

workflow {
    VEP_NORM(channel.of([[id: '${S}'], file('${GD}/${S}/vcf/${S}.vcf.gz'), file('${GD}/${S}/vcf/${S}.vcf.gz.tbi')]),
             channel.value(file('${GD}/reference/ref.fa')), channel.value(file('${GD}/reference/ref.fa.fai')))
}
NF
nf_run vepnorm "process { withName: 'VEP_NORM' { publishDir = [path: '${CASE_TMP}/vepnorm/results', mode: 'copy'] } }"
NORM="${CASE_TMP}/vepnorm/results/${S}.norm.vcf.gz"
check "VEP_NORM wrote and indexed its output" test -s "${NORM}.tbi"
body() { gzip -dc "$1" 2>/dev/null | grep -v '^#'; }
head_no_cmd() { gzip -dc "$1" 2>/dev/null | grep '^#' | grep -v '^##bcftools_'; }
N13=$(body "$SEEN" | wc -l | tr -d ' ')
NNF=$(body "$NORM" | wc -l | tr -d ' ')
echo "records given to VEP: step 13 ${N13}, VEP_NORM ${NNF}"
check_ge "step 13 gave VEP records" "$N13" 100
check "VEP_NORM's records are step 13's, byte for byte" cmp -s <(body "$SEEN") <(body "$NORM")
check "and so is the header but bcftools' command lines" cmp -s <(head_no_cmd "$SEEN") <(head_no_cmd "$NORM")
diff <(body "$SEEN") <(body "$NORM") | head -n 6
PL=$(awk '{print $2}' "${CASE_TMP}/planted.txt" | cut -d: -f2)
SPLIT=$(body "$NORM" | awk -F'\t' -v p="$PL" '$1 == "chr20" && $2 == p {split($9, f, ":"); split($10, v, ":"); for (i in f) if (f[i] == "AD") a = v[i]; print $5 "=" a}' | paste -sd';' -)
echo "the planted record after VEP_NORM: ${SPLIT:-nothing}"
check "VEP_NORM split the planted 1/2 record in two, each AD with the other ALT on REF" \
  has '^[ACGT]=13,15;[ACGT]=15,13$' "${SPLIT:-}"
check_eq "no multiallelic record reaches VEP" "$(body "$NORM" | awk -F'\t' '$5 ~ /,/' | wc -l | tr -d ' ')" 0
check_eq "no 0/0 half reaches VEP" "$(body "$NORM" | awk -F'\t' '{split($10, v, ":")} v[1] ~ /^0[\/|]0$/' | wc -l | tr -d ' ')" 0

# --- 2. NFPAR: steps 20 and 9b, then the two Nextflow processes ------------------------
rm -rf "$SD"
mkdir -p "${SD}/vcf" "${SD}/aligned" "${SD}/expansion_hunter"
ln -f "$HG_VCF" "${SD}/vcf/${S}.vcf.gz"
ln -f "${HG_VCF}.tbi" "${SD}/vcf/${S}.vcf.gz.tbi"
ln -f "$HG_BAM" "${SD}/aligned/${S}_sorted.bam"
ln -f "${HG_BAM}.bai" "${SD}/aligned/${S}_sorted.bam.bai"
GN_PLACED=""
if [ ! -s "$GN" ]; then
  mkdir -p "$(dirname "$GN")"
  curl -fsSL --retry 3 -o "${GN}.part" "$GNOMAD_MITO_URL"
  check_eq "the gnomAD chrM download has the pinned sha256" "$(sha256sum "${GN}.part" 2>/dev/null | awk '{print $1}')" "$GNOMAD_MITO_SHA256"
  mv -f "${GN}.part" "$GN"
  GN_PLACED=1
fi
run_step 20-mtoolbox.sh "$S"
check_step_exit 20-mtoolbox.sh
MT_VCF="${SD}/mito/${S}_chrM_filtered.vcf.gz"
MT_TAB="${SD}/mito/${S}_chrM_annotated.tsv"
check "step 20 wrote the annotation table" test -s "$MT_TAB"
check "step 20's table names the gnomAD file" has $'^# gnomad\tgnomad' "$(head -n 3 "$MT_TAB" 2>/dev/null)"

# Step 9b on made-up reads at RFC1 that Stranger flags (the BAM step 20 read is not needed any more)
rm -f "${SD}/aligned/${S}_sorted.bam" "${SD}/aligned/${S}_sorted.bam.bai"
# The header of the smoke VCF with an ALT line per allele size, one RFC1
# record, and the pathogenic read set of og dev fixture rfc1.
{
  grep '^#' "${REPO}/tests/smoke/stranger_input.vcf.in" | sed '$d'
  for a in 100 150; do printf '##ALT=<ID=STR%s,Description="Allele comprised of %s repeat units">\n' "$a" "$a"; done
  grep '^#' "${REPO}/tests/smoke/stranger_input.vcf.in" | tail -n 1
  printf 'chr4\t39348425\tRFC1\tA\t<STR100>,<STR150>\t.\tPASS\tSVTYPE=STR;END=39348479;REF=11;RL=55;RU=AARRG;VARID=RFC1;REPID=RFC1'
  printf '\tGT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC'
  # The sample's values on their own line: a fraction beside the REPCN key
  # reads as a repeat genotype to scripts/ci/check-personal-data.sh.
  printf '\t%s\n' "1/2:INREPEAT/INREPEAT:100/150:80-120/130-170:0/0:2/2:20/20:30.0"
} > "${SD}/expansion_hunter/${S}_eh.vcf"
"${REPO}/bin/og" dev fixture rfc1 pathogenic "${SD}/aligned/reads.sam"
docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w "/genome/${S}/aligned" "$SAMTOOLS_IMAGE" \
  sh -c "samtools sort -o ${S}_sorted.bam reads.sam && samtools index ${S}_sorted.bam"
rm -f "${SD}/aligned/reads.sam"
run_step 09b-stranger.sh "$S"
check_step_exit 09b-stranger.sh
RFC_TAB="${SD}/expansion_hunter/${S}_rfc1_motif.tsv"
verdict() { awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i} NR == 2 {print $c["verdict"]}' "$1" 2>/dev/null; }
check_eq "step 9b counts the RFC1 motif of the made-up reads" "$(verdict "$RFC_TAB")" pathogenic_motif

cat > "${NP}/main.nf" <<NF
include { STRANGER; STRANGER_RFC1 } from '${REPO}/modules/local/stranger/main'
include { MITO_VARIANTS_ANNOTATE } from '${REPO}/modules/local/mito_variants/main'

workflow {
    STRANGER(channel.of([[id: '${S}'], file('${SD}/expansion_hunter/${S}_eh.vcf')]), [])
    STRANGER_RFC1(STRANGER.out.vcf.map { meta, vcf ->
        [meta, vcf, file('${SD}/aligned/${S}_sorted.bam'), file('${SD}/aligned/${S}_sorted.bam.bai')] }, [])
    MITO_VARIANTS_ANNOTATE(channel.of([[id: '${S}'], file('${MT_VCF}')]),
        channel.value(file('${CV}')), channel.value(file('${CV}.tbi')), channel.value(file('${GN}')))
}
NF
nf_run tools
NF_DIR="${CASE_TMP}/tools/results/${S}"
check_eq "STRANGER_RFC1 counts the same motif" "$(verdict "${NF_DIR}/expansion_hunter/${S}_rfc1_motif.tsv")" pathogenic_motif
check "STRANGER_RFC1 writes step 9b's table, byte for byte" cmp -s "$RFC_TAB" "${NF_DIR}/expansion_hunter/${S}_rfc1_motif.tsv"
check "MITO_VARIANTS_ANNOTATE writes step 20's table, byte for byte" cmp -s "$MT_TAB" "${NF_DIR}/mito/${S}_chrM_annotated.tsv"
diff "$MT_TAB" "${NF_DIR}/mito/${S}_chrM_annotated.tsv" | head -n 6
check_ge "with gnomAD values in it" \
  "$(awk -F'\t' '!/^#/ && $1 != "chrom" && $11 != "." {n++} END {print n + 0}' "${NF_DIR}/mito/${S}_chrM_annotated.tsv" 2>/dev/null)" 1

# --- 3. step 24 and HTML_REPORT on the same tables ---------------------------------------
mkdir -p "${SD}/clinvar" "${SD}/pypgx" "${SD}/slivar" "${SD}/annotsv"
cp "${GENOME_DIR}/${C}/clinvar/${C}_clinvar_hits.vcf" "${SD}/clinvar/${S}_clinvar_hits.vcf"
cp "${GENOME_DIR}/${C}/clinvar/${C}_clinvar_cdna_check.tsv" "${SD}/clinvar/${S}_clinvar_cdna_check.tsv"
cp "${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv" "${SD}/pypgx/${S}_pypgx_summary.tsv"
cp "${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_cyp2d6_depth_check.tsv" "${SD}/pypgx/${S}_cyp2d6_depth_check.tsv"
cp "${GENOME_DIR}/CHPH/slivar/CHPH_prioritized.vcf.gz" "${SD}/slivar/${S}_prioritized.vcf.gz"
cp "${GENOME_DIR}/CHPH/slivar/CHPH_compound_het_pairs.tsv" "${SD}/slivar/${S}_compound_het_pairs.tsv"
printf '%s\n' \
  $'SV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tGene_name\tACMG_class\tCaller\tExon_edge\tExon_edge_detail\tPopulation_AF\tPopulation_source\tCommon' \
  $'chr20\t10100999\t10102999\t-2000\tDEL\tSNAP25\t5\tmanta\tno\t\t\t\tno' \
  $'chr22\t42100000\t42150000\t-50000\tDEL\tCYP2D6;CYP2D7\t4\tdelly\tyes\tNM_000106 exons 2-8\t0.004\tgnomAD-SV_v4.1\tno' \
  > "${SD}/annotsv/${S}_sv_class4_5_checked.tsv"
run_step 24-html-report.sh "$S"
check_step_exit 24-html-report.sh

cat > "${NP}/main.nf" <<NF
include { HTML_REPORT } from '${REPO}/modules/local/html_report/main'

workflow {
    HTML_REPORT(channel.of([[id: '${S}', sex: 'male'], file('${SD}/vcf/${S}.vcf.gz'), [
        file('${SD}/clinvar/${S}_clinvar_hits.vcf'), file('${SD}/clinvar/${S}_clinvar_cdna_check.tsv'),
        file('${SD}/pypgx/${S}_pypgx_summary.tsv'), file('${SD}/pypgx/${S}_cyp2d6_depth_check.tsv'),
        file('${SD}/slivar/${S}_prioritized.vcf.gz'), file('${SD}/slivar/${S}_compound_het_pairs.tsv'),
        file('${SD}/annotsv/${S}_sv_class4_5_checked.tsv'), file('${SD}/expansion_hunter/${S}_eh.vcf'),
        file('${NF_DIR}/expansion_hunter/${S}_eh_stranger.vcf'), file('${NF_DIR}/expansion_hunter/${S}_rfc1_motif.tsv'),
        file('${MT_VCF}'), file('${NF_DIR}/mito/${S}_chrM_annotated.tsv'),
        file('${SD}/mito/${S}_mtdna_copies.tsv'), file('${SD}/mito/${S}_chrM_deletions.tsv'),
        file('${REPO}/assets/mito_confirmed.tsv')]]))
}
NF
nf_run report
NF_SUM="${CASE_TMP}/report/results/${S}/summary.json"
check "HTML_REPORT wrote the report" test -s "${CASE_TMP}/report/results/${S}/${S}_report.html"
# same_sections BASH NF [CHANGE]: the six sections of the two summaries, equal but
# for the file times; CHANGE=1 adds a key to the NF expansions data first (the control)
same_sections() {
  jq -rn --slurpfile b "$1" --slurpfile n "$2" --arg change "${3:-0}" '
    $b[0].sections as $bs
    | ($n[0].sections | if $change == "1" then .expansions.data.control = 1 else . end) as $ns
    | [ "clinvar", "pypgx", "annotsv", "expansions", "mito", "slivar"
        | ($bs[.] | del(.file_date)) as $x | ($ns[.] | del(.file_date)) as $y
        | ($x == $y and $x.state == "ok") as $ok
        | {ok: $ok, line: "  \(.): bash \($x.state), nextflow \($y.state)\(if $ok then "" else ", differ" end)"} ]
    | (.[].line), (if all(.[]; .ok) then empty else "" | halt_error(1) end)'
}
check "the six sections of the Nextflow summary equal the bash ones" same_sections "${SD}/summary.json" "$NF_SUM"
same_sections "${SD}/summary.json" "$NF_SUM" 1 > /dev/null
CTRL=$?
check "control: a changed expansions section is told apart" test "$CTRL" != 0
HTML="${CASE_TMP}/report/results/${S}/${S}_report.html"
check "the Nextflow HTML report shows the RFC1 motif verdict" grep -q 'RFC1 motif count: pathogenic_motif' "$HTML"
check "the Nextflow HTML report has the class 4 and 5 SV card" grep -q 'Structural Variants, Class 4 and 5 (AnnotSV)' "$HTML"
check "the Nextflow HTML report names the gnomAD chrM file" grep -q 'gnomad\.genomes\.v3\.1\.sites\.chrM\.vcf\.bgz' "$HTML"
check "the Nextflow HTML report has the CYP2D6 depth-check line" grep -Eq 'CYP2D6 depth check: (ok|unreliable);' "$HTML"
check "the Nextflow HTML report has the cDNA read check column" grep -q '<th>cDNA read check</th>' "$HTML"
check "the Nextflow HTML report has the pairs by phase" grep -q 'Pairs in trans / cis / unknown' "$HTML"

# --- 4. main.nf passes every table on (stub run of every tool) ----------------------------
# PAIRS: each new table and the step output that says the step ran for the sample
PAIRS="cyp2d6_depth_check.tsv:pypgx_summary.tsv compound_het_pairs.tsv:prioritized.vcf.gz \
sv_class4_5_checked.tsv:sv_consensus.vcf.gz rfc1_motif.tsv:eh_stranger.vcf chrM_annotated.tsv:chrM_filtered.vcf.gz"
# stub_staged DIR: one line per pair, "<new> <HTML_REPORT tasks with the step output> <of them with the table>",
# then "VEP <VEP tasks> <of them reading VEP_NORM's file>"
stub_staged() {
  local d id p new old hits got
  for p in $PAIRS; do
    new=${p%%:*} old=${p##*:} hits=0 got=0
    while IFS= read -r d; do
      d=$(dirname "$d")
      id=$(sed -n "s/^### name: [\"'].*HTML_REPORT (\(.*\))[\"']\$/\1/p" "${d}/.command.run")
      [ -e "${d}/input/${id}_${old}" ] || continue
      hits=$((hits + 1))
      [ ! -e "${d}/input/${id}_${new}" ] || got=$((got + 1))
    done < <(grep -l "^### name: [\"'].*HTML_REPORT (" -r "$1/work" --include=.command.run 2>/dev/null)
    echo "$new $hits $got"
  done
  hits=0 got=0
  while IFS= read -r d; do
    d=$(dirname "$d")
    id=$(sed -n "s/^### name: [\"'].*:VEP (\(.*\))[\"']\$/\1/p" "${d}/.command.run")
    hits=$((hits + 1))
    [ ! -e "${d}/${id}.norm.vcf.gz" ] || got=$((got + 1))
  done < <(grep -l "^### name: [\"'].*:VEP (" -r "$1/work" --include=.command.run 2>/dev/null)
  echo "VEP $hits $got"
}
stub_run() {
  local d="${CASE_TMP}/$1"
  echo "+ nextflow run main.nf -profile test_all -stub ($1)"
  (cd "$d" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${d}/out" -work-dir "${d}/work") \
    > "${d}.log" 2>&1
  check_eq "stub run ($1) exits 0" "$?" 0
}
for n in stub stub-control; do
  rm -rf "${CASE_TMP:?}/${n}"
  mkdir -p "${CASE_TMP}/${n}"
  cp -R "${REPO}/main.nf" "${REPO}/nextflow.config" "${REPO}/nextflow_schema.json" "${REPO}/versions.env" \
    "${REPO}/conf" "${REPO}/modules" "${REPO}/workflows" "${REPO}/bin" "${REPO}/assets" "${REPO}/scripts" "${CASE_TMP}/${n}/"
  mkdir -p "${CASE_TMP}/${n}/tests/fixtures"  # test_all's consumer_sites
  cp -R "${REPO}/tests/fixtures/consumer_export" "${CASE_TMP}/${n}/tests/fixtures/"
done
JOINS='PGX.out.depth_check.map|ANNOTATION.out.comphet_pairs.map|SV.out.annotsv_class4_5.map|BAM_ANALYSIS.out.rfc1_motif.map|BAM_ANALYSIS.out.mito_annotated.map'
check_eq "main.nf joins the five tables into the report inputs" "$(grep -cE "$JOINS" "${CASE_TMP}/stub/main.nf")" 5
sed -i -E "/${JOINS}/d" "${CASE_TMP}/stub-control/main.nf"
sed -i 's/^            VEP_NORM\.out\.vcf,$/            ch_current_vcf,/' "${CASE_TMP}/stub-control/workflows/annotation.nf"
check_eq "control copy: the joins are gone" "$(grep -cE "$JOINS" "${CASE_TMP}/stub-control/main.nf")" 0
check_eq "the annotation workflow gives VEP VEP_NORM's output" "$(grep -c '^            VEP_NORM\.out\.vcf,$' "${CASE_TMP}/stub/workflows/annotation.nf")" 1
check_eq "control copy: VEP reads the VCF itself" "$(grep -c '^            VEP_NORM\.out\.vcf,$' "${CASE_TMP}/stub-control/workflows/annotation.nf")" 0
stub_run stub
stub_run stub-control
while read -r new hits got; do
  check_ge "stub run: tasks with the step before ${new}" "$hits" 1
  check_eq "stub run: each of them gets ${new}" "$got" "$hits"
done < <(stub_staged "${CASE_TMP}/stub")
while read -r new hits got; do
  check_ge "control: tasks with the step before ${new}" "$hits" 1
  check_eq "control: none gets ${new}" "$got" 0
done < <(stub_staged "${CASE_TMP}/stub-control")

[ -z "$GN_PLACED" ] || rm -f "$GN"
rm -rf "$SD" "${CASE_TMP:?}/stub" "${CASE_TMP:?}/stub-control" "${CASE_TMP:?}"/work-* "$GD"
finish
