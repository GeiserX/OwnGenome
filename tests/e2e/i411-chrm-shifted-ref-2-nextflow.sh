#!/usr/bin/env bash
# The shifted-reference pass on the Nextflow path: MITO_VARIANTS on the
# fixture BAM, under the docker profile's own settings, so its tasks start as
# this runner's uid (-u of docker.runOptions) and the module's --user 0:0
# lets BwaSpark's Hadoop login find a name. A pass that failed would leave
# the task green, with a WARNING in its .command.err and the linear calls
# alone, so the traces it leaves are checked here, and the calls are
# compared with step 20's on the same BAM (base case 36).
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G="$GENOME_DIR"
T="$CASE_TMP"
NP="${T}/project"
OUT="${T}/nf/results"
WORK="${T}/nf/work"
rm -rf "$NP" "${T}/nf"
mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
BAM="${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
REF="${G}/reference/GRCh38_no_alt_analysis_set"
cat > "${NP}/main.nf" <<NF
include { MITO_VARIANTS } from '${REPO}/modules/local/mito_variants/main'

workflow {
    MITO_VARIANTS(channel.of([[id: '${SAMPLE}'], file('${BAM}'), file('${BAM}.bai'), [], []]),
                  file('${REF}.fasta'), file('${REF}.fasta.fai'), file('${REF}.dict'))
}
NF
# The docker profile of nextflow.config, as it is there.
cat > "${NP}/nextflow.config" <<'NFCONF'
docker.enabled = true
docker.runOptions = '-u ${UID}:${GROUPS} -e HOME=/tmp'
docker.fixOwnership = true
process.containerOptions = '--network none'
process.cpus = 4
process.memory = '8 GB'
NFCONF
cat >> "${NP}/nextflow.config" <<NFCONF
params.outdir = '${OUT}'
params.publish_dir_mode = 'copy'
NFCONF
echo "+ nextflow run (MITO_VARIANTS) as uid $(id -u)"
(cd "$T" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/conf/containers.config" -work-dir "$WORK") \
  > "${T}/nf.log" 2>&1
NF_RC=$?
check_eq "the Nextflow run exits 0" "$NF_RC" 0
[ "$NF_RC" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 12 "${T}/nf.log" | head -60

VCF="${OUT}/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz"
check "MITO_VARIANTS published the filtered calls" test -s "$VCF"
check "they carry the shifted-pass header line" \
  has '^##og_shifted_pass=chrM:16025-16569 and chrM:1-575 from Mutect2 on chrM rotated by 8000 bp' \
  "$(gzip -dc "$VCF" 2>/dev/null | grep '^##')"
# shellcheck disable=SC2016  # awk program
check_ge "control-region calls in the filtered calls" \
  "$(gzip -dc "$VCF" 2>/dev/null | awk '!/^#/ && ($2 <= 575 || $2 >= 16025)' | wc -l | tr -d ' ')" 1

# The task's log: the reads were kept as pairs, the calls mapped back, and no
# part of the pass failed.
RUN=$(grep -lE '^### name: "([^"]*:)?MITO_VARIANTS \(' "$WORK"/*/*/.command.run 2>/dev/null | head -n 1)
ERR=$(cat "${RUN%.command.run}.command.err" 2>/dev/null)
grep -E 'shifted-reference pass|og vcf shift-back|LoginException' <<< "$ERR"
check "the task's .command.err is found and read" test -n "$ERR"
check "the task kept the chrM reads as pairs" has 'shifted-reference pass: [0-9]+ of [0-9]+ chrM reads kept as complete pairs' "$ERR"
check "the task mapped the rotated calls back" has '^og vcf shift-back: [0-9]+ control-region calls from the rotated pass' "$ERR"
check "the task did not fall back to the linear calls" lacks 'shifted-reference pass failed' "$ERR"
check "BwaSpark found a user name" lacks 'LoginException' "$ERR"
check_eq "files in the work directory left to another user (docker.fixOwnership)" \
  "$(find "$WORK" ! -user "$(id -u)" 2>/dev/null | wc -l | tr -d ' ')" 0

# The same calls as step 20 on the same BAM (base case 36): FilterMutectCalls
# keeps every record, so position, REF and ALT are compared whatever the
# NuMT depth of either run.
calls() { gzip -dc "$1" 2>/dev/null | awk -F'\t' '!/^#/ {print $2, $4, $5}'; }
BASH_CALLS=$(calls "${G}/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz")
check_ge "step 20's calls are there (case 36)" "$(grep -c . <<< "$BASH_CALLS")" 1
NF_CALLS=$(calls "$VCF")
echo "calls: MITO_VARIANTS $(grep -c . <<< "$NF_CALLS"), step 20 $(grep -c . <<< "$BASH_CALLS")"
diff <(echo "$BASH_CALLS") <(echo "$NF_CALLS") | head -20
check_eq "MITO_VARIANTS and step 20 give the same calls (sha256 of POS REF ALT)" \
  "$(sha256sum <<< "$NF_CALLS" | cut -c1-16)" "$(sha256sum <<< "$BASH_CALLS" | cut -c1-16)"

finish
