#!/usr/bin/env bash
# A replaced raw FASTQ is trimmed and aligned again; the old trimmed pair is
# never reused. Step 01b writes fastq_trimmed/<sample>_fastp_inputs.tsv, the
# input_stamp of the raw R1 and R2 it trimmed (name, size, mtime, sha256 of
# the first and last MiB), and skips only while that stamp still matches.
# Step 02 aligns fastq_trimmed/ only while it matches too, and otherwise stops
# before it aligns, saying to run 01b again or to set FASTQ_SUBDIR=fastq.
# Before this, 01b skipped on the trimmed files alone and 02 aligned them: a
# BAM of the old reads.
# The fake fastp copies its input to its output, so the trimmed pair holds the
# raw reads it was made from; the fake samtools answers `view -c` with
# FAKE_PRIMARY, which step 02 compares with the records of the reads it read.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
S="${GENOME_DIR}/s1"
RAW="${S}/fastq" TRIM="${S}/fastq_trimmed"
STAMP="${TRIM}/s1_fastp_inputs.tsv"
BAM="${S}/aligned/s1_sorted.bam"
T="${S}/qc/s1_read_count.tsv"
# fastq N BASES: N four-line records of BASES
fastq() { local i; for ((i = 0; i < $1; i++)); do printf '@r%d\n%s\n+\nIIII\n' "$i" "$2"; done; }
mkdir -p "$RAW"
for r in R1 R2; do fastq 3 ACGT | gzip -c > "${RAW}/s1_${r}.fastq.gz"; done

use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
args=" $* "
last=${!#}
word_after() {
  local -a w
  read -r -a w <<<"$args"
  local i
  for ((i = 0; i < ${#w[@]} - 1; i++)); do
    if [ "${w[i]}" = "$1" ]; then printf '%s' "${w[i + 1]}"; return 0; fi
  done
  return 1
}
case "$args" in
  *" fastp "*)
    cp "$(host_path "$(word_after -i)")" "$(host_path "$(word_after -o)")"
    cp "$(host_path "$(word_after -I)")" "$(host_path "$(word_after -O)")"
    printf '{"summary":{"after_filtering":{"total_reads":0}}}\n' > "$(host_path "$(word_after -j)")"
    : > "$(host_path "$(word_after -h)")" ;;
  *" minimap2 "*" -d "*) : > "$(host_path "$(word_after -d)")" ;;
  *" minimap2 "*) printf '@HD\tVN:1.6\tSO:unsorted\n' ;;
  *markdup*)
    cat > /dev/null
    out=$(host_path "$last")
    printf 'BAI\001\n' > "${out}.bai"
    printf 'BAM\001 run %s\n' "${FAKE_TAG:-1}" > "$out" ;;
  *" samtools view -c "*) echo "${FAKE_PRIMARY:-6}" ;;
  *" samtools quickcheck "*) : ;;
  *) exec "${CASE_WORK}/hook-outputs" image "$@" ;;
esac
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

val() { awk -F'\t' -v k="$1" '$1 == k { print $2; exit }' "$2"; }
fastp_ran() { grep -q ' fastp ' "$FAKE_DOCKER_LOG"; }
aligned() { grep -q 'minimap2 -t' "$FAKE_DOCKER_LOG"; }

# --- first run: 01b trims and records its inputs, 02 aligns the trimmed pair --------
: > "$FAKE_DOCKER_LOG"
run_expect 0 trim-1 "${SCRIPTS}/01b-fastp-qc.sh" s1
fastp_ran || fail "step 01b did not run fastp on the first run"
[ -s "$STAMP" ] || fail "step 01b wrote no ${STAMP}"
[ "$(wc -l < "$STAMP" | tr -d ' ')" = 2 ] || fail "${STAMP} does not hold one line for R1 and one for R2"
[ "$(cut -f1 "$STAMP" | paste -sd, -)" = "s1_R1.fastq.gz,s1_R2.fastq.gz" ] || fail "${STAMP} names '$(cut -f1 "$STAMP" | paste -sd, -)'"
run_expect 0 align-1 env FAKE_PRIMARY=6 "${SCRIPTS}/02-alignment.sh" s1
output_has align-1 'Using trimmed FASTQs from fastp'
[ "$(val input_source "$T")" = "fastq_trimmed/s1_R1.fastq.gz 3 records + fastq_trimmed/s1_R2.fastq.gz 3 records" ] \
  || fail "first alignment: input_source is '$(val input_source "$T")'"

# --- nothing changed: 01b skips ------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 trim-same "${SCRIPTS}/01b-fastp-qc.sh" s1
output_has trim-same 'Trimmed FASTQs of these raw reads already exist'
if fastp_ran; then fail "step 01b trimmed again although the raw reads did not change"; fi

# --- the raw pair is replaced (5 records each): 02 refuses the old trimmed pair ------
for r in R1 R2; do fastq 5 ACGT | gzip -c > "${RAW}/s1_${r}.fastq.gz"; done
bam_before=$(cat "$BAM")
: > "$FAKE_DOCKER_LOG"
run_expect 1 align-stale env FAKE_PRIMARY=6 FAKE_TAG=stale "${SCRIPTS}/02-alignment.sh" s1
output_has align-stale 'fastq_trimmed/ was not trimmed from the raw reads in fastq/ now'
output_has align-stale 'changed since step 01b trimmed them'
output_has align-stale 'Run scripts/01b-fastp-qc.sh s1 first'
if aligned; then fail "step 02 aligned the trimmed pair of the old raw reads"; fi
[ "$(cat "$BAM")" = "$bam_before" ] || fail "step 02 replaced the BAM from the old trimmed pair"

# --- 01b trims again, then 02 aligns the new reads ------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 trim-new "${SCRIPTS}/01b-fastp-qc.sh" s1
output_has trim-new 'The raw reads changed since .* was trimmed from them: trimming again'
fastp_ran || fail "step 01b reused the trimmed pair of the old raw reads"
cmp -s <(gzip -cd "${TRIM}/s1_R1.fastq.gz") <(gzip -cd "${RAW}/s1_R1.fastq.gz") \
  || fail "fastq_trimmed/ does not hold the new reads after step 01b"
: > "$FAKE_DOCKER_LOG"
run_expect 0 align-new env FAKE_PRIMARY=10 FAKE_TAG=new "${SCRIPTS}/02-alignment.sh" s1
aligned || fail "step 02 did not align after step 01b trimmed the new reads"
[ "$(val input_reads "$T")" = 10 ] || fail "second alignment: input_reads is '$(val input_reads "$T")', want 10 (5 + 5)"
grep -q 'run new' "$BAM" || fail "the BAM is not the new alignment"

# --- same size, same time, other bases: the checksum still tells them apart ----------
cp -p "${RAW}/s1_R1.fastq.gz" "${CASE_WORK}/old_R1.fastq.gz"
fastq 5 TGCA | gzip -c > "${RAW}/s1_R1.fastq.gz"
touch -r "${CASE_WORK}/old_R1.fastq.gz" "${RAW}/s1_R1.fastq.gz"
[ "$(wc -c < "${RAW}/s1_R1.fastq.gz")" = "$(wc -c < "${CASE_WORK}/old_R1.fastq.gz")" ] \
  || fail "test setup: the replacement R1 should have the old file's size"
: > "$FAKE_DOCKER_LOG"
run_expect 0 trim-same-size "${SCRIPTS}/01b-fastp-qc.sh" s1
fastp_ran || fail "step 01b skipped although R1 has other content with the same size and time"

# --- a trimmed pair with no record (trimmed before the record existed) ----------------
rm -f "$STAMP"
run_expect 1 align-unstamped "${SCRIPTS}/02-alignment.sh" s1
output_has align-unstamped 'it has no s1_fastp_inputs.tsv'
: > "$FAKE_DOCKER_LOG"
run_expect 0 align-raw env FASTQ_SUBDIR=fastq FAKE_PRIMARY=10 "${SCRIPTS}/02-alignment.sh" s1
output_has align-raw 'Using explicit FASTQ_SUBDIR=fastq'
aligned || fail "FASTQ_SUBDIR=fastq did not align the raw reads"
: > "$FAKE_DOCKER_LOG"
run_expect 0 trim-unstamped "${SCRIPTS}/01b-fastp-qc.sh" s1
output_has trim-unstamped 'has no record of the raw reads it was trimmed from'
fastp_ran || fail "step 01b trusted a trimmed pair with no record of its inputs"
[ -s "$STAMP" ] || fail "step 01b wrote no record after trimming again"
