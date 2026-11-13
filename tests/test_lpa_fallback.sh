#!/usr/bin/env bash
# Step 46 [3/3]: a KIV-2 4925G>A reader that crashes never costs the KIV-2
# count (one og cn call reads both; it writes the count table first). scripts/46-lpa.sh then writes lpa/<id>_lpa_kiv2_4925.tsv saying
# 'not assessed', exits 0 and keeps lpa/<id>_lpa_kiv2.tsv; the LPA Nextflow
# process (modules/local/lpa/main.nf) writes the same fallback table, byte for
# byte. The e2e case genotype-rules-2-lpa-null runs the module's block with a
# crashing reader in its image.
#
# A fake `docker` on PATH stands in for the containers: the og cn call writes
# the count table, then exits 1 at --marker-out (or 0 for the control,
# writing a 'present' table), everything else succeeds.
#
# Controls: with the reader working the table is the reader's ('present'), so
# the not-assessed check can fail; with the fallback removed from a copy of
# the step, the step exits non-zero, so the exit check can fail.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

mkdir -p "${WORK}/bin"
cat > "${WORK}/bin/docker" <<'FAKE'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$FAKE_DOCKER_LOG"
host() { printf '%s' "${GENOME_DIR}${1#/genome}"; }
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
  case "${args[i]}" in
    --table) printf 'metric\tvalue\nstatus\testimated\nkiv2_copies\t38.0\n' > "$(host "${args[i + 1]}")" ;;
    --marker-out)
      if [ "${FAKE_READER:-crash}" = crash ]; then echo "planted crash" >&2; exit 1; fi
      printf 'metric\tvalue\nstatus\tpresent\n' > "$(host "${args[i + 1]}")" ;;
  esac
done
exit 0
FAKE
chmod +x "${WORK}/bin/docker"

GD="${WORK}/genome"
mkdir -p "${GD}/S1/aligned" "${GD}/reference"
: > "${GD}/S1/aligned/S1_sorted.bam"
: > "${GD}/S1/aligned/S1_sorted.bam.bai"
: > "${GD}/reference/GRCh38_no_alt_analysis_set.fasta"
: > "${GD}/reference/GRCh38_no_alt_analysis_set.fasta.fai"
LOG="${WORK}/docker.log"
NULL="${GD}/S1/lpa/S1_lpa_kiv2_4925.tsv"
COUNT="${GD}/S1/lpa/S1_lpa_kiv2.tsv"

# run_step <script> <reader: crash|works>; sets OUT and RC
run_step() {
  : > "$LOG"
  set +e
  OUT=$(PATH="${WORK}/bin:${PATH}" GENOME_DIR="$GD" FAKE_DOCKER_LOG="$LOG" FAKE_READER="$2" bash "$1" S1 2>&1)
  RC=$?
  set -e
}
status() { awk -F'\t' '$1 == "status" {print $2}' "$NULL" 2>/dev/null; }

run_step "${REPO}/scripts/46-lpa.sh" crash
if [ "$RC" -eq 0 ]; then pass "a crashing 4925G>A reader: step 46 exits 0"
else fail "a crashing 4925G>A reader: step 46 exited ${RC}: $(tail -3 <<<"$OUT" | tr '\n' '|')"; fi
if grep -q -- '--marker LPA_KIV2_4925' "$LOG"; then pass "the 4925G>A reader was called"
else fail "the 4925G>A reader was never called"; fi
if [ "$(status)" = "not assessed" ] && grep -q "not assessed, the check failed; see the log of step 46" "$NULL"; then
  pass "a crashing 4925G>A reader: the table says not assessed, the check failed"
else fail "a crashing 4925G>A reader: the table is not the not-assessed one: $(cat "$NULL" 2>/dev/null)"; fi
if [ "$(awk -F'\t' '$1 == "kiv2_copies" {print $2}' "$COUNT" 2>/dev/null)" = "38.0" ]; then
  pass "a crashing 4925G>A reader: the KIV-2 count table stays"
else fail "a crashing 4925G>A reader: the KIV-2 count table is gone: $(cat "$COUNT" 2>/dev/null)"; fi
if grep -q "WARNING: the KIV-2 4925G>A check failed" <<<"$OUT"; then pass "a crashing 4925G>A reader: the log warns"
else fail "a crashing 4925G>A reader: no warning in the log"; fi
STEP_TABLE=$(cat "$NULL" 2>/dev/null || true)

# The module's fallback table: its printf line, cut from main.nf, as the task's bash runs it.
MOD="${REPO}/modules/local/lpa/main.nf"
LINES=$(sed -n "/^ *printf 'metric.*LPA KIV-2 4925G>A/,/> \${meta.id}_lpa_kiv2_4925.tsv\$/p" "$MOD")
if [ "$(wc -l <<<"$LINES" | tr -d ' ')" -eq 2 ]; then
  # shellcheck disable=SC2016  # ${meta.id} is Groovy text to replace, not a shell variable
  MODULE_TABLE=$(cd "$WORK" && sed -e 's/\\\\/\\/g' -e 's/\${meta.id}/S1/g' <<<"$LINES" | bash && cat S1_lpa_kiv2_4925.tsv)
  if [ -n "$STEP_TABLE" ] && [ "$MODULE_TABLE" = "$STEP_TABLE" ]; then
    pass "the LPA process writes step 46's fallback table, byte for byte"
  else fail "the LPA process's fallback table differs from step 46's: $(diff <(echo "$STEP_TABLE") <(echo "$MODULE_TABLE") | tr '\n' '|')"; fi
else fail "the module's fallback printf was not found in ${MOD}: ${LINES}"; fi

# Control 1: a working reader writes its own table; the not-assessed check would fail on it.
run_step "${REPO}/scripts/46-lpa.sh" works
if [ "$RC" -eq 0 ] && [ "$(status)" = "present" ]; then pass "control: a working reader's table is kept (present)"
else fail "control: a working reader gave exit ${RC}, status '$(status)'"; fi

# Control 2: step 46 with the fallback cut out stops on the crash, so the exit check can fail.
mkdir -p "${WORK}/scripts/lib"
cp "${REPO}/scripts/lib/"* "${WORK}/scripts/lib/"
cp "${REPO}/versions.env" "${WORK}/" 2>/dev/null || true
# shellcheck disable=SC2016  # the $(cpath ...) is text of the step, not to expand here
sed -e 's/--marker-out "$(cpath "$NULL")" || rc=\$?$/--marker-out "$(cpath "$NULL")"/' \
  "${REPO}/scripts/46-lpa.sh" > "${WORK}/scripts/46-lpa.sh"
bash -n "${WORK}/scripts/46-lpa.sh" || fail "control: the copy of step 46 without the fallback is not valid bash"
if cmp -s "${REPO}/scripts/46-lpa.sh" "${WORK}/scripts/46-lpa.sh"; then
  fail "control: the fallback could not be cut from the copy of step 46"
else
  run_step "${WORK}/scripts/46-lpa.sh" crash
  if [ "$RC" -ne 0 ] && grep -q "planted crash" <<<"$OUT"; then pass "control: without the fallback a crashing reader stops step 46 (exit ${RC})"
  else fail "control: without the fallback step 46 exited ${RC}, not on the planted crash: $(tail -3 <<<"$OUT" | tr '\n' '|')"; fi
fi

if [ "$FAILS" -ne 0 ]; then echo "RESULT: ${FAILS} FAILED"; exit 1; fi
echo "RESULT: ALL PASS"
