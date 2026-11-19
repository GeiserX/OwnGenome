#!/usr/bin/env bash
# Step 01 (scripts/01-ora-to-fastq.sh) fails when orad exits 0 but writes no
# FASTQ, and keeps nothing it wrote then. A fake orad on ORAD stands in for
# Illumina's: FAKE_ORAD=good writes <name>.fastq.gz with one record (orad keeps
# the ORA file's name), none writes nothing, empty writes a gzip with no
# record, zero writes a 0-byte file, crash exits 3 after a partial write.
#
# The case that matters: step 01 runs once for R1 and once for R2 into the
# same fastq/ folder, so the R1 file from the first call must not make an R2
# call that wrote nothing look like a success (before the fix it printed a
# WARNING only when fastq/ held no .fastq.gz at all, and exited 0).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
STEP=${STEP_SCRIPT:-${REPO}/scripts/01-ora-to-fastq.sh}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

cat > "${WORK}/orad" <<'FAKE'
#!/usr/bin/env bash
# fake orad: --ora-reference DIR --path OUT FILE
set -euo pipefail
out="" file=""
while [ $# -gt 0 ]; do
  case "$1" in
    --ora-reference) shift 2 ;;
    --path) out=$2; shift 2 ;;
    *) file=$1; shift ;;
  esac
done
name=$(basename "${file%.ora}").gz
case "${FAKE_ORAD:-good}" in
  good) printf '@r1\nACGT\n+\nIIII\n' | gzip -c > "${out}/${name}" ;;
  none) : ;;
  empty) gzip -c < /dev/null > "${out}/${name}" ;;
  zero) : > "${out}/${name}" ;;
  crash) printf '@r1\nAC' | gzip -c > "${out}/${name}"; exit 3 ;;
esac
FAKE
chmod +x "${WORK}/orad"

export GENOME_DIR="${WORK}/genome" ORAD="${WORK}/orad"
mkdir -p "${GENOME_DIR}" "${WORK}/oraref" "${WORK}/in"
for r in R1 R2; do printf 'ora\n' > "${WORK}/in/s1_S1_L001_${r}_001.fastq.ora"; done
FQ="${GENOME_DIR}/s1/fastq"

# step RC MODE READ: run step 01 on READ's ORA file with FAKE_ORAD=MODE; true
# when it exited RC. Its output is in ${WORK}/out.
step() {
  local rc=0
  FAKE_ORAD=$2 bash "$STEP" s1 "${WORK}/oraref" "${WORK}/in/s1_S1_L001_${3}_001.fastq.ora" > "${WORK}/out" 2>&1 || rc=$?
  [ "$rc" -eq "$1" ] || { echo "  exit ${rc}, expected $1:"; sed 's/^/    /' "${WORK}/out"; return 1; }
}
leftovers() { find "$FQ" -mindepth 1 -name '.orad.*' 2>/dev/null; }

# R1 decompresses: the file lands in fastq/ under the ORA file's name.
if step 0 good R1 && [ -s "${FQ}/s1_S1_L001_R1_001.fastq.gz" ] && [ -z "$(leftovers)" ]; then
  pass "a good orad run puts s1_S1_L001_R1_001.fastq.gz in fastq/ and leaves no temporary folder"
else fail "a good orad run did not leave exactly the R1 FASTQ in fastq/"; fi

# R2 writes nothing, with the R1 file already in fastq/: the step fails.
if step 1 none R2 && grep -q 'orad exited 0 but wrote no .fastq.gz' "${WORK}/out" \
   && [ ! -e "${FQ}/s1_S1_L001_R2_001.fastq.gz" ] && [ -z "$(leftovers)" ]; then
  pass "an orad that writes nothing fails the step, even with R1 already in fastq/"
else fail "an orad that wrote nothing for R2 did not fail the step"; fi

# A gzip with no record, and a 0-byte file: fail, and neither reaches fastq/.
for m in empty zero; do
  if step 1 "$m" R2 && grep -q 'holds no FASTQ record' "${WORK}/out" \
     && [ ! -e "${FQ}/s1_S1_L001_R2_001.fastq.gz" ] && [ -z "$(leftovers)" ]; then
    pass "FAKE_ORAD=${m}: a .fastq.gz with no record fails the step and is not kept"
  else fail "FAKE_ORAD=${m}: a .fastq.gz with no record did not fail the step, or was kept"; fi
done

# orad fails after a partial write: the step fails and the partial file goes.
if step 3 crash R2 && [ ! -e "${FQ}/s1_S1_L001_R2_001.fastq.gz" ] && [ -z "$(leftovers)" ]; then
  pass "a failing orad fails the step with its exit code and leaves no partial file"
else fail "a failing orad left a partial file or did not fail the step"; fi

# R2 then decompresses: both files are in fastq/, R1 untouched.
if step 0 good R2 && [ -s "${FQ}/s1_S1_L001_R2_001.fastq.gz" ] && [ -s "${FQ}/s1_S1_L001_R1_001.fastq.gz" ]; then
  pass "a good R2 run adds its file next to R1"
else fail "a good R2 run did not add its file"; fi

if [ "$FAILS" -gt 0 ]; then echo "${FAILS} check(s) failed"; exit 1; fi
echo "All step 01 orad checks passed."
