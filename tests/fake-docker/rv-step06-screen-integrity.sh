#!/usr/bin/env bash
# Step 06 (scripts/06-clinvar-screen.sh) never publishes a partial or
# mismatched ClinVar screen:
#   - a failed first bcftools of a nested pipeline (view | norm when ClinVar
#     is normalised, view | norm when the sample is filtered, annotate | view
#     when the hits are written) fails the step and leaves no hits file. The
#     container's sh has no pipefail, so each pipeline must run in bash with
#     pipefail;
#   - two samples screened at the same time against a normalised ClinVar
#     older than its source build it once, and the published VCF and its
#     .tbi come from that one write (with flock, which every Linux runner
#     has; a mkdir lock is not exclusive with the uutils mkdir of Ubuntu
#     26.04);
#   - without flock the step takes lock_acquire's folder and removes it;
#   - a stale copy in a clinvar/ folder the user cannot write stops the step.
# The docker run hook runs each container command for real, with /genome
# mapped to GENOME_DIR and a fake bcftools and og first on PATH. The fake
# tools write plain text and stamp every write with an id, which an index
# copies, so a VCF and a .tbi of two different writes can be told apart.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
C="${GENOME_DIR}/clinvar"
NORM="${C}/clinvar_pathogenic_chr.norm.vcf.gz"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
seed_sample "$GENOME_DIR" sample2

TOOLS="${CASE_WORK}/tools"
mkdir -p "$TOOLS"
export FAKE_TOOL_LOG="${CASE_WORK}/tools.log"

# bcftools: writes a header and a write id to -o (or stdout), indexes by
# copying the data file's id into FILE.tbi, and fails like the real one on a
# missing input. FAKE_FAIL (a regex over the arguments) makes the matching
# call write a header and exit 1. FAKE_NORM_SLEEP slows the norm that
# writes the normalised ClinVar, to hold two runs in that window together.
cat > "${TOOLS}/bcftools" <<'TOOL'
#!/usr/bin/env bash
set -uo pipefail
args="$*"
printf '%s\n' "bcftools ${args}" >> "$FAKE_TOOL_LOG"
if [ -n "${FAKE_FAIL:-}" ] && [[ "$args" =~ $FAKE_FAIL ]]; then
  printf '##fileformat=VCFv4.2\n'
  echo "fake bcftools: failing on purpose: ${args}" >&2
  exit 1
fi
sub=$1
shift
out="" last="" h=false n=false s=false
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out=$2; shift 2; continue ;;
    -h) h=true ;;
    -n) n=true ;;
    -s) s=true ;;
  esac
  last=$1
  shift
done
# No file after the options (or -): bcftools reads stdin.
case "$last" in -*|"") last=-; cat > /dev/null ;; esac
if [ "$last" != - ] && [ ! -e "$last" ]; then echo "fake bcftools: ${last}: No such file" >&2; exit 1; fi
case "$sub" in
  index)
    if $s; then printf 'chr1\t248956422\t3\n'
    elif $n; then echo 3
    else grep -m1 '^##write=' "$last" > "${last}.tbi" || { echo "fake bcftools: ${last} has no write id" >&2; exit 1; }
    fi
    exit 0 ;;
  view)
    if $h; then printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n'; exit 0; fi
    case " ${args} " in *" -H "*) printf 'chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n'; exit 0 ;; esac ;;
  query) [ -z "$out" ] || : > "$out"; exit 0 ;;
esac
if [ "$sub" = norm ] && [[ "$out" == *clinvar_pathogenic_chr.norm* ]]; then
  /bin/sleep "${FAKE_NORM_SLEEP:-0}"
fi
body=$(printf '##fileformat=VCFv4.2\n##write=%s.%s.%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' "$sub" "$$" "$RANDOM")
if [ -n "$out" ]; then printf '%s\n' "$body" > "$out"; else printf '%s\n' "$body"; fi
TOOL

# og: the subcommands step 06 calls, writing what the step reads next.
cat > "${TOOLS}/og" <<'TOOL'
#!/usr/bin/env bash
set -uo pipefail
printf '%s\n' "og $*" >> "$FAKE_TOOL_LOG"
case "$1 $2" in
  "clinvar join-adjacent")
    if [ "$4" = tag ]; then cat "$5"; else cat; fi ;;
  "clinvar deep-intronic")
    cat > /dev/null
    printf 'chrom\tpos\tref\talt\tgenotype\n' ;;
  "cdna check"|"cdna scan")
    while [ $# -gt 0 ]; do
      if [ "$1" = --out ]; then printf 'header\n' > "$2"; fi
      shift
    done ;;
  *) echo "fake og: no answer for: $*" >&2; exit 1 ;;
esac
TOOL
chmod +x "${TOOLS}/bcftools" "${TOOLS}/og"

# The docker run hook: each container path under /genome becomes the host
# path, /pgp-bin/og the fake og; the command runs from the sample's clinvar
# folder (step 06 runs og cdna with -w there and a relative --out).
cat > "${CASE_WORK}/hook-run" <<'HOOK'
#!/usr/bin/env bash
set -uo pipefail
shift   # the image
cmd=()
for a in "$@"; do
  a=${a//\/pgp-bin\/og/og}
  cmd+=("${a//\/genome\//${GENOME_DIR}/}")
done
while IFS= read -r v; do
  case "$v" in
    *:/genome/clinvar*) ;;
    *:/genome/*) [ -d "${v%%:*}/clinvar" ] && cd "${v%%:*}/clinvar" ;;
  esac
done <<<"${FAKE_DOCKER_VOLUMES:-}"
PATH="${FAKE_TOOLS}:${PATH}" exec "${cmd[@]}"
HOOK
chmod +x "${CASE_WORK}/hook-run"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-run" FAKE_TOOLS="$TOOLS"

write_id() { grep -m1 '^##write=' "$1" 2>/dev/null || echo none; }

# --- A full screen with the fake tools succeeds -------------------------------------
# The positive control: the cases below fail only because of the failure they plant.
run_rc clean "${SCRIPTS}/06-clinvar-screen.sh" sample1
expect_rc clean 0
[ -f "${GENOME_DIR}/sample1/clinvar/sample1_clinvar_hits.vcf" ] || fail "a clean screen wrote no hits file"
[ "$(write_id "$NORM")" = "$(write_id "${NORM}.tbi")" ] || fail "a clean screen published a .tbi of another write than the normalised ClinVar"
[ ! -e "${C}/.norm.lock.d" ] || fail "a clean screen left the lock folder clinvar/.norm.lock.d behind"

# --- A failed first stage of each nested pipeline fails the step ---------------------
# name|regex of the first stage's arguments
for stage in \
  "clinvar-view|^view -T .*/clinvar/clinvar_pathogenic_chr\.vcf\.gz$" \
  "sample-view|^view -f .*/sample1\.vcf\.gz$" \
  "hits-annotate|^annotate "; do
  name=${stage%%|*}
  rm -rf "${GENOME_DIR}/sample1/clinvar" "$NORM" "${NORM}.tbi"
  : > "$FAKE_TOOL_LOG"
  FAKE_FAIL=${stage#*|} run_rc "fail-${name}" "${SCRIPTS}/06-clinvar-screen.sh" sample1
  grep -q 'failing on purpose' "${CASE_WORK}/fail-${name}.out" || fail "${name}: the planted failure never ran"
  [ "$RC" -ne 0 ] || fail "${name}: step 06 exited 0 after its first bcftools failed (a partial screen)"
  for f in sample1_clinvar_hits.vcf sample1_clinvar_hits.tsv; do
    [ ! -e "${GENOME_DIR}/sample1/clinvar/${f}" ] || fail "${name}: step 06 left ${f} after a failed bcftools"
  done
  [ ! -e "${C}/.norm.lock.d" ] || fail "${name}: a failed screen left the lock folder clinvar/.norm.lock.d behind"
  echo "${name}: a failed first bcftools fails the step and leaves no hits: PASS"
done

# --- Two samples at once against a normalised ClinVar older than its source ----------
printf '##write=old\n' > "$NORM"
printf '##write=old\n' > "${NORM}.tbi"
touch -t 202001010000 "$NORM" "${NORM}.tbi"
rm -rf "${GENOME_DIR}"/sample[12]/clinvar
: > "$FAKE_TOOL_LOG"
export FAKE_NORM_SLEEP=2
rc1=0 rc2=0
"${SCRIPTS}/06-clinvar-screen.sh" sample1 > "${CASE_WORK}/both-1.out" 2>&1 &
p1=$!
"${SCRIPTS}/06-clinvar-screen.sh" sample2 > "${CASE_WORK}/both-2.out" 2>&1 &
p2=$!
wait "$p1" || rc1=$?
wait "$p2" || rc2=$?
unset FAKE_NORM_SLEEP
cat "${CASE_WORK}/both-1.out" "${CASE_WORK}/both-2.out"
[ "$rc1" -eq 0 ] && [ "$rc2" -eq 0 ] || fail "two screens at once exited ${rc1} and ${rc2}, want 0 and 0"
builds=$(grep -c '^bcftools norm .*clinvar_pathogenic_chr\.norm' "$FAKE_TOOL_LOG" || true)
[ "$builds" -eq 1 ] || fail "two screens at once normalised ClinVar ${builds} times, want once"
vid=$(write_id "$NORM")
tid=$(write_id "${NORM}.tbi")
echo "published VCF ${vid}, .tbi ${tid}"
[ "$vid" != '##write=old' ] && [ "$vid" != none ] || fail "the normalised ClinVar was not rebuilt (${vid})"
[ "$vid" = "$tid" ] || fail "the published .tbi (${tid}) is from another write than the VCF (${vid})"
for s in sample1 sample2; do
  [ -f "${GENOME_DIR}/${s}/clinvar/${s}_clinvar_hits.vcf" ] || fail "${s}: no hits file after two screens at once"
done
[ -z "$(find "$C" -name '*.part*')" ] || fail "two screens at once left part files: $(find "$C" -name '*.part*')"
[ ! -e "${C}/.norm.lock.d" ] || fail "two screens at once left the lock folder clinvar/.norm.lock.d behind"
echo "two screens at once build the normalised ClinVar once, VCF and .tbi of one write: PASS"

# --- Without flock: lock_acquire's folder, taken and removed ---------------------------
# One run only: the folder lock is exclusive with a POSIX mkdir (macOS, GNU),
# but not with the uutils mkdir of Ubuntu 26.04, so two runs here would test
# the runner's mkdir, not the step.
touch -t 202001010000 "$NORM" "${NORM}.tbi"
rm -rf "${GENOME_DIR}/sample1/clinvar"
: > "$FAKE_TOOL_LOG"
(
  hide_commands flock
  run_rc no-flock "${SCRIPTS}/06-clinvar-screen.sh" sample1
  expect_rc no-flock 0
)
grep -q '^bcftools norm .*clinvar_pathogenic_chr\.norm' "$FAKE_TOOL_LOG" || fail "without flock the stale normalised ClinVar was not rebuilt"
[ "$(write_id "$NORM")" = "$(write_id "${NORM}.tbi")" ] || fail "without flock the published .tbi is from another write than the VCF"
[ ! -e "${C}/.norm.lock.d" ] || fail "without flock the lock folder clinvar/.norm.lock.d was left behind"
echo "without flock the step builds under lock_acquire's folder and removes it: PASS"

# --- A stale normalised ClinVar in a folder the user cannot write ----------------------
# Without flock the lock is a folder made in clinvar/, so the step must stop
# with an error there, not wait for a lock it can never take. root writes anywhere, so
# this part needs an ordinary user, as on the CI runner.
if [ "$(id -u)" -eq 0 ]; then
  echo "read-only clinvar/: not checked as root, which can write any folder"
else
  touch -t 202001010000 "$NORM" "${NORM}.tbi"
  rm -rf "${GENOME_DIR}/sample1/clinvar"
  chmod a-w "$C"
  rc=0
  (
    hide_commands flock
    run_rc readonly-clinvar timeout 60 "${SCRIPTS}/06-clinvar-screen.sh" sample1
    exit "$RC"
  ) || rc=$?
  chmod u+w "$C"
  [ "$rc" -ne 124 ] || fail "step 06 was still waiting after 60 s for a lock in a clinvar/ folder it cannot write"
  [ "$rc" -ne 0 ] || fail "step 06 exited 0 with a stale normalised ClinVar it cannot rebuild"
  output_has readonly-clinvar 'clinvar is not writable'
  echo "a stale normalised ClinVar in a read-only folder stops the step with an error: PASS"
fi
