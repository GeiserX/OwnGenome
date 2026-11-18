#!/usr/bin/env bash
# Step 53 (Locityper) in run-all.sh's plan, its switches, and the DONE
# guards of the step script and the module. Moved from
# tests/test_locityper_report.py, cases 7 to 9:
#   7. a default run neither lists step 53 nor passes locityper to
#      Nextflow; TOOLS naming locityper with the database installed lists it
#      as running and passes --tools locityper, --locityper_data and
#      --locityper_bg_region; without the database it is skipped, "data not
#      installed". Control: the opted-in checks fail on the default run's
#      output.
#   8. LOCITYPER_DATA and LOCITYPER_BG_REGION are run-all.sh switches: the
#      command logs/run_status.tsv keeps names both, and 'meta switches'
#      lists them. Control: a run without LOCITYPER_DATA records none.
#   9. scripts/53-locityper.sh and the LOCITYPER module's guard (its script
#      lines before Locityper, run in bash) stop on a missing DONE and on a
#      DONE that names this .fai with another FASTA, and pass on this .fai
#      and FASTA. Control: a DONE that names only the .fai, which the old
#      guard accepted, is refused.
# The table og locityper report writes is tested in
# rust/crates/og-reports/tests/optin_locityper.rs.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

PANEL=$(sed -n 's/^LOCITYPER_PANEL_VERSION="\([^"]*\)".*/\1/p' "${REPO_ROOT}/versions.env")
[ -n "$PANEL" ] || fail "cannot read LOCITYPER_PANEL_VERSION from versions.env"
BG=chr20:10000001-10500000

# run_all NAME TOOLS DONE [VAR=VALUE...]: run-all.sh on placeholder inputs in
# its own genome folder; NFLINE is the nextflow call it made.
run_all() {
  local name=$1 tools=$2 done=$3 g="${CASE_WORK}/$1/g"
  shift 3
  seed_reference "$g"
  seed_sample "$g" s1
  if [ "$done" = yes ]; then
    mkdir -p "${g}/reference/locityper-${PANEL}"
    echo ok > "${g}/reference/locityper-${PANEL}/DONE"
  fi
  : > "$FAKE_DOCKER_LOG"
  run_rc "$name" env HOME="${CASE_WORK}/${name}" GENOME_DIR="$g" SKIP_VALIDATION=true LOCITYPER_BG_REGION="$BG" \
    ${tools:+"TOOLS=$tools"} "$@" bash "${SCRIPTS}/run-all.sh" s1 male
  NFLINE=$(grep '^nextflow :: .*/main\.nf ' "$FAKE_DOCKER_LOG" | tail -n 1 || true)
}

# opted_in NAME: the two checks of an opted-in run, as one line per failure.
opted_in() {
  grep -qE '^  53 Hard loci \(Locityper\) +runs$' "${CASE_WORK}/$1.out" || echo "step 53 is not listed as running"
  grep -qE -- '--tools [^ ]*locityper' <<<"$NFLINE" && grep -qF -- '--locityper_data ' <<<"$NFLINE" \
    && grep -qF -- "--locityper_bg_region ${BG}" <<<"$NFLINE" || echo "nextflow did not get locityper and its data"
}

# 7. the plan
run_all default "" yes
NF_DEFAULT=$NFLINE
[ -n "$NF_DEFAULT" ] || fail "the default run made no nextflow call"
output_lacks default '53 Hard loci'
if grep -qE -- '--tools [^ ]*locityper|--locityper_' <<<"$NF_DEFAULT"; then
  fail "a default run passes locityper to Nextflow: ${NF_DEFAULT}"
fi
missed=$(opted_in default)
[ "$(grep -c . <<<"$missed")" -eq 2 ] || fail "control: the opted-in checks do not both fail on the default run: ${missed}"

run_all on "pharmcat,locityper" yes
missed=$(opted_in on)
[ -z "$missed" ] || fail "TOOLS=pharmcat,locityper: ${missed}; nextflow: ${NFLINE}"

run_all nodata "pharmcat,locityper" no
output_has nodata "53 Hard loci \(Locityper\)     skipped    \(data not installed: reference/locityper-"
if grep -qF -- '--locityper_data' <<<"$NFLINE"; then fail "no database, yet --locityper_data: ${NFLINE}"; fi
echo "PASS run-all.sh plans step 53 only when TOOLS names it and its database is installed"

# 8. the switches run_status.tsv keeps
meta() { awk -F'\t' -v k="$2" '$1 == "meta" && $2 == k { print $3 }' "${CASE_WORK}/$1/g/s1/logs/run_status.tsv" 2>/dev/null || true; }
ALT="${CASE_WORK}/elsewhere/locityper"
mkdir -p "$ALT"
echo ok > "${ALT}/DONE"
run_all alt "pharmcat,locityper" no LOCITYPER_DATA="$ALT"
grep -qF -- "--locityper_data ${ALT}" <<<"$NFLINE" || fail "LOCITYPER_DATA does not reach Nextflow: ${NFLINE}"
cmd=$(meta alt command)
case "$cmd" in *"LOCITYPER_DATA=${ALT} "*) ;; *) fail "the recorded command lacks LOCITYPER_DATA: ${cmd}" ;; esac
case "$cmd" in *"LOCITYPER_BG_REGION=${BG} "*) ;; *) fail "the recorded command lacks LOCITYPER_BG_REGION: ${cmd}" ;; esac
switches=" $(meta alt switches) "
for v in LOCITYPER_DATA LOCITYPER_BG_REGION; do
  case "$switches" in *" ${v} "*) ;; *) fail "'meta switches' lacks ${v}: ${switches}" ;; esac
done
cmd=$(meta default command)
[ -n "$cmd" ] || fail "control: the default run recorded no command"
case "$cmd" in *LOCITYPER_DATA=*) fail "control: a run without LOCITYPER_DATA records it: ${cmd}" ;; esac
echo "PASS LOCITYPER_DATA and LOCITYPER_BG_REGION are recorded switches"

# 9. the DONE guards
G="${CASE_WORK}/guard/g"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
DATA="${G}/reference/locityper-${PANEL}"
seed_sample "$G" s1
mkdir -p "$DATA"
printf '>chr1\nACGT\n' > "$REF"
printf 'chr1\t4\t6\t4\t5\n' > "${REF}.fai"
echo x > "${DATA}/counts.jf"
printf 'chr1\t0\t4\tX\n' > "${DATA}/loci.bed"
command -v sha256sum >/dev/null 2>&1 || fail "no sha256sum on this machine"
sha() { sha256sum "$1" | cut -d' ' -f1; }
GOOD="fai=$(sha "${REF}.fai") fasta=$(sha "$REF") image=x panel=y loci=z"
OTHER="fai=$(sha "${REF}.fai") fasta=$(printf '>chr1 TGCA' | sha256sum | cut -d' ' -f1) image=x panel=y loci=z"
put_done() { if [ -z "$1" ]; then rm -f "${DATA}/DONE"; else printf '%s\n' "$1" > "${DATA}/DONE"; fi; }

step53() {
  put_done "$2"
  run_rc "$1" env HOME="${CASE_WORK}/guard" GENOME_DIR="$G" bash "${SCRIPTS}/53-locityper.sh" s1
}
step53 s53-other "$OTHER"
output_has s53-other 'built from another reference'
output_lacks s53-other 'Database:'
step53 s53-good "$GOOD"
output_lacks s53-good 'another reference'
output_has s53-good 'Database: '
step53 s53-none ""
output_has s53-none 'database not found'
output_lacks s53-none 'Database:'

# The module's guard as the task runs it: its script lines before Locityper,
# with Nextflow's ${...} filled in and \$ unescaped.
GUARD="${CASE_WORK}/guard.sh"
awk '/^    script:/ { s = 1; next } s && /"""/ { g = 1; next } g && /^[[:space:]]*mkdir -p work/ { exit } g' \
  "${REPO_ROOT}/modules/local/locityper/main.nf" > "${GUARD}.in"
grep -q 'DONE' "${GUARD}.in" || fail "module: the guard is not in the script block"
# shellcheck disable=SC2016  # the single-quoted ${...} are the literal text to replace
{
  while IFS= read -r line; do
    line=${line//'${data}'/$DATA}
    line=${line//'${fai}'/${REF}.fai}
    line=${line//'${fasta}'/$REF}
    line=${line//'\$'/'$'}
    printf '%s\n' "$line"
  done < "${GUARD}.in"
  echo 'echo GUARD-PASSED'
} > "$GUARD"
module() { put_done "$2"; run_rc "$1" bash "$GUARD"; }
module mod-other "$OTHER"
output_has mod-other 'built from another reference'
module mod-good "$GOOD"
output_has mod-good 'GUARD-PASSED'
module mod-none ""
output_has mod-none 'not finished \(no DONE\)'
module mod-old "fai=$(sha "${REF}.fai") image=x"
output_lacks mod-old 'GUARD-PASSED'
echo "PASS step 53 and the module stop on a missing DONE or another reference, and pass on this one"
