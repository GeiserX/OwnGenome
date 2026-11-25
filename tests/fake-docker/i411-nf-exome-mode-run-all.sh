#!/usr/bin/env bash
# run-all.sh with --model_type WES (an exome). main.nf leaves the
# exome_skipped_tools() of workflows/exome.nf out of an exome, and indexcov
# (step 16). Nextflow then exits 0 without running them, so run-all.sh must
# not record them as done: they are not passed in --tools and their rows in
# logs/run_status.tsv read 'skipped (not run: exome)', which og-runlog shows
# as skipped with that reason.
#   1. WES with every tool of the list named in TOOLS: each step reads
#      'skipped (not run: exome)', none is in --tools, step 16 the same, no
#      other step says it; survivor_merge keeps Manta and Delly and runs.
#   2. WES with Manta and CNVpytor as the SV callers: the merge (22) goes
#      with CNVpytor, as exome_dropped() decides.
#   3. WES with the default TOOLS: ROH, TelomereHunter and CNVpytor are not
#      run, no step outside the list says so, the merge keeps its two callers.
#   4. The control, a whole genome: the same steps run and are ok, and no
#      step says 'not run: exome'.
# The list is read here from workflows/exome.nf on its own, so run-all.sh
# reading another list, or none, turns 1 red.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
ST="${G}/sample1/logs/run_status.tsv"
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook
BED="${CASE_WORK}/capture.bed"
printf 'chr1\t1000\t2000\n' > "$BED"
NOT_RUN='skipped (not run: exome)'

# The --tools names an exome leaves out, as workflows/exome.nf lists them, and their steps.
mapfile -t LIST < <(awk '/^def exome_skipped_tools\(\)/ {f = 1; next} f && /^}/ {exit} f' "${REPO_ROOT}/workflows/exome.nf" \
  | grep -oE "'[a-z0-9_]+'" | tr -d "'")
declare -A STEP=([cnvpytor]=18 [roh]=11 [telomere_hunter]=10 [archaic]=43 [mocha]=51 [viral]=47 [cyrius]=21
  [parascopy]=35 [locityper]=53 [mei]=50)
[ "${#LIST[@]}" -ge "${#STEP[@]}" ] || fail "read ${#LIST[@]} names from exome_skipped_tools() of workflows/exome.nf: ${LIST[*]}"
for t in "${LIST[@]}"; do [ -n "${STEP[$t]:-}" ] || fail "workflows/exome.nf leaves out '${t}': add its step number to STEP here"; done
for t in "${!STEP[@]}"; do [[ " ${LIST[*]} " == *" ${t} "* ]] || fail "workflows/exome.nf no longer lists '${t}': update STEP and this case"; done

last_nf() { grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -n 1; }
tools_of() { sed -n 's/.* --tools \([^ ]*\) .*/\1/p' <<< "$1" | tr -d "\\\\"; }  # printf %q may escape the commas
row() { awk -F'\t' -v s="$1" '$1 == "step" && $2 == s {r = $3} END {print r}' "$ST"; }
not_run_steps() { awk -F'\t' -v w="$NOT_RUN" '$1 == "step" && $3 == w {print $2}' "$ST" | sort -u | tr '\n' ' '; }

# --- 1. an exome with every tool of the list ------------------------------------------------
ALL=$(IFS=,; echo "${LIST[*]}"),manta,delly,survivor_merge,mosdepth
run_expect 0 wes env TOOLS="$ALL" "${SCRIPTS}/run-all.sh" sample1 male --model_type WES --capture_bed "$BED"
NFLINE=$(last_nf)
grep -qF -- "--model_type WES --capture_bed $(printf '%q' "$BED")" <<< "$NFLINE" || fail "nextflow was not given the exome options: ${NFLINE}"
TL=$(tools_of "$NFLINE")
for s in 16 "${STEP[@]}"; do [ "$(row "$s")" != ok ] || fail "WES: step ${s} reads ok, but main.nf does not run it on an exome"; done
want="16"
for t in "${LIST[@]}"; do
  s=${STEP[$t]}
  got=$(row "$s")
  [ "$got" = "$NOT_RUN" ] || fail "WES: step ${s} (${t}) reads '${got}', want '${NOT_RUN}'"
  [[ ",${TL}," != *",${t},"* ]] || fail "WES: ${t} is in --tools (${TL})"
  want+=" ${s}"
done
got=$(row 16)
[ "$got" = "$NOT_RUN" ] || fail "WES: step 16 (indexcov) reads '${got}', want '${NOT_RUN}'"
want=$(tr ' ' '\n' <<< "$want" | sort -u | tr '\n' ' ')
[ "$(not_run_steps)" = "$want" ] || fail "WES: the steps that read '${NOT_RUN}' are '$(not_run_steps)', want '${want}'"
output_has wes '^  11 ROH +skipped +\(not run: exome\)$'
output_has wes '^  22 SV consensus merge +runs$'
[ "$(row 22)" = ok ] || fail "WES: survivor_merge with Manta and Delly reads '$(row 22)', want ok"
for t in manta delly survivor_merge mosdepth; do [[ ",${TL}," == *",${t},"* ]] || fail "WES: ${t} is not in --tools (${TL})"; done

# --- 2. an exome where CNVpytor was one of the merge's two callers ----------------------------
run_expect 0 wes-two env TOOLS=manta,cnvpytor,survivor_merge,mosdepth "${SCRIPTS}/run-all.sh" sample1 male --model_type=WES --capture_bed "$BED"
TL=$(tools_of "$(last_nf)")
for s in 18 22; do [ "$(row "$s")" = "$NOT_RUN" ] || fail "WES, Manta and CNVpytor: step ${s} reads '$(row "$s")', want '${NOT_RUN}'"; done
[ "$(row 04)" = ok ] || fail "WES, Manta and CNVpytor: Manta reads '$(row 04)', want ok"
[ "$TL" = mosdepth,manta ] || fail "WES, Manta and CNVpytor: --tools is '${TL}', want mosdepth,manta"

# --- 3. an exome with the default TOOLS ------------------------------------------------------
run_expect 0 wes-default "${SCRIPTS}/run-all.sh" sample1 male --model_type WES --capture_bed "$BED"
TL=$(tools_of "$(last_nf)")
for s in 10 11 16 18; do [ "$(row "$s")" = "$NOT_RUN" ] || fail "WES, default: step ${s} reads '$(row "$s")', want '${NOT_RUN}'"; done
for s in $(not_run_steps); do
  [[ " ${want} " == *" ${s} "* ]] || fail "WES, default: step ${s} reads '${NOT_RUN}', but workflows/exome.nf does not leave it out"
done
for t in roh telomere_hunter; do [[ ",${TL}," != *",${t},"* ]] || fail "WES, default: ${t} is in --tools (${TL})"; done
[ "$(row 22)" = ok ] || fail "WES, default: survivor_merge reads '$(row 22)', want ok"

# --- 4. the control: a whole genome ----------------------------------------------------------
run_expect 0 wgs env TOOLS=roh,telomere_hunter,manta,delly,survivor_merge,mosdepth "${SCRIPTS}/run-all.sh" sample1 male
TL=$(tools_of "$(last_nf)")
for s in 10 11 16 22; do [ "$(row "$s")" = ok ] || fail "WGS: step ${s} reads '$(row "$s")', want ok"; done
[ "$TL" = roh,mosdepth,telomere_hunter,manta,delly,survivor_merge ] || fail "WGS: --tools is '${TL}'"
[ -z "$(not_run_steps)" ] || fail "WGS: steps read '${NOT_RUN}': $(not_run_steps)"
output_lacks wgs 'not run: exome'
echo "PASS: run-all.sh records an exome's dropped steps as skipped (not run: exome), never ok"
