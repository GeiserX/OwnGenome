#!/usr/bin/env bash
# Step 06's always-report list reads FILTER and the depth floor, and step 06
# counts the sample records whose REF is not the reference's, in both mirrors:
#   - scripts/06-clinvar-screen.sh: bcftools norm -c w keeps a record whose REF
#     differs from the reference and only warns REF_MISMATCH on stderr. The
#     summary counts those warnings and prints a notice when there are any,
#     so a mismatched record no longer leaves only a log line;
#   - it hands og clinvar deep-intronic the depth floor (CLINVAR_MIN_DP), and
#     its query format ends with %FILTER, so a filtered or shallow reference
#     call no longer reads as plain "ref";
#   - modules/local/clinvar_screen/main.nf builds the same query format (its
#     query_format runs here, Groovy escapes undone), passes ext.min_dp, and
#     counts REF_MISMATCH from its own bcftools norm.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook

# The answers step 06 reads from bcftools, every call's outputs, the two
# REF_MISMATCH warnings bcftools norm writes for the sample when
# ${CASE_WORK}/mismatch exists, and the two cDNA tables og writes relative to
# its -w directory, so the step runs to its summary.
cat > "${CASE_WORK}/hook-step06" <<'HOOK'
#!/usr/bin/env bash
"${CASE_WORK:?}/hook-outputs" "$@" || exit $?
shift
out="${CASE_WORK}/genome/sample1/clinvar"
case "$*" in
  *"index -s "*) printf 'chr1\t248956422\t3\n' ;;
  *"index -n "*) echo 3 ;;
  *"view -H -f PASS "*) echo 1 ;;
  *"view -h "*) printf '##fileformat=VCFv4.2\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' ;;
  *"--multi-overlaps"*)
    if [ -e "${CASE_WORK}/mismatch" ]; then
      printf 'REF_MISMATCH\tchr1\t40\tA\tT\nREF_MISMATCH\tchr1\t150\tGT\tCT\n' >&2
    fi
    printf 'Lines   total/split/joined/realigned/mismatch_removed/dup_removed/skipped:\t4/0/0/0/0/0/0\n' >&2 ;;
  *"cdna check "*) : > "${out}/sample1_clinvar_cdna_check.tsv.part" ;;
  *"cdna scan "*) : > "${out}/sample1_clinvar_cdna_scan.tsv.part" ;;
esac
exit 0
HOOK
chmod +x "${CASE_WORK}/hook-step06"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-step06"

COUNT_LINE='^Records whose REF differs from the reference \(REF_MISMATCH\): '

# --- scripts/06: two mismatched records ---------------------------------------------
touch "${CASE_WORK}/mismatch"
: > "$FAKE_DOCKER_LOG"
run_rc step06-mismatch "${SCRIPTS}/06-clinvar-screen.sh" sample1
expect_rc step06-mismatch 0
output_has step06-mismatch "^REF_MISMATCH	chr1	40	A	T$"
output_has step06-mismatch "${COUNT_LINE}2$"
output_has step06-mismatch '^NOTICE: 2 sample records have a REF that is not the reference'"'"'s'
DEEP_CALL=$(awk '/^run / && /og clinvar deep-intronic /' "$FAKE_DOCKER_LOG")
echo "deep-intronic call: ${DEEP_CALL:-none}"
grep -qE -- '--min-gq 30 --min-dp 10 ' <<<"$DEEP_CALL" \
  || fail "step 06 does not hand og clinvar deep-intronic its depth floor (--min-dp 10): ${DEEP_CALL:-no call}"
QUERY=$(awk '/^run / && /bcftools query -T /' "$FAKE_DOCKER_LOG")
grep -qF -- 'FILTER' <<<"$QUERY" || fail "step 06's always-report query reads no FILTER: ${QUERY:-no query}"

# --- scripts/06: none, and another depth floor -----------------------------------
rm -f "${CASE_WORK}/mismatch"
rm -rf "${GENOME_DIR}/sample1/clinvar"
: > "$FAKE_DOCKER_LOG"
CLINVAR_MIN_DP=8 run_rc step06-clean "${SCRIPTS}/06-clinvar-screen.sh" sample1
expect_rc step06-clean 0
output_has step06-clean "${COUNT_LINE}0$"
output_lacks step06-clean '^NOTICE: .* REF that is not the reference'
grep -qE -- '--min-dp 8 ' <<<"$(awk '/^run / && /og clinvar deep-intronic /' "$FAKE_DOCKER_LOG")" \
  || fail "step 06 does not pass CLINVAR_MIN_DP=8 to og clinvar deep-intronic"

# --- both mirrors build the same query format, ending with FILTER -----------------
WANT='%CHROM\t%POS\t%END\t%REF\t%ALT\t[%GT]\t[%GQ]\t[%DP]\t%FILTER\n'
HEADER="${CASE_WORK}/gq-dp-header.txt"
printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n' > "$HEADER"
# fn_text FILE: the query_format function of FILE, as text.
fn_text() { awk '/^ *query_format\(\) \{/ { on = 1 } on { print } on && /^ *\}$/ { exit }' "$1"; }
BASH_FN=$(fn_text "${SCRIPTS}/06-clinvar-screen.sh")
[ -n "$BASH_FN" ] || fail "scripts/06 has no query_format function"
# shellcheck disable=SC2016  # the function's text runs in the inner shell
GOT=$(HEADER=$HEADER bash -c 'set -euo pipefail; bcftools_run() { cat "$HEADER"; }; cpath() { printf %s "$1"; }; eval "$1"; query_format x.vcf.gz' _ "$BASH_FN")
[ "$GOT" = "$WANT" ] || fail "scripts/06 queries '${GOT}', want '${WANT}'"
MOD="${REPO_ROOT}/modules/local/clinvar_screen/main.nf"
MOD_FN=$(fn_text "$MOD" | sed -e 's/\\\$/$/g' -e 's/\\\\/\\/g')
[ -n "$MOD_FN" ] || fail "the ClinVar screen module has no query_format function"
# shellcheck disable=SC2016  # the function's text runs in the inner shell
GOT=$(cd "$CASE_WORK" && HEADER=$HEADER bash -c 'set -euo pipefail; bcftools() { cat "$HEADER"; }; eval "$1"; query_format x.vcf.gz' _ "$MOD_FN")
[ "$GOT" = "$WANT" ] || fail "the module queries '${GOT}', want '${WANT}'"

# --- the module: depth floor and REF_MISMATCH count ---------------------------------
# shellcheck disable=SC2016  # the module's own text, dollar signs included
grep -qE -- 'og clinvar deep-intronic .*--min-gq \$\{min_gq\} --min-dp \$\{min_dp\} ' "$MOD" \
  || fail "the module does not hand og clinvar deep-intronic ext.min_dp"
# shellcheck disable=SC2016  # the module's own text
grep -qF -- '--multi-overlaps . -c w -f ${reference} -Oz -o ${meta.id}_pass.vcf.gz - 2> norm.log' "$MOD" \
  || fail "the module's sample bcftools norm does not keep its REF_MISMATCH warnings"
grep -qF -- "REF_MISMATCH=\\\$(grep -c '^REF_MISMATCH' norm.log || true)" "$MOD" \
  || fail "the module does not count the REF_MISMATCH warnings"
# shellcheck disable=SC2016  # the module's own text
grep -qF -- 'echo "Records whose REF differs from the reference (REF_MISMATCH): \${REF_MISMATCH}"' "$MOD" \
  || fail "the module's log has no REF_MISMATCH count"
echo "both mirrors read FILTER and the depth floor, and count REF_MISMATCH: PASS"
