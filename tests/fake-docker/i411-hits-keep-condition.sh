#!/usr/bin/env bash
# Step 06 copies ClinVar's condition onto each hit, in both mirrors:
#   - scripts/06-clinvar-screen.sh: the bcftools annotate that writes the hits
#     VCF copies INFO/CLNDN and INFO/CLNDISDB with ID, GENEINFO, CLNSIG and
#     CLNREVSTAT when ClinVar's header defines them (every ClinVar release
#     does). Without them a later step knows the gene but not which of its
#     diseases ClinVar asserted for the allele;
#   - a hand-built ClinVar file without them still runs: bcftools annotate
#     stops on a tag the file does not define, so they are left out, with a
#     notice;
#   - modules/local/clinvar_screen/main.nf builds the same column list from
#     the same header and hands it to its own bcftools annotate.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook

# ClinVar's own header lines for the two tags, from a ClinVar release.
CLNDN_LINE="##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB\">"
CLNDISDB_LINE='##INFO=<ID=CLNDISDB,Number=.,Type=String,Description="Tag-value pairs of disease database name and identifier submitted for germline classifications, e.g. OMIM:NNNNNN">'
BASE_HEADER='##fileformat=VCFv4.1
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene(s)">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Significance">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">'
CHROM_LINE=$(printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')

# The answers step 06 reads from bcftools: one contig with records on both
# sides, a PASS record, record counts, and ClinVar's header from
# ${CASE_WORK}/clinvar-header.txt. Every call also gets its outputs created.
cat > "${CASE_WORK}/hook-step06" <<'HOOK'
#!/usr/bin/env bash
"${CASE_WORK:?}/hook-outputs" "$@" || exit $?
shift
case "$*" in
  *"index -s "*) printf 'chr1\t248956422\t3\n' ;;
  *"index -n "*) echo 3 ;;
  *"view -H -f PASS "*) echo 1 ;;
  *"view -h "*clinvar_pathogenic_chr.norm.vcf.gz*) cat "${CASE_WORK}/clinvar-header.txt" ;;
  *"view -h "*) printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' ;;
esac
exit 0
HOOK
chmod +x "${CASE_WORK}/hook-step06"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-step06"

# The step stops later, at an output the fake docker does not write (an -o
# relative to -w), so its exit code is not checked: the annotate call
# comes before it.
# annotate_cols: the -c list of the hits VCF's bcftools annotate in the docker log.
annotate_cols() {
  awk '/^run / && /bcftools annotate / && /pair-logic exact/' "$FAKE_DOCKER_LOG" \
    | grep -oE -- "-c[^ ]* [^ ]*INFO/CLNREVSTAT[^ ]*" | head -n 1
}

# --- scripts/06: a ClinVar release defines both tags ------------------------------
printf '%s\n%s\n%s\n%s\n' "$BASE_HEADER" "$CLNDN_LINE" "$CLNDISDB_LINE" "$CHROM_LINE" > "${CASE_WORK}/clinvar-header.txt"
: > "$FAKE_DOCKER_LOG"
run_rc step06-release "${SCRIPTS}/06-clinvar-screen.sh" sample1
COLS=$(annotate_cols)
echo "annotate columns: ${COLS:-none}"
[ -n "$COLS" ] || fail "step 06 ran no bcftools annotate for the hits VCF"
grep -qF 'ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB' <<<"$COLS" \
  || fail "step 06 does not copy ClinVar's condition (INFO/CLNDN, INFO/CLNDISDB) onto the hits: ${COLS}"
output_lacks step06-release 'does not define CLNDN'

# --- scripts/06: a hand-built ClinVar file without them --------------------------
printf '%s\n%s\n' "$BASE_HEADER" "$CHROM_LINE" > "${CASE_WORK}/clinvar-header.txt"
rm -rf "${GENOME_DIR}/sample1/clinvar"
: > "$FAKE_DOCKER_LOG"
run_rc step06-handbuilt "${SCRIPTS}/06-clinvar-screen.sh" sample1
COLS=$(annotate_cols)
echo "annotate columns: ${COLS:-none}"
[ -n "$COLS" ] || fail "step 06 ran no bcftools annotate for a ClinVar file without CLNDN"
if grep -q 'CLNDN\|CLNDISDB' <<<"$COLS"; then
  fail "step 06 asks bcftools annotate for CLNDN or CLNDISDB, which this ClinVar file does not define (bcftools stops): ${COLS}"
fi
grep -qF 'ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT' <<<"$COLS" || fail "step 06 lost ClinVar's ID, gene, significance or review status: ${COLS}"
output_has step06-handbuilt 'does not define CLNDN and CLNDISDB'

# --- the Nextflow module: the same columns, from the same header ------------------
# Its clinvar_cols line, with the Groovy escapes undone, runs with a bcftools
# that prints the header; its annotate must take that list.
MOD="${REPO_ROOT}/modules/local/clinvar_screen/main.nf"
LINE=$(awk '/^process CLINVAR_SCREEN \{/ { on = 1 } on && /^ *clinvar_cols=/ { print; exit }' "$MOD")
[ -n "$LINE" ] || fail "the ClinVar screen module builds no clinvar_cols list for bcftools annotate"
LINE=$(sed -e 's/\\\$/$/g' -e 's/\\\\/\\/g' <<<"$LINE")
# shellcheck disable=SC2016  # the module's text, Groovy escape included
awk '/^process CLINVAR_SCREEN \{/ { on = 1 } on && /bcftools annotate .*pair-logic exact/ { getline; print; exit }' "$MOD" \
  | grep -qF -- '-c "\${clinvar_cols}"' || fail "the module's bcftools annotate does not take clinvar_cols"
# module_cols HEADER_FILE: the list the module's line builds from that header.
module_cols() {
  # shellcheck disable=SC2034  # read by the module's line
  CLINVAR_NORM=clinvar.norm.vcf.gz HEADER=$1 bash -c 'set -euo pipefail; bcftools() { cat "$HEADER"; }; eval "$1"; printf "%s" "$clinvar_cols"' _ "$LINE"
}
printf '%s\n%s\n%s\n%s\n' "$BASE_HEADER" "$CLNDISDB_LINE" "$CLNDN_LINE" "$CHROM_LINE" > "${CASE_WORK}/release-header.txt"
printf '%s\n%s\n' "$BASE_HEADER" "$CHROM_LINE" > "${CASE_WORK}/handbuilt-header.txt"
GOT=$(module_cols "${CASE_WORK}/release-header.txt")
[ "$GOT" = 'ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB' ] \
  || fail "the module copies '${GOT}' from a ClinVar release, want ClinVar's condition too (INFO/CLNDN, INFO/CLNDISDB)"
GOT=$(module_cols "${CASE_WORK}/handbuilt-header.txt")
[ "$GOT" = 'ID,INFO/GENEINFO,INFO/CLNSIG,INFO/CLNREVSTAT' ] \
  || fail "the module copies '${GOT}' from a ClinVar file without CLNDN and CLNDISDB, want only ID, GENEINFO, CLNSIG and CLNREVSTAT"
echo "both mirrors copy ClinVar's condition when ClinVar defines it: PASS"
