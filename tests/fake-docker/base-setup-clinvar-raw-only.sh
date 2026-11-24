#!/usr/bin/env bash
# setup.sh when only the raw ClinVar file and its index exist (an earlier
# run stopped before the derived files): it must build clinvar_chr.vcf.gz
# and clinvar_pathogenic_chr.vcf.gz instead of reporting ClinVar as done,
# and, once og is there, the conflicting subset from clinvar_chr.vcf.gz.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
mkdir -p "${GENOME_DIR}/clinvar"
printf 'placeholder\n' > "${GENOME_DIR}/clinvar/clinvar.vcf.gz"
printf 'placeholder\n' > "${GENOME_DIR}/clinvar/clinvar.vcf.gz.tbi"
use_output_hook

# og builds the conflicting subset: a stand-in for each Linux build bin/.og
# lacks satisfies setup.sh's og check, and goes when the case ends.
OGDIR="${REPO_ROOT}/bin/.og"
TREE=$("${REPO_ROOT}/rust/og.sh" tree-hash) || fail "could not compute the rust/ tree hash"
STANDINS=()
trap 'for f in "${STANDINS[@]}"; do rm -f "$f"; done' EXIT
for a in x86_64 aarch64; do
  [ ! -e "${OGDIR}/og-linux-${a}" ] || continue
  mkdir -p "$OGDIR"
  printf 'stand-in og-tree-hash:%s;\n' "$TREE" > "${OGDIR}/og-linux-${a}"; chmod 755 "${OGDIR}/og-linux-${a}"
  STANDINS+=("${OGDIR}/og-linux-${a}")
done
"${REPO_ROOT}/rust/og.sh" check --required > /dev/null || fail "bin/.og does not pass setup.sh's og check"

run_rc setup "${SCRIPTS}/setup.sh" "$GENOME_DIR"
awk '/^run / && /clinvar_chr\.vcf\.gz/ && !/clinvar_pathogenic/ { found = 1 } END { exit !found }' \
  "$FAKE_DOCKER_LOG" || fail "setup.sh did not build clinvar_chr.vcf.gz from the raw ClinVar file"
docker_log_has '^run .*clinvar_pathogenic_chr\.vcf\.gz' \
  "setup.sh did not build clinvar_pathogenic_chr.vcf.gz"
docker_log_has '^run .*conflicting-subset.*/genome/clinvar/clinvar_chr\.vcf\.gz' \
  "setup.sh did not build the conflicting subset from clinvar_chr.vcf.gz"
for f in clinvar_conflicting_plp_chr.vcf.gz clinvar_conflicting_plp_chr.vcf.gz.tbi; do
  [ -f "${GENOME_DIR}/clinvar/${f}" ] || fail "setup.sh did not install ${f}"
done
[ ! -e "${GENOME_DIR}/clinvar/.build-conflicting" ] || fail "setup.sh left clinvar/.build-conflicting behind"
output_has setup 'ClinVar conflicting subset \(records some labs call pathogenic, step 06\) built'
output_lacks setup 'unbound variable'
expect_rc setup 3
