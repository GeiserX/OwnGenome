#!/usr/bin/env bash
# validate-setup.sh reports step 11's allele-frequency panel (roh_af/, the
# file of ROH_AF_PANEL_URL): missing is a [WARN] (step 11 runs with AF 0.4),
# and a file that is not the pinned one is a [FAIL], exit 1.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
seed_clinvar "$GENOME_DIR"
PANEL_URL=$(sed -n 's/^ROH_AF_PANEL_URL="\(.*\)"$/\1/p' "${REPO_ROOT}/versions.env")
[ -n "$PANEL_URL" ] || fail "versions.env has no ROH_AF_PANEL_URL"
PANEL="${GENOME_DIR}/roh_af/${PANEL_URL##*/}"

TOOLS=pharmcat run_expect 0 missing "${SCRIPTS}/validate-setup.sh" sample1
output_has missing '\[WARN\].*ROH allele-frequency panel \(step 11 gives every site AF 0\.4 without it\) not found at: .*/roh_af/'
output_lacks missing '\[FAIL\].*ROH allele-frequency panel'

mkdir -p "${GENOME_DIR}/roh_af"
printf 'chr1\t100\tA,G\t0.1\n' > "$PANEL"
TOOLS=pharmcat run_expect 1 not-pinned "${SCRIPTS}/validate-setup.sh" sample1
output_has not-pinned '\[FAIL\].*ROH allele-frequency panel \(step 11\): .* is not the pinned file'
