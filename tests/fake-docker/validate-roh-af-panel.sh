#!/usr/bin/env bash
# validate-setup.sh reports step 11's allele-frequency panel (roh_af/, the
# file of ROH_AF_PANEL_URL): missing is a [WARN] (step 11 runs with AF 0.4),
# a file that is not the pinned one is a [FAIL], exit 1; the pinned file
# with no index is a [WARN] (step 11 does not use it), and with its .tbi a
# [PASS]. The last two run a checkout copy whose versions.env pins the
# planted file's sha256, so no real panel is needed.
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

# A checkout copy that pins the planted panel: every entry is the real
# checkout's but versions.env, so validate-setup.sh reads that pin.
COPY="${CASE_WORK}/checkout"
mkdir -p "${COPY}/scripts"
for e in "${REPO_ROOT}"/* "${REPO_ROOT}"/.[!.]*; do
  n=${e##*/}
  case "$n" in scripts | versions.env) ;; *) ln -s "$e" "${COPY}/${n}" ;; esac
done
for e in "${REPO_ROOT}"/scripts/*; do ln -s "$e" "${COPY}/scripts/${e##*/}"; done
rm "${COPY}/scripts/validate-setup.sh"
cp "${SCRIPTS}/validate-setup.sh" "${COPY}/scripts/"
SUM=$(sha256sum "$PANEL" | cut -d' ' -f1)
sed "s/^ROH_AF_PANEL_SHA256=.*/ROH_AF_PANEL_SHA256=\"${SUM}\"/" "${REPO_ROOT}/versions.env" > "${COPY}/versions.env"
grep -qx "ROH_AF_PANEL_SHA256=\"${SUM}\"" "${COPY}/versions.env" || fail "the copy does not pin the planted panel"

rm -f "${PANEL}.tbi"
TOOLS=pharmcat run_expect 0 no-index "${COPY}/scripts/validate-setup.sh" sample1
output_has no-index '\[WARN\].*ROH allele-frequency panel \(step 11\) has no index, so step 11 gives every site AF 0\.4: .*\.tbi'
output_lacks no-index '\[(FAIL|PASS)\].*ROH allele-frequency panel'

printf 'index\n' > "${PANEL}.tbi"
TOOLS=pharmcat run_expect 0 indexed "${COPY}/scripts/validate-setup.sh" sample1
output_has indexed '\[PASS\].*ROH allele-frequency panel \(step 11\): present'
output_lacks indexed '\[(WARN|FAIL)\].*ROH allele-frequency panel'
