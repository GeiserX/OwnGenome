#!/usr/bin/env bash
# validate-setup.sh fails on a step-40 site table that is not the pinned
# build, so run-all.sh, which runs it first, stops before such a table is used:
#   1. a table built before og took each site's strand from the liftover chain
#      (no "## table_format: 2" line) is a FAIL that says step 40 refuses it
#      and names the rebuild, or deleting the table when step 40 is not used;
#   2. a table with the format line that differs from the pin is a FAIL that
#      names the same two ways out;
#   3. the pinned table passes (a checkout copy pins the planted table).
# Before the fix, case 1 was a warning and validate-setup.sh exited 0.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
G="${CASE_WORK}/genome" P="${REPO_ROOT}/tests/parity/inputs/consumer"
export GENOME_DIR="$G"
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
SITES="${G}/reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
mkdir -p "$(dirname "$SITES")"

# 1. The old table (the unit table as og wrote it before the format line).
cp "${P}/old-sites.tsv.gz" "$SITES"
run_expect 1 old "${SCRIPTS}/validate-setup.sh" sample1
output_has old "\[FAIL\].*Site table ${SITES} was built before og took each site's strand from the liftover chain"
output_has old 'step 40 refuses it'
output_has old "Build it again: \./scripts/setup\.sh --consumer-sites ${G}, or delete it if you do not use step 40: rm ${SITES}$"
output_lacks old '\[OK\].*Site table'

# 2. A table with the format line that is not the pinned build.
cp "${P}/unit-sites.tsv.gz" "$SITES"
run_expect 1 unpinned "${SCRIPTS}/validate-setup.sh" sample1
output_has unpinned "\[FAIL\].*Site table ${SITES} differs from the pinned build"
output_has unpinned "Build it again: \./scripts/setup\.sh --consumer-sites ${G}, or delete it if you do not use step 40: rm ${SITES}$"
output_lacks unpinned 'built before og took'

# 3. A checkout copy whose pin is this table: it passes.
COPY="${CASE_WORK}/checkout"
mkdir -p "${COPY}/scripts"
for e in "${REPO_ROOT}"/*; do
  case "${e##*/}" in scripts|versions.env) ;; *) ln -s "$e" "${COPY}/${e##*/}" ;; esac
done
for e in "${REPO_ROOT}"/scripts/*; do ln -s "$e" "${COPY}/scripts/${e##*/}"; done
if command -v sha256sum > /dev/null 2>&1; then pin=$(gzip -dc "$SITES" | sha256sum | cut -d' ' -f1)
else pin=$(gzip -dc "$SITES" | shasum -a 256 | cut -d' ' -f1); fi
sed "s/^CONSUMER_SITES_TABLE_SHA256=.*/CONSUMER_SITES_TABLE_SHA256=\"${pin}\"/" "${REPO_ROOT}/versions.env" > "${COPY}/versions.env"
grep -qx "CONSUMER_SITES_TABLE_SHA256=\"${pin}\"" "${COPY}/versions.env" || fail "the checkout copy does not pin the planted table"
run_expect 0 pinned "${COPY}/scripts/validate-setup.sh" sample1
output_has pinned "\[OK\].*Site table of the raw-data export .*the same as the pinned build"
output_lacks pinned '\[FAIL\].*Site table'
echo "PASS validate-setup.sh fails on a site table that is not the pinned build and passes the pinned one"
