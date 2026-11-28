#!/usr/bin/env bash
# Step 40 never exports from a site table built before og took each site's
# strand from the liftover chain. Such a table has no "## table_format: 2"
# line, and it has the wrong bases at sites on inverted blocks; its file name
# is the same as the new table's, so an upgraded install still has it.
#   1. scripts/40-export-consumer.sh refuses it, names the rebuild with this
#      install's folder, and starts no export;
#   2. og consumer export, which the Nextflow module calls directly, refuses
#      it too (exit 2) and names the rebuild;
#   3. the same table with the format line exports.
# Before the fix, case 1 exported the old table.
# og runs on the host for the host's arch (bin/og, which Guard's fake-docker
# job builds).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

A=$(uname -m)
[ -x "${REPO_ROOT}/bin/.og/og-linux-${A}" ] || fail "bin/.og/og-linux-${A} is not built: this case runs og"
# The step's og check asks for the host's arch only; the fake docker runs og here.
export OG_ARCHES=$A

# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
G="${CASE_WORK}/genome" P="${REPO_ROOT}/tests/parity/inputs/consumer"
export GENOME_DIR="$G"
SITES="${G}/reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
OUT="${G}/S/consumer_export"
mkdir -p "${G}/S/vcf" "$(dirname "$SITES")"
cp "${P}/unit.g.vcf.gz" "${G}/S/vcf/S.g.vcf.gz"
# old-sites.tsv.gz is the unit table as og wrote it before the format line;
# unit-sites.tsv.gz is the same table with it.
cp "${P}/old-sites.tsv.gz" "$SITES"
if gzip -dc "$SITES" | grep -q '^## table_format:'; then fail "the planted old table has a format line"; fi
[ "$(gzip -dc "${P}/unit-sites.tsv.gz" | grep -v '^## table_format: 2 ')" = "$(gzip -dc "$SITES")" ] \
  || fail "unit-sites.tsv.gz is not old-sites.tsv.gz with one format line"

use_output_hook
cat > "${CASE_WORK}/hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift  # the image
map() {
  case "$1" in
    /pgp-bin/*) printf '%s' "${REPO_ROOT}/bin/${1#/pgp-bin/}" ;;
    /pgp-assets/*) printf '%s' "${REPO_ROOT}/assets/${1#/pgp-assets/}" ;;
    /genome*) host_path "$1" ;;
    *) printf '%s' "$1" ;;
  esac
}
case "$*" in
  "/pgp-bin/og "*) args=(); for a in "$@"; do args+=("$(map "$a")"); done; exec "${args[@]}" ;;
esac
HOOK
chmod +x "${CASE_WORK}/hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook"

# 1. The step refuses the old table and starts no export.
run_expect 1 old "${SCRIPTS}/40-export-consumer.sh" S
output_has old "the site table ${SITES} was built before og took each site's strand from the liftover chain"
output_has old "run scripts/setup.sh --consumer-sites ${G},"
[ ! -e "${OUT}/S_23andme_v5.zip" ] || fail "step 40 wrote a raw-data file from the old table"
if grep -q 'og consumer export' "$FAKE_DOCKER_LOG"; then fail "step 40 started og consumer export on the old table"; fi

# 2. og alone refuses it, as the Nextflow module runs it.
run_expect 2 og-old "${REPO_ROOT}/bin/og" consumer export --gvcf "${G}/S/vcf/S.g.vcf.gz" --sites "$SITES" \
  --sample S --sex male --par "${REPO_ROOT}/assets/par_grch38.bed" --out-dir "${CASE_WORK}/og-old"
output_has og-old "has no '## table_format: 2' line"
output_has og-old 'rebuild it with \./scripts/setup\.sh --consumer-sites <genome_dir>'
[ ! -e "${CASE_WORK}/og-old/S_23andme_v5.zip" ] || fail "og consumer export wrote a raw-data file from the old table"

# 3. With the format line, the same table exports.
cp "${P}/unit-sites.tsv.gz" "$SITES"
run_expect 0 new "${SCRIPTS}/40-export-consumer.sh" S
[ -s "${OUT}/S_23andme_v5.zip" ] || fail "step 40 wrote no raw-data file from the table with the format line"
grep -q 'og consumer export' "$FAKE_DOCKER_LOG" || fail "step 40 did not run og consumer export"
echo "PASS step 40 and og consumer export refuse a site table built before the strand came from the liftover chain"
