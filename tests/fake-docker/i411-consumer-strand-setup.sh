#!/usr/bin/env bash
# setup.sh --consumer-sites keeps a step-40 site table only when it is the
# pinned build. Any other table, such as one built before og took each site's
# strand from the liftover chain, is built again in its place; the table's
# file name never changes, so an upgraded install keeps its old one otherwise.
# A checkout copy pins small inputs the fake download serves (the parity
# build-sites inputs), and og consumer build-sites runs for real on the host.
#   1. a build that does not give the pinned table exits 1;
#   2. with the pin set to that build, an old table there (no format line,
#      one strand flipped) is built again, and setup exits 0 with the pinned
#      table in its place;
#   3. the pinned table is kept as it is, with no build.
# Before the fix, case 2 warned, kept the old table and exited 0 (and case 1
# exited 0).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

A=$(uname -m)
[ -x "${REPO_ROOT}/bin/.og/og-linux-${A}" ] || fail "bin/.og/og-linux-${A} is not built: this case runs og"
sha() { if command -v sha256sum > /dev/null 2>&1; then sha256sum | cut -d' ' -f1; else shasum -a 256 | cut -d' ' -f1; fi; }

COPY="${CASE_WORK}/checkout" S="${CASE_WORK}/served" G="${CASE_WORK}/genome" P="${REPO_ROOT}/tests/parity/inputs/consumer"
mkdir -p "${COPY}/scripts" "${COPY}/rust" "$S"
for d in scripts/lib scripts/ci bin assets; do ln -s "${REPO_ROOT}/${d}" "${COPY}/${d}"; done
cp "${REPO_ROOT}/scripts/setup.sh" "${COPY}/scripts/setup.sh"
printf '#!/usr/bin/env bash\nexit 0\n' > "${COPY}/rust/og.sh"
chmod +x "${COPY}/rust/og.sh"
# The served inputs, by the base names of their pinned URLs.
cp "${P}/build-tpl.txt.gz" "${S}/23andme_v5_hg19_ref.txt.gz"
cp "${P}/build-dbsnp.vcf.gz" "${S}/GCF_000001405.40.gz"
gzip -nc "${P}/build-chain.txt" > "${S}/hg19ToHg38.over.chain.gz"
# pin_table SHA: the checkout copy's versions.env, pinning the served inputs and SHA as the table.
pin_table() {
  { grep -v '^CONSUMER_SITES_SHA256=\|^CONSUMER_DBSNP_MD5=\|^HG19_TO_HG38_CHAIN_SHA256=\|^CONSUMER_SITES_TABLE_SHA256=' "${REPO_ROOT}/versions.env"
    echo "CONSUMER_SITES_SHA256=\"$(sha < "${S}/23andme_v5_hg19_ref.txt.gz")\""
    echo "CONSUMER_DBSNP_MD5=\"$(fake_md5 < "${S}/GCF_000001405.40.gz")\""
    echo "HG19_TO_HG38_CHAIN_SHA256=\"$(sha < "${S}/hg19ToHg38.over.chain.gz")\""
    echo "CONSUMER_SITES_TABLE_SHA256=\"$1\""
  } > "${COPY}/versions.env"
}
seed_reference "$G"
use_output_hook
cat > "${CASE_WORK}/hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift  # the image
map() {
  case "$1" in
    /pgp-bin/*) printf '%s' "${REPO_ROOT}/bin/${1#/pgp-bin/}" ;;
    /genome*) host_path "$1" ;;
    *) printf '%s' "$1" ;;
  esac
}
case "$*" in
  "/pgp-bin/og "*) args=(); for a in "$@"; do args+=("$(map "$a")"); done; exec "${args[@]}" ;;
esac
HOOK
chmod +x "${CASE_WORK}/hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook" FAKE_DOWNLOAD_DIR="$S" SETUP_DISK_CHECK=false FETCH_TRIES=1 FETCH_WAIT=0 \
  OG="${REPO_ROOT}/bin/og"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
TABLE="${G}/reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
builds() { grep -c 'consumer build-sites' "$FAKE_DOCKER_LOG" || true; }

# 1. A pin no build gives: setup builds the table and fails.
pin_table 0000000000000000000000000000000000000000000000000000000000000000
run_rc first "${COPY}/scripts/setup.sh" --consumer-sites "$G"
FIRST_RC=$RC
[ -s "$TABLE" ] || fail "setup.sh --consumer-sites built no table"
[ "$(builds)" -eq 1 ] || fail "og consumer build-sites ran $(builds) times, not once"
PIN=$(gzip -dc "$TABLE" | sha)

# 2. The pin is that build; the old table there is built again in its place.
pin_table "$PIN"
gzip -dc "$TABLE" | grep -v '^## table_format:' \
  | awk -F'\t' -v OFS='\t' '!/^#/ && !done && $7 == "+" { $7 = "-"; done = 1 } 1' | gzip -c > "${CASE_WORK}/old.tsv.gz"
[ "$(gzip -dc "${CASE_WORK}/old.tsv.gz" | sha)" != "$PIN" ] || fail "the planted old table is the pinned one"
cp "${CASE_WORK}/old.tsv.gz" "$TABLE"
run_expect 0 stale "${COPY}/scripts/setup.sh" --consumer-sites "$G"
output_has stale "\[WARN\] Site table of step 40 already present: ${TABLE} differs from the pinned build"
output_has stale 'building it again in its place'
output_has stale "\[OK\] Site table of step 40 built .*${TABLE}, the same as the pinned build"
[ "$(gzip -dc "$TABLE" | sha)" = "$PIN" ] || fail "the old table was kept: setup.sh did not build it again"
[ "$(builds)" -eq 2 ] || fail "og consumer build-sites ran $(builds) times, not twice"

# 3. The pinned table is kept, with no build.
run_expect 0 pinned "${COPY}/scripts/setup.sh" --consumer-sites "$G"
output_has pinned "\[OK\] Site table of step 40 already present: ${TABLE}, the same as the pinned build"
[ "$(builds)" -eq 2 ] || fail "setup.sh built the pinned table again"

gzip -dc "$TABLE" | grep -q '^## table_format: 2 ' || fail "the built table has no '## table_format: 2' line"
[ "$FIRST_RC" -eq 1 ] || fail "a build that does not give the pinned table exited ${FIRST_RC}, not 1"
output_has first "differs from the pinned build"
output_has first "does not give the pinned table"
echo "PASS setup.sh --consumer-sites builds a table that is not the pinned build again and keeps the pinned one"
