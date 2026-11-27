#!/usr/bin/env bash
# A downloads folder outside GENOME_DIR, named relative to the working
# directory (OG_DOWNLOADS=shared-dl), is made absolute: the download lands in
# it, a container sees it read-only at /downloads and is handed the file at
# its /downloads path, and the file stays there once installed.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

S="${CASE_WORK}/served" G="${CASE_WORK}/genome"
mkdir -p "$S" "$G"
printf 'fake zip\n' > "${S}/hla.dat.zip"
printf 'MD5 (hla.dat.zip) = %s\n' "$(md5sum < "${S}/hla.dat.zip" | cut -d' ' -f1)" > "${S}/md5checksum.txt"
# shellcheck disable=SC2016  # the sed pattern names a shell variable as text
REL=$(sed -n 's/^HLA_DB_RELEASE=\${HLA_DB_RELEASE:-\(.*\)}$/\1/p' "${REPO_ROOT}/scripts/lib/common.sh")
[ -n "$REL" ] || fail "cannot read HLA_DB_RELEASE from scripts/lib/common.sh"
use_output_hook
cat > "${CASE_WORK}/hook" <<'HOOK'
#!/usr/bin/env bash
. "${CASE_WORK:?}/host-path.sh"
{ printf '%s\n' "$FAKE_DOCKER_VOLUMES"; printf 'arg %s\n' "$@"; } > "${CASE_WORK}/call"
h=$(host_path "${!#}") || exit 31
echo 'IPD-IMGT/HLA fake hla.dat' > "$h"
HOOK
chmod +x "${CASE_WORK}/hook"
cd "$CASE_WORK"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook" FAKE_DOWNLOAD_DIR="$S" GENOME_DIR="$G" OG_DOWNLOADS=shared-dl \
  FETCH_TRIES=1 FETCH_WAIT=0 OG="${REPO_ROOT}/bin/og"

# shellcheck disable=SC2016  # $1 belongs to the inner bash
run_expect 0 hla bash -c '. "$1/scripts/lib/common.sh" && install_data_file hla_dat' _ "$REPO_ROOT"
grep -qxF "${CASE_WORK}/shared-dl:/downloads:ro" "${CASE_WORK}/call" \
  || fail "the container did not get ${CASE_WORK}/shared-dl read-only at /downloads: $(cat "${CASE_WORK}/call")"
grep -qxF "arg /downloads/IPD-IMGT-HLA_${REL}_hla.dat.zip" "${CASE_WORK}/call" \
  || fail "the container was not handed the zip at its /downloads path: $(cat "${CASE_WORK}/call")"
[ -s "${CASE_WORK}/shared-dl/IPD-IMGT-HLA_${REL}_hla.dat.zip" ] || fail "the zip is not kept in shared-dl"
[ ! -e "${G}/downloads" ] || fail "a downloads folder was made in GENOME_DIR too"
grep -q 'IPD-IMGT/HLA' "${G}/hla/IPD-IMGT-HLA_${REL}/hla.dat" || fail "hla.dat is not installed"
echo "PASS a relative OG_DOWNLOADS outside GENOME_DIR is mounted at /downloads and keeps the download"
