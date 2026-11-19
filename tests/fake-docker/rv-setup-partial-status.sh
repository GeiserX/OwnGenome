#!/usr/bin/env bash
# setup.sh says "Setup complete!" and exits 0 only when it installed
# everything. When an install fails or is skipped, it prints "Setup finished:
# core ready, N items missing", lists each item with the command that retries
# it, and exits 3.
#
# The fake downloads fail their pinned checksums and og is not built for both
# Linux arches here, so a plain fake setup always misses something. To reach
# the complete case, setup.sh runs from a copy of the checkout (not a git
# checkout, so og.sh has no tree to compare) with two placeholder og binaries,
# and GENOME_DIR holds what the pinned downloads would have installed. A fake
# sha256sum gives the two files setup.sh checks by digest their pinned values.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
seed_clinical_downloads
export FETCH_TRIES=1 FETCH_WAIT=0

COPY="${CASE_WORK}/checkout"
mkdir -p "${COPY}/rust" "${COPY}/bin/.og"
cp -R "${REPO_ROOT}/scripts" "${REPO_ROOT}/assets" "${REPO_ROOT}/versions.env" "$COPY/"
cp "${REPO_ROOT}/rust/og.sh" "${COPY}/rust/"
for f in "${REPO_ROOT}"/bin/*; do
  [ -f "$f" ] && cp "$f" "${COPY}/bin/"
done
for a in x86_64 aarch64; do
  printf '#!/bin/sh\nexit 0\n' > "${COPY}/bin/.og/og-linux-${a}"
  chmod +x "${COPY}/bin/.og/og-linux-${a}"
done
SETUP="${COPY}/scripts/setup.sh"

# A fake sha256sum: the pinned digest for the tandem-repeat map and gnomAD's
# chrM sites, the real one for every other file.
real_sha=$(command -v sha256sum || true)
mkdir -p "${CASE_WORK}/sha-bin"
cat > "${CASE_WORK}/sha-bin/sha256sum" <<SHA
#!/usr/bin/env bash
case "\${1##*/}" in
  ${TANDEM_REPEATS_URL##*/}) echo "${TANDEM_REPEATS_SHA256}  \$1" ;;
  ${GNOMAD_MITO_URL##*/}) echo "${GNOMAD_MITO_SHA256}  \$1" ;;
  *) if [ -n "${real_sha}" ]; then exec "${real_sha}" "\$@"; else exec shasum -a 256 "\$@"; fi ;;
esac
SHA
chmod +x "${CASE_WORK}/sha-bin/sha256sum"
export PATH="${CASE_WORK}/sha-bin:${PATH}"

# seed_installed G: what phases 4 and 5 install, as setup.sh finds it when
# it is already there. AnnotSV's bundle is left out with "no-annotsv".
seed_installed() {
  local g=$1 r="$1/reference" a="$1/annotsv_annotations" u
  mkdir -p "${r}/somalier" "${r}/verifybamid2" "${r}/strchive-${STRCHIVE_VERSION}" "${g}/gnomad_mito"
  printf 'x\n' > "${r}/GRCh38_no_alt_analysis_set.dict"
  # The pinned data files, by the paths common.sh's data_file gives them.
  # shellcheck disable=SC2016  # expanded by the inner bash
  while read -r u; do
    mkdir -p "${u%/*}"
    printf 'x\n' > "$u"
  done < <(GENOME_DIR=$g bash -c '. "$1"; for n in $DATA_FILES; do data_file "$n" || true; done' _ "${COPY}/scripts/lib/common.sh")
  for u in somalier/sites.hg38.vcf.gz \
           verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.UD verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.mu \
           verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.bed "${TANDEM_REPEATS_URL##*/}" \
           "$(basename "${TANDEM_REPEATS_URL}" .gz)" "strchive-${STRCHIVE_VERSION}/catalog.json" \
           "strchive-${STRCHIVE_VERSION}/str_ranges.tsv"; do
    printf 'x\n' > "${r}/${u}"
  done
  printf 'x\n' > "${g}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
  printf 'release=%s eh=%s loci=%s converter=og-str-catalog-1 extra=off\n' \
    "$STRCHIVE_VERSION" "$STRCHIVE_EH_SHA256" "$STRCHIVE_LOCI_SHA256" > "${r}/strchive-${STRCHIVE_VERSION}/catalog.stamp"
  [ "${2:-}" = no-annotsv ] && return 0
  mkdir -p "${a}/Annotations_Human/Genes/GRCh38" "${a}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38"
  printf 'gnomAD-SV %s sha256 %s\n' "$GNOMAD_SV_VERSION" "$GNOMAD_SV_SHA256" \
    > "${a}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38/gnomad-sv.installed"
  printf 'NM_000001.1\n' > "${a}/mane_select_refseq.txt"
  printf 'MANE %s sha256 %s\n' "$MANE_RELEASE" "$MANE_SUMMARY_SHA256" > "${a}/mane_select_refseq.release"
}

# 1. Everything installs: the banner, no missing list, exit 0.
G="${CASE_WORK}/complete"
seed_installed "$G"
run_expect 0 complete "$SETUP" "$G"
output_has complete '^  Setup complete!$'
output_lacks complete 'Setup finished|Missing \(|exits with status 3'

# 2. The AnnotSV download fails: no banner, AnnotSV and its retry, exit 3.
G="${CASE_WORK}/annotsv"
seed_installed "$G" no-annotsv
FAKE_DOWNLOAD_FAIL='lbgi\.fr' run_expect 3 annotsv "$SETUP" "$G"
output_lacks annotsv 'Setup complete'
output_has annotsv '^  Setup finished: core ready, 1 item missing$'
output_has annotsv '^  - AnnotSV data: the bundle, gnomAD-SV .* \(step 05\)$'
output_has annotsv "^      retry: .*/scripts/setup\\.sh --annotsv-data \"${G}\"$"
output_has annotsv 'Setup exits with status 3: a partial setup'
[ "$(grep -c '^  - ' "${CASE_WORK}/annotsv.out")" -eq 1 ] || fail 'annotsv lists more than the AnnotSV data'

# 3. Too little disk for the two VEP caches: each skipped cache is listed
# with its own retry, the bundle still installs, exit 3.
mkdir -p "${CASE_WORK}/low-bin"
cat > "${CASE_WORK}/low-bin/df" <<'DF'
#!/usr/bin/env bash
echo 'Filesystem 1024-blocks Used Available Capacity Mounted on'
echo 'fakefs 41943040 0 41943040 0% /'
DF
chmod +x "${CASE_WORK}/low-bin/df"
G="${CASE_WORK}/low"
seed_installed "$G"
PATH="${CASE_WORK}/low-bin:${PATH}" run_expect 3 low "$SETUP" "$G"
output_lacks low 'Setup complete'
output_has low '^  Setup finished: core ready, 2 items missing$'
output_has low "^  - VEP ${VEP_CACHE_RELEASE} cache \\(steps 13, 30, 23 and 31\\)$"
output_has low "^      retry: .*/scripts/setup\\.sh --vep-cache \"${G}\"$"
output_has low "^  - VEP ${PCGR_VEP_CACHE_RELEASE} cache \\(step 17, CPSR\\)$"
output_has low "^      retry: .*/scripts/setup\\.sh --clinical-data \"${G}\"$"
output_lacks low '^  - PCGR'
[ -d "${G}/pcgr_data/${PCGR_DATA_BUNDLE}/data" ] || fail 'the bundle was not installed beside the skipped caches'

# 4. --no-clinical-data is a choice, not a gap: complete, exit 0.
G="${CASE_WORK}/opt-out"
seed_installed "$G"
run_expect 0 opt_out "$SETUP" --no-clinical-data "$G"
output_has opt_out '^  Setup complete!$'

# 5. A plain fake setup misses several items: each has a retry line.
G="${CASE_WORK}/plain"
mkdir -p "$G"
run_expect 3 plain "$SETUP" "$G"
output_lacks plain 'Setup complete'
items=$(grep -c '^  - ' "${CASE_WORK}/plain.out")
retries=$(grep -c '^      retry: ' "${CASE_WORK}/plain.out")
[ "$items" -ge 5 ] && [ "$items" -eq "$retries" ] || fail "plain lists ${items} items and ${retries} retry lines"
output_has plain "^  Setup finished: core ready, ${items} items missing$"
output_has plain "retry: .*/scripts/setup\\.sh --sample-qc-data \"${G}\"$"
output_has plain "retry: .*/scripts/setup\\.sh --strchive-data \"${G}\"$"
