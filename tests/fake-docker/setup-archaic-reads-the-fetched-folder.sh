#!/usr/bin/env bash
# setup.sh --archaic reads each genome's chromosome from the folder its fetch
# wrote it to: the archaic folder of the downloads folder. The read used to
# name reference/archaic/<version>/download, where nothing is written now.
# A checkout copy pins the served fixtures (versions.env is not read from the
# environment); the hook records the arguments of the genotype read and stops
# the install there.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

COPY="${CASE_WORK}/checkout" S="${CASE_WORK}/served" G="${CASE_WORK}/genome"
mkdir -p "${COPY}/scripts" "${COPY}/rust" "$S"
for d in scripts/lib scripts/ci bin assets; do ln -s "${REPO_ROOT}/${d}" "${COPY}/${d}"; done
cp "${SCRIPTS}/setup.sh" "${COPY}/scripts/setup.sh"
printf '#!/usr/bin/env bash\nexit 0\n' > "${COPY}/rust/og.sh"
chmod +x "${COPY}/rust/og.sh"
sha() { sha256sum "$1" | cut -d' ' -f1; }
printf 'chain\n' | gzip -c > "${S}/hg19ToHg38.over.chain.gz"
printf 'chr22\t0\t100\n' > "${S}/20160622.allChr.mask.bed"
for i in $(seq 1 120); do printf '0\tchr22\t%d\t%d\n' "$((i * 1000))" "$((i * 1000 + 10))"; done | gzip -c > "${S}/genomicSuperDups.txt.gz"
printf '##fileformat=VCFv4.2\n' | gzip -c > "${S}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"
printf '##fileformat=VCFv4.2\n' | gzip -c > "${S}/chr22_mq25_mapab100.vcf.gz"
printf 'mask\n' | gzip -c > "${S}/chr22_mask.bed.gz"
printf '%s  chr22_mq25_mapab100.vcf.gz\n' "$(md5sum < "${S}/chr22_mq25_mapab100.vcf.gz" | cut -d' ' -f1)" > "${S}/md5sums.txt"
{ grep -v '^ARCHAIC_' "${REPO_ROOT}/versions.env"
  grep '^ARCHAIC_DATA_VERSION=' "${REPO_ROOT}/versions.env"
  for g in ALTAI VINDIJA DENISOVA; do
    echo "ARCHAIC_${g}_MD5SUMS_SHA256=\"$(sha "${S}/md5sums.txt")\""
    echo "ARCHAIC_FILTERBED_SHA256_${g}=\"22:$(sha "${S}/chr22_mask.bed.gz")\""
  done
  echo "ARCHAIC_1KG_SITES_MD5=\"$(md5sum < "${S}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz" | cut -d' ' -f1)\""
  echo "ARCHAIC_CHAIN_SHA256=\"$(sha "${S}/hg19ToHg38.over.chain.gz")\""
  echo "ARCHAIC_STRICT_MASK_SHA256=\"$(sha "${S}/20160622.allChr.mask.bed")\""
  echo "ARCHAIC_SEGDUPS_SHA256=\"$(sha "${S}/genomicSuperDups.txt.gz")\""
} > "${COPY}/versions.env"
V=$(sed -n 's/^ARCHAIC_DATA_VERSION="\(.*\)"/\1/p' "${COPY}/versions.env")

seed_reference "$G"
printf '@HD\tVN:1.6\n' > "${G}/reference/GRCh38_no_alt_analysis_set.dict"
printf '22\t0\t50000000\n' > "${CASE_WORK}/regions.bed"
use_output_hook
cat > "${CASE_WORK}/hook" <<'HOOK'
#!/usr/bin/env bash
. "${CASE_WORK:?}/host-path.sh"
case "$*" in
  *nonref.pos*) while [ "$1" != _ ]; do shift; done
    host_path "$2" > "${CASE_WORK}/read-from"; exit 1 ;;
  *indels.tsv*) while [ "$1" != _ ]; do shift; done
    h=$(host_path "$4") && : > "${h}/indels.tsv" && : > "${h}/modern.tsv" ;;
esac
exit 0
HOOK
chmod +x "${CASE_WORK}/hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook" FAKE_DOWNLOAD_DIR="$S" ARCHAIC_REGIONS="${CASE_WORK}/regions.bed" \
  SETUP_DISK_CHECK=false FETCH_TRIES=1 FETCH_WAIT=0 OG="${REPO_ROOT}/bin/og"

run_expect 1 archaic "${COPY}/scripts/setup.sh" --archaic "$G"
output_has archaic '\[WARN\] Could not read the archaic genotypes of chr22'
want="${G}/downloads/archaic-${V}"
[ -s "${want}/Altai.chr22.vcf.gz" ] || fail "the fetch did not write ${want}/Altai.chr22.vcf.gz"
[ "$(cat "${CASE_WORK}/read-from" 2>/dev/null)" = "$want" ] \
  || fail "the genotype read names $(cat "${CASE_WORK}/read-from" 2>/dev/null || echo nothing), not ${want}, the folder the fetch wrote"
echo "PASS setup.sh --archaic reads the archaic genotypes from the folder their fetch wrote"
