#!/usr/bin/env bash
# Exercise the manifest dispatcher and the new indexed-data and BWA installers.
# Stub the older install bodies, which have their own fake-docker cases.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

COPY="${CASE_WORK}/checkout"
mkdir -p "${COPY}/scripts" "${COPY}/rust" "${COPY}/assets"
ln -s "${REPO_ROOT}/scripts/lib" "${COPY}/scripts/lib"
ln -s "${REPO_ROOT}/scripts/ci" "${COPY}/scripts/ci"
ln -s "${REPO_ROOT}/bin" "${COPY}/bin"
cp "${REPO_ROOT}/assets/everything.tsv" "${COPY}/assets/"
cp "${REPO_ROOT}/versions.env" "${COPY}/versions.env"
cat > "${COPY}/rust/og.sh" <<'OG'
#!/usr/bin/env bash
exit 0
OG
chmod +x "${COPY}/rust/og.sh"
SETUP="${COPY}/scripts/setup.sh"

# Keep the dispatch and base setup control flow. The stubs only replace
# installer bodies, and are entered through the same case as a normal run.
legacy='ancestry_panel|ancestry_fine_panel|cyrius|parascopy_data|kir_data|yleaf_data|consumer_sites|ancient_panel|viral_data|exomiser_data|chip_data|strchive_data|megane_kmers|mocha_data|locityper_data|clinical_data|vep_cache_opt_in|annotsv_data|archaic|sample_qc_data|prs|tandem_repeats'
awk -v legacy="$legacy" '
  $0 ~ "^install_(" legacy ")\\(\\) \\{" {
    print
    if ($0 ~ /^install_sample_qc_data/) print "  if $SAMPLE_QC_ONLY; then echo legacy --sample-qc-data >> \"$FAKE_DOCKER_LOG\"; fi"
    else print "  if [ -n \"$OPT_IN\" ]; then printf \"legacy --%s\\n\" \"$OPT_IN\" >> \"$FAKE_DOCKER_LOG\"; fi"

    if ($0 ~ /^install_strchive_data/) print "  if [ \"$OPT_IN\" = strchive-data ]; then [ \"${STR_EXTRA_LOCI:-}\" = yes ] || return 1; fi"
    print "  return 0"
    next
  }
  { print }
' "${SCRIPTS}/setup.sh" > "$SETUP"
chmod +x "$SETUP"

# Pin the exact bytes served by fake wget. Keep the real fetch and digest
# code, including the SpliceAI range checks. All data files stay tiny.
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
export FAKE_DOWNLOAD_DIR="${CASE_WORK}/downloads"
export FETCH_TRIES=1 FETCH_WAIT=0
mkdir -p "$FAKE_DOWNLOAD_DIR"
seed_pins() {
  local base=$1 key file kind pin dest keys=()
  read -r -a keys <<< "$2"
  for key in "${keys[@]}"; do
    read -r file kind _ <<< "${!key}"
    dest="${FAKE_DOWNLOAD_DIR}/${file}"
    fake_body "${base}/${file}" > "${dest}.seed"
    mv "${dest}.seed" "$dest"
    case "$kind" in
      sha256) pin=$(sha256sum "$dest" | awk '{print $1}') ;;
      md5) pin=$(fake_md5 < "$dest") ;;
      gzip-ends)
        pin=$(sha256sum "$dest" | awk '{print $1}')
        pin="$(wc -c < "$dest"):${pin}:${pin}" ;;
      *) fail "unknown fixture pin kind $kind" ;;
    esac
    printf '%s="%s %s %s"\n' "$key" "$file" "$kind" "$pin"
  done
}
for spec in CADD SPLICEAI_MASKED SOMATIC_DATA CNVPYTOR_DATA; do
  base_var="${spec}_BASE_URL" files_var="${spec}_FILES"
  seed_pins "${!base_var}" "${!files_var}" >> "${COPY}/versions.env"
done
for spec in REVEL ALPHAMISSENSE; do
  url_var="${spec}_URL"
  printf '%s_SHA256="%s"\n' "$spec" "$(fake_body "${!url_var}" | sha256sum | awk '{print $1}')" >> "${COPY}/versions.env"
done

G="${CASE_WORK}/data"
seed_reference "$G"
seed_clinvar "$G"
printf 'dict\n' > "${G}/reference/GRCh38_no_alt_analysis_set.dict"
mkdir -p "${G}/gnomad_mito"
printf 'fixture\n' > "${G}/gnomad_mito/${GNOMAD_MITO_URL##*/}"
printf 'GNOMAD_MITO_SHA256="%s"\n' "$(sha256sum "${G}/gnomad_mito/${GNOMAD_MITO_URL##*/}" | awk '{print $1}')" >> "${COPY}/versions.env"

# Reuse the paths the unmodified base setup checks. Each file is a preinstalled fixture.
# shellcheck disable=SC2016  # the inner Bash resolves the installed paths
while read -r file; do
  mkdir -p "${file%/*}"
  printf 'fixture\n' > "$file"
done < <(GENOME_DIR="$G" bash -c '. "$1"; for name in $DATA_FILES; do data_file "$name" || :; done' _ "${COPY}/scripts/lib/common.sh")

ln -s "${REPO_ROOT}/assets/reference" "${COPY}/assets/reference"

export SETUP_FIXTURE_DIR="$G"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tool"
cat > "$FAKE_DOCKER_RUN_HOOK" <<'TOOL'
#!/usr/bin/env bash
set -eu
shift
host() { printf '%s%s' "$SETUP_FIXTURE_DIR" "${1#/genome}"; }
case "$1" in
  /pgp-bin/og)
    [ "$2/$3" = setup/revel ]
    printf '#chr\tpos\tref\talt\tREVEL\nchr1\t10\tA\tT\t0.1\n' > "$(host "$5")" ;;
  bgzip)
    file=$(host "$3")
    gzip -n "$file" ;;
  tabix) printf 'index\n' > "$(host "${@: -1}").tbi" ;;
  bwa)
    [ "$2/$3" = index/-p ]
    prefix=$(host "$4")
    for ext in amb ann bwt pac sa; do printf 'index\n' > "${prefix}.${ext}"; done ;;
esac
TOOL
chmod +x "$FAKE_DOCKER_RUN_HOOK"


# Replace Git's external download with a tagged checkout of the pinned commit.
mkdir -p "${CASE_WORK}/test-bin"
export FAKE_PYPGX_COMMIT="$PYPGX_BUNDLE_COMMIT"
cat > "${CASE_WORK}/test-bin/git" <<'GIT'
#!/usr/bin/env bash
set -eu
case "$1" in
  clone)
    dest=${@: -1}
    mkdir -p "$dest"
    printf '%s\n' "$FAKE_PYPGX_COMMIT" > "${dest}/commit" ;;
  -C)
    [ "$3 $4" = 'rev-parse HEAD' ]
    cat "${2}/commit" ;;
  *) exit 1 ;;
esac
GIT
chmod +x "${CASE_WORK}/test-bin/git"
export PATH="${CASE_WORK}/test-bin:$PATH"

# Use the actual fake wget path, including its resume and failure handling.
hide_commands curl
run_expect 0 plain "$SETUP" "$G"
[ ! -e "${G}/annotations" ] && [ ! -e "${G}/somatic" ] || fail 'plain setup installed opt-in data'
[ ! -e "${G}/reference/GRCh38_no_alt_analysis_set.fasta.sa" ] || fail 'plain setup built a classic BWA index'
output_lacks plain '^  (CADD|SpliceAI precomputed scores|REVEL|AlphaMissense|Cyrius|PGS000014):'

# Force one download to fail. The other score files and later installers
# still finish, and retrying the same directory only fetches the failed pair.
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_FAIL='whole_genome_SNVs\.tsv\.gz$' run_expect 3 partial "$SETUP" --everything "$G"
sed -n 's/^=== Everything installer: \(.*\) ===$/\1/p' "${CASE_WORK}/partial.out" > "${CASE_WORK}/ran"
awk -F'\t' '$1=="setup" {print $2}' "${COPY}/assets/everything.tsv" > "${CASE_WORK}/expected"
diff -u "${CASE_WORK}/expected" "${CASE_WORK}/ran"
sed -n 's/^legacy //p' "$FAKE_DOCKER_LOG" > "${CASE_WORK}/legacy-ran"
awk '$0 !~ /^--(cnvpytor-data|pypgx-data|annotation-scores|somatic-data|bwa-index)$/' "${CASE_WORK}/expected" > "${CASE_WORK}/legacy-expected"
diff -u "${CASE_WORK}/legacy-expected" "${CASE_WORK}/legacy-ran"

for name in CADD 'SpliceAI precomputed scores' REVEL AlphaMissense Cyrius PGS000014 PGS000016 PGS000017 PGS000027; do
  [ "$(grep -c "^  ${name}:" "${CASE_WORK}/partial.out")" -eq 1 ] || fail "expected one licence line for $name"
done
output_has partial '^  - --annotation-scores '
[ "$(grep -c '^  - ' "${CASE_WORK}/partial.out")" -eq 1 ] || fail 'unexpected missing items'
output_lacks partial '^Everything setup complete!$'
[ ! -e "${G}/annotations/whole_genome_SNVs.tsv.gz.tbi" ] || fail 'index published without its data'
for f in gnomad.genomes.r4.0.indel.tsv.gz spliceai_scores.masked.snv.hg38.vcf.gz spliceai_scores.masked.indel.hg38.vcf.gz revel_grch38.tsv.gz AlphaMissense_hg38.tsv.gz; do
  [ -s "${G}/annotations/${f}" ] && [ -s "${G}/annotations/${f}.tbi" ] || fail "missing annotation pair $f"
done
for f in af-only-gnomad 1000g_pon small_exac_common_3; do
  [ -s "${G}/somatic/${f}.hg38.vcf.gz" ] && [ -s "${G}/somatic/${f}.hg38.vcf.gz.tbi" ] || fail "missing somatic pair $f"
done
for ext in amb ann bwt pac sa; do
  [ -s "${G}/reference/GRCh38_no_alt_analysis_set.fasta.${ext}" ] || fail "missing BWA .$ext"
done
[ -s "${G}/reference/cnvpytor/gc_hg38.pytor" ] || fail "CNVpytor did not install"
[ -s "${G}/reference/pypgx-bundle/commit" ] || fail "PyPGx did not install"
docker_log_has '^wget .*/revel-v1.3_all_chromosomes.zip .*\.part$' 'REVEL did not use resumable fake wget'

: > "$FAKE_DOCKER_LOG"
run_expect 0 retry "$SETUP" --everything "$G"
output_has retry '^Everything setup complete!$'
[ "$(grep -c '^wget ' "$FAKE_DOCKER_LOG")" -eq 2 ] || fail 'retry downloaded more than the missing CADD pair'
[ -s "${G}/annotations/whole_genome_SNVs.tsv.gz.tbi" ] || fail 'retry did not finish CADD'
output_has retry 'Classic BWA index already present'

# A manifest row with no handler must fail and be named in the summary.
printf 'setup\t--missing-installer\tvcfanno\t-\tplanted test row\n' >> "${COPY}/assets/everything.tsv"
run_expect 3 unknown "$SETUP" --everything "$G"
output_has unknown '^ERROR: unknown setup flag: --missing-installer$'
output_has unknown '^  - --missing-installer \(planted test row\)$'
