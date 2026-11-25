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
# code, including the SpliceAI range checks. The largest fixture is about two MiB.
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
mkdir -p "${G}/roh_af"
printf 'fixture\n' > "${G}/roh_af/${ROH_AF_PANEL_URL##*/}"
printf 'ROH_AF_PANEL_SHA256="%s"\n' "$(sha256sum "${G}/roh_af/${ROH_AF_PANEL_URL##*/}" | awk '{print $1}')" >> "${COPY}/versions.env"

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
    # Git accepts the checkout created by this setup process. The existing
    # bundle below requires an explicit safe.directory override.
    [[ "$2" = */.pypgx.*/bundle ]] && [ "$3 $4" = 'rev-parse HEAD' ] || exit 1
    cat "${2}/commit" ;;
  -c)
    [ "$2" = "safe.directory=$4" ] && [ "$3 $5 $6" = '-C rev-parse HEAD' ] || exit 1
    cat "${4}/commit" ;;
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

# Force the score preflight to run one KiB short. The same free space
# suffices after a one-KiB partial download is present.
cadd="${G}/annotations/whole_genome_SNVs.tsv.gz"
rm -f "$cadd" "${cadd}.tbi"
printf 'CADD_SNV_BYTES=2048\n' >> "${COPY}/versions.env"
export NORMAL_DF
NORMAL_DF=$(command -v df)
cat > "${CASE_WORK}/test-bin/df" <<'DF'
#!/usr/bin/env bash
if [ -n "${FAKE_SCORE_FREE_KIB:-}" ]; then
  printf 'Filesystem 1024-blocks Used Available Capacity Mounted on\n'
  printf 'fakefs 2147483648 0 %s 0%% /\n' "$FAKE_SCORE_FREE_KIB"
else
  exec "$NORMAL_DF" "$@"
fi
DF
chmod +x "${CASE_WORK}/test-bin/df"
export PATH="${CASE_WORK}/test-bin:$PATH"
: > "$FAKE_DOCKER_LOG"
space_margin=$((5 * 1048576))
space_with_margin=$((space_margin + 1))
FAKE_SCORE_FREE_KIB=$space_with_margin run_expect 1 low-space "$SETUP" --annotation-scores "$G"
output_has low-space '^\[WARN\] Skipping annotation scores: needs 2048 more bytes for CADD and SpliceAI plus a 5 GiB margin\.$'
output_has low-space '^  Install later: .*setup\.sh --annotation-scores '
[ ! -s "$FAKE_DOCKER_LOG" ] || fail 'low-space preflight started a download'
head -c 1024 /dev/zero > "${cadd}.part"
FAKE_SCORE_FREE_KIB=$space_with_margin run_expect 0 resume-with-space "$SETUP" --annotation-scores "$G"
[ -s "$cadd" ] && [ -s "${cadd}.tbi" ] || fail 'resumed download did not publish the CADD pair'
[ ! -e "${cadd}.part" ] || fail 'resumed download left a .part'

# With room for the margin alone, a size-matching score needs only its index.
# Log calls to the real MD5 tool so a duplicate content check fails the case.
rm -f "${cadd}.tbi"
printf 'CADD_SNV_BYTES="%s"\n' "$(stat -c %s "$cadd")" >> "${COPY}/versions.env"
export NORMAL_MD5SUM
NORMAL_MD5SUM=$(command -v md5sum)
cat > "${CASE_WORK}/test-bin/md5sum" <<'MD5SUM'
#!/usr/bin/env bash
printf 'md5sum %s\n' "$@" >> "$FAKE_DOCKER_LOG"
exec "$NORMAL_MD5SUM" "$@"
MD5SUM
chmod +x "${CASE_WORK}/test-bin/md5sum"
: > "$FAKE_DOCKER_LOG"
FAKE_SCORE_FREE_KIB=$space_margin run_expect 0 retry-index-only "$SETUP" --annotation-scores "$G"
[ -s "${cadd}.tbi" ] || fail 'index retry did not publish the CADD index'
docker_log_has '^wget .*/whole_genome_SNVs\.tsv\.gz\.tbi ' 'index retry did not fetch the CADD index'
[ "$(grep -c '^wget ' "$FAKE_DOCKER_LOG")" -eq 1 ] || fail 'index retry fetched more than the CADD index'
cadd_checks=$(awk -v path="$cadd" '$0 == "md5sum " path { count++ } END { print count + 0 }' "$FAKE_DOCKER_LOG")
[ "$cadd_checks" -eq 1 ] || fail "index retry checked CADD data $cadd_checks times, expected 1"
rm -f "${CASE_WORK}/test-bin/md5sum"

# The outer gzip members each span more than one MiB. Recompress each with
# one input byte flipped so it can replace just one end of the stream.
fixture="${CASE_WORK}/spliceai-gzip"
mkdir "$fixture"
for member in first last; do
  input="${fixture}/${member}.raw"
  head -c 1048576 /dev/urandom > "$input"
  gzip -n < "$input" > "${fixture}/${member}.gz"
  member_size=$(stat -c %s "${fixture}/${member}.gz")
  [ "$member_size" -gt 1048576 ] || fail 'gzip member must span more than one MiB'
  byte=$(od -An -tu1 -j 524288 -N 1 "$input")
  printf '%b' "\\$(printf '%03o' "$((byte ^ 1))")" |
    dd of="$input" bs=1 seek=524288 conv=notrunc status=none
  gzip -n < "$input" > "${fixture}/${member}-alt.gz"
  [ "$(stat -c %s "${fixture}/${member}-alt.gz")" -eq "$member_size" ] || fail "$member gzip variant changed size"
done
printf 'SpliceAI test payload\n' | gzip -n > "${fixture}/middle.gz"
cat "${fixture}/first.gz" "${fixture}/middle.gz" "${fixture}/last.gz" > "${fixture}/good.gz"
size=$(stat -c %s "${fixture}/good.gz")
first=$(head -c 1048576 "${fixture}/good.gz" | sha256sum | awk '{print $1}')
last=$(tail -c 1048576 "${fixture}/good.gz" | sha256sum | awk '{print $1}')
read -r snv _ <<< "$SPLICEAI_MASKED_SNV_PIN"
printf 'SPLICEAI_MASKED_SNV_PIN="%s gzip-ends %s:%s:%s"\n' "$snv" "$size" "$first" "$last" >> "${COPY}/versions.env"

cat "${fixture}/first-alt.gz" "${fixture}/middle.gz" "${fixture}/last.gz" > "${fixture}/first-end.gz"
cat "${fixture}/first.gz" "${fixture}/middle.gz" "${fixture}/last-alt.gz" > "${fixture}/last-end.gz"
gzip -t "${fixture}/first-end.gz" "${fixture}/last-end.gz" || fail 'gzip end variants must be valid'

# Corrupt the middle CRC or drop its member to isolate the gzip and size checks.
cp "${fixture}/good.gz" "${fixture}/bad-crc.gz"
crc_offset=$(( $(stat -c %s "${fixture}/first.gz") + $(stat -c %s "${fixture}/middle.gz") - 8 ))
crc_byte=$(od -An -tu1 -j "$crc_offset" -N 1 "${fixture}/good.gz")
printf '%b' "\\$(printf '%03o' "$((crc_byte ^ 1))")" |
  dd of="${fixture}/bad-crc.gz" bs=1 seek="$crc_offset" conv=notrunc status=none
cat "${fixture}/first.gz" "${fixture}/last.gz" > "${fixture}/missing-member.gz"
head -c "$((size - 1))" "${fixture}/good.gz" > "${fixture}/truncated.gz"

masked="${G}/annotations/${snv}"
for body in bad-crc missing-member truncated first-end last-end good; do
  rm -f "$masked" "${masked}.tbi"
  cp "${fixture}/${body}.gz" "${FAKE_DOWNLOAD_DIR}/${snv}"
  if [ "$body" = good ]; then
    run_expect 0 "gzip-${body}" "$SETUP" --annotation-scores "$G"
    cmp "${fixture}/good.gz" "$masked"
    [ -s "${masked}.tbi" ] || fail 'valid gzip is missing its index'
  else
    run_expect 1 "gzip-${body}" "$SETUP" --annotation-scores "$G"
    output_has "gzip-${body}" 'failed its pinned size, SHA-256 range or gzip check'
    [ ! -e "$masked" ] && [ ! -e "${masked}.tbi" ] || fail "published the invalid SpliceAI pair ($body)"
  fi
  [ ! -e "${masked}.part" ] || fail "left a SpliceAI .part ($body)"
done

# A complete raw SNV pair replaces only masked SNVs. Without a raw indel
# index, setup still installs the masked indel pair.
read -r indel _ <<< "$SPLICEAI_MASKED_INDEL_PIN"
masked_indel="${G}/annotations/${indel}"
raw_snv="${masked/.masked./.raw.}"
raw_indel="${masked_indel/.masked./.raw.}"
mv "$masked" "$raw_snv"
mv "${masked}.tbi" "${raw_snv}.tbi"
mv "$masked_indel" "$raw_indel"
rm -f "${masked_indel}.tbi"
FAKE_SCORE_FREE_KIB=$space_with_margin run_expect 0 prefer-raw "$SETUP" --annotation-scores "$G"
output_has prefer-raw '^\[OK\] SpliceAI raw scores already present: .*spliceai_scores\.raw\.snv\.hg38\.vcf\.gz and .*spliceai_scores\.raw\.snv\.hg38\.vcf\.gz\.tbi$'
[ ! -e "$masked" ] && [ ! -e "${masked}.tbi" ] || fail 'downloaded masked SNVs despite a complete raw pair'
[ -s "$masked_indel" ] && [ -s "${masked_indel}.tbi" ] || fail 'incomplete raw indels blocked the masked pair'


# A manifest row with no handler must fail and be named in the summary.
printf 'setup\t--missing-installer\tvcfanno\t-\tplanted test row\n' >> "${COPY}/assets/everything.tsv"
run_expect 3 unknown "$SETUP" --everything "$G"
output_has unknown '^ERROR: unknown setup flag: --missing-installer$'
output_has unknown '^  - --missing-installer \(planted test row\)$'
