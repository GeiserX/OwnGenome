#!/usr/bin/env bash
# bin/write_manifest.sh's run pipeline_version row: the OwnGenome release, the
# version in the manifest { } block of nextflow.config, which both reports
# print beside the pipeline commit.
#   1. The repository's write_manifest.sh writes the row, with the version
#      nextflow.config declares.
#   2. A copy of the tree whose nextflow.config declares 9.8.7 gets 9.8.7: the
#      row is read, not written from a constant.
#   3. Red control: a copy of write_manifest.sh without the row line writes no
#      row, so the check in 1 can fail.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
mkdir -p "${GENOME_DIR}/sample1"
M="${GENOME_DIR}/sample1/run_manifest.tsv"
release() { awk -F'\t' '$1 == "run" && $2 == "pipeline_version" {print $3}' "$M"; }
DECLARED=$(sed -nE "/^manifest[[:space:]]*\{/,/^\}/ s/^[[:space:]]*version[[:space:]]*=[[:space:]]*'([^']+)'.*/\1/p" \
  "${REPO_ROOT}/nextflow.config")
[ -n "$DECLARED" ] || fail "could not read the manifest version of nextflow.config"

# 1. the repository's script
run_expect 0 repo env bash "${REPO_ROOT}/bin/write_manifest.sh" sample1 test
[ "$(release)" = "$DECLARED" ] || fail "run pipeline_version reads '$(release)', nextflow.config declares '${DECLARED}'"
[ "$(grep -c $'^run\tpipeline_version\t' "$M")" = 1 ] || fail "run_manifest.tsv has not exactly one pipeline_version row"

# A copy of the files write_manifest.sh reads: bin/, scripts/lib/, versions.env, nextflow.config
copy_tree() {
  rm -rf "${CASE_WORK}/$1"
  mkdir -p "${CASE_WORK}/$1/scripts"
  cp -R "${REPO_ROOT}/bin" "${CASE_WORK}/$1/"
  cp -R "${REPO_ROOT}/scripts/lib" "${CASE_WORK}/$1/scripts/"
  cp "${REPO_ROOT}/versions.env" "${REPO_ROOT}/nextflow.config" "${CASE_WORK}/$1/"
}

# 2. another declared version
copy_tree other
sed -i.bak -E "/^manifest[[:space:]]*\{/,/^\}/ s/^([[:space:]]*version[[:space:]]*=[[:space:]]*)'[^']+'/\1'9.8.7'/" \
  "${CASE_WORK}/other/nextflow.config"
grep -q "version         = '9.8.7'" "${CASE_WORK}/other/nextflow.config" || fail "could not plant 9.8.7 in the copy"
run_expect 0 other env bash "${CASE_WORK}/other/bin/write_manifest.sh" sample1 test
[ "$(release)" = 9.8.7 ] || fail "with nextflow.config at 9.8.7 the row reads '$(release)'"

# 3. red control: the script without the row line
copy_tree norow
sed -i.bak '/row run pipeline_version/d' "${CASE_WORK}/norow/bin/write_manifest.sh"
! grep -q 'pipeline_version "' "${CASE_WORK}/norow/bin/write_manifest.sh" || fail "the control copy still writes the row"
run_expect 0 norow env bash "${CASE_WORK}/norow/bin/write_manifest.sh" sample1 test
[ -z "$(release)" ] || fail "control: the copy without the row line still wrote '$(release)'"
echo "run pipeline_version is nextflow.config's manifest version (${DECLARED}), and the control without it has none."
