#!/usr/bin/env bash
# An otherwise prepared installation fails validation without the og launcher.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

COPY="${CASE_WORK}/repo"
mkdir -p "${COPY}/rust" "${COPY}/bin"
cp -R "${REPO_ROOT}/scripts" "${REPO_ROOT}/assets" "${REPO_ROOT}/versions.env" "$COPY/"
cp "${REPO_ROOT}/rust/og.sh" "${COPY}/rust/"
ln -s "${REPO_ROOT}/bin/.og" "${COPY}/bin/.og"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"

run_expect 1 no-og "${COPY}/scripts/validate-setup.sh"
output_has no-og '\[FAIL\].*bin/og is missing or not executable'
