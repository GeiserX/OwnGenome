#!/usr/bin/env bash
# ownhealth-facts.sh — Write the file OwnHealth imports from a finished run
# Usage: ./scripts/ownhealth-facts.sh <sample_name>
#
# Writes ${SAMPLE}/facts/${SAMPLE}_genome_facts.json (og facts,
# docs/ownhealth.md) from the sample folder's summary.json and the outputs
# beside it. run-all.sh runs it last, after the reports. The reference is the
# one the run used: the VRS identifiers in the file are worked out on it.
# What og facts leaves out, and why, is in the log.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

# og is the static binary, so the python image runs it. Without summary.json
# (step 24 or generate-report.sh writes it) or the reference and its .fai, og
# facts says which file is missing and exits 1. There is no host-side check:
# the fake-docker cases run run-all.sh with no reference on disk. An earlier
# run's file is removed first, so a failed run leaves no file to import.
rm -f "${GENOME_DIR}/${SAMPLE}/facts/${SAMPLE}_genome_facts.json"
echo "[ownhealth-facts] writing ${GENOME_DIR}/${SAMPLE}/facts/${SAMPLE}_genome_facts.json from ${GENOME_DIR}/${SAMPLE}/summary.json and ${REF_FASTA}"
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" "${PYTHON_IMAGE}" \
  /pgp-bin/og facts \
    --sample "$SAMPLE" \
    --sample-dir "/genome/${SAMPLE}" \
    --reference "$REF_FASTA_C" \
    --assets /pgp-assets
