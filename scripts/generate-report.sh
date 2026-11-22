#!/usr/bin/env bash
# generate-report.sh — Create a plain-text summary report of all pipeline results
# Usage: ./scripts/generate-report.sh <sample_name>
#
# Renders the text report (og report render) from ${SAMPLE}/summary.json,
# the summary og summary collect reads from every step's output. When step 24
# wrote that summary and nothing it reads has changed since (in a run-all.sh
# run, step 24 ran just before), the text report is rendered from it as it is:
# one collect per run, and both reports carry the same values and time.
# Otherwise it collects the summary again first. A result from an earlier run
# whose step was skipped or failed in the latest run-all.sh run is marked
# stale with its date. Run after completing all (or some) pipeline steps.
#
# Then it writes what changed since the previous run of the sample
# (og run delta, docs/delta-report.md): ${SAMPLE}_delta.txt, .html and
# .json, against the newest copy under history/ that holds a summary.json,
# and keeps this run's summary.json and run_manifest.tsv as history/<UTC time>/
# for the next run.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
REPORT="${SAMPLE_DIR}/${SAMPLE}_report.txt"

# Check sample directory exists
if [ ! -d "$SAMPLE_DIR" ]; then
  echo "ERROR: Sample directory not found: ${SAMPLE_DIR}" >&2
  exit 1
fi

if [ ! -f "${SAMPLE_DIR}/run_manifest.tsv" ]; then
  GENOME_DIR="$GENOME_DIR" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" generate-report.sh
fi

# summary_current: 0 when SAMPLE_DIR/summary.json is newer than everything a
# collect reads: every file of the sample folder except what is written from
# the summary or beside it (the reports, the search index, the charts, the
# delta, the logs other than run_status.tsv, and nextflow/, the Nextflow launch
# folder with work/), and bin/ (og and the tables it reads there) and assets/. Anything newer, or a
# summary that is not there, means collect again. A subfolder newer than the
# summary (other than logs/) had a file added or removed since, and so did
# bin/, bin/.og/ or assets/ when it is newer. A section's source that is gone
# (its file or its whole folder removed), or the run status the summary read,
# means the same: collect again.
summary_current() {
  local s="${SAMPLE_DIR}/summary.json" newer src
  [ -s "$s" ] && [ ! -L "$s" ] || return 1
  if grep -q '^  "status_file": true' "$s"; then
    [ -e "${SAMPLE_DIR}/logs/run_status.tsv" ] || return 1
  fi
  # Each section's "source", at the sections' depth of og's one-space indent.
  while IFS= read -r src; do
    case $src in *\\*) return 1 ;; esac
    [ -e "${SAMPLE_DIR}/${src}" ] || return 1
  done < <(sed -n 's/^   "source": "\([^"]*\)",\{0,1\}$/\1/p' "$s")
  newer=$(cd "$SAMPLE_DIR" && find . \( -path ./nextflow -o -path ./facts -o -path ./summary.json \) -prune \
    -o -type d ! -path . ! -path ./logs -newer summary.json -print -quit \
    -o ! -type d -newer summary.json \
    ! \( -path './logs/*' ! -path ./logs/run_status.tsv \) \
    ! -name "${SAMPLE}_report.html" ! -name "${SAMPLE}_report.txt" ! -name "${SAMPLE}_search_index.json" \
    ! -name "${SAMPLE}_charts.json" ! -name "${SAMPLE}_delta.txt" ! -name "${SAMPLE}_delta.html" \
    ! -name "${SAMPLE}_delta.json" -print -quit) || return 1
  [ -z "$newer" ] || return 1
  newer=$(find "${PGP_ROOT}/bin" "${PGP_ROOT}/bin/.og" "${PGP_ROOT}/assets" -maxdepth 0 -newer "$s" -print -quit 2> /dev/null) \
    || return 1
  [ -z "$newer" ] || return 1
  newer=$(find "${PGP_ROOT}/bin" "${PGP_ROOT}/bin/.og" "${PGP_ROOT}/assets" -maxdepth 1 ! -type d -newer "$s" -print -quit 2> /dev/null) \
    || return 1
  [ -z "$newer" ]
}

rm -f "$REPORT"
if summary_current; then
  echo "Rendering from summary.json: step 24 wrote it and nothing it reads has changed since"
  from=(--summary "/genome/${SAMPLE}/summary.json")
else
  # assets/: the ACMG SF rules (acmg_sf.tsv) and the carrier list (carrier_genes.tsv).
  from=(--assets /pgp-assets --json "/genome/${SAMPLE}/summary.json")
fi
# og is the static binary, so the python image runs it.
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" "${OG_IMAGE}" \
  /pgp-bin/og report render \
    "${from[@]}" \
    --sample "$SAMPLE" \
    --sample-dir "/genome/${SAMPLE}" \
    -o "/genome/${SAMPLE}/${SAMPLE}_report.txt"

if [ ! -f "$REPORT" ]; then
  echo "ERROR: the report was not written: ${REPORT}" >&2
  exit 1
fi
cat "$REPORT"
echo ""
echo "Report saved to: ${REPORT}"

# What changed since the previous run, then keep this run as the next baseline.
HISTORY="${SAMPLE_DIR}/history"
# og picks the previous run (--newest-history): the newest history/<stamp>/
# that holds a summary.json, the rule og summary collect's "New since"
# follows too. A newer folder without one, from a copy that stopped half
# way, is passed over. og is the static binary, so the python image runs it.
echo ""
if ! run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${OG_IMAGE}" \
  /pgp-bin/og run delta \
    --sample "$SAMPLE" \
    --summary "/genome/${SAMPLE}/summary.json" \
    --manifest "/genome/${SAMPLE}/run_manifest.tsv" \
    --newest-history "/genome/${SAMPLE}" \
    --out "/genome/${SAMPLE}/${SAMPLE}_delta"; then
  echo "ERROR: the delta against the previous run was not written (og run delta failed)" >&2
  exit 1
fi
# history/<UTC time>[-NN]/: sorts by name in run order.
if [ -f "${SAMPLE_DIR}/summary.json" ]; then
  stamp=$(date -u '+%Y-%m-%dT%H%M%SZ') dest="${HISTORY}/${stamp}" n=1
  while [ -e "$dest" ]; do n=$((n + 1)); dest=$(printf '%s/%s-%02d' "$HISTORY" "$stamp" "$n"); done
  mkdir -p "$dest"
  cp -p "${SAMPLE_DIR}/summary.json" "$dest/"
  [ ! -f "${SAMPLE_DIR}/run_manifest.tsv" ] || cp -p "${SAMPLE_DIR}/run_manifest.tsv" "$dest/"
  echo "Delta saved to: ${SAMPLE_DIR}/${SAMPLE}_delta.txt (this run kept as ${dest#"${SAMPLE_DIR}"/})"
else
  echo "WARNING: no ${SAMPLE_DIR}/summary.json, so this run is not kept in history/" >&2
fi
