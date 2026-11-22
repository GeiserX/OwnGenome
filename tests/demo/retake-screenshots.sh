#!/usr/bin/env bash
# retake-screenshots.sh — rebuild DEMO-001 and retake the three docs pictures
# Usage: tests/demo/retake-screenshots.sh <work_dir> [<pictures_dir>]
#
#   1. `og dev demo-sample` (bin/og, which scripts/setup.sh --og builds)
#      writes the invented sample DEMO-001 in <work_dir>/genome (fixed seed;
#      no real sample is read), with as many records as the live-site demo;
#   2. steps 36, 27, 24 and 28 run on it unchanged, with the images of
#      versions.env;
#   3. `og dev demo-sample check` stops here if a report section is not ok,
#      so no picture shows a card that says "Not run";
#   4. `og dev pictures take` draws demo-html-report.png,
#      demo-cpic-report.png and demo-multiqc-report.png into <pictures_dir>
#      (default <work_dir>/pictures) in a headless Chrome, and checks the
#      text each one shows: DEMO-001 and no private term;
#   5. `og dev pictures check` runs its self-test, then checks that each
#      picture holds image chunks only;
#   6. tests/demo/build_live_demo.py builds the live-site demo as the docs
#      workflow does, into <work_dir>/live-demo, and `og dev pictures
#      same-sample` stops when its DEMO-001 report and the pictured one
#      differ in Total variants: the docs call them the same sample.
#
# Env: DEMO_RECORDS     variant records of the VCF (default 20000, the
#                       live-site demo's count; step 6 fails on another)
#      OG               the og that writes DEMO-001 and takes the pictures
#                       (default bin/og)
#      CHROME           the Chrome or Chromium to run (default: the first
#                       of google-chrome, chromium and their kin on PATH,
#                       else Google Chrome's macOS app)
#      DEMO_TERMS_FILE  your private terms for the text check, one regular
#                       expression per line; keep it outside the repository
#
# Docker must see <work_dir> and this checkout. With Colima on macOS both
# must sit in a folder Colima mounts (the home folder and /private/tmp by
# default): a folder it does not mount looks empty inside the container. The
# script checks this before it starts.
#
# It writes nothing in docs/: look at the pictures, then copy them into
# docs/images/ (docs/testing.md, "Retaking the report pictures").
set -euo pipefail

WORK=${1:?Usage: $0 <work_dir> [<pictures_dir>]}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
SAMPLE="DEMO-001"
DEMO="${REPO}/tests/demo"
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)
OUT=${2:-${WORK}/pictures}
export GENOME_DIR="${WORK}/genome"
# shellcheck source=../../versions.env
. "${REPO}/versions.env"
ENGINE=${CONTAINER_ENGINE:-docker}

# Only a folder this script made is deleted.
if [ -e "$GENOME_DIR" ] && [ ! -f "${GENOME_DIR}/.demo-genome" ]; then
  echo "ERROR: ${GENOME_DIR} exists and was not written by this script; give another <work_dir>." >&2
  exit 1
fi
rm -rf "$GENOME_DIR"
mkdir -p "$GENOME_DIR"
: > "${GENOME_DIR}/.demo-genome"

# The two images the steps run, pulled first: a registry error is reported as
# one, not as a mount problem below. Up to four tries, 30, 60 then 90 seconds
# apart: a shared runner address can hit Docker Hub's anonymous pull limit.
for img in "$OG_IMAGE" "$MULTIQC_IMAGE"; do
  "$ENGINE" image inspect "$img" >/dev/null 2>&1 && continue
  for try in 1 2 3 4; do
    "$ENGINE" pull -q "$img" && break
    if [ "$try" = 4 ]; then
      echo "ERROR: cannot pull ${img}: check the network and the registry, then run again." >&2
      exit 1
    fi
    sleep $((try * 30))
  done
done

echo "=== Docker sees ${GENOME_DIR} and ${REPO}/bin?"
if ! "$ENGINE" run --rm --network none -v "${GENOME_DIR}:/probe:ro" -v "${REPO}/bin:/pgp-bin:ro" "$OG_IMAGE" \
    test -f /probe/.demo-genome -a -f /pgp-bin/og; then
  echo "ERROR: inside a container ${GENOME_DIR} or ${REPO}/bin is empty: Docker does not mount that folder." >&2
  echo "  Put the checkout and <work_dir> under a folder it mounts (Colima: the home folder or /private/tmp)." >&2
  exit 1
fi

echo "=== 1. DEMO-001"
OG=${OG:-${REPO}/bin/og}
"$OG" dev demo-sample --genome-dir "$GENOME_DIR" --records "${DEMO_RECORDS:-20000}" --repo "$REPO"

echo "=== 2. Steps 36, 27, 24 and 28"
bash "${REPO}/scripts/36-pgx-consensus.sh" "$SAMPLE"
bash "${REPO}/scripts/27-cpic-lookup.sh" "$SAMPLE"
bash "${REPO}/scripts/24-html-report.sh" "$SAMPLE"
bash "${REPO}/scripts/28-multiqc.sh" "$SAMPLE"

echo "=== 3. Every section of the summary is ok?"
"$OG" dev demo-sample check --genome-dir "$GENOME_DIR"

echo "=== 4. Pictures"
TERMS=()
[ -z "${DEMO_TERMS_FILE:-}" ] || TERMS+=(--terms-file "$DEMO_TERMS_FILE")
# The HTML picture also stops above the Traits card: its table of genes and
# genotypes would add about 3,700 pixels below the cards.
"$OG" dev pictures take --genome-dir "$GENOME_DIR" --out "$OUT" ${TERMS[@]+"${TERMS[@]}"} \
  --html-stop-before "Traits,Secondary-Findings Genes,Genes With a Non-Normal,ClinVar Hits,Steps Not Run"

echo "=== 5. Checks"
"$OG" dev pictures check --self-test ${TERMS[@]+"${TERMS[@]}"}
"$OG" dev pictures check "${OUT}"/demo-*.png

echo "=== 6. The same sample as the live-site demo?"
LIVE="${WORK}/live-demo"
rm -rf "$LIVE"
env OG="$OG" python3 "${DEMO}/build_live_demo.py" --out "$LIVE" --og "$OG"
"$OG" dev pictures same-sample "${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html" "${LIVE}/runs/${SAMPLE}/report.html"

echo ""
echo "Pictures in ${OUT} (sizes above). Look at them, then: cp ${OUT}/demo-*.png ${REPO}/docs/images/"
