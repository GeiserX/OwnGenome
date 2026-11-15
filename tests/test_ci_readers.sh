#!/usr/bin/env bash
# The migrated callers run on small public/synthetic inputs without Docker.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
OG=${OG:-${REPO}/bin/og}
t=$(mktemp -d)
trap 'rm -rf "$t"' EXIT
export REPO OG
GENOME_DIR=$t SAMPLE=synthetic FIXTURE_DIR=$t E2E_WORK=$t E2E_NOTES=$t/notes \
  bash "$REPO/tests/e2e/rust-z2e-readers.sh"
bash "$REPO/scripts/ci/image-smoke.sh" --self-test
mkdir -p "$t/bash/vcf"
printf '{"genes":{"CYP2C19":{"sourceDiplotypes":[{"allele1":{"name":"*1"}}]}}}\n' > "$t/bash/vcf/s.report.json"
PARITY_ITEMS=pharmcat bash "$REPO/scripts/ci/parity-diff.sh" --self-test "$t/bash" s
