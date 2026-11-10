#!/usr/bin/env bash
# The PharmCAT reports of tests/fixtures/pharmcat and the two tables
# `og pgx cpic-report` reads from assets/, copied as they are, so the pgx-parse
# cases run on the files the pipeline and its tests use.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
cp "${root}"/tests/fixtures/pharmcat/report-3.4.0.json "${root}"/tests/fixtures/pharmcat/report-3.2.0.json \
  "${root}"/tests/fixtures/pharmcat/pharmcat-docs-example.json \
  "${root}"/assets/pgx_inhibitors.tsv "${root}"/assets/pgx_contraindications.tsv .
# PharmCAT's example with VKORC1 called -1639 GA (the warfarin rule's A
# allele), made here rather than committed: the example's own guideline text
# names genotypes, which the personal-data scan refuses outside tests/fixtures.
python3 - <<'PY'
import json
ex = json.load(open("pharmcat-docs-example.json"))
ex["genes"]["VKORC1"]["sourceDiplotypes"][0].update(
    label="rs9923231 reference (C)/rs9923231 variant (T)", phenotypes=["-1639 GA"])
with open("example-vkorc1-ga.json", "w") as f:
    f.write(json.dumps(ex, indent=1) + "\n")
PY
