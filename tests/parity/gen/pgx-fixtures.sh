#!/usr/bin/env bash
# The PharmCAT reports of tests/fixtures/pharmcat and the two tables
# `og pgx cpic-report` reads from assets/, copied as they are, so the pgx-parse
# cases run on the files the pipeline and its tests use.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
cp "${root}"/tests/fixtures/pharmcat/report-3.4.0.json "${root}"/tests/fixtures/pharmcat/report-3.2.0.json \
  "${root}"/tests/fixtures/pharmcat/pharmcat-docs-example.json \
  "${root}"/assets/pgx_inhibitors.tsv "${root}"/assets/pgx_contraindications.tsv .
