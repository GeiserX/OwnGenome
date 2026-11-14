#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
cp "$ROOT"/tests/fixtures/pharmcat/report-3.2.0.json "$ROOT"/tests/fixtures/pharmcat/report-3.4.0.json "$ROOT"/tests/fixtures/pharmcat/pharmcat-docs-example.json .
