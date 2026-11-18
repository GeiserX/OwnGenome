#!/usr/bin/env bash
# The family page keeps a child's closed de novo and compound-het rows off the
# page, the rows of a grandchild read against that child's genome included:
# tests/live/family_hide_closed.js runs assets/live/core.js and
# family-page.js under node on an invented family, with its own controls (the
# trio rule cut out of core.js must fail). No node is a failure, never a skip.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
if ! command -v node >/dev/null 2>&1; then
  echo "FAIL: node is not on PATH; this test needs node 20 or newer"
  exit 1
fi
exec node "${REPO}/tests/live/family_hide_closed.js" "$REPO"
