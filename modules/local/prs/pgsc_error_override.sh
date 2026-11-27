#!/usr/bin/env bash
# pgsc_error_override.sh: the config line that gives 'terminate' to every
# pgsc_calc process that sets its own errorStrategy.
#
# Usage: bash pgsc_error_override.sh <pgsc_calc folder>
#
# PRS (modules/local/prs/main.nf) and step 25 (scripts/25-prs.sh) run
# pgsc_calc with a config of their own, images.config, whose generic
# errorStrategy gives 'terminate' to a failure it does not retry. A process
# that sets errorStrategy itself, a directive in its own .nf file, outranks
# that generic line: at v2.3.0 MATCH_VARIANTS and PLINK2_SCORE set 'finish'.
# Only a withName line outranks a directive. This script finds those
# processes in the pgsc_calc copy the run uses, so a pgsc_calc release that
# adds one is covered with no edit here, and prints
#
#     withName: 'MATCH_VARIANTS|PLINK2_SCORE' { errorStrategy = 'terminate' }
#
# for images.config's process block, or nothing when no process sets it.
# A process is named when a line of its body starts with errorStrategy; a
# process starts at a line that starts with "process NAME". The process of
# each errorStrategy in pgsc_calc's .nf and .config files is checked by
# tests/test_pgsc_error_strategy.sh, on the pinned release.
set -euo pipefail

dir=${1:?usage: pgsc_error_override.sh <pgsc_calc folder>}
[ -f "${dir}/main.nf" ] || { echo "ERROR: ${dir} has no main.nf: not a pgsc_calc folder" >&2; exit 1; }

# shellcheck disable=SC2016  # the dollar signs are awk's
names=$(find -L "$dir" -type f -name '*.nf' -exec awk '
  FNR == 1 { proc = "" }
  /^[[:space:]]*process[[:space:]]+[A-Za-z_][A-Za-z0-9_]*/ {
    proc = $2; sub(/[^A-Za-z0-9_].*$/, "", proc)
  }
  proc != "" && /^[[:space:]]*errorStrategy([^A-Za-z0-9_]|$)/ { print proc }
' {} + | sort -u | paste -sd'|' -)
if [ -n "$names" ]; then
  printf "    withName: '%s' { errorStrategy = 'terminate' }\n" "$names"
fi
