#!/usr/bin/env bash
# Every pgsc_calc task gets 'terminate' from the config PRS writes for it.
#
# PRS (modules/local/prs/main.nf) and step 25 (scripts/25-prs.sh) run
# pgsc_calc with images.config, whose generic errorStrategy gives
# 'terminate' to a failure it does not retry. With pgsc_calc's 'finish' a
# memory kill let it wait for its pending tasks, and in its ancestry branch
# that wait could last for ever, so PRS never ended. A process that sets
# errorStrategy itself outranks a generic config line, and at v2.3.0
# MATCH_VARIANTS and PLINK2_SCORE set 'finish'. So both entry points add a
# withName line, which outranks the directive, for the processes that
# modules/local/prs/pgsc_error_override.sh finds in the pgsc_calc copy.
#
#   bash tests/test_pgsc_error_strategy.sh [PGSC_DIR]
#       PGSC_DIR: a pgsc_calc folder. Default: the pinned release
#       (versions.env PGSC_CALC_VERSION), fetched from GitHub and checked
#       against PGSC_CALC_SHA256.
#
# 1. The audit: every errorStrategy in the release's .nf, .config and
#    .groovy files must be one of
#      - a directive of a process that the withName line names, as each
#        entry point's "own errorStrategy" block writes it for that folder;
#      - a generic line of a config's process block: images.config's
#        generic line is read later, so it wins;
#      - a withLabel or withName line that gives 'retry' or 'ignore' (never
#        'finish'; Nextflow ends a retried task with 'terminate' once its
#        retries are spent).
#    Anything else fails, with its file and line: an errorStrategy the
#    audit cannot place may be a new way for pgsc_calc to wait.
# 2. Each entry point puts the line into images.config's process block.
# 3. Controls on a copy of the release: a planted process with its own
#    errorStrategy 'finish' is named by both blocks and passes the audit;
#    the audit with the v2.3.0 list written out by hand reports it (the
#    audit can fail); a planted withName 'finish' in pgsc_calc's conf/ is
#    reported; a module without the block fails.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
MODULE="${REPO}/modules/local/prs/main.nf"
STEP="${REPO}/scripts/25-prs.sh"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0
ok()  { echo "PASS: $*"; }
bad() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }

# --- the pgsc_calc folder -----------------------------------------------------
if [ $# -ge 1 ]; then
  PGSC=$(cd "$1" && pwd)
else
  VER=$(sed -n 's/^PGSC_CALC_VERSION="\(.*\)"$/\1/p' "${REPO}/versions.env")
  SHA=$(sed -n 's/^PGSC_CALC_SHA256="\(.*\)"$/\1/p' "${REPO}/versions.env")
  [ -n "$VER" ] && [ -n "$SHA" ] || { echo "FAIL: no PGSC_CALC_VERSION or PGSC_CALC_SHA256 in versions.env"; exit 1; }
  curl -sS --fail -L --retry 4 --retry-delay 5 --max-time 120 -o "${WORK}/pgsc.tar.gz" \
    "https://github.com/PGScatalog/pgsc_calc/archive/refs/tags/${VER}.tar.gz"
  if command -v sha256sum >/dev/null; then SUM=(sha256sum); else SUM=(shasum -a 256); fi
  echo "${SHA}  ${WORK}/pgsc.tar.gz" | "${SUM[@]}" -c - >/dev/null \
    || { echo "FAIL: pgsc_calc ${VER} does not match PGSC_CALC_SHA256"; exit 1; }
  mkdir "${WORK}/pgsc"
  tar -xzf "${WORK}/pgsc.tar.gz" -C "${WORK}/pgsc" --strip-components 1
  PGSC="${WORK}/pgsc"
  echo "pgsc_calc ${VER}, sha256 checked"
fi

# --- the two entry points' blocks ---------------------------------------------
# block FILE: the lines between "# >>> own errorStrategy" and
# "# <<< own errorStrategy", with Groovy's \$ and \\ undone (a no-op for the
# bash step). Empty when FILE has no block.
block() {
  awk '
    /# >>> own errorStrategy/ { on = 1; next }
    /# <<< own errorStrategy/ { on = 0 }
    on { gsub(/\\\$/, "$"); gsub(/\\\\/, "\\"); print }' "$1"
}
# module_line DIR [MODULE]: the line the module's block writes for DIR. The
# block names the folder by ${pipeline} (its main.nf) and the script by
# ${moduleDir}, which Groovy fills in; they are filled in here as text.
module_line() {
  # shellcheck disable=SC2016  # the Groovy placeholders, as text
  local b m='${moduleDir}' p='${pipeline}'
  b=$(block "${2:-$MODULE}")
  [ -n "$b" ] || { echo "no own errorStrategy block" >&2; return 1; }
  b=${b//"$m"/${REPO}/modules/local/prs}
  b=${b//"$p"/${1}/main.nf}
  bash -c 'set -euo pipefail; '"$b"$'\nprintf "%s\\n" "$OWN_LINE"'
}
# step_line DIR: the line step 25's block writes for DIR.
step_line() {
  local b
  b=$(block "$STEP")
  [ -n "$b" ] || { echo "no own errorStrategy block" >&2; return 1; }
  PGP_ROOT="$REPO" PGSC_CALC_DIR="$1" \
    bash -c 'set -euo pipefail; '"$b"$'\nprintf "%s\\n" "$OWN_LINE"'
}
# names LINE: the process names of a withName line, one per line.
names() {
  sed -n "s/^ *withName: '\([^']*\)' { errorStrategy = 'terminate' }\$/\1/p" <<< "$1" | tr '|' '\n' | sed '/^$/d'
}

# --- the audit ----------------------------------------------------------------
# audit DIR NAMES: every errorStrategy of DIR placed (see the header); NAMES
# is the covered processes, one per line. Prints one line per problem.
audit() {
  local dir=$1 covered=$2
  # shellcheck disable=SC2016  # the dollar signs are awk's
  find -L "$dir" -type f \( -name '*.nf' -o -name '*.config' -o -name '*.groovy' \) | sort | while IFS= read -r f; do
    awk -v covered="$(paste -sd'|' - <<< "$covered")" -v rel="${f#"$dir"/}" '
      BEGIN { n = split(covered, c, "|"); for (i = 1; i <= n; i++) if (c[i] != "") ok[c[i]] = 1
              nf = rel ~ /\.nf$/ }
      {
        line = $0
        code = line; sub(/\/\/.*$/, "", code)
        if (nf && code ~ /^[[:space:]]*process[[:space:]]+[A-Za-z_][A-Za-z0-9_]*/) {
          split(code, w, /[[:space:]]+/); proc = (w[1] == "" ? w[3] : w[2]); sub(/[^A-Za-z0-9_].*$/, "", proc)
        }
        if (code ~ /errorStrategy/) {
          if (nf) {
            if (code !~ /^[[:space:]]*errorStrategy([^A-Za-z0-9_]|$)/)
              print rel ":" FNR ": an errorStrategy that is not a process directive: " line
            else if (proc == "")
              print rel ":" FNR ": an errorStrategy outside a process: " line
            else if (!(proc in ok))
              print rel ":" FNR ": " proc " sets its own errorStrategy and the withName line does not name it"
          } else {
            top = depth ? head[depth] : ""
            if (top == "process") { }
            else if (top ~ /^with(Label|Name)[[:space:]]*:/ && code ~ /errorStrategy[[:space:]]*=[[:space:]]*["\047](retry|ignore)["\047][[:space:]]*$/) { }
            else print rel ":" FNR ": an errorStrategy under \"" top "\" that may be finish: " line
          }
        }
        if (!nf) {
          s = code; opens = gsub(/\{/, "{", s); closes = gsub(/\}/, "}", s)
          if (opens > closes) {
            h = code; sub(/\{.*$/, "", h); gsub(/^[[:space:]]+|[[:space:]]+$/, "", h)
            depth++; head[depth] = h
          } else if (closes > opens && depth > 0) depth--
        }
      }' "$f"
  done
}

LINE=$(module_line "$PGSC") || LINE=''
SLINE=$(step_line "$PGSC") || SLINE=''
echo "module: ${LINE:-<none>}"
echo "step 25: ${SLINE:-<none>}"
COVERED=$(names "$LINE")
if [ -n "$COVERED" ]; then ok "the module's block names $(paste -sd, - <<< "$COVERED")"; else bad "the module's block writes no withName line"; fi
if [ "$SLINE" = "$LINE" ]; then ok "step 25's block writes the same line"; else bad "step 25 writes '${SLINE}', the module '${LINE}'"; fi
PROBLEMS=$(audit "$PGSC" "$COVERED")
if [ -z "$PROBLEMS" ]; then
  ok "the audit places every errorStrategy of the release"
else
  bad "the audit of the release:"; while IFS= read -r l; do echo "  $l"; done <<< "$PROBLEMS"
fi
# A module whose line misses one process must fail the audit too.
if [ -n "$COVERED" ] && [ -n "$(audit "$PGSC" "$(sed 1d <<< "$COVERED")")" ]; then
  ok "the audit reports a process left out of the line"
else
  bad "the audit passes with a process left out of the line"
fi

# --- 2. images.config gets the line --------------------------------------------
# The module prints its process block with one printf; step 25 echoes the line
# inside its { ... } > images.config group.
MOD_WANT="\"\\\$OWN_LINE\" '}' > images.config"
if grep -F "printf '%s\\\\n' 'process {'" "$MODULE" | grep -qF "$MOD_WANT"; then
  ok "the module puts OWN_LINE into images.config's process block"
else
  bad "the module does not put OWN_LINE into images.config's process block"
fi
# shellcheck disable=SC2016  # the dollar signs are the step's, read as text
if awk -v own='[ -z "$OWN_LINE" ] || echo "$OWN_LINE"' '
    $0 == "{" { g = 1; next }
    g && index($0, own) { f = 1 }
    g && $0 == "  echo '\''}'\''" { found = f; exit }
    END { exit !found }' "$STEP"; then
  ok "step 25 puts OWN_LINE into images.config's process block"
else
  bad "step 25 does not put OWN_LINE into images.config's process block"
fi

# --- 3. controls ----------------------------------------------------------------
P="${WORK}/planted"
mkdir -p "$P"
cp -R "${PGSC}/." "$P/"
cat > "${P}/modules/local/planted_task.nf" <<'NF'
process PLANTED_TASK {
    label 'process_low'
    errorStrategy 'finish'

    script:
    """
    true
    """
}
NF
PL=$(module_line "$P") || PL=''
if names "$PL" | grep -qx PLANTED_TASK && [ "$(step_line "$P" || true)" = "$PL" ]; then
  ok "a planted process with its own 'finish' is named by both blocks"
else
  bad "a planted process with its own 'finish' is not named: '${PL}'"
fi
if [ -z "$(audit "$P" "$(names "$PL")")" ]; then ok "the audit passes with the planted process named"; else bad "the audit fails with the planted process named"; fi
HAND=$(audit "$P" "$(printf 'MATCH_VARIANTS\nPLINK2_SCORE\n')")
if grep -q 'PLANTED_TASK sets its own errorStrategy' <<< "$HAND"; then
  ok "the audit with the v2.3.0 list by hand reports the planted process"
else
  bad "the audit with the v2.3.0 list by hand misses the planted process: ${HAND:-nothing}"
fi
rm -f "${P}/modules/local/planted_task.nf"
printf "process {\n    withName: 'PLANTED_TASK' {\n        errorStrategy = 'finish'\n    }\n}\n" > "${P}/conf/planted.config"
if grep -q 'conf/planted.config:3: an errorStrategy under "withName:' <<< "$(audit "$P" "$(names "$PL")")"; then
  ok "the audit reports a planted withName 'finish' in pgsc_calc's conf/"
else
  bad "the audit misses a planted withName 'finish' in pgsc_calc's conf/"
fi
grep -v 'own errorStrategy' "$MODULE" > "${WORK}/noblock.nf"
if module_line "$PGSC" "${WORK}/noblock.nf" >/dev/null 2>&1; then
  bad "a module without the block still gives a line"
else
  ok "a module without the block fails"
fi

if [ "$FAILS" -ne 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "All checks passed"
