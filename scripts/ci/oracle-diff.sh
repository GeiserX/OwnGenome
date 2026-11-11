#!/usr/bin/env bash
# oracle-diff.sh: run each og port and the Python it replaced (the oracle) on
# the same inputs, and require the same results, by meaning. Every port is
# checked with it before its Python is deleted (tests/parity/README.md).
#
#   scripts/ci/oracle-diff.sh [CASE...]
#       runs the cases of tests/parity/cases.tsv (all, or those named
#       TOOL/CASE or TOOL), prints one line per compared item and a total,
#       exits 1 on any difference that is neither known nor golden
#   scripts/ci/oracle-diff.sh materialise [DIR]
#       writes each oracle of tests/parity/oracle.ref into DIR (default
#       bin/.og/oracle): DIR/<tool>/prefix (the og words it replaces),
#       DIR/<tool>/prog and DIR/<tool>/run. bin/og's shadow mode
#       (OG_SHADOW) runs these
#   scripts/ci/oracle-diff.sh --self-test
#       plants each kind of difference in a fake tool and requires each to
#       be reported, then runs og dev parity-compare --self-test and the
#       shadow mode of bin/og against a fake binary
#
# A case runs in two directories that hold the same staged inputs: og runs
# in one, the oracle in the other, with the same arguments. A staged input
# named `stdin` is each side's standard input; without one it is empty.
# Compared items: the exit code, stdout, stderr, the named outputs and every
# new file (stdout and stderr as bytes, files with the case's compare kind,
# through og dev parity-compare). The expected side of an item is the oracle's,
# unless tests/parity/golden/<tool>/<case>/<item> exists (a deliberate fix:
# og must write the golden instead; <item>.absent means it must write
# nothing). A golden the oracle also matches is stale and fails, and so does
# a known.tsv row whose difference is no longer seen.
#
# Env: OG       the og binary (default: rust/target/release/og, built when
#               missing or older than rust/)
#      ORACLE_RUNNER  docker (default; the case's image from versions.env,
#               no network) or local (the host's interpreter, for a quick
#               look on a laptop; CI never uses it)
#      PARITY_WORK    the work directory (default: a temporary one, removed)
set -euo pipefail
export LC_ALL=C

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
PARITY=${PARITY_DIR:-${ROOT}/tests/parity}
CASES="${PARITY}/cases.tsv"
KNOWN="${PARITY}/known.tsv"
REF="${PARITY}/oracle.ref"
GOLDEN="${PARITY}/golden"
RUNNER=${ORACLE_RUNNER:-docker}

die() { echo "oracle-diff: $*" >&2; exit 2; }
rows() { grep -vE '^[[:space:]]*(#|$)' "$1" || true; }

# --- the oracles --------------------------------------------------------------------

# source_at SHA PATH: PATH as it was at SHA; fetches SHA from origin when a
# shallow checkout lacks it.
source_at() {
  git -C "$ROOT" cat-file -e "$1^{commit}" 2> /dev/null \
    || git -C "$ROOT" fetch -q --depth=1 origin "$1" \
    || die "commit $1 is not here and could not be fetched"
  git -C "$ROOT" show "$1:$2"
}

# materialise DIR: one directory per oracle.ref row (awk reads to the end,
# so git show never dies of SIGPIPE under pipefail). Sources:
#   heredoc:FILE:NAME  the body of NAME=$(cat <<'X' ... X) in FILE at SHA
#   file:FILE          FILE at SHA
materialise() {
  local dir=$1 tool sha src interp prefix body
  mkdir -p "$dir"
  while IFS=$'\t' read -r tool sha src interp prefix; do
    [ -n "$prefix" ] || die "oracle.ref row ${tool}: needs tool, sha, source, interpreter and og prefix"
    case "$src" in
      heredoc:*:*)
        local file=${src#heredoc:} name
        name=${file##*:}
        file=${file%:*}
        body=$(source_at "$sha" "$file" | awk -v n="$name" '
          !on && !done && index($0, n "=$(cat <<") == 1 { on = 1; split($0, a, "<<"); t = a[2]; gsub(/[^A-Za-z0-9_]/, "", t); next }
          on && $0 == t { done = 1; on = 0; next }
          on { print }
          END { if (!done) exit 1 }') || die "${tool}: no heredoc ${name} in ${file} at ${sha}"
        ;;
      file:*) body=$(source_at "$sha" "${src#file:}") ;;
      *) die "${tool}: unknown oracle source '${src}'" ;;
    esac
    rm -rf "${dir:?}/${tool}"
    mkdir -p "${dir}/${tool}"
    printf '%s\n' "$body" > "${dir}/${tool}/prog"
    printf '%s\n' "$prefix" > "${dir}/${tool}/prefix"
    printf '#!/bin/sh\nexec %s "$(dirname "$0")/prog" "$@"\n' "$interp" > "${dir}/${tool}/run"
    chmod 755 "${dir}/${tool}/run"
    echo "oracle ${tool}: ${src} at ${sha:0:12} -> ${dir}/${tool}"
  done < <(rows "$REF")
}

# --- one case -----------------------------------------------------------------------

# stage DIR INPUTS: the case's inputs into DIR. INPUTS is comma separated:
# paths under tests/parity/inputs (copied by base name), or gen:SCRIPT, run
# from tests/parity/gen inside DIR (once per run; later cases copy it).
stage() {
  local dir=$1 item
  local -a items
  IFS=, read -r -a items <<< "$2"
  for item in "${items[@]}"; do
    [ "$item" != - ] || continue
    case "$item" in
      gen:*)
        local cache="${WORK}/gen/${item#gen:}"
        if [ ! -d "$cache" ]; then
          mkdir -p "$cache"
          (cd "$cache" && bash "${PARITY}/gen/${item#gen:}") || die "generator ${item#gen:} failed"
        fi
        cp -pR "$cache"/. "$dir"/ ;;
      *)
        [ -f "${PARITY}/inputs/${item}" ] || die "input ${item} is not in tests/parity/inputs"
        cp -p "${PARITY}/inputs/${item}" "$dir"/ ;;
    esac
  done
}

# new_files DIR: every file under DIR not in DIR/.staged, relative, sorted.
new_files() {
  (cd "$1" && find . -type f ! -name .staged | sed 's|^\./||' | sort | comm -23 - .staged)
}

# run_oracle TOOL IMAGE DIR OUT ARGS...
run_oracle() {
  local tool=$1 image=$2 dir=$3 out=$4
  shift 4
  local in=/dev/null
  [ ! -f "${dir}/stdin" ] || in="${dir}/stdin"
  if [ "$RUNNER" = local ]; then
    (cd "$dir" && "${ORACLES}/${tool}/run" "$@") < "$in" > "${out}/stdout" 2> "${out}/stderr" || echo $? > "${out}/exit"
  else
    # Pulled first: docker run's pull progress would land in the oracle's stderr.
    docker image inspect "$image" > /dev/null 2>&1 || docker pull -q "$image" > /dev/null || die "cannot pull ${image}"
    docker run --rm -i --network none --user "$(id -u):$(id -g)" -e HOME=/tmp \
      -v "${ORACLES}:/oracle:ro" -v "${dir}:/w" -w /w "$image" "/oracle/${tool}/run" "$@" \
      < "$in" > "${out}/stdout" 2> "${out}/stderr" || echo $? > "${out}/exit"
  fi
  [ -s "${out}/exit" ] || echo 0 > "${out}/exit"
}

FAILS=0
ITEMS=0
TIMES=""
report() { # STATUS TOOL CASE ITEM DETAIL
  printf '%-7s %s/%s %s%s\n' "$1" "$2" "$3" "$4" "${5:+: $5}"
  ITEMS=$((ITEMS + 1))
  case "$1" in ok | known | golden) ;; *) FAILS=$((FAILS + 1)) ;; esac
}

# differ KIND A B: empty when A and B are the same by KIND, else why. A path
# that does not exist is "absent"; two absent files are the same.
differ() {
  local kind=$1 a=$2 b=$3 out rc=0
  if [ ! -e "$a" ] && [ ! -e "$b" ]; then return 0; fi
  if [ ! -e "$a" ]; then echo "absent on the expected side, written by og"; return 0; fi
  if [ ! -e "$b" ]; then echo "expected, but og did not write it"; return 0; fi
  out=$("$OG" dev parity-compare "$kind" "$a" "$b" 2>&1) || rc=$?
  case "$rc" in
    0) ;;
    1) echo "$out" ;;
    *) echo "cannot compare ($out)" ;;
  esac
}

run_case() {
  local tool=$1 case=$2 inputs=$3 args=$4 outputs=$5 kind=$6 image_var=$7
  local w="${WORK}/run/${tool}/${case}" image prefix item exp why known g
  local -a argv prefixv
  image=${!image_var:-}
  [ -n "$image" ] || die "${tool}/${case}: ${image_var} is not set in versions.env"
  [ -x "${ORACLES}/${tool}/run" ] || die "${tool}: no oracle (is it in oracle.ref?)"
  read -r -a argv <<< "$args"
  read -r -a prefixv < "${ORACLES}/${tool}/prefix"
  rm -rf "$w"
  mkdir -p "${w}/og" "${w}/oracle" "${w}/og.out" "${w}/oracle.out"
  for d in og oracle; do
    stage "${w}/${d}" "$inputs"
    (cd "${w}/${d}" && find . -type f | sed 's|^\./||' | sort > .staged)
  done
  local t0=${EPOCHREALTIME:-0} t1 t2
  local in=/dev/null
  [ ! -f "${w}/og/stdin" ] || in="${w}/og/stdin"
  (cd "${w}/og" && "$OG" "${prefixv[@]}" "${argv[@]}") < "$in" > "${w}/og.out/stdout" 2> "${w}/og.out/stderr" \
    || echo $? > "${w}/og.out/exit"
  [ -s "${w}/og.out/exit" ] || echo 0 > "${w}/og.out/exit"
  t1=${EPOCHREALTIME:-0}
  run_oracle "$tool" "$image" "${w}/oracle" "${w}/oracle.out" "${argv[@]}"
  t2=${EPOCHREALTIME:-0}
  # Wall seconds of each side (bash 5; the oracle's include starting its container).
  TIMES=$(awk -v a="${t0/,/.}" -v b="${t1/,/.}" -v c="${t2/,/.}" 'BEGIN { printf "og %.2f s, oracle %.2f s", b - a, c - b }')

  # The items: exit, stdout, stderr, the outputs the case names and every
  # other file either side wrote.
  local -a items=(exit stdout stderr)
  local -a named
  IFS=, read -r -a named <<< "$outputs"
  for item in "${named[@]}"; do [ "$item" = - ] || items+=("$item"); done
  while IFS= read -r item; do
    case " ${items[*]} " in *" ${item} "*) ;; *) items+=("$item") ;; esac
  done < <({ new_files "${w}/og"; new_files "${w}/oracle"; } | sort -u)

  for item in "${items[@]}"; do
    local ours oracle k=$kind
    case "$item" in
      exit | stdout | stderr) ours="${w}/og.out/${item}"; oracle="${w}/oracle.out/${item}"; k=bytes ;;
      *) ours="${w}/og/${item}"; oracle="${w}/oracle/${item}" ;;
    esac
    g="${GOLDEN}/${tool}/${case}/${item}"
    known=$(awk -F'\t' -v t="$tool" -v c="$case" -v i="$item" '$1 == t && $2 == c && $3 == i {print $4; exit}' <(rows "$KNOWN"))
    if [ -e "$g" ] || [ -e "${g}.absent" ]; then
      exp=$g
      if [ -z "$(differ "$k" "$oracle" "$exp")" ]; then
        report STALE "$tool" "$case" "$item" "the golden is what the oracle writes: delete it"
        continue
      fi
    else
      exp=$oracle
    fi
    why=$(differ "$k" "$exp" "$ours")
    if [ -n "$known" ]; then
      if [ -n "$why" ]; then report known "$tool" "$case" "$item" "$known"
      else report STALE "$tool" "$case" "$item" "known.tsv lists a difference that is gone: remove the row"; fi
    elif [ -n "$why" ]; then
      report DIFF "$tool" "$case" "$item" "$why"
    elif [ "$exp" = "$g" ]; then
      report golden "$tool" "$case" "$item"
    else
      report ok "$tool" "$case" "$item"
    fi
  done
  # A case must write what it names, or its golden must say it writes nothing.
  for item in "${named[@]}"; do
    [ "$item" = - ] || [ -e "${w}/og/${item}" ] || [ -e "${GOLDEN}/${tool}/${case}/${item}.absent" ] \
      || report MISSING "$tool" "$case" "$item" "the case names it, og did not write it"
  done
}

# selected TOOL CASE: 0 when no CASE argument was given, or one names it.
selected() {
  [ "${#PICK[@]}" -eq 0 ] && return 0
  local p
  for p in "${PICK[@]}"; do [ "$p" = "$1" ] || [ "$p" = "$1/$2" ] && return 0; done
  return 1
}

find_og() {
  if [ -n "${OG:-}" ]; then [ -x "$OG" ] || die "OG=${OG} is not executable"; return; fi
  OG="${CARGO_TARGET_DIR:-${ROOT}/rust/target}/release/og"
  if [ ! -x "$OG" ] || [ -n "$(find "${ROOT}/rust" -path "${ROOT}/rust/target" -prune -o -newer "$OG" -type f -print -quit)" ]; then
    echo "Building og (cargo build --release) ..."
    (cd "${ROOT}/rust" && cargo build -q --release --locked -p og) || die "cargo build failed"
  fi
}

run_all() {
  # shellcheck source=../../versions.env
  . "${ROOT}/versions.env"
  find_og
  materialise "${ORACLES}" > /dev/null
  local tool case inputs args outputs kind image n=0
  while IFS=$'\t' read -r tool case inputs args outputs kind image; do
    selected "$tool" "$case" || continue
    [ -n "$image" ] || die "cases.tsv row ${tool}/${case}: needs 7 columns"
    n=$((n + 1))
    run_case "$tool" "$case" "$inputs" "$args" "$outputs" "$kind" "$image"
    echo "        ${tool}/${case}: ${TIMES}"
  done < <(rows "$CASES")
  [ "$n" -gt 0 ] || die "no case matched: ${PICK[*]:-every case}"
  if [ "$FAILS" -gt 0 ]; then
    echo "oracle-diff: ${FAILS} of ${ITEMS} items in ${n} cases differ"
    return 1
  fi
  echo "oracle-diff: all ${ITEMS} items in ${n} cases agree"
}

# --- the self-test ------------------------------------------------------------------

# A fake tool, "fake echo", whose og side is OG_REAL with OG_FAKE_PLANT set,
# and whose oracle is a shell program: each planted kind of difference must
# be reported, and the agreeing ones must pass.
self_test() {
  find_og
  local real=$OG st="${WORK}/self" bad=0 out rc
  mkdir -p "${st}/parity/inputs" "${st}/parity/golden/fake/golden-out" "${st}/parity/golden/fake/golden-stale" "${st}/bin"
  printf 'a\tb\n' > "${st}/parity/inputs/in.tsv"
  # og side: copies in.tsv to out.tsv, says "copied", changed as PLANT says.
  cat > "${st}/bin/og" <<EOF
#!/usr/bin/env bash
if [ "\$1" = dev ]; then exec "${real}" "\$@"; fi
shift 2
cp "\$1" "\$2"
case "\${OG_FAKE_PLANT:-}" in
  byte) printf 'a\tc\n' > "\$2" ;;
  missing) rm "\$2" ;;
  extra) echo x > extra.txt ;;
  exit) exit 3 ;;
  stdout) echo "copied!"; exit 0 ;;
esac
echo "\${OG_FAKE_SAY:-copied}"
EOF
  chmod 755 "${st}/bin/og"
  # The oracle: the same, written as the "old" program.
  mkdir -p "${st}/oracle/fake"
  printf 'fake echo\n' > "${st}/oracle/fake/prefix"
  printf '#!/bin/sh\ncp "$1" "$2"\necho copied\n' > "${st}/oracle/fake/run"
  chmod 755 "${st}/oracle/fake/run"
  # golden-out: the fix makes og print "fixed" (OG_FAKE_SAY) where the
  # oracle printed "copied"; golden-stale: a golden equal to the oracle.
  printf 'fixed\n' > "${st}/parity/golden/fake/golden-out/stdout"
  printf 'copied\n' > "${st}/parity/golden/fake/golden-stale/stdout"

  # expect LABEL WANT PLANT [KNOWN_ROW] [CASE]: run one case of the fake
  # tool; WANT is the status one item must get, or "pass" for all ok.
  expect() {
    local label=$1 want=$2 plant=$3 known=${4:-} case=${5:-plain}
    printf '' > "${st}/parity/known.tsv"
    [ -z "$known" ] || printf '%s\n' "$known" > "${st}/parity/known.tsv"
    printf 'fake\t%s\tin.tsv\tin.tsv out.tsv\tout.tsv\ttsv\tPYTHON_IMAGE\n' "$case" > "${st}/parity/cases.tsv"
    rc=0
    out=$(OG="${st}/bin/og" OG_FAKE_PLANT=$plant OG_FAKE_SAY=$([ "$case" = golden-out ] && echo fixed || echo copied) \
      PARITY_DIR="${st}/parity" ORACLE_RUNNER=local ORACLE_DIR="${st}/oracle" PYTHON_IMAGE=unused \
      bash "$0" --cases-only 2>&1) || rc=$?
    local ok=false
    case "$want" in
      pass) [ "$rc" -eq 0 ] && ! grep -qE '^(DIFF|STALE|MISSING)' <<< "$out" && ok=true ;;
      *) [ "$rc" -eq 1 ] && grep -qE "^${want} " <<< "$out" && ok=true ;;
    esac
    if $ok; then echo "ok      ${label}"
    else echo "WRONG   ${label} (exit ${rc}, wanted ${want})"; sed 's/^/        /' <<< "$out"; bad=$((bad + 1)); fi
  }
  expect "the same output passes" pass ""
  expect "a changed byte in an output file is reported" DIFF byte
  expect "a missing output file is reported" DIFF missing
  expect "an extra output file is reported" DIFF extra
  expect "a wrong exit code is reported" DIFF exit
  expect "a changed stdout is reported" DIFF stdout
  expect "a known difference passes" pass byte $'fake\tplain\tout.tsv\tplanted in the self-test'
  expect "a stale known row fails" STALE "" $'fake\tplain\tout.tsv\tplanted in the self-test'
  expect "a golden og matches passes" pass "" "" golden-out
  expect "a golden the oracle matches is stale" STALE "" "" golden-stale
  local n=10

  echo "--- og dev parity-compare --self-test"
  if "$real" dev parity-compare --self-test; then :; else bad=$((bad + 1)); fi
  n=$((n + 1))
  echo "--- bin/og shadow mode (OG_SHADOW)"
  if bash "${ROOT}/tests/parity/shadow-self-test.sh" "$real"; then :; else bad=$((bad + 1)); fi
  n=$((n + 1))
  if [ "$bad" -gt 0 ]; then echo "oracle-diff self-test: ${bad} of ${n} checks went wrong"; return 1; fi
  echo "oracle-diff self-test: all ${n} checks came out as they must"
}

# --- main ---------------------------------------------------------------------------

WORK=${PARITY_WORK:-}
if [ -z "$WORK" ]; then
  WORK=$(mktemp -d "${TMPDIR:-/tmp}/oracle-diff.XXXXXX")
  trap 'rm -rf "$WORK"' EXIT
fi
ORACLES=${ORACLE_DIR:-${WORK}/oracle}
PICK=()
case "${1:-}" in
  materialise) materialise "${2:-${ROOT}/bin/.og/oracle}" ;;
  --self-test) self_test ;;
  --cases-only)
    # The self-test's inner run: the oracles are already in ORACLE_DIR.
    [ -n "${ORACLE_DIR:-}" ] || die "--cases-only needs ORACLE_DIR"
    materialise() { :; }
    shift
    PICK=("$@")
    run_all ;;
  -h | --help) sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//' ;;
  *) PICK=("$@"); run_all ;;
esac
