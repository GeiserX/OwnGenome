#!/usr/bin/env bash
# Run the cases of tests/parity/cases.tsv against their frozen goldens.
#
#   scripts/ci/oracle-diff.sh [TOOL[/CASE] ...]
#   scripts/ci/oracle-diff.sh --self-test
#   scripts/ci/oracle-diff.sh --check-goldens
#
# OG selects an existing binary (otherwise builds a release og). PARITY_WORK
# keeps the staged inputs and actual outputs in the given absolute directory.
set -euo pipefail
export LC_ALL=C TZ=UTC

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
PARITY=${PARITY_DIR:-"$ROOT/tests/parity"}
CASES="$PARITY/cases.tsv"
GOLDEN="$PARITY/golden"
CAP=1000000

die() { echo "oracle-diff: $*" >&2; exit 2; }
rows() { awk '!/^[[:space:]]*(#|$)/' "$1"; }

# Count stored bytes across all files in each tool/case directory, including
# absence markers and checksums.
check_goldens() {
  find "$GOLDEN" -type f -exec wc -c {} + | awk -v root="$GOLDEN/" -v cap="$CAP" '
    index($0, root) {
      path = substr($0, index($0, root) + length(root))
      if (split(path, parts, "/") < 3) { print "invalid golden path: " path; bad = 1; next }
      sizes[parts[1] "/" parts[2]] += $1
      total += $1
      files++
    }
    END {
      for (id in sizes) {
        if (sizes[id] > cap) { printf "%s: %d bytes exceeds %d\n", id, sizes[id], cap; bad = 1 }
        if (sizes[id] > max) max = sizes[id]
        cases++
      }
      printf "Goldens: %d bytes in %d files across %d cases (largest case %d bytes, cap %d)\n", total, files, cases, max, cap
      exit (bad || !files)
    }'
}

stage() {
  local dir=$1 item gen cache
  local -a items
  IFS=, read -r -a items <<< "$2"
  for item in "${items[@]}"; do
    case "$item" in
      -) ;;
      gen:*)
        gen=${item#gen:}
        cache="$WORK/gen/$gen"
        if [ ! -d "$cache" ]; then
          mkdir -p "$cache"
          (cd "$cache" && bash "$PARITY/gen/$gen") || die "generator $gen failed"
        fi
        cp -pR "$cache"/. "$dir"/ ;;
      fixture:*)
        tar -xzf "$GOLDEN/fixtures/${item#fixture:}/inputs.tar.gz" -C "$dir" ;;
      *)
        [ -f "$PARITY/inputs/$item" ] || die "missing input $item"
        cp -p "$PARITY/inputs/$item" "$dir"/ ;;
    esac
  done
}

new_files() {
  (cd "$1" && find . -type f -print | sed 's|^\./||' | sort | comm -23 - "$2")
}

digest() {
  "$OG" dev parity-compare --canonical "$1" "$2" | sha256sum | awk 'NR == 1 {print $1}'
}

difference() {
  local kind=$1 expected=$2 actual=$3 want got status=0
  if [ -f "$expected.absent" ]; then
    if [ -e "$expected" ] || [ -e "$expected.sha256" ]; then
      echo "conflicting goldens"
    elif [ -e "$actual" ]; then
      echo "expected absence, but og wrote a file"
    fi
  elif [ -f "$expected.sha256" ]; then
    if [ -e "$expected" ]; then echo "conflicting goldens"; return; fi
    if [ ! -f "$actual" ]; then echo "expected a checksummed file, but og did not write it"; return; fi
    want=$(cat "$expected.sha256")
    if ! [[ "$want" =~ ^[a-f0-9]{64}$ ]]; then echo "invalid golden checksum"; return; fi
    got=$(digest "$kind" "$actual") || { echo "cannot checksum output ($kind)"; return; }
    [ "$got" = "$want" ] || echo "checksum differs (expected $want, got $got)"
  elif [ -f "$expected" ]; then
    if [ ! -f "$actual" ]; then echo "expected a file, but og did not write it"; return; fi
    got=$("$OG" dev parity-compare "$kind" "$expected" "$actual" 2>&1) || status=$?
    case "$status" in
      0) ;;
      1) echo "$got" ;;
      *) echo "cannot compare ($got)" ;;
    esac
  else
    echo "no golden for this item"
  fi
}

# The two Yleaf data-path messages report the staged working directory.
# Preserve the rest of each line, including the name under that directory.
normalize_path() {
  CASE_PATH=$1 awk '{
    s = $0; from = ENVIRON["CASE_PATH"]; result = ""
    while ((i = index(s, from))) {
      result = result substr(s, 1, i - 1) "<case-dir>"
      s = substr(s, i + length(from))
    }
    print result s
  }' "$2" > "$2.normalized"
  mv "$2.normalized" "$2"
}

FAILS=0
ITEMS=0
report() {
  printf '%-7s %s/%s %s%s\n' "$1" "$2" "$3" "$4" "${5:+: $5}"
  ITEMS=$((ITEMS + 1))
  if [ "$1" != ok ]; then FAILS=$((FAILS + 1)); fi
}

# kind_for ITEM SPEC: the compare kind of ITEM. SPEC is one kind, or a
# space-separated list of GLOB=KIND (the first glob ITEM matches, as in
# [[ ITEM == GLOB ]]) ending in the kind of every other item:
# "*.zip=zip-content bytes". A list starts with GLOB=KIND, whose text before
# the first = has no colon, so a kind such as mask:REGEX is never split.
kind_for() {
  local item=$1 spec=$2 el glob first=${2%% *}
  local -a els
  case "${first%%=*}" in *:*) printf '%s\n' "$spec"; return ;; esac
  [[ "$first" == *=* ]] || { printf '%s\n' "$spec"; return; }
  read -r -a els <<< "$spec"   # read, not word splitting: no file name expansion
  for el in "${els[@]}"; do
    if [[ "$el" == *=* && "${el%%=*}" != *:* ]]; then
      glob=${el%%=*}
      # shellcheck disable=SC2053  # the glob is meant to match
      if [[ "$item" == $glob ]]; then printf '%s\n' "${el#*=}"; return; fi
    else
      printf '%s\n' "$el"
      return
    fi
  done
  die "compare kinds '${spec}' end without the kind of every other item (${item})"
}

run_case() {
  local tool=$1 case=$2 inputs=$3 args=$4 outputs=$5 kind=$6 prefix=$7
  local work="$WORK/run/$tool/$case" expected="$GOLDEN/$tool/$case"
  local file actual item k override why input=/dev/null rc=0
  local -a argv prefixv named overrides
  if [ ! -d "$expected" ]; then
    report DIFF "$tool" "$case" goldens "missing case directory"
    return
  fi
  read -r -a argv <<< "$args"
  read -r -a prefixv <<< "$prefix"
  IFS=, read -r -a named <<< "$outputs"
  IFS=';' read -r -a overrides <<< "$kind"
  rm -rf "$work"
  mkdir -p "$work/files" "$work/result"
  stage "$work/files" "$inputs"
  (cd "$work/files" && find . -type f -print | sed 's|^\./||' | sort) > "$work/staged"
  [ ! -f "$work/files/stdin" ] || input="$work/files/stdin"
  (cd "$work/files" && "$OG" "${prefixv[@]}" "${argv[@]}") < "$input" \
    > "$work/result/stdout" 2> "$work/result/stderr" || rc=$?
  printf '%s\n' "$rc" > "$work/result/exit"
  case "$tool/$case" in
    y-yleaf-positions/data) normalize_path "$work/files" "$work/result/stdout" ;;
    y-yleaf-positions/no-tree) normalize_path "$work/files" "$work/result/stderr" ;;
  esac

  {
    printf '%s\n' exit stdout stderr "${named[@]}"
    new_files "$work/files" "$work/staged"
    (cd "$expected" && find . -type f -print | sed -e 's|^\./||' -e 's/\.absent$//' -e 's/\.sha256$//')
  } | sort -u > "$work/items"
  while IFS= read -r item; do
    case "$item" in
      -) continue ;;
      exit | stdout | stderr) actual="$work/result/$item"; k=bytes ;;
      *)
        actual="$work/files/$item"
        k=$(kind_for "$item" "${overrides[0]}")
        for override in "${overrides[@]:1}"; do
          # shellcheck disable=SC2053  # the left side is a glob on purpose
          if [[ "$item" == ${override%%=*} ]]; then k=${override#*=}; fi
        done ;;
    esac
    file="$expected/$item"
    why=$(difference "$k" "$file" "$actual")
    if [ -n "$why" ]; then report DIFF "$tool" "$case" "$item" "$why"
    else report ok "$tool" "$case" "$item"; fi
  done < "$work/items"
  if [ -z "${PARITY_WORK:-}" ]; then rm -rf "$work"; fi
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
  find_og
  export OG OG_BIN_TABLES="$ROOT/bin"
  local tool case inputs args outputs kind prefix extra n=0
  while IFS=$'\t' read -r tool case inputs args outputs kind prefix extra; do
    selected "$tool" "$case" || continue
    [ -n "$prefix" ] && [ -z "$extra" ] || die "cases.tsv row $tool/$case: expected 7 columns"
    n=$((n + 1))
    run_case "$tool" "$case" "$inputs" "$args" "$outputs" "$kind" "$prefix"
  done < <(rows "$CASES")
  [ "$n" -gt 0 ] || die "no case matched: ${PICK[*]:-every case}"
  if [ "$FAILS" -gt 0 ]; then
    echo "oracle-diff: $FAILS of $ITEMS items in $n cases differ"
    return 1
  fi
  echo "oracle-diff: all $ITEMS items in $n cases agree"
}

self_test() {
  find_og
  local real=$OG st="$WORK/self" found=0 n=0 log rc
  local case_dir="$st/parity/golden/fake/copy"
  mkdir -p "$st/parity/inputs" "$case_dir" "$st/bin"
  printf 'a\tb\n' > "$st/parity/inputs/in.tsv"
  cp -p "$st/parity/inputs/in.tsv" "$case_dir/out.tsv"
  printf '0\n' > "$case_dir/exit"
  printf 'copied\n' > "$case_dir/stdout"
  : > "$case_dir/stderr"
  : > "$case_dir/never.absent"
  cat > "$st/bin/og" <<'FAKE'
#!/usr/bin/env bash
set -eu
if [ "$1" = dev ]; then exec "$OG_REAL" "$@"; fi
cp in.tsv out.tsv
case "${OG_FAKE_PLANT:-}" in
  byte) printf 'a\tc\n' > out.tsv ;;
  missing) rm out.tsv ;;
  extra) echo x > extra.txt ;;
  absent) echo x > never ;;
  exit) exit 3 ;;
  stdout) echo changed; exit 0 ;;
  stderr) echo changed >&2 ;;
esac
echo copied
FAKE
  chmod 755 "$st/bin/og"
  printf 'fake\tcopy\tin.tsv\t \tout.tsv,never\t%s\tfake echo\n' '*.zip=zip-content out.*=tsv bytes' > "$st/parity/cases.tsv"

  expect() {
    local label=$1 plant=$2 reject=${3:-}
    rc=0
    log=$(OG="$st/bin/og" OG_REAL="$real" OG_FAKE_PLANT="$plant" PARITY_DIR="$st/parity" \
      PARITY_WORK='' bash "$0" 2>&1) || rc=$?
    n=$((n + 1))
    if { [ -z "$reject" ] && [ "$rc" -eq 0 ]; } ||
      { [ -n "$reject" ] && [ "$rc" -eq 1 ] && grep -qF "fake/copy $reject:" <<< "$log"; }; then
      echo "ok      $label"
    else
      printf 'WRONG   %s (exit %s)\n%s\n' "$label" "$rc" "$log"
      found=$((found + 1))
    fi
  }

  expect "matching outputs and an absent named file pass" ""
  expect "a changed field, using its own comparison kind, fails" byte out.tsv
  expect "a missing file fails" missing out.tsv
  expect "an unexpected file fails" extra extra.txt
  expect "a file with an absence marker fails if written" absent never
  expect "a changed exit code fails" exit exit
  expect "changed stdout fails" stdout stdout
  expect "changed stderr fails" stderr stderr
  printf 'fake\tcopy\tin.tsv\t \tout.tsv,never\tbytes\tfake echo\n' > "$st/parity/cases.tsv"
  digest bytes "$case_dir/out.tsv" > "$case_dir/out.tsv.sha256"
  rm -f "$case_dir/out.tsv"
  expect "a matching checksum passes" ""
  expect "a changed file fails its checksum" byte out.tsv
  printf 'fake\tcopy\tin.tsv\t \tout.tsv,never\t%s\tfake echo\n' 'mask:^masked:' > "$st/parity/cases.tsv"
  printf 'masked: expected\na\tb\n' > "$st/parity/inputs/in.tsv"
  digest 'mask:^masked:' "$st/parity/inputs/in.tsv" > "$case_dir/out.tsv.sha256"
  printf 'masked: actual\na\tb\n' > "$st/parity/inputs/in.tsv"
  expect "a changed masked line keeps the same checksum" ""
  printf 'masked: actual\na\tc\n' > "$st/parity/inputs/in.tsv"
  expect "a changed unmasked line fails its checksum" "" out.tsv
  dd if=/dev/zero of="$case_dir/cap-probe" bs=1000001 count=1 2>/dev/null
  n=$((n + 1))
  rc=0
  log=$(PARITY_DIR="$st/parity" bash "$0" --check-goldens 2>&1) || rc=$?
  if [ "$rc" -eq 1 ] && grep -qF "bytes exceeds 1000000" <<< "$log"; then
    echo "ok      an oversized case fails the cap check"
  else
    printf 'WRONG   oversized case (exit %s)\n%s\n' "$rc" "$log"
    found=$((found + 1))
  fi
  if "$real" dev parity-compare --self-test; then :; else found=$((found + 1)); fi
  n=$((n + 1))
  if [ "$found" -gt 0 ]; then
    echo "oracle-diff self-test: $found of $n checks went wrong"
    return 1
  fi
  echo "oracle-diff self-test: all $n checks came out as they must"
}

WORK=${PARITY_WORK:-}
if [ -z "$WORK" ]; then
  WORK=$(mktemp -d "${TMPDIR:-/tmp}/oracle-diff.XXXXXX")
  trap 'rm -rf "$WORK"' EXIT
fi
PICK=()
case "${1:-}" in
  --self-test) self_test ;;
  --check-goldens) check_goldens ;;
  -h | --help) sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//' ;;
  *) PICK=("$@"); run_all ;;
esac
