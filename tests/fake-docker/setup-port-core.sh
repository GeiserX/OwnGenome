#!/usr/bin/env bash
# og setup install (rust/crates/og-setup/src/install): setup.sh hands each
# flag og has ported to og and keeps every other flag; OG_SETUP=bash keeps
# them all in setup.sh. For every ported flag, og and setup.sh's bash do the
# same thing against the fake downloads and the fake docker: the same exit
# code, the same messages, the same files with the same modes, and the same
# docker and download calls in the same order.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# Each run below says which side it is on (OG_SETUP), whatever the suite was given.
unset OG_SETUP

# An og that notes each `setup install` it is handed, then runs the real one.
mkdir -p "${CASE_WORK}/og-bin"
cat > "${CASE_WORK}/og-bin/og" <<'SHIM'
#!/usr/bin/env bash
if [ "${1:-} ${2:-}" = "setup install" ]; then printf '%s\n' "$*" >> "${CASE_WORK}/og.calls"; fi
exec "${REPO_ROOT}/bin/og" "$@"
SHIM
chmod +x "${CASE_WORK}/og-bin/og"
export OG="${CASE_WORK}/og-bin/og"
: > "${CASE_WORK}/og.calls"
G="${CASE_WORK}/genome"
mkdir -p "$G"

# --- who runs a flag --------------------------------------------------------------
"$OG" setup ported --kir-data > /dev/null || fail "og setup ported: --kir-data is ported"
rc=0; "$OG" setup ported --archaic > /dev/null || rc=$?
[ "$rc" -eq 1 ] || fail "og setup ported --archaic exited ${rc}, expected 1 (bash runs it)"
rc=0; "$OG" setup ported --nope 2> /dev/null || rc=$?
[ "$rc" -eq 2 ] || fail "og setup ported --nope exited ${rc}, expected 2 (no such flag)"

run_expect 2 unknown "${SCRIPTS}/setup.sh" --nope "$G"
[ ! -s "${CASE_WORK}/og.calls" ] || fail "setup.sh handed an unknown flag to og: $(cat "${CASE_WORK}/og.calls")"
output_has unknown '^ERROR: unknown setup flag: --nope$'
run_expect 2 unknown-og "$OG" setup install --nope "$G"
output_has unknown-og '^ERROR: unknown setup flag: --nope$'
run_expect 2 two-flags "$OG" setup install --kir-data "$G" --viral-data
output_has two-flags 'one flag at a time: --viral-data comes after --kir-data'
n=$(wc -l < "${CASE_WORK}/og.calls")
OG_SETUP=rust run_expect 2 bad-og-setup "${SCRIPTS}/setup.sh" --kir-data "$G"
[ "$(wc -l < "${CASE_WORK}/og.calls")" -eq "$n" ] || fail "setup.sh handed --kir-data to og under OG_SETUP=rust"
output_has bad-og-setup "OG_SETUP must be og or bash, got 'rust'"
OG_SETUP=rust run_expect 2 bad-og-setup-og "$OG" setup install --kir-data "$G"
output_has bad-og-setup-og "OG_SETUP must be og or bash, got 'rust'"
run_expect 0 install-help "$OG" setup install --help
output_has install-help 'og setup install'

# --- og against setup.sh's bash, flag by flag ---------------------------------------
# SIDE_TREE DIR: every file and link under DIR with its mode and checksum.
side_tree() {
  (cd "$1" && find . \( -type f -o -type l \) -print | LC_ALL=C sort | while IFS= read -r f; do
    if [ -L "$f" ]; then printf '%s -> %s\n' "$f" "$(readlink "$f")"
    else printf '%s %s %s\n' "$f" "$(stat -c %a "$f" 2>/dev/null || stat -f %Lp "$f")" "$(cksum < "$f")"; fi
  done)
}

# compare NAME ARGS...: setup.sh ARGS under OG_SETUP=bash and under og, each
# in its own install folder (named G in what is compared).
compare() {
  local name=$1 side g start calls want_calls
  shift
  for side in bash og; do
    g="${CASE_WORK}/${name}-${side}/genome"
    mkdir -p "$g"
    start=$(wc -l < "$FAKE_DOCKER_LOG")
    calls=$(wc -l < "${CASE_WORK}/og.calls")
    OG_SETUP=$side run_rc "${name}-${side}" "${SCRIPTS}/setup.sh" "${@//@G@/$g}"
    echo "$RC" > "${CASE_WORK}/${name}-${side}.rc"
    sed "s|${g}|G|g" "${CASE_WORK}/${name}-${side}.out" > "${CASE_WORK}/${name}-${side}.said"
    tail -n "+$((start + 1))" "$FAKE_DOCKER_LOG" | sed "s|${g}|G|g" > "${CASE_WORK}/${name}-${side}.calls"
    if [ -d "$g" ]; then side_tree "$g" > "${CASE_WORK}/${name}-${side}.tree"; else : > "${CASE_WORK}/${name}-${side}.tree"; fi
    want_calls=$calls
    [ "$side" = og ] && want_calls=$((calls + 1))
    [ "$(wc -l < "${CASE_WORK}/og.calls")" -eq "$want_calls" ] \
      || fail "${name}: OG_SETUP=${side} setup.sh $* was$([ "$side" = og ] || echo ' not') to be handed to og; og.calls: $(cat "${CASE_WORK}/og.calls")"
  done
  for what in rc said calls tree; do
    diff -u "${CASE_WORK}/${name}-bash.${what}" "${CASE_WORK}/${name}-og.${what}" \
      || fail "${name}: og and setup.sh's bash differ in ${what} (the diff above, bash first)"
  done
  echo "[OK] ${name}: og and setup.sh's bash give exit $(cat "${CASE_WORK}/${name}-og.rc"), the same messages, $(wc -l < "${CASE_WORK}/${name}-og.calls" | tr -d ' ') calls and $(wc -l < "${CASE_WORK}/${name}-og.tree" | tr -d ' ') files"
}

ported=()
while IFS= read -r flag; do
  if "$OG" setup ported "$flag" > /dev/null 2>&1; then ported+=("$flag"); fi
done < <("$OG" setup flags)
[ "${#ported[@]}" -ge 8 ] || fail "og runs ${#ported[@]} flags, expected at least the eight small installers: ${ported[*]}"

for flag in "${ported[@]}"; do
  case $flag in
    --pull-only)
      compare pull-only --pull-only
      FAKE_DOCKER_MISSING_IMAGES=deepvariant compare pull-missing --pull-only ;;
    *)
      compare "${flag#--}" "$flag" @G@
      compare "${flag#--}-no-keep" --no-keep-downloads "$flag" @G@ ;;
  esac
done
[ "$(grep -c '^setup install' "${CASE_WORK}/og.calls")" -ge $(( 2 * ${#ported[@]} )) ] \
  || fail "og ran fewer installs than compared: $(cat "${CASE_WORK}/og.calls")"
