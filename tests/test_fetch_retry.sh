#!/usr/bin/env bash
# fetch (scripts/lib/common.sh) is `og setup fetch`: its retries, resumes,
# checks and reuse are tested against a local HTTP server in
# rust/crates/og-setup/tests/fetch.rs. This test keeps what needs the real
# curl and wget: a transfer that stalls fails instead of hanging.
#
# Stalls: each try gives up on a transfer that stops sending data, so fetch
# can try again and resume the .part. A local HTTP server (og dev stub-http)
# sends the start of a file and then nothing, holding the connection open.
# With FETCH_STALL=2, fetch with curl and fetch with wget (curl hidden from
# PATH) each fail within seconds, keep the .part, and try FETCH_TRIES times;
# the same server's whole file downloads under the same bounds (at most 20 s).
# Two controls, a curl and a wget on PATH that drop the bounds og passes them,
# still hang at 25 s, when this test stops them: the behaviour without the
# bounds.

set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
SRV=""
# The killed server's wait returns 143; under set -e that would become the
# test's exit status, so the stop never fails the trap.
trap '[ -z "$SRV" ] || { kill "$SRV"; wait "$SRV"; } 2>/dev/null || true; rm -rf "$WORK"' EXIT
export FETCH_WAIT=0
FAILS=0
LIB="${REPO}/scripts/lib/common.sh"

# --- Stalls -------------------------------------------------------------------
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }
command -v curl >/dev/null || { echo "FAIL: curl is needed to run this test"; exit 1; }
command -v wget >/dev/null || { echo "FAIL: wget is needed to run this test"; exit 1; }

# The server, og dev stub-http: /stall.bin sends a 100000-byte
# Content-Length, 1000 bytes, then nothing for 60 s, holding the connection;
# /whole.bin sends the whole file. Guard's unit tests set OG.
OG=${OG:-${REPO}/bin/og}
export OG
"$OG" version > /dev/null 2>&1 || { echo "FAIL: og is needed to run this test (OG, or bin/og: scripts/setup.sh --og)"; exit 1; }
"$OG" dev stub-http /stall.bin=drip:100000:1000:60000 /whole.bin=drip:100000:100000:0 > "${WORK}/port" & SRV=$!
for _ in $(seq 1 50); do [ -s "${WORK}/port" ] && break; sleep 0.1; done
[ -s "${WORK}/port" ] || { echo "FAIL: the test server did not start"; exit 1; }
URL="http://127.0.0.1:$(cat "${WORK}/port")"

# A PATH without curl, so _get takes its wget branch.
mkdir -p "${WORK}/nocurl"
for c in bash wget awk sed mv rm mkdir cat sleep basename dirname gzip md5sum sha256sum wc head tail tr cut date env \
         grep sort stat readlink uname id timeout nproc getconf printf; do
  p=$(command -v "$c") && ln -sf "$p" "${WORK}/nocurl/$c"
done

# unbounded NAME TOOL: a folder holding TOOL, a wrapper that drops the bounds
# og passes it (curl: --connect-timeout, --speed-limit, --speed-time; wget:
# --connect-timeout, --read-timeout) and runs the real one, ahead of PATH.
unbounded() {
  mkdir -p "${WORK}/$1"
  # shellcheck disable=SC2016  # the wrapper's own code, expanded when it runs
  printf '%s\n' '#!/usr/bin/env bash' \
    'a=(); while [ $# -gt 0 ]; do case "$1" in --connect-timeout|--speed-limit|--speed-time) shift 2 ;; --connect-timeout=*|--read-timeout=*) shift ;; *) a+=("$1"); shift ;; esac; done' \
    "exec $(printf '%q' "$(command -v "$2")") \"\${a[@]}\"" > "${WORK}/$1/$2"
  chmod +x "${WORK}/$1/$2"
}
unbounded curl-unbounded curl
unbounded wget-unbounded wget

# fetch_in LIB PATH NAME FILE [LIMIT]: fetch URL/FILE with LIB, PATH and
# FETCH_STALL=2, two tries, under a LIMIT-second timeout (default 30).
# Prints "rc=<exit> secs=<seconds> part=<bytes>".
fetch_in() {
  local start rc=0 dest="${WORK}/$3/$4"
  start=$(date +%s)
  # shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
  PATH="$2" FETCH_STALL=2 FETCH_TRIES=2 FETCH_WAIT=0 timeout -k 2 "${5:-30}" bash -c '. "$1" && fetch "$2" "$3"' \
    _ "$1" "${URL}/$4" "$dest" > "${WORK}/$3.log" 2>&1 || rc=$?
  echo "rc=${rc} secs=$(( $(date +%s) - start )) part=$( { wc -c < "${dest}.part" || echo 0; } 2>/dev/null | tr -d ' ')"
}
# expect_re DESC GOT WANT [LOG]: GOT matches the extended regex WANT.
expect_re() {
  if [[ "$2" =~ $3 ]]; then pass "$1 ($2)"; else fail "$1: got '$2', want /$3/"; sed 's/^/  /' "${WORK}/${4:-none}.log" 2>/dev/null | tail -n 5; fi
}

expect_re "curl: a stalled transfer fails both tries in seconds and keeps the .part" \
  "$(fetch_in "$LIB" "$PATH" curl-stall stall.bin)" '^rc=1 secs=([0-9]|1[0-9]|20) part=1000$' curl-stall
expect_re "curl: a whole file downloads under the same bounds" \
  "$(fetch_in "$LIB" "$PATH" curl-whole whole.bin)" '^rc=0 ' curl-whole
expect_re "wget: a stalled transfer fails both tries in seconds and keeps the .part" \
  "$(fetch_in "$LIB" "${WORK}/nocurl" wget-stall stall.bin)" '^rc=1 secs=([0-9]|1[0-9]|20) part=1000$' wget-stall
expect_re "wget: a whole file downloads under the same bounds" \
  "$(fetch_in "$LIB" "${WORK}/nocurl" wget-whole whole.bin)" '^rc=0 ' wget-whole
expect_re "control: curl without the bounds hangs until the timeout" \
  "$(fetch_in "$LIB" "${WORK}/curl-unbounded:$PATH" curl-mutant stall.bin 25)" '^rc=124 ' curl-mutant
expect_re "control: wget without the bounds hangs until the timeout" \
  "$(fetch_in "$LIB" "${WORK}/wget-unbounded:${WORK}/nocurl" wget-mutant stall.bin 25)" '^rc=124 ' wget-mutant
for z in 0 00 x; do
  # shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
  if FETCH_STALL=$z bash -c '. "$1" && _get "$2" "$3"' _ "$LIB" "${URL}/whole.bin" "${WORK}/zero" > "${WORK}/zero.log" 2>&1; then
    fail "FETCH_STALL=${z} was accepted"
  else
    expect_re "FETCH_STALL=${z} is refused" "$(cat "${WORK}/zero.log")" 'FETCH_STALL must be'
  fi
done

if [ "$FAILS" -gt 0 ]; then echo "${FAILS} check(s) failed"; exit 1; fi
echo "All fetch checks passed."
