#!/usr/bin/env bash
# og-leg: OG_REQUIRED
# live_post (scripts/lib/common.sh) sends run-all.sh's step messages to the
# live site (og live serve). This test starts the site on loopback and checks
# that a message reaches events.jsonl, that without OWNGENOME_LIVE_URL nothing is
# sent, and that with the site down it warns once and returns 0 every time, so
# a run never stops for it. It first proves the store check can fail: a
# message to another run must not show up in this run's file.
# It also checks that scripts/owngenome-live.sh hash-password exits non-zero, with no
# success line, when it cannot write the hash file; the control is the old
# one-line `a && b && c` form, under which set -e let the script go on.
#
# LIVE_POST_FROM=FILE takes the function from another file. The site is og
# built for this machine (scripts/setup.sh --og), through bin/og or OG; without
# it the test is skipped with a note, and with OG_REQUIRED set the skip is a
# failure (the og leg of guard.yml sets it).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FROM=${LIVE_POST_FROM:-${REPO}/scripts/lib/common.sh}
WORK=$(mktemp -d)
PID=""
cleanup() { [ -z "$PID" ] || kill "$PID" 2>/dev/null || true; rm -rf "$WORK"; }
trap cleanup EXIT
fail() { echo "FAIL: $*"; exit 1; }

SRC=$(sed -n '/^live_post() {/,/^}/p' "$FROM")
[ -n "$SRC" ] || fail "no live_post() in ${FROM}"
eval "$SRC"
command -v curl >/dev/null || fail "curl is needed"
OG=${OG:-${REPO}/bin/og}
export OG  # scripts/owngenome-live.sh hash-password below runs the same og
if ! "$OG" version > /dev/null 2>&1; then
  [ -z "${OG_REQUIRED:-}" ] || fail "og is not built (scripts/setup.sh --og)"
  echo "skip: og is not built for this machine, and the live site is og live serve (scripts/setup.sh --og)"
  exit 0
fi

# free_port: a port on 127.0.0.1 that nothing answers on (bash's /dev/tcp
# connect is refused), picked at random from 20000 to 52767.
free_port() {
  local p
  for _ in $(seq 100); do
    p=$((20000 + RANDOM % 32768))
    if ! (exec 3<> "/dev/tcp/127.0.0.1/${p}") 2> /dev/null; then echo "$p"; return; fi
  done
  fail "no free port on 127.0.0.1"
}
WEB=$(free_port) ING=$(free_port)
mkdir -p "${WORK}/genome"
"$OG" live serve --root "${WORK}/genome" --state-dir "${WORK}/state" --port "$WEB" \
  --ingest-port "$ING" --quiet 2> "${WORK}/server.err" &
PID=$!
for _ in $(seq 50); do [ -s "${WORK}/state/ingest-token" ] && curl -s -o /dev/null "http://127.0.0.1:${ING}/" && break; sleep 0.1; done
# og live serve takes a message only under the token of its state folder, the first part of the path.
TOK=$(cat "${WORK}/state/ingest-token") || fail "og live serve wrote no ingest token: $(cat "${WORK}/server.err")"
EV="${WORK}/state/runs/demo1/events.jsonl"
MSG='{"event":"pgp_step","step":"14","state":"started"}'

# Without OWNGENOME_LIVE_URL nothing is sent
unset OWNGENOME_LIVE_URL OWNGENOME_LIVE_WARNED
live_post demo1 "$MSG" || fail "live_post without OWNGENOME_LIVE_URL returned non-zero"
[ ! -e "$EV" ] || fail "a message was stored without OWNGENOME_LIVE_URL"

# Control: without the token the site stores nothing, and live_post still returns 0
export OWNGENOME_LIVE_URL="http://127.0.0.1:${ING}/"
live_post demo1 "$MSG" 2> /dev/null || fail "live_post returned non-zero when the site refused it"
[ ! -e "$EV" ] || fail "control: the site stored a message posted without its token"
unset OWNGENOME_LIVE_WARNED

# Control: a message to another run is not in this run's file
export OWNGENOME_LIVE_URL="http://127.0.0.1:${ING}/${TOK}/"
live_post other1 "$MSG"
if [ -e "$EV" ] && grep -q '"step":"14"' "$EV"; then fail "control: another run's message is in demo1's file"; fi
[ -s "${WORK}/state/runs/other1/events.jsonl" ] || fail "control: the other run's message was not stored"
echo "Control: the store check tells one run's file from another's."

# The site up: the message is stored (a URL ending in '/' works too)
live_post demo1 "$MSG"
grep -q '"event":"pgp_step","step":"14","state":"started"' "$EV" || fail "the message is not in ${EV}: $(cat "$EV" 2>/dev/null)"
echo "PASS: live_post stores the message"

# The site down: one warning, exit 0 each time
kill "$PID"; wait "$PID" 2>/dev/null || true; PID=""
unset OWNGENOME_LIVE_WARNED
rc=0; { live_post demo1 "$MSG"; live_post demo1 "$MSG"; } 2> "${WORK}/warn.txt" || rc=$?
[ "$rc" -eq 0 ] || fail "live_post returned ${rc} with the site down"
[ "$(grep -c 'WARNING: the live site' "${WORK}/warn.txt")" -eq 1 ] || fail "expected one warning: $(cat "${WORK}/warn.txt")"
if grep -q -- "$TOK" "${WORK}/warn.txt"; then fail "the warning shows the ingest token: $(cat "${WORK}/warn.txt")"; fi
echo "PASS: with the site down live_post warns once, without the token, and returns 0"

# Control: a live_post that passes the failure on must be caught by the check above
eval "$(printf '%s\n' "$SRC" | sed -e 's/^live_post()/live_post_bad()/' -e 's/^  return 0$/  return 1/')"
rc=0; live_post_bad demo1 "$MSG" 2>/dev/null || rc=$?
[ "$rc" -ne 0 ] || fail "control: the planted live_post_bad returned 0, so the exit check cannot fail"
echo "Control: a live_post that returns the failure is caught."

# hash-password: a hash file it cannot write is an error, never "Password hash written"
hash_into() {  # SCRIPT DIR: run hash-password into DIR/password; prints the exit code, output in WORK/hash.out
  local rc=0
  printf 'synthetic test password 7\n' | GENOME_DIR="${WORK}/genome" OWNGENOME_LIVE_PASSWORD_FILE="$2/password" \
    "$1" hash-password > "${WORK}/hash.out" 2>&1 || rc=$?
  echo "$rc"
}
if [ "$(id -u)" -eq 0 ]; then
  echo "SKIP: hash-password write failure (root writes into a read-only folder)"
else
  mkdir -p "${WORK}/ok" "${WORK}/ro" && chmod 555 "${WORK}/ro"
  rc=$(hash_into "${REPO}/scripts/owngenome-live.sh" "${WORK}/ok")
  [ "$rc" -eq 0 ] && [ "$(stat -c %a "${WORK}/ok/password" 2>/dev/null || stat -f %Lp "${WORK}/ok/password")" = 600 ] \
    || fail "hash-password into a writable folder: exit ${rc}, $(cat "${WORK}/hash.out")"
  rc=$(hash_into "${REPO}/scripts/owngenome-live.sh" "${WORK}/ro")
  { [ "$rc" -ne 0 ] && ! grep -q 'Password hash written' "${WORK}/hash.out"; } \
    || fail "hash-password into a read-only folder: exit ${rc}, $(cat "${WORK}/hash.out")"
  echo "PASS: hash-password fails, with no success line, when it cannot write the hash"
  # Control: the && chain it replaced. A copy beside links to the rest of the repo.
  mkdir -p "${WORK}/tree/scripts"
  for f in "${REPO}"/* "${REPO}"/.[!.]*; do [ "$(basename "$f")" = scripts ] || ln -s "$f" "${WORK}/tree/"; done
  ln -s "${REPO}/scripts/lib" "${WORK}/tree/scripts/lib"
  # The three lines of the fixed write, put back as the one-line chain.
  # shellcheck disable=SC2016  # the dollar signs are the script's
  awk '
    { line[++n] = $0 }
    END {
      for (i = 1; i <= n; i++) {
        if (i + 2 <= n && line[i] == "  (umask 077; printf '\''%s\\n'\'' \"$h\" > \"${PWF}.tmp\")" \
            && line[i + 1] == "  chmod 600 \"${PWF}.tmp\"" && line[i + 2] == "  mv \"${PWF}.tmp\" \"$PWF\"") {
          print "  (umask 077 && printf '\''%s\\n'\'' \"$h\" > \"${PWF}.tmp\") && chmod 600 \"${PWF}.tmp\" && mv \"${PWF}.tmp\" \"$PWF\""
          i += 2
        } else print line[i]
      }
    }' "${REPO}/scripts/owngenome-live.sh" > "${WORK}/tree/scripts/owngenome-live.sh"
  chmod +x "${WORK}/tree/scripts/owngenome-live.sh"
  grep -q 'umask 077 && printf' "${WORK}/tree/scripts/owngenome-live.sh" || fail "control: the planted copy was not made"
  rc=$(hash_into "${WORK}/tree/scripts/owngenome-live.sh" "${WORK}/ro")
  { [ "$rc" -ne 0 ] && ! grep -q 'Password hash written' "${WORK}/hash.out"; } \
    && fail "control: the && chain also failed (exit ${rc}), so the check cannot tell them apart"
  echo "Control: the old && chain exits ${rc} with a success line; the check catches it."
fi
