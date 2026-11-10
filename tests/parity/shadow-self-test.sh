#!/usr/bin/env bash
# shadow-self-test.sh OG: bin/og's shadow mode (OG_SHADOW) with a fake og
# binary and a fake oracle, so each kind of disagreement must fail the call
# and an agreeing call must pass with og's own output. OG is a real og for
# this machine (its dev parity-compare does the comparing). The shim runs
# only on Linux, so elsewhere this says so and stops; CI runs it on Linux
# (scripts/ci/oracle-diff.sh --self-test).
set -euo pipefail
REAL=$(cd "$(dirname "${1:?usage: shadow-self-test.sh OG}")" && pwd)/$(basename "$1")
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
if [ "$(uname -s)" != Linux ]; then
  echo "shadow self-test skipped: bin/og runs only on Linux (CI runs it)"
  exit 0
fi
case $(uname -m) in x86_64 | amd64) arch=x86_64 ;; aarch64 | arm64) arch=aarch64 ;; *) arch=$(uname -m) ;; esac
T=$(mktemp -d "${TMPDIR:-/tmp}/shadow-self-test.XXXXXX")
trap 'rm -rf "$T"' EXIT
mkdir -p "${T}/bin/.og/oracle/fake" "${T}/task"
cp "${ROOT}/bin/og" "${T}/bin/og"
# The fake og: "fake echo IN OUT" copies IN to OUT and says so.
cat > "${T}/bin/.og/og-linux-${arch}" <<EOF
#!/bin/sh
if [ "\$1" = dev ]; then exec "${REAL}" "\$@"; fi
cp "\$3" "\$4"
echo copied
EOF
# The fake oracle, changed as FAKE_ORACLE says; it leaves a mark when it runs.
cat > "${T}/bin/.og/oracle/fake/run" <<EOF
#!/bin/sh
touch "${T}/oracle-ran"
cp "\$1" "\$2"
case "\${FAKE_ORACLE:-}" in
  byte) printf 'a\tc\n' > "\$2" ;;
  exit) echo copied; exit 3 ;;
  stdout) echo "copied, said otherwise"; exit 0 ;;
  extra) echo x > extra.txt ;;
  through) echo x > in.tsv ;;
esac
echo copied
EOF
chmod 755 "${T}/bin/.og/og-linux-${arch}" "${T}/bin/.og/oracle/fake/run"
echo "fake echo" > "${T}/bin/.og/oracle/fake/prefix"

bad=0
# call LABEL WANT_RC GREP [ENV...]: run "og fake echo in.tsv OUT" in a
# fresh task dir with ENV; WANT_RC is 0 or "fail"; GREP must match the output.
# OUT is $CALL_OUT (out.tsv by default); $CALL_PRE names a file the task dir
# holds before the call, besides in.tsv (a path in a subdirectory makes it).
call() {
  local label=$1 want=$2 pattern=$3 rc=0 out
  shift 3
  rm -rf "${T}/task" "${T}/oracle-ran"
  mkdir "${T}/task"
  printf 'a\tb\n' > "${T}/task/in.tsv"
  if [ -n "${CALL_PRE:-}" ]; then
    mkdir -p "$(dirname "${T}/task/${CALL_PRE}")"
    echo old > "${T}/task/${CALL_PRE}"
  fi
  out=$(cd "${T}/task" && env "$@" "${T}/bin/og" fake echo in.tsv "${CALL_OUT:-out.tsv}" 2>&1) || rc=$?
  if { [ "$want" = 0 ] && [ "$rc" -eq 0 ]; } || { [ "$want" = fail ] && [ "$rc" -ne 0 ]; }; then
    if grep -qE "$pattern" <<< "$out"; then echo "ok      ${label}"; return; fi
  fi
  echo "WRONG   ${label} (exit ${rc}):"
  sed 's/^/        /' <<< "$out"
  bad=$((bad + 1))
}
call "without OG_SHADOW og runs alone" 0 '^copied$' OG_SHADOW=
[ ! -e "${T}/oracle-ran" ] || { echo "WRONG   the oracle ran without OG_SHADOW"; bad=$((bad + 1)); }
call "an agreeing oracle passes, with og's output" 0 '^copied$' OG_SHADOW=fake
[ -e "${T}/oracle-ran" ] || { echo "WRONG   OG_SHADOW=fake did not run the oracle"; bad=$((bad + 1)); }
[ -s "${T}/task/out.tsv" ] || { echo "WRONG   og's output is gone after a shadow run"; bad=$((bad + 1)); }
call "a byte the oracle writes differently fails the call" fail 'out.tsv: byte' OG_SHADOW=fake FAKE_ORACLE=byte
call "another exit code fails the call" fail 'exit code: oracle 3, og 0' OG_SHADOW=fake FAKE_ORACLE=exit
call "another stdout fails the call" fail 'stdout:' OG_SHADOW=fake FAKE_ORACLE=stdout
call "a file only the oracle writes fails the call" fail 'files written' OG_SHADOW=fake FAKE_ORACLE=extra
call "OG_SHADOW naming a tool with no oracle fails" fail 'holds no oracle' OG_SHADOW=nosuch
call "an agreeing call says it agrees" 0 'og shadow \(fake\): og and the Python it replaces agree' OG_SHADOW=fake
CALL_PRE=out.tsv call "og overwriting a file that was there fails the call" fail 'og wrote where the shadow cannot compare' OG_SHADOW=fake
CALL_PRE=sub/keep.txt CALL_OUT=sub/out.tsv call "og writing into a directory that was there fails the call" fail 'og wrote where the shadow cannot compare.*sub/out.tsv' OG_SHADOW=fake
call "the oracle writing through a link into the real tree fails the call" fail 'the oracle wrote into the real tree.*in.tsv' OG_SHADOW=fake FAKE_ORACLE=through
rm -f "${T}/oracle-ran"
out=$(cd "${T}/task" && OG_SHADOW=fake FAKE_ORACLE=byte "${T}/bin/og" fake other in.tsv out2.tsv 2>&1) || true
if [ ! -e "${T}/oracle-ran" ] && [ "$out" = copied ]; then echo "ok      a call that is not the tool's is not shadowed"
else echo "WRONG   a call that is not the tool's ran the oracle: ${out}"; bad=$((bad + 1)); fi

if [ "$bad" -gt 0 ]; then echo "shadow self-test: ${bad} checks went wrong"; exit 1; fi
echo "shadow self-test: every check came out as it must"
