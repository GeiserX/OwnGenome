#!/usr/bin/env bash
# og-leg: OG_REQUIRED
# scripts/owngenome-view.sh, the launcher of the local genome viewer, with
# the og on this machine (OG, else bin/og). The server itself (og view
# serve) is tested in rust/crates/og-live/tests/port_viewer.rs; this checks
# what the script adds, on a 3 kb reference with no reads.
#
# Checks:
#   1. it downloads the igv.js versions.env pins into GENOME_DIR/tools with
#      its sha256, and sets no pin of its own;
#   2. it starts og view serve and prints the address with the locus;
#   3. it picks Manta's diploid SVs before step 22's consensus;
#   4. the server answers on 127.0.0.1, and refuses a foreign Host and a
#      request a proxy forwarded (403);
#   5. it logs each request on stderr and closes each connection;
#   6. a reversed locus stops it with exit 2 before it serves;
#   7. it needs no python3: every run has a python3 that fails first on PATH.
# Needs curl and the network once (igv.js). Without og the test is skipped
# with a note; with OG_REQUIRED set, the skip is a failure.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
OG=${OG:-${REPO}/bin/og}
WORK=$(mktemp -d)
SERVER=""
cleanup() {
  [ -z "$SERVER" ] || kill "$SERVER" 2> /dev/null || true
  rm -rf "$WORK"
}
trap cleanup EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

if ! "$OG" version > /dev/null 2>&1; then
  if [ -n "${OG_REQUIRED:-}" ]; then echo "FAIL: og is not built (scripts/setup.sh --og)"; exit 1; fi
  echo "skip: og is not built for this machine, and the viewer runs it (scripts/setup.sh --og)"
  exit 0
fi
export OG

# A 3 kb chr1, its .fai, and sample S with both SV VCFs and no reads.
G="${WORK}/genome"
mkdir -p "${G}/reference" "${G}/S/manta/results/variants" "${G}/S/sv_merged"
bases=$(head -c 3000 < /dev/zero | tr '\0' 'A')
{ echo ">chr1"; fold -w 60 <<< "$bases"; } > "${G}/reference/synthetic.fasta"
printf 'chr1\t3000\t6\t60\t61\n' > "${G}/reference/synthetic.fasta.fai"
: > "${G}/S/manta/results/variants/diploidSV.vcf.gz"
: > "${G}/S/sv_merged/S_sv_consensus.vcf.gz"
PORT=$((20000 + RANDOM % 20000))

# A python3 that fails, first on PATH: the viewer must not need one.
NOPY="${WORK}/nopy"
mkdir -p "$NOPY"
printf '#!/bin/sh\necho "python3 was called: $*" >&2\nexit 127\n' > "${NOPY}/python3"
chmod 755 "${NOPY}/python3"

view() {
  GENOME_DIR="$G" REF_FASTA="${G}/reference/synthetic.fasta" OWNGENOME_VIEW_PORT="$PORT" READS=none \
    PATH="${NOPY}:${PATH}" bash "${REPO}/scripts/owngenome-view.sh" "$@"
}

# --- 1-2. start it ------------------------------------------------------------------------
view S chr1:451 > "${WORK}/server.log" 2>&1 &
SERVER=$!
up=""
for _ in $(seq 1 240); do
  kill -0 "$SERVER" 2> /dev/null || break
  if grep -qF "Viewer for S: http://127.0.0.1:${PORT}/#chr1:451" "${WORK}/server.log"; then up=1; break; fi
  sleep 0.5
done
if [ -n "$up" ]; then pass "owngenome-view.sh starts og view serve and prints the address"; else
  fail "the viewer did not start: $(cat "${WORK}/server.log")"
  exit 1
fi

# shellcheck source=../versions.env
. "${REPO}/versions.env"
igv="${G}/tools/igv.js-${IGV_JS_VERSION}/igv.min.js"
sum=$( (sha256sum "$igv" 2> /dev/null || shasum -a 256 "$igv") | cut -d' ' -f1)
if [ "$sum" = "$IGV_JS_SHA256" ]; then pass "igv.js ${IGV_JS_VERSION} is downloaded with the pinned sha256"; else
  fail "igv.js has sha256 '${sum}', versions.env pins ${IGV_JS_SHA256}"
fi
if grep -Eq '^[[:space:]]*IGV_JS_(VERSION|SHA256)=' "${REPO}/scripts/owngenome-view.sh"; then
  fail "owngenome-view.sh sets an igv.js pin of its own"
else pass "owngenome-view.sh sets no igv.js pin of its own"; fi

# --- 3. the SV track ----------------------------------------------------------------------
if grep -q 'SV VCF: .*manta/results/variants/diploidSV.vcf.gz' "${WORK}/server.log"; then
  pass "Manta's diploid SVs are picked before step 22's consensus"
else fail "the SV VCF is not Manta's: $(grep 'SV VCF' "${WORK}/server.log")"; fi

# --- 4. what it answers -------------------------------------------------------------------
base="http://127.0.0.1:${PORT}"
code() { curl -s -o /dev/null -w '%{http_code}' "$@"; }
cfg=$(curl -s "${base}/config.json")
if grep -q '"url": "data/2/diploidSV.vcf.gz"' <<< "$cfg"; then pass "the config names Manta's SV VCF"; else
  fail "the config: ${cfg}"
fi
# status DESC WANT CURL-ARGS...: the request answers WANT.
status() {
  local desc=$1 want=$2 got
  shift 2
  got=$(code "$@")
  if [ "$got" = "$want" ]; then pass "${desc} (${got})"; else fail "${desc}: got ${got}, want ${want}"; fi
}
status "the page is served" 200 "${base}/"
status "a foreign Host is refused" 403 -H 'Host: attacker.example' "${base}/config.json"
status "a forwarded request is refused" 403 -H 'X-Forwarded-For: 203.0.113.9' "${base}/config.json"

# --- 4b. it logs each request and closes each connection ----------------------------------
# A handler that blocks on the log keeps its socket open: the client still
# gets the answer (Content-Length), but the server runs out of descriptors.
logged=""
for _ in $(seq 1 20); do
  if grep -qF '"GET /config.json HTTP/1.1" 200 -' "${WORK}/server.log"; then logged=1; break; fi
  sleep 0.5
done
if [ -n "$logged" ]; then pass "each request is logged on stderr"; else
  fail "no request was logged: $(cat "${WORK}/server.log")"
fi
closed=""
if exec 3<> "/dev/tcp/127.0.0.1/${PORT}"; then
  printf 'GET /config.json HTTP/1.1\r\nHost: 127.0.0.1:%s\r\n\r\n' "$PORT" >&3
  # Read to the end of the stream, at most 10 s. A timeout ends the loop too
  # (and with status 1 when a partial line was read, on some bash), so the
  # time it took is what tells the two apart.
  start=$SECONDS
  while IFS= read -r -t 10 _ <&3; do :; done
  [ $((SECONDS - start)) -lt 9 ] && closed=1
  exec 3<&-
fi
if [ -n "$closed" ]; then pass "the server closes the connection after the answer"; else
  fail "the connection was still open 10 s after the answer"
fi
kill "$SERVER" 2> /dev/null || true
wait "$SERVER" 2> /dev/null || true
SERVER=""

# --- 5. a reversed locus ------------------------------------------------------------------
rc=0
out=$(view S chr1:2000-1000 2>&1) || rc=$?
if [ "$rc" -eq 2 ] && grep -q "ends before it starts" <<< "$out"; then pass "a reversed locus stops it (exit 2)"; else
  fail "a reversed locus: exit ${rc}: ${out}"
fi

echo
if [ "$FAILS" -eq 0 ]; then echo "RESULT: ALL PASS"; else echo "RESULT: ${FAILS} FAILED"; exit 1; fi
