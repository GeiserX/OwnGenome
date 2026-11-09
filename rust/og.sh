#!/usr/bin/env bash
# og.sh: how the og binary reaches bin/.og, and the check that it is the one
# this checkout's rust/ makes. bin/og (the shim tasks call) runs
# bin/.og/og-linux-<arch>; both arches are kept, because an Apple Silicon
# host runs amd64 images under emulation, where uname -m says x86_64.
#
#   rust/og.sh install     scripts/setup.sh --og: on a release tag (vX.Y.Z,
#                          clean rust/) download the release's binaries and
#                          check them against SHA256SUMS (and the build
#                          attestation when gh is installed); otherwise, or
#                          when the release has none, build them
#   rust/og.sh build       build both from rust/ in OG_BUILDER (versions.env)
#   rust/og.sh check [--required]
#                          one line, OK, NOTE or ERROR; exit 1 on ERROR: a
#                          binary built from another rust/, or (with
#                          --required, or once REQUIRED below is true) none
#   rust/og.sh tree-hash   the git tree hash of rust/ as it is on disk
#   rust/og.sh built       the tree hash the binaries were built from, for
#                          run_manifest.tsv ("not built" when there are none)
#
# The tree hash is git's: `git rev-parse HEAD:rust` on a clean checkout, and
# for an edited rust/ the tree those edits would commit (a throwaway index),
# so any edit, committed or not, makes the built binary stale. The build
# stamps it into the binary as og-tree-hash:<hash>; so the check
# reads it with grep and never runs a binary the host may not be able to.
set -euo pipefail

# True once a pipeline step calls og: until then a checkout without the
# binary runs as before, and only a stale one is refused.
REQUIRED=false

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BIN="${ROOT}/bin/.og"
ARCHES="x86_64 aarch64"
REPO="GeiserX/OwnGenome"
ENGINE=${CONTAINER_ENGINE:-docker}

triple() { echo "${1}-unknown-linux-musl"; }

tree_hash() {
  local idx rc=0
  git -C "$ROOT" rev-parse --verify -q HEAD:rust > /dev/null 2>&1 || return 3
  if [ -z "$(git -C "$ROOT" status --porcelain --untracked-files=all -- rust 2> /dev/null)" ]; then
    git -C "$ROOT" rev-parse HEAD:rust
    return
  fi
  idx=$(mktemp -d "${TMPDIR:-/tmp}/og-tree.XXXXXX")
  { GIT_INDEX_FILE="${idx}/index" git -C "$ROOT" read-tree HEAD \
    && GIT_INDEX_FILE="${idx}/index" git -C "$ROOT" add -A -- rust \
    && GIT_INDEX_FILE="${idx}/index" git -C "$ROOT" write-tree --prefix=rust/; } || rc=$?
  rm -rf "$idx"
  return "$rc"
}

# embedded FILE: the hash stamped into a binary, or "unknown".
embedded() {
  local h
  h=$(LC_ALL=C grep -ao 'og-tree-hash:[0-9a-f]\{40,64\};' "$1" 2> /dev/null | head -n 1 || true)
  h=${h#og-tree-hash:}
  echo "${h%;}" | grep . || echo unknown
}

check() {
  local required=$REQUIRED expected f h missing="" stale=""
  [ "${1:-}" != --required ] || required=true
  expected=$(tree_hash 2> /dev/null) || expected=""
  for a in $ARCHES; do
    f="${BIN}/og-linux-${a}"
    if [ ! -x "$f" ]; then missing+=" ${a}"; continue; fi
    h=$(embedded "$f")
    [ -z "$expected" ] || [ "$h" = "$expected" ] || stale+=" bin/.og/og-linux-${a} (rust/ tree ${h})"
  done
  if [ -n "$stale" ]; then
    echo "ERROR: og is stale:${stale}, but rust/ here is tree ${expected}: run scripts/setup.sh --og"
    return 1
  fi
  if [ -n "$missing" ] && $required; then
    echo "ERROR: og is not built for${missing}: run scripts/setup.sh --og"
    return 1
  fi
  if [ -n "$missing" ]; then
    echo "NOTE: og is not built for${missing}; no step needs it yet (scripts/setup.sh --og builds it)"
  elif [ -z "$expected" ]; then
    echo "NOTE: og is built, but this is not a git checkout, so its rust/ tree cannot be checked"
  else
    echo "OK: og is built for ${ARCHES// /, } from this checkout's rust/ (tree ${expected})"
  fi
}

built() {
  local hs="" a
  for a in $ARCHES; do
    if [ -x "${BIN}/og-linux-${a}" ]; then hs+="${a}=$(embedded "${BIN}/og-linux-${a}") "; fi
  done
  hs=${hs% }
  if [ -z "$hs" ]; then echo "not built"; return; fi
  # One hash when both arches agree, as they do after setup.sh --og.
  local -a h
  read -r -a h <<< "$hs"
  if [ "${#h[@]}" -eq 2 ] && [ "${h[0]#*=}" = "${h[1]#*=}" ]; then echo "${h[0]#*=}"; else echo "$hs"; fi
}

# put SRC ARCH: move a binary into bin/.og in one rename.
put() {
  mkdir -p "$BIN"
  cp "$1" "${BIN}/.og-linux-${2}.tmp"
  chmod 755 "${BIN}/.og-linux-${2}.tmp"
  mv -f "${BIN}/.og-linux-${2}.tmp" "${BIN}/og-linux-${2}"
}

build() {
  local hash cache t0 a out targets=()
  hash=$(tree_hash) || { echo "ERROR: ${ROOT} is not a git checkout: og is stamped with the git tree of rust/, so it builds only in one" >&2; return 1; }
  # shellcheck source=../versions.env
  . "${ROOT}/versions.env"
  : "${OG_BUILDER:?OG_BUILDER is not set in versions.env}"
  # Build cache and output under rust/target (ignored by git, so outside the
  # tree hash); the sources go in read-only.
  for a in $ARCHES; do targets+=(--target "$(triple "$a")"); done
  cache="${ROOT}/rust/target/og-build"
  mkdir -p "$cache"
  echo "Building og ${hash} for ${ARCHES} in ${OG_BUILDER%@*} ..."
  t0=$(date +%s)
  "$ENGINE" run --rm --user "$(id -u):$(id -g)" \
    -e HOME=/cache/home -e CARGO_HOME=/cache/cargo -e CARGO_TARGET_DIR=/cache/target -e OG_TREE_HASH="$hash" \
    -v "${ROOT}/rust:/src:ro" -v "${cache}:/cache" -w /src "$OG_BUILDER" \
    cargo zigbuild --release --locked -p og "${targets[@]}"
  for a in $ARCHES; do
    out="${cache}/target/$(triple "$a")/release/og"
    [ "$(embedded "$out")" = "$hash" ] || { echo "ERROR: ${out} does not carry tree ${hash}" >&2; return 1; }
    put "$out" "$a"
  done
  echo "Built bin/.og/og-linux-{${ARCHES// /,}} in $(($(date +%s) - t0)) s"
}

sha256() {
  if command -v sha256sum > /dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

get() {  # get URL DEST: 0 when downloaded, 1 when not
  if command -v curl > /dev/null 2>&1; then curl -fsL --retry 3 -o "$2" "$1"; else wget -q -O "$2" "$1"; fi
}

# fetch TAG: 0 when the release's binaries are installed, 2 when the release
# has none (build instead), 1 when a download does not check out.
fetch() {
  local tag=$1 tmp want hash a f rc=0 attest=0
  hash=$(tree_hash)
  # gh attestation verify needs a login: without one (gh missing, or never
  # logged in) only SHA256SUMS is checked. gh auth token, not gh auth status,
  # which can hang on enterprise hosts.
  if command -v gh > /dev/null 2>&1 && gh auth token > /dev/null 2>&1; then attest=1; fi
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/og-fetch.XXXXXX")
  # OG_RELEASE_BASE: another server with the same layout (tests/test_og_check.sh).
  local base="${OG_RELEASE_BASE:-https://github.com/${REPO}/releases/download}/${tag}"
  if ! get "${base}/SHA256SUMS" "${tmp}/SHA256SUMS"; then
    rm -rf "$tmp"
    return 2
  fi
  for a in $ARCHES; do
    f="og-linux-${a}"
    want=$(awk -v f="$f" '$2 == f || $2 == "*" f {print $1}' "${tmp}/SHA256SUMS")
    if [ -z "$want" ] || ! get "${base}/${f}" "${tmp}/${f}"; then echo "ERROR: release ${tag} lists no ${f}" >&2; rc=1; break; fi
    if [ "$(sha256 "${tmp}/${f}")" != "$want" ]; then echo "ERROR: ${f} of release ${tag} does not match its SHA256SUMS" >&2; rc=1; break; fi
    if [ "$attest" = 1 ]; then
      gh attestation verify "${tmp}/${f}" --repo "$REPO" > /dev/null || { echo "ERROR: ${f} of release ${tag} has no valid build attestation" >&2; rc=1; break; }
    fi
    if [ "$(embedded "${tmp}/${f}")" != "$hash" ]; then echo "ERROR: ${f} of release ${tag} was not built from this rust/ (tree ${hash})" >&2; rc=1; break; fi
  done
  if [ "$rc" -eq 0 ]; then
    for a in $ARCHES; do put "${tmp}/og-linux-${a}" "$a"; done
    [ "$attest" = 1 ] || echo "NOTE: gh is not installed or not logged in, so the build attestation was not checked (SHA256SUMS was)"
    echo "Installed og ${hash} from release ${tag}"
  fi
  rm -rf "$tmp"
  return "$rc"
}

install() {
  local tag rc=0 hash
  hash=$(tree_hash) || { echo "ERROR: ${ROOT} is not a git checkout; clone the repository to install og" >&2; return 1; }
  tag=$(git -C "$ROOT" describe --exact-match --tags --match 'v[0-9]*.[0-9]*.[0-9]*' HEAD 2> /dev/null || true)
  if [[ "$tag" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] && [ "$hash" = "$(git -C "$ROOT" rev-parse HEAD:rust)" ]; then
    echo "Release ${tag}: downloading its og binaries"
    fetch "$tag" || rc=$?
    [ "$rc" -ne 0 ] || return 0
    [ "$rc" -eq 2 ] || return 1
    echo "Release ${tag} has no og binaries: building them"
  fi
  build
}

case "${1:-}" in
  install) install ;;
  build) build ;;
  check) shift; check "$@" ;;
  tree-hash) tree_hash ;;
  built) built ;;
  *) sed -n '2,/^set -euo/p' "$0" | sed -e '$d' -e 's/^# \{0,1\}//' >&2; exit 2 ;;
esac
