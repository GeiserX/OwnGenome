#!/usr/bin/env bash
# og.sh: how the og binary reaches bin/.og, and the check that it is the one
# this checkout's rust/ makes. bin/og (the shim tasks call) runs
# bin/.og/og-linux-<arch>; both arches are kept, because an Apple Silicon
# host runs amd64 images under emulation, where uname -m says x86_64. On a
# macOS host (the live site, scripts/owngenome-live.sh) bin/og runs
# bin/.og/og-darwin-<arm64|x86_64> instead: built without the cargo feature
# `bam`, because that host never reads a BAM (htslib is still linked,
# through og-summary and og-setup).
#
#   rust/og.sh install     scripts/setup.sh --og: on a release tag (vX.Y.Z,
#                          clean rust/) download the release's binaries and
#                          check them against SHA256SUMS (and the build
#                          attestation when gh is logged in); otherwise, or
#                          when the release has none, build them
#   rust/og.sh build       build all four from rust/ in OG_BUILDER (versions.env),
#                          and bin/.og/THIRD-PARTY-LICENSES.txt: the licences of
#                          every crate and C library they link, which a
#                          release ships beside them (og dev licences notices)
#   rust/og.sh check [--required]
#                          one line, OK, NOTE or ERROR; exit 1 on ERROR: a
#                          binary built from another rust/, or (with
#                          --required, or once REQUIRED below is true) no
#                          Linux pair. The macOS pair is optional, since no
#                          pipeline step runs it, but a stale one is refused
#   rust/og.sh tree-hash   the source hash of rust/ and its embedded table
#   rust/og.sh built       the tree hash the binaries were built from, for
#                          run_manifest.tsv ("not built" when there are none)
#
# Env: OG_ARCHES         the Linux arches to build, fetch and require (default
#                        "x86_64 aarch64"); one of them must be the arch the
#                        builder runs on, which writes the licence notices
#      OG_DARWIN_ARCHES  the macOS arches to build and fetch (default
#                        "arm64 x86_64"; empty for none)
#      The E2E workflow sets OG_ARCHES=x86_64 and OG_DARWIN_ARCHES= for the
#      one binary its runners run; a release and Guard build all four.
#
# The source stamp hashes two git object IDs: rust/'s tree as it is on
# disk and the blob of assets/imprinted_regions.tsv. The latter is embedded
# at compile time, so an edit to it must also make a binary stale. The
# legacy tree-hash verb and og-tree-hash:<hash>; stamp keep their format.
# Checks read the stamp without executing a binary for another platform.
set -euo pipefail

# Steps 04, 05, 07, 09b, 21, 25, 32 and 33 (and the modules of 09b, 21 and 32:
# STRANGER_RFC1, CYRIUS and PYPGX) call og, and stop with bin/og's message
# without it. A run is not refused up front yet: the fake-docker and Nextflow
# stub tests start runs without og, so they need it first. Until then a
# checkout without the binary starts, the check names the steps that need it,
# and only a stale binary is refused.
REQUIRED=false

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BIN="${ROOT}/bin/.og"
ARCHES=${OG_ARCHES:-x86_64 aarch64}
# macOS arches as uname -m names them, so bin/og finds og-darwin-$(uname -m).
DARWIN_ARCHES=${OG_DARWIN_ARCHES-arm64 x86_64}
for a in $ARCHES; do
  case $a in x86_64|aarch64) ;; *) echo "ERROR: OG_ARCHES may name x86_64 and aarch64, not '${a}'" >&2; exit 2 ;; esac
done
for a in $DARWIN_ARCHES; do
  case $a in arm64|x86_64) ;; *) echo "ERROR: OG_DARWIN_ARCHES may name arm64 and x86_64, not '${a}'" >&2; exit 2 ;; esac
done
REPO="GeiserX/OwnGenome"
ENGINE=${CONTAINER_ENGINE:-docker}

triple() { echo "${1}-unknown-linux-musl"; }
darwin_triple() { case $1 in arm64) echo aarch64-apple-darwin ;; *) echo "${1}-apple-darwin" ;; esac; }

# names: every binary bin/.og holds, Linux first; a release ships the same.
names() {
  local a
  for a in $ARCHES; do echo "og-linux-${a}"; done
  for a in $DARWIN_ARCHES; do echo "og-darwin-${a}"; done
}

# tree_hash: exit 3 when this is not a git checkout, 1 when git fails.
rust_tree_hash() {
  local idx rc=0
  git -C "$ROOT" rev-parse --verify -q HEAD:rust > /dev/null 2>&1 || return 3
  if [ -z "$(git -C "$ROOT" status --porcelain --untracked-files=all -- rust 2> /dev/null)" ]; then
    git -C "$ROOT" rev-parse HEAD:rust
    return
  fi
  idx=$(mktemp -d "${TMPDIR:-/tmp}/og-tree.XXXXXX") || return 1
  { GIT_INDEX_FILE="${idx}/index" git -C "$ROOT" read-tree HEAD \
    && GIT_INDEX_FILE="${idx}/index" git -C "$ROOT" add -A -- rust \
    && GIT_INDEX_FILE="${idx}/index" git -C "$ROOT" write-tree --prefix=rust/; } || rc=$?
  rm -rf "$idx"
  [ "$rc" -eq 0 ] || return 1
}

# Hash the two source object IDs without writing to the checkout's object
# store. Both committed and uncommitted table edits change the stamp.
source_hash() { printf '%s\n%s\n' "$1" "$2" | git -C "$ROOT" hash-object --stdin; }
tree_hash() {
  local rust asset
  rust=$(rust_tree_hash) || return $?
  asset=$(git -C "$ROOT" hash-object -- assets/imprinted_regions.tsv) || return 1
  source_hash "$rust" "$asset"
}
committed_hash() {
  local rust asset
  rust=$(git -C "$ROOT" rev-parse HEAD:rust) || return 1
  asset=$(git -C "$ROOT" rev-parse HEAD:assets/imprinted_regions.tsv) || return 1
  source_hash "$rust" "$asset"
}

# embedded FILE: the hash stamped into a binary, or "unknown".
embedded() {
  local h
  h=$(LC_ALL=C grep -ao 'og-tree-hash:[0-9a-f]\{40,64\};' "$1" 2> /dev/null | head -n 1 || true)
  h=${h#og-tree-hash:}
  echo "${h%;}" | grep . || echo unknown
}

check() {
  local required=$REQUIRED expected f h missing="" stale="" trc=0
  [ "${1:-}" != --required ] || required=true
  expected=$(tree_hash 2> /dev/null) || trc=$?
  if [ "$trc" -eq 3 ]; then
    expected=""
  elif [ "$trc" -ne 0 ]; then
    echo "ERROR: the git tree hash of ${ROOT}/rust could not be computed, so og cannot be checked (is TMPDIR writable?)"
    return 1
  fi
  local mac=""
  for a in $ARCHES; do
    f="${BIN}/og-linux-${a}"
    if [ ! -x "$f" ]; then missing+=" ${a}"; continue; fi
    h=$(embedded "$f")
    [ -z "$expected" ] || [ "$h" = "$expected" ] || stale+=" bin/.og/og-linux-${a} (rust/ tree ${h})"
  done
  # The macOS pair: optional, but never stale.
  for a in $DARWIN_ARCHES; do
    f="${BIN}/og-darwin-${a}"
    [ -x "$f" ] || continue
    mac+=", ${a}"
    h=$(embedded "$f")
    [ -z "$expected" ] || [ "$h" = "$expected" ] || stale+=" bin/.og/og-darwin-${a} (rust/ tree ${h})"
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
    echo "NOTE: og is not built for${missing}; steps 04, 05, 07, 09b, 21, 25, 32 and 33 stop without it: run scripts/setup.sh --og"
  elif [ -z "$expected" ]; then
    echo "NOTE: og is built, but this is not a git checkout, so its rust/ tree cannot be checked"
  elif [ -n "$mac" ]; then
    echo "OK: og is built for ${ARCHES// /, } and macOS ${mac#, } from this checkout's rust/ (tree ${expected})"
  else
    echo "OK: og is built for ${ARCHES// /, } from this checkout's rust/ (tree ${expected})"
  fi
}

built() {
  local hs="" a h n=0 first="" same=true
  for a in $ARCHES; do
    n=$((n + 1))
    [ -x "${BIN}/og-linux-${a}" ] || { same=false; continue; }
    h=$(embedded "${BIN}/og-linux-${a}")
    hs+="${a}=${h} "
    if [ -z "$first" ]; then first=$h; elif [ "$h" != "$first" ]; then same=false; fi
  done
  hs=${hs% }
  if [ -z "$hs" ]; then echo "not built"; return; fi
  # One hash when every arch of ARCHES is built and they agree, as after
  # setup.sh --og.
  if $same && [ "$n" -gt 0 ]; then echo "$first"; else echo "$hs"; fi
}

# put SRC...: one binary per name of names(), in that order, into bin/.og.
# All are copied beside their targets first; only when every copy worked are
# they renamed into place, so a failed copy (a full disk) leaves the old set.
put() {
  local n m i=0 srcs=("$@")
  mkdir -p "$BIN" || return 1
  for n in $(names); do
    if ! { cp "${srcs[i]}" "${BIN}/.${n}.tmp" && chmod 755 "${BIN}/.${n}.tmp"; }; then
      echo "ERROR: could not copy ${n} into ${BIN}" >&2
      for m in $(names); do rm -f "${BIN}/.${m}.tmp"; done
      return 1
    fi
    i=$((i + 1))
  done
  for n in $(names); do
    mv -f "${BIN}/.${n}.tmp" "${BIN}/${n}" || { echo "ERROR: could not install ${BIN}/${n}" >&2; return 1; }
  done
}

build() {
  local hash cache t0 a out targets=() mac_targets=() outs=()
  hash=$(tree_hash) || { echo "ERROR: ${ROOT} is not a git checkout: og is stamped with the git tree of rust/, so it builds only in one" >&2; return 1; }
  # shellcheck source=../versions.env
  . "${ROOT}/versions.env"
  : "${OG_BUILDER:?OG_BUILDER is not set in versions.env}"
  # Build cache and output under rust/target (ignored by git, so outside the
  # tree hash); the sources go in read-only.
  for a in $ARCHES; do targets+=(--target "$(triple "$a")"); done
  for a in $DARWIN_ARCHES; do mac_targets+=(--target "$(darwin_triple "$a")"); done
  cache="${ROOT}/rust/target/og-build"
  mkdir -p "$cache"
  echo "Building og ${hash} for ${ARCHES} and macOS ${DARWIN_ARCHES:-(none)} in ${OG_BUILDER%@*} ..."
  t0=$(date +%s)
  # Two cargo runs, because the feature set is per run: Linux with BAM (the
  # default features), macOS without it (zig links against the image's
  # macOS SDK).
  zigbuild() {
    "$ENGINE" run --rm --user "$(id -u):$(id -g)" \
      -e HOME=/cache/home -e CARGO_HOME=/cache/cargo -e CARGO_TARGET_DIR=/cache/target -e OG_TREE_HASH="$hash" \
      -v "${ROOT}/rust:/src:ro" \
      -v "${ROOT}/assets/imprinted_regions.tsv:/assets/imprinted_regions.tsv:ro" -v "${cache}:/cache" -w /src "$OG_BUILDER" \
      cargo zigbuild --release --locked -p og "$@"
  }
  zigbuild "${targets[@]}" || return
  # No macOS target: skip the run, which would otherwise build for the host.
  if [ "${#mac_targets[@]}" -gt 0 ]; then zigbuild --no-default-features "${mac_targets[@]}" || return; fi
  for a in $ARCHES; do outs+=("${cache}/target/$(triple "$a")/release/og"); done
  for a in $DARWIN_ARCHES; do outs+=("${cache}/target/$(darwin_triple "$a")/release/og"); done
  for out in "${outs[@]}"; do
    [ "$(embedded "$out")" = "$hash" ] || { echo "ERROR: ${out} does not carry tree ${hash}" >&2; return 1; }
  done
  # The notices of everything the four link, written by the new og itself in
  # the builder, where cargo and the crate sources are: the Linux binary of
  # the builder's own arch, with the targets and feature sets built above.
  # Two parts are no crate, so the builder names their notices: the Rust
  # standard library (the toolchain's COPYRIGHT-library.html) and musl, the
  # C library of the Linux pair (musl's COPYRIGHT, which zig ships).
  local lic=() notices="${cache}/THIRD-PARTY-LICENSES.txt"
  for a in $ARCHES; do lic+=(--target "$(triple "$a")"); done
  for a in $DARWIN_ARCHES; do lic+=(--target-no-default "$(darwin_triple "$a")"); done
  rm -f "$notices"
  # shellcheck disable=SC2016 # the sh -c script expands in the builder
  "$ENGINE" run --rm --user "$(id -u):$(id -g)" \
    -e HOME=/cache/home -e CARGO_HOME=/cache/cargo -e CARGO_TARGET_DIR=/cache/target \
    -v "${ROOT}/rust:/src:ro" -v "${cache}:/cache" -w /src "$OG_BUILDER" \
    sh -c 'sr=$(rustc --print sysroot) || exit
      set -- "$@" --extra "Rust standard library=${sr}/share/doc/rust/COPYRIGHT-library.html"
      for m in /usr/local/zig*/lib/libc/musl/COPYRIGHT; do set -- "$@" --extra "musl libc=${m}"; done
      exec "/cache/target/$(uname -m)-unknown-linux-musl/release/og" "$@"' og \
    dev licences notices "${lic[@]}" --out /cache/THIRD-PARTY-LICENSES.txt || return
  [ -s "$notices" ] || { echo "ERROR: og dev licences notices wrote no ${notices}" >&2; return 1; }
  put "${outs[@]}" || return 1
  { cp "$notices" "${BIN}/.THIRD-PARTY-LICENSES.txt.tmp" && mv -f "${BIN}/.THIRD-PARTY-LICENSES.txt.tmp" "${BIN}/THIRD-PARTY-LICENSES.txt"; } \
    || { echo "ERROR: could not install ${BIN}/THIRD-PARTY-LICENSES.txt" >&2; return 1; }
  echo "Built $(names | sed "s|^|bin/.og/|" | paste -sd" " -) and bin/.og/THIRD-PARTY-LICENSES.txt in $(($(date +%s) - t0)) s"
}

sha256() {
  if command -v sha256sum > /dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

get() {  # get URL DEST: 0 when downloaded, 2 when the server says 404, 1 on any other error
  local code
  if command -v curl > /dev/null 2>&1; then
    code=$(curl -sL --retry 3 -o "$2" -w '%{http_code}' "$1") || return 1
    case "$code" in 2??) return 0 ;; 404) return 2 ;; *) return 1 ;; esac
  fi
  wget -q -O "$2" "$1" && return 0
  # wget exits 8 on an error response from the server; that is all it tells.
  [ "$?" -eq 8 ] && return 2 || return 1
}

# fetch TAG: 0 when the release's binaries are installed, 2 when the release
# has none (SHA256SUMS is not there: build instead), 1 when a download fails
# or does not check out.
fetch() {
  local tag=$1 tmp want hash f rc=0 grc=0 attest=0
  hash=$(tree_hash)
  # gh attestation verify needs a login: without one (gh missing, or never
  # logged in) only SHA256SUMS is checked. gh auth token, not gh auth status,
  # which can hang on enterprise hosts.
  if command -v gh > /dev/null 2>&1 && gh auth token > /dev/null 2>&1; then attest=1; fi
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/og-fetch.XXXXXX")
  # OG_RELEASE_BASE: another server with the same layout (tests/test_og_check.sh).
  local base="${OG_RELEASE_BASE:-https://github.com/${REPO}/releases/download}/${tag}"
  get "${base}/SHA256SUMS" "${tmp}/SHA256SUMS" || grc=$?
  if [ "$grc" -ne 0 ]; then
    rm -rf "$tmp"
    [ "$grc" -eq 2 ] && return 2
    echo "ERROR: could not download ${base}/SHA256SUMS; check the network and run scripts/setup.sh --og again" >&2
    return 1
  fi
  for f in $(names); do
    want=$(awk -v f="$f" '$2 == f || $2 == "*" f {print $1}' "${tmp}/SHA256SUMS")
    if [ -z "$want" ]; then echo "ERROR: release ${tag} lists no ${f}" >&2; rc=1; break; fi
    if ! get "${base}/${f}" "${tmp}/${f}"; then echo "ERROR: could not download ${f} of release ${tag}" >&2; rc=1; break; fi
    if [ "$(sha256 "${tmp}/${f}")" != "$want" ]; then echo "ERROR: ${f} of release ${tag} does not match its SHA256SUMS" >&2; rc=1; break; fi
    if [ "$attest" = 1 ]; then
      gh attestation verify "${tmp}/${f}" --repo "$REPO" > /dev/null || { echo "ERROR: ${f} of release ${tag} has no valid build attestation" >&2; rc=1; break; }
    fi
    if [ "$(embedded "${tmp}/${f}")" != "$hash" ]; then echo "ERROR: ${f} of release ${tag} was not built from this rust/ (tree ${hash})" >&2; rc=1; break; fi
  done
  local srcs=()
  for f in $(names); do srcs+=("${tmp}/${f}"); done
  if [ "$rc" -eq 0 ] && ! put "${srcs[@]}"; then rc=1; fi
  if [ "$rc" -eq 0 ]; then
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
  if [[ "$tag" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] && [ "$hash" = "$(committed_hash)" ]; then
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
