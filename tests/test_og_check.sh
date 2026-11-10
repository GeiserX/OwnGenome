#!/usr/bin/env bash
# rust/og.sh check on planted trees: the stale-binary guard of bin/.og.
# A copy of the checkout's scripts, bin/ and rust/ becomes its own git repo;
# fake binaries carry a tree hash the way the real build stamps it
# (og-tree-hash:<hash>;). Checked: no binary is a note (an error with
# --required); binaries of this rust/ pass; an edit to rust/, uncommitted,
# untracked or committed, makes them stale; a file git ignores does not; a
# binary stamped "unknown" is stale; a tree hash that cannot be computed is
# an error; built prints the hash for run_manifest.tsv; outside git the hash
# cannot be checked. Release downloads: checksum, attestation, a network
# error, a release without binaries, and a copy into bin/.og that fails. Then the callers:
# run-all.sh (with SKIP_VALIDATION=true) and validate-setup.sh refuse a stale
# binary with the setup command, and write_manifest.sh records the hash.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fails=0
ok()  { echo "PASS $*"; }
bad() { echo "FAIL $*"; fails=$((fails + 1)); }
expect() {  # expect NAME RC WANT_RC OUTPUT PATTERN
  if [ "$2" -eq "$3" ] && [[ "$4" == $5 ]]; then ok "$1"; else bad "$1: exit $2 (want $3), output: $4"; fi
}

C="${T}/repo"
mkdir -p "$C"
cp -R "${REPO}/scripts" "${REPO}/bin" "${REPO}/rust" "${REPO}/versions.env" "${REPO}/nextflow.config" "${REPO}/.gitignore" "$C/"
rm -rf "${C}/bin/.og" "${C}/rust/target"
g() { git -C "$C" -c user.name=t -c user.email=t@example.invalid "$@"; }
g init -q && g add -A && g commit -qm base
OG="${C}/rust/og.sh"
plant() {  # plant HASH [ARCH...]: fake binaries stamped with HASH
  local h=$1 a; shift
  [ $# -gt 0 ] || set -- x86_64 aarch64
  mkdir -p "${C}/bin/.og"
  for a in "$@"; do
    printf '\177ELF junk og-tree-hash:%s;attempt more junk' "$h" > "${C}/bin/.og/og-linux-${a}"
    chmod 755 "${C}/bin/.og/og-linux-${a}"
  done
}

h0=$(bash "$OG" tree-hash)
[ "$h0" = "$(g rev-parse HEAD:rust)" ] && ok "tree-hash of a clean rust/ is HEAD:rust" || bad "tree-hash ${h0} != HEAD:rust"

out=$(bash "$OG" check); rc=$?
expect "no binary: a note, exit 0" "$rc" 0 "$out" "NOTE: og is not built for x86_64 aarch64*setup.sh --og*"
out=$(bash "$OG" check --required); rc=$?
expect "no binary with --required: refused" "$rc" 1 "$out" "ERROR: og is not built for x86_64 aarch64: run scripts/setup.sh --og"
out=$(bash "$OG" built); [ "$out" = "not built" ] && ok "built without binaries: not built" || bad "built printed '${out}'"

plant "$h0"
out=$(bash "$OG" check); rc=$?
expect "binaries of this rust/: OK" "$rc" 0 "$out" "OK: og is built for x86_64, aarch64*${h0}*"
out=$(bash "$OG" built); [ "$out" = "$h0" ] && ok "built prints the one hash" || bad "built printed '${out}'"

plant "$h0" x86_64
rm -f "${C}/bin/.og/og-linux-aarch64"
out=$(bash "$OG" check --required); rc=$?
expect "one arch missing with --required: refused" "$rc" 1 "$out" "ERROR: og is not built for aarch64*"

plant "$h0"
mkdir -p "${C}/rust/target/debug" && echo junk > "${C}/rust/target/debug/x"
out=$(bash "$OG" check); rc=$?
expect "a file git ignores (rust/target) does not make it stale" "$rc" 0 "$out" "OK:*"

echo "// edited" >> "${C}/rust/og/src/main.rs"
h1=$(bash "$OG" tree-hash)
[ "$h1" != "$h0" ] && ok "an uncommitted edit changes the tree hash" || bad "tree hash did not change after an edit"
out=$(bash "$OG" check); rc=$?
expect "uncommitted edit to rust/: stale, refused" "$rc" 1 "$out" "ERROR: og is stale:*${h0}*${h1}*run scripts/setup.sh --og"
g add rust/og/src/main.rs && g commit -qm edit
[ "$(bash "$OG" tree-hash)" = "$h1" ] && [ "$h1" = "$(g rev-parse HEAD:rust)" ] \
  && ok "committing the edit keeps the hash the build would stamp" || bad "hash after commit differs"
out=$(bash "$OG" check); rc=$?
expect "committed edit to rust/: stale, refused" "$rc" 1 "$out" "ERROR: og is stale:*"
plant "$h1"
echo "new" > "${C}/rust/new-file.txt"
out=$(bash "$OG" check); rc=$?
expect "an untracked file in rust/: stale" "$rc" 1 "$out" "ERROR: og is stale:*"
rm -f "${C}/rust/new-file.txt"
out=$(bash "$OG" check); rc=$?
expect "rebuilt from the edited rust/: OK again" "$rc" 0 "$out" "OK:*${h1}*"

printf 'og-tree-hash:unknown;' > "${C}/bin/.og/og-linux-x86_64"
out=$(bash "$OG" check); rc=$?
expect "a binary stamped unknown (plain cargo build): stale" "$rc" 1 "$out" "ERROR: og is stale:*og-linux-x86_64 (rust/ tree unknown)*"
out=$(bash "$OG" built); [[ "$out" == "x86_64=unknown aarch64=${h1}" ]] && ok "built lists both hashes when they differ" || bad "built printed '${out}'"

# The tree hash of an edited rust/ cannot be computed (TMPDIR not writable):
# an error, never the note of a checkout without git.
echo "// edited" >> "${C}/rust/og/src/main.rs"
out=$(TMPDIR="${T}/no-such-dir" bash "$OG" check 2> /dev/null); rc=$?
expect "a tree hash that fails is an error" "$rc" 1 "$out" "ERROR: the git tree hash of*could not be computed*"
g checkout -q -- rust/og/src/main.rs

# The callers, on a stale binary.
plant "$h0"
mkdir -p "${T}/fake" "${T}/genome"
printf '#!/bin/sh\necho "openjdk version \\"17.0.1\\"" >&2\n' > "${T}/fake/java"
printf '#!/bin/sh\necho "fake nextflow: must not be reached" >&2\nexit 99\n' > "${T}/fake/nextflow"
chmod 755 "${T}/fake/java" "${T}/fake/nextflow"
out=$(PATH="${T}/fake:${PATH}" GENOME_DIR="${T}/genome" SKIP_VALIDATION=true bash "${C}/scripts/run-all.sh" s1 male 2>&1); rc=$?
expect "run-all.sh with SKIP_VALIDATION=true stops on a stale og" "$rc" 2 "$out" "*ERROR: og is stale:*run scripts/setup.sh --og*"
[[ "$out" != *"must not be reached"* ]] && ok "run-all.sh stopped before nextflow" || bad "run-all.sh reached nextflow"
out=$(GENOME_DIR="${T}/genome" bash "${C}/scripts/validate-setup.sh" 2>&1); rc=$?
expect "validate-setup.sh fails a stale og" "$rc" 1 "$out" "*[[]FAIL[]]*og is stale:*run scripts/setup.sh --og*"
plant "$h1"
out=$(GENOME_DIR="${T}/genome" bash "${C}/scripts/validate-setup.sh" 2>&1)
[[ "$out" == *"[OK]"*"og is built for x86_64, aarch64"* ]] && ok "validate-setup.sh passes a fresh og" || bad "validate-setup.sh: $(grep -i ' og ' <<< "$out")"
mkdir -p "${T}/genome/s1"
CONTAINER_ENGINE=true GENOME_DIR="${T}/genome" bash "${C}/bin/write_manifest.sh" s1 test > /dev/null 2>&1
row=$(awk -F'\t' '$1 == "run" && $2 == "og_tree_hash" {print $3}' "${T}/genome/s1/run_manifest.tsv" 2> /dev/null)
[ "$row" = "$h1" ] && ok "run_manifest.tsv records og_tree_hash" || bad "run_manifest.tsv og_tree_hash row: '${row}'"

# A release tag: the binaries come from the release, checked against
# SHA256SUMS and the attestation (a fake gh here), and are refused when
# either does not match; a gh that is not logged in skips the attestation. OG_RELEASE_BASE points at a local server.
if command -v python3 > /dev/null 2>&1 && command -v curl > /dev/null 2>&1; then
  rm -rf "${C}/bin/.og"
  g tag v9.9.9
  S="${T}/srv/v9.9.9"
  mkdir -p "$S"
  for a in x86_64 aarch64; do printf 'og-tree-hash:%s;' "$h1" > "${S}/og-linux-${a}"; done
  (cd "$S" && { sha256sum og-linux-* 2> /dev/null || shasum -a 256 og-linux-*; } > SHA256SUMS)
  port=$((20000 + RANDOM % 20000))
  (cd "${T}/srv" && exec python3 -m http.server "$port" --bind 127.0.0.1 > /dev/null 2>&1) &
  srv=$!
  for _ in $(seq 50); do curl -fs "http://127.0.0.1:${port}/v9.9.9/SHA256SUMS" > /dev/null 2>&1 && break; sleep 0.1; done
  printf '#!/bin/sh
exit 0
' > "${T}/fake/gh"; chmod 755 "${T}/fake/gh"
  out=$(PATH="${T}/fake:${PATH}" OG_RELEASE_BASE="http://127.0.0.1:${port}" bash "$OG" install 2>&1); rc=$?
  expect "a release tag installs the release's binaries" "$rc" 0 "$out" "*Installed og ${h1} from release v9.9.9*"
  [ "$(bash "$OG" built)" = "$h1" ] && ok "the installed binaries carry the tree" || bad "built after install: $(bash "$OG" built)"
  rm -rf "${C}/bin/.og"
  printf '#!/bin/sh
[ "$1" = auth ] && exit 0
echo "no attestation" >&2
exit 1
' > "${T}/fake/gh"
  out=$(PATH="${T}/fake:${PATH}" OG_RELEASE_BASE="http://127.0.0.1:${port}" bash "$OG" install 2>&1); rc=$?
  expect "a binary without a valid attestation is refused" "$rc" 1 "$out" "*has no valid build attestation*"
  # gh installed but not logged in: attestation verify cannot run, so only
  # SHA256SUMS is checked, with a note, never a refusal that reads as tampering.
  printf '#!/bin/sh
echo "To get started with GitHub CLI, please run: gh auth login" >&2
exit 4
' > "${T}/fake/gh"
  out=$(PATH="${T}/fake:${PATH}" OG_RELEASE_BASE="http://127.0.0.1:${port}" bash "$OG" install 2>&1); rc=$?
  expect "gh not logged in: installed, attestation not checked" "$rc" 0 "$out" "*attestation was not checked*Installed og ${h1} from release v9.9.9*"
  rm -rf "${C}/bin/.og"
  printf '#!/bin/sh
exit 0
' > "${T}/fake/gh"
  printf 'tampered' >> "${S}/og-linux-aarch64"
  out=$(PATH="${T}/fake:${PATH}" OG_RELEASE_BASE="http://127.0.0.1:${port}" bash "$OG" install 2>&1); rc=$?
  expect "a binary that does not match SHA256SUMS is refused" "$rc" 1 "$out" "*og-linux-aarch64 of release v9.9.9 does not match its SHA256SUMS*"
  [ ! -e "${C}/bin/.og/og-linux-x86_64" ] && ok "a refused release installs nothing" || bad "a refused release left a binary"
  # A copy that fails (bin/.og not writable) keeps the old pair and fails.
  (cd "$S" && for a in x86_64 aarch64; do printf 'og-tree-hash:%s;' "$h1" > "og-linux-${a}"; done \
    && { sha256sum og-linux-* 2> /dev/null || shasum -a 256 og-linux-*; } > SHA256SUMS)
  plant "$h0"
  chmod 555 "${C}/bin/.og"
  out=$(PATH="${T}/fake:${PATH}" OG_RELEASE_BASE="http://127.0.0.1:${port}" bash "$OG" install 2>&1); rc=$?
  chmod 755 "${C}/bin/.og"
  expect "a copy into bin/.og that fails is an error" "$rc" 1 "$out" "*could not copy og-linux-x86_64*"
  [[ "$out" != *"Installed og"* ]] && [ "$(bash "$OG" built)" = "$h0" ] && [ -z "$(ls -A "${C}/bin/.og" | grep tmp)" ] \
    && ok "a failed copy keeps the old pair and leaves no temporary file" || bad "after a failed copy: $(ls -A "${C}/bin/.og"), built $(bash "$OG" built)"
  rm -rf "${C}/bin/.og"
  # A network error on SHA256SUMS is an error, never a fallback to a build.
  printf '#!/bin/sh\necho "fake engine: must not be reached" >&2\nexit 99\n' > "${T}/fake/engine"; chmod 755 "${T}/fake/engine"
  out=$(PATH="${T}/fake:${PATH}" CONTAINER_ENGINE=engine OG_RELEASE_BASE="http://127.0.0.1:1" bash "$OG" install 2>&1); rc=$?
  expect "a release download that fails is an error, not a build" "$rc" 1 "$out" "*could not download http://127.0.0.1:1/v9.9.9/SHA256SUMS*"
  [[ "$out" != *"must not be reached"* ]] && ok "no build after a failed download" || bad "a failed download started a build"
  # A release without SHA256SUMS (404) has no binaries: build instead.
  out=$(PATH="${T}/fake:${PATH}" CONTAINER_ENGINE=engine OG_RELEASE_BASE="http://127.0.0.1:${port}/none" bash "$OG" install 2>&1); rc=$?
  expect "a release without binaries falls back to a build" "$rc" 99 "$out" "*has no og binaries: building them*must not be reached*"
  kill "$srv" 2> /dev/null; wait "$srv" 2> /dev/null
  g tag -d v9.9.9 > /dev/null
else
  echo "SKIP release download cases: no python3 or curl"
fi

# Outside git the hash cannot be checked: a note, never a silent OK.
N="${T}/nogit"
plant "$h1"
cp -R "$C" "$N" && rm -rf "${N}/.git"
out=$(bash "${N}/rust/og.sh" check); rc=$?
expect "binaries outside a git checkout: a note" "$rc" 0 "$out" "NOTE: og is built, but this is not a git checkout*"

[ "$fails" -eq 0 ] && echo "All og check cases passed" || echo "${fails} og check case(s) failed"
exit "$fails"
