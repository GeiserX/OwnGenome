#!/usr/bin/env bash
# check-release.sh: the checks a release must pass before it is published.
#
#   check-release.sh tag TAG
#       TAG (vMAJOR.MINOR.PATCH) must equal "v" + the manifest version that
#       nextflow.config declares at TAG, so `nextflow run -r TAG` reports the
#       version the release is named after. The release notes are the tag's
#       CHANGELOG.md entry; scripts/ci/release-notes.sh checks and writes them.
#
#   check-release.sh checks RUNS.TSV [--ignore NAME]...
#       RUNS.TSV is the tagged commit's check runs, one per line:
#       name<TAB>status<TAB>conclusion (as `gh api .../commits/SHA/check-runs`
#       gives them). Every run must have finished, and finished as success,
#       skipped or neutral; Guard's ci-ok, which runs on every push to main,
#       must be among them. NAME (the release job itself) is left out.
#
#   check-release.sh declared
#       The working tree: nextflow.config declares exactly one manifest
#       version, of the form MAJOR.MINOR.PATCH. Guard runs this through
#       tests/test_release.sh.
#
#   check-release.sh dist DIR
#       DIR is what a release uploads: the four og binaries,
#       THIRD-PARTY-LICENSES.txt (rust/og.sh build) and SHA256SUMS. Each file
#       must be there, not empty, listed once in SHA256SUMS and match its
#       line. The licence file must list its crates and carry the notices of
#       the Rust standard library, musl, htslib and htscodecs (hts-sys's
#       htslib/LICENSE and htslib/htscodecs/LICENSE.md). Prints each file
#       with its checksum.
#
# --repo DIR runs against another git work tree (the tests plant one); the
# default is the checkout this script is in.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
if [ "${1:-}" = --repo ]; then ROOT=${2:?--repo needs a directory}; shift 2; fi
MODE=${1:-}
shift || true

die() { echo "ERROR: $*" >&2; exit 1; }
semver='^v([0-9]+)\.([0-9]+)\.([0-9]+)$'

# manifest_version: the version = '...' line inside the manifest { } block of
# nextflow.config on stdin. Exactly one, or an error.
manifest_version() {
  local v
  v=$(awk '
    /^[[:space:]]*manifest[[:space:]]*\{/ {inside = 1; next}
    inside && /^[[:space:]]*\}/ {inside = 0}
    inside && /^[[:space:]]*version[[:space:]]*=/ {
      s = $0; sub(/^[^=]*=[[:space:]]*/, "", s); sub(/[[:space:]]*(\/\/.*)?$/, "", s)
      gsub(/^["\047]|["\047]$/, "", s); print s
    }')
  [ -n "$v" ] || die "nextflow.config declares no version in its manifest block"
  [ "$(printf '%s\n' "$v" | wc -l)" -eq 1 ] || die "nextflow.config declares more than one manifest version"
  printf '%s\n' "$v"
}

case "$MODE" in
  tag)
    TAG=${1:?Usage: $0 tag TAG}
    [ $# -eq 1 ] || die "tag takes one argument, TAG; the release notes come from scripts/ci/release-notes.sh"
    [[ "$TAG" =~ $semver ]] || die "the tag is not of the form vMAJOR.MINOR.PATCH"
    git -C "$ROOT" rev-parse -q --verify "refs/tags/${TAG}" >/dev/null || die "tag ${TAG} does not exist in this repository"
    declared=$(git -C "$ROOT" show "${TAG}:nextflow.config" | manifest_version)
    [ "v${declared}" = "$TAG" ] \
      || die "tag ${TAG} does not match the version nextflow.config declares at that commit (${declared}): bump manifest.version in a pull request, then tag the merged commit"
    echo "OK: ${TAG} matches manifest.version"
    ;;
  checks)
    RUNS=${1:?Usage: $0 checks RUNS.TSV [--ignore NAME]...}
    shift
    ignore=()
    while [ $# -gt 0 ]; do
      [ "$1" = --ignore ] || die "unknown option: $1"
      ignore+=("${2:?--ignore needs a name}")
      shift 2
    done
    [ -f "$RUNS" ] || die "no check-run list at ${RUNS}"
    n=0 bad=0 cios=0
    while IFS=$'\t' read -r name status conclusion || [ -n "$name" ]; do
      [ -n "$name" ] || continue
      skip=0
      for i in "${ignore[@]+"${ignore[@]}"}"; do [ "$name" != "$i" ] || skip=1; done
      [ "$skip" = 0 ] || continue
      n=$((n + 1))
      if [ "$status" != completed ]; then
        echo "NOT FINISHED: ${name} (${status})"; bad=1
      else
        case "$conclusion" in
          success|skipped|neutral) [ "$name" != ci-ok ] || [ "$conclusion" != success ] || cios=1 ;;
          *) echo "NOT GREEN: ${name} (${conclusion:-no conclusion})"; bad=1 ;;
        esac
      fi
    done < "$RUNS"
    [ "$n" -gt 0 ] || die "the tagged commit has no check runs: tag a commit that main's CI ran on"
    [ "$cios" = 1 ] || die "Guard's ci-ok did not succeed on the tagged commit: tag a commit that main's CI ran on"
    [ "$bad" = 0 ] || die "the tagged commit's checks are not all green: wait for them or fix main, then run the release again"
    echo "OK: all ${n} check runs of the tagged commit finished green"
    ;;
  declared)
    declared=$(manifest_version < "${ROOT}/nextflow.config")
    [[ "v${declared}" =~ $semver ]] || die "nextflow.config declares version '${declared}', which is not MAJOR.MINOR.PATCH"
    echo "OK: declared version ${declared}"
    ;;
  dist)
    DIST=${1:?Usage: $0 dist DIR}
    sums="${DIST}/SHA256SUMS"
    [ -s "$sums" ] || die "${DIST} has no SHA256SUMS"
    for f in og-linux-x86_64 og-linux-aarch64 og-darwin-arm64 og-darwin-x86_64 THIRD-PARTY-LICENSES.txt; do
      [ -s "${DIST}/${f}" ] || die "${DIST} has no ${f}, or it is empty"
      want=$(awk -v f="$f" '$2 == f || $2 == "*" f {print $1}' "$sums")
      [ "$(printf '%s\n' "$want" | grep -c .)" -eq 1 ] || die "SHA256SUMS must list ${f} once"
      if command -v sha256sum > /dev/null 2>&1; then got=$(sha256sum "${DIST}/${f}"); else got=$(shasum -a 256 "${DIST}/${f}"); fi
      [ "${got%% *}" = "$want" ] || die "${f} does not match its SHA256SUMS line"
      echo "${want}  ${f}"
    done
    lic="${DIST}/THIRD-PARTY-LICENSES.txt"
    grep -qE '^Crates \([1-9][0-9]*\)$' "$lic" || die "THIRD-PARTY-LICENSES.txt lists no crates"
    grep -qE '^For Rust standard library: ' "$lic" || die "THIRD-PARTY-LICENSES.txt does not carry the Rust standard library's notices"
    grep -qE '^For musl libc: COPYRIGHT$' "$lic" || die "THIRD-PARTY-LICENSES.txt does not carry musl's licence"
    grep -qE '^For hts-sys [^:]+: htslib/LICENSE$' "$lic" || die "THIRD-PARTY-LICENSES.txt does not carry htslib's licence (hts-sys htslib/LICENSE)"
    grep -qE '^For hts-sys [^:]+: htslib/htscodecs/LICENSE\.md$' "$lic" || die "THIRD-PARTY-LICENSES.txt does not carry htscodecs' licence (hts-sys htslib/htscodecs/LICENSE.md)"
    echo "OK: the og binaries and THIRD-PARTY-LICENSES.txt match SHA256SUMS, and the notices of the Rust standard library, musl, htslib and htscodecs are in it"
    ;;
  *)
    die "usage: $0 [--repo DIR] {tag TAG | checks RUNS.TSV [--ignore NAME]... | declared | dist DIR}"
    ;;
esac
