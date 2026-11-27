#!/usr/bin/env bash
# scripts/ci/release-notes.sh on planted trees, then on the real UPGRADING.md.
#   A tag's notes are its "## X.Y.Z" CHANGELOG.md entry, as the tagged commit
#   holds it: the entry ends at the next "## " heading, a heading that only
#   starts with the version (## 1.2.30) is another entry, and relative links
#   point into the tagged tree on GitHub.
#   A tag with no entry, or an empty one, fails and prints no notes.
#   A major (X.0.0) needs "## Upgrading to X.0" in UPGRADING.md, and its
#   notes open with a link to that section; without the section it fails.
#   A REF argument reads another commit, so notes can be checked before the
#   tag exists. The real UPGRADING.md has the 1.0 section the 1.0.0 tag needs.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
NOTES="${REPO}/scripts/ci/release-notes.sh"
T=$(mktemp -d)
E=$(mktemp)
trap 'rm -rf "$T" "$E"' EXIT
fails=0
ok()  { echo "PASS $*"; }
bad() { echo "FAIL $*"; fails=$((fails + 1)); }
# run ARGS...: the script on the planted tree; sets out, err and rc.
run() {
  rc=0
  out=$(GITHUB_REPOSITORY=acme/genome bash "$NOTES" --repo "$T" "$@" 2> "$E") || rc=$?
  err=$(cat "$E")
}
# refused WHAT CAUSE ARGS...: the script fails, names CAUSE and prints no notes.
refused() {
  local what=$1 cause=$2
  shift 2
  run "$@"
  if [ "$rc" != 0 ] && [ -z "$out" ] && grep -qF -- "$cause" <<< "$err"; then
    ok "$what"
  else
    bad "$what (exit ${rc}, stdout '${out}', stderr '${err}')"
  fi
}
# has WHAT TEXT: the last notes hold TEXT on a line of its own.
has() { if grep -qxF -- "$2" <<< "$out"; then ok "$1"; else bad "$1: no line '$2' in: ${out}"; fi; }
lacks() { if grep -qF -- "$2" <<< "$out"; then bad "$1: '$2' is in: ${out}"; else ok "$1"; fi; }
g() { git -C "$T" -c user.name=t -c user.email=t@example.invalid -c tag.gpgSign=false -c commit.gpgSign=false -c core.hooksPath=/dev/null "$@"; }

g init -q
cat > "${T}/CHANGELOG.md" <<'C'
# Changelog

## Unreleased

- Not in any release yet.

## 1.2.30 (2026-12-09)

- Only in 1.2.30.

## 1.2.3 (2026-12-01)

- Fixes the [step 6 page](docs/06-clinvar-screen.md#floor) and keeps [the site](https://example.org/x) and [this anchor](#top).

## 1.2.2

## 1.2.1

- Only in 1.2.1.

## 2.0.0 (2026-12-20)

- The second major.

## Releases

- 0.2.0 to 0.8.2: https://github.com/acme/genome/releases
C
g add -A && g commit -qm "docs: the changelog"
for t in v1.2.3 v1.2.2 v1.2.4 v2.0.0; do g tag "$t"; done

# A patch release: its own entry, nothing from the entries around it.
run v1.2.3
if [ "$rc" = 0 ]; then ok "v1.2.3 has notes"; else bad "v1.2.3 has notes (exit ${rc}): ${err}"; fi
has "the notes are the entry's line, links into the tagged tree" \
  "- Fixes the [step 6 page](https://github.com/acme/genome/blob/v1.2.3/docs/06-clinvar-screen.md#floor) and keeps [the site](https://example.org/x) and [this anchor](#top)."
lacks "the entry ends at the next heading" "Only in 1.2.1"
lacks "## 1.2.30 is not the entry of 1.2.3" "Only in 1.2.30"
lacks "the heading itself is not in the notes" "## 1.2.3"
lacks "a minor or patch release has no upgrade link" "Upgrading"

refused "a tag with no CHANGELOG entry fails" "CHANGELOG.md has no '## 1.2.4' entry" v1.2.4
refused "an empty entry fails" "CHANGELOG.md has no '## 1.2.2' entry" v1.2.2
refused "a tag that does not exist fails" "v9.9.9" v9.9.9
refused "a tag not of the form vMAJOR.MINOR.PATCH fails" "vMAJOR.MINOR.PATCH" 1.2.3

# A major release with no UPGRADING.md, then with only another major's section.
refused "a major with no UPGRADING.md fails" "UPGRADING.md has no '## Upgrading to 2.0' section" v2.0.0
printf '# Upgrading\n\n## Upgrading to 3.0\n\nSteps.\n' > "${T}/UPGRADING.md"
g add UPGRADING.md && g commit -qm "docs: a 3.0 section only"
refused "a major whose UPGRADING.md has no '## Upgrading to 2.0' fails" "UPGRADING.md has no '## Upgrading to 2.0' section" v2.0.0 HEAD

# The same major once UPGRADING.md has the section: read through REF before the tag moves.
printf '# Upgrading\n\n## Upgrading to 2.0\n\nWhat changed, then the steps.\n' > "${T}/UPGRADING.md"
g commit -qam "docs: the 2.0 upgrade section"
refused "the tag still reads the tagged commit, not the newer one" "UPGRADING.md has no '## Upgrading to 2.0' section" v2.0.0
run v2.0.0 HEAD
if [ "$rc" = 0 ]; then ok "v2.0.0 at HEAD has notes"; else bad "v2.0.0 at HEAD has notes (exit ${rc}): ${err}"; fi
first=$(sed -n '/[^[:space:]]/{p;q;}' <<< "$out")
want="Upgrading from 1.x? Read [Upgrading to 2.0](https://github.com/acme/genome/blob/v2.0.0/UPGRADING.md#upgrading-to-20): what changed and every step."
if [ "$first" = "$want" ]; then ok "a major's notes open with the link to its UPGRADING.md section"; else bad "first line of the 2.0.0 notes: '${first}'"; fi
has "then the entry" "- The second major."

# The real page: the 1.0.0 cut needs this section, and its anchor is the one the notes link.
if grep -qxF "## Upgrading to 1.0" "${REPO}/UPGRADING.md"; then ok "UPGRADING.md has '## Upgrading to 1.0'"; else bad "UPGRADING.md has no '## Upgrading to 1.0' line"; fi
cp "${REPO}/UPGRADING.md" "${T}/UPGRADING.md"
printf '# Changelog\n\n## 1.0.0 (2026-12-01)\n\n- The first stable release.\n' > "${T}/CHANGELOG.md"
g commit -qam "docs: the real upgrade page"
run v1.0.0 HEAD
if [ "$rc" = 0 ] && grep -qF "UPGRADING.md#upgrading-to-10)" <<< "$out"; then ok "a 1.0.0 entry with the real UPGRADING.md gives notes that link its 1.0 section"; else bad "1.0.0 with the real UPGRADING.md (exit ${rc}): ${out} ${err}"; fi

if [ "$fails" -gt 0 ]; then echo "${fails} check(s) failed"; exit 1; fi
echo "All release-notes checks passed"
