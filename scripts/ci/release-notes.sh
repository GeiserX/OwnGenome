#!/usr/bin/env bash
# release-notes.sh: the GitHub Release notes of a tag, from its CHANGELOG.md
# entry. The release workflow runs it on every v* tag; by hand it prints the
# notes of any tag, or of a commit before it is tagged.
#
#   release-notes.sh [--repo DIR] TAG [REF]
#
# The notes are the body of the "## X.Y.Z" entry of CHANGELOG.md (a space and
# more may follow the version, such as a date), up to the next "## " heading,
# as REF holds the file. REF defaults to TAG, so a tag's notes are what the
# tagged commit says. A major release (X.0.0) must have a "## Upgrading to X.0"
# section in UPGRADING.md at REF, and its notes open with a link to it.
# Relative links point into the tagged tree on GitHub, so they work on the
# Releases page. A missing or empty entry, a major without its section, or
# notes too long for a Release body fail with nothing on stdout, before
# anything is published.
#
# GitHub refuses a Release body over 125,000 characters. The release job adds
# its og binaries paragraph after these notes, so the notes may hold at most
# NOTES_MAX characters (Unicode characters, not bytes): 125,000 less a
# 2,000-character reserve for that paragraph.
#
# --repo DIR reads another git work tree (the tests plant one); the default is
# the checkout this script is in. GITHUB_REPOSITORY names the repository the
# links point to (default GeiserX/OwnGenome).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
if [ "${1:-}" = --repo ]; then ROOT=${2:?--repo needs a directory}; shift 2; fi
TAG=${1:?usage: release-notes.sh [--repo DIR] vX.Y.Z [REF]}
REF=${2:-$TAG}
REPO=${GITHUB_REPOSITORY:-GeiserX/OwnGenome}
BODY_MAX=125000 NOTES_MAX=123000

die() { echo "release-notes: $*" >&2; exit 1; }
[[ "$TAG" =~ ^v([0-9]+)\.([0-9]+)\.([0-9]+)$ ]] || die "the tag is not of the form vMAJOR.MINOR.PATCH"
major=${BASH_REMATCH[1]} minor=${BASH_REMATCH[2]} patch=${BASH_REMATCH[3]}
version=${TAG#v}
git -C "$ROOT" rev-parse -q --verify "${REF}^{commit}" > /dev/null || die "${REF} is not a commit or tag in this repository"

changelog=$(git -C "$ROOT" show "${REF}:CHANGELOG.md" 2> /dev/null) || die "${REF} has no CHANGELOG.md"
# The entry's lines, without the blank lines before its first line.
notes=$(awk -v h="## ${version}" '
  found && /^## / { exit }
  found { print }
  $0 == h || index($0, h " ") == 1 { found = 1 }
' <<< "$changelog" | sed '/[^[:space:]]/,$!d')
[ -n "$(tr -d '[:space:]' <<< "$notes")" ] || die "CHANGELOG.md has no '## ${version}' entry at ${REF}, or it is empty: write the release's changes there first"

if [ "$minor" = 0 ] && [ "$patch" = 0 ]; then
  section="## Upgrading to ${major}.0"
  git -C "$ROOT" show "${REF}:UPGRADING.md" 2> /dev/null | grep -qxF "$section" \
    || die "${version} is a major release and UPGRADING.md has no '${section}' section at ${REF}: write what changed and every step to adapt there first"
  notes="Upgrading from $((major - 1)).x? Read [Upgrading to ${major}.0](UPGRADING.md#upgrading-to-${major}0): what changed and every step.

${notes}"
fi

# A link target with no scheme, no anchor-only form and no mailto: is a path
# in the repository: point it at the tagged tree.
out=$(printf '%s\n' "$notes" | TAG="$TAG" REPO="$REPO" perl -pe \
  's{\]\((?!https?://|#|mailto:)([^)\s]+)\)}{](https://github.com/$ENV{REPO}/blob/$ENV{TAG}/$1)}g')

# The size of what this script prints, newline included, in characters.
size=$(printf '%s\n' "$out" | perl -CSD -0777 -ne 'print length')
[ "$size" -le "$NOTES_MAX" ] \
  || die "the notes of ${TAG} are ${size} characters, over the limit of ${NOTES_MAX} (a GitHub Release body holds at most ${BODY_MAX}, and the release job adds up to 2,000 after the notes): cut the CHANGELOG.md entry down first"
printf '%s\n' "$out"
