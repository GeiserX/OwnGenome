#!/usr/bin/env bash
# scripts/ci/check-release.sh on planted trees, then on the real one.
#   tag:      a tag that differs from manifest.version fails (the state main
#             was in before the release check: 0.8.2 declared, v1.0.0 pushed);
#             a matching tag passes, lightweight or annotated, major, minor or
#             patch, with no notes page (the notes are the CHANGELOG.md entry,
#             tests/test_release_notes.sh); a missing tag fails.
#   checks:   a run still going, a failed or cancelled run, no ci-ok, or no
#             runs at all fail; the ignored release job itself does not count.
#   declared: a declared X.Y.0 needs no notes page; a manifest with no
#             version, or one not of the form MAJOR.MINOR.PATCH, fails.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
CHECK="${REPO}/scripts/ci/check-release.sh"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fails=0
ok()  { echo "PASS $*"; }
bad() { echo "FAIL $*"; fails=$((fails + 1)); }
# expect RC WHAT ARGS...: run the check, compare its exit code (0 or non-0).
expect() {
  local want=$1 what=$2 rc=0
  shift 2
  out=$(bash "$CHECK" "$@" 2>&1) || rc=$?
  if { [ "$want" = 0 ] && [ "$rc" = 0 ]; } || { [ "$want" != 0 ] && [ "$rc" != 0 ]; }; then
    ok "$what"
  else
    bad "$what (exit ${rc}): ${out}"
  fi
}

# cause TEXT: the last refusal names its cause, so an unrelated failure cannot pass.
cause() { grep -qF -- "$1" <<< "$out" || bad "the refusal names its cause ('$1'): ${out}"; }
g() { git -C "$T" -c user.name=t -c user.email=t@example.invalid -c tag.gpgSign=false -c commit.gpgSign=false -c core.hooksPath=/dev/null "$@"; }
config() {
  cat > "${T}/nextflow.config" <<NF
params {
    cyrius_install = null   // tools/cyrius-<version>
    version        = 'not this one'
}
manifest {
    name            = 'OwnGenome'
    nextflowVersion = '!>=25.10.4'
    version         = '$1'
}
NF
}
g init -q
config 0.8.2
g add -A && g commit -qm "feat: the last commit before the cut"
g tag v0.8.2
g commit -q --allow-empty -m "fix: a change after 0.8.2"
g tag v1.0.0

# --- tag ---
expect 1 "tag v1.0.0 on a commit that declares 0.8.2 is refused" --repo "$T" tag v1.0.0
grep -q "does not match the version" <<< "$out" || bad "the mismatch message names the cause: ${out}"
expect 1 "a tag that does not exist is refused" --repo "$T" tag v9.9.9
cause "does not exist"
expect 1 "a tag not of the form vX.Y.Z is refused" --repo "$T" tag 1.0.0
cause "not of the form"

config 1.0.0
g add -A && g commit -qm "chore: declare 1.0.0"
g tag -f v1.0.0 >/dev/null
expect 0 "v1.0.0 with matching manifest.version passes, with no notes page" --repo "$T" tag v1.0.0
g tag -f -a -m "OwnGenome 1.0.0" v1.0.0 >/dev/null
expect 0 "an annotated v1.0.0 tag with matching manifest.version passes" --repo "$T" tag v1.0.0
g tag -a -m "OwnGenome 1.1.0" v1.1.0
expect 1 "an annotated v1.1.0 tag on a commit that declares 1.0.0 is refused" --repo "$T" tag v1.1.0
cause "does not match the version"
g tag -d v1.1.0 >/dev/null
expect 1 "the removed --notes-out option is refused, not ignored" --repo "$T" tag v1.0.0 --notes-out "${T}/body.md"
cause "tag takes one argument"

config 1.0.1
g add -A && g commit -qm "fix: a patch after 1.0.0"
g tag v1.0.1
expect 0 "patch tag v1.0.1 with matching manifest.version passes" --repo "$T" tag v1.0.1

# --- checks ---
tab=$'\t'
runs() { printf '%s\n' "$@" > "${T}/runs.tsv"; }
runs "ci-ok${tab}completed${tab}success" "shellcheck${tab}completed${tab}success" "deploy${tab}completed${tab}skipped" \
     "create-release${tab}in_progress${tab}"
expect 0 "all checks green, the release job itself left out" checks "${T}/runs.tsv" --ignore create-release
expect 1 "the release job counts when it is not left out" checks "${T}/runs.tsv"
cause "NOT FINISHED: create-release"
runs "ci-ok${tab}completed${tab}success" "nextflow (pinned)${tab}in_progress${tab}"
expect 1 "a check still running is refused" checks "${T}/runs.tsv"
cause "NOT FINISHED"
runs "ci-ok${tab}completed${tab}success" "unit-tests${tab}completed${tab}failure"
expect 1 "a failed check is refused" checks "${T}/runs.tsv"
cause "NOT GREEN: unit-tests"
runs "ci-ok${tab}completed${tab}success" "e2e${tab}completed${tab}cancelled"
expect 1 "a cancelled check is refused" checks "${T}/runs.tsv"
cause "NOT GREEN: e2e"
runs "shellcheck${tab}completed${tab}success"
expect 1 "no ci-ok among the checks is refused" checks "${T}/runs.tsv"
cause "ci-ok did not succeed"
runs "ci-ok${tab}completed${tab}failure"
expect 1 "a failed ci-ok is refused" checks "${T}/runs.tsv"
cause "ci-ok did not succeed"
: > "${T}/runs.tsv"
expect 1 "a commit with no check runs is refused" checks "${T}/runs.tsv"
cause "no check runs"

# --- declared ---
D="${T}/declared"
mkdir -p "$D"
cat > "${D}/nextflow.config" <<'NF'
manifest {
    version = '1.0.0'
}
NF
expect 0 "declared 1.0.0 with no notes page passes" --repo "$D" declared
sed -i.bak "s/'1.0.0'/'1.0'/" "${D}/nextflow.config"
expect 1 "a declared version not of the form MAJOR.MINOR.PATCH is refused" --repo "$D" declared
cause "which is not MAJOR.MINOR.PATCH"
printf 'manifest {\n    name = "x"\n}\n' > "${D}/nextflow.config"
expect 1 "a manifest with no version is refused" --repo "$D" declared
cause "declares no version"

# --- the real tree ---
expect 0 "the real tree: the declared version" declared
echo "$out"

[ "$fails" -eq 0 ] || { echo "${fails} release check case(s) failed"; exit 1; }
echo "All release check cases passed"
