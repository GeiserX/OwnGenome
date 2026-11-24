#!/usr/bin/env bash
# og recheck and step 06's conflicting subset, through bin/write_manifest.sh,
# for one ClinVar release that never changes:
#   1. a run with no subset installed records clinvar_conflicting "not
#      installed", and og recheck finds nothing to redo;
#   2. the subset installed later (setup.sh after an upgrade) makes step 06 and
#      the report rerun-ready, and nothing downstream of step 06;
#   3. a manifest written with the subset installed records the value og
#      recheck reads (the subset's ##fileDate and sha256): nothing to redo;
#   4. the subset rebuilt with other records: rerun-ready again;
#   5. a manifest from before the row existed, with the subset here, is
#      rerun-ready (the run had no subset).
# The upstream answers come from an offline file: nothing here asks the network.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
S="${G}/sample1"
M="${S}/run_manifest.tsv"
STATE="${S}/recheck/state.json"
CONF="${G}/clinvar/clinvar_conflicting_plp_chr.vcf.gz"
seed_reference "$G"
mkdir -p "$S" "${G}/clinvar"
printf '##fileformat=VCFv4.1\n##fileDate=20260901\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
  | gzip -c > "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz"
ANSWERS="${REPO_ROOT}/tests/fixtures/recheck/offline-nothing-newer.json"

# subset ID: the conflicting subset of the same release, one record, with its .tbi.
subset() {
  printf '##fileformat=VCFv4.1\n##fileDate=20260901\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t%s\tA\tG\t.\t.\tCLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Pathogenic_(2)|Uncertain_significance_(1)\n' "$1" \
    | gzip -c > "$CONF"
  printf 'tbi\n' > "${CONF}.tbi"
}
manifest() { run_expect 0 "$1" bash "${REPO_ROOT}/bin/write_manifest.sh" sample1 test; }
row() { awk -F'\t' '$1 == "data" && $2 == "clinvar_conflicting" {print $3}' "$M"; }
# recheck NAME: og recheck on the sample; its exit code is not this case's
# subject (an item of another database may have no answer offline).
recheck() {
  rm -f "$STATE"
  "${REPO_ROOT}/bin/og" recheck --sample sample1 --offline "$ANSWERS" > "${CASE_WORK}/$1.out" 2>&1 || true
  [ -s "$STATE" ] || fail "$1: og recheck wrote no recheck/state.json: $(cat "${CASE_WORK}/$1.out")"
}
item() { jq -r '.items[] | select(.item == "clinvar_conflicting") | .'"$1"' // "null"' "$STATE"; }
ready() { jq -r '.rerun_ready_steps | join(" ")' "$STATE"; }
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1"; else shasum -a 256 "$1"; fi | awk '{print $1}'; }

# 1. No subset
manifest none
[ "$(row)" = "not installed" ] || fail "without the subset the manifest row reads '$(row)'"
recheck none
[ "$(item class)" = no-newer-release ] || fail "no subset, none at the run: class $(item class)"
[ -z "$(ready)" ] || fail "no subset, none at the run: rerun-ready steps $(ready)"

# 2. The subset appears for the same release
subset 1
recheck appeared
[ "$(item class)" = rerun-ready ] || fail "the subset appeared: class $(item class)"
[ "$(ready)" = "06 24" ] || fail "the subset appeared: rerun-ready steps '$(ready)', expected '06 24'"
[ "$(item installed)" = "20260901 (sha256 $(sha "$CONF"))" ] || fail "og recheck reads the subset as '$(item installed)'"

# 3. A run that had it
manifest recorded
[ "$(row)" = "20260901 (sha256 $(sha "$CONF"))" ] || fail "with the subset the manifest row reads '$(row)'"
recheck recorded
[ "$(item ran_with)" = "$(item installed)" ] || fail "the manifest row ($(item ran_with)) is not what og recheck reads ($(item installed))"
[ "$(item class)" = no-newer-release ] || fail "the run's own subset: class $(item class)"
[ -z "$(ready)" ] || fail "the run's own subset: rerun-ready steps $(ready)"

# 4. Rebuilt with other records, same release
subset 2
recheck rebuilt
[ "$(item class)" = rerun-ready ] || fail "the subset rebuilt: class $(item class)"
[ "$(ready)" = "06 24" ] || fail "the subset rebuilt: rerun-ready steps '$(ready)'"

# 5. A manifest from before the row existed
grep -v $'^data\tclinvar_conflicting\t' "$M" > "${M}.old" && mv "${M}.old" "$M"
recheck old-manifest
[ "$(item ran_with)" = null ] || fail "a manifest without the row: ran_with $(item ran_with)"
[ "$(item class)" = rerun-ready ] || fail "a manifest without the row, subset here: class $(item class)"
echo "OK: og recheck marks step 06 when the conflicting subset appears or is rebuilt for the same ClinVar release"
