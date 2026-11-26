#!/usr/bin/env bash
# run-all.sh never leaves an earlier run's OwnHealth file (facts/<sample>_genome_facts.json)
# beside a new run's results, where it could be imported by mistake:
#   - step 24 fails: the file is skipped (summary.json may be an earlier
#     run's) and the earlier file is removed;
#   - step 24 works but og facts writes nothing (the fake og here): the
#     earlier file is removed before og facts runs, so none is left;
#   - the earlier file is gone before the pipeline starts, so a run cut
#     short before og facts leaves none beside the new results.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
mkdir -p "${G}/prs_scores"

# The ExpansionHunter catalog comes out of its image; step 24's render fails while fail-24 exists.
use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
case "${*:2}" in
  "cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json") echo '[{"LocusId": "HTT"}]' ;;
  "/pgp-bin/og report render"*)
    [ ! -e "${GENOME_DIR}/sample1/facts/sample1_genome_facts.json" ] || touch "${CASE_WORK}/facts-seen-at-24"
    [ ! -e "${CASE_WORK}/fail-24" ] || { touch "${CASE_WORK}/render-failed"; exit 1; } ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

FACTS="${G}/sample1/facts/sample1_genome_facts.json"
stale() { mkdir -p "${G}/sample1/facts"; echo '{"format":"an earlier run"}' > "$FACTS"; }

stale; touch "${CASE_WORK}/fail-24"
run_expect 1 report-failed "${SCRIPTS}/run-all.sh" sample1 male
[ -e "${CASE_WORK}/render-failed" ] || fail "step 24 failed before its render, not at the planted failure"
output_has report-failed '^\[ownhealth-facts\] skipped: 24-html-report.sh failed'
[ ! -e "$FACTS" ] || fail "step 24 failed and the earlier run's OwnHealth file is still there"
[ ! -e "${CASE_WORK}/facts-seen-at-24" ] || fail "the earlier run's OwnHealth file was still there when step 24 ran"
! grep -q '/pgp-bin/og facts' "$FAKE_DOCKER_LOG" || fail "og facts ran although step 24 failed"

stale; rm -f "${CASE_WORK}/fail-24"
run_expect 0 report-ok "${SCRIPTS}/run-all.sh" sample1 male
grep -q '/pgp-bin/og facts --sample sample1 ' "$FAKE_DOCKER_LOG" || fail "og facts did not run after step 24 worked"
[ ! -e "$FACTS" ] || fail "og facts wrote nothing and the earlier run's OwnHealth file is still there"
output_has report-ok 'the OwnHealth file [^ ]*/sample1/facts/sample1_genome_facts\.json'
echo "run-all.sh removes an earlier run's OwnHealth file before the pipeline starts, and writes none when step 24 fails."
