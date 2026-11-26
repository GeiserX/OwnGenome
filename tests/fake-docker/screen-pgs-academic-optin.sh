#!/usr/bin/env bash
# The four PGS Catalog scores free for academic research use only (licence
# academic-research in assets/pgs_scores.tsv) are opt-in, against the fake
# docker, wget and nextflow:
#   - a plain setup.sh downloads every other score of the list and none of
#     the four; setup.sh --prs-academic downloads exactly the four;
#   - step 25 by default neither downloads them nor passes --academic to
#     og prs format; PRS_ACADEMIC=true downloads them and passes it;
#   - PRS_ACADEMIC with another value than true or false is refused.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

LIST="${REPO_ROOT}/assets/pgs_scores.tsv"
# ids LICENCE: the ids of the list whose licence column is LICENCE, sorted.
ids() {
  awk -F'\t' -v want="$1" '!h && $1 == "pgs_id" { for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next }
    h && $1 ~ /^PGS[0-9]+$/ && $c["licence"] == want {print $1}' "$LIST" | LC_ALL=C sort
}
# fetched: the PGS ids whose scoring file the fake curl or wget was asked for, sorted.
fetched() { grep -oE '^(curl|wget) [^ ]*/PGS[0-9]+_hmPOS_GRCh38\.txt\.gz ' "$FAKE_DOCKER_LOG" | grep -oE 'PGS[0-9]+_' | tr -d _ | LC_ALL=C sort -u; }
ACADEMIC=$(ids academic-research)
AUTHOR=$(ids author)
[ "$(tr '\n' ' ' <<<"$ACADEMIC")" = "PGS000014 PGS000016 PGS000017 PGS000027 " ] || fail "the list does not mark the four academic-only scores: ${ACADEMIC}"
[ "$(grep -c . <<<"$AUTHOR")" -eq 11 ] || fail "the list has not 11 author-terms scores: ${AUTHOR}"

# --- setup.sh ---------------------------------------------------------------------------
seed_clinical_downloads
G="${CASE_WORK}/genome"
mkdir -p "$G"
: > "$FAKE_DOCKER_LOG"
run_expect 3 setup "${SCRIPTS}/setup.sh" "$G"
[ "$(fetched)" = "$AUTHOR" ] || fail "a plain setup did not download exactly the 11 author-terms scores:
$(diff <(echo "$AUTHOR") <(fetched))"
for id in $ACADEMIC; do [ ! -e "${G}/prs_scores/${id}.txt.gz" ] || fail "a plain setup installed ${id}"; done

: > "$FAKE_DOCKER_LOG"
run_expect 0 prs-academic "${SCRIPTS}/setup.sh" --prs-academic "$G"
[ "$(fetched)" = "$ACADEMIC" ] || fail "setup.sh --prs-academic did not download exactly the four:
$(diff <(echo "$ACADEMIC") <(fetched))"
for id in $ACADEMIC; do [ -s "${G}/prs_scores/${id}.txt.gz" ] || fail "setup.sh --prs-academic did not install ${id}"; done
output_has prs-academic '^  PGS000014: freely available to the academic community for research use\.$'
output_has prs-academic '^\[OK\] Academic-only PGS Catalog scores'

# --- step 25 ----------------------------------------------------------------------------
G="${CASE_WORK}/run"
mkdir -p "$G"
export GENOME_DIR="$G" NXF_HOME="${CASE_WORK}/nxf"
seed_reference "$G"
seed_sample "$G" sample1
mkdir -p "${G}/tools/pgsc_calc-${PGSC_CALC_VERSION}" "${NXF_HOME}/plugins/nf-schema-${PGSC_CALC_NF_SCHEMA}" "${G}/prs_scores"
: > "${G}/tools/pgsc_calc-${PGSC_CALC_VERSION}/main.nf"
BODY='#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\n20\t100\tA\t1\n'
for id in $AUTHOR; do printf '%b' "$BODY" | gzip -c > "${G}/prs_scores/${id}.txt.gz"; done
# What the fake wget serves for the four, so their downloads pass the build check.
export FAKE_DOWNLOAD_DIR="${CASE_WORK}/downloads" FETCH_TRIES=1 FETCH_WAIT=0
mkdir -p "$FAKE_DOWNLOAD_DIR"
for id in $ACADEMIC; do printf '%b' "$BODY" | gzip -nc > "${FAKE_DOWNLOAD_DIR}/${id}_hmPOS_GRCh38.txt.gz"; done

use_output_hook
cat > "${CASE_WORK}/hook-prs" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
args=("$@")
val() { local i; for ((i = 0; i < ${#args[@]} - 1; i++)); do [ "${args[i]}" = "$1" ] && { host_path "${args[i + 1]}"; return; }; done; return 1; }
case " $* " in
  *" prs format "*)
    out=$(val --out); alleles=$(val --alleles)
    mkdir -p "$out"; : > "${out}/PGS000018.txt.gz"
    printf 'chr20\t100\tA\nchr20\t100\tG\n' > "$alleles" ;;
  *" prs table "*)
    out=$(val --out)
    printf 'Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\tMatched_Pct\tPercentile\tAncestry_Group\tChrX\tAmbiguous\tApplies_To\tSex_Check\tInput\n' > "$out"
    printf 'Coronary artery disease\tPGS000018\t1\t1\t1\t100.0\tNA\tNA\t0\t0\t\t\tgvcf\n' >> "$out" ;;
  *" bcftools "*)
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\nchr20\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n' ;;
esac
exit 0
HOOK
chmod +x "${CASE_WORK}/hook-prs"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-prs"

: > "$FAKE_DOCKER_LOG"
run_expect 0 default "${SCRIPTS}/25-prs.sh" sample1 female
[ -z "$(fetched)" ] || fail "step 25 downloaded scores by default: $(fetched | tr '\n' ' ')"
for id in $ACADEMIC; do [ ! -e "${G}/prs_scores/${id}.txt.gz" ] || fail "step 25 installed ${id} by default"; done
docker_log_has 'og prs format ' "step 25 did not run og prs format"
if grep -E 'og prs format .*--academic' "$FAKE_DOCKER_LOG"; then fail "og prs format got --academic without PRS_ACADEMIC"; fi
output_has default 'Academic-only scores: left out \(PRS_ACADEMIC=true to include\)'

: > "$FAKE_DOCKER_LOG"
PRS_ACADEMIC=true run_expect 0 academic "${SCRIPTS}/25-prs.sh" sample1 female
[ "$(fetched)" = "$ACADEMIC" ] || fail "PRS_ACADEMIC=true did not download exactly the four:
$(diff <(echo "$ACADEMIC") <(fetched))"
docker_log_has 'og prs format .*--academic' "PRS_ACADEMIC=true did not reach og prs format"
output_has academic 'Academic-only scores: scored \(PRS_ACADEMIC=true\)'

: > "$FAKE_DOCKER_LOG"
PRS_ACADEMIC=yes run_expect 1 bad-value "${SCRIPTS}/25-prs.sh" sample1
output_has bad-value "PRS_ACADEMIC must be 'true' or 'false', got 'yes'"
if grep -qE '^run ' "$FAKE_DOCKER_LOG"; then fail "the refused run started a container"; fi
echo "the four academic-only scores are installed and scored only on request"
