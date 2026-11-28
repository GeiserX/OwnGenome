#!/usr/bin/env bash
# The live site (og live serve) on a real Nextflow run: the pinned Nextflow
# loads the nf-weblog plugin from the config run-all.sh writes when
# OWNGENOME_LIVE_URL is set (the same text, the pin read from versions.env), posts each
# task to the site's ingest port, and the site stores the messages and shows
# step 11 (ROH) done in state.json. A control reads state.json before the run:
# step 11 is 'not planned' then, so the state check can fail. The token goes as
# the weblog's HTTP Basic password, which Nextflow reads from
# OWNGENOME_LIVE_TOKEN: the plugin names its URL in .nextflow.log on every
# post (the check that the log has those lines is the control), and neither
# that log, a rotated one nor live.config holds the token.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
PLUGIN=$(sed -n 's/^NF_WEBLOG_PLUGIN="\([^"]*\)"$/\1/p' "${REPO}/versions.env")
check "versions.env pins nf-weblog with a version" has '^nf-weblog@[0-9]+\.[0-9]+\.[0-9]+$' "$PLUGIN"

G="$GENOME_DIR"
ROOT="${CASE_TMP}/root" STATE="${CASE_TMP}/state"
mkdir -p "${ROOT}/${SAMPLE}/nextflow"   # a run folder the site lists (it has nextflow/)
"${REPO}/bin/og" live serve --root "$ROOT" --state-dir "$STATE" --port 0 --ingest-port 0 --quiet \
  2> "${CASE_TMP}/site.err" &
SITE=$!
trap 'kill "$SITE" 2>/dev/null' EXIT
# Port 0: the site's start line, printed once both listen, names the ports it bound. Each port
# must be followed by '/', so a line still being written is never read.
WEB="" ING=""
for _ in $(seq 100); do
  read -r WEB ING < <(sed -n 's|^og live serve: browser http://127\.0\.0\.1:\([0-9]*\)/ ; ingest http://127\.0\.0\.1:\([0-9]*\)/.*|\1 \2|p' \
    "${CASE_TMP}/site.err") && break
  kill -0 "$SITE" 2>/dev/null || break
  sleep 0.1
done
check "the site printed the ports it bound" has '^[0-9]+ [0-9]+$' "${WEB} ${ING}"
# og live serve takes a message only with its token (docs/live-site.md).
TOK=$(cat "${STATE}/ingest-token")
check "the site made a token of 43 characters" has '^[A-Za-z0-9_-]{43}$' "$TOK"
state11() { curl -s "http://127.0.0.1:${WEB}/runs/${SAMPLE}/state.json" \
  | jq -er 'first(.steps[] | select(.step == "11") | .state // empty)' 2>/dev/null \
  || echo "no state.json"; }
check_eq "control: before the run, step 11 is not planned" "$(state11)" "not planned"

SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$SAMPLE" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
  "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" > "$SHEET"
(umask 077; printf "plugins {\n    id '%s'\n}\nweblog {\n    enabled = true\n    url     = '%s'\n%s}\n" "$PLUGIN" \
  "http://127.0.0.1:${ING}/ingest/${SAMPLE}" "    basicToken = \"og:\${env('OWNGENOME_LIVE_TOKEN')}\""$'\n' > "${CASE_TMP}/live.config")
echo "+ nextflow run main.nf --tools roh -c live.config"
(
  cd "$CASE_TMP" && OWNGENOME_LIVE_TOKEN="$TOK" nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -c "${CASE_TMP}/live.config" \
    -work-dir "${CASE_TMP}/work" --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools roh --outdir "${CASE_TMP}/out" --max_cpus 2 --max_memory 6.GB
) 2>&1 | tee "$STEP_LOG"
check_eq "nextflow run exits 0" "${PIPESTATUS[0]}" 0
check "the run loaded the plugin" grep -q "nf-weblog" "${CASE_TMP}/.nextflow.log"
# Control: the plugin logs the URL of each post it sent, so the scan below reads lines that name it.
check "control: .nextflow.log names the URL of a post the site took" \
  grep -qF "Successfully sent message to http://127.0.0.1:${ING}/ingest/${SAMPLE}" "${CASE_TMP}/.nextflow.log"
for f in "${CASE_TMP}"/.nextflow.log* "${CASE_TMP}/live.config"; do
  check "$(basename "$f") does not hold the ingest token" lacks "$TOK" "$(cat "$f")"
done

EV="${STATE}/runs/${SAMPLE}/events.jsonl"
# count EVENT PROCESS_SUFFIX: the stored messages of that event (and process, unless "").
# jq prints 0 before it fails on a missing file, so its output is kept only when it exits 0.
count() {
  local n
  n=$(jq -n --arg e "$1" --arg p "$2" \
    '[inputs | select(.event == $e and ($p == "" or (((.trace // {}).process // "") | endswith($p))))] | length' "$EV" 2>/dev/null) \
    || n=0
  echo "$n"
}
check_ge "the site stored the run's 'started' message" "$(count started '')" 1
check_ge "the site stored ROH's process_completed" "$(count process_completed ROH)" 1
check_ge "the site stored the run's 'completed' message" "$(count completed '')" 1
check "no stored message holds the samplesheet or the VCF path" lacks "${SHEET}|${SAMPLE}\\.vcf\\.gz" "$(cat "$EV" 2>/dev/null)"
check_eq "state.json shows step 11 done" "$(state11)" "done"
finish
