#!/usr/bin/env bash
# run-all.sh and the live site (docs/live-site.md):
#   1. with OWNGENOME_LIVE_URL set to a running og live serve's ingest port, run-all.sh
#      reads the token of the site's state folder (OWNGENOME_LIVE_STATE_DIR) into
#      OWNGENOME_LIVE_TOKEN; nextflow starts with OWNGENOME_LIVE_RUN=sample1, that token
#      and -c <sample>/nextflow/live.config (mode 600), which loads the pinned nf-weblog
#      plugin, points the weblog at OWNGENOME_LIVE_URL/ingest/sample1 and names the token
#      only as env('OWNGENOME_LIVE_TOKEN'), its HTTP Basic password;
#      the script step 04b (GRIDSS=true) posts 'started' and its result, and
#      the run ends with a pgp_run_finished message carrying the exit code, each with
#      the token in a curl config on fd 3 and never in an argument; the site
#      takes every one of those posts, and refuses a post without the token, with
#      another token and a message it does not know; no line of the output, no
#      argument and neither live.config nor `meta command` hold the token;
#      with OWNGENOME_LIVE_PAGE_URL set it prints the run's page, which the site serves;
#   2. without OWNGENOME_LIVE_URL, nextflow sees no OWNGENOME_LIVE_RUN and no live.config,
#      and nothing is posted; no page is printed, even with OWNGENOME_LIVE_PAGE_URL;
#      a URL that is not a plain http://host:port, with the token or without it,
#      exits 2;
#   3. with the token in the URL (http://host:port/<token>) and the site refusing
#      every post, the exit code is the same as in 1 and one warning says so,
#      naming the site without its token; the posts still carry it on fd 3, and
#      `meta command`, which og recheck --rerun prints and replays, names the
#      plain URL only.
# Wrappers in front of the fake nextflow and the fake curl log what they see.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
use_output_hook
seed_sample "$G" sample1

mkdir -p "${CASE_WORK}/bin"
cat > "${CASE_WORK}/bin/nextflow" <<'NF'
#!/usr/bin/env bash
echo "nextflow OWNGENOME_LIVE_RUN=${OWNGENOME_LIVE_RUN:-unset} OWNGENOME_LIVE_URL=${OWNGENOME_LIVE_URL:-unset} OWNGENOME_LIVE_TOKEN=${OWNGENOME_LIVE_TOKEN:-unset}" >> "${CASE_WORK}/live.log"
exec "${REPO_ROOT}/scripts/ci/fake-docker/nextflow" "$@"
NF
cat > "${CASE_WORK}/bin/curl" <<'CURL'
#!/usr/bin/env bash
echo "argv $*" >> "${CASE_WORK}/argv.log"
for a in "$@"; do case "$a" in */ingest/*)
  echo "post ${a} $(cat)" >> "${CASE_WORK}/live.log"
  for k in "$@"; do [ "$k" != -K/dev/fd/3 ] || echo "auth ${a} $(cat /dev/fd/3)" >> "${CASE_WORK}/live.log"; done ;;
esac; done
exec "${REPO_ROOT}/scripts/ci/fake-docker/curl" "$@"
CURL
chmod +x "${CASE_WORK}/bin/nextflow" "${CASE_WORK}/bin/curl"
P="${CASE_WORK}/bin:${PATH}"
LOG="${CASE_WORK}/live.log"
# A real og live serve: its ingest port, and the token run-all.sh reads from its state folder.
OG="${REPO_ROOT}/bin/og"
"$OG" version > /dev/null || fail "og does not run here: build it into bin/.og (rust/og.sh, the Guard fake-docker job does)"
CURL=$(PATH=/usr/local/bin:/usr/bin:/bin command -v curl) || fail "the real curl is needed"
STATE="${CASE_WORK}/live-state"
WEB=$("$OG" dev fixture free-port) ING=$("$OG" dev fixture free-port)
"$OG" live serve --root "$G" --state-dir "$STATE" --port "$WEB" --ingest-port "$ING" --quiet 2> "${CASE_WORK}/site.err" &
SITE=$!
trap 'kill "$SITE" 2>/dev/null || true' EXIT
for _ in $(seq 100); do [ -s "${STATE}/ingest-token" ] && "$CURL" -s -o /dev/null "http://127.0.0.1:${ING}/" && break; sleep 0.1; done
TOK=$(cat "${STATE}/ingest-token") || fail "og live serve wrote no ingest token: $(cat "${CASE_WORK}/site.err")"
URL="http://127.0.0.1:${ING}"

# The token as a post carries it: user og, in the curl config on fd 3.
AUTH="user = \"og:${TOK}\""
# meta_command: the command run-all.sh recorded for og recheck --rerun.
meta_command() { awk -F'\t' '$1 == "meta" && $2 == "command" { print $3 }' "${G}/sample1/logs/run_status.tsv"; }

# --- 1. the site is set ----------------------------------------------------------------------
: > "$LOG"; : > "${CASE_WORK}/argv.log"
# GRIDSS fails on the fake docker, so run-all.sh exits 1; what counts is what it posted.
run_rc on env PATH="$P" GRIDSS=true OWNGENOME_LIVE_URL="$URL" OWNGENOME_LIVE_STATE_DIR="$STATE" OWNGENOME_LIVE_PAGE_URL="http://127.0.0.1:${WEB}/" "${SCRIPTS}/run-all.sh" sample1 male
ON_RC=$RC
[ "$ON_RC" -le 1 ] || fail "run-all.sh exited ${ON_RC}, expected 0 or 1"
grep -qxF "nextflow OWNGENOME_LIVE_RUN=sample1 OWNGENOME_LIVE_URL=${URL} OWNGENOME_LIVE_TOKEN=${TOK}" "$LOG" \
  || fail "nextflow did not see OWNGENOME_LIVE_RUN=sample1, the plain URL and the token in OWNGENOME_LIVE_TOKEN: $(cat "$LOG")"
grep -q "^post ${URL}/ingest/sample1 {\"event\":\"pgp_step\",\"step\":\"04b\",\"state\":\"started\"," "$LOG" \
  || fail "no 'started' post for step 04b: $(cat "$LOG")"
grep -qE "^post ${URL}/ingest/sample1 \{\"event\":\"pgp_step\",\"step\":\"04b\",\"state\":\"(ok|failed)\"," "$LOG" \
  || fail "no result post for step 04b: $(cat "$LOG")"
grep -q "^post ${URL}/ingest/sample1 {\"event\":\"pgp_run_finished\",\"exit\":${ON_RC}," "$LOG" \
  || fail "no pgp_run_finished post with exit ${ON_RC}: $(cat "$LOG")"
[ "$(grep -c '^post ' "$LOG")" -eq "$(grep -cxF "auth ${URL}/ingest/sample1 ${AUTH}" "$LOG")" ] \
  || fail "a post did not carry the token on fd 3: $(grep '^post \|^auth ' "$LOG")"
if grep -qF -- "$TOK" "${CASE_WORK}/argv.log"; then fail "curl was started with the token in its arguments"; fi
output_has on "Live site: messages go to ${URL}/ingest/sample1"
# The run's page, from OWNGENOME_LIVE_PAGE_URL with its trailing slash: one slash before runs/, and a page the site serves.
output_has on "Live site: this run's page is http://127.0.0.1:${WEB}/runs/sample1/"
output_lacks on "//runs/"
got=$("$CURL" -s -o /dev/null -w '%{http_code}' "http://127.0.0.1:${WEB}/runs/sample1/")
[ "$got" = 200 ] || fail "the run's page answered ${got}"
output_lacks on "$TOK"
LC="${G}/sample1/nextflow/live.config"
grep -q '^nextflow :: .* -c '"$LC" "$FAKE_DOCKER_LOG" || fail "nextflow was not given -c ${LC}: $(grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1)"
grep -qE "^    id 'nf-weblog@[0-9]+\.[0-9]+\.[0-9]+'$" "$LC" || fail "live.config does not load a pinned nf-weblog: $(cat "$LC")"
grep -qxF "    url     = '${URL}/ingest/sample1'" "$LC" || fail "live.config does not point at the site: $(cat "$LC")"
grep -qxF "    basicToken = \"og:\${env('OWNGENOME_LIVE_TOKEN')}\"" "$LC" || fail "live.config does not send the token as HTTP Basic: $(cat "$LC")"
if grep -qF -- "$TOK" "$LC"; then fail "live.config holds the token"; fi
case "$(ls -l "$LC")" in -rw-------*) ;; *) fail "live.config is not mode 600: $(ls -l "$LC")" ;; esac
[ "$(meta_command | grep -cF "OWNGENOME_LIVE_URL=${URL} ")" -eq 1 ] || fail "meta command does not name the plain URL: $(meta_command)"
if meta_command | grep -qF -- "$TOK"; then fail "meta command holds the token"; fi
# Every post is a message the site takes: each body, posted again with the real curl to og live serve.
sed -n 's/^post [^ ]* //p' "$LOG" > "${CASE_WORK}/bodies.jsonl"
n=0
while IFS= read -r body; do
  got=$(printf '%s' "$body" | "$CURL" -s -w ' %{http_code}' -H 'Content-Type: application/json' -d @- -K/dev/fd/3 \
    "${URL}/ingest/sample1" 3< <(printf '%s\n' "$AUTH"))
  [ "$(printf '%s' "$got" | tr -d '\n')" = "ok 200" ] || fail "og live serve did not take a post run-all.sh sent (${got}): ${body}"
  n=$((n + 1))
done < "${CASE_WORK}/bodies.jsonl"
[ "$n" -ge 3 ] || fail "only ${n} posts to replay: $(cat "$LOG")"
[ "$(wc -l < "${STATE}/runs/sample1/events.jsonl")" -eq "$n" ] || fail "events.jsonl does not hold the ${n} posts"
# Controls: the same post without the token, with another token, and a step the site does not know, are refused.
body=$(head -1 "${CASE_WORK}/bodies.jsonl")
got=$(printf '%s' "$body" | "$CURL" -s -o /dev/null -w '%{http_code}' -H 'Content-Type: application/json' -d @- "${URL}/ingest/sample1")
[ "$got" = 403 ] || fail "og live serve answered ${got} to a post without the token"
got=$(printf '%s' "$body" | "$CURL" -s -o /dev/null -w '%{http_code}' -H 'Content-Type: application/json' -d @- -K/dev/fd/3 \
  "${URL}/ingest/sample1" 3< <(printf 'user = "og:%s"\n' "${TOK%?}x"))
[ "$got" = 403 ] || fail "og live serve answered ${got} to a post with another token"
got=$(printf '%s' '{"event":"pgp_step","step":"4","state":"ok"}' | "$CURL" -s -o /dev/null -w '%{http_code}' -H 'Content-Type: application/json' -d @- \
  -K/dev/fd/3 "${URL}/ingest/sample1" 3< <(printf '%s\n' "$AUTH"))
[ "$got" = 400 ] || fail "og live serve answered ${got} to a step it does not know"
[ "$(wc -l < "${STATE}/runs/sample1/events.jsonl")" -eq "$n" ] || fail "a refused post reached events.jsonl"

# --- 2. no site --------------------------------------------------------------------------------
: > "$LOG"
run_rc off env PATH="$P" GRIDSS=true OWNGENOME_LIVE_PAGE_URL="http://127.0.0.1:${WEB}/" "${SCRIPTS}/run-all.sh" sample1 male
grep -q '^nextflow OWNGENOME_LIVE_RUN=unset OWNGENOME_LIVE_URL=unset OWNGENOME_LIVE_TOKEN=unset$' "$LOG" || fail "nextflow saw a live setting without OWNGENOME_LIVE_URL: $(cat "$LOG")"
if grep -q '^post ' "$LOG"; then fail "run-all.sh posted without OWNGENOME_LIVE_URL: $(cat "$LOG")"; fi
output_lacks off 'Live site'
if grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1 | grep -q 'live\.config'; then fail "nextflow was given live.config without OWNGENOME_LIVE_URL"; fi
[ ! -e "${G}/sample1/nextflow/live.config" ] || fail "an earlier run's live.config was left in place"
run_expect 2 badurl env PATH="$P" OWNGENOME_LIVE_URL="http://127.0.0.1:8788/x'y" "${SCRIPTS}/run-all.sh" sample1 male
output_has badurl 'OWNGENOME_LIVE_URL must be a plain'
# The ingest listener is plain HTTP and takes /ingest/<run> with the token: https, another path or a token
# of the wrong shape would lose every post. Each is refused before run_status.tsv is written, without its value.
STATUS_BEFORE=$(cat "${G}/sample1/logs/run_status.tsv")
run_expect 2 httpsurl env PATH="$P" OWNGENOME_LIVE_URL="https://127.0.0.1:8788" "${SCRIPTS}/run-all.sh" sample1 male
output_has httpsurl 'OWNGENOME_LIVE_URL must be a plain'
run_expect 2 pathurl env PATH="$P" OWNGENOME_LIVE_URL="http://127.0.0.1:8788/live" "${SCRIPTS}/run-all.sh" sample1 male
output_has pathurl 'OWNGENOME_LIVE_URL must be a plain'
run_expect 2 shorttoken env PATH="$P" OWNGENOME_LIVE_URL="${URL}/${TOK%?}" "${SCRIPTS}/run-all.sh" sample1 male
output_has shorttoken 'OWNGENOME_LIVE_URL must be a plain'
output_lacks shorttoken "${TOK%?}"
run_expect 2 badtoken env PATH="$P" OWNGENOME_LIVE_URL="$URL" OWNGENOME_LIVE_TOKEN="${TOK}x" "${SCRIPTS}/run-all.sh" sample1 male
output_has badtoken 'OWNGENOME_LIVE_TOKEN is not an ingest token'
output_lacks badtoken "${TOK}x"
[ "$(cat "${G}/sample1/logs/run_status.tsv")" = "$STATUS_BEFORE" ] || fail "a refused live setting rewrote run_status.tsv"

# --- 3. the site refuses every post -------------------------------------------------------------
# The token in the URL itself, as when the site runs as another account: the warning names the site without it.
: > "$LOG"; : > "${CASE_WORK}/argv.log"
run_rc down env PATH="$P" GRIDSS=true OWNGENOME_LIVE_URL="${URL}/${TOK}" FAKE_DOWNLOAD_FAIL='/ingest/' "${SCRIPTS}/run-all.sh" sample1 male
[ "$RC" -eq "$ON_RC" ] || fail "with the site down run-all.sh exited ${RC}, with it up ${ON_RC}"
output_has down "WARNING: the live site at ${URL} did not take a message"
output_has down "Live site: messages go to ${URL}/ingest/sample1"
output_lacks down "$TOK"
output_lacks down "this run's page"
grep -qxF "auth ${URL}/ingest/sample1 ${AUTH}" "$LOG" || fail "the token of the URL did not go with the posts: $(grep '^auth ' "$LOG")"
if grep -qF -- "$TOK" "${CASE_WORK}/argv.log"; then fail "curl was started with the token of the URL in its arguments"; fi
[ "$(meta_command | grep -cF "OWNGENOME_LIVE_URL=${URL} ")" -eq 1 ] || fail "meta command does not name the plain URL: $(meta_command)"
if meta_command | grep -qF -- "$TOK"; then fail "meta command holds the token of the URL"; fi
if grep -qF -- "$TOK" "${G}/sample1/nextflow/live.config"; then fail "live.config holds the token of the URL"; fi
echo "run-all.sh posts its steps when OWNGENOME_LIVE_URL is set, nothing without it, and a site that is down changes nothing."
