#!/usr/bin/env bash
# Public synthetic reports and traces. No containers or downloaded fixture needed.
. "$(dirname "$0")/lib.sh"
printf '{"genes":{"CYP2C19":{"sourceDiplotypes":[{"allele1":{"name":"*1"}}]}}}\n' > "$CASE_TMP/report.json"
check_pharmcat_called "a readable report has a named gene" "$CASE_TMP/report.json"
printf '{"genes":' > "$CASE_TMP/broken.json"
rc=0
pharmcat_called "$CASE_TMP/broken.json" > "$CASE_TMP/out" 2> "$CASE_TMP/err" || rc=$?
check_eq "truncated report exits 2" "$rc" 2
check "truncated report is unreadable, never a zero count" grep -q 'unreadable report:' "$CASE_TMP/err"
check_eq "truncated report has no call output" "$(wc -c < "$CASE_TMP/out" | tr -d ' ')" 0
printf 'task_id\tname\tattempt\texit\tmemory\n1\tPRS (s)\t1\t137\t25769803776\n2\tPRS (s)\t2\t0\t51539607552\n' > "$CASE_TMP/trace.txt"
check "retry memory and exit match" "${OG:-${REPO}/bin/og}" dev nf-trace attempts "$CASE_TMP/trace.txt" --name-glob 'PRS (s)' --want '1:137:24,2:0:48' --nextflow-exit 0
rc=0
"${OG:-${REPO}/bin/og}" dev nf-trace attempts "$CASE_TMP/trace*.txt" --name-glob 'PRS (s)' --want '1:137:24,2:0:48' --nextflow-exit 0 > "$CASE_TMP/out" 2> "$CASE_TMP/err" || rc=$?
check_eq "ambiguous trace glob exits 2" "$rc" 2
check "trace glob refusal explains the required path" grep -q 'trace glob refused' "$CASE_TMP/err"
finish
