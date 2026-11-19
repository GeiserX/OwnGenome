#!/usr/bin/env bash
# bin/write_manifest.sh with two or more writers of one sample at once (a run
# and a report script can both write run_manifest.tsv). A fake container
# engine holds each writer at its first image lookup, after it has started
# its manifest, until every writer is there; then each lookup sleeps a little.
# So the writers always overlap, however the runner schedules them.
# Checked: every writer exits 0; the manifest is whole (three tab-separated
# fields per line, each section/key once, the same keys as a lone writer
# writes); no temporary file is left beside it. Red control: a copy of the
# script that builds the manifest in one shared temporary file
# (run_manifest.tsv.tmp, the old name) must leave a manifest that is not
# whole under the same overlap.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fails=0
ok()  { echo "PASS $*"; }
bad() { echo "FAIL $*"; fails=$((fails + 1)); }

# A container engine that knows no image: write_manifest.sh asks it twice per
# *_IMAGE line of versions.env. The barrier: the first call of writer
# WRITER_ID marks BARRIER/WRITER_ID and waits (30 s at most) until WRITERS
# writers have marked it. A writer's first call comes after it has written
# its run and version rows into its temporary file.
mkdir -p "${T}/fake"
cat > "${T}/fake/engine" <<'ENGINE'
#!/bin/sh
if [ ! -e "${BARRIER}/${WRITER_ID}" ]; then
  : > "${BARRIER}/${WRITER_ID}"
  i=0
  while [ "$(find "$BARRIER" -type f | wc -l)" -lt "$WRITERS" ] && [ "$i" -lt 600 ]; do
    sleep 0.05
    i=$((i + 1))
  done
fi
sleep 0.01
exit 0
ENGINE
chmod 755 "${T}/fake/engine"

# whole MANIFEST REF: MANIFEST has three fields per line and the same
# section/key rows as REF, each as often (a row a writer repeated or lost
# shows in the diff). Prints what is wrong.
whole() {
  local m=$1 ref=$2
  [ -s "$m" ] || { echo "missing or empty"; return 1; }
  awk -F'\t' 'NF != 3 {print "line " NR " has " NF " fields"; e = 1} END {exit e}' "$m" || return 1
  local diff
  diff=$(diff <(cut -f1,2 "$ref" | sort) <(cut -f1,2 "$m" | sort) | grep '^[<>]' | head -3)
  [ -z "$diff" ] || { echo "rows differ from a lone writer's: ${diff//$'\t'/ }"; return 1; }
}

# writers SCRIPT GENOME_DIR N: N writers of sample s1 at once, held together
# at the barrier; prints their exit codes.
writers() {
  local script=$1 gd=$2 n=$3 i pids=() rc codes=""
  mkdir -p "${gd}/s1" "${gd}/barrier"
  for ((i = 1; i <= n; i++)); do
    CONTAINER_ENGINE="${T}/fake/engine" BARRIER="${gd}/barrier" WRITER_ID="$i" WRITERS="$n" \
      GENOME_DIR="$gd" bash "$script" s1 "writer${i}" > /dev/null 2> "${gd}/err${i}" &
    pids+=("$!")
  done
  for i in "${pids[@]}"; do
    rc=0
    wait "$i" || rc=$?
    codes="${codes}${rc} "
  done
  echo "$codes"
}

# A lone writer, for the set of rows.
mkdir -p "${T}/lone/s1"
if CONTAINER_ENGINE=true GENOME_DIR="${T}/lone" bash "${REPO}/bin/write_manifest.sh" s1 lone > /dev/null 2> "${T}/lone/err"; then
  ok "a lone writer writes run_manifest.tsv ($(($(wc -l < "${T}/lone/s1/run_manifest.tsv"))) rows)"
else
  bad "a lone writer failed: $(cat "${T}/lone/err")"
fi
REF="${T}/lone/s1/run_manifest.tsv"

# 1. the repository's script, three writers at once
codes=$(writers "${REPO}/bin/write_manifest.sh" "${T}/g" 3)
if [ "$codes" = "0 0 0 " ]; then
  ok "three writers at once all exit 0"
else
  bad "three writers at once exit ${codes}: $(cat "${T}"/g/err*)"
fi
if why=$(whole "${T}/g/s1/run_manifest.tsv" "$REF"); then
  ok "the manifest of three writers at once is whole"
else
  bad "the manifest of three writers at once is not whole: ${why}"
fi
left=$(find "${T}/g/s1" -mindepth 1 ! -name run_manifest.tsv)
if [ -z "$left" ]; then
  ok "no temporary file is left beside the manifest"
else
  bad "left beside the manifest: ${left}"
fi
written_by=$(awk -F'\t' '$1 == "run" && $2 == "written_by" {print $3}' "${T}/g/s1/run_manifest.tsv")
if [[ "$written_by" =~ ^writer[123]$ ]]; then
  ok "the manifest is one writer's (${written_by})"
else
  bad "written_by is '${written_by}'"
fi

# 2. red control: the old shared temporary name, in a copy of the files the
# script reads (bin/, scripts/lib/, versions.env, nextflow.config).
C="${T}/control"
mkdir -p "${C}/bin" "${C}/scripts"
cp "${REPO}/bin/write_manifest.sh" "${C}/bin/"
cp -R "${REPO}/scripts/lib" "${C}/scripts/"
cp "${REPO}/versions.env" "${REPO}/nextflow.config" "${C}/"
# shellcheck disable=SC2016  # the planted lines are text, never expanded here
sed -i.bak -E -e 's/^TMP=\$\(mktemp .*$/TMP="${OUT}.tmp"; : > "$TMP"/' \
  -e '/^trap .rm -f "\$TMP". EXIT$/d' "${C}/bin/write_manifest.sh"
# shellcheck disable=SC2016  # as above
if grep -qx 'TMP="${OUT}.tmp"; : > "$TMP"' "${C}/bin/write_manifest.sh" \
  && ! grep -q '^trap' "${C}/bin/write_manifest.sh"; then
  codes=$(writers "${C}/bin/write_manifest.sh" "${T}/gc" 3)
  if why=$(whole "${T}/gc/s1/run_manifest.tsv" "$REF"); then
    bad "control: three writers sharing run_manifest.tsv.tmp left a whole manifest (exits ${codes}); this test cannot fail"
  else
    ok "control: three writers sharing run_manifest.tsv.tmp break the manifest (${why})"
  fi
else
  bad "control: could not plant the shared temporary name in the copy"
fi

echo "${fails} failure(s)"
[ "$fails" -eq 0 ]
