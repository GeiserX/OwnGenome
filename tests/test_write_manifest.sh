#!/usr/bin/env bash
# og-leg: OG_REQUIRED
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
#
# The run intervals row (INTERVALS, run-all.sh's --intervals): INTERVALS set
# writes one row, its regions on one line with single spaces; INTERVALS unset,
# blank, genome or chip (step 29's whole genome and CHIP genes) writes none. Then og summary collect reads the manifest the script
# wrote, beside a sample VCF with calls on all 22 autosomes: the row marks the
# run partial ("limited to the run's intervals"). Red controls: a copy of the
# script without the row line writes no row with INTERVALS set, and the same
# folder without the row is not marked. The og part needs og (OG, else
# bin/og); without it that part is skipped with a note, and with OG_REQUIRED
# set the skip is a failure.
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fails=0
# The run intervals row comes from INTERVALS: only section 3 sets it.
unset INTERVALS
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

# 3. the run intervals row. manifest SCRIPT GENOME_DIR [INTERVALS]: the
# manifest SCRIPT writes for sample s1, with INTERVALS set when given.
manifest() {
  local script=$1 gd=$2
  mkdir -p "${gd}/s1"
  if [ "$#" -gt 2 ]; then
    INTERVALS=$3 CONTAINER_ENGINE=true GENOME_DIR="$gd" bash "$script" s1 test > /dev/null 2> "${gd}/err"
  else
    CONTAINER_ENGINE=true GENOME_DIR="$gd" bash "$script" s1 test > /dev/null 2> "${gd}/err"
  fi || bad "${script} failed for ${gd}: $(cat "${gd}/err")"
}
intervals_rows() { awk -F'\t' '$1 == "run" && $2 == "intervals"' "$1/s1/run_manifest.tsv"; }
want=$'run\tintervals\tchr1:1-1000000 chr2:1-1000000'
manifest "${REPO}/bin/write_manifest.sh" "${T}/iv" $' chr1:1-1000000\t chr2:1-1000000\n'
rows=$(intervals_rows "${T}/iv")
if [ "$rows" = "$want" ]; then
  ok "INTERVALS set writes one run intervals row, on one line"
else
  bad "INTERVALS set wrote the rows '${rows//$'\t'/ }', want '${want//$'\t'/ }'"
fi
manifest "${REPO}/bin/write_manifest.sh" "${T}/noiv"
manifest "${REPO}/bin/write_manifest.sh" "${T}/blankiv" "  "
# genome and chip are step 29's own INTERVALS values (the whole genome, the
# CHIP genes), not regions the run was limited to.
manifest "${REPO}/bin/write_manifest.sh" "${T}/genomeiv" " genome "
manifest "${REPO}/bin/write_manifest.sh" "${T}/chipiv" "chip"
for gd in noiv blankiv genomeiv chipiv; do
  rows=$(intervals_rows "${T}/${gd}")
  if [ -z "$rows" ]; then
    ok "${gd}: no run intervals row (a whole-genome run)"
  else
    bad "${gd}: wrote '${rows//$'\t'/ }'"
  fi
done
# Red control: a copy of the script without the row line, INTERVALS set.
N="${T}/norow"
mkdir -p "${N}/bin" "${N}/scripts"
cp "${REPO}/bin/write_manifest.sh" "${N}/bin/"
cp -R "${REPO}/scripts/lib" "${N}/scripts/"
cp "${REPO}/versions.env" "${REPO}/nextflow.config" "${N}/"
sed -i.bak '/row run intervals /d' "${N}/bin/write_manifest.sh"
if grep -q 'row run intervals ' "${N}/bin/write_manifest.sh"; then
  bad "control: could not drop the row line from the copy"
else
  manifest "${N}/bin/write_manifest.sh" "${T}/norowiv" "chr1:1-1000000 chr2:1-1000000"
  if [ -z "$(intervals_rows "${T}/norowiv")" ]; then
    ok "control: the copy without the row line writes no row with INTERVALS set"
  else
    bad "control: the copy without the row line still wrote it; the check above cannot fail"
  fi
fi

# 4. og summary collect reads the row: a sample VCF with calls on all 22
# autosomes and the manifest of 3 is marked partial; without the row it is not.
OG=${OG:-${REPO}/bin/og}
if ! "$OG" version > /dev/null 2>&1; then
  if [ -n "${OG_REQUIRED:-}" ]; then
    bad "og is not built (scripts/setup.sh --og), and OG_REQUIRED is set"
  else
    echo "skip: og is not built for this machine (scripts/setup.sh --og): the run intervals row is not read back"
  fi
else
  vcf="${T}/iv/s1/vcf/s1.vcf.gz"
  mkdir -p "$(dirname "$vcf")"
  {
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
    for i in $(seq 1 22); do printf 'chr%s\t1000000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n' "$i"; done
  } | gzip > "$vcf"
  mark='"mark": "partial: limited to the run'"'"'s intervals (chr1:1-1000000 chr2:1-1000000)"'
  # scope SAMPLE_DIR: og summary collect writes SAMPLE_DIR/summary.json
  # (emptied first, so a failed run leaves nothing to read).
  scope() {
    : > "$1/summary.json"
    "$OG" summary collect --sample s1 --sample-dir "$1" --assets "${REPO}/assets" --out "$1/summary.json" > /dev/null 2> "$1/og.err" \
      || bad "og summary collect failed: $(cat "$1/og.err")"
  }
  scope "${T}/iv/s1"
  json=$(cat "${T}/iv/s1/summary.json")
  if grep -qF "$mark" <<< "$json" && grep -qF '"autosomes_covered": 22' <<< "$json"; then
    ok "og summary collect marks the run partial from the row (22 autosomes covered)"
  else
    bad "og summary collect did not mark the run partial from the row: $(grep -A1 -F '"mark"' <<< "$json" | head -2)"
  fi
  # Red control: the same folder, its manifest without the row.
  awk -F'\t' '!($1 == "run" && $2 == "intervals")' "${T}/iv/s1/run_manifest.tsv" > "${T}/iv/m" \
    && mv "${T}/iv/m" "${T}/iv/s1/run_manifest.tsv"
  [ -z "$(intervals_rows "${T}/iv")" ] || bad "control: could not drop the row from the manifest"
  scope "${T}/iv/s1"
  json=$(cat "${T}/iv/s1/summary.json")
  if grep -qF '"autosomes_covered": 22' <<< "$json" && ! grep -qF "$mark" <<< "$json" \
    && grep -qF '"mark": null' <<< "$json"; then
    ok "control: without the row the same run is not marked partial"
  else
    bad "control: without the row the run is still marked; the check above cannot fail"
  fi
fi

# 5. The Nextflow side: HTML_REPORT (modules/local/html_report/main.nf) writes
# the same row from --intervals. No Nextflow here, so the test does what
# Nextflow would: it evaluates the module's def intervals (toString, then each
# .replace("'", "'\\''"); any other step is refused, so a new step needs a new
# test), puts the result in the task script's intervals lines, Groovy's
# escapes undone, and runs those lines in bash. A value with a tab and a
# newline must give the row section 3 wants; a value with a quote must reach
# the row as it is. Red controls: module copies without the one-line step,
# without the quote escape and without the printf line must not.
# HTML_REPORT_MODULE overrides the module.
MOD=${HTML_REPORT_MODULE:-${REPO}/modules/local/html_report/main.nf}
# nf_def MODULE VALUE: what the module's def intervals gives for --intervals
# VALUE (a Groovy '' is false, and both sides give '').
nf_def() {
  local expr v=$2 q="'" r="'\\''"
  local head='params.intervals ? params.intervals.toString()' tail=" : ''"
  local esc=".replace(\"'\", \"'\\\\''\")"
  expr=$(sed -n -E 's/^[[:space:]]*def intervals = (.*)$/\1/p' "$1")
  [[ "$expr" == "$head"*"$tail" ]] || { echo "def intervals is not '${head}...${tail}': '${expr}'"; return 1; }
  expr=${expr#"$head"}
  expr=${expr%"$tail"}
  while [ -n "$expr" ]; do
    [[ "$expr" == "$esc"* ]] || { echo "def intervals has a step this test cannot evaluate: '${expr}'"; return 1; }
    v=${v//"$q"/"$r"}
    expr=${expr#"$esc"}
  done
  printf '%s' "$v"
}
# nf_rows MODULE VALUE: the rows the module's intervals lines write for VALUE.
nf_rows() {
  local d val
  d=$(mktemp -d "${T}/nf.XXXXXX")
  nf_def "$1" "$2" > "${d}/val" || { cat "${d}/val"; return; }
  val=$(cat "${d}/val"; printf x)
  val=${val%x}
  # The lines from intervals='${intervals}' to the printf of the row, as bash
  # sees them: \$ is $, \\ is \ (any other escape is refused), then the def's
  # value in place of ${intervals}, as data.
  # shellcheck disable=SC2016  # the dollar signs are awk's and the script's
  if ! OG_VAL="$val" awk '
    index($0, "intervals='"'"'${intervals}'"'"'") { on = 1 }
    on {
      out = ""
      for (i = 1; i <= length($0); i++) {
        c = substr($0, i, 1)
        if (c == "\\") {
          nxt = substr($0, i + 1, 1)
          if (nxt != "\\" && nxt != "$") { print "unexpected Groovy escape \\" nxt > "/dev/stderr"; exit 1 }
          out = out nxt; i++
        } else out = out c
      }
      k = index(out, "${intervals}")
      if (k) out = substr(out, 1, k - 1) ENVIRON["OG_VAL"] substr(out, k + 12)
      print out
      if (index($0, "printf") && index($0, "run\\\\tintervals")) { found = 1; exit }
    }
    END { exit !found }' "$1" > "${d}/block.sh"; then
    echo "no intervals lines"
    return
  fi
  (cd "$d" && S="$d" bash -euo pipefail "${d}/block.sh") 2> "${d}/err" || { echo "the lines failed: $(cat "${d}/err")"; return; }
  [ ! -e "${d}/run_manifest.tsv" ] || awk -F'\t' '$1 == "run" && $2 == "intervals"' "${d}/run_manifest.tsv"
}
messy=$' chr1:1-1000000\t chr2:1-1000000\n'
quoted="chr1:1-1000000 it's"
want_quoted=$'run\tintervals\t'"$quoted"
rows=$(nf_rows "$MOD" "$messy")
if [ "$rows" = "$want" ]; then
  ok "HTML_REPORT writes the row write_manifest.sh writes, --intervals with a tab and a newline on one line"
else
  bad "HTML_REPORT wrote '${rows//$'\t'/ }', want '${want//$'\t'/ }'"
fi
rows=$(nf_rows "$MOD" "$quoted")
if [ "$rows" = "$want_quoted" ]; then
  ok "HTML_REPORT writes a quote in --intervals as it is"
else
  bad "HTML_REPORT wrote '${rows//$'\t'/ }' for a quote, want '${want_quoted//$'\t'/ }'"
fi
rows=$(nf_rows "$MOD" "")
if [ -z "$rows" ]; then
  ok "HTML_REPORT without --intervals writes no row"
else
  bad "HTML_REPORT without --intervals wrote '${rows//$'\t'/ }'"
fi
# Red controls: module copies that each break one step.
# shellcheck disable=SC2016  # the module's text, not expansions
sed 's/"\\${regions\[\*\]}" >>/"\\$intervals" >>/' "$MOD" > "${T}/html_report_rawiv.nf"
# shellcheck disable=SC2016  # as above
if ! grep -qF '"\$intervals" >>' "${T}/html_report_rawiv.nf"; then
  bad "control: could not drop the one-line step from the module copy"
elif [ "$(nf_rows "${T}/html_report_rawiv.nf" "$messy")" != "$want" ]; then
  ok "control: the module without the one-line step does not write the row"
else
  bad "control: the module without the one-line step still wrote the row; the check above cannot fail"
fi
sed '/def intervals = /s/\.replace([^)]*)//' "$MOD" > "${T}/html_report_noesc.nf"
if grep -E 'def intervals = .*replace' "${T}/html_report_noesc.nf" > /dev/null; then
  bad "control: could not drop the quote escape from the module copy"
elif [ "$(nf_rows "${T}/html_report_noesc.nf" "$quoted")" != "$want_quoted" ]; then
  ok "control: the module without the quote escape does not write a quote as it is"
else
  bad "control: the module without the quote escape still wrote the quote; the check above cannot fail"
fi
sed '/printf .run\\\\tintervals/d' "$MOD" > "${T}/html_report_norow.nf"
if grep -q 'run\\\\tintervals' "${T}/html_report_norow.nf"; then
  bad "control: could not drop the printf line from the module copy"
elif [ "$(nf_rows "${T}/html_report_norow.nf" "chr1:1-1000000 chr2:1-1000000")" != "$want" ]; then
  ok "control: the module without the printf line does not write the row"
else
  bad "control: the module without the printf line still wrote a row; the check above cannot fail"
fi

echo "${fails} failure(s)"
[ "$fails" -eq 0 ]
