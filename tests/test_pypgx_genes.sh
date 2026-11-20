#!/usr/bin/env bash
# og-leg: OG_REQUIRED
# Step 32 (pypgx) takes its genes from assets/pgx_genes.tsv. This test checks:
#  1. the table against pypgx's own gene table at the pinned tag: every bam
#     row is a gene pypgx marks SV=TRUE, and every SV=TRUE gene is a bam row
#     (a planted SLC28A3 bam row and UGT2B17 set to vcf must both be reported);
#  2. scripts/32-pypgx.sh under a fake docker, with a fake pypgx that logs its
#     calls: the genes the step ran before the table (four bam, 19 vcf), given
#     as a table, give the same calls in the same order; each bam row added to
#     that table adds exactly one BAM-mode call; the shipped table gives one
#     BAM-mode call per bam row;
#  3. the not-callable line: from the callable_without_alt column, printed for
#     a BAM whose header lacks the gene's ALT contig and not for one that has
#     it, and not printed when the column says yes;
#  4. the PYPGX module reads the same table and spells no gene list;
#  5. the reports (og summary collect, then og report render): every bam
#     gene has exactly one row in the PGx copy number block (its call, or why
#     there is none), pypgx's, even beside a CNVpytor section, and a
#     not-callable gene gets its not-assessable line
#     (two planted summaries must be reported);
#  6. one owner per gene: no row of the depth engine's catalog
#     (bin/cn_loci.tsv) is a gene pypgx calls from the BAM (a planted UGT2B17
#     row must be reported);
#  7. no old call survives a rerun: a results.zip left by an earlier run is
#     gone when depth preparation fails, a gene whose pypgx run fails after
#     writing part of its result leaves no results.zip, and a gene that
#     succeeds has its result under <gene>/;
#  8. an ALT contig check that cannot run (the BAM header cannot be read, or
#     og cn check-bam fails) stops the step with exit 2 before pypgx runs;
#  9. docs/32-pypgx.md names exactly the genes pypgx has no phenotype rule
#     for (PhenotypeMethod N/A in its gene table) and, of those, the weakly
#     supported ones no other tool here calls (second_caller none).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
CAT="${REPO}/assets/pgx_genes.tsv"
FAILS=0
bad() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
ok() { echo "ok: $*"; }

# The genes and modes of scripts/32-pypgx.sh before the table, as it spelled
# them (BAM_GENES, VCF_GENES).
OLD_BAM="CYP2D6 CYP2A6 GSTM1 GSTT1"
OLD_VCF="CYP1A2 CYP2B6 CYP2C9 CYP2C19 CYP3A4 CYP3A5 CYP4F2 DPYD TPMT NUDT15 UGT1A1 SLCO1B1 VKORC1 NAT2 COMT MTHFR ABCB1 G6PD IFNL3"

rows() { awk -F'\t' '!/^#/ && $1 != "gene" && NF' "$1"; }
row_of() { awk -F'\t' -v g="$1" '$1 == g' "$CAT"; }
# set_col GENE N VALUE: the shipped table with column N of GENE's row set to VALUE
set_col() { awk -F'\t' -v OFS='\t' -v g="$1" -v n="$2" -v v="$3" '$1 == g {$n = v} {print}' "$CAT"; }

# --- 1. the table against pypgx's gene table ---------------------------------
# shellcheck source=../versions.env
. "${REPO}/versions.env"
PYPGX_VERSION=${PYPGX_IMAGE##*:}
PYPGX_VERSION=${PYPGX_VERSION%%--*}
[ "$PYPGX_VERSION" = "$PYPGX_BUNDLE_VERSION" ] || bad "PYPGX_IMAGE is ${PYPGX_VERSION}, PYPGX_BUNDLE_VERSION ${PYPGX_BUNDLE_VERSION}"
GT="${WORK}/gene-table.csv"
curl -sS --fail --retry 4 --retry-delay 5 --max-time 60 -o "$GT" \
  "https://raw.githubusercontent.com/sbslee/pypgx/${PYPGX_VERSION}/pypgx/api/data/gene-table.csv"

# check_table TABLE: 0 when the bam rows are exactly pypgx's SV=TRUE genes and
# every row is a pypgx target gene
check_table() {
  # pypgx's gene table is CSV; the columns read here (Gene, Target, SV) come
  # before its first quoted field, so a comma split reads them.
  awk -F'\t' '
    NR == FNR {
      n = split($0, f, ",")
      if (FNR == 1) { for (i = 1; i <= n; i++) col[f[i]] = i; if (!col["Gene"] || !col["Target"] || !col["SV"]) { print "gene-table.csv has no Gene, Target or SV column"; bad = 2; exit 2 }; next }
      g = f[col["Gene"]]; target[g] = f[col["Target"]]
      if (f[col["SV"]] == "TRUE") sv[g] = 1; else delete sv[g]
      next
    }
    !/[^ \t\r]/ || /^#/ { next }
    !header++ { next }
    { rows[++nr] = $1; if ($2 == "bam") bam[$1] = 1 }
    END {
      if (bad) exit bad
      for (g in sv) nsv++
      for (g in bam) nbam++
      for (i = 1; i <= nr; i++) if (target[rows[i]] != "TRUE") errs[++ne] = rows[i] " is not a pypgx target gene"
      m = 0; for (g in bam) if (!(g in sv)) x[++m] = g
      sort(x, m); for (i = 1; i <= m; i++) errs[++ne] = x[i] " is a bam row but pypgx has no copy-number model for it (SV=FALSE)"
      m = 0; for (g in sv) if (!(g in bam)) x[++m] = g
      sort(x, m); for (i = 1; i <= m; i++) errs[++ne] = x[i] " has a copy-number model (SV=TRUE) but is not a bam row"
      if (!ne) { print (nbam + 0) " bam rows = pypgx'\''s " (nsv + 0) " SV genes"; exit 0 }
      for (i = 1; i <= ne; i++) print errs[i]
      exit 1
    }
    function sort(a, n,   i, j, t) { for (i = 2; i <= n; i++) { t = a[i]; for (j = i - 1; j >= 1 && a[j] > t; j--) a[j + 1] = a[j]; a[j + 1] = t } }
  ' "$GT" "$1"
}
if out=$(check_table "$CAT"); then ok "$out"; else bad "the shipped table: $out"; fi
{ cat "$CAT"; printf 'SLC28A3\tbam\tyes\t.\tnone\tno\n'; } > "${WORK}/planted-slc28a3.tsv"
if out=$(check_table "${WORK}/planted-slc28a3.tsv"); then bad "a planted SLC28A3 bam row passed"; else ok "control, planted SLC28A3 bam row reported: $out"; fi
set_col UGT2B17 2 vcf > "${WORK}/planted-ugt2b17-vcf.tsv"
if out=$(check_table "${WORK}/planted-ugt2b17-vcf.tsv"); then bad "UGT2B17 as a vcf row passed"; else ok "control, UGT2B17 as vcf reported: $out"; fi

# --- 2 and 3. scripts/32-pypgx.sh under a fake docker -------------------------
G="${WORK}/genome" S=sample1
mkdir -p "${G}/${S}/aligned" "${G}/${S}/vcf" "${WORK}/bin" "${WORK}/pypgx-bin"
for f in "${S}/aligned/${S}_sorted.bam.bai" "${S}/vcf/${S}.vcf.gz" "${S}/vcf/${S}.vcf.gz.tbi"; do echo placeholder > "${G}/${f}"; done
B="${G}/reference/pypgx-bundle"
mkdir -p "$B"
git -C "$B" init -q
git -C "$B" -c core.hooksPath=/dev/null -c user.name=test -c user.email=test@example.invalid commit -q --allow-empty -m bundle
git -C "$B" tag "$PYPGX_BUNDLE_VERSION"
ln -s "${REPO}/scripts/ci/fake-docker/docker" "${WORK}/bin/docker"
# The bash running this test runs the step too (it needs bash 4.4 or later).
BASH_DIR=$(dirname "$BASH")

# A fake pypgx: logs each call, and run-ngs-pipeline writes a results.zip
# in its output directory. FAKE_PYPGX_FAIL=prep fails depth preparation;
# FAKE_PYPGX_FAIL=GENE fails that gene after writing part of its result.
cat > "${WORK}/pypgx-bin/pypgx" <<'SH'
#!/usr/bin/env bash
echo "$*" >> "$PYPGX_LOG"
case "$1" in
  prepare-depth-of-coverage) [ "${FAKE_PYPGX_FAIL:-}" != prep ] || exit 1 ;;
  run-ngs-pipeline)
    mkdir -p "$3"
    if [ "${FAKE_PYPGX_FAIL:-}" = "$2" ]; then echo partial > "$3/results.zip"; exit 1; fi
    echo new > "$3/results.zip" ;;
esac
exit 0
SH
chmod +x "${WORK}/pypgx-bin/pypgx"
# The docker run hook: og bam contigs reads the @SQ lines of the BAM
# placeholder (the header fixture); samtools view -H prints it, or fails
# with FAKE_HEADER_FAIL; og cn check-bam passes, or fails (exit 3) with
# FAKE_CHECK_FAIL; the pypgx container's bash -c runs with the fake pypgx,
# its /genome paths on the host; the summary call logs its arguments and
# writes a header.
cat > "${WORK}/hook" <<'SH'
#!/usr/bin/env bash
img=$1; shift
host() { printf '%s' "${GENOME_DIR}${1#/genome}"; }
if [ "$1" = /pgp-bin/og ] && [ "$2" = bam ] && [ "$3" = contigs ]; then
  awk -F'\t' '$1 == "@SQ" {sub(/^SN:/, "", $2); print $2}' "$(host "$4")"
elif [ "$1" = samtools ] && [ "$2" = view ] && [ "$3" = -H ]; then
  [ -z "${FAKE_HEADER_FAIL:-}" ] || exit 1
  cat "$(host "$4")"
elif [ "$1" = /pgp-bin/og ] && [ "$2" = cn ] && [ "$3" = check-bam ]; then
  [ -z "${FAKE_CHECK_FAIL:-}" ] || exit 3
elif [ "$1" = bash ] && [ "$2" = -c ]; then
  PATH="${PYPGX_BIN}:$PATH" bash -c "${3//\/genome\//${GENOME_DIR}/}"
elif [ "$1" = /pgp-bin/og ] && [ "$2" = pgx ] && [ "$3" = pypgx-summary ]; then
  printf '%s\n' "$*" >> "$SUMMARY_LOG"
  while [ $# -gt 0 ]; do [ "$1" = --out ] && printf 'Gene\n' > "$(host "$2")"; shift; done
fi
exit 0
SH
chmod +x "${WORK}/hook"

# run32 TABLE HEADER NAME: run the step with TABLE and a BAM whose header is
# HEADER (alt or noalt); leaves NAME.out, NAME.pypgx and NAME.summary
run32() {
  local hdr=$'@HD\tVN:1.6\n@SQ\tSN:chr22\tLN:50818468'
  [ "$2" = noalt ] || hdr+=$'\n@SQ\tSN:chr22_KI270879v1_alt\tLN:304548'
  printf '%s\n' "$hdr" > "${G}/${S}/aligned/${S}_sorted.bam"
  : > "${WORK}/$3.pypgx"; : > "${WORK}/$3.summary"
  local rc=0
  env -i HOME="$WORK" PATH="${WORK}/bin:${BASH_DIR}:/usr/bin:/bin" GENOME_DIR="$G" PGX_GENES="$1" \
    FAKE_DOCKER_LOG="${WORK}/docker.log" FAKE_DOCKER_RUN_HOOK="${WORK}/hook" PYPGX_BIN="${WORK}/pypgx-bin" \
    PYPGX_LOG="${WORK}/$3.pypgx" SUMMARY_LOG="${WORK}/$3.summary" FAKE_PYPGX_FAIL="${FAKE_PYPGX_FAIL:-}" \
    FAKE_HEADER_FAIL="${FAKE_HEADER_FAIL:-}" FAKE_CHECK_FAIL="${FAKE_CHECK_FAIL:-}" \
    "$BASH" "${REPO}/scripts/32-pypgx.sh" "$S" > "${WORK}/$3.out" 2>&1 || rc=$?
  RC=$rc
}
# bam_calls NAME, vcf_calls NAME: the genes of the BAM-mode and VCF-only
# run-ngs-pipeline calls, in order; prep NAME: the genes depth was prepared for
bam_calls() { awk '$1 == "run-ngs-pipeline" && /--depth-of-coverage/ {print $2}' "${WORK}/$1.pypgx" | paste -sd ' ' -; }
vcf_calls() { awk '$1 == "run-ngs-pipeline" && !/--depth-of-coverage/ {print $2}' "${WORK}/$1.pypgx" | paste -sd ' ' -; }
prep() { sed -n 's/^prepare-depth-of-coverage .* --genes //p' "${WORK}/$1.pypgx"; }
sum_arg() { sed -n "s/.* --$2 \([^-]*[^ -]\) --.*/\1/p" "${WORK}/$1.summary" | head -n 1; }
expect() { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1: got '$2', expected '$3'"; fi; }

# The table of the genes the step ran before: their rows of the shipped table
# with the mode they had, in the order they had.
old_table() {
  grep '^#' "$CAT"; grep '^gene	' "$CAT"
  local g
  for g in $OLD_BAM; do row_of "$g" | awk -F'\t' -v OFS='\t' '{$2 = "bam"; print}'; done
  for g in $OLD_VCF; do row_of "$g" | awk -F'\t' -v OFS='\t' '{$2 = "vcf"; print}'; done
}
old_table > "${WORK}/old.tsv"
[ "$(rows "${WORK}/old.tsv" | wc -l)" -eq 23 ] || bad "the old genes are not all in the shipped table"

run32 "${WORK}/old.tsv" alt old
expect "old table: exit 0" "$RC" 0
expect "old table: BAM-mode calls as before" "$(bam_calls old)" "$OLD_BAM"
expect "old table: VCF-only calls as before" "$(vcf_calls old)" "$OLD_VCF"
expect "old table: depth prepared for the old bam genes" "$(prep old)" "$OLD_BAM"
expect "old table: summary --bam-genes as before" "$(sum_arg old bam-genes)" "$OLD_BAM"
expect "old table: summary --vcf-genes as before" "$(sum_arg old vcf-genes)" "$OLD_VCF"

# Each bam row the shipped table adds, alone on top of the old table.
for g in $(rows "$CAT" | awk -F'\t' '$2 == "bam" {print $1}'); do
  case " $OLD_BAM " in *" $g "*) continue ;; esac
  { awk -F'\t' -v g="$g" '$1 != g' "${WORK}/old.tsv"; row_of "$g"; } > "${WORK}/plus-${g}.tsv"
  run32 "${WORK}/plus-${g}.tsv" alt "plus-${g}"
  n_old=$(wc -w <<<"$OLD_BAM")
  expect "+${g}: exactly one more BAM-mode call, for ${g}" "$(bam_calls "plus-${g}")" "${OLD_BAM} ${g}"
  if vcf_calls "plus-${g}" | tr ' ' '\n' | grep -qxF "$g"; then bad "+${g}: also called from the VCF only"; fi
  [ "$(wc -w <<<"$(bam_calls "plus-${g}")")" -eq $((n_old + 1)) ] || bad "+${g}: not one more BAM-mode call"
done

# The shipped table: one BAM-mode call per bam row, one VCF-only call per vcf row.
WANT_BAM=$(rows "$CAT" | awk -F'\t' '$2 == "bam" {print $1}' | paste -sd ' ' -)
WANT_VCF=$(rows "$CAT" | awk -F'\t' '$2 == "vcf" {print $1}' | paste -sd ' ' -)
run32 "$CAT" alt shipped
expect "shipped table, ALT header: BAM-mode calls" "$(bam_calls shipped)" "$WANT_BAM"
expect "shipped table, ALT header: VCF-only calls" "$(vcf_calls shipped)" "$WANT_VCF"
expect "shipped table: 13 BAM-mode calls" "$(wc -w <<<"$(bam_calls shipped)" | tr -d ' ')" 13
if grep -q 'not callable on this reference' "${WORK}/shipped.out"; then bad "ALT header: a not-callable line"; else ok "ALT header: no not-callable line"; fi
if grep -q -- '--not-callable' "${WORK}/shipped.summary"; then bad "ALT header: --not-callable passed"; fi

run32 "$CAT" noalt noalt
NOALT_BAM=$(tr ' ' '\n' <<<"$WANT_BAM" | grep -vxF GSTT1 | paste -sd ' ' -)
expect "no-ALT header: exit 0" "$RC" 0
expect "no-ALT header: GSTT1 left out of the BAM-mode calls" "$(bam_calls noalt)" "$NOALT_BAM"
expect "no-ALT header: GSTT1 left out of depth preparation" "$(prep noalt)" "$NOALT_BAM"
expect "no-ALT header: the not-callable line" \
  "$(grep -c '^NOTICE: the BAM has no chr22_KI270879v1_alt contig (reference without ALT contigs); GSTT1 is not callable on this reference and is skipped$' "${WORK}/noalt.out")" 1
if grep -q -- '--not-callable GSTT1' "${WORK}/noalt.summary"; then ok "no-ALT header: --not-callable GSTT1 to the summary"; else bad "no-ALT header: no --not-callable GSTT1"; fi
expect "no-ALT header: the summary still lists GSTT1" "$(sum_arg noalt bam-genes)" "$WANT_BAM"

set_col GSTT1 3 yes > "${WORK}/gstt1-yes.tsv"
run32 "${WORK}/gstt1-yes.tsv" noalt gstt1-yes
if grep -q 'not callable on this reference' "${WORK}/gstt1-yes.out"; then bad "callable_without_alt=yes: a not-callable line"; else ok "control: with callable_without_alt=yes no not-callable line, so it comes from the column"; fi

{ cat "$CAT"; printf 'CYP2C9\tbam\tyes\t.\tnone\tno\n'; } > "${WORK}/dup.tsv"
run32 "${WORK}/dup.tsv" alt dup
expect "a gene twice stops the step" "$RC" 1
{ cat "$CAT"; printf 'ABCG2\tBAM\tyes\t.\tnone\tno\n'; } > "${WORK}/mode.tsv"
run32 "${WORK}/mode.tsv" alt mode
expect "an unknown mode stops the step" "$RC" 1

# --- 4. the module reads the same table ---------------------------------------
NF="${REPO}/modules/local/pypgx/main.nf"
# shellcheck disable=SC2016 # the ${projectDir} is matched, not expanded
grep -qF 'file("${projectDir}/assets/pgx_genes.tsv"' "$NF" || bad "the PYPGX module does not read assets/pgx_genes.tsv"
if grep -qE '"(CYP[0-9A-Z]+ ){2,}' "$NF"; then bad "the PYPGX module spells a gene list"; else ok "the PYPGX module spells no gene list"; fi

# --- 5. the reports ------------------------------------------------------------
# A summary with every gene of the table: GSTT1 not callable, CYP2E1 failed in
# pypgx, CYP2D6 withheld by the depth check, every other bam gene with a
# copy-number call. og summary collect reads it and og report render writes
# the text report, with a CNVpytor section beside pypgx's: every bam gene
# must have exactly one row in its PGx copy number block (its call, "not
# called: pypgx failed" or "not called: the CYP2D6 depth check failed"), and
# GSTT1 its not-assessable line. Two planted
# summaries must be reported: one where GSTT1 reads as called (no
# not-assessable line), and one where a bam gene's row says VCF (it drops out
# of the block).
OG_BIN=${OG:-${REPO}/bin/og}
# pypgx_table GSTT1_STATUS VCF_SOURCE_FOR: step 32's summary table for every gene.
pypgx_table() {
  printf 'Gene\tDiplotype\tPhenotype\tCNV_call\tSource\tStatus\tDefining_variants\tMin_GQ\tWeakest_variant\tDefining_filter\tQuality\n'
  rows "$CAT" | awk -F'\t' -v OFS='\t' -v gs="$1" -v vf="$2" '{
    g = $1; src = (g == vf) ? "VCF" : toupper($2); na5 = "N/A\tN/A\tN/A\tN/A\tN/A"
    if (g == "GSTT1") print g, "N/A", "N/A", "N/A", src, gs, na5
    else if (g == "CYP2E1") print g, "FAILED", "N/A", "N/A", src, "failed", na5
    else if (g == "CYP2D6") print g, "Indeterminate", "Indeterminate (CYP2D6 depth check)", "N/A (CYP2D6 depth check)", src, "no call (CYP2D6 depth check)", na5
    else print g, "*1/*1", "Normal Metabolizer", ($2 == "bam" ? "CNV-" g : "N/A"), src, "called", "none", "N/A", "N/A", "N/A", "N/A"
  }'
}
# report_problems NAME GSTT1_STATUS VCF_SOURCE_FOR: what is wrong with the reports of that summary, one per line.
report_problems() {
  local d="${WORK}/report-$1" g want got n
  mkdir -p "${d}/s1/pypgx"
  pypgx_table "$2" "$3" > "${d}/s1/pypgx/s1_pypgx_summary.tsv"
  "$OG_BIN" summary collect --sample s1 --sample-dir "${d}/s1" --out "${d}/collected.json" > /dev/null 2>&1 \
    || { echo "og summary collect failed"; return; }
  # CNVpytor's section beside pypgx's: pypgx owns every bam gene, so each
  # has one row, and it is pypgx's.
  jq '.sections.cnvpytor = {"title": "Copy number (CNVpytor)", "step": "18", "state": "ok", "source": null,
      "file_date": null, "note": null,
      "data": {"total": 0, "deletions": 0, "duplications": 0, "significant": 0}}' "${d}/collected.json" > "${d}/summary.json"
  "$OG_BIN" report render --summary "${d}/summary.json" -o "${d}/report.txt" > /dev/null 2>&1 \
    || { echo "og report render failed"; return; }
  # The block's rows as GENE<TAB>TEXT.
  sed -n '/^## PGx copy number$/,/^$/p' "${d}/report.txt" | sed -n 's/^  \([^ ]*\)  *\(.*\)$/\1\t\2/p' > "${d}/block.tsv"
  grep -qF 'could not be read' "${d}/report.txt" && echo "a report section could not be read"
  cut -f1 "${d}/block.tsv" | sort | uniq -c | awk '$1 > 1 { print $2 ": " $1 " rows in the copy-number block" }'
  while IFS= read -r g; do
    if [ "$g" = GSTT1 ]; then
      grep -qF 'Not assessable on this reference: GSTT1' "${d}/report.txt" || echo "GSTT1 has no not-assessable line"
      continue
    fi
    case $g in
      CYP2E1) want="not called: pypgx failed" ;;
      CYP2D6) want="not called: the CYP2D6 depth check failed" ;;
      *) want="CNV-${g}" ;;
    esac
    want="${want} (pypgx, read depth)"
    got=$(awk -F'\t' -v g="$g" '$1 == g { print $2; exit }' "${d}/block.tsv")
    [ "$got" = "$want" ] || echo "${g}: copy-number block says '${got}', not '${want}'"
  done < <(rows "$CAT" | awk -F'\t' '$2 == "bam" { print $1 }')
}
if ! "$OG_BIN" version > /dev/null 2>&1; then
  if [ -n "${OG_REQUIRED:-}" ]; then bad "the reports: og is not built (set OG or run scripts/setup.sh --og)"; else
    echo "skip: the reports need og (set OG or run scripts/setup.sh --og)"
  fi
else
  got=$(report_problems shipped "not callable on this reference" "")
  if [ -z "$got" ]; then
    n=$(rows "$CAT" | awk -F'\t' '$2 == "bam"' | wc -l)
    ok "the reports show the copy-number call of $((n - 3)) bam genes, CYP2E1's failed call, CYP2D6's depth-check no call and GSTT1's not-assessable line"
  else
    while IFS= read -r l; do bad "$l"; done <<< "$got"
  fi
  got=$(report_problems gstt1-called "no call" "")
  if [ -n "$got" ]; then ok "control, GSTT1 read as called reported: $(tr '\n' ';' <<< "$got")"; else bad "control, GSTT1 read as called passed"; fi
  got=$(report_problems slc22a2-vcf "not callable on this reference" SLC22A2)
  if [ -n "$got" ]; then ok "control, SLC22A2's row says VCF reported: $(tr '\n' ';' <<< "$got")"; else bad "control, SLC22A2's row says VCF passed"; fi
fi

# --- 6. one owner per gene ------------------------------------------------------
# A gene pypgx calls from the BAM (a bam row) is pypgx's: the depth engine's
# catalog has no row for it. A copy of the catalog with a UGT2B17 row planted
# back must be reported.
CN="${REPO}/bin/cn_loci.tsv"
cn_owner_clash() {   # prints the catalog loci that are bam rows of the table
  awk -F'\t' 'NR == FNR { if ($0 !~ /^#/ && $2 == "bam") bam[$1] = 1; next }
              $0 !~ /^#/ && $1 != "locus" && ($1 in bam) { print $1 }' "$CAT" "$1"
}
if [ ! -s "$CN" ]; then
  bad "bin/cn_loci.tsv is missing or empty"
else
  clash=$(cn_owner_clash "$CN")
  if [ -n "$clash" ]; then bad "bin/cn_loci.tsv has a row for a gene pypgx calls from the BAM: $clash"
  else ok "bin/cn_loci.tsv has no row for any of the $(awk -F'\t' '$2 == "bam"' "$CAT" | wc -l | tr -d ' ') bam genes"; fi
  { cat "$CN"; printf 'UGT2B17\tgene\tchr4\tUGT2B17:68516300-68627039\n'; } > "${WORK}/planted-cn-ugt2b17.tsv"
  clash=$(cn_owner_clash "${WORK}/planted-cn-ugt2b17.tsv")
  if [ "$clash" = "UGT2B17" ]; then ok "control, a planted UGT2B17 catalog row reported"
  else bad "control, a planted UGT2B17 catalog row passed (got '${clash}')"; fi
fi

# --- 7. no old call survives a rerun --------------------------------------------
# Before, step 32 removed only the summary: a gene this run did not call kept
# the last run's <gene>/results.zip, and the summary read it as this run's
# call; a gene pypgx failed on after writing part of its result kept that.
P="${G}/${S}/pypgx"
plant() { mkdir -p "${P}/$1"; echo old > "${P}/$1/results.zip"; }
plant CYP2D6; plant CYP2A6; plant CYP2C9
for z in depth_of_coverage.zip control_statistics.zip; do echo old > "${P}/${z}"; done
FAKE_PYPGX_FAIL=prep run32 "$CAT" alt stale-prep
expect "depth preparation failed: exit 0 (the VCF genes ran)" "$RC" 0
expect "depth preparation failed: no BAM-mode call" "$(bam_calls stale-prep)" ""
for g in CYP2D6 CYP2A6; do
  if [ -e "${P}/${g}/results.zip" ]; then bad "depth preparation failed: the earlier run's ${g}/results.zip is still there"
  else ok "depth preparation failed: no ${g}/results.zip from the earlier run"; fi
done
expect "depth preparation failed: CYP2C9 has this run's result" "$(cat "${P}/CYP2C9/results.zip" 2>/dev/null)" new
for z in depth_of_coverage.zip control_statistics.zip; do
  if [ -e "${P}/${z}" ]; then bad "depth preparation failed: the earlier run's ${z} is still there"
  else ok "depth preparation failed: no ${z} from the earlier run"; fi
done

plant CYP2C9
FAKE_PYPGX_FAIL=CYP2C9 run32 "$CAT" alt half-written
expect "a gene that fails: exit 0 (the other genes ran)" "$RC" 0
if [ -e "${P}/CYP2C9/results.zip" ]; then bad "a gene that failed after writing part of its result has a results.zip: $(cat "${P}/CYP2C9/results.zip")"
else ok "a gene that failed after writing part of its result has no results.zip"; fi
if [ -e "${P}/CYP2C9.part" ]; then bad "a gene that failed leaves CYP2C9.part"; fi
expect "a gene that succeeds: its result is under CYP2D6/" "$(cat "${P}/CYP2D6/results.zip" 2>/dev/null)" new
expect "a gene that succeeds: no CYP2D6.part left" "$(compgen -G "${P}/*.part" || true)" ""

# --- 8. an ALT contig check that cannot run stops the step ------------------------
# Before, the step warned and went on, so a BAM that might have ALT contigs
# reached pypgx's copy-number calls.
FAKE_HEADER_FAIL=1 run32 "$CAT" alt header-fail
expect "the BAM header cannot be read: exit 2" "$RC" 2
expect "the BAM header cannot be read: pypgx never runs" "$(wc -l < "${WORK}/header-fail.pypgx" | tr -d ' ')" 0
if grep -q '^ERROR: could not read the header of .*ALT contig check could not run' "${WORK}/header-fail.out"; then
  ok "the BAM header cannot be read: the ERROR line"
else bad "the BAM header cannot be read: no ERROR line: $(cat "${WORK}/header-fail.out")"; fi
FAKE_CHECK_FAIL=1 run32 "$CAT" alt check-fail
expect "og cn check-bam fails: exit 2" "$RC" 2
expect "og cn check-bam fails: pypgx never runs" "$(wc -l < "${WORK}/check-fail.pypgx" | tr -d ' ')" 0

# --- 9. the weakly supported genes in the doc --------------------------------------
# From pypgx's gene table and assets/pgx_genes.tsv: the genes of the table
# pypgx has no phenotype rule for (PhenotypeMethod N/A), and of those the ones
# no other tool here calls (second_caller none). docs/32-pypgx.md lists each
# set on its own line; a gene that changes either way fails until the doc
# follows.
DOC="${REPO}/docs/32-pypgx.md"
# no_rule [none]: the table's genes with PhenotypeMethod N/A (and second_caller none)
no_rule() {
  awk -v only="${1:-}" '
    NR == FNR {
      n = split($0, f, ",")
      if (FNR == 1) { for (i = 1; i <= n; i++) col[f[i]] = i; next }
      pm[f[col["Gene"]]] = f[col["PhenotypeMethod"]]; next
    }
    { split($0, t, "\t") }
    /^#/ || t[1] == "gene" || !NF { next }
    pm[t[1]] == "N/A" && (only == "" || t[5] == only) { print t[1] }
  ' "$GT" "$CAT" | sort | paste -sd ' ' -
}
# doc_list LEAD: the genes on the doc line that starts with **LEAD (N):**, sorted;
# doc_count LEAD: its N
doc_list() { sed -n "s/^\*\*$1 ([0-9]*):\*\* //p" "$DOC" | head -n 1 | tr -d '.' | tr ',' '\n' | tr -d ' ' | grep . | sort | paste -sd ' ' -; }
doc_count() { sed -n "s/^\*\*$1 (\([0-9]*\)):\*\*.*/\1/p" "$DOC" | head -n 1; }
for spec in "No pypgx phenotype rule|" "Weakly supported, pypgx alone|none"; do
  lead=${spec%|*} only=${spec#*|}
  want=$(no_rule "$only")
  [ -n "$want" ] || bad "no gene of the table lacks a pypgx phenotype rule: did gene-table.csv change shape?"
  expect "docs/32-pypgx.md: ${lead}" "$(doc_list "$lead")" "$want"
  expect "docs/32-pypgx.md: ${lead}, the count" "$(doc_count "$lead")" "$(wc -w <<<"$want" | tr -d ' ')"
done
# A control: a gene table where CYP2C9 has no phenotype rule must change the
# weakly supported set (CYP2C9's second caller is PharmCAT, so only the first).
sed 's/^CYP2C9,\([^,]*,[^,]*,[^,]*,[^,]*,[^,]*\),Score,/CYP2C9,\1,N\/A,/' "$GT" > "${WORK}/gt-cyp2c9.csv"
if [ "$(GT="${WORK}/gt-cyp2c9.csv" no_rule)" != "$(no_rule)" ]; then ok "control, CYP2C9 without a phenotype rule changes the no-rule set"
else bad "control, CYP2C9 without a phenotype rule left the no-rule set as it was"; fi

if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "all checks passed"
