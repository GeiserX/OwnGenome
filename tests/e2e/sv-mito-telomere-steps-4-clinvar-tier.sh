#!/usr/bin/env bash
# Step 23's ClinVar tier follows the ClinVar file step 13 annotates with
# --custom (CSQ field ClinVar_CLNSIG), not the cache release's CLIN_SIG. The
# fixture's VEP output (built with --database, no --custom) gets the field
# added: record A is Pathogenic in ClinVar_CLNSIG and not in CLIN_SIG (newly
# classified), record B is pathogenic in CLIN_SIG only (reclassified since the
# cache release). On a copy of the sample with no step 06 hits, the tier must
# hold A and not B.
. "$(dirname "$0")/lib.sh"

T="${SAMPLE}cv"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/vep"
SRC="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf"
OUT="${GENOME_DIR}/${T}/vep/${T}_vep.vcf"
# The CSQ header gains ClinVar_CLNSIG (and CLIN_SIG when it lacks one); each
# PASS record's transcripts gain the field: Pathogenic on the first (A), empty
# on the rest, with CLIN_SIG set to "" on A and to "pathogenic" on the second
# (B). A and B go to ab.tsv.
awk -F'\t' -v OFS='\t' -v ab="${CASE_TMP}/ab.tsv" '
  function repl(s, x, y,   o, i) {  # every x in s replaced by y, as written
    o = ""
    while ((i = index(s, x)) > 0) { o = o substr(s, 1, i - 1) y; s = substr(s, i + length(x)) }
    return o s
  }
  /^##INFO=<ID=CSQ/ {
    if (match($0, /Format: [^"]+/)) {
      old = substr($0, RSTART + 8, RLENGTH - 8)
      n = split(old, f, "|"); ci = 0
      for (i = 1; i <= n; i++) if (f[i] == "CLIN_SIG") { ci = i; break }
      missing = !ci
      if (missing) ci = n + 1
      $0 = repl($0, old, old (missing ? "|CLIN_SIG" : "") "|ClinVar_CLNSIG")
    }
    print; next
  }
  /^#/ || $7 != "PASS" { print; next }
  {
    m = split($8, info, ";"); k = 0
    for (i = 1; i <= m; i++) if (substr(info[i], 1, 4) == "CSQ=") { k = i; break }
    if (a == "") { a = $1 "\t" $2; add = "Pathogenic"; cache = ""; set = 1 }
    else if (b == "") { b = $1 "\t" $2; add = ""; cache = "pathogenic"; set = 1 }
    else { add = ""; set = 0 }
    t = split(substr(info[k], 5), trs, ","); out = ""
    for (j = 1; j <= t; j++) {
      nv = split(trs[j], v, "|")
      if (missing) v[++nv] = ""
      if (set) v[ci] = cache
      s = v[1]; for (x = 2; x <= nv; x++) s = s "|" v[x]
      out = out (j > 1 ? "," : "") s "|" add
    }
    info[k] = "CSQ=" out
    s = info[1]; for (x = 2; x <= m; x++) s = s ";" info[x]
    $8 = s
    print
  }
  END { printf "%s\n%s\n", a, b > ab }' "$SRC" > "$OUT"
read -r A_CHR A_POS < <(sed -n 1p "${CASE_TMP}/ab.tsv")
read -r B_CHR B_POS < <(sed -n 2p "${CASE_TMP}/ab.tsv")
echo "A ${A_CHR}:${A_POS} (ClinVar_CLNSIG Pathogenic), B ${B_CHR}:${B_POS} (CLIN_SIG pathogenic only)"

run_step 23-clinical-filter.sh "$T"
check_step_exit 23-clinical-filter.sh
check "the log names ClinVar_CLNSIG as the tier's source" has 'ClinVar tier: ClinVar_CLNSIG' "$(cat "$STEP_LOG")"
TIER="${T}/clinical/${T}_clinvar_pathogenic.vcf.gz"
check "the ClinVar tier VCF is readable" vcf_ok "$TIER"
check_ge "record A (Pathogenic in the current ClinVar file) is in the tier" "$(vcf_count -r "${A_CHR}:${A_POS}" "$TIER")" 1
check_eq "record B (pathogenic in the cache only) is not" "$(vcf_count -r "${B_CHR}:${B_POS}" "$TIER")" 0

finish
