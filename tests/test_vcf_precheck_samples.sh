#!/usr/bin/env bash
# One samplesheet row is one sample: a VCF with more than one sample column
# is refused at intake, with the count and the names, and a single-sample VCF
# goes on. The same block refuses a VCF whose ##contig lengths are not
# GRCh38's (a GRCh37 file named the chr way or the Ensembl way), and lets a
# GRCh38 header through. A header without contig lengths gets its REF bases
# compared with the reference instead: more than 5% of the first 1000
# records different stops the run, 5% or less goes on. A GT key after
# another FORMAT key stops the run with exit 2 before a reader can mistake
# depth or quality for a genotype.
#
# Two places apply the rule, and both are run here on the synthetic VCFs in
# tests/fixtures/vcf/:
#   - VCF_PRECHECK (modules/local/vcf_precheck/main.nf): its script: block is
#     rendered as Nextflow would (Groovy escapes, the few ${...} it uses) and
#     run under bash -euo pipefail, with a fake bcftools on PATH that answers
#     the calls the block makes from the plain-text fixture;
#   - vcf_header_samples in scripts/lib/common.sh, which validate-setup.sh
#     uses (its own wiring is the validate-multisample-vcf fake-docker case).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FIX="${REPO}/tests/fixtures/vcf"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

# --- the script: block of VCF_PRECHECK, as a bash file ---------------------------
# Literal dollar signs must survive Groovy before the rendered bash runs.
# A raw regex end anchor parses in bash/awk but breaks a Nextflow GString.
# shellcheck disable=SC2016  # Match literal dollar signs in the Groovy source.
if sed -n '/^    script:/,/^    stub:/p' "${REPO}/modules/local/vcf_precheck/main.nf" \
   | grep -nE '(^|[^\\])\$([^[:alpha:]_{]|$)'; then
  fail "VCF_PRECHECK has an unescaped literal dollar in its Groovy script block"
  exit 1
fi
pass "VCF_PRECHECK escapes literal dollar signs for Groovy"

# The Groovy values the block interpolates. A new one fails here, so this
# test is updated instead of running a script with a hole in it. Groovy's
# escapes in a triple-quoted GString go first (an unknown one is refused),
# then each ${...} that is not escaped; then the common indent comes off.
# shellcheck disable=SC2016  # the dollar signs are awk's and Groovy's
awk '
  function die(m) { print m > "/dev/stderr"; bad = 1; exit 1 }
  { text = text "\n" $0 }
  END {
    if (bad) exit 1
    text = text "\n"
    q = "\n    \"\"\"\n"
    a = index(text, "\n    script:\n")
    if (a) { rest = substr(text, a + 13); o = index(rest, q) }
    if (a && o) { rest = substr(rest, o + length(q) - 1); c = index(rest, q "\n    stub:") }
    if (!a || !o || !c) die("could not find the script: block of VCF_PRECHECK")
    body = substr(rest, 2, c - 2)
    v["vcf"] = "input.vcf.gz"; v["vcf.name"] = "input.vcf.gz"; v["meta.id"] = "S1"
    v["allow_unfiltered"] = "false"; v["task.process"] = "VCF_PRECHECK"
    v["reference"] = "ref.fa"; v["reference_fai"] = "ref.fa.fai"
    esc["\\"] = "\\"; esc["$"] = "$"; esc["n"] = "\n"; esc["t"] = "\t"; esc["\""] = "\""; esc["\047"] = "\047"
    out = ""; n = length(body)
    for (i = 1; i <= n; i++) {
      ch = substr(body, i, 1)
      if (ch == "\\" && i < n && substr(body, i + 1, 1) != "\n") {
        nx = substr(body, i + 1, 1)
        if (!(nx in esc)) die("unexpected Groovy escape \\" nx " in the script: block")
        out = out esc[nx]; i++
      } else if (substr(body, i, 2) == "${" && (j = index(substr(body, i + 2), "}"))) {
        e = substr(body, i + 2, j - 1)
        if (e ~ /^task\.container/) out = out "test"
        else if (e in v) out = out v[e]
        else die("unknown Groovy value in the script: block: ${" e "}")
        i += j + 1
      } else out = out ch
    }
    # textwrap.dedent: the common leading whitespace of the lines that are
    # not blank comes off every line; blank lines become empty.
    k = split(out, L, "\n"); margin = ""; first = 1
    for (i = 1; i <= k; i++) {
      if (L[i] !~ /[^ \t]/) continue
      match(L[i], /^[ \t]*/); w = substr(L[i], 1, RLENGTH)
      if (first) { margin = w; first = 0 }
      else while (substr(w, 1, length(margin)) != margin) margin = substr(margin, 1, length(margin) - 1)
    }
    for (i = 1; i <= k; i++) print (L[i] ~ /[^ \t]/ ? substr(L[i], length(margin) + 1) : "")
  }' "${REPO}/modules/local/vcf_precheck/main.nf" > "${WORK}/precheck.sh"

# --- a fake bcftools, from the plain-text VCF beside the input -----------------
mkdir -p "${WORK}/bin"
cat > "${WORK}/bin/bcftools" <<'FAKE'
#!/usr/bin/env bash
# Answers the calls VCF_PRECHECK makes. The input is gzip-compressed text.
case "$1 ${2:-}" in
  "index -s")
    gzip -dc "$3" | awk -F'\t' '!/^#/ { n[$1]++ } END { for (c in n) printf "%s\t.\t%d\n", c, n[c] }' ;;
  "view -H")
    gzip -dc "$3" | grep -v '^#' ;;
  "view -h")
    gzip -dc "$3" | grep '^#' ;;
  "head -h")
    # head -h 0 -n N FILE: the first N records, no header
    [ "$3 $4" = "0 -n" ] || { echo "fake bcftools: unexpected head call: $*" >&2; exit 2; }
    gzip -dc "$6" | grep -v '^#' | head -n "$5" ;;
  "norm -c")
    # norm -c x -f REF -o OUT IN: keep the records whose REF is the reference's
    # bases there. Like the real one it fails on a record whose contig the
    # header does not declare, or the reference does not have.
    [ "$3 $4 $6" = "x -f -o" ] || { echo "fake bcftools: unexpected norm call: $*" >&2; exit 2; }
    awk -F'\t' -v OFS='\t' -v out="$7" '
      FNR == 1 { f++ }
      f == 1 { if (/^>/) { name = substr($1, 2); sub(/[ \t].*/, "", name) } else seq[name] = seq[name] $0; next }
      /^##contig=<ID=/ { id = $0; sub(/^##contig=<ID=/, "", id); sub(/[,>].*/, "", id); declared[id] = 1 }
      /^#/ { print > out; next }
      !($1 in declared) { print "[E::vcf_format] Invalid BCF, CONTIG id=0 not present in the header" > "/dev/stderr"; exit 255 }
      !($1 in seq) { print "faidx_fetch_seq failed at " $1 ":" $2 > "/dev/stderr"; exit 255 }
      toupper(substr(seq[$1], $2, length($4))) == toupper($4) { print > out }' "$5" "$8" ;;
  "query -f")
    # The former query pass, kept so the pre-fix control runs on the same fake.
    case "$3" in
      *ALT*) gzip -dc "$4" | awk -F'\t' -v OFS='\t' '!/^#/ {
        end = "."
        n = split($8, info, ";")
        for (i = 1; i <= n; i++) if (info[i] ~ /^END=/) end = substr(info[i], 5)
        print $7, $5, end }' ;;
      *) gzip -dc "$4" | awk -F'\t' '!/^#/ { print $7 }' ;;
    esac ;;
  *) echo "fake bcftools: unexpected call: $*" >&2; exit 2 ;;
esac
FAKE
chmod +x "${WORK}/bin/bcftools"

# --- a small reference -------------------------------------------------------
# chr1 and chrM, every base G except the REF bases of the fixture's three
# records: chr1:10000100 A, chr1:10000200 C, chrM:100 T.
awk -v fa="${WORK}/ref.fa" -v fai="${WORK}/ref.fa.fai" '
  # seq NAME LENGTH: NAME of LENGTH bases, all G but those in base[]
  function seq(name, len,   s, p) {
    s = "G"
    while (length(s) < len) s = s s
    s = substr(s, 1, len)
    for (p in base) if (p + 0 <= len) s = substr(s, 1, p - 1) base[p] substr(s, p + 1)
    printf ">%s\n", name > fa
    off += length(name) + 2
    printf "%s\t%d\t%d\t%d\t%d\n", name, len, off, len, len + 1 > fai
    printf "%s\n", s > fa
    off += len + 1
    split("", base)
  }
  BEGIN {
    base[10000100] = "A"; base[10000200] = "C"; seq("chr1", 10001500)
    base[100] = "T"; seq("chrM", 16569)
  }'

# run_precheck VCF_TEXT_FILE: run the block in a fresh task directory; sets OUT, RC
run_precheck() {
  local dir
  dir=$(mktemp -d "${WORK}/task.XXXXXX")
  gzip -c "$1" > "${dir}/input.vcf.gz"
  : > "${dir}/input.vcf.gz.tbi"
  ln -s "${WORK}/ref.fa" "${dir}/ref.fa"
  ln -s "${WORK}/ref.fa.fai" "${dir}/ref.fa.fai"
  set +e
  OUT=$(cd "$dir" && PATH="${WORK}/bin:${PATH}" bash -euo pipefail "${WORK}/precheck.sh" 2>&1)
  RC=$?
  set -e
}

run_precheck "${FIX}/two_samples.vcf"
if [ "$RC" -ne 0 ] && grep -q "Sample 'S1': input.vcf.gz holds 2 samples (SAMPLE_A, SAMPLE_B)" <<<"$OUT" \
   && grep -q 'One samplesheet row is one sample' <<<"$OUT" \
   && grep -q 'bcftools view -s SAMPLE_A -a -c 1 -Oz -o SAMPLE_A.vcf.gz input.vcf.gz' <<<"$OUT"; then
  pass "VCF_PRECHECK refuses a two-sample VCF, naming the count, the samples and the fix"
else
  fail "VCF_PRECHECK on a two-sample VCF: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

awk -F'\t' -v OFS='\t' '/^#CHROM/ { $10 = "S1"; $11 = "S2\tS3\tS4\tS5\tS6\tS7" } { print }' \
  "${FIX}/two_samples.vcf" | grep '^#' > "${WORK}/seven.vcf"
printf 'chr1\t10000100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\n' >> "${WORK}/seven.vcf"
run_precheck "${WORK}/seven.vcf"
if [ "$RC" -ne 0 ] && grep -q 'holds 7 samples (S1, S2, S3, S4, S5, ...)' <<<"$OUT"; then
  pass "VCF_PRECHECK lists the first five of seven samples, then ..."
else
  fail "VCF_PRECHECK on a seven-sample VCF: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# A sample name is free text: the printed command quotes it, so copying it
# cannot run what the name holds.
# shellcheck disable=SC2016  # the $(...) is the literal sample name under test
awk -F'\t' -v OFS='\t' '/^#CHROM/ { $10 = "x$(touch pwned) y" } { print }' \
  "${FIX}/two_samples.vcf" > "${WORK}/odd.vcf"
run_precheck "${WORK}/odd.vcf"
# shellcheck disable=SC2016
if [ "$RC" -ne 0 ] && grep -qF 'bcftools view -s x\$\(touch\ pwned\)\ y -a -c 1 -Oz -o x\$\(touch\ pwned\)\ y.vcf.gz input.vcf.gz' <<<"$OUT"; then
  pass "VCF_PRECHECK shell-quotes the sample name in the command it prints"
else
  fail "VCF_PRECHECK on an odd sample name: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

run_precheck "${FIX}/one_sample.vcf"
if [ "$RC" -eq 0 ] && ! grep -q 'samples (' <<<"$OUT" && grep -q 'FILTER counts PASS=3' <<<"$OUT"; then
  pass "VCF_PRECHECK accepts a one-sample VCF and goes on to count FILTER values"
else
  fail "VCF_PRECHECK on a one-sample VCF: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# --- FORMAT order: never read a depth value as the genotype -------------------
# Only the second record is malformed: checking the header or first record
# alone cannot catch it. Repeated errors must name the first bad record.
awk -F'\t' -v OFS='\t' '
  !/^#/ && ++n >= 2 { $9 = "DP:GT:GQ"; $10 = "30:0/1:50" }
  { print }' "${FIX}/one_sample.vcf" > "${WORK}/gt-not-first.vcf"
run_precheck "${WORK}/gt-not-first.vcf"
if [ "$RC" -eq 2 ] && grep -q "Sample .S1.: input.vcf.gz" <<<"$OUT" \
   && grep -q 'chr1:10000200.*DP:GT:GQ' <<<"$OUT" \
   && grep -q 'GT must be the first FORMAT key' <<<"$OUT" \
   && grep -q 'FORMAT keys and each sample.*values together' <<<"$OUT"; then
  pass "VCF_PRECHECK refuses misplaced GT with exit 2, the first record and the fix"
else
  fail "VCF_PRECHECK on misplaced GT: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# GT may also be the last key, even on a filtered record.
awk -F'\t' -v OFS='\t' '
  !/^#/ && ++n == 3 { $7 = "LowQual"; $9 = "DP:GQ:GT"; $10 = "30:50:0/1" }
  { print }' "${FIX}/one_sample.vcf" > "${WORK}/gt-last.vcf"
run_precheck "${WORK}/gt-last.vcf"
if [ "$RC" -eq 2 ] && grep -q 'chrM:100.*DP:GQ:GT' <<<"$OUT"; then
  pass "VCF_PRECHECK checks the last record and filtered records for misplaced GT"
else
  fail "VCF_PRECHECK on GT last: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# The VCF rule applies when GT is present. A sites-only VCF, FORMAT '.', or
# FORMAT without GT is valid and must not be labelled a FORMAT-order error.
for format in GT GT:DP:GQ DP:GQ . XGT:DP; do
  awk -F'\t' -v OFS='\t' -v fmt="$format" '
    !/^#/ { $9 = fmt; $10 = (fmt == "GT" ? "0/1" : fmt == "GT:DP:GQ" ? "0/1:30:50" : fmt == "." ? "." : "30:50") }
    { print }' "${FIX}/one_sample.vcf" > "${WORK}/format.vcf"
  run_precheck "${WORK}/format.vcf"
  if [ "$RC" -eq 0 ] && grep -q 'FILTER counts PASS=3' <<<"$OUT"; then
    pass "VCF_PRECHECK preserves FORMAT $format"
  else
    fail "VCF_PRECHECK on FORMAT $format: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
  fi
done
cut -f1-8 "${FIX}/one_sample.vcf" > "${WORK}/sites-only.vcf"
run_precheck "${WORK}/sites-only.vcf"
if [ "$RC" -eq 0 ] && grep -q 'FILTER counts PASS=3' <<<"$OUT"; then
  pass "VCF_PRECHECK preserves a sites-only VCF"
else
  fail "VCF_PRECHECK on sites-only VCF: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# Changing the scan must preserve FILTER counts and reference-block detection.
# END may be surrounded by other INFO keys; a variant with END is not a block.
awk -F'\t' -v OFS='\t' '
  /^#CHROM/ { print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">" }
  !/^#/ {
    n++; $7 = (n == 1 ? "PASS" : n == 2 ? "." : "LowQual")
    if (n == 2) { $5 = "<NON_REF>"; $8 = "X=1;END=10000300;Y=2" }
    if (n == 3) $8 = "END=101"
  }
  { print }' "${FIX}/one_sample.vcf" > "${WORK}/blocks.vcf"
run_precheck "${WORK}/blocks.vcf"
if [ "$RC" -eq 0 ] && grep -q 'FILTER counts PASS=1 .=1 other=1' <<<"$OUT" \
   && grep -q 'gVCF blocks, reference-block records 1' <<<"$OUT"; then
  pass "VCF_PRECHECK preserves FILTER and reference-block counts in the shared scan"
else
  fail "VCF_PRECHECK on reference blocks: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# --- build: ##contig lengths that are not GRCh38's --------------------------------------
# GRCh37's chr1 under the chr name: the stop names the contig, both lengths and GRCh37.
sed 's/^##contig=<ID=chr1,length=248956422>/##contig=<ID=chr1,length=249250621>/' \
  "${FIX}/one_sample.vcf" > "${WORK}/grch37.vcf"
grep -q 'length=249250621' "${WORK}/grch37.vcf" || fail "the GRCh37 test input was not built"
run_precheck "${WORK}/grch37.vcf"
if [ "$RC" -ne 0 ] && grep -q "Sample 'S1': input.vcf.gz is not on GRCh38" <<<"$OUT" \
   && grep -q 'chr1 length 249250621 (GRCh38: 248956422)' <<<"$OUT" \
   && grep -q 'A chr1 length of 249250621 is GRCh37 (hg19)' <<<"$OUT"; then
  pass "VCF_PRECHECK stops a chr-named VCF with GRCh37's chr1 length, naming the build"
else
  fail "VCF_PRECHECK on a GRCh37 header: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# The same under Ensembl names: the build stop comes before the rename advice.
awk -F'\t' -v OFS='\t' '/^##contig=<ID=chr1,/ { $0 = "##contig=<ID=1,length=249250621>" } $1 == "chr1" { $1 = "1" } { print }' \
  "${FIX}/one_sample.vcf" > "${WORK}/grch37-ensembl.vcf"
run_precheck "${WORK}/grch37-ensembl.vcf"
if [ "$RC" -ne 0 ] && grep -q '    1 length 249250621 (GRCh38: 248956422)' <<<"$OUT"; then
  pass "VCF_PRECHECK stops an Ensembl-named VCF with GRCh37's chr1 length"
else
  fail "VCF_PRECHECK on an Ensembl GRCh37 header: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# One wrong length on another contig is enough, and the GRCh37 line is left out.
awk '/^#CHROM/ { print "##contig=<ID=chr2,length=243199373,assembly=hg19>" } { print }' \
  "${FIX}/one_sample.vcf" > "${WORK}/chr2.vcf"
run_precheck "${WORK}/chr2.vcf"
if [ "$RC" -ne 0 ] && grep -q 'chr2 length 243199373 (GRCh38: 242193529)' <<<"$OUT" \
   && ! grep -q 'is GRCh37' <<<"$OUT"; then
  pass "VCF_PRECHECK stops on a wrong chr2 length and names it"
else
  fail "VCF_PRECHECK on a wrong chr2 length: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# Every GRCh38 length of chr1-22, X and Y, in any attribute order, goes on;
# so does a header without ##contig lines (validate-setup.sh spot-checks those).
{
  grep '^##fileformat' "${FIX}/one_sample.vcf"
  awk -F'[(", )]+' '/^    \("chr/ { for (i = 2; i < NF; i += 2) print $i, $(i + 1) }' \
    "${REPO}/tests/demo/make_demo_sample.py" \
    | awk '{ printf "##contig=<ID=%s,assembly=GRCh38,length=%s>\n", $1, $2 }'
  grep -v '^##fileformat' "${FIX}/one_sample.vcf" | grep -v '^##contig'
} > "${WORK}/grch38-all.vcf"
[ "$(grep -c '^##contig' "${WORK}/grch38-all.vcf")" -eq 24 ] || fail "the GRCh38 test header does not have 24 contigs"
run_precheck "${WORK}/grch38-all.vcf"
if [ "$RC" -eq 0 ] && grep -q 'FILTER counts PASS=3' <<<"$OUT" && ! grep -q 'REF bases' <<<"$OUT"; then
  pass "VCF_PRECHECK accepts GRCh38's 24 lengths (the demo sample's NCBI table), with no REF spot-check"
else
  fail "VCF_PRECHECK on a full GRCh38 header: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# --- build: REF bases, when the header gives no contig length ------------------------
# The fixture's records match the small reference: the file goes on.
grep -v '^##contig' "${FIX}/one_sample.vcf" > "${WORK}/nocontig.vcf"
run_precheck "${WORK}/nocontig.vcf"
if [ "$RC" -eq 0 ] && ! grep -q 'not on GRCh38' <<<"$OUT" \
   && grep -q 'REF bases differ from the reference in 0 of the first 3 records (0.0%)' <<<"$OUT"; then
  pass "VCF_PRECHECK spot-checks a header without ##contig lines and lets matching REF bases through"
else
  fail "VCF_PRECHECK on a header without ##contig lines: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# spot_vcf FILE N_SAME N_DIFF [CONTIG]: a header without contig lengths and
# N_SAME records whose REF is the reference's G, then N_DIFF whose REF is T.
spot_vcf() {
  local c=${4:-chr1}
  {
    grep -v '^##contig' "${FIX}/one_sample.vcf" | grep '^#'
    awk -v s="$2" -v d="$3" -v c="$c" 'BEGIN { OFS = "\t"
      for (i = 1; i <= s + d; i++) print c, 10000300 + i, ".", (i <= s ? "G" : "T"), "A", 50, "PASS", ".", "GT", "0/1" }'
  } > "$1"
}

# Every REF different: the stop names the sample, the count, the rate and the fix.
spot_vcf "${WORK}/spot_all.vcf" 0 20
run_precheck "${WORK}/spot_all.vcf"
if [ "$RC" -ne 0 ] && grep -q "Sample 'S1': input.vcf.gz is not on GRCh38. Its header gives no contig lengths" <<<"$OUT" \
   && grep -q 'they differ in 20 of the first 20 records (100.0%)' <<<"$OUT" \
   && grep -q 'start from the FASTQ' <<<"$OUT" && ! grep -q 'FILTER counts' <<<"$OUT"; then
  pass "VCF_PRECHECK stops a header without contig lengths whose REF bases are not the reference's"
else
  fail "VCF_PRECHECK on REF bases all different: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# The line is validate-setup.sh's: 1 of 20 (5%) goes on; 1 of 19 (5.3%) and
# 2 of 20 (10%) stop, and the rate printed is above 5.
spot_vcf "${WORK}/spot_5.vcf" 19 1
run_precheck "${WORK}/spot_5.vcf"
if [ "$RC" -eq 0 ] && grep -q 'REF bases differ from the reference in 1 of the first 20 records (5.0%)' <<<"$OUT"; then
  pass "VCF_PRECHECK lets 5% different REF bases through"
else
  fail "VCF_PRECHECK on 1 of 20 REF bases different: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi
spot_vcf "${WORK}/spot_19.vcf" 18 1
run_precheck "${WORK}/spot_19.vcf"
if [ "$RC" -ne 0 ] && grep -q 'they differ in 1 of the first 19 records (5.3%)' <<<"$OUT"; then
  pass "VCF_PRECHECK stops at 1 of 19 different REF bases, and prints 5.3%"
else
  fail "VCF_PRECHECK on 1 of 19 REF bases different: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi
spot_vcf "${WORK}/spot_10.vcf" 18 2
run_precheck "${WORK}/spot_10.vcf"
if [ "$RC" -ne 0 ] && grep -q 'they differ in 2 of the first 20 records (10.0%)' <<<"$OUT"; then
  pass "VCF_PRECHECK stops at 10% different REF bases"
else
  fail "VCF_PRECHECK on 2 of 20 REF bases different: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# Only the first 1000 records count.
spot_vcf "${WORK}/spot_1000.vcf" 1000 200
run_precheck "${WORK}/spot_1000.vcf"
if [ "$RC" -eq 0 ] && grep -q 'differ from the reference in 0 of the first 1000 records' <<<"$OUT"; then
  pass "VCF_PRECHECK compares the first 1000 records only"
else
  fail "VCF_PRECHECK on 1200 records: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# Records on a contig the reference lacks are left out, not compared.
spot_vcf "${WORK}/spot_chr2.vcf" 0 20 chr2
run_precheck "${WORK}/spot_chr2.vcf"
if [ "$RC" -eq 0 ] && grep -q 'no record is on a reference contig: REF bases not compared' <<<"$OUT"; then
  pass "VCF_PRECHECK leaves out records on contigs the reference lacks"
else
  fail "VCF_PRECHECK on records off the reference: rc=${RC}, output: $(tr '\n' '|' <<<"$OUT")"
fi

# --- vcf_header_samples in scripts/lib/common.sh ------------------------------------
# shellcheck source=../scripts/lib/common.sh
. "${REPO}/scripts/lib/common.sh"
got=$(vcf_header_samples < "${FIX}/two_samples.vcf" | paste -sd, -)
if [ "$got" = "SAMPLE_A,SAMPLE_B" ]; then pass "vcf_header_samples: two samples"; else fail "vcf_header_samples two: '${got}'"; fi
got=$(vcf_header_samples < "${FIX}/one_sample.vcf" | paste -sd, -)
if [ "$got" = "SAMPLE_A" ]; then pass "vcf_header_samples: one sample"; else fail "vcf_header_samples one: '${got}'"; fi
got=$(cut -f1-8 "${FIX}/one_sample.vcf" | vcf_header_samples | paste -sd, -)
if [ -z "$got" ]; then pass "vcf_header_samples: a sites-only VCF has none"; else fail "vcf_header_samples sites-only: '${got}'"; fi

echo ""
if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "All checks passed"
